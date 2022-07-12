/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "testing mesh refinement algorithm\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
    ierr = PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
    CHKERRG(ierr);
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
              "*** ERROR -my_file (MESH FILE NEEDED)");
    }

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    const char *option;
    option = "";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    auto get_dim = [&]() {
      int dim;
      int nb_ents_3d;
      CHKERR m_field.get_moab().get_number_entities_by_dimension(
          0, 3, nb_ents_3d, true);
      if (nb_ents_3d > 0) {
        dim = 3;
      } else {
        int nb_ents_2d;
        CHKERR m_field.get_moab().get_number_entities_by_dimension(
            0, 2, nb_ents_2d, true);
        if (nb_ents_2d > 0) {
          dim = 2;
        } else {
          dim = 1;
        }
      }
      return dim;
    };

    auto dim = get_dim();

    auto refine = m_field.getInterface<MeshRefinement>();

    BitRefLevel bit_level0;
    bit_level0.set(0);
    if (dim == 3 || dim == 2) {
      CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
          0, dim, bit_level0);
    } else {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
              "Dimension not handled by test");
    }

    BitRefLevel bit_level1;
    bit_level1.set(1);

    auto refine_edges = [&](auto bit0, auto bit) {
      MoFEMFunctionBegin;
      auto meshset_ptr = get_temp_meshset_ptr(moab);
      CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByRefLevel(
          bit0, BitRefLevel().set(), *meshset_ptr);

      // random mesh refinement
      auto meshset_ref_edges_ptr = get_temp_meshset_ptr(moab);
      Range edges_to_refine;
      CHKERR moab.get_entities_by_type(*meshset_ptr, MBEDGE, edges_to_refine);
      int ii = 0;
      for (Range::iterator eit = edges_to_refine.begin();
           eit != edges_to_refine.end(); eit++, ii++) {
        int numb = ii % 2;
        if (numb == 0) {
          CHKERR moab.add_entities(*meshset_ref_edges_ptr, &*eit, 1);
        }
      }
      CHKERR refine->addVerticesInTheMiddleOfEdges(*meshset_ref_edges_ptr, bit,
                                                   false, QUIET, 10000);
      if (dim == 3) {
        CHKERR refine->refineTets(*meshset_ptr, bit, QUIET, true);
      } else if (dim == 2) {
        CHKERR refine->refineTris(*meshset_ptr, bit, QUIET, true);
      } else {
        SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                "Dimension not handled by test");
      }
      MoFEMFunctionReturn(0);
    };

    auto refine_ents_hanging_nodes = [&](auto bit0, auto bit) {
      MoFEMFunctionBegin;
      auto meshset_ptr = get_temp_meshset_ptr(moab);
      CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByRefLevel(
          bit0, BitRefLevel().set(), *meshset_ptr);

      // random mesh refinement
      auto meshset_ref_edges_ptr = get_temp_meshset_ptr(moab);
      Range ents_dim;
      CHKERR moab.get_entities_by_dimension(*meshset_ptr, dim, ents_dim);
      int ii = 0;
      for (Range::iterator eit = ents_dim.begin(); eit != ents_dim.end();
           eit++, ii++) {
        int numb = ii % 2;
        if (numb == 0) {
          std::vector<EntityHandle> adj_ents;
          CHKERR moab.get_adjacencies(&*eit, 1, 1, false, adj_ents);
          CHKERR moab.add_entities(*meshset_ref_edges_ptr, &*adj_ents.begin(),
                                   adj_ents.size());
        }
      }

      CHKERR refine->addVerticesInTheMiddleOfEdges(*meshset_ref_edges_ptr, bit,
                                                   false, QUIET, 10000);
      if (dim == 3) {
        CHKERR refine->refineTetsHangingNodes(*meshset_ptr, bit, QUIET, true);
      } else if (dim == 2) {
        CHKERR refine->refineTrisHangingNodes(*meshset_ptr, bit, QUIET, true);
      } else {
        SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                "Dimension not handled by test");
      }
      MoFEMFunctionReturn(0);
    };

    auto save_blessed_field = [&](auto bit) {
      MoFEMFunctionBegin;

      std::ofstream myfile;
      myfile.open("mesh_refine.txt");

      auto out_meshset_tet_ptr = get_temp_meshset_ptr(moab);
      CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByDimAndRefLevel(
          bit, BitRefLevel().set(), dim, *out_meshset_tet_ptr);
      Range tets;
      CHKERR moab.get_entities_by_handle(*out_meshset_tet_ptr, tets);
      {
        int ii = 0;
        for (Range::iterator tit = tets.begin(); tit != tets.end(); tit++) {
          int num_nodes;
          const EntityHandle *conn;
          CHKERR moab.get_connectivity(*tit, conn, num_nodes, true);

          for (int nn = 0; nn < num_nodes; nn++) {
            // cout << conn[nn] << " ";
            myfile << conn[nn] << " ";
          }
          // cout << std::endl;
          myfile << std::endl;
          if (ii > 25)
            break;
        }
      }

      myfile.close();

      MoFEMFunctionReturn(0);
    };

    auto save_vtk = [&](auto bit) {
      MoFEMFunctionBegin;
      auto out_meshset_tet_ptr = get_temp_meshset_ptr(moab);
      CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByDimAndRefLevel(
          bit, BitRefLevel().set(), dim, *out_meshset_tet_ptr);
      CHKERR moab.write_file("out_mesh_refine.vtk", "VTK", "",
                             out_meshset_tet_ptr->get_ptr(), 1);
      MoFEMFunctionReturn(0);
    };

    CHKERR refine_edges(BitRefLevel().set(0), BitRefLevel().set(1));
    CHKERR refine_ents_hanging_nodes(BitRefLevel().set(1),
                                     BitRefLevel().set(2));
    CHKERR save_blessed_field(BitRefLevel().set(2));
    CHKERR save_vtk(BitRefLevel().set(2));
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
