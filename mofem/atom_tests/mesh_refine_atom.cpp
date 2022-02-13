/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

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
    if (dim == 3) {
      CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
          0, 3, bit_level0);
    } else if (dim == 2) {
      CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
          0, 2, bit_level0);
    } else {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
              "Dimension not handled by test");
    }

    BitRefLevel bit_level1;
    bit_level1.set(1);

    auto meshset_level0_ptr = get_temp_meshset_ptr(moab);
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByRefLevel(
        bit_level0, BitRefLevel().set(), *meshset_level0_ptr);

    // random mesh refinement
    auto meshset_ref_edges_ptr = get_temp_meshset_ptr(moab);
    Range edges_to_refine;
    CHKERR moab.get_entities_by_type(*meshset_level0_ptr, MBEDGE,
                                     edges_to_refine);
    int ii = 0;
    for (Range::iterator eit = edges_to_refine.begin();
         eit != edges_to_refine.end(); eit++, ii++) {
      int numb = ii % 2;
      if (numb == 0) {
        CHKERR moab.add_entities(*meshset_ref_edges_ptr, &*eit, 1);
      }
    }
    CHKERR refine->addVerticesInTheMiddleOfEdges(
        *meshset_ref_edges_ptr, bit_level1, false, QUIET, 10000);
    if (dim == 3) {
      CHKERR refine->refineTets(*meshset_level0_ptr, bit_level1, QUIET, true);
    } else if (dim == 2) {
      CHKERR refine->refineTris(*meshset_level0_ptr, bit_level1, QUIET, true);
    } else {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
              "Dimension not handled by test");
    }

    std::ofstream myfile;
    myfile.open("mesh_refine.txt");

    auto out_meshset_tet_ptr = get_temp_meshset_ptr(moab);
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByDimAndRefLevel(
        bit_level1, BitRefLevel().set(), dim, *out_meshset_tet_ptr);
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

    CHKERR moab.write_file("out_mesh_refine.vtk", "VTK", "",
                           out_meshset_tet_ptr->get_ptr(), 1);
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
