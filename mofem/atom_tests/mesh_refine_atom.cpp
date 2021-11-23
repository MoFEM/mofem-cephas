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

    MeshRefinement *refine;
    CHKERR m_field.getInterface(refine);

    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);

    BitRefLevel bit_level1;
    bit_level1.set(1);

    EntityHandle meshset_level0;
    CHKERR moab.create_meshset(MESHSET_SET, meshset_level0);
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByRefLevel(
        bit_level0, BitRefLevel().set(), meshset_level0);

    // random mesh refinement
    EntityHandle meshset_ref_edges;
    CHKERR moab.create_meshset(MESHSET_SET, meshset_ref_edges);
    Range edges_to_refine;
    CHKERR moab.get_entities_by_type(meshset_level0, MBEDGE, edges_to_refine);
    int ii = 0;
    for (Range::iterator eit = edges_to_refine.begin();
         eit != edges_to_refine.end(); eit++, ii++) {
      int numb = ii % 2;
      if (numb == 0) {
        CHKERR moab.add_entities(meshset_ref_edges, &*eit, 1);
      }
    }
    CHKERR refine->addVerticesInTheMiddleOfEdges(
        meshset_ref_edges, bit_level1, false, QUIET, 10000);
    CHKERR refine->refineTets(meshset_level0, bit_level1, false, QUIET);

    std::ofstream myfile;
    myfile.open("mesh_refine.txt");

    EntityHandle out_meshset_tet;
    CHKERR moab.create_meshset(MESHSET_SET, out_meshset_tet);
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit_level1, BitRefLevel().set(), MBTET, out_meshset_tet);
    Range tets;
    CHKERR moab.get_entities_by_handle(out_meshset_tet, tets);
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

    CHKERR moab.write_file("out_mesh_refine.vtk", "VTK", "", &out_meshset_tet,
                           1);

  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
