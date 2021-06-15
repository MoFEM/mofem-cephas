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

static char help[] = "...\n\n";
static int debug = 1;

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Read parameters from line command
    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
#else
    CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-my_file",
                                 mesh_file_name, 255, &flg);
#endif
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");
    }

    // Read mesh to MOAB
    const char *option;
    option = ""; //"PARALLEL=BCAST;";//;DEBUG_IO";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // Create MoFEM (Joseph) databas
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);

    NodeMergerInterface *node_merger_iface;
    CHKERR m_field.getInterface(node_merger_iface);

    int ii = 1;
    for (; ii < 2; ii++) {

      BitRefLevel bit_level1;
      bit_level1.set(ii);

      Range edges;
      // CHKERR moab.get_entities_by_type(0,MBEDGE,edges,false);
      CHKERR m_field.getInterface<BitRefManager>()
          ->getEntitiesByTypeAndRefLevel(bit_level0, BitRefLevel().set(),
                                         MBEDGE, edges);
      Range::iterator eit = edges.begin();

      const EntityHandle *conn;
      int num_nodes;
      CHKERR moab.get_connectivity(*eit, conn, num_nodes, true);
      CHKERR node_merger_iface->mergeNodes(conn[0], conn[1], bit_level1,
                                           bit_level0);

      EntityHandle meshset_level1;
      CHKERR moab.create_meshset(MESHSET_SET, meshset_level1);
      CHKERR m_field.getInterface<BitRefManager>()
          ->getEntitiesByTypeAndRefLevel(bit_level1, BitRefLevel().set(), MBTET,
                                         meshset_level1);

      std::ostringstream ss;
      ss << "node_merger_" << ii << ".vtk";

      if (debug)
        CHKERR moab.write_file(ss.str().c_str(), "VTK", "", &meshset_level1, 1);
      bit_level0 = bit_level1;
    }

    Range tets;
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        BitRefLevel().set(ii - 1), BitRefLevel().set(), MBTET, tets);

    std::cout << tets << std::endl;
    if (tets.size() != 10) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
               "diffrent number of tets than expected = %u", tets.size());
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  return 0;
}
