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

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, PETSC_NULL, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Reade parameters from line command
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
    PetscInt order;
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-my_order", &order, &flg);
#else
    CHKERR PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-my_order", &order,
                              &flg);
#endif
    if (flg != PETSC_TRUE) {
      order = 1;
    }

    // Read mesh to MOAB
    const char *option;
    option = ""; //"PARALLEL=BCAST;";//;DEBUG_IO";
    CHKERR moab.load_file(mesh_file_name, 0, option);
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

    // Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // stl::bitset see for more details
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle meshset_level0;
    CHKERR moab.create_meshset(MESHSET_SET, meshset_level0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByRefLevel(
        bit_level0, BitRefLevel().set(), meshset_level0);

    /***/
    // Define problem

    // Coord system
    CoordSystemsManager *cs_manger_ptr;
    CHKERR m_field.getInterface(cs_manger_ptr);
    {
      int cs_dim[] = {0, 3, 0, 3};
      CHKERR cs_manger_ptr->addCoordinateSystem(cs_dim,
                                                "BASE_FOR_TWO_POINT_TENSOR");
    }

    // Fields
    CHKERR m_field.add_field("FIELD_A", H1, AINSWORTH_LEGENDRE_BASE, 9);
    CHKERR cs_manger_ptr->setFieldCoordinateSystem("FIELD_A",
                                                   "BASE_FOR_TWO_POINT_TENSOR");

    CHKERR m_field.add_ents_to_field_by_type(0, MBTET, "FIELD_A");

    CHKERR m_field.set_field_order(0, MBTET, "FIELD_A", order);
    CHKERR m_field.set_field_order(0, MBTRI, "FIELD_A", order);
    CHKERR m_field.set_field_order(0, MBEDGE, "FIELD_A", order);
    CHKERR m_field.set_field_order(0, MBVERTEX, "FIELD_A", 1);

    // build field
    CHKERR m_field.build_fields();

    int cs_dim[4];
    std::string cs_name;

    // Open mesh_file_name.txt for writing
    for (_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field, "FIELD_A", dof_ptr)) {

      for (int alpha = 0; alpha < 4; alpha++) {
        cs_dim[alpha] = (*dof_ptr)->getCoordSysDim(alpha);
      }

      if (cs_dim[1] != 3) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong base dim");
      }
      if (cs_dim[3] != 3) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong base dim");
      }

      break;
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  return 0;
}
