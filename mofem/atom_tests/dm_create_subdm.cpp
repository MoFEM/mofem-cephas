/** \file dm_create_subdm.cpp
  \example dm_create_subdm.cpp
  \brief Atom test for Data Manager Interface and create sub-problem
*/

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

static const bool debug = false;

int main(int argc, char *argv[]) {

  // initialize petsc
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

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

    // register new dm type, i.e. mofem
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);
    DMType dm_name_sub = "DMSUB";
    CHKERR DMRegister_MoFEM(dm_name_sub);

    // craete dm instance
    DM dm;
    CHKERR DMCreate(PETSC_COMM_WORLD, &dm);
    CHKERR DMSetType(dm, dm_name);

    DM subdm0, subdm1;
    CHKERR DMCreate(PETSC_COMM_WORLD, &subdm0);
    CHKERR DMSetType(subdm0, dm_name_sub);
    CHKERR DMCreate(PETSC_COMM_WORLD, &subdm1);
    CHKERR DMSetType(subdm1, dm_name_sub);

    // read mesh and create moab and mofem data structures
    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    auto moab_comm_wrap =
        boost::make_shared<WrapMPIComm>(PETSC_COMM_WORLD, false);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, moab_comm_wrap->get_comm());

    const char *option;
    option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION="
             "PARALLEL_PARTITION;";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    EntityHandle root_set = moab.get_root_set();
    // add all entities to database, all of them will be used
    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        root_set, 3, bit_level0);
    // define & build field
    const int field_rank = 1;
    CHKERR m_field.add_field("FIELD0", H1, AINSWORTH_LEGENDRE_BASE, field_rank);
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD0");
    // set app. order
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD0", 1);
    CHKERR m_field.add_field("FIELD1", H1, AINSWORTH_LEGENDRE_BASE, field_rank);
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD1");
    // set app. order
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD1", 1);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD1", 2);

    // build data structures for fields
    CHKERR m_field.build_fields();

    // define & build finite elements
    CHKERR m_field.add_finite_element("FE00");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("FE00", "FIELD0");
    CHKERR m_field.modify_finite_element_add_field_col("FE00", "FIELD0");
    CHKERR m_field.modify_finite_element_add_field_data("FE00", "FIELD0");
    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "FE00");

    // define & build finite elements
    CHKERR m_field.add_finite_element("FE11");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("FE11", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("FE11", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("FE11", "FIELD1");
    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "FE11");

    // define & build finite elements
    CHKERR m_field.add_finite_element("FE01");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("FE01", "FIELD0");
    CHKERR m_field.modify_finite_element_add_field_col("FE01", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("FE01", "FIELD0");
    CHKERR m_field.modify_finite_element_add_field_data("FE01", "FIELD1");
    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "FE01");

    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);

    // set dm data structure which created mofem data structures
    CHKERR DMMoFEMCreateMoFEM(dm, &m_field, "MAIN_PROBLEM", bit_level0);
    CHKERR DMSetFromOptions(dm);
    CHKERR DMMoFEMAddElement(dm, "FE00");
    CHKERR DMMoFEMAddElement(dm, "FE11");
    CHKERR DMMoFEMAddElement(dm, "FE01");
    CHKERR DMSetUp(dm);
    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>(
            "MAIN_PROBLEM", -1, -1, 1);

    int nf;
    char **field_names;
    IS *fields;
    CHKERR DMCreateFieldIS(dm, &nf, &field_names, &fields);
    for (int f = 0; f != nf; f++) {
      PetscPrintf(PETSC_COMM_WORLD, "%d field %s\n", f, field_names[f]);
      CHKERR PetscFree(field_names[f]);
      CHKERR ISDestroy(&(fields[f]));
    }
    CHKERR PetscFree(field_names);
    CHKERR PetscFree(fields);

    CHKERR DMMoFEMCreateSubDM(subdm0, dm, "SUB0");
    CHKERR DMMoFEMSetSquareProblem(subdm0, PETSC_TRUE);
    CHKERR DMMoFEMAddElement(subdm0, "FE11");
    CHKERR DMMoFEMAddSubFieldRow(subdm0, "FIELD1", MBVERTEX, MBVERTEX);
    CHKERR DMMoFEMAddSubFieldCol(subdm0, "FIELD1", MBVERTEX, MBVERTEX);
    CHKERR DMSetUp(subdm0);
    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>("SUB0", -1,
                                                                   -1, 1);
    if (debug) {
      Mat A;
      CHKERR DMCreateMatrix(subdm0, &A);
      MatView(A, PETSC_VIEWER_DRAW_WORLD);
      std::string wait;
      std::cin >> wait;
      CHKERR MatDestroy(&A);
    }

    CHKERR DMMoFEMCreateSubDM(subdm1, dm, "SUB1");
    CHKERR DMMoFEMSetSquareProblem(subdm1, PETSC_FALSE);
    CHKERR DMMoFEMAddElement(subdm1, "FE01");
    CHKERR DMMoFEMAddSubFieldRow(subdm1, "FIELD0");
    CHKERR DMMoFEMAddSubFieldCol(subdm1, "FIELD1");
    CHKERR DMSetUp(subdm1);
    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>("SUB1", -1,
                                                                   -1, 1);

    if (debug) {
      Mat B;
      CHKERR DMCreateMatrix(subdm1, &B);
      MatView(B, PETSC_VIEWER_DRAW_WORLD);
      std::string wait;
      std::cin >> wait;
      CHKERR MatDestroy(&B);
    }

    // destry dm
    CHKERR DMDestroy(&dm);
    CHKERR DMDestroy(&subdm0);
    CHKERR DMDestroy(&subdm1);
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, ect.
  CHKERR MoFEM::Core::Finalize();

  return 0;
}
