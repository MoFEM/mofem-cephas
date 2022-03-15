/** \file dm_build_partitioned_mesh.cpp
  \example dm_build_partitioned_mesh.cpp
  \brief Testing problem for partitioned mesh

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

struct CountUp : FEMethod {

  CountUp(int &counter) : FEMethod(), cOunter(counter) {}
  MoFEMErrorCode preProcess() { return 0; }
  MoFEMErrorCode operator()() {
    ++cOunter;
    return 0;
  }
  MoFEMErrorCode postProcess() { return 0; }

private:
  int &cOunter;

};

struct CountDown : FEMethod {

  CountDown(int &counter) : FEMethod(), cOunter(counter) {}
    MoFEMErrorCode preProcess() { return 0; }
  MoFEMErrorCode operator()() {
    --cOunter;
    return 0;
  }
  MoFEMErrorCode postProcess() { return 0; }

private:
  int &cOunter;

};

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

    // create dm instance
    DM dm;
    CHKERR DMCreate(PETSC_COMM_WORLD, &dm);
    CHKERR DMSetType(dm, dm_name);

    // read mesh and create moab and mofem data structures
    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    auto moab_comm_wrap =
        boost::make_shared<WrapMPIComm>(PETSC_COMM_WORLD, false);
    if (pcomm == NULL)
      pcomm =
          new ParallelComm(&moab, moab_comm_wrap->get_comm());

    const char *option;
    option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION="
             "PARALLEL_PARTITION;";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    MoFEM::Core core(moab, PETSC_COMM_WORLD);
    MoFEM::Interface &m_field = core;

    EntityHandle root_set = moab.get_root_set();
    // add all entities to database, all of them will be used
    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        root_set, 3, bit_level0);
    // define & build field
    int field_rank = 3;
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-my_field_rank", &field_rank,
                              &flg);
#else
    CHKERR PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-my_field_rank",
                              &field_rank, &flg);
#endif
    CHKERR m_field.add_field("FIELD1", H1, AINSWORTH_LEGENDRE_BASE, field_rank);
    CHKERR m_field.add_field("FIELD2", L2, AINSWORTH_LEGENDRE_BASE, 1);

    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD1");
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD2");

    // set app. order
    int order = 4;
    CHKERR m_field.set_field_order(root_set, MBTET, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD1", 1);
    CHKERR m_field.set_field_order(root_set, MBTET, "FIELD2", order);

    // build data structures for fields
    CHKERR m_field.build_fields();

    // define & build finite elements
    CHKERR m_field.add_finite_element("FE");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("FE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("FE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("FE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_row("FE", "FIELD2");
    CHKERR m_field.modify_finite_element_add_field_col("FE", "FIELD2");
    CHKERR m_field.modify_finite_element_add_field_data("FE", "FIELD2");
    // Only data
    CHKERR m_field.add_finite_element("FE_ONLY_DATA");
    CHKERR m_field.modify_finite_element_add_field_data("FE_ONLY_DATA",
                                                        "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("FE_ONLY_DATA",
                                                        "FIELD2");
    // Add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "FE");
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET,
                                                      "FE_ONLY_DATA");
    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);

    // set dm data structure which created mofem data structures
    CHKERR DMMoFEMCreateMoFEM(dm, &m_field, "TEST_PROBLEM", bit_level0);
    CHKERR DMMoFEMSetIsPartitioned(dm, PETSC_TRUE);
    CHKERR DMMoFEMSetSquareProblem(
        dm, PETSC_FALSE); // this is for testing (this problem has the same rows
                          // and cols)
    CHKERR DMSetFromOptions(dm);
    CHKERR DMMoFEMAddElement(dm, "FE");
    CHKERR DMMoFEMAddElement(dm, "FE_ONLY_DATA");
    CHKERR DMSetUp(dm);

    // dump data to file, just to check if something was changed
    SmartPetscObj<Mat> m;
    CHKERR DMCreateMatrix_MoFEM(dm, m);

    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>(
            "TEST_PROBLEM", -1, -1, 1);

    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJMatrixFillIn<PetscGlobalIdx_mi_tag>("TEST_PROBLEM", -1, -1,
                                                         1);

    std::vector<std::string> fields_list;
    fields_list.push_back("FIELD1");

    // PetscSection section;
    PetscSection section;
    CHKERR m_field.getInterface<ISManager>()->sectionCreate("TEST_PROBLEM",
                                                            &section);
    CHKERR PetscSectionView(section, PETSC_VIEWER_STDOUT_WORLD);
    CHKERR DMSetSection(dm, section);
    CHKERR PetscSectionDestroy(&section);

    PetscBool save_file = PETSC_TRUE;
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-my_save_file", &save_file,
                               &flg);
#else
    CHKERR PetscOptionsGetBool(PETSC_NULL, PETSC_NULL, "-my_save_file",
                               &save_file, &flg);
#endif
    if (save_file) {
      PetscViewer viewer;
      CHKERR PetscViewerASCIIOpen(PETSC_COMM_WORLD,
                                  "dm_build_partitioned_mesh.txt", &viewer);
      CHKERR PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_INFO);
      MatView(m, viewer);
      CHKERR PetscViewerDestroy(&viewer);
    }

    int count = 0;
    CHKERR DMoFEMLoopFiniteElements(dm, "FE",
                                    boost::make_shared<CountUp>(count));
    CHKERR DMoFEMLoopFiniteElements(dm, "FE_ONLY_DATA",
                                    boost::make_shared<CountDown>(count));
    if(count)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Should be zero %d",
               count);

    // destry dm
    CHKERR DMDestroy(&dm);

  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, ect.
  CHKERR MoFEM::Core::Finalize();

  return 0;
}
