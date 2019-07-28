/** \file dm_build_partitioned_mesh.cpp
  \example dm_partitioned_no_field.cpp
  \brief Testing problem for partitioned mesh with NOFIELD field

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

int main(int argc, char *argv[]) {
  // initialize petsc
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");

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
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
    const std::string options = "PARALLEL=READ_PART;"
                                "PARALLEL_RESOLVE_SHARED_ENTS;"
                                "PARTITION=PARALLEL_PARTITION;";
    CHKERR moab.load_file(mesh_file_name, 0, options.c_str());

    MoFEM::Core core(moab, PETSC_COMM_WORLD);
    MoFEM::Interface &m_field = core;

    EntityHandle root_set = moab.get_root_set();
    // add all entities to database, all of them will be used
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        root_set, 3, BitRefLevel().set(0));

    // add field
    CHKERR m_field.add_field("FIELD", H1, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR m_field.add_field("LAMBDA", NOFIELD, NOBASE, 3);

    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD");

    // set app. order
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD", 1);
    // Create vertices for NOFILE
    std::array<double, 6> coords = {0, 0, 0, 0, 0, 0};
    CHKERR m_field.create_vertices_and_add_to_field("LAMBDA", coords.data(), 2);
    CHKERR m_field.make_field_entities_multishared("LAMBDA", 0, NOISY);
    CHKERR m_field.getInterface<BitRefManager>()->setFieldEntitiesBitRefLevel(
        "LAMBDA", BitRefLevel().set(0));

    // build data structures for fields
    CHKERR m_field.build_fields();

    const FieldEntity_multiIndex *field_ents;
    CHKERR m_field.get_field_ents(&field_ents);
    for (auto it = field_ents->get<FieldName_mi_tag>().lower_bound("LAMBDA");
         it != field_ents->get<FieldName_mi_tag>().upper_bound("LAMBDA");
         ++it) {

      std::ostringstream ss;
      ss << "Rank " << m_field.get_comm_rank() << " -> " << **it << endl;
      PetscSynchronizedPrintf(m_field.get_comm(), "%s", ss.str().c_str());
      PetscSynchronizedFlush(m_field.get_comm(), PETSC_STDOUT);

    }

    // define & build finite elements
    CHKERR m_field.add_finite_element("FE");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("FE", "FIELD");
    CHKERR m_field.modify_finite_element_add_field_col("FE", "FIELD");
    CHKERR m_field.modify_finite_element_add_field_data("FE", "FIELD");
    CHKERR m_field.modify_finite_element_add_field_row("FE", "LAMBDA");
    CHKERR m_field.modify_finite_element_add_field_col("FE", "LAMBDA");
    CHKERR m_field.modify_finite_element_add_field_data("FE", "LAMBDA");

    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "FE");
    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(BitRefLevel().set());

    // set dm data structure which created mofem data structures
    CHKERR DMMoFEMCreateMoFEM(dm, &m_field, "TEST_PROBLEM",
                              BitRefLevel().set(0));
    CHKERR DMMoFEMSetSquareProblem(
        dm, PETSC_FALSE); // this is for testing (this problem has the same rows
                          // and cols)
    CHKERR DMSetFromOptions(dm);
    CHKERR DMMoFEMAddElement(dm, "FE");
    CHKERR DMSetUp(dm);

    Mat m;
    CHKERR DMCreateMatrix(dm, &m);

    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>(
            "TEST_PROBLEM", -1, -1, 1);

        CHKERR MatDestroy(&m);
    // destry dm
    CHKERR DMDestroy(&dm);
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, ect.
  CHKERR MoFEM::Core::Finalize();

  return 0;
}
