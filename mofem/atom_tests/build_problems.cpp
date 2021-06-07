/** \file build_large_problem.cpp
  \example build_large_problem.cpp
  \brief Atom test for building problems

  \bug Not verifying what if partitioned mesh is loaded.

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

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

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

    CHKERR moab.load_file(mesh_file_name, 0, "");

    // Create MoFEM database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // set entities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);

    // Fields
    CHKERR m_field.add_field("F1", L2, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR m_field.add_field("F2", HDIV, AINSWORTH_LEGENDRE_BASE, 1);

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "F1");
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "F2");

    // set app. order
    // see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes
    // (Mark Ainsworth & Joe Coyle)
    int order = 2;
    CHKERR m_field.set_field_order(root_set, MBTET, "F1", order - 1);
    CHKERR m_field.set_field_order(root_set, MBTET, "F2", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "F2", order);

    CHKERR m_field.build_fields();

    // add elements
    CHKERR m_field.add_finite_element("E1");
    CHKERR m_field.add_finite_element("E2");
    CHKERR m_field.add_finite_element("E3");

    CHKERR m_field.modify_finite_element_add_field_row("E1", "F1");
    CHKERR m_field.modify_finite_element_add_field_col("E1", "F1");
    CHKERR m_field.modify_finite_element_add_field_row("E2", "F2");
    CHKERR m_field.modify_finite_element_add_field_col("E2", "F2");
    // To build composite problem
    CHKERR m_field.modify_finite_element_add_field_row("E3", "F1");
    CHKERR m_field.modify_finite_element_add_field_col("E3", "F2");

    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "E1");
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "E2");
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "E3");

    CHKERR m_field.build_finite_elements();
    CHKERR m_field.build_adjacencies(bit_level0);

    // Problems
    CHKERR m_field.add_problem("P1");
    CHKERR m_field.add_problem("P2");

    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("P1", bit_level0);
    CHKERR m_field.modify_problem_ref_level_add_bit("P2", bit_level0);

    CHKERR m_field.modify_problem_add_finite_element("P1", "E1");
    CHKERR m_field.modify_problem_add_finite_element("P2", "E2");

    // build problems
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("P1", true);
    CHKERR prb_mng_ptr->buildProblem("P2", true);
    CHKERR prb_mng_ptr->partitionProblem("P1");
    CHKERR prb_mng_ptr->partitionProblem("P2");

    CHKERR prb_mng_ptr->partitionFiniteElements("P1");
    CHKERR prb_mng_ptr->partitionGhostDofs("P1");
    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>("P1", -1, -1,
                                                                   0);
    CHKERR prb_mng_ptr->partitionFiniteElements("P2");
    CHKERR prb_mng_ptr->partitionGhostDofs("P2");
    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>("P2", -1, -1,
                                                                   0);

    // compose problem
    CHKERR m_field.add_problem("P3");
    CHKERR m_field.modify_problem_ref_level_add_bit("P3", bit_level0);
    CHKERR m_field.modify_problem_add_finite_element("P3", "E3");

    CHKERR prb_mng_ptr->buildProblem("P3", false);
    CHKERR prb_mng_ptr->inheritPartition("P3", "P1", false, "P2", true);
    CHKERR prb_mng_ptr->partitionFiniteElements("P3");
    CHKERR prb_mng_ptr->partitionGhostDofs("P3");
    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>("P3", -1, -1,
                                                                   0);
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  CHKERR MoFEM::Core::Finalize();

  return 0;
}
