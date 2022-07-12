/**
 * \file serial_matrix.cpp
 *
 * \brief testing create serial/sequential matrix
 */

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

namespace bio = boost::iostreams;
using bio::stream;
using bio::tee_device;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, PETSC_NULL, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

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
    auto moab_comm_wrap =
        boost::make_shared<WrapMPIComm>(PETSC_COMM_WORLD, false);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, moab_comm_wrap->get_comm());

    int do_for_rank = 0;
    if (pcomm->rank() ==
        (unsigned int)do_for_rank) { // should work only with rank 0

      // Create MoFEM (Joseph) database
      // second argument set communicator for sequential problem
      // last argument make mofem QUIET
      MoFEM::Core core(moab, PETSC_COMM_SELF, -1);
      MoFEM::Interface &m_field = core;

      // stl::bitset see for more details
      BitRefLevel bit_level0;
      bit_level0.set(0);
      EntityHandle meshset_level0;
      CHKERR moab.create_meshset(MESHSET_SET, meshset_level0);
      CHKERRG(rval);
      CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
          0, 3, bit_level0);
      CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByRefLevel(
          bit_level0, BitRefLevel().set(), meshset_level0);

      /***/
      // Define problem

      // Fields
      CHKERR m_field.add_field("FIELD_A", H1, AINSWORTH_LEGENDRE_BASE, 3,
                               MB_TAG_DENSE);
      CHKERR m_field.add_field("FIELD_B", L2, AINSWORTH_LEGENDRE_BASE, 1,
                               MB_TAG_DENSE);

      CHKERR m_field.add_ents_to_field_by_type(0, MBTET, "FIELD_A");
      CHKERR m_field.add_ents_to_field_by_type(0, MBTET, "FIELD_B");

      CHKERR m_field.set_field_order(0, MBTET, "FIELD_A", order);
      CHKERR m_field.set_field_order(0, MBTRI, "FIELD_A", order);
      CHKERR m_field.set_field_order(0, MBEDGE, "FIELD_A", order);
      CHKERR m_field.set_field_order(0, MBVERTEX, "FIELD_A", 1);

      CHKERR m_field.set_field_order(0, MBTET, "FIELD_B", order - 1);
      CHKERR m_field.set_field_order(0, MBTRI, "FIELD_B", order - 1);
      CHKERR m_field.set_field_order(0, MBEDGE, "FIELD_B", order - 1);

      // build field
      CHKERR m_field.build_fields();

      // Element
      CHKERR m_field.add_finite_element("FE1");
      CHKERR m_field.add_finite_element("FE2");

      // Define rows/cols and element data
      CHKERR m_field.modify_finite_element_add_field_row("FE1", "FIELD_A");
      CHKERR m_field.modify_finite_element_add_field_col("FE1", "FIELD_A");
      CHKERR m_field.modify_finite_element_add_field_data("FE1", "FIELD_A");
      CHKERR m_field.modify_finite_element_add_field_col("FE1", "FIELD_B");
      CHKERR m_field.modify_finite_element_add_field_data("FE1", "FIELD_B");

      // Define rows/cols and element data
      CHKERR m_field.modify_finite_element_add_field_row("FE2", "FIELD_B");
      CHKERR m_field.modify_finite_element_add_field_col("FE2", "FIELD_A");
      CHKERR m_field.modify_finite_element_add_field_data("FE2", "FIELD_B");

      // add entities to finite element
      CHKERR m_field.add_ents_to_finite_element_by_type(0, MBTET, "FE1");
      CHKERR m_field.add_ents_to_finite_element_by_type(0, MBTET, "FE2");

      // build finite elemnts
      CHKERR m_field.build_finite_elements();
      // build adjacencies
      CHKERR m_field.build_adjacencies(bit_level0);

      // Problem
      CHKERR m_field.add_problem("TEST_PROBLEM");

      // set refinement level for problem
      CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",
                                                      bit_level0);
      CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "FE1");
      CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "FE2");

      ProblemsManager *prb_mng_ptr;
      CHKERR m_field.getInterface(prb_mng_ptr);
      CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);
      CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");

      // partition finite elements
      CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
      // what are ghost nodes, see Petsc Manual
      CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

      Vec F;
      CHKERR m_field.getInterface<VecManager>()->vecCreateSeq("TEST_PROBLEM",
                                                              ROW, &F);
      Mat A;
      CHKERR m_field.getInterface<MatrixManager>()
          ->createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>("TEST_PROBLEM", &A);

      PetscViewer viewer;
      CHKERR PetscViewerASCIIOpen(PETSC_COMM_SELF, "serial_matrix.txt",
                                  &viewer);
      CHKERR PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_INFO);
      CHKERR MatView(A, viewer);

      CHKERR MatDestroy(&A);
      CHKERR VecDestroy(&F);
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  return 0;
}
