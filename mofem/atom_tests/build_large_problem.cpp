/** \file build_problems.cpp
  \example build_problems.cpp
  \brief atom test build problem
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

    if (flg != PETSC_TRUE)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "file not found -my_file %s",
               mesh_file_name);

    CHKERR moab.load_file(mesh_file_name, 0, "");

    // Create MoFEM database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // set entities bit level
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, BitRefLevel().set(0));

    // fields
    CHKERR m_field.add_field("F1", H1, AINSWORTH_LEGENDRE_BASE, 3,
                             MB_TAG_DENSE);

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBVERTEX, "F1");
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBEDGE, "F1");
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTRI, "F1");
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "F1");

    // set app. order
    constexpr int order = 3;
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "F1", 1);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "F1", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "F1", order);
    CHKERR m_field.set_field_order(root_set, MBTET, "F1", order);

    CHKERR m_field.build_fields();

    // add elements
    CHKERR m_field.add_finite_element("E1");
    CHKERR m_field.modify_finite_element_add_field_row("E1", "F1");
    CHKERR m_field.modify_finite_element_add_field_col("E1", "F1");
    CHKERR m_field.modify_finite_element_add_field_data("E1", "F1");
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "E1");
    CHKERR m_field.build_finite_elements();
    CHKERR m_field.build_adjacencies(BitRefLevel().set(0));

    // Problems
    CHKERR m_field.add_problem("P1");
    CHKERR m_field.modify_problem_ref_level_add_bit("P1", BitRefLevel().set(0));
    CHKERR m_field.modify_problem_add_finite_element("P1", "E1");

    // build problems
    auto prb_mng_ptr = m_field.getInterface<ProblemsManager>();
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("P1", true);
    CHKERR prb_mng_ptr->partitionProblem("P1");
    CHKERR prb_mng_ptr->partitionFiniteElements("P1");
    CHKERR prb_mng_ptr->partitionGhostDofs("P1");

    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_LOG("WORLD", Sev::inform) << "Create matrix";
    SmartPetscObj<Mat> A;
    CHKERR m_field.getInterface<MatrixManager>()
        ->createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>("P1", A);
    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_LOG("WORLD", Sev::inform) << "Done";

    auto fe_ptr = boost::make_shared<FEMethod>();

    auto pre_proc_hook = [&]() {
      MoFEMFunctionBegin;
      MoFEMFunctionReturn(0);
    };

    auto post_proc_hook = [&]() {
      MoFEMFunctionBegin;
      MoFEMFunctionReturn(0);
    };

    auto op_hook = [&]() {
      MoFEMFunctionBegin;

      MOFEM_LOG_CHANNEL("WORLD");
      MOFEM_LOG("WORLD", Sev::verbose)
          << "Iterate finite element " << *(fe_ptr->numeredEntFiniteElementPtr);

      auto row_dofs = fe_ptr->getRowDofsPtr();
      for (auto &dof : *row_dofs)
        MOFEM_LOG("WORLD", Sev::noisy) << *dof << endl;

      MoFEMFunctionReturn(0);
    };

    fe_ptr->preProcessHook = pre_proc_hook;
    fe_ptr->postProcessHook = post_proc_hook;
    fe_ptr->operatorHook = op_hook;

    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_LOG("WORLD", Sev::inform) << "Iterate finite elements";
    CHKERR m_field.loop_finite_elements("P1", "E1", *fe_ptr);
    MOFEM_LOG("WORLD", Sev::inform) << "Done";

    using Vol = VolumeElementForcesAndSourcesCore;
    using VolOp = VolumeElementForcesAndSourcesCore::UserDataOperator;
    using EntData = EntitiesFieldData::EntData;

    auto vol_ptr = boost::make_shared<Vol>(m_field);

    struct Op : public VolOp {

      Op() : VolOp("F1", VolOp::OPROWCOL, true) {}

      MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                            EntityType col_type, EntData &row_data,
                            EntData &col_data) {
        return 0;
      };
    };

    vol_ptr->getOpPtrVector().push_back(new Op());

    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_LOG("WORLD", Sev::inform) << "Iterate volume finite elements";
    CHKERR m_field.loop_finite_elements("P1", "E1", *vol_ptr);
    MOFEM_LOG("WORLD", Sev::inform) << "Done";

  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  CHKERR MoFEM::Core::Finalize();

  return 0;
}
