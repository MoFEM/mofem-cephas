/** \file build_problems.cpp
  \example build_problems.cpp
  \brief atom test build problem
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

    if (flg != PETSC_TRUE)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "file not found -my_file %s",
               mesh_file_name);

    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
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

    auto ref_ents = m_field.get_ref_ents();
    auto ref_fe_ents = m_field.get_ref_finite_elements();
    auto field_ents = m_field.get_field_ents();
    auto dofs = m_field.get_dofs();
    auto finite_elements = m_field.get_finite_elements();
    auto adjacencies = m_field.get_ents_elements_adjacency();

    // This realease data structures not used by the code. Once problem is
    // build, and you not plan create more problems, you can realease those
    // data structures. However, underlying data, on entities are not
    // realeased, since data can be sill acessed form finite element.
    const_cast<RefEntity_multiIndex *>(ref_ents)->clear();
    const_cast<RefElement_multiIndex *>(ref_fe_ents)->clear();
    const_cast<FieldEntity_multiIndex *>(field_ents)->clear();
    const_cast<DofEntity_multiIndex *>(dofs)->clear();
    const_cast<FiniteElement_multiIndex *>(finite_elements)->clear();

    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_LOG("WORLD", Sev::inform) << "Create matrix";
    SmartPetscObj<Mat> A;
    CHKERR m_field.getInterface<MatrixManager>()
        ->createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>("P1", A);
    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_LOG("WORLD", Sev::inform) << "Done";

    // Once matrix is created we do not need adjacencies data structures
    const_cast<FieldEntityEntFiniteElementAdjacencyMap_multiIndex *>(
        adjacencies)
        ->clear();

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

      auto row_dofs = fe_ptr->getRowDofs();
      for (auto &dof : row_dofs)
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

    auto vol_ptr =
        boost::make_shared<VolumeElementForcesAndSourcesCore>(m_field);
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
