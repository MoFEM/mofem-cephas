/** \file remove_entities_from_problem.cpp
  \example remove_entities_from_problem.cpp
  \brief Remove field entities from the problem

*/

#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

constexpr int SPACE_DIM = 3;

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
    if (!flg)
      SETERRQ(PETSC_COMM_WORLD, 1, "*** ERROR -my_file (MESH FILE NEEDED)");
    
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
    const char *option;
    option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION="
             "PARALLEL_PARTITION;";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // Create MoFEM database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // set entities bit level
    const auto bit_level = BitRefLevel().set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level);

    // Fields
    CHKERR m_field.add_field("F1", HCURL, DEMKOWICZ_JACOBI_BASE, 3);
    CHKERR m_field.add_field("F2", HDIV, DEMKOWICZ_JACOBI_BASE, 1);

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_dim(root_set, SPACE_DIM, "F1");
    CHKERR m_field.add_ents_to_field_by_dim(root_set, SPACE_DIM, "F2");

    // set app. order
    // see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes
    // (Mark Ainsworth & Joe Coyle)
    int order = 2;
    for (EntityType t = CN::TypeDimensionMap[1].first;
         t <= CN::TypeDimensionMap[SPACE_DIM].second; ++t) {
      CHKERR m_field.set_field_order(root_set, t, "F1", order);
    }
    for (EntityType t = CN::TypeDimensionMap[2].first;
         t <= CN::TypeDimensionMap[SPACE_DIM].second; ++t) {
      CHKERR m_field.set_field_order(root_set, t, "F2", order);
    }

    CHKERR m_field.build_fields();

    // add elements
    CHKERR m_field.add_finite_element("E1");

    CHKERR m_field.modify_finite_element_add_field_row("E1", "F1");
    CHKERR m_field.modify_finite_element_add_field_row("E1", "F2");
    CHKERR m_field.modify_finite_element_add_field_col("E1", "F2");
    CHKERR m_field.modify_finite_element_add_field_col("E1", "F1");
    CHKERR m_field.modify_finite_element_add_field_data("E1", "F1");
    CHKERR m_field.modify_finite_element_add_field_data("E1", "F2");

    CHKERR m_field.add_ents_to_finite_element_by_dim(root_set, SPACE_DIM, "E1");

    CHKERR m_field.build_finite_elements();
    CHKERR m_field.build_adjacencies(bit_level);

    // Problems
    CHKERR m_field.add_problem("P1");

    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("P1", bit_level);
    CHKERR m_field.modify_problem_add_finite_element("P1", "E1");

    // build problems
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblemOnDistributedMesh("P1", true);
    CHKERR prb_mng_ptr->partitionFiniteElements("P1");
    CHKERR prb_mng_ptr->partitionGhostDofs("P1");

    auto get_triangles_on_skin = [&](Range &tets_skin) {
      MoFEMFunctionBegin;
      Range tets;
      CHKERR m_field.get_moab().get_entities_by_dimension(root_set, SPACE_DIM,
                                                          tets);
      Skinner skin(&m_field.get_moab());
      CHKERR skin.find_skin(0, tets, false, tets_skin);
      Range adj;
      CHKERR m_field.get_moab().get_adjacencies(tets_skin, SPACE_DIM - 2, false,
                                                adj, moab::Interface::UNION);
      tets_skin.merge(adj);
      MoFEMFunctionReturn(0);
    };

    Range tets;
    CHKERR m_field.get_moab().get_entities_by_dimension(root_set, SPACE_DIM,
                                                        tets);
    SmartPetscObj<IS> is_before_remove;
    CHKERR m_field.getInterface<ISManager>()->isCreateProblemFieldAndRank(
        "P1", ROW, "F1", 0, 1, is_before_remove, &tets);

    Range tets_skin;
    CHKERR get_triangles_on_skin(tets_skin);
    CHKERR prb_mng_ptr->removeDofsOnEntities("P1", "F1", tets_skin, 0, 1, 0,
                                             100, NOISY, true);

    SmartPetscObj<IS> is_after_remove;
    CHKERR m_field.getInterface<ISManager>()->isCreateProblemFieldAndRank(
        "P1", ROW, "F1", 0, 1, is_after_remove, &tets);

    auto test_is = [&](auto prb_name, auto is_before, auto is_after) {
      MoFEMFunctionBegin;
      const Problem *prb_ptr = m_field.get_problem(prb_name);
      if (auto sub_data = prb_ptr->getSubData()) {
        auto sub_ao = sub_data->getSmartRowMap();
        if (sub_ao) {
          CHKERR AOApplicationToPetscIS(sub_ao, is_before);
          PetscBool is_the_same;
          CHKERR ISEqual(is_before, is_after, &is_the_same);
          if (is_the_same == PETSC_FALSE) {
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                    "IS should be the same if map is correctly implemented");
          } else {
            MOFEM_LOG("WORLD", Sev::inform) << "Sub data map is correct";
          }
        } else {
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "AO map should exist");
        }
      } else {
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "Sub DM should exist");
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR test_is("P1", is_before_remove, is_after_remove);

    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>("P1", -2, -2,
                                                                   0);

    CHKERR m_field.add_problem("SUB");
     // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("SUB", bit_level);
    CHKERR m_field.modify_problem_add_finite_element("SUB", "E1");
    CHKERR prb_mng_ptr->buildSubProblem("SUB", {"F1"}, {"F1"}, "P1", PETSC_TRUE,
                                        nullptr, nullptr);
    CHKERR prb_mng_ptr->partitionFiniteElements("SUB", true, 0,
                                                m_field.get_comm_size());
    CHKERR prb_mng_ptr->partitionGhostDofsOnDistributedMesh("SUB");

    CHKERR prb_mng_ptr->removeDofsOnEntities("SUB", "F1", tets_skin, 0, 1, 0,
                                             100, NOISY, true);
    SmartPetscObj<IS> is_sub_after_remove;
    CHKERR m_field.getInterface<ISManager>()->isCreateProblemFieldAndRank(
        "SUB", ROW, "F1", 0, MAX_DOFS_ON_ENTITY, is_sub_after_remove, &tets);
    CHKERR test_is("SUB", is_before_remove, is_after_remove);

    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>("SUB", -2,
                                                                   -2, 0);
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  CHKERR MoFEM::Core::Finalize();

  return 0;
}
