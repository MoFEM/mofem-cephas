/**
 * @file schur_test_hdiv_contact.cpp
 * @example schur_test_hdiv_contact.cpp
 * @brief Test block matrix and Schur complement matrix for contact problem with
 * H(div) approximation of LM
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

template <int DIM> struct ElementsAndOps {};

template <> struct ElementsAndOps<2> {
  using DomainEle = PipelineManager::FaceEle;
};

template <> struct ElementsAndOps<3> {
  using DomainEle = VolumeElementForcesAndSourcesCore;
};

//! [Define dimension]
constexpr int SPACE_DIM = 2;
constexpr int FIELD_DIM = SPACE_DIM;

constexpr AssemblyType A = AssemblyType::BLOCK_MAT; //< selected assembly type
constexpr IntegrationType I =
    IntegrationType::GAUSS; //< selected integration type

using EntData = EntitiesFieldData::EntData;
using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainEleOp = DomainEle::UserDataOperator;

SmartPetscObj<Mat> petsc_mat;
SmartPetscObj<Mat> block_mat;

template <>
typename MoFEM::OpBaseImpl<PETSC, DomainEleOp>::MatSetValuesHook
    MoFEM::OpBaseImpl<PETSC, DomainEleOp>::matSetValuesHook =
        [](ForcesAndSourcesCore::UserDataOperator *op_ptr,
           const EntitiesFieldData::EntData &row_data,
           const EntitiesFieldData::EntData &col_data, MatrixDouble &m) {
          return MatSetValues<AssemblyTypeSelector<PETSC>>(
              petsc_mat, row_data, col_data, m, ADD_VALUES);
        };

template <>
typename MoFEM::OpBaseImpl<BLOCK_MAT, DomainEleOp>::MatSetValuesHook
    MoFEM::OpBaseImpl<BLOCK_MAT, DomainEleOp>::matSetValuesHook =
        [](ForcesAndSourcesCore::UserDataOperator *op_ptr,
           const EntitiesFieldData::EntData &row_data,
           const EntitiesFieldData::EntData &col_data, MatrixDouble &m) {
          return MatSetValues<AssemblyTypeSelector<BLOCK_MAT>>(
              block_mat, row_data, col_data, m, ADD_VALUES);
        };

template <>
typename MoFEM::OpBaseImpl<BLOCK_SCHUR, DomainEleOp>::MatSetValuesHook
    MoFEM::OpBaseImpl<BLOCK_SCHUR, DomainEleOp>::matSetValuesHook =
        [](ForcesAndSourcesCore::UserDataOperator *op_ptr,
           const EntitiesFieldData::EntData &row_data,
           const EntitiesFieldData::EntData &col_data, MatrixDouble &m) {
          return MatSetValues<AssemblyTypeSelector<BLOCK_SCHUR>>(
              block_mat, row_data, col_data, m, ADD_VALUES);
        };

template <>
typename MoFEM::OpBaseImpl<BLOCK_PRECONDITIONER_SCHUR,
                           DomainEleOp>::MatSetValuesHook
    MoFEM::OpBaseImpl<BLOCK_PRECONDITIONER_SCHUR,
                      DomainEleOp>::matSetValuesHook =
        [](ForcesAndSourcesCore::UserDataOperator *op_ptr,
           const EntitiesFieldData::EntData &row_data,
           const EntitiesFieldData::EntData &col_data, MatrixDouble &m) {
          return MatSetValues<AssemblyTypeSelector<BLOCK_PRECONDITIONER_SCHUR>>(
              block_mat, row_data, col_data, m, ADD_VALUES);
        };

constexpr bool debug = false;

int main(int argc, char *argv[]) {

  // initialize petsc
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    // Create MoAB database
    moab::Core moab_core;
    moab::Interface &moab = moab_core;

    // Create MoFEM database and link it to MoAB
    MoFEM::Core mofem_core(moab);
    MoFEM::Interface &m_field = mofem_core;

    // Register DM Manager
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Add logging channel for example
    auto core_log = logging::core::get();
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "Timeline"));
    LogManager::setLog("Timeline");
    MOFEM_LOG_TAG("Timeline", "Timeline");

    // Simple interface
    auto simple = m_field.getInterface<Simple>();

    // get options from command line
    CHKERR simple->getOptions();
    // load mesh file
    CHKERR simple->loadFile();

    // CHKERR simple->addDomainField("V", H1, AINSWORTH_LEGENDRE_BASE,
    // FIELD_DIM); CHKERR simple->addDomainField("T", L2,
    // AINSWORTH_LEGENDRE_BASE, FIELD_DIM); CHKERR simple->addDomainField("S",
    // L2, AINSWORTH_LEGENDRE_BASE, FIELD_DIM); CHKERR
    // simple->addDomainField("O", L2, AINSWORTH_LEGENDRE_BASE, FIELD_DIM);

    CHKERR simple->addDomainField("U", H1, AINSWORTH_LEGENDRE_BASE, SPACE_DIM);
    CHKERR simple->addBoundaryField("U", H1, AINSWORTH_LEGENDRE_BASE,
                                    SPACE_DIM);
    CHKERR simple->addDomainField("S", HCURL, DEMKOWICZ_JACOBI_BASE,
                                  SPACE_DIM);
    CHKERR simple->addBoundaryField("S", HCURL, DEMKOWICZ_JACOBI_BASE,
                                    SPACE_DIM);

    // set fields order, i.e. for most first cases order is sufficient.
    int order = 2;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order, PETSC_NULL);
    CHKERR simple->setFieldOrder("U", order);

    auto get_skin = [&]() {
      Range body_ents;
      CHKERR m_field.get_moab().get_entities_by_dimension(0, SPACE_DIM,
                                                          body_ents);
      Skinner skin(&m_field.get_moab());
      Range skin_ents;
      CHKERR skin.find_skin(0, body_ents, false, skin_ents);
      return skin_ents;
    };

    auto filter_blocks = [&](auto skin) {
      bool is_contact_block = false;
      Range contact_range;
      for (auto m : m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
               std::regex(

                   (boost::format("%s(.*)") % "CONTACT").str()

                       ))

      ) {
        is_contact_block = true; ///< blocs interation is collective, so that is
                                 ///< set irrespective if there are entities in
                                 ///< given rank or not in the block
        MOFEM_LOG("Timeline", Sev::inform)
            << "Find contact block set:  " << m->getName();
        auto meshset = m->getMeshset();
        Range contact_meshset_range;
        CHKERR m_field.get_moab().get_entities_by_dimension(
            meshset, SPACE_DIM - 1, contact_meshset_range, true);

        CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(
            contact_meshset_range);
        contact_range.merge(contact_meshset_range);
      }
      if (is_contact_block) {
        MOFEM_LOG("Timeline", Sev::inform)
            << "Nb entities in contact surface: " << contact_range.size();
        MOFEM_LOG_SYNCHRONISE(m_field.get_comm());
        skin = intersect(skin, contact_range);
      }
      return skin;
    };

    auto filter_true_skin = [&](auto skin) {
      Range boundary_ents;
      ParallelComm *pcomm =
          ParallelComm::get_pcomm(&m_field.get_moab(), MYPCOMM_INDEX);
      CHKERR pcomm->filter_pstatus(skin, PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                   PSTATUS_NOT, -1, &boundary_ents);
      return boundary_ents;
    };

    auto boundary_ents = filter_true_skin(filter_blocks(get_skin()));
    CHKERR simple->setFieldOrder("S", 0);
    CHKERR simple->setFieldOrder("S", order - 1, &boundary_ents);

    // setup problem
    CHKERR simple->setUp();

    auto schur_dm = createDM(m_field.get_comm(), "DMMOFEM");
    CHKERR DMMoFEMCreateSubDM(schur_dm, simple->getDM(), "SCHUR");
    CHKERR DMMoFEMSetSquareProblem(schur_dm, PETSC_TRUE);
    CHKERR DMMoFEMAddElement(schur_dm, simple->getDomainFEName());
    CHKERR DMMoFEMAddSubFieldRow(schur_dm, "U");
    CHKERR DMMoFEMAddSubFieldCol(schur_dm, "U");
    CHKERR DMSetUp(schur_dm);

    auto block_dm = createDM(m_field.get_comm(), "DMMOFEM");
    CHKERR DMMoFEMCreateSubDM(block_dm, simple->getDM(), "BLOCK");
    CHKERR DMMoFEMSetSquareProblem(block_dm, PETSC_TRUE);
    CHKERR DMMoFEMAddElement(block_dm, simple->getDomainFEName());
    CHKERR DMMoFEMAddSubFieldRow(block_dm, "S");
    CHKERR DMMoFEMAddSubFieldCol(block_dm, "S");
    // CHKERR DMMoFEMAddSubFieldRow(block_dm, "S");
    // CHKERR DMMoFEMAddSubFieldCol(block_dm, "S");
    // CHKERR DMMoFEMAddSubFieldRow(block_dm, "O");
    // CHKERR DMMoFEMAddSubFieldCol(block_dm, "O");
    CHKERR DMSetUp(block_dm);

    petsc_mat = createDMMatrix(simple->getDM());
    auto S = createDMMatrix(schur_dm);

    auto shell_data =
        createBlockMat(simple->getDM(),

                       createBlockMatStructure(simple->getDM(),

                                               {{simple->getDomainFEName(),

                                                 {{"U", "U"},
                                                  {"S", "S"},
                                                  {"U", "S"},
                                                  {"S", "U"}
                                                 }}}

                                               )

        );

    auto [mat, block_data_ptr] = shell_data;
    block_mat = mat;

    std::vector<std::string> fields{"S"};

    auto [nested_mat, nested_data_ptr] = createSchurNestedMatrix(

        getNestSchurData(

            {schur_dm, block_dm}, block_data_ptr,

            fields, {nullptr, nullptr, nullptr}, true

            )

    );

    using OpMassPETSCAssemble = FormsIntegrators<DomainEleOp>::Assembly<
        PETSC>::BiLinearForm<GAUSS>::OpMass<1, FIELD_DIM>;
    using OpMassBlockAssemble = FormsIntegrators<DomainEleOp>::Assembly<
        BLOCK_SCHUR>::BiLinearForm<GAUSS>::OpMass<1, FIELD_DIM>;
    using OpMassBlockPreconditionerAssemble =
        FormsIntegrators<DomainEleOp>::Assembly<BLOCK_PRECONDITIONER_SCHUR>::
            BiLinearForm<GAUSS>::OpMass<1, FIELD_DIM>;

    // get operators tester
    auto pip_mng = m_field.getInterface<PipelineManager>(); // get interface to
                                                            // pipeline manager

    auto close_zero = [](double, double, double) { return 1; };
    auto beta = [](double, double, double) { return -1. / 2; };
    auto gamma = [](double, double, double) { return -1. / 4; };

    pip_mng->setDomainLhsIntegrationRule([](int, int, int o) { return 2 * o; });
    auto &pip_lhs = pip_mng->getOpDomainLhsPipeline();

    pip_lhs.push_back(new OpMassPETSCAssemble("U", "U"));
    // pip_lhs.push_back(new OpMassPETSCAssemble("T", "T"));
    pip_lhs.push_back(new OpMassPETSCAssemble("S", "S"));
    pip_lhs.push_back(new OpMassPETSCAssemble("U", "S"));
    pip_lhs.push_back(new OpMassPETSCAssemble("S", "U"));
    // pip_lhs.push_back(new OpMassPETSCAssemble("S", "S", close_zero));
    // pip_lhs.push_back(new OpMassPETSCAssemble("S", "T", beta));
    // pip_lhs.push_back(new OpMassPETSCAssemble("T", "S", beta));
    // pip_lhs.push_back(new OpMassPETSCAssemble("O", "O", close_zero));
    // pip_lhs.push_back(new OpMassPETSCAssemble("T", "O", beta));
    // pip_lhs.push_back(new OpMassPETSCAssemble("O", "T", beta));
    // pip_lhs.push_back(new OpMassPETSCAssemble("S", "O", gamma));
    // pip_lhs.push_back(new OpMassPETSCAssemble("O", "S", gamma));

    pip_lhs.push_back(createOpSchurAssembleBegin());
    // pip_lhs.push_back(new OpMassBlockAssemble("V", "V"));
    // pip_lhs.push_back(new OpMassBlockAssemble("T", "T"));
    pip_lhs.push_back(new OpMassBlockAssemble("U", "U"));
    pip_lhs.push_back(new OpMassBlockAssemble("S", "S"));
    pip_lhs.push_back(new OpMassBlockAssemble("U", "S"));
    pip_lhs.push_back(new OpMassBlockAssemble("S", "U"));
    // pip_lhs.push_back(new OpMassBlockAssemble("S", "S", close_zero));
    // pip_lhs.push_back(new OpMassBlockAssemble("S", "T", beta));
    // pip_lhs.push_back(new OpMassBlockAssemble("T", "S", beta));
    // pip_lhs.push_back(new OpMassBlockAssemble("O", "O", close_zero));
    // pip_lhs.push_back(new OpMassBlockAssemble("T", "O", beta));
    // pip_lhs.push_back(new OpMassBlockAssemble("O", "T", beta));
    // pip_lhs.push_back(new OpMassBlockAssemble("S", "O", gamma));
    // pip_lhs.push_back(new OpMassBlockAssemble("O", "S", gamma));
    // pip_lhs.push_back(new OpMassBlockPreconditionerAssemble("T", "T"));

    auto schur_is = getDMSubData(schur_dm)->getSmartRowIs();
    auto ao_up = createAOMappingIS(schur_is, PETSC_NULL);

    pip_lhs.push_back(createOpSchurAssembleEnd(

        fields,

        {nullptr, nullptr, nullptr},

        {nullptr, nullptr, ao_up}, {nullptr, nullptr, S},

        {true, true, true}, true, block_data_ptr)

    );

    {
      MOFEM_LOG_CHANNEL("Timeline");
      MOFEM_LOG_TAG("Timeline", "timer");
      BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
      MOFEM_LOG("Timeline", Sev::inform) << "Assemble start";
      CHKERR DMoFEMLoopFiniteElements(simple->getDM(),
                                      simple->getDomainFEName(),
                                      pip_mng->getDomainLhsFE());
      MOFEM_LOG("Timeline", Sev::inform) << "Assemble end";
    }
    {
      MOFEM_LOG_CHANNEL("Timeline");
      MOFEM_LOG_TAG("Timeline", "timer");
      BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
      MOFEM_LOG("Timeline", Sev::inform) << "Mat assemble start";
      CHKERR MatAssemblyBegin(petsc_mat, MAT_FINAL_ASSEMBLY);
      CHKERR MatAssemblyEnd(petsc_mat, MAT_FINAL_ASSEMBLY);
      MOFEM_LOG("Timeline", Sev::inform) << "Mat assemble end";
    }

    auto get_random_vector = [&](auto dm) {
      auto v = createDMVector(dm);
      PetscRandom rctx;
      PetscRandomCreate(PETSC_COMM_WORLD, &rctx);
      CHK_MOAB_THROW(VecSetRandom(v, rctx), "generate rand vector");
      PetscRandomDestroy(&rctx);
      return v;
    };
    auto v = get_random_vector(simple->getDM());

    auto y_petsc = createDMVector(simple->getDM());
    auto y_block = createDMVector(simple->getDM());

    auto test = [](auto msg, auto y, double norm0) {
      MoFEMFunctionBegin;

      double eps = 1e-10;
      CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-eps", &eps, PETSC_NULL);

      PetscReal norm;
      CHKERR VecNorm(y, NORM_2, &norm);
      norm = norm / norm0;
      MOFEM_LOG_CHANNEL("WORLD");
      MOFEM_TAG_AND_LOG("WORLD", Sev::inform, "TestBlockMat")
          << msg << ": norm of difference: " << norm;
      if (norm > eps || std::isnan(norm) || std::isinf(norm)) {
        SETERRQ(PETSC_COMM_WORLD, 1, "norm of difference is too big");
      }
      MoFEMFunctionReturn(0);
    };

    std::vector<int> zero_rows_and_cols = {
        0, 1, 10, 20,
        500}; // not to remove dofs for TENSOR filed, inverse will not work
    CHKERR MatZeroRowsColumns(petsc_mat, zero_rows_and_cols.size(),
                              &*zero_rows_and_cols.begin(), 1, PETSC_NULL,
                              PETSC_NULL);
    CHKERR MatZeroRowsColumns(block_mat, zero_rows_and_cols.size(),
                              &*zero_rows_and_cols.begin(), 1, PETSC_NULL,
                              PETSC_NULL);

    {
      MOFEM_LOG_CHANNEL("Timeline");
      MOFEM_LOG_TAG("Timeline", "timer");
      BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
      MOFEM_LOG("Timeline", Sev::inform)
          << "MatMult(petsc_mat, v, y_petsc) star";
      CHKERR MatMult(petsc_mat, v, y_petsc);
      MOFEM_LOG("Timeline", Sev::inform)
          << "MatMult(petsc_mat, v, y_petsc) end";
    }
    double nrm0;
    CHKERR VecNorm(y_petsc, NORM_2, &nrm0);

    {
      MOFEM_LOG_CHANNEL("Timeline");
      MOFEM_LOG_TAG("Timeline", "timer");
      BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
      MOFEM_LOG("Timeline", Sev::inform)
          << "MatMult(block_mat, v, y_block) star";
      CHKERR MatMult(block_mat, v, y_block);
      MOFEM_LOG("Timeline", Sev::inform)
          << "MatMult(block_mat, v, y_block) end";
    }

    CHKERR VecAXPY(y_petsc, -1.0, y_block);
    CHKERR test("mult", y_petsc, nrm0);

    MOFEM_LOG_CHANNEL("Timeline");
    MOFEM_LOG_TAG("Timeline", "timer");

    CHKERR MatMult(petsc_mat, v, y_petsc);
    CHKERR MatMult(block_mat, v, y_block);
    CHKERR MatMultAdd(petsc_mat, v, y_petsc, y_petsc);
    CHKERR MatMultAdd(block_mat, v, y_block, y_block);
    CHKERR VecAXPY(y_petsc, -1.0, y_block);

    CHKERR test("mult add", y_petsc, nrm0);
    CHKERR schurSwitchPreconditioner(std::get<1>(*nested_data_ptr)[3]);
    auto y_nested = createDMVector(simple->getDM());
    CHKERR MatMult(petsc_mat, v, y_petsc);

    CHKERR MatMult(nested_mat, v, y_nested);
    CHKERR VecAXPY(y_petsc, -1.0, y_nested);

    CHKERR test("mult nested", y_petsc, nrm0);

    CHKERR schurSwitchPreconditioner(std::get<1>(*nested_data_ptr)[3]);
    auto diag_mat = std::get<0>(*nested_data_ptr)[3];
    auto diag_block_x = get_random_vector(block_dm);
    auto diag_block_f = createDMVector(block_dm);
    auto block_solved_x = createDMVector(block_dm);
    CHKERR MatMult(diag_mat, diag_block_x, diag_block_f);
    // That is if one like to use MatSolve directly, not though PC, as it is
    // below
    // CHKERR MatSolve(diag_mat, diag_block_f, block_solved_x);

    // set matrix type to shell, set data
    CHKERR DMSetMatType(block_dm, MATSHELL);
    CHKERR DMMoFEMSetBlocMatData(block_dm, std::get<1>(*nested_data_ptr)[3]);
    // set empty operator, since block data are already calculated
    CHKERR DMKSPSetComputeOperators(
        block_dm,
        [](KSP, Mat, Mat, void *) {
          MOFEM_LOG("WORLD", Sev::inform) << "empty operator";
          return 0;
        },
        nullptr);

    auto ksp = createKSP(m_field.get_comm());
    CHKERR KSPSetDM(ksp, block_dm);
    CHKERR KSPSetFromOptions(ksp);

    // set preconditioner to block mat
    auto get_pc = [](auto ksp) {
      PC pc_raw;
      CHKERR KSPGetPC(ksp, &pc_raw);
      return SmartPetscObj<PC>(pc_raw, true); // bump reference
    };
    CHKERR setSchurMatSolvePC(get_pc(ksp));
    CHKERR KSPSetUp(ksp);

    CHKERR VecZeroEntries(block_solved_x);
    CHKERR KSPSolve(ksp, diag_block_f, block_solved_x);

    auto diag_block_f_test = createDMVector(block_dm);
    CHKERR MatMult(diag_mat, block_solved_x, diag_block_f_test);
    CHKERR VecAXPY(diag_block_f_test, -1.0, diag_block_f);
    CHKERR test("diag solve", diag_block_f_test, nrm0);

    if (m_field.get_comm_rank() == 0) {
      CHKERR schurSaveBlockMesh(block_data_ptr, "block_mesh.vtk");
    }

    petsc_mat.reset();
    block_mat.reset();
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}