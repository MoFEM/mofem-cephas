/**
 * @file operators_tests.cpp
 * @example operators_tests.cpp
 * @brief Test operators in forms integrators
 * @date 2022-12-11
 *
 * @copyright Copyright (c) 2022
 *
 * TODO: Add more operators.
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
constexpr int SPACE_DIM = 3;
constexpr int FIELD_DIM = SPACE_DIM;

constexpr AssemblyType A =
    AssemblyType::BLOCK_MAT; //< selected assembly type
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
        LogManager::createSink(LogManager::getStrmWorld(), "OpTester"));
    LogManager::setLog("OpTester");
    MOFEM_LOG_TAG("OpTester", "OpTester");

    // Simple interface
    auto simple = m_field.getInterface<Simple>();

    // get options from command line
    CHKERR simple->getOptions();
    // load mesh file
    CHKERR simple->loadFile();

    // Scalar fields and vector field is tested. Add more fields, i.e. vector
    // field if needed.
    CHKERR simple->addDomainField("VECTOR", H1, AINSWORTH_LEGENDRE_BASE,
                                  FIELD_DIM);
    CHKERR simple->addDomainField("TENSOR", L2, AINSWORTH_LEGENDRE_BASE,
                                  FIELD_DIM);

    // set fields order, i.e. for most first cases order is sufficient.
    constexpr int order = 4;
    CHKERR simple->setFieldOrder("VECTOR", order);
    CHKERR simple->setFieldOrder("TENSOR", order);

    // setup problem
    CHKERR simple->setUp();

    auto schur_dm = createDM(m_field.get_comm(), "DMMOFEM");
    CHKERR DMMoFEMCreateSubDM(schur_dm, simple->getDM(), "SCHUR");
    CHKERR DMMoFEMSetSquareProblem(schur_dm, PETSC_TRUE);
    CHKERR DMMoFEMAddElement(schur_dm, simple->getDomainFEName());
    CHKERR DMMoFEMAddSubFieldRow(schur_dm, "VECTOR");
    CHKERR DMMoFEMAddSubFieldCol(schur_dm, "VECTOR");
    CHKERR DMSetUp(schur_dm);

    auto block_dm = createDM(m_field.get_comm(), "DMMOFEM");
    CHKERR DMMoFEMCreateSubDM(block_dm, simple->getDM(), "BLOCK");
    CHKERR DMMoFEMSetSquareProblem(block_dm, PETSC_TRUE);
    CHKERR DMMoFEMAddElement(block_dm, simple->getDomainFEName());
    CHKERR DMMoFEMAddSubFieldRow(block_dm, "TENSOR");
    CHKERR DMMoFEMAddSubFieldCol(block_dm, "TENSOR");
    CHKERR DMSetUp(block_dm);

    petsc_mat = createDMMatrix(simple->getDM());
    auto S = createDMMatrix(schur_dm);

    auto data_struture =
        createSchurBlockMatStructure(simple->getDM(),

                                     {{simple->getDomainFEName(),

                                       {{"VECTOR", "VECTOR"},
                                        {"TENSOR", "TENSOR"}

                                        ,

                                        {"VECTOR", "TENSOR"},
                                        {"TENSOR", "VECTOR"}

                                       }}}

        );
    auto shell_data = createSchurBlockMat(simple->getDM(), data_struture);
    auto [mat, data] = shell_data;
    block_mat = mat;

    using OpMassPETSCAssemble = FormsIntegrators<DomainEleOp>::Assembly<
        PETSC>::BiLinearForm<GAUSS>::OpMass<1, FIELD_DIM>;
    using OpMassBlockAssemble = FormsIntegrators<DomainEleOp>::Assembly<
        BLOCK_SCHUR>::BiLinearForm<GAUSS>::OpMass<1, FIELD_DIM>;

    // get operators tester
    auto pip_mng = m_field.getInterface<PipelineManager>(); // get interface to
                                                        // pipeline manager
    auto &pip_lhs = pip_mng->getOpDomainLhsPipeline();

    pip_lhs.push_back(new OpMassPETSCAssemble("VECTOR", "VECTOR"));
    pip_lhs.push_back(new OpMassPETSCAssemble("TENSOR", "TENSOR"));
    pip_lhs.push_back(new OpMassPETSCAssemble("VECTOR", "TENSOR"));
    pip_lhs.push_back(new OpMassPETSCAssemble("TENSOR", "VECTOR"));
    pip_lhs.push_back(createOpSchurAssembleBegin());
    pip_lhs.push_back(new OpMassBlockAssemble("VECTOR", "VECTOR"));
    pip_lhs.push_back(new OpMassBlockAssemble("TENSOR", "TENSOR"));
    pip_lhs.push_back(new OpMassBlockAssemble("VECTOR", "TENSOR"));
    pip_lhs.push_back(new OpMassBlockAssemble("TENSOR", "VECTOR"));

    auto schur_is = getDMSubData(schur_dm)->getSmartRowIs();
    auto ao_up = createAOMappingIS(schur_is, PETSC_NULL);

    pip_lhs.push_back(createOpSchurAssembleEnd(

        {"TENSOR"}, {nullptr}, {ao_up}, {S}, {false}, true)

    );

    CHKERR DMoFEMLoopFiniteElements(simple->getDM(), simple->getDomainFEName(),
                                    pip_mng->getDomainLhsFE());
    CHKERR MatAssemblyBegin(petsc_mat, MAT_FINAL_ASSEMBLY);
    CHKERR MatAssemblyEnd(petsc_mat, MAT_FINAL_ASSEMBLY);

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

    auto test = [](auto msg, auto y) {
      MoFEMFunctionBegin;
      PetscReal norm;
      CHKERR VecNorm(y, NORM_2, &norm);
      MOFEM_LOG("WORLD", Sev::inform)
          << msg << ": norm of difference: " << norm;
      constexpr double eps = 1e-10;
      if (norm > eps || std::isnan(norm) || std::isinf(norm)) {
        SETERRQ(PETSC_COMM_WORLD, 1, "norm of difference is too big");
      }
      MoFEMFunctionReturn(0);
    };

    std::vector<int> zero_rows_and_cols = {0, 1, 10, 20, 500, 1000};
    CHKERR MatZeroRowsColumns(petsc_mat, zero_rows_and_cols.size(),
                              &*zero_rows_and_cols.begin(), 1, PETSC_NULL,
                              PETSC_NULL);
    CHKERR MatZeroRowsColumns(block_mat, zero_rows_and_cols.size(),
                              &*zero_rows_and_cols.begin(), 1, PETSC_NULL,
                              PETSC_NULL);

    CHKERR MatMult(petsc_mat, v, y_petsc);
    CHKERR MatMult(block_mat, v, y_block);
    CHKERR VecAXPY(y_petsc, -1.0, y_block);

    CHKERR test("mult", y_petsc);

    CHKERR MatMult(petsc_mat, v, y_petsc);
    CHKERR MatMult(block_mat, v, y_block);
    CHKERR MatMultAdd(petsc_mat, v, y_petsc, y_petsc);
    CHKERR MatMultAdd(block_mat, v, y_block, y_block);
    CHKERR VecAXPY(y_petsc, -1.0, y_block);

    CHKERR test("mult add", y_petsc);

    auto [nested_mat, nested_data] = createSchurNestedMatrix(

        {schur_dm, block_dm},

        getSchurNestMatArray(

            {schur_dm, block_dm}, shell_data,

            {"TENSOR"}, {nullptr}

        )

    );

    auto y_nested = createDMVector(simple->getDM());
    CHKERR MatMult(petsc_mat, v, y_petsc);
    CHKERR MatMult(nested_mat, v, y_nested);
    CHKERR VecAXPY(y_petsc, -1.0, y_nested);

    CHKERR test("mult nested", y_petsc);

    auto diag_mat = nested_data.first[3];
    auto diag_block_x = get_random_vector(block_dm);
    auto diag_block_f = vectorDuplicate(diag_block_x);
    CHKERR MatMult(diag_mat, diag_block_x, diag_block_f);
    auto block_solved_x = vectorDuplicate(diag_block_x);
    CHKERR MatSolve(diag_mat, diag_block_f, block_solved_x);

    // auto ksp = createKSP(m_field.get_comm());
    // CHKERR KSPSetOperators(ksp, diag_mat, diag_mat);
    // CHKERR KSPSetFromOptions(ksp);

    // PC pc;
    // CHKERR KSPGetPC(ksp, &pc);
    // CHKERR PCSetType(pc, PCNONE);

    // CHKERR KSPSetFromOptions(ksp);
    // CHKERR KSPSolve(ksp, diag_block_f, block_solved_x);

    CHKERR VecAXPY(diag_block_x, -1.0, block_solved_x);
    CHKERR test("diag solve", diag_block_x);


    petsc_mat.reset();
    block_mat.reset();
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}