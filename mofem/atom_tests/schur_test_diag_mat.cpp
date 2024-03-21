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

    // set fields order, i.e. for most first cases order is sufficient.
    constexpr int order = 4;
    CHKERR simple->setFieldOrder("VECTOR", order);

    // setup problem
    CHKERR simple->setUp();

    petsc_mat = createDMMatrix(simple->getDM());

    auto data_struture = createSchurBlockMatStructure(
        simple->getDM(), {"VECTOR"}, {simple->getDomainFEName()},
        {boost::make_shared<DomainEle>(m_field)});
    auto [mat, data] = createSchurBlockMat(simple->getDM(), data_struture);
    block_mat = mat;

    using OpMassPETSCAssemble = FormsIntegrators<DomainEleOp>::Assembly<
        PETSC>::BiLinearForm<GAUSS>::OpMass<1, FIELD_DIM>;
    using OpMassBlockAssemble = FormsIntegrators<DomainEleOp>::Assembly<
        BLOCK_MAT>::BiLinearForm<GAUSS>::OpMass<1, FIELD_DIM>;

    // get operators tester
    auto pip_mng = m_field.getInterface<PipelineManager>(); // get interface to
                                                        // pipeline manager
    auto &pip_lhs = pip_mng->getOpDomainLhsPipeline();

    pip_lhs.push_back(new OpMassPETSCAssemble("VECTOR", "VECTOR"));
    pip_lhs.push_back(new OpMassBlockAssemble("VECTOR", "VECTOR"));
    CHKERR DMoFEMLoopFiniteElements(simple->getDM(), simple->getDomainFEName(),
                                    pip_mng->getDomainLhsFE());
    CHKERR MatAssemblyBegin(petsc_mat, MAT_FINAL_ASSEMBLY);
    CHKERR MatAssemblyEnd(petsc_mat, MAT_FINAL_ASSEMBLY);

    auto get_random_vector = [&]() {
      auto v = createDMVector(simple->getDM());
      PetscRandom rctx;
      PetscRandomCreate(PETSC_COMM_WORLD, &rctx);
      CHK_MOAB_THROW(VecSetRandom(v, rctx), "generate rand vector");
      PetscRandomDestroy(&rctx);
      return v;
    };
    auto v = get_random_vector();

    auto y_petsc = createDMVector(simple->getDM());
    auto y_block = createDMVector(simple->getDM());

    auto test = [](auto msg, auto y) {
      MoFEMFunctionBegin;
      PetscReal norm;
      CHKERR VecNorm(y, NORM_2, &norm);
      MOFEM_LOG("WORLD", Sev::inform)
          << msg << ": norm of difference: " << norm;
      constexpr double eps = 1e-12;
      if (norm > eps || std::isnan(norm) || std::isinf(norm)) {
        SETERRQ(PETSC_COMM_WORLD, 1, "norm of difference is too big");
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR MatMult(petsc_mat, v, y_petsc);
    CHKERR MatMult(block_mat, v, y_block);
    CHKERR VecAXPY(y_petsc, -1.0, y_block);

    CHKERR test("mult", y_petsc);

    CHKERR MatMultAdd(petsc_mat, v, y_petsc, y_petsc);
    CHKERR MatMultAdd(block_mat, v, y_block, y_block);
    CHKERR VecAXPY(y_petsc, -1.0, y_block);

    CHKERR test("mult", y_petsc);

    petsc_mat.reset();
    block_mat.reset();
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}