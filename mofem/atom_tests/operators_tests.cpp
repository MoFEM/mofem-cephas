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
constexpr AssemblyType A = AssemblyType::PETSC; //< selected assembly type
constexpr IntegrationType I =
    IntegrationType::GAUSS; //< selected integration type

using EntData = EntitiesFieldData::EntData;
using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainEleOp = DomainEle::UserDataOperator;
using PostProcEle = PostProcBrokenMeshInMoab<DomainEle>;


// Specializations for tested operators

template <int FIELD_DIM>
using OpGradGrad = FormsIntegrators<DomainEleOp>::Assembly<A>::BiLinearForm<
    I>::OpGradGrad<1, FIELD_DIM, SPACE_DIM>;

template <int FIELD_DIM>
using OpGradTimesTensor = FormsIntegrators<DomainEleOp>::Assembly<
    A>::LinearForm<I>::OpGradTimesTensor<1, FIELD_DIM, SPACE_DIM>;

template <int FIELD_DIM>
using OpConvectiveTermRhs = FormsIntegrators<DomainEleOp>::Assembly<
    A>::LinearForm<I>::OpConvectiveTermRhs<1, FIELD_DIM, SPACE_DIM>;

template <int FIELD_DIM>
using OpConvectiveTermLhsDu = FormsIntegrators<DomainEleOp>::Assembly<
    A>::BiLinearForm<I>::OpConvectiveTermLhsDu<1, FIELD_DIM, SPACE_DIM>;

template <int FIELD_DIM>
using OpConvectiveTermLhsDy = FormsIntegrators<DomainEleOp>::Assembly<
    A>::BiLinearForm<I>::OpConvectiveTermLhsDy<1, FIELD_DIM, SPACE_DIM>;

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
    CHKERR simple->addDomainField("SCALAR", H1, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR simple->addDomainField("VECTOR", H1, AINSWORTH_LEGENDRE_BASE,
                                  SPACE_DIM);

    // set fields order, i.e. for most first cases order is sufficient.
    CHKERR simple->setFieldOrder("SCALAR", 1);
    CHKERR simple->setFieldOrder("VECTOR", 1);

    // setup problem
    CHKERR simple->setUp();

    // get operators tester
    auto opt = m_field.getInterface<OperatorsTester>(); // get interface to
                                                        // OperatorsTester
    auto pip = m_field.getInterface<PipelineManager>(); // get interface to
                                                        // pipeline manager

    auto post_proc = [&](auto dm, auto f_res, auto out_name) {
      MoFEMFunctionBegin;
      auto post_proc_fe =
          boost::make_shared<PostProcBrokenMeshInMoab<DomainEle>>(m_field);

      using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

      auto scal_vec = boost::make_shared<VectorDouble>();
      auto vec_mat = boost::make_shared<MatrixDouble>();

      post_proc_fe->getOpPtrVector().push_back(
          new OpCalculateScalarFieldValues("SCALAR", scal_vec, f_res));
      post_proc_fe->getOpPtrVector().push_back(
          new OpCalculateVectorFieldValues<SPACE_DIM>("VECTOR", vec_mat,
                                                      f_res));

      post_proc_fe->getOpPtrVector().push_back(

          new OpPPMap(

              post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

              {{"SCALAR", scal_vec}},

              {{"VECTOR", vec_mat}},

              {}, {})

      );

      CHKERR DMoFEMLoopFiniteElements(dm, simple->getDomainFEName(),
                                      post_proc_fe);
      post_proc_fe->writeFile(out_name);
      MoFEMFunctionReturn(0);
    };

    // Note: I do not push OPs transferring base to physical coordinates. That
    // is not needed to test OPs, and not required, since we like to test only
    // operators.

    // Test grad-grad operator for scalar and vector field specialisation
    auto TestOpGradGrad = [&]() {
      MoFEMFunctionBegin;
      MOFEM_LOG("OpTester", Sev::verbose) << "TestOpGradGrad";

      pip->getOpDomainLhsPipeline().clear();
      pip->getOpDomainRhsPipeline().clear();

      pip->setDomainLhsIntegrationRule([](int, int, int o) { return 2 * o; });
      pip->setDomainLhsIntegrationRule([](int, int, int o) { return 2 * o; });

      pip->getOpDomainLhsPipeline().push_back(new OpGradGrad<1>(
          "SCALAR", "SCALAR", [](double, double, double) { return 1; }));
      pip->getOpDomainLhsPipeline().push_back(new OpGradGrad<SPACE_DIM>(
          "VECTOR", "VECTOR", [](double, double, double) { return 1; }));

      auto scl_mat = boost::make_shared<MatrixDouble>();
      auto vec_mat = boost::make_shared<MatrixDouble>();

      pip->getOpDomainRhsPipeline().push_back(
          new OpCalculateScalarFieldGradient<SPACE_DIM>("SCALAR", scl_mat));
      pip->getOpDomainRhsPipeline().push_back(
          new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>("VECTOR",
                                                                   vec_mat));

      pip->getOpDomainRhsPipeline().push_back(new OpGradTimesTensor<1>(
          "SCALAR", scl_mat, [](double, double, double) { return 1; }));
      pip->getOpDomainRhsPipeline().push_back(new OpGradTimesTensor<SPACE_DIM>(
          "VECTOR", vec_mat, [](double, double, double) { return 1; }));

      constexpr double eps = 1e-6;

      auto x = opt->setRandomFields(simple->getDM(),
                                    {{"SCALAR", {-1, 1}}, {"VECTOR", {-1, 1}}});
      auto diff_x = opt->setRandomFields(
          simple->getDM(), {{"SCALAR", {-1, 1}}, {"VECTOR", {-1, 1}}});

      auto diff_res = opt->checkCentralFiniteDifference(
          simple->getDM(), simple->getDomainFEName(), pip->getDomainRhsFE(),
          pip->getDomainLhsFE(), x, SmartPetscObj<Vec>(), SmartPetscObj<Vec>(),
          diff_x, 0, 1, eps);

      // Calculate norm of difference between directive calculated from finite
      // difference, and tangent matrix.
      double fnorm;
      CHKERR VecNorm(diff_res, NORM_2, &fnorm);
      MOFEM_LOG_C("OpTester", Sev::inform, "TestOpGradGrad %3.4e", fnorm);

      constexpr double err = 1e-9;
      if (fnorm > err)
        SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                "Norm of directional derivative too large");

      if (debug) {
        // Example how to plot direction in direction diff_x. If instead
        // directionalCentralFiniteDifference(...) diff_res is used, then error
        // on directive is plotted.
        CHKERR post_proc(simple->getDM(),
                         //
                         opt->directionalCentralFiniteDifference(
                             simple->getDM(), simple->getDomainFEName(),
                             pip->getDomainRhsFE(), x, SmartPetscObj<Vec>(),
                             SmartPetscObj<Vec>(), diff_x, 0, 1, eps),
                         //
                         "out_directional_directive_gradgrad.h5m");
      }

      MoFEMFunctionReturn(0);
    };

    // test operator for convection (part of material time directives)
    auto TestOpConvection = [&]() {
      MoFEMFunctionBegin;
      MOFEM_LOG("OpTester", Sev::verbose) << "TestOpConvection";

      pip->getOpDomainLhsPipeline().clear();
      pip->getOpDomainRhsPipeline().clear();

      pip->setDomainLhsIntegrationRule([](int, int, int o) { return 2 * o; });
      pip->setDomainLhsIntegrationRule([](int, int, int o) { return 2 * o; });

      auto evaluate_field = [&](auto &pipe) {
        auto scl_mat = boost::make_shared<MatrixDouble>();
        auto vec_mat = boost::make_shared<MatrixDouble>();
        auto dot_vec_vel = boost::make_shared<MatrixDouble>();
        pipe.push_back(
            new OpCalculateScalarFieldGradient<SPACE_DIM>("SCALAR", scl_mat));
        pipe.push_back(new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>(
            "VECTOR", vec_mat));
        pipe.push_back(new OpCalculateVectorFieldValuesDot<SPACE_DIM>(
            "VECTOR", dot_vec_vel));
        return std::make_tuple(scl_mat, vec_mat, dot_vec_vel);
      };

      auto [rhs_scl_mat, rhs_vec_mat, rhs_dot_vec_vel] =
          evaluate_field(pip->getOpDomainRhsPipeline());
      pip->getOpDomainRhsPipeline().push_back(
          new OpConvectiveTermRhs<1>("SCALAR", rhs_dot_vec_vel, rhs_scl_mat));
      pip->getOpDomainRhsPipeline().push_back(
          new OpConvectiveTermRhs<SPACE_DIM>("VECTOR", rhs_dot_vec_vel,
                                             rhs_vec_mat));

      auto [lhs_scl_mat, lhs_vec_mat, lhs_dot_vec_vel] =
          evaluate_field(pip->getOpDomainLhsPipeline());

      // I create op, set scaling function to calculate time directive, and add
      // operator pinter to pipeline
      auto op_convective_term_lhs_du_scalar =
          new OpConvectiveTermLhsDu<1>("SCALAR", "VECTOR", lhs_scl_mat);
      op_convective_term_lhs_du_scalar->feScalingFun =
          [](const FEMethod *fe_ptr) { return fe_ptr->ts_a; };
      pip->getOpDomainLhsPipeline().push_back(op_convective_term_lhs_du_scalar);
      pip->getOpDomainLhsPipeline().push_back(
          new OpConvectiveTermLhsDy<1>("SCALAR", "SCALAR", lhs_dot_vec_vel));

      auto op_convective_term_lhs_du_vec =
          new OpConvectiveTermLhsDu<SPACE_DIM>("VECTOR", "VECTOR", lhs_vec_mat);
      op_convective_term_lhs_du_vec->feScalingFun = [](const FEMethod *fe_ptr) {
        return fe_ptr->ts_a;
      };
      pip->getOpDomainLhsPipeline().push_back(op_convective_term_lhs_du_vec);
      pip->getOpDomainLhsPipeline().push_back(
          new OpConvectiveTermLhsDy<SPACE_DIM>("VECTOR", "VECTOR",
                                               lhs_dot_vec_vel));

      constexpr double eps = 1e-6;

      auto x = opt->setRandomFields(simple->getDM(),
                                    {{"SCALAR", {-1, 1}}, {"VECTOR", {-1, 1}}});
      auto x_t = opt->setRandomFields(
          simple->getDM(), {{"SCALAR", {-1, 1}}, {"VECTOR", {-1, 1}}});
      auto diff_x = opt->setRandomFields(
          simple->getDM(), {{"SCALAR", {-1, 1}}, {"VECTOR", {-1, 1}}});

      auto diff_res = opt->checkCentralFiniteDifference(
          simple->getDM(), simple->getDomainFEName(), pip->getDomainRhsFE(),
          pip->getDomainLhsFE(), x, x_t, SmartPetscObj<Vec>(), diff_x, 0, 1,
          eps);

      // Calculate norm of difference between directive calculated from finite
      // difference, and tangent matrix.
      double fnorm;
      CHKERR VecNorm(diff_res, NORM_2, &fnorm);
      MOFEM_LOG_C("OpTester", Sev::inform, "TestOpGradGrad %3.4e", fnorm);

      constexpr double err = 1e-9;
      if (fnorm > err) {
        SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                "Norm of directional derivative too large");
      }

      if (debug) {
        // Example how to plot direction in direction diff_x. If instead
        // directionalCentralFiniteDifference(...) diff_res is used, then error
        // on directive is plotted.
        CHKERR post_proc(simple->getDM(),
                         //
                         opt->directionalCentralFiniteDifference(
                             simple->getDM(), simple->getDomainFEName(),
                             pip->getDomainRhsFE(), x, x_t,
                             SmartPetscObj<Vec>(), diff_x, 0, 1, eps),
                         //
                         "out_directional_directive_convection.h5m");
      }

      MoFEMFunctionReturn(0);
    };

    CHKERR TestOpGradGrad();
    CHKERR TestOpConvection();
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}