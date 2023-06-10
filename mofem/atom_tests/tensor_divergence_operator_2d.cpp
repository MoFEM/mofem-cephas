/**
 * \file hcurl_divergence_operator_2d.cpp
 * \example hcurl_divergence_operator_2d.cpp
 *
 * Testing Hcurl base, transfromed to Hdiv base in 2d using Green theorem.
 *
 * Note 0: This is low-level implementation.
 * Note 1: Generic implementation for Quad/Tri mesh of arbitrary order.
 *
 */

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

constexpr AssemblyType A = AssemblyType::PETSC; //< selected assembly type
constexpr IntegrationType I =
    IntegrationType::GAUSS;                     //< selected integration type

constexpr int SPACE_DIM = EXECUTABLE_DIMENSION;
FTensor::Index<'i', SPACE_DIM> i;

template <int DIM> struct ElementsAndOps {};

template <> struct ElementsAndOps<2> : PipelineManager::ElementsAndOpsByDim<2> {
  static constexpr FieldSpace HDIV_SPACE = HCURL;
};

template <> struct ElementsAndOps<3> : PipelineManager::ElementsAndOpsByDim<3> {
  static constexpr FieldSpace HDIV_SPACE = HDIV;
};

constexpr FieldSpace ElementsAndOps<2>::HDIV_SPACE;
constexpr FieldSpace ElementsAndOps<3>::HDIV_SPACE;

using EntData = EntitiesFieldData::EntData;
using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainEleOp = DomainEle::UserDataOperator;
using BoundaryEle = ElementsAndOps<SPACE_DIM>::BoundaryEle;
using BoundaryEleOp = BoundaryEle::UserDataOperator;

using AssemblyDomainEleOp = FormsIntegrators<DomainEleOp>::Assembly<A>::OpBase;
using AssemblyBoundaryEleOp =
    FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;
constexpr FieldSpace HDIV_SPACE = ElementsAndOps<SPACE_DIM>::HDIV_SPACE;

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    // Create MoFEM instance
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // Declare elements
    enum bases { AINSWORTH, DEMKOWICZ, LASBASETOP };
    const char *list_bases[] = {"ainsworth", "demkowicz"};
    PetscBool flg;
    PetscInt choice_base_value = AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list_bases,
                                LASBASETOP, &choice_base_value, &flg);
    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "base not set");
    FieldApproximationBase base = AINSWORTH_LEGENDRE_BASE;
    if (choice_base_value == AINSWORTH)
      base = AINSWORTH_LEGENDRE_BASE;
    else if (choice_base_value == DEMKOWICZ)
      base = DEMKOWICZ_JACOBI_BASE;
    int order = 5;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order, PETSC_NULL);

    // Register DM Manager
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Add logging channel for example
    auto core_log = logging::core::get();
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "ATOM"));
    LogManager::setLog("ATOM");
    MOFEM_LOG_TAG("ATOM", "ATOM");

    // Simple interface
    auto simple = m_field.getInterface<Simple>();

    // get options from command line
    CHKERR simple->getOptions();
    // load mesh file
    CHKERR simple->loadFile();

    CHKERR simple->addDomainField("U", H1, base, SPACE_DIM);
    CHKERR simple->addDomainField("SIGMA", HDIV_SPACE, DEMKOWICZ_JACOBI_BASE,
                                  SPACE_DIM);
    CHKERR simple->addBoundaryField("U", H1, base, SPACE_DIM);
    CHKERR simple->addBoundaryField("SIGMA", HDIV_SPACE, DEMKOWICZ_JACOBI_BASE,
                                    SPACE_DIM);

    CHKERR simple->addDataField("GEOMETRY", H1, base, SPACE_DIM);

    // set fields order, i.e. for most first cases order is sufficient.
    CHKERR simple->setFieldOrder("U", order);
    CHKERR simple->setFieldOrder("SIGMA", order);
    CHKERR simple->setFieldOrder("GEOMETRY", 2);

    // setup problem
    CHKERR simple->setUp();

    auto project_ho_geometry = [&]() {
      Projection10NodeCoordsOnField ent_method(m_field, "GEOMETRY");
      return m_field.loop_dofs("GEOMETRY", ent_method);
    };
    CHKERR project_ho_geometry();

    // auto bc_mng = mField.getInterface<BcManager>();
    // CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "SYM",
    //                                          "U", 0, 0, true);
    // CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "SYM",
    //                                          "SIGMA", 0, 0, true);


    // get operators tester
    auto opt = m_field.getInterface<OperatorsTester>(); // get interface to
                                                        // OperatorsTester
    auto pip = m_field.getInterface<PipelineManager>(); // get interface to
                                                        // pipeline manager

    pip->getOpDomainLhsPipeline().clear();
    pip->getOpDomainRhsPipeline().clear();

    pip->setDomainLhsIntegrationRule([](int, int, int o) { return 2 * o; });
    pip->setDomainLhsIntegrationRule([](int, int, int o) { return 2 * o; });

    auto x = opt->setRandomFields(simple->getDM(),
                                  {{"U", {-1, 1}}, {"SIGMA", {-1, 1}}});
    CHKERR DMoFEMMeshToLocalVector(simple->getDM(), x, ADD_VALUES,
                                   SCATTER_REVERSE);

    using OpMixDivURhs = FormsIntegrators<DomainEleOp>::Assembly<A>::LinearForm<
        GAUSS>::OpMixDivTimesU<3, SPACE_DIM, SPACE_DIM>;
    using OpMixLambdaGradURhs = FormsIntegrators<DomainEleOp>::Assembly<
        PETSC>::LinearForm<I>::OpMixTensorTimesGradU<SPACE_DIM>;
    using OpMixLambdaGradURhs = FormsIntegrators<DomainEleOp>::Assembly<
        PETSC>::LinearForm<I>::OpMixTensorTimesGradU<SPACE_DIM>;
    using OpMixUTimesDivLambdaRhs = FormsIntegrators<DomainEleOp>::Assembly<
        PETSC>::LinearForm<I>::OpMixVecTimesDivLambda<SPACE_DIM>;

    using OpMixUTimesLambdaRhs = FormsIntegrators<DomainEleOp>::Assembly<
        PETSC>::LinearForm<I>::OpGradTimesTensor<1, SPACE_DIM, SPACE_DIM>;

    using OpMixNormalLambdaURhs = FormsIntegrators<BoundaryEleOp>::Assembly<
        PETSC>::LinearForm<I>::OpNormalMixVecTimesVectorField<SPACE_DIM>;
    using OpUTimeTractionRhs = FormsIntegrators<BoundaryEleOp>::Assembly<
        PETSC>::LinearForm<I>::OpBaseTimesVector<1, SPACE_DIM, 1>;

    auto pip_mng = m_field.getInterface<PipelineManager>();
    // integration rule
    auto rule = [&](int, int, int p) { return 2 * p + 1; };
    CHKERR pip_mng->setDomainRhsIntegrationRule(rule);
    CHKERR pip_mng->setBoundaryRhsIntegrationRule(rule);

    auto beta_domain = [](double x, double y, double z) { return 1.; };
    auto beta_bdy = [](double x, double y, double z) { return -1; };

    auto ops_rhs_interior = [&](auto &pip) {
      MoFEMFunctionBegin;
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1, HDIV},
                                                            "GEOMETRY");
      auto u_ptr = boost::make_shared<MatrixDouble>();
      auto grad_u_ptr = boost::make_shared<MatrixDouble>();
      pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));
      pip.push_back(new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>(
          "U", grad_u_ptr));

      pip.push_back(new OpMixDivURhs("SIGMA", u_ptr, beta_domain));
      pip.push_back(new OpMixLambdaGradURhs("SIGMA", grad_u_ptr, beta_domain));

      auto sigma_ptr = boost::make_shared<MatrixDouble>();
      auto sigma_div_ptr = boost::make_shared<MatrixDouble>();
      pip.push_back(new OpCalculateHVecTensorDivergence<SPACE_DIM, SPACE_DIM>(
          "SIGMA", sigma_div_ptr));
      pip.push_back(new OpCalculateHVecTensorField<SPACE_DIM, SPACE_DIM>(
          "SIGMA", sigma_ptr));
      pip.push_back(new OpMixUTimesDivLambdaRhs("U", sigma_div_ptr));
      pip.push_back(new OpMixUTimesLambdaRhs("U", sigma_ptr));

      MoFEMFunctionReturn(0);
    };

    auto ops_rhs_boundary = [&](auto &pip) {
      MoFEMFunctionBegin;
      CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(pip, {HDIV},
                                                                "GEOMETRY");
      auto u_ptr = boost::make_shared<MatrixDouble>();
      pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));
      auto traction_ptr = boost::make_shared<MatrixDouble>();
      pip.push_back(new OpCalculateHVecTensorTrace<SPACE_DIM, BoundaryEleOp>(
          "SIGMA", traction_ptr));

      // We have to integrate on curved face geometry, thus integration weight
      // have to adjusted.
      pip.push_back(new OpSetHOWeightsOnSubDim<SPACE_DIM>());
      pip.push_back(new OpMixNormalLambdaURhs("SIGMA", u_ptr, beta_bdy));
      pip.push_back(new OpUTimeTractionRhs("U", traction_ptr, beta_bdy));

      MoFEMFunctionReturn(0);
    };

    CHKERR ops_rhs_interior(pip_mng->getOpDomainRhsPipeline());
    CHKERR ops_rhs_boundary(pip_mng->getOpBoundaryRhsPipeline());

    auto f = createDMVector(simple->getDM());
    pip_mng->getDomainRhsFE()->f = f;
    pip_mng->getBoundaryRhsFE()->f = f;

    CHKERR VecZeroEntries(f);

    CHKERR DMoFEMLoopFiniteElements(simple->getDM(), simple->getDomainFEName(),
                                    pip_mng->getDomainRhsFE());
    CHKERR DMoFEMLoopFiniteElements(simple->getDM(),
                                    simple->getBoundaryFEName(),
                                    pip_mng->getBoundaryRhsFE());

    CHKERR VecAssemblyBegin(f);
    CHKERR VecAssemblyEnd(f);

    double f_nrm2;
    CHKERR VecNorm(f, NORM_2, &f_nrm2);

    MOFEM_LOG("ATOM", Sev::inform) << "f_norm2 = " << f_nrm2;
    if (std::fabs(f_nrm2) > 1e-10)
      SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "Test failed");

    // Testing Lhs domain operators

    using OpMixDivULhs = FormsIntegrators<DomainEleOp>::Assembly<
        PETSC>::BiLinearForm<I>::OpMixDivTimesVec<SPACE_DIM>;
    using OpLambdaGraULhs = FormsIntegrators<DomainEleOp>::Assembly<
        PETSC>::BiLinearForm<I>::OpMixTensorTimesGrad<SPACE_DIM>;

    auto op_lhs_domain = [&](auto &pip) {
      MoFEMFunctionBegin;
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1, HDIV},
                                                            "GEOMETRY");

      auto unity = []() { return 1; };
      pip.push_back(new OpMixDivULhs("SIGMA", "U", unity, true));
      pip.push_back(new OpLambdaGraULhs("SIGMA", "U", unity, true));
      MoFEMFunctionReturn(0);
    };

    CHKERR op_lhs_domain(pip_mng->getOpDomainLhsPipeline());

    auto diff_x = opt->setRandomFields(simple->getDM(),
                                       {{"U", {-1, 1}}, {"SIGMA", {-1, 1}}});
    constexpr double eps = 1e-6;
    auto diff_res = opt->checkCentralFiniteDifference(
        simple->getDM(), simple->getDomainFEName(), pip->getDomainRhsFE(),
        pip->getDomainLhsFE(), x, SmartPetscObj<Vec>(), SmartPetscObj<Vec>(),
        diff_x, 0, 1, eps);

    // Calculate norm of difference between directive calculated from finite
    // difference, and tangent matrix.
    double fnorm_res;
    CHKERR VecNorm(diff_res, NORM_2, &fnorm_res);
    MOFEM_LOG_C("ATOM", Sev::inform, "Test Lhs OPs %3.4e", fnorm_res);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
