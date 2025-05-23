/**
 * \file tensor_divergence_operator.cpp
 * \example tensor_divergence_operator.cpp
 * 
 * Unit test for:
 * 1. Integration linear forms, and consistency of boundary integrals
 * 2. Integration integration on high-order geometry (2d and 3d)
 * 3. Integration for axi-symmetric case
 * 4. Consistency of Lhs operators, and Rhs operators (using OperatorsTester)
 * 
 * Note: Two sets are tested, when Sigma or u is variation.
 * 
 * \f[
 * \int_\Gamma n_i \Sigma_{ij} u_j \textrm{d}\Gamma
 * =
 * \int_\Omega \left( \Sigma_{ij} u_j \right)_i \textrm{d}\Omega
 * =
 * \int_\Omega \Sigma_{ij,i} u_j \textrm{d}\Omega
 * +
 * \int_\Omega \Sigma_{ij} u_{j,i} \textrm{d}\Omega
 * \f]
 *
 */

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

constexpr AssemblyType A = AssemblyType::PETSC; //< selected assembly type
constexpr IntegrationType I =
    IntegrationType::GAUSS;                     //< selected integration type

constexpr int SPACE_DIM = EXECUTABLE_DIMENSION;
constexpr CoordinateTypes COORD_TYPE = EXECUTABLE_COORD_TYPE;

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
    int order = 4;
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
    // CHKERR simple->setFieldOrder("SIGMA", order);
    CHKERR simple->setFieldOrder("GEOMETRY", 2);

    auto get_skin = [&]() {
      Range body_ents;
      CHKERR m_field.get_moab().get_entities_by_dimension(0, SPACE_DIM,
                                                         body_ents);
      Skinner skin(&m_field.get_moab());
      Range skin_ents;
      CHKERR skin.find_skin(0, body_ents, false, skin_ents);
      return skin_ents;
    };

    auto filter_true_skin = [&](auto skin) {
      Range boundary_ents;
      ParallelComm *pcomm =
          ParallelComm::get_pcomm(&m_field.get_moab(), MYPCOMM_INDEX);
      CHKERR pcomm->filter_pstatus(skin, PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                   PSTATUS_NOT, -1, &boundary_ents);
      return boundary_ents;
    };

    auto boundary_ents = filter_true_skin(get_skin());
    CHKERR simple->setFieldOrder("SIGMA", 0);
    CHKERR simple->setFieldOrder("SIGMA", order, &boundary_ents);

    // setup problem
    CHKERR simple->setUp();

    auto bc_mng = m_field.getInterface<BcManager>();
    CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "SYM",
                                             "U", 0, MAX_DOFS_ON_ENTITY, true);
    CHKERR bc_mng->removeBlockDOFsOnEntities(
        simple->getProblemName(), "SYM", "SIGMA", 0, MAX_DOFS_ON_ENTITY, true);

    auto project_ho_geometry = [&]() {
      Projection10NodeCoordsOnField ent_method(m_field, "GEOMETRY");
      return m_field.loop_dofs("GEOMETRY", ent_method);
    };
    CHKERR project_ho_geometry();

    // get operators tester
    auto opt = m_field.getInterface<OperatorsTester>(); // get interface to
                                                        // OperatorsTester
    auto pip_mng = m_field.getInterface<PipelineManager>(); // get interface to
                                                            // pipeline manager

    pip_mng->getOpDomainLhsPipeline().clear();
    pip_mng->getOpDomainRhsPipeline().clear();

    // integration rule
    auto rule = [&](int, int, int p) { return 2 * p + 2; };
    CHKERR pip_mng->setDomainRhsIntegrationRule(rule);
    CHKERR pip_mng->setDomainLhsIntegrationRule(rule);
    CHKERR pip_mng->setBoundaryRhsIntegrationRule(rule);

    auto x = opt->setRandomFields(simple->getDM(),
                                  {{"U", {-1, 1}}, {"SIGMA", {-1, 1}}});
    CHKERR DMoFEMMeshToLocalVector(simple->getDM(), x, INSERT_VALUES,
                                   SCATTER_REVERSE);

    // set integration rules

    auto jacobian = [&](double r) {
      if constexpr (COORD_TYPE == CYLINDRICAL)
        return 2 * M_PI * r;
      else
        return 1.;
    };

    auto beta_domain = [jacobian](double r, double, double) {
      return jacobian(r);
    };
    auto beta_bdy = [jacobian](double r, double, double) {
      return -jacobian(r);
    };

    auto post_proc = [&](auto dm, auto f_res, auto out_name) {
      MoFEMFunctionBegin;
      auto post_proc_fe =
          boost::make_shared<PostProcBrokenMeshInMoab<DomainEle>>(m_field);

      using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

      auto sigma_ptr = boost::make_shared<MatrixDouble>();
      post_proc_fe->getOpPtrVector().push_back(
          new OpCalculateHVecTensorField<SPACE_DIM, SPACE_DIM>("SIGMA",
                                                               sigma_ptr));
      auto u_ptr = boost::make_shared<MatrixDouble>();
      post_proc_fe->getOpPtrVector().push_back(
          new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));

      post_proc_fe->getOpPtrVector().push_back(

          new OpPPMap(

              post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),
              {}, {{"U", u_ptr}}, {{"SIGMA", sigma_ptr}}, {})

      );

      CHKERR DMoFEMMeshToLocalVector(simple->getDM(), f_res, INSERT_VALUES,
                                     SCATTER_REVERSE);
      CHKERR DMoFEMLoopFiniteElements(dm, simple->getDomainFEName(),
                                      post_proc_fe);
      post_proc_fe->writeFile(out_name);
      MoFEMFunctionReturn(0);
    };

    auto test_consistency_of_domain_and_bdy_integrals = [&]() {
      MoFEMFunctionBegin;
      using OpMixDivURhs =
          FormsIntegrators<DomainEleOp>::Assembly<A>::LinearForm<
              GAUSS>::OpMixDivTimesU<3, SPACE_DIM, SPACE_DIM, COORD_TYPE>;
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
        pip.push_back(
            new OpMixLambdaGradURhs("SIGMA", grad_u_ptr, beta_domain));

        auto sigma_ptr = boost::make_shared<MatrixDouble>();
        auto sigma_div_ptr = boost::make_shared<MatrixDouble>();
        pip.push_back(new OpCalculateHVecTensorDivergence<SPACE_DIM, SPACE_DIM,
                                                          COORD_TYPE>(
            "SIGMA", sigma_div_ptr));
        pip.push_back(new OpCalculateHVecTensorField<SPACE_DIM, SPACE_DIM>(
            "SIGMA", sigma_ptr));
        pip.push_back(
            new OpMixUTimesDivLambdaRhs("U", sigma_div_ptr, beta_domain));
        pip.push_back(new OpMixUTimesLambdaRhs("U", sigma_ptr, beta_domain));

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

      CHKERR DMoFEMLoopFiniteElements(simple->getDM(),
                                      simple->getDomainFEName(),
                                      pip_mng->getDomainRhsFE());
      CHKERR DMoFEMLoopFiniteElements(simple->getDM(),
                                      simple->getBoundaryFEName(),
                                      pip_mng->getBoundaryRhsFE());

      CHKERR VecAssemblyBegin(f);
      CHKERR VecAssemblyEnd(f);

      double f_nrm2;
      CHKERR VecNorm(f, NORM_2, &f_nrm2);

      MOFEM_LOG("ATOM", Sev::inform) << "f_norm2 = " << f_nrm2;
      if (std::fabs(f_nrm2) > 1e-10) {
        CHKERR post_proc(simple->getDM(), f,
                         "tensor_divergence_operator_res_vec.h5m");
        SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "Test failed");
      }

      MoFEMFunctionReturn(0);
    };

    // Testing Lhs domain operators

    auto test_lhs_ops = [&]() {
      MoFEMFunctionBegin;
      using OpMixDivULhs = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::BiLinearForm<I>::OpMixDivTimesVec<SPACE_DIM, COORD_TYPE>;
      using OpLambdaGraULhs = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::BiLinearForm<I>::OpMixTensorTimesGrad<SPACE_DIM>;

      auto op_lhs_domain = [&](auto &pip) {
        MoFEMFunctionBegin;
        CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1, HDIV},
                                                              "GEOMETRY");

        auto unity = []() { return 1; };
        pip.push_back(
            new OpMixDivULhs("SIGMA", "U", unity, beta_domain, true, false));
        pip.push_back(
            new OpLambdaGraULhs("SIGMA", "U", unity, beta_domain, true, false));
        MoFEMFunctionReturn(0);
      };

      CHKERR op_lhs_domain(pip_mng->getOpDomainLhsPipeline());

      auto diff_x = opt->setRandomFields(simple->getDM(),
                                         {{"U", {-1, 1}}, {"SIGMA", {-1, 1}}});
      constexpr double eps = 1e-5;
      auto diff_res = opt->checkCentralFiniteDifference(
          simple->getDM(), simple->getDomainFEName(), pip_mng->getDomainRhsFE(),
          pip_mng->getDomainLhsFE(), x, SmartPetscObj<Vec>(),
          SmartPetscObj<Vec>(), diff_x, 0, 1, eps);

      // Calculate norm of difference between directive calculated from finite
      // difference, and tangent matrix.
      double fnorm_res;
      CHKERR VecNorm(diff_res, NORM_2, &fnorm_res);
      MOFEM_LOG_C("ATOM", Sev::inform, "Test Lhs OPs %3.4e", fnorm_res);
      if (std::fabs(fnorm_res) > 1e-8)
        SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "Test failed"); 

      MoFEMFunctionReturn(0);
    };

    CHKERR test_consistency_of_domain_and_bdy_integrals();
    CHKERR test_lhs_ops();
}
CATCH_ERRORS;

CHKERR MoFEM::Core::Finalize();
}
