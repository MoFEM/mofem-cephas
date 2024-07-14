/**
 * \example test_broken_space.cpp
 *
 * Testing broken spaces
 *
 */

#include <MoFEM.hpp>
#include <FormsBrokenSpaceConstraintImpl.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

constexpr int SPACE_DIM =
    EXECUTABLE_DIMENSION; //< Space dimension of problem, mesh

using DomainEle = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::DomainEle;
using BoundaryEle =
    PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::BoundaryEle;
using EleOnSide = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::FaceSideEle;

using EntData = EntitiesFieldData::EntData;
using DomainEleOp = DomainEle::UserDataOperator;
using BdyEleOp = BoundaryEle::UserDataOperator;
using SideEleOp = EleOnSide::UserDataOperator;

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    // Add logging channel for example
    auto core_log = logging::core::get();
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "AT"));
    LogManager::setLog("AT");
    MOFEM_LOG_TAG("AT", "atom_test");

    // Create MoFEM instance
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    auto *simple = m_field.getInterface<Simple>();
    CHKERR simple->getOptions();

    simple->getAddSkeletonFE() = true;
    simple->getAddBoundaryFE() = true;

    CHKERR simple->loadFile("", "");

    // Declare elements
    enum bases {
      AINSWORTH,
      AINSWORTH_LOBATTO,
      DEMKOWICZ,
      BERNSTEIN,
      LASBASETOP
    };
    const char *list_bases[] = {"ainsworth", "ainsworth_labatto", "demkowicz",
                                "bernstein"};
    PetscBool flg;
    PetscInt choice_base_value = AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list_bases,
                                LASBASETOP, &choice_base_value, &flg);

    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "base not set");
    FieldApproximationBase base = AINSWORTH_LEGENDRE_BASE;
    if (choice_base_value == AINSWORTH)
      base = AINSWORTH_LEGENDRE_BASE;
    if (choice_base_value == AINSWORTH_LOBATTO)
      base = AINSWORTH_LOBATTO_BASE;
    else if (choice_base_value == DEMKOWICZ)
      base = DEMKOWICZ_JACOBI_BASE;
    else if (choice_base_value == BERNSTEIN)
      base = AINSWORTH_BERNSTEIN_BEZIER_BASE;

    enum spaces { hdiv, hcurl, last_space };
    const char *list_spaces[] = {"hdiv", "hcurl"};
    PetscInt choice_space_value = hdiv;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-space", list_spaces,
                                last_space, &choice_space_value, &flg);
    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "space not set");
    FieldSpace space = HDIV;
    if (choice_space_value == hdiv)
      space = HDIV;
    else if (choice_space_value == hcurl)
      space = HCURL;

    int approx_order = 1;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &approx_order,
                              PETSC_NULL);

    CHKERR simple->addDomainBrokenField("BROKEN", space, base, 1);
    CHKERR simple->addDomainField("U", L2, base, 1);
    CHKERR simple->addSkeletonField("HYBRID", L2, base, 1);

    CHKERR simple->setFieldOrder("BROKEN", approx_order);
    CHKERR simple->setFieldOrder("U", approx_order - 1);
    CHKERR simple->setFieldOrder("HYBRID", approx_order - 1);

    CHKERR simple->setUp();

    auto assemble_domain_lhs = [&](auto &pip) {
      MoFEMFunctionBegin;
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {HDIV});

      using OpHdivHdiv = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::BiLinearForm<GAUSS>::OpMass<3, SPACE_DIM>;
      using OpHdivU = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::BiLinearForm<GAUSS>::OpMixDivTimesScalar<SPACE_DIM>;

      using OpMass = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::BiLinearForm<GAUSS>::OpMass<1, 1>;

      auto beta = [](const double, const double, const double) constexpr {
        return 1;
      };

      pip.push_back(new OpHdivHdiv("BROKEN", "BROKEN", beta));
      auto unity = []() constexpr { return 1; };
      pip.push_back(new OpHdivU("BROKEN", "U", unity, true));

      MoFEMFunctionReturn(0);
    };

    auto assemble_domain_rhs = [&](auto &pip, double alpha) {
      MoFEMFunctionBegin;
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {HDIV});
      using OpDomainSource = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::LinearForm<GAUSS>::OpSource<1, 1>;
      auto source = [&](const double x, const double y,
                        const double z) constexpr { return -1; };
      pip.push_back(new OpDomainSource("U", source));
      MoFEMFunctionReturn(0);
    };

    auto get_broken_ptr = [&]() {
      auto broken_data_ptr =
          boost::make_shared<std::vector<BrokenBaseSideData>>();
      auto flux_mat_ptr = boost::make_shared<MatrixDouble>();
      auto op_loop_side = new OpLoopSide<EleOnSide>(
          m_field, simple->getDomainFEName(), SPACE_DIM, Sev::noisy);
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
          op_loop_side->getOpPtrVector(), {HDIV});
      op_loop_side->getOpPtrVector().push_back(
          new OpGetBrokenBaseSideData<SideEleOp>("BROKEN", broken_data_ptr));
      op_loop_side->getOpPtrVector().push_back(
          new OpCalculateHVecTensorField<1, 3>("BROKEN", flux_mat_ptr));
      op_loop_side->getOpPtrVector().push_back(
          new OpSetFlux<SideEleOp>(broken_data_ptr, flux_mat_ptr));
      return std::make_tuple(op_loop_side, broken_data_ptr);
    };

    auto assemble_skeleton_lhs = [&](auto &pip, auto &&broken_data_tuple) {
      MoFEMFunctionBegin;
      using OpC = FormsIntegrators<BdyEleOp>::Assembly<PETSC>::BiLinearForm<
          GAUSS>::OpBrokenSpaceConstrain<1>;
      auto [op_loop_side, broken_data_ptr] = broken_data_tuple;
      pip.push_back(op_loop_side);
      CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
          op_loop_side->getOpPtrVector(), {});
      pip.push_back(new OpC("HYBRID", broken_data_ptr, 1., true, false));
      MoFEMFunctionReturn(0);
    };

    auto *pip_mng = m_field.getInterface<PipelineManager>();

    CHKERR assemble_domain_lhs(pip_mng->getOpDomainLhsPipeline());
    CHKERR assemble_skeleton_lhs(pip_mng->getOpSkeletonLhsPipeline(),
                                 get_broken_ptr());
    CHKERR assemble_domain_rhs(pip_mng->getOpDomainRhsPipeline(), 1);

    auto integration_rule = [](int, int, int p) { return 2 * p + 1; };
    CHKERR pip_mng->setDomainLhsIntegrationRule(integration_rule);
    CHKERR pip_mng->setDomainRhsIntegrationRule(integration_rule);
    CHKERR pip_mng->setSkeletonLhsIntegrationRule(integration_rule);
    CHKERR pip_mng->setSkeletonRhsIntegrationRule(integration_rule);

    auto ksp = pip_mng->createKSP();

    CHKERR KSPSetFromOptions(ksp);
    CHKERR KSPSetUp(ksp);

    auto x = createDMVector(simple->getDM());
    auto f = vectorDuplicate(x);

    CHKERR KSPSolve(ksp, f, x);
    CHKERR VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    CHKERR VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
    CHKERR DMoFEMMeshToLocalVector(simple->getDM(), x, INSERT_VALUES,
                                   SCATTER_REVERSE);

    auto check_residual = [&](auto x, auto f) {
      MoFEMFunctionBegin;
      auto *simple = m_field.getInterface<Simple>();
      auto *pip_mng = m_field.getInterface<PipelineManager>();

      auto &skeleton_rhs = pip_mng->getOpSkeletonRhsPipeline();
      auto &domain_rhs = pip_mng->getOpDomainRhsPipeline();
      skeleton_rhs.clear();
      domain_rhs.clear();

      CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(skeleton_rhs,
                                                                {});
      using OpC_dHybrid = FormsIntegrators<BdyEleOp>::Assembly<
          PETSC>::LinearForm<GAUSS>::OpBrokenSpaceConstrainDHybrid<1>;
      using OpC_dBroken = FormsIntegrators<BdyEleOp>::Assembly<
          PETSC>::LinearForm<GAUSS>::OpBrokenSpaceConstrainDFlux<1>;
      auto broken_data_tuple = get_broken_ptr();
      auto [op_loop_side, broken_data_ptr] = broken_data_tuple;
      skeleton_rhs.push_back(op_loop_side);
      CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
          op_loop_side->getOpPtrVector(), {});
      skeleton_rhs.push_back(new OpC_dHybrid("HYBRID", broken_data_ptr, 1.));
      auto hybrid_ptr = boost::make_shared<MatrixDouble>();
      skeleton_rhs.push_back(
          new OpCalculateVectorFieldValues<1>("HYBRID", hybrid_ptr));

      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(domain_rhs, {HDIV});
      auto flux_ptr = boost::make_shared<MatrixDouble>();
      domain_rhs.push_back(
          new OpCalculateHVecVectorField<3>("BROKEN", flux_ptr));
      auto div_flux_ptr = boost::make_shared<VectorDouble>();
      domain_rhs.push_back(new OpCalculateHdivVectorDivergence<3, SPACE_DIM>(
          "BROKEN", div_flux_ptr));
      using OpUDivFlux = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::LinearForm<GAUSS>::OpBaseTimesScalarField<1>;
      auto beta = [](double, double, double) constexpr { return 1; };
      domain_rhs.push_back(new OpUDivFlux("U", div_flux_ptr, beta));
      auto source = [&](const double x, const double y,
                        const double z) constexpr { return 1; };
      using OpDomainSource = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::LinearForm<GAUSS>::OpSource<1, 1>;
      domain_rhs.push_back(new OpDomainSource("U", source));

      using OpHDivH = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::LinearForm<GAUSS>::OpMixDivTimesU<3, 1, SPACE_DIM>;
      domain_rhs.push_back(new OpHdivHdiv("BROKEN", "BROKEN", beta));
      

      auto u_ptr = boost::make_shared<VectorDouble>();
      domain_rhs.push_back(new OpCalculateScalarFieldValues("U", u_ptr));

      CHKERR VecZeroEntries(f);
      CHKERR VecGhostUpdateBegin(f, INSERT_VALUES, SCATTER_FORWARD);
      CHKERR VecGhostUpdateEnd(f, INSERT_VALUES, SCATTER_FORWARD);

      pip_mng->getDomainRhsFE()->f = f;
      pip_mng->getSkeletonRhsFE()->f = f;
      pip_mng->getDomainRhsFE()->x = x;
      pip_mng->getSkeletonRhsFE()->x = x;

      CHKERR DMoFEMLoopFiniteElements(simple->getDM(),
                                      simple->getDomainFEName(),
                                      pip_mng->getDomainRhsFE());
      CHKERR DMoFEMLoopFiniteElements(simple->getDM(),
                                      simple->getSkeletonFEName(),
                                      pip_mng->getSkeletonRhsFE());

      CHKERR VecGhostUpdateBegin(f, ADD_VALUES, SCATTER_REVERSE);
      CHKERR VecGhostUpdateEnd(f, ADD_VALUES, SCATTER_REVERSE);
      CHKERR VecAssemblyBegin(f);
      CHKERR VecAssemblyEnd(f);

      double fnrm;
      CHKERR VecNorm(f, NORM_2, &fnrm);
      MOFEM_LOG_C("AT", Sev::inform, "Residual %3.4e", fnrm);

      constexpr double eps = 1e-8;
      if (fnrm > eps)
        SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                "Residual norm larger than accepted");

      MoFEMFunctionReturn(0);
    };

    auto get_post_proc_fe = [&]() {
      using PostProcEle = PostProcBrokenMeshInMoab<DomainEle>;
      using OpPPMap = OpPostProcMapInMoab<3, SPACE_DIM>;
      auto post_proc_fe = boost::make_shared<PostProcEle>(m_field);

      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
          post_proc_fe->getOpPtrVector(), {HDIV});
      auto u_vec_ptr = boost::make_shared<VectorDouble>();
      auto flux_mat_ptr = boost::make_shared<MatrixDouble>();
      post_proc_fe->getOpPtrVector().push_back(
          new OpCalculateScalarFieldValues("U", u_vec_ptr));
      post_proc_fe->getOpPtrVector().push_back(
          new OpCalculateHVecVectorField<3>("BROKEN", flux_mat_ptr));

      post_proc_fe->getOpPtrVector().push_back(

          new OpPPMap(

              post_proc_fe->getPostProcMesh(),

              post_proc_fe->getMapGaussPts(),

              {{"U", u_vec_ptr}},

              {{"BROKEN", flux_mat_ptr}},

              {}, {})

      );

      return post_proc_fe;
    };

    auto post_proc_fe = get_post_proc_fe();
    CHKERR DMoFEMLoopFiniteElements(simple->getDM(), simple->getDomainFEName(),
                                    post_proc_fe);
    CHKERR post_proc_fe->writeFile("out_result.h5m");

    CHKERR check_residual(x, f);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
