/**
 * \example test_broken_space.cpp
 *
 * Testing broken spaces
 *
 */

#include <MoFEM.hpp>

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
    CHKERR simple->addSkeletonField("HYBRID", L2, base, SPACE_DIM);
    CHKERR simple->addBoundaryField("HYBRID", L2, base, SPACE_DIM);
    CHKERR simple->setFieldOrder("BROKEN", approx_order);
    CHKERR simple->setFieldOrder("HYBRID", approx_order);

    CHKERR simple->setUp();

    auto assemble_domain_lhs = [&](auto &pip) {
      MoFEMFunctionBegin;
      using OpMass = FormsIntegrators<DomainEleOp>::Assembly<
          SCHUR>::BiLinearForm<GAUSS>::OpMass<3, 3>;
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {HDIV});
      pip.push_back(new OpMass("BROKEN", "BROKEN",
                               [](double, double, double) { return 1.; }));
      MoFEMFunctionReturn(0);
    };

    auto assemble_domain_rhs = [&](auto &pip) {
      MoFEMFunctionBegin;
      using OpSource = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::LinearForm<GAUSS>::OpSource<3, 3>;
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {HDIV});
      pip.push_back(new OpSource("BROKEN", [](double, double, double) {
        return FTensor::Tensor1<double, 3>{1., 0., 0.};
      }));

      MoFEMFunctionReturn(0);
    };

    auto get_broken_ptr = [&]() {
      auto broken_data_ptr = boost::make_shared<BrokenBaseSideData>();
      auto flux_mat_ptr = boost::make_shared<MatrixDouble>();
      auto op_loop_side = new OpLoopSide<EleOnSide>(
          m_field, simple->getDomainFEName(), SPACE_DIM, Sev::noisy);
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
          op_loop_side->getOpPtrVector(), {HDIV});
      op_loop_side->getOpPtrVector().push_back(
          new OpGetBrokenBaseSideData("BROKEN", broken_data_ptr));
      op_loop_side->getOpPtrVector().push_back(
          new OpCalculateHVecVectorField<3>("BROKEN", flux_mat_ptr));
      return std::make_tuple(op_loop_side, broken_data_ptr, flux_mat_ptr);
    };

    auto assemble_skeleton_lhs = [&](auto &pip, auto &&broken_data_tuple) {
      MoFEMFunctionBegin;
      using OpC = FormsIntegrators<BdyEleOp>::Assembly<PETSC>::BiLinearForm<
          GAUSS>::OpBrokenSpaceConstrain<SPACE_DIM>;
      auto [op_loop_side, broken_data_ptr, flux_mat_ptr] = broken_data_tuple;
      pip.push_back(op_loop_side);
      CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
          op_loop_side->getOpPtrVector(), {});
      pip.push_back(new OpC("HYBRID", broken_data_ptr, 1., true, false));
      MoFEMFunctionReturn(0);
    };

    auto assemble_skeleton_rhs = [&](auto &pip, auto &&broken_data_tuple) {
      MoFEMFunctionBegin;
      using OpC_dHybrid = FormsIntegrators<BdyEleOp>::Assembly<
          PETSC>::LinearForm<GAUSS>::OpBrokenSpaceConstrainDHybrid<SPACE_DIM>;
      using OpC_dBroken = FormsIntegrators<BdyEleOp>::Assembly<
          PETSC>::LinearForm<GAUSS>::OpBrokenSpaceConstrainDFlux<SPACE_DIM>;
      auto [op_loop_side, broken_data_ptr, flux_mat_ptr] = broken_data_tuple;
      pip.push_back(op_loop_side);
      CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
          op_loop_side->getOpPtrVector(), {});
      auto hybrid_ptr = boost::make_shared<MatrixDouble>();
      pip.push_back(
          new OpCalculateVectorFieldValues<SPACE_DIM>("HYBRID", hybrid_ptr));
      // pip.push_back(new OpC_dHybrid("HYBRID", flux_mat_ptr, 1.));
      // pip.push_back(new OpC_dBroken(broken_data_ptr, hybrid_ptr, 1.));
      MoFEMFunctionReturn(0);
    };

    auto *pip_mng = m_field.getInterface<PipelineManager>();

    CHKERR assemble_domain_lhs(pip_mng->getOpDomainLhsPipeline());
    CHKERR assemble_domain_rhs(pip_mng->getOpDomainLhsPipeline());
    CHKERR assemble_skeleton_lhs(pip_mng->getOpSkeletonLhsPipeline(),
                                 get_broken_ptr());
    CHKERR assemble_skeleton_rhs(pip_mng->getOpSkeletonRhsPipeline(),
                                 get_broken_ptr());

    auto integration_rule = [](int, int, int p_data) { return 2 * p_data + 1; };
    CHKERR pip_mng->setDomainRhsIntegrationRule(integration_rule);
    CHKERR pip_mng->setDomainLhsIntegrationRule(integration_rule);
    CHKERR pip_mng->setSkeletonLhsIntegrationRule(integration_rule);
    CHKERR pip_mng->setSkeletonRhsIntegrationRule(integration_rule);

    auto ksp = pip_mng->createKSP();

    CHKERR KSPSetFromOptions(ksp);
    CHKERR KSPSetUp(ksp);

    auto x = createDMVector(simple->getDM());
    auto f = vectorDuplicate(x);

    CHKERR KSPSolve(ksp, x, f);
    CHKERR VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    CHKERR VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
    CHKERR DMoFEMMeshToLocalVector(simple->getDM(), x, INSERT_VALUES,
                                   SCATTER_REVERSE);

    auto get_post_proc_fe = [&]() {
      using PostProcEle = PostProcBrokenMeshInMoab<DomainEle>;
      using OpPPMap = OpPostProcMapInMoab<3, SPACE_DIM>;
      auto post_proc_fe = boost::make_shared<PostProcEle>(m_field);

      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
          post_proc_fe->getOpPtrVector(), {HDIV});
      auto flux_mat_ptr = boost::make_shared<MatrixDouble>();
      post_proc_fe->getOpPtrVector().push_back(
          new OpCalculateHVecVectorField<3>("BROKEN", flux_mat_ptr));

      post_proc_fe->getOpPtrVector().push_back(

          new OpPPMap(

              post_proc_fe->getPostProcMesh(),

              post_proc_fe->getMapGaussPts(),

              {},

              {{"BROKEN", flux_mat_ptr}},

              {}, {})

      );

      return post_proc_fe;
    };

    auto post_proc_fe = get_post_proc_fe();
    CHKERR DMoFEMLoopFiniteElements(simple->getDM(), simple->getDomainFEName(),
                                    post_proc_fe);
    CHKERR post_proc_fe->writeFile("out_result.h5m");
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
