/**
 * \example test_broken_space.cpp
 *
 * Testing broken spaces. Implementations works for 2d and 3d meshes, is aiming
 * to test H-div broken base functions, and L2 base on skeleton.
 *
 * Also, it test block matrix with Schur complement.
 *
 */

#include <MoFEM.hpp>
#include <FormsBrokenSpaceConstraintImpl.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

constexpr bool debug = false;

constexpr AssemblyType AT =
    (SCHUR_ASSEMBLE) ? AssemblyType::BLOCK_SCHUR
                     : AssemblyType::PETSC; //< selected assembly type

constexpr IntegrationType IT =
    IntegrationType::GAUSS; //< selected integration type

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

struct SetUpSchur {
  static boost::shared_ptr<SetUpSchur>
  createSetUpSchur(MoFEM::Interface &m_field);
  virtual MoFEMErrorCode setUp(SmartPetscObj<KSP>) = 0;

protected:
  SetUpSchur() = default;
  virtual ~SetUpSchur() = default;
};

int approx_order = 1;

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    //! [Register MoFEM discrete manager in PETSc]
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);
    DMType dm_name_mg = "DMMOFEM_MG";
    CHKERR DMRegister_MGViaApproxOrders(dm_name_mg);
    //! [Register MoFEM discrete manager in PETSc

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    // Add logging channel for example
    auto core_log = logging::core::get();
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "AT"));
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "TIMER"));
    LogManager::setLog("AT");
    LogManager::setLog("TIMER");
    MOFEM_LOG_TAG("AT", "atom_test");
    MOFEM_LOG_TAG("TIMER", "timer");

    // Create MoFEM instance
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    auto *simple = m_field.getInterface<Simple>();
    CHKERR simple->getOptions();
    simple->getAddBoundaryFE() = true;
    CHKERR simple->loadFile();

    auto add_shared_entities_on_skeleton = [&]() {
      MoFEMFunctionBegin;
      auto boundary_meshset = simple->getBoundaryMeshSet();
      auto skeleton_meshset = simple->getSkeletonMeshSet();
      Range bdy_ents;
      CHKERR m_field.get_moab().get_entities_by_handle(boundary_meshset,
                                                       bdy_ents, true);
      Range skeleton_ents;
      CHKERR m_field.get_moab().get_entities_by_dimension(
          0, simple->getDim() - 1, skeleton_ents, true);
      skeleton_ents = subtract(skeleton_ents, bdy_ents);
      CHKERR m_field.get_moab().clear_meshset(&skeleton_meshset, 1);
      CHKERR m_field.get_moab().add_entities(skeleton_meshset, skeleton_ents);
      MoFEMFunctionReturn(0);
    };

    CHKERR add_shared_entities_on_skeleton();

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

    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &approx_order,
                              PETSC_NULL);

    CHKERR simple->addDomainBrokenField("BROKEN", space, base, 1);
    CHKERR simple->addDomainField("U", L2, base, 1);
    CHKERR simple->addSkeletonField("HYBRID", L2, base, 1);

    CHKERR simple->setFieldOrder("BROKEN", approx_order);
    CHKERR simple->setFieldOrder("U", approx_order - 1);
    CHKERR simple->setFieldOrder("HYBRID", approx_order - 1);

    CHKERR simple->setUp();

    auto integration_rule = [](int, int, int p) { return 2 * p; };

    auto assemble_domain_lhs = [&](auto &pip) {
      MoFEMFunctionBegin;
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {HDIV});

      using OpHdivHdiv = FormsIntegrators<DomainEleOp>::Assembly<
          AT>::BiLinearForm<IT>::OpMass<3, SPACE_DIM>;
      using OpHdivU = FormsIntegrators<DomainEleOp>::Assembly<AT>::BiLinearForm<
          IT>::OpMixDivTimesScalar<SPACE_DIM>;

      auto beta = [](const double, const double, const double) constexpr {
        return 1;
      };

      pip.push_back(new OpHdivHdiv("BROKEN", "BROKEN", beta));
      auto unity = []() constexpr { return 1; };
      pip.push_back(new OpHdivU("BROKEN", "U", unity, true));

      // First: Iterate over skeleton FEs adjacent to Domain FEs
      // Note:  BoundaryEle, i.e. uses skeleton interation rule
      auto op_loop_skeleton_side = new OpLoopSide<BoundaryEle>(
          m_field, simple->getSkeletonFEName(), SPACE_DIM - 1, Sev::noisy);
      op_loop_skeleton_side->getSideFEPtr()->getRuleHook = integration_rule;
      CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
          op_loop_skeleton_side->getOpPtrVector(), {});

      // Second: Iterate over domain FEs adjacent to skelton, particularly one
      // domain element.
      auto broken_data_ptr =
          boost::make_shared<std::vector<BrokenBaseSideData>>();
      // Note: EleOnSide, i.e. uses on domain projected skeleton rule
      auto op_loop_domain_side = new OpBrokenLoopSide<EleOnSide>(
          m_field, simple->getDomainFEName(), SPACE_DIM, Sev::noisy);
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
          op_loop_domain_side->getOpPtrVector(), {HDIV});
      op_loop_domain_side->getOpPtrVector().push_back(
          new OpGetBrokenBaseSideData<SideEleOp>("BROKEN", broken_data_ptr));

      op_loop_skeleton_side->getOpPtrVector().push_back(op_loop_domain_side);
      using OpC = FormsIntegrators<BdyEleOp>::Assembly<AT>::BiLinearForm<
          IT>::OpBrokenSpaceConstrain<1>;
      op_loop_skeleton_side->getOpPtrVector().push_back(
          new OpC("HYBRID", broken_data_ptr, 1., true, false));

      if (debug) {
        // print skeleton elements on partition
        constexpr int partition = 1;
        auto op_print = new BdyEleOp(NOSPACE, BdyEleOp::OPSPACE);
        op_print->doWorkRhsHook = [&](DataOperator *base_op_ptr, int side,
                                      EntityType type,
                                      EntitiesFieldData::EntData &data) {
          MoFEMFunctionBegin;
          if (auto op_ptr = dynamic_cast<BdyEleOp *>(base_op_ptr)) {
            auto fe_method = op_ptr->getFEMethod();
            auto num_fe = fe_method->numeredEntFiniteElementPtr;

            if (m_field.get_comm_rank() == partition) {
              if (num_fe->getPStatus() & PSTATUS_SHARED)
                MOFEM_LOG("SELF", Sev::inform) << "Num FE: " << *num_fe;
            }
          }
          MoFEMFunctionReturn(0);
        };
        op_loop_skeleton_side->getOpPtrVector().push_back(op_print);
      };

      pip.push_back(op_loop_skeleton_side);

      MoFEMFunctionReturn(0);
    };

    auto assemble_domain_rhs = [&](auto &pip) {
      MoFEMFunctionBegin;
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {HDIV});
      using OpDomainSource = FormsIntegrators<DomainEleOp>::Assembly<
          AT>::LinearForm<IT>::OpSource<1, 1>;
      auto source = [&](const double x, const double y,
                        const double z) constexpr {
        return -1;//sin(100 * (x / 10.) * M_PI_2);
      };
      pip.push_back(new OpDomainSource("U", source));
      MoFEMFunctionReturn(0);
    };

    auto *pip_mng = m_field.getInterface<PipelineManager>();

    CHKERR assemble_domain_lhs(pip_mng->getOpDomainLhsPipeline());
    CHKERR assemble_domain_rhs(pip_mng->getOpDomainRhsPipeline());

    CHKERR pip_mng->setDomainLhsIntegrationRule(integration_rule);
    CHKERR pip_mng->setDomainRhsIntegrationRule(integration_rule);
    CHKERR pip_mng->setSkeletonLhsIntegrationRule(integration_rule);
    CHKERR pip_mng->setSkeletonRhsIntegrationRule(integration_rule);

    TetPolynomialBase::switchCacheHDivBaseOn(
        {pip_mng->getDomainLhsFE().get(), pip_mng->getDomainRhsFE().get()});

    auto x = createDMVector(simple->getDM());
    auto f = vectorDuplicate(x);

    if (AT == PETSC) {
      auto ksp = pip_mng->createKSP();

      CHKERR KSPSetFromOptions(ksp);
      BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
      MOFEM_LOG("TIMER", Sev::inform) << "KSPSetUp";
      CHKERR KSPSetUp(ksp);
      MOFEM_LOG("TIMER", Sev::inform) << "KSPSetUp <= Done";

      MOFEM_LOG("TIMER", Sev::inform) << "KSPSolve";
      CHKERR KSPSolve(ksp, f, x);
      MOFEM_LOG("TIMER", Sev::inform) << "KSPSolve <= Done";

      CHKERR VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
      CHKERR VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
      CHKERR DMoFEMMeshToLocalVector(simple->getDM(), x, INSERT_VALUES,
                                     SCATTER_REVERSE);
    } else {
      auto ksp = pip_mng->createKSP();
      auto schur_ptr = SetUpSchur::createSetUpSchur(m_field);
      BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
      MOFEM_LOG("TIMER", Sev::inform) << "KSPSetUp";
      CHKERR schur_ptr->setUp(ksp);
      MOFEM_LOG("TIMER", Sev::inform) << "KSPSetUp <= Done";

      MOFEM_LOG("TIMER", Sev::inform) << "KSPSolve";
      CHKERR KSPSolve(ksp, f, x);
      MOFEM_LOG("TIMER", Sev::inform) << "KSPSolve <= Done";

      CHKERR VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
      CHKERR VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
      CHKERR DMoFEMMeshToLocalVector(simple->getDM(), x, INSERT_VALUES,
                                     SCATTER_REVERSE);
    }

    auto check_residual = [&](auto x, auto f) {
      MoFEMFunctionBegin;
      auto *simple = m_field.getInterface<Simple>();
      auto *pip_mng = m_field.getInterface<PipelineManager>();

      // auto &skeleton_rhs = pip_mng->getOpSkeletonRhsPipeline();
      auto &domain_rhs = pip_mng->getOpDomainRhsPipeline();
      // skeleton_rhs.clear();
      domain_rhs.clear();

      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(domain_rhs, {HDIV});

      auto div_flux_ptr = boost::make_shared<VectorDouble>();
      domain_rhs.push_back(new OpCalculateHdivVectorDivergence<3, SPACE_DIM>(
          "BROKEN", div_flux_ptr));
      using OpUDivFlux = FormsIntegrators<DomainEleOp>::Assembly<
          AT>::LinearForm<IT>::OpBaseTimesScalarField<1>;
      auto beta = [](double, double, double) constexpr { return 1; };
      domain_rhs.push_back(new OpUDivFlux("U", div_flux_ptr, beta));
      auto source = [&](const double x, const double y,
                        const double z) constexpr { return 1; };
      using OpDomainSource = FormsIntegrators<DomainEleOp>::Assembly<
          AT>::LinearForm<IT>::OpSource<1, 1>;
      domain_rhs.push_back(new OpDomainSource("U", source));

      using OpHDivH = FormsIntegrators<DomainEleOp>::Assembly<AT>::LinearForm<
          IT>::OpMixDivTimesU<3, 1, SPACE_DIM>;
      using OpHdivFlux = FormsIntegrators<DomainEleOp>::Assembly<
          AT>::LinearForm<IT>::OpBaseTimesVector<3, 3, 1>;
      auto flux_ptr = boost::make_shared<MatrixDouble>();
      domain_rhs.push_back(
          new OpCalculateHVecVectorField<3>("BROKEN", flux_ptr));
      boost::shared_ptr<VectorDouble> u_ptr =
          boost::make_shared<VectorDouble>();
      domain_rhs.push_back(new OpCalculateScalarFieldValues("U", u_ptr));
      auto minus = [](double, double, double) constexpr { return -1; };
      domain_rhs.push_back(new OpHDivH("BROKEN", u_ptr, beta));
      domain_rhs.push_back(new OpHdivFlux("BROKEN", flux_ptr, beta));

      // First: Iterate over skeleton FEs adjacent to Domain FEs
      // Note:  BoundaryEle, i.e. uses skeleton interation rule
      auto op_loop_skeleton_side = new OpLoopSide<BoundaryEle>(
          m_field, simple->getSkeletonFEName(), SPACE_DIM - 1, Sev::noisy);
      op_loop_skeleton_side->getSideFEPtr()->getRuleHook = integration_rule;
      CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
          op_loop_skeleton_side->getOpPtrVector(), {});

      // Second: Iterate over domain FEs adjacent to skelton, particularly one
      // domain element.
      auto broken_data_ptr =
          boost::make_shared<std::vector<BrokenBaseSideData>>();
      // Note: EleOnSide, i.e. uses on domain projected skeleton rule
      auto op_loop_domain_side = new OpBrokenLoopSide<EleOnSide>(
          m_field, simple->getDomainFEName(), SPACE_DIM, Sev::noisy);
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
          op_loop_domain_side->getOpPtrVector(), {HDIV});
      op_loop_domain_side->getOpPtrVector().push_back(
          new OpGetBrokenBaseSideData<SideEleOp>("BROKEN", broken_data_ptr));
      auto flux_mat_ptr = boost::make_shared<MatrixDouble>();
      op_loop_domain_side->getOpPtrVector().push_back(
          new OpCalculateHVecTensorField<1, 3>("BROKEN", flux_mat_ptr));
      op_loop_domain_side->getOpPtrVector().push_back(
          new OpSetFlux<SideEleOp>(broken_data_ptr, flux_mat_ptr));

      // Assemble on skeleton
      op_loop_skeleton_side->getOpPtrVector().push_back(op_loop_domain_side);
      using OpC_dHybrid = FormsIntegrators<BdyEleOp>::Assembly<AT>::LinearForm<
          IT>::OpBrokenSpaceConstrainDHybrid<1>;
      using OpC_dBroken = FormsIntegrators<BdyEleOp>::Assembly<AT>::LinearForm<
          IT>::OpBrokenSpaceConstrainDFlux<1>;
      op_loop_skeleton_side->getOpPtrVector().push_back(
          new OpC_dHybrid("HYBRID", broken_data_ptr, 1.));
      auto hybrid_ptr = boost::make_shared<MatrixDouble>();
      op_loop_skeleton_side->getOpPtrVector().push_back(
          new OpCalculateVectorFieldValues<1>("HYBRID", hybrid_ptr));
      op_loop_skeleton_side->getOpPtrVector().push_back(
          new OpC_dBroken(broken_data_ptr, hybrid_ptr, 1.));

      // Add skeleton to domain pipeline
      domain_rhs.push_back(op_loop_skeleton_side);

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

    auto calculate_error = [&]() {
      MoFEMFunctionBegin;

      // auto &skeleton_rhs = pip_mng->getOpSkeletonRhsPipeline();
      auto &domain_rhs = pip_mng->getOpDomainRhsPipeline();
      // skeleton_rhs.clear();
      domain_rhs.clear();

      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(domain_rhs, {HDIV});

      auto u_grad_ptr = boost::make_shared<MatrixDouble>();
      auto flux_val_ptr = boost::make_shared<MatrixDouble>();
      auto div_val_ptr = boost::make_shared<VectorDouble>();
      auto source_ptr = boost::make_shared<VectorDouble>();

      domain_rhs.push_back(
          new OpCalculateScalarFieldGradient<SPACE_DIM>("U", u_grad_ptr));
      domain_rhs.push_back(
          new OpCalculateHVecVectorField<3, SPACE_DIM>("BROKEN", flux_val_ptr));
      domain_rhs.push_back(new OpCalculateHdivVectorDivergence<3, SPACE_DIM>(
          "BROKEN", div_val_ptr));
      auto source = [&](const double x, const double y,
                        const double z) constexpr { return -1; };
      domain_rhs.push_back(new OpGetTensor0fromFunc(source_ptr, source));

      enum { DIV, GRAD, LAST};
      auto mpi_vec = createVectorMPI(
          m_field.get_comm(), (!m_field.get_comm_rank()) ? LAST : 0, LAST);
      domain_rhs.push_back(
          new OpCalcNormL2Tensor0(div_val_ptr, mpi_vec, DIV, source_ptr));
      domain_rhs.push_back(new OpCalcNormL2Tensor1<SPACE_DIM>(
          u_grad_ptr, mpi_vec, GRAD, flux_val_ptr));

      CHKERR DMoFEMLoopFiniteElements(simple->getDM(),
                                      simple->getDomainFEName(),
                                      pip_mng->getDomainRhsFE());
      CHKERR VecAssemblyBegin(mpi_vec);
      CHKERR VecAssemblyEnd(mpi_vec);

      if (!m_field.get_comm_rank()) {
        const double *error_ind;
        CHKERR VecGetArrayRead(mpi_vec, &error_ind);
        MOFEM_LOG("AT", Sev::inform)
            << "Approximation error ||div flux - source||: "
            << std::sqrt(error_ind[DIV]);
        MOFEM_LOG("AT", Sev::inform)
            << "Approximation error ||grad-flux||: "
            << std::sqrt(error_ind[GRAD]);
        CHKERR VecRestoreArrayRead(mpi_vec, &error_ind);
      }

      MoFEMFunctionReturn(0);
    };

    auto get_post_proc_fe = [&]() {
      using PostProcEle = PostProcBrokenMeshInMoab<BoundaryEle>;
      using OpPPMap = OpPostProcMapInMoab<3, SPACE_DIM>;
      auto post_proc_fe = boost::make_shared<PostProcEle>(m_field);

      auto op_loop_side = new OpLoopSide<EleOnSide>(
          m_field, simple->getDomainFEName(), SPACE_DIM, Sev::noisy,
          boost::make_shared<
              ForcesAndSourcesCore::UserDataOperator::AdjCache>());
      post_proc_fe->getOpPtrVector().push_back(op_loop_side);

      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
          op_loop_side->getOpPtrVector(), {HDIV});
      auto u_vec_ptr = boost::make_shared<VectorDouble>();
      auto flux_mat_ptr = boost::make_shared<MatrixDouble>();
      op_loop_side->getOpPtrVector().push_back(
          new OpCalculateScalarFieldValues("U", u_vec_ptr));
      op_loop_side->getOpPtrVector().push_back(
          new OpCalculateHVecVectorField<3>("BROKEN", flux_mat_ptr));

      op_loop_side->getOpPtrVector().push_back(

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
    CHKERR DMoFEMLoopFiniteElements(simple->getDM(),
                                    simple->getBoundaryFEName(), post_proc_fe);
    CHKERR post_proc_fe->writeFile("out_result.h5m");

    CHKERR calculate_error();
    CHKERR check_residual(x, f);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}

struct SetUpSchurImpl : public SetUpSchur {

  SetUpSchurImpl(MoFEM::Interface &m_field) : SetUpSchur(), mField(m_field) {}

  virtual ~SetUpSchurImpl() = default;

  MoFEMErrorCode setUp(SmartPetscObj<KSP>);

private:
  MoFEM::Interface &mField;
  SmartPetscObj<Mat> S;
};

MoFEMErrorCode SetUpSchurImpl::setUp(SmartPetscObj<KSP> ksp) {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  auto pip_mng = mField.getInterface<PipelineManager>();

  CHKERR KSPSetFromOptions(ksp);
  PC pc;
  CHKERR KSPGetPC(ksp, &pc);

  PetscBool is_pcfs = PETSC_FALSE;
  PetscObjectTypeCompare((PetscObject)pc, PCFIELDSPLIT, &is_pcfs);
  if (is_pcfs) {

    MOFEM_LOG("AT", Sev::inform) << "Setup Schur pc";

    auto create_sub_dm = [&]() {
      auto simple = mField.getInterface<Simple>();

      auto create_dm = [&](

                           std::string problem_name,
                           std::vector<std::string> fe_names,
                           std::vector<std::string> fields,

                           auto dm_type

                       ) {
        auto dm = createDM(mField.get_comm(), dm_type);
        auto create_dm_imp = [&]() {
          MoFEMFunctionBegin;
          CHKERR DMMoFEMCreateSubDM(dm, simple->getDM(), problem_name.c_str());
          CHKERR DMMoFEMSetSquareProblem(dm, PETSC_TRUE);
          for (auto fe : fe_names) {
            CHKERR DMMoFEMAddElement(dm, fe);
          }
          CHKERR DMMoFEMAddElement(dm, simple->getSkeletonFEName());
          for (auto field : fields) {
            CHKERR DMMoFEMAddSubFieldRow(dm, field);
            CHKERR DMMoFEMAddSubFieldCol(dm, field);
          }
          CHKERR DMSetUp(dm);
          MoFEMFunctionReturn(0);
        };
        CHK_THROW_MESSAGE(
            create_dm_imp(),
            "Error in creating schurDM. It is possible that schurDM is "
            "already created");
        return dm;
      };

      auto schur_dm = create_dm(

          "SCHUR",

          {simple->getDomainFEName(), simple->getSkeletonFEName()},

          {"HYBRID"},

          "DMMOFEM_MG");

      auto block_dm = create_dm(

          "BLOCK",

          {simple->getDomainFEName(), simple->getSkeletonFEName()},

          {"BROKEN", "U"},

          "DMMOFEM");

      return std::make_tuple(schur_dm, block_dm);
    };

    auto get_nested_mat_data = [&](auto schur_dm, auto block_dm) {
      auto block_mat_data = createBlockMatStructure(
          simple->getDM(),

          {

              {

                  simple->getDomainFEName(),

                  {

                      {"BROKEN", "BROKEN"},
                      {"U", "U"},
                      {"BROKEN", "U"},
                      {"U", "BROKEN"}

                  }

              },

              {

                  simple->getSkeletonFEName(),

                  {

                      {"BROKEN", "HYBRID"}, {"HYBRID", "BROKEN"}

                  }

              }

          }

      );

      return getNestSchurData(

          {schur_dm, block_dm}, block_mat_data,

          {"BROKEN", "U"}, {nullptr, nullptr}, true

      );
    };

    auto set_ops = [&](auto schur_dm) {
      MoFEMFunctionBegin;
      auto dm_is = getDMSubData(schur_dm)->getSmartRowIs();
      auto ao_up = createAOMappingIS(dm_is, PETSC_NULL);

      boost::shared_ptr<BlockStructure> block_data;
      CHKERR DMMoFEMGetBlocMatData(simple->getDM(), block_data);

      pip_mng->getOpDomainLhsPipeline().push_front(
          createOpSchurAssembleBegin());
      pip_mng->getOpDomainLhsPipeline().push_back(

          createOpSchurAssembleEnd(
              {"BROKEN", "U"}, {nullptr, nullptr}, {SmartPetscObj<AO>(), ao_up},
              {SmartPetscObj<Mat>(), S}, {true, true}, true, block_data)

      );

      auto pre_proc_schur_lhs_ptr = boost::make_shared<FEMethod>();
      auto post_proc_schur_lhs_ptr = boost::make_shared<FEMethod>();

      pre_proc_schur_lhs_ptr->preProcessHook = [this]() {
        MoFEMFunctionBegin;
        CHKERR MatZeroEntries(S);
        MOFEM_LOG("AT", Sev::verbose) << "Lhs Assemble Begin";
        MoFEMFunctionReturn(0);
      };

      post_proc_schur_lhs_ptr->postProcessHook = [this, ao_up,
                                                  post_proc_schur_lhs_ptr]() {
        MoFEMFunctionBegin;
        MOFEM_LOG("AT", Sev::verbose) << "Lhs Assemble End";
        CHKERR MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
        CHKERR MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);
        MOFEM_LOG("AT", Sev::verbose) << "Lhs Assemble Finish";
        MoFEMFunctionReturn(0);
      };

      auto ksp_ctx_ptr = getDMKspCtx(simple->getDM());
      ksp_ctx_ptr->getPreProcSetOperators().push_front(pre_proc_schur_lhs_ptr);
      ksp_ctx_ptr->getPostProcSetOperators().push_back(post_proc_schur_lhs_ptr);

      MoFEMFunctionReturn(0);
    };

    auto set_pc = [&](auto pc, auto block_dm) {
      MoFEMFunctionBegin;
      auto block_is = getDMSubData(block_dm)->getSmartRowIs();
      CHKERR PCFieldSplitSetIS(pc, NULL, block_is);
      CHKERR PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, S);
      MoFEMFunctionReturn(0);
    };

    auto set_diagonal_pc = [&](auto pc, auto schur_dm) {
      MoFEMFunctionBegin;

      auto A = createDMBlockMat(simple->getDM());
      auto P = createDMNestSchurMat(simple->getDM());
      CHKERR PCSetOperators(pc, A, P);

      KSP *subksp;
      CHKERR PCFieldSplitSchurGetSubKSP(pc, PETSC_NULL, &subksp);
      auto get_pc = [](auto ksp) {
        PC pc_raw;
        CHKERR KSPGetPC(ksp, &pc_raw);
        return pc_raw;
      };
      CHKERR setSchurA00MatSolvePC(SmartPetscObj<PC>(get_pc(subksp[0]), true));

      auto set_pc_p_mg = [&](auto dm, auto pc) {
        MoFEMFunctionBegin;

        CHKERR PCSetDM(pc, dm);
        PetscBool same = PETSC_FALSE;
        PetscObjectTypeCompare((PetscObject)pc, PCMG, &same);
        if (same) {
          // By default do not use shell mg mat. Implementation of SOR is slow.
          CHKERR PCMGSetUpViaApproxOrders(
              pc, createPCMGSetUpViaApproxOrdersCtx(dm, S, false));
          CHKERR PCSetFromOptions(pc);
        }
        MoFEMFunctionReturn(0);
      };

      CHKERR set_pc_p_mg(schur_dm, get_pc(subksp[1]));

      CHKERR PetscFree(subksp);
      MoFEMFunctionReturn(0);
    };

    auto [schur_dm, block_dm] = create_sub_dm();
    auto nested_mat_data = get_nested_mat_data(schur_dm, block_dm);
    CHKERR DMMoFEMSetNestSchurData(simple->getDM(), nested_mat_data);
    S = createDMHybridisedL2Matrix(schur_dm);
    CHKERR MatSetDM(S, PETSC_NULL);
    int bs = (SPACE_DIM == 2) ? NBEDGE_L2(approx_order - 1)
                              : NBFACETRI_L2(approx_order - 1);
    CHKERR MatSetBlockSize(S, bs);

    CHKERR set_ops(schur_dm);
    CHKERR set_pc(pc, block_dm);
    DM solver_dm;
    CHKERR KSPGetDM(ksp, &solver_dm);
    CHKERR DMSetMatType(solver_dm, MATSHELL);

    CHKERR KSPSetUp(ksp);
    CHKERR set_diagonal_pc(pc, schur_dm);

  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
            "PC is not set to PCFIELDSPLIT");
  }
  MoFEMFunctionReturn(0);
}

boost::shared_ptr<SetUpSchur>
SetUpSchur::createSetUpSchur(MoFEM::Interface &m_field) {
  return boost::shared_ptr<SetUpSchur>(new SetUpSchurImpl(m_field));
}