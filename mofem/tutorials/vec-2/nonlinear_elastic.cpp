/**
 * \file nonlinear_elastic.cpp
 * \example nonlinear_elastic.cpp
 *
 * Plane stress elastic dynamic problem
 *
 */

#include <MoFEM.hpp>
#include <MatrixFunction.hpp>

using namespace MoFEM;

constexpr int SPACE_DIM =
    EXECUTABLE_DIMENSION; //< Space dimension of problem, mesh

constexpr AssemblyType A = AssemblyType::PETSC; //< selected assembly type

using EntData = EntitiesFieldData::EntData;
using DomainEle = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::DomainEle;
using BoundaryEle =
    PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::BoundaryEle;
using DomainEleOp = DomainEle::UserDataOperator;
using BoundaryEleOp = BoundaryEle::UserDataOperator;

using DomainNaturalBC =
    NaturalBC<DomainEleOp>::Assembly<PETSC>::LinearForm<GAUSS>;
using OpBodyForce =
    DomainNaturalBC::OpFlux<NaturalMeshsetType<BLOCKSET>, 1, SPACE_DIM>;

using BoundaryNaturalBC =
    NaturalBC<BoundaryEleOp>::Assembly<PETSC>::LinearForm<GAUSS>;
using OpForce = BoundaryNaturalBC::OpFlux<NaturalForceMeshsets, 1, SPACE_DIM>;

template <int DIM> struct PostProcEleByDim;

template <> struct PostProcEleByDim<2> {
  using PostProcEleDomain = PostProcBrokenMeshInMoabBase<DomainEle>;
  using PostProcEleBdy = PostProcBrokenMeshInMoabBase<BoundaryEle>;
  using SideEle = PipelineManager::ElementsAndOpsByDim<2>::FaceSideEle;
};

template <> struct PostProcEleByDim<3> {
  using PostProcEleDomain = PostProcBrokenMeshInMoabBase<BoundaryEle>;
  using PostProcEleBdy = PostProcBrokenMeshInMoabBase<BoundaryEle>;
  using SideEle = PipelineManager::ElementsAndOpsByDim<3>::FaceSideEle;
};

using SideEle = PostProcEleByDim<SPACE_DIM>::SideEle;
using PostProcEleBdy = PostProcEleByDim<SPACE_DIM>::PostProcEleBdy;
using PostProcEdges =
    PostProcBrokenMeshInMoabBase<EdgeElementForcesAndSourcesCore>;

#include <HenckyOps.hpp>
using namespace HenckyOps;
#include <RigidBodyTieConstraint.hpp>

struct Example {

  Example(MoFEM::Interface &m_field) : mField(m_field) {}

  MoFEMErrorCode runProblem();

private:
  MoFEM::Interface &mField;

  MoFEMErrorCode readMesh();
  MoFEMErrorCode setupProblem();
  MoFEMErrorCode boundaryCondition();
  MoFEMErrorCode assembleSystem();
  MoFEMErrorCode solveSystem();
  MoFEMErrorCode gettingNorms();
  MoFEMErrorCode outputResults();
  MoFEMErrorCode checkResults();

  struct TieBlock {
    Range tieFaces; // tie faces
    FTensor::Tensor1<double, 3> tieCoord;
    FTensor::Tensor1<double, 3> tieDirection;
    EntityHandle refPoint;
  };

  std::vector<TieBlock> tieBlocks; //< Store infomation about tie blocks
};

//! [Run problem]
MoFEMErrorCode Example::runProblem() {
  MoFEMFunctionBegin;
  CHKERR readMesh();
  CHKERR setupProblem();
  CHKERR boundaryCondition();
  CHKERR assembleSystem();
  CHKERR solveSystem();
  CHKERR gettingNorms();
  CHKERR outputResults();
  CHKERR checkResults();
  MoFEMFunctionReturn(0);
}
//! [Run problem]

//! [Read mesh]
MoFEMErrorCode Example::readMesh() {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  CHKERR simple->getOptions();
  char meshFileName[255];
  CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-file_name",
                               meshFileName, 255, PETSC_NULL);
  CHKERR simple->loadFile("", meshFileName,
                          EssentialPreProc<MPCsType>::loadFileWithMPCs);
  Range tie_ents;
  for (auto m :
       mField.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

           (boost::format("%s(.*)") % "TIE_MATRIX").str()

               ))

  ) {
    auto meshset = m->getMeshset();
    Range tie_meshset_range;
    CHKERR mField.get_moab().get_entities_by_dimension(meshset, SPACE_DIM - 1,
                                                       tie_meshset_range, true);
    std::vector<double> attributes;
    CHKERR m->getAttributes(attributes);
    if (attributes.size() != 6) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
               "Wrong number of head parameters %d", attributes.size());
    }
    tieBlocks.push_back({tie_meshset_range,
                         FTensor::Tensor1<double, 3>(
                             attributes[0], attributes[1], attributes[2]),
                         FTensor::Tensor1<double, 3>(
                             attributes[3], attributes[4], attributes[5])});
  }
  MoFEMFunctionReturn(0);
}
//! [Read mesh]

//! [Set up problem]
MoFEMErrorCode Example::setupProblem() {
  MoFEMFunctionBegin;
  Simple *simple = mField.getInterface<Simple>();

  enum bases { AINSWORTH, DEMKOWICZ, LASBASETOPT };
  const char *list_bases[LASBASETOPT] = {"ainsworth", "demkowicz"};
  PetscInt choice_base_value = AINSWORTH;
  CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list_bases,
                              LASBASETOPT, &choice_base_value, PETSC_NULL);

  FieldApproximationBase base;
  switch (choice_base_value) {
  case AINSWORTH:
    base = AINSWORTH_LEGENDRE_BASE;
    MOFEM_LOG("WORLD", Sev::inform)
        << "Set AINSWORTH_LEGENDRE_BASE for displacements";
    break;
  case DEMKOWICZ:
    base = DEMKOWICZ_JACOBI_BASE;
    MOFEM_LOG("WORLD", Sev::inform)
        << "Set DEMKOWICZ_JACOBI_BASE for displacements";
    break;
  default:
    base = LASTBASE;
    break;
  }

  // Add field
  CHKERR simple->addDomainField("U", H1, base, SPACE_DIM);
  CHKERR simple->addBoundaryField("U", H1, base, SPACE_DIM);
  int order = 2;
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order, PETSC_NULL);
  CHKERR simple->setFieldOrder("U", order);

  Range rigid_body_surface_ents;
  auto add_tie_lagrange_multiplier = [&]() {
    MoFEMFunctionBegin;
    CHKERR simple->addBoundaryField("LAMBDA", H1, base, SPACE_DIM);
    CHKERR simple->addMeshsetField("RIGID_BODY_LAMBDA", NOFIELD, NOBASE, 3);
    CHKERR simple->addMeshsetField("RIGID_BODY_THETA", NOFIELD, NOBASE, 3);
    CHKERR simple->addMeshsetField("LAMBDA", H1, base, SPACE_DIM);
    for (auto &t : tieBlocks) {
      // add tie nodes to rigid body meshset
      Range tie_nodes;
      CHKERR mField.get_moab().get_adjacencies(t.tieFaces, 0, false, tie_nodes,
                                               moab::Interface::UNION);
      // add tie edges to rigid body meshset
      Range tie_edges;
      CHKERR mField.get_moab().get_adjacencies(t.tieFaces, 1, false, tie_edges,
                                               moab::Interface::UNION);
      // add to range of rigid body entities
      rigid_body_surface_ents.merge(tie_nodes);
      rigid_body_surface_ents.merge(tie_edges);
      rigid_body_surface_ents.merge(t.tieFaces);
    }
    simple->getMeshsetFiniteElementEntities().push_back(
        rigid_body_surface_ents);

    CHKERR simple->setFieldOrder("LAMBDA", 0);
    for (auto &t : tieBlocks) {
      CHKERR simple->setFieldOrder("LAMBDA", order, &t.tieFaces);
    }
    MoFEMFunctionReturn(0);
  };

  CHKERR add_tie_lagrange_multiplier();

  CHKERR simple->setUp();
  MoFEMFunctionReturn(0);
}
//! [Set up problem]

//! [Boundary condition]
MoFEMErrorCode Example::boundaryCondition() {
  MoFEMFunctionBegin;
  auto *pipeline_mng = mField.getInterface<PipelineManager>();
  auto simple = mField.getInterface<Simple>();
  auto bc_mng = mField.getInterface<BcManager>();
  auto time_scale = boost::make_shared<TimeScale>();

  auto integration_rule = [](int, int, int approx_order) {
    return 2 * (approx_order - 1);
  };

  CHKERR pipeline_mng->setDomainRhsIntegrationRule(integration_rule);
  CHKERR pipeline_mng->setDomainLhsIntegrationRule(integration_rule);
  CHKERR pipeline_mng->setBoundaryRhsIntegrationRule(integration_rule);

  CHKERR BoundaryNaturalBC::AddFluxToPipeline<OpForce>::add(
      pipeline_mng->getOpBoundaryRhsPipeline(), mField, "U", {time_scale},
      "FORCE", Sev::inform);

  //! [Define gravity vector]
  CHKERR DomainNaturalBC::AddFluxToPipeline<OpBodyForce>::add(
      pipeline_mng->getOpDomainRhsPipeline(), mField, "U", {time_scale},
      "BODY_FORCE", Sev::inform);

  // Essential BC
  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "REMOVE_X",
                                           "U", 0, 0);
  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "REMOVE_Y",
                                           "U", 1, 1);
  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "REMOVE_Z",
                                           "U", 2, 2);
  CHKERR bc_mng->pushMarkDOFsOnEntities<DisplacementCubitBcData>(
      simple->getProblemName(), "U");

  // Essential MPCs BC
  CHKERR bc_mng->addBlockDOFsToMPCs(simple->getProblemName(), "U");

  auto add_constrain_lhs = [&](auto &pip) {
    MoFEMFunctionBegin;
    //CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(pip, {NOSPACE});
    auto u_ptr = boost::make_shared<MatrixDouble>();
    auto lambda_ptr = boost::make_shared<MatrixDouble>();
    auto translation_ptr = boost::make_shared<VectorDouble>();
    auto theta_ptr = boost::make_shared<VectorDouble>();

    pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));
    pip.push_back(
        new OpCalculateVectorFieldValues<SPACE_DIM>("LAMBDA", lambda_ptr));
    pip.push_back(new OpCalculateNoFieldVectorValues("RIGID_BODY_LAMBDA",
                                                     translation_ptr));
    pip.push_back(
        new OpCalculateNoFieldVectorValues("RIGID_BODY_THETA", theta_ptr));

    for (auto &t : tieBlocks) {
      pip.push_back(new OpTieTermConstrainRigidBodyLhs_dU<SPACE_DIM>(
          "LAMBDA", "U", u_ptr, t.tieCoord,
          boost::make_shared<Range>(t.tieFaces), translation_ptr, theta_ptr));
      pip.push_back(new OpTieTermConstrainRigidBodyLhs_dTranslation<SPACE_DIM>(
          "LAMBDA", "RIGID_BODY_LAMBDA",
          boost::make_shared<Range>(t.tieFaces)));
      pip.push_back(new OpTieTermConstrainRigidBodyLhs_dRotation<SPACE_DIM>(
          "LAMBDA", "RIGID_BODY_THETA", t.tieCoord,
          boost::make_shared<Range>(t.tieFaces)));
    }
    MoFEMFunctionReturn(0);
  };

  auto add_constrain_rhs = [&](auto &pip) {
    MoFEMFunctionBegin;
    auto u_ptr = boost::make_shared<MatrixDouble>();
    auto lambda_ptr = boost::make_shared<MatrixDouble>();
    auto translation_ptr = boost::make_shared<VectorDouble>();
    auto theta_ptr = boost::make_shared<VectorDouble>();

    pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));
    pip.push_back(
        new OpCalculateVectorFieldValues<SPACE_DIM>("LAMBDA", lambda_ptr));

    pip.push_back(new OpCalculateNoFieldVectorValues("RIGID_BODY_LAMBDA",
                                                     translation_ptr));
    pip.push_back(
        new OpCalculateNoFieldVectorValues("RIGID_BODY_THETA", theta_ptr));

    for (auto &t : tieBlocks) {
      pip.push_back(new OpTieTermConstrainRigidBodyRhs<SPACE_DIM>(
          "LAMBDA", u_ptr, t.tieCoord,
          boost::make_shared<Range>(t.tieFaces), translation_ptr, theta_ptr));
    }
    MoFEMFunctionReturn(0);
  };
  

  CHKERR add_constrain_lhs(pipeline_mng->getOpBoundaryLhsPipeline());
  CHKERR add_constrain_rhs(pipeline_mng->getOpBoundaryRhsPipeline());

  auto add_constrain_force_term = [&](auto &pip) {
    MoFEMFunctionBegin;
    auto lambda_ptr = boost::make_shared<MatrixDouble>();
    auto u_ptr = boost::make_shared<MatrixDouble>();
    pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));
    pip.push_back(
        new OpCalculateVectorFieldValues<SPACE_DIM>("LAMBDA", lambda_ptr));

    for (auto &t : tieBlocks) {
      pip.push_back(new OpTieTermConstrainRigidBodyRhs_du<SPACE_DIM>(
          "U", lambda_ptr, boost::make_shared<Range>(t.tieFaces)));
    }
    MoFEMFunctionReturn(0);
  };

  CHKERR add_constrain_force_term(pipeline_mng->getOpBoundaryRhsPipeline());

  auto add_rigid_body_rhs = [&](auto &pip) {
    MoFEMFunctionBegin;

    // add ops
    auto lambda_ptr = boost::make_shared<MatrixDouble>();
    auto translation_ptr = boost::make_shared<VectorDouble>();
    auto u_ptr = boost::make_shared<MatrixDouble>();
    auto int_translation_ptr =
        boost::make_shared<VectorDouble>();
    auto int_rotation_ptr =
        boost::make_shared<MatrixDouble>();

    Range rigid_body_ents;
    auto v_rigid_body_ents = simple->getMeshsetFiniteElementEntities();

    // loop over all rigid body entities
    for (auto ents : v_rigid_body_ents) {
      rigid_body_ents.merge(ents);
    }

    auto rigid_body_ents_ptr = boost::make_shared<Range>(rigid_body_ents);

    for (auto &t : tieBlocks) {
      auto op_loop_side = new OpLoopSide<BoundaryEle>(
          mField, "bFE", SPACE_DIM - 1, rigid_body_ents_ptr);
      CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
          op_loop_side->getOpPtrVector(), {});
      op_loop_side->getOpPtrVector().push_back(
          new OpCalculateVectorFieldValues<SPACE_DIM>("LAMBDA", lambda_ptr));
      op_loop_side->getOpPtrVector().push_back(
          new OpTieTermConstrainRigidBodyGlobalTranslationIntegralRhs<
              SPACE_DIM>("LAMBDA", boost::make_shared<Range>(t.tieFaces),
                         lambda_ptr, int_translation_ptr));
      op_loop_side->getOpPtrVector().push_back(
          new OpTieTermConstrainRigidBodyGlobalRotationIntegralRhs<SPACE_DIM>(
              "LAMBDA", t.tieCoord, lambda_ptr, int_rotation_ptr));
      pip.push_back(op_loop_side);

      pip.push_back(
          new OpTieTermConstrainRigidBodyGlobalTranslationRhs<SPACE_DIM>(
              "RIGID_BODY_LAMBDA", int_translation_ptr));
      pip.push_back(new OpTieTermConstrainRigidBodyGlobalRotationRhs<SPACE_DIM>(
          "RIGID_BODY_THETA", int_rotation_ptr));
    }
    MoFEMFunctionReturn(0);
  };

  // auto add_rigid_body_lhs = [&](auto &pip) {
  //   MoFEMFunctionBegin;

  //   // add ops
  //   auto lambda_ptr = boost::make_shared<MatrixDouble>();
  //   auto translation_ptr = boost::make_shared<VectorDouble>();
  //   auto u_ptr = boost::make_shared<MatrixDouble>();

  //   Range rigid_body_ents;
  //   auto v_rigid_body_ents = simple->getMeshsetFiniteElementEntities();

  //   // loop over all rigid body entities
  //   for (auto ents : v_rigid_body_ents) {
  //     rigid_body_ents.merge(ents);
  //   }

  //   auto rigid_body_ents_ptr = boost::make_shared<Range>(rigid_body_ents);

  //   for (auto &t : tieBlocks) {
  //     auto op_loop_side = new OpLoopSide<BoundaryEle>(
  //         mField, "bFE", SPACE_DIM - 1, rigid_body_ents_ptr);
  //     CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
  //         op_loop_side->getOpPtrVector(), {});
  //     op_loop_side->getOpPtrVector().push_back(
  //         new OpCalculateVectorFieldValues<SPACE_DIM>("LAMBDA", lambda_ptr));
  //     // push operator to gather col_data?

  //     pip.push_back(op_loop_side);

  //     pip.push_back(new OpTieTermConstrainTranslation_dlambda<SPACE_DIM>(
  //         "RIGID_BODY_LAMBDA", "LAMBDA", u_ptr, t.tieCoord, t.tieDirection,
  //         boost::make_shared<Range>(t.tieFaces)));
  //     pip.push_back(new OpTieTermConstrainRotation_dlambda<SPACE_DIM>(
  //         "RIGID_BODY_THETA", "LAMBDA", u_ptr, t.tieCoord, t.tieDirection,
  //         boost::make_shared<Range>(t.tieFaces)));
  //   }
  //   MoFEMFunctionReturn(0);
  // };

  CHKERR add_rigid_body_rhs(pipeline_mng->getOpMeshsetRhsPipeline());


  MoFEMFunctionReturn(0);
}
//! [Boundary condition]

//! [Push operators to pipeline]
MoFEMErrorCode Example::assembleSystem() {
  MoFEMFunctionBegin;
  auto *pipeline_mng = mField.getInterface<PipelineManager>();

  auto add_domain_ops_lhs = [&](auto &pip) {
    MoFEMFunctionBegin;
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1});
    CHKERR HenckyOps::opFactoryDomainLhs<SPACE_DIM, PETSC, GAUSS, DomainEleOp>(
        mField, pip, "U", "MAT_ELASTIC", Sev::inform);
    MoFEMFunctionReturn(0);
  };

  auto add_domain_ops_rhs = [&](auto &pip) {
    MoFEMFunctionBegin;
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1});
    CHKERR HenckyOps::opFactoryDomainRhs<SPACE_DIM, PETSC, GAUSS, DomainEleOp>(
        mField, pip, "U", "MAT_ELASTIC", Sev::inform);
    MoFEMFunctionReturn(0);
  };

  CHKERR add_domain_ops_lhs(pipeline_mng->getOpDomainLhsPipeline());
  CHKERR add_domain_ops_rhs(pipeline_mng->getOpDomainRhsPipeline());

  MoFEMFunctionReturn(0);
}
//! [Push operators to pipeline]

/**
 * @brief Monitor solution
 *
 * This functions is called by TS solver at the end of each step. It is used
 * to output results to the hard drive.
 */
struct Monitor : public FEMethod {
  Monitor(SmartPetscObj<DM> dm, boost::shared_ptr<PostProcEleBdy> post_proc)
      : dM(dm), postProc(post_proc){};
  MoFEMErrorCode postProcess() {
    MoFEMFunctionBegin;
    constexpr int save_every_nth_step = 1;
    if (ts_step % save_every_nth_step == 0) {
      CHKERR DMoFEMLoopFiniteElements(dM, "bFE", postProc);
      CHKERR EssentialPreProc<MPCsType>::addLinksToPostProcMesh(
          postProc->mField, postProc->getPostProcMesh(), {"U"});

      CHKERR postProc->writeFile(
          "out_step_" + boost::lexical_cast<std::string>(ts_step) + ".h5m");
    }
    MoFEMFunctionReturn(0);
  }

private:
  SmartPetscObj<DM> dM;
  boost::shared_ptr<PostProcEleBdy> postProc;
};

//! [Solve]
MoFEMErrorCode Example::solveSystem() {
  MoFEMFunctionBegin;
  auto *simple = mField.getInterface<Simple>();
  auto *pipeline_mng = mField.getInterface<PipelineManager>();

  auto dm = simple->getDM();
  auto ts = pipeline_mng->createTSIM();

  // Setup postprocessing
  auto create_post_proc_fe = [dm, this, simple]() {
    auto post_proc_ele_domain = [dm, this](auto &pip_domain) {
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip_domain, {H1});
      auto common_ptr = commonDataFactory<SPACE_DIM, GAUSS, DomainEleOp>(
          mField, pip_domain, "U", "MAT_ELASTIC", Sev::inform);
      return common_ptr;
    };

    auto post_proc_fe_bdy = boost::make_shared<PostProcEleBdy>(mField);
    auto u_ptr = boost::make_shared<MatrixDouble>();
    post_proc_fe_bdy->getOpPtrVector().push_back(
        new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));
    auto op_loop_side =
        new OpLoopSide<SideEle>(mField, simple->getDomainFEName(), SPACE_DIM);
    auto common_ptr = post_proc_ele_domain(op_loop_side->getOpPtrVector());
    post_proc_fe_bdy->getOpPtrVector().push_back(op_loop_side);

    using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

    post_proc_fe_bdy->getOpPtrVector().push_back(

        new OpPPMap(

            post_proc_fe_bdy->getPostProcMesh(),
            post_proc_fe_bdy->getMapGaussPts(),

            {},

            {{"U", u_ptr}},

            {{"GRAD", common_ptr->matGradPtr},
             {"FIRST_PIOLA", common_ptr->getMatFirstPiolaStress()}},

            {}

            )

    );
    return post_proc_fe_bdy;
  };

  auto add_extra_finite_elements_to_ksp_solver_pipelines = [&]() {
    MoFEMFunctionBegin;

    auto pre_proc_ptr = boost::make_shared<FEMethod>();
    auto post_proc_rhs_ptr = boost::make_shared<FEMethod>();
    auto post_proc_lhs_ptr = boost::make_shared<FEMethod>();

    auto time_scale = boost::make_shared<TimeScale>();

    auto get_bc_hook_rhs = [this, pre_proc_ptr, time_scale]() {
      MoFEMFunctionBeginHot;
      CHKERR EssentialPreProc<DisplacementCubitBcData>(mField, pre_proc_ptr,
                                                       {time_scale}, false)();
      MoFEMFunctionReturnHot(0);
    };

    pre_proc_ptr->preProcessHook = get_bc_hook_rhs;

    auto get_post_proc_hook_rhs = [this, post_proc_rhs_ptr]() {
      MoFEMFunctionBegin;
      CHKERR EssentialPreProcReaction<DisplacementCubitBcData>(
          mField, post_proc_rhs_ptr, nullptr, Sev::verbose)();
      CHKERR EssentialPostProcRhs<DisplacementCubitBcData>(
          mField, post_proc_rhs_ptr, 1.)();
      CHKERR EssentialPostProcRhs<MPCsType>(mField, post_proc_rhs_ptr)();
      MoFEMFunctionReturn(0);
    };
    auto get_post_proc_hook_lhs = [this, post_proc_lhs_ptr]() {
      MoFEMFunctionBegin;
      CHKERR EssentialPostProcLhs<DisplacementCubitBcData>(
          mField, post_proc_lhs_ptr, 1.)();
      CHKERR EssentialPostProcLhs<MPCsType>(mField, post_proc_lhs_ptr)();
      MoFEMFunctionReturn(0);
    };
    post_proc_rhs_ptr->postProcessHook = get_post_proc_hook_rhs;
    post_proc_lhs_ptr->postProcessHook = get_post_proc_hook_lhs;

    auto set_global_bc = [&]() {
      MoFEMFunctionBegin;
      auto ksp_ctx_ptr = getDMKspCtx(dm);
      auto ts_ctx_ptr = getDMTsCtx(dm);
      auto fe_post_proc_rhs = boost::make_shared<FEMethod>();
      auto fe_post_proc_lhs = boost::make_shared<FEMethod>();

      auto set_post_rhs = [this, fe_post_proc_rhs, ts_ctx_ptr]() {
        MoFEMFunctionBeginHot;
        auto is_mng = mField.getInterface<ISManager>();
        auto simple = mField.getInterface<Simple>();
        SmartPetscObj<IS> is;
        CHKERR is_mng->isCreateProblemFieldAndRankLocal(
            simple->getProblemName(), ROW, "RIGID_BODY_LAMBDA", 0,
            MAX_DOFS_ON_ENTITY, is, nullptr);

        const int *i_ptr;
        CHKERR ISGetIndices(is, &i_ptr);
        int size;
        CHKERR ISGetLocalSize(is, &size);
        double *f_ptr;
        double *x_ptr;
        CHKERR VecGetArray(fe_post_proc_rhs->f, &f_ptr);
        CHKERR VecGetArray(fe_post_proc_rhs->x, &x_ptr);
        for (auto i = 0; i != size; ++i) {
          if (i % SPACE_DIM == 0) {
            // FIXME: this is a hack, we should use a proper way to set the tieBlock
            f_ptr[i_ptr[i]] = x_ptr[i_ptr[i]] - (tieBlocks[0].tieDirection(0) * fe_post_proc_rhs->ts_t);
          }
        }
        CHKERR VecRestoreArray(fe_post_proc_rhs->f, &f_ptr);
        CHKERR VecRestoreArray(fe_post_proc_rhs->x, &x_ptr);
        CHKERR ISRestoreIndices(is, &i_ptr);
        MoFEMFunctionReturnHot(0);
      };

      auto set_post_lhs = [this, fe_post_proc_lhs]() {
        MoFEMFunctionBeginHot;
        auto is_mng = mField.getInterface<ISManager>();
        auto simple = mField.getInterface<Simple>();
        SmartPetscObj<IS> is;
        CHKERR is_mng->isCreateProblemFieldAndRankLocal(
            simple->getProblemName(), ROW, "RIGID_BODY_LAMBDA", 0, SPACE_DIM,
            is, nullptr);
        *(fe_post_proc_lhs->matAssembleSwitch) = PETSC_TRUE;
        CHKERR MatAssemblyBegin(fe_post_proc_lhs->B, MAT_FINAL_ASSEMBLY);
        CHKERR MatAssemblyEnd(fe_post_proc_lhs->B, MAT_FINAL_ASSEMBLY);
        int is_size;
        CHKERR ISGetSize(is, &is_size);
        if(is_size != SPACE_DIM)
          SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY, "Wrong size of rigid translation dofs");
        const int *i_ptr;
        CHKERR ISGetIndices(is, &i_ptr);
        const int idx = 0;
        CHKERR MatSetOption(fe_post_proc_lhs->B,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
        CHKERR MatZeroRows(fe_post_proc_lhs->B, 1, &i_ptr[idx], 1, PETSC_NULL,
                           PETSC_NULL);
        CHKERR ISRestoreIndices(is, &i_ptr);

        MoFEMFunctionReturnHot(0);
      };

      fe_post_proc_rhs->postProcessHook = set_post_rhs;
      ts_ctx_ptr->getPostProcessIFunction().push_back(fe_post_proc_rhs);
      fe_post_proc_lhs->postProcessHook = set_post_lhs;
      ts_ctx_ptr->getPostProcessIJacobian().push_back(fe_post_proc_lhs);
      MoFEMFunctionReturn(0);
    };

    CHKERR set_global_bc();

    // This is low level pushing finite elements (pipelines) to solver
    auto ts_ctx_ptr = getDMTsCtx(simple->getDM());
    ts_ctx_ptr->getPreProcessIFunction().push_front(pre_proc_ptr);
    ts_ctx_ptr->getPreProcessIJacobian().push_front(pre_proc_ptr);
    ts_ctx_ptr->getPostProcessIFunction().push_back(post_proc_rhs_ptr);
    ts_ctx_ptr->getPostProcessIJacobian().push_back(post_proc_lhs_ptr);
    MoFEMFunctionReturn(0);
  };

  // Add extra finite elements to SNES solver pipelines to resolve essential
  // boundary conditions
  CHKERR add_extra_finite_elements_to_ksp_solver_pipelines();

  auto create_monitor_fe = [dm](auto &&post_proc_fe) {
    return boost::make_shared<Monitor>(dm, post_proc_fe);
  };

  // Set monitor which postprocessing results and saves them to the hard drive
  boost::shared_ptr<FEMethod> null_fe;
  auto monitor_ptr = create_monitor_fe(create_post_proc_fe());
  CHKERR DMMoFEMTSSetMonitor(dm, ts, simple->getDomainFEName(), null_fe,
                             null_fe, monitor_ptr);

  // Set time solver
  double ftime = 1;
  CHKERR TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);

  auto D = createDMVector(simple->getDM());

  CHKERR TSSetSolution(ts, D);
  CHKERR TSSetFromOptions(ts);

  CHKERR TSSolve(ts, NULL);
  CHKERR TSGetTime(ts, &ftime);

  PetscInt steps, snesfails, rejects, nonlinits, linits;
  CHKERR TSGetStepNumber(ts, &steps);
  CHKERR TSGetSNESFailures(ts, &snesfails);
  CHKERR TSGetStepRejections(ts, &rejects);
  CHKERR TSGetSNESIterations(ts, &nonlinits);
  CHKERR TSGetKSPIterations(ts, &linits);
  MOFEM_LOG_C("EXAMPLE", Sev::inform,
              "steps %d (%d rejected, %d SNES fails), ftime %g, nonlinits "
              "%d, linits %d",
              steps, rejects, snesfails, ftime, nonlinits, linits);

  MoFEMFunctionReturn(0);
}
//! [Solve]

//! [Getting norms]
MoFEMErrorCode Example::gettingNorms() {
  MoFEMFunctionBegin;

  auto simple = mField.getInterface<Simple>();
  auto dm = simple->getDM();

  auto T = createDMVector(simple->getDM());
  CHKERR DMoFEMMeshToLocalVector(simple->getDM(), T, INSERT_VALUES,
                                 SCATTER_FORWARD);
  double nrm2;
  CHKERR VecNorm(T, NORM_2, &nrm2);
  MOFEM_LOG("EXAMPLE", Sev::inform) << "Solution norm " << nrm2;

  auto post_proc_norm_fe = boost::make_shared<DomainEle>(mField);

  auto post_proc_norm_rule_hook = [](int, int, int p) -> int { return 2 * p; };
  post_proc_norm_fe->getRuleHook = post_proc_norm_rule_hook;

  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      post_proc_norm_fe->getOpPtrVector(), {H1});

  enum NORMS { U_NORM_L2 = 0, PIOLA_NORM, LAST_NORM };
  auto norms_vec =
      createVectorMPI(mField.get_comm(),
                      (mField.get_comm_rank() == 0) ? LAST_NORM : 0, LAST_NORM);
  CHKERR VecZeroEntries(norms_vec);

  auto u_ptr = boost::make_shared<MatrixDouble>();
  post_proc_norm_fe->getOpPtrVector().push_back(
      new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));

  post_proc_norm_fe->getOpPtrVector().push_back(
      new OpCalcNormL2Tensor1<SPACE_DIM>(u_ptr, norms_vec, U_NORM_L2));

  auto common_ptr = commonDataFactory<SPACE_DIM, GAUSS, DomainEleOp>(
      mField, post_proc_norm_fe->getOpPtrVector(), "U", "MAT_ELASTIC",
      Sev::inform);

  post_proc_norm_fe->getOpPtrVector().push_back(
      new OpCalcNormL2Tensor2<SPACE_DIM, SPACE_DIM>(
          common_ptr->getMatFirstPiolaStress(), norms_vec, PIOLA_NORM));

  CHKERR DMoFEMLoopFiniteElements(dm, simple->getDomainFEName(),
                                  post_proc_norm_fe);

  CHKERR VecAssemblyBegin(norms_vec);
  CHKERR VecAssemblyEnd(norms_vec);

  MOFEM_LOG_CHANNEL("SELF"); // Clear channel from old tags
  if (mField.get_comm_rank() == 0) {
    const double *norms;
    CHKERR VecGetArrayRead(norms_vec, &norms);
    MOFEM_TAG_AND_LOG("SELF", Sev::inform, "example")
        << "norm_u: " << std::scientific << std::sqrt(norms[U_NORM_L2]);
    MOFEM_TAG_AND_LOG("SELF", Sev::inform, "example")
        << "norm_piola: " << std::scientific << std::sqrt(norms[PIOLA_NORM]);
    CHKERR VecRestoreArrayRead(norms_vec, &norms);
  }

  MoFEMFunctionReturn(0);
}
//! [Getting norms]

//! [Postprocessing results]
MoFEMErrorCode Example::outputResults() {
  MoFEMFunctionBegin;
  PetscInt test_nb = 0;
  PetscBool test_flg = PETSC_FALSE;
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-test", &test_nb, &test_flg);

  if (test_flg) {
    auto simple = mField.getInterface<Simple>();
    auto T = createDMVector(simple->getDM());
    CHKERR DMoFEMMeshToLocalVector(simple->getDM(), T, INSERT_VALUES,
                                   SCATTER_FORWARD);
    double nrm2;
    CHKERR VecNorm(T, NORM_2, &nrm2);
    MOFEM_LOG("EXAMPLE", Sev::verbose) << "Regression norm " << nrm2;
    double regression_value = 0;
    switch (test_nb) {
    case 1:
      regression_value = 1.02789;
      break;
    case 2:
      regression_value = 1.8841e+00;
      break;
    case 3:
      regression_value = 1.8841e+00;
      break;
    case 4: // just links
      regression_value = 4.9625e+00;
      break;
    case 5: // link master 
      regression_value = 6.6394e+00;
      break;
    case 6: // link master swap
      regression_value = 4.98764e+00;
      break;
    case 7: // link Y
      regression_value = 4.9473e+00;
      break;
    case 8: // link_3D_repr
      regression_value = 2.5749e-01;
      break;

    default:
      SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "Wrong test number.");
      break;
    }
    if (fabs(nrm2 - regression_value) > 1e-2)
      SETERRQ2(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
               "Regression test field; wrong norm value. %6.4e != %6.4e", nrm2,
               regression_value);
  }
  MoFEMFunctionReturn(0);
}
//! [Postprocessing results]

//! [Check]
MoFEMErrorCode Example::checkResults() {
  MoFEMFunctionBegin;
  MoFEMFunctionReturn(0);
}
//! [Check]

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  // Initialisation of MoFEM/PETSc and MOAB data structures
  const char param_file[] = "param_file.petsc";
  MoFEM::Core::Initialize(&argc, &argv, param_file, help);

  // Add logging channel for example
  auto core_log = logging::core::get();
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmWorld(), "EXAMPLE"));
  LogManager::setLog("EXAMPLE");
  MOFEM_LOG_TAG("EXAMPLE", "example");

  try {

    //! [Register MoFEM discrete manager in PETSc]
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);
    //! [Register MoFEM discrete manager in PETSc

    //! [Create MoAB]
    moab::Core mb_instance;              ///< mesh database
    moab::Interface &moab = mb_instance; ///< mesh database interface
    //! [Create MoAB]

    //! [Create MoFEM]
    MoFEM::Core core(moab);           ///< finite element database
    MoFEM::Interface &m_field = core; ///< finite element database interface
    //! [Create MoFEM]

    //! [Example]
    Example ex(m_field);
    CHKERR ex.runProblem();
    //! [Example]
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
