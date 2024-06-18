/**
 * \file PostProcContact.hpp
 *
 *
 * @copyright Copyright (c) 2023
 */

namespace ContactOps {

template <int DIM> struct PostProcEleByDim;

template <> struct PostProcEleByDim<2> {
  using PostProcEleDomain = PostProcBrokenMeshInMoabBaseCont<DomainEle>;
  using PostProcEleBdy = PostProcBrokenMeshInMoabBaseCont<BoundaryEle>;
  using SideEle = PipelineManager::ElementsAndOpsByDim<2>::FaceSideEle;
};

template <> struct PostProcEleByDim<3> {
  using PostProcEleDomain = PostProcBrokenMeshInMoabBaseCont<BoundaryEle>;
  using PostProcEleBdy = PostProcBrokenMeshInMoabBaseCont<BoundaryEle>;
  using SideEle = PipelineManager::ElementsAndOpsByDim<3>::FaceSideEle;
};

using PostProcEleDomain = PostProcEleByDim<SPACE_DIM>::PostProcEleDomain;
using SideEle = PostProcEleByDim<SPACE_DIM>::SideEle;
using PostProcEleBdy = PostProcEleByDim<SPACE_DIM>::PostProcEleBdy;

struct Monitor : public FEMethod {

  Monitor(SmartPetscObj<DM> &dm, double scale,
          boost::shared_ptr<GenericElementInterface> mfront_interface = nullptr,
          bool is_axisymmetric = false)
      : dM(dm), moabVertex(mbVertexPostproc), sTEP(0),
        mfrontInterface(mfront_interface) {

    MoFEM::Interface *m_field_ptr;
    CHKERR DMoFEMGetInterfacePtr(dM, &m_field_ptr);

    using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

    struct OpScale : public ForcesAndSourcesCore::UserDataOperator {
      OpScale(boost::shared_ptr<MatrixDouble> m_ptr, double s)
          : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPSPACE),
            mPtr(m_ptr), scale(s) {}
      MoFEMErrorCode doWork(int, EntityType, EntitiesFieldData::EntData &) {
        *mPtr *= 1./scale;
        return 0;
      }

    private:
      boost::shared_ptr<MatrixDouble> mPtr;
      double scale;
    };

    auto push_domain_ops = [&](auto &pip) {
      CHK_THROW_MESSAGE((AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
                            pip, {H1, HDIV}, "GEOMETRY")),
                        "Apply base transform");
      auto henky_common_data_ptr =
          commonDataFactory<SPACE_DIM, GAUSS, DomainEleOp>(
              *m_field_ptr, pip, "U", "MAT_ELASTIC", Sev::inform, scale);
      auto contact_stress_ptr = boost::make_shared<MatrixDouble>();
      pip.push_back(new OpCalculateHVecTensorField<SPACE_DIM, SPACE_DIM>(
          "SIGMA", contact_stress_ptr));
      pip.push_back(new OpScale(contact_stress_ptr, scale));
      return std::make_tuple(henky_common_data_ptr, contact_stress_ptr);
    };

    auto push_bdy_ops = [&](auto &pip) {
      // evaluate traction
      auto common_data_ptr = boost::make_shared<ContactOps::CommonData>();
      pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>(
          "U", common_data_ptr->contactDispPtr()));
      pip.push_back(new OpCalculateHVecTensorTrace<SPACE_DIM, BoundaryEleOp>(
          "SIGMA", common_data_ptr->contactTractionPtr()));
      pip.push_back(new OpScale(common_data_ptr->contactTractionPtr(), scale));
      using C = ContactIntegrators<BoundaryEleOp>;
      pip.push_back(new typename C::template OpEvaluateSDF<SPACE_DIM, GAUSS>(
          common_data_ptr));
      return common_data_ptr;
    };

    auto get_domain_pip = [&](auto &pip)
        -> boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> & {
      if constexpr (SPACE_DIM == 3) {
        auto op_loop_side = new OpLoopSide<SideEle>(
            *m_field_ptr, "dFE", SPACE_DIM, Sev::noisy,
            boost::make_shared<
                ForcesAndSourcesCore::UserDataOperator::AdjCache>());
        pip.push_back(op_loop_side);
        return op_loop_side->getOpPtrVector();
      } else {
        return pip;
      }
    };

    auto get_post_proc_domain_fe = [&]() {
      auto post_proc_fe =
          boost::make_shared<PostProcEleDomain>(*m_field_ptr, postProcMesh);
      auto &pip = post_proc_fe->getOpPtrVector();

      auto [henky_common_data_ptr, contact_stress_ptr] =
          push_domain_ops(get_domain_pip(pip));

      auto u_ptr = boost::make_shared<MatrixDouble>();
      pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));
      auto X_ptr = boost::make_shared<MatrixDouble>();
      pip.push_back(
          new OpCalculateVectorFieldValues<SPACE_DIM>("GEOMETRY", X_ptr));



      pip.push_back(

          new OpPPMap(

              post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

              {},
              {

                  {"U", u_ptr}, {"GEOMETRY", X_ptr}

              },
              {

                  {"SIGMA", contact_stress_ptr},

                  {"G", henky_common_data_ptr->matGradPtr},

                  {"PK1", henky_common_data_ptr->getMatFirstPiolaStress()}

              },
              {}

              )

      );

      if (SPACE_DIM == 3) {

        CHK_THROW_MESSAGE((AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
                              pip, {HDIV}, "GEOMETRY")),
                          "Apply transform");
        auto common_data_ptr = push_bdy_ops(pip);

        pip.push_back(

            new OpPPMap(

                post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

                {{"SDF", common_data_ptr->sdfPtr()},
                 {"CONSTRAINT_CONTACT", common_data_ptr->constraintPtr()}},

                {

                    {"TRACTION_CONTACT", common_data_ptr->contactTractionPtr()},
                    {"GRAD_SDF", common_data_ptr->gradSdfPtr()}

                },

                {},

                {{"HESS_SDF", common_data_ptr->hessSdfPtr()}}

                )

        );
      }

      return post_proc_fe;
    };

    auto get_post_proc_bdy_fe = [&]() {
      auto post_proc_fe =
          boost::make_shared<PostProcEleBdy>(*m_field_ptr, postProcMesh);
      auto &pip = post_proc_fe->getOpPtrVector();

      CHK_THROW_MESSAGE((AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
                            pip, {HDIV}, "GEOMETRY")),
                        "Apply transform");
      auto common_data_ptr = push_bdy_ops(pip);

      // create OP which run element on side
      auto op_loop_side = new OpLoopSide<SideEle>(
          *m_field_ptr, "dFE", SPACE_DIM, Sev::noisy,
          boost::make_shared<
              ForcesAndSourcesCore::UserDataOperator::AdjCache>());
      pip.push_back(op_loop_side);

      auto [henky_common_data_ptr, contact_stress_ptr] =
          push_domain_ops(op_loop_side->getOpPtrVector());

      auto X_ptr = boost::make_shared<MatrixDouble>();
      pip.push_back(
          new OpCalculateVectorFieldValues<SPACE_DIM>("GEOMETRY", X_ptr));

      pip.push_back(

          new OpPPMap(

              post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

              {{"SDF", common_data_ptr->sdfPtr()},
               {"CONSTRAINT_CONTACT", common_data_ptr->constraintPtr()}},

              {{"U", common_data_ptr->contactDispPtr()},
               {"GEOMETRY", X_ptr},
               {"TRACTION_CONTACT", common_data_ptr->contactTractionPtr()},
               {"GRAD_SDF", common_data_ptr->gradSdfPtr()}

              },

              {},

              {{"HESS_SDF", common_data_ptr->hessSdfPtr()}}

              )

      );

      return post_proc_fe;
    };

    auto get_integrate_traction = [&]() {
      auto integrate_traction = boost::make_shared<BoundaryEle>(*m_field_ptr);
      CHK_THROW_MESSAGE(
          (AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
              integrate_traction->getOpPtrVector(), {HDIV}, "GEOMETRY")),
          "Apply transform");
      // We have to integrate on curved face geometry, thus integration weight
      // have to adjusted.
      integrate_traction->getOpPtrVector().push_back(
          new OpSetHOWeightsOnSubDim<SPACE_DIM>());
      integrate_traction->getRuleHook = [](int, int, int approx_order) {
        return 2 * approx_order + geom_order - 1;
      };

      CHK_THROW_MESSAGE(
          (opFactoryCalculateTraction<SPACE_DIM, GAUSS, BoundaryEleOp>(
              integrate_traction->getOpPtrVector(), "SIGMA", is_axisymmetric)),
          "push operators to calculate traction");

      return integrate_traction;
    };

    auto get_integrate_area = [&]() {
      auto integrate_area = boost::make_shared<BoundaryEle>(*m_field_ptr);

      CHK_THROW_MESSAGE(
          (AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
              integrate_area->getOpPtrVector(), {HDIV}, "GEOMETRY")),
          "Apply transform");
      // We have to integrate on curved face geometry, thus integration weight
      // have to adjusted.
      integrate_area->getOpPtrVector().push_back(
          new OpSetHOWeightsOnSubDim<SPACE_DIM>()); // Ask Lukasz
      integrate_area->getRuleHook = [](int, int, int approx_order) {
        return 2 * approx_order + geom_order - 1;
      };
      Range contact_range;
      for (auto m :
           m_field_ptr->getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
               std::regex((boost::format("%s(.*)") % "CONTACT").str()))) {
        auto meshset = m->getMeshset();
        Range contact_meshset_range;
        CHKERR m_field_ptr->get_moab().get_entities_by_dimension(
            meshset, SPACE_DIM - 1, contact_meshset_range, true);

        CHKERR m_field_ptr->getInterface<CommInterface>()->synchroniseEntities(
            contact_meshset_range);
        contact_range.merge(contact_meshset_range);
      }

      auto contact_range_ptr = boost::make_shared<Range>(contact_range);

      auto op_loop_side = new OpLoopSide<SideEle>(
          *m_field_ptr, m_field_ptr->getInterface<Simple>()->getDomainFEName(),
          SPACE_DIM);
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
          op_loop_side->getOpPtrVector(), {H1}, "GEOMETRY");

      CHK_THROW_MESSAGE(
          (opFactoryCalculateArea<SPACE_DIM, GAUSS, BoundaryEleOp>(
              integrate_area->getOpPtrVector(), op_loop_side, "SIGMA", "U",
              is_axisymmetric, contact_range_ptr)),
          "push operators to calculate area");

      return integrate_area;
    };

    postProcDomainFe = get_post_proc_domain_fe();
    if constexpr (SPACE_DIM == 2)
      postProcBdyFe = get_post_proc_bdy_fe();

    integrateTraction = get_integrate_traction();
    integrateArea = get_integrate_area();

    normsVec = createVectorMPI(
        m_field_ptr->get_comm(),
        (m_field_ptr->get_comm_rank() == 0) ? LAST_NORM : 0, LAST_NORM);
  }

  MoFEMErrorCode preProcess() { return 0; }
  MoFEMErrorCode operator()() { return 0; }

  MoFEMErrorCode postProcess() {
    MoFEMFunctionBegin;
    MoFEM::Interface *m_field_ptr;
    CHKERR DMoFEMGetInterfacePtr(dM, &m_field_ptr);

    auto post_proc = [&]() {
      MoFEMFunctionBegin;

      if (!mfrontInterface) {
        auto post_proc_begin =
            boost::make_shared<PostProcBrokenMeshInMoabBaseBegin>(*m_field_ptr,
                                                                  postProcMesh);
        auto post_proc_end =
            boost::make_shared<PostProcBrokenMeshInMoabBaseEnd>(*m_field_ptr,
                                                                postProcMesh);

        CHKERR DMoFEMPreProcessFiniteElements(dM,
                                              post_proc_begin->getFEMethod());
        if (!postProcBdyFe) {
          postProcDomainFe->copyTs(*this); // this here is a Monitor
          CHKERR DMoFEMLoopFiniteElements(dM, "bFE", postProcDomainFe);
        } else {
          postProcDomainFe->copyTs(*this); // this here is a Monitor
          postProcBdyFe->copyTs(*this);
          CHKERR DMoFEMLoopFiniteElements(dM, "dFE", postProcDomainFe);
          CHKERR DMoFEMLoopFiniteElements(dM, "bFE", postProcBdyFe);
        }
        CHKERR DMoFEMPostProcessFiniteElements(dM,
                                               post_proc_end->getFEMethod());

        CHKERR post_proc_end->writeFile(
            "out_contact_" + boost::lexical_cast<std::string>(sTEP) + ".h5m");
      } else {
        CHKERR mfrontInterface->postProcessElement(
            ts_step, dM,
            m_field_ptr->getInterface<Simple>()->getDomainFEName());
      }

      MoFEMFunctionReturn(0);
    };

    auto calculate_force = [&] {
      MoFEMFunctionBegin;
      CHKERR VecZeroEntries(CommonData::totalTraction);
      CHKERR DMoFEMLoopFiniteElements(dM, "bFE", integrateTraction);
      CHKERR VecAssemblyBegin(CommonData::totalTraction);
      CHKERR VecAssemblyEnd(CommonData::totalTraction);
      MoFEMFunctionReturn(0);
    };

    auto calculate_area = [&] {
      MoFEMFunctionBegin;
      integrateArea->copyTs(*this);
      CHKERR DMoFEMLoopFiniteElements(dM, "bFE", integrateArea);
      CHKERR VecAssemblyBegin(CommonData::totalTraction);
      CHKERR VecAssemblyEnd(CommonData::totalTraction);
      MoFEMFunctionReturn(0);
    };

    auto calculate_reactions = [&]() {
      MoFEMFunctionBegin;

      auto res = createDMVector(dM);

      auto assemble_domain = [&]() {
        MoFEMFunctionBegin;
        auto fe_rhs = boost::make_shared<DomainEle>(*m_field_ptr);
        auto &pip = fe_rhs->getOpPtrVector();
        fe_rhs->f = res;

        auto integration_rule = [](int, int, int approx_order) {
          return 2 * approx_order + geom_order - 1;
        };
        fe_rhs->getRuleHook = integration_rule;

        CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1},
                                                              "GEOMETRY");
        CHKERR
        ContactOps::opFactoryDomainRhs<SPACE_DIM, PETSC, IT, DomainEleOp>(
            pip, "SIGMA", "U", is_axisymmetric);

        if (!mfrontInterface) {
          CHKERR
          HenckyOps::opFactoryDomainRhs<SPACE_DIM, PETSC, IT, DomainEleOp>(
              *m_field_ptr, pip, "U", "MAT_ELASTIC", Sev::inform, scale);
        } else {
          CHKERR mfrontInterface->opFactoryDomainRhs(pip);
        }
        CHKERR DMoFEMLoopFiniteElements(dM, "dFE", fe_rhs);

        CHKERR VecAssemblyBegin(res);
        CHKERR VecAssemblyEnd(res);
        CHKERR VecGhostUpdateBegin(res, ADD_VALUES, SCATTER_REVERSE);
        CHKERR VecGhostUpdateEnd(res, ADD_VALUES, SCATTER_REVERSE);

        MoFEMFunctionReturn(0);
      };

      auto assemble_boundary = [&]() {
        MoFEMFunctionBegin;
        auto fe_rhs = boost::make_shared<BoundaryEle>(*m_field_ptr);
        auto &pip = fe_rhs->getOpPtrVector();
        fe_rhs->f = res;

        auto integration_rule = [](int, int, int approx_order) {
          return 2 * approx_order + geom_order - 1;
        };
        fe_rhs->getRuleHook = integration_rule;

        CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(pip, {},
                                                                  "GEOMETRY");
        // We have to integrate on curved face geometry, thus integration weight
        // have to adjusted.
        pip.push_back(new OpSetHOWeightsOnSubDim<SPACE_DIM>());

        auto u_disp = boost::make_shared<MatrixDouble>();
        pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_disp));
        pip.push_back(
            new OpSpringRhs("U", u_disp, [this](double, double, double) {
              return spring_stiffness;
            }));

        CHKERR DMoFEMLoopFiniteElements(dM, "bFE", fe_rhs);

        MoFEMFunctionReturn(0);
      };

      CHKERR assemble_domain();
      CHKERR assemble_boundary();

      auto fe_post_proc_ptr = boost::make_shared<FEMethod>();
      auto get_post_proc_hook_rhs = [this, fe_post_proc_ptr, res,
                                     m_field_ptr]() {
        MoFEMFunctionBegin;
        CHKERR EssentialPreProcReaction<DisplacementCubitBcData>(
            *m_field_ptr, fe_post_proc_ptr, res)();
        MoFEMFunctionReturn(0);
      };
      fe_post_proc_ptr->postProcessHook = get_post_proc_hook_rhs;
      CHKERR DMoFEMPostProcessFiniteElements(dM, fe_post_proc_ptr.get());

      MoFEMFunctionReturn(0);
    };

    auto print_max_min = [&](auto &tuple, const std::string msg) {
      MoFEMFunctionBegin;
      CHKERR VecScatterBegin(std::get<1>(tuple), ts_u, std::get<0>(tuple),
                             INSERT_VALUES, SCATTER_FORWARD);
      CHKERR VecScatterEnd(std::get<1>(tuple), ts_u, std::get<0>(tuple),
                           INSERT_VALUES, SCATTER_FORWARD);
      double max, min;
      CHKERR VecMax(std::get<0>(tuple), PETSC_NULL, &max);
      CHKERR VecMin(std::get<0>(tuple), PETSC_NULL, &min);
      MOFEM_LOG_C("CONTACT", Sev::inform, "%s time %6.4e min %6.4e max %6.4e",
                  msg.c_str(), ts_t, min, max);
      MoFEMFunctionReturn(0);
    };

    auto print_force_and_area = [&]() {
      MoFEMFunctionBegin;
      MoFEM::Interface *m_field_ptr;
      CHKERR DMoFEMGetInterfacePtr(dM, &m_field_ptr);
      if (!m_field_ptr->get_comm_rank()) {
        const double *t_ptr;
        CHKERR VecGetArrayRead(CommonData::totalTraction, &t_ptr);
        MOFEM_LOG_C("CONTACT", Sev::inform,
                    "Contact force: time %6.4e Fx: %6.5e Fy: %6.5e Fz: %6.5e",
                    ts_t, t_ptr[0], t_ptr[1], t_ptr[2]);
        MOFEM_LOG_C("CONTACT", Sev::inform,
                    "Contact area: time %6.4e Active: %6.5e Potential: %6.5e",
                    ts_t, t_ptr[3], t_ptr[4]);
        CHKERR VecRestoreArrayRead(CommonData::totalTraction, &t_ptr);
      }
      MoFEMFunctionReturn(0);
    };

    if (mfrontInterface) {
      CHKERR mfrontInterface->updateElementVariables(
          dM, m_field_ptr->getInterface<Simple>()->getDomainFEName());
    }

    auto calculate_error = [&](MoFEM::ScalarFun &fun) {
      MoFEMFunctionBegin;
      struct OpCalcTractions : public BoundaryEleOp {
        OpCalcTractions(boost::shared_ptr<MatrixDouble> m_ptr,
                        boost::shared_ptr<VectorDouble> p_ptr,
                        boost::shared_ptr<VectorDouble> mag_ptr,
                        boost::shared_ptr<VectorDouble> traction_y_ptr,
                        boost::shared_ptr<MatrixDouble> t_ptr,
                        boost::shared_ptr<MatrixDouble> grad_sdf_ptr)
            : BoundaryEleOp(NOSPACE, OPSPACE), mPtr(m_ptr), pPtr(p_ptr),
              magPtr(mag_ptr), tyPtr(traction_y_ptr), tPtr(t_ptr),
              gradSDFPtr(grad_sdf_ptr) {}
        MoFEMErrorCode doWork(int, EntityType, EntitiesFieldData::EntData &) {
          MoFEMFunctionBegin;
          FTensor::Index<'i', SPACE_DIM> i;
          mPtr->resize(SPACE_DIM, pPtr->size());
          mPtr->clear();
          magPtr->resize(pPtr->size());
          magPtr->clear();
          tyPtr->resize(pPtr->size());
          tyPtr->clear();

          auto t_traction = getFTensor1FromMat<SPACE_DIM>(*mPtr);
          auto t_contact_traction = getFTensor1FromMat<SPACE_DIM>(*tPtr);
          auto t_p = getFTensor0FromVec(*pPtr);
          int nb_gauss_pts = pPtr->size();
          auto t_normal = getFTensor1FromMat<SPACE_DIM>(*gradSDFPtr);
          auto t_normal_at_gauss = getFTensor1NormalsAtGaussPts();
          auto t_mag = getFTensor0FromVec(*magPtr);
          auto t_ty = getFTensor0FromVec(*tyPtr);

          for (int gg = 0; gg != nb_gauss_pts; gg++) {
            FTensor::Tensor1<double, SPACE_DIM> normal;
            t_traction(i) = t_p * (-(t_normal(i) / t_normal.l2()));
            t_mag = t_contact_traction.l2();
            t_ty = t_contact_traction(1);

            ++t_normal;
            ++t_traction;
            ++t_p;
            ++t_mag;
            ++t_contact_traction;
            ++t_ty;
            ++t_normal_at_gauss;
          }
          MoFEMFunctionReturn(0);
        }

      private:
        boost::shared_ptr<MatrixDouble> mPtr;
        boost::shared_ptr<VectorDouble> pPtr;
        boost::shared_ptr<VectorDouble> magPtr;
        boost::shared_ptr<MatrixDouble> tPtr;
        boost::shared_ptr<VectorDouble> tyPtr;
        boost::shared_ptr<MatrixDouble> gradSDFPtr;
      };

      auto post_proc_norm_fe = boost::make_shared<BoundaryEle>(*m_field_ptr);
      auto common_data_ptr = boost::make_shared<ContactOps::CommonData>();
      auto simple = m_field_ptr->getInterface<Simple>();
      Range contact_range;
      for (auto m :
           m_field_ptr->getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
               std::regex((boost::format("%s(.*)") % "CONTACT").str()))) {
        auto meshset = m->getMeshset();
        Range contact_meshset_range;
        CHKERR m_field_ptr->get_moab().get_entities_by_dimension(
            meshset, SPACE_DIM - 1, contact_meshset_range, true);

        CHKERR m_field_ptr->getInterface<CommInterface>()->synchroniseEntities(
            contact_meshset_range);
        contact_range.merge(contact_meshset_range);
      }

      CHK_THROW_MESSAGE(
          (AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
              post_proc_norm_fe->getOpPtrVector(), {HDIV}, "GEOMETRY")),
          "Apply transform");
      // We have to integrate on curved face geometry, thus integration weight
      // have to adjusted.
      post_proc_norm_fe->getOpPtrVector().push_back(
          new OpSetHOWeightsOnSubDim<SPACE_DIM>());
      post_proc_norm_fe->getRuleHook = [](int, int, int approx_order) {
        return 2 * approx_order + geom_order - 1;
      };

      post_proc_norm_fe->getOpPtrVector().push_back(
          new OpCalculateVectorFieldValues<SPACE_DIM>(
              "U", common_data_ptr->contactDispPtr()));
      post_proc_norm_fe->getOpPtrVector().push_back(
          new OpCalculateHVecTensorTrace<SPACE_DIM, BoundaryEleOp>(
              "SIGMA", common_data_ptr->contactTractionPtr()));
      using C = ContactIntegrators<BoundaryEleOp>;
      post_proc_norm_fe->getOpPtrVector().push_back(
          new typename C::template OpEvaluateSDF<SPACE_DIM, GAUSS>(
              common_data_ptr));

      auto analytical_traction_ptr = boost::make_shared<MatrixDouble>();
      auto analytical_pressure_ptr = boost::make_shared<VectorDouble>();
      auto mag_traction_ptr = boost::make_shared<VectorDouble>();
      auto traction_y_ptr = boost::make_shared<VectorDouble>();
      auto contact_range_ptr = boost::make_shared<Range>(contact_range);

      post_proc_norm_fe->getOpPtrVector().push_back(
          new OpGetTensor0fromFunc(analytical_pressure_ptr, fun));

      post_proc_norm_fe->getOpPtrVector().push_back(new OpCalcTractions(
          analytical_traction_ptr, analytical_pressure_ptr, mag_traction_ptr,
          traction_y_ptr, common_data_ptr->contactTractionPtr(),
          common_data_ptr->gradSdfPtr()));

      post_proc_norm_fe->getOpPtrVector().push_back(
          new OpCalcNormL2Tensor1<SPACE_DIM>(
              common_data_ptr->contactTractionPtr(), normsVec, TRACTION_NORM_L2,
              analytical_traction_ptr, contact_range_ptr));

      // calculate magnitude of traction

      post_proc_norm_fe->getOpPtrVector().push_back(new OpCalcNormL2Tensor0(
          mag_traction_ptr, normsVec, MAG_TRACTION_NORM_L2,
          analytical_pressure_ptr, contact_range_ptr));

      post_proc_norm_fe->getOpPtrVector().push_back(
          new OpCalcNormL2Tensor0(traction_y_ptr, normsVec, TRACTION_Y_NORM_L2,
                                  analytical_pressure_ptr, contact_range_ptr));

      CHKERR VecZeroEntries(normsVec);
      post_proc_norm_fe->copyTs(*this); // set time as is in Monitor
      CHKERR DMoFEMLoopFiniteElements(dM, "bFE", post_proc_norm_fe);
      CHKERR VecAssemblyBegin(normsVec);
      CHKERR VecAssemblyEnd(normsVec);

      MOFEM_LOG_CHANNEL("SELF"); // Clear channel from old tags
      if (m_field_ptr->get_comm_rank() == 0) {
        const double *norms;
        CHKERR VecGetArrayRead(normsVec, &norms);
        MOFEM_TAG_AND_LOG("SELF", Sev::inform, "Errors")
            << "norm_traction: " << std::scientific
            << std::sqrt(norms[TRACTION_NORM_L2]);
        MOFEM_TAG_AND_LOG("SELF", Sev::inform, "Errors")
            << "norm_mag_traction: " << std::scientific
            << std::sqrt(norms[MAG_TRACTION_NORM_L2]);
        MOFEM_TAG_AND_LOG("SELF", Sev::inform, "Errors")
            << "norm_traction_y: " << std::scientific
            << std::sqrt(norms[TRACTION_Y_NORM_L2]);
        CHKERR VecRestoreArrayRead(normsVec, &norms);
      }
      MoFEMFunctionReturn(0);
    };

    int se = 1;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-save_every", &se, PETSC_NULL);

    if (!(ts_step % se)) {
      MOFEM_LOG("CONTACT", Sev::inform)
          << "Write file at time " << ts_t << " write step " << sTEP;
      CHKERR post_proc();
    }
    CHKERR calculate_force();
    CHKERR calculate_area();

    CHKERR calculate_reactions();
    if (atom_test == 1 && sTEP > 0)
      CHKERR calculate_error(analyticalHertzPressurePlaneStress);
    if (atom_test == 2 && sTEP > 0)
      CHKERR calculate_error(analyticalHertzPressurePlaneStrain);
    if (atom_test == 5 && sTEP > 0)
      CHKERR calculate_error(analyticalHertzPressureAxisymmetric);
    if (atom_test == 6 && sTEP > 0)
      CHKERR calculate_error(analyticalWavy2DPressure);

    CHKERR print_max_min(uXScatter, "Ux");
    CHKERR print_max_min(uYScatter, "Uy");
    if (SPACE_DIM == 3)
      CHKERR print_max_min(uZScatter, "Uz");
    CHKERR print_force_and_area();
    ++sTEP;

    MoFEMFunctionReturn(0);
  }

  //! [Analytical function]
  MoFEM::ScalarFun analyticalWavy2DPressure = [](double x, double y, double z) {
    // update to atom test values
    double E_star = young_modulus / (1 - poisson_ratio * poisson_ratio);
    // delta
    double delta = 0.0002;
    // lambda
    double lambda = 2;
    // pressure star
    double p_star = M_PI * E_star * delta / lambda;

    // Pressure = p_bar + p_star * cos(2 * pi * x / lambda)
    return p_star + p_star * std::cos(2. * M_PI * x / lambda);
  };

  MoFEM::ScalarFun analyticalHertzPressureAxisymmetric = [](double x, double y,
                                                            double z) {
    // update to atom test values
    double E_star = young_modulus / (1 - poisson_ratio * poisson_ratio);
    // Radius
    double R = 100.;
    // Indentation
    double d = 0.01;
    // Force
    double F = (4. / 3.) * E_star * std::sqrt(R) * std::pow(d, 1.5);
    // Contact area radius
    double a = std::pow((3. * F * R) / (4. * E_star), 1. / 3.);
    // Maximum pressure
    double p_max = (3. * F) / (2. * M_PI * a * a);

    double r = std::sqrt((x * x) + (y * y));

    if (r > a) {
      return 0.;
    }
    // Pressure = p_max * sqrt(1 - (r^2 / a^2))
    return p_max * std::sqrt(1 - ((r * r) / (a * a)));
  };

  MoFEM::ScalarFun analyticalHertzPressurePlaneStrain = [](double x, double y,
                                                           double z) {
    // update to atom test values
    double E_star = young_modulus / (1 - poisson_ratio * poisson_ratio);
    // Radius
    double R = 100.;
    // Indentation
    double d = 0.02745732273553991;
    // Contact area radius
    double a = 1;
    // current radius
    double r = std::sqrt((x * x) + (y * y));

    if (r > a) {
      return 0.;
    }
    // Pressure = p_max * sqrt(1 - (r^2 / a^2))
    return E_star / (2. * R) * std::sqrt(a * a - r * r);
  };
  MoFEM::ScalarFun analyticalHertzPressurePlaneStress = [](double x, double y,
                                                           double z) {
    // update to atom test values
    double E_star = young_modulus;
    // Radius
    double R = 100.;
    // Indentation
    double d = 0.02745732273553991;
    // Contact area radius
    double a = 1;
    // current radius
    double r = std::sqrt((x * x) + (y * y));

    if (r > a) {
      return 0.;
    }
    // Pressure = p_max * sqrt(1 - (r^2 / a^2))
    return E_star / (2. * R) * std::sqrt(a * a - r * r);
  };

  // ***DISPLACMENT NOT TESTED***
  MoFEM::VectorFun<SPACE_DIM> analyticalHertzDisplacement3D = [](double x,
                                                                 double y,
                                                                 double z) {
    // update to atom test values
    double E_star = young_modulus / (1 - poisson_ratio * poisson_ratio);
    // Radius
    double R = 100.;
    // Contact area radius
    double a = 1;
    // max pressure
    double p_0 = (2. * E_star * a) / (M_PI * R);
    // current radius
    double r = std::sqrt((x * x) + (y * y));
    FTensor::Tensor1<double, SPACE_DIM> u;
    std::vector<double> v_u;

    double u_z = 0.;
    double u_r = 0.;
    // outside contact zone
    if (r > a) {
      u_z = (1. - std::pow(poisson_ratio, 2.)) / young_modulus *
            ((p_0) / (2. * a)) *
            ((2. * std::pow(a, 2.) - std::pow(r, 2.)) * asin(a / r) +
             std::pow(r, 2.) * (a / r) *
                 std::pow(1 - (std::pow(a, 2.) / std::pow(r, 2.)), 2.));
      u_r = -((1. - 2. * poisson_ratio) * (1. + poisson_ratio)) /
            (3. * young_modulus) * ((std::pow(a, 2) / r)) * p_0;

      if (SPACE_DIM == 2)
        v_u = {u_r, u_z};
      else
        v_u = {u_r, u_z, u_r};

      for (int i = 0; i < SPACE_DIM; ++i)
        u(i) = v_u[i];

      return u;
    }

    // In contact zone
    u_z = ((1. - std::pow(poisson_ratio, 2.)) / young_modulus) *
          ((M_PI * p_0) / 4. * a) * (2. * std::pow(a, 2.) - std::pow(r, 2.));
    u_r = -((1. - 2. * poisson_ratio) * (1. + poisson_ratio)) /
          (3. * young_modulus) * ((std::pow(a, 2.) / r)) * p_0 *
          (1 - std::pow(1 - (std::pow(r, 2.) / std::pow(a, 2.)), 1.5));

    if (SPACE_DIM == 2)
      v_u = {u_r, u_z};
    else
      v_u = {u_r, u_z, u_r};

    for (int i = 0; i < SPACE_DIM; ++i)
      u(i) = v_u[i];

    return u;
  };

  MoFEMErrorCode setScatterVectors(
      std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> ux_scatter,
      std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uy_scatter,
      std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uz_scatter) {
    MoFEMFunctionBegin;
    uXScatter = ux_scatter;
    uYScatter = uy_scatter;
    uZScatter = uz_scatter;
    MoFEMFunctionReturn(0);
  }

  MoFEMErrorCode getErrorNorm(int normType) {
    const double *norm;
    CHKERR VecGetArrayRead(normsVec, &norm);
    double norm_val = std::sqrt(norm[normType]);
    CHKERR VecRestoreArrayRead(normsVec, &norm);
    return norm_val;
  }

private:
  SmartPetscObj<DM> dM;
  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uXScatter;
  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uYScatter;
  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uZScatter;

  boost::shared_ptr<moab::Core> postProcMesh = boost::make_shared<moab::Core>();

  boost::shared_ptr<PostProcEleDomain> postProcDomainFe;
  boost::shared_ptr<PostProcEleBdy> postProcBdyFe;

  boost::shared_ptr<BoundaryEle> integrateTraction;
  boost::shared_ptr<BoundaryEle> integrateArea;

  enum NORMS {
    TRACTION_NORM_L2 = 0,
    MAG_TRACTION_NORM_L2,
    TRACTION_Y_NORM_L2,
    LAST_NORM
  };
  SmartPetscObj<Vec> normsVec;

  moab::Core mbVertexPostproc;
  moab::Interface &moabVertex;

  double lastTime;
  double deltaTime;
  int sTEP;

  boost::shared_ptr<GenericElementInterface> mfrontInterface;
};

} // namespace ContactOps