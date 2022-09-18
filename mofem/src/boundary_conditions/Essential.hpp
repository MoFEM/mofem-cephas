/**
 * @file Essential.hpp
 * @brief Setting essential boundary conditions
 * @version 0.13.1
 * @date 2022-08-12
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _ESSENTIAL_BC_
#define _ESSENTIAL_BC_

namespace MoFEM {

/**
 * @brief Wrapper on user dtat (element) operator used to select specialization
 * of essential bc.
 *
 * @tparam T
 */
template <typename T> struct EssentialOpType {};

/**
 * @brief Specialisation for b.c. applied by different types of meshsets
 *
 * @tparam BC
 */
template <CubitBC BC> struct EssentialMeshsetType {};

/**
 * @brief Class (Function) to enforce essential constrains
 *
 * Class (Function) to enforce essential constrains for DOFs which were
 * removed from the system
 *
 * @tparam T
 */
template <typename T> struct EssentialPreProc {};

/**
 * @brief Specialization for DisplacementCubitBcData
 *
 * Specialization to enforce blocksets which DisplacementCubitBcData ptr. That
 * is to enforce displacement constraints. set
 *
 * @tparam
 */
template <> struct EssentialPreProc<DisplacementCubitBcData> {
  EssentialPreProc(MoFEM::Interface &m_field,
                   boost::shared_ptr<FEMethod> fe_ptr,
                   std::vector<boost::shared_ptr<ScalingMethod>> smv,
                   bool get_coords = false);

  MoFEMErrorCode operator()();

protected:
  MoFEM::Interface &mField;
  boost::weak_ptr<FEMethod> fePtr;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
  bool getCoords;
};

/**
 * @brief Specialization for TemperatureCubitBcData
 *
 * Specialization to enforce blocksets which TemperatureCubitBcData ptr. That is
 * to enforce constrains on temperature. set
 *
 * @tparam
 */
template <> struct EssentialPreProc<TemperatureCubitBcData> {
  EssentialPreProc(MoFEM::Interface &m_field,
                   boost::shared_ptr<FEMethod> fe_ptr,
                   std::vector<boost::shared_ptr<ScalingMethod>> smv);

  MoFEMErrorCode operator()();

protected:
  MoFEM::Interface &mField;
  boost::weak_ptr<FEMethod> fePtr;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
};

/**
 * @brief Enforce essential constrains on rhs.
 * 
 * This class is used when constrains are enforced by least square method.
 * 
 * @tparam T 
 * @tparam BASE_DIM 
 * @tparam FIELD_DIM 
 * @tparam A 
 * @tparam I 
 * @tparam OpBase 
 */
template <typename T, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct OpEssentialRhsImpl;

/**
 * @brief Enforce essential constrains on lhs
 * 
 * This class is used when constrains are enforced by least square method.
 * 
 * @tparam T 
 * @tparam BASE_DIM 
 * @tparam FIELD_DIM 
 * @tparam A 
 * @tparam I 
 * @tparam OpBase 
 */
template <typename T, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct OpEssentialLhsImpl;

/**
 * @brief Function (factory) for setting operators for rhs pipeline
 * 
 * @tparam T 
 * @tparam A 
 * @tparam I 
 * @tparam OpBase 
 */
template <typename T, AssemblyType A, IntegrationType I, typename OpBase>
struct AddEssentialToRhsPipelineImpl;

/**
 * @brief Function (factory) for setting operators for lhs pipeline
 * 
 * @tparam T 
 * @tparam A 
 * @tparam I 
 * @tparam OpBase 
 */
template <typename T, AssemblyType A, IntegrationType I, typename OpBase>
struct AddEssentialToLhsPipelineImpl;

/**
 * @brief Natural boundary conditions
 * @ingroup mofem_forms
 *
 * @tparam EleOp
 */
template <typename EleOp> struct EssentialBC {

  using EntData = EntitiesFieldData::EntData;
  using OpType = typename EleOp::OpType;

  /**
   * @brief Assembly methods
   * @ingroup mofem_forms
   *
   * @tparam A
   */
  template <AssemblyType A> struct Assembly {

    template <IntegrationType I> struct LinearForm {
      template <typename T, int BASE_DIM, int FIELD_DIM>
      using OpEssentialRhs =
          OpEssentialRhsImpl<T, BASE_DIM, FIELD_DIM, A, I, EleOp>;

      template <typename T>
      static MoFEMErrorCode addEssentialToRhsPipeline(
          EssentialOpType<T>, MoFEM::Interface &m_field,
          boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
          std::string problem_name, std::string field_name,
          boost::shared_ptr<MatrixDouble> field_mat_ptr,
          std::vector<boost::shared_ptr<ScalingMethod>> smv);
    };

    template <IntegrationType I> struct BiLinearForm {
      template <typename T, int BASE_DIM, int FIELD_DIM>
      using OpEssentialLhs =
          OpEssentialLhsImpl<T, BASE_DIM, FIELD_DIM, A, I, EleOp>;
      template <typename T>
      static MoFEMErrorCode addEssentialToLhsPipeline(
          EssentialOpType<T>, MoFEM::Interface &m_field,
          boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
          std::string problem_name, std::string field_name);
    };
  };
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpEssentialRhsImpl<DisplacementCubitBcData, 1, FIELD_DIM, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpSource<1, FIELD_DIM> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpSource<1, FIELD_DIM>;

  OpEssentialRhsImpl(const std::string field_name,
                     boost::shared_ptr<DisplacementCubitBcData> bc_data,
                     boost::shared_ptr<Range> ents_ptr,
                     std::vector<boost::shared_ptr<ScalingMethod>> smv);

private:
  FTensor::Tensor1<double, FIELD_DIM> tVal;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpEssentialRhsImpl<DisplacementCubitBcData, 1, FIELD_DIM, A, I, OpBase>::
    OpEssentialRhsImpl(const std::string field_name,
                       boost::shared_ptr<DisplacementCubitBcData> bc_data,
                       boost::shared_ptr<Range> ents_ptr,
                       std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpSource(
          field_name, [this](double, double, double) { return tVal; },
          ents_ptr),
      vecOfTimeScalingMethods(smv) {
  static_assert(FIELD_DIM > 1, "Is not implemented for scalar field");

  FTensor::Index<'i', FIELD_DIM> i;
  tVal(i) = 0;
  if (bc_data->data.flag1 == 1)
    tVal(0) = -bc_data->data.value1;
  if (bc_data->data.flag2 == 1 && FIELD_DIM > 1)
    tVal(1) = -bc_data->data.value2;
  if (bc_data->data.flag3 == 1 && FIELD_DIM > 2)
    tVal(2) = -bc_data->data.value3;

  this->timeScalingFun = [this](const double t) {
    double s = 1;
    for (auto &o : vecOfTimeScalingMethods) {
      s *= o->getScale(t);
    }
    return s;
  };
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpEssentialRhsImpl<HeatFluxCubitBcData, 3, FIELD_DIM, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<I>::
          template OpSource<3, FIELD_DIM, SourceBoundaryNormalSpecialization> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<A>::
      template LinearForm<I>::template OpSource<
          3, FIELD_DIM, SourceBoundaryNormalSpecialization>;

  OpEssentialRhsImpl(const std::string field_name,
                     boost::shared_ptr<HeatFluxCubitBcData> bc_data,
                     boost::shared_ptr<Range> ents_ptr,
                     std::vector<boost::shared_ptr<ScalingMethod>> smv);

private:
  double heatFlux;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpEssentialRhsImpl<HeatFluxCubitBcData, 3, FIELD_DIM, A, I, OpBase>::
    OpEssentialRhsImpl(const std::string field_name,
                       boost::shared_ptr<HeatFluxCubitBcData> bc_data,
                       boost::shared_ptr<Range> ents_ptr,
                       std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpSource(
          field_name, [this](double, double, double) { return heatFlux; },
          ents_ptr),
      vecOfTimeScalingMethods(smv) {
  heatFlux = -bc_data->data.value1;
  this->timeScalingFun = [this](const double t) {
    double s = 1;
    for (auto &o : vecOfTimeScalingMethods) {
      s *= o->getScale(t);
    }
    return s;
  };
}

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct OpEssentialLhsImpl<DisplacementCubitBcData, BASE_DIM, FIELD_DIM, A, I,
                          OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template BiLinearForm<
          I>::template OpMass<BASE_DIM, FIELD_DIM> {

  using OpMass = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template BiLinearForm<I>::template OpMass<BASE_DIM, FIELD_DIM>;

  OpEssentialLhsImpl(const std::string field_name,
                     boost::shared_ptr<Range> ents_ptr);
};

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
OpEssentialLhsImpl<DisplacementCubitBcData, BASE_DIM, FIELD_DIM, A, I,
                   OpBase>::OpEssentialLhsImpl(const std::string field_name,
                                               boost::shared_ptr<Range>
                                                   ents_ptr)
    : OpMass(
          field_name, field_name,

          [](double, double, double) constexpr { return 1; },

          ents_ptr) {}

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct OpEssentialLhsImpl<HeatFluxCubitBcData, BASE_DIM, FIELD_DIM, A, I,
                          OpBase>
    : OpEssentialLhsImpl<DisplacementCubitBcData, BASE_DIM, FIELD_DIM, A, I,
                         OpBase> {
  using OpEssentialLhsImpl<DisplacementCubitBcData, BASE_DIM, FIELD_DIM, A, I,
                           OpBase>::OpEssentialLhsImpl;
};

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct AddEssentialToRhsPipelineImpl<

    OpEssentialRhsImpl<DisplacementCubitBcData, BASE_DIM, FIELD_DIM, A, I,
                       OpBase>,
    A, I, OpBase

    > {

  AddEssentialToRhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      EssentialOpType<OpEssentialRhsImpl<DisplacementCubitBcData, BASE_DIM,
                                         FIELD_DIM, A, I, OpBase>>,

      MoFEM::Interface &m_field,
      boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      const std::string problem_name, std::string field_name,
      boost::shared_ptr<MatrixDouble> field_mat_ptr,
      std::vector<boost::shared_ptr<ScalingMethod>> smv) {
    MoFEMFunctionBegin;

    using OP =
        typename EssentialBC<OpBase>::template Assembly<A>::template LinearForm<
            I>::template OpEssentialRhs<DisplacementCubitBcData, BASE_DIM,
                                        FIELD_DIM>;
    using OpInternal = typename FormsIntegrators<OpBase>::template Assembly<
        A>::template LinearForm<I>::template OpBaseTimesVector<BASE_DIM,
                                                               FIELD_DIM, 1>;

    auto add_op = [&](auto &bcs) {
      MoFEMFunctionBeginHot;
      for (auto &m : bcs) {
        if (auto bc = m.second->dispBcPtr) {
          auto &bc_id = m.first;
          auto regex_str =
              (boost::format("%s_%s_(.*)") % problem_name % field_name).str();
          if (std::regex_match(bc_id, std::regex(regex_str))) {
            MOFEM_TAG_AND_LOG("SELF", Sev::noisy, "OpEssentialRhs") << *bc;
            pipeline.push_back(
                new OpSetBc(field_name, false, m.second->getBcMarkersPtr()));
            pipeline.push_back(
                new OP(field_name, bc, m.second->getBcEntsPtr(), smv));
            pipeline.push_back(new OpInternal(
                field_name, field_mat_ptr,
                [](double, double, double) constexpr { return 1.; },
                m.second->getBcEntsPtr()));
            pipeline.push_back(new OpUnSetBc(field_name));
          }
        }
      }
      MOFEM_LOG_CHANNEL("SELF");
      MoFEMFunctionReturnHot(0);
    };

    CHKERR add_op(m_field.getInterface<BcManager>()->getBcMapByBlockName());

    MoFEMFunctionReturn(0);
  }
};

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct AddEssentialToLhsPipelineImpl<

    OpEssentialLhsImpl<DisplacementCubitBcData, BASE_DIM, FIELD_DIM, A, I,
                       OpBase>,
    A, I, OpBase

    > {

  AddEssentialToLhsPipelineImpl() = delete;

  static MoFEMErrorCode
  add(EssentialOpType<OpEssentialLhsImpl<DisplacementCubitBcData, BASE_DIM,
                                         FIELD_DIM, A, I, OpBase>>,
      MoFEM::Interface &m_field,
      boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      const std::string problem_name, std::string field_name) {
    MoFEMFunctionBegin;

    using OP = typename EssentialBC<OpBase>::template Assembly<A>::
        template BiLinearForm<I>::template OpEssentialLhs<
            DisplacementCubitBcData, BASE_DIM, FIELD_DIM>;

    auto add_op = [&](auto &bcs) {
      MoFEMFunctionBeginHot;
      for (auto &m : bcs) {
        if (auto bc = m.second->dispBcPtr) {
          auto &bc_id = m.first;
          auto regex_str =
              (boost::format("%s_%s_(.*)") % problem_name % field_name).str();
          if (std::regex_match(bc_id, std::regex(regex_str))) {
            MOFEM_TAG_AND_LOG("SELF", Sev::noisy, "OpEssentialLhs") << *bc;
            pipeline.push_back(
                new OpSetBc(field_name, false, m.second->getBcMarkersPtr()));
            pipeline.push_back(new OP(field_name, m.second->getBcEntsPtr()));
            pipeline.push_back(new OpUnSetBc(field_name));
          }
        }
      }
      MOFEM_LOG_CHANNEL("SELF");
      MoFEMFunctionReturnHot(0);
    };

    CHKERR add_op(m_field.getInterface<BcManager>()->getBcMapByBlockName());

    MoFEMFunctionReturn(0);
  }
};

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct AddEssentialToRhsPipelineImpl<

    OpEssentialRhsImpl<HeatFluxCubitBcData, BASE_DIM, FIELD_DIM, A, I, OpBase>,
    A, I, OpBase

    > {

  AddEssentialToRhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      EssentialOpType<OpEssentialRhsImpl<HeatFluxCubitBcData, BASE_DIM,
                                         FIELD_DIM, A, I, OpBase>>,

      MoFEM::Interface &m_field,
      boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      const std::string problem_name, std::string field_name,
      boost::shared_ptr<MatrixDouble> field_mat_ptr,
      std::vector<boost::shared_ptr<ScalingMethod>> smv) {
    MoFEMFunctionBegin;

    using OP =
        typename EssentialBC<OpBase>::template Assembly<A>::template LinearForm<
            I>::template OpEssentialRhs<HeatFluxCubitBcData, BASE_DIM,
                                        FIELD_DIM>;
    using OpInternal = typename FormsIntegrators<OpBase>::template Assembly<
        A>::template LinearForm<I>::template OpBaseTimesVector<BASE_DIM,
                                                               FIELD_DIM, 1>;

    auto add_op = [&](auto &bcs) {
      MoFEMFunctionBeginHot;
      for (auto &m : bcs) {
        if (auto bc = m.second->heatFluxBcPtr) {
          auto &bc_id = m.first;
          auto regex_str =
              (boost::format("%s_%s_(.*)") % problem_name % field_name).str();
          if (std::regex_match(bc_id, std::regex(regex_str))) {
            MOFEM_TAG_AND_LOG("SELF", Sev::noisy, "OpEssentialRhs") << *bc;
            pipeline.push_back(
                new OpSetBc(field_name, false, m.second->getBcMarkersPtr()));
            pipeline.push_back(
                new OP(field_name, bc, m.second->getBcEntsPtr(), smv));
            pipeline.push_back(new OpInternal(
                field_name, field_mat_ptr,
                [](double, double, double) constexpr { return 1.; },
                m.second->getBcEntsPtr()));
            pipeline.push_back(new OpUnSetBc(field_name));
          }
        }
      }
      MOFEM_LOG_CHANNEL("SELF");
      MoFEMFunctionReturnHot(0);
    };

    CHKERR add_op(m_field.getInterface<BcManager>()->getBcMapByBlockName());

    MoFEMFunctionReturn(0);
  }
};

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct AddEssentialToLhsPipelineImpl<

    OpEssentialLhsImpl<HeatFluxCubitBcData, BASE_DIM, FIELD_DIM, A, I, OpBase>,
    A, I, OpBase

    > {

  AddEssentialToLhsPipelineImpl() = delete;

  static MoFEMErrorCode
  add(EssentialOpType<OpEssentialLhsImpl<HeatFluxCubitBcData, BASE_DIM,
                                         FIELD_DIM, A, I, OpBase>>,
      MoFEM::Interface &m_field,
      boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      const std::string problem_name, std::string field_name) {
    MoFEMFunctionBegin;

    using OP = typename EssentialBC<OpBase>::template Assembly<A>::
        template BiLinearForm<I>::template OpEssentialLhs<HeatFluxCubitBcData,
                                                          BASE_DIM, FIELD_DIM>;

    auto add_op = [&](auto &bcs) {
      MoFEMFunctionBeginHot;
      for (auto &m : bcs) {
        if (auto bc = m.second->heatFluxBcPtr) {
          auto &bc_id = m.first;
          auto regex_str =
              (boost::format("%s_%s_(.*)") % problem_name % field_name).str();
          if (std::regex_match(bc_id, std::regex(regex_str))) {
            MOFEM_TAG_AND_LOG("SELF", Sev::noisy, "OpEssentialLhs") << *bc;
            pipeline.push_back(
                new OpSetBc(field_name, false, m.second->getBcMarkersPtr()));
            pipeline.push_back(new OP(field_name, m.second->getBcEntsPtr()));
            pipeline.push_back(new OpUnSetBc(field_name));
          }
        }
      }
      MOFEM_LOG_CHANNEL("SELF");
      MoFEMFunctionReturnHot(0);
    };

    CHKERR add_op(m_field.getInterface<BcManager>()->getBcMapByBlockName());

    MoFEMFunctionReturn(0);
  }
};

template <typename OpBase>
template <AssemblyType A>
template <IntegrationType I>
template <typename T>
MoFEMErrorCode
EssentialBC<OpBase>::Assembly<A>::LinearForm<I>::addEssentialToRhsPipeline(
    EssentialOpType<T>, MoFEM::Interface &m_field,
    boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::string problem_name, std::string field_name,
    boost::shared_ptr<MatrixDouble> field_mat_ptr,
    std::vector<boost::shared_ptr<ScalingMethod>> smv) {
  return AddEssentialToRhsPipelineImpl<T, A, I, OpBase>::add(
      EssentialOpType<T>(), m_field, pipeline, problem_name, field_name,
      field_mat_ptr, smv);
}

template <typename OpBase>
template <AssemblyType A>
template <IntegrationType I>
template <typename T>
MoFEMErrorCode
EssentialBC<OpBase>::Assembly<A>::BiLinearForm<I>::addEssentialToLhsPipeline(
    EssentialOpType<T>, MoFEM::Interface &m_field,
    boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::string problem_name, std::string field_name) {
  return AddEssentialToLhsPipelineImpl<T, A, I, OpBase>::add(
      EssentialOpType<T>(), m_field, pipeline, problem_name, field_name);
}

} // namespace MoFEM

#endif //_ESSENTIAL_BC_