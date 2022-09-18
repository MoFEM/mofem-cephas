/**
 * @file EssentialHeatFluxCubitBcData.hpp
 * @brief Specialization for essential b.c. with HeatFluxCubitBcData
 * @version 0.13.2
 * @date 2022-09-18
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _ESSENTIAL_HEATFLUXCUBITBCDATA_HPP_
#define _ESSENTIAL_HEATFLUXCUBITBCDATA_HPP_

namespace MoFEM {

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

} // namespace MoFEM

#endif //_ESSENTIAL_HEATFLUXCUBITBCDATA_HPP_