/**
 * @file EssentialDisplacementCubitBcData.hpp
 * @brief Implementations specialisation for DisplacementCubitBcData
 * @version 0.13.2
 * @date 2022-09-18
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _ESSENTIAL_DISPLACEMENTCUBITBCDATA_HPP__
#define _ESSENTIAL_DISPLACEMENTCUBITBCDATA_HPP__

namespace MoFEM {

inline FTensor::Tensor1<double, 3>
_getRotDisp(FTensor::Tensor1<double, 3> t_angles, double x, double y,
            double z) {

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;
  FTensor::Tensor1<double, 3> tRot;
  tRot(i) = 0;

  auto get_rotation = [&](auto &t_omega) {
    FTensor::Tensor2<double, 3, 3> t_R;

    constexpr auto t_kd = FTensor::Kronecker_Delta<int>();
    t_R(i, j) = t_kd(i, j);

    const double angle = sqrt(t_omega(i) * t_omega(i));
    if (std::abs(angle) < std::numeric_limits<double>::epsilon())
      return t_R;

    FTensor::Tensor2<double, 3, 3> t_Omega;
    t_Omega(i, j) = FTensor::levi_civita<double>(i, j, k) * t_omega(k);
    const double a = sin(angle) / angle;
    const double ss_2 = sin(angle / 2.);
    const double b = 2. * ss_2 * ss_2 / (angle * angle);
    t_R(i, j) += a * t_Omega(i, j);
    t_R(i, j) += b * t_Omega(i, k) * t_Omega(k, j);

    return t_R;
  };

  FTensor::Tensor1<double, 3> t_coords(x, y, z);
  auto t_rot = get_rotation(t_angles);
  tRot(i) += t_rot(j, i) * t_coords(j) - t_coords(i);

  return tRot;
}

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
  FTensor::Tensor1<double, 3> tAngles;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
  boost::shared_ptr<DisplacementCubitBcData> bcData;
  VectorFun<FIELD_DIM> dispFunction;
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpEssentialRhsImpl<DisplacementCubitBcData, 1, FIELD_DIM, A, I, OpBase>::
    OpEssentialRhsImpl(const std::string field_name,
                       boost::shared_ptr<DisplacementCubitBcData> bc_data,
                       boost::shared_ptr<Range> ents_ptr,
                       std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpSource(field_name, dispFunction, ents_ptr),
      vecOfTimeScalingMethods(smv), bcData(bc_data) {
  static_assert(FIELD_DIM > 1, "Is not implemented for scalar field");

  FTensor::Index<'i', FIELD_DIM> i;
  tVal(i) = 0;
  tAngles(i) = 0;

  if (bc_data->data.flag1 == 1)
    tVal(0) = -bc_data->data.value1;
  if (bc_data->data.flag2 == 1 && FIELD_DIM > 1)
    tVal(1) = -bc_data->data.value2;
  if (bc_data->data.flag3 == 1 && FIELD_DIM > 2)
    tVal(2) = -bc_data->data.value3;
  if (bc_data->data.flag4 == 1 && FIELD_DIM > 2)
    tAngles(0) = -bc_data->data.value4;
  if (bc_data->data.flag5 == 1 && FIELD_DIM > 2)
    tAngles(1) = -bc_data->data.value5;
  if (bc_data->data.flag6 == 1 && FIELD_DIM > 1)
    tAngles(2) = -bc_data->data.value6;

  this->timeScalingFun = [this](const double t) {
    double s = 1;
    for (auto &o : vecOfTimeScalingMethods) {
      s *= o->getScale(t);
    }
    return s;
  };
  if (bc_data->data.flag4 || bc_data->data.flag5 || bc_data->data.flag6) {
    this->dispFunction = [this](double x, double y, double z) {
      FTensor::Index<'i', FIELD_DIM> i;
      auto rot = _getRotDisp(tAngles, x, y, z);
      tVal(i) += rot(i);
      return tVal;
    };
  } else {
    this->dispFunction = [this](double x, double y, double z) { return tVal; };
  }
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
struct AddEssentialToRhsPipelineImpl<

    OpEssentialRhsImpl<DisplacementCubitBcData, BASE_DIM, FIELD_DIM, A, I,
                       OpBase>,
    A, I, OpBase

    > {

  AddEssentialToRhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      MoFEM::Interface &m_field,
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      const std::string problem_name, std::string field_name,
      boost::shared_ptr<MatrixDouble> field_mat_ptr,
      std::vector<boost::shared_ptr<ScalingMethod>> smv

  ) {
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

  static MoFEMErrorCode add(

      MoFEM::Interface &m_field,
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      const std::string problem_name, std::string field_name

  ) {
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

} // namespace MoFEM

#endif // _ESSENTIAL_DISPLACEMENTCUBITBCDATA_HPP__
