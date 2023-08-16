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

/**
 * @brief A specialized version of DisplacementCubitBcData that includes an
 * additional rotation offset.
 *
 * DisplacementCubitBcDataWithRotation extends the DisplacementCubitBcData
 * structure with a rotation offset and a method to calculate the rotated
 * displacement given the rotation angles and coordinates.
 */
struct DisplacementCubitBcDataWithRotation : public DisplacementCubitBcData {

  std::array<double, 3> rotOffset;
  DisplacementCubitBcDataWithRotation() : rotOffset{0, 0, 0} {}

  /**
   * @brief Calculates the rotated displacement given the rotation angles,
   * coordinates, and an optional offset.
   *
   * @param angles A FTensor::Tensor1 containing the rotation angles.
   * @param coordinates A FTensor::Tensor1 containing the coordinates.
   * @param offset An optional FTensor::Tensor1 containing the offset
   * (default is {0., 0., 0.}).
   * @return FTensor::Tensor1<double, 3> representing the rotated displacement.
   */
  static FTensor::Tensor1<double, 3>
  GetRotDisp(const FTensor::Tensor1<double, 3> &angles,
             FTensor::Tensor1<double, 3> coordinates,
             FTensor::Tensor1<double, 3> offset = FTensor::Tensor1<double, 3>{
                 0., 0., 0.}) {

    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;
    FTensor::Index<'k', 3> k;

    FTensor::Tensor1<double, 3> rotated_displacement;
    rotated_displacement(i) = 0;

    auto get_rotation = [&](auto &omega) {
      FTensor::Tensor2<double, 3, 3> rotation_matrix;

      constexpr auto kronecker_delta = FTensor::Kronecker_Delta<int>();
      rotation_matrix(i, j) = kronecker_delta(i, j);

      const double angle = sqrt(omega(i) * omega(i));
      if (std::abs(angle) < std::numeric_limits<double>::epsilon())
        return rotation_matrix;

      FTensor::Tensor2<double, 3, 3> t_omega;
      t_omega(i, j) = FTensor::levi_civita<double>(i, j, k) * omega(k);
      const double a = sin(angle) / angle;
      const double sin_squared = sin(angle / 2.) * sin(angle / 2.);
      const double b = 2. * sin_squared / (angle * angle);
      rotation_matrix(i, j) += a * t_omega(i, j);
      rotation_matrix(i, j) += b * t_omega(i, k) * t_omega(k, j);

      return rotation_matrix;
    };

    auto rotation = get_rotation(angles);
    FTensor::Tensor1<double, 3> coordinate_distance;
    coordinate_distance(i) = offset(i) - coordinates(i);
    rotated_displacement(i) =
        coordinate_distance(i) - rotation(j, i) * coordinate_distance(j);

    return rotated_displacement;
  }
};

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
 * @brief Specialization for DisplacementCubitBcData
 *
 * @tparam
 */
template <> struct EssentialPostProcRhs<DisplacementCubitBcData> {
  EssentialPostProcRhs(MoFEM::Interface &m_field,
                      boost::shared_ptr<FEMethod> fe_ptr, double diag,
                      SmartPetscObj<Vec> rhs = nullptr);

  MoFEMErrorCode operator()();

protected:
  MoFEM::Interface &mField;
  boost::weak_ptr<FEMethod> fePtr;
  double vDiag;
  SmartPetscObj<Vec> vRhs;
};

/**
 * @brief Specialization for DisplacementCubitBcData
 *
 * @tparam
 */
template <> struct EssentialPostProcLhs<DisplacementCubitBcData> {
  EssentialPostProcLhs(MoFEM::Interface &m_field,
                      boost::shared_ptr<FEMethod> fe_ptr, double diag,
                      SmartPetscObj<Mat> lhs = nullptr,
                      SmartPetscObj<AO> ao = nullptr);

  MoFEMErrorCode operator()();

protected:
  MoFEM::Interface &mField;
  boost::weak_ptr<FEMethod> fePtr;
  double vDiag;
  SmartPetscObj<Mat> vLhs;
  SmartPetscObj<AO> vAO;
};

/**
 * @brief Specialization for DisplacementCubitBcData
 *
 * @tparam
 */
template <> struct EssentialPreProcReaction<DisplacementCubitBcData> {
  EssentialPreProcReaction(MoFEM::Interface &m_field,
                           boost::shared_ptr<FEMethod> fe_ptr,
                           SmartPetscObj<Vec> rhs = nullptr,
                           LogManager::SeverityLevel sev = Sev::inform);

  MoFEMErrorCode operator()();

protected:
  MoFEM::Interface &mField;
  boost::weak_ptr<FEMethod> fePtr;
  SmartPetscObj<Vec> vRhs;
  LogManager::SeverityLevel sevLevel;
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
  FTensor::Tensor1<double, 3> tOffset;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
  VectorFun<FIELD_DIM> dispFunction;

  FTensor::Tensor1<double, FIELD_DIM> rotFunction(double x, double y,
                                                  double z) {
    FTensor::Index<'i', FIELD_DIM> i;
    auto rot = DisplacementCubitBcDataWithRotation::GetRotDisp(
        tAngles, FTensor::Tensor1<double, 3>{x, y, z}, tOffset);
    FTensor::Tensor1<double, FIELD_DIM> t_ret;
    t_ret(i) = tVal(i) + rot(i);
    return t_ret;
  };
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
  tAngles(i) = 0;
  tOffset(i) = 0;

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

  if (auto ext_bc_data =
          dynamic_cast<DisplacementCubitBcDataWithRotation const *>(
              bc_data.get())) {
    for (int a = 0; a != 3; ++a)
      tOffset(a) = ext_bc_data->rotOffset[a];
  }

  this->timeScalingFun = [this](const double t) {
    double s = 1;
    for (auto &o : vecOfTimeScalingMethods) {
      s *= o->getScale(t);
    }
    return s;
  };
  if (bc_data->data.flag4 || bc_data->data.flag5 || bc_data->data.flag6) {
    this->dispFunction = [this](double x, double y, double z) {
      return this->rotFunction(x, y, z);
    };
  } else {
    // use default set at construction
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
