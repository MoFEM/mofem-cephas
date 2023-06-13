/**
 * @file EssentialTemperatureCubitBcData.hpp
 * @brief Specialization for essential b.c. with TemperatureCubitBcData
 * @version 0.13.2
 * @date 2022-09-18
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _ESSENTIAL_TEMPERATURECUBITBCDATA_HPP_
#define _ESSENTIAL_TEMPERATURECUBITBCDATA_HPP_

namespace MoFEM {

// FIXME: OpEssentialRhsImpl<TemperatureCubitBcData, 1, 1, A, I, OpBase> not
// tested.

template <AssemblyType A, IntegrationType I, typename OpBase>
struct OpEssentialRhsImpl<TemperatureCubitBcData, 1, 1, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpSource<1, 1> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpSource<1, 1>;

  OpEssentialRhsImpl(const std::string field_name,
                     boost::shared_ptr<TemperatureCubitBcData> bc_data,
                     boost::shared_ptr<Range> ents_ptr,
                     std::vector<boost::shared_ptr<ScalingMethod>> smv);

private:
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
  double bcVal;
};

template <AssemblyType A, IntegrationType I, typename OpBase>
OpEssentialRhsImpl<TemperatureCubitBcData, 1, 1, A, I, OpBase>::
    OpEssentialRhsImpl(const std::string field_name,
                       boost::shared_ptr<TemperatureCubitBcData> bc_data,
                       boost::shared_ptr<Range> ents_ptr,
                       std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpSource(
          field_name, [this](double, double, double) { return bcVal; },
          ents_ptr),
      vecOfTimeScalingMethods(smv) {

  bcVal = bc_data->data.value1;
  this->timeScalingFun = [this](const double t) {
    double s = 1;
    for (auto &o : vecOfTimeScalingMethods) {
      s *= o->getScale(t);
    }
    return s;
  };

}

template <AssemblyType A, IntegrationType I, typename OpBase>
struct OpEssentialLhsImpl<TemperatureCubitBcData, 1, 1, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template BiLinearForm<
          I>::template OpMass<1, 1> {

  using OpMass = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template BiLinearForm<I>::template OpMass<1, 1>;

  OpEssentialLhsImpl(const std::string field_name,
                     boost::shared_ptr<Range> ents_ptr);
};

template <AssemblyType A, IntegrationType I, typename OpBase>
OpEssentialLhsImpl<TemperatureCubitBcData, 1, 1, A, I,
                   OpBase>::OpEssentialLhsImpl(const std::string field_name,
                                               boost::shared_ptr<Range>
                                                   ents_ptr)
    : OpMass(
          field_name, field_name,

          [](double, double, double) constexpr { return 1; },

          ents_ptr) {}

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

} // namespace MoFEM

#endif //_ESSENTIAL_TEMPERATURECUBITBCDATA_HPP_