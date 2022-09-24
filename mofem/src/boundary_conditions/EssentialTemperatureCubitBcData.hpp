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