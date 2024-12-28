/** \file MFrontInterface.cpp
 * \brief MFrontInterface
 *
 * MFrontInterface
 *
 */

#ifdef WITH_MGIS

#include <MFrontInterface.hpp>

namespace MoFEM {

template struct MFrontInterface<TRIDIMENSIONAL>;
template struct MFrontInterface<AXISYMMETRICAL>;
template struct MFrontInterface<PLANESTRAIN>;

template <ModelHypothesis H>
MFrontInterface<H>::MFrontInterface(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)) {

  // if (!LogManager::checkIfChannelExist("FieldEvaluatorWorld")) {
  //   auto core_log = logging::core::get();

  //   core_log->add_sink(LogManager::createSink(LogManager::getStrmWorld(),
  //                                             "FieldEvaluatorWorld"));
  //   core_log->add_sink(LogManager::createSink(LogManager::getStrmSync(),
  //                                             "FieldEvaluatorSync"));
  //   core_log->add_sink(LogManager::createSink(LogManager::getStrmSelf(),
  //                                             "FieldEvaluatorSelf"));

  //   LogManager::setLog("FieldEvaluatorWorld");
  //   LogManager::setLog("FieldEvaluatorSync");
  //   LogManager::setLog("FieldEvaluatorSelf");

  //   MOFEM_LOG_TAG("FieldEvaluatorWorld", "FieldEvaluator");
  //   MOFEM_LOG_TAG("FieldEvaluatorSync", "FieldEvaluator");
  //   MOFEM_LOG_TAG("FieldEvaluatorSelf", "FieldEvaluator");
  // }

  // MOFEM_LOG("FieldEvaluatorWorld", Sev::noisy) << "Field evaluator
  // intreface";
}

template <ModelHypothesis H>
MoFEMErrorCode
MFrontInterface<H>::query_interface(boost::typeindex::type_index type_index,
                                    UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = const_cast<MFrontInterface<H> *>(this);
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM

#endif // WITH_MGIS
