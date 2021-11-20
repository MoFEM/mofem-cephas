/** \file PipelineManager.hpp
 * \brief Header file for basic interface
 * \ingroup mofem_basic_interface
 *
 * Make basic interface, to speedup problem setup and analysts.
 * See discussion here
 * <a
 * href=https://groups.google.com/d/msg/mofem-group/Vkc00aia4dU/o9RF3ZmPAAAJ>link
 * to google groups</a>
 *
 */

/* MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

#ifndef __BASIC_HPP__
#define __BASIC_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

/**
 * \brief PipelineManager interface
 * \ingroup mofem_basic_interface
 */
struct PipelineManager : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  PipelineManager(const MoFEM::Core &core);

  using UserDataOperator = MoFEM::ForcesAndSourcesCore::UserDataOperator;
  using RuleHookFun = MoFEM::ForcesAndSourcesCore::RuleHookFun;

  using VolEle = MoFEM::VolumeElementForcesAndSourcesCore;
  using FaceEle = MoFEM::FaceElementForcesAndSourcesCore;
  using EdgeEle = MoFEM::EdgeElementForcesAndSourcesCore;

  inline boost::shared_ptr<FEMethod> &getDomainLhsFE();

  inline boost::shared_ptr<FEMethod> &getDomainRhsFE();

  inline boost::shared_ptr<FEMethod> &getBoundaryLhsFE();

  inline boost::shared_ptr<FEMethod> &getBoundaryRhsFE();

  inline boost::shared_ptr<FEMethod> &getSkeletonLhsFE();

  inline boost::shared_ptr<FEMethod> &getSkeletonRhsFE();

  inline boost::shared_ptr<FEMethod> &getDomainExplicitRhsFE();

  inline boost::shared_ptr<FEMethod> &getBoundaryExplicitRhsFE();

  inline boost::shared_ptr<FEMethod> &getSkeletonExplicitRhsFE();

  template <int DIM = -1>
  inline MoFEMErrorCode setDomainLhsIntegrationRule(RuleHookFun rule);

  template <int DIM = -1>
  inline MoFEMErrorCode setDomainRhsIntegrationRule(RuleHookFun rule);

  template <int DIM = -1>
  inline MoFEMErrorCode setBoundaryLhsIntegrationRule(RuleHookFun rule);

  template <int DIM = -1>
  inline MoFEMErrorCode setBoundaryRhsIntegrationRule(RuleHookFun rule);

  template <int DIM = -1>
  inline MoFEMErrorCode setSkeletonLhsIntegrationRule(RuleHookFun rule);

  template <int DIM = -1>
  inline MoFEMErrorCode setSkeletonRhsIntegrationRule(RuleHookFun rule);

  template <int DIM = -1>
  inline MoFEMErrorCode setDomainExplicitRhsIntegrationRule(RuleHookFun rule);

  template <int DIM = -1>
  inline MoFEMErrorCode setBoundaryExplicitRhsIntegrationRule(RuleHookFun rule);

  template <int DIM = -1>
  inline MoFEMErrorCode setSkeletonExplicitRhsIntegrationRule(RuleHookFun rule);

  /**
   * @brief Get the Op Domain Lhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @tparam -1
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &
  getOpDomainLhsPipeline();

  /**
   * @brief Get the Op Domain Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @tparam -1
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &
  getOpDomainRhsPipeline();

  /**
   * @brief Get the Op Boundary Lhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @tparam -1
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &
  getOpBoundaryLhsPipeline();

  /**
   * @brief Get the Op Boundary Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @tparam -1
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &
  getOpBoundaryRhsPipeline();

  /**
   * @brief Get the Op Skeleton Lhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &
  getOpSkeletonLhsPipeline();

  /**
   * @brief Get the Op Skeleton Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @tparam -1
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &
  getOpSkeletonRhsPipeline();

  /**
   * @brief Get the Op Domain Rhs Pipeline object for implicit-explicit G term
   * @ingroup mofem_basic_interface
   *
   * @tparam -1
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &
  getOpDomainExplicitRhsPipeline();

  /**
   * @brief Get the Op Bondary Rhs Pipeline object for implicit-explicit G term
   * @ingroup mofem_basic_interface
   *
   * @tparam -1
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &getOpBoundaryExplicitRhsPipeline();

  /**
   * @brief Get the Op Skeleton Rhs Pipeline object for implicit-explicit G term
   * @ingroup mofem_basic_interface
   *
   * @tparam -1
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &getOpSkeletonExplicitRhsPipeline();

  /**
   * @brief Iterate finite elements
   * @ingroup mofem_basic_interface
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode loopFiniteElements(SmartPetscObj<DM> dm = nullptr);

  /**
   * @brief Create KSP (linear) solver
   * @ingroup mofem_basic_interface
   *
   * @param dm
   * @return SmartPetscObj<KSP>
   */
  SmartPetscObj<KSP> createKSP(SmartPetscObj<DM> dm = nullptr);

  /**
   * @brief Create SNES (nonlinear) solver
   * @ingroup mofem_basic_interface
   *
   * @param dm
   * @return SmartPetscObj<SNES>
   */
  SmartPetscObj<SNES> createSNES(SmartPetscObj<DM> dm = nullptr);

  enum TSType { EX, IM, IM2, IMEX };

  /**
   * @brief reate TS (time) solver
   * 
   * @param type Type of time solver PipelineManager:EX/IM/IM2/IMEX
   * @param dm 
   * @return SmartPetscObj<TS> 
   */
  SmartPetscObj<TS> createTS(const TSType type, SmartPetscObj<DM> dm = nullptr);

  /**
   * @brief Create TS (time) explit solver
   * @ingroup mofem_basic_interface
   *
   * @param dm
   * @return SmartPetscObj<TS>
   */
  SmartPetscObj<TS> createTSEX(SmartPetscObj<DM> dm = nullptr);


  /**
   * @brief Create TS (time) implicit solver
   * @ingroup mofem_basic_interface
   *
   * @param dm
   * @return SmartPetscObj<TS>
   */
  SmartPetscObj<TS> createTSIM(SmartPetscObj<DM> dm = nullptr);

  /**
   * @deprecated  Use version with explicit TS solver type
   */
  DEPRECATED auto createTS(SmartPetscObj<DM> dm = nullptr) {
    return createTSIM(dm);
  }

  /**
   * @brief Create TS (time) solver for second order equation in time
   * @ingroup mofem_basic_interface
   *
   * @param dm
   * @return SmartPetscObj<TS>
   */
  SmartPetscObj<TS> createTSIM2(SmartPetscObj<DM> dm = nullptr);

  /**
   * @deprecated  Change name. Use createTSIM2 instead.
   */
  inline DEPRECATED auto createTS2(SmartPetscObj<DM> dm = nullptr) {
    return createTSIM2(dm);
  }

  /**
   * @brief Create TS (time) implicit-explicit solver
   * @ingroup mofem_basic_interface
   *
   * @param dm
   * @return SmartPetscObj<TS>
   */
  SmartPetscObj<TS> createTSIMEX(SmartPetscObj<DM> dm = nullptr);

private:
  MoFEM::Core &cOre;

  boost::shared_ptr<FEMethod>
      feDomainRhs; ///< Element to assemble RHS side by integrating domain
  boost::shared_ptr<FEMethod>
      feDomainLhs; ///< Element to assemble LHS side by integrating domain
  boost::shared_ptr<FEMethod>
      feBoundaryRhs; ///< Element to assemble RHS side by integrating boundary
  boost::shared_ptr<FEMethod>
      feBoundaryLhs; ///< Element to assemble LHS side by integrating boundary
  boost::shared_ptr<FEMethod>
      feSkeletonRhs; ///< Element to assemble RHS side by integrating skeleton
  boost::shared_ptr<FEMethod>
      feSkeletonLhs; ///< Element to assemble LHS side by integrating skeleton

  boost::shared_ptr<FEMethod>
      feDomainExplicitRhs; ///< Element to assemble explict Rhs for IMEX solver
  boost::shared_ptr<FEMethod>
      feBoundaryExplicitRhs; ///< Element to assemble explict Rhs for IMEX solver
  boost::shared_ptr<FEMethod>
      feSkeletonExplicitRhs; ///< Element to assemble explict Rhs for IMEX solver

  template <int DIM>
  inline boost::shared_ptr<FEMethod> &
  createDomainFEPipeline(boost::shared_ptr<FEMethod> &fe);

  template <int DIM>
  inline boost::shared_ptr<FEMethod> &
  createBoundaryFEPipeline(boost::shared_ptr<FEMethod> &fe);
};

template <int DIM>
boost::shared_ptr<FEMethod> &
PipelineManager::createDomainFEPipeline(boost::shared_ptr<FEMethod> &fe) {
  static_assert(DIM == 1 || DIM == 2 || DIM == 3, "not implemented");
  fe = boost::make_shared<FEMethod>();
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
PipelineManager::createDomainFEPipeline<3>(boost::shared_ptr<FEMethod> &fe) {
  if (!fe)
    fe = boost::make_shared<VolEle>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
PipelineManager::createDomainFEPipeline<2>(boost::shared_ptr<FEMethod> &fe) {
  if (!fe)
    fe = boost::make_shared<FaceEle>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
PipelineManager::createDomainFEPipeline<1>(boost::shared_ptr<FEMethod> &fe) {
  if (!fe)
    fe = boost::make_shared<EdgeEle>(cOre);
  return fe;
}

template <int DIM>
boost::shared_ptr<FEMethod> &
PipelineManager::createBoundaryFEPipeline(boost::shared_ptr<FEMethod> &fe) {
  static_assert(DIM == 1 || DIM == 2 || DIM == 3, "not implemented");
  fe = boost::make_shared<FEMethod>();
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
PipelineManager::createBoundaryFEPipeline<3>(boost::shared_ptr<FEMethod> &fe) {
  if (!fe)
    fe = boost::make_shared<FaceElementForcesAndSourcesCore>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
PipelineManager::createBoundaryFEPipeline<2>(boost::shared_ptr<FEMethod> &fe) {
  if (!fe)
    fe = boost::make_shared<EdgeEle>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
PipelineManager::createBoundaryFEPipeline<1>(boost::shared_ptr<FEMethod> &fe) {
  if (!fe)
    fe = boost::make_shared<VertexElementForcesAndSourcesCore>(cOre);
  return fe;
}

boost::shared_ptr<FEMethod> &PipelineManager::getDomainLhsFE() { return feDomainLhs; }

boost::shared_ptr<FEMethod> &PipelineManager::getDomainRhsFE() { return feDomainRhs; }

boost::shared_ptr<FEMethod> &PipelineManager::getBoundaryLhsFE() { return feBoundaryLhs; }

boost::shared_ptr<FEMethod> &PipelineManager::getBoundaryRhsFE() { return feBoundaryRhs; }

boost::shared_ptr<FEMethod> &PipelineManager::getSkeletonLhsFE() { return feSkeletonLhs; }

boost::shared_ptr<FEMethod> &PipelineManager::getSkeletonRhsFE() { return feSkeletonRhs; }

boost::shared_ptr<FEMethod> &PipelineManager::getDomainExplicitRhsFE() {
  return feDomainExplicitRhs;
}

boost::shared_ptr<FEMethod> &PipelineManager::getBoundaryExplicitRhsFE() {
  return feBoundaryExplicitRhs;
}

boost::shared_ptr<FEMethod> &PipelineManager::getSkeletonExplicitRhsFE() {
  return feSkeletonExplicitRhs;
}

template <int DIM>
MoFEMErrorCode PipelineManager::setDomainLhsIntegrationRule(PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feDomainLhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
PipelineManager::setDomainLhsIntegrationRule<-1>(PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return setDomainLhsIntegrationRule<1>(rule);
  case 2:
    return setDomainLhsIntegrationRule<2>(rule);
  case 3:
    return setDomainLhsIntegrationRule<3>(rule);
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

template <int DIM>
MoFEMErrorCode PipelineManager::setDomainRhsIntegrationRule(PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feDomainRhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
PipelineManager::setDomainRhsIntegrationRule<-1>(PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return setDomainRhsIntegrationRule<1>(rule);
  case 2:
    return setDomainRhsIntegrationRule<2>(rule);
  case 3:
    return setDomainRhsIntegrationRule<3>(rule);
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

template <int DIM>
MoFEMErrorCode PipelineManager::setBoundaryLhsIntegrationRule(PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createBoundaryFEPipeline<DIM>(feBoundaryLhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
PipelineManager::setBoundaryLhsIntegrationRule<-1>(PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return setBoundaryLhsIntegrationRule<1>(rule);
  case 2:
    return setBoundaryLhsIntegrationRule<2>(rule);
  case 3:
    return setBoundaryLhsIntegrationRule<3>(rule);
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

template <int DIM>
MoFEMErrorCode PipelineManager::setBoundaryRhsIntegrationRule(PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createBoundaryFEPipeline<DIM>(feBoundaryRhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
PipelineManager::setBoundaryRhsIntegrationRule<-1>(PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return setBoundaryRhsIntegrationRule<1>(rule);
  case 2:
    return setBoundaryRhsIntegrationRule<2>(rule);
  case 3:
    return setBoundaryRhsIntegrationRule<3>(rule);
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

template <int DIM>
MoFEMErrorCode PipelineManager::setSkeletonLhsIntegrationRule(PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createBoundaryFEPipeline<DIM>(feSkeletonLhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
PipelineManager::setSkeletonLhsIntegrationRule<-1>(PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return setSkeletonLhsIntegrationRule<1>(rule);
  case 2:
    return setSkeletonLhsIntegrationRule<2>(rule);
  case 3:
    return setSkeletonLhsIntegrationRule<3>(rule);
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

template <int DIM>
MoFEMErrorCode PipelineManager::setSkeletonRhsIntegrationRule(PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createBoundaryFEPipeline<DIM>(feSkeletonRhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
PipelineManager::setSkeletonRhsIntegrationRule<-1>(PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return setSkeletonRhsIntegrationRule<1>(rule);
  case 2:
    return setSkeletonRhsIntegrationRule<2>(rule);
  case 3:
    return setSkeletonRhsIntegrationRule<3>(rule);
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

template <int DIM>
MoFEMErrorCode PipelineManager::setDomainExplicitRhsIntegrationRule(
    PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feDomainExplicitRhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <int DIM>
MoFEMErrorCode PipelineManager::setBoundaryExplicitRhsIntegrationRule(
    PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createBoundaryFEPipeline<DIM>(feBoundaryExplicitRhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <int DIM>
MoFEMErrorCode PipelineManager::setSkeletonExplicitRhsIntegrationRule(
    PipelineManager::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createBoundaryFEPipeline<DIM>(feSkeletonExplicitRhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <int DIM>
boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpDomainLhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createDomainFEPipeline<DIM>(feDomainLhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpDomainLhsPipeline<-1>() {
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return getOpDomainLhsPipeline<1>();
  case 2:
    return getOpDomainLhsPipeline<2>();
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpDomainLhsPipeline<3>();
}

template <int DIM>
boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpDomainRhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createDomainFEPipeline<DIM>(feDomainRhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpDomainRhsPipeline<-1>() {
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return getOpDomainRhsPipeline<1>();
  case 2:
    return getOpDomainRhsPipeline<2>();
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpDomainRhsPipeline<3>();
}

template <int DIM>
boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpBoundaryLhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feBoundaryLhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpBoundaryLhsPipeline<-1>() {
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return getOpBoundaryLhsPipeline<1>();
  case 2:
    return getOpBoundaryLhsPipeline<2>();
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpBoundaryLhsPipeline<3>();
}

template <int DIM>
boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpBoundaryRhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feBoundaryRhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpBoundaryRhsPipeline<-1>() {
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return getOpBoundaryRhsPipeline<1>();
  case 2:
    return getOpBoundaryRhsPipeline<2>();
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpBoundaryRhsPipeline<3>();
}

template <int DIM>
boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpSkeletonLhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feSkeletonRhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpSkeletonLhsPipeline<-1>() {
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return getOpSkeletonLhsPipeline<1>();
  case 2:
    return getOpSkeletonLhsPipeline<2>();
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpSkeletonLhsPipeline<3>();
}

template <int DIM>
boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpSkeletonRhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feSkeletonLhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpSkeletonRhsPipeline<-1>() {
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return getOpSkeletonRhsPipeline<1>();
  case 2:
    return getOpSkeletonRhsPipeline<2>();
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpSkeletonRhsPipeline<3>();
}

template <int DIM>
boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpDomainExplicitRhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createDomainFEPipeline<DIM>(feDomainExplicitRhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpDomainExplicitRhsPipeline<-1>() {
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return getOpDomainExplicitRhsPipeline<1>();
  case 2:
    return getOpDomainExplicitRhsPipeline<2>();
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpDomainExplicitRhsPipeline<3>();
}

template <int DIM>
boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpBoundaryExplicitRhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feSkeletonExplicitRhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpBoundaryExplicitRhsPipeline<-1>() {
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return getOpBoundaryExplicitRhsPipeline<1>();
  case 2:
    return getOpBoundaryExplicitRhsPipeline<2>();
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpBoundaryExplicitRhsPipeline<3>();
}

template <int DIM>
boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpSkeletonExplicitRhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feSkeletonExplicitRhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<PipelineManager::UserDataOperator> &
PipelineManager::getOpSkeletonExplicitRhsPipeline<-1>() {
  switch (cOre.getInterface<Simple>()->getDim()) {
  case 1:
    return getOpSkeletonExplicitRhsPipeline<1>();
  case 2:
    return getOpSkeletonExplicitRhsPipeline<2>();
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpSkeletonExplicitRhsPipeline<3>();
}

} // namespace MoFEM

#endif // __BASIC_HPP__

/**
 * \defgroup mofem_basic_interface PipelineManager interface
 * \brief Implementation of basic interface for rapid problem implementation.
 *
 * \ingroup mofem
 **/
