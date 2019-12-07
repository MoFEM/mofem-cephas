/** \file Basic.hpp
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

static const MOFEMuuid IDD_MOFEMBasic =
    MOFEMuuid(BitIntefaceId(BASIC_INTERFACE));

/**
 * \brief Basic interface
 * \ingroup mofem_basic_interface
 */
struct Basic : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  Basic(const MoFEM::Core &core);

  using UserDataOperator = MoFEM::ForcesAndSourcesCore::UserDataOperator;
  using RuleHookFun = MoFEM::ForcesAndSourcesCore::RuleHookFun;

  using FaceEle2D = MoFEM::FaceElementForcesAndSourcesCoreSwitch<
      FaceElementForcesAndSourcesCore::NO_CONTRAVARIANT_TRANSFORM_HDIV |
      FaceElementForcesAndSourcesCore::NO_COVARIANT_TRANSFORM_HCURL>;
  using EdgeEle2D = MoFEM::EdgeElementForcesAndSourcesCoreSwitch<
      EdgeElementForcesAndSourcesCore::NO_COVARIANT_TRANSFORM_HCURL>;
  using EdgeEle1D = EdgeEle2D;

  inline boost::shared_ptr<FEMethod> &getDomainLhsFE();

  inline boost::shared_ptr<FEMethod> &getDomainRhsFE();

  inline boost::shared_ptr<FEMethod> &getBoundaryLhsFE();

  inline boost::shared_ptr<FEMethod> &getBoundaryRhsFE();

  inline boost::shared_ptr<FEMethod> &getSkeletonLhsFE();

  inline boost::shared_ptr<FEMethod> &getSkeletonRhsFE();

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

  /**
   * @brief Create TS (time) solver
   * @ingroup mofem_basic_interface
   *
   * @param dm
   * @return SmartPetscObj<TS>
   */
  SmartPetscObj<TS> createTS(SmartPetscObj<DM> dm = nullptr);

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

  template <int DIM>
  inline boost::shared_ptr<FEMethod> &
  createDomainFEPipeline(boost::shared_ptr<FEMethod> &fe);

  template <int DIM>
  inline boost::shared_ptr<FEMethod> &
  createBoundaryFEPipeline(boost::shared_ptr<FEMethod> &fe);
};

template <int DIM>
boost::shared_ptr<FEMethod> &
Basic::createDomainFEPipeline(boost::shared_ptr<FEMethod> &fe) {
  static_assert(DIM == 1 || DIM == 2 || DIM == 3, "not implemented");
  fe = boost::make_shared<FEMethod>();
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
Basic::createDomainFEPipeline<3>(boost::shared_ptr<FEMethod> &fe) {
  if (!fe)
    fe = boost::make_shared<VolumeElementForcesAndSourcesCore>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
Basic::createDomainFEPipeline<2>(boost::shared_ptr<FEMethod> &fe) {
  if (!fe)
    fe = boost::make_shared<FaceEle2D>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
Basic::createDomainFEPipeline<1>(boost::shared_ptr<FEMethod> &fe) {
  if (!fe)
    fe = boost::make_shared<EdgeEle1D>(cOre);
  return fe;
}

template <int DIM>
boost::shared_ptr<FEMethod> &
Basic::createBoundaryFEPipeline(boost::shared_ptr<FEMethod> &fe) {
  static_assert(DIM == 1 || DIM == 2 || DIM == 3, "not implemented");
  fe = boost::make_shared<FEMethod>();
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
Basic::createBoundaryFEPipeline<3>(boost::shared_ptr<FEMethod> &fe) {
  if (!fe)
    fe = boost::make_shared<FaceElementForcesAndSourcesCore>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
Basic::createBoundaryFEPipeline<2>(boost::shared_ptr<FEMethod> &fe) {
  if (!fe)
    fe = boost::make_shared<EdgeEle2D>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
Basic::createBoundaryFEPipeline<1>(boost::shared_ptr<FEMethod> &fe) {
  if (!fe)
    fe = boost::make_shared<VertexElementForcesAndSourcesCore>(cOre);
  return fe;
}

boost::shared_ptr<FEMethod> &Basic::getDomainLhsFE() { return feDomainLhs; }

boost::shared_ptr<FEMethod> &Basic::getDomainRhsFE() { return feDomainRhs; }

boost::shared_ptr<FEMethod> &Basic::getBoundaryLhsFE() { return feBoundaryLhs; }

boost::shared_ptr<FEMethod> &Basic::getBoundaryRhsFE() { return feBoundaryRhs; }

boost::shared_ptr<FEMethod> &Basic::getSkeletonLhsFE() { return feSkeletonLhs; }

boost::shared_ptr<FEMethod> &Basic::getSkeletonRhsFE() { return feSkeletonRhs; }

template <int DIM>
MoFEMErrorCode Basic::setDomainLhsIntegrationRule(Basic::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feDomainLhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
Basic::setDomainLhsIntegrationRule<-1>(Basic::RuleHookFun rule) {
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
MoFEMErrorCode Basic::setDomainRhsIntegrationRule(Basic::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feDomainRhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
Basic::setDomainRhsIntegrationRule<-1>(Basic::RuleHookFun rule) {
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
MoFEMErrorCode Basic::setBoundaryLhsIntegrationRule(Basic::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feBoundaryLhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
Basic::setBoundaryLhsIntegrationRule<-1>(Basic::RuleHookFun rule) {
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
MoFEMErrorCode Basic::setBoundaryRhsIntegrationRule(Basic::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feBoundaryRhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
Basic::setBoundaryRhsIntegrationRule<-1>(Basic::RuleHookFun rule) {
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
MoFEMErrorCode Basic::setSkeletonLhsIntegrationRule(Basic::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feSkeletonLhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
Basic::setSkeletonLhsIntegrationRule<-1>(Basic::RuleHookFun rule) {
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
MoFEMErrorCode Basic::setSkeletonRhsIntegrationRule(Basic::RuleHookFun rule) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feSkeletonRhs))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
Basic::setSkeletonRhsIntegrationRule<-1>(Basic::RuleHookFun rule) {
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
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainLhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createDomainFEPipeline<DIM>(feDomainLhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainLhsPipeline<-1>() {
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
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainRhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createDomainFEPipeline<DIM>(feDomainRhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainRhsPipeline<-1>() {
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
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryLhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feBoundaryLhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryLhsPipeline<-1>() {
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
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryRhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feBoundaryRhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryRhsPipeline<-1>() {
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
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonLhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feSkeletonRhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonLhsPipeline<-1>() {
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
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonRhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feSkeletonLhs))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonRhsPipeline<-1>() {
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

} // namespace MoFEM

#endif // __BASIC_HPP__

/**
 * \defgroup mofem_basic_interface Simple interface
 * \brief Implementation of basic interface for rapid problem implementation.
 *
 * \ingroup mofem
 **/
