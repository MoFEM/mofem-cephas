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
struct Basic : public MoFEM::Simple {

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

  template <typename T = ForcesAndSourcesCore>
  inline boost::shared_ptr<T> &getOpDomainLhsFE();

  template <typename T = ForcesAndSourcesCore>
  inline boost::shared_ptr<T> &getOpDomainRhsFE();

  template <typename T = ForcesAndSourcesCore>
  inline boost::shared_ptr<T> &getOpBoundaryLhsFE();

  template <typename T = ForcesAndSourcesCore>
  inline boost::shared_ptr<T> &getOpBoundaryRhsFE();

  template <typename T = ForcesAndSourcesCore>
  inline boost::shared_ptr<T> &getOpSkeletonLhsFE();

  template <typename T = ForcesAndSourcesCore>
  inline boost::shared_ptr<T> &getOpSkeletonRhsFE();

  template <int DIM = -1>
  inline MoFEMErrorCode setDomainLhsIntegrationRule(RuleHookFun rule,
                                                    const bool reset = 0);

  template <int DIM = -1>
  inline MoFEMErrorCode setDomainRhsIntegrationRule(RuleHookFun rule,
                                                    const bool reset = 0);

  template <int DIM = -1>
  inline MoFEMErrorCode setBoundaryLhsIntegrationRule(RuleHookFun rule,
                                                      const bool reset = 0);

  template <int DIM = -1>
  inline MoFEMErrorCode setBoundaryRhsIntegrationRule(RuleHookFun rule,
                                                      const bool reset = 0);

  template <int DIM = -1>
  inline MoFEMErrorCode setSkeletonLhsIntegrationRule(RuleHookFun rule,
                                                      const bool reset = 0);

  template <int DIM = -1>
  inline MoFEMErrorCode setSkeletonRhsIntegrationRule(RuleHookFun rule,
                                                      const bool reset = 0);

  /**
   * @brief Get the Op Domain Lhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @tparam -1
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &
  getOpDomainLhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Domain Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @tparam -1
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &
  getOpDomainRhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Boundary Lhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @tparam -1
   * @param reset  If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &
  getOpBoundaryLhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Boundary Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @tparam -1
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &
  getOpBoundaryRhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Skeleton Lhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &
  getOpSkeletonLhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Skeleton Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @tparam -1
   * @param reset  If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM = -1>
  inline boost::ptr_vector<UserDataOperator> &
  getOpSkeletonRhsPipeline(const bool reset = false);

  /**
   * @brief Iterate finite elements
   * @ingroup mofem_basic_interface
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode loopFiniteElements();

  /**
   * @brief Create KSP (linear) solver
   * @ingroup mofem_basic_interface
   *
   * @return SmartPetscObj<KSP>
   */
  SmartPetscObj<KSP> createKSP();

  /**
   * @brief Create SNES (nonlinear) solver
   * @ingroup mofem_basic_interface
   *
   * @return SmartPetscObj<SNES>
   */
  SmartPetscObj<SNES> createSNES();

  /**
   * @brief Create TS (time) solver
   * @ingroup mofem_basic_interface
   *
   * @return SmartPetscObj<TS>
   */
  SmartPetscObj<TS> createTS();

private:
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
  createDomainFEPipeline(boost::shared_ptr<FEMethod> &fe,
                         const bool reset = false);

  template <int DIM>
  inline boost::shared_ptr<FEMethod> &
  createBoundaryFEPipeline(boost::shared_ptr<FEMethod> &fe,
                           const bool reset = false);
};

template <int DIM>
boost::shared_ptr<FEMethod> &
Basic::createDomainFEPipeline(boost::shared_ptr<FEMethod> &fe,
                              const bool reset) {
  static_assert(DIM == 1 || DIM == 2 || DIM == 3, "not implemented");
  fe = boost::make_shared<FEMethod>();
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
Basic::createDomainFEPipeline<3>(boost::shared_ptr<FEMethod> &fe,
                                 const bool reset) {
  if (!fe || reset)
    fe = boost::make_shared<VolumeElementForcesAndSourcesCore>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
Basic::createDomainFEPipeline<2>(boost::shared_ptr<FEMethod> &fe,
                                 const bool reset) {
  if (!fe || reset)
    fe = boost::make_shared<FaceEle2D>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
Basic::createDomainFEPipeline<1>(boost::shared_ptr<FEMethod> &fe,
                                 const bool reset) {
  if (!fe || reset)
    fe = boost::make_shared<EdgeEle1D>(cOre);
  return fe;
}

template <int DIM>
boost::shared_ptr<FEMethod> &
Basic::createBoundaryFEPipeline(boost::shared_ptr<FEMethod> &fe,
                                const bool reset) {
  static_assert(DIM == 1 || DIM == 2 || DIM == 3, "not implemented");
  fe = boost::make_shared<FEMethod>();
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
Basic::createBoundaryFEPipeline<3>(boost::shared_ptr<FEMethod> &fe,
                                   const bool reset) {
  if (!fe || reset)
    fe = boost::make_shared<FaceElementForcesAndSourcesCore>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
Basic::createBoundaryFEPipeline<2>(boost::shared_ptr<FEMethod> &fe,
                                   const bool reset) {
  if (!fe || reset)
    fe = boost::make_shared<EdgeEle2D>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod> &
Basic::createBoundaryFEPipeline<1>(boost::shared_ptr<FEMethod> &fe,
                                   const bool reset) {
  if (!fe || reset)
    fe = boost::make_shared<VertexElementForcesAndSourcesCore>(cOre);
  return fe;
}

template <typename T> boost::shared_ptr<T> &Basic::getOpDomainLhsFE() {
  return boost::dynamic_pointer_cast<T>(feDomainLhs);
}

template <typename T> boost::shared_ptr<T> &Basic::getOpDomainRhsFE() {
  return boost::dynamic_pointer_cast<T>(feDomainRhs);
}

template <typename T> boost::shared_ptr<T> &Basic::getOpBoundaryLhsFE() {
  return boost::dynamic_pointer_cast<T>(feBoundaryLhs);
}

template <typename T> boost::shared_ptr<T> &Basic::getOpBoundaryRhsFE() {
  return boost::dynamic_pointer_cast<T>(feBoundaryRhs);
}

template <typename T> boost::shared_ptr<T> &Basic::getOpSkeletonLhsFE() {
  return boost::dynamic_pointer_cast<T>(feSkeletonLhs);
}

template <typename T> boost::shared_ptr<T> &Basic::getOpSkeletonRhsFE() {
  return boost::dynamic_pointer_cast<T>(feSkeletonRhs);
}

template <int DIM>
MoFEMErrorCode Basic::setDomainLhsIntegrationRule(Basic::RuleHookFun rule,
                                                  const bool reset) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feDomainLhs, reset))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
Basic::setDomainLhsIntegrationRule<-1>(Basic::RuleHookFun rule,
                                       const bool reset) {
  MoFEMFunctionBegin;
  switch (getDim()) {
  case 1:
    return setDomainLhsIntegrationRule<1>(rule, reset);
  case 2:
    return setDomainLhsIntegrationRule<2>(rule, reset);
  case 3:
    return setDomainLhsIntegrationRule<3>(rule, reset);
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

template <int DIM>
MoFEMErrorCode Basic::setDomainRhsIntegrationRule(Basic::RuleHookFun rule,
                                                  const bool reset) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feDomainRhs, reset))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
Basic::setDomainRhsIntegrationRule<-1>(Basic::RuleHookFun rule,
                                       const bool reset) {
  MoFEMFunctionBegin;
  switch (getDim()) {
  case 1:
    return setDomainRhsIntegrationRule<1>(rule, reset);
  case 2:
    return setDomainRhsIntegrationRule<2>(rule, reset);
  case 3:
    return setDomainRhsIntegrationRule<3>(rule, reset);
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

template <int DIM>
MoFEMErrorCode Basic::setBoundaryLhsIntegrationRule(Basic::RuleHookFun rule,
                                                    const bool reset) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feBoundaryLhs, reset))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
Basic::setBoundaryLhsIntegrationRule<-1>(Basic::RuleHookFun rule,
                                         const bool reset) {
  MoFEMFunctionBegin;
  switch (getDim()) {
  case 1:
    return setBoundaryLhsIntegrationRule<1>(rule, reset);
  case 2:
    return setBoundaryLhsIntegrationRule<2>(rule, reset);
  case 3:
    return setBoundaryLhsIntegrationRule<3>(rule, reset);
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

template <int DIM>
MoFEMErrorCode Basic::setBoundaryRhsIntegrationRule(Basic::RuleHookFun rule,
                                                    const bool reset) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feBoundaryRhs, reset))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
Basic::setBoundaryRhsIntegrationRule<-1>(Basic::RuleHookFun rule,
                                         const bool reset) {
  MoFEMFunctionBegin;
  switch (getDim()) {
  case 1:
    return setBoundaryRhsIntegrationRule<1>(rule, reset);
  case 2:
    return setBoundaryRhsIntegrationRule<2>(rule, reset);
  case 3:
    return setBoundaryRhsIntegrationRule<3>(rule, reset);
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

template <int DIM>
MoFEMErrorCode Basic::setSkeletonLhsIntegrationRule(Basic::RuleHookFun rule,
                                                    const bool reset) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feSkeletonLhs, reset))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
Basic::setSkeletonLhsIntegrationRule<-1>(Basic::RuleHookFun rule,
                                         const bool reset) {
  MoFEMFunctionBegin;
  switch (getDim()) {
  case 1:
    return setSkeletonLhsIntegrationRule<1>(rule, reset);
  case 2:
    return setSkeletonLhsIntegrationRule<2>(rule, reset);
  case 3:
    return setSkeletonLhsIntegrationRule<3>(rule, reset);
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

template <int DIM>
MoFEMErrorCode Basic::setSkeletonRhsIntegrationRule(Basic::RuleHookFun rule,
                                                    const bool reset) {
  MoFEMFunctionBegin;
  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
      createDomainFEPipeline<DIM>(feSkeletonRhs, reset))
      ->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <>
inline MoFEMErrorCode
Basic::setSkeletonRhsIntegrationRule<-1>(Basic::RuleHookFun rule,
                                         const bool reset) {
  MoFEMFunctionBegin;
  switch (getDim()) {
  case 1:
    return setSkeletonRhsIntegrationRule<1>(rule, reset);
  case 2:
    return setSkeletonRhsIntegrationRule<2>(rule, reset);
  case 3:
    return setSkeletonRhsIntegrationRule<3>(rule, reset);
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

template <int DIM>
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainLhsPipeline(const bool reset) {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createDomainFEPipeline<DIM>(feDomainLhs, reset))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainLhsPipeline<-1>(const bool reset) {
  switch (getDim()) {
  case 1:
    return getOpDomainLhsPipeline<1>(reset);
  case 2:
    return getOpDomainLhsPipeline<2>(reset);
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpDomainLhsPipeline<3>(reset);
}

template <int DIM>
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainRhsPipeline(const bool reset) {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createDomainFEPipeline<DIM>(feDomainRhs, reset))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainRhsPipeline<-1>(const bool reset) {
  switch (getDim()) {
  case 1:
    return getOpDomainRhsPipeline<1>(reset);
  case 2:
    return getOpDomainRhsPipeline<2>(reset);
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpDomainRhsPipeline<3>(reset);
}

template <int DIM>
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryLhsPipeline(const bool reset) {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feBoundaryRhs, reset))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryLhsPipeline<-1>(const bool reset) {
  switch (getDim()) {
  case 1:
    return getOpBoundaryLhsPipeline<1>(reset);
  case 2:
    return getOpBoundaryLhsPipeline<2>(reset);
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpBoundaryLhsPipeline<3>(reset);
}

template <int DIM>
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryRhsPipeline(const bool reset) {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feBoundaryLhs, reset))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryRhsPipeline<-1>(const bool reset) {
  switch (getDim()) {
  case 1:
    return getOpBoundaryRhsPipeline<1>(reset);
  case 2:
    return getOpBoundaryRhsPipeline<2>(reset);
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpBoundaryRhsPipeline<3>(reset);
}

template <int DIM>
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonLhsPipeline(const bool reset) {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feSkeletonRhs, reset))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonLhsPipeline<-1>(const bool reset) {
  switch (getDim()) {
  case 1:
    return getOpSkeletonLhsPipeline<1>(reset);
  case 2:
    return getOpSkeletonLhsPipeline<2>(reset);
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpSkeletonLhsPipeline<3>(reset);
}

template <int DIM>
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonRhsPipeline(const bool reset) {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createBoundaryFEPipeline<DIM>(feSkeletonLhs, reset))
      ->getOpPtrVector();
}

template <>
inline boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonRhsPipeline<-1>(const bool reset) {
  switch (getDim()) {
  case 1:
    return getOpSkeletonRhsPipeline<1>(reset);
  case 2:
    return getOpSkeletonRhsPipeline<2>(reset);
  case 3:
    break;
  default:
    THROW_MESSAGE("Not implemented");
  }
  return getOpSkeletonRhsPipeline<3>(reset);
}

} // namespace MoFEM

#endif // __BASIC_HPP__

/**
 * \defgroup mofem_basic_interface Simple interface
 * \brief Implementation of basic interface for rapid problem implementation.
 *
 * \ingroup mofem
 **/
