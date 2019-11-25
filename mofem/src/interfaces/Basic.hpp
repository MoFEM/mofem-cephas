/** \file Basic.hpp
 * \brief Header file for basic interface
 * \ingroup mofem_basic_interface
 *
 * Make simplified interface, to speedup problem setup and analysts.
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

  using Simple::Simple::Simple;

  using UserDataOperator = MoFEM::ForcesAndSourcesCore::UserDataOperator;

  /**
   * @brief Get the Op Domain Lhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  boost::ptr_vector<UserDataOperator> &
  inline getOpDomainLhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Domain Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  boost::ptr_vector<UserDataOperator> &
  inline getOpDomainRhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Boundary Lhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset  If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  boost::ptr_vector<UserDataOperator> &
  inline getOpBoundaryLhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Boundary Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  boost::ptr_vector<UserDataOperator> &
  inline getOpBoundaryRhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Skeleton Lhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  boost::ptr_vector<UserDataOperator> &
  inline getOpSkeletonLhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Skeleton Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset  If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  boost::ptr_vector<UserDataOperator> &
  inline getOpSkeletonRhsPipeline(const bool reset = false);

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

  boost::shared_ptr<ForcesAndSourcesCore>
  createDomainFEPipeline(boost::shared_ptr<ForcesAndSourcesCore> &fe,
                         const bool reset = false);

  boost::shared_ptr<ForcesAndSourcesCore>
  createBoundaryFEPipeline(boost::shared_ptr<ForcesAndSourcesCore> &fe,
                           const bool reset = false);

  boost::shared_ptr<ForcesAndSourcesCore>
      feDomainRhs; ///< Element to assemble RHS side by integrating domain
  boost::shared_ptr<ForcesAndSourcesCore>
      feDomainLhs; ///< Element to assemble LHS side by integrating domain
  boost::shared_ptr<ForcesAndSourcesCore>
      feBcRhs; ///< Element to assemble RHS side by integrating boundary
  boost::shared_ptr<ForcesAndSourcesCore>
      feBcLhs; ///< Element to assemble LHS side by integrating boundary
  boost::shared_ptr<ForcesAndSourcesCore>
      feSkeletonRhs; ///< Element to assemble RHS side by integrating skeleton
  boost::shared_ptr<ForcesAndSourcesCore>
      feSkeletonLhs; ///< Element to assemble LHS side by integrating skeleton
};

boost::shared_ptr<ForcesAndSourcesCore>
Basic::createDomainFEPipeline(boost::shared_ptr<ForcesAndSourcesCore> &fe,
                              const bool reset) {
  if (!fe || reset) {
    switch (getDim()) {
    case 2:
      fe = boost::make_shared<FaceElementForcesAndSourcesCore>(cOre);
      break;
    case 3:
      fe = boost::make_shared<VolumeElementForcesAndSourcesCore>(cOre);
      break;
    default:
      THROW_MESSAGE("Dimension not implemented Dim = " +
                    boost::lexical_cast<std::string>(getDim()));
    }
  }
  return fe;
}

boost::shared_ptr<ForcesAndSourcesCore>
Basic::createBoundaryFEPipeline(boost::shared_ptr<ForcesAndSourcesCore> &fe,
                                const bool reset) {
  if (!fe || reset) {
    switch (getDim()) {
    case 2:
      fe = boost::make_shared<EdgeElementForcesAndSourcesCore>(cOre);
      break;
    case 3:
      fe = boost::make_shared<FaceElementForcesAndSourcesCore>(cOre);
      break;
    default:
      THROW_MESSAGE("Dimension not implemented Dim = " +
                    boost::lexical_cast<std::string>(getDim()));
    }
  }
  return fe;
}

boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainLhsPipeline(const bool reset) {
  return createDomainFEPipeline(feDomainLhs, reset)->getOpPtrVector();
}

boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainRhsPipeline(const bool reset) {
  return createDomainFEPipeline(feDomainRhs, reset)->getOpPtrVector();
}

boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryLhsPipeline(const bool reset) {
  return createBoundaryFEPipeline(feBcLhs, reset)->getOpPtrVector();
}

boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryRhsPipeline(const bool reset) {
  return createBoundaryFEPipeline(feBcRhs, reset)->getOpPtrVector();
}

boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonLhsPipeline(const bool reset) {
  return createBoundaryFEPipeline(feSkeletonLhs, reset)->getOpPtrVector();
}

boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonRhsPipeline(const bool reset) {
  return createBoundaryFEPipeline(feSkeletonRhs, reset)->getOpPtrVector();
}

MoFEMErrorCode Basic::loopFiniteElements() {
  MoFEMFunctionBegin;

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getDomainFEName(), feDomainLhs);
  if (feBcLhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getBoundaryFEName(), feBcLhs);
  if (feSkeletonLhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getSkeletonFEName(),
                                    feSkeletonLhs);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getDomainFEName(), feDomainRhs);
  if (feBcRhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getBoundaryFEName(), feBcRhs);
  if (feSkeletonRhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getSkeletonFEName(),
                                    feSkeletonRhs);

  MoFEMFunctionReturn(0);
}

SmartPetscObj<KSP> Basic::createKSP() {
  Interface &m_field = cOre;

  boost::shared_ptr<KspCtx> snes_ctx(new KspCtx(m_field, getProblemName()));
  CHKERR DMMoFEMSetKspCtx(getDM(), snes_ctx);

  boost::shared_ptr<FEMethod> null;

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMMoFEMKSPSetComputeOperators(getDM(), getDomainFEName(),
                                         feDomainLhs, null, null);
  if (feBcLhs)
    CHKERR DMMoFEMKSPSetComputeOperators(getDM(), getBoundaryFEName(), feBcLhs,
                                         null, null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMKSPSetComputeOperators(getDM(), getSkeletonFEName(),
                                         feSkeletonLhs, null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMKSPSetComputeRHS(getDM(), getDomainFEName(), feDomainRhs,
                                   null, null);
  if (feBcRhs)
    CHKERR DMMoFEMKSPSetComputeRHS(getDM(), getBoundaryFEName(), feBcRhs, null,
                                   null);
  if (feSkeletonRhs)
    CHKERR DMMoFEMKSPSetComputeRHS(getDM(), getSkeletonFEName(), feSkeletonRhs,
                                   null, null);

  auto ksp = MoFEM::createKSP(m_field.get_comm());
  CHKERR KSPSetDM(ksp, getDM());
  return ksp;
}

SmartPetscObj<SNES> Basic::createSNES() {
  Interface &m_field = cOre;

  boost::shared_ptr<MoFEM::SnesCtx> snes_ctx(
      new SnesCtx(m_field, getProblemName()));
  CHKERR DMMoFEMSetSnesCtx(getDM(), snes_ctx);

  boost::shared_ptr<FEMethod> null;

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMMoFEMSNESSetJacobian(getDM(), getDomainFEName(), feDomainLhs, null,
                                  null);
  if (feBcLhs)
    CHKERR DMMoFEMSNESSetJacobian(getDM(), getBoundaryFEName(), feBcLhs, null,
                                  null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMSNESSetJacobian(getDM(), getSkeletonFEName(), feSkeletonLhs,
                                  null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMSNESSetFunction(getDM(), getDomainFEName(), feDomainRhs, null,
                                  null);
  if (feBcRhs)
    CHKERR DMMoFEMSNESSetFunction(getDM(), getBoundaryFEName(), feBcRhs, null,
                                  null);
  if (feSkeletonRhs)
    CHKERR DMMoFEMSNESSetFunction(getDM(), getSkeletonFEName(), feSkeletonRhs,
                                  null, null);

  auto snes = MoFEM::createSNES(m_field.get_comm());
  CHKERR SNESSetDM(snes, getDM());
  return snes;
}

SmartPetscObj<TS> Basic::createTS() {
  Interface &m_field = cOre;

  boost::shared_ptr<MoFEM::TsCtx> ts_ctx(new TsCtx(m_field, getProblemName()));
  CHKERR DMMoFEMSetTsCtx(getDM(), ts_ctx);

  boost::shared_ptr<FEMethod> null;

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMMoFEMTSSetIJacobian(getDM(), getDomainFEName(), feDomainLhs, null,
                                 null);
  if (feBcLhs)
    CHKERR DMMoFEMTSSetIJacobian(getDM(), getBoundaryFEName(), feBcLhs, null,
                                 null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMTSSetIJacobian(getDM(), getSkeletonFEName(), feSkeletonLhs,
                                 null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMTSSetIFunction(getDM(), getDomainFEName(), feDomainRhs, null,
                                 null);
  if (feBcRhs)
    CHKERR DMMoFEMTSSetIFunction(getDM(), getBoundaryFEName(), feBcRhs, null,
                                 null);
  if (feSkeletonRhs)
    CHKERR DMMoFEMTSSetIFunction(getDM(), getSkeletonFEName(), feSkeletonRhs,
                                 null, null);

  // Note: More cases for explit, and implicit time ingeration cases can be
  // implemented here.

  auto ts = MoFEM::createTS(m_field.get_comm());
  CHKERR TSSetDM(ts, getDM());
  return ts;
}

} // namespace MoFEM

#endif // __BASIC_HPP__

/**
 * \defgroup mofem_basic_interface Simple interface
 * \brief Implementation of basic interface for rapid problem implementation.
 *
 * \ingroup mofem
 **/
