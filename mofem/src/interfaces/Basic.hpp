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
  template <int DIM>
  inline boost::ptr_vector<UserDataOperator> &
  getOpDomainLhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Domain Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM>
  inline boost::ptr_vector<UserDataOperator> &
  getOpDomainRhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Boundary Lhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset  If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM>
  inline boost::ptr_vector<UserDataOperator> &
  getOpBoundaryLhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Boundary Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM>
  inline boost::ptr_vector<UserDataOperator> &
  getOpBoundaryRhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Skeleton Lhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM>
  inline boost::ptr_vector<UserDataOperator> &
  getOpSkeletonLhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Skeleton Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset  If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  template <int DIM>
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
      feBcRhs; ///< Element to assemble RHS side by integrating boundary
  boost::shared_ptr<FEMethod>
      feBcLhs; ///< Element to assemble LHS side by integrating boundary
  boost::shared_ptr<FEMethod>
      feSkeletonRhs; ///< Element to assemble RHS side by integrating skeleton
  boost::shared_ptr<FEMethod>
      feSkeletonLhs; ///< Element to assemble LHS side by integrating skeleton

  template <int DIM>
  inline boost::shared_ptr<FEMethod>
  createDomainFEPipeline(boost::shared_ptr<FEMethod> &fe,
                         const bool reset = false);

  template <int DIM>
  inline boost::shared_ptr<FEMethod>
  createBoundaryFEPipeline(boost::shared_ptr<FEMethod> &fe,
                           const bool reset = false);
};

template <int DIM>
boost::shared_ptr<FEMethod>
Basic::createDomainFEPipeline(boost::shared_ptr<FEMethod> &fe,
                              const bool reset) {
  static_assert(DIM == 1 || DIM == 2 || DIM == 3, "not implemented");
  return boost::make_shared<FEMethod>();
}

template <>
boost::shared_ptr<FEMethod>
inline Basic::createDomainFEPipeline<3>(boost::shared_ptr<FEMethod> &fe,
                                 const bool reset) {
  if (!fe || reset)
    fe = boost::make_shared<VolumeElementForcesAndSourcesCore>(cOre);
  return fe;
}

template <>
boost::shared_ptr<FEMethod>
inline Basic::createDomainFEPipeline<2>(boost::shared_ptr<FEMethod> &fe,
                                 const bool reset) {
  if (!fe || reset)
    fe = boost::make_shared<FaceElementForcesAndSourcesCore>(cOre);
  return fe;
}

template <>
boost::shared_ptr<FEMethod>
inline Basic::createDomainFEPipeline<1>(boost::shared_ptr<FEMethod> &fe,
                                 const bool reset) {
  if (!fe || reset)
    fe = boost::make_shared<EdgeElementForcesAndSourcesCore>(cOre);
  return fe;
}

template <int DIM>
boost::shared_ptr<FEMethod>
Basic::createBoundaryFEPipeline(boost::shared_ptr<FEMethod> &fe,
                                const bool reset) {
  static_assert(DIM == 1 || DIM == 2 || DIM == 3, "not implemented");
  return boost::make_shared<FEMethod>();
}

template <>
inline boost::shared_ptr<FEMethod>
Basic::createBoundaryFEPipeline<3>(boost::shared_ptr<FEMethod> &fe,
                                   const bool reset) {
  if (!fe || reset)
    fe = boost::make_shared<FaceElementForcesAndSourcesCore>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod>
Basic::createBoundaryFEPipeline<2>(boost::shared_ptr<FEMethod> &fe,
                                   const bool reset) {
  if (!fe || reset)
    fe = boost::make_shared<EdgeElementForcesAndSourcesCore>(cOre);
  return fe;
}

template <>
inline boost::shared_ptr<FEMethod>
Basic::createBoundaryFEPipeline<1>(boost::shared_ptr<FEMethod> &fe,
                                   const bool reset) {
  if (!fe || reset)
    fe = boost::make_shared<VertexElementForcesAndSourcesCore>(cOre);
  return fe;
}

template <int DIM>
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainLhsPipeline(const bool reset) {
  return dynamic_cast<boost::shared_ptr<ForcesAndSourcesCore>>(
             createDomainFEPipeline<DIM>(feDomainLhs, reset))
      ->getOpPtrVector();
}

template <int DIM>
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainRhsPipeline(const bool reset) {
  return dynamic_cast<boost::shared_ptr<ForcesAndSourcesCore>>(
             createDomainFEPipeline<DIM>(feDomainRhs, reset))
      ->getOpPtrVector();
}

template <int DIM>
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryLhsPipeline(const bool reset) {
  return dynamic_cast<boost::shared_ptr<ForcesAndSourcesCore>>(
             createBoundaryFEPipeline<DIM>(feBcRhs, reset))
      ->getOpPtrVector();
}

template <int DIM>
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryRhsPipeline(const bool reset) {
  return dynamic_cast<boost::shared_ptr<ForcesAndSourcesCore>>(
             createBoundaryFEPipeline<DIM>(feBcLhs, reset))
      ->getOpPtrVector();
}

template <int DIM>
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonLhsPipeline(const bool reset) {
  return dynamic_cast<boost::shared_ptr<ForcesAndSourcesCore>>(
             createBoundaryFEPipeline<DIM>(feSkeletonRhs, reset))
      ->getOpPtrVector();
}

template <int DIM>
boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonRhsPipeline(const bool reset) {
  return dynamic_cast<boost::shared_ptr<ForcesAndSourcesCore>>(
             createBoundaryFEPipeline<DIM>(feSkeletonLhs, reset))
      ->getOpPtrVector();
}

} // namespace MoFEM

#endif // __BASIC_HPP__

/**
 * \defgroup mofem_basic_interface Simple interface
 * \brief Implementation of basic interface for rapid problem implementation.
 *
 * \ingroup mofem
 **/
