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
  getOpDomainLhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Domain Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  boost::ptr_vector<UserDataOperator> &
  getOpDomainRhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Boundary Lhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset  If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  boost::ptr_vector<UserDataOperator> &
  getOpBoundaryLhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Boundary Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  boost::ptr_vector<UserDataOperator> &
  getOpBoundaryRhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Skeleton Lhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  boost::ptr_vector<UserDataOperator> &
  getOpSkeletonLhsPipeline(const bool reset = false);

  /**
   * @brief Get the Op Skeleton Rhs Pipeline object
   * @ingroup mofem_basic_interface
   *
   * @param reset  If true reset pipeline (reset finite element pointer)
   * @return boost::ptr_vector<UserDataOperator>&
   */
  boost::ptr_vector<UserDataOperator> &
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

} // namespace MoFEM

#endif // __BASIC_HPP__

/**
 * \defgroup mofem_basic_interface Simple interface
 * \brief Implementation of basic interface for rapid problem implementation.
 *
 * \ingroup mofem
 **/
