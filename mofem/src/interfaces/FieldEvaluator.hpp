/** \file FieldEvaluator.hpp
 * \brief Field Evaluator
 *
 * Evaluate field at given coordinate 
 *
 *
 * \ingroup field_evaluator
 */

/*
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __FIELD_EVALUATOR_HPP__
#define __FIELD_EVALUATOR_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMFieldEvaluator =
    MOFEMuuid(BitIntefaceId(FIELDEVALUATOR_INTERFACE));

/** \brief Field evaluator interface

  * \ingroup field_evaluator
  */
struct FieldEvaluatorInterface : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  FieldEvaluatorInterface(const MoFEM::Core &core)
      : cOre(const_cast<MoFEM::Core &>(core)) {}

  /**
   * @brief Build spatial tree
   * 
   * @param finite_element finite element name
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode buildTree3D(const std::string finite_element);

  /**
   * @brief Default evaluator for setting integration points
   * 
   */
  struct SetGaussPts {

    /**
     * @brief Set the Gauss Pts object
     * 
     * @param fe_method finite element instance
     * @param eval_points pointer to array with evaluation points
     * @param nb_eval_points number of evaluated points
     * @param eps tolerance used to find if point is in the element
     * @param verb 
     */
    SetGaussPts(boost::shared_ptr<MoFEM::ForcesAndSourcesCore> fe_method,
                const double *eval_points, const int nb_eval_points,
                const double eps, VERBOSITY_LEVELS verb = QUIET)
        : feMethod(fe_method), evalPoints(eval_points),
          nbEvalPoints(nb_eval_points), eps(eps), verb(verb) {
      localCoords.resize(nbEvalPoints, 3);
      shapeFunctions.resize(nbEvalPoints, 4);
    }

    MoFEMErrorCode operator()(int order_row, int order_col, int order_data);

    inline void setEvalPoints(const double *ptr, const int nb_eval_points) {
      evalPoints = ptr;
      nbEvalPoints = nb_eval_points;
    }

  private:
    boost::shared_ptr<MoFEM::ForcesAndSourcesCore> feMethod;
    const double *evalPoints;
    int nbEvalPoints;
    double eps;
    VERBOSITY_LEVELS verb;

    MatrixDouble localCoords;
    MatrixDouble shapeFunctions;
  };

  /**
   * @brief Evaluate field at artbitray position
   *
   * \code

    std::array<double, 3> point = {0, 0, 0};
    const double dist = 0.3;
    std::array<double, 6> eval_points = {-1.,  -1., -1., 1., 1., 1. };

    boost::shared_ptr<VolumeElementForcesAndSourcesCore> vol_ele(
          new VolumeElementForcesAndSourcesCore(m_field));

    // push operators to finite element instance, e.g.
    vol_ele->getOpPtrVector().push_back(new MyOp());

    // use default evaluator for gauss points
    boost::shared_ptr<FieldEvaluatorInterface::SetGaussPts> set_gauss_pts(
          new FieldEvaluatorInterface::SetGaussPts(
              vol_ele, &eval_points[0], eval_points.size() / 3, 1e-12));

    // iterate over elemnts with evaluated points
    CHKERR m_field.getInterface<FieldEvaluatorInterface>()
          ->evalFEAtThePoint3D(&point[0], dist, prb_ptr->getName(),
                               simple_interface->getDomainFEName(), vol_ele,
                               set_gauss_pts, m_field.get_comm_rank(),
                               m_field.get_comm_rank(), MF_EXIST, QUIET);

   * \endcode
   *
   * @param point point used to find tetrahedrons
   * @param distance distance from the point where tetrahedrons are searched
   * @param problem problem name
   * @param finite_element finite element name
   * @param fe_method method to evaluate field quantities
   * @param set_gauss_pts  method to set gauss points
   * @param lower_rank lower processor
   * @param upper_rank upper process
   * @param bh control looping over entities, e.g. throwing error if element not
   found
   * @param verb verbosity level
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode
  evalFEAtThePoint3D(const double *const point, const double distance,
                      const std::string problem,
                     const std::string finite_element,
                     boost::shared_ptr<MoFEM::ForcesAndSourcesCore> fe_method, 
                     int lower_rank,
                     int upper_rank, MoFEMTypes bh = MF_EXIST,
                     VERBOSITY_LEVELS verb = QUIET);

private:

  EntityHandle rooTreeSet;
  boost::shared_ptr<BVHTree> treePtr;

};

}

#endif // __FIELD_EVALUATOR_HPP__


/**
 * \defgroup field_evaluator Field Evaluator
 * \brief Evaluate field at the point
 *
 * \ingroup mofem
 */