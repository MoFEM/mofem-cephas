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

/** \brief Field evaluator interface

  * \ingroup field_evaluator
  */
struct FieldEvaluatorInterface : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  FieldEvaluatorInterface(const MoFEM::Core &core);

  struct SetPtsData {

    SetPtsData() = delete;

    /**
     * @brief Set the Gauss Pts data
     *
     * @param fe_method_ptr pointer to finite element instance
     * @param eval_points pointer to array with evaluation points
     * @param nb_eval_points number of evaluated points
     * @param eps tolerance used to find if point is in the element
     * @param verb
     */
    SetPtsData(boost::shared_ptr<MoFEM::ForcesAndSourcesCore> fe_method_ptr,
               const double *eval_points, const int nb_eval_points,
               const double eps, VERBOSITY_LEVELS verb = QUIET)
        : feMethodPtr(fe_method_ptr), evalPoints(eval_points),
          nbEvalPoints(nb_eval_points), eps(eps), verb(verb) {
      localCoords.resize(nbEvalPoints, 3);
      shapeFunctions.resize(nbEvalPoints, 4);
    }

    inline void setEvalPoints(const double *ptr, const int nb_eval_points) {
      evalPoints = ptr;
      nbEvalPoints = nb_eval_points;
      localCoords.resize(nbEvalPoints, 3, false);
      shapeFunctions.resize(nbEvalPoints, 4, false);
    }

    boost::weak_ptr<MoFEM::ForcesAndSourcesCore> feMethodPtr;
    const double *evalPoints;
    int nbEvalPoints;
    double eps;
    VERBOSITY_LEVELS verb;

    MatrixDouble localCoords;
    MatrixDouble shapeFunctions;
    std::vector<EntityHandle> evalPointEntityHandle;

    EntityHandle rooTreeSet;
    boost::scoped_ptr<AdaptiveKDTree> treePtr;
  };

  /**
   * @brief Default evaluator for setting integration points
   *
   */
  struct SetPts {
    SetPts() = delete;
    SetPts(boost::shared_ptr<SetPtsData> data_ptr) : dataPtr(data_ptr) {}
    MoFEMErrorCode operator()(ForcesAndSourcesCore *fe_raw_ptr, int order_row,
                              int order_col, int order_data);

  private:
    boost::weak_ptr<SetPtsData> dataPtr;
  };

  /**
   * @brief Get the Data object
   *
   * Pack pointers with data structures for field evaluator and finite
   * element. Function return shared pointer if returned shared pointer
   * is reset; all data are destroyed. The idea is to pack all data in
   * one structure, create shared pointer to it and return aliased shared
   * pointer to one of the elements of the data structure. It is a bit
   * complicated, but has great flexibility.
   *
   * @tparam VE
   * @tparam SetPtsData
   * @tparam SetPts
   * @param ptr
   * @param nb_eval_points
   * @param eps
   * @param verb
   * @return boost::shared_ptr<SPD>
   */
  template <typename VE, typename SPD = SetPtsData, typename SP = SetPts>
  boost::shared_ptr<SPD>
  getData(const double *ptr = nullptr, const int nb_eval_points = 0,
          const double eps = 1e-12, VERBOSITY_LEVELS verb = QUIET) {
    struct PackData {
      boost::scoped_ptr<VE> elePtr;
      boost::scoped_ptr<SPD> setPtsDataPtr;
      boost::scoped_ptr<SP> setPtsPtr;
    };
    boost::shared_ptr<PackData> pack_data(new PackData());
    MoFEM::Interface &m_field = cOre;
    pack_data->elePtr.reset(new VE(m_field));
    pack_data->setPtsDataPtr.reset(
        new SPD(boost::shared_ptr<VE>(pack_data, pack_data->elePtr.get()), ptr,
                nb_eval_points, eps, verb));
    pack_data->setPtsPtr.reset(new SP(
        boost::shared_ptr<SPD>(pack_data, pack_data->setPtsDataPtr.get())));
    pack_data->elePtr->setRuleHook = boost::ref(*pack_data->setPtsPtr);
    boost::shared_ptr<SPD> data(pack_data, pack_data->setPtsDataPtr.get());
    return data;
  }

  /**
   * @brief Build spatial tree
   *
   * @param finite_element finite element name
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode buildTree3D(boost::shared_ptr<SetPtsData> spd_ptr,
                             const std::string finite_element);

  /**
   * @copydoc buildTree3D
   */
  MoFEMErrorCode buildTree2D(boost::shared_ptr<SetPtsData> spd_ptr,
                             const std::string finite_element);

  /**
   * @brief Evaluate field at artbitray position
   *
   * \code

    std::array<double, 3> point = {0, 0, 0};
    const double dist = 0.3;
    std::array<double, 6> eval_points = {-1.,  -1., -1., 1., 1., 1. };

    using VolEle = VolumeElementForcesAndSourcesCore;
    auto data_ptr =
      m_field.getInterface<FieldEvaluatorInterface>()->
      getData<VolEle>(point.data(), point.size/3);

    if(auto vol_ele = data_ptr->feMethod.lock()) {
      // push operators to finite element instance, e.g.
      vol_ele->getOpPtrVector().push_back(new MyOp());
      // iterate over elemnts with evaluated points
      auto cache_ptr = boost::make_shared<CacheTuple>();
      CHKERR m_field.cache_problem_entities(prb_ptr->getName(), cache_ptr);
      CHKERR m_field.getInterface<FieldEvaluatorInterface>()
        ->evalFEAtThePoint3D(point.data(), dist, prb_ptr->getName(),
                             "FINITE_ELEMENT_NAME",
                             data_ptr, m_field.get_comm_rank(),
                             m_field.get_comm_rank(), cache_ptr);
    }

   * \endcode
   *
   * @param point point used to find tetrahedrons
   * @param distance distance from the point where tetrahedrons are searched
   * @param problem problem name
   * @param finite_element finite element name
   * @param data_ptr pointer to data abut gauss points
   * @param lower_rank lower processor
   * @param upper_rank upper process
   * @param cache_ptr cache or problem entities
   * @param bh control looping over entities, e.g. throwing error if element not
   found
   * @param verb verbosity level
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode
  evalFEAtThePoint3D(const double *const point, const double distance,
                     const std::string problem,
                     const std::string finite_element,
                     boost::shared_ptr<SetPtsData> data_ptr, int lower_rank,
                     int upper_rank, boost::shared_ptr<CacheTuple> cache_ptr,
                     MoFEMTypes bh = MF_EXIST, VERBOSITY_LEVELS verb = QUIET);

  /**
     * @copydoc evalFEAtThePoint3D
     */
  MoFEMErrorCode
  evalFEAtThePoint2D(const double *const point, const double distance,
                     const std::string problem,
                     const std::string finite_element,
                     boost::shared_ptr<SetPtsData> data_ptr, int lower_rank,
                     int upper_rank, boost::shared_ptr<CacheTuple> cache_ptr,
                     MoFEMTypes bh = MF_EXIST, VERBOSITY_LEVELS verb = QUIET);

private:
  /**
   * @copydoc buildTree3D
   */
  template <int D>
  MoFEMErrorCode buildTree(boost::shared_ptr<SetPtsData> spd_ptr,
                           const std::string finite_element);

  /**
     * @copydoc evalFEAtThePoint3D
     */
  template <int D>
  MoFEMErrorCode
  evalFEAtThePoint(const double *const point, const double distance,
                     const std::string problem,
                     const std::string finite_element,
                     boost::shared_ptr<SetPtsData> data_ptr, int lower_rank,
                     int upper_rank, boost::shared_ptr<CacheTuple> cache_ptr,
                     MoFEMTypes bh = MF_EXIST, VERBOSITY_LEVELS verb = QUIET);

};

} // namespace MoFEM

#endif // __FIELD_EVALUATOR_HPP__

/**
 * \defgroup field_evaluator Field Evaluator
 * \brief Evaluate field at the point
 *
 * \ingroup mofem
 */