/** \file FieldBlas.hpp
 * \brief Field basic algebra
 * \ingroup mofem_is_managers
 *
 * Managing problems, build and partitioning.
 *
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

#ifndef __FIELD_BLAS_HPP__
#define __FIELD_BLAS_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

/**
 * \brief Basic algebra on fields
 * \ingroup mofem_field_algebra
 *
 */
struct FieldBlas : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  const MoFEM::Interface &cOre;
  bool dEbug;

  FieldBlas(const MoFEM::Core &core);

  ~FieldBlas() = default;

  /**
   * @brief function to set a field value
   *
   */
  using OneFieldFunctionOnValues =
      boost::function<double(const double field_val)>;

  /**
   * @brief
   *
   */
  using OneFieldFunctionOnEntities = boost::function<MoFEMErrorCode(
      boost::shared_ptr<FieldEntity> field_entity_ptr)>;

  /**
   * @brief
   * \param y_field_value_reference field "y_field" values
   * \param x_field_value
   */
  using TwoFieldFunctionOnValues = boost::function<double(
      const double y_field_value_reference, const double x_field_value)>;

  /**
   * @brief
   * \param y_field_entity_ptr pointer to second "y_field" field entities
   * \param x_field_entity_ptr pointer to first "x_field" field entities
   */
  using TwoFieldFunctionOnEntities = boost::function<MoFEMErrorCode(
      boost::shared_ptr<FieldEntity> y_field_entity_ptr,
      const boost::shared_ptr<FieldEntity> x_field_entity_ptr)>;

  /** \brief filed lambda
     * \ingroup mofem_field_algebra
     *
     * Do calculation on one field using lambda function with entity field input
     *
     * \code
     auto field_abs = [&](
      boost::shared_ptr<FieldEntity> ent_ptr) {
        MoFEMFunctionBeginHot;
        auto field_data = ent_ptr->getEntFieldData();
        for (auto &v : field_data)
          v = std::abs(v);
        MoFEMFunctionReturnHot(0);
      };

     CHKERR
     m_field.getInterface<FieldBlas>()->fieldLambdaOnEntities(field_abs,
     "U", block_ents_ptr);
     * \endcode
     *
     * \param function f(double &x, double)
     * \param field_name name of field
     * \param ents_ptr execute lambda function only on entities given by pointer
     to range
     *
     */
  MoFEMErrorCode fieldLambdaOnEntities(OneFieldFunctionOnEntities lambda,
                                       const std::string field_name,
                                       Range *ents_ptr = nullptr);

  /** \brief filed lambda
   * \ingroup mofem_field_algebra
   *
   * Do calculation on one field using lambda function with field value input
   *
   * \code
     auto scale_field = [&](const double val) {
       double time = ts_t; // global variable
       return val * time;
      };
     CHKERR
     m_field.getInterface<FieldBlas>()->fieldLambdaOnValues(scale_field,
     "U", false, block_ents_ptr);
   * \endcode
   *
   * \param function f(const double)
   * \param field_name name of field
   * \param ents_ptr execute lambda function only on entities given by
   pointer to range
   *
   */
  MoFEMErrorCode fieldLambdaOnValues(OneFieldFunctionOnValues lambda,
                                     const std::string field_name,
                                     Range *ents_ptr = nullptr);

  /** \brief filed lambda
   * \ingroup mofem_field_algebra
   *
   * Do calculation on two fields and save result to field fy
   * input function usees field values
   *   \code
      auto field_axpy = [&](const double val_y,
        const double val_x) {
        // y = 2 * x + y
        return 2 * val_x + val_y;
      };
     CHKERR
     m_field.getInterface<FieldBlas>()->fieldLambdaOnValues(field_axpy,
     "U", "P" false, block_ents_ptr);
   *  \endcode
   *
   * \param function f(double &x, double)
   * \param field_name_x name of field_x
   * \param field_name_y name of field_y
   * \param error_if_missing throw error if missing entities of field y
   * \param ents_ptr execute lambda function only on entities given by pointer
   to range
   *
   */
  MoFEMErrorCode fieldLambdaOnValues(TwoFieldFunctionOnValues lambda,
                                     const std::string field_name_x,
                                     const std::string field_name_y,
                                     bool error_if_missing = false,
                                     Range *ents_ptr = nullptr);

  /** \brief filed lambda
   * \ingroup mofem_field_algebra
   *
   * Do calculation on two fields and save result to field fy
   * input function usees field entities
   *
   *  \code
     auto vector_times_scalar_field = [&](boost::shared_ptr<FieldEntity> ent_ptr_y,
        boost::shared_ptr<FieldEntity> ent_ptr_x) {
        MoFEMFunctionBeginHot;
        auto x_data = ent_ptr_x->getEntFieldData();
        auto y_data = ent_ptr_y->getEntFieldData();
        const auto size_x = x_data.size(); // scalar
        const auto size_y = y_data.size(); // vector

        for (size_t dd = 0;; dd != size_y; ++dd)
          y_data[dd] *= x_data[0];
    
        MoFEMFunctionReturnHot(0);
      };

     CHKERR
     m_field.getInterface<FieldBlas>()->fieldLambdaOnValues(vector_times_scalar_field,
          "U", "P" false, block_ents_ptr);

   *  \endcode
   *
   * \param function f(double &x, double)
   * \param field_name_x name of field_x
   * \param field_name_y name of field_y
   * \param error_if_missing throw error if missing entities of field y
   * \param ents_ptr execute lambda function only on entities given by pointer
   to range
   *
   */
  MoFEMErrorCode fieldLambdaOnEntities(TwoFieldFunctionOnEntities lambda,
                                       const std::string field_name_x,
                                       const std::string field_name_y,
                                       bool error_if_missing = false,
                                       Range *ents_ptr = nullptr);

  /**
   * @deprecated This is with obsolete last option
   *
   * @param lambda
   * @param field_name_x
   * @param field_name_y
   * @param error_if_missing
   * @param create_if_missing
   * @return DEPRECATED
   */
  DEPRECATED MoFEMErrorCode fieldLambda(TwoFieldFunctionOnValues lambda,
                                        const std::string field_name_x,
                                        const std::string field_name_y,
                                        bool error_if_missing,
                                        bool create_if_missing);

  /** \brief axpy fields
   * \ingroup mofem_field_algebra
   *
   * field_y = field_y + alpha*field_x
   *
   * \param alpha
   * \param field_name_x name of field_x
   * \param field_name_y name of field_y
   * \param error_if_missing throw error if entity/dof exist in field_x but not
   * \param ents_ptr execute lambda function only on entities given by pointer
   * on field_y
   *
   */
  MoFEMErrorCode fieldAxpy(const double alpha, const std::string field_name_x,
                           const std::string field_name_y,
                           bool error_if_missing = false,
                           Range *ents_ptr = nullptr);

  /** \deprecated axpy fields
   * \ingroup mofem_field_algebra
   *
   * field_y = field_y + alpha*field_x
   *
   * \param alpha
   * \param field_name_x name of field_x
   * \param field_name_y name of field_y
   * \param error_if_missing throw error if entity/dof exist in field_x but not
   * on field_y \param create_if_missing creat dof in field_y from field_x if it
   * is not database
   *
   */
  DEPRECATED MoFEMErrorCode fieldAxpy(const double alpha,
                                      const std::string field_name_x,
                                      const std::string field_name_y,
                                      bool error_if_missing,
                                      bool create_if_missing);

  /** \brief copy and scale fields
   * \ingroup mofem_field_algebra
   *
   * field_y = alpha*field_x
   *
   * \param alpha
   * \param field_name_x name of field_x
   * \param field_name_y name of field_y
   * \param error_if_missing throw error if entity/dof exist in field_x but not
   * \param ents_ptr execute lambda function only on entities given by pointer
   * on field_y \param create_if_missing creat dof in field_y from field_x if it
   * is not database
   *
   */
  MoFEMErrorCode fieldCopy(const double alpha, const std::string field_name_x,
                           const std::string field_name_y,
                           bool error_if_missing = false,
                           Range *ents_ptr = nullptr);

  DEPRECATED MoFEMErrorCode fieldCopy(const double alpha,
                                      const std::string field_name_x,
                                      const std::string field_name_y,
                                      bool error_if_missing,
                                      bool create_if_missing);

  /** \brief field scale
   * \ingroup mofem_field_algebra
   * 
   * \param scale field by value
   * \param field_name name of field
   * \param error_if_missing throw error if missing entities of field y
   * \param ents_ptr execute lambda function only on entities given by
   pointer to range
   */
  MoFEMErrorCode fieldScale(const double alpha, const std::string field_name,
                            Range *ents_ptr = nullptr);

  using VertexCoordsFunction =
      boost::function<MoFEMErrorCode(VectorAdaptor &&field_data, double *xcoord,
                                     double *ycoord, double *zcoord)>;

  /** \brief Set DOFs on vertices using user function
   * \ingroup mofem_field_algebra
   *
   * Example:
   *
   * \code
   * auto do_something = [&](VectorAdaptor &field_data, double *x,
   *                         double *y, double *z) {
   *   MoFEMFunctionBegin;
   *   field_data[0] = (*x);
   *   field_data[1] = (*y);
   *   field_data[2] = (*z);
   *   MoFEMFunctionReturn(0);
   * };
   * CHKERR m_field.getInterface<FieldBlas>()->setVertexDofs(set_distance,
   * "DISP"); \endcode
   *
   * \note Function works both ways, using it coordinates can be set from
   * field.
   *
   * \param lambda function evaluating field at points
   * \param field_name  is a field name
   * \param verts pointer to vertices if null all vertices in the field are
   * evaluated)
   *
   */
  MoFEMErrorCode setVertexDofs(VertexCoordsFunction lambda,
                               const std::string field_name,
                               Range *verts = nullptr);

  /** \deprecated scale field
   * \ingroup mofem_field_algebra
   *
   * \param val is a set parameter
   * \field_name  is a field name
   *
   */
  MoFEMErrorCode setField(const double val, const EntityType type,
                          const std::string field_name);

  /** \brief set field
   * \ingroup mofem_field_algebra
   *
   * field_y = val
   *
   * \param val
   * \param entity type
   * \param field_name
   *
   */
  MoFEMErrorCode setField(const double val, const EntityType type,
                          const Range &ents, const std::string field_name);

  /** \brief set field
   * \ingroup mofem_field_algebra
   *
   * field_y = val
   *
   * \param val
   * \param field_name
   *
   */
  MoFEMErrorCode setField(const double val, const std::string field_name);

};

} // namespace MoFEM

#endif // __FIELD_BLAS_HPP__
