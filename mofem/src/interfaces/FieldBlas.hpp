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

  /**
   * \brief Destructor
   */
  ~FieldBlas();

  typedef boost::function<MoFEMErrorCode(double &, const double)>
      TwoFieldFunction;

  /** \brief filed lambda
   * \ingroup mofem_field_algebra
   *
   * Do calculation on two fields and save result to field fy
   *
   * \code
   struct Axpy {
    const double aLpha;
    Axpy(const double alpha) : aLpha(alpha) {}
    inline MoFEMErrorCode operator(double &fy, double fx) {
      MoFEMFunctionBeginHot;
      fy += Alpha * fx;
      MoFEMFunctionReturnHot(0);
    }
   };
   CHKERR m_fiel.getInterface<FieldBlas>()->fieldLambda(Axpy(aLpha),
   field_name_x, field_name_y);
   * \endcode
   *
   * \param function f(double &x, double)
   * \param field_name_x name of field_x
   * \param field_name_y name of field_y
   * \param error_if_missing throw error if entity/dof exist in field_x but not
   * on field_y \param create_if_missing creat dof in field_y from field_x if it
   * is not database
   *
   */
  MoFEMErrorCode fieldLambda(TwoFieldFunction lambda,
                             const std::string field_name_x,
                             const std::string field_name_y,
                             bool error_if_missing = false,
                             bool creat_if_missing = false);

  /** \brief axpy fields
   * \ingroup mofem_field_algebra
   *
   * field_y = field_y + alpha*field_x
   *
   * \param alpha
   * \param field_name_x name of field_x
   * \param field_name_y name of field_y
   * \param error_if_missing throw error if entity/dof exist in field_x but not
   * on field_y \param create_if_missing creat dof in field_y from fiedl_x if it
   * is not database
   *
   */
  MoFEMErrorCode fieldAxpy(const double alpha, const std::string field_name_x,
                           const std::string field_name_y,
                           bool error_if_missing = false,
                           bool creat_if_missing = false);

  /** \brief copy and scale fields
   * \ingroup mofem_field_algebra
   *
   * field_y = alpha*field_x
   *
   * \param alpha
   * \param field_name_x name of field_x
   * \param field_name_y name of field_y
   * \param error_if_missing throw error if entity/dof exist in field_x but not
   * on field_y \param create_if_missing creat dof in field_y from fiedl_x if it
   * is not database
   *
   */
  MoFEMErrorCode fieldCopy(const double alpha, const std::string field_name_x,
                           const std::string field_name_y,
                           bool error_if_missing = false,
                           bool creat_if_missing = false);

  typedef boost::function<MoFEMErrorCode(VectorAdaptor &&field_data,
                                         double *xcoord, double *ycoord,
                                         double *zcoord)>
      VertexCoordsFunction;

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
   * \note Function works both ways, using it coordinates can be set from field.
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

  /** \brief scale field
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

  /** \brief set field
   * \ingroup mofem_field_algebra
   *
   * field_y = val
   *
   * \param val
   * \param entity type
   * \param on enties
   * \param field_name
   *
   */
  MoFEMErrorCode fieldScale(const double alpha, const std::string field_name);
};

} // namespace MoFEM

#endif // __FIELD_BLAS_HPP__

/**
 * \defgroup mofem_field_algebra Field Basic Algebra
 * \brief Basic algebraic operation on fields
 *
 * \ingroup mofem
 */
