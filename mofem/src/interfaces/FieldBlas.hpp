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

  static const MOFEMuuid IDD_MOFEMFieldBlas = MOFEMuuid( BitIntefaceId(FIELDBLAS_INTERFACE) );

  /**
   * \brief Section manager is used to create sections
   * \mofem_is_managers
   *
   */
  struct FieldBlas: public UnknownInterface {

    MoFEMErrorCode query_interface(const MOFEMuuid& uuid, UnknownInterface** iface) const;

    const MoFEM::Interface& cOre;
    bool dEbug;

    FieldBlas(const MoFEM::Core& core);

    /**
    * \brief Destructor
    */
    ~FieldBlas();


    /** \brief axpy fields
      * \ingroup mofem_field_algebra
      * \todo should be moved to independent interface, i.e. FieldAlgebra
      *
      * field_y = field_y + alpha*field_x
      *
      * \param alpha
      * \param field_name_x name of field_x
      * \param field_name_y name of field_y
      * \param error_if_missing throw error if entity/dof exist in field_x but not on field_y
      * \param create_if_missing creat dof in field_y from fiedl_x if it is not database
      *
      */
    MoFEMErrorCode fieldAxpy(
      const double alpha,const std::string& field_name_x,const std::string& field_name_y,
      bool error_if_missing = false,bool creat_if_missing = false
    );


    /** \brief scale field
      * \ingroup mofem_field_algebra
      * \todo should be moved to independent interface, i.e. FieldAlgebra
      *
      * \param alpha is a scaling factor
      * \field_name  is a field name
      *
      */
    MoFEMErrorCode setField(const double val,const EntityType type,const std::string& field_name);

    /** \brief set field
      * \ingroup mofem_field_algebra
      * \todo should be moved to independent interface, i.e. FieldAlgebra
      *
      * field_y = val
      *
      * \param val
      * \param entity type
      * \param field_name
      *
      */
    MoFEMErrorCode setField(
      const double val,const EntityType type,const Range &ents,const std::string& field_name
    );

    /** \brief set field
      * \ingroup mofem_field_algebra
      * \todo should be moved to independent interface, i.e. FieldAlgebra
      *
      * field_y = val
      *
      * \param val
      * \param entity type
      * \param on enties
      * \param field_name
      *
      */
    MoFEMErrorCode fieldScale(const double alpha,const std::string& field_name);

  };

}

#endif // __FIELD_BLAS_HPP__

/***************************************************************************//**
 * \defgroup mofem_field_algebra Field Basic Algebra
 * \brief Basic algebraic operation on fields
 *
 * \ingroup mofem
 ******************************************************************************/
