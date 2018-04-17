/** \file ISManager.hpp
 * \brief Interface managing IS
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

#ifndef _ISMANAGER_HPP__
#define _ISMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

  static const MOFEMuuid IDD_MOFEMISManager = MOFEMuuid( BitIntefaceId(ISMANAGER_INTERFACE) );

  /**
   * \brief Section manager is used to create indexes and sections
   * \ingroup mofem_is_managers
   *
   */
  struct ISManager: public UnknownInterface {

    MoFEMErrorCode query_interface(const MOFEMuuid& uuid, UnknownInterface** iface) const;

    const MoFEM::Interface& cOre;
    bool dEbug;

    ISManager(const MoFEM::Core& core);

    /**
    * \brief Destructor
    */
    ~ISManager();

    /**
     * \brief Create global selection
     *
     * Create section for given problem, such that points are sorted by UId,
     * and number of dofs on point is equal to number of dofs on entity
     *
     * It takes all fields
     *
     * @param  problem_name
     * @param  fields_list
     * @param  s
     * @param  row_col      ROE or COL, default is ROW
     * @return              error code
     */
    MoFEMErrorCode sectionCreate(
      const std::string &problem_name,
      PetscSection *s,
      const RowColData row_col = COL
    ) const;

    /**
      * \brief create IS for given order range (collective)
      * \ingroup mofem_is_managers

      * \param problem name
      * \param rc ROW or COL
      * \param min_order
      * \param max_order
      * \retval is out value

      */
    MoFEMErrorCode isCreateProblemOrder(
      const std::string &problem,RowColData rc,int min_order,int max_order,IS *is
    ) const;

    /**
      * \brief create IS for given problem, field and rank range (collective)
      * \ingroup mofem_is_managers

      * \param problem name
      * \param rc ROW or COL
      * \param field name
      * \param min_coeff_idx
      * \param max_coeff_idx
      * \retval is out value

      */
    MoFEMErrorCode isCreateProblemFieldAndRank(
      const std::string &problem,
      RowColData rc,
      const std::string &field,
      int min_coeff_idx,
      int max_coeff_idx,
      IS *is,
      Range *ents = nullptr
    ) const;

    /** \brief create IS for give two problems and field
      * \ingroup mofem_is_managers

      Note that indices are ordered in ascending order of local indices in problem_y

      \param x_problem name of problem
      \param x_field_name name of field in problem_x
      \param x_rc that is ROW or COL
      \param y_problem name of problem
      \param y_field_name name of field in problem_y
      \param y_rc that is ROW or COL

      \retval idx indexes in problem_x
      \retval idy indexes in problem_y

      */
    MoFEMErrorCode isCreateFromProblemFieldToOtherProblemField(
      const std::string &x_problem,
      const std::string &x_field_name,RowColData x_rc,
      const std::string &y_problem,
      const std::string &y_field_name,RowColData y_rc,
      std::vector<int> &idx,std::vector<int> &idy
    ) const;

    /** \brief create IS for give two problems and field
      * \ingroup mofem_is_managers

      Indices are sorted by global PETSc index in problem_x.

      \param x_problem name of problem
      \param x_field_name name of field in problem_x
      \param x_rc that is ROW or COL
      \param y_problem name of problem
      \param y_field_name name of field in problem_y
      \param y_rc that is ROW or COL

      \retval ix IS indexes in problem_x (can be PETSC_NULL)
      \retval iy IS indexes in problem_y

      */
    MoFEMErrorCode isCreateFromProblemFieldToOtherProblemField(
      const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
      const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
      IS *ix,IS *iy
    ) const;

    MoFEMErrorCode isCreateFromProblemToOtherProblem(
      const std::string &x_problem,RowColData x_rc,
      const std::string &y_problem,RowColData y_rc,
      std::vector<int> &idx,std::vector<int> &idy
    ) const;

    MoFEMErrorCode isCreateFromProblemToOtherProblem(
      const std::string &x_problem,RowColData x_rc,
      const std::string &y_problem,RowColData y_rc,
      IS *ix,IS *iy
    ) const;

  };

}

/***************************************************************************//**
 * \defgroup mofem_is_managers Index sets (IS)
 * \brief Construct index sets for MoFEM problems
 *
 * \ingroup mofem
 ******************************************************************************/


#endif // _ISMANAGER_HPP__
