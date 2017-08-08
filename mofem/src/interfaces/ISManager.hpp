/** \file ISManager.hpp
 * \brief Interface managing sections
 * \ingroup mofem_section_manager
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

#ifndef __SECTIONMANAGER_HPP__
#define __SECTIONMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

  static const MOFEMuuid IDD_MOFEMISManager = MOFEMuuid( BitIntefaceId(ISMANAGER_INTERFACE) );

  /**
   * \brief Section manager is used to create sections
   * \mofem_section_manager
   *
   */
  struct ISManager: public UnknownInterface {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

    const MoFEM::Interface& cOre;
    bool dEbug;

    ISManager(const MoFEM::Core& core);

    /**
    * \brief Destructor
    */
    ~ISManager();

    PetscErrorCode sectionCreateFromFieldsList(
      const std::string &problem_name,
      const std::vector<std::string> &fields_list,
      PetscSection *s,
      const RowColData row_col = COL
    ) const;

    /**
      * \brief create IS for given order range (collective)
      * \ingroup mofem_section_manager

      * \param problem name
      * \param rc ROW or COL
      * \param min_order
      * \param max_order
      * \retval is out value

      */
    PetscErrorCode isCreateProblemOrder(
      const std::string &problem,RowColData rc,int min_order,int max_order,IS *is
    ) const;

    /**
      * \brief create IS for given problem, field and rank range (collective)
      * \ingroup mofem_section_manager

      * \param problem name
      * \param rc ROW or COL
      * \param field name
      * \param min_coeff_idx
      * \param max_coeff_idx
      * \retval is out value

      */
    PetscErrorCode isCreateProblemFieldAndRank(
      const std::string &problem,
      RowColData rc,
      const std::string &field,
      int min_coeff_idx,
      int max_coeff_idx,
      IS *is
    ) const;

    /** \brief create IS for give two problems and field
      * \ingroup mofem_section_manager

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
    PetscErrorCode isCreateFromProblemFieldToOtherProblemField(
      const std::string &x_problem,
      const std::string &x_field_name,RowColData x_rc,
      const std::string &y_problem,
      const std::string &y_field_name,RowColData y_rc,
      std::vector<int> &idx,std::vector<int> &idy
    ) const;

    /** \brief create IS for give two problems and field
      * \ingroup mofem_section_manager

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
    PetscErrorCode isCreateFromProblemFieldToOtherProblemField(
      const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
      const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
      IS *ix,IS *iy
    ) const;

    PetscErrorCode isCreateFromProblemToOtherProblem(
      const std::string &x_problem,RowColData x_rc,
      const std::string &y_problem,RowColData y_rc,
      std::vector<int> &idx,std::vector<int> &idy
    ) const;

    PetscErrorCode isCreateFromProblemToOtherProblem(
      const std::string &x_problem,RowColData x_rc,
      const std::string &y_problem,RowColData y_rc,
      IS *ix,IS *iy
    ) const;

  };

}

#endif // __SECTIONMANAGER_HPP__
