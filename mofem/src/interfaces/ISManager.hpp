/** \file ISManager.hpp
 * \brief Interface managing IS
 * \ingroup mofem_is_managers
 *
 * Managing problems, build and partitioning.
 *
 */

#ifndef _ISMANAGER_HPP__
#define _ISMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

/**
 * \brief Section manager is used to create indexes and sections
 * \ingroup mofem_is_managers
 *
 * FIXME: ISManager is not properly tested by atom tests.
 *
 */
struct ISManager : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  const MoFEM::Interface &cOre;
  bool dEbug;

  ISManager(const MoFEM::Core &core);

  /**
   * \brief Destructor
   */
  virtual ~ISManager() = default;

  /**
   * \brief Create global selection
   * \ingroup mofem_is_managers
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
  MoFEMErrorCode sectionCreate(const std::string problem_name, PetscSection *s,
                               const RowColData row_col = COL) const;

  /**
   * \brief Create global selection
   * \ingroup mofem_is_managers
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
  SmartPetscObj<PetscSection>
  sectionCreate(const std::string problem_name,
                const RowColData row_col = COL) const;

  /**
   * @brief Create IS for problem
   *
   * @param problem_name
   * @param rc
   * @param is
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode isCreateProblem(const std::string problem_name, RowColData rc,
                                 IS *is) const;

  /**
   * @brief Create IS for problem
   *
   * @param problem_name
   * @param rc
   * @param is
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode isCreateProblem(const std::string problem_name, RowColData rc,
                                 SmartPetscObj<IS> &is) const;

  /**
    * \brief create IS for given order range (collective)
    * \ingroup mofem_is_managers

    * \param problem name
    * \param rc ROW or COL
    * \param min_order
    * \param max_order
    * \retval is out value

    */
  MoFEMErrorCode isCreateProblemOrder(const std::string problem_name,
                                      RowColData rc, int min_order,
                                      int max_order, IS *is) const;

  /**
   * @copydoc MoFEM::ISManager::isCreateProblemOrder
   */
  MoFEMErrorCode isCreateProblemOrder(const std::string problem_name,
                                      RowColData rc, int min_order,
                                      int max_order,
                                      SmartPetscObj<IS> &is) const;

  /**
    * \brief create IS for given problem, field and rank range (collective)
    * \ingroup mofem_is_managers

    * \param problem name
    * \param rc ROW or COL
    * \param field name
    * \param min_coeff_idx
    * \param max_coeff_idx
    * \param ents if not null get dofs only on given entities
    * \retval is out value

    */
  MoFEMErrorCode isCreateProblemFieldAndRank(const std::string problem_name,
                                             RowColData rc,
                                             const std::string field,
                                             int min_coeff_idx,
                                             int max_coeff_idx, IS *is,
                                             Range *ents = nullptr) const;

  /**
   * \copybrief create IS for given problem, field and rank range (collective)
   * \ingroup mofem_is_managers
   *
   * \param problem name
   * \param rc ROW or COL
   * \param field name
   * \param min_coeff_idx
   * \param max_coeff_idx
   * \param ents if not null get dofs only on given entities
   * \retval is out value
   */
  MoFEMErrorCode
  isCreateProblemFieldAndRank(const std::string problem_name, RowColData rc,
                              const std::string field, int min_coeff_idx,
                              int max_coeff_idx, SmartPetscObj<IS> &smart_is,
                              Range *ents = nullptr) const;

  /**
    * \brief create IS for given problem, field and rank range (collective)
    * \ingroup mofem_is_managers

    * \param problem name
    * \param rc ROW or COL
    * \param field name
    * \param min_coeff_idx
    * \param max_coeff_idx
    * \param ents if not null get dofs only on given entities
    * \retval is out value

    */
  MoFEMErrorCode
  isCreateProblemFieldAndRankLocal(const std::string problem_name,
                                   RowColData rc, const std::string field,
                                   int min_coeff_idx, int max_coeff_idx, IS *is,
                                   Range *ents = nullptr) const;

  /**
   * \copybrief create IS for given problem, field and rank range (collective)
   * \ingroup mofem_is_managers
   *
   * \param problem name
   * \param rc ROW or COL
   * \param field name
   * \param min_coeff_idx
   * \param max_coeff_idx
   * \param ents if not null get dofs only on given entities
   * \retval is out value
   */
  MoFEMErrorCode isCreateProblemFieldAndRankLocal(
      const std::string problem_name, RowColData rc, const std::string field,
      int min_coeff_idx, int max_coeff_idx, SmartPetscObj<IS> &smart_is,
      Range *ents = nullptr) const;

  /**
   * @brief create IS for given problem, field and type range (collective)
   * \ingroup mofem_is_managers
   *
   * @param problem
   * @param rc
   * @param field
   * @param low_type
   * @param hi_type
   * @param min_coeff_idx
   * @param max_coeff_idx
   * @param is
   * @param ents
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode isCreateProblemFieldAndEntityType(
      const std::string problem_name, RowColData rc, const std::string field,
      EntityType low_type, EntityType hi_type, int min_coeff_idx,
      int max_coeff_idx, IS *is, Range *ents = nullptr) const;

  /** \brief create IS for give two problems and field
    * \ingroup mofem_is_managers

    Note that indices are ordered in ascending order of local indices in
    problem_y

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
      const std::string x_problem, const std::string x_field_name,
      RowColData x_rc, const std::string y_problem,
      const std::string y_field_name, RowColData y_rc, std::vector<int> &idx,
      std::vector<int> &idy) const;

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
      const std::string x_problem, const std::string x_field_name,
      RowColData x_rc, const std::string y_problem,
      const std::string y_field_name, RowColData y_rc, IS *ix, IS *iy) const;

  /**
   * @brief Create is from one problem to other problem
   * @ingroup mofem_is_managers
   *
   * @param x_problem
   * @param x_rc
   * @param y_problem
   * @param y_rc
   * @param idx
   * @param idy
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode isCreateFromProblemToOtherProblem(
      const std::string x_problem, RowColData x_rc, const std::string y_problem,
      RowColData y_rc, std::vector<int> &idx, std::vector<int> &idy) const;

  /**
   * @brief Create is from one problem to other problem
   * @ingroup mofem_is_managers
   *
   * @param x_problem
   * @param x_rc
   * @param y_problem
   * @param y_rc
   * @param ix
   * @param iy
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode isCreateFromProblemToOtherProblem(const std::string x_problem,
                                                   RowColData x_rc,
                                                   const std::string y_problem,
                                                   RowColData y_rc, IS *ix,
                                                   IS *iy) const;
};

} // namespace MoFEM

/**
 * \defgroup mofem_is_managers Index sets (IS)
 * \brief Construct index sets for MoFEM problems
 *
 * \ingroup mofem
 */

#endif // _ISMANAGER_HPP__
