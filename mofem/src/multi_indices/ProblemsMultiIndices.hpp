/** \file ProblemsMultiIndices.hpp
 * \brief Mylti-index containers, data structures for problems and other low-level functions
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __PROBLEMSMULTIINDICES_HPP__
#define __PROBLEMSMULTIINDICES_HPP__

namespace MoFEM {

/** \brief keeps data about problem
  * \ingroup problems_multi_indices
  *
  * \todo Implement layouts for DOFs
  *
  */
struct MoFEMProblem {
  EntityHandle meshset;
  BitProblemId* tag_id_data;
  const void* tag_name_data;
  int tag_name_size;
  DofIdx* tag_nbdof_data_row;
  DofIdx* tag_nbdof_data_col;
  DofIdx* tag_local_nbdof_data_row;
  DofIdx* tag_local_nbdof_data_col;
  DofIdx* tag_ghost_nbdof_data_row;
  DofIdx* tag_ghost_nbdof_data_col;
  BitFEId* tag_BitFEId_data;
  BitRefLevel* tag_BitRefLevel;
  BitRefLevel* tag_BitRefLevel_DofMask;
  boost::shared_ptr<NumeredDofEntity_multiIndex> numered_dofs_rows;
  boost::shared_ptr<NumeredDofEntity_multiIndex> numered_dofs_cols;
  NumeredEntFiniteElement_multiIndex numeredFiniteElements;

  /**
    * use with loops to iterate problem fes
    * \ingroup problems_multi_indices
    *
    * for(_IT_NUMEREDFEMOFEMENTITY_FOR_LOOP_(MOFEMPROBLEM,NAME,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDFEMOFEMENTITY_BY_NAME_FOR_LOOP_(MOFEMPROBLEM,NAME,IT) \
    NumeredEntFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_fes_begin(NAME); \
    IT!=MOFEMPROBLEM->get_numered_fes_end(NAME); IT++

  NumeredEntFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator get_numered_fes_begin(std::string fe_name) const {
    return numeredFiniteElements.get<FiniteElement_name_mi_tag>().lower_bound(fe_name);
  }

  NumeredEntFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator get_numered_fes_end(std::string fe_name) const {
    return numeredFiniteElements.get<FiniteElement_name_mi_tag>().upper_bound(fe_name);
  }

  /**
    * \brief use with loops to iterate problem fes
    * \ingroup problems_multi_indices
    *
    * \code
    * for(_IT_NUMEREDFEMOFEMENTITY_BY_NAME_AND_PART_FOR_LOOP_(MOFEMPROBLEM,NAME,PART,IT)) {
    *   ...
    * }
    * \endcode
    *
    */
  #define _IT_NUMEREDFEMOFEMENTITY_BY_NAME_AND_PART_FOR_LOOP_(MOFEMPROBLEM,NAME,PART,IT) \
    NumeredEntFiniteElement_multiIndex::index<Composite_Name_And_Part_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_fes_begin(NAME,PART); \
    IT!=MOFEMPROBLEM->get_numered_fes_end(NAME,PART); IT++

  NumeredEntFiniteElement_multiIndex::index<Composite_Name_And_Part_mi_tag>::type::iterator get_numered_fes_begin(std::string fe_name,int part) const {
    return numeredFiniteElements.get<Composite_Name_And_Part_mi_tag>().lower_bound(boost::make_tuple(fe_name,part));
  }

  NumeredEntFiniteElement_multiIndex::index<Composite_Name_And_Part_mi_tag>::type::iterator get_numered_fes_end(std::string fe_name,int part) const {
    return numeredFiniteElements.get<Composite_Name_And_Part_mi_tag>().upper_bound(boost::make_tuple(fe_name,part));
  }

  /**
    * \brief use with loops to iterate row dofs
    * \ingroup problems_multi_indices
    *
    * \code
    * for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    * \endcode
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(MOFEMPROBLEM,IT) \
    NumeredDofEntity_multiIndex::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_begin(); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_end(); IT++

  /**
    * use with loops to iterate col dofs
    * \ingroup problems_multi_indices
    *
    * \code
    * for(_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    * \endcode
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(MOFEMPROBLEM,IT) \
    NumeredDofEntity_multiIndex::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_begin(); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_end(); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator get_numered_dofs_rows_begin() const { return numered_dofs_rows->begin(); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator get_numered_dofs_rows_end() const { return numered_dofs_rows->end(); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator get_numered_dofs_cols_begin() const { return numered_dofs_cols->begin(); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator get_numered_dofs_cols_end() const { return numered_dofs_cols->end(); }

  /**
    * \brief get iterator of dof in row by uid
    * \ingroup problems_multi_indices
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_UID_(MOFEMPROBLEM,UID,IT) \
    NumeredDofEntity_multiIndex::index<Unique_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_row_dof_by_uid(UID);

  /**
    * \brief get iterator of dof in col by uid
    * \ingroup problems_multi_indices
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_UID_(MOFEMPROBLEM,UID,IT) \
    NumeredDofEntity_multiIndex::index<Unique_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_col_dof_by_uid(UID);

  /// get iterator of dof in row by uid (instead you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_UID_FOR_LOOP_)
  NumeredDofEntity_multiIndex::index<Unique_mi_tag>::type::iterator get_row_dof_by_uid(GlobalUId uid) const { return numered_dofs_rows->get<Unique_mi_tag>().find(uid); };

  /// get iterator of dof in column by uid (instead you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_UID_FOR_LOOP_)
  NumeredDofEntity_multiIndex::index<Unique_mi_tag>::type::iterator get_col_dof_by_uid(GlobalUId uid) const { return numered_dofs_cols->get<Unique_mi_tag>().find(uid); };

  /**
    * \brief use with loops to iterate row dofs
    * \ingroup problems_multi_indices
    *
    * \code
    * for(_IT_NUMEREDDOFMOFEMENTITY_BY_LOCIDX_ROW_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    * \endcode
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_LOCIDX_FOR_LOOP_(MOFEMPROBLEM,IT) \
    NumeredDofEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_by_locidx_begin(0); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_by_locidx_end(MOFEMPROBLEM->get_nb_local_dofs_row()-1); IT++

  /**
    * \brief use with loops to iterate col dofs
    *
    * \code
    * for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_LOCIDX_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    * \endcode
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_LOCIDX_FOR_LOOP_(MOFEMPROBLEM,IT) \
    NumeredDofEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_by_locidx_begin(0); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_by_locidx_end(MOFEMPROBLEM->get_nb_local_dofs_row()-1); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator get_numered_dofs_rows_by_locidx_begin(const DofIdx locidx) const
    { return numered_dofs_rows->get<PetscLocalIdx_mi_tag>().lower_bound(locidx); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator get_numered_dofs_rows_by_locidx_end(const DofIdx locidx) const
    { return numered_dofs_rows->get<PetscLocalIdx_mi_tag>().upper_bound(locidx); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator get_numered_dofs_cols_by_locidx_begin(const DofIdx locidx) const
    { return numered_dofs_cols->get<PetscLocalIdx_mi_tag>().lower_bound(locidx); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator get_numered_dofs_cols_by_locidx_end(const DofIdx locidx) const
    { return numered_dofs_cols->get<PetscLocalIdx_mi_tag>().upper_bound(locidx); }

  /**
    * \brief use with loops to iterate row dofs
    * \ingroup problems_multi_indices
    *
    * \code
    * for(_IT_NUMEREDDOFMOFEMENTITY_BY_ENT_ROW_FOR_LOOP_(MOFEMPROBLEM,ENT,IT)) {
    *   ...
    * }
    * \endcode
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(MOFEMPROBLEM,ENT,IT) \
    NumeredDofEntity_multiIndex::index<Ent_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_by_ent_begin(ENT); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_by_ent_end(ENT); IT++

  /**
    * \brief use with loops to iterate col dofs
    * \ingroup problems_multi_indices
    *
    * \code
    * for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_(MOFEMPROBLEM,ENT,IT)) {
    *   ...
    * }
    * \endcode
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_(MOFEMPROBLEM,ENT,IT) \
    NumeredDofEntity_multiIndex::index<Ent_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_by_ent_begin(ENT); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_by_ent_end(ENT); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<Ent_mi_tag>::type::iterator get_numered_dofs_rows_by_ent_begin(const EntityHandle ent) const
    { return numered_dofs_rows->get<Ent_mi_tag>().lower_bound(ent); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<Ent_mi_tag>::type::iterator get_numered_dofs_rows_by_ent_end(const EntityHandle ent) const
    { return numered_dofs_rows->get<Ent_mi_tag>().upper_bound(ent); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<Ent_mi_tag>::type::iterator get_numered_dofs_cols_by_ent_begin(const EntityHandle ent) const
    { return numered_dofs_cols->get<Ent_mi_tag>().lower_bound(ent); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<Ent_mi_tag>::type::iterator get_numered_dofs_cols_by_ent_end(const EntityHandle ent) const
    { return numered_dofs_cols->get<Ent_mi_tag>().upper_bound(ent); }

  /**
    * use with loops to iterate row dofs
    * \ingroup problems_multi_indices
    *
    * \code
    * for(_IT_NUMEREDDOFMOFEMENTITY_BY_NAME_ROW_FOR_LOOP_(MOFEMPROBLEM,NAME,IT)) {
    *   ...
    * }
    * \endcode
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(MOFEMPROBLEM,NAME,IT) \
    NumeredDofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_begin(NAME); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_end(NAME); IT++

  /**
    * \brief use with loops to iterate col dofs
    * \ingroup problems_multi_indices
    *
    * \code
    * for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_FOR_LOOP_(MOFEMPROBLEM,NAME,IT)) {
    *   ...
    * }
    * \endcode
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_FOR_LOOP_(MOFEMPROBLEM,NAME,IT) \
    NumeredDofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_begin(NAME); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_end(NAME); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_numered_dofs_rows_begin(const std::string& name) const
    { return numered_dofs_rows->get<FieldName_mi_tag>().lower_bound(name); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_numered_dofs_rows_end(const std::string& name) const
    { return numered_dofs_rows->get<FieldName_mi_tag>().upper_bound(name); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_numered_dofs_cols_begin(const std::string& name) const
    { return numered_dofs_cols->get<FieldName_mi_tag>().lower_bound(name); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_numered_dofs_cols_end(const std::string& name) const
    { return numered_dofs_cols->get<FieldName_mi_tag>().upper_bound(name); }

  /**
    * \brief use with loops to iterate row dofs
    * \ingroup problems_multi_indices
    *
    * \code
    * for(_IT_NUMEREDDOFMOFEMENTITY_BY_NAME_ENT_PART_ROW_FOR_LOOP_(MOFEMPROBLEM,NAME,ENT,PART,IT)) {
    *   ...
    * }
    * \endcode
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(MOFEMPROBLEM,NAME,ENT,PART,IT) \
    NumeredDofEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_begin(NAME,ENT,PART); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_end(NAME,ENT,PART); IT++

  /**
    * use with loops to iterate col dofs
    * \ingroup problems_multi_indices
    *
    * \code
    * for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_ENT_PART_FOR_LOOP_(MOFEMPROBLEM,NAME,ENT,PART,IT)) {
    *   ...
    * }
    * \endcode
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_ENT_PART_FOR_LOOP_(MOFEMPROBLEM,NAME,ENT,PART,IT) \
    NumeredDofEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_begin(NAME,ENT,PART); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_end(NAME,ENT,PART); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator get_numered_dofs_rows_begin(const std::string& name,const EntityHandle ent,const int part) const {
    return numered_dofs_rows->get<Composite_Name_Ent_And_Part_mi_tag>().lower_bound(boost::make_tuple(name,ent,part));
  }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator get_numered_dofs_rows_end(const std::string& name,const EntityHandle ent,const int part) const {
    return numered_dofs_rows->get<Composite_Name_Ent_And_Part_mi_tag>().upper_bound(boost::make_tuple(name,ent,part));
  }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_ENT_PART_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator get_numered_dofs_cols_begin(const std::string& name,const EntityHandle ent,const int part) const {
    return numered_dofs_cols->get<Composite_Name_Ent_And_Part_mi_tag>().lower_bound(boost::make_tuple(name,ent,part));
  }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_ENT_PART_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator get_numered_dofs_cols_end(const std::string& name,const EntityHandle ent,const int part) const {
    return numered_dofs_cols->get<Composite_Name_Ent_And_Part_mi_tag>().upper_bound(boost::make_tuple(name,ent,part));
  }

  /**
   * \berief Loop over problem DOFs in row by part
   * @param  MOFEMPROBLEM problem pointer
   * @param  PART         partition number
   * @param  IT           iterator
   * @return              error code
   * \ingroup problems_multi_indices
   *
   * \code
   * for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_PART_FOR_LOOP_(PART,IT)) {
   *   ...
   * }
   * \endcode
   */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_PART_FOR_LOOP_(MOFEMPROBLEM,PART,IT) \
  NumeredDofEntity_multiIndex::index<Part_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_begin(PART); \
  IT!=MOFEMPROBLEM->get_numered_dofs_rows_end(PART); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_PART_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<Part_mi_tag>::type::iterator get_numered_dofs_rows_begin(const int part) const {
    return numered_dofs_rows->get<Part_mi_tag>().lower_bound(part);
  }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_PART_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<Part_mi_tag>::type::iterator get_numered_dofs_rows_end(const int part) const {
    return numered_dofs_rows->get<Part_mi_tag>().upper_bound(part);
  }

  /**
   * \berief Loop over problem DOFs in col by part
   * @param  MOFEMPROBLEM problem pointer
   * @param  PART         partition number
   * @param  IT           iterator
   * @return              error code
   * \ingroup problems_multi_indices
   *
   * \code
   * for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_PART_FOR_LOOP_(PART,IT)) {
   *   ...
   * }
   * \endcode
   *
   */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_PART_FOR_LOOP_(MOFEMPROBLEM,PART,IT) \
  NumeredDofEntity_multiIndex::index<Part_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_begin(PART); \
  IT!=MOFEMPROBLEM->get_numered_dofs_cols_end(PART); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_PART_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<Part_mi_tag>::type::iterator get_numered_dofs_cols_begin(const int part) const {
    return numered_dofs_rows->get<Part_mi_tag>().lower_bound(part);
  }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_PART_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::index<Part_mi_tag>::type::iterator get_numered_dofs_cols_end(const int part) const {
    return numered_dofs_rows->get<Part_mi_tag>().upper_bound(part);
  }

  MoFEMProblem(Interface &moab,const EntityHandle _meshset);
  inline BitProblemId getId() const { return *((BitProblemId*)tag_id_data); }

  // /** \deprecated Use getId() instead
  // */
  // DEPRECATED inline BitProblemId get_id() const { return getId(); }

  inline std::string getName() const { return std::string((char *)tag_name_data,tag_name_size); }

  // /** \deprecated Use getName() instead
  // */
  // DEPRECATED inline std::string get_name() const { return getName(); }


  inline DofIdx getNbDofsRow() const { return *((DofIdx*)tag_nbdof_data_row); }
  inline DofIdx getNbDofsCol() const { return *((DofIdx*)tag_nbdof_data_col); }
  inline DofIdx get_nb_local_dofs_row() const { return *((DofIdx*)tag_local_nbdof_data_row); }
  inline DofIdx get_nb_local_dofs_col() const { return *((DofIdx*)tag_local_nbdof_data_col); }
  inline DofIdx get_nb_ghost_dofs_row() const { return *((DofIdx*)tag_ghost_nbdof_data_row); }
  inline DofIdx get_nb_ghost_dofs_col() const { return *((DofIdx*)tag_ghost_nbdof_data_col); }

  inline BitRefLevel getBitRefLevel() const { return *tag_BitRefLevel; }

  // /** \deprecated Use getBitRefLevel() instead
  // */
  // DEPRECATED inline BitRefLevel get_BitRefLevel() const { return getBitRefLevel(); }


  inline BitRefLevel get_DofMask_BitRefLevel() const { return *tag_BitRefLevel_DofMask; }
  PetscErrorCode getRowDofsByPetscGlobalDofIdx(DofIdx idx,const NumeredDofEntity **dof_ptr) const;
  PetscErrorCode getColDofsByPetscGlobalDofIdx(DofIdx idx,const NumeredDofEntity **dof_ptr) const;
  BitFEId get_BitFEId() const;
  friend std::ostream& operator<<(std::ostream& os,const MoFEMProblem& e);

  /**
   * \brief Get number of finite elements by name on processors
   *
   * Size retuned IS is equal to size of processors.
   *
   * What PetscLayout see, <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/IS/index.html>
   *
   * Example of usage for layout
   * \code

   PetscInt rstart, rend;
   ierr = PetscLayoutGetRange(layout, &rstart, &rend); CHKERRQ(ierr);
   int global_size;
   ierr = PetscLayoutGetSize(layout,&global_size); CHKERRQ(ierr);

   \endcode
   *
   * @param  comm Communicator
   * @param  name Finite element name
   * @param  layout Get number of elements on each processor
   * @return      error code
   */
  PetscErrorCode getNumberOfElementsByNameAndPart(MPI_Comm comm,const std::string name,PetscLayout *layout) const;

  /**
   * \brief Get number of finite elements on processors
   *
   * Size retuned IS is equal to size of processors.
   *
   * What PetscLayout see, <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/IS/index.html>
   *
   * Example of usage for layout
   * \code

   PetscInt rstart, rend;
   ierr = PetscLayoutGetRange(layout, &rstart, &rend); CHKERRQ(ierr);
   int global_size;
   ierr = PetscLayoutGetSize(layout,&global_size); CHKERRQ(ierr);

   \endcode
   *
   * @param  comm Communicator
   * @param  layout Get number of elements on each processor
   * @return      error code
   */
  PetscErrorCode getNumberOfElementsByPart(MPI_Comm comm,PetscLayout *layout) const;

};

/**
 * @relates multi_index_container
 * \brief MultiIndex for entities for MoFEMProblem
 * \ingroup fe_multi_indices
 */
typedef multi_index_container<
  MoFEMProblem,
  indexed_by<
    ordered_unique<
      tag<Meshset_mi_tag>, member<MoFEMProblem,EntityHandle,&MoFEMProblem::meshset> >,
    hashed_unique<
      tag<BitProblemId_mi_tag>, const_mem_fun<MoFEMProblem,BitProblemId,&MoFEMProblem::getId>, HashBit<BitProblemId>, EqBit<BitProblemId> >,
    hashed_unique<
      tag<Problem_mi_tag>, const_mem_fun<MoFEMProblem,std::string,&MoFEMProblem::getName> >
  > > MoFEMProblem_multiIndex;

/** \brief add ref level to problem
  * \ingroup problems_multi_indices
  */
struct ProblemChangeRefLevelBitAdd {
  BitRefLevel bit;
  ProblemChangeRefLevelBitAdd(const BitRefLevel _bit): bit(_bit) {};
  void operator()(MoFEMProblem &p) { *(p.tag_BitRefLevel) |= bit; };
};

/** \brief set ref level to problem
  * \ingroup problems_multi_indices
  */
struct ProblemChangeRefLevelBitSet {
  BitRefLevel bit;
  ProblemChangeRefLevelBitSet(const BitRefLevel _bit): bit(_bit) {};
  void operator()(MoFEMProblem &p) { *(p.tag_BitRefLevel) = bit; };
};

/** \brief set prof dof bit ref mask
  * \ingroup problems_multi_indices
  */
struct ProblemChangeRefLevelBitDofMaskSet {
  BitRefLevel bit;
  ProblemChangeRefLevelBitDofMaskSet(const BitRefLevel _bit): bit(_bit) {};
  void operator()(MoFEMProblem &p) { *(p.tag_BitRefLevel_DofMask) = bit; };
};

/** \brief add finite element to problem
  * \ingroup problems_multi_indices
  */
struct ProblemFiniteElementChangeBitAdd {
  BitFEId f_id;
  ProblemFiniteElementChangeBitAdd(const BitFEId _f_id): f_id(_f_id) {};
  void operator()(MoFEMProblem &p);
};

/** \brief remove finite element from problem
  * \ingroup problems_multi_indices
  */
struct ProblemFiniteElementChangeBitUnSet {
  BitFEId f_id;
  ProblemFiniteElementChangeBitUnSet(const BitFEId _f_id): f_id(_f_id) {};
  void operator()(MoFEMProblem &p);
};

/** \brief increase nb. dof in row
  * \ingroup problems_multi_indices
  */
struct ProblemAddRowDof {
  const boost::shared_ptr<DofEntity> dof_ptr;
  ProblemAddRowDof(const boost::shared_ptr<DofEntity> _dof_ptr);
  std::pair<NumeredDofEntity_multiIndex::iterator,bool> p;
  void operator()(MoFEMProblem &e);
};

/** \brief increase nb. dof in col
  * \ingroup problems_multi_indices
  */
struct ProblemAddColDof {
  const boost::shared_ptr<DofEntity> dof_ptr;
  ProblemAddColDof(const boost::shared_ptr<DofEntity> _dof_ptr);
  std::pair<NumeredDofEntity_multiIndex::iterator,bool> p;
  void operator()(MoFEMProblem &e);
};

/** \brief zero nb. of dofs in row
  * \ingroup problems_multi_indices
  */
struct ProblemZeroNbRowsChange {
  void operator()(MoFEMProblem &e);
};

/** \brief zero nb. of dofs in col
  * \ingroup problems_multi_indices
  */
struct ProblemZeroNbColsChange {
  void operator()(MoFEMProblem &e);
};

/** \brief clear problem finite elements
  * \ingroup problems_multi_indices
  */
struct ProblemClearNumeredFiniteElementsChange {
  void operator()(MoFEMProblem &e);
};

/** \brief number dofs in row
  * \ingroup problems_multi_indices
  */
struct ProblemRowNumberChange {
  ProblemRowNumberChange() {};
  void operator()(MoFEMProblem &e);
};

/** \brief number dofs in col
  * \ingroup problems_multi_indices
  */
struct ProblemColNumberChange {
  ProblemColNumberChange() {};
  void operator()(MoFEMProblem &e);
};

}

#endif //__PROBLEMSMULTIINDICES_HPP__

/***************************************************************************//**
 * \defgroup problems_multi_indices Problems structures and multi-indices
 * \ingroup mofem
 ******************************************************************************/
