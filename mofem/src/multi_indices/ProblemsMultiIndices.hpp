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

struct MoFEMProblem;

/**
 * Data structure created when composite problem is created
 */
struct ComposedProblemsData {

  std::vector<const MoFEMProblem*> rowProblemsAdd;
  std::vector<const MoFEMProblem*> colProblemsAdd;

  std::vector<IS> rowIs;
  std::vector<IS> colIs;

  inline PetscErrorCode getRowIs(IS *is,int pp) {
    PetscFunctionBegin;
    PetscObjectReference((PetscObject)rowIs[pp]);
    if(pp<=rowIs.size()) {
      SETERRQ1(
        PETSC_COMM_WORLD,MOFEM_INVALID_DATA,
        "Exceed size of array pp<%d",rowIs.size()
      );
    }
    *is = rowIs[pp];
    PetscFunctionReturn(0);
  }

  inline PetscErrorCode getColIs(IS *is,int pp) {
    PetscFunctionBegin;
    PetscObjectReference((PetscObject)colIs[pp]);
    if(pp<=colIs.size()) {
      SETERRQ1(
        PETSC_COMM_WORLD,MOFEM_INVALID_DATA,
        "Exceed size of array pp<%d",colIs.size()
      );
    }
    *is = colIs[pp];
    PetscFunctionReturn(0);
  }

  virtual ~ComposedProblemsData() {
    for(int ii = 0;ii!=rowIs.size();ii++) {
      ISDestroy(&rowIs[ii]);
    }
    for(int jj = 0;jj!=colIs.size();jj++) {
      ISDestroy(&colIs[jj]);
    }
  }

};


/** \brief keeps basic data about problem
  * \ingroup problems_multi_indices
  *
  * This is low level structure with information about problem, what elements
  * compose problem and what dofs are on rows and columns.
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
  BitRefLevel* tag_MaskBitRefLevel;

  mutable boost::shared_ptr<NumeredDofEntity_multiIndex> numered_dofs_rows; // FIXME name convention
  mutable boost::shared_ptr<NumeredDofEntity_multiIndex> numered_dofs_cols; // FIXME name convention
  mutable NumeredEntFiniteElement_multiIndex numeredFiniteElements;


  /**
   * \brief Subproblem problem data
   */
  struct SubProblemData {

    IS rowIs; ///< indices of main problem of which sub problem is this
    IS colIs; ///< indices of main problem of which sub problem is this
    AO rowMap; ///< mapping form main problem indices to sub-problem indices
    AO colMap;

    /**
     * get row Is for sub problem
     * @param  is create is
     * @return    error code
     */
    inline PetscErrorCode getRowIs(IS *is) {
      PetscFunctionBegin;
      PetscObjectReference((PetscObject)rowIs);
      *is = rowIs;
      PetscFunctionReturn(0);
    }

    /**
     * get col Is for sub problem
     * @param  is create is
     * @return    error code
     */
    inline PetscErrorCode getColIs(IS *is) {
      PetscFunctionBegin;
      PetscObjectReference((PetscObject)colIs);
      *is = colIs;
      PetscFunctionReturn(0);
    };

    /**
     * get row AO mapping for sub problem
     * @param  ao get mapping
     * @return    error code
     */
    inline PetscErrorCode getRowMap(AO *ao) {
      PetscFunctionBegin;
      PetscObjectReference((PetscObject)rowMap);
      *ao = rowMap;
      PetscFunctionReturn(0);
    }

    /**
     * get col AO mapping for sub problem
     * @param  ao get mapping
     * @return    error code
     */
    inline PetscErrorCode getColMap(AO *ao) {
      PetscFunctionBegin;
      PetscObjectReference((PetscObject)colMap);
      *ao = colMap;
      PetscFunctionReturn(0);
    }

    SubProblemData():
    rowIs(PETSC_NULL),
    colIs(PETSC_NULL),
    rowMap(PETSC_NULL),
    colMap(PETSC_NULL) {}
    virtual ~SubProblemData() {
      // int flg;
      // MPI_Finalized(&flg);
      // if(!flg) {
        ISDestroy(&rowIs);
        ISDestroy(&colIs);
        AODestroy(&rowMap);
        AODestroy(&colMap);
      // }
    }
  };

  /**
   * Pointer to data structure. This pointer has allocated data only for
   * sub problems.
   */
  mutable boost::shared_ptr<SubProblemData> subProblemData;

  /**
   * \brief Get main problem of sub-problem is
   * @return    sub problem data structure
   */
  inline boost::shared_ptr<SubProblemData> getSubData() const {
    return subProblemData;
  }

  /**
   * Pointer to data structure from which this problem is composed
   */
  mutable boost::shared_ptr<ComposedProblemsData> composedProblemsData;

  /**
   * \brief Het composed problems data structure
   */
  inline boost::shared_ptr<ComposedProblemsData> getComposedProblemsData() const {
    return composedProblemsData;
  }

  /**
   * \brief get dof from problem
   *
   * Note that \e ent_dof_idx is not coefficient number, is local number of dof on
   * the entity. The coefficient number and local index of dof or entity are
   * the same on vertices and H1 approximation.
   *
   * @param  name       field name
   * @param  ent        entity handle
   * @param  ent_dof_idx index of dof on entity
   * @param  row_or_col ROW or COL
   * @param  dof_ptr    shared pointer to dof if found
   * @return            error code
   */
  PetscErrorCode getDofByNameEntAndEntDofIdx(
    const string name,
    const EntityHandle ent,
    const int ent_dof_idx,
    const RowColData row_or_col,
    boost::shared_ptr<NumeredDofEntity> &dof_ptr
  ) const;

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
    NumeredEntFiniteElementbyName::iterator IT = MOFEMPROBLEM->getNumeredFEsBegin(NAME); \
    IT!=MOFEMPROBLEM->getNumeredFEsEnd(NAME); IT++

  NumeredEntFiniteElementbyName::iterator getNumeredFEsBegin(std::string fe_name) const {
    return numeredFiniteElements.get<FiniteElement_name_mi_tag>().lower_bound(fe_name);
  }

  NumeredEntFiniteElementbyName::iterator getNumeredFEsEnd(std::string fe_name) const {
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
    NumeredEntFiniteElementbyNameAndPart::iterator IT = MOFEMPROBLEM->getNumeredFEsBegin(NAME,PART); \
    IT!=MOFEMPROBLEM->getNumeredFEsEnd(NAME,PART); IT++

  NumeredEntFiniteElementbyNameAndPart::iterator getNumeredFEsBegin(std::string fe_name,int part) const {
    return numeredFiniteElements.get<Composite_Name_And_Part_mi_tag>().lower_bound(boost::make_tuple(fe_name,part));
  }

  NumeredEntFiniteElementbyNameAndPart::iterator getNumeredFEsEnd(std::string fe_name,int part) const {
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
    NumeredDofEntity_multiIndex::iterator IT = MOFEMPROBLEM->getNumeredDofsRowsBegin(); \
    IT!=MOFEMPROBLEM->getNumeredDofsRowsEnd(); IT++

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
    NumeredDofEntity_multiIndex::iterator IT = MOFEMPROBLEM->getNumeredDofsColsBegin(); \
    IT!=MOFEMPROBLEM->getNumeredDofsColsEnd(); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredDofsRowsBegin() const { return numered_dofs_rows->begin(); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredDofsRowsEnd() const { return numered_dofs_rows->end(); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredDofsColsBegin() const { return numered_dofs_cols->begin(); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredDofsColsEnd() const { return numered_dofs_cols->end(); }

  /**
    * \brief get iterator of dof in row by uid
    * \ingroup problems_multi_indices
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_UID_(MOFEMPROBLEM,UID,IT) \
    NumeredDofEntityByUId::iterator IT = MOFEMPROBLEM->get_row_dof_by_uid(UID);

  /**
    * \brief get iterator of dof in col by uid
    * \ingroup problems_multi_indices
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_UID_(MOFEMPROBLEM,UID,IT) \
    NumeredDofEntityByUId::iterator IT = MOFEMPROBLEM->get_col_dof_by_uid(UID);

  /// get iterator of dof in row by uid (instead you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_UID_FOR_LOOP_)
  NumeredDofEntityByUId::iterator get_row_dof_by_uid(GlobalUId uid) const {
    return numered_dofs_rows->get<Unique_mi_tag>().find(uid);
  };

  /// get iterator of dof in column by uid (instead you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_UID_FOR_LOOP_)
  NumeredDofEntityByUId::iterator get_col_dof_by_uid(GlobalUId uid) const {
    return numered_dofs_cols->get<Unique_mi_tag>().find(uid);
  };

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
    NumeredDofEntityByLocalIdx::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_by_locidx_begin(0); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_by_locidx_end(MOFEMPROBLEM->getNbLocalDofsRow()-1); IT++

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
    NumeredDofEntityByUId::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_by_locidx_begin(0); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_by_locidx_end(MOFEMPROBLEM->getNbLocalDofsRow()-1); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofEntityByLocalIdx::iterator get_numered_dofs_rows_by_locidx_begin(const DofIdx locidx) const
    { return numered_dofs_rows->get<PetscLocalIdx_mi_tag>().lower_bound(locidx); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofEntityByLocalIdx::iterator get_numered_dofs_rows_by_locidx_end(const DofIdx locidx) const
    { return numered_dofs_rows->get<PetscLocalIdx_mi_tag>().upper_bound(locidx); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofEntityByLocalIdx::iterator get_numered_dofs_cols_by_locidx_begin(const DofIdx locidx) const
    { return numered_dofs_cols->get<PetscLocalIdx_mi_tag>().lower_bound(locidx); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofEntityByLocalIdx::iterator get_numered_dofs_cols_by_locidx_end(const DofIdx locidx) const
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
    NumeredDofEntityByEnt::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_by_ent_begin(ENT); \
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
    NumeredDofEntityByEnt::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_by_ent_begin(ENT); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_by_ent_end(ENT); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofEntityByEnt::iterator get_numered_dofs_rows_by_ent_begin(const EntityHandle ent) const
    { return numered_dofs_rows->get<Ent_mi_tag>().lower_bound(ent); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofEntityByEnt::iterator get_numered_dofs_rows_by_ent_end(const EntityHandle ent) const
    { return numered_dofs_rows->get<Ent_mi_tag>().upper_bound(ent); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofEntityByEnt::iterator get_numered_dofs_cols_by_ent_begin(const EntityHandle ent) const
    { return numered_dofs_cols->get<Ent_mi_tag>().lower_bound(ent); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofEntityByEnt::iterator get_numered_dofs_cols_by_ent_end(const EntityHandle ent) const
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
    NumeredDofEntityByFieldName::iterator IT = MOFEMPROBLEM->getNumeredDofsRowsBegin(NAME); \
    IT!=MOFEMPROBLEM->getNumeredDofsRowsEnd(NAME); IT++

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
    NumeredDofEntityByFieldName::iterator IT = MOFEMPROBLEM->getNumeredDofsColsBegin(NAME); \
    IT!=MOFEMPROBLEM->getNumeredDofsColsEnd(NAME); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofEntityByFieldName::iterator getNumeredDofsRowsBegin(const std::string& name) const
    { return numered_dofs_rows->get<FieldName_mi_tag>().lower_bound(name); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofEntityByFieldName::iterator getNumeredDofsRowsEnd(const std::string& name) const
    { return numered_dofs_rows->get<FieldName_mi_tag>().upper_bound(name); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofEntityByFieldName::iterator getNumeredDofsColsBegin(const std::string& name) const
    { return numered_dofs_cols->get<FieldName_mi_tag>().lower_bound(name); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofEntityByFieldName::iterator getNumeredDofsColsEnd(const std::string& name) const
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
    NumeredDofEntityByNameEntAndPart::iterator IT = MOFEMPROBLEM->getNumeredDofsRowsBegin(NAME,ENT,PART); \
    IT!=MOFEMPROBLEM->getNumeredDofsRowsEnd(NAME,ENT,PART); IT++

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
    NumeredDofEntityByNameEntAndPart::iterator IT = MOFEMPROBLEM->getNumeredDofsColsBegin(NAME,ENT,PART); \
    IT!=MOFEMPROBLEM->getNumeredDofsColsEnd(NAME,ENT,PART); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_ for loops)
  NumeredDofEntityByNameEntAndPart::iterator getNumeredDofsRowsBegin(const std::string& name,const EntityHandle ent,const int part) const {
    return numered_dofs_rows->get<Composite_Name_Ent_And_Part_mi_tag>().lower_bound(boost::make_tuple(name,ent,part));
  }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_ for loops)
  NumeredDofEntityByNameEntAndPart::iterator getNumeredDofsRowsEnd(const std::string& name,const EntityHandle ent,const int part) const {
    return numered_dofs_rows->get<Composite_Name_Ent_And_Part_mi_tag>().upper_bound(boost::make_tuple(name,ent,part));
  }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_ENT_PART_FOR_LOOP_ for loops)
  NumeredDofEntityByNameEntAndPart::iterator getNumeredDofsColsBegin(const std::string& name,const EntityHandle ent,const int part) const {
    return numered_dofs_cols->get<Composite_Name_Ent_And_Part_mi_tag>().lower_bound(boost::make_tuple(name,ent,part));
  }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_ENT_PART_FOR_LOOP_ for loops)
  NumeredDofEntityByNameEntAndPart::iterator getNumeredDofsColsEnd(const std::string& name,const EntityHandle ent,const int part) const {
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
  NumeredDofEntity_multiIndex::iterator IT = MOFEMPROBLEM->getNumeredDofsRowsBegin(PART); \
  IT!=MOFEMPROBLEM->getNumeredDofsRowsEnd(PART); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_PART_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredDofsRowsBegin(const int part) const {
    return numered_dofs_rows->get<Unique_mi_tag>().lower_bound(
      DofEntity::getGlobalUniqueIdCalculate_Low_Proc(part)
    );
  }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_PART_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredDofsRowsEnd(const int part) const {
    return numered_dofs_rows->get<Unique_mi_tag>().upper_bound(
      DofEntity::getGlobalUniqueIdCalculate_Hi_Proc(part)
    );
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
  NumeredDofEntity_multiIndex::iterator IT = MOFEMPROBLEM->getNumeredDofsColsBegin(PART); \
  IT!=MOFEMPROBLEM->getNumeredDofsColsEnd(PART); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_PART_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredDofsColsBegin(const int part) const {
    return numered_dofs_rows->get<Unique_mi_tag>().lower_bound(
      DofEntity::getGlobalUniqueIdCalculate_Low_Proc(part)
    );
  }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_PART_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredDofsColsEnd(const int part) const {
    return numered_dofs_rows->get<Unique_mi_tag>().upper_bound(
      DofEntity::getGlobalUniqueIdCalculate_Hi_Proc(part)
    );
  }

  MoFEMProblem(Interface &moab,const EntityHandle meshset);

  virtual ~MoFEMProblem();

  inline BitProblemId getId() const { return *((BitProblemId*)tag_id_data); }

  inline std::string getName() const { return std::string((char *)tag_name_data,tag_name_size); }

  inline DofIdx getNbDofsRow() const { return *((DofIdx*)tag_nbdof_data_row); }
  inline DofIdx getNbDofsCol() const { return *((DofIdx*)tag_nbdof_data_col); }
  inline DofIdx getNbLocalDofsRow() const { return *((DofIdx*)tag_local_nbdof_data_row); }
  inline DofIdx getNbLocalDofsCol() const { return *((DofIdx*)tag_local_nbdof_data_col); }
  inline DofIdx getNbGhostDofsRow() const { return *((DofIdx*)tag_ghost_nbdof_data_row); }
  inline DofIdx getNbGhostDofsCol() const { return *((DofIdx*)tag_ghost_nbdof_data_col); }

  inline BitRefLevel getBitRefLevel() const { return *tag_BitRefLevel; }
  inline BitRefLevel getMaskBitRefLevel() const { return *tag_MaskBitRefLevel; }

  PetscErrorCode getRowDofsByPetscGlobalDofIdx(DofIdx idx,const NumeredDofEntity **dof_ptr) const;
  PetscErrorCode getColDofsByPetscGlobalDofIdx(DofIdx idx,const NumeredDofEntity **dof_ptr) const;

  BitFEId getBitFEId() const;

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

  /**
   * \brief Get weak_ptr reference to sequence/vector storing dofs.
   *
   * Vector is automatically destroy when last DOF in vector os destroyed. Every
   * shared_ptr to the DOF has aliased shared_ptr to vector of DOFs in that vector.
   * That do the trick.
   *
   * \note It is week_ptr, so it is no guaranteed that sequence is there. Check
   * if sequence is there.
   * \code
   * if(boost::shared_ptr<std::vector<NumeredDofEntity> > ptr=fe->getRowDofsSeqence().lock()) {
   *  // use ptr
   * }
   * \endcode
   *
   */
  inline boost::weak_ptr<std::vector<NumeredDofEntity> >& getRowDofsSeqence() const {
    return dofsRowSequence;
  }

  /**
   * \brief Get weak_ptr reference to sequence/vector storing dofs.
   *
   * Vector is automatically destroy when last DOF in vector os destroyed. Every
   * shared_ptr to the DOF has aliased shared_ptr to vector of DOFs in that vector.
   * That do the trick.
   *
   * \note It is week_ptr, so it is no guaranteed that sequence is there. Check
   * if sequence is there.
   * \code
   * if(boost::shared_ptr<std::vector<NumeredDofEntity> > ptr=fe->getColDofsSeqence().lock()) {
   *  // use ptr
   * }
   * \endcode
   *
   */
  inline boost::weak_ptr<std::vector<NumeredDofEntity> >& getColDofsSeqence() const {
    return dofsColSequence;
  }

  /**
   * \brief Get weak_ptr reference to sequence/vector storing finite elements.
   *
   * \note It is week_ptr, so it is no guaranteed that sequence is there. Check
   * if sequence is there.
   * \code
   * if(boost::shared_ptr<std::vector<NumeredDofEntity> > ptr=fe->getFeSeqence().lock()) {
   *  // use ptr
   * }
   * \endcode
   *
   */
  // inline boost::weak_ptr<std::vector<NumeredEntFiniteElement> >& getFeSeqence() const {
  //   return feSequence;
  // }


private:

  // Keep vector of DoFS on entity
  mutable boost::weak_ptr<std::vector<NumeredDofEntity> > dofsRowSequence;
  mutable boost::weak_ptr<std::vector<NumeredDofEntity> > dofsColSequence;

  // // Keeps finite elements on entities
  // mutable boost::weak_ptr<std::vector<NumeredEntFiniteElement> > feSequence;

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
      tag<Meshset_mi_tag>,
      member<MoFEMProblem,EntityHandle,&MoFEMProblem::meshset>
    >,
    hashed_unique<
      tag<BitProblemId_mi_tag>,
      const_mem_fun<MoFEMProblem,BitProblemId,&MoFEMProblem::getId>, HashBit<BitProblemId>, EqBit<BitProblemId>
    >,
    hashed_unique<
      tag<Problem_mi_tag>,
      const_mem_fun<MoFEMProblem,std::string,&MoFEMProblem::getName>
    >
  >
> MoFEMProblem_multiIndex;

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
  void operator()(MoFEMProblem &p) { *(p.tag_MaskBitRefLevel) = bit; };
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

}

#endif //__PROBLEMSMULTIINDICES_HPP__

/***************************************************************************//**
 * \defgroup problems_multi_indices Problems structures and multi-indices
 * \ingroup mofem
 ******************************************************************************/
