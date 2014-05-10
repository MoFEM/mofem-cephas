/** \file ProblemsMultiIndices.hpp
 * \brief Myltindex containes, data structures for problems and other low-level functions 
 * 
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 *
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
 * It can be freely used for educational and research purposes 
 * by other institutions. If you use this softwre pleas cite my work. 
 *
 * MoFEM is free software: you can redistribute it and/or modify it under
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

/// \brief keeps data about problem
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
  NumeredDofMoFEMEntity_multiIndex numered_dofs_rows;
  NumeredDofMoFEMEntity_multiIndex numered_dofs_cols;
  NumeredMoFEMFiniteElement_multiIndex numeredFiniteElements;

  /**
    * use with loops to iterate problem fes 
    *
    * for(_IT_NUMEREDFEMOFEMENTITY_FOR_LOOP_(MOFEMPROBLEM,NAME,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDFEMOFEMENTITY_BY_NAME_FOR_LOOP_(MOFEMPROBLEM,NAME,IT) \
    NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_fes_begin(NAME); \
    IT!=MOFEMPROBLEM->get_numered_fes_end(NAME); IT++

  NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator get_numered_fes_begin(string fe_name) const { 
    return numeredFiniteElements.get<MoFEMFiniteElement_name_mi_tag>().lower_bound(fe_name);
  }

  NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator get_numered_fes_end(string fe_name) const { 
    return numeredFiniteElements.get<MoFEMFiniteElement_name_mi_tag>().upper_bound(fe_name);
  }

  /**
    * use with loops to iterate problem fes 
    *
    * for(_IT_NUMEREDFEMOFEMENTITY_BY_NAME_AND_PART_FOR_LOOP_(MOFEMPROBLEM,NAME,PART,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDFEMOFEMENTITY_BY_NAME_AND_PART_FOR_LOOP_(MOFEMPROBLEM,NAME,PART,IT) \
    NumeredMoFEMFiniteElement_multiIndex::index<Composite_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_fes_begin(NAME,PART); \
    IT!=MOFEMPROBLEM->get_numered_fes_end(NAME,PART); IT++

  NumeredMoFEMFiniteElement_multiIndex::index<Composite_mi_tag>::type::iterator get_numered_fes_begin(string fe_name,int part) const { 
    return numeredFiniteElements.get<Composite_mi_tag>().lower_bound(boost::make_tuple(fe_name,part));
  }

  NumeredMoFEMFiniteElement_multiIndex::index<Composite_mi_tag>::type::iterator get_numered_fes_end(string fe_name,int part) const { 
    return numeredFiniteElements.get<Composite_mi_tag>().upper_bound(boost::make_tuple(fe_name,part));
  }

  /**
    * use with loops to iterate row dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(MOFEMPROBLEM,IT) \
    NumeredDofMoFEMEntity_multiIndex::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_begin(); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_end(); IT++

  /**
    * use with loops to iterate col dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(MOFEMPROBLEM,IT) \
    NumeredDofMoFEMEntity_multiIndex::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_begin(); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_end(); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::iterator get_numered_dofs_rows_begin() const { return numered_dofs_rows.begin(); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::iterator get_numered_dofs_rows_end() const { return numered_dofs_rows.end(); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::iterator get_numered_dofs_cols_begin() const { return numered_dofs_cols.begin(); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::iterator get_numered_dofs_cols_end() const { return numered_dofs_cols.end(); }

  /**
    * get iterator of dof in row by uid
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_UID_(MOFEMPROBLEM,UID,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_row_dof_by_uid(UID);

  /**
    * get iterator of dof in col by uid
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_UID_(MOFEMPROBLEM,UID,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_col_dof_by_uid(UID);

  /// get iterator of dof in row by uid (instead you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_UID_FOR_LOOP_)
  NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator get_row_dof_by_uid(UId uid) const { return numered_dofs_rows.get<Unique_mi_tag>().find(uid); };

  /// get iterator of dof in column by uid (instead you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_UID_FOR_LOOP_)
  NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator get_col_dof_by_uid(UId uid) const { return numered_dofs_cols.get<Unique_mi_tag>().find(uid); };

  /**
    * use with loops to iterate row dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_BY_LOCIDX_ROW_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_LOCIDX_FOR_LOOP_(MOFEMPROBLEM,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_by_locidx_begin(0); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_by_locidx_end(MOFEMPROBLEM->get_nb_local_dofs_row()); IT++

  /**
    * use with loops to iterate col dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_LOCIDX_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_LOCIDX_FOR_LOOP_(MOFEMPROBLEM,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_by_locidx_begin(0); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_by_locidx_end(MOFEMPROBLEM->get_nb_local_dofs_row()); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator get_numered_dofs_rows_by_locidx_begin(const DofIdx locidx) const 
    { return numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(locidx); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator get_numered_dofs_rows_by_locidx_end(const DofIdx locidx) const 
    { return numered_dofs_rows.get<PetscLocalIdx_mi_tag>().upper_bound(locidx); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator get_numered_dofs_cols_by_locidx_begin(const DofIdx locidx) const 
    { return numered_dofs_cols.get<PetscLocalIdx_mi_tag>().lower_bound(locidx); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator get_numered_dofs_cols_by_locidx_end(const DofIdx locidx) const 
    { return numered_dofs_cols.get<PetscLocalIdx_mi_tag>().upper_bound(locidx); }

  /**
    * use with loops to iterate row dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_BY_ENT_ROW_FOR_LOOP_(MOFEMPROBLEM,ENT,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(MOFEMPROBLEM,ENT,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_by_ent_begin(ENT); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_by_ent_end(ENT); IT++

  /**
    * use with loops to iterate col dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_(MOFEMPROBLEM,ENT,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_(MOFEMPROBLEM,ENT,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_by_ent_begin(ENT); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_by_ent_end(ENT); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator get_numered_dofs_rows_by_ent_begin(const EntityHandle ent) const 
    { return numered_dofs_rows.get<MoABEnt_mi_tag>().lower_bound(ent); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator get_numered_dofs_rows_by_ent_end(const EntityHandle ent) const 
    { return numered_dofs_rows.get<MoABEnt_mi_tag>().upper_bound(ent); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator get_numered_dofs_cols_by_ent_begin(const EntityHandle ent) const 
    { return numered_dofs_cols.get<MoABEnt_mi_tag>().lower_bound(ent); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator get_numered_dofs_cols_by_ent_end(const EntityHandle ent) const 
    { return numered_dofs_cols.get<MoABEnt_mi_tag>().upper_bound(ent); }

  /**
    * use with loops to iterate row dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_BY_NAME_ROW_FOR_LOOP_(MOFEMPROBLEM,NAME,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(MOFEMPROBLEM,NAME,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_begin(NAME); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_end(NAME); IT++

  /**
    * use with loops to iterate col dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_FOR_LOOP_(MOFEMPROBLEM,NAME,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_FOR_LOOP_(MOFEMPROBLEM,NAME,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_begin(NAME); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_end(NAME); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_numered_dofs_rows_begin(const string& name) const 
    { return numered_dofs_rows.get<FieldName_mi_tag>().lower_bound(name); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_numered_dofs_rows_end(const string& name) const 
    { return numered_dofs_rows.get<FieldName_mi_tag>().upper_bound(name); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_numered_dofs_cols_begin(const string& name) const 
    { return numered_dofs_cols.get<FieldName_mi_tag>().lower_bound(name); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_numered_dofs_cols_end(const string& name) const 
    { return numered_dofs_cols.get<FieldName_mi_tag>().upper_bound(name); }

  /**
    * use with loops to iterate row dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_BY_NAME_ENT_PART_ROW_FOR_LOOP_(MOFEMPROBLEM,NAME,ENT,PART,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(MOFEMPROBLEM,NAME,ENT,PART,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_begin(NAME,ENT,PART); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_end(NAME,ENT,PART); IT++

  /**
    * use with loops to iterate col dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_ENT_PART_FOR_LOOP_(MOFEMPROBLEM,NAME,ENT,PART,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_ENT_PART_FOR_LOOP_(MOFEMPROBLEM,NAME,ENT,PART,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_begin(NAME,ENT,PART); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_end(NAME,ENT,PART); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator get_numered_dofs_rows_begin(const string& name,const EntityHandle ent,const int part) const 
    { return numered_dofs_rows.get<Composite_Name_Ent_And_Part_mi_tag>().lower_bound(boost::make_tuple(name,ent,part)); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator get_numered_dofs_rows_end(const string& name,const EntityHandle ent,const int part) const 
    { return numered_dofs_rows.get<Composite_Name_Ent_And_Part_mi_tag>().upper_bound(boost::make_tuple(name,ent,part)); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator get_numered_dofs_cols_begin(const string& name,const EntityHandle ent,const int part) const 
    { return numered_dofs_cols.get<Composite_Name_Ent_And_Part_mi_tag>().lower_bound(boost::make_tuple(name,ent,part)); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_NAME_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator get_numered_dofs_cols_end(const string& name,const EntityHandle ent,const int part) const 
    { return numered_dofs_cols.get<Composite_Name_Ent_And_Part_mi_tag>().upper_bound(boost::make_tuple(name,ent,part)); }

  MoFEMProblem(Interface &moab,const EntityHandle _meshset);
  inline BitProblemId get_id() const { return *((BitProblemId*)tag_id_data); }
  inline string get_name() const { return string((char *)tag_name_data,tag_name_size); }
  inline DofIdx get_nb_dofs_row() const { return *((DofIdx*)tag_nbdof_data_row); }
  inline DofIdx get_nb_dofs_col() const { return *((DofIdx*)tag_nbdof_data_col); }
  inline DofIdx get_nb_local_dofs_row() const { return *((DofIdx*)tag_local_nbdof_data_row); }
  inline DofIdx get_nb_local_dofs_col() const { return *((DofIdx*)tag_local_nbdof_data_col); }
  inline DofIdx get_nb_ghost_dofs_row() const { return *((DofIdx*)tag_ghost_nbdof_data_row); }
  inline DofIdx get_nb_ghost_dofs_col() const { return *((DofIdx*)tag_ghost_nbdof_data_col); }
  inline BitRefLevel get_BitRefLevel() const { return *tag_BitRefLevel; }
  inline BitRefLevel get_DofMask_BitRefLevel() const { return *tag_BitRefLevel_DofMask; }
  BitFEId get_BitFEId() const;
  friend ostream& operator<<(ostream& os,const MoFEMProblem& e);
};

typedef multi_index_container<
  MoFEMProblem,
  indexed_by<
    ordered_unique<
      tag<Meshset_mi_tag>, member<MoFEMProblem,EntityHandle,&MoFEMProblem::meshset> >,
    hashed_unique<
      tag<BitProblemId_mi_tag>, const_mem_fun<MoFEMProblem,BitProblemId,&MoFEMProblem::get_id>, hashbit<BitProblemId>, eqbit<BitProblemId> >,
    hashed_unique<
      tag<MoFEMProblem_mi_tag>, const_mem_fun<MoFEMProblem,string,&MoFEMProblem::get_name> >
  > > MoFEMProblem_multiIndex;

/// \brief add ref level to problem
struct problem_change_ref_level_bit_add {
  BitRefLevel bit;
  problem_change_ref_level_bit_add(const BitRefLevel _bit): bit(_bit) {};
  void operator()(MoFEMProblem &p) { *(p.tag_BitRefLevel) |= bit; };
};

/// \brief set ref level to problem
struct problem_change_ref_level_bit_set {
  BitRefLevel bit;
  problem_change_ref_level_bit_set(const BitRefLevel _bit): bit(_bit) {};
  void operator()(MoFEMProblem &p) { *(p.tag_BitRefLevel) = bit; };
};

/// \brief set prof dof bit ref mask
struct problem_change_ref_level_bit_dof_mask_set {
  BitRefLevel bit;
  problem_change_ref_level_bit_dof_mask_set(const BitRefLevel _bit): bit(_bit) {};
  void operator()(MoFEMProblem &p) { *(p.tag_BitRefLevel_DofMask) = bit; };
};

/// \brief add finite element to problem
struct problem_MoFEMFiniteElement_change_bit_add {
  BitFEId f_id;
  problem_MoFEMFiniteElement_change_bit_add(const BitFEId _f_id): f_id(_f_id) {};
  void operator()(MoFEMProblem &p);
};

/// \brief increase nb. dof in row
struct problem_row_change {
  const DofMoFEMEntity *dof_ptr;
  problem_row_change(const DofMoFEMEntity *_dof_ptr);
  void operator()(MoFEMProblem &e);
};

/// \brief increase nb. dof in col
struct problem_col_change {
  const DofMoFEMEntity *dof_ptr;
  problem_col_change(const DofMoFEMEntity *_dof_ptr);
  void operator()(MoFEMProblem &e);
};

/// \brief zero nb. of dofs in row
struct problem_zero_nb_rows_change {
  void operator()(MoFEMProblem &e);
};

/// \brief zero nb. of dofs in col
struct problem_zero_nb_cols_change {
  void operator()(MoFEMProblem &e);
};

/// \brief clear problem finite elements 
struct problem_clear_numered_finite_elements_change {
  void operator()(MoFEMProblem &e);
};

/// \brief number dofs in row
struct problem_row_number_change {
  problem_row_number_change() {};
  void operator()(MoFEMProblem &e);
};

/// \brief number dofs in col
struct problem_col_number_change {
  problem_col_number_change() {};
  void operator()(MoFEMProblem &e);
};

}

#endif //__PROBLEMSMULTIINDICES_HPP__
