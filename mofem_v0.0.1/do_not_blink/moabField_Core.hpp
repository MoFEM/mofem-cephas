/** \file moabField_Core.hpp
 * \brief Core moabField class for user interface
 * 
 * Low level data structures not used directly by user
 *
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
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

#ifndef __MOABFIELD_CORE_HPP__
#define __MOABFIELD_CORE_HPP__

#include "moabField.hpp"
#include "Core_dataStructures.hpp"

namespace MoFEM {

/** \brief Core moabField class
 *
 * This clas is not used directly by the ueser
 */
struct moabField_Core: public moabField {
  ErrorCode rval;
  PetscErrorCode ierr;
  //Data and low level methods 
  Tag th_Part;
  Tag th_RefType,th_RefParentHandle,th_RefBitLevel,th_RefBitEdge;
  Tag th_FieldId,th_FieldName,th_FieldSpace;
  Tag th_FEId,th_FEName;
  Tag th_FEIdCol,th_FEIdRow,th_FEIdData;
  Tag th_ProblemId,th_ProblemName,th_ProblemFEId;
  Tag th_ProblemNbDofsRow,th_ProblemNbDofsCol;
  Tag th_ProblemLocalNbDofRow,th_ProblemGhostNbDofRow;
  Tag th_ProblemLocalNbDofCol,th_ProblemGhostNbDofCol;
  Tag th_ProblemShift,th_FieldShift,th_FEShift;
  Tag nsTag,ssTag,nsTag_data,ssTag_data,bhTag,bhTag_header;
  Tag th_ElemType;

  Interface& moab;
  int *f_shift,*MoFEMFE_shift,*p_shift;
  int verbose;

  //database

  //ref
  RefMoFEMEntity_multiIndex ref_entities;
  RefMoFEMFiniteElement_multiIndex ref_finite_elements;
  //field
  MoFEMField_multiIndex moabfields;
  MoFEMEntity_multiIndex ents_moabfield;
  DofMoFEMEntity_multiIndex dofs_moabfield;
  //finite element
  MoFEMFE_multiIndex finite_elements;
  EntMoFEMFE_multiIndex finite_elements_data;
  //finite elemts and dofs
  MoFEMAdjacencies_multiIndex adjacencies;
  //problems
  MoFEMProblem_multiIndex problems;
  //prism 
  AdjBasicMoFEMEntity_multiIndex Adj_prisms;
  //cubit
  moabBaseMeshSet_multiIndex cubit_meshsets;

  //safty nets
  Tag th_MoFEMBuild;
  int *build_MoFEM;

  //core metgids 
  PetscErrorCode clear_map();
  BitFieldId get_field_shift();
  BitFEId get_BitFEId();
  BitProblemId get_problem_shift();
  PetscErrorCode map_from_mesh(int verb = -1);
  Interface& get_moab();

  //meshsets
  PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const Cubit_BC_bitset CubitBCType, const int dimension,Range &entities,const bool recursive = false);
  PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const Cubit_BC_bitset CubitBCType, Range &entities,const bool recursive = false);
  PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType, const int dimension,Range &entities,const bool recursive = false);
  PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType, Range &entities,const bool recursive = false);
  PetscErrorCode get_msId_meshset(const int msId,const unsigned int CubitBCType,EntityHandle &meshset);
  PetscErrorCode get_CubitBCType_meshsets(const unsigned int CubitBCType,Range &meshsets);

  //refine
  PetscErrorCode seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit);
  PetscErrorCode seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit);
  PetscErrorCode add_verices_in_the_middel_of_edges(
    const EntityHandle meshset,const BitRefLevel &bit,const bool recursive = false,int verb = -1);
  PetscErrorCode refine_TET(const EntityHandle meshset,const BitRefLevel &bit,const bool respect_interface = true);
  PetscErrorCode refine_PRISM(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode refine_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,const bool recursive = false,int verb = -1);
  PetscErrorCode refine_get_finite_elements(const BitRefLevel &bit,const EntityHandle meshset);
  PetscErrorCode refine_get_ents(const BitRefLevel &bit,const EntityHandle meshset);
  PetscErrorCode refine_get_childern(
    const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child, EntityType child_type,
    const bool recursive = false, int verb = -1);

  //field
  PetscErrorCode add_field(const string& name,const BitFieldId id,const FieldSpace space,const ApproximationRank rank,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TETs(const EntityHandle meshset,const BitFieldId id);
  PetscErrorCode add_ents_to_field_by_TETs(const EntityHandle meshset,const string& name);
  PetscErrorCode add_ents_to_field_by_PRISMs(const EntityHandle meshset,const BitFieldId id);
  PetscErrorCode add_ents_to_field_by_PRISMs(const EntityHandle meshset,const string& name);
  PetscErrorCode set_field_order(const EntityHandle meshset,const EntityType type,const BitFieldId id,const ApproximationOrder order);
  PetscErrorCode set_field_order(const EntityHandle meshset,const EntityType type,const string& name,const ApproximationOrder order);
  PetscErrorCode dofs_NoField(const BitFieldId id);
  PetscErrorCode dofs_L2H1HcurlHdiv(const BitFieldId id,int verb = -1);
  PetscErrorCode list_dof_by_id(const BitFieldId id) const;
  PetscErrorCode list_ent_by_id(const BitFieldId id) const;
  PetscErrorCode list_field() const;
  PetscErrorCode add_field(const string& name,const FieldSpace space,const ApproximationRank rank,int verb = -1);
  BitFieldId get_BitFieldId(const string& name) const;
  string get_BitFieldId_name(const BitFieldId id) const;
  EntityHandle get_field_meshset(const BitFieldId id) const;
  EntityHandle get_field_meshset(const string& name) const;

  //MoFEMFE
  PetscErrorCode add_finite_element(const string &MoFEMFE_name);
  PetscErrorCode modify_finite_element_add_field_data(const string &MoFEMFE_name,const string &name_filed);
  PetscErrorCode modify_finite_element_add_field_row(const string &MoFEMFE_name,const string &name_row);
  PetscErrorCode modify_finite_element_add_field_col(const string &MoFEMFE_name,const string &name_col);
  PetscErrorCode add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const BitFEId id);
  PetscErrorCode add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const string &name);
  PetscErrorCode add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const string& name);
  PetscErrorCode add_ents_to_finite_element_by_MESHSETs(const EntityHandle meshset,const string& name);
  PetscErrorCode add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit_ref,const string &name,EntityType type);
  BitFEId get_BitFEId(const string& name) const;
  string get_BitFEId_name(const BitFEId id) const;
  EntityHandle get_meshset_by_BitFEId(const BitFEId id) const;
  EntityHandle get_meshset_by_BitFEId(const string& name) const;
  PetscErrorCode list_finite_elements() const;

  //problem
  PetscErrorCode add_problem(const BitProblemId id,const string& name);
  PetscErrorCode add_problem(const string& name);
  PetscErrorCode modify_problem_add_finite_element(const string &name_problem,const string &MoFEMFE_name);
  PetscErrorCode modify_problem_ref_level_add_bit(const string &name_problem,const BitRefLevel &bit);
  BitProblemId get_BitProblemId(const string& name) const;
  PetscErrorCode list_problem() const;

  //build problem, adjacencies,MoFEMFE and dofs
  PetscErrorCode build_fields(int verb = -1);

  ///add entity EntFe to finite element data databse and resolve dofs on that entity
  PetscErrorCode build_finite_element(const EntMoFEMFE &EntFe,int verb = -1);

  //loop over all finite elements, resolve its meshsets, and resolve dofs on that entities
  PetscErrorCode build_finite_elements(int verb = -1);

  PetscErrorCode build_adjacencies(const BitRefLevel bit);
  PetscErrorCode build_problems(int verb = -1);

  //adjacencies
  PetscErrorCode list_adjacencies() const;

  //problem buildig
  PetscErrorCode partition_problems(const string &name,int verb = -1);
  PetscErrorCode partition_ghost_dofs(const string &name);
  PetscErrorCode partition_finite_elements(const string &name,bool do_skip = true,int verb = -1);

  //clean active
  PetscErrorCode erase_inactive_dofs_moabfield();
  PetscErrorCode erase_inacrive_dofs_numered();

  //save meshsets
  PetscErrorCode problem_get_FE(const string &name,const string &fe_name,const EntityHandle meshset);

  //vector and matrices 
  PetscErrorCode VecCreateGhost(const string &name,RowColData rc,Vec *V);
  PetscErrorCode set_local_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode);
  PetscErrorCode set_global_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode);
  PetscErrorCode MatCreateMPIAIJWithArrays(const string &name,Mat *Aij,int verb = -1);

  //topology
  PetscErrorCode get_msId_3dENTS_sides(const int msId,const Cubit_BC_bitset CubitBCType,
    const bool recursive = false,int verb = -1);
  PetscErrorCode get_msId_3dENTS_sides(const EntityHandle SideSet,
    const bool recursive = false,int verb = -1);
  PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const int msId,const Cubit_BC_bitset CubitBCType,
    const bool add_iterfece_entities,const bool recursive = false,int verb = -1);
  PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const EntityHandle SideSet,const bool add_iterfece_entities,const bool recursive = false,int verb = -1);
  PetscErrorCode add_prism_to_Adj_prisms(const EntityHandle prism,int verb = -1);

  //loops
  PetscErrorCode loop_finite_elements(
    const string &problem_name,const string &fe_name,FEMethod &method,
    int lower_rank,int upper_rank,int verb = -1);
  PetscErrorCode loop_finite_elements(const string &problem_name,const string &fe_name,FEMethod &method,int verb = -1);
  PetscErrorCode loop_dofs(const string &problem_name,const string &field_name,RowColData rc,EntMethod &method,int verb = -1);

  //get multi_index form database
  PetscErrorCode get_problems_database(const string &problem_name,const MoFEMProblem **problem_ptr);
  PetscErrorCode get_dofs_moabfield(const DofMoFEMEntity_multiIndex **dofs_moabfield_ptr);

  //Copy Field to Another
  //NOT TESTED DONT USE PetscErrorCode set_other_filed_values(const string& fiel_name,const string& cpy_field_name,InsertMode mode,ScatterMode scatter_mode);
  
  //Copy Vector of Field to Another
  PetscErrorCode set_other_global_VecCreateGhost(const string &name,const string& fiel_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode);


  //constructor
  moabField_Core(Interface& _moab,int _verbose = 1);
  ~moabField_Core();

  //Templates
  template<typename Tag> 
  PetscErrorCode partition_create_Mat(const string &name,Mat *Adj,Mat *Aij,const bool no_diagonals = true,int verb = -1) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    typedef typename boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Tag>::type NumeredDofMoFEMEntitys_by_idx;
    typedef NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type NumeredDofMoFEMEntitys_by_unique_id;
    typedef MoFEMAdjacencies_multiIndex::index<Composite_mi_tag>::type adj_by_ent;
    //find p_miit
    typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
    problems_by_name &problems_set = problems.get<MoFEMProblem_mi_tag>();
    problems_by_name::iterator p_miit = problems_set.find(name);
    //
    const NumeredDofMoFEMEntitys_by_idx &dofs_row_by_idx = p_miit->numered_dofs_rows.get<Tag>();
    const NumeredDofMoFEMEntitys_by_unique_id &dofs_col_by_id = p_miit->numered_dofs_cols.get<Unique_mi_tag>();
    DofIdx nb_dofs_row = dofs_row_by_idx.size();
    assert(p_miit->get_nb_dofs_row()==nb_dofs_row);
    typename boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Tag>::type::iterator miit_row,hi_miit_row;
    if(Tag::IamNotPartitioned) {
      DofIdx nb_dofs_row_on_proc = (DofIdx)ceil(nb_dofs_row/pcomm->size());
      DofIdx lower_dof_row = nb_dofs_row_on_proc*pcomm->rank();
      miit_row = dofs_row_by_idx.lower_bound(lower_dof_row);
      DofIdx upper_dof_row = pcomm->rank()==pcomm->size()-1 ? nb_dofs_row-1 : nb_dofs_row_on_proc*(pcomm->rank()+1)-1;
      hi_miit_row = dofs_row_by_idx.upper_bound(upper_dof_row);
      if(verb > 1) {
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"partition_create_Mat: row lower %d row upper %d\n",lower_dof_row,upper_dof_row);
	PetscSynchronizedFlush(PETSC_COMM_WORLD); 
      }
    } else {
      miit_row = dofs_row_by_idx.lower_bound(pcomm->rank());
      hi_miit_row = dofs_row_by_idx.upper_bound(pcomm->rank());
    }
    MoFEMEntity *MoFEMEntity_ptr = NULL;
    vector<DofIdx> dofs_vec,dofs_vec2;
    vector<PetscInt> i,j;
    // loop local rows
    for(;miit_row!=hi_miit_row;miit_row++) {
      i.push_back(j.size());
      if( (MoFEMEntity_ptr == NULL) ? 1 : (MoFEMEntity_ptr->get_unique_id() != miit_row->field_ptr->field_ptr->get_unique_id()) ) {
	// get field ptr
	MoFEMEntity_ptr = const_cast<MoFEMEntity*>(miit_row->field_ptr->field_ptr);
	adj_by_ent::iterator adj_miit = adjacencies.get<Composite_mi_tag>().lower_bound(boost::make_tuple(MoFEMEntity_ptr->get_meshset(),MoFEMEntity_ptr->get_ent()));
	adj_by_ent::iterator hi_adj_miit = adjacencies.get<Composite_mi_tag>().upper_bound(boost::make_tuple(MoFEMEntity_ptr->get_meshset(),MoFEMEntity_ptr->get_ent()));
	dofs_vec.resize(0);
	for(;adj_miit!=hi_adj_miit;adj_miit++) {
	  if(!(adj_miit->by_other&by_row)) continue;
	  if((adj_miit->EntMoFEMFE_ptr->get_id()&p_miit->get_BitFEId()).none()) continue;
	  if((adj_miit->EntMoFEMFE_ptr->get_BitRefLevel()&miit_row->get_BitRefLevel()).none()) continue;
	  int size  = adj_miit->EntMoFEMFE_ptr->tag_col_uids_size/sizeof(UId);
	  for(int ii = 0;ii<size;ii++) {
	    UId uid = adj_miit->EntMoFEMFE_ptr->tag_col_uids_data[ii];
	    NumeredDofMoFEMEntitys_by_unique_id::iterator miiit = dofs_col_by_id.find(uid);
	    if(miiit == p_miit->numered_dofs_cols.get<Unique_mi_tag>().end()) continue;
	    dofs_vec.insert(dofs_vec.end(),Tag::get_index(miiit));
	  }
	}
	sort(dofs_vec.begin(),dofs_vec.end());
	vector<DofIdx>::iterator new_end = unique(dofs_vec.begin(),dofs_vec.end());
	dofs_vec.erase(new_end,dofs_vec.end());
      }
      if(!dofs_vec.empty()) {
	if(no_diagonals) {
	  dofs_vec2.resize(0);
	  dofs_vec2.insert(dofs_vec2.end(),dofs_vec.begin(),dofs_vec.end());
	  vector<DofIdx>::iterator vit = find(dofs_vec2.begin(),dofs_vec2.end(),Tag::get_index(miit_row));
	  if(vit==dofs_vec2.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  dofs_vec2.erase(vit);
	  j.insert(j.end(),dofs_vec2.begin(),dofs_vec2.end());
	} else {
	  j.insert(j.end(),dofs_vec.begin(),dofs_vec.end());
	}
      }
    }
    //build adj matrix
    i.push_back(j.size());
    PetscInt *_i,*_j;
    PetscMalloc(i.size()*sizeof(PetscInt),&_i);
    PetscMalloc(j.size()*sizeof(PetscInt),&_j);
    copy(i.begin(),i.end(),_i);
    copy(j.begin(),j.end(),_j);
    PetscInt nb_row_dofs = p_miit->get_nb_dofs_row();
    PetscInt nb_col_dofs = p_miit->get_nb_dofs_col();
    if(Adj!=NULL) {
      ierr = MatCreateMPIAdj(PETSC_COMM_WORLD,i.size()-1,nb_col_dofs,_i,_j,PETSC_NULL,Adj); CHKERRQ(ierr);
    }
    if(Aij!=NULL) {
      PetscInt nb_local_dofs_row = p_miit->get_nb_local_dofs_row();
      if((unsigned int)nb_local_dofs_row!=i.size()-1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      PetscInt nb_local_dofs_col = p_miit->get_nb_local_dofs_col();
      ierr = ::MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD,nb_local_dofs_row,nb_local_dofs_col,nb_row_dofs,nb_col_dofs,_i,_j,PETSC_NULL,Aij); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }
  
  //low level finite element data
  double diffN_TET[12]; 

  //Petsc Logs
  PetscLogEvent USER_EVENT_preProcess;
  PetscLogEvent USER_EVENT_operator;
  PetscLogEvent USER_EVENT_postProcess;
};

}

#endif // __MOABFIELD_CORE_HPP__
