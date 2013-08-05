/** \file Core_dataStructures.hpp
 * \brief Myltindex containes, data structures and other low-level functions 
 * 
 * Low level data structures not used directly by user
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

#ifndef __DATASTRUCTURES_HPP__
#define __DATASTRUCTURES_HPP__

#include "common.hpp"

namespace MoFEM {

const int prism_adj_edges[] = { 6,7,8, -1,-1,-1, 0,1,2 };
const int prism_edges_conn[6][2] = { {0,1},{1,2},{2,0}, {3,4}, {4,5}, {5,3} };

inline int fNBENTITYSET_nofield(int P) { (void)P; return 1; }
//
inline int fNBVERTEX_L2(int P) { (void)P; return 0; }
inline int fNBEDGE_L2(int P) { (void)P; return 0; }
inline int fNBFACE_L2(int P) { (void)P; return 0; }
inline int fNBVOLUME_L2(int P) { return NBVOLUME_L2(P); }
//2D
inline int fNBSURFACE_L2(int P) { return NBSURFACE_L2(P); }

//
/// number of approx. functions for H1 space on vertex
inline int fNBVERTEX_H1(int P) { return (P==1) ? 1 : 0; }
/// number of approx. functions for H1 space on edge
inline int fNBEDGE_H1(int P) { return NBEDGE_H1(P); }
/// number of approx. functions for H1 space on face
inline int fNBFACE_H1(int P) { return NBFACE_H1(P); }
/// number of approx. functions for H1 space on volume
inline int fNBVOLUME_H1(int P) { return NBVOLUME_H1(P); }
//
/// number of approx. functions for Hcurl space on vertex
inline int fNBVERTEX_Hcurl(int P) { (void)P; return 0; }
inline int fNBEDGE_Hcurl(int P) { return NBEDGE_Hcurl(P); }
inline int fNBFACE_Hcurl(int P) { return NBFACE_Hcurl(P); }
inline int fNBVOLUME_Hcurl(int P) { return NBVOLUME_Hcurl(P); }
//
/// \brief number of approx. functions for Hdiv space on vertex
///
/// zero number of digrees of freedom on vertex for that space
inline int fNBVERTEX_Hdiv(int P) { (void)P; return 0; }
/// number of approx. functions for Hdiv space on edge
inline int fNBEDGE_Hdiv(int P) { assert(P==P); (void)P; return NBEDGE_Hdiv(P); }
/// number of approx. functions for Hcurl space on face
inline int fNBFACE_Hdiv(int P) { return NBFACE_Hdiv(P); }
/// number of approx. functions for Hcurl space on voulem
inline int fNBVOLUME_Hdiv(int P) { return NBVOLUME_Hdiv(P); }

//MultiIndex Tags
struct BitFieldId_col_id_mi_tag {};
struct BitFieldId_row_id_mi_tag {};
struct AppOrder_mi_tag {};
struct Active_mi_tag {};
struct BitProblemId_mi_tag {};
struct MoFEMProblem_mi_tag {};
struct MoFEMFE_Part_mi_tag {};
struct BitRefLevel_mi_tag {};
struct ParentEntType_mi_tag {};
struct MeshSetEnt_mi_tag{};
struct CubitMeshSets_mi_tag {};

struct ltstr
{ inline bool operator()(const string &s1, const string& s2) const
  { return strcmp(s1.c_str(), s2.c_str()) < 0; } };

/**
 * \typedef moabBaseMeshSet_multiIndex
 * type multiIndex container for CubitMeshSets
 *
 */
typedef multi_index_container<
  CubitMeshSets,
  indexed_by<
    hashed_unique<
      tag<Meshset_mi_tag>, member<CubitMeshSets,EntityHandle,&CubitMeshSets::meshset> >,
    ordered_non_unique<
      tag<CubitMeshSets_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_ulong> >,
    hashed_unique<
      tag<Composite_mi_tag>,       
      composite_key<
	CubitMeshSets,
	  const_mem_fun<CubitMeshSets,int,&CubitMeshSets::get_msId>,
	  const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_ulong> > >
  > > moabBaseMeshSet_multiIndex;

struct RefMoFEMEntity_change_add_bit {
  BitRefLevel bit;
  RefMoFEMEntity_change_add_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(RefMoFEMEntity &e) { 
    bit |= *(e.tag_BitRefLevel); 
    *e.tag_BitRefLevel = bit;
  }
};

/**
 * \typedef RefMoFEMEntity_multiIndex
 * type multiIndex container for RefMoFEMEntity
 *
 * \param hashed_unique MoABEnt_mi_tag 
 * \param ordered_non_unique Meshset_mi_tag 
 * \param ordered_non_unique MoABEnt_MoABEnt_mi_tag
 * \param ordered_non_unique EntType_mi_tag
 * \param ordered_non_unique ParentEntType_mi_tag
 * \param ordered_non_unique Composite_mi_tag
 */
typedef multi_index_container<
  RefMoFEMEntity,
  indexed_by<
    hashed_unique<
      tag<MoABEnt_mi_tag>, member<RefMoFEMEntity::BasicMoFEMEntity,EntityHandle,&RefMoFEMEntity::ent> >,
    ordered_non_unique<
      tag<MoABEnt_MoABEnt_mi_tag>, const_mem_fun<RefMoFEMEntity,EntityHandle,&RefMoFEMEntity::get_parent_ent> >,
    ordered_non_unique<
      tag<EntType_mi_tag>, const_mem_fun<RefMoFEMEntity::BasicMoFEMEntity,EntityType,&RefMoFEMEntity::get_ent_type> >,
    ordered_non_unique<
      tag<ParentEntType_mi_tag>, const_mem_fun<RefMoFEMEntity,EntityType,&RefMoFEMEntity::get_parent_ent_type> >,
    ordered_non_unique<
      tag<Composite_mi_tag>, 
      composite_key<
	RefMoFEMEntity,
	const_mem_fun<RefMoFEMEntity::BasicMoFEMEntity,EntityType,&RefMoFEMEntity::get_ent_type>,
	const_mem_fun<RefMoFEMEntity,EntityType,&RefMoFEMEntity::get_parent_ent_type> > >,
    ordered_non_unique<
      tag<Composite_mi_tag2>, 
      composite_key<
	RefMoFEMEntity,
	const_mem_fun<RefMoFEMEntity,EntityHandle,&RefMoFEMEntity::get_parent_ent>,
	const_mem_fun<RefMoFEMEntity::BasicMoFEMEntity,EntityType,&RefMoFEMEntity::get_ent_type> > >
  > > RefMoFEMEntity_multiIndex;

typedef multi_index_container<
  const RefMoFEMEntity*,
  indexed_by<
    hashed_unique<
      const_mem_fun<RefMoFEMEntity,EntityHandle,&RefMoFEMEntity::get_parent_ent> >
  > > RefMoFEMEntity_multiIndex_view_by_parent_entity;

struct ptrWrapperRefMoFEMFiniteElement: public interface_RefMoFEMFiniteElement<RefMoFEMFiniteElement> {
  typedef interface_RefMoFEMEntity<RefMoFEMFiniteElement> interface_type_RefMoFEMEntity;
  typedef interface_RefMoFEMFiniteElement<RefMoFEMFiniteElement> interface_type_RefMoFEMFiniteElement;
  int wrapp;
  ptrWrapperRefMoFEMFiniteElement(const RefMoFEMFiniteElement *__ptr): interface_RefMoFEMFiniteElement<RefMoFEMFiniteElement>(__ptr),wrapp(1) {}
  ptrWrapperRefMoFEMFiniteElement(const ptrWrapperRefMoFEMFiniteElement &ref): interface_RefMoFEMFiniteElement<RefMoFEMFiniteElement>(ref) { 
    wrapp = 1;
    assert(ref.wrapp == 1);
    (const_cast<ptrWrapperRefMoFEMFiniteElement&>(ref)).wrapp++;
  }
  ~ptrWrapperRefMoFEMFiniteElement() { 
    if(wrapp == 1) {
      delete interface_RefMoFEMEntity<RefMoFEMFiniteElement>::ref_ptr; 
    }
  }
};

/**
 * \typedef RefMoFEMFiniteElement
 * type multiIndex container for RefMoFEMFiniteElement
 *
 * \param hashed_unique MoABEnt_mi_tag 
 * \param ordered_non_unique Meshset_mi_tag 
 * \param ordered_non_unique MoABEnt_MoABEnt_mi_tag
 */
typedef multi_index_container<
  ptrWrapperRefMoFEMFiniteElement,
  indexed_by<
    hashed_unique<
      tag<MoABEnt_mi_tag>, const_mem_fun<ptrWrapperRefMoFEMFiniteElement::interface_type_RefMoFEMEntity,EntityHandle,&ptrWrapperRefMoFEMFiniteElement::get_ref_ent> >,
    ordered_non_unique<
      tag<MoABEnt_MoABEnt_mi_tag>, const_mem_fun<ptrWrapperRefMoFEMFiniteElement::interface_type_RefMoFEMEntity,EntityHandle,&ptrWrapperRefMoFEMFiniteElement::get_parent_ent> >,
    ordered_non_unique<
      tag<EntType_mi_tag>, const_mem_fun<ptrWrapperRefMoFEMFiniteElement::interface_type_RefMoFEMEntity,EntityType,&ptrWrapperRefMoFEMFiniteElement::get_ent_type> >,
    ordered_non_unique<
      tag<Composite_mi_tag>,
      composite_key<
	ptrWrapperRefMoFEMFiniteElement,
	const_mem_fun<ptrWrapperRefMoFEMFiniteElement::interface_type_RefMoFEMEntity,EntityHandle,&ptrWrapperRefMoFEMFiniteElement::get_parent_ent>,
	const_mem_fun<ptrWrapperRefMoFEMFiniteElement::interface_type_RefMoFEMFiniteElement,int,&ptrWrapperRefMoFEMFiniteElement::get_BitRefEdges_ulong> > >
  > > RefMoFEMFiniteElement_multiIndex;

typedef multi_index_container<
  const MoFEMEntity*,
  indexed_by<
    hashed_non_unique<
      tag<MoABEnt_mi_tag>, const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::get_ent> >
  > > MoFEMEntity_multiIndex_ent_view;

struct DofMoFEMEntity_active_change {
  bool active;
  DofMoFEMEntity_active_change(bool _active);
  void operator()(DofMoFEMEntity &_dof_);
};

typedef multi_index_container<
  const DofMoFEMEntity*,
  indexed_by<
    ordered_non_unique< 
      const_mem_fun<DofMoFEMEntity,int,&DofMoFEMEntity::get_active> >
  > > DofMoFEMEntity_multiIndex_active_view;

typedef multi_index_container<
  const DofMoFEMEntity*,
  indexed_by<
    ordered_non_unique< 
      const_mem_fun<DofMoFEMEntity,ApproximationOrder,&DofMoFEMEntity::get_dof_order> >
  > > DofMoFEMEntity_multiIndex_order_view;

typedef multi_index_container<
  const DofMoFEMEntity*,
  indexed_by<
    ordered_non_unique<
      const_mem_fun<DofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&DofMoFEMEntity::get_ent_type> >
  > > DofMoFEMEntity_multiIndex_ent_type_view;

struct NumeredDofMoFEMEntity_part_change {
  unsigned int part;
  DofIdx petsc_gloabl_dof_idx;
  NumeredDofMoFEMEntity_part_change(const unsigned int _part,const DofIdx _petsc_gloabl_dof_idx): 
    part(_part),
    petsc_gloabl_dof_idx(_petsc_gloabl_dof_idx) {};
  void operator()(NumeredDofMoFEMEntity &dof) { 
    dof.part = part;
    dof.petsc_gloabl_dof_idx = petsc_gloabl_dof_idx; 
    dof.petsc_local_dof_idx = -1;
  }
};

struct NumeredDofMoFEMEntity_local_idx_change {
  DofIdx petsc_local_dof_idx;
  NumeredDofMoFEMEntity_local_idx_change(const DofIdx _petsc_local_dof_idx): 
    petsc_local_dof_idx(_petsc_local_dof_idx) {};
  void operator()(NumeredDofMoFEMEntity &dof) { 
    dof.petsc_local_dof_idx = petsc_local_dof_idx; 
  }
};

typedef multi_index_container<
  const NumeredDofMoFEMEntity*,
  indexed_by<
    ordered_unique< 
      member<NumeredDofMoFEMEntity,const DofIdx,&NumeredDofMoFEMEntity::petsc_gloabl_dof_idx> > 
 > > NumeredDofMoFEMEntity_multiIndex_global_index_view;

struct MoFEMFE_col_change_bit_add {
  BitFieldId f_id_col;
  MoFEMFE_col_change_bit_add(const BitFieldId _f_id_col): f_id_col(_f_id_col) {};
  void operator()(MoFEMFE &MoFEMFE);
};
struct MoFEMFE_row_change_bit_add {
  BitFieldId f_id_row;
  MoFEMFE_row_change_bit_add(const BitFieldId _f_id_row): f_id_row(_f_id_row) {};
  void operator()(MoFEMFE &MoFEMFE);
};
struct EntMoFEMFE_change_bit_add {
  BitFieldId f_id_data;
  EntMoFEMFE_change_bit_add(const BitFieldId _f_id_data): f_id_data(_f_id_data) {};
  void operator()(MoFEMFE &MoFEMFE);
};

/// set uids for finite elements dofs in rows
struct EntMoFEMFE_row_dofs_change {
  Interface &moab;
  const DofMoFEMEntity_multiIndex_uid_view &uids_view;
  EntMoFEMFE_row_dofs_change(Interface &_moab,const DofMoFEMEntity_multiIndex_uid_view &_uids_view): 
    moab(_moab),uids_view(_uids_view) {};
  void operator()(EntMoFEMFE &MoFEMFE);
};

/// set uids for finite elements dofs in rows
struct EntMoFEMFE_col_dofs_change {
  Interface &moab;
  const DofMoFEMEntity_multiIndex_uid_view &uids_view;
  EntMoFEMFE_col_dofs_change(Interface &_moab,const DofMoFEMEntity_multiIndex_uid_view &_uids_view): 
    moab(_moab),uids_view(_uids_view) {};
  void operator()(EntMoFEMFE &MoFEMFE);
};

/// set uids for finite elements dofs need to calulate element matrices and vectors
struct EntMoFEMFE_data_dofs_change {
  Interface &moab;
  const DofMoFEMEntity_multiIndex_uid_view &uids_view;
  EntMoFEMFE_data_dofs_change(Interface &_moab,const DofMoFEMEntity_multiIndex_uid_view &_uids_view):
    moab(_moab),uids_view(_uids_view) {};
  void operator()(EntMoFEMFE &MoFEMFE);
};

/**
  * Internal _NumeredMoFEMFE_
  *
  */
struct _NumeredMoFEMFE_: public NumeredMoFEMFE {
  FENumeredDofMoFEMEntity_multiIndex rows_dofs;
  FENumeredDofMoFEMEntity_multiIndex cols_dofs;
  _NumeredMoFEMFE_(const EntMoFEMFE *_EntMoFEMFE_ptr): NumeredMoFEMFE(_EntMoFEMFE_ptr) {};
};

typedef multi_index_container<
  _NumeredMoFEMFE_,
  indexed_by<
    hashed_unique<
      tag<Composite_unique_mi_tag>,       
      composite_key<
	_NumeredMoFEMFE_,
	const_mem_fun<NumeredMoFEMFE::interface_type_MoFEMFE,EntityHandle,&NumeredMoFEMFE::get_meshset>,
	const_mem_fun<NumeredMoFEMFE::interface_type_EntMoFEMFE,EntityHandle,&NumeredMoFEMFE::get_ent> > >,
    ordered_non_unique<
      tag<MoFEMFE_name_mi_tag>, const_mem_fun<NumeredMoFEMFE::interface_type_MoFEMFE,string,&NumeredMoFEMFE::get_name> >,
    ordered_non_unique<
      tag<MoFEMFE_Part_mi_tag>, member<NumeredMoFEMFE,unsigned int,&NumeredMoFEMFE::part> >,
    ordered_non_unique<
      tag<Composite_mi_tag>,       
      composite_key<
	_NumeredMoFEMFE_,
	const_mem_fun<NumeredMoFEMFE::interface_type_MoFEMFE,string,&NumeredMoFEMFE::get_name>,
	member<NumeredMoFEMFE,unsigned int,&NumeredMoFEMFE::part> > >
  > > NumeredMoFEMFE_multiIndex;

struct NumeredMoFEMFE_change_part {
  unsigned int part;
  NumeredMoFEMFE_change_part(unsigned int _part): part(_part) {};
  void operator()(NumeredMoFEMFE &MoFEMFE) {
    MoFEMFE.part = part;
  }
};

struct MoFEMAdjacencies_change_by_what {
  by_what by;
  MoFEMAdjacencies_change_by_what(const by_what _by): by(_by) {}
  void operator()(MoFEMAdjacencies &e) { e.by_other |= by; }
};

/** 
  * _MoFEMProblem_ hidden from user
  */
struct _MoFEMProblem_: public MoFEMProblem {
  NumeredMoFEMFE_multiIndex numered_finite_elements;
  _MoFEMProblem_(Interface &moab,const EntityHandle _meshset): MoFEMProblem(moab,_meshset) {};
};

typedef multi_index_container<
  _MoFEMProblem_,
  indexed_by<
    ordered_unique<
      tag<Meshset_mi_tag>, member<_MoFEMProblem_::MoFEMProblem,EntityHandle,&MoFEMProblem::meshset> >,
    hashed_unique<
      tag<BitProblemId_mi_tag>, const_mem_fun<_MoFEMProblem_::MoFEMProblem,BitProblemId,&MoFEMProblem::get_id>, hashbit<BitProblemId>, eqbit<BitFieldId> >,
    hashed_unique<
      tag<MoFEMProblem_mi_tag>, const_mem_fun<_MoFEMProblem_::MoFEMProblem,string,&MoFEMProblem::get_name> >
  > > MoFEMProblem_multiIndex;

/// \brief add ref level to problem
struct problem_change_ref_level_bit_add {
  BitRefLevel bit;
  problem_change_ref_level_bit_add(const BitRefLevel _bit): bit(_bit) {};
  void operator()(MoFEMProblem &p) { *(p.tag_BitRefLevel) |= bit; };
};
/// \brief add finite element to problem
struct problem_MoFEMFE_change_bit_add {
  BitFEId f_id;
  problem_MoFEMFE_change_bit_add(const BitFEId _f_id): f_id(_f_id) {};
  void operator()(MoFEMProblem &p);
};
/// \brief increase nb. dof in row
struct problem_row_change {
  const DofMoFEMEntity *dof_ptr;
  problem_row_change(const DofMoFEMEntity *_dof_ptr);
  void operator()(_MoFEMProblem_ &e);
};
/// \brief increase nb. dof in col
struct problem_col_change {
  const DofMoFEMEntity *dof_ptr;
  problem_col_change(const DofMoFEMEntity *_dof_ptr);
  void operator()(_MoFEMProblem_ &e);
};
/// \brief zero nb. of dofs in row
struct problem_zero_nb_rows_change {
  void operator()(_MoFEMProblem_ &e);
};
/// \brief zero nb. of dofs in col
struct problem_zero_nb_cols_change {
  void operator()(_MoFEMProblem_ &e);
};

template<typename Tag> 
void get_vector_by_multi_index_tag(vector<DofMoFEMEntity> &vec_dof,const DofMoFEMEntity_multiIndex &dofs,Tag* = 0) {
  const typename boost::multi_index::index<DofMoFEMEntity_multiIndex,Tag>::type& i = get<Tag>(dofs);
  vec_dof.insert(vec_dof.end(),i.begin(),i.end());
}

template <typename T,typename V>
PetscErrorCode get_MoFEMFE_dof_uid_view(
  const T &dofs_moabfield,V &dofs_view,
  const int operation_type,const void* tag_data,const int tag_size) {
  PetscFunctionBegin;
  typedef typename boost::multi_index::index<T,Unique_mi_tag>::type dofs_by_uid;
  typedef typename boost::multi_index::index<T,Unique_mi_tag>::type::value_type value_type;
  const dofs_by_uid &dofs = dofs_moabfield.get<Unique_mi_tag>();
  const UId *uids = (UId*)tag_data;
  int size = tag_size/sizeof(UId);
  vector<const value_type*> vec;
  int ii = 0;
  for(;ii<size;ii++) {
    UId uid = uids[ii];
    typename dofs_by_uid::iterator miit = dofs.find(uid);
    if(miit==dofs.end()) continue;
    vec.push_back(&*miit);
  }
  if(operation_type==Interface::UNION) {
    dofs_view.insert(vec.begin(),vec.end());
  } else {
    //FIXME not implemented
    assert(0);
  }
  PetscFunctionReturn(0);
}

/**
 * \brief test if MoFEM is comapatibele with linked version of moab
 *
 * test EntityID <br>
 * test EntityType <br>
 */
PetscErrorCode test_moab(Interface &moab,const EntityHandle ent);

/**
 * \briMoFEMFE gives connectivity if all the edges are refined
 *
 * \param moab intefcae
 * \param conn refined tet connectivity
 * \param edge_new_nodes new nodes by edges
 * \param new_tets_conn return new egdge coonectivity
 */
void tet_type_6(Interface& moab,const EntityHandle *conn,const EntityHandle *edge_new_nodes,EntityHandle *new_tets_conn);
/**
 * \briMoFEMFE gives connectivity if 5 out of 6 egses are refined
 *
 * \param moab intefcae
 * \param conn refined tet connectivity
 * \param edge_new_nodes new nodes by edges
 * \param new_tets_conn return new egdge coonectivity
 * \return sub refinment type
 */
int tet_type_5(Interface& moab,const EntityHandle *conn,const EntityHandle *edge_new_nodes,EntityHandle *new_tets_conn);
/**
* \briMoFEMFE gives connectivity if 4 out of 6 egses are refined
*
* \param moab intefcae
* \param conn refined tet connectivity
* \param edge_new_nodes new nodes by edges
* \param new_tets_conn return new egdge coonectivity
* \return sub refinment type
*/
int tet_type_4(const EntityHandle *conn,const int *split_edges,const EntityHandle *edge_new_nodes,EntityHandle *new_tets_conn);
/**
* \briMoFEMFE gives connectivity if 3 out of 6 egses are refined
*
* \param moab intefcae
* \param conn refined tet connectivity
* \param edge_new_nodes new nodes by edges
* \param new_tets_conn return new egdge coonectivity
* \return sub refinment type
*/
int tet_type_3(const EntityHandle *conn,const int *split_edges,const EntityHandle *edge_new_nodes,EntityHandle *new_tets_conn);
/**
* \briMoFEMFE gives connectivity if 2 out of 6 egses are refined
*
* \param moab intefcae
* \param conn refined tet connectivity
* \param edge_new_nodes new nodes by edges
* \param new_tets_conn return new egdge coonectivity
* \return sub refinment type
*/
int tet_type_2(const EntityHandle *conn,const int *split_edges,const EntityHandle *edge_new_nodes,EntityHandle *new_tets_conn);
/**
* \briMoFEMFE gives connectivity if 1 out of 6 egses are refined
*
* \param moab intefcae
* \param conn refined tet connectivity
* \param edge_new_nodes new nodes by edges
* \param new_tets_conn return new egdge coonectivity
* \return sub refinment type
*/
void tet_type_1(const EntityHandle *conn,const int split_edge,const EntityHandle edge_new_node,EntityHandle *new_tets_conn);

/// only two edges are split
PetscErrorCode prism_type_1(const EntityHandle *conn,const BitRefEdges split_edges,const EntityHandle *edge_new_node,EntityHandle *new_prism_conn);
/// only tow edges are split
PetscErrorCode prism_type_2(const EntityHandle *conn,const BitRefEdges split_edges,const EntityHandle *edge_new_node,EntityHandle *new_prism_conn);
/// all edges are split
PetscErrorCode prism_type_3(const EntityHandle *conn,const BitRefEdges split_edges,const EntityHandle *edge_new_node,EntityHandle *new_prism_conn);

}

#endif //__DATASTRUCTURES_HPP__

