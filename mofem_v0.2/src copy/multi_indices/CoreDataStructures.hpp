/** \file CoreDataStructures.hpp
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

namespace MoFEM {

const int prism_adj_edges[] = { 6,7,8, -1,-1,-1, 0,1,2 };
const int prism_edges_conn[6][2] = { {0,1},{1,2},{2,0}, {3,4}, {4,5}, {5,3} };

inline int fNBENTITYSET_nofield(int P) { (void)P; return 1; }
//
inline int fNBVERTEX_L2(int P) { (void)P; return 0; }
inline int fNBVOLUME_L2(int P) { return NBVOLUME_L2(P); }
inline int fNBFACE_L2(int P) { return NBFACE_L2(P); }
inline int fNBEDGE_L2(int P) { return NBEDGE_L2(P); }

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
/// number of approx. functions for HCURL space on vertex
inline int fNBVERTEX_HCURL(int P) { (void)P; return 0; }
inline int fNBEDGE_HCURL(int P) { return NBEDGE_HCURL(P); }
inline int fNBFACE_HCURL(int P) { return NBFACE_HCURL(P); }
inline int fNBVOLUME_HCURL(int P) { return NBVOLUME_HCURL(P); }
//
/// \brief number of approx. functions for HDIV space on vertex
///
/// zero number of digrees of freedom on vertex for that space
inline int fNBVERTEX_HDIV(int P) { (void)P; return 0; }
/// number of approx. functions for HDIV space on edge
inline int fNBEDGE_HDIV(int P) { assert(P==P); (void)P; return NBEDGE_HDIV(P); }
/// number of approx. functions for HDIV space on face
inline int fNBFACE_HDIV(int P) { return NBFACE_HDIV(P); }
/// number of approx. functions for HDIV space on voulem
inline int fNBVOLUME_HDIV(int P) { return NBVOLUME_HDIV(P); }

struct ltstr
{ inline bool operator()(const string &s1, const string& s2) const
  { return strcmp(s1.c_str(), s2.c_str()) < 0; } };

struct CubitMeshSets_change_add_bit_to_CubitBCType {
  CubitBC_BitSet bit;
  CubitMeshSets_change_add_bit_to_CubitBCType(const CubitBC_BitSet &_bit): bit(_bit) {};
  void operator()(CubitMeshSets &e) { 
    e.CubitBCType |= bit;
  }
};

typedef multi_index_container<
  const MoFEMEntity*,
  indexed_by<
    hashed_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::get_ent> >
  > > MoFEMEntity_multiIndex_ent_view;

typedef multi_index_container<
  const NumeredDofMoFEMEntity*,
  indexed_by<
    ordered_unique< 
      member<NumeredDofMoFEMEntity,const DofIdx,&NumeredDofMoFEMEntity::petsc_gloabl_dof_idx> > 
 > > NumeredDofMoFEMEntity_multiIndex_global_index_view;

struct MoFEMFiniteElement_col_change_bit_add {
  BitFieldId f_id_col;
  MoFEMFiniteElement_col_change_bit_add(const BitFieldId _f_id_col): f_id_col(_f_id_col) {};
  void operator()(MoFEMFiniteElement &MoFEMFiniteElement);
};
struct MoFEMFiniteElement_row_change_bit_add {
  BitFieldId f_id_row;
  MoFEMFiniteElement_row_change_bit_add(const BitFieldId _f_id_row): f_id_row(_f_id_row) {};
  void operator()(MoFEMFiniteElement &MoFEMFiniteElement);
};
struct EntMoFEMFiniteElement_change_bit_add {
  BitFieldId f_id_data;
  EntMoFEMFiniteElement_change_bit_add(const BitFieldId _f_id_data): f_id_data(_f_id_data) {};
  void operator()(MoFEMFiniteElement &MoFEMFiniteElement);
};

struct MoFEMFiniteElement_col_change_bit_off {
  BitFieldId f_id_col;
  MoFEMFiniteElement_col_change_bit_off(const BitFieldId _f_id_col): f_id_col(_f_id_col) {};
  void operator()(MoFEMFiniteElement &MoFEMFiniteElement);
};
struct MoFEMFiniteElement_row_change_bit_off {
  BitFieldId f_id_row;
  MoFEMFiniteElement_row_change_bit_off(const BitFieldId _f_id_row): f_id_row(_f_id_row) {};
  void operator()(MoFEMFiniteElement &MoFEMFiniteElement);
};
struct EntMoFEMFiniteElement_change_bit_off {
  BitFieldId f_id_data;
  EntMoFEMFiniteElement_change_bit_off(const BitFieldId _f_id_data): f_id_data(_f_id_data) {};
  void operator()(MoFEMFiniteElement &MoFEMFiniteElement);
};

struct MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_ByWhat {
  ByWhat by;
  MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_ByWhat(const ByWhat _by): by(_by) {}
  void operator()(MoFEMEntityEntMoFEMFiniteElementAdjacencyMap &e) { e.by_other |= by; }
};

template<typename Tag> 
void get_vector_by_multi_index_tag(vector<DofMoFEMEntity> &vec_dof,const DofMoFEMEntity_multiIndex &dofs,Tag* = 0) {
  const typename boost::multi_index::index<DofMoFEMEntity_multiIndex,Tag>::type& i = get<Tag>(dofs);
  vec_dof.insert(vec_dof.end(),i.begin(),i.end());
}

/** \brief if moab entity hanled
  */
PetscErrorCode test_moab(Interface &moab,const EntityHandle ent);

}

#endif //__DATASTRUCTURES_HPP__

