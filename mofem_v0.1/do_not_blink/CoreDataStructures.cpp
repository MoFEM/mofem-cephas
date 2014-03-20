/** \file CoreDataStructures.cpp
 * \brief Myltindex containes, data structures and other low-level functions 
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


#include <CoreDataStructures.hpp>

namespace MoFEM {

const bool Idx_mi_tag::IamNotPartitioned = true;
const bool Part_mi_tag::IamNotPartitioned = false;

//fields 
MoFEMField::MoFEMField(Interface &moab,const EntityHandle _meshset): meshset(_meshset),
  tag_id_data(NULL),tag_space_data(NULL),tag_rank_data(NULL),tag_name_data(NULL),tag_name_size(0) { 
  //Change those tags only by modyfiers
  ErrorCode rval;
  //id
  Tag th_FieldId;
  rval = moab.tag_get_handle("_FieldId",th_FieldId); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_FieldId,&meshset,1,(const void **)&tag_id_data); CHKERR_THROW(rval);
  //space
  Tag th_FieldSpace;
  rval = moab.tag_get_handle("_FieldSpace",th_FieldSpace); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_FieldSpace,&meshset,1,(const void **)&tag_space_data); CHKERR_THROW(rval);
  //name
  Tag th_FieldName;
  rval = moab.tag_get_handle("_FieldName",th_FieldName); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_FieldName,&meshset,1,(const void **)&tag_name_data,&tag_name_size); CHKERR_THROW(rval);
  //name prefix
  Tag th_FieldName_DataNamePrefix;
  rval = moab.tag_get_handle("_FieldName_DataNamePrefix",th_FieldName_DataNamePrefix); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_FieldName_DataNamePrefix,&meshset,1,(const void **)&tag_name_prefix_data,&tag_name_prefix_size); CHKERR_THROW(rval);
  string name_data_prefix((char *)tag_name_prefix_data,tag_name_prefix_size);
  //data
  string Tag_data_name = name_data_prefix+get_name();
  rval = moab.tag_get_handle(Tag_data_name.c_str(),th_FieldData); CHKERR_THROW(rval);
  //order
  string Tag_ApproximationOrder_name = "_App_Order_"+get_name();
  rval = moab.tag_get_handle(Tag_ApproximationOrder_name.c_str(),th_AppOrder); CHKERR_THROW(rval);
  //dof order
  string Tag_dof_ApproximationOrder_name = "_App_Dof_Order"+get_name();
  rval = moab.tag_get_handle(Tag_dof_ApproximationOrder_name.c_str(),th_AppDofOrder); CHKERR_THROW(rval);
  //rank
  Tag th_Rank;
  string Tag_rank_name = "_Field_Rank_"+get_name();
  rval = moab.tag_get_handle(Tag_rank_name.c_str(),th_Rank); CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_Rank,&meshset,1,(const void **)&tag_rank_data); CHKERR_THROW(rval);
  //dof rank
  string Tag_dof_rank_name = "_Field_Dof_Rank_"+get_name();
  rval = moab.tag_get_handle(Tag_dof_rank_name.c_str(),th_DofRank); CHKERR_THROW(rval);
  switch (*tag_space_data) {
    case H1:
      forder_entityset = NULL;
      forder_vertex = fNBVERTEX_H1;
      forder_edge = fNBEDGE_H1;
      forder_face = fNBFACE_H1;
      forder_elem = fNBVOLUME_H1;
      break;
    case Hdiv:
      forder_entityset = NULL;
      forder_vertex = fNBVERTEX_Hdiv;
      forder_edge = fNBEDGE_Hdiv;
      forder_face = fNBFACE_Hdiv;
      forder_elem = fNBVOLUME_Hdiv;
      break;
    case Hcurl:
      forder_entityset = NULL;
      forder_vertex = fNBVERTEX_Hcurl;
      forder_edge = fNBEDGE_Hcurl;
      forder_face = fNBFACE_Hcurl;
      forder_elem = fNBVOLUME_Hcurl;
      break;
    case L2:
      forder_entityset = NULL;
      forder_vertex = fNBVERTEX_L2;
      forder_edge = fNBEDGE_L2;
      forder_face = fNBFACE_L2;
      forder_elem = fNBVOLUME_L2;
      break;
    case NoField:
      forder_entityset = fNBENTITYSET_nofield;
      forder_vertex = NULL;
      forder_edge = NULL;
      forder_face = NULL;
      forder_elem = NULL;
      break;
    default:
      THROW_AT_LINE("not implemented");
  }
}
ostream& operator<<(ostream& os,const MoFEMField& e) {
  os << "name "<<e.get_name_ref()<<" BitFieldId "<< e.get_id().to_ulong() << " bit number " << e.get_bit_number() 
    << " space " << e.get_space() << " rank " << e.get_max_rank() << " meshset " << e.meshset;
  return os;
}

//moab dof
DofMoFEMEntity::DofMoFEMEntity(const MoFEMEntity *_MoFEMEntity_ptr,const ApproximationOrder _dof_order,const ApproximationRank _dof_rank,const DofIdx _dof): 
    interface_MoFEMEntity<MoFEMEntity>(_MoFEMEntity_ptr), dof(_dof),active(false) {
  if(field_ptr->tag_dof_order_data==NULL) {
    ostringstream ss;
    ss << "at " << __LINE__ << " in " << __FILE__;
    ss << " field_ptr->tag_dof_order_data==NULL";
    ss << " (top tip: check if order set to vertices is 1)";
    throw(ss.str().c_str());
  }
  assert(field_ptr->tag_dof_rank_data!=NULL);
  ((ApproximationOrder*)field_ptr->tag_dof_order_data)[dof] = _dof_order;
  ((ApproximationRank*)field_ptr->tag_dof_rank_data)[dof] = _dof_rank;
  uid = get_unique_id_calculate();
}
ostream& operator<<(ostream& os,const DofMoFEMEntity& e) {
  os << "dof_uid " << e.get_unique_id() 
    << " dof_order " << e.get_dof_order()
    << " dof_rank " << e.get_dof_rank()
    << " dof " << e.dof
    << " active " << e.active 
    << " " << *(e.field_ptr)
    << " Data " << e.get_FieldData();
  return os;
}
DofMoFEMEntity_active_change::DofMoFEMEntity_active_change(bool _active): active(_active) {}
void DofMoFEMEntity_active_change::operator()(DofMoFEMEntity &_dof_) {
  _dof_.active = active;
  assert((_dof_.get_dof_order()<=_dof_.get_max_order()));
}

//numered dof
NumeredDofMoFEMEntity::NumeredDofMoFEMEntity(const DofIdx idx,const DofMoFEMEntity* _DofMoFEMEntity_ptr): 
    interface_DofMoFEMEntity<DofMoFEMEntity>(_DofMoFEMEntity_ptr),
    dof_idx(idx),petsc_gloabl_dof_idx(-1),petsc_local_dof_idx(-1),part(-1) {}
ostream& operator<<(ostream& os,const NumeredDofMoFEMEntity& e) {
  os << "idx " << e.dof_idx << " part " << e.part 
    << " petsc idx " << e.petsc_gloabl_dof_idx 
    << " ( " << e.petsc_local_dof_idx <<  " ) "
    << *e.field_ptr;
  return os;
}

FEDofMoFEMEntity::FEDofMoFEMEntity(
  SideNumber *_side_number_ptr,
  const DofMoFEMEntity *_DofMoFEMEntity_ptr): 
  BaseFEDofMoFEMEntity(_side_number_ptr), interface_DofMoFEMEntity<DofMoFEMEntity>(_DofMoFEMEntity_ptr) {}
ostream& operator<<(ostream& os,const FEDofMoFEMEntity& e) {
  os << "local dof MoFEMFiniteElement idx " 
    << "side_number " << e.side_number_ptr->side_number << " "
    << "sense " << e.side_number_ptr->sense << " "
    << *e.field_ptr;
  return os;
}

FENumeredDofMoFEMEntity::FENumeredDofMoFEMEntity(
  SideNumber *_side_number_ptr,
  const NumeredDofMoFEMEntity *_NumeredDofMoFEMEntity_ptr): 
  BaseFEDofMoFEMEntity(_side_number_ptr), interface_NumeredDofMoFEMEntity<NumeredDofMoFEMEntity>(_NumeredDofMoFEMEntity_ptr) {}
ostream& operator<<(ostream& os,const FENumeredDofMoFEMEntity& e) {
  os << "local dof MoFEMFiniteElement idx " 
    << "side_number " << e.side_number_ptr->side_number << " "
    << "sense " << e.side_number_ptr->sense << " "
    << *e.field_ptr;
  return os;
}

//MoFEMEntityEntMoFEMFiniteElementAdjacencyMap
MoFEMEntityEntMoFEMFiniteElementAdjacencyMap::MoFEMEntityEntMoFEMFiniteElementAdjacencyMap(const MoFEMEntity *_MoFEMEntity_ptr,const EntMoFEMFiniteElement *_EntMoFEMFiniteElement_ptr):
  by_other(0),MoFEMEntity_ptr(_MoFEMEntity_ptr),EntMoFEMFiniteElement_ptr(_EntMoFEMFiniteElement_ptr) {}
ostream& operator<<(ostream& os,const MoFEMEntityEntMoFEMFiniteElementAdjacencyMap& e) {
  os << "by_other " << bitset<3>(e.by_other) << " "
    << *e.MoFEMEntity_ptr << endl << *e.EntMoFEMFiniteElement_ptr->fe_ptr;
  return os;
}

//....
PetscErrorCode test_moab(Interface &moab,const EntityHandle ent) {
  PetscFunctionBegin;
  //tets type
  EntityType type = (EntityType)((ent&MB_TYPE_MASK)>>MB_ID_WIDTH);
  if(type != moab.type_from_handle(ent)) SETERRQ(PETSC_COMM_SELF,1,"incosistencies with type_from_handle");
  //tets id
  EntityID id = (EntityType)(ent&MB_ID_MASK);
  if(id != moab.id_from_handle(ent)) SETERRQ(PETSC_COMM_SELF,1,"incosistencies with id_from_handle");
  PetscFunctionReturn(0);
}

}
