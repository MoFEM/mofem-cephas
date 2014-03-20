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

//basic moab ent
BasicMoFEMEntity::BasicMoFEMEntity(const EntityHandle _ent): ent(_ent) {
  switch (get_ent_type()) {
    case MBVERTEX:
    case MBEDGE:
    case MBTRI:
    case MBTET:
    case MBPRISM:
    case MBENTITYSET:
    break;
    default:
      THROW_AT_LINE("this entity type is currently not implemented");
  }
}
//ref moab ent
BitRefEdges MoFEM::RefMoFEMElement::DummyBitRefEdges = BitRefEdges(0);
RefMoFEMEntity::RefMoFEMEntity(
  Interface &moab,const EntityHandle _ent): 
    BasicMoFEMEntity(_ent),tag_parent_ent(NULL),tag_BitRefLevel(NULL) {
  ErrorCode rval;
  Tag th_RefParentHandle,th_RefBitLevel;
  rval = moab.tag_get_handle("_RefParentHandle",th_RefParentHandle); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_RefBitLevel",th_RefBitLevel); CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefParentHandle,&ent,1,(const void **)&tag_parent_ent); CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitLevel,&ent,1,(const void **)&tag_BitRefLevel); CHKERR_THROW(rval);
}
ostream& operator<<(ostream& os,const RefMoFEMEntity& e) {
  os << "ent " << e.ent;
  os << " parent ent " << e.get_parent_ent();
  os << " BitRefLevel " << e.get_BitRefLevel();
  os << " ent type " << e.get_ent_type();
  os << " ent parent type " << e.get_parent_ent_type();
  return os;
}

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

//moab ent
MoFEMEntity::MoFEMEntity(Interface &moab,const MoFEMField *_FieldData,const RefMoFEMEntity *_ref_ent_ptr): 
  interface_MoFEMField<MoFEMField>(_FieldData),interface_RefMoFEMEntity<RefMoFEMEntity>(_ref_ent_ptr),ref_mab_ent_ptr(_ref_ent_ptr),
  tag_order_data(NULL),tag_FieldData(NULL),tag_FieldData_size(0),tag_dof_order_data(NULL),tag_dof_rank_data(NULL) {
  const EntityType type = get_ent_type(); 
  switch (type) {
    case MBVERTEX:
      forder = field_ptr->forder_vertex;
      break;
    case MBEDGE:
      forder = field_ptr->forder_edge;
      break;
    case MBTRI:
      forder = field_ptr->forder_face;
      break;
    case MBTET:
      forder = field_ptr->forder_elem;
      break;
    case MBPRISM:
      forder = field_ptr->forder_elem;
      break;
    case MBENTITYSET:
      forder = field_ptr->forder_entityset;
      break;
    default: {
      THROW_AT_LINE("data inconsistency");
    }
  }
  ErrorCode rval;
  EntityHandle ent = get_ent();
  rval = moab.tag_get_by_ptr(field_ptr->th_AppOrder,&ent,1,(const void **)&tag_order_data); CHKERR_THROW(rval);
  uid = get_unique_id_calculate();
  rval = moab.tag_get_by_ptr(field_ptr->th_FieldData,&ent,1,(const void **)&tag_FieldData,&tag_FieldData_size); 
  if(rval == MB_SUCCESS) {
    if( (unsigned int)tag_FieldData_size != 0 ) {
      int tag_size[1];
      rval = moab.tag_get_by_ptr(field_ptr->th_AppDofOrder,&ent,1,(const void **)&tag_dof_order_data,tag_size); CHKERR_THROW(rval);
      assert(tag_size[0]/sizeof(ApproximationOrder) == tag_FieldData_size/sizeof(FieldData));
      rval = moab.tag_get_by_ptr(field_ptr->th_DofRank,&ent,1,(const void **)&tag_dof_rank_data,tag_size); CHKERR_THROW(rval);
      assert(tag_size[0]/sizeof(ApproximationRank) == tag_FieldData_size/sizeof(FieldData));
    }
  }
}
ostream& operator<<(ostream& os,const MoFEMEntity& e) {
  os << "ent_uid " << e.get_unique_id() << " entity "<< e.get_ent() << " type " << e.get_ent_type()
    <<" order "<<e.get_max_order()<<" "<<*e.field_ptr;
  return os;
}
void MoFEMEntity_change_order::operator()(MoFEMEntity &e) {
  ErrorCode rval;
  int nb_dofs = e.forder(order)*e.get_max_rank();
  ApproximationOrder& ent_order = *((ApproximationOrder*)e.tag_order_data);
  ent_order = order;
  EntityHandle ent = e.get_ent();
  rval = moab.tag_get_by_ptr(e.field_ptr->th_FieldData,&ent,1,(const void **)&e.tag_FieldData,&e.tag_FieldData_size); 
  if(rval == MB_SUCCESS) {
    //data
    if( nb_dofs*sizeof(FieldData) == (unsigned int)e.tag_FieldData_size ) {
      rval = moab.tag_get_by_ptr(e.field_ptr->th_FieldData,&ent,1,(const void **)&e.tag_FieldData,&e.tag_FieldData_size); CHKERR(rval);
      int tag_size[1];
      rval = moab.tag_get_by_ptr(e.field_ptr->th_AppDofOrder,&ent,1,(const void **)&e.tag_dof_order_data,tag_size); CHKERR(rval);
      assert(tag_size[0]/sizeof(ApproximationOrder) == e.tag_FieldData_size/sizeof(FieldData));
      rval = moab.tag_get_by_ptr(e.field_ptr->th_DofRank,&ent,1,(const void **)&e.tag_dof_rank_data,tag_size); CHKERR(rval);
      assert(tag_size[0]/sizeof(ApproximationRank) == e.tag_FieldData_size/sizeof(FieldData));
      return;
    }
    data.resize(e.tag_FieldData_size/sizeof(FieldData));
    FieldData *ptr_begin = (FieldData*)e.tag_FieldData;
    FieldData *ptr_end = (FieldData*)e.tag_FieldData + e.tag_FieldData_size/sizeof(FieldData);
    copy(ptr_begin,ptr_end,data.begin());
    int tag_size[1];
    //order
    rval = moab.tag_get_by_ptr(e.field_ptr->th_AppDofOrder,&ent,1,(const void **)&e.tag_dof_order_data,tag_size); 
    assert(rval == MB_SUCCESS);
    assert(tag_size[0]/sizeof(ApproximationOrder) == e.tag_FieldData_size/sizeof(FieldData));
    data_dof_order.resize(e.tag_FieldData_size/sizeof(FieldData));
    ApproximationOrder *ptr_dof_order_begin = (ApproximationOrder*)e.tag_dof_order_data;
    ApproximationOrder *ptr_dof_order_end = (ApproximationOrder*)e.tag_dof_order_data + e.tag_FieldData_size/sizeof(FieldData);
    copy(ptr_dof_order_begin,ptr_dof_order_end,data_dof_order.begin());
    //rank
    rval = moab.tag_get_by_ptr(e.field_ptr->th_DofRank,&ent,1,(const void **)&e.tag_dof_rank_data,tag_size); 
    assert(rval == MB_SUCCESS);
    assert(tag_size[0]/sizeof(ApproximationRank) == e.tag_FieldData_size/sizeof(FieldData));
    data_dof_rank.resize(e.tag_FieldData_size/sizeof(FieldData));
    ApproximationRank *ptr_dof_rank_begin = (ApproximationRank*)e.tag_dof_rank_data;
    ApproximationRank *ptr_dof_rank_end = (ApproximationRank*)e.tag_dof_rank_data + e.tag_FieldData_size/sizeof(FieldData);
    copy(ptr_dof_rank_begin,ptr_dof_rank_end,data_dof_rank.begin());
  }
  if(nb_dofs>0) {
    data.resize(nb_dofs,0);
    int tag_size[1]; tag_size[0] = data.size()*sizeof(FieldData);
    void const* tag_data[] = { &data[0] };
    rval = moab.tag_set_by_ptr(e.field_ptr->th_FieldData,&ent,1,tag_data,tag_size); CHKERR(rval); assert(rval==MB_SUCCESS);
    rval = moab.tag_get_by_ptr(e.field_ptr->th_FieldData,&ent,1,(const void **)&e.tag_FieldData,&e.tag_FieldData_size); CHKERR(rval);
    //order
    data_dof_order.resize(nb_dofs,0);
    tag_size[0] = data_dof_order.size()*sizeof(ApproximationOrder);
    tag_data[0] = &data_dof_order[0];
    rval = moab.tag_set_by_ptr(e.field_ptr->th_AppDofOrder,&ent,1,tag_data,tag_size); CHKERR(rval); assert(rval==MB_SUCCESS);
    rval = moab.tag_get_by_ptr(e.field_ptr->th_AppDofOrder,&ent,1,(const void **)&e.tag_dof_order_data,tag_size); CHKERR(rval);
    //rank
    data_dof_rank.resize(nb_dofs,0);
    tag_size[0] = data_dof_rank.size()*sizeof(ApproximationRank);
    tag_data[0] = &data_dof_rank[0];
    rval = moab.tag_set_by_ptr(e.field_ptr->th_DofRank,&ent,1,tag_data,tag_size); CHKERR(rval); assert(rval==MB_SUCCESS);
    rval = moab.tag_get_by_ptr(e.field_ptr->th_DofRank,&ent,1,(const void **)&e.tag_dof_rank_data,tag_size); CHKERR(rval);
  }
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
