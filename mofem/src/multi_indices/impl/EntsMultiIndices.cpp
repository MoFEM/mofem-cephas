/** \file CoreDataStructures.cpp
 * \brief Multi-index containers, data structures and other low-level functions
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

#include <Includes.hpp>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <MaterialBlocks.hpp>
#include <CubitBCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
// #include <ProblemsMultiIndices.hpp>
// #include <AdjacencyMultiIndices.hpp>
// #include <BCMultiIndices.hpp>
// #include <CoreDataStructures.hpp>

namespace MoFEM {

//basic moab ent
BasicMoFEMEntity::BasicMoFEMEntity(Interface &moab,const EntityHandle _ent):
  ent(_ent),sharing_procs_ptr(NULL),sharing_handlers_ptr(NULL) {
  switch (get_ent_type()) {
    case MBVERTEX:
    case MBEDGE:
    case MBTRI:
    case MBQUAD:
    case MBTET:
    case MBPRISM:
    case MBENTITYSET:
    break;
    default:
      THROW_MESSAGE("this entity type is currently not implemented");
  }
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  ErrorCode rval;
  rval = pcomm->get_owner_handle(ent,owner_proc,moab_owner_handle); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(pcomm->pstatus_tag(),&ent,1,(const void **)&pstatus_val_ptr); CHKERR_MOAB(rval);

  if(*pstatus_val_ptr & PSTATUS_MULTISHARED) {
    // entity is multi shared
    rval = moab.tag_get_by_ptr(pcomm->sharedps_tag(),&ent,1,(const void **)&sharing_procs_ptr); CHKERR_MOAB(rval);
    rval = moab.tag_get_by_ptr(pcomm->sharedhs_tag(),&ent,1,(const void **)&sharing_handlers_ptr); CHKERR_MOAB(rval);
  } else if(*pstatus_val_ptr & PSTATUS_SHARED) {
    // shared
    rval = moab.tag_get_by_ptr(pcomm->sharedp_tag(),&ent,1,(const void **)&sharing_procs_ptr); CHKERR_MOAB(rval);
    rval = moab.tag_get_by_ptr(pcomm->sharedh_tag(),&ent,1,(const void **)&sharing_handlers_ptr); CHKERR_MOAB(rval);
  }

}
//ref moab ent
BitRefEdges MoFEM::RefMoFEMElement::DummyBitRefEdges = BitRefEdges(0);
RefMoFEMEntity::RefMoFEMEntity(
  Interface &moab,const EntityHandle _ent):
    BasicMoFEMEntity(moab,_ent),tag_parent_ent(NULL),tag_BitRefLevel(NULL) {
  ErrorCode rval;
  Tag th_RefParentHandle,th_RefBitLevel;
  rval = moab.tag_get_handle("_RefParentHandle",th_RefParentHandle); MOAB_THROW(rval);
  rval = moab.tag_get_handle("_RefBitLevel",th_RefBitLevel); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefParentHandle,&ent,1,(const void **)&tag_parent_ent); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitLevel,&ent,1,(const void **)&tag_BitRefLevel); MOAB_THROW(rval);
}
ostream& operator<<(ostream& os,const RefMoFEMEntity& e) {
  os << "ent " << e.ent;
  os << " pstatus "<< bitset<8>(e.get_pstatus());
  os << " owner ent " << e.get_owner_ent();
  os << " owner proc " << e.get_owner_proc();
  os << " parent ent " << e.get_parent_ent();
  //os << " BitRefLevel " << e.get_BitRefLevel();
  os << " ent type " << e.get_ent_type();
  os << " ent parent type " << e.get_parent_ent_type();
  return os;
}

//moab ent
MoFEMEntity::MoFEMEntity(Interface &moab,const MoFEMField *_field_ptr,const RefMoFEMEntity *_ref_ent_ptr):
  interface_MoFEMField<MoFEMField>(_field_ptr),interface_RefMoFEMEntity<RefMoFEMEntity>(_ref_ent_ptr),ref_mab_ent_ptr(_ref_ent_ptr),
  tag_order_data(NULL),tag_FieldData(NULL),tag_FieldData_size(0),tag_dof_order_data(NULL),tag_dof_rank_data(NULL) {
  ErrorCode rval;
  EntityHandle ent = get_ent();
  rval = moab.tag_get_by_ptr(field_ptr->th_AppOrder,&ent,1,(const void **)&tag_order_data); MOAB_THROW(rval);
  local_uid = get_local_unique_id_calculate();
  global_uid = get_global_unique_id_calculate();
  rval = moab.tag_get_by_ptr(field_ptr->th_FieldData,&ent,1,(const void **)&tag_FieldData,&tag_FieldData_size);
  if(rval == MB_SUCCESS) {
    if( (unsigned int)tag_FieldData_size != 0 ) {
      int tag_size[1];
      rval = moab.tag_get_by_ptr(field_ptr->th_AppDofOrder,&ent,1,(const void **)&tag_dof_order_data,tag_size); MOAB_THROW(rval);
      assert(tag_size[0]/sizeof(ApproximationOrder) == tag_FieldData_size/sizeof(FieldData));
      rval = moab.tag_get_by_ptr(field_ptr->th_DofRank,&ent,1,(const void **)&tag_dof_rank_data,tag_size); MOAB_THROW(rval);
      assert(tag_size[0]/sizeof(FieldCoefficientsNumber) == tag_FieldData_size/sizeof(FieldData));
    }
  }
}
MoFEMEntity::~MoFEMEntity() {}
ostream& operator<<(ostream& os,const MoFEMEntity& e) {
  os << "ent_global_uid " << (UId)e.get_global_unique_id()
    << " ent_local_uid " << (UId)e.get_local_unique_id()
    << " entity "<< e.get_ent() << " type " << e.get_ent_type()
    << " pstatus "<< bitset<8>(e.get_pstatus()) << " owner handle " << e.get_owner_ent() << " owner proc " << e.get_owner_proc()
    << " order "<<e.get_max_order()<<" "<<*e.field_ptr;
  return os;
}
void MoFEMEntity_change_order::operator()(MoFEMEntity &e) {
  ErrorCode rval;
  int nb_dofs = e.get_order_nb_dofs(order)*e.get_nb_of_coeffs();
  ApproximationOrder& ent_order = *((ApproximationOrder*)e.tag_order_data);
  ent_order = order;
  EntityHandle ent = e.get_ent();
  rval = moab.tag_get_by_ptr(e.field_ptr->th_FieldData,&ent,1,(const void **)&e.tag_FieldData,&e.tag_FieldData_size);
  if(rval == MB_SUCCESS) {
    //data
    if( nb_dofs*sizeof(FieldData) == (unsigned int)e.tag_FieldData_size ) {
      rval = moab.tag_get_by_ptr(e.field_ptr->th_FieldData,&ent,1,(const void **)&e.tag_FieldData,&e.tag_FieldData_size); CHKERR_MOAB(rval);
      int tag_size[1];
      rval = moab.tag_get_by_ptr(e.field_ptr->th_AppDofOrder,&ent,1,(const void **)&e.tag_dof_order_data,tag_size); CHKERR_MOAB(rval);
      assert(tag_size[0]/sizeof(ApproximationOrder) == e.tag_FieldData_size/sizeof(FieldData));
      rval = moab.tag_get_by_ptr(e.field_ptr->th_DofRank,&ent,1,(const void **)&e.tag_dof_rank_data,tag_size); CHKERR_MOAB(rval);
      assert(tag_size[0]/sizeof(FieldCoefficientsNumber) == e.tag_FieldData_size/sizeof(FieldData));
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
    assert(tag_size[0]/sizeof(FieldCoefficientsNumber) == e.tag_FieldData_size/sizeof(FieldData));
    data_dof_rank.resize(e.tag_FieldData_size/sizeof(FieldData));
    FieldCoefficientsNumber *ptr_dof_rank_begin = (FieldCoefficientsNumber*)e.tag_dof_rank_data;
    FieldCoefficientsNumber *ptr_dof_rank_end = (FieldCoefficientsNumber*)e.tag_dof_rank_data + e.tag_FieldData_size/sizeof(FieldData);
    copy(ptr_dof_rank_begin,ptr_dof_rank_end,data_dof_rank.begin());
  }
  if(nb_dofs>0) {
    data.resize(nb_dofs,0);
    int tag_size[1];
    tag_size[0] = data.size()*sizeof(FieldData);
    void const* tag_data[] = { &data[0] };
    rval = moab.tag_set_by_ptr(e.field_ptr->th_FieldData,&ent,1,tag_data,tag_size); CHKERR_MOAB(rval); assert(rval==MB_SUCCESS);
    rval = moab.tag_get_by_ptr(e.field_ptr->th_FieldData,&ent,1,(const void **)&e.tag_FieldData,&e.tag_FieldData_size); CHKERR_MOAB(rval);
    //order
    data_dof_order.resize(nb_dofs,0);
    tag_size[0] = data_dof_order.size()*sizeof(ApproximationOrder);
    tag_data[0] = &data_dof_order[0];
    rval = moab.tag_set_by_ptr(e.field_ptr->th_AppDofOrder,&ent,1,tag_data,tag_size); CHKERR_MOAB(rval); assert(rval==MB_SUCCESS);
    rval = moab.tag_get_by_ptr(e.field_ptr->th_AppDofOrder,&ent,1,(const void **)&e.tag_dof_order_data,tag_size); CHKERR_MOAB(rval);
    //rank
    data_dof_rank.resize(nb_dofs,0);
    tag_size[0] = data_dof_rank.size()*sizeof(FieldCoefficientsNumber);
    tag_data[0] = &data_dof_rank[0];
    rval = moab.tag_set_by_ptr(e.field_ptr->th_DofRank,&ent,1,tag_data,tag_size); CHKERR_MOAB(rval); assert(rval==MB_SUCCESS);
    rval = moab.tag_get_by_ptr(e.field_ptr->th_DofRank,&ent,1,(const void **)&e.tag_dof_rank_data,tag_size); CHKERR_MOAB(rval);
  }
}

}
