/** \file CoreDataStructures.cpp
 * \brief Multi-index containers, data structures and other low-level functions
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

#include <Includes.hpp>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>

#define IS_BUILDING_MB
#include <moab/Error.hpp>
#include <SparseTag.hpp>
#include <DenseTag.hpp>


namespace MoFEM {

static MoABErrorCode rval;
static moab::Error error;

inline void* get_tag_ptr(SequenceManager *sequence_manager,Tag th,EntityHandle ent,int *tag_size) {
  ApproximationOrder *ret_val;
  if(th->get_storage_type()==MB_TAG_SPARSE) {
    rval = static_cast<SparseTag*>(th)->get_data(
      sequence_manager,&error,&ent,1,(const void**)&ret_val,tag_size
    ); MOAB_THROW(rval);
    return ret_val;
  } else {
    rval = static_cast<DenseTag*>(th)->get_data(
      sequence_manager,&error,&ent,1,(const void**)&ret_val,tag_size
    ); MOAB_THROW(rval);
    return ret_val;
  }
}

BasicEntityData::BasicEntityData(const moab::Interface &moab):
moab(const_cast<moab::Interface&>(moab)),
distributedMesh(true) {
  rval = moab.tag_get_handle("_RefParentHandle",th_RefParentHandle); MOAB_THROW(rval);
  rval = moab.tag_get_handle("_RefBitLevel",th_RefBitLevel); MOAB_THROW(rval);
}
BasicEntityData::~BasicEntityData() {
}

//basic moab ent
BasicEntity::BasicEntity(
  const boost::shared_ptr<BasicEntityData>& basic_data_ptr,
  const EntityHandle ent
):
basicDataPtr(basic_data_ptr),
ent(ent) {
  switch (getEntType()) {
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
  ParallelComm* pcomm = ParallelComm::get_pcomm(&basicDataPtr->moab,MYPCOMM_INDEX);
  if(pcomm == NULL) THROW_MESSAGE("pcomm is null");
  // unsigned char pstatus;
  // rval = basicDataPtr->
  // moab.tag_get_data(pcomm->pstatus_tag(),&ent,1,&pstatus); MOAB_THROW(rval);
  rval = pcomm->get_owner_handle(ent,owner_proc,moab_owner_handle); MOAB_THROW(rval);
  part_proc = owner_proc;
}

unsigned char BasicEntity::getPStatus() const {
  ParallelComm* pcomm = ParallelComm::get_pcomm(&basicDataPtr->moab,MYPCOMM_INDEX);
  return *((unsigned char*)MoFEM::get_tag_ptr(
    static_cast<moab::Core*>(&basicDataPtr->moab)->sequence_manager(),pcomm->pstatus_tag(),ent,NULL
  ));
}

//ref moab ent
BitRefEdges MoFEM::RefElement::DummyBitRefEdges = BitRefEdges(0);
RefEntity::RefEntity(
  const boost::shared_ptr<BasicEntityData>& basic_data_ptr,
  const EntityHandle ent):
BasicEntity(basic_data_ptr,ent) {
}

EntityHandle* RefEntity::getParentEntPtr() const {
  return (EntityHandle*)get_tag_ptr(
    static_cast<moab::Core*>(&basicDataPtr->moab)->sequence_manager(),basicDataPtr->th_RefParentHandle,ent,NULL
  );
}

BitRefLevel* RefEntity::getBitRefLevelPtr() const {
  return (BitRefLevel*)get_tag_ptr(
    static_cast<moab::Core*>(&basicDataPtr->moab)->sequence_manager(),basicDataPtr->th_RefBitLevel,ent,NULL
  );
}

PetscErrorCode getPatentEnt(Interface &moab,Range ents,std::vector<EntityHandle> vec_patent_ent) {
  MoABErrorCode rval;
  PetscFunctionBegin;
  Tag th_ref_parent_handle;
  rval = moab.tag_get_handle("_RefParentHandle",th_ref_parent_handle); CHKERRQ_MOAB(rval);
  vec_patent_ent.resize(ents.size());
  rval = moab.tag_get_data(th_ref_parent_handle,ents,&*vec_patent_ent.begin()); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}

PetscErrorCode RefEntity::getBitRefLevel(Interface &moab,Range ents,std::vector<BitRefLevel> vec_bit_ref_level) {
  MoABErrorCode rval;
  PetscFunctionBegin;
  Tag th_ref_bit_level;
  rval = moab.tag_get_handle("_RefBitLevel",th_ref_bit_level); MOAB_THROW(rval);
  vec_bit_ref_level.resize(ents.size());
  rval = moab.tag_get_data(th_ref_bit_level,ents,&*vec_bit_ref_level.begin()); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}

std::ostream& operator<<(std::ostream& os,const RefEntity& e) {
  os << "ent " << e.ent;
  os << " pstatus "<< std::bitset<8>(e.getPStatus());
  os << " owner ent " << e.getOwnerEnt();
  os << " owner proc " << e.getOwnerProc();
  os << " parent ent " << e.getParentEnt();
  //os << " BitRefLevel " << e.getBitRefLevel();
  os << " ent type " << e.getEntType();
  os << " ent parent type " << e.getParentEntType();
  return os;
}

//moab ent
FieldEntity::FieldEntity(
  const boost::shared_ptr<Field>& field_ptr,
  const boost::shared_ptr<RefEntity>& ref_ent_ptr
):
interface_Field<Field>(field_ptr),
interface_RefEntity<RefEntity>(ref_ent_ptr),
// tag_order_data(NULL),
tag_FieldData(NULL),
tag_FieldData_size(0) {
  //MoABErrorCode rval;
  EntityHandle ent = getEnt();
  moab::Interface &moab = ref_ent_ptr->basicDataPtr->moab;
  rval = moab.tag_get_by_ptr(field_ptr->th_FieldData,&ent,1,(const void **)&tag_FieldData,&tag_FieldData_size);
  globalUid = getGlobalUniqueIdCalculate();
  getDofOrderMap().resize(MAX_DOFS_ON_ENTITY,-1);
}

ApproximationOrder* FieldEntity::getMaxOrderPtr() {
  return (ApproximationOrder*)MoFEM::get_tag_ptr(
    static_cast<moab::Core*>(&sFieldPtr->moab)->sequence_manager(),sFieldPtr->th_AppOrder,sPtr->ent,NULL
  );
}
ApproximationOrder FieldEntity::getMaxOrder() const {
  return *(ApproximationOrder*)MoFEM::get_tag_ptr(
    static_cast<moab::Core*>(&sFieldPtr->moab)->sequence_manager(),sFieldPtr->th_AppOrder,sPtr->ent,NULL
  );

}

FieldEntity::~FieldEntity() {}
std::ostream& operator<<(std::ostream& os,const FieldEntity& e) {
  os << "ent_global_uid " << (UId)e.getGlobalUniqueId()
    // << " ent_local_uid " << (UId)e.get_local_unique_id()
    << " entity "<< e.getEnt() << " type " << e.getEntType()
    << " pstatus "<< std::bitset<8>(e.getPStatus()) << " owner handle " << e.getOwnerEnt() << " owner proc " << e.getOwnerProc()
    << " order "<<e.getMaxOrder()<<" "<< *e.sFieldPtr;
  return os;
}
void FieldEntity_change_order::operator()(FieldEntity *e) {
  MoABErrorCode rval;
  moab::Interface &moab = e->sPtr->basicDataPtr->moab;
  unsigned int nb_dofs = e->getOrderNbDofs(order)*e->getNbOfCoeffs();
  ApproximationOrder& ent_order = *(e->getMaxOrderPtr());
  ent_order = order;
  EntityHandle ent = e->getEnt();
  // Get pointer and size of field values tag
  rval = moab.tag_get_by_ptr(
    e->sFieldPtr->th_FieldData,
    &ent,1,(const void **)&e->tag_FieldData,&e->tag_FieldData_size
  );
  // Tag exist and are some data on it
  if(rval == MB_SUCCESS) {
    // Check if size of filed values tag is correct
    if(nb_dofs*sizeof(FieldData) == (unsigned int)e->tag_FieldData_size) {
      return;
    } else if(nb_dofs == 0) {
      // Delete data on this entity
      rval = moab.tag_delete_data(e->sFieldPtr->th_FieldData,&ent,1); MOAB_THROW(rval);
      e->tag_FieldData_size = 0;
      return;
    }
    // Size of tag is different than new seize, so copy data to new container
    data.resize(e->tag_FieldData_size/sizeof(FieldData));
    FieldData *ptr_begin = (FieldData*)e->tag_FieldData;
    FieldData *ptr_end = (FieldData*)e->tag_FieldData + e->tag_FieldData_size/sizeof(FieldData);
    std::copy(ptr_begin,ptr_end,data.begin());
  }
  // Set new data
  if(nb_dofs>0) {
    // Set field dof data
    data.resize(nb_dofs,0);
    int tag_size[1];
    tag_size[0] = data.size()*sizeof(FieldData);
    void const* tag_data[] = { &data[0] };
    rval = moab.tag_set_by_ptr(e->sFieldPtr->th_FieldData,&ent,1,tag_data,tag_size); MOAB_THROW(rval);
    rval = moab.tag_get_by_ptr(
      e->sFieldPtr->th_FieldData,
      &ent,1,(const void **)&e->tag_FieldData,&e->tag_FieldData_size
    ); MOAB_THROW(rval);
    if(nb_dofs != e->tag_FieldData_size/sizeof(FieldData))
    THROW_MESSAGE("Data inconsistency");
  }
}

}
