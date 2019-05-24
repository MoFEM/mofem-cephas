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

#define IS_BUILDING_MB
#include <moab/Error.hpp>

namespace MoFEM {

static moab::Error error;

inline void *get_tag_ptr(moab::Interface &moab, Tag th, EntityHandle ent,
                         int *tag_size) {
  ApproximationOrder *ret_val;
  rval = moab.tag_get_by_ptr(th, &ent, 1, (const void **)&ret_val, tag_size);
  if (rval != MB_SUCCESS) {
    *tag_size = 0;
    return NULL;
  } else {
    return ret_val;
  }
}

BasicEntityData::BasicEntityData(const moab::Interface &moab,
                                 const int pcomm_id)
    : moab(const_cast<moab::Interface &>(moab)), pcommID(pcomm_id),
      distributedMesh(true) {
  rval = moab.tag_get_handle("_RefParentHandle", th_RefParentHandle);
  MOAB_THROW(rval);
  rval = moab.tag_get_handle("_RefBitLevel", th_RefBitLevel);
  MOAB_THROW(rval);
}
BasicEntityData::~BasicEntityData() {}

// basic moab ent
BasicEntity::BasicEntity(
    const boost::shared_ptr<BasicEntityData> &basic_data_ptr,
    const EntityHandle ent)
    : basicDataPtr(basic_data_ptr), ent(ent) {
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
  ParallelComm *pcomm =
      ParallelComm::get_pcomm(&basicDataPtr->moab, basicDataPtr->pcommID);
  if (pcomm == NULL)
    THROW_MESSAGE("pcomm is null");
  rval = pcomm->get_owner_handle(ent, owner_proc, moab_owner_handle);
  MOAB_THROW(rval);
  part_proc = owner_proc;
}

unsigned char BasicEntity::getPStatus() const {
  ParallelComm *pcomm =
      ParallelComm::get_pcomm(&basicDataPtr->moab, basicDataPtr->pcommID);
  return *((unsigned char *)MoFEM::get_tag_ptr(
      basicDataPtr->moab, pcomm->pstatus_tag(), ent, NULL));
}

// ref moab ent
BitRefEdges MoFEM::RefElement::DummyBitRefEdges = BitRefEdges(0);
RefEntity::RefEntity(const boost::shared_ptr<BasicEntityData> &basic_data_ptr,
                     const EntityHandle ent)
    : BasicEntity(basic_data_ptr, ent) {}

EntityHandle *RefEntity::getParentEntPtr() const {
  return static_cast<EntityHandle *>(get_tag_ptr(
      basicDataPtr->moab, basicDataPtr->th_RefParentHandle, ent, NULL));
}

BitRefLevel *RefEntity::getBitRefLevelPtr() const {
  return static_cast<BitRefLevel *>(
      get_tag_ptr(basicDataPtr->moab, basicDataPtr->th_RefBitLevel, ent, NULL));
}

MoFEMErrorCode getParentEnt(moab::Interface &moab, Range ents,
                            std::vector<EntityHandle> vec_patent_ent) {

  MoFEMFunctionBegin;
  Tag th_ref_parent_handle;
  CHKERR moab.tag_get_handle("_RefParentHandle", th_ref_parent_handle);
  vec_patent_ent.resize(ents.size());
  CHKERR moab.tag_get_data(th_ref_parent_handle, ents,
                           &*vec_patent_ent.begin());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
RefEntity::getBitRefLevel(moab::Interface &moab, Range ents,
                          std::vector<BitRefLevel> &vec_bit_ref_level) {

  MoFEMFunctionBegin;
  Tag th_ref_bit_level;
  CHKERR moab.tag_get_handle("_RefBitLevel", th_ref_bit_level);
  vec_bit_ref_level.resize(ents.size());
  CHKERR moab.tag_get_data(th_ref_bit_level, ents, &*vec_bit_ref_level.begin());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode RefEntity::getBitRefLevel(
    moab::Interface &moab, Range ents,
    std::vector<const BitRefLevel *> &vec_ptr_bit_ref_level) {
  MoFEMFunctionBegin;
  Tag th_ref_bit_level;
  CHKERR moab.tag_get_handle("_RefBitLevel", th_ref_bit_level);
  vec_ptr_bit_ref_level.resize(ents.size());
  CHKERR moab.tag_get_by_ptr(
      th_ref_bit_level, ents,
      reinterpret_cast<const void **>(&*vec_ptr_bit_ref_level.begin()));
  MoFEMFunctionReturn(0);
}

std::ostream &operator<<(std::ostream &os, const RefEntity &e) {
  os << "ent " << e.ent;
  os << " pstatus " << std::bitset<8>(e.getPStatus());
  os << " owner ent " << e.getOwnerEnt();
  os << " owner proc " << e.getOwnerProc();
  os << " parent ent " << e.getParentEnt();
  // os << " BitRefLevel " << e.getBitRefLevel();
  os << " ent type " << e.getEntType();
  os << " ent parent type " << e.getParentEntType();
  return os;
}

// moab ent
FieldEntity::FieldEntity(const boost::shared_ptr<Field> &field_ptr,
                         const boost::shared_ptr<RefEntity> &ref_ent_ptr,
                         boost::shared_ptr<const int> &&t_max_order_ptr)
    : interface_Field<Field>(field_ptr), interface_RefEntity<RefEntity>(
                                             ref_ent_ptr),
      tagMaxOrderPtr(t_max_order_ptr) {
  globalUId = getGlobalUniqueIdCalculate();
  vectorAdaptorPtr = makeSharedFieldDataAdaptorPtr();

  if (!tagMaxOrderPtr)
    THROW_MESSAGE("Pointer to max order not set");
}

boost::shared_ptr<VectorAdaptor>
FieldEntity::makeSharedFieldDataAdaptorPtr() const {
  int size;
  double *ptr;
  switch (getEntType()) {
  case MBVERTEX:
    size = getNbOfCoeffs();
    ptr = static_cast<double *>(MoFEM::get_tag_ptr(
        sFieldPtr->moab, sFieldPtr->th_FieldData, sPtr->ent, &size));
    return boost::make_shared<VectorAdaptor>(
        size, ublas::shallow_array_adaptor<double>(size, ptr));
    break;
  default:
    ptr = static_cast<double *>(MoFEM::get_tag_ptr(
        sFieldPtr->moab, sFieldPtr->th_FieldData, sPtr->ent, &size));
    return boost::make_shared<VectorAdaptor>(
        size, ublas::shallow_array_adaptor<double>(size, ptr));
  }
}

UId *FieldEntity::getEntFieldDataLastUid = NULL;
double *FieldEntity::getEntFieldDataLastPtr = NULL;
int FieldEntity::getEntFieldDataLastSize = 0;
int FieldEntity::getEntFieldDataLastTagSize = 0;

VectorAdaptor FieldEntity::getEntFieldData() const {
  if (getEntFieldDataLastUid != &globalUId) {
    getEntFieldDataLastUid = const_cast<UId *>(&globalUId);
    switch (getEntType()) {
    case MBVERTEX:
      getEntFieldDataLastSize = getNbOfCoeffs();
      getEntFieldDataLastTagSize = getEntFieldDataLastSize;
      getEntFieldDataLastPtr = static_cast<double *>(MoFEM::get_tag_ptr(
          sFieldPtr->moab, sFieldPtr->th_FieldDataVerts, sPtr->ent, NULL));
      break;
    default:
      getEntFieldDataLastSize = getNbDofsOnEnt();
      getEntFieldDataLastPtr = static_cast<double *>(
          MoFEM::get_tag_ptr(sFieldPtr->moab, sFieldPtr->th_FieldData,
                             sPtr->ent, &getEntFieldDataLastTagSize));
      getEntFieldDataLastTagSize /= sizeof(double);
      if (PetscUnlikely(getEntFieldDataLastTagSize < getEntFieldDataLastSize))
        THROW_MESSAGE("Data inconsistency");
    }
  }
  return VectorAdaptor(getEntFieldDataLastSize,
                       ublas::shallow_array_adaptor<double>(
                           getEntFieldDataLastSize, getEntFieldDataLastPtr));
  }

  FieldEntity::~FieldEntity() {}
  std::ostream &operator<<(std::ostream &os, const FieldEntity &e) {
    os << "ent_global_uid "
       << (UId)e.getGlobalUniqueId()
       // << " ent_local_uid " << (UId)e.get_local_unique_id()
       << " entity " << e.getEnt() << " type " << e.getEntType() << " pstatus "
       << std::bitset<8>(e.getPStatus()) << " owner handle " << e.getOwnerEnt()
       << " owner proc " << e.getOwnerProc() << " order " << e.getMaxOrder()
       << " " << *e.sFieldPtr;
    return os;
  }

  void FieldEntity_change_order::operator()(FieldEntity *e) {

    moab::Interface &moab = e->sPtr->basicDataPtr->moab;
    const EntityHandle ent = e->getEnt();
    rval = moab.tag_set_data(e->sFieldPtr->th_AppOrder, &ent, 1, &order);
    MOAB_THROW(rval);
    unsigned int nb_dofs = e->getOrderNbDofs(order) * e->getNbOfCoeffs();

    double *tag_field_data;
    int tag_field_data_size;

    switch (e->getEntType()) {
    case MBVERTEX: {
      if (e->sFieldPtr->th_FieldDataVertsType == MB_TAG_SPARSE) {
        // Get pointer and size of field values tag
        rval = moab.tag_get_by_ptr(e->sFieldPtr->th_FieldDataVerts, &ent, 1,
                                   (const void **)&tag_field_data,
                                   &tag_field_data_size);
        if (nb_dofs) {
          if (nb_dofs != tag_field_data_size) {
            rval = moab.tag_set_data(e->sFieldPtr->th_FieldDataVerts, &ent, 1,
                                     &*data.begin());
            MOAB_THROW(rval);
          }
        } else if (rval == MB_SUCCESS) {
          rval = moab.tag_delete_data(e->sFieldPtr->th_FieldDataVerts, &ent, 1);
          MOAB_THROW(rval);
        }
      } else {
        rval = moab.tag_get_by_ptr(e->sFieldPtr->th_FieldDataVerts, &ent, 1,
                                   (const void **)&tag_field_data);
        MOAB_THROW(rval);
        rval = moab.tag_set_data(e->sFieldPtr->th_FieldDataVerts, &ent, 1,
                                 tag_field_data);
        MOAB_THROW(rval);
      }
    } break;
    default: {
      // Get pointer and size of field values tag
      rval = moab.tag_get_by_ptr(e->sFieldPtr->th_FieldData, &ent, 1,
                                 (const void **)&tag_field_data,
                                 &tag_field_data_size);
      // Tag exist and are some data on it
      if (rval == MB_SUCCESS) {
        // Check if size of filed values tag is correct
        if (nb_dofs * sizeof(FieldData) <= (unsigned int)tag_field_data_size)
          return;
        else if (nb_dofs == 0) {
          // Delete data on this entity
          rval = moab.tag_delete_data(e->sFieldPtr->th_FieldData, &ent, 1);
          MOAB_THROW(rval);
          return;
        }
        // Size of tag is different than new seize, so copy data to new
        // container
        data.resize(tag_field_data_size / sizeof(FieldData));
        FieldData *ptr_begin = (FieldData *)tag_field_data;
        FieldData *ptr_end = (FieldData *)tag_field_data +
                             tag_field_data_size / sizeof(FieldData);
        std::copy(ptr_begin, ptr_end, data.begin());
      }
      // Set new data
      if (nb_dofs > 0) {
        // Set field dof data
        data.resize(nb_dofs, 0);
        int tag_size[1];
        tag_size[0] = data.size() * sizeof(FieldData);
        void const *tag_data[] = {&data[0]};
        rval = moab.tag_set_by_ptr(e->sFieldPtr->th_FieldData, &ent, 1,
                                   tag_data, tag_size);
        MOAB_THROW(rval);
        rval = moab.tag_get_by_ptr(e->sFieldPtr->th_FieldData, &ent, 1,
                                   (const void **)&tag_field_data,
                                   &tag_field_data_size);
        MOAB_THROW(rval);
        if (nb_dofs != tag_field_data_size / sizeof(FieldData))
          THROW_MESSAGE("Data inconsistency");
      }
    }
    }
  }

} // namespace MoFEM
