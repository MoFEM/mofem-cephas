/** \file RefEntsMultiIndices.cpp
 * \brief Multi-index containers for entities
 */


#define IS_BUILDING_MB
#include <moab/Error.hpp>

namespace MoFEM {

BasicEntityData::BasicEntityData(const moab::Interface &moab,
                                 const int pcomm_id)
    : moab(const_cast<moab::Interface &>(moab)), pcommID(pcomm_id) {
  rval = moab.tag_get_handle("_RefParentHandle", th_RefParentHandle);
  MOAB_THROW(rval);
  rval = moab.tag_get_handle("_RefBitLevel", th_RefBitLevel);
  MOAB_THROW(rval);
  rval = moab.tag_get_handle("_MeshsetPartition", th_MeshsetPart);
  MOAB_THROW(rval);
}

boost::weak_ptr<BasicEntityData> RefEntityTmp<0>::basicDataPtr;
boost::weak_ptr<RefElement> RefEntityTmp<0>::refElementPtr;

RefEntityTmp<0>::RefEntityTmp(
    const boost::shared_ptr<BasicEntityData> &basic_data_ptr,
    const EntityHandle ent)
    : ent(ent) {
  entParentTagPtr = static_cast<EntityHandle *>(get_tag_ptr(
      basic_data_ptr->moab, basic_data_ptr->th_RefParentHandle, ent, NULL));
}

RefEntityTmp<0>::RefEntityTmp(
    const boost::shared_ptr<BasicEntityData> &basic_data_ptr,
    const EntityHandle ent, EntityHandle *ent_parent_tag_ptr)
    : ent(ent), entParentTagPtr(ent_parent_tag_ptr) {}

int RefEntityTmp<0>::getSideNumber() const {
  return getRefElementPtr()->getSideNumberPtr(ent)->side_number;
}

boost::shared_ptr<SideNumber> RefEntityTmp<0>::getSideNumberPtr() const {
  return getRefElementPtr()->getSideNumberPtr(ent);
}

std::ostream &operator<<(std::ostream &os, const RefEntity &e) {
  os << "ent " << e.ent;
  os << " pstatus " << std::bitset<8>(e.getPStatus());
  os << " owner ent " << e.getOwnerEnt();
  os << " owner proc " << e.getOwnerProc();
  os << " parent ent " << e.getParentEnt();
  os << " ent type " << e.getEntTypeName();
  os << " ent parent type " << e.getParentEntType();
  return os;
}

} // namespace MoFEM