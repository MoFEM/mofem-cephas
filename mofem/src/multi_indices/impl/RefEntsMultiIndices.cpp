/** \file RefEntsMultiIndices.cpp
 * \brief Multi-index containers for entities
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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