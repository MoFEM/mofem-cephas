/** \file RefEntsMultiIndices.cpp
 * \brief Multi-index containers for entities
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

// static moab::Error error;

BasicEntityData::BasicEntityData(const moab::Interface &moab,
                                 const int pcomm_id)
    : moab(const_cast<moab::Interface &>(moab)), pcommID(pcomm_id),
      distributedMesh(true) {
  rval = moab.tag_get_handle("_RefParentHandle", th_RefParentHandle);
  MOAB_THROW(rval);
  rval = moab.tag_get_handle("_RefBitLevel", th_RefBitLevel);
  MOAB_THROW(rval);
}

MoFEMErrorCode RefEntityTmp<0>::setFESideNumberPtr(const RefElement &fe) {
  MoFEMFunctionBegin;
  if (auto ptr = fe.getSideNumberPtr(ent))
    sideNumberPtr = ptr;
  else
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
            "Side entity on element ot found");
  MoFEMFunctionReturn(0);
}

boost::weak_ptr<BasicEntityData> RefEntityTmp<0>::basicDataPtr;
boost::weak_ptr<SideNumber> RefEntityTmp<0>::sideNumberPtr;

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

} // namespace MoFEM