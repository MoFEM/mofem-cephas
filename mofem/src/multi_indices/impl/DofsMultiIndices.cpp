/** \file CoreDataStructures.cpp
 * \brief Multi-index container, data structures and other low-level functions
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

#include <Common.hpp>
#include <Includes.hpp>
#include <definitions.h>

#include <h1_hdiv_hcurl_l2.h>

#include <BCData.hpp>
#include <CoordSysMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <MaterialBlocks.hpp>
#include <TagMultiIndices.hpp>

namespace MoFEM {

// moab dof
DofEntity::DofEntity(const boost::shared_ptr<FieldEntity> &entity_ptr,
                     const ApproximationOrder dof_order,
                     const FieldCoefficientsNumber dof_rank, const DofIdx dof,
                     const bool is_active)
    : interface_FieldEntity<FieldEntity>(entity_ptr), active(is_active) {

  if (!entity_ptr) {
    THROW_MESSAGE("FieldEntity pinter not initialized");
  }
  if (!sPtr) {
    THROW_MESSAGE("FieldEntity pinter not initialized");
  }
  if (!getFieldEntityPtr()) {
    THROW_MESSAGE("FieldEntity pinter not initialized");
  }

  globalUId = getGlobalUniqueIdCalculate(dof, entity_ptr);

  // set order to DOF
  ApproximationOrder &order = getDofOrderMap()[dof];
  if (order != dof_order) {
    if (order != -1) {
      cerr << dof << " " << dof_order << " " << order;
      THROW_MESSAGE(
          "Order of DOFs inconsistent with order set before for entity");
    }
    order = dof_order;
  }

  // verify data consistency
  if (dof_rank != dof % getNbOfCoeffs()) {
    std::ostringstream ss;
    ss << dof_rank << " " << dof << " " << getNbOfCoeffs() << endl;
    ss << *entity_ptr;
    cerr << ss.str() << endl;
    THROW_MESSAGE("Inconsitent DOFs rank with index of DOF on entity");
  }
}

std::ostream &operator<<(std::ostream &os, const DofEntity &e) {
  os << "dof_uid " << e.getGlobalUniqueId() << " dof_order " << e.getDofOrder()
     << " dof_rank " << e.getDofCoeffIdx() << " dof " << e.getEntDofIdx()
     << " active " << e.active << " " << *(e.sFieldPtr);
  return os;
}

DofEntity_active_change::DofEntity_active_change(bool active)
    : aCtive(active) {}
void DofEntity_active_change::operator()(boost::shared_ptr<DofEntity> &dof) {
  dof->active = aCtive;
  if (aCtive && dof->getDofOrder() > dof->getMaxOrder()) {
    cerr << *dof << endl;
    THROW_MESSAGE("Set DoF active which has order larger than maximal order "
                  "set to entity");
  }
}

// numbered dof
NumeredDofEntity::NumeredDofEntity(
    const boost::shared_ptr<DofEntity> &dof_entity_ptr, const int dof_idx,
    const int petsc_gloabl_dof_idx, const int petsc_local_dof_idx,
    const int part

    )
    : interface_DofEntity<DofEntity>(dof_entity_ptr), dofIdx(dof_idx),
      petscGloablDofIdx(petsc_gloabl_dof_idx),
      petscLocalDofIdx(petsc_local_dof_idx), pArt(part) {}

std::ostream &operator<<(std::ostream &os, const NumeredDofEntity &e) {
  os << "idx " << e.dofIdx << " part " << e.pArt << " petsc idx "
     << e.petscGloablDofIdx << " ( " << e.petscLocalDofIdx << " ) "
     << *e.sFieldPtr;
  return os;
}

FEDofEntity::FEDofEntity(
    const boost::tuple<const boost::shared_ptr<SideNumber> &,
                       const boost::shared_ptr<DofEntity> &> &t)
    : BaseFEEntity(t.get<0>()), interface_DofEntity<DofEntity>(t.get<1>()) {}

FEDofEntity::FEDofEntity(const boost::shared_ptr<SideNumber> &side_number_ptr,
                         const boost::shared_ptr<DofEntity> &dof_ptr)
    : BaseFEEntity(side_number_ptr), interface_DofEntity<DofEntity>(
                                            dof_ptr) {}

std::ostream &operator<<(std::ostream &os, const FEDofEntity &e) {
  os << "local dof FiniteElement idx "
     << "side_number " << (int)e.sideNumberPtr->side_number << " "
     << "sense " << (int)e.sideNumberPtr->sense << " " << *e.sFieldPtr;
  return os;
}

FENumeredDofEntity::FENumeredDofEntity(
    const boost::shared_ptr<SideNumber> &side_number_ptr,
    const boost::shared_ptr<NumeredDofEntity> &dof_ptr)
    : BaseFEEntity(side_number_ptr),
      interface_NumeredDofEntity<NumeredDofEntity>(dof_ptr) {}

FENumeredDofEntity::FENumeredDofEntity(
    const boost::tuple<const boost::shared_ptr<SideNumber> &,
                       const boost::shared_ptr<NumeredDofEntity> &> &t)
    : BaseFEEntity(t.get<0>()), interface_NumeredDofEntity<NumeredDofEntity>(
                                       t.get<1>()) {}

std::ostream &operator<<(std::ostream &os, const FENumeredDofEntity &e) {
  os << "local dof FiniteElement idx "
     << "side_number " << (int)e.sideNumberPtr->side_number << " "
     << "sense " << (int)e.sideNumberPtr->sense << " " << *e.sFieldPtr;
  return os;
}

} // namespace MoFEM
