/** \file DofsMultiIndices.cpp
 * \brief Multi-index containers for DOFs
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

namespace MoFEM {

// moab dof
DofEntity::DofEntity(const boost::shared_ptr<FieldEntity> &entity_ptr,
                     const ApproximationOrder dof_order,
                     const FieldCoefficientsNumber dof_rank, const DofIdx dof)
    : interface_FieldEntity<FieldEntity>(entity_ptr), dof(dof) {

  if (PetscUnlikely(!entity_ptr))
    THROW_MESSAGE("FieldEntity pointer not initialized");
  if (PetscUnlikely(!sPtr))
    THROW_MESSAGE("FieldEntity pointer not initialized");
  if (PetscUnlikely(!getFieldEntityPtr()))
    THROW_MESSAGE("FieldEntity pointer not initialized");
  // verify dof order
  if (PetscUnlikely(dof_order != getDofOrderMap()[dof]))
    THROW_MESSAGE(
        "Inconsistent DOF order with order set before for entity " +
        boost::lexical_cast<std::string>(*entity_ptr) + " dof_order->" +
        boost::lexical_cast<std::string>(dof_order) +
        " != " + boost::lexical_cast<std::string>(getDofOrderMap()[dof]) +
        "<-getDofOrderMap()[dof]");
  // verify dof rank
  if (PetscUnlikely(dof_rank != dof % getNbOfCoeffs()))
    THROW_MESSAGE("Inconsistent DOFs rank with index of DOF on entity");
}

std::ostream &operator<<(std::ostream &os, const DofEntity &e) {
  os << "dof_uid " << e.getLocalUniqueId() << " dof_order " << e.getDofOrder()
     << " dof_rank " << e.getDofCoeffIdx() << " dof " << e.getEntDofIdx()
     << " active " << (e.dof < 0 ? false : true) << " "
     << *e.getFieldEntityPtr();
  return os;
}

DofEntity_active_change::DofEntity_active_change(bool active)
    : aCtive(active) {}
void DofEntity_active_change::operator()(boost::shared_ptr<DofEntity> &dof) {

  if (aCtive)
    dof->dof = std::abs(dof->dof);
  else
    dof->dof = -std::abs(dof->dof);

  if (PetscUnlikely(aCtive && dof->getDofOrder() > dof->getMaxOrder())) {
    MOFEM_LOG_CHANNEL("SELF");
    MOFEM_LOG_TAG("SELF", "DofEntity_active_change");
    MOFEM_LOG("SELF", Sev::error) << *dof;
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
      petscLocalDofIdx(petsc_local_dof_idx), pArt(part) {
}

std::ostream &operator<<(std::ostream &os, const NumeredDofEntity &e) {
  os << "idx " << e.dofIdx << " part " << e.pArt << " petsc idx "
     << e.petscGloablDofIdx << " ( " << e.petscLocalDofIdx << " ) "
     << *e.getFieldEntityPtr();
  return os;
}

std::ostream &operator<<(std::ostream &os, const FEDofEntity &e) {
  os << "local dof FiniteElement idx "
     << "side_number " << static_cast<int>(e.getSideNumberPtr()->side_number)
     << " "
     << "sense " << static_cast<int>(e.getSideNumberPtr()->sense) << " "
     << *e.getFieldEntityPtr();
  return os;
}

std::ostream &operator<<(std::ostream &os, const FENumeredDofEntity &e) {
  os << "local dof FiniteElement idx "
     << "side_number " << (int)e.getSideNumberPtr()->side_number << " "
     << "sense " << (int)e.getSideNumberPtr()->sense << " "
     << *e.getFieldEntityPtr();
  return os;
}

} // namespace MoFEM
