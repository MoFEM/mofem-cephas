/** \file DofsMultiIndices.cpp
 * \brief Multi-index containers for DOFs
 */



namespace MoFEM {

// moab dof
DofEntity::DofEntity(const boost::shared_ptr<FieldEntity> &entity_ptr,
                     const ApproximationOrder dof_order,
                     const FieldCoefficientsNumber dof_rank, const DofIdx dof)
    : interface_FieldEntity<FieldEntity>(entity_ptr), dof(dof) {

#ifndef NDEBUG
  if (PetscUnlikely(!entity_ptr))
    THROW_MESSAGE("FieldEntity pointer not initialized");
  if (PetscUnlikely(!sPtr))
    THROW_MESSAGE("FieldEntity pointer not initialized");
  if (PetscUnlikely(!getFieldEntityPtr()))
    THROW_MESSAGE("FieldEntity pointer not initialized");
  // verify dof order
  if (PetscUnlikely(dof_order != getDofOrderMap()[dof])) {
    THROW_MESSAGE(
        "Inconsistent DOF order with order set before for entity " +
        boost::lexical_cast<std::string>(*entity_ptr) + "\n" + "dof->" +
        boost::lexical_cast<std::string>(dof) + " dof_order->" +
        boost::lexical_cast<std::string>(dof_order) +
        " != " + boost::lexical_cast<std::string>(getDofOrderMap()[dof]) +
        "<-getDofOrderMap()[dof]");
  }

  // verify dof rank
  if (PetscUnlikely(dof_rank != dof % getNbOfCoeffs()))
    THROW_MESSAGE("Inconsistent DOFs rank with index of DOF on entity");
#endif
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
