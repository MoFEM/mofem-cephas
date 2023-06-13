/** \file ProblemsMultiIndices.cpp
 * \brief Multindex contains for problems
 */


namespace MoFEM {

// moab problem
Problem::Problem(moab::Interface &moab, const EntityHandle meshset)
    : meshset(meshset), numeredRowDofsPtr(new NumeredDofEntity_multiIndex()),
      numeredColDofsPtr(new NumeredDofEntity_multiIndex()),
      numeredFiniteElementsPtr(new NumeredEntFiniteElement_multiIndex()),
      sequenceRowDofContainer(new SequenceDofContainer()),
      sequenceColDofContainer(new SequenceDofContainer()) {
  Tag th_ProblemId;
  MOAB_THROW(moab.tag_get_handle("_ProblemId", th_ProblemId));
  MOAB_THROW(
      moab.tag_get_by_ptr(th_ProblemId, &meshset, 1, (const void **)&tagId));
  Tag th_ProblemName;
  MOAB_THROW(moab.tag_get_handle("_ProblemName", th_ProblemName));
  MOAB_THROW(moab.tag_get_by_ptr(th_ProblemName, &meshset, 1,
                                 (const void **)&tagName, &tagNameSize));
  Tag th_ProblemFEId;
  MOAB_THROW(moab.tag_get_handle("_ProblemFEId", th_ProblemFEId));
  MOAB_THROW(moab.tag_get_by_ptr(th_ProblemFEId, &meshset, 1,
                                 (const void **)&tagBitFEId));
  Tag th_RefBitLevel;
  MOAB_THROW(moab.tag_get_handle("_RefBitLevel", th_RefBitLevel));
  MOAB_THROW(moab.tag_get_by_ptr(th_RefBitLevel, &meshset, 1,
                                 (const void **)&tagBitRefLevel));
  Tag th_RefBitLevel_Mask;
  MOAB_THROW(moab.tag_get_handle("_RefBitLevelMask", th_RefBitLevel_Mask));
  MOAB_THROW(moab.tag_get_by_ptr(th_RefBitLevel_Mask, &meshset, 1,
                                 (const void **)&tagBitRefLevelMask));
}

std::ostream &operator<<(std::ostream &os, const Problem &e) {
  os << "Problem id " << e.getId() << " name " << e.getName() << endl;
  os << "FiniteElement id " << e.getBitFEId() << endl;
  os << "BitRefLevel " << e.getBitRefLevel() << endl;
  os << "BitRefLevelMask " << e.getBitRefLevelMask();
  return os;
}

BitFEId Problem::getBitFEId() const { return *tagBitFEId; }

boost::weak_ptr<NumeredDofEntity>
Problem::getRowDofsByPetscGlobalDofIdx(DofIdx idx) const {
  MoFEMFunctionBeginHot;
  boost::weak_ptr<NumeredDofEntity> dof_weak_ptr;
  NumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit;
  dit = numeredRowDofsPtr->get<PetscGlobalIdx_mi_tag>().find(idx);
  if (dit != numeredRowDofsPtr->get<PetscGlobalIdx_mi_tag>().end())
    dof_weak_ptr = *dit;
  return dof_weak_ptr;
}

boost::weak_ptr<NumeredDofEntity>
Problem::getColDofsByPetscGlobalDofIdx(DofIdx idx) const {
  MoFEMFunctionBeginHot;
  boost::weak_ptr<NumeredDofEntity> dof_weak_ptr;
  NumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit;
  dit = numeredColDofsPtr->get<PetscGlobalIdx_mi_tag>().find(idx);
  if (dit != numeredColDofsPtr->get<PetscGlobalIdx_mi_tag>().end())
    dof_weak_ptr = *dit;
  return dof_weak_ptr;
}

MoFEMErrorCode Problem::getRowDofsByPetscGlobalDofIdx(
    DofIdx idx, const NumeredDofEntity **dof_ptr, MoFEMTypes bh) const {
  MoFEMFunctionBeginHot;
  auto weak_dof_ptr = getRowDofsByPetscGlobalDofIdx(idx);
  if (auto shared_dof_ptr = weak_dof_ptr.lock()) {
    *dof_ptr = shared_dof_ptr.get();
  } else {
    if (bh == MF_EXIST)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "row dof <%d> not found", idx);
    *dof_ptr = nullptr;
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Problem::getColDofsByPetscGlobalDofIdx(
    DofIdx idx, const NumeredDofEntity **dof_ptr, MoFEMTypes bh) const {
  MoFEMFunctionBeginHot;
  auto weak_dof_ptr = getColDofsByPetscGlobalDofIdx(idx);
  if (auto shared_dof_ptr = weak_dof_ptr.lock()) {
    *dof_ptr = shared_dof_ptr.get();
  } else {
    if (bh == MF_EXIST)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "row dof <%d> not found", idx);
    *dof_ptr = nullptr;
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Problem::getNumberOfElementsByNameAndPart(MPI_Comm comm, const std::string name,
                                          PetscLayout *layout) const {
  MoFEMFunctionBegin;
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  CHKERR PetscLayoutCreate(comm, layout);
  CHKERR PetscLayoutSetBlockSize(*layout, 1);
  const auto &fe_by_name_and_part =
      numeredFiniteElementsPtr->get<Composite_Name_And_Part_mi_tag>();
  int nb_elems;
  nb_elems = fe_by_name_and_part.count(boost::make_tuple(name, rank));
  CHKERR PetscLayoutSetLocalSize(*layout, nb_elems);
  CHKERR PetscLayoutSetUp(*layout);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Problem::getNumberOfElementsByPart(MPI_Comm comm,
                                                  PetscLayout *layout) const {
  MoFEMFunctionBegin;
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  CHKERR PetscLayoutCreate(comm, layout);
  CHKERR PetscLayoutSetBlockSize(*layout, 1);
  typedef NumeredEntFiniteElement_multiIndex::index<Part_mi_tag>::type FeByPart;
  const FeByPart &fe_by_part = numeredFiniteElementsPtr->get<Part_mi_tag>();
  int nb_elems;
  nb_elems = fe_by_part.count(rank);
  CHKERR PetscLayoutSetLocalSize(*layout, nb_elems);
  CHKERR PetscLayoutSetUp(*layout);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Problem::getDofByNameEntAndEntDofIdx(
    const int field_bit_number, const EntityHandle ent, const int ent_dof_idx,
    const RowColData row_or_col,
    boost::shared_ptr<NumeredDofEntity> &dof_ptr) const {
  MoFEMFunctionBegin;
  decltype(numeredRowDofsPtr) numered_dofs;
  switch (row_or_col) {
  case ROW:
    if (!numeredRowDofsPtr) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
              "Row numbered index in problem not allocated");
    }
    numered_dofs = numeredRowDofsPtr;
    break;
  case COL:
    if (!numeredColDofsPtr) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
              "Col numbered index in problem not allocated");
    }
    numered_dofs = numeredColDofsPtr;
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
            "Only ROW and COL is possible for 3rd argument");
  }

  auto it =
      numered_dofs->get<Unique_mi_tag>().find(DofEntity::getUniqueIdCalculate(
          ent_dof_idx,
          FieldEntity::getLocalUniqueIdCalculate(field_bit_number, ent)));
  if (it != numered_dofs->get<Unique_mi_tag>().end()) {
    dof_ptr = *it;
  } else {
    dof_ptr = boost::shared_ptr<NumeredDofEntity>();
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Problem::eraseElements(Range entities) const {
  MoFEMFunctionBegin;
  for (auto p = entities.pair_begin(); p != entities.pair_end(); ++p) {
    auto lo = numeredFiniteElementsPtr->get<Ent_mi_tag>().lower_bound(p->first);
    auto hi =
        numeredFiniteElementsPtr->get<Ent_mi_tag>().upper_bound(p->second);
    numeredFiniteElementsPtr->get<Ent_mi_tag>().erase(lo, hi);
  }
  MoFEMFunctionReturn(0);
}

void ProblemFiniteElementChangeBitAdd::operator()(Problem &p) {
  *(p.tagBitFEId) |= f_id;
}
void ProblemFiniteElementChangeBitUnSet::operator()(Problem &p) {
  *(p.tagBitFEId) &= ~f_id;
}
void ProblemZeroNbRowsChange::operator()(Problem &e) {
  e.nbDofsRow = 0;
  e.nbLocDofsRow = 0;
  e.nbGhostDofsRow = 0;
  e.numeredRowDofsPtr->clear();
}
void ProblemZeroNbColsChange::operator()(Problem &e) {
  e.nbDofsCol = 0;
  e.nbLocDofsCol = 0;
  e.nbGhostDofsCol = 0;
  e.numeredColDofsPtr->clear();
}
void ProblemClearNumeredFiniteElementsChange::operator()(Problem &e) {
  e.numeredFiniteElementsPtr->clear();
}

} // namespace MoFEM
