/** \file ProblemsMultiIndices.cpp
 * \brief Multindex contains for problems
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

namespace MoFEM {

// moab problem
Problem::Problem(moab::Interface &moab, const EntityHandle meshset)
    : meshset(meshset), numeredRowDofs(new NumeredDofEntity_multiIndex()),
      numeredColDofs(new NumeredDofEntity_multiIndex()),
      numeredFiniteElements(new NumeredEntFiniteElement_multiIndex()),
      sequenceRowDofContainer(new SequenceDofContainer()),
      sequenceColDofContainer(new SequenceDofContainer()) {
  ErrorCode rval;
  Tag th_ProblemId;
  rval = moab.tag_get_handle("_ProblemId", th_ProblemId);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemId, &meshset, 1, (const void **)&tagId);
  MOAB_THROW(rval);
  Tag th_ProblemName;
  rval = moab.tag_get_handle("_ProblemName", th_ProblemName);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemName, &meshset, 1,
                             (const void **)&tagName, &tagNameSize);
  MOAB_THROW(rval);
  Tag th_ProblemFEId;
  rval = moab.tag_get_handle("_ProblemFEId", th_ProblemFEId);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemFEId, &meshset, 1,
                             (const void **)&tagBitFEId);
  MOAB_THROW(rval);
  Tag th_RefBitLevel;
  rval = moab.tag_get_handle("_RefBitLevel", th_RefBitLevel);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitLevel, &meshset, 1,
                             (const void **)&tagBitRefLevel);
  MOAB_THROW(rval);
  Tag th_RefBitLevel_Mask;
  rval = moab.tag_get_handle("_RefBitLevelMask", th_RefBitLevel_Mask);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitLevel_Mask, &meshset, 1,
                             (const void **)&tagMaskBitRefLevel);
  MOAB_THROW(rval);
}

Problem::~Problem() {}

std::ostream &operator<<(std::ostream &os, const Problem &e) {
  os << "problem id " << e.getId() << " FiniteElement id " << e.getBitFEId()
     << " name " << e.getName();
  return os;
}

BitFEId Problem::getBitFEId() const { return *tagBitFEId; }

boost::weak_ptr<NumeredDofEntity>
Problem::getRowDofsByPetscGlobalDofIdx(DofIdx idx) const {
  MoFEMFunctionBeginHot;
  boost::weak_ptr<NumeredDofEntity> dof_weak_ptr;
  NumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit;
  dit = numeredRowDofs->get<PetscGlobalIdx_mi_tag>().find(idx);
  if (dit != numeredRowDofs->get<PetscGlobalIdx_mi_tag>().end())
    dof_weak_ptr = *dit;
  return dof_weak_ptr;
}

boost::weak_ptr<NumeredDofEntity>
Problem::getColDofsByPetscGlobalDofIdx(DofIdx idx) const {
  MoFEMFunctionBeginHot;
  boost::weak_ptr<NumeredDofEntity> dof_weak_ptr;
  NumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit;
  dit = numeredColDofs->get<PetscGlobalIdx_mi_tag>().find(idx);
  if (dit != numeredColDofs->get<PetscGlobalIdx_mi_tag>().end())
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
  const NumeredEntFiniteElementbyNameAndPart &fe_by_name_and_part =
      numeredFiniteElements->get<Composite_Name_And_Part_mi_tag>();
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
  const FeByPart &fe_by_part = numeredFiniteElements->get<Part_mi_tag>();
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
  decltype(numeredRowDofs) numered_dofs;
  switch (row_or_col) {
  case ROW:
    if (!numeredRowDofs) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
              "Row numbered index in problem not allocated");
    }
    numered_dofs = numeredRowDofs;
    break;
  case COL:
    if (!numeredColDofs) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
              "Col numbered index in problem not allocated");
    }
    numered_dofs = numeredColDofs;
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
  e.numeredRowDofs->clear();
}
void ProblemZeroNbColsChange::operator()(Problem &e) {
  e.nbDofsCol = 0;
  e.nbLocDofsCol = 0;
  e.nbGhostDofsCol = 0;
  e.numeredColDofs->clear();
}
void ProblemClearNumeredFiniteElementsChange::operator()(Problem &e) {
  e.numeredFiniteElements->clear();
}

} // namespace MoFEM
