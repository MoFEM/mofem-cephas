/** \file CoreDataStructures.cpp
 * \brief Myltindex contains, data structures and other low-level functions
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
#include <ProblemsMultiIndices.hpp>
#include <CoreDataStructures.hpp>

namespace MoFEM {

//moab problem
Problem::Problem(Interface &moab,const EntityHandle meshset):
  meshset(meshset),
  numeredDofsRows(boost::shared_ptr<NumeredDofEntity_multiIndex>(new NumeredDofEntity_multiIndex())),
  numeredDofsCols(boost::shared_ptr<NumeredDofEntity_multiIndex>(new NumeredDofEntity_multiIndex())),
  // numered_dofs_rows(numeredDofsRows), // this is deprecated
  // numered_dofs_cols(numeredDofsCols), // this is deprecated
  sequenceRowDofContainer(boost::make_shared<SequenceDofContainer>()),
  sequenceColDofContainer(boost::make_shared<SequenceDofContainer>()) {
  ErrorCode rval;
  Tag th_ProblemId;
  rval = moab.tag_get_handle("_ProblemId",th_ProblemId); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemId,&meshset,1,(const void **)&tag_id_data); MOAB_THROW(rval);
  Tag th_ProblemName;
  rval = moab.tag_get_handle("_ProblemName",th_ProblemName); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemName,&meshset,1,(const void **)&tag_name_data,&tag_name_size); MOAB_THROW(rval);
  Tag th_ProblemNbDofsRow;
  rval = moab.tag_get_handle("_ProblemNbDofsRow",th_ProblemNbDofsRow); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemNbDofsRow,&meshset,1,(const void **)&tag_nbdof_data_row); MOAB_THROW(rval);
  Tag th_ProblemNbDofsCol;
  rval = moab.tag_get_handle("_ProblemNbDofsCol",th_ProblemNbDofsCol); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemNbDofsCol,&meshset,1,(const void **)&tag_nbdof_data_col); MOAB_THROW(rval);
  Tag th_ProblemLocalNbDofRow;
  rval = moab.tag_get_handle("_ProblemLocalNbDofsRow",th_ProblemLocalNbDofRow); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemLocalNbDofRow,&meshset,1,(const void **)&tag_local_nbdof_data_row); MOAB_THROW(rval);
  Tag th_ProblemGhostNbDofRow;
  rval = moab.tag_get_handle("_ProblemGhostNbDofsRow",th_ProblemGhostNbDofRow); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemGhostNbDofRow,&meshset,1,(const void **)&tag_ghost_nbdof_data_row); MOAB_THROW(rval);
  Tag th_ProblemLocalNbDofCol;
  rval = moab.tag_get_handle("_ProblemLocalNbDofsCol",th_ProblemLocalNbDofCol); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemLocalNbDofCol,&meshset,1,(const void **)&tag_local_nbdof_data_col); MOAB_THROW(rval);
  Tag th_ProblemGhostNbDofCol;
  rval = moab.tag_get_handle("_ProblemGhostNbDofsCol",th_ProblemGhostNbDofCol); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemGhostNbDofCol,&meshset,1,(const void **)&tag_ghost_nbdof_data_col); MOAB_THROW(rval);
  Tag th_ProblemFEId;
  rval = moab.tag_get_handle("_ProblemFEId",th_ProblemFEId); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemFEId,&meshset,1,(const void **)&tag_BitFEId_data); MOAB_THROW(rval);
  Tag th_RefBitLevel;
  rval = moab.tag_get_handle("_RefBitLevel",th_RefBitLevel); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitLevel,&meshset,1,(const void **)&tag_BitRefLevel); MOAB_THROW(rval);
  Tag th_RefBitLevel_Mask;
  rval = moab.tag_get_handle("_RefBitLevelMask",th_RefBitLevel_Mask); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitLevel_Mask,&meshset,1,(const void **)&tag_MaskBitRefLevel); MOAB_THROW(rval);
}

Problem::~Problem() {
}

std::ostream& operator<<(std::ostream& os,const Problem& e) {
  os << "problem id " << e.getId()
    << " FiniteElement id " << e.getBitFEId()
    << " name "<<e.getName();
  return os;
}

BitFEId Problem::getBitFEId() const {
  return *tag_BitFEId_data;
}

PetscErrorCode Problem::getRowDofsByPetscGlobalDofIdx(DofIdx idx,const NumeredDofEntity **dof_ptr) const {
  PetscFunctionBegin;
  NumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit;
  dit = numeredDofsRows->get<PetscGlobalIdx_mi_tag>().find(idx);
  if(dit==numeredDofsRows->get<PetscGlobalIdx_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"row dof <%d> not found",idx);
  }
  *dof_ptr = &*(*dit);
  PetscFunctionReturn(0);
}

PetscErrorCode Problem::getColDofsByPetscGlobalDofIdx(DofIdx idx,const NumeredDofEntity **dof_ptr) const {
  PetscFunctionBegin;
  NumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit;
  dit = numeredDofsCols->get<PetscGlobalIdx_mi_tag>().find(idx);
  if(dit==numeredDofsCols->get<PetscGlobalIdx_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"row dof <%d> not found",idx);
  }
  *dof_ptr = &*(*dit);
  PetscFunctionReturn(0);
}

PetscErrorCode Problem::getNumberOfElementsByNameAndPart(MPI_Comm comm,const std::string name,PetscLayout *layout) const {
  PetscFunctionBegin;

  int size, rank;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  ierr = PetscLayoutCreate(comm,layout); CHKERRQ(ierr);
  ierr = PetscLayoutSetBlockSize(*layout,1); CHKERRQ(ierr);
  const NumeredEntFiniteElementbyNameAndPart &fe_by_name_and_part = numeredFiniteElements.get<Composite_Name_And_Part_mi_tag>();
  int nb_elems;
  nb_elems = fe_by_name_and_part.count(boost::make_tuple(name,rank));
  ierr = PetscLayoutSetLocalSize(*layout,nb_elems); CHKERRQ(ierr);
  ierr = PetscLayoutSetUp(*layout); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Problem::getNumberOfElementsByPart(MPI_Comm comm,PetscLayout *layout) const {
  PetscFunctionBegin;

  int size, rank;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  ierr = PetscLayoutCreate(comm,layout); CHKERRQ(ierr);
  ierr = PetscLayoutSetBlockSize(*layout,1); CHKERRQ(ierr);
  typedef NumeredEntFiniteElement_multiIndex::index<Part_mi_tag>::type FeByPart;
  const FeByPart &fe_by_part = numeredFiniteElements.get<Part_mi_tag>();
  int nb_elems;
  nb_elems = fe_by_part.count(rank);
  ierr = PetscLayoutSetLocalSize(*layout,nb_elems); CHKERRQ(ierr);
  ierr = PetscLayoutSetUp(*layout); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Problem::getDofByNameEntAndEntDofIdx(
  const string name,
  const EntityHandle ent,
  const int ent_dof_idx,
  const RowColData row_or_col,
  boost::shared_ptr<NumeredDofEntity> &dof_ptr
) const {
  PetscFunctionBegin;
  typedef NumeredDofEntity_multiIndex::index<
  Composite_Name_And_Ent_And_EntDofIdx_mi_tag
  >::type NumberdDofByNameEntAndEndDofIdx;
  NumberdDofByNameEntAndEndDofIdx::iterator it;
  // not use shared pointer is local here, direct pointer is more efficient
  NumeredDofEntity_multiIndex *numered_dofs;
  switch (row_or_col) {
    case ROW:
    if(!numeredDofsRows) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA,"Row numbered index in problem not allocated");
    }
    numered_dofs = numeredDofsRows.get();
    break;
    case COL:
    if(!numeredDofsRows) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA,"Col numbered index in problem not allocated");
    }
    numered_dofs = numeredDofsCols.get();
    break;
    default:
    SETERRQ(
      PETSC_COMM_SELF,MOFEM_INVALID_DATA,"Only ROW and COL is possible for 3rd argument"
    );
  }
  it = numered_dofs->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().find(
    boost::make_tuple(name,ent,ent_dof_idx)
  );
  if(it!=numered_dofs->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().end()) {
    dof_ptr = *it;
  } else {
    dof_ptr = boost::shared_ptr<NumeredDofEntity>();
  }
  PetscFunctionReturn(0);
}


void ProblemFiniteElementChangeBitAdd::operator()(Problem &p) {
  *(p.tag_BitFEId_data) |= f_id;
}
void ProblemFiniteElementChangeBitUnSet::operator()(Problem &p) {
  *(p.tag_BitFEId_data) &= ~f_id;
}
void ProblemZeroNbRowsChange::operator()(Problem &e) {
  (*(DofIdx*)e.tag_nbdof_data_row) = 0;
  (*(DofIdx*)e.tag_local_nbdof_data_row) = 0;
  (*(DofIdx*)e.tag_ghost_nbdof_data_row) = 0;
  e.numeredDofsRows->clear();
}
void ProblemZeroNbColsChange::operator()(Problem &e) {
  (*(DofIdx*)e.tag_nbdof_data_col) = 0;
  (*(DofIdx*)e.tag_local_nbdof_data_col) = 0;
  (*(DofIdx*)e.tag_ghost_nbdof_data_col) = 0;
  e.numeredDofsCols->clear();
}
void ProblemClearNumeredFiniteElementsChange::operator()(Problem &e) {
  e.numeredFiniteElements.clear();
}

}
