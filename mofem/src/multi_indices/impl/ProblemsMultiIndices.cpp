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
#include <CubitBCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <CoreDataStructures.hpp>

namespace MoFEM {

//moab problem
MoFEMProblem::MoFEMProblem(Interface &moab,const EntityHandle _meshset): meshset(_meshset) {
  ErrorCode rval;
  Tag th_ProblemId;
  rval = moab.tag_get_handle("_ProblemId",th_ProblemId); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemId,&meshset,1,(const void **)&tag_id_data); CHKERR(rval);
  Tag th_ProblemName;
  rval = moab.tag_get_handle("_ProblemName",th_ProblemName); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemName,&meshset,1,(const void **)&tag_name_data,&tag_name_size); CHKERR(rval);
  Tag th_ProblemNbDofsRow;
  rval = moab.tag_get_handle("_ProblemNbDofsRow",th_ProblemNbDofsRow); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemNbDofsRow,&meshset,1,(const void **)&tag_nbdof_data_row); CHKERR(rval);
  Tag th_ProblemNbDofsCol;
  rval = moab.tag_get_handle("_ProblemNbDofsCol",th_ProblemNbDofsCol); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemNbDofsCol,&meshset,1,(const void **)&tag_nbdof_data_col); CHKERR(rval);
  Tag th_ProblemLocalNbDofRow;
  rval = moab.tag_get_handle("_ProblemLocalNbDofsRow",th_ProblemLocalNbDofRow); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemLocalNbDofRow,&meshset,1,(const void **)&tag_local_nbdof_data_row); CHKERR(rval);
  Tag th_ProblemGhostNbDofRow;
  rval = moab.tag_get_handle("_ProblemGhostNbDofsRow",th_ProblemGhostNbDofRow); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemGhostNbDofRow,&meshset,1,(const void **)&tag_ghost_nbdof_data_row); CHKERR(rval);
  Tag th_ProblemLocalNbDofCol;
  rval = moab.tag_get_handle("_ProblemLocalNbDofsCol",th_ProblemLocalNbDofCol); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemLocalNbDofCol,&meshset,1,(const void **)&tag_local_nbdof_data_col); CHKERR(rval);
  Tag th_ProblemGhostNbDofCol;
  rval = moab.tag_get_handle("_ProblemGhostNbDofsCol",th_ProblemGhostNbDofCol); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemGhostNbDofCol,&meshset,1,(const void **)&tag_ghost_nbdof_data_col); CHKERR(rval);
  Tag th_ProblemFEId;
  rval = moab.tag_get_handle("_ProblemFEId",th_ProblemFEId); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemFEId,&meshset,1,(const void **)&tag_BitFEId_data); CHKERR(rval);
  Tag th_RefBitLevel;
  rval = moab.tag_get_handle("_RefBitLevel",th_RefBitLevel); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_RefBitLevel,&meshset,1,(const void **)&tag_BitRefLevel); CHKERR(rval);
  Tag th_RefBitLevel_Mask;
  rval = moab.tag_get_handle("_RefBitLevelMask",th_RefBitLevel_Mask); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_RefBitLevel_Mask,&meshset,1,(const void **)&tag_BitRefLevel_DofMask); CHKERR(rval);
}
ostream& operator<<(ostream& os,const MoFEMProblem& e) {
  os << "problem id " << e.get_id()
    << " MoFEMFiniteElement id " << e.get_BitFEId()
    << " name "<<e.get_name();
  return os;
}
BitFEId MoFEMProblem::get_BitFEId() const {
  return *tag_BitFEId_data;
}
PetscErrorCode MoFEMProblem::get_row_dofs_by_petsc_gloabl_dof_idx(DofIdx idx,const NumeredDofMoFEMEntity **dof_ptr) const {
  PetscFunctionBegin;
  NumeredDofMoFEMEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit;
  dit = numered_dofs_rows.get<PetscGlobalIdx_mi_tag>().find(idx);
  if(dit==numered_dofs_rows.get<PetscGlobalIdx_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"row dof <%d> not found",idx);
  }
  *dof_ptr = &*dit;
  PetscFunctionReturn(0);
}
PetscErrorCode MoFEMProblem::get_col_dofs_by_petsc_gloabl_dof_idx(DofIdx idx,const NumeredDofMoFEMEntity **dof_ptr) const {
  PetscFunctionBegin;
  NumeredDofMoFEMEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit;
  dit = numered_dofs_cols.get<PetscGlobalIdx_mi_tag>().find(idx);
  if(dit==numered_dofs_cols.get<PetscGlobalIdx_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"row dof <%d> not found",idx);
  }
  *dof_ptr = &*dit;
  PetscFunctionReturn(0);
}
void ProblemFiniteElementChangeBitAdd::operator()(MoFEMProblem &p) {
  *(p.tag_BitFEId_data) |= f_id;
}
void ProblemFiniteElementChangeBitUnSet::operator()(MoFEMProblem &p) {
  *(p.tag_BitFEId_data) &= ~f_id;
}
ProblemAddRowDof::ProblemAddRowDof(const DofMoFEMEntity *_dof_ptr): dof_ptr(_dof_ptr) {
  assert(dof_ptr->active);
}
void ProblemAddRowDof::operator()(MoFEMProblem &e) {
  p = e.numered_dofs_rows.insert(NumeredDofMoFEMEntity(dof_ptr));
  if(p.second) {
    (*(DofIdx*)e.tag_nbdof_data_row)++;
  }
}
ProblemAddColDof::ProblemAddColDof(const DofMoFEMEntity *_dof_ptr): dof_ptr(_dof_ptr) {}
void ProblemAddColDof::operator()(MoFEMProblem &e) {
  p = e.numered_dofs_cols.insert(NumeredDofMoFEMEntity(dof_ptr));
  if(p.second) {
    (*(DofIdx*)e.tag_nbdof_data_col)++;
  }
}
void ProblemZeroNbRowsChange::operator()(MoFEMProblem &e) {
  (*(DofIdx*)e.tag_nbdof_data_row) = 0;
  (*(DofIdx*)e.tag_local_nbdof_data_row) = 0;
  (*(DofIdx*)e.tag_ghost_nbdof_data_row) = 0;
  e.numered_dofs_rows.clear();
}
void ProblemZeroNbColsChange::operator()(MoFEMProblem &e) {
  (*(DofIdx*)e.tag_nbdof_data_col) = 0;
  (*(DofIdx*)e.tag_local_nbdof_data_col) = 0;
  (*(DofIdx*)e.tag_ghost_nbdof_data_col) = 0;
  e.numered_dofs_cols.clear();
}
void ProblemClearNumeredFiniteElementsChange::operator()(MoFEMProblem &e) {
  e.numeredFiniteElements.clear();
}
void ProblemRowNumberChange::operator()(MoFEMProblem &e) {
  NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator dit;
  dit = e.numered_dofs_rows.get<Unique_mi_tag>().begin();
  int idx = 0;
  for(;dit!=e.numered_dofs_rows.get<Unique_mi_tag>().end();dit++,idx++) {
    bool success =
      e.numered_dofs_rows.modify(dit,NumeredDofMoFEMEntity_mofem_index_change(idx));
    if(!success) {
      throw "modification unsuccessful";
    }
  }
}
void ProblemColNumberChange::operator()(MoFEMProblem &e) {
  NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator dit;
  dit = e.numered_dofs_cols.get<Unique_mi_tag>().begin();
  int idx = 0;
  for(;dit!=e.numered_dofs_cols.get<Unique_mi_tag>().end();dit++,idx++) {
    bool success =
      e.numered_dofs_cols.modify(dit,NumeredDofMoFEMEntity_mofem_index_change(idx));
    if(!success) {
      throw "modification unsuccessful";
    }
  }
}

}
