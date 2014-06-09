
/** \file FieldInterface.cpp
 * \brief Myltindex containes, data structures and other low-level functions 
 * 
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 *
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
 * It can be freely used for educational and research purposes 
 * by other institutions. If you use this softwre pleas cite my work. 
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

#include<FieldInterface.hpp>
#include<FEM.h>

namespace MoFEM {

//SNES
PetscErrorCode FieldInterface::SnesMethod::set_snes_ctx(const snes_context ctx_) {
  PetscFunctionBegin;
  snes_ctx = ctx_;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::SnesMethod::set_snes(SNES _snes) { 
  PetscFunctionBegin;
  snes = _snes;
  PetscFunctionReturn(0);
}
//TS
PetscErrorCode FieldInterface::TSMethod::set_ts_ctx(const ts_context ctx_) {
  PetscFunctionBegin;
  ts_ctx = ctx_;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::TSMethod::set_ts(TS _ts) {
  PetscFunctionBegin;
  ts = _ts;
  PetscFunctionReturn(0);
}
//BasicMethod
FieldInterface::BasicMethod::BasicMethod():
  refinedMoFemEntities(NULL),refinedMoFemElements(NULL),
  problem_ptr(NULL),moabfields(NULL),ents_moabfield(NULL),dofs_moabfield(NULL),
  finite_elements(NULL),finite_elements_moabents(NULL),fem_adjacencies(NULL) {};
PetscErrorCode FieldInterface::BasicMethod::set_problem(const MoFEMProblem *_problem_ptr) {
  PetscFunctionBegin;
  if(_problem_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  problem_ptr = _problem_ptr;
  PetscFunctionReturn(0);
}
FieldInterface::FEMethod::FEMethod(): BasicMethod(),
  fe_ptr(NULL),data_multiIndex(NULL),
  row_multiIndex(NULL),col_multiIndex(NULL) {}

PetscErrorCode FieldInterface::FEMethod::preProcess() {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class (preProcess)");
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::FEMethod::postProcess() {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class (postProcess)");
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::FEMethod::operator()() {   
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class (operator)");
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::FEMethod::set_fe(const NumeredMoFEMFiniteElement *_fe_ptr) {
  PetscFunctionBegin;
  if(_fe_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  fe_ptr = _fe_ptr;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::FEMethod::set_data_multIndex(const FEDofMoFEMEntity_multiIndex *_data_multiIndex) {
  PetscFunctionBegin;
  if(_data_multiIndex == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  data_multiIndex = _data_multiIndex;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::FEMethod::set_row_multIndex(const FENumeredDofMoFEMEntity_multiIndex *_row_multiIndex) {
  PetscFunctionBegin;
  if(_row_multiIndex == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  row_multiIndex = _row_multiIndex;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::FEMethod::set_col_multIndex(const FENumeredDofMoFEMEntity_multiIndex *_col_multiIndex) {
  PetscFunctionBegin;
  if(_col_multiIndex == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  col_multiIndex = _col_multiIndex;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::BasicMethod::set_fields(const MoFEMField_multiIndex *_moabfields) {
  PetscFunctionBegin;
  if(_moabfields == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  moabfields = _moabfields;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::BasicMethod::set_ents_multiIndex(
    const RefMoFEMEntity_multiIndex *_refinedMoFemEntities,
    const MoFEMEntity_multiIndex* _ents_moabfield) {
  PetscFunctionBegin;
  if(_refinedMoFemEntities == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  refinedMoFemEntities = _refinedMoFemEntities;
  if(_ents_moabfield == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  ents_moabfield = _ents_moabfield;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::BasicMethod::set_dofs_multiIndex(const DofMoFEMEntity_multiIndex *_dofs_moabfield) {
  PetscFunctionBegin;
  if(_dofs_moabfield == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  dofs_moabfield = _dofs_moabfield;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::BasicMethod::set_fes_multiIndex(
    const RefMoFEMElement_multiIndex *_refinedMoFemElements,
    const MoFEMFiniteElement_multiIndex *_finite_elements) {
  PetscFunctionBegin;
  if(_refinedMoFemElements == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  refinedMoFemElements = _refinedMoFemElements;
  if(_finite_elements == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  finite_elements = _finite_elements;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::BasicMethod::set_adjacencies(const MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex *_fem_adjacencies) {
  PetscFunctionBegin;
  if(_fem_adjacencies == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  fem_adjacencies = _fem_adjacencies;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::BasicMethod::set_fes_data_multiIndex(const EntMoFEMFiniteElement_multiIndex *_finite_elements_moabents) {
  PetscFunctionBegin;
  if(_finite_elements_moabents == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  finite_elements_moabents = _finite_elements_moabents;
  PetscFunctionReturn(0);
}
FieldInterface::EntMethod::EntMethod(): BasicMethod(),dof_ptr(NULL) {}
PetscErrorCode FieldInterface::EntMethod::preProcess() {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class");
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::EntMethod::postProcess() {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class");
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::EntMethod::operator()() {   
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class");
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::EntMethod::set_dof(const DofMoFEMEntity *_dof_ptr) {
  PetscFunctionBegin;
  if(_dof_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  dof_ptr = _dof_ptr;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::EntMethod::set_numered_dof(const NumeredDofMoFEMEntity *_dof_ptr) {
  PetscFunctionBegin;
  //if(_dof_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  dof_numered_ptr = _dof_ptr;
  PetscFunctionReturn(0);
}


}
