
/** \file moabField.cpp
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

#include<moabField.hpp>
#include<FEM.h>

namespace MoFEM {


PetscErrorCode moabField::SnesMethod::set_ctx(const snes_context ctx_) {
  PetscFunctionBegin;
  ctx = ctx_;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::SnesMethod::set_snes(SNES _snes) { 
  PetscFunctionBegin;
  snes = _snes;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::SnesMethod::set_x(Vec _x) {
  PetscFunctionBegin;
  snes_x = _x;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::SnesMethod::set_f(Vec _f) {
  PetscFunctionBegin;
  snes_f = _f;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::SnesMethod::set_A(Mat *_A) {
  PetscFunctionBegin;
  snes_A = _A;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::SnesMethod::set_B(Mat *_B) {
  PetscFunctionBegin;
  snes_A = _B;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::SnesMethod::set_flag(MatStructure *_flag) {
  PetscFunctionBegin;
  snes_flag = _flag;
  PetscFunctionReturn(0);
}
moabField::BasicMethod::BasicMethod():
  moabfields(NULL),ents_moabfield(NULL),dofs_moabfield(NULL),
  finite_elements(NULL),finite_elements_data(NULL),fem_adjacencies(NULL) {};

moabField::FEMethod::FEMethod(Interface& _moab): BasicMethod(), moab(_moab),
  problem_ptr(NULL),fe_ptr(NULL),data_multiIndex(NULL),
  row_multiIndex(NULL),col_multiIndex(NULL) {}

PetscErrorCode moabField::FEMethod::preProcess() {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class (preProcess)");
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::FEMethod::postProcess() {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class (postProcess)");
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::FEMethod::operator()() {   
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class (operator)");
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::FEMethod::set_problem(const MoFEMProblem *_problem_ptr) {
  PetscFunctionBegin;
  if(_problem_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  problem_ptr = _problem_ptr;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::FEMethod::set_fe(const NumeredMoFEMFE *_fe_ptr) {
  PetscFunctionBegin;
  if(_fe_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  fe_ptr = _fe_ptr;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::FEMethod::set_data_multIndex(const FEDofMoFEMEntity_multiIndex *_data_multiIndex) {
  PetscFunctionBegin;
  if(_data_multiIndex == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  data_multiIndex = _data_multiIndex;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::FEMethod::set_row_multIndex(const FENumeredDofMoFEMEntity_multiIndex *_row_multiIndex) {
  PetscFunctionBegin;
  if(_row_multiIndex == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  row_multiIndex = _row_multiIndex;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::FEMethod::set_col_multIndex(const FENumeredDofMoFEMEntity_multiIndex *_col_multiIndex) {
  PetscFunctionBegin;
  if(_col_multiIndex == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  col_multiIndex = _col_multiIndex;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::BasicMethod::set_moabfields(const MoFEMField_multiIndex *_moabfields) {
  PetscFunctionBegin;
  if(_moabfields == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  moabfields = _moabfields;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::BasicMethod::set_ents_multiIndex(const MoFEMEntity_multiIndex* _ents_moabfield) {
  PetscFunctionBegin;
  if(_ents_moabfield == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  ents_moabfield = _ents_moabfield;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::BasicMethod::set_dofs_multiIndex(const DofMoFEMEntity_multiIndex *_dofs_moabfield) {
  PetscFunctionBegin;
  if(_dofs_moabfield == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  dofs_moabfield = _dofs_moabfield;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::BasicMethod::set_fes_multiIndex(const MoFEMFE_multiIndex *_finite_elements) {
  PetscFunctionBegin;
  if(_finite_elements == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  finite_elements = _finite_elements;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::BasicMethod::set_adjacencies(const MoFEMAdjacencies_multiIndex *_fem_adjacencies) {
  PetscFunctionBegin;
  if(_fem_adjacencies == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  fem_adjacencies = _fem_adjacencies;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::BasicMethod::set_fes_data_multiIndex(const EntMoFEMFE_multiIndex *_finite_elements_data) {
  PetscFunctionBegin;
  if(_finite_elements_data == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  finite_elements_data = _finite_elements_data;
  PetscFunctionReturn(0);
}
moabField::EntMethod::EntMethod(Interface& _moab): BasicMethod(), moab(_moab),problem_ptr(NULL),dof_ptr(NULL) {}
PetscErrorCode moabField::EntMethod::preProcess() {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class");
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::EntMethod::postProcess() {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class");
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::EntMethod::operator()() {   
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class");
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::EntMethod::set_problem(const MoFEMProblem *_problem_ptr) {
  PetscFunctionBegin;
  if(_problem_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  problem_ptr = _problem_ptr;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField::EntMethod::set_dof(const NumeredDofMoFEMEntity *_dof_ptr) {
  PetscFunctionBegin;
  if(_dof_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  dof_ptr = _dof_ptr;
  PetscFunctionReturn(0);
}


}
