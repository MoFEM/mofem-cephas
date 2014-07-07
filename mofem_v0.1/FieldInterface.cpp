
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
PetscErrorCode FieldInterface::SnesMethod::set_snes_ctx(const SNESContext ctx_) {
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
PetscErrorCode FieldInterface::TSMethod::set_ts_ctx(const TSContext ctx_) {
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
  refinedEntitiesPtr(NULL),refinedFiniteElementsPtr(NULL),
  problemPtr(NULL),fieldsPtr(NULL),entitiesPtr(NULL),dofsPtr(NULL),
  finiteElementsPtr(NULL),finiteElementsEntitiesPtr(NULL),adjacenciesPtr(NULL) {};
PetscErrorCode FieldInterface::BasicMethod::setProblem(const MoFEMProblem *_problemPtr) {
  PetscFunctionBegin;
  if(_problemPtr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  problemPtr = _problemPtr;
  PetscFunctionReturn(0);
}
FieldInterface::FEMethod::FEMethod(): BasicMethod(),
  fePtr(NULL),dataPtr(NULL),
  rowPtr(NULL),colPtr(NULL) {}

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
PetscErrorCode FieldInterface::FEMethod::setFE(const NumeredMoFEMFiniteElement *_fe_ptr) {
  PetscFunctionBegin;
  if(_fe_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  fePtr = _fe_ptr;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::FEMethod::setData(const FEDofMoFEMEntity_multiIndex *_data_ptr) {
  PetscFunctionBegin;
  if(_data_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  dataPtr = _data_ptr;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::FEMethod::setRowData(const FENumeredDofMoFEMEntity_multiIndex *_row_ptr) {
  PetscFunctionBegin;
  if(_row_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  rowPtr = _row_ptr;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::FEMethod::setColData(const FENumeredDofMoFEMEntity_multiIndex *_col_ptr) {
  PetscFunctionBegin;
  if(_col_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  colPtr = _col_ptr;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::BasicMethod::setFields(const MoFEMField_multiIndex *_fieldsPtr) {
  PetscFunctionBegin;
  if(_fieldsPtr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  fieldsPtr = _fieldsPtr;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::BasicMethod::setEnts(
    const RefMoFEMEntity_multiIndex *_refinedEntitiesPtr,
    const MoFEMEntity_multiIndex* _entitiesPtr) {
  PetscFunctionBegin;
  if(_refinedEntitiesPtr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  refinedEntitiesPtr = _refinedEntitiesPtr;
  if(_entitiesPtr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  entitiesPtr = _entitiesPtr;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::BasicMethod::setDofs(const DofMoFEMEntity_multiIndex *_dofsPtr) {
  PetscFunctionBegin;
  if(_dofsPtr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  dofsPtr = _dofsPtr;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::BasicMethod::setFiniteElements(
    const RefMoFEMElement_multiIndex *_refinedFiniteElementsPtr,
    const MoFEMFiniteElement_multiIndex *_finiteElementsPtr) {
  PetscFunctionBegin;
  if(_refinedFiniteElementsPtr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  refinedFiniteElementsPtr = _refinedFiniteElementsPtr;
  if(_finiteElementsPtr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  finiteElementsPtr = _finiteElementsPtr;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::BasicMethod::setAdjacencies(const MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex *_adjacenciesPtr) {
  PetscFunctionBegin;
  if(_adjacenciesPtr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  adjacenciesPtr = _adjacenciesPtr;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::BasicMethod::setFiniteElementsEntities(const EntMoFEMFiniteElement_multiIndex *_finiteElementsEntitiesPtr) {
  PetscFunctionBegin;
  if(_finiteElementsEntitiesPtr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  finiteElementsEntitiesPtr = _finiteElementsEntitiesPtr;
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
PetscErrorCode FieldInterface::EntMethod::setDof(const DofMoFEMEntity *_dof_ptr) {
  PetscFunctionBegin;
  if(_dof_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  dof_ptr = _dof_ptr;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldInterface::EntMethod::set_numered_dof(const NumeredDofMoFEMEntity *_dof_ptr) {
  PetscFunctionBegin;
  //if(_dof_ptr == NULL) SETERRQ(PETSC_COMM_SELF,1,"can not be NULL");
  dofPtr = _dof_ptr;
  PetscFunctionReturn(0);
}


}
