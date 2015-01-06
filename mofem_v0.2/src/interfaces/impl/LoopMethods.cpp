
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

#include <petscsys.h>
#include <petscvec.h> 
#include <petscmat.h> 
#include <petscsnes.h> 
#include <petscts.h> 

#include <definitions.h>
#include <Common.hpp>

#include<LoopMethods.hpp>

namespace MoFEM {

//KSP
PetscErrorCode KspMethod::set_ksp_ctx(const KSPContext ctx_) {
  PetscFunctionBegin;
  ksp_ctx = ctx_;
  PetscFunctionReturn(0);
}
PetscErrorCode KspMethod::set_ksp(KSP ksp_) {
  PetscFunctionBegin;
  ksp = ksp_;
  PetscFunctionReturn(0);
}
//SNES
PetscErrorCode SnesMethod::set_snes_ctx(const SNESContext ctx_) {
  PetscFunctionBegin;
  snes_ctx = ctx_;
  PetscFunctionReturn(0);
}
PetscErrorCode SnesMethod::set_snes(SNES _snes) { 
  PetscFunctionBegin;
  snes = _snes;
  PetscFunctionReturn(0);
}
//TS
PetscErrorCode TSMethod::set_ts_ctx(const TSContext ctx_) {
  PetscFunctionBegin;
  ts_ctx = ctx_;
  PetscFunctionReturn(0);
}
PetscErrorCode TSMethod::set_ts(TS _ts) {
  PetscFunctionBegin;
  ts = _ts;
  PetscFunctionReturn(0);
}
//BasicMethod
BasicMethod::BasicMethod():
  refinedEntitiesPtr(NULL),refinedFiniteElementsPtr(NULL),
  problemPtr(NULL),fieldsPtr(NULL),entitiesPtr(NULL),dofsPtr(NULL),
  finiteElementsPtr(NULL),finiteElementsEntitiesPtr(NULL),adjacenciesPtr(NULL) {};
FEMethod::FEMethod(): BasicMethod(),
  fePtr(NULL),dataPtr(NULL),
  rowPtr(NULL),colPtr(NULL) {}

PetscErrorCode FEMethod::preProcess() {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class (preProcess)");
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod::postProcess() {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class (postProcess)");
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod::operator()() {   
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class (operator)");
  PetscFunctionReturn(0);
}
EntMethod::EntMethod(): BasicMethod(),dofPtr(NULL),dofNumeredPtr(NULL) {}
PetscErrorCode EntMethod::preProcess() {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class");
  PetscFunctionReturn(0);
}
PetscErrorCode EntMethod::postProcess() {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class");
  PetscFunctionReturn(0);
}
PetscErrorCode EntMethod::operator()() {   
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"should be implemented by user in derived class");
  PetscFunctionReturn(0);
}

}
