/** \file LoopMethods.cpp
\brief methods for managing loops
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

#include <version.h>
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
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <UnknownInterface.hpp>
#include <LoopMethods.hpp>

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
PetscErrorCode KspMethod::copy_ksp(const KspMethod &ksp) {
  PetscFunctionBegin;
  this->ksp_ctx = ksp.ksp_ctx;
  this->ksp = ksp.ksp;
  this->ksp_f = ksp.ksp_f;
  this->ksp_A = ksp.ksp_A;
  this->ksp_B = ksp.ksp_B;
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
PetscErrorCode SnesMethod::copy_snes(const SnesMethod &snes) {
  PetscFunctionBegin;
  this->snes_ctx = snes.snes_ctx;
  this->snes = snes.snes;
  this->snes_x = snes.snes_x;
  this->snes_f = snes.snes_f;
  this->snes_A = snes.snes_A;
  this->snes_B = snes.snes_B;
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
PetscErrorCode TSMethod::copy_ts(const TSMethod &ts) {
  PetscFunctionBegin;
  this->ts_ctx = ts.ts_ctx;
  this->ts = ts.ts;
  this->ts_u = ts.ts_u;
  this->ts_u_t = ts.ts_u_t;
  this->ts_F = ts.ts_F;
  this->ts_A = ts.ts_A;
  this->ts_B = ts.ts_B;
  this->ts_step = ts.ts_step;
  this->ts_a = ts.ts_a;
  this->ts_t = ts.ts_t;
  PetscFunctionReturn(0);
}

//BasicMethod
BasicMethod::BasicMethod():
KspMethod(),
SnesMethod(),
TSMethod(),
nInTheLoop(0),
loopSize(0),
refinedEntitiesPtr(NULL),
refinedFiniteElementsPtr(NULL),
problemPtr(NULL),
fieldsPtr(NULL),
entitiesPtr(NULL),
dofsPtr(NULL),
finiteElementsPtr(NULL),
finiteElementsEntitiesPtr(NULL),
adjacenciesPtr(NULL) {
}

PetscErrorCode BasicMethod::copy_basic_method(const BasicMethod &basic) {
  PetscFunctionBegin;

  this->nInTheLoop = basic.nInTheLoop;
  this->loopSize = basic.loopSize;
  this->rAnk = basic.rAnk;
  this->sIze = basic.sIze;
  this->refinedEntitiesPtr = basic.refinedEntitiesPtr;
  this->refinedFiniteElementsPtr = basic.refinedFiniteElementsPtr;
  this->problemPtr = basic.problemPtr;
  this->fieldsPtr = basic.fieldsPtr;
  this->entitiesPtr = basic.entitiesPtr;
  this->dofsPtr = basic.dofsPtr;
  this->finiteElementsPtr = basic.finiteElementsPtr;
  this->finiteElementsEntitiesPtr = basic.finiteElementsEntitiesPtr;
  this->adjacenciesPtr = basic.adjacenciesPtr;

  PetscFunctionReturn(0);
}

//FEMethod
FEMethod::FEMethod():
BasicMethod() {
}

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

//EntMethod
EntMethod::EntMethod(): BasicMethod() {}

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
