/* This file is part of MoFEM.
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <UnknownInterface.hpp>

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

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <AuxPETSc.hpp>
#include <TsCtx.hpp>

// #if PETSC_VERSION_GE(3,6,0)
//   #include <petsc/private/tsimpl.h>
// #else
//   #include <petsc-private/tsimpl.h>
// #endif

namespace MoFEM {

PetscErrorCode f_TSSetIFunction(TS ts,PetscReal t,Vec u,Vec u_t,Vec F,void *ctx) {
  // PetscValidHeaderSpecific(ts,TS_CLASSID,1);
  MoFEMFunctionBeginHot;
  TsCtx* ts_ctx = (TsCtx*)ctx;
  PetscLogEventBegin(ts_ctx->MOFEM_EVENT_TsCtxIFunction,0,0,0,0);
  ierr = VecGhostUpdateBegin(u,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(u,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecGhostUpdateBegin(u_t,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(u_t,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = ts_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
    ts_ctx->problemName,COL,u,INSERT_VALUES,SCATTER_REVERSE
  ); CHKERRG(ierr);
  ierr = VecZeroEntries(F); CHKERRG(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  int step;
  ierr = TSGetTimeStepNumber(ts,&step); CHKERRG(ierr);
  //preprocess
  TsCtx::BasicMethodsSequence::iterator bit = ts_ctx->preProcess_IFunction.begin();
  for(;bit!=ts_ctx->preProcess_IFunction.end();bit++) {
    (*bit)->ts_u = u;
    (*bit)->ts_u_t = u_t;
    (*bit)->ts_F = F;
    (*bit)->ts_t = t;
    (*bit)->ts_step = step;
    ierr = (*bit)->setTsCtx(TSMethod::CTX_TSSETIFUNCTION);
    ierr = (*bit)->setTs(ts); CHKERRG(ierr);
    ierr = ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,*(*(bit)));  CHKERRG(ierr);
    ierr = (*bit)->setTsCtx(TSMethod::CTX_TSNONE);
  }
  //fe loops
  TsCtx::FEMethodsSequence::iterator lit = ts_ctx->loops_to_do_IFunction.begin();
  for(;lit!=ts_ctx->loops_to_do_IFunction.end();lit++) {
    lit->second->ts_u = u;
    lit->second->ts_u_t = u_t;
    lit->second->ts_F = F;
    lit->second->ts_t = t;
    lit->second->ts_step = step;
    ierr = lit->second->setTsCtx(TSMethod::CTX_TSSETIFUNCTION);
    ierr = lit->second->setTs(ts); CHKERRG(ierr);
    ierr = ts_ctx->mField.loop_finite_elements(ts_ctx->problemName,lit->first,*(lit->second),ts_ctx->bH); CHKERRG(ierr);
    ierr = lit->second->setTsCtx(TSMethod::CTX_TSNONE);
  }
  //post process
  bit = ts_ctx->postProcess_IFunction.begin();
  for(;bit!=ts_ctx->postProcess_IFunction.end();bit++) {
    (*bit)->ts_u = u;
    (*bit)->ts_u_t = u_t;
    (*bit)->ts_F = F;
    (*bit)->ts_t = t;
    (*bit)->ts_step = step;
    ierr = (*bit)->setTsCtx(TSMethod::CTX_TSSETIFUNCTION);
    ierr = (*bit)->setTs(ts); CHKERRG(ierr);
    ierr = ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,*(*(bit)));  CHKERRG(ierr);
    ierr = (*bit)->setTsCtx(TSMethod::CTX_TSNONE);
  }
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  ierr = VecAssemblyBegin(F); CHKERRG(ierr);
  ierr = VecAssemblyEnd(F); CHKERRG(ierr);
  PetscLogEventEnd(ts_ctx->MOFEM_EVENT_TsCtxIFunction,0,0,0,0);
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode f_TSSetIJacobian(TS ts,PetscReal t,Vec u,Vec u_t,PetscReal a,Mat A,Mat B,void *ctx) {
  // PetscValidHeaderSpecific(ts,TS_CLASSID,1);
  MoFEMFunctionBeginHot;
  TsCtx* ts_ctx = (TsCtx*)ctx;
  PetscLogEventBegin(ts_ctx->MOFEM_EVENT_TsCtxIFunction,0,0,0,0);
  ierr = VecGhostUpdateBegin(u,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(u,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecGhostUpdateBegin(u_t,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(u_t,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = ts_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
    ts_ctx->problemName,COL,u,INSERT_VALUES,SCATTER_REVERSE
  ); CHKERRG(ierr);
  if(ts_ctx->zeroMatrix) {
    ierr = MatZeroEntries(B); CHKERRG(ierr);
  }
  int step;
  ierr = TSGetTimeStepNumber(ts,&step); CHKERRG(ierr);
  //preproces
  TsCtx::BasicMethodsSequence::iterator bit = ts_ctx->preProcess_IJacobian.begin();
  for(;bit!=ts_ctx->preProcess_IJacobian.end();bit++) {
    (*bit)->ts_u = u;
    (*bit)->ts_u_t = u_t;
    (*bit)->ts_A = A;
    (*bit)->ts_B = B;
    (*bit)->ts_t = t;
    (*bit)->ts_a = a;
    (*bit)->ts_step = step;
    ierr = (*bit)->setTsCtx(TSMethod::CTX_TSSETIJACOBIAN);
    ierr = (*bit)->setTs(ts); CHKERRG(ierr);
    ierr = ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,*(*(bit))); CHKERRG(ierr);
    ierr = (*bit)->setTsCtx(TSMethod::CTX_TSNONE); CHKERRG(ierr); CHKERRG(ierr);
  }
  TsCtx::FEMethodsSequence::iterator lit = ts_ctx->loops_to_do_IJacobian.begin();
  for(;lit!=ts_ctx->loops_to_do_IJacobian.end();lit++) {
    lit->second->ts_u = u;
    lit->second->ts_u_t = u_t;
    lit->second->ts_A = A;
    lit->second->ts_B = B;
    lit->second->ts_t = t;
    lit->second->ts_a = a;
    lit->second->ts_step = step;
    ierr = lit->second->setTsCtx(TSMethod::CTX_TSSETIJACOBIAN);
    ierr = lit->second->setTs(ts); CHKERRG(ierr);
    ierr = ts_ctx->mField.loop_finite_elements(ts_ctx->problemName,lit->first,*(lit->second),ts_ctx->bH); CHKERRG(ierr);
    ierr = lit->second->setTsCtx(TSMethod::CTX_TSNONE); CHKERRG(ierr);
  }
  //post process
  bit = ts_ctx->postProcess_IJacobian.begin();
  for(;bit!=ts_ctx->postProcess_IJacobian.end();bit++) {
    (*bit)->ts_u = u;
    (*bit)->ts_u_t = u_t;
    (*bit)->ts_A = A;
    (*bit)->ts_B = B;
    (*bit)->ts_t = t;
    (*bit)->ts_a = a;
    (*bit)->ts_step = step;
    ierr = (*bit)->setTsCtx(TSMethod::CTX_TSSETIJACOBIAN);
    ierr = (*bit)->setTs(ts); CHKERRG(ierr);
    ierr = ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,*(*(bit))); CHKERRG(ierr);
    ierr = (*bit)->setTsCtx(TSMethod::CTX_TSNONE); CHKERRG(ierr);
  }
  if(ts_ctx->zeroMatrix) {
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY); CHKERRG(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY); CHKERRG(ierr);
  }
  PetscLogEventEnd(ts_ctx->MOFEM_EVENT_TsCtxIFunction,0,0,0,0);
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode f_TSMonitorSet(TS ts,PetscInt step,PetscReal t,Vec u,void *ctx) {
  // PetscValidHeaderSpecific(ts,TS_CLASSID,1);
  MoFEMFunctionBeginHot;
  TsCtx* ts_ctx = (TsCtx*)ctx;
  PetscLogEventBegin(ts_ctx->MOFEM_EVENT_TsCtxRHSFunction,0,0,0,0);
  ierr = VecGhostUpdateBegin(u,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(u,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = ts_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
    ts_ctx->problemName,COL,u,INSERT_VALUES,SCATTER_REVERSE
  ); CHKERRG(ierr);
  //preproces
  TsCtx::BasicMethodsSequence::iterator bit = ts_ctx->preProcess_Monitor.begin();
  for(;bit!=ts_ctx->preProcess_Monitor.end();bit++) {
    (*bit)->ts_u = u;
    (*bit)->ts_t = t;
    (*bit)->ts_step = step;
    (*bit)->ts_F = PETSC_NULL;
    ierr = (*bit)->setTsCtx(TSMethod::CTX_TSSETIJACOBIAN);
    ierr = (*bit)->setTs(ts); CHKERRG(ierr);
    ierr = ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,*(*(bit))); CHKERRG(ierr);
    ierr = (*bit)->setTsCtx(TSMethod::CTX_TSNONE); CHKERRG(ierr); CHKERRG(ierr);
  }
  TsCtx::FEMethodsSequence::iterator lit = ts_ctx->loops_to_do_Monitor.begin();
  for(;lit!=ts_ctx->loops_to_do_Monitor.end();lit++) {
    lit->second->ts_u = u;
    lit->second->ts_t = t;
    lit->second->ts_step = step;
    lit->second->ts_F = PETSC_NULL;
    ierr = lit->second->setTsCtx(TSMethod::CTX_TSTSMONITORSET);
    ierr = lit->second->setTs(ts); CHKERRG(ierr);
    ierr = ts_ctx->mField.loop_finite_elements(ts_ctx->problemName,lit->first,*(lit->second),ts_ctx->bH); CHKERRG(ierr);
    ierr = lit->second->setTsCtx(TSMethod::CTX_TSNONE);
  }
  //post process
  bit = ts_ctx->postProcess_Monitor.begin();
  for(;bit!=ts_ctx->postProcess_Monitor.end();bit++) {
    (*bit)->ts_u = u;
    (*bit)->ts_t = t;
    (*bit)->ts_step = step;
    (*bit)->ts_F = PETSC_NULL;
    ierr = (*bit)->setTsCtx(TSMethod::CTX_TSSETIJACOBIAN);
    ierr = (*bit)->setTs(ts); CHKERRG(ierr);
    ierr = ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,*(*(bit))); CHKERRG(ierr);
    ierr = (*bit)->setTsCtx(TSMethod::CTX_TSNONE); CHKERRG(ierr);
  }
  PetscLogEventEnd(ts_ctx->MOFEM_EVENT_TsCtxRHSFunction,0,0,0,0);
  MoFEMFunctionReturnHot(0);
}


}
