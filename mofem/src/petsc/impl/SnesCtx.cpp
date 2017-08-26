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
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <VecManager.hpp>
#include <AuxPTESc.hpp>
#include <SnesCtx.hpp>

namespace MoFEM {

PetscErrorCode SnesRhs(SNES snes,Vec x,Vec f,void *ctx) {
  PetscFunctionBegin;

  SnesCtx* snes_ctx = (SnesCtx*)ctx;
  PetscLogEventBegin(snes_ctx->USER_EVENT_SnesRhs,0,0,0,0);
  ierr = VecGhostUpdateBegin(x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = snes_ctx->mField.query_interface<VecManager>()->setLocalGhostVector(snes_ctx->problemName,COL,x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecZeroEntries(f); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  SnesCtx::BasicMethodsSequence::iterator bit = snes_ctx->preProcess_Rhs.begin();
  for(;bit!=snes_ctx->preProcess_Rhs.end();bit++) {
    ierr = (*bit)->setSnes(snes); CHKERRQ(ierr);
    (*bit)->snes_x = x;
    (*bit)->snes_f = f;
    ierr = (*bit)->setSnesCtx(SnesMethod::CTX_SNESSETFUNCTION);  CHKERRQ(ierr);
    ierr = snes_ctx->mField.problem_basic_method_preProcess(snes_ctx->problemName,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->setSnesCtx(SnesMethod::CTX_SNESNONE);  CHKERRQ(ierr);
  }
  SnesCtx::FEMethodsSequence::iterator lit = snes_ctx->loops_to_do_Rhs.begin();
  for(;lit!=snes_ctx->loops_to_do_Rhs.end();lit++) {
    ierr = lit->second->setSnesCtx(SnesMethod::CTX_SNESSETFUNCTION);  CHKERRQ(ierr);
    ierr = lit->second->setSnes(snes); CHKERRQ(ierr);
    lit->second->snes_x = x;
    lit->second->snes_f = f;
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\t\tLoop FE for Rhs: %s\n",lit->first.c_str());
    //PetscSynchronizedFlush(PETSC_COMM_WORLD);
    ierr = snes_ctx->mField.loop_finite_elements(snes_ctx->problemName,lit->first,*(lit->second),snes_ctx->bH);  CHKERRQ(ierr);
    ierr = lit->second->setSnesCtx(SnesMethod::CTX_SNESNONE);  CHKERRQ(ierr);
  }
  bit = snes_ctx->postProcess_Rhs.begin();
  for(;bit!=snes_ctx->postProcess_Rhs.end();bit++) {
    ierr = (*bit)->setSnes(snes); CHKERRQ(ierr);
    (*bit)->snes_x = x;
    (*bit)->snes_f = f;
    ierr = (*bit)->setSnesCtx(SnesMethod::CTX_SNESSETFUNCTION);  CHKERRQ(ierr);
    ierr = snes_ctx->mField.problem_basic_method_postProcess(snes_ctx->problemName,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->setSnesCtx(SnesMethod::CTX_SNESNONE);  CHKERRQ(ierr);
  }
  ierr = VecGhostUpdateBegin(f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(f); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(f); CHKERRQ(ierr);
  PetscLogEventEnd(snes_ctx->USER_EVENT_SnesRhs,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode SnesMat(SNES snes,Vec x,Mat A,Mat B,void *ctx) {
  PetscFunctionBegin;

  SnesCtx* snes_ctx = (SnesCtx*)ctx;
  PetscLogEventBegin(snes_ctx->USER_EVENT_SnesMat,0,0,0,0);
  if(snes_ctx->zeroPreCondMatrixB) {
    ierr = MatZeroEntries(B); CHKERRQ(ierr);
  }
  SnesCtx::BasicMethodsSequence::iterator bit = snes_ctx->preProcess_Mat.begin();
  for(;bit!=snes_ctx->preProcess_Mat.end();bit++) {
    ierr = (*bit)->setSnes(snes); CHKERRQ(ierr);
    (*bit)->snes_x = x;
    (*bit)->snes_A = A;
    (*bit)->snes_B = B;
    ierr = (*bit)->setSnesCtx(SnesMethod::CTX_SNESSETJACOBIAN); CHKERRQ(ierr);
    ierr = snes_ctx->mField.problem_basic_method_preProcess(snes_ctx->problemName,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->setSnesCtx(SnesMethod::CTX_SNESNONE);  CHKERRQ(ierr);
  }
  SnesCtx::FEMethodsSequence::iterator lit = snes_ctx->loops_to_do_Mat.begin();
  for(;lit!=snes_ctx->loops_to_do_Mat.end();lit++) {
    ierr = lit->second->setSnesCtx(SnesMethod::CTX_SNESSETJACOBIAN); CHKERRQ(ierr);
    ierr = lit->second->setSnes(snes); CHKERRQ(ierr);
    lit->second->snes_x = x;
    lit->second->snes_A = A;
    lit->second->snes_B = B;
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\t\tLoop FE for Mat: %s\n",lit->first.c_str());
    //PetscSynchronizedFlush(PETSC_COMM_WORLD);
    ierr = snes_ctx->mField.loop_finite_elements(snes_ctx->problemName,lit->first,*(lit->second),snes_ctx->bH);  CHKERRQ(ierr);
    ierr = lit->second->setSnesCtx(SnesMethod::CTX_SNESNONE);
  }
  bit = snes_ctx->postProcess_Mat.begin();
  for(;bit!=snes_ctx->postProcess_Mat.end();bit++) {
    ierr = (*bit)->setSnes(snes); CHKERRQ(ierr);
    (*bit)->snes_x = x;
    (*bit)->snes_A = A;
    (*bit)->snes_B = B;
    ierr = (*bit)->setSnesCtx(SnesMethod::CTX_SNESSETJACOBIAN); CHKERRQ(ierr);
    ierr = snes_ctx->mField.problem_basic_method_postProcess(snes_ctx->problemName,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->setSnesCtx(SnesMethod::CTX_SNESNONE);  CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  PetscLogEventEnd(snes_ctx->USER_EVENT_SnesMat,0,0,0,0);
  PetscFunctionReturn(0);
}

}
