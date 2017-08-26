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
#include <Core.hpp>

#include <AuxPTESc.hpp>
#include <KspCtx.hpp>

#if PETSC_VERSION_GE(3,6,0)
  #include <petsc/private/kspimpl.h>
#else
  #include <petsc-private/kspimpl.h>
#endif

namespace MoFEM {

PetscErrorCode KspRhs(KSP ksp,Vec f,void *ctx) {
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,1);
  PetscFunctionBegin;
  KspCtx* ksp_ctx = (KspCtx*)ctx;
  PetscLogEventBegin(ksp_ctx->USER_EVENT_KspRhs,0,0,0,0);
  KspCtx::BasicMethodsSequence::iterator bit = ksp_ctx->preProcess_Rhs.begin();
  for(;bit!=ksp_ctx->preProcess_Rhs.end();bit++) {
    ierr = (*bit)->setKsp(ksp); CHKERRQ(ierr);
    (*bit)->ksp_f = f;
    ierr = (*bit)->setKspCtx(KspMethod::CTX_SETFUNCTION);  CHKERRQ(ierr);
    ierr = ksp_ctx->mField.problem_basic_method_preProcess(ksp_ctx->problemName,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->setKspCtx(KspMethod::CTX_KSPNONE);  CHKERRQ(ierr);
  }
  KspCtx::FEMethodsSequence::iterator lit = ksp_ctx->loops_to_do_Rhs.begin();
  for(;lit!=ksp_ctx->loops_to_do_Rhs.end();lit++) {
    ierr = lit->second->setKspCtx(KspMethod::CTX_SETFUNCTION);  CHKERRQ(ierr);
    ierr = lit->second->setKsp(ksp); CHKERRQ(ierr);
    lit->second->ksp_f = f;
    ierr = ksp_ctx->mField.loop_finite_elements(ksp_ctx->problemName,lit->first,*(lit->second),ksp_ctx->bH);  CHKERRQ(ierr);
    ierr = lit->second->setKspCtx(KspMethod::CTX_KSPNONE);  CHKERRQ(ierr);
  }
  bit = ksp_ctx->postProcess_Rhs.begin();
  for(;bit!=ksp_ctx->postProcess_Rhs.end();bit++) {
    ierr = (*bit)->setKsp(ksp); CHKERRQ(ierr);
    (*bit)->ksp_f = f;
    ierr = (*bit)->setKspCtx(KspMethod::CTX_SETFUNCTION);  CHKERRQ(ierr);
    ierr = ksp_ctx->mField.problem_basic_method_postProcess(ksp_ctx->problemName,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->setKspCtx(KspMethod::CTX_KSPNONE);  CHKERRQ(ierr);
  }
  ierr = VecGhostUpdateBegin(f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(f); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(f); CHKERRQ(ierr);
  PetscLogEventEnd(ksp_ctx->USER_EVENT_KspRhs,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode KspMat(KSP ksp,Mat A,Mat B,void *ctx) {
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,1);
  PetscFunctionBegin;
  KspCtx* ksp_ctx = (KspCtx*)ctx;
  PetscLogEventBegin(ksp_ctx->USER_EVENT_KspMat,0,0,0,0);
  KspCtx::BasicMethodsSequence::iterator bit = ksp_ctx->preProcess_Mat.begin();
  for(;bit!=ksp_ctx->preProcess_Mat.end();bit++) {
    ierr = (*bit)->setKsp(ksp); CHKERRQ(ierr);
    (*bit)->ksp_A = A;
    (*bit)->ksp_B = B;
    ierr = (*bit)->setKspCtx(KspMethod::CTX_OPERATORS); CHKERRQ(ierr);
    ierr = ksp_ctx->mField.problem_basic_method_preProcess(ksp_ctx->problemName,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->setKspCtx(KspMethod::CTX_KSPNONE);  CHKERRQ(ierr);
  }
  KspCtx::FEMethodsSequence::iterator lit = ksp_ctx->loops_to_do_Mat.begin();
  for(;lit!=ksp_ctx->loops_to_do_Mat.end();lit++) {
    lit->second->ksp_A = A;
    lit->second->ksp_B = B;
    ierr = lit->second->setKsp(ksp); CHKERRQ(ierr);
    ierr = lit->second->setKspCtx(KspMethod::CTX_OPERATORS); CHKERRQ(ierr);
    ierr = ksp_ctx->mField.loop_finite_elements(ksp_ctx->problemName,lit->first,*(lit->second),ksp_ctx->bH);  CHKERRQ(ierr);
    ierr = lit->second->setKspCtx(KspMethod::CTX_KSPNONE);
  }
  bit = ksp_ctx->postProcess_Mat.begin();
  for(;bit!=ksp_ctx->postProcess_Mat.end();bit++) {
    ierr = (*bit)->setKsp(ksp); CHKERRQ(ierr);
    (*bit)->ksp_A = A;
    (*bit)->ksp_B = B;
    ierr = (*bit)->setKspCtx(KspMethod::CTX_OPERATORS); CHKERRQ(ierr);
    ierr = ksp_ctx->mField.problem_basic_method_postProcess(ksp_ctx->problemName,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->setKspCtx(KspMethod::CTX_KSPNONE);  CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  // MatView(A,PETSC_VIEWER_DRAW_WORLD);
  // std::string wait;
  // std::cin >> wait;
  PetscLogEventEnd(ksp_ctx->USER_EVENT_KspMat,0,0,0,0);
  PetscFunctionReturn(0);
}

}
