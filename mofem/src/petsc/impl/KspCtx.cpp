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
#include <FEMMultiIndices.hpp>
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

#include <KspCtx.hpp>

namespace MoFEM {

PetscErrorCode KspRhs(KSP ksp,Vec f,void *ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  KspCtx* ksp_ctx = (KspCtx*)ctx;
  PetscLogEventBegin(ksp_ctx->USER_EVENT_KspRhs,0,0,0,0);
  KspCtx::basic_method_to_do::iterator bit = ksp_ctx->preProcess_Rhs.begin();
  for(;bit!=ksp_ctx->preProcess_Rhs.end();bit++) {
    ierr = (*bit)->set_ksp(ksp); CHKERRQ(ierr);
    (*bit)->ksp_f = f;
    ierr = (*bit)->set_ksp_ctx(KspMethod::CTX_SETFUNCTION);  CHKERRQ(ierr);
    ierr = ksp_ctx->mField.problem_basic_method_preProcess(ksp_ctx->problemName,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->set_ksp_ctx(KspMethod::CTX_KSPNONE);  CHKERRQ(ierr);
  }
  KspCtx::loops_to_do_type::iterator lit = ksp_ctx->loops_to_do_Rhs.begin();
  for(;lit!=ksp_ctx->loops_to_do_Rhs.end();lit++) {
    ierr = lit->second->set_ksp_ctx(KspMethod::CTX_SETFUNCTION);  CHKERRQ(ierr);
    ierr = lit->second->set_ksp(ksp); CHKERRQ(ierr);
    lit->second->ksp_f = f;
    ierr = ksp_ctx->mField.loop_finite_elements(ksp_ctx->problemName,lit->first,*(lit->second),ksp_ctx->bH);  CHKERRQ(ierr);
    ierr = lit->second->set_ksp_ctx(KspMethod::CTX_KSPNONE);  CHKERRQ(ierr);
  }
  bit = ksp_ctx->postProcess_Rhs.begin();
  for(;bit!=ksp_ctx->postProcess_Rhs.end();bit++) {
    ierr = (*bit)->set_ksp(ksp); CHKERRQ(ierr);
    (*bit)->ksp_f = f;
    ierr = (*bit)->set_ksp_ctx(KspMethod::CTX_SETFUNCTION);  CHKERRQ(ierr);
    ierr = ksp_ctx->mField.problem_basic_method_postProcess(ksp_ctx->problemName,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->set_ksp_ctx(KspMethod::CTX_KSPNONE);  CHKERRQ(ierr);
  }
  PetscLogEventEnd(ksp_ctx->USER_EVENT_KspRhs,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode KspMat(KSP ksp,Mat A,Mat B,void *ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  KspCtx* ksp_ctx = (KspCtx*)ctx;
  PetscLogEventBegin(ksp_ctx->USER_EVENT_KspMat,0,0,0,0);
  KspCtx::basic_method_to_do::iterator bit = ksp_ctx->preProcess_Mat.begin();
  for(;bit!=ksp_ctx->preProcess_Mat.end();bit++) {
    ierr = (*bit)->set_ksp(ksp); CHKERRQ(ierr);
    (*bit)->ksp_A = A;
    (*bit)->ksp_B = B;
    ierr = (*bit)->set_ksp_ctx(KspMethod::CTX_OPERATORS); CHKERRQ(ierr);
    ierr = ksp_ctx->mField.problem_basic_method_preProcess(ksp_ctx->problemName,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->set_ksp_ctx(KspMethod::CTX_KSPNONE);  CHKERRQ(ierr);
  }
  KspCtx::loops_to_do_type::iterator lit = ksp_ctx->loops_to_do_Mat.begin();
  for(;lit!=ksp_ctx->loops_to_do_Mat.end();lit++) {
    ierr = lit->second->set_ksp_ctx(KspMethod::CTX_OPERATORS); CHKERRQ(ierr);
    ierr = lit->second->set_ksp(ksp); CHKERRQ(ierr);
    lit->second->ksp_A = A;
    lit->second->ksp_B = B;
    ierr = ksp_ctx->mField.loop_finite_elements(ksp_ctx->problemName,lit->first,*(lit->second),ksp_ctx->bH);  CHKERRQ(ierr);
    ierr = lit->second->set_ksp_ctx(KspMethod::CTX_KSPNONE);
  }
  bit = ksp_ctx->postProcess_Mat.begin();
  for(;bit!=ksp_ctx->postProcess_Mat.end();bit++) {
    ierr = (*bit)->set_ksp(ksp); CHKERRQ(ierr);
    (*bit)->ksp_A = A;
    (*bit)->ksp_B = B;
    ierr = (*bit)->set_ksp_ctx(KspMethod::CTX_OPERATORS); CHKERRQ(ierr);
    ierr = ksp_ctx->mField.problem_basic_method_postProcess(ksp_ctx->problemName,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->set_ksp_ctx(KspMethod::CTX_KSPNONE);  CHKERRQ(ierr);
  }
  PetscLogEventEnd(ksp_ctx->USER_EVENT_KspMat,0,0,0,0);
  PetscFunctionReturn(0);
}

}
