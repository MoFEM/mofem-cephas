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

// #if PETSC_VERSION_GE(3,6,0)
//   #include <petsc/private/kspimpl.h>
// #else
//   #include <petsc-private/kspimpl.h>
// #endif

namespace MoFEM {

PetscErrorCode KspRhs(KSP ksp, Vec f, void *ctx) {
  // PetscValidHeaderSpecific(ksp,KSP_CLASSID,1);
  MoFEMFunctionBegin;
  KspCtx *ksp_ctx = static_cast<KspCtx *>(ctx);
  PetscLogEventBegin(ksp_ctx->MOFEM_EVENT_KspRhs, 0, 0, 0, 0);
  KspCtx::BasicMethodsSequence::iterator bit = ksp_ctx->preProcess_Rhs.begin();
  for (; bit != ksp_ctx->preProcess_Rhs.end(); bit++) {
    CHKERR(*bit)->setKsp(ksp);
    (*bit)->ksp_f = f;
    CHKERR(*bit)->setKspCtx(KspMethod::CTX_SETFUNCTION);
    CHKERR ksp_ctx->mField.problem_basic_method_preProcess(ksp_ctx->problemName,
                                                           *(*(bit)));
    CHKERR(*bit)->setKspCtx(KspMethod::CTX_KSPNONE);
  }
  KspCtx::FEMethodsSequence::iterator lit = ksp_ctx->loops_to_do_Rhs.begin();
  for (; lit != ksp_ctx->loops_to_do_Rhs.end(); lit++) {
    CHKERR lit->second->setKspCtx(KspMethod::CTX_SETFUNCTION);
    CHKERR lit->second->setKsp(ksp);
    lit->second->ksp_f = f;
    CHKERR ksp_ctx->mField.loop_finite_elements(
        ksp_ctx->problemName, lit->first, *(lit->second), nullptr, ksp_ctx->bH);
    CHKERR lit->second->setKspCtx(KspMethod::CTX_KSPNONE);
  }
  bit = ksp_ctx->postProcess_Rhs.begin();
  for (; bit != ksp_ctx->postProcess_Rhs.end(); bit++) {
    CHKERR(*bit)->setKsp(ksp);
    (*bit)->ksp_f = f;
    CHKERR(*bit)->setKspCtx(KspMethod::CTX_SETFUNCTION);
    CHKERR ksp_ctx->mField.problem_basic_method_postProcess(
        ksp_ctx->problemName, *(*(bit)));
    CHKERR(*bit)->setKspCtx(KspMethod::CTX_KSPNONE);
  }
  CHKERR VecGhostUpdateBegin(f, ADD_VALUES, SCATTER_REVERSE);
  CHKERR VecGhostUpdateEnd(f, ADD_VALUES, SCATTER_REVERSE);
  CHKERR VecAssemblyBegin(f);
  CHKERR VecAssemblyEnd(f);
  PetscLogEventEnd(ksp_ctx->MOFEM_EVENT_KspRhs, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}
PetscErrorCode KspMat(KSP ksp, Mat A, Mat B, void *ctx) {
  // PetscValidHeaderSpecific(ksp,KSP_CLASSID,1);
  MoFEMFunctionBegin;
  KspCtx *ksp_ctx = static_cast<KspCtx *>(ctx);
  PetscLogEventBegin(ksp_ctx->MOFEM_EVENT_KspMat, 0, 0, 0, 0);
  KspCtx::BasicMethodsSequence::iterator bit = ksp_ctx->preProcess_Mat.begin();
  for (; bit != ksp_ctx->preProcess_Mat.end(); bit++) {
    CHKERR(*bit)->setKsp(ksp);
    (*bit)->ksp_A = A;
    (*bit)->ksp_B = B;
    CHKERR(*bit)->setKspCtx(KspMethod::CTX_OPERATORS);
    CHKERR ksp_ctx->mField.problem_basic_method_preProcess(ksp_ctx->problemName,
                                                           *(*(bit)));
    CHKERR(*bit)->setKspCtx(KspMethod::CTX_KSPNONE);
  }
  KspCtx::FEMethodsSequence::iterator lit = ksp_ctx->loops_to_do_Mat.begin();
  for (; lit != ksp_ctx->loops_to_do_Mat.end(); lit++) {
    lit->second->ksp_A = A;
    lit->second->ksp_B = B;
    CHKERR lit->second->setKsp(ksp);
    CHKERR lit->second->setKspCtx(KspMethod::CTX_OPERATORS);
    CHKERR ksp_ctx->mField.loop_finite_elements(
        ksp_ctx->problemName, lit->first, *(lit->second), nullptr, ksp_ctx->bH);
    CHKERR lit->second->setKspCtx(KspMethod::CTX_KSPNONE);
  }
  bit = ksp_ctx->postProcess_Mat.begin();
  for (; bit != ksp_ctx->postProcess_Mat.end(); bit++) {
    CHKERR(*bit)->setKsp(ksp);
    (*bit)->ksp_A = A;
    (*bit)->ksp_B = B;
    CHKERR(*bit)->setKspCtx(KspMethod::CTX_OPERATORS);
    CHKERR ksp_ctx->mField.problem_basic_method_postProcess(
        ksp_ctx->problemName, *(*(bit)));
    CHKERR(*bit)->setKspCtx(KspMethod::CTX_KSPNONE);
  }
  CHKERR MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
  CHKERR MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
  // MatView(A,PETSC_VIEWER_DRAW_WORLD);
  // std::string wait;
  // std::cin >> wait;
  PetscLogEventEnd(ksp_ctx->MOFEM_EVENT_KspMat, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
