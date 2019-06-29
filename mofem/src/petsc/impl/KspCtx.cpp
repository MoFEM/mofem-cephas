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

  ksp_ctx->vecAssembleSwitch = boost::movelib::make_unique<bool>(true);

  for (auto &bit : ksp_ctx->preProcess_Rhs) {
    bit->vecAssembleSwitch = boost::move(ksp_ctx->vecAssembleSwitch);
    CHKERR bit->setKsp(ksp);
    bit->ksp_f = f;
    CHKERR bit->setKspCtx(KspMethod::CTX_SETFUNCTION);
    CHKERR ksp_ctx->mField.problem_basic_method_preProcess(ksp_ctx->problemName,
                                                           *bit);
    CHKERR bit->setKspCtx(KspMethod::CTX_KSPNONE);
    ksp_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }
  for (auto &lit : ksp_ctx->loops_to_do_Rhs) {
    lit.second->vecAssembleSwitch = boost::move(ksp_ctx->vecAssembleSwitch);
    CHKERR lit.second->setKspCtx(KspMethod::CTX_SETFUNCTION);
    CHKERR lit.second->setKsp(ksp);
    lit.second->ksp_f = f;
    CHKERR ksp_ctx->mField.loop_finite_elements(
        ksp_ctx->problemName, lit.first, *(lit.second), nullptr, ksp_ctx->bH);
    CHKERR lit.second->setKspCtx(KspMethod::CTX_KSPNONE);
    ksp_ctx->vecAssembleSwitch = boost::move(lit.second->vecAssembleSwitch);
  }
  for (auto &bit : ksp_ctx->postProcess_Rhs) {
    bit->vecAssembleSwitch = boost::move(ksp_ctx->vecAssembleSwitch);
    CHKERR bit->setKsp(ksp);
    bit->ksp_f = f;
    CHKERR bit->setKspCtx(KspMethod::CTX_SETFUNCTION);
    CHKERR ksp_ctx->mField.problem_basic_method_postProcess(
        ksp_ctx->problemName, *bit);
    CHKERR bit->setKspCtx(KspMethod::CTX_KSPNONE);
    ksp_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }

  if (*ksp_ctx->vecAssembleSwitch) {
    CHKERR VecGhostUpdateBegin(f, ADD_VALUES, SCATTER_REVERSE);
    CHKERR VecGhostUpdateEnd(f, ADD_VALUES, SCATTER_REVERSE);
    CHKERR VecAssemblyBegin(f);
    CHKERR VecAssemblyEnd(f);
  }
  PetscLogEventEnd(ksp_ctx->MOFEM_EVENT_KspRhs, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}
PetscErrorCode KspMat(KSP ksp, Mat A, Mat B, void *ctx) {
  // PetscValidHeaderSpecific(ksp,KSP_CLASSID,1);
  MoFEMFunctionBegin;
  KspCtx *ksp_ctx = static_cast<KspCtx *>(ctx);
  PetscLogEventBegin(ksp_ctx->MOFEM_EVENT_KspMat, 0, 0, 0, 0);

  ksp_ctx->matAssembleSwitch = boost::movelib::make_unique<bool>(true);

  for (auto &bit : ksp_ctx->preProcess_Mat) {
    bit->matAssembleSwitch = boost::move(ksp_ctx->matAssembleSwitch);
    CHKERR bit->setKsp(ksp);
    bit->ksp_A = A;
    bit->ksp_B = B;
    CHKERR bit->setKspCtx(KspMethod::CTX_OPERATORS);
    CHKERR ksp_ctx->mField.problem_basic_method_preProcess(ksp_ctx->problemName,
                                                           *bit);
    CHKERR bit->setKspCtx(KspMethod::CTX_KSPNONE);
    ksp_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }
  for (auto &lit : ksp_ctx->loops_to_do_Mat) {
    lit.second->matAssembleSwitch = boost::move(ksp_ctx->matAssembleSwitch);
    lit.second->ksp_A = A;
    lit.second->ksp_B = B;
    CHKERR lit.second->setKsp(ksp);
    CHKERR lit.second->setKspCtx(KspMethod::CTX_OPERATORS);
    CHKERR ksp_ctx->mField.loop_finite_elements(
        ksp_ctx->problemName, lit.first, *(lit.second), nullptr, ksp_ctx->bH);
    CHKERR lit.second->setKspCtx(KspMethod::CTX_KSPNONE);
    ksp_ctx->matAssembleSwitch = boost::move(lit.second->matAssembleSwitch);
  }
  for (auto &bit : ksp_ctx->postProcess_Mat) {
    bit->matAssembleSwitch = boost::move(ksp_ctx->matAssembleSwitch);
    CHKERR bit->setKsp(ksp);
    bit->ksp_A = A;
    bit->ksp_B = B;
    CHKERR bit->setKspCtx(KspMethod::CTX_OPERATORS);
    CHKERR ksp_ctx->mField.problem_basic_method_postProcess(
        ksp_ctx->problemName, *bit);
    CHKERR bit->setKspCtx(KspMethod::CTX_KSPNONE);
    ksp_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

  if (ksp_ctx->matAssembleSwitch) {
    CHKERR MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    CHKERR MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
  }
  PetscLogEventEnd(ksp_ctx->MOFEM_EVENT_KspMat, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
