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

namespace MoFEM {

PetscErrorCode TsSetIFunction(TS ts, PetscReal t, Vec u, Vec u_t, Vec F,
                              void *ctx) {
  MoFEMFunctionBegin;
  TsCtx *ts_ctx = static_cast<TsCtx *>(ctx);
  PetscLogEventBegin(ts_ctx->MOFEM_EVENT_TsCtxIFunction, 0, 0, 0, 0);
  CHKERR VecGhostUpdateBegin(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateBegin(u_t, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u_t, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR ts_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
      ts_ctx->problemName, COL, u, INSERT_VALUES, SCATTER_REVERSE);

  auto zero_ghost_vec = [](Vec g) {
    MoFEMFunctionBegin;
    Vec l;
    CHKERR VecGhostGetLocalForm(g, &l);
    double *a;
    CHKERR VecGetArray(l, &a);
    int s;
    CHKERR VecGetLocalSize(l, &s);
    for (int i = 0; i != s; ++i)
      a[i] = 0;
    CHKERR VecRestoreArray(l, &a);
    CHKERR VecGhostRestoreLocalForm(g, &l);
    MoFEMFunctionReturn(0);
  };
  CHKERR zero_ghost_vec(F);

  ts_ctx->vecAssembleSwitch = boost::movelib::make_unique<bool>(true);

  int step;
  CHKERR TSGetTimeStepNumber(ts, &step);
  // preprocess
  for (auto &bit : ts_ctx->preProcess_IFunction) {
    bit->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    bit->ts_u = u;
    bit->ts_u_t = u_t;
    bit->ts_F = F;
    bit->ts_t = t;
    bit->ts_step = step;
    CHKERR bit->setTsCtx(TSMethod::CTX_TSSETIFUNCTION);
    CHKERR bit->setTs(ts);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    CHKERR bit->setTsCtx(TSMethod::CTX_TSNONE);
    ts_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }

  // fe loops
  for (auto &lit : ts_ctx->loops_to_do_IFunction) {
    ts_ctx->vecAssembleSwitch = boost::move(lit.second->vecAssembleSwitch);
    lit.second->ts_u = u;
    lit.second->ts_u_t = u_t;
    lit.second->ts_F = F;
    lit.second->ts_t = t;
    lit.second->ts_step = step;
    CHKERR lit.second->setTsCtx(TSMethod::CTX_TSSETIFUNCTION);
    CHKERR lit.second->setTs(ts);
    CHKERR ts_ctx->mField.loop_finite_elements(
        ts_ctx->problemName, lit.first, *(lit.second), nullptr, ts_ctx->bH);
    CHKERR lit.second->setTsCtx(TSMethod::CTX_TSNONE);
    ts_ctx->vecAssembleSwitch = boost::move(lit.second->vecAssembleSwitch);
  }

  // post process
  for (auto &bit : ts_ctx->postProcess_IFunction) {
    bit->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    bit->ts_u = u;
    bit->ts_u_t = u_t;
    bit->ts_F = F;
    bit->ts_t = t;
    bit->ts_step = step;
    CHKERR bit->setTsCtx(TSMethod::CTX_TSSETIFUNCTION);
    CHKERR bit->setTs(ts);
    CHKERR ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,
                                                           *bit);
    CHKERR bit->setTsCtx(TSMethod::CTX_TSNONE);
    ts_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }

  if (*ts_ctx->vecAssembleSwitch) {
    CHKERR VecGhostUpdateBegin(F, ADD_VALUES, SCATTER_REVERSE);
    CHKERR VecGhostUpdateEnd(F, ADD_VALUES, SCATTER_REVERSE);
    CHKERR VecAssemblyBegin(F);
    CHKERR VecAssemblyEnd(F);
  }

  PetscLogEventEnd(ts_ctx->MOFEM_EVENT_TsCtxIFunction, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

PetscErrorCode TsSetIJacobian(TS ts, PetscReal t, Vec u, Vec u_t, PetscReal a,
                              Mat A, Mat B, void *ctx) {
  MoFEMFunctionBegin;

  TsCtx *ts_ctx = static_cast<TsCtx *>(ctx);
  PetscLogEventBegin(ts_ctx->MOFEM_EVENT_TsCtxIFunction, 0, 0, 0, 0);
  CHKERR VecGhostUpdateBegin(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateBegin(u_t, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u_t, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR ts_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
      ts_ctx->problemName, COL, u, INSERT_VALUES, SCATTER_REVERSE);
  if (ts_ctx->zeroMatrix) {
    CHKERR MatZeroEntries(B);
  }
  int step;
  CHKERR TSGetTimeStepNumber(ts, &step);

  ts_ctx->matAssembleSwitch =
      boost::movelib::make_unique<bool>(ts_ctx->zeroMatrix);

  // preproces
  for (auto &bit : ts_ctx->preProcess_IJacobian) {
    bit->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    bit->ts_u = u;
    bit->ts_u_t = u_t;
    bit->ts_A = A;
    bit->ts_B = B;
    bit->ts_t = t;
    bit->ts_a = a;
    bit->ts_step = step;
    CHKERR bit->setTsCtx(TSMethod::CTX_TSSETIJACOBIAN);
    CHKERR bit->setTs(ts);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    CHKERR bit->setTsCtx(TSMethod::CTX_TSNONE);
    ts_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

  for (auto &lit : ts_ctx->loops_to_do_IJacobian) {
    lit.second->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    lit.second->ts_u = u;
    lit.second->ts_u_t = u_t;
    lit.second->ts_A = A;
    lit.second->ts_B = B;
    lit.second->ts_t = t;
    lit.second->ts_a = a;
    lit.second->ts_step = step;
    CHKERR lit.second->setTsCtx(TSMethod::CTX_TSSETIJACOBIAN);
    CHKERR lit.second->setTs(ts);
    CHKERR ts_ctx->mField.loop_finite_elements(
        ts_ctx->problemName, lit.first, *(lit.second), nullptr, ts_ctx->bH);
    CHKERR lit.second->setTsCtx(TSMethod::CTX_TSNONE);
    ts_ctx->matAssembleSwitch = boost::move(lit.second->matAssembleSwitch);
  }

  // post process
  for (auto &bit : ts_ctx->postProcess_IJacobian) {
    bit->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    bit->ts_u = u;
    bit->ts_u_t = u_t;
    bit->ts_A = A;
    bit->ts_B = B;
    bit->ts_t = t;
    bit->ts_a = a;
    bit->ts_step = step;
    CHKERR bit->setTsCtx(TSMethod::CTX_TSSETIJACOBIAN);
    CHKERR bit->setTs(ts);
    CHKERR ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,
                                                           *bit);
    CHKERR bit->setTsCtx(TSMethod::CTX_TSNONE);
    ts_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

  if (ts_ctx->matAssembleSwitch) {
    CHKERR MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    CHKERR MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
  }
  PetscLogEventEnd(ts_ctx->MOFEM_EVENT_TsCtxIFunction, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

PetscErrorCode TsMonitorSet(TS ts, PetscInt step, PetscReal t, Vec u,
                            void *ctx) {
  MoFEMFunctionBegin;
  TsCtx *ts_ctx = static_cast<TsCtx *>(ctx);
  PetscLogEventBegin(ts_ctx->MOFEM_EVENT_TsCtxRHSFunction, 0, 0, 0, 0);
  CHKERR VecGhostUpdateBegin(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR ts_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
      ts_ctx->problemName, COL, u, INSERT_VALUES, SCATTER_REVERSE);

  // preproces
  for (auto &bit : ts_ctx->preProcess_Monitor) {
    bit->ts_u = u;
    bit->ts_t = t;
    bit->ts_step = step;
    bit->ts_F = PETSC_NULL;
    CHKERR bit->setTsCtx(TSMethod::CTX_TSSETIJACOBIAN);
    CHKERR bit->setTs(ts);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    CHKERR bit->setTsCtx(TSMethod::CTX_TSNONE);
  }

  for (auto &lit : ts_ctx->loops_to_do_Monitor) {
    lit.second->ts_u = u;
    lit.second->ts_t = t;
    lit.second->ts_step = step;
    lit.second->ts_F = PETSC_NULL;
    CHKERR lit.second->setTsCtx(TSMethod::CTX_TSTSMONITORSET);
    CHKERR lit.second->setTs(ts);
    CHKERR ts_ctx->mField.loop_finite_elements(
        ts_ctx->problemName, lit.first, *(lit.second), nullptr, ts_ctx->bH);
    CHKERR lit.second->setTsCtx(TSMethod::CTX_TSNONE);
  }

  // post process
  for (auto &bit : ts_ctx->postProcess_Monitor) {
    bit->ts_u = u;
    bit->ts_t = t;
    bit->ts_step = step;
    bit->ts_F = PETSC_NULL;
    CHKERR bit->setTsCtx(TSMethod::CTX_TSSETIJACOBIAN);
    CHKERR bit->setTs(ts);
    CHKERR ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,
                                                           *bit);
    CHKERR bit->setTsCtx(TSMethod::CTX_TSNONE);
  }
  PetscLogEventEnd(ts_ctx->MOFEM_EVENT_TsCtxRHSFunction, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

PetscErrorCode TSSetRHSFunction(TS ts, PetscReal t, Vec u, Vec F, void *ctx) {
  MoFEMFunctionBegin;
  TsCtx *ts_ctx = static_cast<TsCtx *>(ctx);
  PetscLogEventBegin(ts_ctx->MOFEM_EVENT_TsCtxRHSFunction, 0, 0, 0, 0);
  CHKERR VecGhostUpdateBegin(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR ts_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
      ts_ctx->problemName, COL, u, INSERT_VALUES, SCATTER_REVERSE);

  auto zero_ghost_vec = [](Vec g) {
    MoFEMFunctionBegin;
    Vec l;
    CHKERR VecGhostGetLocalForm(g, &l);
    double *a;
    CHKERR VecGetArray(l, &a);
    int s;
    CHKERR VecGetLocalSize(l, &s);
    for (int i = 0; i != s; ++i)
      a[i] = 0;
    CHKERR VecRestoreArray(l, &a);
    CHKERR VecGhostRestoreLocalForm(g, &l);
    MoFEMFunctionReturn(0);
  };
  CHKERR zero_ghost_vec(F);

  ts_ctx->vecAssembleSwitch = boost::movelib::make_unique<bool>(true);

  int step;
  CHKERR TSGetTimeStepNumber(ts, &step);
  for (auto &bit : ts_ctx->preProcess_RHSJacobian) {
    bit->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    bit->ts_u = u;
    bit->ts_F = F;
    bit->ts_t = t;
    bit->ts_step = step;
    CHKERR bit->setTsCtx(TSMethod::CTX_TSSETRHSFUNCTION);
    CHKERR bit->setTs(ts);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    CHKERR bit->setTsCtx(TSMethod::CTX_TSNONE);
    ts_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }

  // fe loops
  for (auto &lit : ts_ctx->loops_to_do_RHSFunction) {
    lit.second->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    lit.second->ts_u = u;
    lit.second->ts_F = F;
    lit.second->ts_t = t;
    lit.second->ts_step = step;
    CHKERR lit.second->setTsCtx(TSMethod::CTX_TSSETRHSFUNCTION);
    CHKERR lit.second->setTs(ts);
    CHKERR ts_ctx->mField.loop_finite_elements(
        ts_ctx->problemName, lit.first, *(lit.second), nullptr, ts_ctx->bH);
    CHKERR lit.second->setTsCtx(TSMethod::CTX_TSNONE);
    ts_ctx->vecAssembleSwitch = boost::move(lit.second->vecAssembleSwitch);
  }

  // post process
  for (auto &bit : ts_ctx->postProcess_RHSJacobian) {
    bit->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    bit->ts_u = u;
    bit->ts_F = F;
    bit->ts_t = t;
    bit->ts_step = step;
    CHKERR bit->setTsCtx(TSMethod::CTX_TSSETRHSFUNCTION);
    CHKERR bit->setTs(ts);
    CHKERR ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,
                                                           *bit);
    CHKERR bit->setTsCtx(TSMethod::CTX_TSNONE);
    ts_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }

  if (*ts_ctx->vecAssembleSwitch) {
    CHKERR VecGhostUpdateBegin(F, ADD_VALUES, SCATTER_REVERSE);
    CHKERR VecGhostUpdateEnd(F, ADD_VALUES, SCATTER_REVERSE);
    CHKERR VecAssemblyBegin(F);
    CHKERR VecAssemblyEnd(F);
  }

  PetscLogEventEnd(ts_ctx->MOFEM_EVENT_TsCtxRHSFunction, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

PetscErrorCode TSSetRHSJacobian(TS ts, PetscReal t, Vec u, Mat A, Mat B,
                                void *ctx) {
  MoFEMFunctionBegin;
  TsCtx *ts_ctx = static_cast<TsCtx *>(ctx);
  PetscLogEventBegin(ts_ctx->MOFEM_EVENT_TsCtxRHSJacobian, 0, 0, 0, 0);
  CHKERR VecGhostUpdateBegin(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR ts_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
      ts_ctx->problemName, COL, u, INSERT_VALUES, SCATTER_REVERSE);

  if (ts_ctx->zeroMatrix) {
    CHKERR MatZeroEntries(B);
  }

  ts_ctx->matAssembleSwitch =
      boost::movelib::make_unique<bool>(ts_ctx->zeroMatrix);

  int step;
  CHKERR TSGetTimeStepNumber(ts, &step);
  // preprocess
  for (auto &bit : ts_ctx->preProcess_RHSJacobian) {
    bit->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    bit->ts_u = u;
    bit->ts_A = A;
    bit->ts_B = B;
    bit->ts_t = t;
    bit->ts_step = step;
    CHKERR bit->setTsCtx(TSMethod::CTX_TSSETRHSJACOBIAN);
    CHKERR bit->setTs(ts);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    CHKERR bit->setTsCtx(TSMethod::CTX_TSNONE);
    ts_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

  // fe loops
  for (auto &lit : ts_ctx->loops_to_do_RHSJacobian) {
    lit.second->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    lit.second->ts_u = u;
    lit.second->ts_A = A;
    lit.second->ts_B = B;
    lit.second->ts_t = t;
    lit.second->ts_step = step;
    CHKERR lit.second->setTsCtx(TSMethod::CTX_TSSETRHSJACOBIAN);
    CHKERR lit.second->setTs(ts);
    CHKERR ts_ctx->mField.loop_finite_elements(
        ts_ctx->problemName, lit.first, *(lit.second), nullptr, ts_ctx->bH);
    CHKERR lit.second->setTsCtx(TSMethod::CTX_TSNONE);
    ts_ctx->matAssembleSwitch = boost::move(lit.second->matAssembleSwitch);
  }

  // post process
  for (auto &bit : ts_ctx->postProcess_RHSJacobian) {
    bit->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    bit->ts_u = u;
    bit->ts_A = A;
    bit->ts_B = B;
    bit->ts_t = t;
    bit->ts_step = step;
    CHKERR bit->setTsCtx(TSMethod::CTX_TSSETRHSJACOBIAN);
    CHKERR bit->setTs(ts);
    CHKERR ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,
                                                           *bit);
    CHKERR bit->setTsCtx(TSMethod::CTX_TSNONE);
    ts_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

  if (ts_ctx->matAssembleSwitch) {
    CHKERR MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    CHKERR MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
  }

  PetscLogEventEnd(ts_ctx->MOFEM_EVENT_TsCtxRHSJacobian, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
