/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

namespace MoFEM {

MoFEMErrorCode TsCtx::clearLoops() {
  MoFEMFunctionBeginHot;
  loopsIJacobian.clear();
  loopsIFunction.clear();
  loopsMonitor.clear();
  loopsRHSJacobian.clear();
  loopsRHSFunction.clear();
  preProcessIJacobian.clear();
  postProcessIJacobian.clear();
  preProcessIFunction.clear();
  postProcessIFunction.clear();
  preProcessMonitor.clear();
  postProcessMonitor.clear();
  preProcessRHSJacobian.clear();
  preProcessRHSFunction.clear();
  postProcessRHSJacobian.clear();
  postProcessRHSFunction.clear();
  MoFEMFunctionReturnHot(0);
}  

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
#if PETSC_VERSION_GE(3, 8, 0)
  CHKERR TSGetStepNumber(ts, &step);
#else
  CHKERR TSGetTimeStepNumber(ts, &step);
#endif

  auto set = [&](auto &fe) {
    fe.ts = ts;
    fe.ts_u = u;
    fe.ts_u_t = u_t;
    fe.ts_F = F;
    fe.ts_t = t;
    fe.ts_step = step;
    fe.ts_ctx = TSMethod::CTX_TSSETIFUNCTION;
    fe.snes_ctx = SnesMethod::CTX_SNESSETFUNCTION;
    fe.ksp_ctx = KspMethod::CTX_SETFUNCTION;
    fe.data_ctx = PetscData::CtxSetF | PetscData::CtxSetX |
                  PetscData::CtxSetX_T | PetscData::CtxSetTime;
  };

  auto unset = [&](auto &fe) {
    fe.ts_ctx = TSMethod::CTX_TSNONE;
    fe.snes_ctx = SnesMethod::CTX_SNESNONE;
    fe.ksp_ctx = KspMethod::CTX_KSPNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };

  // preprocess
  for (auto &bit : ts_ctx->preProcessIFunction) {
    bit->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    unset(*bit);
    ts_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }

  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

  // fe loops
  for (auto &lit : ts_ctx->loopsIFunction) {
    lit.second->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    set(*lit.second);
    CHKERR ts_ctx->mField.loop_finite_elements(ts_ctx->problemName, lit.first,
                                               *(lit.second), nullptr,
                                               ts_ctx->bH, cache_ptr);
    unset(*lit.second);
    ts_ctx->vecAssembleSwitch = boost::move(lit.second->vecAssembleSwitch);
  }

  // post process
  for (auto &bit : ts_ctx->postProcessIFunction) {
    bit->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,
                                                           *bit);
    unset(*bit);
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
#if PETSC_VERSION_GE(3, 8, 0)
  CHKERR TSGetStepNumber(ts, &step);
#else
  CHKERR TSGetTimeStepNumber(ts, &step);
#endif

  ts_ctx->matAssembleSwitch =
      boost::movelib::make_unique<bool>(ts_ctx->zeroMatrix);

  auto set = [&](auto &fe) {
    fe.ts = ts;
    fe.ts_u = u;
    fe.ts_u_t = u_t;
    fe.ts_A = A;
    fe.ts_B = B;
    fe.ts_t = t;
    fe.ts_a = a;
    fe.ts_step = step;
    fe.ts_ctx = TSMethod::CTX_TSSETIJACOBIAN;
    fe.snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
    fe.ksp_ctx = KspMethod::CTX_OPERATORS;
    fe.data_ctx = PetscData::CtxSetA | PetscData::CtxSetB | PetscData::CtxSetX |
                  PetscData::CtxSetX_T | PetscData::CtxSetTime;
  };

  auto unset = [&](auto &fe) {
    fe.ts_ctx = TSMethod::CTX_TSNONE;
    fe.snes_ctx = SnesMethod::CTX_SNESNONE;
    fe.ksp_ctx = KspMethod::CTX_KSPNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };

  // preproces
  for (auto &bit : ts_ctx->preProcessIJacobian) {
    bit->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    unset(*bit);
    ts_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

  for (auto &lit : ts_ctx->loopsIJacobian) {
    lit.second->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    set(*lit.second);
    CHKERR ts_ctx->mField.loop_finite_elements(ts_ctx->problemName, lit.first,
                                               *(lit.second), nullptr,
                                               ts_ctx->bH, cache_ptr);
    unset(*lit.second);
    ts_ctx->matAssembleSwitch = boost::move(lit.second->matAssembleSwitch);
  }

  // post process
  for (auto &bit : ts_ctx->postProcessIJacobian) {
    bit->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,
                                                           *bit);
    unset(*bit);
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
  PetscLogEventBegin(ts_ctx->MOFEM_EVENT_TsCtxMonitor, 0, 0, 0, 0);
  CHKERR VecGhostUpdateBegin(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR ts_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
      ts_ctx->problemName, COL, u, INSERT_VALUES, SCATTER_REVERSE);

  auto set = [&](auto &fe) {
    fe.ts = ts;
    fe.ts_u = u;
    fe.ts_t = t;
    fe.ts_step = step;
    fe.ts_F = PETSC_NULL;
    fe.ts_ctx = TSMethod::CTX_TSTSMONITORSET;
    fe.snes_ctx = SnesMethod::CTX_SNESNONE;
    fe.ksp_ctx = KspMethod::CTX_KSPNONE;
    fe.data_ctx = PetscData::CtxSetX | PetscData::CtxSetTime;
  };

  auto unset = [&](auto &fe) {
    fe.ts_ctx = TSMethod::CTX_TSNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };

  // preproces
  for (auto &bit : ts_ctx->preProcessMonitor) {
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    unset(*bit);
  }

  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

  for (auto &lit : ts_ctx->loopsMonitor) {
    set(*lit.second);
    CHKERR ts_ctx->mField.loop_finite_elements(ts_ctx->problemName, lit.first,
                                               *(lit.second), nullptr,
                                               ts_ctx->bH, cache_ptr);
    unset(*lit.second);
  }

  // post process
  for (auto &bit : ts_ctx->postProcessMonitor) {
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,
                                                           *bit);
    unset(*bit);
  }

  PetscLogEventEnd(ts_ctx->MOFEM_EVENT_TsCtxMonitor, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

PetscErrorCode TsSetRHSFunction(TS ts, PetscReal t, Vec u, Vec F, void *ctx) {
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
#if PETSC_VERSION_GE(3, 8, 0)
  CHKERR TSGetStepNumber(ts, &step);
#else
  CHKERR TSGetTimeStepNumber(ts, &step);
#endif

  auto set = [&](auto &fe) {
    fe.ts_u = u;
    fe.ts_F = F;
    fe.ts_t = t;
    fe.ts = ts;
    fe.ts_step = step;
    fe.ts_ctx = TSMethod::CTX_TSSETRHSFUNCTION;
    fe.snes_ctx = SnesMethod::CTX_SNESSETFUNCTION;
    fe.ksp_ctx = KspMethod::CTX_SETFUNCTION;
    fe.data_ctx =
        PetscData::CtxSetF | PetscData::CtxSetX | PetscData::CtxSetTime;
  };

  auto unset = [&](auto &fe) {
    fe.ts_ctx = TSMethod::CTX_TSNONE;
    fe.snes_ctx = SnesMethod::CTX_SNESNONE;
    fe.ksp_ctx = KspMethod::CTX_KSPNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };

  for (auto &bit : ts_ctx->preProcessRHSJacobian) {
    bit->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    unset(*bit);
    ts_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }

  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

  // fe loops
  for (auto &lit : ts_ctx->loopsRHSFunction) {
    lit.second->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    set(*lit.second);
    CHKERR ts_ctx->mField.loop_finite_elements(ts_ctx->problemName, lit.first,
                                               *(lit.second), nullptr,
                                               ts_ctx->bH, cache_ptr);
    unset(*lit.second);
    ts_ctx->vecAssembleSwitch = boost::move(lit.second->vecAssembleSwitch);
  }

  // post process
  for (auto &bit : ts_ctx->postProcessRHSJacobian) {
    bit->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,
                                                           *bit);
    unset(*bit);
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

PetscErrorCode TsSetRHSJacobian(TS ts, PetscReal t, Vec u, Mat A, Mat B,
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
#if PETSC_VERSION_GE(3, 8, 0)
  CHKERR TSGetStepNumber(ts, &step);
#else
  CHKERR TSGetTimeStepNumber(ts, &step);
#endif

  auto set = [&](auto &fe) {
    fe.ts_u = u;
    fe.ts_A = A;
    fe.ts_B = B;
    fe.ts_t = t;
    fe.ts_step = step;
    fe.ts_ctx = TSMethod::CTX_TSSETRHSJACOBIAN;
    fe.snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
    fe.ksp_ctx = KspMethod::CTX_OPERATORS;
    fe.data_ctx = PetscData::CtxSetA | PetscData::CtxSetB | PetscData::CtxSetX |
                  PetscData::CtxSetTime;
    fe.ts = ts;
  };

  auto unset = [&](auto &fe) {
    fe.ts_ctx = TSMethod::CTX_TSNONE;
    fe.snes_ctx = SnesMethod::CTX_SNESNONE;
    fe.ksp_ctx = KspMethod::CTX_KSPNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };

  // preprocess
  for (auto &bit : ts_ctx->preProcessRHSJacobian) {
    bit->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    unset(*bit);
    ts_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

  // fe loops
  for (auto &lit : ts_ctx->loopsRHSJacobian) {
    lit.second->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    set(*lit.second);
    CHKERR ts_ctx->mField.loop_finite_elements(ts_ctx->problemName, lit.first,
                                               *(lit.second), nullptr,
                                               ts_ctx->bH, cache_ptr);
    unset(*lit.second);
    ts_ctx->matAssembleSwitch = boost::move(lit.second->matAssembleSwitch);
  }

  // post process
  for (auto &bit : ts_ctx->postProcessRHSJacobian) {
    bit->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,
                                                           *bit);
    unset(*bit);
    ts_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

  if (ts_ctx->matAssembleSwitch) {
    CHKERR MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    CHKERR MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
  }

  PetscLogEventEnd(ts_ctx->MOFEM_EVENT_TsCtxRHSJacobian, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

PetscErrorCode TsSetI2Jacobian(TS ts, PetscReal t, Vec u, Vec u_t, Vec u_tt,
                               PetscReal a, PetscReal aa, Mat A, Mat B,
                               void *ctx) {
  MoFEMFunctionBegin;

  TsCtx *ts_ctx = static_cast<TsCtx *>(ctx);
  PetscLogEventBegin(ts_ctx->MOFEM_EVENT_TsCtxI2Function, 0, 0, 0, 0);
  CHKERR VecGhostUpdateBegin(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateBegin(u_t, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u_t, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateBegin(u_tt, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u_tt, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR ts_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
      ts_ctx->problemName, COL, u, INSERT_VALUES, SCATTER_REVERSE);
  if (ts_ctx->zeroMatrix) {
    CHKERR MatZeroEntries(B);
  }
  int step;
#if PETSC_VERSION_GE(3, 8, 0)
  CHKERR TSGetStepNumber(ts, &step);
#else
  CHKERR TSGetTimeStepNumber(ts, &step);
#endif

  ts_ctx->matAssembleSwitch =
      boost::movelib::make_unique<bool>(ts_ctx->zeroMatrix);

  auto set = [&](auto &fe) {
    fe.ts_u = u;
    fe.ts_u_t = u_t;
    fe.ts_u_tt = u_tt;
    fe.ts_A = A;
    fe.ts_B = B;
    fe.ts_t = t;
    fe.ts_a = a;
    fe.ts_aa = aa;
    fe.ts_step = step;

    fe.ts_ctx = TSMethod::CTX_TSSETIJACOBIAN;
    fe.snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
    fe.ksp_ctx = KspMethod::CTX_OPERATORS;
    fe.data_ctx = PetscData::CtxSetA | PetscData::CtxSetB | PetscData::CtxSetX |
                  PetscData::CtxSetX_T | PetscData::CtxSetX_TT |
                  PetscData::CtxSetTime;
    fe.ts = ts;
  };

  auto unset = [&](auto &fe) {
    fe.ts_ctx = TSMethod::CTX_TSNONE;
    fe.snes_ctx = SnesMethod::CTX_SNESNONE;
    fe.ksp_ctx = KspMethod::CTX_KSPNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };

  // preproces
  for (auto &bit : ts_ctx->preProcessIJacobian) {
    bit->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    unset(*bit);
    ts_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

  for (auto &lit : ts_ctx->loopsIJacobian) {
    lit.second->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    set(*lit.second);
    CHKERR ts_ctx->mField.loop_finite_elements(ts_ctx->problemName, lit.first,
                                               *(lit.second), nullptr,
                                               ts_ctx->bH, cache_ptr);
    unset(*lit.second);
    ts_ctx->matAssembleSwitch = boost::move(lit.second->matAssembleSwitch);
  }

  // post process
  for (auto &bit : ts_ctx->postProcessIJacobian) {
    bit->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,
                                                           *bit);
    unset(*bit);
    ts_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

  if (ts_ctx->matAssembleSwitch) {
    CHKERR MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    CHKERR MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
  }
  PetscLogEventEnd(ts_ctx->MOFEM_EVENT_TsCtxI2Function, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

PetscErrorCode TsSetI2Function(TS ts, PetscReal t, Vec u, Vec u_t, Vec u_tt,
                               Vec F, void *ctx) {
  MoFEMFunctionBegin;
  TsCtx *ts_ctx = static_cast<TsCtx *>(ctx);
  PetscLogEventBegin(ts_ctx->MOFEM_EVENT_TsCtxIFunction, 0, 0, 0, 0);
  CHKERR VecGhostUpdateBegin(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateBegin(u_t, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u_t, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateBegin(u_tt, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(u_tt, INSERT_VALUES, SCATTER_FORWARD);
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
#if PETSC_VERSION_GE(3, 8, 0)
  CHKERR TSGetStepNumber(ts, &step);
#else
  CHKERR TSGetTimeStepNumber(ts, &step);
#endif

  auto set = [&](auto &fe) {
    fe.ts_u = u;
    fe.ts_u_t = u_t;
    fe.ts_u_tt = u_tt;
    fe.ts_F = F;
    fe.ts_t = t;
    fe.ts_step = step;
    fe.ts_ctx = TSMethod::CTX_TSSETIFUNCTION;
    fe.snes_ctx = SnesMethod::CTX_SNESSETFUNCTION;
    fe.ksp_ctx = KspMethod::CTX_SETFUNCTION;
    fe.data_ctx = PetscData::CtxSetF | PetscData::CtxSetX |
                  PetscData::CtxSetX_T | PetscData::CtxSetX_TT |
                  PetscData::CtxSetTime;
    fe.ts = ts;
  };

  auto unset = [&](auto &fe) {
    fe.ts_ctx = TSMethod::CTX_TSNONE;
    fe.snes_ctx = SnesMethod::CTX_SNESNONE;
    fe.ksp_ctx = KspMethod::CTX_KSPNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };

  // preprocess
  for (auto &bit : ts_ctx->preProcessIFunction) {
    bit->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    unset(*bit);
    ts_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }

  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

  // fe loops
  for (auto &lit : ts_ctx->loopsIFunction) {
    lit.second->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    set(*lit.second);
    CHKERR ts_ctx->mField.loop_finite_elements(ts_ctx->problemName, lit.first,
                                               *(lit.second), nullptr,
                                               ts_ctx->bH, cache_ptr);
    unset(*lit.second);
    ts_ctx->vecAssembleSwitch = boost::move(lit.second->vecAssembleSwitch);
  }

  // post process
  for (auto &bit : ts_ctx->postProcessIFunction) {
    bit->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_postProcess(ts_ctx->problemName,
                                                           *bit);
    unset(*bit);
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

} // namespace MoFEM
