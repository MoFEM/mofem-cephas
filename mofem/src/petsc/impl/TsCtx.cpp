

namespace MoFEM {

TsCtx::TsCtx(MoFEM::Interface &m_field, const std::string &problem_name)
    : mField(m_field), moab(m_field.get_moab()), problemName(problem_name),
      bH(MF_EXIST), zeroMatrix(true) {
  PetscLogEventRegister("LoopTsIFunction", 0, &MOFEM_EVENT_TsCtxIFunction);
  PetscLogEventRegister("LoopTsIJacobian", 0, &MOFEM_EVENT_TsCtxIJacobian);
  PetscLogEventRegister("LoopTsRHSFunction", 0, &MOFEM_EVENT_TsCtxRHSFunction);
  PetscLogEventRegister("LoopTsRHSJacobian", 0, &MOFEM_EVENT_TsCtxRHSJacobian);
  PetscLogEventRegister("LoopTsMonitor", 0, &MOFEM_EVENT_TsCtxMonitor);
  PetscLogEventRegister("LoopTsI2Function", 0, &MOFEM_EVENT_TsCtxI2Function);
  PetscLogEventRegister("LoopTsI2Jacobian", 0, &MOFEM_EVENT_TsCtxI2Jacobian);

  if (!LogManager::checkIfChannelExist("TSWORLD")) {
    auto core_log = logging::core::get();

    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "TSWORLD"));
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSync(), "TSSYNC"));
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSelf(), "TSSELF"));

    LogManager::setLog("TSWORLD");
    LogManager::setLog("TSSYNC");
    LogManager::setLog("TSSELF");

    MOFEM_LOG_TAG("TSWORLD", "TS");
    MOFEM_LOG_TAG("TSSYNC", "TS");
    MOFEM_LOG_TAG("TSSELF", "TS");
  }
}

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

  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

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

    CHKERR TSGetSNES(ts, &fe.snes);
    CHKERR SNESGetKSP(fe.snes, &fe.ksp);

    fe.cacheWeakPtr = cache_ptr;
    CHK_THROW_MESSAGE(TSGetTimeStep(ts, &fe.ts_dt), "get time step failed");
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
  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

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

    CHKERR TSGetSNES(ts, &fe.snes);
    CHKERR SNESGetKSP(fe.snes, &fe.ksp);

    fe.cacheWeakPtr = cache_ptr;
    CHK_THROW_MESSAGE(TSGetTimeStep(ts, &fe.ts_dt), "get time step failed");
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

  if (*(ts_ctx->matAssembleSwitch)) {
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

  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

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

    CHKERR TSGetSNES(ts, &fe.snes);
    CHKERR SNESGetKSP(fe.snes, &fe.ksp);

    fe.cacheWeakPtr = cache_ptr;
    CHK_THROW_MESSAGE(TSGetTimeStep(ts, &fe.ts_dt), "get time step failed");
  };

  auto unset = [&](auto &fe) {
    fe.ts_ctx = TSMethod::CTX_TSNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };

  // preprocess
  for (auto &bit : ts_ctx->preProcessMonitor) {
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    unset(*bit);
  }

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
  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

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
    fe.cacheWeakPtr = cache_ptr;
    CHK_THROW_MESSAGE(TSGetTimeStep(ts, &fe.ts_dt), "get time step failed");
  };

  auto unset = [&](auto &fe) {
    fe.ts_ctx = TSMethod::CTX_TSNONE;
    fe.snes_ctx = SnesMethod::CTX_SNESNONE;
    fe.ksp_ctx = KspMethod::CTX_KSPNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };

  for (auto &bit : ts_ctx->preProcessRHSFunction) {
    bit->vecAssembleSwitch = boost::move(ts_ctx->vecAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    unset(*bit);
    ts_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }

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
  for (auto &bit : ts_ctx->postProcessRHSFunction) {
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
  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

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
    fe.cacheWeakPtr = cache_ptr;
    CHK_THROW_MESSAGE(TSGetTimeStep(ts, &fe.ts_dt), "get time step failed");
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

  if (*(ts_ctx->matAssembleSwitch)) {
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
  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

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

    CHKERR TSGetSNES(ts, &fe.snes);
    CHKERR SNESGetKSP(fe.snes, &fe.ksp);
                  
    fe.ts = ts;
    fe.cacheWeakPtr = cache_ptr;
    CHK_THROW_MESSAGE(TSGetTimeStep(ts, &fe.ts_dt), "get time step failed");
  };

  auto unset = [&](auto &fe) {
    fe.ts_ctx = TSMethod::CTX_TSNONE;
    fe.snes_ctx = SnesMethod::CTX_SNESNONE;
    fe.ksp_ctx = KspMethod::CTX_KSPNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };

  // preprocess
  for (auto &bit : ts_ctx->preProcessIJacobian) {
    bit->matAssembleSwitch = boost::move(ts_ctx->matAssembleSwitch);
    set(*bit);
    CHKERR ts_ctx->mField.problem_basic_method_preProcess(ts_ctx->problemName,
                                                          *bit);
    unset(*bit);
    ts_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

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

  if (*(ts_ctx->matAssembleSwitch)) {
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
  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ts_ctx->mField.cache_problem_entities(ts_ctx->problemName, cache_ptr);

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

    CHKERR TSGetSNES(ts, &fe.snes);
    CHKERR SNESGetKSP(fe.snes, &fe.ksp);

    fe.ts = ts;
    fe.cacheWeakPtr = cache_ptr;
    CHK_THROW_MESSAGE(TSGetTimeStep(ts, &fe.ts_dt), "get time step failed");
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

TSAdaptMoFEM::TSAdaptMoFEM()
    : alpha(0.75), gamma(0.5), desiredIt(6), offApat(PETSC_FALSE) {
  CHKERR PetscOptionsGetScalar("", "-ts_mofem_adapt_alpha", &alpha, PETSC_NULL);
  CHKERR PetscOptionsGetScalar("", "-ts_mofem_adapt_gamma", &gamma, PETSC_NULL);
  CHKERR PetscOptionsGetInt("", "-ts_mofem_adapt_desired_it", &desiredIt,
                            PETSC_NULL);
  CHKERR PetscOptionsGetBool("", "-ts_mofem_adapt_off", &offApat, PETSC_NULL);

  MOFEM_LOG("TSWORLD", Sev::inform)
      << "TS adaptivity: alpha = " << alpha << ", gamma = " << gamma
      << ", desiredIt = " << desiredIt << ", offAdapt = " << offApat;
}

PetscErrorCode TSAdaptChooseMoFEM(TSAdapt adapt, TS ts, PetscReal h,
                                  PetscInt *next_sc, PetscReal *next_h,
                                  PetscBool *accept, PetscReal *wlte,
                                  PetscReal *wltea, PetscReal *wlter) {
  MoFEMFunctionBegin;

  auto ts_adapt_mofem = boost::make_shared<TSAdaptMoFEM>();

  *next_sc = 0; /* Reuse the same order scheme */
  *wlte = -1;   /* Weighted local truncation error was not evaluated */
  *wltea = -1;  /* Weighted absolute local truncation error is not used */
  *wlter = -1;  /* Weighted relative local truncation error is not used */

  *accept = PETSC_TRUE;
  *next_h = h; /* Reuse the old step */

  if (!ts_adapt_mofem->offApat) {

    SNES snes;
    CHKERR TSGetSNES(ts, &snes);

    SNESConvergedReason reason;
    CHKERR SNESGetConvergedReason(snes, &reason);

    int it;
    CHKERR SNESGetIterationNumber(snes, &it);

    if (reason < 0) {
      h *= ts_adapt_mofem->alpha;
      *next_h = h;
      MOFEM_LOG_C(
          "TSWORLD", Sev::warning,
          "\tDiverged set step length: it = %d, h = %3.4g set h = %3.4g \n", it,
          h, *next_h);
    } else if (reason > 0) {
      h *= pow((static_cast<double>(ts_adapt_mofem->desiredIt) /
                static_cast<double>(it + 1)),
               static_cast<double>(ts_adapt_mofem->gamma));
      *next_h = PetscClipInterval(h, adapt->dt_min, adapt->dt_max);
      MOFEM_LOG_C(
          "TSWORLD", Sev::inform,
          "\tConverged set step length: it = %d, h = %3.4g set h = %3.4g \n",
          it, h, *next_h);
    }
  }

  MoFEMFunctionReturn(0);
}

PetscErrorCode TSAdaptResetMoFEM(TSAdapt adapt) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

PetscErrorCode TSAdaptDestroyMoFEM(TSAdapt adapt) {
  MoFEMFunctionBegin;
  CHKERR TSAdaptResetMoFEM(adapt);
  MoFEMFunctionReturn(0);
}

PetscErrorCode TSAdaptCreateMoFEM(TSAdapt adapt) {
  PetscFunctionBegin;
  adapt->ops->choose = TSAdaptChooseMoFEM;
  adapt->ops->reset = TSAdaptResetMoFEM;
  adapt->ops->destroy = TSAdaptDestroyMoFEM;
  PetscFunctionReturn(0);
}

} // namespace MoFEM
