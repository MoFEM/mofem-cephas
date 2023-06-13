

// #if PETSC_VERSION_GE(3,6,0)
//   #include <petsc/private/kspimpl.h>
// #else
//   #include <petsc-private/kspimpl.h>
// #endif

namespace MoFEM {

MoFEMErrorCode KspCtx::clearLoops() {
  MoFEMFunctionBegin;
  loops_to_do_Mat.clear();
  loops_to_do_Rhs.clear();
  preProcess_Mat.clear();
  postProcess_Mat.clear();
  preProcess_Rhs.clear();
  MoFEMFunctionReturn(0);
}

PetscErrorCode KspRhs(KSP ksp, Vec f, void *ctx) {
  // PetscValidHeaderSpecific(ksp,KSP_CLASSID,1);
  MoFEMFunctionBegin;
  KspCtx *ksp_ctx = static_cast<KspCtx *>(ctx);
  PetscLogEventBegin(ksp_ctx->MOFEM_EVENT_KspRhs, 0, 0, 0, 0);

  ksp_ctx->vecAssembleSwitch = boost::movelib::make_unique<bool>(true);

  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ksp_ctx->mField.cache_problem_entities(ksp_ctx->problemName,
                                                cache_ptr);

  auto set = [&](auto &fe) {
    fe.ksp = ksp;
    fe.ksp_ctx = KspMethod::CTX_SETFUNCTION;
    fe.data_ctx = PetscData::CtxSetF;
    fe.ksp_f = f;
    fe.cacheWeakPtr = cache_ptr;
  };

  auto unset = [&](auto &fe) {
    fe.ksp_ctx = KspMethod::CTX_KSPNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };

  // pre-process
  for (auto &bit : ksp_ctx->preProcess_Rhs) {
    bit->vecAssembleSwitch = boost::move(ksp_ctx->vecAssembleSwitch);
    set(*bit);
    CHKERR ksp_ctx->mField.problem_basic_method_preProcess(ksp_ctx->problemName,
                                                           *bit);
    unset(*bit);
    ksp_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }

  // operators
  for (auto &lit : ksp_ctx->loops_to_do_Rhs) {
    lit.second->vecAssembleSwitch = boost::move(ksp_ctx->vecAssembleSwitch);
    set(*lit.second);
    CHKERR ksp_ctx->mField.loop_finite_elements(ksp_ctx->problemName, lit.first,
                                                *(lit.second), nullptr,
                                                ksp_ctx->bH, cache_ptr);
    unset(*lit.second);
    ksp_ctx->vecAssembleSwitch = boost::move(lit.second->vecAssembleSwitch);
  }

  // post-process
  for (auto &bit : ksp_ctx->postProcess_Rhs) {
    bit->vecAssembleSwitch = boost::move(ksp_ctx->vecAssembleSwitch);
    set(*bit);
    CHKERR ksp_ctx->mField.problem_basic_method_postProcess(
        ksp_ctx->problemName, *bit);
    unset(*bit);
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
  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR ksp_ctx->mField.cache_problem_entities(ksp_ctx->problemName,
                                                cache_ptr);

  auto set = [&](auto &fe) {
    fe.ksp = ksp;
    fe.ksp_A = A;
    fe.ksp_B = B;
    fe.ksp_ctx = KspMethod::CTX_OPERATORS;
    fe.data_ctx = PetscData::CtxSetA | PetscData::CtxSetB;
    fe.cacheWeakPtr = cache_ptr;
  };

  auto unset = [&](auto &fe) {
    fe.ksp_ctx = KspMethod::CTX_KSPNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };

  auto ent_data_cache = boost::make_shared<std::vector<EntityCacheDofs>>();
  auto ent_row_cache =
      boost::make_shared<std::vector<EntityCacheNumeredDofs>>();
  auto ent_col_cache =
      boost::make_shared<std::vector<EntityCacheNumeredDofs>>();

  // pre-procsess
  for (auto &bit : ksp_ctx->preProcess_Mat) {
    bit->matAssembleSwitch = boost::move(ksp_ctx->matAssembleSwitch);
    set(*bit);
    CHKERR ksp_ctx->mField.problem_basic_method_preProcess(ksp_ctx->problemName,
                                                           *bit);
    unset(*bit);
    ksp_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

  // operators
  for (auto &lit : ksp_ctx->loops_to_do_Mat) {
    lit.second->matAssembleSwitch = boost::move(ksp_ctx->matAssembleSwitch);
    set(*lit.second);
    CHKERR ksp_ctx->mField.loop_finite_elements(ksp_ctx->problemName, lit.first,
                                                *(lit.second), nullptr,
                                                ksp_ctx->bH, cache_ptr);
    unset(*lit.second);
    ksp_ctx->matAssembleSwitch = boost::move(lit.second->matAssembleSwitch);
  }

  // post-process
  for (auto &bit : ksp_ctx->postProcess_Mat) {
    bit->matAssembleSwitch = boost::move(ksp_ctx->matAssembleSwitch);
    set(*bit);
    CHKERR ksp_ctx->mField.problem_basic_method_postProcess(
        ksp_ctx->problemName, *bit);
    unset(*bit);
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
