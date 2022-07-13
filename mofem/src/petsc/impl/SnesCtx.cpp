

namespace MoFEM {

MoFEMErrorCode SnesCtx::copyLoops(const SnesCtx &snes_ctx) {
  MoFEMFunctionBeginHot;
  loops_to_do_Mat = snes_ctx.loops_to_do_Mat;
  loops_to_do_Rhs = snes_ctx.loops_to_do_Rhs;
  preProcess_Mat = snes_ctx.preProcess_Mat;
  postProcess_Mat = snes_ctx.postProcess_Mat;
  preProcess_Rhs = snes_ctx.preProcess_Rhs;
  postProcess_Rhs = snes_ctx.postProcess_Rhs;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode SnesCtx::clearLoops() {
  MoFEMFunctionBeginHot;
  loops_to_do_Mat.clear();
  loops_to_do_Rhs.clear();
  preProcess_Mat.clear();
  postProcess_Mat.clear();
  preProcess_Rhs.clear();
  postProcess_Rhs.clear();
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode SnesRhs(SNES snes, Vec x, Vec f, void *ctx) {
  SnesCtx *snes_ctx = (SnesCtx *)ctx;
  // PetscValidHeaderSpecific(snes,SNES_CLASSID,1);
  MoFEMFunctionBegin;
  PetscLogEventBegin(snes_ctx->MOFEM_EVENT_SnesRhs, 0, 0, 0, 0);
  CHKERR VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
  if (snes_ctx->vErify) {
    // Verify finite elements, check for not a number
    CHKERR VecAssemblyBegin(f);
    CHKERR VecAssemblyEnd(f);
    MPI_Comm comm = PetscObjectComm((PetscObject)f);
    PetscSynchronizedPrintf(comm, "SNES Verify x\n");
    const Problem *prb_ptr;
    CHKERR snes_ctx->mField.get_problem(snes_ctx->problemName, &prb_ptr);
    CHKERR snes_ctx->mField.getInterface<Tools>()->checkVectorForNotANumber(
        prb_ptr, COL, x);
  }
  CHKERR snes_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
      snes_ctx->problemName, COL, x, INSERT_VALUES, SCATTER_REVERSE);

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
  CHKERR zero_ghost_vec(f);

  snes_ctx->vecAssembleSwitch = boost::movelib::make_unique<bool>(true);

  auto set = [&](auto &fe) {
    fe.snes = snes;
    fe.snes_x = x;
    fe.snes_f = f;
    fe.snes_ctx = SnesMethod::CTX_SNESSETFUNCTION;
    fe.ksp_ctx = KspMethod::CTX_SETFUNCTION;
    fe.data_ctx = PetscData::CtxSetF | PetscData::CtxSetX;
  };

  auto unset = [&](auto &fe) {
    fe.snes_ctx = SnesMethod::CTX_SNESNONE;
    fe.ksp_ctx = KspMethod::CTX_KSPNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };

  for (auto &bit : snes_ctx->preProcess_Rhs) {
    bit->vecAssembleSwitch = boost::move(snes_ctx->vecAssembleSwitch);
    set(*bit);
    CHKERR snes_ctx->mField.problem_basic_method_preProcess(
        snes_ctx->problemName, *bit);
    unset(*bit);
    snes_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }

  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR snes_ctx->mField.cache_problem_entities(snes_ctx->problemName,
                                                 cache_ptr);

  for (auto &lit : snes_ctx->loops_to_do_Rhs) {
    lit.second->vecAssembleSwitch = boost::move(snes_ctx->vecAssembleSwitch);
    set(*lit.second);
    CHKERR snes_ctx->mField.loop_finite_elements(
        snes_ctx->problemName, lit.first, *lit.second, nullptr, snes_ctx->bH,
        cache_ptr);
    unset(*lit.second);
    if (snes_ctx->vErify) {
      // Verify finite elements, check for not a number
      CHKERR VecAssemblyBegin(f);
      CHKERR VecAssemblyEnd(f);
      MPI_Comm comm = PetscObjectComm((PetscObject)f);
      PetscSynchronizedPrintf(comm, "SNES Verify f FE < %s >\n",
                              lit.first.c_str());
      const Problem *prb_ptr;
      CHKERR snes_ctx->mField.get_problem(snes_ctx->problemName, &prb_ptr);
      CHKERR snes_ctx->mField.getInterface<Tools>()->checkVectorForNotANumber(
          prb_ptr, ROW, f);
    }

    snes_ctx->vecAssembleSwitch = boost::move(lit.second->vecAssembleSwitch);
  }

  for (auto &bit : snes_ctx->postProcess_Rhs) {
    bit->vecAssembleSwitch = boost::move(snes_ctx->vecAssembleSwitch);
    set(*bit);
    CHKERR snes_ctx->mField.problem_basic_method_postProcess(
        snes_ctx->problemName, *bit);
    unset(*bit);
    snes_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }

  if (snes_ctx->vecAssembleSwitch) {
    CHKERR VecGhostUpdateBegin(f, ADD_VALUES, SCATTER_REVERSE);
    CHKERR VecGhostUpdateEnd(f, ADD_VALUES, SCATTER_REVERSE);
    CHKERR VecAssemblyBegin(f);
    CHKERR VecAssemblyEnd(f);
  }
  PetscLogEventEnd(snes_ctx->MOFEM_EVENT_SnesRhs, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

PetscErrorCode SnesMat(SNES snes, Vec x, Mat A, Mat B, void *ctx) {
  SnesCtx *snes_ctx = (SnesCtx *)ctx;
  // PetscValidHeaderSpecific(snes,SNES_CLASSID,1);
  MoFEMFunctionBegin;
  PetscLogEventBegin(snes_ctx->MOFEM_EVENT_SnesMat, 0, 0, 0, 0);
  if (snes_ctx->zeroPreCondMatrixB)
    CHKERR MatZeroEntries(B);

  snes_ctx->matAssembleSwitch = boost::movelib::make_unique<bool>(true);

  auto set = [&](auto &fe) {
    fe.snes = snes;
    fe.snes_x = x;
    fe.snes_A = A;
    fe.snes_B = B;
    fe.snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
    fe.ksp_ctx = KspMethod::CTX_OPERATORS;
    fe.data_ctx = PetscData::CtxSetA | PetscData::CtxSetB | PetscData::CtxSetX;
  };

  auto unset = [&](auto &fe) {
    fe.snes_ctx = SnesMethod::CTX_SNESNONE;
    fe.ksp_ctx = KspMethod::CTX_KSPNONE;
    fe.data_ctx = PetscData::CtxSetNone;
  };


  CHKERR VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR snes_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
      snes_ctx->problemName, COL, x, INSERT_VALUES, SCATTER_REVERSE);
  for (auto &bit : snes_ctx->preProcess_Mat) {
    bit->matAssembleSwitch = boost::move(snes_ctx->matAssembleSwitch);
    set(*bit);
    CHKERR snes_ctx->mField.problem_basic_method_preProcess(
        snes_ctx->problemName, *bit);
    unset(*bit);
    snes_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

  auto cache_ptr = boost::make_shared<CacheTuple>();
  CHKERR snes_ctx->mField.cache_problem_entities(snes_ctx->problemName,
                                                 cache_ptr);

  for (auto &lit : snes_ctx->loops_to_do_Mat) {
    lit.second->matAssembleSwitch = boost::move(snes_ctx->matAssembleSwitch);
    set(*lit.second);
    CHKERR snes_ctx->mField.loop_finite_elements(
        snes_ctx->problemName, lit.first, *(lit.second), nullptr, snes_ctx->bH,
        cache_ptr);
    unset(*lit.second);
    snes_ctx->matAssembleSwitch = boost::move(lit.second->matAssembleSwitch);
  }

  for (auto &bit : snes_ctx->postProcess_Mat) {
    bit->matAssembleSwitch = boost::move(snes_ctx->matAssembleSwitch);
    set(*bit);
    CHKERR snes_ctx->mField.problem_basic_method_postProcess(
        snes_ctx->problemName, *bit);
    unset(*bit);
    snes_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }

  if (*snes_ctx->matAssembleSwitch) {
    CHKERR MatAssemblyBegin(B, snes_ctx->typeOfAssembly);
    CHKERR MatAssemblyEnd(B, snes_ctx->typeOfAssembly);
  }
  PetscLogEventEnd(snes_ctx->MOFEM_EVENT_SnesMat, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SnesMoFEMSetAssemblyType(SNES snes, MatAssemblyType type) {
  SnesCtx *snes_ctx;
  // PetscValidHeaderSpecific(snes,SNES_CLASSID,1);
  MoFEMFunctionBegin;
  CHKERR SNESGetApplicationContext(snes, &snes_ctx);
  snes_ctx->typeOfAssembly = type;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SnesMoFEMSetBehavior(SNES snes, MoFEMTypes bh) {
  SnesCtx *snes_ctx;
  MoFEMFunctionBegin;
  CHKERR SNESGetApplicationContext(snes, &snes_ctx);
  snes_ctx->bH = bh;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MoFEMSNESMonitorFields(SNES snes, PetscInt its, PetscReal fgnorm,
                                      SnesCtx *snes_ctx) {
  MoFEMFunctionBegin;
  auto &m_field =  snes_ctx->mField;
  auto problem_ptr = m_field.get_problem(snes_ctx->problemName);
  auto fields_ptr = m_field.get_fields();
  auto dofs = problem_ptr->numeredRowDofsPtr;

  std::vector<double> lnorms(fields_ptr->size(), 0),
      norms(fields_ptr->size(), 0);

  Vec res;
  CHKERR SNESGetFunction(snes, &res, NULL, NULL);

  const double *r;
  CHKERR VecGetArrayRead(res, &r);
  {
    int f = 0;
    for (auto fi : *fields_ptr) {
      const auto lo_uid = FieldEntity::getLoBitNumberUId(fi->bitNumber);
      const auto hi_uid = FieldEntity::getHiBitNumberUId(fi->bitNumber);
      const auto hi = dofs->get<Unique_mi_tag>().upper_bound(hi_uid);
      for (auto lo = dofs->get<Unique_mi_tag>().lower_bound(lo_uid); lo != hi;
           ++lo) {
        const DofIdx loc_idx = (*lo)->getPetscLocalDofIdx();
        if (loc_idx >= 0 && loc_idx < problem_ptr->nbLocDofsRow) {
          lnorms[f] += PetscRealPart(PetscSqr(r[loc_idx]));
        }
      }
      ++f;
    }
  }
  CHKERR VecRestoreArrayRead(res, &r);

  MPIU_Allreduce(&*lnorms.begin(), &*norms.begin(), lnorms.size(), MPIU_REAL,
                 MPIU_SUM, PetscObjectComm((PetscObject)snes));

  std::stringstream s;
  int tl;
  CHKERR PetscObjectGetTabLevel((PetscObject)snes, &tl);
  for (auto t = 0; t != tl; ++t)
    s << "  ";
  s << its << " Function norm " << boost::format("%14.12e") % (double)fgnorm
    << " [";
  {
    int f = 0;
    for (auto fi : *fields_ptr) {
      if (f > 0)
        s << ", ";
      s << boost::format("%14.12e") % (double)PetscSqrtReal(norms[f]);
      ++f;
    }
    s << "]";
  }

  MOFEM_LOG("SNES_WORLD", Sev::inform) << s.str();

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
