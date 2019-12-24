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
//   #include <petsc/private/snesimpl.h>
// #else
//   #include <petsc-private/snesimpl.h>
// #endif

namespace MoFEM {

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

  snes_ctx->vecAssembleSwitch =
      boost::movelib::make_unique<bool>(true);

  for (auto &bit : snes_ctx->preProcess_Rhs) {
    bit->vecAssembleSwitch = boost::move(snes_ctx->vecAssembleSwitch);
    bit->snes = snes;
    bit->snes_x = x;
    bit->snes_f = f;
    bit->snes_ctx = SnesMethod::CTX_SNESSETFUNCTION;
    CHKERR snes_ctx->mField.problem_basic_method_preProcess(
        snes_ctx->problemName, *bit);
    bit->snes_ctx = SnesMethod::CTX_SNESNONE;
    snes_ctx->vecAssembleSwitch = boost::move(bit->vecAssembleSwitch);
  }

  for (auto &lit : snes_ctx->loops_to_do_Rhs) {
    lit.second->vecAssembleSwitch = boost::move(snes_ctx->vecAssembleSwitch);
    lit.second->snes_ctx = SnesMethod::CTX_SNESSETFUNCTION;
    lit.second->snes = snes;
    lit.second->snes_x = x;
    lit.second->snes_f = f;
    CHKERR snes_ctx->mField.loop_finite_elements(
        snes_ctx->problemName, lit.first, *lit.second, nullptr, snes_ctx->bH);
    lit.second->snes_ctx = SnesMethod::CTX_SNESNONE;
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
    bit->snes = snes;
    bit->snes_x = x;
    bit->snes_f = f;
    bit->snes_ctx = SnesMethod::CTX_SNESSETFUNCTION;
    CHKERR snes_ctx->mField.problem_basic_method_postProcess(
        snes_ctx->problemName, *bit);
    bit->snes_ctx = SnesMethod::CTX_SNESNONE;
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

  CHKERR VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR snes_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
      snes_ctx->problemName, COL, x, INSERT_VALUES, SCATTER_REVERSE);
  for (auto &bit : snes_ctx->preProcess_Mat) {
    bit->matAssembleSwitch = boost::move(snes_ctx->matAssembleSwitch);
    bit->snes = snes;
    bit->snes_x = x;
    bit->snes_A = A;
    bit->snes_B = B;
    bit->snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
    CHKERR snes_ctx->mField.problem_basic_method_preProcess(
        snes_ctx->problemName, *bit);
    bit->snes_ctx = SnesMethod::CTX_SNESNONE;
    snes_ctx->matAssembleSwitch = boost::move(bit->matAssembleSwitch);
  }
  for (auto &lit : snes_ctx->loops_to_do_Mat) {
    lit.second->matAssembleSwitch = boost::move(snes_ctx->matAssembleSwitch);
    lit.second->snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
    lit.second->snes = snes;
    lit.second->snes_x = x;
    lit.second->snes_A = A;
    lit.second->snes_B = B;
    CHKERR snes_ctx->mField.loop_finite_elements(
        snes_ctx->problemName, lit.first, *(lit.second), nullptr, snes_ctx->bH);
    lit.second->snes_ctx = SnesMethod::CTX_SNESNONE;
    snes_ctx->matAssembleSwitch = boost::move(lit.second->matAssembleSwitch);
  }
  for (auto &bit : snes_ctx->postProcess_Mat) {
    bit->matAssembleSwitch = boost::move(snes_ctx->matAssembleSwitch);
    bit->snes = snes;
    bit->snes_x = x;
    bit->snes_A = A;
    bit->snes_B = B;
    bit->snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
    CHKERR snes_ctx->mField.problem_basic_method_postProcess(
        snes_ctx->problemName, *bit);
    bit->snes_ctx = SnesMethod::CTX_SNESNONE;
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

} // namespace MoFEM
