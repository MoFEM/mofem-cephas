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
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <VecManager.hpp>
#include <AuxPETSc.hpp>
#include <SnesCtx.hpp>

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
  CHKERR snes_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
      snes_ctx->problemName, COL, x, INSERT_VALUES, SCATTER_REVERSE);
  CHKERR VecZeroEntries(f);
  CHKERR VecGhostUpdateBegin(f, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(f, INSERT_VALUES, SCATTER_FORWARD);
  SnesCtx::BasicMethodsSequence::iterator bit =
      snes_ctx->preProcess_Rhs.begin();
  for (; bit != snes_ctx->preProcess_Rhs.end(); bit++) {
    CHKERR (*bit)->setSnes(snes);
    (*bit)->snes_x = x;
    (*bit)->snes_f = f;
    CHKERR (*bit)->setSnesCtx(SnesMethod::CTX_SNESSETFUNCTION);
    CHKERR snes_ctx->mField.problem_basic_method_preProcess(
        snes_ctx->problemName, *(*(bit)));
    CHKERR (*bit)->setSnesCtx(SnesMethod::CTX_SNESNONE);
  }
  SnesCtx::FEMethodsSequence::iterator lit = snes_ctx->loops_to_do_Rhs.begin();
  for (; lit != snes_ctx->loops_to_do_Rhs.end(); lit++) {
    CHKERR lit->second->setSnesCtx(SnesMethod::CTX_SNESSETFUNCTION);
    CHKERR lit->second->setSnes(snes);
    lit->second->snes_x = x;
    lit->second->snes_f = f;
    // PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\t\tLoop FE for Rhs:
    // %s\n",lit->first.c_str());  PetscSynchronizedFlush(PETSC_COMM_WORLD);
    CHKERR snes_ctx->mField.loop_finite_elements(
        snes_ctx->problemName, lit->first, *(lit->second), snes_ctx->bH);
    CHKERR lit->second->setSnesCtx(SnesMethod::CTX_SNESNONE);
  }
  bit = snes_ctx->postProcess_Rhs.begin();
  for (; bit != snes_ctx->postProcess_Rhs.end(); bit++) {
    CHKERR (*bit)->setSnes(snes);
    (*bit)->snes_x = x;
    (*bit)->snes_f = f;
    CHKERR (*bit)->setSnesCtx(SnesMethod::CTX_SNESSETFUNCTION);
    CHKERR snes_ctx->mField.problem_basic_method_postProcess(
        snes_ctx->problemName, *(*(bit)));
    CHKERR (*bit)->setSnesCtx(SnesMethod::CTX_SNESNONE);
  }
  CHKERR VecGhostUpdateBegin(f, ADD_VALUES, SCATTER_REVERSE);
  CHKERR VecGhostUpdateEnd(f, ADD_VALUES, SCATTER_REVERSE);
  CHKERR VecAssemblyBegin(f);
  CHKERR VecAssemblyEnd(f);
  PetscLogEventEnd(snes_ctx->MOFEM_EVENT_SnesRhs, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}
PetscErrorCode SnesMat(SNES snes, Vec x, Mat A, Mat B, void *ctx) {
  SnesCtx *snes_ctx = (SnesCtx *)ctx;
  // PetscValidHeaderSpecific(snes,SNES_CLASSID,1);
  MoFEMFunctionBegin;
  PetscLogEventBegin(snes_ctx->MOFEM_EVENT_SnesMat, 0, 0, 0, 0);
  if (snes_ctx->zeroPreCondMatrixB) {
    CHKERR MatZeroEntries(B);
  }
  CHKERR VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR snes_ctx->mField.getInterface<VecManager>()->setLocalGhostVector(
      snes_ctx->problemName, COL, x, INSERT_VALUES, SCATTER_REVERSE);
  SnesCtx::BasicMethodsSequence::iterator bit =
      snes_ctx->preProcess_Mat.begin();
  for (; bit != snes_ctx->preProcess_Mat.end(); bit++) {
    CHKERR (*bit)->setSnes(snes);
    (*bit)->snes_x = x;
    (*bit)->snes_A = A;
    (*bit)->snes_B = B;
    CHKERR (*bit)->setSnesCtx(SnesMethod::CTX_SNESSETJACOBIAN);
    CHKERR snes_ctx->mField.problem_basic_method_preProcess(
        snes_ctx->problemName, *(*(bit)));
    CHKERR (*bit)->setSnesCtx(SnesMethod::CTX_SNESNONE);
  }
  SnesCtx::FEMethodsSequence::iterator lit = snes_ctx->loops_to_do_Mat.begin();
  for (; lit != snes_ctx->loops_to_do_Mat.end(); lit++) {
    CHKERR lit->second->setSnesCtx(SnesMethod::CTX_SNESSETJACOBIAN);
    CHKERR lit->second->setSnes(snes);
    lit->second->snes_x = x;
    lit->second->snes_A = A;
    lit->second->snes_B = B;
    // PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\t\tLoop FE for Mat:
    // %s\n",lit->first.c_str());  PetscSynchronizedFlush(PETSC_COMM_WORLD);
    CHKERR snes_ctx->mField.loop_finite_elements(
        snes_ctx->problemName, lit->first, *(lit->second), snes_ctx->bH);
    CHKERR lit->second->setSnesCtx(SnesMethod::CTX_SNESNONE);
  }
  bit = snes_ctx->postProcess_Mat.begin();
  for (; bit != snes_ctx->postProcess_Mat.end(); bit++) {
    CHKERR (*bit)->setSnes(snes);
    (*bit)->snes_x = x;
    (*bit)->snes_A = A;
    (*bit)->snes_B = B;
    CHKERR (*bit)->setSnesCtx(SnesMethod::CTX_SNESSETJACOBIAN);
    CHKERR snes_ctx->mField.problem_basic_method_postProcess(
        snes_ctx->problemName, *(*(bit)));
    CHKERR (*bit)->setSnesCtx(SnesMethod::CTX_SNESNONE);
  }
  CHKERR MatAssemblyBegin(B, snes_ctx->typeOfAssembly);
  CHKERR MatAssemblyEnd(B, snes_ctx->typeOfAssembly);
  PetscLogEventEnd(snes_ctx->MOFEM_EVENT_SnesMat, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SNESMoFEMSetAssmblyType(SNES snes, MatAssemblyType type) {
  SnesCtx *snes_ctx;
  // PetscValidHeaderSpecific(snes,SNES_CLASSID,1);
  MoFEMFunctionBeginHot;
  ierr = SNESGetApplicationContext(snes, &snes_ctx);
  CHKERRG(ierr);
  snes_ctx->typeOfAssembly = type;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode SNESMoFEMSetBehavior(SNES snes, MoFEMTypes bh) {
  SnesCtx *snes_ctx;
  // PetscValidHeaderSpecific(snes,SNES_CLASSID,1);
  MoFEMFunctionBeginHot;
  ierr = SNESGetApplicationContext(snes, &snes_ctx);
  CHKERRG(ierr);
  snes_ctx->bH = bh;
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM
