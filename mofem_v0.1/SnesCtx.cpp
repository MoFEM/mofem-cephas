/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
 */

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

#include "SnesCtx.hpp"

namespace MoFEM {

PetscErrorCode SnesRhs(SNES snes,Vec x,Vec f,void *ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  SnesCtx* snes_ctx = (SnesCtx*)ctx;
  PetscLogEventBegin(snes_ctx->USER_EVENT_SnesRhs,0,0,0,0);
  ierr = VecGhostUpdateBegin(x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = snes_ctx->mField.set_local_VecCreateGhost(snes_ctx->problem_name,COL,x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecZeroEntries(f); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  SnesCtx::basic_method_to_do::iterator bit = snes_ctx->preProcess_Rhs.begin();
  for(;bit!=snes_ctx->preProcess_Rhs.end();bit++) {
    ierr = (*bit)->set_snes(snes); CHKERRQ(ierr);
    (*bit)->snes_x = x;
    (*bit)->snes_f = f;
    ierr = (*bit)->set_snes_ctx(FieldInterface::SnesMethod::CTX_SNESSETFUNCTION);  CHKERRQ(ierr);
    ierr = snes_ctx->mField.problem_basic_method_preProcess(snes_ctx->problem_name,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->set_snes_ctx(FieldInterface::SnesMethod::CTX_SNESNONE);  CHKERRQ(ierr);
  }
  SnesCtx::loops_to_do_type::iterator lit = snes_ctx->loops_to_do_Rhs.begin();
  for(;lit!=snes_ctx->loops_to_do_Rhs.end();lit++) {
    ierr = lit->second->set_snes_ctx(FieldInterface::SnesMethod::CTX_SNESSETFUNCTION);  CHKERRQ(ierr);
    ierr = lit->second->set_snes(snes); CHKERRQ(ierr);
    lit->second->snes_x = x;
    lit->second->snes_f = f;
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\t\tLoop FE for Rhs: %s\n",lit->first.c_str());
    //PetscSynchronizedFlush(PETSC_COMM_WORLD);
    ierr = snes_ctx->mField.loop_finite_elements(snes_ctx->problem_name,lit->first,*(lit->second));  CHKERRQ(ierr);
    ierr = lit->second->set_snes_ctx(FieldInterface::SnesMethod::CTX_SNESNONE);  CHKERRQ(ierr);
  }
  bit = snes_ctx->postProcess_Rhs.begin();
  for(;bit!=snes_ctx->postProcess_Rhs.end();bit++) {
    ierr = (*bit)->set_snes(snes); CHKERRQ(ierr);
    (*bit)->snes_x = x;
    (*bit)->snes_f = f;
    ierr = (*bit)->set_snes_ctx(FieldInterface::SnesMethod::CTX_SNESSETFUNCTION);  CHKERRQ(ierr);
    ierr = snes_ctx->mField.problem_basic_method_postProcess(snes_ctx->problem_name,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->set_snes_ctx(FieldInterface::SnesMethod::CTX_SNESNONE);  CHKERRQ(ierr);
  }
  ierr = VecGhostUpdateBegin(f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(f); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(f); CHKERRQ(ierr);
  PetscLogEventEnd(snes_ctx->USER_EVENT_SnesRhs,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode SnesMat(SNES snes,Vec x,Mat *A,Mat *B,MatStructure *flag,void *ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  SnesCtx* snes_ctx = (SnesCtx*)ctx;
  PetscLogEventBegin(snes_ctx->USER_EVENT_SnesMat,0,0,0,0);
  ierr = MatZeroEntries(*B); CHKERRQ(ierr);
  SnesCtx::basic_method_to_do::iterator bit = snes_ctx->preProcess_Mat.begin();
  for(;bit!=snes_ctx->preProcess_Mat.end();bit++) {
    ierr = (*bit)->set_snes(snes); CHKERRQ(ierr);
    (*bit)->snes_x = x;
    (*bit)->snes_A = A;
    (*bit)->snes_B = B;
    (*bit)->snes_flag = flag;
    ierr = (*bit)->set_snes_ctx(FieldInterface::SnesMethod::CTX_SNESSETJACOBIAN); CHKERRQ(ierr);
    ierr = snes_ctx->mField.problem_basic_method_preProcess(snes_ctx->problem_name,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->set_snes_ctx(FieldInterface::SnesMethod::CTX_SNESNONE);  CHKERRQ(ierr);
  }
  SnesCtx::loops_to_do_type::iterator lit = snes_ctx->loops_to_do_Mat.begin();
  for(;lit!=snes_ctx->loops_to_do_Mat.end();lit++) {
    ierr = lit->second->set_snes_ctx(FieldInterface::SnesMethod::CTX_SNESSETJACOBIAN); CHKERRQ(ierr);
    ierr = lit->second->set_snes(snes); CHKERRQ(ierr);
    lit->second->snes_x = x;
    lit->second->snes_A = A;
    lit->second->snes_B = B;
    lit->second->snes_flag = flag;
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\t\tLoop FE for Mat: %s\n",lit->first.c_str());
    //PetscSynchronizedFlush(PETSC_COMM_WORLD);
    ierr = snes_ctx->mField.loop_finite_elements(snes_ctx->problem_name,lit->first,*(lit->second));  CHKERRQ(ierr);
    ierr = lit->second->set_snes_ctx(FieldInterface::SnesMethod::CTX_SNESNONE);
  }
  bit = snes_ctx->postProcess_Mat.begin();
  for(;bit!=snes_ctx->postProcess_Mat.end();bit++) {
    ierr = (*bit)->set_snes(snes); CHKERRQ(ierr);
    (*bit)->snes_x = x;
    (*bit)->snes_A = A;
    (*bit)->snes_B = B;
    (*bit)->snes_flag = flag;
    ierr = (*bit)->set_snes_ctx(FieldInterface::SnesMethod::CTX_SNESSETJACOBIAN); CHKERRQ(ierr);
    ierr = snes_ctx->mField.problem_basic_method_postProcess(snes_ctx->problem_name,*(*(bit)));  CHKERRQ(ierr);
    ierr = (*bit)->set_snes_ctx(FieldInterface::SnesMethod::CTX_SNESNONE);  CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  PetscLogEventEnd(snes_ctx->USER_EVENT_SnesMat,0,0,0,0);
  PetscFunctionReturn(0);
}

}

