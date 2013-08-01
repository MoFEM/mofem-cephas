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

#ifndef __MOABSNES_HPP__
#define __MOABSNES_HPP__

#include "moabField.hpp"
#include <petsc.h>
#include <petscmat.h>
#include <petscsnes.h>

namespace MoFEM {

struct moabSnesCtx {

  ErrorCode rval;
  PetscErrorCode ierr;

  moabField &mField;
  Interface &moab;

  string problem_name;

  typedef pair<string,moabField::FEMethod*> loop_pair_type;
  typedef vector<loop_pair_type > loops_to_do_type;
  loops_to_do_type loops_to_do_Mat;
  loops_to_do_type loops_to_do_Rhs;

  PetscLogEvent USER_EVENT_SnesRhs;
  PetscLogEvent USER_EVENT_SnesMat;

  moabSnesCtx(moabField &_mField,const string &_problem_name): 
    mField(_mField),moab(_mField.get_moab()),problem_name(_problem_name) {
    PetscLogEventRegister("LoopSnesRhs",0,&USER_EVENT_SnesRhs);
    PetscLogEventRegister("LoopSnesMat",0,&USER_EVENT_SnesMat);
  }

  const moabField& get_mField() const { return mField; }
  const Interface& get_moab() const { return moab; }
  loops_to_do_type& get_loops_to_do_Mat() { return loops_to_do_Mat; }
  loops_to_do_type& get_loops_to_do_Rhs() { return loops_to_do_Rhs; }

  friend PetscErrorCode SnesFunc(SNES snes,Vec x,Vec f,moabSnesCtx *);
  friend PetscErrorCode SnesMat(SNES snes,Vec x,Mat *A,Mat *B,MatStructure *flag,void *ctx);

};

PetscErrorCode SnesRhs(SNES snes,Vec x,Vec f,void *ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  moabSnesCtx* snes_ctx = (moabSnesCtx*)ctx;
  PetscLogEventBegin(snes_ctx->USER_EVENT_SnesRhs,0,0,0,0);
  ierr = VecGhostUpdateBegin(x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = snes_ctx->mField.set_local_VecCreateGhost(snes_ctx->problem_name,Row,x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = snes_ctx->mField.set_global_VecCreateGhost(snes_ctx->problem_name,Row,x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecZeroEntries(f); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  moabSnesCtx::loops_to_do_type::iterator lit = snes_ctx->loops_to_do_Rhs.begin();
  for(;lit!=snes_ctx->loops_to_do_Rhs.end();lit++) {
    ierr = lit->second->set_snes_ctx(moabField::SnesMethod::ctx_SNESSetFunction);
    ierr = lit->second->set_snes(snes); CHKERRQ(ierr);
    lit->second->snes_x = x;
    lit->second->snes_f = f;
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\t\tLoop FE for Rhs: %s\n",lit->first.c_str());
    //PetscSynchronizedFlush(PETSC_COMM_WORLD);
    ierr = snes_ctx->mField.loop_finite_elements(snes_ctx->problem_name,lit->first,*(lit->second));  CHKERRQ(ierr);
    ierr = lit->second->set_snes_ctx(moabField::SnesMethod::ctx_SNESNone);
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
  moabSnesCtx* snes_ctx = (moabSnesCtx*)ctx;
  PetscLogEventBegin(snes_ctx->USER_EVENT_SnesMat,0,0,0,0);
  ierr = MatZeroEntries(*B); CHKERRQ(ierr);
  moabSnesCtx::loops_to_do_type::iterator lit = snes_ctx->loops_to_do_Mat.begin();
  for(;lit!=snes_ctx->loops_to_do_Mat.end();lit++) {
    ierr = lit->second->set_snes_ctx(moabField::SnesMethod::ctx_SNESSetJacobian);
    ierr = lit->second->set_snes(snes); CHKERRQ(ierr);
    lit->second->snes_x = x;
    lit->second->snes_A = A;
    lit->second->snes_B = B;
    lit->second->snes_flag = flag;
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\t\tLoop FE for Mat: %s\n",lit->first.c_str());
    //PetscSynchronizedFlush(PETSC_COMM_WORLD);
    ierr = snes_ctx->mField.loop_finite_elements(snes_ctx->problem_name,lit->first,*(lit->second));  CHKERRQ(ierr);
    ierr = lit->second->set_snes_ctx(moabField::SnesMethod::ctx_SNESNone);
  }
  ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  PetscLogEventEnd(snes_ctx->USER_EVENT_SnesMat,0,0,0,0);
  PetscFunctionReturn(0);
}

}

#endif // __MOABSNES_HPP__
