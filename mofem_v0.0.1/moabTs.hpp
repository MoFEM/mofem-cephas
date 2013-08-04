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

#ifndef __MOABTS_HPP__
#define __MOABTS_HPP__

#include "moabField.hpp"
#include <petsc.h>
#include <petscmat.h>
#include <petscsnes.h>

#include<petscts.h>

namespace MoFEM {

struct moabTsCtx {

  ErrorCode rval;
  PetscErrorCode ierr;

  moabField &mField;
  Interface &moab;

  string problem_name;

  typedef pair<string,moabField::FEMethod*> loop_pair_type;
  typedef vector<loop_pair_type > loops_to_do_type;
  loops_to_do_type loops_to_do_IJacobian;
  loops_to_do_type loops_to_do_IFunction;
  loops_to_do_type loops_to_do_RHSFunction;
  loops_to_do_type loops_to_do_RHSJacobian;
  loops_to_do_type loops_to_do_Monitor;

  PetscLogEvent USER_EVENT_moabTsRHSFunction;
  PetscLogEvent USER_EVENT_moabTsRHSJacobian;
  PetscLogEvent USER_EVENT_moabTsIFunction;
  PetscLogEvent USER_EVENT_moabTsIJacobian;
  PetscLogEvent USER_EVENT_moabTsMonitor;

  moabTsCtx(moabField &_mField,const string &_problem_name): 
    mField(_mField),moab(_mField.get_moab()),problem_name(_problem_name) {
    PetscLogEventRegister("LoopTsIFunction",0,&USER_EVENT_moabTsIFunction);
    PetscLogEventRegister("LoopTsIJacobian",0,&USER_EVENT_moabTsIJacobian);
    PetscLogEventRegister("LoopTsRHSFunction",0,&USER_EVENT_moabTsRHSFunction);
    PetscLogEventRegister("LoopTsRHSJacobian",0,&USER_EVENT_moabTsRHSJacobian);
    PetscLogEventRegister("LoopTsMonitor",0,&USER_EVENT_moabTsMonitor);
  }

  const moabField& get_mField() const { return mField; }
  const Interface& get_moab() const { return moab; }
  loops_to_do_type& get_loops_to_do_RHSFunction() { return loops_to_do_RHSFunction; }
  loops_to_do_type& get_loops_to_do_RHSJacobian() { return loops_to_do_RHSJacobian; }
  loops_to_do_type& get_loops_to_do_IFunction() { return loops_to_do_IFunction; }
  loops_to_do_type& get_loops_to_do_IJacobian() { return loops_to_do_IJacobian; }
  loops_to_do_type& get_loops_to_do_Monitor() { return loops_to_do_Monitor; }

  friend PetscErrorCode f_TSSetIFunction(TS ts,PetscReal t,Vec u,Vec u_t,Vec F,void *ctx);
  friend PetscErrorCode f_TSSetIJacobian(TS ts,PetscReal t,Vec u,Vec U_t,PetscReal a,Mat *A,Mat *B,MatStructure *flag,void *ctx);
  friend PetscErrorCode f_TSSetRHSFunction(TS ts,PetscReal t,Vec u,Vec F,void *ctx);
  friend PetscErrorCode f_TSSetRHSJacobian(TS ts,PetscReal t,Vec u,Mat *A,Mat *B,MatStructure *flag,void *ctx);
  friend PetscErrorCode f_TSMonitorSet(TS ts,PetscInt step,PetscReal t,Vec u,void *ctx);

};

PetscErrorCode f_TSSetIFunction(TS ts,PetscReal t,Vec u,Vec u_t,Vec F,void *ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  moabTsCtx* ts_ctx = (moabTsCtx*)ctx;
  PetscLogEventBegin(ts_ctx->USER_EVENT_moabTsIFunction,0,0,0,0);
  ierr = VecGhostUpdateBegin(u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(u_t,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(u_t,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = ts_ctx->mField.set_local_VecCreateGhost(ts_ctx->problem_name,Row,u,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  moabTsCtx::loops_to_do_type::iterator lit = ts_ctx->loops_to_do_IFunction.begin();
  for(;lit!=ts_ctx->loops_to_do_IFunction.end();lit++) {
    lit->second->ts_u = u;
    lit->second->ts_u_t = u_t;
    lit->second->ts_F = F;
    lit->second->ts_t = t;
    ierr = lit->second->set_ts_ctx(moabField::TSMethod::ctx_TSSetIFunction);
    ierr = lit->second->set_ts(ts); CHKERRQ(ierr);
    ierr = ts_ctx->mField.loop_finite_elements(ts_ctx->problem_name,lit->first,*(lit->second)); CHKERRQ(ierr);
    ierr = lit->second->set_ts_ctx(moabField::TSMethod::ctx_TSNone);
  }
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  PetscLogEventEnd(ts_ctx->USER_EVENT_moabTsIFunction,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode f_TSSetIJacobian(TS ts,PetscReal t,Vec u,Vec u_t,PetscReal a,Mat *A,Mat *B,MatStructure *flag,void *ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  moabTsCtx* ts_ctx = (moabTsCtx*)ctx;
  PetscLogEventBegin(ts_ctx->USER_EVENT_moabTsIFunction,0,0,0,0);
  moabTsCtx::loops_to_do_type::iterator lit = ts_ctx->loops_to_do_IJacobian.begin();
  ierr = MatZeroEntries(*B); CHKERRQ(ierr);
  for(;lit!=ts_ctx->loops_to_do_IJacobian.end();lit++) {
    lit->second->ts_u = u;
    lit->second->ts_u_t = u_t;
    lit->second->ts_A = A;
    lit->second->ts_B = B;
    lit->second->ts_flag = flag;
    lit->second->ts_t = t;
    lit->second->ts_a = a;
    ierr = lit->second->set_ts_ctx(moabField::TSMethod::ctx_TSSetIJacobian);
    ierr = lit->second->set_ts(ts); CHKERRQ(ierr);
    ierr = ts_ctx->mField.loop_finite_elements(ts_ctx->problem_name,lit->first,*(lit->second)); CHKERRQ(ierr);
    ierr = lit->second->set_ts_ctx(moabField::TSMethod::ctx_TSNone);
  }
  ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  PetscLogEventEnd(ts_ctx->USER_EVENT_moabTsIFunction,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode f_TSSetRHSFunction(TS ts,PetscReal t,Vec u,Vec F,void *ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  moabTsCtx* ts_ctx = (moabTsCtx*)ctx;
  PetscLogEventBegin(ts_ctx->USER_EVENT_moabTsRHSFunction,0,0,0,0);
  ierr = VecGhostUpdateBegin(u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = ts_ctx->mField.set_local_VecCreateGhost(ts_ctx->problem_name,Row,u,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  moabTsCtx::loops_to_do_type::iterator lit = ts_ctx->loops_to_do_RHSFunction.begin();
  for(;lit!=ts_ctx->loops_to_do_RHSFunction.end();lit++) {
    lit->second->ts_u = u;
    lit->second->ts_F = F;
    lit->second->ts_t = t;
    ierr = lit->second->set_ts_ctx(moabField::TSMethod::ctx_TSSetRHSFunction);
    ierr = lit->second->set_ts(ts); CHKERRQ(ierr);
    ierr = ts_ctx->mField.loop_finite_elements(ts_ctx->problem_name,lit->first,*(lit->second)); CHKERRQ(ierr);
    ierr = lit->second->set_ts_ctx(moabField::TSMethod::ctx_TSNone);
  }
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  PetscLogEventEnd(ts_ctx->USER_EVENT_moabTsRHSFunction,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode f_TSSetRHSJacobian(TS ts,PetscReal t,Vec u,Mat *A,Mat *B,MatStructure *flag,void *ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  moabTsCtx* ts_ctx = (moabTsCtx*)ctx;
  PetscLogEventBegin(ts_ctx->USER_EVENT_moabTsRHSJacobian,0,0,0,0);
  moabTsCtx::loops_to_do_type::iterator lit = ts_ctx->loops_to_do_RHSJacobian.begin();
  ierr = MatZeroEntries(*B); CHKERRQ(ierr);
  for(;lit!=ts_ctx->loops_to_do_RHSJacobian.end();lit++) {
    lit->second->ts_u = u;
    lit->second->ts_A = A;
    lit->second->ts_B = B;
    lit->second->ts_flag = flag;
    lit->second->ts_t = t;
    ierr = lit->second->set_ts_ctx(moabField::TSMethod::ctx_TSSetRHSJacobian);
    ierr = lit->second->set_ts(ts); CHKERRQ(ierr);
    ierr = ts_ctx->mField.loop_finite_elements(ts_ctx->problem_name,lit->first,*(lit->second)); CHKERRQ(ierr);
    ierr = lit->second->set_ts_ctx(moabField::TSMethod::ctx_TSNone);
  }
  ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  PetscLogEventEnd(ts_ctx->USER_EVENT_moabTsRHSJacobian,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode f_TSMonitorSet(TS ts,PetscInt step,PetscReal t,Vec u,void *ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  moabTsCtx* ts_ctx = (moabTsCtx*)ctx;
  PetscLogEventBegin(ts_ctx->USER_EVENT_moabTsRHSFunction,0,0,0,0);
  ierr = VecGhostUpdateBegin(u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = ts_ctx->mField.set_local_VecCreateGhost(ts_ctx->problem_name,Row,u,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  moabTsCtx::loops_to_do_type::iterator lit = ts_ctx->loops_to_do_Monitor.begin();
  for(;lit!=ts_ctx->loops_to_do_Monitor.end();lit++) {
    lit->second->ts_u = u;
    lit->second->ts_t = t;
    lit->second->ts_step = step;
    lit->second->ts_F = PETSC_NULL;
    ierr = lit->second->set_ts_ctx(moabField::TSMethod::ctx_TSTSMonitorSet);
    ierr = lit->second->set_ts(ts); CHKERRQ(ierr);
    ierr = ts_ctx->mField.loop_finite_elements(ts_ctx->problem_name,lit->first,*(lit->second)); CHKERRQ(ierr);
    ierr = lit->second->set_ts_ctx(moabField::TSMethod::ctx_TSNone);
  }
  PetscLogEventEnd(ts_ctx->USER_EVENT_moabTsRHSFunction,0,0,0,0);
  PetscFunctionReturn(0);
}


}

#endif // __MOABTS_HPP__
