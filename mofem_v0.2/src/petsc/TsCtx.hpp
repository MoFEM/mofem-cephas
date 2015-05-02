/** \file TsCtx.hpp 
 * \brief Context for PETSc Time Stepping 
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

namespace MoFEM {

struct TsCtx {

  ErrorCode rval;
  PetscErrorCode ierr;

  FieldInterface &mField;
  Interface &moab;

  string problemName;

  typedef pair<string,FEMethod*> loop_pair_type;
  typedef vector<loop_pair_type > loops_to_do_type;
  loops_to_do_type loops_to_do_IJacobian;
  loops_to_do_type loops_to_do_IFunction;
  loops_to_do_type loops_to_do_Monitor;

  typedef vector<BasicMethod*> basic_method_to_do;
  basic_method_to_do preProcess_IJacobian;
  basic_method_to_do postProcess_IJacobian;
  basic_method_to_do preProcess_IFunction;
  basic_method_to_do postProcess_IFunction;
  basic_method_to_do preProcess_Monitor;
  basic_method_to_do postProcess_Monitor;

  PetscLogEvent USER_EVENT_TsCtxRHSFunction;
  PetscLogEvent USER_EVENT_TsCtxRHSJacobian;
  PetscLogEvent USER_EVENT_TsCtxIFunction;
  PetscLogEvent USER_EVENT_TsCtxIJacobian;
  PetscLogEvent USER_EVENT_TsCtxMonitor;

  bool zero_matrix;
  TsCtx(FieldInterface &_mField,const string &_problem_name): 
    mField(_mField),moab(_mField.get_moab()),problemName(_problem_name),zero_matrix(true) {
    PetscLogEventRegister("LoopTsIFunction",0,&USER_EVENT_TsCtxIFunction);
    PetscLogEventRegister("LoopTsIJacobian",0,&USER_EVENT_TsCtxIJacobian);
    PetscLogEventRegister("LoopTsRHSFunction",0,&USER_EVENT_TsCtxRHSFunction);
    PetscLogEventRegister("LoopTsRHSJacobian",0,&USER_EVENT_TsCtxRHSJacobian);
    PetscLogEventRegister("LoopTsMonitor",0,&USER_EVENT_TsCtxMonitor);
  }

  const FieldInterface& get_mField() const { return mField; }
  const Interface& get_moab() const { return moab; }
  loops_to_do_type& get_loops_to_do_IFunction() { return loops_to_do_IFunction; }
  loops_to_do_type& get_loops_to_do_IJacobian() { return loops_to_do_IJacobian; }
  loops_to_do_type& get_loops_to_do_Monitor() { return loops_to_do_Monitor; }

  basic_method_to_do& get_preProcess_to_do_IFunction() { return preProcess_IFunction; }
  basic_method_to_do& get_postProcess_to_do_IFunction() { return postProcess_IFunction; }
  basic_method_to_do& get_preProcess_to_do_IJacobian() { return preProcess_IJacobian; }
  basic_method_to_do& get_postProcess_to_do_IJacobian() { return postProcess_IJacobian; }
  basic_method_to_do& get_preProcess_to_do_Monitor() { return preProcess_Monitor; }
  basic_method_to_do& get_postProcess_to_do_Monitor() { return postProcess_Monitor; }

  friend PetscErrorCode f_TSSetIFunction(TS ts,PetscReal t,Vec u,Vec u_t,Vec F,void *ctx);
  friend PetscErrorCode f_TSSetIJacobian(TS ts,PetscReal t,Vec u,Vec U_t,PetscReal a,Mat A,Mat B,void *ctx);
  friend PetscErrorCode f_TSMonitorSet(TS ts,PetscInt step,PetscReal t,Vec u,void *ctx);

};

PetscErrorCode f_TSSetIFunction(TS ts,PetscReal t,Vec u,Vec u_t,Vec F,void *ctx);
PetscErrorCode f_TSSetIJacobian(TS ts,PetscReal t,Vec u,Vec u_t,PetscReal a,Mat A,Mat B,void *ctx);
PetscErrorCode f_TSMonitorSet(TS ts,PetscInt step,PetscReal t,Vec u,void *ctx);


}

#endif // __MOABTS_HPP__
