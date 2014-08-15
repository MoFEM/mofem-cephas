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

#include <MoFEM.hpp>

using namespace MoFEM;
#include <ArcLengthTools.hpp>

namespace ObosleteUsersModules {

PetscErrorCode arc_length_mult_shell(Mat A,Vec x,Vec f) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  void *void_ctx;
  ierr = MatShellGetContext(A,&void_ctx); CHKERRQ(ierr);
  ArcLengthMatShell *ctx = (ArcLengthMatShell*)void_ctx;
  ierr = MatMult(ctx->Aij,x,f); CHKERRQ(ierr);
  double lambda;
  ierr = ctx->set_lambda(x,&lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  double db_dot_x;
  ierr = VecDot(ctx->arc_ptr->db,x,&db_dot_x); CHKERRQ(ierr);
  double f_lambda;
  f_lambda = ctx->arc_ptr->diag*lambda + db_dot_x;
  ierr = ctx->set_lambda(f,&f_lambda,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAXPY(f,lambda,ctx->arc_ptr->F_lambda); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode pc_apply_arc_length(PC pc,Vec pc_f,Vec pc_x) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  void *void_ctx;
  ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
  PCShellCtx *PCCtx = (PCShellCtx*)void_ctx;
  void *void_MatCtx;
  MatShellGetContext(PCCtx->ShellAij,&void_MatCtx);
  ArcLengthMatShell *MatCtx = (ArcLengthMatShell*)void_MatCtx;
  ierr = PCApply(PCCtx->pc,pc_f,pc_x); CHKERRQ(ierr);
  ierr = PCApply(PCCtx->pc,PCCtx->arc_ptr->F_lambda,PCCtx->arc_ptr->x_lambda); CHKERRQ(ierr);
  double db_dot_pc_x,db_dot_x_lambda;
  ierr = VecDot(PCCtx->arc_ptr->db,pc_x,&db_dot_pc_x); CHKERRQ(ierr);
  ierr = VecDot(PCCtx->arc_ptr->db,PCCtx->arc_ptr->x_lambda,&db_dot_x_lambda); CHKERRQ(ierr);
  double denominator = PCCtx->arc_ptr->diag+db_dot_x_lambda;
  double res_lambda;
  ierr = MatCtx->set_lambda(pc_f,&res_lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  double ddlambda = (res_lambda - db_dot_pc_x)/denominator;
  if(ddlambda != ddlambda) {
    ostringstream ss;
    ss << "problem with ddlambda: " << res_lambda << " " << ddlambda << " " << db_dot_pc_x << " " << db_dot_x_lambda << " " << PCCtx->arc_ptr->diag;
    //cerr << ss.str() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }
  ierr = VecAXPY(pc_x,ddlambda,PCCtx->arc_ptr->x_lambda); CHKERRQ(ierr);
  ierr = MatCtx->set_lambda(pc_x,&ddlambda,SCATTER_REVERSE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode pc_setup_arc_length(PC pc) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  void *void_ctx;
  ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
  PCShellCtx *ctx = (PCShellCtx*)void_ctx;
  ierr = PCSetFromOptions(ctx->pc); CHKERRQ(ierr);
  ierr = PCGetOperators(pc,&ctx->ShellAij,&ctx->Aij); CHKERRQ(ierr);
  ierr = PCSetOperators(ctx->pc,ctx->ShellAij,ctx->Aij); CHKERRQ(ierr);
  ierr = PCSetUp(ctx->pc); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

}
