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

#ifndef __MOABCONSTRAINSBYMARKAINSWORTH_HPP__
#define __MOABCONSTRAINSBYMARKAINSWORTH_HPP__

struct matPROJ_ctx {
  moabField& mField;
  Mat C,CT,CCT;
  KSP ksp;
  Vec _x_;
  VecScatter scatter;
  Vec Cx,CCTm1_Cx,CT_CCTm1_Cx;
  string x_problem,y_problem;
  bool init;
  bool debug;
  matPROJ_ctx(moabField& _mField,Mat _C,Mat _CT,Mat _CCT,string _x_problem,string _y_problem): mField(_mField),C(_C),CT(_CT),CCT(_CCT),
    x_problem(_x_problem),y_problem(_y_problem),init(true),debug(true) {}
  PetscErrorCode Destroy() {
    PetscFunctionBegin;
    if(init) PetscFunctionReturn(0);
    PetscErrorCode ierr;
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);
    ierr = VecDestroy(&_x_); CHKERRQ(ierr);
    ierr = VecDestroy(&Cx); CHKERRQ(ierr);
    ierr = VecDestroy(&CCTm1_Cx); CHKERRQ(ierr);
    ierr = VecDestroy(&CT_CCTm1_Cx); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  friend PetscErrorCode matQ_mult_shell(Mat Q,Vec x,Vec f);
  friend PetscErrorCode matP_mult_shell(Mat P,Vec x,Vec f);
};

PetscErrorCode matQ_mult_shell(Mat Q,Vec x,Vec f) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  void *void_ctx;
  ierr = MatShellGetContext(Q,&void_ctx); CHKERRQ(ierr);
  matPROJ_ctx *ctx = (matPROJ_ctx*)void_ctx;
  if(ctx->init) {
    ctx->init = false;
    ierr = KSPCreate(PETSC_COMM_WORLD,&(ctx->ksp)); CHKERRQ(ierr);
    ierr = KSPSetOperators(ctx->ksp,ctx->CCT,ctx->CCT,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ctx->ksp); CHKERRQ(ierr);
    ierr = KSPSetUp(ctx->ksp); CHKERRQ(ierr);
    ierr = MatGetVecs(ctx->C,&ctx->_x_,PETSC_NULL); CHKERRQ(ierr);
    ierr = MatGetVecs(ctx->C,PETSC_NULL,&ctx->Cx); CHKERRQ(ierr);
    ierr = MatGetVecs(ctx->CCT,PETSC_NULL,&ctx->CCTm1_Cx); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx->_x_,&ctx->CT_CCTm1_Cx); CHKERRQ(ierr);
    ierr = ctx->mField.VecScatterCreate(x,ctx->x_problem,Row,ctx->_x_,ctx->y_problem,Col,&ctx->scatter); CHKERRQ(ierr);
  }
  ierr = VecCopy(x,f); CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //ierr = VecView(ctx->_x_,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  if(ctx->debug) {
    ierr = VecScatterBegin(ctx->scatter,ctx->_x_,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatter,ctx->_x_,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    PetscBool  flg;
    ierr = VecEqual(x,f,&flg); CHKERRQ(ierr);
    if(flg ==  PETSC_FALSE) SETERRQ(PETSC_COMM_SELF,1,"scatter is not working");
  }
  ierr = MatMult(ctx->C,ctx->_x_,ctx->Cx);  CHKERRQ(ierr);
  ierr = KSPSolve(ctx->ksp,ctx->Cx,ctx->CCTm1_Cx); CHKERRQ(ierr);
  ierr = MatMult(ctx->CT,ctx->CCTm1_Cx,ctx->CT_CCTm1_Cx);  CHKERRQ(ierr);
  ierr = VecScale(ctx->CT_CCTm1_Cx,-1); CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx->scatter,ctx->CT_CCTm1_Cx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,ctx->CT_CCTm1_Cx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode matP_mult_shell(Mat P,Vec x,Vec f) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  void *void_ctx;
  ierr = MatShellGetContext(P,&void_ctx); CHKERRQ(ierr);
  matPROJ_ctx *ctx = (matPROJ_ctx*)void_ctx;
  if(ctx->init) {
    ctx->init = false;
    ierr = KSPCreate(PETSC_COMM_WORLD,&(ctx->ksp)); CHKERRQ(ierr);
    ierr = KSPSetOperators(ctx->ksp,ctx->CCT,ctx->CCT,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ctx->ksp); CHKERRQ(ierr);
    ierr = KSPSetUp(ctx->ksp); CHKERRQ(ierr);
    ierr = MatGetVecs(ctx->C,&ctx->_x_,PETSC_NULL); CHKERRQ(ierr);
    ierr = MatGetVecs(ctx->C,PETSC_NULL,&ctx->Cx); CHKERRQ(ierr);
    ierr = MatGetVecs(ctx->CCT,PETSC_NULL,&ctx->CCTm1_Cx); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx->_x_,&ctx->CT_CCTm1_Cx); CHKERRQ(ierr);
    ierr = ctx->mField.VecScatterCreate(x,ctx->x_problem,Row,ctx->_x_,ctx->y_problem,Col,&ctx->scatter); CHKERRQ(ierr);
  }
  ierr = VecScatterBegin(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatMult(ctx->C,ctx->_x_,ctx->Cx);  CHKERRQ(ierr);
  ierr = KSPSolve(ctx->ksp,ctx->Cx,ctx->CCTm1_Cx); CHKERRQ(ierr);
  ierr = MatMult(ctx->CT,ctx->CCTm1_Cx,ctx->CT_CCTm1_Cx);  CHKERRQ(ierr);
  ierr = VecZeroEntries(f); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx->scatter,ctx->CT_CCTm1_Cx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,ctx->CT_CCTm1_Cx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#endif //__NONLINEAR_ELASTICITY_HPP__
