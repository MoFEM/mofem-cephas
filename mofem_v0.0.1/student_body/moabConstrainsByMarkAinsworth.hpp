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
  Mat C,K;
  KSP ksp;
  Vec _x_;
  VecScatter scatter;
  Mat CT,CCT,CTC;
  Vec Cx,CCTm1_Cx,CT_CCTm1_Cx,CTCx;
  string x_problem,y_problem;
  bool init;
  bool debug;
  PetscLogEvent USER_EVENT_projInit;
  PetscLogEvent USER_EVENT_projQ;
  PetscLogEvent USER_EVENT_projP;
  PetscLogEvent USER_EVENT_projCTC_QTKQ;
  matPROJ_ctx(moabField& _mField,Mat _K,Mat _C,string _x_problem,string _y_problem): 
    mField(_mField),K(_K),C(_C),x_problem(_x_problem),y_problem(_y_problem),init(true),debug(true) {
    PetscLogEventRegister("ProjectionInit",0,&USER_EVENT_projInit);
    PetscLogEventRegister("ProjectionQ",0,&USER_EVENT_projQ);
    PetscLogEventRegister("ProjectionP",0,&USER_EVENT_projP);
    PetscLogEventRegister("ProjectionCTC_QTKQ",0,&USER_EVENT_projCTC_QTKQ);
  }
  PetscErrorCode Init(Vec x) {
    PetscFunctionBegin;
    if(init) {
      init = false;
      PetscErrorCode ierr;
      PetscLogEventBegin(USER_EVENT_projInit,0,0,0,0);
      ierr = MatTranspose(C,MAT_INITIAL_MATRIX,&CT); CHKERRQ(ierr);
      ierr = MatTransposeMatMult(CT,CT,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CCT); CHKERRQ(ierr);
      ierr = MatTransposeMatMult(C,C,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CTC); CHKERRQ(ierr);
      if(debug) {
	//MatView(CCT,PETSC_VIEWER_DRAW_WORLD);
	int m,n;
	MatGetSize(CCT,&m,&n);
	PetscPrintf(PETSC_COMM_WORLD,"CCT_ALL size (%d,%d)\n",m,n);
	//std::string wait;
	//std::cin >> wait;
      }
      ierr = KSPCreate(PETSC_COMM_WORLD,&(ksp)); CHKERRQ(ierr);
      ierr = KSPSetOperators(ksp,CCT,CCT,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
      ierr = KSPSetUp(ksp); CHKERRQ(ierr);
      ierr = MatGetVecs(C,&_x_,PETSC_NULL); CHKERRQ(ierr);
      ierr = MatGetVecs(C,PETSC_NULL,&Cx); CHKERRQ(ierr);
      ierr = MatGetVecs(CCT,PETSC_NULL,&CCTm1_Cx); CHKERRQ(ierr);
      ierr = MatGetVecs(CTC,PETSC_NULL,&CTCx); CHKERRQ(ierr);
      ierr = VecDuplicate(_x_,&CT_CCTm1_Cx); CHKERRQ(ierr);
      ierr = mField.VecScatterCreate(x,x_problem,Row,_x_,y_problem,Col,&scatter); CHKERRQ(ierr);
      PetscLogEventEnd(USER_EVENT_projInit,0,0,0,0);
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode Destroy() {
    PetscFunctionBegin;
    if(init) PetscFunctionReturn(0);
    PetscErrorCode ierr;
    ierr = MatDestroy(&CT); CHKERRQ(ierr);
    ierr = MatDestroy(&CCT); CHKERRQ(ierr);
    ierr = MatDestroy(&CTC); CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);
    ierr = VecDestroy(&_x_); CHKERRQ(ierr);
    ierr = VecDestroy(&Cx); CHKERRQ(ierr);
    ierr = VecDestroy(&CCTm1_Cx); CHKERRQ(ierr);
    ierr = VecDestroy(&CT_CCTm1_Cx); CHKERRQ(ierr);
    ierr = VecDestroy(&CTCx); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  friend PetscErrorCode matQ_mult_shell(Mat Q,Vec x,Vec f);
  friend PetscErrorCode matP_mult_shell(Mat P,Vec x,Vec f);
  friend PetscErrorCode matCTC_QTKQ_mult_shell(Mat CTC_QTKQ,Vec x,Vec f);
};

PetscErrorCode matQ_mult_shell(Mat Q,Vec x,Vec f) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  void *void_ctx;
  ierr = MatShellGetContext(Q,&void_ctx); CHKERRQ(ierr);
  matPROJ_ctx *ctx = (matPROJ_ctx*)void_ctx;
  PetscLogEventBegin(ctx->USER_EVENT_projQ,0,0,0,0);
  ierr = ctx->Init(x); CHKERRQ(ierr);
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
  PetscLogEventEnd(ctx->USER_EVENT_projQ,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode matP_mult_shell(Mat P,Vec x,Vec f) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  void *void_ctx;
  ierr = MatShellGetContext(P,&void_ctx); CHKERRQ(ierr);
  matPROJ_ctx *ctx = (matPROJ_ctx*)void_ctx;
  PetscLogEventBegin(ctx->USER_EVENT_projP,0,0,0,0);
  ierr = ctx->Init(x); CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatMult(ctx->C,ctx->_x_,ctx->Cx);  CHKERRQ(ierr);
  ierr = KSPSolve(ctx->ksp,ctx->Cx,ctx->CCTm1_Cx); CHKERRQ(ierr);
  ierr = MatMult(ctx->CT,ctx->CCTm1_Cx,ctx->CT_CCTm1_Cx);  CHKERRQ(ierr);
  ierr = VecZeroEntries(f); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx->scatter,ctx->CT_CCTm1_Cx,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,ctx->CT_CCTm1_Cx,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PetscLogEventEnd(ctx->USER_EVENT_projP,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode matCTC_QTKQ_mult_shell(Mat CTC_QTKQ,Vec x,Vec f) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  void *void_ctx;
  ierr = MatShellGetContext(CTC_QTKQ,&void_ctx); CHKERRQ(ierr);
  matPROJ_ctx *ctx = (matPROJ_ctx*)void_ctx;
  PetscLogEventBegin(ctx->USER_EVENT_projCTC_QTKQ,0,0,0,0);
  Mat Q;
  int M,N,m,n;
  ierr = MatGetSize(ctx->K,&M,&N); CHKERRQ(ierr);
  ierr = MatGetLocalSize(ctx->K,&m,&n); CHKERRQ(ierr);
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,ctx,&Q); CHKERRQ(ierr);
  ierr = MatMult(Q,x,f); CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatMult(ctx->CCT,ctx->_x_,ctx->CTCx); CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx->scatter,ctx->CTCx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,ctx->CTCx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = MatDestroy(&Q); CHKERRQ(ierr);
  PetscLogEventEnd(ctx->USER_EVENT_projCTC_QTKQ,0,0,0,0);
  PetscFunctionReturn(0);
}

#endif //__NONLINEAR_ELASTICITY_HPP__
