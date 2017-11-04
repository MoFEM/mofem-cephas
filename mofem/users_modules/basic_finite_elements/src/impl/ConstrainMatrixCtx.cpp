/** \file ConstrainMatrixCtx.cpp
 * \brief Implementation of projection matrix
 *
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
#include <ConstrainMatrixCtx.hpp>

const static bool debug = false;

#define INIT_DATA_CONSTRAINMATRIXCTX \
C(PETSC_NULL), \
CT(PETSC_NULL),  \
CCT(PETSC_NULL), \
CTC(PETSC_NULL), \
K(PETSC_NULL), \
Cx(PETSC_NULL), \
CCTm1_Cx(PETSC_NULL), \
CT_CCTm1_Cx(PETSC_NULL), \
CTCx(PETSC_NULL), \
Qx(PETSC_NULL), \
KQx(PETSC_NULL), \
initQorP(true), \
initQTKQ(true), \
createKSP(create_ksp), \
createScatter(true), \
cancelKSPMonitor(true), \
ownConstrainMatrix(own_contrain_matrix)

ConstrainMatrixCtx::ConstrainMatrixCtx(
  MoFEM::Interface& m_field,string x_problem,string y_problem,bool create_ksp,bool own_contrain_matrix
):
mField(m_field),
INIT_DATA_CONSTRAINMATRIXCTX,
xProblem(x_problem),
yProblem(y_problem)
{
  PetscLogEventRegister("ProjectionInit",0,&MOFEM_EVENT_projInit);
  PetscLogEventRegister("ProjectionQ",0,&MOFEM_EVENT_projQ);
  PetscLogEventRegister("ProjectionP",0,&MOFEM_EVENT_projP);
  PetscLogEventRegister("ProjectionR",0,&MOFEM_EVENT_projR);
  PetscLogEventRegister("ProjectionRT",0,&MOFEM_EVENT_projRT);
  PetscLogEventRegister("ProjectionCTC_QTKQ",0,&MOFEM_EVENT_projCTC_QTKQ);
}

ConstrainMatrixCtx::ConstrainMatrixCtx(
  MoFEM::Interface& m_field,
  VecScatter scatter,
  bool create_ksp,
  bool own_contrain_matrix
):
mField(m_field),
INIT_DATA_CONSTRAINMATRIXCTX,
sCatter(scatter) {
  PetscLogEventRegister("ProjectionInit",0,&MOFEM_EVENT_projInit);
  PetscLogEventRegister("ProjectionQ",0,&MOFEM_EVENT_projQ);
  PetscLogEventRegister("ProjectionP",0,&MOFEM_EVENT_projP);
  PetscLogEventRegister("ProjectionR",0,&MOFEM_EVENT_projR);
  PetscLogEventRegister("ProjectionRT",0,&MOFEM_EVENT_projRT);
  PetscLogEventRegister("ProjectionCTC_QTKQ",0,&MOFEM_EVENT_projCTC_QTKQ);
}


MoFEMErrorCode ConstrainMatrixCtx::initializeQorP(Vec x) {
    MoFEMFunctionBeginHot;
    if(initQorP) {
      initQorP = false;
      
      PetscLogEventBegin(MOFEM_EVENT_projInit,0,0,0,0);
      ierr = MatTranspose(C,MAT_INITIAL_MATRIX,&CT); CHKERRG(ierr);
      ierr = MatTransposeMatMult(CT,CT,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CCT); CHKERRG(ierr); // need to be calculated when C is changed
      if(createKSP) {
        ierr = KSPCreate(mField.get_comm(),&kSP); CHKERRG(ierr); // neet to be recalculated when C is changed
        ierr = KSPSetOperators(kSP,CCT,CCT); CHKERRG(ierr);
        ierr = KSPSetFromOptions(kSP); CHKERRG(ierr);
        ierr = KSPSetInitialGuessKnoll(kSP,PETSC_TRUE); CHKERRG(ierr);
        ierr = KSPGetTolerances(kSP,&rTol,&absTol,&dTol,&maxIts); CHKERRG(ierr);
        ierr = KSPSetUp(kSP); CHKERRG(ierr);
        if(cancelKSPMonitor) {
          ierr = KSPMonitorCancel(kSP); CHKERRG(ierr);
        }
      }
      #if PETSC_VERSION_GE(3,5,3)
      ierr = MatCreateVecs(C,&X,PETSC_NULL); CHKERRG(ierr);
      ierr = MatCreateVecs(C,PETSC_NULL,&Cx); CHKERRG(ierr);
      ierr = MatCreateVecs(CCT,PETSC_NULL,&CCTm1_Cx); CHKERRG(ierr);
      #else
      ierr = MatGetVecs(C,&X,PETSC_NULL); CHKERRG(ierr);
      ierr = MatGetVecs(C,PETSC_NULL,&Cx); CHKERRG(ierr);
      ierr = MatGetVecs(CCT,PETSC_NULL,&CCTm1_Cx); CHKERRG(ierr);
      #endif
      ierr = VecDuplicate(X,&CT_CCTm1_Cx); CHKERRG(ierr);
      if(createScatter) {
        ierr = mField.getInterface<VecManager>()->vecScatterCreate(x,xProblem,ROW,X,yProblem,COL,&sCatter); CHKERRG(ierr);
      }
      PetscLogEventEnd(MOFEM_EVENT_projInit,0,0,0,0);
    }
    MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ConstrainMatrixCtx::recalculateCTandCCT() {
    MoFEMFunctionBeginHot;
    if(initQorP) MoFEMFunctionReturnHot(0);
    
    ierr = MatTranspose(C,MAT_REUSE_MATRIX,&CT); CHKERRG(ierr);
    ierr = MatTransposeMatMult(CT,CT,MAT_REUSE_MATRIX,PETSC_DEFAULT,&CCT); CHKERRG(ierr);
    MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ConstrainMatrixCtx::destroyQorP() {
    MoFEMFunctionBeginHot;
    if(initQorP) MoFEMFunctionReturnHot(0);
    
    ierr = MatDestroy(&CT); CHKERRG(ierr);
    ierr = MatDestroy(&CCT); CHKERRG(ierr);
    if(createKSP) {
      ierr = KSPDestroy(&kSP); CHKERRG(ierr);
    }
    ierr = VecDestroy(&X); CHKERRG(ierr);
    ierr = VecDestroy(&Cx); CHKERRG(ierr);
    ierr = VecDestroy(&CCTm1_Cx); CHKERRG(ierr);
    ierr = VecDestroy(&CT_CCTm1_Cx); CHKERRG(ierr);
    if(createScatter) {
      ierr = VecScatterDestroy(&sCatter); CHKERRG(ierr);
    }
    initQorP = true;
    MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ConstrainMatrixCtx::initializeQTKQ() {
    MoFEMFunctionBeginHot;
    if(initQTKQ) {
      initQTKQ = false;
      
      PetscLogEventBegin(MOFEM_EVENT_projInit,0,0,0,0);
      ierr = MatTransposeMatMult(C,C,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CTC); CHKERRG(ierr); // need to be recalculated when C is changed
      if(debug) {
        //MatView(CCT,PETSC_VIEWER_DRAW_WORLD);
        int m,n;
        MatGetSize(CCT,&m,&n);
        PetscPrintf(mField.get_comm(),"CCT size (%d,%d)\n",m,n);
        //std::string wait;
        //std::cin >> wait;
      }
      #if PETSC_VERSION_GE(3,5,3)
      ierr = MatCreateVecs(K,&Qx,PETSC_NULL); CHKERRG(ierr);
      ierr = MatCreateVecs(K,PETSC_NULL,&KQx); CHKERRG(ierr);
      ierr = MatCreateVecs(CTC,PETSC_NULL,&CTCx); CHKERRG(ierr);
      #else
      ierr = MatGetVecs(K,&Qx,PETSC_NULL); CHKERRG(ierr);
      ierr = MatGetVecs(K,PETSC_NULL,&KQx); CHKERRG(ierr);
      ierr = MatGetVecs(CTC,PETSC_NULL,&CTCx); CHKERRG(ierr);
      #endif
      PetscLogEventEnd(MOFEM_EVENT_projInit,0,0,0,0);
    }
    MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ConstrainMatrixCtx::recalculateCTC() {
    MoFEMFunctionBeginHot;
    if(initQTKQ) MoFEMFunctionReturnHot(0);
    
    ierr = MatTransposeMatMult(C,C,MAT_REUSE_MATRIX,PETSC_DEFAULT,&CTC); CHKERRG(ierr);
    MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ConstrainMatrixCtx::destroyQTKQ() {
    MoFEMFunctionBeginHot;
    if(initQTKQ) MoFEMFunctionReturnHot(0);
    
    ierr = MatDestroy(&CTC); CHKERRG(ierr);
    ierr = VecDestroy(&Qx); CHKERRG(ierr);
    ierr = VecDestroy(&KQx); CHKERRG(ierr);
    ierr = VecDestroy(&CTCx); CHKERRG(ierr);
    initQTKQ = true;
    MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ProjectionMatrixMultOpQ(Mat Q,Vec x,Vec f) {
  MoFEMFunctionBeginHot;
  
  void *void_ctx;
  ierr = MatShellGetContext(Q,&void_ctx); CHKERRG(ierr);
  ConstrainMatrixCtx *ctx = (ConstrainMatrixCtx*)void_ctx;
  PetscLogEventBegin(ctx->MOFEM_EVENT_projQ,0,0,0,0);
  ierr = ctx->initializeQorP(x); CHKERRG(ierr);
  ierr = VecCopy(x,f); CHKERRG(ierr);
  ierr = VecScatterBegin(ctx->sCatter,x,ctx->X,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecScatterEnd(ctx->sCatter,x,ctx->X,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  if(debug) {
    //ierr = VecView(ctx->X,PETSC_VIEWER_STDOUT_WORLD); CHKERRG(ierr);
    ierr = VecScatterBegin(ctx->sCatter,ctx->X,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
    ierr = VecScatterEnd(ctx->sCatter,ctx->X,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
    PetscBool  flg;
    ierr = VecEqual(x,f,&flg); CHKERRG(ierr);
    if(flg ==  PETSC_FALSE) SETERRQ(PETSC_COMM_SELF,1,"scatter is not working");
  }
  ierr = MatMult(ctx->C,ctx->X,ctx->Cx);  CHKERRG(ierr);
  ierr = KSPSolve(ctx->kSP,ctx->Cx,ctx->CCTm1_Cx); CHKERRG(ierr);
  ierr = MatMult(ctx->CT,ctx->CCTm1_Cx,ctx->CT_CCTm1_Cx);  CHKERRG(ierr);
  ierr = VecScale(ctx->CT_CCTm1_Cx,-1); CHKERRG(ierr);
  ierr = VecScatterBegin(ctx->sCatter,ctx->CT_CCTm1_Cx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  ierr = VecScatterEnd(ctx->sCatter,ctx->CT_CCTm1_Cx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  PetscLogEventEnd(ctx->MOFEM_EVENT_projQ,0,0,0,0);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ConstrainMatrixMultOpP(Mat P,Vec x,Vec f) {
  MoFEMFunctionBeginHot;
  
  void *void_ctx;
  ierr = MatShellGetContext(P,&void_ctx); CHKERRG(ierr);
  ConstrainMatrixCtx *ctx = (ConstrainMatrixCtx*)void_ctx;
  PetscLogEventBegin(ctx->MOFEM_EVENT_projP,0,0,0,0);
  ierr = ctx->initializeQorP(x); CHKERRG(ierr);
  ierr = VecScatterBegin(ctx->sCatter,x,ctx->X,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecScatterEnd(ctx->sCatter,x,ctx->X,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = MatMult(ctx->C,ctx->X,ctx->Cx);  CHKERRG(ierr);
  ierr = KSPSolve(ctx->kSP,ctx->Cx,ctx->CCTm1_Cx); CHKERRG(ierr);
  ierr = MatMult(ctx->CT,ctx->CCTm1_Cx,ctx->CT_CCTm1_Cx);  CHKERRG(ierr);
  ierr = VecZeroEntries(f); CHKERRG(ierr);
  ierr = VecGhostUpdateBegin(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecScatterBegin(ctx->sCatter,ctx->CT_CCTm1_Cx,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  ierr = VecScatterEnd(ctx->sCatter,ctx->CT_CCTm1_Cx,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  PetscLogEventEnd(ctx->MOFEM_EVENT_projP,0,0,0,0);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ConstrainMatrixMultOpR(Mat R,Vec x,Vec f) {
  MoFEMFunctionBeginHot;
  
  void *void_ctx;
  ierr = MatShellGetContext(R,&void_ctx); CHKERRG(ierr);
  ConstrainMatrixCtx *ctx = (ConstrainMatrixCtx*)void_ctx;
  PetscLogEventBegin(ctx->MOFEM_EVENT_projR,0,0,0,0);
  if(ctx->initQorP) SETERRQ(PETSC_COMM_SELF,1,"you have to call first initQorP or use Q matrix");
  ierr = KSPSolve(ctx->kSP,x,ctx->CCTm1_Cx); CHKERRG(ierr);
  ierr = MatMult(ctx->CT,ctx->CCTm1_Cx,ctx->CT_CCTm1_Cx);  CHKERRG(ierr);
  ierr = VecZeroEntries(f); CHKERRG(ierr);
  ierr = VecGhostUpdateBegin(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecScatterBegin(ctx->sCatter,ctx->CT_CCTm1_Cx,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  ierr = VecScatterEnd(ctx->sCatter,ctx->CT_CCTm1_Cx,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  PetscLogEventEnd(ctx->MOFEM_EVENT_projR,0,0,0,0);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ConstrainMatrixMultOpRT(Mat RT,Vec x,Vec f) {
  MoFEMFunctionBeginHot;
  
  void *void_ctx;
  ierr = MatShellGetContext(RT,&void_ctx); CHKERRG(ierr);
  ConstrainMatrixCtx *ctx = (ConstrainMatrixCtx*)void_ctx;
  PetscLogEventBegin(ctx->MOFEM_EVENT_projRT,0,0,0,0);
  if(ctx->initQorP) SETERRQ(PETSC_COMM_SELF,1,"you have to call first initQorP or use Q matrix");
  ierr = VecScatterBegin(ctx->sCatter,x,ctx->X,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecScatterEnd(ctx->sCatter,x,ctx->X,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = MatMult(ctx->C,ctx->X,ctx->Cx);  CHKERRG(ierr);
  ierr = KSPSolve(ctx->kSP,ctx->Cx,f); CHKERRG(ierr);
  PetscLogEventEnd(ctx->MOFEM_EVENT_projRT,0,0,0,0);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ConstrainMatrixMultOpCTC_QTKQ(Mat CTC_QTKQ,Vec x,Vec f) {
  MoFEMFunctionBeginHot;
  
  void *void_ctx;
  ierr = MatShellGetContext(CTC_QTKQ,&void_ctx); CHKERRG(ierr);
  ConstrainMatrixCtx *ctx = (ConstrainMatrixCtx*)void_ctx;
  PetscLogEventBegin(ctx->MOFEM_EVENT_projCTC_QTKQ,0,0,0,0);
  Mat Q;
  int M,N,m,n;
  ierr = MatGetSize(ctx->K,&M,&N); CHKERRG(ierr);
  ierr = MatGetLocalSize(ctx->K,&m,&n); CHKERRG(ierr);
  ierr = MatCreateShell(ctx->mField.get_comm(),m,n,M,N,ctx,&Q); CHKERRG(ierr);
  ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))ProjectionMatrixMultOpQ); CHKERRG(ierr);
  ierr = ctx->initializeQTKQ(); CHKERRG(ierr);
  ierr = MatMult(Q,x,ctx->Qx); CHKERRG(ierr);
  ierr = MatMult(ctx->K,ctx->Qx,ctx->KQx); CHKERRG(ierr);
  ierr = MatMult(Q,ctx->KQx,f); CHKERRG(ierr);
  ierr = VecScatterBegin(ctx->sCatter,x,ctx->X,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecScatterEnd(ctx->sCatter,x,ctx->X,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = MatMult(ctx->CTC,ctx->X,ctx->CTCx); CHKERRG(ierr);
  ierr = VecScatterBegin(ctx->sCatter,ctx->CTCx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  ierr = VecScatterEnd(ctx->sCatter,ctx->CTCx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  ierr = MatDestroy(&Q); CHKERRG(ierr);
  PetscLogEventEnd(ctx->MOFEM_EVENT_projCTC_QTKQ,0,0,0,0);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ConstrainMatrixDestroyOpPorQ(Mat Q) {
  MoFEMFunctionBeginHot;
  
  void *void_ctx;
  ierr = MatShellGetContext(Q,&void_ctx); CHKERRG(ierr);
  ConstrainMatrixCtx *ctx = (ConstrainMatrixCtx*)void_ctx;
  ierr = ctx->destroyQorP(); CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode ConstrainMatrixDestroyOpQTKQ(Mat QTKQ) {
  MoFEMFunctionBeginHot;
  
  void *void_ctx;
  ierr = MatShellGetContext(QTKQ,&void_ctx); CHKERRG(ierr);
  ConstrainMatrixCtx *ctx = (ConstrainMatrixCtx*)void_ctx;
  ierr = ctx->destroyQTKQ(); CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}
