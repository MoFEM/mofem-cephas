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

#ifndef __PROJECTIONMATRIXCTX_HPP__
#define __PROJECTIONMATRIXCTX_HPP__

namespace ObosleteUsersModules {

/**
  * \brief structure for projection matries
  * \ingroup projection_matrix
  *
  */
struct ProjectionMatrixCtx {
  FieldInterface& mField;
  KSP kSP;
  Vec X;
  VecScatter sCatter;
  Mat CT,CCT,CTC;
  Vec Cx,CCTm1_Cx,CT_CCTm1_Cx,CTCx;
  Vec Qx,KQx;
  string xProblem,yProblem;
  bool initQorP,initQTKQ;

  PetscLogEvent USER_EVENT_projInit;
  PetscLogEvent USER_EVENT_projQ;
  PetscLogEvent USER_EVENT_projP;
  PetscLogEvent USER_EVENT_projR;
  PetscLogEvent USER_EVENT_projRT;
  PetscLogEvent USER_EVENT_projCTC_QTKQ;
  ProjectionMatrixCtx(FieldInterface& _mField,string x_problem,string y_problem): 
    mField(_mField),
    CT(PETSC_NULL),CCT(PETSC_NULL),CTC(PETSC_NULL),
    Cx(PETSC_NULL),CCTm1_Cx(PETSC_NULL),CT_CCTm1_Cx(PETSC_NULL),CTCx(PETSC_NULL),
    Qx(PETSC_NULL),KQx(PETSC_NULL),
    xProblem(x_problem),yProblem(y_problem),
    initQorP(true),initQTKQ(true),
    C(PETSC_NULL),K(PETSC_NULL),g(PETSC_NULL) {
    PetscLogEventRegister("ProjectionInit",0,&USER_EVENT_projInit);
    PetscLogEventRegister("ProjectionQ",0,&USER_EVENT_projQ);
    PetscLogEventRegister("ProjectionP",0,&USER_EVENT_projP);
    PetscLogEventRegister("ProjectionR",0,&USER_EVENT_projR);
    PetscLogEventRegister("ProjectionRT",0,&USER_EVENT_projRT);
    PetscLogEventRegister("ProjectionCTC_QTKQ",0,&USER_EVENT_projCTC_QTKQ);
  }


  PetscReal rTol,absTol,dTol;
  PetscInt maxIts;
  Mat C,K;
  Vec g;

  /**
    * \brief Init vectors and matrices for Q and P shell matrices, stacttering is set based on x_problem and y_problem
    */
  PetscErrorCode initializeQorP(Vec x);
  PetscErrorCode initializeQTKQ();


  PetscErrorCode recalculateCTandCCT();
  PetscErrorCode destroyQorP();
  PetscErrorCode recalculateCTC();
  PetscErrorCode destroyQTKQ();

  friend PetscErrorCode matQ_mult_shell(Mat Q,Vec x,Vec f);
  friend PetscErrorCode matP_mult_shell(Mat P,Vec x,Vec f);
  friend PetscErrorCode matR_mult_shell(Mat R,Vec x,Vec f);
  friend PetscErrorCode matRT_mult_shell(Mat R,Vec x,Vec f);
  friend PetscErrorCode matCTC_QTKQ_mult_shell(Mat CTC_QTKQ,Vec x,Vec f);

};

PetscErrorCode matQ_mult_shell(Mat Q,Vec x,Vec f);
PetscErrorCode matP_mult_shell(Mat P,Vec x,Vec f);
PetscErrorCode matR_mult_shell(Mat R,Vec x,Vec f);
PetscErrorCode matRT_mult_shell(Mat RT,Vec x,Vec f);
PetscErrorCode matCTC_QTKQ_mult_shell(Mat CTC_QTKQ,Vec x,Vec f);

}

#endif //__PROJECTIONMATRIXCTX_HPP__


/***************************************************************************//**
 * \defgroup projection_matrix Projection Matrix
 ******************************************************************************/

