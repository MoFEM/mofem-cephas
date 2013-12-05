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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"

namespace MoFEM {

/**
  * \brief structure for projection matries
  *
  */
struct matPROJ_ctx {
  FieldInterface& mField;
  KSP ksp,ksp_self;
  Vec _x_;
  VecScatter scatter;
  Mat CT,CCT,CTC;
  Vec Cx,CCTm1_Cx,CT_CCTm1_Cx,CTCx;
  Vec Qx,KQx;
  string x_problem,y_problem;
  bool initQorP,initQTKQ;
  bool debug;
  PetscLogEvent USER_EVENT_projInit;
  PetscLogEvent USER_EVENT_projQ;
  PetscLogEvent USER_EVENT_projP;
  PetscLogEvent USER_EVENT_projR;
  PetscLogEvent USER_EVENT_projRT;
  PetscLogEvent USER_EVENT_projCTC_QTKQ;
  matPROJ_ctx(FieldInterface& _mField,string _x_problem,string _y_problem): 
    mField(_mField),
    CT(PETSC_NULL),CCT(PETSC_NULL),CTC(PETSC_NULL),
    Cx(PETSC_NULL),CCTm1_Cx(PETSC_NULL),CT_CCTm1_Cx(PETSC_NULL),CTCx(PETSC_NULL),
    Qx(PETSC_NULL),KQx(PETSC_NULL),
    x_problem(_x_problem),y_problem(_y_problem),
    initQorP(true),initQTKQ(true),debug(true),
    C(PETSC_NULL),K(PETSC_NULL),g(PETSC_NULL) {
    PetscLogEventRegister("ProjectionInit",0,&USER_EVENT_projInit);
    PetscLogEventRegister("ProjectionQ",0,&USER_EVENT_projQ);
    PetscLogEventRegister("ProjectionP",0,&USER_EVENT_projP);
    PetscLogEventRegister("ProjectionR",0,&USER_EVENT_projR);
    PetscLogEventRegister("ProjectionRT",0,&USER_EVENT_projRT);
    PetscLogEventRegister("ProjectionCTC_QTKQ",0,&USER_EVENT_projCTC_QTKQ);
  }


  PetscReal rtol,abstol,dtol;
  PetscInt maxits;
  Mat C,K;
  Vec g;

  /**
    * \brief Init vectors and matrices for Q and P shell matrices, stacttering is set based on x_problem and y_problem
    */
  PetscErrorCode InitQorP(Vec x);

  PetscErrorCode RecalculateCTandCCT();
  PetscErrorCode DestroyQorP();
  PetscErrorCode InitQTKQ();
  PetscErrorCode RecalulateCTC();
  PetscErrorCode DestroyQTKQ();

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

#endif //__MOABCONSTRAINSBYMARKAINSWORTH_HPP__

