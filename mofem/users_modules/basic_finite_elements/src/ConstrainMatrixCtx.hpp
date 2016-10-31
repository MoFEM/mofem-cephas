/** \file ConstrainMatrixCtx.hpp
 *
 * Can be used if constrains are linear, i.e. are not function of time.
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

/**
  * \brief structure for projection matrices
  * \ingroup projection_matrix
  *
  */
struct ConstrainMatrixCtx {
  MoFEM::Interface& mField;
  KSP kSP;
  VecScatter sCatter;
  Mat C,CT,CCT,CTC,K;
  Vec Cx,CCTm1_Cx,CT_CCTm1_Cx,CTCx;
  Vec X,Qx,KQx;
  string xProblem,yProblem;
  bool initQorP,initQTKQ;
  bool createKSP,createScatter;

  PetscLogEvent USER_EVENT_projInit;
  PetscLogEvent USER_EVENT_projQ;
  PetscLogEvent USER_EVENT_projP;
  PetscLogEvent USER_EVENT_projR;
  PetscLogEvent USER_EVENT_projRT;
  PetscLogEvent USER_EVENT_projCTC_QTKQ;

  ConstrainMatrixCtx(
    MoFEM::Interface& m_field,
    string x_problem,
    string y_problem,
    bool create_ksp = true
  );

  ConstrainMatrixCtx(
    MoFEM::Interface& m_field,
    VecScatter scatter,
    bool create_ksp = true
  );

  PetscReal rTol,absTol,dTol;
  PetscInt maxIts;

  /**
    * \brief initialize vectors and matrices for Q and P shell matrices, stacttering is set based on x_problem and y_problem
    */
  PetscErrorCode initializeQorP(Vec x);

  /**
    * \brief initialize vectors and matrices for CTC+QTKQ shell matrices, stacttering is set based on x_problem and y_problem
    */
  PetscErrorCode initializeQTKQ();

  /**
    * \brief recalculete CT and CCT if C matrix has been changed since initalisation
    */
  PetscErrorCode recalculateCTandCCT();

  /**
    * \brief recalculete CTC matrix has been changed since initalisation
    */
  PetscErrorCode recalculateCTC();

  /**
    * \brief destroy submatrices sused for shell marices P, Q, R, RT
    */
  PetscErrorCode destroyQorP();

  /**
    * \brief destroy submatrices sused for shell marix QTKQ
    */
  PetscErrorCode destroyQTKQ();

  friend PetscErrorCode PorjectionMatrixMultOpQ(Mat Q,Vec x,Vec f);
  friend PetscErrorCode ConstrainMatrixMultOpP(Mat P,Vec x,Vec f);
  friend PetscErrorCode ConstrainMatrixMultOpR(Mat R,Vec x,Vec f);
  friend PetscErrorCode ConstrainMatrixMultOpRT(Mat RT,Vec x,Vec f);
  friend PetscErrorCode ConstrainMatrixMultOpCTC_QTKQ(Mat CTC_QTKQ,Vec x,Vec f);

  friend PetscErrorCode ConstrainMatrixDestroyOpPorQ();
  friend PetscErrorCode ConstrainMatrixDestroyOpQTKQ();

};

/**
  * \brief Multiplication operator for Q = I-CTC(CCT)^-1C
  *
  * \code
  * Mat Q; //for poroblem
  * ConstrainMatrixCtx projection_matrix_ctx(m_fiel,problem_name,contrains_porblem_name);
  * ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,&projection_matrix_ctx,&Q); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))PorjectionMatrixMultOpQ); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(Q,(void(*)(void))mat_destroy_PorQ); CHKERRQ(ierr);
  *
  * \endcode

  * \ingroup projection_matrix

  */
PetscErrorCode PorjectionMatrixMultOpQ(Mat Q,Vec x,Vec f);

/**
  * \brief Multiplication operator for P = CT(CCT)^-1C
  *
  * \code
  * Mat P; //for poroblem
  * ConstrainMatrixCtx projection_matrix_ctx(m_fiel,problem_name,contrains_porblem_name);
  * ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,&projection_matrix_ctx,&P); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(P,MATOP_MULT,(void(*)(void))ConstrainMatrixMultOpP); CHKERRQ(ierr);
  *
  * \endcode

  * \ingroup projection_matrix

  */
PetscErrorCode ConstrainMatrixMultOpP(Mat P,Vec x,Vec f);

/**
  * \brief Multiplication operator for R = CT(CCT)^-1
  *
  * \code
  * Mat R; //for poroblem
  * ConstrainMatrixCtx projection_matrix_ctx(m_fiel,problem_name,contrains_porblem_name);
  * ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,&projection_matrix_ctx,&R); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(R,MATOP_MULT,(void(*)(void))ConstrainMatrixMultOpR); CHKERRQ(ierr);
  *
  * \endcode

  * \ingroup projection_matrix

  */
PetscErrorCode ConstrainMatrixMultOpR(Mat R,Vec x,Vec f);

/**
  * \brief Multiplication operator for RT = (CCT)^-TC
  *
  * \code
  * Mat RT; //for poroblem
  * ConstrainMatrixCtx projection_matrix_ctx(m_fiel,problem_name,contrains_porblem_name);
  * ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,&projection_matrix_ctx,&RT); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(RT,MATOP_MULT,(void(*)(void))ConstrainMatrixMultOpRT); CHKERRQ(ierr);
  *
  * \endcode

  * \ingroup projection_matrix

  */
PetscErrorCode ConstrainMatrixMultOpRT(Mat RT,Vec x,Vec f);

/**
  * \brief Multiplication operator for RT = (CCT)^-TC
  *
  * \code
  * Mat CTC_QTKQ; //for poroblem
  * ConstrainMatrixCtx projection_matrix_ctx(m_fiel,problem_name,contrains_porblem_name);
  * ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,&projection_matrix_ctx,&CTC_QTKQ); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(CTC_QTKQ,MATOP_MULT,(void(*)(void))ConstrainMatrixMultOpCTC_QTKQ); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(CTC_QTKQ,MATOP_DESTROY,(void(*)(void))mat_destroy_QTKQ); CHKERRQ(ierr);
  *
  * \endcode
  *

  * \ingroup projection_matrix

  */
PetscErrorCode ConstrainMatrixMultOpCTC_QTKQ(Mat CTC_QTKQ,Vec x,Vec f);

/**
  * \brief Destroy shell matrix Q
  *
  * \code
  * Mat Q; //for poroblem
  * ConstrainMatrixCtx projection_matrix_ctx(m_fiel,problem_name,contrains_porblem_name);
  * ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,&projection_matrix_ctx,&Q); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))PorjectionMatrixMultOpQ); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(Q,(void(*)(void))mat_destroy_PorQ); CHKERRQ(ierr);
  *
  * \endcode

  * \ingroup projection_matrix

  */
PetscErrorCode ConstrainMatrixDestroyOpPorQ(Mat Q);

/**
  * \brief Destroy shell matrix
  *
  * \code
  * Mat CTC_QTKQ; //for poroblem
  * ConstrainMatrixCtx projection_matrix_ctx(m_fiel,problem_name,contrains_porblem_name);
  * ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,&projection_matrix_ctx,&Q); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))ConstrainMatrixMultOpCTC_QTKQ); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(Q,MATOP_DESTROY,(void(*)(void))mat_destroy_QTKQ); CHKERRQ(ierr);
  *
  * \endcode

  * \ingroup projection_matrix

  */
PetscErrorCode ConstrainMatrixDestroyOpQTKQ(Mat QTKQ);

#endif //__PROJECTIONMATRIXCTX_HPP__


/***************************************************************************//**
 \defgroup projection_matrix Constrain Projection Matrix
 \ingroup user_modules
 ******************************************************************************/
