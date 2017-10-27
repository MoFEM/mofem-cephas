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

#ifndef __PROJECTION_MATRIX_CTX_HPP__
#define __PROJECTION_MATRIX_CTX_HPP__

/**
  * \brief structure for projection matrices
  * \ingroup projection_matrix
  *
  */
struct ConstrainMatrixCtx {

  MoFEM::Interface& mField;

  KSP kSP;
  Mat C,CT,CCT,CTC,K;
  Vec Cx,CCTm1_Cx,CT_CCTm1_Cx,CTCx;
  Vec X,Qx,KQx;
  bool initQorP,initQTKQ;
  bool createKSP;
  bool createScatter;
  bool cancelKSPMonitor;
  bool ownConstrainMatrix;

  // Scatter is created form problem_x to problem_y, or scatter is given
  // in the constructor

  VecScatter sCatter;
  string xProblem,yProblem;


  PetscLogEvent MOFEM_EVENT_projInit;
  PetscLogEvent MOFEM_EVENT_projQ;
  PetscLogEvent MOFEM_EVENT_projP;
  PetscLogEvent MOFEM_EVENT_projR;
  PetscLogEvent MOFEM_EVENT_projRT;
  PetscLogEvent MOFEM_EVENT_projCTC_QTKQ;

  /**
   * Construct data structure to build operators for projection matrices
   *
   * User need to set matrix C to make it work
   *
   * \param x_problem problem on which vector is projected
   * \param y_problem problem used to construct projection matrices
   * \param create_ksp create ksp solver otherwise  user need to set it up
   */
  ConstrainMatrixCtx(
    MoFEM::Interface& m_field,
    string x_problem,
    string y_problem,
    bool create_ksp = true,
    bool own_contrain_matrix = false
  );

  ConstrainMatrixCtx(
    MoFEM::Interface& m_field,
    VecScatter scatter,
    bool create_ksp = true,
    bool own_contrain_matrix = false
  );

  virtual ~ConstrainMatrixCtx() {
    ierr = destroyQorP(); CHKERRABORT(mField.get_comm(),ierr);
    ierr = destroyQTKQ(); CHKERRABORT(mField.get_comm(),ierr);
    if(ownConstrainMatrix) {
      ierr = MatDestroy(&C); CHKERRABORT(mField.get_comm(),ierr);
    }
  };

  PetscReal rTol,absTol,dTol;
  PetscInt maxIts;

  /**
    * \brief initialize vectors and matrices for Q and P shell matrices, scattering is set based on x_problem and y_problem
    *
    * \param x is a vector from problem x
    */
  PetscErrorCode initializeQorP(Vec x);

  /**
    * \brief initialize vectors and matrices for CTC+QTKQ shell matrices, scattering is set based on x_problem and y_problem
    */
  PetscErrorCode initializeQTKQ();

  /**
    * \brief re-calculate CT and CCT if C matrix has been changed since initialization
    */
  PetscErrorCode recalculateCTandCCT();

  /**
    * \brief re-calculate CTC matrix has been changed since initialization
    */
  PetscErrorCode recalculateCTC();

  /**
    * \brief destroy sub-matrices used for shell matrices P, Q, R, RT
    */
  PetscErrorCode destroyQorP();

  /**
    * \brief destroy sub-matrices used for shell matrix QTKQ
    */
  PetscErrorCode destroyQTKQ();

  friend PetscErrorCode ProjectionMatrixMultOpQ(Mat Q,Vec x,Vec f);
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
  * Mat Q; //for problem
  * ConstrainMatrixCtx projection_matrix_ctx(m_field,problem_name,contrains_problem_name);
  * ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,&projection_matrix_ctx,&Q); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))ProjectionMatrixMultOpQ); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(Q,MATOP_DESTROY,(void(*)(void))ConstrainMatrixDestroyOpPorQ); CHKERRQ(ierr);
  *
  * \endcode

  * \ingroup projection_matrix

  */
PetscErrorCode ProjectionMatrixMultOpQ(Mat Q,Vec x,Vec f);

/**
 * \deprecated Use ProjectionMatrixMultOpQ
 */
DEPRECATED PetscErrorCode PorjectionMatrixMultOpQ(Mat Q, Vec x, Vec f) {
  return ProjectionMatrixMultOpQ(Q, x, f);
}

/**
  * \brief Multiplication operator for P = CT(CCT)^-1C
  *
  * \code
  * Mat P; //for problem
  * ConstrainMatrixCtx
  projection_matrix_ctx(m_field,problem_name,contrains_problem_name);
  * ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,&projection_matrix_ctx,&P);
  CHKERRQ(ierr);
  * ierr =
  MatShellSetOperation(P,MATOP_MULT,(void(*)(void))ConstrainMatrixMultOpP);
  CHKERRQ(ierr);
  *
  * \endcode

  * \ingroup projection_matrix

  */
PetscErrorCode ConstrainMatrixMultOpP(Mat P, Vec x, Vec f);

/**
  * \brief Multiplication operator for R = CT(CCT)^-1
  *
  * \code
  * Mat R; //for problem
  * ConstrainMatrixCtx projection_matrix_ctx(m_field,problem_name,contrains_problem_name);
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
  * Mat RT; //for problem
  * ConstrainMatrixCtx projection_matrix_ctx(m_field,problem_name,contrains_problem_name);
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
  * Mat CTC_QTKQ; //for problem
  * ConstrainMatrixCtx projection_matrix_ctx(m_field,problem_name,contrains_problem_name);
  * ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,&projection_matrix_ctx,&CTC_QTKQ); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(CTC_QTKQ,MATOP_MULT,(void(*)(void))ConstrainMatrixMultOpCTC_QTKQ); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(CTC_QTKQ,MATOP_DESTROY,(void(*)(void))ConstrainMatrixDestroyOpQTKQ); CHKERRQ(ierr);
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
  * Mat Q; //for problem
  * ConstrainMatrixCtx projection_matrix_ctx(m_field,problem_name,contrains_problem_name);
  * ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,&projection_matrix_ctx,&Q); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))ProjectionMatrixMultOpQ); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(Q,MATOP_DESTROY,(void(*)(void))ConstrainMatrixDestroyOpPorQ); CHKERRQ(ierr);
  *
  * \endcode

  * \ingroup projection_matrix

  */
PetscErrorCode ConstrainMatrixDestroyOpPorQ(Mat Q);

/**
  * \brief Destroy shell matrix
  *
  * \code
  * Mat CTC_QTKQ; //for problem
  * ConstrainMatrixCtx projection_matrix_ctx(m_field,problem_name,contrains_problem_name);
  * ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,&projection_matrix_ctx,&Q); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))ConstrainMatrixMultOpCTC_QTKQ); CHKERRQ(ierr);
  * ierr = MatShellSetOperation(Q,MATOP_DESTROY,(void(*)(void))mat_destroy_QTKQ); CHKERRQ(ierr);
  *
  * \endcode

  * \ingroup projection_matrix

  */
PetscErrorCode ConstrainMatrixDestroyOpQTKQ(Mat QTKQ);

#endif // __PROJECTION_MATRIX_CTX_HPP__

/***************************************************************************//**
 \defgroup projection_matrix Constrain Projection Matrix
 \ingroup user_modules
 ******************************************************************************/
