/** \file ArcLengthTools.cpp
 *
 * Implementation of arc-length control method
 *
 * FIXME: Some variables not comply with naming convention, need to be fixed.
 *
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

#ifndef __ARCLEGHTTOOLS_HPP__
#define __ARCLEGHTTOOLS_HPP__

/**
 * \brief Store variables for ArcLength analysis
 * \ingroup arc_length_control
 *
 * r_lambda = f_lambda - s
 * f_lambda = alpha*f(dx*dx) + beta^2*(dlambda*sqrt(F_lambda*F_lambda)
 *
 * dx = x-x0
 *
 * db*ddx + diag*ddlambda - r_lambda = 0
 *
 * User need to implement functions calculating f_lambda, i.e. f(dx*dx) and
 * derivative, db
 *
 * alpha,beta parameters
 * dlambda is load factor
 * s arc-length radius
 * F_lambda reference load vector
 * F_lambda2 dot product of F_lambda
 * diag value on matrix diagonal
 * x0  displacement vector at beginning of step
 * x current displacement vector
 * dx2 dot product of dx vector
 * db derivative of f(dx*dx), i.e. db = d[ f(dx*dx) ]/dx
 *
 * x_lambda is solution of eq. K*x_lambda = F_lambda
 */
struct ArcLengthCtx {

  FieldInterface &mField;

  double s;	///< arc length radius
  double beta; 	///< force scaling factor
  double alpha; ///< displacement scaling factor

  double dlambda;	///< increment of load factor
  double diag;		///< diagonal value
  double dx2;		///< inner_prod(dX,dX)
  double F_lambda2;	///< inner_prod(F_lambda,F_lambda);
  double res_lambda;	///< f_lambda - s
  Vec F_lambda;		///< F_lambda reference load vector
  Vec db;		///< db derivative of f(dx*dx), i.e. db = d[ f(dx*dx) ]/dx
  Vec x_lambda;		///< solution of eq. K*x_lambda = F_lambda
  Vec x0;		///< displacement vector at beginning of step
  Vec dx;		///< dx = x-x0

  /**
    * \brief set arc radius
    */
  PetscErrorCode setS(double s);

  /**
   * \brief set parameters controlling arc-length equations
   * alpha controls off diagonal therms
   * beta controls diagonal therm
   */
  PetscErrorCode setAlphaBeta(double alpha,double beta);

  ArcLengthCtx(FieldInterface &m_field,const string &problem_name);
  ~ArcLengthCtx();

  NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dIt;
  DofIdx getPetscGloablDofIdx() { return dIt->get_petsc_gloabl_dof_idx(); };
  DofIdx getPetscLocalDofIdx() { return dIt->get_petsc_local_dof_idx(); };
  FieldData& getFieldData() { return dIt->get_FieldData(); }
  int getPart() { return dIt->get_part(); };

};

/**
 * \brief It is ctx structure passed to SNES solver
 * \ingroup arc_length_control
 */
struct ArcLengthSnesCtx: public SnesCtx {
  ArcLengthCtx* arcPtr;
  ArcLengthSnesCtx(FieldInterface &m_field,const string &problem_name,ArcLengthCtx* arc_ptr):
    SnesCtx(m_field,problem_name),arcPtr(arc_ptr) {}
};

/** \brief shell matrix for arc-length method
 * \ingroup arc_length_control
 *
 * Shell matrix which has structure
 * [ K 		-dF_lambda]
 * [ db		 diag	]
 */
struct ArcLengthMatShell {

  Mat Aij;
  ArcLengthCtx* arcPtr;
  string problemName;

  ArcLengthMatShell(Mat aij,ArcLengthCtx *arc_ptr,string problem_name);
  virtual ~ArcLengthMatShell();

  PetscErrorCode setLambda(Vec ksp_x,double *lambda,ScatterMode scattermode);

  friend PetscErrorCode ArcLengthMatMultShellOp(Mat A,Vec x,Vec f);
};

/**
 * mult operator for Arc Length Shell Mat
 */
PetscErrorCode ArcLengthMatMultShellOp(Mat A,Vec x,Vec f);

/**
 * \brief structure for Arc Length pre-conditioner
 * \ingroup arc_length_control
 */
struct PCArcLengthCtx {
  PC pC;
  Mat shellAij,Aij;
  ArcLengthCtx* arcPtr;
  PCArcLengthCtx(Mat shell_Aij,Mat _Aij,ArcLengthCtx* arc_ptr);
  ~PCArcLengthCtx();

  friend PetscErrorCode PCApplyArcLength(PC pc,Vec pc_f,Vec pc_x);
  friend PetscErrorCode PCSetupArcLength(PC pc);
};

/**
 * apply operator for Arc Length pre-conditioner
 * solves K*pc_x = pc_f
 * solves K*x_lambda = -dF_lambda
 * solves ddlambda = ( res_lambda - db*x_lambda )/( diag + db*pc_x )
 * calculate pc_x = pc_x + ddlambda*x_lambda
 */
PetscErrorCode PCApplyArcLength(PC pc,Vec pc_f,Vec pc_x);

/**
 * set up structure for Arc Length pre-conditioner

 * it sets pre-conditioner for matrix K
 */
PetscErrorCode PCSetupArcLength(PC pc);

/**
 * \brief Pre and Post Process for Arc Length
 * preProcess - zero F_lambda
 * postProcess - assembly F_lambda
 * Example: \code
      SnesCtx::basic_method_to_do& preProcess_to_do_Rhs = SnesCtx.get_preProcess_to_do_Rhs();
      SnesCtx::basic_method_to_do& postProcess_to_do_Rhs = SnesCtx.get_postProcess_to_do_Rhs();
      SnesCtx.get_preProcess_to_do_Rhs().push_back(&PrePostFE); //Zero F_lambda before looping over FEs
      loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ELASTIC",&MyFE));
      loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("INTERFACE",&IntMyFE));
      loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&MyArcMethod));
      SnesCtx.get_postProcess_to_do_Rhs().push_back(&PrePostFE); //finally, assemble F_lambda
  \endcode
 */
struct PrePostProcessForArcLength: public FEMethod {

  ArcLengthCtx *arcPtr;
  PrePostProcessForArcLength(ArcLengthCtx *arc_ptr);

  PetscErrorCode preProcess();
  PetscErrorCode postProcess();

};

/** \brief Implementation of cylindrical arc-length method, i.e. alpha*dx*x + dlmabda^2*beta^2*F_lamda*F_lambda = s^2
  * \ingroup arc_length_control
  *
  * This is particular implementation of ArcLength control, i.e. spherical arc
  * length control. If beta is set to 0 and alpha is non-zero it is cylindrical
  * arc-length control. Works well with general problem with non-linear
  * geometry. It not guarantee dissipative loading path in case of physical
  * nonlinearities.
  *
  */
struct SphericalArcLengthControl: public FEMethod {

  ArcLengthCtx* arcPtr;
  Vec ghostDiag;

  SphericalArcLengthControl(ArcLengthCtx *arc_ptr);
  ~SphericalArcLengthControl();

  PetscErrorCode preProcess();
  PetscErrorCode operator()();
  PetscErrorCode postProcess();

  /**
    * alpha*dx*x + dlmabda^2*beta^2*F_lamda*F_lambda
    */
  double calculateLambdaInt();

  /**
    * db(x)/dx  = 2*dx
    */
  PetscErrorCode calculateDb();
  PetscErrorCode calculateDxAndDlambda(Vec x);
  PetscErrorCode calculateInitDlambda(double *dlambda);
  PetscErrorCode setDlambdaToX(Vec x,double dlambda);

};

#endif // __ARCLEGHTTOOLS_HPP__

/***************************************************************************//**
 \defgroup arc_length_control Arc-Length control
 \ingroup user_modules
 ******************************************************************************/
