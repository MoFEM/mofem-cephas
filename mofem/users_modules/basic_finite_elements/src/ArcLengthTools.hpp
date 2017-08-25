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
 *
 * \ingroup arc_length_control

 The constrain function if given by
 \f[
 r_\lambda = f_\lambda(\mathbf{x},\lambda) - s^2
 \f]
 where \f$f_\lambda(\mathbf{x},\lambda)\f$ is some constrain function, which has general form given by
 \f[
 f_\lambda(\mathbf{x},\lambda) = \alpha f(\|\Delta\mathbf{x}\|^2) +
 \beta^2 \Delta\lambda^2 \| \mathbf{F}_{\lambda} \|^2
 \f]
 where \f$f(\|\Delta\mathbf{x}\|^2)\f$ is some user defined function evaluating
 increments vector of degrees of freedom  \f$\Delta\mathbf{x}\f$. The increment
 vector is
 \f[
 \Delta \mathbf{x} = \mathbf{x}-\mathbf{x}_0.
 \f]
 For convenience we assume that
 \f[
 \frac{\partial f}{\partial \mathbf{x}}\Delta \mathbf{x}
 =
 \textrm{d}\mathbf{b} \Delta \mathbf{x},
 \f]
 as result linearized constrain equation takes form
 \f[
 \textrm{d}\mathbf{b} \delta \Delta x +
 D \delta \Delta\lambda - r_{\lambda} = 0
 \f]
 where
 \f[
 D = 2\beta^2 \Delta\lambda \| \mathbf{F}_{\lambda} \|^2.
 \f]

 User need to implement functions calculating \f$f(\mathbf{x},\lambda)\f$, i.e. function
 \f$f(\|\Delta\mathbf{x}\|^2)\f$  and its derivative, \f$\textrm{d}\mathbf{b}\f$.

 */
struct ArcLengthCtx {

  MoFEM::Interface &mField;

  double s;	///< arc length radius
  double beta; 	///< force scaling factor
  double alpha; ///< displacement scaling factor

  Vec ghosTdLambda;
  double dLambda;	///< increment of load factor
  Vec ghostDiag;
  double dIag;		///< diagonal value

  double dx2;		///< inner_prod(dX,dX)
  double F_lambda2;	///< inner_prod(F_lambda,F_lambda);
  double res_lambda;	///< f_lambda - s
  Vec F_lambda;		///< F_lambda reference load vector
  Vec db;		///< db derivative of f(dx*dx), i.e. db = d[ f(dx*dx) ]/dx
  Vec xLambda;		///< solution of eq. K*xLambda = F_lambda
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

  ArcLengthCtx(
    MoFEM::Interface& m_field,
    const std::string& problem_name,
    const std::string& field_name = "LAMBDA"
  );
  virtual ~ArcLengthCtx();

  NumeredDofEntityByFieldName::iterator dIt;

  /** \brief Get global index of load factor
  */
  DofIdx getPetscGlobalDofIdx() { return (*dIt)->getPetscGlobalDofIdx(); };

  /** \brief Get local index of load factor
  */
  DofIdx getPetscLocalDofIdx() { return (*dIt)->getPetscLocalDofIdx(); };

  /** \brief Get value of load factor
  */
  FieldData& getFieldData() { return (*dIt)->getFieldData(); }

  /** \brief Get proc owning lambda dof
  */
  int getPart() { return (*dIt)->getPart(); };

};

#ifdef __SNESCTX_HPP__

/**
 * \brief It is ctx structure passed to SNES solver
 * \ingroup arc_length_control
 */
struct ArcLengthSnesCtx: public SnesCtx {
  ArcLengthCtx* arcPtr;
  ArcLengthSnesCtx(
    MoFEM::Interface &m_field,const std::string &problem_name,ArcLengthCtx* arc_ptr
  ):
  SnesCtx(m_field,problem_name),
  arcPtr(arc_ptr) {
  }
};

#endif //__SNESCTX_HPP__

#ifdef __TSCTX_HPP__

/**
 * \brief It is ctx structure passed to SNES solver
 * \ingroup arc_length_control
 */
struct ArcLengthTsCtx: public TsCtx {
  ArcLengthCtx* arcPtr;
  ArcLengthTsCtx(
    MoFEM::Interface &m_field,const std::string &problem_name,ArcLengthCtx* arc_ptr
  ):
  TsCtx(m_field,problem_name),
  arcPtr(arc_ptr) {
  }
};

#endif // __TSCTX_HPP__

/** \brief shell matrix for arc-length method
 *
 * \ingroup arc_length_control

 Shell matrix which has structure:
 \f[
 \left[
  \begin{array}{cc}
   \mathbf{K} & -\mathbf{F}_\lambda \\
   \textrm{d}\mathbf{b} & D
  \end{array}
  \right]
 \left\{
 \begin{array}{c}
 \delta \Delta \mathbf{x} \\
 \delta \Delta \lambda
 \end{array}
 \right\}
 =
 \left[
  \begin{array}{c}
    -\mathbf{f}_\textrm{int} \\
    -r_\lambda
  \end{array}
  \right]
 \f]

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
  PCArcLengthCtx(
    Mat shell_Aij,Mat _Aij,ArcLengthCtx* arc_ptr
  );
  virtual ~PCArcLengthCtx();

  friend PetscErrorCode PCApplyArcLength(PC pc,Vec pc_f,Vec pc_x);
  friend PetscErrorCode PCSetupArcLength(PC pc);
};

/**
 * apply operator for Arc Length pre-conditioner
 * solves K*pc_x = pc_f
 * solves K*xLambda = -dF_lambda
 * solves ddlambda = ( res_lambda - db*xLambda )/( diag + db*pc_x )
 * calculate pc_x = pc_x + ddlambda*xLambda
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

/** \brief Implementation of spherical arc-length method
  * \ingroup arc_length_control

  \f[
  \alpha \| \Delta\mathbf{x} \|^2
  + \Delta\lambda^2 \beta^2 \| \mathbf{F}_\lambda \|^2
  = s^2
  \f]

  This is particular implementation of ArcLength control, i.e. spherical arc
  length control. If beta is set to 0 and alpha is non-zero it is cylindrical
  arc-length control. Works well with general problem with non-linear
  geometry. It not guarantee dissipative loading path in case of physical
  nonlinearities.

  */
struct SphericalArcLengthControl: public FEMethod {

  ArcLengthCtx* arcPtr;

  SphericalArcLengthControl(ArcLengthCtx *arc_ptr);
  virtual ~SphericalArcLengthControl();

  PetscErrorCode preProcess();
  PetscErrorCode operator()();
  PetscErrorCode postProcess();

  /** \brief Calculate f_lambda(dx,lambda)

  \f[
  f_\lambda(\Delta\mathbf{x},\lambda) =
  \alpha \| \Delta\mathbf{x} \|^2
  + \Delta\lambda^2 \beta^2 \| \mathbf{F}_\lambda \|^2
  \f]

  */
  virtual double calculateLambdaInt();

  /** \brief Calculate db

  \f[
  \textrm{d}\mathbf{b} = 2 \alpha \Delta\mathbf{x}
  \f]

  */
  virtual PetscErrorCode calculateDb();
  virtual PetscErrorCode calculateDxAndDlambda(Vec x);
  virtual PetscErrorCode calculateInitDlambda(double *dlambda);
  virtual PetscErrorCode setDlambdaToX(Vec x,double dlambda);

};

#endif // __ARCLEGHTTOOLS_HPP__

/***************************************************************************//**
 \defgroup arc_length_control Arc-Length control
 \ingroup user_modules
 ******************************************************************************/
