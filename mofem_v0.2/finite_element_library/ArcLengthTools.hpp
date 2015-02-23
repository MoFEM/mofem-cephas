/** \file ArcLengthTools.cpp
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

#ifndef __ARCLEGHTTOOLS_HPP__
#define __ARCLEGHTTOOLS_HPP__

namespace MoFEM {

/**
 * \brief Store variables for ArcLength analaysis
 *
 * r_lambda = f_lambda - s
 * f_lambda = alpha*f(dx*dx) + beta*(dlambda*sqrt(F_lambda*F_lambda) 
 *
 * dx = x-x0
 *
 * db*ddx + diag*ddlambda - r_lambda = 0
 *
 * User need to implement fuctions calulating f_lambda, i.e. f(dx*x) and
 * direvative, db
 * 
 * alpha,beta parameters
 * dlambda is load factor
 * s arc-length radius
 * F_lambda reference load vcetor
 * F_lambda2 dot product of F_lambda
 * diag value on matrix diagonal
 * x0  displacement vetor at begining of step
 * x current displacemengt vector
 * dx2 dot product of dx vector
 * db direvative of f(dx*dx), i.e. db = d[ f(dx*dx) ]/dx
 *
 * x_lambda is solution of eq. K*x_lambda = F_lambda
 */
struct ArcLengthCtx {

  ErrorCode rval;
  PetscErrorCode ierr;

  double s;	///< arc length radius
  double beta; 	///< force scaling factor 
  double alpha; ///< displacement scaling factor

  double dlambda;	///< increment of load factor
  double diag;		///< diagonal value
  double dx2;		///< inner_prod(dX,dX)
  double F_lambda2;	///< inner_prod(F_lambda,F_lambda);
  double res_lambda;	///< f_lambda - s
  Vec F_lambda;		///< F_lambda reference load vcetor
  Vec db;		///< db direvative of f(dx*dx), i.e. db = d[ f(dx*dx) ]/dx
  Vec x_lambda;		///< solution of eq. K*x_lambda = F_lambda
  Vec x0;		///< displacement vetor at begining of step
  Vec dx;		///< dx = x-x0

  /** 
    * \brief set arc radius
    */
  PetscErrorCode setS(double s);

  /**
   * \brief set parematers controling arc-length equaitions
   * alpha controls off diagonal therms
   * beta controls diagonal therm
   */
  PetscErrorCode setAlphaBeta(double alpha,double beta);

  ArcLengthCtx(FieldInterface &mField,const string &problem_name);
  ~ArcLengthCtx();

  NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dIt;
  DofIdx getPetscGloablDofIdx() { return dIt->get_petsc_gloabl_dof_idx(); };
  DofIdx getPetscLocalDofIdx() { return dIt->get_petsc_local_dof_idx(); };
  FieldData& getFieldData() { return dIt->get_FieldData(); }
  int getPart() { return dIt->get_part(); };

};

/**
 * \brief It is ctx structure passed to SNES solver
 */
struct ArcLengthSnesCtx: public SnesCtx {
  ErrorCode rval;
  PetscErrorCode ierr;
  ArcLengthCtx* arcPtr;
  ArcLengthSnesCtx(FieldInterface &m_field,const string &problem_name,ArcLengthCtx* arc_ptr):
    SnesCtx(m_field,problem_name),arcPtr(arc_ptr) {}
};

/** \brief shell matrix for arc-length method
 *
 * Shell matrix which has tructure
 * [ K 		-dF_lambda]
 * [ db		 diag	]
 */
struct ArcLengthMatShell {

  ErrorCode rval;
  PetscErrorCode ierr;
  FieldInterface& mField;
  Mat Aij;
  ArcLengthCtx* arcPtr;
  string problemName;

  ArcLengthMatShell(FieldInterface& _mField,Mat _Aij,ArcLengthCtx *arc_ptr,string _problem_name);
  virtual ~ArcLengthMatShell();

  PetscErrorCode set_lambda(Vec ksp_x,double *lambda,ScatterMode scattermode); 

  friend PetscErrorCode ArcLengthMatMultShellOp(Mat A,Vec x,Vec f);
};

/**
 * mult operator for Arc Length Shell Mat
 */
PetscErrorCode ArcLengthMatMultShellOp(Mat A,Vec x,Vec f);

/**
 * strutture for Arc Length precodnditioner
 */
struct PCShellCtx {
  PC pC;
  Mat shellAij,Aij;
  ArcLengthCtx* arcPtr;
  PCShellCtx(Mat shell_Aij,Mat _Aij,ArcLengthCtx* arc_ptr); 
  ~PCShellCtx();

  friend PetscErrorCode PCApplyArcLength(PC pC,Vec pc_f,Vec pc_x);
  friend PetscErrorCode PCSetupArcLength(PC pC);
};

/**
 * apply oppertor for Arc Length precoditionet
 * solves K*pc_x = pc_f
 * solves K*x_lambda = -dF_lambda
 * solves ddlambda = ( res_lambda - db*x_lambda )/( diag + db*pc_x )
 * calculate pc_x = pc_x + ddlambda*x_lambda
 */
PetscErrorCode PCApplyArcLength(PC pc,Vec pc_f,Vec pc_x);

/**
 * set up struture for Arc Length precoditioner
 * it sets precoditioner for matrix K
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
  
  FieldInterface& mField;
  ArcLengthCtx *arcPtr;
  
  PrePostProcessForArcLength(FieldInterface& _mField, ArcLengthCtx *arc_ptr);

  PetscErrorCode preProcess();    
  PetscErrorCode postProcess();

};

}

#endif // __ARCLEGHTTOOLS_HPP__


