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

#ifndef __ARCLEGHTTOOLS_HPP__
#define __ARCLEGHTTOOLS_HPP__

#include "FieldInterface.hpp"
#include "SnesCtx.hpp"

namespace MoFEM {

/**
 * Store variables for ArcLength analaysis
 *
 * r_lambda = f_lambda - s
 * f_lambda = alpha*f(dx*dx) + beta*(dlambda*sqrt(F_lambda*F_lambda) 
 *
 * dx = x-x0
 *
 * db*ddx + diag*ddlambda - r_lambda = 0
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

  double s,beta,alpha;
  PetscErrorCode set_s(double _s) { 
    PetscFunctionBegin;
    s = _s;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\tSet s = %6.4e\n",s); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
   * set parematers controling arc-length equaitions
   * alpha controls off diagonal therms
   * beta controls diagonal therm
   */
  PetscErrorCode set_alpha_and_beta(double _alpha,double _beta) { 
    PetscFunctionBegin;
    alpha = _alpha;
    beta = _beta;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\tSet alpha = %6.4e beta = %6.4e\n",alpha,beta); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  double dlambda;
  //dx2 - dot product of 
  double diag,dx2,F_lambda2,res_lambda;
  Vec F_lambda,db,x_lambda,x0,dx;

  NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit;
  DofIdx get_petsc_gloabl_dof_idx() { return dit->get_petsc_gloabl_dof_idx(); };
  DofIdx get_petsc_local_dof_idx() { return dit->get_petsc_local_dof_idx(); };
  FieldData& get_FieldData() { return dit->get_FieldData(); }
  int get_part() { return dit->get_part(); };

  ArcLengthCtx(FieldInterface &mField,const string &problem_name):
    dlambda(0),diag(0),dx2(0),F_lambda2(0),res_lambda(0) {
    PetscErrorCode ierr;

    ierr = mField.VecCreateGhost(problem_name,ROW,&F_lambda); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecSetOption(F_lambda,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = mField.VecCreateGhost(problem_name,ROW,&db); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = mField.VecCreateGhost(problem_name,ROW,&x_lambda); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = mField.VecCreateGhost(problem_name,ROW,&x0); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = mField.VecCreateGhost(problem_name,ROW,&dx); CHKERRABORT(PETSC_COMM_WORLD,ierr);

    const MoFEMProblem *problem_ptr;
    ierr = mField.get_problem(problem_name,&problem_ptr); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    NumeredDofMoFEMEntity_multiIndex& dofsPtr_no_const 
	    = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
    NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator hi_dit;
    dit = dofsPtr_no_const.get<FieldName_mi_tag>().lower_bound("LAMBDA");
    hi_dit = dofsPtr_no_const.get<FieldName_mi_tag>().upper_bound("LAMBDA");
    if(distance(dit,hi_dit)!=1) {
      PetscTraceBackErrorHandler(PETSC_COMM_WORLD,__LINE__,PETSC_FUNCTION_NAME,__FILE__,__SDIR__,1,PETSC_ERROR_INITIAL,
	"can not find unique LAMBDA (load factor)",PETSC_NULL);
      PetscMPIAbortErrorHandler(PETSC_COMM_WORLD,__LINE__,PETSC_FUNCTION_NAME,__FILE__,__SDIR__,1,PETSC_ERROR_INITIAL,
	"can not find unique LAMBDA (load factor)",PETSC_NULL);
    }

  }

  ~ArcLengthCtx() {
    PetscErrorCode ierr;
    ierr = VecDestroy(&F_lambda); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecDestroy(&db); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecDestroy(&x_lambda); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecDestroy(&x0); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecDestroy(&dx); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }


};

/**
 * It is ctx structure passed to SNES solver
 */
struct ArcLengthSnesCtx: public SnesCtx {

  ErrorCode rval;
  PetscErrorCode ierr;

  ArcLengthCtx* arc_ptr;
  ArcLengthSnesCtx(FieldInterface &_mField,const string &_problem_name,ArcLengthCtx* _arc_ptr):
    SnesCtx(_mField,_problem_name),arc_ptr(_arc_ptr) {}

};

/**
 * Shell matrix which has tructure
 * [ K 		-dF_lambda]
 * [ db		 diag	]
 */
struct ArcLengthMatShell {

  ErrorCode rval;
  PetscErrorCode ierr;


  FieldInterface& mField;

  Mat Aij;
  ArcLengthCtx* arc_ptr;
  string problem_name;
  ArcLengthMatShell(FieldInterface& _mField,Mat _Aij,ArcLengthCtx *_arc_ptr,string _problem_name): 
    mField(_mField),Aij(_Aij),arc_ptr(_arc_ptr),problem_name(_problem_name) {};
  PetscErrorCode set_lambda(Vec ksp_x,double *lambda,ScatterMode scattermode) {
    PetscFunctionBegin;
    const MoFEMProblem *problem_ptr;
    ierr = mField.get_problem(problem_name,&problem_ptr); CHKERRQ(ierr);
    if(arc_ptr->get_petsc_local_dof_idx()!=-1) {
      PetscScalar *array;
      ierr = VecGetArray(ksp_x,&array); CHKERRQ(ierr);
      switch(scattermode) {
	case SCATTER_FORWARD:
	  *lambda = array[arc_ptr->get_petsc_local_dof_idx()];
	  break;
	case SCATTER_REVERSE:
	  array[arc_ptr->get_petsc_local_dof_idx()] = *lambda;
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      ierr = VecRestoreArray(ksp_x,&array); CHKERRQ(ierr);
    } 
    int part = arc_ptr->get_part();
    MPI_Bcast(lambda,1,MPI_DOUBLE,part,PETSC_COMM_WORLD);

    PetscFunctionReturn(0);
  }
  ~ArcLengthMatShell() { }

  friend PetscErrorCode arc_length_mult_shell(Mat A,Vec x,Vec f);
};

/**
 * mult operator for Arc Length Shell Mat
 */
PetscErrorCode arc_length_mult_shell(Mat A,Vec x,Vec f);

/**
 * strutture for Arc Length precodnditioner
 */
struct PCShellCtx {
  PC pc;
  Mat ShellAij,Aij;
  ArcLengthCtx* arc_ptr;
  PCShellCtx(Mat _ShellAij,Mat _Aij,ArcLengthCtx* _arc_ptr): 
    ShellAij(_ShellAij),Aij(_Aij),arc_ptr(_arc_ptr) {
    PetscErrorCode ierr;
    ierr = PCCreate(PETSC_COMM_WORLD,&pc); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  ~PCShellCtx() {
    PetscErrorCode ierr;
    ierr = PCDestroy(&pc); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  friend PetscErrorCode pc_apply_arc_length(PC pc,Vec pc_f,Vec pc_x);
  friend PetscErrorCode pc_setup_arc_length(PC pc);
};

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
struct PrePostProcessFEMethod_For_F_lambda: public FieldInterface::FEMethod {
  
  FieldInterface& mField;
  ArcLengthCtx *arc_ptr;
  
  PrePostProcessFEMethod_For_F_lambda(FieldInterface& _mField, ArcLengthCtx *_arc_ptr):
    mField(_mField),arc_ptr(_arc_ptr) {}

  PetscErrorCode ierr;
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      
      switch(snes_ctx) {
        case CTX_SNESNONE:
        case CTX_SNESSETFUNCTION: {
          //F_lambda
          ierr = VecZeroEntries(arc_ptr->F_lambda); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        break;
        case CTX_SNESSETJACOBIAN: {
        }
        break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      
      switch(snes_ctx) {
        case CTX_SNESNONE: {
        }
        case CTX_SNESSETFUNCTION: {
          //F_lambda
          ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
          //F_lambda2
          ierr = VecDot(arc_ptr->F_lambda,arc_ptr->F_lambda,&arc_ptr->F_lambda2); CHKERRQ(ierr);
          PetscPrintf(PETSC_COMM_WORLD,"\tFlambda2 = %6.4e\n",arc_ptr->F_lambda2);
        }
        break;
        case CTX_SNESSETJACOBIAN: {
        }
        break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      
      PetscFunctionReturn(0);
    }
};

/**
 * apply oppertor for Arc Length precoditionet
 * solves K*pc_x = pc_f
 * solves K*x_lambda = -dF_lambda
 * solves ddlambda = ( res_lambda - db*x_lambda )/( diag + db*pc_x )
 * calulate pc_x = pc_x + ddlambda*x_lambda
 */
PetscErrorCode pc_apply_arc_length(PC pc,Vec pc_f,Vec pc_x);

/**
 * set up struture for Arc Length precoditioner
 * it sets precoditioner for matrix K
 */
PetscErrorCode pc_setup_arc_length(PC pc);

};

#endif // __ARCLEGHTTOOLS_HPP__


