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

#ifndef __ARC_LENGHT_NONLINEAR_ELASTICITY_HPP__
#define __ARC_LENGHT_NONLINEAR_ELASTICITY_HPP__

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include <petscksp.h>

#include "moabSnes.hpp"
#include "PostProcDisplacementOnMesh.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "moabSnes.hpp"
#include "moabFEMethod_ComplexForLazy.hpp"
#include "moabFEMethod_DriverComplexForLazy.hpp"

#include "complex_for_lazy.h"

#include "nonlinear_elasticity.hpp"

#ifdef __cplusplus
extern "C" {
#endif
#include <petsc-private/snesimpl.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

ErrorCode rval;
PetscErrorCode ierr;

struct ArcLenghtCtx {

  double s,beta,alpha;
  PetscErrorCode set_s(double _s) { 
    PetscFunctionBegin;
    s = _s;
    PetscFunctionReturn(0);
  }

PetscErrorCode set_alpha_and_beta(double _alpha,double _beta) { 
    PetscFunctionBegin;
    alpha = _alpha;
    beta = _beta;
    PetscFunctionReturn(0);
  }


  double dlambda;
  PetscErrorCode set_dlambda(double _dlambda) {
    PetscFunctionBegin;
    dlambda = _dlambda;
    PetscFunctionReturn(0);
  } 


  double diag,dx2,F_lambda2,res_lambda,;
  Vec F_lambda,db,x_lambda,y_residual,x0,dx;
  ArcLenghtCtx(moabField &mField,const string &problem_name) {

    mField.VecCreateGhost(problem_name,Row,&F_lambda);
    mField.VecCreateGhost(problem_name,Row,&db);
    mField.VecCreateGhost(problem_name,Row,&x_lambda);
    mField.VecCreateGhost(problem_name,Row,&y_residual);
    mField.VecCreateGhost(problem_name,Row,&x0);
    mField.VecCreateGhost(problem_name,Row,&dx);

  }

  ~ArcLenghtCtx() {
    VecDestroy(&F_lambda);
    VecDestroy(&db);
    VecDestroy(&x_lambda);
    VecDestroy(&y_residual);
    VecDestroy(&x0);
  }


};

struct MySnesCtx: public moabSnesCtx {

  ArcLenghtCtx* arc_ptr;
  MySnesCtx(moabField &_mField,const string &_problem_name,ArcLenghtCtx* _arc_ptr):
    moabSnesCtx(_mField,_problem_name),arc_ptr(_arc_ptr) {}

};

struct MatShellCtx {

  double scale_lambda;
  moabField& mField;

  Mat Aij;
  ArcLenghtCtx* arc_ptr;
  MatShellCtx(moabField& _mField,Mat _Aij,ArcLenghtCtx *_arc_ptr): scale_lambda(1),mField(_mField),Aij(_Aij),arc_ptr(_arc_ptr) {};
  PetscErrorCode set_lambda(Vec ksp_x,double *lambda,ScatterMode scattermode) {
    PetscFunctionBegin;
    const MoFEMProblem *problem_ptr;
    ierr = mField.get_problems_database("ELASTIC_MECHANICS",&problem_ptr); CHKERRQ(ierr);
    //get problem dofs
    NumeredDofMoFEMEntity_multiIndex &numered_dofs_rows = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
    NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit;
    dit = numered_dofs_rows.get<FieldName_mi_tag>().find("LAMBDA");
    DofIdx lambda_dof_index = dit->get_petsc_local_dof_idx();
    int part = dit->part;
    if(lambda_dof_index!=-1) {
      PetscScalar *array;
      ierr = VecGetArray(ksp_x,&array); CHKERRQ(ierr);
      switch(scattermode) {
	case SCATTER_FORWARD:
	  *lambda = array[lambda_dof_index];
	  break;
	case SCATTER_REVERSE:
	  array[lambda_dof_index] = *lambda;
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      ierr = VecRestoreArray(ksp_x,&array); CHKERRQ(ierr);
    } 
    MPI_Bcast(lambda,1,MPI_DOUBLE,part,PETSC_COMM_WORLD);

    PetscFunctionReturn(0);
  }
  ~MatShellCtx() { }

  friend PetscErrorCode arc_lenght_mult_shell(Mat A,Vec x,Vec f);
};
PetscErrorCode arc_lenght_mult_shell(Mat A,Vec x,Vec f) {
  PetscFunctionBegin;
  void *void_ctx;
  ierr = MatShellGetContext(A,&void_ctx); CHKERRQ(ierr);
  MatShellCtx *ctx = (MatShellCtx*)void_ctx;
  ierr = MatMult(ctx->Aij,x,f); CHKERRQ(ierr);
  double lambda;
  ierr = ctx->set_lambda(x,&lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  double db_dot_x;
  ierr = VecDot(ctx->arc_ptr->db,x,&db_dot_x); CHKERRQ(ierr);
  double f_lambda;
  f_lambda = ctx->scale_lambda*(ctx->arc_ptr->diag*lambda + db_dot_x);
  ierr = ctx->set_lambda(f,&f_lambda,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAXPY(f,-lambda,ctx->arc_ptr->F_lambda); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

struct PCShellCtx {
  PC pc;
  Mat ShellAij,Aij;
  ArcLenghtCtx* arc_ptr;
  PCShellCtx(Mat _ShellAij,Mat _Aij,ArcLenghtCtx* _arc_ptr): 
    ShellAij(_ShellAij),Aij(_Aij),arc_ptr(_arc_ptr) {
    PCCreate(PETSC_COMM_WORLD,&pc);
  }
  ~PCShellCtx() {
    PCDestroy(&pc);
  }
  friend PetscErrorCode pc_apply_arc_length(PC pc,Vec pc_f,Vec pc_x);
  friend PetscErrorCode pc_setup_arc_length(PC pc);
};
PetscErrorCode pc_apply_arc_length(PC pc,Vec pc_f,Vec pc_x) {
  PetscFunctionBegin;
  void *void_ctx;
  ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
  PCShellCtx *PCCtx = (PCShellCtx*)void_ctx;
  void *void_MatCtx;
  MatShellGetContext(PCCtx->ShellAij,&void_MatCtx);
  MatShellCtx *MatCtx = (MatShellCtx*)void_MatCtx;
  ierr = PCApply(PCCtx->pc,pc_f,pc_x); CHKERRQ(ierr);
  ierr = PCApply(PCCtx->pc,PCCtx->arc_ptr->F_lambda,PCCtx->arc_ptr->x_lambda); CHKERRQ(ierr);
  double db_dot_pc_x,db_dot_x_lambda;
  ierr = VecDot(PCCtx->arc_ptr->db,pc_x,&db_dot_pc_x); CHKERRQ(ierr);
  ierr = VecDot(PCCtx->arc_ptr->db,PCCtx->arc_ptr->x_lambda,&db_dot_x_lambda); CHKERRQ(ierr);
  double denominator = MatCtx->scale_lambda*PCCtx->arc_ptr->diag+MatCtx->scale_lambda*db_dot_x_lambda;
  double res_lambda;
  ierr = MatCtx->set_lambda(pc_f,&res_lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  double ddlambda = (res_lambda - MatCtx->scale_lambda*db_dot_pc_x)/denominator;
  //cerr << res_lambda << " " << ddlambda << " " << db_dot_pc_x << " " << db_dot_x_lambda << " " << PCCtx->arc_ptr->diag << endl;
  if(ddlambda != ddlambda) SETERRQ(PETSC_COMM_SELF,1,"problem with constraint");
  ierr = VecAXPY(pc_x,ddlambda,PCCtx->arc_ptr->x_lambda); CHKERRQ(ierr);
  ierr = MatCtx->set_lambda(pc_x,&ddlambda,SCATTER_REVERSE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode pc_setup_arc_length(PC pc) {
  PetscFunctionBegin;
  void *void_ctx;
  ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
  PCShellCtx *ctx = (PCShellCtx*)void_ctx;
  ierr = PCSetFromOptions(ctx->pc); CHKERRQ(ierr);
  MatStructure flag;
  ierr = PCGetOperators(pc,&ctx->ShellAij,&ctx->Aij,&flag); CHKERRQ(ierr);
  ierr = PCSetOperators(ctx->pc,ctx->ShellAij,ctx->Aij,flag); CHKERRQ(ierr);
  ierr = PCSetUp(ctx->pc); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

struct MyElasticFEMethod: public FEMethod_DriverComplexForLazy {

  ArcLenghtCtx* arc_ptr;
  Range& SideSet1;
  Range& SideSet2;
  Range& SideSet3;
  Range& SideSet4;
  Range& NodeSet1;

  Range SideSet1_;
  Range SideSet3_;
  Range SideSet4_;

  MyElasticFEMethod(Interface& _moab,double _lambda,double _mu,
      ArcLenghtCtx *_arc_ptr,
      Range &_SideSet1,Range &_SideSet2,
      Range &_SideSet3,Range &_SideSet4,
      Range &_NodeSet1,
      int _verbose = 0): 
      FEMethod_DriverComplexForLazy(_moab,_lambda,_mu,_verbose), 
      arc_ptr(_arc_ptr),
      SideSet1(_SideSet1),SideSet2(_SideSet2),
      SideSet3(_SideSet3),SideSet4(_SideSet4),
      NodeSet1(_NodeSet1) {

    set_PhysicalEquationNumber(neohookean);

    Range SideSetEdges,SideSetNodes;
    rval = moab.get_connectivity(SideSet1,SideSetNodes,true); CHKERR_THROW(rval);
    SideSet1_.insert(SideSet1.begin(),SideSet1.end());
    SideSet1_.insert(SideSetNodes.begin(),SideSetNodes.end());

    SideSetEdges.clear();
    SideSetNodes.clear();
    rval = moab.get_adjacencies(SideSet3,1,false,SideSetEdges,Interface::UNION); CHKERR_THROW(rval);
    rval = moab.get_connectivity(SideSet3,SideSetNodes,true); CHKERR_THROW(rval);
    SideSet3_.insert(SideSet3.begin(),SideSet3.end());
    SideSet3_.insert(SideSetEdges.begin(),SideSetEdges.end());
    SideSet3_.insert(SideSetNodes.begin(),SideSetNodes.end());

    SideSetEdges.clear();
    SideSetNodes.clear();
    rval = moab.get_adjacencies(SideSet4,1,false,SideSetEdges,Interface::UNION); CHKERR_THROW(rval);
    rval = moab.get_connectivity(SideSet4,SideSetNodes,true); CHKERR_THROW(rval);
    SideSet4_.insert(SideSet4.begin(),SideSet4.end());
    SideSet4_.insert(SideSetEdges.begin(),SideSetEdges.end());
    SideSet4_.insert(SideSetNodes.begin(),SideSetNodes.end());

  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy::preProcess(); CHKERRQ(ierr);
    switch(ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	//F_lambda
      	ierr = VecZeroEntries(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }


  PetscErrorCode operator()() {
    PetscFunctionBegin;

    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndices(); CHKERRQ(ierr);

    Range& DirihletSideSet = SideSet1_;
    Range& NeumannSideSet = SideSet2;
    Range& SymmBC_Y = SideSet3_;
    Range& SymmBC_X = SideSet4_;

    ierr = ApplyDirihletBC(DirihletSideSet,fixed_z,true); CHKERRQ(ierr);
    //cerr << DirihletBC.size() << endl;
    ierr = ApplyDirihletBC(SymmBC_Y,fixed_y,false); CHKERRQ(ierr);
    //cerr << DirihletBC.size() << endl;
    ierr = ApplyDirihletBC(SymmBC_X,fixed_x,false); CHKERRQ(ierr);
    //cerr << DirihletBC.size() << endl;
    if(Diagonal!=PETSC_NULL) {
	if(DirihletBC.size()>0) {
	  DirihletBCDiagVal.resize(DirihletBC.size());
	  fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
	  ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
	}
    }

    double t[] = { 0,0,t_val, 0,0,t_val, 0,0,t_val };

    switch(ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	  ierr = CalculateFint(snes_f); CHKERRQ(ierr);
	  ierr = CaluclateFext(arc_ptr->F_lambda,t,NeumannSideSet); CHKERRQ(ierr);
	}
	break;
      case ctx_SNESSetJacobian: {
	  ierr = CalculateTangent(*snes_B); CHKERRQ(ierr);
	  FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
	  dit = row_multiIndex->get<FieldName_mi_tag>().lower_bound("LAMBDA");
	  hi_dit = row_multiIndex->get<FieldName_mi_tag>().upper_bound("LAMBDA");
	  if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  double _lambda_ = dit->get_FieldData();
	  cblas_dscal(9,_lambda_,t,1);
	  ierr = CalculateTangentExt(*snes_B,t,NeumannSideSet); CHKERRQ(ierr);
	}
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    switch(ctx) {
      case ctx_SNESSetFunction: {  
      }
      break;
      case ctx_SNESSetJacobian: {
      }
      break;
      default:
      break;
    }
 
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(ctx) {
      case ctx_SNESNone:
	//snes_f //assemble only if ctx_SNESNone
	ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      case ctx_SNESSetFunction: { 
	//F_lambda
	ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	//F_lambda2
	ierr = VecDot(arc_ptr->F_lambda,arc_ptr->F_lambda,&arc_ptr->F_lambda2); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tFlambda^2 = %6.4e\n",arc_ptr->F_lambda2);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = VecAssemblyBegin(Diagonal); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(Diagonal); CHKERRQ(ierr);
	ierr = MatDiagonalSet(*snes_B,Diagonal,ADD_VALUES); CHKERRQ(ierr);
	ierr = VecDestroy(&Diagonal); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode potsProcessLoadPath() {
    PetscFunctionBegin;
    NumeredDofMoFEMEntity_multiIndex &numered_dofs_rows = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
    NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator lit;
    lit = numered_dofs_rows.get<FieldName_mi_tag>().find("LAMBDA");
    if(lit == numered_dofs_rows.get<FieldName_mi_tag>().end()) PetscFunctionReturn(0);
    Range::iterator nit = NodeSet1.begin();
    for(;nit!=NodeSet1.end();nit++) {
      NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator dit,hi_dit;
      dit = numered_dofs_rows.get<MoABEnt_mi_tag>().lower_bound(*nit);
      hi_dit = numered_dofs_rows.get<MoABEnt_mi_tag>().upper_bound(*nit);
      for(;dit!=hi_dit;dit++) {
	PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e -> ",lit->get_name().c_str(),lit->get_dof_rank(),lit->get_FieldData());
	PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e\n",dit->get_name().c_str(),dit->get_dof_rank(),dit->get_FieldData());
      }
    }
    PetscFunctionReturn(0);
  }

};

struct ArcLenghtElemFEMethod: public moabField::FEMethod {

  ArcLenghtCtx* arc_ptr;
  Vec GhostDiag;
  ArcLenghtElemFEMethod(Interface& _moab,ArcLenghtCtx *_arc_ptr): FEMethod(_moab),arc_ptr(_arc_ptr) {
    PetscInt ghosts[1] = { 0 };
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm->rank() == 0) {
      VecCreateGhost(PETSC_COMM_WORLD,1,1,1,ghosts,&GhostDiag);
    } else {
      VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&GhostDiag);
    }
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    switch(ctx) {
      case ctx_SNESSetFunction: { 
	ierr = calulate_dx_and_dlambda(snes_x); CHKERRQ(ierr);
	ierr = calulate_db(); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
      }
      break;
      default:
      break;
    }
    PetscFunctionReturn(0);
  }

  double calulate_lambda_int() {
    PetscFunctionBegin;
    return arc_ptr->alpha*arc_ptr->dx2 + pow(arc_ptr->dlambda,2)*pow(arc_ptr->beta,2)*arc_ptr->F_lambda2;
    PetscFunctionReturn(0);
  }

  PetscErrorCode calulate_db() {
    PetscFunctionBegin;
    //db
    ierr = VecCopy(arc_ptr->dx,arc_ptr->db); CHKERRQ(ierr);
    ierr = VecScale(arc_ptr->db,2*arc_ptr->alpha); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;

    FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
    dit = row_multiIndex->get<FieldName_mi_tag>().lower_bound("LAMBDA");
    hi_dit = row_multiIndex->get<FieldName_mi_tag>().upper_bound("LAMBDA");
    //only one LAMBDA
    if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");

    switch(ctx) {
      case ctx_SNESSetFunction: {
	arc_ptr->res_lambda = calulate_lambda_int() - pow(arc_ptr->s,2);
	ierr = VecSetValue(snes_f,dit->get_petsc_gloabl_dof_idx(),arc_ptr->res_lambda,ADD_VALUES); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tres_lambda = %6.4e\n",arc_ptr->res_lambda);
      }
      break; 
      case ctx_SNESSetJacobian: {
	double diag = 2*arc_ptr->dlambda*pow(arc_ptr->beta,2)*arc_ptr->F_lambda2;
	ierr = VecSetValue(GhostDiag,0,diag,INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(*snes_B,dit->get_petsc_gloabl_dof_idx(),dit->get_petsc_gloabl_dof_idx(),1,ADD_VALUES); CHKERRQ(ierr);
      }
      break;
      default:
      break;
    }	
    
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(ctx) {
      case ctx_SNESSetFunction: { 
	//assemble
	ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	//add F_lambda
	NumeredDofMoFEMEntity_multiIndex& dofs_moabfield_no_const 
	    = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
	NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
	dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().lower_bound("LAMBDA");
	hi_dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().upper_bound("LAMBDA");
	if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	ierr = VecAXPY(snes_f,-dit->get_FieldData(),arc_ptr->F_lambda); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",dit->get_FieldData());  
	//snes_f norm
	double fnorm;
	ierr = VecNormBegin(snes_f,NORM_2,&fnorm); CHKERRQ(ierr);	
	ierr = VecNormEnd(snes_f,NORM_2,&fnorm);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tfnorm = %6.4e\n",fnorm);  
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(GhostDiag); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(GhostDiag); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(GhostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(GhostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	double *diag;
	ierr = VecGetArray(GhostDiag,&diag); CHKERRQ(ierr);
	arc_ptr->diag = *diag;
	ierr = VecRestoreArray(GhostDiag,&diag); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tdiag = %6.4e\n",arc_ptr->diag);
	//Matrix View
	//MatView(*snes_B,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
	//std::string wait;
	//std::cin >> wait;
      }
      break;
      default:
      break;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode calulate_dx_and_dlambda(Vec x) {
    PetscFunctionBegin;
    //dx
    ierr = VecCopy(x,arc_ptr->dx); CHKERRQ(ierr);
    ierr = VecAXPY(arc_ptr->dx,-1,arc_ptr->x0); CHKERRQ(ierr);
    //dlambda
    NumeredDofMoFEMEntity_multiIndex& dofs_moabfield_no_const 
	  = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
    NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
    dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().lower_bound("LAMBDA");
    hi_dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().upper_bound("LAMBDA");
    if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    if(dit->get_petsc_local_dof_idx()!=-1) {
	  double *array;
	  ierr = VecGetArray(arc_ptr->dx,&array); CHKERRQ(ierr);
	  arc_ptr->dlambda = array[dit->get_petsc_local_dof_idx()];
	  array[dit->get_petsc_local_dof_idx()] = 0;
	  ierr = VecRestoreArray(arc_ptr->dx,&array); CHKERRQ(ierr);
    }
    int part = dit->part;
    MPI_Bcast(&(arc_ptr->dlambda),1,MPI_DOUBLE,part,PETSC_COMM_WORLD);
    //dx2
    ierr = VecDot(arc_ptr->dx,arc_ptr->dx,&arc_ptr->dx2); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"\tdlambda = %6.4e dx2 = %6.4e\n",arc_ptr->dlambda,arc_ptr->dx2);
    PetscFunctionReturn(0);
  }

  PetscErrorCode calculate_init_dlambda(double *dlambda) {

      PetscFunctionBegin;

      *dlambda = sqrt(pow(arc_ptr->s,2)/(pow(arc_ptr->beta,2)*arc_ptr->F_lambda2));
      if(!(*dlambda == *dlambda)) {
	      ostringstream sss;
	      sss << "s " << arc_ptr->s << " " << arc_ptr->beta << " " << arc_ptr->F_lambda2;
	      SETERRQ(PETSC_COMM_SELF,1,sss.str().c_str());
      }

      PetscFunctionReturn(0);
  }

  PetscErrorCode set_dlambda_to_x(Vec x,double dlambda) {
      PetscFunctionBegin;

      NumeredDofMoFEMEntity_multiIndex& dofs_moabfield_no_const 
	    = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
      NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
      dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().lower_bound("LAMBDA");
      hi_dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().upper_bound("LAMBDA");
      if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");

      if(dit->get_petsc_local_dof_idx()!=-1) {
	    double *array;
	    ierr = VecGetArray(x,&array); CHKERRQ(ierr);
	    double lambda_old = array[dit->get_petsc_local_dof_idx()];
	    if(!(dlambda == dlambda)) {
	      ostringstream sss;
	      sss << "s " << arc_ptr->s << " " << arc_ptr->beta << " " << arc_ptr->F_lambda2;
	      SETERRQ(PETSC_COMM_SELF,1,sss.str().c_str());
	    }
	    array[dit->get_petsc_local_dof_idx()] = lambda_old + dlambda;
	    PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e, %6.4e (%6.4e)\n",
	      lambda_old, array[dit->get_petsc_local_dof_idx()], dlambda);
	    ierr = VecRestoreArray(x,&array); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
  }


};

PetscErrorCode snes_apply_arc_lenght(_p_SNES *snes,Vec x) {
  PetscFunctionBegin;
  
  PetscInt           lits;
  MatStructure       flg = DIFFERENT_NONZERO_PATTERN;
  Vec                Y,F;
  KSPConvergedReason kspreason;

  /*void *void_snes_ctx;
  ierr = SNESShellGetContext(snes,&void_snes_ctx); CHKERRQ(ierr);
  MySnesCtx *SnesCtx = (MySnesCtx*)void_snes_ctx;*/

  ierr = SNESSetUpMatrices(snes);CHKERRQ(ierr);

  snes->numFailures            = 0;
  snes->numLinearSolveFailures = 0;
  snes->reason                 = SNES_CONVERGED_ITERATING;
  snes->iter                   = 0;
  snes->norm                   = 0.0;

  //PC pc;
  PetscBool domainerror;
  PetscReal fnorm,ynorm,xnorm;
  PetscInt  maxits,i;

  xnorm = 0;
  maxits = snes->max_its;	/* maximum number of iterations */

  F = snes->vec_func;
  Y = snes->vec_sol_update;

  ierr = SNESComputeFunction(snes,x,F);CHKERRQ(ierr);
  if (snes->domainerror) {
    snes->reason = SNES_DIVERGED_FUNCTION_DOMAIN;
    PetscFunctionReturn(0);
  }

  if (!snes->norm_init_set) {
    ierr = VecNormBegin(F,NORM_2,&fnorm);CHKERRQ(ierr);	/* fnorm <- ||F||  */
    ierr = VecNormEnd(F,NORM_2,&fnorm);CHKERRQ(ierr);
    if (PetscIsInfOrNanReal(fnorm)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FP,"User provided compute function generated a Not-a-Number");
  } else {
    fnorm = snes->norm_init;
    snes->norm_init_set = PETSC_FALSE;
  }
  ierr = PetscObjectTakeAccess(snes);CHKERRQ(ierr);
  snes->norm = fnorm;
  ierr = PetscObjectGrantAccess(snes);CHKERRQ(ierr);
  SNESLogConvHistory(snes,fnorm,0);
  ierr = SNESMonitor(snes,0,fnorm);CHKERRQ(ierr);

  /* set parameter for default relative tolerance convergence test */
  snes->ttol = fnorm*snes->rtol;
  /* test convergence */
  ierr = (*snes->ops->converged)(snes,0,0.0,0.0,fnorm,&snes->reason,snes->cnvP);CHKERRQ(ierr);
  if (snes->reason) PetscFunctionReturn(0);

  for (i=0; i<maxits; i++) {

    /* Call general purpose update function */
    if (snes->ops->update) {
      ierr = (*snes->ops->update)(snes, snes->iter);CHKERRQ(ierr);
    }

    /* Solve J Y = F, where J is Jacobian matrix */
    ierr = SNESComputeJacobian(snes,x,&snes->jacobian,&snes->jacobian_pre,&flg);CHKERRQ(ierr);
    ierr = KSPSetOperators(snes->ksp,snes->jacobian,snes->jacobian_pre,flg);CHKERRQ(ierr);

    /* PC */
    /*{
      ierr = KSPGetPC(snes->ksp,&pc); CHKERRQ(ierr);
      ierr = PCSetUp(pc); CHKERRQ(ierr);
      void *void_pc_ctx;
      ierr = PCShellGetContext(pc,&void_pc_ctx); CHKERRQ(ierr);
      PCShellCtx *PCCtx = (PCShellCtx*)void_pc_ctx;
      ierr = PCApply(PCCtx->pc,F,SnesCtx->arc_ptr->y_residual); CHKERRQ(ierr);
      ierr = PCApply(PCCtx->pc,SnesCtx->arc_ptr->F_lambda,SnesCtx->arc_ptr->x_lambda); CHKERRQ(ierr);

    }*/

    ierr = KSPSolve(snes->ksp,F,Y);CHKERRQ(ierr);
    ierr = KSPGetConvergedReason(snes->ksp,&kspreason);CHKERRQ(ierr);
    if (kspreason < 0 && ++snes->numLinearSolveFailures >= snes->maxLinearSolveFailures) {
      ierr = PetscInfo2(snes,"iter=%D, number linear solve failures %D greater than current SNES allowed, stopping solve\n",snes->iter,snes->numLinearSolveFailures);CHKERRQ(ierr);
      snes->reason = SNES_DIVERGED_LINEAR_SOLVE;
    } 	else {
      snes->reason = SNES_CONVERGED_ITS;
    }
    ierr = KSPGetIterationNumber(snes->ksp,&lits);CHKERRQ(ierr);
    snes->linear_its += lits;
    ierr = PetscInfo2(snes,"iter=%D, linear solve iterations=%D\n",snes->iter,lits);CHKERRQ(ierr);
    snes->iter++;



    /* Take the computed step. */
    ierr = VecAXPY(x,-1.0,Y);CHKERRQ(ierr);

    ierr = SNESComputeFunction(snes,x,F);CHKERRQ(ierr);
    ierr = SNESGetFunctionDomainError(snes, &domainerror);CHKERRQ(ierr);
    if (domainerror) {
      snes->reason = SNES_DIVERGED_FUNCTION_DOMAIN;
      PetscFunctionReturn(0);
    }

    ierr = VecNormBegin(F,NORM_2,&fnorm);CHKERRQ(ierr);	/* fnorm <- ||F||  */
    ierr = VecNormEnd(F,NORM_2,&fnorm);CHKERRQ(ierr);
    if (PetscIsInfOrNanReal(fnorm)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FP,"User provided compute function generated a Not-a-Number");
    ierr = VecNormBegin(Y,NORM_2,&ynorm);CHKERRQ(ierr);	/* fnorm <- ||F||  */
    ierr = VecNormEnd(Y,NORM_2,&ynorm);CHKERRQ(ierr);
    if (PetscIsInfOrNanReal(ynorm)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FP,"User provided compute function generated a Not-a-Number");

    /* Monitor convergence */
    ierr = PetscObjectTakeAccess(snes);CHKERRQ(ierr);
    snes->iter = i+1;
    snes->norm = fnorm;

    ierr = PetscObjectGrantAccess(snes);CHKERRQ(ierr);
    SNESLogConvHistory(snes,snes->norm,lits);
    ierr = SNESMonitor(snes,snes->iter,snes->norm);CHKERRQ(ierr);
    /* Test for convergence, xnorm = || x || */
    if (snes->ops->converged != SNESSkipConverged) { ierr = VecNorm(x,NORM_2,&xnorm);CHKERRQ(ierr); }
    ierr = (*snes->ops->converged)(snes,snes->iter,xnorm,ynorm,fnorm,&snes->reason,snes->cnvP);CHKERRQ(ierr);
    if (snes->reason) break;
  }
  if (i == maxits) {
    ierr = PetscInfo1(snes,"Maximum number of iterations has been reached: %D\n",maxits);CHKERRQ(ierr);
    if(!snes->reason) snes->reason = SNES_DIVERGED_MAX_IT;
  }

  PetscFunctionReturn(0);
}

/*struct StroreField: public EntMethod {

  string storaga_field_name;
  StroreField(Interface& _moab,string _storaga_field_name): EntMethod(_moab),storaga_field_name(_storaga_field_name) {}

  PetscErrorCode preProcess() {
  }

  PetscErrorCode operator()() {
  }

  PetscErrorCode postProcess() {
  }


};*/

}

#endif //__ARC_LENGHT_NONLINEAR_ELASTICITY_HPP__

