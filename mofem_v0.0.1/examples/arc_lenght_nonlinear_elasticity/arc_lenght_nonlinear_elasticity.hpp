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

namespace MoFEM {

ErrorCode rval;
PetscErrorCode ierr;

struct MatShellCtx {
  moabField& mField;

  Mat Aij;
  Vec F_lambda,db;
  double diag_lambda;
  MatShellCtx(moabField& _mField,Mat _Aij,Vec _F_lambda,Vec _db): mField(_mField),Aij(_Aij),F_lambda(_F_lambda),db(_db) {};
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
  MatShellGetContext(A,&void_ctx);
  MatShellCtx *ctx = (MatShellCtx*)void_ctx;
  ierr = MatMult(ctx->Aij,x,f); CHKERRQ(ierr);
  double db_dot_x;
  ierr = VecDot(ctx->db,x,&db_dot_x); CHKERRQ(ierr);
  double f_lambda;
  ierr = ctx->set_lambda(f,&f_lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  f_lambda += db_dot_x;
  ierr = ctx->set_lambda(f,&f_lambda,SCATTER_REVERSE); CHKERRQ(ierr);
  double lambda;
  ierr = ctx->set_lambda(x,&lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  //PetscPrintf(PETSC_COMM_WORLD,"mat_mult lambda = %6.4e\n",lambda);
  ierr = VecAXPY(f,lambda,ctx->F_lambda); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


struct MyElasticFEMethod: public FEMethod_DriverComplexForLazy {

  Vec F_lambda,b,db;
  Range& SideSet1;
  Range& SideSet2;
  Range& SideSet3;
  Range& SideSet4;
  Range& NodeSet1;

  Range SideSet1_;
  Range SideSet3_;
  Range SideSet4_;

  MyElasticFEMethod(Interface& _moab,double _lambda,double _mu,
      Vec _F_lambda,Vec _b,Vec _db,
      Range &_SideSet1,Range &_SideSet2,
      Range &_SideSet3,Range &_SideSet4,
      Range &_NodeSet1,
      int _verbose = 0): 
      FEMethod_DriverComplexForLazy(_moab,_lambda,_mu,_verbose), 
      F_lambda(_F_lambda),b(_b),db(_db),
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
      case ctx_SNESSetFunction: { 
      	ierr = VecZeroEntries(F_lambda); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      	ierr = VecZeroEntries(b); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(b,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(b,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
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

    switch(ctx) {
      case ctx_SNESSetFunction: { 
	  ierr = CalculateFint(snes_f); CHKERRQ(ierr);
	  double t[] = { 0,0,t_val, 0,0,t_val, 0,0,t_val };
	  ierr = CaluclateFext(F_lambda,t,NeumannSideSet); CHKERRQ(ierr);
	}
	break;
      case ctx_SNESSetJacobian: {
	  ierr = CalculateTangent(*snes_B); CHKERRQ(ierr);
	  FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
	  dit = row_multiIndex->get<FieldName_mi_tag>().lower_bound("LAMBDA");
	  hi_dit = row_multiIndex->get<FieldName_mi_tag>().upper_bound("LAMBDA");
	  //only one LAMBDA
	  if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  double _lambda_ = dit->get_FieldData();
	  //PetscPrintf(PETSC_COMM_WORLD,"snes _lambda_ = %6.4e\n",_lambda_);  
	  double t[] = { 0,0,t_val, 0,0,t_val, 0,0,t_val };
	  cblas_dscal(6,_lambda_,t,1);
	  ierr = CalculateTangentExt(*snes_B,t,NeumannSideSet); CHKERRQ(ierr);
	}
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    switch(ctx) {
      case ctx_SNESSetFunction: { 
  
	FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
	dit = row_multiIndex->get<FieldName_mi_tag>().lower_bound("SPATIAL_POSITION");
	hi_dit = row_multiIndex->get<FieldName_mi_tag>().upper_bound("SPATIAL_POSITION");
	for(;dit!=hi_dit;dit++) {
	  //if(dit->get_ent_type()!=MBVERTEX) continue;
	  //if(find(NodeSet1.begin(),NodeSet1.end(),dit->get_ent())==NodeSet1.end()) continue;
	  //(x0+dx)*(x0+dx) + (lambd+dlambda)^2*beta^2*F_lambda*F_Lambda - s = 0
	  //x0*x0 + 2x0*dx + dx*dx + (lambda^2+2*lambda*dlmabda+dlambda^2)*(....) - s = 0
	  ierr = VecSetValue(b,dit->get_petsc_gloabl_dof_idx(),dit->get_FieldData(),INSERT_VALUES); CHKERRQ(ierr);
	  ierr = VecSetValue(db,dit->get_petsc_gloabl_dof_idx(),2.*dit->get_FieldData(),INSERT_VALUES); CHKERRQ(ierr);
	}
      }
      break;
      default:
      break;
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

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecGhostUpdateBegin(F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(F_lambda); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F_lambda); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(b); CHKERRQ(ierr);
	//ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(db); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(db); CHKERRQ(ierr);
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

};

struct ArcLenghtElemFEMethod: public moabField::FEMethod {

  double s;
  PetscErrorCode set_s(double _s) { 
    PetscFunctionBegin;
    s = _s;
    PetscFunctionReturn(0);
  }

  Vec GhostLambda,GhostDiag;
  Vec F_lambda,b;
  double beta;
  ArcLenghtElemFEMethod(Interface& _moab,Vec _F_lambda,Vec _b): FEMethod(_moab),F_lambda(_F_lambda),b(_b),beta(1e-3) {
    PetscInt ghosts[1] = { 0 };
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm->rank() == 0) {
      VecCreateGhost(PETSC_COMM_WORLD,1,1,1,ghosts,&GhostLambda);
      VecCreateGhost(PETSC_COMM_WORLD,1,1,1,ghosts,&GhostDiag);
    } else {
      VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&GhostLambda);
      VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&GhostDiag);
    }
  }

  PetscInt iter;
  double s0;
  double b_dot_x;
  double F_lambda2;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    switch(ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecDot(snes_x,b,&b_dot_x); CHKERRQ(ierr);
	ierr = SNESGetIterationNumber(snes,&iter); CHKERRQ(ierr);
	if(iter == 0) {
	  ierr = VecDot(snes_x,b,&s0); CHKERRQ(ierr);
	  ierr = VecDot(F_lambda,F_lambda,&F_lambda2); CHKERRQ(ierr);
	  DofMoFEMEntity_multiIndex& dofs_moabfield_no_const = const_cast<DofMoFEMEntity_multiIndex&>(*dofs_moabfield);
	  DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
	  dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().lower_bound("LAMBDA");
	  hi_dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().upper_bound("LAMBDA");
	  if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  double lambda = dit->get_FieldData();
	  s0 += lambda*lambda*beta*F_lambda2;
	}
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

    FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
    dit = row_multiIndex->get<FieldName_mi_tag>().lower_bound("LAMBDA");
    hi_dit = row_multiIndex->get<FieldName_mi_tag>().upper_bound("LAMBDA");
    //only one LAMBDA
    if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");

    switch(ctx) {
      case ctx_SNESSetFunction: {
	double res_lambda,lambda;
	PetscScalar *array;
	//
	ierr = VecGetArray(snes_x,&array); CHKERRQ(ierr);
	lambda = array[dit->get_petsc_local_dof_idx()];
	ierr = VecRestoreArray(snes_x,&array); CHKERRQ(ierr);
	ierr = VecSetValue(GhostLambda,0,lambda,INSERT_VALUES); CHKERRQ(ierr);
	//
	res_lambda = b_dot_x + lambda*lambda*beta*F_lambda2  - (s0 + s);
	ierr = VecSetValue(snes_f,dit->get_petsc_gloabl_dof_idx(),res_lambda,ADD_VALUES); CHKERRQ(ierr);
      }
      break; 
      case ctx_SNESSetJacobian: {
	double lambda;
	PetscScalar *array;
	//
	ierr = VecGetArray(snes_x,&array); CHKERRQ(ierr);
	lambda = array[dit->get_petsc_local_dof_idx()];
	ierr = VecRestoreArray(snes_x,&array); CHKERRQ(ierr);
	const double eps = 1e-12;
	double diag = 2*lambda*beta*F_lambda2;
	if(fabs(diag) < eps) diag = eps;
	ierr = VecSetValue(GhostDiag,0,diag,INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(*snes_B,dit->get_petsc_gloabl_dof_idx(),dit->get_petsc_gloabl_dof_idx(),diag,ADD_VALUES); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"no implemented");
    }
    
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(GhostLambda); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(GhostLambda); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(GhostLambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(GhostLambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	double *lambda;
	ierr = VecGetArray(GhostLambda,&lambda); CHKERRQ(ierr);
	ierr = VecAXPY(snes_f,*lambda,F_lambda); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",*lambda);  
	ierr = VecRestoreArray(GhostLambda,&lambda); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
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
	void *void_ctx;
	MatShellGetContext(*snes_A,&void_ctx);
	MatShellCtx *MatCtx = (MatShellCtx*)void_ctx;
	MatCtx->diag_lambda = *diag;
	ierr = VecRestoreArray(GhostDiag,&diag); CHKERRQ(ierr);
	//Matrix View
	//MatView(*snes_B,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
	//std::string wait;
	//std::cin >> wait;
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

};

struct PCShellCtx {
  PC pc;
  Mat ShellAij,Aij;
  Vec x_lambda;
  PCShellCtx(Mat _ShellAij,Mat _Aij,Vec _x_lambda): 
    ShellAij(_ShellAij),Aij(_Aij),x_lambda(_x_lambda) {
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
  double res_lambda;
  ierr = MatCtx->set_lambda(pc_f,&res_lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = PCApply(PCCtx->pc,pc_f,pc_x); CHKERRQ(ierr);
  ierr = PCApply(PCCtx->pc,MatCtx->F_lambda,PCCtx->x_lambda); CHKERRQ(ierr);
  // db \dot pc_x + dlambda*diag_lambda - res_lambda = 0
  // pc_x = x_int + dlambda*x_lambda
  // db \dot x_int + dlambda*(b \dot x_lambda + diag_lambda) - res_lambda = 0
  // dlambda = (res_lambda - db \dot x_int)/(b \cdot x_lambda + diag_lambda)
  double db_dot_x_int,db_dot_x_lambda;
  ierr = VecDot(MatCtx->db,pc_x,&db_dot_x_int); CHKERRQ(ierr);
  ierr = VecDot(MatCtx->db,PCCtx->x_lambda,&db_dot_x_lambda); CHKERRQ(ierr);
  double dlambda;
  dlambda = (res_lambda - db_dot_x_int)/(db_dot_x_lambda+MatCtx->diag_lambda);
  if(dlambda != dlambda) SETERRQ(PETSC_COMM_SELF,1,"db \\dot x_lambda = 0\nCheck constrint vector, SideSet, ect.");
  //PetscPrintf(PETSC_COMM_WORLD,"pc dlambda = %6.4e\n",dlambda);
  ierr = VecAXPY(pc_x,dlambda,PCCtx->x_lambda); CHKERRQ(ierr);
  ierr = MatCtx->set_lambda(pc_x,&dlambda,SCATTER_REVERSE); CHKERRQ(ierr);
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

}

#endif //__ARC_LENGHT_NONLINEAR_ELASTICITY_HPP__

