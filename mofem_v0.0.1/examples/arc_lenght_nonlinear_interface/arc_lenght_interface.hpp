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

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "ElasticFEMethod.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include "ElasticFEMethodForInterface.hpp"
#include "ArcLeghtTools.hpp"

#ifdef __cplusplus
extern "C" {
#endif
#include <petsc-private/snesimpl.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

struct ArcInterfaceElasticFEMethod: public InterfaceElasticFEMethod {

  ArcInterfaceElasticFEMethod(Interface& _moab): InterfaceElasticFEMethod(_moab) {};

  ArcLenghtCtx *arc_ptr;
  ArcInterfaceElasticFEMethod(
      Interface& _moab,Mat &_Aij,Vec& _F,
      double _lambda,double _mu,Range &_SideSet1,Range &_SideSet2,Range &_SideSet3,
      ArcLenghtCtx *_arc_ptr): 
      InterfaceElasticFEMethod(_moab,_Aij,_F,_lambda,_mu,_SideSet1,_SideSet2,_SideSet3)
      ,arc_ptr(_arc_ptr) {};

  PetscErrorCode NeumannBC() {
      PetscFunctionBegin;
      
      ublas::vector<FieldData,ublas::bounded_array<double,3> > traction2(3);
      traction2[0] = 0;
      traction2[1] = +1;
      traction2[2] = 0;
      ierr = ElasticFEMethod::NeumannBC(arc_ptr->F_lambda,traction2,SideSet2); CHKERRQ(ierr);

      ublas::vector<FieldData,ublas::bounded_array<double,3> > traction3(3);
      traction3[0] = 0;
      traction3[1] = -1;
      traction3[2] = 0;
      ierr = ElasticFEMethod::NeumannBC(arc_ptr->F_lambda,traction3,SideSet3); CHKERRQ(ierr);

      PetscFunctionReturn(0);
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    g_NTET.resize(4*45);
    ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
    g_NTRI.resize(3*13);
    ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
    // See FEAP - - A Finite Element Analysis Program
    D_lambda = ublas::zero_matrix<FieldData>(6,6);
    for(int rr = 0;rr<3;rr++) {
	ublas::matrix_row<ublas::matrix<FieldData> > row_D_lambda(D_lambda,rr);
	for(int cc = 0;cc<3;cc++) {
	  row_D_lambda[cc] = 1;
	}
    }
    D_mu = ublas::zero_matrix<FieldData>(6,6);
    for(int rr = 0;rr<6;rr++) {
	D_mu(rr,rr) = rr<3 ? 2 : 1;
    }
    D = lambda*D_lambda + mu*D_mu;


    switch(ctx) {
      case ctx_SNESNone: {}
      break;
      case ctx_SNESSetFunction: { 
	//F_lambda
	ierr = VecZeroEntries(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	Diagonal = PETSC_NULL;
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = VecDuplicate(F,&Diagonal); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }


  PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      ierr = GetMatrices(); CHKERRQ(ierr);

      //Dirihlet Boundary Condition

      switch(ctx) {
	case ctx_SNESNone: {
	}
	break;
	case ctx_SNESSetJacobian: 
	case ctx_SNESSetFunction: { 
	  ApplyDirihletBC();
	  if(Diagonal!=PETSC_NULL) {
	    if(DirihletBC.size()>0) {
	      DirihletBCDiagVal.resize(DirihletBC.size());
	      fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
	      ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
	    }
	  }
	}
	break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }

      switch(ctx) {
	case ctx_SNESNone: 
	case ctx_SNESSetFunction: {
	  //Assembly  F
	  ierr = Fint(); CHKERRQ(ierr);
	  //Neumann Boundary Conditions
	  ierr = NeumannBC(); CHKERRQ(ierr);
	}
	break;
	case ctx_SNESSetJacobian: {
	  //Assembly  F
	  ierr = Lhs(); CHKERRQ(ierr);
	}
	break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
  }


  PetscErrorCode postProcess() {
    PetscFunctionBegin;

    switch(ctx) {
      case ctx_SNESNone: {
      }
      case ctx_SNESSetFunction: { 
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
      case ctx_SNESSetJacobian: {
	ierr = VecAssemblyBegin(Diagonal); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(Diagonal); CHKERRQ(ierr);
	ierr = MatDiagonalSet(Aij,Diagonal,ADD_VALUES); CHKERRQ(ierr);
	ierr = VecDestroy(&Diagonal); CHKERRQ(ierr);
	//Note MAT_FLUSH_ASSEMBLY
	ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

};

struct ArcInterfaceFEMethod: public InterfaceFEMethod {
  ArcInterfaceFEMethod(
      Interface& _moab,double _YoungModulus): 
      InterfaceFEMethod(_moab,_YoungModulus) {};

  Vec D;
  ArcInterfaceFEMethod(
      Interface& _moab,Mat &_Aij,Vec& _F,Vec& _D,double _YoungModulus,Range &_SideSet1,Range &_SideSet2,Range &_SideSet3): 
      InterfaceFEMethod(_moab,_Aij,_F,_YoungModulus,_SideSet1,_SideSet2,_SideSet3),D(_D) {};

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    ierr = PetscGetTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
    g_NTET.resize(4*45);
    ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
    g_NTRI.resize(3*13);
    ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 

    switch(ctx) {
      case ctx_SNESNone: {}
      break;
      case ctx_SNESSetFunction: { }
      break;
      case ctx_SNESSetJacobian: { }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }


    PetscFunctionReturn(0);
  }


  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = OpStudentStart_PRISM(g_NTRI); CHKERRQ(ierr);

    //Rotation matrix
    ierr = CalcR(); CHKERRQ(ierr);
    //Dglob
    ierr = CalcDglob(); CHKERRQ(ierr);
    //Calculate Matrices
    ierr = Matrices();    CHKERRQ(ierr);

    switch(ctx) {
      case ctx_SNESNone: {
      }
      break;
      case ctx_SNESSetJacobian: 
      case ctx_SNESSetFunction: { 
	//Apply Dirihlet BC
	ApplyDirihletBC();
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    switch(ctx) {
      case ctx_SNESNone: {
	ierr = RhsInt(); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetFunction: { 
	ierr = RhsInt(); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: { 
	ierr = LhsInt(); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }


    ierr = OpStudentEnd(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


  PetscErrorCode postProcess() {
    PetscFunctionBegin;

    switch(ctx) {
      case ctx_SNESNone: {
      }
      break;
      case ctx_SNESSetFunction: {}
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }


    PetscFunctionReturn(0);
  }

};

struct ArcLenghtIntElemFEMethod: public moabField::FEMethod {
  Interface& moab;
  PetscErrorCode ierr;

  Mat Aij;
  Vec F,D;
  ArcLenghtCtx* arc_ptr;
  Vec GhostDiag;
  ArcLenghtIntElemFEMethod(Interface& _moab,Mat &_Aij,Vec& _F,Vec& _D,
    ArcLenghtCtx *_arc_ptr): FEMethod(),moab(_moab),Aij(_Aij),F(_F),D(_D),arc_ptr(_arc_ptr) {
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
	ierr = calulate_dx_and_dlambda(D); CHKERRQ(ierr);
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
	ierr = VecSetValue(F,dit->get_petsc_gloabl_dof_idx(),arc_ptr->res_lambda,ADD_VALUES); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tres_lambda = %6.4e\n",arc_ptr->res_lambda);
      }
      break; 
      case ctx_SNESSetJacobian: {
	double diag = 2*arc_ptr->dlambda*pow(arc_ptr->beta,2)*arc_ptr->F_lambda2;
	ierr = VecSetValue(GhostDiag,0,diag,INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(Aij,dit->get_petsc_gloabl_dof_idx(),dit->get_petsc_gloabl_dof_idx(),1,ADD_VALUES); CHKERRQ(ierr);
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
	//add F_lambda
	NumeredDofMoFEMEntity_multiIndex& dofs_moabfield_no_const 
	    = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
	NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
	dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().lower_bound("LAMBDA");
	hi_dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().upper_bound("LAMBDA");
	if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	ierr = VecAXPY(F,-dit->get_FieldData(),arc_ptr->F_lambda); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",dit->get_FieldData());  
	//snes_f norm
	double fnorm;
	ierr = VecNormBegin(F,NORM_2,&fnorm); CHKERRQ(ierr);	
	ierr = VecNormEnd(F,NORM_2,&fnorm);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tfnorm = %6.4e\n",fnorm);  
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = VecAssemblyBegin(GhostDiag); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(GhostDiag); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(GhostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(GhostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	double *diag;
	ierr = VecGetArray(GhostDiag,&diag); CHKERRQ(ierr);
	arc_ptr->diag = *diag;
	ierr = VecRestoreArray(GhostDiag,&diag); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tdiag = %6.4e\n",arc_ptr->diag);
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
      PetscPrintf(PETSC_COMM_WORLD,"\tInit dlambda = %6.4e s = %6.4e beta = %6.4e F_lambda2 = %6.4e\n",*dlambda,arc_ptr->s,arc_ptr->beta,arc_ptr->F_lambda2);
      double a = *dlambda;
      if(a - a != 0) {
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
      if(dit == dofs_moabfield_no_const.get<FieldName_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
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




}
