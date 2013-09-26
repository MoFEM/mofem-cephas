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
#include "PostProcVertexMethod.hpp"
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

#include "ArcLeghtTools.hpp"

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

struct MyElasticFEMethod: public FEMethod_DriverComplexForLazy_Spatial {

  ArcLenghtCtx* arc_ptr;

  Range& NodeSet1;

  MyElasticFEMethod(
      moabField& _mField,BaseDirihletBC *_dirihlet_ptr,double _lambda,double _mu,
      ArcLenghtCtx *_arc_ptr,Range &_NodeSet1,int _verbose = 0): 
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_ptr,_verbose), 
      FEMethod_DriverComplexForLazy_Spatial(_mField,_dirihlet_ptr,_lambda,_mu,_verbose), 
      arc_ptr(_arc_ptr),NodeSet1(_NodeSet1) {

    set_PhysicalEquationNumber(neohookean);

  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy_Spatial::preProcess(); CHKERRQ(ierr);
    switch(snes_ctx) {
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
    ierr = GetIndicesSpatial(); CHKERRQ(ierr);

    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(this,RowGlobSpatial,ColGlobSpatial,DirihletBC); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	  ierr = CalculateSpatialFint(snes_f); CHKERRQ(ierr);
	}
	break;
      case ctx_SNESSetJacobian: {
	  ierr = CalculateSpatialTangent(*snes_B); CHKERRQ(ierr);
	}
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|PressureSet,it)) {

      pressure_cubit_bc_data mydata;
      ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
      /*ostringstream ss;
      ss << *it << endl;
      ss << mydata;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());*/

      double t_val = *(this->t_val)*mydata.data.value1;
      double t[] = { 0,0,t_val, 0,0,t_val, 0,0,t_val };

      Range NeumannSideSet;
      ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),2,NeumannSideSet,true); CHKERRQ(ierr);

      switch(snes_ctx) {
	case ctx_SNESNone:
	case ctx_SNESSetFunction: { 
	  ierr = CaluclateSpatialFext(arc_ptr->F_lambda,t,NeumannSideSet); CHKERRQ(ierr);
	  }
	  break;
	case ctx_SNESSetJacobian: {
	  FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
	  dit = row_multiIndex->get<FieldName_mi_tag>().lower_bound("LAMBDA");
	  hi_dit = row_multiIndex->get<FieldName_mi_tag>().upper_bound("LAMBDA");
	  if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  double _lambda_ = dit->get_FieldData();
	  cblas_dscal(9,_lambda_,t,1);
	  ierr = CalculateSpatialTangentExt(*snes_B,t,NeumannSideSet); CHKERRQ(ierr);
	  } break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }

    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case ctx_SNESNone:
	//snes_f //assemble only if ctx_SNESNone
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,snes_f); CHKERRQ(ierr);
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
	ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,*snes_B); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
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
  Interface& moab;

  ArcLenghtCtx* arc_ptr;
  Vec GhostDiag;
  ArcLenghtElemFEMethod(Interface& _moab,ArcLenghtCtx *_arc_ptr): FEMethod(),moab(_moab),arc_ptr(_arc_ptr) {
    PetscInt ghosts[1] = { 0 };
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm->rank() == 0) {
      VecCreateGhost(PETSC_COMM_WORLD,1,1,0,ghosts,&GhostDiag);
    } else {
      VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&GhostDiag);
    }
  }
  ~ArcLenghtElemFEMethod() {
    VecDestroy(&GhostDiag);
  }


  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
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

    switch(snes_ctx) {
      case ctx_SNESSetFunction: {
	arc_ptr->res_lambda = calulate_lambda_int() - pow(arc_ptr->s,2);
	ierr = VecSetValue(snes_f,dit->get_petsc_gloabl_dof_idx(),arc_ptr->res_lambda,ADD_VALUES); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_SELF,"\tres_lambda = %6.4e\n",arc_ptr->res_lambda);
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
    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	/*//assemble
	ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);*/
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

}

#endif //__ARC_LENGHT_NONLINEAR_ELASTICITY_HPP__

