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

#ifndef __ARC_LENGTH_NONLINEAR_ELASTICITY_HPP__
#define __ARC_LENGTH_NONLINEAR_ELASTICITY_HPP__

#include "FEMethod_DriverComplexForLazy.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include "ArcLengthTools.hpp"

namespace MoFEM {

ErrorCode rval;
PetscErrorCode ierr;

struct ArcComplexForLazyElasticFEMethod: public FEMethod_DriverComplexForLazy_Spatial {

  ArcLengthCtx* arc_ptr;

  Range& NodeSet1;

  ArcComplexForLazyElasticFEMethod(
      FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,double _lambda,double _mu,
      ArcLengthCtx *_arc_ptr,Range &_NodeSet1,int _verbose = 0): 
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
    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(
      this,RowGlobSpatial,ColGlobSpatial,DirihletBC); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	  ierr = CalculateSpatialFint(snes_f); CHKERRQ(ierr);
	  ierr = CalculateSpatialKFext(PETSC_NULL,arc_ptr->F_lambda,1); CHKERRQ(ierr);
	}
	break;
      case ctx_SNESSetJacobian: {
	  ierr = CalculateSpatialTangent(*snes_B); CHKERRQ(ierr);
	  ierr = CalculateSpatialKFext(*snes_B,PETSC_NULL,arc_ptr->get_FieldData()); CHKERRQ(ierr);
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
      case ctx_SNESNone: {
	//snes_f //assemble only if ctx_SNESNone
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,snes_f); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      case ctx_SNESSetFunction: { 
	//F_lambda
	ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	//F_lambda2
	ierr = VecDot(arc_ptr->F_lambda,arc_ptr->F_lambda,&arc_ptr->F_lambda2); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tFlambda^2 = %6.4e\n",arc_ptr->F_lambda2);
	//add F_lambda
	ierr = VecAXPY(snes_f,-arc_ptr->get_FieldData(),arc_ptr->F_lambda); CHKERRQ(ierr);
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
    //NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator lit;
    //lit = numered_dofs_rows.get<FieldName_mi_tag>().find("LAMBDA");
    //if(lit == numered_dofs_rows.get<FieldName_mi_tag>().end()) PetscFunctionReturn(0);
    Range::iterator nit = NodeSet1.begin();
    for(;nit!=NodeSet1.end();nit++) {
      NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator dit,hi_dit;
      dit = numered_dofs_rows.get<MoABEnt_mi_tag>().lower_bound(*nit);
      hi_dit = numered_dofs_rows.get<MoABEnt_mi_tag>().upper_bound(*nit);
      for(;dit!=hi_dit;dit++) {
	PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e -> ","LAMBDA",0,arc_ptr->get_FieldData());
	PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e\n",dit->get_name().c_str(),dit->get_dof_rank(),dit->get_FieldData());
      }
    }
    PetscFunctionReturn(0);
  }

};


struct ArcElasticThermalFEMethod: public ArcComplexForLazyElasticFEMethod {


  double *t_thermal_load_factor_val;
  //set load factor
  PetscErrorCode set_thermal_load_factor(double t_thermal_load_factor_val_) {
      PetscFunctionBegin;
      *t_thermal_load_factor_val = t_thermal_load_factor_val_;
      PetscFunctionReturn(0);
  }

  Tag th_thermal_load_factor;
  ArcElasticThermalFEMethod(
      FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,double _lambda,double _mu,
      ArcLengthCtx *_arc_ptr,Range &_NodeSet1,int _verbose = 0): 
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_ptr,_verbose), 
      ArcComplexForLazyElasticFEMethod(_mField,_dirihlet_ptr,_lambda,_mu,_arc_ptr,_NodeSet1,_verbose) {

    set_PhysicalEquationNumber(neohookean);
    //set_PhysicalEquationNumber(hooke);
    set_ThermalDeformationEquationNumber(linear_expansion_true_volume);
    //set_ThermalDeformationEquationNumber(linear_expanison);

    double def_t_val = 0;
    const EntityHandle root_meshset = mField.get_moab().get_root_set();
    rval = mField.get_moab().tag_get_handle("_ThermalExpansionFactor_t_alpha_val",1,MB_TYPE_DOUBLE,th_thermal_load_factor,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&def_t_val); 
    if(rval == MB_ALREADY_ALLOCATED) {
      rval = mField.get_moab().tag_get_by_ptr(th_thermal_load_factor,&root_meshset,1,(const void**)&t_thermal_load_factor_val); CHKERR_THROW(rval);
    } else {
      CHKERR_THROW(rval);
      rval = mField.get_moab().tag_set_data(th_thermal_load_factor,&root_meshset,1,&def_t_val); CHKERR_THROW(rval);
      rval = mField.get_moab().tag_get_by_ptr(th_thermal_load_factor,&root_meshset,1,(const void**)&t_thermal_load_factor_val); CHKERR_THROW(rval);
    }

  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    ierr = set_thermal_load_factor(arc_ptr->get_FieldData()); CHKERRQ(ierr);
    thermal_load_factor = *t_thermal_load_factor_val;
    ierr = FEMethod_DriverComplexForLazy_Spatial::preProcess(); CHKERRQ(ierr);
    switch(snes_ctx) {
      case ctx_SNESNone: 
      case ctx_SNESSetFunction: {
	//F_lambda
	ierr = VecZeroEntries(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
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
    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(
      this,RowGlobSpatial,ColGlobSpatial,DirihletBC); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	  ierr = CalculateSpatialFint(snes_f); CHKERRQ(ierr);
	  analysis _type_of_analysis = type_of_analysis;
	  type_of_analysis = scaled_themp_direvative_spatial;
	  ierr = GetFint(); CHKERRQ(ierr);
	  iFint_h /= -eps;
	  ierr = VecSetValues(arc_ptr->F_lambda,RowGlobSpatial[i_nodes].size(),&(RowGlobSpatial[i_nodes][0]),&(iFint_h.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	  for(int ee = 0;ee<6;ee++) {
	    if(RowGlobSpatial[1+ee].size()>0) {
	      iFint_h_edge_data[ee] /= -eps;
	      ierr = VecSetValues(arc_ptr->F_lambda,RowGlobSpatial[1+ee].size(),&(RowGlobSpatial[1+ee][0]),&(iFint_h_edge_data[ee].data()[0]),ADD_VALUES); CHKERRQ(ierr);
	    }
	  }
	  for(int ff = 0;ff<4;ff++) {
	    if(RowGlobSpatial[1+6+ff].size()>0) {
	      iFint_h_face_data[ff] /= -eps;
	      ierr = VecSetValues(arc_ptr->F_lambda,RowGlobSpatial[1+6+ff].size(),&(RowGlobSpatial[1+6+ff][0]),&(iFint_h_face_data[ff].data()[0]),ADD_VALUES); CHKERRQ(ierr);
	    }
	  }
	  if(RowGlobSpatial[i_volume].size()>0) {
	    iFint_h_volume /= -eps;
	    ierr = VecSetValues(arc_ptr->F_lambda,RowGlobSpatial[i_volume].size(),&(RowGlobSpatial[i_volume][0]),&(iFint_h_volume.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	  }
	  type_of_analysis = _type_of_analysis;
	}
	break;
      case ctx_SNESSetJacobian: {
	  ierr = CalculateSpatialTangent(*snes_B); CHKERRQ(ierr);
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
      case ctx_SNESNone: {
	//snes_f //assemble only if ctx_SNESNone
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,snes_f); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      case ctx_SNESSetFunction: {
	//F_lambda2
	ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
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

};

struct ArcLengthElemFEMethod: public FieldInterface::FEMethod {
  Interface& moab;

  ArcLengthCtx* arc_ptr;
  Vec GhostDiag;
  ArcLengthElemFEMethod(Interface& _moab,ArcLengthCtx *_arc_ptr): FEMethod(),moab(_moab),arc_ptr(_arc_ptr) {
    PetscInt ghosts[1] = { 0 };
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm->rank() == 0) {
      VecCreateGhost(PETSC_COMM_WORLD,1,1,0,ghosts,&GhostDiag);
    } else {
      VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&GhostDiag);
    }
  }
  ~ArcLengthElemFEMethod() {
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

    switch(snes_ctx) {
      case ctx_SNESSetFunction: {
	arc_ptr->res_lambda = calulate_lambda_int() - pow(arc_ptr->s,2);
	ierr = VecSetValue(snes_f,arc_ptr->get_petsc_gloabl_dof_idx(),arc_ptr->res_lambda,ADD_VALUES); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_SELF,"\tres_lambda = %6.4e\n",arc_ptr->res_lambda);
      }
      break; 
      case ctx_SNESSetJacobian: {
	double diag = 2*arc_ptr->dlambda*pow(arc_ptr->beta,2)*arc_ptr->F_lambda2;
	ierr = VecSetValue(GhostDiag,0,diag,INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(*snes_B,arc_ptr->get_petsc_gloabl_dof_idx(),arc_ptr->get_petsc_gloabl_dof_idx(),1,ADD_VALUES); CHKERRQ(ierr);
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
	PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",arc_ptr->get_FieldData());  
	//snes_f norm
	double fnorm;
	ierr = VecNorm(snes_f,NORM_2,&fnorm); CHKERRQ(ierr);	
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
    if(arc_ptr->get_petsc_local_dof_idx()!=-1) {
      double *array;
      ierr = VecGetArray(arc_ptr->dx,&array); CHKERRQ(ierr);
      arc_ptr->dlambda = array[arc_ptr->get_petsc_local_dof_idx()];
      array[arc_ptr->get_petsc_local_dof_idx()] = 0;
      ierr = VecRestoreArray(arc_ptr->dx,&array); CHKERRQ(ierr);
    }
    int part = arc_ptr->get_part();
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
    //check if locl dof idx is non zero, i.e. that lambda is acessible from this processor
    if(arc_ptr->get_petsc_local_dof_idx()!=-1) {
      double *array;
      ierr = VecGetArray(x,&array); CHKERRQ(ierr);
      double lambda_old = array[arc_ptr->get_petsc_local_dof_idx()];
      if(!(dlambda == dlambda)) {
	ostringstream sss;
	sss << "s " << arc_ptr->s << " " << arc_ptr->beta << " " << arc_ptr->F_lambda2;
	SETERRQ(PETSC_COMM_SELF,1,sss.str().c_str());
      }
      array[arc_ptr->get_petsc_local_dof_idx()] = lambda_old + dlambda;
      PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e, %6.4e (%6.4e)\n",
	lambda_old, array[arc_ptr->get_petsc_local_dof_idx()], dlambda);
      ierr = VecRestoreArray(x,&array); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

};

}

#endif //__ARC_LENGTH_NONLINEAR_ELASTICITY_HPP__

