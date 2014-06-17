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
  bool nodal_forces_not_added_and_imperfection;

  ArcComplexForLazyElasticFEMethod(
      FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,double _lambda,double _mu,
      ArcLengthCtx *_arc_ptr,Range &_NodeSet1,int _verbose = 0): 
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_ptr,_verbose), 
      FEMethod_DriverComplexForLazy_Spatial(_mField,_dirihlet_ptr,_lambda,_mu,_verbose), 
      arc_ptr(_arc_ptr),NodeSet1(_NodeSet1),nodal_forces_not_added_and_imperfection(true) {

    set_PhysicalEquationNumber(neohookean);

  }

  PetscScalar imperfection_factor;

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscBool flg = PETSC_TRUE;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-my_if",&imperfection_factor,&flg); CHKERRQ(ierr);
    if(flg!=PETSC_TRUE) {
      imperfection_factor = 1;
    }
    ierr = FEMethod_DriverComplexForLazy_Spatial::preProcess(); CHKERRQ(ierr);
    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	nodal_forces_not_added_and_imperfection = true;
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
	  ierr = CalculateSpatialKFext(PETSC_NULL,arc_ptr->F_lambda,1.); CHKERRQ(ierr);
	  if(nodal_forces_not_added_and_imperfection) {
	    nodal_forces_not_added = true;
	  }
	  ierr = CalculateSpatialKFext(PETSC_NULL,snes_f,imperfection_factor,"Imperfection"); CHKERRQ(ierr);
	  nodal_forces_not_added_and_imperfection = false;
	}
	break;
      case ctx_SNESSetJacobian: {
	  ierr = CalculateSpatialTangent(*snes_B); CHKERRQ(ierr);
	  ierr = CalculateSpatialKFext(*snes_B,PETSC_NULL,arc_ptr->get_FieldData()); CHKERRQ(ierr);
	  ierr = CalculateSpatialKFext(*snes_B,PETSC_NULL,imperfection_factor,"Imperfection"); CHKERRQ(ierr);
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




}

#endif //__ARC_LENGTH_NONLINEAR_ELASTICITY_HPP__

