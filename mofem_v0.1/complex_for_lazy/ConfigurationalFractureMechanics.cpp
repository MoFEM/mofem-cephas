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

#include "ConfigurationalFractureMechanics.hpp"
#include "FieldCore.hpp"
#include "FEMethod_ComplexConstArea.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include "petscShellMATs_ConstrainsByMarkAinsworth.hpp"
#include "FaceSplittingTool.hpp"

#include <petscsnes.h>

using namespace MoFEM;

phisical_equation_volume eq_solid = hooke; /*stvenant_kirchhoff;*/

struct NL_ElasticFEMethod: public FEMethod_DriverComplexForLazy_Spatial {
  
    NL_ElasticFEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0): 
        FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
        FEMethod_DriverComplexForLazy_Spatial(_mField,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose)  {
      set_PhysicalEquationNumber(eq_solid);
    }
  
  };
  
struct NL_MaterialFEMethod: public FEMethod_DriverComplexForLazy_Material {
  
  NL_MaterialFEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0): 
        FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
        FEMethod_DriverComplexForLazy_Material(_mField,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose)  {
      set_PhysicalEquationNumber(eq_solid);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,*snes_B); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  
  };
  
struct NL_ElasticFEMethodCoupled: public FEMethod_DriverComplexForLazy_CoupledSpatial {
  
  ArcLengthCtx *arc_ptr;
  NL_ElasticFEMethodCoupled(
    FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,ArcLengthCtx *_arc_ptr,int _verbose = 0):
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
      FEMethod_DriverComplexForLazy_CoupledSpatial(_mField,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose), arc_ptr(_arc_ptr) {
	set_PhysicalEquationNumber(eq_solid);
      }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case ctx_SNESNone: {
	nodal_forces_not_added = true;
	ierr = VecZeroEntries(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
      } break;
      case ctx_SNESSetFunction: { 
        nodal_forces_not_added = true;
	ierr = VecZeroEntries(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
      } break;
      case ctx_SNESSetJacobian: {
        ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      } break;
      default:
        SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndicesSpatial(); CHKERRQ(ierr);
    DirihletBC.resize(0);
    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	ierr = CalculateSpatialKFext(PETSC_NULL,arc_ptr->F_lambda,-1.); CHKERRQ(ierr);
	ierr = CalculateSpatialFint(snes_f); CHKERRQ(ierr);
      } break;
      case ctx_SNESSetJacobian: {
	ierr = CalculateSpatialKFext(*snes_B,PETSC_NULL,-arc_ptr->get_FieldData()); CHKERRQ(ierr);
	ierr = CalculateSpatialTangent(*snes_B); CHKERRQ(ierr);
	ierr = GetIndicesRow(RowGlobMaterial,material_field_name); CHKERRQ(ierr);
	ierr = AssembleSpatialCoupledTangent(*snes_B); CHKERRQ(ierr);
      } break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case ctx_SNESNone: {
      } break;
      case ctx_SNESSetFunction:
      case ctx_SNESSetJacobian: {
	switch(snes_ctx) {
	  case ctx_SNESSetFunction: { 
	    ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	  } break;
	  case ctx_SNESSetJacobian: {
	    ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	    ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  } break;
	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
	}
      } break;
      default:
	break;
    }
    PetscFunctionReturn(0);
  }

};

struct NL_MaterialFEMethodCoupled: public FEMethod_DriverComplexForLazy_CoupledMaterial {
  ArcLengthCtx *arc_ptr;
  NL_MaterialFEMethodCoupled(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,ArcLengthCtx *_arc_ptr,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    FEMethod_DriverComplexForLazy_CoupledMaterial(_mField,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose),arc_ptr(_arc_ptr) {
	set_PhysicalEquationNumber(eq_solid);
    }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case ctx_SNESNone: {
      } break;
      case ctx_SNESSetFunction: { 
        nodal_forces_not_added = true;
        ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
        ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
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
    ierr = GetIndicesMaterial(); CHKERRQ(ierr);
    ierr = GetData(
	order_x_edges,order_x_faces,order_x_volume,
	dofs_x_edge_data,dofs_x_edge,
	dofs_x_face_data,dofs_x_face,
	dofs_x_volume,dofs_x,
	spatial_field_name); CHKERRQ(ierr);
    switch(snes_ctx) {
	case ctx_SNESNone: 
	case ctx_SNESSetFunction: { 
	  ierr = CalculateMaterialFint(snes_f); CHKERRQ(ierr);
	  ierr = CaluclateMaterialKFext(PETSC_NULL,arc_ptr->F_lambda,-1); CHKERRQ(ierr);
	} break;
	case ctx_SNESSetJacobian: {
	  ierr = GetTangent(); CHKERRQ(ierr);
	  ierr = AssembleMaterialTangent(*snes_B); CHKERRQ(ierr);
	  ierr = CaluclateMaterialKFext(*snes_B,PETSC_NULL,-arc_ptr->get_FieldData()); CHKERRQ(ierr);
	  ierr = GetIndicesRow(RowGlobSpatial,spatial_field_name); CHKERRQ(ierr);
	  ierr = AssembleMaterialCoupledTangent(*snes_B); CHKERRQ(ierr);
	} break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionBegin;
    switch(snes_ctx) {
      case ctx_SNESNone: 
      case ctx_SNESSetFunction: {
	  ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	  ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,arc_ptr->F_lambda); CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	  //F_lambda2
	  ierr = VecDot(arc_ptr->F_lambda,arc_ptr->F_lambda,&arc_ptr->F_lambda2); CHKERRQ(ierr);
	  PetscPrintf(PETSC_COMM_WORLD,"\tsquare root of reference load vector FlambdaT*Flambda = %6.4e\n",arc_ptr->F_lambda2);
	  ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	  //Add external forces to RHS
	  ierr = VecAXPY(snes_f,-arc_ptr->get_FieldData(),arc_ptr->F_lambda); CHKERRQ(ierr);
	  PetscPrintf(PETSC_COMM_WORLD,"\tload factor lambda = %6.4e\n",arc_ptr->get_FieldData());  
	  ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	  ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,snes_f); CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	} break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,*snes_B); CHKERRQ(ierr);
      } break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  };

struct NL_MeshSmootherCoupled: public FEMethod_DriverComplexForLazy_MeshSmoothing {
  
  NL_MeshSmootherCoupled(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,int _verbose = 0):
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
      FEMethod_DriverComplexForLazy_MeshSmoothing(_mField,_dirihlet_bc_method_ptr) {
	set_qual_ver(0);
      }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case ctx_SNESNone: {
      } break;
      case ctx_SNESSetFunction: { 
        nodal_forces_not_added = true;
        ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	if(!crackFrontNodes.empty()) {
	  if(front_f == PETSC_NULL) {
	    ierr = VecDuplicate(snes_f,&front_f); CHKERRQ(ierr);
	    ierr = VecSetOption(front_f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); CHKERRQ(ierr);
	    ierr = VecDuplicate(snes_f,&tangent_front_f); CHKERRQ(ierr);
	    ierr = VecSetOption(tangent_front_f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); CHKERRQ(ierr);
	  }
	  ierr = VecZeroEntries(front_f); CHKERRQ(ierr);
	  ierr = VecGhostUpdateBegin(front_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	  ierr = VecGhostUpdateEnd(front_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	}
      } break;
      case ctx_SNESSetJacobian: {
        ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
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
      } break;
      case ctx_SNESSetFunction:
      case ctx_SNESSetJacobian: {
	switch(snes_ctx) {
	  case ctx_SNESSetFunction: { 
	    ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	    if(!crackFrontNodes.empty()) {
	      ierr = VecAssemblyBegin(front_f); CHKERRQ(ierr);
	      ierr = VecAssemblyEnd(front_f); CHKERRQ(ierr);
	      ierr = VecGhostUpdateBegin(front_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	      ierr = VecGhostUpdateEnd(front_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	      ierr = VecGhostUpdateBegin(front_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	      ierr = VecGhostUpdateEnd(front_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    }
	  } break;
	  case ctx_SNESSetJacobian: {
	    ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	    ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  } break;
	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
	}
      } break;
      default:
	break;
    }
    PetscFunctionReturn(0);
  }

  
  };

PetscErrorCode ConfigurationalFractureMechanics::CubitDisplacementDirihletBC_Coupled::SetDirihletBC_to_ElementIndicies(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::CubitDisplacementDirihletBC_Coupled::SetDirihletBC_to_ElementIndiciesRow(
      FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirihletBC) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}
  
PetscErrorCode ConfigurationalFractureMechanics::CubitDisplacementDirihletBC_Coupled::SetDirihletBC_to_ElementIndiciesCol(
  FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::CubitDisplacementDirihletBC_Coupled::SetDirihletBC_to_MatrixDiagonal(FieldInterface::FEMethod *fe_method_ptr,Mat Aij) {
  PetscFunctionBegin;
  set<DofIdx> set_zero_rows;
  if(fixAllSpatialDispacements) {
    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_LOCIDX_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
      if( (dit->get_name() != "SPATIAL_POSITION") ) continue;
      set_zero_rows.insert(dit->get_petsc_gloabl_dof_idx());
    }
  }
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
    if(dit->get_ent_type()!=MBVERTEX) continue;
    if( (dit->get_name() == "SPATIAL_POSITION") ) continue;
    if(find(cornersNodes.begin(),cornersNodes.end(),dit->get_ent()) == cornersNodes.end()) continue;
    set_zero_rows.insert(dit->get_petsc_gloabl_dof_idx());
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|DisplacementSet,it)) {
    for(int dim = 0;dim<3;dim++) {
      Range _ents;
      ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),dim,_ents,true); CHKERRQ(ierr);
      if(dim>1) {
        Range _edges;
        ierr = mField.get_moab().get_adjacencies(_ents,1,false,_edges,Interface::UNION); CHKERRQ(ierr);
        _ents.insert(_edges.begin(),_edges.end());
      }
      if(dim>0) {
        Range _nodes;
        rval = mField.get_moab().get_connectivity(_ents,_nodes,true); CHKERR_PETSC(rval);
        _ents.insert(_nodes.begin(),_nodes.end());
      }
      for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
	if(dit->get_ent_type()!=MBVERTEX) continue;
	if(
	  (dit->get_name() == "SPATIAL_POSITION")
	) continue;
	if(find(_ents.begin(),_ents.end(),dit->get_ent()) == _ents.end()) continue;
	set_zero_rows.insert(dit->get_petsc_gloabl_dof_idx());
      }
    }
  }
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_LOCIDX_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
    if(
      (dit->get_name() != "SPATIAL_POSITION")
      ) continue;
    for(int ss = 0;ss<3;ss++) {
      if(dit->get_dof_rank()==ss) {
	map<int,Range>::iterator bit = bc_map[ss].begin();
	for(;bit!=bc_map[ss].end();bit++) {
	  if(find(bit->second.begin(),bit->second.end(),dit->get_ent()) == bit->second.end()) continue;
	  set_zero_rows.insert(dit->get_petsc_gloabl_dof_idx());
	}
      }
    }
  }
  vector<DofIdx> zero_rows(set_zero_rows.size());
  copy(set_zero_rows.begin(),set_zero_rows.end(),zero_rows.begin());
  ierr = MatZeroRowsColumns(
    Aij,zero_rows.size(),&*zero_rows.begin(),1.,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::CubitDisplacementDirihletBC_Coupled::SetDirihletBC_to_RHS(FieldInterface::FEMethod *fe_method_ptr,Vec F) {
  PetscFunctionBegin;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
  if(fixAllSpatialDispacements) {
    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(fe_method_ptr->problem_ptr,"SPATIAL_POSITION",dit)) {
      ierr = VecSetValue(F,dit->get_petsc_gloabl_dof_idx(),0.,INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
    if(dit->get_part()!=pcomm->rank()) continue;
    if(dit->get_ent_type()!=MBVERTEX) continue;
    if(
      (dit->get_name() == "SPATIAL_POSITION")
      ) continue;
    if(find(cornersNodes.begin(),cornersNodes.end(),dit->get_ent()) == cornersNodes.end()) continue;
    ierr = VecSetValue(F,dit->get_petsc_gloabl_dof_idx(),0.,INSERT_VALUES); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|DisplacementSet,it)) {
    for(int dim = 0;dim<3;dim++) {
      Range _ents;
      ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),dim,_ents,true); CHKERRQ(ierr);
      if(dim>1) {
        Range _edges;
        ierr = mField.get_moab().get_adjacencies(_ents,1,false,_edges,Interface::UNION); CHKERRQ(ierr);
        _ents.insert(_edges.begin(),_edges.end());
      }
      if(dim>0) {
        Range _nodes;
        rval = mField.get_moab().get_connectivity(_ents,_nodes,true); CHKERR_PETSC(rval);
        _ents.insert(_nodes.begin(),_nodes.end());
      }
      for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
	if(dit->get_ent_type()!=MBVERTEX) continue;
	if(
	  (dit->get_name() == "SPATIAL_POSITION")
	) continue;
	if(find(_ents.begin(),_ents.end(),dit->get_ent()) == _ents.end()) continue;
	ierr = VecSetValue(F,dit->get_petsc_gloabl_dof_idx(),0.,INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_LOCIDX_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
    if(dit->get_part()!=pcomm->rank()) continue;
    if(
      (dit->get_name() != "SPATIAL_POSITION")
      ) continue;
    for(int ss = 0;ss<3;ss++) {
      if(dit->get_dof_rank()==ss) {
	map<int,Range>::iterator bit = bc_map[ss].begin();
	for(;bit!=bc_map[ss].end();bit++) {
	  if(find(bit->second.begin(),bit->second.end(),dit->get_ent()) == bit->second.end()) continue;
	  ierr = VecSetValue(F,dit->get_petsc_gloabl_dof_idx(),0.,INSERT_VALUES); CHKERRQ(ierr);
	}
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::CubitDisplacementDirihletBC_Coupled::SetDirihletBC_to_ElementIndiciesFace(FieldInterface::FEMethod *fe_method_ptr,
  vector<DofIdx>& DirihletBC,vector<DofIdx>& FaceNodeGlobalDofs,vector<vector<DofIdx> > &FaceEdgeGlobalDofs,vector<DofIdx> &FaceGlobalDofs) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::set_material_fire_wall(FieldInterface& mField) {
  PetscFunctionBegin;
  ErrorCode rval;
  BitRefLevel def_bit_level = 0;
  rval = mField.get_moab().tag_get_handle("_Materiar_FireWall",sizeof(Material_FirelWall_def),MB_TYPE_OPAQUE,
    th_MaterialFireWall,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level); 
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  rval = mField.get_moab().tag_get_by_ptr(th_MaterialFireWall,&root_meshset,1,(const void**)&material_FirelWall); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::thermal_field(FieldInterface& mField) {
  PetscFunctionBegin;
  if(material_FirelWall->operator[](FW_thermal_field)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_thermal_field);

  PetscErrorCode ierr;

  ierr = mField.add_field("TEMPERATURE",H1,1,MF_ZERO); CHKERRQ(ierr);

  Range level_tets;
  ierr = mField.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(level_tets,"TEMPERATURE"); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTET,"TEMPERATURE",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"TEMPERATURE",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"TEMPERATURE",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"TEMPERATURE",1); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::spatial_problem_definition(FieldInterface& mField) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_spatial_problem_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_spatial_problem_definition);

  PetscErrorCode ierr;
  ErrorCode rval;

  //Fields
  ierr = mField.add_field("SPATIAL_POSITION",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("SPATIAL_DISPLACEMENT",H1,3,MF_ZERO); CHKERRQ(ierr);
  //ierr = mField.add_field("MESH_NODE_DISPLACEMENT",H1,3,MF_ZERO); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC",MF_ZERO); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  if(mField.check_field("TEMPERATURE")) {
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","TEMPERATURE"); CHKERRQ(ierr);
  }

  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS",MF_ZERO); CHKERRQ(ierr);
  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);

  Range level_tets;
  ierr = mField.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(level_tets,"SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(level_tets,"SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(level_tets,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //ierr = mField.add_ents_to_field_by_TETs(level_tets,"MESH_NODE_DISPLACEMENT"); CHKERRQ(ierr);
  Range blocked_tets;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-1);
  ierr = mField.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBTET,blocked_tets); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(blocked_tets,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //ierr = mField.add_ents_to_field_by_TETs(blocked_tets,"MESH_NODE_DISPLACEMENT"); CHKERRQ(ierr);

  PetscInt order;
  PetscBool flg = PETSC_TRUE;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    Tag th_set_order;
    rval = mField.get_moab().tag_get_handle("_SET_ORDER",th_set_order); CHKERR_PETSC(rval);
    const EntityHandle root_meshset = mField.get_moab().get_root_set();
    rval = mField.get_moab().tag_get_data(th_set_order,&root_meshset,1,&order); CHKERR_PETSC(rval);
  }

  //set app. order
  ierr = mField.set_field_order(0,MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"SPATIAL_DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"SPATIAL_DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"SPATIAL_DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_DISPLACEMENT",1); CHKERRQ(ierr);
  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  //ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_DISPLACEMENT",1); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode ConfigurationalFractureMechanics::material_problem_definition(FieldInterface& mField) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_material_problem_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_material_problem_definition);

  PetscErrorCode ierr;

  //Fields
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("MATERIAL_FORCE",H1,3,MF_ZERO); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("MATERIAL",MF_ZERO); CHKERRQ(ierr);

  //Define rows/cols and element data 
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //
  ierr = mField.modify_finite_element_add_field_row("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_off_field_row("MATERIAL","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MATERIAL","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  if(mField.check_field("TEMPERATURE")) {
    ierr = mField.modify_finite_element_add_field_data("MATERIAL","TEMPERATURE"); CHKERRQ(ierr);
  }

  //define problems
  ierr = mField.add_problem("MATERIAL_MECHANICS",MF_ZERO); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("MATERIAL_MECHANICS","MATERIAL"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::mesh_smoothing_problem_definition(FieldInterface& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //define problems
  ierr = mField.add_problem("MESH_SMOOTHING_PROBLEM",MF_ZERO); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("MESH_SMOOTHING_PROBLEM","MATERIAL"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("MESH_SMOOTHING_PROBLEM","CandCT_SURFACE_ELEM"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << msId;
    ierr = mField.modify_problem_add_finite_element("MESH_SMOOTHING_PROBLEM",ss2.str()); CHKERRQ(ierr);
  }

  bool cs = true;
  if(cs) {
    ierr = mField.modify_problem_add_finite_element("MESH_SMOOTHING_PROBLEM","CandCT_CRACK_SURFACE_ELEM_WITH_CRACK_FRONT"); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::coupled_problem_definition(FieldInterface& mField) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_coupled_problem_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_coupled_problem_definition);

  PetscErrorCode ierr;

  //Fields
  //ierr = mField.add_field("SPATIAL_POSITION",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("MATERIAL_FORCE",H1,3,MF_ZERO); CHKERRQ(ierr);
  //FE
  ierr = mField.add_finite_element("ELASTIC_COUPLED",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("MATERIAL_COUPLED",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("MESH_SMOOTHER",MF_ZERO); CHKERRQ(ierr);
  
  //fes definitions
  ierr = mField.modify_finite_element_add_field_row("ELASTIC_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("ELASTIC_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  if(mField.check_field("TEMPERATURE")) {
    ierr = mField.modify_finite_element_add_field_data("ELASTIC_COUPLED","TEMPERATURE"); CHKERRQ(ierr);
  }

  //
  ierr = mField.modify_finite_element_add_field_row("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("MATERIAL_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MATERIAL_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  if(mField.check_field("TEMPERATURE")) {
    ierr = mField.modify_finite_element_add_field_data("MATERIAL_COUPLED","TEMPERATURE"); CHKERRQ(ierr);
  }

  //
  ierr = mField.modify_finite_element_add_field_row("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("COUPLED_PROBLEM",MF_ZERO); CHKERRQ(ierr);
  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","ELASTIC_COUPLED"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","MESH_SMOOTHER"); CHKERRQ(ierr);

  ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","CandCT_SURFACE_ELEM"); CHKERRQ(ierr);
  bool cs = true;
  if(cs) {
    ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","CandCT_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << msId;
    ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM",ss2.str()); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::arclength_problem_definition(FieldInterface& mField) {
  PetscFunctionBegin;

  ErrorCode rval;

  if(material_FirelWall->operator[](FW_arc_lenhghat_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_arc_lenhghat_definition);

  PetscErrorCode ierr;
  
  ierr = mField.add_field("LAMBDA",NoField,1); CHKERRQ(ierr);

  ierr = mField.add_finite_element("ARC_LENGTH"); CHKERRQ(ierr);
  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);
  //elem data
  ierr = mField.modify_finite_element_add_field_data("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("ELASTIC_COUPLED","LAMBDA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC_COUPLED","LAMBDA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC_COUPLED","LAMBDA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("MATERIAL_COUPLED","LAMBDA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MATERIAL_COUPLED","LAMBDA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MATERIAL_COUPLED","LAMBDA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("MESH_SMOOTHER","LAMBDA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MESH_SMOOTHER","LAMBDA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MESH_SMOOTHER","LAMBDA"); CHKERRQ(ierr);

  ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","ARC_LENGTH"); CHKERRQ(ierr);

  //this entity will carray data for this finite element
  EntityHandle meshset_FE_ARC_LENGTH;
  rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_FE_ARC_LENGTH); CHKERR_PETSC(rval);
  //get LAMBDA field meshset
  EntityHandle meshset_field_LAMBDA = mField.get_field_meshset("LAMBDA");
  //add LAMBDA field meshset to finite element ARC_LENGTH
  rval = mField.get_moab().add_entities(meshset_FE_ARC_LENGTH,&meshset_field_LAMBDA,1); CHKERR_PETSC(rval);
  //add finite element ARC_LENGTH meshset to refinment database (all ref bit leveles)
  ierr = mField.seed_ref_level_MESHSET(meshset_FE_ARC_LENGTH,BitRefLevel().set()); CHKERRQ(ierr);
  //finally add created meshset to the ARC_LENGTH finite element
  ierr = mField.add_ents_to_finite_element_by_MESHSET(meshset_FE_ARC_LENGTH,"ARC_LENGTH"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::constrains_problem_definition(FieldInterface& mField) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_constrains_problem_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_constrains_problem_definition);

  ErrorCode rval;
  PetscErrorCode ierr;

  bool cs = true;

  //Fields
  ierr = mField.add_field("LAMBDA_SURFACE",H1,1,MF_ZERO); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    ierr = mField.add_field(ss.str(),H1,1,MF_ZERO); CHKERRQ(ierr);
  }
  //CRACK
  if(cs) {
    ierr = mField.add_field("LAMBDA_CRACK_SURFACE",H1,1,MF_ZERO); CHKERRQ(ierr);
    ierr = mField.add_field("LAMBDA_CRACK_SURFACE_WITH_CRACK_FRONT",H1,1,MF_ZERO); CHKERRQ(ierr);
  }

  //FE
  ierr = mField.add_finite_element("C_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CTC_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CandCT_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss0;
    ss0 << "C_SURFACE_ELEM_msId_" << msId;
    ostringstream ss1;
    ss1 << "CTC_SURFACE_ELEM_msId_" << msId;
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << msId;
    ierr = mField.add_finite_element(ss0.str(),MF_ZERO); CHKERRQ(ierr);
    ierr = mField.add_finite_element(ss1.str(),MF_ZERO); CHKERRQ(ierr);
    ierr = mField.add_finite_element(ss2.str(),MF_ZERO); CHKERRQ(ierr);
  }

  //CRACK
  if(cs) {
    ierr = mField.add_finite_element("C_CRACK_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.add_finite_element("CTC_CRACK_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.add_finite_element("CandCT_CRACK_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.add_finite_element("CandCT_CRACK_SURFACE_ELEM_WITH_CRACK_FRONT",MF_ZERO); CHKERRQ(ierr);
  }

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("C_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("C_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("CandCT_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("CandCT_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("CandCT_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("CandCT_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("CandCT_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("CandCT_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);

  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    ostringstream ss0;
    ss0 << "C_SURFACE_ELEM_msId_" << msId;
    ierr = mField.modify_finite_element_add_field_row(ss0.str(),ss.str()); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(ss0.str(),"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(ss0.str(),"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(ss0.str(),ss.str()); CHKERRQ(ierr);
    ostringstream ss1;
    ss1 << "CTC_SURFACE_ELEM_msId_" << msId;
    ierr = mField.modify_finite_element_add_field_row(ss1.str(),ss.str()); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(ss1.str(),ss.str()); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(ss1.str(),ss.str()); CHKERRQ(ierr);
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << msId;
    ierr = mField.modify_finite_element_add_field_row(ss2.str(),ss.str()); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(ss2.str(),"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(ss2.str(),"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(ss2.str(),ss.str()); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(ss2.str(),"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(ss2.str(),ss.str()); CHKERRQ(ierr);
  }

  //CRACK
  if(cs) {
    ierr = mField.modify_finite_element_add_field_row("C_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("C_CRACK_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("C_CRACK_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("C_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);

    ierr = mField.modify_finite_element_add_field_row("CTC_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("CTC_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("CTC_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);

    ierr = mField.modify_finite_element_add_field_row("CandCT_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("CandCT_CRACK_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("CandCT_CRACK_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("CandCT_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("CandCT_CRACK_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("CandCT_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);

    ierr = mField.modify_finite_element_add_field_row("CandCT_CRACK_SURFACE_ELEM_WITH_CRACK_FRONT","LAMBDA_CRACK_SURFACE_WITH_CRACK_FRONT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("CandCT_CRACK_SURFACE_ELEM_WITH_CRACK_FRONT","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("CandCT_CRACK_SURFACE_ELEM_WITH_CRACK_FRONT","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("CandCT_CRACK_SURFACE_ELEM_WITH_CRACK_FRONT","LAMBDA_CRACK_SURFACE_WITH_CRACK_FRONT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("CandCT_CRACK_SURFACE_ELEM_WITH_CRACK_FRONT","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("CandCT_CRACK_SURFACE_ELEM_WITH_CRACK_FRONT","LAMBDA_CRACK_SURFACE_WITH_CRACK_FRONT"); CHKERRQ(ierr);
  }

  //define problems
  ierr = mField.add_problem("CCT_ALL_MATRIX",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_problem("C_ALL_MATRIX",MF_ZERO); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("C_ALL_MATRIX","C_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_SURFACE_ELEM"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss0;
    ss0 << "C_SURFACE_ELEM_msId_" << msId;
    ostringstream ss1;
    ss1 << "CTC_SURFACE_ELEM_msId_" << msId;
    ierr = mField.modify_problem_add_finite_element("C_ALL_MATRIX",ss0.str()); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX",ss1.str()); CHKERRQ(ierr);
  }

  //CRACK //this for testing only
  /*if(cs) {
    ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("C_ALL_MATRIX","C_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
  }*/

  Range level_tris,level_edges,level_nodes;
  ierr = mField.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBVERTEX,level_nodes); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBEDGE,level_edges); CHKERRQ(ierr);

  Range corners_edges,corners_nodes,surfaces_faces;
  ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,corners_edges,true); CHKERRQ(ierr);
  corners_edges = intersect(corners_edges,level_edges);
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,corners_nodes,true); CHKERRQ(ierr);
  corners_nodes = intersect(corners_nodes,level_nodes);
  ierr = mField.get_Cubit_msId_entities_by_dimension(102,SideSet,2,surfaces_faces,true); CHKERRQ(ierr);
  surfaces_faces = intersect(surfaces_faces,level_tris);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 100 = %d\n",corners_edges.size()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of NodeSet 101 = %d\n",corners_nodes.size()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 102 = %d\n",surfaces_faces.size()); CHKERRQ(ierr);

  Range crack_surfaces_faces,crack_front_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,crack_surfaces_faces,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRQ(ierr);
  crack_surfaces_faces = intersect(crack_surfaces_faces,level_tris);
  crack_front_edges = intersect(crack_front_edges,level_edges);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack Surfaces Faces = %d\n",crack_surfaces_faces.size()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack Front Edges = %d\n",crack_front_edges.size()); CHKERRQ(ierr);

  Interface& moab = mField.get_moab();

  ierr = mField.seed_finite_elements(surfaces_faces); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(surfaces_faces,"C_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(surfaces_faces,"CTC_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(surfaces_faces,"CandCT_SURFACE_ELEM"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss0;
    ss0 << "C_SURFACE_ELEM_msId_" << msId;
    ostringstream ss1;
    ss1 << "CTC_SURFACE_ELEM_msId_" << msId;
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << msId;
    Range surfaces_faces_msId;
    ierr = mField.get_Cubit_msId_entities_by_dimension(msId,SideSet,2,surfaces_faces_msId,true); CHKERRQ(ierr);
    surfaces_faces_msId = intersect(surfaces_faces_msId,level_tris);
    ierr = mField.seed_finite_elements(surfaces_faces_msId); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet ( mdId =  %d ) = %d\n",msId,surfaces_faces_msId.size()); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(surfaces_faces_msId,ss0.str()); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(surfaces_faces_msId,ss1.str()); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(surfaces_faces_msId,ss2.str()); CHKERRQ(ierr);
  }

  if(cs) {
    ierr = mField.seed_finite_elements(crack_surfaces_faces); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(crack_surfaces_faces,"C_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(crack_surfaces_faces,"CTC_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(crack_surfaces_faces,"CandCT_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(crack_surfaces_faces,"CandCT_CRACK_SURFACE_ELEM_WITH_CRACK_FRONT"); CHKERRQ(ierr);
  }

  Range crack_front_edges_nodes;
  rval = moab.get_connectivity(crack_front_edges,crack_front_edges_nodes,true); CHKERR_PETSC(rval);
  Range corners_edges_nodes;
  rval = moab.get_connectivity(corners_edges,corners_edges_nodes,true); CHKERR_PETSC(rval);
  corners_nodes.insert(corners_edges_nodes.begin(),corners_edges_nodes.end());

  Range surfaces_nodes;
  rval = moab.get_connectivity(surfaces_faces,surfaces_nodes,true); CHKERR_PETSC(rval);
  surfaces_nodes = subtract(surfaces_nodes,corners_nodes);
  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_VERTICEs(surfaces_nodes,"LAMBDA_SURFACE"); CHKERRQ(ierr);
  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_SURFACE",1); CHKERRQ(ierr);

  //CRCAK
  if(cs) {
    Range crack_surfaces_nodes;
    rval = moab.get_connectivity(crack_surfaces_faces,crack_surfaces_nodes,true); CHKERR_PETSC(rval);
    crack_surfaces_nodes = subtract(crack_surfaces_nodes,corners_nodes);
    ierr = mField.add_ents_to_field_by_VERTICEs(crack_surfaces_nodes,"LAMBDA_CRACK_SURFACE_WITH_CRACK_FRONT"); CHKERRQ(ierr);
    crack_surfaces_nodes = subtract(crack_surfaces_nodes,crack_front_edges_nodes);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack Surfaces Nodes = %d\n",crack_surfaces_nodes.size()); CHKERRQ(ierr);
    ierr = mField.add_ents_to_field_by_VERTICEs(crack_surfaces_nodes,"LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CRACK_SURFACE",1); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CRACK_SURFACE_WITH_CRACK_FRONT",1); CHKERRQ(ierr);
  }

  
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    Range surfaces_faces_msId;
    ierr = mField.get_Cubit_msId_entities_by_dimension(msId,SideSet,2,surfaces_faces_msId,true); CHKERRQ(ierr);
    surfaces_faces_msId = intersect(surfaces_faces_msId,level_tris);
    Range surfaces_nodes_msId;
    rval = mField.get_moab().get_connectivity(surfaces_faces_msId,surfaces_nodes_msId,true); CHKERR_PETSC(rval);
    surfaces_nodes_msId = subtract(surfaces_nodes_msId,corners_nodes);
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    ierr = mField.add_ents_to_field_by_VERTICEs(surfaces_nodes_msId,ss.str()); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,ss.str(),1); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::constrains_crack_front_problem_definition(FieldInterface& mField,string problem) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_constrains_crack_front_problem_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_constrains_crack_front_problem_definition);

  PetscErrorCode ierr;
  ErrorCode rval;

  Range level_tris,level_edges,level_nodes;
  ierr = mField.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBVERTEX,level_nodes); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBEDGE,level_edges); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_field("LAMBDA_CRACKFRONT_AREA",H1,1,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("LAMBDA_CRACK_TANGENT_CONSTRAIN",H1,1,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("GRIFFITH_FORCE",H1,3,MF_ZERO); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("C_CRACKFRONT_AREA_ELEM",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CTC_CRACKFRONT_AREA_ELEM",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("dCT_CRACKFRONT_AREA_ELEM",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("C_TANGENT_ELEM",MF_ZERO); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("C_CRACKFRONT_AREA_ELEM","LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("C_CRACKFRONT_AREA_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_CRACKFRONT_AREA_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_CRACKFRONT_AREA_ELEM","LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("CTC_CRACKFRONT_AREA_ELEM","LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("CTC_CRACKFRONT_AREA_ELEM","LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("CTC_CRACKFRONT_AREA_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("CTC_CRACKFRONT_AREA_ELEM","LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("dCT_CRACKFRONT_AREA_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("dCT_CRACKFRONT_AREA_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("dCT_CRACKFRONT_AREA_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("dCT_CRACKFRONT_AREA_ELEM","LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("C_TANGENT_ELEM","LAMBDA_CRACK_TANGENT_CONSTRAIN"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("C_TANGENT_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("C_TANGENT_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("C_TANGENT_ELEM","LAMBDA_CRACK_TANGENT_CONSTRAIN"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_TANGENT_ELEM","LAMBDA_CRACK_TANGENT_CONSTRAIN"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_TANGENT_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //Problem
  ierr = mField.add_problem("C_CRACKFRONT_MATRIX",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_problem("CTC_CRACKFRONT_MATRIX",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("C_CRACKFRONT_MATRIX","C_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("CTC_CRACKFRONT_MATRIX","CTC_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);

  {
    Range crack_surfaces_faces,crack_front_edges;
    ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,crack_surfaces_faces,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRQ(ierr);
    crack_surfaces_faces = intersect(crack_surfaces_faces,level_tris);
    crack_front_edges = intersect(crack_front_edges,level_edges);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack Faces = %d\n",crack_surfaces_faces.size()); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack Front = %d\n",crack_front_edges.size()); CHKERRQ(ierr);
    Range crack_front_nodes;
    rval = mField.get_moab().get_connectivity(crack_front_edges,crack_front_nodes,true); CHKERR_PETSC(rval);
    crack_front_nodes = intersect(crack_front_nodes,level_nodes);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Front Nodes = %d\n",crack_front_nodes.size()); CHKERRQ(ierr);

    Range surfaces_faces;
    ierr = mField.get_Cubit_msId_entities_by_dimension(102,SideSet,2,surfaces_faces,true); CHKERRQ(ierr);
    Range surfaces_nodes;
    rval = mField.get_moab().get_connectivity(surfaces_faces,surfaces_nodes,true); CHKERR_PETSC(rval);
    Range tangent_crack_front_nodes = subtract(crack_front_nodes,surfaces_nodes);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Tangent Crack Front Nodes = %d\n",tangent_crack_front_nodes.size()); CHKERRQ(ierr);

    Range crack_surfaces_edge_faces;
    rval = mField.get_moab().get_adjacencies(
      crack_front_nodes,2,false,crack_surfaces_edge_faces,Interface::UNION); CHKERR_PETSC(rval);
    crack_surfaces_edge_faces = crack_surfaces_edge_faces.subset_by_type(MBTRI);
    crack_surfaces_edge_faces = intersect(crack_surfaces_edge_faces,crack_surfaces_faces);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Front Faces = %d\n",crack_surfaces_edge_faces.size()); CHKERRQ(ierr);

    Range level_tris;
    ierr = mField.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
    crack_surfaces_edge_faces = intersect(crack_surfaces_edge_faces,level_tris);

    ierr = mField.seed_finite_elements(crack_surfaces_edge_faces); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(crack_surfaces_edge_faces,"C_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(crack_surfaces_edge_faces,"CTC_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(crack_surfaces_edge_faces,"dCT_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(crack_surfaces_edge_faces,"C_TANGENT_ELEM"); CHKERRQ(ierr);

    Range level_edges;
    ierr = mField.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBEDGE,level_edges); CHKERRQ(ierr);
    crack_front_edges = intersect(crack_front_edges,level_edges);
    ierr = mField.seed_finite_elements(crack_front_edges); CHKERRQ(ierr);

    //add entitities (by tets) to the field
    ierr = mField.add_ents_to_field_by_VERTICEs(crack_front_nodes,"LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_field_by_VERTICEs(tangent_crack_front_nodes,"LAMBDA_CRACK_TANGENT_CONSTRAIN"); CHKERRQ(ierr);
  }

  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CRACKFRONT_AREA",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CRACK_TANGENT_CONSTRAIN",1); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element(problem,"dCT_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);

  if(problem == "COUPLED_PROBLEM") {
    ierr = mField.modify_problem_add_finite_element(problem,"C_TANGENT_ELEM"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("MESH_SMOOTHER","LAMBDA_CRACK_TANGENT_CONSTRAIN"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("MESH_SMOOTHER","LAMBDA_CRACK_TANGENT_CONSTRAIN"); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::spatial_partition_problems(FieldInterface& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS",1); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS",1); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::material_partition_problems(FieldInterface& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition MATERIAL_MECHANICS
  ierr = mField.partition_problem("MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("MATERIAL_MECHANICS"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::coupled_partition_problems(FieldInterface& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = mField.partition_problem("COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("COUPLED_PROBLEM"); CHKERRQ(ierr);

  ierr = mField.partition_problem("MESH_SMOOTHING_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("MESH_SMOOTHING_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("MESH_SMOOTHING_PROBLEM"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::constrains_partition_problems(FieldInterface& mField,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = mField.simple_partition_problem("CCT_ALL_MATRIX",0); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  //partition
  ierr = mField.compose_problem("C_ALL_MATRIX","CCT_ALL_MATRIX",false,problem,true); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("C_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("C_ALL_MATRIX"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::crackfront_partition_problems(FieldInterface& mField,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = mField.simple_partition_problem("CTC_CRACKFRONT_MATRIX",0); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("CTC_CRACKFRONT_MATRIX"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("CTC_CRACKFRONT_MATRIX"); CHKERRQ(ierr);
  //partition
  ierr = mField.compose_problem("C_CRACKFRONT_MATRIX","CTC_CRACKFRONT_MATRIX",false,problem,true); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("C_CRACKFRONT_MATRIX"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("C_CRACKFRONT_MATRIX"); CHKERRQ(ierr);

  PetscFunctionReturn(0);

}

PetscErrorCode ConfigurationalFractureMechanics::set_spatial_positions(FieldInterface& mField) {
  PetscFunctionBegin;
  if(material_FirelWall->operator[](FW_set_spatial_positions)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_set_spatial_positions);

  ErrorCode rval;

  EntityHandle node = 0;
  double coords[3];
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"SPATIAL_POSITION",dof_ptr)) {
      if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
      EntityHandle ent = dof_ptr->get_ent();
      int dof_rank = dof_ptr->get_dof_rank();
      double &fval = dof_ptr->get_FieldData();
      if(node!=ent) {
	rval = mField.get_moab().get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
      }
      fval = coords[dof_rank];
  }

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::set_material_positions(FieldInterface& mField) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_set_material_positions)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_set_material_positions);

  ErrorCode rval;

  EntityHandle node = 0;
  double coords[3];
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"MESH_NODE_POSITIONS",dof_ptr)) {
    if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
    EntityHandle ent = dof_ptr->get_ent();
    int dof_rank = dof_ptr->get_dof_rank();
    double &fval = dof_ptr->get_FieldData();
    if(node!=ent) {
	rval = mField.get_moab().get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
    }
    fval = coords[dof_rank];
  }

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::set_coordinates_from_material_solution(FieldInterface& mField) {
  PetscFunctionBegin;
  //PetscErrorCode ierr;
  ErrorCode rval;
  double coords[3];
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"MESH_NODE_POSITIONS",dof_ptr)) {
    if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
    EntityHandle ent = dof_ptr->get_ent();
    int dof_rank = dof_ptr->get_dof_rank();
    double fval = dof_ptr->get_FieldData();
    rval = mField.get_moab().get_coords(&ent,1,coords); CHKERR_PETSC(rval);
    coords[dof_rank] = fval;
    rval = mField.get_moab().set_coords(&ent,1,coords); CHKERR_PETSC(rval);
  }
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"LAMBDA_SURFACE",dof_ptr)) {
    dof_ptr->get_FieldData() = 0;
  }
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"LAMBDA_CRACK_SURFACE",dof_ptr)) {
    dof_ptr->get_FieldData() = 0;
  }
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"LAMBDA_CRACK_TANGENT_CONSTRAIN",dof_ptr)) {
    dof_ptr->get_FieldData() = 0;
  }
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,ss.str(),dof_ptr)) {
      dof_ptr->get_FieldData() = 0;
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::solve_spatial_problem(FieldInterface& mField,SNES *snes,bool postproc) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //create matrices
  Vec F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  DirihletBCMethod_DriverComplexForLazy myDirihletBCSpatial(mField,"ELASTIC_MECHANICS","SPATIAL_POSITION");
  ierr = myDirihletBCSpatial.Init(); CHKERRQ(ierr);

  const double YoungModulus = 1;
  const double PoissonRatio = 0.;
  NL_ElasticFEMethod MySpatialFE(mField,&myDirihletBCSpatial,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  PetscBool flg = PETSC_TRUE;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_thermal_expansion",&MySpatialFE.thermal_expansion,&flg); CHKERRQ(ierr);


  SnesCtx snes_ctx(mField,"ELASTIC_MECHANICS");
  
  ierr = SNESSetApplicationContext(*snes,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(*snes,F,SnesRhs,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(*snes,Aij,Aij,SnesMat,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(*snes); CHKERRQ(ierr);

  SnesCtx::loops_to_do_type& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ELASTIC",&MySpatialFE));
  SnesCtx::loops_to_do_type& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ELASTIC",&MySpatialFE));

  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = mField.set_local_VecCreateGhost("ELASTIC_MECHANICS",Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = SNESSolve(*snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(*snes,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  if(postproc) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Save tags with nodal postions for spatial problem\n"); CHKERRQ(ierr);
    PostProcVertexMethod ent_method(mField.get_moab(),"SPATIAL_POSITION");
    ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Col,ent_method); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Save tags with residuals for spatial problem\n"); CHKERRQ(ierr);
    PostProcVertexMethod ent_method_res(mField.get_moab(),"SPATIAL_POSITION",F,"SPATIAL_RESIDUAL");
    ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Col,ent_method_res); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Save tags on refined mesh with stresses\n"); CHKERRQ(ierr);
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    if(pcomm->rank()==0) {
      if(fe_post_proc_stresses_method!=NULL) delete fe_post_proc_stresses_method;
      fe_post_proc_stresses_method = new PostProcStressNonLinearElasticity(mField.get_moab(),MySpatialFE);
      fe_post_proc_stresses_method->do_broadcast = false;
      ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",*fe_post_proc_stresses_method,0,pcomm->size());  CHKERRQ(ierr);
    }
  }

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Spatial problem solved\n"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::surface_projection_data(FieldInterface& mField,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  ierr = delete_surface_projection_data(mField); CHKERRQ(ierr);

  //C_ALL_MATRIX is problem used to constrain surface constrain matrices
  if(projSurfaceCtx==NULL) {
    projSurfaceCtx = new matPROJ_ctx(mField,problem,"C_ALL_MATRIX");
    ierr = mField.MatCreateMPIAIJWithArrays("C_ALL_MATRIX",&projSurfaceCtx->C); CHKERRQ(ierr);
  }

  Range corners_edges,corners_nodes;
  ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,corners_edges,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,corners_nodes,true); CHKERRQ(ierr);
  ErrorCode rval;
  Range corners_edges_nodes;
  rval = mField.get_moab().get_connectivity(corners_edges,corners_edges_nodes,true); CHKERR_PETSC(rval);
  corners_nodes.merge(corners_edges_nodes);
  Range nodes_to_block;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-1);
  ierr = mField.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBVERTEX,nodes_to_block); CHKERRQ(ierr);
  corners_nodes.merge(nodes_to_block);
  CubitDisplacementDirihletBC_Coupled myDirihletBC(mField,"C_ALL_MATRIX",corners_nodes);

  //Loops over body surface (CFE_SURFACE) and crack surface (CFE_CRACK_SURFACE)
  C_SURFACE_FEMethod CFE_SURFACE(mField,&myDirihletBC,projSurfaceCtx->C);
  C_SURFACE_FEMethod CFE_CRACK_SURFACE(mField,&myDirihletBC,projSurfaceCtx->C,"LAMBDA_CRACK_SURFACE");

  map<int,C_SURFACE_FEMethod*> CFE_SURFACE_msId_ptr;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    CFE_SURFACE_msId_ptr[msId] = new C_SURFACE_FEMethod(mField,&myDirihletBC,projSurfaceCtx->C,ss.str());
  }

  ierr = MatSetOption(projSurfaceCtx->C,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(projSurfaceCtx->C,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

  ierr = MatZeroEntries(projSurfaceCtx->C); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_SURFACE_ELEM",CFE_SURFACE);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_CRACK_SURFACE_ELEM",CFE_CRACK_SURFACE);  CHKERRQ(ierr);
  for(map<int,C_SURFACE_FEMethod*>::iterator mit = CFE_SURFACE_msId_ptr.begin();
    mit!=CFE_SURFACE_msId_ptr.end();mit++) {
    ostringstream ss0;
    ss0 << "C_SURFACE_ELEM_msId_" << mit->first;
    ierr = mField.loop_finite_elements("C_ALL_MATRIX",ss0.str(),*(mit->second));  CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(projSurfaceCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(projSurfaceCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  for(map<int,C_SURFACE_FEMethod*>::iterator mit = CFE_SURFACE_msId_ptr.begin();
    mit!=CFE_SURFACE_msId_ptr.end();mit++) {
    delete mit->second;
  }

  /*{
    //Matrix View
    MatView(projSurfaceCtx->C,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
    std::string wait;
    std::cin >> wait;
  }*/

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::delete_surface_projection_data(FieldInterface& mField) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  if(projSurfaceCtx!=NULL) {
    ierr = projSurfaceCtx->DestroyQorP(); CHKERRQ(ierr);
    ierr = projSurfaceCtx->DestroyQTKQ(); CHKERRQ(ierr);
    ierr = MatDestroy(&projSurfaceCtx->C); CHKERRQ(ierr);
    delete projSurfaceCtx;
    projSurfaceCtx = NULL;
  }

  PetscFunctionReturn(0);
}


PetscErrorCode ConfigurationalFractureMechanics::project_force_vector(FieldInterface& mField,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  //ErrorCode rval;

  Vec F_Material;
  ierr = mField.VecCreateGhost(problem,Row,&F_Material); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost(problem,"MESH_NODE_POSITIONS","MATERIAL_FORCE",Row,F_Material,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  int M,m;
  ierr = VecGetSize(F_Material,&M); CHKERRQ(ierr);
  ierr = VecGetLocalSize(F_Material,&m); CHKERRQ(ierr);
  Mat Q;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,projSurfaceCtx,&Q); CHKERRQ(ierr);
  ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))matQ_mult_shell); CHKERRQ(ierr);

  Vec QTF_Material;
  ierr = VecDuplicate(F_Material,&QTF_Material); CHKERRQ(ierr);
  ierr = MatMult(Q,F_Material,QTF_Material); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_material(mField.get_moab(),"MESH_NODE_POSITIONS",QTF_Material,"MATERIAL_FORCE_PROJECTED");
  ierr = mField.loop_dofs(problem,"MESH_NODE_POSITIONS",Row,ent_method_material); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost(problem,"MESH_NODE_POSITIONS","MATERIAL_FORCE",Row,QTF_Material,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  double nrm2_material;
  ierr = VecNorm(QTF_Material,NORM_2,&nrm2_material);   CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"nrm2_QTF_Material = %6.4e\n",nrm2_material); CHKERRQ(ierr);

  ierr = VecDestroy(&F_Material); CHKERRQ(ierr);
  ierr = MatDestroy(&Q); CHKERRQ(ierr);
  ierr = VecDestroy(&QTF_Material); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode ConfigurationalFractureMechanics::project_form_th_projection_tag(FieldInterface& mField,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;

  Vec D,dD;
  ierr = mField.VecCreateGhost(problem,Col,&D); CHKERRQ(ierr);
  ierr = VecDuplicate(D,&dD); CHKERRQ(ierr);
  ierr = VecZeroEntries(dD); CHKERRQ(ierr);
  ierr = VecSetOption(dD,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);

  ierr = mField.set_local_VecCreateGhost(problem,Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost(problem,Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  Tag th_projection;
  rval = mField.get_moab().tag_get_handle("PROJECTION_CRACK_SURFACE",th_projection); CHKERR_PETSC(rval);

  const MoFEMProblem *problem_ptr;
  ierr = mField.get_problem(problem,&problem_ptr); CHKERRQ(ierr);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);

  Range crack_front_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRQ(ierr);
  Range crack_front_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_front_edges,crack_front_edges_nodes,true); CHKERR_PETSC(rval);
  Range::iterator nit = crack_front_edges_nodes.begin();
  for(;nit!=crack_front_edges_nodes.end();nit++) {

    NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator dit,hi_dit;
    dit = problem_ptr->numered_dofs_cols.get<Composite_Name_Ent_And_Part_mi_tag>().lower_bound(
      boost::make_tuple("MESH_NODE_POSITIONS",*nit,pcomm->rank()));
    hi_dit = problem_ptr->numered_dofs_cols.get<Composite_Name_Ent_And_Part_mi_tag>().upper_bound(
      boost::make_tuple("MESH_NODE_POSITIONS",*nit,pcomm->rank()));
    if(distance(dit,hi_dit)==0) continue;
    ublas::vector<PetscInt,ublas::bounded_array<PetscInt,3> > indices;
    indices.resize(3);
    fill(indices.begin(),indices.end(),-1);
    ublas::vector<double,ublas::bounded_array<double,3> > new_coords;
    new_coords.resize(3);
    rval = mField.get_moab().tag_get_data(th_projection,&*nit,1,&*new_coords.data().begin()); CHKERR_PETSC(rval);
    ublas::vector<double,ublas::bounded_array<double,3> > delta;
    delta.resize(3);
    for(;dit!=hi_dit;dit++) {
      indices[dit->get_dof_rank()] = dit->get_petsc_gloabl_dof_idx();
      delta[dit->get_dof_rank()] = new_coords[dit->get_dof_rank()] - dit->get_FieldData();
    }
    ierr = VecSetValues(dD,3,&*indices.data().begin(),&*delta.data().begin(),INSERT_VALUES); CHKERRQ(ierr);

  }

  ierr = VecAssemblyBegin(dD); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(dD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  int M,m;
  ierr = VecGetSize(dD,&M); CHKERRQ(ierr);
  ierr = VecGetLocalSize(dD,&m); CHKERRQ(ierr);
  Mat Q;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,projSurfaceCtx,&Q); CHKERRQ(ierr);
  ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))matQ_mult_shell); CHKERRQ(ierr);

  Vec QTdD;
  ierr = VecDuplicate(D,&QTdD); CHKERRQ(ierr);
  ierr = MatMult(Q,dD,QTdD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(QTdD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(QTdD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = VecAXPY(D,1.,QTdD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost(problem,Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_res(mField.get_moab(),"MESH_NODE_POSITIONS",D,"PROJECTION_CRACK_SURFACE");
  ierr = mField.loop_dofs(problem,"MESH_NODE_POSITIONS",Col,ent_method_res); CHKERRQ(ierr);

  VecDestroy(&D);
  VecDestroy(&dD);
  VecDestroy(&QTdD);
  MatDestroy(&Q);

  PetscFunctionReturn(0);
}


PetscErrorCode ConfigurationalFractureMechanics::front_projection_data(FieldInterface& mField,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  
  ierr = delete_front_projection_data(mField); CHKERRQ(ierr);

  if(projFrontCtx==NULL) {
    projFrontCtx = new matPROJ_ctx(mField,problem,"C_CRACKFRONT_MATRIX");
    ierr = mField.MatCreateMPIAIJWithArrays("C_CRACKFRONT_MATRIX",&projFrontCtx->C); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::delete_front_projection_data(FieldInterface& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  if(projFrontCtx!=NULL) {
    ierr = projFrontCtx->DestroyQorP(); CHKERRQ(ierr);
    ierr = projFrontCtx->DestroyQTKQ(); CHKERRQ(ierr);
    ierr = MatDestroy(&projFrontCtx->C); CHKERRQ(ierr);
    delete projFrontCtx;
    projFrontCtx = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::griffith_force_vector(FieldInterface& mField,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  Vec LambdaVec,GriffithForceVec;
  ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",Row,&LambdaVec); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",Col,&GriffithForceVec); CHKERRQ(ierr);

  double gc;
  PetscBool flg;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
  //if(flg != PETSC_TRUE) {
    //SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
  //}
  if(flg == PETSC_TRUE) {
    ierr = VecSet(LambdaVec,gc); CHKERRQ(ierr);
  } else {
    ierr = VecSet(LambdaVec,1); CHKERRQ(ierr);
  }

  Mat Q;
  {
    int M,m;
    ierr = VecGetSize(GriffithForceVec,&M); CHKERRQ(ierr);
    ierr = VecGetLocalSize(GriffithForceVec,&m); CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,projSurfaceCtx,&Q); CHKERRQ(ierr);
    ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))matQ_mult_shell); CHKERRQ(ierr);
  }

  Range corners_edges,corners_nodes;
  ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,corners_edges,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,corners_nodes,true); CHKERRQ(ierr);
  ErrorCode rval;
  Range corners_edges_nodes;
  rval = mField.get_moab().get_connectivity(corners_edges,corners_edges_nodes,true); CHKERR_PETSC(rval);
  corners_nodes.merge(corners_edges_nodes);
  Range blocked_nodes;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-1);
  ierr = mField.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBVERTEX,blocked_nodes); CHKERRQ(ierr);
  corners_nodes.merge(blocked_nodes);
  CubitDisplacementDirihletBC_Coupled myDirihletBC(mField,"C_CRACKFRONT_MATRIX",corners_nodes);

  C_CONSTANT_AREA_FEMethod C_AREA_ELEM(mField,projFrontCtx->C,Q,"LAMBDA_CRACKFRONT_AREA");

  ierr = MatSetOption(projFrontCtx->C,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(projFrontCtx->C,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

  ierr = MatZeroEntries(projFrontCtx->C); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_CRACKFRONT_MATRIX","C_CRACKFRONT_AREA_ELEM",C_AREA_ELEM);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(projFrontCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(projFrontCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = MatMultTranspose(projFrontCtx->C,LambdaVec,GriffithForceVec); CHKERRQ(ierr);

  Vec QTGriffithForceVec;
  ierr = VecDuplicate(GriffithForceVec,&QTGriffithForceVec); CHKERRQ(ierr);
  ierr = MatMult(Q,GriffithForceVec,QTGriffithForceVec); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(QTGriffithForceVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(QTGriffithForceVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost(
    problem,"MESH_NODE_POSITIONS","GRIFFITH_FORCE",Row,QTGriffithForceVec,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_area(mField.get_moab(),"MESH_NODE_POSITIONS",QTGriffithForceVec,"GRIFFITH_FORCE");
  ierr = mField.loop_dofs(problem,"MESH_NODE_POSITIONS",Col,ent_method_area); CHKERRQ(ierr);

  double nrm2_griffith_force;
  ierr = VecNorm(QTGriffithForceVec,NORM_2,&nrm2_griffith_force);   CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"nrm2_QTGriffithForceVec = %6.4e\n",nrm2_griffith_force); CHKERRQ(ierr);

  //Tangent front froce
  C_FRONT_TANGENT_FEMethod C_TANGENT_ELEM(mField,projFrontCtx->C,Q,"LAMBDA_CRACKFRONT_AREA");

  ierr = MatZeroEntries(projFrontCtx->C); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_CRACKFRONT_MATRIX","C_CRACKFRONT_AREA_ELEM",C_TANGENT_ELEM);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(projFrontCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(projFrontCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = MatMultTranspose(projFrontCtx->C,LambdaVec,GriffithForceVec); CHKERRQ(ierr);
  ierr = MatMult(Q,GriffithForceVec,QTGriffithForceVec); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(QTGriffithForceVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(QTGriffithForceVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_tangent(
    mField.get_moab(),"MESH_NODE_POSITIONS",QTGriffithForceVec,"GRIFFITH_TANGENT_FORCE");
  ierr = mField.loop_dofs(problem,"MESH_NODE_POSITIONS",Col,ent_method_tangent); CHKERRQ(ierr);

  ierr = MatDestroy(&Q); CHKERRQ(ierr);
  ierr = VecDestroy(&LambdaVec); CHKERRQ(ierr);
  ierr = VecDestroy(&GriffithForceVec); CHKERRQ(ierr);
  ierr = VecDestroy(&QTGriffithForceVec); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::griffith_g(FieldInterface& mField,string problem) {
  PetscFunctionBegin;
  
  PetscErrorCode ierr;
  //ErrorCode rval;

  Vec F_Material;
  ierr = mField.VecCreateGhost(problem,Row,&F_Material); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost(problem,"MESH_NODE_POSITIONS","MATERIAL_FORCE",Row,F_Material,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  Vec F_Griffith;
  ierr = mField.VecCreateGhost(problem,Row,&F_Griffith); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost(problem,"MESH_NODE_POSITIONS","GRIFFITH_FORCE",Row,F_Griffith,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  Vec LambdaVec;
  ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",Row,&LambdaVec); CHKERRQ(ierr);

  Mat Q;
  {
    int M,m;
    ierr = VecGetSize(F_Material,&M); CHKERRQ(ierr);
    ierr = VecGetLocalSize(F_Material,&m); CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,projSurfaceCtx,&Q); CHKERRQ(ierr);
    ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))matQ_mult_shell); CHKERRQ(ierr);
  }

  ErrorCode rval;
  Range corners_edges,corners_nodes;
  ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,corners_edges,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,corners_nodes,true); CHKERRQ(ierr);
  Range corners_edges_nodes;
  rval = mField.get_moab().get_connectivity(corners_edges,corners_edges_nodes,true); CHKERR_PETSC(rval);
  corners_nodes.insert(corners_edges_nodes.begin(),corners_edges_nodes.end());
  Range blocked_nodes;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-1);
  ierr = mField.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBVERTEX,blocked_nodes); CHKERRQ(ierr);
  corners_nodes.merge(blocked_nodes);

  CubitDisplacementDirihletBC_Coupled myDirihletBC(mField,"C_CRACKFRONT_MATRIX",corners_nodes);
  C_CONSTANT_AREA_FEMethod C_AREA_ELEM(mField,projFrontCtx->C,Q,"LAMBDA_CRACKFRONT_AREA");

  ierr = MatSetOption(projFrontCtx->C,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(projFrontCtx->C,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

  ierr = MatZeroEntries(projFrontCtx->C); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_CRACKFRONT_MATRIX","C_CRACKFRONT_AREA_ELEM",C_AREA_ELEM);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(projFrontCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(projFrontCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  /*{
    //Matrix View
    MatView(projFrontCtx->C,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
    std::string wait;
    std::cin >> wait;
  }*/

  // R = CT*(CC^T)^(-1) [ unit m/m^(-2) = 1/m ] [ 0.5/0.25 = 2 ]
  // R^T = (CC^T)^(-T)C [ unit m/m^(-2) = 1/m ] 
  Mat RT;
  {
    int N,n;
    int M,m;
    ierr = VecGetSize(F_Material,&N); CHKERRQ(ierr);
    ierr = VecGetLocalSize(F_Material,&n); CHKERRQ(ierr);
    ierr = VecGetSize(LambdaVec,&M); CHKERRQ(ierr);
    ierr = VecGetLocalSize(LambdaVec,&m); CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,projFrontCtx,&RT); CHKERRQ(ierr);
    ierr = MatShellSetOperation(RT,MATOP_MULT,(void(*)(void))matRT_mult_shell); CHKERRQ(ierr);
  }
  
  double gc = 1;
  PetscBool flg;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
  //if(flg != PETSC_TRUE) {
  //    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
  //}

  ierr = projFrontCtx->InitQorP(F_Material); CHKERRQ(ierr);
  // unit of LambdaVec [ N * 1/m = N*m/m^2 = J/m^2 ]
  ierr = VecScale(F_Material,-1./gc); CHKERRQ(ierr);
  ierr = MatMult(RT,F_Material,LambdaVec); CHKERRQ(ierr);  
  ierr = VecScale(F_Material,gc); CHKERRQ(ierr);
  ierr = VecScale(LambdaVec,gc); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(LambdaVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(LambdaVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost("C_CRACKFRONT_MATRIX",Row,LambdaVec,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  Vec JVec;
  ierr = VecDuplicate(LambdaVec,&JVec); CHKERRQ(ierr);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);

  ublas::vector<double> coords;
  coords.resize(3);
  PetscPrintf(PETSC_COMM_WORLD,"\n\ngriffith force\n\n");

  Tag th_g;
  const double def_val = 0;
  rval = mField.get_moab().tag_get_handle("G",1,MB_TYPE_DOUBLE,th_g,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERR_THROW(rval);

  map_ent_g.clear();
  map_ent_j.clear();
  const MoFEMProblem *problem_ptr;
  ierr = mField.get_problem("C_CRACKFRONT_MATRIX",&problem_ptr); CHKERRQ(ierr);
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(problem_ptr,"LAMBDA_CRACKFRONT_AREA",diit)) {
    EntityHandle ent = diit->get_ent();
    rval = mField.get_moab().get_coords(&ent,1,&*coords.data().begin()); CHKERR_PETSC(rval);
    int dd = 0;
    ublas::vector<double,ublas::bounded_array<double,9> > material_force(3);
    for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,"MATERIAL_FORCE",diit->get_ent(),diiit)) {
      material_force[diiit->get_dof_rank()] = diiit->get_FieldData();
      dd++;
    }
    if(dd != 3) SETERRQ1(PETSC_COMM_SELF,1,"can not find material force vector at node %ld",diit->get_ent());
    ublas::vector<double,ublas::bounded_array<double,9> > griffith_force(3);
    for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,"GRIFFITH_FORCE",diit->get_ent(),diiit)) {
      griffith_force[diiit->get_dof_rank()] = diiit->get_FieldData();
      dd++;
    }
    if(dd != 6) SETERRQ1(PETSC_COMM_SELF,1,"can not find griffith force vector at node %ld",diit->get_ent());
    double j = norm_2(material_force)/(norm_2(griffith_force)/gc);
    if(diit->get_part()==pcomm->rank()) {
      ierr = VecSetValue(JVec,diit->get_petsc_gloabl_dof_idx(),j,INSERT_VALUES); CHKERRQ(ierr);
    }
    ostringstream ss;
    ss << "griffith force at ";
    ss << "ent " << diit->get_ent();
    ss << "\tcoords";
    ss << " " << setw(10) << setprecision(4) << coords[0];
    ss << " " << setw(10) << setprecision(4) << coords[1];
    ss << " " << setw(10) << setprecision(4) << coords[2];
    ss << "\t\tg " << scientific << setprecision(4) << diit->get_FieldData();
    ss << " / " << scientific << setprecision(4) << j;
    ss << " ( " << scientific << setprecision(4) << diit->get_FieldData()/j << " )";
    if(flg == PETSC_TRUE) {
      ss << "\t relative error (gc-g)/gc " << scientific << setprecision(4) << (gc-diit->get_FieldData())/gc;
      ss << " / " << scientific << setprecision(4) << (gc-j)/gc;
    }
    ss << endl; 
    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
    double val = diit->get_FieldData();
    map_ent_g[ent] = val;
    map_ent_j[ent] = j;
    rval = mField.get_moab().tag_set_data(th_g,&ent,1,&val); CHKERR_PETSC(rval);
  }

  ierr = VecSum(LambdaVec,&ave_g); CHKERRQ(ierr);
  ierr = VecMin(LambdaVec,PETSC_NULL,&min_g); CHKERRQ(ierr);
  ierr = VecMax(LambdaVec,PETSC_NULL,&max_g); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(JVec); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(JVec); CHKERRQ(ierr);
  ierr = VecSum(JVec,&ave_j); CHKERRQ(ierr);
  ierr = VecMin(JVec,PETSC_NULL,&min_j); CHKERRQ(ierr);
  ierr = VecMax(JVec,PETSC_NULL,&max_j); CHKERRQ(ierr);
  ierr = VecDestroy(&JVec); CHKERRQ(ierr);

  Range crackSurfacesFaces;
  ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,crackSurfacesFaces,true); CHKERRQ(ierr);
  Range level_tris;
  ierr = mField.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
  crackSurfacesFaces = intersect(crackSurfacesFaces,level_tris);

  double aRea = 0;
  ublas::vector<double,ublas::bounded_array<double,6> > diffNTRI(6);
  ShapeDiffMBTRI(&diffNTRI.data()[0]);
  for(Range::iterator fit = crackSurfacesFaces.begin();fit!=crackSurfacesFaces.end();fit++) {
    const EntityHandle* conn; 
    int num_nodes; 
    rval = mField.get_moab().get_connectivity(*fit,conn,num_nodes,true); CHKERR_PETSC(rval);
    ublas::vector<double,ublas::bounded_array<double,9> > dofsX(9);
    ublas::vector<double,ublas::bounded_array<double,3> > normal(3);
    for(int nn = 0;nn<num_nodes; nn++) {
      for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,"MESH_NODE_POSITIONS",conn[nn],dit)) {
	dofsX[nn*3+dit->get_dof_rank()] = dit->get_FieldData();
      }
    }
    ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&dofsX.data()[0],&normal.data()[0]); CHKERRQ(ierr);
    //crack surface area is a half of crack top and bottom body surface
    aRea += norm_2(normal)*0.25;
  }

  {
    int N;
    ierr = VecGetSize(LambdaVec,&N); CHKERRQ(ierr);
    ave_g /= N;
    ave_j /= N;
  }
  PetscPrintf(PETSC_COMM_WORLD,"\naverage griffith force %6.4e / %6.4e Crack surface area %6.4e\n",ave_g,ave_j,aRea);
  PetscPrintf(PETSC_COMM_WORLD,"\n\n");

  PostProcVertexMethod ent_method(mField.get_moab(),"LAMBDA_CRACKFRONT_AREA");
  ierr = mField.loop_dofs("C_CRACKFRONT_MATRIX","LAMBDA_CRACKFRONT_AREA",Row,ent_method); CHKERRQ(ierr);

  ierr = MatDestroy(&Q); CHKERRQ(ierr);
  ierr = MatDestroy(&RT); CHKERRQ(ierr);
  ierr = VecDestroy(&F_Griffith); CHKERRQ(ierr);
  ierr = VecDestroy(&F_Material); CHKERRQ(ierr);
  ierr = VecDestroy(&LambdaVec); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode MySnesConvernceTest(SNES snes,int it,double xnorm,double gnorm,double fnorm,SNESConvergedReason *reason,void *void_ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = SNESConvergedDefault(snes,it,xnorm,gnorm,fnorm,reason,PETSC_NULL); CHKERRQ(ierr);
  const PetscReal div = 10e3;
  if(fnorm > div) {
    *reason = SNES_DIVERGED_FUNCTION_DOMAIN;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::solve_mesh_smooting_problem(FieldInterface& mField,SNES *snes) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  Range corners_edges,corners_nodes,crack_front_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,corners_edges,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,corners_nodes,true); CHKERRQ(ierr);
  Range corners_edges_nodes;
  rval = mField.get_moab().get_connectivity(corners_edges,corners_edges_nodes,true); CHKERR_PETSC(rval);
  corners_nodes.merge(corners_edges_nodes);


  Range nodes_to_block;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-1);
  ierr = mField.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBVERTEX,nodes_to_block); CHKERRQ(ierr);
  corners_nodes.merge(nodes_to_block);

  CubitDisplacementDirihletBC_Coupled myDirihletBC(mField,"MESH_SMOOTHING_PROBLEM",corners_nodes);
  ierr = myDirihletBC.Init(); CHKERRQ(ierr);

  NL_MeshSmootherCoupled MyMeshSmoother(mField,&myDirihletBC);
  C_SURFACE_FEMethod_ForSnes MySurfaceConstrainsBodySurface(mField,&myDirihletBC,"LAMBDA_SURFACE");
  C_SURFACE_FEMethod_ForSnes MySurfaceConstrainsCrackSurface(mField,&myDirihletBC,"LAMBDA_CRACK_SURFACE_WITH_CRACK_FRONT");
  MySurfaceConstrainsCrackSurface.use_projection_from_crack_front = true;

  MySurfaceConstrainsBodySurface.nonlinear = true;
  MySurfaceConstrainsCrackSurface.nonlinear = true;


  //add additional surfaces
  map<int,C_SURFACE_FEMethod_ForSnes*> C_SURFACE_FEMethod_ForSnes_msId_ptr;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    C_SURFACE_FEMethod_ForSnes_msId_ptr[msId] = new C_SURFACE_FEMethod_ForSnes(mField,&myDirihletBC,ss.str());
    C_SURFACE_FEMethod_ForSnes_msId_ptr[msId]->nonlinear = true;
  }

  SnesCtx snes_ctx(mField,"MESH_SMOOTHING_PROBLEM");

  struct SNES_set_global_VecCreateGhost: public FieldInterface::FEMethod {
    FieldInterface& mField;
    SNES_set_global_VecCreateGhost(FieldInterface& _mField): mField(_mField) {} 
    PetscErrorCode ierr;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      ierr = mField.set_global_VecCreateGhost(
	"MESH_SMOOTHING_PROBLEM",Col,snes_x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  };
  SNES_set_global_VecCreateGhost my_set_global_VecCreateGhost(mField);

  struct SNES_dirihlet_bc: public FieldInterface::FEMethod {

    FieldInterface& mField;
    BaseDirihletBC *dirihlet_bc_method_ptr;
    SNES_dirihlet_bc(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr):
      mField(_mField),dirihlet_bc_method_ptr(_dirihlet_bc_method_ptr) {} 

    PetscErrorCode ierr;

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      switch(snes_ctx) {
	case ctx_SNESSetFunction: { 
	  ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	  ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,snes_f); CHKERRQ(ierr);
	}
	break;
	case ctx_SNESSetJacobian: {
	  ierr = MatAssemblyBegin(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	  ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,*snes_B); CHKERRQ(ierr);
	}
	break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      PetscFunctionReturn(0);
    }

  };
  SNES_dirihlet_bc my_set_dirihlet_bc(mField,&myDirihletBC);

  SnesCtx::loops_to_do_type& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
  snes_ctx.get_preProcess_to_do_Rhs().push_back(&my_set_global_VecCreateGhost);
  snes_ctx.get_postProcess_to_do_Rhs().push_back(&my_set_dirihlet_bc);
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MATERIAL",&MyMeshSmoother));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("CandCT_SURFACE_ELEM",&MySurfaceConstrainsBodySurface));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM_WITH_CRACK_FRONT",&MySurfaceConstrainsCrackSurface));

  SnesCtx::loops_to_do_type& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
  snes_ctx.get_postProcess_to_do_Mat().push_back(&my_set_dirihlet_bc);
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MATERIAL",&MyMeshSmoother));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("CandCT_SURFACE_ELEM",&MySurfaceConstrainsBodySurface));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM_WITH_CRACK_FRONT",&MySurfaceConstrainsCrackSurface));

  for(map<int,C_SURFACE_FEMethod_ForSnes*>::iterator mit = C_SURFACE_FEMethod_ForSnes_msId_ptr.begin();
    mit!=C_SURFACE_FEMethod_ForSnes_msId_ptr.end();mit++) {
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << mit->first;
    loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type(ss2.str(),mit->second));
    loops_to_do_Mat.push_back(SnesCtx::loop_pair_type(ss2.str(),mit->second));
  }

  Vec F,D;
  ierr = mField.VecCreateGhost("MESH_SMOOTHING_PROBLEM",Row,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("MESH_SMOOTHING_PROBLEM",Row,&D); CHKERRQ(ierr);
  Mat K;
  ierr = mField.MatCreateMPIAIJWithArrays("MESH_SMOOTHING_PROBLEM",&K); CHKERRQ(ierr);

  ierr = SNESSetApplicationContext(*snes,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(*snes,F,SnesRhs,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(*snes,K,K,SnesMat,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetConvergenceTest(*snes,MySnesConvernceTest,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);

  ierr = mField.set_local_VecCreateGhost("MESH_SMOOTHING_PROBLEM",Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = SNESSolve(*snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(*snes,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost(
    "MESH_SMOOTHING_PROBLEM",Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&K); CHKERRQ(ierr);

  for(map<int,C_SURFACE_FEMethod_ForSnes*>::iterator mit = C_SURFACE_FEMethod_ForSnes_msId_ptr.begin();
    mit!=C_SURFACE_FEMethod_ForSnes_msId_ptr.end();mit++) {
    delete mit->second;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::solve_coupled_problem(FieldInterface& mField,SNES *snes,const double da,const double fraction_treshold) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ierr = front_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);

  //create matrices
  Mat K;
  ierr = mField.MatCreateMPIAIJWithArrays("COUPLED_PROBLEM",&K); CHKERRQ(ierr);
  //create vectors
  Vec F;
  ierr = mField.VecCreateGhost("COUPLED_PROBLEM",Row,&F); CHKERRQ(ierr);
  Vec D;
  ierr = mField.VecCreateGhost("COUPLED_PROBLEM",Col,&D); CHKERRQ(ierr);

  if(material_FirelWall->operator[](FW_arc_lenhghat_definition)) {
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"arc length not initialised)");
  }

  ArcLengthCtx arc_ctx(mField,"COUPLED_PROBLEM");
  ArcLengthSnesCtx arc_snes_ctx(mField,"COUPLED_PROBLEM",&arc_ctx);
  ArcLengthMatShell arc_mat_ctx(mField,K,&arc_ctx,"COUPLED_PROBLEM");
  PCShellCtx pc_ctx(K,K,&arc_ctx);
  ArcLengthElemFEMethod arc_elem(mField,this,&arc_ctx);

  double gc;
  PetscBool flg;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
  }

  ErrorCode rval;
  Tag th_freez;
  const int def_order = 0;
  rval = mField.get_moab().tag_get_handle("FROZEN_NODE",1,MB_TYPE_INTEGER,th_freez,MB_TAG_CREAT|MB_TAG_SPARSE,&def_order); CHKERR_PETSC(rval);

  Range corners_edges,corners_nodes,unfreez_nodes;
  ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,corners_edges,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,corners_nodes,true); CHKERRQ(ierr);
  Range corners_edgesNodes;
  rval = mField.get_moab().get_connectivity(corners_edges,corners_edgesNodes,true); CHKERR_PETSC(rval);
  corners_nodes.insert(corners_edgesNodes.begin(),corners_edgesNodes.end());
  Range nodes_to_block;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-1);
  ierr = mField.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBVERTEX,nodes_to_block); CHKERRQ(ierr);
  corners_nodes.merge(nodes_to_block);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"freeze front nodes:\n");
  nb_un_freez_nodes = 0;
  bool freez_the_rest = false;
  for(
    map<EntityHandle,double>::iterator mit = map_ent_j.begin();
    mit!=map_ent_j.end();mit++) {
    double fraction = (max_j-mit->second)/max_j;
    //(fmin(max_j,gc)-mit->second)/fmin(max_j,gc);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "front node = %ld max_g = %6.4e g = %6.4e (%6.4e)",
      mit->first,max_j,mit->second,fraction); CHKERRQ(ierr);
    if(fraction > fraction_treshold || freez_the_rest) {
      ierr = PetscPrintf(PETSC_COMM_WORLD," freeze\n");
      corners_nodes.insert(mit->first);
      int freez = 1;
      EntityHandle node = mit->first;
      rval = mField.get_moab().tag_set_data(th_freez,&node,1,&freez); CHKERR_PETSC(rval);
    } else {
      int freez;
      EntityHandle node = mit->first;
      rval = mField.get_moab().tag_get_data(th_freez,&node,1,&freez); CHKERR_PETSC(rval);
      if(freez == 0) {
	unfreez_nodes.insert(mit->first);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");
	freez_the_rest = freeze_all_but_one;
      } else if(freez == 1) {
	ierr = PetscPrintf(PETSC_COMM_WORLD," unfreeze\n");
	int freez = 0;
	rval = mField.get_moab().tag_set_data(th_freez,&node,1,&freez); CHKERR_PETSC(rval);
	nb_un_freez_nodes++;
      } 
    }
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");
  //if(nb_un_freez_nodes) {
    //corners_nodes.insert(unfreez_nodes.begin(),unfreez_nodes.end());
  //}

  CubitDisplacementDirihletBC_Coupled myDirihletBC(mField,"COUPLED_PROBLEM",corners_nodes);
  ierr = myDirihletBC.Init(); CHKERRQ(ierr);

  const double YoungModulus = 1;
  const double PoissonRatio = 0.;
  NL_ElasticFEMethodCoupled MySpatialFE(
    mField,&myDirihletBC,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),&arc_ctx);
  NL_MaterialFEMethodCoupled MyMaterialFE(
    mField,&myDirihletBC,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),&arc_ctx);
  NL_MeshSmootherCoupled MyMeshSmoother(mField,&myDirihletBC);
  C_SURFACE_FEMethod_ForSnes MySurfaceConstrainsBodySurface(mField,&myDirihletBC,"LAMBDA_SURFACE");
  C_SURFACE_FEMethod_ForSnes MySurfaceConstrainsCrackSurface(mField,&myDirihletBC,"LAMBDA_CRACK_SURFACE");
  Snes_CTgc_CONSTANT_AREA_FEMethod MyCTgc(mField,*projFrontCtx,"COUPLED_PROBLEM","LAMBDA_CRACKFRONT_AREA");
  Snes_dCTgc_CONSTANT_AREA_FEMethod MydCTgc(mField,K,"LAMBDA_CRACKFRONT_AREA");
  TangentWithMeshSmoothingFrontConstrain_FEMethod MyTangentFrontContrain(mField,&MyMeshSmoother,"LAMBDA_CRACK_TANGENT_CONSTRAIN");

  ////******
  ierr = MySpatialFE.initCrackFrontData(mField); CHKERRQ(ierr);
  ierr = MyMaterialFE.initCrackFrontData(mField); CHKERRQ(ierr);
  ierr = MyMeshSmoother.initCrackFrontData(mField); CHKERRQ(ierr);
  MySurfaceConstrainsBodySurface.nonlinear = true;
  MySurfaceConstrainsCrackSurface.nonlinear = true;
  ////******

  map<int,C_SURFACE_FEMethod_ForSnes*> C_SURFACE_FEMethod_ForSnes_msId_ptr;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    C_SURFACE_FEMethod_ForSnes_msId_ptr[msId] = new C_SURFACE_FEMethod_ForSnes(mField,&myDirihletBC,ss.str());
    C_SURFACE_FEMethod_ForSnes_msId_ptr[msId]->nonlinear = true;
  }

  SnesCtx::loops_to_do_type& loops_to_do_Rhs = arc_snes_ctx.get_loops_to_do_Rhs();
  arc_snes_ctx.get_preProcess_to_do_Rhs().push_back(&MyCTgc);
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("CandCT_SURFACE_ELEM",&MySurfaceConstrainsBodySurface));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM",&MySurfaceConstrainsCrackSurface));
  for(map<int,C_SURFACE_FEMethod_ForSnes*>::iterator mit = C_SURFACE_FEMethod_ForSnes_msId_ptr.begin();
    mit!=C_SURFACE_FEMethod_ForSnes_msId_ptr.end();mit++) {
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << mit->first;
    loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type(ss2.str(),mit->second));
  }
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MESH_SMOOTHER",&MyMeshSmoother));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("C_TANGENT_ELEM",&MyTangentFrontContrain));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ELASTIC_COUPLED",&MySpatialFE));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MATERIAL_COUPLED",&MyMaterialFE));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&arc_elem));

  SnesCtx::loops_to_do_type& loops_to_do_Mat = arc_snes_ctx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("dCT_CRACKFRONT_AREA_ELEM",&MydCTgc));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("CandCT_SURFACE_ELEM",&MySurfaceConstrainsBodySurface));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM",&MySurfaceConstrainsCrackSurface));
  for(map<int,C_SURFACE_FEMethod_ForSnes*>::iterator mit = C_SURFACE_FEMethod_ForSnes_msId_ptr.begin();
    mit!=C_SURFACE_FEMethod_ForSnes_msId_ptr.end();mit++) {
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << mit->first;
    loops_to_do_Mat.push_back(SnesCtx::loop_pair_type(ss2.str(),mit->second));
  }
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MESH_SMOOTHER",&MyMeshSmoother));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("C_TANGENT_ELEM",&MyTangentFrontContrain));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ELASTIC_COUPLED",&MySpatialFE));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MATERIAL_COUPLED",&MyMaterialFE));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&arc_elem));

  PetscInt M,N;
  ierr = MatGetSize(K,&M,&N); CHKERRQ(ierr);
  PetscInt m,n;
  MatGetLocalSize(K,&m,&n);
  Mat ShellK;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,(void*)&arc_mat_ctx,&ShellK); CHKERRQ(ierr);
  ierr = MatShellSetOperation(ShellK,MATOP_MULT,(void(*)(void))arc_length_mult_shell); CHKERRQ(ierr);

  ierr = SNESSetApplicationContext(*snes,&arc_snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(*snes,F,SnesRhs,&arc_snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(*snes,ShellK,K,SnesMat,&arc_snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetConvergenceTest(*snes,MySnesConvernceTest,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
  PetscReal my_tol;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_tol",&my_tol,&flg); CHKERRQ(ierr);
  if(flg == PETSC_TRUE) {
    PetscReal atol,rtol,stol;
    PetscInt maxit,maxf;
    ierr = SNESGetTolerances(*snes,&atol,&rtol,&stol,&maxit,&maxf); CHKERRQ(ierr);
    atol = my_tol;
    rtol = atol*1e2;
    ierr = SNESSetTolerances(*snes,atol,rtol,stol,maxit,maxf); CHKERRQ(ierr);
  }

  KSP ksp;
  ierr = SNESGetKSP(*snes,&ksp); CHKERRQ(ierr);
  if(flg == PETSC_TRUE) {
    PetscReal rtol,atol,dtol;
    PetscInt maxits;
    ierr = KSPGetTolerances(ksp,&rtol,&atol,&dtol,&maxits); CHKERRQ(ierr);
    atol = my_tol*1e-2;
    rtol = atol*1e-2;
    ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxits); CHKERRQ(ierr);
  }
  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
  ierr = PCShellSetContext(pc,&pc_ctx); CHKERRQ(ierr);
  ierr = PCShellSetApply(pc,pc_apply_arc_length); CHKERRQ(ierr);
  ierr = PCShellSetSetUp(pc,pc_setup_arc_length); CHKERRQ(ierr);

  arc_ctx.get_FieldData() = *(MySpatialFE.t_val);
  ierr = mField.set_local_VecCreateGhost("COUPLED_PROBLEM",Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.get_problem("COUPLED_PROBLEM",&(arc_elem.problem_ptr)); CHKERRQ(ierr);
  ierr = arc_elem.set_dlambda_to_x(D,0); CHKERRQ(ierr);
  ierr = VecCopy(D,arc_ctx.x0); CHKERRQ(ierr);
  ierr = arc_elem.get_dlambda(D); CHKERRQ(ierr);
  ierr = arc_ctx.set_alpha_and_beta(1,0); CHKERRQ(ierr);
  MySpatialFE.set_snes_ctx(FieldInterface::SnesMethod::ctx_SNESNone);
  MyMaterialFE.set_snes_ctx(FieldInterface::SnesMethod::ctx_SNESNone);
  MySpatialFE.snes_x = D;
  MySpatialFE.snes_f = F;
  MyMaterialFE.snes_x = D;
  MyMaterialFE.snes_f = F;
  ierr = mField.loop_finite_elements("COUPLED_PROBLEM","ELASTIC_COUPLED",MySpatialFE);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("COUPLED_PROBLEM","MATERIAL_COUPLED",MyMaterialFE);  CHKERRQ(ierr);
  //set s
  arc_elem.snes_x = D;
  ierr = arc_elem.calulate_lambda_int(); CHKERRQ(ierr);
  ierr = arc_ctx.set_s(arc_elem.lambda_int+da/arc_elem.aRea0); CHKERRQ(ierr);

  ierr = SNESSolve(*snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(*snes,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  *(MySpatialFE.t_val) = arc_ctx.get_FieldData();

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("COUPLED_PROBLEM",Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = MySpatialFE.set_t_val(arc_ctx.get_FieldData()); CHKERRQ(ierr);
  aRea = arc_elem.aRea;
  lambda = arc_ctx.get_FieldData();

  ierr = mField.set_field(0,MBVERTEX,"SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.set_field(0,MBEDGE,"SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.set_field(0,MBTRI,"SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.set_field(0,MBTET,"SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.field_axpy(+1.,"SPATIAL_POSITION","SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.field_axpy(-1.,"MESH_NODE_POSITIONS","SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);

  Vec DISP;
  ierr = VecDuplicate(D,&DISP); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost(
    "COUPLED_PROBLEM","SPATIAL_POSITION","SPATIAL_DISPLACEMENT",Col,DISP,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecDot(DISP,arc_ctx.F_lambda,&energy); CHKERRQ(ierr);
  energy = 0.5*lambda*energy;

  SNESConvergedReason reason;
  ierr = SNESGetConvergedReason(*snes,&reason); CHKERRQ(ierr);

  int verb = 1;
  if(reason>=0) {

    //Save field on mesh
    if(verb>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Save SPATIAL_POSITION on tags\n"); CHKERRQ(ierr);
    }
    PostProcVertexMethod ent_method_spatial(mField.get_moab(),"SPATIAL_POSITION");
    ierr = mField.loop_dofs("COUPLED_PROBLEM","SPATIAL_POSITION",Col,ent_method_spatial); CHKERRQ(ierr);
    if(verb>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Save MESH_NODE_POSITIONS on tags\n"); CHKERRQ(ierr);
    }
    PostProcVertexMethod ent_method_material(mField.get_moab(),"MESH_NODE_POSITIONS");
    ierr = mField.loop_dofs("COUPLED_PROBLEM","MESH_NODE_POSITIONS",Col,ent_method_material); CHKERRQ(ierr);
    if(verb>2) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Save MATERIAL_RESIDUAL on tags\n"); CHKERRQ(ierr);
      PostProcVertexMethod ent_method_res_mat(mField.get_moab(),"MESH_NODE_POSITIONS",F,"MATERIAL_RESIDUAL");
      ierr = mField.loop_dofs("COUPLED_PROBLEM","MESH_NODE_POSITIONS",Col,ent_method_res_mat); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Save SPATIAL_RESIDUAL on tags\n"); CHKERRQ(ierr);
      PostProcVertexMethod ent_method_res_spat(mField.get_moab(),"SPATIAL_POSITION",F,"SPATIAL_RESIDUAL");
      ierr = mField.loop_dofs("COUPLED_PROBLEM","SPATIAL_POSITION",Col,ent_method_res_spat); CHKERRQ(ierr);
    }
    if(verb>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Save SPATIAL_DISPLACEMENT on tags\n"); CHKERRQ(ierr);
      PostProcVertexMethod ent_method_disp_spat(mField.get_moab(),"SPATIAL_POSITION",DISP,"SPATIAL_DISPLACEMENT");
      ierr = mField.loop_dofs("COUPLED_PROBLEM","SPATIAL_POSITION",Col,ent_method_disp_spat); CHKERRQ(ierr);
    }
 
  }

  ierr = VecDestroy(&DISP); CHKERRQ(ierr);
  ierr = MatDestroy(&ShellK); CHKERRQ(ierr);
  ierr = MatDestroy(&K); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);

  for(map<int,C_SURFACE_FEMethod_ForSnes*>::iterator mit = C_SURFACE_FEMethod_ForSnes_msId_ptr.begin();
    mit!=C_SURFACE_FEMethod_ForSnes_msId_ptr.end();mit++) {
    delete mit->second;
  }

  ierr = delete_front_projection_data(mField); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::calculate_material_forces(FieldInterface& mField,string problem,string fe) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  Range CornersEdges,corners_nodes;
  ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,CornersEdges,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,corners_nodes,true); CHKERRQ(ierr);
  ErrorCode rval;
  Range CornersEdgesNodes;
  rval = mField.get_moab().get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
  corners_nodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
  Range nodes_to_block;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-1);
  ierr = mField.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBVERTEX,nodes_to_block); CHKERRQ(ierr);
  corners_nodes.merge(nodes_to_block);
  CubitDisplacementDirihletBC_Coupled myDirihletBC(mField,problem,corners_nodes);
  ierr = myDirihletBC.Init(); CHKERRQ(ierr);

  const double YoungModulus = 1;
  const double PoissonRatio = 0.;
  NL_MaterialFEMethod MyMaterialFE(mField,&myDirihletBC,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  PetscBool flg = PETSC_TRUE;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_thermal_expansion",&MyMaterialFE.thermal_expansion,&flg); CHKERRQ(ierr);

  Vec F_Material;
  ierr = mField.VecCreateGhost(problem,Col,&F_Material); CHKERRQ(ierr);
  ierr = VecZeroEntries(F_Material); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_Material,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Material,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  MyMaterialFE.snes_f = F_Material;
  MyMaterialFE.set_snes_ctx(FieldInterface::SnesMethod::ctx_SNESSetFunction);
  ierr = mField.loop_finite_elements(problem,fe,MyMaterialFE);  CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_Material,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Material,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F_Material); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F_Material); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_material(mField.get_moab(),"MESH_NODE_POSITIONS",F_Material,"MATERIAL_FORCE");
  ierr = mField.loop_dofs(problem,"MESH_NODE_POSITIONS",Col,ent_method_material); CHKERRQ(ierr);

  double nrm2_material;
  ierr = VecNorm(F_Material,NORM_2,&nrm2_material);   CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"nrm2_F_Material = %6.4e\n",nrm2_material); CHKERRQ(ierr);

  //Fields
  ierr = mField.set_other_global_VecCreateGhost(
    problem,"MESH_NODE_POSITIONS","MATERIAL_FORCE",
    Row,F_Material,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"MATERIAL_FORCE",dof)) {
    //cerr << *dof << endl;
  //}

  //detroy matrices
  ierr = VecDestroy(&F_Material); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
ConfigurationalFractureMechanics::ArcLengthElemFEMethod::ArcLengthElemFEMethod(
  FieldInterface& _mField,ConfigurationalFractureMechanics *_conf_prob,ArcLengthCtx *_arc_ptr): 
    mField(_mField),conf_prob(_conf_prob),arc_ptr(_arc_ptr) {
    PetscErrorCode ierr;
    ErrorCode rval;

    //ghost dofs for diagonal value
    PetscInt ghosts[1] = { 0 };
    Interface &moab = mField.get_moab();
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm->rank() == 0) {
      ierr = VecCreateGhost(PETSC_COMM_WORLD,1,1,0,ghosts,&ghostDiag); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    } else {
      ierr = VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&ghostDiag); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }

    //get crack surcace
    ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,crackSurfacesFaces,true); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    Range level_tris;
    ierr = mField.get_entities_by_type_and_ref_level(*conf_prob->ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    crackSurfacesFaces = intersect(crackSurfacesFaces,level_tris);
    //get crack surface faces adjacenet to crack front
    Range crack_front_edges;
    ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    mField.get_moab().get_connectivity(crack_front_edges,crackFrontNodes,true);

    //get pointer to coupled problem
    ierr = mField.get_problem("COUPLED_PROBLEM",&problem_ptr); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    set<DofIdx> set_idx; //set of global dofs indexes for crack surface
    for(Range::iterator fit = crackSurfacesFaces.begin();fit!=crackSurfacesFaces.end();fit++) {
      const EntityHandle* conn; 
      int num_nodes; 
      rval = mField.get_moab().get_connectivity(*fit,conn,num_nodes,true); 
      if (MB_SUCCESS != rval) { 
	PetscAbortErrorHandler(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,__SDIR__,rval,PETSC_ERROR_INITIAL,
	  "can not get connectibility",PETSC_NULL);
	CHKERRABORT(PETSC_COMM_SELF,rval);
      }
      for(int nn = 0;nn<num_nodes; nn++) {
	for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problem_ptr,conn[nn],dit)) {
	  if(dit->get_name()!="MESH_NODE_POSITIONS") continue;
	  set_idx.insert(dit->get_petsc_gloabl_dof_idx());
	}
      }
    }

    //vector surfaceDofs keeps material positions for crack surface dofs
    ierr = PetscMalloc(set_idx.size()*sizeof(PetscInt),&isIdx); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    copy(set_idx.begin(),set_idx.end(),isIdx);
    ierr = ISCreateGeneral(PETSC_COMM_SELF,set_idx.size(),isIdx,PETSC_USE_POINTER,&isSurface); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF,set_idx.size(),&surfaceDofs); CHKERRABORT(PETSC_COMM_WORLD,ierr);

    ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",Row,&lambdaVec); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecSet(lambdaVec,1); CHKERRABORT(PETSC_COMM_WORLD,ierr);

}
ConfigurationalFractureMechanics::ArcLengthElemFEMethod::~ArcLengthElemFEMethod() {
  PetscErrorCode ierr;
  ierr = VecScatterDestroy(&surfaceScatter); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&ghostDiag); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&surfaceDofs); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = ISDestroy(&isSurface); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = PetscFree(isIdx); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&lambdaVec); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::calulate_area() {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  //scatter from snes_x to surfaceDofs
  ierr = VecScatterCreate(snes_x,isSurface,surfaceDofs,PETSC_NULL,&surfaceScatter); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecScatterBegin(surfaceScatter,snes_x,surfaceDofs,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(surfaceScatter,snes_x,surfaceDofs,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterDestroy(&surfaceScatter); CHKERRABORT(PETSC_COMM_SELF,ierr);

  //get crack front nodal positions form surfaceDofs
  double *array;
  ierr = VecGetArray(surfaceDofs,&array); CHKERRQ(ierr);  
  PetscInt size;
  ISGetSize(isSurface,&size);
  for(int ii = 0;ii<size;ii++) {
    problem_ptr->numered_dofs_rows.get<
      PetscGlobalIdx_mi_tag>().find(isIdx[ii])->get_FieldData() = array[ii];
  }
  ierr = VecRestoreArray(surfaceDofs,&array); CHKERRABORT(PETSC_COMM_SELF,ierr);

  //calulate crack surface area
  aRea = 0;
  aRea0 = 0;
  ublas::vector<double,ublas::bounded_array<double,6> > diffNTRI(6);
  ShapeDiffMBTRI(&diffNTRI.data()[0]);
  for(Range::iterator fit = crackSurfacesFaces.begin();fit!=crackSurfacesFaces.end();fit++) {
    const EntityHandle* conn; 
    int num_nodes; 
    rval = mField.get_moab().get_connectivity(*fit,conn,num_nodes,true); CHKERR_PETSC(rval);
    ublas::vector<double,ublas::bounded_array<double,9> > dofsX(9);
    ierr = mField.get_moab().get_coords(conn,num_nodes,&*dofsX.data().begin()); CHKERRQ(ierr);
    ublas::vector<double,ublas::bounded_array<double,3> > normal(3);
    ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&dofsX.data()[0],&normal.data()[0]); CHKERRQ(ierr);
    aRea0 += norm_2(normal)*0.25;
    for(int nn = 0;nn<num_nodes; nn++) {
      if(find(crackFrontNodes.begin(),crackFrontNodes.end(),conn[nn])!=crackFrontNodes.end()) {
      for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,"MESH_NODE_POSITIONS",conn[nn],dit)) {
	dofsX[nn*3+dit->get_dof_rank()] = dit->get_FieldData();
      }}
    }
    ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&dofsX.data()[0],&normal.data()[0]); CHKERRQ(ierr);
    //crack surface area is a half of crack top and bottom body surface
    aRea += norm_2(normal)*0.25;
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::calulate_lambda_int() {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = calulate_area(); CHKERRQ(ierr);
  lambda_int = arc_ptr->alpha*aRea/aRea0 + arc_ptr->dlambda*arc_ptr->beta*sqrt(arc_ptr->F_lambda2);
  //PetscPrintf(PETSC_COMM_WORLD,"\tsurface crack area = %6.4e lambda_int = %6.4e\n",aRea,lambda_int);
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::calulate_db() {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  if(arc_ptr->alpha!=0) {
    ierr = MatMultTranspose(conf_prob->projFrontCtx->C,lambdaVec,arc_ptr->db); CHKERRQ(ierr);
    ierr = VecScale(arc_ptr->db,1./aRea0); CHKERRQ(ierr);
    if(arc_ptr->alpha!=1) {
      ierr = VecScale(arc_ptr->db,arc_ptr->alpha); CHKERRQ(ierr);
    }
  } else {
    ierr = VecZeroEntries(arc_ptr->db); CHKERRQ(ierr);
  }   
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::set_dlambda_to_x(Vec x,double dlambda) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
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
    PetscPrintf(PETSC_COMM_WORLD,"\tload factor lambda_old,lambda_new/dlambda = %6.4e, %6.4e (%6.4e)\n",lambda_old, array[arc_ptr->get_petsc_local_dof_idx()], dlambda);
    ierr = VecRestoreArray(x,&array); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::get_dlambda(Vec x) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  //set dx
  ierr = VecCopy(x,arc_ptr->dx); CHKERRQ(ierr);
  ierr = VecAXPY(arc_ptr->dx,-1,arc_ptr->x0); CHKERRQ(ierr);
  //if LAMBDA dof is on this partition
  if(arc_ptr->get_petsc_local_dof_idx()!=-1) {
    double *array;
    ierr = VecGetArray(arc_ptr->dx,&array); CHKERRQ(ierr);
    arc_ptr->dlambda = array[arc_ptr->get_petsc_local_dof_idx()];
    array[arc_ptr->get_petsc_local_dof_idx()] = 0;
    ierr = VecRestoreArray(arc_ptr->dx,&array); CHKERRQ(ierr);
  }
  //brodcast dlambda
  int part = arc_ptr->get_part();
  MPI_Bcast(&(arc_ptr->dlambda),1,MPI_DOUBLE,part,PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,"\tload factor increment dlambda = %6.4e\n",arc_ptr->dlambda);
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::preProcess() {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  switch(snes_ctx) {
    case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = get_dlambda(snes_x); CHKERRQ(ierr);
	ierr = calulate_lambda_int(); CHKERRQ(ierr);
      }
      break;
    case ctx_SNESSetJacobian: 
	ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = calulate_db(); CHKERRQ(ierr);
      break;
    default:
      break;
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::operator()() {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  //get dlambda dof 
  switch(snes_ctx) {
    case ctx_SNESSetFunction: {
      //calulate residual for arc length row
      arc_ptr->res_lambda = lambda_int - arc_ptr->s;
      ierr = VecSetValue(snes_f,arc_ptr->get_petsc_gloabl_dof_idx(),arc_ptr->res_lambda,INSERT_VALUES); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_SELF,"\n");
      PetscPrintf(PETSC_COMM_SELF,"\t** Arc-Length residual:\n");
      PetscPrintf(PETSC_COMM_SELF,"\t  residual of arc-length control res_lambda = %6.4e crack area/f_lambda_int = %6.4e ( %6.4e )\n"
	,arc_ptr->res_lambda,aRea,aRea/aRea0);
    } break; 
    case ctx_SNESSetJacobian: {
      //calulate diagonal therm
      double diag = arc_ptr->beta*sqrt(arc_ptr->F_lambda2);
      ierr = VecSetValue(ghostDiag,0,diag,INSERT_VALUES); CHKERRQ(ierr);
      ierr = MatSetValue(*snes_B,arc_ptr->get_petsc_gloabl_dof_idx(),arc_ptr->get_petsc_gloabl_dof_idx(),1,INSERT_VALUES); CHKERRQ(ierr);
    } break;
    default:
    break;
  }	
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::postProcess() {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  switch(snes_ctx) {
    case ctx_SNESSetFunction: { 
	ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
	double res_nrm2[6];
	Vec res_nrm2_vec;
	if(pcomm->rank()==0) {
	  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1,6,6,res_nrm2,&res_nrm2_vec); CHKERRQ(ierr);
	} else {
	  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1,0,6,res_nrm2,&res_nrm2_vec); CHKERRQ(ierr);
	}
	//
	Range CrackFrontEdges;
	ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,CrackFrontEdges,true); CHKERRQ(ierr);
	Range CrackFrontNodes;
	ErrorCode rval;
	rval = mField.get_moab().get_connectivity(CrackFrontEdges,CrackFrontNodes,true); CHKERR_PETSC(rval);
	//
	double *array;
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = VecGetArray(snes_f,&array); CHKERRQ(ierr);
	ierr = VecZeroEntries(res_nrm2_vec); CHKERRQ(ierr);
	for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_LOCIDX_FOR_LOOP_(problem_ptr,dof)) {
	  if(dof->get_part()!=pcomm->rank()) continue;
	  double val = pow(array[dof->get_petsc_local_dof_idx()],2);
	  if(dof->get_name() == "SPATIAL_POSITION") {
	    ierr = VecSetValue(res_nrm2_vec,0,val,ADD_VALUES); CHKERRQ(ierr);
	  }
	  if(dof->get_name() == "MESH_NODE_POSITIONS") {
	    if(find(CrackFrontNodes.begin(),CrackFrontNodes.end(),dof->get_ent())!=CrackFrontNodes.end()) {
	      ierr = VecSetValue(res_nrm2_vec,1,val,ADD_VALUES); CHKERRQ(ierr);
	    } else {
	      ierr = VecSetValue(res_nrm2_vec,2,val,ADD_VALUES); CHKERRQ(ierr);
	    }
	  }
	  if(dof->get_name().compare(0,14,"LAMBDA_SURFACE") == 0) {
	    ierr = VecSetValue(res_nrm2_vec,3,val,ADD_VALUES); CHKERRQ(ierr);
	  }
	  if(dof->get_name() == "LAMBDA_CRACK_SURFACE") {
	    ierr = VecSetValue(res_nrm2_vec,4,val,ADD_VALUES); CHKERRQ(ierr);
	  }
	  if(dof->get_name() == "LAMBDA_CRACK_TANGENT_CONSTRAIN") {
	    ierr = VecSetValue(res_nrm2_vec,5,val,ADD_VALUES); CHKERRQ(ierr);
	  }
	}
	ierr = VecAssemblyBegin(res_nrm2_vec); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(res_nrm2_vec); CHKERRQ(ierr);
	if(pcomm->rank()==0) {
	  PetscPrintf(PETSC_COMM_WORLD,"\t** Equilibrium Residuals:\n");
	  PetscPrintf(PETSC_COMM_WORLD,"\t  equilibrium of physical forces res_spatial_nrm2 = %6.4e\n",sqrt(res_nrm2[0]));
	  PetscPrintf(PETSC_COMM_WORLD,"\t  equilibrium of material forces res_crack_front_nrm2 = %6.4e\n",sqrt(res_nrm2[1]));
	  PetscPrintf(PETSC_COMM_WORLD,"\t** Mesh Quality Residuals:\n");
	  PetscPrintf(PETSC_COMM_WORLD,"\t  residual of mesh smoother res_mesh_smoother_nrm2 = %6.4e\n",sqrt(res_nrm2[2]));
	  PetscPrintf(PETSC_COMM_WORLD,"\t  residual of surface constraint res_surface_constrain_nrm2 = %6.4e\n",sqrt(res_nrm2[3]));
	  PetscPrintf(PETSC_COMM_WORLD,"\t  residual of crack surface constraint res_crack_surface_constrain_nrm2 = %6.4e\n",sqrt(res_nrm2[4]));
	  PetscPrintf(PETSC_COMM_WORLD,"\t  residual of crack front mesh smoother tangential constraint res_crack_front_tangent_constrain_nrm2 = %6.4e\n",sqrt(res_nrm2[5]));
	  PetscPrintf(PETSC_COMM_WORLD,"\n");
	}
	ierr = VecRestoreArray(snes_f,&array); CHKERRQ(ierr);
	ierr = VecDestroy(&res_nrm2_vec); CHKERRQ(ierr);
    } break;
    case ctx_SNESSetJacobian: {
	ierr = VecAssemblyBegin(ghostDiag); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(ghostDiag); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(ghostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(ghostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	double *diag;
	ierr = VecGetArray(ghostDiag,&diag); CHKERRQ(ierr);
	arc_ptr->diag = *diag;
	ierr = VecRestoreArray(ghostDiag,&diag); CHKERRQ(ierr);
	//PetscPrintf(PETSC_COMM_WORLD,"\tdiag = %6.4e\n",arc_ptr->diag);
	/*ierr = MatAssemblyBegin(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
 	//Matrix View
	MatView(*snes_B,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
	std::string wait;
	std::cin >> wait;*/
    } break;
    default:
	break;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode main_spatial_solution(FieldInterface& mField,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  ierr = conf_prob.set_material_fire_wall(mField); CHKERRQ(ierr);

  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); 
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  EntityHandle meshset_level0;
  rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  ierr = conf_prob.spatial_problem_definition(mField); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_set_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = conf_prob.spatial_partition_problems(mField); CHKERRQ(ierr);

  //solve problem
  ierr = conf_prob.set_spatial_positions(mField); CHKERRQ(ierr);
  ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);

  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);  
  ierr = conf_prob.solve_spatial_problem(mField,&snes); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  PetscFunctionReturn(0);

}
PetscErrorCode main_material_forces(FieldInterface& mField,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr; 

  ierr = conf_prob.set_material_fire_wall(mField); CHKERRQ(ierr);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  ierr = conf_prob.spatial_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.constrains_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.material_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.constrains_crack_front_problem_definition(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MATERIAL",MBTET); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_set_bit("MATERIAL_MECHANICS",bit_level0); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_set_bit("CCT_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("C_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("C_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("CTC_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = conf_prob.spatial_partition_problems(mField); CHKERRQ(ierr);
  ierr = conf_prob.material_partition_problems(mField); CHKERRQ(ierr);
  ierr = conf_prob.constrains_partition_problems(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.crackfront_partition_problems(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);

  //caculate material forces
  ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);
  ierr = conf_prob.front_projection_data(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.surface_projection_data(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.calculate_material_forces(mField,"MATERIAL_MECHANICS","MATERIAL"); CHKERRQ(ierr);
  ierr = conf_prob.griffith_force_vector(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.project_force_vector(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.griffith_g(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.delete_surface_projection_data(mField); CHKERRQ(ierr);
  ierr = conf_prob.delete_front_projection_data(mField); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode main_arc_length_setup(FieldInterface& mField,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;
  //PetscBool flg = PETSC_TRUE;

  ierr = conf_prob.set_material_fire_wall(mField); CHKERRQ(ierr);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  ierr = conf_prob.spatial_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.constrains_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.material_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.mesh_smoothing_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.coupled_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.constrains_crack_front_problem_definition(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.arclength_problem_definition(mField); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC_COUPLED",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MATERIAL_COUPLED",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MESH_SMOOTHER",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MATERIAL",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_set_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("MATERIAL_MECHANICS",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("MESH_SMOOTHING_PROBLEM",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("COUPLED_PROBLEM",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("CCT_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("C_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("C_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("CTC_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = conf_prob.spatial_partition_problems(mField); CHKERRQ(ierr);
  ierr = conf_prob.material_partition_problems(mField); CHKERRQ(ierr);
  ierr = conf_prob.coupled_partition_problems(mField); CHKERRQ(ierr);
  ierr = conf_prob.constrains_partition_problems(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.crackfront_partition_problems(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode main_rescale_load_factor(FieldInterface& mField,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  double gc;
  PetscBool flg;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
  }

  const EntityHandle root_meshset = mField.get_moab().get_root_set();

  Tag th_t_val;
  rval = mField.get_moab().tag_get_handle("_LoadFactor_t_val",th_t_val); CHKERR_PETSC(rval);
  double *load_factor_ptr;
  rval = mField.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&load_factor_ptr); CHKERR_THROW(rval);
  double& load_factor = *load_factor_ptr;

  double max_j = conf_prob.max_j;
  double a = fabs(max_j)/pow(load_factor,2);
  double new_load_factor = copysign(sqrt(gc/a),load_factor);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\ncoefficient a = %6.4e\n",a); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"new load factor value = %6.4e\n\n",new_load_factor); CHKERRQ(ierr);
  load_factor = new_load_factor;
  SNES snes;
  //solve spatial problem
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);  
  ierr = conf_prob.solve_spatial_problem(mField,&snes,false); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);
  ierr = conf_prob.calculate_material_forces(mField,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
  ierr = conf_prob.front_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.surface_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.project_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.griffith_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.griffith_g(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode main_arc_length_solve(FieldInterface& mField,ConfigurationalFractureMechanics& conf_prob,bool face_splitting) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscBool flg = PETSC_TRUE;

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  //caculate material forces
  ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);

  double da_0 = 0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_da",&da_0,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_da (what is the crack area increment ?)");
  }
  double da = da_0;

  int def_nb_load_steps = 0;
  Tag th_nb_load_steps;
  rval = mField.get_moab().tag_get_handle("_NB_LOAD_STEPS",1,MB_TYPE_INTEGER,
    th_nb_load_steps,MB_TAG_CREAT|MB_TAG_SPARSE,&def_nb_load_steps); 
  int *ptr_nb_load_steps;
  rval = mField.get_moab().tag_get_by_ptr(th_nb_load_steps,&root_meshset,1,(const void**)&ptr_nb_load_steps); CHKERR_PETSC(rval);
  int &step = *ptr_nb_load_steps;
  if(step!=0) {
    step++;
  }

  PetscInt nb_load_steps = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_load_steps",&nb_load_steps,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_load_steps (what is the number of load_steps ?)");
  }

  PetscInt odd_face_split = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_odd_face_split",&odd_face_split,&flg); CHKERRQ(ierr);

  if(step == 0) {
    ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);
  }

  Tag th_t_val;
  rval = mField.get_moab().tag_get_handle("_LoadFactor_t_val",th_t_val); CHKERR_PETSC(rval);
  double *load_factor_ptr;
  rval = mField.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&load_factor_ptr); CHKERR_THROW(rval);
  double& load_factor = *load_factor_ptr;

  for(int aa = 0;step<nb_load_steps;step++,aa++) {

    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;
    ierr = PetscTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n** number of step = %D\n",step); CHKERRQ(ierr);

    if(aa == 0) {
      ierr = conf_prob.calculate_material_forces(mField,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
      ierr = conf_prob.front_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.surface_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.project_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.griffith_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.griffith_g(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.delete_surface_projection_data(mField); CHKERRQ(ierr);
      ierr = conf_prob.delete_front_projection_data(mField); CHKERRQ(ierr);
    }

    double gc;
    PetscBool flg;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
    }

    //calulate initial load factor
    if(step == 0) {
      ierr = main_rescale_load_factor(mField,conf_prob); CHKERRQ(ierr);
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD,"* da = %6.4e\n",da); CHKERRQ(ierr);
    SNES snes;
    ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
    
    Vec D0;
    ierr = mField.VecCreateGhost("COUPLED_PROBLEM",Col,&D0); CHKERRQ(ierr);
    ierr = mField.set_local_VecCreateGhost("COUPLED_PROBLEM",Col,D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    double load_factor0 = load_factor;
    
    int nb_sub_steps;
    ierr = PetscOptionsGetInt("","-my_nb_sub_steps",&nb_sub_steps,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      nb_sub_steps = 20;
    }

    bool first_step_converged = false;
    bool not_converged_state = false;
    bool line_searcher_set = false;
    double _da_ = 0;
    for(int ii = 0;ii<nb_sub_steps;ii++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n* number of substeps = %D\n\n",ii); CHKERRQ(ierr);
      if(ii == 0) {
	conf_prob.freeze_all_but_one = false;
	ierr = conf_prob.solve_coupled_problem(mField,&snes,(aa == 0) ? 0 : da); CHKERRQ(ierr);
      } else {
	ierr = conf_prob.solve_coupled_problem(mField,&snes,_da_); CHKERRQ(ierr);
      }
      int its;
      ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
      if(its == 0) break;
      SNESConvergedReason reason;
      ierr = SNESGetConvergedReason(snes,&reason); CHKERRQ(ierr);
      if(reason > 0) {
	conf_prob.freeze_all_but_one = false;
	if(aa > 0 && ii == 0) {
	  int its_d;
	  ierr = PetscOptionsGetInt("","-my_its_d",&its_d,&flg); CHKERRQ(ierr);
	  if(flg != PETSC_TRUE) {
	    its_d = 8;
	  }
	  double gamma = 0.5,reduction = 1;
	  reduction = pow((double)its_d/(double)(its+1),gamma);
	  const double max_da_reduction = 10;
	  if(reduction<1 || da < max_da_reduction*da_0) {
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n* change of da = %6.4e\n\n",reduction); CHKERRQ(ierr);
	    da = fmin(da*reduction,max_da_reduction*da_0);
	  }
	}
	_da_ = 0;
	first_step_converged = true;
	not_converged_state = false;
	ierr = mField.set_local_VecCreateGhost("COUPLED_PROBLEM",Col,D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	load_factor0 = load_factor;
	ierr = conf_prob.calculate_material_forces(mField,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
	ierr = conf_prob.front_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.surface_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.project_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.griffith_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.griffith_g(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.delete_surface_projection_data(mField); CHKERRQ(ierr);
	ierr = conf_prob.delete_front_projection_data(mField); CHKERRQ(ierr);
      } else {
	not_converged_state = true;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"* reset unknowns vector\n"); CHKERRQ(ierr);
	ierr = mField.set_global_VecCreateGhost("COUPLED_PROBLEM",Col,D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	load_factor = load_factor0;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"* failed to converge, recalculate spatial positions\n"); CHKERRQ(ierr);
	{
	  SNES snes_spatial;
	  ierr = SNESCreate(PETSC_COMM_WORLD,&snes_spatial); CHKERRQ(ierr);  
	  ierr = conf_prob.solve_spatial_problem(mField,&snes_spatial,false); CHKERRQ(ierr);
	  ierr = SNESDestroy(&snes_spatial); CHKERRQ(ierr);
	  ierr = conf_prob.calculate_material_forces(mField,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
	  ierr = conf_prob.front_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	  ierr = conf_prob.surface_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	  ierr = conf_prob.project_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	  ierr = conf_prob.griffith_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	  ierr = conf_prob.griffith_g(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	}
	if(!first_step_converged) {
	  if(da > 0) {
	    da = 0.5*da;
	    _da_ = da;
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"* failed to converge, set da = %6.4e ( 0.5 )\n",_da_); CHKERRQ(ierr);
	  } 
	} else {
	  if(da>0) {
	    if(not_converged_state) {
	      _da_+= 0.1*da;
	      ierr = PetscPrintf(PETSC_COMM_WORLD,"* failed to converge, set da = %6.4e\n",_da_); CHKERRQ(ierr);
	    }
	  }
	}
      }
      //just split and make animation going
      if(reason < 0 && odd_face_split != 0) {
	if(line_searcher_set) {
	  if(conf_prob.freeze_all_but_one) {
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"* use spatial solution only and split faces\n"); CHKERRQ(ierr);
	    break;
	  } else {
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"* freez all but one\n"); CHKERRQ(ierr);
	    conf_prob.freeze_all_but_one = true;
	  }
	}
      }
      //set line sercher L2 for not converged state
      if(reason < 0) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"* set L2 linesercher\n",_da_); CHKERRQ(ierr);
	//set line sercher L2 for not converged state
	ierr = SNESDestroy(&snes); CHKERRQ(ierr);
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
	ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
	SNESLineSearch linesearch;
	ierr = SNESGetLineSearch(snes,&linesearch); CHKERRQ(ierr);
	ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHL2); CHKERRQ(ierr);
	line_searcher_set = true;
      }
    }
    ierr = VecDestroy(&D0); CHKERRQ(ierr);
    ierr = SNESDestroy(&snes); CHKERRQ(ierr);

    ierr = conf_prob.set_coordinates_from_material_solution(mField); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "load_path: %4D Area %6.4e Lambda %6.4e Energy %6.4e\n",
      step,conf_prob.aRea,conf_prob.lambda,conf_prob.energy); CHKERRQ(ierr);
    if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      ierr = mField.problem_get_FE("COUPLED_PROBLEM","MATERIAL_COUPLED",out_meshset); CHKERRQ(ierr);
      ostringstream ss;
      ss << "out_load_step_" << step << ".vtk";
      rval = mField.get_moab().write_file(ss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);

      {
	Range SurfacesFaces;
	ierr = mField.get_Cubit_msId_entities_by_dimension(102,SideSet,2,SurfacesFaces,true); CHKERRQ(ierr);
	Range CrackSurfacesFaces;
	ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,CrackSurfacesFaces,true); CHKERRQ(ierr);
	for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
	  int msId = it->get_msId();
	  if((msId < 10200)||(msId >= 10300)) continue;
	  Range SurfacesFaces_msId;
	  ierr = mField.get_Cubit_msId_entities_by_dimension(msId,SideSet,2,SurfacesFaces_msId,true); CHKERRQ(ierr);
	  SurfacesFaces.insert(SurfacesFaces_msId.begin(),SurfacesFaces_msId.end());
	}

	Range level_tris;
	ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
	SurfacesFaces = intersect(SurfacesFaces,level_tris);
	CrackSurfacesFaces = intersect(CrackSurfacesFaces,level_tris);

	EntityHandle out_meshset;
	rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	rval = mField.get_moab().add_entities(out_meshset,CrackSurfacesFaces); CHKERR_PETSC(rval);
	ostringstream ss1;
	ss1 << "out_crack_surface_" << step << ".vtk";
	rval = mField.get_moab().write_file(ss1.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	rval = mField.get_moab().add_entities(out_meshset,SurfacesFaces); CHKERR_PETSC(rval);
	ostringstream ss2;
	ss2 << "out_surface_" << step << ".vtk";
	rval = mField.get_moab().write_file(ss2.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      }

      ostringstream ss3;
      ss3 << "restart_" << step << ".h5m";
      rval = mField.get_moab().write_file(ss3.str().c_str()); CHKERR_PETSC(rval);

      const MoFEMProblem *problem_ptr;
      ierr = mField.get_problem("COUPLED_PROBLEM",&problem_ptr); CHKERRQ(ierr);

      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|UnknownCubitName,it)) {
	if(it->get_Cubit_name() != "LoadPath") continue;

	Range nodes;
	rval = mField.get_moab().get_entities_by_type(it->meshset,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
	for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {

	  double coords[3];
	  rval = mField.get_moab().get_coords(&*nit,1,coords); CHKERR_PETSC(rval);
	  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problem_ptr,*nit,dof)) {
	    if(dof->get_name()!="SPATIAL_POSITION") continue;
	    ierr = PetscPrintf(PETSC_COMM_WORLD,
	      "load_path_disp ent %ld dim %d "
	      "coords ( %8.6f %8.6f %8.6f ) "
	      "val %6.4e Lambda %6.4e\n",
	      dof->get_ent(),dof->get_dof_rank(),
	      coords[0],coords[1],coords[2],
	      dof->get_FieldData()-coords[dof->get_dof_rank()],      
	      conf_prob.lambda); CHKERRQ(ierr);
	  }
	}

      }
    }

    ierr = PetscTime(&v2);CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Step Time Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);

    if(odd_face_split != 0) {
      
      if( aa>0 && aa % odd_face_split == 0 ) {

	{ //cat mesh
  
	  FaceSplittingTools face_splitting(mField);
	  ierr = main_refine_and_meshcat(mField,face_splitting,false); CHKERRQ(ierr);
	  ierr = face_splitting.cleanMeshsets(); CHKERRQ(ierr);

	}

	//project and set coords
	conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
	conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_material_positions) = 0;
	ierr = conf_prob.set_spatial_positions(mField); CHKERRQ(ierr);
	ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);

	//find faces for split
	FaceSplittingTools face_splitting(mField);
	ierr = main_select_faces_for_splitting(mField,face_splitting,0); CHKERRQ(ierr);
	//do splittig
	ierr = main_split_faces_and_update_field_and_elements(mField,face_splitting,0); CHKERRQ(ierr);
	
	//rebuild fields, finite elementa and problems
	ierr = main_face_splitting_restart(mField,conf_prob); CHKERRQ(ierr);

	//face spliting job done
	ierr = face_splitting.cleanMeshsets(); CHKERRQ(ierr);

	//solve spatial problem and calulate griffith forces

	//project and set coords
	conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
	conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_material_positions) = 0;
	ierr = conf_prob.set_spatial_positions(mField); CHKERRQ(ierr);
	ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);
 
	/*//debuging
	ierr = mField.partition_check_matrix_fill_in("MESH_SMOOTHING_PROBLEM",1); CHKERRQ(ierr);
	ierr = mField.partition_check_matrix_fill_in("ELASTIC_MECHANICS",1); CHKERRQ(ierr);
	ierr = mField.partition_check_matrix_fill_in("COUPLED_PROBLEM",1); CHKERRQ(ierr);*/

	SNES snes;
	//solve mesh smoothing
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);  
	ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
	SNESLineSearch linesearch;
	ierr = SNESGetLineSearch(snes,&linesearch); CHKERRQ(ierr);
	ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHL2); CHKERRQ(ierr);
	Vec D_tmp_mesh_positions;
	ierr = mField.VecCreateGhost("MESH_SMOOTHING_PROBLEM",Col,&D_tmp_mesh_positions); CHKERRQ(ierr);
	ierr = mField.set_local_VecCreateGhost("MESH_SMOOTHING_PROBLEM",Col,D_tmp_mesh_positions,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = conf_prob.front_projection_data(mField,"MESH_SMOOTHING_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.surface_projection_data(mField,"MESH_SMOOTHING_PROBLEM"); CHKERRQ(ierr);
	Range tets;
	ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);

	struct CheckForNegatieVolume {
	  static PetscErrorCode F(FieldInterface &mField,const Range &tets,bool &flg) {
	    PetscFunctionBegin;
	    ErrorCode rval;
	    PetscErrorCode ierr;
	    double diffNTET[12],coords[12],V;
	    ierr = ShapeDiffMBTET(diffNTET); CHKERRQ(ierr);
	    Range::iterator tit = tets.begin();
	    for(;tit!=tets.end();tit++) {
	      int num_nodes;
	      const EntityHandle *conn;
	      rval = mField.get_moab().get_connectivity(*tit,conn,num_nodes,true); CHKERR_PETSC(rval);
	      ierr = mField.get_FielData("MESH_NODE_POSITIONS",conn,num_nodes,coords); CHKERRQ(ierr);
	      V = Shape_intVolumeMBTET(diffNTET,coords); 
	      if(V<=0) break;	  
	    }
	    if(tit==tets.end()) {
	      flg = true;
	    } else {
	      flg = false;
	    }
	    PetscFunctionReturn(0);
	  }
	};

	int nb_sub_steps = 1;
	int nn;
	do { 

	  nn = 1;
	  for(;nn<=nb_sub_steps;nn++) {

	    ierr = PetscPrintf(PETSC_COMM_WORLD,"Mesh projection substep = %d out of %d\n",nn,nb_sub_steps); CHKERRQ(ierr);
	    double alpha = (double)nn/(double)nb_sub_steps;
	    ierr = face_splitting.calculateDistanceCrackFrontNodesFromCrackSurface(alpha); CHKERRQ(ierr);
	    //project nodes on crack surface
	    ierr = conf_prob.project_form_th_projection_tag(mField,"MESH_SMOOTHING_PROBLEM"); CHKERRQ(ierr);
	    bool flg;
	    ierr = CheckForNegatieVolume::F(mField,tets,flg); CHKERRQ(ierr);
	    if(!flg) {
	      ierr = mField.set_global_VecCreateGhost(
		"MESH_SMOOTHING_PROBLEM",
		Col,D_tmp_mesh_positions,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	      nb_sub_steps++;
	      break;
	    } else {

	      ierr = conf_prob.solve_mesh_smooting_problem(mField,&snes);  CHKERRQ(ierr);
	      SNESConvergedReason reason;
	      ierr = SNESGetConvergedReason(snes,&reason); CHKERRQ(ierr);
	      ierr = CheckForNegatieVolume::F(mField,tets,flg); CHKERRQ(ierr);
	      if(reason < 0 || !flg) {
		ierr = mField.set_global_VecCreateGhost(
		  "MESH_SMOOTHING_PROBLEM",
		  Col,D_tmp_mesh_positions,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
		nb_sub_steps++;
		break;
	      }
	      ierr = mField.set_local_VecCreateGhost(
		"MESH_SMOOTHING_PROBLEM",
		Col,D_tmp_mesh_positions,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	    }
	  }

	  if(nb_sub_steps >= 10) break;

	} while(nn-1 != nb_sub_steps);
	ierr = conf_prob.set_coordinates_from_material_solution(mField); CHKERRQ(ierr);
	ierr = VecDestroy(&D_tmp_mesh_positions); CHKERRQ(ierr);
	ierr = SNESDestroy(&snes); CHKERRQ(ierr);

	ierr = conf_prob.delete_surface_projection_data(mField); CHKERRQ(ierr);
	ierr = conf_prob.delete_front_projection_data(mField); CHKERRQ(ierr);

	//save on mesh spatial positions
	conf_prob.material_FirelWall->
	  operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
	ierr = conf_prob.set_spatial_positions(mField); CHKERRQ(ierr);


	//solve spatial problem
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);  
	ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
	ierr = conf_prob.solve_spatial_problem(mField,&snes,false); CHKERRQ(ierr);
	ierr = SNESDestroy(&snes); CHKERRQ(ierr);

	//calulate Griffth forces
	ierr = conf_prob.calculate_material_forces(mField,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
	ierr = conf_prob.front_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.surface_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.project_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.griffith_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.griffith_g(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.delete_surface_projection_data(mField); CHKERRQ(ierr);
	ierr = conf_prob.delete_front_projection_data(mField); CHKERRQ(ierr);

      }

    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode main_face_splitting_restart(FieldInterface& mField,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  Tag th_griffith_force;
  rval = mField.get_moab().tag_get_handle("GRIFFITH_FORCE",th_griffith_force); CHKERR_PETSC(rval);
  rval = mField.get_moab().tag_delete(th_griffith_force); CHKERR_PETSC(rval);
  Tag th_freez;
  rval = mField.get_moab().tag_get_handle("FROZEN_NODE",th_freez); CHKERR_PETSC(rval);
  rval = mField.get_moab().tag_delete(th_freez); CHKERR_PETSC(rval);

  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_spatial_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_material_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_coupled_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_constrains_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_constrains_crack_front_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_material_positions) = 0;

  //ierr = mField.rebuild_database(0); CHKERRQ(ierr);
  ierr = main_arc_length_setup(mField,conf_prob); CHKERRQ(ierr);
  ierr = mField.check_number_of_ents_in_ents_field(); CHKERRQ(ierr);
  ierr = mField.check_number_of_ents_in_ents_finite_element(); CHKERRQ(ierr);

  //ierr = mField.list_ent_by_field_name("SPATIAL_POSITION"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


