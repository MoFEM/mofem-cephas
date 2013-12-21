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

using namespace MoFEM;

phisical_equation_volume eq_solid = hooke; /*stvenant_kirchhoff;*/

struct materialDirihletBC;

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
    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(this,RowGlobSpatial,ColGlobSpatial,DirihletBC); CHKERRQ(ierr);
    switch(snes_ctx) {
      case ctx_SNESNone: {
	ierr = CalulateSpatialKFext(PETSC_NULL,arc_ptr->F_lambda,-1.); CHKERRQ(ierr);
      } break;
      case ctx_SNESSetFunction: { 
	ierr = CalculateSpatialFint(snes_f); CHKERRQ(ierr);
	ierr = CalulateSpatialKFext(PETSC_NULL,arc_ptr->F_lambda,-1.); CHKERRQ(ierr);
      } break;
      case ctx_SNESSetJacobian: {
	ierr = CalculateSpatialTangent(*snes_B); CHKERRQ(ierr);
	ierr = GetIndicesRow(RowGlobMaterial,material_field_name); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndiciesRow(this,RowGlobMaterial,DirihletBC); CHKERRQ(ierr);
	ierr = AssembleSpatialCoupledTangent(*snes_B); CHKERRQ(ierr);
	ierr = CalulateSpatialKFext(*snes_B,PETSC_NULL,-arc_ptr->get_FieldData()); CHKERRQ(ierr);

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
	    //ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,snes_f); CHKERRQ(ierr);
	  } break;
	  case ctx_SNESSetJacobian: {
	    ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	    ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	    //ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,*snes_B); CHKERRQ(ierr);
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
    ierr = GetData(dofs_x_edge_data,dofs_x_edge,
	dofs_x_face_data,dofs_x_face,
	dofs_x_volume,dofs_x,
	spatial_field_name); CHKERRQ(ierr);
    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(this,RowGlobMaterial,ColGlobMaterial,DirihletBC); CHKERRQ(ierr);
    switch(snes_ctx) {
	case ctx_SNESNone: {
	  ierr = CaluclateMaterialKFext(PETSC_NULL,arc_ptr->F_lambda,-1); CHKERRQ(ierr);
	}  break;
	case ctx_SNESSetFunction: { 
	  ierr = CalculateMaterialFint(snes_f); CHKERRQ(ierr);
	  ierr = CaluclateMaterialKFext(PETSC_NULL,arc_ptr->F_lambda,-1); CHKERRQ(ierr);
	} break;
	case ctx_SNESSetJacobian: {
	  ierr = GetTangent(); CHKERRQ(ierr);
	  ierr = AssembleMaterialTangent(*snes_B); CHKERRQ(ierr);
	  ierr = CaluclateMaterialKFext(*snes_B,PETSC_NULL,-arc_ptr->get_FieldData()); CHKERRQ(ierr);
	  ierr = GetIndicesRow(RowGlobSpatial,spatial_field_name); CHKERRQ(ierr);
	  ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndiciesRow(this,RowGlobSpatial,DirihletBC); CHKERRQ(ierr);
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
      case ctx_SNESNone: {
	ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	//F_lambda2
	ierr = VecDot(arc_ptr->F_lambda,arc_ptr->F_lambda,&arc_ptr->F_lambda2); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tFlambda2 = %6.4e\n",arc_ptr->F_lambda2);
      } break;
      case ctx_SNESSetFunction:
      case ctx_SNESSetJacobian: {
	switch(snes_ctx) {
	  case ctx_SNESSetFunction: { 
	    ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	    //F_lambda2
	    ierr = VecDot(arc_ptr->F_lambda,arc_ptr->F_lambda,&arc_ptr->F_lambda2); CHKERRQ(ierr);
	    PetscPrintf(PETSC_COMM_WORLD,"\tFlambda2 = %6.4e\n",arc_ptr->F_lambda2);
	    ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	    //Add external forces to RHS
	    ierr = VecAXPY(snes_f,-arc_ptr->get_FieldData(),arc_ptr->F_lambda); CHKERRQ(ierr);
	    PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",arc_ptr->get_FieldData());  
	    ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,snes_f); CHKERRQ(ierr);
	    ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	  } break;
	  case ctx_SNESSetJacobian: {
	    ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	    ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,*snes_B); CHKERRQ(ierr);
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

struct NL_MeshSmootherCoupled: public FEMethod_DriverComplexForLazy_MeshSmoothing {
  
    NL_MeshSmootherCoupled(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,int _verbose = 0):
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
      FEMethod_DriverComplexForLazy_MeshSmoothing(_mField,_dirihlet_bc_method_ptr) {
	set_qual_ver(0);
      }
  
  };

struct materialDirihletBC: public BaseDirihletBC {

  Interface& moab;
  Range &CornersNodes;
  string field_name;
  materialDirihletBC(Interface &_moab,Range& _CornerNodes): moab(_moab),CornersNodes(_CornerNodes),field_name("MESH_NODE_POSITIONS") {}

  PetscErrorCode SetDirihletBC_to_ElementIndicies(
      FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlob,vector<vector<DofIdx> > &ColGlob,vector<DofIdx>& DirihletBC) {
      PetscFunctionBegin;
      //Dirihlet form SideSet1
      DirihletBC.resize(0);
      Range::iterator siit1 = CornersNodes.begin();
      for(;siit1!=CornersNodes.end();siit1++) {
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
	  for(;riit!=hi_riit;riit++) {
	    if(riit->get_name()!=field_name) continue;
	    // all fixed
	    DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
	    for(unsigned int rr = 0;rr<RowGlob.size();rr++) {
	      vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
	      if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
	    }
	  }
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator ciit = fe_method_ptr->col_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_ciit = fe_method_ptr->col_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
	  for(;ciit!=hi_ciit;ciit++) {
	    if(ciit->get_name()!=field_name) continue;
	    for(unsigned int cc = 0;cc<ColGlob.size();cc++) {
	      vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),ciit->get_petsc_gloabl_dof_idx());
	      if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
	    }
	  }
      }
      PetscFunctionReturn(0);
  }

  PetscErrorCode SetDirihletBC_to_ElementIndiciesFace(
    FieldInterface::FEMethod *fe_method_ptr,vector<DofIdx>& DirihletBC,
    vector<DofIdx> &FaceNodeIndices,
    vector<vector<DofIdx> > &FaceEdgeIndices,
    vector<DofIdx> &FaceIndices) {
      PetscFunctionBegin;
      vector<DofIdx>::iterator dit = DirihletBC.begin();
      for(;dit!=DirihletBC.end();dit++) {
	vector<DofIdx>::iterator it = find(FaceNodeIndices.begin(),FaceNodeIndices.end(),*dit);
	if(it!=FaceNodeIndices.end()) *it = -1; // of idx is set -1 row is not assembled
	if(!FaceEdgeIndices.empty()) {
	for(int ee = 0;ee<3;ee++) {
	  it = find(FaceEdgeIndices[ee].begin(),FaceEdgeIndices[ee].end(),*dit);
	  if(it!=FaceEdgeIndices[ee].end()) *it = -1; // of idx is set -1 row is not assembled
	}}
	it = find(FaceIndices.begin(),FaceIndices.end(),*dit);
	if(it!=FaceIndices.end()) *it = -1; // of idx is set -1 row is not assembled
      }
      PetscFunctionReturn(0);
  }

  PetscErrorCode SetDirihletBC_to_MatrixDiagonal(
    FieldInterface::FEMethod *fe_method_ptr,Mat Aij) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_LOCIDX_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
      if(find(CornersNodes.begin(),CornersNodes.end(),dit->get_ent()) == CornersNodes.end()) continue;
      ierr = MatSetValue(Aij,dit->get_petsc_gloabl_dof_idx(),dit->get_petsc_gloabl_dof_idx(),1.,INSERT_VALUES); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode SetDirihletBC_to_RHS(FieldInterface::FEMethod *fe_method_ptr,Vec F) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_LOCIDX_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
      if(find(CornersNodes.begin(),CornersNodes.end(),dit->get_ent()) == CornersNodes.end()) continue;
      ierr = VecSetValue(F,dit->get_petsc_gloabl_dof_idx(),0.,INSERT_VALUES); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

};

PetscErrorCode ConfigurationalFractureMechanics::CubitDisplacementDirihletBC_Coupled::SetDirihletBC_to_ElementIndiciesRow(
      FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirihletBC) {
  PetscFunctionBegin;
  field_name = "SPATIAL_POSITION";
  ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndiciesRow(fe_method_ptr,RowGlobDofs,DirihletBC); CHKERRQ(ierr);
  field_name = "MESH_NODE_POSITIONS";
  ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndiciesRow(fe_method_ptr,RowGlobDofs,DirihletBC); CHKERRQ(ierr);
  Range::iterator siit1 = CornersNodes.begin();
  for(;siit1!=CornersNodes.end();siit1++) {
    for(_IT_GET_FEROW_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(fe_method_ptr,"MESH_NODE_POSITIONS",*siit1,riit)) {
      // all fixed
      DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
      for(unsigned int rr = 0;rr<RowGlobDofs.size();rr++) {
	vector<DofIdx>::iterator it = find(RowGlobDofs[rr].begin(),RowGlobDofs[rr].end(),riit->get_petsc_gloabl_dof_idx());
	if( it!=RowGlobDofs[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
      }
    }
  }
  PetscFunctionReturn(0);
}
  
PetscErrorCode ConfigurationalFractureMechanics::CubitDisplacementDirihletBC_Coupled::SetDirihletBC_to_ElementIndiciesCol(
  FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC) {
  PetscFunctionBegin;
  CubitDisplacementDirihletBC::field_name = "SPATIAL_POSITION";
  ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndiciesCol(fe_method_ptr,ColGlobDofs,DirihletBC); CHKERRQ(ierr);
  CubitDisplacementDirihletBC::field_name = "MESH_NODE_POSITIONS";
  ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndiciesCol(fe_method_ptr,ColGlobDofs,DirihletBC); CHKERRQ(ierr);
  Range::iterator siit1 = CornersNodes.begin();
  for(;siit1!=CornersNodes.end();siit1++) {
    for(_IT_GET_FECOL_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(fe_method_ptr,"MESH_NODE_POSITIONS",*siit1,ciit)) {
      for(unsigned int cc = 0;cc<ColGlobDofs.size();cc++) {
	vector<DofIdx>::iterator it = find(ColGlobDofs[cc].begin(),ColGlobDofs[cc].end(),ciit->get_petsc_gloabl_dof_idx());
	if( it!=ColGlobDofs[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::CubitDisplacementDirihletBC_Coupled::SetDirihletBC_to_MatrixDiagonal(FieldInterface::FEMethod *fe_method_ptr,Mat Aij) {
  PetscFunctionBegin;
  ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_MatrixDiagonal(fe_method_ptr,Aij); CHKERRQ(ierr);
  PetscErrorCode ierr;
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(fe_method_ptr->problem_ptr,"MESH_NODE_POSITIONS",dit)) {
    if(dit->get_name()!="MESH_NODE_POSITIONS") {
      SETERRQ(PETSC_COMM_SELF,1,"inconsistency");
    }
    if(find(CornersNodes.begin(),CornersNodes.end(),dit->get_ent()) == CornersNodes.end()) continue;
    ierr = MatSetValue(Aij,dit->get_petsc_gloabl_dof_idx(),dit->get_petsc_gloabl_dof_idx(),1.,INSERT_VALUES); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::CubitDisplacementDirihletBC_Coupled::SetDirihletBC_to_RHS(FieldInterface::FEMethod *fe_method_ptr,Vec F) {
  PetscFunctionBegin;
  ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_RHS(fe_method_ptr,F); CHKERRQ(ierr);
  for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_LOCIDX_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
    if(find(CornersNodes.begin(),CornersNodes.end(),dit->get_ent()) == CornersNodes.end()) continue;
    ierr = VecSetValue(F,dit->get_petsc_gloabl_dof_idx(),0.,INSERT_VALUES); CHKERRQ(ierr);
  }
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

PetscErrorCode ConfigurationalFractureMechanics::spatial_problem_definition(FieldInterface& mField) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_spatial_problem_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_spatial_problem_definition);

  PetscErrorCode ierr;

  //Fields
  ierr = mField.add_field("SPATIAL_POSITION",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC",MF_ZERO); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);
  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  PetscInt order;
  PetscBool flg = PETSC_TRUE;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }
  //set app. order
  ierr = mField.set_field_order(0,MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);
  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);


  /*PetscInt nb_ref_levels;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_ref",&nb_ref_levels,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    nb_ref_levels = 0;
  }

  ErrorCode rval;
  EntityHandle meshset_level;
  rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_level); CHKERR_PETSC(rval);

  for(int ll = 2;ll<nb_ref_levels+2;ll++) {

    rval = mField.get_moab().clear_meshset(&meshset_level,1); CHKERR(rval);

    BitRefLevel bit_level;
    bit_level.set(ll);
    for(int lll = ll+1;lll<nb_ref_levels+2;lll++) {
      bit_level.set(lll);
    }

    ierr = mField.refine_get_ents(bit_level,bit_level,meshset_level); CHKERRQ(ierr);
  
    int ref_order = order > ll ? order : ll;
    ref_order = ref_order > 5 ? 5 : ref_order;

    ierr = mField.set_field_order(meshset_level,MBTET,"SPATIAL_POSITION",ref_order,2); CHKERRQ(ierr);
    ierr = mField.set_field_order(meshset_level,MBTRI,"SPATIAL_POSITION",ref_order,0); CHKERRQ(ierr);
    ierr = mField.set_field_order(meshset_level,MBEDGE,"SPATIAL_POSITION",ref_order,0); CHKERRQ(ierr);

  }
  rval = mField.get_moab().delete_entities(&meshset_level,1); CHKERR(rval);

  Tag th_order;
  const int def_order = -1;
  rval = mField.get_moab().tag_get_handle("ORDER",1,MB_TYPE_INTEGER,th_order,MB_TAG_CREAT|MB_TAG_SPARSE,&def_order); CHKERR_THROW(rval);
  for(_IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(mField,"SPATIAL_POSITION",dit)) {
    EntityHandle ent = dit->get_ent();
    int order = dit->get_max_order();
    rval = mField.get_moab().tag_set_data(th_order,&ent,1,&order); CHKERR(rval);
  }*/

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

  //define problems
  ierr = mField.add_problem("MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.add_problem("MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS"); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("MATERIAL_MECHANICS","MATERIAL"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS","MATERIAL"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS","CandCT_SURFACE_ELEM"); CHKERRQ(ierr);

  bool cs = true;
  if(cs) {
    ierr = mField.modify_problem_add_finite_element("MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS","CandCT_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
  }

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

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
  //
  ierr = mField.modify_finite_element_add_field_row("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("MATERIAL_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MATERIAL_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  //
  ierr = mField.modify_finite_element_add_field_row("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("COUPLED_PROBLEM"); CHKERRQ(ierr);
  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","ELASTIC_COUPLED"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","MESH_SMOOTHER"); CHKERRQ(ierr);

  ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","CandCT_SURFACE_ELEM"); CHKERRQ(ierr);
  bool cs = true;
  if(cs) {
    ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","CandCT_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
  }

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);


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

  //Field for ArcLength
  ierr = mField.add_field("X0_MATERIAL_POSITION",H1,3); CHKERRQ(ierr);

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
  //CRACK
  if(cs) {
    ierr = mField.add_field("LAMBDA_CRACK_SURFACE",H1,1,MF_ZERO); CHKERRQ(ierr);
  }

  //FE
  ierr = mField.add_finite_element("C_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CTC_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CandCT_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);

  //CRACK
  if(cs) {
    ierr = mField.add_finite_element("C_CRACK_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.add_finite_element("CTC_CRACK_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.add_finite_element("CandCT_CRACK_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
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

  }

  //define problems
  ierr = mField.add_problem("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.add_problem("C_ALL_MATRIX"); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("C_ALL_MATRIX","C_SURFACE_ELEM"); CHKERRQ(ierr);
  //CRACK
  if(cs) {
    ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("C_ALL_MATRIX","C_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
  }


  {
    //
    Range CornersEdges,CornersNodes,SurfacesFaces;
    ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,CornersEdges,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,CornersNodes,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(102,SideSet,2,SurfacesFaces,true); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 100 = %d\n",CornersEdges.size()); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of NodeSet 101 = %d\n",CornersNodes.size()); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 102 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
    Range CrackSurfacesFaces,CrackCornersEdges;
    ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,CrackSurfacesFaces,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,CrackCornersEdges,true); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack SideSet 200 = %d\n",CrackSurfacesFaces.size()); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack SideSet 201 = %d\n",CrackCornersEdges.size()); CHKERRQ(ierr);

    Interface& moab = mField.get_moab();

    Range level_tris;
    ierr = mField.refine_get_ents(*ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
    SurfacesFaces = intersect(SurfacesFaces,level_tris);

    ierr = mField.seed_finite_elements(SurfacesFaces); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"C_SURFACE_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"CTC_SURFACE_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"CandCT_SURFACE_ELEM"); CHKERRQ(ierr);

    if(cs) {
      CrackSurfacesFaces = intersect(CrackSurfacesFaces,level_tris);
      ierr = mField.seed_finite_elements(CrackSurfacesFaces); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TRIs(CrackSurfacesFaces,"C_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TRIs(CrackSurfacesFaces,"CTC_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TRIs(CrackSurfacesFaces,"CandCT_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
    }

    Range CrackCornersEdgesNodes;
    rval = moab.get_connectivity(CrackCornersEdges,CrackCornersEdgesNodes,true); CHKERR_PETSC(rval);
    Range CornersEdgesNodes;
    rval = moab.get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
    {
      CornersNodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
      rval = moab.create_meshset(MESHSET_SET,cornersNodesMeshset); CHKERR_PETSC(rval);	
      rval = moab.add_entities(cornersNodesMeshset,CornersNodes); CHKERR_PETSC(rval);
    }
    Range SurfacesNodes;
    {
      rval = moab.get_connectivity(SurfacesFaces,SurfacesNodes,true); CHKERR_PETSC(rval);
      SurfacesNodes = subtract(SurfacesNodes,CornersNodes);
      SurfacesNodes = subtract(SurfacesNodes,CornersEdgesNodes);
      rval = moab.create_meshset(MESHSET_SET,surfacesFacesNodesMeshset); CHKERR_PETSC(rval);	
      rval = moab.add_entities(surfacesFacesNodesMeshset,SurfacesNodes); CHKERR_PETSC(rval);
    }
    //CRCAK
    if(cs) {
      {
	Range CrackSurfacesNodes;
	rval = moab.get_connectivity(CrackSurfacesFaces,CrackSurfacesNodes,true); CHKERR_PETSC(rval);
	CrackSurfacesNodes = subtract(CrackSurfacesNodes,CornersNodes);
	CrackSurfacesNodes = subtract(CrackSurfacesNodes,CornersEdgesNodes);
	CrackSurfacesNodes = subtract(CrackSurfacesNodes,CrackCornersEdgesNodes);
	rval = moab.create_meshset(MESHSET_SET,crackSurfacesFacesNodesMeshset); CHKERR_PETSC(rval);	
	rval = moab.add_entities(crackSurfacesFacesNodesMeshset,CrackSurfacesNodes); CHKERR_PETSC(rval);
      }
    }
  }

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_VERTICEs(surfacesFacesNodesMeshset,"LAMBDA_SURFACE"); CHKERRQ(ierr);
  //NOTE: always order should be 1
  ierr = mField.set_field_order(surfacesFacesNodesMeshset,MBVERTEX,"LAMBDA_SURFACE",1); CHKERRQ(ierr);
  //CRACK
  if(cs) {
    ierr = mField.add_ents_to_field_by_VERTICEs(crackSurfacesFacesNodesMeshset,"LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = mField.set_field_order(crackSurfacesFacesNodesMeshset,MBVERTEX,"LAMBDA_CRACK_SURFACE",1); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::constrains_crack_front_problem_definition(FieldInterface& mField,string problem) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_constrains_crack_front_problem_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_constrains_crack_front_problem_definition);

  PetscErrorCode ierr;
  ErrorCode rval;

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
  ierr = mField.add_problem("C_CRACKFRONT_MATRIX"); CHKERRQ(ierr);
  ierr = mField.add_problem("CTC_CRACKFRONT_MATRIX"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("C_CRACKFRONT_MATRIX","C_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("CTC_CRACKFRONT_MATRIX","CTC_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);

  {
    Range CrackSurfacesFaces,CrackFrontEdges;
    ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,CrackSurfacesFaces,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,CrackFrontEdges,true); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack Faces = %d\n",CrackSurfacesFaces.size()); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack Front = %d\n",CrackFrontEdges.size()); CHKERRQ(ierr);
    Range CrackFrontNodes;
    rval = mField.get_moab().get_connectivity(CrackFrontEdges,CrackFrontNodes,true); CHKERR_PETSC(rval);
    rval = mField.get_moab().create_meshset(MESHSET_SET,crackForntMeshset); CHKERR_PETSC(rval);
    rval = mField.get_moab().add_entities(crackForntMeshset,CrackFrontNodes); CHKERR_PETSC(rval);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Front Nodes = %d\n",CrackFrontNodes.size()); CHKERRQ(ierr);

    Range SurfacesFaces;
    ierr = mField.get_Cubit_msId_entities_by_dimension(102,SideSet,2,SurfacesFaces,true); CHKERRQ(ierr);
    Range SurfacesNodes;
    rval = mField.get_moab().get_connectivity(SurfacesFaces,SurfacesNodes,true); CHKERR_PETSC(rval);
    Range FrontTangentNodes = subtract(CrackFrontNodes,SurfacesNodes);
    rval = mField.get_moab().create_meshset(MESHSET_SET,crackFrontTangentConstrains); CHKERR_PETSC(rval);	
    rval = mField.get_moab().add_entities(crackFrontTangentConstrains,FrontTangentNodes); CHKERR_PETSC(rval);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Front Tangent Nodes = %d\n",FrontTangentNodes.size()); CHKERRQ(ierr);

    Range CrackSurfacesEdgeFaces;
    rval = mField.get_moab().get_adjacencies(CrackFrontNodes,2,false,CrackSurfacesEdgeFaces,Interface::UNION); CHKERR_PETSC(rval);
    CrackSurfacesEdgeFaces = CrackSurfacesEdgeFaces.subset_by_type(MBTRI);
    CrackSurfacesEdgeFaces = intersect(CrackSurfacesEdgeFaces,CrackSurfacesFaces);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Front Faces = %d\n",CrackSurfacesEdgeFaces.size()); CHKERRQ(ierr);

    Range level_tris;
    ierr = mField.refine_get_ents(*ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
    CrackSurfacesEdgeFaces = intersect(CrackSurfacesEdgeFaces,level_tris);

    ierr = mField.seed_finite_elements(CrackSurfacesEdgeFaces); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(CrackSurfacesEdgeFaces,"C_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(CrackSurfacesEdgeFaces,"CTC_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(CrackSurfacesEdgeFaces,"dCT_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(CrackSurfacesEdgeFaces,"C_TANGENT_ELEM"); CHKERRQ(ierr);

    Range level_edges;
    ierr = mField.refine_get_ents(*ptr_bit_level0,BitRefLevel().set(),MBEDGE,level_edges); CHKERRQ(ierr);
    CrackFrontEdges = intersect(CrackFrontEdges,level_edges);
    ierr = mField.seed_finite_elements(CrackFrontEdges); CHKERRQ(ierr);

  }

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_VERTICEs(crackForntMeshset,"LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_VERTICEs(crackFrontTangentConstrains,"LAMBDA_CRACK_TANGENT_CONSTRAIN"); CHKERRQ(ierr);

  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CRACKFRONT_AREA",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CRACK_TANGENT_CONSTRAIN",1); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element(problem,"dCT_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);

  if(problem == "COUPLED_PROBLEM") {
    ierr = mField.modify_problem_add_finite_element(problem,"C_TANGENT_ELEM"); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::spatial_partition_problems(FieldInterface& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
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

  //partition MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS
  ierr = mField.partition_problem("MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::coupled_partition_problems(FieldInterface& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = mField.partition_problem("COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("COUPLED_PROBLEM"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::constrains_partition_problems(FieldInterface& mField,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = mField.partition_problem("CCT_ALL_MATRIX"); CHKERRQ(ierr);
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
  ierr = mField.partition_problem("CTC_CRACKFRONT_MATRIX"); CHKERRQ(ierr);
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

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::solve_spatial_problem(FieldInterface& mField,SNES *snes) {
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

  PostProcVertexMethod ent_method(mField.get_moab(),"SPATIAL_POSITION");
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Col,ent_method); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_res(mField.get_moab(),"SPATIAL_POSITION",F,"SPATIAL_RESIDUAL");
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Col,ent_method_res); CHKERRQ(ierr);

  if(fe_post_proc_stresses_method!=NULL) delete fe_post_proc_stresses_method;
  fe_post_proc_stresses_method = new PostProcStressNonLinearElasticity(mField.get_moab(),MySpatialFE);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",*fe_post_proc_stresses_method);  CHKERRQ(ierr);

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::surface_projection_data(FieldInterface& mField,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  Interface& moab = mField.get_moab();

  ierr = delete_surface_projection_data(mField); CHKERRQ(ierr);

  if(projSurfaceCtx==NULL) {
    projSurfaceCtx = new matPROJ_ctx(mField,problem,"C_ALL_MATRIX");
    ierr = mField.MatCreateMPIAIJWithArrays("C_ALL_MATRIX",&projSurfaceCtx->C); CHKERRQ(ierr);
  }

  Range CornersEdges,CornersNodes;
  ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,CornersEdges,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,CornersNodes,true); CHKERRQ(ierr);
  ErrorCode rval;
  Range CornersEdgesNodes;
  rval = mField.get_moab().get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
  CornersNodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
  materialDirihletBC myDirihletBC(mField.get_moab(),CornersNodes);

  C_SURFACE_FEMethod CFE_SURFACE(moab,&myDirihletBC,projSurfaceCtx->C);
  C_SURFACE_FEMethod CFE_CRACK_SURFACE(moab,&myDirihletBC,projSurfaceCtx->C,"LAMBDA_CRACK_SURFACE");
  CFE_SURFACE.updated = true;
  CFE_CRACK_SURFACE.updated = true;

  ierr = MatSetOption(projSurfaceCtx->C,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(projSurfaceCtx->C,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

  ierr = MatZeroEntries(projSurfaceCtx->C); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_SURFACE_ELEM",CFE_SURFACE);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_CRACK_SURFACE_ELEM",CFE_CRACK_SURFACE);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(projSurfaceCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(projSurfaceCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

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

  double gc;
  PetscBool flg;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is fracture energy ?)");
  }

  Vec LambdaVec,GriffithForceVec;
  ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",Row,&LambdaVec); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",Col,&GriffithForceVec); CHKERRQ(ierr);
  ierr = VecSet(LambdaVec,gc); CHKERRQ(ierr);

  Mat Q;
  {
    int M,m;
    ierr = VecGetSize(GriffithForceVec,&M); CHKERRQ(ierr);
    ierr = VecGetLocalSize(GriffithForceVec,&m); CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,projSurfaceCtx,&Q); CHKERRQ(ierr);
    ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))matQ_mult_shell); CHKERRQ(ierr);
  }

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

  ierr = mField.set_other_global_VecCreateGhost(problem,"MESH_NODE_POSITIONS","GRIFFITH_FORCE",Row,QTGriffithForceVec,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_res(mField.get_moab(),"MESH_NODE_POSITIONS",QTGriffithForceVec,"GRIFFITH_FORCE");
  ierr = mField.loop_dofs(problem,"MESH_NODE_POSITIONS",Col,ent_method_res); CHKERRQ(ierr);

  double nrm2_griffith_force;
  ierr = VecNorm(QTGriffithForceVec,NORM_2,&nrm2_griffith_force);   CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"nrm2_QTGriffithForceVec = %6.4e\n",nrm2_griffith_force); CHKERRQ(ierr);

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

  double gc;
  PetscBool flg;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is fracture energy ?)");
  }

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

  //if(projFrontCtx->C != PETSC_NULL) {
    //ierr = MatDestroy(&projFrontCtx->C); CHKERRQ(ierr);
  //}
  //ierr = mField.MatCreateMPIAIJWithArrays("C_CRACKFRONT_MATRIX",&projFrontCtx->C); CHKERRQ(ierr);

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
  
  ierr = projFrontCtx->InitQorP(F_Material); CHKERRQ(ierr);
  // unit of LambdaVec [ N * 1/m = N*m/m^2 = J/m^2 ]
  ierr = VecScale(F_Material,1./gc); CHKERRQ(ierr);
  ierr = MatMult(RT,F_Material,LambdaVec); CHKERRQ(ierr);  
  ierr = VecScale(F_Material,gc); CHKERRQ(ierr);
  ierr = VecScale(LambdaVec,gc); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(LambdaVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(LambdaVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost("C_CRACKFRONT_MATRIX",Row,LambdaVec,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  Vec JVec;
  ierr = VecDuplicate(LambdaVec,&JVec); CHKERRQ(ierr);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);

  ErrorCode rval;
  ublas::vector<double> coords;
  coords.resize(3);
  PetscPrintf(PETSC_COMM_WORLD,"\n\ngriffith force\n\n");

  Tag th_g;
  const double def_val = 0;
  rval = mField.get_moab().tag_get_handle("G",1,MB_TYPE_DOUBLE,th_g,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERR_THROW(rval);

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
    ss << "\t\tg " << scientific << setprecision(4) << -diit->get_FieldData();
    ss << " / " << scientific << setprecision(4) << j;
    ss << " ( " << scientific << setprecision(4) << -diit->get_FieldData()/j << " )";
    ss << "\t relative error (gc-g)/gc " << scientific << setprecision(4) << (gc+diit->get_FieldData())/gc;
    ss << " / " << scientific << setprecision(4) << (gc-j)/gc;
    ss << endl; 
    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());

    double val = -diit->get_FieldData();
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
  ierr = mField.refine_get_ents(*ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
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
  PetscPrintf(PETSC_COMM_WORLD,"\naverage griffith force %6.4e / %6.4e Crack surface area %6.4e\n",-ave_g,ave_j,aRea);
  PetscPrintf(PETSC_COMM_WORLD,"\n\n");

  PostProcVertexMethod ent_method(mField.get_moab(),"LAMBDA_CRACKFRONT_AREA");
  ierr = mField.loop_dofs("C_CRACKFRONT_MATRIX","LAMBDA_CRACKFRONT_AREA",Row,ent_method); CHKERRQ(ierr);

  /*// unit of F_Griffith [ N/m * m = N ]
  ierr = MatMult(projFrontCtx->CT,LambdaVec,F_Griffith); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_Griffith,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Griffith,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_griffith(mField.get_moab(),"MESH_NODE_POSITIONS",F_Griffith,"G");
  ierr = mField.loop_dofs(problem,"MESH_NODE_POSITIONS",Row,ent_method_griffith); CHKERRQ(ierr);*/

  ierr = MatDestroy(&Q); CHKERRQ(ierr);
  ierr = MatDestroy(&RT); CHKERRQ(ierr);
  ierr = VecDestroy(&F_Griffith); CHKERRQ(ierr);
  ierr = VecDestroy(&F_Material); CHKERRQ(ierr);
  ierr = VecDestroy(&LambdaVec); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::solve_material_problem(FieldInterface& mField,SNES *snes) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  ierr = front_projection_data(mField,"MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS"); CHKERRQ(ierr);

  Vec F,D;
  ierr = mField.VecCreateGhost("MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS",Row,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS",Row,&D); CHKERRQ(ierr);
  Mat K;
  ierr = mField.MatCreateMPIAIJWithArrays("MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS",&K); CHKERRQ(ierr);

  const double YoungModulus = 1;
  const double PoissonRatio = 0.;
  double gc;// = 8.2414e-07;//4e-7;
  PetscBool flg;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is fracture energy ?)");
  }

  Range CornersEdges,CornersNodes;
  ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,CornersEdges,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,CornersNodes,true); CHKERRQ(ierr);
  ErrorCode rval;
  Range CornersEdgesNodes;
  rval = mField.get_moab().get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
  CornersNodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
  materialDirihletBC myDirihletBC(mField.get_moab(),CornersNodes);

  NL_MaterialFEMethod MyMaterialFE(mField,&myDirihletBC,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  C_SURFACE_FEMethod_ForSnes MySurfaceConstrainsBodySurface(mField,&myDirihletBC,"LAMBDA_SURFACE");
  C_SURFACE_FEMethod_ForSnes MySurfaceConstrainsCrackSurface(mField,&myDirihletBC,"LAMBDA_CRACK_SURFACE");
  Snes_CTgc_CONSTANT_AREA_FEMethod MyCTgc(mField,*projFrontCtx,"MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS","LAMBDA_CRACKFRONT_AREA");
  Snes_dCTgc_CONSTANT_AREA_FEMethod MydCTgc(mField,K,"LAMBDA_CRACKFRONT_AREA");

  SnesCtx snes_ctx(mField,"MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS");

  SnesCtx::loops_to_do_type& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
  snes_ctx.get_preProcess_to_do_Rhs().push_back(&MyCTgc);
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MATERIAL",&MyMaterialFE));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("CandCT_SURFACE_ELEM",&MySurfaceConstrainsBodySurface));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM",&MySurfaceConstrainsCrackSurface));

  SnesCtx::loops_to_do_type& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MATERIAL",&MyMaterialFE));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("CandCT_SURFACE_ELEM",&MySurfaceConstrainsBodySurface));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM",&MySurfaceConstrainsCrackSurface));

  ierr = SNESSetApplicationContext(*snes,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(*snes,F,SnesRhs,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(*snes,K,K,SnesMat,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(*snes); CHKERRQ(ierr);

  ierr = mField.set_local_VecCreateGhost("MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS",Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = SNESSolve(*snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(*snes,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&K); CHKERRQ(ierr);

  ierr = delete_front_projection_data(mField); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::solve_coupled_problem(FieldInterface& mField,SNES *snes,double da) {
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

  Range CornersEdges,CornersNodes;
  ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,CornersEdges,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,CornersNodes,true); CHKERRQ(ierr);
  ErrorCode rval;
  Range CornersEdgesNodes;
  rval = mField.get_moab().get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
  CornersNodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
  CubitDisplacementDirihletBC_Coupled myDirihletBC(mField,"COUPLED_PROBLEM",CornersNodes);
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
  TangentFrontConstrain_FEMethod MyTangentFrontContrain(mField,"LAMBDA_CRACK_TANGENT_CONSTRAIN");

  ////******
  ierr = MySpatialFE.initCrackFrontData(mField); CHKERRQ(ierr);
  ierr = MyMaterialFE.initCrackFrontData(mField); CHKERRQ(ierr);
  ierr = MyMeshSmoother.initCrackFrontData(mField); CHKERRQ(ierr);
  MySurfaceConstrainsBodySurface.nonlinear = true;
  MySurfaceConstrainsCrackSurface.nonlinear = true;
  ////******

  SnesCtx::loops_to_do_type& loops_to_do_Rhs = arc_snes_ctx.get_loops_to_do_Rhs();
  arc_snes_ctx.get_preProcess_to_do_Rhs().push_back(&MyCTgc);
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ELASTIC_COUPLED",&MySpatialFE));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MATERIAL_COUPLED",&MyMaterialFE));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MESH_SMOOTHER",&MyMeshSmoother));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("CandCT_SURFACE_ELEM",&MySurfaceConstrainsBodySurface));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM",&MySurfaceConstrainsCrackSurface));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("C_TANGENT_ELEM",&MyTangentFrontContrain));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&arc_elem));

  SnesCtx::loops_to_do_type& loops_to_do_Mat = arc_snes_ctx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ELASTIC_COUPLED",&MySpatialFE));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MATERIAL_COUPLED",&MyMaterialFE));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MESH_SMOOTHER",&MyMeshSmoother));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("CandCT_SURFACE_ELEM",&MySurfaceConstrainsBodySurface));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM",&MySurfaceConstrainsCrackSurface));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("dCT_CRACKFRONT_AREA_ELEM",&MydCTgc));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("C_TANGENT_ELEM",&MyTangentFrontContrain));
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
  ierr = SNESSetFromOptions(*snes); CHKERRQ(ierr);

  KSP ksp;
  ierr = SNESGetKSP(*snes,&ksp); CHKERRQ(ierr);
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
  ierr = arc_ctx.set_s(arc_elem.lambda_int+da); CHKERRQ(ierr);

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

  //Seve field on mesh
  PostProcVertexMethod ent_method_spatial(mField.get_moab(),"SPATIAL_POSITION");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","SPATIAL_POSITION",Col,ent_method_spatial); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_material(mField.get_moab(),"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","MESH_NODE_POSITIONS",Col,ent_method_material); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_res_mat(mField.get_moab(),"MESH_NODE_POSITIONS",F,"MATERIAL_RESIDUAL");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","MESH_NODE_POSITIONS",Col,ent_method_res_mat); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_res_spat(mField.get_moab(),"SPATIAL_POSITION",F,"SPATIAL_RESIDUAL");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","SPATIAL_POSITION",Col,ent_method_res_spat); CHKERRQ(ierr);

  ierr = MatDestroy(&ShellK); CHKERRQ(ierr);
  ierr = MatDestroy(&K); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);

  ierr = delete_front_projection_data(mField); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode ConfigurationalFractureMechanics::calculate_spatial_residual(FieldInterface& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  DirihletBCMethod_DriverComplexForLazy myDirihletBCSpatial(mField,"ELASTIC_MECHANICS","SPATIAL_POSITION");
  ierr = myDirihletBCSpatial.Init(); CHKERRQ(ierr);

  const double YoungModulus = 1;
  const double PoissonRatio = 0.;
  NL_ElasticFEMethod MySpatialFE(mField,&myDirihletBCSpatial,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  Vec F_Spatial;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&F_Spatial); CHKERRQ(ierr);
  ierr = VecZeroEntries(F_Spatial); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_Spatial,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Spatial,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  MySpatialFE.snes_f = F_Spatial;
  MySpatialFE.set_snes_ctx(FieldInterface::SnesMethod::ctx_SNESSetFunction);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MySpatialFE);  CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_Spatial,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Spatial,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F_Spatial); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F_Spatial); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_res_spatial(mField.get_moab(),"SPATIAL_POSITION",F_Spatial,"SPATIAL_RESIDUAL");
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Col,ent_method_res_spatial); CHKERRQ(ierr);

  ierr = VecDestroy(&F_Spatial); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::calculate_material_forces(FieldInterface& mField,string problem,string fe) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  Range CornersEdges,CornersNodes;
  ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,CornersEdges,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,CornersNodes,true); CHKERRQ(ierr);
  ErrorCode rval;
  Range CornersEdgesNodes;
  rval = mField.get_moab().get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
  CornersNodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
  CubitDisplacementDirihletBC_Coupled myDirihletBC(mField,problem,CornersNodes);
  ierr = myDirihletBC.Init(); CHKERRQ(ierr);

  const double YoungModulus = 1;
  const double PoissonRatio = 0.;
  NL_MaterialFEMethod MyMaterialFE(mField,&myDirihletBC,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));

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
    problem,"MESH_NODE_POSITIONS","MATERIAL_FORCE",Row,F_Material,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"MATERIAL_FORCE",dof)) {
    //cerr << *dof << endl;
  //}

  //detroy matrices
  ierr = VecDestroy(&F_Material); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::save_edge_lenght_in_tags(FieldInterface& mField,BitRefLevel mask) {
  PetscFunctionBegin;
  ErrorCode rval;
  PetscErrorCode ierr;
  Range bit_level_edges;
  ierr = mField.refine_get_ents(*ptr_bit_level0,mask,MBEDGE,bit_level_edges); CHKERRQ(ierr);
  double def_VAL[1] = {0};
  Tag th_edge_length;
  rval = mField.get_moab().tag_get_handle("EDGE_LENGTH",1,MB_TYPE_DOUBLE,th_edge_length,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR(rval);
  vector<double> edges_length(bit_level_edges.size());
  Range::iterator eit = bit_level_edges.begin();
  for(int ee = 0;eit!=bit_level_edges.end();eit++,ee++) {
    const EntityHandle* conn; 
    int num_nodes; 
    rval = mField.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
    double coords[6];
    rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    cblas_daxpy(3,-1,&coords[3],1,&coords[0],1);
    edges_length[ee] = cblas_dnrm2(3,coords,1);
    if(edges_length[ee] <= 0) SETERRQ(PETSC_COMM_SELF,1,"edge length is 0 (or negative)");
  }
  rval = mField.get_moab().tag_set_data(th_edge_length,bit_level_edges,&edges_length[0]); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::save_edge_strech_lenght_in_tags(FieldInterface& mField) {
  PetscFunctionBegin;
  ErrorCode rval;
  PetscErrorCode ierr;
  Range bit_level_edges;
  ierr = mField.refine_get_ents(*ptr_bit_level0,BitRefLevel().set(),MBEDGE,bit_level_edges); CHKERRQ(ierr);
  double def_VAL[1] = {0};
  Tag th_edge_strech;
  rval = mField.get_moab().tag_get_handle("EDGE_STRECH",1,MB_TYPE_DOUBLE,th_edge_strech,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR(rval);
  vector<double> edges_strech(bit_level_edges.size());
  Range::iterator eit = bit_level_edges.begin();
  for(int ee = 0;eit!=bit_level_edges.end();eit++,ee++) {
    const EntityHandle* conn; 
    int num_nodes; 
    rval = mField.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
    double coords[6];
    rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    cblas_daxpy(3,-1,&coords[3],1,&coords[0],1);
    edges_strech[ee] = cblas_dnrm2(3,coords,1);
  }
  Tag th_edge_length;
  rval = mField.get_moab().tag_get_handle("EDGE_LENGTH",th_edge_length); CHKERR(rval);
  vector<double> edges_length0(bit_level_edges.size());
  rval = mField.get_moab().tag_get_data(th_edge_length,bit_level_edges,&*edges_length0.begin()); CHKERR_PETSC(rval);
  vector<double>::iterator vit_strech = edges_strech.begin();
  vector<double>::iterator vit_length0 = edges_length0.begin();
  for(;vit_strech != edges_strech.end();vit_strech++,vit_length0++) {
    if(*vit_length0 <= 0) SETERRQ(PETSC_COMM_SELF,1,"edge length is 0 (or negative)");
    *vit_strech /= *vit_length0;
  }
  rval = mField.get_moab().tag_set_data(th_edge_strech,bit_level_edges,&*edges_strech.begin()); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::refine_streched_edges(FieldInterface& mField,double strech_treshold,bool crack_crack_surface_only) {
  PetscFunctionBegin;
  try {
  ErrorCode rval;
  PetscErrorCode ierr;
  int last = 0;
  for(unsigned int ii = 0;ii<ptr_bit_level0->size();ii++) {
    if(ptr_bit_level0->test(ii)) last = ii;
  }
  last++;
  Range bit_level_edges;
  ierr = mField.refine_get_ents(*ptr_bit_level0,BitRefLevel().set(),MBEDGE,bit_level_edges); CHKERRQ(ierr);
  if(crack_crack_surface_only) {
    Range crack_surfaces_faces,adj_cracksurfaces_edges;
    ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,crack_surfaces_faces,true); CHKERRQ(ierr);
    rval = mField.get_moab().get_adjacencies(crack_surfaces_faces,1,false,adj_cracksurfaces_edges,Interface::UNION); CHKERR_PETSC(rval);
    bit_level_edges = intersect(bit_level_edges,adj_cracksurfaces_edges);
  }
  Tag th_edge_strech;
  rval = mField.get_moab().tag_get_handle("EDGE_STRECH",th_edge_strech); CHKERR_PETSC(rval);
  Range edges_to_refine;
  Range::iterator eit = bit_level_edges.begin();
  for(int ee = 0;eit!=bit_level_edges.end();eit++,ee++) {
    double edge_strech;
    rval = mField.get_moab().tag_get_data(th_edge_strech,&*eit,1,&edge_strech); CHKERR_PETSC(rval);
    if(edge_strech<strech_treshold) continue;
    edges_to_refine.insert(*eit);
  }
  BitRefLevel last_ref = BitRefLevel().set(last);
  ierr = mField.add_verices_in_the_middel_of_edges(edges_to_refine,last_ref,2); CHKERRQ(ierr);
  Range level_tets;
  ierr = mField.refine_get_ents(*ptr_bit_level0,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
  ierr = mField.refine_TET(level_tets,last_ref,false); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
    EntityHandle cubit_meshset = cubit_it->meshset; 
    ierr = mField.refine_get_childern(cubit_meshset,last_ref,cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(cubit_meshset,last_ref,cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(cubit_meshset,last_ref,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(cubit_meshset,last_ref,cubit_meshset,MBTET,true); CHKERRQ(ierr);
  }
  *ptr_bit_level0 = last_ref;
  ierr = save_edge_lenght_in_tags(mField,*ptr_bit_level0); CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
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
    ierr = mField.refine_get_ents(*conf_prob->ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    crackSurfacesFaces = intersect(crackSurfacesFaces,level_tris);

    //get pointer to coupled problem
    ierr = mField.get_problem("COUPLED_PROBLEM",&problem_ptr); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    set<DofIdx> set_idx; //set of global dofs indexes for crack surface
    for(Range::iterator fit = crackSurfacesFaces.begin();fit!=crackSurfacesFaces.end();fit++) {
      const EntityHandle* conn; 
      int num_nodes; 
      rval = mField.get_moab().get_connectivity(*fit,conn,num_nodes,true); 
      if (MB_SUCCESS != rval) { 
	PetscAbortErrorHandler(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,__SDIR__,rval,PETSC_ERROR_INITIAL,
	  "can not get connectibity",PETSC_NULL);
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
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::calulate_lambda_int() {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = calulate_area(); CHKERRQ(ierr);
  lambda_int = arc_ptr->alpha*aRea + arc_ptr->dlambda*arc_ptr->beta*sqrt(arc_ptr->F_lambda2);
  PetscPrintf(PETSC_COMM_WORLD,"\tsurface crack area = %6.4e lambda_int = %6.4e\n",aRea,lambda_int);
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::calulate_db() {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  if(arc_ptr->alpha!=0) {
    ierr = MatMultTranspose(conf_prob->projFrontCtx->C,lambdaVec,arc_ptr->db); CHKERRQ(ierr);
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
    PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e, %6.4e (%6.4e)\n",lambda_old, array[arc_ptr->get_petsc_local_dof_idx()], dlambda);
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
  PetscPrintf(PETSC_COMM_WORLD,"\tdlambda = %6.4e\n",arc_ptr->dlambda);
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
	ierr = VecSetValue(snes_f,arc_ptr->get_petsc_gloabl_dof_idx(),arc_ptr->res_lambda,ADD_VALUES); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_SELF,"\tres_lambda = %6.4e lambda_int = %6.4e\n",arc_ptr->res_lambda,lambda_int);
      } break; 
      case ctx_SNESSetJacobian: {
	//calulate diagonal therm
	double diag = arc_ptr->beta*sqrt(arc_ptr->F_lambda2);
	ierr = VecSetValue(ghostDiag,0,diag,INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(*snes_B,arc_ptr->get_petsc_gloabl_dof_idx(),arc_ptr->get_petsc_gloabl_dof_idx(),1,ADD_VALUES); CHKERRQ(ierr);
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
	  if(dof->get_name() == "LAMBDA_SURFACE") {
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
	  PetscPrintf(PETSC_COMM_WORLD,"\tres spatial nrm2 = %6.4e\n",sqrt(res_nrm2[0]));
	  PetscPrintf(PETSC_COMM_WORLD,"\tres crack front nrm2 = %6.4e\n",sqrt(res_nrm2[1]));
	  PetscPrintf(PETSC_COMM_WORLD,"\tres mesh smoother nrm2 = %6.4e\n",sqrt(res_nrm2[2]));
	  PetscPrintf(PETSC_COMM_WORLD,"\tres surface constrain nrm2 = %6.4e\n",sqrt(res_nrm2[3]));
	  PetscPrintf(PETSC_COMM_WORLD,"\tres crack surface constrain nrm2 = %6.4e\n",sqrt(res_nrm2[4]));
	  PetscPrintf(PETSC_COMM_WORLD,"\tres crack front tangent constrain nrm2 = %6.4e\n",sqrt(res_nrm2[5]));
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
	PetscPrintf(PETSC_COMM_WORLD,"\tdiag = %6.4e\n",arc_ptr->diag);
      } break;
      default:
	break;
    }
    PetscFunctionReturn(0);
}

