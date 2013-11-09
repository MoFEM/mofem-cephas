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
  
struct NL_MaterialFEMethodProjected: public FEMethod_DriverComplexForLazy_MaterialProjected {
  
    NL_MaterialFEMethodProjected(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0): 
        FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
        FEMethod_DriverComplexForLazy_MaterialProjected(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose)  {
      set_PhysicalEquationNumber(eq_solid);
    }
  
  };
  
struct NL_ElasticFEMethodCoupled: public FEMethod_DriverComplexForLazy_CoupledSpatial {
  
  ArcLengthCtx *arc_ptr;
  NL_ElasticFEMethodCoupled(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,ArcLengthCtx *_arc_ptr = NULL,int _verbose = 0):
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
      FEMethod_DriverComplexForLazy_CoupledSpatial(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose), arc_ptr(_arc_ptr) {
	set_PhysicalEquationNumber(eq_solid);
      }

  PetscErrorCode operator()() {
    PetscFunctionBegin;

    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndicesSpatial(); CHKERRQ(ierr);

    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(this,RowGlobSpatial,ColGlobSpatial,DirihletBC); CHKERRQ(ierr);
    ierr = calulateKFint(proj_all_ctx.K,snes_f); CHKERRQ(ierr);
    if(arc_ptr==NULL) {
      ierr = calulateKFext(proj_all_ctx.K,snes_f,*(this->t_val)); CHKERRQ(ierr);
    } else {
      ierr = calulateKFext(proj_all_ctx.K,arc_ptr->F_lambda,-1.); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }
  
};

struct NL_ElasticFEMethodCoupled_F_lambda_Only: public FEMethod_DriverComplexForLazy_CoupledSpatial {

  ArcLengthCtx *arc_ptr;
  NL_ElasticFEMethodCoupled_F_lambda_Only(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,ArcLengthCtx *_arc_ptr = NULL,int _verbose = 0):
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
      FEMethod_DriverComplexForLazy_CoupledSpatial(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose), arc_ptr(_arc_ptr) {
	set_PhysicalEquationNumber(eq_solid);
      }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;

    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndicesSpatial(); CHKERRQ(ierr);
    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(this,RowGlobSpatial,ColGlobSpatial,DirihletBC); CHKERRQ(ierr);
    snes_ctx = ctx_SNESSetFunction;
    ierr = calulateKFext(PETSC_NULL,arc_ptr->F_lambda,-1); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};
  
struct NL_MaterialFEMethodCoupled: public FEMethod_DriverComplexForLazy_CoupledMaterial {
  
    ArcLengthCtx *arc_ptr;
    NL_MaterialFEMethodCoupled(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,ArcLengthCtx *_arc_ptr = NULL,int _verbose = 0):
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
      FEMethod_DriverComplexForLazy_CoupledMaterial(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose),arc_ptr(_arc_ptr) {
	set_PhysicalEquationNumber(eq_solid);
      }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      switch(snes_ctx) {
	case ctx_SNESSetFunction:  
	  if(arc_ptr==NULL) {
	    ierr = FEMethod_DriverComplexForLazy_Material::operator()(&(proj_all_ctx.K),snes_f,RHS|LHS|TANGENT|FEXT|FINT); CHKERRQ(ierr);
	  } else {
	    ierr = FEMethod_DriverComplexForLazy_Material::operator()(&(proj_all_ctx.K),snes_f,RHS|LHS|TANGENT|FINT); CHKERRQ(ierr);
	    ierr = FEMethod_DriverComplexForLazy_Material::operator()(&(proj_all_ctx.K),arc_ptr->F_lambda,RHS|FEXT); CHKERRQ(ierr);
	  }
	  ierr = GetIndicesRow(RowGlobSpatial,spatial_field_name); CHKERRQ(ierr);
	  ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndiciesRow(this,RowGlobSpatial,DirihletBC); CHKERRQ(ierr);
	  ierr = AssembleMaterialCoupledTangent(proj_all_ctx.K); CHKERRQ(ierr);
	break;
	case ctx_SNESSetJacobian:
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }

      PetscFunctionReturn(0);
    }

  
  };

struct NL_MaterialFEMethodCoupled_F_lambda_Only: public FEMethod_DriverComplexForLazy_CoupledMaterial {
  
    ArcLengthCtx *arc_ptr;
    NL_MaterialFEMethodCoupled_F_lambda_Only(
      FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,ArcLengthCtx *_arc_ptr = NULL,int _verbose = 0):
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
      FEMethod_DriverComplexForLazy_CoupledMaterial(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose),arc_ptr(_arc_ptr) {
	set_PhysicalEquationNumber(eq_solid);
      }

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = FEMethod_DriverComplexForLazy_Material::operator()(&(proj_all_ctx.K),arc_ptr->F_lambda,RHS|FEXT); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

  
  };


struct NL_MeshSmootherCoupled: public FEMethod_DriverComplexForLazy_CoupledMeshSmoother {
  
    NL_MeshSmootherCoupled(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,double _alpha3,int _verbose = 0):
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
      FEMethod_DriverComplexForLazy_CoupledMeshSmoother(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,_alpha3) {
	set_qual_ver(0);
      }
  
  };

struct NL_Projection: public FEMethod_DriverComplexForLazy_CoupledProjected {

  ArcLengthCtx *arc_ptr;

  NL_Projection(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,string _problem_name,ArcLengthCtx *_arc_ptr = NULL,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    FEMethod_DriverComplexForLazy_CoupledProjected(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,_problem_name,_verbose),arc_ptr(_arc_ptr) {};

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy_CoupledProjected::preProcess(); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	if(arc_ptr != NULL) {
	  ierr = VecZeroEntries(arc_ptr->F_lambda); CHKERRQ(ierr);
	  ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	  ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
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

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy_CoupledProjected::postProcess(); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	if(arc_ptr != NULL) {
	  ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	  ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	  //add F_lambda
	  ierr = VecAXPY(snes_f,-arc_ptr->get_FieldData(),arc_ptr->F_lambda); CHKERRQ(ierr);
	  PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",arc_ptr->get_FieldData());  
	  //F_lambda2
	  ierr = VecDot(arc_ptr->F_lambda,arc_ptr->F_lambda,&arc_ptr->F_lambda2); CHKERRQ(ierr);
	  PetscPrintf(PETSC_COMM_WORLD,"\tFlambda2 = %6.4e\n",arc_ptr->F_lambda2);
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

};
  
PetscErrorCode ConfigurationalFractureMechanics::CubitDisplacementDirihletBC_Coupled::SetDirihletBC_to_ElementIndiciesRow(
      FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirihletBC) {
  PetscFunctionBegin;
  field_name = "SPATIAL_POSITION";
  ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndiciesRow(fe_method_ptr,RowGlobDofs,DirihletBC); CHKERRQ(ierr);
  field_name = "MESH_NODE_POSITIONS";
  ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndiciesRow(fe_method_ptr,RowGlobDofs,DirihletBC); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
  
PetscErrorCode ConfigurationalFractureMechanics::CubitDisplacementDirihletBC_Coupled::SetDirihletBC_to_ElementIndiciesCol(
  FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC) {
  PetscFunctionBegin;
  CubitDisplacementDirihletBC::field_name = "SPATIAL_POSITION";
  ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndiciesCol(fe_method_ptr,ColGlobDofs,DirihletBC); CHKERRQ(ierr);
  CubitDisplacementDirihletBC::field_name = "MESH_NODE_POSITIONS";
  ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndiciesCol(fe_method_ptr,ColGlobDofs,DirihletBC); CHKERRQ(ierr);
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
  ierr = mField.add_field("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);

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

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("MATERIAL_MECHANICS","MATERIAL"); CHKERRQ(ierr);

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
  ierr = mField.add_finite_element("ELASTIC_COUPLED"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("MATERIAL_COUPLED"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("MESH_SMOOTHER"); CHKERRQ(ierr);
  
  //fes definitions
  ierr = mField.modify_finite_element_add_field_row("ELASTIC_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("ELASTIC_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //
  ierr = mField.modify_finite_element_add_field_row("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("MATERIAL_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
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

  /*//add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);
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
  ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);*/

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
  ierr = mField.add_field("LAMBDA_SURFACE",H1,1); CHKERRQ(ierr);
  ierr = mField.add_field("LAMBDA_CORNER",H1,3); CHKERRQ(ierr);
  //CRACK
  if(cs) {
    ierr = mField.add_field("LAMBDA_CRACK_SURFACE",H1,1); CHKERRQ(ierr);
  }

  //FE
  ierr = mField.add_finite_element("C_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CTC_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("C_CORNER_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CTC_CORNER_ELEM"); CHKERRQ(ierr);
  //CRACK
  if(cs) {
    ierr = mField.add_finite_element("C_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("CTC_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
  }

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("C_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("C_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("C_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("C_CORNER_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_CORNER_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("CTC_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("CTC_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("CTC_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);

  //CRACK
  if(cs) {
    ierr = mField.modify_finite_element_add_field_row("C_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("C_CRACK_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("C_CRACK_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("C_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);

    ierr = mField.modify_finite_element_add_field_row("CTC_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("CTC_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("CTC_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
  }

  //define problems
  ierr = mField.add_problem("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.add_problem("C_ALL_MATRIX"); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_CORNER_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("C_ALL_MATRIX","C_CORNER_ELEM"); CHKERRQ(ierr);
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
    ierr = mField.refine_get_ents(bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
    SurfacesFaces = intersect(SurfacesFaces,level_tris);

    ierr = mField.seed_finite_elements(SurfacesFaces); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"C_SURFACE_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"CTC_SURFACE_ELEM"); CHKERRQ(ierr);

    if(cs) {
      Range CrackEdgeNodes;
      rval = moab.get_connectivity(CrackCornersEdges,CrackEdgeNodes,true); CHKERR_PETSC(rval);
      Range AdjCrackFrontFaces;
      rval = moab.get_adjacencies(CrackEdgeNodes,2,false,AdjCrackFrontFaces,Interface::UNION); CHKERR_PETSC(rval);
      CrackSurfacesFaces = intersect(CrackSurfacesFaces,level_tris);
      CrackSurfacesFaces = subtract(CrackSurfacesFaces,AdjCrackFrontFaces);
      ierr = mField.seed_finite_elements(CrackSurfacesFaces); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TRIs(CrackSurfacesFaces,"C_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TRIs(CrackSurfacesFaces,"CTC_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
    }

    Range CrackCornersEdgesNodes;
    rval = moab.get_connectivity(CrackCornersEdges,CrackCornersEdgesNodes,true); CHKERR_PETSC(rval);
    Range CrackCornerTets;
    rval = moab.get_adjacencies(CrackCornersEdgesNodes,3,false,CrackCornerTets,Interface::UNION); CHKERR_PETSC(rval);
    Range CornersEdgesNodes;
    rval = moab.get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
    {
      CornersNodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
      rval = moab.create_meshset(MESHSET_SET,cornersNodesMeshset); CHKERR_PETSC(rval);	
      rval = moab.add_entities(cornersNodesMeshset,CornersNodes); CHKERR_PETSC(rval);
      //add surface elements
      ierr = mField.seed_finite_elements(CornersNodes); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_VERTICEs(CornersNodes,"C_CORNER_ELEM"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_VERTICEs(CornersNodes,"CTC_CORNER_ELEM"); CHKERRQ(ierr);
    }
    {
      Range SurfacesNodes;
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
  ierr = mField.add_ents_to_field_by_VERTICEs(cornersNodesMeshset,"LAMBDA_CORNER"); CHKERRQ(ierr);
  //CRACK
  if(cs) {
    ierr = mField.add_ents_to_field_by_VERTICEs(crackSurfacesFacesNodesMeshset,"LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
  }

  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_SURFACE",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CORNER",1); CHKERRQ(ierr);
  if(cs) {
    ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CRACK_SURFACE",1); CHKERRQ(ierr);
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
  ierr = mField.add_field("LAMBDA_CRACKFRONT_AREA",H1,1); CHKERRQ(ierr);
  ierr = mField.add_field("GRIFFITH_FORCE",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("C_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CTC_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("dCT_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("FRONT_CONSTRAIN"); CHKERRQ(ierr);

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

  ierr = mField.modify_finite_element_add_field_row("FRONT_CONSTRAIN","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("FRONT_CONSTRAIN","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("FRONT_CONSTRAIN","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

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
    Range CrackSurfacesEdgeFaces;
    rval = mField.get_moab().get_adjacencies(CrackFrontNodes,2,false,CrackSurfacesEdgeFaces,Interface::UNION); CHKERR_PETSC(rval);
    CrackSurfacesEdgeFaces = CrackSurfacesEdgeFaces.subset_by_type(MBTRI);
    CrackSurfacesEdgeFaces = intersect(CrackSurfacesEdgeFaces,CrackSurfacesFaces);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Front Faces = %d\n",CrackSurfacesEdgeFaces.size()); CHKERRQ(ierr);

    Range level_tris;
    ierr = mField.refine_get_ents(bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
    CrackSurfacesEdgeFaces = intersect(CrackSurfacesEdgeFaces,level_tris);

    ierr = mField.seed_finite_elements(CrackSurfacesEdgeFaces); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(CrackSurfacesEdgeFaces,"C_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(CrackSurfacesEdgeFaces,"CTC_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(CrackSurfacesEdgeFaces,"dCT_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);

    Range level_edges;
    ierr = mField.refine_get_ents(bit_level0,BitRefLevel().set(),MBEDGE,level_edges); CHKERRQ(ierr);
    CrackFrontEdges = intersect(CrackFrontEdges,level_edges);



    ierr = mField.seed_finite_elements(CrackFrontEdges); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_EDGEs(CrackFrontEdges,"FRONT_CONSTRAIN"); CHKERRQ(ierr);

  }

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_VERTICEs(crackForntMeshset,"LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);
  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CRACKFRONT_AREA",1); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element(problem,"dCT_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element(problem,"FRONT_CONSTRAIN"); CHKERRQ(ierr);

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

  //partition
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

  C_SURFACE_FEMethod CFE_SURFACE(moab,projSurfaceCtx->C);
  C_SURFACE_FEMethod CFE_CRACK_SURFACE(moab,projSurfaceCtx->C,"LAMBDA_CRACK_SURFACE");
  C_CORNER_FEMethod CFE_CORNER(moab,projSurfaceCtx->C);

  ierr = MatSetOption(projSurfaceCtx->C,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(projSurfaceCtx->C,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

  ierr = MatZeroEntries(projSurfaceCtx->C); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_SURFACE_ELEM",CFE_SURFACE);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_CRACK_SURFACE_ELEM",CFE_CRACK_SURFACE);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(projSurfaceCtx->C,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(projSurfaceCtx->C,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_CORNER_ELEM",CFE_CORNER);  CHKERRQ(ierr);
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

  Vec LambdaVec,GriffithForceVec;
  ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",Row,&LambdaVec); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",Col,&GriffithForceVec); CHKERRQ(ierr);
  ierr = VecSet(LambdaVec,1); CHKERRQ(ierr);

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
  ierr = MatMult(RT,F_Material,LambdaVec); CHKERRQ(ierr);  

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

    double j = norm_2(material_force)/norm_2(griffith_force);
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

  }
  ierr = VecSum(LambdaVec,&ave_g); CHKERRQ(ierr);
  ierr = VecMin(LambdaVec,PETSC_NULL,&min_g); CHKERRQ(ierr);
  ierr = VecMax(LambdaVec,PETSC_NULL,&max_g); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(JVec); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(JVec); CHKERRQ(ierr);
  ierr = VecSum(JVec,&ave_j); CHKERRQ(ierr);
  ierr = VecMin(JVec,PETSC_NULL,&min_j); CHKERRQ(ierr);
  ierr = VecMax(JVec,PETSC_NULL,&max_j); CHKERRQ(ierr);

  {
    int N;
    ierr = VecGetSize(LambdaVec,&N); CHKERRQ(ierr);
    ave_g /= N;
    ave_j /= N;
  }
  PetscPrintf(PETSC_COMM_WORLD,"\naverage griffith force %6.4e / %6.4e\n",-ave_g,ave_j);
  PetscPrintf(PETSC_COMM_WORLD,"\n\n");

  PostProcVertexMethod ent_method(mField.get_moab(),"LAMBDA_CRACKFRONT_AREA");
  ierr = mField.loop_dofs("C_CRACKFRONT_MATRIX","LAMBDA_CRACKFRONT_AREA",Row,ent_method); CHKERRQ(ierr);

  // unit of F_Griffith [ N/m * m = N ]
  ierr = MatMult(projFrontCtx->CT,LambdaVec,F_Griffith); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_Griffith,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Griffith,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_griffith(mField.get_moab(),"MESH_NODE_POSITIONS",F_Griffith,"G");
  ierr = mField.loop_dofs(problem,"MESH_NODE_POSITIONS",Row,ent_method_griffith); CHKERRQ(ierr);

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

  DirihletBCMethod_DriverComplexForLazy myDirihletBCMaterial(mField,"MATERIAL_MECHANICS","MESH_POSITION");
  ierr = myDirihletBCMaterial.Init(); CHKERRQ(ierr);

  ierr = front_projection_data(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = surface_projection_data(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);

  Vec F;
  ierr = mField.VecCreateGhost("MATERIAL_MECHANICS",Row,&F); CHKERRQ(ierr);
  Vec D;
  ierr = mField.VecCreateGhost("MATERIAL_MECHANICS",Col,&D); CHKERRQ(ierr);
  ierr = mField.set_local_VecCreateGhost("MATERIAL_MECHANICS",Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  Mat precK;
  ierr = mField.MatCreateMPIAIJWithArrays("MATERIAL_MECHANICS",&precK); CHKERRQ(ierr);
  ierr = MatDuplicate(precK,MAT_DO_NOT_COPY_VALUES,&projSurfaceCtx->K); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("C_ALL_MATRIX",Row,&projSurfaceCtx->g); CHKERRQ(ierr);

  Mat CTC_QTKQ;
  int M,N,m,n;
  ierr = MatGetSize(projSurfaceCtx->K,&M,&N); CHKERRQ(ierr);
  ierr = MatGetLocalSize(projSurfaceCtx->K,&m,&n); CHKERRQ(ierr);
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,projSurfaceCtx,&CTC_QTKQ); CHKERRQ(ierr);
  ierr = MatShellSetOperation(CTC_QTKQ,MATOP_MULT,(void(*)(void))matCTC_QTKQ_mult_shell); CHKERRQ(ierr);

  Mat C_crack_fornt;
  ierr = mField.MatCreateMPIAIJWithArrays("C_CRACKFRONT_MATRIX",&C_crack_fornt); CHKERRQ(ierr);

  const double YoungModulus = 1;
  const double PoissonRatio = 0.;
  double gc;// = 8.2414e-07;//4e-7;
  PetscBool flg;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is fracture energy ?)");
  }

  double alpha3 = 1;
  NL_MaterialFEMethodProjected MyMaterialFE(
    mField,*projSurfaceCtx,&myDirihletBCMaterial,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),alpha3);

  Snes_CTgc_CONSTANT_AREA_FEMethod MyCTgc(mField,*projFrontCtx,"MATERIAL_MECHANICS","LAMBDA_CRACKFRONT_AREA");
  Snes_dCTgc_CONSTANT_AREA_FEMethod MydCTgc(mField,*projSurfaceCtx,"LAMBDA_CRACKFRONT_AREA");

  FEMethod_DriverComplexForLazy_CoupledProjected Projection(mField,*projSurfaceCtx,&myDirihletBCMaterial,"MATERIAL_MECHANICS");
  SnesCtx snes_ctx(mField,"MATERIAL_MECHANICS");
  
  ierr = SNESSetApplicationContext(*snes,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(*snes,F,SnesRhs,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(*snes,CTC_QTKQ,precK,SnesMat,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(*snes); CHKERRQ(ierr);

  SnesCtx::loops_to_do_type& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MATERIAL",&MyMaterialFE));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("dCT_CRACKFRONT_AREA_ELEM",&MydCTgc));
  SnesCtx::loops_to_do_type& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MATERIAL",&MyMaterialFE));

  snes_ctx.get_preProcess_to_do_Rhs().push_back(&Projection);
  snes_ctx.get_preProcess_to_do_Rhs().push_back(&MyCTgc);
  snes_ctx.get_postProcess_to_do_Rhs().push_back(&Projection);
  snes_ctx.get_preProcess_to_do_Mat().push_back(&Projection);
  snes_ctx.get_postProcess_to_do_Mat().push_back(&Projection);

  ierr = SNESSolve(*snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(*snes,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost("MATERIAL_MECHANICS",Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  PostProcVertexMethod ent_method(mField.get_moab(),"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MATERIAL_MECHANICS","MESH_NODE_POSITIONS",Col,ent_method); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_res(mField.get_moab(),"MESH_NODE_POSITIONS",F,"MATERIAL_RESIDUAL");
  ierr = mField.loop_dofs("MATERIAL_MECHANICS","MESH_NODE_POSITIONS",Col,ent_method_res); CHKERRQ(ierr);

  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&precK); CHKERRQ(ierr);
  ierr = MatDestroy(&CTC_QTKQ); CHKERRQ(ierr);
  ierr = MatDestroy(&C_crack_fornt); CHKERRQ(ierr);
  ierr = MatDestroy(&projSurfaceCtx->K); CHKERRQ(ierr);
  ierr = VecDestroy(&projSurfaceCtx->g); CHKERRQ(ierr);

  ierr = delete_surface_projection_data(mField); CHKERRQ(ierr);
  ierr = delete_front_projection_data(mField); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::solve_coupled_problem(FieldInterface& mField,SNES *snes,double alpha3,double da) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  ierr = front_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = surface_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);

  //create matrices
  Mat precK;
  ierr = mField.MatCreateMPIAIJWithArrays("COUPLED_PROBLEM",&precK); CHKERRQ(ierr);
  ierr = MatDuplicate(precK,MAT_DO_NOT_COPY_VALUES,&projSurfaceCtx->K); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("C_ALL_MATRIX",Row,&projSurfaceCtx->g); CHKERRQ(ierr);

  Mat CTC_QTKQ;
  int M,N,m,n;
  ierr = MatGetSize(projSurfaceCtx->K,&M,&N); CHKERRQ(ierr);
  ierr = MatGetLocalSize(projSurfaceCtx->K,&m,&n); CHKERRQ(ierr);
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,projSurfaceCtx,&CTC_QTKQ); CHKERRQ(ierr);
  ierr = MatShellSetOperation(CTC_QTKQ,MATOP_MULT,(void(*)(void))matCTC_QTKQ_mult_shell); CHKERRQ(ierr);

  Vec F;
  ierr = mField.VecCreateGhost("COUPLED_PROBLEM",Row,&F); CHKERRQ(ierr);
  Vec D;
  ierr = mField.VecCreateGhost("COUPLED_PROBLEM",Col,&D); CHKERRQ(ierr);

  Mat C_crack_fornt;
  ierr = mField.MatCreateMPIAIJWithArrays("C_CRACKFRONT_MATRIX",&C_crack_fornt); CHKERRQ(ierr);

  CubitDisplacementDirihletBC_Coupled myDirihletBC(mField,"COUPLED_PROBLEM");
  ierr = myDirihletBC.Init(); CHKERRQ(ierr);

  SnesCtx* snes_ctx;
  Mat Arc_CTC_QTKQ;
  ArcLengthCtx* arc_ctx = NULL;
  ArcLengthSnesCtx* arc_snes_ctx;
  ArcLengthMatShell* arc_mat_ctx;
  PCShellCtx* pc_ctx;
  ArcLengthElemFEMethod *arc_elem;
  if(material_FirelWall->operator[](FW_arc_lenhghat_definition)) {
    arc_ctx = new ArcLengthCtx(mField,"COUPLED_PROBLEM");
    arc_snes_ctx = new ArcLengthSnesCtx(mField,"COUPLED_PROBLEM",arc_ctx);
    snes_ctx = arc_snes_ctx;
    arc_mat_ctx = new ArcLengthMatShell(mField,CTC_QTKQ,arc_ctx,"COUPLED_PROBLEM");
    ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,(void*)arc_mat_ctx,&Arc_CTC_QTKQ); CHKERRQ(ierr);
    ierr = MatShellSetOperation(Arc_CTC_QTKQ,MATOP_MULT,(void(*)(void))arc_length_mult_shell); CHKERRQ(ierr);
    pc_ctx = new PCShellCtx(CTC_QTKQ,precK,arc_ctx);
    arc_elem = new ArcLengthElemFEMethod(mField,this,arc_ctx);
  } else {
    snes_ctx = new SnesCtx(mField,"COUPLED_PROBLEM");
  }

  const double YoungModulus = 1;
  const double PoissonRatio = 0.;
  NL_ElasticFEMethodCoupled MySpatialFE(
    mField,*projSurfaceCtx,&myDirihletBC,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),arc_ctx);
  NL_MaterialFEMethodCoupled MyMaterialFE(
    mField,*projSurfaceCtx,&myDirihletBC,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  NL_MeshSmootherCoupled MyMeshSmoother(mField,*projSurfaceCtx,&myDirihletBC,1.);
  ConstrainCrackForntEdges_FEMethod FrontPenalty(mField,this,alpha3);
 
  ////******
  ierr = MySpatialFE.init_crack_front_data(); CHKERRQ(ierr);
  ierr = MyMaterialFE.init_crack_front_data(); CHKERRQ(ierr);
  ierr = MyMeshSmoother.init_crack_front_data(false,false); CHKERRQ(ierr);
  ////******

  if(material_FirelWall->operator[](FW_arc_lenhghat_definition)) {
    arc_ctx->get_FieldData() = *(MySpatialFE.t_val);
  }
  ierr = mField.set_local_VecCreateGhost("COUPLED_PROBLEM",Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  if(material_FirelWall->operator[](FW_arc_lenhghat_definition)) {
    ierr = mField.get_problem("COUPLED_PROBLEM",&(arc_elem->problem_ptr)); CHKERRQ(ierr);
    ierr = arc_elem->set_dlambda_to_x(D,0);/**(MySpatialFE.t_val));*/ CHKERRQ(ierr);
    ierr = VecCopy(D,arc_ctx->x0); CHKERRQ(ierr);
    ierr = arc_elem->get_dlambda(D); CHKERRQ(ierr);
    ierr = arc_ctx->set_alpha_and_beta(1,0); CHKERRQ(ierr);
    NL_ElasticFEMethodCoupled_F_lambda_Only MySpatialFE_F_lambda_Only(
      mField,*projSurfaceCtx,&myDirihletBC,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),arc_ctx);
    NL_MaterialFEMethodCoupled_F_lambda_Only MyMaterialFE_F_lambda_Only(
      mField,*projSurfaceCtx,&myDirihletBC,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),arc_ctx);
    ierr = MyMaterialFE_F_lambda_Only.init_crack_front_data(); CHKERRQ(ierr);
    ierr = VecZeroEntries(arc_ctx->F_lambda); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(arc_ctx->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(arc_ctx->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(arc_ctx->F_lambda); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(arc_ctx->F_lambda); CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("COUPLED_PROBLEM","ELASTIC_COUPLED",MySpatialFE_F_lambda_Only);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("COUPLED_PROBLEM","MATERIAL_COUPLED",MyMaterialFE_F_lambda_Only);  CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(arc_ctx->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(arc_ctx->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(arc_ctx->F_lambda); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(arc_ctx->F_lambda); CHKERRQ(ierr);
    //F_lambda2
    ierr = VecDot(arc_ctx->F_lambda,arc_ctx->F_lambda,&arc_ctx->F_lambda2); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"\tFlambda2 = %6.4e\n",arc_ctx->F_lambda2);
    //set s
    arc_elem->snes_x = D;
    ierr = arc_elem->calulate_lambda_int(); CHKERRQ(ierr);
    ierr = arc_ctx->set_s(arc_elem->lambda_int+da); CHKERRQ(ierr);
  } else {
  }

  Snes_CTgc_CONSTANT_AREA_FEMethod MyCTgc(mField,*projFrontCtx,"COUPLED_PROBLEM","LAMBDA_CRACKFRONT_AREA");
  Snes_dCTgc_CONSTANT_AREA_FEMethod MydCTgc(mField,*projSurfaceCtx,"LAMBDA_CRACKFRONT_AREA");

  NL_Projection Projection(mField,*projSurfaceCtx,&myDirihletBC,"COUPLED_PROBLEM",arc_ctx);
 
  ierr = SNESSetApplicationContext(*snes,snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(*snes,F,SnesRhs,snes_ctx); CHKERRQ(ierr);
  if(material_FirelWall->operator[](FW_arc_lenhghat_definition)) {

    ierr = SNESSetJacobian(*snes,Arc_CTC_QTKQ,precK,SnesMat,snes_ctx); CHKERRQ(ierr);

    KSP ksp;
    ierr = SNESGetKSP(*snes,&ksp); CHKERRQ(ierr);
    PC pc;
    ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
    ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
    ierr = PCShellSetContext(pc,pc_ctx); CHKERRQ(ierr);
    ierr = PCShellSetApply(pc,pc_apply_arc_length); CHKERRQ(ierr);
    ierr = PCShellSetSetUp(pc,pc_setup_arc_length); CHKERRQ(ierr);

  } else {
    ierr = SNESSetJacobian(*snes,CTC_QTKQ,precK,SnesMat,snes_ctx); CHKERRQ(ierr);
  }
  ierr = SNESSetFromOptions(*snes); CHKERRQ(ierr);

  SnesCtx::loops_to_do_type& loops_to_do_Rhs = snes_ctx->get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ELASTIC_COUPLED",&MySpatialFE));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MATERIAL_COUPLED",&MyMaterialFE));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("dCT_CRACKFRONT_AREA_ELEM",&MydCTgc));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MESH_SMOOTHER",&MyMeshSmoother));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("FRONT_CONSTRAIN",&FrontPenalty));
  if(material_FirelWall->operator[](FW_arc_lenhghat_definition)) {
    loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",arc_elem));
  }
  SnesCtx::loops_to_do_type& loops_to_do_Mat = snes_ctx->get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ELASTIC_COUPLED",&MySpatialFE));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MATERIAL_COUPLED",&MyMaterialFE));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MESH_SMOOTHER",&MyMeshSmoother));

  snes_ctx->get_preProcess_to_do_Rhs().push_back(&Projection);
  snes_ctx->get_preProcess_to_do_Rhs().push_back(&MyCTgc);
  snes_ctx->get_postProcess_to_do_Rhs().push_back(&Projection);
  snes_ctx->get_preProcess_to_do_Mat().push_back(&Projection);
  snes_ctx->get_postProcess_to_do_Mat().push_back(&Projection);

  ierr = SNESSolve(*snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(*snes,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("COUPLED_PROBLEM",Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  if(material_FirelWall->operator[](FW_arc_lenhghat_definition)) {
    ierr = MySpatialFE.set_t_val(arc_ctx->get_FieldData()); CHKERRQ(ierr);
    aRea = arc_elem->aRea;
    lambda = arc_ctx->get_FieldData();
  }

  PostProcVertexMethod ent_method_spatial(mField.get_moab(),"SPATIAL_POSITION");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","SPATIAL_POSITION",Col,ent_method_spatial); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_material(mField.get_moab(),"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","MESH_NODE_POSITIONS",Col,ent_method_material); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_res_mat(mField.get_moab(),"MESH_NODE_POSITIONS",F,"MATERIAL_RESIDUAL");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","MESH_NODE_POSITIONS",Col,ent_method_res_mat); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_res_spat(mField.get_moab(),"SPATIAL_POSITION",F,"SPATIAL_RESIDUAL");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","SPATIAL_POSITION",Col,ent_method_res_spat); CHKERRQ(ierr);

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  MySpatialFE.snes_x = D;
  MySpatialFE.snes_f = F;
  MyMaterialFE.snes_x = D;
  MyMaterialFE.snes_f = F;
  MyCTgc.snes_x = D;
  MyCTgc.snes_f = F;
  Projection.snes_x = D;
  Projection.snes_f = F;
  Projection.set_snes_ctx(FieldInterface::SnesMethod::ctx_SNESSetFunction);
  MyCTgc.set_snes_ctx(FieldInterface::SnesMethod::ctx_SNESSetFunction);
  MySpatialFE.set_snes_ctx(FieldInterface::SnesMethod::ctx_SNESSetFunction);
  MyMaterialFE.set_snes_ctx(FieldInterface::SnesMethod::ctx_SNESSetFunction);
  ierr = mField.problem_basic_method_preProcess("COUPLED_PROBLEM",Projection); CHKERRQ(ierr);
  ierr = mField.problem_basic_method_preProcess("COUPLED_PROBLEM",MyCTgc); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("COUPLED_PROBLEM","ELASTIC_COUPLED",MySpatialFE);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("COUPLED_PROBLEM","MATERIAL_COUPLED",MyMaterialFE); CHKERRQ(ierr);
  ierr = mField.problem_basic_method_postProcess("COUPLED_PROBLEM",Projection); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_front_equlibrium(mField.get_moab(),"MESH_NODE_POSITIONS",F,"FRONT_EQULIBRIUM");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","MESH_NODE_POSITIONS",Col,ent_method_front_equlibrium); CHKERRQ(ierr);

  ierr = VecNorm(F,NORM_2,&nrm2_front_equlibrium); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"nrm2_Front_Equlibrium = %6.4e\n",nrm2_front_equlibrium); CHKERRQ(ierr);

  //detroy matrices
  ierr = MatDestroy(&precK); CHKERRQ(ierr);
  ierr = MatDestroy(&projSurfaceCtx->K); CHKERRQ(ierr);
  ierr = VecDestroy(&projSurfaceCtx->g); CHKERRQ(ierr);

  ierr = MatDestroy(&CTC_QTKQ); CHKERRQ(ierr);

  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);

  ierr = MatDestroy(&C_crack_fornt); CHKERRQ(ierr);

  ierr = delete_surface_projection_data(mField); CHKERRQ(ierr);
  ierr = delete_front_projection_data(mField); CHKERRQ(ierr);

  if(material_FirelWall->operator[](FW_arc_lenhghat_definition)) {
    ierr = MatDestroy(&Arc_CTC_QTKQ); CHKERRQ(ierr);
    delete arc_ctx;
    delete arc_snes_ctx;
    delete arc_mat_ctx;
    delete pc_ctx;
    delete arc_elem;
  } else {
    delete snes_ctx;
  }

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

  DirihletBCMethod_DriverComplexForLazy myDirihletBCMaterial(mField,problem,"MESH_NODE_POSITIONS");
  ierr = myDirihletBCMaterial.Init(); CHKERRQ(ierr);

  const double YoungModulus = 1;
  const double PoissonRatio = 0.;
  NL_MaterialFEMethod MyMaterialFE(mField,&myDirihletBCMaterial,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));

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
PetscErrorCode ConfigurationalFractureMechanics::ConstrainCrackForntEdges_FEMethod::preProcess() {
  PetscFunctionBegin;

    PetscErrorCode ierr;

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = VecDuplicate(snes_f,&tmp_snes_f); CHKERRQ(ierr);
	ierr = VecSwap(snes_f,tmp_snes_f); CHKERRQ(ierr);
	ierr = VecZeroEntries(snes_f); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(conf_prob->projSurfaceCtx->K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(conf_prob->projSurfaceCtx->K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	break;
      case ctx_SNESSetJacobian: 
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ConstrainCrackForntEdges_FEMethod::operator()() {
  PetscFunctionBegin;
  EntityHandle edge = fe_ptr->get_ent();

  ErrorCode rval;
  PetscErrorCode ierr;

  try {

  const EntityHandle* conn; 
  int num_nodes; 
  rval = mField.get_moab().get_connectivity(edge,conn,num_nodes,true); CHKERR_PETSC(rval);
  if(num_nodes!=2) SETERRQ(PETSC_COMM_SELF,1,"this implementation works for edges (2 nodes entities)");
  for(int nn = 0;nn<2; nn++) {
    for(_IT_GET_FEROW_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(this,"MESH_NODE_POSITIONS",conn[nn],dit)) { 
      rowDofs[3*nn+dit->get_EntDofIdx()] = dit->get_petsc_gloabl_dof_idx();
      dofsX[3*nn+dit->get_EntDofIdx()] = dit->get_FieldData();
    }
    for(_IT_GET_FECOL_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(this,"MESH_NODE_POSITIONS",conn[nn],dit)) { 
      colDofs[3*nn+dit->get_EntDofIdx()] = dit->get_petsc_gloabl_dof_idx();
    }
  }
  rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
  for(int dd = 0;dd<3;dd++) {
    delta0[dd] = coords[3+dd] - coords[dd];
    delta[dd] = dofsX[3+dd] - dofsX[dd];
  }
  double l0_2 = pow(delta0[0],2)+pow(delta0[1],2)+pow(delta0[2],2);
  double l0 = sqrt(l0_2);
  double l_2 = pow(delta[0],2)+pow(delta[1],2)+pow(delta[2],2);
  double e = (l_2-l0_2)/(2*l0_2);
  ublas::vector<FieldData> c = delta/l0;
  for(int dd = 0;dd<3;dd++) {
    f[dd+0] = -e*c[dd];
    f[dd+3] = +e*c[dd];
    for(int ddd = 0;ddd<3;ddd++) {
      K(dd+0,ddd+0) = +c[dd]*c[ddd] + e*( dd==ddd ? +1 : 0);
      K(dd+3,ddd+3) = +c[dd]*c[ddd] + e*( dd==ddd ? +1 : 0);
      K(dd+0,ddd+3) = -c[dd]*c[ddd] + e*( dd==ddd ? -1 : 0);
      K(dd+3,ddd+0) = -c[dd]*c[ddd] + e*( dd==ddd ? -1 : 0);
    }
  }
  K = K/l0;
  f = f*alpha3;
  K = K*alpha3;
  ierr = VecSetValues(snes_f,rowDofs.size(),&*rowDofs.data().begin(),&*f.data().begin(),ADD_VALUES); CHKERRQ(ierr);
  ierr = MatSetValues(conf_prob->projSurfaceCtx->K,
    rowDofs.size(),&*rowDofs.data().begin(),
    colDofs.size(),&*colDofs.data().begin(),
    &*K.data().begin(),ADD_VALUES); CHKERRQ(ierr);

  } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ConstrainCrackForntEdges_FEMethod::postProcess() {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  switch(snes_ctx) {
    case ctx_SNESSetFunction: {
      ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      //ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      PetscReal snes_f_nrm2;
      ierr = VecNorm(snes_f, NORM_2,&snes_f_nrm2); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"\ttruss f_nrm2 = %6.4e\n",snes_f_nrm2);
      ierr = VecAXPY(tmp_snes_f,1,snes_f); CHKERRQ(ierr);
      ierr = VecSwap(snes_f,tmp_snes_f); CHKERRQ(ierr);
      ierr = VecDestroy(&tmp_snes_f); CHKERRQ(ierr);
      ierr = MatAssemblyBegin(conf_prob->projSurfaceCtx->K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(conf_prob->projSurfaceCtx->K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    } break;
    case ctx_SNESSetJacobian: 
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }

  PetscFunctionReturn(0);
}
ConfigurationalFractureMechanics::ArcLengthElemFEMethod::ArcLengthElemFEMethod(
  FieldInterface& _mField,ConfigurationalFractureMechanics *_conf_prob,ArcLengthCtx *_arc_ptr): 
    mField(_mField),conf_prob(_conf_prob),arc_ptr(_arc_ptr) {
    PetscErrorCode ierr;
    ErrorCode rval;

    PetscInt ghosts[1] = { 0 };
    Interface &moab = mField.get_moab();
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm->rank() == 0) {
      ierr = VecCreateGhost(PETSC_COMM_WORLD,1,1,0,ghosts,&ghostDiag); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    } else {
      ierr = VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&ghostDiag); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }

    ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,crackSurfacesFaces,true); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    Range level_tris;
    ierr = mField.refine_get_ents(conf_prob->bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    crackSurfacesFaces = intersect(crackSurfacesFaces,level_tris);

    ierr = mField.get_problem("COUPLED_PROBLEM",&problem_ptr); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    set<DofIdx> set_idx;
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
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::calulate_lambda_int() {
  PetscFunctionBegin;
  ErrorCode rval;
  PetscErrorCode ierr;

  ierr = VecScatterCreate(snes_x,isSurface,surfaceDofs,PETSC_NULL,&surfaceScatter); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecScatterBegin(surfaceScatter,snes_x,surfaceDofs,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(surfaceScatter,snes_x,surfaceDofs,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterDestroy(&surfaceScatter); CHKERRABORT(PETSC_COMM_SELF,ierr);

  double *array;
  ierr = VecGetArray(surfaceDofs,&array); CHKERRQ(ierr);  
  PetscInt size;
  ISGetSize(isSurface,&size);
  for(int ii = 0;ii<size;ii++) {
    problem_ptr->numered_dofs_rows.get<
      PetscGlobalIdx_mi_tag>().find(isIdx[ii])->get_FieldData() = array[ii];
  }
  ierr = VecRestoreArray(surfaceDofs,&array); CHKERRABORT(PETSC_COMM_SELF,ierr);

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
	ierr = MatAssemblyBegin(conf_prob->projSurfaceCtx->K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(conf_prob->projSurfaceCtx->K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = get_dlambda(snes_x); CHKERRQ(ierr);
	ierr = calulate_lambda_int(); CHKERRQ(ierr);
	ierr = calulate_db(); CHKERRQ(ierr);
      }
      break;
    case ctx_SNESSetJacobian: 
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
	//calulate diagonal therm
	double diag = arc_ptr->beta*sqrt(arc_ptr->F_lambda2);
	ierr = VecSetValue(ghostDiag,0,diag,INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(conf_prob->projSurfaceCtx->K,arc_ptr->get_petsc_gloabl_dof_idx(),arc_ptr->get_petsc_gloabl_dof_idx(),1,ADD_VALUES); CHKERRQ(ierr);
      }
      break; 
      case ctx_SNESSetJacobian:       
      break;
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
	ierr = VecAssemblyBegin(ghostDiag); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(ghostDiag); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(ghostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(ghostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	double *diag;
	ierr = VecGetArray(ghostDiag,&diag); CHKERRQ(ierr);
	arc_ptr->diag = *diag;
	ierr = VecRestoreArray(ghostDiag,&diag); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tdiag = %6.4e\n",arc_ptr->diag);
      }
      break;
      case ctx_SNESSetJacobian: 
      break;
      default:
      break;
    }
    PetscFunctionReturn(0);
}

