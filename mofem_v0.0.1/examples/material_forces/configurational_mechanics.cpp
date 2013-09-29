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

#include "configurational_mechanics.hpp"

#include "moabField_Core.hpp"
#include <petscksp.h>

#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include "moabFEMethod_DriverComplexForLazy.hpp"

using namespace MoFEM;

struct NL_ElasticFEMethod: public FEMethod_DriverComplexForLazy_Spatial {
  
    NL_ElasticFEMethod(moabField& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0): 
        FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
        FEMethod_DriverComplexForLazy_Spatial(_mField,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose)  {
      set_PhysicalEquationNumber(stvenant_kirchhoff);
    }
  
  };
  
struct NL_MaterialFEMethod: public FEMethod_DriverComplexForLazy_Material {
  
    NL_MaterialFEMethod(moabField& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0): 
        FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
        FEMethod_DriverComplexForLazy_Material(_mField,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose)  {
      set_PhysicalEquationNumber(stvenant_kirchhoff);
    }
  
  };
  
struct NL_MaterialFEMethodProjected: public FEMethod_DriverComplexForLazy_MaterialProjected {
  
    NL_MaterialFEMethodProjected(moabField& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0): 
        FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
        FEMethod_DriverComplexForLazy_MaterialProjected(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose)  {
      set_PhysicalEquationNumber(stvenant_kirchhoff);
    }
  
  };
  
struct NL_ElasticFEMethodCoupled: public FEMethod_DriverComplexForLazy_CoupledSpatial {
  
    NL_ElasticFEMethodCoupled(moabField& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0):
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
      FEMethod_DriverComplexForLazy_CoupledSpatial(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose) {
        set_PhysicalEquationNumber(stvenant_kirchhoff);
      }
  
  };
  
struct NL_MaterialFEMethodCoupled: public FEMethod_DriverComplexForLazy_CoupledMaterial {
  
    NL_MaterialFEMethodCoupled(moabField& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0):
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
      FEMethod_DriverComplexForLazy_CoupledMaterial(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose) {
        set_PhysicalEquationNumber(stvenant_kirchhoff);
      }
  
  };
  
struct NL_MeshSmootherCoupled: public FEMethod_DriverComplexForLazy_CoupledMeshSmoother {
  
    NL_MeshSmootherCoupled(moabField& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,int _verbose = 0):
      FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
      FEMethod_DriverComplexForLazy_CoupledMeshSmoother(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr) {
      set_qual_ver(0);
      }
  
  };
  
struct CubitDisplacementDirihletBC_Coupled: public CubitDisplacementDirihletBC {
  
    CubitDisplacementDirihletBC_Coupled (moabField& _mField,const string _problem_name): 
      CubitDisplacementDirihletBC(_mField,_problem_name,"None") {}
  
    PetscErrorCode SetDirihletBC_to_ElementIndiciesRow(
      moabField::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirihletBC) {
      PetscFunctionBegin;
      field_name = "SPATIAL_POSITION";
      ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndiciesRow(fe_method_ptr,RowGlobDofs,DirihletBC); CHKERRQ(ierr);
      field_name = "MESH_NODE_POSITIONS";
      ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndiciesRow(fe_method_ptr,RowGlobDofs,DirihletBC); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  
    PetscErrorCode SetDirihletBC_to_ElementIndiciesCol(
      moabField::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC) {
      PetscFunctionBegin;
      CubitDisplacementDirihletBC::field_name = "SPATIAL_POSITION";
      ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndiciesCol(fe_method_ptr,ColGlobDofs,DirihletBC); CHKERRQ(ierr);
      CubitDisplacementDirihletBC::field_name = "MESH_NODE_POSITIONS";
      ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndiciesCol(fe_method_ptr,ColGlobDofs,DirihletBC); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  
    PetscErrorCode SetDirihletBC_to_ElementIndicies(
      moabField::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC) {
      PetscFunctionBegin;
      field_name = "SPATIAL_POSITION";
      ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndicies(fe_method_ptr,RowGlobDofs,ColGlobDofs,DirihletBC); CHKERRQ(ierr);
      field_name = "MESH_NODE_POSITIONS";
      ierr = CubitDisplacementDirihletBC::SetDirihletBC_to_ElementIndicies(fe_method_ptr,RowGlobDofs,ColGlobDofs,DirihletBC); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  
  };
  
PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_SetMaterialFireWall(moabField& mField) {
  PetscFunctionBegin;

  ErrorCode rval;

  BitRefLevel def_bit_level = 0;
  rval = mField.get_moab().tag_get_handle("_Materiar_FireWall",sizeof(Material_FirelWall_def),MB_TYPE_OPAQUE,
    th_MaterialFireWall,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level); 
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  rval = mField.get_moab().tag_get_by_ptr(th_MaterialFireWall,&root_meshset,1,(const void**)&Material_FirelWall); CHKERR_PETSC(rval);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_SpatialProblemDefinition(moabField& mField) {
  PetscFunctionBegin;

  if(Material_FirelWall->operator[](0)) PetscFunctionReturn(0);
  Material_FirelWall->set(0);

  PetscErrorCode ierr;

  //Fields
  ierr = mField.add_field("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);


  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);

  //add entitities (by tets) to the field
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
  ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_MaterialProblemDefinition(moabField& mField) {
  PetscFunctionBegin;

  if(Material_FirelWall->operator[](1)) PetscFunctionReturn(0);
  Material_FirelWall->set(1);

  PetscErrorCode ierr;

  //Fields
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("MATERIAL_FORCE",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("MATERIAL"); CHKERRQ(ierr);

  //Define rows/cols and element data 
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //
  ierr = mField.modify_finite_element_add_field_row("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
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

PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_CoupledProblemDefinition(moabField& mField) {
  PetscFunctionBegin;

  if(Material_FirelWall->operator[](2)) PetscFunctionReturn(0);
  Material_FirelWall->set(3);

  PetscErrorCode ierr;

  //Fields
  ierr = mField.add_field("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("MATERIAL_FORCE",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("MATERIAL"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("MESH_SMOOTHER"); CHKERRQ(ierr);
  
  //fes definitions
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //
  ierr = mField.modify_finite_element_add_field_row("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("MATERIAL","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MATERIAL","SPATIAL_POSITION"); CHKERRQ(ierr);
  //
  ierr = mField.modify_finite_element_add_field_row("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("COUPLED_PROBLEM"); CHKERRQ(ierr);
  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","MATERIAL"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_PROBLEM","MESH_SMOOTHER"); CHKERRQ(ierr);

  //add entitities (by tets) to the field
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
  ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);


  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_ConstrainsProblemDefinition(moabField& mField) {
  PetscFunctionBegin;

  if(Material_FirelWall->operator[](4)) PetscFunctionReturn(0);
  Material_FirelWall->set(4);

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

  //add tets on corners
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
    //
    Range SurfacesTets;
    Interface& moab = mField.get_moab();
    rval = moab.get_adjacencies(SurfacesFaces,3,false,SurfacesTets,Interface::UNION); CHKERR_PETSC(rval);
    Range CrackSurfacesTets;
    rval = moab.get_adjacencies(CrackSurfacesFaces,3,false,CrackSurfacesTets,Interface::UNION); CHKERR_PETSC(rval);
    CrackSurfacesTets = CrackSurfacesTets.subset_by_type(MBTET);
    {
      Range CornersEdgesNodes;
      rval = moab.get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
      CornersNodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
      //Range CrackCornersEdgesNodes;
      //rval = moab.get_connectivity(CrackCornersEdges,CrackCornersEdgesNodes,true); CHKERR_PETSC(rval);
      //CornersNodes.insert(CrackCornersEdgesNodes.begin(),CrackCornersEdgesNodes.end());
      rval = moab.create_meshset(MESHSET_SET,CornersNodesMeshset); CHKERR_PETSC(rval);	
      rval = moab.add_entities(CornersNodesMeshset,CornersNodes); CHKERR_PETSC(rval);
      //add surface elements
      Range CornersTets;
      rval = moab.get_adjacencies(CornersNodes,3,false,CornersTets,Interface::UNION); CHKERR_PETSC(rval);
      EntityHandle CornersTetsMeshset;
      rval = moab.create_meshset(MESHSET_SET,CornersTetsMeshset); CHKERR_PETSC(rval);	
      rval = moab.add_entities(CornersTetsMeshset,CornersTets); CHKERR_PETSC(rval);
      CornersTets = intersect(CornersTets,SurfacesTets);
      ierr = mField.add_ents_to_finite_element_by_TETs(CornersTetsMeshset,"C_CORNER_ELEM"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TETs(CornersTetsMeshset,"CTC_CORNER_ELEM"); CHKERRQ(ierr);
      rval = moab.delete_entities(&CornersTetsMeshset,1); CHKERR_PETSC(rval);
    }
    {
      Range SurfacesNodes;
      rval = moab.get_adjacencies(SurfacesFaces,0,false,SurfacesNodes,Interface::UNION); CHKERR_PETSC(rval);
      Range CornersEdgesNodes;
      rval = moab.get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
      //Range CrackCornersEdgesNodes;
      //rval = moab.get_connectivity(CrackCornersEdges,CrackCornersEdgesNodes,true); CHKERR_PETSC(rval);
      SurfacesNodes = subtract(SurfacesNodes,CornersNodes);
      SurfacesNodes = subtract(SurfacesNodes,CornersEdgesNodes);
      //SurfacesNodes = subtract(SurfacesNodes,CrackCornersEdgesNodes);
      rval = moab.create_meshset(MESHSET_SET,SurfacesFacesMeshset); CHKERR_PETSC(rval);	
      rval = moab.add_entities(SurfacesFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
      rval = moab.add_entities(SurfacesFacesMeshset,SurfacesNodes); CHKERR_PETSC(rval);
      //add surface elements
      EntityHandle SurfacesTetsMeshset;
      rval = moab.create_meshset(MESHSET_SET,SurfacesTetsMeshset); CHKERR_PETSC(rval);	
      rval = moab.add_entities(SurfacesTetsMeshset,SurfacesTets); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TETs(SurfacesTetsMeshset,"C_SURFACE_ELEM"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TETs(SurfacesTetsMeshset,"CTC_SURFACE_ELEM"); CHKERRQ(ierr);
      rval = moab.delete_entities(&SurfacesTetsMeshset,1); CHKERR_PETSC(rval);
    }
    //CRCAK
    if(cs) {
      {
	Range CrackSurfacesNodes;
	rval = moab.get_adjacencies(CrackSurfacesFaces,0,false,CrackSurfacesNodes,Interface::UNION); CHKERR_PETSC(rval);
	Range CornersEdgesNodes;
	rval = moab.get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
	Range CrackCornersEdgesNodes;
	rval = moab.get_connectivity(CrackCornersEdges,CrackCornersEdgesNodes,true); CHKERR_PETSC(rval);
	CrackSurfacesNodes = subtract(CrackSurfacesNodes,CornersNodes);
	CrackSurfacesNodes = subtract(CrackSurfacesNodes,CornersEdgesNodes);
	CrackSurfacesNodes = subtract(CrackSurfacesNodes,CrackCornersEdgesNodes);
	rval = moab.create_meshset(MESHSET_SET,CrackSurfacesFacesMeshset); CHKERR_PETSC(rval);	
	rval = moab.add_entities(CrackSurfacesFacesMeshset,CrackSurfacesFaces); CHKERR_PETSC(rval);
	rval = moab.add_entities(CrackSurfacesFacesMeshset,CrackSurfacesNodes); CHKERR_PETSC(rval);
	//add surface elements
	EntityHandle CrackSurfacesTetsMeshset;
	rval = moab.create_meshset(MESHSET_SET,CrackSurfacesTetsMeshset); CHKERR_PETSC(rval);	
	rval = moab.add_entities(CrackSurfacesTetsMeshset,CrackSurfacesTets); CHKERR_PETSC(rval);
	ierr = mField.add_ents_to_finite_element_by_TETs(CrackSurfacesTetsMeshset,"C_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
	ierr = mField.add_ents_to_finite_element_by_TETs(CrackSurfacesTetsMeshset,"CTC_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
	rval = moab.delete_entities(&CrackSurfacesTetsMeshset,1); CHKERR_PETSC(rval);
      }
    }
  }

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_VERTICEs(SurfacesFacesMeshset,"LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_VERTICEs(CornersNodesMeshset,"LAMBDA_CORNER"); CHKERRQ(ierr);
  //CRACK
  if(cs) {
    ierr = mField.add_ents_to_field_by_VERTICEs(CrackSurfacesFacesMeshset,"LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
  }

  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_SURFACE",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CORNER",1); CHKERRQ(ierr);
  if(cs) {
    ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CRACK_SURFACE",1); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}


PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_SpatialPartitionProblems(moabField& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_MaterialPartitionProblems(moabField& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = mField.partition_problem("MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("MATERIAL_MECHANICS"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_CoupledPartitionProblems(moabField& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = mField.partition_problem("COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("COUPLED_PROBLEM"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_ConstrainsPartitionProblems(moabField& mField,string problem) {
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

PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_SetSpatialPositions(moabField& mField) {
  PetscFunctionBegin;
  if(Material_FirelWall->operator[](5)) PetscFunctionReturn(0);
  Material_FirelWall->set(5);

  ErrorCode rval;

  EntityHandle node = 0;
  double coords[3];
  for(_IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP_(mField,"SPATIAL_POSITION",dof_ptr)) {
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
PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_SetMaterialPositions(moabField& mField) {
  PetscFunctionBegin;

  if(Material_FirelWall->operator[](6)) PetscFunctionReturn(0);
  Material_FirelWall->set(6);

  ErrorCode rval;

  EntityHandle node = 0;
  double coords[3];
  for(_IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP_(mField,"MESH_NODE_POSITIONS",dof_ptr)) {
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

PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_SolveSpatialProblem(moabField& mField) {
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
  const double PoissonRatio = 0.25;
  NL_ElasticFEMethod MySpatialFE(mField,&myDirihletBCSpatial,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));

  moabSnesCtx SnesCtx(mField,"ELASTIC_MECHANICS");
  
  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,Aij,Aij,SnesMat,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  moabSnesCtx::loops_to_do_type& loops_to_do_Rhs = SnesCtx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MySpatialFE));
  moabSnesCtx::loops_to_do_type& loops_to_do_Mat = SnesCtx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MySpatialFE));

  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = mField.set_local_VecCreateGhost("ELASTIC_MECHANICS",Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  double step_size = -1e-3;
  for(int step = 1;step<2; step++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Load Setp %D\n",step); CHKERRQ(ierr);
    ierr = MySpatialFE.set_t_val(step_size*step); CHKERRQ(ierr);
    ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
    int its;
    ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  }
  
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PostProcVertexMethod ent_method(mField.get_moab(),"SPATIAL_POSITION");
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Col,ent_method); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_res(mField.get_moab(),"SPATIAL_POSITION",F,"PHYSICAL_RESIDUAL");
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Col,ent_method_res); CHKERRQ(ierr);

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_ProjectForceVector(moabField& mField,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;

  Range CornersEdges,CornersNodes,SurfacesFaces,CrackSurfacesFaces,CrackCornersEdges;

  ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,CornersEdges,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,CornersNodes,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(102,SideSet,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,CrackSurfacesFaces,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,CrackCornersEdges,true); CHKERRQ(ierr);

  Interface& moab = mField.get_moab();

  Range CornersEdgesNodes;
  rval = moab.get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
  CornersNodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
  //Range CrackCornersEdgesNodes;
  //rval = moab.get_connectivity(CrackCornersEdges,CrackCornersEdgesNodes,true); CHKERR_PETSC(rval);
  //CornersNodes.insert(CrackCornersEdgesNodes.begin(),CrackCornersEdgesNodes.end());

  matPROJ_ctx proj_all_ctx(mField,problem,"C_ALL_MATRIX");
  ierr = mField.MatCreateMPIAIJWithArrays("C_ALL_MATRIX",&proj_all_ctx.C); CHKERRQ(ierr);

  C_SURFACE_FEMethod CFE_SURFACE(moab,SurfacesFaces,proj_all_ctx.C);
  C_SURFACE_FEMethod CFE_CRACK_SURFACE(moab,CrackSurfacesFaces,proj_all_ctx.C,"LAMBDA_CRACK_SURFACE");
  C_CORNER_FEMethod CFE_CORNER(moab,CornersNodes,proj_all_ctx.C);

  ierr = MatSetOption(proj_all_ctx.C,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(proj_all_ctx.C,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

  ierr = MatZeroEntries(proj_all_ctx.C); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_SURFACE_ELEM",CFE_SURFACE);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_CRACK_SURFACE_ELEM",CFE_CRACK_SURFACE);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(proj_all_ctx.C,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(proj_all_ctx.C,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_CORNER_ELEM",CFE_CORNER);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(proj_all_ctx.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(proj_all_ctx.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  {
    //Matrix View
    MatView(proj_all_ctx.C,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
    std::string wait;
    std::cin >> wait;
  }
  
  Vec F_Material;
  ierr = mField.VecCreateGhost(problem,Col,&F_Material); CHKERRQ(ierr);
  /*for(_IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP_(mField,"MATERIAL_FORCE",dof)) {
    cerr << *dof << endl;
  }*/
  ierr = mField.set_other_global_VecCreateGhost(problem,"MESH_NODE_POSITIONS","MATERIAL_FORCE",Row,F_Material,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = proj_all_ctx.InitQorP(F_Material); CHKERRQ(ierr);

  int M,m;
  ierr = VecGetSize(F_Material,&M); CHKERRQ(ierr);
  ierr = VecGetLocalSize(F_Material,&m); CHKERRQ(ierr);
  Mat Q;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,&proj_all_ctx,&Q); CHKERRQ(ierr);
  ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))matQ_mult_shell); CHKERRQ(ierr);

  Vec QTF_Material;
  ierr = VecDuplicate(F_Material,&QTF_Material); CHKERRQ(ierr);
  ierr = MatMult(Q,F_Material,QTF_Material); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_material(mField.get_moab(),"MESH_NODE_POSITIONS",QTF_Material,"MATERIAL_FORCE_PROJECTED");
  ierr = mField.loop_dofs(problem,"MESH_NODE_POSITIONS",Col,ent_method_material); CHKERRQ(ierr);

  ierr = VecDestroy(&QTF_Material); CHKERRQ(ierr);
  ierr = VecDestroy(&F_Material); CHKERRQ(ierr);
  ierr = MatDestroy(&proj_all_ctx.C); CHKERRQ(ierr);
  ierr = MatDestroy(&Q); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_SolveMaterialProblem(moabField& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  DirihletBCMethod_DriverComplexForLazy myDirihletBCMaterial(mField,"MATERIAL_MECHANICS","MESH_POSITION");
  ierr = myDirihletBCMaterial.Init(); CHKERRQ(ierr);

  matPROJ_ctx proj_all_ctx(mField,"MATERIAL_MECHANICS","C_ALL_MATRIX");
  Mat precK;
  ierr = mField.MatCreateMPIAIJWithArrays("MATERIAL_MECHANICS",&precK); CHKERRQ(ierr);
  ierr = MatDuplicate(precK,MAT_DO_NOT_COPY_VALUES,&proj_all_ctx.K); CHKERRQ(ierr);
  ierr = mField.MatCreateMPIAIJWithArrays("C_ALL_MATRIX",&proj_all_ctx.C); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("C_ALL_MATRIX",Row,&proj_all_ctx.g); CHKERRQ(ierr);

  Vec F;
  ierr = mField.VecCreateGhost("MATERIAL_MECHANICS",Row,&F); CHKERRQ(ierr);


  Mat CTC_QTKQ;
  int M,N,m,n;
  ierr = MatGetSize(proj_all_ctx.K,&M,&N); CHKERRQ(ierr);
  ierr = MatGetLocalSize(proj_all_ctx.K,&m,&n); CHKERRQ(ierr);
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,&proj_all_ctx,&CTC_QTKQ); CHKERRQ(ierr);
  ierr = MatShellSetOperation(CTC_QTKQ,MATOP_MULT,(void(*)(void))matCTC_QTKQ_mult_shell); CHKERRQ(ierr);

  const double YoungModulus = 1;
  const double PoissonRatio = 0.25;
  NL_MaterialFEMethodProjected MyMaterialFE(mField,proj_all_ctx,&myDirihletBCMaterial,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));

  moabSnesCtx SnesCtx(mField,"MATERIAL_MECHANICS");
  
  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,CTC_QTKQ,precK,SnesMat,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);


  moabSnesCtx::loops_to_do_type& loops_to_do_Rhs = SnesCtx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("MATERIAL",&MyMaterialFE));
  moabSnesCtx::loops_to_do_type& loops_to_do_Mat = SnesCtx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("MATERIAL",&MyMaterialFE));

  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = mField.set_local_VecCreateGhost("MATERIAL_MECHANICS",Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  PostProcVertexMethod ent_method(mField.get_moab(),"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MATERIAL_MECHANICS","MESH_NODE_POSITIONS",Col,ent_method); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_res(mField.get_moab(),"MESH_NODE_POSITIONS",F,"MATERIAL_RESIDUAL");
  ierr = mField.loop_dofs("MATERIAL_MECHANICS","MESH_NODE_POSITIONS",Col,ent_method_res); CHKERRQ(ierr);


  ierr = proj_all_ctx.DestroyQorP(); CHKERRQ(ierr);
  ierr = proj_all_ctx.DestroyQTKQ(); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&proj_all_ctx.K); CHKERRQ(ierr);
  ierr = MatDestroy(&precK); CHKERRQ(ierr);
  ierr = MatDestroy(&proj_all_ctx.C); CHKERRQ(ierr);
  ierr = VecDestroy(&proj_all_ctx.g); CHKERRQ(ierr);
  ierr = MatDestroy(&CTC_QTKQ); CHKERRQ(ierr);

  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_CalculateSpatialResidual(moabField& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  DirihletBCMethod_DriverComplexForLazy myDirihletBCSpatial(mField,"ELASTIC_MECHANICS","SPATIAL_POSITION");
  ierr = myDirihletBCSpatial.Init(); CHKERRQ(ierr);

  const double YoungModulus = 1;
  const double PoissonRatio = 0.25;
  NL_ElasticFEMethod MySpatialFE(mField,&myDirihletBCSpatial,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  Vec F_Spatial;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&F_Spatial); CHKERRQ(ierr);
  ierr = VecZeroEntries(F_Spatial); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_Spatial,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Spatial,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  MySpatialFE.snes_f = F_Spatial;
  MySpatialFE.set_snes_ctx(moabField::SnesMethod::ctx_SNESSetFunction);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MySpatialFE);  CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_Spatial,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Spatial,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F_Spatial); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F_Spatial); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_res_spatial(mField.get_moab(),"SPATIAL_POSITION",F_Spatial,"PHYSICAL_RESIDUAL");
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Col,ent_method_res_spatial); CHKERRQ(ierr);

  ierr = VecDestroy(&F_Spatial); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_CalculateMaterialForces(moabField& mField,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  DirihletBCMethod_DriverComplexForLazy myDirihletBCMaterial(mField,problem,"MESH_NODE_POSITIONS");
  ierr = myDirihletBCMaterial.Init(); CHKERRQ(ierr);

  const double YoungModulus = 1;
  const double PoissonRatio = 0.25;
  NL_MaterialFEMethod MyMaterialFE(mField,&myDirihletBCMaterial,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));

  Vec F_Material;
  ierr = mField.VecCreateGhost(problem,Col,&F_Material); CHKERRQ(ierr);
  ierr = VecZeroEntries(F_Material); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_Material,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Material,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  MyMaterialFE.snes_f = F_Material;
  MyMaterialFE.set_snes_ctx(moabField::SnesMethod::ctx_SNESSetFunction);
  ierr = mField.loop_finite_elements(problem,"MATERIAL",MyMaterialFE);  CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_Material,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Material,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F_Material); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F_Material); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_material(mField.get_moab(),"MESH_NODE_POSITIONS",F_Material,"MATERIAL_FORCE");
  ierr = mField.loop_dofs(problem,"MESH_NODE_POSITIONS",Col,ent_method_material); CHKERRQ(ierr);

  //Fields
  ierr = mField.set_other_global_VecCreateGhost(problem,"MESH_NODE_POSITIONS","MATERIAL_FORCE",Row,F_Material,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //for(_IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP_(mField,"MATERIAL_FORCE",dof)) {
    //cerr << *dof << endl;
  //}


  //detroy matrices
  ierr = VecDestroy(&F_Material); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalMechanics::ConfigurationalMechanics_SolveCoupledProblem(moabField& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  matPROJ_ctx proj_all_ctx(mField,"COUPLED_PROBLEM","C_ALL_MATRIX");
  Mat precK;
  ierr = mField.MatCreateMPIAIJWithArrays("COUPLED_PROBLEM",&precK); CHKERRQ(ierr);
  ierr = MatDuplicate(precK,MAT_DO_NOT_COPY_VALUES,&proj_all_ctx.K); CHKERRQ(ierr);
  ierr = mField.MatCreateMPIAIJWithArrays("C_ALL_MATRIX",&proj_all_ctx.C); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("C_ALL_MATRIX",Row,&proj_all_ctx.g); CHKERRQ(ierr);

  //create matrices
  Mat CTC_QTKQ;
  int M,N,m,n;
  ierr = MatGetSize(proj_all_ctx.K,&M,&N); CHKERRQ(ierr);
  ierr = MatGetLocalSize(proj_all_ctx.K,&m,&n); CHKERRQ(ierr);
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,&proj_all_ctx,&CTC_QTKQ); CHKERRQ(ierr);
  ierr = MatShellSetOperation(CTC_QTKQ,MATOP_MULT,(void(*)(void))matCTC_QTKQ_mult_shell); CHKERRQ(ierr);
  Vec F;
  ierr = mField.VecCreateGhost("COUPLED_PROBLEM",Row,&F); CHKERRQ(ierr);

  CubitDisplacementDirihletBC_Coupled myDirihletBC(mField,"COUPLED_PROBLEM");
  ierr = myDirihletBC.Init(); CHKERRQ(ierr);

  const double YoungModulus = 1;
  const double PoissonRatio = 0.25;
  NL_ElasticFEMethodCoupled MySpatialFE(mField,proj_all_ctx,&myDirihletBC,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  NL_MaterialFEMethodCoupled MyMaterialFE(mField,proj_all_ctx,&myDirihletBC,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  NL_MeshSmootherCoupled MyMeshSmoother(mField,proj_all_ctx,&myDirihletBC);
  FEMethod_DriverComplexForLazy_CoupledProjected Projection(mField,proj_all_ctx,&myDirihletBC,"COUPLED_PROBLEM");

  moabSnesCtx SnesCtx(mField,"COUPLED_PROBLEM");
  
  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,CTC_QTKQ,precK,SnesMat,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  moabSnesCtx::loops_to_do_type& loops_to_do_Rhs = SnesCtx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MySpatialFE));
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("MATERIAL",&MyMaterialFE));
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("MESH_SMOOTHER",&MyMeshSmoother));
  moabSnesCtx::loops_to_do_type& loops_to_do_Mat = SnesCtx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MySpatialFE));
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("MATERIAL",&MyMaterialFE));
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("MESH_SMOOTHER",&MyMeshSmoother));

  SnesCtx.get_preProcess_to_do_Rhs().push_back(&Projection);
  SnesCtx.get_postProcess_to_do_Rhs().push_back(&Projection);
  SnesCtx.get_preProcess_to_do_Mat().push_back(&Projection);
  SnesCtx.get_postProcess_to_do_Mat().push_back(&Projection);

  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = mField.set_local_VecCreateGhost("COUPLED_PROBLEM",Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  double step_size = -1e-3;
  for(int step = 1;step<2; step++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Load Setp %D\n",step); CHKERRQ(ierr);
    ierr = MySpatialFE.set_t_val(step_size*step); CHKERRQ(ierr);
    ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
    int its;
    ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  }
  
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("COUPLED_PROBLEM",Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_spatial(mField.get_moab(),"SPATIAL_POSITION");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","SPATIAL_POSITION",Col,ent_method_spatial); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_material(mField.get_moab(),"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","MESH_NODE_POSITIONS",Col,ent_method_material); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_res_mat(mField.get_moab(),"MESH_NODE_POSITIONS",F,"MATERIAL_RESIDUAL");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","MESH_NODE_POSITIONS",Col,ent_method_res_mat); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_res_spat(mField.get_moab(),"SPATIAL_POSITION",F,"PHYSICAL_RESIDUAL");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","SPATIAL_POSITION",Col,ent_method_res_spat); CHKERRQ(ierr);

  //detroy matrices
  ierr = proj_all_ctx.DestroyQorP(); CHKERRQ(ierr);
  ierr = proj_all_ctx.DestroyQTKQ(); CHKERRQ(ierr);
  ierr = MatDestroy(&proj_all_ctx.K); CHKERRQ(ierr);
  ierr = MatDestroy(&precK); CHKERRQ(ierr);
  ierr = MatDestroy(&proj_all_ctx.C); CHKERRQ(ierr);
  ierr = VecDestroy(&proj_all_ctx.g); CHKERRQ(ierr);
  ierr = MatDestroy(&CTC_QTKQ); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&CTC_QTKQ); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


