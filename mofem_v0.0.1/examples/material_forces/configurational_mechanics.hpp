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
#include <petscksp.h>

#include "moabSnes.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include "moabFEMethod_DriverComplexForLazy.hpp"

using namespace MoFEM;

struct NL_ElasticFEMethod: public FEMethod_DriverComplexForLazy_Spatial {

  NL_ElasticFEMethod(moabField& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0): 
      FEMethod_DriverComplexForLazy_Spatial(_mField,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose)  {
    set_PhysicalEquationNumber(hooke);
  }

};

struct NL_MaterialFEMethod: public FEMethod_DriverComplexForLazy_Material {

  NL_MaterialFEMethod(moabField& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0): 
      FEMethod_DriverComplexForLazy_Material(_mField,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose)  {
    set_PhysicalEquationNumber(hooke);
  }

};

PetscErrorCode ConfigurationalMechanics_PhysicalProblemDefinition(moabField& mField) {
  PetscFunctionBegin;

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


PetscErrorCode ConfigurationalMechanics_MaterialProblemDefinition(moabField& mField) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  //Fields
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("MATERIAL"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
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

PetscErrorCode ConfigurationalMechanics_ConstrainsProblemDefinition(moabField& mField) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  //Fields
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("C_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CTC_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("C_CORNER_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CTC_CORNER_ELEM"); CHKERRQ(ierr);

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

  //define problems
  ierr = mField.add_problem("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.add_problem("C_ALL_MATRIX"); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_CORNER_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("C_ALL_MATRIX","C_CORNER_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("C_ALL_MATRIX","C_SURFACE_ELEM"); CHKERRQ(ierr);

  //add tets on corners
  EntityHandle CornersNodesMeshset,SurfacesFacesMeshset;
  {
    Range CornersEdges,CornersNodes,SurfacesFaces;
    ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,CornersEdges,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,CornersNodes,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(102,SideSet,2,SurfacesFaces,true); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 100 = %d\n",CornersEdges.size()); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of NodeSet 101 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 102 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
    Range SurfacesTets;
    Interface& moab = mField.get_moab();
    rval = moab.get_adjacencies(SurfacesFaces,3,false,SurfacesTets,Interface::UNION); CHKERR_PETSC(rval);
    {
      Range CornersEdgesNodes;
      rval = moab.get_adjacencies(CornersEdges,0,false,CornersEdgesNodes,Interface::UNION); CHKERR_PETSC(rval);
      rval = moab.create_meshset(MESHSET_SET,CornersNodesMeshset); CHKERR_PETSC(rval);	
      CornersNodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
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
      SurfacesNodes = subtract(SurfacesNodes,CornersNodes);
      SurfacesNodes = subtract(SurfacesNodes,CornersEdgesNodes);
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
  }

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_VERTICEs(SurfacesFacesMeshset,"LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_VERTICEs(CornersNodesMeshset,"LAMBDA_CORNER"); CHKERRQ(ierr);

  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_SURFACE",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CORNER",1); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode ConfigurationalMechanics_PhysicalPartitionProblems(moabField& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalMechanics_MaterialPartitionProblems(moabField& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = mField.partition_problem("MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("MATERIAL_MECHANICS"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalMechanics_ConstrainsPartitionProblems(moabField& mField) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = mField.partition_problem("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  //partition
  ierr = mField.compose_problem("C_ALL_MATRIX","CCT_ALL_MATRIX",false,"MATERIAL_MECHANICS",true); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("C_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("C_ALL_MATRIX"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

