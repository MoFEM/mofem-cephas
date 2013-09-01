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

#include "hook_material_forces.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";


int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
 
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  moabField_Core core(moab);
  moabField& mField = core;

  Range CubitSideSets_meshsets;
  ierr = mField.get_CubitBCType_meshsets(SideSet,CubitSideSets_meshsets); CHKERRQ(ierr);

  Range CornersEdges,CornersNodes,SurfacesFaces;
  ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,CornersEdges,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,CornersNodes,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(102,SideSet,2,SurfacesFaces,true); CHKERRQ(ierr);

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_field("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("MATERIAL_FORCE",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("LAMBDA_SURFACE",H1,1); CHKERRQ(ierr);
  ierr = mField.add_field("LAMBDA_EDGE",H1,2); CHKERRQ(ierr);
  ierr = mField.add_field("LAMBDA_CORNER",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("MATERIAL"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("C_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CTC_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("C_EDGE_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CTC_EDGE_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("C_CORNER_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CTC_CORNER_ELEM"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MATERIAL","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("C_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("C_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("C_EDGE_ELEM","LAMBDA_EDGE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("C_EDGE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_EDGE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_EDGE_ELEM","LAMBDA_EDGE"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("CTC_EDGE_ELEM","LAMBDA_EDGE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("CTC_EDGE_ELEM","LAMBDA_EDGE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("CTC_EDGE_ELEM","LAMBDA_EDGE"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("C_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("C_CORNER_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_CORNER_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("CTC_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("CTC_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("CTC_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MATERIAL",MBTET); CHKERRQ(ierr);

  //add tets on corners
  EntityHandle CornersNodesMeshset;
  Range SurfacesTets;
  rval = moab.get_adjacencies(SurfacesFaces,3,false,SurfacesTets,Interface::UNION); CHKERR_PETSC(rval);
  {
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
  //add tets on surface
  EntityHandle CornersEdgesMeshset;
  {
    Range CornersEdgesNodes;
    rval = moab.get_adjacencies(CornersEdges,0,false,CornersEdgesNodes,Interface::UNION); CHKERR_PETSC(rval);
    //CornersEdgesNodes = subtract(CornersEdgesNodes,CornersNodes);
    rval = moab.create_meshset(MESHSET_SET,CornersEdgesMeshset); CHKERR_PETSC(rval);	
    rval = moab.add_entities(CornersEdgesMeshset,CornersEdges); CHKERR_PETSC(rval);
    rval = moab.add_entities(CornersEdgesMeshset,CornersEdgesNodes); CHKERR_PETSC(rval);
    rval = moab.add_entities(CornersEdgesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
    //add surface elements
    Range CornersEdgesTets;
    rval = moab.get_adjacencies(CornersEdges,3,false,CornersEdgesTets,Interface::UNION); CHKERR_PETSC(rval);
    CornersEdgesTets = intersect(CornersEdgesTets,SurfacesTets);
    EntityHandle CornersEdgesTetsMeshset;
    rval = moab.create_meshset(MESHSET_SET,CornersEdgesTetsMeshset); CHKERR_PETSC(rval);	
    rval = moab.add_entities(CornersEdgesTetsMeshset,CornersEdgesTets); CHKERR_PETSC(rval);
    ierr = mField.add_ents_to_finite_element_by_TETs(CornersEdgesTetsMeshset,"C_EDGE_ELEM"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TETs(CornersEdgesTetsMeshset,"CTC_EDGE_ELEM"); CHKERRQ(ierr);
    rval = moab.delete_entities(&CornersEdgesTetsMeshset,1); CHKERR_PETSC(rval);
  }
  //add tets on surface
  EntityHandle SurfacesFacesMeshset;
  {
    Range SurfacesNodes;
    rval = moab.get_adjacencies(SurfacesFaces,0,false,SurfacesNodes,Interface::UNION); CHKERR_PETSC(rval);
    Range CornersEdgesNodes;
    rval = moab.get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
    //SurfacesNodes = subtract(SurfacesNodes,CornersNodes);
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

  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.add_problem("MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.add_problem("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.add_problem("C_ALL_MATRIX"); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("MATERIAL_MECHANICS","MATERIAL"); CHKERRQ(ierr);
  //
  ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_CORNER_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_EDGE_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("C_ALL_MATRIX","C_CORNER_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("C_ALL_MATRIX","C_EDGE_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("C_ALL_MATRIX","C_SURFACE_ELEM"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("MATERIAL_MECHANICS",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("CCT_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("C_ALL_MATRIX",bit_level0); CHKERRQ(ierr);

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_VERTICEs(SurfacesFacesMeshset,"LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_VERTICEs(CornersEdgesMeshset,"LAMBDA_EDGE"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_VERTICEs(CornersNodesMeshset,"LAMBDA_CORNER"); CHKERRQ(ierr);

  //set app. order
  int order = 2;
  ierr = mField.set_field_order(0,MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);
  
  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_SURFACE",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_EDGE",1); CHKERRQ(ierr);
  //ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CORNER",1); CHKERRQ(ierr);

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //partition
  ierr = mField.partition_problem("MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("MATERIAL_MECHANICS"); CHKERRQ(ierr);
  //partition
  ierr = mField.partition_problem("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  //partition
  ierr = mField.compose_problem("C_ALL_MATRIX","CCT_ALL_MATRIX","MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("C_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("C_ALL_MATRIX"); CHKERRQ(ierr);

  //create matrices
  Vec F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  SetPositionsEntMethod set_positions(moab);
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Col,set_positions); CHKERRQ(ierr);
  {
    EntityHandle node = no_handle;
    double coords[3];
    const DofMoFEMEntity_multiIndex *dofs_moabfield_ptr;
    ierr = mField.get_dofs_moabfield(&dofs_moabfield_ptr); CHKERRQ(ierr);
    DofMoFEMEntity_multiIndex::iterator dit = dofs_moabfield_ptr->begin();
    for(;dit!=dofs_moabfield_ptr->end();dit++) {
      if(dit->get_ent_type()!=MBVERTEX) continue;
      if(dit->get_name()!="MESH_NODE_POSITIONS") continue;
      EntityHandle ent = dit->get_ent();
      int dof_rank = dit->get_dof_rank();
      if(node!=ent) {
	rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
      }
      dit->get_FieldData() = coords[dof_rank];
    }
  }

  Range SideSet1,SideSet2;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 2 : %u\n",SideSet2.size());

  const double YoungModulus = 1;
  const double PoissonRatio = 0.25;
  ExampleDiriheltBC myDirihletBC(moab,SideSet1);
  Spatial_ElasticFEMethod MyFE(moab,&myDirihletBC,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),SideSet2);

  moabSnesCtx SnesCtx(mField,"ELASTIC_MECHANICS");
  
  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,Aij,Aij,SnesMat,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  ierr = MatSetOption(Aij,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(Aij,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

  moabSnesCtx::loops_to_do_type& loops_to_do_Rhs = SnesCtx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MyFE));
  moabSnesCtx::loops_to_do_type& loops_to_do_Mat = SnesCtx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MyFE));

  Vec D;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&D); CHKERRQ(ierr);
  ierr = mField.set_local_VecCreateGhost("ELASTIC_MECHANICS",Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  double step_size = -1e-3;
  for(int step = 1;step<2; step++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Load Setp %D\n",step); CHKERRQ(ierr);
    ierr = MyFE.set_t_val(step_size*step); CHKERRQ(ierr);
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

  PostProcVertexMethod ent_method(moab,"SPATIAL_POSITION");
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Col,ent_method); CHKERRQ(ierr);

  Vec F_MATERIAL;
  ierr = mField.VecCreateGhost("MATERIAL_MECHANICS",Row,&F_MATERIAL); CHKERRQ(ierr);
  MaterialForcesFEMethod MyMaterialFE(moab,&myDirihletBC,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),F_MATERIAL);
  ierr = VecZeroEntries(F_MATERIAL); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_MATERIAL,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_MATERIAL,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("MATERIAL_MECHANICS","MATERIAL",MyMaterialFE);  CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_MATERIAL,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_MATERIAL,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F_MATERIAL); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F_MATERIAL); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_MATERIAL,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_MATERIAL,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //ierr = VecView(F_MATERIAL,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost(
    "MATERIAL_MECHANICS","MESH_NODE_POSITIONS","MATERIAL_FORCE",Row,F_MATERIAL,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_material_forces(moab,"MESH_NODE_POSITIONS",F_MATERIAL,"MATERIAL_FORCE");
  ierr = mField.loop_dofs("MATERIAL_MECHANICS","MESH_NODE_POSITIONS",Row,ent_method_material_forces); CHKERRQ(ierr);

  matPROJ_ctx proj_all_ctx(mField,"MATERIAL_MECHANICS","C_ALL_MATRIX");
  ierr = mField.MatCreateMPIAIJWithArrays("MATERIAL_MECHANICS",&proj_all_ctx.K); CHKERRQ(ierr);
  ierr = mField.MatCreateMPIAIJWithArrays("C_ALL_MATRIX",&proj_all_ctx.C); CHKERRQ(ierr);
  C_SURFACE_FEMethod CFE_SURFACE_ALL(moab,SurfacesFacesMeshset,proj_all_ctx.C);
  C_EDGE_FEMethod CFE_EDGE_ALL(moab,CornersEdgesMeshset,proj_all_ctx.C);
  C_CORNER_FEMethod CFE_CORNER_ALL(moab,CornersNodes,proj_all_ctx.C);
  ierr = MatZeroEntries(proj_all_ctx.C); CHKERRQ(ierr);
  ierr = MatSetOption(proj_all_ctx.C,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(proj_all_ctx.C,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_SURFACE_ELEM",CFE_SURFACE_ALL);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_EDGE_ELEM",CFE_EDGE_ALL);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(proj_all_ctx.C,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(proj_all_ctx.C,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  //ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_CORNER_ELEM",CFE_CORNER_ALL);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(proj_all_ctx.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(proj_all_ctx.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  {
    //MatView(proj_all_ctx.C,PETSC_VIEWER_DRAW_WORLD);
    int m,n;
    MatGetSize(proj_all_ctx.C,&m,&n);
    PetscPrintf(PETSC_COMM_WORLD,"C size (%d,%d)\n",m,n);
    //std::string wait;
    //std::cin >> wait;
  }

  int M,N,m,n;
  ierr = MatGetSize(proj_all_ctx.K,&M,&N); CHKERRQ(ierr);
  ierr = MatGetLocalSize(proj_all_ctx.K,&m,&n); CHKERRQ(ierr);
  //
  Mat Q_ALL;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,&proj_all_ctx,&Q_ALL); CHKERRQ(ierr);
  ierr = MatShellSetOperation(Q_ALL,MATOP_MULT,(void(*)(void))matQ_mult_shell); CHKERRQ(ierr);

  Vec QTF_ALL_MATERIAL;
  ierr = VecDuplicate(F_MATERIAL,&QTF_ALL_MATERIAL); CHKERRQ(ierr);
  ierr = MatMult(Q_ALL,F_MATERIAL,QTF_ALL_MATERIAL); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(QTF_ALL_MATERIAL,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(QTF_ALL_MATERIAL,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = proj_all_ctx.DestroyQorP(); CHKERRQ(ierr);
  
  ierr = mField.set_other_global_VecCreateGhost(
    "MATERIAL_MECHANICS","MESH_NODE_POSITIONS","MATERIAL_FORCE",Row,QTF_ALL_MATERIAL,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_qt_all_material_forces(moab,"MESH_NODE_POSITIONS",QTF_ALL_MATERIAL,"QT_ALL_MATERIAL_FORCE");
  ierr = mField.loop_dofs("MATERIAL_MECHANICS","MESH_NODE_POSITIONS",Row,ent_method_qt_all_material_forces); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PostProcFieldsAndGradientOnRefMesh fe_post_proc_method(moab);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  if(pcomm->rank()==0) {
    ostringstream sss;
    sss << "save_space_solution.h5m";
    rval = moab.write_file(sss.str().c_str()); CHKERR_PETSC(rval);
  }

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&F_MATERIAL); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = MatDestroy(&proj_all_ctx.K); CHKERRQ(ierr);
  ierr = MatDestroy(&proj_all_ctx.C); CHKERRQ(ierr);
  ierr = MatDestroy(&Q_ALL); CHKERRQ(ierr);
  ierr = VecDestroy(&QTF_ALL_MATERIAL); CHKERRQ(ierr);

  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

  return 0;
}



