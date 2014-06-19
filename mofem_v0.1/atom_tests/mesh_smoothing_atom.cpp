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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"

#include "SurfacePressureComplexForLazy.hpp"
#include "FEMethod_DriverComplexForLazy.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

//Rounding
#define RND_EPS 1e-6

struct MyMeshSmoothing_ElasticFEMethod_LagnageMultiplaiers: public MeshSmoothingFEMethod  {

  MyMeshSmoothing_ElasticFEMethod_LagnageMultiplaiers(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    MeshSmoothingFEMethod(_mField,_dirihlet_bc_method_ptr,_verbose) {
    set_qual_ver(1);
  }

};

int main(int argc, char *argv[]) {

  try {

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

  FieldCore core(moab);
  FieldInterface& mField = core;

  Range CubitSideSets_meshsets;
  ierr = mField.get_Cubit_meshsets(SideSet,CubitSideSets_meshsets); CHKERRQ(ierr);

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  BitRefLevel problem_level = bit_level0;

  //Fields
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("LAMBDA_SURFACE",H1,1); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("MESH_SMOOTHER"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("C_SURFACE_ELEM"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("C_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("C_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("C_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("C_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);

  ierr = mField.add_problem("MESH_SMOOTHING"); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("MESH_SMOOTHING","MESH_SMOOTHER"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("MESH_SMOOTHING","C_SURFACE_ELEM"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("MESH_SMOOTHING",problem_level); CHKERRQ(ierr);

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_level,"MESH_SMOOTHER",MBTET); CHKERRQ(ierr);

  //add tets on corners
  Range corner_nodes;
  EntityHandle coner_nodes_meshset,surface_faces_meshset;
  {
    Range corner_edges,surface_faces;
    ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,corner_edges,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,corner_nodes,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(102,SideSet,2,surface_faces,true); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 100 = %d\n",corner_edges.size()); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of NodeSet 101 = %d\n",corner_nodes.size()); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 102 = %d\n",surface_faces.size()); CHKERRQ(ierr);
    ierr = mField.seed_finite_elements(surface_faces); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(surface_faces,"C_SURFACE_ELEM"); CHKERRQ(ierr);

    if(surface_faces.empty()) SETERRQ(PETSC_COMM_SELF,1,"no surface elements");
    Range corner_edges_nodes;
    rval = moab.get_connectivity(corner_edges,corner_edges_nodes,true); CHKERR_PETSC(rval);
    corner_nodes.insert(corner_edges_nodes.begin(),corner_edges_nodes.end());
    {
      rval = moab.create_meshset(MESHSET_SET,coner_nodes_meshset); CHKERR_PETSC(rval);	
      rval = moab.add_entities(coner_nodes_meshset,corner_nodes); CHKERR_PETSC(rval);
      //add surface elements
      ierr = mField.seed_finite_elements(corner_nodes); CHKERRQ(ierr);
    }
    {
      Range surface_nodes;
      rval = moab.get_connectivity(surface_faces,surface_nodes,true); CHKERR_PETSC(rval);
      surface_nodes = subtract(surface_nodes,corner_nodes);
      rval = moab.create_meshset(MESHSET_SET,surface_faces_meshset); CHKERR_PETSC(rval);
      rval = moab.add_entities(surface_faces_meshset,surface_nodes); CHKERR_PETSC(rval);
    }
  }

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_VERTICEs(surface_faces_meshset,"LAMBDA_SURFACE"); CHKERRQ(ierr);

  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(surface_faces_meshset,MBVERTEX,"LAMBDA_SURFACE",1); CHKERRQ(ierr);

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(problem_level); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition
  ierr = mField.partition_problem("MESH_SMOOTHING"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("MESH_SMOOTHING"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("MESH_SMOOTHING"); CHKERRQ(ierr);

  {
    EntityHandle node = 0;
    double coords[3];
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"MESH_NODE_POSITIONS",dof_ptr)) {
      if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
      EntityHandle ent = dof_ptr->get_ent();
      int dof_rank = dof_ptr->get_dof_rank();
      double &fval = dof_ptr->get_FieldData();
      if(node!=ent) {
	rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
      }
      fval = coords[dof_rank];
    }
  }

  Mat K;
  ierr = mField.MatCreateMPIAIJWithArrays("MESH_SMOOTHING",&K); CHKERRQ(ierr);
  Vec F;
  ierr = mField.VecCreateGhost("MESH_SMOOTHING",Row,&F); CHKERRQ(ierr);
  Vec D;
  ierr = mField.VecCreateGhost("MESH_SMOOTHING",Col,&D); CHKERRQ(ierr);

  //materialDirihletBC myDirihletBC(moab,corner_nodes);
  BaseDirihletBC myDirihletBC;

  FixMaterialPoints fix_material_pts(mField,"MESH_NODE_POSITIONS",K,D,F,corner_nodes);
  //fix_material_pts.field_names.push_back("LAMBDA_SURFACE");
  MyMeshSmoothing_ElasticFEMethod_LagnageMultiplaiers bulk_fe(mField,&myDirihletBC);
  C_SURFACE_FEMethod_ForSnes surface_fe(mField,&myDirihletBC);
  surface_fe.nonlinear = true;

  SnesCtx snes_ctx(mField,"MESH_SMOOTHING");

  snes_ctx.get_preProcess_to_do_Rhs().push_back(&fix_material_pts);
  SnesCtx::loops_to_do_type& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("C_SURFACE_ELEM",&surface_fe));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MESH_SMOOTHER",&bulk_fe));
  snes_ctx.get_postProcess_to_do_Rhs().push_back(&fix_material_pts);

  snes_ctx.get_preProcess_to_do_Mat().push_back(&fix_material_pts);
  SnesCtx::loops_to_do_type& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("C_SURFACE_ELEM",&surface_fe));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MESH_SMOOTHER",&bulk_fe));
  snes_ctx.get_postProcess_to_do_Mat().push_back(&fix_material_pts);

  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,K,K,SnesMat,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  /*{
    MatStructure flg = SAME_NONZERO_PATTERN;
    ierr = SnesMat(snes,D,&K,&K,&flg,&snes_ctx); CHKERRQ(ierr);
    MatView(K,PETSC_VIEWER_DRAW_WORLD);
    //MatView(K,PETSC_VIEWER_STDOUT_WORLD);
    std::string wait;
    std::cin >> wait;
  }*/

  ierr = mField.set_local_VecCreateGhost("MESH_SMOOTHING",Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  PetscViewer viewer;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mesh_smoothing.txt",&viewer);
  ierr = VecChop(D,RND_EPS); CHKERRQ(ierr);
  ierr = VecView(D,viewer); CHKERRQ(ierr);
  PetscViewerDestroy(&viewer);


  ierr = MatDestroy(&K); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);

  PetscFinalize();

  } catch (const char* msg) {
        SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  return 0;
}



