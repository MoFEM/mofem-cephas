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
#include "FaceSplittinfTool.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

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
 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  BARRIER_RANK_START(pcomm) 
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  BARRIER_RANK_END(pcomm) 

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  FieldCore core(moab);
  FieldInterface& mField = core;

  ConfigurationalFractureMechanics conf_prob(mField);

  ierr = conf_prob.set_material_fire_wall(mField); CHKERRQ(ierr);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  ierr = mField.build_fields(); CHKERRQ(ierr);
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  { //cat mesh
  
    FaceSplittingTools face_splitting(mField);
    //ierr = main_refine_and_meshcat(mField,face_splitting,false,2); CHKERRQ(ierr);

  }

  //find faces for split

  FaceSplittingTools face_splitting(mField);
  ierr = main_select_faces_for_splitting(mField,face_splitting,2); CHKERRQ(ierr);
  //do splittig
  ierr = main_split_faces_and_update_field_and_elements(mField,face_splitting,2); CHKERRQ(ierr);
	
  //rebuild fields, finite elementa and problems
  ierr = main_face_splitting_restart(mField,conf_prob); CHKERRQ(ierr);

  //solve spatial problem and calulate griffith forces

  //project and set coords
  ierr = conf_prob.set_spatial_positions(mField); CHKERRQ(ierr);
  ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);
  
  SNES snes;

  //solve mesh smoothing
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);  
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
  SNESLineSearch linesearch;
  ierr = SNESGetLineSearch(snes,&linesearch); CHKERRQ(ierr);
  ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHL2); CHKERRQ(ierr);
  Vec D_tmp_mesh_positions;
  ierr = mField.VecCreateGhost("MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS",Col,&D_tmp_mesh_positions); CHKERRQ(ierr);
  int nb_sub_steps = 1;
  int nn;
  do { 
    nn = 1;
    for(;nn<=nb_sub_steps;nn++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Mesh projection substep = %D\n",nn); CHKERRQ(ierr);
      ierr = mField.set_local_VecCreateGhost(
	"MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS",
	Col,D_tmp_mesh_positions,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      double alpha = (double)nn/(double)nb_sub_steps;
      ierr = face_splitting.calculateDistanceCrackFrontNodesFromCrackSurface(alpha); CHKERRQ(ierr);
      //project nodes on crack surface
      ierr = face_splitting.projectCrackFrontNodes(); CHKERRQ(ierr);
      ierr = conf_prob.solve_mesh_smooting_problem(mField,&snes); 
      SNESConvergedReason reason;
      SNESGetConvergedReason(snes,&reason); CHKERRQ(ierr);
      if(ierr == 0 && reason > 0 && reason != 5) {
	ierr = mField.set_local_VecCreateGhost(
	  "MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS",
	  Col,D_tmp_mesh_positions,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      } else {
	ierr = mField.set_local_VecCreateGhost(
	  "MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS",
	  Col,D_tmp_mesh_positions,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	nb_sub_steps++;
	break;
      }
    }
    if(nb_sub_steps == 10) break;
  } while(nn-1 != nb_sub_steps);
  ierr = VecDestroy(&D_tmp_mesh_positions); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  ierr = conf_prob.set_coordinates_from_material_solution(mField); CHKERRQ(ierr);
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_material_positions) = 0;
  ierr = conf_prob.set_spatial_positions(mField); CHKERRQ(ierr);
  ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);

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

  if(pcomm->rank()==0) {

    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("cat_out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);

  }

  if(pcomm->rank()==0) {

    EntityHandle meshset200;
    ierr = mField.get_Cubit_msId_meshset(200,SideSet,meshset200); CHKERRQ(ierr);
    rval = moab.write_file("cat_CrackFrontSurface.vtk","VTK","",&meshset200,1); CHKERR_PETSC(rval);
    EntityHandle meshset201;
    ierr = mField.get_Cubit_msId_meshset(201,SideSet,meshset201); CHKERRQ(ierr);
    rval = moab.write_file("cat_CrackFrontEdges.vtk","VTK","",&meshset201,1); CHKERR_PETSC(rval);

  }

  if(pcomm->rank()==0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Save results in out_split_faces.h5m\n"); CHKERRQ(ierr);
    rval = moab.write_file("out_split_faces.h5m"); CHKERR_PETSC(rval);
  }

  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  return 0;
}
