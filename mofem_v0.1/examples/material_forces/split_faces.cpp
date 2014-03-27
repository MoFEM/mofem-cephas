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
    ierr = main_refine_and_meshcat(mField,face_splitting,false,2); CHKERRQ(ierr);

  }

  { //find faces for split

    FaceSplittingTools face_splitting(mField);
    ierr = main_select_faces_for_splitting(mField,face_splitting,2); CHKERRQ(ierr);
    //do splittig
    ierr = main_split_faces_and_update_field_and_elements(mField,face_splitting,2); CHKERRQ(ierr);
    //project and set coords
    ierr = conf_prob.set_spatial_positions(mField); CHKERRQ(ierr);
    ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);
    ierr = face_splitting.projectCrackFrontNodes(); CHKERRQ(ierr);

  }
 

  Tag th_griffith_force;
  rval = moab.tag_get_handle("GRIFFITH_FORCE",th_griffith_force); CHKERR_PETSC(rval);
  rval = moab.tag_delete(th_griffith_force); CHKERR_PETSC(rval);
  Tag th_freez;
  rval = moab.tag_get_handle("FROZEN_NODE",th_freez); CHKERR_PETSC(rval);
  rval = moab.tag_delete(th_freez); CHKERR_PETSC(rval);

  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_material_positions) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_spatial_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_constrains_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_material_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_constrains_crack_front_problem_definition) = 0;


  EntityHandle meshset_level0;
  rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  ierr = conf_prob.spatial_problem_definition(mField); CHKERRQ(ierr);
  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_set_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);
  ierr = main_arc_length_setup(mField,conf_prob); CHKERRQ(ierr);

  if(pcomm->rank()==0) {

    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("cat_out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);

  }

  {

    double gc;
    PetscBool flg;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
    }

    Tag th_t_val;
    rval = mField.get_moab().tag_get_handle("_LoadFactor_t_val",th_t_val); CHKERR_PETSC(rval);
    double *load_factor_ptr;
    rval = mField.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&load_factor_ptr); CHKERR_THROW(rval);
    double& load_factor = *load_factor_ptr;

    {
      //solve spatial problem
      SNES snes;
      ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);  
      ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
      ierr = conf_prob.solve_spatial_problem(mField,&snes,false); CHKERRQ(ierr);
      ierr = conf_prob.calculate_material_forces(mField,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
      ierr = conf_prob.front_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.surface_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.project_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.griffith_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.griffith_g(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
    }

  }

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
