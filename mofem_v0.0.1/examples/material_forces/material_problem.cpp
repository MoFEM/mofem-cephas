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
#include "FieldCore.hpp"

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
 
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  FieldCore core(moab);
  FieldInterface& mField = core;

  ConfigurationalMechanics conf_prob(mField);

  ierr = conf_prob.set_material_fire_wall(mField); CHKERRQ(ierr);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  ierr = conf_prob.spatial_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.material_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.constrains_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.constrains_crack_front_problem_definition(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MATERIAL",MBTET); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("MATERIAL_MECHANICS",bit_level0); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("CCT_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("C_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("C_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("CTC_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);

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

  //solve material problem
  ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);
  ierr = conf_prob.surface_projection_data(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.front_projection_data(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);

  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = conf_prob.solve_material_problem(mField,&snes); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  ierr = conf_prob.calculate_material_forces(mField,"MATERIAL_MECHANICS","MATERIAL"); CHKERRQ(ierr);
  ierr = conf_prob.project_force_vector(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.griffith_force_vector(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.griffith_g(mField,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.calculate_spatial_residual(mField); CHKERRQ(ierr);

  rval = moab.write_file("out_material_problem.h5m"); CHKERR_PETSC(rval);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("MATERIAL_MECHANICS","MATERIAL",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  return 0;
}
