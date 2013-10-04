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

  ierr = mField.printCubitDisplacementSet(); CHKERRQ(ierr);
  ierr = mField.printCubitPressureSet(); CHKERRQ(ierr);

  ConfigurationalMechanics conf_prob;

  ierr = conf_prob.ConfigurationalMechanics_SetMaterialFireWall(mField); CHKERRQ(ierr);

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,1); CHKERRQ(ierr);
  //BitRefLevel bit_level0;
  //bit_level0.set(0);

  //Interface
  EntityHandle meshset_interface;
  ierr = mField.get_msId_meshset(200,SideSet,meshset_interface); CHKERRQ(ierr);
  ierr = mField.get_msId_3dENTS_sides(meshset_interface,true); CHKERRQ(ierr);
  // stl::bitset see for more details
  BitRefLevel bit_level_interface;
  bit_level_interface.set(1);
  ierr = mField.get_msId_3dENTS_split_sides(0,bit_level_interface,meshset_interface,false,true); CHKERRQ(ierr);
  EntityHandle meshset_level_interface;
  rval = moab.create_meshset(MESHSET_SET,meshset_level_interface); CHKERR_PETSC(rval);
  ierr = mField.refine_get_ents(bit_level_interface,meshset_level_interface); CHKERRQ(ierr);

  //add refined ent to cubit meshsets
  for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
    EntityHandle cubit_meshset = cubit_it->meshset; 
    ierr = mField.refine_get_childern(cubit_meshset,bit_level_interface,cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(cubit_meshset,bit_level_interface,cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(cubit_meshset,bit_level_interface,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
  }

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(2);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(meshset_level_interface,bit_level0); CHKERRQ(ierr);
  ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

  ierr = conf_prob.ConfigurationalMechanics_SpatialProblemDefinition(mField); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = conf_prob.ConfigurationalMechanics_SpatialPartitionProblems(mField); CHKERRQ(ierr);

  //solve problem
  ierr = conf_prob.ConfigurationalMechanics_SetSpatialPositions(mField); CHKERRQ(ierr);
  ierr = conf_prob.ConfigurationalMechanics_SolveSpatialProblem(mField); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    rval = moab.write_file("out_spatial.h5m"); CHKERR_PETSC(rval);
  }

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

  return 0;
}




