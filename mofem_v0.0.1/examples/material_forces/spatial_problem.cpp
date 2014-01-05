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

  PetscInt nb_ref_levels;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_ref",&nb_ref_levels,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    nb_ref_levels = 0;
  }
 
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  FieldCore core(moab);
  FieldInterface& mField = core;

  ierr = mField.printCubitDisplacementSet(); CHKERRQ(ierr);
  ierr = mField.printCubitPressureSet(); CHKERRQ(ierr);
  ierr = mField.printCubitMaterials(); CHKERRQ(ierr);

  Tag th_my_ref_level;
  BitRefLevel def_bit_level = 0;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",sizeof(BitRefLevel),MB_TYPE_OPAQUE,
    th_my_ref_level,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level); 
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  ConfigurationalFractureMechanics conf_prob(mField);

  ierr = conf_prob.set_material_fire_wall(mField); CHKERRQ(ierr);

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  //BitRefLevel bit_level0;
  //bit_level0.set(0);

  PetscBool no_add_interface;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-my_no_add_interface",&no_add_interface,&flg); CHKERRQ(ierr);
  if(no_add_interface == PETSC_TRUE) {
    conf_prob.material_FirelWall->set(ConfigurationalFractureMechanics::FW_add_crack);
  }

  BitRefLevel bit_level_interface;
  if(mField.check_msId_meshset(200,SideSet)) {
  if(!conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_add_crack)) {

    conf_prob.material_FirelWall->set(ConfigurationalFractureMechanics::FW_add_crack);

    //Interface
    EntityHandle meshset_interface;
    ierr = mField.get_msId_meshset(200,SideSet,meshset_interface); CHKERRQ(ierr);
    ierr = mField.get_msId_3dENTS_sides(meshset_interface,true); CHKERRQ(ierr);
    // stl::bitset see for more details
    bit_level_interface.set(0);
    ierr = mField.get_msId_3dENTS_split_sides(0,bit_level_interface,meshset_interface,false,true); CHKERRQ(ierr);

    //add refined ent to cubit meshsets
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
      EntityHandle cubit_meshset = cubit_it->meshset; 
      ierr = mField.refine_get_childern(cubit_meshset,bit_level_interface,cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
      ierr = mField.refine_get_childern(cubit_meshset,bit_level_interface,cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
      ierr = mField.refine_get_childern(cubit_meshset,bit_level_interface,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
      ierr = mField.refine_get_childern(cubit_meshset,bit_level_interface,cubit_meshset,MBTET,true); CHKERRQ(ierr);
    }

  }} else {
    bit_level_interface.set(0);
    ierr = mField.seed_ref_level_3D(0,bit_level_interface); CHKERRQ(ierr);
  }

  if(mField.check_msId_meshset(201,SideSet)) {
  if(!conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_refine_near_crack_tip)) {

    conf_prob.material_FirelWall->set(ConfigurationalFractureMechanics::FW_refine_near_crack_tip);

    BitRefLevel last_ref = bit_level_interface;
    for(int ll = 1;ll<nb_ref_levels+1;ll++) {

      Range crack_edges,crack_nodes,edge_tets,level_tets,edges_to_refine;

      ierr = mField.refine_get_ents(last_ref,BitRefLevel().set(),level_tets); CHKERRQ(ierr);

      ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_edges,true); CHKERRQ(ierr);
      rval = moab.get_connectivity(crack_edges,crack_nodes,true); CHKERR_PETSC(rval);
      rval = moab.get_adjacencies(crack_nodes,3,false,edge_tets,Interface::UNION); CHKERR_PETSC(rval);
      edge_tets = intersect(level_tets.subset_by_type(MBTET),edge_tets);
      rval = moab.get_adjacencies(edge_tets,1,false,edges_to_refine,Interface::UNION); CHKERR_PETSC(rval);

      last_ref = BitRefLevel().set(ll);
      ierr = mField.add_verices_in_the_middel_of_edges(edges_to_refine,last_ref,2); CHKERRQ(ierr);
      ierr = mField.refine_TET(level_tets,last_ref,false); CHKERRQ(ierr);

      for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
	EntityHandle cubit_meshset = cubit_it->meshset; 
	ierr = mField.refine_get_childern(cubit_meshset,last_ref,cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
	ierr = mField.refine_get_childern(cubit_meshset,last_ref,cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
	ierr = mField.refine_get_childern(cubit_meshset,last_ref,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
	ierr = mField.refine_get_childern(cubit_meshset,last_ref,cubit_meshset,MBTET,true); CHKERRQ(ierr);
      }
  
    }

    bit_level0 = last_ref;


  }} else {
    bit_level0 = bit_level_interface;
  } 

  double *t_val;
  Tag th_t_val;
  double def_t_val = 0;
  rval = mField.get_moab().tag_get_handle("_LoadFactor_t_val",1,MB_TYPE_DOUBLE,th_t_val,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&def_t_val); 
  if(rval == MB_ALREADY_ALLOCATED) {
    rval = mField.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&t_val); CHKERR_PETSC(rval);
  } else {
    CHKERR_PETSC(rval);
    rval = mField.get_moab().tag_set_data(th_t_val,&root_meshset,1,&def_t_val); CHKERR_PETSC(rval);
    rval = mField.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&t_val); CHKERR_PETSC(rval);
  }

  if(!conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_load_factor)) {

    conf_prob.material_FirelWall->set(ConfigurationalFractureMechanics::FW_set_load_factor);

    ierr = PetscOptionsGetReal(PETSC_NULL,"-my_load",t_val,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_WORLD,1,"*** ERROR -my_load (what is load factor?)");
    }

  }

  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.refine_get_ents(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  ierr = conf_prob.spatial_problem_definition(mField); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_set_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

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

  //solve problem
  ierr = conf_prob.set_spatial_positions(mField); CHKERRQ(ierr);
  ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);

  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);  
  ierr = conf_prob.solve_spatial_problem(mField,&snes); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    rval = moab.write_file("out_spatial.h5m"); CHKERR_PETSC(rval);
  }

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    if(conf_prob.fe_post_proc_stresses_method!=NULL) {
      rval = conf_prob.fe_post_proc_stresses_method->moab_post_proc.write_file("out_stresses.vtk","VTK",""); CHKERR_PETSC(rval);
    }
  }

  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  return 0;
}




