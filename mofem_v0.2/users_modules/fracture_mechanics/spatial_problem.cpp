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

#ifdef WITH_TETGEM

#include <tetgen.h>
#ifdef REAL
  #undef REAL
#endif

#endif
#include <MoFEM.hpp>
using namespace MoFEM;

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>
#include <ArcLengthTools.hpp>
#include <MatShellConstrainsByMarkAinsworth.hpp>

extern "C" {
  #include <complex_for_lazy.h>
}

#include <FEMethod_ComplexForLazy.hpp>
#include <FEMethod_DriverComplexForLazy.hpp>

#include <moab/Skinner.hpp>
#include <moab/AdaptiveKDTree.hpp>

#include <petsctime.h>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>
#include <PostProcNonLinearElasticityStresseOnRefindedMesh.hpp>

#include <FaceSplittingTool.hpp>
#include <ConfigurationalFractureMechanics.hpp>

using namespace ObosleteUsersModules;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  try {

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  const EntityHandle root_meshset = moab.get_root_set();

  PetscInt nb_ref_levels;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_ref",&nb_ref_levels,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    nb_ref_levels = 0;
  }

  int def_set_ref_level = 0;
  Tag th_set_ref_level;
  rval = moab.tag_get_handle("_SET_REF_LEVEL",1,MB_TYPE_INTEGER,
    th_set_ref_level,MB_TAG_CREAT|MB_TAG_MESH,&def_set_ref_level); 
  rval = moab.tag_set_data(th_set_ref_level,&root_meshset,1,&nb_ref_levels); CHKERR_PETSC(rval);

  PetscInt order;
  flg = PETSC_TRUE;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }

  int def_set_order = 1;
  Tag th_set_order;
  rval = moab.tag_get_handle("_SET_ORDER",1,MB_TYPE_INTEGER,
    th_set_order,MB_TAG_CREAT|MB_TAG_MESH,&def_set_order); 
  rval = moab.tag_set_data(th_set_order,&root_meshset,1,&order); CHKERR_PETSC(rval);
 
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

  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  ierr = m_field.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field.print_cubit_pressure_set(); CHKERRQ(ierr);
  ierr = m_field.print_cubit_materials_set(); CHKERRQ(ierr);

  Tag th_my_ref_level;
  BitRefLevel def_bit_level = 0;
  rval = m_field.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",sizeof(BitRefLevel),MB_TYPE_OPAQUE,
    th_my_ref_level,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level); 
  BitRefLevel *ptr_bit_level0;
  rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  ConfigurationalFractureMechanics conf_prob(m_field);
  ierr = conf_prob.set_material_fire_wall(m_field); CHKERRQ(ierr);

  //mesh refine and split faces
  FaceSplittingTools face_splitting_tools(m_field);

  PetscBool no_add_interface = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-my_restart",&no_add_interface,&flg); CHKERRQ(ierr);
  if(no_add_interface == PETSC_TRUE) {
    conf_prob.material_FirelWall->set(ConfigurationalFractureMechanics::FW_add_crack);
  } else {
    bit_level0 = BitRefLevel().set(0);
    ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
    face_splitting_tools.meshRefineBitLevels.resize(0);
    face_splitting_tools.meshRefineBitLevels.push_back(0);
  }

  if(!conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_refine_near_crack_tip)) {
    conf_prob.material_FirelWall->set(ConfigurationalFractureMechanics::FW_refine_near_crack_tip);
    ierr = face_splitting_tools.meshRefine(); CHKERRQ(ierr);
  }

  //if(!conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_add_crack)) {
    //conf_prob.material_FirelWall->set(ConfigurationalFractureMechanics::FW_add_crack);
    //ierr = face_splitting_tools.splitFaces(); CHKERRQ(ierr);
  //}

  #ifdef WITH_TETGEM
    char switches1[] = "pAz";
    ierr = face_splitting_tools.rebuildMeshWithTetGen(switches1,1); CHKERRQ(ierr);	
    //char switches2[] = "rO/1AzFV";
    //ierr = face_splitting_tools.rebuildMeshWithTetGen(switches2,1); CHKERRQ(ierr);	
  #endif
  bit_level0 = BitRefLevel().set(face_splitting_tools.meshIntefaceBitLevels.back());

  //PetscFinalize();
  //return 0;


  //load factor
  double *t_val;
  Tag th_t_val;
  double def_t_val = 0;
  rval = m_field.get_moab().tag_get_handle("_LoadFactor_Scale_",1,MB_TYPE_DOUBLE,th_t_val,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&def_t_val); 
  if(rval == MB_ALREADY_ALLOCATED) {
    rval = m_field.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&t_val); CHKERR_PETSC(rval);
  } else {
    CHKERR_PETSC(rval);
    rval = m_field.get_moab().tag_set_data(th_t_val,&root_meshset,1,&def_t_val); CHKERR_PETSC(rval);
    rval = m_field.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&t_val); CHKERR_PETSC(rval);
  }

  if(!conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_load_factor)) {

    conf_prob.material_FirelWall->set(ConfigurationalFractureMechanics::FW_set_load_factor);

    ierr = PetscOptionsGetReal(PETSC_NULL,"-my_load",t_val,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_WORLD,1,"*** ERROR -my_load (what is the load factor?)");
    }

  }

  ierr = main_spatial_solution(m_field,conf_prob); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Save results in out_spatial.h5m\n"); CHKERRQ(ierr);
    rval = moab.write_file("out_spatial.h5m"); CHKERR_PETSC(rval);
  }

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = m_field.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Save results in VTK mesh out.vtk\n"); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    if(conf_prob.fe_post_proc_stresses_method!=NULL) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Save stresses in post-processing mesh out_stresses.vtk\n"); CHKERRQ(ierr);
      rval = conf_prob.fe_post_proc_stresses_method->moab_post_proc.write_file("out_stresses.vtk","VTK",""); CHKERR_PETSC(rval);
    }
  }

  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);


  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  return 0;
}




