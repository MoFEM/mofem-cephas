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

#ifdef WITH_TETGEN

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
#include <ConstrainMatrixCtx.hpp>

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

#include <adolc/adolc.h> 
#include <ConvectiveMassElement.hpp>
#include <ConfigurationalFractureForDynamics.hpp>

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

  ConfigurationalFracturDynamics conf_prob(m_field);
  
  { //definitions
    ierr = conf_prob.set_material_fire_wall(m_field); CHKERRQ(ierr);
    ierr = conf_prob.constrains_problem_definition(m_field); CHKERRQ(ierr);
    ierr = conf_prob.coupled_dynamic_problem_definition(m_field); CHKERRQ(ierr);
    ierr = conf_prob.constrains_crack_front_problem_definition(m_field,"COUPLED_DYNAMIC"); CHKERRQ(ierr);
  }
  
  { //build
  
    //ref meshset ref level 0
    Tag th_my_ref_level;
    rval = m_field.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
    const EntityHandle root_meshset = m_field.get_moab().get_root_set();
    BitRefLevel *ptr_bit_level0;
    rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
    BitRefLevel& bit_level0 = *ptr_bit_level0;

    ierr = m_field.modify_problem_ref_level_set_bit("C_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
    ierr = m_field.modify_problem_ref_level_set_bit("CCT_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
    ierr = m_field.modify_problem_ref_level_set_bit("C_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);
    ierr = m_field.modify_problem_ref_level_set_bit("CTC_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);
    ierr = m_field.modify_problem_ref_level_set_bit("COUPLED_DYNAMIC",bit_level0); CHKERRQ(ierr);

    //add finite elements entities
    ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC_COUPLED",MBTET); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MATERIAL_COUPLED",MBTET); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MESH_SMOOTHER",MBTET); CHKERRQ(ierr);

    //build field
    ierr = m_field.build_fields(); CHKERRQ(ierr);
    //build finite elemnts
    ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
    //build adjacencies
    ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
    //build problem
    ierr = m_field.build_problems(); CHKERRQ(ierr);

  }

  { //partition
    ierr = conf_prob.coupled_dynamic_partition_problems(m_field); CHKERRQ(ierr);
    ierr = conf_prob.constrains_partition_problems(m_field,"COUPLED_DYNAMIC"); CHKERRQ(ierr);
    ierr = conf_prob.crackfront_partition_problems(m_field,"COUPLED_DYNAMIC"); CHKERRQ(ierr);
  }

  { 
    int zero = 0;
    Tag th_step;
    rval = m_field.get_moab().tag_get_handle("_TsStep_",1,MB_TYPE_INTEGER,th_step,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&zero); 
    if(rval != MB_ALREADY_ALLOCATED) {
      //init data
      //conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
      //conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_material_positions) = 0;
      //ierr = conf_prob.set_spatial_positions(m_field); CHKERRQ(ierr);
      //ierr = conf_prob.set_material_positions(m_field);  CHKERRQ(ierr);
    }
  }

  //create tS
  TS ts;
  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);
  ierr = conf_prob.solve_dynmaic_problem(m_field,ts); CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);

  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  return 0;
}
