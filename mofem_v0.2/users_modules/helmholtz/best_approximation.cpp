/** \file best_approximation.cpp
  \ingroup mofem_helmholtz_elem

  Calculates best approximation for incident wave problem. 

 */

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.

 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. 
*/

#include <MoFEM.hpp>
#include <Projection10NodeCoordsOnField.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <FieldApproximation.hpp>
#include <PotsProcOnRefMesh.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <petsctime.h>
#include <fstream>
#include <iostream>

#include <stdexcept>
#include <cmath>
#include <boost/math/special_functions.hpp>
#include <complex>

using namespace std;
using namespace boost::math;
namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;
using namespace MoFEM;

#include <AnalyticalSolutions.hpp>

static char help[] = "...\n\n";

// argc = argument counts, argv = argument vectors
int main(int argc, char *argv[]) {
  
  ErrorCode rval;
  PetscErrorCode ierr;
  
  PetscInitialize(&argc,&argv,(char *)0,help);
  
  moab::Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  PetscBool is_partitioned = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-my_is_partitioned",&is_partitioned,&flg); CHKERRQ(ierr);
  if(is_partitioned == PETSC_TRUE) {
    //Read mesh to MOAB
    const char *option;
    option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  } else {
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  }

  // create MoFEM database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;
  
  // count the comsumption of time by single run
  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
  
  // set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  
  // define fields
  ierr = m_field.add_field("reEX",H1,1); CHKERRQ(ierr);
  ierr = m_field.add_field("imEX",H1,1); CHKERRQ(ierr);
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  
  // define finite element
  ierr = m_field.add_finite_element("FE1"); CHKERRQ(ierr);
  
  // Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("FE1","reEX"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("FE1","reEX"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("FE1","reEX"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("FE1","imEX"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("FE1","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  if(m_field.check_field("rePRES") && m_field.check_field("imPRES")) {

    ierr = m_field.modify_finite_element_add_field_data("FE1","rePRES"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("FE1","imPRES"); CHKERRQ(ierr);

  }
  
  // meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  // add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"reEX"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"imEX"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  // add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_TETs(root_set,"FE1"); CHKERRQ(ierr);
  
  // set app. order
  // see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 3;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 3;
  }
  ierr = m_field.set_field_order(root_set,MBTET,"reEX",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"reEX",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"reEX",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"reEX",1); CHKERRQ(ierr);
  
  ierr = m_field.set_field_order(root_set,MBTET,"imEX",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"imEX",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"imEX",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"imEX",1); CHKERRQ(ierr);

  ierr = m_field.set_field_order(root_set,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  // build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
  // build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  // build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  // define problem
  ierr = m_field.add_problem("EX1_PROBLEM"); CHKERRQ(ierr);
  // set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("EX1_PROBLEM","FE1"); CHKERRQ(ierr);
  // set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("EX1_PROBLEM",bit_level0); CHKERRQ(ierr);

  // build porblems
  if(is_partitioned) {
    // if mesh is partitioned
    ierr = m_field.build_partitioned_problem("EX1_PROBLEM",true); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("EX1_PROBLEM",true); CHKERRQ(ierr);
  } else {
    // if not partitioned mesh is load to all processes 
    ierr = m_field.build_problem("EX1_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_problem("EX1_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("EX1_PROBLEM"); CHKERRQ(ierr);
  }
  ierr = m_field.partition_ghost_dofs("EX1_PROBLEM"); CHKERRQ(ierr);
  
  PetscBool wavenumber_flg;
  double wavenumber;
  // set wave number from line command, that overwrite numbre form block set
  ierr = PetscOptionsGetScalar(NULL,"-wave_number",&wavenumber,&wavenumber_flg); CHKERRQ(ierr);
  if(!wavenumber_flg) {

    SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA,"wave number not given, set in line command -wave_number to fix problem");

  }

  double power_of_incident_wave = 1;
  ierr = PetscOptionsGetScalar(NULL,"-power_of_incident_wave",&power_of_incident_wave,NULL); CHKERRQ(ierr);

  //wave direction unit vector=[x,y,z]^T
  ublas::vector<double> wave_direction;
  wave_direction.resize(3);
  wave_direction.clear();
  wave_direction[2] = 1; // default:X direction [0,0,1]

  int nmax = 3;
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-wave_direction",&wave_direction[0],&nmax,NULL); CHKERRQ(ierr);
  if(nmax > 0 && nmax != 3) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA,"*** ERROR -wave_direction [3*1 vector] default:X direction [0,0,1]");
  }
  
  PetscInt choise_value = 0;
  // set type of analytical solution  
  ierr = PetscOptionsGetEList(NULL,"-analytical_solution_type",analytical_solution_types,6,&choise_value,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -analytical_solution_type needed, WARNING!!!!!!.");
  }
  double scattering_sphere_radius = 0.5;
  ierr = PetscOptionsGetScalar(NULL,"-scattering_sphere_radius",&scattering_sphere_radius,NULL); CHKERRQ(ierr);
    
  switch((AnalyticalSolutionTypes)choise_value) {
    
    case HARD_SPHERE_SCATTER_WAVE:
    
      {
	HardSphereScatterWave function_evaluator(wavenumber,scattering_sphere_radius);
	ierr = solve_problem(m_field,"EX1_PROBLEM","FE1","reEX","imEX",INSERT_VALUES,function_evaluator,is_partitioned); CHKERRQ(ierr);
      }
    
      break;
      
      
    case SOFT_SPHERE_SCATTER_WAVE:

      {
        SoftSphereScatterWave function_evaluator(wavenumber,scattering_sphere_radius);
        ierr = solve_problem(m_field,"EX1_PROBLEM","FE1","reEX","imEX",INSERT_VALUES,function_evaluator,is_partitioned); CHKERRQ(ierr);
      }

    break;

    case PLANE_WAVE:

      {
        PlaneWave function_evaluator(wavenumber,0.25*M_PI);
        ierr = solve_problem(m_field,"EX1_PROBLEM","FE1","reEX","imEX",INSERT_VALUES,function_evaluator,is_partitioned); CHKERRQ(ierr);
      }

    break;
      
    case HARD_CYLINDER_SCATTER_WAVE:
      
      {  
        HardCylinderScatterWave function_evaluator(wavenumber);
        ierr = solve_problem(m_field,"EX1_PROBLEM","FE1","reEX","imEX",INSERT_VALUES,function_evaluator,is_partitioned); CHKERRQ(ierr);
      }
      
    break;

    case SOFT_CYLINDER_SCATTER_WAVE:

      {  
        SoftCylinderScatterWave function_evaluator(wavenumber);
        ierr = solve_problem(m_field,"EX1_PROBLEM","FE1","reEX","imEX",INSERT_VALUES,function_evaluator,is_partitioned); CHKERRQ(ierr);
      }

    break;

    case INCIDENT_WAVE:

      {  
        IncidentWave function_evaluator(wavenumber,wave_direction);
        ierr = solve_problem(m_field,"EX1_PROBLEM","FE1","reEX","imEX",INSERT_VALUES,function_evaluator,is_partitioned); CHKERRQ(ierr);
      }

    break;

    default:

	SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA,"No analytical solution has been defined");

  }

  PetscBool add_incident_wave = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,"-add_incident_wave",&add_incident_wave,NULL); CHKERRQ(ierr);
  if(add_incident_wave) {

	//// define problem
    //ierr = m_field.add_problem("INCIDENT_WAVE"); CHKERRQ(ierr);
    //// set finite elements for problem
    //ierr = m_field.modify_problem_add_finite_element("INCIDENT_WAVE","FE1"); CHKERRQ(ierr);
    //// set refinment level for problem
    //ierr = m_field.modify_problem_ref_level_add_bit("INCIDENT_WAVE",bit_level0); CHKERRQ(ierr);
	
    //// build porblems
    //if(is_partitioned) {
    //  // if mesh is partitioned
    //  ierr = m_field.build_partitioned_problem("INCIDENT_WAVE",true); CHKERRQ(ierr);
    //  ierr = m_field.partition_finite_elements("INCIDENT_WAVE",true); CHKERRQ(ierr);
    //} else {
    //  // if not partitioned mesh is load to all processes 
    //  ierr = m_field.build_problem("INCIDENT_WAVE"); CHKERRQ(ierr);
    //  ierr = m_field.partition_problem("INCIDENT_WAVE"); CHKERRQ(ierr);
    //  ierr = m_field.partition_finite_elements("INCIDENT_WAVE"); CHKERRQ(ierr);
    //}
    //ierr = m_field.partition_ghost_dofs("INCIDENT_WAVE"); CHKERRQ(ierr);
	
    IncidentWave function_evaluator(wavenumber,wave_direction,1);
    ierr = solve_problem(m_field,"EX1_PROBLEM","FE1","reEX","imEX",ADD_VALUES,function_evaluator,is_partitioned); CHKERRQ(ierr);

  }
 
  if(is_partitioned) {
    rval = moab.write_file("best_solution.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
  } else {
    if(!pcomm->rank()) {
      rval = moab.write_file("analytical_solution.h5m"); CHKERR_PETSC(rval);
    }
  }

  
  PetscBool save_postproc_mesh = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,"-save_postproc_mesh",&save_postproc_mesh,NULL); CHKERRQ(ierr);
  if(save_postproc_mesh) {

    PostPocOnRefinedMesh post_proc(m_field);
    ierr = post_proc.generateRefereneElemenMesh(); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc("reEX"); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc("imEX"); CHKERRQ(ierr);

    if(m_field.check_field("rePRES") && m_field.check_field("imPRES")) {
      ierr = post_proc.addFieldValuesPostProc("rePRES"); CHKERRQ(ierr);
      ierr = post_proc.addFieldValuesPostProc("imPRES"); CHKERRQ(ierr);
    }

    ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.loop_finite_elements("EX1_PROBLEM","FE1",post_proc); CHKERRQ(ierr);
    rval = post_proc.postProcMesh.write_file("best_solution_mesh_post_proc.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);

  }
  
  // calulate total time
  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f S CPU Time = %f S \n",pcomm->rank(),v2-v1,t2-t1);
  
  ierr = PetscFinalize(); CHKERRQ(ierr);
  
  return 0;

}

 
