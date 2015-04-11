/* \file fe_approximation.cpp
 
  Calculates finite element (Galerkin) approximation for incident wave problem. 

  Note: 

  In this implementation, first pressure field is approximated on
  boundary and then finite element problem is solved. 

  For more rigorous convergence study, trace of best approximations on boundary
  can be calculated and then finite element for domain and Neumann/mix boundary.
  That will give exact pollution error.

 */

/* 
 * This file is part of MoFEM.
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. 
 *
 */
 
#include <MoFEM.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <DirichletBC.hpp>
#include <PotsProcOnRefMesh.hpp>

#include <Projection10NodeCoordsOnField.hpp>
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

#include <PCMGSetUpViaApproxOrders.hpp>
#include <AnalyticalSolutions.hpp>
#include <AnalyticalDirihlet.hpp>


#include <HelmholtzElement.hpp>

static char help[] = "...\n\n";
//argc = argument counts, argv = argument vectors
int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  //Core mb_instance;
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
  PetscBool is_partitioned = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-my_is_partitioned",&is_partitioned,&flg); CHKERRQ(ierr);
  if(is_partitioned) {
    //Read mesh to MOAB
    const char *option;
    option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  } else {
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  }

  // Create MoFEM (cephas) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  // Get start time for analyse
  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
  
  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("rePRES",H1,1); CHKERRQ(ierr);  
  ierr = m_field.add_field("imPRES",H1,1); CHKERRQ(ierr);
  
  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"rePRES"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"imPRES"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 2;
  }

  ierr = m_field.set_field_order(root_set,MBTET,"rePRES",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"rePRES",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"rePRES",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"rePRES",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTET,"imPRES",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"imPRES",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"imPRES",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"imPRES",1); CHKERRQ(ierr);
  
  if(!m_field.check_field("MESH_NODE_POSITIONS")) {
    ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

    ierr = m_field.build_fields(); CHKERRQ(ierr);
    Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
    ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  } else {

    ierr = m_field.build_fields(); CHKERRQ(ierr);

  }

  // Finite Elements

  HelmholtzElement helmholtz_element(m_field); 
  ierr = helmholtz_element.addHelmholtzElements("rePRES","imPRES"); CHKERRQ(ierr);

  if(m_field.check_field("reEX") && m_field.check_field("imEX")) {

    ierr = m_field.modify_finite_element_add_field_data("HELMHOLTZ_RERE_FE","reEX"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("HELMHOLTZ_RERE_FE","imEX"); CHKERRQ(ierr);

  }

  Range bc_tris;
  for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field,"ANALYTICAL_BC",it)) {
    rval = moab.get_entities_by_type(it->get_meshset(),MBTRI,bc_tris,true); CHKERR_PETSC(rval);
  }
  AnalyticalDirihletBC analytical_bc_real(m_field,bc_tris);
  AnalyticalDirihletBC analytical_bc_imag(m_field,bc_tris);
  ierr = analytical_bc_real.initializeProblem(m_field,"BCREAL_FE","rePRES"); CHKERRQ(ierr);
  ierr = analytical_bc_imag.initializeProblem(m_field,"BCIMAG_FE","imPRES"); CHKERRQ(ierr);

  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  // Problem
  ierr = m_field.add_problem("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.add_problem("BCREAL_PROBLEM"); CHKERRQ(ierr); //analytical Dirichlet for real field
  ierr = m_field.add_problem("BCIMAG_PROBLEM"); CHKERRQ(ierr); //analytical Dirichlet for imag field

  // Set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("ACOUSTIC_PROBLEM",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_add_bit("BCREAL_PROBLEM",bit_level0); CHKERRQ(ierr);  //analytical Dirichlet
  ierr = m_field.modify_problem_ref_level_add_bit("BCIMAG_PROBLEM",bit_level0); CHKERRQ(ierr);  //analytical Dirichlet

  // Add elements to problems
  ierr = m_field.modify_problem_add_finite_element("ACOUSTIC_PROBLEM","HELMHOLTZ_RERE_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ACOUSTIC_PROBLEM","HELMHOLTZ_IMIM_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ACOUSTIC_PROBLEM","HELMHOLTZ_REIM_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ACOUSTIC_PROBLEM","HELMHOLTZ_IMRE_FE"); CHKERRQ(ierr);

  ierr = m_field.modify_problem_add_finite_element("BCREAL_PROBLEM","BCREAL_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("BCIMAG_PROBLEM","BCIMAG_FE"); CHKERRQ(ierr);

  // Build problems

  // build porblems
  if(is_partitioned) {
    // if mesh is partitioned

    ierr = m_field.build_partitioned_problem("ACOUSTIC_PROBLEM",true); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("ACOUSTIC_PROBLEM",true); CHKERRQ(ierr);

    ierr = m_field.build_partitioned_problem("BCREAL_PROBLEM",true); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("BCREAL_PROBLEM",true); CHKERRQ(ierr);

    ierr = m_field.build_partitioned_problem("BCIMAG_PROBLEM",true); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("BCIMAG_PROBLEM",true); CHKERRQ(ierr);

  } else {
    // if not partitioned mesh is load to all processes 

    ierr = m_field.build_problem("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_problem("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);

    ierr = m_field.build_problem("BCREAL_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_problem("BCREAL_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("BCREAL_PROBLEM"); CHKERRQ(ierr);

    ierr = m_field.build_problem("BCIMAG_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_problem("BCIMAG_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("BCIMAG_PROBLEM"); CHKERRQ(ierr);

  }

  ierr = m_field.partition_ghost_dofs("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("BCREAL_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("BCIMAG_PROBLEM"); CHKERRQ(ierr);

  // Get problem matrices and vectors 

  Vec F;  //Right hand side vector
  ierr = m_field.VecCreateGhost("ACOUSTIC_PROBLEM",ROW,&F); CHKERRQ(ierr);
  Vec T; //Solution vector
  ierr = VecDuplicate(F,&T); CHKERRQ(ierr);
  Mat A; //Left hand side matrix
  ierr = m_field.MatCreateMPIAIJWithArrays("ACOUSTIC_PROBLEM",&A); CHKERRQ(ierr);

  // Solve for analytical Dirichlet bc dofs
  ierr = analytical_bc_real.setProblem(m_field,"BCREAL_PROBLEM"); CHKERRQ(ierr);
  ierr = analytical_bc_imag.setProblem(m_field,"BCIMAG_PROBLEM"); CHKERRQ(ierr);

  PetscInt choise_value = 0;
  // set type of analytical solution  
  ierr = PetscOptionsGetEList(NULL,"-analytical_solution_type",analytical_solution_types,2,&choise_value,NULL); CHKERRQ(ierr);

  double angularfreq = 1;
  double speed = 1; 

  /// this works only for one block 
  int nb_of_blocks = 0; 
  for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field,"MAT_HELMHOLTZ",it)) {

    //  get block attributes
    vector<double> attributes;
    ierr = it->get_Cubit_attributes(attributes); CHKERRQ(ierr);
    if(attributes.size()<2) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_INVALID_DATA,
	"not enough block attributes, expected 2 attributes ( angular freq., speed) , attributes.size() = %d ",attributes.size());
    }
    angularfreq = attributes[0];
    speed = attributes[1];  
    nb_of_blocks++;

  }
  
  if(nb_of_blocks!=1) {
    PetscPrintf(PETSC_COMM_SELF,"Warning: wave number is set to all blocks based on last evaluated block");
  }
  double wavenumber = angularfreq/speed;  

  // set wave number from line command, that overwrite numbre form block set
  ierr = PetscOptionsGetScalar(NULL,"-wave_number",&wavenumber,NULL); CHKERRQ(ierr);
  
  switch((AnalyticalSolutionTypes)choise_value) {

    case SPHERE_INCIDENT_WAVE:

      {
	SphereIncidentWave function_evaluator(wavenumber);
	ierr = analytical_bc_real.setApproxOps(m_field,"rePRES",function_evaluator,GenericAnalyticalSolution::REAL); CHKERRQ(ierr); 
	ierr = analytical_bc_imag.setApproxOps(m_field,"imPRES",function_evaluator,GenericAnalyticalSolution::IMAG); CHKERRQ(ierr);
      }

      break;

    case PLANE_WAVE:

      {
	PlaneWave function_evaluator(wavenumber,0.25*M_PI);
	ierr = analytical_bc_real.setApproxOps(m_field,"rePRES",function_evaluator,GenericAnalyticalSolution::REAL); CHKERRQ(ierr); 
	ierr = analytical_bc_imag.setApproxOps(m_field,"imPRES",function_evaluator,GenericAnalyticalSolution::IMAG); CHKERRQ(ierr);
      }

      break;

    case CYLINDER_INCIDENT_WAVE:

      {	
	CylinderIncidentWave function_evaluator(wavenumber);
	ierr = analytical_bc_real.setApproxOps(m_field,"rePRES",function_evaluator,GenericAnalyticalSolution::REAL); CHKERRQ(ierr); 
	ierr = analytical_bc_imag.setApproxOps(m_field,"imPRES",function_evaluator,GenericAnalyticalSolution::IMAG); CHKERRQ(ierr);
      }

      break;

  }

  // Analytical boundary conditions
  AnalyticalDirihletBC::DirichletBC analytical_ditihlet_bc_real(m_field,"rePRES",A,T,F);
  AnalyticalDirihletBC::DirichletBC analytical_ditihlet_bc_imag(m_field,"imPRES",A,T,F);

  ierr = analytical_bc_real.solveProblem(m_field,"BCREAL_PROBLEM","BCREAL_FE",analytical_ditihlet_bc_real); CHKERRQ(ierr);
  ierr = analytical_bc_imag.solveProblem(m_field,"BCIMAG_PROBLEM","BCIMAG_FE",analytical_ditihlet_bc_imag); CHKERRQ(ierr);  

  ierr = analytical_bc_real.destroyProblem(); CHKERRQ(ierr);
  ierr = analytical_bc_imag.destroyProblem(); CHKERRQ(ierr);

  // Zero vectors
  ierr = VecZeroEntries(T); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  // Assemble problem
  ierr = m_field.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc_real); CHKERRQ(ierr);
  ierr = m_field.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc_imag); CHKERRQ(ierr);

  ierr = helmholtz_element.setOperators(A,F,"rePRES","imPRES"); CHKERRQ(ierr);
  ierr = helmholtz_element.calculateA("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
  ierr = helmholtz_element.calculateF("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);

  ierr = m_field.problem_basic_method_postProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc_real); CHKERRQ(ierr);
  ierr = m_field.problem_basic_method_postProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc_imag); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecScale(F,-1); CHKERRQ(ierr);

  // Solve problem
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,T); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  if(is_partitioned) {

    // no need for global communication
    ierr = m_field.set_local_ghost_vector("ACOUSTIC_PROBLEM",ROW,T,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);  

  } else {

    ierr = m_field.set_global_ghost_vector("ACOUSTIC_PROBLEM",ROW,T,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);  

  }

  // Destroy the KSP solvers
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&T); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);

  if(is_partitioned) {
    rval = moab.write_file("fe_solution.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
  } else {
    if(!pcomm->rank()) {
      rval = moab.write_file("fe_solution.h5m"); CHKERR_PETSC(rval);
    }
  }
 
  PetscBool save_postproc_mesh = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL,"-save_postproc_mesh",&save_postproc_mesh,NULL); CHKERRQ(ierr);
  if(save_postproc_mesh) {

    PostPocOnRefinedMesh post_proc(m_field);
    ierr = post_proc.generateRefereneElemenMesh(); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc("rePRES"); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc("imPRES"); CHKERRQ(ierr);

    if(m_field.check_field("reEX") && m_field.check_field("imEX")) {
      ierr = post_proc.addFieldValuesPostProc("reEX"); CHKERRQ(ierr);
      ierr = post_proc.addFieldValuesPostProc("imEX"); CHKERRQ(ierr);
    }

    ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_RERE_FE",post_proc); CHKERRQ(ierr);
    rval = post_proc.postProcMesh.write_file("fe_solution_mesh_post_proc.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);

  }
  
  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f S CPU Time = %f S \n",pcomm->rank(),v2-v1,t2-t1);
   
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}

