/** \file fe_approximation.cpp

  Calculates finite element (Galerkin) approximation for incident wave problem.

  Note:

  In this implementation, first acoustic potential field is approximated on
  boundary and then finite element problem is solved.  For more rigorous
  convergence study, trace of best approximations on boundary can be calculated
  and then finite element for domain and Neumann/mix boundary.  That will give
  exact pollution error.

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
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <DirichletBC.hpp>
#include <PostProcOnRefMesh.hpp>

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

#include <FieldApproximation.hpp>

using namespace std;
using namespace boost::math;

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;


#include <boost/shared_array.hpp>
#include <kiss_fft.h>
#include <kiss_fft.c>

#include <AnalyticalSolutions.hpp>
#include <AnalyticalDirichlet.hpp>
#include <HelmholtzElement.hpp>
#include <TimeSeries.hpp>

struct PlaneIncidentWaveSacttrerData {

  Range tRis;

};

static char help[] = "...\n\n";

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
    SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  // Create moab parallel communicator
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
  ierr = m_field.add_field("P",H1,1); CHKERRQ(ierr);  // in time domain

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"rePRES"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"imPRES"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"P"); CHKERRQ(ierr);

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

  ierr = m_field.set_field_order(root_set,MBTET,"P",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"P",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"P",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"P",1); CHKERRQ(ierr);

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
  ierr = helmholtz_element.getGlobalParametersFromLineCommandOptions(); CHKERRQ(ierr);
  ierr = helmholtz_element.addHelmholtzElements("rePRES","imPRES"); CHKERRQ(ierr);
  if(m_field.check_field("reEX") && m_field.check_field("imEX")) {
    ierr = m_field.modify_finite_element_add_field_data("HELMHOLTZ_RERE_FE","reEX"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("HELMHOLTZ_RERE_FE","imEX"); CHKERRQ(ierr);
  }

  bool Dirichlet_bc_set = false;
  Range bc_dirichlet_tris,analytical_bc_tris;
  for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field,"ANALYTICAL_BC",it)) {
    rval = moab.get_entities_by_type(it->get_meshset(),MBTRI,analytical_bc_tris,true); CHKERR_PETSC(rval);
    Dirichlet_bc_set = true;
  }
  bc_dirichlet_tris.merge(analytical_bc_tris);
  AnalyticalDirichletBC analytical_bc_real(m_field);
  AnalyticalDirichletBC analytical_bc_imag(m_field);
  ierr = analytical_bc_real.initializeProblem(m_field,"BCREAL_FE","rePRES",analytical_bc_tris); CHKERRQ(ierr);
  ierr = analytical_bc_imag.initializeProblem(m_field,"BCIMAG_FE","imPRES",analytical_bc_tris); CHKERRQ(ierr);

  PetscBool wavenumber_flg;
  double wavenumber;
  // set wave number from line command, that overwrite numbre form block set
  ierr = PetscOptionsGetScalar(NULL,"-wave_number",&wavenumber,&wavenumber_flg); CHKERRQ(ierr);
  if(!wavenumber_flg) {

    SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA,"wave number not given, set in line command -wave_number to fix problem");

  }

  double power_of_incident_wave = 1;
  ierr = PetscOptionsGetScalar(NULL,"-power_of_incident_wave",&power_of_incident_wave,NULL); CHKERRQ(ierr);


  // This is added for a case than on some surface, defined by the user a
  // incident plane wave is scattered.
  map<int,PlaneIncidentWaveSacttrerData> planeWaveScatterData;
  for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field,"SOFT_INCIDENT_WAVE_BC",it)) {

    rval = moab.get_entities_by_type(it->get_meshset(),MBTRI,planeWaveScatterData[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
    ierr = analytical_bc_real.initializeProblem(m_field,"BCREAL_FE","rePRES",planeWaveScatterData[it->get_msId()].tRis); CHKERRQ(ierr);
    ierr = analytical_bc_imag.initializeProblem(m_field,"BCIMAG_FE","imPRES",planeWaveScatterData[it->get_msId()].tRis); CHKERRQ(ierr);
    bc_dirichlet_tris.merge(planeWaveScatterData[it->get_msId()].tRis);

    Dirichlet_bc_set = true;

  }

  //ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //// Build adjacencies
  //ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

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

  ierr = m_field.modify_problem_add_finite_element("BCREAL_PROBLEM","BCREAL_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("BCIMAG_PROBLEM","BCIMAG_FE"); CHKERRQ(ierr);


  // Build problems
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  // Build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  // build porblems
  if(is_partitioned) {
    // if mesh is partitioned

    ierr = m_field.build_partitioned_problem("ACOUSTIC_PROBLEM",true); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("ACOUSTIC_PROBLEM",true); CHKERRQ(ierr);

    if(Dirichlet_bc_set) {
      ierr = m_field.build_partitioned_problem("BCREAL_PROBLEM",true); CHKERRQ(ierr);
      ierr = m_field.partition_finite_elements("BCREAL_PROBLEM",true); CHKERRQ(ierr);

      ierr = m_field.build_partitioned_problem("BCIMAG_PROBLEM",true); CHKERRQ(ierr);
      ierr = m_field.partition_finite_elements("BCIMAG_PROBLEM",true); CHKERRQ(ierr);

    }


  } else {
    // if not partitioned mesh is load to all processes

    ierr = m_field.build_problem("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_problem("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);

    if(Dirichlet_bc_set) {
      ierr = m_field.build_problem("BCREAL_PROBLEM"); CHKERRQ(ierr);
      ierr = m_field.partition_problem("BCREAL_PROBLEM"); CHKERRQ(ierr);
      ierr = m_field.partition_finite_elements("BCREAL_PROBLEM"); CHKERRQ(ierr);



      ierr = m_field.build_problem("BCIMAG_PROBLEM"); CHKERRQ(ierr);
      ierr = m_field.partition_problem("BCIMAG_PROBLEM"); CHKERRQ(ierr);
      ierr = m_field.partition_finite_elements("BCIMAG_PROBLEM"); CHKERRQ(ierr);
    }
  }

  ierr = m_field.partition_ghost_dofs("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
  if(Dirichlet_bc_set) {
    ierr = m_field.partition_ghost_dofs("BCREAL_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_ghost_dofs("BCIMAG_PROBLEM"); CHKERRQ(ierr);
  }

  // Get problem matrices and vectors
  Vec F;  //Right hand side vector
  ierr = m_field.VecCreateGhost("ACOUSTIC_PROBLEM",ROW,&F); CHKERRQ(ierr);
  Vec T; //Solution vector
  ierr = VecDuplicate(F,&T); CHKERRQ(ierr);
  Mat A; //Left hand side matrix
  ierr = m_field.MatCreateMPIAIJWithArrays("ACOUSTIC_PROBLEM",&A); CHKERRQ(ierr);
  ierr = helmholtz_element.setOperators(A,F,"rePRES","imPRES"); CHKERRQ(ierr);

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

  PetscInt choise_value = NO_ANALYTICAL_SOLUTION;
  // set type of analytical solution
  ierr = PetscOptionsGetEList(NULL,"-analytical_solution_type",analytical_solution_types,6,&choise_value,NULL); CHKERRQ(ierr);

  switch((AnalyticalSolutionTypes)choise_value) {

    case HARD_SPHERE_SCATTER_WAVE:

    {
      double scattering_sphere_radius = 0.5;
      ierr = PetscOptionsGetScalar(NULL,"-scattering_sphere_radius",&scattering_sphere_radius,NULL); CHKERRQ(ierr);



      boost::shared_ptr<HardSphereScatterWave> function_evaluator = boost::shared_ptr<HardSphereScatterWave>(
        new HardSphereScatterWave(wavenumber,scattering_sphere_radius)
      );
      ierr = analytical_bc_real.setApproxOps(
        m_field,"rePRES",analytical_bc_tris,function_evaluator,GenericAnalyticalSolution::REAL
      ); CHKERRQ(ierr);
      ierr = analytical_bc_imag.setApproxOps(
        m_field,"imPRES",analytical_bc_tris,function_evaluator,GenericAnalyticalSolution::IMAG
      ); CHKERRQ(ierr);
      Dirichlet_bc_set = true;

    }

    break;

    case SOFT_SPHERE_SCATTER_WAVE:

      {

      double scattering_sphere_radius = 0.5;
      ierr = PetscOptionsGetScalar(NULL,"-scattering_sphere_radius",&scattering_sphere_radius,NULL); CHKERRQ(ierr);


      boost::shared_ptr<SoftSphereScatterWave> function_evaluator = boost::shared_ptr<SoftSphereScatterWave>(new SoftSphereScatterWave(wavenumber,scattering_sphere_radius));
      ierr = analytical_bc_real.setApproxOps(m_field,"rePRES",analytical_bc_tris,function_evaluator,GenericAnalyticalSolution::REAL); CHKERRQ(ierr);
      ierr = analytical_bc_imag.setApproxOps(m_field,"imPRES",analytical_bc_tris,function_evaluator,GenericAnalyticalSolution::IMAG); CHKERRQ(ierr);
      Dirichlet_bc_set = true;

    }

    break;

    case PLANE_WAVE:

    {

      double angle = 0.25;
      // set wave number from line command, that overwrite numbre form block set
      ierr = PetscOptionsGetScalar(NULL,"-wave_guide_angle",&angle,NULL); CHKERRQ(ierr);


      boost::shared_ptr<PlaneWave> function_evaluator = boost::shared_ptr<PlaneWave>(new PlaneWave(wavenumber,angle*M_PI));
      ierr = analytical_bc_real.setApproxOps(m_field,"rePRES",analytical_bc_tris,function_evaluator,GenericAnalyticalSolution::REAL); CHKERRQ(ierr);
      ierr = analytical_bc_imag.setApproxOps(m_field,"imPRES",analytical_bc_tris,function_evaluator,GenericAnalyticalSolution::IMAG); CHKERRQ(ierr);
      Dirichlet_bc_set = true;

    }

    break;

    case HARD_CYLINDER_SCATTER_WAVE:

    {

      boost::shared_ptr<HardCylinderScatterWave> function_evaluator = boost::shared_ptr<HardCylinderScatterWave>(new HardCylinderScatterWave(wavenumber));
      ierr = analytical_bc_real.setApproxOps(m_field,"rePRES",analytical_bc_tris,function_evaluator,GenericAnalyticalSolution::REAL); CHKERRQ(ierr);
      ierr = analytical_bc_imag.setApproxOps(m_field,"imPRES",analytical_bc_tris,function_evaluator,GenericAnalyticalSolution::IMAG); CHKERRQ(ierr);
      Dirichlet_bc_set = true;


    }

    break;

    case SOFT_CYLINDER_SCATTER_WAVE:

    {

      boost::shared_ptr<SoftCylinderScatterWave> function_evaluator = boost::shared_ptr<SoftCylinderScatterWave>(new SoftCylinderScatterWave(wavenumber));
      ierr = analytical_bc_real.setApproxOps(m_field,"rePRES",analytical_bc_tris,function_evaluator,GenericAnalyticalSolution::REAL); CHKERRQ(ierr);
      ierr = analytical_bc_imag.setApproxOps(m_field,"imPRES",analytical_bc_tris,function_evaluator,GenericAnalyticalSolution::IMAG); CHKERRQ(ierr);
      Dirichlet_bc_set = true;

    }

    break;

    case INCIDENT_WAVE:

    {

      boost::shared_ptr<IncidentWave> function_evaluator =
      boost::shared_ptr<IncidentWave>(new IncidentWave(wavenumber,wave_direction,power_of_incident_wave));
      ierr = analytical_bc_real.setApproxOps(m_field,"rePRES",analytical_bc_tris,function_evaluator,GenericAnalyticalSolution::REAL); CHKERRQ(ierr);
      ierr = analytical_bc_imag.setApproxOps(m_field,"imPRES",analytical_bc_tris,function_evaluator,GenericAnalyticalSolution::IMAG); CHKERRQ(ierr);
      Dirichlet_bc_set = true;

    }

    break;

    case NO_ANALYTICAL_SOLUTION:

    {
      Dirichlet_bc_set = false;
    }

    break;

  }

  // Analytical boundary conditions
  AnalyticalDirichletBC::DirichletBC analytical_ditihlet_bc_real(m_field,"rePRES",A,T,F);
  AnalyticalDirichletBC::DirichletBC analytical_ditihlet_bc_imag(m_field,"imPRES",A,T,F);

  if(Dirichlet_bc_set) {

    {

      map<int,PlaneIncidentWaveSacttrerData>::iterator mit = planeWaveScatterData.begin();
      for(;mit!=planeWaveScatterData.end();mit++) {

        // note negative field, scatter field should cancel incident wave
        boost::shared_ptr<IncidentWave> function_evaluator = boost::shared_ptr<IncidentWave>(
          new IncidentWave(wavenumber,wave_direction,-power_of_incident_wave)
        );
        ierr = analytical_bc_real.setApproxOps(
          m_field,"rePRES",mit->second.tRis,function_evaluator,GenericAnalyticalSolution::REAL
        ); CHKERRQ(ierr);
        ierr = analytical_bc_imag.setApproxOps(
          m_field,"imPRES",mit->second.tRis,function_evaluator,GenericAnalyticalSolution::IMAG
        ); CHKERRQ(ierr);

      }

    }
    // Solve for analytical Dirichlet bc dofs
    ierr = analytical_bc_real.setProblem(m_field,"BCREAL_PROBLEM"); CHKERRQ(ierr);
    ierr = analytical_bc_imag.setProblem(m_field,"BCIMAG_PROBLEM"); CHKERRQ(ierr);
    ierr = analytical_bc_real.solveProblem(
      m_field,"BCREAL_PROBLEM","BCREAL_FE",analytical_ditihlet_bc_real,bc_dirichlet_tris
    ); CHKERRQ(ierr);
    ierr = analytical_bc_imag.solveProblem(
      m_field,"BCIMAG_PROBLEM","BCIMAG_FE",analytical_ditihlet_bc_imag,bc_dirichlet_tris
    ); CHKERRQ(ierr);

    ierr = analytical_bc_real.destroyProblem(); CHKERRQ(ierr);
    ierr = analytical_bc_imag.destroyProblem(); CHKERRQ(ierr);

  }


  PetscBool monochromatic_wave = PETSC_TRUE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-monochromatic_wave",&monochromatic_wave,NULL); CHKERRQ(ierr);

  PetscBool add_incident_wave = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,"-add_incident_wave",&add_incident_wave,NULL); CHKERRQ(ierr);
  if(add_incident_wave) {

    // define problem
    ierr = m_field.add_problem("INCIDENT_WAVE"); CHKERRQ(ierr);
    // set finite elements for problem
    ierr = m_field.modify_problem_add_finite_element("INCIDENT_WAVE","HELMHOLTZ_RERE_FE"); CHKERRQ(ierr);
    // set refinment level for problem
    ierr = m_field.modify_problem_ref_level_add_bit("INCIDENT_WAVE",bit_level0); CHKERRQ(ierr);

    // build porblems
    if(is_partitioned) {
      // if mesh is partitioned
      ierr = m_field.build_partitioned_problem("INCIDENT_WAVE",true); CHKERRQ(ierr);
      ierr = m_field.partition_finite_elements("INCIDENT_WAVE",true); CHKERRQ(ierr);
    } else {
      // if not partitioned mesh is load to all processes
      ierr = m_field.build_problem("INCIDENT_WAVE"); CHKERRQ(ierr);
      ierr = m_field.partition_problem("INCIDENT_WAVE"); CHKERRQ(ierr);
      ierr = m_field.partition_finite_elements("INCIDENT_WAVE"); CHKERRQ(ierr);
      }
    ierr = m_field.partition_ghost_dofs("INCIDENT_WAVE"); CHKERRQ(ierr);

  }


  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);

  if(monochromatic_wave) {

    // Zero vectors
    ierr = VecZeroEntries(T); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = MatZeroEntries(A); CHKERRQ(ierr);

    // Assemble problem
    if(Dirichlet_bc_set) {
      ierr = m_field.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc_real); CHKERRQ(ierr);
      ierr = m_field.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc_imag); CHKERRQ(ierr);
    }

    ierr = helmholtz_element.calculateA("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
    ierr = helmholtz_element.calculateF("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);

    if(Dirichlet_bc_set) {
      ierr = m_field.problem_basic_method_postProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc_real); CHKERRQ(ierr);
      ierr = m_field.problem_basic_method_postProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc_imag); CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = VecScale(F,-1); CHKERRQ(ierr);

    // Solve problem
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

  } else {


    // define problem
    ierr = m_field.add_problem("PRESSURE_IN_TIME"); CHKERRQ(ierr);
    // set finite elements for problem
    ierr = m_field.modify_problem_add_finite_element("PRESSURE_IN_TIME","PRESSURE_FE"); CHKERRQ(ierr);
    // set refinment level for problem
    ierr = m_field.modify_problem_ref_level_add_bit("PRESSURE_IN_TIME",bit_level0); CHKERRQ(ierr);

    // build porblems
    if(is_partitioned) {
      // if mesh is partitioned
      ierr = m_field.build_partitioned_problem("PRESSURE_IN_TIME",true); CHKERRQ(ierr);
      ierr = m_field.partition_finite_elements("PRESSURE_IN_TIME",true); CHKERRQ(ierr);
    } else {
      // if not partitioned mesh is load to all processes
      ierr = m_field.build_problem("PRESSURE_IN_TIME"); CHKERRQ(ierr);
      ierr = m_field.partition_problem("PRESSURE_IN_TIME"); CHKERRQ(ierr);
      ierr = m_field.partition_finite_elements("PRESSURE_IN_TIME"); CHKERRQ(ierr);
    }
    ierr = m_field.partition_ghost_dofs("PRESSURE_IN_TIME"); CHKERRQ(ierr);

    TimeSeries time_series(m_field,helmholtz_element,
      analytical_ditihlet_bc_real,analytical_ditihlet_bc_imag,
      Dirichlet_bc_set);

    ierr = time_series.readData(); CHKERRQ(ierr);
    ierr = time_series.createPressureSeries(T); CHKERRQ(ierr);
    ierr = time_series.forwardSpaceDft(); CHKERRQ(ierr);
    ierr = time_series.pressureForwardDft(); CHKERRQ(ierr);
    ierr = time_series.solveForwardDFT(solver,A,F,T); CHKERRQ(ierr);
    ierr = time_series.pressureInverseDft(); CHKERRQ(ierr);
    ierr = time_series.generateReferenceElementMesh(); CHKERRQ(ierr);
    ierr = time_series.saveResults(); CHKERRQ(ierr);
    ierr = time_series.destroyPressureSeries(); CHKERRQ(ierr);

  }


  //Vec P,M;
  //ierr = m_field.VecCreateGhost("EX1_PROBLEM",COL,&M); CHKERRQ(ierr);
  //ierr = VecDuplicate(M,&P); CHKERRQ(ierr);

  //ierr = m_field.set_local_ghost_vector("EX1_PROBLEM",COL,M,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //ierr = m_field.set_other_global_ghost_vector("EX1_PROBLEM","reEX","imEX",COL,P,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //double nrm2_M;
  //ierr = VecNorm(M,NORM_2,&nrm2_M);  CHKERRQ(ierr);

  //Vec V;
  //ierr = m_field.VecCreateGhost("ACOUSTIC_PROBLEM",COL,&V); CHKERRQ(ierr);
  //ierr = m_field.set_local_ghost_vector("ACOUSTIC_PROBLEM",COL,V,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //VecScatter scatter_real,scatter_imag;

  //ierr = m_field.VecScatterCreate(V,"ACOUSTIC_PROBLEM","rePRES",COL,M,"EX1_PROBLEM","reEX",COL,&scatter_real); CHKERRQ(ierr);

  //ierr = m_field.VecScatterCreate(V,"ACOUSTIC_PROBLEM","imPRES",COL,P,"EX1_PROBLEM","reEX",COL,&scatter_imag); CHKERRQ(ierr);

  //VecScale(V,-1);

  //ierr = VecScatterBegin(scatter_real,V,M,ADD_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //ierr = VecScatterEnd(scatter_real,V,M,ADD_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //double nrm2_ErrM;
  //ierr = VecNorm(M,NORM_2,&nrm2_ErrM);  CHKERRQ(ierr);
  //PetscPrintf(PETSC_COMM_WORLD,"L2 relative error on real field of acoustic problem %6.4e\n",(nrm2_ErrM)/(nrm2_M));


  //ierr = VecDestroy(&M); CHKERRQ(ierr);
  //ierr = VecDestroy(&P); CHKERRQ(ierr);
  //ierr = VecDestroy(&V); CHKERRQ(ierr);

  // Destroy the KSP solvers
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&T); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  if(monochromatic_wave) {

    if(add_incident_wave) {

      IncidentWave function_evaluator(wavenumber,wave_direction,power_of_incident_wave);

      ierr = solve_problem(m_field,"INCIDENT_WAVE","HELMHOLTZ_RERE_FE","rePRES","imPRES",
        ADD_VALUES,function_evaluator,is_partitioned); CHKERRQ(ierr);

    }

    PetscBool save_postproc_mesh = PETSC_TRUE;
    ierr = PetscOptionsGetBool(NULL,"-save_postproc_mesh",&save_postproc_mesh,NULL); CHKERRQ(ierr);
    if(save_postproc_mesh) {

      PostPocOnRefinedMesh post_proc(m_field);
      ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
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
    if(is_partitioned) {
      rval = moab.write_file("fe_solution.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
    } else {

      if(!pcomm->rank()) {
        rval = moab.write_file("fe_solution.h5m"); CHKERR_PETSC(rval);
      }

    }

  }

  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f S CPU Time = %f S \n",pcomm->rank(),v2-v1,t2-t1);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
