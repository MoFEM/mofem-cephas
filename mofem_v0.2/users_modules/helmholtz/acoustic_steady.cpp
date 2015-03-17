/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
 */

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *PhD student Thomas Xuan Meng/Users/xuanmeng/Documents/mofem-cephas/mofem_v0.2/users_modules/helmholtz/README
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. 
 * abs = sqrt(rePRES^2+imPRES^2) */
 

#include <MoFEM.hpp>


#include <DirichletBC.hpp>
#include <PotsProcOnRefMesh.hpp>
#include <HelmholtzElement.hpp>

#include <Projection10NodeCoordsOnField.hpp>
#include <AnalyticalDirichletHelmholtz.hpp>
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
void error1( const string& msg )
{
	throw( runtime_error( msg ) );
}


namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

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

  bool useImpedance;
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  
  char impedance[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-use_impedance",impedance,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
	  SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -use_impedance (true of false needed)");
  }
  if (strcmp ("true",impedance ) == 0) {useImpedance = true;}
  else if(strcmp ("false",impedance ) == 0) {useImpedance = false;}
  
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  BARRIER_RANK_START(pcomm) 
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  BARRIER_RANK_END(pcomm) 

  //Create MoFEM (cephas) database
  //FieldCore core(moab);
  MoFEM::Core core(moab);
  FieldInterface& mField = core;

  //count the comsumption of time by single run
  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
  
  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_field("rePRES",H1,1); CHKERRQ(ierr);  //field order distinguish the scalar field and vector field.
  ierr = mField.add_field("imPRES",H1,1); CHKERRQ(ierr);
  
  //Problem
  ierr = mField.add_problem("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.add_problem("BC1_PROBLEM"); CHKERRQ(ierr); //analytical Dirichlet for real field
  ierr = mField.add_problem("BC2_PROBLEM"); CHKERRQ(ierr); //analytical Dirichlet for imag field
  
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ACOUSTIC_PROBLEM",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("BC1_PROBLEM",bit_level0); CHKERRQ(ierr);  //analytical Dirichlet
  ierr = mField.modify_problem_ref_level_add_bit("BC2_PROBLEM",bit_level0); CHKERRQ(ierr);  //analytical Dirichlet
  
  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"rePRES"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"imPRES"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 2;
  }

  ierr = mField.set_field_order(root_set,MBTET,"rePRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"rePRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"rePRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"rePRES",1); CHKERRQ(ierr);
  
  if(!mField.check_field("MESH_NODE_POSITIONS")) {
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  }
  ierr = mField.set_field_order(root_set,MBTET,"imPRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"imPRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"imPRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"imPRES",1); CHKERRQ(ierr);

  HelmholtzElement helmholtz_elements(mField); //Create the HelmholtzElement class in the header-file
  
  ierr = helmholtz_elements.addHelmholtzElements("ACOUSTIC_PROBLEM","rePRES","imPRES"); CHKERRQ(ierr);
  ierr = helmholtz_elements.addHelmholtzFluxElement("ACOUSTIC_PROBLEM","rePRES","imPRES"); CHKERRQ(ierr);
  if(useImpedance) {
  ierr = helmholtz_elements.addHelmholtzImpedanceElement("ACOUSTIC_PROBLEM","rePRES","imPRES"); CHKERRQ(ierr);
  }
  
  
  //Set up the analytical Dirichlet BC
  Range bc_tris;
  for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"ANALYTICAL_BC",it)) {
   rval = moab.get_entities_by_type(it->get_meshset(),MBTRI,bc_tris,true); CHKERR_PETSC(rval);
  }
  
  AnalyticalDirihletBC analytical_bc1(mField,bc_tris);
  AnalyticalDirihletBC analytical_bc2(mField,bc_tris);
  ierr = analytical_bc1.initializeBcProblem(mField,"BC1_PROBLEM","BC1_FE","rePRES"); CHKERRQ(ierr);
  ierr = analytical_bc2.initializeBcProblem(mField,"BC2_PROBLEM","BC2_FE","imPRES"); CHKERRQ(ierr);
  /*** add exact solution data in finite element */
  if(mField.check_field("reEX") && mField.check_field("imEX")) {
	  ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_FE","reEX"); CHKERRQ(ierr);
	  ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_FE","imEX"); CHKERRQ(ierr);
  }
  //End of Dirichlet set up
  

  /****/
  //build database
  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  Projection10NodeCoordsOnField ent_method_material(mField,"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  /****/
  //mesh partitioning 
  //partition
  ierr = mField.partition_problem("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
  
  //mesh partitinoning for analytical Dirichlet
  ierr = mField.simple_partition_problem("BC1_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("BC1_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.simple_partition_problem("BC2_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("BC2_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("BC1_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("BC2_PROBLEM"); CHKERRQ(ierr);

  Vec F;  //Right hand side vector
  ierr = mField.VecCreateGhost("ACOUSTIC_PROBLEM",ROW,&F); CHKERRQ(ierr);
  Vec T; //Solution vector
  ierr = VecDuplicate(F,&T); CHKERRQ(ierr);
  Mat A; //Left hand side matrix
  ierr = mField.MatCreateMPIAIJWithArrays("ACOUSTIC_PROBLEM",&A); CHKERRQ(ierr);
    
  bool useScalar = true;
  //ierr = helmholtz_elements.setHelmholtzFiniteElementRhs_FOperators("rePRES","rePRES",F,useScalar); CHKERRQ(ierr); //The analytical F source vector
  ierr = helmholtz_elements.setHelmholtzFiniteElementRhsOperators("rePRES","imPRES",F,useImpedance); CHKERRQ(ierr); //the Rhs of Dirichlet BC
  //ierr = helmholtz_elements.setHelmholtzFiniteElementRhsOperators_imim("imPRES","imPRES",F); CHKERRQ(ierr); //the Rhs residual of Dirichlet BC
  ierr = helmholtz_elements.setHelmholtzFiniteElementLhsOperators("rePRES","imPRES",(A)); CHKERRQ(ierr);//Stiffness Matrix
  ierr = helmholtz_elements.setHelmholtzMassFiniteElementLhsOperators("rePRES","imPRES",(A)); CHKERRQ(ierr);//Mass Matrix
  //ierr = helmholtz_elements.setHelmholtzFluxFiniteElementRhsOperators("rePRES","rePRES",F); CHKERRQ(ierr);  //real Neumann BC
  //ierr = helmholtz_elements.setHelmholtzFluxFiniteElementRhsOperators("imPRES","imPRES",F); CHKERRQ(ierr);    //Imag Neumann BC
  ierr = helmholtz_elements.setHelmholtzIncidentWaveFiniteElementRhsOperators("rePRES","imPRES",F); CHKERRQ(ierr); // Incident wave flux.
  if(useImpedance) {
	//The boundary Impedance BC.
	ierr = helmholtz_elements.setHelmholtzImpedanceFiniteElementLhsOperators("rePRES","imPRES",(A)); CHKERRQ(ierr);
  }

  
  ierr = VecZeroEntries(T); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);
  
  //analytical dirichlet bc
  AnalyticalDirihletBC::DirichletBC analytical_ditihlet_bc1(mField,"rePRES",A,T,F);
  AnalyticalDirihletBC::DirichletBC analytical_ditihlet_bc2(mField,"imPRES",A,T,F);
  
  //solve for analytical dirichlet bc dofs
  ierr = analytical_bc1.setBcProblem(mField,"BC1_PROBLEM"); CHKERRQ(ierr);
  ierr = analytical_bc2.setBcProblem(mField,"BC2_PROBLEM"); CHKERRQ(ierr);
  
  static double aNgularfreq;
  static double sPeed; //Without static. got error:use of local variable with automatic storage from containing function

  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it))
  {
	  cout << endl << *it << endl;
	  

	  if(it->get_Cubit_name().compare(0,13,"MAT_HELMHOLTZ") == 0) {
		  
		  //get block attributes
		  vector<double> attributes;
		  ierr = it->get_Cubit_attributes(attributes); CHKERRQ(ierr);
		  if(attributes.size()<2) {
			  SETERRQ1(PETSC_COMM_SELF,1,"not enough block attributes to deffine fluid pressure element, attributes.size() = %d ",attributes.size());
		  }
		  
		  aNgularfreq = attributes[0];
		  sPeed = attributes[1];

	  }
  }
  
  //Extract the data output to .txt file.
  std::string filename("scattered_sphere_outputs.txt" );
  static ofstream ofs( filename.c_str() ); //put the data from cpu into file
  if( !ofs ){
  error1( "Error opening file" );}
  ofs.precision( 18 );
  cout.precision( 18 );

  /* this function compute the scattered field of helmholtz operator */
  struct AnaliticalFunction {
	  static double fUN(double x,double y,double z,bool use_real) {
		  
		  bool useReal;
		  
		  const double pi = atan( 1.0 ) * 4.0;
		  double R = sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)); //radius
		  double theta = atan2(y,x)+2*pi; //the arctan of radians (y/x)
		  		 
		  	  if(theta != theta) cerr << "theta\n";
		  
		  const double wAvenumber = aNgularfreq/sPeed;
		  
		  const double k = wAvenumber;  //Wave number
		  const double a = 0.5;         //radius of the sphere,wait to modify by user
		  const double const1 = k * a;
		  double const2 = k * R;
		  
		  
		  const complex< double > i( 0.0, 1.0 );
		  
		  // magnitude of incident wave
		  const double phi_incident_mag = 1.0;
		  
		  const double tol = 1.0e-6;
		  double max = 0.0;
		  double min = 999999.0;
		  
		  complex< double > result = 0.0;
		  complex< double > prev_result;
		  
		  double error = 100.0;
		  unsigned int n = 0; //initialized the infinite series loop
		  
		  while( error > tol )  //finding the acoustic potential in one single point.
		  {
			  double jn_der = n / const1 * sph_bessel( n, const1 ) - sph_bessel( n + 1, const1 );  //The derivative of bessel function
			  if(jn_der != jn_der) cerr << "error jn_der\n";
			  
			  complex< double > hn_der = n / const1 * sph_hankel_1( n, const1 ) - sph_hankel_1( n + 1, const1 );
			  //if(hn_der != hn_der) cerr << "error hn_der\n";
			  //complex< double > hn_der = 0.5 * ( sph_hankel_1( n - 1, const1 ) -
			  //( sph_hankel_1( n, const1 ) + const1 * sph_hankel_1( n + 1, const1 ) ) / const1 );
			  double Pn = legendre_p( n, cos( theta ) );
			  if(Pn != Pn) cerr << "Pn \n";
			  complex< double >hn = sph_hankel_1( n, const2 );  //S Hankel first kind function
			  if(const2 != const2) cerr << "const2 \n";
			  
			  if(n == 0) { complex< double > hn_c = -i*exp(i*const2)*(1/const2); cout << "\n hn_c = \n" << hn_c << endl;}
			  if(hn != hn) {cerr << "hn \n"; cout << hn << "\n n = \n" << n << endl; cout << "\n k * r = \n" << const2 << endl;}
			  prev_result = result;
			  result -= pow( i, n ) * ( 2.0 * n + 1.0 ) * jn_der / hn_der * Pn * hn;
			  error = abs( abs( result ) - abs( prev_result ) );
			  ++n;
		  }
          
		  //const complex< double > inc_field = exp( i * k * R * cos( theta ) );  //incident wave
		  //const complex< double > total_field = inc_field + result;
		  
		  //double val = std::real(result);
		  //ofs << theta << "\t" << abs( result ) << "\t" << abs( inc_field ) << "\t" << abs( total_field ) << "\t" << R << endl; //write the file
		  
		  ///* cube 2D */
		  //double theta = pi/4;
		  //result = exp(i*(k*cos(theta)*x+k*sin(theta)*y));
		  ///* cube 2D */
		  
		  //if(std::real(result)!=std::real(result)) {
			//cerr << "error real\n";  
		  //}
		  
		  if(useReal) {
			  return std::real(result);
		  } else {
			  return std::imag(result);
		  }
		  
		  //return std::real((exp(i*k*x)-1-i*exp(i*k)*sin(k*x))/(pow(k,2.0))); //exact solution of 1D problem
		  //return 0;
	  }

  };
  

  
  ierr = analytical_bc1.setApproxOps(mField,"rePRES",AnaliticalFunction::fUN); CHKERRQ(ierr); //Triangles
  ierr = analytical_bc2.setApproxOps(mField,"imPRES",AnaliticalFunction::fUN); CHKERRQ(ierr);
  
  ierr = analytical_bc1.solveBcProblem(mField,"BC1_PROBLEM","BC1_FE",analytical_ditihlet_bc1); CHKERRQ(ierr);
  ierr = analytical_bc2.solveBcProblem(mField,"BC2_PROBLEM","BC2_FE",analytical_ditihlet_bc2); CHKERRQ(ierr);  
  //wait for confirmation
  
  ierr = analytical_bc1.destroyBcProblem(); CHKERRQ(ierr);
  ierr = analytical_bc2.destroyBcProblem(); CHKERRQ(ierr);
  
  //preproc
  //Preprocess the analytical Dirichlet BC
  ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc1); CHKERRQ(ierr);
  ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc2); CHKERRQ(ierr);
  
  //std::string wait;
  //ierr = MatView(A,PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
  //std::cin >> wait;
 
  ierr = mField.set_global_VecCreateGhost("ACOUSTIC_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	
  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",helmholtz_elements.getLoopFeRhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",helmholtz_elements.getLoopFeLhs()); CHKERRQ(ierr);
  //ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FLUX_FE",helmholtz_elements.getLoopFeFlux()); CHKERRQ(ierr); //scalar flux
  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FLUX_FE",helmholtz_elements.getLoopfeIncidentWave()); CHKERRQ(ierr); //Incident wave flux
  

  if(useImpedance) {
  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_IMPEDANCE_FE",helmholtz_elements.getLoopFeImpedanceLhs()); CHKERRQ(ierr);
  }
  /*above terms related to operators in HelmholtzElement.hpp */
  
  //ierr = MatView(A,PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
  //std::cin >> wait;
  //int ii1,jj1,N3;
  //ierr=MatGetSize(A,&ii1,&jj1);
  //ierr=VecGetSize(F,&N3);
  //std::cout << "\n size of stiffness matrix = \n" << ii1 << " X " << jj1 << "\n size of load vector = \n" << N3 << std::endl;
  //postproc
  //Postprocess the Analytical Dirichlet BC
  ierr = mField.problem_basic_method_postProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc1); CHKERRQ(ierr);
  ierr = mField.problem_basic_method_postProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc2); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecScale(F,-1); CHKERRQ(ierr);

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,T); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //ierr = VecView(T,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc1); CHKERRQ(ierr);
  ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc2); CHKERRQ(ierr);
  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ACOUSTIC_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);  
  
  //Wait to putput the data in format
  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Acoustic_Impining_Sphere.txt",&viewer); CHKERRQ(ierr);
  VecView(T,viewer);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  
  
  
  //if(pcomm->rank()==0) {
  rval = moab.write_file("impinging_numerical.h5m"); CHKERR_PETSC(rval);
  //}
  //destroy the KSP solvers
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&T); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  
  
  PostPocOnRefinedMesh post_proc1(mField);
  ierr = post_proc1.generateRefereneElemenMesh(); CHKERRQ(ierr);
  ierr = post_proc1.addFieldValuesPostProc("rePRES"); CHKERRQ(ierr);
  ierr = post_proc1.addFieldValuesGradientPostProc("rePRES"); CHKERRQ(ierr);
  ierr = post_proc1.addFieldValuesPostProc("imPRES"); CHKERRQ(ierr);
  ierr = post_proc1.addFieldValuesGradientPostProc("imPRES"); CHKERRQ(ierr);
  if(mField.check_field("reEX") && mField.check_field("imEX")) {
	  ierr = post_proc1.addFieldValuesPostProc("reEX"); CHKERRQ(ierr);
	  ierr = post_proc1.addFieldValuesGradientPostProc("reEX"); CHKERRQ(ierr);
	  ierr = post_proc1.addFieldValuesPostProc("imEX"); CHKERRQ(ierr);
	  ierr = post_proc1.addFieldValuesGradientPostProc("imEX"); CHKERRQ(ierr);
	  
	  ierr = post_proc1.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",post_proc1); CHKERRQ(ierr);
	  rval = post_proc1.postProcMesh.write_file("four_fields.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
	  
	  //output the results from Docker
	  char command1[] = "mbconvert ./four_fields.h5m ./four_fields.vtk && cp ./four_fields.vtk ../../../../../mnt/home/Desktop/U_pan/helmholtz_results/";
	  int todo1 = system( command1 );
	  
  } else {

  ierr = post_proc1.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",post_proc1); CHKERRQ(ierr);
  rval = post_proc1.postProcMesh.write_file("acoustic_impinging_out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
  
  //output the results from Docker
  char command1[] = "mbconvert ./acoustic_impinging_out.h5m ./acoustic_impinging_out.vtk && cp ./acoustic_impinging_out.vtk ../../../../../mnt/home/Desktop/U_pan/helmholtz_results/";
  int todo1 = system( command1 );
  }
  
  /** get the time interval **/
  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f S CPU Time = %f S \n",pcomm->rank(),v2-v1,t2-t1);
  
  Vec M;
  
  
  
  
  //wait to calculate the the magnitude of acoustic potential-abs(phi) using real and imag components.
  //ierr = VecView(T,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //std::cin >> wait;

  ////Calculate the L^2 Norm manually below
  //ierr = VecPow(C,2.0); CHKERRQ(ierr); 
  //ierr = VecAbs(C); CHKERRQ(ierr);
  //Vec H;
  //ierr = VecDuplicate(T,&H); CHKERRQ(ierr);
  //ierr = VecCopy(T,H); CHKERRQ(ierr);
  ////ierr = VecPow(H,2.0); CHKERRQ(ierr); 
  //ierr = VecAbs(H); CHKERRQ(ierr);
  ////std::cout << "\n \n \n \n \n \n \n " << std::endl;
  ////ierr = VecView(H,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ////std::cin >> wait;
  //ierr = VecAXPBY(C,-1.0,1.0,H); CHKERRQ(ierr);
  //PetscReal l2norm,pointwisenorm,sum;
  //PetscScalar pow = 2;
  ////Compute the Global L^2 norm error and point wise error
  //Vec G;
  //ierr = VecDuplicate(C,&G); CHKERRQ(ierr);
  //ierr = VecCopy(C,G); CHKERRQ(ierr);
  //ierr = VecPow(G,pow); CHKERRQ(ierr);
  //ierr = VecSum(G,&sum); CHKERRQ(ierr);
  //l2norm = sqrt(sum);
  //ierr = VecDestroy(&G); CHKERRQ(ierr);
  
  ////ierr = VecNorm(C,NORM_FROBENIUS,&l2norm);;
  ////ierr = VecNorm(C,NORM_MAX,&pointwisenorm);
  //ierr = VecMax(C,NULL,&pointwisenorm);
  ////std::cout << "\n The Global L2 Norm of error for this problem is : --\n" << l2norm << std::endl;
  //std::cout << "\n The Global Pointwise Norm of error for this problem is : --\n" << pointwisenorm << std::endl;
  //std::cin >> wait;
  ////ierr = VecSqrtAbs(C); CHKERRQ(ierr);
  //ierr = VecView(C,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //std::cin >> wait;
  //std::cout << "\n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n " << std::endl;
  																//give the right value of Im(Scalar Solution)
  

  
  
  
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


