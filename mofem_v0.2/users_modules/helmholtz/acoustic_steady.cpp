/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
 */

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *PhD student Thomas Xuan Meng
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
//#include <L2normElement.hpp>

#include <Projection10NodeCoordsOnField.hpp>
#include <AnalyticalDirichletHelmholtz.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
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

  //Create MoFEM (Joseph) database
  //FieldCore core(moab);
  MoFEM::Core core(moab);
  FieldInterface& mField = core;

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

  ierr = mField.add_problem("EX1_PROBLEM"); CHKERRQ(ierr); //analytical Dirichlet for real field
  ierr = mField.add_problem("EX2_PROBLEM"); CHKERRQ(ierr); //analytical Dirichlet for imag field

  
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ACOUSTIC_PROBLEM",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("BC1_PROBLEM",bit_level0); CHKERRQ(ierr);  //analytical Dirichlet
  ierr = mField.modify_problem_ref_level_add_bit("BC2_PROBLEM",bit_level0); CHKERRQ(ierr);  //analytical Dirichlet
  
  ierr = mField.modify_problem_ref_level_add_bit("EX1_PROBLEM",bit_level0); CHKERRQ(ierr);  //analytical Dirichlet
  ierr = mField.modify_problem_ref_level_add_bit("EX2_PROBLEM",bit_level0); CHKERRQ(ierr);  //analytical Dirichlet

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

  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(root_set,MBTET,"imPRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"imPRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"imPRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"imPRES",1); CHKERRQ(ierr);

  HelmholtzElement helmholtz_elements(mField);               //Create the HelmholtzElement class in the header-file
  //L2normElement l2norm_element(mField);   //l2norm
  
  ierr = helmholtz_elements.addHelmholtzElements("ACOUSTIC_PROBLEM","rePRES","imPRES"); CHKERRQ(ierr);
  ierr = helmholtz_elements.addHelmholtzFluxElement("ACOUSTIC_PROBLEM","rePRES","imPRES"); CHKERRQ(ierr);
  ierr = helmholtz_elements.addHelmholtzImpedanceElement("ACOUSTIC_PROBLEM","rePRES","imPRES"); CHKERRQ(ierr);
  
  //ierr = l2norm_element.addL2NormElements("ACOUSTIC_PROBLEM","L2_NORM","rePRES","imPRES"); CHKERRQ(ierr); //l2norm
  
  //Set up the analytical Dirichlet BC.And Exact Solution
  Range bc_tris;
  Range bc_tets;
  for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"ANALYTICAL_BC",it)) {
   rval = moab.get_entities_by_type(it->get_meshset(),MBTRI,bc_tris,true); CHKERR_PETSC(rval);
  }
  
  //for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"EXACT_SOL",it)) {
//	  rval = moab.get_entities_by_type(it->get_meshset(),MBTET,bc_tets,true); CHKERR_PETSC(rval);
  //}
  
  AnalyticalDirihletBC analytical_bc1(mField,bc_tris,bc_tets);
  AnalyticalDirihletBC analytical_bc2(mField,bc_tris,bc_tets);
  AnalyticalDirihletBC analytical_bc3(mField,bc_tris,bc_tets);
  AnalyticalDirihletBC analytical_bc4(mField,bc_tris,bc_tets);
  ierr = analytical_bc1.initializeBcProblem(mField,"BC1_PROBLEM","BC1_FE","rePRES"); CHKERRQ(ierr);
  ierr = analytical_bc2.initializeBcProblem(mField,"BC2_PROBLEM","BC2_FE","imPRES"); CHKERRQ(ierr);
  
  //ierr = analytical_bc3.initializeExactProblem(mField,"EX1_PROBLEM","EX1_FE","rePRES"); CHKERRQ(ierr);
  //ierr = analytical_bc4.initializeExactProblem(mField,"EX2_PROBLEM","EX2_FE","imPRES"); CHKERRQ(ierr);
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
  
  //mesh partitinoning for analytical Dirichlet
  ierr = mField.simple_partition_problem("EX1_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("EX1_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.simple_partition_problem("EX2_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("EX2_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("EX1_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("EX2_PROBLEM"); CHKERRQ(ierr);

  
  Vec F;  //Right hand side vector
  ierr = mField.VecCreateGhost("ACOUSTIC_PROBLEM",ROW,&F); CHKERRQ(ierr);
  Vec T; //Solution vector
  ierr = VecDuplicate(F,&T); CHKERRQ(ierr);
  Mat A; //Left hand side matrix
  ierr = mField.MatCreateMPIAIJWithArrays("ACOUSTIC_PROBLEM",&A); CHKERRQ(ierr);
  
  //Vec D; //l2norm vector
  //ierr = VecDuplicate(F,&D); CHKERRQ(ierr); //l2norm vector

  //int ii1,jj1;
  //ierr=MatGetSize(A,&ii1,&jj1);
  std::string wait;
  //MatView(A,PETSC_VIEWER_DRAW_WORLD);
  //std::cin >> wait;
  //PetscInt i,j,N1,N2,N3,N4;


  
  bool useScalar = true;
  //ierr = helmholtz_elements.setHelmholtzFiniteElementRhsOperators_F("rePRES","rePRES",F,useScalar); CHKERRQ(ierr); //The analytical F source vector
  ierr = helmholtz_elements.setHelmholtzFiniteElementRhsOperators_rere("rePRES","rePRES",F); CHKERRQ(ierr); 
  ierr = helmholtz_elements.setHelmholtzFiniteElementRhsOperators_imim("imPRES","imPRES",F); CHKERRQ(ierr); 
  ierr = helmholtz_elements.setHelmholtzFiniteElementLhsOperators("rePRES","imPRES",(A)); CHKERRQ(ierr);
  ierr = helmholtz_elements.setHelmholtzMassFiniteElementLhsOperators("rePRES","imPRES",(A)); CHKERRQ(ierr);
  //ierr = helmholtz_elements.setHelmholtzFluxFiniteElementRhsOperators("rePRES","rePRES",F); CHKERRQ(ierr);  //real Neumann BC
  //ierr = helmholtz_elements.setHelmholtzFluxFiniteElementRhsOperators("imPRES","imPRES",F); CHKERRQ(ierr);    //Imag Neumann BC
  ierr = helmholtz_elements.setHelmholtzIncidentWaveFiniteElementRhsOperators("rePRES","imPRES",F); CHKERRQ(ierr); // Incident wave flux.

  ierr = helmholtz_elements.setHelmholtzImpedanceFiniteElementLhsOperators("rePRES","imPRES",(A)); CHKERRQ(ierr);
  

  
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

  //MoFEM::HelmholtzElement::BlockData blockData;
  //double wavenumber = helmholtz_elements.getwavenumber();
  //// = helmholtz_elements.getblockdata();
  ////MoFEM::HelmholtzElement::getblockdata();  //wait why cant use this to get the attribute? 
  
  static double aNgularfreq;
  static double sPeed; //Without static. I got error:use of local variable with automatic storage from containing function

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

  
  struct AnaliticalFunction {
	  static double fUN_real(double x,double y,double z) {
		  ////return pow(x,1);
		  const double pi = atan( 1.0 ) * 4.0;
		  double R = sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)); //radius
		  double theta = atan2(y,x)+pi; //the arctan of radians (y/x)
		  
		  //if(theta < 0) {theta += 2 * pi;}  //if the return radians is in the lower half plane, add 2 pi to the results.
		  const double wAvenumber = aNgularfreq/sPeed;
		  
		  const double k = wAvenumber;  //Wave number
		  const double a = 0.5;         //radius of the sphere,wait to modify by user
		  const double const1 = k * a;
		  double const2 = k * R;
		  
		  
		  const complex< double > i( 0.0, 1.0 );
		  
		  // magnitude of incident wave
		  const double phi_incident_mag = 1.0;
		  
		  const double tol = 1.0e-10;
		  double max = 0.0;
		  double min = 999999.0;
		  
		  complex< double > result = 0.0;
		  complex< double > prev_result;
		  
		  double error = 100.0;
		  unsigned int n = 0; //initialized the infinite series loop
		  
		  while( error > tol )  //finding the acoustic potential in one single point.
		  {
			  double jn_der = n / const1 * sph_bessel( n, const1 ) - sph_bessel( n + 1, const1 );  //The derivative of bessel function
			  complex< double > hn_der = n / const1 * sph_hankel_1( n, const1 ) - sph_hankel_1( n + 1, const1 );
			  //complex< double > hn_der = 0.5 * ( sph_hankel_1( n - 1, const1 ) -
			  //( sph_hankel_1( n, const1 ) + const1 * sph_hankel_1( n + 1, const1 ) ) / const1 );
			  double Pn = legendre_p( n, cos( theta ) );
			  complex< double >hn = sph_hankel_1( n, const2 );  //S Hankel first kind function
			  prev_result = result;
			  result -= pow( i, n ) * ( 2.0 * n + 1.0 ) * jn_der / hn_der * Pn * hn;
			  error = abs( abs( result ) - abs( prev_result ) );
			  ++n;
		  }
          
		  const complex< double > inc_field = exp( i * k * R * cos( theta ) );  //???? Incident wave
		  const complex< double > total_field = inc_field + result;
		  
		  if ( R == 5 ) {
			  //std::string wait;
			  //std::cout << "\n wo lai le R= " << R << std::endl;
		  }
		  
		  double val = std::real(result);
		  //std::string wait;
		  //std::string wait;
		  //std::cout << "\n theta = \n" << theta << std::endl;
		  //std::cout << "\n abs( result ) = \n" << abs( result ) << std::endl;
		  ofs << theta << "\t" << abs( result ) << "\t" << abs( inc_field ) << "\t" << abs( total_field ) << "\t" << R << endl; //write the file
		  
          return std::real(result);

          //return sin(x);  //test

	  }
	  static double fUN_imag(double x,double y,double z) {
		  const double pi = atan( 1.0 ) * 4.0;
		  double R = sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)); //radius
		  double theta = atan2(y,x)+pi; //the arctan of radians (y/x)
		  //if(theta < 0) {theta += 2 * pi;}
		  //(x > 0 ? x : (2*PI + x)) * 360 / (2*PI)
		  //double theta = atan(y/x);
		  const double wAvenumber = aNgularfreq/sPeed;
		  
		  const double k = wAvenumber;  //Wave number
		  const double a = 0.5;         //radius of the sphere,wait to modify by user
		  const double const1 = k * a;
		  double const2 = k * R;
		  
		  
		  const complex< double > i( 0.0, 1.0 );
		 
		  // magnitude of incident wave
		  const double phi_incident_mag = 1.0;
		  
		  const double tol = 1.0e-10;
		  double max = 0.0;
		  double min = 999999.0;
		  
		  complex< double > result = 0.0;
		  complex< double > prev_result;
		  
		  double error = 100.0;
		  unsigned int n = 0; //initialized the infinite series loop
		  
		  while( error > tol )  //finding the acoustic potential in one single point.
		  {
			  double jn_der = n / const1 * sph_bessel( n, const1 ) - sph_bessel( n + 1, const1 );  //The derivative of bessel function
			  complex< double > hn_der = n / const1 * sph_hankel_1( n, const1 ) - sph_hankel_1( n + 1, const1 );
			  //complex< double > hn_der = 0.5 * ( sph_hankel_1( n - 1, const1 ) -
			  //( sph_hankel_1( n, const1 ) + const1 * sph_hankel_1( n + 1, const1 ) ) / const1 );
			  double Pn = legendre_p( n, cos( theta ) );
			  complex< double >hn = sph_hankel_1( n, const2 );  //S Hankel first kind function
			  prev_result = result;
			  result -= pow( i, n ) * ( 2.0 * n + 1.0 ) * jn_der / hn_der * Pn * hn;
			  error = abs( abs( result ) - abs( prev_result ) );
			  ++n;
		  }

          return std::imag(result);
		  
		    //return cos(x); //test
	  
	  }
  };
   
  //ierr = l2norm_element.setL2NormRelativeError(mField,"rePRES","rePRES",D,AnaliticalFunction::fUN_real); CHKERRQ(ierr);//l2norm
  //ierr = l2norm_element.setL2NormRelativeError(mField,"imPRES","imPRES",D,AnaliticalFunction::fUN_imag); CHKERRQ(ierr);//l2norm
  //ierr = VecZeroEntries(D); CHKERRQ(ierr);//l2norm
  //ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);//l2norm
  
  //solve the analytical Dirichlet BC, results in a solution vector contains the the solution on Dirichlet nodes only. 
  
  //Create the Matrix and vector to calculate the exact solution.
  
  //Mat B;
  //MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&B);
  
  //Vec G; //Exavt solution vector = Int_V eXact dV 
  //Vec C;
  //ierr = VecDuplicate(T,&G); CHKERRQ(ierr); //Exact solution Load Vector
  //ierr = VecDuplicate(T,&C); CHKERRQ(ierr); //Exact solution vector
  

  
  ierr = analytical_bc1.setApproxOps(mField,"rePRES",AnaliticalFunction::fUN_real); CHKERRQ(ierr); //Triangles
  ierr = analytical_bc2.setApproxOps(mField,"imPRES",AnaliticalFunction::fUN_imag); CHKERRQ(ierr);
  
  ierr = analytical_bc1.solveBcProblem(mField,"BC1_PROBLEM","BC1_FE",analytical_ditihlet_bc1); CHKERRQ(ierr);
  ierr = analytical_bc2.solveBcProblem(mField,"BC2_PROBLEM","BC2_FE",analytical_ditihlet_bc2); CHKERRQ(ierr);
  
  ////use the same way as Dirichlet bc
  //AnalyticalDirihletBC::ExactBC exact_bc1(mField,"rePRES",B,C,G);
  //AnalyticalDirihletBC::ExactBC exact_bc2(mField,"imPRES",B,C,G);
  
  ////solve for analytical dirichlet exact dofs
  //ierr = analytical_bc3.setExactProblem(mField,"EX1_PROBLEM"); CHKERRQ(ierr);
  //ierr = analytical_bc4.setExactProblem(mField,"EX2_PROBLEM"); CHKERRQ(ierr);
  
  ////tets
  //ierr = analytical_bc3.setExactSolOp(mField,"rePRES",AnaliticalFunction::fUN_real); CHKERRQ(ierr);
  //ierr = analytical_bc4.setExactSolOp(mField,"imPRES",AnaliticalFunction::fUN_imag); CHKERRQ(ierr);
  

  
 // ierr = analytical_bc3.solveExactProblem(mField,"EX1_PROBLEM","EX1_FE"); CHKERRQ(ierr);
 // std::cout << "\n upto here I am fine. Not lying!!!! \n" << std::endl;
 // std::cin >> wait;
 // ierr = analytical_bc4.solveExactProblem(mField,"EX2_PROBLEM","EX2_FE"); CHKERRQ(ierr);
 ////tets

  



  
  

  
  
  //wait for confirmation
  
  ierr = analytical_bc1.destroyBcProblem(); CHKERRQ(ierr);
  ierr = analytical_bc2.destroyBcProblem(); CHKERRQ(ierr);
  //ierr = analytical_bc3.destroyExactProblem(); CHKERRQ(ierr);
  //ierr = analytical_bc4.destroyExactProblem(); CHKERRQ(ierr);

  
  //preproc
  //Preprocess the analytical Dirichlet BC
  ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc1); CHKERRQ(ierr);
  ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc2); CHKERRQ(ierr);
  

  
  
  ////save data on mesh
  //ierr = mField.set_global_VecCreateGhost("ACOUSTIC_PROBLEM",ROW,C,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  
  ////L^2norm 
  
  //ProjectionFieldOn10NodeTet ent_method_on_10nodeTet3(mField,"rePRES",true,false,"rePRES");
  //ProjectionFieldOn10NodeTet ent_method_on_10nodeTet4(mField,"imPRES",true,false,"imPRES");
  //ierr = mField.loop_dofs("rePRES",ent_method_on_10nodeTet3); CHKERRQ(ierr);
  //ierr = mField.loop_dofs("imPRES",ent_method_on_10nodeTet4); CHKERRQ(ierr);
  //ent_method_on_10nodeTet3.set_nodes = false;
  //ent_method_on_10nodeTet4.set_nodes = false;
  //ierr = mField.loop_dofs("rePRES",ent_method_on_10nodeTet3); CHKERRQ(ierr);
  //ierr = mField.loop_dofs("imPRES",ent_method_on_10nodeTet4); CHKERRQ(ierr);
  
  
  
  //if(pcomm->rank()==0) {
//	  EntityHandle out_meshset2;
//	  rval = moab.create_meshset(MESHSET_SET,out_meshset2); CHKERR_PETSC(rval);
//	  ierr = mField.problem_get_FE("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",out_meshset2); CHKERRQ(ierr);
//	  rval = moab.write_file("L2norm.vtk","VTK","",&out_meshset2,1); CHKERR_PETSC(rval);
//	  rval = moab.delete_entities(&out_meshset2,1); CHKERR_PETSC(rval);
  //}
  
  
  
  //char command2[] = "cp L2norm.vtk ../../../../../mnt/home/Documents/mofem-cephas/mofem_v0.2/users_modules/helmholtz/";
  
  //int status = system( command2 ); 
  


  //clear the Dirichlet BC with volume indices.
  
  ierr = mField.set_global_VecCreateGhost("ACOUSTIC_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	
  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",helmholtz_elements.getLoopFeRhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",helmholtz_elements.getLoopFeLhs()); CHKERRQ(ierr);
  //ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FLUX_FE",helmholtz_elements.getLoopFeFlux()); CHKERRQ(ierr); //scalar flux
  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FLUX_FE",helmholtz_elements.getLoopfeIncidentWave()); CHKERRQ(ierr); //Incident wave flux
  
  /*Wait for confirmation */
  //ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_IMPEDANCE_FE",helmholtz_elements.getLoopFeImpedanceRhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_IMPEDANCE_FE",helmholtz_elements.getLoopFeImpedanceLhs()); CHKERRQ(ierr);
  /*above terms related to operators in HelmholtzElement.hpp */
  
  //ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","L2_NORM",l2norm_element.getLoopFeRhs()); CHKERRQ(ierr);//l2norm
  
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
  
  //MatView(A,PETSC_VIEWER_DRAW_WORLD);
  //std::cin >> wait;
  
  //ierr=VecGetSize(F,&N3);
  //ierr=VecGetSize(T,&N4);
  //std::cout<< "\n Solution Vector T = \n" << N3 << std::endl;
  //std::cin >> wait;
  
  //std::cout<< "\n RHS Vector F = \n" << N4 << std::endl;
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //std::cin >> wait;
  
  //int N1,sizeT;
  //ierr = VecGetSize(F,&N1); CHKERRQ(ierr);
  //ierr = VecGetSize(T,&sizeT); CHKERRQ(ierr);
  //std::cout << "\n N1 = \n" << N1 << std::endl;
  //std::cout << "\n sizeT = \n" << sizeT << std::endl;

  ierr = MatView(A,PETSC_VIEWER_DRAW_WORLD);
  std::cin >> wait;
  
  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,T); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  //there the vector T^tude and F^tude is without the Dirichlet BC value inserted.
  //ierr = VecView(T,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //below fill in the solution vector T with analytical dirichlet solutions.

  //ierr = VecView(T,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc1); CHKERRQ(ierr);
  ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc2); CHKERRQ(ierr);
  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ACOUSTIC_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);


  PetscReal pointwisenorm;
  ierr = VecMax(T,NULL,&pointwisenorm);
  
  std::cout << "\n The Global Pointwise Norm of error for this problem is : --\n" << pointwisenorm << std::endl;
  
  //ierr = VecView(T,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  
  //ierr = VecView(T,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  //Wait to putput the data in format
  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Acoustic_Impining_Sphere.txt",&viewer); CHKERRQ(ierr);
  VecView(T,viewer);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  
  if(pcomm->rank()==0) {
	  rval = moab.write_file("solution.h5m"); CHKERR_PETSC(rval);
  }

  //Open mesh_file_name.txt for writing
  ofstream myfile;
  myfile.open("field_rePRES_test.txt");
  
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"rePRES",dof_ptr))
  {
	  if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
	  
	  if(dof_ptr->get_dof_rank()==0)
	  {
		  //Round and truncate to 3 decimal places
		  double fval = dof_ptr->get_FieldData();
		  cout << boost::format("%.3lf") % fval << "  ";
		  myfile << boost::format("%.3lf") % fval << "  ";
	  }
	  //if(dof_ptr->get_dof_rank()==1)
	  //{
	  //    //Round and truncate to 3 decimal places
	  //    double fval = dof_ptr->get_FieldData();
	  //    cout << boost::format("%.3lf") % roundn(fval) << "  ";
	  //    myfile << boost::format("%.3lf") % roundn(fval) << "  ";
	  //}
	  //if(dof_ptr->get_dof_rank()==2)
	  //{
	  //    //Round and truncate to 3 decimal places
	  //    double fval = dof_ptr->get_FieldData();
	  //    cout << boost::format("%.3lf") % roundn(fval) << endl;
	  //    myfile << boost::format("%.3lf") % roundn(fval) << endl;
	  //}
	  
  }
  myfile.close();
  
  
  
  
  PostPocOnRefinedMesh post_proc1(mField);
  ierr = post_proc1.generateRefereneElemenMesh(); CHKERRQ(ierr);
  ierr = post_proc1.addFieldValuesPostProc("reEX"); CHKERRQ(ierr);
  ierr = post_proc1.addFieldValuesGradientPostProc("reEX"); CHKERRQ(ierr);
  ierr = post_proc1.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("EX1_PROBLEM","FE1",post_proc1); CHKERRQ(ierr);
  rval = post_proc1.postProcMesh.write_file("acoustic_out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
  
  PostPocOnRefinedMesh post_proc2(mField);
  ierr = post_proc2.generateRefereneElemenMesh(); CHKERRQ(ierr);
  ierr = post_proc2.addFieldValuesPostProc("imEX"); CHKERRQ(ierr);
  ierr = post_proc2.addFieldValuesGradientPostProc("imEX"); CHKERRQ(ierr);
  ierr = post_proc2.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("EX2_PROBLEM","FE2",post_proc2); CHKERRQ(ierr);
  rval = post_proc2.postProcMesh.write_file("acoustic_out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
  
  char command1[] = "mbconvert ./acoustic_out.h5m ./acoustic_out.vtk && cp ./acoustic_out.vtk ../../../../../mnt/home/Documents/mofem-cephas/mofem_v0.2/users_modules/helmholtz/";
  //char command2[] = "cp ./new_out.vtk ../../../../../mnt/home/Documents/mofem-cephas/mofem_v0.2/users_modules/helmholtz/";
  
  
  int todo1 = system( command1 );
  //int todo2 = system( command2 ); 
  
  
  
  
  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet1(mField,"rePRES",true,false,"rePRES");
  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet2(mField,"imPRES",true,false,"imPRES");
  ierr = mField.loop_dofs("rePRES",ent_method_on_10nodeTet1); CHKERRQ(ierr);
  ierr = mField.loop_dofs("imPRES",ent_method_on_10nodeTet2); CHKERRQ(ierr);
  ent_method_on_10nodeTet1.set_nodes = false;
  ent_method_on_10nodeTet2.set_nodes = false;
  ierr = mField.loop_dofs("rePRES",ent_method_on_10nodeTet1); CHKERRQ(ierr);
  ierr = mField.loop_dofs("imPRES",ent_method_on_10nodeTet2); CHKERRQ(ierr);
  
  
  
  if(pcomm->rank()==0) {
	  EntityHandle out_meshset1;
	  rval = moab.create_meshset(MESHSET_SET,out_meshset1); CHKERR_PETSC(rval);
	  ierr = mField.problem_get_FE("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",out_meshset1); CHKERRQ(ierr);
	  rval = moab.write_file("out_analytical.vtk","VTK","",&out_meshset1,1); CHKERR_PETSC(rval);
	  rval = moab.delete_entities(&out_meshset1,1); CHKERR_PETSC(rval);
  }
  
  //Copy the output .vtk file to desired location.
  char command[] = "cp out_analytical.vtk ../../../../../mnt/home/Documents/mofem-cephas/mofem_v0.2/users_modules/helmholtz/";
  
  int status = system( command );
  
  
  
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
  

  

  
  
  //ierr = VecView(C,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //std::cin >> wait;
  

  
  //L^2norm 
  
  //Range ref_edges;
  //ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBEDGE,ref_edges); CHKERRQ(ierr);
  //rval = moab.list_entities(ref_edges); CHKERR_PETSC(rval);
  //mField.list_dofs_by_field_name("TEMP");


  /*EntityHandle fe_meshset = mField.get_finite_element_meshset("HELMHOLTZ_FE");
  Range tets;
  rval = moab.get_entities_by_type(fe_meshset,MBTET,tets,true); CHKERR_PETSC(rval);
  Range tets_edges;
  rval = moab.get_adjacencies(tets,1,false,tets_edges,Interface::UNION); CHKERR(rval);
  EntityHandle edges_meshset;
  rval = moab.create_meshset(MESHSET_SET,edges_meshset); CHKERR_PETSC(rval);
  rval = moab.add_entities(edges_meshset,tets); CHKERR_PETSC(rval);
  rval = moab.add_entities(edges_meshset,tets_edges); CHKERR_PETSC(rval);
  rval = moab.convert_entities(edges_meshset,true,false,false); CHKERR_PETSC(rval);*/


  
  //if(pcomm->rank()==0) {
  //  PostProcScalarFieldsAndGradientOnRefMesh fe_post_proc_method1(moab,"rePRES");
//	PostProcScalarFieldsAndGradientOnRefMesh fe_post_proc_method2(moab,"imPRES");
  //  fe_post_proc_method1.do_broadcast = false;
//	fe_post_proc_method2.do_broadcast = false;
  //  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",fe_post_proc_method1,0,pcomm->size());  CHKERRQ(ierr);
//	ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",fe_post_proc_method2,0,pcomm->size());  CHKERRQ(ierr);
  //  rval = fe_post_proc_method1.moab_post_proc.write_file("out_post_proc_re.vtk","VTK",""); CHKERR_PETSC(rval);
//	rval = fe_post_proc_method2.moab_post_proc.write_file("out_post_proc_im.vtk","VTK",""); CHKERR_PETSC(rval); //wait, it seems does not give
  //}                 																					//give the right value of Im(Scalar Solution)
  
	  
  //destroy the KSP solvers
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&T); CHKERRQ(ierr);
  //ierr = VecDestroy(&H); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  //wait, destroy the exact solution solvers
  //ierr = MatDestroy(&B); CHKERRQ(ierr);
  //ierr = VecDestroy(&C); CHKERRQ(ierr);
  //ierr = VecDestroy(&G); CHKERRQ(ierr);`
  //wait for confirmation.
  
  
  
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


