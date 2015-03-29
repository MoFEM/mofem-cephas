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
 * abs(result) = sqrt(reEX^2+imEX^2) */

#include <MoFEM.hpp>


#include <DirichletBC.hpp>
#include <PotsProcOnRefMesh.hpp>
#include <HelmholtzElement.hpp>

#include <Projection10NodeCoordsOnField.hpp>
//#include <AnalyticalDirichletHelmholtz.hpp>
#include <AnalyticalHelmholtz.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

#include <stdexcept>
#include <cmath>
#include <boost/math/special_functions.hpp>
#include <complex>

//#define RND_EPS 1e-6

////Rounding
//double roundn(double n)
//{
//	//break n into fractional part (fract) and integral part (intp)
//	double fract, intp;
//	fract = modf(n,&intp);
//	
////    //round up
////    if (fract>=.5)
////    {
////        n*=10;
////        ceil(n);
////        n/=10;
////    }
////	//round down
////    if (fract<=.5)
////    {
////		n*=10;
////        floor(n);
////        n/=10;
////    }
//	// case where n approximates zero, set n to "positive" zero
//	if (abs(intp)==0)
//	{
//		if(abs(fract)<=RND_EPS)
//		{
//			n=0.000;
//		}
//	}
//	return n;
//}

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
	ierr = mField.add_field("reEX",H1,1); CHKERRQ(ierr);  //field order distinguish the scalar field and vector field.
	ierr = mField.add_field("imEX",H1,1); CHKERRQ(ierr);
	
	//Problem
	ierr = mField.add_problem("EX1_PROBLEM"); CHKERRQ(ierr); //analytical Dirichlet for real field
	ierr = mField.add_problem("EX2_PROBLEM"); CHKERRQ(ierr); //analytical Dirichlet for imag field

	
	ierr = mField.modify_problem_ref_level_add_bit("EX1_PROBLEM",bit_level0); CHKERRQ(ierr);  //analytical Dirichlet
	ierr = mField.modify_problem_ref_level_add_bit("EX2_PROBLEM",bit_level0); CHKERRQ(ierr);  //analytical Dirichlet

	
	
	//meshset consisting all entities in mesh
	EntityHandle root_set = moab.get_root_set(); 
	//add entities to field
	ierr = mField.add_ents_to_field_by_TETs(root_set,"reEX"); CHKERRQ(ierr);
	ierr = mField.add_ents_to_field_by_TETs(root_set,"imEX"); CHKERRQ(ierr);

	
	//set app. order
	//see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
	PetscInt order;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) {
		order = 2;
	}

	ierr = mField.set_field_order(root_set,MBTET,"reEX",order); CHKERRQ(ierr);
	ierr = mField.set_field_order(root_set,MBTRI,"reEX",order); CHKERRQ(ierr);
	ierr = mField.set_field_order(root_set,MBEDGE,"reEX",order); CHKERRQ(ierr);
	ierr = mField.set_field_order(root_set,MBVERTEX,"reEX",1); CHKERRQ(ierr);
	
	ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
	ierr = mField.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
	
	ierr = mField.set_field_order(root_set,MBTET,"imEX",order); CHKERRQ(ierr);
	ierr = mField.set_field_order(root_set,MBTRI,"imEX",order); CHKERRQ(ierr);
	ierr = mField.set_field_order(root_set,MBEDGE,"imEX",order); CHKERRQ(ierr);
	ierr = mField.set_field_order(root_set,MBVERTEX,"imEX",1); CHKERRQ(ierr);

    Range bc_tets;
	for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"EXACT_SOL",it)) {
		rval = moab.get_entities_by_type(it->get_meshset(),MBTET,bc_tets,true); CHKERR_PETSC(rval);
	}

	AnalyticalSolution analytical_bc1(mField,bc_tets);
	AnalyticalSolution analytical_bc2(mField,bc_tets);

	ierr = analytical_bc1.initializeExactProblem("EX1_PROBLEM","EX1_FE","reEX"); CHKERRQ(ierr);
	ierr = analytical_bc2.initializeExactProblem("EX2_PROBLEM","EX2_FE","imEX"); CHKERRQ(ierr);

	////add entities to finite element
	//ierr = mField.add_ents_to_finite_element_by_TETs(root_set,"EX1_FE"); CHKERRQ(ierr);
	//ierr = mField.add_ents_to_finite_element_by_TETs(root_set,"EX2_FE"); CHKERRQ(ierr);

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

	//mesh partitinoning for analytical Dirichlet
	ierr = mField.partition_problem("EX1_PROBLEM"); CHKERRQ(ierr);
	ierr = mField.partition_finite_elements("EX1_PROBLEM"); CHKERRQ(ierr);
	//what are ghost nodes, see Petsc Manual
	ierr = mField.partition_ghost_dofs("EX1_PROBLEM"); CHKERRQ(ierr);
		
	//mesh partitinoning for analytical Dirichlet
	ierr = mField.partition_problem("EX2_PROBLEM"); CHKERRQ(ierr);
	ierr = mField.partition_finite_elements("EX2_PROBLEM"); CHKERRQ(ierr);
	//what are ghost nodes, see Petsc Manual
	ierr = mField.partition_ghost_dofs("EX2_PROBLEM"); CHKERRQ(ierr);


    //Exact solution of Impinging sphere from Acoustic isogeometric boundary element analysis by R.N. Simpson etc.
    static double aNgularfreq;
	static double sPeed; //Without static. I got error:use of local variable with automatic storage from containing function
	
	for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it))
	{
		cout << endl << *it << endl;
		
		//Get block name
		string name = it->get_Cubit_name();
		if (name.compare(0,13,"MAT_HELMHOLTZ") == 0)
		{
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
			//Incident wave in Z direction.
			//double sqrtx2y2 = sqrt(pow(x,2.0)+pow(y,2.0));
			//double theta = atan2(sqrtx2y2,z)+pi;
			//double theta = acos(z/R); 
			//Incident wave in X direction.
			double theta = atan2(y,x); //the arctan of radians (y/x)
			
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
			
			
			//const complex< double > inc_field = exp( i * k * R * cos( theta ) );  //???? Incident wave
			//const complex< double > total_field = inc_field + result;
			
	
			//complex< double > result = inc_field;
			//double val = std::real(result);
			//std::string wait;
			//std::string wait;
			//std::cout << "\n theta = \n" << theta << std::endl;
			//std::cout << "\n abs( result ) = \n" << abs( result ) << std::endl;
			
			//ofs << theta << "\t" << abs( result ) << "\t" << abs( inc_field ) << "\t" << abs( total_field ) <<  "\t" << R << endl; //write the file
			
			ofs << theta << "\t" << R <<  endl;
			//if ( R == 5 ) {
			//	std::string wait;
			//std::cout << "\n wo lai le R= " << R << std::endl;
			//}
			
			
			//return std::real(result);
			double tEmp_theta = theta;
			return tEmp_theta;
			//std::string wait;
			//double rX,rY;
			//rX = roundn(x);
			//rY = y
			//std::cout << " \n roundn(x) = \n" << rX<< "\n roundn(y) = \n " << rY << 
			//		  " \n result =  \n" << rX/rY << std::endl; 

			//double result = rX/rY;

			//std::cout << "\n y/x = \n" << result << std::endl;
		    //return result; //test
	
		}
		static double fUN_imag(double x,double y,double z) {
			//const double pi = atan( 1.0 ) * 4.0;
			//double R = sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)); //radius
			////double theta = atan2(y,x)+pi; //the arctan of radians (y/x)
			//double theta = acos(z/R);

			//const double wAvenumber = aNgularfreq/sPeed;
			
			//const double k = wAvenumber;  //Wave number
			//const double a = 0.5;         //radius of the sphere,wait to modify by user
			//const double const1 = k * a;
			//double const2 = k * R;
			
			
			//const complex< double > i( 0.0, 1.0 );
			
			//// magnitude of incident wave
			//const double phi_incident_mag = 1.0;
			
			//const double tol = 1.0e-10;
			//double max = 0.0;
			//double min = 999999.0;
			
			//complex< double > result = 0.0;
			//complex< double > prev_result;
			
			//double error = 100.0;
			//unsigned int n = 0; //initialized the infinite series loop
			
			//while( error > tol )  //finding the acoustic potential in one single point.
			//{
			//	double jn_der = n / const1 * sph_bessel( n, const1 ) - sph_bessel( n + 1, const1 );  //The derivative of bessel function
			//	complex< double > hn_der = n / const1 * sph_hankel_1( n, const1 ) - sph_hankel_1( n + 1, const1 );
			//	//complex< double > hn_der = 0.5 * ( sph_hankel_1( n - 1, const1 ) -
			//	//( sph_hankel_1( n, const1 ) + const1 * sph_hankel_1( n + 1, const1 ) ) / const1 );
			//	double Pn = legendre_p( n, cos( theta ) );
			//	complex< double >hn = sph_hankel_1( n, const2 );  //S Hankel first kind function
			//	prev_result = result;
			//	result -= pow( i, n ) * ( 2.0 * n + 1.0 ) * jn_der / hn_der * Pn * hn;
			//	error = abs( abs( result ) - abs( prev_result ) );
			//	++n;
			//}
	
			
			//if (R=5) { std::cout << "\n the point is on the boundary of x-y plane ! \n" << std::endl;}
			//const complex< double > inc_field = exp( i * k * R * cos( theta ) );  //???? Incident wave
			//const complex< double > total_field = inc_field + result;
			
			//complex< double > result = inc_field;
			//return std::imag(result);
			//return (theta);
			return 5-x;
			
			
			//double result=theta;
			//if(y < 0 && x != 0) {
			//double result = atan2(y,x)+PI;
			//} else if(y > 0) {
			//double result = atan2(y,x)+PI;

			//} else if (y = 0 && x < 0) {
			//double result = atan2(y,x);
			//} else if (y = 0 && x > 0) {
			//double result = atan2(y,x)+PI;
			//}
				
			
			
		}
	};

    Mat A;
	Mat B;
	
	Vec D;
	Vec G; //Exavt solution vector = Int_V eXact dV 
	
	Vec F;
	Vec C;
	
	ierr = mField.VecCreateGhost("EX1_PROBLEM",ROW,&F); CHKERRQ(ierr);
	ierr = mField.VecCreateGhost("EX2_PROBLEM",ROW,&C); CHKERRQ(ierr);

	ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
	ierr = VecDuplicate(C,&G); CHKERRQ(ierr);

	ierr = mField.MatCreateMPIAIJWithArrays("EX1_PROBLEM",&A); CHKERRQ(ierr);
	ierr = mField.MatCreateMPIAIJWithArrays("EX2_PROBLEM",&B); CHKERRQ(ierr);

	////use the same way as Dirichlet bc
	//AnalyticalDirihletBC::ExactBC exact_bc1(mField,"reEX",B,C,G);
	//AnalyticalDirihletBC::ExactBC exact_bc2(mField,"imEX",B,C,G);
	
	
	////solve for analytical exact dofs
	//ierr = analytical_bc1.setExactProblem(mField,"EX1_PROBLEM"); CHKERRQ(ierr);
	//ierr = analytical_bc2.setExactProblem(mField,"EX2_PROBLEM"); CHKERRQ(ierr);
	
	//tets
	ierr = analytical_bc1.setExactSolRhsOp("reEX",AnaliticalFunction::fUN_real,F); CHKERRQ(ierr);
	ierr = analytical_bc1.setExactSolLhsOp("reEX",AnaliticalFunction::fUN_real,A); CHKERRQ(ierr);
	
    ierr = analytical_bc2.setExactSolRhsOp("imEX",AnaliticalFunction::fUN_imag,C); CHKERRQ(ierr);
	ierr = analytical_bc2.setExactSolLhsOp("imEX",AnaliticalFunction::fUN_imag,B); CHKERRQ(ierr);
	

	
	//ierr = analytical_bc1.solveExactProblem("EX1_PROBLEM","EX1_FE","reEX",A,F,D); CHKERRQ(ierr);

	//std::string wait;
	//std::cout << "I am fine " << std::endl;
	//std::cin >> wait;
	//ierr = analytical_bc2.solveExactProblem("EX2_PROBLEM","EX2_FE","imEX",B,C,G); CHKERRQ(ierr);
	
	ierr = VecZeroEntries(F); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecZeroEntries(D); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = MatZeroEntries(A); CHKERRQ(ierr);
    
	ierr = mField.set_global_VecCreateGhost("EX1_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

	analytical_bc1.getLoopFeRhs().snes_f = F;
	analytical_bc1.getLoopFeLhs().snes_B = A;
	
	ierr = mField.loop_finite_elements("EX1_PROBLEM","EX1_FE",analytical_bc1.getLoopFeRhs()); CHKERRQ(ierr); //wait; where problem occurs
	ierr = mField.loop_finite_elements("EX1_PROBLEM","EX1_FE",analytical_bc1.getLoopFeLhs()); CHKERRQ(ierr); //wait; where problem occurs
	

	
	ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	int ii2,jj2;
	ierr=MatGetSize(A,&ii2,&jj2);
	std::cout << "\n I am Batman @ with size \n" << ii2 << " X " << jj2 <<  std::endl;
	
	//ierr = MatView(B,PETSC_VIEWER_DRAW_WORLD);
	//std::cin >> wait;
	
	KSP solver;
	ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
	ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
	ierr = KSPSetUp(solver); CHKERRQ(ierr);
	
	ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	
	//ierr = mField.set_global_VecCreateGhost("EX1_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = mField.set_local_VecCreateGhost("EX1_PROBLEM",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	//ierr = m_field.set_other_global_VecCreateGhost(problem,re_field_name,re_field_name,ROW,G,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	//ierr = VecView(G,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	
	//mField.list_fields();
	if(mField.check_field("reEX")) { std::cout << "\n reEX is in the database . \n" << std::endl;}
	//mField.get_field_structure("reEX");
	
	PetscReal pointwisenorm;
	ierr = VecMax(D,NULL,&pointwisenorm);
	
	std::cout << "\n The Global Pointwise Norm for the solution Vector is : --\n" << pointwisenorm << std::endl;
	
	
	ierr = VecZeroEntries(C); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecZeroEntries(G); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(G,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(G,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = MatZeroEntries(B); CHKERRQ(ierr);
    
	ierr = mField.set_global_VecCreateGhost("EX2_PROBLEM",ROW,G,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

	ierr = mField.loop_finite_elements("EX2_PROBLEM","EX2_FE",analytical_bc2.getLoopFeRhs()); CHKERRQ(ierr); //wait; where problem occurs
	ierr = mField.loop_finite_elements("EX2_PROBLEM","EX2_FE",analytical_bc2.getLoopFeLhs()); CHKERRQ(ierr); //wait; where problem occurs
	
	ierr = VecGhostUpdateBegin(C,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(C,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(C); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(C); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	int ii1,jj1;
	ierr=MatGetSize(B,&ii1,&jj1);
	std::cout << "\n I am Batman @ with size \n" << ii1 << " X " << jj1 <<  std::endl;
	
	//ierr = MatView(B,PETSC_VIEWER_DRAW_WORLD);
	//std::cin >> wait;
	
	KSP solver1;
	ierr = KSPCreate(PETSC_COMM_WORLD,&solver1); CHKERRQ(ierr);
	ierr = KSPSetOperators(solver1,B,B); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(solver1); CHKERRQ(ierr);
	ierr = KSPSetUp(solver1); CHKERRQ(ierr);
	
	ierr = KSPSolve(solver1,C,G); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(G,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(G,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	
	ierr = mField.set_global_VecCreateGhost("EX2_PROBLEM",ROW,G,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	//ierr = mField.set_other_global_VecCreateGhost("EX1_PROBLEM","reEX","imEX",ROW,G,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	//ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	

	
	//PetscReal pointwisenorm2;
	//ierr = VecMax(G,NULL,&pointwisenorm2);
	
	//std::cout << "\n The Global Pointwise Norm of error for this problem is : --\n" << pointwisenorm2 << std::endl;
	
	
	//Wait to putput the data in format
	PetscViewer viewer;
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Exact_Impining_Sphere.txt",&viewer); CHKERRQ(ierr);
	VecView(D,viewer);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
	
	
	////Open mesh_file_name.txt for writing
	//ofstream myfile;
	//myfile.open("field_reEX_test.txt");
	
	//for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"reEX",dof_ptr))
    //{
    //    if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
    //    
    //    if(dof_ptr->get_dof_rank()==0)
    //    {
    //        //Round and truncate to 3 decimal places
    //        double fval = dof_ptr->get_FieldData();
    //        cout << boost::format("%.3lf") % fval << "  ";
    //        myfile << boost::format("%.3lf") % roundn(fval) << "  ";
    //    }
    //    //if(dof_ptr->get_dof_rank()==1)
    //    //{
    //    //    //Round and truncate to 3 decimal places
    //    //    double fval = dof_ptr->get_FieldData();
    //    //    cout << boost::format("%.3lf") % roundn(fval) << "  ";
    //    //    myfile << boost::format("%.3lf") % roundn(fval) << "  ";
    //    //}
    //    //if(dof_ptr->get_dof_rank()==2)
    //    //{
    //    //    //Round and truncate to 3 decimal places
    //    //    double fval = dof_ptr->get_FieldData();
    //    //    cout << boost::format("%.3lf") % roundn(fval) << endl;
    //    //    myfile << boost::format("%.3lf") % roundn(fval) << endl;
    //    //}
    //    
    //}
	//myfile.close();
	
	
	
	
	
    //tets
	
	//ierr = analytical_bc1.destroyExactProblem(); CHKERRQ(ierr);
	//ierr = analytical_bc2.destroyExactProblem(); CHKERRQ(ierr);
	
	
	
	
	
	

	if(pcomm->rank()==0) {
		rval = moab.write_file("solution.h5m"); CHKERR_PETSC(rval);
	}
	
	
	
	
	// Analytical solution	
	ProjectionFieldOn10NodeTet ent_method_on_10nodeTet1(mField,"reEX",true,false,"reEX");
	ProjectionFieldOn10NodeTet ent_method_on_10nodeTet2(mField,"imEX",true,false,"imEX");
	ierr = mField.loop_dofs("reEX",ent_method_on_10nodeTet1); CHKERRQ(ierr);
	ierr = mField.loop_dofs("imEX",ent_method_on_10nodeTet2); CHKERRQ(ierr);
	ent_method_on_10nodeTet1.set_nodes = false;
	ent_method_on_10nodeTet2.set_nodes = false;
	ierr = mField.loop_dofs("reEX",ent_method_on_10nodeTet1); CHKERRQ(ierr);
	ierr = mField.loop_dofs("imEX",ent_method_on_10nodeTet2); CHKERRQ(ierr);
	
	
	
	if(pcomm->rank()==0) {
		EntityHandle out_meshset1;
		rval = moab.create_meshset(MESHSET_SET,out_meshset1); CHKERR_PETSC(rval);
		ierr = mField.problem_get_FE("EX1_PROBLEM","EX1_FE",out_meshset1); CHKERRQ(ierr);
		//ierr = mField.problem_get_FE("EX2_PROBLEM","EX2_FE",out_meshset1); CHKERRQ(ierr);
		rval = moab.write_file("Exact_solution1.vtk","VTK","",&out_meshset1,1); CHKERR_PETSC(rval);
		rval = moab.delete_entities(&out_meshset1,1); CHKERR_PETSC(rval);
	}
	
	//if(pcomm->rank()==0) {
	//	EntityHandle out_meshset2;
	//	rval = moab.create_meshset(MESHSET_SET,out_meshset2); CHKERR_PETSC(rval);
	//	ierr = mField.problem_get_FE("EX2_PROBLEM","EX2_FE",out_meshset2); CHKERRQ(ierr);
	//	rval = moab.write_file("Exact_solution2.vtk","VTK","",&out_meshset2,1); CHKERR_PETSC(rval);
	//	rval = moab.delete_entities(&out_meshset2,1); CHKERR_PETSC(rval);
	//}
	
	
	
	char command2[] = "cp Exact_solution1.vtk ../../../../../mnt/home/Documents/mofem-cephas/mofem_v0.2/users_modules/helmholtz/";
	
	int status = system( command2 ); 
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	ierr = MatDestroy(&A); CHKERRQ(ierr);
	ierr = VecDestroy(&F); CHKERRQ(ierr);
	ierr = VecDestroy(&D); CHKERRQ(ierr);
	
	ierr = MatDestroy(&B); CHKERRQ(ierr);
	ierr = VecDestroy(&C); CHKERRQ(ierr);
	ierr = VecDestroy(&G); CHKERRQ(ierr);


	ierr = KSPDestroy(&solver); CHKERRQ(ierr);
	ierr = KSPDestroy(&solver1); CHKERRQ(ierr);



	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;

}
