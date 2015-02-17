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
#include <Projection10NodeCoordsOnField.hpp>
#include <HelmholtzElement.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <FiledApproximation.hpp>
#include <PotsProcOnRefMesh.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

#include <stdexcept>
#include <cmath>
#include <boost/math/special_functions.hpp>
#include <complex>

#define HOON

using namespace std;
using namespace boost::math;
namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;
using namespace MoFEM;

static char help[] = "...\n\n";

struct MyFunApprox_re {
	
	 ublas::vector<double> result1;
	 //double wAvenumber;
	 //ublas::vector<double>& operator()(double x, double y, double z) {
	//	result.resize(3);
	//	result[0] = x;
	//	result[1] = y;
	//	result[2] = z*z;
	//	return result;
	//}     
	 
	 //MyFunApprox_re(double wavenumber):
		 //wAvenumber(wavenumber) {}
	 //~MyFunApprox_re() {}
	 
	
	ublas::vector<double>& operator()(double x, double y, double z) {
		////return pow(x,1);
		const double pi = atan( 1.0 ) * 4.0;
		double R = sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)); //radius
		//Incident wave in Z direction.
		//double sqrtx2y2 = sqrt(pow(x,2.0)+pow(y,2.0));
		//double theta = atan2(sqrtx2y2,z)+pi;
		//double theta = acos(z/R); 
		//Incident wave in X direction.
		double theta = atan2(y,x)+2*pi; //the arctan of radians (y/x)
		
		//if(theta < 0) {theta += 2 * pi;}  //if the return radians is in the lower half plane, add 2 pi to the results.
		//const double wAvenumber = aNgularfreq/sPeed;
		double wAvenumber = 2;
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
		//ofs << theta << "\t" << abs( result ) << "\t" << abs( inc_field ) << "\t" << abs( total_field ) <<  "\t" << R << endl; //write the file
		
		//if ( R == 5 ) {
		//	std::string wait;
		//std::cout << "\n wo lai le R= " << R << std::endl;
		//}
		
		result1.resize(1);
		result1[0] = std::real(result);
		//result1[0] = x;
		
		return result1;
		
		
	}
	
	
};

struct MyFunApprox_im {
	
	ublas::vector<double> result1;
	//double wAvenumber;
	//ublas::vector<double>& operator()(double x, double y, double z) {
	//	result.resize(3);
	//	result[0] = x;
	//	result[1] = y;
	//	result[2] = z*z;
	//	return result;
	//}     
	
	//MyFunApprox_re(double wavenumber):
	//wAvenumber(wavenumber) {}
	//~MyFunApprox_re() {}
	
	
	ublas::vector<double>& operator()(double x, double y, double z) {
		////return pow(x,1);
		const double pi = atan( 1.0 ) * 4.0;
		double R = sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)); //radius
		//Incident wave in Z direction.
		//double sqrtx2y2 = sqrt(pow(x,2.0)+pow(y,2.0));
		//double theta = atan2(sqrtx2y2,z)+pi;
		//double theta = acos(z/R); 
		//Incident wave in X direction.
		double theta = atan2(y,x)+2*pi; //the arctan of radians (y/x)
		
		//if(theta < 0) {theta += 2 * pi;}  //if the return radians is in the lower half plane, add 2 pi to the results.
		//const double wAvenumber = aNgularfreq/sPeed;
		double wAvenumber = 2;
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
		//ofs << theta << "\t" << abs( result ) << "\t" << abs( inc_field ) << "\t" << abs( total_field ) <<  "\t" << R << endl; //write the file
		
		//if ( R == 5 ) {
		//	std::string wait;
		//std::cout << "\n wo lai le R= " << R << std::endl;
		//}
		
		result1.resize(1);
		result1[0] = std::imag(result);
		//result1[0] = x;
		
		return result1;
		
		
	}
	
	
};



//argc = argument counts, argv = argument vectors
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
		SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
	}
	
	//Create MoFEM (Joseph) database
	MoFEM::Core core(moab);
	FieldInterface& m_field = core;
	
	ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
	if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
	
	const char *option;
	option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
	BARRIER_RANK_START(pcomm) 
	rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
	BARRIER_RANK_END(pcomm) 
	
	//set entitities bit level
	BitRefLevel bit_level0;
	bit_level0.set(0);
	EntityHandle meshset_level0;
	rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
	ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
	
	//Fields
	ierr = m_field.add_field("reEX",H1,1); CHKERRQ(ierr);
	ierr = m_field.add_field("imEX",H1,1); CHKERRQ(ierr);
	#ifdef HOON
	ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
	#endif
	
	//FE
	ierr = m_field.add_finite_element("FE1"); CHKERRQ(ierr);
	ierr = m_field.add_finite_element("FE2"); CHKERRQ(ierr);
	
	//Define rows/cols and element data
	ierr = m_field.modify_finite_element_add_field_row("FE1","reEX"); CHKERRQ(ierr);
	ierr = m_field.modify_finite_element_add_field_col("FE1","reEX"); CHKERRQ(ierr);
	ierr = m_field.modify_finite_element_add_field_data("FE1","reEX"); CHKERRQ(ierr);
	#ifdef HOON
	ierr = m_field.modify_finite_element_add_field_data("FE1","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	#endif
	
	//Define rows/cols and element data
	ierr = m_field.modify_finite_element_add_field_row("FE2","imEX"); CHKERRQ(ierr);
	ierr = m_field.modify_finite_element_add_field_col("FE2","imEX"); CHKERRQ(ierr);
	ierr = m_field.modify_finite_element_add_field_data("FE2","imEX"); CHKERRQ(ierr);
	#ifdef HOON
	ierr = m_field.modify_finite_element_add_field_data("FE2","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	#endif
	
	//Problem
	ierr = m_field.add_problem("EX1_PROBLEM"); CHKERRQ(ierr);
	ierr = m_field.add_problem("EX2_PROBLEM"); CHKERRQ(ierr);
	
	//set finite elements for problem
	ierr = m_field.modify_problem_add_finite_element("EX1_PROBLEM","FE1"); CHKERRQ(ierr);
	ierr = m_field.modify_problem_add_finite_element("EX2_PROBLEM","FE2"); CHKERRQ(ierr);
	//set refinment level for problem
	ierr = m_field.modify_problem_ref_level_add_bit("EX1_PROBLEM",bit_level0); CHKERRQ(ierr);
	ierr = m_field.modify_problem_ref_level_add_bit("EX2_PROBLEM",bit_level0); CHKERRQ(ierr);
	
	//meshset consisting all entities in mesh
	EntityHandle root_set = moab.get_root_set(); 
	//add entities to field
	ierr = m_field.add_ents_to_field_by_TETs(root_set,"reEX"); CHKERRQ(ierr);
	ierr = m_field.add_ents_to_field_by_TETs(root_set,"imEX"); CHKERRQ(ierr);
	#ifdef HOON
    ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	#endif
	//add entities to finite element
	ierr = m_field.add_ents_to_finite_element_by_TETs(root_set,"FE1"); CHKERRQ(ierr);
	ierr = m_field.add_ents_to_finite_element_by_TETs(root_set,"FE2"); CHKERRQ(ierr);
	
	
	//set app. order
	//see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
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
	#ifdef HOON
	ierr = m_field.set_field_order(root_set,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
	#endif
	
	/****/
	//build database
	//build field
	ierr = m_field.build_fields(); CHKERRQ(ierr);
	#ifdef HOON
	Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
	ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
	#endif
	//build finite elemnts
	ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
	//build adjacencies
	ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
	//build problem
	ierr = m_field.build_problems(); CHKERRQ(ierr);
	
	/****/
	//mesh partitioning 
	//partition
	ierr = m_field.simple_partition_problem("EX1_PROBLEM"); CHKERRQ(ierr);
	ierr = m_field.partition_finite_elements("EX1_PROBLEM"); CHKERRQ(ierr);
	ierr = m_field.simple_partition_problem("EX2_PROBLEM"); CHKERRQ(ierr);
	ierr = m_field.partition_finite_elements("EX2_PROBLEM"); CHKERRQ(ierr);
	//what are ghost nodes, see Petsc Manual
	ierr = m_field.partition_ghost_dofs("EX1_PROBLEM"); CHKERRQ(ierr);
	ierr = m_field.partition_ghost_dofs("EX2_PROBLEM"); CHKERRQ(ierr);
	
	Mat A;
	ierr = m_field.MatCreateMPIAIJWithArrays("EX1_PROBLEM",&A); CHKERRQ(ierr);
	Vec D,F;
	ierr = m_field.VecCreateGhost("EX1_PROBLEM",ROW,&F); CHKERRQ(ierr);
	ierr = m_field.VecCreateGhost("EX1_PROBLEM",COL,&D); CHKERRQ(ierr);
	
	Mat B;
	ierr = m_field.MatCreateMPIAIJWithArrays("EX2_PROBLEM",&B); CHKERRQ(ierr);
	Vec C,G;
	ierr = m_field.VecCreateGhost("EX2_PROBLEM",ROW,&G); CHKERRQ(ierr);
	ierr = m_field.VecCreateGhost("EX2_PROBLEM",COL,&C); CHKERRQ(ierr);
	
	//Extract data from MAT_HELMHOLTZ block
	//Exact solution of Impinging sphere from Acoustic isogeometric boundary element analysis by R.N. Simpson etc.
	static double aNgularfreq;
	static double sPeed; //Without static. I got error:use of local variable with automatic storage from containing function
	
	
	//for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field,"MAT_HELMHOLTZ",it)) {
	for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it))
	{
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
			std::string wait;
			std::cout << "\n sPeed = \n" << sPeed << "\n aNgularfreq = " << aNgularfreq << std::endl;
			
		}
	}
	
	
	double wavenumber = aNgularfreq/sPeed;
	
	std::cout << "\n I am here !!!!!!! WAVE NUMBER = \n" << wavenumber << std::endl;
	
	{
		
		MyFunApprox_re function_evaluator_re;
		FieldApproximationH1<MyFunApprox_re> field_approximation_re(m_field);
		
		field_approximation_re.loopMatrixAndVector(
			"EX1_PROBLEM","FE1","reEX",A,F,function_evaluator_re);
	}
	
	{
		
		MyFunApprox_im function_evaluator_im;
		FieldApproximationH1<MyFunApprox_im> field_approximation_im(m_field);
		
		field_approximation_im.loopMatrixAndVector(
			"EX2_PROBLEM","FE2","imEX",B,G,function_evaluator_im);
	}
	
	//solve real part of the acoustic problem
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	
	KSP solver1;
	ierr = KSPCreate(PETSC_COMM_WORLD,&solver1); CHKERRQ(ierr);
	ierr = KSPSetOperators(solver1,A,A); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(solver1); CHKERRQ(ierr);
	ierr = KSPSetUp(solver1); CHKERRQ(ierr);
	
	ierr = KSPSolve(solver1,F,D); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	
	//solve imagine part of the acoustic problem
	ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(G,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(G,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(G,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(G,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	
	KSP solver2;
	ierr = KSPCreate(PETSC_COMM_WORLD,&solver2); CHKERRQ(ierr);
	ierr = KSPSetOperators(solver2,B,B); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(solver2); CHKERRQ(ierr);
	ierr = KSPSetUp(solver2); CHKERRQ(ierr);
	
	ierr = KSPSolve(solver2,G,C); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	
	ierr = m_field.set_global_VecCreateGhost("EX1_PROBLEM",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = m_field.set_global_VecCreateGhost("EX2_PROBLEM",COL,C,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	//destroy the solvers
	ierr = KSPDestroy(&solver1); CHKERRQ(ierr);
	ierr = VecDestroy(&D); CHKERRQ(ierr);
	ierr = VecDestroy(&F); CHKERRQ(ierr);
	ierr = MatDestroy(&A); CHKERRQ(ierr);
	ierr = KSPDestroy(&solver2); CHKERRQ(ierr);
	ierr = VecDestroy(&C); CHKERRQ(ierr);
	ierr = VecDestroy(&G); CHKERRQ(ierr);
	ierr = MatDestroy(&B); CHKERRQ(ierr);
	
	
	////Define rows/cols and element data
	//ierr = m_field.modify_finite_element_add_field_row("FE2","imEX"); CHKERRQ(ierr);
	//ierr = m_field.modify_finite_element_add_field_col("FE2","imEX"); CHKERRQ(ierr);
	//ierr = m_field.modify_finite_element_add_field_data("FE1","imEX"); CHKERRQ(ierr);
	//#ifdef HOON
	//ierr = m_field.modify_finite_element_add_field_data("FE2","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	//#endif
	
	
	PostPocOnRefinedMesh post_proc1(m_field);
	ierr = post_proc1.generateRefereneElemenMesh(); CHKERRQ(ierr);
	ierr = post_proc1.addFieldValuesPostProc("reEX"); CHKERRQ(ierr);
	ierr = post_proc1.addFieldValuesGradientPostProc("reEX"); CHKERRQ(ierr);
	ierr = post_proc1.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	ierr = m_field.loop_finite_elements("EX1_PROBLEM","FE1",post_proc1); CHKERRQ(ierr);
	rval = post_proc1.postProcMesh.write_file("real_out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
	
	PostPocOnRefinedMesh post_proc2(m_field);
	ierr = post_proc2.generateRefereneElemenMesh(); CHKERRQ(ierr);
	ierr = post_proc2.addFieldValuesPostProc("imEX"); CHKERRQ(ierr);
	ierr = post_proc2.addFieldValuesGradientPostProc("imEX"); CHKERRQ(ierr);
	ierr = post_proc2.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	ierr = m_field.loop_finite_elements("EX2_PROBLEM","FE2",post_proc2); CHKERRQ(ierr);
	rval = post_proc2.postProcMesh.write_file("imag_out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
	//PostPocOnRefinedMesh post_proc2(m_field);
	//ierr = post_proc2.generateRefereneElemenMesh(); CHKERRQ(ierr);
	//ierr = post_proc2.addFieldValuesPostProc("imEX"); CHKERRQ(ierr);
	//ierr = post_proc2.addFieldValuesGradientPostProc("imEX"); CHKERRQ(ierr);
	//ierr = post_proc2.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	//ierr = m_field.loop_finite_elements("EX2_PROBLEM","FE2",post_proc2); CHKERRQ(ierr);
	//rval = post_proc2.postProcMesh.write_file("out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
	
	
    
	char command1[] = "mbconvert ./real_out.h5m ./real_new_out.vtk && cp ./real_new_out.vtk ../../../../../mnt/home/Documents/mofem-cephas/mofem_v0.2/users_modules/helmholtz/ && mbconvert ./imag_out.h5m ./imag_new_out.vtk && cp ./imag_new_out.vtk ../../../../../mnt/home/Documents/mofem-cephas/mofem_v0.2/users_modules/helmholtz/";
	//char command2[] = "mbconvert ./imag_out.h5m ./imag_new_out.vtk && cp ./new_out.vtk ../../../../../mnt/home/Documents/mofem-cephas/mofem_v0.2/users_modules/helmholtz/";


	int todo1 = system( command1 );
	//int todo2 = system( command2 ); 


	
	
	
	
	
	
	
	//typedef tee_device<ostream, ofstream> TeeDevice;
	//typedef stream<TeeDevice> TeeStream;

	//ofstream ofs("acoustic_re_field_testing_field_approximation.txt");
	//TeeDevice tee(cout, ofs); 
	//TeeStream my_split(tee);

	//Range nodes;
	//rval = moab.get_entities_by_type(0,MBVERTEX,nodes,true); CHKERR(rval);
	//ublas::matrix<double> nodes_vals;
	//nodes_vals.resize(nodes.size(),3);  //change the parameter from 3 to 1 ?
	//rval = moab.tag_get_data(
	//		   ent_method_field1_on_10nodeTet.th,nodes,&*nodes_vals.data().begin()); CHKERR(rval);
	

	//const double eps = 1e-4;

	//my_split.precision(3);
	//my_split.setf(std::ios::fixed);
	//for(
	//	ublas::unbounded_array<double>::iterator it = nodes_vals.data().begin();
	//	it!=nodes_vals.data().end();it++) {
	//	*it = fabs(*it)<eps ? 0.0 : *it; //if a < b ?then c, :else d
	//}
	//my_split << nodes_vals << endl;

	//const MoFEMProblem *problemPtr;
	//ierr = m_field.get_problem("PROBLEM1",&problemPtr); CHKERRQ(ierr);
	//map<EntityHandle,double> m0,m1,m2;
	//for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(problemPtr,dit)) {

	//	my_split.precision(3);
	//	my_split.setf(std::ios::fixed);
	//	double val = fabs(dit->get_FieldData())<eps ? 0.0 : dit->get_FieldData();
	//	my_split << dit->get_petsc_gloabl_dof_idx() << " " << val << endl;

	//}

	ierr = PetscFinalize(); CHKERRQ(ierr);
	
	return 0;

}
	
	
	












