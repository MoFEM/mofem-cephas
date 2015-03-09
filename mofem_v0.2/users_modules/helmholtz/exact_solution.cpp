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
#include <ExactElement.hpp>

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
	bool useApproxError;
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
	MoFEM::Core core(moab);
	FieldInterface& m_field = core;
	
	//set entitities bit level
	BitRefLevel bit_level0;
	bit_level0.set(0);
	EntityHandle meshset_level0;
	rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
	ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
	

	ierr = m_field.add_field("reAY",H1,1); CHKERRQ(ierr);
	ierr = m_field.add_field("imAY",H1,1); CHKERRQ(ierr);
	
	//Problem
	ierr = m_field.add_problem("AY1_PROBLEM"); CHKERRQ(ierr); //analytical for real field
	ierr = m_field.add_problem("AY2_PROBLEM"); CHKERRQ(ierr); //analytical for imag field
	
	ierr = m_field.modify_problem_ref_level_add_bit("AY1_PROBLEM",bit_level0); CHKERRQ(ierr);  //analytical Dirichlet
	ierr = m_field.modify_problem_ref_level_add_bit("AY2_PROBLEM",bit_level0); CHKERRQ(ierr);  //analytical Dirichlet
	
	//meshset consisting all entities in mesh
	EntityHandle root_set = moab.get_root_set(); 
	//add entities to field
	ierr = m_field.add_ents_to_field_by_TETs(root_set,"reAY"); CHKERRQ(ierr);
	ierr = m_field.add_ents_to_field_by_TETs(root_set,"imAY"); CHKERRQ(ierr);
	
	
	//set app. order
	//see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
	int order = 3;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) {
		order = 3;
	}


	ierr = m_field.set_field_order(root_set,MBTET,"reAY",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBTRI,"reAY",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBEDGE,"reAY",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBVERTEX,"reAY",1); CHKERRQ(ierr);
	
	ierr = m_field.set_field_order(root_set,MBTET,"imAY",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBTRI,"imAY",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBEDGE,"imAY",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBVERTEX,"imAY",1); CHKERRQ(ierr);
	
	ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
	ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
	
	Range tets;
	for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field,"EXACT_SOL",it)) {
		rval = moab.get_entities_by_type(it->get_meshset(),MBTET,tets,true); CHKERR_PETSC(rval);
	}
	
	AnalyticalSolution analytical_1(m_field,tets);
	AnalyticalSolution analytical_2(m_field,tets);
	
	ierr = analytical_1.initializeExactProblem("AY1_PROBLEM","FE1","reAY"); CHKERRQ(ierr);
	ierr = analytical_2.initializeExactProblem("AY2_PROBLEM","FE2","imAY"); CHKERRQ(ierr);
	ierr = m_field.modify_finite_element_add_field_data("FE1","imAY"); CHKERRQ(ierr);

	/****/
	//build database
	//build field
	ierr = m_field.build_fields(); CHKERRQ(ierr);
	//build finite elemnts
	ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
	//build adjacencies
	ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
	//build problem
	ierr = m_field.build_problems(); CHKERRQ(ierr);
	
	Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
	ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
	
	ierr = m_field.partition_problem("AY1_PROBLEM"); CHKERRQ(ierr);
	ierr = m_field.partition_finite_elements("AY1_PROBLEM"); CHKERRQ(ierr);
	ierr = m_field.partition_ghost_dofs("AY1_PROBLEM"); CHKERRQ(ierr);
	
	ierr = m_field.partition_problem("AY2_PROBLEM"); CHKERRQ(ierr);
	ierr = m_field.partition_finite_elements("AY2_PROBLEM"); CHKERRQ(ierr);
	ierr = m_field.partition_ghost_dofs("AY2_PROBLEM"); CHKERRQ(ierr);
	
	//Exact solution of Impinging sphere from Acoustic isogeometric boundary element analysis by R.N. Simpson etc.
    static double aNgularfreq;
	static double sPeed; //Without static. I got error:use of local variable with automatic storage from containing function
	
	for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it))
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
	
	
	
	/* this function compute the scattered field of helmholtz operator */
	struct AnaliticalFunction {
		static double fUN_real(double x,double y,double z) {
			
			const double pi = atan( 1.0 ) * 4.0;
			//double R = sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)); //radius
			//double theta = atan2(y,x)+2*pi; //the arctan of radians (y/x)		  
			const double wAvenumber = aNgularfreq/sPeed;
			
			const double k = wAvenumber;  //Wave number
			//const double a = 0.5;         //radius of the sphere,wait to modify by user
			//const double const1 = k * a;
			//double const2 = k * R;
			
			
			const complex< double > i( 0.0, 1.0 );
			
			//// magnitude of incident wave
			//const double phi_incident_mag = 1.0;
			
			//const double tol = 1.0e-10;
			//double max = 0.0;
			//double min = 999999.0;
			
			complex< double > result = 0.0;
			//complex< double > prev_result;
			
			//double error = 100.0;
			//unsigned int n = 0; //initialized the infinite series loop
			
			//while( error > tol )  //finding the acoustic potential in one single point.
			//{
			//	  double jn_der = n / const1 * sph_bessel( n, const1 ) - sph_bessel( n + 1, const1 );  //The derivative of bessel function
			//	  complex< double > hn_der = n / const1 * sph_hankel_1( n, const1 ) - sph_hankel_1( n + 1, const1 );
			//	  //complex< double > hn_der = 0.5 * ( sph_hankel_1( n - 1, const1 ) -
			//	  //( sph_hankel_1( n, const1 ) + const1 * sph_hankel_1( n + 1, const1 ) ) / const1 );
			//	  double Pn = legendre_p( n, cos( theta ) );
			//	  complex< double >hn = sph_hankel_1( n, const2 );  //S Hankel first kind function
			//	  prev_result = result;
			//	  result -= pow( i, n ) * ( 2.0 * n + 1.0 ) * jn_der / hn_der * Pn * hn;
			//	  error = abs( abs( result ) - abs( prev_result ) );
			//	  ++n;
			//}
			
			//const complex< double > inc_field = exp( i * k * R * cos( theta ) );  //incident wave
			//const complex< double > total_field = inc_field + result;
			
			//double val = std::real(result);
			//ofs << theta << "\t" << abs( result ) << "\t" << abs( inc_field ) << "\t" << abs( total_field ) << "\t" << R << endl; //write the file
			double theta = pi/4;
			result = exp(i*(k*cos(theta)*x+k*sin(theta)*y));
			
			return std::real(result);
			//return std::real((exp(i*k*x)-1-i*exp(i*k)*sin(k*x))/(pow(k,2.0))); //exact solution of 1D problem
			//return 0;
		}
		static double fUN_imag(double x,double y,double z) {
			const double pi = atan( 1.0 ) * 4.0;
			//double R = sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)); //radius
			//double theta = atan2(y,x) + 2*pi; //the arctan of radians (y/x)
			const double wAvenumber = aNgularfreq/sPeed;
			
			const double k = wAvenumber;  //Wave number
			//  const double a = 0.5;         //radius of the sphere,wait to modify by user
			//  const double const1 = k * a;
			//  double const2 = k * R;
			//  
			//  
			const complex< double > i( 0.0, 1.0 );
			// 
			//  // magnitude of incident wave
			//  const double phi_incident_mag = 1.0;
			//  
			//  const double tol = 1.0e-10;
			//  double max = 0.0;
			//  double min = 999999.0;
			//  
			complex< double > result = 0.0;
			//  complex< double > prev_result;
			//  
			//  double error = 100.0;
			//  unsigned int n = 0; //initialized the infinite series loop
			//  
			//  while( error > tol )  //finding the acoustic potential in one single point.
			//  {
			//	  double jn_der = n / const1 * sph_bessel( n, const1 ) - sph_bessel( n + 1, const1 );  //The derivative of bessel function
			//	  complex< double > hn_der = n / const1 * sph_hankel_1( n, const1 ) - sph_hankel_1( n + 1, const1 );
			//	  double Pn = legendre_p( n, cos( theta ) );
			//	  complex< double >hn = sph_hankel_1( n, const2 );  //S Hankel first kind function
			//	  prev_result = result;
			//	  result -= pow( i, n ) * ( 2.0 * n + 1.0 ) * jn_der / hn_der * Pn * hn;
			//	  error = abs( abs( result ) - abs( prev_result ) );
			//	  ++n;
			//  }
			
			//const complex< double > inc_field = exp( i * k * R * cos( theta ) );  //incident wave
			double theta = pi/4;
			result = exp(i*(k*cos(theta)*x+k*sin(theta)*y));
			return std::imag(result);	  
			//return std::imag((exp(i*k*x)-1-i*exp(i*k)*sin(k*x))/(pow(k,2.0))); //exact solution of 1D problem
			//return 0;
		}
	};
	

	Vec M,P;
	ierr = m_field.VecCreateGhost("AY1_PROBLEM",COL,&M); CHKERRQ(ierr);
	ierr = m_field.VecCreateGhost("AY2_PROBLEM",COL,&P); CHKERRQ(ierr);
	

	

	ierr = analytical_1.setExactSolRhsOp("reAY",AnaliticalFunction::fUN_real,M); CHKERRQ(ierr);
	ierr = analytical_2.setExactSolRhsOp("imAY",AnaliticalFunction::fUN_imag,P); CHKERRQ(ierr);
		
		

	
	//ierr = m_field.set_global_VecCreateGhost("AY1_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);	
	//ierr = m_field.set_global_VecCreateGhost("AY2_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);	

	
	ierr = analytical_1.solveExactProblem("AY1_PROBLEM","FE1","reAY",M); CHKERRQ(ierr);
	ierr = analytical_2.solveExactProblem("AY2_PROBLEM","FE2","imAY",P); CHKERRQ(ierr);
	
	//double nrml2_M,nrml2_P;
	//ierr = VecNorm(M,NORM_2,&nrml2_M); CHKERRQ(ierr);
	//ierr = VecNorm(P,NORM_2,&nrml2_P); CHKERRQ(ierr);

	//std::cout << "\n the L2 of real analytical field is : \n" << nrml2_M << std::endl;
	//std::cout << "\n the L2 of imag analytical field is : \n" << nrml2_P << std::endl;

	
	
	//if(pcomm->rank()==0) {
	rval = moab.write_file("exact_solution_mesh.h5m"); CHKERR_PETSC(rval);
	//}
	//destroy the solvers
	ierr = VecDestroy(&M); CHKERRQ(ierr);
	ierr = VecDestroy(&P); CHKERRQ(ierr);

	
	
	PostPocOnRefinedMesh post_proc1(m_field);
	ierr = post_proc1.generateRefereneElemenMesh(); CHKERRQ(ierr);
	ierr = post_proc1.addFieldValuesPostProc("reAY"); CHKERRQ(ierr);
	//ierr = post_proc1.addFieldValuesGradientPostProc("reAY"); CHKERRQ(ierr);
	ierr = post_proc1.addFieldValuesPostProc("imAY"); CHKERRQ(ierr);
	//ierr = post_proc1.addFieldValuesGradientPostProc("imAY"); CHKERRQ(ierr);
	ierr = post_proc1.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	ierr = m_field.loop_finite_elements("AY1_PROBLEM","FE1",post_proc1); CHKERRQ(ierr);
	rval = post_proc1.postProcMesh.write_file("exact_out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);

	//output the results from Docker
	char command1[] = "mbconvert ./exact_out.h5m ./exact_out.vtk && cp ./exact_out.vtk ../../../../../mnt/home/Desktop/U_pan/helmholtz_results/";
	int todo1 = system( command1 );
	
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
	
