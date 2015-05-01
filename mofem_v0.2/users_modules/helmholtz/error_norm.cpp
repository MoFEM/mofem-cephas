/* 
  \file error_norm.cpp
  \ingroup mofem_helmholtz_elem

  Calculates finite element (Galerkin) approximation for difference between two solutions in L^2 and H_1 norm. 
  \bug work for scalar filed only, NEED further modification. As well as calculate A matrix twice.
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
#include <PotsProcOnRefMesh.hpp>
#include <NormElement.hpp>

#include <Projection10NodeCoordsOnField.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <petsctime.h>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {
	
	ErrorCode rval;
	PetscErrorCode ierr;
	
	PetscInitialize(&argc,&argv,(char *)0,help);
	
	moab::Core mb_instance;
	Interface& moab = mb_instance;
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	//read h5m solution file into moab
	bool usel2; //norm type
	PetscBool flg = PETSC_TRUE;
	
	char mesh_file_name[255];
	ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) {
		SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
	}
	
	
	PetscBool userela = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,"-relative_error",&userela,NULL); CHKERRQ(ierr);
	/* FIX ME, change to enum for better performance */
	char type_error_norm[255];
	ierr = PetscOptionsGetString(PETSC_NULL,"-norm_type",type_error_norm,255,&flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) {
		SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -type_error_norm (L2 or H1 type needed)");
	}

	if (strcmp ("l2",type_error_norm ) == 0) {usel2 = true;}
	else if(strcmp ("h1",type_error_norm ) == 0) {usel2 = false;}
	
	ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
	if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
	
	const char *option;
	option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
	BARRIER_RANK_START(pcomm) 
	/* load the mesh files */
	rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
	//rval = moab.load_file(mesh_file_name2, 0, option); CHKERR_PETSC(rval); 
	BARRIER_RANK_END(pcomm) 

	//Create MoFEM (Joseph) database
	MoFEM::Core core(moab);
	FieldInterface& m_field = core;

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
	ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

	//Fields
	ierr = m_field.add_field("erorNORM_re",H1,1); CHKERRQ(ierr);
	ierr = m_field.add_field("erorNORM_im",H1,1); CHKERRQ(ierr);

	
	//Problem
	ierr = m_field.add_problem("NORM_PROBLEM1"); CHKERRQ(ierr);
	ierr = m_field.add_problem("NORM_PROBLEM2"); CHKERRQ(ierr);


	//set refinment level for problem
	ierr = m_field.modify_problem_ref_level_add_bit("NORM_PROBLEM1",bit_level0); CHKERRQ(ierr);
	ierr = m_field.modify_problem_ref_level_add_bit("NORM_PROBLEM2",bit_level0); CHKERRQ(ierr);



	//meshset consisting all entities in mesh
	EntityHandle root_set = moab.get_root_set(); 
	//add entities to field
	ierr = m_field.add_ents_to_field_by_TETs(root_set,"erorNORM_re"); CHKERRQ(ierr);
	ierr = m_field.add_ents_to_field_by_TETs(root_set,"erorNORM_im"); CHKERRQ(ierr);

	
	//set app. order , approximation error of norm should be as least 1 order higher than numerical spaces.
	//see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
	PetscInt order;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) {
		order = 2;
	}
	
	ierr = m_field.set_field_order(root_set,MBTET,"erorNORM_re",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBTRI,"erorNORM_re",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBEDGE,"erorNORM_re",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBVERTEX,"erorNORM_re",1); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBTET,"erorNORM_im",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBTRI,"erorNORM_im",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBEDGE,"erorNORM_im",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBVERTEX,"erorNORM_im",1); CHKERRQ(ierr);
	
	if(!m_field.check_field("MESH_NODE_POSITIONS")) {
	ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
	ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
	}
	
    double real_error;
    double imag_error;
	NormElement norm_elements_re(m_field,real_error);
	NormElement norm_elements_im(m_field,imag_error);

	if(m_field.check_field("reEX") && m_field.check_field("rePRES")) {
	//&& m_field.check_field("imPRES") && m_field.check_field("imEX") ) {
		norm_elements_re.addNormElements("NORM_PROBLEM1","NORM_FE1","erorNORM_re","reEX","rePRES");
		norm_elements_im.addNormElements("NORM_PROBLEM2","NORM_FE2","erorNORM_im","imEX","imPRES");
		ierr = m_field.modify_finite_element_add_field_data("NORM_FE1","erorNORM_im"); CHKERRQ(ierr);
	}
	

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
	
	
	/****/
	//mesh partitioning 
	//partition
	ierr = m_field.partition_problem("NORM_PROBLEM1"); CHKERRQ(ierr);
	ierr = m_field.partition_finite_elements("NORM_PROBLEM1"); CHKERRQ(ierr);
	//what are ghost nodes, see Petsc Manual
	ierr = m_field.partition_ghost_dofs("NORM_PROBLEM1"); CHKERRQ(ierr);
	
	ierr = m_field.partition_problem("NORM_PROBLEM2"); CHKERRQ(ierr);
	ierr = m_field.partition_finite_elements("NORM_PROBLEM2"); CHKERRQ(ierr);
	ierr = m_field.partition_ghost_dofs("NORM_PROBLEM2"); CHKERRQ(ierr);
    
	if(m_field.check_field("reEX") && m_field.check_field("imEX")) {

      ierr = m_field.build_problem("EX1_PROBLEM"); CHKERRQ(ierr);
      ierr = m_field.partition_problem("EX1_PROBLEM"); CHKERRQ(ierr);
      ierr = m_field.partition_finite_elements("EX1_PROBLEM"); CHKERRQ(ierr);
      ierr = m_field.partition_ghost_dofs("EX1_PROBLEM"); CHKERRQ(ierr);

	}
	
	//print block sets with materials
	//ierr = m_field.print_cubit_materials_set(); CHKERRQ(ierr);
	
	Vec F;
	ierr = m_field.VecCreateGhost("NORM_PROBLEM1",ROW,&F); CHKERRQ(ierr);
	Vec T;
	ierr = VecDuplicate(F,&T); CHKERRQ(ierr);
	Mat A;
	ierr = m_field.MatCreateMPIAIJWithArrays("NORM_PROBLEM1",&A); CHKERRQ(ierr);
	
	Vec G;
	ierr = m_field.VecCreateGhost("NORM_PROBLEM2",ROW,&G); CHKERRQ(ierr);
	Vec D;
	ierr = VecDuplicate(G,&D); CHKERRQ(ierr);


	/*   Set operators  */
    ierr = norm_elements_re.setNormFiniteElementRhsOperator("erorNORM_re","reEX","rePRES",F,usel2,userela); CHKERRQ(ierr);
    ierr = norm_elements_re.setNormFiniteElementLhsOperator("erorNORM_re","reEX","rePRES",A); CHKERRQ(ierr);
    
	ierr = norm_elements_im.setNormFiniteElementRhsOperator("erorNORM_im","imEX","imPRES",G,usel2,userela); CHKERRQ(ierr);

    
    
    /* create matrix and vectors */
	ierr = VecZeroEntries(T); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecZeroEntries(F); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = MatZeroEntries(A); CHKERRQ(ierr);
	ierr = m_field.set_global_ghost_vector("NORM_PROBLEM1",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

	ierr = VecZeroEntries(D); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecZeroEntries(G); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(G,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(G,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	ierr = m_field.set_global_ghost_vector("NORM_PROBLEM2",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

    /* looop over finite elements */
	ierr = m_field.loop_finite_elements("NORM_PROBLEM1","NORM_FE1",norm_elements_re.getLoopFe()); CHKERRQ(ierr);
	ierr = m_field.loop_finite_elements("NORM_PROBLEM2","NORM_FE2",norm_elements_im.getLoopFe()); CHKERRQ(ierr);
      
      
      /* create ghost points */
	ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	
	ierr = VecGhostUpdateBegin(G,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(G,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(G); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(G); CHKERRQ(ierr);

	
	//Solver
	KSP solver1;
	ierr = KSPCreate(PETSC_COMM_WORLD,&solver1); CHKERRQ(ierr);
	ierr = KSPSetOperators(solver1,A,A); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(solver1); CHKERRQ(ierr);
	ierr = KSPSetUp(solver1); CHKERRQ(ierr);
	
	ierr = KSPSolve(solver1,F,T); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = m_field.set_global_ghost_vector("NORM_PROBLEM1",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
   
    
	
	ierr = KSPSolve(solver1,G,D); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = m_field.set_global_ghost_vector("NORM_PROBLEM2",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	
	
	
	/* Global error calculation */
	//PetscScalar sum_T,sum_D;
	PetscReal nrm2_T,nrm2_D,nrm2_P,nrm2_M;
	/* FIXME sum and sqrt of Vector T and D */
	//ierr = VecSum(T,&sum_T);
	//ierr = VecSum(D,&sum_D);
	//nrm2_T = sqrt(abs(sum_T));
	//nrm2_D = sqrt(abs(sum_D));
	ierr = VecNorm(T,NORM_2,&nrm2_T);;
	ierr = VecNorm(D,NORM_2,&nrm2_D); CHKERRQ(ierr);
	//ierr = VecNorm(T,NORM_MAX,&pointwisenormT);
    
    Vec P,M;
	ierr = m_field.VecCreateGhost("EX1_PROBLEM",ROW,&M); CHKERRQ(ierr);
    ierr = VecDuplicate(M,&P); CHKERRQ(ierr);
    
	ierr = m_field.set_local_ghost_vector("EX1_PROBLEM",ROW,M,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = m_field.set_other_global_ghost_vector("EX1_PROBLEM","reEX","imEX",ROW,P,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    if(usel2) {
      ierr = VecNorm(M,NORM_INFINITY,&nrm2_M);
      ierr = VecNorm(P,NORM_INFINITY,&nrm2_P);
      PetscPrintf(PETSC_COMM_WORLD,"\n real part of global l2 error is: %f \n",sqrt(real_error)/nrm2_M);
      PetscPrintf(PETSC_COMM_WORLD,"\n imag part of global l2 error is: %f \n",sqrt(imag_error)/nrm2_P);
    } else {
      ierr = VecNorm(M,NORM_INFINITY,&nrm2_M);
      ierr = VecNorm(P,NORM_INFINITY,&nrm2_P);
      PetscPrintf(PETSC_COMM_WORLD,"\n real part of global H1 error is: %f \n",sqrt(real_error)/nrm2_M);
      PetscPrintf(PETSC_COMM_WORLD,"\n imag part of global H1 error is: %f \n",sqrt(imag_error)/nrm2_P);
    }
    	
	ierr = VecNorm(M,NORM_FROBENIUS,&nrm2_M);
	ierr = VecNorm(P,NORM_2,&nrm2_P); CHKERRQ(ierr);
	//ierr = VecNorm(M,NORM_INFINITY,&nrm2_M);
    
	
	//out stream the global error
	if(usel2 && !userela) {
      PetscPrintf(PETSC_COMM_WORLD,"\n The Global L2 error of real field is : %f \n",nrm2_T);
      PetscPrintf(PETSC_COMM_WORLD,"\n The Global L2 error of imag field is : %f \n",nrm2_D);

      PetscPrintf(PETSC_COMM_WORLD,"\n The Global L2 relative error of real field is : %f \n",nrm2_T/nrm2_M);
      PetscPrintf(PETSC_COMM_WORLD,"\n The Global L2 relative error of imag field is : %f \n",nrm2_D/nrm2_P);
      PetscPrintf(PETSC_COMM_WORLD,"\n Global error  of total potential in l2 norm is : %f \n",sqrt((nrm2_T/nrm2_M)*(nrm2_T/nrm2_M) + (nrm2_D/nrm2_P)*(nrm2_D/nrm2_P)));
      //std::cout << "\n The Global L2 relative error of real field is : --\n" << nrm2_T/nrm2_M  << std::endl;
      //std::cout << "\n The Global L2 relative error of imag field is  : --\n" << nrm2_D/nrm2_P << std::endl;
      //cout << "\n Global error  of total potential in l2 norm  \n" << sqrt((nrm2_T/nrm2_M)*(nrm2_T/nrm2_M) + (nrm2_D/nrm2_P)*(nrm2_D/nrm2_P)) << endl;
	}
	else if(!usel2 && !userela) {

      std::cout << "\n The Global H1 relative error of real field is : --\n" << nrm2_T/nrm2_M  << std::endl;
      std::cout << "\n The Global H1 relative error of imag field is  : --\n" << nrm2_D/nrm2_P << std::endl;
      cout << "\n Global error  of total potential in l2 norm  \n" << sqrt((nrm2_T/nrm2_M)*(nrm2_T/nrm2_M) + (nrm2_D/nrm2_P)*(nrm2_D/nrm2_P)) << endl;
	}
	else if(userela) {
		//NEED TO BE FIXED
		//we need two vector, one is exact solution, one is error in the norm, then find the maximum of 
		//two vectors, and devide second one by first one. it is the global pointwise relative error.
		std::cout << "\n The Global least square of H1 Norm of error real field is : --\n" << nrm2_T << std::endl;
		//std::cout << "\n The Global Pointwise of H1 Norm of error for real field is : --\n" << pointwisenorm << std::endl;
	}


	/*  destroy objects  */
	ierr = MatDestroy(&A); CHKERRQ(ierr);
	ierr = VecDestroy(&F); CHKERRQ(ierr);
	ierr = VecDestroy(&T); CHKERRQ(ierr);
	ierr = KSPDestroy(&solver1); CHKERRQ(ierr);
	
	ierr = MatDestroy(&B); CHKERRQ(ierr);
	ierr = VecDestroy(&G); CHKERRQ(ierr);
	ierr = VecDestroy(&D); CHKERRQ(ierr);
	//ierr = KSPDestroy(&solver2); CHKERRQ(ierr);

	ierr = VecDestroy(&M); CHKERRQ(ierr);
	ierr = VecDestroy(&P); CHKERRQ(ierr);
	
	PostPocOnRefinedMesh post_proc1(m_field);
	ierr = post_proc1.generateRefereneElemenMesh(); CHKERRQ(ierr);
	ierr = post_proc1.addFieldValuesPostProc("erorNORM_re"); CHKERRQ(ierr);
	ierr = post_proc1.addFieldValuesPostProc("erorNORM_im"); CHKERRQ(ierr);
	ierr = post_proc1.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	ierr = m_field.loop_finite_elements("NORM_PROBLEM1","NORM_FE1",post_proc1); CHKERRQ(ierr);
	rval = post_proc1.postProcMesh.write_file("norm_error.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);

	ierr = PetscTime(&v2);CHKERRQ(ierr);
	ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f S CPU Time = %f S \n",pcomm->rank(),v2-v1,t2-t1);
	
	
	//output the results from Docker
	//char command1[] = "mbconvert norm_error.h5m ./norm_error.vtk && cp ./norm_error.vtk ../../../../mnt/home/Desktop/U_pan/helmholtz_results/";
	//int todo1 = system( command1 );
	
	
	
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;

}





