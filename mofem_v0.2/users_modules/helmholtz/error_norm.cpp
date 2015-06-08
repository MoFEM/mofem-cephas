/*
  \file error_norm.cpp
  \ingroup mofem_helmholtz_elem

  Calculates finite element (Galerkin) approximation for difference between two solutions in L^2 and H_1 norm.
  \bug work for scalar filed only, NEED further modification.
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
#include <PostProcOnRefMesh.hpp>
#include <NormElement.hpp>

#include <Projection10NodeCoordsOnField.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <petsctime.h>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>

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


	//PetscBool userela = PETSC_FALSE;
	//ierr = PetscOptionsGetBool(NULL,"-relative_error",&userela,NULL); CHKERRQ(ierr);

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
	ierr = m_field.add_field("erorNORM",H1,1); CHKERRQ(ierr);

	//Problem
	ierr = m_field.add_problem("NORM_PROBLEM"); CHKERRQ(ierr);

	//set refinment level for problem
	ierr = m_field.modify_problem_ref_level_add_bit("NORM_PROBLEM",bit_level0); CHKERRQ(ierr);




	//meshset consisting all entities in mesh
	EntityHandle root_set = moab.get_root_set();
	//add entities to field
	ierr = m_field.add_ents_to_field_by_TETs(root_set,"erorNORM"); CHKERRQ(ierr);


	//set app. order , approximation error of norm should be as least 1 order higher than numerical spaces.
	//see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
	PetscInt order;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) {
		order = 2;
	}

	ierr = m_field.set_field_order(root_set,MBTET,"erorNORM",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBTRI,"erorNORM",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBEDGE,"erorNORM",order); CHKERRQ(ierr);
	ierr = m_field.set_field_order(root_set,MBVERTEX,"erorNORM",1); CHKERRQ(ierr);

	if(!m_field.check_field("MESH_NODE_POSITIONS")) {
	ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
	ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
	}

	/* global error and analytical solution */
  double error;
  double analy;


	NormElement norm_elements(m_field,error,analy);

	if(m_field.check_field("reEX") && m_field.check_field("rePRES")) {
	//&& m_field.check_field("imPRES") && m_field.check_field("imEX") ) {
		norm_elements.addNormElements("NORM_PROBLEM","NORM_FE","erorNORM","reEX","rePRES","imEX","imPRES");
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
	ierr = m_field.partition_problem("NORM_PROBLEM"); CHKERRQ(ierr);
	ierr = m_field.partition_finite_elements("NORM_PROBLEM"); CHKERRQ(ierr);
	//what are ghost nodes, see Petsc Manual
	ierr = m_field.partition_ghost_dofs("NORM_PROBLEM"); CHKERRQ(ierr);

	if(m_field.check_field("reEX") && m_field.check_field("imEX")) {

    ierr = m_field.build_problem("EX1_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_problem("EX1_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("EX1_PROBLEM"); CHKERRQ(ierr);
    ierr = m_field.partition_ghost_dofs("EX1_PROBLEM"); CHKERRQ(ierr);

	}

	//print block sets with materials
	//ierr = m_field.print_cubit_materials_set(); CHKERRQ(ierr);

	Vec F;
	ierr = m_field.VecCreateGhost("NORM_PROBLEM",ROW,&F); CHKERRQ(ierr);
	Vec T;
	ierr = VecDuplicate(F,&T); CHKERRQ(ierr);
	Mat A;
	ierr = m_field.MatCreateMPIAIJWithArrays("NORM_PROBLEM",&A); CHKERRQ(ierr);


	/*   Set operators  */
  ierr = norm_elements.setNormFiniteElementRhsOperator("erorNORM","reEX","rePRES","imEX","imPRES",F,usel2); CHKERRQ(ierr);
  ierr = norm_elements.setNormFiniteElementLhsOperator("erorNORM",A); CHKERRQ(ierr);



    /* create matrix and vectors */
	ierr = VecZeroEntries(T); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecZeroEntries(F); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = MatZeroEntries(A); CHKERRQ(ierr);
	ierr = m_field.set_global_ghost_vector("NORM_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

    /* looop over finite elements */
	ierr = m_field.loop_finite_elements("NORM_PROBLEM","NORM_FE",norm_elements.getLoopFe()); CHKERRQ(ierr);


      /* create ghost points */
	ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	//Solver
	KSP solver1;
	ierr = KSPCreate(PETSC_COMM_WORLD,&solver1); CHKERRQ(ierr);
	ierr = KSPSetOperators(solver1,A,A); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(solver1); CHKERRQ(ierr);
	ierr = KSPSetUp(solver1); CHKERRQ(ierr);

	ierr = KSPSolve(solver1,F,T); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = m_field.set_global_ghost_vector("NORM_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

	/* Global error calculation */

    //PetscPrintf(PETSC_COMM_WORLD,"\n real part of l2 norm of analytical solution is: %f \n",real_analy);
    //PetscPrintf(PETSC_COMM_WORLD,"\n imag part of l2 norm of analytical solution is: %f \n",imag_analy);

  if(usel2) {
		double aa = sqrt(error)/sqrt(analy);
    PetscPrintf(PETSC_COMM_WORLD,"\n global l2 realtive error is: %f \n",aa);
  } else {

    double aa = sqrt(error)/sqrt(analy);
    PetscPrintf(PETSC_COMM_WORLD,"\n global H1 realtive error is: %f \n",aa);

  }

	/*  destroy objects  */
	ierr = MatDestroy(&A); CHKERRQ(ierr);
	ierr = VecDestroy(&F); CHKERRQ(ierr);
	ierr = VecDestroy(&T); CHKERRQ(ierr);
	ierr = KSPDestroy(&solver1); CHKERRQ(ierr);

	/* save local error indicator on mesh */
  PetscBool save_postproc_mesh = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL,"-save_postproc_mesh",&save_postproc_mesh,NULL); CHKERRQ(ierr);
  if(save_postproc_mesh) {

    PostPocOnRefinedMesh post_proc1(m_field);
		ierr = post_proc1.generateReferenceElementMesh(); CHKERRQ(ierr);
    ierr = post_proc1.addFieldValuesPostProc("erorNORM"); CHKERRQ(ierr);
    ierr = post_proc1.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.loop_finite_elements("NORM_PROBLEM","NORM_FE",post_proc1); CHKERRQ(ierr);
    rval = post_proc1.postProcMesh.write_file("norm_error.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);

  }

	ierr = PetscTime(&v2);CHKERRQ(ierr);
	ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f S CPU Time = %f S \n",pcomm->rank(),v2-v1,t2-t1);

	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;

}
