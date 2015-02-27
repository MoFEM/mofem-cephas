/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
 */

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 * PhD student Thomas Xuan Meng
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. 
 * abs(result) = sqrt(reEX^2+imEX^2) */


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
	bool userela; //relative error or pure error
	PetscBool flg = PETSC_TRUE;
	
	char mesh_file_name[255];
	ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) {
		SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
	}
	
	//char mesh_file_name2[255];
	//ierr = PetscOptionsGetString(PETSC_NULL,"-my_file2",mesh_file_name2,255,&flg); CHKERRQ(ierr);
	//if(flg != PETSC_TRUE) {
	//	SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file2 (MESH FILE NEEDED)");
	//}
	
	char type_error_norm[255];
	ierr = PetscOptionsGetString(PETSC_NULL,"-norm_type",type_error_norm,255,&flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) {
		SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -type_error_norm (L2 or H1 type needed)");
	}

	if (strcmp ("l2",type_error_norm ) == 0) {usel2 = true;}
	else if(strcmp ("h1",type_error_norm ) == 0) {usel2 = false;}
	
	char relative_error[255];
	ierr = PetscOptionsGetString(PETSC_NULL,"-relative_error",relative_error,255,&flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) {
		SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -type_error_norm (L2 or H1 type needed)");
	}
	
	if (strcmp ("true",relative_error ) == 0) {userela = true;}
	else if(strcmp ("false",relative_error ) == 0) {userela = false;}
	
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
	//ierr = m_field.add_field("relaNORM",H1,1); CHKERRQ(ierr);

	
	//Problem
	ierr = m_field.add_problem("NORM_PROBLEM"); CHKERRQ(ierr);

	//set refinment level for problem
	ierr = m_field.modify_problem_ref_level_add_bit("NORM_PROBLEM",bit_level0); CHKERRQ(ierr);

	//meshset consisting all entities in mesh
	EntityHandle root_set = moab.get_root_set(); 
	//add entities to field
	ierr = m_field.add_ents_to_field_by_TETs(root_set,"erorNORM"); CHKERRQ(ierr);
	//ierr = m_field.add_ents_to_field_by_TETs(root_set,"relaNORM"); CHKERRQ(ierr);
	
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
	//ierr = m_field.set_field_order(root_set,MBTET,"relaNORM",order); CHKERRQ(ierr);
	//ierr = m_field.set_field_order(root_set,MBTRI,"relaNORM",order); CHKERRQ(ierr);
	//ierr = m_field.set_field_order(root_set,MBEDGE,"relaNORM",order); CHKERRQ(ierr);
	//ierr = m_field.set_field_order(root_set,MBVERTEX,"relaNORM",1); CHKERRQ(ierr);
	
	if(!m_field.check_field("MESH_NODE_POSITIONS")) {
	ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
	ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
	ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
	}
	
	NormElement norm_elements(m_field);
	if(m_field.check_field("reEX") && m_field.check_field("rePRES") && m_field.check_field("imPRES") && m_field.check_field("imEX") ) {
	norm_elements.addNormElements("NORM_PROBLEM","NORM_FE","erorNORM","reEX","rePRES");}
	
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
	//print block sets with materials
	//ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);
	
	Vec F,D;
	ierr = m_field.VecCreateGhost("NORM_PROBLEM",ROW,&F); CHKERRQ(ierr);
	ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
	ierr = norm_elements.setNormFiniteElementRhsOperator("erorNORM","reEX","rePRES",F,usel2,userela); CHKERRQ(ierr);
	//ierr = norm_elements.setRelativeNormFiniteElementRhsOperator("relaNORM","reEX","rePRES",F,D,usel2); CHKERRQ(ierr);
    //Could we use 2 problem 1 element, two fields to calculate two result vectors ? 
	
	ierr = VecZeroEntries(F); CHKERRQ(ierr);
	//ierr = VecZeroEntries(D); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	//ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	//ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = m_field.set_global_VecCreateGhost("NORM_PROBLEM",ROW,F,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	//ierr = m_field.set_global_VecCreateGhost("NORM_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

	
	ierr = m_field.loop_finite_elements("NORM_PROBLEM","NORM_FE",norm_elements.getLoopFeRhs()); CHKERRQ(ierr);
	
	ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
	
	//ierr = VecGhostUpdateBegin(D,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	//ierr = VecGhostUpdateEnd(D,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	//ierr = VecAssemblyBegin(D); CHKERRQ(ierr);
	//ierr = VecAssemblyEnd(D); CHKERRQ(ierr);

	ierr = m_field.set_global_VecCreateGhost("NORM_PROBLEM",ROW,F,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	//ierr = m_field.set_global_VecCreateGhost("NORM_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

	
	/* Global error calculation */
	PetscReal l2norm,pointwisenorm;
	ierr = VecNorm(F,NORM_FROBENIUS,&l2norm);;
	ierr = VecNorm(F,NORM_MAX,&pointwisenorm);
    //ierr = VecMax(F,NULL,&pointwisenorm);
	
	//out stream the global error
	if(usel2) {
	std::cout << "\n The Global least square of l2 Norm of error is : --\n" << l2norm << std::endl;
	std::cout << "\n The Global Pointwise of l2 Norm of error for is : --\n" << pointwisenorm << std::endl;}
	else {
	std::cout << "\n The Global least square of H1 Norm of error is : --\n" << l2norm << std::endl;
	std::cout << "\n The Global Pointwise of H1 Norm of error for is : --\n" << pointwisenorm << std::endl;}
	/*    */
    ierr = VecDestroy(&F); CHKERRQ(ierr);
    //ierr = VecDestroy(&D); CHKERRQ(ierr);
	PostPocOnRefinedMesh post_proc1(m_field);
	ierr = post_proc1.generateRefereneElemenMesh(); CHKERRQ(ierr);
	ierr = post_proc1.addFieldValuesPostProc("erorNORM"); CHKERRQ(ierr);
	ierr = post_proc1.addFieldValuesGradientPostProc("erorNORM"); CHKERRQ(ierr);
	//ierr = post_proc1.addFieldValuesPostProc("relaNORM"); CHKERRQ(ierr);
	//ierr = post_proc1.addFieldValuesGradientPostProc("relaNORM"); CHKERRQ(ierr);
	ierr = post_proc1.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
	ierr = m_field.loop_finite_elements("NORM_PROBLEM","NORM_FE",post_proc1); CHKERRQ(ierr);
	rval = post_proc1.postProcMesh.write_file("norm_error.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);

	ierr = PetscTime(&v2);CHKERRQ(ierr);
	ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
	
	
	//output the results from Docker
	char command1[] = "mbconvert norm_error.h5m ./norm_error.vtk && cp ./norm_error.vtk ../../../../mnt/home/Desktop/U_pan/helmholtz_results/";
	int todo1 = system( command1 );
	
	
	
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;

}





