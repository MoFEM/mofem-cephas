/* Copyright (C) 2014, Guoqiang Xue (GUOQIANG XUE AT glasgow.ac.uk)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
 */

/* This file is part of MoFEM.
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#include <MoFEM.hpp>

#include <DirichletBC.hpp>
#include <PotsProcOnRefMesh.hpp>

#include <Projection10NodeCoordsOnField.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

extern "C" {
#include <gm_rule.h>
}


#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <SurfacePressure.hpp>
#include <NodalForce.hpp>
#include <FluidPressure.hpp>
#include <BodyForce.hpp>
#include <ThermalStressElement.hpp>
#include <ThermalElement.hpp>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

using namespace boost::numeric;
using namespace ObosleteUsersModules;

#include "TDPotentialFlowFEMethod.hpp"


ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";


int main(int argc, char *argv[]) {
	
	try {

/*****/
		PetscInitialize(&argc,&argv,(char *)0,help);
		
		moab::Core mb_instance;
		Interface& moab = mb_instance;
		int rank;
		MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

		//Read parameters from line command
		PetscBool flg = PETSC_TRUE;
		char mesh_file_name[255];
		ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
		if(flg != PETSC_TRUE) {
			SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
		}
				
		ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
		if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
		
		//Read mesh to MOAB
		const char *option;
		option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
		BARRIER_RANK_START(pcomm) 
		rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
		BARRIER_RANK_END(pcomm) 
		
		//Create MoFEM database
		MoFEM::Core core(moab);
		FieldInterface& m_field = core;

		//Set order
		PetscInt order;
		ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
		if(flg != PETSC_TRUE) {
			order = 1;
		}	
		
		//ref meshset ref level 0
//		ierr = m_field.seed_ref_level_3D(0,0); CHKERRQ(ierr);
		
		// stl::bitset see for more details
		BitRefLevel bit_level0;
		bit_level0.set(0);
		EntityHandle meshset_level0;
		rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
		ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
//		ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

/*****/		
		/* define the problem */	
		
		//add fields
		ierr = m_field.add_field("2DPOTENTIAL_FIELD",H1,1); CHKERRQ(ierr);
		
		//Problem
		ierr = m_field.add_problem("2DPOTENTIAL_PROBLEM"); CHKERRQ(ierr);

		//set refinment level for problem
		ierr = m_field.modify_problem_ref_level_add_bit("2DPOTENTIAL_PROBLEM",bit_level0); CHKERRQ(ierr);
		
		Range SurTris;
		ierr = m_field.get_Cubit_msId_entities_by_dimension(101,SIDESET,2,SurTris,true); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET 101 = %d\n",SurTris.size()); CHKERRQ(ierr);

		//add entitities (by tris) to the field
		
		EntityHandle BoundFacesMeshset;
		rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
		rval = moab.add_entities(BoundFacesMeshset,SurTris); CHKERR_PETSC(rval);
		ierr = m_field.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);
		ierr = m_field.add_ents_to_field_by_TRIs(BoundFacesMeshset,"2DPOTENTIAL_FIELD",2); CHKERRQ(ierr); 
		
/*****/
		/* declare the problem */		
		
    	ierr = m_field.set_field_order(0,MBVERTEX,"2DPOTENTIAL_FIELD",1); CHKERRQ(ierr);
		ierr = m_field.set_field_order(0,MBEDGE,"2DPOTENTIAL_FIELD",order); CHKERRQ(ierr);
		ierr = m_field.set_field_order(0,MBTRI,"2DPOTENTIAL_FIELD",order); CHKERRQ(ierr);
		
     
        TDPotentialElement TDPotential_elements(m_field);
		ierr = TDPotential_elements.add2DPotential_Element("2DPOTENTIAL_PROBLEM","2DPOTENTIAL_FIELD"); CHKERRQ(ierr);
		
	
/*****/		
		/* build database */
		//build fields
		ierr = m_field.build_fields(); CHKERRQ(ierr);
		//build finite elements
		ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
		//build adjacencies
		ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
		//build problem
		ierr = m_field.build_problems(); CHKERRQ(ierr);
/*****/
		/* mesh partitioning */
		
		//partition problems
		ierr = m_field.partition_problem("2DPOTENTIAL_PROBLEM"); CHKERRQ(ierr);
		ierr = m_field.partition_finite_elements("2DPOTENTIAL_PROBLEM"); CHKERRQ(ierr);
		ierr = m_field.partition_ghost_dofs("2DPOTENTIAL_PROBLEM"); CHKERRQ(ierr);
		
		//print bcs
//		ierr = m_field.print_cubit_displacement_set(); CHKERRQ(ierr);
//		ierr = m_field.print_cubit_pressure_set(); CHKERRQ(ierr);

/*****/
		//create matrices and vectors
		Vec F,D;
		ierr = m_field.VecCreateGhost("2DPOTENTIAL_PROBLEM",ROW,&F); CHKERRQ(ierr);
		ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
		Mat A;
		ierr = m_field.MatCreateMPIAIJWithArrays("2DPOTENTIAL_PROBLEM",&A); CHKERRQ(ierr);
		
//		Projection10NodeCoordsOnField ent_method_material(mField,"MESH_NODE_POSITIONS");
//		ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
		
		//0-----------------0//
		//boundary conditions
		//0-----------------0//
		
		ierr = TDPotential_elements.setTDPotentialFiniteElementLhsOperators("2DPOTENTIAL_FIELD",A); CHKERRQ(ierr);
//		ierr = TDPotential_elements.setTDPotentialFiniteElementRhsOperators("TEMP",F); CHKERRQ(ierr);

//		ierr = VecZeroEntries(D); CHKERRQ(ierr);
//		ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//		ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//		ierr = VecZeroEntries(F); CHKERRQ(ierr);
//		ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//		ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    	ierr = MatZeroEntries(A); CHKERRQ(ierr);
		ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
		ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	    ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);//PETSC_VIEWER_DRAW_WORLD);
		cout << "3333333333333333"<<endl;
//		//ierr = m_field.loop_finite_elements("THERMAL_PROBLEM","THERMAL_CONVECTION_FE",thermal_elements.getLoopFeConvectionRhs()); CHKERRQ(ierr);
		ierr = m_field.loop_finite_elements("2DPOTENTIAL_PROBLEM","2DPOTENTIAL_ELEM",TDPotential_elements.getLoopFe2DpoetentialLhs()); CHKERRQ(ierr);
		ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
		ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
		ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);

		
		//set matrix possitives define and symetric for cholesky and icc preceonditionser
/*		ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
		
		//Solver
		KSP solver;
		ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
		ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
		ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
		ierr = KSPSetUp(solver); CHKERRQ(ierr);
		
		ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
		ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
		ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
		ierr = mField.set_global_VecCreateGhost("POTENTIAL_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
		
		ierr = KSPDestroy(&solver); CHKERRQ(ierr);
		ierr = VecDestroy(&F); CHKERRQ(ierr);
		ierr = VecDestroy(&D); CHKERRQ(ierr);
		ierr = MatDestroy(&A); CHKERRQ(ierr);
		
		/*Tag th_phi;
		 double def_val = 0;
		 rval = moab.tag_get_handle("PHI",1,MB_TYPE_DOUBLE,th_phi,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERR_PETSC(rval);
		 for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"POTENTIAL_FIELD",dof)) {
		 EntityHandle ent = dof->get_ent();
		 double val = dof->get_FieldData();
		 rval = moab.tag_set_data(th_phi,&ent,1,&val); CHKERR_PETSC(rval);
		 }*/
/*		
		ProjectionFieldOn10NodeTet ent_method_phi_on_10nodeTet(mField,"POTENTIAL_FIELD",true,false,"PHI");
		ierr = mField.loop_dofs("POTENTIAL_FIELD",ent_method_phi_on_10nodeTet); CHKERRQ(ierr);
		ent_method_phi_on_10nodeTet.set_nodes = false;
		ierr = mField.loop_dofs("POTENTIAL_FIELD",ent_method_phi_on_10nodeTet); CHKERRQ(ierr);
		
		if(pcomm->rank()==0) {
			rval = moab.write_file("solution.h5m"); CHKERR_PETSC(rval);
		}
		
		if(pcomm->rank()==0) {
			EntityHandle out_meshset;
			rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
			ierr = mField.problem_get_FE("POTENTIAL_PROBLEM","POTENTIAL_ELEM",out_meshset); CHKERRQ(ierr);
			rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
			rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
		}
		
		PostProcScalarFieldsAndGradientOnRefMesh fe_post_proc_method(moab,"POTENTIAL_FIELD");
		ierr = mField.loop_finite_elements("POTENTIAL_PROBLEM","POTENTIAL_ELEM",fe_post_proc_method);  CHKERRQ(ierr);
		
		if(pcomm->rank()==0) {
			rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
		}
*/		
		ierr = MatDestroy(&A); CHKERRQ(ierr);
		ierr = VecDestroy(&F); CHKERRQ(ierr);
		ierr = VecDestroy(&D); CHKERRQ(ierr);
//		ierr = KSPDestroy(&solver); CHKERRQ(ierr);
		
		PetscFinalize();
		
	} catch (const char* msg) {
		SETERRQ(PETSC_COMM_SELF,1,msg);
	} catch (const std::exception& ex) {
		ostringstream ss;
		ss << "thorw in method: " << ex.what() << endl;
		SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	}
	
}
