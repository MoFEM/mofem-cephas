/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "FEMethod_UpLevelStudent.hpp"
#include "PotentialFlowFEMethod.hpp"
#include "FunctionsForFieldData.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";


int main(int argc, char *argv[]) {
    
    try {
        
        PetscInitialize(&argc,&argv,(char *)0,help);
        
        Core mb_instance;
        Interface& moab = mb_instance;
        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
        
        PetscBool flg = PETSC_TRUE;
        char mesh_file_name[255];
        ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
        if(flg != PETSC_TRUE) {
            SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
        }
        PetscInt order;
        ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
        if(flg != PETSC_TRUE) {
            order = 1;
        }
		
        PetscInt mesh_refinement_level;
        ierr = PetscOptionsGetInt(PETSC_NULL,"-my_mesh_ref_level",&mesh_refinement_level,&flg); CHKERRQ(ierr);
        if(flg != PETSC_TRUE) {
            mesh_refinement_level = 0;
        }
		
        
        const char *option;
        option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
        rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval);
        ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
        if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
        
		
        FieldCore core(moab);
        FieldInterface& mField = core;
        
        //add fields
        ierr = mField.add_field("POTENTIAL_FIELD",H1,1); CHKERRQ(ierr);
        
        ///Getting No. of Fibres to be used for Potential Flow Problem
        int noOfFibres=0;
        for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|UnknownCubitName,it)) {
            
            std::size_t found=it->get_Cubit_name().find("PotentialFlow");
            if (found==std::string::npos) continue;
            noOfFibres += 1;
        }
		cout<<"No. of Fibres for Potential Flow : "<<noOfFibres<<endl;
        
		vector<int> fibreList(noOfFibres,0);  //define a vector of type integer of size noOfFibres and initialize with 0
		for (int aa=0; aa<noOfFibres; aa++) {
			fibreList[aa] = aa + 1;
		}
        
		for (int cc = 0; cc < noOfFibres; cc++) {
            
            ostringstream sss, rrr;
            
            //add finite elements
            sss << "POTENTIAL_ELEM" << fibreList[cc];
            cout<<sss.str().c_str()<<endl;
            ierr = mField.add_finite_element( sss.str().c_str() ); CHKERRQ(ierr);
            ierr = mField.modify_finite_element_add_field_row( sss.str().c_str() ,"POTENTIAL_FIELD"); CHKERRQ(ierr);
            ierr = mField.modify_finite_element_add_field_col( sss.str().c_str() ,"POTENTIAL_FIELD"); CHKERRQ(ierr);
            ierr = mField.modify_finite_element_add_field_data( sss.str().c_str() ,"POTENTIAL_FIELD"); CHKERRQ(ierr);
            //add problems
            rrr << "POTENTIAL_PROBLEM" << fibreList[cc];
            ierr = mField.add_problem( rrr.str().c_str() ); CHKERRQ(ierr);
            //define problems and finite elements
            ierr = mField.modify_problem_add_finite_element( rrr.str().c_str() , sss.str().c_str() ); CHKERRQ(ierr);
            
        }
        
		Tag th_meshset_info;
		int def_meshset_info[2] = {0,0};
		rval = moab.tag_get_handle("MESHSET_INFO",2,MB_TYPE_INTEGER,th_meshset_info,MB_TAG_CREAT|MB_TAG_SPARSE,&def_meshset_info);
		
		int meshset_data[2];
        EntityHandle root = moab.get_root_set();
		rval = moab.tag_get_data(th_meshset_info,&root,1,meshset_data); CHKERR_PETSC(rval);
		
        ierr = mField.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);
        vector<BitRefLevel> bit_levels;
        bit_levels.push_back(BitRefLevel().set(0));
		
		//MESH REFINEMENT
		int ll = 1;
        
		//*****INTERFACE INSERTION******
//        for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,cit)) {
        for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|InterfaceSet,cit)) {
//            std::size_t interfaceFound=cit->get_Cubit_name().find("interface");
//            if (interfaceFound==std::string::npos) continue;
            
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Insert Interface %d\n",cit->get_msId()); CHKERRQ(ierr);
            EntityHandle cubit_meshset = cit->get_meshset();
            int meshset_data_int[2];
            rval = moab.tag_get_data(th_meshset_info,&cubit_meshset,1,meshset_data_int); CHKERR_PETSC(rval);
            if(meshset_data_int[1]==0){
                //get tet enties form back bit_level
                EntityHandle ref_level_meshset = 0;
                rval = moab.create_meshset(MESHSET_SET,ref_level_meshset); CHKERR_PETSC(rval);
                ierr = mField.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBTET,ref_level_meshset); CHKERRQ(ierr);
                ierr = mField.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBPRISM,ref_level_meshset); CHKERRQ(ierr);
                Range ref_level_tets;
                rval = moab.get_entities_by_handle(ref_level_meshset,ref_level_tets,true); CHKERR_PETSC(rval);
                //get faces and test to split
                ierr = mField.get_msId_3dENTS_sides(cubit_meshset,bit_levels.back(),true); CHKERRQ(ierr);
                //set new bit level
                bit_levels.push_back(BitRefLevel().set(ll++));
                //split faces and edges
                ierr = mField.get_msId_3dENTS_split_sides(ref_level_meshset,bit_levels.back(),cubit_meshset,true,true,0); CHKERRQ(ierr);
                //clean meshsets
                rval = moab.delete_entities(&ref_level_meshset,1); CHKERR_PETSC(rval);
                int meshset_data_ins[2] = {ll,1};
                rval = moab.tag_set_data(th_meshset_info,&cubit_meshset,1,meshset_data_ins); CHKERR_PETSC(rval);
            }
            //update cubit meshsets
            for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,ciit)) {
                EntityHandle cubit_meshset = ciit->meshset;
                ierr = mField.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
                ierr = mField.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
                ierr = mField.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBTRI,true); CHKERRQ(ierr);
                ierr = mField.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBTET,true); CHKERRQ(ierr);
            }
        }
		
		//******REFINE MESH*****//
		for (int ref_lev=1; ref_lev<mesh_refinement_level+1; ref_lev++ ) {
			EntityHandle refine_meshset;
			rval = moab.create_meshset(MESHSET_SET,refine_meshset); CHKERR_PETSC(rval);
			ierr = mField.get_entities_by_ref_level(bit_levels.back(),BitRefLevel().set(),refine_meshset); CHKERRQ(ierr);
			
			bit_levels.push_back(BitRefLevel().set(ll++));
			
			ierr = mField.add_verices_in_the_middel_of_edges(refine_meshset,bit_levels.back()); CHKERRQ(ierr);
			ierr = mField.refine_TET(refine_meshset,bit_levels.back(),true); CHKERRQ(ierr);
			ierr = mField.refine_PRISM(refine_meshset,bit_levels.back(),0); CHKERRQ(ierr);
            
			for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,ciit)) {
				EntityHandle cubit_meshset = ciit->meshset;
				ierr = mField.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
				ierr = mField.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
				ierr = mField.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBTRI,true); CHKERRQ(ierr);
				ierr = mField.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBTET,true); CHKERRQ(ierr);
			}
		}
		
		//End of refinement, save level of refinement
		int meshset_data_root[2]={ll,0};
		rval = moab.tag_set_data(th_meshset_info,&root,1,meshset_data_root); CHKERR_PETSC(rval);
		
		/******************TETS TO MESHSET AND SAVING TETS ENTITIES******************/
		EntityHandle out_tet_meshset;
		rval = moab.create_meshset(MESHSET_SET,out_tet_meshset); CHKERR_PETSC(rval);
		ierr = mField.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBTET,out_tet_meshset); CHKERRQ(ierr);
		rval = moab.write_file("out_tets.vtk","VTK","",&out_tet_meshset,1); CHKERR_PETSC(rval);
		/*******************************************************/
		
		/******************PRISMS TO MESHSET AND SAVING PRISMS ENTITIES******************/
		EntityHandle out_meshset_prism;
		rval = moab.create_meshset(MESHSET_SET,out_meshset_prism); CHKERR_PETSC(rval);
		ierr = mField.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBPRISM,out_meshset_prism); CHKERRQ(ierr);
		rval = moab.write_file("out_prism.vtk","VTK","",&out_meshset_prism,1); CHKERR_PETSC(rval);
		/*******************************************************/
		
		EntityHandle out_meshset;
		rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
		ierr = mField.get_entities_by_ref_level(bit_levels.back(),BitRefLevel().set(),out_meshset); CHKERRQ(ierr);
		rval = moab.write_file("out_all_mesh.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
		Range LatestRefinedTets;
		rval = moab.get_entities_by_type(out_meshset, MBTET,LatestRefinedTets,true); CHKERR_PETSC(rval);
        
		
        BitRefLevel problem_bit_level = bit_levels.back();
        
        ///Adding entities to Field and FE for Potential Flow Problem
        for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|UnknownCubitName,it)) {
            
            //      std::size_t found=it->get_Cubit_name().find("PotentialFlow");
            //      if (found==std::string::npos) continue;
            //          cout<<it->get_Cubit_name()<<endl;
            
            for (int cc = 0; cc < noOfFibres; cc++) {
                ostringstream sss,rrr;
                //set problem level
                sss << "POTENTIAL_ELEM" << fibreList[cc];
                rrr << "PotentialFlow" << fibreList[cc];
                
                if(it->get_Cubit_name() ==  rrr.str().c_str() ) {
                    Range TetsInBlock;
                    rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
                    Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);
                    EntityHandle block_meshset;
                    rval = moab.create_meshset(MESHSET_SET,block_meshset); CHKERR_PETSC(rval);
                    rval = moab.add_entities(block_meshset,block_rope_bit_level);CHKERR_PETSC(rval);
                    
                    ierr = mField.add_ents_to_field_by_TETs(0,"POTENTIAL_FIELD"); CHKERRQ(ierr);
                    //add finite elements entities
                    ierr = mField.add_ents_to_finite_element_by_TETs(block_meshset, sss.str().c_str() ,true); CHKERRQ(ierr);
                    rval = moab.delete_entities(&block_meshset,1); CHKERR_PETSC(rval);
                }
            }
        }
        
        ierr = mField.set_field_order(0,MBVERTEX,"POTENTIAL_FIELD",1); CHKERRQ(ierr);
        ierr = mField.set_field_order(0,MBEDGE,"POTENTIAL_FIELD",order); CHKERRQ(ierr);
        ierr = mField.set_field_order(0,MBTRI,"POTENTIAL_FIELD",order); CHKERRQ(ierr);
        ierr = mField.set_field_order(0,MBTET,"POTENTIAL_FIELD",order); CHKERRQ(ierr);
        
        for (int cc = 0; cc < noOfFibres; cc++) {
            ostringstream sss;
            //set problem level
            sss << "POTENTIAL_PROBLEM" << fibreList[cc];
            ierr = mField.modify_problem_ref_level_add_bit( sss.str().c_str() ,problem_bit_level); CHKERRQ(ierr);
        }
        
        //build fields
        ierr = mField.build_fields(); CHKERRQ(ierr);
        //build finite elements
        ierr = mField.build_finite_elements(); CHKERRQ(ierr);
        //build adjacencies
        ierr = mField.build_adjacencies(problem_bit_level); CHKERRQ(ierr);
        //build problem
        ierr = mField.build_problems(); CHKERRQ(ierr);
        
        //partition problems
        for (int cc = 0; cc < noOfFibres; cc++) {
            ostringstream sss;
            sss << "POTENTIAL_PROBLEM" << fibreList[cc];
            ierr = mField.partition_problem( sss.str().c_str() ); CHKERRQ(ierr);
            ierr = mField.partition_finite_elements( sss.str().c_str() ); CHKERRQ(ierr);
            ierr = mField.partition_ghost_dofs( sss.str().c_str() ); CHKERRQ(ierr);
        }
        
        //print bcs
        ierr = mField.printCubitDisplacementSet(); CHKERRQ(ierr);
        ierr = mField.printCubitPressureSet(); CHKERRQ(ierr);
        
        //create matrices and vectors
		vector<Vec> F(noOfFibres);
		vector<Vec> D(noOfFibres);
		vector<Mat> A(noOfFibres);
		
        for (int cc = 0; cc < noOfFibres; cc++) {
            ostringstream sss,rrr;
            sss << "POTENTIAL_PROBLEM" << fibreList[cc];
            rrr << "POTENTIAL_ELEM" << fibreList[cc];
            ierr = mField.VecCreateGhost( sss.str().c_str() ,Row,&F[cc]); CHKERRQ(ierr);
            ierr = mField.VecCreateGhost( sss.str().c_str() ,Col,&D[cc]); CHKERRQ(ierr);
            ierr = mField.MatCreateMPIAIJWithArrays( sss.str().c_str() ,&A[cc]); CHKERRQ(ierr);
            
            LaplacianElem elem(mField,A[cc],F[cc]);
            
            ierr = MatZeroEntries(A[cc]); CHKERRQ(ierr);
            ierr = mField.loop_finite_elements( sss.str().c_str() , rrr.str().c_str() ,elem);  CHKERRQ(ierr);
            
            ierr = MatAssemblyBegin(A[cc],MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A[cc],MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
            
            ierr = VecAssemblyBegin(F[cc]); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(F[cc]); CHKERRQ(ierr);
            
            //		VecView(F[0],PETSC_VIEWER_STDOUT_WORLD);
            
            KSP solver;
            ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
            ierr = KSPSetOperators(solver,A[cc],A[cc],SAME_NONZERO_PATTERN); CHKERRQ(ierr);
            ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
            ierr = KSPSetUp(solver); CHKERRQ(ierr);
            
            ierr = KSPSolve(solver,F[cc],D[cc]); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(D[cc],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(D[cc],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = mField.set_global_VecCreateGhost( sss.str().c_str() ,Row,D[cc],INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            
        }
        
        Tag th_phi;
        double def_val = 0;
        rval = moab.tag_get_handle("PHI",1,MB_TYPE_DOUBLE,th_phi,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERR_PETSC(rval);
        for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"POTENTIAL_FIELD",dof)) {
            EntityHandle ent = dof->get_ent();
            double val = dof->get_FieldData();
            rval = moab.tag_set_data(th_phi,&ent,1,&val); CHKERR_PETSC(rval);
        }
        
        
        if(pcomm->rank()==0) {
            rval = moab.write_file("solution1.h5m"); CHKERR_PETSC(rval);
        }
        
		EntityHandle out_meshset1;
		rval = moab.create_meshset(MESHSET_SET,out_meshset1); CHKERR_PETSC(rval);
		ierr = mField.get_entities_by_type_and_ref_level(bit_levels[0],BitRefLevel().set(),MBTET,out_meshset1); CHKERRQ(ierr);
		rval = moab.write_file("solution2.vtk","VTK","",&out_meshset1,1); CHKERR_PETSC(rval);
		
        if(pcomm->rank()==0) {
            EntityHandle out_meshset;
            rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
            
            for (int cc = 0; cc < noOfFibres; cc++) {
                EntityHandle out_meshset1;
                rval = moab.create_meshset(MESHSET_SET,out_meshset1); CHKERR_PETSC(rval);
                ostringstream sss,rrr,ttt;
                sss << "POTENTIAL_PROBLEM" << fibreList[cc];
                rrr << "POTENTIAL_ELEM" << fibreList[cc];
                ttt << "out_potential_flow" << fibreList[cc] <<".vtk";
                ierr = mField.problem_get_FE( sss.str().c_str() , rrr.str().c_str() ,out_meshset); CHKERRQ(ierr);
                ierr = mField.problem_get_FE( sss.str().c_str() , rrr.str().c_str() ,out_meshset1); CHKERRQ(ierr);
                
                rval = moab.write_file( ttt.str().c_str() ,"VTK","",&out_meshset1,1); CHKERR_PETSC(rval);
                rval = moab.delete_entities(&out_meshset1,1); CHKERR_PETSC(rval);
            }
            
            rval = moab.write_file("out_potential_flow.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
            rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
        }
        
        PetscFinalize();
        
    } catch (const char* msg) {
        SETERRQ(PETSC_COMM_SELF,1,msg);
    } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    
}
