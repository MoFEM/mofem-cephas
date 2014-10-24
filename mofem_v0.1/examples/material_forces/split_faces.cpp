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


#include "ConfigurationalFractureMechanics.hpp"
#include "FieldCore.hpp"
#include "FaceSplittingTool.hpp"
#include "PostProcVertexMethod.hpp"


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
 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  BARRIER_RANK_START(pcomm) 
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  BARRIER_RANK_END(pcomm) 

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  FieldCore core(moab);
  FieldInterface& mField = core;

  ConfigurationalFractureMechanics conf_prob(mField);

  ierr = conf_prob.set_material_fire_wall(mField); CHKERRQ(ierr);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  ierr = mField.seed_ref_level_3D(0,BitRefLevel()); CHKERRQ(ierr);

  ierr = mField.build_fields(); CHKERRQ(ierr);
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out0.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }


        {
	  FaceSplittingTools face_splitting(mField);
	  ierr = main_refine_and_meshcat(mField,face_splitting); CHKERRQ(ierr);
	  ierr = face_splitting.cleanMeshsets(); CHKERRQ(ierr);
	}

	//project and set coords
	conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
	conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_material_positions) = 0;
	ierr = conf_prob.set_spatial_positions(mField); CHKERRQ(ierr);
	ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);

	//find faces for split
	FaceSplittingTools face_splitting(mField);
	ierr = main_select_faces_for_splitting(mField,face_splitting,10); CHKERRQ(ierr);
	//do splittig
	ierr = main_split_faces_and_update_field_and_elements(mField,face_splitting,10); CHKERRQ(ierr);
	
	//rebuild fields, finite elementa and problems
	ierr = main_face_splitting_restart(mField,conf_prob); CHKERRQ(ierr);

	//face spliting job done
	ierr = face_splitting.cleanMeshsets(); CHKERRQ(ierr);

	//solve spatial problem and calculate griffith forces

	//project and set coords
	conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
	conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_material_positions) = 0;
	ierr = conf_prob.set_spatial_positions(mField); CHKERRQ(ierr);
	ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);
 
	/*//debuging
	ierr = mField.partition_check_matrix_fill_in("MESH_SMOOTHING_PROBLEM",1); CHKERRQ(ierr);
	ierr = mField.partition_check_matrix_fill_in("ELASTIC_MECHANICS",1); CHKERRQ(ierr);
	ierr = mField.partition_check_matrix_fill_in("COUPLED_PROBLEM",1); CHKERRQ(ierr);*/

	SNES snes;
	//solve mesh smoothing
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);  
	ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
	SNESLineSearch linesearch;
	ierr = SNESGetLineSearch(snes,&linesearch); CHKERRQ(ierr);
	bool use_l2_instead_of_bt = true;
	if(use_l2_instead_of_bt) {
	  ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHL2); CHKERRQ(ierr);
	} else {
	  ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHBT); CHKERRQ(ierr);
	}
	Vec D_tmp_mesh_positions;
	ierr = mField.VecCreateGhost("MESH_SMOOTHING_PROBLEM",COL,&D_tmp_mesh_positions); CHKERRQ(ierr);
	ierr = mField.set_local_VecCreateGhost("MESH_SMOOTHING_PROBLEM",COL,D_tmp_mesh_positions,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = conf_prob.front_projection_data(mField,"MESH_SMOOTHING_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.surface_projection_data(mField,"MESH_SMOOTHING_PROBLEM"); CHKERRQ(ierr);
	Range tets;
	ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);

	struct CheckForNegatieVolume {
	  static PetscErrorCode F(FieldInterface &mField,const Range &tets,bool &flg) {
	    PetscFunctionBegin;
	    ErrorCode rval;
	    PetscErrorCode ierr;
	    double diffNTET[12],coords[12],V;
	    ierr = ShapeDiffMBTET(diffNTET); CHKERRQ(ierr);
	    Range::iterator tit = tets.begin();
	    for(;tit!=tets.end();tit++) {
	      int num_nodes;
	      const EntityHandle *conn;
	      rval = mField.get_moab().get_connectivity(*tit,conn,num_nodes,true); CHKERR_PETSC(rval);
	      int count = 0;
	      ierr = mField.get_FielData("MESH_NODE_POSITIONS",conn,num_nodes,coords,&count); CHKERRQ(ierr);
	      if(count != 3*num_nodes) {
		SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency %d != %d",count,3*num_nodes);
	      }
	      V = Shape_intVolumeMBTET(diffNTET,coords); 
	      if(V<=0) break;	  
	    }
	    if(tit==tets.end()) {
	      flg = true;
	    } else {
	      flg = false;
	    }
	    PetscFunctionReturn(0);
	  }
	};

	bool do_not_project = true;
	int nb_sub_steps = 1;
	int nn;
	do { 

	  nn = 1;
	  for(;nn<=nb_sub_steps;nn++) {

	    ierr = PetscPrintf(PETSC_COMM_WORLD,
	      "Mesh projection substep = %d out of %d (do not project %d)\n",
	      nn,nb_sub_steps,(int)do_not_project); CHKERRQ(ierr);
	    double alpha = fmin((double)nn/(double)nb_sub_steps,1);
	    ierr = face_splitting.calculateDistanceCrackFrontNodesFromCrackSurface(alpha); CHKERRQ(ierr);
	    //project nodes on crack surface
	    ierr = conf_prob.project_form_th_projection_tag(mField,"MESH_SMOOTHING_PROBLEM",do_not_project); CHKERRQ(ierr);
	    bool flg;
	    ierr = CheckForNegatieVolume::F(mField,tets,flg); CHKERRQ(ierr);
	    if(!flg) {
	      ierr = mField.set_global_VecCreateGhost(
		"MESH_SMOOTHING_PROBLEM",
		COL,D_tmp_mesh_positions,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	      nb_sub_steps++;
	      nn = 1;
	      break;
	    } else {
	      ierr = conf_prob.solve_mesh_smooting_problem(mField,&snes);  CHKERRQ(ierr);
	      SNESConvergedReason reason;
	      ierr = SNESGetConvergedReason(snes,&reason); CHKERRQ(ierr);
	      ierr = CheckForNegatieVolume::F(mField,tets,flg); CHKERRQ(ierr);
	      if(reason < 0 || !flg) {
		ierr = mField.set_global_VecCreateGhost(
		  "MESH_SMOOTHING_PROBLEM",
		  COL,D_tmp_mesh_positions,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
		if(use_l2_instead_of_bt) {
		  ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHBT); CHKERRQ(ierr);
		  use_l2_instead_of_bt = false;
		  nn--;
		} else {	  
		  ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHL2); CHKERRQ(ierr);
		  use_l2_instead_of_bt = true;
		  nb_sub_steps++;
		  nn = 1;
		  break;
		}
	      } else {
		ierr = mField.set_local_VecCreateGhost(
		  "MESH_SMOOTHING_PROBLEM",
		  COL,D_tmp_mesh_positions,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
		ierr = conf_prob.set_coordinates_from_material_solution(mField); CHKERRQ(ierr);
		if(nn == nb_sub_steps && do_not_project) {
		  do_not_project = false;
		  nb_sub_steps = 1;
		  nn = 1;
		  break;
		}
	      }
	    }

	  }

	  if(nb_sub_steps >= 10) break;

	} while(nn-1 != nb_sub_steps);
	ierr = conf_prob.set_coordinates_from_material_solution(mField); CHKERRQ(ierr);
	ierr = VecDestroy(&D_tmp_mesh_positions); CHKERRQ(ierr);
	ierr = SNESDestroy(&snes); CHKERRQ(ierr);

	ierr = conf_prob.delete_surface_projection_data(mField); CHKERRQ(ierr);
	ierr = conf_prob.delete_front_projection_data(mField); CHKERRQ(ierr);

	//save on mesh spatial positions
	conf_prob.material_FirelWall->
	  operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
	ierr = conf_prob.set_spatial_positions(mField); CHKERRQ(ierr);


	//solve spatial problem
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);  
	ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
	ierr = conf_prob.solve_spatial_problem(mField,&snes,false); CHKERRQ(ierr);
	ierr = SNESDestroy(&snes); CHKERRQ(ierr);

	//calculate Griffth forces
	ierr = conf_prob.calculate_material_forces(mField,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
	ierr = conf_prob.front_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.surface_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.project_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.griffith_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.griffith_g(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.delete_surface_projection_data(mField); CHKERRQ(ierr);
	ierr = conf_prob.delete_front_projection_data(mField); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_material(mField.get_moab(),"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("COUPLED_PROBLEM","MESH_NODE_POSITIONS",COL,ent_method_material); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("COUPLED_PROBLEM","MATERIAL_COUPLED",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  return 0;
}
