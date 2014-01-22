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

  ierr = conf_prob.constrains_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.material_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.coupled_problem_definition(mField); CHKERRQ(ierr);
  ierr = conf_prob.constrains_crack_front_problem_definition(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.arclength_problem_definition(mField); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC_COUPLED",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MATERIAL_COUPLED",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MATERIAL",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MESH_SMOOTHER",MBTET); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_set_bit("MATERIAL_MECHANICS",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("MATERIAL_MECHANICS_LAGRANGE_MULTIPLAIERS",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("COUPLED_PROBLEM",bit_level0); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_set_bit("CCT_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("C_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("C_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_set_bit("CTC_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = conf_prob.spatial_partition_problems(mField); CHKERRQ(ierr);
  ierr = conf_prob.material_partition_problems(mField); CHKERRQ(ierr);
  ierr = conf_prob.coupled_partition_problems(mField); CHKERRQ(ierr);
  ierr = conf_prob.constrains_partition_problems(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.crackfront_partition_problems(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);

  //caculate material forces
  ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);

  double da_0 = 0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_da",&da_0,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_da (what is the crack area increment ?)");
  }
  double da = da_0;

  int def_nb_load_steps = 0;
  Tag th_nb_load_steps;
  rval = mField.get_moab().tag_get_handle("_NB_LOAD_STEPS",1,MB_TYPE_INTEGER,
    th_nb_load_steps,MB_TAG_CREAT|MB_TAG_SPARSE,&def_nb_load_steps); 
  int *ptr_nb_load_steps;
  rval = mField.get_moab().tag_get_by_ptr(th_nb_load_steps,&root_meshset,1,(const void**)&ptr_nb_load_steps); CHKERR_PETSC(rval);
  int &step = *ptr_nb_load_steps;

  PetscInt nb_load_steps = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_load_steps",&nb_load_steps,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_load_steps (what is the number of load_steps ?)");
  }

  if(step == 0) {
    ierr = conf_prob.set_material_positions(mField); CHKERRQ(ierr);
  }

  Tag th_t_val;
  rval = moab.tag_get_handle("_LoadFactor_t_val",th_t_val); CHKERR_PETSC(rval);
  double *load_factor_ptr;
  rval = mField.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&load_factor_ptr); CHKERR_THROW(rval);
  double& load_factor = *load_factor_ptr;


  for(int aa = 0;step<nb_load_steps;step++,aa++) {

    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;
    ierr = PetscTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n** number of step = %D\n",step); CHKERRQ(ierr);

    ierr = conf_prob.set_coordinates_from_material_solution(mField); CHKERRQ(ierr);
    if(aa == 0) {
      ierr = conf_prob.calculate_material_forces(mField,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
      ierr = conf_prob.front_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.surface_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.project_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.griffith_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.griffith_g(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.delete_surface_projection_data(mField); CHKERRQ(ierr);
      ierr = conf_prob.delete_front_projection_data(mField); CHKERRQ(ierr);
    }

    double gc;
    PetscBool flg;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
    }

    //calulate initial load factor
    if(step == 0) {
      double max_j = conf_prob.max_j;
      double a = fabs(max_j)/pow(load_factor,2);
      double new_load_factor = copysign(sqrt(gc/a),load_factor);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\ncoefficient a = %6.4e\n",a); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"new load factor value = %6.4e\n\n",new_load_factor); CHKERRQ(ierr);
      load_factor = new_load_factor;
      SNES snes;
      //solve spatial problem
      ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);  
      ierr = conf_prob.solve_spatial_problem(mField,&snes,false); CHKERRQ(ierr);
      ierr = SNESDestroy(&snes); CHKERRQ(ierr);
      ierr = conf_prob.calculate_material_forces(mField,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
      ierr = conf_prob.front_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.surface_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.project_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.griffith_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.griffith_g(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD,"* da = %6.4e\n",da); CHKERRQ(ierr);
    SNES snes;
    ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
    
    Vec D0;
    ierr = mField.VecCreateGhost("COUPLED_PROBLEM",Col,&D0); CHKERRQ(ierr);
    ierr = mField.set_local_VecCreateGhost("COUPLED_PROBLEM",Col,D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    double load_factor0 = load_factor;
    
    int nb_sub_steps;
    ierr = PetscOptionsGetInt("","-my_nb_sub_steps",&nb_sub_steps,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      nb_sub_steps = 20;
    }

    bool first_step_converged = false;
    bool not_converged_state = false;
    double _da_ = 0;
    //int last_converged = -1;
    for(int ii = 0;ii<nb_sub_steps;ii++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n* number of substeps = %D\n\n",ii); CHKERRQ(ierr);
      if(ii == 0) {
	ierr = conf_prob.solve_coupled_problem(mField,&snes,(aa == 0) ? 0 : da); CHKERRQ(ierr);
      } else {
	ierr = conf_prob.solve_coupled_problem(mField,&snes,_da_); CHKERRQ(ierr);
      }
      int its;
      ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
      if(its == 0) break;
      SNESConvergedReason reason;
      ierr = SNESGetConvergedReason(snes,&reason); CHKERRQ(ierr);
      if(reason > 0) {
	if(aa > 0 && ii == 0) {
	  int its_d;
	  ierr = PetscOptionsGetInt("","-my_its_d",&its_d,&flg); CHKERRQ(ierr);
	  if(flg != PETSC_TRUE) {
	    its_d = 8;
	  }
	  double gamma = 0.5,reduction = 1;
	  reduction = pow((double)its_d/(double)(its+1),gamma);
	  const double max_da_reduction = 10;
	  if(reduction<1 || da < max_da_reduction*da_0) {
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n* change of da = %6.4e\n\n",reduction); CHKERRQ(ierr);
	    da = fmin(da*reduction,max_da_reduction*da_0);
	  }
	}
	_da_ = 0;
	//last_converged = ii;
	first_step_converged = true;
	not_converged_state = false;
	ierr = mField.set_local_VecCreateGhost("COUPLED_PROBLEM",Col,D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	load_factor0 = load_factor;
	ierr = conf_prob.calculate_material_forces(mField,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
	ierr = conf_prob.front_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.surface_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.project_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.griffith_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.griffith_g(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.delete_surface_projection_data(mField); CHKERRQ(ierr);
	ierr = conf_prob.delete_front_projection_data(mField); CHKERRQ(ierr);
      } else {
	not_converged_state = true;
	if(!first_step_converged) {
	  ierr = mField.set_global_VecCreateGhost("COUPLED_PROBLEM",Col,D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	  load_factor = load_factor0;
	  da = 0.5*da;
	  _da_ = da;
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"* failed to converge, set da = %6.4e ( 0.5 )\n",_da_); CHKERRQ(ierr);
	} else {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"* reset unknowns vector\n"); CHKERRQ(ierr);
	  ierr = mField.set_global_VecCreateGhost("COUPLED_PROBLEM",Col,D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	  load_factor = load_factor0;
	  ierr = conf_prob.set_coordinates_from_material_solution(mField); CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"* failed to converge, recalculate spatial positions only\n"); CHKERRQ(ierr);
	  SNES snes_spatial;
	  //solve spatial problem
	  ierr = SNESCreate(PETSC_COMM_WORLD,&snes_spatial); CHKERRQ(ierr);  
	  ierr = conf_prob.solve_spatial_problem(mField,&snes_spatial,false); CHKERRQ(ierr);
	  ierr = SNESDestroy(&snes_spatial); CHKERRQ(ierr);
	  ierr = conf_prob.calculate_material_forces(mField,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
	  ierr = conf_prob.front_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	  ierr = conf_prob.surface_projection_data(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	  ierr = conf_prob.project_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	  ierr = conf_prob.griffith_force_vector(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	  ierr = conf_prob.griffith_g(mField,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	  if(not_converged_state) {
	    _da_+= 0.1*da;
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"* failed to converge, set da = %6.4e ( 0.1 ) \n",_da_); CHKERRQ(ierr);
	  }
	}
      }
      //set line sercher L2 for not converged state
      if(reason < 0) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"* set L2 linesercher\n",_da_); CHKERRQ(ierr);
	//set line sercher L2 for not converged state
	ierr = SNESDestroy(&snes); CHKERRQ(ierr);
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
	ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
	SNESLineSearch linesearch;
	ierr = SNESGetLineSearch(snes,&linesearch); CHKERRQ(ierr);
	ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHL2); CHKERRQ(ierr);
      }
    }

    ierr = VecDestroy(&D0); CHKERRQ(ierr);
    ierr = SNESDestroy(&snes); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"load_path: %4D Area %6.4e Lambda %6.4e Energy %6.4e\n",step,conf_prob.aRea,conf_prob.lambda,conf_prob.energy); CHKERRQ(ierr);
    if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      ierr = mField.problem_get_FE("COUPLED_PROBLEM","MATERIAL_COUPLED",out_meshset); CHKERRQ(ierr);
      ostringstream ss;
      ss << "out_load_step_" << step << ".vtk";
      rval = moab.write_file(ss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);

      {
	Range SurfacesFaces;
	ierr = mField.get_Cubit_msId_entities_by_dimension(102,SideSet,2,SurfacesFaces,true); CHKERRQ(ierr);
	Range CrackSurfacesFaces;
	ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,CrackSurfacesFaces,true); CHKERRQ(ierr);
	for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
	  int msId = it->get_msId();
	  if((msId < 10200)||(msId >= 10300)) continue;
	  Range SurfacesFaces_msId;
	  ierr = mField.get_Cubit_msId_entities_by_dimension(msId,SideSet,2,SurfacesFaces_msId,true); CHKERRQ(ierr);
	  SurfacesFaces.insert(SurfacesFaces_msId.begin(),SurfacesFaces_msId.end());
	}

	Range level_tris;
	ierr = mField.refine_get_ents(bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
	SurfacesFaces = intersect(SurfacesFaces,level_tris);
	CrackSurfacesFaces = intersect(CrackSurfacesFaces,level_tris);

	EntityHandle out_meshset;
	rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	rval = moab.add_entities(out_meshset,CrackSurfacesFaces); CHKERR_PETSC(rval);
	ostringstream ss1;
	ss1 << "out_crack_surface_" << step << ".vtk";
	rval = moab.write_file(ss1.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	rval = moab.add_entities(out_meshset,SurfacesFaces); CHKERR_PETSC(rval);
	ostringstream ss2;
	ss2 << "out_surface_" << step << ".vtk";
	rval = moab.write_file(ss2.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      }

      ostringstream ss3;
      ss3 << "restart_" << step << ".h5m";
      rval = moab.write_file(ss3.str().c_str()); CHKERR_PETSC(rval);

      const MoFEMProblem *problem_ptr;
      ierr = mField.get_problem("COUPLED_PROBLEM",&problem_ptr); CHKERRQ(ierr);

      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|UnknownCubitName,it)) {
	if(it->get_Cubit_name() != "LoadPath") continue;

	Range nodes;
	rval = moab.get_entities_by_type(it->meshset,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
	for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {

	  double coords[3];
	  rval = moab.get_coords(&*nit,1,coords); CHKERR_PETSC(rval);
	  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problem_ptr,*nit,dof)) {
	    if(dof->get_name()!="SPATIAL_POSITION") continue;
	    ierr = PetscPrintf(PETSC_COMM_WORLD,
	      "load_path_disp ent %ld dim %d "
	      "coords ( %8.6f %8.6f %8.6f ) "
	      "val %6.4e Lambda %6.4e\n",
	      dof->get_ent(),dof->get_dof_rank(),
	      coords[0],coords[1],coords[2],
	      dof->get_FieldData()-coords[dof->get_dof_rank()],      
	      conf_prob.lambda); CHKERRQ(ierr);
	  }
	}

      }

    }

    ierr = PetscTime(&v2);CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Step Time Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);

  }

  if(pcomm->rank()==0) {
    rval = moab.write_file("out_arc_length.h5m"); CHKERR_PETSC(rval);
  }

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.refine_get_ents(bit_level0,BitRefLevel().set(),MBEDGE,out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out_edges.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

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
