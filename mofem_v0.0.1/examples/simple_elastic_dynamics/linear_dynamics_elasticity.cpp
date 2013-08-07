/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * Test for linar elastic dynamics.
 *
 * This is not exactly procedure for linear elatic dynamics, since jacobian is
 * evaluated at every time step and snes procedure is involved. However it is
 * implemented like that, to test methodology for general nonlinear problem.
 *
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

#include "linear_dynamics_elasticity.hpp"

using namespace MoFEM;

static char help[] = "...\n\n";

struct MyBC: public DynamicElasticFEMethod::BC {

    MyBC(): DynamicElasticFEMethod::BC()  {}
    PetscErrorCode f_CalcTraction(double ts_t,ublas::vector<double,ublas::bounded_array<double,3> >& traction) const {
      PetscFunctionBegin;
      
      //Triangular loading over 10s (maximum at 5)
      double scale = 0;
      if(ts_t < 5.) scale = ts_t/5.;
      if(ts_t > 5.) scale = 1.+(5.-ts_t)/5.;
      if(ts_t > 10.) scale = 0;

      //Set Direction of Traction On SideSet2
      traction[0] = 0; //X
      traction[1] = 0; //Y 
      traction[2] = scale; //Z*/

      PetscFunctionReturn(0);
    }

    MyBC(double _final_time): DynamicElasticFEMethod::BC(_final_time) {};
    PetscErrorCode f_CalcDisp(double ts_t,
      ublas::vector<double,ublas::bounded_array<double,3> >& disp,
      ublas::vector<double,ublas::bounded_array<double,3> >& vel) const {

      PetscFunctionBegin;

      double disp_val = 0;
      double vel_val = 0;
      const double v1 = 0.2;
      const double init_t = 1.;
      if(ts_t < init_t) {
	  disp_val = v1 * ( (ts_t*ts_t/2.)/init_t );
	  vel_val = v1* ( ts_t/init_t );
      } else if(ts_t < final_time) {
	  disp_val = v1*( ((init_t*init_t/2.)/init_t) + (ts_t-init_t) );
	  vel_val = v1;
      } 

      disp[0] = 0;
      vel[0] = 0;
      disp[1] = 0;
      vel[1] = 0;
      disp[2] = disp_val;
      vel[2] = vel_val;

      PetscFunctionReturn(0);
    }

};


int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //Reade parameters from line command
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //We need that for code profiling
  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  //Create MoFEM (Joseph) database
  moabField_Core core(moab);
  moabField& mField = core;

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

  /***/
  //Define problem

  //Fields
  ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("VELOCITIES",L2,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("STIFFNESS"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("MASS"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("COPUPLING_VV"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("COPUPLING_VU"); CHKERRQ(ierr);


  //Define rows/cols and element data
  //STIFFNESS
  ierr = mField.modify_finite_element_add_field_row("STIFFNESS","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("STIFFNESS","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("STIFFNESS","DISPLACEMENT"); CHKERRQ(ierr);

  //MASS
  ierr = mField.modify_finite_element_add_field_row("MASS","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MASS","VELOCITIES"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MASS","DISPLACEMENT"); CHKERRQ(ierr);

  //COPUPLING
  //VV
  ierr = mField.modify_finite_element_add_field_row("COPUPLING_VV","VELOCITIES"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("COPUPLING_VV","VELOCITIES"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("COPUPLING_VV","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("COPUPLING_VV","VELOCITIES"); CHKERRQ(ierr);
  //VU
  ierr = mField.modify_finite_element_add_field_row("COPUPLING_VU","VELOCITIES"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("COPUPLING_VU","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("COPUPLING_VU","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("COPUPLING_VU","VELOCITIES"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","STIFFNESS"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","MASS"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","COPUPLING_VV"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","COPUPLING_VU"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  /***/
  //Declare problem

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"VELOCITIES"); CHKERRQ(ierr);


  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"STIFFNESS",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MASS",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"COPUPLING_VV",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"COPUPLING_VU",MBTET); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 1;
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
  order = 1;
  ierr = mField.set_field_order(0,MBTET,"VELOCITIES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"VELOCITIES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"VELOCITIES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"VELOCITIES",1); CHKERRQ(ierr);

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

  /****/
  //mesh partitioning 

  //partition
  ierr = mField.partition_problems("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //Get SideSet 1 and SideSet 2 defined in CUBIT
  Range SideSet1,SideSet2;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 2 : %u\n",SideSet2.size());

  //create matrices
  Vec D,F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&D); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  //TS
  moabTsCtx TsCtx(mField,"ELASTIC_MECHANICS");

  const double YoungModulus = 1;
  const double PoissonRatio = 0.;
  const double rho = 1;

  MyBC mybc;
  DynamicElasticFEMethod MyFE(moab,mField,Aij,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),rho,SideSet1,SideSet2,&mybc);


  moabTsCtx::loops_to_do_type& loops_to_do_Rhs = TsCtx.get_loops_to_do_IFunction();
  loops_to_do_Rhs.push_back(moabTsCtx::loop_pair_type("STIFFNESS",&MyFE));
  loops_to_do_Rhs.push_back(moabTsCtx::loop_pair_type("MASS",&MyFE));
  loops_to_do_Rhs.push_back(moabTsCtx::loop_pair_type("COPUPLING_VV",&MyFE));
  loops_to_do_Rhs.push_back(moabTsCtx::loop_pair_type("COPUPLING_VU",&MyFE));
  moabTsCtx::loops_to_do_type& loops_to_do_Mat = TsCtx.get_loops_to_do_IJacobian();
  loops_to_do_Mat.push_back(moabTsCtx::loop_pair_type("STIFFNESS",&MyFE));
  loops_to_do_Mat.push_back(moabTsCtx::loop_pair_type("MASS",&MyFE));
  loops_to_do_Mat.push_back(moabTsCtx::loop_pair_type("COPUPLING_VV",&MyFE));
  loops_to_do_Mat.push_back(moabTsCtx::loop_pair_type("COPUPLING_VU",&MyFE));

  //Monitor
  moabTsCtx::loops_to_do_type& loops_to_do_Monitor = TsCtx.get_loops_to_do_Monitor();
  loops_to_do_Monitor.push_back(moabTsCtx::loop_pair_type("STIFFNESS",&MyFE));
  loops_to_do_Monitor.push_back(moabTsCtx::loop_pair_type("COPUPLING_VV",&MyFE));

  TS ts;
  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);

  ierr = TSSetIFunction(ts,F,f_TSSetIFunction,&TsCtx); CHKERRQ(ierr);
  ierr = TSSetIJacobian(ts,Aij,Aij,f_TSSetIJacobian,&TsCtx); CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,f_TSMonitorSet,&TsCtx,PETSC_NULL); CHKERRQ(ierr);

  double ftime = 1;
  ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
  ierr = TSSetSolution(ts,D); CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);

  ierr = TSSolve(ts,D,&ftime); CHKERRQ(ierr);

  PetscInt steps,snesfails,rejects,nonlinits,linits;
  ierr = TSGetTimeStepNumber(ts,&steps); CHKERRQ(ierr);
  ierr = TSGetSNESFailures(ts,&snesfails); CHKERRQ(ierr);
  ierr = TSGetStepRejections(ts,&rejects); CHKERRQ(ierr);
  ierr = TSGetSNESIterations(ts,&nonlinits); CHKERRQ(ierr);
  ierr = TSGetKSPIterations(ts,&linits); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,
    "steps %D (%D rejected, %D SNES fails), ftime %G, nonlinits %D, linits %D\n",steps,rejects,snesfails,ftime,nonlinits,linits);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PostProcVertexMethod ent_method(moab);
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Row,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","STIFFNESS",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","STIFFNESS",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  PostProcL2VelocitiesFieldsAndGradientOnRefMesh fe_post_proc_velocities(moab);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","COPUPLING_VV",fe_post_proc_velocities);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_velocities.moab_post_proc.write_file("out_post_proc_velocities.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  //detroy matrices
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

}



