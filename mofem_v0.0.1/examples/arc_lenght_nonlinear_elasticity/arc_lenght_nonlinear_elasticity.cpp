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

#include "arc_lenght_nonlinear_elasticity.hpp"

static char help[] = "\
-my_file mesh file name\n\
-my_sr reduction of step size\n\
-my_ms maximal number of steps\n\n";

using namespace MoFEM;

int main(int argc, char *argv[]) {

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

  PetscScalar step_size_reduction;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_sr",&step_size_reduction,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    step_size_reduction = 1.;
  }

  PetscInt max_steps;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_ms",&max_steps,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    max_steps = 5;
  }
 
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //data stored on mesh for restart
  Tag th_step_size,th_step;
  double def_step_size = 1;
  rval = moab.tag_get_handle("_STEPSIZE",1,MB_TYPE_DOUBLE,th_step_size,MB_TAG_CREAT|MB_TAG_MESH,&def_step_size);  
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR(rval);
  int def_step = 1;
  rval = moab.tag_get_handle("_STEP",1,MB_TYPE_INTEGER,th_step,MB_TAG_CREAT|MB_TAG_MESH,&def_step);  
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR(rval);
  const void* tag_data_step_size[1];
  EntityHandle root = moab.get_root_set();
  rval = moab.tag_get_by_ptr(th_step_size,&root,1,tag_data_step_size); CHKERR_PETSC(rval);
  double& step_size = *(double *)tag_data_step_size[0];
  const void* tag_data_step[1];
  rval = moab.tag_get_by_ptr(th_step,&root,1,tag_data_step); CHKERR_PETSC(rval);
  int& step = *(int *)tag_data_step[0];
  //end of data stored for restart
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Start step %D and step_size = %6.4e\n",step,step_size); CHKERRQ(ierr);

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  moabField_Core core(moab);
  moabField& mField = core;

  BitRefLevel bit_level0;
  bit_level0.set(0);

  if(step == 1) {
    Range CubitSideSets_meshsets;
    ierr = mField.get_CubitBCType_meshsets(SideSet,CubitSideSets_meshsets); CHKERRQ(ierr);

    //ref meshset ref level 0
    ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
    EntityHandle meshset_level0;
    rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
    ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
    ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

    //Fields
    ierr = mField.add_field("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);
    ierr = mField.add_field("LAMBDA",NoField,1); CHKERRQ(ierr);

    //Field for ArcLenght
    ierr = mField.add_field("X0_SPATIAL_POSITION",H1,3); CHKERRQ(ierr);

    //FE
    ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("ARC_LENGHT"); CHKERRQ(ierr);

    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ELASTIC","LAMBDA"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC","LAMBDA"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","LAMBDA"); CHKERRQ(ierr);

    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ARC_LENGHT","LAMBDA"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ARC_LENGHT","LAMBDA"); CHKERRQ(ierr);
    //elem data
    ierr = mField.modify_finite_element_add_field_data("ARC_LENGHT","LAMBDA"); CHKERRQ(ierr);

    //define problems
    ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

    //set finite elements for problems
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ARC_LENGHT"); CHKERRQ(ierr);

    //set refinment level for problem
    ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

    //add entitities (by tets) to the field
    ierr = mField.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);

    //add finite elements entities
    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  
    //this entity will carray data for this finite element
    EntityHandle meshset_FE_ARC_LENGHT;
    rval = moab.create_meshset(MESHSET_SET,meshset_FE_ARC_LENGHT); CHKERR_PETSC(rval);
    //get LAMBDA field meshset
    EntityHandle meshset_field_LAMBDA = mField.get_field_meshset("LAMBDA");
    //add LAMBDA field meshset to finite element ARC_LENGHT
    rval = moab.add_entities(meshset_FE_ARC_LENGHT,&meshset_field_LAMBDA,1); CHKERR_PETSC(rval);
    //add finite element ARC_LENGHT meshset to refinment database (all ref bit leveles)
    ierr = mField.seed_ref_level_MESHSET(meshset_FE_ARC_LENGHT,BitRefLevel().set()); CHKERRQ(ierr);
    //finally add created meshset to the ARC_LENGHT finite element
    ierr = mField.add_ents_to_finite_element_by_MESHSET(meshset_FE_ARC_LENGHT,"ARC_LENGHT"); CHKERRQ(ierr);

    //set app. order
    ierr = mField.set_field_order(0,MBTET,"SPATIAL_POSITION",4); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"SPATIAL_POSITION",4); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"SPATIAL_POSITION",4); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);
  }

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition
  ierr = mField.partition_problems("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //create matrices
  Vec F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);
  
  ArcLenghtCtx* ArcCtx = new ArcLenghtCtx(mField,"ELASTIC_MECHANICS");

  PetscInt M,N;
  ierr = MatGetSize(Aij,&M,&N); CHKERRQ(ierr);
  PetscInt m,n;
  MatGetLocalSize(Aij,&m,&n);
  MatShellCtx* MatCtx = new MatShellCtx(mField,Aij,ArcCtx);
  Mat ShellAij;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,(void*)MatCtx,&ShellAij); CHKERRQ(ierr);
  ierr = MatShellSetOperation(ShellAij,MATOP_MULT,(void(*)(void))arc_lenght_mult_shell); CHKERRQ(ierr);

  if(step==1) {
    SetPositionsEntMethod set_positions(moab);
    ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Row,set_positions); CHKERRQ(ierr);
  }

  Range SideSet1,SideSet2,SideSet3,SideSet4;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,1,SideSet1,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(3,SideSet,2,SideSet3,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(4,SideSet,2,SideSet4,true); CHKERRQ(ierr);
  Range NodeSet1;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,NodeSet,0,NodeSet1,true); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"Nb. edges in SideSet 1 : %u\n",SideSet1.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 2 : %u\n",SideSet2.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 3 : %u\n",SideSet3.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 4 : %u\n",SideSet4.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. nodes in NodeSet 1 : %u\n",NodeSet1.size());

  const double YoungModulus = 1;
  const double PoissonRatio = 0.25;
  MyElasticFEMethod MyFE(moab,
    LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),
    ArcCtx,SideSet1,SideSet2,SideSet3,SideSet4,NodeSet1);

  ArcLenghtElemFEMethod MyArcMethod(moab,ArcCtx);
  ArcLenghtSnesCtx SnesCtx(mField,"ELASTIC_MECHANICS",ArcCtx);
  
  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,ShellAij,Aij,SnesMat,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  //
  /*ierr = SNESSetType(snes,SNESSHELL); CHKERRQ(ierr);
  ierr = SNESShellSetContext(snes,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESShellSetSolve(snes,snes_apply_arc_lenght); CHKERRQ(ierr);*/
  //

  KSP ksp;
  ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  PCShellCtx* PCCtx = new PCShellCtx(Aij,ShellAij,ArcCtx);
  ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
  ierr = PCShellSetContext(pc,PCCtx); CHKERRQ(ierr);
  ierr = PCShellSetApply(pc,pc_apply_arc_length); CHKERRQ(ierr);
  ierr = PCShellSetSetUp(pc,pc_setup_arc_length); CHKERRQ(ierr);

  moabSnesCtx::loops_to_do_type& loops_to_do_Rhs = SnesCtx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MyFE));
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("ARC_LENGHT",&MyArcMethod));
  moabSnesCtx::loops_to_do_type& loops_to_do_Mat = SnesCtx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MyFE));
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("ARC_LENGHT",&MyArcMethod));

  ierr = MyFE.set_t_val(-1); CHKERRQ(ierr);

  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = mField.set_local_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);


  if(step>1) {
    ierr = mField.set_other_global_VecCreateGhost(
      "ELASTIC_MECHANICS","SPATIAL_POSITION","X0_SPATIAL_POSITION",Row,ArcCtx->x0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    double x0_nrm;
    ierr = VecNorm(ArcCtx->x0,NORM_2,&x0_nrm);  CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\tRead x0_nrm = %6.4e dlambda = %6.4e\n",x0_nrm,ArcCtx->dlambda);
    //
    ierr = ArcCtx->set_alpha_and_beta(1,0); CHKERRQ(ierr);
    ierr = MyFE.set_x(D); CHKERRQ(ierr);
    ierr = MyFE.set_f(F); CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
    ierr = MyArcMethod.set_x(D); CHKERRQ(ierr);
    ierr = MyArcMethod.set_f(F); CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyArcMethod);  CHKERRQ(ierr);
  }

  int its_d = 6;
  double gamma = 0.5,reduction = 1;
  //step = 1;
  if(step == 1) {
    step_size = step_size_reduction;
  } else {
    reduction = step_size_reduction;
  }

  for(;step<max_steps;step++) {

    if(step == 1) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Load Setp %D step_size = %6.4e\n",step,step_size); CHKERRQ(ierr);
      ierr = ArcCtx->set_s(step_size); CHKERRQ(ierr);
      ierr = ArcCtx->set_alpha_and_beta(0,1); CHKERRQ(ierr);
      ierr = VecCopy(D,ArcCtx->x0); CHKERRQ(ierr);
      ierr = MyFE.set_x(D); CHKERRQ(ierr);
      ierr = MyFE.set_f(F); CHKERRQ(ierr);
      ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
      ierr = MyArcMethod.set_x(D); CHKERRQ(ierr);
      ierr = MyArcMethod.set_f(F); CHKERRQ(ierr);
      ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyArcMethod);  CHKERRQ(ierr);
      double dlambda;
      ierr = MyArcMethod.calculate_init_dlambda(&dlambda); CHKERRQ(ierr);
      ierr = MyArcMethod.set_dlambda_to_x(D,dlambda); CHKERRQ(ierr);
    } else if(step == 2) {
      ierr = ArcCtx->set_alpha_and_beta(1,0); CHKERRQ(ierr);
      ierr = MyArcMethod.calulate_dx_and_dlambda(D); CHKERRQ(ierr);
      step_size = sqrt(MyArcMethod.calulate_lambda_int());
      ierr = ArcCtx->set_s(step_size); CHKERRQ(ierr);
      double dlambda = ArcCtx->dlambda;
      double dx_nrm;
      ierr = VecNorm(ArcCtx->dx,NORM_2,&dx_nrm);  CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,
	"Load Setp %D step_size = %6.4e dlambda0 = %6.4e dx_nrm = %6.4e dx2 = %6.4e\n",
	step,step_size,dlambda,dx_nrm,ArcCtx->dx2); CHKERRQ(ierr);
      ierr = VecCopy(D,ArcCtx->x0); CHKERRQ(ierr);
      ierr = VecAXPY(D,1.,ArcCtx->dx); CHKERRQ(ierr);
      ierr = MyArcMethod.set_dlambda_to_x(D,dlambda); CHKERRQ(ierr);
    } else {
      ierr = MyArcMethod.calulate_dx_and_dlambda(D); CHKERRQ(ierr);
      step_size *= reduction;
      ierr = ArcCtx->set_s(step_size); CHKERRQ(ierr);
      double dlambda = reduction*ArcCtx->dlambda;
      double dx_nrm;
      ierr = VecScale(ArcCtx->dx,reduction); CHKERRQ(ierr);
      ierr = VecNorm(ArcCtx->dx,NORM_2,&dx_nrm);  CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,
	"Load Setp %D step_size = %6.4e dlambda0 = %6.4e dx_nrm = %6.4e dx2 = %6.4e\n",
	step,step_size,dlambda,dx_nrm,ArcCtx->dx2); CHKERRQ(ierr);
      ierr = VecCopy(D,ArcCtx->x0); CHKERRQ(ierr);
      ierr = VecAXPY(D,1.,ArcCtx->dx); CHKERRQ(ierr);
      ierr = MyArcMethod.set_dlambda_to_x(D,dlambda); CHKERRQ(ierr);
    }

    ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
    int its;
    ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);

    if(step > 1) {
      reduction = pow((double)its_d/(double)(its+1),gamma);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"reduction step_size = %6.4e\n",reduction); CHKERRQ(ierr);
    }

    //Save data on mesh
    ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = mField.set_other_global_VecCreateGhost(
      "ELASTIC_MECHANICS","SPATIAL_POSITION","X0_SPATIAL_POSITION",Row,ArcCtx->x0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    
    //
    PostProcDisplacementsEntMethod ent_method(moab,"SPATIAL_POSITION");
    ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Row,ent_method); CHKERRQ(ierr);
    //
    ierr = MyFE.potsProcessLoadPath(); CHKERRQ(ierr);
    //
    if(step % 1 == 0) {
      if(pcomm->rank()==0) {
	ostringstream sss;
	sss << "restart_" << step << ".h5m";
	rval = moab.write_file(sss.str().c_str()); CHKERR_PETSC(rval);
      }
    }

  }
  
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PostProcDisplacementsEntMethod ent_method(moab,"SPATIAL_POSITION");
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Row,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PostProcFieldsAndGradientOnRefMesh fe_post_proc_method(moab);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = MatDestroy(&ShellAij); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);
  delete MatCtx;
  delete PCCtx;
  delete ArcCtx;

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}



