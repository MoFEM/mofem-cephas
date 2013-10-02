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

#include "ArcFEMethodForInterface.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

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

  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 3;
  }


  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
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

  //We need that for code profiling
  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  //Create MoFEM (Joseph) database
  moabField_Core core(moab);
  moabField& mField = core;

  BitRefLevel bit_level0;
  bit_level0.set(1);
  BitRefLevel bit_level1;
  bit_level1.set(2);

  BitRefLevel problem_bit_level = bit_level0;

  if(step == 1) {
    //ref meshset ref level 0
    ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);

    //Interface
    EntityHandle meshset_interface;
    ierr = mField.get_msId_meshset(4,SideSet,meshset_interface); CHKERRQ(ierr);
    ierr = mField.get_msId_3dENTS_sides(meshset_interface,true); CHKERRQ(ierr);
    // stl::bitset see for more details
    BitRefLevel bit_level_interface;
    bit_level_interface.set(0);
    ierr = mField.get_msId_3dENTS_split_sides(0,bit_level_interface,meshset_interface,true,true); CHKERRQ(ierr);
    EntityHandle meshset_level_interface;
    rval = moab.create_meshset(MESHSET_SET,meshset_level_interface); CHKERR_PETSC(rval);
    ierr = mField.refine_get_ents(bit_level_interface,meshset_level_interface); CHKERRQ(ierr);

    //add refined ent to cubit meshsets
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
      EntityHandle cubit_meshset = cubit_it->meshset; 
      ierr = mField.refine_get_childern(cubit_meshset,bit_level_interface,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
    }

    // stl::bitset see for more details
    EntityHandle meshset_level0;
    rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
    ierr = mField.seed_ref_level_3D(meshset_level_interface,bit_level0); CHKERRQ(ierr);
    ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

    /*
    ierr = mField.add_verices_in_the_middel_of_edges(meshset_level0,bit_level1); CHKERRQ(ierr);
    ierr = mField.refine_TET(meshset_level0,bit_level1); CHKERRQ(ierr);
    ierr = mField.refine_PRISM(meshset_level0,bit_level1); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
      EntityHandle cubit_meshset = cubit_it->meshset; 
      ierr = mField.refine_get_childern(cubit_meshset,bit_level1,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
    }*/

    /***/
    //Define problem

    //Fields
    ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
    ierr = mField.add_field("LAMBDA",NoField,1); CHKERRQ(ierr);
    //Field for ArcLenght
    ierr = mField.add_field("X0_DISPLACEMENT",H1,3); CHKERRQ(ierr);

    //FE
    ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ELASTIC","LAMBDA"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC","LAMBDA"); CHKERRQ(ierr); //this is for paremtis
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","LAMBDA"); CHKERRQ(ierr);

    //FE Interface
    ierr = mField.add_finite_element("INTERFACE"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);

    //FE ArcLenght
    ierr = mField.add_finite_element("ARC_LENGHT"); CHKERRQ(ierr);
    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ARC_LENGHT","LAMBDA"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ARC_LENGHT","LAMBDA"); CHKERRQ(ierr);
    //elem data
    ierr = mField.modify_finite_element_add_field_data("ARC_LENGHT","LAMBDA"); CHKERRQ(ierr);

    //define problems
    ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

    //set finite elements for problem
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","INTERFACE"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ARC_LENGHT"); CHKERRQ(ierr);

    //set refinment level for problem
    ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);

    /***/
    //Declare problem

    //add entitities (by tets) to the field
    ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

    //add finite elements entities
    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"ELASTIC",MBTET); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"INTERFACE",MBPRISM); CHKERRQ(ierr);


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
    //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
    ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

  }

  /****/
  //build database

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(problem_bit_level); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning 

  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS",false); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //print bcs
  ierr = mField.printDisplacementSet(); CHKERRQ(ierr);
  ierr = mField.printCubitForceSet(); CHKERRQ(ierr);

  //create matrices
  Vec F,D;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&F); CHKERRQ(ierr);
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  //Get SideSet 1 and SideSet 2 defined in CUBIT
  Range SideSet1,SideSet2,SideSet3;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(3,SideSet,2,SideSet3,true); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 2 : %u\n",SideSet2.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 3 : %u\n",SideSet3.size());

  //Assemble F and Aij
  const double YoungModulus = 1;
  const double PoissonRatio = 0.0;
  const double h = 1;
  const double beta = 0;
  const double ft = 1;
  const double Gf = 1;

  struct MyArcLenghtIntElemFEMethod: public ArcLenghtIntElemFEMethod {
    Range PostProcNodes;
    MyArcLenghtIntElemFEMethod(Interface& _moab,Mat &_Aij,Vec& _F,Vec& _D,
      ArcLenghtCtx *_arc_ptr): ArcLenghtIntElemFEMethod(_moab,_Aij,_F,_D,_arc_ptr) {

      Range all_nodes;
      rval = moab.get_entities_by_type(0,MBVERTEX,all_nodes,true); CHKERR_THROW(rval);
      for(Range::iterator nit = all_nodes.begin();nit!=all_nodes.end();nit++) {
	double coords[3];
	rval = moab.get_coords(&*nit,1,coords);  CHKERR_THROW(rval);
	if(fabs(coords[0]-5)<1e-6) {
	  PostProcNodes.insert(*nit);
	}
      }
      PetscPrintf(PETSC_COMM_WORLD,"Nb. PostProcNodes %lu\n",PostProcNodes.size());

    };

    PetscErrorCode potsProcessLoadPath() {
      PetscFunctionBegin;
      NumeredDofMoFEMEntity_multiIndex &numered_dofs_rows = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
      NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator lit;
      lit = numered_dofs_rows.get<FieldName_mi_tag>().find("LAMBDA");
      if(lit == numered_dofs_rows.get<FieldName_mi_tag>().end()) PetscFunctionReturn(0);
      Range::iterator nit = PostProcNodes.begin();
      for(;nit!=PostProcNodes.end();nit++) {
	NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator dit,hi_dit;
	dit = numered_dofs_rows.get<MoABEnt_mi_tag>().lower_bound(*nit);
	hi_dit = numered_dofs_rows.get<MoABEnt_mi_tag>().upper_bound(*nit);
	double coords[3];
	rval = moab.get_coords(&*nit,1,coords);  CHKERR_THROW(rval);
	for(;dit!=hi_dit;dit++) {
	  PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e -> ",lit->get_name().c_str(),lit->get_dof_rank(),lit->get_FieldData());
	  PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e ",dit->get_name().c_str(),dit->get_dof_rank(),dit->get_FieldData());
	  PetscPrintf(PETSC_COMM_WORLD,"-> %3.4f %3.4f %3.4f\n",coords[0],coords[1],coords[2]);
	}
      }
      PetscFunctionReturn(0);
    }

  };

  ArcLenghtCtx* ArcCtx = new ArcLenghtCtx(mField,"ELASTIC_MECHANICS");
  MyArcLenghtIntElemFEMethod* MyArcMethod_ptr = new MyArcLenghtIntElemFEMethod(moab,Aij,F,D,ArcCtx);
  MyArcLenghtIntElemFEMethod& MyArcMethod = *MyArcMethod_ptr;
  ArcLenghtSnesCtx SnesCtx(mField,"ELASTIC_MECHANICS",ArcCtx);

  CubitDisplacementDirihletBC myDirihletBC(mField,"ELASTIC_MECHANICS","DISPLACEMENT");
  ierr = myDirihletBC.Init(); CHKERRQ(ierr);

  ArcElasticFEMethod MyFE(mField,&myDirihletBC,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),ArcCtx);
  ArcInterfaceFEMethod IntMyFE(mField,&myDirihletBC,Aij,D,F,YoungModulus,h,beta,ft,Gf,ArcInterfaceFEMethod::ctx_IntLinearSoftening);

  PetscInt M,N;
  ierr = MatGetSize(Aij,&M,&N); CHKERRQ(ierr);
  PetscInt m,n;
  MatGetLocalSize(Aij,&m,&n);
  MatShellCtx* MatCtx = new MatShellCtx(mField,Aij,ArcCtx);
  Mat ShellAij;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,(void*)MatCtx,&ShellAij); CHKERRQ(ierr);
  ierr = MatShellSetOperation(ShellAij,MATOP_MULT,(void(*)(void))arc_lenght_mult_shell); CHKERRQ(ierr);

  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,ShellAij,Aij,SnesMat,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  KSP ksp;
  ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  PCShellCtx* PCCtx = new PCShellCtx(Aij,ShellAij,ArcCtx);
  ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
  ierr = PCShellSetContext(pc,PCCtx); CHKERRQ(ierr);
  ierr = PCShellSetApply(pc,pc_apply_arc_length); CHKERRQ(ierr);
  ierr = PCShellSetSetUp(pc,pc_setup_arc_length); CHKERRQ(ierr);

  //Rhs
  moabSnesCtx::loops_to_do_type& loops_to_do_Rhs = SnesCtx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MyFE));
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("INTERFACE",&IntMyFE));
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("ARC_LENGHT",&MyArcMethod));
  //Mat
  moabSnesCtx::loops_to_do_type& loops_to_do_Mat = SnesCtx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MyFE));
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("INTERFACE",&IntMyFE));
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("ARC_LENGHT",&MyArcMethod));

  int its_d = 6;
  double gamma = 0.5,reduction = 1;
  //step = 1;
  if(step == 1) {
    step_size = step_size_reduction;
  } else {
    reduction = step_size_reduction;
    step++;
  }

  if(step>1) {
    ierr = mField.set_local_VecCreateGhost("ELASTIC_MECHANICS",Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = mField.set_other_global_VecCreateGhost(
      "ELASTIC_MECHANICS","DISPLACEMENT","X0_DISPLACEMENT",Col,ArcCtx->x0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    double x0_nrm;
    ierr = VecNorm(ArcCtx->x0,NORM_2,&x0_nrm);  CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\tRead x0_nrm = %6.4e dlambda = %6.4e\n",x0_nrm,ArcCtx->dlambda);
    //
    ierr = ArcCtx->set_alpha_and_beta(1,0); CHKERRQ(ierr);
    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",IntMyFE);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyArcMethod);  CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  }


  for(;step<max_steps;step++) {

    if(step == 1) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Load Setp %D step_size = %6.4e\n",step,step_size); CHKERRQ(ierr);
      ierr = ArcCtx->set_s(step_size); CHKERRQ(ierr);
      ierr = ArcCtx->set_alpha_and_beta(0,1); CHKERRQ(ierr);
      ierr = VecCopy(D,ArcCtx->x0); CHKERRQ(ierr);
      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
      ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",IntMyFE);  CHKERRQ(ierr);
      ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ARC_LENGHT",MyArcMethod);  CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
      double dlambda;
      ierr = MyArcMethod.calculate_init_dlambda(&dlambda); CHKERRQ(ierr);
      ierr = MyArcMethod.set_dlambda_to_x(D,dlambda); CHKERRQ(ierr);
    } else if(step == 2) {
      ierr = ArcCtx->set_alpha_and_beta(1,0); CHKERRQ(ierr);
      ierr = MyArcMethod.calulate_dx_and_dlambda(D); CHKERRQ(ierr);
      ierr = MyArcMethod.calulate_lambda_int(step_size); CHKERRQ(ierr);
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

    //Distribute displacements on all processors
    ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

    //Update History and Calulate Residual
    //Tell Interface method that kappa is upadated
    ierr = IntMyFE.set_ctx_int(ArcInterfaceFEMethod::ctx_KappaUpdate); CHKERRQ(ierr);
    //run this on all processors, so we could save history tags on all parts and restart
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",IntMyFE,0,pcomm->size());  CHKERRQ(ierr);
    //Standard procedure
    ierr = IntMyFE.set_ctx_int(ArcInterfaceFEMethod::ctx_InterfaceNone); CHKERRQ(ierr);

    int its;
    ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);

    if(step > 1) {
      if(its>0) {
	reduction = pow((double)its_d/(double)(its+1),gamma);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"reduction step_size = %6.4e\n",reduction); CHKERRQ(ierr);
      } else reduction = 1;
    }

    //Save data on mesh
    ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = mField.set_other_global_VecCreateGhost(
      "ELASTIC_MECHANICS","DISPLACEMENT","X0_DISPLACEMENT",Col,ArcCtx->x0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    
    //
    PostProcVertexMethod ent_method(moab);
    ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Col,ent_method); CHKERRQ(ierr);
    //
    ierr = MyArcMethod.potsProcessLoadPath(); CHKERRQ(ierr);
    //
    if(step % 1 == 0) {
      if(pcomm->rank()==0) {
	ostringstream sss;
	sss << "restart_" << step << ".h5m";
	rval = moab.write_file(sss.str().c_str()); CHKERR_PETSC(rval);
      }
    }

  }

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",IntMyFE);  CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  PostProcVertexMethod ent_method(moab);
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Col,ent_method); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_res(moab,"DISPLACEMENT",F,"RESIDUAL");
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Col,ent_method_res); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","INTERFACE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }


  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);
  ierr = MatDestroy(&ShellAij); CHKERRQ(ierr);
  delete ArcCtx;
  delete MatCtx;
  delete PCCtx;
  delete MyArcMethod_ptr;

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

}

