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

#include <MoFEM.hpp>
using namespace MoFEM;

#include <DirichletBC.hpp>

#include <Projection10NodeCoordsOnField.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <SurfacePressure.hpp>
#include <NodalForce.hpp>
#include <FluidPressure.hpp>
#include <BodyForce.hpp>
#include <ThermalStressElement.hpp>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

#include <ElasticFEMethod.hpp>
#include <ElasticFEMethodDynamics.hpp>

#include <iostream>
#include <string>

using namespace boost::numeric;
using namespace ObosleteUsersModules;
using namespace MoFEM;
static char help[] = "...\n\n";

//See file users_modules/elasticity/TimeForceScale.hpp
#include <TimeForceScale.hpp>

int main(int argc, char *argv[]) {
  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Reade parameters from line command
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
    
  //Read mesh to MOAB
  const char *option;
  //option = "PARALLEL=BCAST_DELETE;"
      //"PARTITION=GEOM_DIMENSION,PARTITION_VAL=3,PARTITION_DISTRIBUTE";//;DEBUG_IO";
  option = "";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  //ref meshset ref level 0
  ierr = m_field.seed_ref_level_3D(0,0); CHKERRQ(ierr);

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  /***/
  //Define problem

  //Fields
  ierr = m_field.add_field("DISPLACEMENT",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("VELOCITIES",L2,3,MF_ZERO); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("STIFFNESS",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("MASS",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("COPUPLING_VV",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("COPUPLING_VU",MF_ZERO); CHKERRQ(ierr);


  //Define rows/cols and element data
  //STIFFNESS
  ierr = m_field.modify_finite_element_add_field_row("STIFFNESS","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("STIFFNESS","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("STIFFNESS","DISPLACEMENT"); CHKERRQ(ierr);

  //MASS
  ierr = m_field.modify_finite_element_add_field_row("MASS","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("MASS","VELOCITIES"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("MASS","DISPLACEMENT"); CHKERRQ(ierr);

  //COPUPLING
  //VV
  ierr = m_field.modify_finite_element_add_field_row("COPUPLING_VV","VELOCITIES"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("COPUPLING_VV","VELOCITIES"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("COPUPLING_VV","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("COPUPLING_VV","VELOCITIES"); CHKERRQ(ierr);
  //VU
  ierr = m_field.modify_finite_element_add_field_row("COPUPLING_VU","VELOCITIES"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("COPUPLING_VU","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("COPUPLING_VU","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("COPUPLING_VU","VELOCITIES"); CHKERRQ(ierr);

  //define problems
  ierr = m_field.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","STIFFNESS"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","MASS"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","COPUPLING_VV"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","COPUPLING_VU"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  /***/
  //Declare problem

  //add entitities (by tets) to the field
  ierr = m_field.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"VELOCITIES"); CHKERRQ(ierr);


  //add finite elements entities
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"STIFFNESS",MBTET); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MASS",MBTET); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"COPUPLING_VV",MBTET); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"COPUPLING_VU",MBTET); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt disp_order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_disp_order",&disp_order,&flg); CHKERRQ(ierr);
  if(flg!=PETSC_TRUE) {
    disp_order = 1;	
  }
  ierr = m_field.set_field_order(0,MBTET,"DISPLACEMENT",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DISPLACEMENT",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DISPLACEMENT",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
  PetscInt vel_order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_vel_order",&vel_order,&flg); CHKERRQ(ierr);
  if(flg!=PETSC_TRUE) {
    vel_order = 1;	
  }
  ierr = m_field.set_field_order(0,MBTET,"VELOCITIES",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"VELOCITIES",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"VELOCITIES",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"VELOCITIES",1); CHKERRQ(ierr);

  ierr = MetaNeummanForces::addNeumannBCElements(m_field,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = MetaNodalForces::addNodalForceElement(m_field,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","FORCE_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","FORCE_FE"); CHKERRQ(ierr);

  //body forces
  ierr = m_field.add_finite_element("BODY_FORCE",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","BODY_FORCE"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|BODYFORCESSET,it)) {
    Range tets;
    rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
    ierr = m_field.add_ents_to_finite_element_by_TETs(tets,"BODY_FORCE"); CHKERRQ(ierr);
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

  /****/
  //mesh partitioning 

  //partition
  ierr = m_field.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //print bcs
  ierr = m_field.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = m_field.print_cubit_materials_set(); CHKERRQ(ierr);

  //create matrices
  Vec D,F;
  ierr = m_field.VecCreateGhost("ELASTIC_MECHANICS",COL,&D); CHKERRQ(ierr);
  ierr = m_field.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = m_field.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  //TS
  TsCtx ts_ctx(m_field,"ELASTIC_MECHANICS");

  const double YoungModulus = 1;
  const double PoissonRatio = 0.;
  const double rho = 1;

  struct DynamicBCFEMethodPreAndPostProc: public DisplacementBCFEMethodPreAndPostProc {
    DynamicBCFEMethodPreAndPostProc(FieldInterface& _m_field,const string &_field_name,
    Mat &_Aij,Vec _X,Vec _F): DisplacementBCFEMethodPreAndPostProc(_m_field,_field_name,_Aij,_X,_F) {}

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      ierr = iNitalize(); CHKERRQ(ierr);
      if(dofsIndices.size()>0) {
	ierr = VecSetValues(ts_u,dofsIndices.size(),&dofsIndices[0],&dofsValues[0],INSERT_VALUES); CHKERRQ(ierr);
      }
      ierr = VecAssemblyBegin(ts_u); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(ts_u); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      switch(ts_ctx) {
	case CTX_TSSETIFUNCTION: {
	  ierr = VecAssemblyBegin(ts_F); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(ts_F); CHKERRQ(ierr);
	  for(vector<int>::iterator vit = dofsIndices.begin();vit!=dofsIndices.end();vit++) {
	    ierr = VecSetValue(ts_F,*vit,0,INSERT_VALUES); CHKERRQ(ierr);
	  }
	  ierr = VecAssemblyBegin(ts_F); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(ts_F); CHKERRQ(ierr);
	}
	break;
	case CTX_TSSETIJACOBIAN: {
	  ierr = MatAssemblyBegin(ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	  ierr = MatZeroRowsColumns(ts_B,dofsIndices.size(),&dofsIndices[0],1,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
	}
	break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"unknown snes stage");
      }
      PetscFunctionReturn(0);
    }

  };

  struct UpdateAndControl: public FEMethod {

    TS tS;
    int jacobianLag;
    UpdateAndControl(TS _ts): tS(_ts),jacobianLag(-1) {}

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      SNES snes;
      ierr = TSGetSNES(tS,&snes); CHKERRQ(ierr);
      ierr = SNESSetLagJacobian(snes,jacobianLag); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  TS ts;
  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);

  DynamicBCFEMethodPreAndPostProc my_dirichlet_bc(m_field,"DISPLACEMENT",Aij,D,F);
  DynamicElasticFEMethod my_fe(m_field,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),rho);
  UpdateAndControl update_and_control(ts);

  //Right hand side
  //preprocess
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&my_dirichlet_bc);
  //fe looops
  TsCtx::loops_to_do_type& loops_to_do_Rhs = ts_ctx.get_loops_to_do_IFunction();
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("STIFFNESS",&my_fe));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("MASS",&my_fe));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("COPUPLING_VV",&my_fe));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("COPUPLING_VU",&my_fe));
  //Neumann boundary conditions
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  string fe_name_str ="FORCE_FE";
  neumann_forces.insert(fe_name_str,new NeummanForcesSurface(m_field));
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = neumann_forces.at(fe_name_str).addForce("DISPLACEMENT",F,it->get_msId());  CHKERRQ(ierr);
    neumann_forces.at(fe_name_str).methodsOp.push_back(new TimeForceScale());
  }
  fe_name_str = "PRESSURE_FE";
  neumann_forces.insert(fe_name_str,new NeummanForcesSurface(m_field));
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
    ierr = neumann_forces.at(fe_name_str).addPreassure("DISPLACEMENT",F,it->get_msId()); CHKERRQ(ierr);
    neumann_forces.at(fe_name_str).methodsOp.push_back(new TimeForceScale());
  }
  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
  for(;mit!=neumann_forces.end();mit++) {
    if(mit->second->methodsOp.size() == 0) continue;
    loops_to_do_Rhs.push_back(TsCtx::loop_pair_type(mit->first,&mit->second->getLoopFe()));
  }
  //nodal forces
  boost::ptr_map<string,NodalForce> nodal_forces;
  fe_name_str ="FORCE_FE";
  nodal_forces.insert(fe_name_str,new NodalForce(m_field));
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = nodal_forces.at(fe_name_str).addForce("DISPLACEMENT",F,it->get_msId());  CHKERRQ(ierr);
    nodal_forces.at(fe_name_str).methodsOp.push_back(new TimeForceScale());
  }
  boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    loops_to_do_Rhs.push_back(TsCtx::loop_pair_type(fit->first,&fit->second->getLoopFe()));
  }
  //body forecs
  BodyFroceConstantField body_forces_methods(m_field);
  bool add_body_forces = false;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|BODYFORCESSET,it)) {
    add_body_forces = true;
    ierr = body_forces_methods.addBlock("DISPLACEMENT",F,it->get_msId()); CHKERRQ(ierr);
  }
  if(add_body_forces) {
    loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("BODY_FORCE",&body_forces_methods.getLoopFe()));
  }

  //postprocess
  ts_ctx.get_postProcess_to_do_IFunction().push_back(&my_dirichlet_bc);

  //Left hand side
  //preprocess
  ts_ctx.get_preProcess_to_do_IJacobian().push_back(&my_dirichlet_bc);
  //loops finire elements
  TsCtx::loops_to_do_type& loops_to_do_Mat = ts_ctx.get_loops_to_do_IJacobian();
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("STIFFNESS",&my_fe));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("MASS",&my_fe));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("COPUPLING_VV",&my_fe));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("COPUPLING_VU",&my_fe));
  //postrocess
  ts_ctx.get_postProcess_to_do_IJacobian().push_back(&my_dirichlet_bc);
  ts_ctx.get_postProcess_to_do_IJacobian().push_back(&update_and_control);

  //Monitor
  TsCtx::loops_to_do_type& loops_to_do_Monitor = ts_ctx.get_loops_to_do_Monitor();
  loops_to_do_Monitor.push_back(TsCtx::loop_pair_type("STIFFNESS",&my_fe));
  loops_to_do_Monitor.push_back(TsCtx::loop_pair_type("COPUPLING_VV",&my_fe));

  ierr = TSSetIFunction(ts,F,f_TSSetIFunction,&ts_ctx); CHKERRQ(ierr);
  ierr = TSSetIJacobian(ts,Aij,Aij,f_TSSetIJacobian,&ts_ctx); CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,f_TSMonitorSet,&ts_ctx,PETSC_NULL); CHKERRQ(ierr);

  double ftime = 1;
  ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
  ierr = TSSetSolution(ts,D); CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);

  //PetscReal alpha_m,alpha_f,gamma;
  //ierr = TSAlphaGetParams(ts,&alpha_m,&alpha_f,&gamma); CHKERRQ(ierr);
  //PetscPrintf(PETSC_COMM_WORLD,"alpha_m = %6.4e alpha_f = %6.4e gamma = %6.4e\n",alpha_m,alpha_f,gamma);

  ierr = TSSolve(ts,D); CHKERRQ(ierr);
  ierr = TSGetTime(ts,&ftime); CHKERRQ(ierr);

  PetscInt steps,snesfails,rejects,nonlinits,linits;
  ierr = TSGetTimeStepNumber(ts,&steps); CHKERRQ(ierr);
  ierr = TSGetSNESFailures(ts,&snesfails); CHKERRQ(ierr);
  ierr = TSGetStepRejections(ts,&rejects); CHKERRQ(ierr);
  ierr = TSGetSNESIterations(ts,&nonlinits); CHKERRQ(ierr);
  ierr = TSGetKSPIterations(ts,&linits); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,
    "steps %D (%D rejected, %D SNES fails), ftime %g, nonlinits %D, linits %D\n",
    steps,rejects,snesfails,ftime,nonlinits,linits);

  //Save data on mesh
  ierr = m_field.set_global_VecCreateGhost("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PostProcVertexMethod ent_method(moab);
  ierr = m_field.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",ROW,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = m_field.problem_get_FE("ELASTIC_MECHANICS","STIFFNESS",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab,"DISPLACEMENT");
  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","STIFFNESS",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  PostProcL2VelocitiesFieldsAndGradientOnRefMesh fe_post_proc_velocities(moab);
  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","COPUPLING_VV",fe_post_proc_velocities);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_velocities.moab_post_proc.write_file("out_post_proc_velocities.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  //detroy matrices
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);

  PetscFinalize();

}



