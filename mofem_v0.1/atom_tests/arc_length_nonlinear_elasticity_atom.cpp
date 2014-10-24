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

static char help[] = "\
-my_file mesh file name\n\
-my_sr reduction of step size\n\
-my_ms maximal number of steps\n\n";

#include "SurfacePressureComplexForLazy.hpp"

#include "FEMethod_DriverComplexForLazy.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcNonLinearElasticityStresseOnRefindedMesh.hpp"
#include "Projection10NodeCoordsOnField.hpp"

//Rounding
#define RND_EPS 1e-6
double roundn(double n) {
  //break n into fractional part (fract) and integral part (intp)
  double fract, intp;
  fract = modf(n,&intp);
  // case where n approximates zero, set n to "positive" zero
  if(abs(intp)==0) {
    if(abs(fract)<=RND_EPS) {
      n=0.000;
    }
  }
  return n;
}

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

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

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  FieldCore core(moab);
  FieldInterface& mField = core;

  BitRefLevel bit_level0;
  bit_level0.set(0);

  if(step == 1) {
    Range CubitSIDESETs_meshsets;
    ierr = mField.get_Cubit_meshsets(SIDESET,CubitSIDESETs_meshsets); CHKERRQ(ierr);

    //ref meshset ref level 0
    ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
    EntityHandle meshset_level0;
    rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
    ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
    ierr = mField.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

    //Fields
    ierr = mField.add_field("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);
    ierr = mField.add_field("LAMBDA",NOFIELD,1); CHKERRQ(ierr);

    //Field for ArcLength
    ierr = mField.add_field("X0_SPATIAL_POSITION",H1,3); CHKERRQ(ierr);

    //FE
    ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("ARC_LENGTH"); CHKERRQ(ierr);

    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ELASTIC","LAMBDA"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC","LAMBDA"); CHKERRQ(ierr); //this is for parmetis
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","LAMBDA"); CHKERRQ(ierr);

    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);
    //elem data
    ierr = mField.modify_finite_element_add_field_data("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);

    //define problems
    ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

    //set finite elements for problems
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ARC_LENGTH"); CHKERRQ(ierr);

    //set refinment level for problem
    ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

    //add entitities (by tets) to the field
    ierr = mField.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);

    //add finite elements entities
    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  
    //this entity will carray data for this finite element
    EntityHandle meshset_FE_ARC_LENGTH;
    rval = moab.create_meshset(MESHSET_SET,meshset_FE_ARC_LENGTH); CHKERR_PETSC(rval);
    //get LAMBDA field meshset
    EntityHandle meshset_field_LAMBDA = mField.get_field_meshset("LAMBDA");
    //add LAMBDA field meshset to finite element ARC_LENGTH
    rval = moab.add_entities(meshset_FE_ARC_LENGTH,&meshset_field_LAMBDA,1); CHKERR_PETSC(rval);
    //add finite element ARC_LENGTH meshset to refinment database (all ref bit leveles)
    ierr = mField.seed_ref_level_MESHSET(meshset_FE_ARC_LENGTH,BitRefLevel().set()); CHKERRQ(ierr);
    //finally add created meshset to the ARC_LENGTH finite element
    ierr = mField.add_ents_to_finite_element_by_MESHSET(meshset_FE_ARC_LENGTH,"ARC_LENGTH",false); CHKERRQ(ierr);

    //set app. order
    int order = 1;
    ierr = mField.set_field_order(0,MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);

    //add neumman finite elemnets to add static boundary conditions
    ierr = mField.add_finite_element("NEUAMNN_FE"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","NEUAMNN_FE"); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
      Range tris;
      rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TRIs(tris,"NEUAMNN_FE"); CHKERRQ(ierr);
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {
      Range tris;
      rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TRIs(tris,"NEUAMNN_FE"); CHKERRQ(ierr);
    }

    #ifdef NONLINEAR_TEMPERATUTE

    ierr = mField.add_field("TEMPERATURE",H1,1); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","TEMPERATURE"); CHKERRQ(ierr);

    ierr = mField.add_ents_to_field_by_TETs(0,"TEMPERATURE"); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTET,"TEMPERATURE",1); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"TEMPERATURE",1); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"TEMPERATURE",1); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"TEMPERATURE",1); CHKERRQ(ierr);

    #endif //NONLINEAR_TEMPERATUTE

  }

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  if(step==1) {
    EntityHandle node = 0;
    double coords[3];
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"SPATIAL_POSITION",dof_ptr)) {
      if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
      EntityHandle ent = dof_ptr->get_ent();
      int dof_rank = dof_ptr->get_dof_rank();
      double &fval = dof_ptr->get_FieldData();
      if(node!=ent) {
	rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
      }
      fval = coords[dof_rank];
    }
  }
  #ifdef NONLINEAR_TEMPERATUTE

  if(step==1) {
    EntityHandle node = 0;
    double coords[3];
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"TEMPERATURE",dof_ptr)) {
      if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
      EntityHandle ent = dof_ptr->get_ent();
      if(node!=ent) {
	rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
      }
      double &fval = dof_ptr->get_FieldData();
      fval = coords[0];
    }
  }
  
  #endif //NONLINEAR_TEMPERATUTE


  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);


  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //print bcs
  ierr = mField.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = mField.print_cubit_pressure_set(); CHKERRQ(ierr);
  ierr = mField.print_cubit_force_set(); CHKERRQ(ierr);

  //print block sets with materials
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);

  //create matrices
  Vec F,F_lambda_for_neumann_forces;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",COL,&F); CHKERRQ(ierr);
  ierr = VecDuplicate(F,&F_lambda_for_neumann_forces); CHKERRQ(ierr);
  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);
    
  ArcLengthCtx* arc_ctx = new ArcLengthCtx(mField,"ELASTIC_MECHANICS");

  PetscInt M,N;
  ierr = MatGetSize(Aij,&M,&N); CHKERRQ(ierr);
  PetscInt m,n;
  MatGetLocalSize(Aij,&m,&n);
  ArcLengthMatShell* mat_ctx = new ArcLengthMatShell(mField,Aij,arc_ctx,"ELASTIC_MECHANICS");
  Mat ShellAij;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,(void*)mat_ctx,&ShellAij); CHKERRQ(ierr);
  ierr = MatShellSetOperation(ShellAij,MATOP_MULT,(void(*)(void))arc_length_mult_shell); CHKERRQ(ierr);

  double young_modulus = 1;
  double poisson_ratio = 0.25;
	
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
    cout << endl << *it << endl;
    //Get block name
    string name = it->get_Cubit_name();
    if (name.compare(0,11,"MAT_ELASTIC") == 0) {
      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      cout << mydata;
      young_modulus=mydata.data.Young;
      poisson_ratio=mydata.data.Poisson;
    }
  }

  ArcLengthElemFEMethod* arc_method_ptr = new ArcLengthElemFEMethod(moab,arc_ctx);
  ArcLengthElemFEMethod& arc_method = *arc_method_ptr;

  #ifndef NONLINEAR_TEMPERATUTE
  NonLinearSpatialElasticFEMthod fe(mField,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  set_PhysicalEquationNumber(neohookean);
  #else 
  NonLinearSpatialElasticFEMthod fe(mField,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),arc_ctx);
  set_PhysicalEquationNumber(neohookean);
  set_ThermalDeformationEquationNumber(linear_expansion_true_volume);
  fe.thermal_expansion = 1;
  #endif

  double scaled_reference_load = 1;
  double *scale_lhs = &(arc_ctx->get_FieldData());
  double *scale_rhs = &(scaled_reference_load);
  NeummanForcesSurfaceComplexForLazy neumann_forces(mField,Aij,F_lambda_for_neumann_forces,scale_lhs,scale_rhs);
  NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE &fe_neumann = neumann_forces.getLoopSpatialFe();
  fe_neumann.uSeF = true;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
    ierr = fe_neumann.addForce(it->get_msId()); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {
    ierr = fe_neumann.addPreassure(it->get_msId()); CHKERRQ(ierr);
  }
  SpatialPositionsBCFEMethodPreAndPostProc my_dirichlet_bc(mField,"SPATIAL_POSITION",Aij,D,F);
  ierr = mField.get_problem("ELASTIC_MECHANICS",&my_dirichlet_bc.problemPtr); CHKERRQ(ierr);
  ierr = my_dirichlet_bc.iNitalize(); CHKERRQ(ierr);

  struct MyPrePostProcessFEMethod: public FieldInterface::FEMethod {
    
    FieldInterface& mField;
    ArcLengthCtx *arc_ptr;
    SpatialPositionsBCFEMethodPreAndPostProc *bC;
    Vec F_lambda_for_neumann_forces;


    MyPrePostProcessFEMethod(FieldInterface& _mField,
      ArcLengthCtx *_arc_ptr,SpatialPositionsBCFEMethodPreAndPostProc *bc): 
      mField(_mField),arc_ptr(_arc_ptr),bC(bc),F_lambda_for_neumann_forces(PETSC_NULL) {}

    MyPrePostProcessFEMethod(FieldInterface& _mField,
      ArcLengthCtx *_arc_ptr,SpatialPositionsBCFEMethodPreAndPostProc *bc,
      Vec _F_lambda_for_neumann_forces): 
      mField(_mField),arc_ptr(_arc_ptr),bC(bc),
      F_lambda_for_neumann_forces(_F_lambda_for_neumann_forces) {}

  
    PetscErrorCode ierr;
      
      PetscErrorCode preProcess() {
        PetscFunctionBegin;
        
	//PetscAttachDebugger();
        switch(snes_ctx) {
          case CTX_SNESSETFUNCTION: {
            ierr = VecZeroEntries(snes_f); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecZeroEntries(arc_ptr->F_lambda); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    if(F_lambda_for_neumann_forces!=PETSC_NULL) {
	      ierr = VecZeroEntries(F_lambda_for_neumann_forces); CHKERRQ(ierr);
	      ierr = VecGhostUpdateBegin(F_lambda_for_neumann_forces,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	      ierr = VecGhostUpdateEnd(F_lambda_for_neumann_forces,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    }
          }
          break;
          default:
            SETERRQ(PETSC_COMM_SELF,1,"not implemented");
        }
        
        PetscFunctionReturn(0);
      }
      
      PetscErrorCode postProcess() {
        PetscFunctionBegin;
        switch(snes_ctx) {
          case CTX_SNESSETFUNCTION: {
	    //snes_f
            ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	    //F_lambda
            ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	    if(F_lambda_for_neumann_forces!=PETSC_NULL) {
	      ierr = VecGhostUpdateBegin(F_lambda_for_neumann_forces,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	      ierr = VecGhostUpdateEnd(F_lambda_for_neumann_forces,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	      ierr = VecAssemblyBegin(F_lambda_for_neumann_forces); CHKERRQ(ierr);
	      ierr = VecAssemblyEnd(F_lambda_for_neumann_forces); CHKERRQ(ierr);
	      ierr = VecAXPY(arc_ptr->F_lambda,1,F_lambda_for_neumann_forces); CHKERRQ(ierr);
	      //add F_lambda
	      ierr = VecAXPY(snes_f,arc_ptr->get_FieldData(),F_lambda_for_neumann_forces); CHKERRQ(ierr);
	    }
	    for(vector<int>::iterator vit = bC->dofsIndices.begin();vit!=bC->dofsIndices.end();vit++) {
	      ierr = VecSetValue(arc_ptr->F_lambda,*vit,0,INSERT_VALUES); CHKERRQ(ierr);
	    }
	    ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	    ierr = VecDot(arc_ptr->F_lambda,arc_ptr->F_lambda,&arc_ptr->F_lambda2); CHKERRQ(ierr);
	    PetscPrintf(PETSC_COMM_WORLD,"\tFlambda2 = %6.4e\n",arc_ptr->F_lambda2);
	    PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",arc_ptr->get_FieldData());  
	    double fnorm;
	    ierr = VecNorm(snes_f,NORM_2,&fnorm); CHKERRQ(ierr);	
	    PetscPrintf(PETSC_COMM_WORLD,"\tfnorm = %6.4e\n",fnorm);  
          }
          break;
          default:
            SETERRQ(PETSC_COMM_SELF,1,"not implemented");
        }
        
        PetscFunctionReturn(0);
      }
  };

  MyPrePostProcessFEMethod pre_post_method(mField,arc_ctx,&my_dirichlet_bc,F_lambda_for_neumann_forces);
  ArcLengthSnesCtx snes_ctx(mField,"ELASTIC_MECHANICS",arc_ctx);
  
  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,ShellAij,Aij,SnesMat,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  PetscReal my_tol;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_tol",&my_tol,&flg); CHKERRQ(ierr);
  if(flg == PETSC_TRUE) {
    PetscReal atol,rtol,stol;
    PetscInt maxit,maxf;
    ierr = SNESGetTolerances(snes,&atol,&rtol,&stol,&maxit,&maxf); CHKERRQ(ierr);
    atol = my_tol;
    rtol = atol*1e2;
    ierr = SNESSetTolerances(snes,atol,rtol,stol,maxit,maxf); CHKERRQ(ierr);
  }

  KSP ksp;
  ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  PCShellCtx* pc_ctx = new PCShellCtx(Aij,ShellAij,arc_ctx);
  ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
  ierr = PCShellSetContext(pc,pc_ctx); CHKERRQ(ierr);
  ierr = PCShellSetApply(pc,pc_apply_arc_length); CHKERRQ(ierr);
  ierr = PCShellSetSetUp(pc,pc_setup_arc_length); CHKERRQ(ierr);

  if(flg == PETSC_TRUE) {
    PetscReal rtol,atol,dtol;
    PetscInt maxits;
    ierr = KSPGetTolerances(ksp,&rtol,&atol,&dtol,&maxits); CHKERRQ(ierr);
    atol = my_tol*1e-2;
    rtol = atol*1e-2;
    ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxits); CHKERRQ(ierr);
  }


  SnesCtx::loops_to_do_type& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
  snes_ctx.get_preProcess_to_do_Rhs().push_back(&my_dirichlet_bc);
  snes_ctx.get_preProcess_to_do_Rhs().push_back(&pre_post_method);
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ELASTIC",&fe));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_neumann));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&arc_method));
  snes_ctx.get_postProcess_to_do_Rhs().push_back(&pre_post_method);
  snes_ctx.get_postProcess_to_do_Rhs().push_back(&my_dirichlet_bc);

  SnesCtx::loops_to_do_type& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
  snes_ctx.get_preProcess_to_do_Mat().push_back(&my_dirichlet_bc);
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ELASTIC",&fe));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_neumann));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&arc_method));
  snes_ctx.get_postProcess_to_do_Mat().push_back(&my_dirichlet_bc);

  ierr = mField.set_local_VecCreateGhost("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  if(step>1) {
    ierr = mField.set_other_global_VecCreateGhost(
      "ELASTIC_MECHANICS","SPATIAL_POSITION","X0_SPATIAL_POSITION",COL,arc_ctx->x0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    double x0_nrm;
    ierr = VecNorm(arc_ctx->x0,NORM_2,&x0_nrm);  CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\tRead x0_nrm = %6.4e dlambda = %6.4e\n",x0_nrm,arc_ctx->dlambda);
    ierr = arc_ctx->set_alpha_and_beta(1,0); CHKERRQ(ierr);
  } else {
    ierr = arc_ctx->set_s(0); CHKERRQ(ierr);
    ierr = arc_ctx->set_alpha_and_beta(0,1); CHKERRQ(ierr);
  }
  ierr = SnesRhs(snes,D,F,&snes_ctx); CHKERRQ(ierr);

  int its_d;
  ierr = PetscOptionsGetInt("","-my_its_d",&its_d,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    its_d = 6;
  }

  double gamma = 0.5,reduction = 1;
  //step = 1;
  if(step == 1) {
    step_size = step_size_reduction;
  } else {
    reduction = step_size_reduction;
    step++;
  }

  Vec D0,x00;
  ierr = VecDuplicate(D,&D0); CHKERRQ(ierr);
  ierr = VecDuplicate(arc_ctx->x0,&x00); CHKERRQ(ierr);
  bool converged_state = false;

  for(;step<max_steps;step++) {

    ierr = VecCopy(D,D0); CHKERRQ(ierr);
    ierr = VecCopy(arc_ctx->x0,x00); CHKERRQ(ierr);

    if(step == 1) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Load Step %D step_size = %6.4e\n",step,step_size); CHKERRQ(ierr);
      ierr = arc_ctx->set_s(step_size); CHKERRQ(ierr);
      ierr = arc_ctx->set_alpha_and_beta(0,1); CHKERRQ(ierr);
      ierr = VecCopy(D,arc_ctx->x0); CHKERRQ(ierr);
      double dlambda;
      ierr = arc_method.calculate_init_dlambda(&dlambda); CHKERRQ(ierr);
      ierr = arc_method.set_dlambda_to_x(D,dlambda); CHKERRQ(ierr);
      //ierr = MyFE.potsProcessLoadPath(); CHKERRQ(ierr);
    } else if(step == 2) {
      ierr = arc_ctx->set_alpha_and_beta(1,0); CHKERRQ(ierr);
      ierr = arc_method.calculate_dx_and_dlambda(D); CHKERRQ(ierr);
      step_size = sqrt(arc_method.calculate_lambda_int());
      ierr = arc_ctx->set_s(step_size); CHKERRQ(ierr);
      double dlambda = arc_ctx->dlambda;
      double dx_nrm;
      ierr = VecNorm(arc_ctx->dx,NORM_2,&dx_nrm);  CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,
	"Load Step %D step_size = %6.4e dlambda0 = %6.4e dx_nrm = %6.4e dx2 = %6.4e\n",
	step,step_size,dlambda,dx_nrm,arc_ctx->dx2); CHKERRQ(ierr);
      ierr = VecCopy(D,arc_ctx->x0); CHKERRQ(ierr);
      ierr = VecAXPY(D,1.,arc_ctx->dx); CHKERRQ(ierr);
      ierr = arc_method.set_dlambda_to_x(D,dlambda); CHKERRQ(ierr);
    } else {
      ierr = arc_method.calculate_dx_and_dlambda(D); CHKERRQ(ierr);
      step_size *= reduction;
      ierr = arc_ctx->set_s(step_size); CHKERRQ(ierr);
      double dlambda = reduction*arc_ctx->dlambda;
      double dx_nrm;
      ierr = VecScale(arc_ctx->dx,reduction); CHKERRQ(ierr);
      ierr = VecNorm(arc_ctx->dx,NORM_2,&dx_nrm);  CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,
	"Load Step %D step_size = %6.4e dlambda0 = %6.4e dx_nrm = %6.4e dx2 = %6.4e\n",
	step,step_size,dlambda,dx_nrm,arc_ctx->dx2); CHKERRQ(ierr);
      ierr = VecCopy(D,arc_ctx->x0); CHKERRQ(ierr);
      ierr = VecAXPY(D,1.,arc_ctx->dx); CHKERRQ(ierr);
      ierr = arc_method.set_dlambda_to_x(D,dlambda); CHKERRQ(ierr);
    }

    ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
    int its;
    ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);

    SNESConvergedReason reason;
    ierr = SNESGetConvergedReason(snes,&reason); CHKERRQ(ierr);
    if(reason < 0) {

      ierr = VecCopy(D0,D); CHKERRQ(ierr);
      ierr = VecCopy(x00,arc_ctx->x0); CHKERRQ(ierr);

      double x0_nrm;
      ierr = VecNorm(arc_ctx->x0,NORM_2,&x0_nrm);  CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\tRead x0_nrm = %6.4e dlambda = %6.4e\n",x0_nrm,arc_ctx->dlambda);
      ierr = arc_ctx->set_alpha_and_beta(1,0); CHKERRQ(ierr);

      
      reduction = 0.1;
      converged_state = false;
      
      continue;

    } else {
      if(step > 1 && converged_state) {
	reduction = pow((double)its_d/(double)(its+1),gamma);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"reduction step_size = %6.4e\n",reduction); CHKERRQ(ierr);
      }
      //Save data on mesh
      ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = mField.set_other_global_VecCreateGhost(
	"ELASTIC_MECHANICS","SPATIAL_POSITION","X0_SPATIAL_POSITION",COL,arc_ctx->x0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      converged_state = true;
      
    }
    
    PostProcVertexMethod ent_method(moab,"SPATIAL_POSITION");
    ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",COL,ent_method); CHKERRQ(ierr);

  }


  {

    //Open mesh_file_name.txt for writing
    ofstream myfile;
    #ifndef NONLINEAR_TEMPERATUTE
      myfile.open("arc_length_nonlinear_elasticity_atom.txt");
    #else 
      myfile.open("arc_length_nonlinear_elasticity_thermal_atom.txt");
    #endif
    
    //Output displacements
    cout << "<<<< Displacements (X-Translation, Y-Translation, Z-Translation) >>>>>" << endl;
    myfile << "<<<< Displacements (X-Translation, Y-Translation, Z-Translation) >>>>>" << endl;
    
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"SPATIAL_POSITION",dof_ptr))
    {
        if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
  
	double coords[3];
	EntityHandle ent = dof_ptr->get_ent();
	rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
        
        if(dof_ptr->get_dof_rank()==0)
        {
            //Round and truncate to 3 decimal places
            double fval = dof_ptr->get_FieldData()-coords[dof_ptr->get_dof_rank()];
            cout << boost::format("%.3lf") % roundn(fval) << "  ";
            myfile << boost::format("%.3lf") % roundn(fval) << "  ";
        }
        if(dof_ptr->get_dof_rank()==1)
        {
            //Round and truncate to 3 decimal places
            double fval = dof_ptr->get_FieldData()-coords[dof_ptr->get_dof_rank()];
            cout << boost::format("%.3lf") % roundn(fval) << "  ";
            myfile << boost::format("%.3lf") % roundn(fval) << "  ";
        }
        if(dof_ptr->get_dof_rank()==2)
        {
            //Round and truncate to 3 decimal places
            double fval = dof_ptr->get_FieldData()-coords[dof_ptr->get_dof_rank()];
            cout << boost::format("%.3lf") % roundn(fval) << endl;
            myfile << boost::format("%.3lf") % roundn(fval) << endl;
        }
        
    }

  }


  ierr = VecDestroy(&D0); CHKERRQ(ierr);
  ierr = VecDestroy(&x00); CHKERRQ(ierr);
  
  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&F_lambda_for_neumann_forces); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = MatDestroy(&ShellAij); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  delete mat_ctx;
  delete pc_ctx;
  delete arc_ctx;
  delete arc_method_ptr;

  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}



