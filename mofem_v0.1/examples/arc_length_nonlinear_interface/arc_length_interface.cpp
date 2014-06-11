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

#include "NonLinearFEMethodInterface.hpp"
#include "SurfacePressure.hpp"
#include "BodyForce.hpp"

#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

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
    order = 2;
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
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;

  Tag th_my_ref_level;
  BitRefLevel def_bit_level = 0;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",sizeof(BitRefLevel),MB_TYPE_OPAQUE,
    th_my_ref_level,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level); 
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;
  BitRefLevel problem_bit_level = bit_level0;

  if(step == 1) {
    //ref meshset ref level 0
    ierr = mField.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);
    vector<BitRefLevel> bit_levels;
    bit_levels.push_back(BitRefLevel().set(0));

    int ll = 1;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|InterfaceSet,cit)) {
    //for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,cit)) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Insert Interface %d\n",cit->get_msId()); CHKERRQ(ierr);
      EntityHandle cubit_meshset = cit->get_meshset();
      {
	//get tet enties form back bit_level
	EntityHandle ref_level_meshset = 0;
	rval = moab.create_meshset(MESHSET_SET,ref_level_meshset); CHKERR_PETSC(rval);
	ierr = mField.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBTET,ref_level_meshset); CHKERRQ(ierr);
	ierr = mField.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBPRISM,ref_level_meshset); CHKERRQ(ierr);
	Range ref_level_tets;
	rval = moab.get_entities_by_handle(ref_level_meshset,ref_level_tets,true); CHKERR_PETSC(rval);
	//get faces and test to split
	ierr = mField.get_msId_3dENTS_sides(cubit_meshset,bit_levels.back(),true,0); CHKERRQ(ierr);
	//set new bit level
	bit_levels.push_back(BitRefLevel().set(ll++));
	//split faces and 
	ierr = mField.get_msId_3dENTS_split_sides(ref_level_meshset,bit_levels.back(),cubit_meshset,true,true,0); CHKERRQ(ierr);
	//clean meshsets
	rval = moab.delete_entities(&ref_level_meshset,1); CHKERR_PETSC(rval);
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
    
    bit_level0 = bit_levels.back();
    problem_bit_level = bit_level0;

    /***/
    //Define problem

    //Fields
    ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
    ierr = mField.add_field("LAMBDA",NoField,1); CHKERRQ(ierr);
    //Field for ArcLength
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

    //FE ArcLength
    ierr = mField.add_finite_element("ARC_LENGTH"); CHKERRQ(ierr);
    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);
    //elem data
    ierr = mField.modify_finite_element_add_field_data("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);

    //define problems
    ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

    //set finite elements for problem
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","INTERFACE"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ARC_LENGTH"); CHKERRQ(ierr);

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
    //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
    ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

    /*//reduce level of approximation for entities on inetrface
    Range prims;
    ierr = mField.get_entities_by_type_and_ref_level(problem_bit_level,BitRefLevel().set(),MBPRISM,prims); CHKERRQ(ierr);
    Range prims_faces;
    rval = mField.get_moab().get_adjacencies(prims,2,false,prims_faces,Interface::UNION); CHKERR_PETSC(rval);
    Range prims_faces_edges;
    rval = mField.get_moab().get_adjacencies(prims_faces,1,false,prims_faces_edges,Interface::UNION); CHKERR_PETSC(rval);
    ierr = mField.set_field_order(prims_faces,"DISPLACEMENT",order>1 ? order-1 : 0); CHKERRQ(ierr);
    ierr = mField.set_field_order(prims_faces_edges,"DISPLACEMENT",order>1 ? order-1 : 0); CHKERRQ(ierr);*/

    //Elements with boundary conditions
    ierr = MetaNeummanForces::addNeumannBCElements(mField,"ELASTIC_MECHANICS","DISPLACEMENT"); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|ForceSet,it)) {
      ostringstream fe_name;
      fe_name << "FORCE_FE_" << it->get_msId();
      ierr = mField.add_finite_element(fe_name.str()); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row(fe_name.str(),"LAMBDA"); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col(fe_name.str(),"LAMBDA"); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data(fe_name.str(),"LAMBDA"); CHKERRQ(ierr);
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|PressureSet,it)) {
      ostringstream fe_name;
      fe_name << "PRESSURE_FE_" << it->get_msId();
      ierr = mField.add_finite_element(fe_name.str()); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row(fe_name.str(),"LAMBDA"); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col(fe_name.str(),"LAMBDA"); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data(fe_name.str(),"LAMBDA"); CHKERRQ(ierr);
    }
    ierr = mField.add_finite_element("BODY_FORCE"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("BODY_FORCE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","BODY_FORCE"); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|Block_BodyForcesSet,it)) {
      Range tets;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TETs(tets,"BODY_FORCE"); CHKERRQ(ierr);
    }
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
  ierr = mField.printCubitDisplacementSet(); CHKERRQ(ierr);
  ierr = mField.printCubitForceSet(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = mField.printCubitMaterials(); CHKERRQ(ierr);

  //create matrices
  Vec F,F_body_force,D;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&F); CHKERRQ(ierr);
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = VecDuplicate(F,&F_body_force); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  //Assemble F and Aij
  const double YoungModulus = 1;
  const double PoissonRatio = 0.0;
  const double h = 1;
  const double beta = 0;
  const double ft = 1;
  const double Gf = 1;

  struct MyArcLengthIntElemFEMethod: public ArcLengthIntElemFEMethod {
    FieldInterface& mField;
    Range PostProcNodes;
    MyArcLengthIntElemFEMethod(FieldInterface& _mField,Mat &_Aij,Vec& _F,Vec& _D,
      ArcLengthCtx *_arc_ptr): ArcLengthIntElemFEMethod(_mField.get_moab(),_Aij,_F,_D,_arc_ptr),mField(_mField) {

      for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"LoadPath",cit)) {
	EntityHandle meshset = cit->get_meshset();
	Range nodes;
	rval = moab.get_entities_by_type(meshset,MBVERTEX,nodes,true); CHKERR_THROW(rval);
	PostProcNodes.merge(nodes);
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

  struct MyPrePostProcessFEMethod: public FieldInterface::FEMethod {
    
    FieldInterface& mField;
    Mat &Aij;
    Vec &F;
    Vec &F_body_force;
    ArcLengthCtx *arc_ptr;
    
    MyPrePostProcessFEMethod(FieldInterface& _mField,
      Mat _Aij,Vec _F,Vec _F_body_force,ArcLengthCtx *_arc_ptr):
      mField(_mField),Aij(_Aij),F(_F),F_body_force(_F_body_force),arc_ptr(_arc_ptr) {}
  
    PetscErrorCode ierr;
      
      PetscErrorCode preProcess() {
        PetscFunctionBegin;
        
        switch(snes_ctx) {
          case ctx_SNESNone:
          case ctx_SNESSetFunction: {
            //F_lambda
            ierr = VecZeroEntries(F); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          }
          break;
          case ctx_SNESSetJacobian: {
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
          case ctx_SNESNone: {}
          case ctx_SNESSetFunction: {
            ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
	    //add F_lambda
	    ierr = VecAXPY(F,-arc_ptr->get_FieldData(),arc_ptr->F_lambda); CHKERRQ(ierr);
	    ierr = VecAXPY(F_body_force,-1.,F_body_force); CHKERRQ(ierr);
	    PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",arc_ptr->get_FieldData());  
	    //snes_f norm
	    double fnorm;
	    ierr = VecNormBegin(F,NORM_2,&fnorm); CHKERRQ(ierr);	
	    ierr = VecNormEnd(F,NORM_2,&fnorm);CHKERRQ(ierr);
	    PetscPrintf(PETSC_COMM_WORLD,"\tfnorm = %6.4e\n",fnorm);  
          }
          break;
          case ctx_SNESSetJacobian: {
          }
          break;
          default:
            SETERRQ(PETSC_COMM_SELF,1,"not implemented");
        }
        
        PetscFunctionReturn(0);
      }
  };

  ArcLengthCtx* ArcCtx = new ArcLengthCtx(mField,"ELASTIC_MECHANICS");
  MyArcLengthIntElemFEMethod* MyArcMethod_ptr = new MyArcLengthIntElemFEMethod(mField,Aij,F,D,ArcCtx);
  MyArcLengthIntElemFEMethod& MyArcMethod = *MyArcMethod_ptr;
  ArcLengthSnesCtx SnesCtx(mField,"ELASTIC_MECHANICS",ArcCtx);
  MyPrePostProcessFEMethod PrePostFE(mField,Aij,F,F_body_force,ArcCtx);

  CubitDisplacementDirihletBC myDirihletBC(mField,"ELASTIC_MECHANICS","DISPLACEMENT");
  ierr = myDirihletBC.Init(); CHKERRQ(ierr);

  ElasticFEMethod MyFE(mField,&myDirihletBC,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  NonLinearInterfaceFEMethod IntMyFE(mField,&myDirihletBC,Aij,D,F,YoungModulus,h,beta,ft,Gf,NonLinearInterfaceFEMethod::ctx_IntLinearSoftening);

  PetscInt M,N;
  ierr = MatGetSize(Aij,&M,&N); CHKERRQ(ierr);
  PetscInt m,n;
  MatGetLocalSize(Aij,&m,&n);
  ArcLengthMatShell* MatCtx = new ArcLengthMatShell(mField,Aij,ArcCtx,"ELASTIC_MECHANICS");
  Mat ShellAij;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,(void*)MatCtx,&ShellAij); CHKERRQ(ierr);
  ierr = MatShellSetOperation(ShellAij,MATOP_MULT,(void(*)(void))arc_length_mult_shell); CHKERRQ(ierr);

  //body forces
  BodyFroceConstantField body_forces_methods(mField);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|Block_BodyForcesSet,it)) {
    ierr = body_forces_methods.addBlock("DISPLACEMENT",F_body_force,it->get_msId()); CHKERRQ(ierr);
  }
  ierr = VecZeroEntries(F_body_force); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_body_force,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_body_force,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","BODY_FORCE",body_forces_methods.getLoopFe()); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_body_force,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_body_force,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F_body_force); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F_body_force); CHKERRQ(ierr);

  //surface forces
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|ForceSet,it)) {
    ostringstream fe_name;
    fe_name << "FORCE_FE_" << it->get_msId();
    string fe_name_str = fe_name.str();
    neumann_forces.insert(fe_name_str,new NeummanForcesSurface(mField));
    ierr = neumann_forces.at(fe_name_str).addForce("DISPLACEMENT",ArcCtx->F_lambda,it->get_msId());  CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|PressureSet,it)) {
    ostringstream fe_name;
    fe_name << "PRESSURE_FE_" << it->get_msId();
    string fe_name_str = fe_name.str();
    neumann_forces.insert(fe_name_str,new NeummanForcesSurface(mField));
    ierr = neumann_forces.at(fe_name_str).addPreassure("DISPLACEMENT",ArcCtx->F_lambda,it->get_msId()); CHKERRQ(ierr);
  }
  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
  ierr = VecZeroEntries(ArcCtx->F_lambda); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ArcCtx->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ArcCtx->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  for(;mit!=neumann_forces.end();mit++) {
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
  }
  ierr = VecGhostUpdateBegin(ArcCtx->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ArcCtx->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ArcCtx->F_lambda); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ArcCtx->F_lambda); CHKERRQ(ierr);

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
  SnesCtx::loops_to_do_type& loops_to_do_Rhs = SnesCtx.get_loops_to_do_Rhs();
  SnesCtx.get_preProcess_to_do_Rhs().push_back(&PrePostFE);
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ELASTIC",&MyFE));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("INTERFACE",&IntMyFE));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&MyArcMethod));
  SnesCtx.get_postProcess_to_do_Rhs().push_back(&PrePostFE);

  //Mat
  SnesCtx::loops_to_do_type& loops_to_do_Mat = SnesCtx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ELASTIC",&MyFE));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("INTERFACE",&IntMyFE));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&MyArcMethod));

  int its_d = 6;
  double gamma = 0.5,reduction = 1;
  double step_size0;
  //step = 1;
  if(step == 1) {
    step_size = step_size_reduction;
  } else {
    reduction = step_size;
    step_size0 = step_size_reduction;
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
  }

  for(;step<max_steps;step++) {

    if(step == 1) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Load Step %D step_size = %6.4e\n",step,step_size); CHKERRQ(ierr);
      ierr = ArcCtx->set_s(step_size); CHKERRQ(ierr);
      ierr = ArcCtx->set_alpha_and_beta(0,1); CHKERRQ(ierr);
      ierr = VecCopy(D,ArcCtx->x0); CHKERRQ(ierr);
      double dlambda;
      ierr = MyArcMethod.calculate_init_dlambda(&dlambda); CHKERRQ(ierr);
      ierr = MyArcMethod.set_dlambda_to_x(D,dlambda); CHKERRQ(ierr);
    } else if(step == 2) {
      ierr = ArcCtx->set_alpha_and_beta(1,0); CHKERRQ(ierr);
      ierr = MyArcMethod.calulate_dx_and_dlambda(D); CHKERRQ(ierr);
      ierr = MyArcMethod.calulate_lambda_int(step_size); CHKERRQ(ierr);
      step_size0 = step_size;
      ierr = ArcCtx->set_s(step_size); CHKERRQ(ierr);
      double dlambda = ArcCtx->dlambda;
      double dx_nrm;
      ierr = VecNorm(ArcCtx->dx,NORM_2,&dx_nrm);  CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,
	"Load Step %D step_size = %6.4e dlambda0 = %6.4e dx_nrm = %6.4e dx2 = %6.4e\n",
	step,step_size,dlambda,dx_nrm,ArcCtx->dx2); CHKERRQ(ierr);
      ierr = VecCopy(D,ArcCtx->x0); CHKERRQ(ierr);
      ierr = VecAXPY(D,1.,ArcCtx->dx); CHKERRQ(ierr);
      ierr = MyArcMethod.set_dlambda_to_x(D,dlambda); CHKERRQ(ierr);
    } else {
      ierr = MyArcMethod.calulate_dx_and_dlambda(D); CHKERRQ(ierr);
      double step_size1 = step_size;
      ierr = MyArcMethod.calulate_lambda_int(step_size); CHKERRQ(ierr);
      //step_size0_1/step_size0 = step_stize1/step_size
      //step_size0_1 = step_size0*(step_stize1/step_size)
      step_size0 = step_size0*(step_size1/step_size);
      step_size *= reduction;
      ierr = ArcCtx->set_s(step_size); CHKERRQ(ierr);
      double dlambda = reduction*ArcCtx->dlambda;
      double dx_nrm;
      ierr = VecScale(ArcCtx->dx,reduction); CHKERRQ(ierr);
      ierr = VecNorm(ArcCtx->dx,NORM_2,&dx_nrm);  CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,
	"Load Step %D step_size = %6.4e dlambda0 = %6.4e dx_nrm = %6.4e dx2 = %6.4e\n",
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
    IntMyFE.snes_ctx = FieldInterface::FEMethod::ctx_SNESNone;
    ierr = IntMyFE.set_ctx_int(NonLinearInterfaceFEMethod::ctx_KappaUpdate); CHKERRQ(ierr);
    //run this on all processors, so we could save history tags on all parts and restart
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",IntMyFE,0,pcomm->size());  CHKERRQ(ierr);
    //Standard procedure
    ierr = IntMyFE.set_ctx_int(NonLinearInterfaceFEMethod::ctx_InterfaceNone); CHKERRQ(ierr);
    //Remove nodes of damaged prisms
    ierr = MyArcMethod.remove_damaged_prisms_nodes(); CHKERRQ(ierr);

    int its;
    ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);

    if(step > 1) {
      if(its>0) {
	reduction = pow((double)its_d/(double)(its+1),gamma);
	if(step_size*reduction > 10*step_size0) {
	  reduction = 1;
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"reduction step_size (step_size,step_size0) = %6.4e (%6.4e,%6.4e)\n",reduction,step_size,step_size0); CHKERRQ(ierr);
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
	EntityHandle out_meshset;
	rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
	ostringstream ss;
	ss << "out_" << step << ".vtk";
	rval = moab.write_file(ss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
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

  PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab,"DISPLACEMENT");
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&F_body_force); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);
  ierr = MatDestroy(&ShellAij); CHKERRQ(ierr);
  delete ArcCtx;
  delete MatCtx;
  delete PCCtx;
  delete MyArcMethod_ptr;

  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

}

