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
#include "SurfacePressure.hpp"
#include "NodalForce.hpp"
#include "ElasticFEMethodInterface.hpp"

#include "FEMethod_DriverComplexForLazy.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcNonLinearElasticityStresseOnRefindedMesh.hpp"
#include "Projection10NodeCoordsOnField.hpp"

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
	
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 3;
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

  /*Range tets;
  rval = moab.get_entities_by_type(0,MBTET,tets,true);  CHKERR_PETSC(rval);
  for(Range::iterator tit = tets.begin();tit!=tets.end();tit++) {
    cout << "tet: " << *tit << " conn: ";
    Range conn;
    rval = moab.get_connectivity(&*tit,1,conn);  CHKERR_PETSC(rval);
    ostream_iterator<EntityHandle> out_it (std::cout,", ");
    copy ( conn.begin(), conn.end(), out_it ); 
    cout << endl;
  }*/

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

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);
  vector<BitRefLevel> bit_levels;
  bit_levels.push_back(BitRefLevel().set(0));
  BitRefLevel problem_bit_level;

  if(step == 1) {
    
    int ll = 1;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|INTERFACESET,cit)) {
      //for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,cit)) {
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
    
    problem_bit_level = bit_levels.back();
    
    Range CubitSideSets_meshsets;
    ierr = mField.get_Cubit_meshsets(SIDESET,CubitSideSets_meshsets); CHKERRQ(ierr);

    //Fields
    ierr = mField.add_field("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);
    ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);

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
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","LAMBDA"); CHKERRQ(ierr);

    //Define FE Interface
    ierr = mField.add_finite_element("INTERFACE"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("INTERFACE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("INTERFACE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("INTERFACE","SPATIAL_POSITION"); CHKERRQ(ierr);
    
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
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","INTERFACE"); CHKERRQ(ierr);

    //set refinment level for problem
    ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);

    //add entitities (by tets) to the field
    ierr = mField.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

    //add finite elements entities
    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"ELASTIC",MBTET); CHKERRQ(ierr);
    //add prisms to interface finite element
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
    ierr = mField.set_field_order(0,MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);
    //
    ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

    //add neumman finite elemnets to add static boundary conditions
    ierr = mField.add_finite_element("NEUAMNN_FE"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("NEUAMNN_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
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
    //add nodal force element
    ierr = MetaNodalForces::addNodalForceElement(mField,"ELASTIC_MECHANICS","SPATIAL_POSITION"); CHKERRQ(ierr);

  }

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finiteElementsPtr(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(problem_bit_level); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  if(step==1) {
    //10 node tets
    Projection10NodeCoordsOnField ent_method_material(mField,"MESH_NODE_POSITIONS");
    ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
    ierr = mField.set_field(0,MBVERTEX,"SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.set_field(0,MBEDGE,"SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.field_axpy(1.,"MESH_NODE_POSITIONS","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.set_field(0,MBTRI,"SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.set_field(0,MBTET,"SPATIAL_POSITION"); CHKERRQ(ierr);
  }

  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finiteElementsPtr("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //print bcs
  ierr = mField.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = mField.print_cubit_pressure_set(); CHKERRQ(ierr);
  ierr = mField.print_cubit_force_set(); CHKERRQ(ierr);

  //print block sets with materials
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);

  //create matrices
  Vec F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",COL,&F); CHKERRQ(ierr);
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

  ArcLengthSnesCtx snes_ctx(mField,"ELASTIC_MECHANICS",arc_ctx);

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

  Range node_set;
  for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"LoadPath",cit)) {
    EntityHandle meshset = cit->get_meshset();
    Range nodes;
    rval = moab.get_entities_by_type(meshset,MBVERTEX,nodes,true); CHKERR_THROW(rval);
    node_set.merge(nodes);
  }
  PetscPrintf(PETSC_COMM_WORLD,"Nb. nodes in load path: %u\n",node_set.size());

  ArcLengthElemFEMethod* arc_method_ptr = new ArcLengthElemFEMethod(moab,arc_ctx);
  ArcLengthElemFEMethod& arc_method = *arc_method_ptr;

  NonLinearSpatialElasticFEMthod fe(mField,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));

  const double alpha = 0.1;
  InterfaceFEMethod int_fe(mField,Aij,D,F,young_modulus*alpha,"SPATIAL_POSITION");

  double scaled_reference_load = 1;
  double *scale_lhs = &(arc_ctx->get_FieldData());
  double *scale_rhs = &(scaled_reference_load);
  NeummanForcesSurfaceComplexForLazy neumann_forces(mField,Aij,arc_ctx->F_lambda,scale_lhs,scale_rhs);
  NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE &fe_neumann = neumann_forces.getLoopSpatialFe();
  fe_neumann.uSeF = true;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
    ierr = fe_neumann.addForce(it->get_msId()); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {
    ierr = fe_neumann.addPreassure(it->get_msId()); CHKERRQ(ierr);
  }
  SpatialPositionsBCFEMethodPreAndPostProc my_dirihlet_bc(mField,"SPATIAL_POSITION",Aij,D,F);
  ierr = mField.get_problem("ELASTIC_MECHANICS",&my_dirihlet_bc.problemPtr); CHKERRQ(ierr);
  ierr = my_dirihlet_bc.iNitalize(); CHKERRQ(ierr);

  struct MyPrePostProcessFEMethod: public FieldInterface::FEMethod {
    
    FieldInterface& mField;
    ArcLengthCtx *arc_ptr;

    SpatialPositionsBCFEMethodPreAndPostProc *bC;
    Range &nodeSet;

    MyPrePostProcessFEMethod(FieldInterface& _mField,
      ArcLengthCtx *_arc_ptr,SpatialPositionsBCFEMethodPreAndPostProc *bc,Range &node_set): 
      mField(_mField),arc_ptr(_arc_ptr),bC(bc),nodeSet(node_set) {}
  
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
        }
        break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode potsProcessLoadPath() {
      PetscFunctionBegin;
      NumeredDofMoFEMEntity_multiIndex &numered_dofs_rows = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problemPtr->numered_dofs_rows);
      Range::iterator nit = nodeSet.begin();
      for(;nit!=nodeSet.end();nit++) {
	NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator dit,hi_dit;
	dit = numered_dofs_rows.get<MoABEnt_mi_tag>().lower_bound(*nit);
	hi_dit = numered_dofs_rows.get<MoABEnt_mi_tag>().upper_bound(*nit);
	for(;dit!=hi_dit;dit++) {
	  PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e -> ","LAMBDA",0,arc_ptr->get_FieldData());
	  PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e\n",dit->get_name().c_str(),dit->get_dof_rank(),dit->get_FieldData());
	}
      }
      PetscFunctionReturn(0);
    }

  };

  struct AssembleLambdaFEMethod: public FieldInterface::FEMethod {
    
    FieldInterface& mField;
    ArcLengthCtx *arc_ptr;

    SpatialPositionsBCFEMethodPreAndPostProc *bC;

    AssembleLambdaFEMethod(FieldInterface& _mField,
      ArcLengthCtx *_arc_ptr,SpatialPositionsBCFEMethodPreAndPostProc *bc): 
      mField(_mField),arc_ptr(_arc_ptr),bC(bc) {}
  
    PetscErrorCode ierr;
      
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      switch(snes_ctx) {
        case CTX_SNESSETFUNCTION: {
	  //F_lambda
          ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	  for(vector<int>::iterator vit = bC->dofsIndices.begin();vit!=bC->dofsIndices.end();vit++) {
	    ierr = VecSetValue(arc_ptr->F_lambda,*vit,0,INSERT_VALUES); CHKERRQ(ierr);
	  }
	  ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	  ierr = VecDot(arc_ptr->F_lambda,arc_ptr->F_lambda,&arc_ptr->F_lambda2); CHKERRQ(ierr);
	  PetscPrintf(PETSC_COMM_WORLD,"\tFlambda2 = %6.4e\n",arc_ptr->F_lambda2);
	  //add F_lambda
	  ierr = VecAXPY(snes_f,arc_ptr->get_FieldData(),arc_ptr->F_lambda); CHKERRQ(ierr);
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

  MyPrePostProcessFEMethod pre_post_method(mField,arc_ctx,&my_dirihlet_bc,node_set);
  AssembleLambdaFEMethod assemble_F_lambda(mField,arc_ctx,&my_dirihlet_bc);
  
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

  //
  /*ierr = SNESSetType(snes,SNESSHELL); CHKERRQ(ierr);
  ierr = SNESShellSetContext(snes,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESShellSetSolve(snes,snes_apply_arc_length); CHKERRQ(ierr);*/
  //

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
  snes_ctx.get_preProcess_to_do_Rhs().push_back(&my_dirihlet_bc);
  snes_ctx.get_preProcess_to_do_Rhs().push_back(&pre_post_method);
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ELASTIC",&fe));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("INTERFACE",&int_fe));
  //surface focres and preassures
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_neumann));
  //nodal forces
  boost::ptr_map<string,NodalForce> nodal_forces;
  string fe_name_str ="FORCE_FE";
  nodal_forces.insert(fe_name_str,new NodalForce(mField));
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
    ierr = nodal_forces.at(fe_name_str).addForce("SPATIAL_POSITION",arc_ctx->F_lambda,it->get_msId());  CHKERRQ(ierr);
  }
  boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type(fit->first,&fit->second->getLoopFe()));
  }
  //arc length
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("NONE",&assemble_F_lambda));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&arc_method));
  snes_ctx.get_postProcess_to_do_Rhs().push_back(&pre_post_method);
  snes_ctx.get_postProcess_to_do_Rhs().push_back(&my_dirihlet_bc);

  SnesCtx::loops_to_do_type& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
  snes_ctx.get_preProcess_to_do_Mat().push_back(&my_dirihlet_bc);
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ELASTIC",&fe));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("INTERFACE",&int_fe));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_neumann));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&arc_method));
  snes_ctx.get_postProcess_to_do_Mat().push_back(&my_dirihlet_bc);

  ierr = mField.set_local_VecCreateGhost("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);


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


  if(step>1) {
    ierr = mField.set_other_global_VecCreateGhost(
      "ELASTIC_MECHANICS","SPATIAL_POSITION","X0_SPATIAL_POSITION",COL,arc_ctx->x0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    double x0_nrm;
    ierr = VecNorm(arc_ctx->x0,NORM_2,&x0_nrm);  CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\tRead x0_nrm = %6.4e dlambda = %6.4e\n",x0_nrm,arc_ctx->dlambda);
    ierr = arc_ctx->set_alpha_and_beta(1,0); CHKERRQ(ierr);
  } else {
    ierr = arc_ctx->set_s(step_size); CHKERRQ(ierr);
    ierr = arc_ctx->set_alpha_and_beta(0,1); CHKERRQ(ierr);
  }
  ierr = SnesRhs(snes,D,F,&snes_ctx); CHKERRQ(ierr);

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
    } else if(step == 2) {
      ierr = arc_ctx->set_alpha_and_beta(1,0); CHKERRQ(ierr);
      ierr = arc_method.calulate_dx_and_dlambda(D); CHKERRQ(ierr);
      step_size = sqrt(arc_method.calulate_lambda_int());
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
      ierr = arc_method.calulate_dx_and_dlambda(D); CHKERRQ(ierr);
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

    if(step % 1 == 0) {
      if(pcomm->rank()==0) {
	ostringstream sss;
	sss << "restart_" << step << ".h5m";
	rval = moab.write_file(sss.str().c_str()); CHKERR_PETSC(rval);
      }
      /*PostProcVertexMethod ent_method(moab,"SPATIAL_POSITION");
      ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",COL,ent_method); CHKERRQ(ierr);
      if(pcomm->rank()==0) {
	EntityHandle out_meshset;
	rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
	ostringstream sss;
	sss << "out_" << step << ".vtk";
	rval = moab.write_file(sss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      }*/
      PostProcStressNonLinearElasticity fe_post_proc_method(moab,fe);
      ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
      if(pcomm->rank()==0) {
	ostringstream sss;
	sss << "out_post_proc_" << step << ".vtk";
	rval = fe_post_proc_method.moab_post_proc.write_file(sss.str().c_str(),"VTK",""); CHKERR_PETSC(rval);
      }
    }

    ierr = pre_post_method.potsProcessLoadPath(); CHKERRQ(ierr);

  }

  ierr = VecDestroy(&D0); CHKERRQ(ierr);
  ierr = VecDestroy(&x00); CHKERRQ(ierr);
  
  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
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



