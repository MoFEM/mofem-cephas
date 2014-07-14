/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
 * FIXME: Interface element is implemented with assumptions of small
 * displacemennts. It has to be changed to new element taking in account large
 * displecements, i.e. change of interface normal.
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

#include "ElasticFEMethodInterface.hpp"

#include "SurfacePressureComplexForLazy.hpp"
#include "SurfacePressure.hpp"
#include "NodalForce.hpp"

#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include "PostProcNonLinearElasticityStresseOnRefindedMesh.hpp"

#include "FEMethod_DriverComplexForLazy.hpp"

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

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
 
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

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
    
  BitRefLevel problem_bit_level = bit_levels.back();

  Range CubitSIDESETs_meshsets;
  ierr = mField.get_Cubit_meshsets(SIDESET,CubitSIDESETs_meshsets); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_field("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);

  //Define FE Interface
  ierr = mField.add_finite_element("INTERFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("INTERFACE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("INTERFACE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("INTERFACE","SPATIAL_POSITION"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","INTERFACE"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"ELASTIC",MBTET); CHKERRQ(ierr);
  //add prisms to interface finite element
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"INTERFACE",MBPRISM); CHKERRQ(ierr);

  //set app. order
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 5;
  }
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
  ierr = MetaNodalForces::addNodalForceElement(mField,"ELASTIC_MECHANICS","SPATIAL_POSITION"); CHKERRQ(ierr);

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = mField.build_finiteElementsPtr(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(problem_bit_level); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finiteElementsPtr("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  ierr = mField.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = mField.print_cubit_pressure_set(); CHKERRQ(ierr);
  ierr = mField.print_cubit_force_set(); CHKERRQ(ierr);
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);

  //create matrices
  Vec F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",COL,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);
  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);

  {
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
  SnesCtx snes_ctx(mField,"ELASTIC_MECHANICS");

  const double young_modulus = 1;
  const double poisson_ratio = 0.25;
  NonLinearSpatialElasticFEMthod my_fe(mField,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));

  const double alpha = 0.1;
  InterfaceFEMethod int_fe(mField,Aij,D,F,young_modulus*alpha,"SPATIAL_POSITION");

  NeummanForcesSurfaceComplexForLazy neumann_forces(mField,Aij,F);
  NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE &fe_spatial = neumann_forces.getLoopSpatialFe();
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
    ierr = fe_spatial.addForce(it->get_msId()); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {
    ierr = fe_spatial.addPreassure(it->get_msId()); CHKERRQ(ierr);
  }
  SpatialPositionsBCFEMethodPreAndPostProc my_dirihlet_bc(mField,"SPATIAL_POSITION",Aij,D,F);
 
  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,Aij,Aij,SnesMat,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  snes_ctx.get_preProcess_to_do_Rhs().push_back(&my_dirihlet_bc);
  SnesCtx::loops_to_do_type& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ELASTIC",&my_fe));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("INTERFACE",&int_fe));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_spatial));
  //nodal forces
  boost::ptr_map<string,NodalForce> nodal_forces;
  string fe_name_str ="FORCE_FE";
  nodal_forces.insert(fe_name_str,new NodalForce(mField));
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
    ierr = nodal_forces.at(fe_name_str).addForce("SPATIAL_POSITION",F,it->get_msId(),true);  CHKERRQ(ierr);
    nodal_forces.at(fe_name_str).methodsOp.push_back(new MetaNodalForces::TagForceScale(mField));
  }
  boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type(fit->first,&fit->second->getLoopFe()));
  }
  //postporc, i.e. dirihlet bcs
  snes_ctx.get_postProcess_to_do_Rhs().push_back(&my_dirihlet_bc);

  snes_ctx.get_preProcess_to_do_Mat().push_back(&my_dirihlet_bc);
  SnesCtx::loops_to_do_type& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ELASTIC",&my_fe));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("INTERFACE",&int_fe));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_spatial));
  snes_ctx.get_postProcess_to_do_Mat().push_back(&my_dirihlet_bc);

  ierr = mField.set_local_VecCreateGhost("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  double step_size = -1e-3;
  for(int step = 1;step<4; step++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Load Step %D\n",step); CHKERRQ(ierr);
    ierr = neumann_forces.setForceScale(step_size*step); CHKERRQ(ierr);
    ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
    int its;
    ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  }
  
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PostProcVertexMethod ent_method(moab,"SPATIAL_POSITION");
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",COL,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PostProcFieldsAndGradientOnRefMesh fe_post_proc_method(moab);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  PostProcStressNonLinearElasticity fe_post_proc_stresses_method(moab,my_fe);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_stresses_method);  CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
    rval = fe_post_proc_stresses_method.moab_post_proc.write_file("out_post_proc_stresses.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

  return 0;
}



