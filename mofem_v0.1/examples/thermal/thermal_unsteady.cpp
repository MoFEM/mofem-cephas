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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "ThermalElement.hpp"
#include "DirihletBC.hpp"

#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include "Projection10NodeCoordsOnField.hpp"

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

#include <petscksp.h>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

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

  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_field("TEMP",H1,1); CHKERRQ(ierr);
  ierr = mField.add_field("TEMP_RATE",H1,1); CHKERRQ(ierr);

  //Problem
  ierr = mField.add_problem("THERMAL_PROBLEM"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("THERMAL_PROBLEM",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"TEMP"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"TEMP_RATE"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }
  ierr = mField.set_field_order(root_set,MBTET,"TEMP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"TEMP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"TEMP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"TEMP",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(root_set,MBTET,"TEMP_RATE",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"TEMP_RATE",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"TEMP_RATE",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"TEMP_RATE",1); CHKERRQ(ierr);

  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  ThermalElement thermal_elements(mField);
  ierr = thermal_elements.addThermalElements("THERMAL_PROBLEM","TEMP"); CHKERRQ(ierr);
  ierr = thermal_elements.addThermalFluxElement("THERMAL_PROBLEM","TEMP"); CHKERRQ(ierr);
  //add rate of temerature to data field of finite element
  ierr = mField.modify_finite_element_add_field_data("THERMAL_FE","TEMP_RATE"); CHKERRQ(ierr);

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

  Projection10NodeCoordsOnField ent_method_material(mField,"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  /****/
  //mesh partitioning 
  //partition
  ierr = mField.simple_partition_problem("THERMAL_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("THERMAL_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("THERMAL_PROBLEM"); CHKERRQ(ierr);

  Vec F;
  ierr = mField.VecCreateGhost("THERMAL_PROBLEM",Row,&F); CHKERRQ(ierr);
  Vec T;
  ierr = VecDuplicate(F,&T); CHKERRQ(ierr);
  Mat A;
  ierr = mField.MatCreateMPIAIJWithArrays("THERMAL_PROBLEM",&A); CHKERRQ(ierr);

  //TS
  TsCtx ts_ctx(mField,"THERMAL_PROBLEM");
  TS ts;
  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);

  TemperatureBCFEMethodPreAndPostProc my_dirihlet_bc(mField,"TEMP",A,T,F);
  ThermalElement::UpdateAndControl update_velocities(mField,ts,"TEMP","TEMP_RATE");
  ThermalElement::TimeSeriesMonitor monitor(mField,"THEMP_SERIES","TEMP");

  //preprocess
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&update_velocities);
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&my_dirihlet_bc);
  ts_ctx.get_preProcess_to_do_IJacobian().push_back(&my_dirihlet_bc);

  //and temperature element functions
  ierr = thermal_elements.setTimeSteppingProblem(ts_ctx,"TEMP","TEMP_RATE"); CHKERRQ(ierr);

  //postprocess
  ts_ctx.get_postProcess_to_do_IFunction().push_back(&my_dirihlet_bc);
  ts_ctx.get_postProcess_to_do_IJacobian().push_back(&my_dirihlet_bc);
  ts_ctx.get_postProcess_to_do_IJacobian().push_back(&update_velocities);
  ts_ctx.get_postProcess_to_do_Monitor().push_back(&monitor);

  ierr = TSSetIFunction(ts,F,f_TSSetIFunction,&ts_ctx); CHKERRQ(ierr);
  ierr = TSSetIJacobian(ts,A,A,f_TSSetIJacobian,&ts_ctx); CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,f_TSMonitorSet,&ts_ctx,PETSC_NULL); CHKERRQ(ierr);

  double ftime = 1;
  ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
  ierr = TSSetSolution(ts,T); CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);

  ierr = mField.add_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);
  ierr = mField.initialize_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);

  ierr = TSSolve(ts,T); CHKERRQ(ierr);
  ierr = TSGetTime(ts,&ftime); CHKERRQ(ierr);

  ierr = mField.finalize_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);

  PetscInt steps,snesfails,rejects,nonlinits,linits;
  ierr = TSGetTimeStepNumber(ts,&steps); CHKERRQ(ierr);
  ierr = TSGetSNESFailures(ts,&snesfails); CHKERRQ(ierr);
  ierr = TSGetStepRejections(ts,&rejects); CHKERRQ(ierr);
  ierr = TSGetSNESIterations(ts,&nonlinits); CHKERRQ(ierr);
  ierr = TSGetKSPIterations(ts,&linits); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,
    "steps %D (%D rejected, %D SNES fails), ftime %G, nonlinits %D, linits %D\n",
    steps,rejects,snesfails,ftime,nonlinits,linits);

  if(pcomm->rank()==0) {
    rval = moab.write_file("solution.h5m"); CHKERR_PETSC(rval);
  }
  
  /*EntityHandle fe_meshset = mField.get_finite_element_meshset("THERMAL_FE");
  Range tets;
  rval = moab.get_entities_by_type(fe_meshset,MBTET,tets,true); CHKERR_PETSC(rval);
  Range tets_edges;
  rval = moab.get_adjacencies(tets,1,false,tets_edges,Interface::UNION); CHKERR(rval);
  EntityHandle edges_meshset;
  rval = moab.create_meshset(MESHSET_SET,edges_meshset); CHKERR_PETSC(rval);
  rval = moab.add_entities(edges_meshset,tets); CHKERR_PETSC(rval);
  rval = moab.add_entities(edges_meshset,tets_edges); CHKERR_PETSC(rval);
  rval = moab.convert_entities(edges_meshset,true,false,false); CHKERR_PETSC(rval);*/

  for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(mField,"THEMP_SERIES",sit)) {

    PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());

    ierr = mField.load_series_data("THEMP_SERIES",sit->get_step_number()); CHKERRQ(ierr);
    ierr = mField.set_local_VecCreateGhost("THERMAL_PROBLEM",Row,T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    ProjectionFieldOn10NodeTet ent_method_on_10nodeTet(mField,"TEMP",true,false,"TEMP");
    ent_method_on_10nodeTet.set_nodes = true;
    ierr = mField.loop_dofs("TEMP",ent_method_on_10nodeTet); CHKERRQ(ierr);
    ent_method_on_10nodeTet.set_nodes = false;
    ierr = mField.loop_dofs("TEMP",ent_method_on_10nodeTet); CHKERRQ(ierr);

    if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      ierr = mField.problem_get_FE("THERMAL_PROBLEM","THERMAL_FE",out_meshset); CHKERRQ(ierr);
      ostringstream ss;
      ss << "out_" << sit->step_number << ".vtk";
      rval = moab.write_file(ss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }

    if(pcomm->rank()==0) {
      PostProcScalarFieldsAndGradientOnRefMesh fe_post_proc_method(moab,"TEMP");
      fe_post_proc_method.do_broadcast = false;
      ostringstream ss;
      ss << "out_post_proc_" << sit->step_number << ".vtk";
      ierr = mField.loop_finite_elements("THERMAL_PROBLEM","THERMAL_FE",fe_post_proc_method,0,pcomm->size());  CHKERRQ(ierr);
      rval = fe_post_proc_method.moab_post_proc.write_file(ss.str().c_str(),"VTK",""); CHKERR_PETSC(rval);
    }

  }

  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&T); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


