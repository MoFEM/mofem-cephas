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

  //Problem
  ierr = mField.add_problem("THERMAL_PROBLEM"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("THERMAL_PROBLEM",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"TEMP"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 2;
  }

  ierr = mField.set_field_order(root_set,MBTET,"TEMP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"TEMP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"TEMP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"TEMP",1); CHKERRQ(ierr);

  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  ThermalElement thermal_elements(mField);
  ierr = thermal_elements.addThermalElements("THERMAL_PROBLEM","TEMP"); CHKERRQ(ierr);
  ierr = thermal_elements.addThermalFluxElement("THERMAL_PROBLEM","TEMP"); CHKERRQ(ierr);

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

  TemperatureBCFEMethodPreAndPostProc my_dirihlet_bc(mField,"TEMP",A,T,F);
  ierr = thermal_elements.setThermalFiniteElementRhsOperators("TEMP",F); CHKERRQ(ierr);
  ierr = thermal_elements.setThermalFiniteElementLhsOperators("TEMP",(&A)); CHKERRQ(ierr);
  ierr = thermal_elements.setThermalFluxFiniteElementLhsOperators("TEMP",F); CHKERRQ(ierr);

  ierr = VecZeroEntries(T); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);
  
  //preproc
  ierr = mField.problem_basic_method_preProcess("THERMAL_PROBLEM",my_dirihlet_bc); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost("THERMAL_PROBLEM",Row,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = mField.loop_finite_elements("THERMAL_PROBLEM","THERMAL_FE",thermal_elements.getLoopFeRhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("THERMAL_PROBLEM","THERMAL_FE",thermal_elements.getLoopFeLhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("THERMAL_PROBLEM","THERMAL_FLUX_FE",thermal_elements.getLoopFeFlux()); CHKERRQ(ierr);

  //postproc
  ierr = mField.problem_basic_method_postProcess("THERMAL_PROBLEM",my_dirihlet_bc); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = VecScale(F,-1); CHKERRQ(ierr);

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,T); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = mField.problem_basic_method_preProcess("THERMAL_PROBLEM",my_dirihlet_bc); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("THERMAL_PROBLEM",Row,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = mField.add_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);
  ierr = mField.initialize_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);

  ThermalElement::TimeSeriesMonitor monitor(mField,"THEMP_SERIES","TEMP");
  monitor.ts_u = T;
  ierr = mField.problem_basic_method_postProcess("THERMAL_PROBLEM",monitor); CHKERRQ(ierr);

  ierr = mField.finalize_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);

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

  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet(mField,"TEMP",true,false,"TEMP");
  ierr = mField.loop_dofs("TEMP",ent_method_on_10nodeTet); CHKERRQ(ierr);
  ent_method_on_10nodeTet.set_nodes = false;
  ierr = mField.loop_dofs("TEMP",ent_method_on_10nodeTet); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("THERMAL_PROBLEM","THERMAL_FE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  if(pcomm->rank()==0) {
    PostProcScalarFieldsAndGradientOnRefMesh fe_post_proc_method(moab,"TEMP");
    fe_post_proc_method.do_broadcast = false;
    ierr = mField.loop_finite_elements("THERMAL_PROBLEM","THERMAL_FE",fe_post_proc_method,0,pcomm->size());  CHKERRQ(ierr);
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&T); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


