/* Copyright (C) 2015, Lukasz Kaczmarczyk (likask AT wp.pl)
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

#include <MoFEM.hpp>

#include <DirichletBC.hpp>
#include <PotsProcOnRefMesh.hpp>
#include <ThermalElement.hpp>

#include <Projection10NodeCoordsOnField.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  const char *option;
  option = "";

  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }

  DMType dm_name = "DMTHERMAL";
  ierr = DMRegister_MoFEM(dm_name); CHKERRQ(ierr);
  //craete dm instance
  DM dm;
  ierr = DMCreate(PETSC_COMM_WORLD,&dm);CHKERRQ(ierr);
  ierr = DMSetType(dm,dm_name);CHKERRQ(ierr);

  //create MoAB
  moab::Core mb_instance;
  Interface& moab = mb_instance;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  //create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  //set entitities bit level (this allow to set refinment levels for h-adaptivity)
  //onlt one level is used in this example
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields H1 space rank 1
  ierr = m_field.add_field("TEMP",H1,1); CHKERRQ(ierr);
  ierr = m_field.add_field("TEMP_RATE",H1,1); CHKERRQ(ierr);

  //Add field H1 space rank 3 to approximate gemetry using heierachical basis
  //For 10 node tets, before use, gemetry is projected on that field (see below)
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field (root_mesh, i.e. on all mesh etities fields are approx.)
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"TEMP"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"TEMP_RATE"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //for simplicity of example to all entities is applied the same order
  ierr = m_field.set_field_order(root_set,MBTET,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"TEMP",1); CHKERRQ(ierr);

  ierr = m_field.set_field_order(root_set,MBTET,"TEMP_RATE",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"TEMP_RATE",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"TEMP_RATE",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"TEMP_RATE",1); CHKERRQ(ierr);

  //gemetry approximation is set to 2nd oreder
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  //this defuault class to calulate thermal elements
  ThermalElement thermal_elements(m_field);
  ierr = thermal_elements.addThermalElements("TEMP"); CHKERRQ(ierr);
  ierr = thermal_elements.addThermalFluxElement("TEMP"); CHKERRQ(ierr);
  ierr = thermal_elements.addThermalConvectionElement("TEMP"); CHKERRQ(ierr);
  ierr = thermal_elements.addThermalRadiationElement("TEMP"); CHKERRQ(ierr);
  //add rate of temerature to data field of finite element
  ierr = m_field.modify_finite_element_add_field_data("THERMAL_FE","TEMP_RATE"); CHKERRQ(ierr);
  //and temperature element default element operators at integration (gauss) points
  ierr = thermal_elements.setTimeSteppingProblem("TEMP","TEMP_RATE"); CHKERRQ(ierr);

  //build database, i.e. declare dofs, elements and ajacencies

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //priject 10 node tet approximation of gemetry on hierarhical basis 
  Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //set dm datastruture whict created mofem datastructures
  ierr = DMMoFEMCreateMoFEM(dm,&m_field,dm_name,bit_level0); CHKERRQ(ierr);
  ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
  //add elements to dm
  ierr = DMMoFEMAddElement(dm,"THERMAL_FE"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"THERMAL_FLUX_FE"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"THERMAL_CONVECTION_FE"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"THERMAL_RADIATION_FE"); CHKERRQ(ierr);
  ierr = DMSetUp(dm); CHKERRQ(ierr);

  TemperatureBCFEMethodPreAndPostProc dirichlet_bc(m_field,"TEMP",PETSC_NULL,PETSC_NULL,PETSC_NULL);
  ThermalElement::UpdateAndControl update_velocities(m_field,"TEMP","TEMP_RATE");
  ThermalElement::TimeSeriesMonitor monitor(m_field,"THEMP_SERIES","TEMP");

  //preprocess
  ierr = DMMoFEMTSSetIFunction(dm,DM_NO_ELEMENT,NULL,&update_velocities,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIFunction(dm,DM_NO_ELEMENT,NULL,&dirichlet_bc,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIJacobian(dm,DM_NO_ELEMENT,NULL,&dirichlet_bc,NULL); CHKERRQ(ierr);

  //loops rhs
  ierr = DMMoFEMTSSetIFunction(dm,"THERMAL_FE",&thermal_elements.feRhs,NULL,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIFunction(dm,"THERMAL_FLUX_FE",&thermal_elements.feFlux,NULL,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIFunction(dm,"THERMAL_CONVECTION_FE",&thermal_elements.feConvectionRhs,NULL,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIFunction(dm,"THERMAL_RADIATION_FE",&thermal_elements.feRadiationRhs,NULL,NULL); CHKERRQ(ierr);

  //loops lhs
  ierr = DMMoFEMTSSetIJacobian(dm,"THERMAL_FE",&thermal_elements.feLhs,NULL,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIJacobian(dm,"THERMAL_CONVECTION_FE",&thermal_elements.feConvectionLhs,NULL,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIJacobian(dm,"THERMAL_RADIATION_FE",&thermal_elements.feRadiationLhs,NULL,NULL); CHKERRQ(ierr);

  //postprocess
  ierr = DMMoFEMTSSetIFunction(dm,DM_NO_ELEMENT,NULL,NULL,&dirichlet_bc); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIJacobian(dm,DM_NO_ELEMENT,NULL,NULL,&dirichlet_bc); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIJacobian(dm,DM_NO_ELEMENT,NULL,NULL,&update_velocities); CHKERRQ(ierr);

  TsCtx *ts_ctx;
  DMMoFEMGetTsCtx(dm,&ts_ctx);
  //add monitor opetator
  ts_ctx->get_postProcess_to_do_Monitor().push_back(&monitor);

  //create matrices
  Vec T,F;
  ierr = DMCreateGlobalVector_MoFEM(dm,&T); CHKERRQ(ierr);
  ierr = VecDuplicate(T,&F); CHKERRQ(ierr);
  Mat A;
  ierr = DMCreateMatrix_MoFEM(dm,&A); CHKERRQ(ierr);

  //create time solver
  TS ts;
  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);

  ierr = TSSetIFunction(ts,F,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
  ierr = TSSetIJacobian(ts,A,A,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
  //add monitor to TS solver
  ierr = TSMonitorSet(ts,f_TSMonitorSet,ts_ctx,PETSC_NULL); CHKERRQ(ierr); // !!!

  SeriesRecorder *recorder_ptr;
  ierr = m_field.query_interface(recorder_ptr); CHKERRQ(ierr);
  ierr = recorder_ptr->add_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);
  ierr = recorder_ptr->initialize_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);

  double ftime = 1;
  ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);
  ierr = TSSetDM(ts,dm); CHKERRQ(ierr);

  ierr = TSSolve(ts,T); CHKERRQ(ierr);
  ierr = TSGetTime(ts,&ftime); CHKERRQ(ierr);

  ierr = recorder_ptr->finalize_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);

  PetscInt steps,snesfails,rejects,nonlinits,linits;
  ierr = TSGetTimeStepNumber(ts,&steps); CHKERRQ(ierr);
  ierr = TSGetSNESFailures(ts,&snesfails); CHKERRQ(ierr);
  ierr = TSGetStepRejections(ts,&rejects); CHKERRQ(ierr);
  ierr = TSGetSNESIterations(ts,&nonlinits); CHKERRQ(ierr);
  ierr = TSGetKSPIterations(ts,&linits); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,
    "steps %D (%D rejected, %D SNES fails), ftime %g, nonlinits %D, linits %D\n",
    steps,rejects,snesfails,ftime,nonlinits,linits);
  
  //m_field.list_dofs_by_field_name("TEMP");
  if(pcomm->rank()==0) {
    rval = moab.write_file("solution.h5m"); CHKERR_PETSC(rval);
  }

  for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder_ptr,"THEMP_SERIES",sit)) {

    PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());

    ierr = recorder_ptr->load_series_data("THEMP_SERIES",sit->get_step_number()); CHKERRQ(ierr);
    ierr = m_field.set_local_VecCreateGhost(dm_name,ROW,T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    ProjectionFieldOn10NodeTet ent_method_on_10nodeTet(m_field,"TEMP",true,false,"TEMP");
    ent_method_on_10nodeTet.set_nodes = true;
    ierr = m_field.loop_dofs("TEMP",ent_method_on_10nodeTet); CHKERRQ(ierr);
    ent_method_on_10nodeTet.set_nodes = false;
    ierr = m_field.loop_dofs("TEMP",ent_method_on_10nodeTet); CHKERRQ(ierr);

    if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      ierr = m_field.problem_get_FE(dm_name,"THERMAL_FE",out_meshset); CHKERRQ(ierr);
      ostringstream ss;
      ss << "out_" << sit->step_number << ".vtk";
      rval = moab.write_file(ss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }

  }

  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&T); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);

  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


