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

#include <MoFEM.hpp>
using namespace MoFEM;

#include <DirichletBC.hpp>
#include <PostProcOnRefMesh.hpp>
#include <calculate_wt_with_derivative.hpp>

#include <Projection10NodeCoordsOnField.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;


static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
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
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("Wt",H1,1); CHKERRQ(ierr);
  ierr = m_field.add_field("Wt_RATE",H1,1); CHKERRQ(ierr);
  ierr = m_field.add_field("dWt_dBeta",H1,1); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("Wt_PROBLEM"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("Wt_PROBLEM",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"Wt"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"Wt_RATE"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"dWt_dBeta"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }
  ierr = m_field.set_field_order(root_set,MBTET,"Wt",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"Wt",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"Wt",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"Wt",1); CHKERRQ(ierr);

  ierr = m_field.set_field_order(root_set,MBTET,"Wt_RATE",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"Wt_RATE",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"Wt_RATE",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"Wt_RATE",1); CHKERRQ(ierr);
  
  ierr = m_field.set_field_order(root_set,MBTET,"dWt_dBeta",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"dWt_dBeta",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"dWt_dBeta",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"dWt_dBeta",1); CHKERRQ(ierr);

  Calculate_wt_with_derivative wt_elements(m_field);
  if(m_field.check_field("TEMP")) {
    cout<<"Temprature field exists "<< endl;
    if(m_field.check_field("CONC")) {
      cout<<"Concentration field exists "<< endl;
      ierr = wt_elements.addWtElement("Wt_PROBLEM","Wt_FE","Wt","Wt_RATE","TEMP","CONC","dWt_dBeta"); CHKERRQ(ierr);
    }
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
//
//  Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
//  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  /****/
  //mesh partitioning 
  //partition
  ierr = m_field.partition_problem("Wt_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("Wt_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("Wt_PROBLEM"); CHKERRQ(ierr);

  Vec F;
  ierr = m_field.VecCreateGhost("Wt_PROBLEM",ROW,&F); CHKERRQ(ierr);
  Vec T;
  ierr = VecDuplicate(F,&T); CHKERRQ(ierr);
  Mat A;
  ierr = m_field.MatCreateMPIAIJWithArrays("Wt_PROBLEM",&A); CHKERRQ(ierr);

  ierr = VecZeroEntries(T); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  
  
  //Setting initial conditions for Wt (Degradaiton parameter) which is 1 at each node (for both field (Wt) and vector (T))
  for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(m_field,"Wt",MBVERTEX,dof)) {
//    EntityHandle ent = dof->get_ent();
//    ublas::vector<double> coords(3);
//    rval = moab.get_coords(&ent,1,&coords[0]); CHKERR_PETSC(rval);
//    cout<<"coords "<<coords<<endl;
    dof->get_FieldData() = 1;
  }
  ierr = m_field.set_local_ghost_vector("Wt_PROBLEM",COL,T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

//  ierr = VecView(T,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  std::string wait;
//  std::cin >> wait;

  
  
  //TS
  TsCtx ts_ctx(m_field,"Wt_PROBLEM");
  TS ts;
  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);
  
  Calculate_wt_with_derivative::LoadTimeSeries load_series_data(m_field,"THEMP_SERIES","CONC_SERIES");
  Calculate_wt_with_derivative::UpdateAndControl update_velocities(m_field,ts,"Wt","Wt_RATE");
  //Calculate_wt_with_derivative::UpdateAndControl update_beta(m_field,ts,"Wt","dWt_dBeta");
  Calculate_wt_with_derivative::TimeSeriesMonitor monitor(m_field,"Wt_SERIES","Wt");

  //preprocess
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&load_series_data);
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&update_velocities);
  //ts_ctx.get_preProcess_to_do_IFunction().push_back(&update_beta);

  //and temperature element functions
  ierr = wt_elements.setTimeSteppingProblem(ts_ctx,"Wt_FE","Wt","Wt_RATE","TEMP","CONC","dWt_dBeta"); CHKERRQ(ierr);

  //postprocess
  ts_ctx.get_postProcess_to_do_IJacobian().push_back(&update_velocities);
  ts_ctx.get_postProcess_to_do_Monitor().push_back(&monitor);

  ierr = TSSetIFunction(ts,F,f_TSSetIFunction,&ts_ctx); CHKERRQ(ierr);
  ierr = TSSetIJacobian(ts,A,A,f_TSSetIJacobian,&ts_ctx); CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,f_TSMonitorSet,&ts_ctx,PETSC_NULL); CHKERRQ(ierr);

  double ftime = 1;
  cout<<"\n\nStart to calculate\n\n"<<endl;
  ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
  ierr = TSSetSolution(ts,T); CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);

  SeriesRecorder *recorder_ptr;
  ierr = m_field.query_interface(recorder_ptr); CHKERRQ(ierr);
  ierr = recorder_ptr->add_series_recorder("Wt_SERIES"); CHKERRQ(ierr);
  ierr = recorder_ptr->initialize_series_recorder("Wt_SERIES"); CHKERRQ(ierr);
  cout<<"\n\nstep time "<<ftime<<endl;
  ierr = TSSolve(ts,T); CHKERRQ(ierr);
  ierr = TSGetTime(ts,&ftime); CHKERRQ(ierr);
  cout<<"\n\nstep time "<<ftime<<endl;
  ierr = recorder_ptr->finalize_series_recorder("Wt_SERIES"); CHKERRQ(ierr);
  
  PetscInt steps,snesfails,rejects,nonlinits,linits;
  ierr = TSGetTimeStepNumber(ts,&steps); CHKERRQ(ierr);
  ierr = TSGetSNESFailures(ts,&snesfails); CHKERRQ(ierr);
  ierr = TSGetStepRejections(ts,&rejects); CHKERRQ(ierr);
  ierr = TSGetSNESIterations(ts,&nonlinits); CHKERRQ(ierr);
  ierr = TSGetKSPIterations(ts,&linits); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,
  "steps %D (%D rejected, %D SNES fails), ftime %g, nonlinits %D, linits %D\n",
  steps,rejects,snesfails,ftime,nonlinits,linits);

  
  //m_field.list_dofs_by_field_name("Wt");
  if(pcomm->rank()==0) {
    rval = moab.write_file("solution_wt.h5m"); CHKERR_PETSC(rval);
  }

  
  for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder_ptr,"Wt_SERIES",sit)) {

    PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());

    //ierr = recorder_ptr->load_series_data("Wt_SERIES",sit->get_step_number()); CHKERRQ(ierr);
    ierr = recorder_ptr->load_series_data("THEMP_SERIES",sit->get_step_number()); CHKERRQ(ierr);
    
    ierr = m_field.set_local_ghost_vector("Wt_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
ierr = VecView(T,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ProjectionFieldOn10NodeTet ent_method_on_10nodeTet(m_field,"Wt",true,false,"Wt");
    ent_method_on_10nodeTet.set_nodes = true;
    ierr = m_field.loop_dofs("Wt",ent_method_on_10nodeTet); CHKERRQ(ierr);
    ent_method_on_10nodeTet.set_nodes = false;
    ierr = m_field.loop_dofs("Wt",ent_method_on_10nodeTet); CHKERRQ(ierr);

    if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      ierr = m_field.get_problem_finite_elements_entities("Wt_PROBLEM","Wt_FE",out_meshset); CHKERRQ(ierr);
      ostringstream ss;
      ss << "Wt_" << sit->step_number << ".vtk";
      rval = moab.write_file(ss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }

  }

  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&T); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


