/* Copyright (C) 2015, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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
#include <Projection10NodeCoordsOnField.hpp>

#include <SurfacePressure.hpp>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>


#include <ElasticFEMethod.hpp>  // Old implementaiton by Lukasz
#include "ElasticFE_RVELagrange_Disp.hpp"
#include "ElasticFE_RVELagrange_Disp_Multi_Rhs.hpp"
#include "RVEVolume.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Disp.hpp"





#include <ElasticElement.hpp>  // New implementation (in progress)
using namespace ObosleteUsersModules;
#include <Calculate_RVE_Dmat.hpp>



#include <PotsProcOnRefMesh.hpp>


using namespace boost::numeric;


static char help[] = "...\n\n";

int main(int argc, char *argv[]) {
  ErrorCode rval;
  PetscErrorCode ierr;
  PetscInitialize(&argc,&argv,(char *)0,help);
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscBool flg = PETSC_TRUE;
  const char *option;
  PetscInt order;

  //====================================================================================================
  //  DEFINING RVE PROBLEM
  //====================================================================================================
  
  moab::Core mb_instance_RVE;
  Interface& moab_RVE = mb_instance_RVE;
  
  char mesh_file_name_RVE[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_RVE",mesh_file_name_RVE,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_RVE (MESH FILE NEEDED)");
  }
  ParallelComm* pcomm_RVE = ParallelComm::get_pcomm(&moab_RVE,MYPCOMM_INDEX);
  if(pcomm_RVE == NULL) pcomm_RVE =  new ParallelComm(&moab_RVE,PETSC_COMM_WORLD);
  
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab_RVE.load_file(mesh_file_name_RVE, 0, option); CHKERR_PETSC(rval);
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core_RVE(moab_RVE);
  FieldInterface& m_field_RVE = core_RVE;
  
  //set entitities bit level
  BitRefLevel bit_level0_RVE;
  bit_level0_RVE.set(0);
  EntityHandle meshset_level0_RVE;
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_level0_RVE); CHKERR_PETSC(rval);
  ierr = m_field_RVE.seed_ref_level_3D(0,bit_level0_RVE); CHKERRQ(ierr);
  
  //Fields
  int field_rank=3;
  ierr = m_field_RVE.add_field("DISP_RVE",H1,3); CHKERRQ(ierr);
  ierr = m_field_RVE.add_field("Lagrange_mul_disp",H1,field_rank); CHKERRQ(ierr);
  
  //FE
  ierr = m_field_RVE.add_finite_element("ELASTIC_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.add_finite_element("Lagrange_FE"); CHKERRQ(ierr);
  
  //Define rows/cols and element data
  ierr = m_field_RVE.modify_finite_element_add_field_row("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  
  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = m_field_RVE.modify_finite_element_add_field_row("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  
  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = m_field_RVE.modify_finite_element_add_field_col("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_row("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  
  //data
  ierr = m_field_RVE.modify_finite_element_add_field_data("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  
  //define problems
  ierr = m_field_RVE.add_problem("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  
  //set finite elements for problem
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","Lagrange_FE"); CHKERRQ(ierr);
  
  
  //set refinment level for problem
  ierr = m_field_RVE.modify_problem_ref_level_add_bit("ELASTIC_PROBLEM_RVE",bit_level0_RVE); CHKERRQ(ierr);
  
  /***/
  //Declare problem
  
  //add entitities (by tets) to the field
  ierr = m_field_RVE.add_ents_to_field_by_TETs(0,"DISP_RVE"); CHKERRQ(ierr);
  
  
  //add finite elements entities
  ierr = m_field_RVE.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0_RVE,"ELASTIC_FE_RVE",MBTET); CHKERRQ(ierr);
  Range SurfacesFaces;
  ierr = m_field_RVE.get_Cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_FE"); CHKERRQ(ierr);
  
  
  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab_RVE.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
	rval = moab_RVE.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = m_field_RVE.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Lagrange_mul_disp",2); CHKERRQ(ierr);
  
  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  ierr = m_field_RVE.set_field_order(0,MBTET,"DISP_RVE",order); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBTRI,"DISP_RVE",order); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,"DISP_RVE",order); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"DISP_RVE",1); CHKERRQ(ierr);
  
  ierr = m_field_RVE.set_field_order(0,MBTRI,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);
  
  /****/
  //build database
  
  //build field
  ierr = m_field_RVE.build_fields(); CHKERRQ(ierr);
  
  //build finite elemnts
  ierr = m_field_RVE.build_finite_elements(); CHKERRQ(ierr);
  
  //build adjacencies
  ierr = m_field_RVE.build_adjacencies(bit_level0_RVE); CHKERRQ(ierr);
  
  //build problem
  ierr = m_field_RVE.build_problems(); CHKERRQ(ierr);
  
  
  /****/
  //mesh partitioning
  
  //partition
  ierr = m_field_RVE.partition_problem("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.partition_finite_elements("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field_RVE.partition_ghost_dofs("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  
  //print bcs
  ierr = m_field_RVE.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field_RVE.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = m_field_RVE.print_cubit_materials_set(); CHKERRQ(ierr);
  
  
  //====================================================================================================
  //  DEFINING MACRO PROBLEM
  //====================================================================================================

  moab::Core mb_instance_Macro;
  Interface& moab_Macro = mb_instance_Macro;
  
  char mesh_file_name_Macro[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_Macro",mesh_file_name_Macro,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_Macro (MESH FILE NEEDED)");
  }
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab_Macro,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab_Macro,PETSC_COMM_WORLD);

  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab_Macro.load_file(mesh_file_name_Macro, 0, option); CHKERR_PETSC(rval);

  //Create MoFEM (Joseph) database
  MoFEM::Core core_Macro(moab_Macro);
  FieldInterface& m_field_Macro = core_Macro;

  //set entitities bit level
  BitRefLevel bit_level0_Macro;
  bit_level0_Macro.set(0);
  EntityHandle meshset_level0;
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field_Macro.seed_ref_level_3D(0,bit_level0_Macro); CHKERRQ(ierr);

  //Fields
  ierr = m_field_Macro.add_field("DISP_MACRO",H1,3); CHKERRQ(ierr);

  //Problem
  ierr = m_field_Macro.add_problem("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = m_field_Macro.modify_problem_ref_level_add_bit("ELASTIC_PROBLEM_MACRO",bit_level0_Macro); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab_Macro.get_root_set(); 
  //add entities to field
  ierr = m_field_Macro.add_ents_to_field_by_TETs(root_set,"DISP_MACRO"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }
  ierr = m_field_Macro.set_field_order(root_set,MBTET,"DISP_MACRO",order); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(root_set,MBTRI,"DISP_MACRO",order); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(root_set,MBEDGE,"DISP_MACRO",order); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(root_set,MBVERTEX,"DISP_MACRO",1); CHKERRQ(ierr);
  
  if(!(m_field_Macro.check_field("MESH_NODE_POSITIONS"))){
    ierr = m_field_Macro.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
    ierr = m_field_Macro.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  }

  ElasticElement elastic_elements(m_field_Macro);

  ierr = elastic_elements.addElasticElements("DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO","Wt"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO"); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = m_field_Macro.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field_Macro.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field_Macro.build_adjacencies(bit_level0_Macro); CHKERRQ(ierr);
  //build problem
  ierr = m_field_Macro.build_problems(); CHKERRQ(ierr);

  Projection10NodeCoordsOnField ent_method_material(m_field_Macro,"MESH_NODE_POSITIONS");
  ierr = m_field_Macro.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  /****/
  //mesh partitioning 
  //partition
  ierr = m_field_Macro.partition_problem("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.partition_finite_elements("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field_Macro.partition_ghost_dofs("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);

//  Vec F;
//  ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&F); CHKERRQ(ierr);
//  Vec D;
//  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
//  Mat A;
//  ierr = m_field_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_MACRO",&A); CHKERRQ(ierr);
//
//  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field_Macro,"ELASTIC_PROBLEM_MACRO",A,D,F);
//
//  ierr = VecZeroEntries(D); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecZeroEntries(F); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  
//  //External forces, This vector is assemble only once (as this is not function of Dmat)
//  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
//  ierr = MetaNeummanForces::setNeumannFiniteElementOperators(m_field_Macro,neumann_forces,F,"ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
//  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
//  for(;mit!=neumann_forces.end();mit++) {
//    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
//  }
  
  
//  Vec Fint;
//  ierr = VecDuplicate(F,&Fint); CHKERRQ(ierr);
//  ierr = elastic_elements.setElasticFiniteElementLhsOperators("ELASTIC_PROBLEM_MACRO","Wt",A); CHKERRQ(ierr);
//  ierr = elastic_elements.setElasticFiniteElementRhsOperators("ELASTIC_PROBLEM_MACRO","Wt",Fint); CHKERRQ(ierr);

  Calculate_RVE_Dmat calculate_rve_dmat(m_field_Macro);
  ierr = calculate_rve_dmat.setRVE_DmatRhsOperators(m_field_RVE, "DISP_MACRO","Wt"); CHKERRQ(ierr);

  
  SeriesRecorder *recorder_ptr;
  ierr = m_field_Macro.query_interface(recorder_ptr); CHKERRQ(ierr);
  if( recorder_ptr->check_series("Wt_SERIES") ) {
    cout<<"============== Wt_SERIES exists =============== "<<endl;
    for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder_ptr,"Wt_SERIES",sit)) {
      PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());
      ierr = recorder_ptr->load_series_data("Wt_SERIES",sit->get_step_number()); CHKERRQ(ierr);

//      //We need to assemble A matrix and internal force vector at each time step as these depends on Dmat, which will change at each time step
//      ierr = MatZeroEntries(A); CHKERRQ(ierr);
//      ierr = VecZeroEntries(Fint); CHKERRQ(ierr);
//      ierr = VecGhostUpdateBegin(Fint,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//      ierr = VecGhostUpdateEnd(Fint,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//
//      //preproc
//      ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
//      ierr = m_field_Macro.set_global_VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//      
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",calculate_rve_dmat.getLoopFeRhs()); CHKERRQ(ierr);
//
//
////      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE",elastic_elements.getLoopFeRhs()); CHKERRQ(ierr);
////      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE",elastic_elements.getLoopFeLhs()); CHKERRQ(ierr);
//
//      //postproc
//      ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);

    }
  }

//  //TS
//  TsCtx ts_ctx(m_field_Macro,"ELASTIC_PROBLEM_MACRO");
//  TS ts;
//  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
//  ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);
//
//  TemperatureBCFEMethodPreAndPostProc my_dirichlet_bc(m_field_Macro,"TEMP",A,T,F);
//  ThermalElement::UpdateAndControl update_velocities(m_field_Macro,"TEMP","TEMP_RATE");
//  ThermalElement::TimeSeriesMonitor monitor(m_field_Macro,"THEMP_SERIES","TEMP");
//
//  ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
//  ierr = m_field_Macro.set_global_VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//
//  //preprocess
//  ts_ctx.get_preProcess_to_do_IFunction().push_back(&update_velocities);
//  ts_ctx.get_preProcess_to_do_IFunction().push_back(&my_dirichlet_bc);
//  ts_ctx.get_preProcess_to_do_IJacobian().push_back(&my_dirichlet_bc);
//
//  //and temperature element functions
//  ierr = thermal_elements.setTimeSteppingProblem(ts_ctx,"TEMP","TEMP_RATE"); CHKERRQ(ierr);
//
//  //postprocess
//  ts_ctx.get_postProcess_to_do_IFunction().push_back(&my_dirichlet_bc);
//  ts_ctx.get_postProcess_to_do_IJacobian().push_back(&my_dirichlet_bc);
//  ts_ctx.get_postProcess_to_do_Monitor().push_back(&monitor);
//
//  ierr = TSSetIFunction(ts,F,f_TSSetIFunction,&ts_ctx); CHKERRQ(ierr);
//  ierr = TSSetIJacobian(ts,A,A,f_TSSetIJacobian,&ts_ctx); CHKERRQ(ierr);
//  ierr = TSMonitorSet(ts,f_TSMonitorSet,&ts_ctx,PETSC_NULL); CHKERRQ(ierr);
//
//  double ftime = 1;
//  ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
//  ierr = TSSetSolution(ts,T); CHKERRQ(ierr);
//  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);
//
//  SeriesRecorder *recorder_ptr;
//  ierr = m_field_Macro.query_interface(recorder_ptr); CHKERRQ(ierr);
//  ierr = recorder_ptr->add_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);
//  ierr = recorder_ptr->initialize_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);
//
//  ierr = TSSolve(ts,T); CHKERRQ(ierr);
//  ierr = TSGetTime(ts,&ftime); CHKERRQ(ierr);
//
//  ierr = recorder_ptr->finalize_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);
//
//  PetscInt steps,snesfails,rejects,nonlinits,linits;
//  ierr = TSGetTimeStepNumber(ts,&steps); CHKERRQ(ierr);
//  ierr = TSGetSNESFailures(ts,&snesfails); CHKERRQ(ierr);
//  ierr = TSGetStepRejections(ts,&rejects); CHKERRQ(ierr);
//  ierr = TSGetSNESIterations(ts,&nonlinits); CHKERRQ(ierr);
//  ierr = TSGetKSPIterations(ts,&linits); CHKERRQ(ierr);
//
//  PetscPrintf(PETSC_COMM_WORLD,
//    "steps %D (%D rejected, %D SNES fails), ftime %g, nonlinits %D, linits %D\n",
//    steps,rejects,snesfails,ftime,nonlinits,linits);
//
//  //m_field_Macro.list_dofs_by_field_name("TEMP");
//  if(pcomm->rank()==0) {
//    rval = moab_Macro.write_file("solution.h5m"); CHKERR_PETSC(rval);
//  }
//
//  for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder_ptr,"THEMP_SERIES",sit)) {
//
//    PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());
//
//    ierr = recorder_ptr->load_series_data("THEMP_SERIES",sit->get_step_number()); CHKERRQ(ierr);
//    ierr = m_field_Macro.set_local_VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//
//    ProjectionFieldOn10NodeTet ent_method_on_10nodeTet(m_field_Macro,"TEMP",true,false,"TEMP");
//    ent_method_on_10nodeTet.set_nodes = true;
//    ierr = m_field_Macro.loop_dofs("TEMP",ent_method_on_10nodeTet); CHKERRQ(ierr);
//    ent_method_on_10nodeTet.set_nodes = false;
//    ierr = m_field_Macro.loop_dofs("TEMP",ent_method_on_10nodeTet); CHKERRQ(ierr);
//
//    if(pcomm->rank()==0) {
//      EntityHandle out_meshset;
//      rval = moab_Macro.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
//      ierr = m_field_Macro.problem_get_FE("ELASTIC_PROBLEM_MACRO","THERMAL_FE",out_meshset); CHKERRQ(ierr);
//      ostringstream ss;
//      ss << "out_" << sit->step_number << ".vtk";
//      rval = moab_Macro.write_file(ss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
//      rval = moab_Macro.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
//    }
//
//  }
//
//  ierr = TSDestroy(&ts);CHKERRQ(ierr);
//  ierr = MatDestroy(&A); CHKERRQ(ierr);
//  ierr = VecDestroy(&F); CHKERRQ(ierr);
//  ierr = VecDestroy(&T); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


