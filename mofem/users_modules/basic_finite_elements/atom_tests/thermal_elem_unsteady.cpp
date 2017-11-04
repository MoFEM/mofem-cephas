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

#include <BasicFiniteElements.hpp>
using namespace MoFEM;

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  
  

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

    moab::Core mb_instance;
    moab::Interface& moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
    ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
    }

    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

    const char *option;
    option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
    BARRIER_RANK_START(pcomm)
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);
    BARRIER_RANK_END(pcomm)

    //Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    //set entitities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle meshset_level0;
    rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRG(rval);
    ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

    //Fields
    ierr = m_field.add_field("TEMP",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
    ierr = m_field.add_field("TEMP_RATE",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);

    //Problem
    ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

    //set refinement level for problem
    ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);

    //meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    //add entities to field
    ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"TEMP"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"TEMP_RATE"); CHKERRQ(ierr);

    //set app. order
    //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
    int order = 2;
    ierr = m_field.set_field_order(root_set,MBTET,"TEMP",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBTRI,"TEMP",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBEDGE,"TEMP",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBVERTEX,"TEMP",1); CHKERRQ(ierr);

    ierr = m_field.set_field_order(root_set,MBTET,"TEMP_RATE",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBTRI,"TEMP_RATE",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBEDGE,"TEMP_RATE",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBVERTEX,"TEMP_RATE",1); CHKERRQ(ierr);

    ThermalElement thermal_elements(m_field);
    ierr = thermal_elements.addThermalElements("TEMP"); CHKERRQ(ierr);
    ierr = thermal_elements.addThermalFluxElement("TEMP"); CHKERRQ(ierr);
    //add rate of temerature to data field of finite element
    ierr = m_field.modify_finite_element_add_field_data("THERMAL_FE","TEMP_RATE"); CHKERRQ(ierr);

    ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","THERMAL_FE"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","THERMAL_FLUX_FE"); CHKERRQ(ierr);

    /****/
    //build database
    //build field
    ierr = m_field.build_fields(); CHKERRQ(ierr);
    //build finite elemnts
    ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
    //build adjacencies
    ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
    //build problem
    ProblemsManager *prb_mng_ptr;
    ierr = m_field.getInterface(prb_mng_ptr); CHKERRQ(ierr);
    ierr = prb_mng_ptr->buildProblem("TEST_PROBLEM",true); CHKERRQ(ierr);

    /****/
    //mesh partitioning
    //partition
    ierr = prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM"); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM"); CHKERRQ(ierr);
    //what are ghost nodes, see Petsc Manual
    ierr = prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM"); CHKERRQ(ierr);

    Vec F;
    ierr = m_field.getInterface<VecManager>()->vecCreateGhost("TEST_PROBLEM",ROW,&F); CHKERRQ(ierr);
    Vec T;
    ierr = VecDuplicate(F,&T); CHKERRQ(ierr);
    Mat A;
    ierr = m_field.MatCreateMPIAIJWithArrays("TEST_PROBLEM",&A); CHKERRQ(ierr);

    //TS
    TsCtx ts_ctx(m_field,"TEST_PROBLEM");
    TS ts;
    ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
    ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);

    DirichletTemperatureBc my_dirichlet_bc(m_field,"TEMP",A,T,F);
    ThermalElement::UpdateAndControl update_velocities(m_field,"TEMP","TEMP_RATE");
    ThermalElement::TimeSeriesMonitor monitor(m_field,"THEMP_SERIES","TEMP");

    //preprocess
    ts_ctx.get_preProcess_to_do_IFunction().push_back(&update_velocities);
    ts_ctx.get_preProcess_to_do_IFunction().push_back(&my_dirichlet_bc);
    ts_ctx.get_preProcess_to_do_IJacobian().push_back(&my_dirichlet_bc);

    //and temperature element functions
    ierr = thermal_elements.setTimeSteppingProblem(ts_ctx,"TEMP","TEMP_RATE"); CHKERRQ(ierr);

    //postprocess
    ts_ctx.get_postProcess_to_do_IFunction().push_back(&my_dirichlet_bc);
    ts_ctx.get_postProcess_to_do_IJacobian().push_back(&my_dirichlet_bc);
    ts_ctx.get_postProcess_to_do_Monitor().push_back(&monitor);

    ierr = TSSetIFunction(ts,F,f_TSSetIFunction,&ts_ctx); CHKERRQ(ierr);
    ierr = TSSetIJacobian(ts,A,A,f_TSSetIJacobian,&ts_ctx); CHKERRQ(ierr);
    ierr = TSMonitorSet(ts,f_TSMonitorSet,&ts_ctx,PETSC_NULL); CHKERRQ(ierr);

    double ftime = 1;
    ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
    ierr = TSSetSolution(ts,T); CHKERRQ(ierr);
    ierr = TSSetFromOptions(ts); CHKERRQ(ierr);

    SeriesRecorder *recorder_ptr;
    ierr = m_field.getInterface(recorder_ptr); CHKERRQ(ierr);
    ierr = recorder_ptr->add_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);
    ierr = recorder_ptr->initialize_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);

    #if PETSC_VERSION_GE(3,7,0)
    ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER); CHKERRQ(ierr);
    #endif
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

      // PetscViewer viewer;
      // PetscViewerASCIIOpen(PETSC_COMM_WORLD,"thermal_elem_unsteady.txt",&viewer);

      double sum = 0;
      double fnorm = 0;

      for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder_ptr,"THEMP_SERIES",sit)) {

        ierr = recorder_ptr->load_series_data("THEMP_SERIES",sit->get_step_number()); CHKERRQ(ierr);
        ierr = m_field.getInterface<VecManager>()->setLocalGhostVector("TEST_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

        double sum0;
        ierr = VecSum(T,&sum0); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"sum0  = %9.8e\n",sum0); CHKERRQ(ierr);
        double fnorm0;
        ierr = VecNorm(T,NORM_2,&fnorm0); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"fnorm0  = %9.8e\n",fnorm0); CHKERRQ(ierr);

        sum += sum0;
        fnorm += fnorm0;

        // ierr = VecChop(T,1e-4); CHKERRQ(ierr);
        // ierr = VecView(T,viewer); CHKERRQ(ierr);

      }
      ierr = PetscPrintf(PETSC_COMM_WORLD,"sum  = %9.8e\n",sum); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"fnorm  = %9.8e\n",fnorm); CHKERRQ(ierr);
      if(fabs(sum+1.32314077e+01)>1e-7) {
        SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Failed to pass test");
      }
      if(fabs(fnorm-4.59664623e+00)>1e-6) {
        SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Failed to pass test");
      }


      // ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);


      /*PostProcVertexMethod ent_method(moab,"TEMP");
      ierr = m_field.loop_dofs("TEST_PROBLEM","TEMP",ROW,ent_method); CHKERRQ(ierr);
      if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERRG(rval);
      ierr = m_field.get_problem_finite_elements_entities("TEST_PROBLEM","THERMAL_FE",out_meshset); CHKERRQ(ierr);
      rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERRG(rval);
      rval = moab.delete_entities(&out_meshset,1); CHKERRG(rval);
    }*/

    ierr = TSDestroy(&ts);CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = VecDestroy(&F); CHKERRQ(ierr);
    ierr = VecDestroy(&T); CHKERRQ(ierr);


  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
