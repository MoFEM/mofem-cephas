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


#include <Projection10NodeCoordsOnField.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <calculate_wt.hpp>

using namespace boost::numeric;
using namespace ObosleteUsersModules;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

const double young_modulus = 1;
const double poisson_ratio = 0.0;

int main(int argc, char *argv[]) {
  
  PetscInitialize(&argc,&argv,(char *)0,help);
  
  moab::Core mb_instance;
  Interface& moab = mb_instance;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
  
  //Reade parameters from line command
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 5;
  }
  
  PetscBool is_partitioned = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-my_is_partitioned",&is_partitioned,&flg); CHKERRQ(ierr);
  
  if(is_partitioned == PETSC_TRUE) {
    //Read mesh to MOAB
    const char *option;
    option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
    rval = pcomm->resolve_shared_ents(0,3,0); CHKERR_PETSC(rval);
    rval = pcomm->resolve_shared_ents(0,3,1); CHKERR_PETSC(rval);
    rval = pcomm->resolve_shared_ents(0,3,2); CHKERR_PETSC(rval);
  } else {
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  }
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;
  
  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(0,bit_level0,1,PETSC_COMM_WORLD); CHKERRQ(ierr);
  Range meshset_level0;
  ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"meshset_level0 %d\n",meshset_level0.size());
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  
  //Define problem
  
  //Fields
  ierr = m_field.add_field("DISPLACEMENT",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  
  //FE
  ierr = m_field.add_finite_element("ELASTIC",MF_ZERO); CHKERRQ(ierr);
  
  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  //define problems
  ierr = m_field.add_problem("ELASTIC_PROB"); CHKERRQ(ierr);
  
  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_PROB","ELASTIC"); CHKERRQ(ierr);
  
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("ELASTIC_PROB",bit_level0); CHKERRQ(ierr);
  
  //Declare problem
  
  //add entitities (by tets) to the field
  ierr = m_field.add_ents_to_field_by_TETs(0,"DISPLACEMENT",PETSC_COMM_WORLD,2); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS",PETSC_COMM_WORLD,2); CHKERRQ(ierr);
  
  //add finite elements entities
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  
  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  ierr = m_field.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  
  //define elements for calculaiton of degradaiton parameter w_t
  Calculate_wt calculate_wt_elem(m_field);
  if(m_field.check_field("TEMP")) {
    cout<<"Temprature field exists "<< endl;
    if(m_field.check_field("CONC")) {
      cout<<"Concentration field exists "<< endl;
    }
    ierr = calculate_wt_elem.addCalculateWtElement("ELASTIC_PROB","ELASTIC","DISPLACEMENT","TEMP","CONC"); CHKERRQ(ierr);
  }
  
  //build database
  
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  
  //build problem
  //ierr = m_field.build_problems(); CHKERRQ(ierr);
  if(is_partitioned) {
    ierr = m_field.build_partitioned_problems(PETSC_COMM_WORLD,1); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("ELASTIC_PROB",true,0,pcomm->size(),1); CHKERRQ(ierr);
  } else {
    ierr = m_field.build_problems(); CHKERRQ(ierr);
    ierr = m_field.partition_problem("ELASTIC_PROB"); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("ELASTIC_PROB"); CHKERRQ(ierr);
  }
  
  ierr = m_field.partition_ghost_dofs("ELASTIC_PROB"); CHKERRQ(ierr);
  
  //m_field.list_dofs_by_field_name("DISPLACEMENT",true);
  
  
  //create matrices
  Vec F,D;
  ierr = m_field.VecCreateGhost("ELASTIC_PROB",ROW,&F); CHKERRQ(ierr);
  ierr = m_field.VecCreateGhost("ELASTIC_PROB",COL,&D); CHKERRQ(ierr);
  
  Mat Aij;
  ierr = m_field.MatCreateMPIAIJWithArrays("ELASTIC_PROB",&Aij); CHKERRQ(ierr);
  
  //Matrix View
  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;
  

  Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
  
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  
  SeriesRecorder& recorder = core;
  if(m_field.check_field("TEMP")) {
    if(m_field.check_field("CONC")) {
      cout<<"Both TEMP and CONC fields exists on the mesh "<<endl;
      //read time series and do thermo elastci analysis
      Vec F_thermal;
      ierr = VecDuplicate(F,&F_thermal); CHKERRQ(ierr);
      
      ierr = calculate_wt_elem.setCalculateWtOperators("DISPLACEMENT","TEMP","CONC",F_thermal); CHKERRQ(ierr);
      
      if( recorder.check_series("THEMP_SERIES") ) {
        if( recorder.check_series("CONC_SERIES") ) {
          cout<<"Both THEMP_SERIES and CONC_SERIES fields exists "<<endl;
          
          typedef SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator SERIES_iterator;
          SERIES_iterator sr_temp_it ,hi_sr_temp_it, sr_conc_it ,hi_sr_conc_it;
          
          sr_temp_it    =recorder.get_series_steps_byName_begin("THEMP_SERIES");
          hi_sr_temp_it =recorder.get_series_steps_byName_end("THEMP_SERIES");

          sr_conc_it    =recorder.get_series_steps_byName_begin("CONC_SERIES");
          hi_sr_conc_it =recorder.get_series_steps_byName_end("CONC_SERIES");
          
          
          //Here steps in both series are the same
          for(;sr_temp_it!=hi_sr_temp_it;  sr_temp_it++, sr_conc_it++){
            PetscPrintf(PETSC_COMM_WORLD,"Process temp series step %d\n",sr_temp_it->get_step_number());
            PetscPrintf(PETSC_COMM_WORLD,"Process conc series step %d\n",sr_conc_it->get_step_number());
            
            ierr = recorder.load_series_data("CONC_SERIES" ,sr_conc_it->get_step_number()); CHKERRQ(ierr);
            ierr = recorder.load_series_data("THEMP_SERIES",sr_temp_it->get_step_number()); CHKERRQ(ierr);

            ierr = m_field.loop_finite_elements("ELASTIC_PROB","ELASTIC",calculate_wt_elem.getLoopDegradation()); CHKERRQ(ierr);

          }
//          PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());
        }
      }
    }
  }
  
      
      
      
//
//      for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder,"THEMP_SERIES",sit))
//        
//        PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());
//        ierr = recorder.load_series_data("THEMP_SERIES",sit->get_step_number()); CHKERRQ(ierr);
//        ierr = VecZeroEntries(F_thermal); CHKERRQ(ierr);
//        ierr = VecGhostUpdateBegin(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//        ierr = VecGhostUpdateEnd(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//        
//        ierr = m_field.loop_finite_elements("ELASTIC_PROB","ELASTIC",thermal_stress_elem.getLoopThermalStressRhs()); CHKERRQ(ierr);
//        ierr = VecGhostUpdateBegin(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//        ierr = VecGhostUpdateEnd(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//        ierr = VecAssemblyBegin(F_thermal); CHKERRQ(ierr);
//        ierr = VecAssemblyEnd(F_thermal); CHKERRQ(ierr);
//        
//        PetscReal nrm_F;
//        ierr = VecNorm(F,NORM_2,&nrm_F); CHKERRQ(ierr);
//        PetscPrintf(PETSC_COMM_WORLD,"norm2 F = %6.4e\n",nrm_F);
//        
//        PetscReal nrm_F_thremal;
//        ierr = VecNorm(F_thermal,NORM_2,&nrm_F_thremal); CHKERRQ(ierr);
//        PetscPrintf(PETSC_COMM_WORLD,"norm2 F_thernal = %6.4e\n",nrm_F_thremal);
//        
//        ierr = VecScale(F_thermal,-1); CHKERRQ(ierr);
//        ierr = VecAXPY(F_thermal,1,F); CHKERRQ(ierr);
//        
//        my_dirichlet_bc.snes_x = D;
//        my_dirichlet_bc.snes_f = F_thermal;
//        ierr = m_field.problem_basic_method_postProcess("ELASTIC_PROB",my_dirichlet_bc); CHKERRQ(ierr);
//        
//        ierr = KSPSolve(solver,F_thermal,D); CHKERRQ(ierr);
//        ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//        ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//        
//        //Save data on mesh
//        ierr = m_field.set_local_VecCreateGhost("ELASTIC_PROB",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//        ierr = m_field.loop_finite_elements("ELASTIC_PROB","ELASTIC",post_proc); CHKERRQ(ierr);
//        ostringstream o1;
//        o1 << "out_" << sit->step_number << ".h5m";
//        rval = post_proc.postProcMesh.write_file(o1.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
//    
//      }
//    }
//    }
  


  
  
//    else {
//
//      ierr = VecZeroEntries(F_thermal); CHKERRQ(ierr);
//      ierr = VecGhostUpdateBegin(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//      ierr = VecGhostUpdateEnd(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//      
//      ierr = m_field.loop_finite_elements("ELASTIC_PROB","ELASTIC",thermal_stress_elem.getLoopThermalStressRhs()); CHKERRQ(ierr);
//      ierr = VecGhostUpdateBegin(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//      ierr = VecGhostUpdateEnd(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//      ierr = VecAssemblyBegin(F_thermal); CHKERRQ(ierr);
//      ierr = VecAssemblyEnd(F_thermal); CHKERRQ(ierr);
//      
//      PetscReal nrm_F;
//      ierr = VecNorm(F,NORM_2,&nrm_F); CHKERRQ(ierr);
//      PetscPrintf(PETSC_COMM_WORLD,"norm2 F = %6.4e\n",nrm_F);
//      
//      PetscReal nrm_F_thremal;
//      ierr = VecNorm(F_thermal,NORM_2,&nrm_F_thremal); CHKERRQ(ierr);
//      PetscPrintf(PETSC_COMM_WORLD,"norm2 F_thernal = %6.4e\n",nrm_F_thremal);
//      
//      ierr = VecScale(F_thermal,-1); CHKERRQ(ierr);
//      ierr = VecAXPY(F_thermal,1,F); CHKERRQ(ierr);
//      
//      my_dirichlet_bc.snes_x = D;
//      my_dirichlet_bc.snes_f = F_thermal;
//      ierr = m_field.problem_basic_method_postProcess("ELASTIC_PROB",my_dirichlet_bc); CHKERRQ(ierr);
//      
//      ierr = KSPSolve(solver,F_thermal,D); CHKERRQ(ierr);
//      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//      
//      //Save data on mesh
//      ierr = m_field.set_local_VecCreateGhost("ELASTIC_PROB",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//      ierr = m_field.loop_finite_elements("ELASTIC_PROB","ELASTIC",post_proc); CHKERRQ(ierr);
//      rval = post_proc.postProcMesh.write_file("out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
//      
//    }
//    
//    ierr = VecDestroy(&F_thermal); CHKERRQ(ierr);
//    
//  } else {
//    
//    // elastic analys
//    
//    ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
//    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//    
//    //Save data on mesh
//    ierr = m_field.set_local_VecCreateGhost("ELASTIC_PROB",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//    ierr = m_field.loop_finite_elements("ELASTIC_PROB","ELASTIC",post_proc); CHKERRQ(ierr);
//    rval = post_proc.postProcMesh.write_file("out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
//    
//  }
//  
//  //Destroy matrices
//  ierr = VecDestroy(&F); CHKERRQ(ierr);
//  ierr = VecDestroy(&D); CHKERRQ(ierr);
//  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
//  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  
  PetscFinalize();
  
}

