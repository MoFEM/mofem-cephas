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
#include "FEMethod_UpLevelStudent.hpp"
#include "PotentialFlowFEMethod.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "ThermalFEMethod.hpp"
#include "PostProcVertexMethodTemp.hpp"
#include "PostProcTemperatureOnRefindedMesh.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;



static char help[] = "...\n\n";

int main(int argc, char *argv[]) {
    

  PetscInitialize(&argc,&argv,(char *)0,help);
    
    //We need that for code profiling
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;
    ierr = PetscTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
  
  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  /*if(rank==0) {
    EntityHandle dummy_meshset;
    rval = moab.create_meshset(MESHSET_SET,dummy_meshset); CHKERR_PETSC(rval);
  }*/

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }


  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  FieldCore core(moab);
  FieldInterface& mField = core;

  //add filds
  ierr = mField.add_field("TEMPERATURE",H1,1); CHKERRQ(ierr);

  //add finite elements
  ierr = mField.add_finite_element("THERMAL"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("THERMAL","TEMPERATURE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("THERMAL","TEMPERATURE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("THERMAL","TEMPERATURE"); CHKERRQ(ierr);

  //add problems 
  ierr = mField.add_problem("THERMAL_PROBLEM"); CHKERRQ(ierr);

  //define problems and finite elements
  ierr = mField.modify_problem_add_finite_element("THERMAL_PROBLEM","THERMAL"); CHKERRQ(ierr);

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("THERMAL_PROBLEM",bit_level0); CHKERRQ(ierr);
      
  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"TEMPERATURE"); CHKERRQ(ierr);
  
  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"THERMAL",MBTET); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  //set app. order
  ierr = mField.set_field_order(0,MBVERTEX,"TEMPERATURE",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"TEMPERATURE",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"TEMPERATURE",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"TEMPERATURE",order); CHKERRQ(ierr);

  //build fields
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elements
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = mField.partition_problem("THERMAL_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("THERMAL_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("THERMAL_PROBLEM"); CHKERRQ(ierr);

  //print bcs
  ierr = mField.printCubitTemperatureSet(); CHKERRQ(ierr);
  ierr = mField.printCubitHeatFluxSet(); CHKERRQ(ierr);

  //print block sets with materials
  ierr = mField.printCubitMaterials(); CHKERRQ(ierr);
    
  //create matrices and vectors
  Vec F,D;
  ierr = mField.VecCreateGhost("THERMAL_PROBLEM",Row,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("THERMAL_PROBLEM",Col,&D); CHKERRQ(ierr);

  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("THERMAL_PROBLEM",&Aij); CHKERRQ(ierr);

  struct MyThermalFEMethod: public ThermalFEMethod {
    MyThermalFEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,
                        Mat &_Aij,Vec &_D,Vec& _F,double _Ther_Cond):   //class constructor
	ThermalFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_Ther_Cond) {};
    //This is for KSP solver, residual is on the RHS
    PetscErrorCode Rhs(Vec F) {
      PetscFunctionBegin;
      ierr = Fint(); CHKERRQ(ierr);
      for(int rr = 0;rr<row_mat;rr++) {
	if(RowGlob[rr].size()==0) continue;
	if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	f_int[rr] *= -1;
	ierr = VecSetValues(F,RowGlob[rr].size(),&(RowGlob[rr])[0],&*f_int[rr].data().begin(),ADD_VALUES); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

   };
    
    
   CubitTemperatureDirihletBC myDirihletBC(mField,"THERMAL_PROBLEM","TEMPERATURE");
   ierr = myDirihletBC.Init(); CHKERRQ(ierr);
      
  //Assemble F and Aij
  double Ther_Cond=100;
    
  MyThermalFEMethod MyFE(mField,&myDirihletBC,Aij,D,F,Ther_Cond);
    
//    cout<<"F "<<endl<<endl<<endl<<endl;
//    PetscViewer    viewer;
//    VecView(F,viewer);
//    cout<<endl<<endl<<endl<<endl;

ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  
  ierr = mField.loop_finite_elements("THERMAL_PROBLEM","THERMAL",MyFE);  CHKERRQ(ierr);
    //cout<< "Hi "<<endl;
    
//    cout<<"Rhs "<<endl<<endl<<endl<<endl;
//    PetscViewer    viewer;
//    VecView(F,viewer);
//    cout<<endl<<endl<<endl<<endl;
    

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  //PetscSynchronizedFlush(PETSC_COMM_WORLD);

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);
    
  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("THERMAL_PROBLEM",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PostProcVertexMethodTemp ent_method(moab);
  ierr = mField.loop_dofs("THERMAL_PROBLEM","TEMPERATURE",Row,ent_method); CHKERRQ(ierr);
    
  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("THERMAL_PROBLEM","THERMAL",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);

  }

  PostProcTemperatureOnRefMesh post_proc_on_ref_mesh(moab);
  ierr = mField.loop_finite_elements("THERMAL_PROBLEM","THERMAL",post_proc_on_ref_mesh); CHKERRQ(ierr);
  if(pcomm->rank()==0) {
    rval = post_proc_on_ref_mesh.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

    //Display the temperature vector and Aij matrices
    //cout<<"D "<<endl<<endl;
    PetscViewer    viewer;
    VecView(D,viewer);
    cout<<endl<<endl;
    //cout<<"Aij "<<endl<<endl;
    //MatView(Aij, viewer);
    //cout<<endl<<endl<<endl<<endl;
 
  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
    
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  
  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  
  PetscFinalize();

}
