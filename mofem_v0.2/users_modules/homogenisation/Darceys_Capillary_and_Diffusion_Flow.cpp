/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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
#include <petsctime.h>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

#include <MoistureElement.hpp>
#include "DarcysElement.hpp"

using namespace boost::numeric;
using namespace ObosleteUsersModules;

#include "ElasticFE_RVELagrange_Disp.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Disp.hpp"
#include "RVEVolume.hpp"

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {
  
  PetscInitialize(&argc,&argv,(char *)0,help);
  
  moab::Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
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
    order = 1;
  }
  
  //Applied concentration gradient on the RVE (vector of length 3) strain=[xx, yy, zz]^T
  double myapplied_ConGrad[3];
  int nmax=3;
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-myapplied_ConGrad",myapplied_ConGrad,&nmax,&flg); CHKERRQ(ierr);
  ublas::vector<FieldData> applied_ConGrad;
  applied_ConGrad.resize(3);
  cblas_dcopy(3, &myapplied_ConGrad[0], 1, &applied_ConGrad(0), 1);
  cout<<"applied_ConGrad ="<<applied_ConGrad<<endl;
  
  //Applied concentration gradient on the RVE (vector of length 3) strain=[xx, yy, zz]^T
  double myapplied_PressGrad[3];
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-myapplied_PressGrad",myapplied_PressGrad,&nmax,&flg); CHKERRQ(ierr);
  ublas::vector<FieldData> applied_PressGrad;
  applied_PressGrad.resize(3);
  cblas_dcopy(3, &myapplied_PressGrad[0], 1, &applied_PressGrad(0), 1);
  cout<<"applied_PressGrad ="<<applied_PressGrad<<endl;
  
  
  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
  
  //We need that for code profiling
  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& mField = core;
  
  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  
  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  
  //Field
  int field_rank=1;
  ierr = mField.add_field("CONC",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("PRESSURE",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("FIELD_LAGRANGE_MUL_CONC",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("FIELD_LAGRANGE_MUL_PRESSURE",H1,field_rank); CHKERRQ(ierr);
  
  //problem
  ierr = mField.add_problem("MOISTURE_PROBLEM"); CHKERRQ(ierr);

 //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  
  MoistureElement diffusion_elements(mField);
  DarcysElement darceys_elements(mField);
  ierr = darceys_elements.addDarceysElements("MOISTURE_PROBLEM","CONC"); CHKERRQ(ierr);
  ierr = darceys_elements.addDarceysElements("MOISTURE_PROBLEM","PRESSURE"); CHKERRQ(ierr);

  //FE diffusion
  ierr = mField.add_finite_element("FE_LAGRANGE_DIFFUSION"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("FE_LAGRANGE_CAPILLARY"); CHKERRQ(ierr);
  //=======================================================================================
  // rows and columns for FE_LAGRANGE_DIFFUSION
  //=======================================================================================
  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("FE_LAGRANGE_DIFFUSION","FIELD_LAGRANGE_MUL_CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("FE_LAGRANGE_DIFFUSION","CONC"); CHKERRQ(ierr);
  
  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("FE_LAGRANGE_DIFFUSION","CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("FE_LAGRANGE_DIFFUSION","FIELD_LAGRANGE_MUL_CONC"); CHKERRQ(ierr);
  
  //data
  ierr = mField.modify_finite_element_add_field_data("FE_LAGRANGE_DIFFUSION","FIELD_LAGRANGE_MUL_CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("FE_LAGRANGE_DIFFUSION","CONC"); CHKERRQ(ierr);
  
  //=======================================================================================
  // rows and columns for FE_LAGRANGE_CAPILLARY
  //=======================================================================================
  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("FE_LAGRANGE_CAPILLARY","FIELD_LAGRANGE_MUL_PRESSURE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("FE_LAGRANGE_CAPILLARY","PRESSURE"); CHKERRQ(ierr);
  
  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("FE_LAGRANGE_CAPILLARY","PRESSURE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("FE_LAGRANGE_CAPILLARY","FIELD_LAGRANGE_MUL_PRESSURE"); CHKERRQ(ierr);
  
  //data
  ierr = mField.modify_finite_element_add_field_data("FE_LAGRANGE_CAPILLARY","FIELD_LAGRANGE_MUL_PRESSURE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("FE_LAGRANGE_CAPILLARY","PRESSURE"); CHKERRQ(ierr);
  //=======================================================================================

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("MOISTURE_PROBLEM","FE_LAGRANGE_DIFFUSION"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("MOISTURE_PROBLEM","FE_LAGRANGE_CAPILLARY"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("MOISTURE_PROBLEM",bit_level0); CHKERRQ(ierr);

  //add entities to field
  Range TetsMatix, TetsFibre;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,BLOCKSET,3,TetsMatix,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,BLOCKSET,3,TetsFibre,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of tets in MAT_MOISTURE_MATRIX = %d\n",TetsMatix.size()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of tets in MAT_MOISTURE_FIBRE = %d\n",TetsFibre.size()); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(TetsMatix,"CONC",2); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(TetsFibre,"PRESSURE",2); CHKERRQ(ierr);
  
  Range BoundFibres, BoundMarix;
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,SIDESET,2,BoundFibres,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 101 = %d\n",BoundFibres.size()); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(102,SIDESET,2,BoundMarix,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 102 = %d\n",BoundMarix.size()); CHKERRQ(ierr);
  EntityHandle BoundFibresMeshset, BoundMarixMeshset;
  rval = moab.create_meshset(MESHSET_SET,BoundFibresMeshset); CHKERR_PETSC(rval);
  rval = moab.create_meshset(MESHSET_SET,BoundMarixMeshset); CHKERR_PETSC(rval);
	rval = moab.add_entities(BoundFibresMeshset,BoundFibres); CHKERR_PETSC(rval);
	rval = moab.add_entities(BoundMarixMeshset,BoundMarix); CHKERR_PETSC(rval);
  ierr = mField.add_ents_to_field_by_TRIs(BoundFibresMeshset,"FIELD_LAGRANGE_MUL_PRESSURE",2); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TRIs(BoundMarixMeshset,"FIELD_LAGRANGE_MUL_CONC",2); CHKERRQ(ierr);
  
  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  ierr = mField.set_field_order(root_set,MBTET,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"CONC",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(root_set,MBTET,"PRESSURE",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"PRESSURE",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"PRESSURE",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"PRESSURE",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTRI,"FIELD_LAGRANGE_MUL_CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"FIELD_LAGRANGE_MUL_CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"FIELD_LAGRANGE_MUL_CONC",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTRI,"FIELD_LAGRANGE_MUL_PRESSURE",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"FIELD_LAGRANGE_MUL_PRESSURE",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"FIELD_LAGRANGE_MUL_PRESSURE",1); CHKERRQ(ierr);

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
  

//  /****/
//  //mesh partitioning
//  
//  //partition
//  ierr = mField.partition_problem("MOISTURE_PROBLEM"); CHKERRQ(ierr);
//  ierr = mField.partition_finite_elements("MOISTURE_PROBLEM"); CHKERRQ(ierr);
//  //what are ghost nodes, see Petsc Manual
//  ierr = mField.partition_ghost_dofs("MOISTURE_PROBLEM"); CHKERRQ(ierr);
//  
//  //print block sets with materials
//  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);
//  
//  //create matrices (here F, D and A are matrices for the full problem)
//  Vec F,D;
//  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",ROW,&F); CHKERRQ(ierr);
//  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",COL,&D); CHKERRQ(ierr);
//  Mat A;
//  ierr = mField.MatCreateMPIAIJWithArrays("MOISTURE_PROBLEM",&A); CHKERRQ(ierr);
//
//  MoistureElement diffusion_elements(mField);
//  DarcysElement darceys_elements(mField);
//
//  ierr = diffusion_elements.setMoistureFiniteElementLhsOperators("CONC",A); CHKERRQ(ierr);
//  ElasticFE_RVELagrange_Disp MyFE_RVELagrange_conc(mField,A,D,F,applied_ConGrad,"CONC","FIELD_LAGRANGE_MUL_CONC",field_rank);
//
//  ierr = darceys_elements.setDarceysFiniteElementLhsOperators("PRESSURE",A); CHKERRQ(ierr);
//  ElasticFE_RVELagrange_Disp MyFE_RVELagrange_pressure(mField,A,D,F,applied_PressGrad,"PRESSURE","FIELD_LAGRANGE_MUL_PRESSURE",field_rank);
//
//  ierr = VecZeroEntries(F); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = MatZeroEntries(A); CHKERRQ(ierr);
//
////  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","FE_DIFFUSION",diffusion_elements.getLoopFeLhs()); CHKERRQ(ierr);
////  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","FE_LAGRANGE_DIFFUSION",MyFE_RVELagrange_conc);  CHKERRQ(ierr);
////
////  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","FE_CAPILLARY",darceys_elements.getLoopFeLhs()); CHKERRQ(ierr);
////  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","FE_LAGRANGE_CAPILLARY",MyFE_RVELagrange_pressure);  CHKERRQ(ierr);
////
//  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
//  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
//  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//
//  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  
//  //Matrix View
//  MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//  std::string wait;
//  std::cin >> wait;
//
//  //Solver
//  KSP solver;
//  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
//  ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
//  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
//  ierr = KSPSetUp(solver); CHKERRQ(ierr);
//  
//  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  
////  ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
////  ierr = VecView(C,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  
//  
//  //Save data on mesh
//  ierr = mField.set_global_VecCreateGhost("MOISTURE_PROBLEM",ROW,C,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  
//  //Calculation of Homogenized stress
//  //=======================================================================================================================================================
//  const double young_modulus = 1;
//  const double poisson_ratio = 0.0;
//
//  double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
//  Vec RVE_volume_Vec;
//  ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
//  ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
//  
//  RVEVolume MyRVEVol(mField,A,C,F,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio), RVE_volume_Vec);
//  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","MOISTURE_FE",MyRVEVol);  CHKERRQ(ierr);
////  ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
//  cout<<"Final RVE_volume = "<< RVE_volume <<endl;
//  
//  
//  //create a vector for 6 components of homogenized stress
//  Vec Stress_Homo;
//  ierr = VecCreateMPI(PETSC_COMM_WORLD, 3, 3*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
//  ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
//  
//  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(mField,A,C,F,&RVE_volume, applied_strain, Stress_Homo,"CONC","LAGRANGE_MUL_FIELD",field_rank);
//  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","LAGRANGE_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
//  
////  if(pcomm->rank()) cout<< " Stress_Homo =  "<<endl;
////  ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  
//  if(pcomm->rank()==0){
//    PetscScalar    *avec;
//    VecGetArray(Stress_Homo, &avec);
//    
//    cout<< "\nStress_Homo = \n\n";
//    for(int ii=0; ii<3; ii++){
//      cout <<*avec<<endl; ;
//      avec++;
//    }
//  }
//  cout<< "\n\n";
//
//  //=======================================================================================================================================================
//
//  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet(mField,"CONC",true,false,"CONC");
//  ierr = mField.loop_dofs("CONC",ent_method_on_10nodeTet); CHKERRQ(ierr);
//  ent_method_on_10nodeTet.set_nodes = false;
//  ierr = mField.loop_dofs("CONC",ent_method_on_10nodeTet); CHKERRQ(ierr);
//
//  
//  if(pcomm->rank()==0) {
//    EntityHandle out_meshset;
//    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
//    ierr = mField.problem_get_FE("MOISTURE_PROBLEM","MOISTURE_FE",out_meshset); CHKERRQ(ierr);
//    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
//    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
//  }
//  
//  //PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
//  PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method(mField,"CONC",LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
//  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","MOISTURE_FE",fe_post_proc_method);  CHKERRQ(ierr);
//  
//  if(pcomm->rank()==0) {
//    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
//  }
//  
//   
//  //Destroy matrices
//  ierr = VecDestroy(&F); CHKERRQ(ierr);
//  ierr = VecDestroy(&C); CHKERRQ(ierr);
//  ierr = MatDestroy(&A); CHKERRQ(ierr);
//  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
//  
//  
//  ierr = PetscTime(&v2);CHKERRQ(ierr);
//  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
//  
//  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
//  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  
  PetscFinalize();
  
}

