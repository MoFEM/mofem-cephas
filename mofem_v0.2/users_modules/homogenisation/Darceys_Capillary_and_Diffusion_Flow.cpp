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

#include <Diffusion_and_Capillary_Element.hpp>

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
  
  /***/
  //Define problem
  //Field
  int field_rank=1;
  ierr = mField.add_field("CONC",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("LAGRANGE_MUL_FIELD_CONC",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("PRESSURE",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("LAGRANGE_MUL_FIELD_PRESSURE",H1,field_rank); CHKERRQ(ierr);

  //Problem
  ierr = mField.add_problem("DIFFUSION_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.add_problem("CAPILLARY_PROBLEM"); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();

  //FE
  MoistureElement moisture_elements(mField);
  ierr = moisture_elements.addDiffusionElements("DIFFUSION_PROBLEM","CONC"); CHKERRQ(ierr);
  ierr = moisture_elements.addCapillaryElements("CAPILLARY_PROBLEM","PRESSURE"); CHKERRQ(ierr);
  
  ierr = mField.add_finite_element("LAGRANGE_FE_CONC"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("LAGRANGE_FE_PRESSURE"); CHKERRQ(ierr);
  
  //Diffusion transport (C and C^T)
  //=====================================================================================================
  //C
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE_CONC","LAGRANGE_MUL_FIELD_CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE_CONC","CONC"); CHKERRQ(ierr);
  //C^T
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE_CONC","CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE_CONC","LAGRANGE_MUL_FIELD_CONC"); CHKERRQ(ierr);
  //Data
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE_CONC","LAGRANGE_MUL_FIELD_CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE_CONC","CONC"); CHKERRQ(ierr);
  //=====================================================================================================

  //Cappillary transport (C and C^T)
  //=====================================================================================================
  //C
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE_PRESSURE","LAGRANGE_MUL_FIELD_PRESSURE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE_PRESSURE","PRESSURE"); CHKERRQ(ierr);
  //C^T
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE_PRESSURE","PRESSURE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE_PRESSURE","LAGRANGE_MUL_FIELD_PRESSURE"); CHKERRQ(ierr);
  //Data
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE_PRESSURE","LAGRANGE_MUL_FIELD_PRESSURE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE_PRESSURE","PRESSURE"); CHKERRQ(ierr);
  //=====================================================================================================
  
  
  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("DIFFUSION_PROBLEM","LAGRANGE_FE_CONC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("CAPILLARY_PROBLEM","LAGRANGE_FE_PRESSURE"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("DIFFUSION_PROBLEM",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("CAPILLARY_PROBLEM",bit_level0); CHKERRQ(ierr);

  /***/
  //Declare problem
  //adding
  Range TetsMatix, TetsFibre;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,BLOCKSET,3,TetsMatix,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,BLOCKSET,3,TetsFibre,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of tets in MAT_MOISTURE_MATRIX = %d\n",TetsMatix.size()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of tets in MAT_MOISTURE_FIBRE = %d\n",TetsFibre.size()); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(TetsMatix,"CONC",2); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(TetsFibre,"PRESSURE",2); CHKERRQ(ierr);

  Range BoundFibres, BoundMatrix;
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,SIDESET,2,BoundFibres,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 101 = %d\n",BoundFibres.size()); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(102,SIDESET,2,BoundMatrix,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 102 = %d\n",BoundMatrix.size()); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(BoundMatrix,"LAGRANGE_FE_CONC"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(BoundFibres,"LAGRANGE_FE_PRESSURE"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TRIs(BoundMatrix,"LAGRANGE_MUL_FIELD_CONC",2); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TRIs(BoundFibres,"LAGRANGE_MUL_FIELD_PRESSURE",2); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
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

  ierr = mField.set_field_order(0,MBTRI,"LAGRANGE_MUL_FIELD_CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"LAGRANGE_MUL_FIELD_CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"LAGRANGE_MUL_FIELD_CONC",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTRI,"LAGRANGE_MUL_FIELD_PRESSURE",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"LAGRANGE_MUL_FIELD_PRESSURE",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"LAGRANGE_MUL_FIELD_PRESSURE",1); CHKERRQ(ierr);

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
  
  
  /****/
  //mesh partitioning
  
  //partition
  ierr = mField.partition_problem("DIFFUSION_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_problem("CAPILLARY_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("DIFFUSION_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("CAPILLARY_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("DIFFUSION_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("CAPILLARY_PROBLEM"); CHKERRQ(ierr);

  //print block sets with materials
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);
  
  //create matrices (here F, D and A are matrices for the full problem)
  Vec F_Dif,D_Dif,F_Cap,D_Cap;   //Here Dif=Diffusivity & Cap=Capillary
  ierr = mField.VecCreateGhost("DIFFUSION_PROBLEM",ROW,&F_Dif); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("CAPILLARY_PROBLEM",ROW,&F_Cap); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("DIFFUSION_PROBLEM",COL,&D_Dif); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("CAPILLARY_PROBLEM",COL,&D_Cap); CHKERRQ(ierr);

  Mat A_Dif,A_Cap;
  ierr = mField.MatCreateMPIAIJWithArrays("DIFFUSION_PROBLEM",&A_Dif); CHKERRQ(ierr);
  ierr = mField.MatCreateMPIAIJWithArrays("CAPILLARY_PROBLEM",&A_Cap); CHKERRQ(ierr);
  
  ierr = moisture_elements.setDiffusionElementLhsOperators("CONC",A_Dif); CHKERRQ(ierr);
  ierr = moisture_elements.setCapillaryElementLhsOperators("PRESSURE",A_Cap); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Disp MyFE_RVELagrange_conc(mField,A_Dif,D_Dif,F_Dif,applied_ConGrad,"CONC","LAGRANGE_MUL_FIELD_CONC",field_rank);
  ElasticFE_RVELagrange_Disp MyFE_RVELagrange_pressure(mField,A_Cap,D_Cap,F_Cap,applied_PressGrad,"PRESSURE","LAGRANGE_MUL_FIELD_PRESSURE",field_rank);

  ierr = VecZeroEntries(F_Dif); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_Dif,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Dif,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(A_Dif); CHKERRQ(ierr);

  ierr = VecZeroEntries(F_Cap); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_Cap,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Cap,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(A_Cap); CHKERRQ(ierr);

  ierr = mField.loop_finite_elements("DIFFUSION_PROBLEM","DIFFUSION_FE",moisture_elements.getLoopFeLhsDiffusion()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("DIFFUSION_PROBLEM","LAGRANGE_FE_CONC",MyFE_RVELagrange_conc);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("CAPILLARY_PROBLEM","CAPILLARY_FE",moisture_elements.getLoopFeLhsCapillary()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("CAPILLARY_PROBLEM","LAGRANGE_FE_PRESSURE",MyFE_RVELagrange_pressure);  CHKERRQ(ierr);

//
  ierr = VecGhostUpdateBegin(F_Dif,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Dif,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F_Dif); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F_Dif); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A_Dif,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_Dif,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  
  ierr = VecGhostUpdateBegin(F_Cap,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Cap,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F_Cap); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F_Cap); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A_Cap,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_Cap,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

//  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  ierr = MatView(A_Dif,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  cout<<"\n\n\n\n\n\n"<<endl;
//  ierr = MatView(A_Cap,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

//  //Matrix View
//  MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//  std::string wait;
//  std::cin >> wait;

  //Solver
  KSP solver_Dif;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver_Dif); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver_Dif,A_Dif,A_Dif); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver_Dif); CHKERRQ(ierr);
  ierr = KSPSetUp(solver_Dif); CHKERRQ(ierr);
  
  ierr = KSPSolve(solver_Dif,F_Dif,D_Dif); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D_Dif,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D_Dif,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("DIFFUSION_PROBLEM",ROW,D_Dif,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  //Solver
  KSP solver_Cap;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver_Cap); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver_Cap,A_Cap,A_Cap); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver_Cap); CHKERRQ(ierr);
  ierr = KSPSetUp(solver_Cap); CHKERRQ(ierr);
  
  ierr = KSPSolve(solver_Cap,F_Cap,D_Cap); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D_Cap,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D_Cap,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("CAPILLARY_PROBLEM",ROW,D_Cap,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);



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

  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet_Dif(mField,"CONC",true,false,"CONC");
  ierr = mField.loop_dofs("CONC",ent_method_on_10nodeTet_Dif); CHKERRQ(ierr);
  ent_method_on_10nodeTet_Dif.set_nodes = false;
  ierr = mField.loop_dofs("CONC",ent_method_on_10nodeTet_Dif); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset_Dif;
    rval = moab.create_meshset(MESHSET_SET,out_meshset_Dif); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("DIFFUSION_PROBLEM","DIFFUSION_FE",out_meshset_Dif); CHKERRQ(ierr);
    rval = moab.write_file("out_Dif.vtk","VTK","",&out_meshset_Dif,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset_Dif,1); CHKERR_PETSC(rval);
  }
  
  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet_Cap(mField,"PRESSURE",true,false,"PRESSURE");
  ierr = mField.loop_dofs("PRESSURE",ent_method_on_10nodeTet_Cap); CHKERRQ(ierr);
  ent_method_on_10nodeTet_Cap.set_nodes = false;
  ierr = mField.loop_dofs("PRESSURE",ent_method_on_10nodeTet_Cap); CHKERRQ(ierr);
  
  if(pcomm->rank()==0) {
    EntityHandle out_meshset_Cap;
    rval = moab.create_meshset(MESHSET_SET,out_meshset_Cap); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("CAPILLARY_PROBLEM","CAPILLARY_FE",out_meshset_Cap); CHKERRQ(ierr);
    rval = moab.write_file("out_Cap.vtk","VTK","",&out_meshset_Cap,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset_Cap,1); CHKERR_PETSC(rval);
  }

  
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

