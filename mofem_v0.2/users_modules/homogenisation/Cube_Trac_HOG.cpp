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

#include <Projection10NodeCoordsOnField.hpp>
#include <petsctime.h>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcOnRefMesh.hpp>
#include <PostProcHookStresses.hpp>

#include <adolc/adolc.h>
#include <NonLienarElasticElement.hpp>
#include <Hooke.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

#include "BCs_RVELagrange_Disp.hpp"
#include "BCs_RVELagrange_Trac.hpp"
#include "BCs_RVELagrange_Trac_Rigid_Trans.hpp"
#include "BCs_RVELagrange_Trac_Rigid_Rot.hpp"
#include "BCs_RVEVolume.hpp"
#include "BCs_RVE_Homogenized_Stress_Trac.hpp"


//#include "ElasticFE_RVELagrange_Homogenized_Stress_Traction.hpp"
//#include "ElasticFE_RVELagrange_RigidBodyTranslation.hpp"
//#include "ElasticFE_RVELagrange_RigidBodyRotation.hpp"
//#include "RVEVolume.hpp"

ErrorCode rval;
PetscErrorCode ierr;


PetscErrorCode Stress_cal(FieldInterface &mField, string field_name,string lagrang_field_name, double RVE_volume, map<int,BCs_RVELagrange_Disp::RVEBC_Data> &setOfRVEBC, Vec Stress_Homo, ublas::matrix<FieldData> &Dmat, int Dmatcol) {
  PetscFunctionBegin;
  
  ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
  BCs_RVE_Homogenized_Stress_Trac rve_homo_stress_trac(mField);
  rve_homo_stress_trac.setRVEBCsHomoStressOperators("DISPLACEMENT","Lagrange_mul_disp", Stress_Homo, setOfRVEBC,"MESH_NODE_POSITIONS");
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",rve_homo_stress_trac.getLoopFeRVEBCRhs()); CHKERRQ(ierr);
  
  ierr = VecAssemblyBegin(Stress_Homo); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Stress_Homo); CHKERRQ(ierr);
  ierr = VecScale(Stress_Homo, 1.0/(RVE_volume)); CHKERRQ(ierr);
//  ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PetscScalar *avec;
  VecGetArray(Stress_Homo, &avec);
  for(int ii=0; ii<6; ii++){
    Dmat(ii,Dmatcol)=*avec;
    avec++;
  }
  
  PetscFunctionReturn(0);
}





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
    order = 5;
  }
  
  //Applied strain on the RVE (vector of length 6) strain=[xx, yy, zz, xy, xz, zy]^T
  double myapplied_strain[6];
  int nmax=6;
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-myapplied_strain",myapplied_strain,&nmax,&flg); CHKERRQ(ierr);
  ublas::vector<FieldData> applied_strain;
  applied_strain.resize(6);
  cblas_dcopy(6, &myapplied_strain[0], 1, &applied_strain(0), 1);
  //    cout<<"applied_strain ="<<applied_strain<<endl;
  
  
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
  
  //Fields
  int field_rank=3;
  ierr = mField.add_field("DISPLACEMENT",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("Lagrange_mul_disp",NOFIELD,6); CHKERRQ(ierr);
  ierr = mField.add_field("Lagrange_mul_disp_rigid_trans",NOFIELD,3); CHKERRQ(ierr);   //Control 3 rigid body translations in x, y and z axis
  ierr = mField.add_field("Lagrange_mul_disp_rigid_rotation",NOFIELD,3); CHKERRQ(ierr); //Controla 3 rigid body rotations about x, y and z axis
  
  
  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);

  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  
  //FE
  //Define FE
  //define eleatic element
  Hooke<adouble> hooke_adouble;
  Hooke<double> hooke_double;
  NonlinearElasticElement elastic(mField,2);
  ierr = elastic.setBlocks(&hooke_double,&hooke_adouble); CHKERRQ(ierr);
  ierr = elastic.addElement("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = elastic.setOperators("DISPLACEMENT","MESH_NODE_POSITIONS",false,true); CHKERRQ(ierr);

  
  BCs_RVELagrange_Trac lagrangian_element(mField);
  lagrangian_element.addLagrangiangElement("Lagrange_elem","DISPLACEMENT","Lagrange_mul_disp","MESH_NODE_POSITIONS");
  lagrangian_element.addLagrangiangElement("Lagrange_elem_rigid_trans","DISPLACEMENT","Lagrange_mul_disp_rigid_trans","MESH_NODE_POSITIONS");
  lagrangian_element.addLagrangiangElement("Lagrange_elem_rigid_rotation","DISPLACEMENT","Lagrange_mul_disp_rigid_rotation","MESH_NODE_POSITIONS");
  
  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  
  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","Lagrange_elem"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","Lagrange_elem_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","Lagrange_elem_rigid_rotation"); CHKERRQ(ierr);
  
  
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);
  
  /***/
  //Declare problem
   /****/
  //build database
  
  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  Projection10NodeCoordsOnField ent_method_material(mField,"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  
  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);
  
  /****/
  //mesh partitioning
  
  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  
  //print bcs
  ierr = mField.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = mField.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);
  
  //create matrices (here F, D and Aij are matrices for the full problem)
  Vec F1,F2,F3,F4,F5,F6,D;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F1); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F2); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F3); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F4); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F5); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F6); CHKERRQ(ierr);
  
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",COL,&D); CHKERRQ(ierr);
  
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);
  
  ierr = VecZeroEntries(F1); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  ierr = VecZeroEntries(F2); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  ierr = VecZeroEntries(F3); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  ierr = VecZeroEntries(F4); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  ierr = VecZeroEntries(F5); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  ierr = VecZeroEntries(F6); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  
  
  //Assemble F and Aij
  elastic.getLoopFeLhs().snes_B = Aij;
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",elastic.getLoopFeLhs()); CHKERRQ(ierr);

  lagrangian_element.setRVEBCsOperators("DISPLACEMENT","Lagrange_mul_disp",Aij,F1,F2,F3,F4,F5,F6,"MESH_NODE_POSITIONS");
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",lagrangian_element.getLoopFeRVEBCLhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",lagrangian_element.getLoopFeRVEBCRhs()); CHKERRQ(ierr);
  
  BCs_RVELagrange_Trac_Rigid_Trans lagrangian_element_rigid_body_trans(mField);
  lagrangian_element_rigid_body_trans.setRVEBCsRigidBodyTranOperators("DISPLACEMENT","Lagrange_mul_disp_rigid_trans",Aij, lagrangian_element.setOfRVEBC);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem_rigid_trans",lagrangian_element_rigid_body_trans.getLoopFeRVEBCLhs()); CHKERRQ(ierr);

  BCs_RVELagrange_Trac_Rigid_Rot lagrangian_element_rigid_body_rot(mField);
  lagrangian_element_rigid_body_rot.setRVEBCsRigidBodyRotOperators("DISPLACEMENT","Lagrange_mul_disp_rigid_rotation",Aij, lagrangian_element.setOfRVEBC);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem_rigid_rotation",lagrangian_element_rigid_body_rot.getLoopFeRVEBCLhs()); CHKERRQ(ierr);

//  //Matrix View
//  MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//  std::string wait;
//  std::cin >> wait;

  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(F1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F1); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F1); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(F2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F2); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F2); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(F3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F3); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F3); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(F4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F4); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F4); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(F5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F5); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F5); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(F6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F6); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F6); CHKERRQ(ierr);


  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //ierr = MatView(Aij,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  
////Matrix View
//MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//std::string wait;
//std::cin >> wait;
  
  //  ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //  ierr = MatView(Aij,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  //=============================================================================================================
  //Calculation of RVE volume for homogenised stress calculaiton
  //=============================================================================================================
  double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
  Vec RVE_volume_Vec;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
  ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
  
  BCs_RVEVolume rvevolume_element(mField);
  rvevolume_element.setRVEVolumeOperators(mField, "DISPLACEMENT", RVE_volume_Vec, elastic.setOfBlocks);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",rvevolume_element.getLoopFeLhs()); CHKERRQ(ierr);
  //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
  cout<<"Final RVE_volume = "<< RVE_volume <<endl;
  //=============================================================================================================
  
  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);
  
  //solve for F1
  ierr = KSPSolve(solver,F1,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  //Save data on mesh
  ierr = mField.set_global_ghost_vector("ELASTIC_MECHANICS",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

//  ierr = VecView(F1,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  //=============================================================================================================
  // homogenised stress for strian [1 0 0 0 0 0]^T
  //=============================================================================================================
  ublas::matrix<FieldData> Dmat;
  Dmat.resize(6,6); Dmat.clear();

  //create a vector for 6 components of homogenized stress
  Vec Stress_Homo;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
  
  int Dmatcol=0;
  ierr = Stress_cal(mField, "DISPLACEMENT","Lagrange_mul_disp", RVE_volume, lagrangian_element.setOfRVEBC, Stress_Homo, Dmat, Dmatcol); CHKERRQ(ierr);
  

  //=============================================================================================================
  // homogenised stress for strian [0 1 0 0 0 0]^T
  //=============================================================================================================
  //solve for F2
  ierr = KSPSolve(solver,F2,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_ghost_vector("ELASTIC_MECHANICS",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  Dmatcol=1;
  ierr = Stress_cal(mField, "DISPLACEMENT","Lagrange_mul_disp", RVE_volume, lagrangian_element.setOfRVEBC, Stress_Homo, Dmat, Dmatcol); CHKERRQ(ierr);

  
  //=============================================================================================================
  // homogenised stress for strian [0 0 1 0 0 0]^T
  //=============================================================================================================
  //solve for F3
  ierr = KSPSolve(solver,F3,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_ghost_vector("ELASTIC_MECHANICS",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  Dmatcol=2;
  ierr = Stress_cal(mField, "DISPLACEMENT","Lagrange_mul_disp", RVE_volume, lagrangian_element.setOfRVEBC, Stress_Homo, Dmat, Dmatcol); CHKERRQ(ierr);
  //=============================================================================================================
  // homogenised stress for strian [0 0 0 1 0 0]^T
  //=============================================================================================================
  //solve for F4
  ierr = KSPSolve(solver,F4,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_ghost_vector("ELASTIC_MECHANICS",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  Dmatcol=3;
  ierr = Stress_cal(mField, "DISPLACEMENT","Lagrange_mul_disp", RVE_volume, lagrangian_element.setOfRVEBC, Stress_Homo, Dmat, Dmatcol); CHKERRQ(ierr);
  //=============================================================================================================
  // homogenised stress for strian [0 0 0 0 1 0]^T
  //=============================================================================================================
  //solve for F5
  ierr = KSPSolve(solver,F5,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_ghost_vector("ELASTIC_MECHANICS",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  Dmatcol=4;
  ierr = Stress_cal(mField, "DISPLACEMENT","Lagrange_mul_disp", RVE_volume, lagrangian_element.setOfRVEBC, Stress_Homo, Dmat, Dmatcol); CHKERRQ(ierr);
  //=============================================================================================================
  // homogenised stress for strian [0 0 0 0 0 1]^T
  //=============================================================================================================
  //solve for F6
  ierr = KSPSolve(solver,F6,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_ghost_vector("ELASTIC_MECHANICS",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  Dmatcol=5;
  ierr = Stress_cal(mField, "DISPLACEMENT","Lagrange_mul_disp", RVE_volume, lagrangian_element.setOfRVEBC, Stress_Homo, Dmat, Dmatcol); CHKERRQ(ierr);
  
  
  if(pcomm->rank()==0){
    cout<< "\nStress_Homo = \n\n";
    for(int ii=0; ii<6; ii++){
      for (int jj=0; jj<6; jj++){
        cout <<Dmat(ii,jj)<<"    ";
      }
      cout<<endl;
    }
  }


  //=======================================================================================================================================================

  PostPocOnRefinedMesh post_proc(mField);
  ierr = post_proc.generateRefereneElemenMesh(); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("DISPLACEMENT"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesGradientPostProc("DISPLACEMENT"); CHKERRQ(ierr);
  //add postpocessing for sresses
  post_proc.getOpPtrVector().push_back(
                                          new PostPorcStress(
                                                             mField,
                                                             post_proc.postProcMesh,
                                                             post_proc.mapGaussPts,
                                                             "DISPLACEMENT",
                                                             post_proc.commonData));
  
  
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",post_proc); CHKERRQ(ierr);
  rval = post_proc.postProcMesh.write_file("out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);

  
  //Destroy matrices
  //Destroy matrices and vectors
  ierr = VecDestroy(&F1); CHKERRQ(ierr);
  ierr = VecDestroy(&F2); CHKERRQ(ierr);
  ierr = VecDestroy(&F3); CHKERRQ(ierr);
  ierr = VecDestroy(&F4); CHKERRQ(ierr);
  ierr = VecDestroy(&F5); CHKERRQ(ierr);
  ierr = VecDestroy(&F6); CHKERRQ(ierr);
  ierr = VecDestroy(&Stress_Homo); CHKERRQ(ierr);
  ierr = VecDestroy(&RVE_volume_Vec); CHKERRQ(ierr);
  
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  
  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  
  PetscFinalize();
  
}

