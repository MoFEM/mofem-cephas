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
#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

#include <ThermalElement.hpp>
#include <MoistureTransportElement.hpp>
#include "ElasticFE_RVELagrange_Disp.hpp"
#include "ElasticFE_RVELagrange_Disp_Multi_Rhs.hpp"
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
  
  char output_file_Dmat[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_output_file_Dmat",output_file_Dmat,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (NAME OF DMAT OUTPUT FILE IS NEEDED)");
  }

  
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }
  
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
  ierr = mField.add_field("LAGRANGE_MUL_FIELD",H1,field_rank); CHKERRQ(ierr);
  
  
  //Problem
  ierr = mField.add_problem("MOISTURE_PROBLEM"); CHKERRQ(ierr);
  

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"CONC"); CHKERRQ(ierr);

  //FE
  MoistureTransportElement moisture_elements(mField);
  ierr = moisture_elements.addDiffusionElement("MOISTURE_PROBLEM","CONC"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("LAGRANGE_FE"); CHKERRQ(ierr);
  
  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE","LAGRANGE_MUL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE","CONC"); CHKERRQ(ierr);
  
  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE","LAGRANGE_MUL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE","CONC"); CHKERRQ(ierr);
  
  //As for stress we need both displacement and temprature (Lukasz)
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE","LAGRANGE_MUL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE","CONC"); CHKERRQ(ierr);
  
  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("MOISTURE_PROBLEM","LAGRANGE_FE"); CHKERRQ(ierr);
  
  
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("MOISTURE_PROBLEM",bit_level0); CHKERRQ(ierr);
  
  /***/
  //Declare problem
  
  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"CONC"); CHKERRQ(ierr);
  
  //add finite elements entities
  Range SurfacesFaces;
  ierr = mField.get_Cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"LAGRANGE_FE"); CHKERRQ(ierr);
  
  
  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
	rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TRIs(BoundFacesMeshset,"LAGRANGE_MUL_FIELD",2); CHKERRQ(ierr);
  
  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  ierr = mField.set_field_order(root_set,MBTET,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"CONC",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTRI,"LAGRANGE_MUL_FIELD",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"LAGRANGE_MUL_FIELD",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"LAGRANGE_MUL_FIELD",1); CHKERRQ(ierr);
  
  
  
  //create these elements to calculate separte volumes for each fibres and matrix.
  //======================================================================================
  //FE
  ierr = mField.add_finite_element("FE_MATRXIX"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("FE_FIBRES"); CHKERRQ(ierr);
  
  //Row and columns
  ierr = mField.modify_finite_element_add_field_row("FE_MATRXIX","CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("FE_MATRXIX","CONC"); CHKERRQ(ierr);
  
  ierr = mField.modify_finite_element_add_field_row("FE_FIBRES","CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("FE_FIBRES","CONC"); CHKERRQ(ierr);
  
  ierr = mField.modify_finite_element_add_field_data("FE_MATRXIX","CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("FE_FIBRES","CONC"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("MOISTURE_PROBLEM","FE_MATRXIX"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("MOISTURE_PROBLEM","FE_FIBRES"); CHKERRQ(ierr);
  
  //add finite elements entities
  EntityHandle meshset_Elastic, meshset_Elastic_stiff_inclusion;
  rval = moab.create_meshset(MESHSET_SET,meshset_Elastic); CHKERR_PETSC(rval);
  rval = moab.create_meshset(MESHSET_SET,meshset_Elastic_stiff_inclusion); CHKERR_PETSC(rval);
  
  
  double density_matrix, density_fibres;
  Range TetsInBlock, TetsInBlock_stiff_inc;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)){
		if(it->get_Cubit_name() == "MAT_MOISTURE_MATRIX") {
			rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
      
      Mat_Moisture mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      density_matrix=mydata.data.Density;

		}
    if(it->get_Cubit_name() == "MAT_MOISTURE_FIBRES") {
			rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock_stiff_inc,true); CHKERR_PETSC(rval);
      
      Mat_Moisture mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      density_fibres=mydata.data.Density;

		}
	}
  cout<<"TetsInBlock "<<TetsInBlock<<endl;
  rval = moab.add_entities(meshset_Elastic,TetsInBlock);CHKERR_PETSC(rval);
  
  cout<<"TetsInBlock "<<TetsInBlock_stiff_inc<<endl;
  rval = moab.add_entities(meshset_Elastic_stiff_inclusion,TetsInBlock_stiff_inc);CHKERR_PETSC(rval);
  
  ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Elastic,"FE_MATRXIX",true); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Elastic_stiff_inclusion,"FE_FIBRES",true); CHKERRQ(ierr);
  
  
  //======================================================================================

  
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
  ierr = mField.partition_problem("MOISTURE_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("MOISTURE_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("MOISTURE_PROBLEM"); CHKERRQ(ierr);
  
  //print block sets with materials
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);
  
  //create matrices (here F, D and Aij are matrices for the full problem)
  Vec F1,F2,F3,C1,C2,C3;
  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",ROW,&F1); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",ROW,&F2); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",ROW,&F3); CHKERRQ(ierr);
  
  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",COL,&C1); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",COL,&C2); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",COL,&C3); CHKERRQ(ierr);
 
  Mat A;
  ierr = mField.MatCreateMPIAIJWithArrays("MOISTURE_PROBLEM",&A); CHKERRQ(ierr);

  
  ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
  applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
  
  ierr = moisture_elements.setThermalFiniteElementLhsOperators("CONC",A); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Disp_Multi_Rhs MyFE_RVELagrange(mField,A,C1,F1,F2,F3,applied_strain,"CONC","LAGRANGE_MUL_FIELD",field_rank);

  ierr = VecZeroEntries(F1); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  ierr = VecZeroEntries(F2); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  ierr = VecZeroEntries(F1); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = MatZeroEntries(A); CHKERRQ(ierr);
  
  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","DIFFUSION_FE",moisture_elements.getLoopFeLhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","LAGRANGE_FE",MyFE_RVELagrange);  CHKERRQ(ierr);

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

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
//  //Matrix View
//  MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//  std::string wait;
//  std::cin >> wait;

  
//Calculation of Homogenized stress
//=========================================================================================================================
  const double young_modulus = 1; //dummy values
  const double poisson_ratio = 0.0;
  
  double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
  Vec RVE_volume_Vec;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
  ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
  
  RVEVolume MyRVEVol(mField,A,C1,F1,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio), RVE_volume_Vec);
//  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","DIFFUSION_FE",MyRVEVol);  CHKERRQ(ierr);
  
  //  ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
//  cout<<"Final RVE_volume = "<< RVE_volume <<endl;
  
  double volume_matrix, volume_fibres;
  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","FE_MATRXIX",MyRVEVol);  CHKERRQ(ierr);
  ierr = VecSum(RVE_volume_Vec, &volume_matrix);  CHKERRQ(ierr);
  cout<<"Matrix volume = "<< volume_matrix <<endl;

  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","FE_FIBRES",MyRVEVol);  CHKERRQ(ierr);
  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
  cout<<"Final RVE_volume = "<< RVE_volume <<endl;
  
  volume_fibres=RVE_volume-volume_matrix;
  cout<<"Fibres volume = "<< volume_fibres <<endl;
  
  cout<<"density_matrix = "<< density_matrix <<endl;
  cout<<"density_fibres = "<< density_fibres <<endl;
  double RVE_density;
  RVE_density=(volume_matrix/RVE_volume)*density_matrix +  (volume_fibres/RVE_volume)*density_fibres;
  
//=========================================================================================================================

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);
  
  ierr = KSPSolve(solver,F1,C1); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(C1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(C1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost("MOISTURE_PROBLEM",ROW,C1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  
  ublas::matrix<FieldData> Dmat;
  Dmat.resize(3,3); Dmat.clear();
  //=============================================================================================================
  // homogenised stress for strian [1 0 0]^T
  //=============================================================================================================
  //create a vector for 3 components of homogenized stress
  Vec Stress_Homo;
  if(pcomm->rank()==0) {
    VecCreateGhost(PETSC_COMM_WORLD,3,3,0,PETSC_NULL,&Stress_Homo);
  } else {
    int ghost[] = {0,1,2};
    VecCreateGhost(PETSC_COMM_WORLD,0,3,3,ghost,&Stress_Homo);
    
  }

  {
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(mField,A,C1,F1,&RVE_volume, applied_strain, Stress_Homo,"CONC","LAGRANGE_MUL_FIELD",field_rank);
    ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","LAGRANGE_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
    VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    PetscScalar *avec;
    VecGetArray(Stress_Homo, &avec);
    for(int ii=0; ii<3; ii++){
      Dmat(ii,0)=*avec;
      avec++;
    }
    
    if(pcomm->rank()==0){
      cout<< "\nStress_Homo = \n\n";
      for(int ii=0; ii<3; ii++){
        cout <<Dmat(ii,0)<<endl;
      }
    }
  }

  //=============================================================================================================
  // homogenised stress for strian [0 1 0]^T
  //=============================================================================================================
  ierr = KSPSolve(solver,F2,C2); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(C2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(C2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost("MOISTURE_PROBLEM",ROW,C2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  {
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(mField,A,C2,F2,&RVE_volume, applied_strain, Stress_Homo,"CONC","LAGRANGE_MUL_FIELD",field_rank);
    ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","LAGRANGE_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
    VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    PetscScalar *avec;
    VecGetArray(Stress_Homo, &avec);
    for(int ii=0; ii<3; ii++){
      Dmat(ii,1)=*avec;
      avec++;
    }
    
    if(pcomm->rank()==0){
      cout<< "\nStress_Homo = \n\n";
      for(int ii=0; ii<3; ii++){
        cout <<Dmat(ii,1)<<endl;
      }
    }
  }

  
  //=============================================================================================================
  // homogenised stress for strian [0 0 1]^T
  //=============================================================================================================
  ierr = KSPSolve(solver,F3,C3); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(C3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(C3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost("MOISTURE_PROBLEM",ROW,C3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  {
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(mField,A,C3,F3,&RVE_volume, applied_strain, Stress_Homo,"CONC","LAGRANGE_MUL_FIELD",field_rank);
    ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","LAGRANGE_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
    VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    PetscScalar *avec;
    VecGetArray(Stress_Homo, &avec);
    for(int ii=0; ii<3; ii++){
      Dmat(ii,2)=*avec;
      avec++;
    }
    
    if(pcomm->rank()==0){
      cout<< "\nStress_Homo = \n\n";
      for(int ii=0; ii<3; ii++){
        cout <<Dmat(ii,2)<<endl;
      }
    }
    
  }
  //====================================================================================================
  //Writing Dmat in to binary file for use in the unsteady diffusion problem
  cout <<"Conductivity matrix="<< Dmat<< endl;
  Dmat=Dmat/RVE_density; //Here we will save Dmat/RVE_density [mm2/s]
  
  cout <<"Effective desnsity RVE ="<< RVE_density<< endl;
  cout <<"Diffusivity matix ="<< Dmat<< endl;
  
//  Dmat=1.0e-6*Dmat; //Here we will save Dmat/RVE_density [m2/s]
  cout <<"Diffusivity After Unit conversoin ="<< Dmat<< endl;

  if(pcomm->rank()==0){
    int fd;
    PetscViewer view_out;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,output_file_Dmat,FILE_MODE_WRITE,&view_out);
    PetscViewerBinaryGetDescriptor(view_out,&fd);
    PetscBinaryWrite(fd,&Dmat(0,0),9,PETSC_DOUBLE,PETSC_FALSE);
    PetscViewerDestroy(&view_out);
  }
  //====================================================================================================


  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet(mField,"CONC",true,false,"CONC");
  ierr = mField.loop_dofs("CONC",ent_method_on_10nodeTet); CHKERRQ(ierr);
  ent_method_on_10nodeTet.set_nodes = false;
  ierr = mField.loop_dofs("CONC",ent_method_on_10nodeTet); CHKERRQ(ierr);
  
  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("MOISTURE_PROBLEM","DIFFUSION_FE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  
   //Destroy matrices
  ierr = VecDestroy(&F1); CHKERRQ(ierr);
  ierr = VecDestroy(&F2); CHKERRQ(ierr);
  ierr = VecDestroy(&F2); CHKERRQ(ierr);

  ierr = VecDestroy(&C1); CHKERRQ(ierr);
  ierr = VecDestroy(&C2); CHKERRQ(ierr);
  ierr = VecDestroy(&C3); CHKERRQ(ierr);

  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  
  
  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  
  PetscFinalize();
  
}

