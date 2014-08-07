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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "FEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "ForcesAndSurcesCore.hpp"
#include "TsCtx.hpp"
#include "MoistureElement.hpp"

#include "ElasticFE_RVELagrange_Traction.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Traction.hpp"
#include "ElasticFE_RVELagrange_RigidBodyTranslation.hpp"
#include "RVEVolume.hpp"

#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include "Projection10NodeCoordsOnField.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {
  
  PetscInitialize(&argc,&argv,(char *)0,help);
  
  Core mb_instance;
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
  
  //Applied strain on the RVE (vector of length 6) strain=[xx, yy, zz, xy, xz, zy]^T
  double myapplied_strain[3];
  int nmax=3;
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-myapplied_strain",myapplied_strain,&nmax,&flg); CHKERRQ(ierr);
  ublas::vector<FieldData> applied_strain;
  applied_strain.resize(3);
  cblas_dcopy(3, &myapplied_strain[0], 1, &applied_strain(0), 1);
  cout<<"applied_strain ="<<applied_strain<<endl;
  
  
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
  FieldCore core(moab);
  FieldInterface& mField = core;
  
  
  //=======================================================================================================
  //Seting nodal coordinates on the surface to make sure they are periodic
  //=======================================================================================================
  
  Range SurTrisNeg, SurTrisPos;
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,SIDESET,2,SurTrisNeg,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 101 = %d\n",SurTrisNeg.size()); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(102,SIDESET,2,SurTrisPos,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 102 = %d\n",SurTrisPos.size()); CHKERRQ(ierr);
  
  Range SurNodesNeg,SurNodesPos;
  rval = moab.get_connectivity(SurTrisNeg,SurNodesNeg,true); CHKERR_PETSC(rval);
  cout<<" All nodes on negative surfaces " << SurNodesNeg.size()<<endl;
  rval = moab.get_connectivity(SurTrisPos,SurNodesPos,true); CHKERR_PETSC(rval);
  cout<<" All nodes on positive surfaces " << SurNodesPos.size()<<endl;
  
  
  double roundfact=1000.0;   double coords_nodes[3];
  //Populating the Multi-index container with nodes on -ve faces
  for(Range::iterator nit = SurNodesNeg.begin(); nit!=SurNodesNeg.end();  nit++) {
    rval = moab.get_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
    //round values to 3 disimal places
    if(coords_nodes[0]>=0) coords_nodes[0]=double(int(coords_nodes[0]*roundfact+0.5))/roundfact;  else coords_nodes[0]=double(int(coords_nodes[0]*roundfact-0.5))/roundfact;
    if(coords_nodes[1]>=0) coords_nodes[1]=double(int(coords_nodes[1]*roundfact+0.5))/roundfact;  else coords_nodes[1]=double(int(coords_nodes[1]*roundfact-0.5))/roundfact;
    if(coords_nodes[2]>=0) coords_nodes[2]=double(int(coords_nodes[2]*roundfact+0.5))/roundfact;  else coords_nodes[2]=double(int(coords_nodes[2]*roundfact-0.5))/roundfact;
    rval = moab.set_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
    //      cout<<"   coords_nodes[0]= "<<coords_nodes[0] << "   coords_nodes[1]= "<< coords_nodes[1] << "   coords_nodes[2]= "<< coords_nodes[2] <<endl;
  }
  
  ///Populating the Multi-index container with nodes on +ve faces
  for(Range::iterator nit = SurNodesPos.begin(); nit!=SurNodesPos.end();  nit++) {
    rval = moab.get_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
    //round values to 3 disimal places
    if(coords_nodes[0]>=0) coords_nodes[0]=double(int(coords_nodes[0]*roundfact+0.5))/roundfact;  else coords_nodes[0]=double(int(coords_nodes[0]*roundfact-0.5))/roundfact;
    if(coords_nodes[1]>=0) coords_nodes[1]=double(int(coords_nodes[1]*roundfact+0.5))/roundfact;  else coords_nodes[1]=double(int(coords_nodes[1]*roundfact-0.5))/roundfact;
    if(coords_nodes[2]>=0) coords_nodes[2]=double(int(coords_nodes[2]*roundfact+0.5))/roundfact;  else coords_nodes[2]=double(int(coords_nodes[2]*roundfact-0.5))/roundfact;
    rval = moab.set_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
    //      cout<<"   coords_nodes[0]= "<<coords_nodes[0] << "   coords_nodes[1]= "<< coords_nodes[1] << "   coords_nodes[2]= "<< coords_nodes[2] <<endl;
  }
  //=======================================================================================================
  
  
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
  int field_rank=1;
  ierr = mField.add_field("CONC",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("LAGRANGE_MUL_FIELD",NOFIELD,3); CHKERRQ(ierr);
  ierr = mField.add_field("LAGRANGE_MUL_FIELD_RIGID_TRANS",NOFIELD,1); CHKERRQ(ierr);   //Control 3 rigid body translations in x, y and z axis
  
  //Problem
  ierr = mField.add_problem("MOISTURE_PROBLEM"); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"CONC"); CHKERRQ(ierr);

  //FE
  MoistureElement moisture_elements(mField);
  ierr = moisture_elements.addMoistureElements("MOISTURE_PROBLEM","CONC"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("LAGRANGE_FE"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("LAGRANGE_FE_RIGID_TRANS"); CHKERRQ(ierr);

  //======================================================================================================
  //C row as LAGRANGE_MUL_FIELD and col as CONC
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE","LAGRANGE_MUL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE","CONC"); CHKERRQ(ierr);
  
  //CT col as LAGRANGE_MUL_FIELD and row as CONC
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE","LAGRANGE_MUL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE","CONC"); CHKERRQ(ierr);
  
  //As for stress we need both displacement and temprature (Lukasz)
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE","LAGRANGE_MUL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE","CONC"); CHKERRQ(ierr);
  //======================================================================================================


  //Define rows/cols and element data for C1 and C1T (for lagrange multipliers to contol the rigid body translations)
  //============================================================================================================
  //C1 row as Lagrange_elem_rigid_trans and col as CONC
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE_RIGID_TRANS","LAGRANGE_MUL_FIELD_RIGID_TRANS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE_RIGID_TRANS","CONC"); CHKERRQ(ierr);
  
  //C1T col as Lagrange_elem_rigid_trans and row as CONC
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE_RIGID_TRANS","LAGRANGE_MUL_FIELD_RIGID_TRANS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE_RIGID_TRANS","CONC"); CHKERRQ(ierr);
  
  //As for stress we need both CONC and temprature (Lukasz)
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE_RIGID_TRANS","LAGRANGE_MUL_FIELD_RIGID_TRANS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE_RIGID_TRANS","CONC"); CHKERRQ(ierr);
  //============================================================================================================
  

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("MOISTURE_PROBLEM","LAGRANGE_FE"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("MOISTURE_PROBLEM","LAGRANGE_FE_RIGID_TRANS"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("MOISTURE_PROBLEM",bit_level0); CHKERRQ(ierr);

  /***/
  //Declare problem
  
  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"CONC",2); CHKERRQ(ierr);
  
  Range SurfacesFaces;
  ierr = mField.get_Cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"LAGRANGE_FE"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"LAGRANGE_FE_RIGID_TRANS"); CHKERRQ(ierr);

  
  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
	rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);
  
  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  ierr = mField.set_field_order(0,MBTET,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"CONC",1); CHKERRQ(ierr);
  
  
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
  Vec F,C;
  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",ROW,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",COL,&C); CHKERRQ(ierr);
  
  Mat A;
  ierr = mField.MatCreateMPIAIJWithArrays("MOISTURE_PROBLEM",&A); CHKERRQ(ierr);
  

  //Assemble F and A
  ierr = moisture_elements.setMoistureFiniteElementLhsOperators("CONC",(&A)); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Traction MyFE_RVELagrange(mField,A,C,F,applied_strain,"CONC","LAGRANGE_MUL_FIELD",field_rank);
  ElasticFE_RVELagrange_RigidBodyTranslation MyFE_RVELagrangeRigidBodyTrans(mField,A,C,F,applied_strain,"CONC","LAGRANGE_MUL_FIELD",field_rank,"LAGRANGE_MUL_FIELD_RIGID_TRANS");

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);
  
  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","MOISTURE_FE",moisture_elements.getLoopFeLhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","LAGRANGE_FE",MyFE_RVELagrange);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","LAGRANGE_FE_RIGID_TRANS",MyFE_RVELagrangeRigidBodyTrans);  CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  
//    //Matrix View
//    MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//    std::string wait;
//    std::cin >> wait;
  
  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);
  
  ierr = KSPSolve(solver,F,C); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  //  ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //  ierr = VecView(C,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  
  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("MOISTURE_PROBLEM",ROW,C,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  //Calculation of Homogenized stress
  //=======================================================================================================================================================
  const double young_modulus = 1;
  const double poisson_ratio = 0.0;
  
  double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
  Vec RVE_volume_Vec;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
  ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
  
  RVEVolume MyRVEVol(mField,A,C,F,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio), RVE_volume_Vec);
  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","MOISTURE_FE",MyRVEVol);  CHKERRQ(ierr);
  //  ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
  cout<<"Final RVE_volume = "<< RVE_volume <<endl;
  
  
  //create a vector for 3 components of homogenized stress
  Vec Stress_Homo;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 3, 3*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
  
  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Traction MyFE_RVEHomoStressDisp(mField,A,C,F,&RVE_volume, applied_strain, Stress_Homo,"CONC","LAGRANGE_MUL_FIELD",field_rank);
  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","LAGRANGE_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
  
  //  if(pcomm->rank()) cout<< " Stress_Homo =  "<<endl;
  //  ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec;
    VecGetArray(Stress_Homo, &avec);
    
    cout<< "\nStress_Homo = \n\n";
    for(int ii=0; ii<3; ii++){
      cout <<*avec<<endl; ;
      avec++;
    }
  }
  cout<< "\n\n";
  
  //=======================================================================================================================================================
  
  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet(mField,"CONC",true,false,"CONC");
  ierr = mField.loop_dofs("CONC",ent_method_on_10nodeTet); CHKERRQ(ierr);
  ent_method_on_10nodeTet.set_nodes = false;
  ierr = mField.loop_dofs("CONC",ent_method_on_10nodeTet); CHKERRQ(ierr);
  
  
  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("MOISTURE_PROBLEM","MOISTURE_FE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  
  //PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
  PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method(mField,"CONC",LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","MOISTURE_FE",fe_post_proc_method);  CHKERRQ(ierr);
  
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }
  
  
  //Destroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&C); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  
  
  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  
  PetscFinalize();
  
}

