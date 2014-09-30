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

#include "Coupled_MechFEMethod.hpp"
#include "MoistureFEMethod.hpp"
#include "Coupled_MechMoistureFEMethod.hpp"

#include "ElasticFE_RVELagrange_Disp.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Disp.hpp"
#include "RVEVolume.hpp"

#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include "FEMethod_DriverComplexForLazy.hpp"

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
  
  
  //Mechanical
  //Applied strain on the RVE (vector of length 6) strain=[xx, yy, zz, xy, xz, zy]^T
  double myapplied_strain[6];
  int nmax=6;
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-myapplied_strain",myapplied_strain,&nmax,&flg); CHKERRQ(ierr);
  ublas::vector<FieldData> applied_strain;
  applied_strain.resize(6);
  cblas_dcopy(6, &myapplied_strain[0], 1, &applied_strain(0), 1);
  cout<<"applied_strain ="<<applied_strain<<endl;
  
  //Moisture-transport
  //Applied concentration gradient on the RVE (vector of length 3) strain=[xx, yy, zz]^T
  double myapplied_congrad[3];
  int nmax1=3;
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-myapplied_congrad",myapplied_congrad,&nmax1,&flg); CHKERRQ(ierr);
  ublas::vector<FieldData> applied_congrad;
  applied_congrad.resize(3);
  cblas_dcopy(3, &myapplied_congrad[0], 1, &applied_congrad(0), 1);
  cout<<"applied_congradient ="<<applied_congrad<<endl;

  
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
  int field_rank_mech=3;
  ierr = mField.add_field("DISPLACEMENT",H1,field_rank_mech); CHKERRQ(ierr);
  ierr = mField.add_field("Lagrange_mul_disp",H1,field_rank_mech); CHKERRQ(ierr);
  
  int field_rank_mois=1;
  ierr = mField.add_field("CONC",H1,field_rank_mois); CHKERRQ(ierr);
  ierr = mField.add_field("Lagrange_mul_conc",H1,field_rank_mois); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("Kuu"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("Kuc"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("Kcc"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("Lagrange_elm_disp"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("Lagrange_elm_conc"); CHKERRQ(ierr);

  //Define rows/cols and element data for Kuu
  ierr = mField.modify_finite_element_add_field_row("Kuu","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("Kuu","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Kuu","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Kuu","CONC"); CHKERRQ(ierr);

  //Define rows/cols and element data for Kcc
  ierr = mField.modify_finite_element_add_field_row("Kcc","CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("Kcc","CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Kcc","CONC"); CHKERRQ(ierr);

//  //Define rows/cols and element data for Kuc
  ierr = mField.modify_finite_element_add_field_row("Kuc","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("Kuc","CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Kuc","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Kuc","CONC"); CHKERRQ(ierr);

  
  //Mechanical (Cu and CuT)
  //=====================================================================================================
  //Cu -> row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elm_disp","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elm_disp","DISPLACEMENT"); CHKERRQ(ierr);
  
  //CuT -> col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elm_disp","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elm_disp","DISPLACEMENT"); CHKERRQ(ierr);
  
  //As for stress we need both displacement and temprature
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elm_disp","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elm_disp","DISPLACEMENT"); CHKERRQ(ierr);
  //=====================================================================================================
  
  //Moisture (Cu and CuT)
  //=====================================================================================================
  //Cc -> row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elm_conc","Lagrange_mul_conc"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elm_conc","CONC"); CHKERRQ(ierr);
  
  //CcT -> col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elm_conc","Lagrange_mul_conc"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elm_conc","CONC"); CHKERRQ(ierr);
  
  //As for stress we need both displacement and temprature
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elm_conc","Lagrange_mul_conc"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elm_conc","CONC"); CHKERRQ(ierr);
  //=====================================================================================================

  
  //define problems
  ierr = mField.add_problem("COUPLED_MECH_MOIS"); CHKERRQ(ierr);
  
  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("COUPLED_MECH_MOIS","Kuu"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_MECH_MOIS","Kuc"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_MECH_MOIS","Kcc"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_MECH_MOIS","Lagrange_elm_disp"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_MECH_MOIS","Lagrange_elm_conc"); CHKERRQ(ierr);

  
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("COUPLED_MECH_MOIS",bit_level0); CHKERRQ(ierr);
  
  /***/
  //Declare problem
  
  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT",2); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"CONC",2); CHKERRQ(ierr);

  
  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"Kuu",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"Kuc",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"Kcc",MBTET); CHKERRQ(ierr);

  Range SurfacesFaces;
  ierr = mField.get_Cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elm_disp"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elm_conc"); CHKERRQ(ierr);

  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
	rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Lagrange_mul_disp",2); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Lagrange_mul_conc",2); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(0,MBTET,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"CONC",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTRI,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(0,MBTRI,"Lagrange_mul_conc",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"Lagrange_mul_conc",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"Lagrange_mul_conc",1); CHKERRQ(ierr);

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
  ierr = mField.partition_problem("COUPLED_MECH_MOIS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("COUPLED_MECH_MOIS"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("COUPLED_MECH_MOIS"); CHKERRQ(ierr);
  
  //print block sets with materials
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);
  
  //create matrices (here F, D and Aij are matrices for the full problem)
  Vec F,D;
  ierr = mField.VecCreateGhost("COUPLED_MECH_MOIS",ROW,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("COUPLED_MECH_MOIS",COL,&D); CHKERRQ(ierr);
  
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("COUPLED_MECH_MOIS",&Aij); CHKERRQ(ierr);
  
  //Assemble F and Aij
//  SpatialPositionsBCFEMethodPreAndPostProc MyDirichletBC(mField,"DISPLACEMENT",Aij,D,F);
  SnesCtx snes_ctx(mField,"COUPLED_MECH_MOIS");
  const double young_modulus = 1;
  const double poisson_ratio = 0.0;
    
  Coupled_MechFEMethod my_fe_mech(mField,Aij,D,F,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  MoistureFEMethod my_fe_mois(mField,Aij,D,F);
  Coupled_MechMoistureFEMethod my_fe_coupled_mechmois(mField,Aij,D,F);
  

  ElasticFE_RVELagrange_Disp MyFE_RVELagrange_mech(mField,Aij,D,F,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank_mech);
  ElasticFE_RVELagrange_Disp MyFE_RVELagrange_mois(mField,Aij,D,F,applied_congrad,"CONC","Lagrange_mul_conc",field_rank_mois);

  
//  //*********************************************************************************************************
//  //to solve linear problem
//  //*********************************************************************************************************
//  ierr = VecZeroEntries(F); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
//  
//  ierr = mField.loop_finite_elements("COUPLED_MECH_MOIS","Kuu",my_fe_mech);  CHKERRQ(ierr);
////  ierr = mField.loop_finite_elements("COUPLED_MECH_MOIS","Kcc",my_fe_mois);  CHKERRQ(ierr);
//  ierr = mField.loop_finite_elements("COUPLED_MECH_MOIS","Lagrange_elm_disp",MyFE_RVELagrange_mech);  CHKERRQ(ierr);
////  ierr = mField.loop_finite_elements("COUPLED_MECH_MOIS","Lagrange_elm_conc",MyFE_RVELagrange_mois);  CHKERRQ(ierr);
//
//  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
//  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
//  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//  
//  PetscSynchronizedFlush(PETSC_COMM_WORLD);
//  //Matrix View
//  MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//  std::string wait;
//  std::cin >> wait;
//
//  //Solver
//  KSP solverM;
//  ierr = KSPCreate(PETSC_COMM_WORLD,&solverM); CHKERRQ(ierr);
//  ierr = KSPSetOperators(solverM,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
//  ierr = KSPSetFromOptions(solverM); CHKERRQ(ierr);
//  ierr = KSPSetUp(solverM); CHKERRQ(ierr);
//  
//  ierr = KSPSolve(solverM,F,D); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  
//  //Save data on mesh
//  ierr = mField.set_global_VecCreateGhost("COUPLED_MECH_MOIS",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  //*********************************************************************************************************
  
  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,Aij,Aij,SnesMat,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  SnesCtx::loops_to_do_type& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
//  snes_ctx.get_preProcess_to_do_Rhs().push_back(&MyDirichletBC);
//  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("Kuu",&my_fe_mech));  // we are calculting this already in Kuc
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("Kcc",&my_fe_mois));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("Kuc",&my_fe_coupled_mechmois));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("Lagrange_elm_disp",&MyFE_RVELagrange_mech));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("Lagrange_elm_conc",&MyFE_RVELagrange_mois));
//  snes_ctx.get_postProcess_to_do_Rhs().push_back(&MyDirichletBC);

  SnesCtx::loops_to_do_type& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
//  snes_ctx.get_preProcess_to_do_Mat().push_back(&MyDirichletBC);
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("Kuu",&my_fe_mech));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("Kcc",&my_fe_mois));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("Kuc",&my_fe_coupled_mechmois));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("Lagrange_elm_disp",&MyFE_RVELagrange_mech));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("Lagrange_elm_conc",&MyFE_RVELagrange_mois));
//  snes_ctx.get_postProcess_to_do_Mat().push_back(&MyDirichletBC);

  ierr = mField.set_local_VecCreateGhost("COUPLED_MECH_MOIS",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("COUPLED_MECH_MOIS",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);


  PostProcVertexMethod ent_method(mField.get_moab(),"DISPLACEMENT");
  ierr = mField.loop_dofs("COUPLED_MECH_MOIS","DISPLACEMENT",ROW,ent_method); CHKERRQ(ierr);
  PostProcVertexMethod ent_method1(mField.get_moab(),"CONC");
  ierr = mField.loop_dofs("CONC",ent_method1); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("COUPLED_MECH_MOIS","Kuu",out_meshset); CHKERRQ(ierr);
    rval = mField.get_moab().write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  
  //Calculation of Homogenized stress
  //=======================================================================================================================================================
  double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
  Vec RVE_volume_Vec;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
  ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
  
  RVEVolume MyRVEVol(mField,Aij,D,F,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio), RVE_volume_Vec);
  ierr = mField.loop_finite_elements("COUPLED_MECH_MOIS","Kuu",MyRVEVol);  CHKERRQ(ierr);
  //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
    cout<<"Final RVE_volume = "<< RVE_volume <<endl;
  
  
  //create a vector for 6 components of homogenized stress
  Vec Stress_Homo;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
  
  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank_mech);
  ierr = mField.loop_finite_elements("COUPLED_MECH_MOIS","Lagrange_elm_disp",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
  
//  if(pcomm->rank()) cout<< " Stress_Homo =  "<<endl;
//  ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec;
    VecGetArray(Stress_Homo, &avec);
    
    cout<< "\nStress_Homo = \n\n";
    for(int ii=0; ii<6; ii++){
      cout <<*avec<<endl; ;
      avec++;
    }
  }
  cout<< "\n\n";

  //=======================================================================================================================================================

  
//  //PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
//  PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method(mField,"DISPLACEMENT",LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
//  ierr = mField.loop_finite_elements("COUPLED_MECH_MOIS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
//  
//  if(pcomm->rank()==0) {
//    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
//  }
//  
//  //End Disp
//  /*Range ents;*/
//  //ierr = mField.get_Cubit_msId_entities_by_dimension(1,NodeSet,0,ents,true); CHKERRQ(ierr);
//  //for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"DISPLACEMENT",dit)) {
//  //if(find(ents.begin(),ents.end(),dit->get_ent())!=ents.end()) {
//  //PetscSynchronizedPrintf(PETSC_COMM_WORLD, "val = %6.7e\n",dit->get_FieldData());
//  //}
//  /*}*/
//  
//  //Support stresses
//  //Range ents;
//  //ierr = mField.get_Cubit_msId_entities_by_dimension(4,NodeSet,0,ents,true); CHKERRQ(ierr);
//  
//  //Destroy matrices
//  ierr = VecDestroy(&F); CHKERRQ(ierr);
//  ierr = VecDestroy(&D); CHKERRQ(ierr);
//  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
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

