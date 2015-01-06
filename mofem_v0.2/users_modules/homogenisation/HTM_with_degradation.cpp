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

using namespace boost::numeric;
using namespace ObosleteUsersModules;

#include <ThermalElement.hpp>
#include <MoistureElement.hpp>
#include <ElasticFEMethod.hpp>

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
  
  
  //Mechanical
  //Applied strain on the RVE (vector of length 6) strain=[xx, yy, zz, xy, xz, zy]^T
  double myapplied_strain[6];
  int nmax=6;
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-myapplied_strain",myapplied_strain,&nmax,&flg); CHKERRQ(ierr);
  ublas::vector<FieldData> applied_strain;
  applied_strain.resize(6);
  cblas_dcopy(6, &myapplied_strain[0], 1, &applied_strain(0), 1);
  cout<<"applied_strain ="<<applied_strain<<endl;
  
  //Thermal
  //Applied temprature gradient on the RVE (vector of length 3) gradT=[xx, yy, zz]^T
  double myapplied_Tgrad[3];
  int nmaxT=3;
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-myapplied_Tgrad",myapplied_Tgrad,&nmaxT,&flg); CHKERRQ(ierr);
  ublas::vector<FieldData> applied_Tgrad;
  applied_Tgrad.resize(3);
  cblas_dcopy(3, &myapplied_Tgrad[0], 1, &applied_Tgrad(0), 1);
  cout<<"applied_Tgrad ="<<applied_Tgrad<<endl;

  //Moisture-transport
  //Applied concentration gradient on the RVE (vector of length 3) strain=[xx, yy, zz]^T
  double myapplied_cgrad[3];
  int nmaxc=3;
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-myapplied_congrad",myapplied_cgrad,&nmaxc,&flg); CHKERRQ(ierr);
  ublas::vector<FieldData> applied_cgrad;
  applied_cgrad.resize(3);
  cblas_dcopy(3, &myapplied_cgrad[0], 1, &applied_cgrad(0), 1);
  cout<<"applied_cgradient ="<<applied_cgrad<<endl;

  
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
  
//  //Fields
  int field_rank_TC=1;
  ierr = mField.add_field("Field_Temp",H1,field_rank_TC); CHKERRQ(ierr);
  ierr = mField.add_field("Field_LagMul_Temp",H1,field_rank_TC); CHKERRQ(ierr);

  ierr = mField.add_field("Field_Conc",H1,field_rank_TC); CHKERRQ(ierr);
  ierr = mField.add_field("Field_LagMul_Conc",H1,field_rank_TC); CHKERRQ(ierr);
  
  int field_rank_mech=3;
  ierr = mField.add_field("Field_Disp",H1,field_rank_mech); CHKERRQ(ierr);
  ierr = mField.add_field("Field_LagMul_Disp",H1,field_rank_mech); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("PROB_THERMAL"); CHKERRQ(ierr);
  ierr = mField.add_problem("PROB_MOIS"); CHKERRQ(ierr);
  ierr = mField.add_problem("PROB_MECH"); CHKERRQ(ierr);
  
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("PROB_THERMAL",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("PROB_MOIS",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("PROB_MECH",bit_level0); CHKERRQ(ierr);
  
  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"Field_Temp"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"Field_Conc"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"Field_Disp"); CHKERRQ(ierr);


  ThermalElement thermal_elements(mField);
  ierr = thermal_elements.addThermalElements("PROB_THERMAL","Field_Temp"); CHKERRQ(ierr);

  MoistureElement moisture_elements(mField);
  ierr = moisture_elements.addMoistureElements("PROB_MOIS","Field_Conc"); CHKERRQ(ierr);

  
  //FE
  ierr = mField.add_finite_element("THERMAL_FE_LagMul"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("MOISTURE_FE_LagMul"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("MECH_FE"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("MECH_FE_LagMul"); CHKERRQ(ierr);

  //Define rows/cols for THERMAL_FE_LagMul element
  //=============================================================================================
   //CT
  ierr = mField.modify_finite_element_add_field_row("THERMAL_FE_LagMul","Field_LagMul_Temp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("THERMAL_FE_LagMul","Field_Temp"); CHKERRQ(ierr);

  //CT^T
  ierr = mField.modify_finite_element_add_field_row("THERMAL_FE_LagMul","Field_Temp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("THERMAL_FE_LagMul","Field_LagMul_Temp"); CHKERRQ(ierr);

  //data
  ierr = mField.modify_finite_element_add_field_data("THERMAL_FE_LagMul","Field_LagMul_Temp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("THERMAL_FE_LagMul","Field_Temp"); CHKERRQ(ierr);
  //=============================================================================================

  
  //Define rows/cols for MOISTURE_FE_LagMul element
  //=============================================================================================
  //Cc
  ierr = mField.modify_finite_element_add_field_row("MOISTURE_FE_LagMul","Field_LagMul_Conc"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MOISTURE_FE_LagMul","Field_Conc"); CHKERRQ(ierr);
  
  //Cc^T
  ierr = mField.modify_finite_element_add_field_row("MOISTURE_FE_LagMul","Field_Conc"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MOISTURE_FE_LagMul","Field_LagMul_Conc"); CHKERRQ(ierr);
  
  //data
  ierr = mField.modify_finite_element_add_field_data("MOISTURE_FE_LagMul","Field_Conc"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MOISTURE_FE_LagMul","Field_LagMul_Conc"); CHKERRQ(ierr);
  //=============================================================================================

  
  //Define rows/cols and data for moisture problem
  //=============================================================================================
  ierr = mField.modify_finite_element_add_field_row("MECH_FE","Field_Disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MECH_FE","Field_Disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MECH_FE","Field_Disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MECH_FE","Field_Temp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MECH_FE","Field_Conc"); CHKERRQ(ierr);

  //Cu
  ierr = mField.modify_finite_element_add_field_row("MECH_FE_LagMul","Field_LagMul_Disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MECH_FE_LagMul","Field_Disp"); CHKERRQ(ierr);
  
  //Cu^T
  ierr = mField.modify_finite_element_add_field_row("MECH_FE_LagMul","Field_Disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MECH_FE_LagMul","Field_LagMul_Disp"); CHKERRQ(ierr);
  
  //data
  ierr = mField.modify_finite_element_add_field_data("MECH_FE_LagMul","Field_Disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MECH_FE_LagMul","Field_LagMul_Disp"); CHKERRQ(ierr);
  //=============================================================================================

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("PROB_THERMAL","THERMAL_FE_LagMul"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("PROB_MOIS","MOISTURE_FE_LagMul"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("PROB_MECH","MECH_FE"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("PROB_MECH","MECH_FE_LagMul"); CHKERRQ(ierr);

//  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MECH_FE",MBTET); CHKERRQ(ierr);

  Range SurfacesFaces;
  ierr = mField.get_Cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"THERMAL_FE_LagMul"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"MOISTURE_FE_LagMul"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"MECH_FE_LagMul"); CHKERRQ(ierr);

  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
	rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Field_LagMul_Temp",2); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Field_LagMul_Conc",2); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Field_LagMul_Disp",2); CHKERRQ(ierr);

  
  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  ierr = mField.set_field_order(0,MBTET,"Field_Temp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"Field_Temp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"Field_Temp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"Field_Temp",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(0,MBTET,"Field_Conc",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"Field_Conc",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"Field_Conc",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"Field_Conc",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(0,MBTET,"Field_Disp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"Field_Disp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"Field_Disp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"Field_Disp",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(0,MBTRI,"Field_LagMul_Temp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"Field_LagMul_Temp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"Field_LagMul_Temp",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(0,MBTRI,"Field_LagMul_Conc",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"Field_LagMul_Conc",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"Field_LagMul_Conc",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(0,MBTRI,"Field_LagMul_Disp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"Field_LagMul_Disp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"Field_LagMul_Disp",1); CHKERRQ(ierr);
  
  

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
  //partition PROB_THERMAL
  ierr = mField.partition_problem("PROB_THERMAL"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("PROB_THERMAL"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("PROB_THERMAL"); CHKERRQ(ierr);
  
  //partition PROB_MOIS
  ierr = mField.partition_problem("PROB_MOIS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("PROB_MOIS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("PROB_MOIS"); CHKERRQ(ierr);
  
  //partition PROB_MECH
  ierr = mField.partition_problem("PROB_MECH"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("PROB_MECH"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("PROB_MECH"); CHKERRQ(ierr);

  //print block sets with materials
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);
  
  //create matrices and vectors
  Vec FT,DT,Fc,Dc,FM,DM;
  ierr = mField.VecCreateGhost("PROB_THERMAL",ROW,&FT); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("PROB_THERMAL",COL,&DT); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("PROB_MOIS",ROW,&Fc); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("PROB_MOIS",COL,&Dc); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("PROB_MECH",ROW,&FM); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("PROB_MECH",COL,&DM); CHKERRQ(ierr);

  Mat AT,Ac,AM;
  ierr = mField.MatCreateMPIAIJWithArrays("PROB_THERMAL",&AT); CHKERRQ(ierr);
  ierr = mField.MatCreateMPIAIJWithArrays("PROB_MOIS",&Ac); CHKERRQ(ierr);
  ierr = mField.MatCreateMPIAIJWithArrays("PROB_MECH",&AM); CHKERRQ(ierr);

  //***************************************************************************************************
  //Solve thermal problem to calculate temprature distribution
  //***************************************************************************************************

  ierr = thermal_elements.setThermalFiniteElementLhsOperators("Field_Temp",AT); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Disp MyFE_RVELagrange_thermal(mField,AT,DT,FT,applied_Tgrad,"Field_Temp","Field_LagMul_Temp",field_rank_TC);
  
  ierr = VecZeroEntries(FT); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(FT,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(FT,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(AT); CHKERRQ(ierr);
  
  ierr = mField.loop_finite_elements("PROB_THERMAL","THERMAL_FE",thermal_elements.getLoopFeLhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("PROB_THERMAL","THERMAL_FE_LagMul",MyFE_RVELagrange_thermal);  CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(FT,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(FT,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(FT); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(FT); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(AT,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(AT,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
//  //Matrix View
//  ierr = VecView(FT,PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
//  ierr = MatView(AT,PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
//  std::string wait;
//  std::cin >> wait;
  
  //Solver
  KSP solver_thermal;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver_thermal); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver_thermal,AT,AT); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver_thermal); CHKERRQ(ierr);
  ierr = KSPSetUp(solver_thermal); CHKERRQ(ierr);
  
  ierr = KSPSolve(solver_thermal,FT,DT); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(DT,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(DT,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("PROB_THERMAL",ROW,DT,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(DT,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  //***************************************************************************************************
  //Solve moisture problem to calculate moisture concentration
  //***************************************************************************************************

  
  ierr = moisture_elements.setMoistureFiniteElementLhsOperators("Field_Conc",Ac); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Disp MyFE_RVELagrange_moisture(mField,Ac,Dc,Fc,applied_cgrad,"Field_Conc","Field_LagMul_Conc",field_rank_TC);
  
  ierr = VecZeroEntries(Fc); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(Fc,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(Fc,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Ac); CHKERRQ(ierr);
  
  ierr = mField.loop_finite_elements("PROB_MOIS","MOISTURE_FE",moisture_elements.getLoopFeLhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("PROB_MOIS","MOISTURE_FE_LagMul",MyFE_RVELagrange_moisture);  CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(Fc,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(Fc,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(Fc); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Fc); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Ac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Ac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
//  //Matrix View
//  ierr = VecView(Fc,PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
//  ierr = MatView(Ac,PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
//  std::string wait;
//  std::cin >> wait;

  //Solver
  KSP solver_moisture;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver_moisture); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver_moisture,Ac,Ac); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver_moisture); CHKERRQ(ierr);
  ierr = KSPSetUp(solver_moisture); CHKERRQ(ierr);
  
  ierr = KSPSolve(solver_moisture,Fc,Dc); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(Dc,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(Dc,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("PROB_MOIS",ROW,Dc,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(Dc,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  
  
  //***************************************************************************************************
  //Mechanical problem
  //***************************************************************************************************
  struct MyElasticFEMethod: public ElasticFEMethod {
    MyElasticFEMethod(FieldInterface& _m_field,Mat _Aij,Vec _D,Vec& _F,double _lambda,double _mu):
    ElasticFEMethod(_m_field,_Aij,_D,_F,_lambda,_mu) {};
    
    PetscErrorCode Fint(Vec F_int) {
      PetscFunctionBegin;
      ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        if(RowGlob[rr].size()==0) continue;
        f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
        ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }
    
  };

  
  MyElasticFEMethod my_fe_mech(mField,AM,DM,FM,0,0);
//  ElasticFE_RVELagrange_Disp MyFE_RVELagrange_mech(mField,AM,DM,FM,applied_strain,"Field_Disp","Field_LagMul_Disp",field_rank);
//  ierr = VecZeroEntries(FM); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(FM,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(FM,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = MatZeroEntries(AM); CHKERRQ(ierr);
//  
//  ierr = mField.loop_finite_elements("PROB_MECH","MECH_FE",my_fe_mech);  CHKERRQ(ierr);
//  ierr = mField.loop_finite_elements("PROB_MECH","MECH_FE_LagMul",MyFE_RVELagrange_mech);  CHKERRQ(ierr);

//  MoistureFEMethod my_fe_mois(mField,Ac,Dc,Fc);
//  ElasticFE_RVELagrange_Disp MyFE_RVELagrange_mois(mField,Ac,Dc,Fc,applied_congrad,"CONC","Lagrange_mul_conc",field_rank_mois);
//  ierr = VecZeroEntries(Fc); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(Fc,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(Fc,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = MatZeroEntries(Ac); CHKERRQ(ierr);
//
//  ierr = mField.loop_finite_elements("PROB_MOIS","Kcc",my_fe_mois);  CHKERRQ(ierr);
//  ierr = mField.loop_finite_elements("PROB_MOIS","Lagrange_elm_conc",MyFE_RVELagrange_mois);  CHKERRQ(ierr);
//
//  ierr = VecGhostUpdateBegin(Fc,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(Fc,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecAssemblyBegin(Fc); CHKERRQ(ierr);
//  ierr = VecAssemblyEnd(Fc); CHKERRQ(ierr);
//  ierr = MatAssemblyBegin(Ac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//  ierr = MatAssemblyEnd(Ac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//
//  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  //ierr = MatView(Aij,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//
////  //Matrix View
////  MatView(Ac,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
////  std::string wait;
////  std::cin >> wait;
//
//  //Solver
//  KSP solverC;
//  ierr = KSPCreate(PETSC_COMM_WORLD,&solverC); CHKERRQ(ierr);
//  ierr = KSPSetOperators(solverC,Ac,Ac); CHKERRQ(ierr);
//  ierr = KSPSetFromOptions(solverC); CHKERRQ(ierr);
//  ierr = KSPSetUp(solverC); CHKERRQ(ierr);
//
//  ierr = KSPSolve(solverC,Fc,Dc); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(Dc,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(Dc,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//
//  //Save data on mesh
//  ierr = mField.set_global_VecCreateGhost("PROB_MOIS",ROW,Dc,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
////  ierr = VecView(Dc,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  
//
//  struct MyElasticFEMethod: public Coupled_MechFEMethod {
//    MyElasticFEMethod(FieldInterface& _mField, Mat &_Aij,Vec &_D,Vec& _F,double _lambda,double _mu):
//    Coupled_MechFEMethod(_mField,_Aij,_D,_F,_lambda,_mu) {};
//    
//    PetscErrorCode Fint(Vec F_int) {
//      PetscFunctionBegin;
//      ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
//      for(int rr = 0;rr<row_mat;rr++) {
//        if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
//        if(RowGlob[rr].size()==0) continue;
//        f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
//        ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
//      }
//      PetscFunctionReturn(0);
//    }
//    
//  };
//
//  //Assemble F and Aij
//  const double young_modulus = 1;
//  const double poisson_ratio = 0.0;
//  MyElasticFEMethod my_fe_mech(mField,Au,Du,Fu,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
//  ElasticFE_RVELagrange_Disp MyFE_RVELagrange_mech(mField,Au,Du,Fu,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank_mech);
//
//  ierr = VecZeroEntries(Fu); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(Fu,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(Fu,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = MatZeroEntries(Au); CHKERRQ(ierr);
//  
//  ierr = mField.loop_finite_elements("PROB_MECH","Kuu",my_fe_mech);  CHKERRQ(ierr);
//  ierr = mField.loop_finite_elements("PROB_MECH","Lagrange_elm_disp",MyFE_RVELagrange_mech);  CHKERRQ(ierr);
//
//  ierr = VecGhostUpdateBegin(Fu,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(Fu,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecAssemblyBegin(Fu); CHKERRQ(ierr);
//  ierr = VecAssemblyEnd(Fu); CHKERRQ(ierr);
//  ierr = MatAssemblyBegin(Au,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//  ierr = MatAssemblyEnd(Au,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//  
//  PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
//  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  //ierr = MatView(Aij,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  
////  //Matrix View
////  MatView(Au,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
////  std::string wait;
////  std::cin >> wait;
//
//  //Solver
//  KSP solverM;
//  ierr = KSPCreate(PETSC_COMM_WORLD,&solverM); CHKERRQ(ierr);
//  ierr = KSPSetOperators(solverM,Au,Au); CHKERRQ(ierr);
//  ierr = KSPSetFromOptions(solverM); CHKERRQ(ierr);
//  ierr = KSPSetUp(solverM); CHKERRQ(ierr);
//  
//  ierr = KSPSolve(solverM,Fu,Du); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(Du,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(Du,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  
//  //Save data on mesh
//  ierr = mField.set_global_VecCreateGhost("PROB_MECH",ROW,Du,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
////  ierr = VecView(Du,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  
//  //=======================================================================================================================================================
//  //Calculation of RVE volume for Homogenized stress
//  //=======================================================================================================================================================
//  double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
//  Vec RVE_volume_Vec;
//  ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
//  ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
//  
//  RVEVolume MyRVEVol(mField,Au,Du,Fu,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio), RVE_volume_Vec);
//  ierr = mField.loop_finite_elements("PROB_MECH","Kuu",MyRVEVol);  CHKERRQ(ierr);
//  //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
//  cout<<"Final RVE_volume = "<< RVE_volume <<endl;
//
//  
//  //=======================================================================================================================================================
//  //Calculation of Homogenized stress mechanical
//  //=======================================================================================================================================================
//  
//  //create a vector for 6 components of homogenized stress
//  Vec Stress_Homo_mech;
//  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_mech);  CHKERRQ(ierr);
//  ierr = VecZeroEntries(Stress_Homo_mech); CHKERRQ(ierr);
//  
////  ierr = VecView(Du,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_disp(mField,Au,Du,Fu,&RVE_volume, applied_strain, Stress_Homo_mech,"DISPLACEMENT","Lagrange_mul_disp",field_rank_mech);
//  ierr = mField.loop_finite_elements("PROB_MECH","Lagrange_elm_disp",MyFE_RVEHomoStressDisp_disp);  CHKERRQ(ierr);
//  
//  //  if(pcomm->rank()) cout<< " Stress_Homo_mech =  "<<endl;
//  //  ierr = VecView(Stress_Homo_mech,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  
//  if(pcomm->rank()==0){
//    PetscScalar    *avec;
//    VecGetArray(Stress_Homo_mech, &avec);
//    
//    cout<< "\nStress_Homo_mech = \n\n";
//    for(int ii=0; ii<6; ii++){
//      cout <<*avec<<endl; ;
//      avec++;
//    }
//  }
//  cout<< "\n\n";
//
//  //=======================================================================================================================================================
//  //Calculation of Homogenized flux moisture
//  //=======================================================================================================================================================
//  
//  //create a vector for 6 components of homogenized stress
//  Vec Stress_Homo_mois;
//  ierr = VecCreateMPI(PETSC_COMM_WORLD, 3, 3*pcomm->size(), &Stress_Homo_mois);  CHKERRQ(ierr);
//  ierr = VecZeroEntries(Stress_Homo_mois); CHKERRQ(ierr);
//  
//  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_mois(mField,Ac,Dc,Fc,&RVE_volume, applied_strain, Stress_Homo_mois,"CONC","Lagrange_mul_conc",field_rank_mois);
//  ierr = mField.loop_finite_elements("PROB_MOIS","Lagrange_elm_conc",MyFE_RVEHomoStressDisp_mois);  CHKERRQ(ierr);
//  
//  //  if(pcomm->rank()) cout<< " Stress_Homo_mois =  "<<endl;
//  //  ierr = VecView(Stress_Homo_mois,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  
//  if(pcomm->rank()==0){
//    PetscScalar    *avec;
//    VecGetArray(Stress_Homo_mois, &avec);
//    
//    cout<< "\nStress_Homo_mois = \n\n";
//    for(int ii=0; ii<6; ii++){
//      cout <<*avec<<endl; ;
//      avec++;
//    }
//  }
//  cout<< "\n\n";
////=======================================================================================================================================================
//
//  
//  PostProcVertexMethod ent_method_mech(moab);
//  ierr = mField.loop_dofs("PROB_MECH","DISPLACEMENT",ROW,ent_method_mech); CHKERRQ(ierr);
//
//  if(pcomm->rank()==0) {
//    EntityHandle out_meshset;
//    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
//    ierr = mField.problem_get_FE("PROB_MECH","Kuu",out_meshset); CHKERRQ(ierr);
//    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
//    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
//  }
//  
//  
//  //Destroy matrices
//  ierr = VecDestroy(&Fu); CHKERRQ(ierr);
//  ierr = VecDestroy(&Du); CHKERRQ(ierr);
//  ierr = MatDestroy(&Au); CHKERRQ(ierr);
//  ierr = VecDestroy(&Fc); CHKERRQ(ierr);
//  ierr = VecDestroy(&Dc); CHKERRQ(ierr);
//  ierr = MatDestroy(&Ac); CHKERRQ(ierr);
//
//  ierr = KSPDestroy(&solverC); CHKERRQ(ierr);
//  ierr = KSPDestroy(&solverM); CHKERRQ(ierr);
//
//  
//  ierr = PetscTime(&v2);CHKERRQ(ierr);
//  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
//  
//  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
//  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  
  PetscFinalize();
  
}

