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

#include <SurfacePressure.hpp>
#include <NodalForce.hpp>
#include <FluidPressure.hpp>
#include <BodyForce.hpp>
#include <ThermalStressElement.hpp>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

#include "Coupled_MechFEMethod.hpp"
#include "MoistureFEMethod.hpp"

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
  int field_rank_mech=3;
  ierr = mField.add_field("DISPLACEMENT",H1,field_rank_mech); CHKERRQ(ierr);
  ierr = mField.add_field("Lagrange_mul_disp",H1,field_rank_mech); CHKERRQ(ierr);
  
  int field_rank_mois=1;
  ierr = mField.add_field("CONC",H1,field_rank_mois); CHKERRQ(ierr);
  ierr = mField.add_field("Lagrange_mul_conc",H1,field_rank_mois); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("Kuu"); CHKERRQ(ierr);
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
  ierr = mField.add_problem("PROB_MECH"); CHKERRQ(ierr);
  ierr = mField.add_problem("PROB_MOIS"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("PROB_MECH","Kuu"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("PROB_MECH","Lagrange_elm_disp"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("PROB_MOIS","Kcc"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("PROB_MOIS","Lagrange_elm_conc"); CHKERRQ(ierr);

  
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("PROB_MECH",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("PROB_MOIS",bit_level0); CHKERRQ(ierr);

  /***/
  //Declare problem
  
  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT",2); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"CONC",2); CHKERRQ(ierr);

  
  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"Kuu",MBTET); CHKERRQ(ierr);
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
  
  //partition PROB_MECH
  ierr = mField.partition_problem("PROB_MECH"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("PROB_MECH"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("PROB_MECH"); CHKERRQ(ierr);
  
  //partition PROB_MOIS
  ierr = mField.partition_problem("PROB_MOIS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("PROB_MOIS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("PROB_MOIS"); CHKERRQ(ierr);

  //print block sets with materials
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);
  
  //create matrices and vectors
  Vec Fu,Du,Fc,Dc;
  ierr = mField.VecCreateGhost("PROB_MECH",ROW,&Fu); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("PROB_MECH",COL,&Du); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("PROB_MOIS",ROW,&Fc); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("PROB_MOIS",COL,&Dc); CHKERRQ(ierr);

  Mat Au,Ac;
  ierr = mField.MatCreateMPIAIJWithArrays("PROB_MECH",&Au); CHKERRQ(ierr);
  ierr = mField.MatCreateMPIAIJWithArrays("PROB_MOIS",&Ac); CHKERRQ(ierr);

  //***************************************************************************************************
  //First solve moisture problem to calculate moisture concentration
  //***************************************************************************************************
  
  MoistureFEMethod my_fe_mois(mField,Ac,Dc,Fc);
  ElasticFE_RVELagrange_Disp MyFE_RVELagrange_mois(mField,Ac,Dc,Fc,applied_congrad,"CONC","Lagrange_mul_conc",field_rank_mois);
  ierr = VecZeroEntries(Fc); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(Fc,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(Fc,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Ac); CHKERRQ(ierr);

  ierr = mField.loop_finite_elements("PROB_MOIS","Kcc",my_fe_mois);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("PROB_MOIS","Lagrange_elm_conc",MyFE_RVELagrange_mois);  CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(Fc,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(Fc,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(Fc); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Fc); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Ac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Ac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //ierr = MatView(Aij,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

//  //Matrix View
//  MatView(Ac,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//  std::string wait;
//  std::cin >> wait;

  //Solver
  KSP solverC;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solverC); CHKERRQ(ierr);
  ierr = KSPSetOperators(solverC,Ac,Ac); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solverC); CHKERRQ(ierr);
  ierr = KSPSetUp(solverC); CHKERRQ(ierr);

  ierr = KSPSolve(solverC,Fc,Dc); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(Dc,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(Dc,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("PROB_MOIS",ROW,Dc,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(Dc,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  

  struct MyElasticFEMethod: public Coupled_MechFEMethod {
    MyElasticFEMethod(FieldInterface& _mField, Mat &_Aij,Vec &_D,Vec& _F,double _lambda,double _mu):
    Coupled_MechFEMethod(_mField,_Aij,_D,_F,_lambda,_mu) {};
    
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

  //Assemble F and Aij
  const double young_modulus = 1;
  const double poisson_ratio = 0.0;
  MyElasticFEMethod my_fe_mech(mField,Au,Du,Fu,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ElasticFE_RVELagrange_Disp MyFE_RVELagrange_mech(mField,Au,Du,Fu,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank_mech);

  ierr = VecZeroEntries(Fu); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(Fu,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(Fu,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Au); CHKERRQ(ierr);
  
  ierr = mField.loop_finite_elements("PROB_MECH","Kuu",my_fe_mech);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("PROB_MECH","Lagrange_elm_disp",MyFE_RVELagrange_mech);  CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(Fu,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(Fu,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(Fu); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Fu); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Au,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Au,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //ierr = MatView(Aij,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
//  //Matrix View
//  MatView(Au,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//  std::string wait;
//  std::cin >> wait;

  //Solver
  KSP solverM;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solverM); CHKERRQ(ierr);
  ierr = KSPSetOperators(solverM,Au,Au); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solverM); CHKERRQ(ierr);
  ierr = KSPSetUp(solverM); CHKERRQ(ierr);
  
  ierr = KSPSolve(solverM,Fu,Du); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(Du,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(Du,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("PROB_MECH",ROW,Du,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(Du,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  //=======================================================================================================================================================
  //Calculation of RVE volume for Homogenized stress
  //=======================================================================================================================================================
  double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
  Vec RVE_volume_Vec;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
  ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
  
  RVEVolume MyRVEVol(mField,Au,Du,Fu,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio), RVE_volume_Vec);
  ierr = mField.loop_finite_elements("PROB_MECH","Kuu",MyRVEVol);  CHKERRQ(ierr);
  //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
  cout<<"Final RVE_volume = "<< RVE_volume <<endl;

  
  //=======================================================================================================================================================
  //Calculation of Homogenized stress mechanical
  //=======================================================================================================================================================
  
  //create a vector for 6 components of homogenized stress
  Vec Stress_Homo_mech;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_mech);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo_mech); CHKERRQ(ierr);
  
//  ierr = VecView(Du,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_disp(mField,Au,Du,Fu,&RVE_volume, applied_strain, Stress_Homo_mech,"DISPLACEMENT","Lagrange_mul_disp",field_rank_mech);
  ierr = mField.loop_finite_elements("PROB_MECH","Lagrange_elm_disp",MyFE_RVEHomoStressDisp_disp);  CHKERRQ(ierr);
  
  //  if(pcomm->rank()) cout<< " Stress_Homo_mech =  "<<endl;
  //  ierr = VecView(Stress_Homo_mech,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec;
    VecGetArray(Stress_Homo_mech, &avec);
    
    cout<< "\nStress_Homo_mech = \n\n";
    for(int ii=0; ii<6; ii++){
      cout <<*avec<<endl; ;
      avec++;
    }
  }
  cout<< "\n\n";

  //=======================================================================================================================================================
  //Calculation of Homogenized flux moisture
  //=======================================================================================================================================================
  
  //create a vector for 6 components of homogenized stress
  Vec Stress_Homo_mois;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 3, 3*pcomm->size(), &Stress_Homo_mois);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo_mois); CHKERRQ(ierr);
  
  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_mois(mField,Ac,Dc,Fc,&RVE_volume, applied_strain, Stress_Homo_mois,"CONC","Lagrange_mul_conc",field_rank_mois);
  ierr = mField.loop_finite_elements("PROB_MOIS","Lagrange_elm_conc",MyFE_RVEHomoStressDisp_mois);  CHKERRQ(ierr);
  
  //  if(pcomm->rank()) cout<< " Stress_Homo_mois =  "<<endl;
  //  ierr = VecView(Stress_Homo_mois,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec;
    VecGetArray(Stress_Homo_mois, &avec);
    
    cout<< "\nStress_Homo_mois = \n\n";
    for(int ii=0; ii<6; ii++){
      cout <<*avec<<endl; ;
      avec++;
    }
  }
  cout<< "\n\n";
//=======================================================================================================================================================

  
  PostProcVertexMethod ent_method_mech(moab);
  ierr = mField.loop_dofs("PROB_MECH","DISPLACEMENT",ROW,ent_method_mech); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("PROB_MECH","Kuu",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  
  
  //Destroy matrices
  ierr = VecDestroy(&Fu); CHKERRQ(ierr);
  ierr = VecDestroy(&Du); CHKERRQ(ierr);
  ierr = MatDestroy(&Au); CHKERRQ(ierr);
  ierr = VecDestroy(&Fc); CHKERRQ(ierr);
  ierr = VecDestroy(&Dc); CHKERRQ(ierr);
  ierr = MatDestroy(&Ac); CHKERRQ(ierr);

  ierr = KSPDestroy(&solverC); CHKERRQ(ierr);
  ierr = KSPDestroy(&solverM); CHKERRQ(ierr);

  
  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  
  PetscFinalize();
  
}

