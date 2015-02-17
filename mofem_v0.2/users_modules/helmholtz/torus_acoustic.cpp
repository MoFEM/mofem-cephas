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

//#include "FieldInterface.hpp"
//#include "FieldCore.hpp"

//#include "ForcesAndSurcesCore.hpp"
//#include "TsCtx.hpp"
#include <MoFEM.hpp>


#include <DirichletBC.hpp>
#include <PotsProcOnRefMesh.hpp>
#include <HelmholtzElement.hpp>
//#include "DirichletBC.hpp"
//#include "FEM.h"
//#include "FEMethod_UpLevelStudent.hpp"
//#include "PostProcVertexMethod.hpp"
//#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include <Projection10NodeCoordsOnField.hpp>
#include <AnalyticalDirichletacoustic.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>
#include <cmath>
#include <boost/math/special_functions.hpp>
#include <complex>
//#include <petscksp.h>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";
//argc = argument counts, argv = argument vectors
int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  //Core mb_instance;
  moab::Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  BARRIER_RANK_START(pcomm) 
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  BARRIER_RANK_END(pcomm) 

  //Create MoFEM (Joseph) database
  //FieldCore core(moab);
  MoFEM::Core core(moab);
  FieldInterface& mField = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_field("rePRES",H1,1); CHKERRQ(ierr);  //field order distinguish the scalar field and vector field.
  ierr = mField.add_field("imPRES",H1,1); CHKERRQ(ierr);
  
  //Problem
  ierr = mField.add_problem("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.add_problem("BC_PROBLEM"); CHKERRQ(ierr); //analytical Dirichlet
  
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ACOUSTIC_PROBLEM",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("BC_PROBLEM",bit_level0); CHKERRQ(ierr);  //analytical Dirichlet

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"rePRES"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"imPRES"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 2;
  }

  ierr = mField.set_field_order(root_set,MBTET,"rePRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"rePRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"rePRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"rePRES",1); CHKERRQ(ierr);

  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(root_set,MBTET,"imPRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"imPRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"imPRES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"imPRES",1); CHKERRQ(ierr);

  HelmholtzElement helmholtz_elements(mField);               //Create the HelmholtzElement class in the header-file
  ierr = helmholtz_elements.addHelmholtzElements("ACOUSTIC_PROBLEM","rePRES","imPRES"); CHKERRQ(ierr);
  ierr = helmholtz_elements.addHelmholtzFluxElement("ACOUSTIC_PROBLEM","rePRES","imPRES"); CHKERRQ(ierr);
  ierr = helmholtz_elements.addHelmholtzImpedanceElement("ACOUSTIC_PROBLEM","rePRES","imPRES"); CHKERRQ(ierr);
  
  
  //Set up the analytical Dirichlet BC.
  Range bc_tris;
  for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"ANALYTICAL_BC",it)) {
   rval = moab.get_entities_by_type(it->get_meshset(),MBTRI,bc_tris,true); CHKERR_PETSC(rval);
  }
  
  AnalyticalDirihletBC analytical_bc1(mField,bc_tris);
  AnalyticalDirihletBC analytical_bc2(mField,bc_tris);
  ierr = analytical_bc1.initializeProblem(mField,"BC_PROBLEM","BC_FE","rePRES"); CHKERRQ(ierr);
  ierr = analytical_bc2.initializeProblem(mField,"BC_PROBLEM","BC_FE","imPRES"); CHKERRQ(ierr);
  
  

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

  Projection10NodeCoordsOnField ent_method_material(mField,"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  /****/
  //mesh partitioning 
  //partition
  ierr = mField.partition_problem("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
  
  //mesh partitinoning for analytical Dirichlet
  ierr = mField.simple_partition_problem("BC_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("BC_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("BC_PROBLEM"); CHKERRQ(ierr);
  
  Vec F;  //Right hand side vector
  ierr = mField.VecCreateGhost("ACOUSTIC_PROBLEM",ROW,&F); CHKERRQ(ierr);
  Vec T; //Solution vector
  ierr = VecDuplicate(F,&T); CHKERRQ(ierr);
  Mat A; //Left hand side matrix
  ierr = mField.MatCreateMPIAIJWithArrays("ACOUSTIC_PROBLEM",&A); CHKERRQ(ierr);
  
  //Create the Imaginary part of the left hand side operator.
  //Mat C;
  //ierr = mField.MatCreateMPIAIJWithArrays("ACOUSTIC_PROBLEM",&C); CHKERRQ(ierr);

  int ii1,jj1;
  ierr=MatGetSize(A,&ii1,&jj1);
  std::string wait;
  //MatView(A,PETSC_VIEWER_DRAW_WORLD);
  std::cin >> wait;
  //PetscInt i,j,N1,N2,N3,N4;

  //TemperatureBCFEMethodPreAndPostProc my_dirichlet_bc1(mField,"rePRES",A,T,F);//scalar Dirichlet
  //TemperatureBCFEMethodPreAndPostProc my_dirichlet_bc2(mField,"imPRES",A,T,F);//scalar Dirichlet
  //TemperatureBCFEMethodPreAndPostProc my_dirichlet_bc(mField,"TEMP",C,T,F);
  //ierr=MatCopy(A,C,DIFFERENT_NONZERO_PATTERN); //Copy the sparse matrix A to C with identical size.
  
  bool useScalar = true;
  //ierr = helmholtz_elements.setHelmholtzFiniteElementRhsOperators("rePRES","rePRES",F,useScalar); CHKERRQ(ierr); //The analytical F source vector
  ierr = helmholtz_elements.setHelmholtzFiniteElementLhsOperators("rePRES","imPRES",(A)); CHKERRQ(ierr);
  ierr = helmholtz_elements.setHelmholtzMassFiniteElementLhsOperators("rePRES","imPRES",(A)); CHKERRQ(ierr);
  ierr = helmholtz_elements.setHelmholtzFluxFiniteElementRhsOperators("rePRES","rePRES",F); CHKERRQ(ierr);  //real Neumann BC
  //ierr = helmholtz_elements.setHelmholtzFluxFiniteElementRhsOperators("imPRES","imPRES",F); CHKERRQ(ierr);    //Imag Neumann BC
  //ierr = helmholtz_elements.setHelmholtzIncidentWaveFiniteElementRhsOperators("imPRES","imPRES",F); CHKERRQ(ierr); // Incident wave flux.
  
  ierr = helmholtz_elements.setHelmholtzImpedanceFiniteElementLhsOperators("rePRES","imPRES",(A)); CHKERRQ(ierr);
  
  
  
  ierr = VecZeroEntries(T); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);
  
  //analytical dirihlet bc
  AnalyticalDirihletBC::DirihletBC analytical_ditihlet_bc1(mField,"rePRES",A,T,F);
  AnalyticalDirihletBC::DirihletBC analytical_ditihlet_bs2(mField,"imPRES",A,T,F);

  //solve for analytical ditihlet bc dofs
  ierr = analytical_bc1.setProblem(mField,"BC_PROBLEM"); CHKERRQ(ierr);
  ierr = analytical_bc2.setProblem(mField,"BC_PROBLEM"); CHKERRQ(ierr);

  //MoFEM::HelmholtzElement::BlockData blockData;
  //double wavenumber = helmholtz_elements.getwavenumber();
  //// = helmholtz_elements.getblockdata();
  ////MoFEM::HelmholtzElement::getblockdata();
  
  static double aNgularfreq;
  static double sPeed; //Without static. I got error:use of local variable with automatic storage from containing function
  //for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_HELMHOLTZSET,it)) {
//	  Mat_Helmholtz pres_data;
//	  ierr = it->get_attribute_data_structure(pres_data); CHKERRQ(ierr);
//	  //wait
//	  aNgularfreq = pres_data.data.AngularFreq;
//	  sPeed = pres_data.data.Speed;
  //}
  
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it))
  {
	  cout << endl << *it << endl;
	  
	  //Get block name
	  string name = it->get_Cubit_name();
	  if (name.compare(0,13,"MAT_HELMHOLTZ") == 0)
	  {
		  Mat_Helmholtz pres_data;
		  ierr = it->get_attribute_data_structure(pres_data); CHKERRQ(ierr);
		  cout << pres_data;
		  aNgularfreq = pres_data.data.AngularFreq;
		  sPeed = pres_data.data.Speed;
	  }
  }
  
  
  struct AnaliticalFunction {
	  static double fUN_real(double x,double y,double z) {
		  //return pow(x,1);
		  
		  const double wAvenumber = aNgularfreq/sPeed;
		  
		  const double k = wAvenumber;  //Wave number
		  const double pi = atan( 1.0 ) * 4.0;  //PI
		  const1 = (k*x)/sqrt(3);
		  const2 = (k*y)/sqrt(3);
		  const3 = (k*z)/sqrt(3);
		  result = sin(const1)*sin(const2)*sin(const3);
		  
          return std::real(result);
      
	  }
	  static double fUN_imag(double x,double y,double z) {
		  const double wAvenumber = aNgularfreq/sPeed;
		  
		  const double k = wAvenumber;  //Wave number
		  const double pi = atan( 1.0 ) * 4.0;  //PI
		  const1 = (k*x)/sqrt(3);
		  const2 = (k*y)/sqrt(3);
		  const3 = (k*z)/sqrt(3);
		  result = sin(const1)*sin(const2)*sin(const3);
		  
          return std::imag(result);
	  
	  }
  };
  
  //solve the analytical Dirichlet BC, results in a solution vector contains the the solution on Dirichlet nodes only. 
  ierr = analytical_bc1.setApproxOps(mField,"rePRES",AnaliticalFunction::fUN_real); CHKERRQ(ierr);
  ierr = analytical_bc2.setApproxOps(mField,"imPRES",AnaliticalFunction::fUN_imag); CHKERRQ(ierr);

  ierr = analytical_bc1.solveProblem(mField,"BC_PROBLEM","BC_FE",analytical_ditihlet_bc1); CHKERRQ(ierr);
  ierr = analytical_bc2.solveProblem(mField,"BC_PROBLEM","BC_FE",analytical_ditihlet_bc2); CHKERRQ(ierr);

  ierr = analytical_bc1.destroyProblem(); CHKERRQ(ierr);
  ierr = analytical_bc2.destroyProblem(); CHKERRQ(ierr);
  
  
  //preproc
  //Preprocess the analytical Dirichlet BC
  ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc1); CHKERRQ(ierr);
  ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc2); CHKERRQ(ierr);

  ierr = mField.set_global_VecCreateGhost("ACOUSTIC_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  
  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",helmholtz_elements.getLoopFeRhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",helmholtz_elements.getLoopFeLhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FLUX_FE",helmholtz_elements.getLoopFeFlux()); CHKERRQ(ierr); //scalar flux
  //ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FLUX_FE",helmholtz_elements.getLoopfeIncidentWave()); CHKERRQ(ierr); //Incident wave flux
  
  /*Wait for confirmation */
  //ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_IMPEDANCE_FE",helmholtz_elements.getLoopFeImpedanceRhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_IMPEDANCE_FE",helmholtz_elements.getLoopFeImpedanceLhs()); CHKERRQ(ierr);
  /*above terms related to operators in HelmholtzElement.hpp */
  
  //postproc
  //Postprocess the Analytical Dirichlet BC
  ierr = mField.problem_basic_method_postProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc1); CHKERRQ(ierr);
  ierr = mField.problem_basic_method_postProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc2); CHKERRQ(ierr);



  
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecScale(F,-1); CHKERRQ(ierr);
  
  //MatView(A,PETSC_VIEWER_DRAW_WORLD);
  std::cin >> wait;
  
  //ierr=VecGetSize(F,&N3);
  //ierr=VecGetSize(T,&N4);
  //std::cout<< "\n Solution Vector T = \n" << N3 << std::endl;
  //ierr = VecView(T,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //std::cin >> wait;
  
  //std::cout<< "\n RHS Vector F = \n" << N4 << std::endl;
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //std::cin >> wait;
  

  ierr=VecGetSize(F,&N1);
  ierr=VecGetSize(T,&N2);
  //std::cout<< "\n Solution Vector T = \n" << N2 << std::endl;
  //ierr = VecView(T,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //std::cin >> wait;
  
  //std::cout<< "\n RHS Vector F = \n" << N1 << std::endl;
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //std::cin >> wait;
  
  ////the way of implement the CR-Transformation is not elegant. Wait for fix
  ////Right hand side vector
  //PetscInt i,j,N1,N2,N3,N4;
  //ierr=VecGetSize(F,&N3);
  //Vec F_CR,T_CR;
  //ierr = VecCreate(PETSC_COMM_WORLD,&F_CR);
  //PetscInt vecsize = 2*N3;
  //ierr = VecSetSizes(F_CR,PETSC_DECIDE,vecsize);
  //ierr = VecSetFromOptions(F_CR);
  //ierr = VecZeroEntries(F_CR); CHKERRQ(ierr);
  //ierr = VecGhostUpdateBegin(F_CR,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //ierr = VecGhostUpdateEnd(F_CR,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ////ierr=VecGetSize(F_CR,&N4);
  
  //PetscInt loc [N3];
  //for (i=0; i<N3; ++ i) {
//	  loc[i] = i;
  //}
  //PetscScalar f_cr [N4];
  //ierr = VecGetValues(F,N3,loc,f_cr); //Currently can only get values on the same processor, this routine will be a problem
  ////ierr = VecNestSetSubVecs(F_CR,N3,loc,&F);  //when using multiple processors. How to do it?
  //VecSetValues(F_CR,N3,loc,f_cr,ADD_VALUES);   //Embed the real part of force vector into the complex force vector
  //ierr = VecGhostUpdateBegin(F_CR,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecGhostUpdateEnd(F_CR,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecAssemblyBegin(F_CR); CHKERRQ(ierr);
  //ierr = VecAssemblyEnd(F_CR); CHKERRQ(ierr);
  
  
  
  //the way of implement the CR-Transformation is not elegant. Wait for fix
  //Now We deal with the Matrix embed
  //PetscInt i1,j1,i2,j2,i3,j3;
  //ierr=MatGetSize(A,&i1,&j1);
  //ierr=MatGetSize(C,&i2,&j2);
  //i3 = i1*2;
  //j3 = j1*2;
  //PetscInt row1 [i1];
  //for (j=0; j<i1; ++ j) {
//	  row1[j] = j;
  //}
  //PetscInt row2 [i1];
  //for (j=i1; j < (i1*2); ++ j) {
//	  row2[j-i1] = j;
  //}
  //PetscInt col1 [j1];
  //for (j=0; j<j1; ++ j) {
//	  col1[j] = j;
  //}
  //PetscInt col2 [j1];
  //for (j=j1; j < (j1*2); ++ j) {
//	  col2[j-j1] = j;
  //}
  ////ierr = PetscIntView(j1,col1,PETSC_VIEWER_STDOUT_WORLD);
  ////  ierr = PetscIntView(j1,col2,PETSC_VIEWER_STDOUT_WORLD);
  //PetscScalar A_array [i1*j1];
  //ierr = MatGetValues(A,i1,row1,j1,col1,A_array);
  ////ierr = PetscScalarView(i1*j1,A_array,PETSC_VIEWER_STDOUT_WORLD);
  //PetscScalar C_array [i2*j2];
  //PetscScalar C_array_scale [i2*j2];
  //ierr = MatGetValues(C,i2,row1,j2,col1,C_array);
  //ierr = MatScale(C,-1);
  //ierr = MatGetValues(C,i2,row1,j2,col1,C_array_scale);
  ////ierr = MatView(C,PETSC_VIEWER_DRAW_WORLD);
  ////std::cin >> wait;
  ////ierr = PetscScalarView(i2*j2,C_array_scale,PETSC_VIEWER_STDOUT_WORLD);
  //Mat A_CR;
  //ierr = MatCreate(PETSC_COMM_WORLD,&A_CR);
  //ierr = MatSetSizes(A_CR,PETSC_DECIDE,PETSC_DECIDE,i1+i2,j1+j2);
  //ierr = MatSetType(A_CR,MATMPIAIJ);
  ////////ierr = MatSetBlockSize(Pmat,2);
  ////////ierr = MatSetType(Pmat,MATAIJ);
  ////////ierr = MatSeqAIJSetPreallocation(Pmat,18,NULL);
  ////////ierr = MatMPIAIJSetPreallocation(Pmat,18,NULL,12,NULL);
  //ierr = MatSetUp(A_CR);
  //ierr = MatZeroEntries(A_CR); CHKERRQ(ierr);
  //ierr = MatSetValues(A_CR,i1,row1,j1,col1,A_array,ADD_VALUES);
  //ierr = MatSetValues(A_CR,i1,row2,j1,col2,A_array,ADD_VALUES);
  //ierr = MatSetValues(A_CR,i2,row1,j2,col2,C_array_scale,ADD_VALUES);
  //ierr = MatSetValues(A_CR,i2,row2,j2,col1,C_array,ADD_VALUES);
  //ierr = MatAssemblyBegin(A_CR,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  //ierr = MatAssemblyEnd(A_CR,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ////ierr = MatView(A_CR,PETSC_VIEWER_DRAW_WORLD);
  ////std::cin >> wait;
  //The way of manually implement the A_CR Matrix
  
  

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,T); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  //there the vector T^tude and F^tude is without the Dirichlet BC value inserted.
  //ierr = VecView(T,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //below fill in the solution vector T with analytical dirichlet solutions.
  ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc1); CHKERRQ(ierr);
  ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analytical_ditihlet_bc2); CHKERRQ(ierr);


  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ACOUSTIC_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  
  //ierr = VecView(T,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  

  
  //Wait to putput the data in format
  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Acoustic_Impining_Sphere.txt",&viewer); CHKERRQ(ierr);
  VecView(T,viewer);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  
  //wait to calculate the the magnitude of acoustic potential-abs(phi) using real and imag value.

  //Range ref_edges;
  //ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBEDGE,ref_edges); CHKERRQ(ierr);
  //rval = moab.list_entities(ref_edges); CHKERR_PETSC(rval);
  //mField.list_dofs_by_field_name("TEMP");

  if(pcomm->rank()==0) {
    rval = moab.write_file("solution.h5m"); CHKERR_PETSC(rval);
  }

  /*EntityHandle fe_meshset = mField.get_finite_element_meshset("HELMHOLTZ_FE");
  Range tets;
  rval = moab.get_entities_by_type(fe_meshset,MBTET,tets,true); CHKERR_PETSC(rval);
  Range tets_edges;
  rval = moab.get_adjacencies(tets,1,false,tets_edges,Interface::UNION); CHKERR(rval);
  EntityHandle edges_meshset;
  rval = moab.create_meshset(MESHSET_SET,edges_meshset); CHKERR_PETSC(rval);
  rval = moab.add_entities(edges_meshset,tets); CHKERR_PETSC(rval);
  rval = moab.add_entities(edges_meshset,tets_edges); CHKERR_PETSC(rval);
  rval = moab.convert_entities(edges_meshset,true,false,false); CHKERR_PETSC(rval);*/

  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet1(mField,"rePRES",true,false,"rePRES");
  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet2(mField,"imPRES",true,false,"imPRES");
  ierr = mField.loop_dofs("rePRES",ent_method_on_10nodeTet1); CHKERRQ(ierr);
  ierr = mField.loop_dofs("imPRES",ent_method_on_10nodeTet2); CHKERRQ(ierr);
  ent_method_on_10nodeTet1.set_nodes = false;
  ent_method_on_10nodeTet2.set_nodes = false;
  ierr = mField.loop_dofs("rePRES",ent_method_on_10nodeTet1); CHKERRQ(ierr);
  ierr = mField.loop_dofs("imPRES",ent_method_on_10nodeTet2); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }


  //Copy the output .vtk file to desired location.
  char command[] = "cp out.vtk ../../../../../mnt/home/Documents/mofem-cephas/mofem_v0.2/users_modules/helmholtz/";
  
  int status = system( command );
  
  //if(pcomm->rank()==0) {
  //  PostProcScalarFieldsAndGradientOnRefMesh fe_post_proc_method1(moab,"rePRES");
//	PostProcScalarFieldsAndGradientOnRefMesh fe_post_proc_method2(moab,"imPRES");
  //  fe_post_proc_method1.do_broadcast = false;
//	fe_post_proc_method2.do_broadcast = false;
  //  ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",fe_post_proc_method1,0,pcomm->size());  CHKERRQ(ierr);
//	ierr = mField.loop_finite_elements("ACOUSTIC_PROBLEM","HELMHOLTZ_FE",fe_post_proc_method2,0,pcomm->size());  CHKERRQ(ierr);
  //  rval = fe_post_proc_method1.moab_post_proc.write_file("out_post_proc_re.vtk","VTK",""); CHKERR_PETSC(rval);
//	rval = fe_post_proc_method2.moab_post_proc.write_file("out_post_proc_im.vtk","VTK",""); CHKERR_PETSC(rval); //wait, it seems does not give
  //}                 																					//give the right value of Im(Scalar Solution)

  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&T); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


