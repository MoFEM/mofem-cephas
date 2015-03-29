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

#include <MoFEM.hpp>
#include <Projection10NodeCoordsOnField.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <FiledApproximation.hpp>
#include <PotsProcOnRefMesh.hpp>

#define HOON

using namespace MoFEM;

static char help[] = "...\n\n";

struct MyFunApprox {

  ublas::vector<double> result;
  ublas::vector<double>& operator()(double x, double y, double z) {
    result.resize(3);
    result[0] = x;
    result[1] = y;
    const int k = 4;
    result[2] = sin(M_PI*z*k);
    return result;
  }     

};

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

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
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("FIELD1",H1,3); CHKERRQ(ierr);
  #ifdef HOON
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  #endif

  //FE
  ierr = m_field.add_finite_element("TEST_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE","FIELD1"); CHKERRQ(ierr);
  #ifdef HOON
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  #endif

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"FIELD1"); CHKERRQ(ierr);
  #ifdef HOON
    ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  #endif
  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_TETs(root_set,"TEST_FE"); CHKERRQ(ierr);


  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 3;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 3;
  }
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD1",1); CHKERRQ(ierr);
  #ifdef HOON
  ierr = m_field.set_field_order(root_set,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  #endif

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  #ifdef HOON
  Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
  #endif
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning 
  //partition
  ierr = m_field.simple_partition_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  Mat A;
  ierr = m_field.MatCreateMPIAIJWithArrays("TEST_PROBLEM",&A); CHKERRQ(ierr);
  Vec D,F;
  ierr = m_field.VecCreateGhost("TEST_PROBLEM",ROW,&F); CHKERRQ(ierr);
  ierr = m_field.VecCreateGhost("TEST_PROBLEM",COL,&D); CHKERRQ(ierr);

  static double aNgularfreq;
  static double sPeed; //Without static. I got error:use of local variable with automatic storage from containing function
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it))
  {
	  //Get block name
	  string name = it->get_Cubit_name();
	  if (name.compare(0,12,"MAT_HELMHOLTZ") == 0)
	  {
		  //get block attributes
		  vector<double> attributes;
		  ierr = it->get_Cubit_attributes(attributes); CHKERRQ(ierr);
		  if(attributes.size()<2) {
			  SETERRQ1(PETSC_COMM_SELF,1,"not enough block attributes to deffine fluid pressure element, attributes.size() = %d ",attributes.size());
		  }
		  aNgularfreq = attributes[0];
		  sPeed = attributes[1];
		  std::string wait;
		  std::cout << "\n sPeed = \n" << sPeed << "\n aNgularfreq = " << aNgularfreq << std::endl;
		  
	  }
  }
  
  
  {
    MyFunApprox function_evaluator;
    FieldApproximationH1<MyFunApprox> field_approximation(m_field);
    field_approximation.loopMatrixAndVector(
      "TEST_PROBLEM","TEST_FE","FIELD1",A,F,function_evaluator);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = m_field.set_global_ghost_vector("TEST_PROBLEM",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);

  PostPocOnRefinedMesh post_proc(m_field);
  ierr = post_proc.generateRefereneElemenMesh(); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("FIELD1"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesGradientPostProc("FIELD1"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE",post_proc); CHKERRQ(ierr);
  rval = post_proc.postProcMesh.write_file("out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


