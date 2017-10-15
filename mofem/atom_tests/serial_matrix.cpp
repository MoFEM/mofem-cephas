/**
 * \file serial_matrix.cpp
 *
 * \brief testing create serial/sequential matrix
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

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,PETSC_NULL,help);

  try {

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;

  //Reade parameters from line command
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetString(PETSC_NULL,"","-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #else
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #endif
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  PetscInt order;
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetInt(PETSC_NULL,"","-my_order",&order,&flg); CHKERRQ(ierr);
  #else
  ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  #endif
  if(flg != PETSC_TRUE) {
    order = 1;
  }

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);


  int do_for_rank = 0;
  if(pcomm->rank()==(unsigned int)do_for_rank) { // should work only with rank 0

  //Create MoFEM (Joseph) database
  //second argument set communicator for sequential problem
  //last argument make mofem quaiet
  MoFEM::Core core(moab,PETSC_COMM_SELF,-1);
  MoFEM::Interface& m_field = core;

  //ref meshset ref level 0
  ierr = m_field.query_interface<BitRefManager>()->setBitRefLevelByDim(0,3,0); CHKERRQ(ierr);

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.query_interface<BitRefManager>()->setBitRefLevelByDim(0,3,bit_level0); CHKERRQ(ierr);
  ierr = m_field.query_interface<Tools>()->getEntitiesByRefLevel(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  /***/
  //Define problem

  //Fields
  ierr = m_field.add_field("FIELD_A",H1,AINSWORTH_LEGENDRE_BASE,3,MB_TAG_DENSE); CHKERRQ(ierr);
  ierr = m_field.add_field("FIELD_B",L2,AINSWORTH_LEGENDRE_BASE,1,MB_TAG_DENSE); CHKERRQ(ierr);

  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"FIELD_A"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"FIELD_B"); CHKERRQ(ierr);

  ierr = m_field.set_field_order(0,MBTET,"FIELD_A",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"FIELD_A",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"FIELD_A",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"FIELD_A",1); CHKERRQ(ierr);

  ierr = m_field.set_field_order(0,MBTET,"FIELD_B",order-1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"FIELD_B",order-1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"FIELD_B",order-1); CHKERRQ(ierr);

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);

  //Element
  ierr = m_field.add_finite_element("FE1"); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("FE2"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("FE1","FIELD_A"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("FE1","FIELD_A"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("FE1","FIELD_A"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("FE1","FIELD_B"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("FE1","FIELD_B"); CHKERRQ(ierr);


  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("FE2","FIELD_B"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("FE2","FIELD_A"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("FE2","FIELD_B"); CHKERRQ(ierr);

  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_type(0,MBTET,"FE1"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_type(0,MBTET,"FE2"); CHKERRQ(ierr);

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);


  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","FE1"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","FE2"); CHKERRQ(ierr);

  ProblemsManager *prb_mng_ptr;
  ierr = m_field.query_interface(prb_mng_ptr); CHKERRQ(ierr);
  ierr = prb_mng_ptr->buildProblem("TEST_PROBLEM",true); CHKERRQ(ierr);
  ierr = prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM"); CHKERRQ(ierr);

  // partition finite elements
  ierr = prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM"); CHKERRQ(ierr);
  // what are ghost nodes, see Petsc Manual
  ierr = prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM"); CHKERRQ(ierr);

  Vec F;
  ierr = m_field.query_interface<VecManager>()->vecCreateSeq("TEST_PROBLEM",ROW,&F); CHKERRQ(ierr);
  Mat A;
  ierr = m_field.MatCreateMPIAIJWithArrays("TEST_PROBLEM",&A); CHKERRQ(ierr);

  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"serial_matrix.txt",&viewer); CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO); CHKERRQ(ierr);
  ierr = MatView(A,viewer); CHKERRQ(ierr);

  //ierr = MatView(A,PETSC_VIEWER_DRAW_SELF); CHKERRQ(ierr);
  //std::string wait;
  //std::cin >> wait;

  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);

  }

  // if(pcomm->rank()!=(unsigned int)do_for_rank) {
  //   std::string wait;
  //   std::cin >> wait;
  // }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFinalize();
  return 0;

}
