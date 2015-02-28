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
#include <BodyForce.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";

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
  BARRIER_RANK_START(pcomm) 
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  BARRIER_RANK_END(pcomm) 

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& mField = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);


  //FE
  ierr = mField.add_finite_element("TEST_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("TEST_FE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("TEST_FE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TEST_FE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TEST_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //Problem
  ierr = mField.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);


  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = mField.add_ents_to_finite_element_by_TETs(root_set,"TEST_FE"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 2;
  ierr = mField.set_field_order(root_set,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(root_set,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

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
  ierr = mField.simple_partition_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  //set from positions of 10 node tets
  Projection10NodeCoordsOnField ent_method(mField,"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method); CHKERRQ(ierr);

  Vec F;
  ierr = mField.VecCreateGhost("TEST_PROBLEM",ROW,&F); CHKERRQ(ierr);

  typedef tee_device<ostream, ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;
  ofstream ofs("body_force_atom_test.txt");
  TeeDevice my_tee(cout, ofs); 
  TeeStream my_split(my_tee);

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  BodyFroceConstantField body_forces_methods(mField);

  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|BODYFORCESSET,it)) {
    Block_BodyForces mydata;
    ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
    my_split << mydata << endl;
    ierr = body_forces_methods.addBlock("DISPLACEMENT",F,it->get_msId()); CHKERRQ(ierr);
  }
  ierr = mField.loop_finite_elements("TEST_PROBLEM","TEST_FE",body_forces_methods.getLoopFe()); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  ierr = mField.set_global_VecCreateGhost("TEST_PROBLEM",COL,F,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  const MoFEMProblem *problemPtr;
  ierr = mField.get_problem("TEST_PROBLEM",&problemPtr); CHKERRQ(ierr);
  map<EntityHandle,double> m0,m1,m2;
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(problemPtr,dit)) {

    if(dit->get_dof_rank()!=1) continue;

    my_split.precision(3);
    my_split.setf(std::ios::fixed);
    my_split << dit->get_petsc_gloabl_dof_idx() << " " << dit->get_FieldData() << endl;

  }

  double sum = 0;
  ierr = VecSum(F,&sum); CHKERRQ(ierr);
  my_split << endl << "Sum : " << setprecision(3) << sum << endl;

  ierr = VecDestroy(&F); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


