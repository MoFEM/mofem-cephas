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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "ForcesAndSurcesCore.hpp"
#include "FiledApproximation.hpp"

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include "Projection10NodeCoordsOnField.hpp"

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";


struct MyFunApprox {

  ublas::vector<double> result;
  ublas::vector<double>& operator()(double x, double y, double z) {
    result.resize(3);
    result[0] = x;
    result[1] = y;
    result[2] = z*z;
    return result;
  }     

};

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  BARRIER_RANK_START(pcomm) 
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  BARRIER_RANK_END(pcomm) 

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_field("FIELD1",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("TEST_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TEST_FE","FIELD1"); CHKERRQ(ierr);

  //Problem
  ierr = mField.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"FIELD1"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = mField.add_ents_to_finite_element_by_TETs(root_set,"TEST_FE"); CHKERRQ(ierr);


  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 3;
  ierr = mField.set_field_order(root_set,MBTET,"FIELD1",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"FIELD1",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"FIELD1",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"FIELD1",1); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = mField.build_finiteElementsPtr(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning 
  //partition
  ierr = mField.simple_partition_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finiteElementsPtr("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  Mat A;
  ierr = mField.MatCreateMPIAIJWithArrays("TEST_PROBLEM",&A); CHKERRQ(ierr);
  Vec D,F;
  ierr = mField.VecCreateGhost("TEST_PROBLEM",ROW,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("TEST_PROBLEM",COL,&D); CHKERRQ(ierr);

  {
  MyFunApprox function_evaluator;
  FieldApproximationH1<MyFunApprox> field_approximation(mField);
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
  ierr = KSPSetOperators(solver,A,A,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = mField.set_global_VecCreateGhost("TEST_PROBLEM",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);

  /*PostProcDisplacementsOnRefMesh fe_postproc(moab,"FIELD1");
  ierr = mField.loop_finiteElementsPtr("TEST_PROBLEM","TEST_FE",fe_postproc);  CHKERRQ(ierr);
  if(pcomm->rank()==0) {
    rval = fe_postproc.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }*/

  EntityHandle fe_meshset = mField.get_finite_element_meshset("TEST_FE");
  Range tets;
  rval = moab.get_entities_by_type(fe_meshset,MBTET,tets,true); CHKERR_PETSC(rval);
  Range tets_edges;
  rval = moab.get_adjacencies(tets,1,false,tets_edges,Interface::UNION); CHKERR(rval);
  EntityHandle edges_meshset;
  rval = moab.create_meshset(MESHSET_SET,edges_meshset); CHKERR_PETSC(rval);
  rval = moab.add_entities(edges_meshset,tets); CHKERR_PETSC(rval);
  rval = moab.add_entities(edges_meshset,tets_edges); CHKERR_PETSC(rval);
  rval = moab.convert_entities(edges_meshset,true,false,false); CHKERR_PETSC(rval);

  ProjectionFieldOn10NodeTet ent_method_field1_on_10nodeTet(mField,"FIELD1",true,false,"FIELD1");
  ierr = mField.loop_dofs("FIELD1",ent_method_field1_on_10nodeTet); CHKERRQ(ierr);
  ent_method_field1_on_10nodeTet.set_nodes = false;
  ierr = mField.loop_dofs("FIELD1",ent_method_field1_on_10nodeTet); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("TEST_PROBLEM","TEST_FE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  typedef tee_device<ostream, ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;

  ofstream ofs("forces_and_sources_testing_field_approximation.txt");
  TeeDevice tee(cout, ofs); 
  TeeStream my_split(tee);

  Range nodes;
  rval = moab.get_entities_by_type(0,MBVERTEX,nodes,true); CHKERR(rval);
  ublas::matrix<double> nodes_vals;
  nodes_vals.resize(nodes.size(),3);
  rval = moab.tag_get_data(
    ent_method_field1_on_10nodeTet.th,nodes,&*nodes_vals.data().begin()); CHKERR(rval);
  

  const double eps = 1e-4;

  my_split.precision(3);
  my_split.setf(std::ios::fixed);
  for(
    ublas::unbounded_array<double>::iterator it = nodes_vals.data().begin();
    it!=nodes_vals.data().end();it++) {
    *it = fabs(*it)<eps ? 0.0 : *it;
  }
  my_split << nodes_vals << endl;

  const MoFEMProblem *problemPtr;
  ierr = mField.get_problem("TEST_PROBLEM",&problemPtr); CHKERRQ(ierr);
  map<EntityHandle,double> m0,m1,m2;
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(problemPtr,dit)) {

    my_split.precision(3);
    my_split.setf(std::ios::fixed);
    double val = fabs(dit->get_FieldData())<eps ? 0.0 : dit->get_FieldData();
    my_split << dit->get_petsc_gloabl_dof_idx() << " " << val << endl;

  }

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


