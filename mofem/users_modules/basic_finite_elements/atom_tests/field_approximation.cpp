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

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

#include <MoFEM.hpp>
#include <Projection10NodeCoordsOnField.hpp>

using namespace MoFEM;
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <FieldApproximation.hpp>

#define HOON

static char help[] = "...\n\n";

/// Example approx. function
struct MyFunApprox {

  std::vector<ublas::vector<double> > result;

  std::vector<ublas::vector<double> >& operator()(double x, double y, double z) {
    result.resize(1);
    result[0].resize(3);
    (result[0])[0] = x;
    (result[0])[1] = y;
    (result[0])[2] = z*z;
    return result;
  }

};

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  BARRIER_RANK_START(pcomm)
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  BARRIER_RANK_END(pcomm)

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
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
  ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
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
  ierr = m_field.partition_simple_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  Mat A;
  ierr = m_field.MatCreateMPIAIJWithArrays("TEST_PROBLEM",&A); CHKERRQ(ierr);
  Vec D,F;
  ierr = m_field.VecCreateGhost("TEST_PROBLEM",ROW,&F); CHKERRQ(ierr);
  ierr = m_field.VecCreateGhost("TEST_PROBLEM",COL,&D); CHKERRQ(ierr);

  std::vector<Vec> vec_F;
  vec_F.push_back(F);

  {
    MyFunApprox function_evaluator;
    FieldApproximationH1 field_approximation(m_field);
    field_approximation.loopMatrixAndVectorVolume(
      "TEST_PROBLEM","TEST_FE","FIELD1",A,vec_F,function_evaluator);
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

  EntityHandle fe_meshset = m_field.get_finite_element_meshset("TEST_FE");
  Range tets;
  rval = moab.get_entities_by_type(fe_meshset,MBTET,tets,true); CHKERRQ_MOAB(rval);
  Range tets_edges;
  rval = moab.get_adjacencies(tets,1,false,tets_edges,moab::Interface::UNION); CHKERR_MOAB(rval);
  EntityHandle edges_meshset;
  rval = moab.create_meshset(MESHSET_SET,edges_meshset); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(edges_meshset,tets); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(edges_meshset,tets_edges); CHKERRQ_MOAB(rval);
  rval = moab.convert_entities(edges_meshset,true,false,false); CHKERRQ_MOAB(rval);

  ProjectionFieldOn10NodeTet ent_method_field1_on_10nodeTet(m_field,"FIELD1",true,false,"FIELD1");
  ierr = m_field.loop_dofs("FIELD1",ent_method_field1_on_10nodeTet); CHKERRQ(ierr);
  ent_method_field1_on_10nodeTet.setNodes = false;
  ierr = m_field.loop_dofs("FIELD1",ent_method_field1_on_10nodeTet); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERRQ_MOAB(rval);
    ierr = m_field.get_problem_finite_elements_entities("TEST_PROBLEM","TEST_FE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERRQ_MOAB(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERRQ_MOAB(rval);
  }

  typedef tee_device<std::ostream, std::ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;

  std::ofstream ofs("field_approximation.txt");
  TeeDevice tee(cout, ofs);
  TeeStream my_split(tee);

  Range nodes;
  rval = moab.get_entities_by_type(0,MBVERTEX,nodes,true); CHKERR_MOAB(rval);
  ublas::matrix<double> nodes_vals;
  nodes_vals.resize(nodes.size(),3);
  rval = moab.tag_get_data(
    ent_method_field1_on_10nodeTet.th,nodes,&*nodes_vals.data().begin()); CHKERR_MOAB(rval);

  const double eps = 1e-4;

  my_split.precision(3);
  my_split.setf(std::ios::fixed);
  for(
    ublas::unbounded_array<double>::iterator it = nodes_vals.data().begin();
    it!=nodes_vals.data().end();it++) {
    *it = fabs(*it)<eps ? 0.0 : *it;
  }
  my_split << nodes_vals << std::endl;

  const MoFEMProblem *problemPtr;
  ierr = m_field.get_problem("TEST_PROBLEM",&problemPtr); CHKERRQ(ierr);
  std::map<EntityHandle,double> m0,m1,m2;
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(problemPtr,dit)) {

    my_split.precision(3);
    my_split.setf(std::ios::fixed);
    double val = fabs(dit->get()->getFieldData())<eps ? 0.0 : dit->get()->getFieldData();
    my_split << dit->get()->getPetscGlobalDofIdx() << " " << val << std::endl;

  }

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
