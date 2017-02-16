/** \file thermal_steady.cpp
  \ingroup mofem_thermal_elem
  \brief Example of steady thermal analysis
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


#include <BasicFiniteElements.hpp>
using namespace MoFEM;

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

static char help[] = "...\n\n";

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

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("TEMP",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("THERMAL_PROBLEM"); CHKERRQ(ierr);

  //set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("THERMAL_PROBLEM",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"TEMP"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 2;
  }

  ierr = m_field.set_field_order(root_set,MBTET,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"TEMP",1); CHKERRQ(ierr);

  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  ThermalElement thermal_elements(m_field);
  ierr = thermal_elements.addThermalElements("TEMP"); CHKERRQ(ierr);
  ierr = thermal_elements.addThermalFluxElement("TEMP"); CHKERRQ(ierr);
  ierr = thermal_elements.addThermalConvectionElement("TEMP"); CHKERRQ(ierr);

  ierr = m_field.modify_problem_add_finite_element("THERMAL_PROBLEM","THERMAL_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("THERMAL_PROBLEM","THERMAL_FLUX_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("THERMAL_PROBLEM","THERMAL_CONVECTION_FE"); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);


  ProblemsManager *prb_mng_ptr;
  ierr = m_field.query_interface(prb_mng_ptr); CHKERRQ(ierr);
  //build problem
  ierr = prb_mng_ptr->buildProblem("THERMAL_PROBLEM",true); CHKERRQ(ierr);

  Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  /****/
  //mesh partitioning
  //partition
  ierr = prb_mng_ptr->partitionProblem("THERMAL_PROBLEM"); CHKERRQ(ierr);
  ierr = prb_mng_ptr->partitionFiniteElements("THERMAL_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = prb_mng_ptr->partitionGhostDofs("THERMAL_PROBLEM"); CHKERRQ(ierr);

  Vec F;
  ierr = m_field.VecCreateGhost("THERMAL_PROBLEM",ROW,&F); CHKERRQ(ierr);
  Vec T;
  ierr = VecDuplicate(F,&T); CHKERRQ(ierr);
  Mat A;
  ierr = m_field.MatCreateMPIAIJWithArrays("THERMAL_PROBLEM",&A); CHKERRQ(ierr);

  TemperatureBCFEMethodPreAndPostProc my_dirichlet_bc(m_field,"TEMP",A,T,F);
  ierr = thermal_elements.setThermalFiniteElementRhsOperators("TEMP",F); CHKERRQ(ierr);
  ierr = thermal_elements.setThermalFiniteElementLhsOperators("TEMP",A); CHKERRQ(ierr);
  ierr = thermal_elements.setThermalFluxFiniteElementRhsOperators("TEMP",F); CHKERRQ(ierr);
  ierr = thermal_elements.setThermalConvectionFiniteElementLhsOperators("TEMP",A); CHKERRQ(ierr);
  ierr = thermal_elements.setThermalConvectionFiniteElementRhsOperators("TEMP",F); CHKERRQ(ierr);

  ierr = VecZeroEntries(T); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  //preproc
  ierr = m_field.problem_basic_method_preProcess("THERMAL_PROBLEM",my_dirichlet_bc); CHKERRQ(ierr);
  ierr = m_field.set_global_ghost_vector("THERMAL_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = m_field.loop_finite_elements("THERMAL_PROBLEM","THERMAL_FE",thermal_elements.getLoopFeRhs()); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("THERMAL_PROBLEM","THERMAL_FE",thermal_elements.getLoopFeLhs()); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("THERMAL_PROBLEM","THERMAL_FLUX_FE",thermal_elements.getLoopFeFlux()); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("THERMAL_PROBLEM","THERMAL_CONVECTION_FE",thermal_elements.getLoopFeConvectionRhs()); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("THERMAL_PROBLEM","THERMAL_CONVECTION_FE",thermal_elements.getLoopFeConvectionLhs()); CHKERRQ(ierr);

  //postproc
  ierr = m_field.problem_basic_method_postProcess("THERMAL_PROBLEM",my_dirichlet_bc); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = VecScale(F,-1); CHKERRQ(ierr);

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,T); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = m_field.problem_basic_method_preProcess("THERMAL_PROBLEM",my_dirichlet_bc); CHKERRQ(ierr);

  //Save data on mesh
  ierr = m_field.set_global_ghost_vector("THERMAL_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //Range ref_edges;
  //ierr = m_field.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBEDGE,ref_edges); CHKERRQ(ierr);
  //rval = moab.list_entities(ref_edges); CHKERRQ_MOAB(rval);
  //m_field.list_dofs_by_field_name("TEMP");

  if(pcomm->rank()==0) {
    rval = moab.write_file("solution.h5m"); CHKERRQ_MOAB(rval);
  }

  /*EntityHandle fe_meshset = m_field.get_finite_element_meshset("THERMAL_FE");
  Range tets;
  rval = moab.get_entities_by_type(fe_meshset,MBTET,tets,true); CHKERRQ_MOAB(rval);
  Range tets_edges;
  rval = moab.get_adjacencies(tets,1,false,tets_edges,moab::Interface::UNION); CHKERRQ_MOAB(rval);
  EntityHandle edges_meshset;
  rval = moab.create_meshset(MESHSET_SET,edges_meshset); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(edges_meshset,tets); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(edges_meshset,tets_edges); CHKERRQ_MOAB(rval);
  rval = moab.convert_entities(edges_meshset,true,false,false); CHKERRQ_MOAB(rval);*/

  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet(m_field,"TEMP",true,false,"TEMP");
  ierr = m_field.loop_dofs("TEMP",ent_method_on_10nodeTet); CHKERRQ(ierr);
  ent_method_on_10nodeTet.setNodes = false;
  ierr = m_field.loop_dofs("TEMP",ent_method_on_10nodeTet); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERRQ_MOAB(rval);
    ierr = m_field.get_problem_finite_elements_entities("THERMAL_PROBLEM","THERMAL_FE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERRQ_MOAB(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERRQ_MOAB(rval);
  }

  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&T); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
