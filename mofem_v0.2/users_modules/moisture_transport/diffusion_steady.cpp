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
#include <PostProcOnRefMesh.hpp>
#include <ThermalElement.hpp>
#include <Projection10NodeCoordsOnField.hpp>

#include <MoistureTransportElement.hpp>

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
  ierr = mField.add_field("CONC",H1,1); CHKERRQ(ierr);

  //Problem
  ierr = mField.add_problem("DIFFUSION_PROBLEM"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("DIFFUSION_PROBLEM",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"CONC"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }

  ierr = mField.set_field_order(root_set,MBTET,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"CONC",1); CHKERRQ(ierr);

  MoistureTransportElement moisture_element(mField);
  ierr = moisture_element.addDiffusionElement("DIFFUSION_PROBLEM","CONC"); CHKERRQ(ierr);
  ierr = moisture_element.addDiffusionFluxElement("DIFFUSION_PROBLEM","CONC"); CHKERRQ(ierr);

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
  ierr = mField.partition_problem("DIFFUSION_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("DIFFUSION_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("DIFFUSION_PROBLEM"); CHKERRQ(ierr);

  Vec F;
  ierr = mField.VecCreateGhost("DIFFUSION_PROBLEM",ROW,&F); CHKERRQ(ierr);
  Vec C;
  ierr = VecDuplicate(F,&C); CHKERRQ(ierr);
  Mat A;
  ierr = mField.MatCreateMPIAIJWithArrays("DIFFUSION_PROBLEM",&A); CHKERRQ(ierr);

  //New way to implement the dirichlet BCs from CUBIT Blockset (MASS_CONC)
  DirichletBCFromBlockSetFEMethodPreAndPostProc my_dirichlet_bc(mField,"CONC","MASS_CONC",A,C,F);

  //These operatore are the same for both thermal and diffusion problem and no need to replace for diffusion problem
  ierr = moisture_element.setThermalFiniteElementRhsOperators("CONC",F); CHKERRQ(ierr);
  ierr = moisture_element.setThermalFiniteElementLhsOperators("CONC",A); CHKERRQ(ierr);
  ierr = moisture_element.setThermalFluxFiniteElementRhsOperators("CONC",F); CHKERRQ(ierr);

  ierr = VecZeroEntries(C); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  //preproc
  ierr = mField.problem_basic_method_preProcess("DIFFUSION_PROBLEM",my_dirichlet_bc); CHKERRQ(ierr);
  ierr = mField.set_global_ghost_vector("DIFFUSION_PROBLEM",ROW,C,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = mField.loop_finite_elements("DIFFUSION_PROBLEM","DIFFUSION_FE",moisture_element.getLoopFeRhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("DIFFUSION_PROBLEM","DIFFUSION_FE",moisture_element.getLoopFeLhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("DIFFUSION_PROBLEM","DIFFUSION_FLUX_FE",moisture_element.getLoopFeFlux()); CHKERRQ(ierr);

  //postproc
  ierr = mField.problem_basic_method_postProcess("DIFFUSION_PROBLEM",my_dirichlet_bc); CHKERRQ(ierr);

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

  ierr = KSPSolve(solver,F,C); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = mField.problem_basic_method_preProcess("DIFFUSION_PROBLEM",my_dirichlet_bc); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_ghost_vector("DIFFUSION_PROBLEM",ROW,C,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //Range ref_edges;
  //ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBEDGE,ref_edges); CHKERRQ(ierr);
  //rval = moab.list_entities(ref_edges); CHKERR_PETSC(rval);
  //mField.list_dofs_by_field_name("TEMP");

  if(pcomm->rank()==0) {
    rval = moab.write_file("solution.h5m"); CHKERR_PETSC(rval);
  }

  /*EntityHandle fe_meshset = mField.get_finite_element_meshset("THERMAL_FE");
  Range tets;
  rval = moab.get_entities_by_type(fe_meshset,MBTET,tets,true); CHKERR_PETSC(rval);
  Range tets_edges;
  rval = moab.get_adjacencies(tets,1,false,tets_edges,Interface::UNION); CHKERR(rval);
  EntityHandle edges_meshset;
  rval = moab.create_meshset(MESHSET_SET,edges_meshset); CHKERR_PETSC(rval);
  rval = moab.add_entities(edges_meshset,tets); CHKERR_PETSC(rval);
  rval = moab.add_entities(edges_meshset,tets_edges); CHKERR_PETSC(rval);
  rval = moab.convert_entities(edges_meshset,true,false,false); CHKERR_PETSC(rval);*/

  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet(mField,"CONC",true,false,"CONC");
  ierr = mField.loop_dofs("CONC",ent_method_on_10nodeTet); CHKERRQ(ierr);
  ent_method_on_10nodeTet.set_nodes = false;
  ierr = mField.loop_dofs("CONC",ent_method_on_10nodeTet); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.get_problem_finite_elements_entities("DIFFUSION_PROBLEM","DIFFUSION_FE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&C); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
