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
#include "FEMethod_UpLevelStudent.hpp"
#include "PotentialFlowFEMethod.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";


int main(int argc, char *argv[]) {

  try {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  /*if(rank==0) {
    EntityHandle dummy_meshset;
    rval = moab.create_meshset(MESHSET_SET,dummy_meshset); CHKERR_PETSC(rval);
  }*/

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  PetscInt order_laplacian;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order_laplacian",&order_laplacian,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order_laplacian = 1;
  }
  PetscInt order_pressure;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order_pressure",&order_pressure,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order_pressure = 1;
  }

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  FieldCore core(moab);
  FieldInterface& mField = core;

  //add filds
  ierr = mField.add_field("POTENTIAL_FIELD",H1,1); CHKERRQ(ierr);
  ierr = mField.add_field("PRESSURE_FIELD",H1,1); CHKERRQ(ierr);

  //add finite elements
  ierr = mField.add_finite_element("LAPLACIAN_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("LAPLACIAN_ELEM","POTENTIAL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("LAPLACIAN_ELEM","POTENTIAL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("LAPLACIAN_ELEM","POTENTIAL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("LAPLACIAN_ELEM","PRESSURE_FIELD"); CHKERRQ(ierr);

  ierr = mField.add_finite_element("PRESSURE_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("PRESSURE_ELEM","PRESSURE_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("PRESSURE_ELEM","PRESSURE_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("PRESSURE_ELEM","POTENTIAL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("PRESSURE_ELEM","PRESSURE_FIELD"); CHKERRQ(ierr);

  ierr = mField.add_finite_element("BERNOULLY_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("BERNOULLY_ELEM","PRESSURE_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("BERNOULLY_ELEM","POTENTIAL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("BERNOULLY_ELEM","PRESSURE_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("BERNOULLY_ELEM","POTENTIAL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("BERNOULLY_ELEM","PRESSURE_FIELD"); CHKERRQ(ierr);

  //add problems 
  ierr = mField.add_problem("LAPLACIAN_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.add_problem("PRESSURE_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.add_problem("BERNOULLY_PROBLEM"); CHKERRQ(ierr);

  //define problems and finite elements
  ierr = mField.modify_problem_add_finite_element("LAPLACIAN_PROBLEM","LAPLACIAN_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("PRESSURE_PROBLEM","PRESSURE_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("BERNOULLY_PROBLEM","BERNOULLY_ELEM"); CHKERRQ(ierr);

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|UnknownCubitName,it)) {
    if(it->get_Cubit_name() == "PotentialFlow") {
 
      //add ents to field and set app. order_laplacian
      ierr = mField.add_ents_to_field_by_TETs(it->meshset,"POTENTIAL_FIELD"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_field_by_TETs(it->meshset,"PRESSURE_FIELD"); CHKERRQ(ierr);

      //add finite elements entities
      ierr = mField.add_ents_to_finite_element_by_TETs(it->meshset,"LAPLACIAN_ELEM",true); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TETs(it->meshset,"PRESSURE_ELEM",true); CHKERRQ(ierr);

    }
  }

  //laplacian
  ierr = mField.set_field_order(0,MBVERTEX,"POTENTIAL_FIELD",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"POTENTIAL_FIELD",order_laplacian); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"POTENTIAL_FIELD",order_laplacian); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"POTENTIAL_FIELD",order_laplacian); CHKERRQ(ierr);
  //pressure
  ierr = mField.set_field_order(0,MBVERTEX,"PRESSURE_FIELD",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"PRESSURE_FIELD",order_pressure); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"PRESSURE_FIELD",order_pressure); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"PRESSURE_FIELD",order_pressure); CHKERRQ(ierr);

  //set problem level
  ierr = mField.modify_problem_ref_level_add_bit("LAPLACIAN_PROBLEM",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("PRESSURE_PROBLEM",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("BERNOULLY_PROBLEM",bit_level0); CHKERRQ(ierr);

  //build fields
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elements
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = mField.partition_problem("LAPLACIAN_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("LAPLACIAN_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("LAPLACIAN_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_problem("PRESSURE_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("PRESSURE_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("PRESSURE_PROBLEM"); CHKERRQ(ierr);

  //print bcs
  ierr = mField.printCubitDisplacementSet(); CHKERRQ(ierr);
  ierr = mField.printCubitPressureSet(); CHKERRQ(ierr);

  //**** solve lapalacian problem ****

  //create matrices and vectors
  Vec F_Laplacian,D_Laplacian;
  ierr = mField.VecCreateGhost("LAPLACIAN_PROBLEM",Row,&F_Laplacian); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("LAPLACIAN_PROBLEM",Col,&D_Laplacian); CHKERRQ(ierr);
  Mat A_Laplacian;
  ierr = mField.MatCreateMPIAIJWithArrays("LAPLACIAN_PROBLEM",&A_Laplacian); CHKERRQ(ierr);

  LaplacianElem elem(mField,A_Laplacian,F_Laplacian);
  ierr = VecZeroEntries(F_Laplacian); CHKERRQ(ierr);
  ierr = MatZeroEntries(A_Laplacian); CHKERRQ(ierr);
  ierr = elem.set_bodyVelocity(1,0,0.0);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("LAPLACIAN_PROBLEM","LAPLACIAN_ELEM",elem);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A_Laplacian,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_Laplacian,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F_Laplacian); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F_Laplacian); CHKERRQ(ierr);

  //Matrix View
  /*{
    MatView(A_Laplacian,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
    std::string wait;
    std::cin >> wait;
  }*/

  //Solver
  KSP solver_laplacian;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver_laplacian); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver_laplacian,A_Laplacian,A_Laplacian,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver_laplacian); CHKERRQ(ierr);
  ierr = KSPSetUp(solver_laplacian); CHKERRQ(ierr);

  ierr = KSPSolve(solver_laplacian,F_Laplacian,D_Laplacian); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D_Laplacian,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D_Laplacian,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost("LAPLACIAN_PROBLEM",Row,D_Laplacian,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = KSPDestroy(&solver_laplacian); CHKERRQ(ierr);
  ierr = VecDestroy(&F_Laplacian); CHKERRQ(ierr);
  ierr = VecDestroy(&D_Laplacian); CHKERRQ(ierr);
  ierr = MatDestroy(&A_Laplacian); CHKERRQ(ierr);

  Tag th_phi;
  double def_val = 0;
  rval = moab.tag_get_handle("PHI",1,MB_TYPE_DOUBLE,th_phi,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERR_PETSC(rval);
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"POTENTIAL_FIELD",dof)) {
    EntityHandle ent = dof->get_ent();
    double val = dof->get_FieldData();
    rval = moab.tag_set_data(th_phi,&ent,1,&val); CHKERR_PETSC(rval);
  }

  /*if(pcomm->rank()==0) {
    rval = moab.write_file("solution_laplacian.h5m"); CHKERR_PETSC(rval);
  }*/

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("LAPLACIAN_PROBLEM","LAPLACIAN_ELEM",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out_laplacian.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  //**** solve calculate pressure ****

  Vec F_Preassure,D_Preassure;
  ierr = mField.VecCreateGhost("PRESSURE_PROBLEM",Row,&F_Preassure); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("PRESSURE_PROBLEM",Col,&D_Preassure); CHKERRQ(ierr);
  Mat A_Preassure;
  ierr = mField.MatCreateMPIAIJWithArrays("PRESSURE_PROBLEM",&A_Preassure); CHKERRQ(ierr);

  SteadyBernoullyElem steady(mField,A_Preassure,F_Preassure);
  ierr = VecZeroEntries(F_Preassure); CHKERRQ(ierr);
  ierr = MatZeroEntries(A_Preassure); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("PRESSURE_PROBLEM","PRESSURE_ELEM",steady);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A_Preassure,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_Preassure,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F_Preassure); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F_Preassure); CHKERRQ(ierr);

  //Solver
  KSP solver_pressure;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver_pressure); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver_pressure,A_Preassure,A_Preassure,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver_pressure); CHKERRQ(ierr);
  ierr = KSPSetUp(solver_pressure); CHKERRQ(ierr);

  ierr = KSPSolve(solver_pressure,F_Preassure,D_Preassure); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D_Preassure,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D_Preassure,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost("PRESSURE_PROBLEM",Row,D_Preassure,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = KSPDestroy(&solver_pressure); CHKERRQ(ierr);
  ierr = VecDestroy(&F_Preassure); CHKERRQ(ierr);
  ierr = VecDestroy(&D_Preassure); CHKERRQ(ierr);
  ierr = MatDestroy(&A_Preassure); CHKERRQ(ierr);

  /*for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|UnknownCubitName,it)) {
    if(it->get_Cubit_name() != "ZeroPressure") continue;
    Range nodes;
    rval = moab.get_entities_by_type(it->meshset,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
    double pressure_shift = 0;
    for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,"PRESSURE_FIELD",*nodes.begin(),dof)) {
      pressure_shift = dof->get_FieldData();
    }
    for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(mField,"PRESSURE_FIELD",MBVERTEX,dof)) {
      dof->get_FieldData() -= pressure_shift;
    }
  }*/

  SurfaceForces forces_from_pressures(mField);
  Vec total_forces;
  ierr = forces_from_pressures.create_totalForceVector(&total_forces); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("PRESSURE_PROBLEM","PRESSURE_ELEM",forces_from_pressures);  CHKERRQ(ierr);
  ierr = VecView(total_forces,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecDestroy(&total_forces); CHKERRQ(ierr);

  Tag th_p;
  rval = moab.tag_get_handle("P",1,MB_TYPE_DOUBLE,th_p,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERR_PETSC(rval);
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"PRESSURE_FIELD",dof)) {
    EntityHandle ent = dof->get_ent();
    double val = dof->get_FieldData();
    rval = moab.tag_set_data(th_p,&ent,1,&val); CHKERR_PETSC(rval);
  }

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("PRESSURE_PROBLEM","PRESSURE_ELEM",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out_pressure.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  /*if(pcomm->rank()==0) {
    rval = moab.write_file("solution_pressure.h5m"); CHKERR_PETSC(rval);
  }*/

  /*PostProcPotentialFlowOnRefMesh post_proc_on_ref_mesh(moab);
  ierr = mField.loop_finite_elements("PRESSURE_PROBLEM","PRESSURE_ELEM",post_proc_on_ref_mesh);  CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    rval = post_proc_on_ref_mesh.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }*/

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

}
