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
#include "ThermalFEMethod.hpp"
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
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }


  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  FieldCore core(moab);
  FieldInterface& mField = core;

  //add filds
  ierr = mField.add_field("THERMAL_FIELD",H1,1); CHKERRQ(ierr);

  //add finite elements
  ierr = mField.add_finite_element("THERMAL_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("THERMAL_ELEM","THERMAL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("THERMAL_ELEM","THERMAL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("THERMAL_ELEM","THERMAL_FIELD"); CHKERRQ(ierr);

  //add problems 
  ierr = mField.add_problem("THERMAL_PROBLEM"); CHKERRQ(ierr);

  //define problems and finite elements
  ierr = mField.modify_problem_add_finite_element("THERMAL_PROBLEM","THERMAL_ELEM"); CHKERRQ(ierr);

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|UnknownCubitName,it)) {
    if(it->get_Cubit_name() == "THER_MAT") {
 
      //add ents to field and set app. order
      ierr = mField.add_ents_to_field_by_TETs(0,"THERMAL_FIELD"); CHKERRQ(ierr);

      //add finite elements entities
      ierr = mField.add_ents_to_finite_element_by_TETs(it->meshset,"THERMAL_ELEM",true); CHKERRQ(ierr);

    }
  }

  ierr = mField.set_field_order(0,MBVERTEX,"THERMAL_FIELD",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"THERMAL_FIELD",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"THERMAL_FIELD",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"THERMAL_FIELD",order); CHKERRQ(ierr);

  //set problem level
  ierr = mField.modify_problem_ref_level_add_bit("THERMAL_PROBLEM",bit_level0); CHKERRQ(ierr);

  //build fields
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elements
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = mField.partition_problem("THERMAL_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("THERMAL_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("THERMAL_PROBLEM"); CHKERRQ(ierr);

  //print bcs
  ierr = mField.printCubitDisplacementSet(); CHKERRQ(ierr);
  ierr = mField.printCubitPressureSet(); CHKERRQ(ierr);

  //create matrices and vectors
  Vec F,D;
  ierr = mField.VecCreateGhost("THERMAL_PROBLEM",Row,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("THERMAL_PROBLEM",Col,&D); CHKERRQ(ierr);
  Mat A;
  ierr = mField.MatCreateMPIAIJWithArrays("THERMAL_PROBLEM",&A); CHKERRQ(ierr);

  ThermalElem elem(mField,A,F);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("THERMAL_PROBLEM","THERMAL_ELEM",elem);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  //Matrix View
  /*{
    MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
    std::string wait;
    std::cin >> wait;
  }*/
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost("THERMAL_PROBLEM",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  Tag th_temp;
  double def_val = 0;
  rval = moab.tag_get_handle("TEMP",1,MB_TYPE_DOUBLE,th_temp,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERR_PETSC(rval);
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"THERMAL_FIELD",dof)) {
    EntityHandle ent = dof->get_ent();
    double val = dof->get_FieldData();
    rval = moab.tag_set_data(th_temp,&ent,1,&val); CHKERR_PETSC(rval);
  }

  if(pcomm->rank()==0) {
    rval = moab.write_file("solution.h5m"); CHKERR_PETSC(rval);
  }

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("THERMAL_PROBLEM","THERMAL_ELEM",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

}
