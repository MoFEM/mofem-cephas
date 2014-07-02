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

using namespace MoFEM;

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,PETSC_NULL,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //Reade parameters from line command
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

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  /***/
  //Define problem

  //Fields
  ierr = mField.add_field("FIELD_A",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("FIELD_B",H1,3); CHKERRQ(ierr);

  ierr = mField.add_ents_to_field_by_TETs(0,"FIELD_A"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"FIELD_B"); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTET,"FIELD_A",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"FIELD_A",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"FIELD_A",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"FIELD_A",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTET,"FIELD_B",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"FIELD_B",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"FIELD_B",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"FIELD_B",1); CHKERRQ(ierr);

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  ierr = mField.set_field(0,MBVERTEX,"FIELD_B"); CHKERRQ(ierr);
  ierr = mField.set_field(1,MBVERTEX,"FIELD_A"); CHKERRQ(ierr);

  ierr = mField.add_series_recorder("TEST_SERIES1"); CHKERRQ(ierr);
  ierr = mField.add_series_recorder("TEST_SERIES2"); CHKERRQ(ierr);

  //initialize
  ierr = mField.initialize_series_recorder("TEST_SERIES1"); CHKERRQ(ierr);

  ierr = mField.record_begin("TEST_SERIES1"); CHKERRQ(ierr);
  ierr = mField.record_field("TEST_SERIES1","FIELD_B",bit_level0,bit_level0); CHKERRQ(ierr);
  ierr = mField.record_end("TEST_SERIES1"); CHKERRQ(ierr);

  ierr = mField.field_axpy(1.,"FIELD_A","FIELD_B"); CHKERRQ(ierr);
  ierr = mField.record_begin("TEST_SERIES1"); CHKERRQ(ierr);
  ierr = mField.record_field("TEST_SERIES1","FIELD_B",bit_level0,bit_level0); CHKERRQ(ierr);

  ierr = mField.initialize_series_recorder("TEST_SERIES2"); CHKERRQ(ierr);
  ierr = mField.record_begin("TEST_SERIES2"); CHKERRQ(ierr);
  ierr = mField.record_field("TEST_SERIES2","FIELD_A",bit_level0,bit_level0); CHKERRQ(ierr);
  ierr = mField.record_field("TEST_SERIES2","FIELD_B",bit_level0,bit_level0); CHKERRQ(ierr);
  ierr = mField.record_end("TEST_SERIES2"); CHKERRQ(ierr);
  ierr = mField.finalize_series_recorder("TEST_SERIES2"); CHKERRQ(ierr);

  ierr = mField.record_end("TEST_SERIES1"); CHKERRQ(ierr);

  //finalize
  ierr = mField.finalize_series_recorder("TEST_SERIES1"); CHKERRQ(ierr);

  ierr = mField.print_series_steps(); CHKERRQ(ierr);

  ierr = mField.field_scale(2,"FIELD_A"); CHKERRQ(ierr);

  FieldCore core2(moab);
  FieldInterface& mField2 = core2;
  ierr = mField2.print_series_steps(); CHKERRQ(ierr);

  //build field
  ierr = mField2.build_fields(); CHKERRQ(ierr);

  typedef tee_device<ostream, ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;
  ofstream ofs("record_series_atom.txt");
  TeeDevice my_tee(cout, ofs); 
  TeeStream my_split(my_tee);

  my_split << "TEST_SERIES1" << endl;
  for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(mField2,"TEST_SERIES1",sit)) {

    ierr = mField.load_series_data("TEST_SERIES1",sit->get_step_number()); CHKERRQ(ierr);

    my_split << "next step:\n";
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField2,"FIELD_B",dof)) {
      my_split << *dof << "\n";
    }

  }

  my_split << "TEST_SERIES2" << endl;
  for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(mField2,"TEST_SERIES2",sit)) {

    ierr = mField.load_series_data("TEST_SERIES2",sit->get_step_number()); CHKERRQ(ierr);

    my_split << "next step:\n";
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField2,"FIELD_A",dof)) {
      my_split << *dof << "\n";
    }
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField2,"FIELD_B",dof)) {
      my_split << *dof << "\n";
    }


  }

  PetscFinalize();
  return 0;

}

