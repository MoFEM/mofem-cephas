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

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,PETSC_NULL,help);

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

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
    order = 3;
  }

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  //ref meshset ref level 0
  ierr = m_field.seed_ref_level_3D(0,0); CHKERRQ(ierr);

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  /***/
  //Define problem

  //Fields
  ierr = m_field.add_field("FIELD_A",H1,3); CHKERRQ(ierr);
  ierr = m_field.add_field("FIELD_B",H1,3); CHKERRQ(ierr);

  ierr = m_field.add_ents_to_field_by_TETs(0,"FIELD_A"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"FIELD_B"); CHKERRQ(ierr);

  ierr = m_field.set_field_order(0,MBTET,"FIELD_A",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"FIELD_A",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"FIELD_A",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"FIELD_A",1); CHKERRQ(ierr);

  ierr = m_field.set_field_order(0,MBTET,"FIELD_B",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"FIELD_B",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"FIELD_B",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"FIELD_B",1); CHKERRQ(ierr);

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);

  ierr = m_field.set_field(0,MBVERTEX,"FIELD_B"); CHKERRQ(ierr);
  ierr = m_field.set_field(1,MBVERTEX,"FIELD_A"); CHKERRQ(ierr);

  SeriesRecorder& recorder = core;

  ierr = recorder.add_series_recorder("TEST_SERIES1"); CHKERRQ(ierr);
  ierr = recorder.add_series_recorder("TEST_SERIES2"); CHKERRQ(ierr);

  //initialize
  ierr = recorder.initialize_series_recorder("TEST_SERIES1"); CHKERRQ(ierr);

  ierr = recorder.record_begin("TEST_SERIES1"); CHKERRQ(ierr);
  ierr = recorder.record_field("TEST_SERIES1","FIELD_B",bit_level0,bit_level0); CHKERRQ(ierr);
  ierr = recorder.record_end("TEST_SERIES1",1); CHKERRQ(ierr);

  ierr = m_field.field_axpy(1.,"FIELD_A","FIELD_B"); CHKERRQ(ierr);
  ierr = recorder.record_begin("TEST_SERIES1"); CHKERRQ(ierr);
  ierr = recorder.record_field("TEST_SERIES1","FIELD_B",bit_level0,bit_level0); CHKERRQ(ierr);

  ierr = recorder.initialize_series_recorder("TEST_SERIES2"); CHKERRQ(ierr);
  ierr = recorder.record_begin("TEST_SERIES2"); CHKERRQ(ierr);
  ierr = recorder.record_field("TEST_SERIES2","FIELD_A",bit_level0,bit_level0); CHKERRQ(ierr);
  ierr = recorder.record_field("TEST_SERIES2","FIELD_B",bit_level0,bit_level0); CHKERRQ(ierr);
  ierr = recorder.record_end("TEST_SERIES2",1); CHKERRQ(ierr);
  ierr = recorder.finalize_series_recorder("TEST_SERIES2"); CHKERRQ(ierr);

  ierr = recorder.record_end("TEST_SERIES1",2); CHKERRQ(ierr);

  //finalize
  ierr = recorder.finalize_series_recorder("TEST_SERIES1"); CHKERRQ(ierr);
  ierr = recorder.print_series_steps(); CHKERRQ(ierr);

  ierr = m_field.field_scale(2,"FIELD_A"); CHKERRQ(ierr);

  MoFEM::Core core2(moab);
  MoFEM::Interface& m_field2 = core2;

  //build field
  ierr = m_field2.build_fields(); CHKERRQ(ierr);

  typedef tee_device<std::ostream, std::ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;
  std::ofstream ofs("record_series_atom.txt");
  TeeDevice my_tee(std::cout, ofs);
  TeeStream my_split(my_tee);

  SeriesRecorder& recorder2 = core2;
  ierr = recorder2.print_series_steps(); CHKERRQ(ierr);

  my_split << "TEST_SERIES1" << std::endl;
  for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_((&recorder2),"TEST_SERIES1",sit)) {

    ierr = recorder2.load_series_data("TEST_SERIES1",sit->get_step_number()); CHKERRQ(ierr);

    my_split << "next step:\n";
    my_split << *sit << std::endl;

    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field2,"FIELD_B",dof)) {
      my_split << *(*dof) << "\n";
    }

  }

  my_split << "TEST_SERIES2" << std::endl;
  for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_((&recorder2),"TEST_SERIES2",sit)) {

    ierr = recorder2.load_series_data("TEST_SERIES2",sit->get_step_number()); CHKERRQ(ierr);

    my_split << "next step:\n";
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field2,"FIELD_A",dof)) {
      my_split << *(*dof) << "\n";
    }
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field2,"FIELD_B",dof)) {
      my_split << *(*dof) << "\n";
    }


  }

  PetscFinalize();
  return 0;

}
