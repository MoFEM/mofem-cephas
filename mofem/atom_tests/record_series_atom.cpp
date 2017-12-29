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

  MoFEM::Core::Initialize(&argc,&argv,PETSC_NULL,help);

  try {

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //Reade parameters from line command
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetString(PETSC_NULL,"","-my_file",mesh_file_name,255,&flg); CHKERRG(ierr);
  #else
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRG(ierr);
  #endif
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  PetscInt order;
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetInt(PETSC_NULL,"","-my_order",&order,&flg); CHKERRG(ierr);
  #else
  ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRG(ierr);
  #endif
  if(flg != PETSC_TRUE) {
    order = 3;
  }

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRG(rval);
  ierr = m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(0,3,bit_level0); CHKERRG(ierr);
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByRefLevel(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRG(ierr);

  /***/
  //Define problem

  //Fields
  ierr = m_field.add_field("FIELD_A",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRG(ierr);
  ierr = m_field.add_field("FIELD_B",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRG(ierr);

  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"FIELD_A"); CHKERRG(ierr);
  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"FIELD_B"); CHKERRG(ierr);

  ierr = m_field.set_field_order(0,MBTET,"FIELD_A",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"FIELD_A",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"FIELD_A",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"FIELD_A",1); CHKERRG(ierr);

  ierr = m_field.set_field_order(0,MBTET,"FIELD_B",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"FIELD_B",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"FIELD_B",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"FIELD_B",1); CHKERRG(ierr);

  //build field
  ierr = m_field.build_fields(); CHKERRG(ierr);

  ierr = m_field.getInterface<FieldBlas>()->setField(0,MBVERTEX,"FIELD_B"); CHKERRG(ierr);
  ierr = m_field.getInterface<FieldBlas>()->setField(1,MBVERTEX,"FIELD_A"); CHKERRG(ierr);

  SeriesRecorder *recorder_ptr;
  ierr = m_field.getInterface(recorder_ptr); CHKERRG(ierr);

  ierr = recorder_ptr->add_series_recorder("TEST_SERIES1"); CHKERRG(ierr);
  ierr = recorder_ptr->add_series_recorder("TEST_SERIES2"); CHKERRG(ierr);

  //initialize
  ierr = recorder_ptr->initialize_series_recorder("TEST_SERIES1"); CHKERRG(ierr);

  ierr = recorder_ptr->record_begin("TEST_SERIES1"); CHKERRG(ierr);
  ierr = recorder_ptr->record_field("TEST_SERIES1","FIELD_B",bit_level0,bit_level0); CHKERRG(ierr);
  ierr = recorder_ptr->record_end("TEST_SERIES1",1); CHKERRG(ierr);

  ierr = m_field.getInterface<FieldBlas>()->fieldAxpy(1.,"FIELD_A","FIELD_B"); CHKERRG(ierr);
  ierr = recorder_ptr->record_begin("TEST_SERIES1"); CHKERRG(ierr);
  ierr = recorder_ptr->record_field("TEST_SERIES1","FIELD_B",bit_level0,bit_level0); CHKERRG(ierr);

  ierr = recorder_ptr->initialize_series_recorder("TEST_SERIES2"); CHKERRG(ierr);
  ierr = recorder_ptr->record_begin("TEST_SERIES2"); CHKERRG(ierr);
  ierr = recorder_ptr->record_field("TEST_SERIES2","FIELD_A",bit_level0,bit_level0); CHKERRG(ierr);
  ierr = recorder_ptr->record_field("TEST_SERIES2","FIELD_B",bit_level0,bit_level0); CHKERRG(ierr);
  ierr = recorder_ptr->record_end("TEST_SERIES2",1); CHKERRG(ierr);
  ierr = recorder_ptr->finalize_series_recorder("TEST_SERIES2"); CHKERRG(ierr);

  ierr = recorder_ptr->record_end("TEST_SERIES1",2); CHKERRG(ierr);

  //finalize
  ierr = recorder_ptr->finalize_series_recorder("TEST_SERIES1"); CHKERRG(ierr);
  ierr = recorder_ptr->print_series_steps(); CHKERRG(ierr);

  ierr = m_field.getInterface<FieldBlas>()->fieldScale(2,"FIELD_A"); CHKERRG(ierr);

  MoFEM::Core core2(moab);
  MoFEM::Interface& m_field2 = core2;

  //build field
  ierr = m_field2.build_fields(); CHKERRG(ierr);

  typedef tee_device<std::ostream, std::ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;
  std::ofstream ofs("record_series_atom.txt");
  TeeDevice my_tee(std::cout, ofs);
  TeeStream my_split(my_tee);

  SeriesRecorder *recorder2_ptr;
  ierr = m_field2.getInterface(recorder2_ptr); CHKERRG(ierr);
  ierr = recorder2_ptr->print_series_steps(); CHKERRG(ierr);

  const DofEntity_multiIndex *dofs_ptr;
  ierr = m_field.get_dofs(&dofs_ptr); CHKERRG(ierr);

  my_split << "TEST_SERIES1" << std::endl;
  for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder2_ptr,"TEST_SERIES1",sit)) {

    ierr = recorder2_ptr->load_series_data("TEST_SERIES1",sit->get_step_number()); CHKERRG(ierr);

    my_split << "next step:\n";
    my_split << *sit << std::endl;

    {
      DofEntity_multiIndex_uid_view dofs_view;
      for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field2,"FIELD_B",dof)) {
        dofs_view.insert(*dof);
      }
      for(
        DofEntity_multiIndex_uid_view::iterator
        dit=dofs_view.begin();dit!=dofs_view.end();dit++
      ) {
        my_split << **dit << endl;
      }
    }

  }

  my_split << "TEST_SERIES2" << std::endl;
  for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder2_ptr,"TEST_SERIES2",sit)) {

    ierr = recorder2_ptr->load_series_data("TEST_SERIES2",sit->get_step_number()); CHKERRG(ierr);

    my_split << "next step:\n";
    {
      DofEntity_multiIndex_uid_view dofs_view;
      for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field2,"FIELD_A",dof)) {
        dofs_view.insert(*dof);
      }
      for(
        DofEntity_multiIndex_uid_view::iterator
        dit=dofs_view.begin();dit!=dofs_view.end();dit++
      ) {
        my_split << **dit << endl;
      }
    }
    {
      DofEntity_multiIndex_uid_view dofs_view;
      for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field2,"FIELD_B",dof)) {
        dofs_view.insert(*dof);
      }
      for(
        DofEntity_multiIndex_uid_view::iterator
        dit=dofs_view.begin();dit!=dofs_view.end();dit++
      ) {
        my_split << **dit << endl;
      }
    }


  }

  } CATCH_ERRORS;

  MoFEM::Core::Finalize();
  return 0;

}
