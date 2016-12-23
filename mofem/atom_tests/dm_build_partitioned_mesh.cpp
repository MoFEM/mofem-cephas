/** \file dm_build_partitioned_mesh.cpp
  \brief Atom test for build mesh which is paragoned

  Data Manager (DM) MoFEM interface is used here for convinience

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


using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  //initialize petsc
  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

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

  //register new dm type, i.e. mofem
  DMType dm_name = "DMMOFEM";
  ierr = DMRegister_MoFEM(dm_name); CHKERRQ(ierr);

  //craete dm instance
  DM dm;
  ierr = DMCreate(PETSC_COMM_WORLD,&dm);CHKERRQ(ierr);
  ierr = DMSetType(dm,dm_name);CHKERRQ(ierr);

  //read mesh and create moab and mofem datastrutures
  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
  const char *option;
  option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);

  MoFEM::Core core(moab,PETSC_COMM_WORLD);
  MoFEM::Interface& m_field = core;

  EntityHandle root_set = moab.get_root_set();
  //add all entities to database, all of them will be used
  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(root_set,bit_level0); CHKERRQ(ierr);
  //define & build field
  int field_rank = 1;
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetInt(PETSC_NULL,"","-my_field_rank",&field_rank,&flg); CHKERRQ(ierr);
  #else
  ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_field_rank",&field_rank,&flg); CHKERRQ(ierr);
  #endif
  ierr = m_field.add_field("FIELD",H1,AINSWORTH_COLE_BASE,field_rank); CHKERRQ(ierr);
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"FIELD"); CHKERRQ(ierr);
  //set app. order
  int order = 4;
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD",1); CHKERRQ(ierr);
  //build data structures for fields
  ierr = m_field.build_fields(); CHKERRQ(ierr);

  //define & build finite elements
  ierr = m_field.add_finite_element("FE"); CHKERRQ(ierr);
  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("FE","FIELD"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("FE","FIELD"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("FE","FIELD"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_TETs(root_set,"FE"); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //set dm data structure which created mofem data structures
  ierr = DMMoFEMCreateMoFEM(dm,&m_field,dm_name,bit_level0); CHKERRQ(ierr);
  ierr = DMMoFEMSetSquareProblem(dm,PETSC_FALSE); CHKERRQ(ierr); // this is for testing (this problem has the same rows and cols)
  ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"FE"); CHKERRQ(ierr);
  ierr = DMSetUp(dm); CHKERRQ(ierr);

  // dump data to file, just to check if something was changed
  Mat m;
  ierr = DMCreateMatrix(dm,&m); CHKERRQ(ierr);

  // if(1) {
  //   MatView(m,PETSC_VIEWER_DRAW_WORLD);
  //   std::string wait;
  //   std::cin >> wait;
  // }

  ierr = m_field.partition_check_matrix_fill_in("DMMOFEM",-1,-1,1); CHKERRQ(ierr);


  PetscBool save_file = PETSC_TRUE;

  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetBool(PETSC_NULL,"","-my_save_fiele",&save_file,&flg); CHKERRQ(ierr);
  #else
  ierr = PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-my_save_fiele",&save_file,&flg); CHKERRQ(ierr);
  #endif
  if(save_file) {

    PetscViewer viewer;
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"dm_build_partitioned_mesh.txt",&viewer); CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO); CHKERRQ(ierr);
    MatView(m,viewer);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  }

  ierr = MatDestroy(&m); CHKERRQ(ierr);
  //destry dm
  ierr = DMDestroy(&dm); CHKERRQ(ierr);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  //finish work cleaning memory, getting statistics, ect.
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
