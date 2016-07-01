/** \file build_problems.cpp

  \brief Atom test for building composite problem

  \bug Not verifying what if partitioned mesh is loaded.

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

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

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

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("F1",L2,1); CHKERRQ(ierr);
  ierr = m_field.add_field("F2",HDIV,1); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"F1"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"F2"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 4;
  ierr = m_field.set_field_order(root_set,MBTET,"F1",order-1,2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTET,"F2",order,2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"F2",order,2); CHKERRQ(ierr);

  ierr = m_field.build_fields(); CHKERRQ(ierr);

  order = 2;
  ierr = m_field.set_field_order(root_set,MBTET,"F1",order-1,2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTET,"F2",order,2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"F2",order,2); CHKERRQ(ierr);

  ierr = m_field.clear_inactive_dofs(); CHKERRQ(ierr);
  ierr = m_field.build_fields(); CHKERRQ(ierr);

  const DofEntity_multiIndex *dofs_ptr;
  ierr = m_field.get_dofs(&dofs_ptr); CHKERRQ(ierr);
  const unsigned int expected_size = 744;
  if(dofs_ptr->size()!=expected_size) {
    SETERRQ1(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "Data inconsistency, should %d dofs",
      expected_size
    );
  }

  //finish work cleaning memory, getting statistics, ect.
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
