/** \file dm_mofem.cpp

  \brief Atom test for Data Manager Interface in MoFEM

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

  
  

  //initialize petsc
  MoFEM::Core::Initialize(&argc,&argv,(char *)0,help);

  try {

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
  const char *option;
  option = "";

  //register new dm type, i.e. mofem
  DMType dm_name = "DMMOFEM";
  ierr = DMRegister_MoFEM(dm_name); CHKERRG(ierr);

  // Create dm instance
  DM dm;
  ierr = DMCreate(PETSC_COMM_WORLD,&dm);CHKERRG(ierr);
  ierr = DMSetType(dm,dm_name);CHKERRG(ierr);

  //read mesh and create moab and mofem datastrutures
  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  EntityHandle root_set = moab.get_root_set();
  //add all entities to database, all of them will be used
  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(root_set,3,bit_level0); CHKERRG(ierr);
  //define & build field
  const int field_rank = 1;
  ierr = m_field.add_field("FIELD",H1,AINSWORTH_LEGENDRE_BASE,field_rank); CHKERRG(ierr);
  //add entities to field
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"FIELD"); CHKERRG(ierr);
  //set app. order
  int order = 4;
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD",1); CHKERRG(ierr);
  //build data structures for fields
  ierr = m_field.build_fields(); CHKERRG(ierr);

  //define & build finite elements
  ierr = m_field.add_finite_element("FE"); CHKERRG(ierr);
  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("FE","FIELD"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_col("FE","FIELD"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("FE","FIELD"); CHKERRG(ierr);
  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"FE"); CHKERRG(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRG(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRG(ierr);

  //set dm data structure which created mofem data structures
  ierr = DMMoFEMCreateMoFEM(dm,&m_field,dm_name,bit_level0); CHKERRG(ierr);
  ierr = DMSetFromOptions(dm); CHKERRG(ierr);
  ierr = DMMoFEMAddElement(dm,"FE"); CHKERRG(ierr);
  ierr = DMSetUp(dm); CHKERRG(ierr);

  Mat m;
  Vec l,g;

  ierr = DMCreateGlobalVector(dm,&g); CHKERRG(ierr);
  ierr = DMCreateLocalVector(dm,&l); CHKERRG(ierr);
  ierr = DMCreateMatrix(dm,&m); CHKERRG(ierr);

  //glob loc
  ierr = VecSet(g,1.1); CHKERRG(ierr);
  ierr = VecGhostUpdateBegin(g,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(g,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);

  ierr = DMGlobalToLocalBegin(dm,g,ADD_VALUES,l); CHKERRG(ierr);
  ierr = DMGlobalToLocalEnd(dm,g,ADD_VALUES,l); CHKERRG(ierr);

  //loc glob
  ierr = DMLocalToGlobalBegin(dm,l,ADD_VALUES,g); CHKERRG(ierr);
  ierr = DMLocalToGlobalEnd(dm,l,ADD_VALUES,g); CHKERRG(ierr);

  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"dm_mofem.txt",&viewer); CHKERRG(ierr);
  const double chop = 1e-4;
  ierr = VecChop(g,chop); CHKERRG(ierr);
  VecView(g,viewer);
  ierr = PetscViewerDestroy(&viewer); CHKERRG(ierr);

  ierr = VecDestroy(&g); CHKERRG(ierr);
  ierr = VecDestroy(&l); CHKERRG(ierr);
  ierr = MatDestroy(&m); CHKERRG(ierr);
  //destry dm
  ierr = DMDestroy(&dm); CHKERRG(ierr);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  //finish work cleaning memory, getting statistics, ect.
  ierr = MoFEM::Core::Finalize(); CHKERRG(ierr);

  return 0;

}
