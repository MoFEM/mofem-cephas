/** \file reading_med_file.cpp

  \brief Testing interface for reading and writing med files

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

  
  

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

    moab::Core mb_instance;
    moab::Interface& moab = mb_instance;
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

    //Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    MedInterface *med_interface_ptr;
    ierr = m_field.getInterface(med_interface_ptr); CHKERRG(ierr);

    ierr = med_interface_ptr->readMed(); CHKERRG(ierr);
    ierr = med_interface_ptr->medGetFieldNames(); CHKERRG(ierr);

    // read field tags
    for(
      std::map<std::string,MedInterface::FieldData>::iterator fit =
      med_interface_ptr->fieldNames.begin();
      fit!=med_interface_ptr->fieldNames.end();
      fit++
    ) {
      ierr = med_interface_ptr->readFields(
        med_interface_ptr->medFileName,fit->first,false,1
      ); CHKERRG(ierr);
    }


    PetscBool check = PETSC_TRUE;
    ierr = PetscOptionsGetBool(PETSC_NULL,"","-check",&check,PETSC_NULL); CHKERRG(ierr);

    int ii = 0;
    const int check_list[] = { 2163, 624, 65, 104};
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,cit)) {
      EntityHandle meshset = cit->getMeshset();
      int nb_ents;
      rval = moab.get_number_entities_by_handle(meshset,nb_ents,true); CHKERRG(rval);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Nb of ents in %s %d\n",cit->getName().c_str(),nb_ents); CHKERRG(ierr);
      if(check && nb_ents!=check_list[ii]) {
        SETERRQ2(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Wrong numbers of entities in meshset %d != %d",nb_ents,check_list[ii]);
      }
      ii++;
    }

    MeshsetsManager *meshset_manager_ptr;
    ierr = m_field.getInterface(meshset_manager_ptr); CHKERRG(ierr);
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_((*meshset_manager_ptr),BLOCKSET,mit)) {
      EntityHandle meshset = mit->getMeshset();
      std::string name = mit->getName();
      PetscPrintf(m_field.get_comm(),"Write mesh %s\n",name.c_str());
      rval = moab.write_file(
        ("out_"+mit->getName()+".vtk").c_str(),
        NULL,
        NULL,
        &meshset,1
      ); CHKERRG(rval);
    }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  ierr = PetscFinalize(); CHKERRG(ierr);

  return 0;
}
