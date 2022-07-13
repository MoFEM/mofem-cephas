/** \file reading_med_file.cpp

  \brief Testing interface for reading and writing med files

*/



#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    // Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    MedInterface *med_interface_ptr;
    CHKERR m_field.getInterface(med_interface_ptr);

    CHKERR med_interface_ptr->readMed();
    CHKERR med_interface_ptr->medGetFieldNames();

    // read field tags
    for (std::map<std::string, MedInterface::FieldData>::iterator fit =
             med_interface_ptr->fieldNames.begin();
         fit != med_interface_ptr->fieldNames.end(); fit++) {
      CHKERR med_interface_ptr->readFields(med_interface_ptr->medFileName,
                                           fit->first, false, 1);
    }

    PetscBool check = PETSC_TRUE;
    CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-check", &check, PETSC_NULL);

    int ii = 0;
    const int check_list[] = {2163, 624, 65, 104};
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, BLOCKSET, cit)) {
      EntityHandle meshset = cit->getMeshset();
      int nb_ents;
      CHKERR moab.get_number_entities_by_handle(meshset, nb_ents, true);
      CHKERR PetscPrintf(PETSC_COMM_WORLD, "Nb of ents in %s %d\n",
                         cit->getName().c_str(), nb_ents);
      CHKERRG(ierr);
      if (check && nb_ents != check_list[ii]) {
        SETERRQ2(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                 "Wrong numbers of entities in meshset %d != %d", nb_ents,
                 check_list[ii]);
      }
      ii++;
    }

    MeshsetsManager *meshset_manager_ptr;
    CHKERR m_field.getInterface(meshset_manager_ptr);
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_((*meshset_manager_ptr),
                                                 BLOCKSET, mit)) {
      EntityHandle meshset = mit->getMeshset();
      std::string name = mit->getName();
      PetscPrintf(m_field.get_comm(), "Write mesh %s\n", name.c_str());
      CHKERR moab.write_file(("out_" + mit->getName() + ".vtk").c_str(), NULL,
                             NULL, &meshset, 1);
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
