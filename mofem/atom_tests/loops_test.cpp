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

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

    moab::Core mb_instance;
    moab::Interface& moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    //Read parameters from line command
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

    //Read mesh to MOAB
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

    //Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    //Open mesh_file_name.txt for writing
    std::ofstream myfile;
    myfile.open ((std::string(mesh_file_name)+".txt").c_str());

    std::cout << "<<<< All BLOCKSETs, SIDESETs and NODESETs >>>>>" << std::endl;
    for(_IT_CUBITMESHSETS_FOR_LOOP_(m_field,it)) {
      std::cout<< it->getName() << std::endl;
      myfile << it->getName() << std::endl;
    }
    std::cout << "<<<< BLOCKSETs >>>>>" << std::endl;
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
      std::cout<< it->getName() << std::endl;
      myfile << it->getName() << std::endl;
    }
    std::cout << "<<<< NODESETs >>>>>" << std::endl;
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,NODESET,it)) {
      std::cout<< it->getName() << std::endl;
      myfile << it->getName() << std::endl;
    }
    std::cout << "<<<< SIDESETs >>>>>" << std::endl;
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
      std::cout<< it->getName() << std::endl;
      myfile << it->getName() << std::endl;
    }
    std::cout <<"<<<< MeshSet of Name Moon >>>>" << std::endl;
    for (_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field,"Moon",it)){
      std::cout << it->getName() << std::endl;
      myfile << it->getName() << std::endl;
      if(it->get_cubit_bc_type_ulong() & BLOCKSET) {
        std::cout << "BLOCKSET" << std::endl;
        myfile << "BLOCKSET" << std::endl;
      }
      if(it->get_cubit_bc_type_ulong() & SIDESET) {
        std::cout << "SIDESET" << std::endl;
        myfile << "SIDESET" << std::endl;
      }
      if(it->get_cubit_bc_type_ulong() & NODESET) {
        std::cout << "NODESET" << std::endl;
        myfile << "NODESET" << std::endl;
      }
    }

    //Close mesh_file_name.txt
    myfile.close();

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFinalize();

}
