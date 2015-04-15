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
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //Read parameters from line command
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Read mesh to MOAB
  const char *option;
  option = "";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;
    
  //Open mesh_file_name.txt for writing
  ofstream myfile;
  myfile.open ((string(mesh_file_name)+".txt").c_str());

  cout << "<<<< All BLOCKSETs, SIDESETs and NODESETs >>>>>" << endl;
  for(_IT_CUBITMESHSETS_FOR_LOOP_(m_field,it)) {
    cout<< it->get_Cubit_name() << endl;
    myfile << it->get_Cubit_name() << endl;
  }
  cout << "<<<< BLOCKSETs >>>>>" << endl;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
    cout<< it->get_Cubit_name() << endl;
    myfile << it->get_Cubit_name() << endl;
  }
  cout << "<<<< NODESETs >>>>>" << endl;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,NODESET,it)) {
    cout<< it->get_Cubit_name() << endl;
    myfile << it->get_Cubit_name() << endl;
  }
  cout << "<<<< SIDESETs >>>>>" << endl;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    cout<< it->get_Cubit_name() << endl;
    myfile << it->get_Cubit_name() << endl;
  }
  cout <<"<<<< MeshSet of Name Moon >>>>" << endl;
  for (_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field,"Moon",it)){
    cout << it->get_Cubit_name() << endl;
    myfile << it->get_Cubit_name() << endl;
    if(it->get_cubit_bc_type_ulong() & BLOCKSET) {
      cout << "BLOCKSET" << endl;
      myfile << "BLOCKSET" << endl;
    }
    if(it->get_cubit_bc_type_ulong() & SIDESET) {
      cout << "SIDESET" << endl;
      myfile << "SIDESET" << endl;
    }
    if(it->get_cubit_bc_type_ulong() & NODESET) {
      cout << "NODESET" << endl;
      myfile << "NODESET" << endl;
    }
  }
  
  //Close mesh_file_name.txt
  myfile.close();
       
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  PetscFinalize();
 
}

