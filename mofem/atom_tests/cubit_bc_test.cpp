/** \file cubit_bc_test.cpp
  \brief Atom test for getting boundary conditions from cubit

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
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetString(PETSC_NULL,"","-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #else
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #endif
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  //Open mesh_file_name.txt for writing
  std::ofstream myfile;
  myfile.open ((std::string(mesh_file_name)+".txt").c_str());

  std::cout << "<<<< NODESETs >>>>>" << std::endl;
  //NODESETs
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,NODESET,it)) {
    std::cout << *it << std::endl;
    ierr = it->print_bc_data(std::cout); CHKERRQ(ierr);
    std::vector<char> bc_data;
    ierr = it->get_bc_data(bc_data); CHKERRQ(ierr);
    if(bc_data.empty()) continue;

      //Displacement
      if (strcmp (&bc_data[0],"Displacement") == 0)
      {
          DisplacementCubitBcData mydata;
          ierr = it->get_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          std::cout << mydata;
          myfile << mydata;
      }

      //Force
      else if (strcmp (&bc_data[0],"Force") == 0)
      {
          ForceCubitBcData mydata;
          ierr = it->get_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          std::cout << mydata;
          myfile << mydata;
      }

      //Velocity
      else if (strcmp (&bc_data[0],"Velocity") == 0)
      {
          VelocityCubitBcData mydata;
          ierr = it->get_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          std::cout << mydata;
          myfile << mydata;
      }

      //Acceleration
      else if (strcmp (&bc_data[0],"Acceleration") == 0)
      {
          AccelerationCubitBcData mydata;
          ierr = it->get_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          std::cout << mydata;
          myfile << mydata;
      }

      //Temperature
      else if (strcmp (&bc_data[0],"Temperature") == 0)
      {
          TemperatureCubitBcData mydata;
          ierr = it->get_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          std::cout << mydata;
          myfile << mydata;
      }

      else SETERRQ(PETSC_COMM_SELF,1,"Error: Unrecognizable BC type");

  }

  std::cout << "<<<< SIDESETs >>>>>" << std::endl;
  //SIDESETs
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    std::cout << *it << std::endl;
    ierr = it->print_bc_data(std::cout); CHKERRQ(ierr);
    std::vector<char> bc_data;
    ierr = it->get_bc_data(bc_data); CHKERRQ(ierr);
    if(bc_data.empty()) continue;

      //Pressure
      if (strcmp (&bc_data[0],"Pressure") == 0)
      {
          PressureCubitBcData mydata;
          ierr = it->get_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          std::cout << mydata;
          myfile << mydata;
      }

      //Heat Flux
      else if (strcmp (&bc_data[0],"HeatFlux") == 0)
      {
          HeatFluxCubitBcData mydata;
          ierr = it->get_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          std::cout << mydata;
          myfile << mydata;
      }

      //cfd_bc
      else if (strcmp (&bc_data[0],"cfd_bc") == 0)
      {
          CfgCubitBcData mydata;
          ierr = it->get_bc_data_structure(mydata); CHKERRQ(ierr);

          //Interface bc (Hex:6 Dec:6)
          if (mydata.data.type == 6) {  // 6 is the decimal value of the corresponding value (hex) in bc_data
              //Print data
              std::cout << std::endl << "Interface" << std::endl;
              myfile << std::endl << "Interface" << std::endl;
              std::cout << mydata;
              myfile << mydata;
          }
          //Pressure inlet (Hex:f Dec:15)
          else if (mydata.data.type == 15) {  // 15 is the decimal value of the corresponding value (hex) in bc_data
              //Print data
              std::cout << std::endl << "Pressure Inlet" << std::endl;
              myfile << std::endl << "Pressure Inlet" << std::endl;
              std::cout << mydata;
              myfile << mydata;
          }
          //Pressure outlet (Hex:10 Dec:16)
          else if (mydata.data.type == 16) {  // 16 is the decimal value of the corresponding value (hex) in bc_data
              //Print data
              std::cout << std::endl << "Pressure Outlet" << std::endl;
              myfile << std::endl << "Pressure Outlet" << std::endl;
              std::cout << mydata;
              myfile << mydata;
          }
      }

      else SETERRQ(PETSC_COMM_SELF,1,"Error: Unrecognizable BC type");
  }


  std::cout << "<<<< BLOCKSETs >>>>>" << std::endl;
  //BLOCKSETs
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it))
  {
      std::cout << std::endl << *it << std::endl;

      //Get and print block name
      ierr = it->print_name(std::cout); CHKERRQ(ierr);
      ierr = it->print_name(myfile); CHKERRQ(ierr);


      //Get and print block attributes
      std::vector<double> attributes;
      ierr = it->get_attributes(attributes); CHKERRQ(ierr);
      ierr = it->print_attributes(std::cout); CHKERRQ(ierr);
      ierr = it->print_attributes(myfile); CHKERRQ(ierr);
  }

  //Get block attributes and assign them as material properties/solution parameters based on the name of each block

  //Conventions:
  //----------------------------------------------------------------------------------------
  //Materials are defined with block names starting with MAT_ e.g. MAT_ELASTIC_abcd,
  //MAT_FRACTcdef etc.
  //Solution procedures are defined with block names starting with SOL_ e.g. SOL_ELASTIC_xx, SOL_NLELASTICxx, SOL_FRACTabcd etc.
  //----------------------------------------------------------------------------------------

  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {

    std::cout << std::endl << *it << std::endl;

    //Get block name
    std::string name = it->get_name();

    //Elastic material
    if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
      Mat_Elastic_TransIso mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      //Print data
      std::cout << mydata;
      myfile << mydata;
    } else if (name.compare(0,11,"MAT_ELASTIC") == 0) {
      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      //Print data
      std::cout << mydata;
      myfile << mydata;
    } else if (name.compare(0,10,"MAT_INTERF") == 0) {
      Mat_Interf mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      //Print data
      std::cout << mydata;
      myfile << mydata;
    } else SETERRQ(PETSC_COMM_SELF,1,"Error: Unrecognizable Material type");

  }

  bool add_block_is_there = false;
  ierr = m_field.add_cubit_msId(BLOCKSET,1000,"ADD_BLOCK_SET"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
    //Get block name
    std::string name = it->get_name();
    if (name.compare(0,13,"ADD_BLOCK_SET") == 0) {
      add_block_is_there = true;
      std::vector<double> attributes;
      it->get_attributes(attributes);
    }
  }

  if(!add_block_is_there) {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_OPERATION_UNSUCCESSFUL,"no added block set");
  }

  //Close mesh_file_name.txt
  myfile.close();


  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFinalize();

}
