/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

    try {

    
  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
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
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;
    
    //Open mesh_file_name.txt for writing
    ofstream myfile;
    myfile.open ((string(mesh_file_name)+".txt").c_str());

  cout << "<<<< NodeSets >>>>>" << endl;
  //NodeSets
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,NodeSet,it)) {
    cout << *it << endl;
    ierr = it->print_Cubit_bc_data(cout); CHKERRQ(ierr);
    vector<char> bc_data;
    ierr = it->get_Cubit_bc_data(bc_data); CHKERRQ(ierr);
    if(bc_data.empty()) continue;
      
      //Displacement
      if (strcmp (&bc_data[0],"Displacement") == 0)
      {
          displacement_cubit_bc_data mydata;
          ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          cout << mydata;
          myfile << mydata;
      }
      
      //Force
      else if (strcmp (&bc_data[0],"Force") == 0)
      {
          force_cubit_bc_data mydata;
          ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          cout << mydata;
          myfile << mydata;
      }
      
      //Velocity
      else if (strcmp (&bc_data[0],"Velocity") == 0)
      {
          velocity_cubit_bc_data mydata;
          ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          cout << mydata;
          myfile << mydata;
      }
      
      //Acceleration
      else if (strcmp (&bc_data[0],"Acceleration") == 0)
      {
          acceleration_cubit_bc_data mydata;
          ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          cout << mydata;
          myfile << mydata;
      }
      
      //Temperature
      else if (strcmp (&bc_data[0],"Temperature") == 0)
      {
          temperature_cubit_bc_data mydata;
          ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          cout << mydata;
          myfile << mydata;
      }
      
      else SETERRQ(PETSC_COMM_SELF,1,"Error: Unrecognizable BC type");
      
  }

  cout << "<<<< SideSets >>>>>" << endl;
  //SideSets
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {
    cout << *it << endl;
    ierr = it->print_Cubit_bc_data(cout); CHKERRQ(ierr);
    vector<char> bc_data;
    ierr = it->get_Cubit_bc_data(bc_data); CHKERRQ(ierr);
    if(bc_data.empty()) continue;
      
      //Pressure
      if (strcmp (&bc_data[0],"Pressure") == 0)
      {
          pressure_cubit_bc_data mydata;
          ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          cout << mydata;
          myfile << mydata;
      }

      //Heat Flux
      else if (strcmp (&bc_data[0],"HeatFlux") == 0)
      {
          heatflux_cubit_bc_data mydata;
          ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          cout << mydata;
          myfile << mydata;
      }
      
      //cfd_bc
      else if (strcmp (&bc_data[0],"cfd_bc") == 0)
      {
          cfd_cubit_bc_data mydata;
          ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
          //Print data
          cout << mydata;
          myfile << mydata;
      }
      
      else SETERRQ(PETSC_COMM_SELF,1,"Error: Unrecognizable BC type");
  }

  cout << "<<<< BlockSets >>>>>" << endl;
  //BlockSets
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BlockSet,it))
  {
      cout << endl << *it << endl;

      //Get and print block name
      ierr = it->print_Cubit_name(cout); CHKERRQ(ierr);
      ierr = it->print_Cubit_name(myfile); CHKERRQ(ierr);
      
      
      //Get and print block attributes
      vector<double> attributes;
      ierr = it->get_Cubit_attributes(attributes); CHKERRQ(ierr);
      ierr = it->print_Cubit_attributes(cout); CHKERRQ(ierr);
      ierr = it->print_Cubit_attributes(myfile); CHKERRQ(ierr);
  }
  
        //Get block attributes and assign them as material properties/solution parameters based on the name of each block
        
        //Conventions:
        //----------------------------------------------------------------------------------------
        //Materials are defined with block names starting with MAT_ e.g. MAT_ELASTIC_abcd,
        //MAT_FRACTcdef etc.
        //Solution procedures are defined with block names starting with SOL_ e.g. SOL_ELASTIC_xx, SOL_NLELASTICxx, SOL_FRACTabcd etc.
        //----------------------------------------------------------------------------------------
        
        for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BlockSet,it))
            {
                cout << endl << *it << endl;
                
                //Get block name
                string name = it->get_Cubit_name();

                //Elastic material
                if (name.compare(0,11,"MAT_ELASTIC") == 0)
                {
                    Mat_Elastic mydata;
                    ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
                    //Print data
                    cout << mydata;
                    myfile << mydata;
                }
                else if (name.compare(0,12,"MAT_TRANSISO") == 0)
                {
                    Mat_TransIso mydata;
                    ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
                    //Print data
                    cout << mydata;
                    myfile << mydata;
                }

                else SETERRQ(PETSC_COMM_SELF,1,"Error: Unrecognizable Material type");
            }
        
  //Close mesh_file_name.txt
  myfile.close();

  PetscFinalize();
        
    } catch (const char* msg) {
        SETERRQ(PETSC_COMM_SELF,1,msg);
    }


}

