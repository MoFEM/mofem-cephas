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

#include "moabField.hpp"
#include "moabField_Core.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

/*! \struct xxxxxxxxx
 *  \brief xxxxxxxxxx
 *
 *   Details ...
 *
 */


// generic bc data structure
struct generic_cubit_bc_data {
    PetscErrorCode ierr;
    
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        SETERRQ(PETSC_COMM_SELF,1,"It makes no sense for the generic bc type");
        PetscFunctionReturn(0);
    }
    
};

struct displacement_cubit_bc_data: public generic_cubit_bc_data {
    char name[12]; // 12 characters for "Displacement"
    char zero;
    char four;
    char flag1; // Flag for X-Translation (0: not used, 1: applied)
    char flag2; // Flag for Y-Translation (0: not used, 1: applied)
    char flag3; // Flag for Z-Translation (0: not used, 1: applied)
    char flag4; // Flag for X-Rotation (0: not used, 1: applied)
    char flag5; // Flag for Y-Rotation (0: not used, 1: applied)
    char flag6; // Flag for Z-Rotation (0: not used, 1: applied)
    double value1; // Value of X-Translation
    double value2; // Value of Y-Translation
    double value3; // Value of Z-Translation
    double value4; // Value of X-Rotation
    double value5; // Value of Y-Rotation
    double value6; // Value of Z-Rotation
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
};

struct force_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) data{
    char name[5]; // 5 characters for "Force"
    char zero[3]; //
    double value1; // Value of Force
    double value2; // Value of Moment
    double value3; // X-component of force direction vector
    double value4; // Y-component of force direction vector
    double value5; // Z-component of force direction vector
    double value6; // X-component of moment direction vector
    double value7; // Y-component of moment direction vector
    double value8; // Z-component of moment direction vector
    char zero2; // 0
    };
    
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
};

struct velocity_cubit_bc_data: public generic_cubit_bc_data {
    char name[8]; // 8 characters for "Velocity"
    char zero[2];
    char flag1; // Flag for X-Translation (0: not used, 1: applied)
    char flag2; // Flag for Y-Translation (0: not used, 1: applied)
    char flag3; // Flag for Z-Translation (0: not used, 1: applied)
    char flag4; // Flag for X-Rotation (0: not used, 1: applied)
    char flag5; // Flag for Y-Rotation (0: not used, 1: applied)
    char flag6; // Flag for Z-Rotation (0: not used, 1: applied)
    double value1; // Value of X-Translation
    double value2; // Value of Y-Translation
    double value3; // Value of Z-Translation
    double value4; // Value of X-Rotation
    double value5; // Value of Y-Rotation
    double value6; // Value of Z-Rotation
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
};

struct acceleration_cubit_bc_data: public generic_cubit_bc_data {
    char name[12]; // 12 characters for "Acceleration"
    char zero[2]; //
    char flag1; // Flag for X-Translation (0: not used, 1: applied)
    char flag2; // Flag for Y-Translation (0: not used, 1: applied)
    char flag3; // Flag for Z-Translation (0: not used, 1: applied)
    char flag4; // Flag for X-Rotation (0: not used, 1: applied)
    char flag5; // Flag for Y-Rotation (0: not used, 1: applied)
    char flag6; // Flag for Z-Rotation (0: not used, 1: applied)
    double value1; // Value of X-Translation
    double value2; // Value of Y-Translation
    double value3; // Value of Z-Translation
    double value4; // Value of X-Rotation
    double value5; // Value of Y-Rotation
    double value6; // Value of Z-Rotation
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
};

struct temperature_cubit_bc_data: public generic_cubit_bc_data {
    char name[11]; // 11 characters for "Temperature"
    char zero[2]; //
    char flag1; // Flag for X-Translation (0: not used, 1: applied)
    char flag2; // Flag for Y-Translation (0: not used, 1: applied)
    char flag3; // Flag for Z-Translation (0: not used, 1: applied)
    char flag4; // Flag for X-Rotation (0: not used, 1: applied)
    char flag5; // Flag for Y-Rotation (0: not used, 1: applied)
    char flag6; // Flag for Z-Rotation (0: not used, 1: applied)
    double value1; // Value of X-Translation
    double value2; // Value of Y-Translation
    double value3; // Value of Z-Translation
    double value4; // Value of X-Rotation
    double value5; // Value of Y-Rotation
    double value6; // Value of Z-Rotation
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
};

struct pressure_cubit_bc_data: public generic_cubit_bc_data {
    char name[8]; // 8 characters for "Pressure"
    char flag1; //
    char flag2; //
    double value1; // Pressure value
    char zero; //
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
};

struct heatflux_cubit_bc_data: public generic_cubit_bc_data {
    char name[8]; // 8 characters for "HeatFlux" (no space)
    char flag1; //
    char flag2; //
    char flag3; //
    char flag4; //
    char flag5; //
    double value1; //
    double value2; //
    double value3; //
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
};


PetscErrorCode func(const vector<char> &bc_data,generic_cubit_bc_data *ptr_cubit_bc_data) {
    PetscFunctionBegin;
    
    PetscErrorCode ierr;
    
    int scom;
    if (strcmp (&bc_data[0],"Displacement") == 0)
        scom = 1;
    else if (strcmp (&bc_data[0],"Force") == 0)
        scom = 2;
    else if (strcmp (&bc_data[0],"Velocity") == 0)
        scom = 3;
    else if (strcmp (&bc_data[0],"Acceleration") == 0)
        scom = 4;
    else if (strcmp (&bc_data[0],"Temperature") == 0)
        scom = 5;
    else if (strcmp (&bc_data[0],"Pressure") == 0)
        scom = 6;
    else if (strcmp (&bc_data[0],"HeatFlux") == 0)
        scom = 7;
    else SETERRQ(PETSC_COMM_SELF,1,"Error: Unrecognizable BC type");
    
    switch (scom)
    {
        case 1: {
            //Displacement
            //FILL displacement_cubit_bc_data
            cout << "<Displacement>";
            ierr = ptr_cubit_bc_data->fill_data(bc_data); CHKERRQ(ierr);
            
        }
            break;
        case 2: {
            //Force
            //FILL force_cubit_bc_data
            ierr = ptr_cubit_bc_data->fill_data(bc_data); CHKERRQ(ierr);
            
        }
            break;
        case 3: {
            //Velocity
            //FILL velocity_cubit_bc_data
            cout << "Velocity" << endl;
            ierr = ptr_cubit_bc_data->fill_data(bc_data); CHKERRQ(ierr);
            
        }
            break;
        case 4: {
            //Acceleration
            //FILL acceleration_cubit_bc_data
            cout << "Acceleration" << endl;
            ierr = ptr_cubit_bc_data->fill_data(bc_data); CHKERRQ(ierr);
        }
            break;
        case 5: {
            //Temperature
            //FILL temperature_cubit_bc_data
            cout << "Temperature" << endl;
            ierr = ptr_cubit_bc_data->fill_data(bc_data); CHKERRQ(ierr);
        }
            break;
        case 6: {
            //Pressure
            //FILL pressure_cubit_bc_data
            cout << "Pressure" << endl;
            ierr = ptr_cubit_bc_data->fill_data(bc_data); CHKERRQ(ierr);
        }
            break;
        case 7: {
            //Heat Flux
            //FILL heatflux_cubit_bc_data
            cout << "Heat Flux" << endl;
            ierr = ptr_cubit_bc_data->fill_data(bc_data); CHKERRQ(ierr);
        }
            break;
    
        default:
            SETERRQ(PETSC_COMM_SELF,1,"Something went wrong filling cubit_bc_data");
    }
    
    PetscFunctionReturn(0);
}


int main(int argc, char *argv[]) {

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

  //We need that for code profiling
  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  //Create MoFEM (Joseph) database
  moabField_Core core(moab);
  moabField& mField = core;

  cout << "<<<< NodeSets >>>>>" << endl;
  //NodeSets
  for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,NodeSet,it)) {
    cout << *it << endl;
    ierr = it->print_Cubit_bc_data(cout); CHKERRQ(ierr);
    vector<char> bc_data;
    ierr = it->get_Cubit_bc_data(bc_data); CHKERRQ(ierr);
      
      
      //displacement_cubit_bc_data disp_cubit_bc_struct;
      force_cubit_bc_data force_cubit_bc_struct;
      
      ierr = func(bc_data,&force_cubit_bc_struct); CHKERRQ(ierr);
      
      
      //testing
      int z;
      z=sizeof(force_cubit_bc_struct); cout << z;
      
      //memcpy(&force_cubit_bc_struct,&bc_data[0],sizeof(z));


       //printf("---> value1 %f\n",force_cubit_bc_struct.value1);

      
  }

  cout << "<<<< SideSets >>>>>" << endl;
  //SideSets
  for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,SideSet,it)) {
    cout << *it << endl;
    ierr = it->print_Cubit_bc_data(cout); CHKERRQ(ierr);
    vector<char> bc_data;
    ierr = it->get_Cubit_bc_data(bc_data); CHKERRQ(ierr);
  } 

  cout << "<<<< BlockSets >>>>>" << endl;
  //BlockSets
  for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,BlockSet,it)) {
    cout << *it << endl;
    ierr = it->print_Cubit_bc_data(cout); CHKERRQ(ierr);
    vector<char> bc_data;
    ierr = it->get_Cubit_bc_data(bc_data); CHKERRQ(ierr);
  } 


  PetscFinalize();

}

