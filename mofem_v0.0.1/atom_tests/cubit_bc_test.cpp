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


// Temperature (NodeSets)
struct __attribute__((packed)) bc_temp {  // pack the structure in order to get rid of padding
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
};

// Velocity (NodeSets)
struct __attribute__((packed)) bc_veloc {  // pack the structure in order to get rid of padding
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
};

// Acceleration (NodeSets)
struct __attribute__((packed)) bc_accel {  // pack the structure in order to get rid of padding
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
};

// Force (NodeSets)
struct __attribute__((packed)) bc_force {  // pack the structure in order to get rid of padding
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

// Displacement (NodeSets)
struct __attribute__((packed)) bc_disp {  // need to pack the structure in order to get rid of padding
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
};

// Pressure (SideSets)
struct __attribute__((packed)) bc_press {  // pack the structure in order to get rid of padding
    char name[8]; // 8 characters for "Pressure"
    char flag1; //
    char flag2; //    
    double value1; // Pressure value
    char zero; //
    
};

// Heat Flux (SideSets)
struct __attribute__((packed)) bc_heatflux {  // need to pack the structure in order to get rid of padding
    char name[8]; // 8 characters for "HeatFlux" (no space)
    char flag1; //
    char flag2; //
    char flag3; //
    char flag4; //
    char flag5; //

    double value1; //
    double value2; //
    double value3; //

};


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

