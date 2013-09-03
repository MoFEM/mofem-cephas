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

/*! \struct generic_cubit_bc_data
 *  \brief Generic bc data structure
 */
struct generic_cubit_bc_data {
    PetscErrorCode ierr;
    
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        SETERRQ(PETSC_COMM_SELF,1,"It makes no sense for the generic bc type");
        PetscFunctionReturn(0);
    }
    
};

/*! \struct displacement_cubit_bc_data
 *  \brief Definition of the displacement bc data structure
 */
struct displacement_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[12]; // 12 characters for "Displacement"
    char pre1; // Always zero
    char pre2; // pre-processing flags for modification of displacement bcs. They should not affect analysis, i.e. safe to ignore; 1: smallest combine, 2: average, 3: largest combine, 4: overwrite or no combination defined (default)
    char flag1; // Flag for X-Translation (0: N/A, 1: specified)
    char flag2; // Flag for Y-Translation (0: N/A, 1: specified)
    char flag3; // Flag for Z-Translation (0: N/A, 1: specified)
    char flag4; // Flag for X-Rotation (0: N/A, 1: specified)
    char flag5; // Flag for Y-Rotation (0: N/A, 1: specified)
    char flag6; // Flag for Z-Rotation (0: N/A, 1: specified)
    double value1; // Value for X-Translation
    double value2; // Value for Y-Translation
    double value3; // Value for Z-Translation
    double value4; // Value for X-Rotation
    double value5; // Value for Y-Rotation
    double value6; // Value for Z-Rotation
    };
    
    _data_ data;
    
        virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
            //Fill data
            memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
    
};

/*! \struct force_cubit_bc_data
 *  \brief Definition of the force bc data structure
 */
struct force_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[5]; // 5 characters for "Force"
    char zero[3]; // 3 zeros
    double value1; // Force magnitude
    double value2; // Moment magnitude
    double value3; // X-component of force direction vector
    double value4; // Y-component of force direction vector
    double value5; // Z-component of force direction vector
    double value6; // X-component of moment direction vector
    double value7; // Y-component of moment direction vector
    double value8; // Z-component of moment direction vector
    char zero2; // 0
    };
    
    _data_ data;
        virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
            //Fill data
            memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
    
};

/*! \struct velocity_cubit_bc_data
 *  \brief Definition of the velocity bc data structure
 */
struct velocity_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[8]; // 8 characters for "Velocity"
    char pre1; // Always zero
    char pre2; // pre-processing flags for modification of displacement bcs. They should not affect analysis, i.e. safe to ignore; 1: smallest combine, 2: average, 3: largest combine, 4: overwrite or no combination defined (default)
    char flag1; // Flag for X-Translation (0: N/A, 1: specified)
    char flag2; // Flag for Y-Translation (0: N/A, 1: specified)
    char flag3; // Flag for Z-Translation (0: N/A, 1: specified)
    char flag4; // Flag for X-Rotation (0: N/A, 1: specified)
    char flag5; // Flag for Y-Rotation (0: N/A, 1: specified)
    char flag6; // Flag for Z-Rotation (0: N/A, 1: specified)
    double value1; // Value for X-Translation
    double value2; // Value for Y-Translation
    double value3; // Value for Z-Translation
    double value4; // Value for X-Rotation
    double value5; // Value for Y-Rotation
    double value6; // Value for Z-Rotation
    };
    
    _data_ data;
    
        virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
            //Fill data
            memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
    
};

/*! \struct acceleration_cubit_bc_data
 *  \brief Definition of the acceleration bc data structure
 */
struct acceleration_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[12]; // 12 characters for "Acceleration"
    char pre1; // Always zero
    char pre2; // pre-processing flags for modification of displacement bcs. They should not affect analysis, i.e. safe to ignore; 1: smallest combine, 2: average, 3: largest combine, 4: overwrite or no combination defined (default)
    char flag1; // Flag for X-Translation (0: N/A, 1: specified)
    char flag2; // Flag for Y-Translation (0: N/A, 1: specified)
    char flag3; // Flag for Z-Translation (0: N/A, 1: specified)
    char flag4; // Flag for X-Rotation (0: N/A, 1: specified)
    char flag5; // Flag for Y-Rotation (0: N/A, 1: specified)
    char flag6; // Flag for Z-Rotation (0: N/A, 1: specified)
    double value1; // Value for X-Translation
    double value2; // Value for Y-Translation
    double value3; // Value for Z-Translation
    double value4; // Value for X-Rotation
    double value5; // Value for Y-Rotation
    double value6; // Value for Z-Rotation
    };
    
    _data_ data;
    
        virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
            //Fill data
            memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
    
};

/*! \struct temperature_cubit_bc_data
 *  \brief Definition of the temperature bc data structure
 */
struct temperature_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[11]; // 11 characters for "Temperature"
    char pre1; // This is always zero
    char pre2; // 0: temperature is not applied on thin shells (default); 1: temperature is applied on thin shells
    char flag1; // 0: N/A, 1: temperature value applied (not on thin shells)
    char flag2; // 0: N/A, 1: temperature applied on thin shell middle
    char flag3; // 0: N/A, 1: thin shell temperature gradient specified
    char flag4; // 0: N/A, 1: top thin shell temperature
    char flag5; // 0: N/A, 1: bottom thin shell temperature
    char flag6; // This is always zero
    double value1; // Temperature (default case - no thin shells)
    double value2; // Temperature for middle of thin shells
    double value3; // Temperature gradient for thin shells
    double value4; // Temperature for top of thin shells
    double value5; // Temperature for bottom of thin shells
    double value6; // This is always zero, i.e. ignore
    };
    
    _data_ data;

        virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
            //Fill data
            memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
    
};

struct pressure_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[8]; // 8 characters for "Pressure"
    char flag1; //
    char flag2; //
    double value1; // Pressure value
    char zero; //
    };
    
    _data_ data;
    
        virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
            memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
    
};

struct heatflux_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
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
    
    _data_ data;

        virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
            memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
    
};


PetscErrorCode func(const vector<char> &bc_data,generic_cubit_bc_data *ptr_cubit_bc_data) {
    PetscFunctionBegin;
    
    PetscErrorCode ierr;
    
    ierr = ptr_cubit_bc_data->fill_data(bc_data); CHKERRQ(ierr);
    
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
      
      //Displacement
      if (strcmp (&bc_data[0],"Displacement") == 0)
      {
          displacement_cubit_bc_data displacement_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&displacement_cubit_bc_struct); CHKERRQ(ierr);
                    
          //Print bc data
          printf(" \n");
          printf("BC type: ");
          for(int uu = 0;uu<12;uu++)
          {
              printf("%c ",displacement_cubit_bc_struct.data.name[uu]);
          }
          printf("\n \n");

          printf("Flag for X-Translation (0/1): %d\n",displacement_cubit_bc_struct.data.flag1);
          printf("Flag for Y-Translation (0/1): %d\n",displacement_cubit_bc_struct.data.flag2);
          printf("Flag for Z-Translation (0/1): %d\n",displacement_cubit_bc_struct.data.flag3);
          printf("Flag for X-Rotation (0/1): %d\n",displacement_cubit_bc_struct.data.flag4);
          printf("Flag for Y-Rotation (0/1): %d\n",displacement_cubit_bc_struct.data.flag5);
          printf("Flag for Z-Rotation (0/1): %d\n",displacement_cubit_bc_struct.data.flag6);
          printf("\n");
          
          
          if (displacement_cubit_bc_struct.data.flag1 == 1)
              printf("Displacement magnitude (X-Translation): %f\n",displacement_cubit_bc_struct.data.value1);
          else
              printf("Displacement magnitude (X-Translation): N/A \n");
          
          if (displacement_cubit_bc_struct.data.flag2 == 1)
              printf("Displacement magnitude (Y-Translation): %f\n",displacement_cubit_bc_struct.data.value2);
          else
              printf("Displacement magnitude (Y-Translation): N/A \n");
          
          if (displacement_cubit_bc_struct.data.flag3 == 1)
              printf("Displacement magnitude (Z-Translation): %f\n",displacement_cubit_bc_struct.data.value3);
          else
              printf("Displacement magnitude (Z-Translation): N/A \n");
          
          if (displacement_cubit_bc_struct.data.flag4 == 1)
              printf("Displacement magnitude (X-Rotation): %f\n",displacement_cubit_bc_struct.data.value4);
          else
              printf("Displacement magnitude (X-Rotation): N/A \n");
          
          if (displacement_cubit_bc_struct.data.flag5 == 1)
              printf("Displacement magnitude (Y-Rotation): %f\n",displacement_cubit_bc_struct.data.value5);
          else
              printf("Displacement magnitude (Y-Rotation): N/A \n");
          
          if (displacement_cubit_bc_struct.data.flag6 == 1)
              printf("Displacement magnitude (Z-Rotation): %f\n",displacement_cubit_bc_struct.data.value6);
          else
              printf("Displacement magnitude (Z-Rotation): N/A \n");
          
          printf("\n");
      }
      
      //Force
      else if (strcmp (&bc_data[0],"Force") == 0)
      {
          force_cubit_bc_data force_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&force_cubit_bc_struct); CHKERRQ(ierr);
          
          //Print bc data
          printf("\n");
          printf("BC type: ");
          for(int uu = 0;uu<5;uu++)
          {
              printf("%c ",force_cubit_bc_struct.data.name[uu]);
          }
          printf("\n \n");
          printf("Force magnitude: %f\n",force_cubit_bc_struct.data.value1);
          printf("Moment magnitude: %f\n",force_cubit_bc_struct.data.value2);
          printf("Force direction vector (X-component): %f\n",force_cubit_bc_struct.data.value3);
          printf("Force direction vector (Y-component): %f\n",force_cubit_bc_struct.data.value4);
          printf("Force direction vector (Z-component): %f\n",force_cubit_bc_struct.data.value5);
          printf("Moment direction vector (X-component): %f\n",force_cubit_bc_struct.data.value6);
          printf("Moment direction vector (Y-component): %f\n",force_cubit_bc_struct.data.value7);
          printf("Moment direction vector (Z-component): %f\n",force_cubit_bc_struct.data.value8);
          printf("\n");
      }
      
      //Velocity
      else if (strcmp (&bc_data[0],"Velocity") == 0)
      {
          velocity_cubit_bc_data velocity_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&velocity_cubit_bc_struct); CHKERRQ(ierr);
          
          //Print bc data
          printf("\n");
          printf("BC type: ");
          for(int uu = 0;uu<8;uu++)
          {
              printf("%c ",velocity_cubit_bc_struct.data.name[uu]);
          }
          printf("\n \n");
          
          if (velocity_cubit_bc_struct.data.flag1 == 1)
              printf("Velocity magnitude (X-Translation): %f\n",velocity_cubit_bc_struct.data.value1);
          else
              printf("Velocity magnitude (X-Translation): N/A \n");

          if (velocity_cubit_bc_struct.data.flag2 == 1)
              printf("Velocity magnitude (Y-Translation): %f\n",velocity_cubit_bc_struct.data.value2);
          else
              printf("Velocity magnitude (Y-Translation): N/A \n");
          
          if (velocity_cubit_bc_struct.data.flag3 == 1)
              printf("Velocity magnitude (Z-Translation): %f\n",velocity_cubit_bc_struct.data.value3);
          else
              printf("Velocity magnitude (Z-Translation): N/A \n");
          
          if (velocity_cubit_bc_struct.data.flag4 == 1)
              printf("Velocity magnitude (X-Rotation): %f\n",velocity_cubit_bc_struct.data.value4);
          else
              printf("Velocity magnitude (X-Rotation): N/A \n");
          
          if (velocity_cubit_bc_struct.data.flag5 == 1)
              printf("Velocity magnitude (Y-Rotation): %f\n",velocity_cubit_bc_struct.data.value5);
          else
              printf("Velocity magnitude (Y-Rotation): N/A \n");
          
          if (velocity_cubit_bc_struct.data.flag6 == 1)
              printf("Velocity magnitude (Z-Rotation): %f\n",velocity_cubit_bc_struct.data.value6);
          else
              printf("Velocity magnitude (Z-Rotation): N/A \n");
          
          printf("\n");
      }
      
      //Acceleration
      else if (strcmp (&bc_data[0],"Acceleration") == 0)
      {
          acceleration_cubit_bc_data acceleration_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&acceleration_cubit_bc_struct); CHKERRQ(ierr);
          
          //Print bc data
          printf("\n");
          printf("BC type: ");
          for(int uu = 0;uu<12;uu++)
          {
              printf("%c ",acceleration_cubit_bc_struct.data.name[uu]);
          }
          printf("\n \n");
          
          if (acceleration_cubit_bc_struct.data.flag1 == 1)
              printf("Acceleration magnitude (X-Translation): %f\n",acceleration_cubit_bc_struct.data.value1);
          else
              printf("Acceleration magnitude (X-Translation): N/A \n");
          
          if (acceleration_cubit_bc_struct.data.flag2 == 1)
              printf("Acceleration magnitude (Y-Translation): %f\n",acceleration_cubit_bc_struct.data.value2);
          else
              printf("Acceleration magnitude (Y-Translation): N/A \n");
          
          if (acceleration_cubit_bc_struct.data.flag3 == 1)
              printf("Acceleration magnitude (Z-Translation): %f\n",acceleration_cubit_bc_struct.data.value3);
          else
              printf("Acceleration magnitude (Z-Translation): N/A \n");
          
          if (acceleration_cubit_bc_struct.data.flag4 == 1)
              printf("Acceleration magnitude (X-Rotation): %f\n",acceleration_cubit_bc_struct.data.value4);
          else
              printf("Acceleration magnitude (X-Rotation): N/A \n");
          
          if (acceleration_cubit_bc_struct.data.flag5 == 1)
              printf("Acceleration magnitude (Y-Rotation): %f\n",acceleration_cubit_bc_struct.data.value5);
          else
              printf("Acceleration magnitude (Y-Rotation): N/A \n");
          
          if (acceleration_cubit_bc_struct.data.flag6 == 1)
              printf("Acceleration magnitude (Z-Rotation): %f\n",acceleration_cubit_bc_struct.data.value6);
          else
              printf("Acceleration magnitude (Z-Rotation): N/A \n");
          
          printf("\n");
      }
      
      //Temperature
      else if (strcmp (&bc_data[0],"Temperature") == 0)
      {
          temperature_cubit_bc_data temperature_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&temperature_cubit_bc_struct); CHKERRQ(ierr);
          
          //Print bc data
          printf("\n");
          printf("BC type: ");
          for(int uu = 0;uu<11;uu++)
          {
              printf("%c ",temperature_cubit_bc_struct.data.name[uu]);
          }
          printf("\n \n");
          
          if (temperature_cubit_bc_struct.data.flag1 == 1)
              printf("Temperature: %f\n",temperature_cubit_bc_struct.data.value1);
          else
              printf("Temperature (default case): N/A \n");
          
          if (temperature_cubit_bc_struct.data.flag2 == 1)
              printf("Temperature (thin shell middle): %f\n",temperature_cubit_bc_struct.data.value2);
          else
              printf("Temperature (thin shell middle): N/A \n");

          if (temperature_cubit_bc_struct.data.flag3 == 1)
              printf("Temperature (thin shell gradient): %f\n",temperature_cubit_bc_struct.data.value3);
          else
              printf("Temperature (thin shell gradient): N/A \n");
          
          if (temperature_cubit_bc_struct.data.flag4 == 1)
              printf("Temperature (thin shell top): %f\n",temperature_cubit_bc_struct.data.value4);
          else
              printf("Temperature (thin shell top): N/A \n");
          
          if (temperature_cubit_bc_struct.data.flag5 == 1)
              printf("Temperature (thin shell bottom): %f\n",temperature_cubit_bc_struct.data.value5);
          else
              printf("Temperature (thin shell bottom): N/A \n");
 
          printf("\n");          
      }
      
      else SETERRQ(PETSC_COMM_SELF,1,"Error: Unrecognizable BC type");
      
  }

  cout << "<<<< SideSets >>>>>" << endl;
  //SideSets
  for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,SideSet,it)) {
    cout << *it << endl;
    ierr = it->print_Cubit_bc_data(cout); CHKERRQ(ierr);
    vector<char> bc_data;
    ierr = it->get_Cubit_bc_data(bc_data); CHKERRQ(ierr);
      
      //Pressure
      if (strcmp (&bc_data[0],"Pressure") == 0)
      {
          pressure_cubit_bc_data pressure_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&pressure_cubit_bc_struct); CHKERRQ(ierr);
          
          //Print bc data
          printf("\n");
          printf("BC type: ");
          for(int uu = 0;uu<8;uu++)
          {
              printf("%c ",pressure_cubit_bc_struct.data.name[uu]);
          }
          printf("\n");                    
          printf("Pressure value: %f\n",pressure_cubit_bc_struct.data.value1);
          printf("\n");
      }

      //Heat Flux
      else if (strcmp (&bc_data[0],"HeatFlux") == 0)
      {
          heatflux_cubit_bc_data heatflux_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&heatflux_cubit_bc_struct); CHKERRQ(ierr);
          
          //Print bc data
          printf("\n");
          printf("BC type: ");
          for(int uu = 0;uu<8;uu++)
          {
              printf("%c ",heatflux_cubit_bc_struct.data.name[uu]);
          }
          printf("\n");
          printf("Heat flux applied on thin shells (Yes=1, No=0): %c\n",heatflux_cubit_bc_struct.data.flag2);
          printf("Heat flux value (no thin shells): %f\n",heatflux_cubit_bc_struct.data.value1);
          printf("Heat flux value (top of thin shells): %f\n",heatflux_cubit_bc_struct.data.value2);
          printf("Heat flux value (bottom of thin shells): %f\n",heatflux_cubit_bc_struct.data.value3);
          printf("\n");
      }
                 
      else SETERRQ(PETSC_COMM_SELF,1,"Error: Unrecognizable BC type");
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

