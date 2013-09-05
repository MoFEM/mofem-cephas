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
    
    friend ostream& operator<<(ostream& os,const displacement_cubit_bc_data& e);
    
};

/*! \brief Print displacement bc data
 */
ostream& operator<<(ostream& os,const displacement_cubit_bc_data& e)
{
    os << "\n";
    os << "D i s p l a c e m e n t \n \n";
    os << "Flag for X-Translation (0/1): " << (int)e.data.flag1 << "\n";
    os << "Flag for Y-Translation (0/1): " << (int)e.data.flag2 << "\n";
    os << "Flag for Z-Translation (0/1): " << (int)e.data.flag3 << "\n";
    os << "Flag for X-Rotation (0/1): " << (int)e.data.flag4 << "\n";
    os << "Flag for Y-Rotation (0/1): " << (int)e.data.flag5 << "\n";
    os << "Flag for Z-Rotation (0/1): " << (int)e.data.flag6 << "\n \n";
    
    if (e.data.flag1 == 1)
        os << "Displacement magnitude (X-Translation): " << e.data.value1 << "\n";
    else os << "Displacement magnitude (X-Translation): N/A" << "\n";
    if (e.data.flag2 == 1)
        os << "Displacement magnitude (Y-Translation): " << e.data.value2 << "\n";
    else os << "Displacement magnitude (Y-Translation): N/A" << "\n";
    if (e.data.flag3 == 1)
        os << "Displacement magnitude (Z-Translation): " << e.data.value3 << "\n";
    else os << "Displacement magnitude (Z-Translation): N/A" << "\n";
    if (e.data.flag4 == 1)
        os << "Displacement magnitude (X-Rotation): " << e.data.value4 << "\n";
    else os << "Displacement magnitude (X-Rotation): N/A" << "\n";
    if (e.data.flag5 == 1)
        os << "Displacement magnitude (Y-Rotation): " << e.data.value5 << "\n";
    else os << "Displacement magnitude (Y-Rotation): N/A" << "\n";
    if (e.data.flag6 == 1)
        os << "Displacement magnitude (Z-Rotation): " << e.data.value6 << "\n";
    else os << "Displacement magnitude (Z-Rotation): N/A" << "\n";
}


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
    
    friend ostream& operator<<(ostream& os,const force_cubit_bc_data& e);
    
};

/*! \brief Print force bc data
 */
ostream& operator<<(ostream& os,const force_cubit_bc_data& e)
{
    os << "\n";
    os << "F o r c e \n \n";
    os << "Force magnitude: " << e.data.value1 << "\n";
    os << "Moment magnitude: " << e.data.value2 << "\n";
    os << "Force direction vector (X-component): " << e.data.value3 << "\n";
    os << "Force direction vector (Y-component): " << e.data.value4 << "\n";
    os << "Force direction vector (Z-component): " << e.data.value5 << "\n";
    os << "Moment direction vector (X-component): " << e.data.value6 << "\n";
    os << "Moment direction vector (Y-component): " << e.data.value7 << "\n";
    os << "Moment direction vector (Z-component): " << e.data.value8 << "\n";
}


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
    
    friend ostream& operator<<(ostream& os,const velocity_cubit_bc_data& e);
    
};

/*! \brief Print velocity bc data
 */
ostream& operator<<(ostream& os,const velocity_cubit_bc_data& e)
{
    os << "\n";
    os << "V e l o c i t y \n \n";
    if (e.data.flag1 == 1)
        os << "Velocity magnitude (X-Translation): " << e.data.value1 << "\n";
    else os << "Velocity magnitude (X-Translation): N/A" << "\n";
    if (e.data.flag2 == 1)
        os << "Velocity magnitude (Y-Translation): " << e.data.value2 << "\n";
    else os << "Velocity magnitude (Y-Translation): N/A" << "\n";
    if (e.data.flag3 == 1)
        os << "Velocity magnitude (Z-Translation): " << e.data.value3 << "\n";
    else os << "Velocity magnitude (Z-Translation): N/A" << "\n";
    if (e.data.flag4 == 1)
        os << "Velocity magnitude (X-Rotation): " << e.data.value4 << "\n";
    else os << "Velocity magnitude (X-Rotation): N/A" << "\n";
    if (e.data.flag5 == 1)
        os << "Velocity magnitude (Y-Rotation): " << e.data.value5 << "\n";
    else os << "Velocity magnitude (Y-Rotation): N/A" << "\n";
    if (e.data.flag6 == 1)
        os << "Velocity magnitude (Z-Rotation): " << e.data.value6 << "\n";
    else os << "Velocity magnitude (Z-Rotation): N/A" << "\n";
}
    

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
    
    friend ostream& operator<<(ostream& os,const acceleration_cubit_bc_data& e);
    
};

/*! \brief Print acceleration bc data
 */
ostream& operator<<(ostream& os,const acceleration_cubit_bc_data& e)
{
    os << "\n";
    os << "A c c e l e r a t i o n \n \n";
    if (e.data.flag1 == 1)
        os << "Acceleration magnitude (X-Translation): " << e.data.value1 << "\n";
    else os << "Acceleration magnitude (X-Translation): N/A" << "\n";
    if (e.data.flag2 == 1)
        os << "Acceleration magnitude (Y-Translation): " << e.data.value2 << "\n";
    else os << "Acceleration magnitude (Y-Translation): N/A" << "\n";
    if (e.data.flag3 == 1)
        os << "Acceleration magnitude (Z-Translation): " << e.data.value3 << "\n";
    else os << "Acceleration magnitude (Z-Translation): N/A" << "\n";
    if (e.data.flag4 == 1)
        os << "Acceleration magnitude (X-Rotation): " << e.data.value4 << "\n";
    else os << "Acceleration magnitude (X-Rotation): N/A" << "\n";
    if (e.data.flag5 == 1)
        os << "Acceleration magnitude (Y-Rotation): " << e.data.value5 << "\n";
    else os << "Acceleration magnitude (Y-Rotation): N/A" << "\n";
    if (e.data.flag6 == 1)
        os << "Acceleration magnitude (Z-Rotation): " << e.data.value6 << "\n";
    else os << "Acceleration magnitude (Z-Rotation): N/A" << "\n";
}



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
    
    friend ostream& operator<<(ostream& os,const temperature_cubit_bc_data& e);
};

/*! \brief Print temperature bc data
 */
ostream& operator<<(ostream& os,const temperature_cubit_bc_data& e)
{
    os << "\n";
    os << "T e m p e r a t u r e \n \n";
    if (e.data.flag1 == 1)
        os << "Temperature: " << e.data.value1 << "\n";
    else os << "Temperature (default case): N/A" << "\n";
    if (e.data.flag2 == 1)
        os << "Temperature (thin shell middle): " << e.data.value2 << "\n";
    else os << "Temperature (thin shell middle): N/A" << "\n";
    if (e.data.flag3 == 1)
        os << "Temperature (thin shell gradient): " << e.data.value3 << "\n";
    else os << "Temperature (thin shell gradient): N/A" << "\n";
    if (e.data.flag4 == 1)
        os << "Temperature (thin shell top): " << e.data.value4 << "\n";
    else os << "Temperature (thin shell top): N/A" << "\n";
    if (e.data.flag5 == 1)
        os << "Temperature (thin shell bottom): " << e.data.value5 << "\n";
    else os << "Temperature (thin shell bottom): N/A" << "\n";
}


/*! \struct pressure_cubit_bc_data
 *  \brief Definition of the pressure bc data structure
 */
struct pressure_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[8]; // 8 characters for "Pressure"
    char flag1; // This is always zero
    char flag2; // 0: Pressure is interpeted as pure pressure 1: pressure is interpreted as total force
    double value1; // Pressure value
    char zero; // This is always zero
    };
    
    _data_ data;
    
        virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
            //Fill data
            memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
    
    friend ostream& operator<<(ostream& os,const pressure_cubit_bc_data& e);
    
};

/*! \brief Print pressure bc data
 */
ostream& operator<<(ostream& os,const pressure_cubit_bc_data& e)
{
    os << "\n";
    os << "P r e s s u r e \n \n";
    os << "Pressure value: " << e.data.value1 << "\n";
}

/*! \struct heatflux_cubit_bc_data
 *  \brief Definition of the heat flux bc data structure
 */
struct heatflux_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[8]; // 8 characters for "HeatFlux" (no space)
    char pre1; // This is always zero
    char pre2; // 0: heat flux is not applied on thin shells (default); 1: heat flux is applied on thin shells
    char flag1; // 0: N/A, 1: normal heat flux case (i.e. single value, case without thin shells)
    char flag2; // 0: N/A, 1: Thin shell top heat flux specified
    char flag3; // 0: N/A, 1: Thin shell bottom heat flux specidied
    double value1; // Heat flux value for default case (no thin shells)
    double value2; // Heat flux (thin shell top)
    double value3; // Heat flux (thin shell bottom)
    };
    
    _data_ data;

        virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
            //Fill data
            memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
    
    friend ostream& operator<<(ostream& os,const heatflux_cubit_bc_data& e);
    
};

/*! \brief Print heat flux bc data
 */
ostream& operator<<(ostream& os,const heatflux_cubit_bc_data& e)
{
    os << "\n";
    os << "H e a t  F l u x \n \n";
    if (e.data.flag1 == 1)
        os << "Heat flux value: " << e.data.value1 << "\n";
    else os << "Heat flux is applied on thin shells" << "\n";
    if (e.data.flag2 == 1)
        os << "Heat flux value (thin shell top): " << e.data.value2 << "\n";
    else os << "Heat flux value (thin shell top): N/A" << "\n";
    if (e.data.flag3 == 1)
        os << "Heat flux value (thin shell bottom): " << e.data.value3 << "\n";
    else os << "Heat flux value (thin shell bottom): N/A" << "\n";
}

/*! \fn func
 *  \brief Function that fills the input struct with bc_data
 */
PetscErrorCode func(const vector<char> &bc_data,generic_cubit_bc_data *ptr_cubit_bc_data) {
    PetscFunctionBegin;
    
    PetscErrorCode ierr;
    
    ierr = ptr_cubit_bc_data->fill_data(bc_data); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

/*! \fn get_type_from_bc_data
 *  \brief Function that returns the Cubit_BC_bitset type of the contents of bc_data
 */
PetscErrorCode get_type_from_bc_data(const vector<char> &bc_data,int &type)
{
    PetscFunctionBegin;
    
    //See Cubit_BC_bitset in common.hpp
    
    if (strcmp (&bc_data[0],"Displacement") == 0)
        type = 4;
    else if (strcmp (&bc_data[0],"Force") == 0)
        type = 5;
    else if (strcmp (&bc_data[0],"Velocity") == 0)
        type = 7;
    else if (strcmp (&bc_data[0],"Acceleration") == 0)
        type = 8;
    else if (strcmp (&bc_data[0],"Temperature") == 0)
        type = 9;
    else if (strcmp (&bc_data[0],"Pressure") == 0)
        type = 6;
    else if (strcmp (&bc_data[0],"HeatFlux") == 0)
        type = 10;
    else type = 0; //UnknownSet
    
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
    
    //Open mesh_file_name.txt for writing
    ofstream myfile;
    myfile.open ((string(mesh_file_name)+".txt").c_str());

  cout << "<<<< NodeSets >>>>>" << endl;
  //NodeSets
  for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,NodeSet,it)) {
    cout << *it << endl;
    ierr = it->print_Cubit_bc_data(cout); CHKERRQ(ierr);
    vector<char> bc_data;
    ierr = it->get_Cubit_bc_data(bc_data); CHKERRQ(ierr);
    if(bc_data.empty()) continue;
      
      //Displacement
      if (strcmp (&bc_data[0],"Displacement") == 0)
      {
          displacement_cubit_bc_data displacement_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&displacement_cubit_bc_struct); CHKERRQ(ierr);
          //Print data
          cout << displacement_cubit_bc_struct;
          myfile << displacement_cubit_bc_struct;
      }
      
      //Force
      else if (strcmp (&bc_data[0],"Force") == 0)
      {
          force_cubit_bc_data force_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&force_cubit_bc_struct); CHKERRQ(ierr);
          //Print data
          cout << force_cubit_bc_struct;
          myfile << force_cubit_bc_struct;
      }
      
      //Velocity
      else if (strcmp (&bc_data[0],"Velocity") == 0)
      {
          velocity_cubit_bc_data velocity_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&velocity_cubit_bc_struct); CHKERRQ(ierr);
          //Print data
          cout << velocity_cubit_bc_struct;
          myfile << velocity_cubit_bc_struct;
      }
      
      //Acceleration
      else if (strcmp (&bc_data[0],"Acceleration") == 0)
      {
          acceleration_cubit_bc_data acceleration_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&acceleration_cubit_bc_struct); CHKERRQ(ierr);
          //Print data
          cout << acceleration_cubit_bc_struct;
          myfile << acceleration_cubit_bc_struct;
      }
      
      //Temperature
      else if (strcmp (&bc_data[0],"Temperature") == 0)
      {
          temperature_cubit_bc_data temperature_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&temperature_cubit_bc_struct); CHKERRQ(ierr);
          //Print data
          cout << temperature_cubit_bc_struct;
          myfile << temperature_cubit_bc_struct;
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
    if(bc_data.empty()) continue;
      
      //Pressure
      if (strcmp (&bc_data[0],"Pressure") == 0)
      {
          pressure_cubit_bc_data pressure_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&pressure_cubit_bc_struct); CHKERRQ(ierr);
          //Print data
          cout << pressure_cubit_bc_struct;
          myfile << pressure_cubit_bc_struct;
      }

      //Heat Flux
      else if (strcmp (&bc_data[0],"HeatFlux") == 0)
      {
          heatflux_cubit_bc_data heatflux_cubit_bc_struct;
          //Fill bc data
          ierr = func(bc_data,&heatflux_cubit_bc_struct); CHKERRQ(ierr);
          //Print data
          cout << heatflux_cubit_bc_struct;
          myfile << heatflux_cubit_bc_struct;
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
    if(bc_data.empty()) continue;
  } 

    //Close mesh_file_name.txt
    myfile.close();

  PetscFinalize();

}

