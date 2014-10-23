/** \file BCMultiIndices.hpp
 * \brief Multi-index containers, data boundary data structures and other low-level functions
 *
 */ 

/*
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 *
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
 * It can be freely used for educational and research purposes 
 * by other institutions. If you use this softwre pleas cite my work. 
 *
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __BCMULTIINDICES_HPP__
#define __BCMULTIINDICES_HPP__

namespace MoFEM {

/** 
 * \typedef CubitBC_BitSet
 * bc & material meshsets
 *
 */
typedef bitset<32> CubitBC_BitSet;

/**
  * Tyeps of sets and boundary conditions
  *
  */
enum CubitBC {
  UNKNOWNSET = 0,
  NODESET = 1<<0,
  SIDESET = 1<<1,
  BLOCKSET = 1<<2,
  MATERIALSET = 1<<3,
  DISPLACEMENTSET = 1<<4,
  FORCESET = 1<<5,
  PRESSURESET = 1<<6,
  VELOCITYSET = 1<<7,
  ACCELERATIONSET = 1<<8,
  TEMPERATURESET = 1<<9,
  HEATFLUXSET = 1<<10,
  INTERFACESET = 1<<11,
  UNKNOWNCUBITNAME = 1<< 12,
  MAT_ELASTICSET = 1<<13,	///< block name is "MAT_ELASTIC"
  MAT_INTERFSET = 1 <<14,
  MAT_THERMALSET = 1<<15,	///< block name is "MAT_THERMAL"
  BLOCK_BODYFORCESSET = 1<<16,	///< block name is "BODY_FORCES"
  MAT_MOISTURESET = 1<<17, ///< block name is "MAT_MOISTURE"
  LASTCUBITSET
};

/*! \struct GenericAttributeData
 *  \brief Generic attribute data structure
 */
struct GenericAttributeData {
    PetscErrorCode ierr;
    
    virtual PetscErrorCode fill_data(const vector<double>& attributes) {
      PetscFunctionBegin;
      SETERRQ(PETSC_COMM_SELF,1,"It makes no sense for the generic attribute type");
      PetscFunctionReturn(0);
    }
    virtual PetscErrorCode set_data(void *tag_ptr,unsigned int size) {
      PetscFunctionBegin;
      SETERRQ(PETSC_COMM_SELF,1,"It makes no sense for the generic attribute type");
      PetscFunctionReturn(0);
    }

};

/** \brief Arbitrary block atributes  data structure
  */
struct BlockSetAttributes: public GenericAttributeData {

    struct __attribute__ ((packed)) _data_{
        double User1; // User attribute 1
        double User2; // User attribute 2
        double User3; // User attribute 3
        double User4; // User attribute 4
        double User5; // User attribute 5
        double User6; // User attribute 6
        double User7; // User attribute 7
        double User8; // User attribute 8
        double User9; // User attribute 9
        double User10; // User attribute 10
    };
    
    _data_ data;
    
    const CubitBC_BitSet type;
    const unsigned int min_number_of_atributes;
    BlockSetAttributes(): type(BLOCKSET),min_number_of_atributes(0) {};
    
    virtual PetscErrorCode fill_data(const vector<double>& attributes) {
      PetscFunctionBegin;
      if(8*attributes.size()>sizeof(data)) {
	SETERRQ(PETSC_COMM_SELF,1,
	  "data inconsistency, please review the number of material properties defined");
      }
      bzero(&data,sizeof(data));
      memcpy(&data, &attributes[0],8*attributes.size());
      PetscFunctionReturn(0);
    }
    virtual PetscErrorCode set_data(void *tag_ptr,unsigned int size) {
      PetscFunctionBegin;
      if(size>sizeof(data)) {
	SETERRQ(PETSC_COMM_SELF,1,
	  "data inconsistency, please review the number of material properties defined");
      }
      memcpy(tag_ptr,&data,size);
      PetscFunctionReturn(0);
    }
    
    /*! \brief Print data
     */
    friend ostream& operator<<(ostream& os,const BlockSetAttributes& e);
    
};

/*! \struct Mat_Elastic
 *  \brief Elastic material data structure
 */
struct Mat_Elastic: public GenericAttributeData {
    struct __attribute__ ((packed)) _data_{
        double Young; 			// Young's modulus
        double Poisson; 		// Poisson's ratio
        double ThermalExpansion;	// Thermal expansion
        double User1; // User attribute 2
        double User2; // User attribute 3
        double User3; // User attribute 4
        double User4; // User attribute 5
        double User5; // User attribute 6
        double User6; // User attribute 7
        double User7; // User attribute 8
    };
    
    _data_ data;
    
    const CubitBC_BitSet type;
    const unsigned int min_number_of_atributes;
    Mat_Elastic(): type(MAT_ELASTICSET),min_number_of_atributes(2) {};
    
    virtual PetscErrorCode fill_data(const vector<double>& attributes) {
        PetscFunctionBegin;
        if(attributes.size()<min_number_of_atributes) {
	  SETERRQ(PETSC_COMM_SELF,1,"Young modulus and/or Poisson ratio is not defined. (top tip: check number of ELASTIC block atributes)");
	}
        if(8*attributes.size()>sizeof(data)) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, please review the number of material properties defined");
	}
	bzero(&data,sizeof(data));
        memcpy(&data, &attributes[0],8*attributes.size());
        PetscFunctionReturn(0);
    }
    virtual PetscErrorCode set_data(void *tag_ptr,unsigned int size) {
      PetscFunctionBegin;
      if(size>sizeof(data)) {
	SETERRQ(PETSC_COMM_SELF,1,
	  "data inconsistency, please review the number of material properties defined");
      }
      memcpy(tag_ptr,&data,size);
      PetscFunctionReturn(0);
    }
    
    /*! \brief Print Mat_Elastic data
     */
    friend ostream& operator<<(ostream& os,const Mat_Elastic& e);
    
};
    
    
/*! \struct Mat_Thermal
*  \brief Thermal material data structure
*/
struct Mat_Thermal: public GenericAttributeData {
  struct __attribute__ ((packed)) _data_{
    double Conductivity; // Thermal conductivity
    double HeatCapacity; // Heat Capacity
    double User2; // User attribute 2
    double User3; // User attribute 3
    double User4; // User attribute 4
    double User5; // User attribute 5
    double User6; // User attribute 6
    double User7; // User attribute 7
    double User8; // User attribute 8
    double User9; // User attribute 9
  };

  _data_ data;
        
  const CubitBC_BitSet type;
  const unsigned int min_number_of_atributes;
  Mat_Thermal(): type(MAT_THERMALSET),min_number_of_atributes(2) {};
        
  virtual PetscErrorCode fill_data(const vector<double>& attributes) {
    PetscFunctionBegin;
    if(attributes.size()<min_number_of_atributes) {
      SETERRQ(PETSC_COMM_SELF,1,"Thermal conductivity is not defined. (top tip: check number of THERMAL block atributes)");
    }
    if(8*attributes.size()>sizeof(data)) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, please review the number of material properties defined");
    }
    bzero(&data,sizeof(data));
    memcpy(&data, &attributes[0],8*attributes.size());
    PetscFunctionReturn(0);
  }
  virtual PetscErrorCode set_data(void *tag_ptr,unsigned int size) {
    PetscFunctionBegin;
    if(size>sizeof(data)) {
      SETERRQ(PETSC_COMM_SELF,1,
	"data inconsistency, please review the number of material properties defined");
    }
    memcpy(tag_ptr,&data,size);
    PetscFunctionReturn(0);
  }

  /*! \brief Print Mat_Elastic data
  */
  friend ostream& operator<<(ostream& os,const Mat_Thermal& e);
};
  
 
  /*! \struct Mat_Moisture
   *  \brief moisture transport material data structure
   */
  
  struct Mat_Moisture: public GenericAttributeData {
    struct __attribute__ ((packed)) _data_{
      double Diffusivity; // moisture diffusivity
      double Viscosity; // Viscosity of water
      double Permeability; // Permeability of material
      double User3; // User attribute 3
      double User4; // User attribute 4
      double User5; // User attribute 5
      double User6; // User attribute 6
      double User7; // User attribute 7
      double User8; // User attribute 8
      double User9; // User attribute 9
    };
    
    _data_ data;
    
    const CubitBC_BitSet type;
    const unsigned int min_number_of_atributes;
    Mat_Moisture(): type(MAT_MOISTURESET),min_number_of_atributes(1) {};
    
    virtual PetscErrorCode fill_data(const vector<double>& attributes) {
      PetscFunctionBegin;
      if(attributes.size()<min_number_of_atributes) {
        SETERRQ(PETSC_COMM_SELF,1,"moisture diffusivity is not defined. (top tip: check number of MOISTURE block atributes)");
      }
      if(8*attributes.size()>sizeof(data)) {
        SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, please review the number of material properties defined");
      }
      bzero(&data,sizeof(data));
      memcpy(&data, &attributes[0],8*attributes.size());
      PetscFunctionReturn(0);
    }
    
    /*! \brief Print Mat_Elastic data
     */
    friend ostream& operator<<(ostream& os,const Mat_Moisture& e);
  };

  
  
/** \brief Body force data structure
  */
struct Block_BodyForces: public GenericAttributeData {
  struct __attribute__ ((packed)) _data_{
    double density; // Thermal conductivity
    double acceleration_x; // User attribute 1
    double acceleration_y; // User attribute 2
    double acceleration_z; // User attribute 3
    double User4; // User attribute 4
    double User5; // User attribute 5
    double User6; // User attribute 6
    double User7; // User attribute 7
    double User8; // User attribute 8
  };

  _data_ data;
        
  const CubitBC_BitSet type;
  const unsigned int min_number_of_atributes;
  Block_BodyForces(): type(BLOCK_BODYFORCESSET),min_number_of_atributes(4) {};
        
  virtual PetscErrorCode fill_data(const vector<double>& attributes) {
    PetscFunctionBegin;
    if(attributes.size()<min_number_of_atributes) {
      SETERRQ(PETSC_COMM_SELF,1,"Material density and/or acceleration is not defined. (top tip: check number of THERMAL block atributes)");
    }
    if(8*attributes.size()>sizeof(data)) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, please review the number of material properties defined");
    }
    bzero(&data,sizeof(data));
    memcpy(&data, &attributes[0],8*attributes.size());
    PetscFunctionReturn(0);
  }
  virtual PetscErrorCode set_data(void *tag_ptr,unsigned int size) {
    PetscFunctionBegin;
    if(size>sizeof(data)) {
      SETERRQ(PETSC_COMM_SELF,1,
	"data inconsistency, please review the number of material properties defined");
    }
    memcpy(tag_ptr,&data,size);
    PetscFunctionReturn(0);
  }
        
  /*! \brief Print Mat_Elastic data
  */
  friend ostream& operator<<(ostream& os,const Block_BodyForces& e);
};
    
    
    
/*! \struct Mat_Elastic_TransIso
 *  \brief Transverse Isotropic material data structure
 */
  struct Mat_Elastic_TransIso: public Mat_Elastic {
    struct __attribute__ ((packed)) _data_{
      double Youngp; // Young's modulus in xy plane (Ep)
      double Youngz; // Young's modulus in z-direction (Ez)
      double Poissonp; // Poisson's ratio in xy plane (vp)
      double Poissonpz; // Poisson's ratio in z-direction (vpz)
      double Shearzp; // Shear modulus in z-direction (Gzp)
    };
    
    _data_ data;
    
    const unsigned int min_number_of_atributes;
    Mat_Elastic_TransIso(): Mat_Elastic(),min_number_of_atributes(5) {};
    
    virtual PetscErrorCode fill_data(const vector<double>& attributes) {
      PetscFunctionBegin;
      //Fill data
      if(attributes.size()<min_number_of_atributes) {
        SETERRQ(PETSC_COMM_SELF,1,"All material data not defined");
      }
      if(8*attributes.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, please review the number of material properties defined");
      memcpy(&data, &attributes[0], sizeof(data));
      bzero(&data,sizeof(data));
      memcpy(&data, &attributes[0],8*attributes.size());
      
      PetscFunctionReturn(0);
    }
    virtual PetscErrorCode set_data(void *tag_ptr,unsigned int size) {
      PetscFunctionBegin;
      if(size>sizeof(data)) {
	SETERRQ(PETSC_COMM_SELF,1,
	  "data inconsistency, please review the number of material properties defined");
      }
      memcpy(tag_ptr,&data,size);
      PetscFunctionReturn(0);
    }
   
    /*! \brief Print Mat_Elastic_TransIso data
     */
    friend ostream& operator<<(ostream& os,const Mat_Elastic_TransIso& e);
    
  };

/*! \struct Mat_Interf
 *  \brief Linear interface data structure
 */
struct Mat_Interf: public GenericAttributeData {
  struct __attribute__ ((packed)) _data_{
    double alpha; // Elastic modulus multiplier
    double beta;  // Damage Coupling multiplier between normal and shear (g=sqrt(gn^2 + beta(gt1^2 + gt2^2)))
    double ft;    // Maximum stress of crack
    double Gf;    // Fracture Energy
  };
      
  _data_ data;
      
  const CubitBC_BitSet type;
  Mat_Interf(): type(MAT_INTERFSET) {};
      
  virtual PetscErrorCode fill_data(const vector<double>& attributes) {
    PetscFunctionBegin;
    //Fill data
    if(8*attributes.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, please review the number of material properties defined");
    memcpy(&data, &attributes[0], sizeof(data));
    PetscFunctionReturn(0);
  }
  virtual PetscErrorCode set_data(void *tag_ptr,unsigned int size) {
    PetscFunctionBegin;
    if(size>sizeof(data)) {
      SETERRQ(PETSC_COMM_SELF,1,
	"data inconsistency, please review the number of material properties defined");
    }
    memcpy(tag_ptr,&data,size);
    PetscFunctionReturn(0);
  }
      
  /*! \brief Print Mat_Interf data
    */
  friend ostream& operator<<(ostream& os,const Mat_Interf& e);
};

/*! \struct Mat_Elastic with Fibres
 *  \brief Elastic material data structure
 */
struct Mat_Elastic_EberleinHolzapfel1: public Mat_Elastic {
    struct __attribute__ ((packed)) _data_{
        double Young; // Young's modulus
        double Poisson; // Poisson's ratio
        double k1; // User attribute 1
        double k2; // User attribute 2
        double a0x; // User attribute 3
        double a0y; // User attribute 4
        double a0z; // User attribute 5
        double a1x; // User attribute 6
        double a1y; // User attribute 7
        double a1z; // User attribute 8
    };
    
    _data_ data;
    
    const unsigned int min_number_of_atributes;
    Mat_Elastic_EberleinHolzapfel1(): Mat_Elastic(),min_number_of_atributes(10) {};
    
    virtual PetscErrorCode fill_data(const vector<double>& attributes) {
        PetscFunctionBegin;
        if(attributes.size()<min_number_of_atributes) {
	  SETERRQ(PETSC_COMM_SELF,1,"All material data not defined");
	}
        if(8*attributes.size()>sizeof(data)) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, please review the number of material properties defined");
	}
	bzero(&data,sizeof(data));
        memcpy(&data, &attributes[0],8*attributes.size());
        PetscFunctionReturn(0);
    }

    
    /*! \brief Print Mat_Elastic data
     */
    friend ostream& operator<<(ostream& os,const Mat_Elastic_EberleinHolzapfel1& e);
    
};


    
/*! \struct GenericCubitBcData
 *  \brief Generic bc data structure
 */
struct GenericCubitBcData {
    PetscErrorCode ierr;
    
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        SETERRQ(PETSC_COMM_SELF,1,"It makes no sense for the generic bc type");
        PetscFunctionReturn(0);
    }
    
};

/*! \struct DisplacementCubitBcData
 *  \brief Definition of the displacement bc data structure
 */
struct DisplacementCubitBcData: public GenericCubitBcData {
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

    const CubitBC_BitSet type;
    DisplacementCubitBcData(): type(DISPLACEMENTSET) {};
    
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
	if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
 
    /*! \brief Print displacement bc data
     */
    friend ostream& operator<<(ostream& os,const DisplacementCubitBcData& e);
    
};

/*! \struct ForceCubitBcData
 *  \brief Definition of the force bc data structure
 */
struct ForceCubitBcData: public GenericCubitBcData {
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
    const CubitBC_BitSet type;
    ForceCubitBcData(): type(FORCESET) {};

    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
	if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
 
    /*! \brief Print force bc data
    */
    friend ostream& operator<<(ostream& os,const ForceCubitBcData& e);
    
};

/*! \struct VelocityCubitBcData
 *  \brief Definition of the velocity bc data structure
 */
struct VelocityCubitBcData: public GenericCubitBcData {
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
    const CubitBC_BitSet type;
    VelocityCubitBcData(): type(VELOCITYSET) {};
   
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
	if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
 
    /*! \brief Print velocity bc data
    */
    friend ostream& operator<<(ostream& os,const VelocityCubitBcData& e);
    
};  

/*! \struct AccelerationCubitBcData
 *  \brief Definition of the acceleration bc data structure
 */    
struct AccelerationCubitBcData: public GenericCubitBcData {
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
    const CubitBC_BitSet type;
    AccelerationCubitBcData(): type(ACCELERATIONSET) {};

    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
	if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
    
    /*! \brief Print acceleration bc data
    */
    friend ostream& operator<<(ostream& os,const AccelerationCubitBcData& e);
    
};

/*! \struct TemperatureCubitBcData
 *  \brief Definition of the temperature bc data structure
 */
struct TemperatureCubitBcData: public GenericCubitBcData {
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
    const CubitBC_BitSet type;
    TemperatureCubitBcData(): type(TEMPERATURESET) {};

    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
      PetscFunctionBegin;
      //Fill data
      if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      memcpy(&data, &bc_data[0], sizeof(data));
      PetscFunctionReturn(0);
    }
    
    /*! \brief Print temperature bc data
    */
    friend ostream& operator<<(ostream& os,const TemperatureCubitBcData& e);
};

/*! \struct PressureCubitBcData
 *  \brief Definition of the pressure bc data structure
 */
struct PressureCubitBcData: public GenericCubitBcData {
    struct __attribute__ ((packed)) _data_{
    char name[8]; // 8 characters for "Pressure"
    char flag1; // This is always zero
    char flag2; // 0: Pressure is interpeted as pure pressure 1: pressure is interpreted as total force
    double value1; // Pressure value
    char zero; // This is always zero
    };
    
    _data_ data;
    const CubitBC_BitSet type;
    PressureCubitBcData(): type(PRESSURESET) {};
   
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
	if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
 
    /*! \brief Print pressure bc data
    */
    friend ostream& operator<<(ostream& os,const PressureCubitBcData& e);
    
};

/*! \struct HeatfluxCubitBcData
 *  \brief Definition of the heat flux bc data structure
 */
struct HeatfluxCubitBcData: public GenericCubitBcData {
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
    const CubitBC_BitSet type;
    HeatfluxCubitBcData(): type(HEATFLUXSET) {};

    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
	if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
 
    /*! \brief Print heat flux bc data
    */
    friend ostream& operator<<(ostream& os,const HeatfluxCubitBcData& e);
    
};

/*! \struct CfgCubitBcData
 *  \brief Definition of the cfd_bc data structure
 */
struct CfgCubitBcData: public GenericCubitBcData {
    struct __attribute__ ((packed)) _data_{
        char name[6]; // 6 characters for "cfd_bc"
        char zero; // This is always zero
        char type; // This is the type of cfd_bc
    };
    
    _data_ data;
    const CubitBC_BitSet type;
    CfgCubitBcData(): type(INTERFACESET) {};
    
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
        if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
    
    /*! \brief Print cfd_bc data
     */
    friend ostream& operator<<(ostream& os,const CfgCubitBcData& e);

};
    
/** 
 * \brief this struct keeps basic methods for moab meshset about material and boundary conditions
 */
struct CubitMeshSets {
  EntityHandle meshset;
  CubitBC_BitSet CubitBCType;
  vector<Tag> tag_handles;
  int *msId;
  char* tag_bc_data;
  int tag_bc_size;
  unsigned int *tag_block_header_data;
  double* tag_block_attributes;
  int tag_block_attributes_size;
  char* tag_name_data;
  const CubitBC_BitSet meshsets_mask;
  CubitMeshSets(Interface &moab,const EntityHandle _meshset);
  CubitMeshSets(Interface &moab,const CubitBC_BitSet _CubitBCType,const int _msId);
  inline int get_msId() const { return *msId; }
  inline CubitBC_BitSet get_CubitBCType() const { return CubitBCType; }

  inline EntityHandle get_meshset() const { return meshset; }
  inline unsigned long int get_CubitBCType_ulong() const { return CubitBCType.to_ulong(); }
  inline unsigned long int get_CubitBCType_mask_meshset_types_ulong() const { return (CubitBCType&meshsets_mask).to_ulong(); }
  inline unsigned long int get_CubitBCType_bc_data_types_ulong() const { return (CubitBCType&(~meshsets_mask)).to_ulong(); }

  PetscErrorCode get_Cubit_msId_entities_by_dimension(Interface &moab,const int dimension,Range &entities,const bool recursive = false) const;
  PetscErrorCode get_Cubit_msId_entities_by_dimension(Interface &moab,Range &entities,const bool recursive = false)  const;

  /** 
   *  \brief Function that returns the CubitBC_BitSet type of the contents of bc_data
   */
  PetscErrorCode get_type_from_bc_data(const vector<char> &bc_data,CubitBC_BitSet &type) const;

  /** 
   *  \brief Function that returns the CubitBC_BitSet type of the contents of bc_data
  */
  PetscErrorCode get_type_from_bc_data(CubitBC_BitSet &type) const;
    
  /**
   * \brief get bc_data vector from MoFEM database
   * 
   * \param bc_data is the in/out vector were bc_data will be stored
   */
  PetscErrorCode get_Cubit_bc_data(vector<char>& bc_data) const;
    
  /**
  * \brief get block_headers vector from MoFEM database
  *
  * \param material_data is the in/out vector were the material data will be stored
  */
  PetscErrorCode get_Cubit_block_header_data(vector<unsigned int>& material_data) const;

  /**
  * \brief print material_data int stream given by os
  *
  * f.e. it->print_Cubit_material_data(cout), i.e. printing to standard output
  * f.e. it->print_Cubit_material_data(cerr), i.e. printing to standard error output
  */
  PetscErrorCode print_Cubit_block_header_data(ostream& os) const;
    
  /**
   * \brief print bc_data int stream given by os
   *
   * f.e. it->print_Cubit_bc_data(cout), i.e. printing to standard output
   * f.e. it->print_Cubit_bc_data(cerr), i.e. printing to standard error output
   */
  PetscErrorCode print_Cubit_bc_data(ostream& os) const;

  template<class _CUBIT_BC_DATA_TYPE_>
  PetscErrorCode get_cubit_bc_data_structure(_CUBIT_BC_DATA_TYPE_& data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if((CubitBCType&data.type).none()) {
      SETERRQ(PETSC_COMM_SELF,1,"bc_data are not for _CUBIT_BC_DATA_TYPE_ structure");  
    }
    vector<char> bc_data;
    get_Cubit_bc_data(bc_data);
    ierr = data.fill_data(bc_data); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
   *  \brief Function that returns the CubitBC_BitSet type of the block name, sideset name etc.
   */
  PetscErrorCode get_type_from_Cubit_name(const string &name,CubitBC_BitSet &type) const;

  /**
   *  \brief Function that returns the CubitBC_BitSet type of the block name, sideset name etc.
   */
  PetscErrorCode get_type_from_Cubit_name(CubitBC_BitSet &type) const;
    
  /**
   * \brief get Cubit block attributes
   *
   * \param attributes is the vector where the block attribute data will be stored
   */
  PetscErrorCode get_Cubit_attributes(vector<double> &attributes) const;

  /**
   * \brief print the attributes vector
   *
   * f.e. it->print_Cubit_attributes(cout), i.e. printing to standard output
   * f.e. it->print_Cubit_attributes(cerr), i.e. printing to standard error output
   */
  PetscErrorCode print_Cubit_attributes(ostream& os) const;

  /**
   * \brief get name of block, sideset etc. (this is set in Cubit block properties)
   *
   * Block Name Conventions:
   * -----------------------
   * Materials are defined with block names starting with MAT_
   * e.g. MAT_ELASTIC_abcd, MAT_FRACTcdef etc.
   * Solution procedures are defined with block names starting with SOL_ e.g.
   * SOL_ELASTIC_xx, SOL_NLELASTICxx, SOL_FRACTabcd etc.
   *
   * List of materials/solution procedures
   * ---------------------------------------------------------------------------
   * Block name /  Number of attributes  / (1) Attribute 1, (2) Attribute 2 etc.
   * ---------------------------------------------------------------------------
   *
   * MAT_ELASTIC / 10 /  (1) Young's  modulus
   *                    (2) Poisson's ratio
   *                    (3) User attribute 8
   *                    ...
   *                    (10) User attribute 8
   *
   * MAT_ELASTIC_TRANSISO / 5 / (1) Young's modulus in xy plane (Ep)
   *                    (2) Young's modulus in z-direction (Ez)
   *                    (3) Poisson's ratio in xy plane (vp)
   *                    (4) Poisson's ratio in z-direction (vpz)
   *                    (5) Shear modulus in z-direction (Gzp)
   *
   * MAT_INTERF / 1 /   (1) Elastic modulus multiplier
   *
   * To be extended as appropriate
   */
   string get_Cubit_name() const;

  /**
   * \brief print name of block, sideset etc. (this is set in Cubit setting properties)
   *
   * e.g. it->print_Cubit_name(cout), i.e. printing to standard output
   * e.g it->print_Cubit_name(cerr), i.e. printing to standard error output
   */
  PetscErrorCode print_Cubit_name(ostream& os) const;
    
  template<class _ATTRIBUTE_TYPE_>
  PetscErrorCode get_attribute_data_structure(_ATTRIBUTE_TYPE_ &data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if((CubitBCType&data.type).none()) {
        SETERRQ(PETSC_COMM_SELF,1,"attributes are not for _ATTRIBUTE_TYPE_ structure");
    }
    vector<double> attributes;
    ierr = get_Cubit_attributes(attributes); CHKERRQ(ierr);
    ierr = data.fill_data(attributes); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  template<class _ATTRIBUTE_TYPE_>
  PetscErrorCode set_attribute_data_structure(_ATTRIBUTE_TYPE_ &data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if((CubitBCType&data.type).none()) {
        SETERRQ(PETSC_COMM_SELF,1,"attributes are not for _ATTRIBUTE_TYPE_ structure");
    }
    double *ptr = const_cast<double*>(tag_block_attributes);
    ierr = data.set_data(ptr,8*tag_block_attributes_size); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
     
  friend ostream& operator<<(ostream& os,const CubitMeshSets& e);

  private:
  Tag nsTag,ssTag,nsTag_data,ssTag_data,bhTag,bhTag_header,block_attribs,entityNameTag;
  PetscErrorCode get_tags_hanlders(Interface &moab);
    
};
    
/**
 * @relates multi_index_container
 * \brief CubitMeshSet_multiIndex
 *
 * \param hashed_unique<
      tag<Meshset_mi_tag>, member<CubitMeshSets,EntityHandle,&CubitMeshSets::meshset> >,
 * \param ordered_non_unique<
      tag<CubitMeshSets_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_ulong> >,
 * \param ordered_non_unique<
      tag<CubitMeshSets_mask_meshset_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_mask_meshset_types_ulong> >,
 * \param ordered_non_unique<
      tag<CubitMeshSets_bc_data_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_bc_data_types_ulong> >,
 * \param ordered_non_unique<
      tag<CubitMeshSets_name>, const_mem_fun<CubitMeshSets,string,&CubitMeshSets::get_Cubit_name> >,
 *
 * \param    hashed_unique<
      tag<Composite_mi_tag>,       
      composite_key<
	CubitMeshSets, <br>
	  const_mem_fun<CubitMeshSets,int,&CubitMeshSets::get_msId>,
	  const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_ulong> > >
 *
 */
typedef multi_index_container<
  CubitMeshSets,
  indexed_by<
    hashed_unique<
      tag<Meshset_mi_tag>, member<CubitMeshSets,EntityHandle,&CubitMeshSets::meshset> >,
    ordered_non_unique<
      tag<CubitMeshSets_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_ulong> >,
    ordered_non_unique<
      tag<CubitMeshSets_mask_meshset_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_mask_meshset_types_ulong> >,
    ordered_non_unique<
      tag<CubitMeshSets_bc_data_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_bc_data_types_ulong> >,
    ordered_non_unique<
      tag<CubitMeshSets_name>, const_mem_fun<CubitMeshSets,string,&CubitMeshSets::get_Cubit_name> >,
    hashed_unique<
      tag<Composite_Cubit_msId_and_MeshSetType_mi_tag>,       
      composite_key<
	CubitMeshSets,
	  const_mem_fun<CubitMeshSets,int,&CubitMeshSets::get_msId>,
	  const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_mask_meshset_types_ulong> > >
  > > CubitMeshSet_multiIndex;

}

#endif // __BCMULTIINDICES_HPP__
