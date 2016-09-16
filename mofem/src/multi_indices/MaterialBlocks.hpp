/** \file MaterialBlocks.hpp
 * \brief Data structures for Meshset/Blocsk with material data
 *
 * Notes:
 * - use BLOCK_ATTRIBUTES tag to store data structures
 * - data structures are tags of meshsets

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

#ifndef __MATERIALBLOCKS_HPP__
#define __MATERIALBLOCKS_HPP__

namespace MoFEM {

/*! \struct GenericAttributeData
 *  \brief Generic attribute data structure
 */
struct GenericAttributeData {
    PetscErrorCode ierr;

    /**
     * \brief get data from structure
     * @param  attributes vector of doubles
     * @return            error code
     */
    virtual PetscErrorCode fill_data(const std::vector<double>& attributes) {
      PetscFunctionBegin;
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"It makes no sense for the generic attribute type");
      PetscFunctionReturn(0);
    }

    /**
     * \brief set data on structure
     * @param  tag_ptr pointer to tag on meshset
     * @param  size    size of data in bytes
     * @return         error code
     */
    virtual PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
      PetscFunctionBegin;
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"It makes no sense for the generic attribute type");
      PetscFunctionReturn(0);
    }

    const CubitBCType tYpe; ///< Type of data (f.e. MAT_ELATIC)

    /**
     * \brief get data type
     * @return data type, see CubitBC
     */
    virtual const CubitBCType& getType() const { return tYpe; }

    unsigned int minNumberOfAtributes; ///< minimal number of attributes

    /**
     * \brief get minimal number of attributes which blockset has to have
     * @return number of minimal data attributes
     */
    virtual unsigned int getMinMumberOfAtributes() const {
      return minNumberOfAtributes;
    }

    /**
     * \brief get data structure size
     * @return size of structure in bytes
     */
    virtual std::size_t getSizeOfData() const = 0;

    /**
     * \brief get pointer to data structure
     * @return pointer
     */
    virtual const void * getDataPtr() const = 0;

    GenericAttributeData(const CubitBCType type,const unsigned int min_number_of_atributes):
    tYpe(type),
    minNumberOfAtributes(min_number_of_atributes) {}

};

/** \brief Arbitrary block attributes  data structure
  */
struct BlockSetAttributes: public GenericAttributeData {

    /** \brief generic block attributes
      *
      */
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

    std::size_t getSizeOfData() const { return sizeof(_data_); }
    const void * getDataPtr() const { return &data; }

    BlockSetAttributes():
    GenericAttributeData(BLOCKSET,0) {
      bzero(&data,sizeof(data));
    }

    PetscErrorCode fill_data(const std::vector<double>& attributes) {
      PetscFunctionBegin;
      if(8*attributes.size()>sizeof(data)) {
        SETERRQ(
          PETSC_COMM_SELF,
          MOFEM_DATA_INCONSISTENCY,
          "data inconsistency, please review the number of material properties defined"
        );
      }
      bzero(&data,sizeof(data));
      memcpy(&data, &attributes[0],8*attributes.size());
      PetscFunctionReturn(0);
    }
    PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
      PetscFunctionBegin;
      if(size>sizeof(data)) {
        SETERRQ(
          PETSC_COMM_SELF,
          MOFEM_DATA_INCONSISTENCY,
          "data inconsistency, please review the number of material properties defined"
        );
      }
      memcpy(tag_ptr,&data,size);
      PetscFunctionReturn(0);
    }

    /*! \brief Print data
     */
    friend std::ostream& operator<<(std::ostream& os,const BlockSetAttributes& e);

};

/*! \struct Mat_Elastic
 *  \brief Elastic material data structure
 */
struct Mat_Elastic: public GenericAttributeData {

    /** \brief block tag data structute
      *
      */
    struct __attribute__ ((packed)) _data_{
        double Young; 			///< Young's modulus
        double Poisson; 		///< Poisson's ratio
        double ThermalExpansion;	///< Thermal expansion
        double User1; // User attribute 2 // For some models is reserved for density
        double User2; // User attribute 3
        double User3; // User attribute 4
        double User4; // User attribute 5
        double User5; // User attribute 6
        double User6; // User attribute 7
        double User7; // User attribute 8
    };

    _data_ data;
    std::size_t getSizeOfData() const { return sizeof(data); };
    const void * getDataPtr() const { return   &data; }

    Mat_Elastic():
    GenericAttributeData(MAT_ELASTICSET,2) {
      bzero(&data,sizeof(data));
    };

    PetscErrorCode fill_data(const std::vector<double>& attributes) {
      PetscFunctionBegin;
      if(attributes.size()<minNumberOfAtributes) {
        SETERRQ2(
          PETSC_COMM_SELF,
          MOFEM_DATA_INCONSISTENCY,
          "Young modulus and/or Poisson ratio is not defined. (top tip: check number of ELASTIC block attributes) %d !< %d",
          attributes.size(),minNumberOfAtributes
        );
      }
      if(8*attributes.size()>sizeof(data)) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency, please review the number of material properties defined");
      }
      bzero(&data,sizeof(data));
      memcpy(&data, &attributes[0],8*attributes.size());
      PetscFunctionReturn(0);
    }
    PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
      PetscFunctionBegin;
      if(size>sizeof(data)) {
        SETERRQ(
          PETSC_COMM_SELF,
          MOFEM_DATA_INCONSISTENCY,
          "data inconsistency, please review the number of material properties defined"
        );
      }
      memcpy(tag_ptr,&data,size);
      PetscFunctionReturn(0);
    }

    /*! \brief Print Mat_Elastic data
     */
    friend std::ostream& operator<<(std::ostream& os,const Mat_Elastic& e);

};


/*! \struct Mat_Thermal
 *  \brief Thermal material data structure
 */
struct Mat_Thermal: public GenericAttributeData {

  /** \brief thermal block attributes
    *
    */
  struct __attribute__ ((packed)) _data_{
    double Conductivity; ///< Thermal conductivity
    double HeatCapacity; ///< Heat Capacity
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
  std::size_t getSizeOfData() const { return sizeof(data); }
  const void * getDataPtr() const { return &data; }

  Mat_Thermal():
  GenericAttributeData(MAT_THERMALSET,2) {
    bzero(&data,sizeof(data));
  }

  PetscErrorCode fill_data(const std::vector<double>& attributes) {
    PetscFunctionBegin;
    if(attributes.size()<minNumberOfAtributes) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Thermal conductivity is not defined. (top tip: check number of THERMAL block atributes)");
    }
    if(8*attributes.size()>sizeof(data)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency, please review the number of material properties defined");
    }
    bzero(&data,sizeof(data));
    memcpy(&data, &attributes[0],8*attributes.size());
    PetscFunctionReturn(0);
  }
  PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
    PetscFunctionBegin;
    if(size>sizeof(data)) {
      SETERRQ(
        PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency, please review the number of material properties defined");
    }
    memcpy(tag_ptr,&data,size);
    PetscFunctionReturn(0);
  }

  /*! \brief Print Mat_Elastic data
  */
  friend std::ostream& operator<<(std::ostream& os,const Mat_Thermal& e);
};


/*! \struct Mat_Moisture
 *  \brief moisture transport material data structure
 */

struct Mat_Moisture: public GenericAttributeData {

  /** \brief moisture block attributes
    *
    */
  struct __attribute__ ((packed)) _data_{
    double Diffusivity; 	///< moisture diffusivity
    double Viscosity;  		///< Viscosity of water
    double Permeability; 	///< Permeability of material
    double User3; // User attribute 3
    double User4; // User attribute 4
    double User5; // User attribute 5
    double User6; // User attribute 6
    double User7; // User attribute 7
    double User8; // User attribute 8
    double User9; // User attribute 9
  };

  _data_ data;
  std::size_t getSizeOfData() const { return sizeof(data); }
  const void * getDataPtr() const { return &data; }

  Mat_Moisture():
  GenericAttributeData(MAT_MOISTURESET,1) {
    bzero(&data,sizeof(data));
  }

  PetscErrorCode fill_data(const std::vector<double>& attributes) {
    PetscFunctionBegin;
    if(attributes.size()<minNumberOfAtributes) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"moisture diffusivity is not defined. (top tip: check number of MOISTURE block atributes)");
    }
    if(8*attributes.size()>sizeof(data)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency, please review the number of material properties defined");
    }
    bzero(&data,sizeof(data));
    memcpy(&data, &attributes[0],8*attributes.size());
    PetscFunctionReturn(0);
  }

  /*! \brief Print Mat_Elastic data
   */
  friend std::ostream& operator<<(std::ostream& os,const Mat_Moisture& e);
};



/** \brief Body force data structure
  */
struct Block_BodyForces: public GenericAttributeData {

  /** \brief body forces
    *
    */
  struct __attribute__ ((packed)) _data_{
    double density; 		///< material density
    double acceleration_x; 	///< acceleration X
    double acceleration_y; 	///< acceleration Y
    double acceleration_z; 	///< acceleration Z
    double User4; // User attribute 4
    double User5; // User attribute 5
    double User6; // User attribute 6
    double User7; // User attribute 7
    double User8; // User attribute 8
  };

  _data_ data;
  std::size_t getSizeOfData() const { return sizeof(data); }
  const void * getDataPtr() const { return &data; }

  Block_BodyForces():
  GenericAttributeData(BODYFORCESSET,4) {}

  PetscErrorCode fill_data(const std::vector<double>& attributes) {
    PetscFunctionBegin;
    if(attributes.size()<minNumberOfAtributes) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Material density and/or acceleration is not defined. (top tip: check number of THERMAL block atributes)");
    }
    if(8*attributes.size()>sizeof(data)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency, please review the number of material properties defined");
    }
    bzero(&data,sizeof(data));
    memcpy(&data, &attributes[0],8*attributes.size());
    PetscFunctionReturn(0);
  }
  PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
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
  friend std::ostream& operator<<(std::ostream& os,const Block_BodyForces& e);
};



/*! \struct Mat_Elastic_TransIso
 *  \brief Transverse Isotropic material data structure
 */
  struct Mat_Elastic_TransIso: public GenericAttributeData {

    /** \brief transverse isotropic
      *
      */
    struct __attribute__ ((packed)) _data_{
      double Youngp; 	///< Young's modulus in xy plane (Ep)
      double Youngz; 	///< Young's modulus in z-direction (Ez)
      double Poissonp; 	///< Poisson's ratio in xy plane (vp)
      double Poissonpz; ///< Poisson's ratio in z-direction (vpz)
      double Shearzp; 	///< Shear modulus in z-direction (Gzp)
    };

    _data_ data;
    std::size_t getSizeOfData() const { return sizeof(data); }
    const void * getDataPtr() const { return &data; }

    Mat_Elastic_TransIso():
    GenericAttributeData(MAT_ELASTICSET,5) {
      bzero(&data,sizeof(data));
    }

    PetscErrorCode fill_data(const std::vector<double>& attributes) {
      PetscFunctionBegin;
      //Fill data
      if(attributes.size()<minNumberOfAtributes) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"All material data not defined");
      }
      if(8*attributes.size()!=sizeof(data)) {
        SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, please review the number of material properties defined");
      }
      memcpy(&data, &attributes[0], sizeof(data));
      bzero(&data,sizeof(data));
      memcpy(&data, &attributes[0],8*attributes.size());

      PetscFunctionReturn(0);
    }
    PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
      PetscFunctionBegin;
      if(size>sizeof(data)) {
        SETERRQ(
          PETSC_COMM_SELF,
          MOFEM_DATA_INCONSISTENCY,
          "data inconsistency, please review the number of material properties defined"
        );
      }
      memcpy(tag_ptr,&data,size);
      PetscFunctionReturn(0);
    }

    /*! \brief Print Mat_Elastic_TransIso data
     */
    friend std::ostream& operator<<(std::ostream& os,const Mat_Elastic_TransIso& e);

  };

/*! \struct Mat_Interf
 *  \brief Linear interface data structure
 */
struct Mat_Interf: public GenericAttributeData {

  /** \brief inteface
    *
    */
  struct __attribute__ ((packed)) _data_{
    double alpha; ///< Elastic modulus multiplier
    double beta;  ///< Damage Coupling multiplier between normal and shear (g=sqrt(gn^2 + beta(gt1^2 + gt2^2)))
    double ft;    ///< Maximum stress of crack
    double Gf;    ///< Fracture Energy
  };

  _data_ data;
  std::size_t getSizeOfData() const { return sizeof(data); }
  const void * getDataPtr() const { return &data; }

  Mat_Interf():
  GenericAttributeData(MAT_INTERFSET,4) {
    bzero(&data,sizeof(data));
  }

  virtual PetscErrorCode fill_data(const std::vector<double>& attributes) {
    PetscFunctionBegin;
    //Fill data
    if(8*attributes.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, please review the number of material properties defined");
    memcpy(&data, &attributes[0], sizeof(data));
    PetscFunctionReturn(0);
  }
  virtual PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
    PetscFunctionBegin;
    if(size>sizeof(data)) {
      SETERRQ(
        PETSC_COMM_SELF,
        MOFEM_DATA_INCONSISTENCY,
        "data inconsistency, please review the number of material properties defined"
      );
    }
    memcpy(tag_ptr,&data,size);
    PetscFunctionReturn(0);
  }

  /*! \brief Print Mat_Interf data
    */
  friend std::ostream& operator<<(std::ostream& os,const Mat_Interf& e);
};

/** \brief Mat_Elastic with Fibres
 *  \brief Elastic material data structure
 */
struct Mat_Elastic_EberleinHolzapfel1: public GenericAttributeData {

    /** \brief Hotzapler soft tissue
      *
      */
    struct __attribute__ ((packed)) _data_{
      double Young; 	///< Young's modulus
      double Poisson; 	///< Poisson's ratio
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
    std::size_t getSizeOfData() const { return sizeof(data); }
    const void * getDataPtr() const { return &data; }

    Mat_Elastic_EberleinHolzapfel1():
    GenericAttributeData(MAT_ELASTICSET,10) {
      bzero(&data,sizeof(data));
    }

    PetscErrorCode fill_data(const std::vector<double>& attributes) {
      PetscFunctionBegin;
      if(attributes.size()<minNumberOfAtributes) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"All material data not defined");
      }
      if(8*attributes.size()>sizeof(data)) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency, please review the number of material properties defined");
      }
      bzero(&data,sizeof(data));
      memcpy(&data, &attributes[0],8*attributes.size());
      PetscFunctionReturn(0);
    }


    /*! \brief Print Mat_Elastic data
     */
    friend std::ostream& operator<<(std::ostream& os,const Mat_Elastic_EberleinHolzapfel1& e);

};

}

#endif // __MATERIALBLOCKS_HPP__
