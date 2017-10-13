/** \file BCData.hpp
 * \brief Data strucures with Cubit native blocks/meshets with boundary conditions
 *
 */

/*
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

#ifndef __BCDATA_HPP__
#define __BCDATA_HPP__

namespace MoFEM {

/*! \struct GenericCubitBcData
 * \brief Generic bc data structure
 * \ingroup mofem_bc
 */
struct GenericCubitBcData {
    

    /**
     * \brief get data from structure
     * @param  attributes vector of doubles
     * @return            error code
     */
    virtual PetscErrorCode fill_data(const std::vector<char>& bc_data) {
      MoFEMFunctionBeginHot;
      SETERRQ(PETSC_COMM_SELF,1,"It makes no sense for the generic bc type");
      MoFEMFunctionReturnHot(0);
    }

    /**
     * \brief set data on structure
     * @param  tag_ptr pointer to tag on meshset
     * @param  size    size of data in bytes
     * @return         error code
     */
    virtual PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
      MoFEMFunctionBeginHot;
      SETERRQ(PETSC_COMM_SELF,1,"It makes no sense for the generic bc type");
      MoFEMFunctionReturnHot(0);
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

    const CubitBCType tYpe; ///< Type of boundary condition

    /**
     * \brief get data type
     * @return data type, see CubitBC
     */
    virtual const CubitBCType& getType() const { return tYpe; }

    GenericCubitBcData(const CubitBCType type):
    tYpe(type) {}

    virtual ~GenericCubitBcData() {}

};

/*! \struct DisplacementCubitBcData
 *  \brief Definition of the displacement bc data structure
 * \ingroup mofem_bc
 */
struct DisplacementCubitBcData: public GenericCubitBcData {

    /*! \brief attributes of DisplacementCubitBcData
    */
    struct __attribute__ ((packed)) _data_{
    char name[12]; //< 12 characters for "Displacement"
    char pre1; //< Always zero
    char pre2; //< pre-processing flags for modification of displacement bcs. They should not affect analysis, i.e. safe to ignore; 1: smallest combine, 2: average, 3: largest combine, 4: overwrite or no combination defined (default)
    char flag1; //< Flag for X-Translation (0: N/A, 1: specified)
    char flag2; //< Flag for Y-Translation (0: N/A, 1: specified)
    char flag3; //< Flag for Z-Translation (0: N/A, 1: specified)
    char flag4; //< Flag for X-Rotation (0: N/A, 1: specified)
    char flag5; //< Flag for Y-Rotation (0: N/A, 1: specified)
    char flag6; //< Flag for Z-Rotation (0: N/A, 1: specified)
    double value1; //< Value for X-Translation
    double value2; //< Value for Y-Translation
    double value3; //< Value for Z-Translation
    double value4; //< Value for X-Rotation
    double value5; //< Value for Y-Rotation
    double value6; //< Value for Z-Rotation
    };

    _data_ data;

    std::size_t getSizeOfData() const { return sizeof(_data_); }
    const void * getDataPtr() const { return &data; }

    DisplacementCubitBcData():
    GenericCubitBcData(DISPLACEMENTSET) {
      bzero(&data,sizeof(data));
    }

    PetscErrorCode fill_data(const std::vector<char>& bc_data) {
      MoFEMFunctionBeginHot;
      //Fill data
      if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      memcpy(&data, &bc_data[0], sizeof(data));
      MoFEMFunctionReturnHot(0);
    }

    PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
      MoFEMFunctionBeginHot;
      if(size!=sizeof(data)) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      memcpy(tag_ptr, &data, size);
      MoFEMFunctionReturnHot(0);
    }


    /*! \brief Print displacement bc data
     */
    friend std::ostream& operator<<(std::ostream& os,const DisplacementCubitBcData& e);

};

/*! \struct ForceCubitBcData
 *  \brief Definition of the force bc data structure
 * \ingroup mofem_bc
 */
struct ForceCubitBcData: public GenericCubitBcData {
    struct __attribute__ ((packed)) _data_{
    char name[5]; //< 5 characters for "Force"
    char zero[3]; //< 3 zeros
    double value1; //< Force magnitude
    double value2; //< Moment magnitude
    double value3; //< X-component of force direction vector
    double value4; //< Y-component of force direction vector
    double value5; //< Z-component of force direction vector
    double value6; //< X-component of moment direction vector
    double value7; //< Y-component of moment direction vector
    double value8; //< Z-component of moment direction vector
    char zero2; // 0
    };

    _data_ data;
    std::size_t getSizeOfData() const { return sizeof(_data_); }
    const void * getDataPtr() const { return &data; }

    ForceCubitBcData():
    GenericCubitBcData(FORCESET) {
      bzero(&data,sizeof(data));
    }

    PetscErrorCode fill_data(const std::vector<char>& bc_data) {
      MoFEMFunctionBeginHot;
      //Fill data
      if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      memcpy(&data, &bc_data[0], sizeof(data));
      MoFEMFunctionReturnHot(0);
    }

    PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
      MoFEMFunctionBeginHot;
      if(size!=sizeof(data)) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      memcpy(tag_ptr, &data, size);
      MoFEMFunctionReturnHot(0);
    }


    /*! \brief Print force bc data
    */
    friend std::ostream& operator<<(std::ostream& os,const ForceCubitBcData& e);

};

/*! \struct VelocityCubitBcData
 * \brief Definition of the velocity bc data structure
 * \ingroup mofem_bc
 */
struct VelocityCubitBcData: public GenericCubitBcData {
    struct __attribute__ ((packed)) _data_{
    char name[8]; //< 8 characters for "Velocity"
    char pre1; //< Always zero
    char pre2; //< pre-processing flags for modification of displacement bcs. They should not affect analysis, i.e. safe to ignore; 1: smallest combine, 2: average, 3: largest combine, 4: overwrite or no combination defined (default)
    char flag1; //< Flag for X-Translation (0: N/A, 1: specified)
    char flag2; //< Flag for Y-Translation (0: N/A, 1: specified)
    char flag3; //< Flag for Z-Translation (0: N/A, 1: specified)
    char flag4; //< Flag for X-Rotation (0: N/A, 1: specified)
    char flag5; //< Flag for Y-Rotation (0: N/A, 1: specified)
    char flag6; //< Flag for Z-Rotation (0: N/A, 1: specified)
    double value1; //< Value for X-Translation
    double value2; //< Value for Y-Translation
    double value3; //< Value for Z-Translation
    double value4; //< Value for X-Rotation
    double value5; //< Value for Y-Rotation
    double value6; //< Value for Z-Rotation
    };

    _data_ data;
    std::size_t getSizeOfData() const { return sizeof(_data_); }
    const void * getDataPtr() const { return &data; }

    VelocityCubitBcData():
    GenericCubitBcData(VELOCITYSET) {
      bzero(&data,sizeof(data));
    }

    PetscErrorCode fill_data(const std::vector<char>& bc_data) {
      MoFEMFunctionBeginHot;
      //Fill data
      if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      memcpy(&data, &bc_data[0], sizeof(data));
      MoFEMFunctionReturnHot(0);
    }

    PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
      MoFEMFunctionBeginHot;
      if(size!=sizeof(data)) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      memcpy(tag_ptr, &data, size);
      MoFEMFunctionReturnHot(0);
    }

    /*! \brief Print velocity bc data
    */
    friend std::ostream& operator<<(std::ostream& os,const VelocityCubitBcData& e);

};

/*! \struct AccelerationCubitBcData
 * \brief Definition of the acceleration bc data structure
 * \ingroup mofem_bc
 */
struct AccelerationCubitBcData: public GenericCubitBcData {
    struct __attribute__ ((packed)) _data_{
    char name[12]; //< 12 characters for "Acceleration"
    char pre1; //< Always zero
    char pre2; //< pre-processing flags for modification of displacement bcs. They should not affect analysis, i.e. safe to ignore; 1: smallest combine, 2: average, 3: largest combine, 4: overwrite or no combination defined (default)
    char flag1; //< Flag for X-Translation (0: N/A, 1: specified)
    char flag2; //< Flag for Y-Translation (0: N/A, 1: specified)
    char flag3; //< Flag for Z-Translation (0: N/A, 1: specified)
    char flag4; //< Flag for X-Rotation (0: N/A, 1: specified)
    char flag5; //< Flag for Y-Rotation (0: N/A, 1: specified)
    char flag6; //< Flag for Z-Rotation (0: N/A, 1: specified)
    double value1; //< Value for X-Translation
    double value2; //< Value for Y-Translation
    double value3; //< Value for Z-Translation
    double value4; //< Value for X-Rotation
    double value5; //< Value for Y-Rotation
    double value6; //< Value for Z-Rotation
    };

    _data_ data;

    std::size_t getSizeOfData() const { return sizeof(_data_); }
    const void * getDataPtr() const { return &data; }

    AccelerationCubitBcData():
    GenericCubitBcData(ACCELERATIONSET) {
      bzero(&data,sizeof(data));
    }

    PetscErrorCode fill_data(const std::vector<char>& bc_data) {
      MoFEMFunctionBeginHot;
      //Fill data
      if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      memcpy(&data, &bc_data[0], sizeof(data));
      MoFEMFunctionReturnHot(0);
    }

    PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
      MoFEMFunctionBeginHot;
      if(size!=sizeof(data)) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      memcpy(tag_ptr, &data, size);
      MoFEMFunctionReturnHot(0);
    }


    /*! \brief Print acceleration bc data
    */
    friend std::ostream& operator<<(std::ostream& os,const AccelerationCubitBcData& e);

};

/*! \struct TemperatureCubitBcData
 *  \brief Definition of the temperature bc data structure
 * \ingroup mofem_bc
 */
 struct TemperatureCubitBcData: public GenericCubitBcData {
   struct __attribute__ ((packed)) _data_{
     char name[11]; //< 11 characters for "Temperature"
     char pre1; //< This is always zero
     char pre2; //< 0: temperature is not applied on thin shells (default); 1: temperature is applied on thin shells
     char flag1; //< 0: N/A, 1: temperature value applied (not on thin shells)
     char flag2; //< 0: N/A, 1: temperature applied on thin shell middle
     char flag3; //< 0: N/A, 1: thin shell temperature gradient specified
     char flag4; //< 0: N/A, 1: top thin shell temperature
     char flag5; //< 0: N/A, 1: bottom thin shell temperature
     char flag6; //< This is always zero
     double value1; //< Temperature (default case - no thin shells)
     double value2; //< Temperature for middle of thin shells
     double value3; //< Temperature gradient for thin shells
     double value4; //< Temperature for top of thin shells
     double value5; //< Temperature for bottom of thin shells
     double value6; //< This is always zero, i.e. ignore
   };

   _data_ data;

   std::size_t getSizeOfData() const { return sizeof(_data_); }
   const void * getDataPtr() const { return &data; }

   TemperatureCubitBcData():
   GenericCubitBcData(TEMPERATURESET) {
     bzero(&data,sizeof(data));
   }

   PetscErrorCode fill_data(const std::vector<char>& bc_data) {
     MoFEMFunctionBeginHot;
     //Fill data
     if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
     memcpy(&data, &bc_data[0], sizeof(data));
     MoFEMFunctionReturnHot(0);
   }

   PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
     MoFEMFunctionBeginHot;
     if(size!=sizeof(data)) {
       SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
     }
     memcpy(tag_ptr, &data, size);
     MoFEMFunctionReturnHot(0);
   }


   /*! \brief Print temperature bc data
   */
   friend std::ostream& operator<<(std::ostream& os,const TemperatureCubitBcData& e);
 };

/*! \struct PressureCubitBcData
 * \brief Definition of the pressure bc data structure
 * \ingroup mofem_bc
 */
struct PressureCubitBcData: public GenericCubitBcData {
    struct __attribute__ ((packed)) _data_{
    char name[8];   //< 8 characters for "Pressure"
    char flag1;     //< This is always zero
    char flag2;     //< 0: Pressure is interpreted as pure pressure 1: pressure is interpreted as total force
    double value1;  //< Pressure value
    char zero;      //< This is always zero
    };

    _data_ data;

    std::size_t getSizeOfData() const { return sizeof(_data_); }
    const void * getDataPtr() const { return &data; }

    PressureCubitBcData():
    GenericCubitBcData(PRESSURESET) {
      bzero(&data,sizeof(data));
    }

    PetscErrorCode fill_data(const std::vector<char>& bc_data) {
      MoFEMFunctionBeginHot;
      //Fill data
      if(bc_data.size()!=sizeof(data)) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      memcpy(&data, &bc_data[0], sizeof(data));
      MoFEMFunctionReturnHot(0);
    }

    PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
      MoFEMFunctionBeginHot;
      if(size!=sizeof(data)) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      memcpy(tag_ptr, &data, size);
      MoFEMFunctionReturnHot(0);
    }

    /*! \brief Print pressure bc data
    */
    friend std::ostream& operator<<(std::ostream& os,const PressureCubitBcData& e);

};

/*! \struct HeatFluxCubitBcData
 * \brief Definition of the heat flux bc data structure
 * \ingroup mofem_bc
 */
struct HeatFluxCubitBcData: public GenericCubitBcData {

  struct __attribute__ ((packed)) _data_ {
    char name[8]; //< 8 characters for "HeatFlux" (no space)
    char pre1; //< This is always zero
    char pre2; //< 0: heat flux is not applied on thin shells (default); 1: heat flux is applied on thin shells
    char flag1; //< 0: N/A, 1: normal heat flux case (i.e. single value, case without thin shells)
    char flag2; //< 0: N/A, 1: Thin shell top heat flux specified
    char flag3; //< 0: N/A, 1: Thin shell bottom heat flux specidied
    double value1; //< Heat flux value for default case (no thin shells)
    double value2; //< Heat flux (thin shell top)
    double value3; //< Heat flux (thin shell bottom)
  };

  _data_ data;

  std::size_t getSizeOfData() const { return sizeof(_data_); }
  const void * getDataPtr() const { return &data; }

  HeatFluxCubitBcData():
  GenericCubitBcData(HEATFLUXSET) {
    bzero(&data,sizeof(data));
  }

  PetscErrorCode fill_data(const std::vector<char>& bc_data) {
    MoFEMFunctionBeginHot;
    //Fill data
    if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    memcpy(&data, &bc_data[0], sizeof(data));
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
    MoFEMFunctionBeginHot;
    if(size!=sizeof(data)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    memcpy(tag_ptr, &data, size);
    MoFEMFunctionReturnHot(0);
  }

  /*! \brief Print heat flux bc data
  */
  friend std::ostream& operator<<(std::ostream& os,const HeatFluxCubitBcData& e);

};

/*! \struct CfgCubitBcData
 * \brief Definition of the cfd_bc data structure
 * \ingroup mofem_bc
 */
struct CfgCubitBcData: public GenericCubitBcData {
    struct __attribute__ ((packed)) _data_{
        char name[6]; //< 6 characters for "cfd_bc"
        char zero; //< This is always zero
        char type; //< This is the type of cfd_bc
    };

    _data_ data;

    std::size_t getSizeOfData() const { return sizeof(_data_); }
    const void * getDataPtr() const { return &data; }

    CfgCubitBcData():
    GenericCubitBcData(INTERFACESET) {
      bzero(&data,sizeof(data));
    }

    PetscErrorCode fill_data(const std::vector<char>& bc_data) {
      MoFEMFunctionBeginHot;
      //Fill data
      if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      memcpy(&data, &bc_data[0], sizeof(data));
      MoFEMFunctionReturnHot(0);
    }

    PetscErrorCode set_data(void *tag_ptr,unsigned int size) const {
      MoFEMFunctionBeginHot;
      if(size!=sizeof(data)) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      memcpy(tag_ptr, &data, size);
      MoFEMFunctionReturnHot(0);
    }

    /*! \brief Print cfd_bc data
     */
    friend std::ostream& operator<<(std::ostream& os,const CfgCubitBcData& e);

};

}

#endif // __BCMULTIINDICES_HPP__
