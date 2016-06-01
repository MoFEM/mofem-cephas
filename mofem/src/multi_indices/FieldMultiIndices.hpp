/** \file FieldMultiIndices.hpp
 * \brief Myltindex containes, for mofem fields data structures and other low-level functions
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

#ifndef __FIELDMULTIINDICES_HPP__
#define __FIELDMULTIINDICES_HPP__

namespace MoFEM {

/** \brief user adjacency function table
  * \ingroup dof_multi_indices
  */
typedef int (*FieldOrderTable[MBMAXTYPE])(const int order);

/** \brief user adjacency function
  * \ingroup fe_multi_indices
  */
typedef int (*FieldOrderFunct)(const int order);

/**
  * \brief Provide data structure for (tensor) field approximation.
  * \ingroup dof_multi_indices

  The Field is intended to provide support for fields, with a strong bias
  towards supporting first and best the capabilities required for scientific
  computing applications. Since we work with discrete spaces, data structure
  has to carry information about type of approximation space, its regularity

  Note: Some concepts and ideas are taken from iFiedl Interface specification
  <https://redmine.scorec.rpi.edu/anonsvn/itaps/software/trunk/tools/doxygen/html/ifield.html>

  */
struct Field {

  moab::Interface &moab;

  EntityHandle meshSet; 		///< keeps entities for this meshset
  boost::shared_ptr<CoordSys> coordSysPtr;

  Tag th_FieldData,th_AppOrder;
  Tag th_AppDofOrder,th_DofRank;

  BitFieldId* tag_id_data; 		///< tag keeps field id
  FieldSpace* tag_space_data;		///< tag keeps field space
  FieldApproximationBase* tag_base_data;		///< tag keeps field space
  FieldCoefficientsNumber* tag_nb_coeff_data; 	///< tag keeps field rank (dimension, f.e. Temperature field has rank 1, displacements field in 3d has rank 3)
  const void* tag_name_data; 		///< tag keeps name of the field
  int tag_name_size; 			///< number of bits necessary to keep field name
  const void* tag_name_prefix_data; 	///< tag keeps name prefix of the field
  int tag_name_prefix_size; 		///< number of bits necessary to keep field name prefix
  FieldOrderTable forder_table;		///< nb. dofs table for entities
  unsigned int bit_number;


  /**
    * \brief constructor for moab field
    *
    * \param meshset which keeps entities for this field
    */
  Field(Interface &moab,const EntityHandle meshset,const boost::shared_ptr<CoordSys> coord_sys_ptr);

  inline EntityHandle get_meshset() const { return meshSet; };

  /**
    * \brief Get dimension of general two-point tensor \ref MoFEM::CoordSys::getDim

    See details here \ref MoFEM::CoordSys::getDim

    */
  inline int getCoordSysDim(const int d = 0) const { return coordSysPtr->getDim(d); }
  inline PetscErrorCode get_E_Base(const double m[]) const {
    PetscFunctionBegin;
    PetscFunctionReturn(coordSysPtr->get_E_Base(m));
  }
  inline PetscErrorCode get_E_DualBase(const double m[]) const {
    PetscFunctionBegin;
    PetscFunctionReturn(coordSysPtr->get_E_DualBase(m));
  }
  inline PetscErrorCode get_e_Base(const double m[]) const {
    PetscFunctionBegin;
    PetscFunctionReturn(coordSysPtr->get_e_Base(m));
  }
  inline PetscErrorCode get_e_DualBase(const double m[]) const {
    PetscFunctionBegin;
    PetscFunctionReturn(coordSysPtr->get_e_DualBase(m));
  }

  inline EntityHandle getCoordSysMeshSet() const { return coordSysPtr->getMeshSet(); }
  inline std::string getCoordSysName() const { return coordSysPtr->getName(); };
  inline boost::string_ref getCoordSysNameRef() const {
    return coordSysPtr->getNameRef();
  };

  inline const BitFieldId& get_id() const { return *((BitFieldId*)tag_id_data); };
  inline boost::string_ref get_name_ref() const { return boost::string_ref((char *)tag_name_data,tag_name_size); };
  inline std::string get_name() const { return std::string((char *)tag_name_data,tag_name_size); };
  inline FieldSpace get_space() const { return *tag_space_data; };
  inline FieldApproximationBase get_approx_base() const { return *tag_base_data; };

  DEPRECATED inline FieldCoefficientsNumber get_max_rank() const { return *tag_nb_coeff_data; };

  /* \brief get number of field coefficients
  */
  inline FieldCoefficientsNumber get_nb_of_coeffs() const { return *tag_nb_coeff_data; };


  /**
    * \brief Get number of set bit in Field ID.
    * Each field has uid, get get_bit_number get number of bit set for given field. Field ID has only one bit set for each field.
    */
  inline unsigned int get_bit_number() const {
    return bit_number;
  }

  /**
    * \brief Calculate number of set bit in Field ID.
    * Each field has uid, get get_bit_number get number of bit set for given field. Field ID has only one bit set for each field.
    */
  inline unsigned int get_bit_number_calculate() const {
    int b = ffsl(((BitFieldId*)tag_id_data)->to_ulong());
    if(b != 0) return b;
    for(int ll = 1;ll<BITFIELDID_SIZE/32;ll++) {
      BitFieldId id;
      id = (*tag_id_data)>>ll*32;
      b = ll*32+ffsl(id.to_ulong());
      if(b!=0) return b;
    }
    return 0;
  }

  friend std::ostream& operator<<(std::ostream& os,const Field& e);
};

/**
 * \brief interface for Field
 * \ingroup dof_multi_indices
 */
template <typename T>
struct interface_Field {
  const boost::shared_ptr<T> sFieldPtr;

  interface_Field(const boost::shared_ptr<T> field_ptr): sFieldPtr(field_ptr) {};

  inline EntityHandle get_meshset() const { return this->sFieldPtr->get_meshset(); };

  inline int getCoordSysId() const { return this->sFieldPtr->getCoordSysId(); }

  /**
    * \brief Get dimension of general two-point tensor \ref MoFEM::CoordSys::getDim

    See details here \ref MoFEM::CoordSys::getDim

    */
  inline int getCoordSysDim(const int d = 0) const { return this->sFieldPtr->getCoordSysDim(d); }
  inline PetscErrorCode get_E_Base(const double m[]) const {
    PetscFunctionBegin;
    PetscFunctionReturn(this->sFieldPtr->get_E_Base(m));
  }
  inline PetscErrorCode get_E_DualBase(const double m[]) const {
    PetscFunctionBegin;
    PetscFunctionReturn(this->sFieldPtr->get_E_DualBase(m));
  }
  inline PetscErrorCode get_e_Base(const double m[]) const {
    PetscFunctionBegin;
    PetscFunctionReturn(this->sFieldPtr->get_e_Base(m));
  }
  inline PetscErrorCode get_e_DualBase(const double m[]) const {
    PetscFunctionBegin;
    PetscFunctionReturn(this->sFieldPtr->get_e_DualBase(m));
  }
  inline EntityHandle getCoordSysMeshSet() const { return this->sFieldPtr->getCoordSysMeshSet(); }
  inline std::string getCoordSysName() const { return this->sFieldPtr->getCoordSysName(); };
  inline boost::string_ref getCoordSysNameRef() const {
    return this->sFieldPtr->getCoordSysNameRef();
  }

  inline const BitFieldId& get_id() const { return this->sFieldPtr->get_id(); }
  inline unsigned int get_bit_number() const { return this->sFieldPtr->get_bit_number(); }
  inline boost::string_ref get_name_ref() const { return this->sFieldPtr->get_name_ref(); }
  inline std::string get_name() const { return this->sFieldPtr->get_name(); }
  inline FieldSpace get_space() const { return this->sFieldPtr->get_space(); }
  inline FieldApproximationBase get_approx_base() const { return this->sFieldPtr->get_approx_base(); }

  DEPRECATED inline FieldCoefficientsNumber get_max_rank() const { return this->sFieldPtr->get_nb_of_coeffs(); }

  /* \brief get number of field coefficients
  */
  inline FieldCoefficientsNumber get_nb_of_coeffs() const { return this->sFieldPtr->get_nb_of_coeffs(); }

  inline const boost::shared_ptr<T> get_Field_ptr() const { return this->sFieldPtr; }

};

/**
 * @relates multi_index_container
 * \brief Field_multiIndex for Field
 *
 */
typedef multi_index_container<
  boost::shared_ptr<Field>,
  indexed_by<
    hashed_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<Field,const BitFieldId&,&Field::get_id>, HashBit<BitFieldId>, EqBit<BitFieldId> >,
    ordered_unique<
      tag<Meshset_mi_tag>, member<Field,EntityHandle,&Field::meshSet> >,
    ordered_unique<
      tag<FieldName_mi_tag>, const_mem_fun<Field,boost::string_ref,&Field::get_name_ref> >,
    ordered_non_unique<
      tag<BitFieldId_space_mi_tag>, const_mem_fun<Field,FieldSpace,&Field::get_space> >
  > > Field_multiIndex;

typedef multi_index_container<
  boost::shared_ptr<Field>,
  indexed_by<
    ordered_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<Field,const BitFieldId&,&Field::get_id>, LtBit<BitFieldId>
    >
> > Field_multiIndex_view;

/** \brief Set field coordinate system
 * \ingroup ent_multi_indices
  */
struct FieldChangeCoordinateSystem {
  boost::shared_ptr<CoordSys> csPtr;
  FieldChangeCoordinateSystem(boost::shared_ptr<CoordSys> cs_ptr):
  csPtr(cs_ptr) {
  }
  void operator()(boost::shared_ptr<Field> &e) {
    e->coordSysPtr = csPtr;
  }
};

}

#endif // __FIELDMULTIINDICES_HPP__
