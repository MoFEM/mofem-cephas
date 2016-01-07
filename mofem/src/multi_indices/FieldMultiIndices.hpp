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

  The MoFEMField is intended to provide support for fields, with a strong bias
  towards supporting first and best the capabilities required for scientific
  computing applications. Since we work with discrete spaces, data structure
  has to carry information about type of approximation space, its regularity

  Note: Some concepts and ideas are taken from iFiedl Interface specification
  <https://redmine.scorec.rpi.edu/anonsvn/itaps/software/trunk/tools/doxygen/html/ifield.html>

  */
struct MoFEMField {

  EntityHandle meshSet; 		///< keeps entities for this meshset
  const CoordSys *coordSysPtr;

  Tag th_FieldData,th_AppOrder;
  Tag th_AppDofOrder,th_DofRank;

  BitFieldId* tag_id_data; 		///< tag keeps field id
  FieldSpace* tag_space_data;		///< tag keeps field space
  ApproximationRank* tag_rank_data; 	///< tag keeps field rank (dimension, f.e. Temperature field has rank 1, displacements field in 3d has rank 3)
  const void* tag_name_data; 		///< tag keeps name of the field
  int tag_name_size; 			///< number of bits necessary to keep field name
  const void* tag_name_prefix_data; 	///< tag keeps name prefix of the field
  int tag_name_prefix_size; 		///< number of bits necessary to keep field name prefix
  FieldOrderTable forder_table;		///< nb. dofs table for entities

  /**
    * \brief constructor for moab field
    *
    * \param meshset which keeps entities for this field
    */
  MoFEMField(Interface &moab,const EntityHandle meshset,const CoordSys *coord_sys_ptr);

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
  inline string getCoordSysName() const { return coordSysPtr->getName(); };
  inline boost::string_ref getCoordSysNameRef() const {
    return coordSysPtr->getNameRef();
  };

  inline const BitFieldId& get_id() const { return *((BitFieldId*)tag_id_data); };
  inline boost::string_ref get_name_ref() const { return boost::string_ref((char *)tag_name_data,tag_name_size); };
  inline string get_name() const { return string((char *)tag_name_data,tag_name_size); };
  inline FieldSpace get_space() const { return *tag_space_data; };
  DEPRECATED inline ApproximationRank get_max_rank() const { return *tag_rank_data; };

  /* \brief get number of field coefficients
  */
  inline ApproximationRank get_nb_of_coeffs() const { return *tag_rank_data; };


  /**
    * \brief get number of set bit in Field ID.
    * Each field has uid, get get_bit_number get number of bit set for given field. Field ID has only one bit set for each field.
    */
  inline unsigned int get_bit_number() const {
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
  const MoFEMField* get_MoFEMField_ptr() const { return this; };
  friend ostream& operator<<(ostream& os,const MoFEMField& e);
};

/**
 * \brief interface for MoFEMField
 * \ingroup dof_multi_indices
 */
template <typename T>
struct interface_MoFEMField {
  const T *field_ptr;
  interface_MoFEMField(const T *_field_ptr): field_ptr(_field_ptr) {};
  inline EntityHandle get_meshset() const { return field_ptr->get_meshset(); };

  inline int getCoordSysId() const { return field_ptr->getCoordSysId(); }

  /**
    * \brief Get dimension of general two-point tensor \ref MoFEM::CoordSys::getDim

    See details here \ref MoFEM::CoordSys::getDim

    */
  inline int getCoordSysDim(const int d = 0) const { return field_ptr->getCoordSysDim(d); }
  inline PetscErrorCode get_E_Base(const double m[]) const {
    PetscFunctionBegin;
    PetscFunctionReturn(field_ptr->get_E_Base(m));
  }
  inline PetscErrorCode get_E_DualBase(const double m[]) const {
    PetscFunctionBegin;
    PetscFunctionReturn(field_ptr->get_E_DualBase(m));
  }
  inline PetscErrorCode get_e_Base(const double m[]) const {
    PetscFunctionBegin;
    PetscFunctionReturn(field_ptr->get_e_Base(m));
  }
  inline PetscErrorCode get_e_DualBase(const double m[]) const {
    PetscFunctionBegin;
    PetscFunctionReturn(field_ptr->get_e_DualBase(m));
  }
  inline EntityHandle getCoordSysMeshSet() const { return field_ptr->getCoordSysMeshSet(); }
  inline string getCoordSysName() const { return field_ptr->getCoordSysName(); };
  inline boost::string_ref getCoordSysNameRef() const {
    return field_ptr->getCoordSysNameRef();
  };

  inline const BitFieldId& get_id() const { return field_ptr->get_id(); };
  inline unsigned int get_bit_number() const { return field_ptr->get_bit_number(); }
  inline boost::string_ref get_name_ref() const { return field_ptr->get_name_ref(); };
  inline string get_name() const { return field_ptr->get_name(); };
  inline FieldSpace get_space() const { return field_ptr->get_space(); };

  DEPRECATED inline ApproximationRank get_max_rank() const { return field_ptr->get_nb_of_coeffs(); };

  /* \brief get number of field coefficients
  */
  inline ApproximationRank get_nb_of_coeffs() const { return field_ptr->get_nb_of_coeffs(); };


  inline const MoFEMField* get_MoFEMField_ptr() const { return field_ptr->get_MoFEMField_ptr(); };
};

/**
 * @relates multi_index_container
 * \brief MoFEMField_multiIndex for MoFEMField
 *
 */
typedef multi_index_container<
  MoFEMField,
  indexed_by<
    hashed_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<MoFEMField,const BitFieldId&,&MoFEMField::get_id>, HashBit<BitFieldId>, EqBit<BitFieldId> >,
    ordered_unique<
      tag<Meshset_mi_tag>, member<MoFEMField,EntityHandle,&MoFEMField::meshSet> >,
    ordered_unique<
      tag<FieldName_mi_tag>, const_mem_fun<MoFEMField,boost::string_ref,&MoFEMField::get_name_ref> >,
    ordered_non_unique<
      tag<BitFieldId_space_mi_tag>, const_mem_fun<MoFEMField,FieldSpace,&MoFEMField::get_space> >
  > > MoFEMField_multiIndex;

typedef multi_index_container<
  const MoFEMField*,
  indexed_by<
    ordered_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<MoFEMField,const BitFieldId&,&MoFEMField::get_id>, LtBit<BitFieldId>
    >
> > MoFEMField_multiIndex_view;

/** \brief Set field coordinate system
 * \ingroup ent_multi_indices
  */
struct MoFEMFieldChangeCoordinateSystem {
  const CoordSys *csPtr;
  MoFEMFieldChangeCoordinateSystem(const CoordSys *cs_ptr):
  csPtr(cs_ptr) {
  }
  void operator()(MoFEMField &e) {
    e.coordSysPtr = csPtr;
  }
};

}

#endif // __FIELDMULTIINDICES_HPP__
