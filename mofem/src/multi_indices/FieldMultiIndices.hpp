/** \file FieldMultiIndices.hpp
 * \brief Field data structure storing information about space, approximation
 * base, coordinate systems, etc.
 *
 * Also, it stores data needed for book keeping, like tags to data on the
 * mesh.
 *
 * Each filed has unique ID and name. This data structure is shared between
 * entities on which is spanning and DOFs on those entities.
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

#ifndef __FIELDMULTIINDICES_HPP__
#define __FIELDMULTIINDICES_HPP__

namespace MoFEM {

/** \brief user adjacency function
 * \ingroup fe_multi_indices
 */
typedef boost::function<int(const int order)> FieldOrderFunct;

/** \brief user adjacency function table
 * \ingroup dof_multi_indices
 */
typedef FieldOrderFunct FieldOrderTable[MBMAXTYPE];

struct FieldEntity;
struct DofEntity;

/**
 * \brief Provide data structure for (tensor) field approximation.
 * \ingroup dof_multi_indices
 *
 * The Field is intended to provide support for fields, with a strong bias
 * towards supporting first and best the capabilities required for scientific
 * computing applications. Since we work with discrete spaces, data structure
 * has to carry information about type of approximation space, its regularity.
 *
 *
 * Field data structure storing information about space, approximation base,
 * coordinate systems, etc. It stores additional data needed for book keeping,
 * like tags to data on the mesh.
 *
 * Each filed has unique ID and name. This
 * data structure is shared between entities on which is spanning and DOFs on
 * those entities.
 *
 */
struct Field {

  virtual ~Field() = default;

  typedef multi_index_container<boost::weak_ptr<std::vector<FieldEntity>>,
                                indexed_by<sequenced<>>>
      SequenceEntContainer;

  typedef multi_index_container<boost::weak_ptr<std::vector<DofEntity>>,
                                indexed_by<sequenced<>>>
      SequenceDofContainer;

  typedef std::array<std::array<int, MAX_DOFS_ON_ENTITY>, MBMAXTYPE>
      DofsOrderMap;

  moab::Interface &moab;

  EntityHandle meshSet; ///< keeps entities for this meshset
  boost::shared_ptr<CoordSys>
      coordSysPtr; ///< Pointer to field coordinate system data structure

  TagType th_FieldDataVertsType; // Tag type for storing data on vertices
  Tag th_FieldDataVerts; ///< Tag storing field values on vertices in the field
  Tag th_FieldData;      ///< Tag storing field values on entity in the field
  Tag th_AppOrder;       ///< Tag storing approximation order on entity

  BitFieldId *tagId;                   ///< tag keeps field id
  FieldSpace *tagSpaceData;            ///< tag keeps field space
  FieldApproximationBase *tagBaseData; ///< tag keeps field base

  /// tag keeps field rank (dimension, f.e. Temperature field has rank 1,
  /// displacements field in 3d has rank 3)
  FieldCoefficientsNumber *tagNbCoeffData;
  const void *tagName; ///< tag keeps name of the field
  int tagNameSize;     ///< number of bits necessary to keep field name
  const void *tagNamePrefixData; ///< tag keeps name prefix of the field
  int tagNamePrefixSize;       ///< number of bits necessary to keep field name
                               ///< prefix
  FieldOrderTable forderTable; ///< nb. DOFs table for entities

  /**
   * @brief Get the Field Order Table
   *
   * @return FieldOrderTable&
   */
  inline FieldOrderTable &getFieldOrderTable() { return forderTable; }

  /**
   * Field Id is bit set. Each field has only one bit on, bitNumber stores
   * number of set bit
   */
  unsigned int bitNumber;

  static UId generateGlobalUniqueIdForTypeLo(const char bit_number,
                                             const EntityType type,
                                             const int owner_proc) {
    constexpr int ent_shift = 8 * sizeof(EntityHandle);
    return (static_cast<UId>(type) << MB_ID_WIDTH |
            static_cast<UId>(bit_number) << 8 * sizeof(EntityHandle) |
            static_cast<UId>(owner_proc) << 5 + ent_shift)
           << 9;
  }

  UId generateGlobalUniqueIdForTypeLo(const EntityType type,
                                      const int owner_proc) const {
    return generateGlobalUniqueIdForTypeLo(bitNumber, type, owner_proc);
  }

  static UId generateGlobalUniqueIdForTypeHi(const char bit_number,
                                             const EntityType type,
                                             const int owner_proc) {
    constexpr int ent_shift = 8 * sizeof(EntityHandle);
    return (static_cast<UId>(type) << MB_ID_WIDTH |
            static_cast<UId>(bit_number) << ent_shift |
            static_cast<UId>(owner_proc) << 5 + ent_shift)
           << 9;
  }

  UId generateGlobalUniqueIdForTypeHi(const EntityType type,
                                      const int owner_proc) const {
    return generateGlobalUniqueIdForTypeHi(bitNumber, type, owner_proc);
  }

  /**
   * \brief constructor for moab field
   *
   * \param meshset which keeps entities for this field
   */
  Field(const moab::Interface &moab, const EntityHandle meshset,
        const boost::shared_ptr<CoordSys> coord_sys_ptr);

  /**
   * \brief Get field meshset
   *

   * To meshsets entity are attached Tags which keeps basic information about
   * field. Those information is field name, approximation base, approximation
   * space, id, etc.

   * In meshset contains entities on which given filed is sparing. Type of
   entities
   * depended on approximations space.

   * @return EntityHandle
   */
  inline EntityHandle getMeshset() const { return meshSet; }

  /**
    * \brief Get dimension of general two-point tensor \ref
    MoFEM::CoordSys::getDim

    See details here \ref MoFEM::CoordSys::getDim

    */
  inline int getCoordSysDim(const int d = 0) const {
    return coordSysPtr->getDim(d);
  }

  /**
   * \brief   Get reference base vectors
   * @param   Array where coefficients (covariant) are returned
   * @return  Error code
   */
  inline MoFEMErrorCode get_E_Base(const double m[]) const {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(coordSysPtr->get_E_Base(m));
  }

  /**
   * \brief   Get reference dual base vectors
   * @param   Array where coefficients (contravariant) are returned
   * @return  Error code
   */
  inline MoFEMErrorCode get_E_DualBase(const double m[]) const {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(coordSysPtr->get_E_DualBase(m));
  }

  /**
   * \brief   Get current dual base vectors
   * @param   Array where coefficients (covariant) are returned
   * @return  Error code
   */
  inline MoFEMErrorCode get_e_Base(const double m[]) const {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(coordSysPtr->get_e_Base(m));
  }

  /**
   * \brief   Get current dual base vectors
   * @param   Array where coefficients (covariant) are returned
   * @return  Error code
   */
  inline MoFEMErrorCode get_e_DualBase(const double m[]) const {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(coordSysPtr->get_e_DualBase(m));
  }

  /**
   * \brief Returns meshset on which Tags defining coordinate system are stored
   * @return Coordinate system EntityHandle
   */
  inline EntityHandle getCoordSysMeshSet() const {
    return coordSysPtr->getMeshset();
  }

  /**
   * \brief   Get coordinate system name
   * @return  Coordinate system name
   */
  inline std::string getCoordSysName() const { return coordSysPtr->getName(); }

  /**
   * \brief   Get coordinate system name
   * @return Return string_ref with name.
   */
  inline boost::string_ref getCoordSysNameRef() const {
    return coordSysPtr->getNameRef();
  }

  /**
   * \brief Get unique field id.
   * @return Filed ID
   */
  inline const BitFieldId &getId() const { return *((BitFieldId *)tagId); }

  /**
   * \brief Get string reference to field name
   * @return Field name
   */
  inline boost::string_ref getNameRef() const {
    return boost::string_ref((char *)tagName, tagNameSize);
  }

  /**
   * \brief   Get field name
   * @return  Field name
   */
  inline std::string getName() const {
    return std::string((char *)tagName, tagNameSize);
  }

  /**
   * \brief   Get field approximation space
   * @return  approximation space
   */
  inline FieldSpace getSpace() const { return *tagSpaceData; }

  /**
   * \brief   Get approximation base
   * @return  Approximation base
   */
  inline FieldApproximationBase getApproxBase() const { return *tagBaseData; }

  /** \brief Get number of field coefficients
    *

    * Scalar field has only one coefficient, vector field in 3D has three. In
    * general number determine space needed to keep data on entities. What
    coefficient
    * means depend on interpretation and associated coordinate system. For
    example
    * 3 coefficient means could be covariant or contravariant, or mean three
    temperatures
    * for mixture of solid, air and water, etc.


  */
  inline FieldCoefficientsNumber getNbOfCoeffs() const {
    return *tagNbCoeffData;
  };

  /**
   * \brief Get number of set bit in Field ID.
   * Each field has uid, get getBitNumber get number of bit set for given field.
   * Field ID has only one bit set for each field.
   */
  inline unsigned int getBitNumber() const { return bitNumber; }

  /**
   * \brief Calculate number of set bit in Field ID.
   * Each field has uid, get getBitNumber get number of bit set for given field.
   * Field ID has only one bit set for each field.
   */
  inline unsigned int getBitNumberCalculate() const {
    int b = ffsl(((BitFieldId *)tagId)->to_ulong());
    if (b != 0)
      return b;
    for (int ll = 1; ll < BITFIELDID_SIZE / 32; ll++) {
      BitFieldId id;
      id = (*tagId) >> ll * 32;
      b = ll * 32 + ffsl(id.to_ulong());
      if (b != 0)
        return b;
    }
    return 0;
  }

  friend std::ostream &operator<<(std::ostream &os, const Field &e);

  /**
   * \brief Get reference to sequence data container
   *
   * In sequence data container data are physically stored. The purpose of this
   * is to allocate MoFEMEntities data in bulk, having only one allocation
   * instead each time entity is inserted. That makes code efficient.
   *
   * The vector in sequence is destroyed if last entity inside that vector is
   * destroyed. All MoFEM::MoFEMEntities have aliased shared_ptr which points to
   * the vector.
   *
   * @return MoFEM::Field::SequenceEntContainer
   */
  inline SequenceEntContainer &getEntSequenceContainer() const {
    return sequenceEntContainer;
  }

  /**
   * \brief Get reference to sequence data container
   *
   * In sequence data container data are physically stored. The purpose of this
   * is to allocate DofEntity data in bulk, having only one allocation instead
   * each time entity is inserted. That makes code efficient.
   *
   * The vector in sequence is destroyed if last entity inside that vector is
   * destroyed. All MoFEM::MoFEMEntities have aliased shared_ptr which points to
   the vector.

   * Not all DOFs are starred in this way, currently such cases are considered;
   * - DOFs on vertices. That is exploited that for H1 space, there is some
   * fixed number of DOFs on each vertex

   * For other cases, DOFs are stored locally in each MoFEM::MoFEMEntities.

   * @return MoFEM::Field::SequenceDofContainer
   */
  inline SequenceDofContainer &getDofSequenceContainer() const {
    return sequenceDofContainer;
  }

  /**
   * \brief get hash-map relating dof index on entity with its order
   *
   * Dofs of given field are indexed on entity
   * of the same type, same space, approximation base and number of
   * coefficients, are sorted in the way.
   *
   */
  inline std::array<int, MAX_DOFS_ON_ENTITY> &
  getDofOrderMap(const EntityType type) const {
    return dofOrderMap[type];
  }

  /**
   * \brief get hash-map relating dof index on entity with its order
   *
   * Dofs of given field are indexed on entity
   * of the same type, same space, approximation base and number of
   * coefficients, are sorted in the way.
   *
   */
  inline DofsOrderMap &getDofOrderMap() const {
    return const_cast<DofsOrderMap &>(dofOrderMap);
  }

  MoFEMErrorCode rebuildDofsOrderMap() const;

private:
  mutable SequenceEntContainer sequenceEntContainer;
  mutable SequenceDofContainer sequenceDofContainer;
  mutable DofsOrderMap dofOrderMap;

};

/**
 * \brief Pointer interface for MoFEM::Field
 *
 * MoFEM::Field class is keeps data and methods. This class is interface to that
 * class, and all other classes, like MoFEMEntities, DofEntity and derived form
 * them inherits pointer interface, not MoFEM::Field class directly.
 *
 * \ingroup dof_multi_indices
 */
template <typename T> struct interface_Field {

  mutable boost::shared_ptr<T> sFieldPtr;

  interface_Field(const boost::shared_ptr<T> &field_ptr)
      : sFieldPtr(field_ptr) {}

  interface_Field(const interface_Field<T> &interface)
      : sFieldPtr(interface.getFieldPtr()) {}

  virtual ~interface_Field() {}

  inline EntityHandle getMeshset() const {
    return this->sFieldPtr->getMeshset();
  }

  inline int getCoordSysId() const { return this->sFieldPtr->getCoordSysId(); }

  /**
    * \brief Get dimension of general two-point tensor \ref
    MoFEM::CoordSys::getDim

    See details here \ref MoFEM::CoordSys::getDim

    */
  inline int getCoordSysDim(const int d = 0) const {
    return this->sFieldPtr->getCoordSysDim(d);
  }

  inline MoFEMErrorCode get_E_Base(const double m[]) const {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(this->sFieldPtr->get_E_Base(m));
  }
  inline MoFEMErrorCode get_E_DualBase(const double m[]) const {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(this->sFieldPtr->get_E_DualBase(m));
  }
  inline MoFEMErrorCode get_e_Base(const double m[]) const {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(this->sFieldPtr->get_e_Base(m));
  }

  inline MoFEMErrorCode get_e_DualBase(const double m[]) const {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(this->sFieldPtr->get_e_DualBase(m));
  }

  /// @return return meshset for coordinate system
  inline EntityHandle getCoordSysMeshSet() const {
    return this->sFieldPtr->getCoordSysMeshSet();
  }

  /// @return return coordinate system name for field
  inline std::string getCoordSysName() const {
    return this->sFieldPtr->getCoordSysName();
  }

  /// @return return coordinate system name for field
  inline boost::string_ref getCoordSysNameRef() const {
    return this->sFieldPtr->getCoordSysNameRef();
  }

  /// @return get field Id
  inline const BitFieldId &getId() const { return this->sFieldPtr->getId(); }

  /// @return get field name
  inline boost::string_ref getNameRef() const {
    return this->sFieldPtr->getNameRef();
  }

  /// @return get field name
  inline std::string getName() const { return this->sFieldPtr->getName(); }

  /// @return get approximation space
  inline FieldSpace getSpace() const { return this->sFieldPtr->getSpace(); }

  /// @return get approximation base
  inline FieldApproximationBase getApproxBase() const {
    return this->sFieldPtr->getApproxBase();
  }

  /// @return get number of coefficients for DOF
  inline FieldCoefficientsNumber getNbOfCoeffs() const {
    return this->sFieldPtr->getNbOfCoeffs();
  }

  /// @return get bit number if filed Id
  inline unsigned int getBitNumber() const {
    return this->sFieldPtr->getBitNumber();
  }

  /// @return get pointer to the field data structure
  inline boost::shared_ptr<T> &getFieldPtr() const { return this->sFieldPtr; }

  /**
   * \brief get hash-map relating dof index on entity with its order
   *
   * Dofs of given field are indexed on entity
   * of the same type, same space, approximation base and number of
   * coefficients, are sorted in the way.
   *
   */
  inline std::vector<ApproximationOrder> &
  getDofOrderMap(const EntityType type) const {
    return this->sFieldPtr->getDofOrderMap(type);
  }
};

/**
 * @relates multi_index_container
 * \brief Field_multiIndex for Field
 *
 */
typedef multi_index_container<
    boost::shared_ptr<Field>,
    indexed_by<
        hashed_unique<tag<BitFieldId_mi_tag>,
                      const_mem_fun<Field, const BitFieldId &, &Field::getId>,
                      HashBit<BitFieldId>, EqBit<BitFieldId>>,
        ordered_unique<tag<Meshset_mi_tag>,
                       member<Field, EntityHandle, &Field::meshSet>>,
        ordered_unique<
            tag<FieldName_mi_tag>,
            const_mem_fun<Field, boost::string_ref, &Field::getNameRef>>,
        ordered_non_unique<tag<BitFieldId_space_mi_tag>,
                           const_mem_fun<Field, FieldSpace, &Field::getSpace>>>>
    Field_multiIndex;

typedef multi_index_container<
    boost::shared_ptr<Field>,
    indexed_by<
        ordered_unique<tag<BitFieldId_mi_tag>,
                       const_mem_fun<Field, const BitFieldId &, &Field::getId>,
                       LtBit<BitFieldId>>>>
    Field_multiIndex_view;

/** \brief Set field coordinate system
 * \ingroup ent_multi_indices
 */
struct FieldChangeCoordinateSystem {
  boost::shared_ptr<CoordSys> csPtr;
  FieldChangeCoordinateSystem(const boost::shared_ptr<CoordSys> &cs_ptr)
      : csPtr(cs_ptr) {}
  void operator()(boost::shared_ptr<Field> &e) { e->coordSysPtr = csPtr; }
};


} // namespace MoFEM

#endif // __FIELDMULTIINDICES_HPP__
  