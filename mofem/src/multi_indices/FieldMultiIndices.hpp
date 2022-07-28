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



#ifndef __FIELDMULTIINDICES_HPP__
#define __FIELDMULTIINDICES_HPP__

namespace MoFEM {

template <typename T> struct interface_RefEntity;
struct DofEntity;

/** \brief user adjacency function
 * \ingroup fe_multi_indices
 */
typedef boost::function<int(const int order)> FieldOrderFunct;

/** \brief user adjacency function table
 * \ingroup dof_multi_indices
 */
typedef FieldOrderFunct FieldOrderTable[MBMAXTYPE];

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

  /**
   * \brief constructor for moab field
   *
   * \param meshset which keeps entities for this field
   */
  Field(moab::Interface &moab, const EntityHandle meshset);

  virtual ~Field() = default;

  using SequenceDofContainer = multi_index_container<

      boost::weak_ptr<std::vector<DofEntity>>,

      indexed_by<sequenced<>>>;

  typedef std::array<std::array<int, MAX_DOFS_ON_ENTITY>, MBMAXTYPE>
      DofsOrderMap;

  moab::Interface &moab;

  EntityHandle meshSet; ///< keeps entities for this meshset

  TagType tagFieldDataVertsType; // Tag type for storing data on vertices
  Tag th_FieldDataVerts; ///< Tag storing field values on vertices in the field
  Tag th_FieldData;      ///< Tag storing field values on entity in the field
  Tag th_AppOrder;       ///< Tag storing approximation order on entity
  Tag th_FieldRank;      /// Tag field rank

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
   * \brief   Get field approximation space
   * @return  approximation space name
   */
  inline auto getSpaceName() const {
    return std::string(FieldSpaceNames[getSpace()]);
  }

  /**
   * \brief   Get approximation base
   * @return  Approximation base
   */
  inline FieldApproximationBase getApproxBase() const { return *tagBaseData; }

  /**
   * \brief   Get approximation base
   * @return  Approximation base name
   */
  inline auto getApproxBaseName() const {
    return std::string(ApproximationBaseNames[getApproxBase()]);
  }

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
  inline FieldBitNumber getBitNumber() const { return bitNumber; }

  /**
   * \brief Calculate number of set bit in Field ID.
   * Each field has uid, get getBitNumber get number of bit set for given field.
   * Field ID has only one bit set for each field.
   */
  static inline FieldBitNumber getBitNumberCalculate(const BitFieldId &id) {
    static_assert(BITFIELDID_SIZE >= 32,
                  "Too many fields allowed, can be more but ...");
    FieldBitNumber b = ffsl(id.to_ulong());
    if (b != 0)
      return b;
    return 0;
  }

  /**
   * \brief Calculate number of set bit in Field ID.
   * Each field has uid, get getBitNumber get number of bit set for given field.
   * Field ID has only one bit set for each field.
   */
  inline FieldBitNumber getBitNumberCalculate() const {
    return getBitNumberCalculate(static_cast<BitFieldId &>(*tagId));
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
  inline const std::array<ApproximationOrder, MAX_DOFS_ON_ENTITY> &
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
  inline const DofsOrderMap &getDofOrderMap() const { return dofOrderMap; }

  MoFEMErrorCode rebuildDofsOrderMap();

  friend std::ostream &operator<<(std::ostream &os, const Field &e);

  inline const Field *getFieldRawPtr() const { return this; };

private:
  mutable SequenceDofContainer sequenceDofContainer;
  mutable DofsOrderMap dofOrderMap;
};

/**
 * \brief Pointer interface for MoFEM::Field
 *
 * MoFEM::Field class is keeps data and methods. This class is interface to
 * that class, and all other classes, like MoFEMEntities, DofEntity and
 * derived form them inherits pointer interface, not MoFEM::Field class
 * directly.
 *
 * \ingroup dof_multi_indices
 */
template <typename FIELD, typename REFENT>
struct interface_FieldImpl : public interface_RefEntity<REFENT> {

  using interface_type_RefEntity = interface_RefEntity<REFENT>;

  interface_FieldImpl(const boost::shared_ptr<FIELD> &field_ptr,
                      const boost::shared_ptr<REFENT> &ref_ents_ptr)
      : interface_RefEntity<REFENT>(ref_ents_ptr) {}
  virtual ~interface_FieldImpl() = default;
};

template <typename FIELD, typename REFENT>
struct interface_Field : public interface_FieldImpl<FIELD, REFENT> {

  interface_Field(const boost::shared_ptr<FIELD> &field_ptr,
                  const boost::shared_ptr<REFENT> &ref_ents_ptr)
      : interface_FieldImpl<FIELD, REFENT>(field_ptr, ref_ents_ptr),
        sFieldPtr(field_ptr) {}

  inline EntityHandle getMeshset() const {
    return getFieldRawPtr()->getMeshset();
  }

  inline int getCoordSysDim(const int d = 0) const {
    return getFieldRawPtr()->getCoordSysDim(d);
  }

  /// @return get field Id
  inline const BitFieldId &getId() const {
    return getFieldRawPtr()->getId();
  }

  /// @return get field name
  inline boost::string_ref getNameRef() const {
    return getFieldRawPtr()->getNameRef();
  }

  /// @return get field name
  inline std::string getName() const {
    return getFieldRawPtr()->getName();
  }

  /// @return get approximation space
  inline FieldSpace getSpace() const {
    return getFieldRawPtr()->getSpace();
  }

  /// @return get approximation base
  inline FieldApproximationBase getApproxBase() const {
    return getFieldRawPtr()->getApproxBase();
  }

  /// @return get number of coefficients for DOF
  inline FieldCoefficientsNumber getNbOfCoeffs() const {
    return getFieldRawPtr()->getNbOfCoeffs();
  }

  /// @return get bit number if filed Id
  inline FieldBitNumber getBitNumber() const {
    return getFieldRawPtr()->getBitNumber();
  }

  /**
   * \brief get hash-map relating dof index on entity with its order
   *
   * Dofs of given field are indexed on entity
   * of the same type, same space, approximation base and number of
   * coefficients, are sorted in the way.
   *
   */
  inline std::array<ApproximationOrder, MAX_DOFS_ON_ENTITY> &
  getDofOrderMap(const EntityType type) const {
    return getFieldRawPtr()->getDofOrderMap(type);
  }

  inline const Field *getFieldRawPtr() const {
    return sFieldPtr->getFieldRawPtr();
  };

  inline FieldOrderTable &getFieldOrderTable() {
    return sFieldPtr->getFieldOrderTable();
  };

private:
  mutable boost::shared_ptr<FIELD> sFieldPtr;
};

template <typename T>
struct interface_Field<T, T> : public interface_FieldImpl<T, T> {
  interface_Field(const boost::shared_ptr<T> &ptr)
      : interface_FieldImpl<T, T>(ptr, ptr) {}

  using interface_type_FieldImpl = interface_FieldImpl<T,T>;

  inline EntityHandle getMeshset() const {
    return getFieldRawPtr()->getMeshset();
  }

  /// @return get field Id
  inline const BitFieldId &getId() const {
    return getFieldRawPtr()->getId();
  }

  /// @return get field name
  inline boost::string_ref getNameRef() const {
    return getFieldRawPtr()->getNameRef();
  }

  /// @return get field name
  inline std::string getName() const {
    return getFieldRawPtr()->getName();
  }

  /// @return get approximation space
  inline FieldSpace getSpace() const {
    return getFieldRawPtr()->getSpace();
  }

  /// @return get approximation base
  inline auto getSpaceName() const { return getFieldRawPtr()->getSpaceName(); }

  /// @return get approximation base
  inline FieldApproximationBase getApproxBase() const {
    return getFieldRawPtr()->getApproxBase();
  }

  /// @return get approximation base
  inline auto getApproxBaseName() const {
    return getFieldRawPtr()->getApproxBaseName();
  }

  /// @return get number of coefficients for DOF
  inline FieldCoefficientsNumber getNbOfCoeffs() const {
    return getFieldRawPtr()->getNbOfCoeffs();
  }

  /// @return get bit number if filed Id
  inline FieldBitNumber getBitNumber() const {
    return getFieldRawPtr()->getBitNumber();
  }

  /**
   * \brief get hash-map relating dof index on entity with its order
   *
   * Dofs of given field are indexed on entity
   * of the same type, same space, approximation base and number of
   * coefficients, are sorted in the way.
   *
   */
  inline std::array<ApproximationOrder, MAX_DOFS_ON_ENTITY> &
  getDofOrderMap(const EntityType type) const {
    return getFieldRawPtr()->getDofOrderMap(type);
  }

  inline const Field *getFieldRawPtr() const {
    return boost::static_pointer_cast<T>(this->getRefEntityPtr())
        ->getFieldRawPtr();
  };
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

} // namespace MoFEM

#endif // __FIELDMULTIINDICES_HPP__
