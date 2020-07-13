/** \file FieldEntsMultiIndices.hpp
 * \brief Multi-index contains, for mofem entities data structures and other
 * low-level functions
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

#ifndef __FIELD_ENTSMULTIINDICES_HPP__
#define __FIELD_ENTSMULTIINDICES_HPP__

namespace MoFEM {

struct DofEntity;

template <int N, int F>
struct FieldEntityTmp : public FieldEntityTmp<N, F - 1> {

  using FieldEntityTmp<N, F - 1>::FieldEntityTmp;

  virtual boost::shared_ptr<const FieldTmp<0, 0>> &getFieldPtr() const {
    return sFieldPtr;
  }

  static boost::shared_ptr<const FieldTmp<0, 0>> sFieldPtr;
};

template <int N>
struct FieldEntityTmp<N, 0>
    : public FieldEntityTmp<N - 1, BITFIELDID_SIZE - 1> {

  using FieldEntityTmp<N - 1, BITFIELDID_SIZE - 1>::FieldEntityTmp;

  virtual boost::shared_ptr<const FieldTmp<0, 0>> &getFieldPtr() const {
    return sFieldPtr;
  }

  static boost::shared_ptr<const FieldTmp<0, 0>> sFieldPtr;
};

template <int N, int F>
boost::shared_ptr<const FieldTmp<0, 0>> FieldEntityTmp<N, F>::sFieldPtr;

template <int N>
boost::shared_ptr<const FieldTmp<0, 0>> FieldEntityTmp<N, 0>::sFieldPtr;

/**
 * \brief Struct keeps handle to entity in the field.
 * \ingroup ent_multi_indices
 */
template <>
struct FieldEntityTmp<0, 0>
    : public interface_FieldImpl<FieldTmp<0, 0>, RefEntity> {

  using interface_type_Field = interface_FieldImpl<FieldTmp<0, 0>, RefEntity>;
  using interface_type_RefEntity = interface_RefEntity<RefEntity>;

  UId globalUId; ///< Global unique id for this entity

  FieldEntityTmp(const boost::shared_ptr<FieldTmp<0, 0>> field_ptr,
                 const boost::shared_ptr<RefEntity> ref_ents_ptr,
                 boost::shared_ptr<double *const> field_data_adaptor_ptr,
                 boost::shared_ptr<const int> t_max_order_ptr);

  virtual ~FieldEntityTmp() = default;

  virtual boost::shared_ptr<const FieldTmp<0, 0>> &getFieldPtr() const {
    return sFieldPtr;
  }

  /**
   * \brief Get entity handle
   * @return EntityHandle
   */
  inline EntityHandle getEnt() const { return this->getRefEnt(); }

  /**
   * \brief Get number of active DOFs on entity
   * @return Number of DOFs
   */
  inline int getNbDofsOnEnt() const {
    return getOrderNbDofs(getMaxOrder()) * this->getNbOfCoeffs();
  }

  /**
   * @brief Return shared pointer to entity field data vector adaptor
   *
   * @return boost::shared_ptr<VectorAdaptor>
   */
  static boost::shared_ptr<FieldData *const> makeSharedFieldDataAdaptorPtr(
      const boost::shared_ptr<Field> &field_ptr,
      const boost::shared_ptr<RefEntity> &ref_ents_ptr);

  /**
   * @brief Get shared ptr to vector adaptor pointing to the field tag data on
   * entity
   *
   * @return boost::shared_ptr<VectorAdaptor>&
   */
  inline boost::shared_ptr<FieldData *const> &getEntFieldDataPtr() const {
    return fieldDataAdaptorPtr;
  }

  /**
   * \brief Get vector of DOFs active values on entity
   * @return Vector of DOFs values
   */
  inline VectorAdaptor getEntFieldData() const {
    return getVectorAdaptor(*fieldDataAdaptorPtr, getNbDofsOnEnt());
  }

  /**
   * \brief Get number of DOFs on entity for given order of approximation
   * @param  order Order of approximation
   * @return       Number of DOFs
   */
  inline int getOrderNbDofs(ApproximationOrder order) const {
    return (this->getFieldPtr()->forderTable[this->getEntType()])(order);
  }

  /**
   * \brief Get difference of number of DOFs between order and order-1
   * @param  order Approximation order
   * @return       Difference number of DOFs
   */
  inline int getOrderNbDofsDiff(ApproximationOrder order) const {
    return getOrderNbDofs(order) - getOrderNbDofs(order - 1);
  }

  /**
   * \brief Get pinter to Tag keeping approximation order
   * @return Pointer to Tag
   */
  inline const ApproximationOrder *getMaxOrderPtr() const {
    return tagMaxOrderPtr.get();
  }

  /**
   * \brief Get order set to the entity (Allocated tag size for such number)
   * @return Approximation order
   */
  inline ApproximationOrder getMaxOrder() const {
    return *tagMaxOrderPtr.get();
  }

  /**
   * \brief Get global unique id
   * @return Global UId
   */
  const UId &getGlobalUniqueId() const { return globalUId; }

  /**
   * \brief Calculate UId for field entity
   *
   * UId is constructed such that all DOFs are ordered by processor, entity,
   * field.
   * 
   * UId is 128 bit 
   *
   * @param  owner_proc               owning processor
   * @param  bit_number               field bit number
   * @param  moab_owner_handle        entity handle on owning processor
   * @param  true_if_distributed_mesh if true UId is constructed for distributed
   * meshes
   * @return                          UId
   */
  static inline UId
  getGlobalUniqueIdCalculate(const int owner_proc, const char bit_number,
                             const EntityHandle moab_owner_handle,
                             const bool true_if_distributed_mesh) {
    constexpr int dof_shift = 9; // Maximal number of DOFs on entity
    constexpr int ent_shift = 64;  // EntityHandle size
    constexpr int proc_shift = 10; // Maximal number of 1024 processors
    if (true_if_distributed_mesh)
      return

          (static_cast<UId>(moab_owner_handle) |
           static_cast<UId>(owner_proc) << ent_shift |
           static_cast<UId>(bit_number) << proc_shift + ent_shift)
          << dof_shift;
    else
      return

          (static_cast<UId>(moab_owner_handle) | static_cast<UId>(bit_number)
                                                     << proc_shift + ent_shift)
          << dof_shift;
  }

  static inline UId getLoBitNumberUId(const char bit_number) {
    constexpr int dof_shift = 9;   // Maximal number of DOFs on entity
    constexpr int ent_shift = 64;  // EntityHandle size
    constexpr int proc_shift = 10; // Maximal number of 1024 processors
    return

        static_cast<UId>(bit_number) << dof_shift + ent_shift + proc_shift;
  }

  static inline UId getHiBitNumberUId(const char bit_number) {
    constexpr int dof_shift = 9;   // Maximal number of DOFs on entity
    constexpr int ent_shift = 64;  // EntityHandle size
    constexpr int proc_shift = 10; // Maximal number of 1024 processors

    return static_cast<UId>(MAX_DOFS_ON_ENTITY - 1) |
           static_cast<UId>(std::numeric_limits<EntityHandle>::max())
               << dof_shift |
           static_cast<UId>(MAX_PROCESSORS_NUMBER - 1)
               << dof_shift + ent_shift |
           static_cast<UId>(bit_number) << dof_shift + ent_shift + proc_shift;
  }

  static inline UId getLoEntBitNumberUId(const UId &uid) {
    return ((~static_cast<UId>(MAX_DOFS_ON_ENTITY - 1)) & uid);
  }

  static inline UId getHiEntBitNumberUId(const UId &uid) {
    return getLoEntBitNumberUId(uid) | static_cast<UId>(MAX_DOFS_ON_ENTITY - 1);
  }

  /**
   * \brief Calculate global UId
   * @return Global UId
   */
  inline UId getGlobalUniqueIdCalculate() const {
    return getGlobalUniqueIdCalculate(
        this->getRefEntityPtr()->getOwnerProc(), this->getBitNumber(),
        this->getRefEntityPtr()->getOwnerEnt(),
        this->getBasicDataPtr()->trueIfDistributedMesh());
  }

  /**
   * \brief get hash-map relating dof index on entity with its order
   *
   * DOFs of given field are indexed on entity
   * of the same type, same space, approximation base and number of
   * coefficients, are sorted in the way.
   *
   */
  inline std::array<int, MAX_DOFS_ON_ENTITY> &getDofOrderMap() const {
    return this->getFieldPtr()->getDofOrderMap(this->getEntType());
  }

  friend std::ostream &operator<<(std::ostream &os, const FieldEntity &e);

  static boost::shared_ptr<const FieldTmp<0, 0>> sFieldPtr;

private:
  mutable boost::shared_ptr<const ApproximationOrder> tagMaxOrderPtr;
  mutable boost::shared_ptr<FieldData *const> fieldDataAdaptorPtr;
};

template <> struct FieldEntityTmp<-1, -1> : public FieldEntityTmp<0, 0> {

  FieldEntityTmp(const boost::shared_ptr<FieldTmp<-1, -1>> field_ptr,
                 const boost::shared_ptr<RefEntity> ref_ents_ptr,
                 boost::shared_ptr<double *const> field_data_adaptor_ptr,
                 boost::shared_ptr<const int> t_max_order_ptr);

  virtual boost::shared_ptr<const FieldTmp<0, 0>> &getFieldPtr() const {
    return sFieldPtr;
  }

  mutable boost::shared_ptr<const FieldTmp<0, 0>> sFieldPtr;
};

using FieldEntity = FieldEntityTmp<0, 0>;

/**
 * \brief Interface to FieldEntity
 * \ingroup ent_multi_indices
 *
 * interface to FieldEntity
 */
template <typename T>
struct interface_FieldEntity : public interface_Field<T, T> {

  interface_FieldEntity(const boost::shared_ptr<T> &sptr)
      : interface_Field<T, T>(sptr) {}

  /// @return get entity handle
  inline EntityHandle getEnt() const { return this->sPtr->getEnt(); }

  /// @return get number of dofs on entity
  inline int getNbDofsOnEnt() const { return this->sPtr->getNbDofsOnEnt(); }

  /// @return get field data on entity
  inline VectorAdaptor getEntFieldData() const {
    return this->sPtr->getEntFieldData();
  }

  /// @return get number of DOFs for given order
  inline int getOrderNbDofs(ApproximationOrder order) const {
    return this->sPtr->getOrderNbDofs(order);
  }

  /// @return get increase of DOFs by increase to this order
  inline int getOrderNbDofsDiff(ApproximationOrder order) const {
    return this->sPtr->getOrderNbDofsDiff(order);
  }

  /// @return get maximal order on entity
  inline ApproximationOrder getMaxOrder() const {
    return this->sPtr->getMaxOrder();
  }

  /// @return get entity UId
  inline UId getGlobalUniqueId() const {
    return this->sPtr->getGlobalUniqueId();
  }

  /// @return return pointer to reference entity data structure
  inline boost::shared_ptr<RefEntity> &getRefEntityPtr() const {
    return this->sPtr->getRefEntityPtr();
  }

  /// @return get pointer to mofem entity data structure
  inline boost::shared_ptr<FieldEntity> &getFieldEntityPtr() const {
    return this->sPtr;
  };

  /**
   * \brief get hash-map relating dof index on entity with its order
   *
   * DOFs of given field are indexed on entity
   * of the same type, same space, approximation base and number of
   * coefficients, are sorted in the way.
   *
   */
  inline std::array<int, MAX_DOFS_ON_ENTITY> &getDofOrderMap() const {
    return this->sPtr->getDofOrderMap();
  }
};

/**
 * \brief structure to change FieldEntity order
 * \ingroup ent_multi_indices
 */
struct FieldEntity_change_order {

  FieldEntity_change_order(const ApproximationOrder order,
                           const bool reduce_tag_size = false)
      : order(order), reduceTagSize(reduce_tag_size) {}
  inline void operator()(boost::shared_ptr<FieldEntity> &e) {
    (*this)(e.get());
  }
  void operator()(FieldEntity *e);

private:
  const ApproximationOrder order;
  const bool reduceTagSize;
  std::vector<FieldData> data;
};

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps FieldEntity
 * \ingroup ent_multi_indices
 *
 */
typedef multi_index_container<
    boost::shared_ptr<FieldEntity>,
    indexed_by<
        ordered_unique<tag<Unique_mi_tag>,
                       member<FieldEntity, UId, &FieldEntity::globalUId>>,
        ordered_non_unique<
            tag<Ent_mi_tag>,
            const_mem_fun<FieldEntity, EntityHandle, &FieldEntity::getEnt>>,
        ordered_non_unique<
            tag<Composite_Name_And_Ent_mi_tag>,
            composite_key<
                FieldEntity,
                const_mem_fun<FieldEntity::interface_type_Field,
                              boost::string_ref, &FieldEntity::getNameRef>,
                const_mem_fun<FieldEntity, EntityHandle,
                              &FieldEntity::getEnt>>>>>
    FieldEntity_multiIndex;

/** \brief Entity index by field name
 *
 * \ingroup ent_multi_indices
 */
using FieldEntityByUId = FieldEntity_multiIndex::index<Unique_mi_tag>::type;

typedef multi_index_container<
    boost::shared_ptr<FieldEntity>,
    indexed_by<

        sequenced<>,

        ordered_non_unique<
            tag<Ent_mi_tag>,
            const_mem_fun<FieldEntity, EntityHandle, &FieldEntity::getEnt>>

        >>
    FieldEntity_multiIndex_ent_view;

typedef multi_index_container<
    boost::shared_ptr<FieldEntity>,
    indexed_by<

        sequenced<>,

        ordered_non_unique<
            tag<Composite_EntType_and_Space_mi_tag>,
            composite_key<FieldEntity,

                          const_mem_fun<FieldEntity::interface_type_RefEntity,
                                        EntityType, &FieldEntity::getEntType>,

                          const_mem_fun<FieldEntity::interface_type_Field,
                                        FieldSpace, &FieldEntity::getSpace>

                          >>

        >>
    FieldEntity_multiIndex_spaceType_view;

typedef std::vector<boost::weak_ptr<FieldEntity>> FieldEntity_vector_view;

} // namespace MoFEM

#endif // __FIELD_ENTSMULTIINDICES_HPP__

/**
 * \defgroup ent_multi_indices Entities structures and multi-indices
 * \ingroup mofem
 **/
