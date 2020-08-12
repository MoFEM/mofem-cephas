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

struct EntityCacheDofs;

template <int N, int F>
struct FieldEntityTmp : public FieldEntityTmp<N, F - 1> {

  using FieldEntityTmp<N, F - 1>::FieldEntityTmp;

  virtual const FieldTmp<0, 0> *getFieldRawPtr() const { return sFieldRawPtr; }
  static FieldTmp<0, 0> *sFieldRawPtr;

private:
  template <int S> friend struct CoreTmp;
};

template <int N>
struct FieldEntityTmp<N, 0>
    : public FieldEntityTmp<N - 1, BITFIELDID_SIZE - 1> {

  using FieldEntityTmp<N - 1, BITFIELDID_SIZE - 1>::FieldEntityTmp;

  virtual const FieldTmp<0, 0> *getFieldRawPtr() const { return sFieldRawPtr; }
  static FieldTmp<0, 0> *sFieldRawPtr;

private:
};

template <int N, int F> FieldTmp<0, 0> *FieldEntityTmp<N, F>::sFieldRawPtr;
template <int N> FieldTmp<0, 0> *FieldEntityTmp<N, 0>::sFieldRawPtr;

/**
 * \brief Struct keeps handle to entity in the field.
 * \ingroup ent_multi_indices
 */
template <>
struct FieldEntityTmp<0, 0>
    : public interface_FieldImpl<FieldTmp<0, 0>, RefEntity> {

  using interface_type_Field = interface_FieldImpl<FieldTmp<0, 0>, RefEntity>;
  using interface_type_RefEntity = interface_RefEntity<RefEntity>;

  UId localUId; ///< Global unique id for this entity

  FieldEntityTmp(const boost::shared_ptr<FieldTmp<0, 0>> field_ptr,
                 const boost::shared_ptr<RefEntity> ref_ents_ptr,
                 boost::shared_ptr<double *const> field_data_adaptor_ptr,
                 boost::shared_ptr<const int> t_max_order_ptr,
                 boost::shared_ptr<EntityCacheDofs *> ents_cache_dofs);

  virtual ~FieldEntityTmp() = default;

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
    return (this->getFieldRawPtr()->forderTable[this->getEntType()])(order);
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
   * @brief Get the Local Unique Id Calculate
   *
   * @param bit_number
   * @param handle
   * @return UId
   */
  static inline UId getLocalUniqueIdCalculate(const char bit_number,
                                              const EntityHandle handle) {
    constexpr int dof_shift = 9;   // Maximal number of DOFs on entity
    constexpr int ent_shift = 64;  // EntityHandle size
    constexpr int proc_shift = 10; // Maximal number of 1024 processors
    return

        (static_cast<UId>(handle) << proc_shift | static_cast<UId>(bit_number)
                                                      << proc_shift + ent_shift)
        << dof_shift;
  }

  /**
   * @brief Get the Local Unique Id Calculate object
   *
   * @return UId
   */
  inline UId getLocalUniqueIdCalculate() {
    return getLocalUniqueIdCalculate(this->getBitNumber(), this->getEnt());
  }

  /**
   * \brief Get global unique id
   * @return Global UId
   */
  const UId getGlobalUniqueId() const { return getGlobalUniqueIdCalculate(); }

  /**
   * \brief Get global unique id
   * @return Global UId
   */
  const UId &getLocalUniqueId() const { return localUId; }

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
   * meshes
   * @return                          UId
   */
  static inline UId
  getGlobalUniqueIdCalculate(const int owner_proc, const char bit_number,
                             const EntityHandle moab_owner_handle) {
    constexpr int dof_shift = 9;   // Maximal number of DOFs on entity
    constexpr int ent_shift = 64;  // EntityHandle size
    constexpr int proc_shift = 10; // Maximal number of 1024 processors
    return

        (static_cast<UId>(owner_proc) |
         static_cast<UId>(moab_owner_handle) << proc_shift |
         static_cast<UId>(bit_number) << proc_shift + ent_shift)
        << dof_shift;
  }

  /**
   * \brief Calculate global UId
   * @return Global UId
   */
  inline UId getGlobalUniqueIdCalculate() const {
    return getGlobalUniqueIdCalculate(this->getRefEntityPtr()->getOwnerProc(),
                                      this->getBitNumber(),
                                      this->getRefEntityPtr()->getOwnerEnt());
  }

  static inline UId getLoBitNumberUId(const FieldBitNumber bit_number) {
    constexpr int dof_shift = 9;   // Maximal number of DOFs on entity
    constexpr int ent_shift = 64;  // EntityHandle size
    constexpr int proc_shift = 10; // Maximal number of 1024 processors
    return

        static_cast<UId>(bit_number) << dof_shift + ent_shift + proc_shift;
  }

  static inline UId getHiBitNumberUId(const FieldBitNumber bit_number) {
    constexpr int dof_shift = 9;   // Maximal number of DOFs on entity
    constexpr int ent_shift = 64;  // EntityHandle size
    constexpr int proc_shift = 10; // Maximal number of 1024 processors

    return static_cast<UId>(MAX_DOFS_ON_ENTITY - 1) |
           static_cast<UId>(MAX_PROCESSORS_NUMBER - 1) << dof_shift |
           static_cast<UId>(std::numeric_limits<EntityHandle>::max())
               << dof_shift + proc_shift |
           static_cast<UId>(bit_number) << dof_shift + ent_shift + proc_shift;
  }

  static inline UId getLoFieldEntityUId(const UId &uid) {
    return ((~static_cast<UId>(MAX_DOFS_ON_ENTITY - 1)) & uid);
  }

  static inline UId getHiFieldEntityUId(const UId &uid) {
    return getLoFieldEntityUId(uid) | static_cast<UId>(MAX_DOFS_ON_ENTITY - 1);
  }

  template <typename T> void getLoFieldEntityUId(T &uid) = delete;
  template <typename T> void getHiFieldEntityUId(T &uid) = delete;

  static inline UId getLoLocalEntityBitNumber(const char bit_number,
                                              const EntityHandle ent) {
    return getLocalUniqueIdCalculate(

        bit_number,

        ent

    );
  }

  static inline UId getHiLocalEntityBitNumber(const char bit_number,
                                              const EntityHandle ent) {
    return getLoLocalEntityBitNumber(bit_number, ent) |
           static_cast<UId>(MAX_DOFS_ON_ENTITY - 1);
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
    return this->getFieldRawPtr()->getDofOrderMap(this->getEntType());
  }

  friend std::ostream &operator<<(std::ostream &os, const FieldEntity &e);

  virtual const FieldTmp<0, 0> *getFieldRawPtr() const { return sFieldRawPtr; }
  static FieldTmp<0, 0> *sFieldRawPtr;

private:
  mutable boost::shared_ptr<const ApproximationOrder> tagMaxOrderPtr;
  mutable boost::shared_ptr<FieldData *const> fieldDataAdaptorPtr;

  friend struct EntFiniteElement;
  mutable boost::shared_ptr<EntityCacheDofs *> entityCacheDofs;
};

template <> struct FieldEntityTmp<-1, -1> : public FieldEntityTmp<0, 0> {

  FieldEntityTmp(const boost::shared_ptr<FieldTmp<-1, -1>> field_ptr,
                 const boost::shared_ptr<RefEntity> ref_ents_ptr,
                 boost::shared_ptr<double *const> field_data_adaptor_ptr,
                 boost::shared_ptr<const int> t_max_order_ptr,
                 boost::shared_ptr<EntityCacheDofs *> ents_cache_dofs);

  virtual const FieldTmp<0, 0> *getFieldRawPtr() const {
    return sFieldPtr.get();
  }

private:
  template <int S> friend struct CoreTmp;
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

  /// @return get entity UId
  inline UId &getLocalUniqueId() const {
    return this->sPtr->getLocalUniqueId();
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
                       member<FieldEntity, UId, &FieldEntity::localUId>>,
        ordered_non_unique<tag<Ent_mi_tag>,
                           const_mem_fun<FieldEntity::interface_type_RefEntity,
                                         EntityHandle, &FieldEntity::getEnt>>

        >>

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

        ordered_non_unique<tag<Ent_mi_tag>,
                           const_mem_fun<FieldEntity::interface_type_RefEntity,
                                         EntityHandle, &FieldEntity::getEnt>>

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
