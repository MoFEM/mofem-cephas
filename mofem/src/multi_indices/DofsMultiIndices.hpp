/** \file DofsMultiIndices.hpp
 * \ingroup dof_multi_indices
 * \brief Multi-Index contains, data structures for mofem dofs and other
 * low-level functions
 */

#ifndef __DOFSMULTIINDICES_HPP__
#define __DOFSMULTIINDICES_HPP__

namespace MoFEM {

/**
 * \brief keeps information about DOF on the entity
 * \ingroup dof_multi_indices
 */
struct DofEntity : public interface_FieldEntity<FieldEntity> {

  virtual ~DofEntity() = default;

  using interface_type_Field = interface_FieldImpl<FieldEntity, FieldEntity>;
  using interface_type_FieldEntity = interface_FieldEntity<FieldEntity>;
  using interface_type_RefEntity = interface_RefEntity<FieldEntity>;

  static inline ShortId
  getNonNonuniqueShortId(const DofIdx dof,
                         const boost::shared_ptr<FieldEntity> &ent_ptr) {
    return static_cast<ShortId>(dof) |
           (static_cast<ShortId>(ent_ptr->getBitNumber()) << 9);
  }

  DofEntity(const boost::shared_ptr<FieldEntity> &entity_ptr,
            const ApproximationOrder dof_order,
            const FieldCoefficientsNumber dof_rank, const DofIdx dof);

  /// @return get dof index on entity
  inline DofIdx getEntDofIdx() const { return std::abs(dof); }

  /// @return get field data on dof
  inline FieldData &getFieldData() const {
    return const_cast<FieldData &>(
        (*this->sPtr->getEntFieldDataPtr())[getEntDofIdx()]);
  }

  /// @return get entity unique dof id
  inline const UId &getEntLocalUniqueId() const {
    return this->sPtr->getLocalUniqueId();
  }

  /// @return get entity unique dof id
  inline const UId getEntGlobalUniqueId() const {
    return this->sPtr->getGlobalUniqueId();
  }

  static inline UId getUniqueIdCalculate(const DofIdx dof, UId ent_uid) {
    return ent_uid | dof;
  }

  static inline UId
  getLocalUniqueIdCalculate(const DofIdx dof,
                            const boost::shared_ptr<FieldEntity> &ent_ptr) {
    return getUniqueIdCalculate(dof, ent_ptr->getLocalUniqueId());
  }

  /// @return get unique dof id
  inline UId getLocalUniqueId() const {
    return getUniqueIdCalculate(std::abs(dof), this->sPtr->getLocalUniqueId());
  }

  static inline UId getLoFieldEntityUId(const FieldBitNumber bit,
                                        const EntityHandle ent) {
    return FieldEntity::getLoFieldEntityUId(
        FieldEntity::getLocalUniqueIdCalculate(bit, ent));
  }

  static inline UId getHiFieldEntityUId(const FieldBitNumber bit,
                                        const EntityHandle ent) {
    return FieldEntity::getHiFieldEntityUId(
        FieldEntity::getLocalUniqueIdCalculate(bit, ent));
  }

  /**
   * \brief Calculate UId for DOF
   *
   * UId is constructed such that all DOFs are ordered by processor, entity,
   * field and dof index on entity, On entity dofs index is constructed such
   * that coefficient number and dofs increase with dofs index on entity.
   *
   * @param  dof     dof index on entity
   * @param  ent_ptr pointer to field entity
   * @return         UId
   */
  static inline UId
  getGlobalUniqueIdCalculate(const DofIdx dof,
                             const boost::shared_ptr<FieldEntity> &ent_ptr) {
    return getUniqueIdCalculate(dof, ent_ptr->getGlobalUniqueId());
  }

  /// @return get unique dof id
  inline UId getGlobalUniqueId() const {
    return getUniqueIdCalculate(std::abs(dof), this->sPtr->getGlobalUniqueId());
  }

  /** \brief get short uid it is unique in combination with entity handle
   *
   * EntityHandle are controlled by MOAB, which is unique in
   * MOAB instance. However two MOAB instances, can have attached different
   * EntityHandles to the same entity.
   *
   * Relation between MoAB EntityHandle can be handled by saving entity handle
   * data into tag, see MB_TYPE_HANDLE. MOAB at time of file reading or
   * creating new MOAB instance, substitute tag value by approbate entity
   * handle.
   *
   * ShortId is created to handle problems related to saving data series, and
   * reading those data using different MoAB instances.
   *
   */
  inline ShortId getNonNonuniqueShortId() const {
    return getNonNonuniqueShortId(getEntDofIdx(), getFieldEntityPtr());
  }

  /// @return get dof entity handle
  inline EntityHandle getEnt() const { return this->sPtr->getEnt(); }

  /// @return get dof approximation order
  inline ApproximationOrder getDofOrder() const {
    return getDofOrderMap()[getEntDofIdx()];
  }

  /// @return get dof coefficient index
  inline FieldCoefficientsNumber getDofCoeffIdx() const {
    return getEntDofIdx() % getNbOfCoeffs();
  }

  /// @return return true if dof us active
  inline char getActive() const { return dof < 0 ? 0 : 1; }

  friend std::ostream &operator<<(std::ostream &os, const DofEntity &e);

private:
  DofIdx dof;

  friend struct DofEntity_active_change;
};

/**
 * \brief Interface to DofEntity
 *
 * In MoFEM DOFs classes (and Ent and Finite Element classes) are derived by
 * interface, i.e. not class is derived but interface to it.
 *
 * \ingroup dof_multi_indices
 */
template <typename T>
struct interface_DofEntity : public interface_FieldEntity<T> {

  virtual ~interface_DofEntity() = default;

  interface_DofEntity(const boost::shared_ptr<T> &sptr)
      : interface_FieldEntity<T>(sptr) {}

  /// @return return dof unique id
  inline UId getGlobalUniqueId() const {
    return this->sPtr->getGlobalUniqueId();
  }

  /// @return return entity unique id
  inline const UId &getEntGlobalUniqueId() const {
    return this->sPtr->getEntGlobalUniqueId();
  }

  /// @return return dof unique id
  inline UId getLocalUniqueId() const { return this->sPtr->getLocalUniqueId(); }

  /// @return return entity unique id
  inline const UId &getEntLocalUniqueId() const {
    return this->sPtr->getEntLocalUniqueId();
  }

  /// @return return short id (used by data recorder)
  inline ShortId getNonNonuniqueShortId() const {
    return this->sPtr->getNonNonuniqueShortId();
  }

  /// @return return index of dof on the entity
  inline DofIdx getEntDofIdx() const { return this->sPtr->getEntDofIdx(); }

  /// @return return data on dof
  inline FieldData &getFieldData() const { return this->sPtr->getFieldData(); }

  /// @return return entity handle
  inline EntityHandle getEnt() const { return this->sPtr->getEnt(); }

  /// @return get dof approximation order
  inline ApproximationOrder getDofOrder() const {
    return this->sPtr->getDofOrder();
  }

  /// @return get dof coefficient index
  inline FieldCoefficientsNumber getDofCoeffIdx() const {
    return this->sPtr->getDofCoeffIdx();
  }

  /// @return return true if dof is active
  inline char getActive() const { return this->sPtr->getActive(); }

  /// @return get pointer to dof data structure
  inline boost::shared_ptr<DofEntity> &getDofEntityPtr() const {
    return this->sPtr;
  }

  /// @return get pioneer do dof's entity data structure
  inline boost::shared_ptr<FieldEntity> &getFieldEntityPtr() const {
    return this->sPtr->getFieldEntityPtr();
  }
};

/**
 * \brief keeps information about indexed dofs for the problem
 * \ingroup dof_multi_indices
 *
 * FIXME: Is too many iterator, this has to be manage more efficiently, some
 * iterators could be moved to multi_indices views.
 *
 */
struct NumeredDofEntity : public interface_DofEntity<DofEntity> {

  virtual ~NumeredDofEntity() = default;

  using interface_type_Field = interface_FieldImpl<DofEntity, DofEntity>;
  using interface_type_FieldEntity = interface_FieldEntity<DofEntity>;
  using interface_type_DofEntity = interface_DofEntity<DofEntity>;
  using interface_type_RefEntity = interface_RefEntity<DofEntity>;

  DofIdx dofIdx;
  DofIdx petscGloablDofIdx;
  DofIdx petscLocalDofIdx;
  unsigned int pArt;

  /// @return MoFEM DoF index
  inline DofIdx getDofIdx() const { return dofIdx; }

  /// @return PETSc global DoF index
  inline DofIdx getPetscGlobalDofIdx() const { return petscGloablDofIdx; }

  /// @return PETSc local DoF index
  inline DofIdx getPetscLocalDofIdx() const { return petscLocalDofIdx; }

  /// @return Owning partition (i.e. process/processor)
  inline unsigned int getPart() const { return pArt; }

  /// @return True if local index is set
  inline bool getHasLocalIndex() const {
    return !std::signbit(petscLocalDofIdx);
  }

  NumeredDofEntity(const boost::shared_ptr<DofEntity> &dof_entity_ptr,
                   const int dof_idx = -1, const int petsc_gloabl_dof_idx = -1,
                   const int petsc_local_dof_idx = -1, const int part = -1);
  friend std::ostream &operator<<(std::ostream &os, const NumeredDofEntity &e);
};

/**
 * \brief keeps information about indexed dofs for the finite element
 * \ingroup dof_multi_indices
 */
struct FEDofEntity : public DofEntity {

  FEDofEntity() = delete;

private:
  using DofEntity::DofEntity;

  // DO NOT MAKE ANY MAMBER DATA HARE !!!

  friend std::ostream &operator<<(std::ostream &os, const FEDofEntity &e);
};

/**
 * \brief keeps information about indexed dofs for the finite element
 * \ingroup dof_multi_indices
 */
struct FENumeredDofEntity : public NumeredDofEntity {

  FENumeredDofEntity() = delete;

private:
  using NumeredDofEntity::NumeredDofEntity;

  // DO NOT MAKE ANY MAMBER DATA HARE !!!

  friend std::ostream &operator<<(std::ostream &os,
                                  const FENumeredDofEntity &e);
};

std::ostream &operator<<(std::ostream &os, const FENumeredDofEntity &e);

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps DofEntity
 * \ingroup dof_multi_indices
 */
typedef multi_index_container<
    boost::shared_ptr<DofEntity>,
    indexed_by<
        // unique
        ordered_unique<
            tag<Unique_mi_tag>,
            const_mem_fun<DofEntity, UId, &DofEntity::getLocalUniqueId>>,

        // non_unique
        ordered_non_unique<
            tag<Ent_mi_tag>,
            const_mem_fun<DofEntity, EntityHandle, &DofEntity::getEnt>>

        >>
    DofEntity_multiIndex;

using DofEntityByUId = DofEntity_multiIndex::index<Unique_mi_tag>::type;

/** \brief Dof multi-index by entity
 *
 * \ingroup dof_multi_indices
 */
using DofEntityByEnt = DofEntity_multiIndex::index<Ent_mi_tag>::type;

/** \brief multi-index view on DofEntity by uid
  \ingroup dof_multi_indices
*/
using DofEntity_multiIndex_uid_view =
    multi_index_container<boost::shared_ptr<DofEntity>,
                          indexed_by<

                              ordered_unique<const_mem_fun<
                                  DofEntity, UId, &DofEntity::getLocalUniqueId>>

                              >>;

/** \brief multi-index view on DofEntity by uid
  \ingroup dof_multi_indices
*/
using DofEntity_multiIndex_global_uid_view = multi_index_container<
    boost::shared_ptr<DofEntity>,
    indexed_by<

        ordered_unique<
            const_mem_fun<DofEntity, UId, &DofEntity::getGlobalUniqueId>>

        >>;

/** \brief vector view on DofEntity by uid
  \ingroup dof_multi_indices
*/
typedef std::vector<boost::weak_ptr<DofEntity>> DofEntity_vector_view;

/** \brief multi-index view on DofEntity activity
  \ingroup dof_multi_indices
*/
typedef multi_index_container<
    boost::shared_ptr<DofEntity>,
    indexed_by<

        ordered_unique<
            const_mem_fun<DofEntity, UId, &DofEntity::getLocalUniqueId>>,

        ordered_non_unique<
            const_mem_fun<DofEntity, char, &DofEntity::getActive>>

        >>
    DofEntity_multiIndex_active_view;

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps FEDofEntity
 * \ingroup dof_multi_indices

 */
typedef multi_index_container<
    boost::shared_ptr<FEDofEntity>,
    indexed_by<

        ordered_unique<tag<Unique_mi_tag>,
                       const_mem_fun<FEDofEntity::DofEntity, UId,
                                     &FEDofEntity::getLocalUniqueId>>,

        ordered_non_unique<tag<Ent_mi_tag>,
                           const_mem_fun<FEDofEntity::DofEntity, EntityHandle,
                                         &FEDofEntity::getEnt>>

        >>
    FEDofEntity_multiIndex;

/** \brief Dof entity multi-index by UId and entity
 *
 * \ingroup dof_multi_indices
 */
using FEDofEntityByUId = FEDofEntity_multiIndex::index<Unique_mi_tag>::type;

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps FENumeredDofEntity

 * \ingroup dof_multi_indices
 */
typedef multi_index_container<
    boost::shared_ptr<FENumeredDofEntity>,
    indexed_by<ordered_unique<
                   tag<Unique_mi_tag>,
                   const_mem_fun<FENumeredDofEntity::interface_type_DofEntity,
                                 UId, &FENumeredDofEntity::getLocalUniqueId>>,
               ordered_non_unique<
                   tag<Ent_mi_tag>,
                   const_mem_fun<FENumeredDofEntity::interface_type_DofEntity,
                                 EntityHandle, &FENumeredDofEntity::getEnt>>

               >>
    FENumeredDofEntity_multiIndex;

/** \brief Dof entity multi-index by UId
 *
 * \ingroup dof_multi_indices
 */
using FENumeredDofEntityByUId =
    FENumeredDofEntity_multiIndex::index<Unique_mi_tag>::type;

/** \brief Numbered DoF multi-index by entity
 *
 * \ingroup dof_multi_indices
 */
using FENumeredDofEntityByEnt =
    FENumeredDofEntity_multiIndex::index<Ent_mi_tag>::type;

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps NumeredDofEntity
 * \ingroup dof_multi_indices
 */
typedef multi_index_container<
    boost::shared_ptr<NumeredDofEntity>,

    // unique
    indexed_by<
        ordered_unique<tag<Unique_mi_tag>,
                       const_mem_fun<NumeredDofEntity::interface_type_DofEntity,
                                     UId, &NumeredDofEntity::getLocalUniqueId>>,

        // non unique
        ordered_non_unique<
            tag<Part_mi_tag>,
            member<NumeredDofEntity, unsigned int, &NumeredDofEntity::pArt>>,

        ordered_non_unique<tag<Idx_mi_tag>, member<NumeredDofEntity, DofIdx,
                                                   &NumeredDofEntity::dofIdx>>,

        ordered_non_unique<tag<PetscGlobalIdx_mi_tag>,
                           member<NumeredDofEntity, DofIdx,
                                  &NumeredDofEntity::petscGloablDofIdx>>,

        ordered_non_unique<tag<PetscLocalIdx_mi_tag>,
                           member<NumeredDofEntity, DofIdx,
                                  &NumeredDofEntity::petscLocalDofIdx>>,

        ordered_non_unique<
            tag<Ent_mi_tag>,
            const_mem_fun<NumeredDofEntity::interface_type_DofEntity,
                          EntityHandle, &NumeredDofEntity::getEnt>>

        >>
    NumeredDofEntity_multiIndex;

/** \brief Numbered DoF multi-index by UId
 *
 * \ingroup dof_multi_indices
 */
using NumeredDofEntityByUId =
    NumeredDofEntity_multiIndex::index<Unique_mi_tag>::type;

/** \brief Numbered DoF multi-index by local index
 *
 * \ingroup dof_multi_indices
 */
using NumeredDofEntityByLocalIdx =
    NumeredDofEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type;

/** \brief Numbered DoF multi-index by entity
 *
 * \ingroup dof_multi_indices
 */
using NumeredDofEntityByEnt =
    NumeredDofEntity_multiIndex::index<Ent_mi_tag>::type;

using NumeredDofEntity_multiIndex_uid_view_ordered =
    multi_index_container<boost::shared_ptr<NumeredDofEntity>,
                          indexed_by<

                              ordered_unique<const_mem_fun<

                                  NumeredDofEntity::interface_type_DofEntity,
                                  UId, &NumeredDofEntity::getLocalUniqueId

                                  >>

                              >>;

using NumeredDofEntity_multiIndex_idx_view_hashed = multi_index_container<
    boost::shared_ptr<NumeredDofEntity>,
    indexed_by<hashed_unique<const_mem_fun<NumeredDofEntity, DofIdx,
                                           &NumeredDofEntity::getDofIdx>>>>;

using NumeredDofEntity_multiIndex_petsc_local_dof_view_ordered_non_unique =
    multi_index_container<boost::shared_ptr<NumeredDofEntity>,
                          indexed_by<ordered_non_unique<const_mem_fun<
                              NumeredDofEntity, DofIdx,
                              &NumeredDofEntity::getPetscLocalDofIdx>>>>;

using NumeredDofEntity_multiIndex_coeff_idx_ordered_non_unique =
    multi_index_container<
        boost::shared_ptr<NumeredDofEntity>,
        indexed_by<ordered_non_unique<const_mem_fun<
            NumeredDofEntity::interface_type_DofEntity, FieldCoefficientsNumber,
            &NumeredDofEntity::getDofCoeffIdx>>>>;

/**
 * Activate or deactivate dofs (multi-index modifier)
 * \ingroup dof_multi_indices
 */
struct DofEntity_active_change {
  bool aCtive;
  DofEntity_active_change(bool active);
  void operator()(boost::shared_ptr<DofEntity> &dof);
};

/**
 * Change part and global pestc index (multi-index modifier)
 * \ingroup dof_multi_indices
 */
struct NumeredDofEntity_part_and_glob_idx_change {
  const unsigned int pArt;
  const DofIdx petscGloablDofIdx;
  NumeredDofEntity_part_and_glob_idx_change(const unsigned int part,
                               const DofIdx petsc_gloabl_dof_idx)
      : pArt(part), petscGloablDofIdx(petsc_gloabl_dof_idx){};
  inline void operator()(boost::shared_ptr<NumeredDofEntity> &dof) const {
    dof->pArt = pArt;
    dof->petscGloablDofIdx = petscGloablDofIdx;
  }
};

struct NumeredDofEntity_part_and_mofem_glob_idx_change {
  const unsigned int pArt;
  const DofIdx mofemDofIdx;
  const DofIdx petscGloablDofIdx;
  NumeredDofEntity_part_and_mofem_glob_idx_change(
      const unsigned int part, const DofIdx mofem_dof_idx,
      const DofIdx petsc_gloabl_dof_idx)
      : pArt(part), mofemDofIdx(mofem_dof_idx),
        petscGloablDofIdx(petsc_gloabl_dof_idx){};
  inline void operator()(boost::shared_ptr<NumeredDofEntity> &dof) const {
    dof->pArt = pArt;
    dof->petscGloablDofIdx = petscGloablDofIdx;
  }
};

struct NumeredDofEntity_part_and_indices_change {
  const unsigned int pArt;
  const DofIdx petscGloablDofIdx;
  const DofIdx petscLocalDofIdx;
  NumeredDofEntity_part_and_indices_change(const unsigned int part,
                                           const DofIdx petsc_gloabl_dof_idx,
                                           const DofIdx petsc_local_dof_idx)
      : pArt(part), petscGloablDofIdx(petsc_gloabl_dof_idx),
        petscLocalDofIdx(petsc_local_dof_idx){};
  inline void operator()(boost::shared_ptr<NumeredDofEntity> &dof) const {
    dof->pArt = pArt;
    dof->petscGloablDofIdx = petscGloablDofIdx;
    dof->petscLocalDofIdx = petscLocalDofIdx;
  }
};

/**
 * Change part and local pestc index (multi-index modifier)
 * \ingroup dof_multi_indices
 */
struct NumeredDofEntity_local_idx_change {
  const DofIdx petscLocalDofIdx;
  NumeredDofEntity_local_idx_change(const DofIdx petsc_local_dof_idx)
      : petscLocalDofIdx(petsc_local_dof_idx){};
  inline void operator()(boost::shared_ptr<NumeredDofEntity> &dof) const {
    dof->petscLocalDofIdx = petscLocalDofIdx;
  }
};

/**
 * Change part and mofem index (multi-index modifier)
 * \ingroup dof_multi_indices
 */
struct NumeredDofEntity_mofem_index_change {
  const DofIdx mofemIdx;
  NumeredDofEntity_mofem_index_change(const DofIdx mofem_idx)
      : mofemIdx(mofem_idx){};
  inline void operator()(boost::shared_ptr<NumeredDofEntity> &dof) const {
    dof->dofIdx = mofemIdx;
  }
};

/**
 * Change part and mofem/pestc global and local index (multi-index modifier)
 * \ingroup dof_multi_indices
 */
struct NumeredDofEntity_part_and_all_indices_change {
  const unsigned int pArt;
  const DofIdx mofemIdx;
  const DofIdx petscGloablDofIdx;
  const DofIdx petscLocalDofIdx;
  NumeredDofEntity_part_and_all_indices_change(
      const unsigned int part, const DofIdx mofem_idx,
      const DofIdx petsc_gloabl_dof_idx, const DofIdx petsc_local_dof_idx)
      : pArt(part), mofemIdx(mofem_idx),
        petscGloablDofIdx(petsc_gloabl_dof_idx),
        petscLocalDofIdx(petsc_local_dof_idx){};
  inline void operator()(boost::shared_ptr<NumeredDofEntity> &dof) const {
    dof->pArt = pArt;
    dof->dofIdx = mofemIdx;
    dof->petscGloablDofIdx = petscGloablDofIdx;
    dof->petscLocalDofIdx = petscLocalDofIdx;
  }
};

} // namespace MoFEM
#endif // __DOFSMULTIINDICES_HPP__

/**
 * \defgroup dof_multi_indices Dofs structures and multi-indices
 * \ingroup mofem
 **/
