/** \file EntsMultiIndices.hpp
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

#ifndef __ENTSMULTIINDICES_HPP__
#define __ENTSMULTIINDICES_HPP__

namespace MoFEM {

/**
 * @brief Get the tag ptr object
 *
 * @param moab
 * @param th
 * @param ent
 * @param tag_size
 * @return void*
 */
inline void *get_tag_ptr(moab::Interface &moab, Tag th, EntityHandle ent,
                         int *tag_size) {
  void *ret_val;
  rval = moab.tag_get_by_ptr(th, &ent, 1, (const void **)&ret_val, tag_size);
  if (rval != MB_SUCCESS) {
    *tag_size = 0;
    return NULL;
  } else {
    return ret_val;
  }
}

/**
 * \brief keeps information about side number for the finite element
 * \ingroup ent_multi_indices
 */
struct __attribute__((__packed__)) SideNumber {
  EntityHandle ent;
  char side_number;
  char sense;
  char offset;
  char brother_side_number;
  inline EntityType getEntType() const {
    return (EntityType)((ent & MB_TYPE_MASK) >> MB_ID_WIDTH);
  }

  SideNumber(EntityHandle _ent, int _side_number, int _sense, int _offset)
      : ent(_ent), side_number(_side_number), sense(_sense), offset(_offset),
        brother_side_number(-1) {}
  virtual ~SideNumber() = default;
};

/**
 * @relates multi_index_container
 * \brief SideNumber_multiIndex for SideNumber
 * \ingroup ent_multi_indices
 *
 */
typedef multi_index_container<
    boost::shared_ptr<SideNumber>,
    indexed_by<
        hashed_unique<member<SideNumber, EntityHandle, &SideNumber::ent>>,
        ordered_non_unique<

            composite_key<
                SideNumber,
                const_mem_fun<SideNumber, EntityType, &SideNumber::getEntType>,
                member<SideNumber, char, &SideNumber::side_number>

                >>,
        ordered_non_unique<
            const_mem_fun<SideNumber, EntityType, &SideNumber::getEntType>>

        >>
    SideNumber_multiIndex;

/**
 * \brief PipelineManager data. like access to moab interface and basic tag
 * handlers.
 */
struct BasicEntityData {
  moab::Interface &moab;
  int pcommID;
  Tag th_RefParentHandle;
  Tag th_RefBitLevel;
  BasicEntityData(const moab::Interface &mfield,
                  const int pcomm_id = MYPCOMM_INDEX);
  virtual ~BasicEntityData() = default;
  inline void setDistributedMesh() { distributedMesh = true; }
  inline void unSetDistributedMesh() { distributedMesh = false; }
  inline bool trueIfDistributedMesh() const { return distributedMesh; }

private:
  bool distributedMesh;
};

template <int N> struct RefEntityTmp : public RefEntityTmp<N - 1> {

  using RefEntityTmp<N - 1>::RefEntityTmp;

  virtual boost::shared_ptr<BasicEntityData> &getBasicDataPtr() {
    return basicDataPtr;
  }

  virtual const boost::shared_ptr<BasicEntityData> &getBasicDataPtr() const {
    return basicDataPtr;
  }

  static boost::shared_ptr<BasicEntityData> basicDataPtr;
};

template <int N>
boost::shared_ptr<BasicEntityData> RefEntityTmp<N>::basicDataPtr;

/**
 * \brief Struct keeps handle to refined handle.
 * \ingroup ent_multi_indices

  \todo th_RefType "_RefType" is set as two integers, need to be fixed, it is
  waste of space.

 */
template <> struct RefEntityTmp<0> {

  RefEntityTmp(const boost::shared_ptr<BasicEntityData> &basic_data_ptr,
               const EntityHandle ent)
      : ent(ent) {}

  virtual ~RefEntityTmp() = default;

  virtual boost::shared_ptr<BasicEntityData> &getBasicDataPtr() {
    return basicDataPtr;
  }

  virtual const boost::shared_ptr<BasicEntityData> &getBasicDataPtr() const {
    return basicDataPtr;
  }

  /**
   * @brief Get the entity handle
   *
   * @return EntityHandle
   */
  inline EntityHandle getRefEnt() const { return this->ent; }

  /** \brief Get entity type
   */
  inline EntityType getEntType() const {
    return (EntityType)((this->ent & MB_TYPE_MASK) >> MB_ID_WIDTH);
  }

  /** \brief get entity id
   */
  inline EntityID getEntId() const {
    return (EntityID)(this->ent & MB_ID_MASK);
  };

  /** \brief Owner handle on this or other processors
   */
  inline EntityHandle getOwnerEnt() const {
    ParallelComm *pcomm = ParallelComm::get_pcomm(
        &this->getBasicDataPtr()->moab, this->getBasicDataPtr()->pcommID);
    auto pstat = *static_cast<unsigned char *>(MoFEM::get_tag_ptr(
        this->getBasicDataPtr()->moab, pcomm->pstatus_tag(), this->ent, NULL));
    if (!(pstat & PSTATUS_NOT_OWNED)) {
      return this->ent;
    } else if (pstat & PSTATUS_MULTISHARED) {
      return static_cast<EntityHandle *>(
          MoFEM::get_tag_ptr(this->getBasicDataPtr()->moab,
                             pcomm->sharedhs_tag(), this->ent, NULL))[0];
    } else if (pstat & PSTATUS_SHARED) {
      return static_cast<EntityHandle *>(
          MoFEM::get_tag_ptr(this->getBasicDataPtr()->moab,
                             pcomm->sharedh_tag(), this->ent, NULL))[0];
    } else {
      return 0;
    }
  }

  /** \brief Get processor owning entity
   */
  inline int getOwnerProc() const {
    ParallelComm *pcomm = ParallelComm::get_pcomm(
        &this->getBasicDataPtr()->moab, this->getBasicDataPtr()->pcommID);
    auto pstat = *static_cast<unsigned char *>(MoFEM::get_tag_ptr(
        this->getBasicDataPtr()->moab, pcomm->pstatus_tag(), this->ent, NULL));
    if (!(pstat & PSTATUS_NOT_OWNED)) {
      return pcomm->rank();
    } else if (pstat & PSTATUS_MULTISHARED) {
      return static_cast<int *>(
          MoFEM::get_tag_ptr(this->getBasicDataPtr()->moab,
                             pcomm->sharedps_tag(), this->ent, NULL))[0];
    } else if (pstat & PSTATUS_SHARED) {
      return static_cast<int *>(
          MoFEM::get_tag_ptr(this->getBasicDataPtr()->moab,
                             pcomm->sharedp_tag(), this->ent, NULL))[0];
    } else {
      return -1;
    }
  }

  /** \brief Get processor
   */
  inline int getPartProc() const {
    ParallelComm *pcomm = ParallelComm::get_pcomm(
        &this->getBasicDataPtr()->moab, this->getBasicDataPtr()->pcommID);
    return *static_cast<int *>(MoFEM::get_tag_ptr(this->getBasicDataPtr()->moab,
                                                  pcomm->partition_tag(),
                                                  this->ent, NULL));
  }

  /** \brief Get processor owning entity
   */
  int &getPartProc() {
    ParallelComm *pcomm = ParallelComm::get_pcomm(
        &this->getBasicDataPtr()->moab, this->getBasicDataPtr()->pcommID);
    return *static_cast<int *>(MoFEM::get_tag_ptr(this->getBasicDataPtr()->moab,
                                                  pcomm->partition_tag(),
                                                  this->ent, NULL));
  }

  /** \brief get pstatus
   * This tag stores various aspects of parallel status in bits; see also
   * define following, to be used in bit mask operations.  If an entity is
   * not shared with any other processors, the pstatus is 0, otherwise it's > 0
   *
   * bit 0: !owned (0=owned, 1=not owned)
   * bit 1: shared (0=not shared, 1=shared)
   * bit 2: multishared (shared by > 2 procs; 0=not shared, 1=shared)
   * bit 3: interface (0=not interface, 1=interface)
   * bit 4: ghost (0=not ghost, 1=ghost)
   *
   */
  inline unsigned char getPStatus() const {
    ParallelComm *pcomm = ParallelComm::get_pcomm(
        &this->getBasicDataPtr()->moab, this->getBasicDataPtr()->pcommID);
    return *static_cast<unsigned char *>(MoFEM::get_tag_ptr(
        this->getBasicDataPtr()->moab, pcomm->pstatus_tag(), this->ent, NULL));
  }

  /** \brief get shared processors

  Returning list to shared processors. Lists end with -1. Returns NULL if not
  sharing processors.

  DO NOT MODIFY LIST.

  \code
    BasicEntity *ent_ptr = BasicEntity(moab,entity_handle);
    for(int proc = 0; proc<MAX_SHARING_PROCS && -1 !=
  ent_ptr->getSharingProcsPtr[proc]; proc++) {
        if(ent_ptr->getSharingProcsPtr[proc] == -1) {
        // End of the list
        break;
        }
        int sharing_proc = ent_ptr->getSharingProcsPtr[proc];
        EntityHandle sharing_ent = ent_ptr->getSharingHandlersPtr[proc];
        if(!(ent_ptr->getPStatus()&PSTATUS_MULTISHARED)) {
        break;
        }
      }
  \endcode

  */
  int *getSharingProcsPtr() const {
    auto &moab = this->getBasicDataPtr()->moab;
    int *sharing_procs_ptr = NULL;
    ParallelComm *pcomm =
        ParallelComm::get_pcomm(&moab, this->getBasicDataPtr()->pcommID);
    if (getPStatus() & PSTATUS_MULTISHARED) {
      // entity is multi shared
      rval = moab.tag_get_by_ptr(pcomm->sharedps_tag(), &this->ent, 1,
                                 (const void **)&sharing_procs_ptr);
      MOAB_THROW(rval);
    } else if (getPStatus() & PSTATUS_SHARED) {
      // shared
      rval = moab.tag_get_by_ptr(pcomm->sharedp_tag(), &this->ent, 1,
                                 (const void **)&sharing_procs_ptr);
      MOAB_THROW(rval);
    }
    return sharing_procs_ptr;
  }

  /** \brief get sharid entity handlers

  Returning list to shared entity handlers. Use it with getSharingProcsPtr()

  DO NOT MODIFY LIST.

  \code
    BasicEntity *ent_ptr = BasicEntity(moab,entity_handle);
    for(int proc = 0; proc<MAX_SHARING_PROCS && -1 !=
     ent_ptr->getSharingProcsPtr[proc]; proc++) {
      if(ent_ptr->getSharingProcsPtr[proc] == -1) {
        // End of the list
        break;
      }
      int sharing_proc = ent_ptr->getSharingProcsPtr[proc];
      EntityHandle sharing_ent = ent_ptr->getSharingHandlersPtr[proc];
      if(!(ent_ptr->getPStatus()&PSTATUS_MULTISHARED)) {
        break;
      }
    }
  \endcode

    */
  inline EntityHandle *getSharingHandlersPtr() const {
    EntityHandle *sharing_handlers_ptr = NULL;
    auto &moab = this->getBasicDataPtr()->moab;
    ParallelComm *pcomm =
        ParallelComm::get_pcomm(&moab, this->getBasicDataPtr()->pcommID);
    if (getPStatus() & PSTATUS_MULTISHARED) {
      // entity is multi shared
      rval = moab.tag_get_by_ptr(pcomm->sharedhs_tag(), &this->ent, 1,
                                 (const void **)&sharing_handlers_ptr);
      MOAB_THROW(rval);
    } else if (getPStatus() & PSTATUS_SHARED) {
      // shared
      rval = moab.tag_get_by_ptr(pcomm->sharedh_tag(), &this->ent, 1,
                                 (const void **)&sharing_handlers_ptr);
      MOAB_THROW(rval);
    }
    return sharing_handlers_ptr;
  }

  static MoFEMErrorCode getParentEnt(Interface &moab, Range ents,
                                     std::vector<EntityHandle> vec_patent_ent);

  static MoFEMErrorCode
  getBitRefLevel(Interface &moab, Range ents,
                 std::vector<BitRefLevel> &vec_bit_ref_level) {
    MoFEMFunctionBegin;
    Tag th_ref_bit_level;
    CHKERR moab.tag_get_handle("_RefBitLevel", th_ref_bit_level);
    vec_bit_ref_level.resize(ents.size());
    CHKERR moab.tag_get_data(th_ref_bit_level, ents,
                             &*vec_bit_ref_level.begin());
    MoFEMFunctionReturn(0);
  }

  static MoFEMErrorCode
  getBitRefLevel(Interface &moab, Range ents,
                 std::vector<const BitRefLevel *> &vec_ptr_bit_ref_level) {
    MoFEMFunctionBegin;
    Tag th_ref_bit_level;
    CHKERR moab.tag_get_handle("_RefBitLevel", th_ref_bit_level);
    vec_ptr_bit_ref_level.resize(ents.size());
    CHKERR moab.tag_get_by_ptr(
        th_ref_bit_level, ents,
        reinterpret_cast<const void **>(&*vec_ptr_bit_ref_level.begin()));
    MoFEMFunctionReturn(0);
  }

  /**
   * \brief Get pointer to parent entity tag.
   *
   * Each refined entity has his parent. Such information is stored on tags.
   * This function get pinter to tag.
   *
   * @return Pointer to tag on entity
   */
  EntityHandle *getParentEntPtr() const {
    return static_cast<EntityHandle *>(get_tag_ptr(
        this->getBasicDataPtr()->moab,
        this->getBasicDataPtr()->th_RefParentHandle, this->ent, NULL));
  }

  /**
   * \brief Get pointer to bit ref level tag

   * Every entity belongs to some refinement level or levels. Each level is
   marked
   * by bit set in BitRefLevel() (bitset) structure.
   *
   * See \ref mix_mesh_refinement for explanation.

   * @return Return pointer to tag.
   */
  BitRefLevel *getBitRefLevelPtr() const {
    return static_cast<BitRefLevel *>(
        get_tag_ptr(this->getBasicDataPtr()->moab,
                    this->getBasicDataPtr()->th_RefBitLevel, this->ent, NULL));
  }

  /** \brief Get patent entity
   */
  inline EntityType getParentEntType() const {
    EntityHandle *tag_parent_ent = getParentEntPtr();
    if (*tag_parent_ent == 0)
      return MBMAXTYPE;
    return (EntityType)((*tag_parent_ent & MB_TYPE_MASK) >> MB_ID_WIDTH);
  }

  /** \brief Get parent entity, i.e. entity form one refinement level up
   */
  inline EntityHandle getParentEnt() const { return *(getParentEntPtr()); }

  /** \brief Get entity ref bit refinement signature
   */
  inline const BitRefLevel &getBitRefLevel() const {
    return *getBitRefLevelPtr();
  }

  /** \brief Get entity ref bit refinement as ulong
   */
  inline unsigned long int getBitRefLevelULong() const {
    return getBitRefLevel().to_ulong();
  }

  friend std::ostream &operator<<(std::ostream &os, const RefEntityTmp &e);

  EntityHandle ent;
  static boost::shared_ptr<BasicEntityData> basicDataPtr;
};

template <> struct RefEntityTmp<-1> : public RefEntityTmp<0> {

  RefEntityTmp(const boost::shared_ptr<BasicEntityData> &basic_data_ptr,
               const EntityHandle ent)
      : RefEntityTmp<0>(basic_data_ptr, ent), basicDataPtr(basic_data_ptr) {}

  virtual boost::shared_ptr<BasicEntityData> &getBasicDataPtr() {
    return basicDataPtr;
  }

  virtual const boost::shared_ptr<BasicEntityData> &getBasicDataPtr() const {
    return basicDataPtr;
  }

private:
  mutable boost::shared_ptr<BasicEntityData> basicDataPtr;
};

using RefEntity = RefEntityTmp<0>;

/**
 * \brief interface to RefEntity
 * \ingroup ent_multi_indices
 */
template <typename T> struct interface_RefEntity {

  mutable boost::shared_ptr<T> sPtr;

  interface_RefEntity(const boost::shared_ptr<T> &sptr) : sPtr(sptr) {}

  interface_RefEntity(const interface_RefEntity<T> &interface)
      : sPtr(interface.getRefEntityPtr()) {}

  virtual ~interface_RefEntity() = default;

  inline boost::shared_ptr<BasicEntityData> &getBasicDataPtr() {
    return this->sPtr->getBasicDataPtr();
  }

  inline const boost::shared_ptr<BasicEntityData> &getBasicDataPtr() const {
    return this->sPtr->getBasicDataPtr();
  }

  inline EntityHandle getRefEnt() const { return this->sPtr->getRefEnt(); }

  inline EntityType getParentEntType() const {
    return this->sPtr->getParentEntType();
  };

  inline EntityHandle getParentEnt() const {
    return this->sPtr->getParentEnt();
  }

  inline BitRefLevel *getBitRefLevelPtr() const {
    return this->sPtr->getBitRefLevelPtr();
  }

  inline const BitRefLevel &getBitRefLevel() const {
    return this->sPtr->getBitRefLevel();
  }

  inline unsigned long int getBitRefLevelULong() const {
    return this->sPtr->getBitRefLevelULong();
  }

  inline EntityType getEntType() const { return this->sPtr->getEntType(); };

  inline EntityID getEntId() const { return this->sPtr->getEntId(); };

  inline EntityHandle getOwnerEnt() const { return this->sPtr->getOwnerEnt(); }

  inline EntityHandle &getOwnerEnt() { return this->sPtr->getOwnerEnt(); }

  inline int getOwnerProc() const { return this->sPtr->getOwnerProc(); }

  inline int getPartProc() const { return this->sPtr->getPartProc(); }

  inline unsigned char getPStatus() const { return this->sPtr->getPStatus(); }

  inline int *getSharingProcsPtr() const {
    return this->sPtr->getSharingProcsPtr();
  }

  inline EntityHandle *getSharingHandlersPtr() const {
    return this->sPtr->getSharingHandlersPtr();
  }

  inline boost::shared_ptr<T> &getRefEntityPtr() const { return this->sPtr; }
};

/**
 * \typedef RefEntity_multiIndex
 * type multiIndex container for RefEntity
 * \ingroup ent_multi_indices
 *
 * \param hashed_unique Ent_mi_tag
 * \param ordered_non_unique Meshset_mi_tag
 * \param ordered_non_unique Ent_Ent_mi_tag
 * \param ordered_non_unique EntType_mi_tag
 * \param ordered_non_unique ParentEntType_mi_tag
 * \param ordered_non_unique Composite_EntType_And_ParentEntType_mi_tag
 * \param ordered_non_unique Composite_ParentEnt_And_EntType_mi_tag
 */
using RefEntity_multiIndex = multi_index_container<
    boost::shared_ptr<RefEntity>,
    indexed_by<
        ordered_unique<

            tag<Ent_mi_tag>,
            const_mem_fun<RefEntity, EntityHandle, &RefEntity::getRefEnt>

            >,

        ordered_non_unique<

            tag<Ent_Ent_mi_tag>,
            const_mem_fun<RefEntity, EntityHandle, &RefEntity::getParentEnt>

            >,

        ordered_non_unique<

            tag<EntType_mi_tag>,
            const_mem_fun<RefEntity, EntityType, &RefEntity::getEntType>

            >,

        ordered_non_unique<

            tag<ParentEntType_mi_tag>,
            const_mem_fun<RefEntity, EntityType, &RefEntity::getParentEntType>

            >,

        ordered_non_unique<
            tag<Composite_EntType_and_ParentEntType_mi_tag>,
            composite_key<
                RefEntity,
                const_mem_fun<RefEntity, EntityType, &RefEntity::getEntType>,
                const_mem_fun<RefEntity, EntityType,
                              &RefEntity::getParentEntType>>

            >,

        ordered_non_unique<
            tag<Composite_ParentEnt_And_EntType_mi_tag>,
            composite_key<
                RefEntity,
                const_mem_fun<RefEntity, EntityType, &RefEntity::getEntType>,
                const_mem_fun<RefEntity, EntityHandle,
                              &RefEntity::getParentEnt>>>>

    >;

/** \brief multi-index view of RefEntity by parent entity
  \ingroup ent_multi_indices
*/
using RefEntity_multiIndex_view_by_hashed_parent_entity = multi_index_container<
    boost::shared_ptr<RefEntity>,
    indexed_by<
        hashed_non_unique<
            const_mem_fun<RefEntity, EntityHandle, &RefEntity::getParentEnt>>,
        hashed_unique<tag<Composite_EntType_and_ParentEntType_mi_tag>,
                      composite_key<boost::shared_ptr<RefEntity>,
                                    const_mem_fun<RefEntity, EntityHandle,
                                                  &RefEntity::getRefEnt>,
                                    const_mem_fun<RefEntity, EntityHandle,
                                                  &RefEntity::getParentEnt>>>>

    >;

using RefEntity_multiIndex_view_by_ordered_parent_entity =
    multi_index_container<
        boost::shared_ptr<RefEntity>,
        indexed_by<ordered_non_unique<const_mem_fun<RefEntity, EntityHandle,
                                                    &RefEntity::getParentEnt>>,
                   hashed_unique<
                       tag<Composite_EntType_and_ParentEntType_mi_tag>,
                       composite_key<boost::shared_ptr<RefEntity>,
                                     const_mem_fun<RefEntity, EntityHandle,
                                                   &RefEntity::getRefEnt>,
                                     const_mem_fun<RefEntity, EntityHandle,
                                                   &RefEntity::getParentEnt>>>>

        >;

using RefEntity_multiIndex_view_sequence_ordered_view = multi_index_container<
    boost::shared_ptr<RefEntity>,
    indexed_by<
        sequenced<>,
        ordered_unique<tag<Ent_mi_tag>, const_mem_fun<RefEntity, EntityHandle,
                                                      &RefEntity::getRefEnt>>>>;

template <class T> struct Entity_update_pcomm_data {
  const int pcommID;
  Entity_update_pcomm_data(const int pcomm_id = MYPCOMM_INDEX)
      : pcommID(pcomm_id) {}
  void operator()(boost::shared_ptr<T> &e) {
    e->getBasicDataPtr()->pcommID = pcommID;
    ParallelComm *pcomm =
        ParallelComm::get_pcomm(&e->getBasicDataPtr()->moab, pcommID);
    if (pcomm == NULL)
      THROW_MESSAGE("pcomm is null");
    if (e->getBasicDataPtr()->trueIfDistributedMesh()) {
      THROW_MESSAGE("Can not change owner proc if distributed mesh, this will "
                    "make undetermined behavior");
    }
    rval = pcomm->get_owner_handle(e->getRefEnt(), e->getOwnerProc(),
                                   e->getOwnerEnt());
    MOAB_THROW(rval);
    EntityHandle ent = e->getRefEnt();
    rval = e->getBasicDataPtr()->moab.tag_get_data(pcomm->part_tag(), &ent, 1,
                                                   &e->getPartProc());
    MOAB_THROW(rval);
  }
};

/** \brief change parent
  * \ingroup ent_multi_indices
  *
  * Use this function with care. Some other multi-indices can deponent on this.

  Known dependent multi-indices (verify if that list is full):
  - RefEntity_multiIndex
  - RefElement_multiIndex

  */
struct RefEntity_change_parent {
  EntityHandle pArent;
  RefEntity_change_parent(EntityHandle parent) : pArent(parent) {}
  inline void operator()(boost::shared_ptr<RefEntity> &e) {
    rval = e->getBasicDataPtr()->moab.tag_set_data(
        e->getBasicDataPtr()->th_RefParentHandle, &e->ent, 1, &pArent);
    MOAB_THROW(rval);
  }
};

/** \brief ref mofem entity, left shift
 * \ingroup ent_multi_indices
 */
struct RefEntity_change_left_shift {
  int shift;
  BitRefLevel mask;
  RefEntity_change_left_shift(const int _shift,
                              const BitRefLevel _mask = BitRefLevel().set())
      : shift(_shift), mask(_mask) {}
  inline void operator()(boost::shared_ptr<RefEntity> &e) {
    BitRefLevel bit = *(e->getBitRefLevelPtr());
    (*e->getBitRefLevelPtr()) = ((bit & mask) << shift) | (bit & ~mask);
  };
};

/** \brief ref mofem entity, right shift
 * \ingroup ent_multi_indices
 */
struct RefEntity_change_right_shift {
  int shift;
  BitRefLevel mask;
  RefEntity_change_right_shift(const int _shift,
                               const BitRefLevel _mask = BitRefLevel().set())
      : shift(_shift), mask(_mask) {}
  inline void operator()(boost::shared_ptr<RefEntity> &e) {
    BitRefLevel bit = *(e->getBitRefLevelPtr());
    *(e->getBitRefLevelPtr()) = ((bit & mask) >> shift) | (bit & ~mask);
  };
};

struct DofEntity;

/**
 * \brief Struct keeps handle to entity in the field.
 * \ingroup ent_multi_indices
 */
struct FieldEntity : public interface_Field<Field>,
                     interface_RefEntity<RefEntity> {

  using interface_type_Field = interface_Field<Field>;
  using interface_type_RefEntity = interface_RefEntity<RefEntity>;
  
  UId globalUId; ///< Global unique id for this entity

  FieldEntity(const boost::shared_ptr<Field> &field_ptr,
              const boost::shared_ptr<RefEntity> &ref_ents_ptr,
              boost::shared_ptr<double *const> &&field_data_adaptor_ptr,
              boost::shared_ptr<const int> &&t_max_order_ptr);

  virtual ~FieldEntity() = default;

  /**
   * \brief Get entity handle
   * @return EntityHandle
   */
  inline EntityHandle getEnt() const { return getRefEnt(); }

  /**
   * \brief Get number of active DOFs on entity
   * @return Number of DOFs
   */
  inline int getNbDofsOnEnt() const {
    return getOrderNbDofs(getMaxOrder()) * getNbOfCoeffs();
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
    return (this->getFieldPtr()->forderTable[getEntType()])(order);
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
    // assert(bit_number < 32);
    // assert(owner_proc < 1024);
    constexpr int ent_shift = 8 * sizeof(EntityHandle);
    if (true_if_distributed_mesh)
      return (static_cast<UId>(moab_owner_handle) |
              static_cast<UId>(bit_number) << ent_shift |
              static_cast<UId>(owner_proc) << 5 + ent_shift)
             << 9;
    else
      return (static_cast<UId>(moab_owner_handle) | static_cast<UId>(bit_number)
                                                        << ent_shift)
             << 9;
  }

  static inline UId getGlobalUniqueIdCalculate_Low_Proc(const int owner_proc) {
    return getGlobalUniqueIdCalculate(owner_proc, 0, 0, true);
  }

  static inline UId getGlobalUniqueIdCalculate_Hi_Proc(const int owner_proc) {
    return getGlobalUniqueIdCalculate(owner_proc, BITFIELDID_SIZE - 1,
                                      MBMAXTYPE, true);
  }

  /**
   * \brief Calculate global UId
   * @return Global UId
   */
  inline UId getGlobalUniqueIdCalculate() const {
    return getGlobalUniqueIdCalculate(
        sPtr->getOwnerProc(), getBitNumber(), sPtr->getOwnerEnt(),
        getBasicDataPtr()->trueIfDistributedMesh());
  }

  /**
   * \brief Get pointer to RefEntity
   */
  inline boost::shared_ptr<RefEntity> &getRefEntityPtr() { return this->sPtr; }

  /**
   * \brief Get pointer to Field data structure associated with this entity
   */
  inline boost::shared_ptr<Field> &getFieldPtr() const {
    return this->sFieldPtr;
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
    return getFieldPtr()->getDofOrderMap(getEntType());
  }

  friend std::ostream &operator<<(std::ostream &os, const FieldEntity &e);

private:
  mutable boost::shared_ptr<const ApproximationOrder> tagMaxOrderPtr;
  mutable boost::shared_ptr<FieldData *const> fieldDataAdaptorPtr;
};

/**
 * \brief Interface to FieldEntity
 * \ingroup ent_multi_indices
 *
 * interface to FieldEntity
 */
template <typename T>
struct interface_FieldEntity : public interface_Field<T>,
                               interface_RefEntity<T> {

  interface_FieldEntity(const boost::shared_ptr<T> &sptr)
      : interface_Field<T>(sptr), interface_RefEntity<T>(sptr) {}

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

  /// @return get pointer to field data structure
  inline boost::shared_ptr<Field> &getFieldPtr() const {
    return this->sFieldPtr->getFieldPtr();
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
            tag<FieldName_mi_tag>,
            const_mem_fun<FieldEntity::interface_type_Field, boost::string_ref,
                          &FieldEntity::getNameRef>>,
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
typedef FieldEntity_multiIndex::index<FieldName_mi_tag>::type
    FieldEntityByFieldName;

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

/**
 * \brief Keeps basic information about entity on the finite element
 */
struct BaseFEEntity {
  BaseFEEntity(const boost::shared_ptr<SideNumber> &side_number_ptr)
      : sideNumberPtr(side_number_ptr){};
  virtual ~BaseFEEntity() = default;
  boost::shared_ptr<SideNumber> sideNumberPtr;
  inline int getSideNumber() { return sideNumberPtr->side_number; }
};

} // namespace MoFEM

#endif // __ENTSMULTIINDICES_HPP__

/**
 * \defgroup ent_multi_indices Entities structures and multi-indices
 * \ingroup mofem
 **/
