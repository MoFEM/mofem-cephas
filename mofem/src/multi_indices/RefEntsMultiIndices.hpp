/** \file RefEntsMultiIndices.hpp
 * \brief Multi-index contains, for mofem entities data structures and other
 * low-level functions
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __REF_ENTSMULTIINDICES_HPP__
#define __REF_ENTSMULTIINDICES_HPP__

namespace MoFEM {

template <EntityType TYPE> inline EntityHandle get_id_for_max_type() {
  return (static_cast<EntityHandle>(TYPE) << MB_ID_WIDTH) |
         (~static_cast<EntityHandle>(MB_TYPE_MASK));
};

template <EntityType TYPE> inline EntityHandle get_id_for_min_type() {
  return (static_cast<EntityHandle>(TYPE) << MB_ID_WIDTH);
};

inline EntityHandle get_id_for_max_type(const EntityType type) {
  return (static_cast<EntityHandle>(type) << MB_ID_WIDTH) |
         (~static_cast<EntityHandle>(MB_TYPE_MASK));
};

inline EntityHandle get_id_for_min_type(const EntityType type) {
  return (static_cast<EntityHandle>(type) << MB_ID_WIDTH);
};

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
    if (tag_size)
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
struct /*__attribute__((__packed__))*/ SideNumber {
  EntityHandle ent;
  signed char side_number;
  signed char sense;
  signed char offset;
  signed char brother_side_number;

  inline EntityType getEntType() const {
    return static_cast<EntityType>((ent & MB_TYPE_MASK) >> MB_ID_WIDTH);
  }

  inline int getEntDimension() const {
    return moab::CN::Dimension(getEntType());
  }

  inline const char *getEntTypeName() const {
    return moab::CN::EntityTypeName(getEntType());
  }

  SideNumber(EntityHandle ent, int side_number, int sense, int offset)
      : ent(ent), side_number(side_number), sense(sense), offset(offset),
        brother_side_number(-1) {}
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
        ordered_unique<member<SideNumber, EntityHandle, &SideNumber::ent>>,
        ordered_non_unique<

            composite_key<
                SideNumber,
                const_mem_fun<SideNumber, EntityType, &SideNumber::getEntType>,
                member<SideNumber, signed char, &SideNumber::side_number>

                >>

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
  Tag th_MeshsetPart;
  BasicEntityData(const moab::Interface &mfield,
                  const int pcomm_id = MYPCOMM_INDEX);
  virtual ~BasicEntityData() = default;
};

template <int N> struct RefEntityTmp : public RefEntityTmp<N - 1> {

  using RefEntityTmp<N - 1>::RefEntityTmp;

  virtual const boost::shared_ptr<BasicEntityData> getBasicDataPtr() const {
    if (auto ptr = basicDataPtr.lock())
      return ptr;
    else
      return nullptr;
  }

  static boost::weak_ptr<BasicEntityData> basicDataPtr;
};

template <int N> boost::weak_ptr<BasicEntityData> RefEntityTmp<N>::basicDataPtr;

struct RefElement;

/**
* \brief Struct keeps handle to refined handle.
* \ingroup ent_multi_indices

\todo th_RefType "_RefType" is set as two integers, need to be fixed, it is
waste of space.

*/
template <> struct RefEntityTmp<0> {

  RefEntityTmp(const boost::shared_ptr<BasicEntityData> &basic_data_ptr,
               const EntityHandle ent);

  RefEntityTmp(const boost::shared_ptr<BasicEntityData> &basic_data_ptr,
               const EntityHandle ent, EntityHandle *ent_parent_tag_ptr);

  virtual ~RefEntityTmp() = default;

  /**
   * @brief Get pointer to basic data struture
   *
   * @return const boost::shared_ptr<BasicEntityData>
   */
  virtual const boost::shared_ptr<BasicEntityData> getBasicDataPtr() const {
    if (auto ptr = basicDataPtr.lock())
      return ptr;
    else
      return nullptr;
  }

  int getSideNumber() const;

  /**
   * @brief Get the Side number
   *
   * @return int
   */
  boost::shared_ptr<SideNumber> getSideNumberPtr() const;

  /**
   * @brief Get the entity handle
   *
   * @return EntityHandle
   */
  inline EntityHandle getEnt() const { return this->ent; }

  /** \brief Get entity type
   */
  inline EntityType getEntType() const {
    return (EntityType)((this->ent & MB_TYPE_MASK) >> MB_ID_WIDTH);
  }

  /**
   * @brief Get the Ent Dimension 
   * 
   * @return int 
   */
  inline int getEntDimension() const {
    return moab::CN::Dimension(getEntType());
  }

  /**
   * @brief Get the Ent Type Name 
   * 
   * @return const char* 
   */
  inline const char *getEntTypeName() const {
    return moab::CN::EntityTypeName(getEntType());
  }

  /** \brief get entity id
   */
  inline EntityID getEntId() const {
    return (EntityID)(this->ent & MB_ID_MASK);
  };

  static inline EntityHandle
  getOwnerEnt(const EntityHandle ent,
              boost::shared_ptr<BasicEntityData> basic_ent_data) {
    auto *pcomm =
        ParallelComm::get_pcomm(&basic_ent_data->moab, basic_ent_data->pcommID);

    auto get_owner_ent = [&](auto tag) {
      if (auto owner_ent_ptr = static_cast<EntityHandle *>(
              MoFEM::get_tag_ptr(basic_ent_data->moab, tag, ent, NULL))) {
        return owner_ent_ptr[0];
      } else {
        return ent;
      }
    };

    if (auto pstat_ptr = static_cast<unsigned char *>(MoFEM::get_tag_ptr(
            basic_ent_data->moab, pcomm->pstatus_tag(), ent, NULL))) {
      if (!(*pstat_ptr & PSTATUS_NOT_OWNED)) {
        return ent;
      } else if (*pstat_ptr & PSTATUS_MULTISHARED) {
        return get_owner_ent(pcomm->sharedhs_tag());
      } else if (*pstat_ptr & PSTATUS_SHARED) {
        return get_owner_ent(pcomm->sharedh_tag());
      } else {
        return ent;
      }
    } else {
      return ent;
    }
  }

  /** \brief Owner handle on this or other processors
   */
  inline EntityHandle getOwnerEnt() const {
    return getOwnerEnt(this->ent, this->getBasicDataPtr());
  }

  static inline int
  getOwnerProc(const EntityHandle ent,
               boost::shared_ptr<BasicEntityData> basic_ent_data) {
    auto *pcomm =
        ParallelComm::get_pcomm(&basic_ent_data->moab, basic_ent_data->pcommID);

    auto get_owner_proc = [&](auto tag) {
      if (auto owner_proc = static_cast<int *>(
              MoFEM::get_tag_ptr(basic_ent_data->moab, tag, ent, NULL))) {
        return owner_proc[0];
      } else {
        return 0;
      }
    };

    if (auto pstat_ptr = static_cast<unsigned char *>(MoFEM::get_tag_ptr(
            basic_ent_data->moab, pcomm->pstatus_tag(), ent, NULL))) {
      if (!(*pstat_ptr & PSTATUS_NOT_OWNED)) {
        return pcomm->rank();
      } else if (*pstat_ptr & PSTATUS_MULTISHARED) {
        return get_owner_proc(pcomm->sharedps_tag());
      } else if (*pstat_ptr & PSTATUS_SHARED) {
        return get_owner_proc(pcomm->sharedp_tag());
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  }

  inline int getOwnerProc() const {
    return getOwnerProc(this->ent, this->getBasicDataPtr());
  }

  static inline int
  getPartProc(const EntityHandle ent, const Tag th,
              const boost::shared_ptr<BasicEntityData> &basic_ptr) {
    if (auto part_ptr = static_cast<int *>(
            MoFEM::get_tag_ptr(basic_ptr->moab, th, ent, NULL))) {
      return *part_ptr >= 0 ? *part_ptr : 0;
    } else
      return 0;
  };

  /** \brief Get processor
   */
  inline int getPartProc() const {
    if (PetscUnlikely(getEntType() == MBENTITYSET)) {
      return getPartProc(ent, this->getBasicDataPtr()->th_MeshsetPart,
                         this->getBasicDataPtr());
    } else {
      auto *pcomm = ParallelComm::get_pcomm(&this->getBasicDataPtr()->moab,
                                            this->getBasicDataPtr()->pcommID);
      return getPartProc(ent, pcomm->partition_tag(), this->getBasicDataPtr());
    }
  }

  /**
   * @brief Get pointer to partition tag
   * 
   * @return int* 
   */
  inline int *getPartProcPtr() const {
    if (PetscUnlikely(getEntType() == MBENTITYSET)) {
      return static_cast<int *>(MoFEM::get_tag_ptr(
          this->getBasicDataPtr()->moab,
          this->getBasicDataPtr()->th_MeshsetPart, ent, NULL));
    } else {
      auto *pcomm = ParallelComm::get_pcomm(&this->getBasicDataPtr()->moab,
                                            this->getBasicDataPtr()->pcommID);
      return static_cast<int *>(MoFEM::get_tag_ptr(
          this->getBasicDataPtr()->moab, pcomm->partition_tag(), ent, NULL));
    }
  }

  /** \brief get pstatus
   * This tag stores various aspects of parallel status in bits; see also
   * define following, to be used in bit mask operations.  If an entity is
   * not shared with any other processors, the pstatus is 0, otherwise it's >
   * 0
   *
   * bit 0: !owned (0=owned, 1=not owned)
   * bit 1: shared (0=not shared, 1=shared)
   * bit 2: multishared (shared by > 2 procs; 0=not shared, 1=shared)
   * bit 3: interface (0=not interface, 1=interface)
   * bit 4: ghost (0=not ghost, 1=ghost)
   *
   */
  inline unsigned char getPStatus() const {
    auto *pcomm = ParallelComm::get_pcomm(&this->getBasicDataPtr()->moab,
                                          this->getBasicDataPtr()->pcommID);
    if (auto pstat_ptr = static_cast<unsigned char *>(MoFEM::get_tag_ptr(
            this->getBasicDataPtr()->moab, pcomm->pstatus_tag(), ent, NULL))) {
      return *pstat_ptr;
    } else {
      return 0;
    }
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
  inline EntityHandle *getParentEntPtr() const { return entParentTagPtr; }

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
    return (EntityType)((getParentEnt() & MB_TYPE_MASK) >> MB_ID_WIDTH);
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
  static boost::weak_ptr<BasicEntityData> basicDataPtr;

private:
  friend struct EntFiniteElement;
  friend struct NumeredEntFiniteElement;

  /**
   * @brief Get the pointer to reference element
   *
   * @return const boost::shared_ptr<RefElement>
   */
  inline const boost::shared_ptr<RefElement> getRefElementPtr() const {
    if (auto ptr = refElementPtr.lock())
      return ptr;
    else
      return nullptr;
  }

  static boost::weak_ptr<RefElement> refElementPtr;
  EntityHandle *entParentTagPtr; ///< Tag ptr to parent entity handle
};

template <> struct RefEntityTmp<-1> : public RefEntityTmp<0> {

  RefEntityTmp(const boost::shared_ptr<BasicEntityData> &basic_data_ptr,
               const EntityHandle ent)
      : RefEntityTmp<0>(basic_data_ptr, ent), basicDataPtr(basic_data_ptr) {}

  RefEntityTmp(const boost::shared_ptr<BasicEntityData> &basic_data_ptr,
               const EntityHandle ent, EntityHandle *ent_parent_tag_ptr)
      : RefEntityTmp<0>(basic_data_ptr, ent, ent_parent_tag_ptr),
        basicDataPtr(basic_data_ptr) {}

  virtual const boost::shared_ptr<BasicEntityData> getBasicDataPtr() const {
    if (auto ptr = basicDataPtr.lock())
      return ptr;
    else
      return nullptr;
  }

private:
  mutable boost::weak_ptr<BasicEntityData> basicDataPtr;
};

using RefEntity = RefEntityTmp<0>;

/**
 * \brief interface to RefEntity
 * \ingroup ent_multi_indices
 */
template <typename T> struct interface_RefEntity {

  interface_RefEntity(const boost::shared_ptr<T> &sptr) : sPtr(sptr) {}

  interface_RefEntity(const interface_RefEntity<T> &interface)
      : sPtr(interface.getRefEntityPtr()) {}

  virtual ~interface_RefEntity() = default;

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getSideNumber
   */
  inline int getSideNumber() const { return this->sPtr->getSideNumber(); }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getSideNumberPtr
   */
  inline boost::shared_ptr<SideNumber> getSideNumberPtr() const {
    return this->sPtr->getSideNumberPtr();
  }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getBasicDataPtr
   */
  inline const boost::shared_ptr<BasicEntityData> getBasicDataPtr() const {
    return this->sPtr->getBasicDataPtr();
  }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getEnt
   */
  inline EntityHandle getEnt() const { return this->sPtr->getEnt(); }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getParentEntType
   */
  inline EntityType getParentEntType() const {
    return this->sPtr->getParentEntType();
  };

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getParentEnt
   */
  inline EntityHandle getParentEnt() const {
    return this->sPtr->getParentEnt();
  }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getBitRefLevelPtr
   */
  inline BitRefLevel *getBitRefLevelPtr() const {
    return this->sPtr->getBitRefLevelPtr();
  }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getBitRefLevel
   */
  inline const BitRefLevel &getBitRefLevel() const {
    return this->sPtr->getBitRefLevel();
  }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getBitRefLevelULong
   */
  inline unsigned long int getBitRefLevelULong() const {
    return this->sPtr->getBitRefLevelULong();
  }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getEntType
   */
  inline EntityType getEntType() const { return this->sPtr->getEntType(); };

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getEntDimension
   */
  inline int getEntDimension() const { return this->sPtr->getEntDimension(); }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getEntTypeName
   */
  inline const char *getEntTypeName() const {
    return this->sPtr->getEntTypeName();
  }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getEntId
   */
  inline EntityID getEntId() const { return this->sPtr->getEntId(); };

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getOwnerEnt
   */
  inline EntityHandle getOwnerEnt() const { return this->sPtr->getOwnerEnt(); }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getOwnerEnt
   */
  inline EntityHandle &getOwnerEnt() { return this->sPtr->getOwnerEnt(); }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getOwnerProc
   */
  inline int getOwnerProc() const { return this->sPtr->getOwnerProc(); }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getPartProc
   */
  inline int getPartProc() const { return this->sPtr->getPartProc(); }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getPartProcPtr
   */
  inline int *getPartProcPtr() const { return this->sPtr->getPartProcPtr(); }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getPartProc
   */
  inline unsigned char getPStatus() const { return this->sPtr->getPStatus(); }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getSharingProcsPtr
   */
  inline int *getSharingProcsPtr() const {
    return this->sPtr->getSharingProcsPtr();
  }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getSharingHandlersPtr
   */
  inline EntityHandle *getSharingHandlersPtr() const {
    return this->sPtr->getSharingHandlersPtr();
  }

  /**
   * @copydoc MoFEM::RefEntityTmp<0>::getRefEntityPtr
   */
  inline boost::shared_ptr<T> &getRefEntityPtr() const { return this->sPtr; }

  mutable boost::shared_ptr<T> sPtr;
};

/**
 * \typedef RefEntity_multiIndex
 * type multiIndex container for RefEntity
 * \ingroup ent_multi_indices
 *
 * \param hashed_unique Ent_mi_tag
 * \param ordered_non_unique Ent_Ent_mi_tag
 * \param ordered_non_unique Composite_EntType_And_ParentEntType_mi_tag
 * \param ordered_non_unique Composite_ParentEnt_And_EntType_mi_tag
 */
using RefEntity_multiIndex = multi_index_container<
    boost::shared_ptr<RefEntity>,
    indexed_by<ordered_unique<

                   tag<Ent_mi_tag>,
                   const_mem_fun<RefEntity, EntityHandle, &RefEntity::getEnt>

                   >,

               ordered_non_unique<

                   tag<Ent_Ent_mi_tag>, const_mem_fun<RefEntity, EntityHandle,
                                                      &RefEntity::getParentEnt>

                   >,

               ordered_non_unique<
                   tag<Composite_EntType_and_ParentEntType_mi_tag>,
                   composite_key<RefEntity,
                                 const_mem_fun<RefEntity, EntityType,
                                               &RefEntity::getEntType>,
                                 const_mem_fun<RefEntity, EntityType,
                                               &RefEntity::getParentEntType>>

                   >,

               ordered_non_unique<
                   tag<Composite_ParentEnt_And_EntType_mi_tag>,
                   composite_key<RefEntity,
                                 const_mem_fun<RefEntity, EntityType,
                                               &RefEntity::getEntType>,
                                 const_mem_fun<RefEntity, EntityHandle,
                                               &RefEntity::getParentEnt>>

                   >>

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
                                                  &RefEntity::getEnt>,
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
                                                   &RefEntity::getEnt>,
                                     const_mem_fun<RefEntity, EntityHandle,
                                                   &RefEntity::getParentEnt>>>>

        >;

using RefEntity_multiIndex_view_sequence_ordered_view = multi_index_container<
    boost::shared_ptr<RefEntity>,
    indexed_by<
        sequenced<>,
        ordered_unique<tag<Ent_mi_tag>, const_mem_fun<RefEntity, EntityHandle,
                                                      &RefEntity::getEnt>>>>;

/** \brief change parent
* \ingroup ent_multi_indices
*
* Use this function with care. Some other multi-indices can deponent on
this.

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
} // namespace MoFEM

#endif // __REF_ENTSMULTIINDICES_HPP__

/**
 * \defgroup ent_multi_indices Entities structures and multi-indices
 * \ingroup mofem
 **/