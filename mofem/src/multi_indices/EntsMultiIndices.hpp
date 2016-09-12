/** \file EntsMultiIndices.hpp
 * \brief Multi-index contains, for mofem entities data structures and other low-level functions
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
 * \brief keeps information about side number for the finite element
 * \ingroup ent_multi_indices
 */
struct __attribute__((__packed__))  SideNumber {
  EntityHandle ent;
  char side_number;
  char sense;
  char offset;
  char brother_side_number;
  inline EntityType getEntType() const {
    return (EntityType)((ent&MB_TYPE_MASK)>>MB_ID_WIDTH);
  }

  // /** \deprecated Use getEntType() instead
  // */
  // DEPRECATED EntityType get_ent_type() const { return getEntType(); }

  SideNumber(EntityHandle _ent,int _side_number,int _sense,int _offset):
  ent(_ent),
  side_number(_side_number),
  sense(_sense),
  offset(_offset),
  brother_side_number(-1) {
  }
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
    hashed_unique<
      member<SideNumber,EntityHandle,&SideNumber::ent>
    >,
    ordered_non_unique<
      composite_key<
      SideNumber,
      const_mem_fun<SideNumber,EntityType,&SideNumber::getEntType>,
      member<SideNumber,char,&SideNumber::side_number> >
    >,
    ordered_non_unique<
      const_mem_fun<SideNumber,EntityType,&SideNumber::getEntType>
    >
  > > SideNumber_multiIndex;

/**
 * \brief Basic data. like access to moab interface and basic tag handlers.
 */
struct BasicEntityData {
  moab::Interface &moab;
  Tag th_RefParentHandle;
  Tag th_RefBitLevel;
  BasicEntityData(moab::Interface &mfield);
  virtual ~BasicEntityData();
};

/**
 * \brief this struct keeps basic methods for moab entity
 * \ingroup ent_multi_indices

  \todo BasicEntity in should be linked to directly to MoAB data structures
  such that connectivity and nodal coordinates could be quickly accessed,
  without need of using native MoAB functions.

 */
struct BasicEntity {

  boost::shared_ptr<BasicEntityData> basicDataPtr;

  EntityHandle ent;
  int owner_proc;
  EntityHandle moab_owner_handle;

  BasicEntity(
    boost::shared_ptr<BasicEntityData> basic_data_ptr,const EntityHandle ent
  );

  /** \brief Get entity type
  */
  inline EntityType getEntType() const { return (EntityType)((ent&MB_TYPE_MASK)>>MB_ID_WIDTH); }

  // /** \deprecated Name changed to getEntType()
  // */
  // DEPRECATED inline EntityType get_ent_type() const { return getEntType(); }

  /** \brief get entity id
  */
  inline EntityID getEntId() const { return (EntityID)(ent&MB_ID_MASK); };

  // /** \deprecated Name changed to getEntId()
  // */
  // inline DEPRECATED EntityID get_ent_id() const { return getEntId(); };

  /** \brief Owner handle on this or other processors
    */
  inline EntityHandle getOwnerEnt() const { return moab_owner_handle; }

  // /** \deprecated Name changed to getOwnerEnt()
  // */
  // DEPRECATED inline EntityHandle get_owner_ent() const { return getOwnerEnt(); }


  /** \brief Get processor owning entity
    */
  inline EntityHandle getOwnerProc() const { return owner_proc; }

  // /** \deprecated Name changed to getOwnerProc()
  // */
  // DEPRECATED inline EntityHandle get_owner_proc() const { return getOwnerProc(); }


  /** \brief get pstatus
    * This tag stores various aspects of parallel status in bits; see also
    * #define's following, to be used in bit mask operations.  If an entity is
    * not shared with any other processors, the pstatus is 0, otherwise it's > 0
    *
    * bit 0: !owned (0=owned, 1=not owned)
    * bit 1: shared (0=not shared, 1=shared)
    * bit 2: multishared (shared by > 2 procs; 0=not shared, 1=shared)
    * bit 3: interface (0=not interface, 1=interface)
    * bit 4: ghost (0=not ghost, 1=ghost)
    *
    */
  unsigned char getPStatus() const;

  // /** \deprecated Name changed to getPStatus()
  // */
  // DEPRECATED inline unsigned char get_pstatus() const { return getPStatus(); }

  /** \berief get shared processors

  Returning list to shared processors. Lists end with -1. Returns NULL if not
  sharing processors.

  DO NOT MODIFY LIST.

  \code
    BasicEntity *ent_ptr = BasicEntity(moan,entity_handle);
    for(int proc = 0; proc<MAX_SHARING_PROCS && -1 != ent_ptr->getSharingProcsPtr[proc]; proc++) {
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
  int* getSharingProcsPtr() const {
    MoABErrorCode rval;
    moab::Interface &moab = basicDataPtr->moab;
    int *sharing_procs_ptr = NULL;
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(getPStatus() & PSTATUS_MULTISHARED) {
      // entity is multi shared
      rval = moab.tag_get_by_ptr(pcomm->sharedps_tag(),&ent,1,(const void **)&sharing_procs_ptr); CHKERR_MOAB(rval);
    } else if(getPStatus() & PSTATUS_SHARED) {
      // shared
      rval = moab.tag_get_by_ptr(pcomm->sharedp_tag(),&ent,1,(const void **)&sharing_procs_ptr); CHKERR_MOAB(rval);
    }
    return sharing_procs_ptr;
  }

  // /** \deprecated Name changed to getSharingProcsPtr()
  // */
  // DEPRECATED int* get_sharing_procs_ptr() const { return getSharingProcsPtr(); }


  /** \berief get sharid entity handlers

  Returning list to shared entity hanlders. Use it with getSharingProcsPtr()

  DO NOT MODIFY LIST.

\code
  BasicEntity *ent_ptr = BasicEntity(moan,entity_handle);
  for(int proc = 0; proc<MAX_SHARING_PROCS && -1 != ent_ptr->getSharingProcsPtr[proc]; proc++) {
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
  inline EntityHandle* getSharingHandlersPtr() const {
    MoABErrorCode rval;
    EntityHandle *sharing_handlers_ptr = NULL;
    moab::Interface &moab = basicDataPtr->moab;
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(getPStatus() & PSTATUS_MULTISHARED) {
      // entity is multi shared
      rval = moab.tag_get_by_ptr(pcomm->sharedhs_tag(),&ent,1,(const void **)&sharing_handlers_ptr); CHKERR_MOAB(rval);
    } else if(getPStatus() & PSTATUS_SHARED) {
      // shared
      rval = moab.tag_get_by_ptr(pcomm->sharedh_tag(),&ent,1,(const void **)&sharing_handlers_ptr); CHKERR_MOAB(rval);
    }
    return sharing_handlers_ptr;
  }

  // /** \deprecated Name changed to getSharingHandlersPtr()
  // */
  // DEPRECATED inline EntityHandle* get_sharing_handlers_ptr() const { return getSharingHandlersPtr(); }

};

/**
 * \brief Struct keeps handle to refined handle.
 * \ingroup ent_multi_indices

  \todo th_RefType "_RefType" is set as two integers, need to be fixed, it is
  waste of space.

 */
struct RefEntity: public BasicEntity {

  RefEntity(boost::shared_ptr<BasicEntityData> basic_data_ptr,const EntityHandle _ent);

  static PetscErrorCode getPatentEnt(Interface &moab,Range ents,std::vector<EntityHandle> vec_patent_ent);
  static PetscErrorCode getBitRefLevel(Interface &moab,Range ents,std::vector<BitRefLevel> vec_bit_ref_level);

  /**
   * \brief Get pointer to parent entity tag.
   *
   * Each refined entity has his parent. Such information is stored on tags.
   * This function get pinter to tag.
   *
   * @return Pointer to tag on entity
   */
  EntityHandle* getParentEntPtr() const;

  /**
   * \brief Get pointer to bit ref level tag

   * Every entity belongs to some refinement level or levels. Each level is marked
   * by bit set in BitRefLevel() (bitset) structure.
   * FIXME: Better explanation here needed.

   * @return Return pointer to tag.
   */
  BitRefLevel* getBitRefLevelPtr() const;

  /** \brief Get entity
  */
  inline EntityHandle getRefEnt() const { return ent; }

  // /** \deprecated Name changed to getRefEnt()
  // */
  // DEPRECATED inline EntityHandle get_ref_ent() const { return getRefEnt(); }


  /** \brief Get patent entity
  */
  inline EntityType getParentEntType() const {
    EntityHandle* tag_parent_ent = getParentEntPtr();
    if(*tag_parent_ent == 0)  return MBMAXTYPE;
    return (EntityType)((*tag_parent_ent&MB_TYPE_MASK)>>MB_ID_WIDTH);
  }

  // /** \deprecated Name changed to getParentEntType()
  // */
  // DEPRECATED inline EntityType get_parent_ent_type() const { return getParentEntType(); }

  /** \brief Get parent entity, i.e. entity form one refinement level up
  */
  inline EntityHandle getParentEnt() const {
    EntityHandle* tag_parent_ent = getParentEntPtr();
    return *tag_parent_ent;
  }

  // /** \deprecated Name changed to getParentEnt()
  // */
  // DEPRECATED inline EntityHandle get_parent_ent() const { return getParentEnt(); }

  /** \brief Get entity ref bit refinement signature
  */
  inline const BitRefLevel& getBitRefLevel() const { return *getBitRefLevelPtr(); }

  // /** \deprecated Name changed to getBitRefLevel()
  // */
  // DEPRECATED inline const BitRefLevel& get_BitRefLevel() const { return getBitRefLevel(); }


  friend std::ostream& operator<<(std::ostream& os,const RefEntity& e);

};


/**
 * \brief interface to RefEntity
 * \ingroup ent_multi_indices
 */
template <typename T>
struct interface_RefEntity {

  const boost::shared_ptr<T> sPtr;
  interface_RefEntity(const boost::shared_ptr<T> sptr):
  sPtr(sptr) {}

  inline EntityHandle getRefEnt() const { return this->sPtr->getRefEnt(); }

  // /** \deprecated Name changed to getRefEnt()
  // */
  // DEPRECATED inline EntityHandle get_ref_ent() const { return this->sPtr->getRefEnt(); }

  inline EntityType getParentEntType() const { return this->sPtr->getParentEntType(); };

  // /** \deprecated Name changed to getParentEntType()
  // */
  // DEPRECATED inline EntityType get_parent_ent_type() const { return this->sPtr->getParentEntType(); };

  inline EntityHandle getParentEnt() const { return this->sPtr->getParentEnt(); }

  // /** \deprecated Name changed to getParentEnt()
  // */
  // DEPRECATED inline EntityHandle get_parent_ent() const { return this->sPtr->getParentEnt(); }
  inline const BitRefLevel& getBitRefLevel() const { return this->sPtr->getBitRefLevel(); }

  // /** \deprecated Name changed to getBitRefLevel()
  // */
  // DEPRECATED inline const BitRefLevel& get_BitRefLevel() const { return this->sPtr->getBitRefLevel(); }

  inline EntityType getEntType() const { return this->sPtr->getEntType(); };

  // /** \deprecated Name changed to getEntType()
  // */
  // DEPRECATED inline EntityType get_ent_type() const { return this->sPtr->getEntType(); };

  inline EntityID getEntId() const { return this->sPtr->getEntId(); };

  // /** \deprecated Name changed to getEntId()
  // */
  // DEPRECATED inline EntityID get_ent_id() const { return this->sPtr->getEntId(); };

  inline EntityHandle getOwnerEnt() const { return this->sPtr->getOwnerEnt(); }

  // /** \deprecated Name changed to getOwnerEnt()
  // */
  // DEPRECATED inline EntityHandle get_owner_ent() const { return this->sPtr->getOwnerEnt(); }

  inline EntityHandle getOwnerProc() const { return this->sPtr->getOwnerProc(); }

  // /** \deprecated Name changed to getOwnerProc()
  // */
  // DEPRECATED EntityHandle get_owner_proc() const { return this->sPtr->getOwnerProc(); }

  inline unsigned char getPStatus() const { return this->sPtr->getPStatus(); }

  // /** \deprecated Name changed to getPStatus()
  // */
  // DEPRECATED inline unsigned char get_pstatus() const { return this->sPtr->getPStatus(); }

  inline int* getSharingProcsPtr() const { return this->sPtr->getSharingProcsPtr(); }

  // /** \deprecated Name changed to getSharingProcsPtr()
  // */
  // DEPRECATED inline int* get_sharing_procs_ptr() const { return this->sPtr->getSharingProcsPtr(); }

  inline EntityHandle* getSharingHandlersPtr() const { return this->sPtr->getSharingHandlersPtr(); }

  // /** \deprecated Name changed to getSharingHandlersPtr()
  // */
  // DEPRECATED inline EntityHandle* get_sharing_handlers_ptr() const { return this->sPtr->getSharingHandlersPtr(); }

  virtual ~interface_RefEntity() {}

  inline const boost::shared_ptr<T> getRefEntityPtr() {
    return this->sPtr;
  }

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
 * \param ordered_non_unique Composite_Ent_And_ParentEntType_mi_tag
 */
typedef multi_index_container<
  boost::shared_ptr<RefEntity>,
  indexed_by<
    hashed_unique<
      tag<Ent_mi_tag>, member<RefEntity::BasicEntity,EntityHandle,&RefEntity::ent> >,
    hashed_non_unique<
      tag<Ent_Owner_mi_tag>, member<RefEntity::BasicEntity,EntityHandle,&RefEntity::moab_owner_handle> >,
    ordered_non_unique<
      tag<Proc_mi_tag>, member<RefEntity::BasicEntity,int,&RefEntity::owner_proc> >,
    ordered_non_unique<
      tag<Ent_ParallelStatus>, const_mem_fun<RefEntity::BasicEntity,unsigned char,&RefEntity::getPStatus> >,
    ordered_non_unique<
      tag<Ent_Ent_mi_tag>, const_mem_fun<RefEntity,EntityHandle,&RefEntity::getParentEnt> >,
    ordered_non_unique<
      tag<EntType_mi_tag>, const_mem_fun<RefEntity::BasicEntity,EntityType,&RefEntity::getEntType> >,
    ordered_non_unique<
      tag<ParentEntType_mi_tag>, const_mem_fun<RefEntity,EntityType,&RefEntity::getParentEntType> >,
    ordered_non_unique<
      tag<Composite_EntType_and_ParentEntType_mi_tag>,
      composite_key<
      	RefEntity,
      	const_mem_fun<RefEntity::BasicEntity,EntityType,&RefEntity::getEntType>,
      	const_mem_fun<RefEntity,EntityType,&RefEntity::getParentEntType> > >,
    ordered_non_unique<
      tag<Composite_Ent_And_ParentEntType_mi_tag>,
      composite_key<
      	RefEntity,
      	const_mem_fun<RefEntity,EntityHandle,&RefEntity::getParentEnt>,
      	const_mem_fun<RefEntity::BasicEntity,EntityType,&RefEntity::getEntType> > >
  > > RefEntity_multiIndex;

/** \brief multi-index view of RefEntity by parent entity
  \ingroup ent_multi_indices
*/
typedef multi_index_container<
  boost::shared_ptr<RefEntity>,
  indexed_by<
    hashed_unique<
      const_mem_fun<RefEntity,EntityHandle,&RefEntity::getParentEnt> >,
    hashed_unique<
      tag<Composite_EntType_and_ParentEntType_mi_tag>,
    composite_key<
	    boost::shared_ptr<RefEntity>,
	    const_mem_fun<RefEntity,EntityHandle,&RefEntity::getRefEnt>,
	    const_mem_fun<RefEntity,EntityHandle,&RefEntity::getParentEnt> > >
  > > RefEntity_multiIndex_view_by_parent_entity;

/** \brief ref mofem entity, remove parent
 * \ingroup ent_multi_indices
 */
struct RefEntity_change_remove_parent {
  ErrorCode rval;
  RefEntity_change_remove_parent() {
  }
  void operator()(boost::shared_ptr<RefEntity> &e) {
    rval = e->basicDataPtr->moab.tag_delete_data(
      e->basicDataPtr->th_RefParentHandle,&e->ent,1
    ); MOAB_THROW(rval);
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
  ErrorCode rval;
  RefEntity_change_parent(EntityHandle parent): pArent(parent) {}
  void operator()(boost::shared_ptr<RefEntity> &e) {
    *(e->getParentEntPtr()) = pArent;
  }
};

/** \brief ref mofem entity, left shift
  * \ingroup ent_multi_indices
  */
struct RefEntity_change_left_shift {
  int shift;
  RefEntity_change_left_shift(const int _shift): shift(_shift) {};
  void operator()(boost::shared_ptr<RefEntity> &e) { (*e->getBitRefLevelPtr())<<=shift;  };
};

/** \brief ref mofem entity, right shift
 * \ingroup ent_multi_indices
  */
struct RefEntity_change_right_shift {
  int shift;
  RefEntity_change_right_shift(const int _shift): shift(_shift) {};
  void operator()(boost::shared_ptr<RefEntity> &e) { *(e->getBitRefLevelPtr())>>=shift;  };
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices
  */
struct RefEntity_change_add_bit {
  BitRefLevel bit;
  RefEntity_change_add_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(boost::shared_ptr<RefEntity> &e) {
    bit |= *(e->getBitRefLevelPtr());
    *(e->getBitRefLevelPtr()) = bit;
  }
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices
  */
struct RefEntity_change_and_bit {
  BitRefLevel bit;
  RefEntity_change_and_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(boost::shared_ptr<RefEntity> &e) {
    bit &= *(e->getBitRefLevelPtr());
    *(e->getBitRefLevelPtr()) = bit;
  }
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices
  */
struct RefEntity_change_xor_bit {
  BitRefLevel bit;
  RefEntity_change_xor_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(boost::shared_ptr<RefEntity> &e) {
    bit ^= *(e->getBitRefLevelPtr());
    *(e->getBitRefLevelPtr()) = bit;
  }
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices
  */
struct RefEntity_change_set_bit {
  BitRefLevel bit;
  RefEntity_change_set_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(boost::shared_ptr<RefEntity> &e) {
    *(e->getBitRefLevelPtr()) = bit;
  }
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices
  */
struct RefEntity_change_set_nth_bit {
  int n;
  bool b;
  RefEntity_change_set_nth_bit(const int _n,bool _b): n(_n),b(_b) {};
  void operator()(boost::shared_ptr<RefEntity> &e) {
    (*(e->getBitRefLevelPtr()))[n] = b;
  }
};

/**
  * \brief Struct keeps handle to entity in the field.
  * \ingroup ent_multi_indices
  */
struct MoFEMEntity:
  public
  interface_Field<Field>,
  interface_RefEntity<RefEntity> {

  typedef interface_Field<Field> interface_type_Field;
  typedef interface_RefEntity<RefEntity> interface_type_RefEntity;
  const FieldData* tag_FieldData;
  int tag_FieldData_size;
  const ApproximationOrder* tag_dof_order_data;
  const FieldCoefficientsNumber* tag_dof_rank_data;
  MoFEMEntity(
    const boost::shared_ptr<Field> field_ptr,
    const boost::shared_ptr<RefEntity> ref_ent_ptr
  );
  ~MoFEMEntity();

  /**
   * \brief Get entity handle
   * @return EntityHandle
   */
  inline EntityHandle getEnt() const { return getRefEnt(); }

  // /** \deprecated Name change to getEnt()
  // */
  // DEPRECATED EntityHandle get_ent() const { return getEnt(); }

  /**
   * \brief Get number of DOFs on entity
   * @return Number of DOFs
   */
  inline int getNbDofsOnEnt() const { return tag_FieldData_size/sizeof(FieldData); }

  // /** \deprecated Name change to getNbDofsOnEnt()
  // */
  // DEPRECATED inline int get_nb_dofs_on_ent() const { return getNbDofsOnEnt(); }

  /**
   * \brief Get Vector of DOFs values on entity
   * @return Vector of DOFs values
   */
  inline VectorAdaptor getEntFieldData() const {
    int size = getNbDofsOnEnt();
    double* ptr = const_cast<FieldData*>(tag_FieldData);
    return VectorAdaptor(size,ublas::shallow_array_adaptor<FieldData>(size,ptr));
  }

  // /** \deprecated Name changed to getEntFieldData()
  // */
  // DEPRECATED inline VectorAdaptor get_ent_FieldData() const { return getEntFieldData(); }

  /**
   * \brief Get number of DOFs on entity for given order of approximation
   * @param  order Order of approximation
   * @return       Number of DOFs
   */
  inline int getOrderNbDofs(int order) const { return (this->sFieldPtr->forder_table[getEntType()])(order); }

  // /** \deprecated Name changed to getOrderNbDofs()
  // */
  // DEPRECATED inline int get_order_nb_dofs(int order) const { return getOrderNbDofs(order); }

  /**
   * \brief Get difference of number of DOFs between order and order-1
   * @param  order Approximation order
   * @return       Difference number of DOFs
   */
  inline int getOrderNbDofsDiff(int order) const { return getOrderNbDofs(order)-getOrderNbDofs(order-1); }

  // /** \deprecated Name changed to getOrderNbDofsDiff()
  // */
  // DEPRECATED inline int get_order_nb_dofs_diff(int order) const { return getOrderNbDofsDiff(order); }

  /**
   * \brief Get pinter to Tag keeping approximation order
   * @return Pointer to Tag
   */
  ApproximationOrder* getMaxOrderPtr();

  /**
   * \brief Get order set to the entity (Allocated tag size for such number)
   * @return Approximation order
   */
  ApproximationOrder getMaxOrder() const;

  // /** \deprecated Name changed to getMaxOrder()
  // */
  // DEPRECATED inline ApproximationOrder get_max_order() const { return getMaxOrder(); }

  GlobalUId global_uid; ///< Global unique id for this entity

  /**
   * \brief Get global unique id
   * @return Global UId
   */
  const GlobalUId& getGlobalUniqueId() const { return global_uid; }

  // /** \deprecated Name changed getGlobalUniqueId()
  // */
  // DEPRECATED const GlobalUId& get_global_unique_id() const { return getGlobalUniqueId(); }

  /**
   * \brief Calculate global UId
   * @return Global UId
   */
  inline GlobalUId getGlobalUniqueIdCalculate() const {
    const char bit_number = getBitNumber();
    assert(bit_number<32);
    assert(sPtr->owner_proc<1024);
    GlobalUId _uid_ = (UId)0;
    _uid_ |= (UId)sPtr->moab_owner_handle;
    _uid_ |= (UId)bit_number << 8*sizeof(EntityHandle);
    _uid_ |= (UId)sPtr->owner_proc << 5+8*sizeof(EntityHandle);
    return _uid_;
  }

  /**
   * \brief Get pointer to RefEntity
   */
  inline const boost::shared_ptr<RefEntity> getRefEntityPtr() { return this->sPtr; }

  // /** \deprecated Name changed to getRefEntityPtr()
  // */
  // DEPRECATED inline const boost::shared_ptr<RefEntity> get_RefEntity_ptr() { return getRefEntityPtr(); }

  /**
   * \brief Get pointer to Field data structure associated with this entity
   */
  inline const boost::shared_ptr<Field> getFieldPtr() const { return this->sFieldPtr; }

  // /** \deprecated Name changed to getFieldPtr()
  // */
  // DEPRECATED inline const boost::shared_ptr<Field> get_Field_ptr() const { return getFieldPtr(); }

  friend std::ostream& operator<<(std::ostream& os,const MoFEMEntity& e);

};

/**
 * \brief Interface to MoFEMEntity
 * \ingroup ent_multi_indices
 *
 * interface to MoFEMEntity
 */
template <typename T>
struct interface_MoFEMEntity:
public
interface_Field<T>,
interface_RefEntity<T> {

  interface_MoFEMEntity(const boost::shared_ptr<T> sptr):
  interface_Field<T>(sptr),
  interface_RefEntity<T>(sptr) {
  };
  inline EntityHandle getEnt() const { return this->sPtr->getEnt(); }

  // /** \deprecated Use getEnt() instead
  // */
  // DEPRECATED inline EntityHandle get_ent() const { return this->sPtr->getEnt(); }

  inline int getNbDofsOnEnt() const { return this->sPtr->getNbDofsOnEnt(); }

  // /** \deprecated Chnage name to getNbDofsOnEnt()
  // */
  // DEPRECATED inline int get_nb_dofs_on_ent() const { return this->sPtr->getNbDofsOnEnt(); }

  inline VectorAdaptor getEntFieldData() const { return this->sPtr->getEntFieldData(); }

  // /** \deprecated Change name getEntFieldData()
  // */
  // DEPRECATED inline VectorAdaptor get_ent_FieldData() const { return this->sPtr->getEntFieldData(); }

  inline int getOrderNbDofs(int order) const { return this->sFieldPtr->getOrderNbDofs(order); }

  // /** \deprecated Change name to getOrderNbDofs()
  // */
  // DEPRECATED inline int get_order_nb_dofs(int order) const { return this->sFieldPtr->getOrderNbDofs(order); }

  inline int getOrderNbDofsDiff(int order) const { return this->sPtr->getOrderNbDofsDiff(order); }

  // /** \deprecated Chnage name to getOrderNbDofsDiff()
  // */
  // DEPRECATED inline int get_order_nb_dofs_diff(int order) const { return this->sPtr->getOrderNbDofsDiff(order); }

  inline ApproximationOrder getMaxOrder() const { return this->sPtr->getMaxOrder(); }

  // /** \deprecated Change name to getMaxOrder()
  // */
  // DEPRECATED inline ApproximationOrder get_max_order() const { return this->sPtr->getMaxOrder(); }

  inline GlobalUId getGlobalUniqueId() const { return this->sPtr->getGlobalUniqueId(); }

  // /** \deprecated Change name to getGlobalUniqueId()
  // */
  // DEPRECATED inline GlobalUId get_global_unique_id() const { return this->sPtr->getGlobalUniqueId(); }

  inline const boost::shared_ptr<RefEntity> getRefEntityPtr() { return this->sPtr->getRefEntityPtr(); }

  // /** \deprecated Change name to getRefEntityPtr()
  // */
  // DEPRECATED inline const boost::shared_ptr<RefEntity> get_RefEntity_ptr() { return this->sPtr->getRefEntityPtr(); }

  inline const boost::shared_ptr<Field> getFieldPtr() const { return this->sFieldPtr->getFieldPtr(); }

  // /** \deprecated Change name to getFieldPtr()
  // */
  // DEPRECATED inline const boost::shared_ptr<Field> get_Field_ptr() const { return this->sFieldPtr->getFieldPtr(); }

  inline const boost::shared_ptr<MoFEMEntity> getMoFEMEntityPtr() const { return this->sPtr; };

  // /** \deprecated Name is deprecated, use getMoFEMEntityPtr instead.
  // */
  // DEPRECATED inline const boost::shared_ptr<MoFEMEntity> get_MoFEMEntity_ptr() const { return getMoFEMEntityPtr(); };

};

/**
 * \brief structure to chane MoFEMEntity order
 * \ingroup ent_multi_indices
 */
struct MoFEMEntity_change_order {
  ApproximationOrder order;
  std::vector<FieldData> data;
  std::vector<ApproximationOrder> data_dof_order;
  std::vector<FieldCoefficientsNumber> data_dof_rank;
  MoFEMEntity_change_order(ApproximationOrder _order):
  order(_order) {};
  void operator()(boost::shared_ptr<MoFEMEntity> &e);
};

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps MoFEMEntity
 * \ingroup ent_multi_indices
 *
 */
typedef multi_index_container<
  boost::shared_ptr<MoFEMEntity>,
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>, member<MoFEMEntity,GlobalUId,&MoFEMEntity::global_uid> >,
    ordered_non_unique<
      tag<Ent_ParallelStatus>, const_mem_fun<MoFEMEntity::interface_type_RefEntity,unsigned char,&MoFEMEntity::getPStatus> >,
    ordered_non_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_Field,const BitFieldId&,&MoFEMEntity::getId>, LtBit<BitFieldId> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_Field,boost::string_ref,&MoFEMEntity::getNameRef> >,
    hashed_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::getEnt> >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_mi_tag>,
      composite_key<
      	MoFEMEntity,
      	const_mem_fun<MoFEMEntity::interface_type_Field,boost::string_ref,&MoFEMEntity::getNameRef>,
      	const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::getEnt>
      > >
  > > MoFEMEntity_multiIndex;

  typedef multi_index_container<
    boost::shared_ptr<MoFEMEntity>,
    indexed_by<
      hashed_non_unique<
        tag<Ent_mi_tag>, const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::getEnt> >
  > > MoFEMEntity_multiIndex_ent_view;

}

#endif // __ENTSMULTIINDICES_HPP__

/***************************************************************************//**
 * \defgroup ent_multi_indices Entities structures and multi-indices
 * \ingroup mofem
 ******************************************************************************/
