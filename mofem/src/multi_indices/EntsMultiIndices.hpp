/** \file EntsMultiIndices.hpp
 * \brief Myltindex contains, for mofem entities data structures and other low-level functions
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
struct SideNumber {
  EntityHandle ent;
  int side_number;
  int sense;
  int offset;
  int brother_side_number;
  inline EntityType get_ent_type() const { return (EntityType)((ent&MB_TYPE_MASK)>>MB_ID_WIDTH); }
  SideNumber(EntityHandle _ent,int _side_number,int _sense,int _offset):
    ent(_ent),side_number(_side_number),sense(_sense),offset(_offset),brother_side_number(-1) {};
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
      const_mem_fun<SideNumber,EntityType,&SideNumber::get_ent_type>,
      member<SideNumber,int,&SideNumber::side_number> >
    >,
    ordered_non_unique<
      const_mem_fun<SideNumber,EntityType,&SideNumber::get_ent_type>
    >
  > > SideNumber_multiIndex;

/**
 * \brief this struct keeps basic methods for moab entity
 * \ingroup ent_multi_indices

  \todo BasicEntity in should be linked to directly to MoAB data structures
  such that connectivity and nodal coordinates could be quickly accessed,
  without need of using native MoAB functions.

 */
struct BasicEntity {
  EntityHandle ent;
  int owner_proc;
  EntityHandle moab_owner_handle;
  unsigned char *pstatus_val_ptr;
  int *sharing_procs_ptr;
  EntityHandle *sharing_handlers_ptr;

  BasicEntity();
  BasicEntity(Interface &moab,const EntityHandle _ent);

  PetscErrorCode iterateBasicEntity(
    EntityHandle _ent,
    int _owner_proc,
    EntityHandle _moab_owner_handle,
    unsigned char *_pstatus_val_ptr,
    int *_sharing_procs_ptr,
    EntityHandle *_sharing_handlers_ptr
  );

  /// get entity type
  inline EntityType get_ent_type() const { return (EntityType)((ent&MB_TYPE_MASK)>>MB_ID_WIDTH); }

  /// get entity id
  inline EntityID get_ent_id() const { return (EntityID)(ent&MB_ID_MASK); };

  /** \brief maob partitioning owner handle
    */
  inline EntityHandle get_owner_ent() const { return moab_owner_handle; }

  /** \brife moab get owner proc
    */
  inline EntityHandle get_owner_proc() const { return owner_proc; }

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
  inline unsigned char get_pstatus() const { return *pstatus_val_ptr; }

  /** \berief get shared processors

  Returning list to shared processors. Lists end with -1. Returns NULL if not
  sharing processors.

  DO NOT MODIFY LIST.

\code
  BasicEntity *ent_ptr = BasicEntity(moan,entity_handle);
  for(int proc = 0; proc<MAX_SHARING_PROCS && -1 != ent_ptr->get_sharing_procs_ptr[proc]; proc++) {
      if(ent_ptr->get_sharing_procs_ptr[proc] == -1) {
	// End of the list
	break;
      }
      int sharing_proc = ent_ptr->get_sharing_procs_ptr[proc];
      EntityHandle sharing_ent = ent_ptr->get_sharing_handlers_ptr[proc];
      if(!(ent_ptr->get_pstatus()&PSTATUS_MULTISHARED)) {
	break;
      }
    }
\endcode

    */
  inline int* get_sharing_procs_ptr() const { return sharing_procs_ptr; }

  /** \berief get sharid entity handlers

  Returning list to shared entity hanlders. Use it with get_sharing_procs_ptr()

  DO NOT MODIFY LIST.

\code
  BasicEntity *ent_ptr = BasicEntity(moan,entity_handle);
  for(int proc = 0; proc<MAX_SHARING_PROCS && -1 != ent_ptr->get_sharing_procs_ptr[proc]; proc++) {
      if(ent_ptr->get_sharing_procs_ptr[proc] == -1) {
	// End of the list
	break;
      }
      int sharing_proc = ent_ptr->get_sharing_procs_ptr[proc];
      EntityHandle sharing_ent = ent_ptr->get_sharing_handlers_ptr[proc];
      if(!(ent_ptr->get_pstatus()&PSTATUS_MULTISHARED)) {
	break;
      }
    }
\endcode

    */
  inline EntityHandle* get_sharing_handlers_ptr() const { return sharing_handlers_ptr; }

};

/**
 * \brief struct keeps handle to refined handle.
 * \ingroup ent_multi_indices

  \todo th_RefType "_RefType" is set as two integers, need to be fixed, it is
  waste of space.

 */
struct RefEntity: public BasicEntity {
  EntityHandle *tag_parent_ent;
  BitRefLevel *tag_BitRefLevel;

  RefEntity();
  RefEntity(Interface &moab,const EntityHandle _ent);

  PetscErrorCode iterateRefEntity(
    EntityHandle _ent,
    int _owner_proc,
    EntityHandle _moab_owner_handle,
    unsigned char *_pstatus_val_ptr,
    int *_sharing_procs_ptr,
    EntityHandle *_sharing_handlers_ptr,
    EntityHandle *_tag_parent_ent,
    BitRefLevel *_tag_BitRefLevel
  );

  static PetscErrorCode getPatentEnt(Interface &moab,Range ents,std::vector<EntityHandle> vec_patent_ent);
  static PetscErrorCode getBitRefLevel(Interface &moab,Range ents,std::vector<BitRefLevel> vec_bit_ref_level);

  /// get entity
  inline EntityHandle get_ref_ent() const { return ent; }
  /// get patent entity
  inline EntityType get_parent_ent_type() const {
    if(tag_parent_ent == NULL) return MBMAXTYPE;
    if(*tag_parent_ent == 0)  return MBMAXTYPE;
    return (EntityType)((*tag_parent_ent&MB_TYPE_MASK)>>MB_ID_WIDTH);
  }
  /// get entity ref bit refinment signature
  inline const BitRefLevel& get_BitRefLevel() const { return *tag_BitRefLevel; }
  /// get parent entity, i.e. entity form one refinment level up
  inline EntityHandle get_parent_ent() const {
    if(tag_parent_ent == NULL) return 0;
    return *tag_parent_ent;
  }

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

  inline EntityHandle get_ref_ent() const { return this->sPtr->get_ref_ent(); }
  inline EntityHandle get_parent_ent() const { return this->sPtr->get_parent_ent(); }
  inline const BitRefLevel& get_BitRefLevel() const { return this->sPtr->get_BitRefLevel(); }
  inline EntityType get_ent_type() const { return this->sPtr->get_ent_type(); };
  inline EntityType get_parent_ent_type() const { return this->sPtr->get_parent_ent_type(); };
  inline EntityID get_ent_id() const { return this->sPtr->get_ent_id(); };
  inline unsigned char get_pstatus() const { return this->sPtr->get_pstatus(); }
  inline EntityHandle get_owner_ent() const { return this->sPtr->get_owner_ent(); }
  inline EntityHandle get_owner_proc() const { return this->sPtr->get_owner_proc(); }
  inline int* get_sharing_procs_ptr() const { return this->sPtr->get_sharing_procs_ptr(); }
  inline EntityHandle* get_sharing_handlers_ptr() const { return this->sPtr->get_sharing_handlers_ptr(); }
  virtual ~interface_RefEntity() {}

  inline const boost::shared_ptr<T> get_RefEntity_ptr() {
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
      tag<Ent_ParallelStatus>, const_mem_fun<RefEntity::BasicEntity,unsigned char,&RefEntity::get_pstatus> >,
    ordered_non_unique<
      tag<Ent_Ent_mi_tag>, const_mem_fun<RefEntity,EntityHandle,&RefEntity::get_parent_ent> >,
    ordered_non_unique<
      tag<EntType_mi_tag>, const_mem_fun<RefEntity::BasicEntity,EntityType,&RefEntity::get_ent_type> >,
    ordered_non_unique<
      tag<ParentEntType_mi_tag>, const_mem_fun<RefEntity,EntityType,&RefEntity::get_parent_ent_type> >,
    ordered_non_unique<
      tag<Composite_EntType_and_ParentEntType_mi_tag>,
      composite_key<
      	RefEntity,
      	const_mem_fun<RefEntity::BasicEntity,EntityType,&RefEntity::get_ent_type>,
      	const_mem_fun<RefEntity,EntityType,&RefEntity::get_parent_ent_type> > >,
    ordered_non_unique<
      tag<Composite_Ent_And_ParentEntType_mi_tag>,
      composite_key<
      	RefEntity,
      	const_mem_fun<RefEntity,EntityHandle,&RefEntity::get_parent_ent>,
      	const_mem_fun<RefEntity::BasicEntity,EntityType,&RefEntity::get_ent_type> > >
  > > RefEntity_multiIndex;

/** \brief multi-index view of RefEntity by parent entity
  \ingroup ent_multi_indices
*/
typedef multi_index_container<
  boost::shared_ptr<RefEntity>,
  indexed_by<
    hashed_unique<
      const_mem_fun<RefEntity,EntityHandle,&RefEntity::get_parent_ent> >,
    hashed_unique<
      tag<Composite_EntType_and_ParentEntType_mi_tag>,
    composite_key<
	    boost::shared_ptr<RefEntity>,
	    const_mem_fun<RefEntity,EntityHandle,&RefEntity::get_ref_ent>,
	    const_mem_fun<RefEntity,EntityHandle,&RefEntity::get_parent_ent> > >
  > > RefEntity_multiIndex_view_by_parent_entity;

/** \brief ref mofem entity, remove parent
 * \ingroup ent_multi_indices
 */
struct RefEntity_change_remove_parent {
  Interface &mOab;
  Tag th_RefParentHandle;
  ErrorCode rval;
  RefEntity_change_remove_parent(Interface &moab): mOab(moab) {
    rval = mOab.tag_get_handle("_RefParentHandle",th_RefParentHandle); MOAB_THROW(rval);
  }
  void operator()(boost::shared_ptr<RefEntity> &e) {
    rval = mOab.tag_delete_data(th_RefParentHandle,&e->ent,1); MOAB_THROW(rval);
    rval = mOab.tag_get_by_ptr(
      th_RefParentHandle,&e->ent,1,(const void **)&(e->tag_parent_ent)
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
  Interface &mOab;
  EntityHandle pArent;
  Tag th_RefParentHandle;
  ErrorCode rval;
  RefEntity_change_parent(Interface &moab,EntityHandle parent): mOab(moab),pArent(parent) {
    rval = mOab.tag_get_handle("_RefParentHandle",th_RefParentHandle); MOAB_THROW(rval);
  }
  void operator()(boost::shared_ptr<RefEntity> &e) {
    rval = mOab.tag_get_by_ptr(
      th_RefParentHandle,&e->ent,1,(const void **)&(e->tag_parent_ent)
    ); MOAB_THROW(rval);
    *(e->tag_parent_ent) = pArent;
  }
};

/** \brief ref mofem entity, left shift
  * \ingroup ent_multi_indices
  */
struct RefEntity_change_left_shift {
  int shift;
  RefEntity_change_left_shift(const int _shift): shift(_shift) {};
  void operator()(boost::shared_ptr<RefEntity> &e) { (*e->tag_BitRefLevel)<<=shift;  };
};

/** \brief ref mofem entity, right shift
 * \ingroup ent_multi_indices
  */
struct RefEntity_change_right_shift {
  int shift;
  RefEntity_change_right_shift(const int _shift): shift(_shift) {};
  void operator()(boost::shared_ptr<RefEntity> &e) { *(e->tag_BitRefLevel)>>=shift;  };
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices
  */
struct RefEntity_change_add_bit {
  BitRefLevel bit;
  RefEntity_change_add_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(boost::shared_ptr<RefEntity> &e) {
    bit |= *(e->tag_BitRefLevel);
    *(e->tag_BitRefLevel) = bit;
  }
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices
  */
struct RefEntity_change_and_bit {
  BitRefLevel bit;
  RefEntity_change_and_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(boost::shared_ptr<RefEntity> &e) {
    bit &= *(e->tag_BitRefLevel);
    *(e->tag_BitRefLevel) = bit;
  }
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices
  */
struct RefEntity_change_xor_bit {
  BitRefLevel bit;
  RefEntity_change_xor_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(boost::shared_ptr<RefEntity> &e) {
    bit ^= *(e->tag_BitRefLevel);
    *(e->tag_BitRefLevel) = bit;
  }
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices
  */
struct RefEntity_change_set_bit {
  BitRefLevel bit;
  RefEntity_change_set_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(boost::shared_ptr<RefEntity> &e) {
    *(e->tag_BitRefLevel) = bit;
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
    (*(e->tag_BitRefLevel))[n] = b;
  }
};

/**
  * \brief struct keeps handle to entity in the field.
  * \ingroup ent_multi_indices
  */
struct MoFEMEntity:
  public
  interface_Field<Field>,
  interface_RefEntity<RefEntity> {

  typedef interface_Field<Field> interface_type_Field;
  typedef interface_RefEntity<RefEntity> interface_type_RefEntity;
  const ApproximationOrder* tag_order_data;
  const FieldData* tag_FieldData;
  int tag_FieldData_size;
  const ApproximationOrder* tag_dof_order_data;
  const FieldCoefficientsNumber* tag_dof_rank_data;
  LocalUId local_uid;
  GlobalUId global_uid;
  MoFEMEntity(
    Interface &moab,
    const boost::shared_ptr<Field> field_ptr,
    const boost::shared_ptr<RefEntity> ref_ent_ptr
  );
  ~MoFEMEntity();
  inline EntityHandle get_ent() const { return get_ref_ent(); }
  inline int get_nb_dofs_on_ent() const { return tag_FieldData_size/sizeof(FieldData); }
  // inline FieldData* get_ent_FieldData() const { return const_cast<FieldData*>(tag_FieldData); }

  inline VectorAdaptor get_ent_FieldData() const {
    int size = get_nb_dofs_on_ent();
    double* ptr = const_cast<FieldData*>(tag_FieldData);
    return VectorAdaptor(size,ublas::shallow_array_adaptor<FieldData>(size,ptr));
  }

  inline int get_order_nb_dofs(int order) const { return (this->sFieldPtr->forder_table[get_ent_type()])(order); }
  inline int get_order_nb_dofs_diff(int order) const { return get_order_nb_dofs(order)-get_order_nb_dofs(order-1); }

  inline ApproximationOrder get_max_order() const { return *((ApproximationOrder*)tag_order_data); }
  const LocalUId& get_local_unique_id() const { return local_uid; }
  LocalUId get_local_unique_id_calculate() const {
    char bit_number = get_bit_number();
    assert(bit_number<32);
    LocalUId _uid_ = (UId)0;
    _uid_ |= (UId)sPtr->ent;
    _uid_ |= (UId)bit_number << 8*sizeof(EntityHandle);
    return _uid_;
  }
  const GlobalUId& get_global_unique_id() const { return global_uid; }
  GlobalUId get_global_unique_id_calculate() const {
    char bit_number = get_bit_number();
    assert(bit_number<32);
    assert(sPtr->owner_proc<1024);
    GlobalUId _uid_ = (UId)0;
    _uid_ |= (UId)sPtr->moab_owner_handle;
    _uid_ |= (UId)bit_number << 8*sizeof(EntityHandle);
    _uid_ |= (UId)sPtr->owner_proc << 5+8*sizeof(EntityHandle);
    return _uid_;
  }
  friend std::ostream& operator<<(std::ostream& os,const MoFEMEntity& e);

  inline const boost::shared_ptr<RefEntity> get_RefEntity_ptr() {
    return this->sPtr;
  }
  inline const boost::shared_ptr<Field> get_Field_ptr() const {
    return this->sFieldPtr;
  }


};

/**
 * \brief interface to MoFEMEntity
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
  inline EntityHandle get_ent() const { return this->sPtr->get_ent(); }

  inline int get_nb_dofs_on_ent() const { return this->sPtr->get_nb_dofs_on_ent(); }
  inline VectorAdaptor get_ent_FieldData() const { return this->sPtr->get_FieldData(); }
  inline int get_order_nb_dofs(int order) const { return this->sFieldPtr->get_order_nb_dofs(order); }
  inline int get_order_nb_dofs_diff(int order) const { return this->sPtr->get_order_nb_dofs_diff(order); }
  inline ApproximationOrder get_max_order() const { return this->sPtr->get_max_order(); }
  inline const LocalUId& get_local_unique_id() const { return this->sPtr->get_local_unique_id(); }
  inline const GlobalUId& get_global_unique_id() const { return this->sPtr->get_global_unique_id(); }

  inline const boost::shared_ptr<MoFEMEntity> get_MoFEMEntity_ptr() const { return this->sPtr; };
  inline const boost::shared_ptr<RefEntity> get_RefEntity_ptr() {
    return this->sPtr->get_RefEntity_ptr();
  }
  inline const boost::shared_ptr<Field> get_Field_ptr() const {
    return this->sFieldPtr->get_Field_ptr();
  }

};

/**
 * \brief structure to chane MoFEMEntity order
 * \ingroup ent_multi_indices
 */
struct MoFEMEntity_change_order {
  Interface& moab;
  ApproximationOrder order;
  std::vector<FieldData> data;
  std::vector<ApproximationOrder> data_dof_order;
  std::vector<FieldCoefficientsNumber> data_dof_rank;
  MoFEMEntity_change_order(Interface& _moab,ApproximationOrder _order):
  moab(_moab),
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
      tag<Ent_ParallelStatus>, const_mem_fun<MoFEMEntity::interface_type_RefEntity,unsigned char,&MoFEMEntity::get_pstatus> >,
    ordered_non_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_Field,const BitFieldId&,&MoFEMEntity::get_id>, LtBit<BitFieldId> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_Field,boost::string_ref,&MoFEMEntity::get_name_ref> >,
    hashed_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::get_ent> >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_mi_tag>,
      composite_key<
      	MoFEMEntity,
      	const_mem_fun<MoFEMEntity::interface_type_Field,boost::string_ref,&MoFEMEntity::get_name_ref>,
      	const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::get_ent>
      > >
  > > MoFEMEntity_multiIndex;

  typedef multi_index_container<
    boost::shared_ptr<MoFEMEntity>,
    indexed_by<
      hashed_non_unique<
        tag<Ent_mi_tag>, const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::get_ent> >
  > > MoFEMEntity_multiIndex_ent_view;

}

#endif // __ENTSMULTIINDICES_HPP__

/***************************************************************************//**
 * \defgroup ent_multi_indices Entities structures and multi-indices
 * \ingroup mofem
 ******************************************************************************/
