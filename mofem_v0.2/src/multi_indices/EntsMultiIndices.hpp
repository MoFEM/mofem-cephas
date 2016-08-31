/** \file EntsMultiIndices.hpp
 * \brief Myltindex containes, for mofem enitities data structures and other low-level functions 
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
  SideNumber,
  indexed_by<
    hashed_unique<
      member<SideNumber,EntityHandle,&SideNumber::ent> >,
    ordered_non_unique<
      composite_key<
	SideNumber,
	const_mem_fun<SideNumber,EntityType,&SideNumber::get_ent_type>,
	member<SideNumber,int,&SideNumber::side_number> > >,
    ordered_non_unique<
      const_mem_fun<SideNumber,EntityType,&SideNumber::get_ent_type> >
  > > SideNumber_multiIndex; 

/** 
 * \brief this struct keeps basic methods for moab entity
 * \ingroup ent_multi_indices 

  \bug BasicMoFEMEntity in should be linked to directly to MoAB data structures such
  that connectivity and nodal coordinates could be quickly accessed, without
  need of using native MoAB functions.

 */
struct BasicMoFEMEntity {
  EntityHandle ent;
  int owner_proc;
  EntityHandle moab_owner_handle;
  unsigned char *pstatus_val_ptr;
  int *sharing_procs_ptr;
  EntityHandle *sharing_handlers_ptr;

  BasicMoFEMEntity(Interface &moab,const EntityHandle _ent);

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

  /** \berief get sharid processors

  Returning list to shared processors. Lists end with -1. Returns NULL if not
  sharing processors.

  DO NOT MODIFY LIST. 

\code
  BasicMoFEMEntity *ent_ptr = BasicMoFEMEntity(moan,entity_handle);
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
  BasicMoFEMEntity *ent_ptr = BasicMoFEMEntity(moan,entity_handle);
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

  \bug th_RefType "_RefType" is set as two integers, need to be fixed, it is waset of space.

 */
struct RefMoFEMEntity: public BasicMoFEMEntity {
  EntityHandle *tag_parent_ent;
  int tag_parent_ent_size;
  BitRefLevel *tag_BitRefLevel;
  RefMoFEMEntity(Interface &moab,const EntityHandle _ent);
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
  const RefMoFEMEntity* get_RefMoFEMEntity_ptr() { return this; }
  friend ostream& operator<<(ostream& os,const RefMoFEMEntity& e);
};


/** 
 * \brief interface to RefMoFEMEntity
 * \ingroup ent_multi_indices 
 */
template <typename T>
struct interface_RefMoFEMEntity {
  const T *ref_ptr;
  interface_RefMoFEMEntity(const T *_ref_ptr): ref_ptr(_ref_ptr) {}
  inline EntityHandle get_ref_ent() const { return ref_ptr->get_ref_ent(); }
  inline EntityHandle get_parent_ent() const { return ref_ptr->get_parent_ent(); }
  inline const BitRefLevel& get_BitRefLevel() const { return ref_ptr->get_BitRefLevel(); }
  inline EntityType get_ent_type() const { return ref_ptr->get_ent_type(); };
  inline EntityType get_parent_ent_type() const { return ref_ptr->get_parent_ent_type(); };
  inline EntityID get_ent_id() const { return ref_ptr->get_ent_id(); };
  inline const RefMoFEMEntity* get_RefMoFEMEntity_ptr() { return ref_ptr->get_RefMoFEMEntity_ptr(); }
  inline unsigned char get_pstatus() const { return ref_ptr->get_pstatus(); }
  inline EntityHandle get_owner_ent() const { return ref_ptr->get_owner_ent(); }
  inline EntityHandle get_owner_proc() const { return ref_ptr->get_owner_proc(); }
  inline int* get_sharing_procs_ptr() const { return ref_ptr->get_sharing_procs_ptr(); }
  inline EntityHandle* get_sharing_handlers_ptr() const { return ref_ptr->get_sharing_handlers_ptr(); }
  virtual ~interface_RefMoFEMEntity() {}
};

/**
 * \typedef RefMoFEMEntity_multiIndex
 * type multiIndex container for RefMoFEMEntity
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
  RefMoFEMEntity,
  indexed_by<
    hashed_unique<
      tag<Ent_mi_tag>, member<RefMoFEMEntity::BasicMoFEMEntity,EntityHandle,&RefMoFEMEntity::ent> >,
    hashed_non_unique<
      tag<Ent_Owner_mi_tag>, member<RefMoFEMEntity::BasicMoFEMEntity,EntityHandle,&RefMoFEMEntity::moab_owner_handle> >,
    ordered_non_unique<
      tag<Proc_mi_tag>, member<RefMoFEMEntity::BasicMoFEMEntity,int,&RefMoFEMEntity::owner_proc> >,
    ordered_non_unique<
      tag<Ent_ParallelStatus>, const_mem_fun<RefMoFEMEntity::BasicMoFEMEntity,unsigned char,&RefMoFEMEntity::get_pstatus> >,
    ordered_non_unique<
      tag<Ent_Ent_mi_tag>, const_mem_fun<RefMoFEMEntity,EntityHandle,&RefMoFEMEntity::get_parent_ent> >,
    ordered_non_unique<
      tag<EntType_mi_tag>, const_mem_fun<RefMoFEMEntity::BasicMoFEMEntity,EntityType,&RefMoFEMEntity::get_ent_type> >,
    ordered_non_unique<
      tag<ParentEntType_mi_tag>, const_mem_fun<RefMoFEMEntity,EntityType,&RefMoFEMEntity::get_parent_ent_type> >,
    ordered_non_unique<
      tag<Composite_EntType_and_ParentEntType_mi_tag>, 
      composite_key<
	RefMoFEMEntity,
	const_mem_fun<RefMoFEMEntity::BasicMoFEMEntity,EntityType,&RefMoFEMEntity::get_ent_type>,
	const_mem_fun<RefMoFEMEntity,EntityType,&RefMoFEMEntity::get_parent_ent_type> > >,
    ordered_non_unique<
      tag<Composite_Ent_And_ParentEntType_mi_tag>, 
      composite_key<
	RefMoFEMEntity,
	const_mem_fun<RefMoFEMEntity,EntityHandle,&RefMoFEMEntity::get_parent_ent>,
	const_mem_fun<RefMoFEMEntity::BasicMoFEMEntity,EntityType,&RefMoFEMEntity::get_ent_type> > >
  > > RefMoFEMEntity_multiIndex;

/** \brief multi-index view of RefMoFEMEntity by parent entity
  \ingroup ent_multi_indices
*/
typedef multi_index_container<
  const RefMoFEMEntity*,
  indexed_by<
    hashed_unique<
      const_mem_fun<RefMoFEMEntity,EntityHandle,&RefMoFEMEntity::get_parent_ent> >,
    hashed_unique<
      tag<Composite_EntType_and_ParentEntType_mi_tag>,
      composite_key<
	const RefMoFEMEntity*,
	const_mem_fun<RefMoFEMEntity,EntityHandle,&RefMoFEMEntity::get_ref_ent>,
	const_mem_fun<RefMoFEMEntity,EntityHandle,&RefMoFEMEntity::get_parent_ent> > >
  > > RefMoFEMEntity_multiIndex_view_by_parent_entity;

/** \brief ref mofem entity, remove parent
 * \ingroup ent_multi_indices 
 */
struct RefMoFEMEntity_change_remove_parent {
  Interface &mOab;
  Tag th_RefParentHandle;
  ErrorCode rval;
  RefMoFEMEntity_change_remove_parent(Interface &moab): mOab(moab) {
    rval = mOab.tag_get_handle("_RefParentHandle",th_RefParentHandle); CHKERR_THROW(rval);
  }
  void operator()(RefMoFEMEntity &e) { 
    rval = mOab.tag_delete_data(th_RefParentHandle,&e.ent,1); CHKERR_THROW(rval);
    rval = mOab.tag_get_by_ptr(th_RefParentHandle,&e.ent,1,(const void **)&(e.tag_parent_ent)); CHKERR_THROW(rval);
  }
};

/** \brief change parent
  * \ingroup ent_multi_indices 
  *
  * Use this function with care. Some other multi-indices can deponent on this.

  Known dependent multi-indices (verify if that list is full): 
  - RefMoFEMEntity_multiIndex
  - RefMoFEMElement_multiIndex

  */
struct RefMoFEMEntity_change_parent {
  Interface &mOab;
  EntityHandle pArent;
  Tag th_RefParentHandle;
  ErrorCode rval;
  RefMoFEMEntity_change_parent(Interface &moab,EntityHandle parent): mOab(moab),pArent(parent) {
    rval = mOab.tag_get_handle("_RefParentHandle",th_RefParentHandle); CHKERR_THROW(rval);
  }
  void operator()(RefMoFEMEntity &e) { 
    rval = mOab.tag_get_by_ptr(th_RefParentHandle,&e.ent,1,(const void **)&(e.tag_parent_ent)); CHKERR_THROW(rval);
    *(e.tag_parent_ent) = pArent;
  }
};

/** \brief ref mofem entity, left shift
  * \ingroup ent_multi_indices 
  */
struct RefMoFEMEntity_change_left_shift {
  int shift;
  RefMoFEMEntity_change_left_shift(const int _shift): shift(_shift) {};
  void operator()(RefMoFEMEntity &e) { (*e.tag_BitRefLevel)<<=shift;  };
};

/** \brief ref mofem entity, right shift
 * \ingroup ent_multi_indices 
  */
struct RefMoFEMEntity_change_right_shift {
  int shift;
  RefMoFEMEntity_change_right_shift(const int _shift): shift(_shift) {};
  void operator()(RefMoFEMEntity &e) { (*e.tag_BitRefLevel)>>=shift;  };
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices 
  */
struct RefMoFEMEntity_change_add_bit {
  BitRefLevel bit;
  RefMoFEMEntity_change_add_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(RefMoFEMEntity &e) { 
    bit |= *(e.tag_BitRefLevel); 
    *e.tag_BitRefLevel = bit;
  }
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices 
  */
struct RefMoFEMEntity_change_and_bit {
  BitRefLevel bit;
  RefMoFEMEntity_change_and_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(RefMoFEMEntity &e) { 
    bit &= *(e.tag_BitRefLevel); 
    *e.tag_BitRefLevel = bit;
  }
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices 
  */
struct RefMoFEMEntity_change_xor_bit {
  BitRefLevel bit;
  RefMoFEMEntity_change_xor_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(RefMoFEMEntity &e) { 
    bit ^= *(e.tag_BitRefLevel); 
    *e.tag_BitRefLevel = bit;
  }
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices 
  */
struct RefMoFEMEntity_change_set_bit {
  BitRefLevel bit;
  RefMoFEMEntity_change_set_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(RefMoFEMEntity &e) { 
    *e.tag_BitRefLevel = bit;
  }
};

/** \brief ref mofem entity, change bit
  * \ingroup ent_multi_indices 
  */
struct RefMoFEMEntity_change_set_nth_bit {
  int n;
  bool b;
  RefMoFEMEntity_change_set_nth_bit(const int _n,bool _b): n(_n),b(_b) {};
  void operator()(RefMoFEMEntity &e) { 
    (*e.tag_BitRefLevel)[n] = b;
  }
};

/**
  * \brief struct keeps handle to entity in the field.
  * \ingroup ent_multi_indices 
  */
struct MoFEMEntity: public interface_MoFEMField<MoFEMField>, interface_RefMoFEMEntity<RefMoFEMEntity> {
  typedef interface_MoFEMField<MoFEMField> interface_type_MoFEMField;
  typedef interface_RefMoFEMEntity<RefMoFEMEntity> interface_type_RefMoFEMEntity;
  const RefMoFEMEntity *ref_mab_ent_ptr;
  const ApproximationOrder* tag_order_data;
  const FieldData* tag_FieldData;
  int tag_FieldData_size;
  const ApproximationOrder* tag_dof_order_data;
  const ApproximationRank* tag_dof_rank_data;
  LocalUId local_uid;
  GlobalUId global_uid;
  MoFEMEntity(Interface &moab,const MoFEMField *_field_ptr,const RefMoFEMEntity *_ref_mab_ent_ptr);
  ~MoFEMEntity();
  inline EntityHandle get_ent() const { return get_ref_ent(); }
  inline int get_nb_dofs_on_ent() const { return tag_FieldData_size/sizeof(FieldData); }
  inline FieldData* get_ent_FieldData() const { return const_cast<FieldData*>(tag_FieldData); }
  inline int get_order_nb_dofs(int order) const { return (interface_MoFEMField<MoFEMField>::field_ptr->forder_table[get_ent_type()])(order); }
  inline int get_order_nb_dofs_diff(int order) const { return get_order_nb_dofs(order)-get_order_nb_dofs(order-1); }
  inline ApproximationOrder get_max_order() const { return *((ApproximationOrder*)tag_order_data); }
  inline const RefMoFEMEntity* get_RefMoFEMEntity_ptr() const { return ref_mab_ent_ptr; }
  const LocalUId& get_local_unique_id() const { return local_uid; }
  LocalUId get_local_unique_id_calculate() const {
    char bit_number = get_bit_number();
    assert(bit_number<64); //assert(bit_number<32);
    LocalUId _uid_ = (UId)0;
    _uid_ |= (UId)ref_ptr->ent;
    _uid_ |= (UId)bit_number << 8*sizeof(EntityHandle);
    return _uid_;
  }
  const GlobalUId& get_global_unique_id() const { return global_uid; }
  GlobalUId get_global_unique_id_calculate() const {
    char bit_number = get_bit_number();
    assert(bit_number<64); //assert(bit_number<32);
    assert(ref_ptr->owner_proc<1024);
    GlobalUId _uid_ = (UId)0;
    _uid_ |= (UId)ref_ptr->moab_owner_handle;
    _uid_ |= (UId)bit_number << 8*sizeof(EntityHandle);
    _uid_ |= (UId)ref_ptr->owner_proc << 5+8*sizeof(EntityHandle);
    return _uid_;
  }
  const MoFEMEntity* get_MoFEMEntity_ptr() const { return this; };
  friend ostream& operator<<(ostream& os,const MoFEMEntity& e);
};

/**
 * \brief interface to MoFEMEntity
 * \ingroup ent_multi_indices 
 *
 * interface to MoFEMEntity
 */
template <typename T>
struct interface_MoFEMEntity: public interface_MoFEMField<T>,interface_RefMoFEMEntity<RefMoFEMEntity> {
  interface_MoFEMEntity(const T *_ptr): interface_MoFEMField<T>(_ptr),interface_RefMoFEMEntity<RefMoFEMEntity>(_ptr->get_RefMoFEMEntity_ptr()) {};
  inline EntityHandle get_ent() const { return interface_MoFEMField<T>::get_ent(); }
  inline int get_nb_dofs_on_ent() const { return interface_MoFEMField<T>::field_ptr->get_nb_dofs_on_ent(); }
  inline FieldData* get_ent_FieldData() const { return interface_MoFEMField<T>::field_ptr->get_FieldData(); }
  inline int get_order_nb_dofs(int order) const { return interface_MoFEMField<T>::field_ptr->get_order_nb_dofs(order); }
  inline int get_order_nb_dofs_diff(int order) const { return interface_MoFEMField<T>::field_ptr->get_order_nb_dofs_diff(order); }
  inline ApproximationOrder get_max_order() const { return interface_MoFEMField<T>::field_ptr->get_max_order(); }
  inline const LocalUId& get_local_unique_id() const { return interface_MoFEMField<T>::field_ptr->get_local_unique_id(); }
  inline const LocalUId& get_global_unique_id() const { return interface_MoFEMField<T>::field_ptr->get_global_unique_id(); }
  inline const MoFEMEntity* get_MoFEMEntity_ptr() const { return interface_MoFEMField<T>::field_ptr->get_MoFEMEntity_ptr(); };
  inline const RefMoFEMEntity* get_RefMoFEMEntity_ptr() const { return interface_MoFEMField<T>::field_ptr->get_RefMoFEMEntity_ptr(); }
};

/**
 * \brief structure to chane MoFEMEntity order
 * \ingroup ent_multi_indices 
 */
struct MoFEMEntity_change_order {
  Interface& moab;
  ApproximationOrder order;
  vector<FieldData> data;
  vector<ApproximationOrder> data_dof_order;
  vector<ApproximationRank> data_dof_rank;
  MoFEMEntity_change_order(Interface& _moab,ApproximationOrder _order): moab(_moab),order(_order) {};
  void operator()(MoFEMEntity &e);
};

/** 
 * @relates multi_index_container
 * \brief MultiIndex container keeps MoFEMEntity
 * \ingroup ent_multi_indices 
 *
 */
typedef multi_index_container<
  MoFEMEntity,
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>, member<MoFEMEntity,GlobalUId,&MoFEMEntity::global_uid> >,
    ordered_non_unique<
      tag<Ent_ParallelStatus>, const_mem_fun<MoFEMEntity::interface_type_RefMoFEMEntity,unsigned char,&MoFEMEntity::get_pstatus> >,
    ordered_non_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_MoFEMField,const BitFieldId&,&MoFEMEntity::get_id>, LtBit<BitFieldId> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_MoFEMField,boost::string_ref,&MoFEMEntity::get_name_ref> >,
    hashed_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::get_ent> >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_mi_tag>, 
      composite_key<
	MoFEMEntity,
	const_mem_fun<MoFEMEntity::interface_type_MoFEMField,boost::string_ref,&MoFEMEntity::get_name_ref>,
	const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::get_ent>
      > >
  > > MoFEMEntity_multiIndex;

}

#endif // __ENTSMULTIINDICES_HPP__

/***************************************************************************//**
 * \defgroup ent_multi_indices Entities structures and multi-indices
 * \ingroup mofem
 ******************************************************************************/


