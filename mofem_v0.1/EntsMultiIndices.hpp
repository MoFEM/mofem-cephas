/** \file EntsMultiIndices.hpp
 * \brief Myltindex containes, for mofem enitities data structures and other low-level functions 
 * 
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 *
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
 * It can be freely used for educational and research purposes 
 * by other institutions. If you use this softwre pleas cite my work. 
 *
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
 *
 *  \param hashed_unique<
 *     member<SideNumber,EntityHandle,&SideNumber::ent> >,
 *  \param ordered_non_unique<
 *     composite_key<
 *	SideNumber,
 *	const_mem_fun<SideNumber,EntityType,&SideNumber::get_ent_type>,
 *	member<SideNumber,int,&SideNumber::side_number> > >,
 *  \param ordered_non_unique<
 *     const_mem_fun<SideNumber,EntityType,&SideNumber::get_ent_type> >
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
 */
struct BasicMoFEMEntity {
  EntityHandle ent;
  /// \param ent handle to moab entity
  BasicMoFEMEntity(const EntityHandle _ent);
  /// get entity type
  inline EntityType get_ent_type() const { return (EntityType)((ent&MB_TYPE_MASK)>>MB_ID_WIDTH); }
  /// get entity id
  inline EntityID get_ent_id() const { return (EntityID)(ent&MB_ID_MASK); };
};

/** 
 * \brief struct keeps handle to refined handle.
 */
struct RefMoFEMEntity: public BasicMoFEMEntity {
  EntityHandle *tag_parent_ent;
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
 * \typedef RefMoFEMEntity_multiIndex
 * type multiIndex container for RefMoFEMEntity
 *
 * \param hashed_unique MoABEnt_mi_tag 
 * \param ordered_non_unique Meshset_mi_tag 
 * \param ordered_non_unique MoABEnt_MoABEnt_mi_tag
 * \param ordered_non_unique EntType_mi_tag
 * \param ordered_non_unique ParentEntType_mi_tag
 * \param ordered_non_unique Composite_EntityType_And_ParentEntityType_mi_tag
 * \param ordered_non_unique Composite_EntityHandle_And_ParentEntityType_mi_tag
 */
typedef multi_index_container<
  RefMoFEMEntity,
  indexed_by<
    hashed_unique<
      tag<MoABEnt_mi_tag>, member<RefMoFEMEntity::BasicMoFEMEntity,EntityHandle,&RefMoFEMEntity::ent> >,
    ordered_non_unique<
      tag<MoABEnt_MoABEnt_mi_tag>, const_mem_fun<RefMoFEMEntity,EntityHandle,&RefMoFEMEntity::get_parent_ent> >,
    ordered_non_unique<
      tag<EntType_mi_tag>, const_mem_fun<RefMoFEMEntity::BasicMoFEMEntity,EntityType,&RefMoFEMEntity::get_ent_type> >,
    ordered_non_unique<
      tag<ParentEntType_mi_tag>, const_mem_fun<RefMoFEMEntity,EntityType,&RefMoFEMEntity::get_parent_ent_type> >,
    ordered_non_unique<
      tag<Composite_EntType_mi_tag_and_ParentEntType_mi_tag>, 
      composite_key<
	RefMoFEMEntity,
	const_mem_fun<RefMoFEMEntity::BasicMoFEMEntity,EntityType,&RefMoFEMEntity::get_ent_type>,
	const_mem_fun<RefMoFEMEntity,EntityType,&RefMoFEMEntity::get_parent_ent_type> > >,
    ordered_non_unique<
      tag<Composite_EntityType_And_ParentEntityType_mi_tag>, 
      composite_key<
	RefMoFEMEntity,
	const_mem_fun<RefMoFEMEntity::BasicMoFEMEntity,EntityType,&RefMoFEMEntity::get_ent_type>,
	const_mem_fun<RefMoFEMEntity,EntityType,&RefMoFEMEntity::get_parent_ent_type> > >,
    ordered_non_unique<
      tag<Composite_EntityHandle_And_ParentEntityType_mi_tag>, 
      composite_key<
	RefMoFEMEntity,
	const_mem_fun<RefMoFEMEntity,EntityHandle,&RefMoFEMEntity::get_parent_ent>,
	const_mem_fun<RefMoFEMEntity::BasicMoFEMEntity,EntityType,&RefMoFEMEntity::get_ent_type> > >
  > > RefMoFEMEntity_multiIndex;


/// \brief ref mofem entity, remove parent
struct RefMoFEMEntity_change_remove_parent {
  Interface &moab;
  Tag th_RefParentHandle;
  ErrorCode rval;
  RefMoFEMEntity_change_remove_parent(Interface &_moab):moab(_moab) {
    rval = moab.tag_get_handle("_RefParentHandle",th_RefParentHandle); CHKERR_THROW(rval);
  };
  void operator()(RefMoFEMEntity &e) { 
    rval = moab.tag_delete_data(th_RefParentHandle,&e.ent,1); CHKERR_THROW(rval);
    rval = moab.tag_get_by_ptr(th_RefParentHandle,&e.ent,1,(const void **)&(e.tag_parent_ent)); CHKERR_THROW(rval);
  };
};

/// \brief ref mofem entity, left shift
struct RefMoFEMEntity_change_left_shift {
  int shift;
  RefMoFEMEntity_change_left_shift(const int _shift): shift(_shift) {};
  void operator()(RefMoFEMEntity &e) { (*e.tag_BitRefLevel)<<=shift;  };
};

/// \brief ref mofem entity, right shift
struct RefMoFEMEntity_change_right_shift {
  int shift;
  RefMoFEMEntity_change_right_shift(const int _shift): shift(_shift) {};
  void operator()(RefMoFEMEntity &e) { (*e.tag_BitRefLevel)>>=shift;  };
};

/// \brief ref mofem entity, change bit
struct RefMoFEMEntity_change_add_bit {
  BitRefLevel bit;
  RefMoFEMEntity_change_add_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(RefMoFEMEntity &e) { 
    bit |= *(e.tag_BitRefLevel); 
    *e.tag_BitRefLevel = bit;
  }
};

/// \brief ref mofem entity, change bit
struct RefMoFEMEntity_change_and_bit {
  BitRefLevel bit;
  RefMoFEMEntity_change_and_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(RefMoFEMEntity &e) { 
    bit &= *(e.tag_BitRefLevel); 
    *e.tag_BitRefLevel = bit;
  }
};

/// \brief ref mofem entity, change bit
struct RefMoFEMEntity_change_xor_bit {
  BitRefLevel bit;
  RefMoFEMEntity_change_xor_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(RefMoFEMEntity &e) { 
    bit ^= *(e.tag_BitRefLevel); 
    *e.tag_BitRefLevel = bit;
  }
};

/// \brief ref mofem entity, change bit
struct RefMoFEMEntity_change_set_bit {
  BitRefLevel bit;
  RefMoFEMEntity_change_set_bit(const BitRefLevel &_bit): bit(_bit) {};
  void operator()(RefMoFEMEntity &e) { 
    *e.tag_BitRefLevel = bit;
  }
};

/// \brief ref mofem entity, change bit
struct RefMoFEMEntity_change_set_nth_bit {
  int n;
  bool b;
  RefMoFEMEntity_change_set_nth_bit(const int _n,bool _b): n(_n),b(_b) {};
  void operator()(RefMoFEMEntity &e) { 
    (*e.tag_BitRefLevel)[n] = b;
  }
};

/** 
 * \brief interface to RefMoFEMEntity
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
  virtual ~interface_RefMoFEMEntity() {}
};

/**
 * \brief struct keeps handle to entity in the field.
 */
struct MoFEMEntity: public interface_MoFEMField<MoFEMField>, interface_RefMoFEMEntity<RefMoFEMEntity> {
  typedef interface_MoFEMField<MoFEMField> interface_type_MoFEMField;
  const RefMoFEMEntity *ref_mab_ent_ptr;
  const ApproximationOrder* tag_order_data;
  const FieldData* tag_FieldData;
  int tag_FieldData_size;
  const ApproximationOrder* tag_dof_order_data;
  const ApproximationRank* tag_dof_rank_data;
  int (*forder)(int);
  UId uid;
  MoFEMEntity(Interface &moab,const MoFEMField *_FieldData,const RefMoFEMEntity *_ref_mab_ent_ptr);
  ~MoFEMEntity();
  inline EntityHandle get_ent() const { return get_ref_ent(); }
  inline int get_nb_dofs_on_ent() const { return tag_FieldData_size/sizeof(FieldData); }
  inline FieldData* get_ent_FieldData() const { return const_cast<FieldData*>(tag_FieldData); }
  inline int get_order_nb_dofs(int order) const { return forder(order); }
  inline int get_order_nb_dofs_diff(int order) const { return forder(order)-forder(order-1); }
  inline ApproximationOrder get_max_order() const { return *((ApproximationOrder*)tag_order_data); }
  inline const RefMoFEMEntity* get_RefMoFEMEntity_ptr() const { return ref_mab_ent_ptr; }
  const UId& get_unique_id() const { return uid; }
  UId get_unique_id_calculate() const {
    char bit_number = get_bit_number();
    assert(bit_number<=32);
    UId _uid_ = (ref_ptr->ent)|(((UId)bit_number)<<(8*sizeof(EntityHandle)));
    return _uid_;
  }
  const MoFEMEntity* get_MoFEMEntity_ptr() const { return this; };
  friend ostream& operator<<(ostream& os,const MoFEMEntity& e);
};

/**
 * \brief interface to MoFEMEntity
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
  inline const UId& get_unique_id() const { return interface_MoFEMField<T>::field_ptr->get_unique_id(); }
  inline const MoFEMEntity* get_MoFEMEntity_ptr() const { return interface_MoFEMField<T>::field_ptr->get_MoFEMEntity_ptr(); };
  inline const RefMoFEMEntity* get_RefMoFEMEntity_ptr() const { return interface_MoFEMField<T>::field_ptr->get_RefMoFEMEntity_ptr(); }
};

/**
 * \brief structure to chane MoFEMEntity order
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
 *
 * \param ordered_unique<
 *    tag<Unique_mi_tag>, member<MoFEMEntity,UId,&MoFEMEntity::uid> >,
 * \param ordered_non_unique<
 *    tag<BitFieldId_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_MoFEMField,const BitFieldId&,&MoFEMEntity::get_id>, ltbit<BitFieldId> >,
 * \param ordered_non_unique<
 *    tag<FieldName_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_MoFEMField,boost::string_ref,&MoFEMEntity::get_name_ref> >,
 * \param hashed_non_unique<
 *    tag<MoABEnt_mi_tag>, const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::get_ent> >,
 * \param ordered_non_unique<
 *   tag<Composite_Name_And_Ent_mi_tag>, 
 *     composite_key<
 *	MoFEMEntity,
 *	const_mem_fun<MoFEMEntity::interface_type_MoFEMField,boost::string_ref,&MoFEMEntity::get_name_ref>,
 *	const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::get_ent>
 *     > >
 */
typedef multi_index_container<
  MoFEMEntity,
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>, member<MoFEMEntity,UId,&MoFEMEntity::uid> >,
    ordered_non_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_MoFEMField,const BitFieldId&,&MoFEMEntity::get_id>, ltbit<BitFieldId> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_MoFEMField,boost::string_ref,&MoFEMEntity::get_name_ref> >,
    hashed_non_unique<
      tag<MoABEnt_mi_tag>, const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::get_ent> >,
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
