/** \file DofsMultiIndices.hpp
 * \ingroup dof_multi_indices
 * \brief Multi-Index contains, data structures for mofem dofs and other low-level functions
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.

 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.

 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __DOFSMULTIINDICES_HPP__
#define __DOFSMULTIINDICES_HPP__

namespace MoFEM {

/**
 * \brief keeps information about indexed dofs
 * \ingroup dof_multi_indices

 */
struct DofMoFEMEntity: public interface_MoFEMEntity<MoFEMEntity> {

  typedef interface_MoFEMField<MoFEMEntity> interface_type_MoFEMField;
  typedef interface_MoFEMEntity<MoFEMEntity> interface_type_MoFEMEntity;
  typedef interface_RefMoFEMEntity<RefMoFEMEntity> interface_type_RefMoFEMEntity;

  static LocalUId get_local_unique_id_calculate(const DofIdx _dof_,const MoFEMEntity *_ent_ptr_) {
    if(_dof_>=512) THROW_MESSAGE("_dof>=512");
    LocalUId _uid_ = ((UId)_dof_)|((_ent_ptr_->get_local_unique_id())<<9);
    return _uid_;
  }

  static GlobalUId get_global_unique_id_calculate(const DofIdx _dof_,const MoFEMEntity *_ent_ptr_) {
    if(_dof_>=512) THROW_MESSAGE("_dof>=512");
    GlobalUId _uid_ = ((UId)_dof_)|((_ent_ptr_->get_global_unique_id())<<9);
    return _uid_;
  }

  static ShortId get_non_nonunique_short_id(const DofIdx _dof_,const MoFEMEntity *_ent_ptr_) {
    if(_dof_>=512) THROW_MESSAGE("_dof>=512")
    if(sizeof(ShortId) < sizeof(char)+2) THROW_MESSAGE("sizeof(ShortId)< sizeof(char)+2")
    char bit_number = _ent_ptr_->get_bit_number();
    ShortId _uid_ = ((ShortId)_dof_)|(((ShortId)bit_number)<<9);
    return _uid_;
  }

  DofIdx dof;
  bool active;
  LocalUId local_uid;
  GlobalUId global_uid;
  ShortId short_uid;

  DofMoFEMEntity(
    const MoFEMEntity *entity_ptr,
    const ApproximationOrder dof_order,
    const FieldCoefficientsNumber dof_rank,
    const DofIdx _dof
  );

  inline DofIdx get_EntDofIdx() const { return dof; }
  inline FieldData& get_FieldData() const { return const_cast<FieldData&>(field_ptr->tag_FieldData[dof]); }

  /** \brief unique dof id
    */
  inline const LocalUId& get_local_unique_id() const { return local_uid; };
  inline LocalUId get_local_unique_id_calculate() const { return get_local_unique_id_calculate(dof,get_MoFEMEntity_ptr()); }

  inline const GlobalUId& get_global_unique_id() const { return global_uid; };
  // inline GlobalUId get_global_unique_id() const { return global_uid; };
  inline GlobalUId get_global_unique_id_calculate() const { return get_global_unique_id_calculate(dof,get_MoFEMEntity_ptr()); }

  /** \brief get short uid it is unique in combination with entity handle
    *
    * EntityHandle are controlled by MOAB, which is unique in
    * MOAB instance. However two MOAB instances, can have attached different
    * EntityHandles to the same entity.
    *
    * Relation between MoAB EntityHandle can be handled by saving entity handle
    * data into tag, see MB_TYPE_HANDLE. MOAB at time of reading file or
    * creating new MOAB instance, substitute tag value by approbate entity
    * handle.
    *
    * ShortId is created to handle problems related to saving data series, and
    * reading those data using different MoAB instances.
    *
    */
  inline ShortId get_non_nonunique_short_id() const  { return short_uid; }
  inline ShortId get_non_nonunique_short_id_calculate() const { return get_non_nonunique_short_id(dof,get_MoFEMEntity_ptr()); }
  inline EntityHandle get_ent() const { return field_ptr->get_ent(); };
  //inline EntityType get_ent_type() const { return field_ptr->get_ent_type(); };
  inline ApproximationOrder get_dof_order() const {
    return ((ApproximationOrder*)field_ptr->tag_dof_order_data)[dof];
  };

  DEPRECATED inline FieldCoefficientsNumber get_dof_rank() const {
    return ((FieldCoefficientsNumber*)field_ptr->tag_dof_rank_data)[dof];
  };

  /** \brief Get dof coefficient
  */
  inline FieldCoefficientsNumber get_dof_coeff_idx() const {
    return ((FieldCoefficientsNumber*)field_ptr->tag_dof_rank_data)[dof];
  };

  //check if node is active
  inline char get_active() const { return active ? 1 : 0; }
  friend ostream& operator<<(ostream& os,const DofMoFEMEntity& e);

};

/**
 * \brief interface to DofMoFEMEntitys
 * \ingroup dof_multi_indices
 */
template <typename T>
struct interface_DofMoFEMEntity: public interface_MoFEMEntity<T> {
  const boost::shared_ptr<T> sPtr; // FIXME: This will be moved down inherence tree

  interface_DofMoFEMEntity(const boost::shared_ptr<T> sptr):
  interface_MoFEMEntity<T>(&*sptr),
  sPtr(sptr) {

  };

  inline const LocalUId& get_local_unique_id() const { return this->sPtr->get_local_unique_id(); }
  inline const GlobalUId& get_global_unique_id() const { return this->sPtr->get_global_unique_id(); }
  inline ShortId get_non_nonunique_short_id() const { return this->sPtr->get_non_nonunique_short_id(); }
  inline DofIdx get_EntDofIdx() const { return this->sPtr->get_EntDofIdx(); }
  inline FieldData& get_FieldData() const { return this->sPtr->get_FieldData(); }
  inline EntityHandle get_ent() const { return this->sPtr->get_ent(); };
  inline ApproximationOrder get_dof_order() const { return this->sPtr->get_dof_order(); };

  DEPRECATED inline FieldCoefficientsNumber get_dof_rank() const {
    return this->sPtr->get_dof_coeff_idx();
  };

  inline FieldCoefficientsNumber get_dof_coeff_idx() const {
    return this->sPtr->get_dof_coeff_idx();
  };

  inline char get_active() const { return this->sPtr->get_active(); }
  inline const boost::shared_ptr<DofMoFEMEntity> get_DofMoFEMEntity_ptr() const {
    return this->sPtr;
  };

};

/**
 * \brief keeps information about indexed dofs for the problem
 * \ingroup dof_multi_indices
 */
struct NumeredDofMoFEMEntity: public interface_DofMoFEMEntity<DofMoFEMEntity> {
  typedef interface_MoFEMField<DofMoFEMEntity> interface_type_MoFEMField;
  typedef interface_MoFEMEntity<DofMoFEMEntity> interface_type_MoFEMEntity;
  typedef interface_DofMoFEMEntity<DofMoFEMEntity> interface_type_DofMoFEMEntity;
  DofIdx dof_idx;
  DofIdx petsc_gloabl_dof_idx;
  DofIdx petsc_local_dof_idx;
  unsigned int part;
  inline DofIdx get_dof_idx() const { return dof_idx; }
  inline DofIdx get_petsc_gloabl_dof_idx() const { return petsc_gloabl_dof_idx;  }
  inline DofIdx get_petsc_local_dof_idx() const { return petsc_local_dof_idx; }
  inline unsigned int get_part() const { return part;  }
  inline bool get_has_local_index() const { return !signbit(petsc_local_dof_idx); }
  NumeredDofMoFEMEntity(const boost::shared_ptr<DofMoFEMEntity> _DofMoFEMEntity_ptr);
  inline bool operator<(const NumeredDofMoFEMEntity& _dof) const { return (UId)get_global_unique_id()<(UId)_dof.get_global_unique_id(); }
  friend ostream& operator<<(ostream& os,const NumeredDofMoFEMEntity& e);
};

/**
 * \brief interface to NumeredDofMoFEMEntity
 * \ingroup dof_multi_indices
 */
template <typename T>
struct interface_NumeredDofMoFEMEntity: public interface_DofMoFEMEntity<T> {
  interface_NumeredDofMoFEMEntity(const boost::shared_ptr<T> sptr): interface_DofMoFEMEntity<T>(sptr) {};
  inline DofIdx get_dof_idx() const { return this->sPtr->get_dof_idx(); }
  inline DofIdx get_petsc_gloabl_dof_idx() const { return this->sPtr->get_petsc_gloabl_dof_idx();  }
  inline DofIdx get_petsc_local_dof_idx() const { return this->sPtr->get_petsc_local_dof_idx(); }
  inline unsigned int get_part() const { return this->sPtr->get_part();  }
  inline bool get_has_local_index() const { return this->sPtr->get_has_local_index(); }
  inline boost::shared_ptr<NumeredDofMoFEMEntity> get_NumeredDofMoFEMEntity_ptr() const { return this->sPtr; };
};

/**
 * \brief keeps basic information about indexed dofs for the finite element
 */
struct BaseFEDofMoFEMEntity {
  BaseFEDofMoFEMEntity(SideNumber *_side_number_ptr): side_number_ptr(_side_number_ptr) {};
  SideNumber *side_number_ptr;
};

/**
 * \brief keeps information about indexed dofs for the finite element
 * \ingroup dof_multi_indices
 */
struct FEDofMoFEMEntity: public BaseFEDofMoFEMEntity,interface_DofMoFEMEntity<DofMoFEMEntity> {
  typedef interface_MoFEMField<DofMoFEMEntity> interface_type_MoFEMField;
  typedef interface_DofMoFEMEntity<DofMoFEMEntity> interface_type_DofMoFEMEntity;
  typedef interface_RefMoFEMEntity<RefMoFEMEntity> interface_type_RefMoFEMEntity;
  FEDofMoFEMEntity(
    SideNumber *_side_number_ptr,
    const boost::shared_ptr<DofMoFEMEntity> _DofMoFEMEntity_ptr
  );
  FEDofMoFEMEntity(boost::tuple<SideNumber *,const boost::shared_ptr<DofMoFEMEntity> > t);
  friend ostream& operator<<(ostream& os,const FEDofMoFEMEntity& e);
};

/**
 * \brief keeps information about indexed dofs for the finite element
 * \ingroup dof_multi_indices
 */
 struct FENumeredDofMoFEMEntity:
 public
 BaseFEDofMoFEMEntity,
 interface_NumeredDofMoFEMEntity<NumeredDofMoFEMEntity> {
   typedef interface_MoFEMField<NumeredDofMoFEMEntity> interface_type_MoFEMField;
   typedef interface_DofMoFEMEntity<NumeredDofMoFEMEntity> interface_type_DofMoFEMEntity;
   typedef interface_RefMoFEMEntity<RefMoFEMEntity> interface_type_RefMoFEMEntity;
   typedef interface_NumeredDofMoFEMEntity<NumeredDofMoFEMEntity> interface_type_NumeredDofMoFEMEntity;
   FENumeredDofMoFEMEntity(
     SideNumber *_side_number_ptr,
     const boost::shared_ptr<NumeredDofMoFEMEntity> _NumeredDofMoFEMEntity_ptr
   );
   FENumeredDofMoFEMEntity(
     boost::tuple<SideNumber *,const boost::shared_ptr<NumeredDofMoFEMEntity> > t
   );
   friend ostream& operator<<(ostream& os,const FENumeredDofMoFEMEntity& e);
 };

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps DofMoFEMEntity
 * \ingroup dof_multi_indices
 */
typedef multi_index_container<
  boost::shared_ptr<DofMoFEMEntity>,
  indexed_by<
    //uniqe
    ordered_unique<
      tag<Unique_mi_tag>, const_mem_fun<DofMoFEMEntity,const GlobalUId&,&DofMoFEMEntity::get_global_unique_id> >,
    ordered_unique<
      tag<Composite_Ent_and_ShortId_mi_tag>,
        composite_key<
        DofMoFEMEntity,
          const_mem_fun<DofMoFEMEntity,EntityHandle,&DofMoFEMEntity::get_ent>,
          const_mem_fun<DofMoFEMEntity,ShortId,&DofMoFEMEntity::get_non_nonunique_short_id>
        > >,
    ordered_unique<
      tag<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>,
      composite_key<
        DofMoFEMEntity,
          const_mem_fun<DofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&DofMoFEMEntity::get_name_ref>,
          const_mem_fun<DofMoFEMEntity,EntityHandle,&DofMoFEMEntity::get_ent>,
          const_mem_fun<DofMoFEMEntity,DofIdx,&DofMoFEMEntity::get_EntDofIdx>
    > >,
    //non_unique
    ordered_non_unique<
      const_mem_fun<DofMoFEMEntity,char,&DofMoFEMEntity::get_active> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<DofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&DofMoFEMEntity::get_name_ref> >,
    ordered_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<DofMoFEMEntity,EntityHandle,&DofMoFEMEntity::get_ent> >,
    ordered_non_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<DofMoFEMEntity::interface_type_MoFEMField,const BitFieldId&,&DofMoFEMEntity::get_id>, LtBit<BitFieldId> >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_mi_tag>,
      composite_key<
        DofMoFEMEntity,
        const_mem_fun<DofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&DofMoFEMEntity::get_name_ref>,
        const_mem_fun<DofMoFEMEntity,EntityHandle,&DofMoFEMEntity::get_ent>
      > >,
      ordered_non_unique<
        tag<Composite_Name_And_Type_mi_tag>,
      composite_key<
        DofMoFEMEntity,
        const_mem_fun<DofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&DofMoFEMEntity::get_name_ref>,
        const_mem_fun<DofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&DofMoFEMEntity::get_ent_type>
      > >,
        ordered_non_unique<
        tag<Composite_Name_Ent_Order_And_CoeffIdx_mi_tag>,
        composite_key<
        DofMoFEMEntity,
          const_mem_fun<DofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&DofMoFEMEntity::get_name_ref>,
          const_mem_fun<DofMoFEMEntity,EntityHandle,&DofMoFEMEntity::get_ent>,
          const_mem_fun<DofMoFEMEntity,ApproximationOrder,&DofMoFEMEntity::get_dof_order>,
          const_mem_fun<DofMoFEMEntity,FieldCoefficientsNumber,&DofMoFEMEntity::get_dof_coeff_idx>
        > >
  > > DofMoFEMEntity_multiIndex;

/** \brief multi-index view on DofMoFEMEntity by uid
  \ingroup dof_multi_indices
*/
typedef multi_index_container<
  boost::shared_ptr<DofMoFEMEntity>,
  indexed_by<
    ordered_unique<
      member<DofMoFEMEntity,const GlobalUId,&DofMoFEMEntity::global_uid> >
  > > DofMoFEMEntity_multiIndex_uid_view;

/** \brief multi-index view on DofMoFEMEntity activity
  \ingroup dof_multi_indices
*/
typedef multi_index_container<
  boost::shared_ptr<DofMoFEMEntity>,
  indexed_by<
    ordered_unique<
      const_mem_fun<DofMoFEMEntity,const GlobalUId&,&DofMoFEMEntity::get_global_unique_id> >,
    ordered_non_unique<
      const_mem_fun<DofMoFEMEntity,char,&DofMoFEMEntity::get_active> >
  > > DofMoFEMEntity_multiIndex_active_view;

/** \brief multi-index view on DofMoFEMEntity order
  \ingroup dof_multi_indices
*/
typedef multi_index_container<
  boost::shared_ptr<DofMoFEMEntity>,
  indexed_by<
    ordered_non_unique<
      const_mem_fun<DofMoFEMEntity,ApproximationOrder,&DofMoFEMEntity::get_dof_order> >
  > > DofMoFEMEntity_multiIndex_order_view;

/** \brief multi-index view on DofMoFEMEntity type
  \ingroup dof_multi_indices
*/
typedef multi_index_container<
  boost::shared_ptr<DofMoFEMEntity>,
  indexed_by<
    ordered_non_unique<
      const_mem_fun<DofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&DofMoFEMEntity::get_ent_type> >
  > > DofMoFEMEntity_multiIndex_ent_type_view;

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps FEDofMoFEMEntity
 * \ingroup dof_multi_indices

 */
typedef multi_index_container<
  boost::shared_ptr<FEDofMoFEMEntity>,
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>, const_mem_fun<FEDofMoFEMEntity::interface_type_DofMoFEMEntity,const GlobalUId&,&FEDofMoFEMEntity::get_global_unique_id> >,
    ordered_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<FEDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&FEDofMoFEMEntity::get_ent> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FEDofMoFEMEntity::get_name_ref> >,
    ordered_non_unique<
      tag<EntType_mi_tag>, const_mem_fun<FEDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FEDofMoFEMEntity::get_ent_type> >,
    ordered_non_unique<
      tag<Composite_Name_Type_And_Side_Number_mi_tag>,
      composite_key<
	FEDofMoFEMEntity,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FEDofMoFEMEntity::get_name_ref>,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FEDofMoFEMEntity::get_ent_type>,
	  KeyFromKey<
	    member<SideNumber,int,&SideNumber::side_number>,
	    member<FEDofMoFEMEntity::BaseFEDofMoFEMEntity,SideNumber *,&FEDofMoFEMEntity::side_number_ptr>
	  >
      > >,
    ordered_non_unique<
      tag<Composite_Name_And_Type_mi_tag>,
      composite_key<
	FEDofMoFEMEntity,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FEDofMoFEMEntity::get_name_ref>,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FEDofMoFEMEntity::get_ent_type>
	> >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_mi_tag>,
      composite_key<
	FEDofMoFEMEntity,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FEDofMoFEMEntity::get_name_ref>,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&FEDofMoFEMEntity::get_ent>
	> >,
    ordered_non_unique<
      tag<Composite_EntType_and_Space_mi_tag>,
      composite_key<
	FEDofMoFEMEntity,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FEDofMoFEMEntity::get_ent_type>,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,FieldSpace,&FEDofMoFEMEntity::get_space>
	> >
  > > FEDofMoFEMEntity_multiIndex;

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps FENumeredDofMoFEMEntity
 * \ingroup dof_multi_indices
 *
 */
typedef multi_index_container<
  boost::shared_ptr<FENumeredDofMoFEMEntity>,
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>, const_mem_fun<FENumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,const GlobalUId&,&FENumeredDofMoFEMEntity::get_global_unique_id> >,
    ordered_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<FENumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&FENumeredDofMoFEMEntity::get_ent> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<FENumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FENumeredDofMoFEMEntity::get_name_ref> >,
    ordered_non_unique<
      tag<PetscGlobalIdx_mi_tag>, const_mem_fun<FENumeredDofMoFEMEntity::interface_type_NumeredDofMoFEMEntity,DofIdx,&FENumeredDofMoFEMEntity::get_petsc_gloabl_dof_idx> >,
    ordered_non_unique<
      tag<Composite_Name_Type_And_Side_Number_mi_tag>,
      composite_key<
	FENumeredDofMoFEMEntity,
	  const_mem_fun<FENumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FENumeredDofMoFEMEntity::get_name_ref>,
	  const_mem_fun<FENumeredDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FENumeredDofMoFEMEntity::get_ent_type>,
	  KeyFromKey<
	    member<SideNumber,int,&SideNumber::side_number>,
	    member<FENumeredDofMoFEMEntity::BaseFEDofMoFEMEntity,SideNumber*,&FENumeredDofMoFEMEntity::side_number_ptr>
	  >
      > >,
    ordered_non_unique<
      tag<Composite_Name_And_Type_mi_tag>,
      composite_key<
	FENumeredDofMoFEMEntity,
	  const_mem_fun<FENumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FENumeredDofMoFEMEntity::get_name_ref>,
	  const_mem_fun<FENumeredDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FENumeredDofMoFEMEntity::get_ent_type>
	> >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_mi_tag>,
      composite_key<
	FENumeredDofMoFEMEntity,
	  const_mem_fun<FENumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FENumeredDofMoFEMEntity::get_name_ref>,
	  const_mem_fun<FENumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&FENumeredDofMoFEMEntity::get_ent>
	> >
  > > FENumeredDofMoFEMEntity_multiIndex;

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps NumeredDofMoFEMEntity
 * \ingroup dof_multi_indices
 */
typedef multi_index_container<
  boost::shared_ptr<NumeredDofMoFEMEntity>,
  //unique
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>,
      const_mem_fun<NumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,const GlobalUId&,&NumeredDofMoFEMEntity::get_global_unique_id> >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>,
      composite_key<
	      NumeredDofMoFEMEntity,
	      const_mem_fun<NumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&NumeredDofMoFEMEntity::get_name_ref>,
	      const_mem_fun<NumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&NumeredDofMoFEMEntity::get_ent>,
	      const_mem_fun<NumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,DofIdx,&NumeredDofMoFEMEntity::get_EntDofIdx>
    > >,
    //non unique
    ordered_non_unique<
      tag<Idx_mi_tag>, member<NumeredDofMoFEMEntity,DofIdx,&NumeredDofMoFEMEntity::dof_idx> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<NumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&NumeredDofMoFEMEntity::get_name_ref> >,
    ordered_non_unique<
      tag<PetscGlobalIdx_mi_tag>, member<NumeredDofMoFEMEntity,DofIdx,&NumeredDofMoFEMEntity::petsc_gloabl_dof_idx> >,
    ordered_non_unique<
      tag<PetscLocalIdx_mi_tag>, member<NumeredDofMoFEMEntity,DofIdx,&NumeredDofMoFEMEntity::petsc_local_dof_idx> >,
    ordered_non_unique<
      tag<Part_mi_tag>, member<NumeredDofMoFEMEntity,unsigned int,&NumeredDofMoFEMEntity::part> >,
    ordered_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<NumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&NumeredDofMoFEMEntity::get_ent> >,
    ordered_non_unique<
      tag<Order_mi_tag>, const_mem_fun<NumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,ApproximationOrder,&NumeredDofMoFEMEntity::get_dof_order> >,
    ordered_non_unique<
      tag<Composite_Part_And_Oder_mi_tag>,
      composite_key<
	      NumeredDofMoFEMEntity,
	       member<NumeredDofMoFEMEntity,unsigned int,&NumeredDofMoFEMEntity::part>,
	      const_mem_fun<NumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,ApproximationOrder,&NumeredDofMoFEMEntity::get_dof_order>
	    > >,
    ordered_non_unique<
      tag<Composite_Name_Part_And_CoeffIdx_mi_tag>,
      composite_key<
	     NumeredDofMoFEMEntity,
	      const_mem_fun<NumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&NumeredDofMoFEMEntity::get_name_ref>,
	      member<NumeredDofMoFEMEntity,unsigned int,&NumeredDofMoFEMEntity::part>,
	      const_mem_fun<NumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,FieldCoefficientsNumber,&NumeredDofMoFEMEntity::get_dof_coeff_idx>
	    > >,
    ordered_non_unique<
      tag<Composite_Name_And_Part_mi_tag>,
      composite_key<
	     NumeredDofMoFEMEntity,
	     const_mem_fun<NumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&NumeredDofMoFEMEntity::get_name_ref>,
	     member<NumeredDofMoFEMEntity,unsigned int,&NumeredDofMoFEMEntity::part>
	  > >,
    ordered_non_unique<
      tag<Composite_Name_Ent_And_Part_mi_tag>,
      composite_key<
	     NumeredDofMoFEMEntity,
	     const_mem_fun<NumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&NumeredDofMoFEMEntity::get_name_ref>,
	     const_mem_fun<NumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&NumeredDofMoFEMEntity::get_ent>,
	     member<NumeredDofMoFEMEntity,unsigned int,&NumeredDofMoFEMEntity::part>
	  > >,
    ordered_non_unique<
      tag<Composite_Name_And_HasLocalIdx_mi_tag>,
      composite_key<
	     NumeredDofMoFEMEntity,
	     const_mem_fun<NumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&NumeredDofMoFEMEntity::get_name_ref>,
	     const_mem_fun<NumeredDofMoFEMEntity,bool,&NumeredDofMoFEMEntity::get_has_local_index>
	  > >
  > > NumeredDofMoFEMEntity_multiIndex;

typedef multi_index_container<
  boost::shared_ptr<NumeredDofMoFEMEntity>,
  indexed_by<
    ordered_unique<
      const_mem_fun<NumeredDofMoFEMEntity,DofIdx,&NumeredDofMoFEMEntity::get_dof_idx> >
  > > NumeredDofMoFEMEntity_multiIndex_uid_view_ordered;

typedef multi_index_container<
  boost::shared_ptr<NumeredDofMoFEMEntity>,
  indexed_by<
    hashed_unique<
      const_mem_fun<NumeredDofMoFEMEntity,DofIdx,&NumeredDofMoFEMEntity::get_dof_idx> >
  > > NumeredDofMoFEMEntity_multiIndex_uid_view_hashed;

typedef multi_index_container<
  boost::shared_ptr<NumeredDofMoFEMEntity>,
  indexed_by<
    ordered_non_unique<
      const_mem_fun<NumeredDofMoFEMEntity,DofIdx,&NumeredDofMoFEMEntity::get_petsc_local_dof_idx> >
  > > NumeredDofMoFEMEntity_multiIndex_petsc_local_dof_view_ordered_non_unique;

struct DofMoFEMEntity_active_change {
  bool active;
  DofMoFEMEntity_active_change(bool _active);
  void operator()(boost::shared_ptr<DofMoFEMEntity> &_dof_);
};

struct NumeredDofMoFEMEntity_part_change {
  unsigned int part;
  DofIdx petsc_gloabl_dof_idx;
  NumeredDofMoFEMEntity_part_change(const unsigned int _part,const DofIdx _petsc_gloabl_dof_idx):
    part(_part),
    petsc_gloabl_dof_idx(_petsc_gloabl_dof_idx) {};
  void operator()(boost::shared_ptr<NumeredDofMoFEMEntity> &dof) {
    dof->part = part;
    dof->petsc_gloabl_dof_idx = petsc_gloabl_dof_idx;
  }
};

struct NumeredDofMoFEMEntity_local_idx_change {
  DofIdx petsc_local_dof_idx;
  NumeredDofMoFEMEntity_local_idx_change(const DofIdx _petsc_local_dof_idx):
    petsc_local_dof_idx(_petsc_local_dof_idx) {};
  void operator()(boost::shared_ptr<NumeredDofMoFEMEntity> &dof) {
    dof->petsc_local_dof_idx = petsc_local_dof_idx;
  }
};

struct NumeredDofMoFEMEntity_mofem_index_change {
  DofIdx mofem_idx;
  NumeredDofMoFEMEntity_mofem_index_change(const DofIdx _mofem_idx):
    mofem_idx(_mofem_idx) {};
  void operator()(boost::shared_ptr<NumeredDofMoFEMEntity> &dof) {
    dof->dof_idx = mofem_idx;
  }
};

typedef multi_index_container<
  boost::shared_ptr<NumeredDofMoFEMEntity>,
  indexed_by<
    ordered_unique<
      member<NumeredDofMoFEMEntity,const DofIdx,&NumeredDofMoFEMEntity::petsc_gloabl_dof_idx> >
 > > NumeredDofMoFEMEntity_multiIndex_global_index_view;

}
#endif // __DOFSMULTIINDICES_HPP__

/***************************************************************************//**
 * \defgroup dof_multi_indices Dofs structures and multi-indices
 * \ingroup mofem
 ******************************************************************************/
