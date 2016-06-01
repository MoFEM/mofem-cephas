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
struct DofEntity: public interface_MoFEMEntity<MoFEMEntity> {

  typedef interface_Field<MoFEMEntity> interface_type_Field;
  typedef interface_MoFEMEntity<MoFEMEntity> interface_type_MoFEMEntity;
  typedef interface_RefEntity<MoFEMEntity> interface_type_RefEntity;

  static inline GlobalUId get_global_unique_id_calculate(const DofIdx dof,const boost::shared_ptr<MoFEMEntity> ent_ptr) {
    if(dof>=512) THROW_MESSAGE("_dof>=512");
    GlobalUId _uid_ = ((UId)dof)|((ent_ptr->get_global_unique_id())<<9);
    return _uid_;
  }

  static inline ShortId get_non_nonunique_short_id(const DofIdx dof,const boost::shared_ptr<MoFEMEntity> ent_ptr) {
    if(dof>=512) THROW_MESSAGE("_dof>=512")
    if(sizeof(ShortId) < sizeof(char)+2) THROW_MESSAGE("sizeof(ShortId)< sizeof(char)+2")
    char bit_number = ent_ptr->get_bit_number();
    ShortId _uid_ = ((ShortId)dof)|(((ShortId)bit_number)<<9);
    return _uid_;
  }

  bool active;
  ShortId short_uid;

  DofEntity(
    const boost::shared_ptr<MoFEMEntity> entity_ptr,
    const ApproximationOrder dof_order,
    const FieldCoefficientsNumber dof_rank,
    const DofIdx dof
  );

  inline DofIdx get_EntDofIdx() const { return (DofIdx)(short_uid&UID_DOF_MAK); }
  inline FieldData& get_FieldData() const { return const_cast<FieldData&>(this->sPtr->tag_FieldData[get_EntDofIdx()]); }

  /** \brief Get unique dof id
    */
  inline const GlobalUId get_global_unique_id() const { return get_global_unique_id_calculate(get_EntDofIdx(),get_MoFEMEntity_ptr()); };

  // inline GlobalUId get_global_unique_id_calculate(const int dof) const { return get_global_unique_id_calculate(dof,get_MoFEMEntity_ptr()); }

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

  /**
   * \brief Calculate short_uid
   * @param  dof DOF number on the entity
   * @return     short_uid
   */
  inline ShortId get_non_nonunique_short_id_calculate(const int dof) const { return get_non_nonunique_short_id(dof,get_MoFEMEntity_ptr()); }

  inline EntityHandle get_ent() const { return this->sPtr->get_ent(); };
  inline ApproximationOrder get_dof_order() const {
    return ((ApproximationOrder*)this->sPtr->tag_dof_order_data)[get_EntDofIdx()];
  };

  DEPRECATED inline FieldCoefficientsNumber get_dof_rank() const {
    return ((FieldCoefficientsNumber*)this->sPtr->tag_dof_rank_data)[get_EntDofIdx()];
  };

  /** \brief Get dof coefficient
  */
  inline FieldCoefficientsNumber get_dof_coeff_idx() const {
    return ((FieldCoefficientsNumber*)this->sPtr->tag_dof_rank_data)[get_EntDofIdx()];
  };

  //check if node is active
  inline char get_active() const { return active ? 1 : 0; }
  friend std::ostream& operator<<(std::ostream& os,const DofEntity& e);

};

/**
 * \brief interface to DofEntitys
 * \ingroup dof_multi_indices
 */
template <typename T>
struct interface_DofEntity: public interface_MoFEMEntity<T> {

  interface_DofEntity(const boost::shared_ptr<T> sptr):
  interface_MoFEMEntity<T>(sptr) {
  };

  // inline const LocalUId& get_local_unique_id() const { return this->sPtr->get_local_unique_id(); }
  inline const GlobalUId get_global_unique_id() const { return this->sPtr->get_global_unique_id(); }
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
  inline const boost::shared_ptr<DofEntity> get_DofEntity_ptr() const {
    return this->sPtr;
  };

  inline const boost::shared_ptr<MoFEMEntity> get_MoFEMEntity_ptr() const {
    return this->sPtr->get_MoFEMEntity_ptr();
  };


};

/**
 * \brief keeps information about indexed dofs for the problem
 * \ingroup dof_multi_indices
 */
struct NumeredDofEntity: public interface_DofEntity<DofEntity> {
  typedef interface_Field<DofEntity> interface_type_Field;
  typedef interface_MoFEMEntity<DofEntity> interface_type_MoFEMEntity;
  typedef interface_DofEntity<DofEntity> interface_type_DofEntity;
  DofIdx dof_idx;
  DofIdx petsc_gloabl_dof_idx;
  DofIdx petsc_local_dof_idx;
  unsigned int part;
  inline DofIdx get_dof_idx() const { return dof_idx; }
  inline DofIdx get_petsc_gloabl_dof_idx() const { return petsc_gloabl_dof_idx;  }
  inline DofIdx get_petsc_local_dof_idx() const { return petsc_local_dof_idx; }
  inline unsigned int get_part() const { return part;  }
  inline bool get_has_local_index() const { return !signbit(petsc_local_dof_idx); }
  NumeredDofEntity(const boost::shared_ptr<DofEntity> _DofEntity_ptr);
  inline bool operator<(const NumeredDofEntity& _dof) const { return (UId)get_global_unique_id()<(UId)_dof.get_global_unique_id(); }
  friend std::ostream& operator<<(std::ostream& os,const NumeredDofEntity& e);
};

/**
 * \brief interface to NumeredDofEntity
 * \ingroup dof_multi_indices
 */
template <typename T>
struct interface_NumeredDofEntity: public interface_DofEntity<T> {
  interface_NumeredDofEntity(const boost::shared_ptr<T> sptr): interface_DofEntity<T>(sptr) {};
  inline DofIdx get_dof_idx() const { return this->sPtr->get_dof_idx(); }
  inline DofIdx get_petsc_gloabl_dof_idx() const { return this->sPtr->get_petsc_gloabl_dof_idx();  }
  inline DofIdx get_petsc_local_dof_idx() const { return this->sPtr->get_petsc_local_dof_idx(); }
  inline unsigned int get_part() const { return this->sPtr->get_part();  }
  inline bool get_has_local_index() const { return this->sPtr->get_has_local_index(); }
  inline boost::shared_ptr<NumeredDofEntity> get_NumeredDofEntity_ptr() const { return this->sPtr; };
};

/**
 * \brief keeps basic information about indexed dofs for the finite element
 */
struct BaseFEDofEntity {
  BaseFEDofEntity(boost::shared_ptr<SideNumber> side_number_ptr):
  sideNumberPtr(side_number_ptr) {};
  boost::shared_ptr<SideNumber> sideNumberPtr;
};

/**
 * \brief keeps information about indexed dofs for the finite element
 * \ingroup dof_multi_indices
 */
struct FEDofEntity: public BaseFEDofEntity,interface_DofEntity<DofEntity> {
  typedef interface_Field<DofEntity> interface_type_Field;
  typedef interface_DofEntity<DofEntity> interface_type_DofEntity;
  typedef interface_RefEntity<DofEntity> interface_type_RefEntity;
  FEDofEntity(
    boost::shared_ptr<SideNumber> side_number_ptr,
    const boost::shared_ptr<DofEntity> dof_ptr
  );
  FEDofEntity(
    boost::tuple<boost::shared_ptr<SideNumber>,const boost::shared_ptr<DofEntity> > t
  );
  friend std::ostream& operator<<(std::ostream& os,const FEDofEntity& e);
};

/**
 * \brief keeps information about indexed dofs for the finite element
 * \ingroup dof_multi_indices
 */
 struct FENumeredDofEntity:
 public
 BaseFEDofEntity,
 interface_NumeredDofEntity<NumeredDofEntity> {
   typedef interface_Field<NumeredDofEntity> interface_type_Field;
   typedef interface_DofEntity<NumeredDofEntity> interface_type_DofEntity;
   typedef interface_RefEntity<NumeredDofEntity> interface_type_RefEntity;
   typedef interface_NumeredDofEntity<NumeredDofEntity> interface_type_NumeredDofEntity;
   FENumeredDofEntity(
     boost::shared_ptr<SideNumber> side_number_ptr,
     const boost::shared_ptr<NumeredDofEntity> dof_ptr
   );
   FENumeredDofEntity(
     boost::tuple<boost::shared_ptr<SideNumber>,const boost::shared_ptr<NumeredDofEntity> > t
   );
   friend std::ostream& operator<<(std::ostream& os,const FENumeredDofEntity& e);
 };

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps DofEntity
 * \ingroup dof_multi_indices
 */
typedef multi_index_container<
  boost::shared_ptr<DofEntity>,
  indexed_by<
    //uniqe
    ordered_unique<
      tag<Unique_mi_tag>, const_mem_fun<DofEntity,const GlobalUId,&DofEntity::get_global_unique_id> >,
    ordered_unique<
      tag<Composite_Ent_and_ShortId_mi_tag>,
        composite_key<
        DofEntity,
          const_mem_fun<DofEntity,EntityHandle,&DofEntity::get_ent>,
          const_mem_fun<DofEntity,ShortId,&DofEntity::get_non_nonunique_short_id>
        > >,
    ordered_unique<
      tag<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>,
      composite_key<
        DofEntity,
          const_mem_fun<DofEntity::interface_type_Field,boost::string_ref,&DofEntity::get_name_ref>,
          const_mem_fun<DofEntity,EntityHandle,&DofEntity::get_ent>,
          const_mem_fun<DofEntity,DofIdx,&DofEntity::get_EntDofIdx>
    > >,
    //non_unique
    ordered_non_unique<
      const_mem_fun<DofEntity,char,&DofEntity::get_active> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<DofEntity::interface_type_Field,boost::string_ref,&DofEntity::get_name_ref> >,
    ordered_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<DofEntity,EntityHandle,&DofEntity::get_ent> >,
    ordered_non_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<DofEntity::interface_type_Field,const BitFieldId&,&DofEntity::get_id>, LtBit<BitFieldId> >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_mi_tag>,
      composite_key<
        DofEntity,
        const_mem_fun<DofEntity::interface_type_Field,boost::string_ref,&DofEntity::get_name_ref>,
        const_mem_fun<DofEntity,EntityHandle,&DofEntity::get_ent>
      > >,
      ordered_non_unique<
        tag<Composite_Name_And_Type_mi_tag>,
      composite_key<
        DofEntity,
        const_mem_fun<DofEntity::interface_type_Field,boost::string_ref,&DofEntity::get_name_ref>,
        const_mem_fun<DofEntity::interface_type_RefEntity,EntityType,&DofEntity::get_ent_type>
      > >,
        ordered_non_unique<
        tag<Composite_Name_Ent_Order_And_CoeffIdx_mi_tag>,
        composite_key<
        DofEntity,
          const_mem_fun<DofEntity::interface_type_Field,boost::string_ref,&DofEntity::get_name_ref>,
          const_mem_fun<DofEntity,EntityHandle,&DofEntity::get_ent>,
          const_mem_fun<DofEntity,ApproximationOrder,&DofEntity::get_dof_order>,
          const_mem_fun<DofEntity,FieldCoefficientsNumber,&DofEntity::get_dof_coeff_idx>
        > >
  > > DofEntity_multiIndex;

/** \brief multi-index view on DofEntity by uid
  \ingroup dof_multi_indices
*/
typedef multi_index_container<
  boost::shared_ptr<DofEntity>,
  indexed_by<
    ordered_unique<
      const_mem_fun<DofEntity,const GlobalUId,&DofEntity::get_global_unique_id>
    >
  > > DofEntity_multiIndex_uid_view;

/** \brief multi-index view on DofEntity activity
  \ingroup dof_multi_indices
*/
typedef multi_index_container<
  boost::shared_ptr<DofEntity>,
  indexed_by<
    ordered_unique<
      const_mem_fun<DofEntity,const GlobalUId,&DofEntity::get_global_unique_id> >,
    ordered_non_unique<
      const_mem_fun<DofEntity,char,&DofEntity::get_active> >
  > > DofEntity_multiIndex_active_view;

/** \brief multi-index view on DofEntity order
  \ingroup dof_multi_indices
*/
typedef multi_index_container<
  boost::shared_ptr<DofEntity>,
  indexed_by<
    ordered_non_unique<
      const_mem_fun<DofEntity,ApproximationOrder,&DofEntity::get_dof_order> >
  > > DofEntity_multiIndex_order_view;

/** \brief multi-index view on DofEntity type
  \ingroup dof_multi_indices
*/
typedef multi_index_container<
  boost::shared_ptr<DofEntity>,
  indexed_by<
    ordered_non_unique<
      const_mem_fun<DofEntity::interface_type_RefEntity,EntityType,&DofEntity::get_ent_type> >
  > > DofEntity_multiIndex_ent_type_view;

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps FEDofEntity
 * \ingroup dof_multi_indices

 */
typedef multi_index_container<
  boost::shared_ptr<FEDofEntity>,
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>, const_mem_fun<FEDofEntity::interface_type_DofEntity,const GlobalUId,&FEDofEntity::get_global_unique_id> >,
    ordered_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<FEDofEntity::interface_type_DofEntity,EntityHandle,&FEDofEntity::get_ent> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<FEDofEntity::interface_type_Field,boost::string_ref,&FEDofEntity::get_name_ref> >,
    ordered_non_unique<
      tag<EntType_mi_tag>, const_mem_fun<FEDofEntity::interface_type_RefEntity,EntityType,&FEDofEntity::get_ent_type> >,
    ordered_non_unique<
      tag<Composite_Name_Type_And_Side_Number_mi_tag>,
      composite_key<
	FEDofEntity,
	  const_mem_fun<FEDofEntity::interface_type_Field,boost::string_ref,&FEDofEntity::get_name_ref>,
	  const_mem_fun<FEDofEntity::interface_type_RefEntity,EntityType,&FEDofEntity::get_ent_type>,
	  KeyFromKey<
	    member<SideNumber,char,&SideNumber::side_number>,
	    member<FEDofEntity::BaseFEDofEntity,boost::shared_ptr<SideNumber>,&FEDofEntity::sideNumberPtr>
	  >
      > >,
    ordered_non_unique<
      tag<Composite_Name_And_Type_mi_tag>,
      composite_key<
	FEDofEntity,
	  const_mem_fun<FEDofEntity::interface_type_Field,boost::string_ref,&FEDofEntity::get_name_ref>,
	  const_mem_fun<FEDofEntity::interface_type_RefEntity,EntityType,&FEDofEntity::get_ent_type>
	> >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_mi_tag>,
      composite_key<
	FEDofEntity,
	  const_mem_fun<FEDofEntity::interface_type_Field,boost::string_ref,&FEDofEntity::get_name_ref>,
	  const_mem_fun<FEDofEntity::interface_type_DofEntity,EntityHandle,&FEDofEntity::get_ent>
	> >,
    ordered_non_unique<
      tag<Composite_EntType_and_Space_mi_tag>,
      composite_key<
	FEDofEntity,
	  const_mem_fun<FEDofEntity::interface_type_RefEntity,EntityType,&FEDofEntity::get_ent_type>,
	  const_mem_fun<FEDofEntity::interface_type_Field,FieldSpace,&FEDofEntity::get_space>
	> >
  > > FEDofEntity_multiIndex;

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps FENumeredDofEntity
 * \ingroup dof_multi_indices
 *
 */
typedef multi_index_container<
  boost::shared_ptr<FENumeredDofEntity>,
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>, const_mem_fun<FENumeredDofEntity::interface_type_DofEntity,const GlobalUId,&FENumeredDofEntity::get_global_unique_id> >,
    ordered_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<FENumeredDofEntity::interface_type_DofEntity,EntityHandle,&FENumeredDofEntity::get_ent> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<FENumeredDofEntity::interface_type_Field,boost::string_ref,&FENumeredDofEntity::get_name_ref> >,
    ordered_non_unique<
      tag<PetscGlobalIdx_mi_tag>, const_mem_fun<FENumeredDofEntity::interface_type_NumeredDofEntity,DofIdx,&FENumeredDofEntity::get_petsc_gloabl_dof_idx> >,
    ordered_non_unique<
      tag<Composite_Name_Type_And_Side_Number_mi_tag>,
      composite_key<
	FENumeredDofEntity,
	  const_mem_fun<FENumeredDofEntity::interface_type_Field,boost::string_ref,&FENumeredDofEntity::get_name_ref>,
	  const_mem_fun<FENumeredDofEntity::interface_type_RefEntity,EntityType,&FENumeredDofEntity::get_ent_type>,
	  KeyFromKey<
	    member<SideNumber,char,&SideNumber::side_number>,
	    member<FENumeredDofEntity::BaseFEDofEntity,boost::shared_ptr<SideNumber>,&FENumeredDofEntity::sideNumberPtr>
	  >
      > >,
    ordered_non_unique<
      tag<Composite_Name_And_Type_mi_tag>,
      composite_key<
	FENumeredDofEntity,
	  const_mem_fun<FENumeredDofEntity::interface_type_Field,boost::string_ref,&FENumeredDofEntity::get_name_ref>,
	  const_mem_fun<FENumeredDofEntity::interface_type_RefEntity,EntityType,&FENumeredDofEntity::get_ent_type>
	> >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_mi_tag>,
      composite_key<
	FENumeredDofEntity,
	  const_mem_fun<FENumeredDofEntity::interface_type_Field,boost::string_ref,&FENumeredDofEntity::get_name_ref>,
	  const_mem_fun<FENumeredDofEntity::interface_type_DofEntity,EntityHandle,&FENumeredDofEntity::get_ent>
	> >
  > > FENumeredDofEntity_multiIndex;

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps NumeredDofEntity
 * \ingroup dof_multi_indices
 */
typedef multi_index_container<
  boost::shared_ptr<NumeredDofEntity>,
  //unique
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>,
      const_mem_fun<NumeredDofEntity::interface_type_DofEntity,const GlobalUId,&NumeredDofEntity::get_global_unique_id> >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>,
      composite_key<
	      NumeredDofEntity,
	      const_mem_fun<NumeredDofEntity::interface_type_Field,boost::string_ref,&NumeredDofEntity::get_name_ref>,
	      const_mem_fun<NumeredDofEntity::interface_type_DofEntity,EntityHandle,&NumeredDofEntity::get_ent>,
	      const_mem_fun<NumeredDofEntity::interface_type_DofEntity,DofIdx,&NumeredDofEntity::get_EntDofIdx>
    > >,
    //non unique
    ordered_non_unique<
      tag<Idx_mi_tag>, member<NumeredDofEntity,DofIdx,&NumeredDofEntity::dof_idx> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<NumeredDofEntity::interface_type_Field,boost::string_ref,&NumeredDofEntity::get_name_ref> >,
    ordered_non_unique<
      tag<PetscGlobalIdx_mi_tag>, member<NumeredDofEntity,DofIdx,&NumeredDofEntity::petsc_gloabl_dof_idx> >,
    ordered_non_unique<
      tag<PetscLocalIdx_mi_tag>, member<NumeredDofEntity,DofIdx,&NumeredDofEntity::petsc_local_dof_idx> >,
    ordered_non_unique<
      tag<Part_mi_tag>, member<NumeredDofEntity,unsigned int,&NumeredDofEntity::part> >,
    ordered_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<NumeredDofEntity::interface_type_DofEntity,EntityHandle,&NumeredDofEntity::get_ent> >,
    ordered_non_unique<
      tag<Order_mi_tag>, const_mem_fun<NumeredDofEntity::interface_type_DofEntity,ApproximationOrder,&NumeredDofEntity::get_dof_order> >,
    ordered_non_unique<
      tag<Composite_Part_And_Oder_mi_tag>,
      composite_key<
	      NumeredDofEntity,
	       member<NumeredDofEntity,unsigned int,&NumeredDofEntity::part>,
	      const_mem_fun<NumeredDofEntity::interface_type_DofEntity,ApproximationOrder,&NumeredDofEntity::get_dof_order>
	    > >,
    ordered_non_unique<
      tag<Composite_Name_Part_And_CoeffIdx_mi_tag>,
      composite_key<
	     NumeredDofEntity,
	      const_mem_fun<NumeredDofEntity::interface_type_Field,boost::string_ref,&NumeredDofEntity::get_name_ref>,
	      member<NumeredDofEntity,unsigned int,&NumeredDofEntity::part>,
	      const_mem_fun<NumeredDofEntity::interface_type_DofEntity,FieldCoefficientsNumber,&NumeredDofEntity::get_dof_coeff_idx>
	    > >,
    ordered_non_unique<
      tag<Composite_Name_And_Part_mi_tag>,
      composite_key<
	     NumeredDofEntity,
	     const_mem_fun<NumeredDofEntity::interface_type_Field,boost::string_ref,&NumeredDofEntity::get_name_ref>,
	     member<NumeredDofEntity,unsigned int,&NumeredDofEntity::part>
	  > >,
    ordered_non_unique<
      tag<Composite_Name_Ent_And_Part_mi_tag>,
      composite_key<
	     NumeredDofEntity,
	     const_mem_fun<NumeredDofEntity::interface_type_Field,boost::string_ref,&NumeredDofEntity::get_name_ref>,
	     const_mem_fun<NumeredDofEntity::interface_type_DofEntity,EntityHandle,&NumeredDofEntity::get_ent>,
	     member<NumeredDofEntity,unsigned int,&NumeredDofEntity::part>
	  > >,
    ordered_non_unique<
      tag<Composite_Name_And_HasLocalIdx_mi_tag>,
      composite_key<
	     NumeredDofEntity,
	     const_mem_fun<NumeredDofEntity::interface_type_Field,boost::string_ref,&NumeredDofEntity::get_name_ref>,
	     const_mem_fun<NumeredDofEntity,bool,&NumeredDofEntity::get_has_local_index>
	  > >
  > > NumeredDofEntity_multiIndex;

typedef multi_index_container<
  boost::shared_ptr<NumeredDofEntity>,
  indexed_by<
    ordered_unique<
      const_mem_fun<NumeredDofEntity,DofIdx,&NumeredDofEntity::get_dof_idx> >
  > > NumeredDofEntity_multiIndex_uid_view_ordered;

typedef multi_index_container<
  boost::shared_ptr<NumeredDofEntity>,
  indexed_by<
    hashed_unique<
      const_mem_fun<NumeredDofEntity,DofIdx,&NumeredDofEntity::get_dof_idx> >
  > > NumeredDofEntity_multiIndex_uid_view_hashed;

typedef multi_index_container<
  boost::shared_ptr<NumeredDofEntity>,
  indexed_by<
    ordered_non_unique<
      const_mem_fun<NumeredDofEntity,DofIdx,&NumeredDofEntity::get_petsc_local_dof_idx> >
  > > NumeredDofEntity_multiIndex_petsc_local_dof_view_ordered_non_unique;

struct DofEntity_active_change {
  bool active;
  DofEntity_active_change(bool _active);
  void operator()(boost::shared_ptr<DofEntity> &_dof_);
};

struct NumeredDofEntity_part_change {
  unsigned int pArt;
  DofIdx petscGloablDofIdx;
  NumeredDofEntity_part_change(const unsigned int part,const DofIdx petsc_gloabl_dof_idx):
  pArt(part),
  petscGloablDofIdx(petsc_gloabl_dof_idx) {};
  void operator()(boost::shared_ptr<NumeredDofEntity> &dof) {
    dof->part = pArt;
    dof->petsc_gloabl_dof_idx = petscGloablDofIdx;
  }
};

struct NumeredDofEntity_local_idx_change {
  DofIdx petscLocalDofIdx;
  NumeredDofEntity_local_idx_change(const DofIdx petsc_local_dof_idx):
  petscLocalDofIdx(petsc_local_dof_idx) {};
  void operator()(boost::shared_ptr<NumeredDofEntity> &dof) {
    dof->petsc_local_dof_idx = petscLocalDofIdx;
  }
};

struct NumeredDofEntity_mofem_index_change {
  DofIdx mofemIdx;
  NumeredDofEntity_mofem_index_change(const DofIdx mofem_idx):
  mofemIdx(mofem_idx) {};
  void operator()(boost::shared_ptr<NumeredDofEntity> &dof) {
    dof->dof_idx = mofemIdx;
  }
};

typedef multi_index_container<
  boost::shared_ptr<NumeredDofEntity>,
  indexed_by<
    ordered_unique<
      member<NumeredDofEntity,const DofIdx,&NumeredDofEntity::petsc_gloabl_dof_idx> >
 > > NumeredDofEntity_multiIndex_global_index_view;

}
#endif // __DOFSMULTIINDICES_HPP__

/***************************************************************************//**
 * \defgroup dof_multi_indices Dofs structures and multi-indices
 * \ingroup mofem
 ******************************************************************************/
