/** \file FEMMultiIndices.hpp
 * \brief Myltindex contains, data structures for mofem finite elements and other low-level functions
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
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

#ifndef __FEMMULTIINDICES_HPP__
#define __FEMMULTIINDICES_HPP__

namespace MoFEM {

/**
 * \brief keeps data about abstract refined finite element
 * \ingroup fe_multi_indices
 */
struct RefMoFEMElement: public interface_RefMoFEMEntity<RefMoFEMEntity> {
  typedef interface_RefMoFEMEntity<RefMoFEMEntity> interface_type_RefMoFEMEntity;

  static BitRefEdges DummyBitRefEdges;

  SideNumber_multiIndex side_number_table;
  RefMoFEMElement(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr);
  virtual const BitRefEdges& get_BitRefEdges() const { return DummyBitRefEdges; }
  virtual int get_BitRefEdges_ulong() const { return 0; }
  SideNumber_multiIndex &get_side_number_table() const { return const_cast<SideNumber_multiIndex&>(side_number_table); };
  virtual SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const {
    NOT_USED(moab);
    NOT_USED(ent);
    return NULL;
  };
  const RefMoFEMElement* get_RefMoFEMElement() const { return this; }
  friend ostream& operator<<(ostream& os,const RefMoFEMElement& e);
};

/**
 * \brief keeps data about abstract MESHSET finite element
 * \ingroup fe_multi_indices
 */
struct RefMoFEMElement_MESHSET: public RefMoFEMElement {
  RefMoFEMElement_MESHSET(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr);
  const RefMoFEMElement* get_RefMoFEMElement() const { return this; }
  SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const;
};

/**
 * \brief keeps data about abstract PRISM finite element
 * \ingroup fe_multi_indices
 */
struct RefMoFEMElement_PRISM: public RefMoFEMElement {
  BitRefEdges *tag_BitRefEdges;
  RefMoFEMElement_PRISM(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr);
  const RefMoFEMElement* get_RefMoFEMElement() const { return this; }
  SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const;
  const BitRefEdges& get_BitRefEdges() const { return *tag_BitRefEdges; }
  int get_BitRefEdges_ulong() const { return get_BitRefEdges().to_ulong(); }
};

/**
 * \brief keeps data about abstract TET finite element
 * \ingroup fe_multi_indices
 */
struct RefMoFEMElement_TET: public RefMoFEMElement {
  BitRefEdges *tag_BitRefEdges;
  const int* tag_type_data;
  RefMoFEMElement_TET(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr);
  const RefMoFEMElement* get_RefMoFEMElement() const { return this; }
  SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const;
  SideNumber_multiIndex &get_side_number_table() const { return const_cast<SideNumber_multiIndex&>(side_number_table); };
  const BitRefEdges& get_BitRefEdges() const { return *tag_BitRefEdges; }
  int get_BitRefEdges_ulong() const { return get_BitRefEdges().to_ulong(); }
  inline int get_ref_type() const { return tag_type_data[0]; }
  inline int get_ref_sub_type() const { return tag_type_data[1]; }
  friend ostream& operator<<(ostream& os,const RefMoFEMElement_TET& e);
};

/**
 * \brief keeps data about abstract TRI finite element
 * \ingroup fe_multi_indices
 */
struct RefMoFEMElement_TRI: public RefMoFEMElement {
  RefMoFEMElement_TRI(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr);
  const RefMoFEMElement* get_RefMoFEMElement() const { return this; }
  SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const;
  friend ostream& operator<<(ostream& os,const RefMoFEMElement_TRI& e);
};

/**
 * \brief keeps data about abstract EDGE finite element
 * \ingroup fe_multi_indices
 */
struct RefMoFEMElement_EDGE: public RefMoFEMElement {
  RefMoFEMElement_EDGE(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr);
  const RefMoFEMElement* get_RefMoFEMElement() const { return this; }
  SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const;
  friend ostream& operator<<(ostream& os,const RefMoFEMElement_EDGE& e);
};

/**
 * \brief keeps data about abstract VERTEX finite element
 * \ingroup fe_multi_indices
 */
struct RefMoFEMElement_VERTEX: public RefMoFEMElement {
  RefMoFEMElement_VERTEX(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr);
  const RefMoFEMElement* get_RefMoFEMElement() const { return this; }
  SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const;
  friend ostream& operator<<(ostream& os,const RefMoFEMElement_VERTEX& e);
};

/**
 * \brief intrface to RefMoFEMElement
 * \ingroup fe_multi_indices
 */
template<typename T>
struct interface_RefMoFEMElement: interface_RefMoFEMEntity<T> {
  interface_RefMoFEMElement(const T *_ref_ptr): interface_RefMoFEMEntity<T>(_ref_ptr) {}
  int get_BitRefEdges_ulong() const { return interface_RefMoFEMEntity<T>::ref_ptr->get_BitRefEdges_ulong(); }
  SideNumber_multiIndex &get_side_number_table() const { return interface_RefMoFEMEntity<T>::ref_ptr->get_side_number_table(); }
  SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const { return interface_RefMoFEMEntity<T>::ref_ptr->get_side_number_ptr(moab,ent); }
  inline const RefMoFEMElement* get_RefMoFEMElement() const { return interface_RefMoFEMEntity<T>::ref_ptr->get_RefMoFEMElement(); }
  virtual ~interface_RefMoFEMElement() {}
};

struct ptrWrapperRefMoFEMElement: public interface_RefMoFEMElement<RefMoFEMElement> {
  typedef interface_RefMoFEMEntity<RefMoFEMElement> interface_type_RefMoFEMEntity;
  typedef interface_RefMoFEMElement<RefMoFEMElement> interface_type_RefMoFEMElement;
  int wrapp;
  ptrWrapperRefMoFEMElement(const RefMoFEMElement *__ptr): interface_RefMoFEMElement<RefMoFEMElement>(__ptr),wrapp(1) {}
  ptrWrapperRefMoFEMElement(const ptrWrapperRefMoFEMElement &ref): interface_RefMoFEMElement<RefMoFEMElement>(ref) {
    wrapp = 1;
    assert(ref.wrapp == 1);
    (const_cast<ptrWrapperRefMoFEMElement&>(ref)).wrapp++;
  }
  virtual ~ptrWrapperRefMoFEMElement() {
    if(wrapp == 1) {
      delete interface_RefMoFEMEntity<RefMoFEMElement>::ref_ptr;
    }
  }
};

/**
 * \typedef RefMoFEMElement_multiIndex
 * type multiIndex container for RefMoFEMElement
 * \ingroup fe_multi_indices
 *
 * \param hashed_unique Ent_mi_tag
 * \param ordered_non_unique Meshset_mi_tag
 * \param ordered_non_unique Ent_Ent_mi_tag
 * \param ordered_non_unique Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag
 */
typedef multi_index_container<
  ptrWrapperRefMoFEMElement,
  indexed_by<
    hashed_unique<
      tag<Ent_mi_tag>, const_mem_fun<ptrWrapperRefMoFEMElement::interface_type_RefMoFEMEntity,EntityHandle,&ptrWrapperRefMoFEMElement::get_ref_ent> >,
    ordered_non_unique<
      tag<Ent_Ent_mi_tag>, const_mem_fun<ptrWrapperRefMoFEMElement::interface_type_RefMoFEMEntity,EntityHandle,&ptrWrapperRefMoFEMElement::get_parent_ent> >,
    ordered_non_unique<
      tag<EntType_mi_tag>, const_mem_fun<ptrWrapperRefMoFEMElement::interface_type_RefMoFEMEntity,EntityType,&ptrWrapperRefMoFEMElement::get_ent_type> >,
    ordered_non_unique<
      tag<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>,
      composite_key<
	ptrWrapperRefMoFEMElement,
	const_mem_fun<ptrWrapperRefMoFEMElement::interface_type_RefMoFEMEntity,EntityHandle,&ptrWrapperRefMoFEMElement::get_parent_ent>,
	const_mem_fun<ptrWrapperRefMoFEMElement::interface_type_RefMoFEMElement,int,&ptrWrapperRefMoFEMElement::get_BitRefEdges_ulong> > >,
    hashed_unique<
      tag<Composite_EntType_and_ParentEntType_mi_tag>,
      composite_key<
	ptrWrapperRefMoFEMElement,
	const_mem_fun<ptrWrapperRefMoFEMElement::interface_type_RefMoFEMEntity,EntityHandle,&ptrWrapperRefMoFEMElement::get_ref_ent>,
	const_mem_fun<ptrWrapperRefMoFEMElement::interface_type_RefMoFEMEntity,EntityHandle,&ptrWrapperRefMoFEMElement::get_parent_ent> > >
  > > RefMoFEMElement_multiIndex;

/** \brief change parent
  * \ingroup  fe_multi_indices
  *
  * Using this function with care. Some other multi-indices can deponent on this.

  Known dependent multi-indices (verify if that list is full):
  - RefMoFEMEntity_multiIndex
  - RefMoFEMElement_multiIndex

  */
struct RefMoFEMElement_change_parent {
  Interface &mOab;
  const RefMoFEMEntity_multiIndex *refEntPtr;
  RefMoFEMEntity_multiIndex::iterator refEntIt;
  EntityHandle pArent;
  ErrorCode rval;
  RefMoFEMElement_change_parent(Interface &moab,
    const RefMoFEMEntity_multiIndex *ref_ent_ptr,
    RefMoFEMEntity_multiIndex::iterator ref_ent_it,
    EntityHandle parent):
    mOab(moab),
    refEntPtr(ref_ent_ptr),
    refEntIt(ref_ent_it),
    pArent(parent) {}
  void operator()(ptrWrapperRefMoFEMElement &e) {
    const_cast<RefMoFEMEntity_multiIndex*>(refEntPtr)->modify(refEntIt,RefMoFEMEntity_change_parent(mOab,pArent));
  }
};

struct EntMoFEMFiniteElement;

/** \brief user adjacency function table
  * \ingroup fe_multi_indices
  */
typedef PetscErrorCode (*ElementAdjacencyTable[MBMAXTYPE])(
  Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency);

/** \brief user adjacency function
  * \ingroup fe_multi_indices
  */
typedef PetscErrorCode (*ElementAdjacencyFunct)(
  Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency);

/**
 * \brief Finite element definition
 * \ingroup fe_multi_indices
 */
struct MoFEMFiniteElement {
  EntityHandle meshset;     ///< meshset stores FE ents
  BitFEId* tag_id_data;     ///< ptr to tag storing FE id
  void* tag_name_data;      ///< ptr to tag storing FE name
  int tag_name_size;        ///< numer of characters in FE name
  BitFieldId* tag_BitFieldId_col_data;  ///< tag stores col id_id for fields
  BitFieldId* tag_BitFieldId_row_data;  ///< tag stores row id_id for fields
  BitFieldId* tag_BitFieldId_data;      ///< tag stores data id_id for fields
  MoFEMFiniteElement(Interface &moab,const EntityHandle _meshset);
  inline BitFEId get_id() const { return *tag_id_data; };
  /// get meshset
  inline EntityHandle get_meshset() const { return meshset; }
  /// get FE name
  inline boost::string_ref get_name_ref() const { return boost::string_ref((char *)tag_name_data,tag_name_size); }
  inline string get_name() const { return string((char *)tag_name_data,tag_name_size); }
  /// get BitFieldId col
  inline BitFieldId get_BitFieldId_col() const { return *((BitFieldId*)tag_BitFieldId_col_data); }
  /// get BitFieldId row
  inline BitFieldId get_BitFieldId_row() const { return *((BitFieldId*)tag_BitFieldId_row_data); }
  /// get BitFieldId data
  inline BitFieldId get_BitFieldId_data() const { return *((BitFieldId*)tag_BitFieldId_data); }
  /// get bit number
  inline unsigned int get_bit_number() const { return ffsl(((BitFieldId*)tag_id_data)->to_ulong()); }

  ElementAdjacencyTable element_adjacency_table;  //<- allow to add user specific adjacency map

  friend ostream& operator<<(ostream& os, const MoFEMFiniteElement& e);
};

/** \brief default adjacency map
  * \ingroup fe_multi_indices
  */
struct DefaultElementAdjacency {

  static PetscErrorCode defaultVertex(Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency);
  static PetscErrorCode defaultEdge(Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency);
  static PetscErrorCode defaultTri(Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency);
  static PetscErrorCode defaultTet(Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency);
  static PetscErrorCode defaultPrism(Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency);
  static PetscErrorCode defaultMeshset(Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency);

};

/**
 * \brief Inetface for FE
 * \ingroup fe_multi_indices
 */
template <typename T>
struct interface_MoFEMFiniteElement {
  const T *fe_ptr;
  interface_MoFEMFiniteElement(const T *_ptr): fe_ptr(_ptr) {};
  inline BitFEId get_id() const { return fe_ptr->get_id(); }
  inline EntityHandle get_meshset() const { return fe_ptr->get_meshset(); }
  inline boost::string_ref get_name_ref() const { return fe_ptr->get_name_ref(); }
  inline string get_name() const { return fe_ptr->get_name(); }
  inline BitFieldId get_BitFieldId_col() const { return fe_ptr->get_BitFieldId_col(); }
  inline BitFieldId get_BitFieldId_row() const { return fe_ptr->get_BitFieldId_row(); }
  inline BitFieldId get_BitFieldId_data() const { return fe_ptr->get_BitFieldId_data(); }
  inline unsigned int get_bit_number() const { return fe_ptr->get_bit_number(); }
};

/**
 * \brief Finite element data for entitiy
 * \ingroup fe_multi_indices
 */
struct EntMoFEMFiniteElement: public interface_MoFEMFiniteElement<MoFEMFiniteElement>,interface_RefMoFEMElement<RefMoFEMElement> {
  typedef interface_RefMoFEMEntity<RefMoFEMElement> interface_type_RefMoFEMEntity;
  typedef interface_RefMoFEMElement<RefMoFEMElement> interface_type_RefMoFEMElement;
  typedef interface_MoFEMFiniteElement<MoFEMFiniteElement> interface_type_MoFEMFiniteElement;
  DofMoFEMEntity_multiIndex_uid_view row_dof_view;
  DofMoFEMEntity_multiIndex_uid_view col_dof_view;
  DofMoFEMEntity_multiIndex_uid_view data_dof_view;
  FEDofMoFEMEntity_multiIndex data_dofs;
  GlobalUId global_uid;
  EntMoFEMFiniteElement(Interface &moab,const RefMoFEMElement *_ref_MoFEMFiniteElement,const MoFEMFiniteElement *_MoFEMFiniteElement_ptr);
  inline const MoFEMFiniteElement* get_MoFEMFiniteElementPtr() { return interface_MoFEMFiniteElement<MoFEMFiniteElement>::fe_ptr; };
  const GlobalUId& get_global_unique_id() const { return global_uid; }
  GlobalUId get_global_unique_id_calculate() const {
    char bit_number = get_bit_number();
    assert(bit_number<=32);
    GlobalUId _uid_ = (ref_ptr->get_ref_ent())|(((GlobalUId)bit_number)<<(8*sizeof(EntityHandle)));
    return _uid_;
  }
  inline EntityHandle get_ent() const { return get_ref_ent(); }
  inline DofIdx get_nb_dofs_row() const { return row_dof_view.size(); }
  inline DofIdx get_nb_dofs_col() const { return col_dof_view.size(); }
  inline DofIdx get_nb_dofs_data() const { return data_dof_view.size(); }
  inline const FEDofMoFEMEntity_multiIndex& get_data_dofs() const { return data_dofs; };
  friend ostream& operator<<(ostream& os,const EntMoFEMFiniteElement& e);
  PetscErrorCode get_MoFEMFiniteElement_row_dof_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_active_view &dofs_view,
    const int operation_type = Interface::UNION) const;
  PetscErrorCode get_MoFEMFiniteElement_col_dof_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_active_view &dofs_view,
    const int operation_type = Interface::UNION) const;
  PetscErrorCode get_MoFEMFiniteElement_data_dof_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_active_view &dofs_view,
    const int operation_type = Interface::UNION) const;
  PetscErrorCode get_MoFEMFiniteElement_row_dof_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_uid_view &dofs_view,
    const int operation_type = Interface::UNION) const;
  PetscErrorCode get_MoFEMFiniteElement_col_dof_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_uid_view &dofs_view,
    const int operation_type = Interface::UNION) const;
  PetscErrorCode get_MoFEMFiniteElement_row_dof_view(
    const NumeredDofMoFEMEntity_multiIndex &dofs,NumeredDofMoFEMEntity_multiIndex_uid_view_ordered &dofs_view,
    const int operation_type = Interface::UNION) const;
  PetscErrorCode get_MoFEMFiniteElement_col_dof_view(
    const NumeredDofMoFEMEntity_multiIndex &dofs,NumeredDofMoFEMEntity_multiIndex_uid_view_ordered &dofs_view,
    const int operation_type = Interface::UNION) const;
  PetscErrorCode get_MoFEMFiniteElement_row_dof_view(
    const NumeredDofMoFEMEntity_multiIndex &dofs,NumeredDofMoFEMEntity_multiIndex_uid_view_hashed &dofs_view,
    const int operation_type = Interface::UNION) const;
  PetscErrorCode get_MoFEMFiniteElement_col_dof_view(
    const NumeredDofMoFEMEntity_multiIndex &dofs,NumeredDofMoFEMEntity_multiIndex_uid_view_hashed &dofs_view,
    const int operation_type = Interface::UNION) const;

  PetscErrorCode get_element_adjacency(Interface &moab,const MoFEMField *field_ptr,Range &adjacency) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    const EntMoFEMFiniteElement *this_fe_ptr = this;
    if(get_MoFEMFiniteElementPtr()->element_adjacency_table[get_ent_type()] == NULL) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }
    ierr = (get_MoFEMFiniteElementPtr()->element_adjacency_table[get_ent_type()])(
      moab,field_ptr,this_fe_ptr,adjacency
    ); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};

/**
 * \brief interface to EntMoFEMFiniteElement
 * \ingroup fe_multi_indices
 */
template <typename T>
struct interface_EntMoFEMFiniteElement:public interface_MoFEMFiniteElement<T>,interface_RefMoFEMElement<T> {
  interface_EntMoFEMFiniteElement(const T *_ptr): interface_MoFEMFiniteElement<T>(_ptr),interface_RefMoFEMElement<T>(_ptr) {};
  inline const MoFEMFiniteElement* get_MoFEMFiniteElementPtr() { return interface_MoFEMFiniteElement<T>::get_MoFEMFiniteElementPtr(); };
  inline EntityID get_ent_id() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_ent_id(); }
  inline EntityType get_ent_type() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_ent_type(); }
  //
  inline const FEDofMoFEMEntity_multiIndex& get_data_dofs() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_data_dofs(); };
  inline DofIdx get_nb_dofs_row() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_nb_dofs_row(); }
  inline DofIdx get_nb_dofs_col() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_nb_dofs_col(); }
  inline DofIdx get_nb_dofs_data() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_nb_dofs_data(); }
  inline EntityHandle get_ent() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_ref_ent(); };
  inline GlobalUId get_global_unique_id() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_global_unique_id(); }
  //
  SideNumber_multiIndex &get_side_number_table() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_side_number_table(); }
  SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_side_number_ptr(moab,ent); }
  //
  inline PetscErrorCode get_element_adjacency(Interface &moab,const MoFEMField *field_ptr,Range &adjacency) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = interface_MoFEMFiniteElement<T>::get_element_adjacency(moab,field_ptr,adjacency); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
};

/** \brief Partitioned Finite Element in Problem
 * \ingroup fe_multi_indices
 */
struct NumeredMoFEMFiniteElement: public interface_EntMoFEMFiniteElement<EntMoFEMFiniteElement> {
  typedef interface_MoFEMFiniteElement<EntMoFEMFiniteElement> interface_type_MoFEMFiniteElement;
  typedef interface_EntMoFEMFiniteElement<EntMoFEMFiniteElement> interface_type_EntMoFEMFiniteElement;
  unsigned int part;
  FENumeredDofMoFEMEntity_multiIndex rows_dofs;
  FENumeredDofMoFEMEntity_multiIndex cols_dofs;
  NumeredMoFEMFiniteElement(const EntMoFEMFiniteElement *EntMoFEMFiniteElement_ptr): interface_EntMoFEMFiniteElement<EntMoFEMFiniteElement>(EntMoFEMFiniteElement_ptr),part(-1) {};
  inline unsigned int get_part() const { return part; };

  /** \brief get FE dof
    * \ingroup mofem_dofs
    */
  inline const FENumeredDofMoFEMEntity_multiIndex& get_rows_dofs() const { return rows_dofs; };

  /** \brief get FE dof
    * \ingroup mofem_dofs
    */
  inline const FENumeredDofMoFEMEntity_multiIndex& get_cols_dofs() const { return cols_dofs; };

  /** \brief get FE dof by petsc index
    * \ingroup mofem_dofs
    */
  PetscErrorCode get_row_dofs_by_petsc_gloabl_dof_idx(DofIdx idx,const FENumeredDofMoFEMEntity **dof_ptr) const;

  /** \brief get FE dof by petsc index
    * \ingroup mofem_dofs
    */
  PetscErrorCode get_col_dofs_by_petsc_gloabl_dof_idx(DofIdx idx,const FENumeredDofMoFEMEntity **dof_ptr) const;

  friend ostream& operator<<(ostream& os,const NumeredMoFEMFiniteElement& e) {
    os << "part " << e.part << " " << *(e.fe_ptr);
    return os;
  }

  /**
   * Loop over DOFs in row on element
   * @param  FEPTR pointer to element structure \ref NumeredMoFEMFiniteElement
   * @param  IT    iterator
   * @return       user return in for(_IT_FENUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(FEPTR,IT))
   * \ingroup fe_multi_indices
   */
  #define _IT_FENUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(FEPTR,IT) \
  FENumeredDofMoFEMEntity_multiIndex::iterator IT = FEPTR->rows_dofs.begin(); \
  IT!=FEPTR->rows_dofs.end(); IT++

  /**
   * Loop over DOFs in col on element
   * @param  FEPTR pointer to element structure \ref NumeredMoFEMFiniteElement
   * @param  IT    iterator
   * @return       user return in for(_IT_FENUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(FEPTR,IT))
   * \ingroup fe_multi_indices
   */
  #define _IT_FENUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(FEPTR,IT) \
  FENumeredDofMoFEMEntity_multiIndex::iterator IT = FEPTR->cols_dofs.begin(); \
  IT!=FEPTR->cols_dofs.end(); IT++

  /**
   * Loop over DOFs in row on element for particular filed
   * @param  FEPTR pointer to element structure \ref NumeredMoFEMFiniteElement
   * @param  NAME  name of filed
   * @param  IT    iterator
   * @return       user return in for(_IT_FENUMEREDDOFMOFEMENTITY_BY_NAME_ROW_FOR_LOOP_(FEPTR,NAME,IT))
   * \ingroup fe_multi_indices
   */
  #define _IT_FENUMEREDDOFMOFEMENTITY_BY_NAME_ROW_FOR_LOOP_(FEPTR,NAME,IT) \
  FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator IT = FEPTR->rows_dofs.get<FieldName_mi_tag>().lower_bound(NAME); \
  IT!=FEPTR->rows_dofs.get<FieldName_mi_tag>().upper_bound(NAME); IT++

  /**
   * Loop over DOFs in col on element for particular filed
   * @param  FEPTR pointer to element structure \ref NumeredMoFEMFiniteElement
   * @param  NAME  name of filed
   * @param  IT    iterator
   * @return       user return in for(_IT_FENUMEREDDOFMOFEMENTITY_BY_NAME_COL_FOR_LOOP_(FEPTR,NAME,IT))
   * \ingroup fe_multi_indices
   */
  #define _IT_FENUMEREDDOFMOFEMENTITY_BY_NAME_COL_FOR_LOOP_(FEPTR,NAME,IT) \
  FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator IT = FEPTR->cols_dofs.get<FieldName_mi_tag>().lower_bound(NAME); \
  IT!=FEPTR->cols_dofs.get<FieldName_mi_tag>().upper_bound(NAME); IT++



};

/** \brief interface for NumeredMoFEMFiniteElement
 * \ingroup fe_multi_indices
 */
template <typename T>
struct interface_NumeredMoFEMFiniteElement: public interface_EntMoFEMFiniteElement<T> {
  const T *ptr;
  interface_NumeredMoFEMFiniteElement(const T *_ptr): interface_EntMoFEMFiniteElement<T>(_ptr) {};
  inline unsigned int get_part() const { return ptr->get_part(); }
  inline const FENumeredDofMoFEMEntity_multiIndex& get_rows_dofs() const { return ptr->get_rows_dofs(); };
  inline const FENumeredDofMoFEMEntity_multiIndex& get_cols_dofs() const { return ptr->get_cols_dofs(); };
};

/**
 * @relates multi_index_container
 * \brief MultiIndex container for EntMoFEMFiniteElement
 * \ingroup fe_multi_indices
 *
 */
typedef multi_index_container<
  EntMoFEMFiniteElement,
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>, member<EntMoFEMFiniteElement,GlobalUId,&EntMoFEMFiniteElement::global_uid> >,
    ordered_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<EntMoFEMFiniteElement,EntityHandle,&EntMoFEMFiniteElement::get_ent> >,
    ordered_non_unique<
      tag<FiniteElement_name_mi_tag>, const_mem_fun<EntMoFEMFiniteElement::interface_type_MoFEMFiniteElement,boost::string_ref,&EntMoFEMFiniteElement::get_name_ref> >,
    ordered_non_unique<
      tag<BitFEId_mi_tag>, const_mem_fun<EntMoFEMFiniteElement::interface_type_MoFEMFiniteElement,BitFEId,&EntMoFEMFiniteElement::get_id>, LtBit<BitFEId> >,
    ordered_non_unique<
      tag<EntType_mi_tag>, const_mem_fun<EntMoFEMFiniteElement::interface_type_RefMoFEMEntity,EntityType,&EntMoFEMFiniteElement::get_ent_type> >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_mi_tag>,
      composite_key<
	EntMoFEMFiniteElement,
	const_mem_fun<EntMoFEMFiniteElement::interface_type_MoFEMFiniteElement,boost::string_ref,&EntMoFEMFiniteElement::get_name_ref>,
	const_mem_fun<EntMoFEMFiniteElement,EntityHandle,&EntMoFEMFiniteElement::get_ent> > >
  > > EntMoFEMFiniteElement_multiIndex;

/**
  @relates multi_index_container
  \brief MultiIndex for entities for NumeredMoFEMFiniteElement
  \ingroup fe_multi_indices
 */
typedef multi_index_container<
  NumeredMoFEMFiniteElement,
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>, const_mem_fun<NumeredMoFEMFiniteElement::interface_type_EntMoFEMFiniteElement,GlobalUId,&NumeredMoFEMFiniteElement::get_global_unique_id> >,
    ordered_non_unique<
      tag<FiniteElement_name_mi_tag>, const_mem_fun<NumeredMoFEMFiniteElement::interface_type_MoFEMFiniteElement,boost::string_ref,&NumeredMoFEMFiniteElement::get_name_ref> >,
    ordered_non_unique<
      tag<FiniteElement_Part_mi_tag>, member<NumeredMoFEMFiniteElement,unsigned int,&NumeredMoFEMFiniteElement::part> >,
    ordered_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<NumeredMoFEMFiniteElement::interface_type_EntMoFEMFiniteElement,EntityHandle,&NumeredMoFEMFiniteElement::get_ent> >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_mi_tag>,
      composite_key<
      NumeredMoFEMFiniteElement,
      const_mem_fun<NumeredMoFEMFiniteElement::interface_type_MoFEMFiniteElement,boost::string_ref,&NumeredMoFEMFiniteElement::get_name_ref>,
      const_mem_fun<NumeredMoFEMFiniteElement::interface_type_EntMoFEMFiniteElement,EntityHandle,&NumeredMoFEMFiniteElement::get_ent> > >,
    ordered_non_unique<
      tag<Composite_Name_And_Part_mi_tag>,
      composite_key<
      NumeredMoFEMFiniteElement,
      const_mem_fun<NumeredMoFEMFiniteElement::interface_type_MoFEMFiniteElement,boost::string_ref,&NumeredMoFEMFiniteElement::get_name_ref>,
      member<NumeredMoFEMFiniteElement,unsigned int,&NumeredMoFEMFiniteElement::part> > >
  > > NumeredMoFEMFiniteElement_multiIndex;

/**
  @relates multi_index_container
  \brief MultiIndex for entities for MoFEMFiniteElement
  \ingroup fe_multi_indices
 */
typedef multi_index_container<
  MoFEMFiniteElement,
  indexed_by<
    hashed_unique<
      tag<FiniteElement_Meshset_mi_tag>, member<MoFEMFiniteElement,EntityHandle,&MoFEMFiniteElement::meshset> >,
    hashed_unique<
      tag<BitFEId_mi_tag>, const_mem_fun<MoFEMFiniteElement,BitFEId,&MoFEMFiniteElement::get_id>, HashBit<BitFEId>, EqBit<BitFEId> >,
    ordered_unique<
      tag<FiniteElement_name_mi_tag>, const_mem_fun<MoFEMFiniteElement,boost::string_ref,&MoFEMFiniteElement::get_name_ref> >
  > > MoFEMFiniteElement_multiIndex;

// modificators

struct NumeredMoFEMFiniteElement_change_part {
  unsigned int part;
  NumeredMoFEMFiniteElement_change_part(unsigned int _part): part(_part) {};
  void operator()(NumeredMoFEMFiniteElement &MoFEMFiniteElement) {
    MoFEMFiniteElement.part = part;
  }
};

struct MoFEMFiniteElement_col_change_bit_add {
  BitFieldId f_id_col;
  MoFEMFiniteElement_col_change_bit_add(const BitFieldId _f_id_col): f_id_col(_f_id_col) {};
  void operator()(MoFEMFiniteElement &MoFEMFiniteElement);
};
struct MoFEMFiniteElement_row_change_bit_add {
  BitFieldId f_id_row;
  MoFEMFiniteElement_row_change_bit_add(const BitFieldId _f_id_row): f_id_row(_f_id_row) {};
  void operator()(MoFEMFiniteElement &MoFEMFiniteElement);
};
struct EntMoFEMFiniteElement_change_bit_add {
  BitFieldId f_id_data;
  EntMoFEMFiniteElement_change_bit_add(const BitFieldId _f_id_data): f_id_data(_f_id_data) {};
  void operator()(MoFEMFiniteElement &MoFEMFiniteElement);
};

struct MoFEMFiniteElement_col_change_bit_off {
  BitFieldId f_id_col;
  MoFEMFiniteElement_col_change_bit_off(const BitFieldId _f_id_col): f_id_col(_f_id_col) {};
  void operator()(MoFEMFiniteElement &MoFEMFiniteElement);
};
struct MoFEMFiniteElement_row_change_bit_off {
  BitFieldId f_id_row;
  MoFEMFiniteElement_row_change_bit_off(const BitFieldId _f_id_row): f_id_row(_f_id_row) {};
  void operator()(MoFEMFiniteElement &MoFEMFiniteElement);
};
struct EntMoFEMFiniteElement_change_bit_off {
  BitFieldId f_id_data;
  EntMoFEMFiniteElement_change_bit_off(const BitFieldId _f_id_data): f_id_data(_f_id_data) {};
  void operator()(MoFEMFiniteElement &MoFEMFiniteElement);
};

}

#endif // __FEMMULTIINDICES_HPP__

/***************************************************************************//**
 * \defgroup fe_multi_indices Finite elements structures and multi-indices
 * \ingroup mofem
 ******************************************************************************/
