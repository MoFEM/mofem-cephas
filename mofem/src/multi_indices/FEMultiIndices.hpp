/** \file FEMultiIndices.hpp
 * \brief Myltiindex contains, data structures for mofem finite elements and other low-level functions
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
struct RefElement: public interface_RefEntity<RefEntity> {

  typedef interface_RefEntity<RefEntity> interface_type_RefEntity;

  static BitRefEdges DummyBitRefEdges;

  SideNumber_multiIndex side_number_table;
  RefElement(const boost::shared_ptr<RefEntity> ref_ent_ptr);
  virtual const BitRefEdges& getBitRefEdges() const { return DummyBitRefEdges; }

  virtual int getBitRefEdgesUlong() const { return 0; }

  SideNumber_multiIndex &getSideNumberTable() const {
    return const_cast<SideNumber_multiIndex&>(side_number_table);
  }

  /** \deprecated Use getSideNumberTable() instead
  */
  DEPRECATED SideNumber_multiIndex &get_side_number_table() const {
    return getSideNumberTable();
  }

  virtual boost::shared_ptr<SideNumber> getSideNumberPtr(const EntityHandle ent) const {
    NOT_USED(ent);
    return boost::shared_ptr<SideNumber>();
  };

  /**
   * \deprecated First argument is no longer needed
   */
  virtual DEPRECATED boost::shared_ptr<SideNumber> getSideNumberPtr(
    const moab::Interface &moab,const EntityHandle ent
  ) const {
    NOT_USED(moab);
    return getSideNumberPtr(ent);
  }

  /**
   * \brief Get pointer to RefEntity
   */
  inline boost::shared_ptr<RefEntity>& getRefEntityPtr() const { return this->sPtr; }

  friend std::ostream& operator<<(std::ostream& os,const RefElement& e);

};

/**
 * \brief keeps data about abstract MESHSET finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_MESHSET: public RefElement {
  RefElement_MESHSET(const boost::shared_ptr<RefEntity> ref_ent_ptr);
  boost::shared_ptr<SideNumber> getSideNumberPtr(const EntityHandle ent) const;
};
/**
 * \brief keeps data about abstract PRISM finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_PRISM: public RefElement {
  BitRefEdges *tag_BitRefEdges;
  RefElement_PRISM(const boost::shared_ptr<RefEntity> ref_ent_ptr);
  boost::shared_ptr<SideNumber> getSideNumberPtr(const EntityHandle ent) const;
  const BitRefEdges& getBitRefEdges() const { return *tag_BitRefEdges; }
  int getBitRefEdgesUlong() const { return getBitRefEdges().to_ulong(); }
};

/**
 * \brief keeps data about abstract TET finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_TET: public RefElement {
  BitRefEdges *tag_BitRefEdges;
  const int* tag_type_data;
  RefElement_TET(const boost::shared_ptr<RefEntity> ref_ent_ptr);
  boost::shared_ptr<SideNumber> getSideNumberPtr(const EntityHandle ent) const;
  SideNumber_multiIndex &getSideNumberTable() const { return const_cast<SideNumber_multiIndex&>(side_number_table); };
  const BitRefEdges& getBitRefEdges() const { return *tag_BitRefEdges; }
  int getBitRefEdgesUlong() const { return getBitRefEdges().to_ulong(); }
  inline int getRefType() const { return tag_type_data[0]; }
};

/**
 * \brief keeps data about abstract TRI finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_TRI: public RefElement {
  RefElement_TRI(const boost::shared_ptr<RefEntity> ref_ent_ptr);
  boost::shared_ptr<SideNumber> getSideNumberPtr(const EntityHandle ent) const;
  friend std::ostream& operator<<(std::ostream& os,const RefElement_TRI& e);
};

/**
 * \brief keeps data about abstract EDGE finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_EDGE: public RefElement {
  RefElement_EDGE(const boost::shared_ptr<RefEntity> ref_ent_ptr);
  boost::shared_ptr<SideNumber> getSideNumberPtr(const EntityHandle ent) const;
  friend std::ostream& operator<<(std::ostream& os,const RefElement_EDGE& e);
};

/**
 * \brief keeps data about abstract VERTEX finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_VERTEX: public RefElement {
  RefElement_VERTEX(const boost::shared_ptr<RefEntity> ref_ent_ptr);
  boost::shared_ptr<SideNumber> getSideNumberPtr(const EntityHandle ent) const;
  friend std::ostream& operator<<(std::ostream& os,const RefElement_VERTEX& e);
};

/**
 * \brief intrface to RefElement
 * \ingroup fe_multi_indices
 */
template<typename T>
struct interface_RefElement: interface_RefEntity<T> {

  typedef interface_RefEntity<T> interface_type_RefEntity;
  typedef interface_RefElement<T> interface_type_RefElement;

  const boost::shared_ptr<T> sPtr;

  interface_RefElement(const boost::shared_ptr<T> sptr):
  interface_RefEntity<T>(sptr),
  sPtr(sptr) {}

  inline int getBitRefEdgesUlong() const
  { return this->sPtr->getBitRefEdgesUlong(); }

  inline SideNumber_multiIndex &getSideNumberTable() const
  { return this->sPtr->getSideNumberTable(); }

  DEPRECATED inline SideNumber_multiIndex &get_side_number_table() const
  { return this->sPtr->getSideNumberTable(); }

  inline boost::shared_ptr<SideNumber> getSideNumberPtr(const EntityHandle ent) const
  { return this->sPtr->getSideNumberPtr(ent); }

  /**
   * \deprecated First argument is no longer needed
   */
  virtual DEPRECATED boost::shared_ptr<SideNumber> getSideNumberPtr(
    const moab::Interface &moab,const EntityHandle ent
  ) const {
    NOT_USED(moab);
    return getSideNumberPtr(ent);
  }

  inline boost::shared_ptr<RefEntity>& getRefEntityPtr() const
  { return this->sPtr->getRefEntityPtr(); }

  inline const boost::shared_ptr<T> getRefElement() const { return this->sPtr; }

  virtual ~interface_RefElement() {}

};

typedef interface_RefElement<RefElement> ptrWrapperRefElement;

/**
 * \typedef RefElement_multiIndex
 * type multiIndex container for RefElement
 * \ingroup fe_multi_indices
 *
 * \param hashed_unique Ent_mi_tag
 * \param ordered_non_unique Meshset_mi_tag
 * \param ordered_non_unique Ent_Ent_mi_tag
 * \param ordered_non_unique Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag
 */
typedef multi_index_container<
  ptrWrapperRefElement,
  indexed_by<
    hashed_unique<
      tag<Ent_mi_tag>,
      const_mem_fun<
        ptrWrapperRefElement::interface_type_RefEntity,
        EntityHandle,
        &ptrWrapperRefElement::getRefEnt
      >
    >,
    ordered_non_unique<
      tag<Ent_Ent_mi_tag>,
        const_mem_fun<ptrWrapperRefElement::interface_type_RefEntity,
        EntityHandle,
        &ptrWrapperRefElement::getParentEnt
      >
    >,
    ordered_non_unique<
      tag<EntType_mi_tag>,
      const_mem_fun<
        ptrWrapperRefElement::interface_type_RefEntity,
        EntityType,
        &ptrWrapperRefElement::getEntType
      >
    >,
    ordered_non_unique<
      tag<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>,
      composite_key<
	      ptrWrapperRefElement,
	      const_mem_fun<
          ptrWrapperRefElement::interface_type_RefEntity,
          EntityHandle,
          &ptrWrapperRefElement::getParentEnt
        >,
	      const_mem_fun<
          ptrWrapperRefElement::interface_type_RefElement,
          int,
          &ptrWrapperRefElement::getBitRefEdgesUlong
        >
      >
    >,
    hashed_unique<
      tag<Composite_EntType_and_ParentEntType_mi_tag>,
      composite_key<
	     ptrWrapperRefElement,
	     const_mem_fun<
        ptrWrapperRefElement::interface_type_RefEntity,
        EntityHandle,
        &ptrWrapperRefElement::getRefEnt
       >,
	     const_mem_fun<
        ptrWrapperRefElement::interface_type_RefEntity,
        EntityHandle,
        &ptrWrapperRefElement::getParentEnt
       >
      >
    >
  >
> RefElement_multiIndex;

/** \brief change parent
  * \ingroup  fe_multi_indices
  *
  * Using this function with care. Some other multi-indices can deponent on this.

  Known dependent multi-indices (verify if that list is full):
  - RefEntity_multiIndex
  - RefElement_multiIndex

  */
struct RefElement_change_parent {
  const RefEntity_multiIndex *refEntPtr;
  RefEntity_multiIndex::iterator refEntIt;
  EntityHandle pArent;
  ErrorCode rval;
  RefElement_change_parent(
    const RefEntity_multiIndex *ref_ent_ptr,
    RefEntity_multiIndex::iterator ref_ent_it,
    EntityHandle parent
  ):
  refEntPtr(ref_ent_ptr),
  refEntIt(ref_ent_it),
  pArent(parent) {}
  void operator()(ptrWrapperRefElement &e) {
    const_cast<RefEntity_multiIndex*>(refEntPtr)->modify(refEntIt,RefEntity_change_parent(pArent));
  }
};

struct EntFiniteElement;

/** \brief user adjacency function table
  * \ingroup fe_multi_indices
  */
typedef PetscErrorCode (*ElementAdjacencyTable[MBMAXTYPE])(
  Interface &moab,const Field &field_ptr,const EntFiniteElement &fe_ptr,Range &adjacency);

/** \brief user adjacency function
  * \ingroup fe_multi_indices
  */
typedef PetscErrorCode (*ElementAdjacencyFunct)(
  Interface &moab,const Field &field_ptr,const EntFiniteElement &fe_ptr,Range &adjacency);

/**
 * \brief Finite element definition
 * \ingroup fe_multi_indices
 */
struct FiniteElement {
  EntityHandle meshset;     ///< meshset stores FE ents
  BitFEId* tag_id_data;     ///< ptr to tag storing FE id
  void* tag_name_data;      ///< ptr to tag storing FE name
  int tag_name_size;        ///< numer of characters in FE name
  BitFieldId* tag_BitFieldId_col_data;  ///< tag stores col id_id for fields
  BitFieldId* tag_BitFieldId_row_data;  ///< tag stores row id_id for fields
  BitFieldId* tag_BitFieldId_data;      ///< tag stores data id_id for fields
  FiniteElement(Interface &moab,const EntityHandle _meshset);

  /**
   * \brief Get finite element id
   * @return Finite element Id
   */
  inline BitFEId getId() const { return *tag_id_data; };

  /**
   * \brief Get meshset containing element entities
   * @return Meshset
   */
  inline EntityHandle getMeshset() const { return meshset; }

  /**
   * \brief Get finite element name
   * @return string_ref
   */
  inline boost::string_ref getNameRef() const { return boost::string_ref((char *)tag_name_data,tag_name_size); }

  /**
   * \brief Get finite element name
   * @return string
   */
  inline std::string getName() const { return std::string((char *)tag_name_data,tag_name_size); }

  /**
   * \brief Get field ids on columns
   * @return Bit field ids
   */
  inline BitFieldId getBitFieldIdCol() const { return *((BitFieldId*)tag_BitFieldId_col_data); }

  /**
   * \brief Get field ids on rows
   * @return Bit field ids
   */
  inline BitFieldId getBitFieldIdRow() const { return *((BitFieldId*)tag_BitFieldId_row_data); }

  /**
   * \brief Get field ids on data
   * @return Bit field ids
   */
  inline BitFieldId getBitFieldIdData() const { return *((BitFieldId*)tag_BitFieldId_data); }

  /**
   * \brief Get bit identifying this element
   *
   * Each element like field is identified by bit set. Each element has unique bit set,
   * this function returns number of that bit.
   *
   * @return Bit number
   */
  inline unsigned int getBitNumber() const { return ffsl(((BitFieldId*)tag_id_data)->to_ulong()); }

  ElementAdjacencyTable element_adjacency_table;  //<- allow to add user specific adjacency map

  friend std::ostream& operator<<(std::ostream& os, const FiniteElement& e);

};

/** \brief default adjacency map
  * \ingroup fe_multi_indices
  */
struct DefaultElementAdjacency {

  static PetscErrorCode defaultVertex(
    Interface &moab,const Field& field_ptr,const EntFiniteElement& fe_ptr,Range &adjacency
  );
  static PetscErrorCode defaultEdge(
    Interface &moab,const Field& field_ptr,const EntFiniteElement& fe_ptr,Range &adjacency
  );
  static PetscErrorCode defaultTri(
    Interface &moab,const Field& field_ptr,const EntFiniteElement& fe_ptr,Range &adjacency
  );
  static PetscErrorCode defaultTet(
    Interface &moab,const Field& field_ptr,const EntFiniteElement& fe_ptr,Range &adjacency
  );
  static PetscErrorCode defaultPrism(
    Interface &moab,const Field& field_ptr,const EntFiniteElement& fe_ptr,Range &adjacency
  );
  static PetscErrorCode defaultMeshset(
    Interface &moab,const Field& field_ptr,const EntFiniteElement& fe_ptr,Range &adjacency
  );

};

/**
 * \brief Inetface for FE
 * \ingroup fe_multi_indices
 */
template <typename T>
struct interface_FiniteElement {

  const boost::shared_ptr<T> sFePtr;
  interface_FiniteElement(const boost::shared_ptr<T> ptr): sFePtr(ptr) {};

  inline const boost::shared_ptr<FiniteElement> get_MoFEMFiniteElementPtr() { return this->sFePtr; };

  /**
   * \brief Get finite element id
   * @return Finite element Id
   */
  inline BitFEId getId() const { return this->sFePtr->getId(); }

  /**
   * \brief Get meshset containing element entities
   * @return Meshset
   */
  inline EntityHandle getMeshset() const { return this->sFePtr->getMeshset(); }

  /**
   * \brief Get finite element name
   * @return string_ref
   */
  inline boost::string_ref getNameRef() const { return this->sFePtr->getNameRef(); }

  /**
   * \brief Get finite element name
   * @return string_ref
   */
  inline std::string getName() const { return this->sFePtr->getName(); }

  /**
   * \brief Get field ids on columns
   * @return Bit field ids
   */
  inline BitFieldId getBitFieldIdCol() const { return this->sFePtr->getBitFieldIdCol(); }

  /**
   * \brief Get field ids on rows
   * @return Bit field ids
   */
  inline BitFieldId getBitFieldIdRow() const { return this->sFePtr->getBitFieldIdRow(); }

  /**
   * \brief Get field ids on data
   * @return Bit field ids
   */
  inline BitFieldId getBitFieldIdData() const { return this->sFePtr->getBitFieldIdData(); }

  /**
   * \brief Get bit identifying this element
   *
   * Each element like field is identified by bit set. Each element has unique bit set,
   * this function returns number of that bit.
   *
   * @return Bit number
   */
  inline unsigned int getBitNumber() const { return this->sFePtr->getBitNumber(); }

};

/**
 * \brief Finite element data for entitiy
 * \ingroup fe_multi_indices
 */
struct EntFiniteElement:
public
interface_FiniteElement<FiniteElement>,
interface_RefElement<RefElement> {

  typedef interface_RefEntity<RefElement> interface_type_RefEntity;
  typedef interface_RefElement<RefElement> interface_type_RefElement;
  typedef interface_FiniteElement<FiniteElement> interface_type_MoFEMFiniteElement;
  boost::shared_ptr<DofEntity_multiIndex_uid_view> row_dof_view;
  boost::shared_ptr<DofEntity_multiIndex_uid_view> col_dof_view;
  FEDofEntity_multiIndex data_dofs;
  GlobalUId global_uid;

  EntFiniteElement(
    Interface &moab,
    const boost::shared_ptr<RefElement> ref_finite_element,
    const boost::shared_ptr<FiniteElement> fe_ptr
  );

  /**
   * \brief Get unique UId for finite element entity
   * @return UId
   */
  const GlobalUId& getGlobalUniqueId() const { return global_uid; }

  /**
   * \brief Generaye UId for finite element entity
   * @return [description]
   */
  GlobalUId getGlobalUniqueIdCalculate() const {
    char bit_number = getBitNumber();
    assert(bit_number<=32);
    GlobalUId _uid_ = (sPtr->getRefEnt())|(((GlobalUId)bit_number)<<(8*sizeof(EntityHandle)));
    return _uid_;
  }

  /**
   * \brief Get element entity
   * @return Element entity handle
   */
  inline EntityHandle getEnt() const { return getRefEnt(); }

  /** \deprecated Use getEnt() instead
  */
  DEPRECATED inline EntityHandle get_ent() const { return getEnt(); }

  /**
   * \brief Get number of DOFs on row
   * @return Number of dofs on row
   */
  inline DofIdx getNbDofsRow() const { return row_dof_view->size(); }

  /**
   * \brief Get number of DOFs on col
   * @return Number of dofs on col
   */
  inline DofIdx getNbDofsCol() const { return col_dof_view->size(); }

  /**
   * \brief Get number of DOFs on data
   * @return Number of dofs on data
   */
  inline DofIdx getNbDofsData() const { return data_dofs.size(); }

  /**
   * \brief Get data data dos multi-index structure
   * @return Reference multi-index FEDofEntity_multiIndex
   */
  inline const FEDofEntity_multiIndex& getDataDofs() const { return data_dofs; };

  friend std::ostream& operator<<(std::ostream& os,const EntFiniteElement& e);

  PetscErrorCode getRowDofView(
    const DofEntity_multiIndex &dofs,DofEntity_multiIndex_active_view &dofs_view,
    const int operation_type = moab::Interface::UNION
  ) const;

  PetscErrorCode getColDofView(
    const DofEntity_multiIndex &dofs,DofEntity_multiIndex_active_view &dofs_view,
    const int operation_type = moab::Interface::UNION
  ) const;

  PetscErrorCode getDataDofView(
    const DofEntity_multiIndex &dofs,DofEntity_multiIndex_active_view &dofs_view,
    const int operation_type = moab::Interface::UNION
  ) const;

  PetscErrorCode getRowDofView(
    const DofEntity_multiIndex &dofs,DofEntity_multiIndex_uid_view &dofs_view,
    const int operation_type = moab::Interface::UNION
  ) const;

  PetscErrorCode getColDofView(
    const DofEntity_multiIndex &dofs,DofEntity_multiIndex_uid_view &dofs_view,
    const int operation_type = moab::Interface::UNION
  ) const;

  PetscErrorCode getRowDofView(
    const NumeredDofEntity_multiIndex &dofs,NumeredDofEntity_multiIndex_uid_view_ordered &dofs_view,
    const int operation_type = moab::Interface::UNION
  ) const;

  PetscErrorCode getColDofView(
    const NumeredDofEntity_multiIndex &dofs,NumeredDofEntity_multiIndex_uid_view_ordered &dofs_view,
    const int operation_type = moab::Interface::UNION) const;

  PetscErrorCode getRowDofView(
    const NumeredDofEntity_multiIndex &dofs,NumeredDofEntity_multiIndex_idx_view_hashed &dofs_view,
    const int operation_type = moab::Interface::UNION
  ) const;

  PetscErrorCode getColDofView(
    const NumeredDofEntity_multiIndex &dofs,NumeredDofEntity_multiIndex_idx_view_hashed &dofs_view,
    const int operation_type = moab::Interface::UNION
  ) const;

  PetscErrorCode getElementAdjacency(
    const boost::shared_ptr<Field> field_ptr,Range &adjacency
  );

  inline const boost::shared_ptr<RefElement> getRefElement() const { return this->sPtr; }

  /**
   * \brief Get weak_ptr reference to sequence/vector storing dofs on entity.
   *
   * Vector is automatically destroy when last DOF in vector os destroyed. Every
   * shared_ptr to the DOF has aliased shared_ptr to vector of DOFs in that vector.
   * That do the trick.
   *
   */
  inline boost::weak_ptr<std::vector<FEDofEntity> >& getDofsSeqence() const {
    return dofsSequce;
  }

private:

  // Keep vector of DoFS on entity
  mutable boost::weak_ptr<std::vector<FEDofEntity> > dofsSequce;

};

/**
 * \brief interface to EntFiniteElement
 * \ingroup fe_multi_indices
 */
template <typename T>
struct interface_EntFiniteElement:
public
interface_FiniteElement<T>,
interface_RefElement<T> {

  interface_EntFiniteElement(const boost::shared_ptr<T> sptr):
  interface_FiniteElement<T>(sptr),
  interface_RefElement<T>(sptr) {
  };

  inline const FEDofEntity_multiIndex& getDataDofs() const { return this->sPtr->getDataDofs(); };

  /**
   * \brief Get number of DOFs on row
   * @return Number of dofs on row
   */
  inline DofIdx getNbDofsRow() const { return this->sPtr->getNbDofsRow(); }

  /**
   * \brief Get number of DOFs on col
   * @return Number of dofs on col
   */
  inline DofIdx getNbDofsCol() const { return this->sPtr->getNbDofsCol(); }

  /**
   * \brief Get number of DOFs on data
   * @return Number of dofs on data
   */
  inline DofIdx getNbDofsData() const { return this->sPtr->getNbDofsData(); }

  /**
   * \brief Get element entity
   * @return Element entity handle
   */
  inline EntityHandle getEnt() const { return this->sPtr->getRefEnt(); }

  /** \deprecated Use getEnt() instead
  */
  DEPRECATED inline EntityHandle get_ent() const { return getEnt(); }

  /**
   * \brief Get unique UId for finite element entity
   * @return UId
   */
  inline GlobalUId getGlobalUniqueId() const { return this->sPtr->getGlobalUniqueId(); }


  SideNumber_multiIndex &getSideNumberTable() const { return this->sPtr->getSideNumberTable(); }

  /** \deprecated Use getSideNumberTable() instead
  */
  DEPRECATED SideNumber_multiIndex &get_side_number_table() const {
    return this->sPtr->getSideNumberTable();
  }

  inline PetscErrorCode getElementAdjacency(const Field *field_ptr,Range &adjacency) {
    return this->getElementAdjacency(field_ptr,adjacency);
  }

  /** \deprecated Use getElementAdjacency() instead
  */
  DEPRECATED inline PetscErrorCode get_element_adjacency(
    const Field *field_ptr,Range &adjacency
  ) {
    return this->getElementAdjacency(field_ptr,adjacency);
  }

  inline const boost::shared_ptr<T> getRefElement() const { return this->sPtr->getRefElement(); }

};

/** \brief Partitioned (Indexed) Finite Element in Problem

  * This type of structure is used to compose problem. Problem is build from
  * indexed finite elements. This data structure carry information about
  * partition, which is specific to problem.


  * \ingroup fe_multi_indices
 */
struct NumeredEntFiniteElement: public interface_EntFiniteElement<EntFiniteElement> {

  typedef interface_FiniteElement<EntFiniteElement> interface_type_MoFEMFiniteElement;
  typedef interface_EntFiniteElement<EntFiniteElement> interface_type_EntFiniteElement;

  unsigned int part; ///< Partition number
  boost::shared_ptr<FENumeredDofEntity_multiIndex> rows_dofs;  ///< indexed dofs on rows
  boost::shared_ptr<FENumeredDofEntity_multiIndex> cols_dofs;  ///< indexed dofs on columns

  /**
   * \Construct indexed finite element
   */
  NumeredEntFiniteElement(const boost::shared_ptr<EntFiniteElement> sptr):
  interface_EntFiniteElement<EntFiniteElement>(sptr),
  part(-1),
  rows_dofs(boost::shared_ptr<FENumeredDofEntity_multiIndex>(new FENumeredDofEntity_multiIndex())),
  cols_dofs(boost::shared_ptr<FENumeredDofEntity_multiIndex>(new FENumeredDofEntity_multiIndex())) {
  };

  /**
   * \brief Get partition number
   * @return [description]
   */
  inline unsigned int getPart() const { return part; };

  /** \brief get FE dof on row
    * \ingroup mofem_dofs
    */
  inline const FENumeredDofEntity_multiIndex& getRowsDofs() const { return *(rows_dofs.get()); };

  /** \brief get FE dof on column
    * \ingroup mofem_dofs
    */
  inline const FENumeredDofEntity_multiIndex& getColsDofs() const { return *(cols_dofs.get()); };

  /** \brief get FE dof by petsc index
    * \ingroup mofem_dofs
    */
  PetscErrorCode getRowDofsByPetscGlobalDofIdx(DofIdx idx,const FENumeredDofEntity **dof_ptr) const;

  /** \deprecated Use getRowDofsByPetscGlobalDofIdx() instead
  */
  inline DEPRECATED PetscErrorCode get_row_dofs_by_petsc_gloabl_dof_idx(
    DofIdx idx,const FENumeredDofEntity **dof_ptr
  ) const {
    return getRowDofsByPetscGlobalDofIdx(idx,dof_ptr);
  }

  /** \brief get FE dof by petsc index
    * \ingroup mofem_dofs
    */
  PetscErrorCode getColDofsByPetscGlobalDofIdx(DofIdx idx,const FENumeredDofEntity **dof_ptr) const;

  /** \deprecated Use getColDofsByPetscGlobalDofIdx() instead
  */
  inline DEPRECATED PetscErrorCode get_col_dofs_by_petsc_gloabl_dof_idx(
    DofIdx idx,const FENumeredDofEntity **dof_ptr
  ) const {
    return getColDofsByPetscGlobalDofIdx(idx,dof_ptr);
  }

  friend std::ostream& operator<<(std::ostream& os,const NumeredEntFiniteElement& e) {
    os << "part " << e.part << " " << *(e.sFePtr);
    return os;
  }

  /**
   * \brief Get weak_ptr reference to sequence/vector storing dofs on entity.
   *
   * Vector is automatically destroy when last DOF in vector os destroyed. Every
   * shared_ptr to the DOF has aliased shared_ptr to vector of DOFs in that vector.
   * That do the trick.
   *
   */
  inline boost::weak_ptr<std::vector<FENumeredDofEntity> >& getRowDofsSeqence() const {
    return dofsRowSequce;
  }

  /**
   * \brief Get weak_ptr reference to sequence/vector storing dofs on entity.
   *
   * Vector is automatically destroy when last DOF in vector os destroyed. Every
   * shared_ptr to the DOF has aliased shared_ptr to vector of DOFs in that vector.
   * That do the trick.
   *
   */
  inline boost::weak_ptr<std::vector<FENumeredDofEntity> >& getColDofsSeqence() const {
    return dofsColSequce;
  }

private:

  // Keep vector of DoFS on entity
  mutable boost::weak_ptr<std::vector<FENumeredDofEntity> > dofsRowSequce;
  mutable boost::weak_ptr<std::vector<FENumeredDofEntity> > dofsColSequce;

};

/** \brief interface for NumeredEntFiniteElement
 * \ingroup fe_multi_indices
 */
template <typename T>
struct interface_NumeredEntFiniteElement: public interface_EntFiniteElement<T> {

  interface_NumeredEntFiniteElement(const boost::shared_ptr<T> sptr): interface_EntFiniteElement<T>(sptr) {};

  /**
   * \brief Get partition number
   * @return Partition number
   */
  inline unsigned int getPart() const { return this->sPtr->getPart(); }

  /** \brief get FE dof on row
    * \ingroup mofem_dofs
    */
  inline const FENumeredDofEntity_multiIndex& getRowsDofs() const {
    return this->sPtr->getRowsDofs();
  };

  /** \brief get FE dof on column
    * \ingroup mofem_dofs
    */
  inline const FENumeredDofEntity_multiIndex& getColsDofs() const {
    return this->sPtr->getColsDofs();
  };

};

/**
 * @relates multi_index_container
 * \brief MultiIndex container for EntFiniteElement
 * \ingroup fe_multi_indices
 *
 */
typedef multi_index_container<
  boost::shared_ptr<EntFiniteElement>,
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>,
      member<EntFiniteElement,GlobalUId,&EntFiniteElement::global_uid>
    >,
    ordered_non_unique<
      tag<Ent_mi_tag>,
      const_mem_fun<EntFiniteElement,EntityHandle,&EntFiniteElement::getEnt>
    >,
    ordered_non_unique<
      tag<FiniteElement_name_mi_tag>,
      const_mem_fun<
        EntFiniteElement::interface_type_MoFEMFiniteElement,
        boost::string_ref,
        &EntFiniteElement::getNameRef
      >
    >,
    ordered_non_unique<
      tag<BitFEId_mi_tag>,
      const_mem_fun<
        EntFiniteElement::interface_type_MoFEMFiniteElement,
        BitFEId,
        &EntFiniteElement::getId
      >,
      LtBit<BitFEId>
    >,
    ordered_non_unique<
      tag<EntType_mi_tag>,
      const_mem_fun<
        EntFiniteElement::interface_type_RefEntity,
        EntityType,
        &EntFiniteElement::getEntType
      >
    >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_mi_tag>,
      composite_key<
	      EntFiniteElement,
	      const_mem_fun<
          EntFiniteElement::interface_type_MoFEMFiniteElement,
          boost::string_ref,
          &EntFiniteElement::getNameRef
        >,
	      const_mem_fun<EntFiniteElement,EntityHandle,&EntFiniteElement::getEnt>
      >
    >
  >
> EntFiniteElement_multiIndex;

/**
 *  \brief Entity finite element multi-index by finite element name
 *
 *  \ingroup fe_multi_indices
 */
typedef EntFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type EntFiniteElementbyName;

/**
  @relates multi_index_container
  \brief MultiIndex for entities for NumeredEntFiniteElement
  \ingroup fe_multi_indices
 */
typedef multi_index_container<
  boost::shared_ptr<NumeredEntFiniteElement>,
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>,
      const_mem_fun<
        NumeredEntFiniteElement::interface_type_EntFiniteElement,
        GlobalUId,
        &NumeredEntFiniteElement::getGlobalUniqueId
      >
    >,
    ordered_non_unique<
      tag<Part_mi_tag>,
      member<NumeredEntFiniteElement,unsigned int,&NumeredEntFiniteElement::part>
    >,
    ordered_non_unique<
      tag<FiniteElement_name_mi_tag>,
      const_mem_fun<
        NumeredEntFiniteElement::interface_type_MoFEMFiniteElement,
        boost::string_ref,
        &NumeredEntFiniteElement::getNameRef
      >
    >,
    ordered_non_unique<
      tag<FiniteElement_Part_mi_tag>,
      member<NumeredEntFiniteElement,unsigned int,&NumeredEntFiniteElement::part>
    >,
    ordered_non_unique<
      tag<Ent_mi_tag>,
      const_mem_fun<
        NumeredEntFiniteElement::interface_type_EntFiniteElement,
        EntityHandle,
        &NumeredEntFiniteElement::getEnt
      >
    >,
    ordered_non_unique<
      tag<Composite_Name_And_Ent_mi_tag>,
      composite_key<
        NumeredEntFiniteElement,
        const_mem_fun<
          NumeredEntFiniteElement::interface_type_MoFEMFiniteElement,
          boost::string_ref,
          &NumeredEntFiniteElement::getNameRef
        >,
        const_mem_fun<
          NumeredEntFiniteElement::interface_type_EntFiniteElement,
          EntityHandle,
          &NumeredEntFiniteElement::getEnt
        >
      >
    >,
    ordered_non_unique<
      tag<Composite_Name_And_Part_mi_tag>,
      composite_key<
        NumeredEntFiniteElement,
        const_mem_fun<
          NumeredEntFiniteElement::interface_type_MoFEMFiniteElement,
          boost::string_ref,
          &NumeredEntFiniteElement::getNameRef
        >,
        member<NumeredEntFiniteElement,unsigned int,&NumeredEntFiniteElement::part>
      >
    >
  >
> NumeredEntFiniteElement_multiIndex;

/**
 *  \brief Entity finite element multi-index by finite element name
 *
 *  \ingroup fe_multi_indices
 */
typedef NumeredEntFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type
NumeredEntFiniteElementbyName;


/**
  *  \brief Entity finite element multi-index by finite element name and partition
  *
  *  \ingroup fe_multi_indices
 */
typedef NumeredEntFiniteElement_multiIndex::index<Composite_Name_And_Part_mi_tag>::type
NumeredEntFiniteElementbyNameAndPart;


/**
  @relates multi_index_container
  \brief MultiIndex for entities for FiniteElement
  \ingroup fe_multi_indices
 */
typedef multi_index_container<
  boost::shared_ptr<FiniteElement>,
  indexed_by<
    hashed_unique<
      tag<FiniteElement_Meshset_mi_tag>,
      member<FiniteElement,EntityHandle,&FiniteElement::meshset>
    >,
    hashed_unique<
      tag<BitFEId_mi_tag>,
      const_mem_fun<FiniteElement,BitFEId,&FiniteElement::getId>, HashBit<BitFEId>, EqBit<BitFEId>
    >,
    ordered_unique<
      tag<FiniteElement_name_mi_tag>,
      const_mem_fun<FiniteElement,boost::string_ref,&FiniteElement::getNameRef>
    >
  >
> FiniteElement_multiIndex;

// modificators

/**
 * \brief Change finite element part
 *
 * \ingroup fe_multi_indices
 */
struct NumeredEntFiniteElement_change_part {
  unsigned int pArt;
  NumeredEntFiniteElement_change_part(unsigned int part): pArt(part) {};
  void operator()(boost::shared_ptr<NumeredEntFiniteElement> &fe) {
    fe->part = pArt;
  }
};

/**
 * \brief Add field to column
 *
 * \ingroup fe_multi_indices
 */
struct MoFEMFiniteElement_col_change_bit_add {
  BitFieldId fIdCol;
  MoFEMFiniteElement_col_change_bit_add(const BitFieldId f_id_col): fIdCol(f_id_col) {};
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

/**
 * \brief Add field to row
 *
 * \ingroup fe_multi_indices
 */
struct MoFEMFiniteElement_row_change_bit_add {
  BitFieldId fIdRow;
  MoFEMFiniteElement_row_change_bit_add(const BitFieldId f_id_row): fIdRow(f_id_row) {};
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

/**
 * \brief Add field to data
 *
 * \ingroup fe_multi_indices
 */
struct MoFEMFiniteElement_change_bit_add {
  BitFieldId fIdData;
  MoFEMFiniteElement_change_bit_add(const BitFieldId f_id_data): fIdData(f_id_data) {};
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

/**
 * \brief Unset field from column
 *
 * \ingroup fe_multi_indices
 */
struct MoFEMFiniteElement_col_change_bit_off {
  BitFieldId fIdCol;
  MoFEMFiniteElement_col_change_bit_off(const BitFieldId f_id_col): fIdCol(f_id_col) {};
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

/**
 * \brief Unset field from row
 *
 * \ingroup fe_multi_indices
 */
struct MoFEMFiniteElement_row_change_bit_off {
  BitFieldId fIdRow;
  MoFEMFiniteElement_row_change_bit_off(const BitFieldId f_id_row): fIdRow(f_id_row) {};
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

/**
 * \brief Unset field from data
 *
 * \ingroup fe_multi_indices
 */
struct MoFEMFiniteElement_change_bit_off {
  BitFieldId fIdData;
  MoFEMFiniteElement_change_bit_off(const BitFieldId f_id_data): fIdData(f_id_data) {};
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

}

/**
 * Loop over DOFs in row on element
 * @param  FEPTR pointer to element structure \ref NumeredEntFiniteElement
 * @param  IT    iterator
 * @return       user return in for(_IT_FENUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(FEPTR,IT))
 * \ingroup fe_multi_indices
 */
#define _IT_FENUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(FEPTR,IT) \
FENumeredDofEntity_multiIndex::iterator IT = FEPTR->rows_dofs->begin(); \
IT!=FEPTR->rows_dofs->end(); IT++

/**
 * Loop over DOFs in col on element
 * @param  FEPTR pointer to element structure \ref NumeredEntFiniteElement
 * @param  IT    iterator
 * @return       user return in for(_IT_FENUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(FEPTR,IT))
 * \ingroup fe_multi_indices
 */
#define _IT_FENUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(FEPTR,IT) \
FENumeredDofEntity_multiIndex::iterator IT = FEPTR->cols_dofs->begin(); \
IT!=FEPTR->cols_dofs->end(); IT++

/**
 * Loop over DOFs in row on element for particular filed
 * @param  FEPTR pointer to element structure \ref NumeredEntFiniteElement
 * @param  NAME  name of filed
 * @param  IT    iterator
 * @return       user return in for(_IT_FENUMEREDDOFMOFEMENTITY_BY_NAME_ROW_FOR_LOOP_(FEPTR,NAME,IT))
 * \ingroup fe_multi_indices
 */
#define _IT_FENUMEREDDOFMOFEMENTITY_BY_NAME_ROW_FOR_LOOP_(FEPTR,NAME,IT) \
FENumeredDofEntityByFieldName::iterator IT = FEPTR->rows_dofs->get<FieldName_mi_tag>().lower_bound(NAME); \
IT!=FEPTR->rows_dofs->get<FieldName_mi_tag>().upper_bound(NAME); IT++

/**
 * Loop over DOFs in col on element for particular filed
 * @param  FEPTR pointer to element structure \ref NumeredEntFiniteElement
 * @param  NAME  name of filed
 * @param  IT    iterator
 * @return       user return in for(_IT_FENUMEREDDOFMOFEMENTITY_BY_NAME_COL_FOR_LOOP_(FEPTR,NAME,IT))
 * \ingroup fe_multi_indices
 */
#define _IT_FENUMEREDDOFMOFEMENTITY_BY_NAME_COL_FOR_LOOP_(FEPTR,NAME,IT) \
FENumeredDofEntityByFieldName::iterator IT = FEPTR->cols_dofs->get<FieldName_mi_tag>().lower_bound(NAME); \
IT!=FEPTR->cols_dofs->get<FieldName_mi_tag>().upper_bound(NAME); IT++

#endif // __FEMMULTIINDICES_HPP__

/***************************************************************************//**
 * \defgroup fe_multi_indices Finite elements structures and multi-indices
 * \ingroup mofem
 ******************************************************************************/
