/** \file FEMultiIndices.hpp
 * \brief Multi-index contains, data structures for mofem finite elements and
 * other low-level functions
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
struct RefElement : public interface_RefEntity<RefEntity> {

  typedef interface_RefEntity<RefEntity> interface_type_RefEntity;

  static BitRefEdges DummyBitRefEdges;

  SideNumber_multiIndex side_number_table;
  RefElement(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElement() = default;

  virtual const BitRefEdges &getBitRefEdges() const { return DummyBitRefEdges; }

  virtual int getBitRefEdgesUlong() const { return 0; }

  SideNumber_multiIndex &getSideNumberTable() const {
    return const_cast<SideNumber_multiIndex &>(side_number_table);
  }

  static const boost::shared_ptr<SideNumber> nullSideNumber;

  virtual const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const {
    NOT_USED(ent);
    return nullSideNumber;
  };

  /**
   * \brief Get pointer to RefEntity
   */
  inline boost::shared_ptr<RefEntity> &getRefEntityPtr() const {
    return this->sPtr;
  }

  friend std::ostream &operator<<(std::ostream &os, const RefElement &e);
};

/**
 * \brief keeps data about abstract MESHSET finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_MESHSET : public RefElement {
  RefElement_MESHSET(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElement_MESHSET() = default;
  const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const;
};
/**
 * \brief keeps data about abstract PRISM finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_PRISM : public RefElement {
  BitRefEdges *tag_BitRefEdges;
  RefElement_PRISM(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElement_PRISM() = default;

  const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const;
  const BitRefEdges &getBitRefEdges() const { return *tag_BitRefEdges; }
  int getBitRefEdgesUlong() const { return getBitRefEdges().to_ulong(); }
};

/**
 * \brief keeps data about abstract TET finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_TET : public RefElement {
  BitRefEdges *tag_BitRefEdges;
  const int *tag_type_data;
  RefElement_TET(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElement_TET() = default;

  const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const;
  SideNumber_multiIndex &getSideNumberTable() const {
    return const_cast<SideNumber_multiIndex &>(side_number_table);
  };
  const BitRefEdges &getBitRefEdges() const { return *tag_BitRefEdges; }
  int getBitRefEdgesUlong() const { return getBitRefEdges().to_ulong(); }
};

/**
 * \brief keeps data about abstract TRI finite element
 * \ingroup fe_multi_indices
 */
struct RefElementFace : public RefElement {
  RefElementFace(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElementFace() = default;
  const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const;
  friend std::ostream &operator<<(std::ostream &os, const RefElementFace &e);
};

/**
 * \brief keeps data about abstract EDGE finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_EDGE : public RefElement {
  RefElement_EDGE(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElement_EDGE() = default;
  const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const;
  friend std::ostream &operator<<(std::ostream &os, const RefElement_EDGE &e);
};

/**
 * \brief keeps data about abstract VERTEX finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_VERTEX : public RefElement {
  RefElement_VERTEX(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElement_VERTEX() = default;
  const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const;
  friend std::ostream &operator<<(std::ostream &os, const RefElement_VERTEX &e);
};

/**
 * \brief intrface to RefElement
 * \ingroup fe_multi_indices
 */
template <typename T> struct interface_RefElement : interface_RefEntity<T> {

  typedef interface_RefEntity<T> interface_type_RefEntity;
  typedef interface_RefElement<T> interface_type_RefElement;

  interface_RefElement(const boost::shared_ptr<T> &sptr)
      : interface_RefEntity<T>(sptr) {}
  virtual ~interface_RefElement() = default;

  inline int getBitRefEdgesUlong() const {
    return this->sPtr->getBitRefEdgesUlong();
  }

  inline SideNumber_multiIndex &getSideNumberTable() const {
    return this->sPtr->getSideNumberTable();
  }

  inline const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const {
    return this->sPtr->getSideNumberPtr(ent);
  }

  inline boost::shared_ptr<RefEntity> &getRefEntityPtr() const {
    return this->sPtr->getRefEntityPtr();
  }

  inline const boost::shared_ptr<T> &getRefElement() const {
    return this->sPtr;
  }
};

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
    boost::shared_ptr<RefElement>,
    // ptrWrapperRefElement,
    indexed_by<
        ordered_unique<tag<Ent_mi_tag>,
                       const_mem_fun<RefElement::interface_type_RefEntity,
                                     EntityHandle, &RefElement::getRefEnt>>,
        ordered_non_unique<tag<EntType_mi_tag>,
                           const_mem_fun<RefElement::interface_type_RefEntity,
                                         EntityType, &RefElement::getEntType>>>>
    RefElement_multiIndex;

typedef multi_index_container<
    boost::shared_ptr<RefElement>,
    // ptrWrapperRefElement,
    indexed_by<
        ordered_unique<tag<Ent_mi_tag>,
                       const_mem_fun<RefElement::interface_type_RefEntity,
                                     EntityHandle, &RefElement::getRefEnt>>,
        ordered_non_unique<
            tag<Ent_Ent_mi_tag>,
            const_mem_fun<RefElement::interface_type_RefEntity, EntityHandle,
                          &RefElement::getParentEnt>>,
        ordered_non_unique<
            tag<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>,
            composite_key<
                RefElement,
                const_mem_fun<RefElement::interface_type_RefEntity,
                              EntityHandle, &RefElement::getParentEnt>,
                const_mem_fun<RefElement, int,
                              &RefElement::getBitRefEdgesUlong>>>>>
    RefElement_multiIndex_parents_view;

struct EntFiniteElement;

/** \brief user adjacency function
 * \ingroup fe_multi_indices
 */
typedef boost::function<MoFEMErrorCode(Interface &moab, const Field &field,
                                       const EntFiniteElement &fe,
                                       Range &adjacency)>
    ElementAdjacencyFunct;

/**
 * \brief Finite element definition
 * \ingroup fe_multi_indices
 */
struct FiniteElement {

  EntityHandle meshset;                ///< meshset stores FE ents
  BitFEId *tagId;                      ///< ptr to tag storing FE id
  void *tagName;                       ///< ptr to tag storing FE name
  int tagNameSize;                     ///< numer of characters in FE name
  BitFieldId *tag_BitFieldId_col_data; ///< tag stores col id_id for fields
  BitFieldId *tag_BitFieldId_row_data; ///< tag stores row id_id for fields
  BitFieldId *tag_BitFieldId_data;     ///< tag stores data id_id for fields
  UId feUId;

  FiniteElement(Interface &moab, const EntityHandle _meshset);

  /**
   * @brief Get finite element uid
   *
   * @return const UId&
   */
  inline const UId &getFEUId() const { return feUId; }

  /**
   * \brief Get finite element id
   * @return Finite element Id
   */
  inline BitFEId getId() const { return *tagId; }

  /**
   * \brief Get meshset containing element entities
   * @return Meshset
   */
  inline EntityHandle getMeshset() const { return meshset; }

  /**
   * \brief Get finite element name
   * @return string_ref
   */
  inline boost::string_ref getNameRef() const {
    return boost::string_ref((char *)tagName, tagNameSize);
  }

  /**
   * \brief Get finite element name
   * @return string
   */
  inline std::string getName() const {
    return std::string((char *)tagName, tagNameSize);
  }

  /**
   * \brief Get field ids on columns
   * @return Bit field ids
   */
  inline BitFieldId getBitFieldIdCol() const {
    return *((BitFieldId *)tag_BitFieldId_col_data);
  }

  /**
   * \brief Get field ids on rows
   * @return Bit field ids
   */
  inline BitFieldId getBitFieldIdRow() const {
    return *((BitFieldId *)tag_BitFieldId_row_data);
  }

  /**
   * \brief Get field ids on data
   * @return Bit field ids
   */
  inline BitFieldId getBitFieldIdData() const {
    return *((BitFieldId *)tag_BitFieldId_data);
  }

  /**
   * \brief Get bit identifying this element
   *
   * Each element like field is identified by bit set. Each element has unique
   * bit set, this function returns number of that bit.
   *
   * @return Bit number
   */
  inline unsigned int getBitNumber() const {
    return ffsl(((BitFieldId *)tagId)->to_ulong());
  }

  /**
   * \brief Table of functions retrieving adjacencies for finite element
   * User can alter and change default behavior
   */
  std::array<ElementAdjacencyFunct, MBMAXTYPE> elementAdjacencyTable;

  /**
   * \brief print finite element
   */
  friend std::ostream &operator<<(std::ostream &os, const FiniteElement &e);
};

/** \brief default adjacency map
 * \ingroup fe_multi_indices
 */
struct DefaultElementAdjacency {

  static MoFEMErrorCode defaultVertex(Interface &moab, const Field &field,
                                      const EntFiniteElement &fe,
                                      Range &adjacency);
  static MoFEMErrorCode defaultEdge(Interface &moab, const Field &field,
                                    const EntFiniteElement &fe,
                                    Range &adjacency);
  static MoFEMErrorCode defaultFace(Interface &moab, const Field &field,
                                    const EntFiniteElement &fe,
                                    Range &adjacency);
  static MoFEMErrorCode defaultTet(Interface &moab, const Field &field,
                                   const EntFiniteElement &fe,
                                   Range &adjacency);
  static MoFEMErrorCode defaultPrism(Interface &moab, const Field &field,
                                     const EntFiniteElement &fe,
                                     Range &adjacency);
  static MoFEMErrorCode defaultMeshset(Interface &moab, const Field &field,
                                       const EntFiniteElement &fe,
                                       Range &adjacency);
};

/**
 * \brief Inetface for FE
 * \ingroup fe_multi_indices
 */
template <typename T> struct interface_FiniteElement {

  mutable boost::shared_ptr<T> sFePtr;

  interface_FiniteElement(const boost::shared_ptr<T> &ptr) : sFePtr(ptr){};
  virtual ~interface_FiniteElement() = default;

  inline const boost::shared_ptr<FiniteElement> &get_MoFEMFiniteElementPtr() {
    return this->sFePtr;
  };

  /**
   * @copydoc MoFEM::FiniteElement::getFEUId
   */
  inline const UId &getFEUId() const { return this->sFePtr->getFEUId(); }

  /**
   * @copydoc MoFEM::FiniteElement::getId
   */
  inline BitFEId getId() const { return this->sFePtr->getId(); }

  /**
   * @copydoc MoFEM::FiniteElement::getMeshset
   */
  inline EntityHandle getMeshset() const { return this->sFePtr->getMeshset(); }

  /**
   * @copydoc MoFEM::FiniteElement::getNameRef
   */
  inline boost::string_ref getNameRef() const {
    return this->sFePtr->getNameRef();
  }

  /**
   * @copydoc MoFEM::FiniteElement::getName
   */
  inline std::string getName() const { return this->sFePtr->getName(); }

  /**
   * @copydoc MoFEM::FiniteElement::getBitFieldIdCol
   */
  inline BitFieldId getBitFieldIdCol() const {
    return this->sFePtr->getBitFieldIdCol();
  }

  /**
   * @copydoc MoFEM::FiniteElement::getBitFieldIdRow
   */
  inline BitFieldId getBitFieldIdRow() const {
    return this->sFePtr->getBitFieldIdRow();
  }

  /**
   * @copydoc MoFEM::FiniteElement::getBitFieldIdData
   */
  inline BitFieldId getBitFieldIdData() const {
    return this->sFePtr->getBitFieldIdData();
  }

  /**
   * @copydoc MoFEM::FiniteElement::getBitNumber
   */
  inline unsigned int getBitNumber() const {
    return this->sFePtr->getBitNumber();
  }
};

/**
 * \brief Finite element data for entity
 * \ingroup fe_multi_indices
 */
struct EntFiniteElement : public interface_FiniteElement<FiniteElement>,
                          interface_RefElement<RefElement> {

  typedef interface_RefEntity<RefElement> interface_type_RefEntity;
  typedef interface_RefElement<RefElement> interface_type_RefElement;
  typedef interface_FiniteElement<FiniteElement> interface_type_FiniteElement;

  EntFiniteElement(const boost::shared_ptr<RefElement> &ref_finite_element,
                   const boost::shared_ptr<FiniteElement> &fe_ptr);
  virtual ~EntFiniteElement() = default;

  /**
   * \brief Get unique UId for finite element entity
   * @return UId
   */
  inline UId getGlobalUniqueId() const { return getGlobalUniqueIdCalculate(); }

  static inline UId getGlobalUniqueIdCalculate(const EntityHandle ent,
                                               UId fe_uid) {
    return fe_uid |= ent;
  }

  /**
   * \brief Generate UId for finite element entity
   * @return finite element entity unique Id
   */
  inline UId getGlobalUniqueIdCalculate() const {
    return getGlobalUniqueIdCalculate(getEnt(), getFEUId());
  }

  /**
   * \brief Get element entity
   * @return Element entity handle
   */
  inline EntityHandle getEnt() const { return getRefEnt(); }

  /**
   * \brief Get number of DOFs on data
   * @return Number of dofs on data
   */
  inline DofIdx getNbDofsData() const { return dataDofs->size(); }

  /**
   * \brief Get data data dos multi-index structure
   * @return Reference multi-index FEDofEntity_multiIndex
   */
  inline const FEDofEntity_multiIndex &getDataDofs() const {
    RefEntityTmp<0>::refElementPtr = this->getRefElement();
    return *dataDofs;
  };

  inline boost::shared_ptr<FEDofEntity_multiIndex> &getDataDofsPtr() {
    RefEntityTmp<0>::refElementPtr = this->getRefElement();
    return dataDofs;
  };

  inline const FieldEntity_multiIndex_spaceType_view &
  getDataFieldEntsView() const {
    return *dataFieldEntsView;
  };

  inline boost::shared_ptr<FieldEntity_multiIndex_spaceType_view> &
  getDataFieldEntsViewPtr() {
    return dataFieldEntsView;
  };

  inline const FieldEntity_vector_view &getRowFieldEntsView() const {
    return *rowFieldEntsView;
  };

  inline boost::shared_ptr<FieldEntity_vector_view> &getRowFieldEntsViewPtr() {
    return rowFieldEntsView;
  };

  inline const FieldEntity_vector_view &getColFieldEntsView() const {
    return *colFieldEntsView;
  };

  inline boost::shared_ptr<FieldEntity_vector_view> &getColFieldEntsViewPtr() {
    return colFieldEntsView;
  };

  friend std::ostream &operator<<(std::ostream &os, const EntFiniteElement &e);

  template <typename FE_ENTS, typename MOFEM_DOFS, typename MOFEM_DOFS_VIEW>
  inline MoFEMErrorCode
  getDofView(const FE_ENTS &fe_ents_view, const MOFEM_DOFS &mofem_dofs,
             MOFEM_DOFS_VIEW &dofs_view, const int operation_type) {
    MoFEMFunctionBeginHot;
    if (operation_type == moab::Interface::UNION) {
      for (auto &it : fe_ents_view) {
        if (auto e = it.lock()) {
          auto r = mofem_dofs.template get<Unique_Ent_mi_tag>().equal_range(
              e->getGlobalUniqueId());
          dofs_view.insert(r.first, r.second);
        }
      }
    } else
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
    MoFEMFunctionReturnHot(0);
  }

  template <typename MOFEM_DOFS, typename MOFEM_DOFS_VIEW>
  inline MoFEMErrorCode
  getRowDofView(const MOFEM_DOFS &mofem_dofs, MOFEM_DOFS_VIEW &dofs_view,
                const int operation_type = moab::Interface::UNION) {
    return getDofView(getRowFieldEntsView(), mofem_dofs, dofs_view,
                      operation_type);
  }

  template <typename MOFEM_DOFS, typename MOFEM_DOFS_VIEW>
  inline MoFEMErrorCode
  getColDofView(const MOFEM_DOFS &mofem_dofs, MOFEM_DOFS_VIEW &dofs_view,
                const int operation_type = moab::Interface::UNION) {
    return getDofView(getColFieldEntsView(), mofem_dofs, dofs_view,
                      operation_type);
  }

  MoFEMErrorCode getElementAdjacency(const boost::shared_ptr<Field> field_ptr,
                                     Range &adjacency);

  inline boost::shared_ptr<RefElement> &getRefElement() const {
    return this->sPtr;
  }

  /**
   * \brief Get weak_ptr reference to sequence/vector storing dofs on entity.
   *
   * Vector is automatically destroy when last DOF in vector os destroyed. Every
   * shared_ptr to the DOF has aliased shared_ptr to vector of DOFs in that
   * vector. That do the trick.
   *
   */
  inline boost::weak_ptr<std::vector<FEDofEntity>> &getDofsSequence() const {
    return dofsSequce;
  }

protected:
  boost::shared_ptr<FEDofEntity_multiIndex> dataDofs;
  boost::shared_ptr<FieldEntity_multiIndex_spaceType_view> dataFieldEntsView;
  boost::shared_ptr<FieldEntity_vector_view> rowFieldEntsView;
  boost::shared_ptr<FieldEntity_vector_view> colFieldEntsView;

private:
  // Keep vector of DoFS on entity
  mutable boost::weak_ptr<std::vector<FEDofEntity>> dofsSequce;
};

/**
 * \brief interface to EntFiniteElement
 * \ingroup fe_multi_indices
 */
template <typename T>
struct interface_EntFiniteElement : public interface_FiniteElement<T>,
                                    interface_RefElement<T> {

  interface_EntFiniteElement(const boost::shared_ptr<T> &sptr)
      : interface_FiniteElement<T>(sptr), interface_RefElement<T>(sptr) {}
  virtual ~interface_EntFiniteElement() = default;

  inline const FEDofEntity_multiIndex &getDataDofs() const {
    return this->sPtr->getDataDofs();
  }

  inline boost::shared_ptr<FEDofEntity_multiIndex> &getDataDofsPtr() {
    return this->sPtr->getDataDofsPtr();
  }

  inline const FieldEntity_multiIndex_spaceType_view &
  getDataFieldEntsView() const {
    return this->sPtr->getDataFieldEntsView();
  };

  inline boost::shared_ptr<FieldEntity_multiIndex_spaceType_view> &
  getDataFieldEntsViewPtr() {
    return this->sPtr->getDataFieldEntsViewPtr();
  };

  inline const FieldEntity_vector_view &getRowFieldEntsView() const {
    return this->sPtr->getRowFieldEntsView();
  };

  inline boost::shared_ptr<FieldEntity_vector_view> &getRowFieldEntsViewPtr() {
    return this->sPtr->getRowFieldEntsViewPtr();
  }

  inline const FieldEntity_vector_view &getColFieldEntsView() const {
    return this->sPtr->getColFieldEntsView();
  };

  inline boost::shared_ptr<FieldEntity_vector_view> &getColFieldEntsViewPtr() {
    return this->sPtr->getColFieldEntsViewPtr();
  };

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

  // /** \deprecated Use getEnt() instead
  // */
  // DEPRECATED inline EntityHandle get_ent() const { return getEnt(); }

  /**
   * \brief Get unique UId for finite element entity
   * @return UId
   */
  inline UId getGlobalUniqueId() const {
    return this->sPtr->getGlobalUniqueId();
  }

  SideNumber_multiIndex &getSideNumberTable() const {
    return this->sPtr->getSideNumberTable();
  }

  inline MoFEMErrorCode getElementAdjacency(const Field *field_ptr,
                                            Range &adjacency) {
    return this->getElementAdjacency(field_ptr, adjacency);
  }

  inline boost::shared_ptr<RefElement> &getRefElement() const {
    return this->sPtr->getRefElement();
  }
};

/** \brief Partitioned (Indexed) Finite Element in Problem

  * This type of structure is used to compose problem. Problem is build from
  * indexed finite elements. This data structure carry information about
  * partition, which is specific to problem.


  * \ingroup fe_multi_indices
 */
struct NumeredEntFiniteElement
    : public interface_EntFiniteElement<EntFiniteElement> {

  virtual ~NumeredEntFiniteElement() = default;

  typedef interface_FiniteElement<EntFiniteElement>
      interface_type_FiniteElement;
  typedef interface_EntFiniteElement<EntFiniteElement>
      interface_type_EntFiniteElement;

  unsigned int part; ///< Partition number

  inline boost::shared_ptr<EntFiniteElement> &getEntFiniteElement() const {
    return this->sPtr;
  }

  /**
   * \Construct indexed finite element
   */
  NumeredEntFiniteElement(const boost::shared_ptr<EntFiniteElement> &sptr);

  /**
   * \brief Get partition number
   * @return [description]
   */
  inline unsigned int getPart() const { return part; };

  /** \brief get FE dof on row
   * \ingroup mofem_dofs
   */
  inline const FENumeredDofEntity_multiIndex &getRowDofs() const {
    RefEntityTmp<0>::refElementPtr = this->getRefElement();
    return *rowDofs;
  }

  inline boost::shared_ptr<FENumeredDofEntity_multiIndex> &getRowDofsPtr() {
    RefEntityTmp<0>::refElementPtr = this->getRefElement();
    return rowDofs;
  }

  /** \brief get FE dof on column
   * \ingroup mofem_dofs
   */
  inline const FENumeredDofEntity_multiIndex &getColDofs() const {
    RefEntityTmp<0>::refElementPtr = this->getRefElement();
    return *colDofs;
  }

  inline boost::shared_ptr<FENumeredDofEntity_multiIndex> &getColDofsPtr() {
    RefEntityTmp<0>::refElementPtr = this->getRefElement();
    return colDofs;
  }

  /** \brief get FE dof by petsc index
   * \ingroup mofem_dofs
   */
  boost::weak_ptr<FENumeredDofEntity>
  getRowDofsByPetscGlobalDofIdx(const int idx) const;

  /** \brief get FE dof by petsc index
   * \ingroup mofem_dofs
   */
  boost::weak_ptr<FENumeredDofEntity>
  getColDofsByPetscGlobalDofIdx(const int idx) const;

  /**
   * @deprecated Unsafe to use, will be removed in future releases.
   *
   * Get the Row Dofs By Petsc Global Dof Idx object
   *
   * @param idx
   * @param dof_raw_ptr
   * @return MoFEMErrorCode
   */
  DEPRECATED inline MoFEMErrorCode
  getRowDofsByPetscGlobalDofIdx(const int idx,
                                const FENumeredDofEntity **dof_raw_ptr) const {
    MoFEMFunctionBegin;
    if (auto r = getRowDofsByPetscGlobalDofIdx(idx).lock())
      *dof_raw_ptr = r.get();
    else
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
               "dof which index < %d > not found", idx);
    MoFEMFunctionReturn(0);
  }

  /**
   * @deprecated Unsafe to use, will be removed in future releases.
   *
   * Get the Row Dofs By Petsc Global Dof Idx object
   *
   * @param idx
   * @param dof_raw_ptr
   * @return MoFEMErrorCode
   */
  DEPRECATED inline MoFEMErrorCode
  getColDofsByPetscGlobalDofIdx(const int idx,
                                const FENumeredDofEntity **dof_raw_ptr) const {
    MoFEMFunctionBegin;
    if (auto r = getColDofsByPetscGlobalDofIdx(idx).lock())
      *dof_raw_ptr = r.get();
    else
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
               "dof which index < %d > not found", idx);
    MoFEMFunctionReturn(0);
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const NumeredEntFiniteElement &e) {
    os << "part " << e.part << " " << *(e.sFePtr);
    return os;
  }

  /**
   * \brief Get weak_ptr reference to sequence/vector storing dofs on entity.
   *
   * Vector is automatically destroy when last DOF in vector os destroyed. Every
   * shared_ptr to the DOF has aliased shared_ptr to vector of DOFs in that
   * vector. That do the trick.
   *
   */
  inline boost::weak_ptr<std::vector<FENumeredDofEntity>> &
  getRowDofsSequence() const {
    return dofsRowSequce;
  }

  /**
   * \brief Get weak_ptr reference to sequence/vector storing dofs on entity.
   *
   * Vector is automatically destroy when last DOF in vector os destroyed. Every
   * shared_ptr to the DOF has aliased shared_ptr to vector of DOFs in that
   * vector. That do the trick.
   *
   */
  inline boost::weak_ptr<std::vector<FENumeredDofEntity>> &
  getColDofsSequence() const {
    return dofsColSequce;
  }

protected:
  boost::shared_ptr<FENumeredDofEntity_multiIndex>
      rowDofs; ///< indexed dofs on rows
  boost::shared_ptr<FENumeredDofEntity_multiIndex>
      colDofs; ///< indexed dofs on columns

private:
  // Keep vector of DoFS on entity
  mutable boost::weak_ptr<std::vector<FENumeredDofEntity>> dofsRowSequce;
  mutable boost::weak_ptr<std::vector<FENumeredDofEntity>> dofsColSequce;
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

        ordered_unique<tag<Unique_mi_tag>,
                       const_mem_fun<EntFiniteElement, UId,
                                     &EntFiniteElement::getGlobalUniqueId>>,

        ordered_non_unique<tag<Ent_mi_tag>,
                           const_mem_fun<EntFiniteElement, EntityHandle,
                                         &EntFiniteElement::getEnt>>,

        ordered_non_unique<
            tag<FiniteElement_name_mi_tag>,
            const_mem_fun<EntFiniteElement::interface_type_FiniteElement,
                          boost::string_ref, &EntFiniteElement::getNameRef>>,

        ordered_non_unique<
            tag<BitFEId_mi_tag>,
            const_mem_fun<EntFiniteElement::interface_type_FiniteElement,
                          BitFEId, &EntFiniteElement::getId>,
            LtBit<BitFEId>>,

        ordered_non_unique<
            tag<EntType_mi_tag>,
            const_mem_fun<EntFiniteElement::interface_type_RefEntity,
                          EntityType, &EntFiniteElement::getEntType>>,

        ordered_non_unique<
            tag<Composite_Name_And_Ent_mi_tag>,
            composite_key<
                EntFiniteElement,
                const_mem_fun<EntFiniteElement::interface_type_FiniteElement,
                              boost::string_ref, &EntFiniteElement::getNameRef>,
                const_mem_fun<EntFiniteElement, EntityHandle,
                              &EntFiniteElement::getEnt>>>

        >>
    EntFiniteElement_multiIndex;

/**
 *  \brief Entity finite element multi-index by finite element name
 *
 *  \ingroup fe_multi_indices
 */
typedef EntFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type
    EntFiniteElementByName;

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
                NumeredEntFiniteElement::interface_type_EntFiniteElement, UId,
                &NumeredEntFiniteElement::getGlobalUniqueId>>,
        ordered_non_unique<tag<Part_mi_tag>,
                           member<NumeredEntFiniteElement, unsigned int,
                                  &NumeredEntFiniteElement::part>>,
        ordered_non_unique<
            tag<FiniteElement_name_mi_tag>,
            const_mem_fun<NumeredEntFiniteElement::interface_type_FiniteElement,
                          boost::string_ref,
                          &NumeredEntFiniteElement::getNameRef>>,
        ordered_non_unique<
            tag<Ent_mi_tag>,
            const_mem_fun<
                NumeredEntFiniteElement::interface_type_EntFiniteElement,
                EntityHandle, &NumeredEntFiniteElement::getEnt>>,
        ordered_non_unique<
            tag<Composite_Name_And_Ent_mi_tag>,
            composite_key<
                NumeredEntFiniteElement,
                const_mem_fun<
                    NumeredEntFiniteElement::interface_type_FiniteElement,
                    boost::string_ref, &NumeredEntFiniteElement::getNameRef>,
                const_mem_fun<
                    NumeredEntFiniteElement::interface_type_EntFiniteElement,
                    EntityHandle, &NumeredEntFiniteElement::getEnt>>>,
        ordered_non_unique<
            tag<Composite_Name_And_Part_mi_tag>,
            composite_key<
                NumeredEntFiniteElement,
                const_mem_fun<
                    NumeredEntFiniteElement::interface_type_FiniteElement,
                    boost::string_ref, &NumeredEntFiniteElement::getNameRef>,
                member<NumeredEntFiniteElement, unsigned int,
                       &NumeredEntFiniteElement::part>>>>>
    NumeredEntFiniteElement_multiIndex;

/**
 *  \brief Entity finite element multi-index by finite element name
 *
 *  \ingroup fe_multi_indices
 */
typedef NumeredEntFiniteElement_multiIndex::index<
    FiniteElement_name_mi_tag>::type NumeredEntFiniteElementbyName;

/**
 *  \brief Entity finite element multi-index by finite element name and
 * partition
 *
 *  \ingroup fe_multi_indices
 */
typedef NumeredEntFiniteElement_multiIndex::index<
    Composite_Name_And_Part_mi_tag>::type NumeredEntFiniteElementbyNameAndPart;

/**
  @relates multi_index_container
  \brief MultiIndex for entities for FiniteElement
  \ingroup fe_multi_indices
 */
typedef multi_index_container<
    boost::shared_ptr<FiniteElement>,
    indexed_by<hashed_unique<tag<FiniteElement_Meshset_mi_tag>,
                             member<FiniteElement, EntityHandle,
                                    &FiniteElement::meshset>>,
               hashed_unique<
                   tag<BitFEId_mi_tag>,
                   const_mem_fun<FiniteElement, BitFEId, &FiniteElement::getId>,
                   HashBit<BitFEId>, EqBit<BitFEId>>,
               ordered_unique<tag<FiniteElement_name_mi_tag>,
                              const_mem_fun<FiniteElement, boost::string_ref,
                                            &FiniteElement::getNameRef>>>>
    FiniteElement_multiIndex;

// modificators

/**
 * \brief Change finite element part
 *
 * \ingroup fe_multi_indices
 */
struct NumeredEntFiniteElement_change_part {
  unsigned int pArt;
  NumeredEntFiniteElement_change_part(unsigned int part) : pArt(part){};
  void operator()(boost::shared_ptr<NumeredEntFiniteElement> &fe) {
    fe->part = pArt;
  }
  void operator()(NumeredEntFiniteElement &fe) { fe.part = pArt; }
};

/**
 * \brief Add field to column
 *
 * \ingroup fe_multi_indices
 */
struct FiniteElement_col_change_bit_add {
  BitFieldId fIdCol;
  FiniteElement_col_change_bit_add(const BitFieldId f_id_col)
      : fIdCol(f_id_col){};
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

/**
 * \brief Add field to row
 *
 * \ingroup fe_multi_indices
 */
struct FiniteElement_row_change_bit_add {
  BitFieldId fIdRow;
  FiniteElement_row_change_bit_add(const BitFieldId f_id_row)
      : fIdRow(f_id_row){};
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

/**
 * \brief Add field to data
 *
 * \ingroup fe_multi_indices
 */
struct FiniteElement_change_bit_add {
  BitFieldId fIdData;
  FiniteElement_change_bit_add(const BitFieldId f_id_data)
      : fIdData(f_id_data){};
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

/**
 * \brief Unset field from column
 *
 * \ingroup fe_multi_indices
 */
struct FiniteElement_col_change_bit_off {
  BitFieldId fIdCol;
  FiniteElement_col_change_bit_off(const BitFieldId f_id_col)
      : fIdCol(f_id_col){};
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

/**
 * \brief Unset field from row
 *
 * \ingroup fe_multi_indices
 */
struct FiniteElement_row_change_bit_off {
  BitFieldId fIdRow;
  FiniteElement_row_change_bit_off(const BitFieldId f_id_row)
      : fIdRow(f_id_row){};
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

/**
 * \brief Unset field from data
 *
 * \ingroup fe_multi_indices
 */
struct FiniteElement_change_bit_off {
  BitFieldId fIdData;
  FiniteElement_change_bit_off(const BitFieldId f_id_data)
      : fIdData(f_id_data){};
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

/**
 * \brief Reset field from column
 *
 * \ingroup fe_multi_indices
 */
struct FiniteElement_col_change_bit_reset {
  FiniteElement_col_change_bit_reset() = default;
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

/**
 * \brief Reset field from row
 *
 * \ingroup fe_multi_indices
 */
struct FiniteElement_row_change_bit_reset {
  FiniteElement_row_change_bit_reset() = default;
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

/**
 * \brief Reset field from data
 *
 * \ingroup fe_multi_indices
 */
struct FiniteElement_change_bit_reset {
  FiniteElement_change_bit_reset() = default;
  void operator()(boost::shared_ptr<FiniteElement> &fe);
};

} // namespace MoFEM

/**
 * Loop over DOFs in row on element
 * @param  FEPTR pointer to element structure \ref NumeredEntFiniteElement
 * @param  IT    iterator
 * @return       user return in for(_IT_FENUMEREDDOF_ROW_FOR_LOOP_(FEPTR,IT))
 * \ingroup fe_multi_indices
 */
#define _IT_FENUMEREDDOF_ROW_FOR_LOOP_(FEPTR, IT)                              \
  auto IT = FEPTR->getRowDofsPtr()->begin();                                   \
  IT != FEPTR->getRowDofsPtr()->end();                                         \
  IT++

/// \deprecated use _IT_FENUMEREDDOF_ROW_FOR_LOOP_
#define _IT_FENUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(FEPTR, IT)                   \
  _IT_FENUMEREDDOF_ROW_FOR_LOOP_(FEPTR, IT)

/**
 * Loop over DOFs in col on element
 * @param  FEPTR pointer to element structure \ref NumeredEntFiniteElement
 * @param  IT    iterator
 * @return       user return in for(_IT_FENUMEREDDOF_COL_FOR_LOOP_(FEPTR,IT))
 * \ingroup fe_multi_indices
 */
#define _IT_FENUMEREDDOF_COL_FOR_LOOP_(FEPTR, IT)                              \
  auto IT = FEPTR->getColDofsPtr()->begin();                                   \
  IT != FEPTR->getColDofsPtr()->end();                                         \
  IT++

/// \deprecated use _IT_FENUMEREDDOF_COL_FOR_LOOP_ instead
#define _IT_FENUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(FEPTR, IT)                   \
  _IT_FENUMEREDDOF_COL_FOR_LOOP_(FEPTR, IT)

/**
 * Loop over DOFs in row on element for particular filed
 * @param  FEPTR pointer to element structure \ref NumeredEntFiniteElement
 * @param  NAME  name of filed
 * @param  IT    iterator
 * @return       user return in
 * for(_IT_FENUMEREDDOF_BY_NAME_ROW_FOR_LOOP_(FEPTR,NAME,IT)) \ingroup
 * fe_multi_indices
 */
#define _IT_FENUMEREDDOF_BY_NAME_ROW_FOR_LOOP_(FEPTR, NAME, IT)                \
  auto IT = FEPTR->getRowDofsPtr()->get<FieldName_mi_tag>().lower_bound(NAME); \
  IT != FEPTR->getRowDofsPtr()->get<FieldName_mi_tag>().upper_bound(NAME);     \
  IT++

/// \deprecated use _IT_FENUMEREDDOF_BY_NAME_ROW_FOR_LOOP_ instead
#define _IT_FENUMEREDDOFMOFEMENTITY_BY_NAME_ROW_FOR_LOOP_(FEPTR, NAME, IT)     \
  _IT_FENUMEREDDOF_BY_NAME_ROW_FOR_LOOP_(FEPTR, NAME, IT)

/**
 * Loop over DOFs in col on element for particular filed
 * @param  FEPTR pointer to element structure \ref NumeredEntFiniteElement
 * @param  NAME  name of filed
 * @param  IT    iterator
 * @return       user return in
 * for(_IT_FENUMEREDDOF_BY_NAME_COL_FOR_LOOP_(FEPTR,NAME,IT)) \ingroup
 * fe_multi_indices
 */
#define _IT_FENUMEREDDOF_BY_NAME_COL_FOR_LOOP_(FEPTR, NAME, IT)                \
  auto IT = FEPTR->getColDofsPtr()->get<FieldName_mi_tag>().lower_bound(NAME); \
  IT != FEPTR->getColDofsPtr()->get<FieldName_mi_tag>().upper_bound(NAME);     \
  IT++

/// \deprecated use _IT_FENUMEREDDOF_BY_NAME_COL_FOR_LOOP_ instead
#define _IT_FENUMEREDDOFMOFEMENTITY_BY_NAME_COL_FOR_LOOP_(FEPTR, NAME, IT)     \
  _IT_FENUMEREDDOF_BY_NAME_COL_FOR_LOOP_(FEPTR, NAME, IT)

#endif // __FEMMULTIINDICES_HPP__

/***************************************************************************/ /**
                                                                               * \defgroup fe_multi_indices Finite elements structures and multi-indices
                                                                               * \ingroup mofem
                                                                               ******************************************************************************/
