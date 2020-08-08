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
template <typename FE, typename REFENT>
struct interface_FiniteElementImpl : public interface_RefElement<REFENT> {

  interface_FiniteElementImpl(const boost::shared_ptr<FE> fe_ptr,
                              const boost::shared_ptr<REFENT> ref_ents_ptr)
      : interface_RefElement<REFENT>(ref_ents_ptr){};

  virtual ~interface_FiniteElementImpl() = default;

  virtual boost::shared_ptr<const FiniteElement> &
  getFiniteElementPtr() const = 0;

  /**
   * @copydoc MoFEM::FiniteElement::getFEUId
   */
  inline const UId &getFEUId() const {
    return this->getFiniteElementPtr()->getFEUId();
  }

  /**
   * @copydoc MoFEM::FiniteElement::getId
   */
  inline BitFEId getId() const { return this->getFiniteElementPtr()->getId(); }

  /**
   * @copydoc MoFEM::FiniteElement::getMeshset
   */
  inline EntityHandle getMeshset() const {
    return this->getFiniteElementPtr()->getMeshset();
  }

  /**
   * @copydoc MoFEM::FiniteElement::getNameRef
   */
  inline boost::string_ref getNameRef() const {
    return this->getFiniteElementPtr()->getNameRef();
  }

  /**
   * @copydoc MoFEM::FiniteElement::getName
   */
  inline std::string getName() const {
    return this->getFiniteElementPtr()->getName();
  }

  /**
   * @copydoc MoFEM::FiniteElement::getBitFieldIdCol
   */
  inline BitFieldId getBitFieldIdCol() const {
    return this->getFiniteElementPtr()->getBitFieldIdCol();
  }

  /**
   * @copydoc MoFEM::FiniteElement::getBitFieldIdRow
   */
  inline BitFieldId getBitFieldIdRow() const {
    return this->getFiniteElementPtr()->getBitFieldIdRow();
  }

  /**
   * @copydoc MoFEM::FiniteElement::getBitFieldIdData
   */
  inline BitFieldId getBitFieldIdData() const {
    return this->getFiniteElementPtr()->getBitFieldIdData();
  }

  /**
   * @copydoc MoFEM::FiniteElement::getBitNumber
   */
  inline unsigned int getBitNumber() const {
    return this->getFiniteElementPtr()->getBitNumber();
  }
};

template <typename FE, typename REFENT>
struct interface_FiniteElement : public interface_RefElement<REFENT> {
  using interface_FiniteElementImpl<FE, REFENT>::interface_FiniteElementImpl;
  virtual ~interface_FiniteElement() = default;
};

template <typename T>
struct interface_FiniteElement<T, T>
    : public interface_FiniteElementImpl<T, T> {

  interface_FiniteElement(const boost::shared_ptr<T> fe_ptr,
                          const boost::shared_ptr<T> ref_ents_ptr)
      : interface_FiniteElementImpl<T, T>(fe_ptr, ref_ents_ptr),
        sFiniteElementPtr(fe_ptr) {}

  inline boost::shared_ptr<const FiniteElement> &getFiniteElementPtr() const {
    return this->sFiniteElementPtr->getFiniteElementPtr();
  }

  mutable boost::shared_ptr<T> sFiniteElementPtr;
};

/**
 * \brief Finite element data for entity
 * \ingroup fe_multi_indices
 */
struct EntFiniteElement
    : public interface_FiniteElementImpl<FiniteElement, RefElement> {

  using interface_type_RefEntity = interface_RefEntity<RefElement>;
  using interface_type_RefElement = interface_RefElement<RefElement>;
  using interface_type_FiniteElement =
      interface_FiniteElementImpl<FiniteElement, RefElement>;

  EntFiniteElement(const boost::shared_ptr<RefElement> &ref_finite_element,
                   const boost::shared_ptr<FiniteElement> &fe_ptr);
  virtual ~EntFiniteElement();

  virtual boost::shared_ptr<const FiniteElement> &getFiniteElementPtr() const {
    return finiteElementPtr;
  }

  /**
   * \brief Get unique UId for finite element entity
   * @return UId
   */
  inline UId getLocalUniqueId() const { return getLocalUniqueIdCalculate(); }

  static inline UId getLocalUniqueIdCalculate(const EntityHandle ent,
                                              UId fe_uid) {
    return fe_uid |= ent;
  }

  /**
   * \brief Generate UId for finite element entity
   * @return finite element entity unique Id
   */
  inline UId getLocalUniqueIdCalculate() const {
    return getLocalUniqueIdCalculate(getEnt(), getFEUId());
  }

  /**
   * \brief Get data data dos multi-index structure
   * @return Reference multi-index FEDofEntity_multiIndex
   */
  inline const FEDofEntity_multiIndex &
  getDataDofs(const DofEntity_multiIndex &dofs_field) const {
    RefEntityTmp<0>::refElementPtr = this->getRefElement();
    if (lastSeenEntFiniteElement != this) {
      dataDofs = boost::make_shared<FEDofEntity_multiIndex>();
      if (getDofView(getDataFieldEnts(), dofs_field, *dataDofs,
                     moab::Interface::UNION))
        THROW_MESSAGE("dataDofs can not be created");
      lastSeenEntFiniteElement = this;
    }
    return *dataDofs;
  };

  inline boost::shared_ptr<FEDofEntity_multiIndex> &
  getDataDofsPtr(const DofEntity_multiIndex &dofs_field) {
    RefEntityTmp<0>::refElementPtr = this->getRefElement();
    if (lastSeenEntFiniteElement != this) {
      dataDofs = boost::make_shared<FEDofEntity_multiIndex>();
      if (getDofView(getDataFieldEnts(), dofs_field, *dataDofs,
                     moab::Interface::UNION))
        THROW_MESSAGE("dataDofs can not be created");
      lastSeenEntFiniteElement = this;
    }
    return dataDofs;
  };

  inline const FieldEntity_vector_view &getDataFieldEnts() const {
    return *dataFieldEnts;
  };

  inline boost::shared_ptr<FieldEntity_vector_view> &getDataFieldEntsPtr() {
    return dataFieldEnts;
  };

  inline const FieldEntity_vector_view &getRowFieldEnts() const {
    return *rowFieldEnts;
  };

  inline boost::shared_ptr<FieldEntity_vector_view> &getRowFieldEntsPtr() {
    return rowFieldEnts;
  };

  inline const FieldEntity_vector_view &getColFieldEnts() const {
    return *colFieldEnts;
  };

  inline boost::shared_ptr<FieldEntity_vector_view> &getColFieldEntsPtr() {
    return colFieldEnts;
  };

  friend std::ostream &operator<<(std::ostream &os, const EntFiniteElement &e);

  template <typename FE_ENTS, typename MOFEM_DOFS, typename MOFEM_DOFS_VIEW>
  inline MoFEMErrorCode
  getDofView(const FE_ENTS &fe_ents_view, const MOFEM_DOFS &mofem_dofs,
             MOFEM_DOFS_VIEW &dofs_view, const int operation_type) const {
    MoFEMFunctionBeginHot;

    if (operation_type == moab::Interface::UNION) {

      auto hint = dofs_view.begin();
      using ValType = typename std::remove_reference<decltype(**hint)>::type;

      for (auto &it : fe_ents_view) {
        if (auto e = it.lock()) {
          const auto &uid = e->getLocalUniqueId();
          auto dit = mofem_dofs.lower_bound(uid);
          if (dit != mofem_dofs.end()) {
            const auto hi_dit = mofem_dofs.upper_bound(
                uid | static_cast<UId>(MAX_DOFS_ON_ENTITY - 1));
            for (; dit != hi_dit; ++dit)
              hint = dofs_view.emplace_hint(
                  hint, boost::reinterpret_pointer_cast<ValType>(*dit));
          }
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
    return getDofView(getRowFieldEnts(), mofem_dofs, dofs_view, operation_type);
  }

  template <typename MOFEM_DOFS, typename MOFEM_DOFS_VIEW>
  inline MoFEMErrorCode
  getColDofView(const MOFEM_DOFS &mofem_dofs, MOFEM_DOFS_VIEW &dofs_view,
                const int operation_type = moab::Interface::UNION) {
    return getDofView(getColFieldEnts(), mofem_dofs, dofs_view, operation_type);
  }

  MoFEMErrorCode getElementAdjacency(const boost::shared_ptr<Field> field_ptr,
                                     Range &adjacency);

private:
  mutable boost::shared_ptr<const FiniteElement> finiteElementPtr;
  mutable boost::shared_ptr<FieldEntity_vector_view> dataFieldEnts;
  mutable boost::shared_ptr<FieldEntity_vector_view> rowFieldEnts;
  mutable boost::shared_ptr<FieldEntity_vector_view> colFieldEnts;

  static const EntFiniteElement *lastSeenEntFiniteElement;
  static boost::shared_ptr<FEDofEntity_multiIndex> dataDofs;
};

/**
 * \brief interface to EntFiniteElement
 * \ingroup fe_multi_indices
 */
template <typename T>
struct interface_EntFiniteElement : public interface_FiniteElement<T, T> {

  interface_EntFiniteElement(const boost::shared_ptr<T> &sptr)
      : interface_FiniteElement<T, T>(sptr, sptr) {}
  virtual ~interface_EntFiniteElement() = default;

  inline const FEDofEntity_multiIndex &
  getDataDofs(const DofEntity_multiIndex &dofs_field) const {
    return this->sPtr->getDataDofs(dofs_field);
  }

  inline boost::shared_ptr<FEDofEntity_multiIndex> &
  getDataDofsPtr(const DofEntity_multiIndex &dofs_field) {
    return this->sPtr->getDataDofsPtr(dofs_field);
  }

  inline const FieldEntity_vector_view &getDataFieldEnts() const {
    return this->sPtr->getDataFieldEnts();
  };

  inline boost::shared_ptr<FieldEntity_vector_view> &getDataFieldEntsPtr() {
    return this->sPtr->getDataFieldEntsPtr();
  };

  inline const FieldEntity_vector_view &getRowFieldEnts() const {
    return this->sPtr->getRowFieldEnts();
  };

  inline boost::shared_ptr<FieldEntity_vector_view> &getRowFieldEntsPtr() {
    return this->sPtr->getRowFieldEntsPtr();
  }

  inline const FieldEntity_vector_view &getColFieldEnts() const {
    return this->sPtr->getColFieldEnts();
  };

  inline boost::shared_ptr<FieldEntity_vector_view> &getColFieldEntsPtr() {
    return this->sPtr->getColFieldEntsPtr();
  };

  /**
   * \brief Get unique UId for finite element entity
   * @return UId
   */
  inline UId getLocalUniqueId() const { return this->sPtr->getLocalUniqueId(); }

  SideNumber_multiIndex &getSideNumberTable() const {
    return this->sPtr->getSideNumberTable();
  }

  inline MoFEMErrorCode getElementAdjacency(const Field *field_ptr,
                                            Range &adjacency) {
    return this->getElementAdjacency(field_ptr, adjacency);
  }

  inline const boost::shared_ptr<RefElement> &getRefElement() const {
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

  using interface_type_FiniteElement =
      interface_FiniteElementImpl<EntFiniteElement, EntFiniteElement>;
  using interface_type_EntFiniteElement =
      interface_EntFiniteElement<EntFiniteElement>;

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

  friend std::ostream &operator<<(std::ostream &os,
                                  const NumeredEntFiniteElement &e);

private:
  boost::shared_ptr<FENumeredDofEntity_multiIndex>
      rowDofs; ///< indexed dofs on rows
  boost::shared_ptr<FENumeredDofEntity_multiIndex>
      colDofs; ///< indexed dofs on columns
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
                                     &EntFiniteElement::getLocalUniqueId>>,

        ordered_non_unique<
            tag<Ent_mi_tag>,
            const_mem_fun<EntFiniteElement::interface_type_RefEntity,
                          EntityHandle, &EntFiniteElement::getEnt>>,

        ordered_non_unique<
            tag<FiniteElement_name_mi_tag>,
            const_mem_fun<EntFiniteElement::interface_type_FiniteElement,
                          boost::string_ref, &EntFiniteElement::getNameRef>>,

        ordered_non_unique<
            tag<Composite_Name_And_Ent_mi_tag>,
            composite_key<
                EntFiniteElement,
                const_mem_fun<EntFiniteElement::interface_type_FiniteElement,
                              boost::string_ref, &EntFiniteElement::getNameRef>,
                const_mem_fun<EntFiniteElement::interface_type_RefEntity,
                              EntityHandle, &EntFiniteElement::getEnt>>>

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
                &NumeredEntFiniteElement::getLocalUniqueId>>,
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
            const_mem_fun<NumeredEntFiniteElement::interface_type_RefEntity,
                          EntityHandle, &NumeredEntFiniteElement::getEnt>>,
        ordered_non_unique<
            tag<Composite_Name_And_Ent_mi_tag>,
            composite_key<
                NumeredEntFiniteElement,
                const_mem_fun<
                    NumeredEntFiniteElement::interface_type_FiniteElement,
                    boost::string_ref, &NumeredEntFiniteElement::getNameRef>,
                const_mem_fun<NumeredEntFiniteElement::interface_type_RefEntity,
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

#endif // __FEMMULTIINDICES_HPP__

/**
 * \defgroup fe_multi_indices Finite elements structures and multi-indices
 * \ingroup mofem
 **/
