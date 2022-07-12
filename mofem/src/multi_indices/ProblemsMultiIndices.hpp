/** \file ProblemsMultiIndices.hpp
 * \brief Multi-index containers, data structures for problems and other
 * low-level functions
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef __PROBLEMSMULTIINDICES_HPP__
#define __PROBLEMSMULTIINDICES_HPP__

namespace MoFEM {

struct Problem;

/**
 * Data structure created when composite problem is created
 */
struct ComposedProblemsData {

  std::vector<const Problem *> rowProblemsAdd;
  std::vector<const Problem *> colProblemsAdd;

  std::vector<IS> rowIs;
  std::vector<IS> colIs;

  inline MoFEMErrorCode getRowIs(IS *is, const unsigned int pp) const;

  inline MoFEMErrorCode getColIs(IS *is, const unsigned int pp) const;

  virtual ~ComposedProblemsData();
};

/** \brief keeps basic data about problem
 * \ingroup problems_multi_indices
 *
 * This is low level structure with information about problem, what elements
 * compose problem and what DOFs are on rows and columns.
 *
 * \todo fix names following name convention
 *
 */
struct Problem {

  EntityHandle meshset; ///< Problem meshset (on tag of this meshset all data
                        ///< related to problem are stored)
  BitProblemId *tagId;  ///< Unique problem ID
  const char *tagName;  ///< Problem name
  int tagNameSize;      ///< Size of problem name
  BitFEId *tagBitFEId;  ///< IDs of finite elements in problem
  BitRefLevel *tagBitRefLevel; ///< BitRef level of finite elements in problem
  BitRefLevel *tagBitRefLevelMask; ///< BItRefMask of elements in problem

  mutable DofIdx nbDofsRow;      ///< Global number of DOFs in  row
  mutable DofIdx nbDofsCol;      ///< Global number of DOFs in col
  mutable DofIdx nbLocDofsRow;   ///< Local number of DOFs in row
  mutable DofIdx nbLocDofsCol;   ///< Local number of DOFs in colIs
  mutable DofIdx nbGhostDofsRow; ///< Number of ghost DOFs in row
  mutable DofIdx nbGhostDofsCol; ///< Number of ghost DOFs in col

  mutable boost::shared_ptr<NumeredDofEntity_multiIndex>
      numeredRowDofsPtr; ///< store DOFs on rows for this problem
  mutable boost::shared_ptr<NumeredDofEntity_multiIndex>
      numeredColDofsPtr; ///< store DOFs on columns for this problem
  mutable boost::shared_ptr<NumeredEntFiniteElement_multiIndex>
      numeredFiniteElementsPtr; ///< store finite elements

  /**
   * \brief get access to numeredRowDofsPtr storing DOFs on rows
   */
  inline auto &getNumeredRowDofsPtr() const { return numeredRowDofsPtr; }

  /**
   * \brief get access to numeredColDofsPtr storing DOFs on cols
   */
  inline auto &getNumeredColDofsPtr() const { return numeredColDofsPtr; }

  /**
   * \brief get access to reference for multi-index storing finite elements
   */
  inline const auto &getNumeredFiniteElementsPtr() const {
    return numeredFiniteElementsPtr;
  }

  /**
   * \brief Subproblem problem data
   */
  struct SubProblemData;

  /**
   * Pointer to data structure. This pointer has allocated data only for
   * sub problems.
   */
  mutable boost::shared_ptr<SubProblemData> subProblemData;

  /**
   * \brief Get main problem of sub-problem is
   * @return    sub problem data structure
   */
  inline boost::shared_ptr<SubProblemData> &getSubData() const {
    return subProblemData;
  }

  /**
   * Pointer to data structure from which this problem is composed
   */
  mutable boost::shared_ptr<ComposedProblemsData> composedProblemsData;

  /**
   * \brief Het composed problems data structure
   */
  inline auto &getComposedProblemsData() const { return composedProblemsData; }

  /**
   * \brief get DOFs from problem
   *
   * Note that \e ent_dof_idx is not coefficient number, is local number of DOFs
   * on the entity. The coefficient number and local index of DOFs or entity are
   * the same on vertices and H1 approximation.
   *
   * @param  field_bit_number       field name field_bit_number = (use
   * m_field.get_field_bit_number(field_name);
   * @param  ent        entity handle
   * @param  ent_dof_idx index of DOFs on entity
   * @param  row_or_col ROW or COL
   * @param  dof_ptr    shared pointer to DOFs if found
   * @return            error code
   */
  MoFEMErrorCode getDofByNameEntAndEntDofIdx(
      const int field_bit_number, const EntityHandle ent, const int ent_dof_idx,
      const RowColData row_or_col,
      boost::shared_ptr<NumeredDofEntity> &dof_ptr) const;

/**
 * use with loops to iterate problem fes
 * \ingroup problems_multi_indices
 *
 * for(_IT_NUMEREDFE_BY_NAME_FOR_LOOP_(PROBLEMPTR,NAME,IT)) {
 *   ...
 * }
 *
 */
#define _IT_NUMEREDFE_BY_NAME_FOR_LOOP_(PROBLEMPTR, NAME, IT)                  \
  auto IT = PROBLEMPTR->getNumeredFEsBegin(NAME);                              \
  IT != PROBLEMPTR->getNumeredFEsEnd(NAME);                                    \
  IT++

  auto getNumeredFEsBegin(std::string fe_name) const {
    return numeredFiniteElementsPtr->get<FiniteElement_name_mi_tag>()
        .lower_bound(fe_name);
  }

  auto getNumeredFEsEnd(std::string fe_name) const {
    return numeredFiniteElementsPtr->get<FiniteElement_name_mi_tag>()
        .upper_bound(fe_name);
  }

/**
 * \brief use with loops to iterate row DOFs
 * \ingroup problems_multi_indices
 *
 * \code
 * for(_IT_NUMEREDDOF_ROW_FOR_LOOP_(PROBLEMPTR,IT)) {
 *   ...
 * }
 * \endcode
 *
 */
#define _IT_NUMEREDDOF_ROW_FOR_LOOP_(PROBLEMPTR, IT)                           \
  NumeredDofEntity_multiIndex::iterator IT =                                   \
      PROBLEMPTR->getNumeredRowDofsBegin();                                    \
  IT != PROBLEMPTR->getNumeredRowDofsEnd();                                    \
  IT++

/**
 * use with loops to iterate col DOFs
 * \ingroup problems_multi_indices
 *
 * \code
 * for(_IT_NUMEREDDOF_COL_FOR_LOOP_(PROBLEMPTR,IT)) {
 *   ...
 * }
 * \endcode
 *
 */
#define _IT_NUMEREDDOF_COL_FOR_LOOP_(PROBLEMPTR, IT)                           \
  NumeredDofEntity_multiIndex::iterator IT =                                   \
      PROBLEMPTR->getNumeredColDofsBegin();                                    \
  IT != PROBLEMPTR->getNumeredColDofsEnd();                                    \
  IT++

  /// get begin iterator for numeredRowDofsPtr (insted you can use
  /// #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredRowDofsBegin() const {
    return numeredRowDofsPtr->begin();
  }

  /// get end iterator for numeredRowDofsPtr (insted you can use
  /// #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredRowDofsEnd() const {
    return numeredRowDofsPtr->end();
  }

  /// get begin iterator for numeredColDofsPtr (insted you can use
  /// #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredColDofsBegin() const {
    return numeredColDofsPtr->begin();
  }

  /// get end iterator for numeredColDofsPtr (insted you can use
  /// #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredColDofsEnd() const {
    return numeredColDofsPtr->end();
  }

/**
 * \brief use with loops to iterate row DOFs
 * \ingroup problems_multi_indices
 *
 * \code
 * for(_IT_NUMEREDDOF_ROW_BY_LOCIDX_FOR_LOOP_(PROBLEMPTR,IT)) {
 *   ...
 * }
 * \endcode
 *
 */
#define _IT_NUMEREDDOF_ROW_BY_LOCIDX_FOR_LOOP_(PROBLEMPTR, IT)                 \
  NumeredDofEntityByLocalIdx::iterator IT =                                    \
      PROBLEMPTR->getNumeredRowDofsByLocIdxBegin(0);                           \
  IT != PROBLEMPTR->getNumeredRowDofsByLocIdxEnd(                              \
            PROBLEMPTR->getNbLocalDofsRow() - 1);                              \
  IT++

/**
 * \brief use with loops to iterate col DOFs
 *
 * \code
 * for(_IT_NUMEREDDOF_COL_BY_LOCIDX_FOR_LOOP_(PROBLEMPTR,IT)) {
 *   ...
 * }
 * \endcode
 *
 */
#define _IT_NUMEREDDOF_COL_BY_LOCIDX_FOR_LOOP_(PROBLEMPTR, IT)                 \
  NumeredDofEntityByUId::iterator IT =                                         \
      PROBLEMPTR->getNumeredColDofsByLocIdxBegin(0);                           \
  IT != PROBLEMPTR->getNumeredColDofsByLocIdxEnd(                              \
            PROBLEMPTR->getNbLocalDofsRow() - 1);                              \
  IT++

  /// get begin iterator for numeredRowDofsPtr (insted you can use
  /// #_IT_NUMEREDDOF_ROW_FOR_LOOP_ for loops)
  auto getNumeredRowDofsByLocIdxBegin(const DofIdx locidx) const {
    return numeredRowDofsPtr->get<PetscLocalIdx_mi_tag>().lower_bound(locidx);
  }

  /// get end iterator for numeredRowDofsPtr (insted you can use
  /// #_IT_NUMEREDDOF_ROW_FOR_LOOP_ for loops)
  auto getNumeredRowDofsByLocIdxEnd(const DofIdx locidx) const {
    return numeredRowDofsPtr->get<PetscLocalIdx_mi_tag>().upper_bound(locidx);
  }

  /// get begin iterator for numeredColDofsPtr (insted you can use
  /// #_IT_NUMEREDDOF_COL_FOR_LOOP_ for loops)
  auto getNumeredColDofsByLocIdxBegin(const DofIdx locidx) const {
    return numeredColDofsPtr->get<PetscLocalIdx_mi_tag>().lower_bound(locidx);
  }

  /// get end iterator for numeredColDofsPtr (insted you can use
  /// #_IT_NUMEREDDOF_COL_FOR_LOOP_ for loops)
  auto getNumeredColDofsByLocIdxEnd(const DofIdx locidx) const {
    return numeredColDofsPtr->get<PetscLocalIdx_mi_tag>().upper_bound(locidx);
  }

/**
 * \brief use with loops to iterate row DOFs
 * \ingroup problems_multi_indices
 *
 * \code
 * for(_IT_NUMEREDDOF_BY_ENT_ROW_FOR_LOOP_(PROBLEMPTR,ENT,IT)) {
 *   ...
 * }
 * \endcode
 *
 */
#define _IT_NUMEREDDOF_ROW_BY_ENT_FOR_LOOP_(PROBLEMPTR, ENT, IT)               \
  auto IT = PROBLEMPTR->getNumeredRowDofsByEntBegin(ENT);                      \
  IT != PROBLEMPTR->getNumeredRowDofsByEntEnd(ENT);                            \
  IT++

/**
 * \brief use with loops to iterate col DOFs
 * \ingroup problems_multi_indices
 *
 * \code
 * for(_IT_NUMEREDDOF_COL_BY_ENT_FOR_LOOP_(PROBLEMPTR,ENT,IT)) {
 *   ...
 * }
 * \endcode
 *
 */
#define _IT_NUMEREDDOF_COL_BY_ENT_FOR_LOOP_(PROBLEMPTR, ENT, IT)               \
  auto IT = PROBLEMPTR->getNumeredColDofsByEntBegin(ENT);                      \
  IT != PROBLEMPTR->getNumeredColDofsByEntEnd(ENT);                            \
  IT++

  /// get begin iterator for numeredRowDofsPtr (insted you can use
  /// #_IT_NUMEREDDOF_ROW_BY_ENT_FOR_LOOP_ for loops)
  auto getNumeredRowDofsByEntBegin(const EntityHandle ent) const {
    return numeredRowDofsPtr->get<Ent_mi_tag>().lower_bound(ent);
  }

  /// get end iterator for numeredRowDofsPtr (insted you can use
  /// #_IT_NUMEREDDOF_ROW_BY_ENT_FOR_LOOP_ for loops)
  auto getNumeredRowDofsByEntEnd(const EntityHandle ent) const {
    return numeredRowDofsPtr->get<Ent_mi_tag>().upper_bound(ent);
  }

  /// get begin iterator for numeredColDofsPtr (insted you can use
  /// #_IT_NUMEREDDOF_COL_BY_ENT_FOR_LOOP_ for loops)
  auto getNumeredColDofsByEntBegin(const EntityHandle ent) const {
    return numeredColDofsPtr->get<Ent_mi_tag>().lower_bound(ent);
  }

  /// get end iterator for numeredColDofsPtr (insted you can use
  /// #_IT_NUMEREDDOF_COL_BY_ENT_FOR_LOOP_ for loops)
  auto getNumeredColDofsByEntEnd(const EntityHandle ent) const {
    return numeredColDofsPtr->get<Ent_mi_tag>().upper_bound(ent);
  }

/**
 * use with loops to iterate row DOFs
 * \ingroup problems_multi_indices
 *
 * \code
 * for(_IT_NUMEREDDOF_BY_NAME_ROW_FOR_LOOP_(PROBLEMPTR,m_field.get_field_bit_number(FIELD_BIT_NUMBER),IT))
 * {
 *   ...
 * }
 * \endcode
 *
 */
#define _IT_NUMEREDDOF_ROW_BY_BITNUMBER_FOR_LOOP_(PROBLEMPTR,                  \
                                                  FIELD_BIT_NUMBER, IT)        \
  auto IT = PROBLEMPTR->numeredRowDofsPtr->lower_bound(                        \
      FieldEntity::getLoBitNumberUId(FIELD_BIT_NUMBER));                       \
  IT != PROBLEMPTR->numeredRowDofsPtr->upper_bound(                            \
            FieldEntity::getHiBitNumberUId(FIELD_BIT_NUMBER));                 \
  IT++

/**
 * \brief use with loops to iterate col DOFs
 * \ingroup problems_multi_indices
 *
 * \code
 * for(_IT_NUMEREDDOF_COL_BY_BITNUMBER_FOR_LOOP_(PROBLEMPTR,m_field.get_field_bit_number(FIELD_BIT_NUMBER),IT))
 * {
 *   ...
 * }
 * \endcode
 *
 */
#define _IT_NUMEREDDOF_COL_BY_BITNUMBER_FOR_LOOP_(PROBLEMPTR,                  \
                                                  FIELD_BIT_NUMBER, IT)        \
  auto IT = PROBLEMPTR->numeredColDofsPtr->lower_bound(                        \
      FieldEntity::getLoBitNumberUId(FIELD_BIT_NUMBER));                       \
  IT != PROBLEMPTR->numeredColDofsPtr->upper_bound(                            \
            FieldEntity::getHiBitNumberUId(FIELD_BIT_NUMBER));                 \
  IT++

  Problem(Interface &moab, const EntityHandle meshset);

  virtual ~Problem() = default;

  inline BitProblemId getId() const { return *((BitProblemId *)tagId); }

  inline auto getName() const {
    return std::string((char *)tagName, tagNameSize);
  }

  inline DofIdx getNbDofsRow() const { return nbDofsRow; }
  inline DofIdx getNbDofsCol() const { return nbDofsCol; }
  inline DofIdx getNbLocalDofsRow() const { return nbLocDofsRow; }
  inline DofIdx getNbLocalDofsCol() const { return nbLocDofsCol; }
  inline DofIdx getNbGhostDofsRow() const { return nbGhostDofsRow; }
  inline DofIdx getNbGhostDofsCol() const { return nbGhostDofsCol; }

  inline BitRefLevel getBitRefLevel() const { return *tagBitRefLevel; }
  inline BitRefLevel getBitRefLevelMask() const { return *tagBitRefLevelMask; }

  // DEPRECATED inline BitRefLevel getMaskBitRefLevel() const {
  //   return *tagBitRefLevelMask;
  // }

  /**
   * @brief Get the Row Dofs By Petsc Global Dof Idx object
   *
   * @param idx
   * @param dof_ptr
   * @param bh
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode getRowDofsByPetscGlobalDofIdx(DofIdx idx,
                                               const NumeredDofEntity **dof_ptr,
                                               MoFEMTypes bh = MF_EXIST) const;

  /**
   * @brief Get the Col Dofs By Petsc Global Dof Idx object
   *
   * @param idx
   * @param dof_ptr
   * @param bh
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode getColDofsByPetscGlobalDofIdx(DofIdx idx,
                                               const NumeredDofEntity **dof_ptr,
                                               MoFEMTypes bh = MF_EXIST) const;

  /**
   * @brief Get the Row Dofs By Petsc Global Dof Idx object
   *
   * @param idx
   * @return boost::weak_ptr<NumeredDofEntity>
   */
  boost::weak_ptr<NumeredDofEntity>
  getRowDofsByPetscGlobalDofIdx(DofIdx idx) const;

  /**
   * @brief Get the Col Dofs By Petsc Global Dof Idx object
   *
   * @param idx
   * @return boost::weak_ptr<NumeredDofEntity>
   */
  boost::weak_ptr<NumeredDofEntity>
  getColDofsByPetscGlobalDofIdx(DofIdx idx) const;

  BitFEId getBitFEId() const;

  friend std::ostream &operator<<(std::ostream &os, const Problem &e);

  /**
   * \brief Get number of finite elements by name on processors
   *
   * Size retuned IS is equal to size of processors.
   *
   * What PetscLayout see,
   <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/IS/index.html>
   *
   * Example of usage for layout
   * \code

   PetscInt rstart, rend;
   ierr = PetscLayoutGetRange(layout, &rstart, &rend); CHKERRG(ierr);
   int global_size;
   ierr = PetscLayoutGetSize(layout,&global_size); CHKERRG(ierr);

   \endcode
   *
   * @param  comm Communicator
   * @param  name Finite element name
   * @param  layout Get number of elements on each processor
   * @return      error code
   */
  MoFEMErrorCode getNumberOfElementsByNameAndPart(MPI_Comm comm,
                                                  const std::string name,
                                                  PetscLayout *layout) const;

  /**
   * \brief Get number of finite elements on processors
   *
   * Size retuned IS is equal to size of processors.
   *
   * What PetscLayout see,
   <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/IS/index.html>
   *
   * Example of usage for layout
   * \code

   PetscInt rstart, rend;
   ierr = PetscLayoutGetRange(layout, &rstart, &rend); CHKERRG(ierr);
   int global_size;
   ierr = PetscLayoutGetSize(layout,&global_size); CHKERRG(ierr);

   \endcode
   *
   * @param  comm Communicator
   * @param  layout Get number of elements on each processor
   * @return      error code
   */
  MoFEMErrorCode getNumberOfElementsByPart(MPI_Comm comm,
                                           PetscLayout *layout) const;

  typedef multi_index_container<boost::weak_ptr<std::vector<NumeredDofEntity>>,
                                indexed_by<sequenced<>>>
      SequenceDofContainer;

  /**
   * \brief Get reference to sequence data numbered dof container
   *
   * In sequence data container data are physically stored. The purpose of this
   * is to allocate NumeredDofEntity data in bulk, having only one allocation
   instead
   * each time entity is inserted. That makes code efficient.
   *
   * The vector in sequence is destroyed if last entity inside that vector is
   * destroyed. All MoFEM::NumeredDofEntity have aliased shared_ptr which points
   to the vector.

   * @return MoFEM::Problem::SequenceDofContainer
   */
  inline auto &getRowDofsSequence() const { return sequenceRowDofContainer; }

  /**
   * \brief Get reference to sequence data numbered dof container
   *
   * In sequence data container data are physically stored. The purpose of this
   * is to allocate NumeredDofEntity data in bulk, having only one allocation
   instead
   * each time entity is inserted. That makes code efficient.
   *
   * The vector in sequence is destroyed if last entity inside that vector is
   * destroyed. All MoFEM::NumeredDofEntity have aliased shared_ptr which points
   to the vector.

   * @return MoFEM::Problem::SequenceDofContainer
   */
  inline auto &getColDofsSequence() const { return sequenceColDofContainer; }

  using EmptyFieldBlocks = std::pair<BitFieldId, BitFieldId>;

  /**
   * @brief Get the empty field blocks
   *
   * Emtpy field blocks is a pair contains IDs of the fields for which matrix
   * has zero entries.
   *
   * @return EmptyFieldBlocks&
   */
  inline EmptyFieldBlocks &getEmptyFieldBlocks() const {
    return emptyFieldBlocks;
  }

  /**
   * @brief Add fields to the empty field blocks
   *
   * Emtpy field blocks is a pair contains IDs of the fields for which matrix
   * has zero entries.
   *
   * @param add_fields
   * @return EmptyFieldBlocks&
   */
  inline EmptyFieldBlocks &
  addFieldToEmptyFieldBlocks(const EmptyFieldBlocks add_fields) const {
    emptyFieldBlocks.first |= add_fields.first;
    emptyFieldBlocks.second |= add_fields.second;
    return emptyFieldBlocks;
  }

private:
  // Keep vector of DoFS on entity
  mutable boost::shared_ptr<SequenceDofContainer> sequenceRowDofContainer;
  mutable boost::shared_ptr<SequenceDofContainer> sequenceColDofContainer;

  mutable EmptyFieldBlocks emptyFieldBlocks;
};

using EmptyFieldBlocks = Problem::EmptyFieldBlocks;

/**
 * \brief Subproblem problem data
 */
struct Problem::SubProblemData {

  SmartPetscObj<IS>
      rowIs; ///< indices of main problem of which sub problem is this
  SmartPetscObj<IS>
      colIs; ///< indices of main problem of which sub problem is this
  SmartPetscObj<AO>
      rowMap; ///< mapping form main problem indices to sub-problem indices
  SmartPetscObj<AO>
      colMap; ///< mapping form main problem indices to sub-problem indices

  inline auto getSmartRowIs() { return rowIs; }
  inline auto getSmartColIs() { return colIs; }
  inline auto getSmartRowMap() { return rowMap; }
  inline auto getSmartColMap() { return colMap; }

  /**
   * get row Is for sub problem
   * @param  is create is
   * @return    error code
   */
  inline MoFEMErrorCode getRowIs(IS *is) {
    MoFEMFunctionBeginHot;
    PetscObjectReference((PetscObject)rowIs);
    *is = rowIs;
    MoFEMFunctionReturnHot(0);
  }

  /**
   * get col Is for sub problem
   * @param  is create is
   * @return    error code
   */
  inline MoFEMErrorCode getColIs(IS *is) {
    MoFEMFunctionBeginHot;
    PetscObjectReference((PetscObject)colIs);
    *is = colIs;
    MoFEMFunctionReturnHot(0);
  };

  /**
   * get row AO mapping for sub problem
   * @param  ao get mapping
   * @return    error code
   */
  inline MoFEMErrorCode getRowMap(AO *ao) {
    MoFEMFunctionBeginHot;
    PetscObjectReference((PetscObject)rowMap);
    *ao = rowMap;
    MoFEMFunctionReturnHot(0);
  }

  /**
   * get col AO mapping for sub problem
   * @param  ao get mapping
   * @return    error code
   */
  inline MoFEMErrorCode getColMap(AO *ao) {
    MoFEMFunctionBeginHot;
    PetscObjectReference((PetscObject)colMap);
    *ao = colMap;
    MoFEMFunctionReturnHot(0);
  }

  SubProblemData() = default;
  virtual ~SubProblemData() = default;
};

using EmptyFieldBlocks = Problem::EmptyFieldBlocks;

/**
 * @relates multi_index_container
 * \brief MultiIndex for entities for Problem
 * \ingroup fe_multi_indices
 */
typedef multi_index_container<
    Problem,
    indexed_by<
        ordered_unique<tag<Meshset_mi_tag>,
                       member<Problem, EntityHandle, &Problem::meshset>>,
        hashed_unique<tag<BitProblemId_mi_tag>,
                      const_mem_fun<Problem, BitProblemId, &Problem::getId>,
                      HashBit<BitProblemId>, EqBit<BitProblemId>>,
        hashed_unique<tag<Problem_mi_tag>,
                      const_mem_fun<Problem, std::string, &Problem::getName>>>>
    Problem_multiIndex;

/** \brief add ref level to problem
 * \ingroup problems_multi_indices
 */
struct ProblemChangeRefLevelBitAdd {
  BitRefLevel bit;
  ProblemChangeRefLevelBitAdd(const BitRefLevel _bit) : bit(_bit){};
  void operator()(Problem &p) { *(p.tagBitRefLevel) |= bit; };
};

/** \brief set ref level to problem
 * \ingroup problems_multi_indices
 */
struct ProblemChangeRefLevelBitSet {
  BitRefLevel bit;
  ProblemChangeRefLevelBitSet(const BitRefLevel _bit) : bit(_bit){};
  void operator()(Problem &p) { *(p.tagBitRefLevel) = bit; };
};

/** \brief set prof dof bit ref mask
 * \ingroup problems_multi_indices
 */
struct ProblemChangeRefLevelBitDofMaskSet {
  BitRefLevel bit;
  ProblemChangeRefLevelBitDofMaskSet(const BitRefLevel _bit) : bit(_bit){};
  void operator()(Problem &p) { *(p.tagBitRefLevelMask) = bit; };
};

/** \brief add finite element to problem
 * \ingroup problems_multi_indices
 */
struct ProblemFiniteElementChangeBitAdd {
  BitFEId f_id;
  ProblemFiniteElementChangeBitAdd(const BitFEId _f_id) : f_id(_f_id){};
  void operator()(Problem &p);
};

/** \brief set prof dof bit ref mask
 * \ingroup problems_multi_indices
 */
struct ProblemChangeRefLevelBitDofMaskAdd {
  BitRefLevel bit;
  ProblemChangeRefLevelBitDofMaskAdd(const BitRefLevel _bit) : bit(_bit){};
  void operator()(Problem &p) { *(p.tagBitRefLevelMask) |= bit; };
};

/** \brief remove finite element from problem
 * \ingroup problems_multi_indices
 */
struct ProblemFiniteElementChangeBitUnSet {
  BitFEId f_id;
  ProblemFiniteElementChangeBitUnSet(const BitFEId _f_id) : f_id(_f_id){};
  void operator()(Problem &p);
};

/** \brief zero nb. of DOFs in row
 * \ingroup problems_multi_indices
 */
struct ProblemZeroNbRowsChange {
  void operator()(Problem &e);
};

/** \brief zero nb. of DOFs in col
 * \ingroup problems_multi_indices
 */
struct ProblemZeroNbColsChange {
  void operator()(Problem &e);
};

/** \brief clear problem finite elements
 * \ingroup problems_multi_indices
 */
struct ProblemClearNumeredFiniteElementsChange {
  void operator()(Problem &e);
};

/**
 * \brief Clear sub-problem data structure
 */
struct ProblemClearSubProblemData {
  void operator()(Problem &e) { e.subProblemData.reset(); }
};

/**
 * \brief Clear composed problem data structure
 */
struct ProblemClearComposedProblemData {
  void operator()(Problem &e) { e.composedProblemsData.reset(); }
};

inline MoFEMErrorCode
ComposedProblemsData::getRowIs(IS *is, const unsigned int pp) const {
  MoFEMFunctionBeginHot;
  PetscObjectReference((PetscObject)rowIs[pp]);
  if (pp <= rowIs.size()) {
    SETERRQ1(PETSC_COMM_WORLD, MOFEM_INVALID_DATA, "Exceed size of array pp<%d",
             rowIs.size());
  }
  *is = rowIs[pp];
  MoFEMFunctionReturnHot(0);
}

inline MoFEMErrorCode
ComposedProblemsData::getColIs(IS *is, const unsigned int pp) const {
  MoFEMFunctionBeginHot;
  PetscObjectReference((PetscObject)colIs[pp]);
  if (pp <= colIs.size()) {
    SETERRQ1(PETSC_COMM_WORLD, MOFEM_INVALID_DATA, "Exceed size of array pp<%d",
             colIs.size());
  }
  *is = colIs[pp];
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM

#endif //__PROBLEMSMULTIINDICES_HPP__

/**
 * \defgroup problems_multi_indices Problems structures and multi-indices
 * \ingroup mofem
 ******************************************************************************/
