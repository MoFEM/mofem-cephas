/** \file ProblemsMultiIndices.hpp
 * \brief Multi-index containers, data structures for problems and other
 * low-level functions
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

  inline MoFEMErrorCode getRowIs(IS *is, const unsigned int pp) const {
    MoFEMFunctionBeginHot;
    PetscObjectReference((PetscObject)rowIs[pp]);
    if (pp <= rowIs.size()) {
      SETERRQ1(PETSC_COMM_WORLD, MOFEM_INVALID_DATA,
               "Exceed size of array pp<%d", rowIs.size());
    }
    *is = rowIs[pp];
    MoFEMFunctionReturnHot(0);
  }

  inline MoFEMErrorCode getColIs(IS *is, const unsigned int pp) const {
    MoFEMFunctionBeginHot;
    PetscObjectReference((PetscObject)colIs[pp]);
    if (pp <= colIs.size()) {
      SETERRQ1(PETSC_COMM_WORLD, MOFEM_INVALID_DATA,
               "Exceed size of array pp<%d", colIs.size());
    }
    *is = colIs[pp];
    MoFEMFunctionReturnHot(0);
  }

  virtual ~ComposedProblemsData() {
    for (unsigned int ii = 0; ii != rowIs.size(); ii++) {
      ISDestroy(&rowIs[ii]);
    }
    for (unsigned int jj = 0; jj != colIs.size(); jj++) {
      ISDestroy(&colIs[jj]);
    }
  }
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
  BitRefLevel *tagMaskBitRefLevel; ///< BItRefMask of elements in problem

  mutable DofIdx nbDofsRow;      ///< Global number of DOFs in  row
  mutable DofIdx nbDofsCol;      ///< Global number of DOFs in col
  mutable DofIdx nbLocDofsRow;   ///< Local number of DOFs in row
  mutable DofIdx nbLocDofsCol;   ///< Local number of DOFs in colIs
  mutable DofIdx nbGhostDofsRow; ///< Number of ghost DOFs in row
  mutable DofIdx nbGhostDofsCol; ///< Number of ghost DOFs in col

  mutable boost::shared_ptr<NumeredDofEntity_multiIndex>
      numeredDofsRows; ///< store DOFs on rows for this problem
  mutable boost::shared_ptr<NumeredDofEntity_multiIndex>
      numeredDofsCols; ///< store DOFs on columns for this problem
  mutable boost::shared_ptr<NumeredEntFiniteElement_multiIndex>
      numeredFiniteElements; ///< store finite elements

  /**
   * \brief get access to numeredDofsRows storing DOFs on rows
   */
  const boost::shared_ptr<NumeredDofEntity_multiIndex> &
  getNumeredDofsRows() const {
    return numeredDofsRows;
  }

  /**
   * \brief get access to numeredDofsCols storing DOFs on cols
   */
  const boost::shared_ptr<NumeredDofEntity_multiIndex> &
  getNumeredDofsCols() const {
    return numeredDofsCols;
  }

  /**
   * \brief get access to reference for multi-index storing finite elements
   */
  const boost::shared_ptr<NumeredEntFiniteElement_multiIndex> &
  getNumeredFiniteElements() const {
    return numeredFiniteElements;
  }

  /**
   * \brief Subproblem problem data
   */
  struct SubProblemData {

    IS rowIs;  ///< indices of main problem of which sub problem is this
    IS colIs;  ///< indices of main problem of which sub problem is this
    AO rowMap; ///< mapping form main problem indices to sub-problem indices
    AO colMap;

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

    SubProblemData()
        : rowIs(PETSC_NULL), colIs(PETSC_NULL), rowMap(PETSC_NULL),
          colMap(PETSC_NULL) {}
    virtual ~SubProblemData() {
      // int flg;
      // MPI_Finalized(&flg);
      // if(!flg) {
      ISDestroy(&rowIs);
      ISDestroy(&colIs);
      AODestroy(&rowMap);
      AODestroy(&colMap);
      // }
    }
  };

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
  inline boost::shared_ptr<ComposedProblemsData> &
  getComposedProblemsData() const {
    return composedProblemsData;
  }

  /**
   * \brief get DOFs from problem
   *
   * Note that \e ent_dof_idx is not coefficient number, is local number of DOFs
   * on the entity. The coefficient number and local index of DOFs or entity are
   * the same on vertices and H1 approximation.
   *
   * @param  name       field name
   * @param  ent        entity handle
   * @param  ent_dof_idx index of DOFs on entity
   * @param  row_or_col ROW or COL
   * @param  dof_ptr    shared pointer to DOFs if found
   * @return            error code
   */
  MoFEMErrorCode getDofByNameEntAndEntDofIdx(
      const string name, const EntityHandle ent, const int ent_dof_idx,
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
  NumeredEntFiniteElementbyName::iterator IT =                                 \
      PROBLEMPTR->getNumeredFEsBegin(NAME);                                    \
  IT != PROBLEMPTR->getNumeredFEsEnd(NAME);                                    \
  IT++

  NumeredEntFiniteElementbyName::iterator
  getNumeredFEsBegin(std::string fe_name) const {
    return numeredFiniteElements->get<FiniteElement_name_mi_tag>().lower_bound(
        fe_name);
  }

  NumeredEntFiniteElementbyName::iterator
  getNumeredFEsEnd(std::string fe_name) const {
    return numeredFiniteElements->get<FiniteElement_name_mi_tag>().upper_bound(
        fe_name);
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
      PROBLEMPTR->getNumeredDofsRowsBegin();                                   \
  IT != PROBLEMPTR->getNumeredDofsRowsEnd();                                   \
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
      PROBLEMPTR->getNumeredDofsColsBegin();                                   \
  IT != PROBLEMPTR->getNumeredDofsColsEnd();                                   \
  IT++

  /// get begin iterator for numeredDofsRows (insted you can use
  /// #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredDofsRowsBegin() const {
    return numeredDofsRows->begin();
  }

  /// get end iterator for numeredDofsRows (insted you can use
  /// #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredDofsRowsEnd() const {
    return numeredDofsRows->end();
  }

  /// get begin iterator for numeredDofsCols (insted you can use
  /// #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredDofsColsBegin() const {
    return numeredDofsCols->begin();
  }

  /// get end iterator for numeredDofsCols (insted you can use
  /// #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofEntity_multiIndex::iterator getNumeredDofsColsEnd() const {
    return numeredDofsCols->end();
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
      PROBLEMPTR->getNumeredDofsRowsByLocIdxBegin(0);                          \
  IT != PROBLEMPTR->getNumeredDofsRowsByLocIdxEnd(                             \
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
      PROBLEMPTR->getNumeredDofsColsByLocIdxBegin(0);                          \
  IT != PROBLEMPTR->getNumeredDofsColsByLocIdxEnd(                             \
            PROBLEMPTR->getNbLocalDofsRow() - 1);                              \
  IT++

  /// get begin iterator for numeredDofsRows (insted you can use
  /// #_IT_NUMEREDDOF_ROW_FOR_LOOP_ for loops)
  NumeredDofEntityByLocalIdx::iterator
  getNumeredDofsRowsByLocIdxBegin(const DofIdx locidx) const {
    return numeredDofsRows->get<PetscLocalIdx_mi_tag>().lower_bound(locidx);
  }

  /// get end iterator for numeredDofsRows (insted you can use
  /// #_IT_NUMEREDDOF_ROW_FOR_LOOP_ for loops)
  NumeredDofEntityByLocalIdx::iterator
  getNumeredDofsRowsByLocIdxEnd(const DofIdx locidx) const {
    return numeredDofsRows->get<PetscLocalIdx_mi_tag>().upper_bound(locidx);
  }

  /// get begin iterator for numeredDofsCols (insted you can use
  /// #_IT_NUMEREDDOF_COL_FOR_LOOP_ for loops)
  NumeredDofEntityByLocalIdx::iterator
  getNumeredDofsColsByLocIdxBegin(const DofIdx locidx) const {
    return numeredDofsCols->get<PetscLocalIdx_mi_tag>().lower_bound(locidx);
  }

  /// get end iterator for numeredDofsCols (insted you can use
  /// #_IT_NUMEREDDOF_COL_FOR_LOOP_ for loops)
  NumeredDofEntityByLocalIdx::iterator
  getNumeredDofsColsByLocIdxEnd(const DofIdx locidx) const {
    return numeredDofsCols->get<PetscLocalIdx_mi_tag>().upper_bound(locidx);
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
  NumeredDofEntityByEnt::iterator IT =                                         \
      PROBLEMPTR->getNumeredDofsRowsByEntBegin(ENT);                           \
  IT != PROBLEMPTR->getNumeredDofsRowsByEntEnd(ENT);                           \
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
  NumeredDofEntityByEnt::iterator IT =                                         \
      PROBLEMPTR->getNumeredDofsColsByEntBegin(ENT);                           \
  IT != PROBLEMPTR->getNumeredDofsColsByEntEnd(ENT);                           \
  IT++

  /// get begin iterator for numeredDofsRows (insted you can use
  /// #_IT_NUMEREDDOF_ROW_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofEntityByEnt::iterator
  getNumeredDofsRowsByEntBegin(const EntityHandle ent) const {
    return numeredDofsRows->get<Ent_mi_tag>().lower_bound(ent);
  }

  /// get end iterator for numeredDofsRows (insted you can use
  /// #_IT_NUMEREDDOF_ROW_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofEntityByEnt::iterator
  getNumeredDofsRowsByEntEnd(const EntityHandle ent) const {
    return numeredDofsRows->get<Ent_mi_tag>().upper_bound(ent);
  }

  /// get begin iterator for numeredDofsCols (insted you can use
  /// #_IT_NUMEREDDOF_COL_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofEntityByEnt::iterator
  getNumeredDofsColsByEntBegin(const EntityHandle ent) const {
    return numeredDofsCols->get<Ent_mi_tag>().lower_bound(ent);
  }

  /// get end iterator for numeredDofsCols (insted you can use
  /// #_IT_NUMEREDDOF_COL_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofEntityByEnt::iterator
  getNumeredDofsColsByEntEnd(const EntityHandle ent) const {
    return numeredDofsCols->get<Ent_mi_tag>().upper_bound(ent);
  }

/**
 * use with loops to iterate row DOFs
 * \ingroup problems_multi_indices
 *
 * \code
 * for(_IT_NUMEREDDOF_BY_NAME_ROW_FOR_LOOP_(PROBLEMPTR,NAME,IT)) {
 *   ...
 * }
 * \endcode
 *
 */
#define _IT_NUMEREDDOF_ROW_BY_NAME_FOR_LOOP_(PROBLEMPTR, NAME, IT)             \
  NumeredDofEntityByFieldName::iterator IT =                                   \
      PROBLEMPTR->getNumeredDofsRowsBegin(NAME);                               \
  IT != PROBLEMPTR->getNumeredDofsRowsEnd(NAME);                               \
  IT++

/**
 * \brief use with loops to iterate col DOFs
 * \ingroup problems_multi_indices
 *
 * \code
 * for(_IT_NUMEREDDOF_COL_BY_NAME_FOR_LOOP_(PROBLEMPTR,NAME,IT)) {
 *   ...
 * }
 * \endcode
 *
 */
#define _IT_NUMEREDDOF_COL_BY_NAME_FOR_LOOP_(PROBLEMPTR, NAME, IT)             \
  NumeredDofEntityByFieldName::iterator IT =                                   \
      PROBLEMPTR->getNumeredDofsColsBegin(NAME);                               \
  IT != PROBLEMPTR->getNumeredDofsColsEnd(NAME);                               \
  IT++

  Problem(Interface &moab, const EntityHandle meshset);

  virtual ~Problem();

  inline BitProblemId getId() const { return *((BitProblemId *)tagId); }

  inline std::string getName() const {
    return std::string((char *)tagName, tagNameSize);
  }

  inline DofIdx getNbDofsRow() const { return nbDofsRow; }
  inline DofIdx getNbDofsCol() const { return nbDofsCol; }
  inline DofIdx getNbLocalDofsRow() const { return nbLocDofsRow; }
  inline DofIdx getNbLocalDofsCol() const { return nbLocDofsCol; }
  inline DofIdx getNbGhostDofsRow() const { return nbGhostDofsRow; }
  inline DofIdx getNbGhostDofsCol() const { return nbGhostDofsCol; }

  inline BitRefLevel getBitRefLevel() const { return *tagBitRefLevel; }
  inline BitRefLevel getMaskBitRefLevel() const { return *tagMaskBitRefLevel; }

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

  typedef multi_index_container<boost::weak_ptr<std::vector<NumeredDofEntity> >,
                                indexed_by<sequenced<> > >
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
  inline boost::shared_ptr<SequenceDofContainer> &getRowDofsSequence() const {
    return sequenceRowDofContainer;
  }

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
  inline boost::shared_ptr<SequenceDofContainer> &getColDofsSequence() const {
    return sequenceColDofContainer;
  }

private:
  // Keep vector of DoFS on entity
  mutable boost::shared_ptr<SequenceDofContainer> sequenceRowDofContainer;
  mutable boost::shared_ptr<SequenceDofContainer> sequenceColDofContainer;
};

/// \deprecated use just Problem
DEPRECATED typedef Problem MoFEMProblem;

/**
 * @relates multi_index_container
 * \brief MultiIndex for entities for Problem
 * \ingroup fe_multi_indices
 */
typedef multi_index_container<
    Problem,
    indexed_by<
        ordered_unique<tag<Meshset_mi_tag>,
                       member<Problem, EntityHandle, &Problem::meshset> >,
        hashed_unique<tag<BitProblemId_mi_tag>,
                      const_mem_fun<Problem, BitProblemId, &Problem::getId>,
                      HashBit<BitProblemId>, EqBit<BitProblemId> >,
        hashed_unique<
            tag<Problem_mi_tag>,
            const_mem_fun<Problem, std::string, &Problem::getName> > > >
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
  void operator()(Problem &p) { *(p.tagMaskBitRefLevel) = bit; };
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
  void operator()(Problem &p) { *(p.tagMaskBitRefLevel) |= bit; };
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

} // namespace MoFEM

#endif //__PROBLEMSMULTIINDICES_HPP__

/**
 * \defgroup problems_multi_indices Problems structures and multi-indices
 * \ingroup mofem
 ******************************************************************************/
