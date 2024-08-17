/**
 * @file Schur.cpp
 * @brief Implementation of Schur Complement
 * @date 2023-02-03
 *
 * @copyright Copyright (c) 2023
 *
 */

namespace MoFEM {

constexpr bool debug_schur = false;

struct OpSchurAssembleBaseImpl : public OpSchurAssembleBase {
  OpSchurAssembleBaseImpl() = delete;

protected:
  using OpSchurAssembleBase::OpSchurAssembleBase;

  /**
   * @brief Assemble Schur complement
   *
   */
  inline MoFEMErrorCode assembleSchurMat(

      Mat S,

      const UId &uid_row, VectorInt &row_ind,

      const UId &uid_col, VectorInt &col_ind,

      MatrixDouble &m, InsertMode iora

  );
};

// Assemble specialisations
MoFEMErrorCode OpSchurAssembleBaseImpl::assembleSchurMat(

    Mat S,

    const UId &uid_row, VectorInt &row_ind,

    const UId &uid_col, VectorInt &col_ind,

    MatrixDouble &m, InsertMode iora

) {
  return ::MatSetValues(

      S, row_ind.size(), &*row_ind.begin(), col_ind.size(), &*col_ind.begin(),
      &*m.data().begin(), iora

  );
}

/**
 * @brief Clear Schur complement internal data
 *
 */
struct OpSchurAssembleBegin : public OpSchurAssembleBaseImpl {

  OpSchurAssembleBegin();

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

/**
 * @brief Assemble Schur complement (Implementation)
 *
 */
template <typename OP_SCHUR_ASSEMBLE_BASE>
struct OpSchurAssembleEndImpl : public OP_SCHUR_ASSEMBLE_BASE {

  using OP = OP_SCHUR_ASSEMBLE_BASE;

  /**
   * @brief Construct a new Op Schur Assemble End object
   *
   * @param fields_name list of fields
   * @param field_ents list of entities on which schur complement is applied
   * (can be empty)
   * @param schur_ao map indices to schur matrix indices
   * @param schur_mat schur matrix
   * @param sym_schur true if schur is symmetric
   * @param symm_op true if block diagonal is symmetric
   */
  OpSchurAssembleEndImpl(

      std::vector<std::string> fields_name,
      std::vector<boost::shared_ptr<Range>> field_ents,

      SmartPetscObj<AO> schur_ao, SmartPetscObj<Mat> schur_mat,

      bool sym_schur, bool symm_op,

      boost::shared_ptr<BlockStructure> diag_blocks);

protected:
  template <typename I>
  MoFEMErrorCode doWorkImpl(

      int side, EntityType type, EntitiesFieldData::EntData &data

  );

  std::vector<std::string> fieldsName;
  std::vector<boost::shared_ptr<Range>> fieldEnts;
  SmartPetscObj<AO> schurAO;
  SmartPetscObj<Mat> schurMat;
  bool symSchur;
  boost::shared_ptr<BlockStructure> diagBlocks;

  MatrixDouble blockMat;
  MatrixDouble invMat;
  MatrixDouble bM, abM, abcM;
};

struct SchurDSYSV; ///< SY	symmetric
struct SchurDGESV; ///< GE	general (i.e., nonsymmetric, in some cases
                   ///< rectangular)

/**
 * @brief Assemble Schur complement
 *
 */
template <typename I> struct OpSchurAssembleEnd;

template <>
struct OpSchurAssembleEnd<SchurDSYSV>
    : public OpSchurAssembleEndImpl<OpSchurAssembleBaseImpl> {
  using OpSchurAssembleEndImpl<OpSchurAssembleBaseImpl>::OpSchurAssembleEndImpl;
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

template <>
struct OpSchurAssembleEnd<SchurDGESV>
    : public OpSchurAssembleEndImpl<OpSchurAssembleBaseImpl> {
  using OpSchurAssembleEndImpl<OpSchurAssembleBaseImpl>::OpSchurAssembleEndImpl;
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

/**
 * @brief Schur complement data storage
 *
 */
struct SchurElemMats : public boost::enable_shared_from_this<SchurElemMats> {

  SchurElemMats(const size_t idx, const UId uid_row, const UId uid_col);
  virtual ~SchurElemMats() = default;

  const size_t iDX;
  UId uidRow;
  UId uidCol;

  inline auto &getMat() const { return locMats[iDX]; }
  inline auto &getRowInd() const { return rowIndices[iDX]; }
  inline auto &getColInd() const { return colIndices[iDX]; }

  static MoFEMErrorCode MatSetValues(Mat M,
                                     const EntitiesFieldData::EntData &row_data,
                                     const EntitiesFieldData::EntData &col_data,
                                     const MatrixDouble &mat, InsertMode iora);

protected:
  static MoFEMErrorCode
  assembleStorage(const EntitiesFieldData::EntData &row_data,
                  const EntitiesFieldData::EntData &col_data,
                  const MatrixDouble &mat, InsertMode iora);

  friend struct OpSchurAssembleBegin;
  template <typename OP_SCHUR_ASSEMBLE_BASE>
  friend struct OpSchurAssembleEndImpl;

  struct uid_mi_tag {};
  struct col_mi_tag {};

  using SchurElemStorage = multi_index_container<
      const SchurElemMats *,
      indexed_by<

          ordered_unique<
              tag<uid_mi_tag>,
              composite_key<
                  SchurElemMats,

                  member<SchurElemMats, const UId, &SchurElemMats::uidRow>,
                  member<SchurElemMats, const UId, &SchurElemMats::uidCol>

                  >>,

          ordered_non_unique<tag<col_mi_tag>, member<SchurElemMats, const UId,
                                                     &SchurElemMats::uidCol>>

          >>;

  static boost::ptr_vector<MatrixDouble> locMats;
  static boost::ptr_vector<VectorInt> rowIndices;
  static boost::ptr_vector<VectorInt> colIndices;
  static boost::ptr_vector<SchurElemMats> schurElemMats;
  static size_t maxIndexCounter;

  static SchurElemStorage schurL2Storage;
};

struct DiagBlockIndex {

  virtual ~DiagBlockIndex() = default;

  /**
   * @brief block data indexes
   *
   */
  struct Indexes {

    Indexes(UId uid_row, UId uid_col, int row, int col, int nb_rows,
            int nb_cols, int loc_row, int loc_col, int mat_shift, int inv_shift)
        : uid_row(uid_row), uid_col(uid_col), row(row), col(col),
          nb_rows(nb_rows), nb_cols(nb_cols), loc_row(loc_row),
          loc_col(loc_col), mat_shift(mat_shift), inv_shift(inv_shift) {}

    inline UId getRowUId() const { return uid_row; }
    inline UId getColUId() const { return uid_col; }
    inline int getRow() const { return row; }
    inline int getCol() const { return col; }
    inline int getNbRows() const { return nb_rows; }
    inline int getNbCols() const { return nb_cols; }
    inline int getLocRow() const { return loc_row; }
    inline int getLocCol() const { return loc_col; }
    inline int &getMatShift() const { return mat_shift; }
    inline int &getInvShift() const { return inv_shift; }

    inline int rowShift() const {
      return getRow() + getNbRows();
    } // shift such that lower bound is included

    inline int colShift() const {
      return getCol() + getNbCols();
    } // shift such that lower bound is included

  private:
    UId uid_row;
    UId uid_col;
    int row;
    int col;
    int nb_rows;
    int nb_cols;
    int loc_row;
    int loc_col;
    mutable int mat_shift;
    mutable int inv_shift;
  };

  using BlockIndex = multi_index_container<

      Indexes,

      indexed_by<

          ordered_unique<

              composite_key<Indexes,

                            const_mem_fun<Indexes, UId, &Indexes::getRowUId>,
                            const_mem_fun<Indexes, UId, &Indexes::getColUId>>

              >,

          ordered_non_unique<

              const_mem_fun<Indexes, int, &Indexes::rowShift>

              >,

          ordered_non_unique<

              const_mem_fun<Indexes, int, &Indexes::colShift>

              >

          >>;

  BlockIndex blockIndex; ///< blocks indexes storage
};

struct BlockStructure : public DiagBlockIndex {

  SmartPetscObj<Vec> ghostX;
  SmartPetscObj<Vec> ghostY;

  boost::shared_ptr<std::vector<double>> dataBlocksPtr;
  boost::shared_ptr<std::vector<double>> dataInvBlocksPtr;
  boost::shared_ptr<std::vector<double>> preconditionerBlocksPtr;
  boost::shared_ptr<std::vector<double>> parentBlockStructurePtr;

  bool multiplyByPreconditioner = false;
};

PetscLogEvent SchurEvents::MOFEM_EVENT_schurMatSetValues;
PetscLogEvent SchurEvents::MOFEM_EVENT_opSchurAssembleEnd;
PetscLogEvent SchurEvents::MOFEM_EVENT_BlockStructureSetValues;
PetscLogEvent SchurEvents::MOFEM_EVENT_BlockStructureMult;
PetscLogEvent SchurEvents::MOFEM_EVENT_BlockStructureSolve;
PetscLogEvent SchurEvents::MOFEM_EVENT_zeroRowsAndCols;

SchurEvents::SchurEvents() {
  PetscLogEventRegister("schurMatSetVal", 0, &MOFEM_EVENT_schurMatSetValues);
  PetscLogEventRegister("opSchurAsmEnd", 0, &MOFEM_EVENT_opSchurAssembleEnd);
  PetscLogEventRegister("blockSetVal", 0, &MOFEM_EVENT_BlockStructureSetValues);
  PetscLogEventRegister("blockMult", 0, &MOFEM_EVENT_BlockStructureMult);
  PetscLogEventRegister("blockSolve", 0, &MOFEM_EVENT_BlockStructureSolve);
  PetscLogEventRegister("schurZeroRandC", 0, &MOFEM_EVENT_zeroRowsAndCols);
}

SchurElemMats::SchurElemStorage SchurElemMats::schurL2Storage;
boost::ptr_vector<MatrixDouble> SchurElemMats::locMats;
boost::ptr_vector<VectorInt> SchurElemMats::rowIndices;
boost::ptr_vector<VectorInt> SchurElemMats::colIndices;
boost::ptr_vector<SchurElemMats> SchurElemMats::schurElemMats;
size_t SchurElemMats::maxIndexCounter = 0;

SchurElemMats::SchurElemMats(const size_t idx, const UId uid_row,
                             const UId uid_col)
    : iDX(idx), uidRow(uid_row), uidCol(uid_col) {}

MoFEMErrorCode
SchurElemMats::assembleStorage(const EntitiesFieldData::EntData &row_data,
                               const EntitiesFieldData::EntData &col_data,
                               const MatrixDouble &mat, InsertMode iora) {
  MoFEMFunctionBegin;

#ifndef NDEBUG
  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_schurMatSetValues, 0, 0, 0, 0);
#endif // NDEBUG

  // get row indices, in case of store, get indices from storage
  // storage keeps marked indices to manage boundary conditions
  auto get_row_indices = [&]() -> const VectorInt & {
    if (auto e_ptr = row_data.getFieldEntities()[0]) {
      if (auto stored_data_ptr =
              e_ptr->getSharedStoragePtr<EssentialBcStorage>()) {
        return stored_data_ptr->entityIndices;
      }
    }
    return row_data.getIndices();
  };

  const auto &row_ind = get_row_indices();
  const auto &col_ind = col_data.getIndices();

  const auto nb_rows = row_ind.size();
  const auto nb_cols = col_ind.size();

#ifndef NDEBUG
  if (mat.size1() != nb_rows) {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Wrong mat size %d != %d", mat.size1(), nb_rows);
  }
  if (mat.size2() != nb_cols) {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Wrong mat size %d != %d", mat.size2(), nb_cols);
  }
#endif // NDEBUG

  // get entity uid
  auto get_uid = [](auto &data) {
    if (data.getFieldEntities().size() == 1) {

      return data.getFieldEntities()[0]->getLocalUniqueId();

    } else {

      // Is assumed that sum of entities ids gives unique id, that is not true,
      // but corner case is improbable.

      // @todo: debug version should have test

      auto &uid0 = data.getFieldEntities()[0]->getLocalUniqueId();
      auto field_id0 = FieldEntity::getFieldBitNumberFromUniqueId(uid0);
      auto ent0 = FieldEntity::getHandleFromUniqueId(uid0);
      auto type0 = type_from_handle(ent0);
      auto id = id_from_handle(ent0);

      for (auto i = 1; i < data.getFieldEntities().size(); ++i) {

        // get entity id from ent
        id += id_from_handle(

            // get entity handle from unique uid
            FieldEntity::getHandleFromUniqueId(
                data.getFieldEntities()[i]->getLocalUniqueId())

        );
      }

      return FieldEntity::getLocalUniqueIdCalculate(
          field_id0,

          ent_form_type_and_id(type0, id)

      );
    }
  };

  auto uid_row = get_uid(row_data);
  auto uid_col = get_uid(col_data);

  auto it =
      SchurElemMats::schurL2Storage.template get<SchurElemMats::uid_mi_tag>()
          .find(boost::make_tuple(uid_row, uid_col));

  if (it ==
      SchurElemMats::schurL2Storage.template get<SchurElemMats::uid_mi_tag>()
          .end()) {

    // get size of arrays of matrices
    const auto size = SchurElemMats::locMats.size();

    // expand memory allocation
    if (SchurElemMats::maxIndexCounter == size) {
      SchurElemMats::locMats.push_back(new MatrixDouble());
      SchurElemMats::rowIndices.push_back(new VectorInt());
      SchurElemMats::colIndices.push_back(new VectorInt());
      SchurElemMats::schurElemMats.push_back(
          new SchurElemMats(SchurElemMats::maxIndexCounter, uid_row, uid_col));
    } else {
      SchurElemMats::schurElemMats[SchurElemMats::maxIndexCounter].uidRow =
          uid_row;
      SchurElemMats::schurElemMats[SchurElemMats::maxIndexCounter].uidCol =
          uid_col;
    }

    // add matrix to storage
    auto p = SchurElemMats::schurL2Storage.emplace(
        &SchurElemMats::schurElemMats[SchurElemMats::maxIndexCounter++]);
#ifndef NDEBUG
    if (!p.second) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Failed to insert");
    }
#endif // NDEBUG

    auto asmb = [&](auto &sm) {
      sm.resize(nb_rows, nb_cols, false);
      noalias(sm) = mat;
    };

    asmb((*p.first)->getMat());

    auto add_indices = [](auto &storage, auto &ind) {
      storage.resize(ind.size(), false);
      noalias(storage) = ind;
    };

    add_indices((*p.first)->getRowInd(), row_ind);
    add_indices((*p.first)->getColInd(), col_ind);

  } else {
    // entry (submatrix) already exists

    auto asmb = [&](auto &sm) {
      MoFEMFunctionBeginHot;

#ifndef NDEBUG
      if (sm.size1() != nb_rows) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong mat or storage size %d != %d", sm.size1(), nb_rows);
      }
      if (sm.size2() != nb_cols) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong mat or storage size %d != %d", sm.size2(), nb_cols);
      }
#endif // NDEBUG

      switch (iora) {
      case ADD_VALUES:
        sm += mat;
        break;
      case INSERT_VALUES:
        noalias(sm) = mat;
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                "Assembly type not implemented");
      }
      MoFEMFunctionReturnHot(0);
    };

    CHKERR asmb((*it)->getMat());

    // no need to set indices
  }

#ifndef NDEBUG
  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_schurMatSetValues, 0, 0, 0, 0);
#endif // NDEBUG

  MoFEMFunctionReturn(0);
}

OpSchurAssembleBegin::OpSchurAssembleBegin()
    : OpSchurAssembleBaseImpl(NOSPACE, OPSPACE) {}

MoFEMErrorCode OpSchurAssembleBegin::doWork(int side, EntityType type,
                                            EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
#ifndef NDEBUG
  if constexpr (debug_schur)
    MOFEM_LOG("SELF", Sev::noisy) << "Schur assemble begin";
#endif
  SchurElemMats::schurL2Storage.clear();
  SchurElemMats::maxIndexCounter = 0;

  MoFEMFunctionReturn(0);
}

template <typename OP_SCHUR_ASSEMBLE_BASE>
OpSchurAssembleEndImpl<OP_SCHUR_ASSEMBLE_BASE>::OpSchurAssembleEndImpl(
    std::vector<std::string> fields_name,
    std::vector<boost::shared_ptr<Range>> field_ents,
    SmartPetscObj<AO> schur_ao, SmartPetscObj<Mat> schur_mat, bool sym_schur,
    bool symm_op, boost::shared_ptr<BlockStructure> diag_blocks)
    : OP(NOSPACE, OP::OPSPACE, symm_op), fieldsName(fields_name),
      fieldEnts(field_ents), schurAO(schur_ao), schurMat(schur_mat),
      symSchur(sym_schur), diagBlocks(diag_blocks) {}

template <typename OP_SCHUR_ASSEMBLE_BASE>
template <typename I>
MoFEMErrorCode OpSchurAssembleEndImpl<OP_SCHUR_ASSEMBLE_BASE>::doWorkImpl(
    int side, EntityType type, EntitiesFieldData::EntData &data) {

  MoFEMFunctionBegin;

#ifndef NDEBUG
  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_opSchurAssembleEnd, 0, 0, 0, 0);
#endif

#ifndef NDEBUG
  if constexpr (debug_schur)
    MOFEM_LOG("SELF", Sev::noisy) << "Schur assemble begin -> end";
#endif

  auto get_field_name = [&](auto uid) {
    return OP::getPtrFE()->mField.get_field_name(field_bit_from_bit_number(
        FieldEntity::getFieldBitNumberFromUniqueId(uid)));
  };

  auto get_a00_uids = [&]() {
    auto get_field_bit = [&](auto &name) {
      return OP::getPtrFE()->mField.get_field_bit_number(name);
    };

    std::vector<std::pair<UId, UId>> a00_uids;
    a00_uids.reserve(fieldsName.size());
    for (auto ss = 0; ss != fieldsName.size(); ++ss) {
      auto field_bit = get_field_bit(fieldsName[ss]);
      auto row_ents = fieldEnts[ss];
      if (row_ents) {
        for (auto p = row_ents->pair_begin(); p != row_ents->pair_end(); ++p) {
          auto lo_uid =
              FieldEntity::getLoLocalEntityBitNumber(field_bit, p->first);
          auto hi_uid =
              FieldEntity::getHiLocalEntityBitNumber(field_bit, p->second);
          a00_uids.push_back(std::make_pair(lo_uid, hi_uid));
        }
      } else {
        auto lo_uid = FieldEntity::getLoLocalEntityBitNumber(
            field_bit, get_id_for_min_type<MBVERTEX>());
        auto hi_uid = FieldEntity::getHiLocalEntityBitNumber(
            field_bit, get_id_for_max_type<MBENTITYSET>());
        a00_uids.push_back(std::make_pair(lo_uid, hi_uid));
      }
    }
    return a00_uids;
  };

#ifndef NDEBUG
  auto list_storage = [&](auto &storage) {
    MoFEMFunctionBegin;
    int i = 0;
    for (auto &p : storage) {
      MOFEM_LOG("SELF", Sev::noisy)
          << "List schur storage: " << i << " " << p->iDX << ": "
          << get_field_name(p->uidRow) << " " << get_field_name(p->uidCol)
          << " : " << p->getMat().size1() << " " << p->getMat().size2();
      ++i;
    }
    MoFEMFunctionReturn(0);
  };
#endif // NDEBUG

  auto assemble_dense_blocks = [&]() {
    using matrix_range = ublas::matrix_range<MatrixDouble>;
    using range = ublas::range;
    MoFEMFunctionBegin;
    auto &storage = SchurElemMats::schurL2Storage;

    auto assemble_schur = [this](auto &m, auto &uid_row, auto &uid_col,
                                 auto *row_ind_ptr, auto *col_ind_ptr) {
      MoFEMFunctionBegin;

      if (schurAO) {
        CHKERR AOApplicationToPetsc(schurAO, row_ind_ptr->size(),
                                    &*row_ind_ptr->begin());
        CHKERR AOApplicationToPetsc(schurAO, col_ind_ptr->size(),
                                    &*col_ind_ptr->begin());
      }

      if (schurMat) {

        if (auto ierr = this->assembleSchurMat(

                schurMat, uid_row, *row_ind_ptr, uid_col, *col_ind_ptr, m,
                ADD_VALUES

                )) {
#ifndef NDEBUG
          auto field_ents = OP::getPtrFE()->mField.get_field_ents();
          auto row_ent_it = field_ents->find(uid_row);
          auto col_ent_it = field_ents->find(uid_col);
          MOFEM_LOG_CHANNEL("SELF");
          if (row_ent_it != field_ents->end())
            MOFEM_LOG("SELF", Sev::error)
                << "Assemble row entity: " << (*row_ent_it)->getName() << " "
                << (*col_ent_it)->getEntTypeName() << " side "
                << (*row_ent_it)->getSideNumber();
          if (col_ent_it != field_ents->end())
            MOFEM_LOG("SELF", Sev::error)
                << "Assemble col entity: " << (*col_ent_it)->getName() << " "
                << (*col_ent_it)->getEntTypeName() << " side "
                << (*col_ent_it)->getSideNumber();
#endif // NDEBUG
          CHK_THROW_MESSAGE(ierr, "MatSetValues");
        }
      }

      MoFEMFunctionReturn(0);
    };

    auto assemble_a00 = [&](auto &m, auto &row_uid, auto &col_uid,
                            auto *row_ind_ptr, auto *col_ind_ptr) {
      MoFEMFunctionBegin;

#ifndef NDEBUG
      if (!diagBlocks->dataInvBlocksPtr)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "No dataInvBlocksPtr");
#endif // NDEBUG

      if (m.size1() && m.size2()) {
        auto it = diagBlocks->blockIndex.get<0>().find(
            boost::make_tuple(row_uid, col_uid));
        if (it != diagBlocks->blockIndex.get<0>().end()) {
          auto inv_shift = it->getInvShift();
          if (inv_shift != -1) {
            auto *ptr = &((*diagBlocks->dataInvBlocksPtr)[inv_shift]);
            if (m.size1() != it->getNbRows() || m.size2() != it->getNbCols()) {
              for (auto i = 0; i < m.size1(); ++i) {
                for (auto j = 0; j < m.size2(); ++j) {
                  if ((*row_ind_ptr)[i] != -1 && (*col_ind_ptr)[j] != -1) {
                    *ptr = m(i, j);
                    ++ptr;
                  }
                }
              }
            } else {
              // assemble of diag terms, witch might be changed by Schur
              // complement
              std::copy(m.data().begin(), m.data().end(), ptr);
            }
          } else {
            MOFEM_LOG("SELF", Sev::error)
                << "No blockIndex for " << get_field_name(row_uid) << " "
                << get_field_name(col_uid);
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "No inv_shift");
          }
        } else {
          MOFEM_LOG("SELF", Sev::error)
              << "No blockIndex for " << get_field_name(row_uid) << " "
              << get_field_name(col_uid);
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "No blockIndex");
        }
      }

      MoFEMFunctionReturn(0);
    };

    auto a00_uids = get_a00_uids();

    auto get_block_indexing = [&](auto &a00_uids) {
      // iterate over a00 uids and find blocks
      std::vector<const SchurElemMats *> block_list;
      block_list.reserve(storage.size());
      for (auto &rp_uid : a00_uids) {
        auto [rlo_uid, rhi_uid] = rp_uid;
        for (auto &cp_uid : a00_uids) {
          auto [clo_uid, chi_uid] = cp_uid;

          auto it =
              storage.template get<SchurElemMats::uid_mi_tag>().lower_bound(
                  boost::make_tuple(rlo_uid, clo_uid));
          auto hi_it =
              storage.template get<SchurElemMats::uid_mi_tag>().upper_bound(
                  boost::make_tuple(rhi_uid, chi_uid));

          for (; it != hi_it; ++it) {
            if ((*it)->uidRow >= rlo_uid && (*it)->uidRow < rhi_uid &&
                (*it)->uidCol >= clo_uid && (*it)->uidCol < chi_uid) {
              block_list.push_back(*it);
            }
          }
        }
      }

      // create block indexes map for blockMat
      std::map<UId, std::pair<size_t, const VectorInt *>>
          block_indexing; // uid block map
      for (auto d : block_list) {
        if (block_indexing.find(d->uidRow) == block_indexing.end()) {
          block_indexing[d->uidRow] =
              std::make_pair(d->getRowInd().size(), &(d->getRowInd()));
        }
        if (block_indexing.find(d->uidCol) == block_indexing.end()) {
          block_indexing[d->uidCol] =
              std::make_pair(d->getColInd().size(), &(d->getColInd()));
        }
      }

      // set indexes to block
      int mat_block_size = 0; // size of block matrix
      for (auto &p_uid : a00_uids) {
        auto [lo_uid, hi_uid] = p_uid;
        auto lo = block_indexing.lower_bound(lo_uid);
        auto up = block_indexing.upper_bound(hi_uid);
        for (; lo != up; ++lo) {
          lo->second.first = mat_block_size;
          mat_block_size += lo->second.second->size();
        }
      }

      return std::make_tuple(block_list, block_indexing, mat_block_size);
    };

    auto [block_list, block_indexing, block_mat_size] =
        get_block_indexing(a00_uids);

    for (auto &s : storage) {
      auto &m = s->getMat();
      VectorInt row_ind, col_ind;
      row_ind = s->getRowInd();
      col_ind = s->getColInd();
      CHKERR assemble_schur(m, s->uidRow, s->uidCol, &(row_ind), &(col_ind));
    }

    if (block_mat_size == 0) {
      MoFEMFunctionReturnHot(0);
    }

    blockMat.resize(block_mat_size, block_mat_size, false);
    blockMat.clear();

    auto get_range = [](auto &bi) {
      return range(bi.first, bi.first + bi.second->size());
    };

    for (auto &s : block_list) {
      auto &m = s->getMat();
#ifndef NDEBUG
      if (block_indexing.find(s->uidRow) == block_indexing.end())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Wrong rlo_uid");
      if (block_indexing.find(s->uidCol) == block_indexing.end())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Wrong clo_uid");
#endif // NDEBUG

      auto &rbi = block_indexing.at(s->uidRow);
      auto &cbi = block_indexing.at(s->uidCol);

      auto sub_mat = matrix_range(blockMat, get_range(rbi), get_range(cbi));
      sub_mat = m;
    }

    auto get_zeroed_indices = [&](auto extractor_uid, auto extractor_ind) {
      std::vector<int> zeroed_indices;
      zeroed_indices.reserve(block_mat_size);
      for (auto &s : block_list) {
        auto &bi = block_indexing.at(extractor_uid(s));
        auto &ind = extractor_ind(s);
        for (auto it = ind.begin(); it != ind.end(); ++it) {
          if (*it < 0) {
            auto idx = bi.first + std::distance(ind.begin(), it);
            zeroed_indices.push_back(idx);
          }
        }
      }
      std::sort(zeroed_indices.begin(), zeroed_indices.end());
      auto it = std::unique(zeroed_indices.begin(), zeroed_indices.end());
      zeroed_indices.resize(std::distance(zeroed_indices.begin(), it));
      return zeroed_indices;
    };
    auto zero_rows = get_zeroed_indices(
        [](auto &s) { return s->uidRow; },
        [](auto &s) -> VectorInt & { return s->getRowInd(); });
    auto zero_cols = get_zeroed_indices(
        [](auto &s) { return s->uidCol; },
        [](auto &s) -> VectorInt & { return s->getColInd(); });

    for (auto i : zero_rows) {
      for (auto j = 0; j != blockMat.size2(); ++j) {
        blockMat(i, j) = 0;
      }
    }
    for (auto j : zero_cols) {
      for (auto i = 0; i != blockMat.size1(); ++i) {
        blockMat(i, j) = 0;
      }
    }
    for (auto i : zero_rows) {
      blockMat(i, i) = 1;
    }

    CHKERR I::invertMat(blockMat, invMat);

    // clear storage and block list from a00 blocks, no more needed
    for (auto &s : block_list) {
      auto it = storage.template get<SchurElemMats::uid_mi_tag>().find(
          boost::make_tuple(s->uidRow, s->uidCol));
      storage.template get<SchurElemMats::uid_mi_tag>().erase(it);
    }
    block_list.clear();

    std::vector<const SchurElemMats *> schur_block_list;

    for (

        auto rp_uid_it = a00_uids.begin(); rp_uid_it != a00_uids.end();
        ++rp_uid_it

    ) {
      auto [rlo_uid, rhi_uid] = *rp_uid_it;
      for (auto rm = block_indexing.lower_bound(rlo_uid);
           rm != block_indexing.upper_bound(rhi_uid); ++rm) {
        auto &rbi = rm->second;

        auto a_lo_tmp =
            storage.template get<SchurElemMats::col_mi_tag>().lower_bound(
                rm->first);
        auto a_hi =
            storage.template get<SchurElemMats::col_mi_tag>().upper_bound(
                rm->first);

        for (

            auto cp_uid_it = (symSchur) ? rp_uid_it : a00_uids.begin();
            cp_uid_it != a00_uids.end(); ++cp_uid_it

        ) {
          auto [clo_uid, chi_uid] = *cp_uid_it;
          for (auto cm = block_indexing.lower_bound(clo_uid);
               cm != block_indexing.upper_bound(chi_uid); ++cm) {
            auto &cbi = cm->second;

            auto c_lo_tmp =
                storage.template get<SchurElemMats::uid_mi_tag>().lower_bound(
                    boost::make_tuple(cm->first, 0));
            auto c_hi =
                storage.template get<SchurElemMats::uid_mi_tag>().upper_bound(
                    boost::make_tuple(cm->first, FieldEntity::getHiBitNumberUId(
                                                     BITFIELDID_SIZE - 1)));

            auto sub_inv_mat =
                matrix_range(invMat, get_range(rbi), get_range(cbi));
            bM.resize(sub_inv_mat.size1(), sub_inv_mat.size2());
            noalias(bM) = sub_inv_mat;

            if (diagBlocks) {
              CHKERR assemble_a00(

                  bM, rm->first, cm->first, rbi.second, cbi.second

              );
            }

            for (auto a_lo = a_lo_tmp; a_lo != a_hi; ++a_lo) {
#ifndef NDEBUG
              if (block_indexing.find((*a_lo)->uidRow) != block_indexing.end())
                SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                        "Wrong a_lo->uidRow");
#endif

              auto &a = (*a_lo)->getMat();
              abM.resize(a.size1(), bM.size2(), false);
              abM.clear();
              cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,

                          a.size1(), bM.size2(), a.size2(), 1.,

                          &*a.data().begin(), a.size2(),

                          &*bM.data().begin(), bM.size2(), 0.,

                          &*abM.data().begin(), abM.size2());

              for (auto c_lo = c_lo_tmp; c_lo != c_hi; ++c_lo) {
#ifndef NDEBUG
                if (block_indexing.find((*c_lo)->uidCol) !=
                    block_indexing.end())
                  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                          "Wrong a_lo->uidRow");
#endif

                auto &c = (*c_lo)->getMat();

                abcM.resize(abM.size1(), c.size2(), false);
                abcM.clear();

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            abM.size1(), c.size2(), abM.size2(), -1.,

                            &*abM.data().begin(), abM.size2(),

                            &*c.data().begin(), c.size2(), 0.,

                            &*abcM.data().begin(), abcM.size2());

                VectorInt row_ind, col_ind;
                row_ind = (*a_lo)->getRowInd();
                col_ind = (*c_lo)->getColInd();
                CHKERR assemble_schur(abcM, (*a_lo)->uidRow, (*c_lo)->uidCol,
                                      &(row_ind), &(col_ind));
                if (symSchur && rp_uid_it != cp_uid_it) {
                  abcM = trans(abcM);
                  CHKERR assemble_schur(abcM, (*c_lo)->uidCol, (*a_lo)->uidRow,
                                        &(col_ind), &(row_ind));
                }
              }
            }

            if (diagBlocks) {
              if (symSchur && rp_uid_it != cp_uid_it) {
                bM = trans(bM);
                CHKERR assemble_a00(

                    bM, cm->first, rm->first, cbi.second, rbi.second

                );
              }
            }
          }
        }
      }
    }

    MoFEMFunctionReturn(0);
  };

  // Assemble Schur complements
  CHKERR assemble_dense_blocks();

#ifndef NDEBUG
  if constexpr (debug_schur)
    MOFEM_LOG("SELF", Sev::noisy) << "Schur assemble done";
#endif

#ifndef NDEBUG
  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_opSchurAssembleEnd, 0, 0, 0, 0);
#endif

  MoFEMFunctionReturn(0);
}

struct SchurDSYSV {

  static auto invertMat(MatrixDouble &m, MatrixDouble &inv) {
    MoFEMFunctionBegin;

    const auto nb = m.size1();
#ifndef NDEBUG
    if (nb != m.size2()) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "It should be square matrix %d != %d", nb, m.size2());
    }
#endif

    inv.resize(nb, nb, false);
    inv.swap(m);
    m.resize(2 * nb, 2 * nb, false);

    VectorInt ipiv(nb);
    int info;

    // FIXME: add lapack_dsytrf and lapack_dsytri

    info = lapack_dgetrf(nb, nb, &*inv.data().begin(), nb, &*ipiv.begin());
    if (info)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "lapack error info = %d", info);
    info = lapack_dgetri(nb, &*inv.data().begin(), nb, &*ipiv.begin(),
                         &*m.data().begin(), m.data().size());
    if (info)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "lapack error info = %d", info);

    MoFEMFunctionReturn(0);
  }
};

struct SchurDGESV {

  static auto invertMat(MatrixDouble &m, MatrixDouble &inv) {
    MoFEMFunctionBeginHot;

    const auto nb = m.size1();
#ifndef NDEBUG
    if (nb != m.size2()) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "It should be square matrix %d != %d", nb, m.size2());
    }
#endif

    inv.resize(nb, nb, false);
    inv.swap(m);

    VectorInt ipiv(nb);
    int info;

    info = lapack_dgetrf(nb, nb, &*inv.data().begin(), nb, &*ipiv.begin());
    if (info)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "lapack error info = %d", info);
    info = lapack_dgetri(nb, &*inv.data().begin(), nb, &*ipiv.begin(),
                         &*m.data().begin(), m.data().size());
    if (info)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "lapack error info = %d", info);

    MoFEMFunctionReturnHot(0);
  }
};

MoFEMErrorCode
OpSchurAssembleEnd<SchurDSYSV>::doWork(int side, EntityType type,
                                       EntitiesFieldData::EntData &data) {
  return doWorkImpl<SchurDSYSV>(side, type, data);
}

MoFEMErrorCode
OpSchurAssembleEnd<SchurDGESV>::doWork(int side, EntityType type,
                                       EntitiesFieldData::EntData &data) {
  return doWorkImpl<SchurDGESV>(side, type, data);
}

boost::shared_ptr<BlockStructure> createBlockMatStructure(

    DM dm,                                //< dm
    SchurFEOpsFEandFields schur_fe_op_vec //< block elements

) {

  auto cmp_uid_lo = [](const boost::weak_ptr<FieldEntity> &a, const UId &b) {
    if (auto a_ptr = a.lock()) {
      if (a_ptr->getLocalUniqueId() < b)
        return true;
      else
        return false;
    } else {
      return false;
    }
  };

  auto cmp_uid_hi = [](const UId &b, const boost::weak_ptr<FieldEntity> &a) {
    if (auto a_ptr = a.lock()) {
      if (b < a_ptr->getLocalUniqueId())
        return true;
      else
        return false;
    } else {
      return true;
    }
  };

  // get uids for fields
  auto get_uid_pair = [](const auto &field_id) {
    auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
        field_id, get_id_for_min_type<MBVERTEX>());
    auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
        field_id, get_id_for_max_type<MBENTITYSET>());
    return std::make_pair(lo_uid, hi_uid);
  };

  // get uids pair
  auto get_it_pair = [cmp_uid_lo, cmp_uid_hi](auto &&field_ents, auto &&p_uid) {
    auto lo = std::lower_bound(field_ents.begin(), field_ents.end(),
                               p_uid.first, cmp_uid_lo);
    auto hi = std::upper_bound(field_ents.begin(), field_ents.end(),
                               p_uid.second, cmp_uid_hi);
    return std::make_pair(lo, hi);
  };

  // extract DOFs for rows/columns. DOFs are associated with fields entities
  // for given problem.
  auto row_extractor = [](auto &e) { return e->entityCacheRowDofs; };
  auto col_extractor = [](auto &e) { return e->entityCacheColDofs; };

  auto extract_data = [](auto &&its, auto extractor) {
    std::vector<std::tuple<UId, int, int, int>> data;
    data.reserve(std::distance(its.first, its.second));
    // iterate field dofs
    for (; its.first != its.second; ++its.first) {
      if (auto e = its.first->lock()) {
        if (auto cache = extractor(e).lock()) {
          int nb_dofs = 0;
          for (auto dof = cache->loHi[0]; dof != cache->loHi[1]; ++dof) {
            if ((*dof)->getPetscGlobalDofIdx() != -1)
              ++nb_dofs;
          }
          auto uid = e->getLocalUniqueId();
          if (nb_dofs) {

            auto glob = (*cache->loHi[0])->getPetscGlobalDofIdx();
            auto loc = (*cache->loHi[0])->getPetscLocalDofIdx();
            while (glob == -1 && cache->loHi[0] != cache->loHi[1]) {
              ++cache->loHi[0];
              glob = (*cache->loHi[0])->getPetscGlobalDofIdx();
              loc = (*cache->loHi[0])->getPetscLocalDofIdx();
            }
            data.emplace_back(uid, glob, nb_dofs, loc);

#ifndef NDEBUG

            for (auto lo = cache->loHi[0]; lo != cache->loHi[1]; ++lo) {
              auto glob = (*lo)->getPetscGlobalDofIdx();
              if (glob == -1) {
                CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                                  "Wrong global index");
              }
            }

#endif
          } else {
            data.emplace_back(uid, -1, 0, -1);
          }
        }
      }
    }
    return data;
  };

  auto data_ptr = boost::make_shared<BlockStructure>();
  auto m_field_ptr = getInterfacePtr(dm);

  // create element to extract data
  auto fe_method = boost::shared_ptr<MoFEM::FEMethod>(new MoFEM::FEMethod());

  for (auto &d : schur_fe_op_vec) {

    // extract bit numbers  for row and column
    auto get_bit_numbers = [&d](auto op) {
      std::vector<FieldBitNumber> bit_numbers(d.second.size());
      std::transform(d.second.begin(), d.second.end(), bit_numbers.begin(), op);
      return bit_numbers;
    };

    // extract bit numbers  for row
    auto row_bit_numbers = get_bit_numbers(
        [&](auto &p) { return m_field_ptr->get_field_bit_number(p.first); });
    // extract bit numbers  for row
    auto col_bit_numbers = get_bit_numbers(
        [&](auto &p) { return m_field_ptr->get_field_bit_number(p.second); });

    fe_method->preProcessHook = []() { return 0; };
    fe_method->postProcessHook = []() { return 0; };
    fe_method->operatorHook = [&]() {
      MoFEMFunctionBegin;

      for (auto f = 0; f != row_bit_numbers.size(); ++f) {

        auto row_data =
            extract_data(get_it_pair(fe_method->getRowFieldEnts(),
                                     get_uid_pair(row_bit_numbers[f])),
                         row_extractor);
        auto col_data =
            extract_data(get_it_pair(fe_method->getColFieldEnts(),
                                     get_uid_pair(col_bit_numbers[f])),
                         col_extractor);

        for (auto &r : row_data) {
          auto [r_uid, r_glob, r_nb_dofs, r_loc] = r;
          for (auto &c : col_data) {
            auto [c_uid, c_glob, c_nb_dofs, c_loc] = c;
            data_ptr->blockIndex.insert(
                BlockStructure::Indexes{r_uid, c_uid, r_glob, c_glob, r_nb_dofs,
                                        c_nb_dofs, r_loc, c_loc, -1, -1});
          }
        }
      }

      MoFEMFunctionReturn(0);
    };

    CHKERR DMoFEMLoopFiniteElementsUpAndLowRank(dm, d.first, fe_method, 0,
                                                m_field_ptr->get_comm_size());
  };

  // order by column (that is for matrix multiplication)
  auto mem_size = 0;
  for (auto &v : data_ptr->blockIndex.get<0>()) {
    v.getMatShift() = mem_size;
    mem_size += v.getNbCols() * v.getNbRows();
  }

  std::vector<double> tmp;
  if (mem_size > tmp.max_size())
    CHK_THROW_MESSAGE(MOFEM_OPERATION_UNSUCCESSFUL, "Vector too big");

  data_ptr->dataBlocksPtr =
      boost::make_shared<std::vector<double>>(mem_size, 0);

  auto ghost_x = createDMVector(dm);
  auto ghost_y = createDMVector(dm);
  CHKERR VecSetDM(ghost_x, PETSC_NULL);
  CHKERR VecSetDM(ghost_y, PETSC_NULL);

  data_ptr->ghostX = ghost_x;
  data_ptr->ghostY = ghost_y;

  return data_ptr;
}

static MoFEMErrorCode mult_schur_block_shell(
    Mat mat, Vec x, Vec y, InsertMode iora,
    boost::function<
        int(DiagBlockIndex::BlockIndex::nth_index<0>::type::iterator)>
        shift_extractor,
    boost::shared_ptr<std::vector<double>> data_blocks_ptr,
    bool multiply_by_preconditioner);

static MoFEMErrorCode solve_schur_block_shell(Mat mat, Vec x, Vec y,
                                              InsertMode iora);

static PetscErrorCode mult(Mat mat, Vec x, Vec y) {
  BlockStructure *ctx;
  CHKERR MatShellGetContext(mat, (void **)&ctx);
  return mult_schur_block_shell(
      mat, x, y, INSERT_VALUES,

      [](DiagBlockIndex::BlockIndex::nth_index<0>::type::iterator it) {
        return it->getMatShift();
      },

      ctx->dataBlocksPtr, true);
}
static PetscErrorCode mult_add(Mat mat, Vec x, Vec y) {
  BlockStructure *ctx;
  CHKERR MatShellGetContext(mat, (void **)&ctx);
  return mult_schur_block_shell(
      mat, x, y, ADD_VALUES,

      [](DiagBlockIndex::BlockIndex::nth_index<0>::type::iterator it) {
        return it->getMatShift();
      },

      ctx->dataBlocksPtr, true);
}
static PetscErrorCode solve(Mat mat, Vec x, Vec y) {
  return solve_schur_block_shell(mat, x, y, INSERT_VALUES);
}
static PetscErrorCode solve_add(Mat mat, Vec x, Vec y) {
  return solve_schur_block_shell(mat, x, y, ADD_VALUES);
}

static PetscErrorCode zero_rows_columns(Mat A, PetscInt N,
                                        const PetscInt rows[], PetscScalar diag,
                                        Vec x, Vec b) {

  MoFEMFunctionBeginHot;
  BlockStructure *ctx;
  CHKERR MatShellGetContext(A, (void **)&ctx);

  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_zeroRowsAndCols, 0, 0, 0, 0);

  const int *ptr = &rows[0];
  int is_nb_rows = N;
  SmartPetscObj<IS> is_local;

  MPI_Comm comm;
  CHKERR PetscObjectGetComm((PetscObject)A, &comm);
  int size;
  MPI_Comm_size(comm, &size);
  if (size > 1) {
    auto is = createISGeneral(comm, N, rows, PETSC_USE_POINTER);
    is_local = isAllGather(is);
  }
  if (is_local) {
    CHKERR ISGetSize(is_local, &is_nb_rows);
#ifndef NDEBUG
    if constexpr (debug_schur) {
      CHKERR ISView(is_local, PETSC_VIEWER_STDOUT_WORLD);
    }
#endif
    CHKERR ISGetIndices(is_local, &ptr);
  }

  int loc_m, loc_n;
  CHKERR MatGetLocalSize(A, &loc_m, &loc_n);

  for (auto n = 0; n != is_nb_rows; ++n) {
    auto row = ptr[n];
    auto rlo = ctx->blockIndex.get<1>().lower_bound(row);
    auto rhi = ctx->blockIndex.get<1>().end();
    for (; rlo != rhi; ++rlo) {
      auto r_shift = row - rlo->getRow();
      if (r_shift >= 0 && r_shift < rlo->getNbRows()) {
        auto *ptr = &(*ctx->dataBlocksPtr)[rlo->getMatShift()];
        for (auto i = 0; i != rlo->getNbCols(); ++i) {
          ptr[i + r_shift * rlo->getNbCols()] = 0;
        }
      } else if (rlo->getRow() + rlo->getNbRows() > row) {
        break;
      }
    }
  }

  for (auto n = 0; n != is_nb_rows; ++n) {
    auto col = ptr[n];
    auto clo = ctx->blockIndex.get<2>().lower_bound(col);
    auto chi = ctx->blockIndex.get<2>().end();
    for (; clo != chi; ++clo) {
      auto c_shift = col - clo->getCol();
      if (c_shift >= 0 && c_shift < clo->getNbCols()) {

        auto *ptr = &(*ctx->dataBlocksPtr)[clo->getMatShift()];
        for (auto i = 0; i != clo->getNbRows(); ++i) {
          ptr[c_shift + i * clo->getNbCols()] = 0;
        }

        // diagonal
        if (

            clo->getRow() == clo->getCol() && clo->getLocRow() < loc_m &&
            clo->getLocCol() < loc_n

        ) {
          auto r_shift = col - clo->getCol();
          if (r_shift >= 0 && r_shift < clo->getNbRows()) {
            ptr[c_shift + r_shift * clo->getNbCols()] = diag;
          }
        }
      } else if (clo->getCol() + clo->getNbCols() > col) {
        break;
      }
    }
  }

  if (is_local) {
    CHKERR ISRestoreIndices(is_local, &ptr);
  }

  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_zeroRowsAndCols, 0, 0, 0, 0);

  MoFEMFunctionReturnHot(0);
}

static PetscErrorCode mat_zero(Mat m) {
  MoFEMFunctionBegin;
  BlockStructure *ctx;
  CHKERR MatShellGetContext(m, (void **)&ctx);
  if (ctx->dataBlocksPtr)
    std::fill(ctx->dataBlocksPtr->begin(), ctx->dataBlocksPtr->end(), 0.);
  if (ctx->dataInvBlocksPtr)
    std::fill(ctx->dataInvBlocksPtr->begin(), ctx->dataInvBlocksPtr->end(), 0.);
  if (ctx->preconditionerBlocksPtr)
    std::fill(ctx->preconditionerBlocksPtr->begin(),
              ctx->preconditionerBlocksPtr->end(), 0.);
  MoFEMFunctionReturn(0);
}

static MoFEMErrorCode setSchurBlockMatOps(Mat mat_raw) {
  MoFEMFunctionBegin;
  CHKERR MatShellSetManageScalingShifts(mat_raw);
  CHKERR MatShellSetOperation(mat_raw, MATOP_MULT, (void (*)(void))mult);
  CHKERR MatShellSetOperation(mat_raw, MATOP_MULT_ADD,
                              (void (*)(void))mult_add);
  CHKERR MatShellSetOperation(mat_raw, MATOP_SOLVE, (void (*)(void))solve);
  CHKERR MatShellSetOperation(mat_raw, MATOP_SOLVE_ADD,
                              (void (*)(void))solve_add);
  CHKERR MatShellSetOperation(mat_raw, MATOP_ZERO_ENTRIES,
                              (void (*)(void))mat_zero);
  CHKERR MatShellSetOperation(mat_raw, MATOP_ZERO_ROWS_COLUMNS,
                              (void (*)(void))zero_rows_columns);

  MoFEMFunctionReturn(0);
};

SchurShellMatData createBlockMat(DM dm,
                                 boost::shared_ptr<BlockStructure> data) {

  auto problem_ptr = getProblemPtr(dm);
  auto nb_local = problem_ptr->nbLocDofsRow;
  auto nb_global = problem_ptr->nbDofsRow;

  // check in nb, rows is equal to nb. columns.
  if (nb_local != problem_ptr->nbLocDofsCol) {
    MOFEM_LOG("SELF", Sev::error)
        << "Wrong size " << nb_local << " != " << problem_ptr->nbLocDofsCol;
    CHK_MOAB_THROW(MOFEM_DATA_INCONSISTENCY,
                   "nb. cols is inconsistent with nb. rows");
  }
  if (nb_global != problem_ptr->nbDofsCol) {
    MOFEM_LOG("SELF", Sev::error)
        << "Wrong size " << nb_global << " != " << problem_ptr->nbDofsCol;
    CHK_MOAB_THROW(MOFEM_DATA_INCONSISTENCY,
                   "nb. cols is inconsistent with nb. rows");
  }

  // get comm from DM
  MPI_Comm comm;
  CHKERR PetscObjectGetComm((PetscObject)dm, &comm);

  Mat mat_raw;
  CHKERR MatCreateShell(comm, nb_local, nb_local, nb_global, nb_global,
                        (void *)data.get(), &mat_raw);
  CHKERR setSchurBlockMatOps(mat_raw);
  // CHKERR PetscObjectSetName((PetscObject)mat_raw, MoFEM_BLOCK_MAT);

  return std::make_pair(SmartPetscObj<Mat>(mat_raw), data);
}

constexpr int max_gemv_size = 2;

static MoFEMErrorCode mult_schur_block_shell(
    Mat mat, Vec x, Vec y, InsertMode iora,

    boost::function<
        int(DiagBlockIndex::BlockIndex::nth_index<0>::type::iterator)>
        shift_extractor,

    boost::shared_ptr<std::vector<double>> data_blocks_ptr,
    bool multiply_by_preconditioner) {
  MoFEMFunctionBegin;
  BlockStructure *ctx;
  CHKERR MatShellGetContext(mat, (void **)&ctx);

  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_BlockStructureMult, 0, 0, 0, 0);

  int x_loc_size;
  CHKERR VecGetLocalSize(x, &x_loc_size);
  int y_loc_size;
  CHKERR VecGetLocalSize(y, &y_loc_size);

  Vec ghost_x = ctx->ghostX;
  Vec ghost_y = ctx->ghostY;

  CHKERR VecCopy(x, ghost_x);

  double *y_array;
  Vec loc_ghost_y;
  CHKERR VecGhostGetLocalForm(ghost_y, &loc_ghost_y);
  int nb_y;
  CHKERR VecGetLocalSize(loc_ghost_y, &nb_y);
  CHKERR VecGetArray(loc_ghost_y, &y_array);
  for (auto i = 0; i != nb_y; ++i)
    y_array[i] = 0.;
  CHKERR VecRestoreArray(loc_ghost_y, &y_array);
  CHKERR VecGhostRestoreLocalForm(ghost_y, &loc_ghost_y);

  auto mult = [&](int low_x, int hi_x, int low_y, int hi_y) {
    MoFEMFunctionBegin;

    double *x_array;
    Vec loc_ghost_x;
    CHKERR VecGhostGetLocalForm(ghost_x, &loc_ghost_x);
    CHKERR VecGetArray(loc_ghost_x, &x_array);

    double *y_array;
    Vec loc_ghost_y;
    CHKERR VecGhostGetLocalForm(ghost_y, &loc_ghost_y);
    int nb_y;
    CHKERR VecGetLocalSize(loc_ghost_y, &nb_y);
    CHKERR VecGetArray(loc_ghost_y, &y_array);

    double *block_ptr = &*data_blocks_ptr->begin();
    auto it = ctx->blockIndex.get<0>().begin();
    auto hi = ctx->blockIndex.get<0>().end();

    while (it != hi) {
      if (it->getLocRow() < low_y || it->getLocRow() >= hi_y ||
          it->getLocCol() < low_x || it->getLocCol() >= hi_x) {
        ++it;
        continue;
      }

      auto nb_rows = it->getNbRows();
      auto nb_cols = it->getNbCols();
      auto x_ptr = &x_array[it->getLocCol()];
      auto y_ptr = &y_array[it->getLocRow()];
      auto ptr = &block_ptr[shift_extractor(it)];

      if (std::min(nb_rows, nb_cols) > max_gemv_size) {
        cblas_dgemv(CblasRowMajor, CblasNoTrans, nb_rows, nb_cols, 1.0, ptr,
                    nb_cols, x_ptr, 1, 1.0, y_ptr, 1);
      } else {
        for (auto r = 0; r != nb_rows; ++r) {
          for (auto c = 0; c != nb_cols; ++c) {
            y_ptr[r] += ptr[r * nb_cols + c] * x_ptr[c];
          }
        }
      }
      ++it;
    }

    if (multiply_by_preconditioner && ctx->multiplyByPreconditioner) {

      if (!ctx->preconditionerBlocksPtr)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "No parentBlockStructurePtr");

      auto preconditioner_ptr = &*ctx->preconditionerBlocksPtr->begin();

      auto it = ctx->blockIndex.get<0>().begin();
      auto hi = ctx->blockIndex.get<0>().end();

      while (it != hi) {
        if (it->getLocRow() < low_y || it->getLocRow() >= hi_y ||
            it->getLocCol() < low_x || it->getLocCol() >= hi_x) {
          ++it;
          continue;
        }

        if (it->getInvShift() != -1) {
          auto nb_rows = it->getNbRows();
          auto nb_cols = it->getNbCols();
          auto x_ptr = &x_array[it->getLocCol()];
          auto y_ptr = &y_array[it->getLocRow()];
          auto ptr = &preconditioner_ptr[it->getInvShift()];
          if (std::min(nb_rows, nb_cols) > max_gemv_size) {
            cblas_dgemv(CblasRowMajor, CblasNoTrans, nb_rows, nb_cols, 1.0, ptr,
                        nb_cols, x_ptr, 1, 1.0, y_ptr, 1);
          } else {
            for (auto r = 0; r != nb_rows; ++r) {
              for (auto c = 0; c != nb_cols; ++c) {
                y_ptr[r] += ptr[r * nb_cols + c] * x_ptr[c];
              }
            }
          }
        }

        ++it;
      }
    }

    CHKERR VecRestoreArray(loc_ghost_x, &x_array);
    CHKERR VecRestoreArray(loc_ghost_y, &y_array);
    CHKERR VecGhostRestoreLocalForm(ghost_x, &loc_ghost_x);
    CHKERR VecGhostRestoreLocalForm(ghost_y, &loc_ghost_y);
    MoFEMFunctionReturn(0);
  };

  constexpr auto max_int = std::numeric_limits<int>::max();
  CHKERR VecGhostUpdateBegin(ghost_x, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR mult(0, x_loc_size, 0, max_int);
  CHKERR VecGhostUpdateEnd(ghost_x, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR mult(x_loc_size, max_int, 0, max_int);

  CHKERR VecGhostUpdateBegin(ghost_y, ADD_VALUES, SCATTER_REVERSE);
  CHKERR VecGhostUpdateEnd(ghost_y, ADD_VALUES, SCATTER_REVERSE);

  switch (iora) {
  case INSERT_VALUES:
    CHKERR VecCopy(ghost_y, y);
    break;
  case ADD_VALUES:
    CHKERR VecAXPY(y, 1., ghost_y);
    break;
  default:
    CHK_MOAB_THROW(MOFEM_NOT_IMPLEMENTED, "Wrong InsertMode");
  }

#ifndef NDEBUG

  auto print_norm = [&](auto msg, auto y) {
    MoFEMFunctionBegin;
    PetscReal norm;
    CHKERR VecNorm(y, NORM_2, &norm);
    int nb_loc_y;
    CHKERR VecGetLocalSize(y, &nb_loc_y);
    int nb_y;
    CHKERR VecGetSize(y, &nb_y);
    MOFEM_LOG("WORLD", Sev::noisy)
        << msg << " " << nb_y << " " << nb_loc_y << " norm " << norm;
    MoFEMFunctionReturn(0);
  };

  switch (iora) {
  case INSERT_VALUES:
    print_norm("mult_schur_block_shell insert x", x);
    print_norm("mult_schur_block_shell insert y", y);
    break;
  case ADD_VALUES:
    print_norm("mult_schur_block_shell add x", x);
    print_norm("mult_schur_block_shell add y", y);
    break;
  default:
    CHK_MOAB_THROW(MOFEM_NOT_IMPLEMENTED, "Wrong InsertMode");
  }

#endif // NDEBUG

  // PetscLogFlops(xxx)
  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_BlockStructureMult, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

static MoFEMErrorCode solve_schur_block_shell(Mat mat, Vec y, Vec x,
                                              InsertMode iora) {
  MoFEMFunctionBegin;
  BlockStructure *ctx;
  CHKERR MatShellGetContext(mat, (void **)&ctx);

  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_BlockStructureSolve, 0, 0, 0, 0);

  if (!ctx->dataInvBlocksPtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "No dataInvBlocksPtr");
  CHKERR mult_schur_block_shell(
      mat, y, x, iora,

      [](DiagBlockIndex::BlockIndex::nth_index<0>::type::iterator it) {
        return it->getInvShift();
      },

      ctx->dataInvBlocksPtr, false);

  // PetscLogFlops(xxx)
  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_BlockStructureSolve, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

inline MoFEMErrorCode shell_block_mat_asmb_wrap_impl(
    BlockStructure *ctx, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora,
    boost::function<int(const DiagBlockIndex::Indexes *)> shift_extractor,
    boost::shared_ptr<std::vector<double>> data_blocks_ptr) {

  MoFEMFunctionBegin;

  if (row_data.getIndices().empty())
    MoFEMFunctionReturnHot(0);
  if (col_data.getIndices().empty())
    MoFEMFunctionReturnHot(0);

#ifndef NDEBUG

  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_BlockStructureSetValues, 0, 0, 0,
                     0);

#endif // NDEBUG

  switch (iora) {
  case ADD_VALUES:
  case INSERT_VALUES:
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Wrong InsertMode");
  }

  auto get_row_data = [&]() -> std::pair<bool, VectorInt> {
    if (auto e_ptr = row_data.getFieldEntities()[0]) {
      if (auto stored_data_ptr =
              e_ptr->getSharedStoragePtr<EssentialBcStorage>()) {
        MOFEM_LOG("SELF", Sev::warning) << "Can lead to unhandled behaviour";
        return std::make_pair(true, stored_data_ptr->entityIndices);
      }
    }
    return std::make_pair(false, row_data.getIndices());
  };

  auto row_indices = get_row_data();

  std::vector<int> ent_row_indices;
  std::vector<int> ent_col_indices;

  for (auto &rent : row_data.getFieldEntities()) {
    if (auto r_cache = rent->entityCacheRowDofs.lock()) {

      auto &row_uid = rent->getLocalUniqueId();
      auto &row_ind = row_indices.second;

      for (auto &cent : col_data.getFieldEntities()) {
        if (auto c_cache = cent->entityCacheColDofs.lock()) {

          auto &col_uid = cent->getLocalUniqueId();
          auto &col_ind = col_data.getIndices();

          auto it = ctx->blockIndex.get<0>().find(

              boost::make_tuple(row_uid, col_uid)

          );

#ifndef NDEBUG

          if (it == ctx->blockIndex.get<0>().end()) {
            MOFEM_LOG_CHANNEL("SELF");
            MOFEM_TAG_AND_LOG("SELF", Sev::error, "BlockMat")
                << "missing block: row "
                << row_data.getFieldDofs()[0]->getName() << " col "
                << col_data.getFieldDofs()[0]->getName();
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "Block not allocated");
          }

#endif

          auto ri = row_ind.begin();
          auto rie = row_ind.end();

          auto shift = shift_extractor(&*it);
          auto s_mat = &(*data_blocks_ptr)[shift];

          auto get_ent_indices = [](auto &cache, auto &ind) {
            ind.clear();
            ind.reserve(std::distance(cache->loHi[0], cache->loHi[1]));
            for (auto e = cache->loHi[0]; e != cache->loHi[1]; ++e) {
              auto glob = (*e)->getPetscGlobalDofIdx();
              if (glob != -1)
                ind.push_back(glob);
            }
          };

          get_ent_indices(r_cache, ent_row_indices);
          if (ent_row_indices.empty())
            continue;
          get_ent_indices(c_cache, ent_col_indices);
          if (ent_col_indices.empty())
            continue;

          if (mat.size1() == ent_row_indices.size() &&
              mat.size2() == ent_col_indices.size()) {

            if (iora == ADD_VALUES) {
              cblas_daxpy(mat.data().size(), 1.0, mat.data().data(), 1, s_mat,
                          1);
            } else {
              cblas_dcopy(mat.data().size(), mat.data().data(), 1, s_mat, 1);
            }

          } else {

            int row = 0;
            for (auto re : ent_row_indices) {
              ri = std::find(ri, rie, re);
              if (!(ri == rie && *ri != -1)) {

                auto ci = col_ind.begin();
                auto cie = col_ind.end();
                auto ce = c_cache->loHi[0];

                int col = 0;
                for (auto ce : ent_col_indices) {
                  ci = std::find(ci, cie, ce);
                  if (!(ci == cie && *ci != -1)) {
                    auto &m = s_mat[row * ent_col_indices.size() + col];
                    if (iora == ADD_VALUES) {
                      m += mat(std::distance(row_ind.begin(), ri),
                               std::distance(col_ind.begin(), ci));
                    } else {
                      m = mat(std::distance(row_ind.begin(), ri),
                              std::distance(col_ind.begin(), ci));
                    }
                  }
                  ++col;
                } // cols
              }
              ++row;
            } // rows
          }

        } // iterate entities
      }
    }
  }

#ifndef NDEBUG
  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_BlockStructureSetValues, 0, 0, 0,
                   0);
#endif // NDEBUG

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
shell_block_mat_asmb_wrap(Mat M, const EntitiesFieldData::EntData &row_data,
                          const EntitiesFieldData::EntData &col_data,
                          const MatrixDouble &mat, InsertMode iora) {
  MoFEMFunctionBegin;
  PetscBool is_mat_shell = PETSC_FALSE;
  PetscObjectTypeCompare((PetscObject)M, MATSHELL, &is_mat_shell);
  if (is_mat_shell) {
    BlockStructure *ctx;
    CHKERR MatShellGetContext(M, (void **)&ctx);
    CHKERR shell_block_mat_asmb_wrap_impl(
        ctx, row_data, col_data, mat, iora,
        [](const DiagBlockIndex::Indexes *idx) { return idx->getMatShift(); },
        ctx->dataBlocksPtr);
  } else {
    CHKERR MatSetValues<AssemblyTypeSelector<PETSC>>(M, row_data, col_data, mat,
                                                     iora);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode shell_block_preconditioner_mat_asmb_wrap(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora) {
  MoFEMFunctionBegin;
  PetscBool is_mat_shell = PETSC_FALSE;
  PetscObjectTypeCompare((PetscObject)M, MATSHELL, &is_mat_shell);
  if (is_mat_shell) {
    BlockStructure *ctx;
    CHKERR MatShellGetContext(M, (void **)&ctx);
    if (!ctx->preconditionerBlocksPtr)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "No preconditionerBlocksPtr");
    CHKERR shell_block_mat_asmb_wrap_impl(
        ctx, row_data, col_data, mat, iora,
        [](const DiagBlockIndex::Indexes *idx) { return idx->getInvShift(); },
        ctx->preconditionerBlocksPtr);
  } else {
    CHKERR MatSetValues<AssemblyTypeSelector<PETSC>>(M, row_data, col_data, mat,
                                                     iora);
  }
  MoFEMFunctionReturn(0);
}

boost::shared_ptr<NestSchurData> getNestSchurData(

    std::pair<SmartPetscObj<DM>, SmartPetscObj<DM>> dms,
    boost::shared_ptr<BlockStructure> block_mat_data_ptr,

    std::vector<std::string> fields_names, //< a00 fields
    std::vector<boost::shared_ptr<Range>>
        field_ents, //< a00 ranges (can be null),
    bool add_preconditioner_block) {

  if (!block_mat_data_ptr) {
    CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY, "Block data not set");
  }

  if (fields_names.size() != field_ents.size())
    CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                      "fields_names.size() != field_ents.size()");

  auto [schur_dm, block_dm] = dms;
  auto schur_prb = getProblemPtr(schur_dm);
  auto block_prb = getProblemPtr(block_dm);
  // auto m_field_ptr = getInterfacePtr(block_dm);

  auto schur_dofs_row = schur_prb->getNumeredRowDofsPtr();
  auto schur_dofs_col = schur_prb->getNumeredColDofsPtr();
  auto block_dofs_row = block_prb->getNumeredRowDofsPtr();
  auto block_dofs_col = block_prb->getNumeredColDofsPtr();

  auto ao_schur_row = schur_prb->getSubData()->getSmartRowMap();
  auto ao_schur_col = schur_prb->getSubData()->getSmartColMap();
  auto ao_block_row = block_prb->getSubData()->getSmartRowMap();
  auto ao_block_col = block_prb->getSubData()->getSmartColMap();

  auto schur_vec_x = createDMVector(schur_dm);
  auto block_vec_x = createDMVector(block_dm);
  auto schur_vec_y = vectorDuplicate(schur_vec_x);
  auto block_vec_y = vectorDuplicate(block_vec_x);
  CHKERR VecSetDM(schur_vec_x, PETSC_NULL);
  CHKERR VecSetDM(block_vec_x, PETSC_NULL);
  CHKERR VecSetDM(schur_vec_y, PETSC_NULL);
  CHKERR VecSetDM(block_vec_y, PETSC_NULL);

  auto find_field_ent = [&](auto uid, auto prb, auto rc) {
    boost::shared_ptr<NumeredDofEntity_multiIndex> dofs;

    switch (rc) {
    case ROW:
      dofs = prb->getNumeredRowDofsPtr();
      break;
    case COL:
      dofs = prb->getNumeredColDofsPtr();
      break;
    default:
      CHK_MOAB_THROW(MOFEM_NOT_IMPLEMENTED, "Wrong RowCol");
      break;
    }

    auto lo = dofs->get<Unique_mi_tag>().lower_bound(uid);
    if (lo == dofs->get<Unique_mi_tag>().end())
      return boost::shared_ptr<NumeredDofEntity>();
    auto hi = dofs->get<Unique_mi_tag>().upper_bound(
        DofEntity::getUniqueIdCalculate(MAX_DOFS_ON_ENTITY - 1, uid));
    if (lo != hi)
      return *lo;

    return boost::shared_ptr<NumeredDofEntity>();
  };

  std::array<boost::shared_ptr<BlockStructure>, 4> data_ptrs;

  for (auto r = 0; r != 3; ++r) {
    data_ptrs[r] = boost::make_shared<BlockStructure>();
    data_ptrs[r]->dataBlocksPtr = block_mat_data_ptr->dataBlocksPtr;
  }
  data_ptrs[3] = boost::make_shared<BlockStructure>();
  data_ptrs[3]->dataBlocksPtr = block_mat_data_ptr->dataBlocksPtr;

  data_ptrs[0]->ghostX = schur_vec_x;
  data_ptrs[0]->ghostY = schur_vec_y;
  data_ptrs[1]->ghostX = block_vec_x;
  data_ptrs[1]->ghostY = schur_vec_y;
  data_ptrs[2]->ghostX = schur_vec_x;
  data_ptrs[2]->ghostY = block_vec_y;
  data_ptrs[3]->ghostX = block_vec_x;
  data_ptrs[3]->ghostY = block_vec_y;

  int idx = 0;
  for (auto &d : block_mat_data_ptr->blockIndex.get<0>()) {

    auto insert = [&](auto &m, auto &dof_r, auto &dof_c, auto &d) {
      m.insert(

          BlockStructure::Indexes{
              d.getRowUId(), d.getColUId(),

              dof_r->getPetscGlobalDofIdx(), dof_c->getPetscGlobalDofIdx(),

              d.getNbRows(), d.getNbCols(),

              dof_r->getPetscLocalDofIdx(), dof_c->getPetscLocalDofIdx(),

              d.getMatShift(), d.getInvShift()}

      );
    };

    auto dof_schur_row = find_field_ent(d.getRowUId(), schur_prb, ROW);
    auto dof_schur_col = find_field_ent(d.getColUId(), schur_prb, COL);
    auto dof_block_row = find_field_ent(d.getRowUId(), block_prb, ROW);
    auto dof_block_col = find_field_ent(d.getColUId(), block_prb, COL);

    if (dof_schur_row && dof_schur_col) {
      insert(data_ptrs[0]->blockIndex, dof_schur_row, dof_schur_col, d);
    }

    if (dof_schur_row && dof_block_col) {
      insert(data_ptrs[1]->blockIndex, dof_schur_row, dof_block_col, d);
    }

    if (dof_block_row && dof_schur_col) {
      insert(data_ptrs[2]->blockIndex, dof_block_row, dof_schur_col, d);
    }

    if (dof_block_row && dof_block_col) {
      insert(data_ptrs[3]->blockIndex, dof_block_row, dof_block_col, d);
    }

    ++idx;
  }

  // set data for a00 solve (inverse blocks)
  auto set_up_a00_data = [&](auto inv_block_data) {
    MoFEMFunctionBegin;

    auto size = inv_block_data->dataBlocksPtr->size();
    inv_block_data->dataInvBlocksPtr =
        boost::make_shared<std::vector<double>>(size, 0.0);
    block_mat_data_ptr->dataInvBlocksPtr = inv_block_data->dataInvBlocksPtr;

    for (auto &s : inv_block_data->blockIndex) {
      s.getInvShift() = s.getMatShift();
    }

    for (auto &s : inv_block_data->blockIndex) {
      auto it = block_mat_data_ptr->blockIndex.find(

          boost::make_tuple(s.getRowUId(), s.getColUId())

      );
      if (it == block_mat_data_ptr->blockIndex.end())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Missing block");
      it->getInvShift() = s.getInvShift();
    }

    if (add_preconditioner_block) {
      auto preconditioned_block = boost::make_shared<std::vector<double>>(
          inv_block_data->dataInvBlocksPtr->size(), 0);
      inv_block_data->preconditionerBlocksPtr = preconditioned_block;
      inv_block_data->multiplyByPreconditioner = true;
      block_mat_data_ptr->preconditionerBlocksPtr =
          inv_block_data->preconditionerBlocksPtr;
      block_mat_data_ptr->multiplyByPreconditioner = false;
    }

    MoFEMFunctionReturn(0);
  };

  CHKERR set_up_a00_data(data_ptrs[3]);

  MPI_Comm comm;
  CHKERR PetscObjectGetComm((PetscObject)schur_dm, &comm);

  auto create_shell_mat = [&](auto nb_r_loc, auto nb_c_loc, auto nb_r_glob,
                              auto nb_c_glob, auto data_ptr) {
    Mat mat_raw;
    CHKERR MatCreateShell(comm, nb_r_loc, nb_c_loc, nb_r_glob, nb_c_glob,
                          (void *)data_ptr.get(), &mat_raw);
    CHKERR setSchurBlockMatOps(mat_raw);
    return SmartPetscObj<Mat>(mat_raw);
  };

  auto schur_nb_global = schur_prb->getNbDofsRow();
  auto block_nb_global = block_prb->getNbDofsRow();
  auto schur_nb_local = schur_prb->getNbLocalDofsRow();
  auto block_nb_local = block_prb->getNbLocalDofsRow();

  std::array<SmartPetscObj<Mat>, 4> mats_array;
  mats_array[0] =
      create_shell_mat(schur_nb_local, schur_nb_local, schur_nb_global,
                       schur_nb_global, data_ptrs[0]);
  mats_array[1] =
      create_shell_mat(schur_nb_local, block_nb_local, schur_nb_global,
                       block_nb_global, data_ptrs[1]);
  mats_array[2] =
      create_shell_mat(block_nb_local, schur_nb_local, block_nb_global,
                       schur_nb_global, data_ptrs[2]);
  mats_array[3] =
      create_shell_mat(block_nb_local, block_nb_local, block_nb_global,
                       block_nb_global, data_ptrs[3]);

  MOFEM_TAG_AND_LOG("SYNC", Sev::verbose, "NestedSchur")
      << "(0, 0) " << schur_nb_local << " " << schur_nb_global << " "
      << data_ptrs[0]->blockIndex.size();
  MOFEM_TAG_AND_LOG("SYNC", Sev::verbose, "NestedSchur")
      << "(0, 1) " << schur_nb_local << " " << block_nb_local << " "
      << schur_nb_global << " " << block_nb_global << " "
      << data_ptrs[1]->blockIndex.size();
  MOFEM_TAG_AND_LOG("SYNC", Sev::verbose, "NestedSchur")
      << "(1, 0) " << block_nb_local << " " << schur_nb_local << " "
      << block_nb_global << " " << schur_nb_global << " "
      << data_ptrs[2]->blockIndex.size();
  MOFEM_TAG_AND_LOG("SYNC", Sev::verbose, "NestedSchur")
      << "(1, 1) " << block_nb_local << " " << block_nb_global << " "
      << data_ptrs[3]->blockIndex.size();

  MOFEM_LOG_SEVERITY_SYNC(comm, Sev::verbose);

  auto schur_is = schur_prb->getSubData()->getSmartRowIs();
  auto block_is = block_prb->getSubData()->getSmartRowIs();

  return boost::make_shared<NestSchurData>(

      mats_array, data_ptrs, block_mat_data_ptr,
      std::make_pair(schur_is, block_is)

  );
}

std::pair<SmartPetscObj<Mat>, boost::shared_ptr<NestSchurData>>
createSchurNestedMatrix(boost::shared_ptr<NestSchurData> schur_net_data_ptr) {

  if (!schur_net_data_ptr)
    CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY, "No data");

  auto [mat_arrays, data_ptrs, block_mat_data_ptr, iss] = *schur_net_data_ptr;
  auto [schur_is, block_is] = iss;

  std::array<IS, 2> is_a = {schur_is, block_is};
  std::array<Mat, 4> mats_a = {

      mat_arrays[0], mat_arrays[1], mat_arrays[2], mat_arrays[3]

  };

  MPI_Comm comm;
  CHKERR PetscObjectGetComm((PetscObject)mat_arrays[0], &comm);

  Mat mat_raw;
  CHKERR MatCreateNest(

      comm, 2, is_a.data(), 2, is_a.data(), mats_a.data(), &mat_raw

  );

  return std::make_pair(SmartPetscObj<Mat>(mat_raw), schur_net_data_ptr);
}

OpSchurAssembleBase *createOpSchurAssembleBegin() {
  return new OpSchurAssembleBegin();
}

OpSchurAssembleBase *
createOpSchurAssembleEnd(std::vector<std::string> fields_name,
                         std::vector<boost::shared_ptr<Range>> field_ents,
                         SmartPetscObj<AO> ao, SmartPetscObj<Mat> schur,
                         bool sym_schur, bool symm_op,
                         boost::shared_ptr<BlockStructure> diag_blocks) {
  if (symm_op)
    return new OpSchurAssembleEnd<SchurDSYSV>(
        fields_name, field_ents, ao, schur, sym_schur, symm_op, diag_blocks);
  else
    return new OpSchurAssembleEnd<SchurDGESV>(
        fields_name, field_ents, ao, schur, sym_schur, symm_op, diag_blocks);
}

OpSchurAssembleBase *
createOpSchurAssembleEnd(std::vector<std::string> fields_name,
                         std::vector<boost::shared_ptr<Range>> field_ents,
                         std::vector<SmartPetscObj<AO>> sequence_of_aos,
                         std::vector<SmartPetscObj<Mat>> sequence_of_mats,
                         std::vector<bool> sym_schur, bool symm_op,
                         boost::shared_ptr<BlockStructure> diag_blocks) {
  return createOpSchurAssembleEnd(
      fields_name, field_ents, sequence_of_aos.back(), sequence_of_mats.back(),
      sym_schur.back(), symm_op, diag_blocks);
}

OpSchurAssembleBase *
createOpSchurAssembleEnd(std::vector<std::string> fields_name,
                         std::vector<boost::shared_ptr<Range>> field_ents,
                         std::vector<SmartPetscObj<AO>> sequence_of_aos,
                         std::vector<SmartPetscObj<Mat>> sequence_of_mats,
                         std::vector<bool> sym_schur,
                         std::vector<double> diag_eps, bool symm_op,
                         boost::shared_ptr<BlockStructure> diag_blocks) {
  return createOpSchurAssembleEnd(
      fields_name, field_ents, sequence_of_aos.back(), sequence_of_mats.back(),
      sym_schur.back(), symm_op, diag_blocks);
}

MoFEMErrorCode setSchurA00MatSolvePC(SmartPetscObj<PC> pc) {
  MoFEMFunctionBegin;

  auto apply = [](PC pc, Vec f, Vec x) {
    MoFEMFunctionBeginHot;
    Mat A, P;
    CHKERR PCGetOperators(pc, &A, &P);
    CHKERR MatSolve(P, f, x);
    MoFEMFunctionReturnHot(0);
  };

  CHKERR PCSetType(pc, PCSHELL);
  CHKERR PCShellSetApply(pc, apply);
  CHKERR PCShellSetName(pc, "MoFEMSchurBlockPC");

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode setSchurMatSolvePC(SmartPetscObj<PC> pc) {
  return setSchurA00MatSolvePC(pc);
}

// This is used to assemble local element matrices for Schur complement
// and matrix for KSP
template <>
MoFEMErrorCode
MatSetValues<SchurElemMats>(Mat M, const EntitiesFieldData::EntData &row_data,
                            const EntitiesFieldData::EntData &col_data,
                            const MatrixDouble &mat, InsertMode iora) {
  return SchurElemMats::MatSetValues(M, row_data, col_data, mat, iora);
}

MoFEMErrorCode
schur_mat_set_values_wrap(Mat M, const EntitiesFieldData::EntData &row_data,
                          const EntitiesFieldData::EntData &col_data,
                          const MatrixDouble &mat, InsertMode iora) {
  return MatSetValues<AssemblyTypeSelector<PETSC>>(M, row_data, col_data, mat,
                                                   iora);
}

// Is assumed that standard PETSc assembly works for matrices used by KSP
SchurBackendMatSetValuesPtr::MatSetValuesPtr
    SchurBackendMatSetValuesPtr::matSetValuesPtr = schur_mat_set_values_wrap;

// We assemble matrix for KSP and store local matrices for Schur complement.
// Schur complement is calculated and assembled in OpSchurAssembleEnd.
MoFEMErrorCode
SchurElemMats::MatSetValues(Mat M, const EntitiesFieldData::EntData &row_data,
                            const EntitiesFieldData::EntData &col_data,
                            const MatrixDouble &mat, InsertMode iora) {
  MoFEMFunctionBegin;
  CHKERR assembleStorage(row_data, col_data, mat, iora);
  CHKERR SchurBackendMatSetValuesPtr::matSetValuesPtr(M, row_data, col_data,
                                                      mat, iora);
  MoFEMFunctionReturn(0);
}

// All is now wrapped in specialisation of
// MatSetValues<AssemblyTypeSelector<SCHUR>>
template <>
MoFEMErrorCode MatSetValues<AssemblyTypeSelector<SCHUR>>(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora) {
  return MatSetValues<SchurElemMats>(M, row_data, col_data, mat, iora);
}

// Standard (PETSc) assembly block matrices do noe work
template <>
MoFEMErrorCode
MatSetValues<BlockStructure>(Mat M, const EntitiesFieldData::EntData &row_data,
                             const EntitiesFieldData::EntData &col_data,
                             const MatrixDouble &mat, InsertMode iora) {
  return shell_block_mat_asmb_wrap(M, row_data, col_data, mat, iora);
}

struct SchurElemMatsBlock : public SchurElemMats {

  static MoFEMErrorCode MatSetValues(Mat M,
                                     const EntitiesFieldData::EntData &row_data,
                                     const EntitiesFieldData::EntData &col_data,
                                     const MatrixDouble &mat, InsertMode iora);
};

SchurBackendMatSetValuesPtr::MatSetValuesPtr
    SchurBackendMatSetValuesPtr::matSetValuesBlockPtr =
        shell_block_mat_asmb_wrap;

MoFEMErrorCode
SchurElemMatsBlock::MatSetValues(Mat M,
                                 const EntitiesFieldData::EntData &row_data,
                                 const EntitiesFieldData::EntData &col_data,
                                 const MatrixDouble &mat, InsertMode iora) {
  MoFEMFunctionBegin;
  CHKERR assembleStorage(row_data, col_data, mat, iora);
  CHKERR SchurBackendMatSetValuesPtr::matSetValuesBlockPtr(M, row_data,
                                                           col_data, mat, iora);
  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode
MatSetValues<SchurElemMatsBlock>(Mat M,
                                 const EntitiesFieldData::EntData &row_data,
                                 const EntitiesFieldData::EntData &col_data,
                                 const MatrixDouble &mat, InsertMode iora) {
  return SchurElemMatsBlock::MatSetValues(M, row_data, col_data, mat, iora);
}

struct SchurElemMatsPreconditionedBlock : public SchurElemMats {

  static MoFEMErrorCode MatSetValues(Mat M,
                                     const EntitiesFieldData::EntData &row_data,
                                     const EntitiesFieldData::EntData &col_data,
                                     const MatrixDouble &mat, InsertMode iora);
};

SchurBackendMatSetValuesPtr::MatSetValuesPtr
    SchurBackendMatSetValuesPtr::matSetValuesPreconditionedBlockPtr =
        shell_block_preconditioner_mat_asmb_wrap;

MoFEMErrorCode SchurElemMatsPreconditionedBlock::MatSetValues(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora) {
  MoFEMFunctionBegin;
  CHKERR assembleStorage(row_data, col_data, mat, iora);
  CHKERR SchurBackendMatSetValuesPtr::matSetValuesPreconditionedBlockPtr(
      M, row_data, col_data, mat, iora);
  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode MatSetValues<SchurElemMatsPreconditionedBlock>(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora) {
  return SchurElemMatsPreconditionedBlock::MatSetValues(M, row_data, col_data,
                                                        mat, iora);
}

MoFEMErrorCode
schurSwitchPreconditioner(boost::shared_ptr<BlockStructure> block_mat_data) {
  MoFEMFunctionBegin;
  if (block_mat_data->multiplyByPreconditioner) {
    block_mat_data->multiplyByPreconditioner = false;
  } else {
    if (!block_mat_data->preconditionerBlocksPtr) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "preconditionerBlocksPtr not set");
    }
    block_mat_data->multiplyByPreconditioner = true;
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
schurSaveBlockMesh(boost::shared_ptr<BlockStructure> block_mat_data,
                   std::string filename) {
  MoFEMFunctionBegin;

  moab::Core core;
  moab::Interface &moab = core;

  ReadUtilIface *iface;
  CHKERR moab.query_interface(iface);

  auto size = block_mat_data->blockIndex.size();

  EntityHandle starting_vertex;
  vector<double *> arrays_coord;
  CHKERR iface->get_node_coords(3, 4 * size, 0, starting_vertex, arrays_coord);
  EntityHandle starting_handle;
  EntityHandle *array;
  CHKERR iface->get_element_connect(size, 4, MBQUAD, 0, starting_handle, array);

  double def_val_nrm2 = 0;
  Tag th_nrm2;
  CHKERR moab.tag_get_handle("nrm2", 1, MB_TYPE_DOUBLE, th_nrm2,
                             MB_TAG_CREAT | MB_TAG_DENSE, &def_val_nrm2);

  int def_val_mat_shift = 0;
  Tag th_mat_shift;
  CHKERR moab.tag_get_handle("mat_shift", 1, MB_TYPE_INTEGER, th_mat_shift,
                             MB_TAG_CREAT | MB_TAG_DENSE, &def_val_mat_shift);

  int i = 0;
  for (auto &d : block_mat_data->blockIndex) {
    auto row = d.getRow();
    auto col = d.getCol();
    auto nb_rows = d.getNbRows();
    auto nb_cols = d.getNbCols();
    auto mat_shift = d.getMatShift();

    // q0
    arrays_coord[0][4 * i + 0] = row;
    arrays_coord[1][4 * i + 0] = col;
    arrays_coord[2][4 * i + 0] = 0;

    // q1
    arrays_coord[0][4 * i + 1] = row + nb_rows;
    arrays_coord[1][4 * i + 1] = col;
    arrays_coord[2][4 * i + 1] = 0;

    // q2
    arrays_coord[0][4 * i + 2] = row + nb_rows;
    arrays_coord[1][4 * i + 2] = col + nb_cols;
    arrays_coord[2][4 * i + 2] = 0;

    // q3
    arrays_coord[0][4 * i + 3] = row;
    arrays_coord[1][4 * i + 3] = col + nb_cols;
    arrays_coord[2][4 * i + 3] = 0;

    // ele conn
    array[4 * i + 0] = starting_vertex + 4 * i + 0;
    array[4 * i + 1] = starting_vertex + 4 * i + 1;
    array[4 * i + 2] = starting_vertex + 4 * i + 2;
    array[4 * i + 3] = starting_vertex + 4 * i + 3;

    auto prt = &(*block_mat_data->dataBlocksPtr)[d.getMatShift()];
    auto nrm2 = cblas_dnrm2(nb_rows * nb_cols, prt, 1);
    EntityHandle ele = starting_handle + i;
    CHKERR moab.tag_set_data(th_nrm2, &ele, 1, &nrm2);
    CHKERR moab.tag_set_data(th_mat_shift, &ele, 1, &mat_shift);

    ++i;
  }

  CHKERR iface->update_adjacencies(starting_handle, size, 4, array);
  CHKERR moab.write_file(filename.c_str(), "VTK", "");

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM