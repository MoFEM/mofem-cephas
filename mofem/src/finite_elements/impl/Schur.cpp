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

/**
 * @brief Clear Schur complement internal data
 *
 */
struct OpSchurAssembleBegin : public OpSchurAssembleBase {

  OpSchurAssembleBegin();

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

/**
 * @brief Assemble Schur complement (Implementation)
 *
 */
struct OpSchurAssembleEndImpl : public OpSchurAssembleBase {

  /**
   * @brief Construct a new Op Schur Assemble End object
   *
   * @param fields_name list of fields
   * @param field_ents list of entities on which schur complement is applied
   * (can be empty)
   * @param sequence_of_aos list of maps from base problem to Schur complement
   * matrix
   * @param sequence_of_mats list of Schur complement matrices
   * @param sym_schur true if Schur complement is symmetric
   * @param symm_op true if block diagonal is symmetric
   */
  OpSchurAssembleEndImpl(
      std::vector<std::string> fields_name,
      std::vector<boost::shared_ptr<Range>> field_ents,
      std::vector<SmartPetscObj<AO>> sequence_of_aos,
      std::vector<SmartPetscObj<Mat>> sequence_of_mats,
      std::vector<bool> sym_schur, bool symm_op = true,
      boost::shared_ptr<BlockStruture> diag_blocks = nullptr);

  /**
   * @brief Construct a new Op Schur Assemble End object
   *
   * @param fields_name list of fields
   * @param field_ents list of entities on which schur complement is applied
   * (can be empty)
   * @param sequence_of_aos list of maps from base problem to Schur complement
   * matrix
   * @param sequence_of_mats list of Schur complement matrices
   * @param sym_schur true if Schur complement is symmetric
   * @param diag_eps add epsilon on diagonal of inverted matrix
   * @param symm_op true if block diagonal is symmetric
   */
  OpSchurAssembleEndImpl(
      std::vector<std::string> fields_name,
      std::vector<boost::shared_ptr<Range>> field_ents,
      std::vector<SmartPetscObj<AO>> sequence_of_aos,
      std::vector<SmartPetscObj<Mat>> sequence_of_mats,
      std::vector<bool> sym_schur, std::vector<double> diag_eps,
      bool symm_op = true,
      boost::shared_ptr<BlockStruture> diag_blocks = nullptr);

protected:
  template <typename I>
  MoFEMErrorCode doWorkImpl(int side, EntityType type,
                            EntitiesFieldData::EntData &data);

  std::vector<std::string> fieldsName;
  std::vector<boost::shared_ptr<Range>> fieldEnts;
  std::vector<SmartPetscObj<AO>> sequenceOfAOs;
  std::vector<SmartPetscObj<Mat>> sequenceOfMats;
  std::vector<bool> symSchur;
  std::vector<double> diagEps;

  MatrixDouble invMat;
  MatrixDouble invDiagOffMat;
  MatrixDouble offMatInvDiagOffMat;
  MatrixDouble transOffMatInvDiagOffMat;

  boost::shared_ptr<BlockStruture> diagBlocks;
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
struct OpSchurAssembleEnd<SchurDSYSV> : public OpSchurAssembleEndImpl {
  using OpSchurAssembleEndImpl::OpSchurAssembleEndImpl;
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

template <>
struct OpSchurAssembleEnd<SchurDGESV> : public OpSchurAssembleEndImpl {
  using OpSchurAssembleEndImpl::OpSchurAssembleEndImpl;
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

  const UId uidRow;
  const UId uidCol;

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

  const size_t iDX;

  friend OpSchurAssembleBegin;
  friend OpSchurAssembleEndImpl;

  struct idx_mi_tag {};
  struct uid_mi_tag {};
  struct row_mi_tag {};
  struct col_mi_tag {};

  using SchurElemStorage = multi_index_container<
      SchurElemMats,
      indexed_by<

          ordered_unique<
              tag<uid_mi_tag>,
              composite_key<
                  SchurElemMats,

                  member<SchurElemMats, const UId, &SchurElemMats::uidRow>,
                  member<SchurElemMats, const UId, &SchurElemMats::uidCol>

                  >>,

          ordered_unique<tag<idx_mi_tag>, member<SchurElemMats, const size_t,
                                                 &SchurElemMats::iDX>>,

          ordered_non_unique<tag<row_mi_tag>, member<SchurElemMats, const UId,
                                                     &SchurElemMats::uidRow>>,

          ordered_non_unique<tag<col_mi_tag>, member<SchurElemMats, const UId,
                                                     &SchurElemMats::uidCol>>

          >>;

  static boost::ptr_vector<MatrixDouble> locMats;
  static boost::ptr_vector<VectorInt> rowIndices;
  static boost::ptr_vector<VectorInt> colIndices;
  static SchurElemStorage schurL2Storage;
};

struct DiagBlockIndex {

  virtual ~DiagBlockIndex() = default;

  /**
   * @brief block data indexes
   *
   */
  struct Indexes {
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

          ordered_non_unique<

              composite_key<Indexes,

                            member<Indexes, int, &Indexes::col>,
                            member<Indexes, int, &Indexes::nb_cols>,
                            member<Indexes, int, &Indexes::row>,
                            member<Indexes, int, &Indexes::nb_rows>>

              >,

          ordered_unique<

              composite_key<Indexes,

                            member<Indexes, int, &Indexes::row>,
                            member<Indexes, int, &Indexes::col>,
                            member<Indexes, int, &Indexes::nb_rows>,
                            member<Indexes, int, &Indexes::nb_cols>

                            >>

          >>;

  BlockIndex blockIndex; ///< blocks indexes storage
};

struct BlockStruture : public DiagBlockIndex {

  SmartPetscObj<Vec> ghostX;
  SmartPetscObj<Vec> ghostY;

  boost::shared_ptr<std::vector<double>> dataBlocksPtr;
  boost::shared_ptr<std::vector<double>> dataInvBlocksPtr;
  boost::shared_ptr<std::vector<double>> preconditionerBlocksPtr;
  boost::shared_ptr<std::vector<double>> parentBlockStruturePtr;

  bool multiplyByPreconditioner = false;
};

struct DiagBlockInvStruture : public BlockStruture {

  using SchurSolverView =
      std::pair<std::vector<const Indexes *>, std::vector<int>>;

  SchurSolverView indexView;
};

PetscLogEvent SchurEvents::MOFEM_EVENT_schurMatSetValues;
PetscLogEvent SchurEvents::MOFEM_EVENT_opSchurAssembleEnd;
PetscLogEvent SchurEvents::MOFEM_EVENT_blockStrutureSetValues;
PetscLogEvent SchurEvents::MOFEM_EVENT_blockStrutureMult;
PetscLogEvent SchurEvents::MOFEM_EVENT_blockStrutureSolve;
PetscLogEvent SchurEvents::MOFEM_EVENT_zeroRowsAndCols;

SchurEvents::SchurEvents() {
  PetscLogEventRegister("schurMatSetVal", 0, &MOFEM_EVENT_schurMatSetValues);
  PetscLogEventRegister("opSchurAsmEnd", 0, &MOFEM_EVENT_opSchurAssembleEnd);
  PetscLogEventRegister("blockSetVal", 0, &MOFEM_EVENT_blockStrutureSetValues);
  PetscLogEventRegister("blockMult", 0, &MOFEM_EVENT_blockStrutureMult);
  PetscLogEventRegister("blockSolve", 0, &MOFEM_EVENT_blockStrutureSolve);
  PetscLogEventRegister("schurZeroRandC", 0, &MOFEM_EVENT_zeroRowsAndCols);
}

SchurElemMats::SchurElemStorage SchurElemMats::schurL2Storage;
boost::ptr_vector<MatrixDouble> SchurElemMats::locMats;
boost::ptr_vector<VectorInt> SchurElemMats::rowIndices;
boost::ptr_vector<VectorInt> SchurElemMats::colIndices;

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

  // get size of storage
  const auto idx = SchurElemMats::schurL2Storage.size();
  // get size of arrays of matrices
  const auto size = SchurElemMats::locMats.size();

  // expand memory allocation
  if (idx == size) {
    SchurElemMats::locMats.push_back(new MatrixDouble());
    SchurElemMats::rowIndices.push_back(new VectorInt());
    SchurElemMats::colIndices.push_back(new VectorInt());
  } else if (idx > size) {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Wrong size %d != %d",
             idx, size);
  }

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

  // add matrix to storage
  auto p = SchurElemMats::schurL2Storage.emplace(idx, uid_row, uid_col);
  auto get_storage = [&p]() { return const_cast<SchurElemMats &>(*p.first); };

  if (p.second) {
    // new entry is created

    auto asmb = [&](auto &sm) {
      sm.resize(nb_rows, nb_cols, false);
      noalias(sm) = mat;
    };

    asmb(get_storage().getMat());

    auto add_indices = [](auto &storage, auto &ind) {
      storage.resize(ind.size(), false);
      noalias(storage) = ind;
    };

    add_indices(get_storage().getRowInd(), row_ind);
    add_indices(get_storage().getColInd(), col_ind);

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

    CHKERR asmb(get_storage().getMat());
    // no need to set indices
  }

#ifndef NDEBUG
  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_schurMatSetValues, 0, 0, 0, 0);
#endif // NDEBUG

  MoFEMFunctionReturn(0);
}

OpSchurAssembleBegin::OpSchurAssembleBegin()
    : OpSchurAssembleBase(NOSPACE, OPSPACE) {}

MoFEMErrorCode OpSchurAssembleBegin::doWork(int side, EntityType type,
                                            EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
#ifndef NDEBUG
  MOFEM_LOG("SELF", Sev::noisy) << "Schur assemble begin";
#endif
  SchurElemMats::schurL2Storage.clear();

  MoFEMFunctionReturn(0);
}

OpSchurAssembleEndImpl::OpSchurAssembleEndImpl(
    std::vector<std::string> fields_name,
    std::vector<boost::shared_ptr<Range>> field_ents,
    std::vector<SmartPetscObj<AO>> sequence_of_aos,
    std::vector<SmartPetscObj<Mat>> sequence_of_mats,
    std::vector<bool> sym_schur, std::vector<double> diag_eps, bool symm_op,
    boost::shared_ptr<BlockStruture> diag_blocks)
    : OpSchurAssembleBase(NOSPACE, OPSPACE, symm_op), fieldsName(fields_name),
      fieldEnts(field_ents), sequenceOfAOs(sequence_of_aos),
      sequenceOfMats(sequence_of_mats), symSchur(sym_schur), diagEps(diag_eps),
      diagBlocks(diag_blocks) {}

OpSchurAssembleEndImpl::OpSchurAssembleEndImpl(
    std::vector<std::string> fields_name,
    std::vector<boost::shared_ptr<Range>> field_ents,
    std::vector<SmartPetscObj<AO>> sequence_of_aos,
    std::vector<SmartPetscObj<Mat>> sequence_of_mats,
    std::vector<bool> sym_schur, bool symm_op,
    boost::shared_ptr<BlockStruture> diag_blocks)
    : OpSchurAssembleEndImpl(
          fields_name, field_ents, sequence_of_aos, sequence_of_mats, sym_schur,
          std::vector<double>(fields_name.size(), 0), symm_op, diag_blocks) {}

template <typename I>
MoFEMErrorCode
OpSchurAssembleEndImpl::doWorkImpl(int side, EntityType type,
                                   EntitiesFieldData::EntData &data) {

  MoFEMFunctionBegin;

  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_opSchurAssembleEnd, 0, 0, 0, 0);

#ifndef NDEBUG
  MOFEM_LOG("SELF", Sev::noisy) << "Schur assemble begin -> end";
#endif

  auto get_field_name = [&](auto uid) {
    return getPtrFE()->mField.get_field_name(field_bit_from_bit_number(
        FieldEntity::getFieldBitNumberFromUniqueId(uid)));
  };

  // Assemble Schur complement
  auto assemble_mat = [&](SmartPetscObj<Mat> M, MatSetValuesRaw mat_set_values,
                          auto &storage) {
    MoFEMFunctionBegin;
    if (M) {
      for (auto &s : storage) {
        auto &m = s.getMat();
        if (m.size1()) {
          auto &row_ind = s.getRowInd();
          auto &col_ind = s.getColInd();

          if (auto ierr = mat_set_values(M, row_ind.size(), &*row_ind.begin(),
                                         col_ind.size(), &*col_ind.begin(),
                                         &*m.data().begin(), ADD_VALUES)) {
#ifndef NDEBUG
            auto field_ents = getPtrFE()->mField.get_field_ents();
            auto row_ent_it = field_ents->find(s.uidRow);
            auto col_ent_it = field_ents->find(s.uidCol);
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
      }
    }
    MoFEMFunctionReturn(0);
  };

  auto apply_schur = [&](auto &storage,

                         auto lo_uid, auto hi_uid,

                         double eps, bool symm) {
    MoFEMFunctionBegin;

    // add off diagonal matrix, i.e. schur complement
    auto add_off_mat = [&](auto row_uid, auto col_uid, auto &row_ind,
                           auto &col_ind, auto &offMatInvDiagOffMat) {
      MoFEMFunctionBegin;

      const auto idx = SchurElemMats::schurL2Storage.size();
      const auto size = SchurElemMats::locMats.size();

      if (idx == size) {
        SchurElemMats::locMats.push_back(new MatrixDouble());
        SchurElemMats::rowIndices.push_back(new VectorInt());
        SchurElemMats::colIndices.push_back(new VectorInt());
      } else if (idx > size) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong size %d != %d", idx, size);
      }

      auto it = storage.template get<SchurElemMats::uid_mi_tag>().find(
          boost::make_tuple(row_uid, col_uid));

      if (it == storage.template get<SchurElemMats::uid_mi_tag>().end()) {

        auto get_index = [&](auto storage) {
          size_t tmp_idx = 0;
          for (auto &s : storage.template get<SchurElemMats::idx_mi_tag>()) {
            if (tmp_idx != s.iDX) {
              return tmp_idx;
            }
            ++tmp_idx;
          }
          return storage.size();
        };

        auto p = SchurElemMats::schurL2Storage.emplace(
            get_index(SchurElemMats::schurL2Storage), row_uid, col_uid);
        auto &mat = p.first->getMat();
        auto &set_row_ind = p.first->getRowInd();
        auto &set_col_ind = p.first->getColInd();

        set_row_ind.resize(row_ind.size(), false);
        noalias(set_row_ind) = row_ind;
        set_col_ind.resize(col_ind.size(), false);
        noalias(set_col_ind) = col_ind;

        mat.resize(offMatInvDiagOffMat.size1(), offMatInvDiagOffMat.size2(),
                   false);
        noalias(mat) = offMatInvDiagOffMat;

#ifndef NDEBUG
        if constexpr (debug_schur) {
          MOFEM_LOG("SELF", Sev::noisy) << "insert: " << get_field_name(row_uid)
                                        << " " << get_field_name(col_uid) << " "
                                        << mat.size1() << " " << mat.size2();
        }
#endif // NDEBUG

      } else {

        auto &mat = it->getMat();
#ifndef NDEBUG
        if constexpr (debug_schur) {
          MOFEM_LOG("SELF", Sev::noisy) << "add: " << get_field_name(row_uid)
                                        << " " << get_field_name(col_uid) << " "
                                        << mat.size1() << " " << mat.size2();
        }
#endif // NDEBUG

#ifndef NDEBUG
        if (mat.size1() != offMatInvDiagOffMat.size1()) {
          SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "Wrong size %d != %d", mat.size1(),
                   offMatInvDiagOffMat.size1());
        }
        if (mat.size2() != offMatInvDiagOffMat.size2()) {
          SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "Wrong size %d != %d", mat.size2(),
                   offMatInvDiagOffMat.size2());
        }
#endif // NDEBUG
        mat += offMatInvDiagOffMat;
      }
      MoFEMFunctionReturn(0);
    };

    auto get_row_view = [&]() {
      auto row_it =
          storage.template get<SchurElemMats::row_mi_tag>().lower_bound(lo_uid);
      auto hi_row_it =
          storage.template get<SchurElemMats::row_mi_tag>().upper_bound(hi_uid);
      std::vector<const SchurElemMats *> schur_row_ptr_view;
      schur_row_ptr_view.reserve(std::distance(row_it, hi_row_it));
      for (; row_it != hi_row_it; ++row_it) {
        schur_row_ptr_view.push_back(&*row_it);
      }
      return schur_row_ptr_view;
    };

    auto get_col_view = [&]() {
      auto col_it =
          storage.template get<SchurElemMats::col_mi_tag>().lower_bound(lo_uid);
      auto hi_col_it =
          storage.template get<SchurElemMats::col_mi_tag>().upper_bound(hi_uid);
      std::vector<const SchurElemMats *> schur_col_ptr_view;
      schur_col_ptr_view.reserve(std::distance(col_it, hi_col_it));
      for (; col_it != hi_col_it; ++col_it) {
        schur_col_ptr_view.push_back(&*col_it);
      }
      return schur_col_ptr_view;
    };

    auto schur_row_ptr_view = get_row_view();
    auto schur_col_ptr_view = get_col_view();

    // iterate row entities
    for (auto row_it : schur_row_ptr_view) {
      // only diagonals to get inverted diagonal
      if (row_it->uidRow == row_it->uidCol) {
#ifndef NDEBUG
        if constexpr (debug_schur) {
          MOFEM_LOG("SELF", Sev::noisy)
              << "invert: row_uid " << get_field_name(row_it->uidRow)
              << " row uid " << get_field_name(row_it->uidCol) << " : "
              << row_it->getMat().size1() << " " << row_it->getMat().size2();
        }
#endif // NDEBUG

        // note invMat is multiplied by -1
        CHKERR I::invertMat(row_it, invMat, eps);

        // iterate column entities
        for (auto c_lo : schur_col_ptr_view) {
#ifndef NDEBUG
          if (c_lo->uidCol != row_it->uidRow)
            SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                     "Wrong size %d != %d", c_lo->uidCol, row_it->uidRow);
#endif // NDEBUG

          auto &row_uid = c_lo->uidRow;
          if (row_uid == row_it->uidRow) {
            continue;
          }

          auto &cc_off_mat = c_lo->getMat();
          invDiagOffMat.resize(cc_off_mat.size1(), invMat.size2(), false);
#ifndef NDEBUG
          if (invMat.size1() != cc_off_mat.size2()) {
            MOFEM_LOG("SELF", Sev::error)
                << "row_uid " << get_field_name(row_it->uidRow) << " row uid "
                << get_field_name(row_uid);
            SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                     "Wrong size %d != %d", invMat.size1(), cc_off_mat.size2());
          }
#endif // NDEBUG

          // noalias(invDiagOffMat) = prod(cc_off_mat, invMat);
          // A10*inv(A00)
          cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                      cc_off_mat.size1(), invMat.size2(), cc_off_mat.size2(),
                      1., &*cc_off_mat.data().begin(), cc_off_mat.size2(),
                      &*invMat.data().begin(), invMat.size2(), 0.,
                      &*invDiagOffMat.data().begin(), invDiagOffMat.size2());

          // iterate row entities
          for (auto r_lo : schur_row_ptr_view) {
#ifndef NDEBUG
            if (c_lo->uidCol != row_it->uidRow)
              SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                       "Wrong size %d != %d", c_lo->uidCol, row_it->uidRow);
#endif // NDEBUG

            auto &col_uid = r_lo->uidCol;

            // Skip diagonal
            if (col_uid == row_it->uidRow) {
              continue;
            }

            // If symmetry only run upper off diagonal terms
            if (symm && row_uid > col_uid) {
              continue;
            }

            auto &rr_off_mat = r_lo->getMat();
            offMatInvDiagOffMat.resize(invDiagOffMat.size1(),
                                       rr_off_mat.size2(), false);
#ifndef NDEBUG
            if (rr_off_mat.size1() != invDiagOffMat.size2()) {
              MOFEM_LOG("SELF", Sev::error)
                  << "row_uid " << get_field_name(row_it->uidRow)
                  << ": col uid " << get_field_name(col_uid) << " row uid "
                  << get_field_name(row_uid);
              SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                       "Wrong size %d != %d", rr_off_mat.size1(),
                       invDiagOffMat.size2());
            }
#endif // NDEBUG

            // noalias(offMatInvDiagOffMat) = prod(invDiagOffMat, rr_off_mat);
            // A10*inv(A00)*A01
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        invDiagOffMat.size1(), rr_off_mat.size2(),
                        invDiagOffMat.size2(), 1.,
                        &*invDiagOffMat.data().begin(), invDiagOffMat.size2(),
                        &*rr_off_mat.data().begin(), rr_off_mat.size2(), 0.,
                        &*offMatInvDiagOffMat.data().begin(),
                        offMatInvDiagOffMat.size2());

            // Make on diagonal A11 have Schur complement S
            CHKERR add_off_mat(row_uid, col_uid, c_lo->getRowInd(),
                               r_lo->getColInd(), offMatInvDiagOffMat);

            if (symm && row_uid != col_uid) {
              transOffMatInvDiagOffMat.resize(offMatInvDiagOffMat.size2(),
                                              offMatInvDiagOffMat.size1(),
                                              false);
              noalias(transOffMatInvDiagOffMat) = trans(offMatInvDiagOffMat);
              CHKERR add_off_mat(col_uid, row_uid, r_lo->getColInd(),
                                 c_lo->getRowInd(), transOffMatInvDiagOffMat);
            }
          }
        }
      }
    }

    MoFEMFunctionReturn(0);
  };

  auto erase_factored = [&](auto &storage, auto lo_uid, auto hi_uid) {
    MoFEMFunctionBegin;

    auto r_lo =
        storage.template get<SchurElemMats::row_mi_tag>().lower_bound(lo_uid);
    auto r_hi =
        storage.template get<SchurElemMats::row_mi_tag>().upper_bound(hi_uid);
    storage.template get<SchurElemMats::row_mi_tag>().erase(r_lo, r_hi);

    auto c_lo =
        storage.template get<SchurElemMats::col_mi_tag>().lower_bound(lo_uid);
    auto c_hi =
        storage.template get<SchurElemMats::col_mi_tag>().upper_bound(hi_uid);
    storage.template get<SchurElemMats::col_mi_tag>().erase(c_lo, c_hi);

    MoFEMFunctionReturn(0);
  };

  auto assemble_S = [&](auto &storage, auto ao, auto mat) {
    MoFEMFunctionBegin;

    // apply AO
    if (ao) {
      for (auto &m : storage) {
        auto &ind_row = m.getRowInd();
        CHKERR AOApplicationToPetsc(ao, ind_row.size(), &*ind_row.begin());
        auto &ind_col = m.getColInd();
        CHKERR AOApplicationToPetsc(ao, ind_col.size(), &*ind_col.begin());
      }
    }

    // assemble matrix
    if (mat) {
      CHKERR assemble_mat(mat, matSetValuesSchurRaw, storage);
    }

    MoFEMFunctionReturn(0);
  };

  auto assemble_A00 = [&](auto &storage, auto lo_uid, auto hi_uid) {
    MoFEMFunctionBegin;

    auto add = [&](auto lo, auto hi) {
      MoFEMFunctionBegin;
      for (; lo != hi; ++lo) {
        auto &m = lo->getMat();
        if (m.size1() && m.size2()) {
          auto row_ind = lo->getRowInd()[0];
          auto col_ind = lo->getColInd()[0];
          auto nb_rows = m.size1();
          auto nb_cols = m.size2();
          auto it = diagBlocks->blockIndex.get<1>().find(
              boost::make_tuple(row_ind, col_ind, nb_rows, nb_cols));
          if (it != diagBlocks->blockIndex.get<1>().end()) {
            auto inv_shift = it->inv_shift;
            if (inv_shift != -1) {
#ifndef NDEBUG
              if (!diagBlocks->dataInvBlocksPtr)
                SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                        "No dataInvBlocksPtr");
#endif // NDEBUG
              auto *ptr = &((*diagBlocks->dataInvBlocksPtr)[inv_shift]);
              if (row_ind == col_ind && nb_rows == nb_cols) {
                // assemble inverted diag
                std::copy(invMat.data().begin(), invMat.data().end(), ptr);
              } else {
                // assemble of diag terms, witch might be changed by Schur
                // complement
                std::copy(m.data().begin(), m.data().end(), ptr);
              }
            }
          }
        }
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR add(

        storage.template get<SchurElemMats::row_mi_tag>().lower_bound(lo_uid),
        storage.template get<SchurElemMats::row_mi_tag>().upper_bound(hi_uid)

    );

    MoFEMFunctionReturn(0);
  };

  auto get_a00_uids = [&]() {
    auto get_field_bit = [&](auto &name) {
      return getPtrFE()->mField.get_field_bit_number(name);
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
          << "List schur storage: " << i << " " << p.iDX << ": "
          << get_field_name(p.uidRow) << " " << get_field_name(p.uidCol)
          << " : " << p.getMat().size1() << " " << p.getMat().size2();
      ++i;
    }
    MoFEMFunctionReturn(0);
  };
#endif // NDEBUG

  auto assemble = [&](auto &&a00_uids) {
    MoFEMFunctionBegin;
    auto &storage = SchurElemMats::schurL2Storage;
    int ss = 0;
    for (auto &p : a00_uids) {
      auto [lo_uid, hi_uid] = p;
#ifndef NDEBUG
      if constexpr (debug_schur) {
        list_storage(storage);
        MOFEM_LOG("SELF", Sev::noisy)
            << "Schur assemble: " << get_field_name(lo_uid) << " - "
            << get_field_name(hi_uid);
      }
#endif
      CHKERR apply_schur(storage, lo_uid, hi_uid, diagEps[ss], symSchur[ss]);
      if (diagBlocks)
        CHKERR assemble_A00(storage, lo_uid, hi_uid);
      CHKERR erase_factored(storage, lo_uid, hi_uid);
      CHKERR assemble_S(storage, sequenceOfAOs[ss], sequenceOfMats[ss]);
      ++ss;
    }
    MoFEMFunctionReturn(0);
  };

  // Assemble Schur complements
  CHKERR assemble(get_a00_uids());

#ifndef NDEBUG
  MOFEM_LOG("SELF", Sev::noisy) << "Schur assemble done";
#endif

  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_opSchurAssembleEnd, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

struct SchurDSYSV {
  static auto invertMat(const SchurElemMats *row_ptr, MatrixDouble &inv,
                        double eps) {
    MoFEMFunctionBeginHot;

    auto &m = row_ptr->getMat();

    VectorInt ipiv;
    VectorDouble lapack_work;
    const int nb = m.size1();

    if (eps) {
      for (int cc = 0; cc != nb; ++cc) {
        m(cc, cc) += eps;
      }
    }

    inv.resize(nb, nb, false);
    inv.clear();
    auto ptr = &*inv.data().begin();
    for (int c = 0; c != nb; ++c, ptr += nb + 1)
      *ptr = -1;
    ipiv.resize(nb, false);
    lapack_work.resize(nb * nb, false);
    const auto info =
        lapack_dsysv('L', nb, nb, &*m.data().begin(), nb, &*ipiv.begin(),
                     &*inv.data().begin(), nb, &*lapack_work.begin(), nb * nb);
    if (info != 0)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
               "Can not invert matrix info = %d", info);
    MoFEMFunctionReturnHot(0);
  };
};

struct SchurDGESV {
  static auto invertMat(const SchurElemMats *row_ptr, MatrixDouble &inv,
                        double eps) {
    MoFEMFunctionBeginHot;

    auto &m = row_ptr->getMat();

    VectorInt ipiv;
    const auto nb = m.size1();
#ifndef NDEBUG
    if (nb != m.size2()) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "It should be square matrix %d != %d", nb, m.size2());
    }
#endif

    if (eps) {
      for (int cc = 0; cc != nb; ++cc) {
        m(cc, cc) += eps;
      }
    }

    inv.resize(nb, nb, false);
    inv.clear();
    auto ptr = &*inv.data().begin();
    for (int c = 0; c != nb; ++c, ptr += nb + 1)
      *ptr = -1;
    ipiv.resize(nb, false);
    const auto info = lapack_dgesv(nb, nb, &*m.data().begin(), nb,
                                   &*ipiv.begin(), &*inv.data().begin(), nb);
    if (info != 0)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
               "Can not invert matrix info = %d", info);
    MoFEMFunctionReturnHot(0);
  };
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

boost::shared_ptr<BlockStruture> createBlockMatStructure(

    DM dm,                                //< dm
    SchurFEOpsFEandFields schur_fe_op_vec //< block elements

) {
  auto data_ptr = boost::make_shared<BlockStruture>();

  auto m_field_ptr = getInterfacePtr(dm);

  for (auto &d : schur_fe_op_vec) {
    auto fe_method = boost::shared_ptr<MoFEM::FEMethod>(new MoFEM::FEMethod());

    auto get_bit_numbers = [&d](auto op) {
      std::vector<FieldBitNumber> bit_numbers(d.second.size());
      std::transform(d.second.begin(), d.second.end(), bit_numbers.begin(), op);
      return bit_numbers;
    };

    auto get_uid_pair = [](const auto &field_id) {
      auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
          field_id, get_id_for_min_type<MBVERTEX>());
      auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
          field_id, get_id_for_max_type<MBENTITYSET>());
      return std::make_pair(lo_uid, hi_uid);
    };

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

    auto get_it_pair = [cmp_uid_lo, cmp_uid_hi](auto &&field_ents,
                                                auto &&p_uid) {
      auto lo = std::lower_bound(field_ents.begin(), field_ents.end(),
                                 p_uid.first, cmp_uid_lo);
      auto hi = std::upper_bound(field_ents.begin(), field_ents.end(),
                                 p_uid.second, cmp_uid_hi);
      return std::make_pair(lo, hi);
    };

    auto row_extractor = [](auto &e) { return e->entityCacheRowDofs; };
    auto col_extractor = [](auto &e) { return e->entityCacheColDofs; };

    auto extract_data = [](auto &&its, auto extractor) {
      std::vector<std::tuple<int, int, int>> data;
      data.reserve(std::distance(its.first, its.second));
      for (; its.first != its.second; ++its.first) {
        if (auto e = its.first->lock()) {
          if (auto cache = extractor(e).lock()) {
            auto nb_dofs = std::distance(cache->loHi[0], cache->loHi[1]);
            if (nb_dofs) {
              auto glob = (*cache->loHi[0])->getPetscGlobalDofIdx();
              auto loc = (*cache->loHi[0])->getPetscLocalDofIdx();
              data.emplace_back(glob, nb_dofs, loc);
            }
          }
        }
      }
      return data;
    };

    auto row_bit_numbers = get_bit_numbers(
        [&](auto &p) { return m_field_ptr->get_field_bit_number(p.first); });
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
          auto [r_glob, r_nb_dofs, r_loc] = r;
          for (auto &c : col_data) {
            auto [c_glob, c_nb_dofs, c_loc] = c;
            if (r_glob != -1 && c_glob != -1) {
              data_ptr->blockIndex.insert(BlockStruture::Indexes{
                  r_glob, c_glob, r_nb_dofs, c_nb_dofs, r_loc, c_loc, -1, -1});
            }
          }
        }
      }

      MoFEMFunctionReturn(0);
    };

    CHKERR DMoFEMLoopFiniteElements(dm, d.first, fe_method);
  };

  // order by column (that is for matrix multiplication)
  auto mem_size = 0;
  for (auto &v : data_ptr->blockIndex.get<0>()) {
    v.mat_shift = mem_size;
    mem_size += v.nb_cols * v.nb_rows;
  }

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

static MoFEMErrorCode mult_schur_block_shell(Mat mat, Vec x, Vec y,
                                             InsertMode iora);

static MoFEMErrorCode solve_schur_block_shell(Mat mat, Vec x, Vec y,
                                              InsertMode iora);

static PetscErrorCode mult(Mat mat, Vec x, Vec y) {
  return mult_schur_block_shell(mat, x, y, INSERT_VALUES);
}
static PetscErrorCode mult_add(Mat mat, Vec x, Vec y) {
  return mult_schur_block_shell(mat, x, y, ADD_VALUES);
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
  BlockStruture *ctx;
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

  struct ShiftedBlockView : BlockStruture::Indexes {
    ShiftedBlockView() = delete;
    inline auto rowShift() const {
      return row + nb_rows;
    } // shift such that lower bound is included
    inline auto colShift() const {
      return col + nb_cols;
    } // shift such that lower bound is included
  };

  // this enable esrch by varying ranges
  using BlockIndexView = multi_index_container<

      const ShiftedBlockView *,

      indexed_by<

          ordered_non_unique<

              const_mem_fun<ShiftedBlockView, int, &ShiftedBlockView::rowShift>

              >,

          ordered_non_unique<
              const_mem_fun<ShiftedBlockView, int, &ShiftedBlockView::colShift>>

          >>;

  int max_dofs_on_ent = 0;
  BlockIndexView view;
  for (auto &v : ctx->blockIndex.get<0>()) {
    view.insert(static_cast<const ShiftedBlockView *>(&v));
    max_dofs_on_ent = std::max(max_dofs_on_ent, v.nb_rows);
  }

  for (auto n = 0; n != is_nb_rows; ++n) {
    auto row = ptr[n];
    auto rlo = view.get<0>().lower_bound(row);
    auto rhi = view.get<0>().upper_bound(
        row + max_dofs_on_ent); // add such that upper bound is included
    for (; rlo != rhi; ++rlo) {
      auto r_shift = row - (*rlo)->row;
      if (r_shift >= 0 && r_shift < (*rlo)->nb_rows) {
        auto *ptr = &(*ctx->dataBlocksPtr)[(*rlo)->mat_shift];
        for (auto i = 0; i != (*rlo)->nb_cols; ++i) {
          ptr[i + r_shift * (*rlo)->nb_cols] = 0;
        }
      } else if ((*rlo)->row > row) {
        break;
      }
    }
  }

  for (auto n = 0; n != is_nb_rows; ++n) {
    auto col = ptr[n];
    auto clo = view.get<1>().lower_bound(col);
    auto chi = view.get<1>().upper_bound(
        col + max_dofs_on_ent); // add such that upper bound is included
    for (; clo != chi; ++clo) {
      auto c_shift = col - (*clo)->col;
      if (c_shift >= 0 && c_shift < (*clo)->nb_cols) {

        auto *ptr = &(*ctx->dataBlocksPtr)[(*clo)->mat_shift];
        for (auto i = 0; i != (*clo)->nb_rows; ++i) {
          ptr[c_shift + i * (*clo)->nb_cols] = 0;
        }

        // diagonal
        if ((*clo)->row == (*clo)->col && (*clo)->loc_row < loc_m &&
            (*clo)->loc_col < loc_n) {
          auto r_shift = col - (*clo)->col;
          if (r_shift >= 0 && r_shift < (*clo)->nb_rows) {
            ptr[c_shift + r_shift * (*clo)->nb_cols] = diag;
          }
        }
      } else if ((*clo)->row > col) {
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
  BlockStruture *ctx;
  CHKERR MatShellGetContext(m, (void **)&ctx);
  ctx->dataBlocksPtr->clear();
  if (ctx->dataInvBlocksPtr)
    ctx->dataInvBlocksPtr->clear();
  if (ctx->preconditionerBlocksPtr)
    ctx->preconditionerBlocksPtr->clear();
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

SchurShellMatData createBlockMat(DM dm, boost::shared_ptr<BlockStruture> data) {

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

static MoFEMErrorCode mult_schur_block_shell(Mat mat, Vec x, Vec y,
                                             InsertMode iora) {
  MoFEMFunctionBegin;
  BlockStruture *ctx;
  CHKERR MatShellGetContext(mat, (void **)&ctx);

  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_blockStrutureMult, 0, 0, 0, 0);

  Vec ghost_x = ctx->ghostX;
  Vec ghost_y = ctx->ghostY;

  CHKERR VecCopy(x, ghost_x);

  CHKERR VecGhostUpdateBegin(ghost_x, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(ghost_x, INSERT_VALUES, SCATTER_FORWARD);

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
  for (auto i = 0; i != nb_y; ++i)
    y_array[i] = 0.;

  double *block_ptr = &*ctx->dataBlocksPtr->begin();

  auto it = ctx->blockIndex.get<0>().lower_bound(0);
  auto hi = ctx->blockIndex.get<0>().end();

  while (it != hi) {
    auto nb_rows = it->nb_rows;
    auto nb_cols = it->nb_cols;
    auto x_ptr = &x_array[it->loc_col];
    auto y_ptr = &y_array[it->loc_row];
    auto ptr = &block_ptr[it->mat_shift];
    for (auto r = 0; r != nb_rows; ++r) {
      for (auto c = 0; c != nb_cols; ++c) {
        y_ptr[r] += ptr[r * nb_cols + c] * x_ptr[c];
      }
    }
    ++it;
  }

  if (ctx->multiplyByPreconditioner) {

    if (!ctx->preconditionerBlocksPtr)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "No parentBlockStruturePtr");

    auto preconditioner_ptr = &*ctx->preconditionerBlocksPtr->begin();

    auto it = ctx->blockIndex.get<0>().lower_bound(0);
    auto hi = ctx->blockIndex.get<0>().end();

    while (it != hi) {
      if (it->row == it->col && it->inv_shift != -1) {
        auto nb_rows = it->nb_rows;
        auto nb_cols = it->nb_cols;
        auto x_ptr = &x_array[it->loc_col];
        auto y_ptr = &y_array[it->loc_row];
        auto ptr = &preconditioner_ptr[it->inv_shift];
        for (auto r = 0; r != nb_rows; ++r) {
          for (auto c = 0; c != nb_cols; ++c) {
            y_ptr[r] += ptr[r * nb_cols + c] * x_ptr[c];
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
  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_blockStrutureMult, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

static MoFEMErrorCode solve_schur_block_shell(Mat mat, Vec y, Vec x,
                                              InsertMode iora) {
  MoFEMFunctionBegin;
  BlockStruture *ctx;
  CHKERR MatShellGetContext(mat, (void **)&ctx);

  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_blockStrutureSolve, 0, 0, 0, 0);

  // Note that for solver those two are swapped
  Vec ghost_x = ctx->ghostY;
  Vec ghost_y = ctx->ghostX;

  CHKERR VecCopy(y, ghost_y);
  CHKERR VecZeroEntries(ghost_x);

  CHKERR VecGhostUpdateBegin(ghost_y, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(ghost_y, INSERT_VALUES, SCATTER_FORWARD);

  double *x_array;
  Vec loc_ghost_x;
  CHKERR VecGhostGetLocalForm(ghost_x, &loc_ghost_x);
  CHKERR VecGetArray(loc_ghost_x, &x_array);

  double *y_array;
  Vec loc_ghost_y;
  CHKERR VecGhostGetLocalForm(ghost_y, &loc_ghost_y);
  CHKERR VecGetArray(loc_ghost_y, &y_array);

  auto data_inv_blocks = ctx->dataInvBlocksPtr;
  auto inv_block_ptr = &*data_inv_blocks->begin();
  auto data_blocks = ctx->dataBlocksPtr;
  auto block_ptr = &*data_blocks->begin();

  auto *data_inv = dynamic_cast<DiagBlockInvStruture *>(ctx);
  auto index_view = &data_inv->indexView;

  std::vector<double> f;

  for (auto s1 = 0; s1 != index_view->second.size() - 1; ++s1) {
    auto lo = index_view->second[s1];
    auto hi = index_view->second[s1 + 1];

    auto diag_index_ptr = index_view->first[lo];
    ++lo;

    auto row = diag_index_ptr->loc_row;
    auto col = diag_index_ptr->loc_col;
    auto nb_rows = diag_index_ptr->nb_rows;
    auto nb_cols = diag_index_ptr->nb_cols;
    auto inv_shift = diag_index_ptr->inv_shift;

    f.resize(nb_cols);
    std::copy(&y_array[col], &y_array[col + nb_cols], f.begin());

    for (; lo != hi; ++lo) {
      auto off_index_ptr = index_view->first[lo];
      auto off_col = off_index_ptr->loc_col;
      auto off_nb_cols = off_index_ptr->nb_cols;
      auto off_shift = off_index_ptr->mat_shift;
      auto x_ptr = &x_array[off_col];
      auto ptr = &block_ptr[off_shift];
      for (auto r = 0; r != nb_rows; ++r) {
        for (auto c = 0; c != off_nb_cols; ++c) {
          f[r] -= ptr[r * off_nb_cols + c] * x_ptr[c];
        }
      }
    }

#ifndef NDEBUG
    if (inv_shift == -1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "inv_shift == -1");
    if (inv_shift + nb_rows * nb_cols > data_inv_blocks->size())
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "inv_shift out of range %d > %d", inv_shift + nb_rows * nb_cols,
               data_inv_blocks->size());
#endif // NDEBUG

    auto ptr = &inv_block_ptr[inv_shift];
    for (auto r = 0; r != nb_rows; ++r) {
      for (auto c = 0; c != nb_cols; ++c) {
        x_array[row + r] += inv_block_ptr[inv_shift + r * nb_cols + c] * f[c];
      }
    }
  }

  CHKERR VecRestoreArray(loc_ghost_x, &x_array);
  CHKERR VecRestoreArray(loc_ghost_y, &y_array);
  CHKERR VecGhostRestoreLocalForm(ghost_x, &loc_ghost_x);
  CHKERR VecGhostRestoreLocalForm(ghost_y, &loc_ghost_y);

  CHKERR VecGhostUpdateBegin(ghost_y, ADD_VALUES, SCATTER_REVERSE);
  CHKERR VecGhostUpdateEnd(ghost_y, ADD_VALUES, SCATTER_REVERSE);

  switch (iora) {
  case INSERT_VALUES:
    CHKERR VecCopy(ghost_x, x);
    break;
  case ADD_VALUES:
    CHKERR VecAXPY(x, 1., ghost_x);
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
    print_norm("solve_schur_block_shell insert x", x);
    print_norm("solve_schur_block_shell insert y", y);
    break;
  case ADD_VALUES:
    print_norm("solve_schur_block_shell add x", x);
    print_norm("solve_schur_block_shell add y", y);
    break;
  default:
    CHK_MOAB_THROW(MOFEM_NOT_IMPLEMENTED, "Wrong InsertMode");
  }

#endif // NDEBUG

  // PetscLogFlops(xxx)
  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_blockStrutureSolve, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
shell_block_mat_asmb_wrap_impl(BlockStruture *ctx,
                               const EntitiesFieldData::EntData &row_data,
                               const EntitiesFieldData::EntData &col_data,
                               const MatrixDouble &mat, InsertMode iora) {

  MatrixDouble tmp_mat;
  MoFEMFunctionBegin;

#ifndef NDEBUG

  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_blockStrutureSetValues, 0, 0, 0,
                     0);

#endif // NDEBUG

  auto get_row_indices_ptr = [&row_data]() -> const VectorInt * {
    boost::shared_ptr<EssentialBcStorage> stored_data_ptr;
    if (!row_data.getFieldEntities().empty()) {
      if (auto e_ptr = row_data.getFieldEntities()[0]) {
        if (auto stored_data_ptr =
                e_ptr->getSharedStoragePtr<EssentialBcStorage>()) {
          return &(stored_data_ptr->entityIndices);
        }
      }
    }
    return &row_data.getIndices();
  };

  auto get_negative_indices = [](auto &&row_indices_ptr) {
    std::vector<int> negative;
    int i = 0;
    for (auto &v : *row_indices_ptr) {
      if (v == -1)
        negative.push_back(i);
      ++i;
    }
    return negative;
  };

  // that is not tested (fix it)
  auto skip_negative_indices_mat =
      [&mat, &tmp_mat](auto &&negative) -> const MatrixDouble * {
    if (negative.size()) {
      tmp_mat.resize(mat.size1(), mat.size2());
      noalias(tmp_mat) = mat;
      for (auto v : negative) {
        ublas::matrix_row<decltype(tmp_mat)> mr(tmp_mat, v);
        std::fill(mr.begin(), mr.end(), 0);
      }
      return &tmp_mat;
    } else {
      return &mat;
    }
  };

  auto get_row_first_index = [](auto *row_index_ptr, int nb_r) {
    for (auto i = 0; i != nb_r; ++i) {
      if (row_index_ptr[i] != -1)
        return row_index_ptr[i] - i;
      return -1;
    }
    return -1;
  };

  auto m_ptr =
      skip_negative_indices_mat(get_negative_indices(get_row_indices_ptr()));

  auto nb_rows = row_data.getIndices().size();
  auto nb_cols = col_data.getIndices().size();

  auto set = [&](auto r, auto c, auto nb_r, auto nb_c, auto &mat) {
    MoFEMFunctionBegin;

    auto row_first_index =
        get_row_first_index(&(*get_row_indices_ptr())[r], nb_r);
    if (row_first_index != -1) {

      auto size = nb_r * nb_c;
      auto it = ctx->blockIndex.get<1>().find(boost::make_tuple(
          row_first_index, col_data.getIndices()[c], nb_r, nb_c));

#ifndef NDEBUG

      if (it == ctx->blockIndex.get<1>().end()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Block not allocated");
      }
      if (it->nb_rows != nb_r) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong size %d != %d", it->nb_rows, nb_r);
      }
      if (it->nb_cols != nb_c) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong size %d != %d", it->nb_cols, nb_c);
      }
      if (nb_r != mat.size1()) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong size %d != %d", nb_r, mat.size1());
      }
      if (nb_c != mat.size2()) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong size %d != %d", nb_c, mat.size2());
      }
      if (nb_r * nb_c != mat.data().size()) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong size %d != %d", nb_r * nb_c, mat.data().size());
      }

#endif // NDEBUG

      auto shift = it->mat_shift;
      if (shift != -1) {

        if (iora == ADD_VALUES) {
          cblas_daxpy(size, 1., &*mat.data().begin(), 1,
                      &(*ctx->dataBlocksPtr)[shift], 1);
        } else {
          std::copy(mat.data().begin(), mat.data().end(),
                    &(*ctx->dataBlocksPtr)[shift]);
        }
      } else {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "shift == -1");
      }
    }
    MoFEMFunctionReturn(0);
  };

  if (nb_rows && nb_cols) {
    auto row_type = row_data.getFieldDofs()[0]->getEntType();
    auto col_type = col_data.getFieldDofs()[0]->getEntType();

    if (row_type == MBVERTEX && col_type == MBVERTEX) {
      auto row_nb_coeff = row_data.getFieldDofs()[0]->getNbOfCoeffs();
      auto col_nb_coeff = col_data.getFieldDofs()[0]->getNbOfCoeffs();
      MatrixDouble m(row_nb_coeff, col_nb_coeff);
      for (auto i = 0; i != nb_rows / row_nb_coeff; ++i) {
        for (auto j = 0; j != nb_cols / col_nb_coeff; ++j) {

          for (auto k = 0; k != row_nb_coeff; ++k) {
            for (auto l = 0; l != col_nb_coeff; ++l) {
              m(k, l) = (*m_ptr)(i * row_nb_coeff + k, j * col_nb_coeff + l);
            }
          }

          CHKERR set(i * row_nb_coeff, j * col_nb_coeff, row_nb_coeff,
                     col_nb_coeff, m);
        }
      }

    } else if (row_type == MBVERTEX) {
      auto row_nb_coeff = row_data.getFieldDofs()[0]->getNbOfCoeffs();
      MatrixDouble m(row_nb_coeff, nb_cols);
      for (auto i = 0; i != nb_rows / row_nb_coeff; ++i) {

        for (auto k = 0; k != row_nb_coeff; ++k) {
          for (auto l = 0; l != nb_cols; ++l) {
            m(k, l) = (*m_ptr)(i * row_nb_coeff + k, l);
          }
        }

        CHKERR set(i * row_nb_coeff, 0, row_nb_coeff, nb_cols, m);
      }

    } else if (col_type == MBVERTEX) {
      auto col_nb_coeff = col_data.getFieldDofs()[0]->getNbOfCoeffs();
      MatrixDouble m(nb_rows, col_nb_coeff);
      for (auto j = 0; j != nb_cols / col_nb_coeff; ++j) {
        for (auto k = 0; k != nb_rows; ++k) {
          for (auto l = 0; l != col_nb_coeff; ++l) {
            m(k, l) = (*m_ptr)(k, j * col_nb_coeff + l);
          }
        }
        CHKERR set(0, j * col_nb_coeff, nb_rows, col_nb_coeff, m);
      }

    } else {

      CHKERR set(0, 0, nb_rows, nb_cols, (*m_ptr));
    }
  }

#ifndef NDEBUG
  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_blockStrutureSetValues, 0, 0, 0, 0);
#endif // NDEBUG

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
shell_block_mat_asmb_wrap(Mat M, const EntitiesFieldData::EntData &row_data,
                          const EntitiesFieldData::EntData &col_data,
                          const MatrixDouble &mat, InsertMode iora) {
  MoFEMFunctionBegin;
  BlockStruture *ctx;
  CHKERR MatShellGetContext(M, (void **)&ctx);
  CHKERR shell_block_mat_asmb_wrap_impl(ctx, row_data, col_data, mat, iora);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode shell_block_preconditioner_mat_asmb_wrap_impl(
    BlockStruture *ctx, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora) {

  MoFEMFunctionBegin;

#ifndef NDEBUG
  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_blockStrutureSetValues, 0, 0, 0,
                     0);
#endif // NDEBUG

  if (!ctx->preconditionerBlocksPtr) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "preconditionerBlocksPtr not set");
  }

  auto nb_rows = row_data.getIndices().size();
  auto nb_cols = col_data.getIndices().size();

  auto set = [&](auto r, auto c, auto nb_r, auto nb_c, auto &mat) {
    MoFEMFunctionBegin;

    auto row_first_index = row_data.getIndices()[r];
    if (row_first_index != -1) {

      auto size = nb_r * nb_c;
      auto it = ctx->blockIndex.get<1>().find(boost::make_tuple(
          row_first_index, col_data.getIndices()[c], nb_r, nb_c));

#ifndef NDEBUG

      if (it == ctx->blockIndex.get<1>().end()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Block not allocated");
      }
      if (it->nb_rows != nb_r) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong size %d != %d", it->nb_rows, nb_r);
      }
      if (it->nb_cols != nb_c) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong size %d != %d", it->nb_cols, nb_c);
      }
      if (nb_r != mat.size1()) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong size %d != %d", nb_r, mat.size1());
      }
      if (nb_c != mat.size2()) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong size %d != %d", nb_c, mat.size2());
      }
      if (nb_r * nb_c != mat.data().size()) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong size %d != %d", nb_r * nb_c, mat.data().size());
      }

#endif // NDEBUG

      auto inv_shift = it->inv_shift;
      if (inv_shift != -1) {
        if (iora == ADD_VALUES) {
          cblas_daxpy(size, 1., &*mat.data().begin(), 1,
                      &(*ctx->preconditionerBlocksPtr)[inv_shift], 1);
        } else {
          std::copy(mat.data().begin(), mat.data().end(),
                    &(*ctx->preconditionerBlocksPtr)[inv_shift]);
        }
      } else {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "inv_shift == -1");
      }
    }
    MoFEMFunctionReturn(0);
  };

  if (nb_rows && nb_cols) {
    CHKERR set(0, 0, nb_rows, nb_cols, mat);
  }

#ifndef NDEBUG
  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_blockStrutureSetValues, 0, 0, 0, 0);
#endif // NDEBUG

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode shell_block_preconditioner_mat_asmb_wrap(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora) {
  MoFEMFunctionBegin;
  BlockStruture *ctx;
  CHKERR MatShellGetContext(M, (void **)&ctx);
  CHKERR shell_block_preconditioner_mat_asmb_wrap_impl(ctx, row_data, col_data,
                                                       mat, iora);
  MoFEMFunctionReturn(0);
}

boost::shared_ptr<NestSchurData> getNestSchurData(

    std::pair<SmartPetscObj<DM>, SmartPetscObj<DM>> dms,
    boost::shared_ptr<BlockStruture> block_mat_data_ptr,

    std::vector<std::string> fields_names, //< a00 fields
    std::vector<boost::shared_ptr<Range>>
        field_ents, //< a00 ranges (can be null),
    bool add_preconditioner_block) {

  if (!block_mat_data_ptr) {
    CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY, "Block data not set");
  }

  auto [schur_dm, block_dm] = dms;
  auto schur_prb = getProblemPtr(schur_dm);
  auto block_prb = getProblemPtr(block_dm);
  auto m_field_ptr = getInterfacePtr(block_dm);

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

  auto get_vec = [&](auto schur_data) {
    std::vector<int> vec_r, vec_c;
    vec_r.reserve(schur_data->blockIndex.size());
    vec_c.reserve(schur_data->blockIndex.size());
    for (auto &d : schur_data->blockIndex.template get<1>()) {
      vec_r.push_back(d.row);
      vec_c.push_back(d.col);
    }
    return std::make_pair(vec_r, vec_c);
  };

  auto [vec_r_schur, vec_c_schur] = get_vec(block_mat_data_ptr);
  CHKERR AOApplicationToPetsc(ao_schur_row, vec_r_schur.size(),
                              &*vec_r_schur.begin());
  CHKERR AOApplicationToPetsc(ao_schur_col, vec_c_schur.size(),
                              &*vec_c_schur.begin());
  auto [vec_r_block, vec_c_block] = get_vec(block_mat_data_ptr);
  CHKERR AOApplicationToPetsc(ao_block_row, vec_r_block.size(),
                              &*vec_r_block.begin());
  CHKERR AOApplicationToPetsc(ao_block_col, vec_c_block.size(),
                              &*vec_c_block.begin());

  std::array<boost::shared_ptr<BlockStruture>, 4> data_ptrs;

  for (auto r = 0; r != 3; ++r) {
    data_ptrs[r] = boost::make_shared<BlockStruture>();
    data_ptrs[r]->dataBlocksPtr = block_mat_data_ptr->dataBlocksPtr;
  }
  data_ptrs[3] = boost::make_shared<DiagBlockInvStruture>();
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
  int inv_mem_size = 0;
  for (auto &d : block_mat_data_ptr->blockIndex.get<1>()) {

    auto insert = [&](auto &m, auto &dof_r, auto &dof_c, auto &d) {
      m.insert(

          BlockStruture::Indexes{
              dof_r->getPetscGlobalDofIdx(), dof_c->getPetscGlobalDofIdx(),

              d.nb_rows, d.nb_cols,

              dof_r->getPetscLocalDofIdx(), dof_c->getPetscLocalDofIdx(),

              d.mat_shift, d.inv_shift}

      );
    };

    if (vec_r_schur[idx] != -1 && vec_c_schur[idx] != -1) {
      auto schur_dof_r =
          schur_dofs_row->get<PetscGlobalIdx_mi_tag>().find(vec_r_schur[idx]);
      auto schur_dof_c =
          schur_dofs_col->get<PetscGlobalIdx_mi_tag>().find(vec_c_schur[idx]);
#ifndef NDEBUG
      if (schur_dof_r == schur_dofs_row->get<PetscGlobalIdx_mi_tag>().end()) {
        CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                          "Block <Schur> not found");
      }
      if (schur_dof_c == schur_dofs_col->get<PetscGlobalIdx_mi_tag>().end()) {
        CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                          "Block <Schur> not found");
      }
#endif // NDEBUG
      insert(data_ptrs[0]->blockIndex, *schur_dof_r, *schur_dof_c, d);
    }

    if (vec_r_schur[idx] != -1 && vec_c_block[idx] != -1) {
      auto schur_dof_r =
          schur_dofs_row->get<PetscGlobalIdx_mi_tag>().find(vec_r_schur[idx]);
      auto block_dof_c =
          block_dofs_col->get<PetscGlobalIdx_mi_tag>().find(vec_c_block[idx]);
#ifndef NDEBUG
      if (schur_dof_r == schur_dofs_row->get<PetscGlobalIdx_mi_tag>().end()) {
        CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                          "Block <Schur> not found");
      }
      if (block_dof_c == block_dofs_col->get<PetscGlobalIdx_mi_tag>().end()) {
        CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                          "Block <Block> not found");
      }
#endif
      insert(data_ptrs[1]->blockIndex, *schur_dof_r, *block_dof_c, d);
    }

    if (vec_r_block[idx] != -1 && vec_c_schur[idx] != -1) {
      auto block_dof_r =
          block_dofs_row->get<PetscGlobalIdx_mi_tag>().find(vec_r_block[idx]);
      auto schur_dof_c =
          schur_dofs_col->get<PetscGlobalIdx_mi_tag>().find(vec_c_schur[idx]);
#ifndef NDEBUG
      if (block_dof_r == block_dofs_row->get<PetscGlobalIdx_mi_tag>().end()) {
        CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                          "Block <Block> not found");
      }
      if (schur_dof_c == schur_dofs_col->get<PetscGlobalIdx_mi_tag>().end()) {
        CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                          "Block <Schur> not found");
      }
#endif // NDEBUG
      insert(data_ptrs[2]->blockIndex, *block_dof_r, *schur_dof_c, d);
    }

    if (vec_r_block[idx] != -1 && vec_c_block[idx] != -1) {
      auto block_dof_r =
          block_dofs_row->get<PetscGlobalIdx_mi_tag>().find(vec_r_block[idx]);
      auto block_dof_c =
          block_dofs_col->get<PetscGlobalIdx_mi_tag>().find(vec_c_block[idx]);
      if (d.row == d.col && d.nb_rows == d.nb_cols) {
        // Only store inverse of diaconal blocks
        d.inv_shift = inv_mem_size;
        inv_mem_size += d.nb_cols * d.nb_rows;
      }
      insert(data_ptrs[3]->blockIndex, *block_dof_r, *block_dof_c, d);
    }

    ++idx;
  }

  auto set_up_a00_data = [&](auto inv_block_data) {
    MoFEMFunctionBegin;

    auto get_a00_uids = [&]() {
      auto get_field_bit = [&](auto field_name) {
        return m_field_ptr->get_field_bit_number(field_name);
      };

      std::vector<std::pair<UId, UId>> a00_uids;
      a00_uids.reserve(fields_names.size());
      for (auto ss = 0; ss != fields_names.size(); ++ss) {
        auto field_bit = get_field_bit(fields_names[ss]);
        auto row_ents = field_ents[ss];
        if (row_ents) {
          for (auto p = row_ents->pair_begin(); p != row_ents->pair_end();
               ++p) {
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

    auto get_glob_idex_pairs = [&](auto &&uid_pairs) {
      std::vector<std::pair<int, int>> glob_idex_pairs;
      glob_idex_pairs.reserve(uid_pairs.size());
      auto dofs = block_prb->getNumeredRowDofsPtr();

      auto it = uid_pairs.rbegin();
      auto hi = uid_pairs.rend();

      for (; it != hi; ++it) {
        auto [lo_uid, hi_uid] = *it;
        auto lo_it = dofs->lower_bound(lo_uid);
        auto hi_it = dofs->upper_bound(hi_uid);
        if (lo_it != hi_it) {
          auto lo_idx = (*lo_it)->getPetscGlobalDofIdx();
          glob_idex_pairs.emplace_back(lo_idx, std::distance(lo_it, hi_it));
        }
      }
      return glob_idex_pairs;
    };

    auto &index_view =
        boost::dynamic_pointer_cast<DiagBlockInvStruture>(inv_block_data)
            ->indexView;

    index_view.first.resize(0);
    index_view.second.resize(0);
    index_view.first.reserve(inv_block_data->blockIndex.size());
    index_view.second.reserve(inv_block_data->blockIndex.size());
    index_view.second.push_back(0);

    // this enable esrch by varying ranges
    using BlockIndexView = multi_index_container<

        const DiagBlockIndex::Indexes *,

        indexed_by<

            ordered_non_unique<member<DiagBlockIndex::Indexes, const int,
                                      &DiagBlockIndex::Indexes::row>>

            >>;

    BlockIndexView block_index_view;
    for (auto it = inv_block_data->blockIndex.begin();
         it != inv_block_data->blockIndex.end(); ++it) {
      block_index_view.insert(&*it);
    }

    auto glob_idx_pairs = get_glob_idex_pairs(get_a00_uids());

    // iterate field & entities pairs
    for (auto s1 = 0; s1 != glob_idx_pairs.size(); ++s1) {

      // iterate matrix indexes
      auto [lo_idx, nb] = glob_idx_pairs[s1];
      auto it = block_index_view.lower_bound(lo_idx);
      auto hi = block_index_view.end();

      // iterate rows
      while (it != hi && ((*it)->row + (*it)->nb_rows) <= lo_idx + nb) {

        auto row = (*it)->row;

        auto get_diag_index =
            [&](auto it) -> const MoFEM::DiagBlockIndex::Indexes * {
          while (it != hi && (*it)->row == row) {
            if ((*it)->col == row && (*it)->nb_cols == (*it)->nb_rows) {
              auto ptr = *it;
              return ptr;
            }
            ++it;
          }
          CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY, "Diagonal not found");
          return nullptr;
        };

        auto push_off_diag = [&](auto it, auto s1) {
          while (it != hi && (*it)->row == row) {
            for (int s2 = 0; s2 < s1; ++s2) {
              auto [col_lo_idx, col_nb] = glob_idx_pairs[s2];
              if ((*it)->col >= col_lo_idx &&
                  (*it)->col + (*it)->nb_cols <= col_lo_idx + col_nb) {
                index_view.first.push_back(*it);
              }
            }
            ++it;
          }
        };

        index_view.first.push_back(get_diag_index(it));
        push_off_diag(it, s1);
        index_view.second.push_back(index_view.first.size());

        while (it != hi && (*it)->row == row) {
          ++it;
        }
      }
    }

    auto get_vec = [&]() {
      std::vector<int> vec_r, vec_c;
      vec_r.reserve(index_view.first.size());
      vec_c.reserve(index_view.first.size());
      for (auto &i : index_view.first) {
        vec_r.push_back(i->row);
        vec_c.push_back(i->col);
      }
      return std::make_pair(vec_r, vec_c);
    };

    auto [vec_r_block, vec_c_block] = get_vec();
    CHKERR AOPetscToApplication(ao_block_row, vec_r_block.size(),
                                &*vec_r_block.begin());
    CHKERR AOPetscToApplication(ao_block_col, vec_c_block.size(),
                                &*vec_c_block.begin());

    int inv_mem_size = 0;
    for (auto s1 = 0; s1 != index_view.second.size() - 1; ++s1) {
      auto lo = index_view.second[s1];
      // auto hi = index_view.second[s1 + 1];

      auto diag_index_ptr = index_view.first[lo];
      auto row = vec_r_block[lo];
      auto col = vec_c_block[lo];
      auto nb_rows = diag_index_ptr->nb_rows;
      auto nb_cols = diag_index_ptr->nb_cols;
      ++lo;

      auto it = block_mat_data_ptr->blockIndex.get<1>().find(
          boost::make_tuple(row, col, nb_rows, nb_cols));
      if (it == block_mat_data_ptr->blockIndex.get<1>().end()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Block not found");
      }

      it->inv_shift = inv_mem_size;
      diag_index_ptr->inv_shift = it->inv_shift;
      inv_mem_size += nb_rows * nb_cols;

      // for (; lo != hi; ++lo) {
      // }
    }

    auto inv_data_ptr =
        boost::make_shared<std::vector<double>>(inv_mem_size, 0);
    data_ptrs[3]->dataInvBlocksPtr = inv_data_ptr;
    block_mat_data_ptr->dataInvBlocksPtr = data_ptrs[3]->dataInvBlocksPtr;

    if (add_preconditioner_block) {
      auto preconditioned_block =
          boost::make_shared<std::vector<double>>(inv_mem_size, 0);
      data_ptrs[3]->preconditionerBlocksPtr = preconditioned_block;
      data_ptrs[3]->multiplyByPreconditioner = true;
      block_mat_data_ptr->preconditionerBlocksPtr =
          data_ptrs[3]->preconditionerBlocksPtr;
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

  MOFEM_TAG_AND_LOG("SYNC", Sev::verbose, "NetsetSchur")
      << "(0, 0) " << schur_nb_local << " " << schur_nb_global << " "
      << data_ptrs[0]->blockIndex.size();
  MOFEM_TAG_AND_LOG("SYNC", Sev::verbose, "NetsetSchur")
      << "(0, 1) " << schur_nb_local << " " << block_nb_local << " "
      << schur_nb_global << " " << block_nb_global << " "
      << data_ptrs[1]->blockIndex.size();
  MOFEM_TAG_AND_LOG("SYNC", Sev::verbose, "NetsetSchur")
      << "(1, 0) " << block_nb_local << " " << schur_nb_local << " "
      << block_nb_global << " " << schur_nb_global << " "
      << data_ptrs[2]->blockIndex.size();
  MOFEM_TAG_AND_LOG("SYNC", Sev::verbose, "NetsetSchur")
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
                         std::vector<SmartPetscObj<AO>> sequence_of_aos,
                         std::vector<SmartPetscObj<Mat>> sequence_of_mats,
                         std::vector<bool> sym_schur, bool symm_op,
                         boost::shared_ptr<BlockStruture> diag_blocks) {
  if (symm_op)
    return new OpSchurAssembleEnd<SchurDSYSV>(fields_name, field_ents,
                                              sequence_of_aos, sequence_of_mats,
                                              sym_schur, symm_op, diag_blocks);
  else
    return new OpSchurAssembleEnd<SchurDGESV>(fields_name, field_ents,
                                              sequence_of_aos, sequence_of_mats,
                                              sym_schur, symm_op, diag_blocks);
}

OpSchurAssembleBase *
createOpSchurAssembleEnd(std::vector<std::string> fields_name,
                         std::vector<boost::shared_ptr<Range>> field_ents,
                         std::vector<SmartPetscObj<AO>> sequence_of_aos,
                         std::vector<SmartPetscObj<Mat>> sequence_of_mats,
                         std::vector<bool> sym_schur,
                         std::vector<double> diag_eps, bool symm_op,
                         boost::shared_ptr<BlockStruture> diag_blocks) {
  if (symm_op)
    return new OpSchurAssembleEnd<SchurDSYSV>(
        fields_name, field_ents, sequence_of_aos, sequence_of_mats, sym_schur,
        diag_eps, symm_op, diag_blocks);
  else
    return new OpSchurAssembleEnd<SchurDGESV>(
        fields_name, field_ents, sequence_of_aos, sequence_of_mats, sym_schur,
        diag_eps, symm_op, diag_blocks);
}

MoFEMErrorCode setSchurMatSolvePC(SmartPetscObj<PC> pc) {
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

// Assemble specialisations

// This is used to assemble Schur complement "S"
OpSchurAssembleBase::MatSetValuesRaw OpSchurAssembleBase::matSetValuesSchurRaw =
    ::MatSetValues;

// This used to assemble element assemble local matrices for Schur complement,
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
  CHKERR SchurBackendMatSetValuesPtr::matSetValuesPtr(M, row_data, col_data,
                                                      mat, iora);
  CHKERR assembleStorage(row_data, col_data, mat, iora);
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
MatSetValues<BlockStruture>(Mat M, const EntitiesFieldData::EntData &row_data,
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
  CHKERR SchurBackendMatSetValuesPtr::matSetValuesBlockPtr(M, row_data,
                                                           col_data, mat, iora);
  CHKERR assembleStorage(row_data, col_data, mat, iora);
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
  CHKERR SchurBackendMatSetValuesPtr::matSetValuesPreconditionedBlockPtr(
      M, row_data, col_data, mat, iora);
  CHKERR assembleStorage(row_data, col_data, mat, iora);
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
schurSwitchPreconditioner(boost::shared_ptr<BlockStruture> block_mat_data) {
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
schurSaveBlockMesh(boost::shared_ptr<BlockStruture> block_mat_data,
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

  double def_val_mat_shift = 0;
  Tag th_mat_shift;
  CHKERR moab.tag_get_handle("mat_shift", 1, MB_TYPE_INTEGER, th_mat_shift,
                             MB_TAG_CREAT | MB_TAG_DENSE, &def_val_mat_shift);

  int i = 0;
  for (auto &d : block_mat_data->blockIndex) {
    auto row = d.row;
    auto col = d.col;
    auto nb_rows = d.nb_rows;
    auto nb_cols = d.nb_cols;
    auto mat_shift = d.mat_shift;

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

    auto prt = &(*block_mat_data->dataBlocksPtr)[d.mat_shift];
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