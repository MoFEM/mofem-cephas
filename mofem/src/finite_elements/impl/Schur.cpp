/**
 * @file Schur.cpp
 * @brief Implementation of Schur Complement
 * @date 2023-02-03
 *
 * @copyright Copyright (c) 2023
 *
 */

namespace MoFEM {

PetscLogEvent SchurEvents::MOFEM_EVENT_schurL2MatsMatSetValues;
PetscLogEvent SchurEvents::MOFEM_EVENT_opSchurAssembleEnd;
PetscLogEvent SchurEvents::MOFEM_EVENT_diagBlockStrutureSetValues;
PetscLogEvent SchurEvents::MOFEM_EVENT_diagBlockStrutureMult;
PetscLogEvent SchurEvents::MOFEM_EVENT_zeroRowsAndCols;

SchurEvents::SchurEvents() {
  PetscLogEventRegister("schurMatSetVal", 0,
                        &MOFEM_EVENT_schurL2MatsMatSetValues);
  PetscLogEventRegister("opSchurAsmbEnd", 0, &MOFEM_EVENT_opSchurAssembleEnd);
  PetscLogEventRegister("diagBlockSetVal", 0,
                        &MOFEM_EVENT_diagBlockStrutureSetValues);
  PetscLogEventRegister("diagBlockMult", 0, &MOFEM_EVENT_diagBlockStrutureMult);
  PetscLogEventRegister("schurZeroRandC", 0, &MOFEM_EVENT_zeroRowsAndCols);
}

/**
 * @brief Schur complement data storage
 *
 */
struct SchurL2Mats : public boost::enable_shared_from_this<SchurL2Mats> {

  SchurL2Mats(const size_t idx, const UId uid_row, const UId uid_col);
  virtual ~SchurL2Mats() = default;

  const UId uidRow;
  const UId uidCol;

  inline auto &getMat() const { return locMats[iDX]; }
  inline auto &getRowInd() const { return rowIndices[iDX]; }
  inline auto &getColInd() const { return colIndices[iDX]; }

  using MatSetValuesPtr = boost::function<MoFEMErrorCode(
      Mat M, const EntitiesFieldData::EntData &row_data,
      const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
      InsertMode iora)>;

  static MatSetValuesPtr matSetValuesPtr; ///< backend assembly function

  static MoFEMErrorCode MatSetValues(Mat M,
                                     const EntitiesFieldData::EntData &row_data,
                                     const EntitiesFieldData::EntData &col_data,
                                     const MatrixDouble &mat, InsertMode iora);

private:
  const size_t iDX;

  friend OpSchurAssembleBegin;
  friend OpSchurAssembleEndImpl;

  struct idx_mi_tag {};
  struct uid_mi_tag {};
  struct row_mi_tag {};
  struct col_mi_tag {};

  using SchurL2Storage = multi_index_container<
      SchurL2Mats,
      indexed_by<

          ordered_unique<
              tag<uid_mi_tag>,
              composite_key<
                  SchurL2Mats,

                  member<SchurL2Mats, const UId, &SchurL2Mats::uidRow>,
                  member<SchurL2Mats, const UId, &SchurL2Mats::uidCol>

                  >>,

          ordered_non_unique<tag<row_mi_tag>, member<SchurL2Mats, const UId,
                                                     &SchurL2Mats::uidRow>>,

          ordered_non_unique<tag<col_mi_tag>, member<SchurL2Mats, const UId,
                                                     &SchurL2Mats::uidCol>>

          >>;

  static boost::ptr_vector<MatrixDouble> locMats;
  static boost::ptr_vector<VectorInt> rowIndices;
  static boost::ptr_vector<VectorInt> colIndices;
  static SchurL2Storage schurL2Storage;

  static SmartPetscObj<Mat> a00;
  static SmartPetscObj<Mat> a01;
  static SmartPetscObj<Mat> a10;
  static SmartPetscObj<Mat> a11;

  static SmartPetscObj<Mat> ao0;
  static SmartPetscObj<Mat> ao1;

  static std::vector<int> rowIndices0;
  static std::vector<int> colIndices0;
  static std::vector<int> rowIndices1;
  static std::vector<int> colIndices1;
};

template <>
inline MoFEMErrorCode
MatSetValues<SchurL2Mats>(Mat M, const EntitiesFieldData::EntData &row_data,
                          const EntitiesFieldData::EntData &col_data,
                          const MatrixDouble &mat, InsertMode iora) {
  return SchurL2Mats::MatSetValues(M, row_data, col_data, mat, iora);
}

template <>
MoFEMErrorCode MatSetValues<AssemblyTypeSelector<SCHUR>>(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora) {
  return MatSetValues<SchurL2Mats>(M, row_data, col_data, mat, iora);
}

SchurL2Mats::SchurL2Storage SchurL2Mats::schurL2Storage;
boost::ptr_vector<MatrixDouble> SchurL2Mats::locMats;
boost::ptr_vector<VectorInt> SchurL2Mats::rowIndices;
boost::ptr_vector<VectorInt> SchurL2Mats::colIndices;

SmartPetscObj<Mat> SchurL2Mats::a00;
SmartPetscObj<Mat> SchurL2Mats::a01;
SmartPetscObj<Mat> SchurL2Mats::a10;
SmartPetscObj<Mat> SchurL2Mats::a11;

SmartPetscObj<Mat> SchurL2Mats::ao0;
SmartPetscObj<Mat> SchurL2Mats::ao1;

std::vector<int> SchurL2Mats::rowIndices0;
std::vector<int> SchurL2Mats::colIndices0;
std::vector<int> SchurL2Mats::rowIndices1;
std::vector<int> SchurL2Mats::colIndices1;

SchurL2Mats::SchurL2Mats(const size_t idx, const UId uid_row, const UId uid_col)
    : iDX(idx), uidRow(uid_row), uidCol(uid_col) {}

OpSchurAssembleEndImpl::MatSetValuesRaw
    OpSchurAssembleEndImpl::matSetValuesSchurRaw = ::MatSetValues;

MoFEMErrorCode
schur_mat_set_values_wrap(Mat M, const EntitiesFieldData::EntData &row_data,
                          const EntitiesFieldData::EntData &col_data,
                          const MatrixDouble &mat, InsertMode iora) {
  return MatSetValues<AssemblyTypeSelector<PETSC>>(M, row_data, col_data, mat,
                                                   iora);
}

SchurL2Mats::MatSetValuesPtr SchurL2Mats::matSetValuesPtr =
    schur_mat_set_values_wrap;

MoFEMErrorCode
SchurL2Mats::MatSetValues(Mat M, const EntitiesFieldData::EntData &row_data,
                          const EntitiesFieldData::EntData &col_data,
                          const MatrixDouble &mat, InsertMode iora) {
  MoFEMFunctionBegin;

  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_schurL2MatsMatSetValues, 0, 0, 0,
                     0);

  // PetscBool is_nested = PETSC_FALSE;
  // PetscObjectTypeCompare((PetscObject)M, MATNEST, &is_nested);
  // if (is_nested) {
  //   CHKERR matSetValuesNestedPtr(M, row_data, col_data, mat, iora);
  // } else {
  CHKERR matSetValuesPtr(M, row_data, col_data, mat, iora);
  // }

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

  const auto idx = SchurL2Mats::schurL2Storage.size();
  const auto size = SchurL2Mats::locMats.size();

  if (idx >= size) {
    SchurL2Mats::locMats.push_back(new MatrixDouble());
    SchurL2Mats::rowIndices.push_back(new VectorInt());
    SchurL2Mats::colIndices.push_back(new VectorInt());
  }

  // insert index
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
  auto p = SchurL2Mats::schurL2Storage.emplace(idx, uid_row, uid_col);

  auto get_storage = [&p]() { return const_cast<SchurL2Mats &>(*p.first); };

  if (p.second) {

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
  }

  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_schurL2MatsMatSetValues, 0, 0, 0,
                   0);

  MoFEMFunctionReturn(0);
}

OpSchurAssembleBegin::OpSchurAssembleBegin()
    : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPSPACE) {}

MoFEMErrorCode OpSchurAssembleBegin::doWork(int side, EntityType type,
                                            EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
#ifndef NDEBUG
  MOFEM_LOG("SELF", Sev::noisy) << "Schur assemble begin";
#endif
  SchurL2Mats::schurL2Storage.clear();

  MoFEMFunctionReturn(0);
}

OpSchurAssembleEndImpl::OpSchurAssembleEndImpl(
    std::vector<std::string> fields_name,
    std::vector<boost::shared_ptr<Range>> field_ents,
    std::vector<SmartPetscObj<AO>> sequence_of_aos,
    std::vector<SmartPetscObj<Mat>> sequence_of_mats,
    std::vector<bool> sym_schur, std::vector<double> diag_eps, bool symm_op)
    : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPSPACE, symm_op),
      fieldsName(fields_name), fieldEnts(field_ents),
      sequenceOfAOs(sequence_of_aos), sequenceOfMats(sequence_of_mats),
      symSchur(sym_schur), diagEps(diag_eps) {}

OpSchurAssembleEndImpl::OpSchurAssembleEndImpl(
    std::vector<std::string> fields_name,
    std::vector<boost::shared_ptr<Range>> field_ents,
    std::vector<SmartPetscObj<AO>> sequence_of_aos,
    std::vector<SmartPetscObj<Mat>> sequence_of_mats,
    std::vector<bool> sym_schur, bool symm_op)
    : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPSPACE, symm_op),
      fieldsName(fields_name), fieldEnts(field_ents),
      sequenceOfAOs(sequence_of_aos), sequenceOfMats(sequence_of_mats),
      symSchur(sym_schur), diagEps(fields_name.size(), 0) {}

template <typename I>
MoFEMErrorCode
OpSchurAssembleEndImpl::doWorkImpl(int side, EntityType type,
                                   EntitiesFieldData::EntData &data) {

  MoFEMFunctionBegin;

  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_opSchurAssembleEnd, 0, 0, 0, 0);

#ifndef NDEBUG
  MOFEM_LOG("SELF", Sev::noisy) << "Schur assemble begin -> end";
#endif

  auto assemble = [&](SmartPetscObj<Mat> M, MatSetValuesRaw mat_set_values,
                      auto &storage) {
    MoFEMFunctionBegin;
    if (M) {
      for (auto &s : storage) {
        auto &m = s.getMat();
        if (m.size1()) {
          auto &row_ind = s.getRowInd();
          auto &col_ind = s.getColInd();

          if (mat_set_values(M, row_ind.size(), &*row_ind.begin(),
                             col_ind.size(), &*col_ind.begin(),
                             &*m.data().begin(), ADD_VALUES)) {
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
            CHK_THROW_MESSAGE(ierr, "MatSetValues");
          }
        }
      }
    }
    MoFEMFunctionReturn(0);
  };

  auto apply_schur = [&](auto &storage,

                         const auto ss,

                         auto lo_uid, auto hi_uid, double eps) {
    MoFEMFunctionBegin;

    auto add_off_mat = [&](auto row_uid, auto col_uid, auto &row_ind,
                           auto &col_ind, auto &offMatInvDiagOffMat) {
      MoFEMFunctionBegin;

      const auto idx = SchurL2Mats::schurL2Storage.size();
      const auto size = SchurL2Mats::locMats.size();

      if (idx >= size) {
        SchurL2Mats::locMats.push_back(new MatrixDouble());
        SchurL2Mats::rowIndices.push_back(new VectorInt());
        SchurL2Mats::colIndices.push_back(new VectorInt());
      }

      auto it = storage.template get<SchurL2Mats::uid_mi_tag>().find(
          boost::make_tuple(row_uid, col_uid));

      if (it == storage.template get<SchurL2Mats::uid_mi_tag>().end()) {

        auto p = SchurL2Mats::schurL2Storage.emplace(idx, row_uid, col_uid);
        auto &mat = p.first->getMat();
        auto &set_row_ind = p.first->getRowInd();
        auto &set_col_ind = p.first->getColInd();

        set_row_ind.resize(row_ind.size(), false);
        noalias(set_row_ind) = row_ind;
        set_col_ind.resize(col_ind.size(), false);
        noalias(set_col_ind) = col_ind;

        mat.swap(offMatInvDiagOffMat);

      } else {
        auto &mat = it->getMat();
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
          storage.template get<SchurL2Mats::row_mi_tag>().lower_bound(lo_uid);
      auto hi_row_it =
          storage.template get<SchurL2Mats::row_mi_tag>().upper_bound(hi_uid);
      std::vector<const SchurL2Mats *> schur_row_ptr_view;
      schur_row_ptr_view.reserve(std::distance(row_it, hi_row_it));
      for (; row_it != hi_row_it; ++row_it) {
        schur_row_ptr_view.push_back(&*row_it);
      }
      return schur_row_ptr_view;
    };

    auto get_col_view = [&]() {
      auto col_it =
          storage.template get<SchurL2Mats::col_mi_tag>().lower_bound(lo_uid);
      auto hi_col_it =
          storage.template get<SchurL2Mats::col_mi_tag>().upper_bound(hi_uid);
      std::vector<const SchurL2Mats *> schur_col_ptr_view;
      schur_col_ptr_view.reserve(std::distance(col_it, hi_col_it));
      for (; col_it != hi_col_it; ++col_it) {
        schur_col_ptr_view.push_back(&*col_it);
      }
      return schur_col_ptr_view;
    };

    auto schur_row_ptr_view = get_row_view();
    auto schur_col_ptr_view = get_col_view();

    for (auto row_it : schur_row_ptr_view) {
      if (row_it->uidRow == row_it->uidCol) {

        CHKERR I::invertMat(row_it->getMat(), invMat, eps);

        for (auto c_lo : schur_col_ptr_view) {

          auto &row_uid = c_lo->uidRow;
          if (row_uid == row_it->uidRow)
            continue;

          auto &cc_off_mat = c_lo->getMat();
          invDiagOffMat.resize(cc_off_mat.size1(), invMat.size2(), false);
#ifndef NDEBUG
          if (invMat.size1() != cc_off_mat.size2()) {
            SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                     "Wrong size %d != %d", invMat.size1(), cc_off_mat.size2());
          }
#endif // NDEBUG

          // noalias(invDiagOffMat) = prod(cc_off_mat, invMat);
          cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                      cc_off_mat.size1(), invMat.size2(), cc_off_mat.size2(),
                      1., &*cc_off_mat.data().begin(), cc_off_mat.size2(),
                      &*invMat.data().begin(), invMat.size2(), 0.,
                      &*invDiagOffMat.data().begin(), invDiagOffMat.size2());

          for (auto r_lo : schur_row_ptr_view) {
            auto &col_uid = r_lo->uidCol;

            // Skip diagonal
            if (col_uid == row_it->uidRow)
              continue;

            // If symmetry only run upper off diagonal terms
            if (symSchur[ss] && row_uid > col_uid)
              continue;

            auto &rr_off_mat = r_lo->getMat();
            offMatInvDiagOffMat.resize(invDiagOffMat.size1(),
                                       rr_off_mat.size2(), false);
#ifndef NDEBUG
            if (rr_off_mat.size1() != invDiagOffMat.size2()) {
              SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                       "Wrong size %d != %d", rr_off_mat.size1(),
                       invDiagOffMat.size2());
            }
#endif // NDEBUG

            // noalias(offMatInvDiagOffMat) = prod(invDiagOffMat, rr_off_mat);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        invDiagOffMat.size1(), rr_off_mat.size2(),
                        invDiagOffMat.size2(), 1.,
                        &*invDiagOffMat.data().begin(), invDiagOffMat.size2(),
                        &*rr_off_mat.data().begin(), rr_off_mat.size2(), 0.,
                        &*offMatInvDiagOffMat.data().begin(),
                        offMatInvDiagOffMat.size2());

            CHKERR add_off_mat(row_uid, col_uid, c_lo->getRowInd(),
                               r_lo->getColInd(), offMatInvDiagOffMat);

            if (symSchur[ss] && row_uid != col_uid) {
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

  auto erase_factored = [&](auto &storage,

                            const auto ss, auto lo_uid, auto hi_uid) {
    MoFEMFunctionBegin;

    auto r_lo =
        storage.template get<SchurL2Mats::row_mi_tag>().lower_bound(lo_uid);
    auto r_hi =
        storage.template get<SchurL2Mats::row_mi_tag>().upper_bound(hi_uid);
    storage.template get<SchurL2Mats::row_mi_tag>().erase(r_lo, r_hi);

    auto c_lo =
        storage.template get<SchurL2Mats::col_mi_tag>().lower_bound(lo_uid);
    auto c_hi =
        storage.template get<SchurL2Mats::col_mi_tag>().upper_bound(hi_uid);
    storage.template get<SchurL2Mats::col_mi_tag>().erase(c_lo, c_hi);

    MoFEMFunctionReturn(0);
  };

  auto assemble_off_diagonal = [&](auto &storage, const auto ss) {
    MoFEMFunctionBegin;

    // apply AO
    for (auto &m : storage) {
      if (auto ao = sequenceOfAOs[ss]) {
        auto &ind_row = m.getRowInd();
        CHKERR AOApplicationToPetsc(ao, ind_row.size(), &*ind_row.begin());
        auto &ind_col = m.getColInd();
        CHKERR AOApplicationToPetsc(ao, ind_col.size(), &*ind_col.begin());
      }
    }

    // assemble Schur
    if (sequenceOfMats[ss]) {
      CHKERR assemble(sequenceOfMats[ss], matSetValuesSchurRaw, storage);
    }

    MoFEMFunctionReturn(0);
  };

  auto &storage = SchurL2Mats::schurL2Storage;
  // Assemble to global matrix

  // Assemble Schur complements
  // Iterate complements
  for (auto ss = 0; ss != fieldsName.size(); ++ss) {

    auto assemble = [&](auto &storage, auto ss, auto lo_uid, auto hi_uid) {
      MoFEMFunctionBegin;

      CHKERR apply_schur(storage, ss, lo_uid, hi_uid, diagEps[ss]);
      CHKERR erase_factored(storage, ss, lo_uid, hi_uid);
      CHKERR assemble_off_diagonal(storage, ss);
      MoFEMFunctionReturn(0);
    };

    auto field_bit = getPtrFE()->mField.get_field_bit_number(fieldsName[ss]);
    auto row_ents = fieldEnts[ss];

    if (row_ents) {
      for (auto p = row_ents->pair_begin(); p != row_ents->pair_end(); ++p) {
        auto lo_uid =
            FieldEntity::getLoLocalEntityBitNumber(field_bit, p->first);
        auto hi_uid =
            FieldEntity::getHiLocalEntityBitNumber(field_bit, p->second);
        CHKERR assemble(storage, ss, lo_uid, hi_uid);
      }
    } else {
      // Iterate all entities (typically L2)
      auto lo_uid = FieldEntity::getLoLocalEntityBitNumber(
          field_bit, get_id_for_min_type<MBVERTEX>());
      auto hi_uid = FieldEntity::getHiLocalEntityBitNumber(
          field_bit, get_id_for_max_type<MBENTITYSET>());
      CHKERR assemble(storage, ss, lo_uid, hi_uid);
    }
  }

#ifndef NDEBUG
  MOFEM_LOG("SELF", Sev::noisy) << "Schur assemble done";
#endif

  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_opSchurAssembleEnd, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

struct SCHUR_DSYSV {
  static auto invertMat(MatrixDouble &m, MatrixDouble &inv, double eps) {
    MoFEMFunctionBeginHot;
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

struct SCHUR_DGESV {
  static auto invertMat(MatrixDouble &m, MatrixDouble &inv, double eps) {
    MoFEMFunctionBeginHot;
    VectorInt ipiv;
    const auto nb = m.size1();
#ifndef NDEBUG
    if (nb != m.size2()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "It should be square matrix");
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
OpSchurAssembleEnd<SCHUR_DSYSV>::doWork(int side, EntityType type,
                                        EntitiesFieldData::EntData &data) {
  return doWorkImpl<SCHUR_DSYSV>(side, type, data);
}

MoFEMErrorCode
OpSchurAssembleEnd<SCHUR_DGESV>::doWork(int side, EntityType type,
                                        EntitiesFieldData::EntData &data) {
  return doWorkImpl<SCHUR_DGESV>(side, type, data);
}

struct DiagBlockStruture {

  DiagBlockStruture() = default;
  virtual ~DiagBlockStruture() = default;

  SmartPetscObj<Vec> ghostX;
  SmartPetscObj<Vec> ghostY;

  boost::shared_ptr<std::vector<double>> dataBlocksPtr;
  boost::shared_ptr<std::vector<double>> dataInvBlocksPtr;

  struct Indexes {
    int row;
    int col;
    int nb_rows;
    int nb_cols;
    int loc_row;
    int loc_col;
    mutable int shift;
  };

  using BlockIndex = multi_index_container<

      Indexes,

      indexed_by<

          ordered_non_unique<member<Indexes, int, &Indexes::col>>,

          hashed_unique<

              composite_key<Indexes,

                            member<Indexes, int, &Indexes::row>,
                            member<Indexes, int, &Indexes::col>,
                            member<Indexes, int, &Indexes::nb_rows>,
                            member<Indexes, int, &Indexes::nb_cols>

                            >>

          >>;

  boost::shared_ptr<BlockIndex> blockIndexPtr;
};

struct CountBlocks : public ForcesAndSourcesCore::UserDataOperator {
  using OP = ForcesAndSourcesCore::UserDataOperator;
  CountBlocks(std::string fr, std::string fc,
              boost::shared_ptr<DiagBlockStruture> data_ptr)
      : OP(fr, fc, OP::OPROWCOL), dataPtr(data_ptr) {
    sYmm = false;
  }
  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        EntitiesFieldData::EntData &row_data,
                        EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;

    auto nb_rows = row_data.getIndices().size();
    auto nb_cols = col_data.getIndices().size();

    auto add = [&](auto r, auto c, int nb_r, int nb_c) {
      dataPtr->blockIndexPtr->insert(

          DiagBlockStruture::Indexes{
              row_data.getIndices()[r], col_data.getIndices()[c], nb_r, nb_c,
              row_data.getLocalIndices()[r], col_data.getLocalIndices()[c], -1}

      );
    };

    if (nb_rows && nb_cols) {
      if (row_type == MBVERTEX && col_type == MBVERTEX) {
        auto row_nb_coeff = row_data.getFieldDofs()[0]->getNbOfCoeffs();
        auto col_nb_coeff = col_data.getFieldDofs()[0]->getNbOfCoeffs();
        for (auto i = 0; i < nb_rows / row_nb_coeff; ++i) {
          for (auto j = 0; j < nb_cols / col_nb_coeff; ++j) {
            add(i * row_nb_coeff, j * col_nb_coeff, row_nb_coeff, col_nb_coeff);
          }
        }
      } else if (row_type == MBVERTEX) {
        auto row_nb_coeff = row_data.getFieldDofs()[0]->getNbOfCoeffs();
        for (auto i = 0; i < nb_rows / row_nb_coeff; ++i) {
          add(i * row_nb_coeff, 0, row_nb_coeff, nb_cols);
        }
      } else if (col_type == MBVERTEX) {
        auto col_nb_coeff = col_data.getFieldDofs()[0]->getNbOfCoeffs();
        for (auto j = 0; j < nb_cols / col_nb_coeff; ++j) {
          add(0, j * col_nb_coeff, nb_rows, col_nb_coeff);
        }
      } else {
        add(0, 0, nb_rows, nb_cols);
      }
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<DiagBlockStruture> dataPtr;
};

boost::shared_ptr<DiagBlockStruture> createSchurBlockMatStructure(

    DM dm, //< dm

    SchurFEOpVec schur_fe_op_vec //< block elements

) {
  auto data_ptr = boost::make_shared<DiagBlockStruture>();
  data_ptr->blockIndexPtr = boost::make_shared<DiagBlockStruture::BlockIndex>();

  for (auto &d : schur_fe_op_vec) {

    auto fe_name = d.first.first;
    auto fe_ptr = d.first.second;

    for (auto &f : d.second) {
      fe_ptr->getOpPtrVector().push_back(
          new CountBlocks(f.first, f.second, data_ptr));
    }

    CHKERR DMoFEMLoopFiniteElements(dm, fe_name, fe_ptr);
  }

  auto mem_size = 0;
  for (auto &v : data_ptr->blockIndexPtr->get<0>()) {
    v.shift = mem_size;
    mem_size += v.nb_cols * v.nb_rows;
  }
  data_ptr->dataBlocksPtr =
      boost::make_shared<std::vector<double>>(mem_size, 0);
  data_ptr->dataInvBlocksPtr =
      boost::make_shared<std::vector<double>>(mem_size, 0);

  auto ghost_x = createDMVector(dm);
  auto ghost_y = createDMVector(dm);

  data_ptr->ghostX = ghost_x;
  data_ptr->ghostY = ghost_y;

  return data_ptr;
}

static MoFEMErrorCode mult_schur_block_shell(Mat mat, Vec x, Vec y,
                                             InsertMode iora, bool solve);

static PetscErrorCode mult(Mat mat, Vec x, Vec y) {
  return mult_schur_block_shell(mat, x, y, INSERT_VALUES, false);
}
static PetscErrorCode mult_add(Mat mat, Vec x, Vec y) {
  return mult_schur_block_shell(mat, x, y, ADD_VALUES, false);
}
static PetscErrorCode solve(Mat mat, Vec x, Vec y) {
  return mult_schur_block_shell(mat, x, y, INSERT_VALUES, true);
}
static PetscErrorCode solve_add(Mat mat, Vec x, Vec y) {
  return mult_schur_block_shell(mat, x, y, ADD_VALUES, true);
}

static PetscErrorCode zero_rows_columns(Mat A, PetscInt N,
                                        const PetscInt rows[], PetscScalar diag,
                                        Vec x, Vec b) {

  MoFEMFunctionBeginHot;
  DiagBlockStruture *ctx;
  CHKERR MatShellGetContext(A, (void **)&ctx);

  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_zeroRowsAndCols, 0, 0, 0, 0);

  int loc_m, loc_n;
  CHKERR MatGetLocalSize(A, &loc_m, &loc_n);

  struct ShiftedBlockView : DiagBlockStruture::Indexes {
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

  BlockIndexView view;
  for (auto &v : ctx->blockIndexPtr->get<0>()) {
    view.insert(static_cast<const ShiftedBlockView *>(&v));
  }

  for (auto n = 0; n != N; ++n) {
    auto row = rows[n];
    auto rlo = view.get<0>().lower_bound(row);
    auto rhi = view.get<0>().upper_bound(
        row + MAX_DOFS_ON_ENTITY); // add such that upper bound is included
    for (; rlo != rhi; ++rlo) {
      auto r_shift = row - (*rlo)->row;
      if (r_shift >= 0 && r_shift < (*rlo)->nb_rows) {
        auto *ptr = &(*ctx->dataBlocksPtr)[(*rlo)->shift];
        for (auto i = 0; i != (*rlo)->nb_cols; ++i) {
          ptr[i + r_shift * (*rlo)->nb_cols] = 0;
        }
      }
    }
  }

  for (auto n = 0; n != N; ++n) {
    auto col = rows[n];
    auto clo = view.get<1>().lower_bound(col);
    auto chi = view.get<1>().upper_bound(
        col + MAX_DOFS_ON_ENTITY); // add such that upper bound is included
    for (; clo != chi; ++clo) {
      auto c_shift = col - (*clo)->col;
      if (c_shift >= 0 && c_shift < (*clo)->nb_cols) {

        auto *ptr = &(*ctx->dataBlocksPtr)[(*clo)->shift];
        for (auto i = 0; i != (*clo)->nb_rows; ++i) {
          ptr[c_shift + i * (*clo)->nb_cols] = 0;
        }

        // diagonal
        if ((*clo)->row == (*clo)->col && (*clo)->loc_row < loc_m) {
          auto r_shift = col - (*clo)->row;
          if (r_shift >= 0 && r_shift < (*clo)->nb_rows) {
            ptr[c_shift + r_shift * (*clo)->nb_cols] = diag;
          }
        }
      }
    }
  }

  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_zeroRowsAndCols, 0, 0, 0, 0);

  MoFEMFunctionReturnHot(0);
}

static PetscErrorCode mat_zero(Mat m) {
  MoFEMFunctionBegin;
  DiagBlockStruture *ctx;
  CHKERR MatShellGetContext(m, (void **)&ctx);
  ctx->dataBlocksPtr->clear();
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

SchurShellMatData
createSchurBlockMat(DM dm, boost::shared_ptr<DiagBlockStruture> data) {

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

  return std::make_pair(SmartPetscObj<Mat>(mat_raw), data);
}

static MoFEMErrorCode mult_schur_block_shell(Mat mat, Vec x, Vec y,
                                             InsertMode iora, bool solve) {
  MoFEMFunctionBegin;
  DiagBlockStruture *ctx;
  CHKERR MatShellGetContext(mat, (void **)&ctx);

  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_diagBlockStrutureMult, 0, 0, 0,
                     0);

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

  double *block_ptr;
  if (solve)
    block_ptr = &*ctx->dataInvBlocksPtr->begin();
  else
    block_ptr = &*ctx->dataBlocksPtr->begin();

  std::vector<int> loc_rows;
  std::vector<int> loc_nb_rows;
  std::vector<double> loc_y;

  auto it = ctx->blockIndexPtr->get<0>().lower_bound(0);
  auto hi = ctx->blockIndexPtr->get<0>().end();

  for (; it != hi;) {

    auto shift = it->shift;

    auto loc_col = it->loc_col;
    auto nb_cols = it->nb_cols;
    auto x_ptr = &x_array[loc_col];

    int it_shift = shift;
    loc_rows.resize(0);
    loc_nb_rows.resize(0);
    loc_nb_rows.push_back(0);
    while (it != hi && (it->loc_col == loc_col && it->nb_cols == nb_cols &&
                        it->shift == it_shift)) {
      loc_rows.push_back(it->loc_row);
      loc_nb_rows.push_back(loc_nb_rows.back() + it->nb_rows);
      it_shift += it->nb_cols * it->nb_rows;
      ++it;
    }
    loc_y.resize(loc_nb_rows.back());

    cblas_dgemv(

        CblasRowMajor, CblasNoTrans,

        loc_nb_rows.back(), nb_cols,

        1., &(block_ptr[shift]), nb_cols,

        x_ptr, 1,

        0., &*loc_y.begin(), 1

    );

    for (auto r = 0; r != loc_rows.size(); ++r) {
      cblas_daxpy(loc_nb_rows[r + 1] - loc_nb_rows[r], 1.,
                  &loc_y[loc_nb_rows[r]], 1, &y_array[loc_rows[r]], 1);
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
  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_diagBlockStrutureMult, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
shell_schur_mat_set_values_wrap(Mat M,
                                const EntitiesFieldData::EntData &row_data,
                                const EntitiesFieldData::EntData &col_data,
                                const MatrixDouble &mat, InsertMode iora) {

  MatrixDouble tmp_mat;
  MoFEMFunctionBegin;

  DiagBlockStruture *ctx;
  CHKERR MatShellGetContext(M, (void **)&ctx);

  PetscLogEventBegin(SchurEvents::MOFEM_EVENT_diagBlockStrutureSetValues, 0, 0,
                     0, 0);

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
      auto it = ctx->blockIndexPtr->get<1>().find(boost::make_tuple(
          row_first_index, col_data.getIndices()[c], nb_r, nb_c));

#ifndef NDEBUG

      if (it == ctx->blockIndexPtr->get<1>().end()) {
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

      auto shift = it->shift;
      if (shift != -1) {

        if (iora == ADD_VALUES) {
          cblas_daxpy(size, 1., &*mat.data().begin(), 1,
                      &(*ctx->dataBlocksPtr)[shift], 1);
        } else {
          std::copy(mat.data().begin(), mat.data().end(),
                    &(*ctx->dataBlocksPtr)[shift]);
        }
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

  PetscLogEventEnd(SchurEvents::MOFEM_EVENT_diagBlockStrutureSetValues, 0, 0, 0,
                   0);

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode
MatSetValues<DiagBlockStruture>(Mat M,
                                const EntitiesFieldData::EntData &row_data,
                                const EntitiesFieldData::EntData &col_data,
                                const MatrixDouble &mat, InsertMode iora) {
  return shell_schur_mat_set_values_wrap(M, row_data, col_data, mat, iora);
}

SchurNestMatrixData
getSchurNestMatArray(std::pair<SmartPetscObj<DM>, SmartPetscObj<DM>> dms,
                     SchurShellMatData A) {
  auto [schur_dm, block_dm] = dms;
  auto schur_prb = getProblemPtr(schur_dm);
  auto block_prb = getProblemPtr(block_dm);

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

  auto [schur_mat, schur_data] = A;

  auto get_vec = [&]() {
    std::vector<int> vec_r, vec_c;
    vec_r.reserve(schur_data->blockIndexPtr->size());
    vec_c.reserve(schur_data->blockIndexPtr->size());
    for (auto &d : *schur_data->blockIndexPtr) {
      vec_r.push_back(d.row);
      vec_c.push_back(d.col);
    }
    return std::make_pair(vec_r, vec_c);
  };

  auto [vec_r_schur, vec_c_schur] = get_vec();
  CHKERR AOApplicationToPetsc(ao_schur_row, vec_r_schur.size(),
                              &*vec_r_schur.begin());
  CHKERR AOApplicationToPetsc(ao_schur_col, vec_c_schur.size(),
                              &*vec_c_schur.begin());
  auto [vec_r_block, vec_c_block] = get_vec();
  CHKERR AOApplicationToPetsc(ao_block_row, vec_r_block.size(),
                              &*vec_r_block.begin());
  CHKERR AOApplicationToPetsc(ao_block_col, vec_c_block.size(),
                              &*vec_c_block.begin());

  std::array<boost::shared_ptr<DiagBlockStruture>, 4> data_ptrs;

  for (auto r = 0; r != 4; ++r) {
    data_ptrs[r] = boost::make_shared<DiagBlockStruture>();
    data_ptrs[r]->blockIndexPtr =
        boost::make_shared<DiagBlockStruture::BlockIndex>();
    data_ptrs[r]->dataBlocksPtr = schur_data->dataBlocksPtr;
  }

  data_ptrs[0]->ghostX = schur_vec_x;
  data_ptrs[0]->ghostY = schur_vec_y;
  data_ptrs[1]->ghostX = block_vec_x;
  data_ptrs[1]->ghostY = schur_vec_y;
  data_ptrs[2]->ghostX = schur_vec_x;
  data_ptrs[2]->ghostY = block_vec_y;
  data_ptrs[3]->ghostX = block_vec_x;
  data_ptrs[3]->ghostY = block_vec_y;

  int idx = 0;
  for (auto &d : *schur_data->blockIndexPtr) {

    auto insert = [&](auto &m_ptr, auto &dof_r, auto &dof_c, auto &d) {
      m_ptr->insert(

          DiagBlockStruture::Indexes{dof_r->getPetscGlobalDofIdx(),
                                     dof_c->getPetscGlobalDofIdx(),

                                     d.nb_rows, d.nb_cols,

                                     dof_r->getPetscLocalDofIdx(),
                                     dof_c->getPetscLocalDofIdx(), d.shift}

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
                          "Block <Schur> not not fund");
      }
      if (schur_dof_c == schur_dofs_col->get<PetscGlobalIdx_mi_tag>().end()) {
        CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                          "Block <Schur> not not fund");
      }
#endif // NDEBUG
      insert(data_ptrs[0]->blockIndexPtr, *schur_dof_r, *schur_dof_c, d);
    }

    if (vec_r_schur[idx] != -1 && vec_c_block[idx] != -1) {
      auto schur_dof_r =
          schur_dofs_row->get<PetscGlobalIdx_mi_tag>().find(vec_r_schur[idx]);
      auto block_dof_c =
          block_dofs_col->get<PetscGlobalIdx_mi_tag>().find(vec_c_block[idx]);
#ifndef NDEBUG
      if (schur_dof_r == schur_dofs_row->get<PetscGlobalIdx_mi_tag>().end()) {
        CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                          "Block <Schur> not not fund");
      }
      if (block_dof_c == block_dofs_col->get<PetscGlobalIdx_mi_tag>().end()) {
        CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                          "Block <Block> not not fund");
      }
#endif
      insert(data_ptrs[1]->blockIndexPtr, *schur_dof_r, *block_dof_c, d);
    }

    if (vec_r_block[idx] != -1 && vec_c_schur[idx] != -1) {
      auto block_dof_r =
          block_dofs_row->get<PetscGlobalIdx_mi_tag>().find(vec_r_block[idx]);
      auto schur_dof_c =
          schur_dofs_col->get<PetscGlobalIdx_mi_tag>().find(vec_c_schur[idx]);
#ifndef NDEBUG
      if (block_dof_r == block_dofs_row->get<PetscGlobalIdx_mi_tag>().end()) {
        CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                          "Block <Block> not not fund");
      }
      if (schur_dof_c == schur_dofs_col->get<PetscGlobalIdx_mi_tag>().end()) {
        CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                          "Block <Schur> not not fund");
      }
#endif // NDEBUG
      insert(data_ptrs[2]->blockIndexPtr, *block_dof_r, *schur_dof_c, d);
    }

    if (vec_r_block[idx] != -1 && vec_c_block[idx] != -1) {
      auto block_dof_r =
          block_dofs_row->get<PetscGlobalIdx_mi_tag>().find(vec_r_block[idx]);
      auto block_dof_c =
          block_dofs_col->get<PetscGlobalIdx_mi_tag>().find(vec_c_block[idx]);
      insert(data_ptrs[3]->blockIndexPtr, *block_dof_r, *block_dof_c, d);
    }

    ++idx;
  }

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
      << data_ptrs[0]->blockIndexPtr->size();
  MOFEM_TAG_AND_LOG("SYNC", Sev::verbose, "NetsetSchur")
      << "(0, 1) " << schur_nb_local << " " << block_nb_local << " "
      << schur_nb_global << " " << block_nb_global << " "
      << data_ptrs[1]->blockIndexPtr->size();
  MOFEM_TAG_AND_LOG("SYNC", Sev::verbose, "NetsetSchur")
      << "(1, 0) " << block_nb_local << " " << schur_nb_local << " "
      << block_nb_global << " " << schur_nb_global << " "
      << data_ptrs[2]->blockIndexPtr->size();
  MOFEM_TAG_AND_LOG("SYNC", Sev::verbose, "NetsetSchur")
      << "(1, 1) " << block_nb_local << " " << block_nb_global << " "
      << data_ptrs[3]->blockIndexPtr->size();

  MOFEM_LOG_SEVERITY_SYNC(comm, Sev::verbose);

  return std::make_pair(mats_array, data_ptrs);
}

std::pair<SmartPetscObj<Mat>, SchurNestMatrixData>
createSchurNestedMatrix(std::pair<SmartPetscObj<DM>, SmartPetscObj<DM>> dms,
                        SchurNestMatrixData schur_net_data) {

  auto [schur_dm, block_dm] = dms;
  auto [mat_arrays, data_ptrs] = schur_net_data;

  auto schur_prb = getProblemPtr(schur_dm);
  auto block_prb = getProblemPtr(block_dm);

  auto schur_is = schur_prb->getSubData()->getSmartRowIs();
  auto block_is = block_prb->getSubData()->getSmartRowIs();

  std::array<IS, 2> is_a = {schur_is, block_is};
  std::array<Mat, 4> mats_a = {

      mat_arrays[0], mat_arrays[1], mat_arrays[2], mat_arrays[3]

  };

  MPI_Comm comm;
  CHKERR PetscObjectGetComm((PetscObject)schur_dm, &comm);

  Mat mat_raw;
  CHKERR MatCreateNest(

      comm, 2, is_a.data(), 2, is_a.data(), mats_a.data(), &mat_raw

  );

  return std::make_pair(SmartPetscObj<Mat>(mat_raw), schur_net_data);
}

} // namespace MoFEM