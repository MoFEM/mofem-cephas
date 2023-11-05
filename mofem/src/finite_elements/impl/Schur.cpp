/**
 * @file Schur.cpp
 * @brief Implementation of Schur Complement
 * @date 2023-02-03
 *
 * @copyright Copyright (c) 2023
 *
 */

namespace MoFEM {

SchurL2Mats::SchurL2Storage SchurL2Mats::schurL2Storage;
boost::ptr_vector<MatrixDouble> SchurL2Mats::locMats;
boost::ptr_vector<VectorInt> SchurL2Mats::rowIndices;
boost::ptr_vector<VectorInt> SchurL2Mats::colIndices;

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

  CHKERR matSetValuesPtr(M, row_data, col_data, mat, iora);

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

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSchurAssembleBegin::doWork(int side, EntityType type,
                                            EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
#ifndef NDEBUG
  MOFEM_LOG("SELF", Sev::noisy) << "Schur assemble begin";
#endif
  SchurL2Mats::schurL2Storage.clear();

  // auto zero_mats = [&](auto &mats) {
  //   for (auto &m : mats) {
  //     m.resize(0, 0, false);
  //   }
  // };

  // zero_mats(SchurL2Mats::locMats);

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

} // namespace MoFEM