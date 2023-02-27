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

MoFEMErrorCode
SchurL2Mats::MatSetValues(Mat M, const EntitiesFieldData::EntData &row_data,
                          const EntitiesFieldData::EntData &col_data,
                          const MatrixDouble &mat, InsertMode iora) {
  MoFEMFunctionBegin;

  if (iora != ADD_VALUES) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
            "Other case of assembly not implemented");
  }

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

  const auto nb_rows = row_ind.size();
  const auto nb_cols = col_data.getIndices().size();

  const auto idx = SchurL2Mats::schurL2Storage.size();
  const auto size = SchurL2Mats::locMats.size();

  if (idx >= size) {
    SchurL2Mats::locMats.push_back(new MatrixDouble(nb_rows, nb_cols));
    SchurL2Mats::rowIndices.push_back(new VectorInt(nb_rows));
    SchurL2Mats::colIndices.push_back(new VectorInt(nb_cols));
  }

  // insert index
  auto get_field_name = [](auto &data) {
    return data.getFieldEntities()[0]->getName();
  };
  auto get_type = [](auto &data) {
    return data.getFieldEntities()[0]->getEntType();
  };

  auto get_uid = [](auto &data) {
    return data.getFieldEntities()[0]->getLocalUniqueId();
  };

  auto p = SchurL2Mats::schurL2Storage.emplace(idx, get_uid(row_data),
                                               get_uid(col_data));
#ifndef NDEBUG
  if (!p.second) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Is expected that data are inserted to storage");
  }
#endif // NDEBUG

  auto get_storage = [&p]() { return const_cast<SchurL2Mats &>(*p.first); };

  if (p.second) {
    get_storage().getMat().resize(nb_rows, nb_cols, false);
    noalias(get_storage().getMat()) = mat;

    auto add_indices = [](auto &storage, auto &ind) {
      storage.resize(ind.size(), false);
      noalias(storage) = ind;
    };

    const auto &col_ind = col_data.getIndices();
    add_indices(get_storage().getRowInd(), row_ind);
    add_indices(get_storage().getColInd(), col_ind);

  } else {

    switch (iora) {
    case ADD_VALUES:
      get_storage().getMat() += mat;
      break;
    case INSERT_VALUES:
      noalias(get_storage().getMat()) = mat;
      break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
              "Assembly type not implemented");
    };
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSchurAssembleBegin::doWork(int side, EntityType type,
                                            EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
  SchurL2Mats::schurL2Storage.clear();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSchurAssembleEnd::doWork(int side, EntityType type,
                                          EntitiesFieldData::EntData &data) {
  VectorInt ipiv;
  VectorDouble lapack_work;
  MoFEMFunctionBegin;

  auto invert_symm_mat = [&](MatrixDouble &m, auto &inv) {
    MoFEMFunctionBeginHot;
    const int nb = m.size1();
    inv.resize(nb, nb, false);
    inv.clear();
    for (int cc = 0; cc != nb; ++cc)
      inv(cc, cc) = -1;
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

  auto invert_nonsymm_mat = [&](MatrixDouble &m, auto &inv) {
    MoFEMFunctionBeginHot;
    const auto nb = m.size1();
    inv.resize(nb, nb, false);
    inv.clear();
    for (int c = 0; c != nb; ++c)
      inv(c, c) = -1;
    ipiv.resize(nb, false);
    const auto info = lapack_dgesv(nb, nb, &*m.data().begin(), nb,
                                   &*ipiv.begin(), &*inv.data().begin(), nb);
    if (info != 0)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
               "Can not invert matrix info = %d", info);
    MoFEMFunctionReturnHot(0);
  };

  auto assemble = [&](SmartPetscObj<Mat> M, auto &storage) {
    MoFEMFunctionBegin;
    if (M) {
      for (auto &s : storage) {
        auto &row_ind = s.getRowInd();
        auto &col_ind = s.getColInd();
        auto &m = s.getMat();
        CHKERR ::MatSetValues(M, row_ind.size(), &*row_ind.begin(),
                              col_ind.size(), &*col_ind.begin(),
                              &*m.data().begin(), ADD_VALUES);
      }
    }
    MoFEMFunctionReturn(0);
  };

  auto apply_schur = [&](auto &storage,

                         const auto ss,

                         auto lo_uid, auto hi_uid, auto row_it,
                         auto hi_row_it) {
    MoFEMFunctionBegin;

    auto add_off_mat = [&](auto row_uid, auto col_uid,
                           auto &offMatInvDiagOffMat) {
      MoFEMFunctionBegin;
      auto it = storage.template get<SchurL2Mats::uid_mi_tag>().find(
          boost::make_tuple(row_uid, col_uid));
      if (it == storage.template get<SchurL2Mats::uid_mi_tag>().end()) {
        const auto idx = SchurL2Mats::schurL2Storage.size();
        const auto size = SchurL2Mats::locMats.size();
        const auto nb_rows = offMatInvDiagOffMat.size1();
        const auto nb_cols = offMatInvDiagOffMat.size2();
        if (idx >= size) {
          SchurL2Mats::locMats.push_back(new MatrixDouble(nb_rows, nb_cols));
          SchurL2Mats::rowIndices.push_back(new VectorInt(nb_rows));
          SchurL2Mats::colIndices.push_back(new VectorInt(nb_cols));
        }
        auto p = SchurL2Mats::schurL2Storage.emplace(idx, row_uid, col_uid);
        auto &mat = p.first->getMat();
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

    for (; row_it != hi_row_it; ++row_it) {
      if (row_it->uidRow == row_it->uidCol) {

        if (getSymm()) {
          CHKERR invert_symm_mat(row_it->getMat(), invMat);
        } else {
          CHKERR invert_nonsymm_mat(row_it->getMat(), invMat);
        }

        const auto row_idx = row_it->iDX;

        auto c_lo =
            storage.template get<SchurL2Mats::col_mi_tag>().lower_bound(lo_uid);
        auto c_hi =
            storage.template get<SchurL2Mats::col_mi_tag>().lower_bound(hi_uid);

        for (; c_lo != c_hi; ++c_lo) {

          auto &row_uid = c_lo->uidRow;
          if (row_uid == row_it->uidRow)
            continue;

          auto &cc_off_mat = c_lo->getMat();
          invDiagOffMat.resize(cc_off_mat.size1(), invMat.size2(), false);
#ifndef NDEBUG
          if (invMat.size1() != cc_off_mat.size2()) {
            SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                     "Wrong size %d != %d", invMat.size1(),
                     cc_off_mat.size2());
          }
#endif // NDEBUG
          noalias(invDiagOffMat) = prod(cc_off_mat, invMat);

          auto r_lo =
              storage.template get<SchurL2Mats::row_mi_tag>().lower_bound(
                  lo_uid);
          auto r_hi =
              storage.template get<SchurL2Mats::row_mi_tag>().upper_bound(
                  hi_uid);

          for (; r_lo != r_hi; ++r_lo) {
            auto &col_uid = r_lo->uidCol;

            // Skip dagonal
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
            noalias(offMatInvDiagOffMat) =
                prod(invDiagOffMat, rr_off_mat);

            CHKERR add_off_mat(row_uid, col_uid, offMatInvDiagOffMat);

            if (symSchur[ss] && row_uid != col_uid) {
              transOffMatInvDiagOffMat.resize(
                  offMatInvDiagOffMat.size2(),
                  offMatInvDiagOffMat.size1(), false);
              noalias(transOffMatInvDiagOffMat) =
                  trans(offMatInvDiagOffMat);
              CHKERR add_off_mat(col_uid, row_uid,
                                 transOffMatInvDiagOffMat);
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
      CHKERR assemble(sequenceOfMats[ss], storage);
    }

    MoFEMFunctionReturn(0);
  };

  auto &storage = SchurL2Mats::schurL2Storage;
  CHKERR assemble(SmartPetscObj<Mat>(getKSPB(), true), storage);

  for (auto ss = 0; ss != fieldsName.size(); ++ss) {

    auto assemble = [&](auto &storage, auto ss, auto lo_uid, auto hi_uid) {
      MoFEMFunctionBegin;
      auto row_it =
          storage.template get<SchurL2Mats::row_mi_tag>().lower_bound(lo_uid);
      auto hi_row_it =
          storage.template get<SchurL2Mats::row_mi_tag>().upper_bound(hi_uid);
      CHKERR apply_schur(storage, ss, lo_uid, hi_uid, row_it, hi_row_it);
      CHKERR erase_factored(storage, ss, lo_uid, hi_uid);
      CHKERR assemble_off_diagonal(storage, ss);
      MoFEMFunctionReturn(0);
    };

    auto field_bit = getPtrFE()->mField.get_field_bit_number(fieldsName[ss]);
    auto row_ents = fieldEnts[ss];

    if (row_ents) {
      auto end = storage.template get<SchurL2Mats::row_mi_tag>().end();
      for (auto p = row_ents->pair_begin(); p != row_ents->pair_end(); ++p) {
        auto lo_uid =
            FieldEntity::getLoLocalEntityBitNumber(field_bit, p->first);
        auto hi_uid =
            FieldEntity::getHiLocalEntityBitNumber(field_bit, p->second);
        CHKERR assemble(storage, ss, lo_uid, hi_uid);
      }
    } else {
      auto lo_uid =
          FieldEntity::getLoLocalEntityBitNumberByType(field_bit, MBVERTEX);
      auto hi_uid =
          FieldEntity::getHiLocalEntityBitNumberByType(field_bit, MBENTITYSET);
      CHKERR assemble(storage, ss, lo_uid, hi_uid);
    }
  }

  MoFEMFunctionReturn(0);
}
} // namespace MoFEM