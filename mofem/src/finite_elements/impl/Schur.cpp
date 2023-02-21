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

SchurL2Mats::SchurL2Mats(const size_t idx, const UId uid_row,
                         const EntityType row_type, const std::string row_field,
                         const UId uid_col, const EntityType col_type,
                         const std::string col_field)
    : iDX(idx), uidRow(uid_row), rowType(row_type), rowField(row_field),
      uidCol(uid_col), colType(col_type), colField(col_field) {}

MoFEMErrorCode
SchurL2Mats::MatSetValues(Mat M, const EntitiesFieldData::EntData &row_data,
                          const EntitiesFieldData::EntData &col_data,
                          const MatrixDouble &mat, InsertMode iora) {
  MoFEMFunctionBegin;

  if(iora != ADD_VALUES) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
            "Other case of assembly not implemented");
  }

  const auto nb_rows = row_data.getIndices().size();
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

  auto p = SchurL2Mats::schurL2Storage.emplace(
      idx, get_uid(row_data), get_type(row_data), get_field_name(row_data),
      get_uid(col_data), get_type(col_data), get_field_name(col_data));
#ifndef NDEBUG
  if (!p.second) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Is expected that data are inserted to storage");
  }
#endif // NDEBUG

  // set data to storage
  auto get_storage = [&p]() { return const_cast<SchurL2Mats &>(*p.first); };
  auto add_indices = [](auto &storage, auto &ind) {
    storage.resize(ind.size(), false);
    noalias(storage) = ind;
  };

  get_storage().getMat().resize(nb_rows, nb_cols, false);
  noalias(get_storage().getMat()) = mat;

  const auto &row_ind = row_data.getIndices();
  const auto &col_ind = col_data.getIndices();
  add_indices(get_storage().getRowInd(), row_ind);
  add_indices(get_storage().getColInd(), col_ind);

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

  MatrixDouble inv_mat;
  MatrixDouble inv_diag_off_mat;
  MatrixDouble off_mat_inv_diag_off_mat;

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

  auto apply_schur = [&](auto &storage, const auto ss) {
    MoFEMFunctionBegin;

    const auto row_field = fieldsName[ss];
    const auto row_type = fieldEntTypes[ss];

    auto row_it = storage.template get<SchurL2Mats::row_mi_tag>().find(
        boost::make_tuple(row_type, row_field));
    if (row_it != storage.template get<SchurL2Mats::row_mi_tag>().end()) {

      if (getSymm()) {
        CHKERR invert_symm_mat(row_it->getMat(), inv_mat);
      } else {
        CHKERR invert_nonsymm_mat(row_it->getMat(), inv_mat);
      }

      const auto row_idx = row_it->iDX;
      auto c_lo = storage.template get<SchurL2Mats::col_mi_tag>().lower_bound(
          boost::make_tuple(row_type, row_field));
      auto c_hi = storage.template get<SchurL2Mats::col_mi_tag>().upper_bound(
          boost::make_tuple(row_type, row_field));

      for (; c_lo != c_hi; ++c_lo) {
        if (row_idx != c_lo->iDX) {

          auto &col_uid = c_lo->uidCol;
          auto &cc_off_mat = c_lo->getMat();
          inv_diag_off_mat.resize(inv_mat.size1(), cc_off_mat.size2(), false);
          noalias(inv_diag_off_mat) = prod(inv_mat, cc_off_mat);

          auto r_lo =
              storage.template get<SchurL2Mats::row_mi_tag>().lower_bound(
                  boost::make_tuple(row_type, row_field));
          auto r_hi =
              storage.template get<SchurL2Mats::row_mi_tag>().upper_bound(
                  boost::make_tuple(row_type, row_field));
          for (; r_lo != r_hi; ++r_hi) {
            if (row_idx != r_lo->iDX) {
              auto &row_uid = r_lo->uidRow;
              auto &rr_off_mat = r_lo->getMat();
              off_mat_inv_diag_off_mat.resize(rr_off_mat.size1(),
                                              inv_diag_off_mat.size2(), false);
              noalias(off_mat_inv_diag_off_mat) =
                  prod(rr_off_mat, inv_diag_off_mat);

              auto it = storage.template get<SchurL2Mats::uid_mi_tag>().find(
                  boost::make_tuple(row_uid, col_uid));
              if (it == storage.template get<SchurL2Mats::uid_mi_tag>().end()) {
                const auto idx = SchurL2Mats::schurL2Storage.size();
                const auto size = SchurL2Mats::locMats.size();
                const auto nb_rows = off_mat_inv_diag_off_mat.size1();
                const auto nb_cols = off_mat_inv_diag_off_mat.size2();
                if (idx > size) {
                  SchurL2Mats::locMats.push_back(
                      new MatrixDouble(nb_rows, nb_cols));
                  SchurL2Mats::rowIndices.push_back(new VectorInt(nb_rows));
                  SchurL2Mats::colIndices.push_back(new VectorInt(nb_cols));
                }
                auto p = SchurL2Mats::schurL2Storage.emplace(
                    idx, row_uid, r_lo->rowType, r_lo->rowField, col_uid,
                    c_lo->colType, c_lo->colField);
                auto &mat = p.first->getMat();
                mat.swap(off_mat_inv_diag_off_mat);
              } else {
                auto &mat = it->getMat();
                mat -= off_mat_inv_diag_off_mat;
              }
            }
          }
        }
      }

    } else {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Diag not found");
    }

    MoFEMFunctionReturn(0);
  };

  auto assemble_diagonal = [&](auto &storage, const auto ss) {
    MoFEMFunctionBegin;

    auto row_type = fieldEntTypes[ss];
    auto &row_field = fieldsName[ss];

    auto row_it = storage.template get<SchurL2Mats::row_mi_tag>().find(
        boost::make_tuple(row_type, row_field));
    if (row_it != storage.template get<SchurL2Mats::row_mi_tag>().end()) {
      auto &ind_row = row_it->getRowInd();
      auto &ind_col = row_it->getColInd();
      if (auto ao = sequenceOfAOs[ss]) {
        CHKERR AOApplicationToPetsc(ao, ind_row.size(), &*ind_row.begin());
        CHKERR AOApplicationToPetsc(ao, ind_col.size(), &*ind_col.begin());
      }
      auto &m = row_it->getMat();
      CHKERR ::MatSetValues(sequenceOfMats[ss], ind_row.size(),
                            &*ind_row.begin(), ind_col.size(),
                            &*ind_col.begin(), &*m.data().begin(), ADD_VALUES);

    } else {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Diag not found");
    }

    MoFEMFunctionReturn(0);
  };

  auto erase_factored = [&](auto &storage, const auto ss) {
    MoFEMFunctionBegin;

    auto row_type = fieldEntTypes[ss];
    auto &row_field = fieldsName[ss];

    auto r_lo = storage.template get<SchurL2Mats::row_mi_tag>().lower_bound(
        boost::make_tuple(row_type, row_field));
    auto r_hi = storage.template get<SchurL2Mats::row_mi_tag>().upper_bound(
        boost::make_tuple(row_type, row_field));
    storage.template get<SchurL2Mats::row_mi_tag>().erase(r_lo, r_hi);

    auto c_lo = storage.template get<SchurL2Mats::col_mi_tag>().lower_bound(
        boost::make_tuple(row_type, row_field));
    auto c_hi = storage.template get<SchurL2Mats::col_mi_tag>().upper_bound(
        boost::make_tuple(row_type, row_field));
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
    CHKERR assemble(sequenceOfMats[ss], storage);

    MoFEMFunctionReturn(0);
  };

  auto &storage = SchurL2Mats::schurL2Storage;
  CHKERR assemble(SmartPetscObj<Mat>(getKSPB(), true), storage);

  for (auto ss = 0; ss != fieldsName.size(); ++ss) {
    CHKERR apply_schur(storage, ss);
    CHKERR assemble_diagonal(storage, ss);
    CHKERR erase_factored(storage, ss);
    CHKERR assemble_off_diagonal(storage, ss);
  }

  MoFEMFunctionReturn(0);
}
} // namespace MoFEM