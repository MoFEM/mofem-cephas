/** file DataOperators.cpp

  \brief implementation of Data Operators for Forces and Sources

*/

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#ifdef __cplusplus
extern "C" {
#endif
#include <cblas.h>
#include <lapack_wrap.h>
#include <gm_rule.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

DataOperator::DataOperator(const bool symm)
    :

      sYmm(symm),

      doEntities{true, true, true, true, true, true,
                 true, true, true, true, true, true},

      doVertices(doEntities[MBVERTEX]), doEdges(doEntities[MBEDGE]),
      doQuads(doEntities[MBQUAD]), doTris(doEntities[MBTRI]),
      doTets(doEntities[MBTET]), doPrisms(doEntities[MBPRISM]) {

  /// This not yet implemented, switch off.
  doEntities[MBPOLYGON] = false;
  doEntities[MBPYRAMID] = false;
  doEntities[MBKNIFE] = false;
  doEntities[MBPOLYHEDRON] = false;
}

template <bool Symm>
MoFEMErrorCode DataOperator::opLhs(DataForcesAndSourcesCore &row_data,
                                   DataForcesAndSourcesCore &col_data) {
  MoFEMFunctionBeginHot;

  auto do_col_entity =
      [&](boost::ptr_vector<DataForcesAndSourcesCore::EntData> &row_ent_data,
          const int ss, const EntityType row_type, const EntityType low_type,
          const EntityType hi_type) {
        MoFEMFunctionBegin;
        for (EntityType col_type = low_type; col_type != hi_type; ++col_type) {
          auto &col_ent_data = col_data.dataOnEntities[col_type];
          for (size_t SS = 0; SS != col_ent_data.size(); SS++) {
            if (col_ent_data[SS].getFieldData().size())
              CHKERR doWork(ss, SS, row_type, col_type, row_ent_data[ss],
                            col_ent_data[SS]);
          }
        }
        MoFEMFunctionReturn(0);
      };

  auto do_row_entity = [&](const EntityType type) {
    MoFEMFunctionBegin;
    auto &row_ent_data = row_data.dataOnEntities[type];
    for (size_t ss = 0; ss != row_ent_data.size(); ++ss) {
      size_t SS = 0;
      if (Symm)
        SS = ss;
      for (; SS < col_data.dataOnEntities[type].size(); SS++) {
        CHKERR doWork(ss, SS, type, type, row_ent_data[ss],
                      col_data.dataOnEntities[type][SS]);
      }
      if (!Symm)
        CHKERR do_col_entity(row_ent_data, ss, type, MBVERTEX, type);
      CHKERR do_col_entity(row_ent_data, ss, type,
                           static_cast<EntityType>(type + 1), MBMAXTYPE);
    }
    MoFEMFunctionReturn(0);
  };

  for (EntityType row_type = MBVERTEX; row_type != MBENTITYSET; ++row_type) {
    if (doEntities[row_type]) {
      CHKERR do_row_entity(row_type);
    }
  }

  if (doEntities[MBENTITYSET]) {
    for (unsigned int mm = 0; mm != row_data.dataOnEntities[MBENTITYSET].size();
         ++mm) {
      if (!row_data.dataOnEntities[MBENTITYSET][mm].getFieldData().empty()) {
        CHKERR do_row_entity(MBENTITYSET);
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DataOperator::opLhs(DataForcesAndSourcesCore &row_data,
                                   DataForcesAndSourcesCore &col_data) {
  if (getSymm())
    return opLhs<true>(row_data, col_data);
  else
    return opLhs<false>(row_data, col_data);
}

template <bool ErrorIfNoBase>
MoFEMErrorCode
DataOperator::opRhs(DataForcesAndSourcesCore &data,
                    const std::array<bool, MBMAXTYPE> &do_entities) {
  MoFEMFunctionBegin;

  auto do_entity = [&](auto type) {
    MoFEMFunctionBegin;

    auto &ent_data = data.dataOnEntities[type];
    const size_t size = ent_data.size();
    for (size_t ss = 0; ss != size; ++ss) {

      auto &side_data = ent_data[ss];

      if (ErrorIfNoBase) {
        if (side_data.getFieldData().size() &&
            (side_data.getBase() == NOBASE ||
             side_data.getBase() == LASTBASE)) {
          for (VectorDofs::iterator it = side_data.getFieldDofs().begin();
               it != side_data.getFieldDofs().end(); it++)
            if ((*it) && (*it)->getActive())
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "No base on");
        }
      }

      CHKERR doWork(ss, type, side_data);
    }

    MoFEMFunctionReturn(0);
  };

  for (EntityType row_type = MBVERTEX; row_type != MBENTITYSET; ++row_type) {
    if (do_entities[row_type]) {
      CHKERR do_entity(row_type);
    }
  }

  if (do_entities[MBENTITYSET]) {
    // This is odd behaviour, diffrent than for other entities. Should be
    // changed that behaviour is consistent,
    for (unsigned int mm = 0; mm != data.dataOnEntities[MBENTITYSET].size();
         ++mm) {
      if (!data.dataOnEntities[MBENTITYSET][mm].getFieldData().empty()) {
        CHKERR doWork(mm, MBENTITYSET, data.dataOnEntities[MBENTITYSET][mm]);
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DataOperator::opRhs(DataForcesAndSourcesCore &data,
                                   const bool error_if_no_base) {
  if (error_if_no_base)
    return opRhs<true>(data, doEntities);
  else
    return opRhs<false>(data, doEntities);
}

template <>
MoFEMErrorCode invertTensor3by3<3, double, ublas::row_major, DoubleAllocator>(
    MatrixDouble &jac_data, VectorDouble &det_data,
    MatrixDouble &inv_jac_data) {
  MoFEMFunctionBegin;
  auto A = getFTensor2FromMat<3, 3>(jac_data);
  int nb_gauss_pts = jac_data.size2();
  det_data.resize(nb_gauss_pts, false);
  inv_jac_data.resize(3, nb_gauss_pts, false);
  auto det = getFTensor0FromVec(det_data);
  auto I = getFTensor2FromMat<3, 3>(inv_jac_data);
  for (int gg = 0; gg != nb_gauss_pts; ++gg) {
    CHKERR determinantTensor3by3(A, det);
    CHKERR invertTensor3by3(A, det, I);
    ++A;
    ++det;
    ++I;
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSetInvJacH1::doWork(int side, EntityType type,
                                     DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  auto transform_base = [&](MatrixDouble &diff_n) {
    MoFEMFunctionBeginHot;

    if (!diff_n.size1())
      MoFEMFunctionReturnHot(0);
    if (!diff_n.size2())
      MoFEMFunctionReturnHot(0);

    const int nb_base_functions = diff_n.size2() / 3;
    const int nb_gauss_pts = diff_n.size1();
    diffNinvJac.resize(diff_n.size1(), diff_n.size2(), false);

    double *t_diff_n_ptr = &*diff_n.data().begin();
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_diff_n(
        t_diff_n_ptr, &t_diff_n_ptr[1], &t_diff_n_ptr[2]);
    double *t_inv_n_ptr = &*diffNinvJac.data().begin();
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_inv_diff_n(
        t_inv_n_ptr, &t_inv_n_ptr[1], &t_inv_n_ptr[2]);

    for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
      for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
        t_inv_diff_n(i) = t_diff_n(j) * tInvJac(j, i);
        ++t_diff_n;
        ++t_inv_diff_n;
      }
    }

    diff_n.swap(diffNinvJac);

    MoFEMFunctionReturnHot(0);
  };

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
    const auto base = static_cast<FieldApproximationBase>(b);
    CHKERR transform_base(data.getDiffN(base));
  }

  switch (type) {
  case MBVERTEX:
    for (auto &m : data.getBBDiffNMap())
      if (m.second)
        CHKERR transform_base(*(m.second));
    break;
  default:
    for (auto &ptr : data.getBBDiffNByOrderArray())
      if (ptr)
        CHKERR transform_base(*ptr);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetInvJacHdivAndHcurl::doWork(int side, EntityType type,
                                DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (type == MBVERTEX)
    MoFEMFunctionReturnHot(0);

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    const unsigned int nb_gauss_pts = data.getDiffN(base).size1();
    const unsigned int nb_base_functions = data.getDiffN(base).size2() / 9;
    if (!nb_base_functions)
      continue;

    diffHdivInvJac.resize(nb_gauss_pts, data.getDiffN(base).size2(), false);

    auto t_diff_n = data.getFTensor2DiffN<3, 3>(base);
    double *t_inv_diff_n_ptr = &*diffHdivInvJac.data().begin();
    FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_inv_diff_n(
        t_inv_diff_n_ptr, &t_inv_diff_n_ptr[HVEC0_1],
        &t_inv_diff_n_ptr[HVEC0_2],

        &t_inv_diff_n_ptr[HVEC1_0], &t_inv_diff_n_ptr[HVEC1_1],
        &t_inv_diff_n_ptr[HVEC1_2],

        &t_inv_diff_n_ptr[HVEC2_0], &t_inv_diff_n_ptr[HVEC2_1],
        &t_inv_diff_n_ptr[HVEC2_2]);

    for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
      for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
        t_inv_diff_n(k, i) = t_diff_n(k, j) * tInvJac(j, i);
        ++t_diff_n;
        ++t_inv_diff_n;
      }
    }

    data.getDiffN(base).swap(diffHdivInvJac);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSetContravariantPiolaTransform::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (CN::Dimension(type) > 1) {

    for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {

      FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

      const unsigned int nb_base_functions = data.getN(base).size2() / 3;
      if (!nb_base_functions)
        continue;

      const unsigned int nb_gauss_pts = data.getN(base).size1();
      double const a = 1. / vOlume;

      piolaN.resize(nb_gauss_pts, data.getN(base).size2(), false);
      if (data.getN(base).size2() > 0) {
        auto t_n = data.getFTensor1N<3>(base);
        double *t_transformed_n_ptr = &*piolaN.data().begin();
        FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_transformed_n(
            t_transformed_n_ptr, // HVEC0
            &t_transformed_n_ptr[HVEC1], &t_transformed_n_ptr[HVEC2]);
        for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
          for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
            t_transformed_n(i) = a * (tJac(i, k) * t_n(k));
            ++t_n;
            ++t_transformed_n;
          }
        }
        data.getN(base).swap(piolaN);
      }

      piolaDiffN.resize(nb_gauss_pts, data.getDiffN(base).size2(), false);
      if (data.getDiffN(base).size2() > 0) {
        auto t_diff_n = data.getFTensor2DiffN<3, 3>(base);
        double *t_transformed_diff_n_ptr = &*piolaDiffN.data().begin();
        FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>
            t_transformed_diff_n(t_transformed_diff_n_ptr,
                                 &t_transformed_diff_n_ptr[HVEC0_1],
                                 &t_transformed_diff_n_ptr[HVEC0_2],
                                 &t_transformed_diff_n_ptr[HVEC1_0],
                                 &t_transformed_diff_n_ptr[HVEC1_1],
                                 &t_transformed_diff_n_ptr[HVEC1_2],
                                 &t_transformed_diff_n_ptr[HVEC2_0],
                                 &t_transformed_diff_n_ptr[HVEC2_1],
                                 &t_transformed_diff_n_ptr[HVEC2_2]);
        for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
          for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
            t_transformed_diff_n(i, k) = a * tJac(i, j) * t_diff_n(j, k);
            ++t_diff_n;
            ++t_transformed_diff_n;
          }
        }
        data.getDiffN(base).swap(piolaDiffN);
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetCovariantPiolaTransform::doWork(int side, EntityType type,
                                     DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (type == MBVERTEX)
    MoFEMFunctionReturnHot(0);

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    const unsigned int nb_base_functions = data.getN(base).size2() / 3;
    if (!nb_base_functions)
      continue;

    const unsigned int nb_gauss_pts = data.getN(base).size1();
    piolaN.resize(nb_gauss_pts, data.getN(base).size2(), false);
    piolaDiffN.resize(nb_gauss_pts, data.getDiffN(base).size2(), false);

    auto t_n = data.getFTensor1N<3>(base);
    double *t_transformed_n_ptr = &*piolaN.data().begin();
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_transformed_n(
        t_transformed_n_ptr, &t_transformed_n_ptr[HVEC1],
        &t_transformed_n_ptr[HVEC2]);
    auto t_diff_n = data.getFTensor2DiffN<3, 3>(base);
    double *t_transformed_diff_n_ptr = &*piolaDiffN.data().begin();
    FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_transformed_diff_n(
        t_transformed_diff_n_ptr, &t_transformed_diff_n_ptr[HVEC0_1],
        &t_transformed_diff_n_ptr[HVEC0_2], &t_transformed_diff_n_ptr[HVEC1_0],
        &t_transformed_diff_n_ptr[HVEC1_1], &t_transformed_diff_n_ptr[HVEC1_2],
        &t_transformed_diff_n_ptr[HVEC2_0], &t_transformed_diff_n_ptr[HVEC2_1],
        &t_transformed_diff_n_ptr[HVEC2_2]);

    for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
      for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
        t_transformed_n(i) = tInvJac(k, i) * t_n(k);
        t_transformed_diff_n(i, k) = tInvJac(j, i) * t_diff_n(j, k);
        ++t_n;
        ++t_transformed_n;
        ++t_diff_n;
        ++t_transformed_diff_n;
      }
    }
    data.getN(base).swap(piolaN);
    data.getDiffN(base).swap(piolaDiffN);
  }

  // data.getBase() = base;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpGetCoordsAndNormalsOnPrism::doWork(int side, EntityType type,
                                     DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (data.getFieldData().size() == 0)
    MoFEMFunctionReturnHot(0);
  const int valid_edges3[] = {1, 1, 1, 0, 0, 0, 0, 0, 0};
  const int valid_faces3[] = {0, 0, 0, 1, 0, 0, 0, 0, 0};
  const int valid_edges4[] = {0, 0, 0, 0, 0, 0, 1, 1, 1};
  const int valid_faces4[] = {0, 0, 0, 0, 1, 0, 0, 0, 0};

  if (type == MBEDGE) {
    if (!valid_edges3[side] || valid_edges4[side])
      MoFEMFunctionReturnHot(0);
  } else if (type == MBTRI) {
    if (!valid_faces3[side] || valid_faces4[side])
      MoFEMFunctionReturnHot(0);
  }

  switch (type) {
  case MBVERTEX: {
    for (unsigned int gg = 0; gg < data.getN().size1(); ++gg) {
      for (int dd = 0; dd < 3; dd++) {
        cOords_at_GaussPtF3(gg, dd) =
            cblas_ddot(3, &data.getN(gg)[0], 1, &data.getFieldData()[dd], 3);
        tAngent1_at_GaussPtF3(gg, dd) = cblas_ddot(
            3, &data.getDiffN()(gg, 0), 2, &data.getFieldData()[dd], 3);
        tAngent2_at_GaussPtF3(gg, dd) = cblas_ddot(
            3, &data.getDiffN()(gg, 1), 2, &data.getFieldData()[dd], 3);
        cOords_at_GaussPtF4(gg, dd) = cblas_ddot(
            3, &data.getN(gg)[0], 1, &data.getFieldData()[9 + dd], 3);
        tAngent1_at_GaussPtF4(gg, dd) = cblas_ddot(
            3, &data.getDiffN()(gg, 6 + 0), 2, &data.getFieldData()[9 + dd], 3);
        tAngent2_at_GaussPtF4(gg, dd) = cblas_ddot(
            3, &data.getDiffN()(gg, 6 + 1), 2, &data.getFieldData()[9 + dd], 3);
      }
    }
  } break;
  case MBEDGE:
  case MBTRI: {
    if (2 * data.getN().size2() != data.getDiffN().size2()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    unsigned int nb_dofs = data.getFieldData().size();
    if (nb_dofs % 3 != 0) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    if (nb_dofs > 3 * data.getN().size2()) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "data inconsistency, side %d type %d", side, type);
    }
    for (unsigned int gg = 0; gg < data.getN().size1(); ++gg) {
      for (int dd = 0; dd < 3; dd++) {
        if ((type == MBTRI && valid_faces3[side]) ||
            (type == MBEDGE && valid_edges3[side])) {
          cOords_at_GaussPtF3(gg, dd) += cblas_ddot(
              nb_dofs / 3, &data.getN(gg)[0], 1, &data.getFieldData()[dd], 3);
          tAngent1_at_GaussPtF3(gg, dd) +=
              cblas_ddot(nb_dofs / 3, &data.getDiffN()(gg, 0), 2,
                         &data.getFieldData()[dd], 3);
          tAngent2_at_GaussPtF3(gg, dd) +=
              cblas_ddot(nb_dofs / 3, &data.getDiffN()(gg, 1), 2,
                         &data.getFieldData()[dd], 3);
        } else if ((type == MBTRI && valid_faces4[side]) ||
                   (type == MBEDGE && valid_edges4[side])) {
          cOords_at_GaussPtF4(gg, dd) += cblas_ddot(
              nb_dofs / 3, &data.getN(gg)[0], 1, &data.getFieldData()[dd], 3);
          tAngent1_at_GaussPtF4(gg, dd) +=
              cblas_ddot(nb_dofs / 3, &data.getDiffN()(gg, 0), 2,
                         &data.getFieldData()[dd], 3);
          tAngent2_at_GaussPtF4(gg, dd) +=
              cblas_ddot(nb_dofs / 3, &data.getDiffN()(gg, 1), 2,
                         &data.getFieldData()[dd], 3);
        }
      }
    }
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpGetCoordsAndNormalsOnPrism::calculateNormals() {
  MoFEMFunctionBegin;

  sPin.resize(3, 3);
  sPin.clear();
  nOrmals_at_GaussPtF3.resize(tAngent1_at_GaussPtF3.size1(), 3, false);
  for (unsigned int gg = 0; gg < tAngent1_at_GaussPtF3.size1(); ++gg) {
    ierr = Spin(&*sPin.data().begin(), &tAngent1_at_GaussPtF3(gg, 0));
    CHKERRG(ierr);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1., &*sPin.data().begin(), 3,
                &tAngent2_at_GaussPtF3(gg, 0), 1, 0.,
                &nOrmals_at_GaussPtF3(gg, 0), 1);
  }
  sPin.clear();
  nOrmals_at_GaussPtF4.resize(tAngent1_at_GaussPtF4.size1(), 3, false);
  for (unsigned int gg = 0; gg < tAngent1_at_GaussPtF4.size1(); ++gg) {
    ierr = Spin(&*sPin.data().begin(), &tAngent1_at_GaussPtF4(gg, 0));
    CHKERRG(ierr);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1., &*sPin.data().begin(), 3,
                &tAngent2_at_GaussPtF4(gg, 0), 1, 0.,
                &nOrmals_at_GaussPtF4(gg, 0), 1);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSetContravariantPiolaTransformOnFace::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  FTensor::Index<'i', 3> i;
  MoFEMFunctionBegin;

  if (moab::CN::Dimension(type) != 2)
    MoFEMFunctionReturnHot(0);

  if (normalRawPtr == nullptr && normalsAtGaussPtsRawPtr == nullptr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
            "Pointer to normal/normals not set");

  bool normal_is_at_gauss_pts = (normalsAtGaussPtsRawPtr != nullptr);
  if (normal_is_at_gauss_pts)
    normal_is_at_gauss_pts = (normalsAtGaussPtsRawPtr->size1() != 0);

  auto apply_transform_linear_geometry = [&](auto base, auto nb_gauss_pts,
                                             auto nb_base_functions) {
    MoFEMFunctionBegin;
    const auto &normal = *normalRawPtr;
    auto t_normal = FTensor::Tensor1<double, 3>{normal[normalShift + 0],
                                                normal[normalShift + 1],
                                                normal[normalShift + 2]};
    const auto l02 = t_normal(i) * t_normal(i);
    auto t_base = data.getFTensor1N<3>(base);
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      for (int bb = 0; bb != nb_base_functions; ++bb) {
        const auto v = t_base(0);
        t_base(i) = (v / l02) * t_normal(i);
        ++t_base;
      }
    }
    MoFEMFunctionReturn(0);
  };

  auto apply_transform_nonlinear_geometry = [&](auto base, auto nb_gauss_pts,
                                                auto nb_base_functions) {
    MoFEMFunctionBegin;
    const MatrixDouble &normals_at_pts = *normalsAtGaussPtsRawPtr;
    auto t_normal = FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3>(
        &normals_at_pts(0, 0), &normals_at_pts(0, 1), &normals_at_pts(0, 2));

    auto t_base = data.getFTensor1N<3>(base);
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      const auto l2 = t_normal(i) * t_normal(i);
      for (int bb = 0; bb != nb_base_functions; ++bb) {
        const auto v = t_base(0);
        t_base(i) = (v / l2) * t_normal(i);
        ++t_base;
      }
      ++t_normal;
    }
    MoFEMFunctionReturn(0);
  };

  if (normal_is_at_gauss_pts) {
    for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {

      FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
      const auto &base_functions = data.getN(base);
      const auto nb_gauss_pts = base_functions.size1();

      if (nb_gauss_pts) {

        if (normalsAtGaussPtsRawPtr->size1() != nb_gauss_pts)
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "normalsAtGaussPtsRawPtr has inconsistent number of "
                  "integration "
                  "points");

        const auto nb_base_functions = base_functions.size2() / 3;
        CHKERR apply_transform_nonlinear_geometry(base, nb_gauss_pts,
                                                  nb_base_functions);
      }
    }
  } else {
    for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {

      FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
      const auto &base_functions = data.getN(base);
      const auto nb_gauss_pts = base_functions.size1();

      if (nb_gauss_pts) {
        const auto nb_base_functions = base_functions.size2() / 3;
        CHKERR apply_transform_linear_geometry(base, nb_gauss_pts,
                                               nb_base_functions);
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSetCovariantPiolaTransformOnFace::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  const auto type_dim = moab::CN::Dimension(type);
  if (type_dim != 1 && type_dim != 2)
    MoFEMFunctionReturnHot(0);

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 2> k;

  FTensor::Tensor2<FTensor::PackPtr<const double *, 3>, 3, 3> t_m(
      &tAngent0[0], &tAngent1[0], &nOrmal[0],

      &tAngent0[1], &tAngent1[1], &nOrmal[1],

      &tAngent0[2], &tAngent1[2], &nOrmal[2]);
  double det;
  FTensor::Tensor2<double, 3, 3> t_inv_m;
  CHKERR determinantTensor3by3(t_m, det);
  CHKERR invertTensor3by3(t_m, det, t_inv_m);

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; ++b) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    auto &baseN = data.getN(base);
    auto &diffBaseN = data.getDiffN(base);

    int nb_dofs = baseN.size2() / 3;
    int nb_gauss_pts = baseN.size1();

    MatrixDouble piola_n(baseN.size1(), baseN.size2());
    MatrixDouble diff_piola_n(diffBaseN.size1(), diffBaseN.size2());

    if (nb_dofs > 0 && nb_gauss_pts > 0) {

      FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_h_curl(
          &baseN(0, HVEC0), &baseN(0, HVEC1), &baseN(0, HVEC2));
      FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2> t_diff_h_curl(
          &diffBaseN(0, HVEC0_0), &diffBaseN(0, HVEC0_1),
          &diffBaseN(0, HVEC1_0), &diffBaseN(0, HVEC1_1),
          &diffBaseN(0, HVEC2_0), &diffBaseN(0, HVEC2_1));
      FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_transformed_h_curl(
          &piola_n(0, HVEC0), &piola_n(0, HVEC1), &piola_n(0, HVEC2));
      FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>
          t_transformed_diff_h_curl(
              &diff_piola_n(0, HVEC0_0), &diff_piola_n(0, HVEC0_1),
              &diff_piola_n(0, HVEC1_0), &diff_piola_n(0, HVEC1_1),
              &diff_piola_n(0, HVEC2_0), &diff_piola_n(0, HVEC2_1));

      int cc = 0;
      if (normalsAtGaussPts.size1() == (unsigned int)nb_gauss_pts) {
        // HO geometry is set, so jacobian is different at each gauss point
        FTensor::Tensor2<FTensor::PackPtr<const double *, 3>, 3, 3> t_m_at_pts(
            &tangent0AtGaussPt(0, 0), &tangent1AtGaussPt(0, 0),
            &normalsAtGaussPts(0, 0), &tangent0AtGaussPt(0, 1),
            &tangent1AtGaussPt(0, 1), &normalsAtGaussPts(0, 1),
            &tangent0AtGaussPt(0, 2), &tangent1AtGaussPt(0, 2),
            &normalsAtGaussPts(0, 2));
        for (int gg = 0; gg < nb_gauss_pts; ++gg) {
          CHKERR determinantTensor3by3(t_m_at_pts, det);
          CHKERR invertTensor3by3(t_m_at_pts, det, t_inv_m);
          for (int ll = 0; ll != nb_dofs; ll++) {
            t_transformed_h_curl(i) = t_inv_m(j, i) * t_h_curl(j);
            t_transformed_diff_h_curl(i, k) =
                t_inv_m(j, i) * t_diff_h_curl(j, k);
            ++t_h_curl;
            ++t_transformed_h_curl;
            ++t_diff_h_curl;
            ++t_transformed_diff_h_curl;
            ++cc;
          }
          ++t_m_at_pts;
        }
      } else {
        for (int gg = 0; gg < nb_gauss_pts; ++gg) {
          for (int ll = 0; ll != nb_dofs; ll++) {
            t_transformed_h_curl(i) = t_inv_m(j, i) * t_h_curl(j);
            t_transformed_diff_h_curl(i, k) =
                t_inv_m(j, i) * t_diff_h_curl(j, k);
            ++t_h_curl;
            ++t_transformed_h_curl;
            ++t_diff_h_curl;
            ++t_transformed_diff_h_curl;
            ++cc;
          }
        }
      }
      if (cc != nb_gauss_pts * nb_dofs)
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "Data inconsistency");

      baseN.swap(piola_n);
      diffBaseN.swap(diff_piola_n);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpGetHOTangentOnEdge::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  int nb_dofs = data.getFieldData().size();
  if (nb_dofs == 0)
    MoFEMFunctionReturnHot(0);

  int nb_gauss_pts = data.getN().size1();
  tAngent.resize(nb_gauss_pts, 3, false);

  int nb_approx_fun = data.getN().size2();
  double *diff = &*data.getDiffN().data().begin();
  double *dofs[] = {&data.getFieldData()[0], &data.getFieldData()[1],
                    &data.getFieldData()[2]};

  tAngent.resize(nb_gauss_pts, 3, false);

  switch (type) {
  case MBVERTEX:
    for (int dd = 0; dd != 3; dd++) {
      for (int gg = 0; gg != nb_gauss_pts; ++gg) {
        tAngent(gg, dd) = cblas_ddot(2, diff, 1, dofs[dd], 3);
      }
    }
    break;
  case MBEDGE:
    if (nb_dofs % 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
              "Approximated field should be rank 3, i.e. vector in 3d space");
    }
    for (int dd = 0; dd != 3; dd++) {
      for (int gg = 0; gg != nb_gauss_pts; ++gg) {
        tAngent(gg, dd) +=
            cblas_ddot(nb_dofs / 3, &diff[gg * nb_approx_fun], 1, dofs[dd], 3);
      }
    }
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
            "This operator can calculate tangent vector only on edge");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSetCovariantPiolaTransformOnEdge::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBeginHot;

  if (type != MBEDGE)
    MoFEMFunctionReturnHot(0);

  FTensor::Index<'i', 3> i;
  FTensor::Tensor1<FTensor::PackPtr<const double *, 0>, 3> t_m(
      &tAngent[0], &tAngent[1], &tAngent[2]);
  const double l0 = t_m(i) * t_m(i);

  auto get_base_at_pts = [&](auto base) {
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_h_curl(
        &data.getN(base)(0, HVEC0), &data.getN(base)(0, HVEC1),
        &data.getN(base)(0, HVEC2));
    return t_h_curl;
  };

  auto get_tangent_at_pts = [&]() {
    FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3> t_m_at_pts(
        &tangentAtGaussPt(0, 0), &tangentAtGaussPt(0, 1),
        &tangentAtGaussPt(0, 2));
    return t_m_at_pts;
  };

  auto calculate_squared_edge_length = [&]() {
    std::vector<double> l1;
    int nb_gauss_pts = tangentAtGaussPt.size1();
    if (nb_gauss_pts) {
      l1.resize(nb_gauss_pts);
      auto t_m_at_pts = get_tangent_at_pts();
      for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
        l1[gg] = t_m_at_pts(i) * t_m_at_pts(i);
        ++t_m_at_pts;
      }
    }
    return l1;
  };

  auto l1 = calculate_squared_edge_length();

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
    const size_t nb_gauss_pts = data.getN(base).size1();
    const size_t nb_dofs = data.getN(base).size2() / 3;
    if (nb_gauss_pts && nb_dofs) {
      auto t_h_curl = get_base_at_pts(base);
      int cc = 0;
      if (tangentAtGaussPt.size1() == nb_gauss_pts) {
        auto t_m_at_pts = get_tangent_at_pts();
        for (int gg = 0; gg != nb_gauss_pts; ++gg) {
          const double l0 = l1[gg];
          for (int ll = 0; ll != nb_dofs; ll++) {
            const double val = t_h_curl(0);
            const double a = val / l0;
            t_h_curl(i) = t_m_at_pts(i) * a;
            ++t_h_curl;
            ++cc;
          }
          ++t_m_at_pts;
        }
      } else {
        for (int gg = 0; gg != nb_gauss_pts; ++gg) {
          for (int ll = 0; ll != nb_dofs; ll++) {
            const double val = t_h_curl(0);
            const double a = val / l0;
            t_h_curl(i) = t_m(i) * a;
            ++t_h_curl;
            ++cc;
          }
        }
      }

      if (cc != nb_gauss_pts * nb_dofs)
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "Data inconsistency");
    }
  }

  MoFEMFunctionReturnHot(0);
}

template <>
template <>
FTensor::Tensor1<double *, 3>
OpGetDataAndGradient<3, 3>::getValAtGaussPtsTensor<3>(MatrixDouble &data) {
  double *ptr = &*data.data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2], 3);
}

template <>
template <>
FTensor::Tensor2<double *, 3, 3>
OpGetDataAndGradient<3, 3>::getGradAtGaussPtsTensor<3, 3>(MatrixDouble &data) {
  double *ptr = &*data.data().begin();
  return FTensor::Tensor2<double *, 3, 3>(ptr, &ptr[1], &ptr[2], &ptr[3],
                                          &ptr[4], &ptr[5], &ptr[6], &ptr[7],
                                          &ptr[8], 9);
}

template <>
MoFEMErrorCode OpGetDataAndGradient<3, 3>::calculateValAndGrad(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBeginHot;
  if (data.getBase() == NOBASE)
    MoFEMFunctionReturnHot(0);
  const unsigned int nb_gauss_pts = data.getN().size1();
  const unsigned int nb_base_functions = data.getN().size2();
  const unsigned int nb_dofs = data.getFieldData().size();
  if (!nb_dofs)
    MoFEMFunctionReturnHot(0);
  auto t_n = data.getFTensor0N();
  auto t_val = getValAtGaussPtsTensor<3>(dataAtGaussPts);
  auto t_grad = getGradAtGaussPtsTensor<3, 3>(dataGradAtGaussPts);
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  if (type == MBVERTEX &&
      data.getDiffN().data().size() == 3 * nb_base_functions) {
    for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
      auto t_data = data.getFTensor1FieldData<3>();
      auto t_diff_n = data.getFTensor1DiffN<3>();
      unsigned int bb = 0;
      for (; bb != nb_dofs / 3; ++bb) {
        t_val(i) += t_data(i) * t_n;
        t_grad(i, j) += t_data(i) * t_diff_n(j);
        ++t_n;
        ++t_diff_n;
        ++t_data;
      }
      ++t_val;
      ++t_grad;
      for (; bb != nb_base_functions; ++bb) {
        ++t_n;
      }
    }
  } else {
    auto t_diff_n = data.getFTensor1DiffN<3>();
    for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
      auto t_data = data.getFTensor1FieldData<3>();
      unsigned int bb = 0;
      for (; bb != nb_dofs / 3; ++bb) {
        t_val(i) += t_data(i) * t_n;
        t_grad(i, j) += t_data(i) * t_diff_n(j);
        ++t_n;
        ++t_diff_n;
        ++t_data;
      }
      ++t_val;
      ++t_grad;
      for (; bb != nb_base_functions; ++bb) {
        ++t_n;
        ++t_diff_n;
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

template <>
MoFEMErrorCode OpGetDataAndGradient<1, 3>::calculateValAndGrad(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBeginHot;
  const unsigned int nb_gauss_pts = data.getN().size1();
  const unsigned int nb_base_functions = data.getN().size2();
  // bool constant_diff = false;
  const unsigned int nb_dofs = data.getFieldData().size();
  auto t_n = data.getFTensor0N();
  FTensor::Tensor0<double *> t_val =
      FTensor::Tensor0<double *>(&*dataAtGaussPts.data().begin(), 1);
  double *ptr = &*dataGradAtGaussPts.data().begin();
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_grad(ptr, &ptr[1],
                                                            &ptr[2]);
  FTensor::Index<'i', 3> i;
  if (type == MBVERTEX &&
      data.getDiffN().data().size() == 3 * nb_base_functions) {
    for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
      auto t_data = data.getFTensor0FieldData();
      auto t_diff_n = data.getFTensor1DiffN<3>();
      unsigned int bb = 0;
      for (; bb != nb_dofs / 3; ++bb) {
        t_val += t_data * t_n;
        t_grad(i) += t_data * t_diff_n(i);
        ++t_n;
        ++t_diff_n;
        ++t_data;
      }
      ++t_val;
      ++t_grad;
      for (; bb != nb_base_functions; ++bb) {
        ++t_n;
      }
    }
  } else {
    auto t_diff_n = data.getFTensor1DiffN<3>();
    for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
      auto t_data = data.getFTensor0FieldData();
      unsigned int bb = 0;
      for (; bb != nb_dofs / 3; ++bb) {
        t_val = t_data * t_n;
        t_grad(i) += t_data * t_diff_n(i);
        ++t_n;
        ++t_diff_n;
        ++t_data;
      }
      ++t_val;
      ++t_grad;
      for (; bb != nb_base_functions; ++bb) {
        ++t_n;
        ++t_diff_n;
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}
} // namespace MoFEM
