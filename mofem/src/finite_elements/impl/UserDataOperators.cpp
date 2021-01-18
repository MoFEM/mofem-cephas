/** \file UserDataOperators.cpp

\brief Generic user data operators for evaluate fields, and other common
purposes.

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

namespace MoFEM {

MoFEMErrorCode
OpCalculateJacForFace::doWork(int side, EntityType type,
                              DataForcesAndSourcesCore::EntData &data) {

  MoFEMFunctionBegin;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  size_t nb_gauss_pts = getGaussPts().size2();
  auto &coords = getCoords();
  double *coords_ptr = &*coords.data().begin();
  jac.resize(9, nb_gauss_pts, false);
  jac.clear();

  auto cal_jac_on_tri = [&]() {
    MoFEMFunctionBeginHot;

    double *coords_ptr = &*coords.data().begin();

    double j00 = 0, j01 = 0, j10 = 0, j11 = 0, j20 = 0, j21 = 0;
    // this is triangle, derivative of nodal shape functions is constant.
    // So only need to do one node.
    for (auto n : {0, 1, 2}) {

      j00 += coords_ptr[3 * n + 0] * Tools::diffShapeFunMBTRI[2 * n + 0];
      j01 += coords_ptr[3 * n + 0] * Tools::diffShapeFunMBTRI[2 * n + 1];
      j10 += coords_ptr[3 * n + 1] * Tools::diffShapeFunMBTRI[2 * n + 0];
      j11 += coords_ptr[3 * n + 1] * Tools::diffShapeFunMBTRI[2 * n + 1];

      // 3d
      j20 += coords_ptr[3 * n + 2] * Tools::diffShapeFunMBTRI[2 * n + 0];
      j21 += coords_ptr[3 * n + 2] * Tools::diffShapeFunMBTRI[2 * n + 1];
    }

    double normal_ptr[3];
    CHKERR Tools::getTriNormal(coords_ptr, normal_ptr);
    const double l =
        sqrt(normal_ptr[0] * normal_ptr[0] + normal_ptr[1] * normal_ptr[1] +
             normal_ptr[2] * normal_ptr[2]);
    for (auto d : {0, 1, 2})
      normal_ptr[d] /= l;

    auto t_jac = getFaceJac(jac, FTensor::Number<3>());
    for (size_t gg = 0; gg != nb_gauss_pts; ++gg, ++t_jac) {

      t_jac(0, 0) = j00;
      t_jac(0, 1) = j01;
      t_jac(1, 0) = j10;
      t_jac(1, 1) = j11;

      // 3d
      t_jac(2, 0) = j20;
      t_jac(2, 1) = j21;
      for (auto d : {0, 1, 2})
        t_jac(d, d) += normal_ptr[d];
    }
    MoFEMFunctionReturnHot(0);
  };

  auto cal_jac_on_quad = [&]() {
    MoFEMFunctionBeginHot;

    auto t_jac = getFaceJac(jac, FTensor::Number<3>());
    double *ksi_ptr = &getGaussPts()(0, 0);
    double *zeta_ptr = &getGaussPts()(1, 0);
    for (size_t gg = 0; gg != nb_gauss_pts;
         ++gg, ++t_jac, ++ksi_ptr, ++zeta_ptr) {
      const double &ksi = *ksi_ptr;
      const double &zeta = *zeta_ptr;
      t_jac(0, 0) = coords_ptr[3 * 0 + 0] * diffN_MBQUAD0x(zeta) +
                    coords_ptr[3 * 1 + 0] * diffN_MBQUAD1x(zeta) +
                    coords_ptr[3 * 2 + 0] * diffN_MBQUAD2x(zeta) +
                    coords_ptr[3 * 3 + 0] * diffN_MBQUAD3x(zeta);
      t_jac(0, 1) = coords_ptr[3 * 0 + 0] * diffN_MBQUAD0y(ksi) +
                    coords_ptr[3 * 1 + 0] * diffN_MBQUAD1y(ksi) +
                    coords_ptr[3 * 2 + 0] * diffN_MBQUAD2y(ksi) +
                    coords_ptr[3 * 3 + 0] * diffN_MBQUAD3y(ksi);
      t_jac(1, 0) = coords_ptr[3 * 0 + 1] * diffN_MBQUAD0x(zeta) +
                    coords_ptr[3 * 1 + 1] * diffN_MBQUAD1x(zeta) +
                    coords_ptr[3 * 2 + 1] * diffN_MBQUAD2x(zeta) +
                    coords_ptr[3 * 3 + 1] * diffN_MBQUAD3x(zeta);
      t_jac(1, 1) = coords_ptr[3 * 0 + 1] * diffN_MBQUAD0y(ksi) +
                    coords_ptr[3 * 1 + 1] * diffN_MBQUAD1y(ksi) +
                    coords_ptr[3 * 2 + 1] * diffN_MBQUAD2y(ksi) +
                    coords_ptr[3 * 3 + 1] * diffN_MBQUAD3y(ksi);

      t_jac(2, 0) = coords_ptr[3 * 0 + 2] * diffN_MBQUAD0x(zeta) +
                    coords_ptr[3 * 1 + 2] * diffN_MBQUAD1x(zeta) +
                    coords_ptr[3 * 2 + 2] * diffN_MBQUAD2x(zeta) +
                    coords_ptr[3 * 3 + 2] * diffN_MBQUAD3x(zeta);
      t_jac(2, 1) = coords_ptr[3 * 0 + 2] * diffN_MBQUAD0y(ksi) +
                    coords_ptr[3 * 1 + 2] * diffN_MBQUAD1y(ksi) +
                    coords_ptr[3 * 2 + 2] * diffN_MBQUAD2y(ksi) +
                    coords_ptr[3 * 3 + 2] * diffN_MBQUAD3y(ksi);

      FTensor::Tensor1<double, 3> t_t1{t_jac(0, 0), t_jac(1, 0), t_jac(2, 0)};
      FTensor::Tensor1<double, 3> t_t2{t_jac(0, 1), t_jac(1, 1), t_jac(2, 1)};

      FTensor::Tensor1<double, 3> t_normal;
      t_normal(j) = FTensor::levi_civita(i, j, k) * t_t1(k) * t_t2(i);
      t_normal(i) /= sqrt(t_normal(j) * t_normal(j));
      for (auto d : {0, 1, 2})
        t_jac(d, d) += t_normal(d);
    }
    MoFEMFunctionReturnHot(0);
  };

  switch (getNumeredEntFiniteElementPtr()->getEntType()) {
  case MBTRI:
    CHKERR cal_jac_on_tri();
    break;
  case MBQUAD:
    CHKERR cal_jac_on_quad();
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Operator not implemented for this entity type");
  };

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpCalculateInvJacForFace::doWork(int side, EntityType type,
                                 DataForcesAndSourcesCore::EntData &data) {

  MoFEMFunctionBegin;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  size_t nb_gauss_pts = getGaussPts().size2();
  auto &coords = getCoords();
  double *coords_ptr = &*coords.data().begin();
  invJac.resize(9, nb_gauss_pts, false);
  invJac.clear();

  auto cal_inv_jac_on_tri = [&]() {
    MoFEMFunctionBeginHot;
    VectorDouble &coords = getCoords();
    double *coords_ptr = &*coords.data().begin();
    double j00 = 0, j01 = 0, j10 = 0, j11 = 0, j20 = 0, j21 = 0;

    // this is triangle, derivative of nodal shape functions is constant.
    // So only need to do one node.
    for (auto n : {0, 1, 2}) {
      j00 += coords_ptr[3 * n + 0] * Tools::diffShapeFunMBTRI[2 * n + 0];
      j01 += coords_ptr[3 * n + 0] * Tools::diffShapeFunMBTRI[2 * n + 1];
      j10 += coords_ptr[3 * n + 1] * Tools::diffShapeFunMBTRI[2 * n + 0];
      j11 += coords_ptr[3 * n + 1] * Tools::diffShapeFunMBTRI[2 * n + 1];

      // 3d
      j20 += coords_ptr[3 * n + 2] * Tools::diffShapeFunMBTRI[2 * n + 0];
      j21 += coords_ptr[3 * n + 2] * Tools::diffShapeFunMBTRI[2 * n + 1];
    }

    double normal_ptr[3];
    CHKERR Tools::getTriNormal(coords_ptr, normal_ptr);
    const double l =
        sqrt(normal_ptr[0] * normal_ptr[0] + normal_ptr[1] * normal_ptr[1] +
             normal_ptr[2] * normal_ptr[2]);
    for (auto d : {0, 1, 2})
      normal_ptr[d] /= l;

    FTensor::Tensor2<double, 3, 3> t_jac, t_inv_jac;
    t_jac(0, 0) = j00;
    t_jac(0, 1) = j01;
    t_jac(1, 0) = j10;
    t_jac(1, 1) = j11;
    // // 3d
    // t_jac(2, 0) = j20;
    // t_jac(2, 1) = j21;
    // for (auto d : {0, 1, 2})
    //   t_jac(d, d) += normal_ptr[d];

    // double det;
    // CHKERR determinantTensor3by3(t_jac, det);
    // CHKERR invertTensor3by3(t_jac, det, t_inv_jac);

    auto t_inv_jac_at_pts = getFaceJac(invJac, FTensor::Number<2>());

    for (size_t gg = 0; gg != nb_gauss_pts; ++gg, ++t_inv_jac_at_pts) {
      t_inv_jac_at_pts(i, j) = t_inv_jac(i, j);
    }
    MoFEMFunctionReturnHot(0);
  };

  auto cal_inv_jac_on_quad = [&]() {
    MoFEMFunctionBeginHot;
    VectorDouble &coords = getCoords();
    double *coords_ptr = &*coords.data().begin();
    double *ksi_ptr = &getGaussPts()(0, 0);
    double *zeta_ptr = &getGaussPts()(1, 0);
    FTensor::Tensor2<double, 3, 3> t_jac;

    auto t_inv_jac = getFaceJac(invJac, FTensor::Number<3>());
    for (size_t gg = 0; gg != nb_gauss_pts;
         ++gg, ++t_inv_jac, ++ksi_ptr, ++zeta_ptr) {
      const double &ksi = *ksi_ptr;
      const double &zeta = *zeta_ptr;

      t_jac(0, 0) = coords_ptr[3 * 0 + 0] * diffN_MBQUAD0x(zeta) +
                    coords_ptr[3 * 1 + 0] * diffN_MBQUAD1x(zeta) +
                    coords_ptr[3 * 2 + 0] * diffN_MBQUAD2x(zeta) +
                    coords_ptr[3 * 3 + 0] * diffN_MBQUAD3x(zeta);
      t_jac(0, 1) = coords_ptr[3 * 0 + 0] * diffN_MBQUAD0y(ksi) +
                    coords_ptr[3 * 1 + 0] * diffN_MBQUAD1y(ksi) +
                    coords_ptr[3 * 2 + 0] * diffN_MBQUAD2y(ksi) +
                    coords_ptr[3 * 3 + 0] * diffN_MBQUAD3y(ksi);
      t_jac(1, 0) = coords_ptr[3 * 0 + 1] * diffN_MBQUAD0x(zeta) +
                    coords_ptr[3 * 1 + 1] * diffN_MBQUAD1x(zeta) +
                    coords_ptr[3 * 2 + 1] * diffN_MBQUAD2x(zeta) +
                    coords_ptr[3 * 3 + 1] * diffN_MBQUAD3x(zeta);
      t_jac(1, 1) = coords_ptr[3 * 0 + 1] * diffN_MBQUAD0y(ksi) +
                    coords_ptr[3 * 1 + 1] * diffN_MBQUAD1y(ksi) +
                    coords_ptr[3 * 2 + 1] * diffN_MBQUAD2y(ksi) +
                    coords_ptr[3 * 3 + 1] * diffN_MBQUAD3y(ksi);

      t_jac(2, 0) = coords_ptr[3 * 0 + 2] * diffN_MBQUAD0x(zeta) +
                    coords_ptr[3 * 1 + 2] * diffN_MBQUAD1x(zeta) +
                    coords_ptr[3 * 2 + 2] * diffN_MBQUAD2x(zeta) +
                    coords_ptr[3 * 3 + 2] * diffN_MBQUAD3x(zeta);
      t_jac(2, 1) = coords_ptr[3 * 0 + 2] * diffN_MBQUAD0y(ksi) +
                    coords_ptr[3 * 1 + 2] * diffN_MBQUAD1y(ksi) +
                    coords_ptr[3 * 2 + 2] * diffN_MBQUAD2y(ksi) +
                    coords_ptr[3 * 3 + 2] * diffN_MBQUAD3y(ksi);

      FTensor::Tensor1<double, 3> t_t1{t_jac(0, 0), t_jac(1, 0), t_jac(2, 0)};
      FTensor::Tensor1<double, 3> t_t2{t_jac(0, 1), t_jac(1, 1), t_jac(2, 1)};

      FTensor::Tensor1<double, 3> t_normal;
      t_normal(j) = FTensor::levi_civita(i, j, k) * t_t1(k) * t_t2(i);
      t_normal(i) /= sqrt(t_normal(j) * t_normal(j));
      for (auto d : {0, 1, 2})
        t_jac(d, d) += t_normal(d);

      double det;
      CHKERR determinantTensor3by3(t_jac, det);
      CHKERR invertTensor3by3(t_jac, det, t_inv_jac);
   
    }
    MoFEMFunctionReturnHot(0);
  };

  switch (getNumeredEntFiniteElementPtr()->getEntType()) {
  case MBTRI:
    CHKERR cal_inv_jac_on_tri();
    break;
  case MBQUAD:
    CHKERR cal_inv_jac_on_quad();
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Operator not implemented for this entity type");
  };

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetInvJacSpaceForFace::doWork(int side, EntityType type,
                                DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (getNumeredEntFiniteElementPtr()->getEntType() != MBTRI &&
      getNumeredEntFiniteElementPtr()->getEntType() != MBQUAD)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which is triangle");

  auto apply_transform = [&](MatrixDouble &diff_n) {
    MoFEMFunctionBegin;
    size_t nb_functions = diff_n.size2() / 2;
    if (nb_functions) {
      size_t nb_gauss_pts = diff_n.size1();
      diffNinvJac.resize(nb_gauss_pts, 2 * nb_functions, false);

      switch (type) {
      case MBVERTEX:
      case MBEDGE:
      case MBTRI:
      case MBQUAD: {
        FTensor::Index<'i', 2> i;
        FTensor::Index<'j', 2> j;
        FTensor::Index<'k', 2> k;
        FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2> t_diff_n(
            &diffNinvJac(0, 0), &diffNinvJac(0, 1));
        FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2> t_diff_n_ref(
            &diff_n(0, 0), &diff_n(0, 1));
        FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 2, 2> t_inv_jac(
            &invJac(0, 0), &invJac(1, 0), &invJac(2, 0), &invJac(3, 0));
        for (size_t gg = 0; gg != nb_gauss_pts; ++gg, ++t_inv_jac) {
          for (size_t dd = 0; dd != nb_functions; ++dd) {
            t_diff_n(i) = t_inv_jac(k, i) * t_diff_n_ref(k);
            ++t_diff_n;
            ++t_diff_n_ref;
          }
        }
        diff_n.data().swap(diffNinvJac.data());
      } break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
      }
    }
    MoFEMFunctionReturn(0);
  };

  for (int b = AINSWORTH_LEGENDRE_BASE; b != USER_BASE; b++) {
    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
    CHKERR apply_transform(data.getDiffN(base));
  }

  switch (type) {
  case MBVERTEX:
    for (auto &m : data.getBBDiffNMap())
      CHKERR apply_transform(*(m.second));
    break;
  default:
    for (auto &ptr : data.getBBDiffNByOrderArray())
      if (ptr)
        CHKERR apply_transform(*ptr);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetInvJacHcurlFace::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (type != MBEDGE && type != MBTRI && type != MBQUAD)
    MoFEMFunctionReturnHot(0);

  if (getNumeredEntFiniteElementPtr()->getEntType() != MBTRI &&
      getNumeredEntFiniteElementPtr()->getEntType() != MBQUAD)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which is triangle");

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 2> j;
  FTensor::Index<'k', 2> k;

  for (int b = AINSWORTH_LEGENDRE_BASE; b != USER_BASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    const unsigned int nb_base_functions = data.getDiffN(base).size2() / 6;
    if (nb_base_functions) {
      const unsigned int nb_gauss_pts = data.getDiffN(base).size1();

      diffHcurlInvJac.resize(nb_gauss_pts, data.getDiffN(base).size2(), false);

      auto t_diff_n = data.getFTensor2DiffN<3, 2>(base);
      double *t_inv_diff_n_ptr = &*diffHcurlInvJac.data().begin();
      FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2> t_inv_diff_n(
          t_inv_diff_n_ptr, &t_inv_diff_n_ptr[HVEC0_1],

          &t_inv_diff_n_ptr[HVEC1_0], &t_inv_diff_n_ptr[HVEC1_1],

          &t_inv_diff_n_ptr[HVEC2_0], &t_inv_diff_n_ptr[HVEC2_1]);

      FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 2, 2> t_inv_jac(
          &invJac(0, 0), &invJac(1, 0), &invJac(2, 0), &invJac(3, 0));

      for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg, ++t_inv_jac) {
        for (unsigned int bb = 0; bb != nb_base_functions; bb++) {
          t_inv_diff_n(i, j) = t_diff_n(i, k) * t_inv_jac(k, j);
          ++t_diff_n;
          ++t_inv_diff_n;
        }
      }

      data.getDiffN(base).data().swap(diffHcurlInvJac.data());
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpMakeHdivFromHcurl::doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (type != MBEDGE && type != MBTRI && type != MBQUAD)
    MoFEMFunctionReturnHot(0);

  if (getNumeredEntFiniteElementPtr()->getEntType() != MBTRI &&
      getNumeredEntFiniteElementPtr()->getEntType() != MBQUAD)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which is face");

  for (int b = AINSWORTH_LEGENDRE_BASE; b != USER_BASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    const size_t nb_base_functions = data.getN(base).size2() / 3;
    if (nb_base_functions) {

      const size_t nb_gauss_pts = data.getN(base).size1();

      auto t_n = data.getFTensor1N<3>(base);
      auto t_diff_n = data.getFTensor2DiffN<3, 2>(base);
      for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
        for (size_t bb = 0; bb != nb_base_functions; ++bb) {

          const double a = t_n(0);
          t_n(0) = -t_n(1);
          t_n(1) = a;

          for (auto n : {0, 1}) {
            const double b = t_diff_n(0, n);
            t_diff_n(0, n) = -t_diff_n(1, n);
            t_diff_n(1, n) = b;
          }

          ++t_n;
          ++t_diff_n;
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpMakeHighOrderGeometryWeightsOnFace::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  const size_t nb_int_pts = getGaussPts().size2();
  if (getNormalsAtGaussPts().size1()) {
    if (getNormalsAtGaussPts().size1() == nb_int_pts) {
      const double a = getMeasure();
      auto t_w = getFTensor0IntegrationWeight();
      auto t_normal = getFTensor1NormalsAtGaussPts();
      FTensor::Index<'i', 3> i;
      for (size_t gg = 0; gg != nb_int_pts; ++gg) {
        t_w *= sqrt(t_normal(i) * t_normal(i)) / a;
        ++t_w;
        ++t_normal;
      }
    } else {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
               "Number of rows in getNormalsAtGaussPts should be equal to "
               "number of integration points, but is not, i.e. %d != %d",
               getNormalsAtGaussPts().size1(), nb_int_pts);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSetContravariantPiolaTransformFace ::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {

  MoFEMFunctionBegin;

  if (type != MBEDGE && type != MBTRI && type != MBQUAD)
    MoFEMFunctionReturnHot(0);

  FTensor::Index<'i', 2> i;
  FTensor::Index<'j', 2> j;
  FTensor::Index<'k', 2> k;

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    const size_t nb_base_functions = data.getN(base).size2() / 3;
    if (nb_base_functions) {

      const size_t nb_gauss_pts = data.getN(base).size1();
      piolaN.resize(nb_gauss_pts, data.getN(base).size2(), false);
      if (data.getN(base).size2() > 0) {
        auto t_n = data.getFTensor1N<3>(base);
        double *t_transformed_n_ptr = &*piolaN.data().begin();
        FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_transformed_n(
            t_transformed_n_ptr, // HVEC0
            &t_transformed_n_ptr[HVEC1], &t_transformed_n_ptr[HVEC2]);
        FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 2, 2> t_jac(
            &jAc(0, 0), &jAc(1, 0), &jAc(2, 0), &jAc(3, 0));
        for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg, ++t_jac) {
          double det;
          CHKERR determinantTensor2by2(t_jac, det);
          for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
            t_transformed_n(i) = t_jac(i, k) * t_n(k) / det;
            ++t_n;
            ++t_transformed_n;
          }
        }
        data.getN(base).data().swap(piolaN.data());
      }

      piolaDiffN.resize(nb_gauss_pts, data.getDiffN(base).size2(), false);
      if (data.getDiffN(base).size2() > 0) {
        auto t_diff_n = data.getFTensor2DiffN<3, 2>(base);
        double *t_transformed_diff_n_ptr = &*piolaDiffN.data().begin();
        FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>
            t_transformed_diff_n(t_transformed_diff_n_ptr,
                                 &t_transformed_diff_n_ptr[HVEC0_1],

                                 &t_transformed_diff_n_ptr[HVEC1_0],
                                 &t_transformed_diff_n_ptr[HVEC1_1],

                                 &t_transformed_diff_n_ptr[HVEC2_0],
                                 &t_transformed_diff_n_ptr[HVEC2_1]);
        FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 2, 2> t_jac(
            &jAc(0, 0), &jAc(1, 0), &jAc(2, 0), &jAc(3, 0));
        for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg, ++t_jac) {
          double det;
          CHKERR determinantTensor2by2(t_jac, det);
          for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
            t_transformed_diff_n(i, k) = t_jac(i, j) * t_diff_n(j, k) / det;
            ++t_diff_n;
            ++t_transformed_diff_n;
          }
        }
        data.getDiffN(base).data().swap(piolaDiffN.data());
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSetContravariantPiolaTransformOnEdge::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBeginHot;

  if (type != MBEDGE)
    MoFEMFunctionReturnHot(0);

  const auto &edge_direction = getDirection();

  FTensor::Index<'i', 3> i;
  FTensor::Tensor1<double, 3> t_m(-edge_direction[1], edge_direction[0],
                                  edge_direction[2]);
  const double l0 = t_m(i) * t_m(i);

  std::vector<double> l1;
  {
    int nb_gauss_pts = getTangetAtGaussPts().size1();
    if (nb_gauss_pts) {
      l1.resize(nb_gauss_pts);
      const auto &edge_direction = getTangetAtGaussPts();
      FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3> t_m_at_pts(
          &edge_direction(0, 0), &edge_direction(0, 1), &edge_direction(0, 2));
      for (int gg = 0; gg < nb_gauss_pts; ++gg) {
        l1[gg] = t_m_at_pts(i) * t_m_at_pts(i);
        ++t_m_at_pts;
      }
    }
  }

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
    size_t nb_gauss_pts = data.getN(base).size1();
    size_t nb_dofs = data.getN(base).size2() / 3;
    if (nb_gauss_pts > 0 && nb_dofs > 0) {

      auto t_h_div = data.getFTensor1N<3>(base);

      size_t cc = 0;
      const auto &edge_direction_at_gauss_pts = getTangetAtGaussPts();
      if (edge_direction_at_gauss_pts.size1() == nb_gauss_pts) {

        FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3> t_m_at_pts(
            &edge_direction_at_gauss_pts(0, 1),
            &edge_direction_at_gauss_pts(0, 0),
            &edge_direction_at_gauss_pts(0, 2));

        for (int gg = 0; gg != nb_gauss_pts; ++gg) {
          const double l0 = l1[gg];
          for (int ll = 0; ll != nb_dofs; ++ll) {
            const double val = t_h_div(0);
            const double a = val / l0;
            t_h_div(i) = t_m_at_pts(i) * a;
            t_h_div(0) *= -1;
            ++t_h_div;
            ++cc;
          }
          ++t_m_at_pts;
        }

      } else {

        for (int gg = 0; gg != nb_gauss_pts; ++gg) {
          for (int ll = 0; ll != nb_dofs; ll++) {
            const double val = t_h_div(0);
            const double a = val / l0;
            t_h_div(i) = t_m(i) * a;
            ++t_h_div;
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

MoFEMErrorCode OpMultiplyDeterminantOfJacobianAndWeightsForFatPrisms::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (type == MBVERTEX) {

    VectorDouble &coords = getCoords();
    double *coords_ptr = &*coords.data().begin();

    const int nb_gauss_pts = data.getN(NOBASE).size1();
    auto t_diff_n = data.getFTensor1DiffN<3>(NOBASE);

    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;
    FTensor::Tensor2<double, 3, 3> t_jac;

    auto t_w = getFTensor0IntegrationWeight();
    for (int gg = 0; gg != nb_gauss_pts; gg++) {

      FTensor::Tensor1<double *, 3> t_coords(coords_ptr, &coords_ptr[1],
                                             &coords_ptr[2], 3);
      t_jac(i, j) = 0;
      for (int bb = 0; bb != 6; bb++) {
        t_jac(i, j) += t_coords(i) * t_diff_n(j);
        ++t_diff_n;
        ++t_coords;
      }

      double det;
      CHKERR determinantTensor3by3(t_jac, det);
      t_w *= det / 2.;

      ++t_w;
    }

    double &vol = getVolume();
    auto t_w_scaled = getFTensor0IntegrationWeight();
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      t_w_scaled /= vol;
      ++t_w_scaled;
    }
  }

  doEntities[MBVERTEX] = true;
  std::fill(&doEntities[MBEDGE], &doEntities[MBMAXTYPE], false);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpCalculateInvJacForFatPrism::doWork(int side, EntityType type,
                                     DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (type == MBVERTEX) {

    VectorDouble &coords = getCoords();
    double *coords_ptr = &*coords.data().begin();

    const int nb_gauss_pts = data.getN(NOBASE).size1();
    auto t_diff_n = data.getFTensor1DiffN<3>(NOBASE);
    invJac.resize(9, nb_gauss_pts, false);
    invJac.clear();
    auto t_inv_jac = getFTensor2FromMat<3, 3>(invJac);

    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;
    FTensor::Tensor2<double, 3, 3> t_jac;

    auto t_w = getFTensor0IntegrationWeight();
    for (int gg = 0; gg != nb_gauss_pts; gg++) {

      FTensor::Tensor1<double *, 3> t_coords(coords_ptr, &coords_ptr[1],
                                             &coords_ptr[2], 3);
      t_jac(i, j) = 0;
      for (int bb = 0; bb != 6; bb++) {
        t_jac(i, j) += t_coords(i) * t_diff_n(j);
        ++t_diff_n;
        ++t_coords;
      }

      double det;
      CHKERR determinantTensor3by3(t_jac, det);
      CHKERR invertTensor3by3(t_jac, det, t_inv_jac);
      ++t_inv_jac;
    }
  }

  doEntities[MBVERTEX] = true;
  std::fill(&doEntities[MBEDGE], &doEntities[MBMAXTYPE], false);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetInvJacH1ForFatPrism::doWork(int side, EntityType type,
                                 DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  for (int b = AINSWORTH_LEGENDRE_BASE; b != USER_BASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
    if (data.getN(base).size2() == 0)
      continue;

    const int nb_gauss_pts = data.getN(base).size1();
    auto t_diff_n = data.getFTensor1DiffN<3>(base);
    diffNinvJac.resize(data.getDiffN(base).size1(), data.getDiffN(base).size2(),
                       false);

    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;

    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_inv_diff_n(
        &diffNinvJac(0, 0), &diffNinvJac(0, 1), &diffNinvJac(0, 2));
    auto t_inv_jac = getFTensor2FromMat<3, 3>(invJac);

    const int nb_dofs = data.getN(base).size2();
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      for (int bb = 0; bb != nb_dofs; bb++) {
        t_inv_diff_n(i) = t_diff_n(j) * t_inv_jac(j, i);
        ++t_inv_diff_n;
        ++t_diff_n;
      }
      ++t_inv_jac;
    }

    data.getDiffN(base).data().swap(diffNinvJac.data());
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpCalculateInvJacForFlatPrism::doWork(int side, EntityType type,
                                      DataForcesAndSourcesCore::EntData &data) {

  MoFEMFunctionBegin;

  if (type == MBVERTEX) {

    VectorDouble &coords = getCoords();
    double *coords_ptr = &*coords.data().begin();
    double diff_n[6];
    CHKERR ShapeDiffMBTRI(diff_n);
    double j00_f3, j01_f3, j10_f3, j11_f3;
    for (int gg = 0; gg < 1; gg++) {
      // this is triangle, derivative of nodal shape functions is constant.
      // So only need to do one node.
      j00_f3 = cblas_ddot(3, &coords_ptr[0], 3, &diff_n[0], 2);
      j01_f3 = cblas_ddot(3, &coords_ptr[0], 3, &diff_n[1], 2);
      j10_f3 = cblas_ddot(3, &coords_ptr[1], 3, &diff_n[0], 2);
      j11_f3 = cblas_ddot(3, &coords_ptr[1], 3, &diff_n[1], 2);
    }
    double det_f3 = j00_f3 * j11_f3 - j01_f3 * j10_f3;
    invJacF3.resize(2, 2, false);
    invJacF3(0, 0) = j11_f3 / det_f3;
    invJacF3(0, 1) = -j01_f3 / det_f3;
    invJacF3(1, 0) = -j10_f3 / det_f3;
    invJacF3(1, 1) = j00_f3 / det_f3;
  }

  doEntities[MBVERTEX] = true;
  std::fill(&doEntities[MBEDGE], &doEntities[MBMAXTYPE], false);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetInvJacH1ForFlatPrism::doWork(int side, EntityType type,
                                  DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  for (int b = AINSWORTH_LEGENDRE_BASE; b != USER_BASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    unsigned int nb_dofs = data.getN(base).size2();
    if (nb_dofs == 0)
      MoFEMFunctionReturnHot(0);
    unsigned int nb_gauss_pts = data.getN(base).size1();
    diffNinvJac.resize(nb_gauss_pts, 2 * nb_dofs, false);

    if (type != MBVERTEX) {
      if (nb_dofs != data.getDiffN(base).size2() / 2) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "data inconsistency nb_dofs != data.diffN.size2()/2 ( %u != "
                 "%u/2 )",
                 nb_dofs, data.getDiffN(base).size2());
      }
    }

    switch (type) {
    case MBVERTEX:
    case MBEDGE:
    case MBTRI: {
      for (unsigned int gg = 0; gg < nb_gauss_pts; gg++) {
        for (unsigned int dd = 0; dd < nb_dofs; dd++) {
          cblas_dgemv(CblasRowMajor, CblasTrans, 2, 2, 1,
                      &*invJacF3.data().begin(), 2,
                      &data.getDiffN(base)(gg, 2 * dd), 1, 0,
                      &diffNinvJac(gg, 2 * dd), 1);
        }
      }
      data.getDiffN(base).data().swap(diffNinvJac.data());
    } break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
    }
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
