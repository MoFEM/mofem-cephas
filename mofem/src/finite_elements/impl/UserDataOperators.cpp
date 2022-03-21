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

MoFEMErrorCode OpSetInvJacSpaceForFaceImpl<2>::doWork(
    int side, EntityType type, EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  if (getFEDim() != 2)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which faces");

  auto apply_transform = [&](MatrixDouble &diff_n) {
    MoFEMFunctionBegin;
    size_t nb_functions = diff_n.size2() / 2;
    if (nb_functions) {
      size_t nb_gauss_pts = diff_n.size1();

#ifndef NDEBUG
      if (nb_gauss_pts != getGaussPts().size2())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Wrong number of Gauss Pts");
#endif

      diffNinvJac.resize(diff_n.size1(), diff_n.size2(), false);
      FTensor::Index<'i', 2> i;
      FTensor::Index<'k', 2> k;
      auto t_diff_n = getFTensor1FromPtr<2>(&*diffNinvJac.data().begin());
      auto t_diff_n_ref = getFTensor1FromPtr<2>(&*diff_n.data().begin());
      auto t_inv_jac = getFTensor2FromMat<2, 2>(*invJacPtr);
      for (size_t gg = 0; gg != nb_gauss_pts; ++gg, ++t_inv_jac) {
        for (size_t dd = 0; dd != nb_functions; ++dd) {
          t_diff_n(i) = t_inv_jac(k, i) * t_diff_n_ref(k);
          ++t_diff_n;
          ++t_diff_n_ref;
        }
      }

      diff_n.swap(diffNinvJac);
    }
    MoFEMFunctionReturn(0);
  };

  if (!(type == MBVERTEX && sPace == L2)) {

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
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSetInvJacSpaceForFaceImpl<3>::doWork(
    int side, EntityType type, EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  if (getFEDim() != 2)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which face");

  auto apply_transform = [&](MatrixDouble &diff_n) {
    MoFEMFunctionBegin;
    size_t nb_functions = diff_n.size2() / 2;
    if (nb_functions) {
      size_t nb_gauss_pts = diff_n.size1();

#ifndef NDEBUG
      if (nb_gauss_pts != getGaussPts().size2())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Wrong number of Gauss Pts");
#endif

      diffNinvJac.resize(diff_n.size1(), nb_functions * 3, false);

      FTensor::Index<'i', 3> i;
      FTensor::Index<'K', 2> K;
      auto t_diff_n = getFTensor1FromPtr<3>(&*diffNinvJac.data().begin());
      auto t_diff_n_ref = getFTensor1FromPtr<2>(&*diff_n.data().begin());
      auto t_inv_jac = getFTensor2FromMat<3, 3>(*invJacPtr);
      for (size_t gg = 0; gg != nb_gauss_pts; ++gg, ++t_inv_jac) {
        for (size_t dd = 0; dd != nb_functions; ++dd) {
          t_diff_n(i) = t_inv_jac(K, i) * t_diff_n_ref(K);
          ++t_diff_n;
          ++t_diff_n_ref;
        }
      }

      diff_n.swap(diffNinvJac);
      
    }
    MoFEMFunctionReturn(0);
  };

  if (!(type == MBVERTEX && sPace == L2)) {
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
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetInvJacHcurlFaceImpl<2>::doWork(int side, EntityType type,
                                    EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  if (type != MBEDGE && type != MBTRI && type != MBQUAD)
    MoFEMFunctionReturnHot(0);

  if (getFEDim() != 2)
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

      auto t_inv_jac = getFTensor2FromMat<2, 2>(*invJacPtr);
      for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg, ++t_inv_jac) {
        for (unsigned int bb = 0; bb != nb_base_functions; bb++) {
          t_inv_diff_n(i, j) = t_diff_n(i, k) * t_inv_jac(k, j);
          ++t_diff_n;
          ++t_inv_diff_n;
        }
      }

      data.getDiffN(base).swap(diffHcurlInvJac);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetInvJacHcurlFaceImpl<3>::doWork(int side, EntityType type,
                                    EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  if (type != MBEDGE && type != MBTRI && type != MBQUAD)
    MoFEMFunctionReturnHot(0);

  if (getFEDim() != 2)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which is triangle");

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'K', 2> K;

  for (int b = AINSWORTH_LEGENDRE_BASE; b != USER_BASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    const unsigned int nb_base_functions = data.getDiffN(base).size2() / 6;
    if (nb_base_functions) {
      const unsigned int nb_gauss_pts = data.getDiffN(base).size1();

      diffHcurlInvJac.resize(nb_gauss_pts, nb_base_functions * 9, false);

      auto t_diff_n = data.getFTensor2DiffN<3, 2>(base);
      double *t_inv_diff_n_ptr = &*diffHcurlInvJac.data().begin();
      FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_inv_diff_n(
          t_inv_diff_n_ptr, &t_inv_diff_n_ptr[HVEC0_1],
          &t_inv_diff_n_ptr[HVEC0_2],

          &t_inv_diff_n_ptr[HVEC1_0], &t_inv_diff_n_ptr[HVEC1_1],
          &t_inv_diff_n_ptr[HVEC1_2],

          &t_inv_diff_n_ptr[HVEC2_0], &t_inv_diff_n_ptr[HVEC2_1],
          &t_inv_diff_n_ptr[HVEC2_2]);

      auto t_inv_jac = getFTensor2FromMat<3, 3>(*invJacPtr);
      for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg, ++t_inv_jac) {
        for (unsigned int bb = 0; bb != nb_base_functions; bb++) {
          t_inv_diff_n(i, j) = t_diff_n(i, K) * t_inv_jac(K, j);
          ++t_diff_n;
          ++t_inv_diff_n;
        }
      }

      data.getDiffN(base).swap(diffHcurlInvJac);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpMakeHdivFromHcurl::doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  if (type != MBEDGE && type != MBTRI && type != MBQUAD)
    MoFEMFunctionReturnHot(0);

  if (getFEDim() != 2)
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

MoFEMErrorCode OpSetContravariantPiolaTransformOnFace2DImpl<2>::doWork(
    int side, EntityType type, EntitiesFieldData::EntData &data) {

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
        auto t_jac = getFTensor2FromMat<2, 2>(*jacPtr);
        for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg, ++t_jac) {
          double det;
          CHKERR determinantTensor2by2(t_jac, det);
          for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
            t_transformed_n(i) = t_jac(i, k) * t_n(k) / det;
            ++t_n;
            ++t_transformed_n;
          }
        }
        data.getN(base).swap(piolaN);
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
        auto t_jac = getFTensor2FromMat<2, 2>(*jacPtr);
        for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg, ++t_jac) {
          double det;
          CHKERR determinantTensor2by2(t_jac, det);
          for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
            t_transformed_diff_n(i, k) = t_jac(i, j) * t_diff_n(j, k) / det;
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

MoFEMErrorCode OpSetContravariantPiolaTransformOnFace2DImpl<3>::doWork(
    int side, EntityType type, EntitiesFieldData::EntData &data) {

  MoFEMFunctionBegin;

  if (type != MBEDGE && type != MBTRI && type != MBQUAD)
    MoFEMFunctionReturnHot(0);

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'K', 2> K;

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
        auto t_jac = getFTensor2FromMat<3, 3>(*jacPtr);
        for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg, ++t_jac) {
          double det;
          CHKERR determinantTensor3by3(t_jac, det);
          for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
            t_transformed_n(i) = t_jac(i, j) * t_n(j) / det;
            ++t_n;
            ++t_transformed_n;
          }
        }
        data.getN(base).swap(piolaN);
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

        auto t_jac = getFTensor2FromMat<3, 3>(*jacPtr);
        for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg, ++t_jac) {
          double det;
          CHKERR determinantTensor3by3(t_jac, det);
          for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
            t_transformed_diff_n(i, K) = t_jac(i, j) * t_diff_n(j, K) / det;
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

MoFEMErrorCode OpSetContravariantPiolaTransformOnEdge2D::doWork(
    int side, EntityType type, EntitiesFieldData::EntData &data) {
  MoFEMFunctionBeginHot;

  if (type != MBEDGE)
    MoFEMFunctionReturnHot(0);

  FTensor::Index<'i', 3> i;

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
    if (nb_dofs > 0) {

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
      }

      if (cc != nb_gauss_pts * nb_dofs)
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "Data inconsistency");
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode OpMultiplyDeterminantOfJacobianAndWeightsForFatPrisms::doWork(
    int side, EntityType type, EntitiesFieldData::EntData &data) {
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
                                     EntitiesFieldData::EntData &data) {
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
                                 EntitiesFieldData::EntData &data) {
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

    data.getDiffN(base).swap(diffNinvJac);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpCalculateInvJacForFlatPrism::doWork(int side, EntityType type,
                                      EntitiesFieldData::EntData &data) {

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
                                  EntitiesFieldData::EntData &data) {
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
      data.getDiffN(base).swap(diffNinvJac);
    } break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
    }
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
