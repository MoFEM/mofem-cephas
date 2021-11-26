/** \file HODataOperators.cpp

\brief Set of operators for high-order geometry approximation.

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

OpCalculateHOJacVolume::OpCalculateHOJacVolume(
    boost::shared_ptr<MatrixDouble> jac_ptr)
    : VolumeElementForcesAndSourcesCoreBase::UserDataOperator(H1, OPLAST),
      jacPtr(jac_ptr) {

  for (auto t = MBEDGE; t != MBMAXTYPE; ++t)
    doEntities[t] = false;

  if (!jacPtr)
    THROW_MESSAGE("Jac pointer not allocated");
}

MoFEMErrorCode
OpCalculateHOJacVolume::doWork(int side, EntityType type,
                               DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  const auto nb_base_functions = data.getN(NOBASE).size2();
  if (nb_base_functions) {

    const auto nb_gauss_pts = getGaussPts().size2();
    auto t_diff_base = data.getFTensor1DiffN<3>(NOBASE);
    auto &coords = getCoords();

#ifndef NDEBUG
    if (nb_gauss_pts != data.getDiffN().size1())
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Inconsistent number base functions and gauss points");
    if (nb_base_functions != data.getDiffN().size2() / 3)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Inconsistent number of base functions");
    if (coords.size() != 3 * nb_base_functions)
      SETERRQ2(
          PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
          "Number of vertex coordinates and base functions is inconsistent "
          "%d != %d",
          coords.size() / 3, nb_base_functions);
#endif

    jacPtr->resize(9, nb_gauss_pts, false);
    jacPtr->clear();
    auto t_jac = getFTensor2FromMat<3, 3>(*jacPtr);
    auto t_vol_inv_jac = getInvJac();

    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;
    FTensor::Index<'k', 3> k;

    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
      FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_coords(
          &coords[0], &coords[1], &coords[2]);
      for (size_t bb = 0; bb != nb_base_functions; ++bb) {
        t_jac(i, j) += t_coords(i) * (t_vol_inv_jac(k, j) * t_diff_base(k));
        ++t_diff_base;
        ++t_coords;
      }
      ++t_jac;
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpCalculateHOCoords::doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
  FTensor::Index<'i', 3> i;
  MoFEMFunctionBegin;
  const auto nb_dofs = data.getFieldData().size() / 3;
  if (nb_dofs) {
    if (type == MBVERTEX)
      getCoordsAtGaussPts().clear();
    auto t_base = data.getFTensor0N();
    auto t_coords = getFTensor1CoordsAtGaussPts();
    const auto nb_integration_pts = data.getN().size1();
    const auto nb_base_functions = data.getN().size2();
    for (auto gg = 0; gg != nb_integration_pts; ++gg) {
      auto t_dof = data.getFTensor1FieldData<3>();
      size_t bb = 0;
      for (; bb != nb_dofs; ++bb) {
        t_coords(i) += t_base * t_dof(i);
        ++t_dof;
        ++t_base;
      }
      for (; bb < nb_base_functions; ++bb)
        ++t_base;
      ++t_coords;
    }
  }
  MoFEMFunctionReturn(0);
};

MoFEMErrorCode
OpSetHOInvJacToScalarBases::doWork(int side, EntityType type,
                                   DataForcesAndSourcesCore::EntData &data) {
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  MoFEMFunctionBegin;

  auto transform_base = [&](MatrixDouble &diff_n) {
    MoFEMFunctionBeginHot;

    unsigned int nb_gauss_pts = diff_n.size1();
    if (nb_gauss_pts == 0)
      MoFEMFunctionReturnHot(0);

    if (invJacPtr->size2() == nb_gauss_pts) {

      unsigned int nb_base_functions = diff_n.size2() / 3;
      if (nb_base_functions == 0)
        MoFEMFunctionReturnHot(0);

      auto t_inv_jac = getFTensor2FromMat<3, 3>(*invJacPtr);

      diffNinvJac.resize(nb_gauss_pts, 3 * nb_base_functions, false);

      double *t_diff_n_ptr = &*diff_n.data().begin();
      FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_diff_n(
          t_diff_n_ptr, &t_diff_n_ptr[1], &t_diff_n_ptr[2]);
      double *t_inv_n_ptr = &*diffNinvJac.data().begin();
      FTensor::Tensor1<double *, 3> t_inv_diff_n(t_inv_n_ptr, &t_inv_n_ptr[1],
                                                 &t_inv_n_ptr[2], 3);

      for (unsigned int gg = 0; gg < nb_gauss_pts; ++gg) {
        for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
          t_inv_diff_n(i) = t_diff_n(j) * t_inv_jac(j, i);
          ++t_diff_n;
          ++t_inv_diff_n;
        }
        ++t_inv_jac;
      }

      diff_n.swap(diffNinvJac);
    } else {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Wrong number of gauss pts in invJacPtr, is %d but should be %d",
               invJacPtr->size1(), nb_gauss_pts);
    }

    MoFEMFunctionReturnHot(0);
  };

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
    CHKERR transform_base(data.getDiffN(base));
  }

  switch (type) {
  case MBVERTEX:
    for (auto &m : data.getBBDiffNMap())
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
OpSetHOInvJacVectorBase::doWork(int side, EntityType type,
                                DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    diffHdivInvJac.resize(data.getDiffN(base).size1(),
                          data.getDiffN(base).size2(), false);

    unsigned int nb_gauss_pts = data.getDiffN(base).size1();
    unsigned int nb_base_functions = data.getDiffN(base).size2() / 9;
    if (nb_base_functions == 0)
      continue;

    auto t_diff_n = data.getFTensor2DiffN<3, 3>(base);
    double *t_inv_diff_n_ptr = &*diffHdivInvJac.data().begin();
    FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_inv_diff_n(
        t_inv_diff_n_ptr, &t_inv_diff_n_ptr[HVEC0_1],
        &t_inv_diff_n_ptr[HVEC0_2], &t_inv_diff_n_ptr[HVEC1_0],
        &t_inv_diff_n_ptr[HVEC1_1], &t_inv_diff_n_ptr[HVEC1_2],
        &t_inv_diff_n_ptr[HVEC2_0], &t_inv_diff_n_ptr[HVEC2_1],
        &t_inv_diff_n_ptr[HVEC2_2]);
    auto t_inv_jac = getFTensor2FromMat<3, 3>(*invJacPtr);
    for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
      for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
        t_inv_diff_n(i, j) = t_inv_jac(k, j) * t_diff_n(i, k);
        ++t_diff_n;
        ++t_inv_diff_n;
      }
      ++t_inv_jac;
    }

    data.getDiffN(base).swap(diffHdivInvJac);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetHOWeightsOnFace::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  const size_t nb_int_pts = getGaussPts().size2();
  if (getNormalsAtGaussPts().size1()) {
    if (getNormalsAtGaussPts().size1() == nb_int_pts) {
      double a = getMeasure();
      if (getNumeredEntFiniteElementPtr()->getEntType() == MBTRI)
        a *= 2;
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

MoFEMErrorCode OpSetHOWeights::doWork(int side, EntityType type,
                                      DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  const auto nb_integration_pts = detPtr->size();

#ifndef NDEBUG
  if (nb_integration_pts != getGaussPts().size2())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Inconsistent number of data points");
#endif

  auto t_w = getFTensor0IntegrationWeight();
  auto t_det = getFTensor0FromVec(*detPtr);
  for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
    t_w *= t_det;
    ++t_w;
    ++t_det;
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSetHOContravariantPiolaTransform::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

#ifndef NDEBUG
  if (!detPtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer for detPtr not allocated");
  if (!jacPtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer for jacPtr not allocated");
#endif

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; ++b) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    auto nb_gauss_pts = data.getN(base).size1();
    auto nb_base_functions = data.getN(base).size2() / 3;

#ifndef NDEBUG
    if (data.getDiffN(base).size1() != nb_gauss_pts)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong number integration points");

    if (data.getDiffN(base).size2() / 9 != nb_base_functions)
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Wrong number base functions %d != %d",
               data.getDiffN(base).size2(), nb_base_functions);
#endif

    if (nb_gauss_pts && nb_base_functions) {

      piolaN.resize(nb_gauss_pts, data.getN(base).size2(), false);
      piolaDiffN.resize(nb_gauss_pts, data.getDiffN(base).size2(), false);

      auto t_n = data.getFTensor1N<3>(base);
      double *t_transformed_n_ptr = &*piolaN.data().begin();
      FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_transformed_n(
          t_transformed_n_ptr, // HVEC0
          &t_transformed_n_ptr[HVEC1], &t_transformed_n_ptr[HVEC2]);
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

      auto t_det = getFTensor0FromVec(*detPtr);
      auto t_jac = getFTensor2FromMat<3, 3>(*jacPtr);

      for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
        for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
          const double a = 1. / t_det;
          t_transformed_n(i) = a * t_jac(i, k) * t_n(k);
          t_transformed_diff_n(i, k) = a * t_jac(i, j) * t_diff_n(j, k);
          ++t_n;
          ++t_transformed_n;
          ++t_diff_n;
          ++t_transformed_diff_n;
        }
        ++t_det;
        ++t_jac;
      }

      data.getN(base).swap(piolaN);
      data.getDiffN(base).swap(piolaDiffN);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSetHOCovariantPiolaTransform::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    unsigned int nb_gauss_pts = data.getN(base).size1();
    unsigned int nb_base_functions = data.getN(base).size2() / 3;
    piolaN.resize(nb_gauss_pts, data.getN(base).size2(), false);
    piolaDiffN.resize(nb_gauss_pts, data.getDiffN(base).size2(), false);

    auto t_n = data.getFTensor1N<3>(base);
    double *t_transformed_n_ptr = &*piolaN.data().begin();
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_transformed_n(
        t_transformed_n_ptr, // HVEC0
        &t_transformed_n_ptr[HVEC1], &t_transformed_n_ptr[HVEC2]);
    auto t_diff_n = data.getFTensor2DiffN<3, 3>(base);
    double *t_transformed_diff_n_ptr = &*piolaDiffN.data().begin();
    FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_transformed_diff_n(
        t_transformed_diff_n_ptr, &t_transformed_diff_n_ptr[HVEC0_1],
        &t_transformed_diff_n_ptr[HVEC0_2], &t_transformed_diff_n_ptr[HVEC1_0],
        &t_transformed_diff_n_ptr[HVEC1_1], &t_transformed_diff_n_ptr[HVEC1_2],
        &t_transformed_diff_n_ptr[HVEC2_0], &t_transformed_diff_n_ptr[HVEC2_1],
        &t_transformed_diff_n_ptr[HVEC2_2]);

    auto t_inv_jac = getFTensor2FromMat<3, 3>(*jacInvPtr);

    for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
      for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
        t_transformed_n(i) = t_inv_jac(k, i) * t_n(k);
        t_transformed_diff_n(i, k) = t_inv_jac(j, i) * t_diff_n(j, k);
        ++t_n;
        ++t_transformed_n;
        ++t_diff_n;
        ++t_transformed_diff_n;
      }
      ++t_inv_jac;
    }

    data.getN(base).swap(piolaN);
    data.getDiffN(base).swap(piolaDiffN);
  }

  MoFEMFunctionReturn(0);
}

OpCalculateHOJacForFaceImpl<2>::OpCalculateHOJacForFaceImpl(
    boost::shared_ptr<MatrixDouble> jac_ptr)
    : FaceElementForcesAndSourcesCoreBase::UserDataOperator(NOSPACE),
      jacPtr(jac_ptr) {}

MoFEMErrorCode OpCalculateHOJacForFaceImpl<2>::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {

  MoFEMFunctionBegin;

  const auto nb_gauss_pts = getGaussPts().size2();

  jacPtr->resize(4, nb_gauss_pts, false);
  auto t_jac = getFTensor2FromMat<2, 2>(*jacPtr);

  FTensor::Index<'i', 2> i;
  FTensor::Index<'j', 2> j;

  auto t_t1 = getFTensor1Tangent1AtGaussPts();
  auto t_t2 = getFTensor1Tangent2AtGaussPts();
  FTensor::Number<0> N0;
  FTensor::Number<1> N1;

  for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
    t_jac(i, N0) = t_t1(i);
    t_jac(i, N1) = t_t2(i);
    ++t_t1;
    ++t_t2;
    ++t_jac;
  }

  MoFEMFunctionReturn(0);
};

MoFEMErrorCode OpCalculateHOJacForFaceImpl<3>::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  size_t nb_gauss_pts = getGaussPts().size2();
  jacPtr->resize(9, nb_gauss_pts, false);

  auto t_t1 = getFTensor1Tangent1AtGaussPts();
  auto t_t2 = getFTensor1Tangent2AtGaussPts();
  auto t_normal = getFTensor1NormalsAtGaussPts();

  FTensor::Number<0> N0;
  FTensor::Number<1> N1;
  FTensor::Number<2> N2;

  auto t_jac = getFTensor2FromMat<3, 3>(*jacPtr);
  for (size_t gg = 0; gg != nb_gauss_pts; ++gg, ++t_jac) {

    t_jac(i, N0) = t_t1(i);
    t_jac(i, N1) = t_t2(i);

    const double l = sqrt(t_normal(j) * t_normal(j));

    t_jac(i, N2) = t_normal(i) / l;

    ++t_t1;
    ++t_t2;
    ++t_normal;
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpGetHONormalsOnFace::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBeginHot;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  FTensor::Number<0> N0;
  FTensor::Number<1> N1;

  auto get_ftensor1 = [](MatrixDouble &m) {
    return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
        &m(0, 0), &m(0, 1), &m(0, 2));
  };

  unsigned int nb_dofs = data.getFieldData().size();
  if (nb_dofs != 0) {

    int nb_gauss_pts = data.getN().size1();
    auto &tangent1_at_gauss_pts = getTangent1AtGaussPts();
    auto &tangent2_at_gauss_pts = getTangent2AtGaussPts();
    tangent1_at_gauss_pts.resize(nb_gauss_pts, 3, false);
    tangent2_at_gauss_pts.resize(nb_gauss_pts, 3, false);

    switch (type) {
    case MBVERTEX: {
      tangent1_at_gauss_pts.clear();
      tangent2_at_gauss_pts.clear();
    }
    case MBEDGE:
    case MBTRI:
    case MBQUAD: {

#ifndef NDEBUG
      if (2 * data.getN().size2() != data.getDiffN().size2()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
      if (nb_dofs % 3 != 0) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
#endif

      if (nb_dofs > 3 * data.getN().size2()) {
        unsigned int nn = 0;
        for (; nn != nb_dofs; nn++) {
          if (!data.getFieldDofs()[nn]->getActive())
            break;
        }
        if (nn > 3 * data.getN().size2()) {
          SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "data inconsistency for base %s",
                   ApproximationBaseNames[data.getBase()]);
        } else {
          nb_dofs = nn;
          if (!nb_dofs)
            MoFEMFunctionReturnHot(0);
        }
      }
      const int nb_base_functions = data.getN().size2();
      auto t_base = data.getFTensor0N();
      auto t_diff_base = data.getFTensor1DiffN<2>();
      auto t_t1 = get_ftensor1(tangent1_at_gauss_pts);
      auto t_t2 = get_ftensor1(tangent2_at_gauss_pts);
      for (int gg = 0; gg != nb_gauss_pts; ++gg) {
        auto t_data = data.getFTensor1FieldData<3>();
        int bb = 0;
        for (; bb != nb_dofs / 3; ++bb) {
          t_t1(i) += t_data(i) * t_diff_base(N0);
          t_t2(i) += t_data(i) * t_diff_base(N1);
          ++t_data;
          ++t_base;
          ++t_diff_base;
        }
        for (; bb != nb_base_functions; ++bb) {
          ++t_base;
          ++t_diff_base;
        }
        ++t_t1;
        ++t_t2;
      }
    } break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
    }
  }

  if (moab::CN::Dimension(type) == 2) {

    auto &normal_at_gauss_pts = getNormalsAtGaussPts();
    auto &tangent1_at_gauss_pts = getTangent1AtGaussPts();
    auto &tangent2_at_gauss_pts = getTangent2AtGaussPts();

    const auto nb_gauss_pts = tangent1_at_gauss_pts.size1();
    normal_at_gauss_pts.resize(nb_gauss_pts, 3, false);

    auto t_normal = get_ftensor1(normal_at_gauss_pts);
    auto t_t1 = get_ftensor1(tangent1_at_gauss_pts);
    auto t_t2 = get_ftensor1(tangent2_at_gauss_pts);
    for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
      t_normal(j) = FTensor::levi_civita(i, j, k) * t_t1(k) * t_t2(i);
      ++t_normal;
      ++t_t1;
      ++t_t2;
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode OpHOSetContravariantPiolaTransformOnFace3D::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  FTensor::Index<'i', 3> i;
  MoFEMFunctionBegin;

  if (moab::CN::Dimension(type) != 2)
    MoFEMFunctionReturnHot(0);

  auto get_normals_ptr = [&]() {
    if (normalsAtGaussPts)
      return &*normalsAtGaussPts->data().begin();
    else
      return &*getNormalsAtGaussPts().data().begin();
  };

  auto apply_transform_nonlinear_geometry = [&](auto base, auto nb_gauss_pts,
                                                auto nb_base_functions) {
    MoFEMFunctionBegin;

    auto ptr = get_normals_ptr();
    auto t_normal = FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
        &ptr[0], &ptr[1], &ptr[2]);

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

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
    const auto &base_functions = data.getN(base);
    const auto nb_gauss_pts = base_functions.size1();

    if (nb_gauss_pts) {

      const auto nb_base_functions = base_functions.size2() / 3;
      CHKERR apply_transform_nonlinear_geometry(base, nb_gauss_pts,
                                                nb_base_functions);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpHOSetCovariantPiolaTransformOnFace3D::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  const auto type_dim = moab::CN::Dimension(type);
  if (type_dim != 1 && type_dim != 2)
    MoFEMFunctionReturnHot(0);

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 2> k;

  auto get_jac = [&]() {
    if (normalsAtPts && tangent1AtPts && tangent2AtPts) {
      double *ptr_n = &*normalsAtPts->data().begin();
      double *ptr_t1 = &*tangent1AtPts->data().begin();
      double *ptr_t2 = &*tangent2AtPts->data().begin();
      return FTensor::Tensor2<FTensor::PackPtr<double *, 3>, 3, 3>(
          &ptr_t1[0], &ptr_t2[0], &ptr_n[0],

          &ptr_t1[1], &ptr_t2[1], &ptr_n[1],

          &ptr_t1[2], &ptr_t2[2], &ptr_n[2]);
    } else {
      double *ptr_n = &*getNormalsAtGaussPts().data().begin();
      double *ptr_t1 = &*getTangent1AtGaussPts().data().begin();
      double *ptr_t2 = &*getTangent2AtGaussPts().data().begin();
      return FTensor::Tensor2<FTensor::PackPtr<double *, 3>, 3, 3>(
          &ptr_t1[0], &ptr_t2[0], &ptr_n[0],

          &ptr_t1[1], &ptr_t2[1], &ptr_n[1],

          &ptr_t1[2], &ptr_t2[2], &ptr_n[2]);
    }
  };

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; ++b) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    auto &baseN = data.getN(base);
    auto &diffBaseN = data.getDiffN(base);

    int nb_dofs = baseN.size2() / 3;
    int nb_gauss_pts = baseN.size1();

    piolaN.resize(baseN.size1(), baseN.size2());
    diffPiolaN.resize(diffBaseN.size1(), diffBaseN.size2());

    if (nb_dofs > 0 && nb_gauss_pts > 0) {

      auto t_h_curl = data.getFTensor1N<3>(base);
      auto t_diff_h_curl = data.getFTensor2DiffN<3, 2>(base);

      FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_transformed_h_curl(
          &piolaN(0, HVEC0), &piolaN(0, HVEC1), &piolaN(0, HVEC2));

      FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>
          t_transformed_diff_h_curl(
              &diffPiolaN(0, HVEC0_0), &diffPiolaN(0, HVEC0_1),
              &diffPiolaN(0, HVEC1_0), &diffPiolaN(0, HVEC1_1),
              &diffPiolaN(0, HVEC2_0), &diffPiolaN(0, HVEC2_1));

      int cc = 0;
      double det;
      FTensor::Tensor2<double, 3, 3> t_inv_m;

      // HO geometry is set, so jacobian is different at each gauss point
      auto t_m_at_pts = get_jac();
      for (int gg = 0; gg != nb_gauss_pts; ++gg) {
        CHKERR determinantTensor3by3(t_m_at_pts, det);
        CHKERR invertTensor3by3(t_m_at_pts, det, t_inv_m);
        for (int ll = 0; ll != nb_dofs; ll++) {
          t_transformed_h_curl(i) = t_inv_m(j, i) * t_h_curl(j);
          t_transformed_diff_h_curl(i, k) = t_inv_m(j, i) * t_diff_h_curl(j, k);
          ++t_h_curl;
          ++t_transformed_h_curl;
          ++t_diff_h_curl;
          ++t_transformed_diff_h_curl;
          ++cc;
        }
        ++t_m_at_pts;
      }

#ifndef NDEBUG
      if (cc != nb_gauss_pts * nb_dofs)
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "Data inconsistency");
#endif

      baseN.swap(piolaN);
      diffBaseN.swap(diffPiolaN);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpHOSetContravariantPiolaTransformOnEdge3D::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  FTensor::Index<'i', 3> i;
  MoFEMFunctionBeginHot;

  if (type != MBEDGE)
    MoFEMFunctionReturnHot(0);

  const auto nb_gauss_pts = getGaussPts().size2();

  if (tangentAtGaussPts)
    if (tangentAtGaussPts->size1() != nb_gauss_pts)
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Wrong number of integration points %d!=%d",
               tangentAtGaussPts->size1(), nb_gauss_pts);

  auto get_base_at_pts = [&](auto base) {
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_h_curl(
        &data.getN(base)(0, HVEC0), &data.getN(base)(0, HVEC1),
        &data.getN(base)(0, HVEC2));
    return t_h_curl;
  };

  auto get_tangent_ptr = [&]() {
    if (tangentAtGaussPts) {
      return &*tangentAtGaussPts->data().begin();
    } else {
      return &*getTangetAtGaussPts().data().begin();
    }
  };

  auto get_tangent_at_pts = [&]() {
    auto ptr = get_tangent_ptr();
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_m_at_pts(
        &ptr[0], &ptr[1], &ptr[2]);
    return t_m_at_pts;
  };

  auto calculate_squared_edge_length = [&]() {
    std::vector<double> l1;
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
    const auto nb_dofs = data.getN(base).size2() / 3;

    if (nb_gauss_pts && nb_dofs) {
      auto t_h_curl = get_base_at_pts(base);
      int cc = 0;
      auto t_m_at_pts = get_tangent_at_pts();
      for (int gg = 0; gg != nb_gauss_pts; ++gg) {
        const double l0 = l1[gg];
        for (int ll = 0; ll != nb_dofs; ++ll) {
          const double val = t_h_curl(0);
          const double a = val / l0;
          t_h_curl(i) = t_m_at_pts(i) * a;
          ++t_h_curl;
          ++cc;
        }
        ++t_m_at_pts;
      }

#ifndef NDEBUG
      if (cc != nb_gauss_pts * nb_dofs)
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "Data inconsistency");
#endif
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
OpGetHOTangentsOnEdge::doWork(int side, EntityType type,
                              DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  int nb_dofs = data.getFieldData().size();
  if (nb_dofs == 0)
    MoFEMFunctionReturnHot(0);

  auto get_tangent = [&]() -> MatrixDouble & {
    if (tangentsAtPts)
      return *tangentsAtPts;
    else
      return getTangetAtGaussPts();
  };

  auto &tangent = get_tangent();

  int nb_gauss_pts = data.getN().size1();
  tangent.resize(nb_gauss_pts, 3, false);

  int nb_approx_fun = data.getN().size2();
  double *diff = &*data.getDiffN().data().begin();
  double *dofs[] = {&data.getFieldData()[0], &data.getFieldData()[1],
                    &data.getFieldData()[2]};

  tangent.resize(nb_gauss_pts, 3, false);

  switch (type) {
  case MBVERTEX:
    for (int dd = 0; dd != 3; dd++) {
      for (int gg = 0; gg != nb_gauss_pts; ++gg) {
        tangent(gg, dd) = cblas_ddot(2, diff, 1, dofs[dd], 3);
      }
    }
    break;
  case MBEDGE:
#ifndef NDEBUG
    if (nb_dofs % 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
              "Approximated field should be rank 3, i.e. vector in 3d space");
    }
#endif
    for (int dd = 0; dd != 3; dd++) {
      for (int gg = 0; gg != nb_gauss_pts; ++gg) {
        tangent(gg, dd) +=
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

} // namespace MoFEM