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

  if (!invJacPtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "invJacPtr not allocated");

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

      switch (type) {
      case MBVERTEX:
      case MBEDGE:
      case MBTRI:
      case MBTET: {
        for (unsigned int gg = 0; gg < nb_gauss_pts; ++gg) {
          for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
            t_inv_diff_n(i) = t_diff_n(j) * t_inv_jac(j, i);
            ++t_diff_n;
            ++t_inv_diff_n;
          }
          ++t_inv_jac;
        }
      } break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
      }

      diff_n.data().swap(diffNinvJac.data());
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

  if (!invJacPtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer for invJacPtr not allocated");

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

    data.getDiffN(base).data().swap(diffHdivInvJac.data());
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetHOWeigthsOnFace::doWork(int side, EntityType type,
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

  if (!detPtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer for detPtr not allocated");

  const auto nb_integration_pts = detPtr->size();
  if (nb_integration_pts != getGaussPts().size2())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Inconsistent number of data points");

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

  if (!detPtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer for detPtr not allocated");

  if (!jacPtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer for jacPtr not allocated");

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; ++b) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    auto nb_gauss_pts = data.getN(base).size1();
    auto nb_base_functions = data.getN(base).size2() / 3;

    if (data.getDiffN(base).size1() != nb_gauss_pts)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong number integration points");

    if (data.getDiffN(base).size2() / 9 != nb_base_functions)
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Wrong number base functions %d != %d",
               data.getDiffN(base).size2(), nb_base_functions);

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

      data.getN(base).data().swap(piolaN.data());
      data.getDiffN(base).data().swap(piolaDiffN.data());
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

  if (!jacInvPtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer for jacPtr not allocated");

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

    data.getN(base).data().swap(piolaN.data());
    data.getDiffN(base).data().swap(piolaDiffN.data());
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
      if (2 * data.getN().size2() != data.getDiffN().size2()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
      if (nb_dofs % 3 != 0) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
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

      if (cc != nb_gauss_pts * nb_dofs)
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "Data inconsistency");

      baseN.data().swap(piolaN.data());
      diffBaseN.data().swap(diffPiolaN.data());
    }
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM