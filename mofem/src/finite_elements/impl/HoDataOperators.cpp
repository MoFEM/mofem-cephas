/** \file HoDataOperators.cpp

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

MoFEMErrorCode OpSetHOWeigthsOnFace::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
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

} // namespace MoFEM