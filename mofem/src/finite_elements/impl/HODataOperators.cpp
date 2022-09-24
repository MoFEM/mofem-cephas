/** \file HODataOperators.cpp

\brief Set of operators for high-order geometry approximation.

*/

namespace MoFEM {

OpCalculateHOJacForVolume::OpCalculateHOJacForVolume(
    boost::shared_ptr<MatrixDouble> jac_ptr)
    : VolumeElementForcesAndSourcesCore::UserDataOperator(H1, OPSPACE),
      jacPtr(jac_ptr) {

  for (auto t = MBEDGE; t != MBMAXTYPE; ++t)
    doEntities[t] = false;

  if (!jacPtr)
    THROW_MESSAGE("Jac pointer not allocated");
}

MoFEMErrorCode
OpCalculateHOJacForVolume::doWork(int side, EntityType type,
                                  EntitiesFieldData::EntData &data) {
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
OpSetHOInvJacToScalarBasesImpl::doWork(int side, EntityType type,
                                       EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  if (getFEDim() == 3) {

    auto transform_base = [&](MatrixDouble &diff_n) {
      MoFEMFunctionBeginHot;
      return applyTransform<3, 3, 3, 3>(diff_n);
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

  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Use different operator");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetHOInvJacVectorBase::doWork(int side, EntityType type,
                                EntitiesFieldData::EntData &data) {
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

MoFEMErrorCode OpSetHOWeightsOnFace::doWork(int side, EntityType type,
                                            EntitiesFieldData::EntData &data) {
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
                                      EntitiesFieldData::EntData &data) {
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

MoFEMErrorCode
OpSetHOContravariantPiolaTransform::doWork(int side, EntityType type,
                                           EntitiesFieldData::EntData &data) {
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

MoFEMErrorCode
OpSetHOCovariantPiolaTransform::doWork(int side, EntityType type,
                                       EntitiesFieldData::EntData &data) {
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
    : FaceElementForcesAndSourcesCore::UserDataOperator(NOSPACE),
      jacPtr(jac_ptr) {}

MoFEMErrorCode
OpCalculateHOJacForFaceImpl<2>::doWork(int side, EntityType type,
                                       EntitiesFieldData::EntData &data) {

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

MoFEMErrorCode
OpCalculateHOJacForFaceImpl<3>::doWork(int side, EntityType type,
                                       EntitiesFieldData::EntData &data) {
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

MoFEMErrorCode OpHOSetContravariantPiolaTransformOnFace3D::doWork(
    int side, EntityType type, EntitiesFieldData::EntData &data) {
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
    int side, EntityType type, EntitiesFieldData::EntData &data) {
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
    int side, EntityType type, EntitiesFieldData::EntData &data) {
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
      return &*getTangentAtGaussPts().data().begin();
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

MoFEMErrorCode AddHOOps<2, 2, 2>::add(
    boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::vector<FieldSpace> spaces, std::string geom_field_name) {
  MoFEMFunctionBegin;

  auto det_ptr = boost::make_shared<VectorDouble>();
  auto jac_ptr = boost::make_shared<MatrixDouble>();
  auto inv_jac_ptr = boost::make_shared<MatrixDouble>();

  if (geom_field_name.empty()) {

    pipeline.push_back(new OpCalculateHOJac<2>(jac_ptr));

  } else {

    pipeline.push_back(new OpCalculateHOCoords<2>(geom_field_name));
    pipeline.push_back(
        new OpCalculateVectorFieldGradient<2, 2>(geom_field_name, jac_ptr));
    pipeline.push_back(new OpGetHONormalsOnFace<2>(geom_field_name));
  }

  pipeline.push_back(new OpInvertMatrix<2>(jac_ptr, det_ptr, inv_jac_ptr));
  pipeline.push_back(new OpSetHOWeightsOnFace());

  for (auto s : spaces) {
    switch (s) {
    case NOSPACE:
      break;
    case H1:
    case L2:
      pipeline.push_back(new OpSetHOInvJacToScalarBases<2>(s, inv_jac_ptr));
      break;
    case HCURL:
      MOFEM_TAG_AND_LOG("WORLD", Sev::warning, "AddHOOps<2, 2, 2>")
          << "Missing covariant Piola transform";
      pipeline.push_back(new OpSetInvJacHcurlFace(inv_jac_ptr));
      MOFEM_LOG_CHANNEL("WORLD");
      break;
    case HDIV:
      pipeline.push_back(new OpMakeHdivFromHcurl());
      pipeline.push_back(new OpSetContravariantPiolaTransformOnFace2D(jac_ptr));
      pipeline.push_back(new OpSetInvJacHcurlFace(inv_jac_ptr));
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
               "Space %s not yet implemented", FieldSpaceNames[s]);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode AddHOOps<1, 2, 2>::add(
    boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::vector<FieldSpace> spaces, std::string geom_field_name) {
  MoFEMFunctionBegin;

  if (geom_field_name.empty()) {

  } else {

    pipeline.push_back(new OpCalculateHOCoords<2>(geom_field_name));
    pipeline.push_back(new OpGetHOTangentsOnEdge<2>(geom_field_name));
  }

  for (auto s : spaces) {
    switch (s) {
    case NOSPACE:
      break;
    case HCURL:
      pipeline.push_back(new OpHOSetContravariantPiolaTransformOnEdge3D(HCURL));
      break;
    case HDIV:
      pipeline.push_back(new OpSetContravariantPiolaTransformOnEdge2D());
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
               "Space %s not yet implemented", FieldSpaceNames[s]);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode AddHOOps<1, 3, 3>::add(
    boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::vector<FieldSpace> spaces, std::string geom_field_name) {
  MoFEMFunctionBegin;

  if (geom_field_name.empty()) {

  } else {

    pipeline.push_back(new OpCalculateHOCoords<2>(geom_field_name));
    pipeline.push_back(new OpGetHOTangentsOnEdge<2>(geom_field_name));
  }

  for (auto s : spaces) {
    switch (s) {
    case HCURL:
      pipeline.push_back(new OpHOSetContravariantPiolaTransformOnEdge3D(HCURL));
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
               "Space %s not yet implemented", FieldSpaceNames[s]);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode AddHOOps<2, 3, 3>::add(
    boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::vector<FieldSpace> spaces, std::string geom_field_name) {
  MoFEMFunctionBegin;

  if (geom_field_name.empty()) {
  } else {

    pipeline.push_back(new OpCalculateHOCoords<3>(geom_field_name));
    pipeline.push_back(new OpGetHONormalsOnFace<3>(geom_field_name));
  }

  for (auto s : spaces) {
    switch (s) {
    case HCURL:
      pipeline.push_back(new OpHOSetContravariantPiolaTransformOnFace3D(HDIV));
      break;
    case HDIV:
      pipeline.push_back(new OpHOSetCovariantPiolaTransformOnFace3D(HDIV));
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
               "Space %s not yet implemented", FieldSpaceNames[s]);
      break;
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode AddHOOps<3, 3, 3>::add(
    boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::vector<FieldSpace> spaces, std::string geom_field_name) {
  MoFEMFunctionBegin;

  auto jac = boost::make_shared<MatrixDouble>();
  auto det = boost::make_shared<VectorDouble>();
  auto inv_jac = boost::make_shared<MatrixDouble>();

  if (geom_field_name.empty()) {

    pipeline.push_back(new OpCalculateHOJac<3>(jac));

  } else {

    pipeline.push_back(new OpCalculateHOCoords<3>(geom_field_name));
    pipeline.push_back(
        new OpCalculateVectorFieldGradient<3, 3>(geom_field_name, jac));
  }

  pipeline.push_back(new OpInvertMatrix<3>(jac, det, inv_jac));
  pipeline.push_back(new OpSetHOWeights(det));

  for (auto s : spaces) {
    switch (s) {
    case NOSPACE:
      break;
    case H1:
      pipeline.push_back(new OpSetHOInvJacToScalarBases<3>(H1, inv_jac));
      break;
    case HCURL:
      pipeline.push_back(new OpSetHOCovariantPiolaTransform(HCURL, inv_jac));
      pipeline.push_back(new OpSetHOInvJacVectorBase(HCURL, inv_jac));
      break;
    case HDIV:
      pipeline.push_back(
          new OpSetHOContravariantPiolaTransform(HDIV, det, jac));
      pipeline.push_back(new OpSetHOInvJacVectorBase(HDIV, inv_jac));
      break;
    case L2:
      pipeline.push_back(new OpSetHOInvJacToScalarBases<3>(L2, inv_jac));
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
               "Space %s not yet implemented", FieldSpaceNames[s]);
      break;
    }
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM