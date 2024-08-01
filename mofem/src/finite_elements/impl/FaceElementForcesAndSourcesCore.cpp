/** \file FaceElementForcesAndSourcesCore.cpp

\brief Implementation of face element

*/

namespace MoFEM {

FaceElementForcesAndSourcesCore::FaceElementForcesAndSourcesCore(
    Interface &m_field)
    : ForcesAndSourcesCore(m_field),
      meshPositionsFieldName("MESH_NODE_POSITIONS"), aRea(elementMeasure) {}

MoFEMErrorCode
FaceElementForcesAndSourcesCore::calculateAreaAndNormalAtIntegrationPts() {
  MoFEMFunctionBegin;

  auto type = numeredEntFiniteElementPtr->getEntType();

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  auto get_ftensor_from_vec_3d = [](VectorDouble &v) {
    return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(&v[0], &v[1],
                                                              &v[2]);
  };

  auto get_ftensor_n_diff = [&]() {
    const auto &m = dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE);
    return FTensor::Tensor1<FTensor::PackPtr<const double *, 2>, 2>(&m(0, 0),
                                                                    &m(0, 1));
  };

  auto get_ftensor_from_mat_3d = [](MatrixDouble &m) {
    return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
        &m(0, 0), &m(0, 1), &m(0, 2));
  };

  if (type == MBTRI) {

    const size_t nb_gauss_pts = gaussPts.size2();
    normalsAtGaussPts.resize(nb_gauss_pts, 3);
    tangentOneAtGaussPts.resize(nb_gauss_pts, 3);
    tangentTwoAtGaussPts.resize(nb_gauss_pts, 3);

    auto t_tan1 = get_ftensor_from_mat_3d(tangentOneAtGaussPts);
    auto t_tan2 = get_ftensor_from_mat_3d(tangentTwoAtGaussPts);
    auto t_normal = get_ftensor_from_mat_3d(normalsAtGaussPts);

    auto t_n =
        FTensor::Tensor1<double *, 3>(&nOrmal[0], &nOrmal[1], &nOrmal[2]);
    auto t_t1 = FTensor::Tensor1<double *, 3>(&tangentOne[0], &tangentOne[1],
                                              &tangentOne[2]);
    auto t_t2 = FTensor::Tensor1<double *, 3>(&tangentTwo[0], &tangentTwo[1],
                                              &tangentTwo[2]);

    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      t_normal(i) = t_n(i);
      t_tan1(i) = t_t1(i);
      t_tan2(i) = t_t2(i);
      ++t_tan1;
      ++t_tan2;
      ++t_normal;
    }

  } else if (type == MBQUAD) {

    EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
    CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
    coords.resize(num_nodes * 3, false);
    CHKERR mField.get_moab().get_coords(conn, num_nodes,
                                        &*coords.data().begin());

    const size_t nb_gauss_pts = gaussPts.size2();
    normalsAtGaussPts.resize(nb_gauss_pts, 3);
    tangentOneAtGaussPts.resize(nb_gauss_pts, 3);
    tangentTwoAtGaussPts.resize(nb_gauss_pts, 3);
    normalsAtGaussPts.clear();
    tangentOneAtGaussPts.clear();
    tangentTwoAtGaussPts.clear();

    auto t_t1 = get_ftensor_from_mat_3d(tangentOneAtGaussPts);
    auto t_t2 = get_ftensor_from_mat_3d(tangentTwoAtGaussPts);
    auto t_normal = get_ftensor_from_mat_3d(normalsAtGaussPts);

    FTensor::Number<0> N0;
    FTensor::Number<1> N1;

    auto t_diff = get_ftensor_n_diff();
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      auto t_coords = get_ftensor_from_vec_3d(coords);
      for (int nn = 0; nn != num_nodes; ++nn) {
        t_t1(i) += t_coords(i) * t_diff(N0);
        t_t2(i) += t_coords(i) * t_diff(N1);
        ++t_diff;
        ++t_coords;
      }
      t_normal(j) = FTensor::levi_civita(i, j, k) * t_t1(k) * t_t2(i);

      ++t_t1;
      ++t_t2;
      ++t_normal;
    }
  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Element type not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FaceElementForcesAndSourcesCore::calculateAreaAndNormal() {
  MoFEMFunctionBegin;

  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
  coords.resize(num_nodes * 3, false);
  CHKERR mField.get_moab().get_coords(conn, num_nodes, &*coords.data().begin());
  nOrmal.resize(3, false);
  tangentOne.resize(3, false);
  tangentTwo.resize(3, false);

  auto calc_normal = [&](const double *diff_ptr) {
    MoFEMFunctionBegin;
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_coords(
        &coords[0], &coords[1], &coords[2]);
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_normal(
        &nOrmal[0], &nOrmal[1], &nOrmal[2]);
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_t1(
        &tangentOne[0], &tangentOne[1], &tangentOne[2]);
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_t2(
        &tangentTwo[0], &tangentTwo[1], &tangentTwo[2]);
    FTensor::Tensor1<FTensor::PackPtr<const double *, 2>, 2> t_diff(
        diff_ptr, &diff_ptr[1]);

    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;
    FTensor::Index<'k', 3> k;

    FTensor::Number<0> N0;
    FTensor::Number<1> N1;
    t_t1(i) = 0;
    t_t2(i) = 0;

    for (int nn = 0; nn != num_nodes; ++nn) {
      t_t1(i) += t_coords(i) * t_diff(N0);
      t_t2(i) += t_coords(i) * t_diff(N1);
      ++t_coords;
      ++t_diff;
    }
    t_normal(j) = FTensor::levi_civita(i, j, k) * t_t1(k) * t_t2(i);
    aRea = sqrt(t_normal(i) * t_normal(i));
    MoFEMFunctionReturn(0);
  };

  const double *diff_ptr;
  switch (numeredEntFiniteElementPtr->getEntType()) {
  case MBTRI:
    diff_ptr = Tools::diffShapeFunMBTRI.data();
    CHKERR calc_normal(diff_ptr);
    // FIXME: Normal should be divided not the area for triangle!!
    aRea /= 2;
    break;
  case MBQUAD:
    diff_ptr = Tools::diffShapeFunMBQUADAtCenter.data();
    CHKERR calc_normal(diff_ptr);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Element type not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FaceElementForcesAndSourcesCore::setIntegrationPts() {
  MoFEMFunctionBegin;
  // Set integration points
  int order_data = getMaxDataOrder();
  int order_row = getMaxRowOrder();
  int order_col = getMaxColOrder();

  const auto type = numeredEntFiniteElementPtr->getEntType();

  auto get_rule_by_type = [&]() {
    switch (type) {
    case MBQUAD:
      return getRule(order_row + 1, order_col + 1, order_data + 1);
    default:
      return getRule(order_row, order_col, order_data);
    }
  };

  const int rule = get_rule_by_type();

  auto set_integration_pts_for_tri = [&]() {
    MoFEMFunctionBegin;
    if (rule < QUAD_2D_TABLE_SIZE) {
      if (QUAD_2D_TABLE[rule]->dim != 2) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong dimension");
      }
      if (QUAD_2D_TABLE[rule]->order < rule) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong order %d != %d", QUAD_2D_TABLE[rule]->order, rule);
      }
      const size_t nb_gauss_pts = QUAD_2D_TABLE[rule]->npoints;
      gaussPts.resize(3, nb_gauss_pts, false);
      cblas_dcopy(nb_gauss_pts, &QUAD_2D_TABLE[rule]->points[1], 3,
                  &gaussPts(0, 0), 1);
      cblas_dcopy(nb_gauss_pts, &QUAD_2D_TABLE[rule]->points[2], 3,
                  &gaussPts(1, 0), 1);
      cblas_dcopy(nb_gauss_pts, QUAD_2D_TABLE[rule]->weights, 1,
                  &gaussPts(2, 0), 1);
      dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 3,
                                                             false);
      double *shape_ptr =
          &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      cblas_dcopy(3 * nb_gauss_pts, QUAD_2D_TABLE[rule]->points, 1, shape_ptr,
                  1);
      dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).resize(3, 2, false);
      std::copy(
          Tools::diffShapeFunMBTRI.begin(), Tools::diffShapeFunMBTRI.end(),
          dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).data().begin());
    } else {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "rule > quadrature order %d < %d", rule, QUAD_2D_TABLE_SIZE);
    }
    MoFEMFunctionReturn(0);
  };

  auto calc_base_for_tri = [&]() {
    MoFEMFunctionBegin;
    const size_t nb_gauss_pts = gaussPts.size2();
    auto &base = dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE);
    auto &diff_base = dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE);
    base.resize(nb_gauss_pts, 3, false);
    diff_base.resize(3, 2, false);
    CHKERR ShapeMBTRI(&*base.data().begin(), &gaussPts(0, 0), &gaussPts(1, 0),
                      nb_gauss_pts);
    std::copy(
        Tools::diffShapeFunMBTRI.begin(), Tools::diffShapeFunMBTRI.end(),
        dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).data().begin());
    MoFEMFunctionReturn(0);
  };

  auto calc_base_for_quad = [&]() {
    MoFEMFunctionBegin;
    const size_t nb_gauss_pts = gaussPts.size2();
    auto &base = dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE);
    auto &diff_base = dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE);
    base.resize(nb_gauss_pts, 4, false);
    diff_base.resize(nb_gauss_pts, 8, false);
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      const double ksi = gaussPts(0, gg);
      const double zeta = gaussPts(1, gg);
      base(gg, 0) = N_MBQUAD0(ksi, zeta);
      base(gg, 1) = N_MBQUAD1(ksi, zeta);
      base(gg, 2) = N_MBQUAD2(ksi, zeta);
      base(gg, 3) = N_MBQUAD3(ksi, zeta);
      diff_base(gg, 0) = diffN_MBQUAD0x(zeta);
      diff_base(gg, 1) = diffN_MBQUAD0y(ksi);
      diff_base(gg, 2) = diffN_MBQUAD1x(zeta);
      diff_base(gg, 3) = diffN_MBQUAD1y(ksi);
      diff_base(gg, 4) = diffN_MBQUAD2x(zeta);
      diff_base(gg, 5) = diffN_MBQUAD2y(ksi);
      diff_base(gg, 6) = diffN_MBQUAD3x(zeta);
      diff_base(gg, 7) = diffN_MBQUAD3y(ksi);
    }
    MoFEMFunctionReturn(0);
  };

  if (rule >= 0) {
    switch (type) {
    case MBTRI:
      CHKERR set_integration_pts_for_tri();
      break;
    case MBQUAD:
      CHKERR Tools::outerProductOfEdgeIntegrationPtsForQuad(gaussPts, rule,
                                                            rule);
      CHKERR calc_base_for_quad();
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
               "Element type not implemented: %d", type);
    }

  } else {
    // If rule is negative, set user defined integration points
    CHKERR setGaussPts(order_row, order_col, order_data);
    const size_t nb_gauss_pts = gaussPts.size2();
    if (nb_gauss_pts) {
      switch (type) {
      case MBTRI:
        CHKERR calc_base_for_tri();
        break;
      case MBQUAD:
        CHKERR calc_base_for_quad();
        break;
      default:
        SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                 "Element type not implemented: %d", type);
      }
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
FaceElementForcesAndSourcesCore::getSpaceBaseAndOrderOnElement() {
  MoFEMFunctionBegin;
  // Get spaces order/base and sense of entities.

  CHKERR getSpacesAndBaseOnEntities(dataH1);

  auto type = numeredEntFiniteElementPtr->getEntType();
  auto dim_type = CN::Dimension(type);

  auto get_data_on_ents = [&](auto lower_dim, auto space) {
    MoFEMFunctionBeginHot;
    auto data = dataOnElement[space];
    for (auto dd = dim_type; dd >= lower_dim; --dd) {
      int nb_ents = moab::CN::NumSubEntities(type, dd);
      for (int ii = 0; ii != nb_ents; ++ii) {
        auto sub_ent_type = moab::CN::SubEntityType(type, dd, ii);
        if ((dataH1.spacesOnEntities[sub_ent_type]).test(space)) {
          auto &data_on_ent = data->dataOnEntities[sub_ent_type];
          CHKERR getEntitySense(sub_ent_type, data_on_ent);
          CHKERR getEntityDataOrder(sub_ent_type, space, data_on_ent);
          data->spacesOnEntities[sub_ent_type].set(space);
        }
      }
    }
    MoFEMFunctionReturnHot(0);
  };

  CHKERR get_data_on_ents(1, H1);    // H1
  CHKERR get_data_on_ents(1, HCURL); // Hcurl
  CHKERR get_data_on_ents(2, HDIV);  // Hdiv
  CHKERR get_data_on_ents(2, L2);    // L2

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
FaceElementForcesAndSourcesCore::calculateCoordinatesAtGaussPts() {
  MoFEMFunctionBeginHot;

  const size_t nb_nodes =
      dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).size2();
  double *shape_functions =
      &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
  const size_t nb_gauss_pts = gaussPts.size2();
  coordsAtGaussPts.resize(nb_gauss_pts, 3, false);
  for (int gg = 0; gg != nb_gauss_pts; ++gg)
    for (int dd = 0; dd != 3; ++dd)
      coordsAtGaussPts(gg, dd) = cblas_ddot(
          nb_nodes, &shape_functions[nb_nodes * gg], 1, &coords[dd], 3);

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode FaceElementForcesAndSourcesCore::UserDataOperator::setPtrFE(
    ForcesAndSourcesCore *ptr) {
  MoFEMFunctionBeginHot;
  if (!(ptrFE = dynamic_cast<FaceElementForcesAndSourcesCore *>(ptr)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "User operator and finite element do not work together");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode FaceElementForcesAndSourcesCore::operator()() {
  MoFEMFunctionBegin;

  const auto type = numeredEntFiniteElementPtr->getEntType();
  if (type != lastEvaluatedElementEntityType) {
    switch (type) {
    case MBTRI:
      getElementPolynomialBase() =
          boost::shared_ptr<BaseFunction>(new TriPolynomialBase());
      break;
    case MBQUAD:
      getElementPolynomialBase() =
          boost::shared_ptr<BaseFunction>(new QuadPolynomialBase());
      break;
    default:
      MoFEMFunctionReturnHot(0);
    }
    CHKERR createDataOnElement(type);
    lastEvaluatedElementEntityType = type;
  }

  // Calculate normal and tangent vectors for face geometry
  CHKERR calculateAreaAndNormal();
  CHKERR getSpaceBaseAndOrderOnElement();

  CHKERR setIntegrationPts();
  if (gaussPts.size2() == 0)
    MoFEMFunctionReturnHot(0);

  CHKERR calculateCoordinatesAtGaussPts();
  CHKERR calHierarchicalBaseFunctionsOnElement();
  CHKERR calBernsteinBezierBaseFunctionsOnElement();
  CHKERR calculateAreaAndNormalAtIntegrationPts();

  // Iterate over operators
  CHKERR loopOverOperators();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
FaceElementForcesAndSourcesCore::UserDataOperator::loopSideVolumes(
    const string fe_name, VolumeElementForcesAndSourcesCoreOnSide &fe_method) {
  return loopSide(fe_name, &fe_method, 3);
}

MoFEMErrorCode OpCopyGeomDataToE<2>::doWork(int side, EntityType type,
                                            EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

#ifndef NDEBUG
  if (toElePtr->gaussPts.size1() != getGaussPts().size1()) {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Inconsistent numer of weights %d != %d",
             toElePtr->gaussPts.size1(), getGaussPts().size1());
  }
  if (toElePtr->gaussPts.size2() != getGaussPts().size2()) {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Inconsistent numer of integtaion pts %d != %d",
             toElePtr->gaussPts.size2(), getGaussPts().size2());
  }
#endif

  // TODO: add support for quad element types
  switch (getFEType()) {
  case MBTRI:
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Element type not implemented");
  }

  auto get_ftensor_from_mat_3d = [](MatrixDouble &m) {
    return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
        &m(0, 0), &m(0, 1), &m(0, 2));
  };

  // get local coordinates, i.e. local coordinates on child element using parent
  // local coordinates
  auto get_local_coords_triangle = [&]() {
    std::array<double, 3> ksi0 = {0, 1, 0};
    std::array<double, 3> ksi1 = {0, 0, 1};
    std::array<double, 9> ref_shapes;
    CHKERR Tools::shapeFunMBTRI<1>(ref_shapes.data(), ksi0.data(), ksi1.data(),
                                   3);
    auto &node_coords = getCoords();
    auto &glob_coords = toElePtr->coords;
    std::array<double, 6> local_coords;
    CHKERR Tools::getLocalCoordinatesOnReferenceThreeNodeTri(
        &*node_coords.begin(), &*glob_coords.begin(), 3, local_coords.data());
    return local_coords;
  };

  // get derivative of shape functions
  auto get_diff_triangle = [&]() {
    auto diff_ptr = Tools::diffShapeFunMBTRI.data();
    return FTensor::Tensor1<FTensor::PackPtr<const double *, 2>, 2>(
        diff_ptr, &diff_ptr[1]);
  };

  // get jacobian to map local coordinates of parent to child element
  auto get_jac = [&](auto &&local_coords, auto &&t_diff) {
    FTensor::Index<'I', 2> I;
    FTensor::Index<'J', 2> J;
    FTensor::Tensor2<double, 2, 2> t_jac;
    auto t_local_coords = getFTensor1FromPtr<2>(local_coords.data());
    t_jac(I, J) = 0;
    for (int nn = 0; nn != 3; ++nn) {
      t_jac(I, J) += t_local_coords(I) * t_diff(J);
      ++t_local_coords;
      ++t_diff;
    }
    return t_jac;
  };

  // get tangent vectors tensor
  auto t_mat_tangent = [&](auto &t1, auto &t2) {
    return FTensor::Tensor2<FTensor::PackPtr<double *, 3>, 2, 3>{
        &t1(0), &t1(1), &t1(2), &t2(0), &t2(1), &t2(2)};
  };

  // transform tangent vectors to child element tangents
  auto transform = [&](auto &&t_mat_t, auto &&t_mat_out_t, auto &&t_inv_jac) {
    FTensor::Index<'I', 2> I;
    FTensor::Index<'J', 2> J;
    FTensor::Index<'i', 3> i;
    FTensor::Number<0> N0;
    FTensor::Number<1> N1;
    for (auto gg = 0; gg != getGaussPts().size2(); ++gg) {
      FTensor::Tensor2<double, 2, 3> t_tmp;
      t_mat_out_t(J, i) = t_mat_t(I, i) * t_inv_jac(I, J);
      ++t_mat_t;
      ++t_mat_out_t;
    }
  };

  // calculate normal vector on child element
  auto calc_normal = [&](auto &n, auto &t1, auto &t2) {
    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;
    FTensor::Index<'k', 3> k;
    auto t_t1 = get_ftensor_from_mat_3d(t1);
    auto t_t2 = get_ftensor_from_mat_3d(t2);
    auto t_n = get_ftensor_from_mat_3d(n);
    for (auto gg = 0; gg != getGaussPts().size2(); ++gg) {
      t_n(j) = FTensor::levi_civita(i, j, k) * t_t1(k) * t_t2(i);
      ++t_t1;
      ++t_t2;
      ++t_n;
    }
  };

  transform(

      t_mat_tangent(getTangent1AtGaussPts(), getTangent2AtGaussPts()),
      t_mat_tangent(toElePtr->tangentOneAtGaussPts,
                    toElePtr->tangentTwoAtGaussPts),
      get_jac(get_local_coords_triangle(), get_diff_triangle())

  );
  calc_normal(toElePtr->normalsAtGaussPts, toElePtr->tangentOneAtGaussPts,
              toElePtr->tangentTwoAtGaussPts);

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
