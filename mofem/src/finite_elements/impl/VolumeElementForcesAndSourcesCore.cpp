/** \file VolumeElementForcesAndSourcesCore.cpp

\brief Implementation of volume element

*/

namespace MoFEM {

VolumeElementForcesAndSourcesCore::VolumeElementForcesAndSourcesCore(
    Interface &m_field)
    : ForcesAndSourcesCore(m_field), vOlume(elementMeasure),
      meshPositionsFieldName("MESH_NODE_POSITIONS"), coords(24), jAc(3, 3),
      invJac(3, 3), opSetInvJacH1(invJac),
      opContravariantPiolaTransform(elementMeasure, jAc),
      opCovariantPiolaTransform(invJac), opSetInvJacHdivAndHcurl(invJac),
      tJac(&jAc(0, 0), &jAc(0, 1), &jAc(0, 2), &jAc(1, 0), &jAc(1, 1),
           &jAc(1, 2), &jAc(2, 0), &jAc(2, 1), &jAc(2, 2)),
      tInvJac(&invJac(0, 0), &invJac(0, 1), &invJac(0, 2), &invJac(1, 0),
              &invJac(1, 1), &invJac(1, 2), &invJac(2, 0), &invJac(2, 1),
              &invJac(2, 2)) {}

VolumeElementForcesAndSourcesCore::VolumeElementForcesAndSourcesCore(
    Interface &m_field, const EntityType type)
    : VolumeElementForcesAndSourcesCore(m_field) {}

MoFEMErrorCode VolumeElementForcesAndSourcesCore::setIntegrationPts() {
  MoFEMFunctionBegin;

  int order_data = getMaxDataOrder();
  int order_row = getMaxRowOrder();
  int order_col = getMaxColOrder();
  const auto type = numeredEntFiniteElementPtr->getEntType();

  auto get_rule_by_type = [&]() {
    switch (type) {
    case MBHEX:
      return getRule(order_row + 1, order_col + 1, order_data + 1);
    default:
      return getRule(order_row, order_col, order_data);
    }
  };

  const int rule = get_rule_by_type();

  auto calc_base_for_tet = [&]() {
    MoFEMFunctionBegin;
    const size_t nb_gauss_pts = gaussPts.size2();
    auto &base = dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE);
    auto &diff_base = dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE);
    base.resize(nb_gauss_pts, 4, false);
    diff_base.resize(nb_gauss_pts, 12, false);
    CHKERR Tools::shapeFunMBTET(&*base.data().begin(), &gaussPts(0, 0),
                                &gaussPts(1, 0), &gaussPts(2, 0), nb_gauss_pts);
    double *diff_shape_ptr = &*diff_base.data().begin();
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      for (int nn = 0; nn != 4; ++nn) {
        for (int dd = 0; dd != 3; ++dd, ++diff_shape_ptr) {
          *diff_shape_ptr = Tools::diffShapeFunMBTET[3 * nn + dd];
        }
      }
    }
    MoFEMFunctionReturn(0);
  };

  auto calc_base_for_hex = [&]() {
    MoFEMFunctionBegin;
    const size_t nb_gauss_pts = gaussPts.size2();
    auto &base = dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE);
    auto &diff_base = dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE);
    base.resize(nb_gauss_pts, 8, false);
    diff_base.resize(nb_gauss_pts, 24, false);
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      const double ksi = gaussPts(0, gg);
      const double zeta = gaussPts(1, gg);
      const double eta = gaussPts(2, gg);
      base(gg, 0) = N_MBHEX0(ksi, zeta, eta);
      base(gg, 1) = N_MBHEX1(ksi, zeta, eta);
      base(gg, 2) = N_MBHEX2(ksi, zeta, eta);
      base(gg, 3) = N_MBHEX3(ksi, zeta, eta);
      base(gg, 4) = N_MBHEX4(ksi, zeta, eta);
      base(gg, 5) = N_MBHEX5(ksi, zeta, eta);
      base(gg, 6) = N_MBHEX6(ksi, zeta, eta);
      base(gg, 7) = N_MBHEX7(ksi, zeta, eta);
      diff_base(gg, 0 * 3 + 0) = diffN_MBHEX0x(zeta, eta);
      diff_base(gg, 1 * 3 + 0) = diffN_MBHEX1x(zeta, eta);
      diff_base(gg, 2 * 3 + 0) = diffN_MBHEX2x(zeta, eta);
      diff_base(gg, 3 * 3 + 0) = diffN_MBHEX3x(zeta, eta);
      diff_base(gg, 4 * 3 + 0) = diffN_MBHEX4x(zeta, eta);
      diff_base(gg, 5 * 3 + 0) = diffN_MBHEX5x(zeta, eta);
      diff_base(gg, 6 * 3 + 0) = diffN_MBHEX6x(zeta, eta);
      diff_base(gg, 7 * 3 + 0) = diffN_MBHEX7x(zeta, eta);
      diff_base(gg, 0 * 3 + 1) = diffN_MBHEX0y(ksi, eta);
      diff_base(gg, 1 * 3 + 1) = diffN_MBHEX1y(ksi, eta);
      diff_base(gg, 2 * 3 + 1) = diffN_MBHEX2y(ksi, eta);
      diff_base(gg, 3 * 3 + 1) = diffN_MBHEX3y(ksi, eta);
      diff_base(gg, 4 * 3 + 1) = diffN_MBHEX4y(ksi, eta);
      diff_base(gg, 5 * 3 + 1) = diffN_MBHEX5y(ksi, eta);
      diff_base(gg, 6 * 3 + 1) = diffN_MBHEX6y(ksi, eta);
      diff_base(gg, 7 * 3 + 1) = diffN_MBHEX7y(ksi, eta);
      diff_base(gg, 0 * 3 + 2) = diffN_MBHEX0z(ksi, zeta);
      diff_base(gg, 1 * 3 + 2) = diffN_MBHEX1z(ksi, zeta);
      diff_base(gg, 2 * 3 + 2) = diffN_MBHEX2z(ksi, zeta);
      diff_base(gg, 3 * 3 + 2) = diffN_MBHEX3z(ksi, zeta);
      diff_base(gg, 4 * 3 + 2) = diffN_MBHEX4z(ksi, zeta);
      diff_base(gg, 5 * 3 + 2) = diffN_MBHEX5z(ksi, zeta);
      diff_base(gg, 6 * 3 + 2) = diffN_MBHEX6z(ksi, zeta);
      diff_base(gg, 7 * 3 + 2) = diffN_MBHEX7z(ksi, zeta);
    }
    MoFEMFunctionReturn(0);
  };

  auto set_integration_pts_for_hex = [&]() {
    MoFEMFunctionBegin;
    CHKERR Tools::outerProductOfEdgeIntegrationPtsForHex(gaussPts, rule, rule,
                                                         rule);
    MoFEMFunctionReturn(0);
  };

  auto set_integration_pts_for_tet = [&]() {
    MoFEMFunctionBegin;
    if (rule < QUAD_3D_TABLE_SIZE) {
      if (QUAD_3D_TABLE[rule]->dim != 3) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong dimension");
      }
      if (QUAD_3D_TABLE[rule]->order < rule) {
        SETERRQ2(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                 "wrong order %d != %d", QUAD_3D_TABLE[rule]->order, rule);
      }
      size_t nb_gauss_pts = QUAD_3D_TABLE[rule]->npoints;
      gaussPts.resize(4, nb_gauss_pts, false);
      cblas_dcopy(nb_gauss_pts, &QUAD_3D_TABLE[rule]->points[1], 4,
                  &gaussPts(0, 0), 1);
      cblas_dcopy(nb_gauss_pts, &QUAD_3D_TABLE[rule]->points[2], 4,
                  &gaussPts(1, 0), 1);
      cblas_dcopy(nb_gauss_pts, &QUAD_3D_TABLE[rule]->points[3], 4,
                  &gaussPts(2, 0), 1);
      cblas_dcopy(nb_gauss_pts, QUAD_3D_TABLE[rule]->weights, 1,
                  &gaussPts(3, 0), 1);

      // CHKERR calc_base_for_tet();

      auto &base = dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE);
      auto &diff_base = dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE);
      base.resize(nb_gauss_pts, 4, false);
      diff_base.resize(nb_gauss_pts, 12, false);
      double *shape_ptr = &*base.data().begin();
      cblas_dcopy(4 * nb_gauss_pts, QUAD_3D_TABLE[rule]->points, 1, shape_ptr,
                  1);
      double *diff_shape_ptr = &*diff_base.data().begin();
      for (int gg = 0; gg != nb_gauss_pts; ++gg) {
        for (int nn = 0; nn != 4; ++nn) {
          for (int dd = 0; dd != 3; ++dd, ++diff_shape_ptr) {
            *diff_shape_ptr = Tools::diffShapeFunMBTET[3 * nn + dd];
          }
        }
      }

    } else {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "rule > quadrature order %d < %d", rule, QUAD_3D_TABLE_SIZE);
    }
    MoFEMFunctionReturn(0);
  };

  if (rule >= 0) {
    switch (type) {
    case MBTET:
      CHKERR set_integration_pts_for_tet();
      break;
    case MBHEX:
      CHKERR set_integration_pts_for_hex();
      CHKERR calc_base_for_hex();
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
      case MBTET: {
        CHKERR calc_base_for_tet();
      } break;
      case MBHEX:
        CHKERR calc_base_for_hex();
        break;
      default:
        SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                 "Element type not implemented: %d", type);
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VolumeElementForcesAndSourcesCore::calculateVolumeAndJacobian() {
  MoFEMFunctionBegin;
  const auto ent = numeredEntFiniteElementPtr->getEnt();
  const auto type = numeredEntFiniteElementPtr->getEntType();
  CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
  coords.resize(3 * num_nodes, false);
  CHKERR mField.get_moab().get_coords(conn, num_nodes, &*coords.data().begin());

  auto get_tet_t_diff_n = [&]() {
    return FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3>(
        &Tools::diffShapeFunMBTET[0], &Tools::diffShapeFunMBTET[1],
        &Tools::diffShapeFunMBTET[2]);
  };

  auto get_hex_t_diff_n = [&]() {
    return FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3>(
        &Tools::diffShapeFunMBHEXAtCenter[0],
        &Tools::diffShapeFunMBHEXAtCenter[1],
        &Tools::diffShapeFunMBHEXAtCenter[2]);
  };

  auto get_t_diff_n = [&]() {
    if (type == MBTET)
      return get_tet_t_diff_n();
    return get_hex_t_diff_n();
  };

  auto t_diff_n = get_t_diff_n();

  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_coords(
      &coords[0], &coords[1], &coords[2]);
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  jAc.clear();

  for (size_t n = 0; n != num_nodes; ++n) {
    tJac(i, j) += t_coords(i) * t_diff_n(j);
    ++t_coords;
    ++t_diff_n;
  }
  CHKERR determinantTensor3by3(tJac, vOlume);
  CHKERR invertTensor3by3(tJac, vOlume, tInvJac);

  if (type == MBTET)
    vOlume *= G_TET_W1[0] / 6.;

#ifndef NDEBUG
  if (vOlume <= std::numeric_limits<double>::epsilon() || vOlume != vOlume) {
    MOFEM_LOG("SELF", Sev::warning) << "Volume is zero " << vOlume;
    MOFEM_LOG("SELF", Sev::warning) << "Coords " << coords << endl;
  }
#endif

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VolumeElementForcesAndSourcesCore::calculateCoordinatesAtGaussPts() {
  MoFEMFunctionBegin;
  // Get coords at Gauss points
  FTensor::Index<'i', 3> i;

  auto &shape_functions = dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE);
  double *shape_functions_ptr = &*shape_functions.data().begin();
  const size_t nb_base_functions = shape_functions.size2();
  const size_t nb_gauss_pts = gaussPts.size2();
  coordsAtGaussPts.resize(nb_gauss_pts, 3, false);
  coordsAtGaussPts.clear();
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_coords_at_gauss_ptr(
      &coordsAtGaussPts(0, 0), &coordsAtGaussPts(0, 1),
      &coordsAtGaussPts(0, 2));
  FTensor::Tensor0<FTensor::PackPtr<double *, 1>> t_shape_functions(
      shape_functions_ptr);
  for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_coords(
        &coords[0], &coords[1], &coords[2]);
    for (int bb = 0; bb != nb_base_functions; ++bb) {
      t_coords_at_gauss_ptr(i) += t_coords(i) * t_shape_functions;
      ++t_coords;
      ++t_shape_functions;
    };
    ++t_coords_at_gauss_ptr;
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VolumeElementForcesAndSourcesCore::getSpaceBaseAndOrderOnElement() {
  MoFEMFunctionBegin;

  CHKERR getSpacesAndBaseOnEntities(dataH1);
  CHKERR getFaceNodes(dataH1);

  auto type = numeredEntFiniteElementPtr->getEntType();
  auto dim_type = CN::Dimension(type);

  auto get_data_on_ents = [&](auto lower_dim, auto space) {
    MoFEMFunctionBeginHot;
    auto data = dataOnElement[space];
    data->facesNodes = dataH1.facesNodes;
    data->facesNodesOrder = dataH1.facesNodesOrder;
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
  CHKERR get_data_on_ents(3, L2);    // L2

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VolumeElementForcesAndSourcesCore::transformBaseFunctions() {
  MoFEMFunctionBegin;

  if (numeredEntFiniteElementPtr->getEntType() == MBTET) {
    // OK that code is not nice.
    MatrixDouble new_diff_n;
    for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
      FTensor::Index<'i', 3> i;
      FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
      EntitiesFieldData::EntData &data = dataH1.dataOnEntities[MBVERTEX][0];
      if ((data.getDiffN(base).size1() == 4) &&
          (data.getDiffN(base).size2() == 3)) {
        const size_t nb_gauss_pts = gaussPts.size2();
        const size_t nb_base_functions = 4;
        new_diff_n.resize(nb_gauss_pts, 3 * nb_base_functions, false);
        double *new_diff_n_ptr = &*new_diff_n.data().begin();
        FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_new_diff_n(
            new_diff_n_ptr, &new_diff_n_ptr[1], &new_diff_n_ptr[2]);
        double *t_diff_n_ptr = &*data.getDiffN(base).data().begin();
        for (unsigned int gg = 0; gg != nb_gauss_pts; gg++) {
          FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_diff_n(
              t_diff_n_ptr, &t_diff_n_ptr[1], &t_diff_n_ptr[2]);
          for (unsigned int bb = 0; bb != nb_base_functions; bb++) {
            t_new_diff_n(i) = t_diff_n(i);
            ++t_new_diff_n;
            ++t_diff_n;
          }
        }
        data.getDiffN(base).resize(new_diff_n.size1(), new_diff_n.size2(),
                                   false);
        data.getDiffN(base).swap(new_diff_n);
      }
    }
  }

  if (dataH1.spacesOnEntities[MBVERTEX].test(H1))
    CHKERR opSetInvJacH1.opRhs(dataH1);

  std::array<std::bitset<LASTSPACE>, 4> spaces_by_dim{

      std::bitset<LASTSPACE>(0), //
      std::bitset<LASTSPACE>(0), //
      std::bitset<LASTSPACE>(0), //
      std::bitset<LASTSPACE>(0)

  };
  for (auto type = MBEDGE; type != MBENTITYSET; ++type) {
    spaces_by_dim[CN::Dimension(type)] |= dataH1.spacesOnEntities[type];
  }

  if (

      spaces_by_dim[1].test(HCURL) || //
      spaces_by_dim[2].test(HCURL) || //
      spaces_by_dim[3].test(HCURL)

  ) {
    CHKERR opCovariantPiolaTransform.opRhs(dataHcurl);
    CHKERR opSetInvJacHdivAndHcurl.opRhs(dataHcurl);
  }

  if (spaces_by_dim[2].test(HDIV) || spaces_by_dim[3].test(HDIV)) {
    CHKERR opContravariantPiolaTransform.opRhs(dataHdiv);
    // Fix for tetrahedrons
    if (numeredEntFiniteElementPtr->getEntType() == MBTET) {
      for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
        FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
        for (auto t : {MBTRI, MBTET}) {
          for (auto &d : dataHdiv.dataOnEntities[t]) {
            d.getN(base) /= 6;
            d.getDiffN(base) /= 6;
          }
        }
      }
    }
    CHKERR opSetInvJacHdivAndHcurl.opRhs(dataHdiv);
  }

  if (spaces_by_dim[3].test(L2)) {
    CHKERR opSetInvJacH1.opRhs(dataL2);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VolumeElementForcesAndSourcesCore::UserDataOperator::setPtrFE(
    ForcesAndSourcesCore *ptr) {
  MoFEMFunctionBeginHot;
  if (!(ptrFE = dynamic_cast<VolumeElementForcesAndSourcesCore *>(ptr)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "User operator and finite element do not work together");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode VolumeElementForcesAndSourcesCore::operator()() {
  MoFEMFunctionBegin;

  const auto type = numeredEntFiniteElementPtr->getEntType();
  if (type != lastEvaluatedElementEntityType) {
    switch (type) {
    case MBTET:
      getElementPolynomialBase() =
          boost::shared_ptr<BaseFunction>(new TetPolynomialBase(this));
      break;
    case MBHEX:
      getElementPolynomialBase() =
          boost::shared_ptr<BaseFunction>(new HexPolynomialBase());
      break;
    default:
      MoFEMFunctionReturnHot(0);
    }
    CHKERR createDataOnElement(type);
    lastEvaluatedElementEntityType = type;
  }

  CHKERR calculateVolumeAndJacobian();
  CHKERR getSpaceBaseAndOrderOnElement();
  CHKERR setIntegrationPts();
  if (gaussPts.size2() == 0)
    MoFEMFunctionReturnHot(0);
  CHKERR calculateCoordinatesAtGaussPts();
  CHKERR calHierarchicalBaseFunctionsOnElement();
  CHKERR calBernsteinBezierBaseFunctionsOnElement();

  CHKERR transformBaseFunctions();

  // Iterate over operators
  CHKERR loopOverOperators();

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
