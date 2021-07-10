/** \file VolumeElementForcesAndSourcesCore.cpp

\brief Implementation of volume element

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

VolumeElementForcesAndSourcesCoreBase::VolumeElementForcesAndSourcesCoreBase(
    Interface &m_field, const EntityType type)
    : ForcesAndSourcesCore(m_field), vOlume(elementMeasure),
      meshPositionsFieldName("MESH_NODE_POSITIONS"), coords(12), jAc(3, 3),
      invJac(3, 3), opSetInvJacH1(invJac),
      opContravariantPiolaTransform(elementMeasure, jAc),
      opCovariantPiolaTransform(invJac), opSetInvJacHdivAndHcurl(invJac),
      tJac(&jAc(0, 0), &jAc(0, 1), &jAc(0, 2), &jAc(1, 0), &jAc(1, 1),
           &jAc(1, 2), &jAc(2, 0), &jAc(2, 1), &jAc(2, 2)),
      tInvJac(&invJac(0, 0), &invJac(0, 1), &invJac(0, 2), &invJac(1, 0),
              &invJac(1, 1), &invJac(1, 2), &invJac(2, 0), &invJac(2, 1),
              &invJac(2, 2)) {
  getElementPolynomialBase() =
      boost::shared_ptr<BaseFunction>(new TetPolynomialBase());
}

MoFEMErrorCode VolumeElementForcesAndSourcesCoreBase::setIntegrationPts() {
  MoFEMFunctionBegin;
  int order_data = getMaxDataOrder();
  int order_row = getMaxRowOrder();
  int order_col = getMaxColOrder();
  int rule = getRule(order_row, order_col, order_data);

  if (rule >= 0) {
    if (rule < QUAD_3D_TABLE_SIZE) {
      if (QUAD_3D_TABLE[rule]->dim != 3) {
        SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY, "wrong dimension");
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
      cblas_dcopy(nb_gauss_pts, QUAD_3D_TABLE[rule]->weights, 1, &gaussPts(3, 0),
                  1);
      dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 4,
                                                             false);
      double *shape_ptr =
          &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      cblas_dcopy(4 * nb_gauss_pts, QUAD_3D_TABLE[rule]->points, 1, shape_ptr, 1);
    } else {
      SETERRQ2(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "rule > quadrature order %d < %d", rule, QUAD_3D_TABLE_SIZE);
    }
  } else {
    CHKERR setGaussPts(order_row, order_col, order_data);
    const size_t nb_gauss_pts = gaussPts.size2();
    dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 4,
                                                           false);
    if (nb_gauss_pts > 0) {
      CHKERR Tools::shapeFunMBTET(
          &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
          &gaussPts(0, 0), &gaussPts(1, 0), &gaussPts(2, 0), nb_gauss_pts);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VolumeElementForcesAndSourcesCoreBase::calculateVolumeAndJacobian() {
  MoFEMFunctionBegin;
  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
  CHKERR mField.get_moab().get_coords(conn, num_nodes, &*coords.data().begin());
  FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3> t_diff_n(
      &Tools::diffShapeFunMBTET[0], &Tools::diffShapeFunMBTET[1],
      &Tools::diffShapeFunMBTET[2]);
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_coords(
      &coords[0], &coords[1], &coords[2]);
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  jAc.clear();
  for (auto n : {0, 1, 2, 3}) {
    tJac(i, j) += t_coords(i) * t_diff_n(j);
    ++t_coords;
    ++t_diff_n;
  }
  CHKERR determinantTensor3by3(tJac, vOlume);
  CHKERR invertTensor3by3(tJac, vOlume, tInvJac);
  vOlume *= G_TET_W1[0] / 6.;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VolumeElementForcesAndSourcesCoreBase::calculateCoordinatesAtGaussPts() {
  MoFEMFunctionBegin;
  // Get coords at Gauss points
  FTensor::Index<'i', 3> i;

  double *shape_functions_ptr =
      &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
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
    for (int bb = 0; bb != 4; ++bb) {
      t_coords_at_gauss_ptr(i) += t_coords(i) * t_shape_functions;
      ++t_coords;
      ++t_shape_functions;
    };
    ++t_coords_at_gauss_ptr;
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VolumeElementForcesAndSourcesCoreBase::getSpaceBaseAndOrderOnElement() {
  MoFEMFunctionBegin;

  CHKERR getSpacesAndBaseOnEntities(dataH1);
  CHKERR getFaceTriNodes(dataH1);
  // H1
  if ((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
    CHKERR getEntitySense<MBEDGE>(dataH1);
    CHKERR getEntityDataOrder<MBEDGE>(dataH1, H1);
  }
  if ((dataH1.spacesOnEntities[MBTRI]).test(H1)) {
    CHKERR getEntitySense<MBTRI>(dataH1);
    CHKERR getEntityDataOrder<MBTRI>(dataH1, H1);
  }
  if ((dataH1.spacesOnEntities[MBTET]).test(H1)) {
    CHKERR getEntityDataOrder<MBTET>(dataH1, H1);
  }
  // Hcurl
  if ((dataH1.spacesOnEntities[MBEDGE]).test(HCURL)) {
    CHKERR getEntitySense<MBEDGE>(dataHcurl);
    CHKERR getEntityDataOrder<MBEDGE>(dataHcurl, HCURL);
    dataHcurl.spacesOnEntities[MBEDGE].set(HCURL);
  }
  if ((dataH1.spacesOnEntities[MBTRI]).test(HCURL)) {
    dataHcurl.facesNodes = dataH1.facesNodes;
    CHKERR getEntitySense<MBTRI>(dataHcurl);
    CHKERR getEntityDataOrder<MBTRI>(dataHcurl, HCURL);
    dataHcurl.spacesOnEntities[MBTRI].set(HCURL);
  }
  if ((dataH1.spacesOnEntities[MBTET]).test(HCURL)) {
    CHKERR getEntityDataOrder<MBTET>(dataHcurl, HCURL);
    dataHcurl.spacesOnEntities[MBTET].set(HCURL);
  }
  // Hdiv
  if ((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
    dataHdiv.facesNodes = dataH1.facesNodes;
    CHKERR getEntitySense<MBTRI>(dataHdiv);
    CHKERR getEntityDataOrder<MBTRI>(dataHdiv, HDIV);
    dataHdiv.spacesOnEntities[MBTRI].set(HDIV);
  }
  if ((dataH1.spacesOnEntities[MBTET]).test(HDIV)) {
    CHKERR getEntityDataOrder<MBTET>(dataHdiv, HDIV);
    dataHdiv.spacesOnEntities[MBTET].set(HDIV);
  }
  // L2
  if ((dataH1.spacesOnEntities[MBTET]).test(L2)) {
    CHKERR getEntityDataOrder<MBTET>(dataL2, L2);
    dataL2.spacesOnEntities[MBTET].set(L2);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VolumeElementForcesAndSourcesCoreBase::transformBaseFunctions() {
  MoFEMFunctionBegin;

  CHKERR opSetInvJacH1.opRhs(dataH1);
  if (dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
    CHKERR opCovariantPiolaTransform.opRhs(dataHcurl);
    CHKERR opSetInvJacHdivAndHcurl.opRhs(dataHcurl);
  }
  if (dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    CHKERR opContravariantPiolaTransform.opRhs(dataHdiv);
    CHKERR opSetInvJacHdivAndHcurl.opRhs(dataHdiv);
  }
  if (dataH1.spacesOnEntities[MBTET].test(L2)) {
    CHKERR opSetInvJacH1.opRhs(dataL2);
  }

  MatrixDouble new_diff_n;
  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
    FTensor::Index<'i', 3> i;
    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
    DataForcesAndSourcesCore::EntData &data =
        dataH1.dataOnEntities[MBVERTEX][0];
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
      data.getDiffN(base).resize(new_diff_n.size1(), new_diff_n.size2(), false);
      data.getDiffN(base).data().swap(new_diff_n.data());
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VolumeElementForcesAndSourcesCoreBase::UserDataOperator::setPtrFE(
    ForcesAndSourcesCore *ptr) {
  MoFEMFunctionBeginHot;
  if(!(ptrFE = dynamic_cast<VolumeElementForcesAndSourcesCoreBase *>(ptr))) 
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "User operator and finite element do not work together");
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM
