/** \file FlatPrismElementForcesAndSourcesCore.cpp

\brief Implementation of flat prism element

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

FlatPrismElementForcesAndSourcesCore::FlatPrismElementForcesAndSourcesCore(
    Interface &m_field)
    : ForcesAndSourcesCore(m_field),
      meshPositionsFieldName("MESH_NODE_POSITIONS"),
      opHOCoordsAndNormals(hoCoordsAtGaussPtsF3, nOrmals_at_GaussPtF3,
                           tAngent1_at_GaussPtF3, tAngent2_at_GaussPtF3,
                           hoCoordsAtGaussPtsF4, nOrmals_at_GaussPtF4,
                           tAngent1_at_GaussPtF4, tAngent2_at_GaussPtF4) {
  getElementPolynomialBase() =
      boost::shared_ptr<BaseFunction>(new FlatPrismPolynomialBase());
}

MoFEMErrorCode FlatPrismElementForcesAndSourcesCore::operator()() {
  MoFEMFunctionBegin;

  if (numeredEntFiniteElementPtr->getEntType() != MBPRISM)
    MoFEMFunctionReturnHot(0);
  CHKERR createDataOnElement();

  DataForcesAndSourcesCore &data_div = *dataOnElement[HDIV];
  DataForcesAndSourcesCore &data_curl = *dataOnElement[HCURL];
  DataForcesAndSourcesCore &data_l2 = *dataOnElement[HCURL];

  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  int num_nodes;
  const EntityHandle *conn;
  CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
  {
    coords.resize(num_nodes * 3, false);
    CHKERR mField.get_moab().get_coords(conn, num_nodes,
                                        &*coords.data().begin());

    double diff_n[6];
    CHKERR ShapeDiffMBTRI(diff_n);
    normal.resize(6, false);
    CHKERR ShapeFaceNormalMBTRI(diff_n, &coords[0], &normal[0]);
    CHKERR ShapeFaceNormalMBTRI(diff_n, &coords[9], &normal[3]);
    aRea[0] = cblas_dnrm2(3, &normal[0], 1) * 0.5;
    aRea[1] = cblas_dnrm2(3, &normal[3], 1) * 0.5;
  }

  CHKERR getSpacesAndBaseOnEntities(dataH1);

  // H1
  if ((dataH1.spacesOnEntities[MBVERTEX]).test(H1)) {
    CHKERR getEntitySense<MBEDGE>(dataH1);
    CHKERR getEntityDataOrder<MBEDGE>(dataH1, H1);
    CHKERR getEntitySense<MBTRI>(dataH1);
    CHKERR getEntityDataOrder<MBTRI>(dataH1, H1);
  }

  // H1
  if ((dataH1.spacesOnEntities[MBEDGE]).test(HCURL)) {
    CHKERR getEntitySense<MBEDGE>(data_curl);
    CHKERR getEntityDataOrder<MBEDGE>(data_curl, HCURL);
    CHKERR getEntitySense<MBTRI>(data_curl);
    CHKERR getEntityDataOrder<MBTRI>(data_curl, HCURL);
  }

  // Hdiv
  if ((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
    CHKERR getEntitySense<MBTRI>(data_div);
    CHKERR getEntityDataOrder<MBTRI>(data_div, HDIV);
  }

  // L2
  if ((dataH1.spacesOnEntities[MBTRI]).test(L2)) {
    CHKERR getEntitySense<MBTRI>(data_l2);
    CHKERR getEntityDataOrder<MBTRI>(data_l2, L2);
  }

  int order_data = getMaxDataOrder();
  int order_row = getMaxRowOrder();
  int order_col = getMaxColOrder();
  int rule = getRule(order_row, order_col, order_data);
  int nb_gauss_pts;
  // int rule = getRule(order);
  if (rule >= 0) {
    if (rule < QUAD_2D_TABLE_SIZE) {
      if (QUAD_2D_TABLE[rule]->dim != 2) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong dimension");
      }
      if (QUAD_2D_TABLE[rule]->order < rule) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong order %d != %d", QUAD_2D_TABLE[rule]->order, rule);
      }
      nb_gauss_pts = QUAD_2D_TABLE[rule]->npoints;
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
    } else {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "rule > quadrature order %d < %d", rule, QUAD_2D_TABLE_SIZE);
      nb_gauss_pts = 0;
    }
  } else {
    CHKERR setGaussPts(order_row, order_col, order_data);
    nb_gauss_pts = gaussPts.size2();
    dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 3,
                                                           false);
    if (nb_gauss_pts) {
      CHKERR ShapeMBTRI(
          &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
          &gaussPts(0, 0), &gaussPts(1, 0), nb_gauss_pts);
    }
  }
  if (nb_gauss_pts == 0)
    MoFEMFunctionReturnHot(0);

  {
    coordsAtGaussPts.resize(nb_gauss_pts, 6, false);
    for (int gg = 0; gg < nb_gauss_pts; gg++) {
      for (int dd = 0; dd < 3; dd++) {
        coordsAtGaussPts(gg, dd) = cblas_ddot(
            3, &dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE)(gg, 0), 1,
            &coords[dd], 3);
        coordsAtGaussPts(gg, 3 + dd) = cblas_ddot(
            3, &dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE)(gg, 0), 1,
            &coords[9 + dd], 3);
      }
    }
  }

  for (int space = HCURL; space != LASTSPACE; ++space)
    if (dataOnElement[space]) {
      dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) =
          dataOnElement[H1]->dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
    }

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
    if (dataH1.bAse.test(b)) {
      switch (static_cast<FieldApproximationBase>(b)) {
      case AINSWORTH_LEGENDRE_BASE:
      case AINSWORTH_LOBATTO_BASE:
        if (dataH1.spacesOnEntities[MBVERTEX].test(H1)) {
          CHKERR getElementPolynomialBase()->getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(new FlatPrismPolynomialBaseCtx(
                  dataH1, mField.get_moab(), numeredEntFiniteElementPtr.get(),
                  H1, static_cast<FieldApproximationBase>(b), NOBASE)));
        }
        if (dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Not yet implemented");
        }
        if (dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Not yet implemented");
        }
        if (dataH1.spacesOnEntities[MBTET].test(L2)) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Not yet implemented");
        }
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Not yet implemented");
      }
    }
  }

  auto check_field = [&]() {
    auto field_it =
        fieldsPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName);
    if (field_it != fieldsPtr->get<FieldName_mi_tag>().end())
      if ((numeredEntFiniteElementPtr->getBitFieldIdData() &
           (*field_it)->getId())
              .any())
        return true;
    return false;
  };

  // Check if field meshPositionsFieldName exist
  if (check_field()) {

    hoCoordsAtGaussPtsF3.resize(nb_gauss_pts, 3, false);
    nOrmals_at_GaussPtF3.resize(nb_gauss_pts, 3, false);
    tAngent1_at_GaussPtF3.resize(nb_gauss_pts, 3, false);
    tAngent2_at_GaussPtF3.resize(nb_gauss_pts, 3, false);
    hoCoordsAtGaussPtsF4.resize(nb_gauss_pts, 3, false);
    nOrmals_at_GaussPtF4.resize(nb_gauss_pts, 3, false);
    tAngent1_at_GaussPtF4.resize(nb_gauss_pts, 3, false);
    tAngent2_at_GaussPtF4.resize(nb_gauss_pts, 3, false);
    CHKERR getNodesFieldData(dataH1, meshPositionsFieldName);
    CHKERR getEntityFieldData(dataH1, meshPositionsFieldName, MBEDGE);
    CHKERR getEntityFieldData(dataH1, meshPositionsFieldName, MBEDGE);
    CHKERR opHOCoordsAndNormals.opRhs(dataH1);
    CHKERR opHOCoordsAndNormals.calculateNormals();
  } else {
    hoCoordsAtGaussPtsF3.resize(0, 0, false);
    nOrmals_at_GaussPtF3.resize(0, 0, false);
    tAngent1_at_GaussPtF3.resize(0, 0, false);
    tAngent2_at_GaussPtF3.resize(0, 0, false);
    hoCoordsAtGaussPtsF4.resize(0, 0, false);
    nOrmals_at_GaussPtF4.resize(0, 0, false);
    tAngent1_at_GaussPtF4.resize(0, 0, false);
    tAngent2_at_GaussPtF4.resize(0, 0, false);
  }

  if (dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented yet");
  }

  // Iterate over operators
  CHKERR loopOverOperators();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FlatPrismElementForcesAndSourcesCore::UserDataOperator::setPtrFE(
    ForcesAndSourcesCore *ptr) {
  MoFEMFunctionBeginHot;
  if (!(ptrFE = dynamic_cast<FlatPrismElementForcesAndSourcesCore *>(ptr)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "User operator and finite element do not work together");
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM
