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

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <base_functions.h>
#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <FTensor.hpp>
#include <DataStructures.hpp>
#include <DataOperators.hpp>
#include <ForcesAndSourcesCore.hpp>
#include <BaseFunction.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <FlatPrismPolynomialBase.hpp>
#include <FlatPrismElementForcesAndSourcesCore.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  #include <cblas.h>
  #include <lapack_wrap.h>
  // #include <gm_rule.h>
  #include <quad.h>
#ifdef __cplusplus
}
#endif

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
          dataOnElement[H1]->dataOnEntities[MBVERTEX][0].getNSharedPtr(
              NOBASE);
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

  if (dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName) !=
      dataPtr->get<FieldName_mi_tag>().end()) {
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

  doVertices = true;
  doEdges = false;
  doQuads = false;
  doTris = false;
  doTets = false;
  doPrisms = false;

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
