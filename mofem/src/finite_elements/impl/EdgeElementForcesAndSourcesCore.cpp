/** \file EdgeElementForcesAndSourcesCore.cpp

\brief Implementation of edge element

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

#include <BaseFunction.hpp>
#include <LegendrePolynomial.hpp>
#include <LobattoPolynomial.hpp>

#include <FTensor.hpp>
#include <DataStructures.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <EdgePolynomialBase.hpp> // Base functions on tet
#include <DataOperators.hpp>
#include <ForcesAndSourcesCore.hpp>
#include <EdgeElementForcesAndSourcesCore.hpp>

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

EdgeElementForcesAndSourcesCore::EdgeElementForcesAndSourcesCore(
    Interface &m_field)
    : ForcesAndSourcesCore(m_field),
      meshPositionsFieldName("MESH_NODE_POSITIONS"),
      opGetHoTangentOnEdge(tangentAtGaussPts),
      opCovariantTransform(dIrection, tangentAtGaussPts) {
  getElementPolynomialBase() =
      boost::shared_ptr<BaseFunction>(new EdgePolynomialBase());
}

MoFEMErrorCode EdgeElementForcesAndSourcesCore::calculateEdgeDirection() {
  MoFEMFunctionBegin;
  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  CHKERR mField.get_moab().get_connectivity(ent, cOnn, numNodes, true);
  cOords.resize(numNodes * 3, false);
  CHKERR mField.get_moab().get_coords(cOnn, numNodes, &*cOords.data().begin());
  dIrection.resize(3, false);
  cblas_dcopy(3, &cOords[3], 1, &*dIrection.data().begin(), 1);
  cblas_daxpy(3, -1., &cOords[0], 1, &*dIrection.data().begin(), 1);
  lEngth = cblas_dnrm2(3, &*dIrection.data().begin(), 1);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode EdgeElementForcesAndSourcesCore::setIntegrationPts() {
  MoFEMFunctionBegin;
  int order_data = getMaxDataOrder();
  int order_row = getMaxRowOrder();
  int order_col = getMaxColOrder();
  int rule = getRule(order_row, order_col, order_data);
  
  int nb_gauss_pts;
  if (rule >= 0) {
    if (rule < QUAD_1D_TABLE_SIZE) {
      if (QUAD_1D_TABLE[rule]->dim != 1) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong dimension");
      }
      if (QUAD_1D_TABLE[rule]->order < rule) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong order %d != %d", QUAD_1D_TABLE[rule]->order, rule);
      }
      nb_gauss_pts = QUAD_1D_TABLE[rule]->npoints;
      gaussPts.resize(2, nb_gauss_pts, false);
      cblas_dcopy(nb_gauss_pts, &QUAD_1D_TABLE[rule]->points[1], 2,
                  &gaussPts(0, 0), 1);
      cblas_dcopy(nb_gauss_pts, QUAD_1D_TABLE[rule]->weights, 1,
                  &gaussPts(1, 0), 1);
      dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 2,
                                                             false);
      double *shape_ptr =
          &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      cblas_dcopy(2 * nb_gauss_pts, QUAD_1D_TABLE[rule]->points, 1, shape_ptr,
                  1);
    } else {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "rule > quadrature order %d < %d", rule, QUAD_1D_TABLE_SIZE);
      nb_gauss_pts = 0;
    }
  } else {
    // If rule is negative, set user defined integration points
    CHKERR setGaussPts(order_row, order_col, order_data);
    nb_gauss_pts = gaussPts.size2();
    dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 2,
                                                           false);
    if (nb_gauss_pts) {
      for (int gg = 0; gg != nb_gauss_pts; gg++) {
        dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE)(gg, 0) =
            N_MBEDGE0(gaussPts(0, gg));
        dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE)(gg, 1) =
            N_MBEDGE1(gaussPts(0, gg));
      }
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
EdgeElementForcesAndSourcesCore::calculateCoordsAtIntegrationPts() {
  MoFEMFunctionBeginHot;
  const int nb_gauss_pts = gaussPts.size2();
  coordsAtGaussPts.resize(nb_gauss_pts, 3, false);
  if (nb_gauss_pts) {
    FTensor::Tensor0<FTensor::PackPtr<double *, 1> > t_ksi(
        &*gaussPts.data().begin());
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_coords(
        &coordsAtGaussPts(0, 0), &coordsAtGaussPts(0, 1),
        &coordsAtGaussPts(0, 2));
    FTensor::Tensor1<double, 3> t_coords_node0(cOords[0], cOords[1], cOords[2]);
    FTensor::Tensor1<double, 3> t_coords_node1(cOords[3], cOords[4], cOords[5]);
    FTensor::Index<'i', 3> i;
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      t_coords(i) = N_MBEDGE0(t_ksi) * t_coords_node0(i) +
                    N_MBEDGE1(t_ksi) * t_coords_node1(i);
      ++t_coords;
      ++t_ksi;
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
EdgeElementForcesAndSourcesCore::calculateHoCoordsAtIntegrationPts() {
  MoFEMFunctionBegin;
  
  if (dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName) !=
      dataPtr->get<FieldName_mi_tag>().end()) {
    CHKERR getNodesFieldData(dataH1, meshPositionsFieldName);
    CHKERR getEntityFieldData(dataH1, meshPositionsFieldName, MBEDGE);
    CHKERR opGetHoTangentOnEdge.opRhs(dataH1);
  } else {
    tangentAtGaussPts.resize(0, 3, false);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode EdgeElementForcesAndSourcesCore::operator()() {
  MoFEMFunctionBegin;

  if (numeredEntFiniteElementPtr->getEntType() != MBEDGE)
    MoFEMFunctionReturnHot(0);
 
  CHKERR createDataOnElement();
  
  DataForcesAndSourcesCore &data_curl = *dataOnElement[HCURL];

  CHKERR calculateEdgeDirection();
  CHKERR getSpacesAndBaseOnEntities(dataH1);
  CHKERR getEntityDataOrder<MBEDGE>(dataH1, H1);
  dataH1.dataOnEntities[MBEDGE][0].getSense() =
      1; // set sense to 1, this is this entity

  // Hcurl
  if (dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
    CHKERR getEntityDataOrder<MBEDGE>(data_curl, HCURL);
    data_curl.dataOnEntities[MBEDGE][0].getSense() =
        1; // set sense to 1, this is this entity
    data_curl.spacesOnEntities[MBEDGE].set(HCURL);
  }


  CHKERR setIntegrationPts();
  CHKERR calculateCoordsAtIntegrationPts();
  CHKERR calculateBaseFunctionsOnElement();
  CHKERR calculateHoCoordsAtIntegrationPts();

  if (dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
    CHKERR opCovariantTransform.opRhs(data_curl);
  }

  // Iterate over operators
  CHKERR loopOverOperators();

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
