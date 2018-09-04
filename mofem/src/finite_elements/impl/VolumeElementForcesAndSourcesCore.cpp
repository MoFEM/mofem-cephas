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
#include <TetPolynomialBase.hpp> // Base functions on tet
#include <DataOperators.hpp>
#include <ForcesAndSourcesCore.hpp>
#include <VolumeElementForcesAndSourcesCore.hpp>
#include <FaceElementForcesAndSourcesCore.hpp>

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

VolumeElementForcesAndSourcesCore::VolumeElementForcesAndSourcesCore(
    Interface &m_field, const EntityType type)
    : ForcesAndSourcesCore(m_field), coords(12), jAc(3, 3), invJac(3, 3),
      dataH1(type), derivedDataH1(dataH1), dataL2(type), derivedDataL2(dataL2),
      dataHdiv(type), derivedDataHdiv(dataHdiv), dataHcurl(type),
      derivedDataHcurl(dataHcurl), dataNoField(type), dataNoFieldCol(type),
      opSetInvJacH1(invJac), opContravariantPiolaTransform(vOlume, jAc),
      opCovariantPiolaTransform(invJac), opSetInvJacHdivAndHcurl(invJac),
      meshPositionsFieldName("MESH_NODE_POSITIONS"),
      opHOatGaussPoints(hoCoordsAtGaussPts, hoGaussPtsJac),
      opSetHoInvJacH1(hoGaussPtsInvJac),
      opHoContravariantTransform(hoGaussPtsDetJac, hoGaussPtsJac),
      opHoCovariantTransform(hoGaussPtsInvJac),
      opSetHoInvJacHdivAndHcurl(hoGaussPtsInvJac),
      tJac(&jAc(0, 0), &jAc(0, 1), &jAc(0, 2), &jAc(1, 0), &jAc(1, 1),
           &jAc(1, 2), &jAc(2, 0), &jAc(2, 1), &jAc(2, 2)),
      tInvJac(&invJac(0, 0), &invJac(0, 1), &invJac(0, 2), &invJac(1, 0),
              &invJac(1, 1), &invJac(1, 2), &invJac(2, 0), &invJac(2, 1),
              &invJac(2, 2)) {}

MoFEMErrorCode VolumeElementForcesAndSourcesCore::setIntegrationPts() {
  MoFEMFunctionBegin;
  int order_data = getMaxDataOrder();
  int order_row = getMaxRowOrder();
  int order_col = getMaxColOrder();
  int rule = getRule(order_row, order_col, order_data);
  // std::cerr << order_data << " " << order_row << " " << order_col << " " <<
  // rule << std::endl;
  if (rule >= 0) {
    if (rule < QUAD_3D_TABLE_SIZE) {
      if (QUAD_3D_TABLE[rule]->dim != 3) {
        SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY, "wrong dimension");
      }
      if (QUAD_3D_TABLE[rule]->order < rule) {
        SETERRQ2(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                 "wrong order %d != %d", QUAD_3D_TABLE[rule]->order, rule);
      }
      nbGaussPts = QUAD_3D_TABLE[rule]->npoints;
      gaussPts.resize(4, nbGaussPts, false);
      cblas_dcopy(nbGaussPts, &QUAD_3D_TABLE[rule]->points[1], 4,
                  &gaussPts(0, 0), 1);
      cblas_dcopy(nbGaussPts, &QUAD_3D_TABLE[rule]->points[2], 4,
                  &gaussPts(1, 0), 1);
      cblas_dcopy(nbGaussPts, &QUAD_3D_TABLE[rule]->points[3], 4,
                  &gaussPts(2, 0), 1);
      cblas_dcopy(nbGaussPts, QUAD_3D_TABLE[rule]->weights, 1, &gaussPts(3, 0),
                  1);
      dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nbGaussPts, 4,
                                                             false);
      double *shape_ptr =
          &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      cblas_dcopy(4 * nbGaussPts, QUAD_3D_TABLE[rule]->points, 1, shape_ptr, 1);
    } else {
      SETERRQ2(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "rule > quadrature order %d < %d", rule, QUAD_3D_TABLE_SIZE);
      nbGaussPts = 0;
    }
  } else {
    CHKERR setGaussPts(order_row, order_col, order_data);
    nbGaussPts = gaussPts.size2();
    dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nbGaussPts, 4,
                                                           false);
    if (nbGaussPts > 0) {
      CHKERR ShapeMBTET(
          &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
          &gaussPts(0, 0), &gaussPts(1, 0), &gaussPts(2, 0), nbGaussPts);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VolumeElementForcesAndSourcesCore::calculateVolumeAndJacobian() {
  MoFEMFunctionBegin;
  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
  CHKERR mField.get_moab().get_coords(conn, num_nodes, &*coords.data().begin());
  double diff_n[12];
  CHKERR ShapeDiffMBTET(diff_n);
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_diff_n(
      &diff_n[0], &diff_n[1], &diff_n[2]);
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_coords(
      &coords[0], &coords[1], &coords[2]);
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  jAc.clear();
  for (int nn = 0; nn != 4; nn++) {
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
VolumeElementForcesAndSourcesCore::calculateCoordinatesAtGaussPts() {
  MoFEMFunctionBegin;
  // Get coords at Gauss points
  FTensor::Index<'i', 3> i;
  double *shape_functions_ptr =
      &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
  coordsAtGaussPts.resize(nbGaussPts, 3, false);
  coordsAtGaussPts.clear();
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_coords_at_gauss_ptr(
      &coordsAtGaussPts(0, 0), &coordsAtGaussPts(0, 1),
      &coordsAtGaussPts(0, 2));
  FTensor::Tensor0<FTensor::PackPtr<double *, 1>> t_shape_functions(
      shape_functions_ptr);
  for (unsigned int gg = 0; gg < nbGaussPts; gg++) {
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_coords(
        &coords[0], &coords[1], &coords[2]);
    for (int bb = 0; bb < 4; bb++) {
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
  CHKERR getFaceTriNodes(dataH1);
  // H1
  if ((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
    CHKERR getEdgesSense(dataH1);
    CHKERR getEdgesDataOrder(dataH1, H1);
  }
  if ((dataH1.spacesOnEntities[MBTRI]).test(H1)) {
    CHKERR getTrisSense(dataH1);
    CHKERR getTrisDataOrder(dataH1, H1);
  }
  if ((dataH1.spacesOnEntities[MBTET]).test(H1)) {
    CHKERR getTetDataOrder(dataH1, H1);
  }
  // Hcurl
  if ((dataH1.spacesOnEntities[MBEDGE]).test(HCURL)) {
    CHKERR getEdgesSense(dataHcurl);
    CHKERR getEdgesDataOrder(dataHcurl, HCURL);
    dataHcurl.spacesOnEntities[MBEDGE].set(HCURL);
  }
  if ((dataH1.spacesOnEntities[MBTRI]).test(HCURL)) {
    CHKERR getTrisSense(dataHcurl);
    CHKERR getFaceTriNodes(dataHcurl);
    CHKERR getTrisDataOrder(dataHcurl, HCURL);
    dataHcurl.spacesOnEntities[MBTRI].set(HCURL);
  }
  if ((dataH1.spacesOnEntities[MBTET]).test(HCURL)) {
    CHKERR getTetDataOrder(dataHcurl, HCURL);
    dataHcurl.spacesOnEntities[MBTET].set(HCURL);
  }
  // Hdiv
  if ((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
    CHKERR getTrisSense(dataHdiv);
    CHKERR getFaceTriNodes(dataHdiv);
    CHKERR getTrisDataOrder(dataHdiv, HDIV);
    dataHdiv.spacesOnEntities[MBTRI].set(HDIV);
  }
  if ((dataH1.spacesOnEntities[MBTET]).test(HDIV)) {
    CHKERR getTetDataOrder(dataHdiv, HDIV);
    dataHdiv.spacesOnEntities[MBTET].set(HDIV);
  }
  // L2
  if ((dataH1.spacesOnEntities[MBTET]).test(L2)) {
    CHKERR getTetDataOrder(dataL2, L2);
    dataL2.spacesOnEntities[MBTET].set(L2);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VolumeElementForcesAndSourcesCore::calculateBaseFunctionsOnElement(
    const int b) {
  MoFEMFunctionBegin;
  if (dataH1.bAse.test(b)) {
    switch (ApproximationBaseArray[b]) {
    case NOBASE:
      break;
    case AINSWORTH_LEGENDRE_BASE:
    case AINSWORTH_LOBATTO_BASE:
      if (dataH1.spacesOnEntities[MBTET].test(L2) &&
          dataH1.basesOnEntities[MBTET].test(b) &&
            dataH1.basesOnSpaces[L2].test(b)) {
        CHKERR TetPolynomialBase().getValue(
            gaussPts,
            boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                dataL2, L2, ApproximationBaseArray[b], NOBASE)));
      }
      if (dataH1.spacesOnEntities[MBVERTEX].test(H1) &&
          dataH1.basesOnEntities[MBVERTEX].test(b) &&
            dataH1.basesOnSpaces[H1].test(b)) {
        CHKERR TetPolynomialBase().getValue(
            gaussPts,
            boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                dataH1, H1, ApproximationBaseArray[b], NOBASE)));
      }
    case DEMKOWICZ_JACOBI_BASE:
      if (dataH1.spacesOnEntities[MBEDGE].test(HCURL) &&
          dataH1.basesOnEntities[MBEDGE].test(b) &&
            dataH1.basesOnSpaces[HCURL].test(b)) {
        CHKERR TetPolynomialBase().getValue(
            gaussPts,
            boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                dataHcurl, HCURL, ApproximationBaseArray[b], NOBASE)));
      }
      if (dataH1.spacesOnEntities[MBTRI].test(HDIV) &&
          dataH1.basesOnEntities[MBTRI].test(b) &&
            dataH1.basesOnSpaces[HDIV].test(b)) {
        CHKERR TetPolynomialBase().getValue(
            gaussPts,
            boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                dataHdiv, HDIV, ApproximationBaseArray[b], NOBASE)));
      }
      break;
    default:
      SETERRQ1(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "Base <%s> not yet implemented",
               ApproximationBaseNames[ApproximationBaseArray[b]]);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VolumeElementForcesAndSourcesCore::calculateBaseFunctionsOnElement() {
  MoFEMFunctionBegin;
  /// Use the some node base
  dataHdiv.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) =
      dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
  dataHcurl.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) =
      dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
  dataL2.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) =
      dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
  std::vector<FieldApproximationBase> shape_functions_for_bases;
  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
    CHKERR calculateBaseFunctionsOnElement(b);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VolumeElementForcesAndSourcesCore::transformBaseFunctions() {
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
  MoFEMFunctionReturn(0);
  }

MoFEMErrorCode VolumeElementForcesAndSourcesCore::calculateHoJacobian() {
  MoFEMFunctionBegin;
  if (dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName) !=
      dataPtr->get<FieldName_mi_tag>().end()) {
    const Field *field_struture =
        mField.get_field_structure(meshPositionsFieldName);
    BitFieldId id = field_struture->getId();
    if ((numeredEntFiniteElementPtr->getBitFieldIdData() & id).none()) {
      SETERRQ(mField.get_comm(), MOFEM_NOT_FOUND,
              "no MESH_NODE_POSITIONS in element data");
    }
    CHKERR getEdgesDataOrderSpaceAndBase(dataH1, meshPositionsFieldName);
    CHKERR getTrisDataOrderSpaceAndBase(dataH1, meshPositionsFieldName);
    CHKERR getTetDataOrderSpaceAndBase(dataH1, meshPositionsFieldName);
    CHKERR getNodesFieldData(dataH1, meshPositionsFieldName);
    if (dataH1.dataOnEntities[MBVERTEX][0].getFieldData().size() != 12) {
      SETERRQ(mField.get_comm(), MOFEM_NOT_FOUND,
              "no MESH_NODE_POSITIONS in element data or field has wrong "
              "number of coefficients");
    }
    CHKERR getEdgesFieldData(dataH1, meshPositionsFieldName);
    CHKERR getTrisFieldData(dataH1, meshPositionsFieldName);
    CHKERR getTetsFieldData(dataH1, meshPositionsFieldName);
    CHKERR opHOatGaussPoints.opRhs(dataH1);
      hoGaussPtsInvJac.resize(hoGaussPtsJac.size1(), hoGaussPtsJac.size2(),
                              false);
      // Express Jacobian as tensor
      FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> jac(
          &hoGaussPtsJac(0, 0), &hoGaussPtsJac(0, 1), &hoGaussPtsJac(0, 2),
          &hoGaussPtsJac(0, 3), &hoGaussPtsJac(0, 4), &hoGaussPtsJac(0, 5),
          &hoGaussPtsJac(0, 6), &hoGaussPtsJac(0, 7), &hoGaussPtsJac(0, 8));
      FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> inv_jac(
          &hoGaussPtsInvJac(0, 0), &hoGaussPtsInvJac(0, 1),
          &hoGaussPtsInvJac(0, 2), &hoGaussPtsInvJac(0, 3),
          &hoGaussPtsInvJac(0, 4), &hoGaussPtsInvJac(0, 5),
          &hoGaussPtsInvJac(0, 6), &hoGaussPtsInvJac(0, 7),
          &hoGaussPtsInvJac(0, 8));
      hoGaussPtsDetJac.resize(nbGaussPts, false);
      FTensor::Tensor0<double *> det(&hoGaussPtsDetJac[0]);
      // Calculate inverse and determinant
      for (unsigned int gg = 0; gg != nbGaussPts; gg++) {
      CHKERR determinantTensor3by3(jac, det);
        // if(det<0) {
        //   SETERRQ(mField.get_comm(),MOFEM_DATA_INCONSISTENCY,"Negative
        //   volume");
        // }
      CHKERR invertTensor3by3(jac, det, inv_jac);
        ++jac;
        ++inv_jac;
        ++det;
      }
  } else {
    hoCoordsAtGaussPts.resize(0, 0, false);
    hoGaussPtsInvJac.resize(0, 0, false);
    hoGaussPtsDetJac.resize(0, false);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VolumeElementForcesAndSourcesCore::transformHoBaseFunctions() {
  MoFEMFunctionBegin;

  if (hoCoordsAtGaussPts.size1() > 0) {
      // Transform derivatives of base functions and apply Piola transformation
      // if needed.
    CHKERR opSetHoInvJacH1.opRhs(dataH1);
      if (dataH1.spacesOnEntities[MBTET].test(L2)) {
      CHKERR opSetHoInvJacH1.opRhs(dataL2);
      }
      if (dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
      CHKERR opHoContravariantTransform.opRhs(dataHdiv);
      CHKERR opSetHoInvJacHdivAndHcurl.opRhs(dataHdiv);
      }
      if (dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
      CHKERR opHoCovariantTransform.opRhs(dataHcurl);
      CHKERR opSetHoInvJacHdivAndHcurl.opRhs(dataHcurl);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VolumeElementForcesAndSourcesCore::operator()() {
  MoFEMFunctionBegin;

  if (numeredEntFiniteElementPtr->getEntType() != MBTET)
    MoFEMFunctionReturnHot(0);

  CHKERR calculateVolumeAndJacobian();
  CHKERR getSpaceBaseAndOrderOnElement();
  CHKERR setIntegrationPts();
  if (nbGaussPts == 0)
    MoFEMFunctionReturnHot(0);
  CHKERR calculateCoordinatesAtGaussPts();
  CHKERR calculateBaseFunctionsOnElement();
  CHKERR transformBaseFunctions();

  try {
    MatrixDouble new_diff_n;
    for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
      FTensor::Index<'i', 3> i;
      FieldApproximationBase base = ApproximationBaseArray[b];
      DataForcesAndSourcesCore::EntData &data =
          dataH1.dataOnEntities[MBVERTEX][0];
      if ((data.getDiffN(base).size1() == 4) &&
          (data.getDiffN(base).size2() == 3)) {
        const unsigned int nb_base_functions = 4;
        new_diff_n.resize(nbGaussPts, 3 * nb_base_functions, false);
        double *new_diff_n_ptr = &*new_diff_n.data().begin();
        FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_new_diff_n(
            new_diff_n_ptr, &new_diff_n_ptr[1], &new_diff_n_ptr[2]);
        double *t_diff_n_ptr = &*data.getDiffN(base).data().begin();
        for (unsigned int gg = 0; gg != nbGaussPts; gg++) {
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
        data.getDiffN(base).data().swap(new_diff_n.data());
      }
    }
  } catch (std::exception &ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__
       << " in file " << __FILE__;
    SETERRQ(mField.get_comm(), MOFEM_STD_EXCEPTION_THROW, ss.str().c_str());
  }

  CHKERR calculateHoJacobian();
  CHKERR transformHoBaseFunctions();

  const UserDataOperator::OpType types[2] = {UserDataOperator::OPROW,
                                             UserDataOperator::OPCOL};
  std::vector<std::string> last_eval_field_name(2);
  DataForcesAndSourcesCore *op_data[2];
  FieldSpace space[2];
  FieldApproximationBase base[2];

  boost::ptr_vector<UserDataOperator>::iterator oit, hi_oit;
  oit = opPtrVector.begin();
  hi_oit = opPtrVector.end();

  for (; oit != hi_oit; oit++) {

    oit->setPtrFE(this);

    if (oit->opType == UserDataOperator::OPLAST) {

      // Set field
      switch (oit->sPace) {
      case NOSPACE:
        SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY, "unknown space");
      case H1 :
        op_data[0] = &dataH1;
        break;
      case HCURL:
        op_data[0] = &dataHcurl;
        break;
      case HDIV:
        op_data[0] = &dataHdiv;
        break;
      case L2:
        op_data[0] = &dataL2;
        break;
      case NOFIELD:
        op_data[0] = &dataNoField;
        break;
      case LASTSPACE:
        SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY, "unknown space");
        break;
      }

      // Reseat all data which all field dependent
      op_data[0]->resetFieldDependentData();
      last_eval_field_name[0] = "";
      last_eval_field_name[1] = "";

      // Run operator
      CHKERR oit->opRhs(*op_data[0], oit->doVertices, oit->doEdges,
                        oit->doQuads, oit->doTris, oit->doTets, false);

    } else {

      for (int ss = 0; ss != 2; ss++) {

        std::string field_name = !ss ? oit->rowFieldName : oit->colFieldName;
        const Field *field_struture = mField.get_field_structure(field_name);
        BitFieldId data_id = field_struture->getId();

        if ((oit->getNumeredEntFiniteElementPtr()->getBitFieldIdData() &
             data_id)
                .none()) {
          SETERRQ2(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                   "no data field < %s > on finite element < %s >",
                   field_name.c_str(), feName.c_str());
        }

        if (oit->getOpType() & types[ss] ||
            oit->getOpType() & UserDataOperator::OPROWCOL) {

          space[ss] = field_struture->getSpace();
          switch (space[ss]) {
          case NOSPACE:
            SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "unknown space");
            break;
          case H1:
            op_data[ss] = !ss ? &dataH1 : &derivedDataH1;
            break;
          case HCURL:
            op_data[ss] = !ss ? &dataHcurl : &derivedDataHcurl;
            break;
          case HDIV:
            op_data[ss] = !ss ? &dataHdiv : &derivedDataHdiv;
            break;
          case L2:
            op_data[ss] = !ss ? &dataL2 : &derivedDataL2;
            break;
          case NOFIELD:
            op_data[ss] = !ss ? &dataNoField : &dataNoFieldCol;
            break;
          case LASTSPACE:
            SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "unknown space");
            break;
          }

          base[ss] = field_struture->getApproxBase();
          switch (base[ss]) {
          case NOBASE:
          case AINSWORTH_LEGENDRE_BASE:
          case AINSWORTH_LOBATTO_BASE:
          case DEMKOWICZ_JACOBI_BASE:
            break;
          default:
            SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "unknown or not implemented base");
            break;
          }

          if (last_eval_field_name[ss] != field_name) {

            switch (space[ss]) {
            case NOSPACE:
              SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                      "unknown space");
              break;
            case H1:
              if (!ss) {
                CHKERR getRowNodesIndices(*op_data[ss], field_name);
              } else {
                CHKERR getColNodesIndices(*op_data[ss], field_name);
              }
              CHKERR getNodesFieldData(*op_data[ss], field_name);
            case HCURL:
              if (!ss) {
                CHKERR getEdgesRowIndices(*op_data[ss], field_name);
              } else {
                CHKERR getEdgesColIndices(*op_data[ss], field_name);
              }
              CHKERR getEdgesDataOrderSpaceAndBase(*op_data[ss], field_name);
              CHKERR getEdgesFieldData(*op_data[ss], field_name);
            case HDIV:
              if (!ss) {
                CHKERR getTrisRowIndices(*op_data[ss], field_name);
              } else {
                CHKERR getTrisColIndices(*op_data[ss], field_name);
              }
              CHKERR getTrisDataOrderSpaceAndBase(*op_data[ss], field_name);
              CHKERR getTrisFieldData(*op_data[ss], field_name);
            case L2:
              if (!ss) {
                CHKERR getTetsRowIndices(*op_data[ss], field_name);
              } else {
                CHKERR getTetsColIndices(*op_data[ss], field_name);
              }
              CHKERR getTetDataOrderSpaceAndBase(*op_data[ss], field_name);
              CHKERR getTetsFieldData(*op_data[ss], field_name);
              break;
            case NOFIELD:
              if (!getNinTheLoop()) {
                // NOFIELD data are the same for each element, can be
                // retrieved only once
                if (!ss) {
                  CHKERR getNoFieldRowIndices(*op_data[ss], field_name);
                } else {
                  CHKERR getNoFieldColIndices(*op_data[ss], field_name);
                }
                CHKERR getNoFieldFieldData(*op_data[ss], field_name);
              }
              break;
            case LASTSPACE:
              SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                      "unknown space");
              break;
            }
            last_eval_field_name[ss] = field_name;
          }
        }
      }

      if (oit->getOpType() & UserDataOperator::OPROW) {
        try {
          ierr = oit->opRhs(*op_data[0], false);
          if (PetscUnlikely(ierr)) {
            std::ostringstream ss;
            ss << "(Calling user data operator "
               << boost::typeindex::type_id_runtime(*oit).pretty_name() << ") "
               << PETSC_FUNCTION_NAME;
            return PetscError(PETSC_COMM_SELF, __LINE__, ss.str().c_str(),
                              __FILE__, ierr, PETSC_ERROR_REPEAT, " ");
          }
        } catch (std::exception &ex) {
          std::ostringstream ss;
          ss << "Operator "
             << boost::typeindex::type_id_runtime(*oit).pretty_name()
             << " operator number "
             << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(
                    opPtrVector.begin(), oit)
             << " thorw in method: " << ex.what() << " at line " << __LINE__
             << " in file " << __FILE__;
          SETERRQ(mField.get_comm(), MOFEM_STD_EXCEPTION_THROW,
                  ss.str().c_str());
        }
      }

      if (oit->getOpType() & UserDataOperator::OPCOL) {
        try {
          ierr = oit->opRhs(*op_data[1], false);
          if (PetscUnlikely(ierr)) {
            std::ostringstream ss;
            ss << "(Calling user data operator "
               << boost::typeindex::type_id_runtime(*oit).pretty_name() << ") "
               << PETSC_FUNCTION_NAME;
            return PetscError(PETSC_COMM_SELF, __LINE__, ss.str().c_str(),
                              __FILE__, ierr, PETSC_ERROR_REPEAT, " ");
          }
        } catch (std::exception &ex) {
          std::ostringstream ss;
          ss << "Operator "
             << boost::typeindex::type_id_runtime(*oit).pretty_name()
             << " operator number "
             << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(
                    opPtrVector.begin(), oit)
             << " thorw in method: " << ex.what() << " at line " << __LINE__
             << " in file " << __FILE__;
          SETERRQ(mField.get_comm(), MOFEM_STD_EXCEPTION_THROW,
                  ss.str().c_str());
        }
      }

      if (oit->getOpType() & UserDataOperator::OPROWCOL) {
        try {
          ierr = oit->opLhs(*op_data[0], *op_data[1], oit->sYmm);
          if (PetscUnlikely(ierr)) {
            std::ostringstream ss;
            ss << "(Calling user data operator "
               << boost::typeindex::type_id_runtime(*oit).pretty_name() << ") "
               << PETSC_FUNCTION_NAME;
            return PetscError(PETSC_COMM_SELF, __LINE__, ss.str().c_str(),
                              __FILE__, ierr, PETSC_ERROR_REPEAT, " ");
          }
        } catch (std::exception &ex) {
          std::ostringstream ss;
          ss << "Operator "
             << boost::typeindex::type_id_runtime(*oit).pretty_name()
             << " operator number "
             << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(
                    opPtrVector.begin(), oit)
             << " thorw in method: " << ex.what() << " at line " << __LINE__
             << " in file " << __FILE__;
          SETERRQ(mField.get_comm(), MOFEM_STD_EXCEPTION_THROW,
                  ss.str().c_str());
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VolumeElementForcesAndSourcesCore::UserDataOperator::
    getDivergenceOfHDivBaseFunctions(int side, EntityType type,
                                     DataForcesAndSourcesCore::EntData &data,
                                     int gg, VectorDouble &div) {
  MoFEMFunctionBegin;

  int nb_dofs = data.getFieldData().size();
  if (nb_dofs == 0)
    MoFEMFunctionReturnHot(0);

  if (data.getSpace() != HDIV && data.getSpace() != HCURL) {
    SETERRQ1(getVolumeFE()->mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "This function should be used for HDIV used but is used with %s",
             FieldSpaceNames[data.getSpace()]);
  }

  if ((unsigned int)nb_dofs != data.getDiffHdivN().size2() / 9) {
    SETERRQ3(getVolumeFE()->mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "Data inositency, wrong number of dofs  = %s "
             "%d != %d/9",
             FieldSpaceNames[data.getSpace()], nb_dofs,
             data.getDiffHdivN().size2());
  }

  div.resize(nb_dofs, false);

  FTensor::Tensor0<double *> t_div(&*div.data().begin());
  const double *grad_ptr = &data.getDiffHdivN()(gg, 0);
  FTensor::Tensor1<FTensor::PackPtr<const double *, 9>, 3> t_grad_base(
      grad_ptr, &grad_ptr[HDIV1_1], &grad_ptr[HDIV2_2]);
  for (int dd = 0; dd < nb_dofs; dd++) {
    t_div = t_grad_base(0) + t_grad_base(1) + t_grad_base(2);
    ++t_div;
    ++t_grad_base;
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VolumeElementForcesAndSourcesCore::UserDataOperator::
    getCurlOfHCurlBaseFunctions(int side, EntityType type,
                                DataForcesAndSourcesCore::EntData &data, int gg,
                                MatrixDouble &curl) {
  MoFEMFunctionBegin;

  int nb_dofs = data.getFieldData().size();
  if (nb_dofs == 0)
    MoFEMFunctionReturnHot(0);

  if (data.getSpace() != HDIV && data.getSpace() != HCURL) {
    SETERRQ1(getVolumeFE()->mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "This function should be used for primarily for HCURL"
             " but will work with HDIV used but is used with %s",
             FieldSpaceNames[data.getSpace()]);
  }

  if ((unsigned int)nb_dofs != data.getDiffHcurlN().size2() / 9) {
    SETERRQ3(getVolumeFE()->mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "Data insistency, wrong number of dofs  = %s "
             "%d != %d/9",
             FieldSpaceNames[data.getSpace()], nb_dofs,
             data.getDiffHcurlN().size2());
  }

  curl.resize(nb_dofs, 3, false);
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_curl(
      &curl(0, 0), &curl(0, 1), &curl(0, 2));
  const double *grad_ptr = &data.getDiffHcurlN()(gg, 0);

  FTensor::Tensor2<FTensor::PackPtr<const double *, 9>, 3, 3> t_grad_base(
      grad_ptr, &grad_ptr[HCURL0_1], &grad_ptr[HCURL0_2], &grad_ptr[HCURL1_0],
      &grad_ptr[HCURL1_1], &grad_ptr[HCURL1_2], &grad_ptr[HCURL2_0],
      &grad_ptr[HCURL2_1], &grad_ptr[HCURL2_2]);
  for (int dd = 0; dd != nb_dofs; ++dd) {
    t_curl(0) = t_grad_base(2, 1) - t_grad_base(1, 2);
    t_curl(1) = t_grad_base(0, 2) - t_grad_base(2, 0);
    t_curl(2) = t_grad_base(1, 0) - t_grad_base(0, 1);
    ++t_curl;
    ++t_grad_base;
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VolumeElementForcesAndSourcesCoreOnSide::setGaussPts(int order) {
  MoFEMFunctionBegin;
  if (faceFEPtr == NULL) {
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
            "Pointer to face element is not set");
  }
  const EntityHandle face_entity =
      faceFEPtr->numeredEntFiniteElementPtr->getEnt();
  SideNumber_multiIndex &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());
  SideNumber_multiIndex::nth_index<0>::type::iterator sit =
      side_table.get<0>().find(face_entity);
  if (sit == side_table.get<0>().end()) {
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
            "Face can not be found on volume element");
  }
  faceSense = (*sit)->sense;
  faceSideNumber = (*sit)->side_number;
  fill(tetConnMap, &tetConnMap[4], -1);
  for (int nn = 0; nn != 3; nn++) {
    faceConnMap[nn] =
        std::distance(conn, find(conn, &conn[4], faceFEPtr->conn[nn]));
    tetConnMap[faceConnMap[nn]] = nn;
    if (faceConnMap[nn] > 3) {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
              "No common node on face and element can not be found");
    }
  }
  oppositeNode =
      std::distance(tetConnMap, find(tetConnMap, &tetConnMap[4], -1));
  const int nb_gauss_pts = faceFEPtr->gaussPts.size2();
  gaussPts.resize(4, nb_gauss_pts, false);
  gaussPts.clear();
  const MatrixDouble &face_shape_funtions =
      faceFEPtr->dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE);
  const double tet_coords[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
  for (int gg = 0; gg != nb_gauss_pts; gg++) {
    gaussPts(0, gg) =
        face_shape_funtions(gg, 0) * tet_coords[3 * faceConnMap[0] + 0] +
        face_shape_funtions(gg, 1) * tet_coords[3 * faceConnMap[1] + 0] +
        face_shape_funtions(gg, 2) * tet_coords[3 * faceConnMap[2] + 0];
    gaussPts(1, gg) =
        face_shape_funtions(gg, 0) * tet_coords[3 * faceConnMap[0] + 1] +
        face_shape_funtions(gg, 1) * tet_coords[3 * faceConnMap[1] + 1] +
        face_shape_funtions(gg, 2) * tet_coords[3 * faceConnMap[2] + 1];
    gaussPts(2, gg) =
        face_shape_funtions(gg, 0) * tet_coords[3 * faceConnMap[0] + 2] +
        face_shape_funtions(gg, 1) * tet_coords[3 * faceConnMap[1] + 2] +
        face_shape_funtions(gg, 2) * tet_coords[3 * faceConnMap[2] + 2];
    gaussPts(3, gg) = faceFEPtr->gaussPts(2, gg);
  }
  MoFEMFunctionReturn(0);
}

VectorDouble &
VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::getNormal() {
  return getFaceFEPtr()->nOrmal;
}

MatrixDouble &VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::
    getNormalsAtGaussPts() {
  return getFaceFEPtr()->normalsAtGaussPts;
}

MatrixDouble &VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::
    getFaceCoordsAtGaussPts() {
  return getFaceFEPtr()->coordsAtGaussPts;
}

/** \brief if higher order geometry return normals at Gauss pts.
 *
 * \param gg gauss point number
 */
ublas::matrix_row<MatrixDouble>
VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::getNormalsAtGaussPts(
    const int gg) {
  return ublas::matrix_row<MatrixDouble>(getNormalsAtGaussPts(), gg);
}

} // namespace MoFEM
