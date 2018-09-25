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
              &invJac(2, 2)) {
  getElementPolynomialBase() =
      boost::shared_ptr<BaseFunction>(new TetPolynomialBase());
}

MoFEMErrorCode VolumeElementForcesAndSourcesCore::setIntegrationPts() {
  MoFEMFunctionBegin;
  int order_data = getMaxDataOrder();
  int order_row = getMaxRowOrder();
  int order_col = getMaxColOrder();
  int rule = getRule(order_row, order_col, order_data);

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
    CHKERR getEntitySense<MBTRI>(dataHcurl);
    CHKERR getFaceTriNodes(dataHcurl);
    CHKERR getEntityDataOrder<MBTRI>(dataHcurl, HCURL);
    dataHcurl.spacesOnEntities[MBTRI].set(HCURL);
  }
  if ((dataH1.spacesOnEntities[MBTET]).test(HCURL)) {
    CHKERR getEntityDataOrder<MBTET>(dataHcurl, HCURL);
    dataHcurl.spacesOnEntities[MBTET].set(HCURL);
  }
  // Hdiv
  if ((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
    CHKERR getEntitySense<MBTRI>(dataHdiv);
    CHKERR getFaceTriNodes(dataHdiv);
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

    CHKERR getEntityDataOrderSpaceAndBase<MBEDGE>(dataH1,
                                                  meshPositionsFieldName);
    CHKERR getEntityDataOrderSpaceAndBase<MBTRI>(dataH1,
                                                 meshPositionsFieldName);
    CHKERR getEntityDataOrderSpaceAndBase<MBTET>(dataH1,
                                                 meshPositionsFieldName);
    CHKERR getNodesFieldData(dataH1, meshPositionsFieldName);
    if (dataH1.dataOnEntities[MBVERTEX][0].getFieldData().size() != 12) {
      SETERRQ(mField.get_comm(), MOFEM_NOT_FOUND,
              "no MESH_NODE_POSITIONS in element data or field has wrong "
              "number of coefficients");
    }
    CHKERR getEntityFieldData<MBEDGE>(dataH1, meshPositionsFieldName);
    CHKERR getEntityFieldData<MBTRI>(dataH1, meshPositionsFieldName);
    CHKERR getEntityFieldData<MBTET>(dataH1, meshPositionsFieldName);
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
  CHKERR createDataOnElement();

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
  }
  CATCH_ERRORS;

  CHKERR calculateHoJacobian();
  CHKERR transformHoBaseFunctions();

  // Iterate over operators
  CHKERR loopOverOperators();

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
  DataForcesAndSourcesCore &dataH1_on_face = *faceFEPtr->dataOnElement[H1];
  const MatrixDouble &face_shape_funtions =
      dataH1_on_face.dataOnEntities[MBVERTEX][0].getN(NOBASE);
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
