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
EdgeElementForcesAndSourcesCore::calculateBaseFunctionsOnElement() {
  MoFEMFunctionBegin;
  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
    if (dataH1.bAse.test(b)) {
      switch (ApproximationBaseArray[b]) {
      case AINSWORTH_LEGENDRE_BASE:
      case AINSWORTH_LOBATTO_BASE:
      case DEMKOWICZ_JACOBI_BASE:
        if (dataH1.spacesOnEntities[MBVERTEX].test(H1) &&
            dataH1.basesOnEntities[MBVERTEX].test(b)) {
          CHKERR EdgePolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                  dataH1, H1, ApproximationBaseArray[b], NOBASE)));
        }
        if (dataH1.spacesOnEntities[MBEDGE].test(HCURL) &&
            dataH1.basesOnEntities[MBEDGE].test(b)) {
          CHKERR EdgePolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                  dataHcurl, HCURL, ApproximationBaseArray[b], NOBASE)));
        }
        if (dataH1.spacesOnEntities[MBEDGE].test(HDIV)) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Make no sense to have Hdiv space on edge (If you need h-div "
                  "space in 2d use h-curl space intead and use orthogonal "
                  "vector");
        }
        if (dataH1.spacesOnEntities[MBEDGE].test(L2)) {
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
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
EdgeElementForcesAndSourcesCore::calculateHoCoordsAtIntegrationPts() {
  MoFEMFunctionBegin;
  if (dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName) !=
      dataPtr->get<FieldName_mi_tag>().end()) {
    CHKERR getEdgesDataOrderSpaceAndBase(dataH1, meshPositionsFieldName);
    CHKERR getEdgesFieldData(dataH1, meshPositionsFieldName);
    CHKERR getNodesFieldData(dataH1, meshPositionsFieldName);
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

  CHKERR calculateEdgeDirection();
  CHKERR getSpacesAndBaseOnEntities(dataH1);
  CHKERR getEdgesDataOrder(dataH1, H1);
  dataH1.dataOnEntities[MBEDGE][0].getSense() =
      1; // set sense to 1, this is this entity

  // Hcurl
  if (dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
    CHKERR getEdgesDataOrder(dataHcurl, HCURL);
    dataHcurl.dataOnEntities[MBEDGE][0].getSense() =
        1; // set sense to 1, this is this entity
    dataHcurl.spacesOnEntities[MBEDGE].set(HCURL);
  }

  /// Use the some node base
  dataHcurl.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) =
      dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);

  CHKERR setIntegrationPts();
  CHKERR calculateCoordsAtIntegrationPts();
  CHKERR calculateBaseFunctionsOnElement();
  CHKERR calculateHoCoordsAtIntegrationPts();

  if (dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
    // cerr << dataHcurl.dataOnEntities[MBEDGE][0].getN(AINSWORTH_LEGENDRE_BASE)
    // << endl;
    CHKERR opCovariantTransform.opRhs(dataHcurl);
  }

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

    for (int ss = 0; ss != 2; ss++) {

      std::string field_name = !ss ? oit->rowFieldName : oit->colFieldName;
      const Field *field_struture = mField.get_field_structure(field_name);
      BitFieldId data_id = field_struture->getId();

      if ((oit->getNumeredEntFiniteElementPtr()->getBitFieldIdData() & data_id)
              .none()) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "no data field < %s > on finite element < %s >",
                 field_name.c_str(), feName.c_str());
      }

      if (oit->getOpType() & types[ss] ||
          oit->getOpType() & UserDataOperator::OPROWCOL) {

        space[ss] = field_struture->getSpace();
        switch (space[ss]) {
        case NOSPACE:
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown space");
        case H1:
          op_data[ss] = !ss ? &dataH1 : &derivedDataH1;
          break;
        case HCURL:
          op_data[ss] = !ss ? &dataHcurl : &derivedDataHcurl;
          break;
        case HDIV:
          SETERRQ(
              PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
              "not make sanes on edge in 3d space (for 1d/2d not implemented)");
          break;
        case L2:
          SETERRQ(
              PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
              "not make sanes on edge in 3d space (for 1d/2d not implemented)");
          break;
        case NOFIELD:
          op_data[ss] = !ss ? &dataNoField : &dataNoFieldCol;
          break;
        case LASTSPACE:
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown space");
          break;
        }

        base[ss] = field_struture->getApproxBase();
        switch (base[ss]) {
        case AINSWORTH_LEGENDRE_BASE:
        case AINSWORTH_LOBATTO_BASE:
        case DEMKOWICZ_JACOBI_BASE:
          break;
        default:
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "unknown or not implemented base");
          break;
        }

        if (last_eval_field_name[ss] != field_name) {
          switch (space[ss]) {
          case NOSPACE:
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown space");
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
            break;
          case HDIV:
            SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                    "not make sanes on edge");
            break;
          case L2:
            SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                    "not make sanes on edge");
            break;
          case NOFIELD:
            if (!getNinTheLoop()) {
              // NOFIELD data are the same for each element, can be retrieved
              // only once
              if (!ss) {
                CHKERR getNoFieldRowIndices(*op_data[ss], field_name);
              } else {
                CHKERR getNoFieldColIndices(*op_data[ss], field_name);
              }
              CHKERR getNoFieldFieldData(*op_data[ss], field_name);
            }
            break;
          case LASTSPACE:
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown space");
            break;
          }
          last_eval_field_name[ss] = field_name;
        }
      }
    }

    if (oit->getOpType() & UserDataOperator::OPROW) {
      try {
        ierr = oit->opRhs(*op_data[0], oit->doVertices, oit->doEdges, false,
                          false, false, false);
        CHKERRG(ierr);
      } catch (std::exception &ex) {
        std::ostringstream ss;
        ss << "Operator "
           << boost::typeindex::type_id_runtime(*oit).pretty_name()
           << " operator number "
           << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(
                  opPtrVector.begin(), oit)
           << " throw in method: " << ex.what() << " at line " << __LINE__
           << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF, MOFEM_STD_EXCEPTION_THROW, ss.str().c_str());
      }
    }

    if (oit->getOpType() & UserDataOperator::OPCOL) {
      try {
        ierr = oit->opRhs(*op_data[1], oit->doVertices, oit->doEdges, false,
                          false, false, false);
        CHKERRG(ierr);
      } catch (std::exception &ex) {
        std::ostringstream ss;
        ss << "Operator "
           << boost::typeindex::type_id_runtime(*oit).pretty_name()
           << " operator number "
           << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(
                  opPtrVector.begin(), oit)
           << " throw in method: " << ex.what() << " at line " << __LINE__
           << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF, MOFEM_STD_EXCEPTION_THROW, ss.str().c_str());
      }
    }

    if (oit->getOpType() & UserDataOperator::OPROWCOL) {
      try {
        ierr = oit->opLhs(*op_data[0], *op_data[1], oit->sYmm);
        CHKERRG(ierr);
      } catch (std::exception &ex) {
        std::ostringstream ss;
        ss << "Operator "
           << boost::typeindex::type_id_runtime(*oit).pretty_name()
           << " operator number "
           << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(
                  opPtrVector.begin(), oit)
           << " throw in method: " << ex.what() << " at line " << __LINE__
           << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF, MOFEM_STD_EXCEPTION_THROW, ss.str().c_str());
      } 
    }
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
