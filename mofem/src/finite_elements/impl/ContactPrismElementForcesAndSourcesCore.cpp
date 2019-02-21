/** \file ContactPrismElementForcesAndSourcesCore.cpp

\brief Implementation of contact prism element

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
#include <ContactPrismElementForcesAndSourcesCore.hpp>

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

ContactPrismElementForcesAndSourcesCore::
    ContactPrismElementForcesAndSourcesCore(Interface &m_field)
    : ForcesAndSourcesCore(m_field),
      meshPositionsFieldName("MESH_NODE_POSITIONS"),
                           dataOnMaster{

          nullptr,
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)), // NOFIELD
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)), // H1
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)), // HCURL
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)), // HDIV
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)) // L2

      },
      derivedDataOnMaster{

          nullptr,
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnMaster[NOFIELD])),
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnMaster[H1])),
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnMaster[HCURL])),
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnMaster[HDIV])),
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnMaster[L2]))

      },dataOnSlave{

          nullptr,
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)), // NOFIELD
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)), // H1
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)), // HCURL
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)), // HDIV
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)) // L2

      },
      derivedDataOnSlave{

          nullptr,
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnSlave[NOFIELD])),
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnSlave[H1])),
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnSlave[HCURL])),
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnSlave[HDIV])),
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnSlave[L2]))

      },
      dataH1Master(*dataOnMaster[H1].get()),
      dataH1Slave(*dataOnSlave[H1].get()),
      dataNoFieldSlave(*dataOnSlave[NOFIELD].get()),
      dataNoFieldMaster(*dataOnMaster[NOFIELD].get()),
      dataHcurlMaster(*dataOnMaster[HCURL].get()),dataHcurlSlave(*dataOnSlave[HCURL].get()),
      dataHdivMaster(*dataOnMaster[HDIV].get()), dataL2Master(*dataOnMaster[L2].get()),
      dataHdivSlave(*dataOnSlave[HDIV].get()), dataL2Slave(*dataOnSlave[L2].get()) {

  getUserPolynomialBase() =
      boost::shared_ptr<BaseFunction>(new TriPolynomialBase());

  dataH1Master.dataOnEntities[MBVERTEX].push_back(
      new DataForcesAndSourcesCore::EntData());
  dataH1Slave.dataOnEntities[MBVERTEX].push_back(
      new DataForcesAndSourcesCore::EntData());

  for (int ee = 0; ee != 3; ++ee) {
    dataH1Master.dataOnEntities[MBEDGE].push_back(
        new DataForcesAndSourcesCore::EntData());
    dataH1Slave.dataOnEntities[MBEDGE].push_back(
        new DataForcesAndSourcesCore::EntData());
  }
 
    dataH1Master.dataOnEntities[MBTRI].push_back(
        new DataForcesAndSourcesCore::EntData());
    dataH1Slave.dataOnEntities[MBTRI].push_back(
        new DataForcesAndSourcesCore::EntData());
}

MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::operator()() {
  MoFEMFunctionBegin;

  if (numeredEntFiniteElementPtr->getEntType() != MBPRISM)
    MoFEMFunctionReturnHot(0);
  CHKERR createDataOnElement();

  // Data on elements for proper spaces
   dataOnMaster[H1]->setElementType(MBTRI);
   derivedDataOnMaster[H1]->setElementType(MBTRI);
   dataOnSlave[H1]->setElementType(MBTRI);
   derivedDataOnSlave[H1]->setElementType(MBTRI);

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



auto clean_data = [](DataForcesAndSourcesCore &data) {
  MoFEMFunctionBegin;
    data.sPace.reset();
    data.bAse.reset();
    for (EntityType t = MBVERTEX; t != MBMAXTYPE; ++t) {
      data.spacesOnEntities[t].reset();
      data.basesOnEntities[t].reset();
    }
    for (int s = 0; s != LASTSPACE; ++s) {
      data.basesOnSpaces[s].reset();
    }

    MoFEMFunctionReturn(0);
};

auto copy_data = [](DataForcesAndSourcesCore &data, DataForcesAndSourcesCore &copy_data, const int &shift) {
  MoFEMFunctionBegin;

    if (shift != 0 && shift != 6) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Wrong shift for contact prism element");
    }

      data.sPace = copy_data.sPace;
      data.bAse = copy_data.bAse;
      data.spacesOnEntities[MBVERTEX] = copy_data.spacesOnEntities[MBVERTEX];
      data.spacesOnEntities[MBEDGE] = copy_data.spacesOnEntities[MBEDGE];
      data.spacesOnEntities[MBTRI] = copy_data.spacesOnEntities[MBTRI];
      
      data.basesOnEntities[MBVERTEX] = copy_data.basesOnEntities[MBVERTEX];
      data.basesOnEntities[MBEDGE] = copy_data.basesOnEntities[MBEDGE];
      data.basesOnEntities[MBTRI] = copy_data.basesOnEntities[MBTRI];
      
      data.basesOnSpaces[MBVERTEX] = copy_data.basesOnSpaces[MBVERTEX];
      data.basesOnSpaces[MBEDGE] = copy_data.basesOnSpaces[MBEDGE];
      data.basesOnSpaces[MBTRI] = copy_data.basesOnSpaces[MBTRI];


    for(int ii = 0; ii != 3; ++ii ){
     data.dataOnEntities[MBEDGE][ii].getSense() = copy_data.dataOnEntities[MBEDGE][ii + shift].getSense();
    }

    if(shift == 0){
    data.dataOnEntities[MBTRI][0].getSense() = copy_data.dataOnEntities[MBTRI][3].getSense();
    } else {
    data.dataOnEntities[MBTRI][0].getSense() = copy_data.dataOnEntities[MBTRI][4].getSense();
    }

    MoFEMFunctionReturn(0);
};

  CHKERR clean_data(dataH1Slave);
  CHKERR copy_data(dataH1Slave, dataH1, 6);
  CHKERR clean_data(dataH1Master);
  CHKERR copy_data(dataH1Master, dataH1, 0);


  int order_data = getMaxDataOrder();
  int order_row = getMaxRowOrder();
  int order_col = getMaxColOrder(); // maybe two different rules?
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

      // For master and slave
      gaussPtsMaster.resize(3, nb_gauss_pts / 2, false);
      gaussPtsSlave.resize(3, nb_gauss_pts / 2, false);

      cblas_dcopy(nb_gauss_pts / 2, &QUAD_2D_TABLE[rule]->points[1], 3,
                  &gaussPtsMaster(0, 0), 1);
      cblas_dcopy(nb_gauss_pts / 2, &QUAD_2D_TABLE[rule]->points[2], 3,
                  &gaussPtsMaster(1, 0), 1);
      cblas_dcopy(nb_gauss_pts / 2, QUAD_2D_TABLE[rule]->weights, 1,
                  &gaussPtsMaster(2, 0), 1);

      cblas_dcopy(nb_gauss_pts / 2,
                  &QUAD_2D_TABLE[rule]->points[3 * nb_gauss_pts / 2 + 1], 3,
                  &gaussPtsSlave(0, 0), 1);
      cblas_dcopy(nb_gauss_pts / 2,
                  &QUAD_2D_TABLE[rule]->points[3 * nb_gauss_pts / 2 + 2], 3,
                  &gaussPtsSlave(1, 0), 1);
      cblas_dcopy(nb_gauss_pts / 2,
                  &QUAD_2D_TABLE[rule]->weights[nb_gauss_pts / 2], 1,
                  &gaussPtsSlave(2, 0), 1);

      dataH1Master.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(
          nb_gauss_pts / 2, 2, false);

      dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(
          nb_gauss_pts / 2, 2, false);

      double *shape_ptr_master = &*dataH1Master.dataOnEntities[MBVERTEX][0]
                                       .getN(NOBASE)
                                       .data()
                                       .begin();
      cblas_dcopy(3 * nb_gauss_pts / 2, &QUAD_2D_TABLE[rule]->points[1], 1,
                  shape_ptr_master, 1);
      double *shape_ptr_slave =
          &*dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      cblas_dcopy(3 * nb_gauss_pts / 2, &QUAD_2D_TABLE[rule]->points[2], 1,
                  shape_ptr_slave, 1);
    } else {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "rule > quadrature order %d < %d", rule, QUAD_2D_TABLE_SIZE);
      nb_gauss_pts = 0;
    }
  } else {
    // Master-Slave
    CHKERR setGaussPts(order_row, order_col, order_data);

    dataH1Master.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(
        nb_gauss_pts / 2, 3, false);
    dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(
        nb_gauss_pts / 2, 3, false);

    if (nb_gauss_pts) {
      CHKERR ShapeMBTRI(&*dataH1Master.dataOnEntities[MBVERTEX][0]
                              .getN(NOBASE)
                              .data()
                              .begin(),
                        &gaussPtsMaster(0, 0), &gaussPtsMaster(1, 0),
                        nb_gauss_pts / 2);

      CHKERR ShapeMBTRI(
          &*dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
          &gaussPtsSlave(0, 0), &gaussPtsSlave(1, 0), nb_gauss_pts / 2);
    }

    if (nb_gauss_pts) {
      CHKERR ShapeMBTRI(&*dataH1Master.dataOnEntities[MBVERTEX][0]
                              .getN(NOBASE)
                              .data()
                              .begin(),
                        &gaussPtsMaster(0, 0), &gaussPts(1, 0),
                        nb_gauss_pts / 2);

      CHKERR ShapeMBTRI(
          &*dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
          &gaussPtsSlave(0, 0), &gaussPts(1, 0), nb_gauss_pts / 2);
    }
  }
  if (nb_gauss_pts == 0)
    MoFEMFunctionReturnHot(0);

  {
    coordsAtGaussPtsMaster.resize(nb_gauss_pts / 2, 3, false);
    coordsAtGaussPtsSlave.resize(nb_gauss_pts / 2, 3, false);
    coordsAtGaussPts.resize(nb_gauss_pts, 6, false);
    for (int gg = 0; gg < nb_gauss_pts / 2; gg++) {
      for (int dd = 0; dd < 3; dd++) {
        coordsAtGaussPtsMaster(gg, dd) = cblas_ddot(
            3, &dataH1Master.dataOnEntities[MBVERTEX][0].getN(NOBASE)(gg, 0), 1,
            &coords[dd], 3);
        coordsAtGaussPtsSlave(gg, dd) = cblas_ddot(
            3, &dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE)(gg, 0), 1,
            &coords[9 + dd], 3);
      }
    }
  }

  for (int space = HCURL; space != LASTSPACE; ++space)
    if (dataOnElement[space]) {
      // dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) =
      //     dataOnElement[H1]->dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
      dataH1Master.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) =
          dataOnMaster[H1]->dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
      dataH1Slave.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) =
          dataOnSlave[H1]->dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
    }

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
    if (dataH1.bAse.test(b)) {
      switch (static_cast<FieldApproximationBase>(b)) {
      case AINSWORTH_LEGENDRE_BASE:
      case AINSWORTH_LOBATTO_BASE:
        if (dataH1.spacesOnEntities[MBVERTEX].test(H1)) {
   
          CHKERR getUserPolynomialBase()->getValue(
              gaussPtsMaster,
              boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                  dataH1Master, H1, static_cast<FieldApproximationBase>(b),
                  NOBASE)));

          CHKERR getUserPolynomialBase()->getValue(
              gaussPtsSlave,
              boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                  dataH1Slave, H1, static_cast<FieldApproximationBase>(b),
                  NOBASE)));
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

 

  if (dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented yet");
  }

  // Iterate over operators
  CHKERR loopOverOperators();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpCalculateInvJacForContactPrism::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {

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
    invJacMaster.resize(2, 2, false);
    invJacMaster(0, 0) = j11_f3 / det_f3;
    invJacMaster(0, 1) = -j01_f3 / det_f3;
    invJacMaster(1, 0) = -j10_f3 / det_f3;
    invJacMaster(1, 1) = j00_f3 / det_f3;
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
OpSetInvJacH1ForContactPrism::doWork(int side, EntityType type,
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
                      &*invJacMaster.data().begin(), 2,
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
