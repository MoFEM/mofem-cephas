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

namespace MoFEM {

ContactPrismElementForcesAndSourcesCore::
    ContactPrismElementForcesAndSourcesCore(Interface &m_field)
    : ForcesAndSourcesCore(m_field),
      dataOnMaster{

          nullptr,
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET), // NOFIELD
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET), // H1
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET), // HCURL
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET), // HDIV
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET)  // L2

      },
      dataOnSlave{

          nullptr,
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET), // NOFIELD
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET), // H1
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET), // HCURL
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET), // HDIV
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET)  // L2

      },
      derivedDataOnMaster{

          nullptr,
          boost::make_shared<DerivedDataForcesAndSourcesCore>(
              dataOnMaster[NOFIELD]),
          boost::make_shared<DerivedDataForcesAndSourcesCore>(dataOnMaster[H1]),
          boost::make_shared<DerivedDataForcesAndSourcesCore>(
              dataOnMaster[HCURL]),
          boost::make_shared<DerivedDataForcesAndSourcesCore>(
              dataOnMaster[HDIV]),
          boost::make_shared<DerivedDataForcesAndSourcesCore>(dataOnMaster[L2])

      },
      derivedDataOnSlave{

          nullptr,
          boost::make_shared<DerivedDataForcesAndSourcesCore>(
              dataOnSlave[NOFIELD]),
          boost::make_shared<DerivedDataForcesAndSourcesCore>(dataOnSlave[H1]),
          boost::make_shared<DerivedDataForcesAndSourcesCore>(
              dataOnSlave[HCURL]),
          boost::make_shared<DerivedDataForcesAndSourcesCore>(
              dataOnSlave[HDIV]),
          boost::make_shared<DerivedDataForcesAndSourcesCore>(dataOnSlave[L2])

      },
      dataH1Master(*dataOnMaster[H1].get()),
      dataH1Slave(*dataOnSlave[H1].get()),
      dataNoFieldSlave(*dataOnSlave[NOFIELD].get()),
      dataNoFieldMaster(*dataOnMaster[NOFIELD].get()),
      dataHcurlMaster(*dataOnMaster[HCURL].get()),
      dataHcurlSlave(*dataOnSlave[HCURL].get()),
      dataHdivMaster(*dataOnMaster[HDIV].get()),
      dataL2Master(*dataOnMaster[L2].get()),
      dataHdivSlave(*dataOnSlave[HDIV].get()),
      dataL2Slave(*dataOnSlave[L2].get()) {

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

  // Data on elements for proper spaces
  dataOnMaster[H1]->setElementType(MBTRI);
  derivedDataOnMaster[H1]->setElementType(MBTRI);
  dataOnSlave[H1]->setElementType(MBTRI);
  derivedDataOnSlave[H1]->setElementType(MBTRI);
}

MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::operator()() {
  MoFEMFunctionBegin;

  if (numeredEntFiniteElementPtr->getEntType() != MBPRISM)
    MoFEMFunctionReturnHot(0);
  CHKERR createDataOnElement();

  DataForcesAndSourcesCore &data_div = *dataOnElement[HDIV];
  DataForcesAndSourcesCore &data_curl = *dataOnElement[HCURL];
  DataForcesAndSourcesCore &data_l2 = *dataOnElement[HCURL];

  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  auto get_coord_and_normal = [&]() {
    MoFEMFunctionBegin;
    int num_nodes;
    const EntityHandle *conn;
    CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
    coords.resize(num_nodes * 3, false);
    CHKERR mField.get_moab().get_coords(conn, num_nodes,
                                        &*coords.data().begin());
    normal.resize(6, false);
    CHKERR Tools::getTriNormal(&coords[0], &normal[0]);
    CHKERR Tools::getTriNormal(&coords[9], &normal[3]);
    aRea[0] = cblas_dnrm2(3, &normal[0], 1) * 0.5;
    aRea[1] = cblas_dnrm2(3, &normal[3], 1) * 0.5;
    MoFEMFunctionReturn(0);
  };
  CHKERR get_coord_and_normal();

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
    for (int s = 0; s != LASTSPACE; ++s)
      data.basesOnSpaces[s].reset();

    MoFEMFunctionReturn(0);
  };

  auto copy_data = [](DataForcesAndSourcesCore &data,
                      DataForcesAndSourcesCore &copy_data, const int shift) {
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

    for (int ii = 0; ii != 3; ++ii) {
      data.dataOnEntities[MBEDGE][ii].getSense() =
          copy_data.dataOnEntities[MBEDGE][ii + shift].getSense();
      data.dataOnEntities[MBEDGE][ii].getDataOrder() =
          copy_data.dataOnEntities[MBEDGE][ii + shift].getDataOrder();
    }

    if (shift == 0) {
      data.dataOnEntities[MBTRI][0].getSense() =
          copy_data.dataOnEntities[MBTRI][3].getSense();
      data.dataOnEntities[MBTRI][0].getDataOrder() =
          copy_data.dataOnEntities[MBTRI][3].getDataOrder();
    } else {
      data.dataOnEntities[MBTRI][0].getSense() =
          copy_data.dataOnEntities[MBTRI][4].getSense();
      data.dataOnEntities[MBTRI][0].getDataOrder() =
          copy_data.dataOnEntities[MBTRI][4].getDataOrder();
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
      gaussPtsMaster.resize(3, nb_gauss_pts, false);
      gaussPtsSlave.resize(3, nb_gauss_pts, false);

      cblas_dcopy(nb_gauss_pts, &QUAD_2D_TABLE[rule]->points[1], 3,
                  &gaussPtsMaster(0, 0), 1);
      cblas_dcopy(nb_gauss_pts, &QUAD_2D_TABLE[rule]->points[2], 3,
                  &gaussPtsMaster(1, 0), 1);
      cblas_dcopy(nb_gauss_pts, QUAD_2D_TABLE[rule]->weights, 1,
                  &gaussPtsMaster(2, 0), 1);

      gaussPtsSlave = gaussPtsMaster;

      dataH1Master.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts,
                                                                   3, false);

      dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts,
                                                                  3, false);

      double *shape_ptr_master = &*dataH1Master.dataOnEntities[MBVERTEX][0]
                                       .getN(NOBASE)
                                       .data()
                                       .begin();
      cblas_dcopy(3 * nb_gauss_pts, QUAD_2D_TABLE[rule]->points, 1,
                  shape_ptr_master, 1);
      double *shape_ptr_slave =
          &*dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      cblas_dcopy(3 * nb_gauss_pts, QUAD_2D_TABLE[rule]->points, 1,
                  shape_ptr_slave, 1);
    } else {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "rule > quadrature order %d < %d", rule, QUAD_2D_TABLE_SIZE);
      nb_gauss_pts = 0;
    }

  } else {

    // Master-Slave
    if (gaussPtsMaster.size2() != gaussPtsSlave.size2())
      SETERRQ(
          PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
          "Number of Gauss Points at Master triangle is different than slave");

    CHKERR setGaussPts(order_row, order_col, order_data);
    nb_gauss_pts = gaussPtsMaster.size2();
    dataH1Master.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts,
                                                                 3, false);
    dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 3,
                                                                false);

    if (nb_gauss_pts) {
      CHKERR ShapeMBTRI(&*dataH1Master.dataOnEntities[MBVERTEX][0]
                              .getN(NOBASE)
                              .data()
                              .begin(),
                        &gaussPtsMaster(0, 0), &gaussPtsMaster(1, 0),
                        nb_gauss_pts);

      CHKERR ShapeMBTRI(
          &*dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
          &gaussPtsSlave(0, 0), &gaussPtsSlave(1, 0), nb_gauss_pts);
    }
  }
  if (nb_gauss_pts == 0)
    MoFEMFunctionReturnHot(0);

  {
    coordsAtGaussPtsMaster.resize(nb_gauss_pts, 3, false);
    coordsAtGaussPtsSlave.resize(nb_gauss_pts, 3, false);
    for (int gg = 0; gg < nb_gauss_pts; gg++) {
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

MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::loopOverOperators() {
  MoFEMFunctionBegin;

  const EntityType type = numeredEntFiniteElementPtr->getEntType();
  const UserDataOperator::OpType types[2] = {UserDataOperator::OPROW,
                                             UserDataOperator::OPCOL};
  std::vector<std::string> last_eval_field_name(2);

  auto oit = opPtrVector.begin();
  auto hi_oit = opPtrVector.end();

  for (; oit != hi_oit; oit++) {

    oit->setPtrFE(this);

    if (oit->opType == UserDataOperator::OPLAST) {

      // Set field
      switch (oit->sPace) {
      case NOSPACE:
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Unknown space");
      case NOFIELD:
      case H1:
      case HCURL:
      case HDIV:
      case L2:
        break;
      default:
        SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                 "Not implemented for this space", oit->sPace);
      }

      // Reseat all data which all field dependent
      dataOnMaster[oit->sPace]->resetFieldDependentData();
      dataOnSlave[oit->sPace]->resetFieldDependentData();
      last_eval_field_name[0] = "";

      // Run operator
      try {
        CHKERR oit->opRhs(*dataOnMaster[oit->sPace], oit->doVertices,
                          oit->doEdges, oit->doQuads, oit->doTris, oit->doTets,
                          false);

        CHKERR oit->opRhs(*dataOnSlave[oit->sPace], oit->doVertices,
                          oit->doEdges, oit->doQuads, oit->doTris, oit->doTets,
                          false);
      }
      CATCH_OP_ERRORS(*oit);

    } else {
      boost::shared_ptr<DataForcesAndSourcesCore> op_master_data[2];
      boost::shared_ptr<DataForcesAndSourcesCore> op_slave_data[2];

      for (int ss = 0; ss != 2; ss++) {

        const std::string field_name =
            !ss ? oit->rowFieldName : oit->colFieldName;
        const Field *field_struture = mField.get_field_structure(field_name);
        const BitFieldId data_id = field_struture->getId();
        const FieldSpace space = field_struture->getSpace();
        op_master_data[ss] =
            !ss ? dataOnMaster[space] : derivedDataOnMaster[space];
        op_slave_data[ss] =
            !ss ? dataOnSlave[space] : derivedDataOnSlave[space];

        if ((oit->getNumeredEntFiniteElementPtr()->getBitFieldIdData() &
             data_id)
                .none()) {
          SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "no data field < %s > on finite element < %s >",
                   field_name.c_str(), feName.c_str());
        }

        if (oit->getOpType() & types[ss] ||
            oit->getOpType() & UserDataOperator::OPROWCOL) {

          if (UserDataOperator *cast_oit =
                  dynamic_cast<UserDataOperator *>(&*oit)) {
          } else {
            printf("Check\n");
          }
          switch (space) {
          case NOSPACE:
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown space");
            break;
          case NOFIELD:
          case H1:
          case HCURL:
          case HDIV:
          case L2:
            break;
          default:
            SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                     "Not implemented for this space", space);
          }

          if (last_eval_field_name[ss] != field_name) {
            CHKERR getEntityFieldData(*op_master_data[ss], field_name, MBEDGE);
            CHKERR getEntityFieldData(*op_slave_data[ss], field_name, MBEDGE,
                                      MBPOLYHEDRON, false);
            if (!ss) {
              CHKERR getEntityRowIndices<true>(*op_master_data[ss], field_name,
                                         MBEDGE);
              CHKERR getEntityRowIndices<false>(*op_slave_data[ss], field_name,
                                                MBEDGE, MBPOLYHEDRON);
            } else {
              CHKERR getEntityColIndices<true>(*op_master_data[ss], field_name,
                                         MBEDGE);
              CHKERR getEntityColIndices<false>(*op_slave_data[ss], field_name,
                                                MBEDGE, MBPOLYHEDRON);
            }

            switch (space) {
            case NOSPACE:
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                      "unknown space");
              break;
            case H1:
              if (!ss) {
                CHKERR getRowNodesIndices<true>(*op_master_data[ss],
                                                field_name);
                CHKERR getRowNodesIndices<false>(*op_slave_data[ss],
                                                 field_name);
              } else {
                CHKERR getColNodesIndices<true>(*op_master_data[ss],
                                                field_name);
                CHKERR getColNodesIndices<false>(*op_slave_data[ss],
                                                 field_name);
              }
              CHKERR getNodesFieldData(*op_master_data[ss], field_name, true);
              CHKERR getNodesFieldData(*op_slave_data[ss], field_name, false);
              break;
            case HCURL:
            case HDIV:
              break;
            case L2:
              switch (type) {
              case MBVERTEX:
                CHKERR getNodesFieldData(*op_master_data[ss], field_name, true);
                CHKERR getNodesFieldData(*op_slave_data[ss], field_name, false);
                break;
              default:
                break;
              }
              break;
            case NOFIELD:
              if (!getNinTheLoop()) {
                // NOFIELD data are the same for each element, can be
                // retrieved only once
                if (!ss) {
                  CHKERR getNoFieldRowIndices(*op_master_data[ss], field_name);
                  CHKERR getNoFieldRowIndices(*op_slave_data[ss], field_name);
                } else {
                  CHKERR getNoFieldColIndices(*op_master_data[ss], field_name);
                  CHKERR getNoFieldColIndices(*op_slave_data[ss], field_name);
                }
                CHKERR getNoFieldFieldData(*op_master_data[ss], field_name);
                CHKERR getNoFieldFieldData(*op_slave_data[ss], field_name);
              }
              break;
            default:
              SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                       "Not implemented for this space", space);
            }
            last_eval_field_name[ss] = field_name;
          }
        }
      }

      int type;

      if (UserDataOperator *cast_oit =
              dynamic_cast<UserDataOperator *>(&*oit)) {
        type = cast_oit->getFaceType();
        if (((oit->getOpType() & UserDataOperator::OPROW) ||
             (oit->getOpType() & UserDataOperator::OPCOL)) &&
            ((type & UserDataOperator::FACEMASTERMASTER) ||
             (type & UserDataOperator::FACEMASTERSLAVE) ||
             (type & UserDataOperator::FACESLAVEMASTER) ||
             (type & UserDataOperator::FACESLAVESLAVE))) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Wrong combination of FaceType and OpType, OPROW or OPCOL "
                  "combined with face-face OpType");
        }
        if (!type) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Face type is not set");
        }
      } else {
        type = UserDataOperator::FACEMASTER | UserDataOperator::FACESLAVE |
               UserDataOperator::FACEMASTERMASTER |
               UserDataOperator::FACEMASTERSLAVE |
               UserDataOperator::FACESLAVEMASTER |
               UserDataOperator::FACESLAVESLAVE;
      }

      if (oit->getOpType() & UserDataOperator::OPROW &&
          (type & UserDataOperator::FACEMASTER)) {
        try {
          CHKERR oit->opRhs(*op_master_data[0], false);
        }
        CATCH_OP_ERRORS(*oit);
      }

      if (oit->getOpType() & UserDataOperator::OPROW &&
          (type & UserDataOperator::FACESLAVE)) {
        try {
          CHKERR oit->opRhs(*op_slave_data[0], false);
        }
        CATCH_OP_ERRORS(*oit);
      }

      if (oit->getOpType() & UserDataOperator::OPCOL &&
          (type & UserDataOperator::FACEMASTER)) {
        try {
          CHKERR oit->opRhs(*op_master_data[1], false);
        }
        CATCH_OP_ERRORS(*oit);
      }

      if (oit->getOpType() & UserDataOperator::OPCOL &&
          (type & UserDataOperator::FACESLAVE)) {
        try {
          CHKERR oit->opRhs(*op_slave_data[1], false);
        }
        CATCH_OP_ERRORS(*oit);
      }

      if (oit->getOpType() & UserDataOperator::OPROWCOL &&
          (type & UserDataOperator::FACEMASTERMASTER)) {
        try {
          CHKERR oit->opLhs(*op_master_data[0], *op_master_data[1], oit->sYmm);
        }
        CATCH_OP_ERRORS(*oit);
      }

      if (oit->getOpType() & UserDataOperator::OPROWCOL &&
          (type & UserDataOperator::FACEMASTERSLAVE)) {
        try {
          CHKERR oit->opLhs(*op_master_data[0], *op_slave_data[1], oit->sYmm);
        }
        CATCH_OP_ERRORS(*oit);
      }

      if (oit->getOpType() & UserDataOperator::OPROWCOL &&
          (type & UserDataOperator::FACESLAVEMASTER)) {
        try {
          CHKERR oit->opLhs(*op_slave_data[0], *op_master_data[1], oit->sYmm);
        }
        CATCH_OP_ERRORS(*oit);
      }

      if (oit->getOpType() & UserDataOperator::OPROWCOL &&
          (type & UserDataOperator::FACESLAVESLAVE)) {
        try {
          CHKERR oit->opLhs(*op_slave_data[0], *op_slave_data[1], oit->sYmm);
        }
        CATCH_OP_ERRORS(*oit);
      }
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::getEntityFieldData(
    DataForcesAndSourcesCore &data, const std::string &field_name,
    const EntityType type_lo, const EntityType type_hi,
    const bool master_flag) const {
  MoFEMFunctionBegin;

  for (EntityType t = type_lo; t != type_hi; ++t) {
    for (auto &dat : data.dataOnEntities[t]) {
      dat.getDataOrder() = 0;
      dat.getBase() = NOBASE;
      dat.getSpace() = NOSPACE;
      dat.getFieldData().resize(0, false);
      dat.getFieldDofs().resize(0, false);
    }
  }

  auto &dofs = const_cast<FEDofEntity_multiIndex &>(
      numeredEntFiniteElementPtr->getDataDofs());
  auto &dofs_by_type = dofs.get<Composite_Name_And_Type_mi_tag>();
  auto dit = dofs_by_type.lower_bound(boost::make_tuple(field_name, type_lo));
  if (dit == dofs_by_type.end())
    MoFEMFunctionReturnHot(0);
  auto hi_dit =
      dofs_by_type.lower_bound(boost::make_tuple(field_name, type_hi));
  std::vector<boost::weak_ptr<FEDofEntity>> brother_dofs_vec;
  for (; dit != hi_dit;) {

    auto &dof = **dit;
    const int nb_dofs_on_ent = dof.getNbDofsOnEnt();
    if (nb_dofs_on_ent) {

      const EntityType type = dof.getEntType();
      // must be non const since side has to change since it must be renumbered
      int side = dof.sideNumberPtr->side_number;

      if ((master_flag == 1 && type == MBEDGE && side > 2) ||
          (master_flag == 1 && type == MBTRI && side != 3)) {
        for (int ii = 0; ii != nb_dofs_on_ent; ++ii) {
          ++dit;
        }
        continue;
      } else if ((master_flag == 0 && type == MBEDGE && side < 6) ||
                 (master_flag == 0 && type == MBTRI && side != 4)) {
        for (int ii = 0; ii != nb_dofs_on_ent; ++ii) {
          ++dit;
        }
        continue;
      }

      if (type == MBTRI) {
        side = 0;
      }

      if (master_flag == 0 && type == MBEDGE) {
        side = side - 6;
      }

      auto &dat = data.dataOnEntities[type][side];

      auto &ent_field_dofs = dat.getFieldDofs();
      auto &ent_field_data = dat.getFieldData();
      const int brother_side = dof.sideNumberPtr->brother_side_number;
      if (brother_side != -1)
        brother_dofs_vec.emplace_back(*dit);

      if (ent_field_data.empty()) {

        dat.getBase() = dof.getApproxBase();
        dat.getSpace() = dof.getSpace();
        const int ent_order = dof.getMaxOrder();
        dat.getDataOrder() =
            dat.getDataOrder() > ent_order ? dat.getDataOrder() : ent_order;
        const auto dof_ent_field_data = dof.getEntFieldData();
        ent_field_data.resize(nb_dofs_on_ent, false);
        noalias(ent_field_data) = dof.getEntFieldData();
        ent_field_dofs.resize(nb_dofs_on_ent, false);
        for (int ii = 0; ii != nb_dofs_on_ent; ++ii) {
          ent_field_dofs[ii] = *dit;
          ++dit;
        }
      }
    }
  }

  for (auto &dof_ptr : brother_dofs_vec) {
    if (auto d = dof_ptr.lock()) {
      const EntityType type = d->getEntType();
      const int side = d->sideNumberPtr->side_number;
      const int brother_side = d->sideNumberPtr->brother_side_number;
      auto &dat = data.dataOnEntities[type][side];
      auto &dat_brother = data.dataOnEntities[type][brother_side];
      dat_brother.getBase() = dat.getBase();
      dat_brother.getSpace() = dat.getSpace();
      dat_brother.getDataOrder() = dat.getDataOrder();
      dat_brother.getFieldData() = dat.getFieldData();
      dat_brother.getFieldDofs() = dat.getFieldDofs();
    }
  }

  MoFEMFunctionReturn(0);
}


// ** Indices **

MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::getNodesFieldData(
    DataForcesAndSourcesCore &data, const std::string &field_name,
    const bool &master_flag) const {
  return getNodesFieldData(field_name,
                           const_cast<FEDofEntity_multiIndex &>(
                               numeredEntFiniteElementPtr->getDataDofs()),
                           data.dataOnEntities[MBVERTEX][0].getFieldData(),
                           data.dataOnEntities[MBVERTEX][0].getFieldDofs(),
                           data.dataOnEntities[MBVERTEX][0].getSpace(),
                           data.dataOnEntities[MBVERTEX][0].getBase(),
                           master_flag);
}

// ** Data **

MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::getNodesFieldData(
    const boost::string_ref field_name, FEDofEntity_multiIndex &dofs,
    VectorDouble &nodes_data, VectorDofs &nodes_dofs, FieldSpace &space,
    FieldApproximationBase &base, const bool &master_flag) const {
  MoFEMFunctionBegin;
  auto &dofs_by_name_and_type = dofs.get<Composite_Name_And_Type_mi_tag>();
  auto tuple = boost::make_tuple(field_name, MBVERTEX);
  auto dit = dofs_by_name_and_type.lower_bound(tuple);
  if (dit == dofs_by_name_and_type.end())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "No nodal dofs on element");
  auto hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(tuple);

  if (dit != hi_dit) {
    auto &first_dof = **dit;
    space = first_dof.getSpace();
    base = first_dof.getApproxBase();
    constexpr int num_nodes = 3;
    const int nb_dof_idx = first_dof.getNbOfCoeffs();
    const int max_nb_dofs = nb_dof_idx * num_nodes;
    nodes_data.resize(max_nb_dofs, false);
    nodes_dofs.resize(max_nb_dofs, false);
    nodes_data.clear();

    std::vector<boost::weak_ptr<FEDofEntity>> brother_dofs_vec;
    for (; dit != hi_dit;) {
      const auto &dof_ptr = *dit;
      const auto &dof = *dof_ptr;
      const auto &sn = *dof.sideNumberPtr;
      int side_number = sn.side_number;
      const int brother_side_number = sn.brother_side_number;

      if (master_flag == 1 && side_number > 2) {
        for (int ii = 0; ii != nb_dof_idx; ++ii) {
          ++dit;
        }
        continue;
      } else if (master_flag == 0 && side_number < 3) {
        for (int ii = 0; ii != nb_dof_idx; ++ii) {
          ++dit;
        }
        continue;
      }

      if (master_flag == 0) {
        side_number = side_number - 3;
      }

      if (brother_side_number != -1)
        brother_dofs_vec.emplace_back(dof_ptr);

      int pos = side_number * nb_dof_idx;
      auto ent_filed_data_vec = dof.getEntFieldData();
      for (int ii = 0; ii != nb_dof_idx; ++ii) {
        nodes_data[pos] = ent_filed_data_vec[ii];
        nodes_dofs[pos] = *dit;
        ++pos;
        ++dit;
      }
    }

    for (auto &dof_ptr : brother_dofs_vec) {
      if (const auto d = dof_ptr.lock()) {
        const auto &sn = d->sideNumberPtr;
        int side_number = sn->side_number;
        const int brother_side_number = sn->brother_side_number;

        if (master_flag == 1 && side_number > 2) {
          continue;
        } else if (master_flag == 0 && side_number < 3) {
          continue;
        }

        if (master_flag == 0) {
          side_number = side_number - 3;
        }

        int pos = side_number * nb_dof_idx;
        int brother_pos = brother_side_number * nb_dof_idx;
        for (int ii = 0; ii != nb_dof_idx; ++ii) {
          nodes_data[brother_pos] = nodes_data[pos];
          nodes_dofs[brother_pos] = nodes_dofs[pos];
          ++pos;
          ++brother_pos;
        }
      }
    }

  } else {
    nodes_data.resize(0, false);
    nodes_dofs.resize(0, false);
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
