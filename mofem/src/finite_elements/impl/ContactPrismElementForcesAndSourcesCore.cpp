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
  constexpr std::array<UserDataOperator::OpType, 2> types{
      UserDataOperator::OPROW, UserDataOperator::OPCOL};
  std::array<std::string, 2> last_eval_field_name{std::string(), std::string()};

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

      // Run operator
      try {
        CHKERR oit->opRhs(*dataOnMaster[oit->sPace], false);
        CHKERR oit->opRhs(*dataOnSlave[oit->sPace], false);
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
                .none())
          SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "no data field < %s > on finite element < %s >",
                   field_name.c_str(), feName.c_str());

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
            CHKERR getEntityFieldData(*op_master_data[ss], *op_slave_data[ss],
                                      field_name, MBEDGE);

            if (!ss)
              CHKERR getEntityIndices(
                  *op_master_data[ss], *op_slave_data[ss], field_name,
                  const_cast<FENumeredDofEntity_multiIndex &>(
                      numeredEntFiniteElementPtr->getRowsDofs()),
                  MBEDGE);
            else
              CHKERR getEntityIndices(
                  *op_master_data[ss], *op_slave_data[ss], field_name,
                  const_cast<FENumeredDofEntity_multiIndex &>(
                      numeredEntFiniteElementPtr->getColsDofs()),
                  MBEDGE);

            switch (space) {
            case NOSPACE:
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                      "unknown space");
              break;
            case H1: {

              auto get_indices = [&](auto &master, auto &slave, auto &dofs) {
                return getNodesIndices(
                    field_name, dofs,
                    master.dataOnEntities[MBVERTEX][0].getIndices(),
                    master.dataOnEntities[MBVERTEX][0].getLocalIndices(),
                    slave.dataOnEntities[MBVERTEX][0].getIndices(),
                    slave.dataOnEntities[MBVERTEX][0].getLocalIndices());
              };

              auto get_data = [&](DataForcesAndSourcesCore &master_data,
                                  DataForcesAndSourcesCore &slave_data) {
                return getNodesFieldData(
                    field_name,
                    const_cast<FEDofEntity_multiIndex &>(
                        numeredEntFiniteElementPtr->getDataDofs()),
                    master_data.dataOnEntities[MBVERTEX][0].getFieldData(),
                    slave_data.dataOnEntities[MBVERTEX][0].getFieldData(),
                    master_data.dataOnEntities[MBVERTEX][0].getFieldDofs(),
                    slave_data.dataOnEntities[MBVERTEX][0].getFieldDofs(),
                    master_data.dataOnEntities[MBVERTEX][0].getSpace(),
                    slave_data.dataOnEntities[MBVERTEX][0].getSpace(),
                    master_data.dataOnEntities[MBVERTEX][0].getBase(),
                    slave_data.dataOnEntities[MBVERTEX][0].getBase());
              };

              if (!ss)
                CHKERR get_indices(
                    *op_master_data[ss], *op_slave_data[ss],
                    const_cast<FENumeredDofEntity_multiIndex &>(
                        numeredEntFiniteElementPtr->getRowsDofs()));
              else
                CHKERR get_indices(
                    *op_master_data[ss], *op_slave_data[ss],
                    const_cast<FENumeredDofEntity_multiIndex &>(
                        numeredEntFiniteElementPtr->getColsDofs()));

              CHKERR get_data(*op_master_data[ss], *op_slave_data[ss]);

            } break;
            case HCURL:
            case HDIV:
            case L2:
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
          CHKERR oit->opLhs(*op_master_data[0], *op_master_data[1]);
        }
        CATCH_OP_ERRORS(*oit);
      }

      if (oit->getOpType() & UserDataOperator::OPROWCOL &&
          (type & UserDataOperator::FACEMASTERSLAVE)) {
        try {
          CHKERR oit->opLhs(*op_master_data[0], *op_slave_data[1]);
        }
        CATCH_OP_ERRORS(*oit);
      }

      if (oit->getOpType() & UserDataOperator::OPROWCOL &&
          (type & UserDataOperator::FACESLAVEMASTER)) {
        try {
          CHKERR oit->opLhs(*op_slave_data[0], *op_master_data[1]);
        }
        CATCH_OP_ERRORS(*oit);
      }

      if (oit->getOpType() & UserDataOperator::OPROWCOL &&
          (type & UserDataOperator::FACESLAVESLAVE)) {
        try {
          CHKERR oit->opLhs(*op_slave_data[0], *op_slave_data[1]);
        }
        CATCH_OP_ERRORS(*oit);
      }
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::getEntityFieldData(
    DataForcesAndSourcesCore &master_data, DataForcesAndSourcesCore &slave_data,
    const std::string &field_name, const EntityType type_lo,
    const EntityType type_hi) const {
  MoFEMFunctionBegin;

  auto reset_data = [type_lo, type_hi](auto &data) {
    for (EntityType t = type_lo; t != type_hi; ++t) {
      for (auto &dat : data.dataOnEntities[t]) {
        dat.getDataOrder() = 0;
        dat.getBase() = NOBASE;
        dat.getSpace() = NOSPACE;
        dat.getFieldData().resize(0, false);
        dat.getFieldDofs().resize(0, false);
      }
    }
  };
  reset_data(master_data);
  reset_data(slave_data);

  auto &dofs = const_cast<FEDofEntity_multiIndex &>(
      numeredEntFiniteElementPtr->getDataDofs());
  auto &dofs_by_type = dofs.get<Composite_Name_And_Type_mi_tag>();
  auto dit = dofs_by_type.lower_bound(boost::make_tuple(field_name, type_lo));
  if (dit == dofs_by_type.end())
    MoFEMFunctionReturnHot(0);
  auto hi_dit =
      dofs_by_type.lower_bound(boost::make_tuple(field_name, type_hi));

  auto get_data = [&](auto &data, auto &dof, auto type, auto side) {
    auto &dat = data.dataOnEntities[type][side];
    auto &ent_field_dofs = dat.getFieldDofs();
    auto &ent_field_data = dat.getFieldData();
    if (ent_field_data.empty()) {
      dat.getBase() = dof.getApproxBase();
      dat.getSpace() = dof.getSpace();
      const int ent_order = dof.getMaxOrder();
      dat.getDataOrder() =
          dat.getDataOrder() > ent_order ? dat.getDataOrder() : ent_order;
      const auto dof_ent_field_data = dof.getEntFieldData();
      const int nb_dofs_on_ent = dof.getNbDofsOnEnt();
      ent_field_data.resize(nb_dofs_on_ent, false);
      noalias(ent_field_data) = dof.getEntFieldData();
      ent_field_dofs.resize(nb_dofs_on_ent, false);
      for (int ii = 0; ii != nb_dofs_on_ent; ++ii) {
        ent_field_dofs[ii] = *dit;
        ++dit;
      }
    }
  };

  for (; dit != hi_dit;) {

    auto &dof = **dit;
    if (dof.getNbDofsOnEnt()) {

      const EntityType type = dof.getEntType();
      const int side = dof.sideNumberPtr->side_number;

      switch (type) {
      case MBEDGE:

        if (side < 3)
          get_data(master_data, dof, type, side);
        else if (side > 5)
          get_data(slave_data, dof, type, side - 6);

        break;
      case MBTRI:

        if (side == 3)
          get_data(master_data, dof, type, 0);
        if (side == 4)
          get_data(slave_data, dof, type, 0);

        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Entity type not implemented");
      };

      const int brother_side = dof.sideNumberPtr->brother_side_number;
      if (brother_side != -1)
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                "Case with brother side not implemented");
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::getNodesFieldData(
    const boost::string_ref field_name, FEDofEntity_multiIndex &dofs,
    VectorDouble &master_nodes_data, VectorDouble &slave_nodes_data,
    VectorDofs &master_nodes_dofs, VectorDofs &slave_nodes_dofs,
    FieldSpace &master_space, FieldSpace &slave_space,
    FieldApproximationBase &master_base,
    FieldApproximationBase &slave_base) const {
  MoFEMFunctionBegin;

  auto set_zero = [](auto &nodes_data, auto &nodes_dofs) {
    nodes_data.resize(0, false);
    nodes_dofs.resize(0, false);
  };
  set_zero(master_nodes_data, master_nodes_dofs);
  set_zero(slave_nodes_data, slave_nodes_dofs);

  auto &dofs_by_name_and_type = dofs.get<Composite_Name_And_Type_mi_tag>();
  auto tuple = boost::make_tuple(field_name, MBVERTEX);
  auto dit = dofs_by_name_and_type.lower_bound(tuple);
  if (dit == dofs_by_name_and_type.end())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "No nodal dofs on element");
  auto hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(tuple);

  if (dit != hi_dit) {
    auto &first_dof = **dit;
    const int nb_dof_idx = first_dof.getNbOfCoeffs();

    auto init_set = [&](auto &nodes_data, auto &nodes_dofs, auto &space,
                        auto &base) {
      constexpr int num_nodes = 3;
      const int max_nb_dofs = nb_dof_idx * num_nodes;
      space = first_dof.getSpace();
      base = first_dof.getApproxBase();
      nodes_data.resize(max_nb_dofs, false);
      nodes_dofs.resize(max_nb_dofs, false);
      nodes_data.clear();
    };

    init_set(master_nodes_data, master_nodes_dofs, master_space, master_base);
    init_set(slave_nodes_data, slave_nodes_dofs, slave_space, slave_base);

    std::vector<boost::weak_ptr<FEDofEntity>> brother_dofs_vec;
    for (; dit != hi_dit;) {
      const auto &dof_ptr = *dit;
      const auto &dof = *dof_ptr;
      const auto &sn = *dof.sideNumberPtr;
      int side = sn.side_number;

      auto set_data = [&](auto &nodes_data, auto &nodes_dofs, int pos) {
        auto ent_filed_data_vec = dof.getEntFieldData();
        for (int ii = 0; ii != nb_dof_idx; ++ii) {
          nodes_data[pos] = ent_filed_data_vec[ii];
          nodes_dofs[pos] = *dit;
          ++pos;
          ++dit;
        }
      };

      if (side < 3)
        set_data(master_nodes_data, master_nodes_dofs, side * nb_dof_idx);
      else
        set_data(slave_nodes_data, slave_nodes_dofs, (side - 3) * nb_dof_idx);

      const int brother_side = sn.brother_side_number;
      if (brother_side != -1)
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::getEntityIndices(
    DataForcesAndSourcesCore &master_data, DataForcesAndSourcesCore &slave_data,
    const std::string &field_name, FENumeredDofEntity_multiIndex &dofs,
    const EntityType type_lo, const EntityType type_hi) const {
  MoFEMFunctionBegin;

  auto clear_data = [type_lo, type_hi](auto &data) {
    for (EntityType t = type_lo; t != type_hi; ++t) {
      for (auto &d : data.dataOnEntities[t]) {
        d.getIndices().resize(0, false);
        d.getLocalIndices().resize(0, false);
      }
    }
  };

  clear_data(master_data);
  clear_data(slave_data);

  auto &dofs_by_type = dofs.get<Composite_Name_And_Type_mi_tag>();
  auto dit = dofs_by_type.lower_bound(boost::make_tuple(field_name, type_lo));
  if (dit == dofs_by_type.end())
    MoFEMFunctionReturnHot(0);
  auto hi_dit =
      dofs_by_type.lower_bound(boost::make_tuple(field_name, type_hi));

  auto get_indices = [&](auto &data, auto &dof, const auto type,
                         const auto side) {
    auto &dat = data.dataOnEntities[type][side];
    auto &ent_field_indices = dat.getIndices();
    auto &ent_field_local_indices = dat.getLocalIndices();
    if (ent_field_indices.empty()) {
      const int nb_dofs_on_ent = dof.getNbDofsOnEnt();
      ent_field_indices.resize(nb_dofs_on_ent, false);
      ent_field_local_indices.resize(nb_dofs_on_ent, false);
      std::fill(ent_field_indices.data().begin(),
                ent_field_indices.data().end(), -1);
      std::fill(ent_field_local_indices.data().begin(),
                ent_field_local_indices.data().end(), -1);
    }
    const int idx = dof.getEntDofIdx();
    ent_field_indices[idx] = dof.getPetscGlobalDofIdx();
    ent_field_local_indices[idx] = dof.getPetscLocalDofIdx();
  };

  for (; dit != hi_dit; ++dit) {

    auto &dof = **dit;

    if (dof.getNbDofsOnEnt()) {

      const EntityType type = dof.getEntType();
      const int side = dof.sideNumberPtr->side_number;

      switch (type) {
      case MBEDGE:

        if (side < 3)
          get_indices(master_data, dof, type, side);
        else if (side > 5)
          get_indices(slave_data, dof, type, side - 6);

        break;
      case MBTRI:

        if (side == 3)
          get_indices(master_data, dof, type, 0);
        if (side == 4)
          get_indices(slave_data, dof, type, 0);

        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Entity type not implemented");
      };

      const int brother_side = dof.sideNumberPtr->brother_side_number;
      if (brother_side != -1)
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented case");
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::getNodesIndices(
    const boost::string_ref field_name, FENumeredDofEntity_multiIndex &dofs,
    VectorInt &master_nodes_indices, VectorInt &master_local_nodes_indices,
    VectorInt &slave_nodes_indices,
    VectorInt &slave_local_nodes_indices) const {
  MoFEMFunctionBegin;

  auto &dofs_by_type = dofs.get<Composite_Name_And_Type_mi_tag>();
  auto tuple = boost::make_tuple(field_name, MBVERTEX);
  auto dit = dofs_by_type.lower_bound(tuple);
  auto hi_dit = dofs_by_type.upper_bound(tuple);

  master_nodes_indices.resize(0, false);
  master_local_nodes_indices.resize(0, false);
  slave_nodes_indices.resize(0, false);
  slave_local_nodes_indices.resize(0, false);

  if (dit != hi_dit) {

    const int num_nodes = 3;
    int max_nb_dofs = 0;
    const int nb_dofs_on_vert = (*dit)->getNbOfCoeffs();
    max_nb_dofs = nb_dofs_on_vert * num_nodes;

    auto set_vec_size = [&](auto &nodes_indices, auto &local_nodes_indices) {
      nodes_indices.resize(max_nb_dofs, false);
      local_nodes_indices.resize(max_nb_dofs, false);
      if (std::distance(dit, hi_dit) != max_nb_dofs) {
        std::fill(nodes_indices.begin(), nodes_indices.end(), -1);
        std::fill(local_nodes_indices.begin(), local_nodes_indices.end(), -1);
      }
    };

    set_vec_size(master_nodes_indices, master_local_nodes_indices);
    set_vec_size(slave_nodes_indices, slave_local_nodes_indices);

    auto get_indices = [&](auto &nodes_indices, auto &local_nodes_indices,
                           auto &dof, const auto side) {
      const int idx = dof.getPetscGlobalDofIdx();
      const int local_idx = dof.getPetscLocalDofIdx();
      const int pos = side * nb_dofs_on_vert + dof.getDofCoeffIdx();
      nodes_indices[pos] = idx;
      local_nodes_indices[pos] = local_idx;
    };

    for (; dit != hi_dit; dit++) {
      auto &dof = **dit;
      const int side = dof.sideNumberPtr->side_number;

      if (side < 3)
        get_indices(master_nodes_indices, master_local_nodes_indices, dof,
                    side);
      else if (side > 2)
        get_indices(slave_nodes_indices, slave_local_nodes_indices, dof,
                    side - 3);
      else
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Impossible case");

      const int brother_side = (*dit)->sideNumberPtr->brother_side_number;
      if (brother_side != -1)
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented case");
    }
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
