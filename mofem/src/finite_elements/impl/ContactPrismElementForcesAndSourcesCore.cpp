/** \file ContactPrismElementForcesAndSourcesCore.cpp

\brief Implementation of contact prism element

*/

namespace MoFEM {

ContactPrismElementForcesAndSourcesCore::
    ContactPrismElementForcesAndSourcesCore(Interface &m_field)
    : ForcesAndSourcesCore(m_field),
      dataOnMaster{

          nullptr,
          boost::make_shared<EntitiesFieldData>(MBENTITYSET), // NOFIELD
          boost::make_shared<EntitiesFieldData>(MBENTITYSET), // H1
          boost::make_shared<EntitiesFieldData>(MBENTITYSET), // HCURL
          boost::make_shared<EntitiesFieldData>(MBENTITYSET), // HDIV
          boost::make_shared<EntitiesFieldData>(MBENTITYSET)  // L2

      },
      dataOnSlave{

          nullptr,
          boost::make_shared<EntitiesFieldData>(MBENTITYSET), // NOFIELD
          boost::make_shared<EntitiesFieldData>(MBENTITYSET), // H1
          boost::make_shared<EntitiesFieldData>(MBENTITYSET), // HCURL
          boost::make_shared<EntitiesFieldData>(MBENTITYSET), // HDIV
          boost::make_shared<EntitiesFieldData>(MBENTITYSET)  // L2

      },
      derivedDataOnMaster{

          nullptr,
          boost::make_shared<DerivedEntitiesFieldData>(dataOnMaster[NOFIELD]),
          boost::make_shared<DerivedEntitiesFieldData>(dataOnMaster[H1]),
          boost::make_shared<DerivedEntitiesFieldData>(dataOnMaster[HCURL]),
          boost::make_shared<DerivedEntitiesFieldData>(dataOnMaster[HDIV]),
          boost::make_shared<DerivedEntitiesFieldData>(dataOnMaster[L2])

      },
      derivedDataOnSlave{

          nullptr,
          boost::make_shared<DerivedEntitiesFieldData>(dataOnSlave[NOFIELD]),
          boost::make_shared<DerivedEntitiesFieldData>(dataOnSlave[H1]),
          boost::make_shared<DerivedEntitiesFieldData>(dataOnSlave[HCURL]),
          boost::make_shared<DerivedEntitiesFieldData>(dataOnSlave[HDIV]),
          boost::make_shared<DerivedEntitiesFieldData>(dataOnSlave[L2])

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
      dataL2Slave(*dataOnSlave[L2].get()), opContravariantTransform() {

  getUserPolynomialBase() =
      boost::shared_ptr<BaseFunction>(new TriPolynomialBase());

  // Data on elements for proper spaces
  dataOnMaster[H1]->setElementType(MBTRI);
  derivedDataOnMaster[H1]->setElementType(MBTRI);
  dataOnSlave[H1]->setElementType(MBTRI);
  derivedDataOnSlave[H1]->setElementType(MBTRI);

  dataOnMaster[HDIV]->setElementType(MBTRI);
  derivedDataOnMaster[HDIV]->setElementType(MBTRI);
  dataOnSlave[HDIV]->setElementType(MBTRI);
  derivedDataOnSlave[HDIV]->setElementType(MBTRI);

  dataOnMaster[NOFIELD]->dataOnEntities[MBENTITYSET].push_back(
      new EntitiesFieldData::EntData());
  dataOnSlave[NOFIELD]->dataOnEntities[MBENTITYSET].push_back(
      new EntitiesFieldData::EntData());

  derivedDataOnMaster[NOFIELD]->dataOnEntities[MBENTITYSET].push_back(
      new EntitiesFieldData::EntData());
  derivedDataOnSlave[NOFIELD]->dataOnEntities[MBENTITYSET].push_back(
      new EntitiesFieldData::EntData());

  CHK_THROW_MESSAGE(createDataOnElement(MBPRISM),
                    "Problem with creation data on element");
}

MoFEMErrorCode
ContactPrismElementForcesAndSourcesCore::setDefaultGaussPts(const int rule) {
  MoFEMFunctionBegin;

  if (rule < QUAD_2D_TABLE_SIZE) {
    if (QUAD_2D_TABLE[rule]->dim != 2) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong dimension");
    }
    if (QUAD_2D_TABLE[rule]->order < rule) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "wrong order %d != %d", QUAD_2D_TABLE[rule]->order, rule);
    }
    nbGaussPts = QUAD_2D_TABLE[rule]->npoints;
    // For master and slave
    gaussPtsMaster.resize(3, nbGaussPts, false);
    gaussPtsSlave.resize(3, nbGaussPts, false);

    cblas_dcopy(nbGaussPts, &QUAD_2D_TABLE[rule]->points[1], 3,
                &gaussPtsMaster(0, 0), 1);
    cblas_dcopy(nbGaussPts, &QUAD_2D_TABLE[rule]->points[2], 3,
                &gaussPtsMaster(1, 0), 1);
    cblas_dcopy(nbGaussPts, QUAD_2D_TABLE[rule]->weights, 1,
                &gaussPtsMaster(2, 0), 1);

    gaussPtsSlave = gaussPtsMaster;

    dataH1Master.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nbGaussPts, 3,
                                                                 false);

    dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nbGaussPts, 3,
                                                                false);

    double *shape_ptr_master =
        &*dataH1Master.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
    cblas_dcopy(3 * nbGaussPts, QUAD_2D_TABLE[rule]->points, 1,
                shape_ptr_master, 1);
    double *shape_ptr_slave =
        &*dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
    cblas_dcopy(3 * nbGaussPts, QUAD_2D_TABLE[rule]->points, 1, shape_ptr_slave,
                1);
  } else {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "rule > quadrature order %d < %d", rule, QUAD_2D_TABLE_SIZE);
    nbGaussPts = 0;
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::operator()() {
  MoFEMFunctionBegin;

  if (numeredEntFiniteElementPtr->getEntType() != MBPRISM)
    MoFEMFunctionReturnHot(0);

  EntitiesFieldData &data_div = *dataOnElement[HDIV];
  EntitiesFieldData &data_curl = *dataOnElement[HCURL];
  EntitiesFieldData &data_l2 = *dataOnElement[HCURL];

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

    for (auto &v : {&tangentMasterOne, &tangentMasterTwo, &tangentSlaveOne,
                    &tangentSlaveTwo}) {
      v->resize(3);
      v->clear();
    }

    CHKERR Tools::getTriNormal(&coords[0], &normal[0]);
    CHKERR Tools::getTriNormal(&coords[9], &normal[3]);

    auto get_vec_ptr = [](VectorDouble &vec_double, int r = 0) {
      return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
          &vec_double(r + 0), &vec_double(r + 1), &vec_double(r + 2));
    };

    auto t_coords_master = get_vec_ptr(coords);
    auto t_coords_slave = get_vec_ptr(coords, 9);
    auto t_normal_master = get_vec_ptr(normal);
    auto t_normal_slave = get_vec_ptr(normal, 3);

    auto t_t1_master = get_vec_ptr(tangentMasterOne);
    auto t_t2_master = get_vec_ptr(tangentMasterTwo);
    auto t_t1_slave = get_vec_ptr(tangentSlaveOne);
    auto t_t2_slave = get_vec_ptr(tangentSlaveTwo);

    const double *diff_ptr = Tools::diffShapeFunMBTRI.data();
    FTensor::Tensor1<FTensor::PackPtr<const double *, 2>, 2> t_diff(
        &diff_ptr[0], &diff_ptr[1]);

    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;
    FTensor::Index<'k', 3> k;

    FTensor::Number<0> N0;

    for (int nn = 0; nn != 3; ++nn) {
      t_t1_master(i) += t_coords_master(i) * t_diff(N0);
      t_t1_slave(i) += t_coords_slave(i) * t_diff(N0);
      ++t_coords_master;
      ++t_coords_slave;
      ++t_diff;
    }

    aRea[0] = sqrt(t_normal_master(i) * t_normal_master(i)) / 2.;
    aRea[1] = sqrt(t_normal_slave(i) * t_normal_slave(i)) / 2.;

    t_t2_master(j) =
        FTensor::levi_civita(i, j, k) * t_normal_master(k) * t_t1_master(i);
    t_t2_slave(j) =
        FTensor::levi_civita(i, j, k) * t_normal_slave(k) * t_t1_slave(i);

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

  // Hcurl
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
    dataHdivSlave.spacesOnEntities[MBTRI].set(HDIV);
    dataHdivMaster.spacesOnEntities[MBTRI].set(HDIV);
    data_div.spacesOnEntities[MBTRI].set(HDIV);
  }

  // L2
  if ((dataH1.spacesOnEntities[MBTRI]).test(L2)) {
    CHKERR getEntitySense<MBTRI>(data_l2);
    CHKERR getEntityDataOrder<MBTRI>(data_l2, L2);
  }

  auto clean_data = [](EntitiesFieldData &data) {
    MoFEMFunctionBegin;
    data.bAse.reset();
    for (EntityType t = MBVERTEX; t != MBMAXTYPE; ++t) {
      data.spacesOnEntities[t].reset();
      data.basesOnEntities[t].reset();
    }
    for (int s = 0; s != LASTSPACE; ++s)
      data.basesOnSpaces[s].reset();

    MoFEMFunctionReturn(0);
  };

  auto copy_data = [](EntitiesFieldData &data, EntitiesFieldData &copy_data,
                      const int shift) {
    MoFEMFunctionBegin;

    if (shift != 0 && shift != 6) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong shift for contact prism element");
    }

    data.bAse = copy_data.bAse;
    for (auto t : {MBVERTEX, MBEDGE, MBTRI}) {
      data.spacesOnEntities[t] = copy_data.spacesOnEntities[t];
      data.basesOnEntities[t] = copy_data.basesOnEntities[t];
      data.basesOnSpaces[t] = copy_data.basesOnSpaces[t];
      data.brokenBasesOnSpaces[t] = copy_data.brokenBasesOnSpaces[t];
    }

    for (int ii = 0; ii != 3; ++ii) {
      data.dataOnEntities[MBEDGE][ii].getSense() =
          copy_data.dataOnEntities[MBEDGE][ii + shift].getSense();
      data.dataOnEntities[MBEDGE][ii].getOrder() =
          copy_data.dataOnEntities[MBEDGE][ii + shift].getOrder();
    }

    if (shift == 0) {
      data.dataOnEntities[MBTRI][0].getSense() =
          copy_data.dataOnEntities[MBTRI][3].getSense();
      data.dataOnEntities[MBTRI][0].getOrder() =
          copy_data.dataOnEntities[MBTRI][3].getOrder();
    } else {
      data.dataOnEntities[MBTRI][0].getSense() =
          copy_data.dataOnEntities[MBTRI][4].getSense();
      data.dataOnEntities[MBTRI][0].getOrder() =
          copy_data.dataOnEntities[MBTRI][4].getOrder();
    }

    MoFEMFunctionReturn(0);
  };

  auto copy_data_hdiv = [](EntitiesFieldData &data,
                           EntitiesFieldData &copy_data, const int shift) {
    MoFEMFunctionBegin;

    if (shift != 3 && shift != 4) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong shift for contact prism element");
    }

    data.bAse = copy_data.bAse;

    for (auto t : {MBVERTEX, MBTRI}) {
      data.spacesOnEntities[t] = copy_data.spacesOnEntities[MBVERTEX];
      data.basesOnEntities[t] = copy_data.basesOnEntities[MBVERTEX];
      data.basesOnSpaces[t] = copy_data.basesOnSpaces[MBVERTEX];
    }

    auto &cpy_ent_dat = copy_data.dataOnEntities[MBTRI][shift];
    auto &ent_dat = data.dataOnEntities[MBTRI][0];
    ent_dat.getBase() = cpy_ent_dat.getBase();
    ent_dat.getSpace() = cpy_ent_dat.getSpace();
    ent_dat.getSense() = ent_dat.getSense();
    ent_dat.getOrder() = cpy_ent_dat.getOrder();

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

  if (rule >= 0) {

    CHKERR setDefaultGaussPts(rule);

  } else {

    // Master-Slave
    if (gaussPtsMaster.size2() != gaussPtsSlave.size2())
      SETERRQ(
          PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
          "Number of Gauss Points at Master triangle is different than slave");

    CHKERR setGaussPts(order_row, order_col, order_data);
    nbGaussPts = gaussPtsMaster.size2();
    dataH1Master.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nbGaussPts, 3,
                                                                 false);
    dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nbGaussPts, 3,
                                                                false);

    if (nbGaussPts) {
      CHKERR Tools::shapeFunMBTRI<1>(&*dataH1Master.dataOnEntities[MBVERTEX][0]
                                           .getN(NOBASE)
                                           .data()
                                           .begin(),
                                     &gaussPtsMaster(0, 0),
                                     &gaussPtsMaster(1, 0), nbGaussPts);

      CHKERR Tools::shapeFunMBTRI<1>(
          &*dataH1Slave.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
          &gaussPtsSlave(0, 0), &gaussPtsSlave(1, 0), nbGaussPts);
    }
  }

  if (nbGaussPts == 0)
    MoFEMFunctionReturnHot(0);

  // Get coordinates on slave and master
  {
    coordsAtGaussPtsMaster.resize(nbGaussPts, 3, false);
    coordsAtGaussPtsSlave.resize(nbGaussPts, 3, false);
    for (int gg = 0; gg < nbGaussPts; gg++) {
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

      if (space == HDIV) {

        if (dataH1.spacesOnEntities[MBTRI].test(HDIV)) {

          CHKERR clean_data(dataHdivSlave);
          CHKERR copy_data_hdiv(dataHdivSlave, data_div, 4);
          CHKERR clean_data(dataHdivMaster);
          CHKERR copy_data_hdiv(dataHdivMaster, data_div, 3);

          dataHdivMaster.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) =
              dataOnMaster[H1]->dataOnEntities[MBVERTEX][0].getNSharedPtr(
                  NOBASE);
          dataHdivSlave.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) =
              dataOnSlave[H1]->dataOnEntities[MBVERTEX][0].getNSharedPtr(
                  NOBASE);
        }
      }
    }

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
    if (dataH1.bAse.test(b)) {
      switch (static_cast<FieldApproximationBase>(b)) {
      case AINSWORTH_LEGENDRE_BASE:
      case AINSWORTH_LOBATTO_BASE:
        if (dataH1.spacesOnEntities[MBVERTEX].test(H1)) {
          CHKERR getUserPolynomialBase()
              -> getValue(
                  gaussPtsMaster,
                  boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                      dataH1Master, H1, CONTINUOUS,
                      static_cast<FieldApproximationBase>(b), NOBASE)));

          CHKERR getUserPolynomialBase()
              -> getValue(
                  gaussPtsSlave,
                  boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                      dataH1Slave, H1, CONTINUOUS,
                      static_cast<FieldApproximationBase>(b), NOBASE)));
        }
        break;
      case DEMKOWICZ_JACOBI_BASE:
        if (dataH1.spacesOnEntities[MBTRI].test(HDIV)) {

          CHKERR getUserPolynomialBase()
              -> getValue(
                  gaussPtsMaster,
                  boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                      dataHdivMaster, HDIV, CONTINUOUS,
                      static_cast<FieldApproximationBase>(b), NOBASE)));
          CHKERR getUserPolynomialBase()
              -> getValue(
                  gaussPtsSlave,
                  boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                      dataHdivSlave, HDIV, CONTINUOUS,
                      static_cast<FieldApproximationBase>(b), NOBASE)));

          opContravariantTransform.normalRawPtr = &normal;
          opContravariantTransform.normalShift = 0;
          CHKERR opContravariantTransform.opRhs(dataHdivMaster);
          opContravariantTransform.normalShift = 3;
          CHKERR opContravariantTransform.opRhs(dataHdivSlave);
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

  // Iterate over operators
  CHKERR loopOverOperators();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::loopOverOperators() {
  MoFEMFunctionBegin;

  constexpr std::array<UserDataOperator::OpType, 2> types{
      UserDataOperator::OPROW, UserDataOperator::OPCOL};
  std::array<std::string, 2> last_eval_field_name{std::string(), std::string()};

  auto oit = opPtrVector.begin();
  auto hi_oit = opPtrVector.end();

  for (; oit != hi_oit; oit++) {

    oit->setPtrFE(this);

    if (oit->opType == UserDataOperator::OPSPACE) {

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
      boost::shared_ptr<EntitiesFieldData> op_master_data[2];
      boost::shared_ptr<EntitiesFieldData> op_slave_data[2];

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
                   field_name.c_str(), getFEName().c_str());

        if (oit->getOpType() & types[ss] ||
            oit->getOpType() & UserDataOperator::OPROWCOL) {

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

            if (!ss) {
              struct Extractor {
                boost::weak_ptr<EntityCacheNumeredDofs>
                operator()(boost::shared_ptr<FieldEntity> &e) {
                  return e->entityCacheRowDofs;
                }
              };

              CHKERR getEntityIndices(*op_master_data[ss], *op_slave_data[ss],
                                      field_name, getRowFieldEnts(), MBEDGE,
                                      MBPRISM, Extractor());
            } else {
              struct Extractor {
                boost::weak_ptr<EntityCacheNumeredDofs>
                operator()(boost::shared_ptr<FieldEntity> &e) {
                  return e->entityCacheColDofs;
                }
              };
              CHKERR getEntityIndices(*op_master_data[ss], *op_slave_data[ss],
                                      field_name, getRowFieldEnts(), MBEDGE,
                                      MBPRISM, Extractor());
            }

            switch (space) {
            case NOSPACE:
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                      "unknown space");
              break;
            case H1: {

              auto get_indices = [&](auto &master, auto &slave, auto &ents,
                                     auto &&ex) {
                return getNodesIndices(
                    field_name, ents,
                    master.dataOnEntities[MBVERTEX][0].getIndices(),
                    master.dataOnEntities[MBVERTEX][0].getLocalIndices(),
                    slave.dataOnEntities[MBVERTEX][0].getIndices(),
                    slave.dataOnEntities[MBVERTEX][0].getLocalIndices(), ex);
              };

              auto get_data = [&](EntitiesFieldData &master_data,
                                  EntitiesFieldData &slave_data) {
                return getNodesFieldData(
                    field_name,
                    master_data.dataOnEntities[MBVERTEX][0].getFieldData(),
                    slave_data.dataOnEntities[MBVERTEX][0].getFieldData(),
                    master_data.dataOnEntities[MBVERTEX][0].getFieldDofs(),
                    slave_data.dataOnEntities[MBVERTEX][0].getFieldDofs(),
                    master_data.dataOnEntities[MBVERTEX][0].getFieldEntities(),
                    slave_data.dataOnEntities[MBVERTEX][0].getFieldEntities(),
                    master_data.dataOnEntities[MBVERTEX][0].getSpace(),
                    slave_data.dataOnEntities[MBVERTEX][0].getSpace(),
                    master_data.dataOnEntities[MBVERTEX][0].getBase(),
                    slave_data.dataOnEntities[MBVERTEX][0].getBase());
              };

              if (!ss) {

                struct Extractor {
                  boost::weak_ptr<EntityCacheNumeredDofs>
                  operator()(boost::shared_ptr<FieldEntity> &e) {
                    return e->entityCacheRowDofs;
                  }
                };

                CHKERR get_indices(*op_master_data[ss], *op_slave_data[ss],
                                   getRowFieldEnts(), Extractor());
              } else {
                struct Extractor {
                  boost::weak_ptr<EntityCacheNumeredDofs>
                  operator()(boost::shared_ptr<FieldEntity> &e) {
                    return e->entityCacheColDofs;
                  }
                };

                CHKERR get_indices(*op_master_data[ss], *op_slave_data[ss],
                                   getColFieldEnts(), Extractor());
              }

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
                const auto bit_number = mField.get_field_bit_number(field_name);
                if (!ss) {
                  CHKERR getNoFieldRowIndices(*op_master_data[ss], bit_number);
                  CHKERR getNoFieldRowIndices(*op_slave_data[ss], bit_number);
                } else {
                  CHKERR getNoFieldColIndices(*op_master_data[ss], bit_number);
                  CHKERR getNoFieldColIndices(*op_slave_data[ss], bit_number);
                }
                CHKERR getNoFieldFieldData(*op_master_data[ss], bit_number);
                CHKERR getNoFieldFieldData(*op_slave_data[ss], bit_number);
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

        if ((oit->getOpType() & UserDataOperator::OPROWCOL) &&
            ((type & UserDataOperator::FACEMASTER) ||
             (type & UserDataOperator::FACESLAVE))) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Wrong combination of FaceType and OpType, OPROWCOL "
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
    EntitiesFieldData &master_data, EntitiesFieldData &slave_data,
    const std::string &field_name, const EntityType type_lo,
    const EntityType type_hi) const {
  MoFEMFunctionBegin;

  auto reset_data = [type_lo, type_hi](auto &data) {
    for (EntityType t = type_lo; t != type_hi; ++t) {
      for (auto &dat : data.dataOnEntities[t]) {
        dat.getOrder() = 0;
        dat.getBase() = NOBASE;
        dat.getSpace() = NOSPACE;
        dat.getFieldData().resize(0, false);
        dat.getFieldDofs().resize(0, false);
        dat.getFieldEntities().resize(0, false);
      }
    }
  };
  reset_data(master_data);
  reset_data(slave_data);

  auto &field_ents = getDataFieldEnts();
  auto bit_number = mField.get_field_bit_number(field_name);
  const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
      bit_number, get_id_for_min_type(type_lo));
  auto lo = std::lower_bound(field_ents.begin(), field_ents.end(), lo_uid,
                             cmp_uid_lo);
  if (lo != field_ents.end()) {
    const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
        bit_number, get_id_for_max_type(type_hi));
    auto hi = std::upper_bound(lo, field_ents.end(), hi_uid, cmp_uid_hi);
    if (lo != hi) {
      for (auto it = lo; it != hi; ++it)
        if (auto e = it->lock()) {

          auto get_data = [&](auto &data, auto type, auto side) {
            auto &dat = data.dataOnEntities[type][side];
            auto &ent_field_dofs = dat.getFieldDofs();
            auto &ent_field_data = dat.getFieldData();
            dat.getFieldEntities().resize(1, false);
            dat.getFieldEntities()[0] = e.get();
            dat.getBase() = e->getApproxBase();
            dat.getSpace() = e->getSpace();
            const int ent_order = e->getMaxOrder();
            dat.getOrder() =
                dat.getOrder() > ent_order ? dat.getOrder() : ent_order;
            const auto dof_ent_field_data = e->getEntFieldData();
            const int nb_dofs_on_ent = e->getNbDofsOnEnt();
            ent_field_data.resize(nb_dofs_on_ent, false);
            noalias(ent_field_data) = e->getEntFieldData();
            ent_field_dofs.resize(nb_dofs_on_ent, false);
            std::fill(ent_field_dofs.begin(), ent_field_dofs.end(), nullptr);
            if (auto cache = e->entityCacheDataDofs.lock()) {
              for (auto dit = cache->loHi[0]; dit != cache->loHi[1]; ++dit) {
                ent_field_dofs[(*dit)->getEntDofIdx()] =
                    reinterpret_cast<FEDofEntity *>((*dit).get());
              }
            }
          };

          const EntityType type = e->getEntType();
          const int side = e->getSideNumberPtr()->side_number;

          switch (type) {
          case MBEDGE:

            if (side < 3)
              get_data(master_data, type, side);
            else if (side > 5)
              get_data(slave_data, type, side - 6);

            break;
          case MBTRI:

            if (side == 3)
              get_data(master_data, type, 0);
            if (side == 4)
              get_data(slave_data, type, 0);

            break;
          default:
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "Entity type not implemented (FIXME)");
          };

          const int brother_side = e->getSideNumberPtr()->brother_side_number;
          if (brother_side != -1)
            SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                    "Case with brother side not implemented (FIXME)");
        }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::getNodesFieldData(
    const std::string field_name, VectorDouble &master_nodes_data,
    VectorDouble &slave_nodes_data, VectorDofs &master_nodes_dofs,
    VectorDofs &slave_nodes_dofs,

    VectorFieldEntities &master_field_entities,
    VectorFieldEntities &slave_field_entities,

    FieldSpace &master_space, FieldSpace &slave_space,
    FieldApproximationBase &master_base,
    FieldApproximationBase &slave_base) const {
  MoFEMFunctionBegin;

  auto set_zero = [](auto &nodes_data, auto &nodes_dofs, auto &field_entities) {
    nodes_data.resize(0, false);
    nodes_dofs.resize(0, false);
    field_entities.resize(0, false);
  };
  set_zero(master_nodes_data, master_nodes_dofs, master_field_entities);
  set_zero(slave_nodes_data, slave_nodes_dofs, slave_field_entities);

  auto field_it = fieldsPtr->get<FieldName_mi_tag>().find(field_name);
  if (field_it != fieldsPtr->get<FieldName_mi_tag>().end()) {

    auto bit_number = (*field_it)->getBitNumber();
    const int nb_dofs_on_vert = (*field_it)->getNbOfCoeffs();
    master_space = slave_space = (*field_it)->getSpace();
    master_base = slave_base = (*field_it)->getApproxBase();

    auto &field_ents = getDataFieldEnts();
    const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
        bit_number, get_id_for_min_type<MBVERTEX>());
    auto lo = std::lower_bound(field_ents.begin(), field_ents.end(), lo_uid,
                               cmp_uid_lo);
    if (lo != field_ents.end()) {
      const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
          bit_number, get_id_for_max_type<MBVERTEX>());
      auto hi = std::upper_bound(lo, field_ents.end(), hi_uid, cmp_uid_hi);
      if (lo != hi) {

        int nb_dofs = 0;
        for (auto it = lo; it != hi; ++it) {
          if (auto e = it->lock()) {
            nb_dofs += e->getEntFieldData().size();
          }
        }

        if (nb_dofs) {

          auto init_set = [&](auto &nodes_data, auto &nodes_dofs,
                              auto &field_entities) {
            constexpr int num_nodes = 3;
            const int max_nb_dofs = nb_dofs_on_vert * num_nodes;
            nodes_data.resize(max_nb_dofs, false);
            nodes_dofs.resize(max_nb_dofs, false);
            field_entities.resize(num_nodes, false);
            nodes_data.clear();
            fill(nodes_dofs.begin(), nodes_dofs.end(), nullptr);
            fill(field_entities.begin(), field_entities.end(), nullptr);
          };

          init_set(master_nodes_data, master_nodes_dofs, master_field_entities);
          init_set(slave_nodes_data, slave_nodes_dofs, slave_field_entities);

          for (auto it = lo; it != hi; ++it) {
            if (auto e = it->lock()) {

              const auto &sn = e->getSideNumberPtr();
              int side = sn->side_number;

              auto set_data = [&](auto &nodes_data, auto &nodes_dofs,
                                  auto &field_entities, int side, int pos) {
                field_entities[side] = e.get();
                if (auto cache = e->entityCacheDataDofs.lock()) {
                  for (auto dit = cache->loHi[0]; dit != cache->loHi[1];
                       ++dit) {
                    const auto dof_idx = (*dit)->getEntDofIdx();
                    nodes_data[pos + dof_idx] = (*dit)->getFieldData();
                    nodes_dofs[pos + dof_idx] =
                        reinterpret_cast<FEDofEntity *>((*dit).get());
                  }
                }
              };

              if (side < 3)
                set_data(master_nodes_data, master_nodes_dofs,
                         master_field_entities, side, side * nb_dofs_on_vert);
              else
                set_data(slave_nodes_data, slave_nodes_dofs,
                         slave_field_entities, (side - 3),
                         (side - 3) * nb_dofs_on_vert);

              const int brother_side = sn->brother_side_number;
              if (brother_side != -1)
                SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                        "Not implemented (FIXME please)");
            }
          }
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

template <typename EXTRACTOR>
MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::getEntityIndices(
    EntitiesFieldData &master_data, EntitiesFieldData &slave_data,
    const std::string &field_name, FieldEntity_vector_view &ents_field,
    const EntityType type_lo, const EntityType type_hi,
    EXTRACTOR &&extractor) const {
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

  auto bit_number = mField.get_field_bit_number(field_name);
  const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
      bit_number, get_id_for_min_type(type_lo));
  auto lo = std::lower_bound(ents_field.begin(), ents_field.end(), lo_uid,
                             cmp_uid_lo);
  if (lo != ents_field.end()) {
    const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
        bit_number, get_id_for_max_type(type_hi));
    auto hi = std::upper_bound(lo, ents_field.end(), hi_uid, cmp_uid_hi);
    if (lo != hi) {

      std::vector<boost::weak_ptr<FieldEntity>> brother_ents_vec;

      for (auto it = lo; it != hi; ++it)
        if (auto e = it->lock()) {

          const EntityType type = e->getEntType();
          const int side = e->getSideNumberPtr()->side_number;
          const int brother_side = e->getSideNumberPtr()->brother_side_number;
          if (brother_side != -1)
            SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                    "Not implemented case");

          auto get_indices = [&](auto &data, const auto type, const auto side) {
            if (auto cache = extractor(e).lock()) {
              for (auto dit = cache->loHi[0]; dit != cache->loHi[1]; ++dit) {
                auto &dof = (**dit);
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
              }
            }
          };

          switch (type) {
          case MBEDGE:

            if (side < 3)
              get_indices(master_data, type, side);
            else if (side > 5)
              get_indices(slave_data, type, side - 6);

            break;
          case MBTRI:

            if (side == 3)
              get_indices(master_data, type, 0);
            if (side == 4)
              get_indices(slave_data, type, 0);

            break;
          default:
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "Entity type not implemented");
          }
        }
    }
  }

  MoFEMFunctionReturn(0);
}

template <typename EXTRACTOR>
MoFEMErrorCode ContactPrismElementForcesAndSourcesCore::getNodesIndices(
    const std::string field_name, FieldEntity_vector_view &ents_field,
    VectorInt &master_nodes_indices, VectorInt &master_local_nodes_indices,
    VectorInt &slave_nodes_indices, VectorInt &slave_local_nodes_indices,
    EXTRACTOR &&extractor) const {
  MoFEMFunctionBegin;

  master_nodes_indices.resize(0, false);
  master_local_nodes_indices.resize(0, false);
  slave_nodes_indices.resize(0, false);
  slave_local_nodes_indices.resize(0, false);

  auto field_it = fieldsPtr->get<FieldName_mi_tag>().find(field_name);
  if (field_it != fieldsPtr->get<FieldName_mi_tag>().end()) {

    auto bit_number = (*field_it)->getBitNumber();
    const int nb_dofs_on_vert = (*field_it)->getNbOfCoeffs();

    const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
        bit_number, get_id_for_min_type<MBVERTEX>());
    auto lo = std::lower_bound(ents_field.begin(), ents_field.end(), lo_uid,
                               cmp_uid_lo);
    if (lo != ents_field.end()) {
      const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
          bit_number, get_id_for_max_type<MBVERTEX>());
      auto hi = std::upper_bound(lo, ents_field.end(), hi_uid, cmp_uid_hi);
      if (lo != hi) {

        int nb_dofs = 0;
        for (auto it = lo; it != hi; ++it) {
          if (auto e = it->lock()) {
            if (auto cache = extractor(e).lock()) {
              nb_dofs += std::distance(cache->loHi[0], cache->loHi[1]);
            }
          }
        }

        if (nb_dofs) {

          constexpr int num_nodes = 3;
          const int max_nb_dofs = nb_dofs_on_vert * num_nodes;

          auto set_vec_size = [&](auto &nodes_indices,
                                  auto &local_nodes_indices) {
            nodes_indices.resize(max_nb_dofs, false);
            local_nodes_indices.resize(max_nb_dofs, false);
            std::fill(nodes_indices.begin(), nodes_indices.end(), -1);
            std::fill(local_nodes_indices.begin(), local_nodes_indices.end(),
                      -1);
          };

          set_vec_size(master_nodes_indices, master_local_nodes_indices);
          set_vec_size(slave_nodes_indices, slave_local_nodes_indices);

          for (auto it = lo; it != hi; ++it) {
            if (auto e = it->lock()) {

              const int side = e->getSideNumberPtr()->side_number;
              const int brother_side =
                  e->getSideNumberPtr()->brother_side_number;
              if (brother_side != -1)
                SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                        "Not implemented case");

              auto get_indices = [&](auto &nodes_indices,
                                     auto &local_nodes_indices,
                                     const auto side) {
                if (auto cache = extractor(e).lock()) {
                  for (auto dit = cache->loHi[0]; dit != cache->loHi[1];
                       ++dit) {
                    const int idx = (*dit)->getPetscGlobalDofIdx();
                    const int local_idx = (*dit)->getPetscLocalDofIdx();
                    const int pos =
                        side * nb_dofs_on_vert + (*dit)->getDofCoeffIdx();
                    nodes_indices[pos] = idx;
                    local_nodes_indices[pos] = local_idx;
                  }
                }
              };

              if (side < 3)
                get_indices(master_nodes_indices, master_local_nodes_indices,
                            side);
              else if (side > 2)
                get_indices(slave_nodes_indices, slave_local_nodes_indices,
                            side - 3);
              else
                SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                        "Impossible case");
            }
          }
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ContactPrismElementForcesAndSourcesCore::UserDataOperator::loopSideVolumes(
    const string fe_name,
    VolumeElementForcesAndSourcesCoreOnContactPrismSide &fe_method,
    const int side_type, const EntityHandle ent_for_side) {
  return loopSide(fe_name, &fe_method, side_type, ent_for_side);
}

} // namespace MoFEM
