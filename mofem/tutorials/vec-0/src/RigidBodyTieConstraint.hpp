/**
 * @file RigidBodyTieConstraint.hpp
 * @brief Tie Constraint Implementation
 * @date 2024-11-14
 *
 * @copyright Copyright (c) 2024
 *
 */

// template <int DIM>
struct RigidBodyTieConstraintData {

  RigidBodyTieConstraintData(MoFEM::Interface &m_field) : mField(m_field) {}

  struct TieBlock {
    Range tieFaces;                                 //< Range of faces in block
    FTensor::Tensor1<double, 3> tieCoord;           //< Reference coordinate
    FTensor::Tensor1<double, 3> tieDirection;       //< Translation values
    FTensor::Tensor1<double, 3> tieRotation;        //< Rotation values
    FTensor::Tensor1<double, 3> tieTranslationFlag; //
    FTensor::Tensor1<double, 3> tieRotationFlag;    //
  };

  MoFEMErrorCode getTieBlocks(std::vector<TieBlock> &tieBlocks);

  MoFEM::Interface &mField;
  // std::vector<TieBlock> tieBlocks; //< Store all tie blocks
};

MoFEMErrorCode
RigidBodyTieConstraintData::getTieBlocks(std::vector<TieBlock> &tieBlocks) {
  // current
  MoFEMFunctionBegin;
  Range tie_ents;
  for (auto m :
       mField.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

           (boost::format("%s(.*)") % "TIE_MATRIX").str()

               ))

  ) {
    auto meshset = m->getMeshset();
    Range tie_meshset_range;
    CHKERR mField.get_moab().get_entities_by_dimension(meshset, SPACE_DIM - 1,
                                                       tie_meshset_range, true);
    std::vector<double> attributes;
    CHKERR m->getAttributes(attributes);
    if (attributes.size() != 15) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
               "Wrong number of tie parameters %d", attributes.size());
    }
    tieBlocks.push_back(
        {tie_meshset_range,
         FTensor::Tensor1<double, 3>(attributes[0], attributes[1],
                                     attributes[2]),
         FTensor::Tensor1<double, 3>(attributes[3], attributes[4],
                                     attributes[5]),
         FTensor::Tensor1<double, 3>(attributes[6], attributes[7],
                                     attributes[8]),
         FTensor::Tensor1<double, 3>(attributes[9], attributes[10],
                                     attributes[11]),
         FTensor::Tensor1<double, 3>(attributes[12], attributes[13],
                                     attributes[14])});
    MOFEM_LOG("SYNC", Sev::inform) << "Nb tie blocks " << tieBlocks.size();
    MOFEM_LOG("SYNC", Sev::inform) << "Nb tie faces " << tie_meshset_range.size();
    MOFEM_LOG("SYNC", Sev::inform)
        << "Tie reference position " << attributes[0] << " " << attributes[1]
        << " " << attributes[2];
    MOFEM_LOG("SYNC", Sev::inform)
        << "Tie translation " << attributes[3] << " " << attributes[4] << " "
        << attributes[5];
    MOFEM_LOG("SYNC", Sev::inform) << "Tie rotation " << attributes[6] << " "
                                    << attributes[7] << " " << attributes[8];
    MOFEM_LOG("SYNC", Sev::inform)
        << "Tie translation flag " << attributes[9] << " " << attributes[10]
        << " " << attributes[11];
    MOFEM_LOG("SYNC", Sev::inform)
        << "Tie rotation flag " << attributes[12] << " " << attributes[13]
        << " " << attributes[14];
    MOFEM_LOG("SYNC", Sev::inform) << "Tie meshset added";
    MOFEM_LOG_SYNCHRONISE(mField.get_comm());
  }
  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType AT>
struct OpTieTermConstrainRigidBodyRhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase {

  using OpBase = typename FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase;

  OpTieTermConstrainRigidBodyRhs(
      std::string lambda_name, boost::shared_ptr<MatrixDouble> u_ptr,
      FTensor::Tensor1<double, 3> tie_coord,
      boost::shared_ptr<Range> tie_faces_ptr,
      boost::shared_ptr<VectorDouble> translation_ptr,
      boost::shared_ptr<VectorDouble> theta_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), translationPtr(translation_ptr),
        thetaPtr(theta_ptr) {
    std::fill(&(OpBase::doEntities[MBEDGE]), &(OpBase::doEntities[MBMAXTYPE]), false);
    OpBase::doEntities[MBEDGE] = true;
    OpBase::doEntities[MBTRI] = true;
    OpBase::doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);

    auto nb_integration_pts = OpBase::getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function gradient on rows
    auto t_row_base = data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get coordinate at integration points
    auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
    // get displacement
    auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);

    FTensor::Tensor1<double, 3> t_translation;
    t_translation(0) = (*translationPtr)(0);
    t_translation(1) = (*translationPtr)(1);
    t_translation(2) = (*translationPtr)(2);

    FTensor::Tensor1<double, 3> t_theta;
    t_theta(0) = (*thetaPtr)(0);
    t_theta(1) = (*thetaPtr)(1);
    t_theta(2) = (*thetaPtr)(2);

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      // define rotation matrix
      FTensor::Tensor2<double, 3, 3> t_omega;
      t_omega(i, j) = levi_civita(i, j, k) * t_theta(k);

      // define reference position
      FTensor::Tensor1<double, 3> t_x_ref;
      t_x_ref(i) = t_u(i) - t_translation(i) -
                   t_omega(i, j) * (t_coords(j) - tieCoord(j));

      FTensor::Tensor1<double, SPACE_DIM> g;
      g(i) = alpha * (t_x_ref(i));

      ++t_coords;
      ++t_u;

      auto t_nf = getFTensor1FromArray<DIM, DIM>(OpBase::locF);

      int rr = 0;
      for (; rr != OpBase::nbRows / DIM; ++rr) {
        t_nf(i) += t_row_base * g(i);
        ++t_row_base;
        ++t_nf;
      }
      for (; rr < OpBase::nbRowBaseFunctions; ++rr)
        ++t_row_base;
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> uPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  boost::shared_ptr<VectorDouble> translationPtr;
  boost::shared_ptr<VectorDouble> thetaPtr;
};
template <int DIM, AssemblyType AT>
struct OpTieTermConstrainRigidBodyLhs_dU
    : public FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase {
  using OpBase = typename FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase;

  OpTieTermConstrainRigidBodyLhs_dU(
      std::string lambda_name, std::string col_field_name,
      boost::shared_ptr<MatrixDouble> u_ptr,
      FTensor::Tensor1<double, 3> tie_coord,
      boost::shared_ptr<Range> tie_faces_ptr,
      boost::shared_ptr<VectorDouble> translation_ptr,
      boost::shared_ptr<VectorDouble> theta_ptr)
      : OpBase(lambda_name, col_field_name, OpBase::OPROWCOL, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), translationPtr(translation_ptr),
        thetaPtr(theta_ptr) {
    std::fill(&(OpBase::doEntities[MBEDGE]), &(OpBase::doEntities[MBMAXTYPE]), false);
    OpBase::doEntities[MBEDGE] = true;
    OpBase::doEntities[MBTRI] = true;
    OpBase::doEntities[MBQUAD] = true;
    this->assembleTranspose = true;
    this->sYmm = false;
  }

  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    auto nb_integration_pts = OpBase::getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function gradient on rows
    auto t_row_base = row_data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get coordinate at integration points
    auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
    // get displacement
    auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);

    FTensor::Tensor1<double, 3> t_translation;
    t_translation(0) = (*translationPtr)(0);
    t_translation(1) = (*translationPtr)(1);
    t_translation(2) = (*translationPtr)(2);

    FTensor::Tensor1<double, 3> t_theta;
    t_theta(0) = (*thetaPtr)(0);
    t_theta(1) = (*thetaPtr)(1);
    t_theta(2) = (*thetaPtr)(2);

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
      t_tangent(i, j) = alpha * t_kd(i, j);

      int rr = 0;
      for (; rr != OpBase::nbRows / SPACE_DIM; ++rr) {
        auto t_col_base = col_data.getFTensor0N(gg, 0);
        auto t_mat =
            getFTensor2FromArray<DIM, DIM, DIM>(OpBase::locMat, DIM * rr);
        for (int cc = 0; cc != OpBase::nbCols / SPACE_DIM; cc++) {
          t_mat(i, j) += (t_row_base * t_col_base) * t_tangent(i, j);
          ++t_mat;
          ++t_col_base;
        }
        ++t_row_base;
      }
      for (; rr < OpBase::nbRowBaseFunctions; ++rr)
        ++t_row_base;
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> uPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  boost::shared_ptr<VectorDouble> translationPtr;
  boost::shared_ptr<VectorDouble> thetaPtr;
};

template <int DIM, AssemblyType AT>
struct OpTieTermConstrainRigidBodyLhs_dTranslation
    : public FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase {
  using OpBase = typename FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase;

  OpTieTermConstrainRigidBodyLhs_dTranslation(
      std::string lambda_name, std::string col_field_name,
      boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, col_field_name, OpBase::OPROWCOL, tie_faces_ptr) {
    std::fill(&(OpBase::doEntities[MBEDGE]), &(OpBase::doEntities[MBMAXTYPE]), false);
    OpBase::doEntities[MBEDGE] = true;
    OpBase::doEntities[MBTRI] = true;
    OpBase::doEntities[MBQUAD] = true;
    this->assembleTranspose = true;
    this->sYmm = false;
  }

  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;

    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    auto nb_integration_pts = OpBase::getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function gradient on rows
    auto t_row_base = row_data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
      t_tangent(i, j) = alpha * t_kd(i, j);

      int rr = 0;
      for (; rr != OpBase::nbRows / SPACE_DIM; ++rr) {
        auto t_mat = getFTensor2FromPtr<SPACE_DIM, SPACE_DIM>(
            &OpBase::locMat(SPACE_DIM * rr, 0));
        t_mat(i, j) -= (t_row_base * t_tangent(i, j));
        ++t_row_base;
      }
      for (; rr < OpBase::nbRowBaseFunctions; ++rr)
        ++t_row_base;
    }
    MoFEMFunctionReturn(0);
  }

private:
};

template <int DIM, AssemblyType AT>
struct OpTieTermConstrainRigidBodyLhs_dRotation
    : public FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase {
  using OpBase = typename FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase;

  OpTieTermConstrainRigidBodyLhs_dRotation(
      std::string lambda_name, std::string col_field_name,
      FTensor::Tensor1<double, 3> tie_coord,
      boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, col_field_name, OpBase::OPROWCOL, tie_faces_ptr),
        tieCoord(tie_coord) {
    std::fill(&(OpBase::doEntities[MBEDGE]), &(OpBase::doEntities[MBMAXTYPE]), false);
    OpBase::doEntities[MBEDGE] = true;
    OpBase::doEntities[MBTRI] = true;
    OpBase::doEntities[MBQUAD] = true;
    // doEntities[MBENTITYSET] = true;
    this->assembleTranspose = true;
    this->sYmm = false;
  }

  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    auto nb_integration_pts = OpBase::getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function gradient on rows
    auto t_row_base = row_data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get coordinate at integration points
    auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
      t_tangent(i, k) =
          alpha * levi_civita(i, j, k) * (t_coords(j) - tieCoord(j));
      ++t_coords;

      int rr = 0;
      for (; rr != OpBase::nbRows / SPACE_DIM; ++rr) {
        auto t_mat = getFTensor2FromPtr<SPACE_DIM, SPACE_DIM>(
            &OpBase::locMat(SPACE_DIM * rr, 0));
        t_mat(i, j) -= (t_row_base * t_tangent(i, j));
        ++t_row_base;
      }
      for (; rr < OpBase::nbRowBaseFunctions; ++rr)
        ++t_row_base;
    }
    MoFEMFunctionReturn(0);
  }

private:
  FTensor::Tensor1<double, 3> tieCoord;
};

template <int DIM, AssemblyType AT>
struct OpTieTermConstrainRigidBodyRhs_du
    : public FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase {

  using OpBase = typename FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase;

  OpTieTermConstrainRigidBodyRhs_du(std::string field_name,
                                    boost::shared_ptr<MatrixDouble> lambda_ptr,
                                    boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(field_name, field_name, OpBase::OPROW, tie_faces_ptr),
        lambdaPtr(lambda_ptr) {
    std::fill(&(OpBase::doEntities[MBEDGE]), &(OpBase::doEntities[MBMAXTYPE]), false);
    OpBase::doEntities[MBEDGE] = true;
    OpBase::doEntities[MBTRI] = true;
    OpBase::doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;

    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    auto nb_integration_pts = OpBase::getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function on rows
    auto t_row_base = data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get lambda
    auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

    auto &nf = OpBase::locF;

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
      t_tangent(i, j) = t_kd(i, j);

      auto t_nf = getFTensor1FromPtr<DIM>(&nf[0]);
      FTensor::Tensor1<double, SPACE_DIM> g;
      g(i) = alpha * (t_tangent(j, i) * t_lambda(j));
      ++t_lambda;

      int rr = 0;
      for (; rr != OpBase::nbRows / DIM; ++rr) {
        t_nf(i) += t_row_base * g(i);
        ++t_row_base;
        ++t_nf;
      }

      for (; rr < OpBase::nbRowBaseFunctions; ++rr)
        ++t_row_base;
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> lambdaPtr;
};

template <int DIM>
struct OpTieTermConstrainRigidBodyGlobalTranslationRhs
    : public ForcesAndSourcesCore::UserDataOperator {

  using OpUserDataOp = ForcesAndSourcesCore::UserDataOperator;

  OpTieTermConstrainRigidBodyGlobalTranslationRhs(
      std::string lambda_name,
      boost::shared_ptr<VectorDouble> int_translation_ptr, boost::shared_ptr<Range> rigid_body_ents_ptr)
      : OpUserDataOp(lambda_name, OpUserDataOp::OPROW),
        intTranslationPtr(int_translation_ptr), rigidBodyEntsPtr(rigid_body_ents_ptr) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBENTITYSET] = true;
  }
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;

    if (data.getIndices().empty() || rigidBodyEntsPtr->empty())
      MoFEMFunctionReturnHot(0);

    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);

    auto t_dofs = getFTensor1FromPtr<SPACE_DIM>(
        &data.getFieldEntities()[0]->getEntFieldData()[0]);

    double *intTranslationRawPtr = &(*intTranslationPtr->data().begin());

    if (intTranslationPtr->size() == 0)
      intTranslationPtr->resize(DIM);
    
    auto t_int_translation =
        getFTensor1FromPtr<SPACE_DIM>(intTranslationRawPtr);

    // create a "dummy" base function
    MatrixDouble m_kd(3, 3);
    m_kd(0, 0) = 1.0;
    m_kd(0, 1) = 0.0;
    m_kd(0, 2) = 0.0;
    m_kd(1, 0) = 0.0;
    m_kd(1, 1) = 1.0;
    m_kd(1, 2) = 0.0;
    m_kd(2, 0) = 0.0;
    m_kd(2, 1) = 0.0;
    m_kd(2, 2) = 1.0;

    auto t_dummy_base = getFTensor1FromMat<3>(m_kd);
    int rr = 0;
    for (; rr != 3; ++rr) {
      t_dofs(j) -= (t_dummy_base(i) * t_int_translation(i)) * t_dummy_base(j);
      ++t_dummy_base;
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<VectorDouble> intTranslationPtr;
  boost::shared_ptr<Range> rigidBodyEntsPtr;
};

template <int DIM, AssemblyType AT>
struct OpTieTermConstrainRigidBodyGlobalTranslationIntegralRhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase {

  using OpBase = typename FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase;

  OpTieTermConstrainRigidBodyGlobalTranslationIntegralRhs(
      std::string lambda_name, boost::shared_ptr<Range> tie_faces_ptr,
      boost::shared_ptr<MatrixDouble> lambda_ptr,
      boost::shared_ptr<VectorDouble> int_translation_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        lambdaPtr(lambda_ptr), intTranslationPtr(int_translation_ptr) {
    std::fill(&(OpBase::doEntities[MBEDGE]), &(OpBase::doEntities[MBMAXTYPE]), false);
    OpBase::doEntities[MBEDGE] = true;
    OpBase::doEntities[MBTRI] = true;
    OpBase::doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;

    if (intTranslationPtr->size() == 0)
      intTranslationPtr->resize(DIM);

    FTENSOR_INDEX(SPACE_DIM, i);

    auto nb_integration_pts = OpBase::getGaussPts().size2();

    // get element volume
    const double vol = OpBase::getMeasure();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get lambda
    auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

    double *intTranslationRawPtr = &(*intTranslationPtr->data().begin());

    auto t_int_translation = getFTensor1FromPtr<3>(intTranslationRawPtr);

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor1<double, 3> g;
      g(i) = alpha * t_lambda(i);
      ++t_lambda;

      t_int_translation(i) += g(i);
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<VectorDouble> intTranslationPtr;
  boost::shared_ptr<MatrixDouble> lambdaPtr;
};

template <int DIM, AssemblyType AT>
struct OpTieTermConstrainRigidBodyGlobalRotationIntegralRhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase {

  using OpBase = typename FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase;

  OpTieTermConstrainRigidBodyGlobalRotationIntegralRhs(
      std::string lambda_name, FTensor::Tensor1<double, 3> tie_coord,
      boost::shared_ptr<MatrixDouble> lambda_ptr,
      boost::shared_ptr<MatrixDouble> int_rotation_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW), tieCoord(tie_coord),
        lambdaPtr(lambda_ptr), intRotationPtr(int_rotation_ptr) {
    std::fill(&(OpBase::doEntities[MBEDGE]), &(OpBase::doEntities[MBMAXTYPE]), false);
    OpBase::doEntities[MBEDGE] = true;
    OpBase::doEntities[MBTRI] = true;
    OpBase::doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;

    if (intRotationPtr->size1() == 0 || intRotationPtr->size2() == 0)
      intRotationPtr->resize(DIM, DIM);

    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);

    auto nb_integration_pts = OpBase::getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function gradient on rows
    auto t_row_base = data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get coordinate at integration points
    auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
    // get lambda
    auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

    double *intRotationRawPtr = &(*intRotationPtr->data().begin());

    auto t_int_rotation =
        getFTensor2FromPtr<SPACE_DIM, SPACE_DIM>(intRotationRawPtr);

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> g;
      g(i, j) = alpha * t_lambda(i) * (t_coords(j) - tieCoord(j));
      ++t_lambda;
      ++t_coords;

      t_int_rotation(i, j) += g(i, j);
    }
    MoFEMFunctionReturn(0);
  }

private:
  FTensor::Tensor1<double, 3> tieCoord;
  boost::shared_ptr<MatrixDouble> lambdaPtr;
  boost::shared_ptr<MatrixDouble> intRotationPtr;
};

// template <int DIM>
// struct OpTieTermConstrainRigidBodyGlobalTranslationLhs_dLambda
//     : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

//   using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

//   OpTieTermConstrainRigidBodyGlobalTranslationLhs_dLambda(
//       std::string lambda_name, std::string field_name,
//        boost::shared_ptr<MatrixDouble> u_ptr,
//       FTensor::Tensor1<double, 3> tie_coord,
//       FTensor::Tensor1<double, 3> tie_direction,
//       boost::shared_ptr<Range> tie_faces_ptr)
//       : OpBase(lambda_name, field_name, OpBase::OPROWCOL, tie_faces_ptr),
//         uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction){
//     std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
//     // doEntities[MBEDGE] = true;
//     // doEntities[MBTRI] = true;
//     // doEntities[MBQUAD] = true;
//     // doEntities[MBENTITYSET] = true;
//     this->assembleTranspose = true;
//     this->sYmm = false;
//   }
//   MoFEMErrorCode iNtegrate(EntData &data) {
//   // MoFEMErrorCode doWork(int side, EntityType type,
//   //                       EntitiesFieldData::EntData &data) {
//     MoFEMFunctionBegin;
//     FTENSOR_INDEX(SPACE_DIM, i);
//     FTENSOR_INDEX(SPACE_DIM, j);
//     FTENSOR_INDEX(SPACE_DIM, k);
//     FTENSOR_INDEX(SPACE_DIM, l);
//     FTENSOR_INDEX(SPACE_DIM, m);

//     constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

//     // double time = getTStime();
//     double time = 1.0;

//     auto nb_integration_pts = getGaussPts().size2();

//     auto t_w = OpBase::getFTensor0IntegrationWeight();

//     const double vol = OpBase::getMeasure();

//     auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

//     // create a "dummy" base function
//     MatrixDouble m_kd(3, 3);
//     m_kd(0, 0) = 1.0;
//     m_kd(0, 1) = 0.0;
//     m_kd(0, 2) = 0.0;
//     m_kd(1, 0) = 0.0;
//     m_kd(1, 1) = 1.0;
//     m_kd(1, 2) = 0.0;
//     m_kd(2, 0) = 0.0;
//     m_kd(2, 1) = 0.0;
//     m_kd(2, 2) = 1.0;

//     for (int gg = 0; gg != nb_integration_pts; gg++) {
//       const auto alpha = t_w * vol;
//       ++t_w;

//       FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
//       t_tangent(i, j) = alpha * t_kd(i, j);

//       auto t_dummy_base = getFTensor1FromMat<3>(m_kd);

//       int rr = 0;
//       for (; rr != 3; ++rr) {
//         auto t_mat = getFTensor1FromPtr<SPACE_DIM>(&OpBase::locMat(rr, 0));
//         auto t_col_base = data.getFTensor0N(gg, 0);
//         for (int cc = 0; cc != OpBase::nbCols/SPACE_DIM; cc++) {
//           t_mat(j) += t_col_base * t_dummy_base(i) * t_tangent(i, j);
//           ++t_col_base;
//           ++t_mat;
//         }
//         ++t_dummy_base;
//       }
//     }
//     MoFEMFunctionReturn(0);
//   }

// private:
//   boost::shared_ptr<MatrixDouble> uPtr;
//   FTensor::Tensor1<double, 3> tieCoord;
//   FTensor::Tensor1<double, 3> tieDirection;
//   boost::shared_ptr<MatrixDouble> lambdaPtr;
// };

template <int DIM>
struct OpTieTermConstrainRigidBodyGlobalRotationRhs
    : public ForcesAndSourcesCore::UserDataOperator {

  using OpUserDataOp = ForcesAndSourcesCore::UserDataOperator;

  OpTieTermConstrainRigidBodyGlobalRotationRhs(
      std::string lambda_name, boost::shared_ptr<MatrixDouble> int_rotation_ptr, boost::shared_ptr<Range> rigid_body_ents_ptr)
      : OpUserDataOp(lambda_name, OpUserDataOp::OPROW),
        intRotationPtr(int_rotation_ptr), rigidBodyEntsPtr(rigid_body_ents_ptr) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBENTITYSET] = true;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {

    MoFEMFunctionBegin;

    if (rigidBodyEntsPtr->empty())
      MoFEMFunctionReturnHot(0);

    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);

    auto t_dofs = getFTensor1FromPtr<SPACE_DIM>(
        &data.getFieldEntities()[0]->getEntFieldData()[0]);

    double *intRotationRawPtr = &(*intRotationPtr->data().begin());
    auto t_int_rotation =
        getFTensor2FromPtr<SPACE_DIM, SPACE_DIM>(intRotationRawPtr);

    // create a "dummy" base function
    MatrixDouble m_kd(3, 3);
    m_kd(0, 0) = 1.0;
    m_kd(0, 1) = 0.0;
    m_kd(0, 2) = 0.0;
    m_kd(1, 0) = 0.0;
    m_kd(1, 1) = 1.0;
    m_kd(1, 2) = 0.0;
    m_kd(2, 0) = 0.0;
    m_kd(2, 1) = 0.0;
    m_kd(2, 2) = 1.0;

    auto t_dummy_base = getFTensor1FromMat<3>(m_kd);

    int rr = 0;
    for (; rr != 3; ++rr) {
      t_dofs(l) -=
          (t_dummy_base(k) * levi_civita(i, j, k) * t_int_rotation(i, j)) *
          t_dummy_base(l);
      ++t_dummy_base;
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> intRotationPtr;
  boost::shared_ptr<Range> rigidBodyEntsPtr;
};

// template <int DIM>
// struct OpTieTermConstrainRigidBodyGlobalRotationLhs_dLambda
//     : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

//   using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

//   OpTieTermConstrainRigidBodyGlobalRotationLhs_dLambda(
//       std::string lambda_name, std::string field_name,
//       boost::shared_ptr<MatrixDouble> u_ptr,
//       FTensor::Tensor1<double, 3> tie_coord,
//       FTensor::Tensor1<double, 3> tie_direction,
//       boost::shared_ptr<Range> tie_faces_ptr,
//       boost::shared_ptr<MatrixDouble> lambda_ptr,
//       boost::shared_ptr<VectorDouble> translation_ptr,
//       boost::shared_ptr<VectorDouble> theta_ptr)
//       : OpBase(lambda_name, field_name, OpBase::OPROWCOL, tie_faces_ptr),
//         uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction),
//         lambdaPtr(lambda_ptr), translationPtr(translation_ptr),
//         thetaPtr(theta_ptr) {
//     std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
//     // doEntities[MBEDGE] = true;
//     // doEntities[MBTRI] = true;
//     // doEntities[MBQUAD] = true;
//     // doEntities[MBENTITYSET] = true;
//     this->assembleTranspose = true;
//     this->sYmm = false;
//   }

//   // MoFEMErrorCode doWork(int side, EntityType type,
//   //                       EntitiesFieldData::EntData &data) {
//   MoFEMErrorCode iNtegrate(EntData &data) {
//     MoFEMFunctionBegin;
//     FTENSOR_INDEX(SPACE_DIM, i);
//     FTENSOR_INDEX(SPACE_DIM, j);
//     FTENSOR_INDEX(SPACE_DIM, k);
//     FTENSOR_INDEX(SPACE_DIM, l);
//     FTENSOR_INDEX(SPACE_DIM, m);

//     // double time = getTStime();
//     double time = 1.0;

//     auto nb_integration_pts = getGaussPts().size2();
//     // get element volume
//     const double vol = OpBase::getMeasure();
//     // get base function gradient on rows
//     auto t_row_base = data.getFTensor0N();
//     // get integration weights
//     auto t_w = OpBase::getFTensor0IntegrationWeight();
//     // get coordinate at integration points
//     auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
//     // get displacement
//     auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
//     // get lambda
//     auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

//     FTensor::Tensor1<double, 3> t_translation;
//     t_translation(0) = (*translationPtr)(0);
//     t_translation(1) = (*translationPtr)(1);
//     t_translation(2) = (*translationPtr)(2);

//     FTensor::Tensor1<double, 3> t_theta;
//     t_theta(0) = (*thetaPtr)(0);
//     t_theta(1) = (*thetaPtr)(1);
//     t_theta(2) = (*thetaPtr)(2);

//     // create a "dummy" base function
//     MatrixDouble m_kd(3, 3);
//     m_kd(0, 0) = 1.0;
//     m_kd(0, 1) = 0.0;
//     m_kd(0, 2) = 0.0;
//     m_kd(1, 0) = 0.0;
//     m_kd(1, 1) = 1.0;
//     m_kd(1, 2) = 0.0;
//     m_kd(2, 0) = 0.0;
//     m_kd(2, 1) = 0.0;
//     m_kd(2, 2) = 1.0;

//     auto t_dummy_base = getFTensor1FromMat<3>(m_kd);

//     for (int gg = 0; gg != nb_integration_pts; gg++) {
//       const auto alpha = t_w * vol;
//       ++t_w;

//       FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
//       t_tangent(i, k) = alpha * levi_civita(i,j,k) * (t_coords(j) -
//       tieCoord(j));

//       ++ t_lambda;
//       ++t_u;
//       ++t_coords;

//       int rr = 0;
//       for (; rr != 3; ++rr) {
//         auto t_mat = getFTensor1FromPtr<SPACE_DIM>(&OpBase::locMat(rr, 0));
//         auto t_col_base = data.getFTensor0N(gg, 0);
//         for (int cc = 0; cc != OpBase::nbCols/SPACE_DIM; cc++) {
//           t_mat(j) -= t_col_base * t_dummy_base(i) * t_tangent(i, j);
//           ++t_col_base;
//           ++t_mat;
//         }
//         ++t_dummy_base;
//       }
//     }
//     MoFEMFunctionReturn(0);
//   }

// private:
//   boost::shared_ptr<MatrixDouble> uPtr;
//   FTensor::Tensor1<double, 3> tieCoord;
//   FTensor::Tensor1<double, 3> tieDirection;
//   boost::shared_ptr<MatrixDouble> lambdaPtr;
//   boost::shared_ptr<VectorDouble> translationPtr;
//   boost::shared_ptr<VectorDouble> thetaPtr;
// };

// template <int Tensor_Dim>
struct OpCalculateNoFieldVectorValues
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateNoFieldVectorValues(const std::string field_name,
                                 boost::shared_ptr<VectorDouble> data_ptr,
                                 const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    VectorDouble &vec = *dataPtr;
    const size_t nb_gauss_pts = getGaussPts().size2();
    if (type == zeroType || vec.size() != nb_gauss_pts) {
      vec.resize(nb_gauss_pts, false);
      vec.clear();
    }

    const size_t nb_dofs = data.getFieldData().size();

    if (nb_dofs) {

      if (dataVec.use_count()) {
        dotVector.resize(nb_dofs, false);
        const double *array;
        CHKERR VecGetArrayRead(dataVec, &array);
        const auto &local_indices = data.getLocalIndices();
        for (int i = 0; i != local_indices.size(); ++i)
          if (local_indices[i] != -1)
            dotVector[i] = array[local_indices[i]];
          else
            dotVector[i] = 0;
        CHKERR VecRestoreArrayRead(dataVec, &array);
        data.getFieldData().swap(dotVector);
      }

      vec = data.getFieldData();

      // auto values_at_gauss_pts = getFTensor0FromVec(vec);

      // for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
      //   auto field_data = data.getFTensor1FieldData<SPACE_DIM>();
      //   size_t bb = 0;
      //   for (; bb != nb_dofs; ++bb) {
      //     values_at_gauss_pts += field_data;
      //     //++field_data;
      //     //++base_function;
      //   }
      //   // It is possible to have more base functions than dofs
      //   for (; bb < nb_base_functions; ++bb)
      //     ++base_function;
      //   ++values_at_gauss_pts;
      // }

      if (dataVec.use_count()) {
        data.getFieldData().swap(dotVector);
      }
    }

    MoFEMFunctionReturn(0);
  }

protected:
  boost::shared_ptr<ublas::vector<double, DoubleAllocator>> dataPtr;
  const EntityHandle zeroType;
  SmartPetscObj<Vec> dataVec;
  VectorDouble dotVector;
};

template <int DIM, AssemblyType AT>
struct OpCalculateRigidBodyTieReaction
    : public FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase {

  using OpBase = typename FormsIntegrators<BoundaryEleOp>::Assembly<AT>::OpBase;

  OpCalculateRigidBodyTieReaction(std::string lambda_name,
                                  boost::shared_ptr<MatrixDouble> lambda_ptr,
                                  SmartPetscObj<Vec> total_reaction,
                                  boost::shared_ptr<Range> tie_faces_ptr,
                                  FTensor::Tensor1<double, 3> tie_coord)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        lambdaPtr(lambda_ptr), totalReaction(total_reaction),
        tieFacesPtr(tie_faces_ptr), tieCoord(tie_coord) {
    std::fill(&(OpBase::doEntities[MBEDGE]), &(OpBase::doEntities[MBMAXTYPE]), false);
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    if (tieFacesPtr->find(OpBase::getNumeredEntFiniteElementPtr()->getEnt()) ==
        tieFacesPtr->end()) {
      MoFEMFunctionReturnHot(0);
    }
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);

    FTensor::Tensor1<double, 3> t_sum_reaction{0.0, 0.0, 0.0};
    FTensor::Tensor1<double, 3> t_sum_moment{0.0, 0.0, 0.0};

    auto nb_integration_pts = OpBase::getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function on rows
    auto t_row_base = data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get coordinate at integration points
    auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
    // get lambda
    auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      t_sum_reaction(i) += alpha * t_lambda(i);

      t_sum_moment(i) += alpha * levi_civita(i, j, k) *
                         (t_coords(k) - tieCoord(k)) * t_lambda(j);
      ++t_lambda;
      ++t_coords;
    }

    constexpr int ind_force[] = {0, 1, 2};
    constexpr int ind_moment[] = {3, 4, 5};

    // set reaction
    CHKERR VecSetValues(totalReaction, 3, ind_force, &t_sum_reaction(0),
                        ADD_VALUES);
    CHKERR VecSetValues(totalReaction, 3, ind_moment, &t_sum_moment(0),
                        ADD_VALUES);

    MoFEMFunctionReturn(0);
  }
  SmartPetscObj<Vec> totalReaction;
  boost::shared_ptr<MatrixDouble> lambdaPtr;
  boost::shared_ptr<Range> tieFacesPtr;
  FTensor::Tensor1<double, 3> tieCoord;
};

template <int DIM, AssemblyType AT>
MoFEMErrorCode OpFactoryCalculateTieConstraintForceTermRhs(
    MoFEM::Interface &m_field,
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string field_name,
    boost::shared_ptr<std::vector<RigidBodyTieConstraintData::TieBlock>>
        tie_blocks_ptr,
    Sev sev) {
  MoFEMFunctionBegin;
  auto lambda_ptr = boost::make_shared<MatrixDouble>();
  auto u_ptr = boost::make_shared<MatrixDouble>();
  pip.push_back(new OpCalculateVectorFieldValues<DIM>("U", u_ptr));
  pip.push_back(new OpCalculateVectorFieldValues<DIM>(field_name, lambda_ptr));

  for (auto &t : *tie_blocks_ptr) {
    pip.push_back(new OpTieTermConstrainRigidBodyRhs_du<DIM, AT>(
        "U", lambda_ptr, boost::make_shared<Range>(t.tieFaces)));
  }
  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType AT>
MoFEMErrorCode OpFactoryCalculateRigidBodyConstraintRhs(
    MoFEM::Interface &m_field,
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string field_name,
    boost::shared_ptr<std::vector<RigidBodyTieConstraintData::TieBlock>>
        tie_blocks_ptr,
    Sev sev) {
  MoFEMFunctionBegin;

  auto simple = m_field.getInterface<Simple>();

  auto lambda_ptr = boost::make_shared<MatrixDouble>();
  auto translation_ptr = boost::make_shared<VectorDouble>();
  auto u_ptr = boost::make_shared<MatrixDouble>();
  auto int_translation_ptr = boost::make_shared<VectorDouble>();
  auto int_rotation_ptr = boost::make_shared<MatrixDouble>();

  Range rigid_body_ents;
  auto v_rigid_body_ents = simple->getMeshsetFiniteElementEntities();

  // loop over all rigid body entities
  for (auto ents : v_rigid_body_ents) {
    rigid_body_ents.merge(ents);
  }

  auto rigid_body_ents_ptr = boost::make_shared<Range>(rigid_body_ents);

  for (auto &t : *tie_blocks_ptr) {
    auto op_loop_side = new OpLoopSide<BoundaryEle>(m_field, "bFE", DIM - 1,
                                                    rigid_body_ents_ptr);
    CHKERR AddHOOps<DIM - 1, DIM, DIM>::add(op_loop_side->getOpPtrVector(), {});
    op_loop_side->getOpPtrVector().push_back(
        new OpCalculateVectorFieldValues<DIM>(field_name, lambda_ptr));
    op_loop_side->getOpPtrVector().push_back(
        new OpTieTermConstrainRigidBodyGlobalTranslationIntegralRhs<DIM, AT>(
            field_name, boost::make_shared<Range>(t.tieFaces), lambda_ptr,
            int_translation_ptr));
    op_loop_side->getOpPtrVector().push_back(
        new OpTieTermConstrainRigidBodyGlobalRotationIntegralRhs<DIM, AT>(
            field_name, t.tieCoord, lambda_ptr, int_rotation_ptr));
    pip.push_back(op_loop_side);

    pip.push_back(new OpTieTermConstrainRigidBodyGlobalTranslationRhs<DIM>(
        "RIGID_BODY_LAMBDA", int_translation_ptr, rigid_body_ents_ptr));
    pip.push_back(new OpTieTermConstrainRigidBodyGlobalRotationRhs<DIM>(
        "RIGID_BODY_THETA", int_rotation_ptr, rigid_body_ents_ptr));
  }

  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType AT>
MoFEMErrorCode OpFactoryCalculateRigidBodyConstraintLhs(
    MoFEM::Interface &m_field,
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string field_name,
    boost::shared_ptr<std::vector<RigidBodyTieConstraintData::TieBlock>>
        tie_blocks_ptr,
    Sev sev) {
  MoFEMFunctionBegin;

  auto simple = m_field.getInterface<Simple>();
  auto lambda_ptr = boost::make_shared<MatrixDouble>();
  auto translation_ptr = boost::make_shared<VectorDouble>();
  auto u_ptr = boost::make_shared<MatrixDouble>();

  Range rigid_body_ents;
  auto v_rigid_body_ents = simple->getMeshsetFiniteElementEntities();

  // loop over all rigid body entities
  for (auto ents : v_rigid_body_ents) {
    rigid_body_ents.merge(ents);
  }

  auto rigid_body_ents_ptr = boost::make_shared<Range>(rigid_body_ents);

  for (auto &t : *tie_blocks_ptr) {
    auto op_loop_side = new OpLoopSide<BoundaryEle>(m_field, "bFE", DIM - 1,
                                                    rigid_body_ents_ptr);
    CHKERR AddHOOps<DIM - 1, DIM, DIM>::add(op_loop_side->getOpPtrVector(), {});
    op_loop_side->getOpPtrVector().push_back(
        new OpCalculateVectorFieldValues<DIM>(field_name, lambda_ptr));
    op_loop_side->getOpPtrVector().push_back(
        new OpTieTermConstrainRigidBodyLhs_dTranslation<DIM, AT>(
            field_name, "RIGID_BODY_LAMBDA",
            boost::make_shared<Range>(t.tieFaces)));
    op_loop_side->getOpPtrVector().push_back(
        new OpTieTermConstrainRigidBodyLhs_dRotation<DIM, AT>(
            field_name, "RIGID_BODY_THETA", t.tieCoord,
            boost::make_shared<Range>(t.tieFaces)));

    pip.push_back(op_loop_side);
  }

  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType AT>
MoFEMErrorCode OpFactoryCalculateTieConstraintRhs(
    MoFEM::Interface &m_field,
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string field_name,
    boost::shared_ptr<std::vector<RigidBodyTieConstraintData::TieBlock>>
        tie_blocks_ptr,
    Sev sev) {
  MoFEMFunctionBegin;
  auto u_ptr = boost::make_shared<MatrixDouble>();
  auto lambda_ptr = boost::make_shared<MatrixDouble>();
  auto translation_ptr = boost::make_shared<VectorDouble>();
  auto theta_ptr = boost::make_shared<VectorDouble>();

  pip.push_back(new OpCalculateVectorFieldValues<DIM>("U", u_ptr));
  pip.push_back(new OpCalculateVectorFieldValues<DIM>(field_name, lambda_ptr));

  pip.push_back(
      new OpCalculateNoFieldVectorValues("RIGID_BODY_LAMBDA", translation_ptr));
  pip.push_back(
      new OpCalculateNoFieldVectorValues("RIGID_BODY_THETA", theta_ptr));

  for (auto &t : *tie_blocks_ptr) {
    pip.push_back(new OpTieTermConstrainRigidBodyRhs<DIM, AT>(
        field_name, u_ptr, t.tieCoord, boost::make_shared<Range>(t.tieFaces),
        translation_ptr, theta_ptr));
  }
  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType AT>
MoFEMErrorCode OpFactoryCalculateTieConstraintLhs(
    MoFEM::Interface &m_field,
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string field_name,
    boost::shared_ptr<std::vector<RigidBodyTieConstraintData::TieBlock>>
        tie_blocks_ptr,
    Sev sev) {
  MoFEMFunctionBegin;

  auto u_ptr = boost::make_shared<MatrixDouble>();
  auto lambda_ptr = boost::make_shared<MatrixDouble>();
  auto translation_ptr = boost::make_shared<VectorDouble>();
  auto theta_ptr = boost::make_shared<VectorDouble>();

  pip.push_back(new OpCalculateVectorFieldValues<DIM>("U", u_ptr));
  pip.push_back(new OpCalculateVectorFieldValues<DIM>(field_name, lambda_ptr));
  pip.push_back(
      new OpCalculateNoFieldVectorValues("RIGID_BODY_LAMBDA", translation_ptr));
  pip.push_back(
      new OpCalculateNoFieldVectorValues("RIGID_BODY_THETA", theta_ptr));

  for (auto &t : *tie_blocks_ptr) {
    pip.push_back(new OpTieTermConstrainRigidBodyLhs_dU<DIM, AT>(
        field_name, "U", u_ptr, t.tieCoord,
        boost::make_shared<Range>(t.tieFaces), translation_ptr, theta_ptr));
  }
  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType AT>
MoFEMErrorCode OpFactoryPipelineTieConstraint(
    MoFEM::Interface &m_field, std::string field_name,
    boost::shared_ptr<std::vector<RigidBodyTieConstraintData::TieBlock>>
        tie_block_ptr,
    Sev sev) {
  MoFEMFunctionBegin;
  auto *pipeline_mng = m_field.getInterface<PipelineManager>();

  CHKERR OpFactoryCalculateTieConstraintLhs<SPACE_DIM, AT>(
      m_field, pipeline_mng->getOpBoundaryLhsPipeline(), field_name,
      tie_block_ptr, Sev::inform);

  CHKERR OpFactoryCalculateTieConstraintRhs<SPACE_DIM, AT>(
      m_field, pipeline_mng->getOpBoundaryRhsPipeline(), field_name,
      tie_block_ptr, Sev::inform);

  CHKERR OpFactoryCalculateTieConstraintForceTermRhs<SPACE_DIM, AT>(
      m_field, pipeline_mng->getOpBoundaryRhsPipeline(), field_name,
      tie_block_ptr, Sev::inform);

  CHKERR OpFactoryCalculateRigidBodyConstraintRhs<SPACE_DIM, AT>(
      m_field, pipeline_mng->getOpMeshsetRhsPipeline(), field_name,
      tie_block_ptr, Sev::inform);

  CHKERR OpFactoryCalculateRigidBodyConstraintLhs<SPACE_DIM, AT>(
      m_field, pipeline_mng->getOpMeshsetLhsPipeline(), field_name,
      tie_block_ptr, Sev::inform);

  MoFEMFunctionReturn(0);
}

struct EssentialPostProcRigidBodyTieRhs {
  EssentialPostProcRigidBodyTieRhs(
      MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr,
      boost::shared_ptr<std::vector<RigidBodyTieConstraintData::TieBlock>>
          tie_blocks_ptr);

  virtual ~EssentialPostProcRigidBodyTieRhs() = default;

  MoFEMErrorCode operator()();

protected:
  MoFEM::Interface &mField;
  boost::weak_ptr<FEMethod> fePtr;
  boost::shared_ptr<std::vector<RigidBodyTieConstraintData::TieBlock>>
      tieBlocksPtr;
};

EssentialPostProcRigidBodyTieRhs::EssentialPostProcRigidBodyTieRhs(
    MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr,
    boost::shared_ptr<std::vector<RigidBodyTieConstraintData::TieBlock>>
        tie_blocks_ptr)
    : mField(m_field), fePtr(fe_ptr), tieBlocksPtr(tie_blocks_ptr) {}

MoFEMErrorCode EssentialPostProcRigidBodyTieRhs::operator()() {
  MoFEMFunctionBegin;
  auto is_mng = mField.getInterface<ISManager>();
  auto simple = mField.getInterface<Simple>();

  auto set_rhs = [&](auto index, auto field_name) {
    MoFEMFunctionBegin;
    SmartPetscObj<IS> is;
    CHKERR is_mng->isCreateProblemFieldAndRankLocal(
        simple->getProblemName(), ROW, field_name, 0, MAX_DOFS_ON_ENTITY, is,
        nullptr);

    if (*(fePtr.lock()->vecAssembleSwitch)) {
      CHKERR VecAssemblyBegin(fePtr.lock()->f);
      CHKERR VecAssemblyEnd(fePtr.lock()->f);
      CHKERR VecGhostUpdateBegin(fePtr.lock()->f, ADD_VALUES, SCATTER_REVERSE);
      CHKERR VecGhostUpdateEnd(fePtr.lock()->f, ADD_VALUES, SCATTER_REVERSE);
      *(fePtr.lock()->vecAssembleSwitch) = false;
    }

    CHKERR VecGhostUpdateBegin(fePtr.lock()->x, INSERT_VALUES, SCATTER_FORWARD);
    CHKERR VecGhostUpdateEnd(fePtr.lock()->x, INSERT_VALUES, SCATTER_FORWARD);

    const int *i_ptr;
    CHKERR ISGetIndices(is, &i_ptr);
    int is_size;
    CHKERR ISGetLocalSize(is, &is_size);
    if (is_size != SPACE_DIM)
      SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "Wrong size of rigid dofs");
    double *f_ptr;
    double *x_ptr;
    auto fePtrLocked = fePtr.lock();
    CHKERR VecGetArray(fePtrLocked->f, &f_ptr);
    CHKERR VecGetArray(fePtrLocked->x, &x_ptr);

    if (mField.get_comm_rank() == 0) {
      if (field_name == "RIGID_BODY_LAMBDA") {
        f_ptr[i_ptr[index]] =
            x_ptr[i_ptr[index]] -
            (tieBlocksPtr->at(0).tieDirection(index) * fePtrLocked->ts_t);
      } else if (field_name == "RIGID_BODY_THETA") {
        f_ptr[i_ptr[index]] =
            x_ptr[i_ptr[index]] -
            (tieBlocksPtr->at(0).tieRotation(index) * fePtrLocked->ts_t);
      } else {
        SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "Wrong field name");
      }
    }
    CHKERR VecRestoreArray(fePtrLocked->f, &f_ptr);
    CHKERR VecRestoreArray(fePtrLocked->x, &x_ptr);
    CHKERR ISRestoreIndices(is, &i_ptr);
    MoFEMFunctionReturn(0);
  };

  // set block bcs
  for (int i = 0; i < 3; ++i) {
    if (tieBlocksPtr->at(0).tieTranslationFlag(i)) {
      set_rhs(i, "RIGID_BODY_LAMBDA");
    }
    if (tieBlocksPtr->at(0).tieRotationFlag(i)) {
      set_rhs(i, "RIGID_BODY_THETA");
    }
  }
  MoFEMFunctionReturn(0);
}

struct EssentialPostProcRigidBodyTieLhs {
  EssentialPostProcRigidBodyTieLhs(
      MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr,
      boost::shared_ptr<std::vector<RigidBodyTieConstraintData::TieBlock>>
          tie_blocks_ptr);

  virtual ~EssentialPostProcRigidBodyTieLhs() = default;

  MoFEMErrorCode operator()();

protected:
  MoFEM::Interface &mField;
  boost::weak_ptr<FEMethod> fePtr;
  boost::shared_ptr<std::vector<RigidBodyTieConstraintData::TieBlock>>
      tieBlocksPtr;
};

EssentialPostProcRigidBodyTieLhs::EssentialPostProcRigidBodyTieLhs(
    MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr,
    boost::shared_ptr<std::vector<RigidBodyTieConstraintData::TieBlock>>
        tie_blocks_ptr)
    : mField(m_field), fePtr(fe_ptr), tieBlocksPtr(tie_blocks_ptr) {}

MoFEMErrorCode EssentialPostProcRigidBodyTieLhs::operator()() {
  MoFEMFunctionBegin;
  auto is_mng = mField.getInterface<ISManager>();
  auto simple = mField.getInterface<Simple>();

  auto set_lhs = [&](auto index, auto field_name) {
    MoFEMFunctionBegin;
    auto fePtrLocked = fePtr.lock();
    SmartPetscObj<IS> is;
    CHKERR is_mng->isCreateProblemFieldAndRank(
        simple->getProblemName(), ROW, field_name, 0, SPACE_DIM, is, nullptr);
    //*(fePtrLocked->matAssembleSwitch) = PETSC_TRUE;

    if (*(fePtr.lock()->matAssembleSwitch)) {
      CHKERR MatAssemblyBegin(fePtrLocked->B, MAT_FINAL_ASSEMBLY);
      CHKERR MatAssemblyEnd(fePtrLocked->B, MAT_FINAL_ASSEMBLY);
      *(fePtr.lock()->matAssembleSwitch) = false;
    }
    // double norm;
    // CHKERR MatNorm(fePtrLocked->B, NORM_FROBENIUS, &norm);
    // MOFEM_LOG("SYNC", Sev::inform) << "norm before = " << norm;
    CHKERR MatSetOption(fePtrLocked->B, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);

    if (mField.get_comm_rank() == 0) {
      int is_size;
      CHKERR ISGetSize(is, &is_size);
      if (is_size != SPACE_DIM)
        SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "Wrong size of rigid dofs");
      const int *i_ptr;
      CHKERR ISGetIndices(is, &i_ptr);
      const int idx = index;

      CHKERR MatZeroRows(fePtrLocked->B, 1, &i_ptr[idx], 1, PETSC_NULL,
                         PETSC_NULL);

      CHKERR ISRestoreIndices(is, &i_ptr);
    } else {
      CHKERR MatZeroRows(fePtrLocked->B, 0, PETSC_NULL, 1, PETSC_NULL,
                         PETSC_NULL);
    }
    // double norm2;
    // CHKERR MatNorm(fePtrLocked->B, NORM_FROBENIUS, &norm2);
    // MOFEM_LOG("SYNC", Sev::inform) << "norm after = " << norm2;
    MoFEMFunctionReturn(0);
  };

  // set block bcs
  for (int i = 0; i < 3; ++i) {
    if (tieBlocksPtr->at(0).tieTranslationFlag(i)) {
      set_lhs(i, "RIGID_BODY_LAMBDA");
    }
    if (tieBlocksPtr->at(0).tieRotationFlag(i)) {
      set_lhs(i, "RIGID_BODY_THETA");
    }
  }

  MoFEMFunctionReturn(0);
}