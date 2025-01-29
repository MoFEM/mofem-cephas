/**
 * @file TieConstraint.hpp
 * @brief Tie Constraint Implementation
 * @date 2024-11-14
 *
 * @copyright Copyright (c) 2024
 *
 */
template <int DIM>
struct OpTieTermConstrainDistanceRhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainDistanceRhs(std::string lambda_name,
                                boost::shared_ptr<MatrixDouble> u_ptr,
                                FTensor::Tensor1<double, 3> tie_coord,
                                FTensor::Tensor1<double, 3> tie_direction,
                                boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor1<double, SPACE_DIM> t_delta_current;
      t_delta_current(i) =
          (t_coords(i) + t_u(i)) - (tieCoord(i) + (tieDirection(i) * time));

      FTensor::Tensor1<double, SPACE_DIM> t_delta_initial;
      t_delta_initial(i) = t_coords(i) - tieCoord(i);
      ++t_u;
      ++t_coords;
      auto g = alpha * (t_delta_current.l2() - t_delta_initial.l2());

      int rr = 0;
      for (; rr != OpBase::nbRows; ++rr) {
        OpBase::locF[rr] += t_row_base * g;
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
  FTensor::Tensor1<double, 3> tieDirection;
};

template <int DIM>
struct OpTieTermConstrainDistanceLhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {
  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainDistanceLhs(std::string lambda_name,
                                std::string col_field_name,
                                boost::shared_ptr<MatrixDouble> u_ptr,
                                FTensor::Tensor1<double, 3> tie_coord,
                                FTensor::Tensor1<double, 3> tie_direction,
                                boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, col_field_name, OpBase::OPROWCOL, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
    this->assembleTranspose = true;
    this->sYmm = false;
  }

  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;
      FTensor::Tensor1<double, SPACE_DIM> t_delta_current;
      t_delta_current(i) =
          (t_coords(i) + t_u(i)) - (tieCoord(i) + (tieDirection(i) * time));
      ++t_u;
      ++t_coords;
      FTensor::Tensor1<double, SPACE_DIM> t_tangent;
      t_tangent(i) = alpha * (t_delta_current(i) / t_delta_current.l2());

      int rr = 0;
      for (; rr != OpBase::nbRows; ++rr) {
        auto t_col_base = col_data.getFTensor0N(gg, 0);
        auto t_mat = getFTensor1FromPtr<SPACE_DIM>(&OpBase::locMat(rr, 0));
        for (int cc = 0; cc != OpBase::nbCols / SPACE_DIM; cc++) {
          t_mat(i) += (t_row_base * t_col_base) * t_tangent(i);
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
  FTensor::Tensor1<double, 3> tieDirection;
  double tsTime;
};

template <int DIM>
struct OpTieTermConstrainDistanceRhs_du
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainDistanceRhs_du(std::string field_name,
                                   std::string lambda_name,
                                   boost::shared_ptr<MatrixDouble> u_ptr,
                                   boost::shared_ptr<VectorDouble> lambda_ptr,
                                   FTensor::Tensor1<double, 3> tie_coord,
                                   FTensor::Tensor1<double, 3> tie_direction,
                                   boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(field_name, field_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), lambdaPtr(lambda_ptr), tieCoord(tie_coord),
        tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function on rows
    auto t_row_base = data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get coordinate at integration points
    auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
    // get displacement
    auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
    // get lambda
    auto t_lambda = getFTensor0FromVec(*lambdaPtr);

    auto &nf = OpBase::locF;

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;
      FTensor::Tensor1<double, SPACE_DIM> t_du;
      t_du(i) = t_u(i) + (t_coords(i) - tieCoord(i) - (tieDirection(i) * time));
      ++t_u;
      ++t_coords;
      
      auto t_nf = getFTensor1FromPtr<DIM>(&nf[0]);
      FTensor::Tensor1<double, SPACE_DIM> g;
      g(i) = alpha * t_lambda * (t_du(i) / t_du.l2());
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
  boost::shared_ptr<MatrixDouble> uPtr;
  boost::shared_ptr<VectorDouble> lambdaPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  FTensor::Tensor1<double, 3> tieDirection;
};

// operator to calculate total reaction force
template <int DIM>
struct OpCalculateTieDistanceReactionForce
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpCalculateTieDistanceReactionForce(
      std::string lambda_name, boost::shared_ptr<VectorDouble> lambda_ptr,
      boost::shared_ptr<SmartPetscObj<Vec>> total_reaction_ptr,
      boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        lambdaPtr(lambda_ptr), totalReactionPtr(total_reaction_ptr),
        tieFacesPtr(tie_faces_ptr) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    // doEntities[MBEDGE] = true;
    // doEntities[MBTRI] = true;
    // doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    if (tieFacesPtr->find(getNumeredEntFiniteElementPtr()->getEnt()) ==
        tieFacesPtr->end()) {
      MoFEMFunctionReturnHot(0);
    }
    // std::cout<< "Entity id = "<< getNumeredEntFiniteElementPtr()->getEnt() <<
    // std::endl;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);

    double sum_reaction = 0.0;
    FTensor::Tensor1<double, 3> t_sum_reaction{0.0, 0.0, 0.0};

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
    auto t_lambda = getFTensor0FromVec(*lambdaPtr);
    // get stress
    // auto t_stress = getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(*stressPtr);
    // get normal
    // auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      auto g = alpha * t_lambda;
      ++t_lambda;

      // t_sum_reaction(i) += alpha * t_stress(i, j) * t_normal(i);

      // int rr = 0;
      // for (; rr != OpBase::nbRows / DIM; ++rr) {
      sum_reaction += g;
      //++t_row_base;
      //}

      // for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      //  ++t_row_base;
    }

    constexpr int ind[] = {1, 2, 3};
    constexpr int ind0[] = {0};
    // std::cout << "lambdaPtr: " << *lambdaPtr << std::endl;
    // std::cout << "sum_reaction: " << sum_reaction << std::endl;
    // std::cout << "t_sum_reaction: " << t_sum_reaction(0) << " " <<
    // t_sum_reaction(1) << " " << t_sum_reaction(2) << std::endl; std::cout <<
    // "total_reaction_ptr: " << *totalReactionPtr << std::endl; set reaction
    // from LM
    CHKERR VecSetValues(*totalReactionPtr, 1, ind0, &sum_reaction, ADD_VALUES);
    // set reaction from stress
    CHKERR VecSetValues(*totalReactionPtr, 3, ind, &t_sum_reaction(0),
                        ADD_VALUES);

    MoFEMFunctionReturn(0);
  }
  boost::shared_ptr<SmartPetscObj<Vec>> totalReactionPtr;
  boost::shared_ptr<VectorDouble> lambdaPtr;
  boost::shared_ptr<MatrixDouble> stressPtr;
  boost::shared_ptr<Range> tieFacesPtr;
};

template <int DIM>
struct OpTieTermConstrainRotationRhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRotationRhs(std::string lambda_name,
                                boost::shared_ptr<MatrixDouble> u_ptr,
                                FTensor::Tensor1<double, 3> tie_coord,
                                FTensor::Tensor1<double, 3> tie_direction,
                                boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor1<double, SPACE_DIM> t_rotation;
      t_rotation(i) = (t_coords(i) - tieCoord(i));

      FTensor::Tensor1<double, SPACE_DIM> t_delta_disp;
      t_delta_disp(i) = t_u(i)  - (tieDirection(i) * time);

      FTensor::Tensor1<double, SPACE_DIM> g;
      g(i) = alpha * (levi_civita(i, j, k) * t_rotation(j) * t_delta_disp(k));

      ++t_coords;
      ++t_u;

      auto t_nf = getFTensor1FromArray<DIM, DIM>(OpBase::locF);

      int rr = 0;
      for (; rr != OpBase::nbRows/DIM; ++rr) {
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
  FTensor::Tensor1<double, 3> tieDirection;
};

template <int DIM>
struct OpTieTermConstrainRotationLhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {
  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRotationLhs(std::string lambda_name, std::string col_field_name,
                         boost::shared_ptr<MatrixDouble> u_ptr,
                         FTensor::Tensor1<double, 3> tie_coord,
                         FTensor::Tensor1<double, 3> tie_direction,
                         boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, col_field_name, OpBase::OPROWCOL, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
    this->assembleTranspose = true;
    this->sYmm = false;
  }

  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor1<double, SPACE_DIM> t_rotation;
      t_rotation(i) = (t_coords(i) - tieCoord(i));
      FTensor::Tensor1<double, SPACE_DIM> t_delta_disp;
      t_delta_disp(i) = t_u(i) - (tieDirection(i) * time);

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
      t_tangent(i, l) = alpha * (levi_civita(i, j, k) * t_rotation(j) * t_kd(k, l));
      ++t_u;
      ++t_coords;

      int rr = 0;
      for (; rr != OpBase::nbRows/SPACE_DIM; ++rr) {
        auto t_col_base = col_data.getFTensor0N(gg, 0);
        auto t_mat = getFTensor2FromArray<DIM, DIM, DIM>(OpBase::locMat, DIM * rr);
        for (int cc = 0; cc != OpBase::nbCols/SPACE_DIM; cc++) {
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
  FTensor::Tensor1<double, 3> tieDirection;
};

template <int DIM> 
struct OpTieTermConstrainRotationRhs_du
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRotationRhs_du(std::string field_name, std::string lambda_name,
                            boost::shared_ptr<MatrixDouble> u_ptr,
                            boost::shared_ptr<MatrixDouble> lambda_ptr,
                            FTensor::Tensor1<double, 3> tie_coord,
                            FTensor::Tensor1<double, 3> tie_direction,
                            boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(field_name, field_name, OpBase::OPROW, tie_faces_ptr), uPtr(u_ptr),
        lambdaPtr(lambda_ptr), tieCoord(tie_coord),
        tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function on rows
    auto t_row_base = data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get coordinate at integration points
    auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
    // get displacement
    auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
    // get lambda
    auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

    auto &nf = OpBase::locF;

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor1<double, SPACE_DIM> t_rotation;
      t_rotation(i) = (t_coords(i) - tieCoord(i));
      FTensor::Tensor1<double, SPACE_DIM> t_delta_disp;
      t_delta_disp(i) = t_u(i)  - (tieDirection(i) * time);

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_du;
      t_du(i,l) = (levi_civita(i,j,k) * t_rotation(j) * t_kd(k, l));
      ++t_u;
      ++t_coords;
      
      auto t_nf = getFTensor1FromPtr<DIM>(&nf[0]);
      FTensor::Tensor1<double, SPACE_DIM> g;
      g(l) = alpha * t_lambda(i) * (t_du(i, l));
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
  boost::shared_ptr<MatrixDouble> uPtr;
  boost::shared_ptr<MatrixDouble> lambdaPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  FTensor::Tensor1<double, 3> tieDirection;
};

// operator to calculate total reaction force
template <int DIM>
struct OpCalculateTieRotationReactionForce
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpCalculateTieRotationReactionForce(
      std::string lambda_name, boost::shared_ptr<VectorDouble> lambda_ptr,
      boost::shared_ptr<SmartPetscObj<Vec>> total_reaction_ptr,
      boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        lambdaPtr(lambda_ptr), totalReactionPtr(total_reaction_ptr),
        tieFacesPtr(tie_faces_ptr) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    // doEntities[MBEDGE] = true;
    // doEntities[MBTRI] = true;
    // doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    if (tieFacesPtr->find(getNumeredEntFiniteElementPtr()->getEnt()) ==
        tieFacesPtr->end()) {
      MoFEMFunctionReturnHot(0);
    }
    // std::cout<< "Entity id = "<< getNumeredEntFiniteElementPtr()->getEnt() <<
    // std::endl;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);

    double sum_reaction = 0.0;
    FTensor::Tensor1<double, 3> t_sum_reaction{0.0, 0.0, 0.0};

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
    // get stress
    // auto t_stress = getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(*stressPtr);
    // get normal
    auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      // auto g = alpha * t_lambda;
      // ++t_lambda;

      // t_sum_reaction(i) += g * (t_normal(i)/t_normal.l2());
      // ++t_normal;

      // int rr = 0;
      // for (; rr != OpBase::nbRows / DIM; ++rr) {
      // sum_reaction += g;
      //++t_row_base;
      //}

      // for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      //  ++t_row_base;
    }

    constexpr int ind[] = {1, 2, 3};
    constexpr int ind0[] = {0};
    // std::cout << "lambdaPtr: " << *lambdaPtr << std::endl;
    // std::cout << "sum_reaction: " << sum_reaction << std::endl;
    // std::cout << "t_sum_reaction: " << t_sum_reaction(0) << " " <<
    // t_sum_reaction(1) << " " << t_sum_reaction(2) << std::endl; std::cout <<
    // "total_reaction_ptr: " << *totalReactionPtr << std::endl; set reaction
    // from LM
    CHKERR VecSetValues(*totalReactionPtr, 1, ind0, &sum_reaction, ADD_VALUES);
    // set reaction from stress
    CHKERR VecSetValues(*totalReactionPtr, 3, ind, &t_sum_reaction(0),
                        ADD_VALUES);

    MoFEMFunctionReturn(0);
  }
  boost::shared_ptr<SmartPetscObj<Vec>> totalReactionPtr;
  boost::shared_ptr<MatrixDouble> lambdaPtr;
  boost::shared_ptr<MatrixDouble> stressPtr;
  boost::shared_ptr<Range> tieFacesPtr;
};



template <int DIM>
struct OpTieTermConstrainRotationNormalRhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRotationNormalRhs(std::string lambda_name,
                                boost::shared_ptr<MatrixDouble> u_ptr,
                                FTensor::Tensor1<double, 3> tie_coord,
                                FTensor::Tensor1<double, 3> tie_direction,
                                boost::shared_ptr<Range> tie_faces_ptr,
                                std::string rotation_plane)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction), rotationPlane(rotation_plane) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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
    //get normal
    //auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();
    FTensor::Tensor1<double, 3> t_normal;
    if (rotationPlane == "x")
    {
      t_normal(0) = 1.0;
      t_normal(1) = 0.0;
      t_normal(2) = 0.0;
    }
    else if (rotationPlane == "y")
    {
      t_normal(0) = 0.0;
      t_normal(1) = 1.0;
      t_normal(2) = 0.0;
    }
    else if (rotationPlane == "z")
    {
      t_normal(0) = 0.0;
      t_normal(1) = 0.0;
      t_normal(2) = 1.0;
    }

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor1<double, SPACE_DIM> t_rotation;
      t_rotation(i) = (t_coords(i) - tieCoord(i));

      FTensor::Tensor1<double, SPACE_DIM> t_delta_disp;
      t_delta_disp(i) = t_u(i)  - (tieDirection(i) * time);

      auto g = alpha * (levi_civita(i, j, k) * t_rotation(j) * t_delta_disp(k)) * t_normal(i);

      ++t_coords;
      ++t_u;
      //++t_normal;

      int rr = 0;
      for (; rr != OpBase::nbRows; ++rr) {
        OpBase::locF[rr] += t_row_base * g;
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
  FTensor::Tensor1<double, 3> tieDirection;
  std::string rotationPlane;
};

template <int DIM>
struct OpTieTermConstrainRotationNormalLhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {
  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRotationNormalLhs(std::string lambda_name, std::string col_field_name,
                         boost::shared_ptr<MatrixDouble> u_ptr,
                         FTensor::Tensor1<double, 3> tie_coord,
                         FTensor::Tensor1<double, 3> tie_direction,
                         boost::shared_ptr<Range> tie_faces_ptr,
                         std::string rotation_plane)
      : OpBase(lambda_name, col_field_name, OpBase::OPROWCOL, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction), rotationPlane(rotation_plane) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
    this->assembleTranspose = true;
    this->sYmm = false;
  }

  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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
    // get normal
    //auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();
    FTensor::Tensor1<double, 3> t_normal;
    if (rotationPlane == "x")
    {
      t_normal(0) = 1.0;
      t_normal(1) = 0.0;
      t_normal(2) = 0.0;
    }
    else if (rotationPlane == "y")
    {
      t_normal(0) = 0.0;
      t_normal(1) = 1.0;
      t_normal(2) = 0.0;
    }
    else if (rotationPlane == "z")
    {
      t_normal(0) = 0.0;
      t_normal(1) = 0.0;
      t_normal(2) = 1.0;
    }

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor1<double, SPACE_DIM> t_rotation;
      t_rotation(i) = (t_coords(i) - tieCoord(i));
      FTensor::Tensor1<double, SPACE_DIM> t_delta_disp;
      t_delta_disp(i) = t_u(i) - (tieDirection(i) * time);

      FTensor::Tensor1<double, SPACE_DIM> t_tangent;
      t_tangent(k) = alpha * (levi_civita(i, j, k) * t_rotation(j) * t_normal(i));

      ++t_u;
      ++t_coords;
      //++t_normal;

      int rr = 0;
      for (; rr != OpBase::nbRows; ++rr) {
        auto t_col_base = col_data.getFTensor0N(gg, 0);
        auto t_mat = getFTensor1FromPtr<SPACE_DIM>(&OpBase::locMat(rr, 0));
        for (int cc = 0; cc != OpBase::nbCols / SPACE_DIM; cc++) {
          t_mat(i) += (t_row_base * t_col_base) * t_tangent(i);
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
  FTensor::Tensor1<double, 3> tieDirection;
  std::string rotationPlane;
};

template <int DIM> 
struct OpTieTermConstrainRotationNormalRhs_du
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRotationNormalRhs_du(std::string field_name, std::string lambda_name,
                            boost::shared_ptr<MatrixDouble> u_ptr,
                            boost::shared_ptr<VectorDouble> lambda_ptr,
                            FTensor::Tensor1<double, 3> tie_coord,
                            FTensor::Tensor1<double, 3> tie_direction,
                            boost::shared_ptr<Range> tie_faces_ptr,
                            std::string rotation_plane)
      : OpBase(field_name, field_name, OpBase::OPROW, tie_faces_ptr), uPtr(u_ptr),
        lambdaPtr(lambda_ptr), tieCoord(tie_coord),
        tieDirection(tie_direction), rotationPlane(rotation_plane) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function on rows
    auto t_row_base = data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get coordinate at integration points
    auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
    // get displacement
    auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
    // get lambda
    auto t_lambda = getFTensor0FromVec(*lambdaPtr);
    // get normal
    //auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();
    // rotation plane to constrain

    FTensor::Tensor1<double, 3> t_normal;
    if (rotationPlane == "x")
    {
      t_normal(0) = 1.0;
      t_normal(1) = 0.0;
      t_normal(2) = 0.0;
    }
    else if (rotationPlane == "y")
    {
      t_normal(0) = 0.0;
      t_normal(1) = 1.0;
      t_normal(2) = 0.0;
    }
    else if (rotationPlane == "z")
    {
      t_normal(0) = 0.0;
      t_normal(1) = 0.0;
      t_normal(2) = 1.0;
    }

    auto &nf = OpBase::locF;


    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor1<double, SPACE_DIM> t_rotation;
      t_rotation(i) = (t_coords(i) - tieCoord(i));
      FTensor::Tensor1<double, SPACE_DIM> t_delta_disp;
      t_delta_disp(i) = t_u(i)  - (tieDirection(i) * time);

      FTensor::Tensor1<double, SPACE_DIM> t_du;
      //t_du(i) = (levi_civita(i,j,k) * t_rotation(j) * t_kd(k, l)) * t_normal(l);
      t_du(k) = levi_civita(i,j,k) * t_rotation(j) * t_normal(i);
      ++t_u;
      ++t_coords;
      //++t_normal;
      
      auto t_nf = getFTensor1FromPtr<DIM>(&nf[0]);
      FTensor::Tensor1<double, SPACE_DIM> g;
      g(i) = alpha * t_lambda * (t_du(i));
      //std::cout << "t_lambda: " << t_lambda << std::endl;
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
  boost::shared_ptr<MatrixDouble> uPtr;
  boost::shared_ptr<VectorDouble> lambdaPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  FTensor::Tensor1<double, 3> tieDirection;
  std::string rotationPlane;
};

// operator to calculate total reaction force
template <int DIM>
struct OpCalculateTieRotationNormalReactionForce
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpCalculateTieRotationNormalReactionForce(
      std::string lambda_name, boost::shared_ptr<VectorDouble> lambda_ptr,
      boost::shared_ptr<SmartPetscObj<Vec>> total_reaction_ptr,
      boost::shared_ptr<Range> tie_faces_ptr,
      FTensor::Tensor1<double, 3> tie_coord,
      std::string rotation_plane)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        lambdaPtr(lambda_ptr), totalReactionPtr(total_reaction_ptr),
        tieFacesPtr(tie_faces_ptr), tieCoord(tie_coord), rotationPlane(rotation_plane) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    // doEntities[MBEDGE] = true;
    // doEntities[MBTRI] = true;
    // doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    if (tieFacesPtr->find(getNumeredEntFiniteElementPtr()->getEnt()) ==
        tieFacesPtr->end()) {
      MoFEMFunctionReturnHot(0);
    }
    // std::cout<< "Entity id = "<< getNumeredEntFiniteElementPtr()->getEnt() <<
    // std::endl;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);

    double sum_reaction = 0.0;
    FTensor::Tensor1<double, 3> t_sum_reaction{0.0, 0.0, 0.0};

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
    auto t_lambda = getFTensor0FromVec(*lambdaPtr);
    // get stress
    // auto t_stress = getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(*stressPtr);
    // get normal
    //auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();

    FTensor::Tensor1<double, 3> t_rotation;
    t_rotation(i) = (t_coords(i) - tieCoord(i));
    ++t_coords;


    FTensor::Tensor1<double, 3> t_normal;
    if (rotationPlane == "x")
    {
      t_normal(0) = 1.0;
      t_normal(1) = 0.0;
      t_normal(2) = 0.0;
    }
    else if (rotationPlane == "y")
    {
      t_normal(0) = 0.0;
      t_normal(1) = 1.0;
      t_normal(2) = 0.0;
    }
    else if (rotationPlane == "z")
    {
      t_normal(0) = 0.0;
      t_normal(1) = 0.0;
      t_normal(2) = 1.0;
    }

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      // auto g = alpha * t_lambda;
      // ++t_lambda;

      t_sum_reaction(k) += alpha * t_lambda * (levi_civita(i, j, k) * t_rotation(j)) * t_normal(i);
      // ++t_normal;
      ++t_lambda;

      // int rr = 0;
      // for (; rr != OpBase::nbRows / DIM; ++rr) {
      // sum_reaction += g;
      //++t_row_base;
      //}

      // for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      //  ++t_row_base;
    }

    constexpr int ind[] = {1, 2, 3};
    constexpr int ind0[] = {0};
    // std::cout << "lambdaPtr: " << *lambdaPtr << std::endl;
    // std::cout << "sum_reaction: " << sum_reaction << std::endl;
    // std::cout << "t_sum_reaction: " << t_sum_reaction(0) << " " <<
    // t_sum_reaction(1) << " " << t_sum_reaction(2) << std::endl; std::cout <<
    // "total_reaction_ptr: " << *totalReactionPtr << std::endl; set reaction
    // from LM
    CHKERR VecSetValues(*totalReactionPtr, 1, ind0, &sum_reaction, ADD_VALUES);
    // set reaction from stress
    CHKERR VecSetValues(*totalReactionPtr, 3, ind, &t_sum_reaction(0),
                        ADD_VALUES);

    MoFEMFunctionReturn(0);
  }
  boost::shared_ptr<SmartPetscObj<Vec>> totalReactionPtr;
  boost::shared_ptr<VectorDouble> lambdaPtr;
  boost::shared_ptr<MatrixDouble> stressPtr;
  boost::shared_ptr<Range> tieFacesPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  std::string rotationPlane;
};

template <int DIM>
struct OpTieTermConstrainDisplacementRhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainDisplacementRhs(std::string lambda_name,
                                    boost::shared_ptr<MatrixDouble> u_ptr,
                                    FTensor::Tensor1<double, 3> tie_coord,
                                    FTensor::Tensor1<double, 3> tie_direction,
                                    boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);
    FTENSOR_INDEX(SPACE_DIM, m);

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      
      double t_rotation = pow(pow(t_coords(1) - tieCoord(1), 2.0) + pow(t_coords(2) - tieCoord(2),2.0), 0.5);


      FTensor::Tensor1<double, 3> t_delta_current;
      t_delta_current(i) =
          ((t_u(i)) - ((tieDirection(i) * time)));

      FTensor::Tensor1<double, SPACE_DIM> g;
      g(0) = alpha * (t_u(0) - ((tieDirection(0) * time) - t_rotation * t_delta_current(2) + t_rotation * t_delta_current(1)));
      g(1) = alpha * (t_u(1) - ((tieDirection(1) * time) - t_rotation * t_delta_current(0) + t_rotation * t_delta_current(2)));
      g(2) = alpha * (t_u(2) - ((tieDirection(2) * time) - t_rotation * t_delta_current(1) + t_rotation * t_delta_current(0)));
                               

      ++t_coords;
      ++t_u;

      auto t_nf = getFTensor1FromArray<DIM, DIM>(OpBase::locF);

      int rr = 0;
      for (; rr != OpBase::nbRows / DIM; ++rr) {
        t_nf(0) += t_row_base * g(0);
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
  FTensor::Tensor1<double, 3> tieDirection;
};

template <int DIM>
struct OpTieTermConstrainDisplacementLhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {
  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainDisplacementLhs(std::string lambda_name,
                                std::string col_field_name,
                                boost::shared_ptr<MatrixDouble> u_ptr,
                                FTensor::Tensor1<double, 3> tie_coord,
                                FTensor::Tensor1<double, 3> tie_direction,
                                boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, col_field_name, OpBase::OPROWCOL, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
    this->assembleTranspose = true;
    this->sYmm = false;
  }

  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      double t_rotation = pow(pow(t_coords(1) - tieCoord(1),2.0) + pow(t_coords(2) - tieCoord(2),2.0), 0.5);

      ++t_coords;

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
      t_tangent(0, 0) = 1.0;
      t_tangent(0, 1) = -1.0 * t_rotation;
      t_tangent(0, 2) = 1.0 * t_rotation;

      t_tangent(1, 0) = 1.0 * t_rotation;
      t_tangent(1, 1) = 1.0;
      t_tangent(1, 2) = -1.0 * t_rotation;

      t_tangent(2, 0) = -1.0 * t_rotation;
      t_tangent(2, 1) = 1.0 * t_rotation;
      t_tangent(2, 2) = 1.0;

      t_tangent(i, j) = alpha * t_tangent(i,j);



      int rr = 0;
      for (; rr != OpBase::nbRows/SPACE_DIM; ++rr) {
        auto t_col_base = col_data.getFTensor0N(gg, 0);
        auto t_mat = getFTensor2FromArray<DIM, DIM, DIM>(OpBase::locMat, DIM * rr);
        for (int cc = 0; cc != OpBase::nbCols/SPACE_DIM; cc++) {
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
  FTensor::Tensor1<double, 3> tieDirection;
  double tsTime;
};

template <int DIM>
struct OpTieTermConstrainDisplacementRhs_du
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainDisplacementRhs_du(std::string field_name,
                                   std::string lambda_name,
                                   boost::shared_ptr<MatrixDouble> u_ptr,
                                   boost::shared_ptr<MatrixDouble> lambda_ptr,
                                   FTensor::Tensor1<double, 3> tie_coord,
                                   FTensor::Tensor1<double, 3> tie_direction,
                                   boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(field_name, field_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), lambdaPtr(lambda_ptr), tieCoord(tie_coord),
        tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function on rows
    auto t_row_base = data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get coordinate at integration points
    auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
    // get displacement
    auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
    // get lambda
    auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

    auto &nf = OpBase::locF;

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      double t_rotation = pow(pow(t_coords(1) - tieCoord(1),2.0) + pow(t_coords(2) - tieCoord(2),2.0), 0.5);
      ++t_coords;

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
      t_tangent(0, 0) = 1.0;
      t_tangent(0, 1) = -1.0 * t_rotation;
      t_tangent(0, 2) = 1.0 * t_rotation;

      t_tangent(1, 0) = 1.0 * t_rotation;
      t_tangent(1, 1) = 1.0;
      t_tangent(1, 2) = -1.0 * t_rotation;

      t_tangent(2, 0) = -1.0 * t_rotation;
      t_tangent(2, 1) = 1.0 * t_rotation;
      t_tangent(2, 2) = 1.0;

      //t_tangent(i, j) = t_tangent(i, j);


      auto t_nf = getFTensor1FromPtr<DIM>(&nf[0]);
      FTensor::Tensor1<double, SPACE_DIM> g;
      //g(i) = alpha * (t_tangent(j, i) * t_lambda(j));
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
  boost::shared_ptr<MatrixDouble> uPtr;
  boost::shared_ptr<MatrixDouble> lambdaPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  FTensor::Tensor1<double, 3> tieDirection;
};


// template <int DIM>
// struct OpTieTermConstrainDisplacementRhs
//     : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

//   using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

//   OpTieTermConstrainDisplacementRhs(std::string lambda_name,
//                                     boost::shared_ptr<MatrixDouble> u_ptr,
//                                     FTensor::Tensor1<double, 3> tie_coord,
//                                     FTensor::Tensor1<double, 3> tie_direction,
//                                     boost::shared_ptr<Range> tie_faces_ptr)
//       : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
//         uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction) {
//     std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
//     doEntities[MBEDGE] = true;
//     doEntities[MBTRI] = true;
//     doEntities[MBQUAD] = true;
//   }

//   MoFEMErrorCode iNtegrate(EntData &data) {
//     MoFEMFunctionBegin;
//     FTENSOR_INDEX(SPACE_DIM, i);
//     FTENSOR_INDEX(SPACE_DIM, j);
//     FTENSOR_INDEX(SPACE_DIM, k);
//     FTENSOR_INDEX(SPACE_DIM, l);
//     FTENSOR_INDEX(SPACE_DIM, m);

//     //double time = getTStime();
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
//     // get gradient of displacement
//     //auto t_grad_u = getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(*gradUPtr);

//     for (int gg = 0; gg != nb_integration_pts; gg++) {
//       const auto alpha = t_w * vol;
//       ++t_w;

      
//       FTensor::Tensor1<double, 3> t_rotation;
//       t_rotation(i) = (t_coords(i) - tieCoord(i));


//       FTensor::Tensor1<double, SPACE_DIM> g;
//       g(i) = alpha * (t_u(i) - (tieDirection(i) * time));
                               

//       ++t_coords;
//       ++t_u;

//       auto t_nf = getFTensor1FromArray<DIM, DIM>(OpBase::locF);

//       int rr = 0;
//       for (; rr != OpBase::nbRows / DIM; ++rr) {
//         t_nf(i) += t_row_base * g(i);
//         ++t_row_base;
//         ++t_nf;
//       }
//       for (; rr < OpBase::nbRowBaseFunctions; ++rr)
//         ++t_row_base;
//     }
//     MoFEMFunctionReturn(0);
//   }

// private:
//   boost::shared_ptr<MatrixDouble> uPtr;
//   FTensor::Tensor1<double, 3> tieCoord;
//   FTensor::Tensor1<double, 3> tieDirection;
// };

// template <int DIM>
// struct OpTieTermConstrainDisplacementLhs
//     : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {
//   using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

//   OpTieTermConstrainDisplacementLhs(std::string lambda_name,
//                                 std::string col_field_name,
//                                 boost::shared_ptr<MatrixDouble> u_ptr,
//                                 FTensor::Tensor1<double, 3> tie_coord,
//                                 FTensor::Tensor1<double, 3> tie_direction,
//                                 boost::shared_ptr<Range> tie_faces_ptr)
//       : OpBase(lambda_name, col_field_name, OpBase::OPROWCOL, tie_faces_ptr),
//         uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction) {
//     std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
//     doEntities[MBEDGE] = true;
//     doEntities[MBTRI] = true;
//     doEntities[MBQUAD] = true;
//     this->assembleTranspose = true;
//     this->sYmm = false;
//   }

//   MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
//                            EntitiesFieldData::EntData &col_data) {
//     MoFEMFunctionBegin;
//     FTENSOR_INDEX(SPACE_DIM, i);
//     FTENSOR_INDEX(SPACE_DIM, j);

//     constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

//     //double time = getTStime();
//     double time = 1.0;

//     auto nb_integration_pts = getGaussPts().size2();
//     // get element volume
//     const double vol = OpBase::getMeasure();
//     // get base function gradient on rows
//     auto t_row_base = row_data.getFTensor0N();
//     // get integration weights
//     auto t_w = OpBase::getFTensor0IntegrationWeight();
//     // get coordinate at integration points
//     auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
//     // get displacement
//     auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);

//     for (int gg = 0; gg != nb_integration_pts; gg++) {
//       const auto alpha = t_w * vol;
//       ++t_w;

//       FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
//       t_tangent(i, j) = alpha * t_kd(i , j);



//       int rr = 0;
//       for (; rr != OpBase::nbRows/SPACE_DIM; ++rr) {
//         auto t_col_base = col_data.getFTensor0N(gg, 0);
//         auto t_mat = getFTensor2FromArray<DIM, DIM, DIM>(OpBase::locMat, DIM * rr);
//         for (int cc = 0; cc != OpBase::nbCols/SPACE_DIM; cc++) {
//           t_mat(i, j) += (t_row_base * t_col_base) * t_tangent(i, j);
//           ++t_mat;
// 					++t_col_base;	
//         }
//         ++t_row_base;
//       }
//       for (; rr < OpBase::nbRowBaseFunctions; ++rr)
//         ++t_row_base;
//     }
//     MoFEMFunctionReturn(0);
//   }

// private:
//   boost::shared_ptr<MatrixDouble> uPtr;
//   FTensor::Tensor1<double, 3> tieCoord;
//   FTensor::Tensor1<double, 3> tieDirection;
//   double tsTime;
// };

// template <int DIM>
// struct OpTieTermConstrainDisplacementRhs_du
//     : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

//   using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

//   OpTieTermConstrainDisplacementRhs_du(std::string field_name,
//                                    std::string lambda_name,
//                                    boost::shared_ptr<MatrixDouble> u_ptr,
//                                    boost::shared_ptr<MatrixDouble> lambda_ptr,
//                                    FTensor::Tensor1<double, 3> tie_coord,
//                                    FTensor::Tensor1<double, 3> tie_direction,
//                                    boost::shared_ptr<Range> tie_faces_ptr)
//       : OpBase(field_name, field_name, OpBase::OPROW, tie_faces_ptr),
//         uPtr(u_ptr), lambdaPtr(lambda_ptr), tieCoord(tie_coord),
//         tieDirection(tie_direction) {
//     std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
//     doEntities[MBEDGE] = true;
//     doEntities[MBTRI] = true;
//     doEntities[MBQUAD] = true;
//   }

//   MoFEMErrorCode iNtegrate(EntData &data) {
//     MoFEMFunctionBegin;
//     FTENSOR_INDEX(SPACE_DIM, i);
//     FTENSOR_INDEX(SPACE_DIM, j);

//     constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

//     //double time = getTStime();
//     double time = 1.0;

//     auto nb_integration_pts = getGaussPts().size2();
//     // get element volume
//     const double vol = OpBase::getMeasure();
//     // get base function on rows
//     auto t_row_base = data.getFTensor0N();
//     // get integration weights
//     auto t_w = OpBase::getFTensor0IntegrationWeight();
//     // get coordinate at integration points
//     auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
//     // get displacement
//     auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
//     // get lambda
//     auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

//     auto &nf = OpBase::locF;

//     for (int gg = 0; gg != nb_integration_pts; gg++) {
//       const auto alpha = t_w * vol;
//       ++t_w;

//       FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
//       t_tangent(i, j) = t_kd(i , j);


//       auto t_nf = getFTensor1FromPtr<DIM>(&nf[0]);
//       FTensor::Tensor1<double, SPACE_DIM> g;
//       g(j) = alpha * (t_lambda(i) * t_tangent(i, j) );
//       ++t_lambda;
//       int rr = 0;
//       for (; rr != OpBase::nbRows / DIM; ++rr) {
//         t_nf(i) += t_row_base * g(i);
//         ++t_row_base;
//         ++t_nf;
//       }

//       for (; rr < OpBase::nbRowBaseFunctions; ++rr)
//         ++t_row_base;
//     }
//     MoFEMFunctionReturn(0);
//   }

// private:
//   boost::shared_ptr<MatrixDouble> uPtr;
//   boost::shared_ptr<MatrixDouble> lambdaPtr;
//   FTensor::Tensor1<double, 3> tieCoord;
//   FTensor::Tensor1<double, 3> tieDirection;
// };


template <int DIM>
struct OpTieTermConstrainDisplacementXRhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainDisplacementXRhs(std::string lambda_name,
                                     boost::shared_ptr<MatrixDouble> u_ptr,
                                     FTensor::Tensor1<double, 3> tie_coord,
                                     FTensor::Tensor1<double, 3> tie_direction,
                                     boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);
    FTENSOR_INDEX(SPACE_DIM, m);

    double time = getTStime();
    // double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      
      double t_rotation = pow(pow(t_coords(1) - tieCoord(1),2.0) + pow(t_coords(2) - tieCoord(2),2.0), 0.5);


      FTensor::Tensor1<double, 3> t_delta_current;
      t_delta_current(i) =
          ((t_u(i)) - ((tieDirection(i) * time)));

      
      double g = alpha * (t_u(0) - ((tieDirection(0) * time) - t_rotation * t_delta_current(2) + t_rotation * t_delta_current(1)));
      //double g = alpha * (t_u(0) - ((tieDirection(0) * time) - t_delta_current(2) + t_delta_current(1)));

      ++t_coords;
      ++t_u;

      int rr = 0;
      for (; rr != OpBase::nbRows; ++rr) {
        OpBase::locF[rr] += t_row_base * g;
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
  FTensor::Tensor1<double, 3> tieDirection;
};

template <int DIM>
struct OpTieTermConstrainDisplacementXLhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {
  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainDisplacementXLhs(std::string lambda_name,
                                std::string col_field_name,
                                boost::shared_ptr<MatrixDouble> u_ptr,
                                FTensor::Tensor1<double, 3> tie_coord,
                                FTensor::Tensor1<double, 3> tie_direction,
                                boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, col_field_name, OpBase::OPROWCOL, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
    this->assembleTranspose = true;
    this->sYmm = false;
  }

  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    double time = getTStime();
    // double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      double t_rotation = pow(pow(t_coords(1) - tieCoord(1),2.0) + pow(t_coords(2) - tieCoord(2),2.0), 0.5);
      ++t_coords; 

      FTensor::Tensor1<double, SPACE_DIM> t_tangent;
      t_tangent(0) = alpha;
      t_tangent(1) = -alpha * t_rotation;
      t_tangent(2) = alpha * t_rotation;

      int rr = 0;
      for (; rr != OpBase::nbRows; ++rr) {
        auto t_col_base = col_data.getFTensor0N(gg, 0);
        auto t_mat = getFTensor1FromPtr<SPACE_DIM>(&OpBase::locMat(rr, 0));
        for (int cc = 0; cc != OpBase::nbCols / SPACE_DIM; cc++) {
          t_mat(i) += (t_row_base * t_col_base) * t_tangent(i);
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
  FTensor::Tensor1<double, 3> tieDirection;
  double tsTime;
};

template <int DIM>
struct OpTieTermConstrainDisplacementXRhs_du
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainDisplacementXRhs_du(std::string field_name,
                                   std::string lambda_name,
                                   boost::shared_ptr<MatrixDouble> u_ptr,
                                   boost::shared_ptr<VectorDouble> lambda_ptr,
                                   FTensor::Tensor1<double, 3> tie_coord,
                                   FTensor::Tensor1<double, 3> tie_direction,
                                   boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(field_name, field_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), lambdaPtr(lambda_ptr), tieCoord(tie_coord),
        tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    double time = getTStime();
    // double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function on rows
    auto t_row_base = data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get coordinate at integration points
    auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
    // get displacement
    auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
    // get lambda
    auto t_lambda = getFTensor0FromVec(*lambdaPtr);

    auto &nf = OpBase::locF;

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      double t_rotation = pow(pow(t_coords(1) - tieCoord(1),2.0) + pow(t_coords(2) - tieCoord(2),2.0), 0.5);
      ++t_coords;

      FTensor::Tensor1<double, SPACE_DIM> t_tangent;
      t_tangent(0) = 1;
      t_tangent(1) = -1 * t_rotation;
      t_tangent(2) = 1 * t_rotation;

      auto t_nf = getFTensor1FromPtr<DIM>(&nf[0]);
      FTensor::Tensor1<double, SPACE_DIM> g;
      g(i) = alpha * t_lambda * (t_tangent(i));
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
  boost::shared_ptr<MatrixDouble> uPtr;
  boost::shared_ptr<VectorDouble> lambdaPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  FTensor::Tensor1<double, 3> tieDirection;
};

// template <int DIM>
// struct OpTieTermConstrainDisplacementGradXRhs
//     : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

//   using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

//   OpTieTermConstrainDisplacementGradXRhs(std::string lambda_name,
//                                      boost::shared_ptr<MatrixDouble> u_ptr,
//                                      FTensor::Tensor1<double, 3> tie_coord,
//                                      FTensor::Tensor1<double, 3> tie_direction,
//                                      boost::shared_ptr<Range> tie_faces_ptr,
//                                      boost::shared_ptr<MatrixDouble> grad_u_ptr)
//       : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
//         uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction), gradUPtr(grad_u_ptr) {
//     std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
//     //doEntities[MBEDGE] = true;
//     doEntities[MBTRI] = true;
//     doEntities[MBQUAD] = true;
//   }

//   MoFEMErrorCode iNtegrate(EntData &data) {
//     MoFEMFunctionBegin;
//     FTENSOR_INDEX(SPACE_DIM, i);
//     FTENSOR_INDEX(SPACE_DIM, j);
//     FTENSOR_INDEX(SPACE_DIM, k);
//     FTENSOR_INDEX(SPACE_DIM, l);
//     FTENSOR_INDEX(SPACE_DIM, m);

//     double time = getTStime();
//     // double time = 1.0;

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
//     // get gradient of displacement
//     auto t_grad_u = getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(*gradUPtr);

//     for (int gg = 0; gg != nb_integration_pts; gg++) {
//       const auto alpha = t_w * vol;
//       ++t_w;

//       FTensor::Tensor1<double, 3> t_r;
//       t_r(i) = (t_coords(i) - tieCoord(i));

//       FTensor::Tensor1<double, SPACE_DIM> t_theta;
//       //t_theta(0) = 0.5 * (t_grad_u(2, 1) - t_grad_u(1, 2));
//       //t_theta(1) = 0.5 * (t_grad_u(0, 2) - t_grad_u(2, 0));
//       //t_theta(2) = 0.5 * (t_grad_u(1, 0) - t_grad_u(0, 1));
//       t_theta(k) = 0.25 * levi_civita(i , j, k) * (t_grad_u(i , j) - t_grad_u(j, i));

//       FTensor::Tensor1<double, SPACE_DIM> g;
//       g(i) = alpha * (t_u(i) - ((tieDirection(i) * time) + (levi_civita(i,j,k) * t_r(j) * t_theta(k))));

//       ++t_coords;
//       ++t_u;
//       ++t_grad_u;

//       int rr = 0;
//       for (; rr != OpBase::nbRows; ++rr) {
//         OpBase::locF[rr] += t_row_base * g(0);
//         ++t_row_base;
//       }
//       for (; rr < OpBase::nbRowBaseFunctions; ++rr)
//         ++t_row_base;

//     }
//     MoFEMFunctionReturn(0);
//   }

// private:
//   boost::shared_ptr<MatrixDouble> uPtr;
//   FTensor::Tensor1<double, 3> tieCoord;
//   FTensor::Tensor1<double, 3> tieDirection;
//   boost::shared_ptr<MatrixDouble> gradUPtr;
// };

// template <int DIM>
// struct OpTieTermConstrainDisplacementGradXLhs
//     : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {
//   using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

//   OpTieTermConstrainDisplacementGradXLhs(std::string lambda_name,
//                                 std::string col_field_name,
//                                 boost::shared_ptr<MatrixDouble> u_ptr,
//                                 FTensor::Tensor1<double, 3> tie_coord,
//                                 FTensor::Tensor1<double, 3> tie_direction,
//                                 boost::shared_ptr<Range> tie_faces_ptr)
//       : OpBase(lambda_name, col_field_name, OpBase::OPROWCOL, tie_faces_ptr),
//         uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction) {
//     std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
//     //doEntities[MBEDGE] = true;
//     doEntities[MBTRI] = true;
//     doEntities[MBQUAD] = true;
//     this->assembleTranspose = true;
//     this->sYmm = false;
//   }

//   MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
//                            EntitiesFieldData::EntData &col_data) {
//     MoFEMFunctionBegin;
//     FTENSOR_INDEX(SPACE_DIM, i);
//     FTENSOR_INDEX(SPACE_DIM, j);
//     FTENSOR_INDEX(SPACE_DIM, k);

//     constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

//     double time = getTStime();
//     // double time = 1.0;

//     auto nb_integration_pts = getGaussPts().size2();
//     // get element volume
//     const double vol = OpBase::getMeasure();
//     // get base function gradient on rows
//     auto t_row_base = row_data.getFTensor0N();
//     // get integration weights
//     auto t_w = OpBase::getFTensor0IntegrationWeight();
//     // get coordinate at integration points
//     auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
//     // get displacement
//     auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);

//     for (int gg = 0; gg != nb_integration_pts; gg++) {
//       const auto alpha = t_w * vol;
//       ++t_w;

//       FTensor::Tensor1<double, 3> t_r;
//       t_r(i) = (t_coords(i) - tieCoord(i));
//       ++t_coords;
     
//       int rr = 0;
//       for (; rr != OpBase::nbRows; ++rr) {
//         auto t_col_base = col_data.getFTensor0N(gg, 0);
//         auto t_col_diff_base = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);
//         auto t_mat = getFTensor1FromPtr<SPACE_DIM>(&OpBase::locMat(rr, 0));
//         for (int cc = 0; cc != OpBase::nbCols / SPACE_DIM; cc++) {
//           FTensor::Tensor1<double, SPACE_DIM> t_theta;
//           t_theta(i) = 0.25 * levi_civita(i, j, k) * (t_col_diff_base(j) - t_col_diff_base(i));
//           t_mat(i) += (t_row_base * t_col_base) * (alpha * (t_kd(0,1) - (levi_civita(i, j, k) * t_r(j) * t_theta(k))));
//           ++t_mat;
//           ++t_col_base;
//           ++t_col_diff_base;
//         }
//         ++t_row_base;
//       }
//       for (; rr < OpBase::nbRowBaseFunctions; ++rr)
//         ++t_row_base;

//     }
//     MoFEMFunctionReturn(0);
//   }

// private:
//   boost::shared_ptr<MatrixDouble> uPtr;
//   FTensor::Tensor1<double, 3> tieCoord;
//   FTensor::Tensor1<double, 3> tieDirection;
//   double tsTime;
// };

// template <int DIM>
// struct OpTieTermConstrainDisplacementGradXRhs_du
//     : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

//   using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

//   OpTieTermConstrainDisplacementGradXRhs_du(std::string field_name,
//                                    std::string lambda_name,
//                                    boost::shared_ptr<MatrixDouble> u_ptr,
//                                    boost::shared_ptr<VectorDouble> lambda_ptr,
//                                    FTensor::Tensor1<double, 3> tie_coord,
//                                    FTensor::Tensor1<double, 3> tie_direction,
//                                    boost::shared_ptr<Range> tie_faces_ptr)
//       : OpBase(field_name, field_name, OpBase::OPROW, tie_faces_ptr),
//         uPtr(u_ptr), lambdaPtr(lambda_ptr), tieCoord(tie_coord),
//         tieDirection(tie_direction) {
//     std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
//     //doEntities[MBEDGE] = true;
//     doEntities[MBTRI] = true;
//     doEntities[MBQUAD] = true;
//   }

//   MoFEMErrorCode iNtegrate(EntData &data) {
//     MoFEMFunctionBegin;
//     FTENSOR_INDEX(SPACE_DIM, i);
//     FTENSOR_INDEX(SPACE_DIM, j);
//     FTENSOR_INDEX(SPACE_DIM, k);

//     constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

//     double time = getTStime();
//     // double time = 1.0;

//     auto nb_integration_pts = getGaussPts().size2();
//     // get element volume
//     const double vol = OpBase::getMeasure();
//     // get base function on rows
//     auto t_row_base = data.getFTensor0N();
//     // get integration weights
//     auto t_w = OpBase::getFTensor0IntegrationWeight();
//     // get coordinate at integration points
//     auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
//     // get displacement
//     auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
//     // get lambda
//     auto t_lambda = getFTensor0FromVec(*lambdaPtr);

    

//     auto &nf = OpBase::locF;

//     for (int gg = 0; gg != nb_integration_pts; gg++) {
//       const auto alpha = t_w * vol;
//       ++t_w;

//       FTensor::Tensor1<double, SPACE_DIM, SPACE_DIM> t_rot;
//       t_r(i) = (t_coords(i) - tieCoord(i));




//       auto t_nf = getFTensor1FromPtr<DIM>(&nf[0]);
//       auto t_row_diff_base = data.getFTensor1DiffN<SPACE_DIM>(gg, 0);

//       int rr = 0;
//       for (; rr != OpBase::nbRows / DIM; ++rr) {
//         FTensor::Tensor1<double, SPACE_DIM> t_theta;
//         t_theta(i) = 0.25 * levi_civita(i, j, k) * (t_row_diff_base(j) - t_row_diff_base(i));
//         t_nf(i) += t_row_base * (alpha * (t_kd(0, i) - (levi_civita(i, j, k) * t_r(j) * t_theta(k))));
//         ++t_row_base;
//         ++t_nf;
//         ++t_row_diff_base;
//       }

//       for (; rr < OpBase::nbRowBaseFunctions; ++rr) {
//         ++t_row_base;
//         ++t_row_diff_base;
//       }
//     }
//     MoFEMFunctionReturn(0);
//   }

// private:
//   boost::shared_ptr<MatrixDouble> uPtr;
//   boost::shared_ptr<VectorDouble> lambdaPtr;
//   FTensor::Tensor1<double, 3> tieCoord;
//   FTensor::Tensor1<double, 3> tieDirection;
// };



template <int DIM>
struct OpCalculateTieDisplacementReactionForce
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpCalculateTieDisplacementReactionForce(
      std::string lambda_name, boost::shared_ptr<MatrixDouble> lambda_ptr,
      boost::shared_ptr<SmartPetscObj<Vec>> total_reaction_ptr,
      boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        lambdaPtr(lambda_ptr), totalReactionPtr(total_reaction_ptr),
        tieFacesPtr(tie_faces_ptr) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    // doEntities[MBEDGE] = true;
    // doEntities[MBTRI] = true;
    // doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    if (tieFacesPtr->find(getNumeredEntFiniteElementPtr()->getEnt()) ==
        tieFacesPtr->end()) {
      MoFEMFunctionReturnHot(0);
    }
    // std::cout<< "Entity id = "<< getNumeredEntFiniteElementPtr()->getEnt() <<
    // std::endl;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);

    double sum_reaction = 0.0;
    FTensor::Tensor1<double, 3> t_sum_reaction{0.0, 0.0, 0.0};

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
    // get stress
    // auto t_stress = getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(*stressPtr);
    // get normal
    // auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;  

      t_sum_reaction(i) += alpha * t_lambda(i);
      ++t_lambda;
      // int rr = 0;
      // for (; rr != OpBase::nbRows / DIM; ++rr) {
      //sum_reaction += g;
      //++t_row_base;
      //}

      // for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      //  ++t_row_base;
    }

    constexpr int ind[] = {1, 2, 3};
    constexpr int ind0[] = {0};
    // std::cout << "lambdaPtr: " << *lambdaPtr << std::endl;
    // std::cout << "sum_reaction: " << sum_reaction << std::endl;
    // std::cout << "t_sum_reaction: " << t_sum_reaction(0) << " " <<
    // t_sum_reaction(1) << " " << t_sum_reaction(2) << std::endl; std::cout <<
    // "total_reaction_ptr: " << *totalReactionPtr << std::endl; set reaction
    // from LM
    CHKERR VecSetValues(*totalReactionPtr, 1, ind0, &sum_reaction, ADD_VALUES);
    // set reaction from stress
    CHKERR VecSetValues(*totalReactionPtr, 3, ind, &t_sum_reaction(0),
                        ADD_VALUES);

    MoFEMFunctionReturn(0);
  }
  boost::shared_ptr<SmartPetscObj<Vec>> totalReactionPtr;
  boost::shared_ptr<MatrixDouble> lambdaPtr;
  boost::shared_ptr<MatrixDouble> stressPtr;
  boost::shared_ptr<Range> tieFacesPtr;
};

template <int DIM>
struct OpExtractReferencePoint
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpExtractReferencePoint(std::string field_name,
                          boost::shared_ptr<MatrixDouble> u_ptr,
                          boost::shared_ptr<MatrixDouble> vel_ptr,
                          boost::shared_ptr<VectorDouble> u_ref_ptr,
                          boost::shared_ptr<VectorDouble> vel_ref_ptr,
                          FTensor::Tensor1<double, 3> tie_coord,
                          boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(field_name, field_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), velPtr(vel_ptr), uRefPtr(u_ref_ptr),
        velRefPtr(vel_ref_ptr), tieCoord(tie_coord),
        tieFacesPtr(tie_faces_ptr) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    // doEntities[MBEDGE] = true;
    // doEntities[MBTRI] = true;
    // doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    if (tieFacesPtr->find(getNumeredEntFiniteElementPtr()->getEnt()) ==
        tieFacesPtr->end()) {
      MoFEMFunctionReturnHot(0);
    }

    std::cout << "Entity type = " << type << std::endl;
    // get entity id

    std::cout<< "Entity id = "<< getNumeredEntFiniteElementPtr()->getEnt()<< std::endl;

    // std::cout<< "Entity id = "<< getNumeredEntFiniteElementPtr()->getEnt() <<
    // std::endl;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);

    auto nb_integration_pts = OpBase::getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function on rows
    auto t_row_base = data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get coordinate at integration points
    auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
    // get displacement
    auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
    // get velocity
    auto t_vel = getFTensor1FromMat<SPACE_DIM>(*velPtr);



    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      (*uRefPtr)(0) = t_u(0);
      (*uRefPtr)(1) = t_u(1);
      (*uRefPtr)(2) = t_u(2);

      (*velRefPtr)(0) = t_vel(0);
      (*velRefPtr)(1) = t_vel(1);
      (*velRefPtr)(2) = t_vel(2);

      ++t_u;
      ++t_vel;

    }
    MoFEMFunctionReturn(0);
  }
  boost::shared_ptr<MatrixDouble> uPtr;
  boost::shared_ptr<MatrixDouble> velPtr;
  boost::shared_ptr<VectorDouble> uRefPtr;
  boost::shared_ptr<VectorDouble> velRefPtr;
  boost::shared_ptr<Range> tieFacesPtr;
  FTensor::Tensor1<double, 3> tieCoord;

};


// Testing Global Rigid body constraint

template <int DIM>
struct OpTieTermConstrainRigidBodyRhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRigidBodyRhs(
      std::string lambda_name, boost::shared_ptr<MatrixDouble> u_ptr,
      FTensor::Tensor1<double, 3> tie_coord,
      FTensor::Tensor1<double, 3> tie_direction,
      boost::shared_ptr<Range> tie_faces_ptr,
      boost::shared_ptr<VectorDouble> translation_ptr,
      boost::shared_ptr<VectorDouble> theta_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction),
        translationPtr(translation_ptr), thetaPtr(theta_ptr) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;

    //std::cout << "getNumberedEntities: " << *getNumeredEntFiniteElementPtr() << std::endl;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);
    FTENSOR_INDEX(SPACE_DIM, m);

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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

    //std::cout << "t_translation: " << t_translation(0) << " " << t_translation(1) << " " << t_translation(2) << std::endl;

    FTensor::Tensor1<double, 3> t_theta;
    t_theta(0) = (*thetaPtr)(0);
    t_theta(1) = (*thetaPtr)(1);
    t_theta(2) = (*thetaPtr)(2);

    //std::cout << "t_theta: " << t_theta(0) << " " << t_theta(1) << " " << t_theta(2) << std::endl;

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      // define rotation matrix
      FTensor::Tensor2<double, 3, 3> t_omega;
      t_omega(i, j) = levi_civita(i, j, k) * t_theta(k);

      // define reference position
      FTensor::Tensor1<double, 3> t_x_ref;
      t_x_ref(i) = t_u(i) - t_translation(i) - t_omega(i, j) * (t_coords(j) - tieCoord(j));

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
  FTensor::Tensor1<double, 3> tieDirection;
  boost::shared_ptr<VectorDouble> translationPtr;
  boost::shared_ptr<VectorDouble> thetaPtr;
};
template <int DIM>
struct OpTieTermConstrainRigidBodyLhs_dU
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {
  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRigidBodyLhs_dU(
      std::string lambda_name, std::string col_field_name,
      boost::shared_ptr<MatrixDouble> u_ptr,
      FTensor::Tensor1<double, 3> tie_coord,
      FTensor::Tensor1<double, 3> tie_direction,
      boost::shared_ptr<Range> tie_faces_ptr,
      boost::shared_ptr<VectorDouble> translation_ptr,
      boost::shared_ptr<VectorDouble> theta_ptr)
      : OpBase(lambda_name, col_field_name, OpBase::OPROWCOL, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction),
        translationPtr(translation_ptr), thetaPtr(theta_ptr) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
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

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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

      // define rotation matrix
      FTensor::Tensor2<double, 3, 3> t_omega;
      t_omega(i, j) = levi_civita(i, j, k) * t_theta(k);

      // define reference position
      FTensor::Tensor1<double, 3> t_x_ref;
      t_x_ref(i) = t_u(i) - t_translation(i) - t_omega(i, j) * (t_coords(j) - tieCoord(j));

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
      t_tangent(i, j) = alpha * t_kd(i, j);

      int rr = 0;
      for (; rr != OpBase::nbRows/SPACE_DIM; ++rr) {
        auto t_col_base = col_data.getFTensor0N(gg, 0);
        auto t_mat = getFTensor2FromArray<DIM, DIM, DIM>(OpBase::locMat, DIM * rr);
        for (int cc = 0; cc != OpBase::nbCols/SPACE_DIM; cc++) {
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
  FTensor::Tensor1<double, 3> tieDirection;
  boost::shared_ptr<VectorDouble> translationPtr;
  boost::shared_ptr<VectorDouble> thetaPtr;
};

template <int DIM>
struct OpTieTermConstrainRigidBodyLhs_dTranslation
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {
  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRigidBodyLhs_dTranslation(
      std::string lambda_name, std::string col_field_name,
      boost::shared_ptr<MatrixDouble> u_ptr,
      FTensor::Tensor1<double, 3> tie_coord,
      FTensor::Tensor1<double, 3> tie_direction,
      boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, col_field_name, OpBase::OPROWCOL, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
    //doEntities[MBENTITYSET] = true;
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
    //std::cout << "type in dTrans: = " << getNumeredEntFiniteElementPtr()->getEntType() << std::endl;

    auto &locMat = OpBase::locMat;

    // double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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

    // dummy base
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
        //++t_mat;
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
  FTensor::Tensor1<double, 3> tieDirection;
};

template <int DIM>
struct OpTieTermConstrainRigidBodyLhs_dRotation
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {
  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRigidBodyLhs_dRotation(
      std::string lambda_name, std::string col_field_name,
      boost::shared_ptr<MatrixDouble> u_ptr,
      FTensor::Tensor1<double, 3> tie_coord,
      FTensor::Tensor1<double, 3> tie_direction,
      boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, col_field_name, OpBase::OPROWCOL, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
    doEntities[MBENTITYSET] = true;
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

    //double time = getTStime();
    double time = 1.0;
    //std::cout << "type: = " << getNumeredEntFiniteElementPtr()->getEntType() << std::endl;

    auto nb_integration_pts = getGaussPts().size2();
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
    

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
      t_tangent(i, k) = alpha * levi_civita(i,j,k) * (t_coords(j) - tieCoord(j));
      ++t_coords;

      
      int rr = 0;
      for (; rr != OpBase::nbRows / SPACE_DIM; ++rr) {
        auto t_mat = getFTensor2FromPtr<SPACE_DIM, SPACE_DIM>(
            &OpBase::locMat(SPACE_DIM * rr, 0));
        t_mat(i, j) -= (t_row_base * t_tangent(i, j));
        //++t_mat;
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
  FTensor::Tensor1<double, 3> tieDirection;
};

template <int DIM>
struct OpTieTermConstrainRigidBodyRhs_du
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRigidBodyRhs_du(std::string field_name,
                                   std::string lambda_name,
                                   boost::shared_ptr<MatrixDouble> u_ptr,
                                   boost::shared_ptr<MatrixDouble> lambda_ptr,
                                   FTensor::Tensor1<double, 3> tie_coord,
                                   FTensor::Tensor1<double, 3> tie_direction,
                                   boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(field_name, field_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), lambdaPtr(lambda_ptr), tieCoord(tie_coord),
        tieDirection(tie_direction) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
    // get element volume
    const double vol = OpBase::getMeasure();
    // get base function on rows
    auto t_row_base = data.getFTensor0N();
    // get integration weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();
    // get coordinate at integration points
    auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
    // get displacement
    auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
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
      // std::cout << "g: " << g(0) << " " << g(1) << " " << g(2) << std::endl;
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
  boost::shared_ptr<MatrixDouble> lambdaPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  FTensor::Tensor1<double, 3> tieDirection;
};

template <int DIM>
struct OpTieTermConstrainRigidBodyGlobalTranslationRhs
    : public ForcesAndSourcesCore::UserDataOperator{

      using OpUserDataOp = ForcesAndSourcesCore::UserDataOperator;

      OpTieTermConstrainRigidBodyGlobalTranslationRhs(
          std::string lambda_name, FTensor::Tensor1<double, 3> tie_coord,
          FTensor::Tensor1<double, 3> tie_direction,
          boost::shared_ptr<Range> tie_faces_ptr,
          boost::shared_ptr<MatrixDouble> lambda_ptr,
          boost::shared_ptr<VectorDouble> int_translation_ptr)
          : OpUserDataOp(lambda_name, OpUserDataOp::OPROW),
            tieCoord(tie_coord), tieDirection(tie_direction),
            lambdaPtr(lambda_ptr), intTranslationPtr(int_translation_ptr) {
        std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
        // doEntities[MBEDGE] = true;
        // doEntities[MBTRI] = true;
        // doEntities[MBQUAD] = true;
        doEntities[MBENTITYSET] = true;
      }
  MoFEMErrorCode doWork(int side, EntityType type,
                         EntitiesFieldData::EntData &data) {
  // MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;

    if(data.getIndices().empty())
      MoFEMFunctionReturnHot(0);

    //std::cout << "getNumberedEntities: " << *getNumeredEntFiniteElementPtr() << std::endl;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);
    FTENSOR_INDEX(SPACE_DIM, m);

    // double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();

    //auto t_w = OpUserDataOp::getFTensor0IntegrationWeight();

    //const double vol = OpUserDataOp::getMeasure();
    
    //auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

    auto v_dofs = data.getFieldData();
    auto &vec = data.getFieldData();
    //std::cout << "v_dofs: " << v_dofs(0) << " " << v_dofs(1) << " " << v_dofs(2) << std::endl;
    auto ind = data.getIndices(); 
    auto local_ind = data.getLocalIndices();
    //auto t_dofs = data.getFTensor1FieldData<3>();
    auto t_dofs = getFTensor1FromPtr<SPACE_DIM>(&data.getFieldEntities()[0]->getEntFieldData()[0]);

    double* intTranslationRawPtr = &(*intTranslationPtr->data().begin());

    auto t_int_translation = getFTensor1FromPtr<SPACE_DIM>(intTranslationRawPtr);

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

    
    //for (int gg = 0; gg != nb_integration_pts; gg++) {
      //const auto alpha = t_w * vol;
      //++t_w;

      //FTensor::Tensor1<double, SPACE_DIM> t_dofs;
      //g(i) = t_lambda(i);
      //g(i) = 1.0;

      //++t_lambda;

      //auto t_nf = getFTensor1FromArray<3, 0>(OpUserDataOp::locF);

      auto t_dummy_base = getFTensor1FromMat<3>(m_kd);
      int rr = 0;
      for (; rr != 3; ++rr) {
        // std::cout << "t_dofs: " << t_dofs(0) << " " << t_dofs(1) << " " << t_dofs(2) << std::endl;
        // std::cout << "t_dummy_base: " << t_dummy_base(0) << " " << t_dummy_base(1) << " " << t_dummy_base(2) << std::endl;
        // std::cout << "t_int_translation: " << t_int_translation(0) << " " << t_int_translation(1) << " " << t_int_translation(2) << std::endl;
        t_dofs(j) -= (t_dummy_base(i) * t_int_translation(i)) * t_dummy_base(j);
        ++t_dummy_base;
        //++t_nf;
      }      
      //std::cout << "Translation: " << t_dofs(0) << " " << t_dofs(1) << " " << t_dofs(2) << std::endl;    
    //}
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> uPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  FTensor::Tensor1<double, 3> tieDirection;
  boost::shared_ptr<MatrixDouble> lambdaPtr;
  boost::shared_ptr<VectorDouble> intTranslationPtr;
};

template <int DIM>
struct OpTieTermConstrainRigidBodyGlobalTranslationIntegralRhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

      using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRigidBodyGlobalTranslationIntegralRhs(
      std::string lambda_name, boost::shared_ptr<MatrixDouble> u_ptr,
      FTensor::Tensor1<double, 3> tie_coord,
      FTensor::Tensor1<double, 3> tie_direction,
      boost::shared_ptr<Range> tie_faces_ptr,
      boost::shared_ptr<MatrixDouble> lambda_ptr,
      boost::shared_ptr<VectorDouble> int_translation_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction), lambdaPtr(lambda_ptr),
        intTranslationPtr(int_translation_ptr){
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
    //doEntities[MBENTITYSET] = true;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                         EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;

    // if (data.getIndices().empty())
    //   MoFEMFunctionReturn(0);

    //std::cout << "getNumberedEntities: " << *getNumeredEntFiniteElementPtr() << std::endl;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);
    FTENSOR_INDEX(SPACE_DIM, m);

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
    
    auto &vec = *intTranslationPtr;
    
    
    // if (data.getNumeredEntFiniteElementPtr() == 0) {
      // vec.clear();
    // }
    //std::cout << "type: "<< type << std::endl;

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

    auto nb_dofs = data.getIndices().size();

    double* intTranslationRawPtr = &(*intTranslationPtr->data().begin());

    auto t_int_translation = getFTensor1FromPtr<3>(intTranslationRawPtr);


    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor1<double, 3> g;
      g(i) = alpha * t_lambda(i);
                               
      ++t_lambda;

      t_int_translation(i) += g(i);

      // int rr = 0;
      // for (; rr != OpBase::nbRows/DIM; ++rr) {
      //   std::cout << "t_row_base: " << t_row_base << std::endl;
      //   t_int_translation(i) += t_row_base * g(i);
      //   ++t_row_base;
      // }
      // for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      //   ++t_row_base;
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> uPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  FTensor::Tensor1<double, 3> tieDirection;
  boost::shared_ptr<VectorDouble> intTranslationPtr;
  boost::shared_ptr<MatrixDouble> lambdaPtr;
  boost::shared_ptr<VectorDouble> thetaPtr;
};

template <int DIM>
struct OpTieTermConstrainRigidBodyGlobalRotationIntegralRhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {
      
      using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRigidBodyGlobalRotationIntegralRhs(
      std::string lambda_name, boost::shared_ptr<MatrixDouble> u_ptr,
      FTensor::Tensor1<double, 3> tie_coord,
      FTensor::Tensor1<double, 3> tie_direction,
      boost::shared_ptr<Range> tie_faces_ptr,
      boost::shared_ptr<MatrixDouble> lambda_ptr,
      boost::shared_ptr<MatrixDouble> int_rotation_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction), lambdaPtr(lambda_ptr),
        intRotationPtr(int_rotation_ptr){
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;

    //std::cout << "getNumberedEntities: " << *getNumeredEntFiniteElementPtr() << std::endl;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);
    FTENSOR_INDEX(SPACE_DIM, m);

    //double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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
    double* intRotationRawPtr = &(*intRotationPtr->data().begin());

    auto t_int_rotation = getFTensor2FromPtr<SPACE_DIM, SPACE_DIM>(intRotationRawPtr);


    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> g;
      g(i, j) = alpha * t_lambda(i) * (t_coords(j) - tieCoord(j));

      ++t_lambda;
      ++t_coords;
      t_int_rotation(i, j) += g(i, j);

      // int rr = 0;
      // for (; rr != OpBase::nbRows / DIM; ++rr) {
      //   t_int_rotation(i, j) += t_row_base * g(i, j);
      //   ++t_row_base;
      // }
      // for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      //   ++t_row_base;
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> uPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  FTensor::Tensor1<double, 3> tieDirection;
  boost::shared_ptr<MatrixDouble> lambdaPtr;
  boost::shared_ptr<VectorDouble> translationPtr;
  boost::shared_ptr<MatrixDouble> intRotationPtr;
};

template <int DIM>
struct OpTieTermConstrainRigidBodyGlobalTranslationLhs_dLambda
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRigidBodyGlobalTranslationLhs_dLambda(
      std::string lambda_name, std::string field_name,
       boost::shared_ptr<MatrixDouble> u_ptr,
      FTensor::Tensor1<double, 3> tie_coord,
      FTensor::Tensor1<double, 3> tie_direction,
      boost::shared_ptr<Range> tie_faces_ptr)
      : OpBase(lambda_name, field_name, OpBase::OPROWCOL, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction){
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    // doEntities[MBEDGE] = true;
    // doEntities[MBTRI] = true;
    // doEntities[MBQUAD] = true;
    // doEntities[MBENTITYSET] = true;
    this->assembleTranspose = true;
    this->sYmm = false;
  }
  MoFEMErrorCode iNtegrate(EntData &data) {
  // MoFEMErrorCode doWork(int side, EntityType type,
  //                       EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);
    FTENSOR_INDEX(SPACE_DIM, m);

    constexpr auto t_kd = FTensor::Kronecker_Delta<double>();

    // double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();

    auto t_w = OpBase::getFTensor0IntegrationWeight();

    const double vol = OpBase::getMeasure();
    
    auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

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

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
      t_tangent(i, j) = alpha * t_kd(i, j);

      auto t_dummy_base = getFTensor1FromMat<3>(m_kd);

      int rr = 0;
      for (; rr != 3; ++rr) {
        auto t_mat = getFTensor1FromPtr<SPACE_DIM>(&OpBase::locMat(rr, 0));
        auto t_col_base = data.getFTensor0N(gg, 0);
        for (int cc = 0; cc != OpBase::nbCols/SPACE_DIM; cc++) {
          t_mat(j) += t_col_base * t_dummy_base(i) * t_tangent(i, j);
          ++t_col_base;
          ++t_mat;
        }
        ++t_dummy_base;
      }
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> uPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  FTensor::Tensor1<double, 3> tieDirection;
  boost::shared_ptr<MatrixDouble> lambdaPtr;
};

template <int DIM>
struct OpTieTermConstrainRigidBodyGlobalRotationRhs
    : public ForcesAndSourcesCore::UserDataOperator {

  using OpUserDataOp = ForcesAndSourcesCore::UserDataOperator;

  OpTieTermConstrainRigidBodyGlobalRotationRhs(
      std::string lambda_name, FTensor::Tensor1<double, 3> tie_coord,
      FTensor::Tensor1<double, 3> tie_direction,
      boost::shared_ptr<Range> tie_faces_ptr,
      boost::shared_ptr<MatrixDouble> lambda_ptr,
      boost::shared_ptr<MatrixDouble> int_rotation_ptr)
      : OpUserDataOp(lambda_name, OpUserDataOp::OPROW),
        tieCoord(tie_coord), tieDirection(tie_direction),
        lambdaPtr(lambda_ptr), intRotationPtr(int_rotation_ptr) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    // doEntities[MBEDGE] = true;
    // doEntities[MBTRI] = true;
    // doEntities[MBQUAD] = true;
    doEntities[MBENTITYSET] = true;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    // MoFEMErrorCode iNtegrate(EntData &data) {

    // check if entity is first in the range
    //std::cout << "getNumberedEntities: " << getNumeredEntFiniteElementPtr() << std::endl;

    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);

    // double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
    // get element volume
    // const double vol = OpUserDataOp::getMeasure();
    // // get base function gradient on rows
    // auto t_row_base = data.getFTensor0N();
    // // get integration weights
    // auto t_w = OpUserDataOp::getFTensor0IntegrationWeight();
    // get coordinate at integration points
    //auto t_coords = OpUserDataOp::getFTensor1CoordsAtGaussPts();
    //auto t_coords = OpUserDataOp::getFTensor1Coords();
    // get displacement
    // auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
    // // get lambda
    // auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

    //auto t_dofs = data.getFTensor1FieldData<3>();
    //auto &vec = data.get();
    auto t_dofs = getFTensor1FromPtr<SPACE_DIM>(&data.getFieldEntities()[0]->getEntFieldData()[0]);

    double* intRotationRawPtr = &(*intRotationPtr->data().begin());
    auto t_int_rotation = getFTensor2FromPtr<SPACE_DIM, SPACE_DIM>(intRotationRawPtr);
  

    //data_ptr = &*data.getFieldData().data().begin();

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

    //FTensor::Tensor1<double, 3> t_dofs;

    auto t_dummy_base = getFTensor1FromMat<3>(m_kd);

    // for (int gg = 0; gg != nb_integration_pts; gg++) {
    //   const auto alpha = t_w * vol;
    //   ++t_w;

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> g;
      // g(i, j) = alpha * t_lambda(i) *
      //           (t_coords(j) - tieCoord(j));
      g(i, j) = 1.0;

      //++ t_lambda;
      //++t_coords;
      //++t_u;
      

      //auto t_nf = getFTensor1FromArray<DIM, DIM>(OpUserDataOp::locF);
      int rr = 0;
      for (; rr != 3; ++rr) {
        t_dofs(l) -= (t_dummy_base(k) * levi_civita(i, j, k) * t_int_rotation(i, j)) * t_dummy_base(l);

        //t_nf(k) += t_dummy_base(k) * levi_civita(i, j, k) * g(i, j);
        ++t_dummy_base;
        //++t_nf;
      }
      // vec(0) = t_dofs(0);
      // vec(1) = t_dofs(1);
      // vec(2) = t_dofs(2);
      //std::cout << "Rotation: " << t_dofs(0) << " " << t_dofs(1) << " " << t_dofs(2) << std::endl;    
    // }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> uPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  FTensor::Tensor1<double, 3> tieDirection;
  boost::shared_ptr<MatrixDouble> lambdaPtr;
  boost::shared_ptr<VectorDouble> translationPtr;
  boost::shared_ptr<VectorDouble> thetaPtr;
  boost::shared_ptr<MatrixDouble> intRotationPtr;
};

template <int DIM>
struct OpTieTermConstrainRigidBodyGlobalRotationLhs_dLambda
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstrainRigidBodyGlobalRotationLhs_dLambda(
      std::string lambda_name, std::string field_name,
      boost::shared_ptr<MatrixDouble> u_ptr,
      FTensor::Tensor1<double, 3> tie_coord,
      FTensor::Tensor1<double, 3> tie_direction,
      boost::shared_ptr<Range> tie_faces_ptr,
      boost::shared_ptr<MatrixDouble> lambda_ptr,
      boost::shared_ptr<VectorDouble> translation_ptr,
      boost::shared_ptr<VectorDouble> theta_ptr)
      : OpBase(lambda_name, field_name, OpBase::OPROWCOL, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction),
        lambdaPtr(lambda_ptr), translationPtr(translation_ptr),
        thetaPtr(theta_ptr) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    // doEntities[MBEDGE] = true;
    // doEntities[MBTRI] = true;
    // doEntities[MBQUAD] = true;
    // doEntities[MBENTITYSET] = true;
    this->assembleTranspose = true;
    this->sYmm = false;
  }

  // MoFEMErrorCode doWork(int side, EntityType type,
  //                       EntitiesFieldData::EntData &data) {
  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);
    FTENSOR_INDEX(SPACE_DIM, m);

    // double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();
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
    // get lambda
    auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

    FTensor::Tensor1<double, 3> t_translation;
    t_translation(0) = (*translationPtr)(0);
    t_translation(1) = (*translationPtr)(1);
    t_translation(2) = (*translationPtr)(2);

    FTensor::Tensor1<double, 3> t_theta;
    t_theta(0) = (*thetaPtr)(0);
    t_theta(1) = (*thetaPtr)(1);
    t_theta(2) = (*thetaPtr)(2);

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

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
      t_tangent(i, k) = alpha * levi_civita(i,j,k) * (t_coords(j) - tieCoord(j));

      ++ t_lambda;
      ++t_u;
      ++t_coords;

      int rr = 0;
      for (; rr != 3; ++rr) {
        auto t_mat = getFTensor1FromPtr<SPACE_DIM>(&OpBase::locMat(rr, 0));
        auto t_col_base = data.getFTensor0N(gg, 0);
        for (int cc = 0; cc != OpBase::nbCols/SPACE_DIM; cc++) {
          t_mat(j) -= t_col_base * t_dummy_base(i) * t_tangent(i, j);
          ++t_col_base;
          ++t_mat;
        }
        ++t_dummy_base;
      }
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> uPtr;
  FTensor::Tensor1<double, 3> tieCoord;
  FTensor::Tensor1<double, 3> tieDirection;
  boost::shared_ptr<MatrixDouble> lambdaPtr;
  boost::shared_ptr<VectorDouble> translationPtr;
  boost::shared_ptr<VectorDouble> thetaPtr;
};

// template <int Tensor_Dim>
struct OpCalculateNoFieldVectorValues : public ForcesAndSourcesCore::UserDataOperator{
  
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

      //auto values_at_gauss_pts = getFTensor0FromVec(vec);
      
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

struct OpTieSetDofs : public ForcesAndSourcesCore::UserDataOperator {

  OpTieSetDofs(const std::string field_name,
               boost::shared_ptr<VectorDouble> data_ptr,
               FTensor::Tensor1<double, 3> tie_direction,
               const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), tieDirection(tie_direction), zeroType(zero_type) {
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
    }

    VectorDouble &FieldData = data.getFieldData();

    FieldData(0) = tieDirection(0);
    std::cout << "FieldData(0): " << FieldData(0) << std::endl;

    //vec(0) = tieDirection(0);
    // std::cout << "vec(0): " << vec(0) << std::endl;
    // std::cout << "vec(1): " << vec(1) << std::endl;
    // std::cout << "vec(2): " << vec(2) << std::endl;

    MoFEMFunctionReturn(0);
  }

protected:
  boost::shared_ptr<ublas::vector<double, DoubleAllocator>> dataPtr;
  const EntityHandle zeroType;
  FTensor::Tensor1<double, 3> tieDirection;
};

template <int DIM>
struct OpTieTermSetBc : public ForcesAndSourcesCore::UserDataOperator {

  OpTieTermSetBc(std::string lambda_name, std::string field_name,
                 boost::shared_ptr<Range> tie_faces_ptr)
      : ForcesAndSourcesCore::UserDataOperator(
            lambda_name, field_name,
            ForcesAndSourcesCore::UserDataOperator::OPROWCOL) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    // doEntities[MBEDGE] = true;
    // doEntities[MBTRI] = true;
    // doEntities[MBQUAD] = true;
    doEntities[MBENTITYSET] = true;
    this->sYmm = false;
  }

  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        EntitiesFieldData::EntData &row_data,
                        EntitiesFieldData::EntData &col_data) {
    // MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);
    FTENSOR_INDEX(SPACE_DIM, l);
    FTENSOR_INDEX(SPACE_DIM, m);

    // double time = getTStime();
    double time = 1.0;

    auto nb_integration_pts = getGaussPts().size2();

    VectorDouble &vec = row_data.getFieldData();


      FTensor::Tensor1<double, 3> t_set_bc;
      t_set_bc(0) = 1.0;
      t_set_bc(1) = 0.0;
      t_set_bc(2) = 0.0;

      vec(0) = t_set_bc(0);
      vec(1) = t_set_bc(1);
      vec(2) = t_set_bc(2);



    MoFEMFunctionReturn(0);
  }
};

struct SetTieBcPreProc{
  SetTieBcPreProc(MoFEM::Interface &m_field,
                   boost::shared_ptr<FEMethod> fe_ptr,
                  FTensor::Tensor1<double, 3> tie_direction);

  virtual ~SetTieBcPreProc() = default;

  MoFEMErrorCode operator()();

protected:
  MoFEM::Interface &mField;
  boost::weak_ptr<FEMethod> fePtr;
  FTensor::Tensor1<double, 3> tieDirection;
};

SetTieBcPreProc::SetTieBcPreProc(MoFEM::Interface &m_field,
                                 boost::shared_ptr<FEMethod> fe_ptr,
                                 FTensor::Tensor1<double, 3> tie_direction)
    : mField(m_field), fePtr(fe_ptr), tieDirection(tie_direction) {}

MoFEMErrorCode SetTieBcPreProc::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  if (auto fe_method_ptr = fePtr.lock()) {

    auto fb = mField.getInterface<FieldBlas>();
    const auto problem_name = fe_method_ptr->problemPtr->getName();
    const auto field_name = "RIGID_BODY_LAMBDA";

    auto get_field_coeffs = [&](auto field_name) {
      auto field_ptr = mField.get_field_structure(field_name);
      return field_ptr->getNbOfCoeffs();
    };
    const auto nb_field_coeffs = get_field_coeffs(field_name);

    MOFEM_LOG("WORLD", Sev::noisy) << "Apply Tie Translation: " << problem_name
                                   << "_" << field_name << std::endl;

    FTensor::Tensor1<double, 3> t_vals{0., 0., 0.};

    // auto scale_value = [&](const double &c) {
    //   double val = c;
    //   for (auto s : vecOfTimeScalingMethods) {
    //     val *= s->getScale(fe_method_ptr->ts_t);
    //   }
    //   return val;
    // };

    //t_vals(0) = scale_value(tieDirection(0));
    t_vals(0) = tieDirection(0);
    int idx = 0;
    int coeff = 0;

    auto lambda = [&](boost::shared_ptr<FieldEntity> field_entity_ptr) {
      MoFEMFunctionBegin;

      auto v = t_vals(coeff);

      field_entity_ptr->getEntFieldData()[coeff] = v;
      //++idx;

      MoFEMFunctionReturn(0);
    };

    // auto zero_lambda = [&](boost::shared_ptr<FieldEntity> field_entity_ptr) {
    //   MoFEMFunctionBegin;
    //   auto size = field_entity_ptr->getEntFieldData().size();
    //   for (int i = coeff; i < size; i += nb_field_coeffs)
    //     field_entity_ptr->getEntFieldData()[i] = 0;
    //   MoFEMFunctionReturn(0);
    // };

    //auto verts = bc.second->bcEnts.subset_by_type(MBVERTEX);
    auto meshset_set = mField.get_field_meshset(field_name);
    Range tie_ents;
    CHKERR mField.get_moab().get_entities_by_handle(meshset_set, tie_ents, true);
    //auto not_verts = subtract(bc.second->bcEnts, verts);

    // idx = 0;
    // coeff = 0;
    CHKERR fb->fieldLambdaOnEntities(lambda, field_name, &tie_ents);
    //CHKERR fb->fieldLambdaOnEntities(zero_lambda, field_name, &not_verts);

  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Can not lock shared pointer");
  }

  MoFEMFunctionReturn(0);
}
