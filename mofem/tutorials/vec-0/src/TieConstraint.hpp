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

    double time = getTStime();
    //double time = 1.0;

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

    double time = getTStime();
    //double time = 1.0;

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

    double time = getTStime();
    //double time = 1.0;

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

    double time = getTStime();
    //double time = 1.0;

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
                                    boost::shared_ptr<VectorDouble> u_ref_ptr,
                                    FTensor::Tensor1<double, 3> tie_coord,
                                    FTensor::Tensor1<double, 3> tie_direction,
                                    boost::shared_ptr<Range> tie_faces_ptr,
                                    boost::shared_ptr<MatrixDouble> vel_ptr,
                                    boost::shared_ptr<VectorDouble> vel_ref_ptr)
      : OpBase(lambda_name, lambda_name, OpBase::OPROW, tie_faces_ptr),
        uPtr(u_ptr), tieCoord(tie_coord), tieDirection(tie_direction),
        velPtr(vel_ptr), velRefPtr(vel_ref_ptr), uRefPtr(u_ref_ptr) {
    std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
    //doEntities[MBEDGE] = true;
    doEntities[MBTRI] = true;
    doEntities[MBQUAD] = true;
  }

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;
    FTENSOR_INDEX(SPACE_DIM, i);
    FTENSOR_INDEX(SPACE_DIM, j);
    FTENSOR_INDEX(SPACE_DIM, k);

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
    // get velocity
    auto t_vel = getFTensor1FromMat<SPACE_DIM>(*velPtr);  


    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor1<double, 3> t_delta_initial;
      t_delta_initial(i) = (t_coords(i) - tieCoord(i));

      FTensor::Tensor1<double, 3> t_disp_ref;
      t_disp_ref(1) = (*uRefPtr)(1);
      t_disp_ref(2) = (*uRefPtr)(2);
      t_disp_ref(0) = tieDirection(0) * time;

      FTensor::Tensor1<double, 3> t_delta_disp;
      t_delta_disp(i) = t_u(i) - t_disp_ref(i);



      FTensor::Tensor1<double, 3> t_delta_vel;
      t_delta_vel(0) = t_vel(0) - (*velRefPtr)(0);
      t_delta_vel(1) = t_vel(1) - (*velRefPtr)(1);
      t_delta_vel(2) = t_vel(2) - (*velRefPtr)(2);

      FTensor::Tensor1<double, 3> t_theta;
      t_theta(i) = t_delta_vel(i) / t_delta_initial.l2();


      ++t_coords;
      ++t_u;

      FTensor::Tensor1<double, SPACE_DIM> g;
      g(i) = alpha * (t_delta_disp(i) + levi_civita(i, j, k) * t_theta(j) * t_delta_initial(k));
      

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
  boost::shared_ptr<MatrixDouble> velPtr;
  boost::shared_ptr<VectorDouble> velRefPtr;
  boost::shared_ptr<VectorDouble> uRefPtr;
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
    //doEntities[MBEDGE] = true;
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

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_tangent;
      t_tangent(i, j) = alpha * (t_kd(i, j));


      int rr = 0;
      for (; rr != OpBase::nbRows/SPACE_DIM; ++rr) {
        auto t_col_base = col_data.getFTensor0N(gg, 0);
        auto t_mat = getFTensor2FromArray<DIM, DIM, DIM>(OpBase::locMat, DIM * rr);
        for (int cc = 0; cc != OpBase::nbCols/SPACE_DIM; cc++) {
          t_mat(0, 0) += (t_row_base * t_col_base) * t_tangent(0, 0);
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
    //doEntities[MBEDGE] = true;
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
    auto t_lambda = getFTensor1FromMat<SPACE_DIM>(*lambdaPtr);

    auto &nf = OpBase::locF;

    for (int gg = 0; gg != nb_integration_pts; gg++) {
      const auto alpha = t_w * vol;
      ++t_w;

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_du;
      t_du(i,j) = (t_kd(i,j));

      //std::cout << "t_lambda: " << t_lambda << std::endl;
  
      auto t_nf = getFTensor1FromPtr<DIM>(&nf[0]);
      FTensor::Tensor1<double, SPACE_DIM> g;
      g(j) = alpha * t_lambda(i) * (t_du(i, j));
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