/**
 * @file TieConstraint.hpp
 * @brief Tie Constraint Implementation
 * @date 2024-11-14
 *
 * @copyright Copyright (c) 2024
 *
 */

template <int DIM>
struct OpTieTermConstraintRhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstraintRhs(std::string lambda_name,
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
          (t_coords(i) + t_u(i)) - (tieCoord(i) + tieDirection(i));
      ++t_u;
      FTensor::Tensor1<double, SPACE_DIM> t_delta_initial;
      t_delta_initial(i) = t_coords(i) - tieCoord(i);
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
struct OpTieTermConstraintLhs
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {
  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstraintLhs(std::string lambda_name, std::string col_field_name,
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
          (t_coords(i) + t_u(i)) - (tieCoord(i) + tieDirection(i));
      ++t_u;
      ++t_coords;
      FTensor::Tensor1<double, SPACE_DIM> t_tangent;
      t_tangent(i) = alpha * (t_delta_current(i) / t_delta_current.l2());
      int rr = 0;
      for (; rr != OpBase::nbRows; ++rr) {
        auto t_col_base = col_data.getFTensor0N(gg, 0);
        auto t_mat = getFTensor1FromPtr<SPACE_DIM>(&OpBase::locMat(rr, 0));
        for (int cc = 0; cc != OpBase::nbCols/SPACE_DIM; cc++) {
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
};


template <int DIM> 
struct OpTieTermConstraintRhs_du
    : public FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase {

  using OpBase = FormsIntegrators<BoundaryEleOp>::Assembly<A>::OpBase;

  OpTieTermConstraintRhs_du(std::string field_name, std::string lambda_name,
                            boost::shared_ptr<MatrixDouble> u_ptr,
                            boost::shared_ptr<VectorDouble> lambda_ptr,
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
      t_du(i) = t_u(i) + (t_coords(i) - tieCoord(i) - tieDirection(i));
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