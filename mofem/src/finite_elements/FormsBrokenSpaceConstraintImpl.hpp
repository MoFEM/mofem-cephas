/**
 * @file FormsBrokenSpaceConstraintImpl.hpp
 * @brief Integrator for broken space constraints
 * @date 2024-07-01
 *
 * @copyright Copyright (c) 2024
 *
 */

namespace MoFEM {

#ifdef BILINEAR_FORMS_INTEGRATORS_HPP

template <int FIELD_DIM, IntegrationType I, typename OpBase>
struct OpBrokenSpaceConstrainImpl;

template <int FIELD_DIM, typename OpBase>
struct OpBrokenSpaceConstrainImpl<FIELD_DIM, GAUSS, OpBase> : public OpBase {

  using OP = OpBase;

  OpBrokenSpaceConstrainImpl(
      const std::string row_field,
      boost::shared_ptr<BrokenBaseSideData> broken_base_side_data,
      const double beta, const bool assmb_transpose, const bool only_transpose)
      : OpBase(row_field, row_field, OpBase::OPROW),
        brokenBaseSideData(broken_base_side_data), scalarBeta(beta) {
    this->assembleTranspose = assmb_transpose;
    this->onlyTranspose = only_transpose;
  }

  OpBrokenSpaceConstrainImpl(const std::string row_field,
                             const std::string col_field, const double beta,
                             const bool assmb_transpose,
                             const bool only_transpose)
      : OpBase(row_field, col_field, OpBase::OPROWCOL),
        brokenBaseSideData(nullptr), scalarBeta(beta) {
    this->assembleTranspose = assmb_transpose;
    this->onlyTranspose = only_transpose;
  }

  MoFEMErrorCode doWork(int row_side, EntityType row_type,
                        EntitiesFieldData::EntData &row_data);

protected:
  double scalarBeta;
  boost::shared_ptr<BrokenBaseSideData> brokenBaseSideData;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);

  MatrixDouble transposedC;
};

template <int FIELD_DIM, typename OpBase>
MoFEMErrorCode OpBrokenSpaceConstrainImpl<FIELD_DIM, GAUSS, OpBase>::doWork(
    int row_side, EntityType row_type, EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  if (OP::entsPtr) {
    if (OP::entsPtr->find(this->getFEEntityHandle()) == OP::entsPtr->end())
      MoFEMFunctionReturnHot(0);
  }

#ifndef NDEBUG
  if (!brokenBaseSideData) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "space not set");
  }
  if (brokenBaseSideData->getData().getSpace() != HDIV &&
      brokenBaseSideData->getData().getSpace() != HCURL) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            (std::string("Expect Hdiv or Hcurl space but received ") +
             FieldSpaceNames[brokenBaseSideData->getData().getSpace()])
                .c_str());
  }
#endif // NDEBUG

  OP::nbRows = row_data.getIndices().size();

  // if no dofs on row, exit that work, nothing to do here

  if (!OP::nbRows)
    MoFEMFunctionReturnHot(0);

  OP::rowSide = row_side;
  OP::rowType = row_type;

  // get number of dofs on column
  auto &col_data = brokenBaseSideData->getData();
  OP::nbCols = col_data.getIndices().size();
  // if no dofs on column, exit nothing to do here
  if (!OP::nbCols)
    MoFEMFunctionReturnHot(0);

  OP::colSide = brokenBaseSideData->getSide();
  OP::colType = brokenBaseSideData->getType();

  // get number of integration points
  OP::nbIntegrationPts = OP::getGaussPts().size2();
  // get row base functions
  OP::nbRowBaseFunctions = OP::getNbOfBaseFunctions(row_data);

  // set size of local entity bock
  OP::locMat.resize(OP::nbRows, OP::nbCols, false);
  // clear matrix
  OP::locMat.clear();
  // integrate local matrix for entity block
  CHKERR this->iNtegrate(row_data, col_data);

  // assemble local matrix
  auto check_if_assemble_transpose = [&] {
    if (this->sYmm) {
      if (OP::rowSide != OP::colSide || OP::rowType != OP::colType)
        return true;
      else
        return false;
    } else if (OP::assembleTranspose) {
      return true;
    }

    return false;
  };
  CHKERR this->aSsemble(row_data, col_data, check_if_assemble_transpose());

  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, typename OpBase>
MoFEMErrorCode OpBrokenSpaceConstrainImpl<FIELD_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data,
    EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  auto nb_row_dofs = row_data.getIndices().size();
  auto nb_col_dofs = col_data.getIndices().size();
  if (!nb_row_dofs || !nb_col_dofs)
    MoFEMFunctionReturnHot(0);

#ifndef NDEBUG
  if (brokenBaseSideData->getData().getSpace() != HDIV &&
      brokenBaseSideData->getData().getSpace() == HCURL) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Expect Hdiv or Hcurl space");
  }
#endif

  FTENSOR_INDEX(FIELD_DIM, i);
  FTENSOR_INDEX(FIELD_DIM, j);

  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();
  size_t nb_base_functions = col_data.getN().size2() / 3;

  auto t_col_base = col_data.getFTensor1N<3>();

  OP::locMat.resize(nb_row_dofs, nb_col_dofs, false);

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    auto t_m = getFTensor1FromPtr<FIELD_DIM>(&*transposedC.data().begin());

    int cc = 0;
    for (; cc != nb_col_dofs / FIELD_DIM; cc++) {
      auto t_row_base = row_data.getFTensor0N(gg, 0);
      for (auto rr = 0; rr != nb_row_dofs / FIELD_DIM; ++rr) {
        t_m(i) += t_w * t_row_base * t_normal(j) * t_col_base(j);
        ++t_row_base;
        ++t_m;
      }

      ++t_col_base;
    }

    for (; cc < nb_base_functions; ++cc)
      ++t_col_base;

    ++t_w;
    ++t_normal;
  }

  transposedC *= scalarBeta;
  noalias(OP::locMat) = trans(transposedC);

  MoFEMFunctionReturn(0);
}

#endif // BILINEAR_FORMS_INTEGRATORS_HPP

#ifdef LINER_FORMS_INTEGRATORS_HPP

struct BrokenBaseSideData {
  BrokenBaseSideData() {}
  inline auto &getSide() { return eleSide; }
  inline auto &getType() { return eleType; }
  inline auto &getData() { return entData; }

private:
  int eleSide = 1;
  EntityType eleType = MBENTITYSET;
  EntitiesFieldData::EntData entData;
};

struct OpGetBrokenBaseSideData : public ForcesAndSourcesCore::UserDataOperator {

  using OP = ForcesAndSourcesCore::UserDataOperator;

  OpGetBrokenBaseSideData(
      const std::string field_name,
      boost::shared_ptr<BrokenBaseSideData> broken_base_side_data);
  MoFEMErrorCode doWork(int row_side, EntityType row_type,
                        EntitiesFieldData::EntData &row_data);

private:
  boost::shared_ptr<BrokenBaseSideData> brokenBaseSideData;
};

template <int FIELD_DIM, IntegrationType I, typename OpBase>
struct OpBrokenSpaceConstrainDHybridImpl;

template <int FIELD_DIM, IntegrationType I, typename OpBase>
struct OpBrokenSpaceConstrainDFluxImpl;

template <int FIELD_DIM, typename OpBase>
struct OpBrokenSpaceConstrainDFluxImpl<FIELD_DIM, GAUSS, OpBase>
    : public OpBase {

  using OP = OpBase;

  OpBrokenSpaceConstrainDFluxImpl(
      boost::shared_ptr<BrokenBaseSideData> broken_base_side_data,
      boost::shared_ptr<MatrixDouble> lagrange_ptr, const double beta)
      : OpBase(NOSPACE, OpBase::OPSPACE), scalarBeta(beta),
        brokenBaseSideData(broken_base_side_data), lagrangePtr(lagrange_ptr) {}

  OpBrokenSpaceConstrainDFluxImpl(std::string row_field,
                                  boost::shared_ptr<MatrixDouble> lagrange_ptr,
                                  const double beta)
      : OpBase(row_field, row_field, OpBase::OPROW), scalarBeta(beta),
        brokenBaseSideData(nullptr), lagrangePtr(lagrange_ptr) {}

  MoFEMErrorCode doWork(int row_side, EntityType row_type,
                        EntitiesFieldData::EntData &row_data);

private:
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data);

  double scalarBeta;
  boost::shared_ptr<BrokenBaseSideData> brokenBaseSideData;
  boost::shared_ptr<MatrixDouble> lagrangePtr;
};

template <int FIELD_DIM, typename OpBase>
MoFEMErrorCode
OpBrokenSpaceConstrainDFluxImpl<FIELD_DIM, GAUSS, OpBase>::doWork(
    int side, EntityType type, EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  if (OP::entsPtr) {
    if (OP::entsPtr->find(this->getFEEntityHandle()) == OP::entsPtr->end())
      MoFEMFunctionReturnHot(0);
  }

  EntitiesFieldData::EntData *raw_data_ptr = nullptr;
  switch (OP::opType) {
  case OpBaseDerivativesBase::OPROW:
    OP::rowSide = side;
    OP::rowType = type;
    raw_data_ptr = &data;
    break;
  case OpBaseDerivativesBase::OPSPACE:
    OP::rowSide = brokenBaseSideData->getSide();
    OP::rowType = brokenBaseSideData->getType();
    raw_data_ptr = &brokenBaseSideData->getData();
    break;
  default:
    CHK_MOAB_THROW(MOFEM_IMPOSSIBLE_CASE,
                   (std::string("wrong op type ") +
                    OpBaseDerivativesBase::OpTypeNames[OP::opType])
                       .c_str());
  }

  // get number of dofs on row
  OP::nbRows = raw_data_ptr->getIndices().size();

  if (!OP::nbRows)
    MoFEMFunctionReturnHot(0);
  // get number of integration points
  OP::nbIntegrationPts = OP::getGaussPts().size2();
  // get row base functions
  OP::nbRowBaseFunctions = OP::getNbOfBaseFunctions(*raw_data_ptr);
  // resize and clear the right hand side vector
  OP::locF.resize(OP::nbRows);
  OP::locF.clear();
  // integrate local vector
  CHKERR this->iNtegrate(*raw_data_ptr);
  // assemble local vector
  CHKERR this->aSsemble(*raw_data_ptr);

  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, typename OpBase>
MoFEMErrorCode
OpBrokenSpaceConstrainDFluxImpl<FIELD_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  FTENSOR_INDEX(FIELD_DIM, i);
  FTENSOR_INDEX(FIELD_DIM, j);

  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();
  auto t_lagrange = getFTensor1FromMat<FIELD_DIM>(*lagrangePtr);

  OP::locF.resize(row_data.getIndices().size(), false);
  OP::locF.clear();

  auto t_row_base = row_data.getFTensor1N<3>();
  auto nb_base_functions = row_data.getN().size2() / 3;

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {
    auto t_vec = getFTensor1FromPtr<FIELD_DIM>(&*OP::locF.data().begin());
    size_t rr = 0;
    for (; rr != row_data.getIndices().size() / FIELD_DIM; ++rr) {
      t_vec(i) += t_w * t_row_base(j) * t_normal(j) * t_lagrange(i);
      ++t_row_base;
      ++t_vec;
    }
    for (; rr < nb_base_functions; ++rr)
      ++t_row_base;
    ++t_w;
    ++t_normal;
    ++t_lagrange;
  }

  OP::locF *= scalarBeta;

  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, typename OpBase>
struct OpBrokenSpaceConstrainDHybridImpl<FIELD_DIM, GAUSS, OpBase>
    : public OpBase {

  using OP = OpBase;

  OpBrokenSpaceConstrainDHybridImpl(const std::string row_field,
                                    boost::shared_ptr<MatrixDouble> flux_ptr,
                                    const double beta)
      : OpBase(row_field, row_field, OpBase::OPROW), fluxPtr(flux_ptr),
        scalarBeta(beta) {}

private:
  double scalarBeta;
  boost::shared_ptr<MatrixDouble> fluxPtr;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data);
};

template <int FIELD_DIM, typename OpBase>
MoFEMErrorCode
OpBrokenSpaceConstrainDHybridImpl<FIELD_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  FTENSOR_INDEX(FIELD_DIM, i);
  FTENSOR_INDEX(FIELD_DIM, j);

  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();

  OP::locF.resize(row_data.getIndices().size(), false);
  OP::locF.clear();

  auto t_row_base = row_data.getFTensor0N();
  auto t_flux = getFTensor1FromMat<FIELD_DIM>(*fluxPtr);
  auto nb_base_functions = row_data.getN().size2() / 3;

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {
    auto t_vec = getFTensor1FromPtr<FIELD_DIM>(&*OP::locF.data().begin());
    size_t rr = 0;
    for (; rr != row_data.getIndices().size() / FIELD_DIM; ++rr) {
      t_flux(i) += t_w * t_row_base * t_normal(j) * t_flux(j);
      ++t_row_base;
      ++t_vec;
    }
    for (; rr < nb_base_functions; ++rr)
      ++t_row_base;
    ++t_w;
    ++t_normal;
    ++t_flux;
  }

  OP::locF *= scalarBeta;

  MoFEMFunctionReturn(0);
}

#endif // LINER_FORMS_INTEGRATORS_HPP

} // namespace MoFEM