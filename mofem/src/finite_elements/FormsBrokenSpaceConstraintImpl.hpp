/**
 * @file FormsBrokenSpaceConstraintImpl.hpp
 * @brief Integrator for broken space constraints
 * @date 2024-07-01
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef __FORMSBROKENSPACECONSTRAINTIMPL_HPP__
#define __FORMSBROKENSPACECONSTRAINTIMPL_HPP__

namespace MoFEM {

template <typename E> struct OpBrokenLoopSide : public OpLoopSide<E> {

  using OP = OpLoopSide<E>;
  using OpLoopSide<E>::OpLoopSide;

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;

    auto prev_side_fe_ptr = OP::getSidePtrFE();
    if (OP::sideFEName == prev_side_fe_ptr->getFEName()) {
      auto prev_fe_uid =
          prev_side_fe_ptr->numeredEntFiniteElementPtr->getFEUId();
      OP::sideFEPtr->feName = OP::sideFEName;
      CHKERR OP::sideFEPtr->setSideFEPtr(OP::getPtrFE());
      CHKERR OP::sideFEPtr->copyBasicMethod(*OP::getFEMethod());
      CHKERR OP::sideFEPtr->copyPetscData(*OP::getFEMethod());
      CHKERR OP::sideFEPtr->copyKsp(*OP::getFEMethod());
      CHKERR OP::sideFEPtr->copySnes(*OP::getFEMethod());
      CHKERR OP::sideFEPtr->copyTs(*OP::getFEMethod());
      OP::sideFEPtr->cacheWeakPtr = prev_side_fe_ptr->cacheWeakPtr;
      OP::sideFEPtr->loopSize = 1;
      CHKERR OP::sideFEPtr->preProcess();
      OP::sideFEPtr->nInTheLoop = 0;
      OP::sideFEPtr->numeredEntFiniteElementPtr =
          prev_side_fe_ptr->numeredEntFiniteElementPtr;
      CHKERR (*OP::sideFEPtr)();
      CHKERR OP::sideFEPtr->postProcess();
    } else {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "sideFEName is different");
    }

    MoFEMFunctionReturn(0);
  };
};

struct BrokenBaseSideData {
  BrokenBaseSideData() : entData(false) {}
  inline auto &getSense() { return eleSense; }
  inline auto &getSide() { return eleSide; }
  inline auto &getType() { return eleType; }
  inline auto &getData() { return entData; }
  inline auto &getFlux() { return fluxMat; }

private:
  int eleSense = 0;
  int eleSide = 1;
  EntityType eleType = MBENTITYSET;
  EntitiesFieldData::EntData entData;
  MatrixDouble fluxMat;
};

template <typename OpBase> struct OpGetBrokenBaseSideData : public OpBase {

  using OP = OpBase;

  OpGetBrokenBaseSideData(
      const std::string field_name,
      boost::shared_ptr<std::vector<BrokenBaseSideData>> broken_base_side_data);

  MoFEMErrorCode doWork(int row_side, EntityType row_type,
                        EntitiesFieldData::EntData &row_data);

private:
  boost::shared_ptr<std::vector<BrokenBaseSideData>> brokenBaseSideData;
};

template <typename OpBase>
OpGetBrokenBaseSideData<OpBase>::OpGetBrokenBaseSideData(
    const std::string field_name,
    boost::shared_ptr<std::vector<BrokenBaseSideData>> broken_base_side_data)
    : OP(field_name, field_name, OP::OPROW),
      brokenBaseSideData(broken_base_side_data) {}

template <typename OpBase>
MoFEMErrorCode
OpGetBrokenBaseSideData<OpBase>::doWork(int row_side, EntityType row_type,
                                        EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  if (row_data.getIndices().size() == 0)
    MoFEMFunctionReturnHot(0);

  brokenBaseSideData->resize(OP::getLoopSize());

  const auto n_in_the_loop = OP::getNinTheLoop();
  const auto face_sense = OP::getSkeletonSense();

#ifndef NDEBUG
  if (face_sense != -1 && face_sense != 1)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "face sense not set");
#endif // NDEBUG

  auto set_data = [&](auto &side_data) {
    side_data.getSide() = row_side;
    side_data.getType() = row_type;
    side_data.getSense() = face_sense;
    side_data.getData().sEnse = row_data.sEnse;
    side_data.getData().sPace = row_data.sPace;
    side_data.getData().bAse = row_data.bAse;
    side_data.getData().iNdices = row_data.iNdices;
    side_data.getData().localIndices = row_data.localIndices;
    side_data.getData().dOfs = row_data.dOfs;
    side_data.getData().fieldEntities = row_data.fieldEntities;
    side_data.getData().fieldData = row_data.fieldData;
  };

  auto set_base = [&](auto &side_data) {
    auto base = side_data.getData().getBase();
    for (auto dd = 0; dd != BaseDerivatives::LastDerivative; ++dd) {
      if (auto base_ptr = row_data.baseFunctionsAndBaseDerivatives[dd][base]) {
        side_data.getData().baseFunctionsAndBaseDerivatives[dd][base] =
            boost::make_shared<MatrixDouble>(*base_ptr);
      }
    }
  };

  set_data((*brokenBaseSideData)[n_in_the_loop]);
  set_base((*brokenBaseSideData)[n_in_the_loop]);

  MoFEMFunctionReturn(0);
}

template <typename OpBase> struct OpSetFlux : public OpBase {

  using OP = OpBase;

  OpSetFlux(
      boost::shared_ptr<std::vector<BrokenBaseSideData>> broken_base_side_data,
      boost::shared_ptr<MatrixDouble> flux_ptr);

  MoFEMErrorCode doWork(int row_side, EntityType row_type,
                        EntitiesFieldData::EntData &row_data);

private:
  boost::shared_ptr<std::vector<BrokenBaseSideData>> brokenBaseSideData;
  boost::shared_ptr<MatrixDouble> fluxPtr;
};

template <typename OpBase>
OpSetFlux<OpBase>::OpSetFlux(
    boost::shared_ptr<std::vector<BrokenBaseSideData>> broken_base_side_data,
    boost::shared_ptr<MatrixDouble> flux_ptr)
    : OP(NOSPACE, OP::OPSPACE), brokenBaseSideData(broken_base_side_data),
      fluxPtr(flux_ptr) {}

template <typename OpBase>
MoFEMErrorCode OpSetFlux<OpBase>::doWork(int side, EntityType type,
                                         EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
  auto swap_flux = [&](auto &side_data) { side_data.getFlux().swap(*fluxPtr); };
  swap_flux((*brokenBaseSideData)[OP::getNinTheLoop()]);
  MoFEMFunctionReturn(0);
}

template <typename OpBase> struct OpBrokenBaseImpl : public OpBase {

  using OP = OpBase;

  OpBrokenBaseImpl(
      boost::shared_ptr<std::vector<BrokenBaseSideData>> broken_base_side_data,
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OP(NOSPACE, OP::OPSPACE), brokenBaseSideData(broken_base_side_data) {
    OP::entsPtr = ents_ptr;
  }

  OpBrokenBaseImpl(
      const std::string row_field,
      boost::shared_ptr<std::vector<BrokenBaseSideData>> broken_base_side_data,
      const bool assmb_transpose, const bool only_transpose,
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OP(row_field, row_field, OP::OPROW, ents_ptr),
        brokenBaseSideData(broken_base_side_data) {
    OP::entsPtr = ents_ptr;
    OP::assembleTranspose = assmb_transpose;
    OP::onlyTranspose = only_transpose;
    OP::sYmm = false;
  }

  MoFEMErrorCode doWork(int row_side, EntityType row_type,
                        EntitiesFieldData::EntData &row_data);

protected:
  boost::shared_ptr<std::vector<BrokenBaseSideData>> brokenBaseSideData;
};

template <typename OpBase>
MoFEMErrorCode
OpBrokenBaseImpl<OpBase>::doWork(int row_side, EntityType row_type,
                                 EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  if (OP::entsPtr) {
    if (OP::entsPtr->find(this->getFEEntityHandle()) == OP::entsPtr->end())
      MoFEMFunctionReturnHot(0);
  }

#ifndef NDEBUG
  if (!brokenBaseSideData) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "space not set");
  }
#endif // NDEBUG

  auto do_work_rhs = [this](int row_side, EntityType row_type,
                            EntitiesFieldData::EntData &row_data, int sense) {
    MoFEMFunctionBegin;
    // get number of dofs on row
    OP::nbRows = row_data.getIndices().size();
    if (!OP::nbRows)
      MoFEMFunctionReturnHot(0);
    // get number of integration points
    OP::nbIntegrationPts = OP::getGaussPts().size2();
    // get row base functions
    OP::nbRowBaseFunctions = OP::getNbOfBaseFunctions(row_data);
    // resize and clear the right hand side vector
    OP::locF.resize(OP::nbRows, false);
    OP::locF.clear();
    // integrate local vector
    CHKERR this->iNtegrate(row_data);
    // assemble local vector
    OP::locF *= sense;
    CHKERR this->aSsemble(row_data);
    MoFEMFunctionReturn(0);
  };

  auto do_work_lhs = [this](int row_side, int col_side, EntityType row_type,
                            EntityType col_type,
                            EntitiesFieldData::EntData &row_data,
                            EntitiesFieldData::EntData &col_data, int sense) {
    MoFEMFunctionBegin;

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

    OP::rowSide = row_side;
    OP::rowType = row_type;
    OP::colSide = col_side;
    OP::colType = col_type;
    OP::nbCols = col_data.getIndices().size();
    OP::locMat.resize(OP::nbRows, OP::nbCols, false);
    OP::locMat.clear();
    CHKERR this->iNtegrate(row_data, col_data);
    OP::locMat *= sense;
    CHKERR this->aSsemble(row_data, col_data, check_if_assemble_transpose());
    MoFEMFunctionReturn(0);
  };

  switch (OP::opType) {
  case OP::OPROW:

    OP::nbRows = row_data.getIndices().size();
    if (!OP::nbRows)
      MoFEMFunctionReturnHot(0);
    OP::nbIntegrationPts = OP::getGaussPts().size2();
    OP::nbRowBaseFunctions = OP::getNbOfBaseFunctions(row_data);

    if (!OP::nbRows)
      MoFEMFunctionReturnHot(0);

    for (auto &bd : *brokenBaseSideData) {

#ifndef NDEBUG
      if (bd.getData().getSpace() != HDIV && bd.getData().getSpace() != HCURL) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                (std::string("Expect Hdiv or Hcurl space but received ") +
                 FieldSpaceNames[bd.getData().getSpace()])
                    .c_str());
      }
      if (!bd.getData().getNSharedPtr(bd.getData().getBase())) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "base functions not set");
      }
#endif

      CHKERR do_work_lhs(

          // side
          row_side, bd.getSide(),

          // type
          row_type, bd.getType(),

          // row_data
          row_data, bd.getData(),

          // sense
          bd.getSense()

      );
    }

    break;
  case OP::OPSPACE:
    for (auto &bd : *brokenBaseSideData) {
      CHKERR do_work_rhs(bd.getSide(), bd.getType(), bd.getData(),
                         bd.getSense());
    }
    break;
  default:
    CHK_MOAB_THROW(MOFEM_IMPOSSIBLE_CASE,
                   (std::string("wrong op type ") +
                    OpBaseDerivativesBase::OpTypeNames[OP::opType])
                       .c_str());
  }

  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, IntegrationType I, typename OpBrokenBase>
struct OpBrokenSpaceConstrainImpl;

template <int FIELD_DIM, typename OpBrokenBase>
struct OpBrokenSpaceConstrainImpl<FIELD_DIM, GAUSS, OpBrokenBase>
    : public OpBrokenBase {

  using OP = OpBrokenBase;

  OpBrokenSpaceConstrainImpl(
      const std::string row_field,
      boost::shared_ptr<std::vector<BrokenBaseSideData>> broken_base_side_data,
      boost::shared_ptr<double> beta_ptr, const bool assmb_transpose,
      const bool only_transpose, boost::shared_ptr<Range> ents_ptr = nullptr)
      : OP(row_field, broken_base_side_data, assmb_transpose, only_transpose,
           ents_ptr),
        scalarBetaPtr(beta_ptr) {}

  OpBrokenSpaceConstrainImpl(
      const std::string row_field,
      boost::shared_ptr<std::vector<BrokenBaseSideData>> broken_base_side_data,
      double beta, const bool assmb_transpose, const bool only_transpose,
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBrokenSpaceConstrainImpl(row_field, broken_base_side_data,
                                   boost::make_shared<double>(beta),
                                   assmb_transpose, only_transpose, ents_ptr) {}

protected:
  boost::shared_ptr<double> scalarBetaPtr;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);
};

template <int FIELD_DIM, typename OpBase>
MoFEMErrorCode OpBrokenSpaceConstrainImpl<FIELD_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data,
    EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  auto nb_row_dofs = row_data.getIndices().size();
  auto nb_col_dofs = col_data.getIndices().size();
  if (!nb_row_dofs || !nb_col_dofs)
    MoFEMFunctionReturnHot(0);

  FTENSOR_INDEX(FIELD_DIM, i);
  FTENSOR_INDEX(3, J);

  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();
  size_t nb_base_functions = col_data.getN().size2() / 3;

#ifndef NDEBUG
  if (nb_row_dofs % FIELD_DIM != 0 || nb_col_dofs % FIELD_DIM != 0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "number of dofs not divisible by field dimension");
  }
  if (nb_row_dofs > row_data.getN().size2() * FIELD_DIM ||
      nb_col_dofs > col_data.getN().size2() * FIELD_DIM) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "number of dofs exceeds number of base functions");
  }
#endif // NDEBUG

  auto t_col_base = col_data.getFTensor1N<3>();

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    int cc = 0;
    for (; cc != nb_col_dofs / FIELD_DIM; cc++) {
      auto t_row_base = row_data.getFTensor0N(gg, 0);
      for (auto rr = 0; rr != nb_row_dofs / FIELD_DIM; ++rr) {
        OP::locMat(FIELD_DIM * rr, FIELD_DIM * cc) +=
            (t_w * t_row_base) * (t_normal(J) * t_col_base(J));
        ++t_row_base;
      }

      ++t_col_base;
    }

    for (; cc < nb_base_functions; ++cc)
      ++t_col_base;

    ++t_w;
    ++t_normal;
  }

  for (auto rr = 0; rr != nb_row_dofs / FIELD_DIM; ++rr) {
    for (auto cc = 0; cc != nb_col_dofs / FIELD_DIM; ++cc) {
      for (auto dd = 1; dd < FIELD_DIM; ++dd) {
        OP::locMat(FIELD_DIM * rr + dd, FIELD_DIM * cc + dd) =
            OP::locMat(FIELD_DIM * rr, FIELD_DIM * cc);
      }
    }
  }

  if (scalarBetaPtr)
    OP::locMat *= *scalarBetaPtr;

  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, IntegrationType I, typename OpBase>
struct OpBrokenSpaceConstrainDHybridImpl;

template <int FIELD_DIM, IntegrationType I, typename OpBase>
struct OpBrokenSpaceConstrainDFluxImpl;

template <int FIELD_DIM, typename OpBrokenBase>
struct OpBrokenSpaceConstrainDFluxImpl<FIELD_DIM, GAUSS, OpBrokenBase>
    : public OpBrokenBase {

  using OP = OpBrokenBase;

  OpBrokenSpaceConstrainDFluxImpl(
      boost::shared_ptr<std::vector<BrokenBaseSideData>> broken_base_side_data,
      boost::shared_ptr<MatrixDouble> lagrange_ptr,
      boost::shared_ptr<double> beta_ptr,
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OP(broken_base_side_data, ents_ptr), scalarBetaPtr(beta_ptr),
        lagrangePtr(lagrange_ptr) {}

  OpBrokenSpaceConstrainDFluxImpl(
      boost::shared_ptr<std::vector<BrokenBaseSideData>> broken_base_side_data,
      boost::shared_ptr<MatrixDouble> lagrange_ptr, double beta,
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBrokenSpaceConstrainDFluxImpl(broken_base_side_data, lagrange_ptr,
                                        boost::make_shared<double>(beta),
                                        ents_ptr) {}

private:
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data);

  boost::shared_ptr<double> scalarBetaPtr;
  boost::shared_ptr<MatrixDouble> lagrangePtr;
};

template <int FIELD_DIM, typename OpBase>
MoFEMErrorCode
OpBrokenSpaceConstrainDFluxImpl<FIELD_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  FTENSOR_INDEX(FIELD_DIM, i);
  FTENSOR_INDEX(3, J);

  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();
  auto t_lagrange = getFTensor1FromMat<FIELD_DIM>(*lagrangePtr);

  auto t_row_base = row_data.getFTensor1N<3>();
  auto nb_base_functions = row_data.getN().size2() / 3;

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {
    auto t_vec = getFTensor1FromPtr<FIELD_DIM>(&*OP::locF.data().begin());
    size_t rr = 0;
    for (; rr != row_data.getIndices().size() / FIELD_DIM; ++rr) {
      t_vec(i) += (t_w * (t_row_base(J) * t_normal(J))) * t_lagrange(i);
      ++t_row_base;
      ++t_vec;
    }
    for (; rr < nb_base_functions; ++rr)
      ++t_row_base;
    ++t_w;
    ++t_normal;
    ++t_lagrange;
  }

  if(scalarBetaPtr)
    OP::locF *= *scalarBetaPtr;

  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, typename OpBase>
struct OpBrokenSpaceConstrainDHybridImpl<FIELD_DIM, GAUSS, OpBase>
    : public OpBase {

  using OP = OpBase;

  OpBrokenSpaceConstrainDHybridImpl(
      const std::string row_field,
      boost::shared_ptr<std::vector<BrokenBaseSideData>> broken_side_data_ptr,
      boost::shared_ptr<double> beta_ptr,
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(row_field, row_field, OpBase::OPROW, ents_ptr),
        brokenSideDataPtr(broken_side_data_ptr), scalarBetaPtr(beta_ptr) {}

  OpBrokenSpaceConstrainDHybridImpl(
      const std::string row_field,
      boost::shared_ptr<std::vector<BrokenBaseSideData>> broken_side_data_ptr,
      double beta, boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBrokenSpaceConstrainDHybridImpl(row_field, broken_side_data_ptr,
                                          boost::make_shared<double>(beta),
                                          ents_ptr) {}

private:
  boost::shared_ptr<double> scalarBetaPtr;
  boost::shared_ptr<std::vector<BrokenBaseSideData>> brokenSideDataPtr;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data);
};

template <int FIELD_DIM, typename OpBase>
MoFEMErrorCode
OpBrokenSpaceConstrainDHybridImpl<FIELD_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  FTENSOR_INDEX(FIELD_DIM, i);
  FTENSOR_INDEX(3, J);

  OP::locF.resize(row_data.getIndices().size(), false);
  OP::locF.clear();

  for (auto &bd : *brokenSideDataPtr) {
    auto t_w = this->getFTensor0IntegrationWeight();
    auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();
    auto t_row_base = row_data.getFTensor0N();
    auto t_flux = getFTensor2FromMat<FIELD_DIM, 3>(bd.getFlux());
    auto nb_base_functions = row_data.getN().size2() / 3;
    auto sense = bd.getSense();
    for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {
      auto t_vec = getFTensor1FromPtr<FIELD_DIM>(&*OP::locF.data().begin());
      size_t rr = 0;
      for (; rr != row_data.getIndices().size() / FIELD_DIM; ++rr) {
        t_vec(i) += (sense * t_w) * t_row_base * t_normal(J) * t_flux(i, J);
        ++t_row_base;
        ++t_vec;
      }
      for (; rr < nb_base_functions; ++rr)
        ++t_row_base;
      ++t_w;
      ++t_normal;
      ++t_flux;
    }
  }

  if(scalarBetaPtr)
    OP::locF *= *scalarBetaPtr;

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM

#endif // __FORMSBROKENSPACECONSTRAINTIMPL_HPP__
