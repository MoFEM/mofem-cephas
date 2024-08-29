/**
 * \file FiniteThermalOps.hpp
 * \example FiniteThermalOps.hpp
 *
 * @copyright Copyright (c) 2023
 */

#ifndef __FINITE_THERMAL_OPS_HPP__
#define __FINITE_THERMAL_OPS_HPP__

namespace FiniteThermalOps {

struct CommonData : public boost::enable_shared_from_this<CommonData> {
  boost::shared_ptr<double> thermalConductivityPtr;
  boost::shared_ptr<double> heatCapacityPtr;
  boost::shared_ptr<double> coeffExpansionPtr;
  boost::shared_ptr<double> refTempPtr;
};

template <int DIM, IntegrationType I, typename DomainEleOp>
struct OpCalculateQdotQRhs;

// Calculate 1/J * K^-1 * Q * deltaQ
// Calculate 1/J * (F_Ji^T * k_im^-1 * FmN) * Q_N * deltaQ_J
template <int DIM, typename DomainEleOp>
struct OpCalculateQdotQRhs<DIM, GAUSS, DomainEleOp> : public DomainEleOp {

  OpCalculateQdotQRhs(const std::string field_name,
                      boost::shared_ptr<MatrixDouble> vec,
                      ScalarFun resistance_function,
                      boost::shared_ptr<MatrixDouble> mat_Grad_Ptr,
                      boost::shared_ptr<Range> ents_ptr = nullptr)
      : DomainEleOp(field_name, field_name, DomainEleOp::OPROW), fluxVec(vec),
        resistanceFunction(resistance_function), matGradPtr(mat_Grad_Ptr) {}

protected:
  boost::shared_ptr<MatrixDouble> fluxVec;
  ScalarFun resistanceFunction;
  boost::shared_ptr<MatrixDouble> matGradPtr;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int DIM, typename DomainEleOp>
MoFEMErrorCode OpCalculateQdotQRhs<DIM, GAUSS, DomainEleOp>::iNtegrate(
    EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  FTensor::Index<'i', DIM> i;
  FTensor::Index<'J', DIM> J;
  FTensor::Index<'m', DIM> m;
  FTensor::Index<'N', DIM> N;

  const auto nb_integration_pts = DomainEleOp::getGaussPts().size2();

  FTensor::Tensor2<double, DIM, DIM> t_F;
  FTensor::Tensor1<double, DIM> t_K_inv_Q_over_J;

  constexpr auto t_kd = FTensor::Kronecker_Delta<int>();

  auto t_grad = getFTensor2FromMat<DIM, DIM>(*(matGradPtr));

  // get element volume
  const double vol = DomainEleOp::getMeasure();
  // get integration weights
  auto t_w = DomainEleOp::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = data.getFTensor1N<3>();
  // get flux values
  auto t_flux = getFTensor1FromMat<DIM>(*fluxVec);
  // get coordinate at integration points
  auto t_coords = AssemblyDomainEleOp::getFTensor1CoordsAtGaussPts();

  // loop over integration points
  for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
    // take into account Jacobian
    const double alpha = t_w * vol;

    t_F(i, J) = t_grad(i, J) + t_kd(i, J); // Deformation gradient
    auto t_J = determinantTensor(t_F);     // Volume jacobian

    t_K_inv_Q_over_J(J) =
        resistanceFunction(t_coords(0), t_coords(1), t_coords(2)) * t_F(m, J) *
        t_F(m, N) * t_flux(N) / t_J; // Value which multiplies delta Q

    // loop over rows base functions
    size_t rr = 0;
    for (; rr != DomainEleOp::nbRows; ++rr) {
      DomainEleOp::locF[rr] += t_K_inv_Q_over_J(J) * alpha * (t_row_base(J));
      ++t_row_base;
    }
    for (; rr < DomainEleOp::nbRowBaseFunctions; ++rr)
      ++t_row_base;
    ++t_w; // move to another integration weight
    ++t_coords;
    ++t_flux;
    ++t_grad;
  };
  MoFEMFunctionReturn(0);
}

template <int DIM, IntegrationType I, typename AssemblyDomainEleOp>
struct OpCalculateQdotQLhs_dQ;

// Calculate LHS of OpCalculateQdotQRhs wrt Q
template <int DIM, typename AssemblyDomainEleOp>
struct OpCalculateQdotQLhs_dQ<DIM, GAUSS, AssemblyDomainEleOp>
    : public AssemblyDomainEleOp {

  OpCalculateQdotQLhs_dQ(const std::string row_field_name,
                         const std::string col_field_name,
                         ScalarFun resistance_function,
                         boost::shared_ptr<MatrixDouble> mat_Grad_Ptr,
                         boost::shared_ptr<Range> ents_ptr = nullptr)
      : AssemblyDomainEleOp(row_field_name, col_field_name,
                            AssemblyDomainEleOp::OPROWCOL),
        resistanceFunction(resistance_function), matGradPtr(mat_Grad_Ptr) {}

protected:
  ScalarFun resistanceFunction;
  boost::shared_ptr<MatrixDouble> matGradPtr;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);
};

template <int DIM, typename AssemblyDomainEleOp>
MoFEMErrorCode
OpCalculateQdotQLhs_dQ<DIM, GAUSS, AssemblyDomainEleOp>::iNtegrate(
    EntitiesFieldData::EntData &row_data,
    EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  FTensor::Index<'i', DIM> i;
  FTensor::Index<'J', DIM> J;
  FTensor::Index<'M', DIM> M;
  FTensor::Index<'L', DIM> L;

  auto t_w = this->getFTensor0IntegrationWeight();

  auto t_row_base = row_data.getFTensor1N<3>();

  FTensor::Tensor2<double, DIM, DIM> t_F;
  FTensor::Tensor2<double, DIM, DIM> t_right_Cauchy_Green;
  FTensor::Tensor2<double, DIM, DIM> t_Kinv_over_J;
  FTensor::Tensor1<double, DIM> t_Kinv_over_J_row_base;

  constexpr auto t_kd = FTensor::Kronecker_Delta<int>();

  auto t_grad = getFTensor2FromMat<DIM, DIM>(*(matGradPtr));

  // get coordinate at integration points
  auto t_coords = AssemblyDomainEleOp::getFTensor1CoordsAtGaussPts();

  for (size_t gg = 0; gg != AssemblyDomainEleOp::nbIntegrationPts; ++gg) {

    t_F(i, J) = t_grad(i, J) + t_kd(i, J); // Deformation gradient
    auto t_J = determinantTensor(t_F);     // Volume jacobian

    t_right_Cauchy_Green(M, L) =
        t_F(i, M) *
        t_F(i, L); // The inverse of the material thermal conductivity tensor

    const double alpha = this->getMeasure() * t_w;

    t_Kinv_over_J(M, L) =
        (alpha * resistanceFunction(t_coords(0), t_coords(1), t_coords(2)) /
         t_J) *
        t_right_Cauchy_Green(M, L);

    size_t rr = 0;
    for (; rr != AssemblyDomainEleOp::nbRows; ++rr) {
      t_Kinv_over_J_row_base(L) = t_Kinv_over_J(M, L) * t_row_base(M);
      auto t_col_base = col_data.getFTensor1N<3>(gg, 0);
      for (size_t cc = 0; cc != AssemblyDomainEleOp::nbCols; ++cc) {
        this->locMat(rr, cc) += t_Kinv_over_J_row_base(L) * t_col_base(L);
        ++t_col_base;
      }
      ++t_row_base;
    }
    for (; rr < AssemblyDomainEleOp::nbRowBaseFunctions; ++rr)
      ++t_row_base;

    ++t_w; // move to another integration weight
    ++t_coords;
    ++t_grad;
  }

  MoFEMFunctionReturn(0);
}

template <int DIM, IntegrationType I, typename AssemblyDomainEleOp,
          bool IS_LARGE_STRAINS>
struct OpCalculateQdotQLhs_dU;

// Calculate LHS of OpCalculateQdotQRhs wrt U
template <int DIM, typename AssemblyDomainEleOp>
struct OpCalculateQdotQLhs_dU<DIM, GAUSS, AssemblyDomainEleOp, true>
    : public AssemblyDomainEleOp {

  OpCalculateQdotQLhs_dU(const std::string row_field_name,
                         const std::string col_field_name,
                         boost::shared_ptr<MatrixDouble> vec,
                         ScalarFun resistance_function,
                         boost::shared_ptr<MatrixDouble> mat_Grad_Ptr,
                         boost::shared_ptr<Range> ents_ptr = nullptr)
      : AssemblyDomainEleOp(row_field_name, col_field_name,
                            AssemblyDomainEleOp::OPROWCOL),
        fluxVec(vec), resistanceFunction(resistance_function),
        matGradPtr(mat_Grad_Ptr) {
    AssemblyDomainEleOp::sYmm = false;
  }

protected:
  boost::shared_ptr<MatrixDouble> fluxVec;
  ScalarFun resistanceFunction;
  boost::shared_ptr<MatrixDouble> matGradPtr;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);
};

template <int DIM, typename AssemblyDomainEleOp>
MoFEMErrorCode
OpCalculateQdotQLhs_dU<DIM, GAUSS, AssemblyDomainEleOp, true>::iNtegrate(
    EntitiesFieldData::EntData &row_data,
    EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  FTensor::Index<'i', DIM> i;
  FTensor::Index<'J', DIM> J;
  FTensor::Index<'m', DIM> m;
  FTensor::Index<'M', DIM> M;
  FTensor::Index<'l', DIM> l;
  FTensor::Index<'L', DIM> L;

  auto t_w = this->getFTensor0IntegrationWeight();

  // Always need three components for vectorial basis
  auto t_row_base = row_data.getFTensor1N<3>();

  FTensor::Tensor2<double, DIM, DIM> t_F;
  FTensor::Tensor2<double, DIM, DIM> t_right_Cauchy_Green;
  FTensor::Tensor2<double, DIM, DIM> t_F_inv;
  FTensor::Tensor1<double, DIM> t_F_flux;
  FTensor::Tensor1<double, DIM> t_F_row_base;
  FTensor::Tensor1<double, DIM> t_FtF_Q;
  FTensor::Tensor2<double, DIM, DIM> t_FtF_Finv_Q_row_base;
  FTensor::Tensor2<double, DIM, DIM> t_F_Q_row_base;

  constexpr auto t_kd = FTensor::Kronecker_Delta<int>();

  auto t_grad = getFTensor2FromMat<DIM, DIM>(*(matGradPtr));

  // get flux values
  auto t_flux = getFTensor1FromMat<DIM>(*fluxVec);

  // for each integration point
  for (size_t gg = 0; gg != AssemblyDomainEleOp::nbIntegrationPts; ++gg) {

    t_F(i, J) = t_grad(i, J) + t_kd(i, J); // Deformation gradient
    auto t_J = determinantTensor(t_F);     // Volume jacobian
    CHKERR invertTensor(t_F, t_J, t_F_inv);

    t_right_Cauchy_Green(M, L) =
        t_F(i, M) *
        t_F(i, L); // The inverse of the material thermal conductivity tensor

    const double alpha = this->getMeasure() * t_w;

    // get the vector of the local matrix for the derivative wrt a vectorial
    // field with a scalar basis
    auto t_vec = getFTensor1FromPtr<DIM>(&this->locMat(0, 0));

    // get coordinate at integration points
    auto t_coords = AssemblyDomainEleOp::getFTensor1CoordsAtGaussPts();

    auto alpha_resistance_Jinv =
        alpha * resistanceFunction(t_coords(0), t_coords(1), t_coords(2)) / t_J;
    t_F_flux(l) = t_F(l, L) * t_flux(L);

    t_FtF_Q(M) = t_right_Cauchy_Green(M, L) * t_flux(L);

    // for each row coefficient
    size_t rr = 0;
    for (; rr != AssemblyDomainEleOp::nbRows; ++rr) {
      auto t_col_grad_base = col_data.getFTensor1DiffN<DIM>(gg, 0);

      t_FtF_Finv_Q_row_base(J, l) = t_FtF_Q(M) * t_row_base(M) * t_F_inv(J, l);
      t_F_row_base(l) = t_F(l, M) * t_row_base(M);
      t_F_Q_row_base(l, M) = t_F_flux(l) * t_row_base(M);

      // for each coupled coefficient in the columns
      // divied by DIM since scalar basis for vectorial field
      for (size_t cc = 0; cc != AssemblyDomainEleOp::nbCols / DIM; ++cc) {
        t_vec(l) += alpha_resistance_Jinv *
                    (-t_FtF_Finv_Q_row_base(J, l) * t_col_grad_base(J) +
                     t_F_Q_row_base(l, M) * t_col_grad_base(M) +
                     t_F_row_base(l) * t_col_grad_base(L) * t_flux(L));
        ++t_col_grad_base;
        ++t_vec;
      }
      ++t_row_base;
    }
    for (; rr < AssemblyDomainEleOp::nbRowBaseFunctions; ++rr)
      ++t_row_base;

    ++t_w; // move to another integration weight
    ++t_coords;
    ++t_flux;
    ++t_grad;
  }

  MoFEMFunctionReturn(0);
}

template <int DIM, typename AssemblyDomainEleOp>
struct OpCalculateQdotQLhs_dU<DIM, GAUSS, AssemblyDomainEleOp, false>
    : public AssemblyDomainEleOp {

  OpCalculateQdotQLhs_dU(const std::string row_field_name,
                         const std::string col_field_name,
                         boost::shared_ptr<MatrixDouble> vec,
                         ScalarFun resistance_function,
                         boost::shared_ptr<MatrixDouble> mat_Grad_Ptr,
                         boost::shared_ptr<Range> ents_ptr = nullptr)
      : AssemblyDomainEleOp(row_field_name, col_field_name,
                            AssemblyDomainEleOp::OPROWCOL),
        fluxVec(vec), resistanceFunction(resistance_function),
        matGradPtr(mat_Grad_Ptr) {
    AssemblyDomainEleOp::sYmm = false;
  }

protected:
  boost::shared_ptr<MatrixDouble> fluxVec;
  ScalarFun resistanceFunction;
  boost::shared_ptr<MatrixDouble> matGradPtr;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);
};

template <int DIM, typename AssemblyDomainEleOp>
MoFEMErrorCode
OpCalculateQdotQLhs_dU<DIM, GAUSS, AssemblyDomainEleOp, false>::iNtegrate(
    EntitiesFieldData::EntData &row_data,
    EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  // Empty implementation of OpCalculateQdotQLhs_dU for small the case of small
  // strains. This allows the operator to be called without adding any
  // contribution to the stiffness matrix.

  MoFEMFunctionReturn(0);
}

// Templated on IS_LARGE_STRAINS to allow different implementations with
// different numbers of arguments at compile time
template <int SPACE_DIM, bool IS_LARGE_STRAINS> struct OpHdivHdivImpl;

template <int SPACE_DIM>
struct OpHdivHdivImpl<SPACE_DIM, false>
    : public FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
          GAUSS>::OpMass<3, SPACE_DIM> {
  OpHdivHdivImpl(const std::string row_field_name,
                 const std::string col_field_name,
                 ScalarFun resistance_function,
                 boost::shared_ptr<MatrixDouble> mat_Grad_Ptr,
                 boost::shared_ptr<Range> ents_ptr = nullptr)
      : FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
            GAUSS>::OpMass<3, SPACE_DIM>(row_field_name, col_field_name,
                                         resistance_function, ents_ptr) {}
};

template <int SPACE_DIM>
struct OpHdivHdivImpl<SPACE_DIM, true>
    : public OpCalculateQdotQLhs_dQ<SPACE_DIM, GAUSS, AssemblyDomainEleOp> {
  OpHdivHdivImpl(const std::string row_field_name,
                 const std::string col_field_name,
                 ScalarFun resistance_function,
                 boost::shared_ptr<MatrixDouble> mat_Grad_Ptr,
                 boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpCalculateQdotQLhs_dQ<SPACE_DIM, GAUSS, AssemblyDomainEleOp>(
            row_field_name, col_field_name, resistance_function, mat_Grad_Ptr,
            ents_ptr) {}
};

// Templated on IS_LARGE_STRAINS to allow different implementations with
// different numbers of arguments at compile time
template <int SPACE_DIM, bool IS_LARGE_STRAINS> struct OpHdivFluxImpl;

template <int SPACE_DIM>
struct OpHdivFluxImpl<SPACE_DIM, false>
    : public FormsIntegrators<DomainEleOp>::Assembly<PETSC>::LinearForm<
          GAUSS>::OpBaseTimesVector<3, SPACE_DIM, 1> {
  OpHdivFluxImpl(const std::string field_name,
                 boost::shared_ptr<MatrixDouble> vec,
                 ScalarFun resistance_function,
                 boost::shared_ptr<MatrixDouble> mat_Grad_Ptr,
                 boost::shared_ptr<Range> ents_ptr = nullptr)
      : FormsIntegrators<DomainEleOp>::Assembly<PETSC>::LinearForm<
            GAUSS>::OpBaseTimesVector<3, SPACE_DIM, 1>(field_name, vec,
                                                       resistance_function,
                                                       ents_ptr) {}
};
template <int SPACE_DIM>
struct OpHdivFluxImpl<SPACE_DIM, true>
    : public OpCalculateQdotQRhs<SPACE_DIM, GAUSS, AssemblyDomainEleOp> {
  OpHdivFluxImpl(const std::string field_name,
                 boost::shared_ptr<MatrixDouble> vec,
                 ScalarFun resistance_function,
                 boost::shared_ptr<MatrixDouble> mat_Grad_Ptr,
                 boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpCalculateQdotQRhs<SPACE_DIM, GAUSS, AssemblyDomainEleOp>(
            field_name, vec, resistance_function, mat_Grad_Ptr, ents_ptr) {}
};

} // namespace FiniteThermalOps

#endif // __FINITE_THERMAL_OPS_HPP__