/** \file TriLinearFormsIntegrators.hpp
 * \brief Trilinear forms integrators
 * \ingroup mofem_form
 */

#ifndef __TRILINEAR_FORMS_INTEGRATORS_HPP__
#define __TRILINEAR_FORMS_INTEGRATORS_HPP__

namespace MoFEM {

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, IntegrationType I,
          typename OpBase>
struct OpConvectiveTermLhsDuImpl;

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, IntegrationType I,
          typename OpBase>
struct OpConvectiveTermLhsDyImpl;

template <int SPACE_DIM, typename OpBase>
struct OpConvectiveTermLhsDuImpl<1, 1, SPACE_DIM, GAUSS, OpBase>
    : public OpBase {
  OpConvectiveTermLhsDuImpl(
      const std::string field_name_row, const std::string field_name_col,
      boost::shared_ptr<MatrixDouble> y_grad_ptr,
      ConstantFun alpha_fun = []() { return 1; })
      : OpBase(field_name_row, field_name_col, OpBase::OPROWCOL),
        yGradPtr(y_grad_ptr), alphaConstant(alpha_fun) {

    this->assembleTranspose = false;
    this->onlyTranspose = false;
    this->sYmm = false;
  }

protected:
  boost::shared_ptr<MatrixDouble> yGradPtr;
  ConstantFun alphaConstant;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);
};

template <int SPACE_DIM, typename OpBase>
struct OpConvectiveTermLhsDyImpl<1, 1, SPACE_DIM, GAUSS, OpBase>
    : public OpBase {
  OpConvectiveTermLhsDyImpl(
      const std::string field_name_row, const std::string field_name_col,
      boost::shared_ptr<MatrixDouble> u_ptr,
      ConstantFun alpha_fun = []() { return 1; })
      : OpBase(field_name_row, field_name_col, OpBase::OPROWCOL), uPtr(u_ptr),
        alphaConstant(alpha_fun) {

    this->assembleTranspose = false;
    this->onlyTranspose = false;
    this->sYmm = false;
  }

protected:
  ConstantFun alphaConstant;
  boost::shared_ptr<MatrixDouble> uPtr;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);
};

template <int FIELD_DIM, int SPACE_DIM, typename OpBase>
struct OpConvectiveTermLhsDuImpl<1, FIELD_DIM, SPACE_DIM, GAUSS, OpBase>
    : public OpBase {
  OpConvectiveTermLhsDuImpl(
      const std::string field_name_row, const std::string field_name_col,
      boost::shared_ptr<MatrixDouble> y_grad_ptr,
      ConstantFun alpha_fun = []() { return 1; })
      : OpBase(field_name_row, field_name_col, OpBase::OPROWCOL),
        yGradPtr(y_grad_ptr), alphaConstant(alpha_fun) {

    this->assembleTranspose = false;
    this->onlyTranspose = false;
    this->sYmm = false;
  }

protected:
  ConstantFun alphaConstant;
  boost::shared_ptr<MatrixDouble> yGradPtr;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);
};

template <int FIELD_DIM, int SPACE_DIM, typename OpBase>
struct OpConvectiveTermLhsDyImpl<1, FIELD_DIM, SPACE_DIM, GAUSS, OpBase>
    : public OpBase {
  OpConvectiveTermLhsDyImpl(
      const std::string field_name_row, const std::string field_name_col,
      boost::shared_ptr<MatrixDouble> u_ptr,
      ConstantFun alpha_fun = []() { return 1; })
      : OpBase(field_name_row, field_name_col, OpBase::OPROWCOL), uPtr(u_ptr),
        alphaConstant(alpha_fun) {

    this->assembleTranspose = false;
    this->onlyTranspose = false;
    this->sYmm = false;
  }

protected:
  ConstantFun alphaConstant;
  boost::shared_ptr<MatrixDouble> uPtr;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);
};

/**
 * @brief Trilinear integrator form
 * @ingroup mofem_forms
 *
 * @tparam EleOp
 * @tparam A
 * @tparam I
 */
template <typename EleOp>
template <AssemblyType A>
template <IntegrationType I>
struct FormsIntegrators<EleOp>::Assembly<A>::TriLinearForm {

  /**
   * @brief Integrate \f$(\lambda_{ij},u_{i,j})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM>
  using OpMixTensorTimesGrad = OpMixTensorTimesGradImpl<SPACE_DIM, I, OpBase>;

  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM>
  using OpConvectiveTermLhsDu =
      OpConvectiveTermLhsDuImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I, OpBase>;

  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM>
  using OpConvectiveTermLhsDy =
      OpConvectiveTermLhsDyImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I, OpBase>;
};

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode
OpConvectiveTermLhsDuImpl<1, 1, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data,
    EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();

  auto get_t_vec = [&](const int rr) {
    std::array<double *, SPACE_DIM> ptrs;
    for (auto i = 0; i != SPACE_DIM; ++i)
      ptrs[i] = &OpBase::locMat(rr, i);
    return FTensor::Tensor1<FTensor::PackPtr<double *, SPACE_DIM>, SPACE_DIM>(
        ptrs);
  };

  auto t_grad_y = getFTensor1FromMat<SPACE_DIM>(*yGradPtr);
  FTensor::Index<'i', SPACE_DIM> i;

  const double alpha_constant = alphaConstant();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha = t_w * vol * alpha_constant;
    // access local matrix
    auto t_vec = get_t_vec(0);
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows; rr++) {
      // get column base functions gradient at gauss point gg
      auto t_col_base = col_data.getFTensor0N(gg, 0);
      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols / SPACE_DIM; cc++) {
        // calculate element of local matrix
        t_vec(i) += alpha * t_row_base * t_col_base * t_grad_y(i);
        ++t_col_base;
        ++t_vec;
      }
      ++t_row_base;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;

    ++t_grad_y;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
};

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode
OpConvectiveTermLhsDyImpl<1, 1, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data,
    EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();

  auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
  FTensor::Index<'i', SPACE_DIM> i;
  const double alpha_constant = alphaConstant();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha = t_w * vol * alpha_constant;
    // loop over rows base functions
    auto a_mat_ptr = &*OpBase::locMat.data().begin();
    int rr = 0;
    for (; rr != OpBase::nbRows; rr++) {
      // get column base functions gradient at gauss point gg
      auto t_diff_col_base = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);
      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols; cc++) {
        // calculate element of local matrix
        (*a_mat_ptr) += alpha * t_row_base * t_diff_col_base(i) * t_u(i);
        ++t_diff_col_base;
        ++a_mat_ptr;
      }
      ++t_row_base;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;

    ++t_u;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
};

template <int FIELD_DIM, int SPACE_DIM, typename OpBase>
MoFEMErrorCode
OpConvectiveTermLhsDuImpl<1, FIELD_DIM, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data,
    EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();

  auto t_grad_y = getFTensor2FromMat<FIELD_DIM, SPACE_DIM>(*yGradPtr);
  FTensor::Index<'I', FIELD_DIM> I;
  FTensor::Index<'k', SPACE_DIM> k;

  auto get_t_mat = [&](const int rr) {
    std::array<double *, FIELD_DIM * SPACE_DIM> ptrs;
    int s = 0;
    for (int j = 0; j != FIELD_DIM; ++j)
      for (auto i = 0; i != SPACE_DIM; ++i, ++s)
        ptrs[s] = &OpBase::locMat(rr + j, i);
    return FTensor::Tensor2<FTensor::PackPtr<double *, SPACE_DIM>, FIELD_DIM,
                            SPACE_DIM>(ptrs);
  };

  const double alpha_constant = alphaConstant();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha = t_w * vol * alpha_constant;

    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows / FIELD_DIM; rr++) {
      // get column base functions gradient at gauss point gg
      auto t_col_base = col_data.getFTensor0N(gg, 0);
      // get mat
      auto t_mat = get_t_mat(FIELD_DIM * rr);
      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols / SPACE_DIM; cc++) {
        // calculate element of local matrix
        t_mat(I, k) += (alpha * t_row_base * t_col_base) * t_grad_y(I, k);
        ++t_col_base;
        ++t_mat;
      }
      ++t_row_base;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;

    ++t_grad_y;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
};

template <int FIELD_DIM, int SPACE_DIM, typename OpBase>
MoFEMErrorCode
OpConvectiveTermLhsDyImpl<1, FIELD_DIM, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data,
    EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();

  auto get_t_mat = [&](const int rr) {
    std::array<double *, FIELD_DIM * FIELD_DIM> ptrs;
    int s = 0;
    for (int i = 0; i != FIELD_DIM; ++i)
      for (int j = 0; j != FIELD_DIM; ++j, ++s)
        ptrs[s] = &(OpBase::locMat(rr + i, j));
    return FTensor::Tensor2<FTensor::PackPtr<double *, FIELD_DIM>, FIELD_DIM,
                            FIELD_DIM>(ptrs);
  };

  auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
  constexpr auto t_kd = FTensor::Kronecker_Delta_symmetric<int>();

  FTensor::Index<'I', FIELD_DIM> I;
  FTensor::Index<'L', FIELD_DIM> L;
  FTensor::Index<'k', SPACE_DIM> k;

  const double alpha_constant = alphaConstant();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha = t_w * vol * alpha_constant;

    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows / FIELD_DIM; ++rr) {
      // get matrix vec
      auto t_mat = get_t_mat(FIELD_DIM * rr);
      // get column base functions gradient at gauss point gg
      auto t_diff_col_base = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);
      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols / FIELD_DIM; ++cc) {
        t_mat(I, L) +=
            alpha * t_row_base * t_kd(I, L) * (t_diff_col_base(k) * t_u(k));
        ++t_mat;
        ++t_diff_col_base;
      }
      ++t_row_base;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;

    ++t_u;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
};

} // namespace MoFEM

#endif //__TRILINEAR_FORMS_INTEGRATORS_HPP__