/** \file LinearFormsIntegrators.hpp
  * \brief Linear forms integrators
  * \ingroup mofem_form

*/

#ifndef __LINEAR_FORMS_INTEGRATORS_HPP__
#define __LINEAR_FORMS_INTEGRATORS_HPP__

namespace MoFEM {

struct SourceFunctionSpecialization {
  template <typename OpBase> struct S { S() = delete; };
  SourceFunctionSpecialization() = delete;
};

struct SourceBoundaryNormalSpecialization {
  template <typename OpBase> struct S { S() = delete; };
  SourceBoundaryNormalSpecialization() = delete;
};

template <int BASE_DIM, int FIELD_DIM, IntegrationType I, typename OpBase>
struct OpSourceImpl;

/**
 * @brief Integrate source
 *
 * @tparam OpBase
 */
template <typename OpBase>
struct OpSourceImpl<1, 1, GAUSS, SourceFunctionSpecialization::S<OpBase>>
    : public OpBase {

  /**
   * @brief Construct a new Op Source Impl object
   *
   * @param field_name
   * @param source_fun
   * @param ents_ptr
   */
  OpSourceImpl(const std::string field_name, ScalarFun source_fun,
               boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, ents_ptr),
        sourceFun(source_fun) {}

  /**
   * @brief Construct a new Op Source Impl object
   *
   * @param field_name
   * @param time_fun
   * @param source_fun
   * @param ents_ptr
   */
  OpSourceImpl(const std::string field_name, TimeFun time_fun,
               ScalarFun source_fun,
               boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, time_fun, ents_ptr),
        sourceFun(source_fun) {}

protected:
  ScalarFun sourceFun;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int FIELD_DIM, typename OpBase>
struct OpSourceImpl<1, FIELD_DIM, GAUSS,
                    SourceFunctionSpecialization::S<OpBase>> : public OpBase {

  /**
   * @brief Construct a new Op Source Impl object
   *
   * @param field_name
   * @param time_fun
   * @param source_fun
   * @param ents_ptr
   */
  OpSourceImpl(const std::string field_name, TimeFun time_fun,
               VectorFun<FIELD_DIM> source_fun,
               boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, time_fun, ents_ptr) {}

  /**
   * @brief Construct a new Op Source Impl object
   *
   * @param field_name
   * @param source_fun
   * @param ents_ptr
   */
  OpSourceImpl(const std::string field_name, VectorFun<FIELD_DIM> source_fun,
               boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, ents_ptr),
        sourceFun(source_fun) {}

protected:
  VectorFun<FIELD_DIM> sourceFun;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int FIELD_DIM, typename OpBase>
struct OpSourceImpl<3, FIELD_DIM, GAUSS,
                    SourceFunctionSpecialization::S<OpBase>> : public OpBase {

  OpSourceImpl(const std::string field_name, TimeFun time_fun,
               VectorFun<FIELD_DIM> source_fun,
               boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, time_fun, ents_ptr),
        sourceFun(source_fun) {}

  OpSourceImpl(const std::string field_name, VectorFun<FIELD_DIM> source_fun,
               boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, ents_ptr),
        sourceFun(source_fun) {}

protected:
  VectorFun<FIELD_DIM> sourceFun;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int BASE_DIM, int S, IntegrationType I, typename OpBase>
struct OpBaseTimesScalarImpl;

template <int S, typename OpBase>
struct OpBaseTimesScalarImpl<1, S, GAUSS, OpBase> : public OpBase {

  OpBaseTimesScalarImpl(
      const std::string field_name, boost::shared_ptr<VectorDouble> vec,
      ScalarFun beta_coeff = [](double, double, double) constexpr { return 1; },
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, ents_ptr), sourceVec(vec),
        betaCoeff(beta_coeff) {}

protected:
  ScalarFun betaCoeff;
  boost::shared_ptr<VectorDouble> sourceVec;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int BASE_DIM, int FIELD_DIM, int S, IntegrationType I,
          typename OpBase>
struct OpBaseTimesVectorImpl;

template <int FIELD_DIM, int S, typename OpBase>
struct OpBaseTimesVectorImpl<1, FIELD_DIM, S, GAUSS, OpBase> : public OpBase {

  OpBaseTimesVectorImpl(
      const std::string field_name, boost::shared_ptr<MatrixDouble> vec,
      ScalarFun beta_coeff = [](double, double, double) constexpr { return 1; },
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, ents_ptr), sourceVec(vec),
        betaCoeff(beta_coeff) {}

protected:
  ScalarFun betaCoeff;
  boost::shared_ptr<MatrixDouble> sourceVec;
  FTensor::Index<'i', FIELD_DIM> i;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int FIELD_DIM, int S, typename OpBase>
struct OpBaseTimesVectorImpl<3, FIELD_DIM, S, GAUSS, OpBase> : public OpBase {

  OpBaseTimesVectorImpl(
      const std::string field_name, boost::shared_ptr<MatrixDouble> vec,
      ScalarFun beta_coeff = [](double, double, double) constexpr { return 1; },
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, ents_ptr), sourceVec(vec),
        betaCoeff(beta_coeff) {}

protected:
  ScalarFun betaCoeff;
  boost::shared_ptr<MatrixDouble> sourceVec;
  FTensor::Index<'i', FIELD_DIM> i;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S, IntegrationType I,
          typename OpBase>
struct OpGradTimesTensorImpl;

template <int SPACE_DIM, int S, typename OpBase>
struct OpGradTimesTensorImpl<1, 1, SPACE_DIM, S, GAUSS, OpBase>
    : public OpBase {

  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index

  OpGradTimesTensorImpl(
      const std::string field_name, boost::shared_ptr<MatrixDouble> mat_vals,
      ScalarFun beta_coeff = [](double, double, double) constexpr { return 1; },
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, ents_ptr),
        matVals(mat_vals), betaCoeff(beta_coeff) {}

protected:
  boost::shared_ptr<MatrixDouble> matVals;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
  ScalarFun betaCoeff;
};

template <int SPACE_DIM, int S, typename OpBase>
struct OpGradTimesTensorImpl<1, SPACE_DIM, SPACE_DIM, S, GAUSS, OpBase>
    : public OpBase {

  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  FTensor::Index<'j', SPACE_DIM> j; ///< summit Index

  OpGradTimesTensorImpl(
      const std::string field_name, boost::shared_ptr<MatrixDouble> mat_vals,
      ScalarFun beta_coeff = [](double, double, double) constexpr { return 1; },
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, ents_ptr),
        matVals(mat_vals), betaCoeff(beta_coeff) {}

protected:
  boost::shared_ptr<MatrixDouble> matVals;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
  ScalarFun betaCoeff;
};

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S, IntegrationType I,
          typename OpBase>
struct OpGradTimesSymTensorImpl {};

template <int SPACE_DIM, int S, typename OpBase>
struct OpGradTimesSymTensorImpl<1, SPACE_DIM, SPACE_DIM, S, GAUSS, OpBase>
    : public OpBase {

  OpGradTimesSymTensorImpl(
      const std::string field_name, boost::shared_ptr<MatrixDouble> mat_vals,
      ScalarFun beta_coeff = [](double, double, double) constexpr { return 1; })
      : OpBase(field_name, field_name, OpBase::OPROW), matVals(mat_vals),
        betaCoeff(beta_coeff) {}

protected:
  boost::shared_ptr<MatrixDouble> matVals;
  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'j', SPACE_DIM> j;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
  ScalarFun betaCoeff;
};

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, IntegrationType I,
          typename OpBase, CoordinateTypes CoordSys>
struct OpMixDivTimesUImpl {};

template <int FIELD_DIM, int SPACE_DIM, typename OpBase,
          CoordinateTypes CoordSys>
struct OpMixDivTimesUImpl<3, FIELD_DIM, SPACE_DIM, GAUSS, OpBase, CoordSys>
    : public OpBase {
  OpMixDivTimesUImpl(
      const std::string field_name, boost::shared_ptr<MatrixDouble> mat_vals,
      ScalarFun beta = [](double, double, double) { return 1; },
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, ents_ptr),
        matVals(mat_vals), betaConst(beta) {}

protected:
  ScalarFun betaConst;
  boost::shared_ptr<MatrixDouble> matVals;
  FTensor::Index<'i', FIELD_DIM> i;
  FTensor::Index<'j', SPACE_DIM> j;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int SPACE_DIM, typename OpBase, CoordinateTypes CoordSys>
struct OpMixDivTimesUImpl<3, 1, SPACE_DIM, GAUSS, OpBase, CoordSys>
    : public OpBase {
  OpMixDivTimesUImpl(
      const std::string field_name, boost::shared_ptr<VectorDouble> vec_vals,
      ScalarFun beta = [](double, double, double) constexpr { return 1; },
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, ents_ptr),
        vecVals(vec_vals), betaConst(beta) {}

protected:
  ScalarFun betaConst;
  boost::shared_ptr<VectorDouble> vecVals;
  FTensor::Index<'j', SPACE_DIM> j;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int FIELD_DIM, typename OpBase, CoordinateTypes CoordSys>
struct OpMixDivTimesUImpl<1, FIELD_DIM, FIELD_DIM, GAUSS, OpBase, CoordSys>
    : public OpBase {

  OpMixDivTimesUImpl(
      const std::string field_name, boost::shared_ptr<VectorDouble> vec,
      ScalarFun beta = [](double, double, double) constexpr { return 1; },
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, ents_ptr), sourceVec(vec),
        betaCoeff(beta) {}

protected:
  ScalarFun betaCoeff;
  boost::shared_ptr<VectorDouble> sourceVec;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int SPACE_DIM, IntegrationType I, typename OpBase>
struct OpMixVecTimesDivLambdaImpl {};

template <int SPACE_DIM, typename OpBase>
struct OpMixVecTimesDivLambdaImpl<SPACE_DIM, GAUSS, OpBase> : public OpBase {
  OpMixVecTimesDivLambdaImpl(const std::string field_name,
                             boost::shared_ptr<MatrixDouble> mat_vals)
      : OpBase(field_name, field_name, OpBase::OPROW), matVals(mat_vals) {}

  OpMixVecTimesDivLambdaImpl(const std::string field_name,
                             boost::shared_ptr<MatrixDouble> mat_vals,
                             ScalarFun beta_fun)
      : OpBase(field_name, field_name, OpBase::OPROW), matVals(mat_vals),
        betaCoeff(beta_fun) {}

protected:
  boost::shared_ptr<MatrixDouble> matVals;
  ScalarFun betaCoeff = [](double, double, double) constexpr { return 1; };

  FTensor::Index<'i', SPACE_DIM> i;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int SPACE_DIM, IntegrationType I, typename OpBase>
struct OpMixTensorTimesGradUImpl {};

template <int SPACE_DIM, typename OpBase>
struct OpMixTensorTimesGradUImpl<SPACE_DIM, GAUSS, OpBase> : public OpBase {
  OpMixTensorTimesGradUImpl(const std::string field_name,
                            boost::shared_ptr<MatrixDouble> mat_vals)
      : OpBase(field_name, field_name, OpBase::OPROW), matVals(mat_vals) {}

  OpMixTensorTimesGradUImpl(const std::string field_name,
                            boost::shared_ptr<MatrixDouble> mat_vals,
                            ScalarFun beta_fun)
      : OpBase(field_name, field_name, OpBase::OPROW), matVals(mat_vals),
        betaCoeff(beta_fun) {}

protected:
  boost::shared_ptr<MatrixDouble> matVals;
  ScalarFun betaCoeff = [](double, double, double) constexpr { return 1; };

  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'j', SPACE_DIM> j;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

/**
 * @brief Multiply vector times normal on the face times scalar function
 *
 * This operator typically will be used to evaluate natural boundary conditions
 * for mixed formulation.
 *
 * @tparam BASE_DIM
 * @tparam SPACE_DIM
 * @tparam OpBase
 */
template <int SPACE_DIM, IntegrationType I, typename OpBase>
struct OpNormalMixVecTimesScalarImpl;

/**
 * @brief This is specialisation for sources on boundary which depends on normal
 *
 * @tparam BASE_DIM
 * @tparam OpBase
 */
template <int FIELD_DIM, IntegrationType I, typename OpBase>
struct OpSourceImpl<3, FIELD_DIM, I,
                    SourceBoundaryNormalSpecialization::S<OpBase>>
    : public OpNormalMixVecTimesScalarImpl<FIELD_DIM, I, OpBase> {
  using OpNormalMixVecTimesScalarImpl<FIELD_DIM, I,
                                      OpBase>::OpNormalMixVecTimesScalarImpl;
};

template <typename OpBase>
struct OpNormalMixVecTimesScalarImpl<3, GAUSS, OpBase> : public OpBase {
  OpNormalMixVecTimesScalarImpl(
      const std::string field_name,
      ScalarFun source_fun = [](double, double, double) constexpr { return 1; },
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, ents_ptr),
        sourceFun(source_fun) {}

protected:
  ScalarFun sourceFun;
  FTensor::Index<'i', 3> i;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <typename OpBase>
struct OpNormalMixVecTimesScalarImpl<2, GAUSS, OpBase> : public OpBase {
  OpNormalMixVecTimesScalarImpl(
      const std::string field_name,
      ScalarFun source_fun = [](double, double, double) constexpr { return 1; },
      boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW, ents_ptr),
        sourceFun(source_fun) {}

protected:
  ScalarFun sourceFun;
  FTensor::Index<'i', 3> i;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, IntegrationType I,
          typename OpBase>
struct OpConvectiveTermRhsImpl;

template <int SPACE_DIM, typename OpBase>
struct OpConvectiveTermRhsImpl<1, 1, SPACE_DIM, GAUSS, OpBase> : public OpBase {
  OpConvectiveTermRhsImpl(
      const std::string field_name, boost::shared_ptr<MatrixDouble> u_ptr,
      boost::shared_ptr<MatrixDouble> y_grad_ptr,
      ConstantFun source_fun = []() constexpr { return 1; })
      : OpBase(field_name, field_name, OpBase::OPROW), uPtr(u_ptr),
        yGradPtr(y_grad_ptr), alphaConstant(source_fun) {}

protected:
  boost::shared_ptr<MatrixDouble> uPtr;
  boost::shared_ptr<MatrixDouble> yGradPtr;
  ConstantFun alphaConstant;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int FIELD_DIM, int SPACE_DIM, typename OpBase>
struct OpConvectiveTermRhsImpl<1, FIELD_DIM, SPACE_DIM, GAUSS, OpBase>
    : public OpBase {
  OpConvectiveTermRhsImpl(
      const std::string field_name, boost::shared_ptr<MatrixDouble> u_ptr,
      boost::shared_ptr<MatrixDouble> y_grad_ptr,
      ConstantFun source_fun = []() constexpr { return 1; })
      : OpBase(field_name, field_name, OpBase::OPROW), uPtr(u_ptr),
        yGradPtr(y_grad_ptr), alphaConstant(source_fun) {}

protected:
  boost::shared_ptr<MatrixDouble> uPtr;
  boost::shared_ptr<MatrixDouble> yGradPtr;
  ConstantFun alphaConstant;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

/**
 * @brief Linear integrator form
 * @ingroup mofem_forms
 *
 * @tparam EleOp
 * @tparam A
 * @tparam I
 */
template <typename EleOp>
template <AssemblyType A>
template <IntegrationType I>
struct FormsIntegrators<EleOp>::Assembly<A>::LinearForm {

  /**
   * @brief Integrate \f$(v,f(\mathbf{x}))_\Omega\f$, f is a scalar
   * @ingroup mofem_forms
   *
   * @note \f$f(\mathbf{x})\$ is scalar function of coordinates
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   */
  template <int BASE_DIM, int FIELD_DIM,
            typename S = SourceFunctionSpecialization>
  using OpSource =
      OpSourceImpl<BASE_DIM, FIELD_DIM, I, typename S::template S<OpBase>>;

  /**
   * @brief Vector field integrator \f$(v_i,f)_\Omega\f$, f is a vector
   * @ingroup mofem_forms
   *
   * @tparam BASE_DIM
   */
  template <int BASE_DIM, int S = 1>
  using OpBaseTimesScalar = OpBaseTimesScalarImpl<BASE_DIM, S, I, OpBase>;

  /** @deprecated use instead OpBaseTimesScalar
   */
  template <int BASE_DIM, int S = 1>
  using OpBaseTimesScalarField = OpBaseTimesScalar<BASE_DIM, S>;

  /**
   * @brief Vector field integrator \f$(v,f_i)_\Omega\f$, f is a vector
   * @ingroup mofem_forms
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   * @tparam 0
   */
  template <int BASE_DIM, int FIELD_DIM, int S>
  using OpBaseTimesVector =
      OpBaseTimesVectorImpl<BASE_DIM, FIELD_DIM, S, I, OpBase>;
  //! [Source operator]

  //! [Grad times tensor]

  /**
   * @brief Integrate \f$(v_{,i},f_{ij})_\Omega\f$, f is a vector
   * @ingroup mofem_forms
   *
   * @note \f$f_{ij}\$ is tensor at integration points
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   * @tparam SPACE_DIM
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S = 1>
  using OpGradTimesTensor =
      OpGradTimesTensorImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I, OpBase>;

  /**
   * @brief Integrate \f$(v_{,i},f_{ij})_\Omega\f$, f is symmetric tensor
   * @ingroup mofem_forms
   *
   * @note \f$f_{ij}\$ is tensor at integration points
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   * @tparam SPACE_DIM
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S = 1>
  using OpGradTimesSymTensor =
      OpGradTimesSymTensorImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I, OpBase>;

  /**
   * @brief Integrate \f$(\lambda_{ij,j},u_{i})_\Omega\f$
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   * @tparam SPACE_DIM
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM,
            CoordinateTypes CoordSys = CARTESIAN>
  using OpMixDivTimesU =
      OpMixDivTimesUImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I, OpBase, CoordSys>;

  /**
   * @brief Integrate \f$(\lambda_{ij},u_{i,j})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM>
  using OpMixTensorTimesGradU = OpMixTensorTimesGradUImpl<SPACE_DIM, I, OpBase>;

  /**
   * @brief Integrate \f$(u_{i},\lambda_{ij,j})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM>
  using OpMixVecTimesDivLambda =
      OpMixVecTimesDivLambdaImpl<SPACE_DIM, I, OpBase>;

  /**
   * @brief Multiply vector times normal on the face times scalar function
   *
   * This operator typically will be used to evaluate natural boundary
   * conditions for mixed formulation.
   *
   * @tparam BASE_DIM
   * @tparam SPACE_DIM
   * @tparam OpBase
   */
  template <int SPACE_DIM>
  using OpNormalMixVecTimesScalar =
      OpNormalMixVecTimesScalarImpl<SPACE_DIM, I, OpBase>;

  /**
   * @brief Convective term
   *
   * \f[
   * (v, u_i \mathbf{y}_{,i})
   * \f]
   * where \f$\mathbf{y}\$ can be scalar or vector.
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   * @tparam SPACE_DIM
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM>
  using OpConvectiveTermRhs =
      OpConvectiveTermRhsImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I, OpBase>;
};

template <typename OpBase>
MoFEMErrorCode
OpSourceImpl<1, 1, GAUSS, SourceFunctionSpecialization::S<OpBase>>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha =
        t_w * vol * sourceFun(t_coords(0), t_coords(1), t_coords(2));
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows; ++rr) {
      OpBase::locF[rr] += alpha * t_row_base;
      ++t_row_base;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;
    ++t_coords;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, typename OpBase>
MoFEMErrorCode
OpSourceImpl<1, FIELD_DIM, GAUSS, SourceFunctionSpecialization::S<OpBase>>::
    iNtegrate(EntitiesFieldData::EntData &row_data) {
  FTensor::Index<'i', FIELD_DIM> i;
  MoFEMFunctionBegin;

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // source file
    auto t_source = sourceFun(t_coords(0), t_coords(1), t_coords(2));
    // take into account Jacobian
    const double alpha = t_w * vol;
    // loop over rows base functions
    auto t_nf = getFTensor1FromArray<FIELD_DIM, FIELD_DIM>(OpBase::locF);
    int rr = 0;
    for (; rr != OpBase::nbRows / FIELD_DIM; ++rr) {
      t_nf(i) += alpha * t_row_base * t_source(i);
      ++t_row_base;
      ++t_nf;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;
    ++t_coords;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, typename OpBase>
MoFEMErrorCode
OpSourceImpl<3, FIELD_DIM, GAUSS, SourceFunctionSpecialization::S<OpBase>>::
    iNtegrate(EntitiesFieldData::EntData &row_data) {
  FTensor::Index<'i', FIELD_DIM> i;
  MoFEMFunctionBegin;

  const size_t nb_base_functions = row_data.getN().size2() / 3;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor1N<3>();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // source file
    auto t_source = sourceFun(t_coords(0), t_coords(1), t_coords(2));
    // take into account Jacobian
    const double alpha = t_w * vol;
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows; ++rr) {
      OpBase::locF[rr] += alpha * t_row_base(i) * t_source(i);
      ++t_row_base;
    }
    for (; rr < nb_base_functions; ++rr)
      ++t_row_base;
    ++t_coords;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
}

template <int S, typename OpBase>
MoFEMErrorCode OpBaseTimesScalarImpl<1, S, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();
  // get vector values
  auto t_vec = getFTensor0FromVec<S>(*sourceVec);

#ifndef NDEBUG
  if (sourceVec->size() != OpBase::nbIntegrationPts) {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Wrong number of integration points %d != %d",
             OpBase::nbIntegrationPts, sourceVec->size());
  }
#endif

  // get coords
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha =
        t_w * vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows; ++rr) {
      OpBase::locF[rr] += alpha * t_row_base * t_vec;
      ++t_row_base;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;
    ++t_w; // move to another integration weight
    ++t_vec;
    ++t_coords;
  }
  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, int S, typename OpBase>
MoFEMErrorCode OpBaseTimesVectorImpl<1, FIELD_DIM, S, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();
  // get coords
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // get vector values
  auto t_vec = getFTensor1FromMat<FIELD_DIM, S>(*sourceVec);
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha =
        t_w * vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // get loc vector tensor
    auto t_nf = OpBase::template getNf<FIELD_DIM>();
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows / FIELD_DIM; ++rr) {
      t_nf(i) += alpha * t_row_base * t_vec(i);
      ++t_row_base;
      ++t_nf;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;
    ++t_w; // move to another integration weight
    ++t_vec;
    ++t_coords;
  }
  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, int S, typename OpBase>
MoFEMErrorCode OpBaseTimesVectorImpl<3, FIELD_DIM, S, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  const size_t nb_base_functions = row_data.getN().size2() / 3;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor1N<3>();
  // get coords
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // get vector values
  auto t_vec = getFTensor1FromMat<FIELD_DIM, S>(*sourceVec);
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha =
        t_w * vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows; ++rr) {
      OpBase::locF[rr] += alpha * t_row_base(i) * t_vec(i);
      ++t_row_base;
    }
    for (; rr < nb_base_functions; ++rr)
      ++t_row_base;
    ++t_w; // move to another integration weight
    ++t_vec;
    ++t_coords;
  }
  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, int S, typename OpBase>
MoFEMErrorCode
OpGradTimesTensorImpl<1, 1, SPACE_DIM, S, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_grad = row_data.getFTensor1DiffN<SPACE_DIM>();
  // get filed gradient values
  auto t_val_grad = getFTensor1FromMat<SPACE_DIM, S>(*(matVals));
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    const double beta = vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // take into account Jacobian
    const double alpha = t_w * beta;
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows; rr++) {
      // calculate element of local matrix
      OpBase::locF[rr] += alpha * (t_row_grad(i) * t_val_grad(i));
      ++t_row_grad; // move to another element of gradient of base
                    // function on row
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_grad;

    ++t_coords;
    ++t_val_grad;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, int S, typename OpBase>
MoFEMErrorCode
OpGradTimesTensorImpl<1, SPACE_DIM, SPACE_DIM, S, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_grad = row_data.getFTensor1DiffN<SPACE_DIM>();
  // get filed gradient values
  auto t_val_grad = getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(*(matVals));
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha =
        t_w * vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // get rhs vector
    auto t_nf = OpBase::template getNf<SPACE_DIM>();
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows / SPACE_DIM; rr++) {
      // calculate element of local matrix
      t_nf(i) += alpha * (t_row_grad(j) * t_val_grad(i, j));
      ++t_row_grad; // move to another element of gradient of base
      // function on row
      ++t_nf;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_grad;

    ++t_coords;
    ++t_val_grad;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, int S, typename OpBase>
MoFEMErrorCode
OpGradTimesSymTensorImpl<1, SPACE_DIM, SPACE_DIM, S, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get coords
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_grad = row_data.getFTensor1DiffN<SPACE_DIM>();
  // get filed gradient values
  auto t_val_mat = getFTensor2SymmetricFromMat<SPACE_DIM, S>(*(matVals));
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha =
        t_w * vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // get rhs vector
    auto t_nf = OpBase::template getNf<SPACE_DIM>();
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows / SPACE_DIM; rr++) {
      // calculate element of local matrix
      t_nf(j) += alpha * (t_row_grad(i) * t_val_mat(i, j));
      ++t_row_grad; // move to another element of gradient of base
                    // function on row
      ++t_nf;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_grad;
    ++t_val_mat;
    ++t_coords;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, int SPACE_DIM, typename OpBase,
          CoordinateTypes CoordSys>
MoFEMErrorCode
OpMixDivTimesUImpl<3, FIELD_DIM, SPACE_DIM, GAUSS, OpBase, CoordSys>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  const size_t nb_base_functions = row_data.getN().size2() / 3;
  auto t_w = this->getFTensor0IntegrationWeight();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  auto t_diff_base = row_data.getFTensor2DiffN<3, SPACE_DIM>();
  auto t_base = row_data.getFTensor1N<3>();
  auto t_u = getFTensor1FromMat<FIELD_DIM>(*(matVals));

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    const double alpha = this->getMeasure() * t_w *
                         betaConst(t_coords(0), t_coords(1), t_coords(2));
    auto t_nf = OpBase::template getNf<FIELD_DIM>();

    size_t bb = 0;
    for (; bb != this->nbRows / FIELD_DIM; ++bb) {
      const double t_div_base = t_diff_base(j, j);
      t_nf(i) += alpha * t_div_base * t_u(i);
      if constexpr (CoordSys == CYLINDRICAL) {
        t_nf(i) += alpha * (t_base(0) / t_coords(0)) * t_u(i);
      }
      ++t_nf;
      ++t_diff_base;
      ++t_base;
    }
    for (; bb < nb_base_functions; ++bb) {
      ++t_diff_base;
      ++t_base;
    }

    ++t_u;
    ++t_w;
    ++t_coords;
  }

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase, CoordinateTypes CoordSys>
MoFEMErrorCode
OpMixDivTimesUImpl<3, 1, SPACE_DIM, GAUSS, OpBase, CoordSys>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  const size_t nb_base_functions = row_data.getN().size2() / 3;
  auto t_w = this->getFTensor0IntegrationWeight();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  auto t_diff_base = row_data.getFTensor2DiffN<3, SPACE_DIM>();
  auto t_u = getFTensor0FromVec(*(vecVals));

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    const double alpha = this->getMeasure() * t_w *
                         betaConst(t_coords(0), t_coords(1), t_coords(2));
    ;

    size_t bb = 0;
    for (; bb != this->nbRows; ++bb) {
      const double t_div_base = t_diff_base(j, j);
      OpBase::locF[bb] += alpha * t_div_base * t_u;
      ++t_diff_base;
    }
    for (; bb < nb_base_functions; ++bb)
      ++t_diff_base;

    ++t_u;
    ++t_w;
    ++t_coords;
  }

  MoFEMFunctionReturn(0);
}

/**
 * @brief div U times vector
 *
 * \f[
 * \delta u_j = \phi^m\delta\overline{u}^m_j\\
 * \delta u_{j,i} = \phi^m_{,i}\delta\overline{u}^m_j\\
 * \textrm{tr}[\delta u_{j,i}] = \delta u_{j,i}\delta_{ji}\\
 * (\textrm{tr}[\delta u_{j,i}], v) =\\
 * (\delta u_{j,i} \delta_{ij}, v) =\\
 * (\delta u_{j,i}, \delta_{ij} v) =\\
 * (\phi^m_{,i}\delta\overline{u}^m_j, \delta_{ij} v) \\
 * f_i^m=(\phi^m_{,i}, v)
 * \f]
 *
 * @tparam FIELD_DIM
 * @tparam SPACE_DIM
 * @tparam OpBase
 * @param row_data
 * @return MoFEMErrorCode
 */
template <int FIELD_DIM, typename OpBase, CoordinateTypes CoordSys>
MoFEMErrorCode
OpMixDivTimesUImpl<1, FIELD_DIM, FIELD_DIM, GAUSS, OpBase, CoordSys>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  FTensor::Index<'i', FIELD_DIM> i;
  MoFEMFunctionBegin;

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // get base function gradient on rows
  auto t_row_grad = row_data.getFTensor1DiffN<FIELD_DIM>();
  // get vector values
  auto t_vec = getFTensor0FromVec(*sourceVec);
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha =
        t_w * vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    auto t_nf = getFTensor1FromArray<FIELD_DIM, FIELD_DIM>(OpBase::locF);
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows / FIELD_DIM; ++rr) {
      t_nf(i) += alpha * t_row_grad(i) * t_vec;
      ++t_row_grad;
      ++t_nf;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_grad;
    ++t_w; // move to another integration weight
    ++t_vec;
    ++t_coords;
  }
  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode OpMixTensorTimesGradUImpl<SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  const size_t nb_base_functions = row_data.getN().size2() / 3;
  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_coords = this->getFTensor1CoordsAtGaussPts();
  auto t_base = row_data.getFTensor1N<3>();
  auto t_grad = getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(*(matVals));

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    const double alpha = this->getMeasure() * t_w;
    auto t_nf = OpBase::template getNf<SPACE_DIM>();

    size_t bb = 0;
    for (; bb != this->nbRows / SPACE_DIM; ++bb) {
      t_nf(i) += alpha * betaCoeff(t_coords(0), t_coords(1), t_coords(2)) *
                 t_base(j) * t_grad(i, j);
      ++t_nf;
      ++t_base;
    }
    for (; bb < nb_base_functions; ++bb)
      ++t_base;

    ++t_grad;
    ++t_coords;
    ++t_w;
  }

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode OpMixVecTimesDivLambdaImpl<SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  const size_t nb_base_functions = row_data.getN().size2();
  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_coords = this->getFTensor1CoordsAtGaussPts();
  auto t_base = row_data.getFTensor0N();
  auto t_div = getFTensor1FromMat<SPACE_DIM>(*(matVals));

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {
    const double alpha = this->getMeasure() * t_w;
    auto t_nf = OpBase::template getNf<SPACE_DIM>();

    size_t bb = 0;
    for (; bb != this->nbRows / SPACE_DIM; ++bb) {
      t_nf(i) += alpha * t_base *
                 betaCoeff(t_coords(0), t_coords(1), t_coords(2)) * t_div(i);
      ++t_nf;
      ++t_base;
    }
    for (; bb < nb_base_functions; ++bb)
      ++t_base;

    ++t_coords;
    ++t_div;
    ++t_w;
  }

  MoFEMFunctionReturn(0);
}

template <typename OpBase>
MoFEMErrorCode OpNormalMixVecTimesScalarImpl<3, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  const size_t nb_base_functions = row_data.getN().size2() / 3;
  // get element volume
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor1N<3>();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // get normal
  auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha = t_w * sourceFun(t_coords(0), t_coords(1), t_coords(2));
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows; ++rr) {
      OpBase::locF[rr] += alpha * t_row_base(i) * t_normal(i);
      ++t_row_base;
    }
    for (; rr < nb_base_functions; ++rr)
      ++t_row_base;
    ++t_coords;
    ++t_w; // move to another integration weight
    ++t_normal;
  }
  MoFEMFunctionReturn(0);
}

template <typename OpBase>
MoFEMErrorCode OpNormalMixVecTimesScalarImpl<2, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  const size_t nb_base_functions = row_data.getN().size2() / 3;
  FTensor::Tensor1<double, 3> t_z{0., 0., 1.};
  // get element volume
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor1N<3>();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // get normal
  auto t_tangent = OpBase::getFTensor1TangentAtGaussPts();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha = t_w * sourceFun(t_coords(0), t_coords(1), t_coords(2));
    FTensor::Tensor1<double, 3> t_normal;
    FTensor::Index<'j', 3> j;
    FTensor::Index<'k', 3> k;
    t_normal(i) = FTensor::levi_civita(i, j, k) * t_tangent(j) * t_z(k);
    int rr = 0;
    for (; rr != OpBase::nbRows; ++rr) {
      OpBase::locF[rr] += alpha * t_row_base(i) * t_normal(i);
      ++t_row_base;
    }
    for (; rr < nb_base_functions; ++rr)
      ++t_row_base;
    ++t_coords;
    ++t_tangent;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode
OpConvectiveTermRhsImpl<1, 1, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_base = row_data.getFTensor0N();

  auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
  auto t_grad_y = getFTensor1FromMat<SPACE_DIM>(*yGradPtr);

  FTensor::Index<'i', SPACE_DIM> i;
  const double alpha_constant = alphaConstant();
  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    // get element volume
    const double vol = OpBase::getMeasure();
    const double c = (t_grad_y(i) * t_u(i)) * (t_w * vol * alpha_constant);

    // get element volume
    int rr = 0;
    for (; rr != OpBase::nbRows; ++rr) {
      OpBase::locF[rr] += c * t_base;
      ++t_base;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_base;

    ++t_w; // move to another integration weight
    ++t_u;
    ++t_grad_y;
  }

  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, int SPACE_DIM, typename OpBase>
MoFEMErrorCode
OpConvectiveTermRhsImpl<1, FIELD_DIM, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_base = row_data.getFTensor0N();

  auto t_u = getFTensor1FromMat<SPACE_DIM>(*uPtr);
  auto t_grad_y = getFTensor2FromMat<FIELD_DIM, SPACE_DIM>(*yGradPtr);

  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'J', FIELD_DIM> J;
  const double alpha_constant = alphaConstant();
  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    // get element volume
    const double vol = OpBase::getMeasure();

    FTensor::Tensor1<double, FIELD_DIM> t_c;
    t_c(J) = (t_grad_y(J, i) * t_u(i)) * (t_w * vol * alpha_constant);

    auto t_nf = getFTensor1FromArray<FIELD_DIM, FIELD_DIM>(OpBase::locF);
    int rr = 0;
    for (; rr != OpBase::nbRows / FIELD_DIM; ++rr) {
      t_nf(J) += t_base * t_c(J);
      ++t_base;
      ++t_nf;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_base;
    ++t_w; // move to another integration weight

    ++t_u;
    ++t_grad_y;
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM

#endif // __LINEAR_FORMS_INTEGRATORS_HPP__