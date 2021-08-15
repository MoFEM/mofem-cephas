/** \file LinearFormsIntegrators.hpp
  * \brief Linear forms integrators
  * \ingroup mofem_form

*/

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#ifndef __LINEAR_FORMS_INTEGRATORS_HPP__
#define __LINEAR_FORMS_INTEGRATORS_HPP__

namespace MoFEM {

template <int BASE_DIM, int FIELD_DIM, IntegrationType I, typename OpBase>
struct OpSourceImpl;

/**
 * @brief Integrate source  
 * 
 * @tparam OpBase 
 */
template <typename OpBase>
struct OpSourceImpl<1, 1, GAUSS, OpBase> : public OpBase {

  /**
   * @brief Construct a new Op Source Impl object
   * 
   * @param field_name 
   * @param source_fun 
   * @param ents_ptr 
   */
  OpSourceImpl(const std::string field_name, ScalarFun source_fun,
               boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW), sourceFun(source_fun),
        entsPtr(ents_ptr) {}

protected:
  boost::shared_ptr<Range> entsPtr;
  ScalarFun sourceFun;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
};

template <int FIELD_DIM, typename OpBase>
struct OpSourceImpl<1, FIELD_DIM, GAUSS, OpBase> : public OpBase {
  OpSourceImpl(const std::string field_name, VectorFun<FIELD_DIM> source_fun)
      : OpBase(field_name, field_name, OpBase::OPROW), sourceFun(source_fun) {}

protected:
  VectorFun<FIELD_DIM> sourceFun;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
};

template <int BASE_DIM, typename OpBase>
struct OpSourceImpl<BASE_DIM, BASE_DIM, GAUSS, OpBase> : public OpBase {
  OpSourceImpl(const std::string field_name, VectorFun<BASE_DIM> source_fun)
      : OpBase(field_name, field_name, OpBase::OPROW), sourceFun(source_fun) {}

protected:
  VectorFun<BASE_DIM> sourceFun;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
};

template <int BASE_DIM, int S, IntegrationType I, typename OpBase>
struct OpBaseTimesScalarFieldImpl;

template <int S, typename OpBase>
struct OpBaseTimesScalarFieldImpl<1, S, GAUSS, OpBase> : public OpBase {

  OpBaseTimesScalarFieldImpl(const std::string field_name,
                             boost::shared_ptr<VectorDouble> vec,
                             const double beta_coeff)
      : OpBase(field_name, field_name, OpBase::OPROW), sourceVec(vec),
        betaCoeff(beta_coeff) {}

protected:
  const double betaCoeff;
  boost::shared_ptr<VectorDouble> sourceVec;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
};

template <int BASE_DIM, int FIELD_DIM, int S, IntegrationType I,
          typename OpBase>
struct OpBaseTimesVectorImpl;

template <int FIELD_DIM, int S, typename OpBase>
struct OpBaseTimesVectorImpl<1, FIELD_DIM, S, GAUSS, OpBase> : public OpBase {

  OpBaseTimesVectorImpl(const std::string field_name,
                        boost::shared_ptr<MatrixDouble> vec,
                        ScalarFun beta_coeff,
                        boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW), sourceVec(vec),
        betaCoeff(beta_coeff), entsPtr(ents_ptr) {}

protected:
  ScalarFun betaCoeff;
  boost::shared_ptr<Range> entsPtr;
  boost::shared_ptr<MatrixDouble> sourceVec;
  FTensor::Index<'i', FIELD_DIM> i;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
};

template <int BASE_DIM, int S, typename OpBase>
struct OpBaseTimesVectorImpl<BASE_DIM, BASE_DIM, S, GAUSS, OpBase>
    : public OpBase {

  OpBaseTimesVectorImpl(const std::string field_name,
                        boost::shared_ptr<MatrixDouble> vec,
                        ScalarFun beta_coeff,
                        boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(field_name, field_name, OpBase::OPROW), sourceVec(vec),
        betaCoeff(beta_coeff), entsPtr(ents_ptr) {}

protected:
  ScalarFun betaCoeff;
  boost::shared_ptr<Range> entsPtr;
  boost::shared_ptr<MatrixDouble> sourceVec;
  FTensor::Index<'i', BASE_DIM> i;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
};

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S, IntegrationType I,
          typename OpBase>
struct OpGradTimesTensorImpl;

template <int SPACE_DIM, int S, typename OpBase>
struct OpGradTimesTensorImpl<1, 1, SPACE_DIM, S, GAUSS, OpBase>
    : public OpBase {

  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index

  OpGradTimesTensorImpl(const std::string field_name,
                        boost::shared_ptr<MatrixDouble> mat_vals,
                        ScalarFun beta_coeff)
      : OpBase(field_name, field_name, OpBase::OPROW), matVals(mat_vals),
        betaCoeff(beta_coeff) {}

protected:
  boost::shared_ptr<MatrixDouble> matVals;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
  ScalarFun betaCoeff;
};

template <int SPACE_DIM, int S, typename OpBase>
struct OpGradTimesTensorImpl<1, SPACE_DIM, SPACE_DIM, S, GAUSS, OpBase>
    : public OpBase {

  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  FTensor::Index<'j', SPACE_DIM> j; ///< summit Index

  OpGradTimesTensorImpl(
      const std::string field_name, boost::shared_ptr<MatrixDouble> mat_vals,
      ScalarFun beta_coeff = [](double, double, double) { return 1; })
      : OpBase(field_name, field_name, OpBase::OPROW), matVals(mat_vals),
        betaCoeff(beta_coeff) {}

protected:
  boost::shared_ptr<MatrixDouble> matVals;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
  ScalarFun betaCoeff;
};

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S, IntegrationType I,
          typename OpBase>
struct OpGradTimesSymTensorImpl {};

template <int SPACE_DIM, int S, typename OpBase>
struct OpGradTimesSymTensorImpl<1, SPACE_DIM, SPACE_DIM, S, GAUSS, OpBase>
    : public OpBase {

  OpGradTimesSymTensorImpl(
      const std::string field_name, boost::shared_ptr<MatrixDouble> &mat_vals,
      ScalarFun beta_coeff = [](double, double, double) { return 1; })
      : OpBase(field_name, field_name, OpBase::OPROW), matVals(mat_vals),
        betaCoeff(beta_coeff) {}

protected:
  boost::shared_ptr<MatrixDouble> matVals;
  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'j', SPACE_DIM> j;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
  ScalarFun betaCoeff;
};

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, IntegrationType I,
          typename OpBase>
struct OpMixDivTimesUImpl {};

template <int FIELD_DIM, int SPACE_DIM, typename OpBase>
struct OpMixDivTimesUImpl<3, FIELD_DIM, SPACE_DIM, GAUSS, OpBase>
    : public OpBase {
  OpMixDivTimesUImpl(const std::string field_name,
                     boost::shared_ptr<MatrixDouble> mat_vals,
                     const double beta)
      : OpBase(field_name, field_name, OpBase::OPROW), matVals(mat_vals),
        betaConst(beta) {}

protected:
  const double betaConst;
  boost::shared_ptr<MatrixDouble> matVals;
  FTensor::Index<'i', FIELD_DIM> i;
  FTensor::Index<'j', SPACE_DIM> j;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
};

template <int SPACE_DIM, typename OpBase>
struct OpMixDivTimesUImpl<3, 1, SPACE_DIM, GAUSS, OpBase> : public OpBase {
  OpMixDivTimesUImpl(const std::string field_name,
                     boost::shared_ptr<VectorDouble> vec_vals,
                     const double beta)
      : OpBase(field_name, field_name, OpBase::OPROW), vecVals(vec_vals),
        betaConst(beta) {}

protected:
  const double betaConst;
  boost::shared_ptr<VectorDouble> vecVals;
  FTensor::Index<'j', SPACE_DIM> j;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
};

template <int FIELD_DIM, typename OpBase>
struct OpMixDivTimesUImpl<1, FIELD_DIM, FIELD_DIM, GAUSS, OpBase>
    : public OpBase {

  OpMixDivTimesUImpl(const std::string field_name,
                             boost::shared_ptr<VectorDouble> vec,
                             const double beta_coeff)
      : OpBase(field_name, field_name, OpBase::OPROW), sourceVec(vec),
        betaCoeff(beta_coeff) {}

protected:
  const double betaCoeff;
  boost::shared_ptr<VectorDouble> sourceVec;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
};

template <int SPACE_DIM, IntegrationType I, typename OpBase>
struct OpMixVecTimesDivLambdaImpl {};

template <int SPACE_DIM, typename OpBase>
struct OpMixVecTimesDivLambdaImpl<SPACE_DIM, GAUSS, OpBase> : public OpBase {
  OpMixVecTimesDivLambdaImpl(const std::string field_name,
                             boost::shared_ptr<MatrixDouble> &mat_vals)
      : OpBase(field_name, field_name, OpBase::OPROW), matVals(mat_vals) {}

protected:
  boost::shared_ptr<MatrixDouble> matVals;
  FTensor::Index<'i', SPACE_DIM> i;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
};

template <int SPACE_DIM, IntegrationType I, typename OpBase>
struct OpMixTensorTimesGradUImpl {};

template <int SPACE_DIM, typename OpBase>
struct OpMixTensorTimesGradUImpl<SPACE_DIM, GAUSS, OpBase> : public OpBase {
  OpMixTensorTimesGradUImpl(const std::string field_name,
                            boost::shared_ptr<MatrixDouble> &mat_vals)
      : OpBase(field_name, field_name, OpBase::OPROW), matVals(mat_vals) {}

protected:
  boost::shared_ptr<MatrixDouble> matVals;
  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'j', SPACE_DIM> j;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
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
  template <int BASE_DIM, int FIELD_DIM>
  struct OpSource : public OpSourceImpl<BASE_DIM, FIELD_DIM, I, OpBase> {
    using OpSourceImpl<BASE_DIM, FIELD_DIM, I, OpBase>::OpSourceImpl;
  };

  /**
   * @brief Vector field integrator \f$(v_i,f)_\Omega\f$, f is a vector
   * @ingroup mofem_forms
   *
   * @tparam BASE_DIM
   */
  template <int BASE_DIM, int S = 1>
  struct OpBaseTimesScalarField
      : public OpBaseTimesScalarFieldImpl<BASE_DIM, S, I, OpBase> {
    using OpBaseTimesScalarFieldImpl<BASE_DIM, S, I,
                                     OpBase>::OpBaseTimesScalarFieldImpl;
  };


  /**
   * @brief Vector field integrator \f$(v,f_i)_\Omega\f$, f is a vector
   * @ingroup mofem_forms
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   * @tparam 0
   */
  template <int BASE_DIM, int FIELD_DIM, int S>
  struct OpBaseTimesVector
      : public OpBaseTimesVectorImpl<BASE_DIM, FIELD_DIM, S, I, OpBase> {
    using OpBaseTimesVectorImpl<BASE_DIM, FIELD_DIM, S, I,
                                OpBase>::OpBaseTimesVectorImpl;
  };

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
  struct OpGradTimesTensor
      : public OpGradTimesTensorImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I,
                                     OpBase> {
    using OpGradTimesTensorImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I,
                                OpBase>::OpGradTimesTensorImpl;
  };

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
  struct OpGradTimesSymTensor
      : public OpGradTimesSymTensorImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I,
                                        OpBase> {
    using OpGradTimesSymTensorImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I,
                                   OpBase>::OpGradTimesSymTensorImpl;
  };

  /**
   * @brief Integrate \f$(\lambda_{ij,j},u_{i})_\Omega\f$
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   * @tparam SPACE_DIM
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM>
  struct OpMixDivTimesU
      : public OpMixDivTimesUImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I, OpBase> {
    using OpMixDivTimesUImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I,
                             OpBase>::OpMixDivTimesUImpl;
  };

  /**
   * @brief Integrate \f$(\lambda_{ij},u_{i,j})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM>
  struct OpMixTensorTimesGradU
      : public OpMixTensorTimesGradUImpl<SPACE_DIM, I, OpBase> {
    using OpMixTensorTimesGradUImpl<SPACE_DIM, I,
                                    OpBase>::OpMixTensorTimesGradUImpl;
  };

  /**
   * @brief Integrate \f$(u_{i},\lambda_{ij,j})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM>
  struct OpMixVecTimesDivLambda
      : public OpMixVecTimesDivLambdaImpl<SPACE_DIM, I, OpBase> {
    using OpMixVecTimesDivLambdaImpl<SPACE_DIM, I,
                                     OpBase>::OpMixVecTimesDivLambdaImpl;
  };
};

template <typename OpBase>
MoFEMErrorCode OpSourceImpl<1, 1, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
  MoFEMFunctionBegin;

  if (entsPtr) {
    if (entsPtr->find(OpBase::getFEEntityHandle()) == entsPtr->end())
      MoFEMFunctionReturnHot(0);
  }
  
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
MoFEMErrorCode OpSourceImpl<1, FIELD_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
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

template <int BASE_DIM, typename OpBase>
MoFEMErrorCode OpSourceImpl<BASE_DIM, BASE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
  FTensor::Index<'i', BASE_DIM> i;
  MoFEMFunctionBegin;
  const size_t nb_base_functions = row_data.getN().size2() / BASE_DIM;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor1N<BASE_DIM>();
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
MoFEMErrorCode OpBaseTimesScalarFieldImpl<1, S, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
  MoFEMFunctionBegin;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();
  // get vector values
  auto t_vec = getFTensor0FromVec<S>(*sourceVec);
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha = t_w * vol * betaCoeff;
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
  }
  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, int S, typename OpBase>
MoFEMErrorCode OpBaseTimesVectorImpl<1, FIELD_DIM, S, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
  MoFEMFunctionBegin;
  if (entsPtr) {
    if (entsPtr->find(OpBase::getFEEntityHandle()) == entsPtr->end())
      MoFEMFunctionReturnHot(0);
  }
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

template <int BASE_DIM, int S, typename OpBase>
MoFEMErrorCode
OpBaseTimesVectorImpl<BASE_DIM, BASE_DIM, S, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
  MoFEMFunctionBegin;
  if (entsPtr) {
    if (entsPtr->find(OpBase::getFEEntityHandle()) == entsPtr->end())
      MoFEMFunctionReturnHot(0);
  }
  const size_t nb_base_functions = row_data.getN().size2() / BASE_DIM;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor1N<BASE_DIM>();
  // get coords
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // get vector values
  auto t_vec = getFTensor1FromMat<BASE_DIM, S>(*sourceVec);
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
    DataForcesAndSourcesCore::EntData &row_data) {
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
    DataForcesAndSourcesCore::EntData &row_data) {
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
    DataForcesAndSourcesCore::EntData &row_data) {
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
    const double alpha = t_w * vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
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

template <int FIELD_DIM, int SPACE_DIM, typename OpBase>
MoFEMErrorCode
OpMixDivTimesUImpl<3, FIELD_DIM, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
  MoFEMFunctionBegin;

  const size_t nb_base_functions = row_data.getN().size2() / 3;
  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_diff_base = row_data.getFTensor2DiffN<3, SPACE_DIM>();
  auto t_u = getFTensor1FromMat<FIELD_DIM>(*(matVals));

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    const double alpha = this->getMeasure() * t_w * betaConst;
    auto t_nf = OpBase::template getNf<FIELD_DIM>();

    size_t bb = 0;
    for (; bb != this->nbRows / FIELD_DIM; ++bb) {
      const double t_div_base = t_diff_base(j, j);
      t_nf(i) += alpha * t_div_base * t_u(i);
      ++t_nf;
      ++t_diff_base;
    }
    for (; bb < nb_base_functions; ++bb)
      ++t_diff_base;

    ++t_u;
    ++t_w;
  }

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode OpMixDivTimesUImpl<3, 1, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
  MoFEMFunctionBegin;

  const size_t nb_base_functions = row_data.getN().size2() / 3;
  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_diff_base = row_data.getFTensor2DiffN<3, SPACE_DIM>();
  auto t_u = getFTensor0FromVec(*(vecVals));

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    const double alpha = this->getMeasure() * t_w * betaConst;

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
template <int FIELD_DIM, typename OpBase>
MoFEMErrorCode
OpMixDivTimesUImpl<1, FIELD_DIM, FIELD_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
  FTensor::Index<'i', FIELD_DIM> i;
  MoFEMFunctionBegin;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_grad = row_data.getFTensor1DiffN<FIELD_DIM>();
  // get vector values
  auto t_vec = getFTensor0FromVec(*sourceVec);
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobian
    const double alpha = t_w * vol * betaCoeff;
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
  }
  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode OpMixTensorTimesGradUImpl<SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
  MoFEMFunctionBegin;

  const size_t nb_base_functions = row_data.getN().size2() / 3;
  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_base = row_data.getFTensor1N<3>();
  auto t_grad = getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(*(matVals));

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    const double alpha = this->getMeasure() * t_w;
    auto t_nf = OpBase::template getNf<SPACE_DIM>();

    size_t bb = 0;
    for (; bb != this->nbRows / SPACE_DIM; ++bb) {
      t_nf(i) += alpha * t_base(j) * t_grad(i, j);
      ++t_nf;
      ++t_base;
    }
    for (; bb < nb_base_functions; ++bb)
      ++t_base;

    ++t_grad;
    ++t_w;
  }

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode OpMixVecTimesDivLambdaImpl<SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
  MoFEMFunctionBegin;

  const size_t nb_base_functions = row_data.getN().size2();
  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_base = row_data.getFTensor0N();
  auto t_div = getFTensor1FromMat<SPACE_DIM>(*(matVals));

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {
    const double alpha = this->getMeasure() * t_w;
    auto t_nf = OpBase::template getNf<SPACE_DIM>();

    size_t bb = 0;
    for (; bb != this->nbRows / SPACE_DIM; ++bb) {
      t_nf(i) += alpha * t_base * t_div(i);
      ++t_nf;
      ++t_base;
    }
    for (; bb < nb_base_functions; ++bb)
      ++t_base;

    ++t_div;
    ++t_w;
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM

#endif // __LINEAR_FORMS_INTEGRATORS_HPP__