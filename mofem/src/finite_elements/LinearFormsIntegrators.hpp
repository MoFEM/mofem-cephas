/** \file LinearFormsIntegrators.hpp
  * \brief Linear forms inteegrators
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

template <typename OpBase>
struct OpSourceImpl<1, 1, GAUSS, OpBase> : public OpBase {
  OpSourceImpl(const std::string field_name, ScalarFun source_fun)
      : OpBase(field_name, field_name, OpBase::OPROW), sourceFun(source_fun) {}

protected:
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

template <int BASE_DIM, int FIELD_DIM, int S, IntegrationType I,
          typename OpBase>
struct OpBaseTimesVectorImpl;

template <int FIELD_DIM, int S, typename OpBase>
struct OpBaseTimesVectorImpl<1, FIELD_DIM, S, GAUSS, OpBase>
    : public OpBase {

  OpBaseTimesVectorImpl(const std::string field_name,
                        boost::shared_ptr<MatrixDouble> vec)
      : OpBase(field_name, field_name, OpBase::OPROW), sourceVec(vec) {}

protected:
  boost::shared_ptr<MatrixDouble> sourceVec;
  FTensor::Index<'i', FIELD_DIM> i;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
};

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S, IntegrationType I,
          typename OpBase>
struct OpGradTimesTensorImpl;

template <int SPACE_DIM, typename OpBase>
struct OpGradTimesTensorImpl<1, 1, SPACE_DIM, 1, GAUSS, OpBase>
    : public OpBase {

  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index

  OpGradTimesTensorImpl(const std::string field_name,
                        boost::shared_ptr<MatrixDouble> mat_vals)
      : OpBase(field_name, field_name, OpBase::OPROW), matVals(mat_vals) {}

protected:
  boost::shared_ptr<MatrixDouble> matVals;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data);
};

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S, IntegrationType I,
          typename OpBase>
struct OpGradTimesSymTensorImpl {};

template <int SPACE_DIM, int S, typename OpBase>
struct OpGradTimesSymTensorImpl<1, SPACE_DIM, SPACE_DIM, S, GAUSS, OpBase>
    : public OpBase {

  OpGradTimesSymTensorImpl(const std::string field_name,
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
   * @brief Vector field integrator \f$(v,f_i)_\Omega\f$, f is a vector
   * @ingroup mofem_forms
   *
   * @note \f$f(\mathbf{x})\$ is scalar function of coordinates
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   * @tparam 0
   */
  template <int BASE_DIM, int FIELD_DIM, int S = 0>
  struct OpBaseTimesVector
      : public OpBaseTimesVectorImpl<BASE_DIM, FIELD_DIM, S, I,
                                     OpBase> {
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
   * @brief Integrate \f$(v_{i},f_{ij})_\Omega\f$, f is symmetric tensor
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
      : public OpGradTimesSymTensorImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S,
                                              I, OpBase> {
    using OpGradTimesSymTensorImpl<
        BASE_DIM, FIELD_DIM, SPACE_DIM, S, I,
        OpBase>::OpGradTimesSymTensorImpl;
  };

};

template <typename OpBase>
MoFEMErrorCode OpSourceImpl<1, 1, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
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
    // take into account Jacobean
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
    // take into account Jacobean
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

template <int FIELD_DIM, int S, typename OpBase>
MoFEMErrorCode OpBaseTimesVectorImpl<1, FIELD_DIM, S, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
  MoFEMFunctionBegin;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();
  // get vector values
  auto t_vec = getFTensor1FromMat<FIELD_DIM, S>(*sourceVec);
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobean
    const double alpha = t_w * vol;
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
  }
  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode
OpGradTimesTensorImpl<1, 1, SPACE_DIM, 1, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data) {
  MoFEMFunctionBegin;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_grad = row_data.getFTensor1DiffN<SPACE_DIM>();
  // get filed gradient values
  auto t_val_grad = getFTensor1FromMat<SPACE_DIM>(*(matVals));
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobean
    const double alpha = t_w * vol;
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

    ++t_val_grad;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, int S, typename OpBase>
MoFEMErrorCode
OpGradTimesSymTensorImpl<1, SPACE_DIM, SPACE_DIM, S, GAUSS, OpBase>::
    iNtegrate(DataForcesAndSourcesCore::EntData &row_data) {
  MoFEMFunctionBegin;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_grad = row_data.getFTensor1DiffN<SPACE_DIM>();
  // get filed gradient values
  auto t_val_mat = getFTensor2SymmetricFromMat<SPACE_DIM, S>(*(matVals));
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    // take into account Jacobean
    const double alpha = t_w * vol;
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
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM

#endif // __LINEAR_FORMS_INTEGRATORS_HPP__