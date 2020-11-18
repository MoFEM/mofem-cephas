/** \file BiLinearFormsIntegrators.hpp
  * \brief Bilinear forms integrators
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

#ifndef __BILINEAR_FORMS_INTEGRATORS_HPP__
#define __BILINEAR_FORMS_INTEGRATORS_HPP__

namespace MoFEM {

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, IntegrationType I,
          typename OpBase>
struct OpGradGradImpl {};

template <int SPACE_DIM, typename OpBase>
struct OpGradGradImpl<1, 1, SPACE_DIM, GAUSS, OpBase> : public OpBase {
  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  OpGradGradImpl(const std::string row_field_name,
                 const std::string col_field_name, ScalarFun beta)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        betaCoeff(beta) {
    if(row_field_name == col_field_name)
      this->sYmm = true;
  }

protected:
  ScalarFun betaCoeff;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int BASE_DIM, int FIELD_DIM, IntegrationType I, typename OpBase>
struct OpMassImpl {};

template <typename OpBase>
struct OpMassImpl<1, 1, GAUSS, OpBase> : public OpBase {

  OpMassImpl(const std::string row_field_name, const std::string col_field_name,
             ScalarFun beta)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        betaCoeff(beta) {
    if(row_field_name == col_field_name)
      this->sYmm = true;
  }

protected:
  ScalarFun betaCoeff;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int FIELD_DIM, typename OpBase>
struct OpMassImpl<1, FIELD_DIM, GAUSS, OpBase>
    : public OpMassImpl<1, 1, GAUSS, OpBase> {
  using OpMassImpl<1, 1, GAUSS, OpBase>::OpMassImpl;

protected:
  using OpMassImpl<1, 1, GAUSS, OpBase>::betaCoeff;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S, IntegrationType I,
          typename OpBase>
struct OpGradSymTensorGradImpl {};

template <int SPACE_DIM, int S, typename OpBase>
struct OpGradSymTensorGradImpl<1, SPACE_DIM, SPACE_DIM, S, GAUSS, OpBase>
    : public OpBase {
  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  OpGradSymTensorGradImpl(const std::string row_field_name,
                          const std::string col_field_name,
                          boost::shared_ptr<MatrixDouble> mat_D)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL), matD(mat_D) {
    if (row_field_name == col_field_name)
      this->sYmm = true;
  }

protected:
  boost::shared_ptr<MatrixDouble> matD;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S, IntegrationType I,
          typename OpBase>
struct OpGradTensorGradImpl {};

template <int SPACE_DIM, int S, typename OpBase>
struct OpGradTensorGradImpl<1, SPACE_DIM, SPACE_DIM, S, GAUSS, OpBase>
    : public OpBase {
  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  OpGradTensorGradImpl(const std::string row_field_name,
                       const std::string col_field_name,
                       boost::shared_ptr<MatrixDouble> mat_D)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL), matD(mat_D) {}

protected:
  boost::shared_ptr<MatrixDouble> matD;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int SPACE_DIM, IntegrationType I, typename OpBase>
struct OpMixDivLambdaTimesUImpl {};

template <int SPACE_DIM, typename OpBase>
struct OpMixDivLambdaTimesUImpl<SPACE_DIM, GAUSS, OpBase> : public OpBase {
  OpMixDivLambdaTimesUImpl(const std::string row_field_name,
                            const std::string col_field_name,
                            const double alpha = 1,
                            bool assemble_transpose = false)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        alphaConstant(alpha) {
    this->assembleTranspose = assemble_transpose;
  }

protected:
  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  const double alphaConstant;
  boost::shared_ptr<MatrixDouble> matLoc;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int SPACE_DIM, IntegrationType I, typename OpBase>
struct OpMixLambdaTimesGradUImpl {};

template <int SPACE_DIM, typename OpBase>
struct OpMixLambdaTimesGradUImpl<SPACE_DIM, GAUSS, OpBase> : public OpBase {
  OpMixLambdaTimesGradUImpl(const std::string row_field_name,
                            const std::string col_field_name,
                            const double alpha = 1,
                            bool assemble_transpose = false)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        alphaConstant(alpha) {
    this->assembleTranspose = assemble_transpose;
  }

protected:
  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  FTensor::Index<'j', SPACE_DIM> j; ///< summit Index
  const double alphaConstant;
  boost::shared_ptr<MatrixDouble> matLoc;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

/**
 * @brief Bilinear integrator form
 * @ingroup mofem_forms
 *
 * @tparam EleOp
 * @tparam A
 * @tparam I
 */
template <typename EleOp>
template <AssemblyType A>
template <IntegrationType I>
struct FormsIntegrators<EleOp>::Assembly<A>::BiLinearForm {

  /**
   * @brief Integrate \f$(v_{,i},\beta(\mathbf{x}) u_{,j}))_\Omega\f$
   * @ingroup mofem_forms
   *
   * @tparam SPACE_DIM
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM>
  struct OpGradGrad
      : public OpGradGradImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I, OpBase> {
    using OpGradGradImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I,
                         OpBase>::OpGradGradImpl;
  };

  /**
   * @brief Integrate \f$(v_i,\beta(\mathbf{x}) u_j)_\Omega\f$
   * @ingroup mofem_forms
   *
   * @tparam
   */
  template <int BASE_DIM, int FIELD_DIM>
  struct OpMass : public OpMassImpl<BASE_DIM, FIELD_DIM, I, OpBase> {
    using OpMassImpl<BASE_DIM, FIELD_DIM, I, OpBase>::OpMassImpl;
  };

  /**
   * @brief Integrate \f$(v_k,D_{ijkl} u_{,l})_\Omega\f$
   *
   * \note \f$D_{ijkl}\f$ is * tensor with minor & major symmetry
   *
   * @ingroup mofem_forms
   *
   * @tparam SPACE_DIM
   * @tparam S
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S = 0>
  struct OpGradSymTensorGrad
      : public OpGradSymTensorGradImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I,
                                       OpBase> {
    using OpGradSymTensorGradImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I,
                                  OpBase>::OpGradSymTensorGradImpl;
  };

  /**
   * @brief Integrate \f$(v_k,D_{ijkl} u_{,l})_\Omega\f$
   *
   * \note \f$D_{ijkl}\f$ is * tensor with no symmetries
   *
   * @ingroup mofem_forms
   *
   * @tparam SPACE_DIM
   * @tparam S
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S = 0>
  struct OpGradTensorGrad
      : public OpGradTensorGradImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I,
                                    OpBase> {
    using OpGradTensorGradImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I,
                               OpBase>::OpGradTensorGradImpl;
  };

  template <int SPACE_DIM>
  struct OpMixDivLambdaTimesU
      : public OpMixDivLambdaTimesUImpl<SPACE_DIM, I, OpBase> {
    using OpMixDivLambdaTimesUImpl<SPACE_DIM, I,
                                   OpBase>::OpMixDivLambdaTimesUImpl;
  };

  template <int SPACE_DIM>
  struct OpMixLambdaTimesGradU
      : public OpMixLambdaTimesGradUImpl<SPACE_DIM, I, OpBase> {
    using OpMixLambdaTimesGradUImpl<SPACE_DIM, I,
                                    OpBase>::OpMixLambdaTimesGradUImpl;
  };
};

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode OpGradGradImpl<1, 1, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_grad = row_data.getFTensor1DiffN<SPACE_DIM>();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();

  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    const double beta = vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // take into account Jacobean
    const double alpha = t_w * beta;
    // loop over ros base functions
    int rr = 0;
    for (; rr != OpBase::nbRows; rr++) {
      // get column base functions gradient at gauss point gg
      auto t_col_grad = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);
      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols; cc++) {
        // calculate element of local matrix
        OpBase::locMat(rr, cc) += alpha * (t_row_grad(i) * t_col_grad(i));
        ++t_col_grad; // move to another gradient of base function
                      // on column
      }
      ++t_row_grad; // move to another element of gradient of base
                    // function on row
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_grad;

    ++t_coords;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
}

template <typename OpBase>
MoFEMErrorCode OpMassImpl<1, 1, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
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
    const double beta = vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // take into account Jacobean
    const double alpha = t_w * beta;
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows; rr++) {
      // get column base functions gradient at gauss point gg
      auto t_col_base = col_data.getFTensor0N(gg, 0);
      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols; cc++) {
        // calculate element of local matrix
        OpBase::locMat(rr, cc) += alpha * (t_row_base * t_col_base);
        ++t_col_base;
      }
      ++t_row_base;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;
    ++t_coords;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
};

template <int FIELD_DIM, typename OpBase>
MoFEMErrorCode OpMassImpl<1, FIELD_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();

  FTensor::Index<'i', FIELD_DIM> i;
  auto get_t_vec = [&](const int rr) {
    std::array<double *, FIELD_DIM> ptrs;
    for (auto i = 0; i != FIELD_DIM; ++i)
      ptrs[i] = &OpBase::locMat(rr + i, i);
    FTensor::Tensor1<FTensor::PackPtr<double *, FIELD_DIM>, FIELD_DIM> t_vec(
        ptrs);
    return t_vec;
  };

  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    const double beta = vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // take into account Jacobean
    const double alpha = t_w * beta;
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows / FIELD_DIM; rr++) {
      // get column base functions gradient at gauss point gg
      auto t_col_base = col_data.getFTensor0N(gg, 0);
      // get mat vec
      auto t_vec = get_t_vec(FIELD_DIM * rr);
      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols / FIELD_DIM; cc++) {
        // calculate element of local matrix
        t_vec(i) += alpha * (t_row_base * t_col_base);
        ++t_col_base;
        ++t_vec;
      }
      ++t_row_base;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;
    ++t_coords;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
};

template <int SPACE_DIM, int S, typename OpBase>
MoFEMErrorCode 
OpGradSymTensorGradImpl<1, SPACE_DIM, SPACE_DIM, S, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;

  const size_t nb_row_dofs = row_data.getIndices().size();
  const size_t nb_col_dofs = col_data.getIndices().size();

  if (nb_row_dofs && nb_col_dofs) {

    FTensor::Index<'i', SPACE_DIM> i;
    FTensor::Index<'j', SPACE_DIM> j;
    FTensor::Index<'k', SPACE_DIM> k;
    FTensor::Index<'l', SPACE_DIM> l;

    // get element volume
    double vol = OpBase::getMeasure();

    // get intergrayion weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();

    // get derivatives of base functions on rows
    auto t_row_diff_base = row_data.getFTensor1DiffN<SPACE_DIM>();

    // Elastic stiffness tensor (4th rank tensor with minor and major
    // symmetry)
    auto t_D = getFTensor4DdgFromMat<SPACE_DIM, SPACE_DIM, S>(*matD);

    // iterate over integration points
    for (int gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

      // calculate scalar weight times element volume
      double a = t_w * vol;

      // iterate over row base functions
      int rr = 0;
      for (; rr != OpBase::nbRows / SPACE_DIM; ++rr) {

        // get sub matrix for the row
        auto t_m = OpBase::template getLocMat<SPACE_DIM>(SPACE_DIM * rr);

        FTensor::Christof<double, SPACE_DIM, SPACE_DIM> t_rowD;
        // I mix up the indices here so that it behaves like a
        // Dg.  That way I don't have to have a separate wrapper
        // class Christof_Expr, which simplifies things.
        t_rowD(l, j, k) = t_D(i, j, k, l) * (a * t_row_diff_base(i));

        // get derivatives of base functions for columns
        auto t_col_diff_base = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);

        // iterate column base functions
        for (int cc = 0; cc != OpBase::nbCols / SPACE_DIM; ++cc) {

          // integrate block local stiffens matrix
          t_m(i, j) += t_rowD(i, j, k) * t_col_diff_base(k);

          // move to next column base function
          ++t_col_diff_base;

          // move to next block of local stiffens matrix
          ++t_m;
        }

        // move to next row base function
        ++t_row_diff_base;
      }

      for (; rr < OpBase::nbRowBaseFunctions; ++rr)
        ++t_row_diff_base;

      // move to next integration weight
      ++t_w;
      ++t_D;
    }
  }

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, int S, typename OpBase>
MoFEMErrorCode 
OpGradTensorGradImpl<1, SPACE_DIM, SPACE_DIM, S, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;

  const size_t nb_row_dofs = row_data.getIndices().size();
  const size_t nb_col_dofs = col_data.getIndices().size();

  if (nb_row_dofs && nb_col_dofs) {

    FTensor::Index<'i', SPACE_DIM> i;
    FTensor::Index<'j', SPACE_DIM> j;
    FTensor::Index<'k', SPACE_DIM> k;
    FTensor::Index<'l', SPACE_DIM> l;

    // get element volume
    double vol = OpBase::getMeasure();

    // get intergrayion weights
    auto t_w = OpBase::getFTensor0IntegrationWeight();

    // get derivatives of base functions on rows
    auto t_row_diff_base = row_data.getFTensor1DiffN<SPACE_DIM>();

    //stiffness tensor (4th rank tensor)
    auto t_D = getFTensor4FromMat<SPACE_DIM, SPACE_DIM, SPACE_DIM, SPACE_DIM, S>(*matD);

    // iterate over integration points
    for (int gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

      // calculate scalar weight times element volume
      double a = t_w * vol;

      // iterate over row base functions
      int rr = 0;
      for (; rr != OpBase::nbRows / SPACE_DIM; ++rr) {

        // get sub matrix for the row
        auto t_m = OpBase::template getLocMat<SPACE_DIM>(SPACE_DIM * rr);

        // get derivatives of base functions for columns
        auto t_col_diff_base = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);

        // iterate column base functions
        for (int cc = 0; cc != OpBase::nbCols / SPACE_DIM; ++cc) {

          // integrate block local stiffens matrix
          t_m(i, k) +=
              a * (t_D(i, j, k, l) * (t_row_diff_base(j) * t_col_diff_base(l)));

          // move to next column base function
          ++t_col_diff_base;

          // move to next block of local stiffens matrix
          ++t_m;
        }

        // move to next row base function
        ++t_row_diff_base;
      }

      for (; rr < OpBase::nbRowBaseFunctions; ++rr)
        ++t_row_diff_base;

      // move to next integration weight
      ++t_w;
      ++t_D;
    }
  }

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode OpMixDivLambdaTimesUImpl<SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;

  const size_t row_nb_dofs = row_data.getIndices().size();
  const size_t col_nb_dofs = col_data.getIndices().size();

  if (row_nb_dofs && col_nb_dofs) {

    auto t_w = this->getFTensor0IntegrationWeight();

    this->locMat.resize(row_nb_dofs, col_nb_dofs, false);
    this->locMat.clear();

    size_t nb_base_functions = row_data.getN().size2() / 3;
    auto t_row_diff_base = row_data.getFTensor2DiffN<3, SPACE_DIM>();

    for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

      const double alpha = alphaConstant * this->getMeasure() * t_w;

      size_t rr = 0;
      for (; rr != row_nb_dofs / SPACE_DIM; ++rr) {
        auto t_mat_diag = getFTensor1FromArrayDiag<SPACE_DIM, SPACE_DIM>(
            this->locMat, SPACE_DIM * rr);
        const double t_row_div_base = t_row_diff_base(i, i);
        auto t_col_base = col_data.getFTensor0N(gg, 0);

        for (size_t cc = 0; cc != col_nb_dofs / SPACE_DIM; ++cc) {
          t_mat_diag(i) += alpha * t_row_div_base * t_col_base;
          ++t_col_base;
          ++t_mat_diag;
        }

        ++t_row_diff_base;
      }
      for (; rr < nb_base_functions; ++rr) 
        ++t_row_diff_base;

      ++t_w;
    }
  }

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode OpMixLambdaTimesGradUImpl<SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;

  const size_t row_nb_dofs = row_data.getIndices().size();
  const size_t col_nb_dofs = col_data.getIndices().size();

  if (row_nb_dofs && col_nb_dofs) {

    auto t_w = this->getFTensor0IntegrationWeight();

    this->locMat.resize(row_nb_dofs, col_nb_dofs, false);
    this->locMat.clear();

    size_t nb_base_functions = row_data.getN().size2() / 3;
    auto t_row_base = row_data.getFTensor1N<3>();

    for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

      const double alpha = alphaConstant * this->getMeasure() * t_w;

      size_t rr = 0;
      for (; rr != row_nb_dofs / SPACE_DIM; ++rr) {
        auto t_mat_diag = getFTensor1FromArrayDiag<SPACE_DIM, SPACE_DIM>(
            this->locMat, SPACE_DIM * rr);
        auto t_col_diff_base = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);

        for (size_t cc = 0; cc != col_nb_dofs / SPACE_DIM; ++cc) {
          t_mat_diag(i) += alpha * t_row_base(j) * t_col_diff_base(j);
          ++t_col_diff_base;
          ++t_mat_diag;
        }

        ++t_row_base;
      }
      for (; rr < nb_base_functions; ++rr)
        ++t_row_base;

      ++t_w;
    }
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM

#endif //__BILINEAR_FORMS_INTEGRATORS_HPP__