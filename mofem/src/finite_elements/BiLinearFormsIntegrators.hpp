/** \file BiLinearFormsIntegrators.hpp
  * \brief Bilinear forms integrators
  * \ingroup mofem_form

  \todo SSome operators could be optimised. To do that, we need to write tests
  and use Valgrind to profile code, shaking cache misses. For example, some
  operators should have iteration first over columns, then rows. ome operators.
  Since those operators are used in many problems, an implementation must be
  efficient.

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
                 const std::string col_field_name, ScalarFun beta,
                 boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        betaCoeff(beta) {
    if (row_field_name == col_field_name)
      this->sYmm = true;
  }

protected:
  ScalarFun betaCoeff;
  boost::shared_ptr<Range> entsPtr;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int FIELD_DIM, int SPACE_DIM, typename OpBase>
struct OpGradGradImpl<1, FIELD_DIM, SPACE_DIM, GAUSS, OpBase> : public OpBase {
  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  OpGradGradImpl(const std::string row_field_name,
                 const std::string col_field_name, ScalarFun beta,
                 boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        betaCoeff(beta) {
    if (row_field_name == col_field_name)
      this->sYmm = true;
  }

protected:
  ScalarFun betaCoeff;
  boost::shared_ptr<Range> entsPtr;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int BASE_DIM, int FIELD_DIM, IntegrationType I, typename OpBase>
struct OpMassImpl {};

template <typename OpBase>
struct OpMassImpl<1, 1, GAUSS, OpBase> : public OpBase {

  OpMassImpl(const std::string row_field_name, const std::string col_field_name,
             ScalarFun beta, boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        betaCoeff(beta), entsPtr(ents_ptr) {
    if (row_field_name == col_field_name)
      this->sYmm = true;
  }

protected:
  ScalarFun betaCoeff;
  boost::shared_ptr<Range> entsPtr;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int FIELD_DIM, typename OpBase>
struct OpMassImpl<1, FIELD_DIM, GAUSS, OpBase>
    : public OpMassImpl<1, 1, GAUSS, OpBase> {
  using OpMassImpl<1, 1, GAUSS, OpBase>::OpMassImpl;

protected:
  using OpMassImpl<1, 1, GAUSS, OpBase>::betaCoeff;
  using OpMassImpl<1, 1, GAUSS, OpBase>::entsPtr;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int BASE_DIM, typename OpBase>
struct OpMassImpl<BASE_DIM, BASE_DIM, GAUSS, OpBase> : public OpBase {

  OpMassImpl(const std::string row_field_name, const std::string col_field_name,
             ScalarFun beta, boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        betaCoeff(beta), entsPtr(ents_ptr) {
    if (row_field_name == col_field_name)
      this->sYmm = true;
  }

protected:
  ScalarFun betaCoeff;
  boost::shared_ptr<Range> entsPtr;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <typename OpBase>
struct OpMassImpl<3, 9, GAUSS, OpBase> : public OpBase {
  OpMassImpl(const std::string row_field_name, const std::string col_field_name,
             ScalarFun beta, boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        betaCoeff(beta), entsPtr(ents_ptr) {
    if (row_field_name == col_field_name)
      this->sYmm = true;
  }

protected:
  ScalarFun betaCoeff;
  boost::shared_ptr<Range> entsPtr;
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
struct OpMixDivTimesScalarImpl {};

template <int SPACE_DIM, typename OpBase>
struct OpMixDivTimesScalarImpl<SPACE_DIM, GAUSS, OpBase> : public OpBase {
  OpMixDivTimesScalarImpl(const std::string row_field_name,
                          const std::string col_field_name,
                          ConstantFun alpha_fun,
                          const bool assemble_transpose = false,
                          const bool only_transpose = false)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        alphaConstant(alpha_fun) {
    this->assembleTranspose = assemble_transpose;
    this->onlyTranspose = only_transpose;
  }

protected:
  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  ConstantFun alphaConstant;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int SPACE_DIM, IntegrationType I, typename OpBase>
struct OpMixDivTimesVecImpl {};

template <int SPACE_DIM, typename OpBase>
struct OpMixDivTimesVecImpl<SPACE_DIM, GAUSS, OpBase> : public OpBase {
  OpMixDivTimesVecImpl(const std::string row_field_name,
                       const std::string col_field_name, ConstantFun alpha_fun,
                       const bool assemble_transpose = false,
                       const bool only_transpose = false)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        alphaConstant(alpha_fun) {
    this->assembleTranspose = assemble_transpose;
    this->onlyTranspose = only_transpose;
  }

protected:
  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  ConstantFun alphaConstant;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int SPACE_DIM, IntegrationType I, typename OpBase,
          CoordinateTypes COORDINATE_SYSTEM>
struct OpMixScalarTimesDivImpl {};

template <int SPACE_DIM, typename OpBase, CoordinateTypes COORDINATE_SYSTEM>
struct OpMixScalarTimesDivImpl<SPACE_DIM, GAUSS, OpBase, COORDINATE_SYSTEM>
    : public OpBase {
  OpMixScalarTimesDivImpl(const std::string row_field_name,
                          const std::string col_field_name,
                          ScalarFun alpha_fun,
                          const bool assemble_transpose = false,
                          const bool only_transpose = false)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        alphaConstant(alpha_fun) {
    this->assembleTranspose = assemble_transpose;
    this->onlyTranspose = only_transpose;
    this->sYmm = false;
  }

protected:
  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  ScalarFun alphaConstant;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, IntegrationType I,
          typename OpBase>
struct OpMixVectorTimesGradImpl;

template <int SPACE_DIM, typename OpBase>
struct OpMixVectorTimesGradImpl<3, SPACE_DIM, SPACE_DIM, GAUSS, OpBase>
    : public OpBase {
  OpMixVectorTimesGradImpl(const std::string row_field_name,
                           const std::string col_field_name,
                           ConstantFun alpha_fun,
                           const bool assemble_transpose = false,
                           const bool only_transpose = false)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        alphaConstant(alpha_fun) {
    this->assembleTranspose = assemble_transpose;
    this->onlyTranspose = only_transpose;
  }

protected:
  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  ConstantFun alphaConstant;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int SPACE_DIM, typename OpBase>
struct OpMixVectorTimesGradImpl<1, SPACE_DIM, SPACE_DIM, GAUSS, OpBase>
    : public OpBase {
  OpMixVectorTimesGradImpl(const std::string row_field_name,
                           const std::string col_field_name,
                           ConstantFun alpha_fun,
                           const bool assemble_transpose = false,
                           const bool only_transpose = false)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        alphaConstant(alpha_fun) {
    this->assembleTranspose = assemble_transpose;
    this->onlyTranspose = only_transpose;
  }

protected:
  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  ConstantFun alphaConstant;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int SPACE_DIM, IntegrationType I, typename OpBase>
struct OpMixTensorTimesGradImpl {};

template <int SPACE_DIM, typename OpBase>
struct OpMixTensorTimesGradImpl<SPACE_DIM, GAUSS, OpBase> : public OpBase {
  OpMixTensorTimesGradImpl(const std::string row_field_name,
                           const std::string col_field_name,
                           ConstantFun alpha_fun,
                           const bool assemble_transpose = false,
                           const bool only_transpose = false)
      : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
        alphaConstant(alpha_fun) {
    this->assembleTranspose = assemble_transpose;
    this->onlyTranspose = only_transpose;
  }

protected:
  FTensor::Index<'i', SPACE_DIM> i; ///< summit Index
  FTensor::Index<'j', SPACE_DIM> j; ///< summit Index
  ConstantFun alphaConstant;
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, IntegrationType I,
          typename OpBase>
struct OpConvectiveTermLhsDuImpl;

template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, IntegrationType I,
          typename OpBase>
struct OpConvectiveTermLhsDyImpl;

template <int SPACE_DIM, typename OpBase>
struct OpConvectiveTermLhsDuImpl<1, 1, SPACE_DIM, GAUSS, OpBase>
    : public OpBase {
  OpConvectiveTermLhsDuImpl(const std::string field_name_row,
                            const std::string field_name_col,
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
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int SPACE_DIM, typename OpBase>
struct OpConvectiveTermLhsDyImpl<1, 1, SPACE_DIM, GAUSS, OpBase>
    : public OpBase {
  OpConvectiveTermLhsDyImpl(const std::string field_name_row,
                            const std::string field_name_col,
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
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int FIELD_DIM, int SPACE_DIM, typename OpBase>
struct OpConvectiveTermLhsDuImpl<1, FIELD_DIM, SPACE_DIM, GAUSS, OpBase>
    : public OpBase {
  OpConvectiveTermLhsDuImpl(const std::string field_name_row,
                            const std::string field_name_col,
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
  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
                           DataForcesAndSourcesCore::EntData &col_data);
};

template <int FIELD_DIM, int SPACE_DIM, typename OpBase>
struct OpConvectiveTermLhsDyImpl<1, FIELD_DIM, SPACE_DIM, GAUSS, OpBase>
    : public OpBase {
  OpConvectiveTermLhsDyImpl(const std::string field_name_row,
                            const std::string field_name_col,
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

  /**
   * @brief Integrate \f$(\lambda_{i,i},u)_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM>
  struct OpMixDivTimesScalar
      : public OpMixDivTimesScalarImpl<SPACE_DIM, I, OpBase> {
    using OpMixDivTimesScalarImpl<SPACE_DIM, I,
                                  OpBase>::OpMixDivTimesScalarImpl;
  };

  /**
   * @brief Integrate \f$(\lambda_{ij,j},u_{i})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM>
  struct OpMixDivTimesVec : public OpMixDivTimesVecImpl<SPACE_DIM, I, OpBase> {
    using OpMixDivTimesVecImpl<SPACE_DIM, I, OpBase>::OpMixDivTimesVecImpl;
  };

  /**
   * @brief Integrate \f$(\lambda,u_{i,i})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM, CoordinateTypes COORDINATE_SYSTEM = CARTESIAN>
  struct OpMixScalarTimesDiv
      : public OpMixScalarTimesDivImpl<SPACE_DIM, I, OpBase,
                                       COORDINATE_SYSTEM> {
    using OpMixScalarTimesDivImpl<SPACE_DIM, I, OpBase,
                                  COORDINATE_SYSTEM>::OpMixScalarTimesDivImpl;
  };

  /**
   * @brief Integrate \f$(\lambda_{i},u_{,j})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM>
  struct OpMixVectorTimesGrad
      : public OpMixVectorTimesGradImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I,
                                        OpBase> {
    using OpMixVectorTimesGradImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I,
                                   OpBase>::OpMixVectorTimesGradImpl;
  };

  /**
   * @brief Integrate \f$(\lambda_{ij},u_{i,j})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM>
  struct OpMixTensorTimesGrad
      : public OpMixTensorTimesGradImpl<SPACE_DIM, I, OpBase> {
    using OpMixTensorTimesGradImpl<SPACE_DIM, I,
                                   OpBase>::OpMixTensorTimesGradImpl;
  };

  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM>
  struct OpConvectiveTermLhsDu
      : public OpConvectiveTermLhsDuImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I,
                                         OpBase> {
    using OpConvectiveTermLhsDuImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I,
                                    OpBase>::OpConvectiveTermLhsDuImpl;
  };

  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM>
  struct OpConvectiveTermLhsDy
      : public OpConvectiveTermLhsDyImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I,
                                         OpBase> {
    using OpConvectiveTermLhsDyImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I,
                                    OpBase>::OpConvectiveTermLhsDyImpl;
  };
};

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode OpGradGradImpl<1, 1, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
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
  auto t_row_grad = row_data.getFTensor1DiffN<SPACE_DIM>();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();

  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    const double beta = vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // take into account Jacobian
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

template <int FIELD_DIM, int SPACE_DIM, typename OpBase>
MoFEMErrorCode
OpGradGradImpl<1, FIELD_DIM, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
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
  auto t_row_grad = row_data.getFTensor1DiffN<SPACE_DIM>();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();

  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'j', FIELD_DIM> j;

  auto get_t_vec = [&](const int rr) {
    std::array<double *, FIELD_DIM> ptrs;
    for (auto i = 0; i != FIELD_DIM; ++i)
      ptrs[i] = &OpBase::locMat(rr + i, i);
    return FTensor::Tensor1<FTensor::PackPtr<double *, FIELD_DIM>, FIELD_DIM>(
        ptrs);
  };

  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    const double beta = vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // take into account Jacobian
    const double alpha = t_w * beta;
    // loop over ros base functions
    int rr = 0;
    for (; rr != OpBase::nbRows / FIELD_DIM; rr++) {
      // get diag vec
      auto t_vec = get_t_vec(rr * FIELD_DIM);
      // get column base functions gradient at gauss point gg
      auto t_col_grad = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);
      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols / FIELD_DIM; cc++) {
        // calculate element of local matrix
        t_vec(j) += alpha * (t_row_grad(i) * t_col_grad(i));
        ++t_col_grad; // move to another gradient of base function
                      // on column
        ++t_vec;
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
    const double beta = vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // take into account Jacobian
    const double alpha = t_w * beta;
    // loop over rows base functions
    auto a_mat_ptr = &*OpBase::locMat.data().begin();
    int rr = 0;
    for (; rr != OpBase::nbRows; rr++) {
      // get column base functions gradient at gauss point gg
      auto t_col_base = col_data.getFTensor0N(gg, 0);
      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols; cc++) {
        // calculate element of local matrix
        *a_mat_ptr += alpha * (t_row_base * t_col_base);
        ++t_col_base;
        ++a_mat_ptr;
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

  FTensor::Index<'i', FIELD_DIM> i;
  auto get_t_vec = [&](const int rr) {
    std::array<double *, FIELD_DIM> ptrs;
    for (auto i = 0; i != FIELD_DIM; ++i)
      ptrs[i] = &OpBase::locMat(rr + i, i);
    return FTensor::Tensor1<FTensor::PackPtr<double *, FIELD_DIM>, FIELD_DIM>(
        ptrs);
  };

  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    const double beta = vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // take into account Jacobian
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

template <int BASE_DIM, typename OpBase>
MoFEMErrorCode OpMassImpl<BASE_DIM, BASE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  FTensor::Index<'i', BASE_DIM> i;
  MoFEMFunctionBegin;
  if (entsPtr) {
    if (entsPtr->find(OpBase::getFEEntityHandle()) == entsPtr->end())
      MoFEMFunctionReturnHot(0);
  }
  size_t nb_base_functions = row_data.getN().size2() / BASE_DIM;
  // // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor1N<BASE_DIM>();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    const double beta = vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // take into account Jacobian
    const double alpha = t_w * beta;
    // loop over rows base functions
    auto a_mat_ptr = &*OpBase::locMat.data().begin();
    int rr = 0;
    for (; rr != OpBase::nbRows; rr++) {
      // get column base functions gradient at gauss point gg
      auto t_col_base = col_data.getFTensor1N<BASE_DIM>(gg, 0);
      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols; cc++) {
        // calculate element of local matrix
        (*a_mat_ptr) += alpha * (t_row_base(i) * t_col_base(i));
        ++t_col_base;
        ++a_mat_ptr;
      }
      ++t_row_base;
    }
    for (; rr < nb_base_functions; ++rr)
      ++t_row_base;
    ++t_coords;
    ++t_w; // move to another integration weight
  }
  MoFEMFunctionReturn(0);
};

template <typename OpBase>
MoFEMErrorCode OpMassImpl<3, 9, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'k', 3> k;
  if (entsPtr) {
    if (entsPtr->find(OpBase::getFEEntityHandle()) == entsPtr->end())
      MoFEMFunctionReturnHot(0);
  }
  auto get_t_vec = [&](const int rr) {
    std::array<double *, 3> ptrs;
    for (auto i = 0; i != 3; ++i)
      ptrs[i] = &OpBase::locMat(rr + i, i);
    return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptrs);
  };
  size_t nb_base_functions = row_data.getN().size2() / 3;
  // // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor1N<3>();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {
    const double beta = vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
    // take into account Jacobian
    const double alpha = t_w * beta;
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows / 3; rr++) {
      // get column base functions gradient at gauss point gg
      auto t_col_base = col_data.getFTensor1N<3>(gg, 0);
      auto t_vec = get_t_vec(3 * rr);
      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols / 3; cc++) {
        // calculate element of local matrix
        t_vec(i) += alpha * (t_row_base(k) * t_col_base(k));
        ++t_col_base;
        ++t_vec;
      }
      ++t_row_base;
    }
    for (; rr < nb_base_functions; ++rr)
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

    // stiffness tensor (4th rank tensor)
    auto t_D =
        getFTensor4FromMat<SPACE_DIM, SPACE_DIM, SPACE_DIM, SPACE_DIM, S>(
            *matD);

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
MoFEMErrorCode OpMixDivTimesScalarImpl<SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;

  auto t_w = this->getFTensor0IntegrationWeight();

  size_t nb_base_functions = row_data.getN().size2() / 3;
  auto t_row_diff_base = row_data.getFTensor2DiffN<3, SPACE_DIM>();

  const double alpha_constant = alphaConstant();

  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {
    const double alpha = alpha_constant * this->getMeasure() * t_w;

    size_t rr = 0;
    for (; rr != OpBase::nbRows; ++rr) {
      const double t_row_div_base = t_row_diff_base(i, i);
      auto t_col_base = col_data.getFTensor0N(gg, 0);
      for (size_t cc = 0; cc != OpBase::nbCols; ++cc) {
        this->locMat(rr, cc) += alpha * t_row_div_base * t_col_base;
        ++t_col_base;
      }
      ++t_row_diff_base;
    }
    for (; rr < nb_base_functions; ++rr)
      ++t_row_diff_base;

    ++t_w;
  }

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode OpMixDivTimesVecImpl<SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;

  auto t_w = this->getFTensor0IntegrationWeight();

  size_t nb_base_functions = row_data.getN().size2() / 3;
  auto t_row_diff_base = row_data.getFTensor2DiffN<3, SPACE_DIM>();
  const double alpha_constant = alphaConstant();
  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    const double alpha = alpha_constant * this->getMeasure() * t_w;

    size_t rr = 0;
    for (; rr != OpBase::nbRows / SPACE_DIM; ++rr) {
      auto t_mat_diag = getFTensor1FromArrayDiag<SPACE_DIM, SPACE_DIM>(
          this->locMat, SPACE_DIM * rr);
      const double t_row_div_base = t_row_diff_base(i, i);
      auto t_col_base = col_data.getFTensor0N(gg, 0);

      for (size_t cc = 0; cc != OpBase::nbCols / SPACE_DIM; ++cc) {
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

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase, CoordinateTypes COORDINATE_SYSTEM>
MoFEMErrorCode
OpMixScalarTimesDivImpl<SPACE_DIM, GAUSS, OpBase, COORDINATE_SYSTEM>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;

#ifndef NDEBUG
  if (OpBase::locMat.size2() % SPACE_DIM)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Number of rows in matrix should be multiple of space dimensions");
#endif

  // When we move to C++17 add if constexpr()
  if constexpr (COORDINATE_SYSTEM == POLAR || COORDINATE_SYSTEM == SPHERICAL)
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
             "%s coordiante not implemented",
             CoordinateTypesNames[COORDINATE_SYSTEM]);

  auto t_w = this->getFTensor0IntegrationWeight();
  auto t_coords = this->getFTensor1CoordsAtGaussPts();
  size_t nb_base_functions_row = row_data.getN().size2();
  auto t_row_base = row_data.getFTensor0N();
  const double vol = this->getMeasure();
  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    const double alpha =
        alphaConstant(t_coords(0), t_coords(1), t_coords(2)) * t_w * vol;        

    size_t rr = 0;
    auto t_m = getFTensor1FromPtr<SPACE_DIM>(OpBase::locMat.data().data());

    // When we move to C++17 add if constexpr()
    if constexpr (COORDINATE_SYSTEM == CARTESIAN) {
      for (; rr != OpBase::nbRows; ++rr) {
        const double r_val = alpha * t_row_base;
        auto t_col_diff_base = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);
        for (size_t cc = 0; cc != OpBase::nbCols / SPACE_DIM; ++cc) {
          t_m(i) += r_val * t_col_diff_base(i);
          ++t_col_diff_base;
          ++t_m;
        }
        ++t_row_base;
      }
    }

    // When we move to C++17 add if constexpr()
    if constexpr (COORDINATE_SYSTEM == CYLINDRICAL) {
      for (; rr != OpBase::nbRows; ++rr) {
        const double r_val = alpha * t_row_base;
        auto t_col_base = col_data.getFTensor0N(gg, 0);
        auto t_col_diff_base = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);
        for (size_t cc = 0; cc != OpBase::nbCols / SPACE_DIM; ++cc) {
          t_m(i) += r_val * t_col_diff_base(i);
          t_m(0) += (r_val / t_coords(0)) * t_col_base;
          ++t_col_base;
          ++t_col_diff_base;
          ++t_m;
        }
        ++t_row_base;
      }
    }

    for (; rr < nb_base_functions_row; ++rr)
      ++t_row_base;

    ++t_w;
    ++t_coords;
  }

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode
OpMixVectorTimesGradImpl<3, SPACE_DIM, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;

  auto t_w = this->getFTensor0IntegrationWeight();

  size_t nb_base_functions = row_data.getN().size2() / 3;
  auto t_row_base = row_data.getFTensor1N<3>();
  auto &mat = this->locMat;
  const double alpha_constant = alphaConstant();
  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    const double alpha = alpha_constant * this->getMeasure() * t_w;

    size_t rr = 0;
    for (; rr != OpBase::nbRows; ++rr) {
      auto t_col_diff_base = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);
      for (size_t cc = 0; cc != OpBase::nbCols; ++cc) {
        mat(rr, cc) += alpha * t_row_base(i) * t_col_diff_base(i);
        ++t_col_diff_base;
      }
      ++t_row_base;
    }
    for (; rr < nb_base_functions; ++rr)
      ++t_row_base;

    ++t_w;
  }

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode
OpMixVectorTimesGradImpl<1, SPACE_DIM, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;

  auto t_w = this->getFTensor0IntegrationWeight();

  size_t nb_base_functions = row_data.getN().size2();
  auto t_row_base = row_data.getFTensor0N();

  auto get_t_vec = [&](const int rr) {
    std::array<double *, SPACE_DIM> ptrs;
    for (auto i = 0; i != SPACE_DIM; ++i)
      ptrs[i] = &OpBase::locMat(rr + i, 0);
    return FTensor::Tensor1<FTensor::PackPtr<double *, 1>, SPACE_DIM>(ptrs);
  };

  const double alpha_constant = alphaConstant();
  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    const double alpha = alpha_constant * this->getMeasure() * t_w;

    size_t rr = 0;
    for (; rr != OpBase::nbRows / SPACE_DIM; ++rr) {
      auto t_vec = get_t_vec(SPACE_DIM * rr);
      auto t_col_diff_base = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);
      for (size_t cc = 0; cc != OpBase::nbCols; ++cc) {
        t_vec(i) += alpha * t_row_base * t_col_diff_base(i);
        ++t_col_diff_base;
        ++t_vec;
      }
      ++t_row_base;
    }
    for (; rr < nb_base_functions; ++rr)
      ++t_row_base;

    ++t_w;
  }

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode OpMixTensorTimesGradImpl<SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;

  auto t_w = this->getFTensor0IntegrationWeight();

  size_t nb_base_functions = row_data.getN().size2() / 3;
  auto t_row_base = row_data.getFTensor1N<3>();
  const double alpha_constant = alphaConstant();
  for (size_t gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

    const double alpha = alpha_constant * this->getMeasure() * t_w;

    size_t rr = 0;
    for (; rr != OpBase::nbRows / SPACE_DIM; ++rr) {
      auto t_mat_diag = getFTensor1FromArrayDiag<SPACE_DIM, SPACE_DIM>(
          this->locMat, SPACE_DIM * rr);
      auto t_col_diff_base = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);

      for (size_t cc = 0; cc != OpBase::nbCols / SPACE_DIM; ++cc) {
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

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, typename OpBase>
MoFEMErrorCode
OpConvectiveTermLhsDuImpl<1, 1, SPACE_DIM, GAUSS, OpBase>::iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
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
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
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
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();

  auto t_grad_y = getFTensor2FromMat<FIELD_DIM, SPACE_DIM>(*yGradPtr);
  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'J', FIELD_DIM> j;

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
        t_mat(j, i) += (alpha * t_row_base * t_col_base) * t_grad_y(j, i);
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
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data) {
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
  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'j', FIELD_DIM> j;
  FTensor::Index<'k', FIELD_DIM> k;

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
        t_mat(j, k) += alpha * t_row_base * (t_diff_col_base(i) * t_u(i));
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

#endif //__BILINEAR_FORMS_INTEGRATORS_HPP__