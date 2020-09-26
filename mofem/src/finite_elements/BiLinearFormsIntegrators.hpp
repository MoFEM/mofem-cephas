/** \file BiLinearFormsIntegrators.hpp
  * \brief Bilinear forms integrators

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

/**
 * @brief Bilinear integrator form
 * @ingroup mofem_form
 * 
 * @tparam EleOp 
 * @tparam A 
 * @tparam I 
 */
template <typename EleOp>
template <AssemblyType A>
template <IntegrationType I>
struct FormsIntegrators<EleOp>::Assembly<A>::BiLinearForm {

  //! [Grad operator]
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM> struct OpGradGrad;

  /**
   * @brief Integrate \f$(v_{,i},\beta(\mathbf{x}) u_{,j}))_\Omega\f$
   * @ingroup mofem_forms
   * 
   * @tparam SPACE_DIM 
   */
  template <int SPACE_DIM> struct OpGradGrad<1, 1, SPACE_DIM> : public OpBase {

    FTensor::Index<'i', SPACE_DIM> i; ///< summit Index

    OpGradGrad(const std::string row_field_name,
               const std::string col_field_name, ScalarFun beta)
        : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
          betaCoeff(beta) {}

  protected:
    ScalarFun betaCoeff;

    virtual MoFEMErrorCode iNtegrate(EntData &row_data, EntData &col_data) {
      return integrateImpl<GAUSS>(row_data, col_data);
    }

    template <IntegrationType T>
    MoFEMErrorCode integrateImpl(EntData &row_data, EntData &col_data) {
      return MOFEM_NOT_IMPLEMENTED;
    }

    template <>
    MoFEMErrorCode integrateImpl<GAUSS>(EntData &row_data, EntData &col_data) {
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
        const double beta =
            vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
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
  };
  //! [Grad operator]

  //! [Mass operator]
  template <int BASE_DIM, int FIELD_DIM> struct OpMass;

  /**
   * @brief Integrate \f$(v_i,\beta(\mathbf{x}) u_j)_\Omega\f$
   * @ingroup mofem_forms
   * 
   * @tparam  
   */
  template <> struct OpMass<1, 1> : public OpBase {

    OpMass(const std::string row_field_name, const std::string col_field_name,
           ScalarFun beta)
        : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
          betaCoeff(beta) {}

  protected:
    ScalarFun betaCoeff;

    MoFEMErrorCode iNtegrate(EntData &row_data, EntData &col_data) {
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
        const double beta =
            vol * betaCoeff(t_coords(0), t_coords(1), t_coords(2));
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
          ++t_coords;
        ++t_w; // move to another integration weight
      }
      MoFEMFunctionReturn(0);
    };
  };
  //! [Mass operator]

  //! [Grad SymD Grad]

  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S = 0>
  struct OpGradSymTensorGrad;

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
  template <int SPACE_DIM, int S>
  struct OpGradSymTensorGrad<1, SPACE_DIM, SPACE_DIM, S> : public OpBase {

    FTensor::Index<'i', SPACE_DIM> i; ///< summit Index

    OpGradSymTensorGrad(const std::string row_field_name,
                        const std::string col_field_name,
                        boost::shared_ptr<MatrixDouble> mat_D)
        : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL),
          matD(mat_D) {}

  protected:
    boost::shared_ptr<MatrixDouble> matD;

    template <>
    inline FTensor::Ddg<FTensor::PackPtr<double *, S>, 1, 1>
    getFTensor4DdgFromMat(MatrixDouble &data) {
      if (data.size1() != 1)
        THROW_MESSAGE(
            "getFTensor4DdgFromMat<3, 3>: wrong size of data matrix, number "
            "of rows should be 36 but is " +
            boost::lexical_cast<std::string>(data.size1()));

      return FTensor::Ddg<FTensor::PackPtr<double *, S>, 1, 1>{&data(0, 0)};
    }

    /**
     * @brief Get symmetric tensor rank 4  on first two and last indices from
     * form data matrix
     *
     * @param data matrix container which has 36 rows
     * @return FTensor::Ddg<FTensor::PackPtr<T *, 1>, Tensor_Dim01, TensorDim23>
     */
    template <>
    inline FTensor::Ddg<FTensor::PackPtr<double *, S>, 2, 2>
    getFTensor4DdgFromMat(MatrixDouble &data) {
      if(data.size1() != 9) THROW_MESSAGE(
          "getFTensor4DdgFromMat<3, 3>: wrong size of data matrix, number "
          "of rows should be 36 but is " +
          boost::lexical_cast<std::string>(data.size1()));

      return FTensor::Ddg<FTensor::PackPtr<double *, S>, 2, 2>{
          &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0), &data(4, 0),
          &data(5, 0), &data(6, 0), &data(7, 0), &data(8, 0)};
    }

    /**
     * @brief Get symmetric tensor rank 4  on first two and last indices from
     * form data matrix
     *
     * @param data matrix container which has 36 rows
     * @return FTensor::Ddg<FTensor::PackPtr<T *, 1>, Tensor_Dim01, TensorDim23>
     */
    template <>
    inline FTensor::Ddg<FTensor::PackPtr<double *, S>, 3, 3>
    getFTensor4DdgFromMat(MatrixDouble &data) {
      if (data.size1() != 36)
        THROW_MESSAGE(
            "getFTensor4DdgFromMat<3, 3>: wrong size of data matrix, number "
            "of rows should be 36 but is " +
            boost::lexical_cast<std::string>(data.size1()));

      return FTensor::Ddg<FTensor::PackPtr<double *, S>, 3, 3>{
          &data(0, 0),  &data(1, 0),  &data(2, 0),  &data(3, 0),  &data(4, 0),
          &data(5, 0),  &data(6, 0),  &data(7, 0),  &data(8, 0),  &data(9, 0),
          &data(10, 0), &data(11, 0), &data(12, 0), &data(13, 0), &data(14, 0),
          &data(15, 0), &data(16, 0), &data(17, 0), &data(18, 0), &data(19, 0),
          &data(20, 0), &data(21, 0), &data(22, 0), &data(23, 0), &data(24, 0),
          &data(25, 0), &data(26, 0), &data(27, 0), &data(28, 0), &data(29, 0),
          &data(30, 0), &data(31, 0), &data(32, 0), &data(33, 0), &data(34, 0),
          &data(35, 0)};
    }

    virtual MoFEMErrorCode iNtegrate(EntData &row_data, EntData &col_data) {
      return integrateImpl<GAUSS>(row_data, col_data);
    }

    template <IntegrationType T>
    MoFEMErrorCode integrateImpl(EntData &row_data, EntData &col_data) {
      return MOFEM_NOT_IMPLEMENTED;
    }

    template <>
    MoFEMErrorCode integrateImpl<GAUSS>(EntData &row_data, EntData &col_data) {
      MoFEMFunctionBegin;

      const size_t nb_row_dofs = row_data.getIndices().size();
      const size_t nb_col_dofs = col_data.getIndices().size();

      if (nb_row_dofs && nb_col_dofs) {

        // get sub-block (3x3) of local stiffens matrix, here represented by
        // second order tensor
        auto get_tensor2 = [](MatrixDouble &m, const int r, const int c) {
          return FTensor::Tensor2<FTensor::PackPtr<double *, 3>, 3, 3>(
              &m(r + 0, c + 0), &m(r + 0, c + 1), &m(r + 0, c + 2),
              &m(r + 1, c + 0), &m(r + 1, c + 1), &m(r + 1, c + 2),
              &m(r + 2, c + 0), &m(r + 2, c + 1), &m(r + 2, c + 2));
        };

        FTensor::Index<'i', 3> i;
        FTensor::Index<'j', 3> j;
        FTensor::Index<'k', 3> k;
        FTensor::Index<'l', 3> l;

        // get element volume
        double vol = OpBase::getVolume();

        // get intergrayion weights
        auto t_w = OpBase::getFTensor0IntegrationWeight();

        // get derivatives of base functions on rows
        auto t_row_diff_base = row_data.getFTensor1DiffN<3>();

        // Elastic stiffness tensor (4th rank tensor with minor and major
        // symmetry)
        auto t_D = getFTensor4DdgFromMat<SPACE_DIM, SPACE_DIM>(*matD);

        // iterate over integration points
        for (int gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {

          // calculate scalar weight times element volume
          double a = t_w * vol;

          // iterate over row base functions
          int rr = 0;
          for (; rr != OpBase::nbRows / 3; ++rr) {

            // get sub matrix for the row
            auto t_m = OpBase::template getLocMat<SPACE_DIM>(SPACE_DIM * rr);

            FTensor::Christof<double, 3, 3> t_rowD;
            // I mix up the indices here so that it behaves like a
            // Dg.  That way I don't have to have a separate wrapper
            // class Christof_Expr, which simplifies things.
            t_rowD(l, j, k) = t_D(i, j, k, l) * (a * t_row_diff_base(i));

            // get derivatives of base functions for columns
            auto t_col_diff_base = col_data.getFTensor1DiffN<3>(gg, 0);

            // iterate column base functions
            for (int cc = 0; cc != OpBase::nbCols / 3; ++cc) {

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
  };
  //! [Grad Symm D Grad]
};

} // namespace MoFEM

#endif //__BILINEAR_FORMS_INTEGRATORS_HPP__