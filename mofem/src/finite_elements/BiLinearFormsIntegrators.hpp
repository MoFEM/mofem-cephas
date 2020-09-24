/** \file BiLinearFormsIntegrators.hpp
  * \brief Bilinear forms inteegrators

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

template <typename EleOp>
template <AssemblyType A>
template <IntegrationType I>
struct FormsIntegrators<EleOp>::Assembly<A>::BiLinearForm {

  //! [Grad operator]
  template <int FIELD_DIM, int BASE_DIM, int SPACE_DIM>
  struct OpGradGrad;

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
        // loop over rows base functions
        for (int rr = 0; rr != OpBase::nbRows; rr++) {
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
        ++t_coords;
        ++t_w; // move to another integration weight
      }
      MoFEMFunctionReturn(0);
    }
  };
  //! [Grad operator]

  //! [Mass operator]
  template <int FIELD_DIM, int BASE_DIM> struct OpMass;

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
        for (int rr = 0; rr != OpBase::nbRows; rr++) {
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
        ++t_coords;
        ++t_w; // move to another integration weight
      }
      MoFEMFunctionReturn(0);
    };
  };
  //! [Mass operator]
};

} // namespace MoFEM

#endif //__BILINEAR_FORMS_INTEGRATORS_HPP__