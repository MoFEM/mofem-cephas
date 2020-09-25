/** \file LinearFormsIntegrators.hpp
  * \brief Linear forms inteegrators

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

template <typename EleOp>
template <AssemblyType A>
template <IntegrationType I>
struct FormsIntegrators<EleOp>::Assembly<A>::LinearForm {

  //! [Source operator]
  template <int BASE_DIM, int FIELD_DIM> struct OpSource;

  template <> struct OpSource<1, 1> : public OpBase {

    OpSource(const std::string field_name, ScalarFun source_fun)
        : OpBase(field_name, field_name, OpBase::OPROW), sourceFun(source_fun) {
    }

  protected:
    ScalarFun sourceFun;

    virtual MoFEMErrorCode iNtegrate(EntData &data) {
      return integrateImpl<GAUSS>(data);
    }

    template <IntegrationType T>
    MoFEMErrorCode integrateImpl(EntData &row_data) {
      return MOFEM_NOT_IMPLEMENTED;
    }

    template <> inline MoFEMErrorCode integrateImpl<GAUSS>(EntData &row_data) {
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
        for (int rr = 0; rr != OpBase::nbRows; ++rr) {
          OpBase::locF[rr] += alpha * t_row_base;
          ++t_row_base;
        }
        ++t_coords;
        ++t_w; // move to another integration weight
      }
      MoFEMFunctionReturn(0);
    }
  };
  //! [Source operator]

  //! [Grad times tensor]
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM>
  struct OpGradTimesTensor;

  template <int SPACE_DIM>
  struct OpGradTimesTensor<1, 1, SPACE_DIM> : public OpBase {

    FTensor::Index<'i', SPACE_DIM> i; ///< summit Index

    OpGradTimesTensor(const std::string field_name, ScalarFun beta,
                      boost::shared_ptr<MatrixDouble> &vector_vals)
        : OpBase(field_name, field_name, OpBase::OPROW),
          vectorVals(vector_vals), betaCoeff(beta) {}

  protected:
    ScalarFun betaCoeff;
    boost::shared_ptr<MatrixDouble> vectorVals;

    virtual MoFEMErrorCode iNtegrate(EntData &data) {
      return integrateImpl<GAUSS>(data);
    }

    template <IntegrationType T>
    MoFEMErrorCode integrateImpl(EntData &row_data) {
      return MOFEM_NOT_IMPLEMENTED;
    }

    template <> inline MoFEMErrorCode integrateImpl<GAUSS>(EntData &row_data) {
      MoFEMFunctionBegin;
      // get element volume
      const double vol = OpBase::getMeasure();
      // get integration weights
      auto t_w = OpBase::getFTensor0IntegrationWeight();
      // get base function gradient on rows
      auto t_row_grad = row_data.getFTensor1DiffN<SPACE_DIM>();
      // get filed gradient values
      auto t_val_grad = getFTensor1FromMat<SPACE_DIM>(*(vectorVals));
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
          // calculate element of local matrix
          OpBase::locF[rr] += alpha * (t_row_grad(i) * t_val_grad(i));
          ++t_row_grad; // move to another element of gradient of base
                        // function on row
        }
        ++t_coords;
        ++t_val_grad;
        ++t_w; // move to another integration weight
      }
      MoFEMFunctionReturn(0);
    }
  };
  //! [Grad times tensor]

  //! [Grad times symmetric tensor]
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM>
  struct OpGradTimesSymmetricTensor;

  template <int SPACE_DIM>
  struct OpGradTimesSymmetricTensor<1, SPACE_DIM, SPACE_DIM> : public OpBase {

    FTensor::Index<'i', SPACE_DIM> i; ///< summit Index

    OpGradTimesSymmetricTensor(const std::string field_name, ScalarFun beta,
                               boost::shared_ptr<MatrixDouble> &mat_vals)
        : OpBase(field_name, field_name, OpBase::OPROW), matVals(mat_vals),
          betaCoeff(beta) {}

  protected:
    ScalarFun betaCoeff;
    boost::shared_ptr<MatrixDouble> matVals;

    virtual MoFEMErrorCode iNtegrate(EntData &data) {
      return integrateImpl<GAUSS>(data);
    }

    template <IntegrationType T>
    MoFEMErrorCode integrateImpl(EntData &row_data) {
      return MOFEM_NOT_IMPLEMENTED;
    }

    template <> inline MoFEMErrorCode integrateImpl<GAUSS>(EntData &row_data) {
      MoFEMFunctionBegin;
      // get element volume
      const double vol = OpBase::getMeasure();
      // get integration weights
      auto t_w = OpBase::getFTensor0IntegrationWeight();
      // get base function gradient on rows
      auto t_row_grad = row_data.getFTensor1DiffN<SPACE_DIM>();
      // get filed gradient values
      auto t_val_mat = getFTensor2SymmetricFromMat<SPACE_DIM>(*(matVals));
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
          // calculate element of local matrix
          OpBase::locF[rr] += alpha * (t_row_grad(i) * t_val_mat(i));
          ++t_row_grad; // move to another element of gradient of base
                        // function on row
        }
        ++t_coords;
        ++t_val_mat;
        ++t_w; // move to another integration weight
      }
      MoFEMFunctionReturn(0);
    }
  };
  //! [Grad times tensor]
};

} // namespace MoFEM

#endif // __LINEAR_FORMS_INTEGRATORS_HPP__