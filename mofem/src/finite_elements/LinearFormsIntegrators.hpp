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

  template <int FIELD_DIM, int BASE_DIM> struct OpSource;

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
};

} // namespace MoFEM

#endif // __LINEAR_FORMS_INTEGRATORS_HPP__