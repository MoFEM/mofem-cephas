/** \file BaseDerivativesDataOperators.hpp
  * \brief Base derivatives data operators

*/

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __BASE_DIREVATIVES_DATA_OPERATORS_HPP__
#define __BASE_DIREVATIVES_DATA_OPERATORS_HPP__

namespace MoFEM {

struct OpBaseDerivativesBase : public ForcesAndSourcesCore::UserDataOperator {

  OpBaseDerivativesBase(boost::shared_ptr<MatrixDouble> base_mass_ptr,
                        boost::shared_ptr<EntitiesFieldData> data_l2,
                        const FieldApproximationBase b, const FieldSpace s,
                        int verb = QUIET, Sev sev = Sev::verbose);

protected:
  FieldApproximationBase base;
  int verbosity;
  Sev severityLevel;

  boost::shared_ptr<MatrixDouble> baseMassPtr;
  boost::shared_ptr<EntitiesFieldData> dataL2;
};

template <int BASE_DIM>
struct OpBaseDerivativesMass;

template <> struct OpBaseDerivativesMass<1> : public OpBaseDerivativesBase  {

  using OpBaseDerivativesBase::OpBaseDerivativesBase;

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

template <> struct OpBaseDerivativesMass<3> : public OpBaseDerivativesMass<1> {
  using OpBaseDerivativesMass<1>::OpBaseDerivativesMass;
};

template <int DIM> struct OpBaseDerivativesSetHOInvJacobian;

template <>
struct OpBaseDerivativesSetHOInvJacobian<2>
    : public OpSetInvJacSpaceForFaceImpl<2, 1> {

  OpBaseDerivativesSetHOInvJacobian(
      boost::shared_ptr<EntitiesFieldData> data_l2,
      boost::shared_ptr<MatrixDouble> inv_jac_ptr);

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  boost::shared_ptr<EntitiesFieldData> dataL2;
};

template <int BASE_DIM>
struct OpBaseDerivativesNext;

/**
 * @brief Specialisation for calculate directives for scalar base functions
 * 
 * @tparam  
 */
template <> struct OpBaseDerivativesNext<1> : public OpBaseDerivativesBase {

  OpBaseDerivativesNext(int derivative,
                        boost::shared_ptr<MatrixDouble> base_mass_ptr,
                        boost::shared_ptr<EntitiesFieldData> data_l2,
                        const FieldApproximationBase b, const FieldSpace s,
                        int verb = QUIET, Sev sev = Sev::verbose);

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  int calcBaseDerivative;
  MatrixDouble nF;

  template <int BASE_DIM>
  MoFEMErrorCode doWorkImpl(int side, EntityType type,
                            EntitiesFieldData::EntData &data);

  template <int BASE_DIM, int SPACE_DIM>
  MoFEMErrorCode setBaseImpl(EntitiesFieldData::EntData &data,
                             EntitiesFieldData::EntData &ent_data);
};

/**
 * @brief Specialisation for calculate directives for scalar base functions
 * 
 * @tparam  
 */
template <> struct OpBaseDerivativesNext<3> : public OpBaseDerivativesNext<1> {

  using OpBaseDerivativesNext<1>::OpBaseDerivativesNext;
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

} // namespace MoFEM

#endif //__BASE_DIREVATIVES_DATA_OPERATORS_HPP__