/** \file BaseDerivativesDataOperators.hpp
  * \brief Base derivatives data operators

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
                        EntitiesFieldData::EntData &data) {
    return doWorkImpl<1>(side, type, data);
  }

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