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

template <int BASE_DIM>
struct OpBaseDerivativesNext;

template <> struct OpBaseDerivativesNext<1> : public OpBaseDerivativesBase {

  OpBaseDerivativesNext(int derivative,
                        boost::shared_ptr<MatrixDouble> base_mass_ptr,
                        boost::shared_ptr<EntitiesFieldData> data_l2,
                        const FieldApproximationBase b, const FieldSpace s,
                        int verb = QUIET, Sev sev = Sev::verbose);

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  int calcBaseDerivative;
  MatrixDouble nF;

  template <int SPACE_DIM>
  MoFEMErrorCode setBase(EntitiesFieldData::EntData &data,
                         EntitiesFieldData::EntData &ent_data);
};

} // namespace MoFEM

#endif //__BASE_DIREVATIVES_DATA_OPERATORS_HPP__