/** \file HoDataOperators.hpp
  * \brief Operators managing HO geometry

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

#ifndef __HO_DATA_OPERATORS_HPP__
#define __HO_DATA_OPERATORS_HPP__

namespace MoFEM {

/**
 * @brief Calculate HO coordinates at gauss points
 * 
 */
struct OpCalculateHoCoords : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateHoCoords(const std::string field_name)
      : ForcesAndSourcesCore::UserDataOperator(field_name, OPROW) {}


 MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

};

/**
 * @brief Modify integration weights on volume to take in account higher-order
 * geometry
 * @ingroup mofem_forces_and_sources_user_data_operators
 *
 */
struct OpMakeHighOrderGeometryWeightsOnVolume
    : public VolumeElementForcesAndSourcesCoreBase::UserDataOperator {
  OpMakeHighOrderGeometryWeightsOnVolume()
      : VolumeElementForcesAndSourcesCoreBase::UserDataOperator(NOSPACE) {}
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};


};

#endif