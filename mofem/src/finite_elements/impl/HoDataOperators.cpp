/** \file UserDataOperators.cpp

\brief Generic user data operators for evaluate fields, and other common
purposes.

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

namespace MoFEM {

MoFEMErrorCode
OpCalculateHoCoords::doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
  FTensor::Index<'i', 3> i;
  MoFEMFunctionBegin;
  const auto nb_dofs = data.getFieldData().size() / 3;
  if (nb_dofs) {
    if (type == MBVERTEX)
      getCoordsAtGaussPts().clear();
    auto t_base = data.getFTensor0N();
    auto t_coords = getFTensor1CoordsAtGaussPts();
    const auto nb_integration_pts = data.getN().size1();
    const auto nb_base_functions = data.getN().size2();
    for (auto gg = 0; gg != nb_integration_pts; ++gg) {
      auto t_dof = data.getFTensor1FieldData<3>();
      size_t bb = 0;
      for (; bb != nb_dofs; ++bb) {
        t_coords(i) += t_base * t_dof(i);
        ++t_base;
      }
      for (; bb != nb_base_functions; ++bb)
        ++t_base;
      ++t_coords;
    }
  }
  MoFEMFunctionReturn(0);
};

MoFEMErrorCode OpMakeHighOrderGeometryWeightsOnVolume::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  const size_t nb_int_pts = getGaussPts().size2();

  if (getHoGaussPtsDetJac().size() == nb_int_pts) {
    const double a = getMeasure();
    auto t_w = getFTensor0IntegrationWeight();
    auto t_w_ho = getFTensor0FromVec(getHoGaussPtsDetJac());

    for (size_t gg = 0; gg != nb_int_pts; ++gg) {
      t_w *= t_w_ho;
      ++t_w;
      ++t_w_ho;
    }
  } else {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
             "Number of rows in getHoGaussPtsDetJac should be equal to "
             "number of integration points, but is not, i.e. %d != %d",
             getHoGaussPtsDetJac().size(), nb_int_pts);
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM