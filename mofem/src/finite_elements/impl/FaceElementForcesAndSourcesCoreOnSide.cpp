/** \file FaceElementForcesAndSourcesCoreOnSide.cpp

\brief Implementation of face element

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

int FaceElementForcesAndSourcesCoreOnSide::getRule(int order) {
  return -1;
};

MoFEMErrorCode
FaceElementForcesAndSourcesCoreOnSide::setGaussPts(int order) {
  MoFEMFunctionBegin;
  if (sidePtrFE == nullptr)
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
            "Pointer to face element is not set");

  const EntityHandle edge_entity =
      sidePtrFE->numeredEntFiniteElementPtr->getEnt();
  SideNumber_multiIndex &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());
  SideNumber_multiIndex::nth_index<0>::type::iterator sit =
      side_table.get<0>().find(edge_entity);
  if (sit == side_table.get<0>().end())
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
            "Edge can not be found on face element");

  auto edge_ptr_fe =
      static_cast<EdgeElementForcesAndSourcesCoreBase *>(sidePtrFE);

  edgeSense = (*sit)->sense;
  edgeSideNumber = (*sit)->side_number;
  fill(faceConnMap.begin(), faceConnMap.end(), -1);
  for (int nn = 0; nn != 2; ++nn) {
    edgeConnMap[nn] = std::distance(
        conn, find(conn, &conn[num_nodes], edge_ptr_fe->cOnn[nn]));
    faceConnMap[edgeConnMap[nn]] = nn;
    if (faceConnMap[nn] >= num_nodes)
      SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
              "No common node on face and element can not be found");
  }

  const int nb_gauss_pts = sidePtrFE->gaussPts.size2();
  gaussPts.resize(3, nb_gauss_pts, false);
  gaussPts.clear();
  EntitiesFieldData &data_h1_on_edge = *edge_ptr_fe->dataOnElement[H1];
  const MatrixDouble &edge_shape_funtions =
      data_h1_on_edge.dataOnEntities[MBVERTEX][0].getN(NOBASE);

  auto set_integration_pts_for_tri = [&]() {
    MoFEMFunctionBegin;
    oppositeNode = std::distance(
        faceConnMap.begin(), find(faceConnMap.begin(), faceConnMap.end(), -1));
    constexpr double face_coords[] = {0, 0, 1, 0, 0, 1};
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      gaussPts(0, gg) =
          edge_shape_funtions(gg, 0) * face_coords[2 * edgeConnMap[0] + 0] +
          edge_shape_funtions(gg, 1) * face_coords[2 * edgeConnMap[1] + 0];
      gaussPts(1, gg) =
          edge_shape_funtions(gg, 0) * face_coords[2 * edgeConnMap[0] + 1] +
          edge_shape_funtions(gg, 1) * face_coords[2 * edgeConnMap[1] + 1];
      gaussPts(2, gg) = edge_ptr_fe->gaussPts(1, gg);
    }
    MoFEMFunctionReturn(0);
  };

  auto set_integration_pts_for_quad = [&]() {
    MoFEMFunctionBegin;
    oppositeNode = -1;
    constexpr double face_coords[] = {0, 0, 1, 0, 1, 1, 0, 1};
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      gaussPts(0, gg) =
          edge_shape_funtions(gg, 0) * face_coords[2 * edgeConnMap[0] + 0] +
          edge_shape_funtions(gg, 1) * face_coords[2 * edgeConnMap[1] + 0];
      gaussPts(1, gg) =
          edge_shape_funtions(gg, 0) * face_coords[2 * edgeConnMap[0] + 1] +
          edge_shape_funtions(gg, 1) * face_coords[2 * edgeConnMap[1] + 1];
      gaussPts(2, gg) = edge_ptr_fe->gaussPts(1, gg);
    }
    MoFEMFunctionReturn(0);
  };

  const auto type = numeredEntFiniteElementPtr->getEntType();

  switch (type) {
  case MBTRI:
    CHKERR set_integration_pts_for_tri();
    break;
  case MBQUAD:
    CHKERR set_integration_pts_for_quad();
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
             "Element type not implemented: %d", type);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
FaceElementForcesAndSourcesCoreOnSide::UserDataOperator::setPtrFE(
    ForcesAndSourcesCore *ptr) {
  MoFEMFunctionBeginHot;
  if (!(ptrFE = dynamic_cast<FaceElementForcesAndSourcesCoreOnSide *>(ptr)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "User operator and finite element do not work together");
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM