/** \file VolumeElementForcesAndSourcesCoreOnSide.cpp

\brief Implementation of volume element on side

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

int VolumeElementForcesAndSourcesCoreOnSideBase::getRule(int order) {
  return -1;
};

MoFEMErrorCode
VolumeElementForcesAndSourcesCoreOnSideBase::setGaussPts(int order) {
  MoFEMFunctionBegin;

  if (!sidePtrFE)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Side element not set");

  const EntityHandle face_entity =
      sidePtrFE->numeredEntFiniteElementPtr->getEnt();
  SideNumber_multiIndex &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());
  SideNumber_multiIndex::nth_index<0>::type::iterator sit =
      side_table.get<0>().find(face_entity);
  if (sit == side_table.get<0>().end())
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
            "Face can not be found on volume element");

  auto face_ptr_fe =
      static_cast<FaceElementForcesAndSourcesCoreBase *>(sidePtrFE);

  faceSense = (*sit)->sense;
  faceSideNumber = (*sit)->side_number;
  fill(tetConnMap.begin(), tetConnMap.end(), -1);
  for (int nn = 0; nn != 3; ++nn) {
    faceConnMap[nn] =
        std::distance(conn, find(conn, &conn[4], face_ptr_fe->conn[nn]));
    tetConnMap[faceConnMap[nn]] = nn;
    if (faceConnMap[nn] > 3)
      SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
              "No common node on face and element can not be found");
  }

  oppositeNode = std::distance(tetConnMap.begin(),
                               find(tetConnMap.begin(), tetConnMap.end(), -1));

  const int nb_gauss_pts = face_ptr_fe->gaussPts.size2();
  gaussPts.resize(4, nb_gauss_pts, false);
  gaussPts.clear();
  DataForcesAndSourcesCore &dataH1_on_face = *face_ptr_fe->dataOnElement[H1];
  const MatrixDouble &face_shape_funtions =
      dataH1_on_face.dataOnEntities[MBVERTEX][0].getN(NOBASE);
  const double tet_coords[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
  for (int gg = 0; gg != nb_gauss_pts; ++gg) {
    gaussPts(0, gg) =
        face_shape_funtions(gg, 0) * tet_coords[3 * faceConnMap[0] + 0] +
        face_shape_funtions(gg, 1) * tet_coords[3 * faceConnMap[1] + 0] +
        face_shape_funtions(gg, 2) * tet_coords[3 * faceConnMap[2] + 0];
    gaussPts(1, gg) =
        face_shape_funtions(gg, 0) * tet_coords[3 * faceConnMap[0] + 1] +
        face_shape_funtions(gg, 1) * tet_coords[3 * faceConnMap[1] + 1] +
        face_shape_funtions(gg, 2) * tet_coords[3 * faceConnMap[2] + 1];
    gaussPts(2, gg) =
        face_shape_funtions(gg, 0) * tet_coords[3 * faceConnMap[0] + 2] +
        face_shape_funtions(gg, 1) * tet_coords[3 * faceConnMap[1] + 2] +
        face_shape_funtions(gg, 2) * tet_coords[3 * faceConnMap[2] + 2];
    gaussPts(3, gg) = face_ptr_fe->gaussPts(2, gg);
  }
  MoFEMFunctionReturn(0);
}

bool VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::
    checkFECast() const {
  return dynamic_cast<VolumeElementForcesAndSourcesCoreOnSideBase *>(ptrFE) !=
         nullptr;
}

} // namespace MoFEM
