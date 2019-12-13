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

int FaceElementForcesAndSourcesCoreOnVolumeSideBase::getRule(int order) {
  return -1;
};

MoFEMErrorCode
FaceElementForcesAndSourcesCoreOnVolumeSideBase::setGaussPts(int order) {
  MoFEMFunctionBegin;
  if (sidePtrFE == nullptr)
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
            "Pointer to face element is not set");

  const EntityHandle volume_entity =
      sidePtrFE->numeredEntFiniteElementPtr->getEnt();
  Range adj_faces;
  CHKERR sidePtrFE->mField.getInterface<BitRefManager>()->getAdjacenciesAny(
      volume_entity, 2, adj_faces);
  SideNumber_multiIndex &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());

  SideNumber_multiIndex &side_volume_table =
      const_cast<SideNumber_multiIndex &>(
          sidePtrFE->numeredEntFiniteElementPtr->getSideNumberTable());

  adj_faces = adj_faces.subset_by_type(MBTRI);
  bool face_common = false;
  SideNumber_multiIndex::nth_index<0>::type::iterator sit;
  SideNumber_multiIndex::nth_index<0>::type::iterator sit_volume;
  EntityHandle face_ent;
  for (Range::iterator it_tris = adj_faces.begin(); it_tris != adj_faces.end();
       it_tris++) {
    sit = side_table.get<0>().find(*it_tris);
    sit_volume = side_volume_table.get<0>().find(*it_tris);
    if (sit != side_table.get<0>().end()) {
      face_ent = *it_tris;
      face_common = true;
      break;
    }
  }
  if (!face_common) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Error no common tets");
  }

  auto contact_prism_ptr_fe =
      static_cast<ContactPrismElementForcesAndSourcesCore *>(sidePtrFE);

  faceSense = (*sit)->sense;
  faceSideNumber = (*sit)->side_number;
  int side_of_vol_number = (*sit_volume)->side_number;
  const EntityHandle *ent_on_vol;
  int num_nodes = 3;
  CHKERR contact_prism_ptr_fe->mField.get_moab().get_connectivity(
      face_ent, ent_on_vol, num_nodes, true);

  // const EntityHandle edge_entity =
  //     sidePtrFE->numeredEntFiniteElementPtr->getEnt();
  // SideNumber_multiIndex &side_table = const_cast<SideNumber_multiIndex &>(
  //     numeredEntFiniteElementPtr->getSideNumberTable());
  // SideNumber_multiIndex::nth_index<0>::type::iterator sit =
  //     side_table.get<0>().find(edge_entity);
  // if (sit == side_table.get<0>().end())
  //   SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
  //           "Edge can not be found on face element");

  // auto face_ptr_fe =
  //     static_cast<FaceElementForcesAndSourcesCoreBase *>(sidePtrFE);

  fill(faceConnMap.begin(), faceConnMap.end(), -1);
  for (int nn = 0; nn != 3; ++nn) {
    faceConnMap[nn] = std::distance(conn, find(conn, &conn[3], ent_on_vol[nn]));

    faceConnMap[faceConnMap[nn]] = nn;
    if (faceConnMap[nn] > 2)
      SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
              "No common node on face and element can not be found");
  }

  const int nb_gauss_pts =
      contact_prism_ptr_fe->getGaussPtsMasterFromEleSide().size2();
  gaussPts.resize(3, nb_gauss_pts, false);
  gaussPts.clear();

  if (nb_gauss_pts == 0)
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
            "Error No Gauss points created on face side element on Volume ");

  for (int gg = 0; gg != nb_gauss_pts; ++gg) {
    if (side_of_vol_number == 3) {
      gaussPts(0, gg) =
          contact_prism_ptr_fe->getGaussPtsMasterFromEleSide()(0, gg);
      gaussPts(1, gg) =
          contact_prism_ptr_fe->getGaussPtsMasterFromEleSide()(1, gg);
      gaussPts(2, gg) =
          contact_prism_ptr_fe->getGaussPtsMasterFromEleSide()(2, gg);
    } else {
      gaussPts(0, gg) =
          contact_prism_ptr_fe->getGaussPtsSlaveFromEleSide()(0, gg);
      gaussPts(1, gg) =
          contact_prism_ptr_fe->getGaussPtsSlaveFromEleSide()(1, gg);
      gaussPts(2, gg) =
          contact_prism_ptr_fe->getGaussPtsSlaveFromEleSide()(2, gg);
    }
  }

  
    MoFEMFunctionReturn(0);
}
} // namespace MoFEM
