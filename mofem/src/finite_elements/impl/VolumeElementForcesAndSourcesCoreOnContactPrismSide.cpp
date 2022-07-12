/** \file VolumeElementForcesAndSourcesCoreOnContactPrismSide.cpp

\brief Implementation of volume element on side

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

namespace MoFEM {

int VolumeElementForcesAndSourcesCoreOnContactPrismSide::getRule(int order) {
  return -1;
};

MoFEMErrorCode
VolumeElementForcesAndSourcesCoreOnContactPrismSide::setGaussPts(int order) {
  MoFEMFunctionBegin;

  if (!sidePtrFE)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Side volume element not set");

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

  fill(tetConnMap.begin(), tetConnMap.end(), -1);
  for (int nn = 0; nn != 3; ++nn) {
    faceConnMap[nn] = std::distance(conn, find(conn, &conn[4], ent_on_vol[nn]));

    tetConnMap[faceConnMap[nn]] = nn;
    if (faceConnMap[nn] > 3) {
      cerr << " faceConnMap[nn] " << faceConnMap[nn] << "\n";
      SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
              "No common node on face and element can not be found");
    }
  }

  oppositeNode = std::distance(tetConnMap.begin(),
                               find(tetConnMap.begin(), tetConnMap.end(), -1));

  const int nb_gauss_pts =
      contact_prism_ptr_fe->getGaussPtsMasterFromEleSide().size2();
  gaussPts.resize(4, nb_gauss_pts, false);
  gaussPts.clear();

  const EntityType tri_type = MBTRI;
  boost::shared_ptr<EntitiesFieldData> dataH1_on_face;

  if (side_of_vol_number == 3) {
    dataH1_on_face = boost::shared_ptr<EntitiesFieldData>(
        contact_prism_ptr_fe->getDataOnMasterFromEleSide()[H1]);
  } else if (side_of_vol_number == 4) {
    dataH1_on_face = boost::shared_ptr<EntitiesFieldData>(
        contact_prism_ptr_fe->getDataOnSlaveFromEleSide()[H1]);
  } else {
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
            "Side element for contact: Wrong face for contact!");
  }

  const MatrixDouble &face_shape_funtions =
      dataH1_on_face->dataOnEntities[MBVERTEX][0].getN(NOBASE);

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
    if (side_of_vol_number == 3) {
      gaussPts(3, gg) =
          contact_prism_ptr_fe->getGaussPtsMasterFromEleSide()(2, gg);
    } else {
      gaussPts(3, gg) =
          contact_prism_ptr_fe->getGaussPtsSlaveFromEleSide()(2, gg);
    }
  }
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
