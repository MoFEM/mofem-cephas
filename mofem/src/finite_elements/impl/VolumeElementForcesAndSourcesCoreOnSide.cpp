/** \file VolumeElementForcesAndSourcesCoreOnSide.cpp

\brief Implementation of volume element on side

*/



namespace MoFEM {

int VolumeElementForcesAndSourcesCoreOnSide::getRule(int order) {
  return -1;
};

MoFEMErrorCode
VolumeElementForcesAndSourcesCoreOnSide::setGaussPts(int order) {
  MoFEMFunctionBegin;

  if (!sidePtrFE)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Side element not set");

  const auto type = numeredEntFiniteElementPtr->getEntType();
  const auto nb_nodes_on_ele = CN::VerticesPerEntity(type);

  auto face_ptr_fe =
      static_cast<FaceElementForcesAndSourcesCore *>(sidePtrFE);

  const auto face_entity = sidePtrFE->numeredEntFiniteElementPtr->getEnt();
  auto &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());
  auto sit = side_table.get<0>().find(face_entity);
  if (sit == side_table.get<0>().end())
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
            "Face can not be found on volume element");

  const auto nb_nodes_on_face = CN::VerticesPerEntity((*sit)->getEntType());

  faceSense = (*sit)->sense;
  faceSideNumber = (*sit)->side_number;
  fill(&tetConnMap[0], &tetConnMap[nb_nodes_on_ele], -1);
  for (int nn = 0; nn != nb_nodes_on_face; ++nn) {
    faceConnMap[nn] = std::distance(
        conn, find(conn, &conn[nb_nodes_on_ele], face_ptr_fe->conn[nn]));
    tetConnMap[faceConnMap[nn]] = nn;
    if (faceConnMap[nn] >= nb_nodes_on_ele)
      SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
              "No common node on face and element can not be found");
  }

  oppositeNode = std::distance(
      &tetConnMap[0], find(&tetConnMap[0], &tetConnMap[nb_nodes_on_ele], -1));

  const int nb_gauss_pts = face_ptr_fe->gaussPts.size2();
  gaussPts.resize(4, nb_gauss_pts, false);
  gaussPts.clear();

  constexpr double tet_coords[12] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
  constexpr double hex_coords[24] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                     0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};

  auto set_gauss = [&](auto &coords) {
    MoFEMFunctionBegin;
    FTensor::Index<'i', 3> i;
    auto &data = *face_ptr_fe->dataOnElement[H1];
    auto t_base = data.dataOnEntities[MBVERTEX][0].getFTensor0N(NOBASE);
    FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3> t_gauss_coords{
        &gaussPts(0, 0), &gaussPts(1, 0), &gaussPts(2, 0)};
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      for (int bb = 0; bb != nb_nodes_on_face; ++bb) {
        const int shift = 3 * faceConnMap[bb];
        FTensor::Tensor1<const double, 3> t_coords{
            coords[shift + 0], coords[shift + 1], coords[shift + 2]};
        t_gauss_coords(i) += t_base * t_coords(i);
        ++t_base;
      }
      ++t_gauss_coords;
      gaussPts(3, gg) = face_ptr_fe->gaussPts(2, gg);
    }
    MoFEMFunctionReturn(0);
  };

  switch (type) {
  case MBTET:
    CHKERR set_gauss(tet_coords);
    break;
  case MBHEX:
    CHKERR set_gauss(hex_coords);
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
             "Element type not implemented: %d", type);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::setPtrFE(
    ForcesAndSourcesCore *ptr) {
  MoFEMFunctionBeginHot;
  if (!(ptrFE =
            dynamic_cast<VolumeElementForcesAndSourcesCoreOnSide *>(ptr)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "User operator and finite element do not work together");
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM
