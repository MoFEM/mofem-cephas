/** \file EdgeElementForcesAndSourcesCoreOnParent.cpp

\brief Implementation of edge element for parent

*/



namespace MoFEM {

int EdgeElementForcesAndSourcesCoreOnChildParent::getRule(int order) {
  return -1;
};

MoFEMErrorCode
EdgeElementForcesAndSourcesCoreOnChildParent::setGaussPts(int order) {
  MoFEMFunctionBegin;

  auto ref_fe = refinePtrFE;
  if (ref_fe == nullptr)
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
            "Pointer to face element is not set");

  const auto ref_ent = ref_fe->getFEEntityHandle();
  const auto ent = getFEEntityHandle();

  const auto &ref_gauss_pts = ref_fe->gaussPts;
  const auto nb_integration_points = ref_gauss_pts.size2();

  auto set_integration_pts_for_edge = [&]() {
    MoFEMFunctionBeginHot;

    auto get_coords = [&](const auto ent) {
      int num_nodes;
      const EntityHandle *conn;
      CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
      std::array<double, 6> node_coords;
      CHKERR mField.get_moab().get_coords(conn, num_nodes, node_coords.data());
      return node_coords;
    };

    auto ref_node_coords = get_coords(ref_ent);
    auto node_coords = get_coords(ent);

    MatrixDouble ref_shapes(nb_integration_points, 2);
    CHKERR Tools::shapeFunMBEDGE<1>(&ref_shapes(0, 0), &ref_gauss_pts(0, 0),
                                    nb_integration_points);

#ifndef NDEBUG
    if (ref_shapes.size1() * ref_shapes.size2() !=
        nb_integration_points * ref_node_coords.size() / 3)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong number of shape functions");
#endif

    auto get_glob_coords = [&]() {
      MatrixDouble glob_coords(nb_integration_points, 3);

      auto t_glob_coords = FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>{
          &glob_coords(0, 0), &glob_coords(0, 1), &glob_coords(0, 2)};
      auto t_shape =
          FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(&ref_shapes(0, 0));

      for (auto gg = 0; gg != nb_integration_points; ++gg) {
        FTensor::Index<'i', 3> i;

        auto t_ref_node_coords =
            FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>{
                &ref_node_coords[0], &ref_node_coords[1], &ref_node_coords[2]};

        t_glob_coords(i) = 0;
        for (int nn = 0; nn != ref_node_coords.size() / 3; ++nn) {
          t_glob_coords(i) += t_shape * t_ref_node_coords(i);
          ++t_ref_node_coords;
          ++t_shape;
        }

        ++t_glob_coords;
      }
      return glob_coords;
    };

    auto glob_coords = get_glob_coords();
    VectorDouble local_coords(nb_integration_points);

    CHKERR Tools::getLocalCoordinatesOnReferenceEdgeNodeEdge(
        node_coords.data(), &glob_coords(0, 0), nb_integration_points,
        &local_coords[0]);

    gaussPts.resize(2, nb_integration_points, false);
    auto t_gauss_pts = FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 2>{
        &gaussPts(0, 0), &gaussPts(1, 0)};
    auto t_local_coords =
        FTensor::Tensor0<FTensor::PackPtr<double *, 1>>{&local_coords[0]};

    for (auto gg = 0; gg != nb_integration_points; ++gg) {
      t_gauss_pts(0) = t_local_coords;
      t_gauss_pts(1) = ref_gauss_pts(1, 0);
      ++t_gauss_pts;
      ++t_local_coords;
    }
    MoFEMFunctionReturnHot(0);
  };

  const auto type = numeredEntFiniteElementPtr->getEntType();

  switch (type) {
  case MBEDGE:
    CHKERR set_integration_pts_for_edge();
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
             "Element type not implemented: %d", type);
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM