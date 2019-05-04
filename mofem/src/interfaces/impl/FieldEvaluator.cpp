/** \file FieldEvaluator.cpp
 * \brief TetGen interface for meshing/re-meshing and on the fly mesh creation
 *
 */

/**
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

namespace MoFEM {

MoFEMErrorCode
FieldEvaluatorInterface::query_interface(const MOFEMuuid &uuid,
                                         UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMFieldEvaluator) {
    *iface = const_cast<FieldEvaluatorInterface *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
FieldEvaluatorInterface::buildTree3D(const std::string finite_element) {
  CoreInterface &m_field = cOre;
  MoFEMFunctionBegin;
  EntityHandle fe_meshset;
  fe_meshset = m_field.get_finite_element_meshset(finite_element);
  Range entities_3d;
  CHKERR m_field.get_moab().get_entities_by_dimension(fe_meshset, 3,
                                                      entities_3d, true);
  treePtr = boost::make_shared<BVHTree>(&m_field.get_moab());
  CHKERR treePtr->build_tree(entities_3d, &rooTreeSet);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldEvaluatorInterface::SetGaussPts::
operator()(int order_row, int order_col, int order_data) {
  MoFEMFunctionBegin;

  auto shape_check = [&](const auto &t_shape) {
    if (

        t_shape(0) >= 0 - eps && t_shape(0) <= 1 + eps &&

        t_shape(1) >= 0 - eps && t_shape(1) <= 1 + eps &&

        t_shape(2) >= 0 - eps && t_shape(2) <= 1 + eps &&

        t_shape(3) >= 0 - eps && t_shape(3) <= 1 + eps

    )
      return true;
    else
      return false;
  };

  if (verb >= VERY_NOISY)
    std::cout << std::endl << "Next" << std::endl;

  const auto &elem_coords =
      static_cast<VolumeElementForcesAndSourcesCore &>(*feMethod).coords;

  if (verb >= VERY_NOISY)
    std::cout << "elem coords: " << elem_coords << std::endl;

  if (verb >= VERY_NOISY)
    std::cout << "nbEvalPoints " << nbEvalPoints << std::endl;

  CHKERR Tools::getLocalCoordinatesOnReferenceFourNodeTet(
      &elem_coords[0], evalPoints, nbEvalPoints,
      &*localCoords.data().begin());

  if (verb >= VERY_NOISY)
    std::cout << "local_coords: " << localCoords << endl;

  if (verb >= VERY_NOISY)
    std::cout << "shape: " << shapeFunctions << endl;

  MatrixDouble &gauss_pts = feMethod->gaussPts;
  gauss_pts.resize(4, nbEvalPoints, false);
  gauss_pts.clear();

  CHKERR Tools::shapeFunMBTET<3>(&shapeFunctions(0, 0), &localCoords(0, 0),
                                 &localCoords(0, 1), &localCoords(0, 2),
                                 nbEvalPoints);

  FTensor::Tensor1<FTensor::PackPtr<double *, 4>, 4> t_shape = {
      &shapeFunctions(0, 0), &shapeFunctions(0, 1), &shapeFunctions(0, 2),
      &shapeFunctions(0, 3)};
  FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3> t_gauss_pts = {
      &gauss_pts(0, 0), &gauss_pts(1, 0), &gauss_pts(2, 0)};

  int nb_gauss_pts = 0;
  for (int nn = 0; nn != nbEvalPoints; ++nn) {
    if (shape_check(t_shape)) {
      for (const int i : {0, 1, 2}) {
        t_gauss_pts(i) = t_shape(i + 1);
        gauss_pts(3, nb_gauss_pts) = nn;
      }
      ++t_gauss_pts;
      ++nb_gauss_pts;
    }
    ++t_shape;
  }

  if (verb >= VERY_NOISY)
    std::cout << "nbEvalOPoints / nbGaussPts: " << nbEvalPoints << " / "
              << nb_gauss_pts << std::endl;
  gauss_pts.resize(4, nb_gauss_pts, true);
  static_cast<VolumeElementForcesAndSourcesCore &>(*feMethod).nbGaussPts =
      nb_gauss_pts;

  if (verb >= VERY_NOISY)
    std::cout << "gauss pts: " << gauss_pts << std::endl;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldEvaluatorInterface::evalFEAtThePoint3D(
    const double *const point, const double distance, const std::string problem,
    const std::string finite_element,
    boost::shared_ptr<MoFEM::ForcesAndSourcesCore> fe_method, int lower_rank,
    int upper_rank, MoFEMTypes bh, VERBOSITY_LEVELS verb) {
  CoreInterface &m_field = cOre;
  MoFEMFunctionBegin;

  std::vector<EntityHandle> leafs_out;
  CHKERR treePtr->distance_search(point, distance, leafs_out);
  Range tree_ents;
  for (auto lit : leafs_out)
    CHKERR m_field.get_moab().get_entities_by_dimension(lit, 3, tree_ents,
                                                        false);

  if (verb >= VERY_NOISY)
    std::cout << "tree entities: " << tree_ents << endl;


  const Problem *prb_ptr;
  CHKERR m_field.get_problem(problem, &prb_ptr);
  boost::shared_ptr<NumeredEntFiniteElement_multiIndex> numered_fes(
      new NumeredEntFiniteElement_multiIndex());

  for (auto peit = tree_ents.pair_begin(); peit != tree_ents.pair_end();
       ++peit) {
    auto lo =
        prb_ptr->numeredFiniteElements->get<Composite_Name_And_Ent_mi_tag>()
            .lower_bound(boost::make_tuple(finite_element, peit->first));
    auto hi =
        prb_ptr->numeredFiniteElements->get<Composite_Name_And_Ent_mi_tag>()
            .lower_bound(boost::make_tuple(finite_element, peit->second));
    numered_fes->insert(lo, hi);

    if (verb >= VERY_NOISY)
      std::cout << "numered elements:" << std::endl;
    for (; lo != hi; ++lo)
      if (verb >= VERY_NOISY)
        std::cout << **lo << endl;
  }
  if (verb >= VERY_NOISY)
    std::cout << std::endl;

  CHKERR m_field.loop_finite_elements(prb_ptr, finite_element, *fe_method,
                                      lower_rank, upper_rank, numered_fes, bh,
                                      verb);

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM