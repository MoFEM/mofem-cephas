/** \file FieldEvaluator.cpp
 * \brief Field evaluator
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

FieldEvaluatorInterface::FieldEvaluatorInterface(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)) {

  if (!LogManager::checkIfChannelExist("FieldEvaluatorWorld")) {
    auto core_log = logging::core::get();

    core_log->add_sink(LogManager::createSink(LogManager::getStrmWorld(),
                                              "FieldEvaluatorWorld"));
    core_log->add_sink(LogManager::createSink(LogManager::getStrmSync(),
                                              "FieldEvaluatorSync"));
    core_log->add_sink(LogManager::createSink(LogManager::getStrmSelf(),
                                              "FieldEvaluatorSelf"));

    LogManager::setLog("FieldEvaluatorWorld");
    LogManager::setLog("FieldEvaluatorSync");
    LogManager::setLog("FieldEvaluatorSelf");

    MOFEM_LOG_TAG("FieldEvaluatorWorld", "FieldEvaluator");
    MOFEM_LOG_TAG("FieldEvaluatorSync", "FieldEvaluator");
    MOFEM_LOG_TAG("FieldEvaluatorSelf", "FieldEvaluator");
  }

  MOFEM_LOG("FieldEvaluatorWorld", Sev::noisy) << "Field evaluator intreface";
}

MoFEMErrorCode FieldEvaluatorInterface::query_interface(
    boost::typeindex::type_index type_index, UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = const_cast<FieldEvaluatorInterface *>(this);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
FieldEvaluatorInterface::buildTree3D(boost::shared_ptr<SetPtsData> spd_ptr,
                                     const std::string finite_element) {
  CoreInterface &m_field = cOre;
  MoFEMFunctionBegin;
  EntityHandle fe_meshset;
  fe_meshset = m_field.get_finite_element_meshset(finite_element);
  Range entities_3d;
  CHKERR m_field.get_moab().get_entities_by_dimension(fe_meshset, 3,
                                                      entities_3d, true);
  spd_ptr->treePtr.reset(new AdaptiveKDTree(&m_field.get_moab()));
  CHKERR spd_ptr->treePtr->build_tree(entities_3d, &spd_ptr->rooTreeSet);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
FieldEvaluatorInterface::SetPts::operator()(ForcesAndSourcesCore *fe_raw_ptr,
                                            int order_row, int order_col,
                                            int order_data) {
  MoFEMFunctionBegin;

  if (auto data_ptr = dataPtr.lock()) {

    decltype(data_ptr->nbEvalPoints) nbEvalPoints = data_ptr->nbEvalPoints;
    decltype(data_ptr->verb) verb = data_ptr->verb;

    decltype(data_ptr->shapeFunctions) &shapeFunctions =
        data_ptr->shapeFunctions;
    decltype(data_ptr->evalPointEntityHandle) &evalPointEntityHandle =
        data_ptr->evalPointEntityHandle;

    if (auto fe_ptr = data_ptr->feMethodPtr.lock()) {

      VolumeElementForcesAndSourcesCore &fe =
          static_cast<VolumeElementForcesAndSourcesCore &>(*fe_ptr);

      if (fe_ptr.get() != fe_raw_ptr)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Wrong FE ptr");

      MatrixDouble &gauss_pts = fe.gaussPts;
      gauss_pts.resize(4, nbEvalPoints, false);
      gauss_pts.clear();

      FTensor::Tensor1<FTensor::PackPtr<double *, 4>, 4> t_shape = {
          &shapeFunctions(0, 0), &shapeFunctions(0, 1), &shapeFunctions(0, 2),
          &shapeFunctions(0, 3)};
      FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3> t_gauss_pts = {
          &gauss_pts(0, 0), &gauss_pts(1, 0), &gauss_pts(2, 0)};

      int nb_gauss_pts = 0;
      for (int nn = 0; nn != nbEvalPoints; ++nn) {
        if (evalPointEntityHandle[nn] ==
            fe.numeredEntFiniteElementPtr->getEnt()) {
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

      if (verb >= VERY_NOISY)
        std::cout << "gauss pts: " << gauss_pts << std::endl;

    } else
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Pointer to element does not exists");

  } else
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer to data does not exists");

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldEvaluatorInterface::evalFEAtThePoint3D(
    const double *const point, const double distance, const std::string problem,
    const std::string finite_element, boost::shared_ptr<SetPtsData> data_ptr,
    int lower_rank, int upper_rank, boost::shared_ptr<CacheTuple> cache_ptr,
    MoFEMTypes bh, VERBOSITY_LEVELS verb) {
  CoreInterface &m_field = cOre;
  MoFEMFunctionBegin;

  std::vector<EntityHandle> leafs_out;
  CHKERR data_ptr->treePtr->distance_search(point, distance, leafs_out);
  Range tree_ents;
  for (auto lit : leafs_out)
    CHKERR m_field.get_moab().get_entities_by_dimension(lit, 3, tree_ents,
                                                        true);

  if (verb >= VERY_NOISY)
    std::cout << "tree entities: " << tree_ents << endl;

  data_ptr->evalPointEntityHandle.resize(data_ptr->nbEvalPoints);
  std::fill(data_ptr->evalPointEntityHandle.begin(),
            data_ptr->evalPointEntityHandle.end(), 0);

  for (auto tet : tree_ents) {

    const EntityHandle *conn;
    int num_nodes;
    CHKERR m_field.get_moab().get_connectivity(tet, conn, num_nodes, true);
    std::array<double, 12> coords;
    CHKERR m_field.get_moab().get_coords(conn, num_nodes, coords.data());

    for (int n = 0; n != data_ptr->nbEvalPoints; ++n) {

      std::array<double, 3> local_coords;
      CHKERR Tools::getLocalCoordinatesOnReferenceFourNodeTet(
          coords.data(), &data_ptr->evalPoints[3 * n], 1, local_coords.data());

      std::array<double, 4> shape;
      CHKERR Tools::shapeFunMBTET<3>(shape.data(), &local_coords[0],
                                     &local_coords[1], &local_coords[2], 1);

      const double eps = data_ptr->eps;
      if (shape[0] >= 0 - eps && shape[0] <= 1 + eps &&

          shape[1] >= 0 - eps && shape[1] <= 1 + eps &&

          shape[2] >= 0 - eps && shape[2] <= 1 + eps &&

          shape[3] >= 0 - eps && shape[3] <= 1 + eps) {

        std::copy(shape.begin(), shape.end(), &data_ptr->shapeFunctions(n, 0));
        std::copy(local_coords.begin(), local_coords.end(),
                  &data_ptr->localCoords(n, 0));
        data_ptr->evalPointEntityHandle[n] = tet;
      }
    }
  }

  const Problem *prb_ptr;
  CHKERR m_field.get_problem(problem, &prb_ptr);
  boost::shared_ptr<NumeredEntFiniteElement_multiIndex> numered_fes(
      new NumeredEntFiniteElement_multiIndex());

  Range in_tets;
  in_tets.insert_list(data_ptr->evalPointEntityHandle.begin(),
                      data_ptr->evalPointEntityHandle.end());
  in_tets = in_tets.subset_by_dimension(3);

  if (verb >= VERY_NOISY)
    std::cout << "in tets: " << in_tets << endl;

  for (auto peit = in_tets.pair_begin(); peit != in_tets.pair_end(); ++peit) {
    auto lo =
        prb_ptr->numeredFiniteElementsPtr->get<Composite_Name_And_Ent_mi_tag>()
            .lower_bound(boost::make_tuple(finite_element, peit->first));
    auto hi =
        prb_ptr->numeredFiniteElementsPtr->get<Composite_Name_And_Ent_mi_tag>()
            .upper_bound(boost::make_tuple(finite_element, peit->second));
    numered_fes->insert(lo, hi);

    if (verb >= VERY_NOISY)
      std::cout << "numered elements:" << std::endl;
    for (; lo != hi; ++lo)
      if (verb >= VERY_NOISY)
        std::cout << **lo << endl;
  }
  if (verb >= VERY_NOISY)
    std::cout << std::endl;

  if (auto fe_ptr = data_ptr->feMethodPtr.lock()) {

    MOFEM_LOG_C("FieldEvaluatorSync", Sev::verbose,
                "Number elements %d to evaluate at proc %d",
                numered_fes->size(), m_field.get_comm_rank());
    MOFEM_LOG_SYNCHRONISE(m_field.get_comm());

    if(!cache_ptr) {
      cache_ptr = boost::make_shared<CacheTuple>();
      CHKERR m_field.cache_problem_entities(prb_ptr->getName(), cache_ptr);

      MOFEM_LOG("FieldEvaluatorSync", Sev::noisy)
          << "If you call that function many times in the loop consider to set "
             "cache_ptr outside of the loop. Otherwise code can be slow.";
    }

    CHKERR m_field.loop_finite_elements(prb_ptr, finite_element, *fe_ptr,
                                        lower_rank, upper_rank, numered_fes, bh,
                                        cache_ptr, verb);

    } else
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer to element does not exists");

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM