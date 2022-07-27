/** \file FieldEvaluator.cpp
 * \brief Field evaluator
 *
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

template <int D>
MoFEMErrorCode
FieldEvaluatorInterface::buildTree(boost::shared_ptr<SetPtsData> spd_ptr,
                                   const std::string finite_element) {
  CoreInterface &m_field = cOre;
  MoFEMFunctionBegin;
  EntityHandle fe_meshset;
  fe_meshset = m_field.get_finite_element_meshset(finite_element);
  Range entities;
  CHKERR m_field.get_moab().get_entities_by_dimension(fe_meshset, D, entities,
                                                      true);
  spd_ptr->treePtr.reset(new AdaptiveKDTree(&m_field.get_moab()));
  CHKERR spd_ptr->treePtr->build_tree(entities, &spd_ptr->rooTreeSet);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
FieldEvaluatorInterface::buildTree3D(boost::shared_ptr<SetPtsData> spd_ptr,
                                     const std::string finite_element) {
  return buildTree<3>(spd_ptr, finite_element);
}

MoFEMErrorCode
FieldEvaluatorInterface::buildTree2D(boost::shared_ptr<SetPtsData> spd_ptr,
                                     const std::string finite_element) {
  return buildTree<2>(spd_ptr, finite_element);
}

MoFEMErrorCode
FieldEvaluatorInterface::SetPts::operator()(ForcesAndSourcesCore *fe_raw_ptr,
                                            int order_row, int order_col,
                                            int order_data) {
  MoFEMFunctionBegin;

  if (auto data_ptr = dataPtr.lock()) {

    const auto nb_eval_points = data_ptr->nbEvalPoints;
    const auto &shape_functions = data_ptr->shapeFunctions;
    const auto &eval_pointentity_handle = data_ptr->evalPointEntityHandle;

    if (auto fe_ptr = data_ptr->feMethodPtr.lock()) {

      auto &fe = static_cast<ForcesAndSourcesCore &>(*fe_ptr);

#ifndef NDEBUG
      if (fe_ptr.get() != fe_raw_ptr)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Wrong FE ptr");
#endif

      const auto fe_ent = fe.numeredEntFiniteElementPtr->getEnt();
      const auto fe_type = type_from_handle(fe_ent);
      const auto fe_dim = moab::CN::Dimension(fe_type);

      auto &gauss_pts = fe.gaussPts;
      int nb_gauss_pts;

      if (fe_dim == 3) {
        gauss_pts.resize(4, nb_eval_points, false);
        FTensor::Tensor1<FTensor::PackPtr<const double *, 4>, 4> t_shape{
            &shape_functions(0, 0), &shape_functions(0, 1),
            &shape_functions(0, 2), &shape_functions(0, 3)};
        FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3> t_gauss_pts{
            &gauss_pts(0, 0), &gauss_pts(1, 0), &gauss_pts(2, 0)};

        nb_gauss_pts = 0;
        for (int nn = 0; nn != nb_eval_points; ++nn) {
          if (eval_pointentity_handle[nn] == fe_ent) {
            for (const int i : {0, 1, 2}) {
              t_gauss_pts(i) = t_shape(i + 1);
            }
            gauss_pts(3, nb_gauss_pts) = nn;
            ++t_gauss_pts;
            ++nb_gauss_pts;
          }
          ++t_shape;
        }
        gauss_pts.resize(4, nb_gauss_pts, true);
      } else if (fe_dim == 2) {
        gauss_pts.resize(3, nb_eval_points, false);
        FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3> t_shape{
            &shape_functions(0, 0), &shape_functions(0, 1),
            &shape_functions(0, 2)};
        FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 2> t_gauss_pts{
            &gauss_pts(0, 0), &gauss_pts(1, 0)};
        nb_gauss_pts = 0;
        for (int nn = 0; nn != nb_eval_points; ++nn) {
          if (eval_pointentity_handle[nn] == fe_ent) {
            for (const int i : {0, 1}) {
              t_gauss_pts(i) = t_shape(i + 1);
            }
            gauss_pts(2, nb_gauss_pts) = nn;
            ++t_gauss_pts;
            ++nb_gauss_pts;
          }
          ++t_shape;
        }
        gauss_pts.resize(3, nb_gauss_pts, true);
      } else {
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                "Dimension not implemented");
      }

#ifndef NDEBUG
      MOFEM_LOG("SELF", Sev::noisy)
          << "nbEvalOPoints / nbGau_sspt_s: " << nb_eval_points << " / "
          << nb_gauss_pts;
      MOFEM_LOG("SELF", Sev::noisy) << "gauss pts: " << gauss_pts;
#endif

    } else
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Pointer to element does not exists");

  } else
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer to data does not exists");

  MoFEMFunctionReturn(0);
}

template <int D>
MoFEMErrorCode FieldEvaluatorInterface::evalFEAtThePoint(
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
    CHKERR m_field.get_moab().get_entities_by_dimension(lit, D, tree_ents,
                                                        true);
  if (verb >= NOISY)
    MOFEM_LOG("SELF", Sev::noisy) << "tree entities: " << tree_ents;

  data_ptr->evalPointEntityHandle.resize(data_ptr->nbEvalPoints);
  std::fill(data_ptr->evalPointEntityHandle.begin(),
            data_ptr->evalPointEntityHandle.end(), 0);

  std::vector<double> local_coords;
  std::vector<double> shape;

  auto nb_eval_points = data_ptr->nbEvalPoints;

  for (auto tet : tree_ents) {
    const EntityHandle *conn;
    int num_nodes;
    CHKERR m_field.get_moab().get_connectivity(tet, conn, num_nodes, true);

    if constexpr (D == 3) {

      local_coords.resize(3 * nb_eval_points);
      shape.resize(4 * nb_eval_points);

      std::array<double, 12> coords;
      CHKERR m_field.get_moab().get_coords(conn, num_nodes, coords.data());

      CHKERR Tools::getLocalCoordinatesOnReferenceFourNodeTet(
          coords.data(), &data_ptr->evalPoints[0], nb_eval_points,
          &local_coords[0]);
      CHKERR Tools::shapeFunMBTET<3>(&shape[0], &local_coords[0],
                                     &local_coords[1], &local_coords[2],
                                     nb_eval_points);

      FTensor::Index<'i', 4> i4;
      FTensor::Tensor1<FTensor::PackPtr<double *, 4>, 4> t_shape{
          &shape[0], &shape[1], &shape[2], &shape[3]};
      FTensor::Tensor1<FTensor::PackPtr<double *, 4>, 4> t_shape_data{
          &data_ptr->shapeFunctions(0, 0), &data_ptr->shapeFunctions(0, 1),
          &data_ptr->shapeFunctions(0, 2), &data_ptr->shapeFunctions(0, 3)};

      FTensor::Index<'j', 3> j3;
      FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_local{
          &local_coords[0], &local_coords[1], &local_coords[2]};
      FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_local_data{
          &(data_ptr->localCoords(0, 0)), &(data_ptr->localCoords(0, 1)),
          &(data_ptr->localCoords(0, 2))};


      for (int n = 0; n != nb_eval_points; ++n) {

        const double eps = data_ptr->eps;
        if (t_shape(0) >= 0 - eps && t_shape(0) <= 1 + eps &&

            t_shape(1) >= 0 - eps && t_shape(1) <= 1 + eps &&

            t_shape(2) >= 0 - eps && t_shape(2) <= 1 + eps &&

            t_shape(3) >= 0 - eps && t_shape(3) <= 1 + eps) {

          data_ptr->evalPointEntityHandle[n] = tet;
          t_shape_data(i4) = t_shape(i4);
          t_local_data(j3) = t_local(j3);
        }

        ++t_shape;
        ++t_shape_data;
        ++t_local;
        ++t_local_data; 

      }

    }

    if constexpr (D == 2) {

      local_coords.resize(2 * nb_eval_points);
      shape.resize(3 * nb_eval_points);

      std::array<double, 9> coords;
      CHKERR m_field.get_moab().get_coords(conn, num_nodes, coords.data());

      CHKERR Tools::getLocalCoordinatesOnReferenceTriNodeTri(
          coords.data(), &data_ptr->evalPoints[0], nb_eval_points,
          &local_coords[0]);
      CHKERR Tools::shapeFunMBTRI<2>(&shape[0], &local_coords[0],
                                     &local_coords[1], nb_eval_points);

      FTensor::Index<'i', 3> i3;
      FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_shape{
          &shape[0], &shape[1], &shape[2]};
      FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_shape_data{
          &data_ptr->shapeFunctions(0, 0), &data_ptr->shapeFunctions(0, 1),
          &data_ptr->shapeFunctions(0, 2)};

      FTensor::Index<'j', 2> j2;
      FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2> t_local{
          &local_coords[0], &local_coords[1]};
      data_ptr->localCoords.resize(nb_eval_points, 2);
      FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2> t_local_data{
          &(data_ptr->localCoords(0, 0)), &(data_ptr->localCoords(0, 1))};

      for (int n = 0; n != nb_eval_points; ++n) {

        const double eps = data_ptr->eps;
        if (t_shape(0) >= 0 - eps && t_shape(0) <= 1 + eps &&

            t_shape(1) >= 0 - eps && t_shape(1) <= 1 + eps &&

            t_shape(2) >= 0 - eps && t_shape(2) <= 1 + eps) {

          data_ptr->evalPointEntityHandle[n] = tet;
          t_shape_data(i3) = t_shape(i3);
          t_local_data(j2) = t_local(j2);
        }

        ++t_shape;
        ++t_shape_data;
        ++t_local;
        ++t_local_data;
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
  in_tets = in_tets.subset_by_dimension(D);

  if (verb >= NOISY)
    MOFEM_LOG("SELF", Sev::noisy) << "in tets: " << in_tets << endl;

  for (auto peit = in_tets.pair_begin(); peit != in_tets.pair_end(); ++peit) {
    auto lo =
        prb_ptr->numeredFiniteElementsPtr->get<Composite_Name_And_Ent_mi_tag>()
            .lower_bound(boost::make_tuple(finite_element, peit->first));
    auto hi =
        prb_ptr->numeredFiniteElementsPtr->get<Composite_Name_And_Ent_mi_tag>()
            .upper_bound(boost::make_tuple(finite_element, peit->second));
    numered_fes->insert(lo, hi);

    if (verb >= VERY_NOISY)
      std::cerr << "numered elements:" << std::endl;
    for (; lo != hi; ++lo)
      if (verb >= VERY_NOISY)
        std::cerr << **lo << endl;
  }
  if (verb >= VERY_NOISY)
    std::cerr << std::endl;

  if (auto fe_ptr = data_ptr->feMethodPtr.lock()) {

    MOFEM_LOG_C("FieldEvaluatorSync", Sev::verbose,
                "Number elements %d to evaluate at proc %d",
                numered_fes->size(), m_field.get_comm_rank());
    MOFEM_LOG_SYNCHRONISE(m_field.get_comm());

    if (!cache_ptr) {
      cache_ptr = boost::make_shared<CacheTuple>();
      CHKERR m_field.cache_problem_entities(prb_ptr->getName(), cache_ptr);

      MOFEM_LOG("FieldEvaluatorSync", Sev::noisy)
          << "If you call that function many times in the loop consider to "
             "set "
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

MoFEMErrorCode FieldEvaluatorInterface::evalFEAtThePoint3D(
    const double *const point, const double distance, const std::string problem,
    const std::string finite_element, boost::shared_ptr<SetPtsData> data_ptr,
    int lower_rank, int upper_rank, boost::shared_ptr<CacheTuple> cache_ptr,
    MoFEMTypes bh, VERBOSITY_LEVELS verb) {
  return evalFEAtThePoint<3>(point, distance, problem, finite_element, data_ptr,
                             lower_rank, upper_rank, cache_ptr, bh, verb);
}

MoFEMErrorCode FieldEvaluatorInterface::evalFEAtThePoint2D(
    const double *const point, const double distance, const std::string problem,
    const std::string finite_element, boost::shared_ptr<SetPtsData> data_ptr,
    int lower_rank, int upper_rank, boost::shared_ptr<CacheTuple> cache_ptr,
    MoFEMTypes bh, VERBOSITY_LEVELS verb) {
  return evalFEAtThePoint<2>(point, distance, problem, finite_element, data_ptr,
                             lower_rank, upper_rank, cache_ptr, bh, verb);
}

} // namespace MoFEM