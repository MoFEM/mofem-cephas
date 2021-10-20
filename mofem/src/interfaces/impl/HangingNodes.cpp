/** \file HangingNodes.hpp
 * \brief Refine mesh
 *
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

HangingNodes::HangingNodes(const Core &core) : cOre(const_cast<Core &>(core)) {

  if (!LogManager::checkIfChannelExist("HangingNodesSelf")) {
    auto core_log = logging::core::get();

    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSelf(), "HangingNodesWorld"));
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSelf(), "HangingNodesSelf"));

    LogManager::setLog("HangingNodesWOrld");
    LogManager::setLog("HangingNodesSelf");

    MOFEM_LOG_TAG("HangingNodesWorld", "HangingWorld");
    MOFEM_LOG_TAG("HangingNodesSelf", "HangingNodes");
  }

  MOFEM_LOG("HangingNodesWorld", Sev::noisy) << "Hanging nodes interface";
}

MoFEMErrorCoda HangingNodes::addVerticesInTheMiddleOfEdges(
    const Range &edges, const BitRefLevel &bit, EntityHandle start_v = 0) {

  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  auto refined_ents_ptr = m_field.get_ref_ents();

  MoFEMFunctionBegin;
  auto miit =
      refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>()
          .lower_bound(boost::make_tuple(MBVERTEX, MBEDGE));
  auto hi_miit =
      refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>()
          .upper_bound(boost::make_tuple(MBVERTEX, MBEDGE));
  RefEntity_multiIndex_view_by_ordered_parent_entity ref_parent_ents_view;
  ref_parent_ents_view.insert(miit, hi_miit);

  Range edges = ents.subset_by_type(MBEDGE);
  MOFEM_LOG("HangingNodesSelf", Sev::verbose)
      << "ref level " << bit << " nb. edges to refine " << edges.size();

  std::array<std::vector<double>, 3> vert_coords;
  for (auto &vc : vert_coords)
    vc.reserve(edges.size());

  std::vector<EntityHandle> parent_edge;
  parent_edge.reserve(edges.size());

  std::array<double, 6> coords;
  FTensor::Tensor1<FTensor::PackPtr<double *, 0>, 3> t_0 = {
      &coords[0], &coords[1], &coords[2]};
  FTensor::Tensor1<FTensor::PackPtr<double *, 0>, 3> t_1 = {
      &coords[3], &coords[4], &coords[5]};

  FTensor::Index<'i', 3> i;

  Range add_bit;
  for (auto p_eit = edges.pair_begin(); p_eit != edges.pair_end(); ++p_eit) {
    auto miit_view = ref_parent_ents_view.lower_bound(p_eit->first);

    Range edge_having_parent_vertex;
    if (miit_view != ref_parent_ents_view.end()) {
      auto hi_miit_view = ref_parent_ents_view.upper_bound(p_eit->second);
      for (; miit_view != hi_miit_view; ++miit_view) {
        edge_having_parent_vertex.insert(edge_having_parent_vertex.end(),
                                         miit_view->get()->getParentEnt());
        add_bit.insert(add_bit.end(), miit_view->get()->getEnt());
      }
    }

    Range add_vertex_edges =
        subtract(Range(p_eit->first, p_eit->second), edge_having_parent_vertex);

    for (auto e : add_vertex_edges)
      parent_edge.emplace_back(e);

    auto e_it = add_vertex_edges.begin();
    while (eit != add_vertex_edges) {
      EntityHandle *conn_ptr;
      int vpere, count;
      CHKERR moab.connect_iterate(e_it, add_vertex_edges.end(), conn_ptr, vpere,
                                  count);
      for (int i = 0; i != count; ++i) {
        CHKERR moab.get_coords(conn_ptr, vpere, coords.data());
        t_0(i) += t_1(i);
        t_0(i) *= 0.5;
        conn_ptr += vpere;
      }
      e_eit += count;
    }
  }

  CHKERR m_field.getInterface<BitRefManager>()->addBitRefLevel(add_bit, bit);

  if (!vert_coords[0].empty()) {
    ReadUtilIface *read_util;
    CHKERR moab.query_interface(read_util);
    int num_nodes = vert_coords[0].size();
    vector<double *> arrays_coord;
    CHKERR read_util->get_node_coords(3, num_nodes, 0, start_v, arrays_coord);
    Range verts(start_v, start_v + num_nodes - 1);
    for (auto dd : {0, 1, 2}) {
      std::copy(vert_coords[dd].begin(), vert_coords[dd].end(),
                arrays_coord[dd]);
    }
    CHKERR moab.tag_set_data(cOre.get_th_RefParentHandle(), verts,
                             &*parent_edge.begin());
    CHKERR m_field.getInterface<BitRefManager>()->setEntitiesBitRefLevel(
        verts, bit, verb);
  }

  MoFEMFunctionReturn(0);
}