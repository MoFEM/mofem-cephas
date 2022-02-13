/** \file MeshRefinement.cpp
 * \brief Mesh refinement interface
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

#include <EntityRefine.hpp>

namespace MoFEM {

typedef multi_index_container<
    boost::shared_ptr<RefElement>,
    // ptrWrapperRefElement,
    indexed_by<
        ordered_unique<tag<Ent_mi_tag>,
                       const_mem_fun<RefElement::interface_type_RefEntity,
                                     EntityHandle, &RefElement::getEnt>>,
        ordered_non_unique<
            tag<Ent_Ent_mi_tag>,
            const_mem_fun<RefElement::interface_type_RefEntity, EntityHandle,
                          &RefElement::getParentEnt>>,
        ordered_non_unique<
            tag<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>,
            composite_key<
                RefElement,
                const_mem_fun<RefElement::interface_type_RefEntity,
                              EntityHandle, &RefElement::getParentEnt>,
                const_mem_fun<RefElement, int,
                              &RefElement::getBitRefEdgesUlong>>>>>
    RefElement_multiIndex_parents_view;

MoFEMErrorCode
MeshRefinement::query_interface(boost::typeindex::type_index type_index,
                                UnknownInterface **iface) const {
  *iface = const_cast<MeshRefinement *>(this);
  return 0;
}

MeshRefinement::MeshRefinement(const Core &core)
    : cOre(const_cast<Core &>(core)) {}

MoFEMErrorCode MeshRefinement::addVerticesInTheMiddleOfEdges(
    const EntityHandle meshset, const BitRefLevel &bit, const bool recursive,
    int verb, EntityHandle start_v) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  Range edges;
  CHKERR moab.get_entities_by_type(meshset, MBEDGE, edges, recursive);
  if (edges.empty()) {
    Range tets;
    CHKERR moab.get_entities_by_type(meshset, MBTET, tets, recursive);
    CHKERR moab.get_adjacencies(tets, 1, true, edges, moab::Interface::UNION);
    if (tets.empty()) {
      Range prisms;
      CHKERR moab.get_entities_by_type(meshset, MBPRISM, prisms, recursive);
      for (Range::iterator pit = prisms.begin(); pit != prisms.end(); pit++) {
        const EntityHandle *conn;
        int num_nodes;
        CHKERR moab.get_connectivity(*pit, conn, num_nodes, true);
        assert(num_nodes == 6);
        //
        Range edge;
        CHKERR moab.get_adjacencies(&conn[0], 2, 1, true, edge);
        assert(edge.size() == 1);
        edges.insert(edge[0]);
        edge.clear();
        CHKERR moab.get_adjacencies(&conn[1], 2, 1, true, edge);
        assert(edge.size() == 1);
        edges.insert(edge[0]);
        EntityHandle conn_edge2[] = {conn[2], conn[0]};
        edge.clear();
        CHKERR moab.get_adjacencies(conn_edge2, 2, 1, true, edge);
        assert(edge.size() == 1);
        edges.insert(edge[0]);
        //
        edge.clear();
        CHKERR moab.get_adjacencies(&conn[3], 2, 1, true, edge);
        assert(edge.size() == 1);
        edges.insert(edge[0]);
        edge.clear();
        CHKERR moab.get_adjacencies(&conn[4], 2, 1, true, edge);
        assert(edge.size() == 1);
        edges.insert(edge[0]);
        EntityHandle conn_edge8[] = {conn[5], conn[3]};
        edge.clear();
        CHKERR moab.get_adjacencies(conn_edge8, 2, 1, true, edge);
        assert(edge.size() == 1);
        edges.insert(edge[0]);
      }
    }
  }
  CHKERR addVerticesInTheMiddleOfEdges(edges, bit, verb, start_v);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode MeshRefinement::addVerticesInTheMiddleOfEdges(
    const Range &ents, const BitRefLevel &bit, int verb, EntityHandle start_v) {
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
  if (verb >= VERBOSE) {
    MOFEM_TAG_AND_LOG("WORLD", Sev::verbose, "MeshRefinemnt")
        << "ref level " << bit << " nb. edges to refine " << edges.size();
  }

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
      for (auto hi_miit_view = ref_parent_ents_view.upper_bound(p_eit->second);
           miit_view != hi_miit_view; ++miit_view) {
        edge_having_parent_vertex.insert(edge_having_parent_vertex.end(),
                                         miit_view->get()->getParentEnt());
        add_bit.insert(add_bit.end(), miit_view->get()->getEnt());
      }
    }

    Range add_vertex_edges =
        subtract(Range(p_eit->first, p_eit->second), edge_having_parent_vertex);

    for (auto e : add_vertex_edges)
      parent_edge.emplace_back(e);

    for (auto e : add_vertex_edges) {

      const EntityHandle *conn;
      int num_nodes;
      CHKERR moab.get_connectivity(e, conn, num_nodes, true);
      if (PetscUnlikely(num_nodes != 2)) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "edge should have 2 edges");
      }
      CHKERR moab.get_coords(conn, num_nodes, coords.data());
      t_0(i) += t_1(i);
      t_0(i) *= 0.5;

      for (auto j : {0, 1, 2})
        vert_coords[j].emplace_back(t_0(j));
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

MoFEMErrorCode MeshRefinement::refineTets(const EntityHandle meshset,
                                          const BitRefLevel &bit, int verb,
                                          const bool debug) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  Range tets;
  CHKERR moab.get_entities_by_type(meshset, MBTET, tets, false);
  CHKERR refineTets(tets, bit, verb, debug);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshRefinement::SetParent::operator()(
    const EntityHandle ent, const EntityHandle parent,
    const RefEntity_multiIndex *ref_ents_ptr, MoFEM::Core &cOre) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  auto it = ref_ents_ptr->find(ent);
  if (it != ref_ents_ptr->end()) {
    if (it->get()->getParentEnt() != parent && ent != parent) {
      parentsToChange[ent] = parent;
    }
  } else {
    if (ent != parent) {
      CHKERR m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(),
                                             &ent, 1, &parent);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshRefinement::SetParent::operator()(
    const RefEntity_multiIndex *ref_ents_ptr) {
  MoFEMFunctionBegin;
  for (auto mit = parentsToChange.begin(); mit != parentsToChange.end();
       ++mit) {
    auto it = ref_ents_ptr->find(mit->first);
    bool success = const_cast<RefEntity_multiIndex *>(ref_ents_ptr)
                       ->modify(it, RefEntity_change_parent(mit->second));
    if (!success) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "impossible to set parent");
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshRefinement::refineTets(const Range &_tets,
                                          const BitRefLevel &bit, int verb,
                                          const bool debug) {

  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  auto refined_ents_ptr = m_field.get_ref_ents();
  auto refined_finite_elements_ptr = m_field.get_ref_finite_elements();
  ReadUtilIface *read_util;
  MoFEMFunctionBegin;

  CHKERR m_field.get_moab().query_interface(read_util);

  Range ents_to_set_bit;

  // FIXME: refinement is based on entity handlers, should work on global ids of
  // nodes, this will allow parallelise algorithm in the future

  auto get_parent_ents_view = [&](const auto type_child,
                                  const auto type_parent) {
    auto &ref_ents =
        refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>();
    auto miit =
        ref_ents.lower_bound(boost::make_tuple(type_child, type_parent));
    auto hi_miit =
        ref_ents.upper_bound(boost::make_tuple(type_child, type_parent));
    RefEntity_multiIndex_view_by_ordered_parent_entity ref_parent_ents_view;

    decltype(miit) first;
    decltype(miit) second;

    do {

      while (miit != hi_miit) {
        first = miit;
        if (((*miit)->getBitRefLevel() & bit).any()) {
          ++miit;
          second = miit;
          break;
        } else {
          ++miit;
        }
      }

      while (miit != hi_miit) {
        if (((*miit)->getBitRefLevel() & bit).any()) {
          second = miit;
          break;
        } else {
          ++miit;
        }
      }

      if (first != miit) {
        ref_parent_ents_view.insert(first, second);
      }

    } while (miit != hi_miit);

    return ref_parent_ents_view;
  };

  auto ref_parent_ents_view = get_parent_ents_view(MBVERTEX, MBEDGE);

  auto tets = _tets.subset_by_type(MBTET);

  auto get_ele_parent_view = [&]() {
    auto &ref_ents = refined_finite_elements_ptr->get<Ent_mi_tag>();
    RefElement_multiIndex_parents_view ref_ele_parent_view;
    for (auto p = tets.pair_begin(); p != tets.pair_end(); ++p) {
      auto tmi = ref_ents.lower_bound(p->first);
      auto tme = ref_ents.upper_bound(p->second);
      ref_ele_parent_view.insert(tmi, tme);
    }
    return ref_ele_parent_view;
  };

  auto ref_ele_parent_view = get_ele_parent_view();

  std::vector<EntityHandle> parent_tets_refinded;
  std::vector<EntityHandle> vertices_of_created_tets;
  std::vector<BitRefEdges> parent_edges_bit_vec;
  std::vector<int> nb_new_tets_vec;
  std::vector<int> sub_type_vec;

  parent_tets_refinded.reserve(tets.size());
  vertices_of_created_tets.reserve(4 * tets.size());
  parent_edges_bit_vec.reserve(tets.size());
  nb_new_tets_vec.reserve(tets.size());
  sub_type_vec.reserve(tets.size());

  // make loop over all tets which going to be refined
  for (auto tit : tets) {

    // get tet connectivity
    const EntityHandle *conn;
    int num_nodes;
    CHKERR moab.get_connectivity(tit, conn, num_nodes, true);

    // get edges
    BitRefEdges parent_edges_bit(0);
    EntityHandle edge_new_nodes[] = {0, 0, 0, 0, 0, 0};
    int split_edges[] = {-1, -1, -1, -1, -1, -1};

    for (int ee = 0; ee != 6; ++ee) {
      EntityHandle edge;
      CHKERR moab.side_element(tit, 1, ee, edge);
      auto eit = ref_parent_ents_view.find(edge);
      if (eit != ref_parent_ents_view.end()) {
        edge_new_nodes[ee] = (*eit)->getEnt();
        split_edges[parent_edges_bit.count()] = ee;
        parent_edges_bit.set(ee, 1);
      }
    }

    if (parent_edges_bit.count()) {

      // swap nodes forward
      EntityHandle _conn_[4];
      std::copy(&conn[0], &conn[4], &_conn_[0]);

      // build connectivity for ref tets
      EntityHandle new_tets_conns[8 * 4];
      std::fill(&new_tets_conns[0], &new_tets_conns[8 * 4], no_handle);

      int sub_type = -1, nb_new_tets = 0;
      switch (parent_edges_bit.count()) {
      case 1:
        sub_type = 0;
        tet_type_1(_conn_, split_edges[0], edge_new_nodes[split_edges[0]],
                   new_tets_conns);
        nb_new_tets = 2;
        break;
      case 2:
        sub_type =
            tet_type_2(_conn_, split_edges, edge_new_nodes, new_tets_conns);
        if (sub_type & (4 | 8 | 16)) {
          nb_new_tets = 3;
          break;
        } else if (sub_type == 1) {
          nb_new_tets = 4;
          break;
        };
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "Imposible case");
        break;
      case 3:
        sub_type =
            tet_type_3(_conn_, split_edges, edge_new_nodes, new_tets_conns);
        if (sub_type <= 4) {
          nb_new_tets = 4;
          break;
        } else if (sub_type <= 7) {
          nb_new_tets = 5;
          break;
        }
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "Imposible case");
      case 4:
        sub_type =
            tet_type_4(_conn_, split_edges, edge_new_nodes, new_tets_conns);
        if (sub_type == 0) {
          nb_new_tets = 5;
          break;
        } else if (sub_type <= 7) {
          nb_new_tets = 6;
          break;
        }
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "Imposible case");
      case 5:
        sub_type = tet_type_5(moab, _conn_, edge_new_nodes, new_tets_conns);
        nb_new_tets = 7;
        break;
      case 6:
        sub_type = 0;
        tet_type_6(moab, _conn_, edge_new_nodes, new_tets_conns);
        nb_new_tets = 8;
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "Imposible case");
      }

      // find that tets
      auto &ref_ele_by_parent_and_ref_edges =
          ref_ele_parent_view
              .get<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>();
      auto it_by_ref_edges = ref_ele_by_parent_and_ref_edges.equal_range(
          boost::make_tuple(tit, parent_edges_bit.to_ulong()));
      // check if tet with this refinement scheme already exits
      std::vector<EntityHandle> ents_to_set_bit_vec;
      if (std::distance(it_by_ref_edges.first, it_by_ref_edges.second) ==
          static_cast<size_t>(nb_new_tets)) {
        ents_to_set_bit_vec.reserve(nb_new_tets);
        for (int tt = 0; it_by_ref_edges.first != it_by_ref_edges.second;
             it_by_ref_edges.first++, tt++) {
          auto tet = (*it_by_ref_edges.first)->getEnt();
          ents_to_set_bit_vec.emplace_back(tet);
          // set ref tets entities
          if (debug) {
            // add this tet if exist to this ref level
            auto ref_tet_it = refined_ents_ptr->find(tet);
            if (ref_tet_it == refined_ents_ptr->end()) {
              SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                      "Tetrahedron should be in database");
            }
          }
        }
        ents_to_set_bit.insert_list(ents_to_set_bit_vec.begin(),
                                    ents_to_set_bit_vec.end());

      } else {
        // if this element was not refined or was refined with different
        // patterns of split edges create new elements

        parent_tets_refinded.emplace_back(tit);
        for (int tt = 0; tt != nb_new_tets; ++tt) {
          for (auto nn : {0, 1, 2, 3})
            vertices_of_created_tets.emplace_back(new_tets_conns[4 * tt + nn]);
        }
        parent_edges_bit_vec.emplace_back(parent_edges_bit);
        nb_new_tets_vec.emplace_back(nb_new_tets);
        sub_type_vec.emplace_back(sub_type);
      }
    } else {

      if (debug) {
        RefEntity_multiIndex::iterator tit_miit;
        tit_miit = refined_ents_ptr->find(tit);
        if (tit_miit == refined_ents_ptr->end())
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "can not find this tet");
      }

      ents_to_set_bit.insert(tit);
    }
  }

  // Create tets
  EntityHandle start_e = 0;
  EntityHandle *conn_moab;
  const int nb_new_tets = vertices_of_created_tets.size() / 4;
  read_util->get_element_connect(nb_new_tets, 4, MBTET, 0, start_e, conn_moab);
  std::copy(vertices_of_created_tets.begin(), vertices_of_created_tets.end(),
            conn_moab);
  CHKERR read_util->update_adjacencies(start_e, nb_new_tets, 4, conn_moab);
  ents_to_set_bit.insert(start_e, start_e + nb_new_tets - 1);

  // Create adj entities
  for (auto d : {1, 2}) {
    Range ents;
    CHKERR moab.get_adjacencies(ents_to_set_bit, d, true, ents,
                                moab::Interface::UNION);
  }

  SetParent set_parent;

  // Set parrents and adjacencies
  for (int idx = 0; idx != parent_tets_refinded.size(); ++idx) {

    const EntityHandle tit = parent_tets_refinded[idx];
    const BitRefEdges &parent_edges_bit = parent_edges_bit_vec[idx];
    const int nb_new_tets = nb_new_tets_vec[idx];
    const int sub_type = sub_type_vec[idx];

    std::array<EntityHandle, 8> ref_tets;
    for (int tt = 0; tt != nb_new_tets; ++tt, ++start_e)
      ref_tets[tt] = start_e;

    if (nb_new_tets) {

      int ref_type[2];
      ref_type[0] = parent_edges_bit.count();
      ref_type[1] = sub_type;

      CHKERR moab.tag_clear_data(cOre.get_th_RefType(), &ref_tets[0],
                                 nb_new_tets, ref_type);
      CHKERR moab.tag_clear_data(cOre.get_th_RefBitEdge(), &ref_tets[0],
                                 nb_new_tets, &parent_edges_bit);
      CHKERR moab.tag_clear_data(cOre.get_th_RefParentHandle(), &ref_tets[0],
                                 nb_new_tets, &tit);

      // hash map of nodes (RefEntity) by edges (EntityHandle)
      std::map<EntityHandle /*edge*/, EntityHandle /*node*/>
          map_ref_nodes_by_edges;

      Range tit_edges;
      CHKERR moab.get_adjacencies(&tit, 1, 1, false, tit_edges);
      for (auto edge : tit_edges) {
        auto eit = ref_parent_ents_view.find(edge);
        if (eit != ref_parent_ents_view.end()) {
          map_ref_nodes_by_edges[(*eit)->getParentEnt()] = eit->get()->getEnt();
        }
      }

      // find parents for new edges and faces
      // get tet edges and faces
      Range tit_faces;
      CHKERR moab.get_adjacencies(&tit, 1, 2, false, tit_faces);
      Range edges_nodes[6], faces_nodes[4];
      // for edges - add ref nodes
      // edges_nodes[ee] - contains all nodes on edge ee including mid nodes if
      // exist
      Range::iterator eit = tit_edges.begin();
      for (int ee = 0; eit != tit_edges.end(); eit++, ee++) {
        CHKERR moab.get_connectivity(&*eit, 1, edges_nodes[ee], true);
        auto map_miit = map_ref_nodes_by_edges.find(*eit);
        if (map_miit != map_ref_nodes_by_edges.end()) {
          edges_nodes[ee].insert(map_miit->second);
        }
      }
      // for faces - add ref nodes
      // faces_nodes[ff] - contains all nodes on face ff including mid nodes if
      // exist
      Range::iterator fit = tit_faces.begin();
      for (int ff = 0; fit != tit_faces.end(); fit++, ff++) {
        CHKERR moab.get_connectivity(&*fit, 1, faces_nodes[ff], true);
        // Get edges on face and loop over those edges to add mid-nodes to range
        Range fit_edges;
        CHKERR moab.get_adjacencies(&*fit, 1, 1, false, fit_edges);
        for (Range::iterator eit2 = fit_edges.begin(); eit2 != fit_edges.end();
             eit2++) {
          auto map_miit = map_ref_nodes_by_edges.find(*eit2);
          if (map_miit != map_ref_nodes_by_edges.end()) {
            faces_nodes[ff].insert(map_miit->second);
          }
        }
      }
      // add ref nodes to tet
      // tet_nodes contains all nodes on tet including mid edge nodes
      Range tet_nodes;
      CHKERR moab.get_connectivity(&tit, 1, tet_nodes, true);
      for (std::map<EntityHandle, EntityHandle>::iterator map_miit =
               map_ref_nodes_by_edges.begin();
           map_miit != map_ref_nodes_by_edges.end(); map_miit++) {
        tet_nodes.insert(map_miit->second);
      }
      Range ref_edges;
      // Get all all edges of refined tets
      CHKERR moab.get_adjacencies(ref_tets.data(), nb_new_tets, 1, false,
                                  ref_edges, moab::Interface::UNION);
      // Check for all ref edge and set parents
      for (Range::iterator reit = ref_edges.begin(); reit != ref_edges.end();
           reit++) {
        Range ref_edges_nodes;
        CHKERR moab.get_connectivity(&*reit, 1, ref_edges_nodes, true);

        // Check if ref edge is an coarse edge (loop over coarse tet edges)
        int ee = 0;
        for (; ee < 6; ee++) {
          // Two nodes are common (node[0],node[1],ref_node (if exist))
          // this tests if given edge is contained by edge of refined
          // tetrahedral
          if (intersect(edges_nodes[ee], ref_edges_nodes).size() == 2) {
            EntityHandle edge = tit_edges[ee];
            CHKERR set_parent(*reit, edge, refined_ents_ptr, cOre);
            break;
          }
        }
        if (ee < 6)
          continue; // this refined edge is contained by edge of tetrahedral
        // check if ref edge is in coarse face
        int ff = 0;
        for (; ff < 4; ff++) {
          // two nodes are common (node[0],node[1],ref_node (if exist))
          // this tests if given face is contained by face of  tetrahedral
          if (intersect(faces_nodes[ff], ref_edges_nodes).size() == 2) {
            EntityHandle face = tit_faces[ff];
            CHKERR set_parent(*reit, face, refined_ents_ptr, cOre);
            break;
          }
        }
        if (ff < 4)
          continue; // this refined edge is contained by face of tetrahedral

        // check if ref edge is in coarse tetrahedral (i.e. that is internal
        // edge of refined tetrahedral)
        if (intersect(tet_nodes, ref_edges_nodes).size() == 2) {
          CHKERR set_parent(*reit, tit, refined_ents_ptr, cOre);
          continue;
        }

        // Refined edge is not child of any edge, face or tetrahedral, this is
        // imposible edge
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "impossible refined edge");
      }

      Range ref_faces;
      CHKERR moab.get_adjacencies(ref_tets.data(), nb_new_tets, 2, false,
                                  ref_faces, moab::Interface::UNION);
      Tag th_interface_side;
      const int def_side[] = {0};
      CHKERR moab.tag_get_handle("INTERFACE_SIDE", 1, MB_TYPE_INTEGER,
                                 th_interface_side,
                                 MB_TAG_CREAT | MB_TAG_SPARSE, def_side);
      // Check for all ref faces
      for (auto rfit : ref_faces) {
        Range ref_faces_nodes;
        CHKERR moab.get_connectivity(&rfit, 1, ref_faces_nodes, true);
        // Check if ref face is in coarse face
        int ff = 0;
        for (; ff < 4; ff++) {
          // Check if refined triangle is contained by face of tetrahedral
          if (intersect(faces_nodes[ff], ref_faces_nodes).size() == 3) {
            EntityHandle face = tit_faces[ff];
            CHKERR set_parent(rfit, face, refined_ents_ptr, cOre);
            int side = 0;
            // Set face side if it is on interface
            CHKERR moab.tag_get_data(th_interface_side, &face, 1, &side);
            CHKERR moab.tag_set_data(th_interface_side, &rfit, 1, &side);
            break;
          }
        }
        if (ff < 4)
          continue; // this face is contained by one of tetrahedrons
        // check if ref face is in coarse tetrahedral
        // this is ref face which is contained by tetrahedral volume
        if (intersect(tet_nodes, ref_faces_nodes).size() == 3) {
          CHKERR set_parent(rfit, tit, refined_ents_ptr, cOre);
          continue;
        }
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "impossible refined face");
      }
    }
  }

  CHKERR set_parent(refined_ents_ptr);
  CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevel(ents_to_set_bit,
                                                               bit, true, verb);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode MeshRefinement::refinePrisms(const EntityHandle meshset,
                                            const BitRefLevel &bit, int verb) {

  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  auto refined_ents_ptr = m_field.get_ref_ents();
  auto refined_finite_elements_ptr = m_field.get_ref_finite_elements();

  // FIXME: refinement is based on entity handlers, should work on global ids of
  // nodes, this will allow parallelise algorithm in the future

  MoFEMFunctionBegin;

  RefElement_multiIndex_parents_view ref_ele_parent_view;
  ref_ele_parent_view.insert(
      refined_finite_elements_ptr->get<Ent_mi_tag>().lower_bound(
          get_id_for_min_type<MBPRISM>()),
      refined_finite_elements_ptr->get<Ent_mi_tag>().upper_bound(
          get_id_for_max_type<MBPRISM>()));
  auto &ref_ele_by_parent_and_ref_edges =
      ref_ele_parent_view
          .get<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>();
  // find all vertices which parent is edge
  auto &ref_ents_by_comp =
      refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>();
  RefEntity_multiIndex_view_by_hashed_parent_entity ref_parent_ents_view;
  ref_parent_ents_view.insert(
      ref_ents_by_comp.lower_bound(boost::make_tuple(MBVERTEX, MBEDGE)),
      ref_ents_by_comp.upper_bound(boost::make_tuple(MBVERTEX, MBEDGE)));
  Range prisms;
  CHKERR moab.get_entities_by_type(meshset, MBPRISM, prisms, false);
  Range::iterator pit = prisms.begin();
  for (; pit != prisms.end(); pit++) {
    auto miit_prism = refined_ents_ptr->get<Ent_mi_tag>().find(*pit);
    if (miit_prism == refined_ents_ptr->end()) {
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "this prism is not in ref database");
    }
    if (verb >= NOISY) {
      MOFEM_TAG_AND_LOG("WORLD", Sev::noisy, "MeshRefinment")
          << "ref prism " << **miit_prism;
    }
    // prism connectivity
    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(*pit, conn, num_nodes, true);
    assert(num_nodes == 6);
    // edges connectivity
    EntityHandle edges[6];
    for (int ee = 0; ee < 3; ee++) {
      CHKERR moab.side_element(*pit, 1, ee, edges[ee]);
    }
    for (int ee = 6; ee < 9; ee++) {
      CHKERR moab.side_element(*pit, 1, ee, edges[ee - 3]);
    }
    // detect split edges
    BitRefEdges split_edges(0);
    EntityHandle edge_nodes[6];
    std::fill(&edge_nodes[0], &edge_nodes[6], no_handle);
    for (int ee = 0; ee < 6; ee++) {
      RefEntity_multiIndex_view_by_hashed_parent_entity::iterator miit_view =
          ref_parent_ents_view.find(edges[ee]);
      if (miit_view != ref_parent_ents_view.end()) {
        if (((*miit_view)->getBitRefLevel() & bit).any()) {
          edge_nodes[ee] = (*miit_view)->getEnt();
          split_edges.set(ee);
        }
      }
    }
    if (split_edges.count() == 0) {
      *(const_cast<RefEntity *>(miit_prism->get())->getBitRefLevelPtr()) |= bit;
      if (verb >= VERY_NOISY)
        MOFEM_TAG_AND_LOG("WORLD", Sev::noisy, "MeshRefinment")
            << "no refinement";
      continue;
    }
    // check consistency
    if (verb >= NOISY) {
      std::ostringstream ss;
      ss << "prism split edges " << split_edges << " count "
         << split_edges.count() << std::endl;
      PetscPrintf(m_field.get_comm(), ss.str().c_str());
    }
    // prism ref
    EntityHandle new_prism_conn[4 * 6];
    std::fill(&new_prism_conn[0], &new_prism_conn[4 * 6], no_handle);
    int nb_new_prisms = 0;
    switch (split_edges.count()) {
    case 0:
      break;
    case 2:
      CHKERR prism_type_1(conn, split_edges, edge_nodes, new_prism_conn);
      nb_new_prisms = 2;
      break;
    case 4:
      CHKERR prism_type_2(conn, split_edges, edge_nodes, new_prism_conn);
      nb_new_prisms = 3;
      break;
    case 6:
      CHKERR prism_type_3(conn, split_edges, edge_nodes, new_prism_conn);
      nb_new_prisms = 4;
      break;
    default:
      std::ostringstream ss;
      ss << split_edges << " : [ " << conn[0] << " " << conn[1] << " "
         << conn[2] << " " << conn[3] << " " << conn[4] << " " << conn[5]
         << " ]";
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY, ss.str().c_str());
    }
    // find that prism
    std::bitset<4> ref_prism_bit(0);
    auto it_by_ref_edges = ref_ele_by_parent_and_ref_edges.lower_bound(
        boost::make_tuple(*pit, split_edges.to_ulong()));
    auto hi_it_by_ref_edges = ref_ele_by_parent_and_ref_edges.upper_bound(
        boost::make_tuple(*pit, split_edges.to_ulong()));
    auto it_by_ref_edges2 = it_by_ref_edges;
    for (int pp = 0; it_by_ref_edges2 != hi_it_by_ref_edges;
         it_by_ref_edges2++, pp++) {
      // Add this tet to this ref
      *(const_cast<RefElement *>(it_by_ref_edges2->get())
            ->getBitRefLevelPtr()) |= bit;
      ref_prism_bit.set(pp, 1);
      if (verb > 2) {
        std::ostringstream ss;
        ss << "is refined " << *it_by_ref_edges2 << std::endl;
        PetscPrintf(m_field.get_comm(), ss.str().c_str());
      }
    }
    if (it_by_ref_edges != hi_it_by_ref_edges) {
      if (ref_prism_bit.count() != (unsigned int)nb_new_prisms)
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
    } else {
      EntityHandle ref_prisms[4];
      // create prism
      for (int pp = 0; pp < nb_new_prisms; pp++) {
        if (verb > 3) {
          std::ostringstream ss;
          ss << "ref prism " << ref_prism_bit << std::endl;
          PetscPrintf(m_field.get_comm(), ss.str().c_str());
        }
        if (!ref_prism_bit.test(pp)) {
          CHKERR moab.create_element(MBPRISM, &new_prism_conn[6 * pp], 6,
                                     ref_prisms[pp]);
          CHKERR moab.tag_set_data(cOre.get_th_RefParentHandle(),
                                   &ref_prisms[pp], 1, &*pit);
          CHKERR moab.tag_set_data(cOre.get_th_RefBitLevel(), &ref_prisms[pp],
                                   1, &bit);
          CHKERR moab.tag_set_data(cOre.get_th_RefBitEdge(), &ref_prisms[pp], 1,
                                   &split_edges);
          auto p_ent =
              const_cast<RefEntity_multiIndex *>(refined_ents_ptr)
                  ->insert(boost::shared_ptr<RefEntity>(new RefEntity(
                      m_field.get_basic_entity_data_ptr(), ref_prisms[pp])));
          std::pair<RefElement_multiIndex::iterator, bool> p_fe;
          try {
            p_fe =
                const_cast<RefElement_multiIndex *>(refined_finite_elements_ptr)
                    ->insert(boost::shared_ptr<RefElement>(
                        new RefElement_PRISM(*p_ent.first)));
          } catch (MoFEMException const &e) {
            SETERRQ(PETSC_COMM_SELF, e.errorCode, e.errorMessage);
          }
          ref_prism_bit.set(pp);
          CHKERR cOre.addPrismToDatabase(ref_prisms[pp]);
          if (verb > 2) {
            std::ostringstream ss;
            ss << "add prism: " << **p_fe.first << std::endl;
            if (verb > 7) {
              for (int nn = 0; nn < 6; nn++) {
                ss << new_prism_conn[nn] << " ";
              }
              ss << std::endl;
            }
            PetscPrintf(m_field.get_comm(), ss.str().c_str());
          }
        }
      }
    }
  }
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode MeshRefinement::refineMeshset(const EntityHandle meshset,
                                             const BitRefLevel &bit,
                                             const bool recursive, int verb) {
  Interface &m_field = cOre;
  auto refined_ents_ptr = m_field.get_ref_ents();
  MoFEMFunctionBegin;
  auto miit = refined_ents_ptr->find(meshset);
  if (miit == refined_ents_ptr->end()) {
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "this meshset is not in ref database");
  }
  CHKERR m_field.getInterface<BitRefManager>()->updateMeshsetByEntitiesChildren(
      meshset, bit, meshset, MBEDGE, recursive, verb);
  CHKERR m_field.getInterface<BitRefManager>()->updateMeshsetByEntitiesChildren(
      meshset, bit, meshset, MBTRI, recursive, verb);
  CHKERR m_field.getInterface<BitRefManager>()->updateMeshsetByEntitiesChildren(
      meshset, bit, meshset, MBTET, recursive, verb);
  *(const_cast<RefEntity *>(miit->get())->getBitRefLevelPtr()) |= bit;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshRefinement::refineTris(const EntityHandle meshset,
                                          const BitRefLevel &bit, int verb,
                                          const bool debug) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  Range tris;
  CHKERR moab.get_entities_by_type(meshset, MBTRI, tris, false);
  CHKERR refineTris(tris, bit, verb, debug);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshRefinement::refineTris(const Range &_tris,
                                          const BitRefLevel &bit, int verb,
                                          const bool debug) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  auto refined_ents_ptr = m_field.get_ref_ents();
  auto refined_finite_elements_ptr = m_field.get_ref_finite_elements();
  ReadUtilIface *read_util;

  MoFEMFunctionBegin;

  auto get_parent_ents_view = [&](const auto type_child,
                                  const auto type_parent) {
    auto &ref_ents =
        refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>();
    auto miit =
        ref_ents.lower_bound(boost::make_tuple(type_child, type_parent));
    auto hi_miit =
        ref_ents.upper_bound(boost::make_tuple(type_child, type_parent));
    RefEntity_multiIndex_view_by_ordered_parent_entity ref_parent_ents_view;

    decltype(miit) first;
    decltype(miit) second;

    do {

      while (miit != hi_miit) {
        first = miit;
        if (((*miit)->getBitRefLevel() & bit).any()) {
          ++miit;
          second = miit;
          break;
        } else {
          ++miit;
        }
      }

      while (miit != hi_miit) {
        if (((*miit)->getBitRefLevel() & bit).any()) {
          second = miit;
          break;
        } else {
          ++miit;
        }
      }

      if (first != miit) {
        ref_parent_ents_view.insert(first, second);
      }

    } while (miit != hi_miit);

    return ref_parent_ents_view;
  };

  auto ref_parent_ents_view = get_parent_ents_view(MBVERTEX, MBEDGE);

  auto tris = _tris.subset_by_type(MBTRI);

  auto get_ele_parent_view = [&]() {
    auto &ref_ents = refined_finite_elements_ptr->get<Ent_mi_tag>();
    RefElement_multiIndex_parents_view ref_ele_parent_view;
    for (auto p = tris.pair_begin(); p != tris.pair_end(); ++p) {
      auto tmi = ref_ents.lower_bound(p->first);
      auto tme = ref_ents.upper_bound(p->second);
      ref_ele_parent_view.insert(tmi, tme);
    }
    return ref_ele_parent_view;
  };

  auto ref_ele_parent_view = get_ele_parent_view();

  Range ents_to_set_bit;

  std::vector<EntityHandle> parent_tris_refined; // entities refined
  std::vector<EntityHandle> vertices_of_created_tris;
  std::vector<BitRefEdges> parent_edges_bit_vec;
  std::vector<int> nb_new_tris_vec;

  parent_tris_refined.reserve(tris.size());
  vertices_of_created_tris.reserve(3 * tris.size());
  parent_edges_bit_vec.reserve(tris.size());
  nb_new_tris_vec.reserve(tris.size());

  // make loop over all tets which going to be refined
  for (auto tit : tris) {

    // get tet connectivity
    const EntityHandle *conn;
    int num_nodes;
    CHKERR moab.get_connectivity(tit, conn, num_nodes, true);

    // get edges
    BitRefEdges parent_edges_bit(0);
    EntityHandle edge_new_nodes[] = {0, 0, 0};
    int split_edges[] = {-1, -1, -1};

    for (int ee = 0; ee != 3; ++ee) {
      EntityHandle edge;
      CHKERR moab.side_element(tit, 1, ee, edge);
      auto eit = ref_parent_ents_view.find(edge);
      if (eit != ref_parent_ents_view.end()) {
        edge_new_nodes[ee] = (*eit)->getEnt();
        split_edges[parent_edges_bit.count()] = ee;
        parent_edges_bit.set(ee, 1);
      }
    }

    if (parent_edges_bit.count()) {

      // swap nodes forward
      EntityHandle _conn_[3];
      std::copy(&conn[0], &conn[3], &_conn_[0]);

      // build connectivity for ref tets
      EntityHandle new_tets_conns[4 * 3];
      std::fill(&new_tets_conns[0], &new_tets_conns[4 * 3], no_handle);

      int nb_new_tris = 0;
      switch (parent_edges_bit.count()) {
      case 1:
        CHKERR tri_type_1(conn, split_edges[0], edge_new_nodes[split_edges[0]],
                          new_tets_conns);
        nb_new_tris += 2;
        break;
      case 2:
        CHKERR tri_type_2(conn, split_edges, edge_new_nodes, new_tets_conns);
        nb_new_tris += 3;
        break;
      case 3:
        CHKERR tri_type_3(conn, edge_new_nodes, new_tets_conns);
        nb_new_tris += 4;
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "Imposible case");
      }

      // find that tets
      auto &ref_ele_by_parent_and_ref_edges =
          ref_ele_parent_view
              .get<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>();
      auto it_by_ref_edges = ref_ele_by_parent_and_ref_edges.equal_range(
          boost::make_tuple(tit, parent_edges_bit.to_ulong()));
      // check if tet with this refinement scheme already exits
      std::vector<EntityHandle> ents_to_set_bit_vec;
      if (std::distance(it_by_ref_edges.first, it_by_ref_edges.second) ==
          static_cast<size_t>(nb_new_tris)) {
        ents_to_set_bit_vec.reserve(nb_new_tris);
        for (int tt = 0; it_by_ref_edges.first != it_by_ref_edges.second;
             it_by_ref_edges.first++, tt++) {
          auto tet = (*it_by_ref_edges.first)->getEnt();
          ents_to_set_bit_vec.emplace_back(tet);
          // set ref tets entities
          if (debug) {
            // add this tet if exist to this ref level
            auto ref_tet_it = refined_ents_ptr->find(tet);
            if (ref_tet_it == refined_ents_ptr->end()) {
              SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                      "Tetrahedron should be in database");
            }
          }
        }
        ents_to_set_bit.insert_list(ents_to_set_bit_vec.begin(),
                                    ents_to_set_bit_vec.end());

      } else {
        // if this element was not refined or was refined with different
        // patterns of split edges create new elements

        parent_tris_refined.emplace_back(tit);
        for (int tt = 0; tt != nb_new_tris; ++tt) {
          for (auto nn : {0, 1, 2})
            vertices_of_created_tris.emplace_back(new_tets_conns[3 * tt + nn]);
        }
        parent_edges_bit_vec.emplace_back(parent_edges_bit);
        nb_new_tris_vec.emplace_back(nb_new_tris);
      }
    } else {

      if (debug) {
        RefEntity_multiIndex::iterator tit_miit;
        tit_miit = refined_ents_ptr->find(tit);
        if (tit_miit == refined_ents_ptr->end())
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "can not find this tet");
      }

      ents_to_set_bit.insert(tit);
    }
  }

  // Create tets
  EntityHandle start_e = 0;
  EntityHandle *conn_moab;
  const int nb_new_tris = vertices_of_created_tris.size() / 3;
  CHKERR m_field.get_moab().query_interface(read_util);
  CHKERR read_util->get_element_connect(nb_new_tris, 3, MBTRI, 0, start_e,
                                        conn_moab);
  std::copy(vertices_of_created_tris.begin(), vertices_of_created_tris.end(),
            conn_moab);
  CHKERR read_util->update_adjacencies(start_e, nb_new_tris, 3, conn_moab);
  ents_to_set_bit.insert(start_e, start_e + nb_new_tris - 1);

  if (debug) {
    auto meshset_ptr = get_temp_meshset_ptr(moab);
    CHKERR moab.add_entities(*meshset_ptr, ents_to_set_bit);
    CHKERR moab.write_file("out_refinded_debug.vtk", "VTK", "",
                           meshset_ptr->get_ptr(), 1);
  }

  // Create adj entities
  for (auto d : {1}) {
    Range ents;
    CHKERR moab.get_adjacencies(ents_to_set_bit, d, true, ents,
                                moab::Interface::UNION);
  }

  SetParent set_parent;

  // Set parrents and adjacencies. Iterate over all parent triangles.
  for (int idx = 0; idx != parent_tris_refined.size(); ++idx) {

    const auto tit = parent_tris_refined[idx]; // parent triangle
    const auto &parent_edges_bit =
        parent_edges_bit_vec[idx]; // marks what edges are refined
    const auto nb_new_tris =
        nb_new_tris_vec[idx]; // number of triangles created by refinement

    // Iterate over refined triangles
    std::array<EntityHandle, 4> ref_tris;
    for (int tt = 0; tt != nb_new_tris; ++tt, ++start_e)
      ref_tris[tt] = start_e;

    if (nb_new_tris) {

      int ref_type[2];
      ref_type[0] = parent_edges_bit.count();
      ref_type[1] = 0;

      // Set tags for reference dtriangles
      CHKERR moab.tag_clear_data(cOre.get_th_RefType(), &ref_tris[0],
                                 nb_new_tris, ref_type);
      CHKERR moab.tag_clear_data(cOre.get_th_RefBitEdge(), &ref_tris[0],
                                 nb_new_tris, &parent_edges_bit);
      CHKERR moab.tag_clear_data(cOre.get_th_RefParentHandle(), &ref_tris[0],
                                 nb_new_tris, &tit);

      // Hash map of nodes (RefEntity) by edges (EntityHandle)
      std::map<EntityHandle /*edge*/, EntityHandle /*node*/>
          map_ref_nodes_by_edges;

      // find parents for new edges and faces
      // get tet edges and faces
      Range tit_edges; // parent edges
      CHKERR moab.get_adjacencies(&tit, 1, 1, false, tit_edges);
      for (auto edge : tit_edges) {
        auto eit = ref_parent_ents_view.find(edge);
        if (eit != ref_parent_ents_view.end()) {
          map_ref_nodes_by_edges[(*eit)->getParentEnt()] = eit->get()->getEnt();
        }
      }

      Range edges_nodes[3];
      // for edges - add ref nodes
      // edges_nodes[ee] - contains all nodes on edge ee including mid nodes if
      // exist
      Range::iterator eit = tit_edges.begin();
      for (int ee = 0; eit != tit_edges.end(); eit++, ee++) {
        CHKERR moab.get_connectivity(&*eit, 1, edges_nodes[ee], true);
        auto map_miit = map_ref_nodes_by_edges.find(*eit);
        if (map_miit != map_ref_nodes_by_edges.end()) {
          edges_nodes[ee].insert(map_miit->second);
        }
      }

      // add ref nodes to tri
      // tet_nodes contains all nodes on tet including mid edge nodes
      Range tri_nodes;
      CHKERR moab.get_connectivity(&tit, 1, tri_nodes, true);
      for (auto map_miit = map_ref_nodes_by_edges.begin();
           map_miit != map_ref_nodes_by_edges.end(); map_miit++) {
        tri_nodes.insert(map_miit->second);
      }

      Range ref_edges; // edges on all new tets
      // Get all all edges of refined tets
      CHKERR moab.get_adjacencies(ref_tris.data(), nb_new_tris, 1, false,
                                  ref_edges, moab::Interface::UNION);

      // Check for all ref edge and set parents
      for (auto reit = ref_edges.begin(); reit != ref_edges.end(); ++reit) {
        Range ref_edges_nodes;
        CHKERR moab.get_connectivity(&*reit, 1, ref_edges_nodes, true);
        // Check if ref edge is an coarse edge (loop over coarse tet edges)
        int ee = 0;
        for (; ee < 3; ee++) {
          // Two nodes are common (node[0],node[1],ref_node (if exist))
          // this tests if given edge is contained by edge of refined
          // triangle
          if (intersect(edges_nodes[ee], ref_edges_nodes).size() == 2) {
            EntityHandle edge = tit_edges[ee];
            CHKERR set_parent(*reit, edge, refined_ents_ptr, cOre);
            break;
          }
        }
        if (ee < 3)
          continue; // this refined edge is contained by edge of triangel

        // check if ref edge is in coarse tetrahedral (i.e. that is internal
        // edge of refined tetrahedral)
        if (intersect(tri_nodes, ref_edges_nodes).size() == 2) {
          CHKERR set_parent(*reit, tit, refined_ents_ptr, cOre);
          continue;
        }

        // Refined edge is not child of any edge, face or tetrahedral, this is
        // imposible edge
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "impossible refined edge");
      }
    }
  }

  CHKERR set_parent(refined_ents_ptr);
  CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevel(ents_to_set_bit,
                                                               bit, true, verb);

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
