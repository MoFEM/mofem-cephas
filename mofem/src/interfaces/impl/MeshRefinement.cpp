/** \file MeshRefinementCore.cpp
 * \brief FIXME this is not so good implementation
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

MoFEMErrorCode MeshRefinement::query_interface(const MOFEMuuid &uuid,
                                               UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMMeshRefine) {
    *iface = const_cast<MeshRefinement *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

MeshRefinement::MeshRefinement(const Core &core)
    : cOre(const_cast<Core &>(core)) {}

MoFEMErrorCode MeshRefinement::add_vertices_in_the_middel_of_edges(
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
  CHKERR add_vertices_in_the_middel_of_edges(edges, bit, verb, start_v);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode MeshRefinement::add_vertices_in_the_middel_of_edges(
    const Range &ents, const BitRefLevel &bit, int verb, EntityHandle start_v) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  const RefEntity_multiIndex *refined_ents_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_ref_ents(&refined_ents_ptr);
  Range edges = ents.subset_by_type(MBEDGE);
  auto miit =
      refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>()
          .lower_bound(boost::make_tuple(MBVERTEX, MBEDGE));
  auto hi_miit =
      refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>()
          .upper_bound(boost::make_tuple(MBVERTEX, MBEDGE));
  RefEntity_multiIndex_view_by_ordered_parent_entity ref_parent_ents_view;
  ref_parent_ents_view.insert(miit, hi_miit);
  if (verb >= VERBOSE) {
    std::ostringstream ss;
    ss << "ref level " << bit << " nb. edges to refine " << edges.size()
       << std::endl;
    PetscPrintf(m_field.get_comm(), ss.str().c_str());
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
      auto hi_miit_view = ref_parent_ents_view.upper_bound(p_eit->second);
      for (; miit_view != hi_miit_view; ++miit_view)
        edge_having_parent_vertex.insert(miit_view->get()->getParentEnt());
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
    int num_nodes = vert_coords[0].size();
    vector<double *> arrays_coord;
    EntityHandle startv;
    {
      ReadUtilIface *iface;
      CHKERR moab.query_interface(iface);
      CHKERR iface->get_node_coords(3, num_nodes, start_v, startv,
                                    arrays_coord);
    }
    Range verts(startv, startv + num_nodes - 1);
    for (int dd = 0; dd != 3; ++dd) {
      std::copy(vert_coords[dd].begin(), vert_coords[dd].end(),
                arrays_coord[dd]);
    }
    CHKERR moab.tag_set_data(cOre.get_th_RefParentHandle(), verts,
                             &*parent_edge.begin());
    CHKERR m_field.getInterface<BitRefManager>()->setEntitiesBitRefLevel(
        verts, bit, verb);
}
MoFEMFunctionReturn(0);
} // namespace MoFEM

MoFEMErrorCode MeshRefinement::refine_TET(const EntityHandle meshset,
                                          const BitRefLevel &bit,
                                          const bool respect_interface,
                                          int verb, Range *ref_edges_ptr) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  Range tets;
  CHKERR moab.get_entities_by_type(meshset, MBTET, tets, false);
  CHKERR refine_TET(tets, bit, respect_interface, verb, ref_edges_ptr);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshRefinement::refine_TET(const Range &_tets,
                                          const BitRefLevel &bit,
                                          const bool respect_interface,
                                          int verb, Range *ref_edges_ptr) {

  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  const RefEntity_multiIndex *refined_ents_ptr;
  const RefElement_multiIndex *refined_finite_elements_ptr;
  MoFEMFunctionBegin;

  // Check if refinement is correct
  struct Check {
    map<EntityHandle, EntityHandle> entParentMap;
    MoFEMErrorCode operator()(Range *ref_edges,
                              const RefEntity_multiIndex *ref_ents_ptr,
                              MoFEM::Core &core) {
      MoFEMFunctionBegin;
      if (!ref_edges)
        MoFEMFunctionReturnHot(0);
      for (Range::iterator eit = ref_edges->begin(); eit != ref_edges->end();
           ++eit) {
        RefEntity_multiIndex::index<
            Composite_ParentEnt_And_EntType_mi_tag>::type::iterator vit =
            ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>().find(
                boost::make_tuple(MBVERTEX, *eit));
        if (vit ==
            ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>().end()) {
          RefEntity_multiIndex::iterator e_eit = ref_ents_ptr->find(*eit);
          if (e_eit == ref_ents_ptr->end()) {
            SETERRQ1(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
                     "Edge not found %ld", *eit);
          }
          cerr << "Parent edge" << endl << **e_eit << endl;
          if (entParentMap.find(*eit) != entParentMap.end()) {
            RefEntity_multiIndex::iterator v_eit =
                ref_ents_ptr->find(entParentMap[*eit]);
            if (v_eit == ref_ents_ptr->end()) {
              SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
                      "Vertex not found");
            }
            cerr << "Vertex " << **v_eit << endl;
          }
          RefEntity_multiIndex::index<Ent_Ent_mi_tag>::type::iterator ee_it,
              ee_hi_it;
          ee_it = ref_ents_ptr->get<Ent_Ent_mi_tag>().lower_bound(*eit);
          ee_hi_it = ref_ents_ptr->get<Ent_Ent_mi_tag>().upper_bound(*eit);
          for (; ee_it != ee_hi_it; ++ee_it) {
            cerr << "Ent having edge parent by parent " << **ee_it << endl;
          }
          RefEntity_multiIndex tmp_index;
          tmp_index.insert(ref_ents_ptr->begin(), ref_ents_ptr->end());
          RefEntity_multiIndex::index<
              Composite_ParentEnt_And_EntType_mi_tag>::type::iterator vvit =
              tmp_index.get<Composite_ParentEnt_And_EntType_mi_tag>().find(
                  boost::make_tuple(MBVERTEX, *eit));
          if (vvit !=
              tmp_index.get<Composite_ParentEnt_And_EntType_mi_tag>().end()) {
            cerr << "Tmp idx Vertex " << **vvit << endl;
          }
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "No vertex on trim edges, that make no sense");
        } else {
          entParentMap[vit->get()->getParentEnt()] = vit->get()->getRefEnt();
        }
      }
      MoFEMFunctionReturn(0);
    }
  };

  struct SetParent {
    map<EntityHandle, EntityHandle> parentsToChange;
    MoFEMErrorCode operator()(const EntityHandle ent, const EntityHandle parent,
                              const RefEntity_multiIndex *ref_ents_ptr,
                              MoFEM::Core &cOre) {
      MoFEM::Interface &m_field = cOre;
      MoFEMFunctionBegin;
      RefEntity_multiIndex::iterator it = ref_ents_ptr->find(ent);
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
    MoFEMErrorCode operator()(const RefEntity_multiIndex *ref_ents_ptr) {
      MoFEMFunctionBegin;
      for (map<EntityHandle, EntityHandle>::iterator mit =
               parentsToChange.begin();
           mit != parentsToChange.end(); ++mit) {
        RefEntity_multiIndex::iterator it = ref_ents_ptr->find(mit->first);
        bool success = const_cast<RefEntity_multiIndex *>(ref_ents_ptr)
                           ->modify(it, RefEntity_change_parent(mit->second));
        if (!success) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "impossible to set parent");
        }
      }
      MoFEMFunctionReturn(0);
    }
  };
  SetParent set_parent;

  Range ents_to_set_bit;

  CHKERR m_field.get_ref_ents(&refined_ents_ptr);
  CHKERR m_field.get_ref_finite_elements(&refined_finite_elements_ptr);

  Check check;
  CHKERR check(ref_edges_ptr, refined_ents_ptr, cOre);

  // FIXME: refinement is based on entity handlers, should work on global ids of
  // nodes, this will allow parallelise algorithm in the future

  // Find all vertices which parent is edge
  typedef const RefEntity_multiIndex::index<
      Composite_EntType_and_ParentEntType_mi_tag>::type RefEntsByComposite;
  RefEntsByComposite &ref_ents =
      refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>();
  RefEntity_multiIndex_view_by_hashed_parent_entity ref_parent_ents_view;
  ref_parent_ents_view.insert(
      ref_ents.lower_bound(boost::make_tuple(MBVERTEX, MBEDGE)),
      ref_ents.upper_bound(boost::make_tuple(MBVERTEX, MBEDGE)));
  typedef const RefElement_multiIndex::index<Ent_mi_tag>::type RefElementByEnt;
  RefElementByEnt &ref_finite_element =
      refined_finite_elements_ptr->get<Ent_mi_tag>();
  typedef const RefElement_multiIndex_parents_view::index<
      Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>::type
      RefEntByParentAndRefEdges;
  RefElement_multiIndex_parents_view ref_ele_parent_view;
  ref_ele_parent_view.insert(
      refined_finite_elements_ptr->get<EntType_mi_tag>().lower_bound(MBTET),
      refined_finite_elements_ptr->get<EntType_mi_tag>().upper_bound(MBTET));
  RefEntByParentAndRefEdges &ref_ele_by_parent_and_ref_edges =
      ref_ele_parent_view
          .get<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>();

  if (respect_interface) {
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "not implemented, set last parameter in refine_TET to false");
  }

  // make loop over all tets which going to be refined
  Range tets = _tets.subset_by_type(MBTET);
  Range::iterator tit = tets.begin();
  for (; tit != tets.end(); ++tit) {

    if (ref_finite_element.find(*tit) == ref_finite_element.end()) {
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "this tet is not in refinedFiniteElements");
    }

    // get tet connectivity
    const EntityHandle *conn;
    int num_nodes;
    CHKERR moab.get_connectivity(*tit, conn, num_nodes, true);

    // get edges
    BitRefEdges parent_edges_bit(0);
    EntityHandle edge_new_nodes[] = {0, 0, 0, 0, 0, 0};
    int split_edges[] = {-1, -1, -1, -1, -1, -1};
    // hash map of nodes (RefEntity) by edges (EntityHandle)
    std::map<EntityHandle /*edge*/, EntityHandle /*node*/>
        map_ref_nodes_by_edges;
    for (int ee = 0; ee != 6; ++ee) {
      EntityHandle edge;
      CHKERR moab.side_element(*tit, 1, ee, edge);
      RefEntity_multiIndex_view_by_hashed_parent_entity::iterator eit =
          ref_parent_ents_view.find(edge);
      if (eit != ref_parent_ents_view.end()) {
        if (((*eit)->getBitRefLevel() & bit).any()) {
          edge_new_nodes[ee] = (*eit)->getRefEnt();
          map_ref_nodes_by_edges[(*eit)->getParentEnt()] =
              eit->get()->getRefEnt();
          {
            const EntityHandle *conn_edge;
            int num_nodes;
            moab.get_connectivity(edge, conn_edge, num_nodes, true);
            if (conn_edge[0] == edge_new_nodes[ee]) {
              SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                      "node 0 on the edges is mid node, that make no sense");
            }
            if (conn_edge[1] == edge_new_nodes[ee]) {
              SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                      "node 1 on the edges is mid node, that make no sense");
            }
          }
          split_edges[parent_edges_bit.count()] = ee;
          parent_edges_bit.set(ee, 1);
        }
      }
    }

    // test if nodes used to refine are not part of tet
    for (int ee = 0; ee != 6; ee++) {
      if (edge_new_nodes[ee] == no_handle)
        continue;
      for (int nn = 0; nn != 4; nn++) {
        if (conn[nn] == edge_new_nodes[ee]) {
          std::cerr << "problem on edge: " << ee << endl;
          std::cerr << "tet conn : " << conn[0] << " " << conn[1] << " "
                    << conn[2] << " " << conn[3] << " : " << edge_new_nodes[ee]
                    << std::endl;
          for (int eee = 0; eee != 6; ++eee) {
            EntityHandle edge;
            CHKERR moab.side_element(*tit, 1, eee, edge);
            const EntityHandle *conn_edge;
            int num_nodes;
            CHKERR moab.get_connectivity(edge, conn_edge, num_nodes, true);
            std::cerr << eee << " : " << conn_edge[0] << " " << conn_edge[1]
                      << std::endl;
          }
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "nodes used to refine are not part of tet");
        }
      }
    }

    // swap nodes forward
    EntityHandle _conn_[4];
    std::copy(&conn[0], &conn[4], &_conn_[0]);
    // build connectivity for rf tets
    EntityHandle new_tets_conns[8 * 4];
    std::fill(&new_tets_conns[0], &new_tets_conns[8 * 4], no_handle);
    int sub_type = -1, nb_new_tets = 0;
    switch (parent_edges_bit.count()) {
    case 0: {
      RefEntity_multiIndex::iterator tit_miit;
      tit_miit = refined_ents_ptr->find(*tit);
      if (tit_miit == refined_ents_ptr->end())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "can not find this tet");
      ents_to_set_bit.insert(*tit);
      continue;
    } break;
    case 1:
      sub_type = 0;
      // for(int nn = 0;nn<4;nn++) {
      //   if(_conn_[nn] == edge_new_nodes[split_edges[0]]) {
      //     std::cerr << _conn_[0] << " "
      //     << _conn_[1] << " "
      //     << _conn_[2] << " "
      //     << _conn_[3] << " : "
      //     << edge_new_nodes[split_edges[0]] << std::endl;
      //     SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data
      //     inconsistency");
      //   }
      // }
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
      assert(0);
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
      assert(0);
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
      assert(0);
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
      assert(0);
    }
    // find that tets
    EntityHandle ref_tets[8];
    std::bitset<8> ref_tets_bit(0);
    RefEntByParentAndRefEdges::iterator it_by_ref_edges =
        ref_ele_by_parent_and_ref_edges.lower_bound(
            boost::make_tuple(*tit, parent_edges_bit.to_ulong()));
    RefEntByParentAndRefEdges::iterator hi_it_by_ref_edges =
        ref_ele_by_parent_and_ref_edges.upper_bound(
            boost::make_tuple(*tit, parent_edges_bit.to_ulong()));
    // check if tet with this refinement shame already exits
    if (std::distance(it_by_ref_edges, hi_it_by_ref_edges) ==
        (unsigned int)nb_new_tets) {
      for (int tt = 0; it_by_ref_edges != hi_it_by_ref_edges;
           it_by_ref_edges++, tt++) {
        EntityHandle tet = it_by_ref_edges->get()->getRefEnt();
        // set ref tets entities
        ref_tets[tt] = tet;
        ref_tets_bit.set(tt, 1);
        // add this tet if exist to this ref level
        RefEntity_multiIndex::iterator ref_tet_it;
        ref_tet_it = refined_ents_ptr->find(tet);
        if (ref_tet_it == refined_ents_ptr->end()) {
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "data inconsistency");
        }
        ents_to_set_bit.insert(tet);
      }
    } else {
      // if this element was not refined or was refined with different patterns
      // of split edges create new elements
      for (int tt = 0; tt < nb_new_tets; tt++) {
        if (!ref_tets_bit.test(tt)) {
          CHKERR moab.create_element(MBTET, &new_tets_conns[4 * tt], 4,
                                     ref_tets[tt]);
          ents_to_set_bit.insert(ref_tets[tt]);
          Range ref_tet_nodes;
          CHKERR moab.get_connectivity(&ref_tets[tt], 1, ref_tet_nodes, true);
          int ref_type[2];
          ref_type[0] = parent_edges_bit.count();
          ref_type[1] = sub_type;
          CHKERR moab.tag_set_data(cOre.get_th_RefType(), &ref_tets[tt], 1,
                                   ref_type);
          CHKERR moab.tag_set_data(cOre.get_th_RefBitEdge(), &ref_tets[tt], 1,
                                   &parent_edges_bit);
          CHKERR moab.tag_set_data(cOre.get_th_RefParentHandle(), &ref_tets[tt],
                                   1, &*tit);
          // set bit that this element is now in database
          ref_tets_bit.set(tt);
        }
      }
    }
    // find parents for new edges and faces
    // get tet edges and faces
    Range tit_edges, tit_faces;
    CHKERR moab.get_adjacencies(&*tit, 1, 1, false, tit_edges);
    CHKERR moab.get_adjacencies(&*tit, 1, 2, false, tit_faces);
    Range edges_nodes[6], faces_nodes[4];
    // for edges - add ref nodes
    // edges_nodes[ee] - contains all nodes on edge ee including mid nodes if
    // exist
    Range::iterator eit = tit_edges.begin();
    for (int ee = 0; eit != tit_edges.end(); eit++, ee++) {
      CHKERR moab.get_connectivity(&*eit, 1, edges_nodes[ee], true);
      std::map<EntityHandle, EntityHandle>::iterator map_miit =
          map_ref_nodes_by_edges.find(*eit);
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
        std::map<EntityHandle, EntityHandle>::iterator map_miit =
            map_ref_nodes_by_edges.find(*eit2);
        if (map_miit != map_ref_nodes_by_edges.end()) {
          faces_nodes[ff].insert(map_miit->second);
        }
      }
    }
    // add ref nodes to tet
    // tet_nodes contains all nodes on tet including mid edge nodes
    Range tet_nodes;
    CHKERR moab.get_connectivity(&*tit, 1, tet_nodes, true);
    for (std::map<EntityHandle, EntityHandle>::iterator map_miit =
             map_ref_nodes_by_edges.begin();
         map_miit != map_ref_nodes_by_edges.end(); map_miit++) {
      tet_nodes.insert(map_miit->second);
    }
    Range ref_edges;
    // Get all all edges of refined tets
    CHKERR moab.get_adjacencies(ref_tets, nb_new_tets, 1, true, ref_edges,
                                moab::Interface::UNION);
    // Check for all ref edge and set parents
    for (Range::iterator reit = ref_edges.begin(); reit != ref_edges.end();
         reit++) {
      Range ref_edges_nodes;
      CHKERR moab.get_connectivity(&*reit, 1, ref_edges_nodes, true);
      if (ref_edges_nodes.size() != 2) {
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "data inconsistency, edge should have 2 nodes");
      }
      // Check if ref edge is an coarse edge (loop over coarse tet edges)
      int ee = 0;
      for (; ee < 6; ee++) {
        // Two nodes are common (node[0],node[1],ref_node (if exist))
        // this tests if given edge is contained by edge of refined tetrahedral
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

      // check if ref edge is in coarse tetrahedral (i.e. that is internal edge
      // of refined tetrahedral)
      if (intersect(tet_nodes, ref_edges_nodes).size() == 2) {
        CHKERR set_parent(*reit, *tit, refined_ents_ptr, cOre);
        continue;
      }

      // Refined edge is not child of any edge, face or tetrahedral, this is
      // imposible edge
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "impossible refined edge");
    }

    Range ref_faces;
    CHKERR moab.get_adjacencies(ref_tets, nb_new_tets, 2, true, ref_faces,
                                moab::Interface::UNION);
    Tag th_interface_side;
    const int def_side[] = {0};
    CHKERR moab.tag_get_handle("INTERFACE_SIDE", 1, MB_TYPE_INTEGER,
                               th_interface_side, MB_TAG_CREAT | MB_TAG_SPARSE,
                               def_side);
    // Check for all ref faces
    for (Range::iterator rfit = ref_faces.begin(); rfit != ref_faces.end();
         rfit++) {
      Range ref_faces_nodes;
      CHKERR moab.get_connectivity(&*rfit, 1, ref_faces_nodes, true);
      // Check if ref face is in coarse face
      int ff = 0;
      for (; ff < 4; ff++) {
        // Check if refined triangle is contained by face of tetrahedral
        if (intersect(faces_nodes[ff], ref_faces_nodes).size() == 3) {
          EntityHandle face = tit_faces[ff];
          CHKERR set_parent(*rfit, face, refined_ents_ptr, cOre);
          int side = 0;
          // Set face side if it is on interface
          CHKERR moab.tag_get_data(th_interface_side, &face, 1, &side);
          CHKERR moab.tag_set_data(th_interface_side, &*rfit, 1, &side);
          break;
        }
      }
      if (ff < 4)
        continue; // this face is contained by one of tetrahedrons
      // check if ref face is in coarse tetrahedral
      // this is ref face which is contained by tetrahedral volume
      if (intersect(tet_nodes, ref_faces_nodes).size() == 3) {
        CHKERR set_parent(*rfit, *tit, refined_ents_ptr, cOre);
        continue;
      }
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "impossible refined face");
    }
  }

  CHKERR check(ref_edges_ptr, refined_ents_ptr, cOre);
  CHKERR set_parent(refined_ents_ptr);
  CHKERR check(ref_edges_ptr, refined_ents_ptr, cOre);
  CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevel(ents_to_set_bit,
                                                               bit, true, verb);
  CHKERR check(ref_edges_ptr, refined_ents_ptr, cOre);

  MoFEMFunctionReturn(0);
}
MoFEMErrorCode MeshRefinement::refine_PRISM(const EntityHandle meshset,
                                            const BitRefLevel &bit, int verb) {

  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  const RefEntity_multiIndex *refined_ents_ptr;
  const RefElement_multiIndex *refined_finite_elements_ptr;

  // FIXME: refinement is based on entity handlers, should work on global ids of
  // nodes, this will allow parallelise algorithm in the future

  MoFEMFunctionBeginHot;

  ierr = m_field.get_ref_ents(&refined_ents_ptr);
  CHKERRG(ierr);
  ierr = m_field.get_ref_finite_elements(&refined_finite_elements_ptr);

  typedef const RefEntity_multiIndex::index<Ent_mi_tag>::type RefEntsByEnt;

  typedef const RefElement_multiIndex_parents_view::index<
      Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>::type
      RefEntByParentAndRefEdges;
  RefElement_multiIndex_parents_view ref_ele_parent_view;
  ref_ele_parent_view.insert(
      refined_finite_elements_ptr->get<EntType_mi_tag>().lower_bound(MBPRISM),
      refined_finite_elements_ptr->get<EntType_mi_tag>().upper_bound(MBPRISM));
  RefEntByParentAndRefEdges &ref_ele_by_parent_and_ref_edges =
      ref_ele_parent_view
          .get<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>();
  // find all vertices which parent is edge
  typedef const RefEntity_multiIndex::index<
      Composite_EntType_and_ParentEntType_mi_tag>::type RefEntsByComposite;
  RefEntsByComposite &ref_ents_by_comp =
      refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>();
  RefEntity_multiIndex_view_by_hashed_parent_entity ref_parent_ents_view;
  ref_parent_ents_view.insert(
      ref_ents_by_comp.lower_bound(boost::make_tuple(MBVERTEX, MBEDGE)),
      ref_ents_by_comp.upper_bound(boost::make_tuple(MBVERTEX, MBEDGE)));
  Range prisms;
  rval = moab.get_entities_by_type(meshset, MBPRISM, prisms, false);
  CHKERRQ_MOAB(rval);
  Range::iterator pit = prisms.begin();
  for (; pit != prisms.end(); pit++) {
    RefEntsByEnt::iterator miit_prism =
        refined_ents_ptr->get<Ent_mi_tag>().find(*pit);
    if (miit_prism == refined_ents_ptr->end()) {
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "this prism is not in ref database");
    }
    if (verb >= NOISY) {
      std::ostringstream ss;
      ss << "ref prism " << **miit_prism << std::endl;
      PetscPrintf(m_field.get_comm(), ss.str().c_str());
    }
    // prism connectivity
    int num_nodes;
    const EntityHandle *conn;
    rval = moab.get_connectivity(*pit, conn, num_nodes, true);
    CHKERRQ_MOAB(rval);
    assert(num_nodes == 6);
    // edges connectivity
    EntityHandle edges[6];
    for (int ee = 0; ee < 3; ee++) {
      rval = moab.side_element(*pit, 1, ee, edges[ee]);
      CHKERRQ_MOAB(rval);
    }
    for (int ee = 6; ee < 9; ee++) {
      rval = moab.side_element(*pit, 1, ee, edges[ee - 3]);
      CHKERRQ_MOAB(rval);
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
          edge_nodes[ee] = (*miit_view)->getRefEnt();
          split_edges.set(ee);
        }
      }
    }
    if (split_edges.count() == 0) {
      *(const_cast<RefEntity *>(miit_prism->get())->getBitRefLevelPtr()) |= bit;
      if (verb >= VERY_NOISY)
        PetscPrintf(m_field.get_comm(), "no refinement");
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
      ierr = prism_type_1(conn, split_edges, edge_nodes, new_prism_conn);
      CHKERRG(ierr);
      nb_new_prisms = 2;
      break;
    case 4:
      ierr = prism_type_2(conn, split_edges, edge_nodes, new_prism_conn);
      CHKERRG(ierr);
      nb_new_prisms = 3;
      break;
    case 6:
      ierr = prism_type_3(conn, split_edges, edge_nodes, new_prism_conn);
      CHKERRG(ierr);
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
    RefEntByParentAndRefEdges::iterator it_by_ref_edges =
        ref_ele_by_parent_and_ref_edges.lower_bound(
            boost::make_tuple(*pit, split_edges.to_ulong()));
    RefEntByParentAndRefEdges::iterator hi_it_by_ref_edges =
        ref_ele_by_parent_and_ref_edges.upper_bound(
            boost::make_tuple(*pit, split_edges.to_ulong()));
    RefEntByParentAndRefEdges::iterator it_by_ref_edges2 = it_by_ref_edges;
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
          rval = moab.create_element(MBPRISM, &new_prism_conn[6 * pp], 6,
                                     ref_prisms[pp]);
          CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(cOre.get_th_RefParentHandle(),
                                   &ref_prisms[pp], 1, &*pit);
          CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(cOre.get_th_RefBitLevel(), &ref_prisms[pp],
                                   1, &bit);
          CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(cOre.get_th_RefBitEdge(), &ref_prisms[pp], 1,
                                   &split_edges);
          CHKERRQ_MOAB(rval);
          std::pair<RefEntity_multiIndex::iterator, bool> p_ent =
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
          ierr = cOre.addPrismToDatabase(ref_prisms[pp]);
          CHKERRG(ierr);
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
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode MeshRefinement::refine_MESHSET(const EntityHandle meshset,
                                              const BitRefLevel &bit,
                                              const bool recursive, int verb) {
  Interface &m_field = cOre;
  const RefEntity_multiIndex *refined_ents_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_ref_ents(&refined_ents_ptr);
  typedef const RefEntity_multiIndex::index<Ent_mi_tag>::type RefEntsByEnt;
  RefEntsByEnt::iterator miit = refined_ents_ptr->find(meshset);
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

} // namespace MoFEM
