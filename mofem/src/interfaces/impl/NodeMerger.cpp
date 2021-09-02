/** \file NodeMerger.cpp
 * \brief Interface for merging nodes
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

namespace MoFEM {

MoFEMErrorCode
NodeMergerInterface::query_interface(boost::typeindex::type_index type_index,
                                     UnknownInterface **iface) const {
  *iface = const_cast<NodeMergerInterface *>(this);
  return 0;
}

static auto min_non_abs(const double a, const double b) {
  return std::min(a, b);
};

NodeMergerInterface::NodeMergerInterface(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)), minQualityFunction(min_non_abs),
      successMerge(false), errorIfNoCommonEdge(false) {}

MoFEMErrorCode NodeMergerInterface::getSubInterfaceOptions() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  ierr = PetscOptionsBegin(m_field.get_comm(), "node_merge",
                           "Node merge options", "none");
  CHKERRQ(ierr);
  
  ierr = PetscOptionsEnd();
  CHKERRQ(ierr);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode NodeMergerInterface::mergeNodes(
    EntityHandle father, EntityHandle mother, Range &out_tets, Range *tets_ptr,
    const bool only_if_improve_quality, const double move,
    const int line_search, Tag th, int verb) {
  Interface &m_field = cOre;
  FTensor::Index<'i', 3> i;
  MoFEMFunctionBegin;

  // Get edges adjacent to father and mother, i.e. mother is merged to father.
  Range father_edges;
  CHKERR m_field.get_moab().get_adjacencies(&father, 1, 1, false, father_edges);
  Range mother_edges;
  CHKERR m_field.get_moab().get_adjacencies(&mother, 1, 1, false, mother_edges);

  // Get tets adjacent to mother and father
  Range father_tets;
  CHKERR m_field.get_moab().get_adjacencies(&father, 1, 3, false, father_tets);
  Range mother_tets;
  CHKERR m_field.get_moab().get_adjacencies(&mother, 1, 3, false, mother_tets);
  if (tets_ptr != NULL) {
    mother_tets = intersect(mother_tets, *tets_ptr);
    father_tets = intersect(father_tets, *tets_ptr);
  }

  // Find common edge
  Range common_edge;
  common_edge = intersect(father_edges, mother_edges);
  if (tets_ptr != NULL) {
    Range tets = unite(father_tets, mother_tets);
    Range tets_edges;
    CHKERR m_field.get_moab().get_adjacencies(tets, 1, false, tets_edges,
                                              moab::Interface::UNION);
    common_edge = intersect(common_edge, tets_edges);
    father_edges = intersect(father_edges, tets_edges);
    mother_edges = intersect(mother_edges, tets_edges);
  }

  // No common edge, merge no possible
  if (errorIfNoCommonEdge && common_edge.empty()) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "no common edge between nodes");
  } else if (common_edge.empty()) {
    Range seed_tets;
    if (tets_ptr != NULL) {
      seed_tets.merge(*tets_ptr);
    }
    out_tets = seed_tets;
    successMerge = false;
    MoFEMFunctionReturnHot(0);
  }

  // Common edge tets, that tests will be squashed
  Range edge_tets;
  CHKERR m_field.get_moab().get_adjacencies(common_edge, 3, true, edge_tets);
  // Intersect with ptr_tets (usually associated with some bit level)
  if (tets_ptr != NULL) {
    edge_tets = intersect(edge_tets, *tets_ptr);
  }
  // Mother tets, has only one mother vertex and no father vertex.
  mother_tets = subtract(mother_tets, edge_tets);
  father_tets = subtract(father_tets, edge_tets);

  auto get_coords = [&m_field](Tag th, const EntityHandle *conn,
                               const int num_nodes) {
    VectorDouble coords(3 * num_nodes);
    if (th == NULL) {
      CHKERR m_field.get_moab().get_coords(conn, num_nodes, &coords[0]);
    } else {
      CHKERR m_field.get_moab().tag_get_data(th, conn, num_nodes, &coords[0]);
    }
    return coords;
  };

  auto get_tensor = [](VectorDouble &coords, const int shift) {
    return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
        &coords[shift], &coords[shift + 1], &coords[shift + 2]);
  };

  // move father coord is move > 0
  FTensor::Tensor1<double, 3> t_move;
  if (move > 0) {
    EntityHandle conn[] = {father, mother};
    VectorDouble coords = get_coords(th, conn, 2);
    auto t_n0 = get_tensor(coords, 0);
    auto t_n1 = get_tensor(coords, 3);
    t_move(i) = t_n0(i) + move * (t_n1(i) - t_n0(i));
  }

  if (line_search > 0) {
    Range check_tests = unite(father_tets, mother_tets);
    CHKERR lineSearch(check_tests, father, mother, line_search, t_move, th);
  }

  if (only_if_improve_quality) {
    Range check_tests = mother_tets;
    // no need to check father tets since no change of quality for them
    if (move > 0 || line_search) {
      check_tests.merge(father_tets);
    }

    double min_quality0 = 1;
    CHKERR minQuality(edge_tets, 0, 0, NULL, min_quality0, th,
                      minQualityFunction);
    CHKERR minQuality(check_tests, 0, 0, NULL, min_quality0, th,
                      minQualityFunction);
    double min_quality = 1;
    CHKERR minQuality(check_tests, father, mother,
                      ((move > 0) || line_search) ? &t_move(0) : NULL,
                      min_quality, th, minQualityFunction);
    if (min_quality < min_quality0) {
      if (tets_ptr != NULL) {
        out_tets = *tets_ptr;
      }
      successMerge = false;
      MoFEMFunctionReturnHot(0);
    }
  }

  // Move node
  if (move > 0 || line_search) {
    if (th == NULL) {
      CHKERR m_field.get_moab().set_coords(&father, 1, &t_move(0));
    } else {
      CHKERR m_field.get_moab().tag_set_data(th, &father, 1, &t_move(0));
    }
  }

  auto get_conn = [&m_field](const EntityHandle ent,
                             int *ret_num_nodes = nullptr) {
    int num_nodes;
    const EntityHandle *conn;
    CHKERR m_field.get_moab().get_connectivity(ent, conn, num_nodes, true);
    if (ret_num_nodes)
      *ret_num_nodes = num_nodes;
    return conn;
  };

  auto create_tet = [this, &m_field](const EntityHandle *new_conn,
                                     const EntityHandle parent) {
    EntityHandle tet;
    Range tets;
    CHKERR m_field.get_moab().get_adjacencies(new_conn, 4, 3, false, tets);
    bool tet_found = false;
    for (auto it_tet : tets) {
      const EntityHandle *tet_conn;
      int num_nodes;
      CHKERR m_field.get_moab().get_connectivity(it_tet, tet_conn, num_nodes,
                                                 true);
      const EntityHandle *p = std::find(tet_conn, &tet_conn[4], new_conn[0]);
      if (p != tet_conn + 4) {
        int s = std::distance(tet_conn, p);
        int n = 0;
        for (; n != 4; ++n) {
          const int idx[] = {0, 1, 2, 3, 0, 1, 2, 3};
          if (tet_conn[idx[s + n]] != new_conn[n])
            break;
        }
        if (n == 4 && !tet_found) {
          tet = it_tet;
          tet_found = true;
        } else if (n == 4) {
          THROW_MESSAGE("More that one tet with the same connectivity");
        }
      }
    }
    if (!tet_found) {
      // Create tet with new connectivity
      CHKERR m_field.get_moab().create_element(MBTET, new_conn, 4, tet);
      CHKERR m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(),
                                             &tet, 1, &parent);
      parentChildMap.insert(ParentChild(parent, tet));
    }
    return tet;
  };

  // clear map
  parentChildMap.clear();

  Range created_tets;
  for (auto m_tet : mother_tets) {
    const EntityHandle *conn = get_conn(m_tet);
    EntityHandle new_conn[4];
    // Replace mother vertices by father vertices
    int nb_mother_verts = 0;
    for (int nn = 0; nn != 4; ++nn) {
      if (conn[nn] == father) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Tet has father vertex, impossible but here it is");
      }
      if (conn[nn] == mother) {
        new_conn[nn] = father;
        nb_mother_verts++;
      } else {
        new_conn[nn] = conn[nn];
      }
    }
    if (nb_mother_verts != 1) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Tet should have only one vertex but have %d", nb_mother_verts);
    }

    VectorDouble new_coords = get_coords(th, new_conn, 4);

    // add tet to range
    created_tets.insert(create_tet(new_conn, m_tet));
  }

  // Loop over mother adjacent entities to use them as parents
  auto get_adj_ents = [&](const Range &ents, const bool create) {
    Range adj;
    for (int dd = 1; dd <= 2; dd++)
      CHKERR m_field.get_moab().get_adjacencies(ents, dd, create, adj,
                                                moab::Interface::UNION);
    return adj;
  };
  auto adj_crated_ents = get_adj_ents(created_tets, true);
  adj_crated_ents.erase(common_edge[0]);

  FaceMapIdx face_map;
  for (auto ent : adj_crated_ents) {
    int num_nodes;
    const EntityHandle *conn = get_conn(ent, &num_nodes);
    EntityHandle small_conn[num_nodes];
    int ii = 0;
    bool father_node = false;
    for (int nn = 0; nn != num_nodes; nn++) {
      if (conn[nn] == father)
        father_node = true;
      else
        small_conn[ii++] = conn[nn];
    }
    if (father_node) {
      if (ii > 1) 
        std::sort(&small_conn[0], &small_conn[ii]);
      if (ii == 2) 
        face_map.insert(FaceMap(ent, small_conn[0], small_conn[1]));
       else 
        face_map.insert(FaceMap(ent, small_conn[0], 0));
    }
  }

  auto adj_mother_ents = get_adj_ents(mother_tets, false);
  adj_mother_ents.erase(common_edge[0]);
  for (auto ent : adj_mother_ents) {
    int num_nodes;
    const EntityHandle *conn = get_conn(ent, &num_nodes);
    EntityHandle small_conn[num_nodes];
    int nb_new_node = 0;
    int nn = 0;
    int ii = 0;
    for (; nn != num_nodes; ++nn) {
      if (conn[nn] == mother) {
        nb_new_node++;
      } else {
        small_conn[ii++] = conn[nn];
      }
    }
    if (nb_new_node > 0) {
      if (ii > 1) 
        std::sort(&small_conn[0], &small_conn[ii]);
      
      EntityHandle n0 = small_conn[0], n1 = 0;
      if (ii == 2) 
        n1 = small_conn[1];
      
      FaceMapIdx::iterator fit = face_map.find(boost::make_tuple(n0, n1));
      if (fit != face_map.end()) {
        const EntityHandle child = fit->e;
        const EntityHandle parent = ent;
        if (m_field.get_moab().dimension_from_handle(parent) !=
            m_field.get_moab().dimension_from_handle(child))
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Huston we have a problem!");

        // Set parent child relation
        CHKERR m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(),
                                               &child, 1, &parent);
        // Create map
        parentChildMap.insert(ParentChild(parent, child));
      }
    }
  }

  // Seed tets to given bit level
  Range seed_tets;
  if (tets_ptr != NULL) {
    seed_tets.merge(*tets_ptr);
    seed_tets = subtract(seed_tets, edge_tets);
    seed_tets = subtract(seed_tets, mother_tets);
  }
  seed_tets.merge(created_tets);
  out_tets.swap(seed_tets);

  successMerge = true;

  if (verb > VERY_VERBOSE) 
    std::cerr << "Nodes merged" << endl;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
NodeMergerInterface::minQuality(Range &check_tests, EntityHandle father,
                                EntityHandle mother, double *coords_move,
                                double &min_quality, Tag th,
                                boost::function<double(double, double)> f) {
  Interface &m_field = cOre;
  double coords[12];
  FTensor::Index<'i', 3> i;
  MoFEMFunctionBegin;
  for (Range::iterator tit = check_tests.begin(); tit != check_tests.end();
       tit++) {
    const EntityHandle *conn;
    int num_nodes;
    CHKERR m_field.get_moab().get_connectivity(*tit, conn, num_nodes, true);
    if (mother > 0) {
      EntityHandle new_conn[4];
      // Replace mother vertices by father vertices
      int nb_mother_verts = 0;
      int father_nn = 0;
      for (int nn = 0; nn < 4; nn++) {
        if (conn[nn] == father) {
          father_nn = nn;
        }
        if (conn[nn] == mother) {
          new_conn[nn] = father;
          father_nn = nn;
          nb_mother_verts++;
        } else {
          new_conn[nn] = conn[nn];
        }
      }
      if (nb_mother_verts > 1) {
        SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Tet should have no more than one mother vertex but have %d",
                 nb_mother_verts);
      }
      if (th == NULL) {
        CHKERR m_field.get_moab().get_coords(new_conn, num_nodes, coords);
      } else {
        CHKERR m_field.get_moab().tag_get_data(th, new_conn, num_nodes, coords);
      }
      if (coords_move) {
        int shift = 3 * father_nn;
        for (int nn = 0; nn != 3; nn++) {
          coords[shift + nn] = coords_move[nn];
        }
      }
    } else {
      if (th == NULL) {
        CHKERR m_field.get_moab().get_coords(conn, num_nodes, coords);
      } else {
        CHKERR m_field.get_moab().tag_get_data(th, conn, num_nodes, coords);
      }
    }
    double quality = Tools::volumeLengthQuality(coords);
    if (!std::isnormal(quality)) {
      min_quality = 0;
      MoFEMFunctionReturnHot(0);
    }
    min_quality = f(min_quality, quality);
  }
  MoFEMFunctionReturn(0);
};

MoFEMErrorCode
NodeMergerInterface::lineSearch(Range &check_tests, EntityHandle father,
                                EntityHandle mother, int line_search,
                                FTensor::Tensor1<double, 3> &t_move, Tag th) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  EntityHandle conn[] = {father, mother};

  double coords[6];
  if (th == NULL) {
    CHKERR m_field.get_moab().get_coords(conn, 2, coords);
  } else {
    CHKERR m_field.get_moab().tag_get_data(th, conn, 2, coords);
  }

  FTensor::Index<'i', 3> i;
  FTensor::Tensor1<double, 3> t_coords(coords[0], coords[1], coords[2]);
  FTensor::Tensor1<double, 3> t_delta;
  for (int nn = 0; nn != 3; nn++) {
    t_delta(nn) = coords[3 + nn] - t_coords(nn);
  }

  t_move(i) = t_coords(i);
  double min_quality_i = 1;
  CHKERR minQuality(check_tests, father, mother, &t_move(0), min_quality_i, th,
                    minQualityFunction);

  t_move(i) = t_coords(i) + t_delta(i);
  double min_quality_k = 1;
  CHKERR minQuality(check_tests, father, mother, &t_move(0), min_quality_k, th,
                    minQualityFunction);

  double alpha_i = 0;
  double alpha_k = 1;

  for (int ii = 0; ii != line_search; ii++) {
    double min_quality = 1;
    double alpha = (alpha_i + alpha_k) * 0.5;
    t_move(i) = t_coords(i) + alpha * t_delta(i);
    CHKERR minQuality(check_tests, father, mother, &t_move(0), min_quality, th,
                      minQualityFunction);
    if (min_quality_i >= min_quality_k) {
      min_quality_k = min_quality;
      alpha_k = alpha;
    } else {
      min_quality_i = min_quality;
      alpha_i = alpha;
    }
  }

  if (min_quality_i > min_quality_k) {
    t_move(i) = t_coords(i) + alpha_i * t_delta(i);
  } else {
    t_move(i) = t_coords(i) + alpha_k * t_delta(i);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode NodeMergerInterface::mergeNodes(
    EntityHandle father, EntityHandle mother, BitRefLevel bit, Range *tets_ptr,
    const bool only_if_improve_quality, const double move, Tag th) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  Range out_tets;
  CHKERR mergeNodes(father, mother, out_tets, tets_ptr, only_if_improve_quality,
                    move, 0, th);
  CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevel(out_tets, bit);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode NodeMergerInterface::mergeNodes(
    EntityHandle father, EntityHandle mother, BitRefLevel bit,
    BitRefLevel tets_from_bit_ref_level, const bool only_if_improve_quality,
    const double move, Tag th) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  Range level_tets;
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      tets_from_bit_ref_level, BitRefLevel().set(), MBTET, level_tets);
  CHKERR mergeNodes(father, mother, bit, &level_tets, only_if_improve_quality,
                    move, th);
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
