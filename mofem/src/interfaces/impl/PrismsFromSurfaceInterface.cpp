/** \file PrismsFromSurfaceInterface.cpp
 * \brief Interface for creating prisms from surface elements
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

MoFEMErrorCode PrismsFromSurfaceInterface::query_interface(
    boost::typeindex::type_index type_index, UnknownInterface **iface) const {
  *iface = const_cast<PrismsFromSurfaceInterface *>(this);
  return 0;
}

DEPRECATED MoFEMErrorCode PrismsFromSurfaceInterface::createPrisms(
    const Range &ents, Range &prisms, int verb) {
  MoFEMFunctionBegin;
  CHKERR createPrisms(ents, NO_SWAP, prisms, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PrismsFromSurfaceInterface::createPrisms(
    const Range &ents, const SwapType swap_type, Range &prisms, int verb) {
  MoFEMFunctionBegin;
  Interface &m_field = cOre;
  Range tris = ents.subset_by_type(MBTRI);
  for (Range::iterator tit = tris.begin(); tit != tris.end(); tit++) {
    const EntityHandle *conn;
    int number_nodes = 0;
    CHKERR m_field.get_moab().get_connectivity(*tit, conn, number_nodes, false);
    double coords[3 * number_nodes];
    CHKERR m_field.get_moab().get_coords(conn, number_nodes, coords);
    EntityHandle prism_nodes[6];
    for (int nn = 0; nn < 3; nn++) {
      prism_nodes[nn] = conn[nn];
      if (createdVertices.find(conn[nn]) != createdVertices.end()) {
        prism_nodes[3 + nn] = createdVertices[prism_nodes[nn]];
      } else {
        CHKERR m_field.get_moab().create_vertex(&coords[3 * nn],
                                                prism_nodes[3 + nn]);
        createdVertices[conn[nn]] = prism_nodes[3 + nn];
        CHKERR m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(),
                                               &prism_nodes[3 + nn], 1,
                                               &prism_nodes[nn]);
      }
    }
    switch (swap_type) {
    case SWAP_TRI_NODE_ORDER:
      std::swap(prism_nodes[1], prism_nodes[2]);
      std::swap(prism_nodes[4], prism_nodes[5]);
      break;
    case SWAP_TOP_AND_BOT_TRI:
      std::swap(prism_nodes[0], prism_nodes[3]);
      std::swap(prism_nodes[1], prism_nodes[4]);
      std::swap(prism_nodes[2], prism_nodes[5]);
      break;
    case NO_SWAP:
    default:
      break;
    }
    EntityHandle prism;
    CHKERR m_field.get_moab().create_element(MBPRISM, prism_nodes, 6, prism);
    Range edges;
    CHKERR m_field.get_moab().get_adjacencies(&prism, 1, 1, true, edges,
                                              moab::Interface::UNION);
    Range faces;
    CHKERR m_field.get_moab().get_adjacencies(&prism, 1, 2, true, faces,
                                              moab::Interface::UNION);
    prisms.insert(prism);
    for (int ee = 0; ee <= 2; ee++) {
      EntityHandle e1;
      CHKERR m_field.get_moab().side_element(prism, 1, ee, e1);
      EntityHandle e2;
      CHKERR m_field.get_moab().side_element(prism, 1, ee + 6, e2);
      CHKERR m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(), &e2,
                                             1, &e1);
    }
    EntityHandle f3, f4;
    {
      CHKERR m_field.get_moab().side_element(prism, 2, 3, f3);
      CHKERR m_field.get_moab().side_element(prism, 2, 4, f4);
      CHKERR m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(), &f4,
                                             1, &f3);
    }
    if (number_nodes > 3) {
      EntityHandle meshset;
      CHKERR m_field.get_moab().create_meshset(MESHSET_SET, meshset);
      CHKERR m_field.get_moab().add_entities(meshset, &f4, 1);
      for (int ee = 0; ee <= 2; ee++) {
        EntityHandle e2;
        CHKERR m_field.get_moab().side_element(prism, 1, ee + 6, e2);
        CHKERR m_field.get_moab().add_entities(meshset, &e2, 1);
      }
      CHKERR m_field.get_moab().convert_entities(meshset, true, false, false);
      CHKERR m_field.get_moab().delete_entities(&meshset, 1);
      const EntityHandle *conn_f4;
      int number_nodes_f4 = 0;
      CHKERR m_field.get_moab().get_connectivity(f4, conn_f4, number_nodes_f4,
                                                 false);
      if (number_nodes_f4 != number_nodes) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
      CHKERR m_field.get_moab().set_coords(&conn_f4[3], 3, &coords[9]);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PrismsFromSurfaceInterface::seedPrismsEntities(
    Range &prisms, const BitRefLevel &bit, int verb) {
  Interface &m_field = cOre;
  auto ref_ents_ptr = m_field.get_ref_ents();
  MoFEMFunctionBegin;
  MPI_Comm comm = m_field.get_comm();
  RefEntity_multiIndex *refined_entities_ptr;
  refined_entities_ptr =
      const_cast<RefEntity_multiIndex *>(ref_ents_ptr);
  if (!prisms.empty()) {
    int dim = m_field.get_moab().dimension_from_handle(prisms[0]);
    for (int dd = 0; dd <= dim; dd++) {
      Range ents;
      CHKERR m_field.get_moab().get_adjacencies(prisms, dd, true, ents,
                                                moab::Interface::UNION);
      Range::iterator eit = ents.begin();
      for (; eit != ents.end(); eit++) {
        std::pair<RefEntity_multiIndex::iterator, bool> p_ent =
            refined_entities_ptr->insert(boost::shared_ptr<RefEntity>(
                new RefEntity(m_field.get_basic_entity_data_ptr(), *eit)));
        *(const_cast<RefEntity *>(p_ent.first->get())->getBitRefLevelPtr()) |=
            bit;
        if (verb >= VERY_VERBOSE) {
          std::ostringstream ss;
          ss << *(p_ent.first);
          PetscSynchronizedPrintf(comm, "%s\n", ss.str().c_str());
        }
      }
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PrismsFromSurfaceInterface::createPrismsFromPrisms(
    const Range &prisms, bool from_down, Range &out_prisms, int verb) {
  MoFEMFunctionBegin;
  Interface &m_field = cOre;
  Range tris;
  for (Range::iterator pit = prisms.begin(); pit != prisms.end(); pit++) {
    EntityHandle face;
    if (from_down) {
      CHKERR m_field.get_moab().side_element(*pit, 2, 3, face);
    } else {
      CHKERR m_field.get_moab().side_element(*pit, 2, 4, face);
    }
    tris.insert(face);
  }
  CHKERR createPrisms(tris, out_prisms, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PrismsFromSurfaceInterface::setThickness(
    const Range &prisms, const double director3[], const double director4[]) {
  MoFEMFunctionBegin;
  Interface &m_field = cOre;
  Range nodes_f3, nodes_f4;
  for (Range::iterator pit = prisms.begin(); pit != prisms.end(); pit++) {
    for (int ff = 3; ff <= 4; ff++) {
      EntityHandle face;
      CHKERR m_field.get_moab().side_element(*pit, 2, ff, face);
      const EntityHandle *conn;
      int number_nodes = 0;
      CHKERR m_field.get_moab().get_connectivity(face, conn, number_nodes,
                                                 false);
      if (ff == 3) {
        nodes_f3.insert(&conn[0], &conn[number_nodes]);
      } else {
        nodes_f4.insert(&conn[0], &conn[number_nodes]);
      }
    }
  }
  double coords[3];
  for (Range::iterator nit = nodes_f3.begin(); nit != nodes_f3.end(); nit++) {
    CHKERR m_field.get_moab().get_coords(&*nit, 1, coords);
    cblas_daxpy(3, 1, director3, 1, coords, 1);
    CHKERR m_field.get_moab().set_coords(&*nit, 1, coords);
  }
  for (Range::iterator nit = nodes_f4.begin(); nit != nodes_f4.end(); nit++) {
    CHKERR m_field.get_moab().get_coords(&*nit, 1, coords);
    cblas_daxpy(3, 1, director4, 1, coords, 1);
    CHKERR m_field.get_moab().set_coords(&*nit, 1, coords);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PrismsFromSurfaceInterface::setNormalThickness(
    const Range &prisms, double thickness3, double thickness4) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  auto add_normal = [&](std::map<EntityHandle, std::array<double, 3>> &nodes,
                        EntityHandle face) {
    MoFEMFunctionBeginHot;
    const EntityHandle *conn;
    int number_nodes;
    CHKERR m_field.get_moab().get_connectivity(face, conn, number_nodes, false);
    std::array<double, 9> coords;
    CHKERR m_field.get_moab().get_coords(conn, number_nodes, coords.data());
    std::array<double, 3> normal;
    CHKERR Tools::getTriNormal(coords.data(), normal.data());
    double a = sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                    normal[2] * normal[2]);
    for (auto d : {0, 1, 2})
      normal[d] /= a;
    for (auto n : {0, 1, 2}) {
      try {
        for (auto d : {0, 1, 2})
          nodes.at(conn[n])[d] += normal[d];
      } catch (...) {
        nodes.insert(
            std::pair<EntityHandle, std::array<double, 3>>(conn[n], normal));
      }
    }
    MoFEMFunctionReturnHot(0);
  };

  auto apply_map = [&](auto &nodes, double t) {
    MoFEMFunctionBeginHot;
    for (auto &m : nodes) {
      std::array<double, 3> coords;
      CHKERR m_field.get_moab().get_coords(&m.first, 1, coords.data());
      auto &normal = m.second;
      double a = sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                      normal[2] * normal[2]);
      for (auto d : {0, 1, 2})
        coords[d] += (normal[d] / a) * t;
      CHKERR m_field.get_moab().set_coords(&m.first, 1, coords.data());
    }
    MoFEMFunctionReturnHot(0);
  };

  map<EntityHandle, std::array<double, 3>> nodes_f3, nodes_f4;
  for (Range::iterator pit = prisms.begin(); pit != prisms.end(); pit++) {
    for (int ff = 3; ff <= 4; ff++) {
      EntityHandle face;
      CHKERR m_field.get_moab().side_element(*pit, 2, ff, face);
      if (ff == 3)
        CHKERR add_normal(nodes_f3, face);
      else
        CHKERR add_normal(nodes_f4, face);
    }
  }

  CHKERR apply_map(nodes_f3, thickness3);
  CHKERR apply_map(nodes_f4, thickness4);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
PrismsFromSurfaceInterface::updateMeshestByEdgeBlock(const Range &prisms) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  Range prisms_edges;
  CHKERR m_field.get_moab().get_adjacencies(prisms, 1, true, prisms_edges,
                                            moab::Interface::UNION);
  Range prisms_faces;
  CHKERR m_field.get_moab().get_adjacencies(prisms, 2, true, prisms_faces,
                                            moab::Interface::UNION);
  for (_IT_CUBITMESHSETS_FOR_LOOP_(m_field, it)) {
    Range edges;
    CHKERR m_field.get_moab().get_entities_by_type(it->meshset, MBEDGE, edges,
                                                   true);
    edges = intersect(edges, prisms_edges);
    if (!edges.empty()) {
      Range edges_faces;
      CHKERR m_field.get_moab().get_adjacencies(edges, 2, false, edges_faces,
                                                moab::Interface::UNION);
      edges_faces = intersect(edges_faces, prisms_faces.subset_by_type(MBQUAD));
      EntityHandle meshset = it->getMeshset();
      CHKERR m_field.get_moab().add_entities(meshset, edges_faces);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
PrismsFromSurfaceInterface::updateMeshestByTriBlock(const Range &prisms) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  Range prisms_tris;
  CHKERR m_field.get_moab().get_adjacencies(prisms, 2, true, prisms_tris,
                                            moab::Interface::UNION);
  prisms_tris = prisms_tris.subset_by_type(MBTRI);
  for (_IT_CUBITMESHSETS_FOR_LOOP_(m_field, it)) {
    Range tris;
    CHKERR m_field.get_moab().get_entities_by_type(it->meshset, MBTRI, tris,
                                                   true);
    tris = intersect(tris, prisms_tris);
    if (!tris.empty()) {
      Range tris_ents;
      CHKERR m_field.get_moab().get_adjacencies(tris, 3, false, tris_ents,
                                                moab::Interface::UNION);
      tris_ents = intersect(tris_ents, prisms);
      EntityHandle meshset = it->getMeshset();
      CHKERR m_field.get_moab().add_entities(meshset, tris_ents);
    }
  }
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
