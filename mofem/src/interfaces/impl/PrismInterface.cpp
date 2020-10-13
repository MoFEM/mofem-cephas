/** \file PrismInterface.cpp
 * \brief Insert prisms in the interface between two surfaces
 * \todo FIXME this is not so good implementation
 *
 * \ingroup mofem_prism_interface
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.

 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.

 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

namespace MoFEM {

MoFEMErrorCode PrismInterface::query_interface(const MOFEMuuid &uuid,
                                               UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMPrismInterface) {
    *iface = const_cast<PrismInterface *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

PrismInterface::PrismInterface(const Core &core)
    : cOre(const_cast<Core &>(core)) {

  try {
    MoFEM::LogManager::getLog("PRISM_INTERFACE");
  } catch (...) {
    auto core_log = logging::core::get();
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "PRISM_INTERFACE"));
    LogManager::setLog("PRISM_INTERFACE");
    MOFEM_LOG_TAG("PRISM_INTERFACE", "PrismInterface");
  }

  MOFEM_LOG("PRISM_INTERFACE", Sev::noisy) << "Prism interface created";
}

MoFEMErrorCode PrismInterface::getSides(const int msId,
                                        const CubitBCType cubit_bc_type,
                                        const BitRefLevel mesh_bit_level,
                                        const bool recursive, int verb) {

  Interface &m_field = cOre;
  MeshsetsManager *meshsets_manager_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.getInterface(meshsets_manager_ptr);
  CubitMeshSet_multiIndex::index<
      Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator miit =
      meshsets_manager_ptr->getMeshsetsMultindex()
          .get<Composite_Cubit_msId_And_MeshSetType_mi_tag>()
          .find(boost::make_tuple(msId, cubit_bc_type.to_ulong()));
  if (miit != meshsets_manager_ptr->getMeshsetsMultindex()
                  .get<Composite_Cubit_msId_And_MeshSetType_mi_tag>()
                  .end()) {
    CHKERR getSides(miit->meshset, mesh_bit_level, recursive, verb);
  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "msId is not there");
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PrismInterface::getSides(const EntityHandle sideset,
                                        const BitRefLevel mesh_bit_level,
                                        const bool recursive, int verb) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  Skinner skin(&moab);

  auto save_range = [&moab](const std::string name, const Range r) {
    MoFEMFunctionBegin;
    EntityHandle out_meshset;
    CHKERR moab.create_meshset(MESHSET_SET, out_meshset);
    CHKERR moab.add_entities(out_meshset, r);
    CHKERR moab.write_file(name.c_str(), "VTK", "", &out_meshset, 1);
    CHKERR moab.delete_entities(&out_meshset, 1);
    MoFEMFunctionReturn(0);
  };

  auto get_adj = [&moab](const Range r, const int dim) {
    Range a;
    if (dim)
      CHKERR moab.get_adjacencies(r, dim, false, a, moab::Interface::UNION);
    else
      CHKERR moab.get_connectivity(r, a, true);
    return a;
  };

  auto get_skin = [&skin](const auto r) {
    Range s;
    CHKERR skin.find_skin(0, r, false, s);
    return s;
  };

  MoFEMFunctionBegin;

  Range triangles;
  CHKERR moab.get_entities_by_type(sideset, MBTRI, triangles, recursive);

  Range mesh_level_ents3d, mesh_level_ents3d_tris;
  Range mesh_level_tris;
  Range mesh_level_nodes;
  Range mesh_level_prisms;

  if (mesh_bit_level.any()) {
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBTET, mesh_level_ents3d);
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBTRI, mesh_level_tris);
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBVERTEX, mesh_level_nodes);
    mesh_level_ents3d_tris = get_adj(mesh_level_ents3d, 2);
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBPRISM, mesh_level_prisms);
    mesh_level_ents3d.merge(mesh_level_prisms);
  }

  // get interface triangles from side set
  if (mesh_bit_level.any())
    triangles = intersect(triangles, mesh_level_ents3d_tris);
  if (triangles.empty())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Range of triangles set to split is emtpy. Nothing to split.");

  MOFEM_LOG_C("PRISM_INTERFACE", Sev::verbose, "Nb. of triangles in set %u",
              triangles.size());

  auto get_skin_edges_boundary = [&]() {
    // get nodes, edges and 3d ents (i.e. tets and prisms)
    auto ents3d_with_prisms = get_adj(get_adj(triangles, 0), 3);
    if (mesh_bit_level.any())
      ents3d_with_prisms = intersect(ents3d_with_prisms, mesh_level_ents3d);
    auto ents3d = ents3d_with_prisms.subset_by_type(
        MBTET); // take only tets, add prism later

    // note: that skin faces edges do not contain internal boundary
    // note: that prisms are not included in ents3d, so if ents3d have border
    // with other inteface is like external boundary skin edges boundary are
    // internal edge <- skin_faces_edges contains edges which are on the body
    // boundary <- that is the trick
    auto skin_edges_boundary =
        subtract(get_skin(triangles),
                 get_adj(get_skin(ents3d),
                         1)); // from skin edges subtract edges from skin
                              // faces of 3d ents (only internal edges)

    skin_edges_boundary =
        subtract(skin_edges_boundary,
                 get_adj(ents3d_with_prisms.subset_by_type(MBPRISM),
                         1)); // from skin edges subtract edges from prism
                              // edges, that create internal boundary

    return skin_edges_boundary;
  };

  auto skin_edges_boundary = get_skin_edges_boundary();
  auto skin_nodes_boundary = get_adj(skin_edges_boundary, 0);

  auto get_edges_without_boundary = [&]() {
    // Get front edges
    return subtract(get_adj(triangles, 1), skin_edges_boundary);
  };

  auto get_nodes_without_front = [&]() {
    // use nodes on body boundary and interface (without internal boundary) to
    // find adjacent tets
    return subtract(get_adj(triangles, 0),
                    skin_nodes_boundary); // nodes_without_front adjacent to
                                          // all split face edges except
                                          // those on internal edge
  };

  auto get_ents3d_with_prisms = [&](auto edges_without_boundary,
                                    auto nodes_without_front) {
    // ents3 that are adjacent to front nodes on split faces but not those which
    // are on the front nodes on internal edge
    Range ents3d_with_prisms =
        get_adj(unite(edges_without_boundary, nodes_without_front), 3);

    auto find_triangles_on_front_and_adjacent_tet = [&]() {
      MoFEMFunctionBegin;
      // get all triangles adjacent to front
      auto skin_nodes_boundary_tris = get_adj(skin_nodes_boundary, 2);

      // get nodes of triangles adjacent to front nodes
      // get hanging nodes, i.e. nodes which are not on the front but adjacent
      // to triangles adjacent to crack front
      auto skin_nodes_boundary_tris_nodes =
          subtract(get_adj(skin_nodes_boundary_tris, 2), skin_nodes_boundary);

      // get triangles adjacent to hanging nodes
      auto skin_nodes_boundary_tris_nodes_tris =
          get_adj(skin_nodes_boundary_tris_nodes, 2);

      // triangles which have tree nodes on front boundary
      skin_nodes_boundary_tris =
          intersect(triangles, subtract(skin_nodes_boundary_tris,
                                        skin_nodes_boundary_tris_nodes_tris));
      if (!skin_nodes_boundary_tris.empty()) {
        // Get internal edges of triangle which has three nodes on boundary
        auto skin_nodes_boundary_tris_edges =
            get_adj(skin_nodes_boundary_tris, 1);

        skin_nodes_boundary_tris_edges =
            subtract(skin_nodes_boundary_tris_edges, skin_edges_boundary);
        // Get 3d elements adjacent to internal edge which has three nodes on
        // boundary
        ents3d_with_prisms.merge(get_adj(skin_nodes_boundary_tris_edges, 3));
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR find_triangles_on_front_and_adjacent_tet();

    // prism and tets on both side of interface
    if (mesh_bit_level.any())
      ents3d_with_prisms = intersect(ents3d_with_prisms, mesh_level_ents3d);

    if (verb >= NOISY) {
      CHKERR save_range("skin_edges_boundary.vtk", skin_edges_boundary);
      CHKERR save_range("nodes_without_front.vtk", nodes_without_front);
    }

    return ents3d_with_prisms;
  };

  auto nodes_without_front = get_nodes_without_front();
  auto ents3d_with_prisms =
      get_ents3d_with_prisms(get_edges_without_boundary(), nodes_without_front);
  auto ents3d = ents3d_with_prisms.subset_by_type(MBTET);

  MOFEM_LOG_C("PRISM_INTERFACE", Sev::noisy,
              "Number of adjacents 3d entities to front nodes %u",
              ents3d.size());

  if (verb >= NOISY) {
    CHKERR save_range("triangles.vtk", triangles);
    CHKERR save_range("ents3d.vtk", ents3d);
    CHKERR save_range("skin_nodes_boundary.vtk", skin_nodes_boundary);
  }

  auto find_tetrahedrons_on_the_side = [&]() {
    auto seed = intersect(get_adj(triangles, 3), ents3d);
    Range side_ents3d;
    if (!seed.empty())
      side_ents3d.insert(seed[0]);
    unsigned int nb_side_ents3d = side_ents3d.size();
    Range side_ents3d_tris_on_surface;

    // get all tets adjacent to crack surface, but only on one side of it
    do {

      Range adj_ents3d;
      do {

        Range adj_tris;
        nb_side_ents3d = side_ents3d.size();
        MOFEM_LOG_C("PRISM_INTERFACE", Sev::noisy,
                    "Number of entities on side %u", nb_side_ents3d);

        // get faces
        // subtrace from faces interface
        adj_tris = get_skin(side_ents3d.subset_by_type(MBTET));
        adj_tris = subtract(adj_tris, triangles);
        if (mesh_bit_level.any())
          adj_tris = intersect(adj_tris, mesh_level_tris);

        // get tets adjacent to faces
        CHKERR moab.get_adjacencies(adj_tris, 3, false, adj_ents3d,
                                    moab::Interface::UNION);
        // intersect tets with tets adjacent to inetface
        adj_ents3d = intersect(adj_ents3d, ents3d_with_prisms);

        // add tets to side
        side_ents3d.insert(adj_ents3d.begin(), adj_ents3d.end());
        if (verb >= VERY_NOISY) {
          CHKERR save_range(
              "side_ents3d_" +
                  boost::lexical_cast<std::string>(nb_side_ents3d) + ".vtk",
              side_ents3d);
        }

      } while (nb_side_ents3d != side_ents3d.size());
      Range side_ents3d_tris;
      CHKERR moab.get_adjacencies(side_ents3d, 2, false, side_ents3d_tris,
                                  moab::Interface::UNION);
      side_ents3d_tris_on_surface = intersect(side_ents3d_tris, triangles);

      if (verb >= VERY_NOISY) {
        Range left_triangles = subtract(triangles, side_ents3d_tris_on_surface);
        if (!left_triangles.empty()) {
          CHKERR save_range(
              "left_triangles_" +
                  boost::lexical_cast<std::string>(nb_side_ents3d) + ".vtk",
              left_triangles);
        }
      }

      // This is a case when separate sub-domains are split, so we need
      // additional tetrahedron for seed process
      if (side_ents3d_tris_on_surface.size() != triangles.size()) {
        auto left_triangles = subtract(triangles, side_ents3d_tris_on_surface);
        Range tets;
        CHKERR moab.get_adjacencies(&*left_triangles.begin(), 1, 3, false,
                                    tets);
        tets = intersect(tets, ents3d_with_prisms);

        if (tets.empty()) {
          CHKERR save_range("error.vtk", left_triangles);
          THROW_MESSAGE(
              "Not all faces on surface going to be split, see error.vtk for "
              "problematic triangle. "
              "It could be a case where triangle on front (part boundary of "
              "surface in interior) "
              "has three nodes front.");
        }
        side_ents3d.insert(*tets.begin());
      }

    } while (side_ents3d_tris_on_surface.size() != triangles.size());

    return side_ents3d;
  };

  auto side_ents3d = find_tetrahedrons_on_the_side();
  if (ents3d_with_prisms.size() == side_ents3d.size())
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "All tetrahedrons are on one side of split surface and that is "
            "wrong. Algorithm can not distinguish (find) sides of interface.");

  // other side ents
  auto other_side = subtract(ents3d_with_prisms, side_ents3d);
  // side nodes
  auto side_nodes = get_adj(side_ents3d.subset_by_type(MBTET), 0);
  // nodes on crack surface without front
  nodes_without_front = intersect(nodes_without_front, side_nodes);
  auto side_edges = get_adj(side_ents3d.subset_by_type(MBTET), 1);
  skin_edges_boundary = intersect(skin_edges_boundary, side_edges);

  // make child meshsets
  std::vector<EntityHandle> children;
  CHKERR moab.get_child_meshsets(sideset, children);
  if (children.empty()) {
    children.resize(3);
    CHKERR moab.create_meshset(MESHSET_SET, children[0]);
    CHKERR moab.create_meshset(MESHSET_SET, children[1]);
    CHKERR moab.create_meshset(MESHSET_SET, children[2]);
    CHKERR moab.add_child_meshset(sideset, children[0]);
    CHKERR moab.add_child_meshset(sideset, children[1]);
    CHKERR moab.add_child_meshset(sideset, children[2]);
  } else {
    if (children.size() != 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "this meshset should have 3 children meshsets");
    }
    children.resize(3);
    CHKERR moab.clear_meshset(&children[0], 3);
  }

  EntityHandle &child_side = children[0];
  EntityHandle &child_other_side = children[1];
  EntityHandle &child_nodes_and_skin_edges = children[2];
  CHKERR moab.add_entities(child_side, side_ents3d);
  CHKERR moab.add_entities(child_other_side, other_side);
  CHKERR moab.add_entities(child_nodes_and_skin_edges, nodes_without_front);
  CHKERR moab.add_entities(child_nodes_and_skin_edges, skin_edges_boundary);

  MOFEM_LOG_C("PRISM_INTERFACE", Sev::verbose,
              "Nb. of side 3d elements in set %u", side_ents3d.size());
  MOFEM_LOG_C("PRISM_INTERFACE", Sev::verbose,
              "Nb. of other side 3d elements in set %u", other_side.size());
  MOFEM_LOG_C("PRISM_INTERFACE", Sev::verbose, "Nb. of boundary edges %u",
              skin_edges_boundary.size());

  if (verb >= NOISY) {
    CHKERR moab.write_file("side.vtk", "VTK", "", &children[0], 1);
    CHKERR moab.write_file("other_side.vtk", "VTK", "", &children[1], 1);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PrismInterface::findFacesWithThreeNodesOnInternalSurfaceSkin(
    const EntityHandle sideset, const BitRefLevel mesh_bit_level,
    const bool recursive, Range &faces_with_three_nodes_on_front, int verb) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  Range mesh_level_ents3d;
  Range mesh_level_edges, mesh_level_tris;
  if (mesh_bit_level.any()) {
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBTET, mesh_level_ents3d);

    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBTRI, mesh_level_tris);

    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBEDGE, mesh_level_edges);
  }

  Skinner skin(&moab);
  // get interface triangles from side set
  Range triangles;
  CHKERR moab.get_entities_by_type(sideset, MBTRI, triangles, recursive);

  if (mesh_bit_level.any()) {
    triangles = intersect(triangles, mesh_level_tris);
  }

  MOFEM_LOG_C("PRISM_INTERFACE", Sev::verbose, "Nb. of triangles in set %u",
              triangles.size());

  // get nodes, edges and 3d ents (i.e. tets and prisms)
  Range nodes; // nodes from triangles
  CHKERR moab.get_connectivity(triangles, nodes, true);

  Range ents3d; // 3d ents from nodes
  CHKERR moab.get_adjacencies(nodes, 3, false, ents3d, moab::Interface::UNION);

  if (mesh_bit_level.any()) {
    ents3d = intersect(ents3d, mesh_level_ents3d);
  }
  // take skin faces
  Range skin_faces; // skin faces from 3d ents
  CHKERR skin.find_skin(0, ents3d.subset_by_type(MBTET), false, skin_faces);

  // take skin edges (boundary of surface if there is any)
  Range skin_edges_boundary; // skin edges from triangles
  CHKERR skin.find_skin(0, triangles, false, skin_edges_boundary);

  // take all edges on skin faces (i.e. skin surface)
  Range skin_faces_edges; // edges from skin faces of 3d ents
  CHKERR moab.get_adjacencies(skin_faces, 1, false, skin_faces_edges,
                              moab::Interface::UNION);

  if (mesh_bit_level.any()) {
    skin_faces_edges = intersect(skin_faces_edges, mesh_level_edges);
  }

  // note: that skin faces edges do not contain internal boundary
  // note: that prisms are not included in ents3d, so if ents3d have border with
  // other inteface is like external boundary skin edges boundary are internal
  // edge <- skin_faces_edges contains edges which are on the body boundary <-
  // that is the trick
  skin_edges_boundary =
      subtract(skin_edges_boundary,
               skin_faces_edges); // from skin edges subtract edges from skin
                                  // faces of 3d ents (only internal edges)

  if (verb >= VERY_VERBOSE) {
    EntityHandle out_meshset;
    CHKERR moab.create_meshset(MESHSET_SET, out_meshset);
    CHKERR moab.add_entities(out_meshset, triangles);
    CHKERR moab.write_file("triangles.vtk", "VTK", "", &out_meshset, 1);
    CHKERR moab.delete_entities(&out_meshset, 1);
    CHKERR moab.create_meshset(MESHSET_SET, out_meshset);
    CHKERR moab.add_entities(out_meshset, ents3d);
    CHKERR moab.write_file("ents3d.vtk", "VTK", "", &out_meshset, 1);
    CHKERR moab.delete_entities(&out_meshset, 1);
    CHKERR moab.create_meshset(MESHSET_SET, out_meshset);
    CHKERR moab.add_entities(out_meshset, skin_edges_boundary);
    CHKERR moab.write_file("skin_edges_boundary.vtk", "VTK", "", &out_meshset,
                           1);
    CHKERR moab.delete_entities(&out_meshset, 1);
  }

  // Get nodes on boundary edge
  Range skin_nodes_boundary;
  CHKERR moab.get_connectivity(skin_edges_boundary, skin_nodes_boundary, true);

  // Remove node which are boundary with other existing interface
  Range prisms_nodes;
  CHKERR
  moab.get_connectivity(ents3d.subset_by_type(MBPRISM), prisms_nodes, true);

  skin_nodes_boundary = subtract(skin_nodes_boundary, prisms_nodes);

  // use nodes on body boundary and interface (without internal boundary) to
  // find adjacent tets
  Range nodes_without_front = subtract(
      nodes, skin_nodes_boundary); // nodes_without_front adjacent to all split
                                   // face edges except those on internal edge

  Range skin_nodes_boundary_tris;
  CHKERR moab.get_adjacencies(skin_nodes_boundary, 2, false,
                              skin_nodes_boundary_tris, moab::Interface::UNION);
  Range skin_nodes_boundary_tris_nodes;
  CHKERR moab.get_connectivity(skin_nodes_boundary_tris,
                               skin_nodes_boundary_tris_nodes, true);

  skin_nodes_boundary_tris_nodes =
      subtract(skin_nodes_boundary_tris_nodes, skin_nodes_boundary);
  Range skin_nodes_boundary_tris_nodes_tris;
  CHKERR moab.get_adjacencies(skin_nodes_boundary_tris_nodes, 2, false,
                              skin_nodes_boundary_tris_nodes_tris,
                              moab::Interface::UNION);

  // Triangle which has tree nodes on front boundary
  skin_nodes_boundary_tris =
      intersect(triangles, subtract(skin_nodes_boundary_tris,
                                    skin_nodes_boundary_tris_nodes_tris));
  faces_with_three_nodes_on_front.swap(skin_nodes_boundary_tris);

  MoFEMFunctionReturn(0);
} // namespace MoFEM

MoFEMErrorCode PrismInterface::splitSides(const EntityHandle meshset,
                                          const BitRefLevel &bit,
                                          const int msId,
                                          const CubitBCType cubit_bc_type,
                                          const bool add_interface_entities,
                                          const bool recursive, int verb) {

  Interface &m_field = cOre;
  MeshsetsManager *meshsets_manager_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.getInterface(meshsets_manager_ptr);
  CubitMeshSet_multiIndex::index<
      Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator miit =
      meshsets_manager_ptr->getMeshsetsMultindex()
          .get<Composite_Cubit_msId_And_MeshSetType_mi_tag>()
          .find(boost::make_tuple(msId, cubit_bc_type.to_ulong()));
  if (miit != meshsets_manager_ptr->getMeshsetsMultindex()
                  .get<Composite_Cubit_msId_And_MeshSetType_mi_tag>()
                  .end()) {
    CHKERR splitSides(meshset, bit, miit->meshset, add_interface_entities,
                      recursive, verb);
  } else {
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY, "msId is not there");
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PrismInterface::splitSides(const EntityHandle meshset,
                                          const BitRefLevel &bit,
                                          const EntityHandle sideset,
                                          const bool add_interface_entities,
                                          const bool recursive, int verb) {

  MoFEMFunctionBegin;
  CHKERR splitSides(meshset, bit, BitRefLevel(), BitRefLevel(), sideset,
                    add_interface_entities, recursive, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PrismInterface::splitSides(
    const EntityHandle meshset, const BitRefLevel &bit,
    const BitRefLevel &inhered_from_bit_level,
    const BitRefLevel &inhered_from_bit_level_mask, const EntityHandle sideset,
    const bool add_interface_entities, const bool recursive, int verb) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  auto refined_ents_ptr = m_field.get_ref_ents();
  MoFEMFunctionBegin;

  std::vector<EntityHandle> children;
  // get children meshsets
  CHKERR moab.get_child_meshsets(sideset, children);
  if (children.size() != 3)
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "should be 3 child meshsets, each of them contains tets on two "
            "sides of interface");

  // 3d ents on "father" side
  Range side_ents3d;
  CHKERR moab.get_entities_by_handle(children[0], side_ents3d, false);
  // 3d ents on "mother" side
  Range other_ents3d;
  CHKERR moab.get_entities_by_handle(children[1], other_ents3d, false);
  // faces of interface
  Range triangles;
  CHKERR moab.get_entities_by_type(sideset, MBTRI, triangles, recursive);
  Range side_ents3d_tris;
  CHKERR moab.get_adjacencies(side_ents3d, 2, true, side_ents3d_tris,
                              moab::Interface::UNION);
  triangles = intersect(triangles, side_ents3d_tris);
  // nodes on interface but not on crack front (those should not be splitted)
  Range nodes;
  CHKERR moab.get_entities_by_type(children[2], MBVERTEX, nodes, false);
  Range meshset_3d_ents, meshset_2d_ents;
  CHKERR moab.get_entities_by_dimension(meshset, 3, meshset_3d_ents, true);
  Range meshset_tets = meshset_3d_ents.subset_by_type(MBTET);
  CHKERR moab.get_adjacencies(meshset_tets, 2, false, meshset_2d_ents,
                              moab::Interface::UNION);
  side_ents3d = intersect(meshset_3d_ents, side_ents3d);
  other_ents3d = intersect(meshset_3d_ents, other_ents3d);
  triangles = intersect(meshset_2d_ents, triangles);

  MOFEM_LOG_C("PRISM_INTERFACE", Sev::verbose, "Split sides triangles %u",
              triangles.size());
  MOFEM_LOG_C("PRISM_INTERFACE", Sev::verbose, "Split sides 3d entities %u",
              side_ents3d.size());
  MOFEM_LOG_C("PRISM_INTERFACE", Sev::verbose, "split sides nodes %u",
              nodes.size());

  struct PartentAndChild {
    EntityHandle parent;
    EntityHandle child;
  };

  typedef multi_index_container<
      PartentAndChild,
      indexed_by<

          hashed_unique<
              member<PartentAndChild, EntityHandle, &PartentAndChild::parent>>,

          hashed_unique<
              member<PartentAndChild, EntityHandle, &PartentAndChild::child>>

          >>
      ParentChildMI;

  ParentChildMI parent_child;

  typedef std::map<EntityHandle, /*node on "mother" side*/
                   EntityHandle  /*node on "father" side*/
                   >
      MapNodes;
  MapNodes map_nodes, reverse_map_nodes;

  // Map nodes on sides, set parent node and set bit ref level
  {
    struct CreateSideNodes {
      MoFEM::Core &cOre;
      MoFEM::Interface &m_field;
      std::vector<EntityHandle> splitNodes;
      std::vector<double> splitCoords[3];

      RefEntity_multiIndex_view_by_hashed_parent_entity refParentEntsView;

      CreateSideNodes(MoFEM::Core &core, int split_size = 0)
          : cOre(core), m_field(core) {
        splitNodes.reserve(split_size);
        for (auto dd : {0, 1, 2})
          splitCoords[dd].reserve(split_size);
      }
      MoFEMErrorCode operator()(const double coords[], const EntityHandle n) {
        MoFEMFunctionBegin;
        splitNodes.emplace_back(n);
        for (auto dd : {0, 1, 2})
          splitCoords[dd].emplace_back(coords[dd]);
        MoFEMFunctionReturn(0);
      }

      MoFEMErrorCode operator()(const BitRefLevel &bit, MapNodes &map_nodes,
                                MapNodes &reverse_map_nodes) {
        ReadUtilIface *iface;
        MoFEMFunctionBegin;
        int num_nodes = splitNodes.size();
        std::vector<double *> arrays_coord;
        EntityHandle startv;
        CHKERR m_field.get_moab().query_interface(iface);
        CHKERR iface->get_node_coords(3, num_nodes, 0, startv, arrays_coord);
        Range verts(startv, startv + num_nodes - 1);
        for (int dd = 0; dd != 3; ++dd)
          std::copy(splitCoords[dd].begin(), splitCoords[dd].end(),
                    arrays_coord[dd]);
        for (int nn = 0; nn != num_nodes; ++nn) {
          map_nodes[splitNodes[nn]] = verts[nn];
          reverse_map_nodes[verts[nn]] = splitNodes[nn];
        }
        CHKERR m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(),
                                               verts, &*splitNodes.begin());
        CHKERR m_field.getInterface<BitRefManager>()->setEntitiesBitRefLevel(
            verts, bit, QUIET);
        MoFEMFunctionReturn(0);
      }
    };
    CreateSideNodes create_side_nodes(cOre, nodes.size());

    RefEntity_multiIndex_view_by_hashed_parent_entity ref_parent_ents_view;
    struct CreateParentEntView {
      MoFEMErrorCode
      operator()(const BitRefLevel &bit, const BitRefLevel &mask,
                 const RefEntity_multiIndex *refined_ents_ptr,
                 RefEntity_multiIndex_view_by_hashed_parent_entity
                     &ref_parent_ents_view) const {
        MoFEMFunctionBegin;
        auto &ref_ents =
            refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>();
        // view by parent type (VERTEX)
        auto miit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX, MBVERTEX));
        auto hi_miit =
            ref_ents.upper_bound(boost::make_tuple(MBVERTEX, MBVERTEX));
        for (; miit != hi_miit; miit++) {
          const auto &ent_bit = (*miit)->getBitRefLevel();
          if ((ent_bit & bit).any() && (ent_bit & mask) == ent_bit) {
            auto p_ref_ent_view = ref_parent_ents_view.insert(*miit);
            if (!p_ref_ent_view.second)
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                      "non unique insertion");
          }
        }
        MoFEMFunctionReturn(0);
      }
    };
    if (inhered_from_bit_level.any() && inhered_from_bit_level_mask.any())
      CHKERR CreateParentEntView()(inhered_from_bit_level,
                                   inhered_from_bit_level_mask,
                                   refined_ents_ptr, ref_parent_ents_view);

    // add new nodes on interface and create map
    Range add_bit_nodes;
    for (auto pnit = nodes.pair_begin(); pnit != nodes.pair_end(); ++pnit) {
      auto lo = refined_ents_ptr->lower_bound(pnit->first);
      auto hi = refined_ents_ptr->upper_bound(pnit->second);
      if (std::distance(lo, hi) != (pnit->second - pnit->first + 1))
        SETERRQ(
            PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Can not find some nodes in database that are split on interface");
      Range nodes_in_range;
      insertOrdered(nodes_in_range, RefEntExtractor(), lo, hi);
      std::vector<double> coords_range(nodes_in_range.size() * 3);
      CHKERR moab.get_coords(nodes_in_range, &*coords_range.begin());
      int pos = 0;
      for (; lo != hi; ++lo, pos += 3) {
        const EntityHandle node = (*lo)->getEnt();
        EntityHandle child_entity = 0;
        auto child_it = ref_parent_ents_view.find(node);
        if (child_it != ref_parent_ents_view.end())
          child_entity = (*child_it)->getEnt();
        if (child_entity == 0) {
          CHKERR create_side_nodes(&coords_range[pos], node);
        } else {
          map_nodes[node] = child_entity;
          add_bit_nodes.insert(child_entity);
        }
      }
    }
    add_bit_nodes.merge(nodes);
    CHKERR m_field.getInterface<BitRefManager>()->addBitRefLevel(add_bit_nodes,
                                                                 bit);
    CHKERR create_side_nodes(bit, map_nodes, reverse_map_nodes);
  }

  // crete meshset for new mesh bit level
  EntityHandle meshset_for_bit_level;
  CHKERR moab.create_meshset(MESHSET_SET, meshset_for_bit_level);

  // subtract those elements which will be refined, i.e. disconnected from other
  // side elements, and connected to new prisms, if they are created
  meshset_3d_ents = subtract(meshset_3d_ents, side_ents3d);
  CHKERR moab.add_entities(meshset_for_bit_level, meshset_3d_ents);

  // create new 3d ents on "father" side
  Range new_3d_ents;
  for (Range::iterator eit3d = side_ents3d.begin(); eit3d != side_ents3d.end();
       eit3d++) {
    auto miit_ref_ent = refined_ents_ptr->find(*eit3d);
    if (miit_ref_ent == refined_ents_ptr->end())
      SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
              "Tetrahedron not in database");

    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(*eit3d, conn, num_nodes, true);
    EntityHandle new_conn[num_nodes];
    int nb_new_conn = 0;
    for (int ii = 0; ii < num_nodes; ii++) {
      std::map<EntityHandle, EntityHandle>::iterator mit =
          map_nodes.find(conn[ii]);
      if (mit != map_nodes.end()) {
        new_conn[ii] = mit->second;
        nb_new_conn++;
        if (verb >= VERY_NOISY) {
          MOFEM_LOG_C("PRISM_INTERFACE", Sev::noisy, "nodes %u -> %d",
                      conn[ii], new_conn[ii]);
          MOFEM_LOG_C("PRISM_INTERFACE", Sev::noisy, "nodes %u -> %d",
                      conn[ii], new_conn[ii]);
        }
      } else {
        new_conn[ii] = conn[ii];
      }
    }
    if (nb_new_conn == 0) {
      // Add this tet to bit ref level
      CHKERR moab.add_entities(meshset_for_bit_level, &*eit3d, 1);
      continue;
    }

    // here is created new or prism is on interface
    EntityHandle existing_ent = 0;
    /* check if tet element with new connectivity is in database*/
    auto child_iit =
        refined_ents_ptr->get<Ent_Ent_mi_tag>().lower_bound(*eit3d);
    auto hi_child_iit =
        refined_ents_ptr->get<Ent_Ent_mi_tag>().upper_bound(*eit3d);

    // Check if child entity has the same connectivity
    for (; child_iit != hi_child_iit; child_iit++) {
      const EntityHandle *conn_ref_tet;
      CHKERR moab.get_connectivity(child_iit->get()->getEnt(), conn_ref_tet,
                                   num_nodes, true);
      int nn = 0;
      for (; nn < num_nodes; nn++) {
        if (conn_ref_tet[nn] != new_conn[nn]) {
          break;
        }
      }
      if (nn == num_nodes) {
        if (existing_ent != 0)
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "Should be only one child entity with the same connectivity");
        existing_ent = child_iit->get()->getEnt();
      }
    }

    switch (moab.type_from_handle(*eit3d)) {
    case MBTET: {

      RefEntity_multiIndex::iterator child_it;
      EntityHandle tet;
      if (existing_ent == 0) {
        Range new_conn_tet;
        CHKERR moab.get_adjacencies(new_conn, 4, 3, false, new_conn_tet);
        if (new_conn_tet.empty()) {
          CHKERR moab.create_element(MBTET, new_conn, 4, tet);
          CHKERR moab.tag_set_data(cOre.get_th_RefParentHandle(), &tet, 1,
                                   &*eit3d);

        } else {

          auto new_rit = refined_ents_ptr->get<Ent_mi_tag>().equal_range(
              *new_conn_tet.begin());

          size_t nb_elems = std::distance(new_rit.first, new_rit.second);
          if (nb_elems != 1)
            SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                     "Can't find entity in database, size is %d", nb_elems);
          tet = *new_conn_tet.begin();
        }
      } else {
        tet = existing_ent;
      }

      CHKERR moab.add_entities(meshset_for_bit_level, &tet, 1);
      new_3d_ents.insert(tet);

    } break;
    case MBPRISM: {
      EntityHandle prism;
      if (existing_ent == 0) {
        Range new_conn_prism;
        CHKERR moab.get_adjacencies(new_conn, 6, 3, false, new_conn_prism);
        if (new_conn_prism.empty()) {
          CHKERR moab.create_element(MBPRISM, new_conn, 6, prism);
          CHKERR moab.tag_set_data(cOre.get_th_RefParentHandle(), &prism, 1,
                                   &*eit3d);
        } else {
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "It is prism with such connectivity, that case has to be "
                  "handled but this is not implemented");
        }
      } else {
        prism = existing_ent;
      }
      CHKERR moab.add_entities(meshset_for_bit_level, &prism, 1);
      new_3d_ents.insert(prism);
    } break;
    default:
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "Not implemented element");
    }
  }

  struct SetParent {

    SetParent(MoFEM::Core &core) : cOre(core), mField(core) {}

    MoFEMErrorCode operator()(const EntityHandle ent, const EntityHandle parent,
                              const RefEntity_multiIndex *ref_ents_ptr) {
      MoFEMFunctionBegin;
      if (ent != parent) {
        auto it = ref_ents_ptr->find(ent);
        if (it != ref_ents_ptr->end()) {
          parentsToChange[ent] = parent;
        } else {
          CHKERR mField.get_moab().tag_set_data(cOre.get_th_RefParentHandle(),
                                                &ent, 1, &parent);
        }
      }
      MoFEMFunctionReturn(0);
    }

    MoFEMErrorCode override_parents(const RefEntity_multiIndex *ref_ents_ptr) {
      MoFEMFunctionBegin;
      for (auto &m : parentsToChange) {
        auto it = ref_ents_ptr->find(m.first);
        if (it != ref_ents_ptr->end()) {

          bool success = const_cast<RefEntity_multiIndex *>(ref_ents_ptr)
                             ->modify(it, RefEntity_change_parent(m.second));
          if (!success)
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "Impossible to set parent");
        } else
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Entity not in database");
      }
      MoFEMFunctionReturn(0);
    }

  private:
    MoFEM::Core &cOre;
    MoFEM::Interface &mField;
    map<EntityHandle, EntityHandle> parentsToChange;
  };

  SetParent set_parent(cOre);

  auto get_adj_ents = [&](const bool create) {
    Range adj;
    for (auto d : {1, 2}) {
      // create new entities by adjacencies form new tets
      CHKERR moab.get_adjacencies(new_3d_ents.subset_by_type(MBTET), d, create,
                                  adj, moab::Interface::UNION);
    }
    return adj;
  };

  // Create entities
  get_adj_ents(true);

  auto get_conn = [&](const auto e) {
    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(e, conn, num_nodes, true);
    return std::make_pair(conn, num_nodes);
  };

  auto get_new_conn = [&](auto conn) {
    std::array<EntityHandle, 8> new_conn;
    int nb_new_conn = 0;
    for (int ii = 0; ii != conn.second; ++ii) {
      auto mit = map_nodes.find(conn.first[ii]);
      if (mit != map_nodes.end()) {
        new_conn[ii] = mit->second;
        nb_new_conn++;
        if (verb >= VERY_NOISY)
          PetscPrintf(m_field.get_comm(), "nodes %u -> %d", conn.first[ii],
                      new_conn[ii]);

      } else {
        new_conn[ii] = conn.first[ii];
      }
    }
    return std::make_pair(new_conn, nb_new_conn);
  };

  auto get_reverse_conn = [&](auto conn) {
    std::array<EntityHandle, 8> rev_conn;
    int nb_new_conn = 0;
    for (int ii = 0; ii != conn.second; ++ii) {
      auto mit = reverse_map_nodes.find(conn.first[ii]);
      if (mit != reverse_map_nodes.end()) {
        rev_conn[ii] = mit->second;
        nb_new_conn++;
        if (verb >= VERY_NOISY)
          PetscPrintf(m_field.get_comm(), "nodes %u -> %d", conn.first[ii],
                      rev_conn[ii]);

      } else {
        rev_conn[ii] = conn.first[ii];
      }
    }
    return std::make_pair(rev_conn, nb_new_conn);
  };

  auto get_new_ent = [&](auto new_conn, auto nb_nodes, int dim) {
    Range new_ent;
    CHKERR moab.get_adjacencies(&(new_conn.first[0]), nb_nodes, dim, false,
                                new_ent);
    if (new_ent.size() != 1)
      THROW_MESSAGE("this entity should be in moab database");
    return new_ent.front();
  };

  auto create_prisms = [&]() {
    MoFEMFunctionBegin;

    // Tags for setting side
    Tag th_interface_side;
    const int def_side[] = {0};
    CHKERR moab.tag_get_handle("INTERFACE_SIDE", 1, MB_TYPE_INTEGER,
                               th_interface_side, MB_TAG_CREAT | MB_TAG_SPARSE,
                               def_side);

    for (auto p = triangles.pair_begin(); p != triangles.pair_end(); ++p) {
      auto f = p->first;
      auto s = p->second;

      auto lo = refined_ents_ptr->lower_bound(f);
      auto hi = refined_ents_ptr->upper_bound(s);
      if (std::distance(lo, hi) != (s - f + 1))
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "Some triangles are not in database");

      for (; f <= s; ++f) {

        auto conn = get_conn(f);
        auto new_conn = get_new_conn(conn);

        if (new_conn.second) {

          auto set_side_tag = [&](auto new_triangle) {
            int new_side = 1;
            CHKERR moab.tag_set_data(th_interface_side, &new_triangle, 1,
                                     &new_side);
          };

          auto get_ent3d = [&](auto e) {
            Range ents_3d;
            CHKERR moab.get_adjacencies(&e, 1, 3, false, ents_3d);
            ents_3d = intersect(ents_3d, side_ents3d);

            switch (ents_3d.size()) {
            case 0:
              THROW_MESSAGE(
                  "Did not find adjacent tets on one side of the interface; if "
                  "this error appears for a contact interface, try creating "
                  "separate blocksets for each contact surface");
            case 2:
              THROW_MESSAGE(
                  "Found both adjacent tets on one side of the interface, if "
                  "this error appears for a contact interface, try creating "
                  "separate blocksets for each contact surface");
            default:
              break;
            }

            return ents_3d.front();
          };

          auto get_sense = [&](auto e, auto ent3d) {
            int sense, side, offset;
            CHKERR moab.side_number(ent3d, e, side, sense, offset);
            if (sense != 1 && sense != -1) {
              THROW_MESSAGE(
                  "Undefined sense of a triangle; if this error appears for a "
                  "contact interface, try creating separate blocksets for each "
                  "contact surface");
            }
            return sense;
          };

          auto new_triangle = get_new_ent(new_conn, 3, 2);
          set_side_tag(new_triangle);

          if (add_interface_entities) {

            if (inhered_from_bit_level.any())
              SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                      "not implemented for inhered_from_bit_level");

            // set prism connectivity
            EntityHandle prism_conn[6] = {
                conn.first[0],     conn.first[1],     conn.first[2],

                new_conn.first[0], new_conn.first[1], new_conn.first[2]};
            if (get_sense(f, get_ent3d(f)) == -1) {
              // swap nodes in triangles for correct prism creation
              std::swap(prism_conn[1], prism_conn[2]);
              std::swap(prism_conn[4], prism_conn[5]);
            }

            EntityHandle prism;
            CHKERR moab.create_element(MBPRISM, prism_conn, 6, prism);
            CHKERR moab.add_entities(meshset_for_bit_level, &prism, 1);
          }
        }
      }
    }

    MoFEMFunctionReturn(0);
  };

  auto set_parnets = [&](auto side_adj_faces_and_edges) {
    MoFEMFunctionBegin;

    for (auto p = side_adj_faces_and_edges.pair_begin();
         p != side_adj_faces_and_edges.pair_end(); ++p) {
      auto f = p->first;
      auto s = p->second;

      for (; f <= s; ++f) {
        auto conn = get_conn(f);
        auto rev_conn = get_reverse_conn(conn);
        if (rev_conn.second) {
          auto rev_ent = get_new_ent(rev_conn, conn.second, conn.second - 1);
          CHKERR set_parent(f, rev_ent, refined_ents_ptr);
        }
      }
    };

    MoFEMFunctionReturn(0);
  };

  auto all_new_adj_entities = [&](const bool create) {
    Range adj;
    for (auto d : {1, 2})
      CHKERR moab.get_adjacencies(new_3d_ents.subset_by_type(MBTET), d, create,
                                  adj, moab::Interface::UNION);
    return adj;
  };

  auto add_new_prisms_which_parents_are_part_of_other_intefaces = [&]() {
    MoFEMFunctionBegin;

    Tag th_interface_side;
    CHKERR moab.tag_get_handle("INTERFACE_SIDE", th_interface_side);

    Range new_3d_prims = new_3d_ents.subset_by_type(MBPRISM);
    for (Range::iterator pit = new_3d_prims.begin(); pit != new_3d_prims.end();
         ++pit) {

      // get parent entity
      EntityHandle parent_prism;
      CHKERR moab.tag_get_data(cOre.get_th_RefParentHandle(), &*pit, 1,
                               &parent_prism);
      if (moab.type_from_handle(parent_prism) != MBPRISM)
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "this prism should have parent which is prism as well");

      int num_nodes;
      // parent prism
      const EntityHandle *conn_parent;
      CHKERR moab.get_connectivity(parent_prism, conn_parent, num_nodes, true);
      Range face_side3_parent, face_side4_parent;
      CHKERR moab.get_adjacencies(conn_parent, 3, 2, false, face_side3_parent);
      CHKERR moab.get_adjacencies(&conn_parent[3], 3, 2, false,
                                  face_side4_parent);
      if (face_side3_parent.size() != 1)
        SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                 "parent face3.size() = %u", face_side3_parent.size());

      if (face_side4_parent.size() != 1)
        SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                 "parent face4.size() = %u", face_side4_parent.size());

      // new prism
      const EntityHandle *conn;
      CHKERR moab.get_connectivity(*pit, conn, num_nodes, true);
      Range face_side3, face_side4;
      CHKERR moab.get_adjacencies(conn, 3, 2, false, face_side3);
      CHKERR moab.get_adjacencies(&conn[3], 3, 2, false, face_side4);
      if (face_side3.size() != 1)
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "face3 is missing");

      if (face_side4.size() != 1)
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "face4 is missing");

      std::vector<EntityHandle> face(2), parent_face(2);
      face[0] = *face_side3.begin();
      face[1] = *face_side4.begin();
      parent_face[0] = *face_side3_parent.begin();
      parent_face[1] = *face_side4_parent.begin();
      for (int ff = 0; ff != 2; ++ff) {
        if (parent_face[ff] == face[ff])
          continue;
        int interface_side;
        CHKERR moab.tag_get_data(th_interface_side, &parent_face[ff], 1,
                                 &interface_side);
        CHKERR moab.tag_set_data(th_interface_side, &face[ff], 1,
                                 &interface_side);
        EntityHandle parent_tri;
        CHKERR moab.tag_get_data(cOre.get_th_RefParentHandle(), &face[ff], 1,
                                 &parent_tri);
        if (parent_tri != parent_face[ff]) {
          SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                   "Wrong parent %lu", parent_tri);
        }
      }
    }
    MoFEMFunctionReturn(0);
  };

  CHKERR create_prisms();

  auto set_parents_ents = all_new_adj_entities(true);

  CHKERR set_parnets(unite(all_new_adj_entities(true), side_ents3d));
  CHKERR add_new_prisms_which_parents_are_part_of_other_intefaces();

  CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
      meshset_for_bit_level, 3, bit);
  CHKERR moab.delete_entities(&meshset_for_bit_level, 1);
  CHKERR moab.clear_meshset(&children[0], 3);

  auto reconstruct_refined_ents = [&]() {
    MoFEMFunctionBegin;
    CHKERR reconstructMultiIndex(*m_field.get_ref_ents());
    MoFEMFunctionReturn(0);
  };

  // Finalise by adding new tets and prism ti bit level
  CHKERR set_parent.override_parents(refined_ents_ptr);

  // Add function which reconstruct core multi-index. Node merging is messy
  // process and entity parent could be changed without notification to
  // multi-index. TODO Issue has to be tracked down better, however in principle
  // is better not to modify multi-index each time parent is changed, that makes
  // code slow. Is better to do it in the bulk as below.
  CHKERR reconstruct_refined_ents();

  MoFEMFunctionReturn(0);
}
} // namespace MoFEM
