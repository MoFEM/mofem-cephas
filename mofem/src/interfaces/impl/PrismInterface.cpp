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
    : cOre(const_cast<Core &>(core)) {}

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
  MoFEMFunctionBegin;
  Range mesh_level_ents3d, mesh_level_ents3d_tris;
  Range mesh_level_tris;
  Range mesh_level_edges;
  Range mesh_level_nodes;
  if (mesh_bit_level.any()) {
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBTET, mesh_level_ents3d);
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBTRI, mesh_level_tris);
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBEDGE, mesh_level_edges);
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBVERTEX, mesh_level_nodes);
    CHKERR moab.get_adjacencies(mesh_level_ents3d, 2, false,
                                mesh_level_ents3d_tris, moab::Interface::UNION);
  }
  Range mesh_level_prisms;
  if (mesh_bit_level.any()) {
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBPRISM, mesh_level_prisms);
    mesh_level_ents3d.merge(mesh_level_prisms);
  }
  Skinner skin(&moab);
  // get interface triangles from side set
  Range triangles;
  CHKERR moab.get_entities_by_type(sideset, MBTRI, triangles, recursive);
  if (mesh_bit_level.any()) {
    triangles = intersect(triangles, mesh_level_ents3d_tris);
  }
  if (verb >= VERBOSE) {
    PetscPrintf(m_field.get_comm(), "Nb. of triangles in set %u\n",
                triangles.size());
  }
  // get nodes, edges and 3d ents (i.e. tets and prisms)
  Range nodes; // nodes from triangles
  CHKERR moab.get_connectivity(triangles, nodes, true);
  Range ents3d, ents3d_with_prisms; // 3d ents form nodes
  CHKERR moab.get_adjacencies(nodes, 3, false, ents3d_with_prisms,
                              moab::Interface::UNION);
  if (mesh_bit_level.any()) {
    ents3d_with_prisms = intersect(ents3d_with_prisms, mesh_level_ents3d);
  }
  ents3d = ents3d_with_prisms.subset_by_type(
      MBTET); // take only tets, add prism later
  // take skin faces
  Range skin_faces; // skin faces from 3d ents
  CHKERR skin.find_skin(0, ents3d, false, skin_faces);
  // take skin edges (boundary of surface if there is any)
  Range skin_edges_boundary; // skin edges from triangles
  CHKERR skin.find_skin(0, triangles, false, skin_edges_boundary);
  if (verb >= VERY_VERBOSE)
    PetscPrintf(m_field.get_comm(), "skin_edges_boundary %u\n",
                skin_edges_boundary.size());
  // take all edges on skin faces (i.e. skin surface)
  Range skin_faces_edges; // edges from skin faces of 3d ents
  CHKERR moab.get_adjacencies(skin_faces, 1, false, skin_faces_edges,
                              moab::Interface::UNION);
  if (mesh_bit_level.any()) {
    skin_faces_edges = intersect(skin_faces_edges, mesh_level_edges);
  }
  if (verb >= VERY_VERBOSE)
    PetscPrintf(m_field.get_comm(), "skin_faces_edges %u\n",
                skin_faces_edges.size());
  // note: that skin faces edges do not contain internal boundary
  // note: that prisms are not included in ents3d, so if ents3d have border with
  // other inteface is like external boundary skin edges boundary are internal
  // edge <- skin_faces_edges contains edges which are on the body boundary <-
  // that is the trick
  skin_edges_boundary =
      subtract(skin_edges_boundary,
               skin_faces_edges); // from skin edges subtract edges from skin
                                  // faces of 3d ents (only internal edges)
  if (verb >= NOISY) {
    EntityHandle out_meshset;
    CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, out_meshset);
    CHKERR moab.add_entities(out_meshset, triangles);
    CHKERR moab.write_file("triangles.vtk", "VTK", "", &out_meshset, 1);
    CHKERR moab.delete_entities(&out_meshset, 1);
    CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, out_meshset);
    CHKERR moab.add_entities(out_meshset, ents3d);
    CHKERR moab.write_file("ents3d.vtk", "VTK", "", &out_meshset, 1);
    CHKERR moab.delete_entities(&out_meshset, 1);
    CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, out_meshset);
    CHKERR moab.add_entities(out_meshset, skin_edges_boundary);
    CHKERR moab.write_file("skin_edges_boundary.vtk", "VTK", "", &out_meshset,
                           1);
    CHKERR moab.delete_entities(&out_meshset, 1);
  }
  if (verb >= VERY_VERBOSE)
    PetscPrintf(m_field.get_comm(), "subtract skin_edges_boundary %u\n",
                skin_edges_boundary.size());
  // Get nodes on boundary edge
  Range skin_nodes_boundary;
  CHKERR moab.get_connectivity(skin_edges_boundary, skin_nodes_boundary, true);
  // Remove node which are boundary with other existing interface
  Range prisms_nodes;
  CHKERR moab.get_connectivity(ents3d_with_prisms.subset_by_type(MBPRISM),
                               prisms_nodes, true);
  skin_nodes_boundary = subtract(skin_nodes_boundary, prisms_nodes);
  if (verb >= VERY_VERBOSE)
    PetscPrintf(m_field.get_comm(), "subtract skin_nodes_boundary %u\n",
                skin_nodes_boundary.size());
  // use nodes on body boundary and interface (without internal boundary) to
  // find adjacent tets
  Range nodes_without_front = subtract(
      nodes, skin_nodes_boundary); // nodes_without_front adjacent to all split
                                   // face edges except those on internal edge
  if (verb >= VERY_VERBOSE)
    PetscPrintf(m_field.get_comm(),
                "adj. node if ents3d but not on the internal edge %u\n",
                nodes_without_front.size());
  // ents3 that are adjacent to front nodes on split faces but not those which
  // are on the front nodes on internal edge
  ents3d.clear();
  ents3d_with_prisms.clear();
  CHKERR moab.get_adjacencies(nodes_without_front, 3, false, ents3d_with_prisms,
                              moab::Interface::UNION);

  // Add tets adjacent to front and triangle which has all nodes on crack front
  struct FindTrianglesOnFrontAndAdjacentTet {
    static MoFEMErrorCode
    fUN(moab::Interface &moab, const Range &triangles,
        const Range &skin_nodes_boundary, ///< nodes on front
        const Range &skin_edges_boundary, ///< edges on front
        Range &ents3d_with_prisms) {
      MoFEMFunctionBegin;
      // get all triangles adjacent to front
      Range skin_nodes_boundary_tris;
      CHKERR moab.get_adjacencies(skin_nodes_boundary, 2, false,
                                  skin_nodes_boundary_tris,
                                  moab::Interface::UNION);
      // get nodes of triangles adjacent to front nodes
      Range skin_nodes_boundary_tris_nodes;
      CHKERR moab.get_connectivity(skin_nodes_boundary_tris,
                                   skin_nodes_boundary_tris_nodes, true);
      // get hanging nodes, i.e. nodes which are not on the front but adjacent
      // to triangles adjacent to crack front
      skin_nodes_boundary_tris_nodes =
          subtract(skin_nodes_boundary_tris_nodes, skin_nodes_boundary);
      // get triangles adjacent to hanging nodes
      Range skin_nodes_boundary_tris_nodes_tris;
      CHKERR moab.get_adjacencies(skin_nodes_boundary_tris_nodes, 2, false,
                                  skin_nodes_boundary_tris_nodes_tris,
                                  moab::Interface::UNION);
      // triangles which have tree nodes on front boundary
      skin_nodes_boundary_tris =
          intersect(triangles, subtract(skin_nodes_boundary_tris,
                                        skin_nodes_boundary_tris_nodes_tris));
      if (!skin_nodes_boundary_tris.empty()) {
        // Get internal edges of triangle which has three nodes on boundary
        Range skin_nodes_boundary_tris_edges;
        CHKERR moab.get_adjacencies(skin_nodes_boundary_tris, 1, false,
                                    skin_nodes_boundary_tris_edges,
                                    moab::Interface::UNION);
        skin_nodes_boundary_tris_edges =
            subtract(skin_nodes_boundary_tris_edges, skin_edges_boundary);
        // Get 3d elements adjacent to internal edge which has two nodes on
        // boundary
        CHKERR moab.get_adjacencies(skin_nodes_boundary_tris_edges, 3, false,
                                    ents3d_with_prisms, moab::Interface::UNION);
      }
      MoFEMFunctionReturn(0);
    }
  };
  CHKERR FindTrianglesOnFrontAndAdjacentTet::fUN(
      moab, triangles, skin_nodes_boundary, skin_edges_boundary,
      ents3d_with_prisms);

  // prism and tets on both side of interface
  if (mesh_bit_level.any()) {
    ents3d_with_prisms = intersect(ents3d_with_prisms, mesh_level_ents3d);
  }
  ents3d = ents3d_with_prisms.subset_by_type(MBTET);
  if (verb >= VERY_VERBOSE)
    PetscPrintf(m_field.get_comm(), "adj. ents3d to front nodes %u\n",
                ents3d.size());

  Range side_ents3d;
  unsigned int nb_side_ents3d = side_ents3d.size();
  side_ents3d.insert(*ents3d.begin());
  Range side_ents3d_tris_on_surface;

  // get all tets adjacent to crack surface, but only on one side of it
  do {

    do {
      Range adj_tris, adj_ents3d;
      nb_side_ents3d = side_ents3d.size();
      if (verb >= VERBOSE)
        PetscPrintf(m_field.get_comm(), "nb_side_ents3d %u\n", nb_side_ents3d);
      // get faces
      CHKERR moab.get_adjacencies(side_ents3d.subset_by_type(MBTET), 2, false,
                                  adj_tris, moab::Interface::UNION);
      if (mesh_bit_level.any()) {
        adj_tris = intersect(adj_tris, mesh_level_tris);
      }
      // subtrace from faces interface
      adj_tris = subtract(adj_tris, triangles);
      if (verb >= VERBOSE)
        PetscPrintf(m_field.get_comm(), "adj_tris %u\n", adj_tris.size());
      // get tets adjacent to faces
      CHKERR moab.get_adjacencies(adj_tris, 3, true, adj_ents3d,
                                  moab::Interface::UNION);
      // intersect tets with tets adjacent to inetface
      adj_ents3d = intersect(adj_ents3d, ents3d_with_prisms);
      if (verb >= VERBOSE)
        PetscPrintf(m_field.get_comm(), "adj_ents3d %u\n", adj_ents3d.size());
      // add tets to side
      side_ents3d.insert(adj_ents3d.begin(), adj_ents3d.end());
      if (verb >= VERY_NOISY) {
        EntityHandle out_meshset;
        CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER,
                                   out_meshset);
        CHKERR moab.add_entities(out_meshset, side_ents3d);
        std::string file = "side_ents3d_" +
                           boost::lexical_cast<std::string>(nb_side_ents3d) +
                           ".vtk";
        CHKERR moab.write_file(file.c_str(), "VTK", "", &out_meshset, 1);
        CHKERR moab.delete_entities(&out_meshset, 1);
      }
    } while (nb_side_ents3d != side_ents3d.size());

    Range side_ents3d_tris;
    CHKERR moab.get_adjacencies(side_ents3d, 2, false, side_ents3d_tris,
                                moab::Interface::UNION);
    side_ents3d_tris_on_surface = intersect(side_ents3d_tris, triangles);

    if (verb >= VERY_NOISY) {
      Range left_triangles = subtract(triangles, side_ents3d_tris_on_surface);
      if (!left_triangles.empty()) {
        EntityHandle out_meshset;
        CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER,
                                   out_meshset);
        CHKERR moab.add_entities(out_meshset, left_triangles);
        std::string file = "left_triangles_" +
                           boost::lexical_cast<std::string>(nb_side_ents3d) +
                           ".vtk";
        CHKERR moab.write_file(file.c_str(), "VTK", "", &out_meshset, 1);
        CHKERR moab.delete_entities(&out_meshset, 1);
      }
    }

    // This is a case when separate sub-domains are split, so we need
    // additional tetrahedron for seed process
    if (side_ents3d_tris_on_surface.size() != triangles.size()) {
      Range left_triangles = subtract(triangles, side_ents3d_tris_on_surface);
      Range tets;
      CHKERR moab.get_adjacencies(&*left_triangles.begin(), 1, 3, false, tets);

      tets = intersect(tets, ents3d_with_prisms);
      if (tets.empty()) {
        Range left_triangles_nodes;
        CHKERR moab.get_connectivity(&*left_triangles.begin(), 1,
                                     left_triangles_nodes, true);
        EntityHandle meshset;
        CHKERR moab.create_meshset(MESHSET_SET, meshset);
        CHKERR moab.add_entities(meshset, left_triangles);
        CHKERR moab.write_file("error.vtk", "VTK", "", &meshset, 1);
        CHKERR moab.delete_entities(&meshset, 1);
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "Not all faces on surface going to be split, see error.vtk for "
                "problematic triangle. "
                "It could be a case where triangle on front (part boundary of "
                "surface in interior) "
                "has three nodes front.");
      }
      side_ents3d.insert(*tets.begin());
    }

  } while (side_ents3d_tris_on_surface.size() != triangles.size());

  if (ents3d_with_prisms.size() == side_ents3d.size()) {
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "All tets on one side, no-interface");
  }
  // other side ents
  Range other_side = subtract(ents3d_with_prisms, side_ents3d);
  // side nodes
  Range side_nodes;
  CHKERR moab.get_connectivity(side_ents3d.subset_by_type(MBTET), side_nodes,
                               true);
  // nodes on crack surface without front
  nodes_without_front = intersect(nodes_without_front, side_nodes);
  Range side_edges;
  CHKERR moab.get_adjacencies(side_ents3d.subset_by_type(MBTET), 1, false,
                              side_edges, moab::Interface::UNION);
  skin_edges_boundary = intersect(skin_edges_boundary, side_edges);
  // make child meshsets
  std::vector<EntityHandle> children;
  CHKERR moab.get_child_meshsets(sideset, children);
  if (children.empty()) {
    children.resize(3);
    CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, children[0]);
    CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, children[1]);
    CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, children[2]);
    CHKERR moab.add_child_meshset(sideset, children[0]);
    CHKERR moab.add_child_meshset(sideset, children[1]);
    CHKERR moab.add_child_meshset(sideset, children[2]);
  } else {
    if (children.size() != 3) {
      SETERRQ(PETSC_COMM_SELF, 1,
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
  if (verb >= VERBOSE) {
    PetscPrintf(m_field.get_comm(), "Nb. of side ents3d in set %u\n",
                side_ents3d.size());
    PetscPrintf(m_field.get_comm(), "Nb. of other side ents3d in set %u\n",
                other_side.size());
    PetscPrintf(m_field.get_comm(), "Nb. of boundary edges %u\n",
                skin_edges_boundary.size());
  }
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
  if (verb >= VERBOSE) {
    PetscPrintf(m_field.get_comm(), "Nb. of triangles in set %u\n",
                triangles.size());
  }
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

  if (verb >= VERY_VERBOSE) {
    PetscPrintf(m_field.get_comm(), "skin_edges_boundary %u\n",
                skin_edges_boundary.size());
  }
  // take all edges on skin faces (i.e. skin surface)
  Range skin_faces_edges; // edges from skin faces of 3d ents
  CHKERR moab.get_adjacencies(skin_faces, 1, false, skin_faces_edges,
                              moab::Interface::UNION);

  if (mesh_bit_level.any()) {
    skin_faces_edges = intersect(skin_faces_edges, mesh_level_edges);
  }
  if (verb >= VERY_VERBOSE) {
    PetscPrintf(m_field.get_comm(), "skin_faces_edges %u\n",
                skin_faces_edges.size());
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
    CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, out_meshset);
    CHKERR moab.add_entities(out_meshset, triangles);
    CHKERR moab.write_file("triangles.vtk", "VTK", "", &out_meshset, 1);
    CHKERR moab.delete_entities(&out_meshset, 1);
    CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, out_meshset);
    CHKERR moab.add_entities(out_meshset, ents3d);
    CHKERR moab.write_file("ents3d.vtk", "VTK", "", &out_meshset, 1);
    CHKERR moab.delete_entities(&out_meshset, 1);
    CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, out_meshset);
    CHKERR moab.add_entities(out_meshset, skin_edges_boundary);
    CHKERR moab.write_file("skin_edges_boundary.vtk", "VTK", "", &out_meshset,
                           1);
    CHKERR moab.delete_entities(&out_meshset, 1);
  }
  if (verb >= VERY_VERBOSE) {
    PetscPrintf(m_field.get_comm(), "subtract skin_edges_boundary %u\n",
                skin_edges_boundary.size());
  }

  // Get nodes on boundary edge
  Range skin_nodes_boundary;
  CHKERR moab.get_connectivity(skin_edges_boundary, skin_nodes_boundary, true);

  // Remove node which are boundary with other existing interface
  Range prisms_nodes;
  CHKERR
  moab.get_connectivity(ents3d.subset_by_type(MBPRISM), prisms_nodes, true);

  skin_nodes_boundary = subtract(skin_nodes_boundary, prisms_nodes);
  if (verb >= VERY_VERBOSE) {
    PetscPrintf(m_field.get_comm(), "subtract skin_nodes_boundary %u\n",
                skin_nodes_boundary.size());
  }
  // use nodes on body boundary and interface (without internal boundary) to
  // find adjacent tets
  Range nodes_without_front = subtract(
      nodes, skin_nodes_boundary); // nodes_without_front adjacent to all split
                                   // face edges except those on internal edge
  if (verb >= VERY_VERBOSE) {
    PetscPrintf(m_field.get_comm(),
                "adj. node if ents3d but not on the internal edge %u\n",
                nodes_without_front.size());
  }

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
}

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
  const RefEntity_multiIndex *refined_ents_ptr;
  MoFEMFunctionBegin;

  CHKERR m_field.get_ref_ents(&refined_ents_ptr);

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
  if (verb >= VERY_VERBOSE) {
    PetscPrintf(m_field.get_comm(), "split sides triangles %u\n",
                triangles.size());
    PetscPrintf(m_field.get_comm(), "split sides side_ents3d %u\n",
                side_ents3d.size());
    PetscPrintf(m_field.get_comm(), "split sides nodes %u\n", nodes.size());
  }

  typedef std::map<EntityHandle, /*node on "mother" side*/
                   EntityHandle  /*node on "father" side*/
                   >
      MapNodes;
  MapNodes map_nodes;

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

      MoFEMErrorCode operator()(const BitRefLevel &bit, MapNodes &map_nodes) {
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
        for (int nn = 0; nn != num_nodes; ++nn)
          map_nodes[splitNodes[nn]] = verts[nn];
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
    CHKERR create_side_nodes(bit, map_nodes);
  }

  // crete meshset for new mesh bit level
  EntityHandle meshset_for_bit_level;
  CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER,
                             meshset_for_bit_level);
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
          PetscPrintf(m_field.get_comm(), "nodes %u -> %d\n", conn[ii],
                      new_conn[ii]);
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
          // FIXME: That takes firs element form the list. Should throw error
          // if is more than one or handle it properly.
          RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator new_rit =
              refined_ents_ptr->get<Ent_mi_tag>().find(*new_conn_tet.begin());
          if (new_rit == refined_ents_ptr->get<Ent_mi_tag>().end())
            SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "Can't find entity in database");
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

  auto get_adj_ents = [&](const bool create) {
    Range adj;
    // create new entities by adjacencies form new tets
    CHKERR moab.get_adjacencies(new_3d_ents.subset_by_type(MBTET), 2, create,
                                adj, moab::Interface::UNION);
    CHKERR moab.get_adjacencies(new_3d_ents.subset_by_type(MBTET), 1, create,
                                adj, moab::Interface::UNION);
    return adj;
  };

  Range new_ents_existing = get_adj_ents(false);
  Range new_ents = subtract(get_adj_ents(true), new_ents_existing);

  // Tags for setting side
  Tag th_interface_side;
  const int def_side[] = {0};
  CHKERR moab.tag_get_handle("INTERFACE_SIDE", 1, MB_TYPE_INTEGER,
                             th_interface_side, MB_TAG_CREAT | MB_TAG_SPARSE,
                             def_side);

  struct SetParent {
    map<EntityHandle, EntityHandle> parentsToChange;
    MoFEMErrorCode operator()(const EntityHandle ent, const EntityHandle parent,
                              const RefEntity_multiIndex *ref_ents_ptr,
                              MoFEM::Core &cOre) {
      MoFEM::Interface &m_field = cOre;
      MoFEMFunctionBegin;
      auto it = ref_ents_ptr->find(ent);
      if (it != ref_ents_ptr->end()) {
        if (it->get()->getParentEnt() != parent && ent != parent)
          parentsToChange[ent] = parent;
      } else {
        if (ent != parent)
          CHKERR m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(),
                                                 &ent, 1, &parent);
      }
      MoFEMFunctionReturn(0);
    }
    MoFEMErrorCode operator()(const RefEntity_multiIndex *ref_ents_ptr) {
      MoFEMFunctionBegin;
      for (auto mit = parentsToChange.begin(); mit != parentsToChange.end();
           ++mit) {
        auto it = ref_ents_ptr->find(mit->first);
        bool success = const_cast<RefEntity_multiIndex *>(ref_ents_ptr)
                           ->modify(it, RefEntity_change_parent(mit->second));
        if (!success)
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Impossible to set parent");
      }
      MoFEMFunctionReturn(0);
    }
  };
  SetParent set_parent;

  // add new edges and triangles to mofem database
  Range ents;
  CHKERR moab.get_adjacencies(triangles, 1, false, ents,
                              moab::Interface::UNION);
  ents.insert(triangles.begin(), triangles.end());
  for (Range::iterator eit = ents.begin(); eit != ents.end(); eit++) {
    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(*eit, conn, num_nodes, true);
    int sense = 0; ///< sense of the triangle used to create a prism
    if (moab.type_from_handle(*eit) == MBTRI) {
      Range ents_3d;
      CHKERR moab.get_adjacencies(&*eit, 1, 3, false, ents_3d);
      ents_3d = intersect(ents_3d, side_ents3d);
      switch (ents_3d.size()) {
      case 0:
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "Did not find adjacent tets on one side of the interface; if "
                "this error appears for a contact interface, try creating "
                "separate blocksets for each contact surface");
      case 2:
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "Found both adjacent tets on one side of the interface, if "
                "this error appears for a contact interface, try creating "
                "separate blocksets for each contact surface");
      default:
        break;
      }
      int side, offset;
      CHKERR moab.side_number(ents_3d.front(), *eit, side, sense, offset);
    }
    EntityHandle new_conn[num_nodes];
    int nb_new_conn = 0;
    for (int ii = 0; ii != num_nodes; ++ii) {
      std::map<EntityHandle, EntityHandle>::iterator mit =
          map_nodes.find(conn[ii]);
      if (mit != map_nodes.end()) {
        new_conn[ii] = mit->second;
        nb_new_conn++;
        if (verb >= VERY_NOISY) {
          PetscPrintf(m_field.get_comm(), "nodes %u -> %d\n", conn[ii],
                      new_conn[ii]);
        }
      } else {
        new_conn[ii] = conn[ii];
      }
    }
    if (nb_new_conn == 0)
      continue;
    RefEntity_multiIndex::iterator miit_ref_ent = refined_ents_ptr->find(*eit);
    if (miit_ref_ent == refined_ents_ptr->end())
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "this entity (edge or tri) should be already in database");

    Range new_ent; // contains all entities (edges or triangles) added to mofem
                   // database
    switch (moab.type_from_handle(*eit)) {
    case MBTRI: {
      // get entity based on its connectivity
      CHKERR moab.get_adjacencies(new_conn, 3, 2, false, new_ent);
      if (new_ent.size() != 1)
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "this tri should be in moab database");
      int new_side = 1;
      CHKERR moab.tag_set_data(th_interface_side, &*new_ent.begin(), 1,
                               &new_side);
      if (verb >= VERY_VERBOSE)
        PetscPrintf(m_field.get_comm(), "new_ent %u\n", new_ent.size());
      // add prism element
      if (add_interface_entities) {
        if (inhered_from_bit_level.any()) {
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "not implemented for inhered_from_bit_level");
        }
        // set prism connectivity
        EntityHandle prism_conn[6] = {conn[0],     conn[1],     conn[2],
                                      new_conn[0], new_conn[1], new_conn[2]};
        if (sense != 1 && sense != -1) {
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "Undefined sense of a triangle; if this error appears for a "
                  "contact interface, try creating separate blocksets for each "
                  "contact surface");
        }
        if (sense == -1) {
          // swap nodes in triangles for correct prism creation
          std::swap(prism_conn[1], prism_conn[2]);
          std::swap(prism_conn[4], prism_conn[5]);
        }
        EntityHandle prism;
        CHKERR moab.create_element(MBPRISM, prism_conn, 6, prism);
        CHKERR moab.add_entities(meshset_for_bit_level, &prism, 1);
      }
    } break;
    case MBEDGE: {
      CHKERR moab.get_adjacencies(new_conn, 2, 1, false, new_ent);
      if (new_ent.size() != 1) {
        SETERRQ2(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                 "this edge should be in moab database "
                 "new_ent.size() = %u nb_new_conn = %d",
                 new_ent.size(), nb_new_conn);
      }
    } break;
    default:
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "Houston we have a problem !!!");
    }
    if (new_ent.size() != 1) {
      SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "new_ent.size() = %u, size always should be 1", new_ent.size());
    }
    CHKERR set_parent(new_ent[0], *eit, refined_ents_ptr, cOre);
  }

  // all other entities, some ents like triangles and faces on the side of tets
  auto all_others_adj_entities = [&](const bool create) {
    Range adj;
    for (auto d : {1, 2})
      CHKERR moab.get_adjacencies(side_ents3d.subset_by_type(MBTET), d, create,
                                  adj, moab::Interface::UNION);
    return adj;
  };
  Range side_adj_faces_and_edges = all_others_adj_entities(true);
  
  for (Range::iterator eit = side_adj_faces_and_edges.begin();
       eit != side_adj_faces_and_edges.end(); ++eit) {
    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(*eit, conn, num_nodes, true);
    EntityHandle new_conn[num_nodes];
    int nb_new_conn = 0;
    for (int ii = 0; ii < num_nodes; ii++) {
      std::map<EntityHandle, EntityHandle>::iterator mit =
          map_nodes.find(conn[ii]);
      if (mit != map_nodes.end()) {
        new_conn[ii] = mit->second;
        nb_new_conn++;
        if (verb >= VERY_NOISY) {
          PetscPrintf(m_field.get_comm(), "nodes %u -> %d\n", conn[ii],
                      new_conn[ii]);
        }
      } else {
        new_conn[ii] = conn[ii];
      }
    }
    if (nb_new_conn == 0)
      continue;
    RefEntity_multiIndex::iterator miit_ref_ent = refined_ents_ptr->find(*eit);
    if (miit_ref_ent == refined_ents_ptr->end()) {
      SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "entity should be in MoFem database, num_nodes = %d", num_nodes);
    }
    Range new_ent;
    switch (moab.type_from_handle(*eit)) {
    case MBEDGE:
      CHKERR moab.get_adjacencies(new_conn, 2, 1, false, new_ent);
      break;
    case MBTRI:
      CHKERR moab.get_adjacencies(new_conn, 3, 2, false, new_ent);
      break;
    default:
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "Houston we have a problem");
    }
    if (new_ent.size() != 1) {
      SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "database inconsistency, new_ent.size() = %u", new_ent.size());
    }
    CHKERR set_parent(new_ent[0], *eit, refined_ents_ptr, cOre);
  }

  // add new prisms which parents are part of other intefaces
  Range new_3d_prims = new_3d_ents.subset_by_type(MBPRISM);
  for (Range::iterator pit = new_3d_prims.begin(); pit != new_3d_prims.end();
       ++pit) {
    // get parent entity
    EntityHandle parent_prism;
    CHKERR moab.tag_get_data(cOre.get_th_RefParentHandle(), &*pit, 1,
                             &parent_prism);
    if (moab.type_from_handle(parent_prism) != MBPRISM) {
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "this prism should have parent which is prism as well");
    }
    int num_nodes;
    // parent prism
    const EntityHandle *conn_parent;
    CHKERR moab.get_connectivity(parent_prism, conn_parent, num_nodes, true);
    Range face_side3_parent, face_side4_parent;
    CHKERR moab.get_adjacencies(conn_parent, 3, 2, false, face_side3_parent);
    CHKERR moab.get_adjacencies(&conn_parent[3], 3, 2, false,
                                face_side4_parent);
    if (face_side3_parent.size() != 1) {
      SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "parent face3.size() = %u", face_side3_parent.size());
    }
    if (face_side4_parent.size() != 1) {
      SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "parent face4.size() = %u", face_side4_parent.size());
    }
    // new prism
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(*pit, conn, num_nodes, true);
    Range face_side3, face_side4;
    CHKERR moab.get_adjacencies(conn, 3, 2, false, face_side3);
    CHKERR moab.get_adjacencies(&conn[3], 3, 2, false, face_side4);
    if (face_side3.size() != 1) {
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY, "face3 is missing");
    }
    if (face_side4.size() != 1) {
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY, "face4 is missing");
    }
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
                 "wrong parent %lu", parent_tri);
      }
    }
  }

  // finalise by adding new tets and prism ti bit level
  // FIXME: This is switch of, you can not change parent. 
  // CHKERR set_parent(refined_ents_ptr);

  CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
      meshset_for_bit_level, 3, bit);
  CHKERR moab.delete_entities(&meshset_for_bit_level, 1);
  CHKERR moab.clear_meshset(&children[0], 3);
  MoFEMFunctionReturn(0);
}
} // namespace MoFEM
