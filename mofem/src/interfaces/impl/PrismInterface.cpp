/** \file PrismInterface.cpp
 * \brief Inserting prims interface elements
 * \todo FIXME this is no so good implementation
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
  if(uuid == IDD_MOFEMPrismInterface) {
    *iface = const_cast<PrismInterface*>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
  MoFEMFunctionReturnHot(0);
}

PrismInterface::PrismInterface(const Core &core):
cOre(const_cast<Core&>(core)) {
}

MoFEMErrorCode PrismInterface::getSides(const int msId,
                                        const CubitBCType cubit_bc_type,
                                        const BitRefLevel mesh_bit_level,
                                        const bool recursive, int verb) {

  Interface &m_field = cOre;
  MeshsetsManager *meshsets_manager_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.getInterface(meshsets_manager_ptr);
  CHKERRQ(ierr);
  CubitMeshSet_multiIndex::index<
      Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator miit =
      meshsets_manager_ptr->getMeshsetsMultindex()
          .get<Composite_Cubit_msId_And_MeshSetType_mi_tag>()
          .find(boost::make_tuple(msId, cubit_bc_type.to_ulong()));
  if (miit !=
      meshsets_manager_ptr->getMeshsetsMultindex()
          .get<Composite_Cubit_msId_And_MeshSetType_mi_tag>()
          .end()) {
    ierr = getSides(miit->meshset, mesh_bit_level, recursive, verb);
    CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"msId is not there");
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode PrismInterface::getSides(const EntityHandle sideset,
                                        const BitRefLevel mesh_bit_level,
                                        const bool recursive, int verb) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  Range mesh_level_ents3d,mesh_level_ents3d_tris;
  Range mesh_level_tris;
  Range mesh_level_edges;
  Range mesh_level_nodes;
  if (mesh_bit_level.any()) {
    ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBTET, mesh_level_ents3d);
    CHKERRQ(ierr);
    ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBTRI, mesh_level_tris);
    CHKERRQ(ierr);
    ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBEDGE, mesh_level_edges);
    CHKERRQ(ierr);
    ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBVERTEX, mesh_level_nodes);
    CHKERRQ(ierr);
    rval = moab.get_adjacencies(mesh_level_ents3d, 2, false,
                                mesh_level_ents3d_tris, moab::Interface::UNION);
    CHKERRQ_MOAB(rval);
  }
  Range mesh_level_prisms;
  if(mesh_bit_level.any()) {
    ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBPRISM, mesh_level_prisms);
    CHKERRQ(ierr);
    mesh_level_ents3d.merge(mesh_level_prisms);
  }
  Skinner skin(&moab);
  //get interface triangles from side set
  Range triangles;
  rval = moab.get_entities_by_type(sideset, MBTRI, triangles, recursive);
  CHKERRQ_MOAB(rval);
  if(mesh_bit_level.any()) {
    triangles = intersect(triangles,mesh_level_ents3d_tris);
  }
  if(verb>=VERBOSE) {
    PetscPrintf(m_field.get_comm(), "Nb. of triangles in set %u\n",
                triangles.size());
  }
  //get nodes, edges and 3d ents (i.e. tets and prisms)
  Range nodes; // nodes from triangles
  rval = moab.get_connectivity(triangles,nodes,true); CHKERRQ_MOAB(rval);
  Range ents3d,ents3d_with_prisms; // 3d ents form nodes
  rval = moab.get_adjacencies(nodes, 3, false, ents3d_with_prisms,
                              moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  if(mesh_bit_level.any()) {
    ents3d_with_prisms = intersect(ents3d_with_prisms,mesh_level_ents3d);
  }
  ents3d = ents3d_with_prisms.subset_by_type(MBTET); // take only tets, add prism later
  //take skin faces
  Range skin_faces; // skin faces from 3d ents
  rval = skin.find_skin(0,ents3d,false,skin_faces); CHKERRQ_MOAB(rval);
  //take skin edges (boundary of surface if there is any)
  Range skin_edges_boundary; //skin edges from triangles
  rval = skin.find_skin(0, triangles, false, skin_edges_boundary);
  CHKERRQ_MOAB(rval);
  if (verb >= VERY_VERBOSE)
    PetscPrintf(m_field.get_comm(), "skin_edges_boundary %u\n",
                skin_edges_boundary.size());
  //take all edges on skin faces (i.e. skin surface)
  Range skin_faces_edges; //edges from skin faces of 3d ents
  rval = moab.get_adjacencies(skin_faces, 1, false, skin_faces_edges,
                              moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  if(mesh_bit_level.any()) {
    skin_faces_edges = intersect(skin_faces_edges,mesh_level_edges);
  }
  if (verb >= VERY_VERBOSE)
    PetscPrintf(m_field.get_comm(), "skin_faces_edges %u\n",
                skin_faces_edges.size());
  //note: that skin faces edges do not contain internal boundary
  //note: that prisms are not included in ents3d, so if ents3d have border with other inteface is like external boundary
  //skin edges boundary are internal edge <- skin_faces_edges contains edges which are on the body boundary <- that is the trick
  skin_edges_boundary = subtract(skin_edges_boundary,skin_faces_edges); // from skin edges subtract edges from skin faces of 3d ents (only internal edges)
  if (verb >= NOISY) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, out_meshset);
    CHKERRQ_MOAB(rval);
    rval = moab.add_entities(out_meshset, triangles);
    CHKERRQ_MOAB(rval);
    rval = moab.write_file("triangles.vtk", "VTK", "", &out_meshset, 1);
    CHKERRQ_MOAB(rval);
    rval = moab.delete_entities(&out_meshset, 1);
    CHKERRQ_MOAB(rval);
    rval = moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, out_meshset);
    CHKERRQ_MOAB(rval);
    rval = moab.add_entities(out_meshset, ents3d);
    CHKERRQ_MOAB(rval);
    rval = moab.write_file("ents3d.vtk", "VTK", "", &out_meshset, 1);
    CHKERRQ_MOAB(rval);
    rval = moab.delete_entities(&out_meshset, 1);
    CHKERRQ_MOAB(rval);
    rval = moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, out_meshset);
    CHKERRQ_MOAB(rval);
    rval = moab.add_entities(out_meshset, skin_edges_boundary);
    CHKERRQ_MOAB(rval);
    rval =
        moab.write_file("skin_edges_boundary.vtk", "VTK", "", &out_meshset, 1);
    CHKERRQ_MOAB(rval);
    rval = moab.delete_entities(&out_meshset, 1);
    CHKERRQ_MOAB(rval);
  }
  if (verb >= VERY_VERBOSE)
    PetscPrintf(m_field.get_comm(), "subtract skin_edges_boundary %u\n",
                skin_edges_boundary.size());
  //Get nodes on boundary edge
  Range skin_nodes_boundary;
  rval = moab.get_connectivity(skin_edges_boundary, skin_nodes_boundary, true);
  CHKERRQ_MOAB(rval);
  //Remove node which are boundary with other existing interface
  Range prisms_nodes;
  rval = moab.get_connectivity(ents3d_with_prisms.subset_by_type(MBPRISM),
                               prisms_nodes, true);
  CHKERRQ_MOAB(rval);
  skin_nodes_boundary = subtract(skin_nodes_boundary,prisms_nodes);
  if (verb >= VERY_VERBOSE)
    PetscPrintf(m_field.get_comm(), "subtract skin_nodes_boundary %u\n",
                skin_nodes_boundary.size());
  //use nodes on body boundary and interface (without internal boundary) to find adjacent tets
  Range nodes_without_front = subtract(nodes,skin_nodes_boundary); // nodes_without_front adjacent to all split face edges except those on internal edge
  if (verb >= VERY_VERBOSE)
    PetscPrintf(m_field.get_comm(),
                "adj. node if ents3d but not on the internal edge %u\n",
                nodes_without_front.size());
  //ents3 that are adjacent to front nodes on split faces but not those which are on the front nodes on internal edge
  ents3d.clear();
  ents3d_with_prisms.clear();
  rval = moab.get_adjacencies(
    nodes_without_front,3,false,ents3d_with_prisms,moab::Interface::UNION
  ); CHKERRQ_MOAB(rval);

  // Add tets adjacent to front and triangle which has all nodes on crack front
  struct FindTrianglesOnFrontAndAdjacentTet {
    static MoFEMErrorCode fUN(moab::Interface &moab, const Range &triangles,
                              const Range &skin_nodes_boundary, ///< nodes on front
                              const Range &skin_edges_boundary, ///< edges on front
                              Range &ents3d_with_prisms) {
      MoFEMFunctionBeginHot;
      // get all triangles adjacent to front
      Range skin_nodes_boundary_tris;
      rval = moab.get_adjacencies(skin_nodes_boundary, 2, false,
                                  skin_nodes_boundary_tris,
                                  moab::Interface::UNION);
      CHKERRQ(ierr);
      // get nodes of triangles adjacent to front nodes
      Range skin_nodes_boundary_tris_nodes;
      rval = moab.get_connectivity(skin_nodes_boundary_tris,
                                   skin_nodes_boundary_tris_nodes, true);
      CHKERRQ_MOAB(rval);
      // get hanging nodes, i.e. nodes which are not on the front but adjacent
      // to triangles adjacent to crack front
      skin_nodes_boundary_tris_nodes =
          subtract(skin_nodes_boundary_tris_nodes, skin_nodes_boundary);
      // get triangles adjacent to hanging nodes
      Range skin_nodes_boundary_tris_nodes_tris;
      rval = moab.get_adjacencies(skin_nodes_boundary_tris_nodes, 2, false,
                                  skin_nodes_boundary_tris_nodes_tris,
                                  moab::Interface::UNION);
      CHKERRQ(ierr);
      // triangles which have tree nodes on front boundary
      skin_nodes_boundary_tris =
          intersect(triangles, subtract(skin_nodes_boundary_tris,
                                        skin_nodes_boundary_tris_nodes_tris));
      if (!skin_nodes_boundary_tris.empty()) {
        // Get internal edges of triangle which has three nodes on boundary
        Range skin_nodes_boundary_tris_edges;
        rval = moab.get_adjacencies(skin_nodes_boundary_tris, 1, false,
                                    skin_nodes_boundary_tris_edges,
                                    moab::Interface::UNION);
        CHKERRQ(rval);
        skin_nodes_boundary_tris_edges =
            subtract(skin_nodes_boundary_tris_edges, skin_edges_boundary);
        // Get 3d elements adjacent to internal edge which has two nodes on
        // boundary
        rval = moab.get_adjacencies(skin_nodes_boundary_tris_edges, 3, false,
                                    ents3d_with_prisms, moab::Interface::UNION);
        CHKERRQ(rval);
      }
      MoFEMFunctionReturnHot(0);
    }
  };
  ierr = FindTrianglesOnFrontAndAdjacentTet::fUN(
    moab,triangles,skin_nodes_boundary,skin_edges_boundary,ents3d_with_prisms
  ); CHKERRQ(ierr);

  // prism and tets on both side of interface
  if(mesh_bit_level.any()) {
    ents3d_with_prisms = intersect(ents3d_with_prisms,mesh_level_ents3d);
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
      Range adj_tris,adj_ents3d;
      nb_side_ents3d = side_ents3d.size();
      if (verb >= VERBOSE)
        PetscPrintf(m_field.get_comm(), "nb_side_ents3d %u\n", nb_side_ents3d);
      //get faces
      rval = moab.get_adjacencies(side_ents3d.subset_by_type(MBTET), 2, false,
                                  adj_tris, moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      if(mesh_bit_level.any()) {
        adj_tris = intersect(adj_tris,mesh_level_tris);
      }
      //subtrace from faces interface
      adj_tris = subtract(adj_tris,triangles);
      if (verb >= VERBOSE)
        PetscPrintf(m_field.get_comm(), "adj_tris %u\n", adj_tris.size());
      //get tets adjacent to faces
      rval = moab.get_adjacencies(
        adj_tris,3,true,adj_ents3d,moab::Interface::UNION
      ); CHKERRQ_MOAB(rval);
      //intersect tets with tets adjacent to inetface
      adj_ents3d = intersect(adj_ents3d,ents3d_with_prisms);
      if (verb >= VERBOSE)
        PetscPrintf(m_field.get_comm(), "adj_ents3d %u\n", adj_ents3d.size());
      //add tets to side
      side_ents3d.insert(adj_ents3d.begin(),adj_ents3d.end());
      if (verb >= VERY_NOISY) {
        EntityHandle out_meshset;
        rval =
            moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, out_meshset);
        CHKERRQ_MOAB(rval);
        rval = moab.add_entities(out_meshset,side_ents3d);
        CHKERRQ_MOAB(rval);
        std::string file = "side_ents3d_" +
                           boost::lexical_cast<std::string>(nb_side_ents3d) +
                           ".vtk";
        rval = moab.write_file(file.c_str(), "VTK", "", &out_meshset, 1);
        CHKERRQ_MOAB(rval);
        rval = moab.delete_entities(&out_meshset, 1);
        CHKERRQ_MOAB(rval);
      }
    } while (nb_side_ents3d != side_ents3d.size());

    Range side_ents3d_tris;
    rval = moab.get_adjacencies(
      side_ents3d,2,false,side_ents3d_tris,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    side_ents3d_tris_on_surface = intersect(side_ents3d_tris,triangles);

    if (verb >= VERY_NOISY) {
      Range left_triangles = subtract(triangles, side_ents3d_tris_on_surface);
      if(!left_triangles.empty()) {
        EntityHandle out_meshset;
        rval =
            moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, out_meshset);
        CHKERRQ_MOAB(rval);
        rval = moab.add_entities(out_meshset, left_triangles);
        CHKERRQ_MOAB(rval);
        std::string file = "left_triangles_" +
                           boost::lexical_cast<std::string>(nb_side_ents3d) +
                           ".vtk";
        rval = moab.write_file(file.c_str(), "VTK", "", &out_meshset, 1);
        CHKERRQ_MOAB(rval);
        rval = moab.delete_entities(&out_meshset, 1);
        CHKERRQ_MOAB(rval);
      }
    }

      // This is a case when separate sub-domains are split, so wee need
      // additional
      // tetrahedron for seed process
    if (side_ents3d_tris_on_surface.size() != triangles.size()) {
      Range left_triangles = subtract(triangles, side_ents3d_tris_on_surface);
      Range tets;
      rval = moab.get_adjacencies(&*left_triangles.begin(), 1, 3, false, tets);
      CHKERRQ_MOAB(rval);
      tets = intersect(tets, ents3d_with_prisms);
      if (tets.empty()) {
        Range left_triangles_nodes;
        rval = moab.get_connectivity(&*left_triangles.begin(), 1,
                                     left_triangles_nodes, true);
        CHKERRQ(ierr);
        EntityHandle meshset;
        rval = moab.create_meshset(MESHSET_SET, meshset);
        CHKERRQ_MOAB(rval);
        rval = moab.add_entities(meshset, left_triangles);
        CHKERRQ_MOAB(rval);
        rval = moab.write_file("error.vtk", "VTK", "", &meshset, 1);
        CHKERRQ_MOAB(rval);
        rval = moab.delete_entities(&meshset, 1);
        CHKERRQ_MOAB(rval);
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "Not all faces on surface going to be split, see error.vtk for "
                "problematic triangle. "
                "It could be a case where triangle on front (part boundary of "
                "surface in interior) "
                "has three nodes front.");
      }
      side_ents3d.insert(*tets.begin());
      }

  } while (side_ents3d_tris_on_surface.size()!=triangles.size());

  if(ents3d_with_prisms.size() == side_ents3d.size()) {
    SETERRQ(
      m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,
      "All tets on one side, no-interface"
    );
  }
  //other side ents
  Range other_side = subtract(ents3d_with_prisms,side_ents3d);
  //side nodes
  Range side_nodes;
  rval = moab.get_connectivity(side_ents3d.subset_by_type(MBTET), side_nodes,
                               true);
  CHKERRQ_MOAB(rval);
  //nodes on crack surface without front
  nodes_without_front = intersect(nodes_without_front,side_nodes);
  Range side_edges;
  rval = moab.get_adjacencies(side_ents3d.subset_by_type(MBTET), 1, false,
                              side_edges, moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  skin_edges_boundary = intersect(skin_edges_boundary,side_edges);
  //make child meshsets
  std::vector<EntityHandle> children;
  rval = moab.get_child_meshsets(sideset,children);  CHKERRQ_MOAB(rval);
  if(children.empty()) {
    children.resize(3);
    rval = moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, children[0]);
    CHKERRQ_MOAB(rval);
    rval = moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, children[1]);
    CHKERRQ_MOAB(rval);
    rval = moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, children[2]);
    CHKERRQ_MOAB(rval);
    rval = moab.add_child_meshset(sideset,children[0]); CHKERRQ_MOAB(rval);
    rval = moab.add_child_meshset(sideset,children[1]); CHKERRQ_MOAB(rval);
    rval = moab.add_child_meshset(sideset,children[2]); CHKERRQ_MOAB(rval);
  } else {
    if(children.size()!=3) {
      SETERRQ(PETSC_COMM_SELF,1,"this meshset should have 3 children meshsets");
    }
    children.resize(3);
    ierr = moab.clear_meshset(&children[0],3); CHKERRQ(ierr);
  }
  EntityHandle &child_side = children[0];
  EntityHandle &child_other_side = children[1];
  EntityHandle &child_nodes_and_skin_edges = children[2];
  rval = moab.add_entities(child_side,side_ents3d); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(child_other_side,other_side); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(child_nodes_and_skin_edges, nodes_without_front);
  CHKERRQ_MOAB(rval);
  rval = moab.add_entities(child_nodes_and_skin_edges, skin_edges_boundary);
  CHKERRQ_MOAB(rval);
  if(verb>=VERBOSE) {
    PetscPrintf(m_field.get_comm(), "Nb. of side ents3d in set %u\n",
                side_ents3d.size());
    PetscPrintf(m_field.get_comm(), "Nb. of other side ents3d in set %u\n",
                other_side.size());
    PetscPrintf(m_field.get_comm(), "Nb. of boundary edges %u\n",
                skin_edges_boundary.size());
  }
  if(verb>=NOISY) {
    ierr = moab.write_file("side.vtk", "VTK", "", &children[0], 1);
    CHKERRQ(ierr);
    ierr = moab.write_file("other_side.vtk", "VTK", "", &children[1], 1);
    CHKERRQ(ierr);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode PrismInterface::findIfTringleHasThreeNodesOnInternalSurfaceSkin(
    const EntityHandle sideset, const BitRefLevel mesh_bit_level,
    const bool recursive, Range &faces_with_three_nodes_on_front, int verb) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;

  Range mesh_level_ents3d;
  Range mesh_level_edges,mesh_level_tris;
  if(mesh_bit_level.any()) {
    ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBTET, mesh_level_ents3d);
    CHKERRQ(ierr);
    ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBTRI, mesh_level_tris);
    CHKERRQ(ierr);
    ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        mesh_bit_level, BitRefLevel().set(), MBEDGE, mesh_level_edges);
    CHKERRQ(ierr);
  }

  Skinner skin(&moab);
  //get interface triangles from side set
  Range triangles;
  rval = moab.get_entities_by_type(sideset, MBTRI, triangles, recursive);
  CHKERRQ_MOAB(rval);
  if(mesh_bit_level.any()) {
    triangles = intersect(triangles,mesh_level_tris);
  }
  if(verb>=VERBOSE) {
    PetscPrintf(m_field.get_comm(),"Nb. of triangles in set %u\n",triangles.size());
  }
  //get nodes, edges and 3d ents (i.e. tets and prisms)
  Range nodes; // nodes from triangles
  rval = moab.get_connectivity(triangles, nodes, true);
  CHKERRQ_MOAB(rval);
  Range ents3d; // 3d ents form nodes
  rval = moab.get_adjacencies(nodes, 3, false, ents3d, moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  if(mesh_bit_level.any()) {
    ents3d = intersect(ents3d,mesh_level_ents3d);
  }
  //take skin faces
  Range skin_faces; // skin faces from 3d ents
  rval = skin.find_skin(0, ents3d.subset_by_type(MBTET), false, skin_faces);
  CHKERRQ_MOAB(rval);
  //take skin edges (boundary of surface if there is any)
  Range skin_edges_boundary; //skin edges from triangles
  rval = skin.find_skin(0, triangles, false, skin_edges_boundary);
  CHKERRQ_MOAB(rval);
  if(verb>=VERY_VERBOSE) {
    PetscPrintf(m_field.get_comm(), "skin_edges_boundary %u\n",
                skin_edges_boundary.size());
  }
  //take all edges on skin faces (i.e. skin surface)
  Range skin_faces_edges; //edges from skin faces of 3d ents
  rval = moab.get_adjacencies(skin_faces, 1, false, skin_faces_edges,
                              moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  if(mesh_bit_level.any()) {
    skin_faces_edges = intersect(skin_faces_edges,mesh_level_edges);
  }
  if(verb>=VERY_VERBOSE) {
    PetscPrintf(m_field.get_comm(),"skin_faces_edges %u\n",skin_faces_edges.size());
  }
  //note: that skin faces edges do not contain internal boundary
  //note: that prisms are not included in ents3d, so if ents3d have border with other inteface is like external boundary
  //skin edges boundary are internal edge <- skin_faces_edges contains edges which are on the body boundary <- that is the trick
  skin_edges_boundary = subtract(skin_edges_boundary,skin_faces_edges); // from skin edges subtract edges from skin faces of 3d ents (only internal edges)

  if(verb>=VERY_VERBOSE) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,out_meshset); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(out_meshset,triangles); CHKERRQ_MOAB(rval);
    rval = moab.write_file("triangles.vtk","VTK","",&out_meshset,1); CHKERRQ_MOAB(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERRQ_MOAB(rval);
    rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,out_meshset); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(out_meshset,ents3d); CHKERRQ_MOAB(rval);
    rval = moab.write_file("ents3d.vtk","VTK","",&out_meshset,1); CHKERRQ_MOAB(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERRQ_MOAB(rval);
    rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,out_meshset); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(out_meshset,skin_edges_boundary); CHKERRQ_MOAB(rval);
    rval = moab.write_file("skin_edges_boundary.vtk","VTK","",&out_meshset,1); CHKERRQ_MOAB(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERRQ_MOAB(rval);
  }
  if(verb>=VERY_VERBOSE) {
    PetscPrintf(
      m_field.get_comm(),"subtract skin_edges_boundary %u\n",skin_edges_boundary.size()
    );
  }

  //Get nodes on boundary edge
  Range skin_nodes_boundary;
  rval = moab.get_connectivity(skin_edges_boundary, skin_nodes_boundary, true);
  CHKERRQ_MOAB(rval);
  //Remove node which are boundary with other existing interface
  Range prisms_nodes;
  rval =
      moab.get_connectivity(ents3d.subset_by_type(MBPRISM), prisms_nodes, true);
  CHKERRQ_MOAB(rval);
  skin_nodes_boundary = subtract(skin_nodes_boundary,prisms_nodes);
  if(verb>=VERY_VERBOSE) {
    PetscPrintf(m_field.get_comm(), "subtract skin_nodes_boundary %u\n",
                skin_nodes_boundary.size());
  }
  //use nodes on body boundary and interface (without internal boundary) to find adjacent tets
  Range nodes_without_front = subtract(nodes,skin_nodes_boundary); // nodes_without_front adjacent to all split face edges except those on internal edge
  if(verb>=VERY_VERBOSE) {
    PetscPrintf(m_field.get_comm(),
                "adj. node if ents3d but not on the internal edge %u\n",
                nodes_without_front.size());
  }

  Range skin_nodes_boundary_tris;
  rval = moab.get_adjacencies(
    skin_nodes_boundary,2,false,skin_nodes_boundary_tris,moab::Interface::UNION
  ); CHKERRQ(ierr);
  Range skin_nodes_boundary_tris_nodes;
  rval = moab.get_connectivity(skin_nodes_boundary_tris,
                               skin_nodes_boundary_tris_nodes, true);
  CHKERRQ_MOAB(rval);
  skin_nodes_boundary_tris_nodes =
      subtract(skin_nodes_boundary_tris_nodes, skin_nodes_boundary);
  Range skin_nodes_boundary_tris_nodes_tris;
  rval = moab.get_adjacencies(skin_nodes_boundary_tris_nodes, 2, false,
                              skin_nodes_boundary_tris_nodes_tris,
                              moab::Interface::UNION);
  CHKERRQ(ierr);
  // Triangle which has tree nodes on front boundary
  skin_nodes_boundary_tris = intersect(
    triangles,
    subtract(skin_nodes_boundary_tris,skin_nodes_boundary_tris_nodes_tris)
  );
  faces_with_three_nodes_on_front.swap(skin_nodes_boundary_tris);

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode PrismInterface::splitSides(const EntityHandle meshset,
                                          const BitRefLevel &bit,
                                          const int msId,
                                          const CubitBCType cubit_bc_type,
                                          const bool add_interface_entities,
                                          const bool recursive, int verb) {

  Interface &m_field = cOre;
  MeshsetsManager *meshsets_manager_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.getInterface(meshsets_manager_ptr); CHKERRQ(ierr);
  CubitMeshSet_multiIndex::index<
    Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator miit =
      meshsets_manager_ptr->getMeshsetsMultindex()
          .get<Composite_Cubit_msId_And_MeshSetType_mi_tag>()
          .find(boost::make_tuple(msId, cubit_bc_type.to_ulong()));
  if (miit !=
      meshsets_manager_ptr->getMeshsetsMultindex()
          .get<Composite_Cubit_msId_And_MeshSetType_mi_tag>()
          .end()) {
    ierr = splitSides(meshset, bit, miit->meshset, add_interface_entities,
                      recursive, verb);
    CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"msId is not there");
  }
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode PrismInterface::splitSides(const EntityHandle meshset,
                                          const BitRefLevel &bit,
                                          const EntityHandle sideset,
                                          const bool add_interface_entities,
                                          const bool recursive, int verb) {

  MoFEMFunctionBeginHot;
  ierr = splitSides(meshset, bit, BitRefLevel(), BitRefLevel(), sideset,
                    add_interface_entities, recursive, verb);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode PrismInterface::splitSides(
    const EntityHandle meshset, const BitRefLevel &bit,
    const BitRefLevel &inhered_from_bit_level,
    const BitRefLevel &inhered_from_bit_level_mask, const EntityHandle sideset,
    const bool add_interface_entities, const bool recursive, int verb) {

  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  std::vector<EntityHandle> children;
  //get children meshsets
  rval = moab.get_child_meshsets(sideset,children);  CHKERRQ_MOAB(rval);
  if(children.size()!=3) {
    SETERRQ(
      m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,
      "should be 3 child meshsets, each of them contains tets on two sides of interface"
    );
  }
  //3d ents on "father" side
  Range side_ents3d;
  rval = moab.get_entities_by_handle(children[0], side_ents3d, false);
  CHKERRQ_MOAB(rval);
  //3d ents on "mather" side
  Range other_ents3d;
  rval = moab.get_entities_by_handle(children[1], other_ents3d, false);
  CHKERRQ_MOAB(rval);
  //faces of interface
  Range triangles;
  rval = moab.get_entities_by_type(sideset, MBTRI, triangles, recursive);
  CHKERRQ_MOAB(rval);
  Range side_ents3d_tris;
  rval = moab.get_adjacencies(side_ents3d, 2, true, side_ents3d_tris,
                              moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  triangles = intersect(triangles,side_ents3d_tris);
  //nodes on interface but not on crack front (those should not be splitted)
  Range nodes;
  rval = moab.get_entities_by_type(children[2], MBVERTEX, nodes, false);
  CHKERRQ_MOAB(rval);
  Range meshset_3d_ents,meshset_2d_ents;
  rval = moab.get_entities_by_dimension(meshset, 3, meshset_3d_ents, true);
  CHKERRQ_MOAB(rval);
  Range meshset_tets = meshset_3d_ents.subset_by_type(MBTET);
  rval = moab.get_adjacencies(meshset_tets, 2, false, meshset_2d_ents,
                              moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  side_ents3d = intersect(meshset_3d_ents,side_ents3d);
  other_ents3d = intersect(meshset_3d_ents,other_ents3d);
  triangles = intersect(meshset_2d_ents,triangles);
  if(verb>3) {
    PetscPrintf(m_field.get_comm(),"triangles %u\n",triangles.size());
    PetscPrintf(m_field.get_comm(),"side_ents3d %u\n",side_ents3d.size());
    PetscPrintf(m_field.get_comm(),"nodes %u\n",nodes.size());
  }

  const RefEntity_multiIndex *refined_ents_ptr;
  ierr = m_field.get_ref_ents(&refined_ents_ptr); CHKERRQ(ierr);
  typedef RefEntity_multiIndex::index<Ent_mi_tag>::type RefEntsByEntType;
  const RefEntsByEntType &ref_ents_by_type = refined_ents_ptr->get<Ent_mi_tag>();
  RefEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  //create view index by parent entity
  {
    typedef RefEntity_multiIndex::index<
        Composite_EntType_and_ParentEntType_mi_tag>::type RefEntsByComposite;
    const RefEntsByComposite &ref_ents =
    refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>();

    RefEntsByComposite::iterator miit;
    RefEntsByComposite::iterator hi_miit;
    //view by parent type (VERTEX)
    miit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX,MBVERTEX));
    hi_miit = ref_ents.upper_bound(boost::make_tuple(MBVERTEX,MBVERTEX));
    for(;miit!=hi_miit;miit++) {
      if(((*miit)->getBitRefLevel()&inhered_from_bit_level_mask) == (*miit)->getBitRefLevel()) {
        if(((*miit)->getBitRefLevel()&inhered_from_bit_level).any()) {
          std::pair<RefEntity_multiIndex_view_by_parent_entity::iterator,bool> p_ref_ent_view;
          p_ref_ent_view = ref_parent_ents_view.insert(*miit);
          if(!p_ref_ent_view.second) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"non unique insertion");
          }
        }
      }
    }
  }
  //maps nodes on "father" and "mather" side
  std::map<
    EntityHandle, /*node on "mather" side*/
    EntityHandle /*node on "father" side*/
  > map_nodes;
  //add new nodes on interface and create map
  Range::iterator nit = nodes.begin();
  double coord[3];
  for(;nit!=nodes.end();nit++) {
    //find ref enet
    RefEntsByEntType::iterator miit_ref_ent = ref_ents_by_type.find(*nit);
    if(miit_ref_ent == ref_ents_by_type.end()) {
      SETERRQ(PETSC_COMM_SELF,1,"can not find node in MoFEM database");
    }
    EntityHandle child_entity = 0;
    RefEntity_multiIndex::iterator child_it;
    RefEntity_multiIndex_view_by_parent_entity::iterator child_iit;
    child_iit = ref_parent_ents_view.find(*nit);
    if(child_iit != ref_parent_ents_view.end()) {
      child_it = refined_ents_ptr->find((*child_iit)->getRefEnt());
      BitRefLevel bit_child = (*child_it)->getBitRefLevel();
      if( (inhered_from_bit_level&bit_child).none() ) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      child_entity = (*child_it)->getRefEnt();
    }
    //
    bool success;
    if(child_entity == 0) {
      rval = moab.get_coords(&*nit,1,coord); CHKERRQ_MOAB(rval);
      EntityHandle new_node;
      rval = moab.create_vertex(coord,new_node); CHKERRQ_MOAB(rval);
      map_nodes[*nit] = new_node;
      //create new node on "father" side
      //parent is node on "mather" side
      rval = moab.tag_set_data(cOre.get_th_RefParentHandle(),&new_node,1,&*nit); CHKERRQ_MOAB(rval);
      std::pair<RefEntity_multiIndex::iterator,bool> p_ref_ent =
      const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->insert(
        boost::shared_ptr<RefEntity>
        (new RefEntity(m_field.get_basic_entity_data_ptr(),new_node))
      );
      //set ref bit level to node on "father" side
      success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->
      modify(p_ref_ent.first,RefEntity_change_add_bit(bit));
      if (!success)
        SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                "modification unsuccessful");
    } else {
      map_nodes[*nit] = child_entity;
      //set ref bit level to node on "father" side
      success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->
      modify(child_it,RefEntity_change_add_bit(bit));
      if (!success)
        SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                "modification unsuccessful");
    }
    //set ref bit level to node on "mather" side
    success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->
    modify(miit_ref_ent,RefEntity_change_add_bit(bit));
    if (!success)
      SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
              "modification unsuccessful");
  }
  //crete meshset for new mesh bit level
  EntityHandle meshset_for_bit_level;
  rval = moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER,
                             meshset_for_bit_level);
  CHKERRQ_MOAB(rval);
  //subtract those elements which will be refined, i.e. disconnected form other side elements, and connected to new prisms, if they area created
  meshset_3d_ents = subtract(meshset_3d_ents,side_ents3d);
  rval = moab.add_entities(meshset_for_bit_level, meshset_3d_ents);
  CHKERRQ_MOAB(rval);
  for(int dd = 0;dd<3;dd++) {
    Range ents_dd;
    rval = moab.get_adjacencies(meshset_3d_ents, dd, false, ents_dd,
                                moab::Interface::UNION);
    CHKERRQ_MOAB(rval);
    rval = moab.add_entities(meshset_for_bit_level, ents_dd);
    CHKERRQ_MOAB(rval);
  }
  //
  //typedef RefEntity_multiIndex::index<Composite_ParentEnt_And_EntType_mi_tag>::type ref_ent_by_composite;
  //ref_ent_by_composite &by_composite = refinedEntities.get<Composite_ParentEnt_And_EntType_mi_tag>();
  //create new 3d ents on "father" side
  Range new_3d_ents;
  for (Range::iterator eit3d = side_ents3d.begin(); eit3d != side_ents3d.end();
       eit3d++) {
    RefEntsByEntType::iterator miit_ref_ent = ref_ents_by_type.find(*eit3d);
    if(miit_ref_ent==ref_ents_by_type.end()) {
      SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
              "tet not in database");
    }
    int num_nodes;
    const EntityHandle* conn;
    rval = moab.get_connectivity(*eit3d, conn, num_nodes, true);
    CHKERRQ_MOAB(rval);
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
      rval = moab.add_entities(meshset_for_bit_level, &*eit3d, 1);
      rval = moab.add_entities(meshset_for_bit_level, conn, num_nodes);
      continue;
    }

    const RefElement_multiIndex *refined_finite_elements_ptr;
    ierr = m_field.get_ref_finite_elements(&refined_finite_elements_ptr);
    CHKERRQ(ierr);

    //here is created new or prism is on interface
    EntityHandle existing_ent = 0;
    /* check if tet element with new connectivity is in database*/
    RefElement_multiIndex::index<Ent_Ent_mi_tag>::type::iterator child_iit,hi_child_iit;
    child_iit =
        refined_finite_elements_ptr->get<Ent_Ent_mi_tag>().lower_bound(*eit3d);
    hi_child_iit =
        refined_finite_elements_ptr->get<Ent_Ent_mi_tag>().upper_bound(*eit3d);
    for(;child_iit!=hi_child_iit;child_iit++) {
      const EntityHandle* conn_ref_tet;
      rval = moab.get_connectivity(child_iit->getRefEnt(), conn_ref_tet,
                                   num_nodes, true);
      CHKERRQ_MOAB(rval);
      int nn = 0;
      for(;nn<num_nodes;nn++) {
        if(conn_ref_tet[nn]!=new_conn[nn]) {
          break;
        }
      }
      if(nn == num_nodes) {
        if(existing_ent != 0) {
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "database inconsistency");
        }
        existing_ent = child_iit->getRefEnt();
      }
    }
    switch (moab.type_from_handle(*eit3d)) {
      case MBTET: {
        RefEntsByEntType::iterator child_it;
        EntityHandle tet;
        if(existing_ent == 0) {
          Range new_conn_tet;
          rval = moab.get_adjacencies(new_conn, 4, 3, false, new_conn_tet);
          CHKERRQ_MOAB(rval);
          if(new_conn_tet.empty()) {
            rval = moab.create_element(MBTET, new_conn, 4, tet);
            CHKERRQ_MOAB(rval);
            rval = moab.tag_set_data(cOre.get_th_RefParentHandle(), &tet, 1,
                                     &*eit3d);
            CHKERRQ_MOAB(rval);
          } else {
            RefElement_multiIndex::index<Ent_mi_tag>::type::iterator rit,
                new_rit;
            rit = refined_finite_elements_ptr->get<Ent_mi_tag>().find(*eit3d);
            if(rit==refined_finite_elements_ptr->get<Ent_mi_tag>().end()) {
              SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                      "can't find this in database");
            }
            new_rit = refined_finite_elements_ptr->get<Ent_mi_tag>().find(
                *new_conn_tet.begin());
            if(new_rit==refined_finite_elements_ptr->get<Ent_mi_tag>().end()) {
              SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                      "can't find this in database");
            }
            tet = *new_conn_tet.begin();
            /*std::ostringstream ss;
            ss << "nb new conns: " << nb_new_conn << std::endl;
            ss << "new_conn_tets.size() " << new_conn_tet.size() << std::endl;
            ss << "data inconsistency\n";
            ss << "this ent:\n";
            ss << *rit->ref_ptr << std::endl;
            ss << "found this ent:\n";
            ss << *new_rit->ref_ptr << std::endl;
            SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());*/
          }
        } else {
          tet = existing_ent;
        }
        rval = moab.add_entities(meshset_for_bit_level, &tet, 1);
        CHKERRQ_MOAB(rval);
        rval = moab.add_entities(meshset_for_bit_level, new_conn, 4);
        CHKERRQ_MOAB(rval);
        new_3d_ents.insert(tet);
      } break;
      case MBPRISM: {
        EntityHandle prism;
        if(verb>3) {
          PetscPrintf(m_field.get_comm(), "prims nb_new_nodes %d\n",
                      nb_new_conn);
        }
        if(existing_ent == 0) {
          Range new_conn_prism;
          rval = moab.get_adjacencies(new_conn, 6, 3, false, new_conn_prism);
          CHKERRQ_MOAB(rval);
          if(new_conn_prism.empty()) {
            rval = moab.create_element(MBPRISM, new_conn, 6, prism);
            CHKERRQ_MOAB(rval);
            rval = moab.tag_set_data(cOre.get_th_RefParentHandle(), &prism, 1,
                                     &*eit3d);
            CHKERRQ_MOAB(rval);
          } else {
            SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency");
          }
        } else {
          prism = existing_ent;
        }
        rval = moab.add_entities(meshset_for_bit_level, &prism, 1);
        CHKERRQ_MOAB(rval);
        rval = moab.add_entities(meshset_for_bit_level, new_conn, 4);
        CHKERRQ_MOAB(rval);
        new_3d_ents.insert(prism);
      } break;
      default:
      SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"not implemented");
    }
  }
  Range new_ents;
  //create new entities by adjacencies form new tets
  rval = moab.get_adjacencies(new_3d_ents.subset_by_type(MBTET), 2, true,
                              new_ents, moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  rval = moab.get_adjacencies(new_3d_ents.subset_by_type(MBTET), 1, true,
                              new_ents, moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  //Tags for setting side
  Tag th_interface_side;
  const int def_side[] = {0};
  rval = moab.tag_get_handle("INTERFACE_SIDE", 1, MB_TYPE_INTEGER,
                             th_interface_side, MB_TAG_CREAT | MB_TAG_SPARSE,
                             def_side);
  CHKERRQ_MOAB(rval);
  //add new edges and triangles to mofem database
  Range ents;
  rval =
      moab.get_adjacencies(triangles, 1, false, ents, moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  ents.insert(triangles.begin(),triangles.end());
  Range new_ents_in_database; //this range contains all new entities
  for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
    int num_nodes;
    const EntityHandle* conn;
    rval = moab.get_connectivity(*eit, conn, num_nodes, true);
    CHKERRQ_MOAB(rval);
    EntityHandle new_conn[num_nodes];
    int nb_new_conn = 0;
    int ii = 0;
    for(;ii<num_nodes; ii++) {
      std::map<EntityHandle, EntityHandle>::iterator mit =
          map_nodes.find(conn[ii]);
      if(mit != map_nodes.end()) {
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
    if(nb_new_conn==0) continue;
    RefEntsByEntType::iterator miit_ref_ent = ref_ents_by_type.find(*eit);
    if(miit_ref_ent == ref_ents_by_type.end()) {
      SETERRQ(PETSC_COMM_SELF, 1,
              "this entity (edge or tri) should be already in database");
    }
    Range new_ent; //contains all entities (edges or triangles) added to mofem database
    switch (moab.type_from_handle(*eit)) {
      case MBTRI: {
        //get entity based on its connectivity
        rval = moab.get_adjacencies(new_conn, 3, 2, false, new_ent);
        CHKERRQ_MOAB(rval);
        if (new_ent.size() != 1)
          SETERRQ(PETSC_COMM_SELF, 1, "this tri should be in moab database");
        int new_side = 1;
        rval = moab.tag_set_data(th_interface_side, &*new_ent.begin(), 1,
                                 &new_side);
        CHKERRQ_MOAB(rval);
        if (verb >= VERY_VERBOSE)
          PetscPrintf(m_field.get_comm(), "new_ent %u\n", new_ent.size());
        //add prism element
        if(add_interface_entities) {
          if(inhered_from_bit_level.any()) {
            SETERRQ(PETSC_COMM_SELF, 1,
                    "not implemented for inhered_from_bit_level");
          }
          //set prism connectivity
          EntityHandle prism_conn[6] = {
            conn[0],conn[1],conn[2],
            new_conn[0],new_conn[1],new_conn[2]
          };
          EntityHandle prism;
          rval = moab.create_element(MBPRISM, prism_conn, 6, prism);
          CHKERRQ_MOAB(rval);
          ierr = cOre.addPrismToDatabase(prism, verb);
          CHKERRQ(ierr);
          rval = moab.add_entities(meshset_for_bit_level, &prism, 1);
          CHKERRQ_MOAB(rval);
        }
      } break;
      case MBEDGE: {
        rval = moab.get_adjacencies(new_conn, 2, 1, false, new_ent);
        CHKERRQ_MOAB(rval);
        if (new_ent.size() != 1) {
          if (m_field.get_comm_rank() == 0) {
            EntityHandle out_meshset;
            rval = moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER,
                                       out_meshset);
            CHKERRQ_MOAB(rval);
            rval = moab.add_entities(out_meshset, &*eit, 1);
            CHKERRQ_MOAB(rval);
            rval = moab.write_file("debug_splitSides.vtk", "VTK", "",
                                   &out_meshset, 1);
            CHKERRQ_MOAB(rval);
            rval = moab.delete_entities(&out_meshset, 1);
            CHKERRQ_MOAB(rval);
            rval = moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER,
                                       out_meshset);
            CHKERRQ_MOAB(rval);
            rval = moab.add_entities(out_meshset, side_ents3d);
            CHKERRQ_MOAB(rval);
            rval = moab.write_file("debug_splitSides_side_ents3d.vtk", "VTK",
                                   "", &out_meshset, 1);
            CHKERRQ_MOAB(rval);
            rval = moab.delete_entities(&out_meshset, 1);
            CHKERRQ_MOAB(rval);
            rval = moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER,
                                       out_meshset);
            CHKERRQ_MOAB(rval);
            rval = moab.add_entities(out_meshset, other_ents3d);
            CHKERRQ_MOAB(rval);
            rval = moab.write_file("debug_splitSides_other_ents3d.vtk", "VTK",
                                   "", &out_meshset, 1);
            CHKERRQ_MOAB(rval);
            rval = moab.delete_entities(&out_meshset, 1);
            CHKERRQ_MOAB(rval);
            rval = moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER,
                                       out_meshset);
            CHKERRQ_MOAB(rval);
            rval = moab.add_entities(out_meshset, triangles);
            CHKERRQ_MOAB(rval);
            rval = moab.write_file("debug_get_msId_3dENTS_split_triangles.vtk",
                                   "VTK", "", &out_meshset, 1);
            CHKERRQ_MOAB(rval);
            rval = moab.delete_entities(&out_meshset, 1);
            CHKERRQ_MOAB(rval);
          }
          SETERRQ2(PETSC_COMM_SELF, 1, "this edge should be in moab database "
                                       "new_ent.size() = %u nb_new_conn = %d",
                   new_ent.size(), nb_new_conn);
        }
      } break;
      default:
        SETERRQ(PETSC_COMM_SELF, 1, "Houston we have a problem !!!");
      }
      if (new_ent.size() != 1) {
        SETERRQ1(PETSC_COMM_SELF, 1,
                 "new_ent.size() = %u, size always should be 1",
                 new_ent.size());
    }
    //set parent
    rval = moab.tag_set_data(cOre.get_th_RefParentHandle(), &*new_ent.begin(),
                             1, &*eit);
    CHKERRQ_MOAB(rval);
    // add to database
    std::pair<RefEntity_multiIndex::iterator, bool> p_ref_ent =
        const_cast<RefEntity_multiIndex *>(refined_ents_ptr)
            ->insert(boost::shared_ptr<RefEntity>(new RefEntity(
                m_field.get_basic_entity_data_ptr(), new_ent[0])));
    const_cast<RefEntity_multiIndex *>(refined_ents_ptr)
        ->modify(p_ref_ent.first, RefEntity_change_add_bit(bit));
    new_ents_in_database.insert(new_ent.begin(), new_ent.end());
  }
  //all other entities, some ents like triangles and faces on the side of tets
  Range side_adj_faces_and_edges;
  rval = moab.get_adjacencies(side_ents3d.subset_by_type(MBTET), 1, true,
                              side_adj_faces_and_edges, moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  rval = moab.get_adjacencies(side_ents3d.subset_by_type(MBTET), 2, true,
                              side_adj_faces_and_edges, moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  //subtract entities already added to mofem database
  side_adj_faces_and_edges =
      subtract(side_adj_faces_and_edges, new_ents_in_database);
  for (Range::iterator eit = side_adj_faces_and_edges.begin();
       eit != side_adj_faces_and_edges.end(); eit++) {
    int num_nodes;
    const EntityHandle* conn;
    rval = moab.get_connectivity(*eit, conn, num_nodes, true);
    CHKERRQ_MOAB(rval);
    EntityHandle new_conn[num_nodes];
    int nb_new_conn = 0;
    int ii = 0;
    for (; ii < num_nodes; ii++) {
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
    RefEntsByEntType::iterator miit_ref_ent = ref_ents_by_type.find(*eit);
    if (miit_ref_ent == ref_ents_by_type.end()) {
      SETERRQ1(PETSC_COMM_SELF, 1,
               "entity should be in MoFem database, num_nodes = %d", num_nodes);
    }
    Range new_ent;
    switch (moab.type_from_handle(*eit)) {
    case MBTRI: {
      rval = moab.get_adjacencies(new_conn, 3, 2, false, new_ent);
      CHKERRQ_MOAB(rval);
    } break;
    case MBEDGE: {
      rval = moab.get_adjacencies(new_conn, 2, 1, false, new_ent);
      CHKERRQ_MOAB(rval);
    } break;
    default:
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "Houston we have a problem");
    }
    if (new_ent.size() != 1) {
      SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "database inconsistency, new_ent.size() = %u", new_ent.size());
    }
    // add entity to mofem database
    rval = moab.tag_set_data(cOre.get_th_RefParentHandle(), &*new_ent.begin(),
                             1, &*eit);
    CHKERRQ_MOAB(rval);
    std::pair<RefEntity_multiIndex::iterator, bool> p_ref_ent =
        const_cast<RefEntity_multiIndex *>(refined_ents_ptr)
            ->insert(boost::shared_ptr<RefEntity>(new RefEntity(
                m_field.get_basic_entity_data_ptr(), new_ent[0])));
    const_cast<RefEntity_multiIndex *>(refined_ents_ptr)
        ->modify(p_ref_ent.first, RefEntity_change_add_bit(bit));
    if (verb > VERY_VERBOSE)
      PetscPrintf(m_field.get_comm(), "new_ent %u\n", new_ent.size());
    new_ents_in_database.insert(new_ent.begin(), new_ent.end());
  }
  // add new prisms which parents are part of other intefaces
  Range new_3d_prims = new_3d_ents.subset_by_type(MBPRISM);
  for (Range::iterator pit = new_3d_prims.begin(); pit != new_3d_prims.end();
       pit++) {
    ierr = cOre.addPrismToDatabase(*pit, verb);
    CHKERRQ(ierr);
    // get parent entity
    EntityHandle parent_prism;
    rval = moab.tag_get_data(cOre.get_th_RefParentHandle(), &*pit, 1,
                             &parent_prism);
    CHKERRQ_MOAB(rval);
    const EntityHandle root_meshset = moab.get_root_set();
    if (parent_prism == root_meshset) {
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "this prism should have parent");
    }
    if (moab.type_from_handle(parent_prism) != MBPRISM) {
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "this prism should have parent which is prism as well");
    }
    int num_nodes;
    // parent prism
    const EntityHandle *conn_parent;
    rval = moab.get_connectivity(parent_prism, conn_parent, num_nodes, true);
    MOAB_THROW(rval);
    Range face_side3_parent, face_side4_parent;
    rval = moab.get_adjacencies(conn_parent, 3, 2, false, face_side3_parent);
    CHKERRQ_MOAB(rval);
    rval =
        moab.get_adjacencies(&conn_parent[3], 3, 2, false, face_side4_parent);
    CHKERRQ_MOAB(rval);
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
    rval = moab.get_connectivity(*pit, conn, num_nodes, true);
    MOAB_THROW(rval);
    Range face_side3, face_side4;
    rval = moab.get_adjacencies(conn, 3, 2, false, face_side3);
    CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(&conn[3], 3, 2, false, face_side4);
    CHKERRQ_MOAB(rval);
    if (face_side3.size() != 1) {
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY, "face3 is missing");
    }
    if (face_side4.size() != 1) {
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY, "face4 is missing");
    }
    //
    std::vector<EntityHandle> face(2), parent_face(2);
    face[0] = *face_side3.begin();
    face[1] = *face_side4.begin();
    parent_face[0] = *face_side3_parent.begin();
    parent_face[1] = *face_side4_parent.begin();
    for (int ff = 0; ff < 2; ff++) {
      if (parent_face[ff] == face[ff])
        continue;
      int interface_side;
      rval = moab.tag_get_data(th_interface_side, &parent_face[ff], 1,
                               &interface_side);
      CHKERRQ_MOAB(rval);
      rval =
          moab.tag_set_data(th_interface_side, &face[ff], 1, &interface_side);
      CHKERRQ_MOAB(rval);
      EntityHandle parent_tri;
      rval = moab.tag_get_data(cOre.get_th_RefParentHandle(), &face[ff], 1,
                               &parent_tri);
      CHKERRQ_MOAB(rval);
      if (parent_tri != parent_face[ff]) {
        SETERRQ1(PETSC_COMM_SELF, 1, "wrong parent %lu", parent_tri);
      }
      if (new_ents_in_database.find(face[ff]) == new_ents_in_database.end()) {
        RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent =
            refined_ents_ptr->get<Ent_mi_tag>().find(face[ff]);
        if (miit_ref_ent == refined_ents_ptr->get<Ent_mi_tag>().end()) {
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "this is not in database, but should not be");
        }
      }
    }
  }
  // finalise by adding new tets and prism ti bit level
  ierr = m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
      meshset_for_bit_level, 3, bit);
  CHKERRQ(ierr);
  rval = moab.delete_entities(&meshset_for_bit_level, 1);
  CHKERRQ_MOAB(rval);
  ierr = moab.clear_meshset(&children[0], 3);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}
}
