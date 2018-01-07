/** \file CutMeshInterface.cpp
 * \brief Cut mesh by surface
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
CutMeshInterface::query_interface(const MOFEMuuid &uuid,
                                  UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMCutMesh) {
    *iface = const_cast<CutMeshInterface *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

CutMeshInterface::CutMeshInterface(const Core &core)
    : cOre(const_cast<Core &>(core)) {
  lineSearchSteps = 10;
  nbMaxMergingCycles = 200;
  nbMaxTrimSearchIterations = 20;
}

MoFEMErrorCode CutMeshInterface::setSurface(const Range &surface) {
  MoFEMFunctionBeginHot;
  sUrface = surface;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CutMeshInterface::copySurface(const Range &surface, Tag th,
                                             double *shift, double *origin,
                                             double *transform,
                                             const std::string save_mesh) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  std::map<EntityHandle,EntityHandle> verts_map;
  for (Range::const_iterator tit = surface.begin(); tit != surface.end();
       tit++) {
    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(*tit, conn, num_nodes, true);
    MatrixDouble coords(num_nodes, 3);
    if (th) {
      CHKERR moab.tag_get_data(th, conn, num_nodes, &coords(0, 0));
    } else {
      CHKERR moab.get_coords(conn, num_nodes, &coords(0, 0));
    }
    EntityHandle new_verts[num_nodes];
    for (int nn = 0; nn != num_nodes; nn++) {
      if(verts_map.find(conn[nn])!=verts_map.end()) {
        new_verts[nn] = verts_map[conn[nn]];
      } else {
        if (transform) {
          ublas::matrix_row<MatrixDouble> mr(coords, nn);
          if (origin) {
            VectorAdaptor vec_origin(
                3, ublas::shallow_array_adaptor<double>(3, origin));
            mr = mr - vec_origin;
          }
          MatrixAdaptor mat_transform = MatrixAdaptor(
              3, 3, ublas::shallow_array_adaptor<double>(9, transform));
          mr = prod(mat_transform, mr);
          if (origin) {
            VectorAdaptor vec_origin(
                3, ublas::shallow_array_adaptor<double>(3, origin));
            mr = mr + vec_origin;
          }
        }
        if (shift) {
          ublas::matrix_row<MatrixDouble> mr(coords, nn);
          VectorAdaptor vec_shift(
              3, ublas::shallow_array_adaptor<double>(3, shift));
          mr = mr + vec_shift;
        }
        CHKERR moab.create_vertex(&coords(nn, 0), new_verts[nn]);
        verts_map[conn[nn]] = new_verts[nn];
      }
    }
    EntityHandle ele;
    CHKERR moab.create_element(MBTRI, new_verts, num_nodes, ele);
    sUrface.insert(ele);
  }
  if (!save_mesh.empty()) {
    EntityHandle meshset;
    CHKERR m_field.get_moab().create_meshset(MESHSET_SET, meshset);
    CHKERR m_field.get_moab().add_entities(meshset, sUrface);
    CHKERR m_field.get_moab().write_file(save_mesh.c_str(), "VTK", "", &meshset,
                                         1);
    CHKERR m_field.get_moab().delete_entities(&meshset, 1);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::setVolume(const Range &volume) {
  MoFEMFunctionBeginHot;
  vOlume = volume;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CutMeshInterface::mergeSurface(const Range &surface) {
  MoFEMFunctionBeginHot;
  sUrface.merge(surface);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CutMeshInterface::mergeVolumes(const Range &volume) {
  MoFEMFunctionBeginHot;
  vOlume.merge(volume);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CutMeshInterface::buildTree() {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  treeSurfPtr = boost::shared_ptr<OrientedBoxTreeTool>(
      new OrientedBoxTreeTool(&moab, "ROOTSETSURF", true));
  CHKERR treeSurfPtr->build(sUrface, rootSetSurf);
  MoFEMFunctionReturn(0);
}

struct UpdateMeshsets {
  MoFEMErrorCode operator()(Core &core, const BitRefLevel &bit) const {
    MoFEMFunctionBeginHot;
    ierr = core.getInterface<MeshsetsManager>()
               ->updateAllMeshsetsByEntitiesChildren(bit);
    CHKERRG(ierr);

    MoFEMFunctionReturnHot(0);
  }
};

MoFEMErrorCode CutMeshInterface::cutAndTrim(
    const BitRefLevel &bit_level1, const BitRefLevel &bit_level2, Tag th,
    const double tol_cut, const double tol_cut_close, const double tol_trim,
    const double tol_trim_close, Range *fixed_edges, Range *corner_nodes,
    const bool update_meshsets,const bool debug) {
  CoreInterface &m_field = cOre;
  MoFEMFunctionBegin;
  // cut mesh
  CHKERR findEdgesToCut(tol_cut,QUIET,debug);
  CHKERR projectZeroDistanceEnts(tol_cut_close);
  CHKERR cutEdgesInMiddle(bit_level1, debug);
  if (fixed_edges) {
    CHKERR cOre.getInterface<BitRefManager>()->updateRange(*fixed_edges,
                                                           *fixed_edges);
  }
  if (corner_nodes) {
    CHKERR cOre.getInterface<BitRefManager>()->updateRange(*corner_nodes,
                                                           *corner_nodes);
  }
  if (update_meshsets) {
    CHKERR UpdateMeshsets()(cOre, bit_level1);
  }
  CHKERR moveMidNodesOnCutEdges(th);

  if(debug) {
    CHKERR saveCutEdges();
  }

  // trim mesh
  CHKERR findEdgesToTrim(th, tol_trim);
  CHKERR trimEdgesInTheMiddle(bit_level2, th, tol_trim_close, debug);
  if (fixed_edges) {
    CHKERR cOre.getInterface<BitRefManager>()->updateRange(*fixed_edges,
                                                           *fixed_edges);
  }
  if (corner_nodes) {
    CHKERR cOre.getInterface<BitRefManager>()->updateRange(*corner_nodes,
                                                           *corner_nodes);
  }
  if (update_meshsets) {
    CHKERR UpdateMeshsets()(cOre, bit_level2);
  }
  CHKERR moveMidNodesOnTrimmedEdges(th);

  if(debug) {
    CHKERR cOre.getInterface<BitRefManager>()->writeBitLevelByType(
        bit_level2, BitRefLevel().set(), MBTET, "out_tets_trim.vtk", "VTK", "");
    {
      EntityHandle meshset;
      CHKERR m_field.get_moab().create_meshset(MESHSET_SET, meshset);
      CHKERR m_field.get_moab().add_entities(meshset, trimNewSurfaces);
      CHKERR m_field.get_moab().write_file("trim_new_surfaces.vtk", "VTK", "",
                                           &meshset, 1);
      CHKERR m_field.get_moab().delete_entities(&meshset, 1);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::cutTrimAndMerge(
    const int fraction_level, const BitRefLevel &bit_level1,
    const BitRefLevel &bit_level2, const BitRefLevel &bit_level3, Tag th,
    const double tol_cut, const double tol_cut_close, const double tol_trim,
    const double tol_trim_close, Range &fixed_edges, Range &corner_nodes,
    const bool update_meshsets, const bool debug) {
  MoFEMFunctionBegin;
  if(debug) {
    CHKERR cOre.getInterface<BitRefManager>()->writeEntitiesNotInDatabase(
        "ents_not_in_database.vtk", "VTK", "");
  }
  CHKERR cutAndTrim(bit_level1, bit_level2, th, tol_cut, tol_cut_close,
                    tol_trim, tol_trim_close, &fixed_edges, &corner_nodes,
                    update_meshsets, debug);
  if(debug) {
    CHKERR cOre.getInterface<BitRefManager>()->writeEntitiesNotInDatabase(
        "cut_trim_ents_not_in_database.vtk", "VTK", "");
  }

  CHKERR mergeBadEdges(fraction_level, bit_level2, bit_level1, bit_level3,
                       getNewTrimSurfaces(), fixed_edges, corner_nodes, th,
                       update_meshsets, debug);
  CHKERR removePathologicalFrontTris(bit_level3,
                                     const_cast<Range &>(getMergedSurfaces()));

  if(debug) {
    CHKERR cOre.getInterface<BitRefManager>()->writeBitLevelByType(
        bit_level3, BitRefLevel().set(), MBTET, "out_tets_merged.vtk", "VTK",
        "");
    CHKERR cOre.getInterface<BitRefManager>()->writeEntitiesNotInDatabase(
        "cut_trim_merge_ents_not_in_database.vtk", "VTK", "");
  }

  CHKERR
      cOre.getInterface<BitRefManager>()->updateRange(fixed_edges, fixed_edges);
  CHKERR cOre.getInterface<BitRefManager>()->updateRange(corner_nodes,
                                                         corner_nodes);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::findEdgesToCut(const double low_tol,
                                                int verb,const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  // Range vol_edges;
  // rval = moab.get_adjacencies(
  //   vOlume,1,true,vol_edges,moab::Interface::UNION
  // ); CHKERRG(rval);
  edgesToCut.clear();
  cutEdges.clear();
  double ray_length;
  double ray_point[3], unit_ray_dir[3];
  VectorAdaptor vec_unit_ray_dir(
      3, ublas::shallow_array_adaptor<double>(3, unit_ray_dir));
  VectorAdaptor vec_ray_point(
      3, ublas::shallow_array_adaptor<double>(3, ray_point));

  Tag th_dist;
  rval = moab.tag_get_handle("DIST",th_dist);
  if(rval == MB_SUCCESS) {
    CHKERR moab.tag_delete(th_dist);
  } else {
    rval = MB_SUCCESS;
  }
  Tag th_dist_normal;
  rval = moab.tag_get_handle("DIST_NORMAL", th_dist_normal); 
  if(rval == MB_SUCCESS) {
     CHKERR moab.tag_delete(th_dist_normal);
  } else {
    rval = MB_SUCCESS;
  }
  double def_val[] = {0};
  CHKERR moab.tag_get_handle("DIST", 1, MB_TYPE_DOUBLE, th_dist,
                             MB_TAG_CREAT | MB_TAG_SPARSE, def_val);
  CHKERR moab.tag_get_handle("DIST_NORMAL", 1, MB_TYPE_DOUBLE, th_dist_normal,
                             MB_TAG_CREAT | MB_TAG_SPARSE, def_val);
  Range vol_vertices;
  CHKERR moab.get_connectivity(vOlume, vol_vertices, true);
  for (Range::iterator vit = vol_vertices.begin(); vit != vol_vertices.end();
       vit++) {
    double coords[3];
    CHKERR moab.get_coords(&*vit, 1, coords);
    VectorAdaptor point_in(3, ublas::shallow_array_adaptor<double>(3, coords));
    double p_out[3];
    EntityHandle facets_out;
    CHKERR treeSurfPtr->closest_to_location(&coords[0], rootSetSurf, p_out,
                                            facets_out);
    VectorAdaptor point_out(3, ublas::shallow_array_adaptor<double>(3, p_out));
    double normal[3];
    Util::normal(&moab, facets_out, normal[0], normal[1], normal[2]);
    VectorAdaptor n(3, ublas::shallow_array_adaptor<double>(3, normal));
    VectorDouble3 delta = point_out - point_in;
    double dist = norm_2(delta);
    double dist_normal = inner_prod(delta, n) / norm_2(n);
    CHKERR moab.tag_set_data(th_dist, &*vit, 1, &dist);
    CHKERR moab.tag_set_data(th_dist_normal, &*vit, 1, &dist_normal);
  }

  Range vol_edges;
  CHKERR moab.get_adjacencies(vOlume, 1, true, vol_edges,
                              moab::Interface::UNION);
  aveLength = 0;
  maxLength = 0;
  int nb_ave_length = 0;
  for (Range::iterator eit = vol_edges.begin(); eit != vol_edges.end(); eit++) {
    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(*eit, conn, num_nodes, true);
    double dist[num_nodes];
    CHKERR moab.tag_get_data(th_dist, conn, num_nodes, dist);
    double dist_normal[num_nodes];
    rval = moab.tag_get_data(th_dist_normal, conn, num_nodes, dist_normal);
    CHKERR getRayForEdge(*eit, vec_ray_point, vec_unit_ray_dir, ray_length);
    const double tol = ray_length * low_tol;
    if ((dist_normal[0] * dist_normal[1] < 0) ||
        (dist_normal[0] * dist_normal[1] == 0 &&
         (dist_normal[0] + dist_normal[1]) > 0)) {
      std::vector<double> distances_out;
      std::vector<EntityHandle> facets_out;
      CHKERR treeSurfPtr->ray_intersect_triangles(distances_out, facets_out,
                                                  rootSetSurf, tol, ray_point,
                                                  unit_ray_dir, &ray_length);
      if (!distances_out.empty()) {
        aveLength += ray_length;
        maxLength = fmax(maxLength, ray_length);
        nb_ave_length++;
        edgesToCut[*eit].dIst = distances_out[0];
        edgesToCut[*eit].lEngth = ray_length;
        edgesToCut[*eit].unitRayDir = vec_unit_ray_dir;
        edgesToCut[*eit].rayPoint = vec_ray_point;
        cutEdges.insert(*eit);
      } else if(dist_normal[0] * dist_normal[1] != 0) {
        const double b = dist_normal[0];
        const double a = dist_normal[1] - b;
        const double s = (-b/a)*ray_length;
        double point[3] = {0,0,0};
        VectorAdaptor vec_point(3,
                                ublas::shallow_array_adaptor<double>(3, point));
        noalias(vec_point) = vec_ray_point + s * vec_unit_ray_dir;
        double p_out[3];
        EntityHandle facets_out;
        CHKERR treeSurfPtr->closest_to_location(&point[0], rootSetSurf, p_out,
                                                facets_out);
        VectorAdaptor point_out(3,
                                ublas::shallow_array_adaptor<double>(3, p_out));
        VectorDouble3 delta = point_out - vec_point;
        double dist = norm_2(delta);
        if (dist < tol) {
          aveLength += ray_length;
          maxLength = fmax(maxLength, ray_length);
          nb_ave_length++;
          edgesToCut[*eit].dIst = dist;
          edgesToCut[*eit].lEngth = ray_length;
          edgesToCut[*eit].unitRayDir = vec_unit_ray_dir;
          edgesToCut[*eit].rayPoint = vec_ray_point;
          cutEdges.insert(*eit);
        }
      }
    }
  }
  aveLength /= nb_ave_length;

  cutVolumes.clear();
  // take all volumes adjacent to cut edges
  CHKERR moab.get_adjacencies(cutEdges, 3, false, cutVolumes,
                              moab::Interface::UNION);
  cutVolumes = intersect(cutVolumes,vOlume);

  // get edges on the cut volumes
  Range edges;
  CHKERR moab.get_adjacencies(cutVolumes, 1, false, edges,
                              moab::Interface::UNION);
  edges = subtract(edges, cutEdges);

  // add to cut set edges which are cutted by extension of cutting surface
  for (Range::iterator eit = edges.begin(); eit != edges.end(); eit++) {
    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(*eit, conn, num_nodes, true);
    double dist[num_nodes];
    CHKERR moab.tag_get_data(th_dist, conn, num_nodes, dist);
    double dist_normal[num_nodes];
    CHKERR moab.tag_get_data(th_dist_normal, conn, num_nodes, dist_normal);
    if (dist_normal[0] * dist_normal[1] < 0 ||
        (dist_normal[0] * dist_normal[1] == 0 &&
         (dist_normal[0] + dist_normal[1]) > 0)) {
      CHKERR getRayForEdge(*eit, vec_ray_point, vec_unit_ray_dir, ray_length);
      double s =
          fabs(dist_normal[0]) / (fabs(dist_normal[0]) + fabs(dist_normal[1]));
      edgesToCut[*eit].dIst = s * ray_length;
      edgesToCut[*eit].lEngth = ray_length;
      edgesToCut[*eit].unitRayDir = vec_unit_ray_dir;
      edgesToCut[*eit].rayPoint = vec_ray_point;
      cutEdges.insert(*eit);
    }
  }

  // if(debug) {
  //   EntityHandle out_meshset;
  //   CHKERR m_field.get_moab().create_meshset(MESHSET_SET, out_meshset);
  //   CHKERR m_field.get_moab().add_entities(out_meshset, cutEdges);
  //   CHKERR m_field.get_moab().write_file("cut_edges.vtk", "VTK", "",
  //                                        &out_meshset, 1);
  //   CHKERR m_field.get_moab().delete_entities(&out_meshset, 1);
  // }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::projectZeroDistanceEnts(const double low_tol,
                                                     int verb) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  Skinner skin(&moab);
  Range tets_skin;
  rval = skin.find_skin(0, vOlume, false, tets_skin);
  Range vol_faces, vol_edges, vol_nodes;
  CHKERR moab.get_adjacencies(vOlume, 2, false, vol_faces,
                              moab::Interface::UNION);
  CHKERR moab.get_adjacencies(vOlume, 1, false, vol_edges,
                              moab::Interface::UNION);
  CHKERR moab.get_connectivity(vOlume, vol_nodes, true);
  Range tets_skin_edges;
  CHKERR moab.get_adjacencies(tets_skin, 1, false, tets_skin_edges,
                              moab::Interface::UNION);
  Range tets_skin_verts;
  CHKERR moab.get_connectivity(tets_skin, tets_skin_verts, true);
  Tag th_dist;
  CHKERR moab.tag_get_handle("DIST_NORMAL", th_dist);
  Range cut_edge_verts;
  CHKERR moab.get_connectivity(cutEdges, cut_edge_verts, false);
  // get faces and edges
  Range cut_edges_faces;
  CHKERR moab.get_adjacencies(cut_edge_verts, 1, true, cut_edges_faces,
                              moab::Interface::UNION);
  CHKERR moab.get_adjacencies(cut_edge_verts, 2, true, cut_edges_faces,
                              moab::Interface::UNION);
  cut_edges_faces = intersect(cut_edges_faces,unite(vol_edges,vol_faces));
  zeroDistanceEnts.clear();
  verticesOnCutEdges.clear();
  for (Range::iterator fit = cut_edges_faces.begin();
       fit != cut_edges_faces.end(); fit++) {
    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(*fit, conn, num_nodes, true);
    double dist[] = {0, 0, 0};
    CHKERR moab.tag_get_data(th_dist, conn, num_nodes, dist);
    if (fabs(dist[0]) < low_tol * aveLength &&
        fabs(dist[1]) < low_tol * aveLength &&
        fabs(dist[2]) < low_tol * aveLength) {
      zeroDistanceEnts.insert(*fit);
      Range adj_edges;
      CHKERR moab.get_adjacencies(conn, num_nodes, 1, false, adj_edges,
                                  moab::Interface::UNION);
      for (Range::iterator eit = adj_edges.begin(); eit != adj_edges.end();
           eit++) {
        cutEdges.erase(*eit);
        edgesToCut.erase(*eit);
      }
      double coords[9];
      CHKERR moab.get_coords(conn, num_nodes, coords);
      CHKERRG(ierr);
      for (int nn = 0; nn != num_nodes; nn++) {
        // check if entity is on the skin, if it nothing can be done here
        // is not part of the surface in the body
        if (tets_skin_verts.find(conn[nn]) == tets_skin_verts.end()) {
          VectorAdaptor s0(
              3, ublas::shallow_array_adaptor<double>(3, &coords[3 * nn]));
          double p_out[3];
          EntityHandle facets_out;
          CHKERR treeSurfPtr->closest_to_location(&s0[0], rootSetSurf, p_out,
                                                  facets_out);
          VectorAdaptor point_out(
              3, ublas::shallow_array_adaptor<double>(3, p_out));
          VectorDouble3 ray = point_out - s0;
          double dist0 = norm_2(ray);
          verticesOnCutEdges[conn[nn]].dIst = dist0;
          verticesOnCutEdges[conn[nn]].lEngth = dist0;
          verticesOnCutEdges[conn[nn]].unitRayDir =
              dist0 > 0 ? ray / dist0 : ray;
          verticesOnCutEdges[conn[nn]].rayPoint = s0;
        }
      }
    }
  }
  CHKERR moab.tag_get_handle("DIST", th_dist);
  for (Range::iterator vit = cut_edge_verts.begin();
       vit != cut_edge_verts.end(); vit++) {
    double dist[] = {0};
    CHKERR moab.tag_get_data(th_dist, &*vit, 1, dist);
    if (fabs(dist[0]) < low_tol * aveLength) {
      zeroDistanceVerts.insert(*vit);
      Range adj_edges;
      CHKERR moab.get_adjacencies(&*vit, 1, 1, false, adj_edges);
      adj_edges = intersect(adj_edges,vol_edges);
      for (Range::iterator eit = adj_edges.begin(); eit != adj_edges.end();
           eit++) {
        cutEdges.erase(*eit);
        edgesToCut.erase(*eit);
      }
      double coords[3];
      CHKERR moab.get_coords(&*vit, 1, coords);
      VectorAdaptor s0(3, ublas::shallow_array_adaptor<double>(3, &coords[0]));
      // Vertex is in the body, just project in on the surface
      if (tets_skin_verts.find(*vit) == tets_skin_verts.end()) {
        double p_out[3];
        EntityHandle facets_out;
        CHKERR treeSurfPtr->closest_to_location(&s0[0], rootSetSurf, p_out,
                                                facets_out);
        VectorAdaptor point_out(3,
                                ublas::shallow_array_adaptor<double>(3, p_out));
        VectorDouble3 ray                   = point_out - s0;
        double dist0                        = norm_2(ray);
        verticesOnCutEdges[*vit].dIst       = dist0;
        verticesOnCutEdges[*vit].lEngth     = dist0;
        verticesOnCutEdges[*vit].unitRayDir = dist0 > 0 ? ray / dist0 : ray;
        verticesOnCutEdges[*vit].rayPoint   = s0;
      } else {
        double ray_length;
        double ray_point[3], unit_ray_dir[3];
        VectorAdaptor vec_unit_ray_dir(
            3, ublas::shallow_array_adaptor<double>(3, unit_ray_dir));
        VectorAdaptor vec_ray_point(
            3, ublas::shallow_array_adaptor<double>(3, ray_point));
        // Vertex is on the skin
        adj_edges = intersect(adj_edges,tets_skin_edges);
        if(adj_edges.empty()) {
          SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
                  "Huston we have a problem, how vertex can be on surface but "
                  "not have adjacent surface edges to it?");
        }
        // Make a loop over all adjacent surface edges, to find direction of the
        // ray. 
        for (Range::iterator eit = adj_edges.begin(); eit != adj_edges.end();
             eit++) {
          CHKERR getRayForEdge(*eit, vec_ray_point, vec_unit_ray_dir,
                               ray_length);
          // FIXME: Pick one edge on the skin and send ray to move vertex. Edge
          // should be pick one from corner edges, now is arbitrary
          std::vector<double> distances_out;
          std::vector<EntityHandle> facets_out;
          CHKERR treeSurfPtr->ray_intersect_triangles(
              distances_out, facets_out, rootSetSurf, low_tol, ray_point,
              unit_ray_dir, &ray_length);
          if (!distances_out.empty()) {
            VectorDouble point_out(vec_ray_point +
                                   distances_out[0] * vec_unit_ray_dir);
            VectorDouble3 ray                   = point_out - s0;
            double dist0                        = norm_2(ray);
            verticesOnCutEdges[*vit].dIst       = dist0;
            verticesOnCutEdges[*vit].lEngth     = dist0;
            verticesOnCutEdges[*vit].unitRayDir = dist0 > 0 ? ray / dist0 : ray;
            verticesOnCutEdges[*vit].rayPoint   = s0;
            break;
          }
        }
      }
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::cutEdgesInMiddle(const BitRefLevel bit,
                                                  const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MeshRefinement *refiner;
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBegin;
  if(cutEdges.size() != edgesToCut.size()) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Data inconsistency");
  }
  CHKERR m_field.getInterface(refiner);
  CHKERR m_field.get_ref_ents(&ref_ents_ptr);
  CHKERR refiner->add_verices_in_the_middel_of_edges(cutEdges, bit);
  CHKERR refiner->refine_TET(vOlume, bit, false, QUIET,
                             debug ? &cutEdges : NULL);
  // Tag th_ray_dir;
  // double def_val[] = {0,0,0};
  // CHKERR moab.tag_get_handle(
  //   "RAY_DIR",3,MB_TYPE_DOUBLE,th_ray_dir,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  // ); 
  cutNewVolumes.clear();
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, BitRefLevel().set(), MBTET, cutNewVolumes);
  cutNewSurfaces.clear();
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, bit, MBTRI, cutNewSurfaces);
  // Find new vertices on cut edges
  cutNewVertices.clear();
  CHKERR moab.get_connectivity(zeroDistanceEnts, cutNewVertices, true);
  cutNewVertices.merge(zeroDistanceVerts);
  for (map<EntityHandle, TreeData>::iterator mit = edgesToCut.begin();
       mit != edgesToCut.end(); ++mit) {
    RefEntity_multiIndex::index<
        Composite_ParentEnt_And_EntType_mi_tag>::type::iterator vit =
        ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>().find(
            boost::make_tuple(MBVERTEX, mit->first));
    if (vit ==
        ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>().end()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "No vertex on cut edges, that make no sense");
    }
    const boost::shared_ptr<RefEntity> &ref_ent = *vit;
    if ((ref_ent->getBitRefLevel() & bit).any()) {
      EntityHandle vert = ref_ent->getRefEnt();
      cutNewVertices.insert(vert);
      verticesOnCutEdges[vert] = mit->second;
    } else {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Vertex has wrong bit ref level");
    }
  }
  // Add zero distance entities faces
  Range tets_skin;
  Skinner skin(&moab);
  CHKERR skin.find_skin(0, cutNewVolumes, false, tets_skin);
  cutNewSurfaces.merge(
      subtract(zeroDistanceEnts.subset_by_type(MBTRI), tets_skin));
  // At that point cutNewSurfaces has all newly created faces, now take all
  // nodes on those faces and subtract nodes on catted edges. Faces adjacent to
  // nodes which left are not part of surface.
  Range diff_verts;
  CHKERR moab.get_connectivity(cutNewSurfaces, diff_verts, true);
  diff_verts = subtract(diff_verts, cutNewVertices);
  Range subtract_faces;
  CHKERR moab.get_adjacencies(diff_verts, 2, false, subtract_faces,
                              moab::Interface::UNION);
  cutNewSurfaces = subtract(cutNewSurfaces, subtract_faces);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::moveMidNodesOnCutEdges(Tag th) {
  MoFEMFunctionBeginHot;

  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;

  // Tag th_ray_dir;
  // double def_val[] = {0,0,0};
  // rval = moab.tag_get_handle(
  //   "RAY_DIR_2",3,MB_TYPE_DOUBLE,th_ray_dir,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  // ); CHKERRG(rval);

  // Range out_side_vertices;
  for (map<EntityHandle, TreeData>::iterator mit = verticesOnCutEdges.begin();
       mit != verticesOnCutEdges.end(); mit++) {
    double dist = mit->second.dIst;
    // cout << s << " " << mit->second.dIst << " " << mit->second.lEngth <<
    // endl;
    VectorDouble3 new_coors =
        mit->second.rayPoint + dist * mit->second.unitRayDir;
    if (th) {
      rval = moab.tag_set_data(th, &mit->first, 1, &new_coors[0]);
      CHKERRG(rval);
    } else {
      rval = moab.set_coords(&mit->first, 1, &new_coors[0]);
      CHKERRG(rval);
    }
    // rval =
    // moab.tag_set_data(th_ray_dir,&mit->first,1,&mit->second.unitRayDir[0]);
    // CHKERRG(rval);
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CutMeshInterface::moveMidNodesOnTrimmedEdges(Tag th) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  // Range out_side_vertices;
  for (map<EntityHandle, TreeData>::iterator mit = verticesOnTrimEdges.begin();
       mit != verticesOnTrimEdges.end(); mit++) {
    double dist = mit->second.dIst;
    // cout << s << " " << mit->second.dIst << " " << mit->second.lEngth <<
    // endl;
    VectorDouble3 new_coors =
        mit->second.rayPoint + dist * mit->second.unitRayDir;
    if (th) {
      rval = moab.tag_set_data(th, &mit->first, 1, &new_coors[0]);
      CHKERRG(rval);
    } else {
      rval = moab.set_coords(&mit->first, 1, &new_coors[0]);
      CHKERRG(rval);
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CutMeshInterface::findEdgesToTrim(Tag th, const double tol,
                                                 int verb) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  // takes edges on body skin
  Skinner skin(&moab);
  Range tets_skin;
  CHKERR skin.find_skin(0, cutNewVolumes, false, tets_skin);
  // vertives on the skin
  Range tets_skin_verts;
  CHKERR moab.get_connectivity(tets_skin,tets_skin_verts,true);
  // edges on the skin
  Range adj_edges_tets_skin;
  CHKERR moab.get_adjacencies(tets_skin, 1, false, adj_edges_tets_skin,
                              moab::Interface::UNION);
  // get edges on new surface
  Range edges;
  CHKERR moab.get_adjacencies(cutNewSurfaces, 1, false, edges,
                              moab::Interface::UNION);

  // clear data ranges
  trimEdges.clear();
  edgesToTrim.clear();
  verticesOnTrimEdges.clear();
  trimNewVertices.clear();

  // Solve nonlinera problem of finding point on surface front
  struct ClosestPointProjection {
    boost::shared_ptr<OrientedBoxTreeTool> treeSurfPtr;
    EntityHandle rootSetSurf;
    VectorDouble3 S0;
    VectorDouble3 rAY;
    double pointOut[3];
    VectorAdaptor vecPointOut;
    ClosestPointProjection(boost::shared_ptr<OrientedBoxTreeTool> &tree,
                            EntityHandle root_set, VectorDouble3 &s0,
                            VectorDouble3 &ray)
        : treeSurfPtr(tree), rootSetSurf(root_set), S0(s0), rAY(ray),
          vecPointOut(3,
                      ublas::shallow_array_adaptor<double>(3, &pointOut[0])) {}
    VectorDouble3 &operator()(const int max_it, const double tol) {
      // cerr << "\n\n" << endl;
      VectorDouble3 w;
      double length = norm_2(rAY);
      rAY /= length;
      for (int ii = 0; ii != max_it; ii++) {
        EntityHandle facets_out;
        treeSurfPtr->closest_to_location(&S0[0], rootSetSurf, pointOut,
                                         facets_out);
        w = vecPointOut - S0;
        double s = inner_prod(rAY, w);
        S0 += s * rAY;
        // cerr << "s " << ii << " " << s << " " << norm_2(w) << endl;
        if (s / length < tol)
          break;
      }
      return S0;
    }
  };

  // iterate over entities on new cut surface
  for (Range::iterator eit = edges.begin(); eit != edges.end(); eit++) {
    // Get edge connectivity and coords
    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(*eit, conn, num_nodes, true);
    double coords[3 * num_nodes];
    if (th) {
      CHKERR moab.tag_get_data(th, conn, num_nodes, coords);
    } else {
      CHKERR moab.get_coords(conn, num_nodes, coords);
    }
    // Put edges coords into boost vectors
    VectorAdaptor s0(3, ublas::shallow_array_adaptor<double>(3, &coords[0]));
    VectorAdaptor s1(3, ublas::shallow_array_adaptor<double>(3, &coords[3]));
    // get edge length
    double length = norm_2(s1 - s0);
    // Find point on surface closet to surface
    double point_out0[3];
    EntityHandle facets_out0;
    // find closet point on the surface from first node
    CHKERR treeSurfPtr->closest_to_location(&coords[0], rootSetSurf, point_out0,
                                            facets_out0);
    // find closest point on the surface from second node
    double point_out1[3];
    EntityHandle facets_out1;
    CHKERR treeSurfPtr->closest_to_location(&coords[3], rootSetSurf, point_out1,
                                            facets_out1);
    // Put closest point in boost vectors
    VectorAdaptor p0(3, ublas::shallow_array_adaptor<double>(3, point_out0));
    VectorAdaptor p1(3, ublas::shallow_array_adaptor<double>(3, point_out1));
    // Calculate deltas, i.e. vectors from edges to closet point on surface
    VectorDouble3 delta0, delta1;
    delta0 = p0 - s0;
    delta1 = p1 - s1;
    // moab.tag_set_data(th,&conn[0],1,&delta0[0]);
    // moab.tag_set_data(th,&conn[1],1,&delta1[0]);
    // Calculate distances
    double dist0 = norm_2(delta0);
    double dist1 = norm_2(delta1);
    double min_dist = fmin(dist0, dist1);
    double max_dist = fmax(dist0, dist1);

    // add edge to trim
    double dist;
    VectorDouble3 ray;
    VectorDouble3 trimmed_end;
    VectorDouble3 itersection_point;
    EntityHandle vert;

    if (min_dist / length < tol && max_dist / length > 0) {
      if (max_dist == dist0) {
        // move mid node in reference to node 0
        vert        = conn[1];
        trimmed_end = s0;
        ray         = s1 - trimmed_end;
      } else {
        // move node in reference to node 1
        vert        = conn[0];
        trimmed_end = s1;
        ray         = s0 - trimmed_end;
      }

      itersection_point =
          ClosestPointProjection(treeSurfPtr, rootSetSurf, trimmed_end,
                                 ray)(nbMaxTrimSearchIterations, tol * 1e-4);
      ray  = itersection_point - trimmed_end;
      dist = norm_2(ray);

      // If one of nodes is on the surface and other is not, that edge is to
      // trim
      if (min_dist / length < tol && max_dist / length > tol) {
        // check if edges should be trimmed, i.e. if edge is trimmed at very end
        // simply move closed node rather than trim
        if ((1 - dist / length) > tol) {
          edgesToTrim[*eit].dIst       = dist;
          edgesToTrim[*eit].lEngth     = dist;
          edgesToTrim[*eit].unitRayDir = ray / dist;
          edgesToTrim[*eit].rayPoint   = trimmed_end;
          trimEdges.insert(*eit);
        } else if ((1 - dist / length) > 0) {
          trimNewVertices.insert(vert);
          verticesOnTrimEdges[vert].dIst       = dist;
          verticesOnTrimEdges[vert].lEngth     = dist;
          verticesOnTrimEdges[vert].unitRayDir = ray / dist;
          verticesOnTrimEdges[vert].rayPoint   = trimmed_end;
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::trimEdgesInTheMiddle(const BitRefLevel bit,
                                                      Tag th, const double tol,
                                                      const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab  = m_field.get_moab();
  MeshRefinement *refiner;
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBegin;

  CHKERR m_field.getInterface(refiner);
  CHKERR m_field.get_ref_ents(&ref_ents_ptr);
  CHKERR refiner->add_verices_in_the_middel_of_edges(trimEdges, bit);
  CHKERR refiner->refine_TET(cutNewVolumes, bit, false, QUIET,
                             debug ? &trimEdges : NULL);

  trimNewVolumes.clear();
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, bit, MBTET, trimNewVolumes);
  // Get vertices which are on trim edges
  for (map<EntityHandle, TreeData>::iterator mit = edgesToTrim.begin();
       mit != edgesToTrim.end(); mit++) {
    RefEntity_multiIndex::index<
        Composite_ParentEnt_And_EntType_mi_tag>::type::iterator vit =
        ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>().find(
            boost::make_tuple(MBVERTEX, mit->first));
    if (vit ==
        ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>().end()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "No vertex on trim edges, that make no sense");
    }
    const boost::shared_ptr<RefEntity> &ref_ent = *vit;
    if ((ref_ent->getBitRefLevel() & bit).any()) {
      EntityHandle vert = ref_ent->getRefEnt();
      trimNewVertices.insert(vert);
      verticesOnTrimEdges[vert] = mit->second;
    }
  }

  // Get faces which are trimmed
  trimNewSurfaces.clear();
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, bit, MBTRI, trimNewSurfaces);
  Range trim_new_surfaces_nodes;
  CHKERR moab.get_connectivity(trimNewSurfaces, trim_new_surfaces_nodes, true);
  trim_new_surfaces_nodes = subtract(trim_new_surfaces_nodes, trimNewVertices);
  trim_new_surfaces_nodes = subtract(trim_new_surfaces_nodes, cutNewVertices);
  Range faces_not_on_surface;
  CHKERR moab.get_adjacencies(trim_new_surfaces_nodes, 2, false,
                              faces_not_on_surface, moab::Interface::UNION);
  trimNewSurfaces = subtract(trimNewSurfaces, faces_not_on_surface);

  // Get surfaces which are not trimmed and add them to surface
  Range all_surfaces_on_bit_level;
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, BitRefLevel().set(), MBTRI, all_surfaces_on_bit_level);
  all_surfaces_on_bit_level =
      intersect(all_surfaces_on_bit_level, cutNewSurfaces);
  trimNewSurfaces = unite(trimNewSurfaces, all_surfaces_on_bit_level);

  // check of nodes are outside surface and if it are remove adjacent faces to
  // those nodes.
  Range check_verts;
  CHKERR moab.get_connectivity(trimNewSurfaces, check_verts, true);
  check_verts = subtract(check_verts, trimNewVertices);
  for (Range::iterator vit = check_verts.begin(); vit != check_verts.end();
       vit++) {
    double coords[3];
    if (th) {
      CHKERR moab.tag_get_data(th, &*vit, 1, coords);
    } else {
      CHKERR moab.get_coords(&*vit, 1, coords);
    }
    double point_out[3];
    EntityHandle facets_out;
    CHKERR treeSurfPtr->closest_to_location(coords, rootSetSurf, point_out,
                                            facets_out);
    VectorAdaptor s(3, ublas::shallow_array_adaptor<double>(3, coords));
    VectorAdaptor p(3, ublas::shallow_array_adaptor<double>(3, point_out));
    if (norm_2(s - p) / aveLength > tol) {
      Range adj;
      CHKERR moab.get_adjacencies(&*vit, 1, 2, false, adj);
      trimNewSurfaces = subtract(trimNewSurfaces, adj);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::getRayForEdge(const EntityHandle ent,
                                               VectorAdaptor ray_point,
                                               VectorAdaptor unit_ray_dir,
                                               double &ray_length) const {
  const CoreInterface &m_field = cOre;
  const moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  int num_nodes;
  const EntityHandle *conn;
  rval = moab.get_connectivity(ent, conn, num_nodes, true);
  CHKERRG(rval);
  double coords[6];
  rval = moab.get_coords(conn, num_nodes, coords);
  CHKERRG(rval);
  VectorAdaptor s0(3, ublas::shallow_array_adaptor<double>(3, &coords[0]));
  VectorAdaptor s1(3, ublas::shallow_array_adaptor<double>(3, &coords[3]));
  noalias(ray_point) = s0;
  noalias(unit_ray_dir) = s1 - s0;
  ray_length = norm_2(unit_ray_dir);
  unit_ray_dir /= ray_length;
  MoFEMFunctionReturnHot(0);
}

// int CutMeshInterface::segmentPlane(
//   VectorAdaptor s0,
//   VectorAdaptor s1,
//   VectorAdaptor x0,
//   VectorAdaptor n,
//   double &s
// ) const {
//   VectorDouble3 u = s1 - s0;
//   VectorDouble3 w = s0 - x0;
//   double nu = inner_prod(n,u);
//   double nw = -inner_prod(n,w);
//   const double tol = 1e-4;
//   if (fabs(nu) < tol) {           // segment is parallel to plane
//       if (nw == 0)                      // segment lies in plane
//           return 2;
//       else
//           return 0;                    // no intersection
//   }
//   // they are not parallel
//   // compute intersect param
//   s = nw / nu;
//   if (s < 0 || s > 1)
//       return 0;                        // no intersection
//   return 1;
// }

MoFEMErrorCode
CutMeshInterface::removePathologicalFrontTris(const BitRefLevel split_bit,
                                              Range &ents) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  PrismInterface *interface;
  MoFEMFunctionBegin;
  CHKERR m_field.getInterface(interface);
  // Remove tris on surface front 
  {
    Range front_tris;
    EntityHandle meshset;
    CHKERR moab.create_meshset(MESHSET_SET, meshset);
    CHKERR moab.add_entities(meshset, ents);
    CHKERR interface->findIfTringleHasThreeNodesOnInternalSurfaceSkin(
        meshset, split_bit, true, front_tris);
    CHKERR moab.delete_entities(&meshset, 1);
    ents = subtract(ents, front_tris);
  }
  // Remove entities on skin
  Range tets;
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
    split_bit,BitRefLevel().set(),MBTET,tets
  );
  // Remove entities on skin
  Skinner skin(&moab);
  Range tets_skin;
  rval = skin.find_skin(0, tets, false, tets_skin);
  ents = subtract(ents, tets_skin);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::splitSides(const BitRefLevel split_bit,
                                            const BitRefLevel bit,
                                            const Range &ents, Tag th) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  PrismInterface *interface;
  MoFEMFunctionBegin;
  ierr = m_field.getInterface(interface);
  CHKERRG(ierr);
  EntityHandle meshset_volume;
  rval = moab.create_meshset(MESHSET_SET, meshset_volume);
  CHKERRG(rval);
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      split_bit, BitRefLevel().set(), MBTET, meshset_volume);
  CHKERRG(ierr);
  EntityHandle meshset_surface;
  rval = moab.create_meshset(MESHSET_SET, meshset_surface);
  CHKERRG(rval);
  rval = moab.add_entities(meshset_surface, ents);
  CHKERRG(rval);
  ierr = interface->getSides(meshset_surface, split_bit, true);
  CHKERRG(ierr);
  ierr =
      interface->splitSides(meshset_volume, bit, meshset_surface, true, true);
  CHKERRG(ierr);
  rval = moab.delete_entities(&meshset_volume, 1);
  CHKERRG(rval);
  rval = moab.delete_entities(&meshset_surface, 1);
  CHKERRG(rval);
  if (th) {
    Range prisms;
    ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit, BitRefLevel().set(), MBPRISM, prisms);
    CHKERRG(ierr);
    for (Range::iterator pit = prisms.begin(); pit != prisms.end(); pit++) {
      int num_nodes;
      const EntityHandle *conn;
      rval = moab.get_connectivity(*pit, conn, num_nodes, true);
      CHKERRG(rval);
      MatrixDouble data(3, 3);
      rval = moab.tag_get_data(th, conn, 3, &data(0, 0));
      CHKERRG(rval);
      // cerr << data << endl;
      rval = moab.tag_set_data(th, &conn[3], 3, &data(0, 0));
      CHKERRG(rval);
    }
  }
  MoFEMFunctionReturn(0);
}

struct LengthMapData {
  double lEngth;
  EntityHandle eDge;
  bool skip;
  LengthMapData(const double l, const EntityHandle e)
      : lEngth(l), eDge(e), skip(false) {}
};

typedef multi_index_container<
  LengthMapData, 
  indexed_by<
    ordered_non_unique<
      member<LengthMapData, double, &LengthMapData::lEngth>
    >,
    hashed_unique<
      member<LengthMapData, EntityHandle, &LengthMapData::eDge>
    >
  >
> LengthMapData_multi_index;



MoFEMErrorCode CutMeshInterface::mergeBadEdges(
    const int fraction_level, const Range &tets, const Range &surface,
    const Range &fixed_edges, const Range &corner_nodes, Range &edges_to_merge,
    Range &out_tets, Range &new_surf, Tag th, const bool update_meshsets,
    const BitRefLevel *bit_ptr,const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  /**
   * \brief Merge nodes
   */
  struct MergeNodes {
    CoreInterface &mField;
    const bool onlyIfImproveQuality;
    const int lineSearch;
    Tag tH;
    bool updateMehsets;

    MergeNodes(CoreInterface &m_field,
               const bool only_if_improve_quality, const int line_search,
               Tag th, bool update_mehsets)
        : mField(m_field), onlyIfImproveQuality(only_if_improve_quality),
          lineSearch(line_search), tH(th), updateMehsets(update_mehsets) {
      mField.getInterface(nodeMergerPtr);
    }
    NodeMergerInterface *nodeMergerPtr;
    MoFEMErrorCode operator()(EntityHandle father, EntityHandle mother,
                              Range &proc_tets, Range &new_surf,
                              Range &edges_to_merge, Range &not_merged_edges,
                              bool add_child = true) const {
      moab::Interface &moab = mField.get_moab();
      MoFEMFunctionBeginHot;
      const EntityHandle conn[] = {father, mother};
      Range vert_tets;
      rval = moab.get_adjacencies(conn, 2, 3, false, vert_tets,
                                  moab::Interface::UNION);
      CHKERRG(rval);
      vert_tets = intersect(vert_tets, proc_tets);
      Range out_tets;
      ierr = nodeMergerPtr->mergeNodes(father, mother, out_tets, &vert_tets,
                                       onlyIfImproveQuality, 0, lineSearch, tH);
      CHKERRG(ierr);
      out_tets.merge(subtract(proc_tets, vert_tets));
      proc_tets.swap(out_tets);

      if (add_child && nodeMergerPtr->getSucessMerge()) {

        NodeMergerInterface::ParentChildMap &parent_child_map =
            nodeMergerPtr->getParentChildMap();

        Range child_ents;
        NodeMergerInterface::ParentChildMap::iterator it;
        for (it = parent_child_map.begin(); it != parent_child_map.end();
             it++) {
          child_ents.insert(it->pArent);
        }

        Range new_surf_child_ents = intersect(new_surf, child_ents);
        new_surf = subtract(new_surf, new_surf_child_ents);
        Range child_surf_ents;
        ierr = updateRangeByChilds(parent_child_map, new_surf_child_ents,
                                   child_surf_ents);
        CHKERRG(ierr);
        new_surf.merge(child_surf_ents);

        Range edges_to_merge_child_ents = intersect(edges_to_merge, child_ents);
        edges_to_merge = subtract(edges_to_merge, edges_to_merge_child_ents);
        Range merged_child_edge_ents;
        ierr = updateRangeByChilds(parent_child_map, edges_to_merge_child_ents,
                                   merged_child_edge_ents);
        CHKERRG(ierr);

        Range not_merged_edges_child_ents =
            intersect(not_merged_edges, child_ents);
        not_merged_edges =
            subtract(not_merged_edges, not_merged_edges_child_ents);
        Range not_merged_child_edge_ents;
        ierr =
            updateRangeByChilds(parent_child_map, not_merged_edges_child_ents,
                                not_merged_child_edge_ents);
        CHKERRG(ierr);

        merged_child_edge_ents =
            subtract(merged_child_edge_ents, not_merged_child_edge_ents);
        edges_to_merge.merge(merged_child_edge_ents);
        not_merged_edges.merge(not_merged_child_edge_ents);

        if (updateMehsets) {
          for (_IT_CUBITMESHSETS_FOR_LOOP_(
                   (*mField.getInterface<MeshsetsManager>()), cubit_it)) {
            EntityHandle cubit_meshset = cubit_it->meshset;
            Range parent_ents;
            rval =
                moab.get_entities_by_handle(cubit_meshset, parent_ents, true);
            CHKERRG(rval);
            Range child_ents;
            ierr =
                updateRangeByChilds(parent_child_map, parent_ents, child_ents);
            CHKERRG(ierr);
            rval = moab.add_entities(cubit_meshset, child_ents);
            CHKERRG(rval);
          }
        }
      }
      MoFEMFunctionReturnHot(0);
    }

  private:
    MoFEMErrorCode updateRangeByChilds(
        const NodeMergerInterface::ParentChildMap &parent_child_map,
        const Range &parents, Range &childs) const {
      MoFEMFunctionBeginHot;
      NodeMergerInterface::ParentChildMap::nth_index<0>::type::iterator it;
      for (Range::const_iterator eit = parents.begin(); eit != parents.end();
           eit++) {
        it = parent_child_map.get<0>().find(*eit);
        if (it == parent_child_map.get<0>().end())
          continue;
        childs.insert(it->cHild);
      }
      MoFEMFunctionReturnHot(0);
    }
  };

  /**
   * \brief Calculate edge length
   */
  struct LengthMap {
    Tag tH;
    CoreInterface &mField;
    moab::Interface &moab;
    const double maxLength;
    LengthMap(CoreInterface &m_field, Tag th, double max_length)
        : tH(th), mField(m_field), moab(m_field.get_moab()),
          maxLength(max_length) {}

    MoFEMErrorCode operator()(const Range &tets, const Range &edges,
                              LengthMapData_multi_index &length_map) const {
      int num_nodes;
      const EntityHandle *conn;
      double coords[6];
      MoFEMFunctionBegin;
      VectorAdaptor s0(3, ublas::shallow_array_adaptor<double>(3, &coords[0]));
      VectorAdaptor s1(3, ublas::shallow_array_adaptor<double>(3, &coords[3]));
      VectorDouble3 delta(3);
      for (Range::const_iterator eit = edges.begin(); eit != edges.end();
           eit++) {
        Range adj_tet;
        CHKERR moab.get_adjacencies(&*eit, 1, 3, false, adj_tet);
        adj_tet = intersect(adj_tet, tets);
        CHKERR moab.get_connectivity(*eit, conn, num_nodes, true);
        if (tH) {
          CHKERR moab.tag_get_data(tH, conn, num_nodes, coords);
        } else {
          CHKERR moab.get_coords(conn, num_nodes, coords);
        }
        double q = 1;
        CHKERR mField.getInterface<Tools>()->minTetsQuality(adj_tet, q, tH);
        if (q != q)
          q = -1;
        noalias(delta) = (s0 - s1) / maxLength;
        double dot = inner_prod(delta, delta);
        length_map.insert(LengthMapData(pow(q, 2) * sqrt(dot), *eit));
      }
      MoFEMFunctionReturn(0);
    }
  };

  /**
   * \brief Topological relations
   */
  struct Toplogy {

    CoreInterface &mField;
    Tag tH;
    const double tOL;
    Toplogy(CoreInterface &m_field, Tag th, const double tol)
        : mField(m_field), tH(th), tOL(tol) {}

    enum TYPE {
      FREE         = 0,
      SKIN         = 1 << 0,
      SURFACE      = 1 << 1,
      SURFACE_SKIN = 1 << 2,
      FRONT_ENDS   = 1 << 3,
      FIX_EDGES    = 1 << 4,
      FIX_CORNERS  = 1 << 5
    };

    typedef map<int, Range> SetsMap;

    MoFEMErrorCode classifyVerts(const Range &surface, const Range &tets,
                                 const Range &fixed_edges,
                                 const Range &corner_nodes,
                                 SetsMap &sets_map) const {
      moab::Interface &moab(mField.get_moab());
      Skinner skin(&moab);
      MoFEMFunctionBeginHot;

      sets_map[FIX_CORNERS].merge(corner_nodes);
      Range fixed_verts;
      rval = moab.get_connectivity(fixed_edges, fixed_verts, true);
      CHKERRG(rval);
      sets_map[FIX_EDGES].swap(fixed_verts);

      Range tets_skin;
      rval = skin.find_skin(0, tets, false, tets_skin);
      CHKERRG(rval);
      Range tets_skin_edges;
      rval = moab.get_adjacencies(tets_skin, 1, false, tets_skin_edges,
                                  moab::Interface::UNION);
      CHKERRG(rval);

      // surface skin
      Range surface_skin;
      rval = skin.find_skin(0, surface, false, surface_skin);
      CHKERRG(rval);
      Range front_in_the_body;
      front_in_the_body = subtract(surface_skin, tets_skin_edges);
      Range front_ends;
      rval = skin.find_skin(0, front_in_the_body, false, front_ends);
      CHKERRG(rval);
      sets_map[FRONT_ENDS].swap(front_ends);

      Range surface_skin_verts;
      rval = moab.get_connectivity(surface_skin, surface_skin_verts, true);
      CHKERRG(rval);
      sets_map[SURFACE_SKIN].swap(surface_skin_verts);

      // surface
      Range surface_verts;
      rval = moab.get_connectivity(surface, surface_verts, true);
      CHKERRG(rval);
      sets_map[SURFACE].swap(surface_verts);

      // skin
      Range tets_skin_verts;
      rval = moab.get_connectivity(tets_skin, tets_skin_verts, true);
      CHKERRG(rval);
      sets_map[SKIN].swap(tets_skin_verts);

      Range tets_verts;
      rval = moab.get_connectivity(tets, tets_verts, true);
      CHKERRG(rval);
      sets_map[FREE].swap(tets_verts);

      MoFEMFunctionReturnHot(0);
    }

    MoFEMErrorCode getProcTets(const Range &tets, const Range &edges_to_merge,
                               Range &proc_tets) const {
      moab::Interface &moab(mField.get_moab());
      MoFEMFunctionBeginHot;
      Range edges_to_merge_verts;
      rval = moab.get_connectivity(edges_to_merge, edges_to_merge_verts, true);
      CHKERRG(rval);
      Range edges_to_merge_verts_tets;
      rval = moab.get_adjacencies(edges_to_merge_verts, 3, false,
                                  edges_to_merge_verts_tets,
                                  moab::Interface::UNION);
      CHKERRG(rval);
      edges_to_merge_verts_tets = intersect(edges_to_merge_verts_tets, tets);
      proc_tets.swap(edges_to_merge_verts_tets);
      MoFEMFunctionReturnHot(0);
    }

    MoFEMErrorCode edgesToMerge(const Range &surface, const Range &tets,
                                Range &edges_to_merge) const {
      moab::Interface &moab(mField.get_moab());
      MoFEMFunctionBeginHot;

      Range surface_verts;
      rval = moab.get_connectivity(surface, surface_verts, true);
      CHKERRG(rval);
      Range surface_verts_edges;
      rval = moab.get_adjacencies(surface_verts, 1, false, surface_verts_edges,
                                  moab::Interface::UNION);
      CHKERRG(rval);
      edges_to_merge.merge(surface_verts_edges);
      Range tets_edges;
      rval = moab.get_adjacencies(tets, 1, false, tets_edges,
                                  moab::Interface::UNION);
      CHKERRG(rval);
      edges_to_merge = intersect(edges_to_merge, tets_edges);
      MoFEMFunctionReturnHot(0);
    }

    MoFEMErrorCode removeBadEdges(const Range &surface, const Range &tets,
                                  const Range &fixed_edges,
                                  const Range &corner_nodes,
                                  Range &edges_to_merge,
                                  Range &not_merged_edges) {
      moab::Interface &moab(mField.get_moab());
      MoFEMFunctionBeginHot;

      // find skin
      Skinner skin(&moab);
      Range tets_skin;
      rval = skin.find_skin(0, tets, false, tets_skin);
      CHKERRG(rval);
      Range surface_skin;
      rval = skin.find_skin(0, surface, false, surface_skin);
      CHKERRG(rval);

      // end nodes
      Range tets_skin_edges;
      rval = moab.get_adjacencies(tets_skin, 1, false, tets_skin_edges,
                                  moab::Interface::UNION);
      CHKERRG(rval);
      Range surface_front;
      surface_front = subtract(surface_skin, tets_skin_edges);
      Range ends_nodes;
      rval = skin.find_skin(0, surface_front, false, ends_nodes);
      CHKERRG(rval);

      // remove bad merges

      // get surface and body skin verts
      Range surface_edges;
      rval = moab.get_adjacencies(surface, 1, false, surface_edges,
                                  moab::Interface::UNION);
      CHKERRG(rval);
      // get nodes on the surface
      Range surface_edges_verts;
      rval = moab.get_connectivity(surface_edges, surface_edges_verts, true);
      CHKERRG(rval);
      // get vertices on the body skin
      Range tets_skin_edges_verts;
      rval =
          moab.get_connectivity(tets_skin_edges, tets_skin_edges_verts, true);
      CHKERRG(rval);

      Range edges_to_remove;

      // remove edges self connected to body skin
      {
        Range ents_nodes_and_edges;
        ents_nodes_and_edges.merge(tets_skin_edges_verts);
        ents_nodes_and_edges.merge(tets_skin_edges);
        ierr = removeSelfConectingEdges(ents_nodes_and_edges, edges_to_remove,
                                        0, false);
        CHKERRG(ierr);
      }
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      // remove edges self connected to surface
      {
        Range ents_nodes_and_edges;
        ents_nodes_and_edges.merge(surface_edges_verts);
        ents_nodes_and_edges.merge(surface_edges);
        ents_nodes_and_edges.merge(tets_skin_edges_verts);
        ents_nodes_and_edges.merge(tets_skin_edges);
        ierr = removeSelfConectingEdges(ents_nodes_and_edges, edges_to_remove,
                                        0, false);
        CHKERRG(ierr);
      }
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      // remove edges adjacent corner_nodes execpt those on fixed edges
      Range fixed_edges_nodes;
      rval = moab.get_connectivity(fixed_edges, fixed_edges_nodes, true);
      {
        Range ents_nodes_and_edges;
        ents_nodes_and_edges.merge(fixed_edges_nodes);
        ents_nodes_and_edges.merge(ends_nodes);
        ents_nodes_and_edges.merge(corner_nodes);
        ents_nodes_and_edges.merge(fixed_edges);
        ierr = removeSelfConectingEdges(ents_nodes_and_edges, edges_to_remove,
                                        0, false);
        CHKERRG(ierr);
      }
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      // remove edges self connected to surface
      ierr = removeSelfConectingEdges(surface_edges, edges_to_remove, 0, false);
      CHKERRG(ierr);
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      // remove edges self contented on surface skin
      {
        Range ents_nodes_and_edges;
        ents_nodes_and_edges.merge(surface_skin);
        ents_nodes_and_edges.merge(fixed_edges_nodes);
        ierr = removeSelfConectingEdges(ents_nodes_and_edges, edges_to_remove,
                                        0, false);
        CHKERRG(ierr);
      }
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      // remove edges connecting crack front and fixed edges, except those short
      {
        Range ents_nodes_and_edges;
        ents_nodes_and_edges.merge(surface_skin.subset_by_type(MBEDGE));
        ents_nodes_and_edges.merge(fixed_edges.subset_by_type(MBEDGE));
        ierr = removeSelfConectingEdges(ents_nodes_and_edges, edges_to_remove,
                                        tOL, false);
        CHKERRG(ierr);
      }
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      MoFEMFunctionReturnHot(0);
    }

  private:
    MoFEMErrorCode removeSelfConectingEdges(const Range &ents,
                                            Range &edges_to_remove,
                                            const bool length,
                                            bool debug) const {
      moab::Interface &moab(mField.get_moab());
      MoFEMFunctionBegin;
      // get nodes
      Range ents_nodes = ents.subset_by_type(MBVERTEX);
      if (ents_nodes.empty()) {
        CHKERR moab.get_connectivity(ents, ents_nodes, true);
      }
      // edges adj. to nodes
      Range ents_nodes_edges;
      CHKERR moab.get_adjacencies(ents_nodes, 1, false, ents_nodes_edges,
                                  moab::Interface::UNION);
      // nodes of adj. edges
      Range ents_nodes_edges_nodes;
      CHKERR moab.get_connectivity(ents_nodes_edges, ents_nodes_edges_nodes,
                                   true);
      // hanging nodes
      ents_nodes_edges_nodes = subtract(ents_nodes_edges_nodes, ents_nodes);
      Range ents_nodes_edges_nodes_edges;
      CHKERR moab.get_adjacencies(ents_nodes_edges_nodes, 1, false,
                                  ents_nodes_edges_nodes_edges,
                                  moab::Interface::UNION);
      // remove edges adj. to hanging edges
      ents_nodes_edges =
          subtract(ents_nodes_edges, ents_nodes_edges_nodes_edges);
      ents_nodes_edges =
          subtract(ents_nodes_edges, ents.subset_by_type(MBEDGE));
      if(length>0) {
        Range::iterator eit = ents_nodes_edges.begin();
        for (; eit != ents_nodes_edges.end();) {
          int num_nodes;
          const EntityHandle *conn;
          rval = moab.get_connectivity(*eit, conn, num_nodes, true);
          double coords[6];
          if(tH) {
            CHKERR moab.tag_get_data(tH, conn, num_nodes, coords);
          } else {
            CHKERR moab.get_coords(conn, num_nodes, coords);
          }
          VectorAdaptor s0(3,
                           ublas::shallow_array_adaptor<double>(3, &coords[0]));
          VectorAdaptor s1(3,
                           ublas::shallow_array_adaptor<double>(3, &coords[3]));
          const double edge_length = norm_2(s0-s1);
          if(edge_length<tOL) {
            eit = ents_nodes_edges.erase(eit);
          } else {
            eit++;
          }
        }
      }
      edges_to_remove.swap(ents_nodes_edges);
      if (debug) {
        EntityHandle meshset;
        CHKERR moab.create_meshset(MESHSET_SET, meshset);
        CHKERR moab.add_entities(meshset, edges_to_remove);
        CHKERR moab.write_file("edges_to_remove.vtk", "VTK", "", &meshset, 1);
        CHKERR moab.delete_entities(&meshset, 1);
      }
      MoFEMFunctionReturn(0);
    }
  };

  Range not_merged_edges;
  CHKERR Toplogy(m_field, th, 0.1 * aveLength)
      .edgesToMerge(surface, tets, edges_to_merge);
  CHKERR Toplogy(m_field, th, 0.1 * aveLength)
      .removeBadEdges(surface, tets, fixed_edges, corner_nodes, edges_to_merge,
                      not_merged_edges);
  Toplogy::SetsMap sets_map;
  CHKERR Toplogy(m_field, th, 0.1 * aveLength)
      .classifyVerts(surface, tets, fixed_edges, corner_nodes, sets_map);
  Range proc_tets;
  CHKERR Toplogy(m_field, th, 0.1 * aveLength)
      .getProcTets(tets, edges_to_merge, proc_tets);
  out_tets = subtract(tets, proc_tets);
  if (bit_ptr) {
    for (int dd = 2; dd >= 0; dd--) {
      CHKERR moab.get_adjacencies(out_tets.subset_by_dimension(3), dd, false,
                                  out_tets, moab::Interface::UNION);
    }
    CHKERR m_field.getInterface<BitRefManager>()->addBitRefLevel(out_tets,
                                                                 *bit_ptr);
  }

  int nb_nodes_merged = 0;
  LengthMapData_multi_index length_map;
  new_surf = surface;

  for (int pp = 0; pp != nbMaxMergingCycles; pp++) {

    int nb_nodes_merged_0 = nb_nodes_merged;
    length_map.clear();
    CHKERR LengthMap(m_field, th, maxLength)(proc_tets, edges_to_merge,
                                             length_map);

    int nn = 0;
    Range collapsed_edges;
    for (LengthMapData_multi_index::nth_index<0>::type::iterator
             mit = length_map.get<0>().begin();
         mit != length_map.get<0>().end(); mit++, nn++) {
      // cerr << mit->lEngth << endl; //" " << mit->eDge << endl;
      if(mit->skip) continue;
      int num_nodes;
      const EntityHandle *conn;
      CHKERR moab.get_connectivity(mit->eDge, conn, num_nodes, true);
      int conn_type[2] = {0, 0};
      for (int nn = 0; nn != 2; nn++) {
        conn_type[nn] = 0;
        for (Toplogy::SetsMap::reverse_iterator sit = sets_map.rbegin();
             sit != sets_map.rend(); sit++) {
          if (sit->second.find(conn[nn]) != sit->second.end()) {
            conn_type[nn] |= sit->first;
          }
        }
      }
      int type_father, type_mother;
      EntityHandle father, mother;
      if (conn_type[0] > conn_type[1]) {
        father = conn[0];
        mother = conn[1];
        type_father = conn_type[0];
        type_mother = conn_type[1];
      } else {
        father = conn[1];
        mother = conn[0];
        type_father = conn_type[1];
        type_mother = conn_type[0];
      }
      int line_search = 0;
      if (type_father == type_mother) {
        line_search = lineSearchSteps;
      }
      
      CHKERR MergeNodes(m_field, true, line_search, th,
                        update_meshsets)(father, mother, proc_tets, new_surf,
                                         edges_to_merge, not_merged_edges);

      if (m_field.getInterface<NodeMergerInterface>()->getSucessMerge()) {
        Range adj_edges;
        CHKERR moab.get_adjacencies(conn, 2, 1, false, adj_edges,
                                    moab::Interface::UNION);
        for (Range::iterator ait = adj_edges.begin(); ait != adj_edges.end();
             ait++) {
          LengthMapData_multi_index::nth_index<1>::type::iterator miit =
              length_map.get<1>().find(*ait);
          if (miit != length_map.get<1>().end()) {
            (const_cast<LengthMapData &>(*miit)).skip = true;
          }
        }
        nb_nodes_merged++;
        collapsed_edges.insert(mit->eDge);
      }

      if (nn > length_map.size() / fraction_level)
        break;
    }

    Range adj_faces, adj_edges;
    CHKERR moab.get_adjacencies(proc_tets, 2, false, adj_faces,
                                moab::Interface::UNION);
    new_surf = intersect(new_surf, adj_faces);

    PetscPrintf(m_field.get_comm(), "(%d) Number of nodes merged %d\n", pp,
                nb_nodes_merged);

    if (debug) {
      // current surface and merged edges in step
      EntityHandle meshset;
      std::string name;
      CHKERR moab.create_meshset(MESHSET_SET, meshset);
      CHKERR moab.add_entities(meshset, new_surf);
      CHKERR moab.add_entities(meshset, collapsed_edges);
      name = "node_merge_step_" + boost::lexical_cast<std::string>(pp) + ".vtk";
      CHKERR moab.write_file(name.c_str(), "VTK", "", &meshset, 1);
      CHKERR moab.delete_entities(&meshset, 1);
      // edges to merge
      CHKERR moab.create_meshset(MESHSET_SET, meshset);
      CHKERR moab.add_entities(meshset, edges_to_merge);
      name = "edges_to_merge_step_" + boost::lexical_cast<std::string>(pp) +
             ".vtk";
      CHKERR moab.write_file(name.c_str(), "VTK", "", &meshset, 1);
      CHKERR moab.delete_entities(&meshset, 1);
    }

    if (nb_nodes_merged == nb_nodes_merged_0)
      break;

    CHKERR moab.get_adjacencies(proc_tets, 1, false, adj_edges,
                                moab::Interface::UNION);
    edges_to_merge = intersect(edges_to_merge, adj_edges);
    CHKERR Toplogy(m_field, th, 0.1 * aveLength)
        .removeBadEdges(new_surf, proc_tets, fixed_edges, corner_nodes,
                        edges_to_merge, not_merged_edges);
  }

  if (bit_ptr) {
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevel(proc_tets,
                                                                 *bit_ptr);
  }
  out_tets.merge(proc_tets);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::mergeBadEdges(
    const int fraction_level, const BitRefLevel trim_bit,
    const BitRefLevel cut_bit, const BitRefLevel bit, const Range &surface,
    const Range &fixed_edges, const Range &corner_nodes, Tag th,
    const bool update_meshsets, const bool debug) {
  CoreInterface &m_field = cOre;
  MoFEMFunctionBeginHot;
  Range tets_level;
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      trim_bit, BitRefLevel().set(), MBTET, tets_level);
  CHKERRG(ierr);

  Range edges_to_merge;
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByParentType(
      trim_bit, trim_bit | cut_bit, MBEDGE, edges_to_merge);
  CHKERRG(ierr);
  edges_to_merge = edges_to_merge.subset_by_type(MBEDGE);

  // get all entities not in database
  Range all_ents_not_in_database_before;
  ierr = cOre.getInterface<BitRefManager>()->getAllEntitiesNotInDatabase(
      all_ents_not_in_database_before);
  CHKERRG(ierr);

  Range out_new_tets, out_new_surf;
  ierr = mergeBadEdges(fraction_level, tets_level, surface, fixed_edges,
                       corner_nodes, edges_to_merge, out_new_tets, out_new_surf,
                       th, update_meshsets, &bit, debug);
  CHKERRG(ierr);

  // get all entities not in database after merge
  Range all_ents_not_in_database_after;
  ierr = cOre.getInterface<BitRefManager>()->getAllEntitiesNotInDatabase(
      all_ents_not_in_database_after);
  CHKERRG(ierr);
  // delete hanging entities
  all_ents_not_in_database_after =
      subtract(all_ents_not_in_database_after, all_ents_not_in_database_before);
  Range meshsets;
  CHKERR m_field.get_moab().get_entities_by_type(0, MBENTITYSET, meshsets,
                                                 true);
  for (Range::iterator mit = meshsets.begin(); mit != meshsets.end(); mit++) {
    CHKERR m_field.get_moab().remove_entities(*mit,
                                              all_ents_not_in_database_after);
  }
  m_field.get_moab().delete_entities(all_ents_not_in_database_after);

  mergedVolumes.swap(out_new_tets);
  mergedSurfaces.swap(out_new_surf);
  MoFEMFunctionReturnHot(0);
}

#ifdef WITH_TETGEN

MoFEMErrorCode CutMeshInterface::rebuildMeshWithTetGen(
    vector<string> &switches, const BitRefLevel &mesh_bit,
    const BitRefLevel &bit, const Range &surface, const Range &fixed_edges,
    const Range &corner_nodes, Tag th, const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  TetGenInterface *tetgen_iface;
  MoFEMFunctionBegin;
  CHKERR m_field.getInterface(tetgen_iface);

  tetGenData.clear();
  moabTetGenMap.clear();
  tetGenMoabMap.clear();

  if (tetGenData.size() < 1) {
    tetGenData.push_back(new tetgenio);
  }
  tetgenio &in = tetGenData.back();

  struct BitEnts {

    CoreInterface &mField;
    const BitRefLevel &bIt;
    BitEnts(CoreInterface &m_field, const BitRefLevel &bit)
        : mField(m_field), bIt(bit) {}

    Range mTets;
    Range mTris;
    Range mEdges;
    Range mNodes;

    MoFEMErrorCode getLevelEnts() {
      MoFEMFunctionBeginHot;
      CHKERR mField.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
          bIt, BitRefLevel().set(), MBTET, mTets);
      CHKERR mField.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
          bIt, BitRefLevel().set(), MBTRI, mTris);
      CHKERR mField.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
          bIt, BitRefLevel().set(), MBEDGE, mEdges);
      CHKERR mField.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
          bIt, BitRefLevel().set(), MBVERTEX, mNodes);
      MoFEMFunctionReturnHot(0);
    }

    Range mSkin;
    Range mSkinNodes;
    Range mSkinEdges;

    MoFEMErrorCode getSkin() {
      moab::Interface &moab = mField.get_moab();
      MoFEMFunctionBeginHot;
      Skinner skin(&moab);
      CHKERR skin.find_skin(0, mTets, false, mSkin);
      CHKERR mField.get_moab().get_connectivity(mSkin, mSkinNodes, true);
      CHKERR mField.get_moab().get_adjacencies(mSkin, 1, false, mSkinEdges,
                                               moab::Interface::UNION);
      MoFEMFunctionReturnHot(0);
    }
  };

  struct SurfaceEnts {

    CoreInterface &mField;
    SurfaceEnts(CoreInterface &m_field) : mField(m_field) {}

    Range sNodes;
    Range sEdges;
    Range sVols;
    Range vNodes;

    MoFEMErrorCode getVolume(const BitEnts &bit_ents, const Range &tris) {
      moab::Interface &moab = mField.get_moab();
      MoFEMFunctionBeginHot;
      CHKERR moab.get_connectivity(tris, sNodes, true);
      CHKERR moab.get_adjacencies(tris, 1, false, sEdges,
                                  moab::Interface::UNION);
      CHKERR moab.get_adjacencies(sNodes, 3, false, sVols,
                                  moab::Interface::UNION);
      sVols = intersect(sVols, bit_ents.mTets);
      CHKERR moab.get_connectivity(sVols, vNodes, true);
      MoFEMFunctionReturnHot(0);
    }

    Range sSkin;
    Range sSkinNodes;
    Range vSkin;
    Range vSkinNodes;
    Range vSkinWithoutBodySkin;
    Range vSkinNodesWithoutBodySkin;
    Range vSkinOnBodySkin;
    Range vSkinOnBodySkinNodes;

    MoFEMErrorCode getSkin(const BitEnts &bit_ents, const Range &tris,
                           const int levels = 3) {
      moab::Interface &moab = mField.get_moab();
      MoFEMFunctionBeginHot;
      Skinner skin(&moab);
      rval = skin.find_skin(0, sVols, false, vSkin);
      for (int ll = 0; ll != levels; ll++) {
        CHKERR moab.get_adjacencies(vSkin, 3, false, sVols,
                                    moab::Interface::UNION);
        sVols = intersect(sVols, bit_ents.mTets);
        vSkin.clear();
        CHKERR skin.find_skin(0, sVols, false, vSkin);
      }
      vSkinWithoutBodySkin = subtract(vSkin, bit_ents.mSkin);
      vSkinOnBodySkin = intersect(vSkin, bit_ents.mSkin);
      CHKERR moab.get_connectivity(vSkinOnBodySkin, vSkinOnBodySkinNodes, true);
      CHKERR moab.get_connectivity(sVols, vNodes, true);
      CHKERR moab.get_connectivity(vSkin, vSkinNodes, true);
      vSkinNodesWithoutBodySkin = subtract(vSkinNodes, bit_ents.mSkinNodes);
      CHKERR skin.find_skin(0, tris, false, sSkin);
      CHKERR moab.get_connectivity(sSkin, sSkinNodes, true);
      tVols = sVols;
      MoFEMFunctionReturnHot(0);
    }

    Range tVols;

    MoFEMErrorCode getTetsForRemesh(const BitEnts &bit_ents, Tag th = NULL) {
      moab::Interface &moab = mField.get_moab();
      MoFEMFunctionBeginHot;

      Range tets_with_four_nodes_on_skin;
      rval = moab.get_adjacencies(vSkinOnBodySkinNodes, 3, false,
                                  tets_with_four_nodes_on_skin,
                                  moab::Interface::UNION);
      Range tets_nodes;
      CHKERR moab.get_connectivity(tets_with_four_nodes_on_skin, tets_nodes,
                                   true);
      tets_nodes = subtract(tets_nodes, vSkinOnBodySkinNodes);
      Range other_tets;
      CHKERR moab.get_adjacencies(tets_nodes, 3, false, other_tets,
                                  moab::Interface::UNION);
      tets_with_four_nodes_on_skin =
          subtract(tets_with_four_nodes_on_skin, other_tets);
      Range to_remove;
      for (Range::iterator tit = tets_with_four_nodes_on_skin.begin();
           tit != tets_with_four_nodes_on_skin.end(); tit++) {
        int num_nodes;
        const EntityHandle *conn;
        CHKERR moab.get_connectivity(*tit, conn, num_nodes, true);
        double coords[12];
        if (th) {
          CHKERR moab.tag_get_data(th, conn, num_nodes, coords);
        } else {
          CHKERR moab.get_coords(conn, num_nodes, coords);
        }
        double quality = Tools::volumeLengthQuality(coords);
        if (quality < 1e-2) {
          to_remove.insert(*tit);
        }
      }

      sVols = subtract(sVols, to_remove);

      Skinner skin(&moab);
      vSkin.clear();
      CHKERR skin.find_skin(0, sVols, false, vSkin);
      Range m_skin;
      CHKERR
          skin.find_skin(0, subtract(bit_ents.mSkin, to_remove), false, m_skin);
      vSkinWithoutBodySkin = subtract(vSkin, m_skin);
      vSkinOnBodySkin = intersect(vSkin, m_skin);
      vNodes.clear();
      vSkinNodes.clear();
      vSkinOnBodySkinNodes.clear();
      CHKERR moab.get_connectivity(sVols, vNodes, true);
      CHKERR moab.get_connectivity(vSkinOnBodySkin, vSkinOnBodySkinNodes, true);
      CHKERR moab.get_connectivity(vSkin, vSkinNodes, true);
      MoFEMFunctionReturnHot(0);
    }
  };

  BitEnts bit_ents(m_field, mesh_bit);
  CHKERR bit_ents.getLevelEnts();
  CHKERR bit_ents.getSkin();
  SurfaceEnts surf_ents(m_field);
  CHKERR surf_ents.getVolume(bit_ents, surface);
  CHKERR surf_ents.getSkin(bit_ents, surface);
  CHKERR surf_ents.getTetsForRemesh(bit_ents);

  map<int, Range> types_ents;
  types_ents[TetGenInterface::RIDGEVERTEX].merge(
      surf_ents.vSkinNodesWithoutBodySkin);
  // FREESEGVERTEX
  types_ents[TetGenInterface::FREESEGVERTEX].merge(surf_ents.sSkinNodes);
  types_ents[TetGenInterface::FREESEGVERTEX] =
      subtract(types_ents[TetGenInterface::FREESEGVERTEX],
               types_ents[TetGenInterface::RIDGEVERTEX]);
  // FREEFACETVERTEX
  types_ents[TetGenInterface::FREEFACETVERTEX].merge(surf_ents.sNodes);
  types_ents[TetGenInterface::FREEFACETVERTEX] =
      subtract(types_ents[TetGenInterface::FREEFACETVERTEX],
               types_ents[TetGenInterface::RIDGEVERTEX]);
  types_ents[TetGenInterface::FREEFACETVERTEX] =
      subtract(types_ents[TetGenInterface::FREEFACETVERTEX],
               types_ents[TetGenInterface::FREESEGVERTEX]);
  // FREEVOLVERTEX
  types_ents[TetGenInterface::FREEVOLVERTEX].merge(surf_ents.vNodes);
  types_ents[TetGenInterface::FREEVOLVERTEX] =
      subtract(types_ents[TetGenInterface::FREEVOLVERTEX],
               types_ents[TetGenInterface::RIDGEVERTEX]);
  types_ents[TetGenInterface::FREEVOLVERTEX] =
      subtract(types_ents[TetGenInterface::FREEVOLVERTEX],
               types_ents[TetGenInterface::FREESEGVERTEX]);
  types_ents[TetGenInterface::FREEVOLVERTEX] =
      subtract(types_ents[TetGenInterface::FREEVOLVERTEX],
               types_ents[TetGenInterface::FREEFACETVERTEX]);

  Tag th_marker;
  // Clean markers
  rval = m_field.get_moab().tag_get_handle("TETGEN_MARKER", th_marker);
  if(rval == MB_SUCCESS) {
    CHKERR m_field.get_moab().tag_delete(th_marker);
    rval = MB_SUCCESS;
  }

  int def_marker = 0;
  CHKERR m_field.get_moab().tag_get_handle(
      "TETGEN_MARKER", 1, MB_TYPE_INTEGER, th_marker,
      MB_TAG_CREAT | MB_TAG_SPARSE, &def_marker);

  // Mark surface with id = 1
  vector<int> markers(surface.size(), 1);
  CHKERR m_field.get_moab().tag_set_data(th_marker, surface, &*markers.begin());
  // Mark all side sets
  int shift = 1;
  map<int, int> id_shift_map; // each meshset has set unique bit
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(
           (*cOre.getInterface<MeshsetsManager>()), SIDESET, it)) {
    int ms_id = it->getMeshsetId();
    id_shift_map[ms_id] = 1 << shift; // shift bit
    ++shift;
    Range sideset_faces;
    CHKERR m_field.getInterface<MeshsetsManager>()->getEntitiesByDimension(
        ms_id, SIDESET, 2, sideset_faces, true);
    markers.resize(sideset_faces.size());
    CHKERR m_field.get_moab().tag_get_data(th_marker, sideset_faces,
                                           &*markers.begin());
    for (unsigned int ii = 0; ii < markers.size(); ii++) {
      markers[ii] |= id_shift_map[ms_id]; // add bit to marker
    }
    CHKERR m_field.get_moab().tag_set_data(th_marker, sideset_faces,
                                           &*markers.begin());
  }
  Range nodes_to_remove; // none
  markers.resize(nodes_to_remove.size());
  fill(markers.begin(), markers.end(), -1);
  CHKERR m_field.get_moab().tag_set_data(th_marker, nodes_to_remove,
                                         &*markers.begin());

  // nodes
  if (tetGenData.size() == 1) {

    Range ents_to_tetgen = surf_ents.sVols;
    CHKERR m_field.get_moab().get_connectivity(surf_ents.sVols, ents_to_tetgen,
                                               true);
    for (int dd = 2; dd >= 1; dd--) {
      CHKERR m_field.get_moab().get_adjacencies(
          surf_ents.sVols, dd, false, ents_to_tetgen, moab::Interface::UNION);
    }

    // Load mesh to TetGen data structures
    CHKERR tetgen_iface->inData(ents_to_tetgen, in, moabTetGenMap,
                                tetGenMoabMap, th);
    CHKERR tetgen_iface->setGeomData(in, moabTetGenMap, tetGenMoabMap,
                                     types_ents);
    std::vector<pair<Range, int> > markers;
    for (Range::iterator tit = surface.begin(); tit != surface.end(); tit++) {
      Range facet;
      facet.insert(*tit);
      markers.push_back(pair<Range, int>(facet, 2));
    }
    for (Range::iterator tit = surf_ents.vSkinWithoutBodySkin.begin();
         tit != surf_ents.vSkinWithoutBodySkin.end(); tit++) {
      Range facet;
      facet.insert(*tit);
      markers.push_back(pair<Range, int>(facet, 1));
    }
    Range other_facets;
    other_facets = subtract(surf_ents.vSkin, surf_ents.vSkinWithoutBodySkin);
    for (Range::iterator tit = other_facets.begin(); tit != other_facets.end();
         tit++) {
      Range facet;
      facet.insert(*tit);
      markers.push_back(pair<Range, int>(facet, 0));
    }
    CHKERR tetgen_iface->setFaceData(markers, in, moabTetGenMap, tetGenMoabMap);
  }

  if (debug) {
    if (m_field.get_comm_rank() == 0) {
      char tetgen_in_file_name[] = "in";
      in.save_nodes(tetgen_in_file_name);
      in.save_elements(tetgen_in_file_name);
      in.save_faces(tetgen_in_file_name);
      in.save_edges(tetgen_in_file_name);
      in.save_poly(tetgen_in_file_name);
    }
  }

  // generate new mesh
  {
    vector<string>::iterator sw = switches.begin();
    for (int ii = 0; sw != switches.end(); sw++, ii++) {
      tetgenio &_in_ = tetGenData.back();
      tetGenData.push_back(new tetgenio);
      tetgenio &_out_ = tetGenData.back();
      char *s = const_cast<char *>(sw->c_str());
      CHKERR tetgen_iface->tetRahedralize(s, _in_, _out_);
    }
  }
  tetgenio &out = tetGenData.back();
  // save elems
  if (debug) {
    char tetgen_out_file_name[] = "out";
    out.save_nodes(tetgen_out_file_name);
    out.save_elements(tetgen_out_file_name);
    out.save_faces(tetgen_out_file_name);
    out.save_edges(tetgen_out_file_name);
    out.save_poly(tetgen_out_file_name);
  }

  CHKERR tetgen_iface->outData(in, out, moabTetGenMap, tetGenMoabMap, bit,
                               false, false);

  Range rest_of_ents = subtract(bit_ents.mTets, surf_ents.tVols);
  for (int dd = 2; dd >= 0; dd--) {
    CHKERR moab.get_adjacencies(rest_of_ents.subset_by_dimension(3), dd, false,
                                rest_of_ents, moab::Interface::UNION);
  }
  CHKERR m_field.getInterface<BitRefManager>()->addBitRefLevel(rest_of_ents,
                                                               bit);

  Range tetgen_faces;
  map<int, Range> face_markers_map;
  CHKERR tetgen_iface->getTriangleMarkers(tetGenMoabMap, out, &tetgen_faces,
                                          &face_markers_map);

  tetgenSurfaces = face_markers_map[1];
  for (map<int, Range>::iterator mit = face_markers_map.begin();
       mit != face_markers_map.end(); mit++) {
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(
             (*cOre.getInterface<MeshsetsManager>()), SIDESET, it)) {
      int msId = it->getMeshsetId();
      if (id_shift_map[msId] & mit->first) {
        EntityHandle meshset = it->getMeshset();
        CHKERR m_field.get_moab().add_entities(
            meshset, mit->second.subset_by_type(MBTRI));
      }
    }
  }

  MoFEMFunctionReturn(0);
}

#endif // WITH_TETGEN

MoFEMErrorCode CutMeshInterface::setTagData(Tag th,const BitRefLevel bit) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  Range nodes;
  if(bit.none()) {
    CHKERR moab.get_entities_by_type(0, MBVERTEX, nodes);
  } else {
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit,BitRefLevel().set(),MBVERTEX,nodes
    );
  }
  std::vector<double> coords(3 * nodes.size());
  CHKERR moab.get_coords(nodes, &coords[0]);
  CHKERR moab.tag_set_data(th, nodes, &coords[0]);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::setCoords(Tag th, const BitRefLevel bit,
                                           const BitRefLevel mask) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  Range nodes;
   if(bit.none()) {
    CHKERR moab.get_entities_by_type(0, MBVERTEX, nodes);
  } else {
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit, mask, MBVERTEX, nodes);
  }
  std::vector<double> coords(3 * nodes.size());
  CHKERR moab.tag_get_data(th, nodes, &coords[0]);
  CHKERR moab.set_coords(nodes, &coords[0]);
  MoFEMFunctionReturn(0);
}

struct SaveData {
  moab::Interface &moab;
  SaveData(moab::Interface &moab) : moab(moab) {}
  MoFEMErrorCode operator()(const std::string name, const Range &ents) {
    MoFEMFunctionBeginHot;
    EntityHandle meshset;
    rval = moab.create_meshset(MESHSET_SET, meshset);
    CHKERRG(rval);
    rval = moab.add_entities(meshset, ents);
    CHKERRG(rval);
    rval = moab.write_file(name.c_str(), "VTK", "", &meshset, 1);
    CHKERRG(rval);
    rval = moab.delete_entities(&meshset, 1);
    CHKERRG(rval);
    MoFEMFunctionReturnHot(0);
  }
};

MoFEMErrorCode CutMeshInterface::saveCutEdges() {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;

  ierr = SaveData(moab)("out_vol.vtk", vOlume);
  CHKERRG(ierr);
  ierr = SaveData(moab)("out_surface.vtk", sUrface);
  CHKERRG(ierr);
  ierr = SaveData(moab)("out_cut_edges.vtk", cutEdges);
  CHKERRG(ierr);
  ierr = SaveData(moab)("out_cut_volumes.vtk", cutVolumes);
  CHKERRG(ierr);
  ierr = SaveData(moab)("out_cut_new_volumes.vtk", cutNewVolumes);
  CHKERRG(ierr);
  ierr = SaveData(moab)("out_cut_new_surfaces.vtk", cutNewSurfaces);
  CHKERRG(ierr);

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CutMeshInterface::saveTrimEdges() {
  moab::Interface &moab = cOre.getInterface<CoreInterface>()->get_moab();
  MoFEMFunctionBeginHot;

  ierr = SaveData(moab)("out_trim_new_volumes.vtk", trimNewVolumes);
  CHKERRG(ierr);
  ierr = SaveData(moab)("out_trim_new_surfaces.vtk", trimNewSurfaces);
  CHKERRG(ierr);
  ierr = SaveData(moab)("out_trim_edges.vtk", trimEdges);
  CHKERRG(ierr);
 
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM
