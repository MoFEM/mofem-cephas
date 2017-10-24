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

PetscErrorCode
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

CutMeshInterface::CutMeshInterface(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)) {
  lineSearchSteps = 10;
  nbMaxMergingCycles = 200;
  nbMaxTrimSearchIterations = 20;
}

PetscErrorCode CutMeshInterface::setSurface(const Range &surface) {
  MoFEMFunctionBeginHot;
  sUrface = surface;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::copySurface(const Range &surface, Tag th,
                                             double *shift, double *origin,
                                             double *transform) {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  for (Range::const_iterator tit = surface.begin(); tit != surface.end();
       tit++) {
    int num_nodes;
    const EntityHandle *conn;
    rval = moab.get_connectivity(*tit, conn, num_nodes, true);
    CHKERRQ_MOAB(rval);
    MatrixDouble coords(num_nodes, 3);
    if (th) {
      rval = moab.tag_get_data(th, conn, num_nodes, &coords(0, 0));
      CHKERRQ_MOAB(rval);
    } else {
      rval = moab.get_coords(conn, num_nodes, &coords(0, 0));
      CHKERRQ_MOAB(rval);
    }
    EntityHandle new_verts[num_nodes];
    for (int nn = 0; nn != num_nodes; nn++) {
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
        VectorAdaptor vec_shift(3,
                                ublas::shallow_array_adaptor<double>(3, shift));
        mr = mr + vec_shift;
      }
      rval = moab.create_vertex(&coords(nn, 0), new_verts[nn]);
      CHKERRQ_MOAB(rval);
    }
    EntityHandle ele;
    rval = moab.create_element(MBTRI, new_verts, num_nodes, ele);
    CHKERRQ_MOAB(rval);
    sUrface.insert(ele);
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::setVolume(const Range &volume) {
  MoFEMFunctionBeginHot;
  vOlume = volume;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::mergeSurface(const Range &surface) {
  MoFEMFunctionBeginHot;
  sUrface.merge(surface);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::mergeVolumes(const Range &volume) {
  MoFEMFunctionBeginHot;
  vOlume.merge(volume);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::buildTree() {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  treeSurfPtr = boost::shared_ptr<OrientedBoxTreeTool>(
      new OrientedBoxTreeTool(&moab, "ROOTSETSURF", true));
  rval = treeSurfPtr->build(sUrface, rootSetSurf);
  CHKERRQ_MOAB(rval);
  Range faces;
  rval = moab.get_adjacencies(vOlume, 2, false, faces, moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  MoFEMFunctionReturnHot(0);
}

struct UpdateMeshsets {
  PetscErrorCode operator()(MoFEM::Core &core, const BitRefLevel &bit) const {
    MoFEMFunctionBeginHot;
    ierr = core.getInterface<MeshsetsManager>()
               ->updateAllMeshsetsByEntitiesChildren(bit);
    CHKERRQ(ierr);

    MoFEMFunctionReturnHot(0);
  }
};

PetscErrorCode CutMeshInterface::cutAndTrim(
    const BitRefLevel &bit_level1, const BitRefLevel &bit_level2, Tag th,
    const double tol_cut, const double tol_cut_close, const double tol_trim,
    const double tol_trim_close, Range *fixed_edges, Range *corner_nodes,
    const bool update_meshsets) {

  MoFEMFunctionBeginHot;
  // cut mesh
  ierr = findEdgesToCut(tol_cut);
  CHKERRQ(ierr);
  ierr = getEntsOnCutSurface(tol_cut_close);
  CHKERRQ(ierr);
  ierr = cutEdgesInMiddle(bit_level1);
  CHKERRQ(ierr);
  if (fixed_edges) {
    ierr = cOre.getInterface<BitRefManager>()->updateRange(*fixed_edges,
                                                           *fixed_edges);
    CHKERRQ(ierr);
  }
  if (corner_nodes) {
    ierr = cOre.getInterface<BitRefManager>()->updateRange(*corner_nodes,
                                                           *corner_nodes);
    CHKERRQ(ierr);
  }
  if (update_meshsets) {
    ierr = UpdateMeshsets()(cOre, bit_level1);
    CHKERRQ(ierr);
  }
  ierr = moveMidNodesOnCutEdges(th);
  CHKERRQ(ierr);

  if(debug) {
    ierr = cOre.getInterface<BitRefManager>()->writeBitLevelByType(
        bit_level1, BitRefLevel().set(), MBTET, "out_tets_cut.vtk", "VTK", "");
    CHKERRQ(ierr);
  }

  // trim mesh
  ierr = findEdgesToTrim(th, tol_trim);
  CHKERRQ(ierr);
  ierr = trimEdgesInTheMiddle(bit_level2, th, tol_trim_close);
  CHKERRQ(ierr);
  if (fixed_edges) {
    ierr = cOre.getInterface<BitRefManager>()->updateRange(*fixed_edges,
                                                           *fixed_edges);
    CHKERRQ(ierr);
  }
  if (corner_nodes) {
    ierr = cOre.getInterface<BitRefManager>()->updateRange(*corner_nodes,
                                                           *corner_nodes);
    CHKERRQ(ierr);
  }
  if (update_meshsets) {
    ierr = UpdateMeshsets()(cOre, bit_level2);
    CHKERRQ(ierr);
  }
  ierr = moveMidNodesOnTrimedEdges(th);
  CHKERRQ(ierr);

  if(debug) {
    ierr = cOre.getInterface<BitRefManager>()->writeBitLevelByType(
        bit_level2, BitRefLevel().set(), MBTET, "out_tets_trim.vtk", "VTK", "");
    CHKERRQ(ierr);
  }

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::cutTrimAndMerge(
    const int fraction_level, const BitRefLevel &bit_level1,
    const BitRefLevel &bit_level2, const BitRefLevel &bit_level3, Tag th,
    const double tol_cut, const double tol_cut_close, const double tol_trim,
    const double tol_trim_close, Range &fixed_edges, Range &corner_nodes,
    const bool update_meshsets, const bool debug) {
  MoFEMFunctionBeginHot;
  if(debug) {
    ierr = cOre.getInterface<BitRefManager>()->writeEntitiesNotInDatabse(
        "ents_not_in_databse.vtk", "VTK", "");
    CHKERRQ(ierr);
  }
  ierr =
      cutAndTrim(bit_level1, bit_level2, th, tol_cut, tol_cut_close, tol_trim,
                 tol_trim_close, &fixed_edges, &corner_nodes, update_meshsets);
  CHKERRQ(ierr);
  if(debug) {
    ierr = cOre.getInterface<BitRefManager>()->writeEntitiesNotInDatabse(
        "cut_trim_ents_not_in_databse.vtk", "VTK", "");
    CHKERRQ(ierr);
  }

  ierr =
      mergeBadEdges(fraction_level, bit_level2, bit_level1, bit_level3,
                    getNewTrimSurfaces(), fixed_edges, corner_nodes, th, debug);
  CHKERRQ(ierr);
  ierr = removePathologicalFrontTris(bit_level3,
                                     const_cast<Range &>(getMergedSurfaces()));
  CHKERRQ(ierr);

  if(debug) {
    ierr = cOre.getInterface<BitRefManager>()->writeBitLevelByType(
        bit_level3, BitRefLevel().set(), MBTET, "out_tets_merged.vtk", "VTK",
        "");
    CHKERRQ(ierr);
    ierr = cOre.getInterface<BitRefManager>()->writeEntitiesNotInDatabse(
        "cut_trim_merge_ents_not_in_databse.vtk", "VTK", "");
    CHKERRQ(ierr);
  }

  ierr =
      cOre.getInterface<BitRefManager>()->updateRange(fixed_edges, fixed_edges);
  CHKERRQ(ierr);
  ierr = cOre.getInterface<BitRefManager>()->updateRange(corner_nodes,
                                                         corner_nodes);
  CHKERRQ(ierr);
  if (update_meshsets) {
    ierr = UpdateMeshsets()(cOre, bit_level3);
    CHKERRQ(ierr);
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::findEdgesToCut(const double low_tol,
                                                int verb) {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  // Range vol_edges;
  // rval = moab.get_adjacencies(
  //   vOlume,1,true,vol_edges,moab::Interface::UNION
  // ); CHKERRQ_MOAB(rval);
  edgesToCut.clear();
  cutEdges.clear();
  double ray_length;
  double ray_point[3], unit_ray_dir[3];
  VectorAdaptor vec_unit_ray_dir(
      3, ublas::shallow_array_adaptor<double>(3, unit_ray_dir));
  VectorAdaptor vec_ray_point(
      3, ublas::shallow_array_adaptor<double>(3, ray_point));

  Tag th_dist;
  double def_val[] = {0};
  rval = moab.tag_get_handle("DIST", 1, MB_TYPE_DOUBLE, th_dist,
                             MB_TAG_CREAT | MB_TAG_SPARSE, def_val);
  CHKERRQ_MOAB(rval);
  Tag th_dist_normal;
  rval = moab.tag_get_handle("DIST_NORMAL", 1, MB_TYPE_DOUBLE, th_dist_normal,
                             MB_TAG_CREAT | MB_TAG_SPARSE, def_val);
  CHKERRQ_MOAB(rval);

  Range vol_vertices;
  rval = moab.get_connectivity(vOlume, vol_vertices, true);
  CHKERRQ_MOAB(rval);
  for (Range::iterator vit = vol_vertices.begin(); vit != vol_vertices.end();
       vit++) {
    double coords[3];
    rval = moab.get_coords(&*vit, 1, coords);
    CHKERRQ_MOAB(rval);
    VectorAdaptor point_in(3, ublas::shallow_array_adaptor<double>(3, coords));
    double p_out[3];
    EntityHandle facets_out;
    rval = treeSurfPtr->closest_to_location(&coords[0], rootSetSurf, p_out,
                                            facets_out);
    CHKERRQ_MOAB(rval);
    VectorAdaptor point_out(3, ublas::shallow_array_adaptor<double>(3, p_out));
    double normal[3];
    Util::normal(&moab, facets_out, normal[0], normal[1], normal[2]);
    VectorAdaptor n(3, ublas::shallow_array_adaptor<double>(3, normal));
    VectorDouble3 delta = point_out - point_in;
    double dist = norm_2(delta);
    double dist_normal = inner_prod(delta, n) / norm_2(n);
    rval = moab.tag_set_data(th_dist, &*vit, 1, &dist);
    CHKERRQ_MOAB(rval);
    rval = moab.tag_set_data(th_dist_normal, &*vit, 1, &dist_normal);
    CHKERRQ_MOAB(rval);
  }

  Range vol_edges;
  rval =
      moab.get_adjacencies(vOlume, 1, true, vol_edges, moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  aveLength = 0;
  int nb_ave_length = 0;
  for (Range::iterator eit = vol_edges.begin(); eit != vol_edges.end(); eit++) {
    int num_nodes;
    const EntityHandle *conn;
    rval = moab.get_connectivity(*eit, conn, num_nodes, true);
    CHKERRQ_MOAB(rval);
    double dist[num_nodes];
    rval = moab.tag_get_data(th_dist, conn, num_nodes, dist);
    CHKERRQ_MOAB(rval);
    double dist_normal[num_nodes];
    rval = moab.tag_get_data(th_dist_normal, conn, num_nodes, dist_normal);
    CHKERRQ_MOAB(rval);
    ierr = getRayForEdge(*eit, vec_ray_point, vec_unit_ray_dir, ray_length);
    CHKERRQ(ierr);
    const double tol = ray_length * low_tol;
    if ((dist_normal[0] * dist_normal[1] < 0) ||
        (dist_normal[0] * dist_normal[1] == 0 &&
         (dist_normal[0] + dist_normal[1]) > 0)) {
      std::vector<double> distances_out;
      std::vector<EntityHandle> facets_out;
      rval = treeSurfPtr->ray_intersect_triangles(distances_out, facets_out,
                                                  rootSetSurf, tol, ray_point,
                                                  unit_ray_dir, &ray_length);
      CHKERRQ_MOAB(rval);
      if (!distances_out.empty()) {
        aveLength += ray_length;
        nb_ave_length++;
        edgesToCut[*eit].dIst = distances_out[0];
        edgesToCut[*eit].lEngth = ray_length;
        edgesToCut[*eit].unitRayDir = vec_unit_ray_dir;
        edgesToCut[*eit].rayPoint = vec_ray_point;
        cutEdges.insert(*eit);
      }
    }
  }
  aveLength /= nb_ave_length;

  // tak all volumes adjacent to cut edges
  cutVolumes.clear();
  rval = moab.get_adjacencies(cutEdges, 3, false, cutVolumes,
                              moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  // get edges on the cut volumes
  Range edges;
  rval =
      moab.get_adjacencies(cutVolumes, 1, false, edges, moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  edges = subtract(edges, cutEdges);

  // add to cut set edges which are cutted by extension of cuttting surfeca
  for (Range::iterator eit = edges.begin(); eit != edges.end(); eit++) {
    int num_nodes;
    const EntityHandle *conn;
    rval = moab.get_connectivity(*eit, conn, num_nodes, true);
    CHKERRQ_MOAB(rval);
    double dist[num_nodes];
    rval = moab.tag_get_data(th_dist, conn, num_nodes, dist);
    CHKERRQ_MOAB(rval);
    double dist_normal[num_nodes];
    rval = moab.tag_get_data(th_dist_normal, conn, num_nodes, dist_normal);
    CHKERRQ_MOAB(rval);
    if (dist_normal[0] * dist_normal[1] < 0 ||
        (dist_normal[0] * dist_normal[1] == 0 &&
         (dist_normal[0] + dist_normal[1]) > 0)) {
      ierr = getRayForEdge(*eit, vec_ray_point, vec_unit_ray_dir, ray_length);
      CHKERRQ(ierr);
      double s =
          fabs(dist_normal[0]) / (fabs(dist_normal[0]) + fabs(dist_normal[1]));
      edgesToCut[*eit].dIst = s * ray_length;
      edgesToCut[*eit].lEngth = ray_length;
      edgesToCut[*eit].unitRayDir = vec_unit_ray_dir;
      edgesToCut[*eit].rayPoint = vec_ray_point;
      cutEdges.insert(*eit);
    }
  }

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::getZeroDistanceEnts(const double low_tol,
                                                     int verb) {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  Skinner skin(&moab);
  Range tets_skin;
  rval = skin.find_skin(0, vOlume, false, tets_skin);
  Range tets_skin_verts;
  rval = moab.get_connectivity(tets_skin, tets_skin_verts, true);
  CHKERRQ_MOAB(rval);
  Tag th_dist;
  rval = moab.tag_get_handle("DIST_NORMAL", th_dist);
  CHKERRQ_MOAB(rval);
  Range cut_edge_verts;
  rval = moab.get_connectivity(cutEdges, cut_edge_verts, false);
  CHKERRQ_MOAB(rval);
  // get faces and edges
  Range cut_edges_faces;
  rval = moab.get_adjacencies(cut_edge_verts, 1, true, cut_edges_faces,
                              moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  rval = moab.get_adjacencies(cut_edge_verts, 2, true, cut_edges_faces,
                              moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  zeroDistanseEnts.clear();
  verticecOnCutEdges.clear();
  for (Range::iterator fit = cut_edges_faces.begin();
       fit != cut_edges_faces.end(); fit++) {
    int num_nodes;
    const EntityHandle *conn;
    rval = moab.get_connectivity(*fit, conn, num_nodes, true);
    CHKERRQ_MOAB(rval);
    double dist[] = {0, 0, 0};
    rval = moab.tag_get_data(th_dist, conn, num_nodes, dist);
    CHKERRQ_MOAB(rval);
    if (fabs(dist[0]) < low_tol * aveLength &&
        fabs(dist[1]) < low_tol * aveLength &&
        fabs(dist[2]) < low_tol * aveLength) {
      zeroDistanseEnts.insert(*fit);
      Range adj_edges;
      rval = moab.get_adjacencies(conn, num_nodes, 1, false, adj_edges,
                                  moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      for (Range::iterator eit = adj_edges.begin(); eit != adj_edges.end();
           eit++) {
        cutEdges.erase(*eit);
        edgesToCut.erase(*eit);
      }
      double coords[9];
      rval = moab.get_coords(conn, num_nodes, coords);
      CHKERRQ(ierr);
      for (int nn = 0; nn != num_nodes; nn++) {
        if (tets_skin_verts.find(conn[nn]) == tets_skin_verts.end()) {
        VectorAdaptor s0(
            3, ublas::shallow_array_adaptor<double>(3, &coords[3 * nn]));
        double p_out[3];
        EntityHandle facets_out;
        rval = treeSurfPtr->closest_to_location(&s0[0], rootSetSurf, p_out,
                                                facets_out);
        CHKERRQ_MOAB(rval);
          VectorAdaptor point_out(
              3, ublas::shallow_array_adaptor<double>(3, p_out));
        VectorDouble3 ray = point_out - s0;
        double dist0 = norm_2(ray);
        verticecOnCutEdges[conn[nn]].dIst = dist0;
        verticecOnCutEdges[conn[nn]].lEngth = dist0;
          verticecOnCutEdges[conn[nn]].unitRayDir =
              dist0 > 0 ? ray / dist0 : ray;
        verticecOnCutEdges[conn[nn]].rayPoint = s0;
      }
    }
  }
  }
  rval = moab.tag_get_handle("DIST", th_dist);
  for (Range::iterator vit = cut_edge_verts.begin();
       vit != cut_edge_verts.end(); vit++) {
    double dist[] = {0};
    rval = moab.tag_get_data(th_dist, &*vit, 1, dist);
    CHKERRQ_MOAB(rval);
    if (fabs(dist[0]) < low_tol * aveLength) {
      zeroDistanseVerts.insert(*vit);
      Range adj_edges;
      rval = moab.get_adjacencies(&*vit, 1, 1, false, adj_edges);
      CHKERRQ_MOAB(rval);
      for (Range::iterator eit = adj_edges.begin(); eit != adj_edges.end();
           eit++) {
        cutEdges.erase(*eit);
        edgesToCut.erase(*eit);
      }
      double coords[3];
      rval = moab.get_coords(&*vit, 1, coords);
      CHKERRQ(ierr);
      VectorAdaptor s0(3, ublas::shallow_array_adaptor<double>(3, &coords[0]));
      double p_out[3];
      EntityHandle facets_out;
      rval = treeSurfPtr->closest_to_location(&s0[0], rootSetSurf, p_out,
                                              facets_out);
      CHKERRQ_MOAB(rval);
      VectorAdaptor point_out(3,
                              ublas::shallow_array_adaptor<double>(3, p_out));
      VectorDouble3 ray = point_out - s0;
      double dist0 = norm_2(ray);
      verticecOnCutEdges[*vit].dIst = dist0;
      verticecOnCutEdges[*vit].lEngth = dist0;
      verticecOnCutEdges[*vit].unitRayDir = dist0 > 0 ? ray / dist0 : ray;
      verticecOnCutEdges[*vit].rayPoint = s0;
    }
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::cutEdgesInMiddle(const BitRefLevel bit) {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MeshRefinement *refiner;
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.getInterface(refiner);
  CHKERRQ(ierr);
  ierr = m_field.get_ref_ents(&ref_ents_ptr);
  CHKERRQ(ierr);
  ierr = refiner->add_verices_in_the_middel_of_edges(cutEdges, bit);
  CHKERRQ(ierr);
  ierr = refiner->refine_TET(vOlume, bit, false);
  CHKERRQ(ierr);
  // Tag th_ray_dir;
  // double def_val[] = {0,0,0};
  // rval = moab.tag_get_handle(
  //   "RAY_DIR",3,MB_TYPE_DOUBLE,th_ray_dir,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  // ); CHKERRQ_MOAB(rval);
  cutNewVolumes.clear();
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, BitRefLevel().set(), MBTET, cutNewVolumes);
  CHKERRQ(ierr);
  cutNewSurfaces.clear();
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, bit, MBTRI, cutNewSurfaces);
  CHKERRQ(ierr);
  // Find new vertices on catted edges
  cutNewVertices.clear();
  rval = moab.get_connectivity(zeroDistanseEnts, cutNewVertices, true);
  CHKERRQ_MOAB(rval);
  cutNewVertices.merge(zeroDistanseVerts);
  for (map<EntityHandle, TreeData>::iterator mit = edgesToCut.begin();
       mit != edgesToCut.end(); mit++) {
    boost::shared_ptr<RefEntity> ref_ent =
        *(ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>().find(
            boost::make_tuple(mit->first, MBVERTEX)));
    if ((ref_ent->getBitRefLevel() & bit).any()) {
      EntityHandle vert = ref_ent->getRefEnt();
      cutNewVertices.insert(vert);
      verticecOnCutEdges[vert] = mit->second;
      // rval =
      // moab.tag_set_data(th_ray_dir,&vert,1,&mit->second.unitRayDir[0]);
      // CHKERRQ_MOAB(rval);
    }
  }
  // Add zero distance entities faces
  Range tets_skin;
  Skinner skin(&moab);
  rval = skin.find_skin(0, cutNewVolumes, false, tets_skin);
  CHKERRQ_MOAB(rval);
  cutNewSurfaces.merge(
      subtract(zeroDistanseEnts.subset_by_type(MBTRI), tets_skin));
  // At that point cutNewSurfaces has all newly created faces, now take all
  // nodes on those faces and subtract nodes on catted edges. Faces adjacent to
  // nodes which left are not part of surface.
  Range diff_verts;
  rval = moab.get_connectivity(cutNewSurfaces, diff_verts, true);
  CHKERRQ_MOAB(rval);
  diff_verts = subtract(diff_verts, cutNewVertices);
  Range subtract_faces;
  rval = moab.get_adjacencies(diff_verts, 2, false, subtract_faces,
                              moab::Interface::UNION);
  CHKERRQ_MOAB(rval);
  cutNewSurfaces = subtract(cutNewSurfaces, subtract_faces);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::trimEdgesInTheMiddle(const BitRefLevel bit,
                                                      Tag th,
                                                      const double tol) {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MeshRefinement *refiner;
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBeginHot;

  ierr = m_field.getInterface(refiner);
  CHKERRQ(ierr);
  ierr = m_field.get_ref_ents(&ref_ents_ptr);
  CHKERRQ(ierr);
  ierr = refiner->add_verices_in_the_middel_of_edges(trimEdges, bit);
  CHKERRQ(ierr);
  ierr = refiner->refine_TET(cutNewVolumes, bit, false);
  CHKERRQ(ierr);
  trimNewVolumes.clear();
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, bit, MBTET, trimNewVolumes);
  CHKERRQ(ierr);
  // Get vertices which are on trim edges
  verticecOnTrimEdges.clear();
  trimNewVertices.clear();
  for (map<EntityHandle, TreeData>::iterator mit = edgesToTrim.begin();
       mit != edgesToTrim.end(); mit++) {
    boost::shared_ptr<RefEntity> ref_ent =
        *(ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>().find(
            boost::make_tuple(mit->first, MBVERTEX)));
    if ((ref_ent->getBitRefLevel() & bit).any()) {
      EntityHandle vert = ref_ent->getRefEnt();
      trimNewVertices.insert(vert);
      verticecOnTrimEdges[vert] = mit->second;
    }
  }

  // Get faces which are trimmed
  trimNewSurfaces.clear();
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, bit, MBTRI, trimNewSurfaces);
  CHKERRQ(ierr);
  Range trim_new_surfaces_nodes;
  rval = moab.get_connectivity(trimNewSurfaces, trim_new_surfaces_nodes, true);
  CHKERRQ_MOAB(rval);
  trim_new_surfaces_nodes = subtract(trim_new_surfaces_nodes, trimNewVertices);
  trim_new_surfaces_nodes = subtract(trim_new_surfaces_nodes, cutNewVertices);
  Range faces_not_on_surface;
  rval = moab.get_adjacencies(trim_new_surfaces_nodes, 2, false,
                              faces_not_on_surface, moab::Interface::UNION);
  trimNewSurfaces = subtract(trimNewSurfaces, faces_not_on_surface);

  // Get surfaces which are not trimmed and add them to surface
  Range all_surfaces_on_bit_level;
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, BitRefLevel().set(), MBTRI, all_surfaces_on_bit_level);
  CHKERRQ(ierr);
  all_surfaces_on_bit_level =
      intersect(all_surfaces_on_bit_level, cutNewSurfaces);
  trimNewSurfaces = unite(trimNewSurfaces, all_surfaces_on_bit_level);

  Range check_verts;
  rval = moab.get_connectivity(trimNewSurfaces, check_verts, true);
  CHKERRQ_MOAB(rval);
  check_verts = subtract(check_verts, trimNewVertices);
  for (Range::iterator vit = check_verts.begin(); vit != check_verts.end();
       vit++) {
    double coords[3];
    if (th) {
      rval = moab.tag_get_data(th, &*vit, 1, coords);
      CHKERRQ_MOAB(rval);
    } else {
      rval = moab.get_coords(&*vit, 1, coords);
      CHKERRQ_MOAB(rval);
    }
    double point_out[3];
    EntityHandle facets_out;
    rval = treeSurfPtr->closest_to_location(coords, rootSetSurf, point_out,
                                            facets_out);
    CHKERRQ_MOAB(rval);
    VectorAdaptor s(3, ublas::shallow_array_adaptor<double>(3, coords));
    VectorAdaptor p(3, ublas::shallow_array_adaptor<double>(3, point_out));
    if (norm_2(s - p) / aveLength > tol) {
      Range adj;
      rval = moab.get_adjacencies(&*vit, 1, 2, false, adj);
      CHKERRQ_MOAB(rval);
      trimNewSurfaces = subtract(trimNewSurfaces, adj);
    }
  }

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::moveMidNodesOnCutEdges(Tag th) {
  MoFEMFunctionBeginHot;

  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;

  // Tag th_ray_dir;
  // double def_val[] = {0,0,0};
  // rval = moab.tag_get_handle(
  //   "RAY_DIR_2",3,MB_TYPE_DOUBLE,th_ray_dir,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  // ); CHKERRQ_MOAB(rval);

  // Range out_side_vertices;
  for (map<EntityHandle, TreeData>::iterator mit = verticecOnCutEdges.begin();
       mit != verticecOnCutEdges.end(); mit++) {
    double dist = mit->second.dIst;
    // cout << s << " " << mit->second.dIst << " " << mit->second.lEngth <<
    // endl;
    VectorDouble3 new_coors =
        mit->second.rayPoint + dist * mit->second.unitRayDir;
    if (th) {
      rval = moab.tag_set_data(th, &mit->first, 1, &new_coors[0]);
      CHKERRQ_MOAB(rval);
    } else {
      rval = moab.set_coords(&mit->first, 1, &new_coors[0]);
      CHKERRQ_MOAB(rval);
    }
    // rval =
    // moab.tag_set_data(th_ray_dir,&mit->first,1,&mit->second.unitRayDir[0]);
    // CHKERRQ_MOAB(rval);
  }

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::moveMidNodesOnTrimedEdges(Tag th) {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  // Range out_side_vertices;
  for (map<EntityHandle, TreeData>::iterator mit = verticecOnTrimEdges.begin();
       mit != verticecOnTrimEdges.end(); mit++) {
    double dist = mit->second.dIst;
    // cout << s << " " << mit->second.dIst << " " << mit->second.lEngth <<
    // endl;
    VectorDouble3 new_coors =
        mit->second.rayPoint + dist * mit->second.unitRayDir;
    if (th) {
      rval = moab.tag_set_data(th, &mit->first, 1, &new_coors[0]);
      CHKERRQ_MOAB(rval);
    } else {
      rval = moab.set_coords(&mit->first, 1, &new_coors[0]);
      CHKERRQ_MOAB(rval);
    }
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::findEdgesToTrim(Tag th, const double tol,
                                                 int verb) {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;

  // takes edges on body skin
  Skinner skin(&moab);
  Range tets_skin;
  rval = skin.find_skin(0, cutNewVolumes, false, tets_skin);
  CHKERRQ_MOAB(rval);
  Range adj_edges_tets_skin;
  rval = moab.get_adjacencies(tets_skin, 1, false, adj_edges_tets_skin,
                              moab::Interface::UNION);
  CHKERRQ_MOAB(rval);

  // get edges on new surface
  Range edges;
  rval = moab.get_adjacencies(cutNewSurfaces, 1, false, edges,
                              moab::Interface::UNION);
  CHKERRQ_MOAB(rval);

  // clear data ranges
  trimEdges.clear();
  edgesToTrim.clear();

  struct CloasestPointProjection {

    boost::shared_ptr<OrientedBoxTreeTool> treeSurfPtr;
    EntityHandle rootSetSurf;
    VectorDouble3 S0;
    VectorDouble3 rAY;
    double aveLength;
    double pointOut[3];
    VectorAdaptor vecPointOut;
    CloasestPointProjection(boost::shared_ptr<OrientedBoxTreeTool> &tree,
                            EntityHandle root_set, VectorDouble3 &s0,
                            VectorDouble3 &ray, double ave_length)
        : treeSurfPtr(tree), rootSetSurf(root_set), S0(s0), rAY(ray),
          aveLength(ave_length),
          vecPointOut(3,
                      ublas::shallow_array_adaptor<double>(3, &pointOut[0])) {}

    VectorDouble3 &operator()(const int max_it, const double tol) {
      MoFEMFunctionBeginHot;
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
        if (s / aveLength < tol)
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
    rval = moab.get_connectivity(*eit, conn, num_nodes, true);
    CHKERRQ_MOAB(rval);
    double coords[3 * num_nodes];
    if (th) {
      rval = moab.tag_get_data(th, conn, num_nodes, coords);
      CHKERRQ_MOAB(rval);
    } else {
      rval = moab.get_coords(conn, num_nodes, coords);
      CHKERRQ_MOAB(rval);
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
    rval = treeSurfPtr->closest_to_location(&coords[0], rootSetSurf, point_out0,
                                            facets_out0);
    CHKERRQ_MOAB(rval);
    // find closest point on the surface from second node
    double point_out1[3];
    EntityHandle facets_out1;
    rval = treeSurfPtr->closest_to_location(&coords[3], rootSetSurf, point_out1,
                                            facets_out1);
    CHKERRQ_MOAB(rval);
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
    // If one of nodes is on the surface and other is not, that edge is to trim
    if (min_dist / length < tol && max_dist / length > tol) {
      // add ege to trim
      double dist;
      VectorDouble3 ray;
      VectorDouble3 trimed_end;
      VectorDouble3 itersection_point;
      if (max_dist == dist0) {
        // move mid node in reference to node 0
        trimed_end = s0;
        ray = s1 - trimed_end;
        itersection_point =
            CloasestPointProjection(treeSurfPtr, rootSetSurf, trimed_end, ray,
                                    aveLength)(nbMaxTrimSearchIterations, tol);
        ray = itersection_point - trimed_end;
        dist = norm_2(ray);
      } else {
        // move node in reference to node 1
        trimed_end = s1;
        ray = s0 - trimed_end;
        itersection_point =
            CloasestPointProjection(treeSurfPtr, rootSetSurf, trimed_end, ray,
                                    aveLength)(nbMaxTrimSearchIterations, tol);
        ray = itersection_point - trimed_end;
        dist = norm_2(ray);
      }
      if (fabs(dist - length) / length > tol) {
        edgesToTrim[*eit].dIst = dist;
        edgesToTrim[*eit].lEngth = dist;
        edgesToTrim[*eit].unitRayDir = ray / dist;
        edgesToTrim[*eit].rayPoint = trimed_end;
        trimEdges.insert(*eit);
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::getRayForEdge(const EntityHandle ent,
                                               VectorAdaptor ray_point,
                                               VectorAdaptor unit_ray_dir,
                                               double &ray_length) const {
  const MoFEM::CoreInterface &m_field = cOre;
  const moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  int num_nodes;
  const EntityHandle *conn;
  rval = moab.get_connectivity(ent, conn, num_nodes, true);
  CHKERRQ_MOAB(rval);
  double coords[6];
  rval = moab.get_coords(conn, num_nodes, coords);
  CHKERRQ_MOAB(rval);
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

PetscErrorCode
CutMeshInterface::removePathologicalFrontTris(const BitRefLevel split_bit,
                                              Range &ents) {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  PrismInterface *interface;
  MoFEMFunctionBeginHot;
  ierr = m_field.getInterface(interface);
  CHKERRQ(ierr);
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET, meshset);
  CHKERRQ_MOAB(rval);
  rval = moab.add_entities(meshset, ents);
  CHKERRQ_MOAB(rval);
  Range front_tris;
  ierr = interface->findIfTringleHasThreeNodesOnInternalSurfaceSkin(
      meshset, split_bit, true, front_tris);
  CHKERRQ(ierr);
  ents = subtract(ents, front_tris);
  rval = moab.delete_entities(&meshset, 1);
  CHKERRQ_MOAB(rval);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::splitSides(const BitRefLevel split_bit,
                                            const BitRefLevel bit,
                                            const Range &ents, Tag th) {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  PrismInterface *interface;
  MoFEMFunctionBegin;
  ierr = m_field.getInterface(interface);
  CHKERRQ(ierr);
  EntityHandle meshset_volume;
  rval = moab.create_meshset(MESHSET_SET, meshset_volume);
  CHKERRQ_MOAB(rval);
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      split_bit, BitRefLevel().set(), MBTET, meshset_volume);
  CHKERRQ(ierr);
  EntityHandle meshset_surface;
  rval = moab.create_meshset(MESHSET_SET, meshset_surface);
  CHKERRQ_MOAB(rval);
  rval = moab.add_entities(meshset_surface, ents);
  CHKERRQ_MOAB(rval);
  ierr = interface->getSides(meshset_surface, split_bit, true);
  CHKERRQ(ierr);
  ierr =
      interface->splitSides(meshset_volume, bit, meshset_surface, true, true);
  CHKERRQ(ierr);
  rval = moab.delete_entities(&meshset_volume, 1);
  CHKERRQ_MOAB(rval);
  rval = moab.delete_entities(&meshset_surface, 1);
  CHKERRQ_MOAB(rval);
  if (th) {
    Range prisms;
    ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit, BitRefLevel().set(), MBPRISM, prisms);
    CHKERRQ(ierr);
    for (Range::iterator pit = prisms.begin(); pit != prisms.end(); pit++) {
      int num_nodes;
      const EntityHandle *conn;
      rval = moab.get_connectivity(*pit, conn, num_nodes, true);
      CHKERRQ_MOAB(rval);
      MatrixDouble data(3, 3);
      rval = moab.tag_get_data(th, conn, 3, &data(0, 0));
      CHKERRQ_MOAB(rval);
      // cerr << data << endl;
      rval = moab.tag_set_data(th, &conn[3], 3, &data(0, 0));
      CHKERRQ_MOAB(rval);
    }
  }
  MoFEMFunctionReturn(0);
}

PetscErrorCode CutMeshInterface::mergeBadEdges(
    const int fraction_level, const Range &tets, const Range &surface,
    const Range &fixed_edges, const Range &corner_nodes, Range &edges_to_merge,
    Range &out_tets, Range &new_surf, Tag th, const bool update_meshsets,
    const BitRefLevel *bit_ptr,const bool debug) {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;

  /**
   * \brief Merge nodes
   */
  struct MergeNodes {
    MoFEM::CoreInterface &mField;
    const bool onlyIfImproveQuality;
    const int lineSearch;
    Tag tH;
    bool updateMehsets;

    MergeNodes(MoFEM::CoreInterface &m_field,
               const bool only_if_improve_quality, const int line_search,
               Tag th, bool update_mehsets)
        : mField(m_field), onlyIfImproveQuality(only_if_improve_quality),
          lineSearch(line_search), tH(th), updateMehsets(update_mehsets) {
      mField.getInterface(nodeMergerPtr);
    }
    NodeMergerInterface *nodeMergerPtr;
    PetscErrorCode operator()(EntityHandle father, EntityHandle mother,
                              Range &proc_tets, Range &new_surf,
                              Range &edges_to_merge, Range &not_merged_edges,
                              bool add_child = true) const {
      moab::Interface &moab = mField.get_moab();
      MoFEMFunctionBeginHot;
      const EntityHandle conn[] = {father, mother};
      Range vert_tets;
      rval = moab.get_adjacencies(conn, 2, 3, false, vert_tets,
                                  moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      vert_tets = intersect(vert_tets, proc_tets);
      Range out_tets;
      ierr = nodeMergerPtr->mergeNodes(father, mother, out_tets, &vert_tets,
                                       onlyIfImproveQuality, 0, lineSearch, tH);
      CHKERRQ(ierr);
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
        CHKERRQ(ierr);
        new_surf.merge(child_surf_ents);

        Range edges_to_merge_child_ents = intersect(edges_to_merge, child_ents);
        edges_to_merge = subtract(edges_to_merge, edges_to_merge_child_ents);
        Range merged_child_edge_ents;
        ierr = updateRangeByChilds(parent_child_map, edges_to_merge_child_ents,
                                   merged_child_edge_ents);
        CHKERRQ(ierr);

        Range not_merged_edges_child_ents =
            intersect(not_merged_edges, child_ents);
        not_merged_edges =
            subtract(not_merged_edges, not_merged_edges_child_ents);
        Range not_merged_child_edge_ents;
        ierr =
            updateRangeByChilds(parent_child_map, not_merged_edges_child_ents,
                                not_merged_child_edge_ents);
        CHKERRQ(ierr);

        merged_child_edge_ents =
            subtract(merged_child_edge_ents, not_merged_child_edge_ents);
        edges_to_merge.merge(merged_child_edge_ents);
        not_merged_edges.merge(not_merged_child_edge_ents);

        if (updateMehsets) {
          Range child_ents;
          for (_IT_CUBITMESHSETS_FOR_LOOP_(mField, cubit_it)) {
            EntityHandle cubit_meshset = cubit_it->meshset;
            Range parent_ents;
            child_ents.clear();
            rval =
                moab.get_entities_by_handle(cubit_meshset, parent_ents, true);
            CHKERRQ_MOAB(rval);
            parent_ents = intersect(parent_ents, child_ents);
            ierr =
                updateRangeByChilds(parent_child_map, parent_ents, child_ents);
            CHKERRQ(ierr);
            rval = moab.add_entities(cubit_meshset, child_ents);
            CHKERRQ_MOAB(rval);
          }
        }
      }
      MoFEMFunctionReturnHot(0);
    }

  private:
    PetscErrorCode updateRangeByChilds(
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
    MoFEM::CoreInterface &mField;
    moab::Interface &moab;
    LengthMap(MoFEM::CoreInterface &m_field, Tag th)
        : tH(th), mField(m_field), moab(m_field.get_moab()) {}

    PetscErrorCode
    operator()(const Range &tets, const Range &edges,
               std::map<double, EntityHandle> &length_map) const {
      int num_nodes;
      const EntityHandle *conn;
      double coords[6];
      MoFEMFunctionBeginHot;
      VectorAdaptor s0(3, ublas::shallow_array_adaptor<double>(3, &coords[0]));
      VectorAdaptor s1(3, ublas::shallow_array_adaptor<double>(3, &coords[3]));
      for (Range::const_iterator eit = edges.begin(); eit != edges.end();
           eit++) {
        Range adj_tet;
        rval = moab.get_adjacencies(&*eit, 1, 3, false, adj_tet);
        CHKERRQ_MOAB(rval);
        adj_tet = intersect(adj_tet, tets);
        rval = moab.get_connectivity(*eit, conn, num_nodes, true);
        CHKERRQ_MOAB(rval);
        if (tH) {
          rval = moab.tag_get_data(tH, conn, num_nodes, coords);
          CHKERRQ_MOAB(rval);
        } else {
          rval = moab.get_coords(conn, num_nodes, coords);
          CHKERRQ_MOAB(rval);
        }
        double q = 1;
        ierr = mField.getInterface<Tools>()->minTetsQuality(adj_tet, q);
        CHKERRQ(ierr);
        length_map[q * norm_2(s0 - s1)] = *eit;
      }
      MoFEMFunctionReturnHot(0);
    }
  };

  /**
   * \brief Topological relations
   */
  struct Toplogy {

    MoFEM::CoreInterface &mField;
    Toplogy(MoFEM::CoreInterface &m_field) : mField(m_field) {}

    enum TYPE {
      FREE = 0,
      SKIN = 1 << 0,
      SURFACE = 1 << 1,
      SURFACE_SKIN = 1 << 2,
      FRONT_ENDS = 1 << 3,
      FIX_EDGES = 1 << 4,
      FIX_CORNERS = 1 << 5
    };

    typedef map<int, Range> SetsMap;

    PetscErrorCode classifyVerts(const Range &surface, const Range &tets,
                                 const Range &fixed_edges,
                                 const Range &corner_nodes,
                                 SetsMap &sets_map) const {
      moab::Interface &moab(mField.get_moab());
      Skinner skin(&moab);
      MoFEMFunctionBeginHot;

      sets_map[FIX_CORNERS].merge(corner_nodes);
      Range fixed_verts;
      rval = moab.get_connectivity(fixed_edges, fixed_verts, true);
      CHKERRQ_MOAB(rval);
      sets_map[FIX_EDGES].swap(fixed_verts);

      Range tets_skin;
      rval = skin.find_skin(0, tets, false, tets_skin);
      CHKERRQ_MOAB(rval);
      Range tets_skin_edges;
      rval = moab.get_adjacencies(tets_skin, 1, false, tets_skin_edges,
                                  moab::Interface::UNION);
      CHKERRQ_MOAB(rval);

      // surface skin
      Range surface_skin;
      rval = skin.find_skin(0, surface, false, surface_skin);
      CHKERRQ_MOAB(rval);
      Range front_in_the_body;
      front_in_the_body = subtract(surface_skin, tets_skin_edges);
      Range front_ends;
      rval = skin.find_skin(0, front_in_the_body, false, front_ends);
      CHKERRQ_MOAB(rval);
      sets_map[FRONT_ENDS].swap(front_ends);

      Range surface_skin_verts;
      rval = moab.get_connectivity(surface_skin, surface_skin_verts, true);
      CHKERRQ_MOAB(rval);
      sets_map[SURFACE_SKIN].swap(surface_skin_verts);

      // surface
      Range surface_verts;
      rval = moab.get_connectivity(surface, surface_verts, true);
      CHKERRQ_MOAB(rval);
      sets_map[SURFACE].swap(surface_verts);

      // skin
      Range tets_skin_verts;
      rval = moab.get_connectivity(tets_skin, tets_skin_verts, true);
      CHKERRQ_MOAB(rval);
      sets_map[SKIN].swap(tets_skin_verts);

      Range tets_verts;
      rval = moab.get_connectivity(tets, tets_verts, true);
      CHKERRQ_MOAB(rval);
      sets_map[FREE].swap(tets_verts);

      MoFEMFunctionReturnHot(0);
    }

    PetscErrorCode getProcTets(const Range &tets, const Range &edges_to_merge,
                               Range &proc_tets) const {
      moab::Interface &moab(mField.get_moab());
      MoFEMFunctionBeginHot;
      Range edges_to_merge_verts;
      rval = moab.get_connectivity(edges_to_merge, edges_to_merge_verts, true);
      CHKERRQ_MOAB(rval);
      Range edges_to_merge_verts_tets;
      rval = moab.get_adjacencies(edges_to_merge_verts, 3, false,
                                  edges_to_merge_verts_tets,
                                  moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      edges_to_merge_verts_tets = intersect(edges_to_merge_verts_tets, tets);
      proc_tets.swap(edges_to_merge_verts_tets);
      MoFEMFunctionReturnHot(0);
    }

    PetscErrorCode edgesToMerge(const Range &surface, const Range &tets,
                                Range &edges_to_merge) const {
      moab::Interface &moab(mField.get_moab());
      MoFEMFunctionBeginHot;

      Range surface_verts;
      rval = moab.get_connectivity(surface, surface_verts, true);
      CHKERRQ_MOAB(rval);
      Range surface_verts_edges;
      rval = moab.get_adjacencies(surface_verts, 1, false, surface_verts_edges,
                                  moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      edges_to_merge.merge(surface_verts_edges);
      Range tets_edges;
      rval = moab.get_adjacencies(tets, 1, false, tets_edges,
                                  moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      edges_to_merge = intersect(edges_to_merge, tets_edges);
      MoFEMFunctionReturnHot(0);
    }

    PetscErrorCode removeBadEdges(const Range &surface, const Range &tets,
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
      CHKERRQ_MOAB(rval);
      Range surface_skin;
      rval = skin.find_skin(0, surface, false, surface_skin);
      CHKERRQ_MOAB(rval);

      // end nodes
      Range tets_skin_edges;
      rval = moab.get_adjacencies(tets_skin, 1, false, tets_skin_edges,
                                  moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      Range surface_front;
      surface_front = subtract(surface_skin, tets_skin_edges);
      Range ends_nodes;
      rval = skin.find_skin(0, surface_front, false, ends_nodes);
      CHKERRQ_MOAB(rval);

      // remove bad merges

      // get surface and body skin verts
      Range surface_edges;
      rval = moab.get_adjacencies(surface, 1, false, surface_edges,
                                  moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      Range surface_edges_verts;
      rval = moab.get_connectivity(surface_edges, surface_edges_verts, true);
      CHKERRQ_MOAB(rval);
      Range tets_skin_edges_verts;
      rval =
          moab.get_connectivity(tets_skin_edges, tets_skin_edges_verts, true);
      CHKERRQ_MOAB(rval);

      Range edges_to_remove;

      // romove edges self connected to body skin
      {
        Range ents_nodes_and_edges;
        ents_nodes_and_edges.merge(tets_skin_edges_verts);
        ents_nodes_and_edges.merge(tets_skin_edges);
        ierr = removeSelfConectingEdges(ents_nodes_and_edges, edges_to_remove,
                                        false);
        CHKERRQ(ierr);
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
                                        false);
        CHKERRQ(ierr);
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
                                        false);
        CHKERRQ(ierr);
      }
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      // remove edges self conected to surface
      ierr = removeSelfConectingEdges(surface_edges, edges_to_remove, false);
      CHKERRQ(ierr);
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      // remove edges self contected on surface skin
      {
        Range ents_nodes_and_edges;
        ents_nodes_and_edges.merge(surface_skin);
        ents_nodes_and_edges.merge(fixed_edges_nodes);
        ents_nodes_and_edges.merge(edges_to_remove);
        ierr = removeSelfConectingEdges(ents_nodes_and_edges, edges_to_remove,
                                        false);
        CHKERRQ(ierr);
      }
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      MoFEMFunctionReturnHot(0);
    }

  private:
    PetscErrorCode removeSelfConectingEdges(const Range &ents,
                                            Range &edges_to_remove,
                                            bool debug) const {
      moab::Interface &moab(mField.get_moab());
      MoFEMFunctionBeginHot;
      // get nodes
      Range ents_nodes = ents.subset_by_type(MBVERTEX);
      if (ents_nodes.empty()) {
        rval = moab.get_connectivity(ents, ents_nodes, true);
        CHKERRQ_MOAB(rval);
      }
      // edges adj. to nodes
      Range ents_nodes_edges;
      rval = moab.get_adjacencies(ents_nodes, 1, false, ents_nodes_edges,
                                  moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      // nodes of adj. edges
      Range ents_nodes_edges_nodes;
      rval =
          moab.get_connectivity(ents_nodes_edges, ents_nodes_edges_nodes, true);
      CHKERRQ_MOAB(rval);
      // hanging nodes
      ents_nodes_edges_nodes = subtract(ents_nodes_edges_nodes, ents_nodes);
      Range ents_nodes_edges_nodes_edges;
      rval = moab.get_adjacencies(ents_nodes_edges_nodes, 1, false,
                                  ents_nodes_edges_nodes_edges,
                                  moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      // remove eges adj. to hanging edges
      ents_nodes_edges =
          subtract(ents_nodes_edges, ents_nodes_edges_nodes_edges);
      ents_nodes_edges =
          subtract(ents_nodes_edges, ents.subset_by_type(MBEDGE));
      edges_to_remove.swap(ents_nodes_edges);
      if (debug) {
        EntityHandle meshset;
        rval = moab.create_meshset(MESHSET_SET, meshset);
        CHKERRQ_MOAB(rval);
        rval = moab.add_entities(meshset, edges_to_remove);
        CHKERRQ_MOAB(rval);
        rval = moab.write_file("edges_to_remove.vtk", "VTK", "", &meshset, 1);
        CHKERRQ_MOAB(rval);
        rval = moab.delete_entities(&meshset, 1);
        CHKERRQ_MOAB(rval);
      }
      MoFEMFunctionReturnHot(0);
    }
  };

  Range not_merged_edges;
  ierr = Toplogy(m_field).edgesToMerge(surface, tets, edges_to_merge);
  CHKERRQ(ierr);
  ierr =
      Toplogy(m_field).removeBadEdges(surface, tets, fixed_edges, corner_nodes,
                                      edges_to_merge, not_merged_edges);
  CHKERRQ(ierr);
  Toplogy::SetsMap sets_map;
  ierr = Toplogy(m_field).classifyVerts(surface, tets, fixed_edges,
                                        corner_nodes, sets_map);
  CHKERRQ(ierr);
  Range proc_tets;
  ierr = Toplogy(m_field).getProcTets(tets, edges_to_merge, proc_tets);
  CHKERRQ(ierr);
  out_tets = subtract(tets, proc_tets);
  if (bit_ptr) {
    for (int dd = 2; dd >= 0; dd--) {
      rval = moab.get_adjacencies(out_tets.subset_by_dimension(3), dd, false,
                                  out_tets, moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
    }
    ierr = m_field.getInterface<BitRefManager>()->addBitRefLevel(out_tets,
                                                                 *bit_ptr);
    CHKERRQ(ierr);
  }

  int nb_nodes_merged = 0;
  std::map<double, EntityHandle> length_map;
  new_surf = surface;

  // Range all_ents,sum_of_all_ents;
  // rval = moab.get_entities_by_handle(0, all_ents, true);
  // CHKERRQ_MOAB(rval);

  for (int pp = 0; pp != nbMaxMergingCycles; pp++) {

    int nb_nodes_merged_0 = nb_nodes_merged;
    length_map.clear();
    ierr = LengthMap(m_field, th)(proc_tets, edges_to_merge, length_map);
    CHKERRQ(ierr);

    int nn = 0;
    Range collapsed_edges;
    for (std::map<double, EntityHandle>::iterator mit = length_map.begin();
         mit != length_map.end(); mit++, nn++) {
      int num_nodes;
      const EntityHandle *conn;
      rval = moab.get_connectivity(mit->second, conn, num_nodes, true);
      CHKERRQ_MOAB(rval);
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
      
      ierr = MergeNodes(m_field, true, line_search, th,
                        update_meshsets)(father, mother, proc_tets, new_surf,
                                         edges_to_merge, not_merged_edges);
      CHKERRQ(ierr);
      // sum_of_all_ents.merge(proc_tets);

      if (m_field.getInterface<NodeMergerInterface>()->getSucessMerge()) {
        Range adj_edges;
        rval = moab.get_adjacencies(conn, 2, 1, false, adj_edges,
                                    moab::Interface::UNION);
        CHKERRQ_MOAB(rval);
        std::map<double, EntityHandle>::iterator miit = mit;
        if ((++miit) != length_map.end()) {
          std::map<double, EntityHandle> tmp_length_map = length_map;
          miit = tmp_length_map.find(miit->first);
          for (; miit != tmp_length_map.end(); miit++) {
            if (adj_edges.find(miit->second) != adj_edges.end()) {
              length_map.erase(miit->first);
            }
          }
        }
        // cerr << length_map.size() << " " << adj_edges << endl;
        nb_nodes_merged++;
        collapsed_edges.insert(mit->second);
      }
      if (nn > length_map.size() / fraction_level)
        break;
    }

    Range adj_faces, adj_edges;
    rval = moab.get_adjacencies(proc_tets, 2, false, adj_faces,
                                moab::Interface::UNION);
    CHKERRQ_MOAB(rval);
    new_surf = intersect(new_surf, adj_faces);

    rval = moab.get_adjacencies(proc_tets, 1, false, adj_edges,
                                moab::Interface::UNION);
    CHKERRQ_MOAB(rval);
    edges_to_merge = intersect(edges_to_merge, adj_edges);
    ierr = Toplogy(m_field).removeBadEdges(new_surf, proc_tets, fixed_edges,
                                           corner_nodes, edges_to_merge,
                                           not_merged_edges);
    CHKERRQ(ierr);

    PetscPrintf(m_field.get_comm(), "(%d) Number of nodes merged %d\n", pp,
                nb_nodes_merged);

    if(debug) {
      EntityHandle meshset_new_faces;
      rval = moab.create_meshset(MESHSET_SET, meshset_new_faces);
      CHKERRQ_MOAB(rval);
      rval = moab.add_entities(meshset_new_faces, new_surf);
      CHKERRQ_MOAB(rval);
      rval = moab.add_entities(meshset_new_faces, collapsed_edges);
      CHKERRQ_MOAB(rval);
      std::string name = "node_merge_step_" + boost::lexical_cast<std::string>(pp) + ".vtk";
      rval = moab.write_file(name.c_str(), "VTK", "", &meshset_new_faces, 1);
      CHKERRQ_MOAB(rval);
      rval = moab.delete_entities(&meshset_new_faces, 1);
      CHKERRQ_MOAB(rval);
    }

    if (nb_nodes_merged == nb_nodes_merged_0)
      break;
  }

  if (bit_ptr) {
    ierr = m_field.getInterface<BitRefManager>()->setBitRefLevel(proc_tets,
                                                                 *bit_ptr);
    CHKERRQ(ierr);
  }
  out_tets.merge(proc_tets);

  // // delete left created and no needed ents
  // {
  //   Range adj;
  //   for (int dd = 2; dd >= 0; dd--) {
  //     rval = moab.get_adjacencies(sum_of_all_ents, dd,
  //                                 false, adj, moab::Interface::UNION);
  //     CHKERRQ_MOAB(rval);
  //   }
  //   sum_of_all_ents.merge(adj);
  // }

  // sum_of_all_ents = subtract(sum_of_all_ents,all_ents);
  // sum_of_all_ents = subtract(sum_of_all_ents,out_tets);
  // {
  //   Range adj;
  //   for (int dd = 2; dd >= 0; dd--) {
  //     rval = moab.get_adjacencies(out_tets, dd, false, adj,
  //                                 moab::Interface::UNION);
  //     CHKERRQ_MOAB(rval);
  //   }
  //   sum_of_all_ents = subtract(sum_of_all_ents, adj);
  // }
  // rval = moab.delete_entities(sum_of_all_ents); CHKERRQ_MOAB(rval);

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::mergeBadEdges(
    const int fraction_level, const BitRefLevel trim_bit,
    const BitRefLevel cut_bit, const BitRefLevel bit, const Range &surface,
    const Range &fixed_edges, const Range &corner_nodes, Tag th,
    const bool debug) {
  MoFEM::CoreInterface &m_field = cOre;
  MoFEMFunctionBeginHot;
  Range tets_level;
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      trim_bit, BitRefLevel().set(), MBTET, tets_level);
  CHKERRQ(ierr);

  Range edges_to_merge;
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByParentType(
      trim_bit, trim_bit | cut_bit, MBEDGE, edges_to_merge);
  CHKERRQ(ierr);
  edges_to_merge = edges_to_merge.subset_by_type(MBEDGE);

  // get all entities not in databse
  Range all_ents_not_in_databse_before;
  ierr = cOre.getInterface<BitRefManager>()->getAllEntitiesNotInDatabase(
      all_ents_not_in_databse_before);
  CHKERRQ(ierr);

  Range out_new_tets, out_new_surf;
  ierr = mergeBadEdges(fraction_level, tets_level, surface, fixed_edges,
                       corner_nodes, edges_to_merge, out_new_tets, out_new_surf,
                       th, true, &bit, debug);
  CHKERRQ(ierr);

  // get all entities not in database after merge
  Range all_ents_not_in_databse_after;
  ierr = cOre.getInterface<BitRefManager>()->getAllEntitiesNotInDatabase(
      all_ents_not_in_databse_after);
  CHKERRQ(ierr);
  // delete hanging entities
  all_ents_not_in_databse_after =
      subtract(all_ents_not_in_databse_after, all_ents_not_in_databse_before);
  m_field.get_moab().delete_entities(all_ents_not_in_databse_after);

  mergedVolumes.swap(out_new_tets);
  mergedSurfaces.swap(out_new_surf);
  MoFEMFunctionReturnHot(0);
}

#ifdef WITH_TETGEN

PetscErrorCode CutMeshInterface::rebuildMeshWithTetGen(
    vector<string> &switches, const BitRefLevel &mesh_bit,
    const BitRefLevel &bit, const Range &surface, const Range &fixed_edges,
    const Range &corner_nodes, Tag th, const bool debug) {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  TetGenInterface *tetgen_iface;
  MoFEMFunctionBeginHot;
  ierr = m_field.getInterface(tetgen_iface);
  CHKERRQ(ierr);

  tetGenData.clear();
  moabTetGenMap.clear();
  tetGenMoabMap.clear();

  if (tetGenData.size() < 1) {
    tetGenData.push_back(new tetgenio);
  }
  tetgenio &in = tetGenData.back();

  struct BitEnts {

    MoFEM::CoreInterface &mField;
    const BitRefLevel &bIt;
    BitEnts(MoFEM::CoreInterface &m_field, const BitRefLevel &bit)
        : mField(m_field), bIt(bit) {}

    Range mTets;
    Range mTris;
    Range mEdges;
    Range mNodes;

    PetscErrorCode getLevelEnts() {
      MoFEMFunctionBeginHot;
      ierr = mField.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
          bIt, BitRefLevel().set(), MBTET, mTets);
      CHKERRQ(ierr);
      ierr = mField.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
          bIt, BitRefLevel().set(), MBTRI, mTris);
      CHKERRQ(ierr);
      ierr = mField.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
          bIt, BitRefLevel().set(), MBEDGE, mEdges);
      CHKERRQ(ierr);
      ierr = mField.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
          bIt, BitRefLevel().set(), MBVERTEX, mNodes);
      CHKERRQ(ierr);
      MoFEMFunctionReturnHot(0);
    }

    Range mSkin;
    Range mSkinNodes;
    Range mSkinEdges;

    PetscErrorCode getSkin() {
      moab::Interface &moab = mField.get_moab();
      MoFEMFunctionBeginHot;
      Skinner skin(&moab);
      rval = skin.find_skin(0, mTets, false, mSkin);
      CHKERRQ_MOAB(rval);
      rval = mField.get_moab().get_connectivity(mSkin, mSkinNodes, true);
      CHKERRQ_MOAB(rval);
      rval = mField.get_moab().get_adjacencies(mSkin, 1, false, mSkinEdges,
                                               moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      MoFEMFunctionReturnHot(0);
    }
  };

  struct SurfaceEnts {

    MoFEM::CoreInterface &mField;
    SurfaceEnts(MoFEM::CoreInterface &m_field) : mField(m_field) {}

    Range sNodes;
    Range sEdges;
    Range sVols;
    Range vNodes;

    PetscErrorCode getVolume(const BitEnts &bit_ents, const Range &tris) {
      moab::Interface &moab = mField.get_moab();
      MoFEMFunctionBeginHot;
      rval = moab.get_connectivity(tris, sNodes, true);
      CHKERRQ_MOAB(rval);
      rval =
          moab.get_adjacencies(tris, 1, false, sEdges, moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      rval =
          moab.get_adjacencies(sNodes, 3, false, sVols, moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      sVols = intersect(sVols, bit_ents.mTets);
      rval = moab.get_connectivity(sVols, vNodes, true);
      CHKERRQ_MOAB(rval);
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

    PetscErrorCode getSkin(const BitEnts &bit_ents, const Range &tris,
                           const int levels = 3) {
      moab::Interface &moab = mField.get_moab();
      MoFEMFunctionBeginHot;
      Skinner skin(&moab);
      rval = skin.find_skin(0, sVols, false, vSkin);
      CHKERRQ_MOAB(rval);
      for (int ll = 0; ll != levels; ll++) {
        rval = moab.get_adjacencies(vSkin, 3, false, sVols,
                                    moab::Interface::UNION);
        CHKERRQ_MOAB(rval);
        sVols = intersect(sVols, bit_ents.mTets);
        vSkin.clear();
        rval = skin.find_skin(0, sVols, false, vSkin);
        CHKERRQ_MOAB(rval);
      }
      vSkinWithoutBodySkin = subtract(vSkin, bit_ents.mSkin);
      vSkinOnBodySkin = intersect(vSkin, bit_ents.mSkin);
      rval = moab.get_connectivity(vSkinOnBodySkin, vSkinOnBodySkinNodes, true);
      CHKERRQ_MOAB(rval);
      rval = moab.get_connectivity(sVols, vNodes, true);
      CHKERRQ_MOAB(rval);
      rval = moab.get_connectivity(vSkin, vSkinNodes, true);
      CHKERRQ_MOAB(rval);
      vSkinNodesWithoutBodySkin = subtract(vSkinNodes, bit_ents.mSkinNodes);
      rval = skin.find_skin(0, tris, false, sSkin);
      CHKERRQ_MOAB(rval);
      rval = moab.get_connectivity(sSkin, sSkinNodes, true);
      CHKERRQ_MOAB(rval);
      tVols = sVols;
      MoFEMFunctionReturnHot(0);
    }

    Range tVols;

    PetscErrorCode getTetsForRemesh(const BitEnts &bit_ents, Tag th = NULL) {
      moab::Interface &moab = mField.get_moab();
      MoFEMFunctionBeginHot;

      Range tets_with_four_nodes_on_skin;
      rval = moab.get_adjacencies(vSkinOnBodySkinNodes, 3, false,
                                  tets_with_four_nodes_on_skin,
                                  moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      Range tets_nodes;
      rval =
          moab.get_connectivity(tets_with_four_nodes_on_skin, tets_nodes, true);
      CHKERRQ_MOAB(rval);
      tets_nodes = subtract(tets_nodes, vSkinOnBodySkinNodes);
      Range other_tets;
      rval = moab.get_adjacencies(tets_nodes, 3, false, other_tets,
                                  moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
      tets_with_four_nodes_on_skin =
          subtract(tets_with_four_nodes_on_skin, other_tets);
      Range to_remove;
      for (Range::iterator tit = tets_with_four_nodes_on_skin.begin();
           tit != tets_with_four_nodes_on_skin.end(); tit++) {
        int num_nodes;
        const EntityHandle *conn;
        rval = moab.get_connectivity(*tit, conn, num_nodes, true);
        CHKERRQ_MOAB(rval);
        double coords[12];
        if (th) {
          rval = moab.tag_get_data(th, conn, num_nodes, coords);
          CHKERRQ_MOAB(rval);
        } else {
          rval = moab.get_coords(conn, num_nodes, coords);
          CHKERRQ_MOAB(rval);
        }
        double quality = Tools::volumeLengthQuality(coords);
        if (quality < 1e-2) {
          to_remove.insert(*tit);
        }
      }

      sVols = subtract(sVols, to_remove);

      Skinner skin(&moab);
      vSkin.clear();
      rval = skin.find_skin(0, sVols, false, vSkin);
      CHKERRQ_MOAB(rval);
      Range m_skin;
      rval =
          skin.find_skin(0, subtract(bit_ents.mSkin, to_remove), false, m_skin);
      CHKERRQ_MOAB(rval);

      vSkinWithoutBodySkin = subtract(vSkin, m_skin);
      vSkinOnBodySkin = intersect(vSkin, m_skin);
      vNodes.clear();
      vSkinNodes.clear();
      vSkinOnBodySkinNodes.clear();
      rval = moab.get_connectivity(sVols, vNodes, true);
      CHKERRQ_MOAB(rval);
      rval = moab.get_connectivity(vSkinOnBodySkin, vSkinOnBodySkinNodes, true);
      CHKERRQ_MOAB(rval);
      rval = moab.get_connectivity(vSkin, vSkinNodes, true);
      CHKERRQ_MOAB(rval);

      MoFEMFunctionReturnHot(0);
    }
  };

  BitEnts bit_ents(m_field, mesh_bit);
  ierr = bit_ents.getLevelEnts();
  CHKERRQ(ierr);
  ierr = bit_ents.getSkin();
  CHKERRQ(ierr);
  SurfaceEnts surf_ents(m_field);
  ierr = surf_ents.getVolume(bit_ents, surface);
  CHKERRQ(ierr);
  ierr = surf_ents.getSkin(bit_ents, surface);
  CHKERRQ(ierr);
  ierr = surf_ents.getTetsForRemesh(bit_ents);
  CHKERRQ(ierr);

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
  int def_marker = 0;
  rval = m_field.get_moab().tag_get_handle(
      "TETGEN_MARKER", 1, MB_TYPE_INTEGER, th_marker,
      MB_TAG_CREAT | MB_TAG_SPARSE, &def_marker);
  CHKERRQ_MOAB(rval);

  vector<int> markers(surf_ents.vNodes.size(), 0);
  rval = moab.tag_set_data(th_marker, surf_ents.vNodes, &*markers.begin());
  CHKERRQ_MOAB(rval);
  {
    rval =
        m_field.get_moab().tag_get_data(th_marker, surface, &*markers.begin());
    CHKERRQ_MOAB(rval);
    fill(markers.begin(), markers.end(), 1);
    rval =
        m_field.get_moab().tag_set_data(th_marker, surface, &*markers.begin());
    CHKERRQ_MOAB(rval);
  }
  int shift = 3;
  map<int, int> id_shift_map; // each meshset has set unique bit
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, SIDESET, it)) {
    int ms_id = it->getMeshsetId();
    id_shift_map[ms_id] = 1 << shift; // shift bit
    shift++;
    Range sideset_faces;
    ierr = m_field.getInterface<MeshsetsManager>()->getEntitiesByDimension(
        ms_id, SIDESET, 2, sideset_faces, true);
    CHKERRQ(ierr);
    sideset_faces = intersect(sideset_faces, surf_ents.vNodes);
    markers.resize(sideset_faces.size());
    rval = m_field.get_moab().tag_get_data(th_marker, sideset_faces,
                                           &*markers.begin());
    CHKERRQ_MOAB(rval);
    for (unsigned int ii = 0; ii < markers.size(); ii++) {
      markers[ii] |= id_shift_map[ms_id]; // add bit to marker
    }
    rval = m_field.get_moab().tag_set_data(th_marker, sideset_faces,
                                           &*markers.begin());
    CHKERRQ_MOAB(rval);
  }
  Range nodes_to_remove; // none
  markers.resize(nodes_to_remove.size());
  fill(markers.begin(), markers.end(), -1);
  rval = m_field.get_moab().tag_set_data(th_marker, nodes_to_remove,
                                         &*markers.begin());
  CHKERRQ_MOAB(rval);

  // nodes
  if (tetGenData.size() == 1) {

    Range ents_to_tetgen = surf_ents.sVols;
    for (int dd = 2; dd >= 0; dd--) {
      rval = m_field.get_moab().get_adjacencies(
          surf_ents.sVols, dd, false, ents_to_tetgen, moab::Interface::UNION);
      CHKERRQ_MOAB(rval);
    }

    // Load mesh to TetGen data structures
    ierr = tetgen_iface->inData(ents_to_tetgen, in, moabTetGenMap,
                                tetGenMoabMap, th);
    CHKERRQ(ierr);
    ierr =
        tetgen_iface->setGeomData(in, moabTetGenMap, tetGenMoabMap, types_ents);
    CHKERRQ(ierr);

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
    ierr = tetgen_iface->setFaceData(markers, in, moabTetGenMap, tetGenMoabMap);
    CHKERRQ(ierr);
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
      ierr = tetgen_iface->tetRahedralize(s, _in_, _out_);
      CHKERRQ(ierr);
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

  ierr = tetgen_iface->outData(in, out, moabTetGenMap, tetGenMoabMap, bit,
                               false, false);
  CHKERRQ(ierr);

  Range rest_of_ents = subtract(bit_ents.mTets, surf_ents.tVols);
  for (int dd = 2; dd >= 0; dd--) {
    rval = moab.get_adjacencies(rest_of_ents.subset_by_dimension(3), dd, false,
                                rest_of_ents, moab::Interface::UNION);
    CHKERRQ_MOAB(rval);
  }
  ierr =
      m_field.getInterface<BitRefManager>()->addBitRefLevel(rest_of_ents, bit);
  CHKERRQ(ierr);

  Range tetgen_faces;
  map<int, Range> face_markers_map;
  ierr = tetgen_iface->getTriangleMarkers(tetGenMoabMap, out, &tetgen_faces,
                                          &face_markers_map);
  CHKERRQ(ierr);

  tetgenSurfaces = face_markers_map[1];
  for (map<int, Range>::iterator mit = face_markers_map.begin();
       mit != face_markers_map.end(); mit++) {
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, SIDESET, it)) {
      int msId = it->getMeshsetId();
      if (id_shift_map[msId] & mit->first) {
        EntityHandle meshset = it->getMeshset();
        ierr = m_field.get_moab().add_entities(
            meshset, mit->second.subset_by_type(MBTRI));
        CHKERRQ(ierr);
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

#endif // WITH_TETGEN

PetscErrorCode CutMeshInterface::setTagData(Tag th) {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  Range nodes;
  rval = moab.get_entities_by_type(0, MBVERTEX, nodes);
  CHKERRQ_MOAB(rval);
  std::vector<double> coords(3 * nodes.size());
  rval = moab.get_coords(nodes, &coords[0]);
  CHKERRQ_MOAB(rval);
  rval = moab.tag_set_data(th, nodes, &coords[0]);
  CHKERRQ_MOAB(rval);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::setCoords(Tag th) {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  Range nodes;
  rval = moab.get_entities_by_type(0, MBVERTEX, nodes);
  CHKERRQ_MOAB(rval);
  std::vector<double> coords(3 * nodes.size());
  rval = moab.tag_get_data(th, nodes, &coords[0]);
  CHKERRQ_MOAB(rval);
  rval = moab.set_coords(nodes, &coords[0]);
  CHKERRQ_MOAB(rval);
  MoFEMFunctionReturnHot(0);
}

struct SaveData {
  moab::Interface &moab;
  SaveData(moab::Interface &moab) : moab(moab) {}
  PetscErrorCode operator()(const std::string name, const Range &ents) {
    MoFEMFunctionBeginHot;
    EntityHandle meshset;
    rval = moab.create_meshset(MESHSET_SET, meshset);
    CHKERRQ_MOAB(rval);
    rval = moab.add_entities(meshset, ents);
    CHKERRQ_MOAB(rval);
    rval = moab.write_file(name.c_str(), "VTK", "", &meshset, 1);
    CHKERRQ_MOAB(rval);
    rval = moab.delete_entities(&meshset, 1);
    CHKERRQ_MOAB(rval);
    MoFEMFunctionReturnHot(0);
  }
};

PetscErrorCode CutMeshInterface::saveCutEdges() {
  MoFEM::CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;

  ierr = SaveData(moab)("out_vol.vtk", vOlume);
  CHKERRQ(ierr);
  ierr = SaveData(moab)("out_surface.vtk", sUrface);
  CHKERRQ(ierr);
  ierr = SaveData(moab)("out_cut_edges.vtk", cutEdges);
  CHKERRQ(ierr);
  ierr = SaveData(moab)("out_cut_volumes.vtk", cutVolumes);
  CHKERRQ(ierr);
  ierr = SaveData(moab)("out_cut_new_volumes.vtk", cutNewVolumes);
  CHKERRQ(ierr);
  ierr = SaveData(moab)("out_cut_new_surfaces.vtk", cutNewSurfaces);
  CHKERRQ(ierr);

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode CutMeshInterface::saveTrimEdges() {
  moab::Interface &moab = cOre.getInterface<CoreInterface>()->get_moab();
  MoFEMFunctionBeginHot;

  ierr = SaveData(moab)("out_trim_new_volumes.vtk", trimNewVolumes);
  CHKERRQ(ierr);
  ierr = SaveData(moab)("out_trim_new_surfaces.vtk", trimNewSurfaces);
  CHKERRQ(ierr);
  ierr = SaveData(moab)("out_trim_edges.vtk", trimEdges);
  CHKERRQ(ierr);
 
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM
