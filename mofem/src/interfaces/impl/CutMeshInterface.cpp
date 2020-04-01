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
  projectEntitiesQualityTrashold = 0.5;
  trimFixedEdges = true;
}

MoFEMErrorCode CutMeshInterface::clearMap() {
  MoFEMFunctionBegin;
  treeSurfPtr.reset();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::setFront(const Range front) {
  MoFEMFunctionBeginHot;
  fRont = front;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CutMeshInterface::setSurface(const Range surface) {
  MoFEMFunctionBeginHot;
  sUrface = surface;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CutMeshInterface::copySurface(const Range surface, Tag th,
                                             double *shift, double *origin,
                                             double *transform,
                                             const std::string save_mesh) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  sUrface.clear();
  std::map<EntityHandle, EntityHandle> verts_map;
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
      if (verts_map.find(conn[nn]) != verts_map.end()) {
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
  if (!save_mesh.empty())
    CHKERR SaveData(m_field.get_moab())(save_mesh, sUrface);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::setVolume(const Range volume) {
  MoFEMFunctionBeginHot;
  vOlume = volume;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CutMeshInterface::setConstrainSurface(const Range surf) {
  MoFEMFunctionBeginHot;
  constrainSurface = surf;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CutMeshInterface::mergeSurface(const Range surface) {
  MoFEMFunctionBeginHot;
  sUrface.merge(surface);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CutMeshInterface::mergeVolumes(const Range volume) {
  MoFEMFunctionBeginHot;
  vOlume.merge(volume);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CutMeshInterface::snapSurfaceSkinToEdges(const Range fixed_edges,
                                                        const double rel_tol,
                                                        const double abs_tol,
                                                        Tag th,
                                                        const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  // Get cutting surface skin
  Skinner skin(&moab);
  Range surface_skin;
  CHKERR skin.find_skin(0, sUrface, false, surface_skin);

  CHKERR snapSurfaceToEdges(surface_skin, fixed_edges, rel_tol, abs_tol, th,
                            debug);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::snapSurfaceToEdges(const Range surface_edges,
                                                    const Range fixed_edges,
                                                    const double rel_tol,
                                                    const double abs_tol,
                                                    Tag th, const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  FTensor::Index<'i', 3> i;
  MoFEMFunctionBegin;

  map<EntityHandle, double> map_verts_length;

  for (auto f : surface_edges) {
    int num_nodes;
    const EntityHandle *conn_skin;
    CHKERR moab.get_connectivity(f, conn_skin, num_nodes, true);
    VectorDouble6 coords_skin(6);
    if (th)
      CHKERR moab.tag_get_data(th, conn_skin, num_nodes, &coords_skin[0]);
    else
      CHKERR moab.get_coords(conn_skin, num_nodes, &coords_skin[0]);
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_n0(
        &coords_skin[0], &coords_skin[1], &coords_skin[2]);
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_n1(
        &coords_skin[3], &coords_skin[4], &coords_skin[5]);
    FTensor::Tensor1<double, 3> t_skin_delta;
    t_skin_delta(i) = t_n1(i) - t_n0(i);
    const double skin_edge_length = sqrt(t_skin_delta(i) * t_skin_delta(i));
    for (int nn = 0; nn != 2; ++nn) {
      auto it = map_verts_length.find(conn_skin[nn]);
      if (it != map_verts_length.end())
        it->second = std::min(it->second, skin_edge_length);
      else
        map_verts_length[conn_skin[nn]] = skin_edge_length;
    }
  }

  for (auto m : map_verts_length) {

    FTensor::Tensor1<double, 3> t_n;
    if (th)
      CHKERR moab.tag_get_data(th, &m.first, 1, &t_n(0));
    else
      CHKERR moab.get_coords(&m.first, 1, &t_n(0));

    double min_dist = rel_tol * m.second;
    FTensor::Tensor1<double, 3> t_min_coords;
    CHKERR cOre.getInterface<Tools>()->findMinDistanceFromTheEdges(
        &t_n(0), 1, fixed_edges, &min_dist, &t_min_coords(0));

    if (min_dist < rel_tol * m.second || min_dist < abs_tol) {
      if (debug)
        cerr << "Snap " << min_dist << endl;
      if (th)
        CHKERR moab.tag_set_data(th, &m.first, 1, &t_min_coords(0));
      else
        CHKERR moab.set_coords(&m.first, 1, &t_min_coords(0));
    }
  }

  MoFEMFunctionReturn(0);
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

MoFEMErrorCode
CutMeshInterface::cutOnly(Range vol, const BitRefLevel cut_bit, Tag th,
                          const double tol_geometry, const double tol_cut_close,
                          Range *fixed_edges, Range *corner_nodes,
                          const bool update_meshsets, const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  // cut mesh
  CHKERR findEdgesToCut(vol, QUIET, debug);
  CHKERR projectZeroDistanceEnts(fixed_edges, corner_nodes, tol_geometry,
                                 tol_cut_close, QUIET, debug);
  CHKERR cutEdgesInMiddle(cut_bit, cutNewVolumes, cutNewSurfaces,
                          cutNewVertices, debug);

  CHKERR cOre.getInterface<BitRefManager>()->updateRange(constrainSurface,
                                                         constrainSurface);
  if (fixed_edges)
    CHKERR cOre.getInterface<BitRefManager>()->updateRange(*fixed_edges,
                                                           *fixed_edges);
  if (corner_nodes)
    CHKERR cOre.getInterface<BitRefManager>()->updateRange(*corner_nodes,
                                                           *corner_nodes);
  if (update_meshsets)
    CHKERR m_field.getInterface<MeshsetsManager>()
        ->updateAllMeshsetsByEntitiesChildren(cut_bit);
  CHKERR moveMidNodesOnCutEdges(th);

  if (debug) {
    CHKERR saveCutEdges();
    if (fixed_edges)
      CHKERR SaveData(moab)("out_cut_new_fixed_edges.vtk", *fixed_edges);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::trimOnly(const BitRefLevel trim_bit, Tag th,
                                          const double tol_geometry,
                                          const double tol_trim_close,
                                          Range *fixed_edges,
                                          Range *corner_nodes,
                                          const bool update_meshsets,
                                          const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  // trim mesh
  CHKERR findEdgesToTrim(fixed_edges, corner_nodes, th, tol_geometry,
                         tol_trim_close, debug);
  CHKERR trimEdgesInTheMiddle(trim_bit, debug);

  CHKERR cOre.getInterface<BitRefManager>()->updateRange(constrainSurface,
                                                         constrainSurface);
  if (fixed_edges)
    CHKERR cOre.getInterface<BitRefManager>()->updateRange(*fixed_edges,
                                                           *fixed_edges);

  if (corner_nodes)
    CHKERR cOre.getInterface<BitRefManager>()->updateRange(*corner_nodes,
                                                           *corner_nodes);

  if (update_meshsets)
    CHKERR m_field.getInterface<MeshsetsManager>()
        ->updateAllMeshsetsByEntitiesChildren(trim_bit);

  // move nodes
  CHKERR moveMidNodesOnTrimmedEdges(th);

  // remove faces
  CHKERR trimSurface(fixed_edges, corner_nodes, debug);

  if (debug) {
    CHKERR saveTrimEdges();
    Range bit2_edges;
    CHKERR cOre.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        trim_bit, BitRefLevel().set(), MBEDGE, bit2_edges);
    CHKERR SaveData(moab)("trim_fixed_edges.vtk",
                          intersect(*fixed_edges, bit2_edges));
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::cutAndTrim(
    int &first_bit, Tag th, const double tol_geometry, const double tol_cut_close,
    const double tol_trim_close, Range *fixed_edges, Range *corner_nodes,
    const bool update_meshsets, const bool debug) {
  CoreInterface &m_field = cOre;
  MoFEMFunctionBegin;

  std::vector<BitRefLevel> bit_levels;

  auto get_back_bit_levels = [&]() {
    bit_levels.push_back(BitRefLevel());
    bit_levels.back().set(first_bit + bit_levels.size() - 1);
    return bit_levels.back();
  };

  BitRefLevel cut_bit = get_back_bit_levels();

  CHKERR cutOnly(unite(cutSurfaceVolumes, cutFrontVolumes), cut_bit, th,
                 tol_geometry, tol_cut_close, fixed_edges, corner_nodes,
                 update_meshsets, debug);

  auto get_min_quality = [&m_field](const BitRefLevel bit, Tag th) {
    Range tets_level; // test at level
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit, BitRefLevel().set(), MBTET, tets_level);
    double min_q = 1;
    CHKERR m_field.getInterface<Tools>()->minTetsQuality(tets_level, min_q, th);
    return min_q;
  };

  PetscPrintf(PETSC_COMM_WORLD, "Min quality cut %6.4g\n",
              get_min_quality(cut_bit, th));

  Range starting_volume = cutNewVolumes;
  Range starting_volume_nodes;
  CHKERR m_field.get_moab().get_connectivity(starting_volume,
                                             starting_volume_nodes, true);
  std::vector<double> staring_volume_coords(3 * starting_volume_nodes.size());
  CHKERR m_field.get_moab().get_coords(starting_volume_nodes,
                                       &*staring_volume_coords.begin());

  if (th) {
    std::vector<double> staring_volume_th_coords(3 *
                                                 starting_volume_nodes.size());
    CHKERR m_field.get_moab().tag_get_data(th, starting_volume_nodes,
                                           &*staring_volume_th_coords.begin());
    CHKERR m_field.get_moab().set_coords(starting_volume_nodes,
                                         &*staring_volume_th_coords.begin());
  }

  if (th)
    CHKERR setTagData(th);

  BitRefLevel trim_bit = get_back_bit_levels();

  CHKERR trimOnly(trim_bit, th, tol_geometry, tol_trim_close, fixed_edges,
                  corner_nodes, update_meshsets, debug);

  PetscPrintf(PETSC_COMM_WORLD, "Min quality trim %3.2g\n",
              get_min_quality(trim_bit, th));

  first_bit += bit_levels.size() - 1;

  if (th)
    CHKERR m_field.get_moab().set_coords(starting_volume_nodes,
                                         &*staring_volume_coords.begin());

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::cutTrimAndMerge(
    int &first_bit, const int fraction_level, Tag th, const double tol_geometry,
    const double tol_cut_close, const double tol_trim_close, Range &fixed_edges,
    Range &corner_nodes, const bool update_meshsets, const bool debug) {
  CoreInterface &m_field = cOre;
  MoFEMFunctionBegin;

  std::vector<BitRefLevel> bit_levels;

  auto get_back_bit_levels = [&]() {
    bit_levels.push_back(BitRefLevel());
    bit_levels.back().set(first_bit + bit_levels.size() - 1);
    return bit_levels.back();
  };

  if (debug) {
    CHKERR cOre.getInterface<BitRefManager>()->writeEntitiesNotInDatabase(
        "ents_not_in_database.vtk", "VTK", "");
  }

  CHKERR cutAndTrim(first_bit, th, tol_geometry, tol_cut_close, tol_trim_close,
                    &fixed_edges, &corner_nodes, update_meshsets, debug);
  if (debug)
    CHKERR cOre.getInterface<BitRefManager>()->writeEntitiesNotInDatabase(
        "cut_trim_ents_not_in_database.vtk", "VTK", "");

  BitRefLevel bit_level1 = BitRefLevel().set(first_bit - 1);
  BitRefLevel bit_level2 = get_back_bit_levels();
  BitRefLevel bit_level3 = get_back_bit_levels();

  CHKERR mergeBadEdges(fraction_level, bit_level1, bit_level2, bit_level3,
                       getNewTrimSurfaces(), fixed_edges, corner_nodes, th,
                       update_meshsets, debug);

  auto get_min_quality = [&m_field, debug](const BitRefLevel bit, Tag th) {
    Range tets_level; // test at level
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit, BitRefLevel().set(), MBTET, tets_level);
    double min_q = 1;
    CHKERR m_field.getInterface<Tools>()->minTetsQuality(tets_level, min_q, th);
    if (min_q < 0 && debug) {
      CHKERR m_field.getInterface<Tools>()->writeTetsWithQuality(
          "negative_tets.vtk", "VTK", "", tets_level, th);
    }
    return min_q;
  };

  PetscPrintf(PETSC_COMM_WORLD, "Min quality node merge %6.4g\n",
              get_min_quality(bit_level3, th));

  CHKERR cOre.getInterface<BitRefManager>()->updateRange(constrainSurface,
                                                         constrainSurface);
  CHKERR cOre.getInterface<BitRefManager>()->updateRange(fixed_edges,
                                                         fixed_edges);
  CHKERR cOre.getInterface<BitRefManager>()->updateRange(corner_nodes,
                                                         corner_nodes);

  first_bit += bit_levels.size() - 1;

  if (debug) {
    CHKERR cOre.getInterface<BitRefManager>()->writeBitLevelByType(
        bit_level3, BitRefLevel().set(), MBTET, "out_tets_merged.vtk", "VTK",
        "");
    CHKERR cOre.getInterface<BitRefManager>()->writeEntitiesNotInDatabase(
        "cut_trim_merge_ents_not_in_database.vtk", "VTK", "");
    CHKERR SaveData(m_field.get_moab())("merged_surface.vtk", mergedSurfaces);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::makeFront(const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  Skinner skin(&moab);
  Range tets_skin;
  CHKERR skin.find_skin(0, vOlume, false, tets_skin);
  Range tets_skin_edges;
  CHKERR moab.get_adjacencies(tets_skin, 1, false, tets_skin_edges,
                              moab::Interface::UNION);
  Range surface_skin;
  CHKERR skin.find_skin(0, sUrface, false, surface_skin);
  fRont = subtract(surface_skin, tets_skin_edges);
  if (debug)
    CHKERR SaveData(moab)("front.vtk", fRont);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::createSurfaceLevelSets(int verb,
                                                        const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  auto tools_interface = m_field.getInterface<Tools>();

  auto create_tag = [&](const std::string name, const int dim) {
    Tag th;
    rval = moab.tag_get_handle(name.c_str(), th);
    if (rval == MB_SUCCESS)
      return th;
    std::vector<double> def_val(dim, 0);
    CHKERR moab.tag_get_handle(name.c_str(), dim, MB_TYPE_DOUBLE, th,
                               MB_TAG_CREAT | MB_TAG_SPARSE, &*def_val.begin());

    return th;
  };

  auto set_vol = [&](const Range &vol_verts, std::vector<double> &coords,
                     std::vector<double> &dist_surface_vec,
                     std::vector<double> &dist_surface_normal_vec,
                     std::vector<double> &dist_surface_digned_dist_vec) {
    MoFEMFunctionBegin;

    coords.resize(3 * vol_verts.size());
    dist_surface_vec.resize(3 * vol_verts.size());
    dist_surface_normal_vec.resize(3 * vol_verts.size());
    dist_surface_digned_dist_vec.resize(vol_verts.size());

    CHKERR moab.get_coords(vol_verts, &*coords.begin());
    std::srand(0);

    for (auto v : vol_verts) {

      const int index = vol_verts.index(v);
      auto point_in = getVectorAdaptor(&coords[3 * index], 3);
      VectorDouble3 point_out(3);
      EntityHandle facets_out;
      CHKERR treeSurfPtr->closest_to_location(&point_in[0], rootSetSurf,
                                              &point_out[0], facets_out);

      VectorDouble3 n(3);
      CHKERR tools_interface->getTriNormal(facets_out, &*n.begin());
      n /= norm_2(n);

      VectorDouble3 delta = point_out - point_in;
      const double dot = inner_prod(delta, n);

      if (std::abs(dot) < std::numeric_limits<double>::epsilon()) {
        if (std::rand() % 2 == 0)
          delta += n * std::numeric_limits<double>::epsilon();
        else
          delta -= n * std::numeric_limits<double>::epsilon();
      }

      auto dist_vec = getVectorAdaptor(&dist_surface_vec[3 * index], 3);
      noalias(dist_vec) = delta;

      auto dist_normal_vec =
          getVectorAdaptor(&dist_surface_normal_vec[3 * index], 3);

      dist_surface_digned_dist_vec[index] = dot;
      noalias(dist_normal_vec) = dot * n;
    }

    MoFEMFunctionReturn(0);
  };

  std::vector<double> coords;
  std::vector<double> dist_surface_vec;
  std::vector<double> dist_surface_normal_vec;
  std::vector<double> dist_surface_signed_dist_vec;
  Range vol_verts;
  CHKERR moab.get_connectivity(vOlume, vol_verts, true);

  CHKERR set_vol(vol_verts, coords, dist_surface_vec, dist_surface_normal_vec,
                 dist_surface_signed_dist_vec);

  CHKERR moab.tag_set_data(create_tag("DIST_SURFACE_VECTOR", 3), vol_verts,
                           &*dist_surface_vec.begin());
  CHKERR moab.tag_set_data(create_tag("DIST_SURFACE_NORMAL_VECTOR", 3),
                           vol_verts, &*dist_surface_normal_vec.begin());
  CHKERR moab.tag_set_data(create_tag("DIST_SURFACE_NORMAL_SIGNED", 1),
                           vol_verts, &*dist_surface_signed_dist_vec.begin());

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::createFrontLevelSets(Range vol, Tag th,
                                                      int verb,
                                                      const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  auto create_tag = [&](const std::string name, const int dim) {
    Tag th;
    rval = moab.tag_get_handle(name.c_str(), th);
    if (rval == MB_SUCCESS)
      return th;
    std::vector<double> def_val(dim, 0);
    CHKERR moab.tag_get_handle(name.c_str(), dim, MB_TYPE_DOUBLE, th,
                               MB_TAG_CREAT | MB_TAG_SPARSE, &*def_val.begin());
    return th;
  };

  Range vol_vertices;
  CHKERR moab.get_connectivity(vol, vol_vertices, true);

  std::vector<double> min_distances_from_front(vol_vertices.size(), -1);
  std::vector<double> points_on_edges(3 * vol_vertices.size(), 0);
  std::vector<EntityHandle> closest_edges(vol_vertices.size(), 0);

  std::vector<double> coords(3 * vol_vertices.size());
  if (th)
    CHKERR moab.tag_get_data(th, vol_vertices, &*coords.begin());
  else
    CHKERR moab.get_coords(vol_vertices, &*coords.begin());

  CHKERR cOre.getInterface<Tools>()->findMinDistanceFromTheEdges(
      &*coords.begin(), vol_vertices.size(), fRont,
      &*min_distances_from_front.begin(), &*points_on_edges.begin(),
      &*closest_edges.begin());

  if (!points_on_edges.empty()) {
    for (int i = 0; i != min_distances_from_front.size(); ++i) {
      Range faces;
      CHKERR moab.get_adjacencies(&closest_edges[0], 1, 2, false, faces);
      auto point_in = getVectorAdaptor(&coords[3 * i], 3);
      auto point_out = getVectorAdaptor(&points_on_edges[3 * i], 3);
      point_out -= point_in;
    }
  }

  auto th_dist_front_vec = create_tag("DIST_FRONT_VECTOR", 3);
  CHKERR moab.tag_set_data(th_dist_front_vec, vol_vertices,
                           &*points_on_edges.begin());

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::findLevelSetVolumes(
    Tag th, Range &vol_edges, const bool remove_adj_prims_edges, int verb,
    const bool debug, const std::string edges_file_name) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  auto get_tag_dim = [&](auto th) {
    int dim;
    CHKERR moab.tag_get_length(th, dim);
    return dim;
  };
  auto dim = get_tag_dim(th);

  auto get_tag_data = [&](auto th, auto conn) {
    const void *ptr;
    CHKERR moab.tag_get_by_ptr(th, &conn, 1, &ptr);
    return getVectorAdaptor(
        const_cast<double *>(static_cast<const double *>(ptr)), 3);
  };

  auto get_edge_ray = [&](auto conn) {
    VectorDouble6 coords(6);
    CHKERR moab.get_coords(conn, 2, &coords[0]);
    VectorAdaptor n0 = getVectorAdaptor(&coords[0], 3);
    VectorAdaptor n1 = getVectorAdaptor(&coords[3], 3);
    VectorDouble3 ray = n1 - n0;
    return ray;
  };

  Range edges;
  CHKERR moab.get_adjacencies(vOlume, 1, true, edges, moab::Interface::UNION);

  auto remove_prisms_edges = [&](Range &edges) {
    MoFEMFunctionBegin;
    Range prisms;
    CHKERR moab.get_adjacencies(edges, 3, false, prisms,
                                moab::Interface::UNION);
    prisms = prisms.subset_by_type(MBPRISM);
    Range prisms_verts;
    CHKERR moab.get_connectivity(prisms, prisms_verts, true);
    Range prism_edges;
    CHKERR moab.get_adjacencies(prisms_verts, 1, false, prism_edges,
                                moab::Interface::UNION);
    edges = subtract(edges, prism_edges);
    MoFEMFunctionReturn(0);
  };
  if (remove_adj_prims_edges)
    CHKERR remove_prisms_edges(edges);

  std::vector<EntityHandle> cut_edges;
  cut_edges.reserve(edges.size());

  auto get_cut_edges_vec = [&](auto th, auto &cut_edges, auto e, auto &&conn) {
    MoFEMFunctionBegin;
    auto ray = get_edge_ray(conn);
    const double length = norm_2(ray);
    ray /= length;
    auto signed_norm = [&](const auto &v) { return inner_prod(ray, v); };
    const auto dist0 = get_tag_data(th, conn[0]);
    const auto dist1 = get_tag_data(th, conn[1]);
    // const double min_dist = std::min(norm_2(dist0), norm_2(dist1));
    // if (min_dist < 2 * length) {
    auto opposite = inner_prod(dist0, dist1);
    if (opposite <= std::numeric_limits<double>::epsilon()) {
      const double sign_dist0 = signed_norm(dist0);
      const double sign_dist1 = signed_norm(dist1);
      if (sign_dist0 > -std::numeric_limits<double>::epsilon() &&
          sign_dist1 < std::numeric_limits<double>::epsilon())
        cut_edges.push_back(e);
    }
    // }
    MoFEMFunctionReturn(0);
  };

  auto get_cut_edges_signed_dist = [&](auto th, auto &cut_edges, auto e,
                                       auto &&conn) {
    MoFEMFunctionBegin;
    auto get_tag_signed_dist = [&](auto conn) {
      double dist;
      CHKERR moab.tag_get_data(th, &conn, 1, &dist);
      return dist;
    };
    const auto dist0 = get_tag_signed_dist(conn[0]);
    const auto dist1 = get_tag_signed_dist(conn[1]);
    const double opposite = dist0 * dist1;
    if (opposite < 0)
      cut_edges.push_back(e);
    MoFEMFunctionReturn(0);
  };

  auto get_conn = [&](auto e) {
    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(e, conn, num_nodes, true);
    return conn;
  };

  if (dim == 3)
    for (auto e : edges)
      CHKERR get_cut_edges_vec(th, cut_edges, e, get_conn(e));
  else
    for (auto e : edges)
      CHKERR get_cut_edges_signed_dist(th, cut_edges, e, get_conn(e));

  CHKERR moab.get_adjacencies(&*cut_edges.begin(), cut_edges.size(), 3, false,
                              vol_edges, moab::Interface::UNION);
  vol_edges = intersect(vol_edges, vOlume);

  auto convert_to_ranger = [](auto &v) {
    Range r;
    r.insert_list(v.begin(), v.end());
    return r;
  };

  if (debug && !edges_file_name.empty())
    CHKERR SaveData(m_field.get_moab())(edges_file_name,
                                        convert_to_ranger(cut_edges));

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::findLevelSetVolumes(int verb,
                                                     const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  CHKERR createFrontLevelSets(vOlume, nullptr, verb, debug);
  Tag th_dist_front_vec;
  CHKERR moab.tag_get_handle("DIST_FRONT_VECTOR", th_dist_front_vec);
  CHKERR findLevelSetVolumes(th_dist_front_vec, cutFrontVolumes, true, verb,
                             debug, "cutFrontEdges.vtk");

  CHKERR createSurfaceLevelSets(verb, debug);

  Tag th_dist_surface_vec;
  CHKERR moab.tag_get_handle("DIST_SURFACE_VECTOR", th_dist_surface_vec);
  cutSurfaceVolumes.clear();
  CHKERR findLevelSetVolumes(th_dist_surface_vec, cutSurfaceVolumes, true, verb,
                             debug, "cutSurfaceEdges.vtk");
  Tag th_dist_signed_normal;
  CHKERR moab.tag_get_handle("DIST_SURFACE_NORMAL_SIGNED",
                             th_dist_signed_normal);
  Range signed_vol;
  CHKERR findLevelSetVolumes(th_dist_signed_normal, signed_vol, true, verb,
                             debug, "cutSurfaceEdges.vtk");
  cutSurfaceVolumes = intersect(cutSurfaceVolumes, signed_vol);

  if (debug)
    CHKERR SaveData(m_field.get_moab())("level_sets.vtk", vOlume);
  if (debug)
    CHKERR SaveData(m_field.get_moab())("cutSurfaceVolumes.vtk",
                                        cutSurfaceVolumes);
  if (debug)
    CHKERR SaveData(m_field.get_moab())("cutFrontVolumes.vtk", cutFrontVolumes);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::refineMesh(const int init_bit_level,
                                            const int surf_levels,
                                            const int front_levels,
                                            Range *fixed_edges, int verb,
                                            const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MeshRefinement *refiner;
  BitRefManager *bit_ref_manager;
  MoFEMFunctionBegin;
  CHKERR m_field.getInterface(refiner);
  CHKERR m_field.getInterface(bit_ref_manager);

  auto add_bit = [&](const int bit) {
    MoFEMFunctionBegin;
    CHKERR bit_ref_manager->addBitRefLevel(vOlume, BitRefLevel().set(bit),
                                           verb);
    Range adj_ents;
    for (auto d : {2, 1, 0})
      CHKERR moab.get_adjacencies(vOlume, d, true, adj_ents,
                                  moab::Interface::UNION);
    CHKERR bit_ref_manager->addBitRefLevel(vOlume, BitRefLevel().set(bit),
                                           verb);
    MoFEMFunctionReturn(0);
  };
  CHKERR add_bit(init_bit_level);

  auto update_range = [&](Range *r_ptr) {
    MoFEMFunctionBegin;
    if (r_ptr) {
      Range childs;
      CHKERR bit_ref_manager->updateRange(*r_ptr, childs);
      r_ptr->merge(childs);
    }
    MoFEMFunctionReturn(0);
  };

  auto refine = [&](const BitRefLevel bit, const Range tets) {
    MoFEMFunctionBegin;
    Range verts;
    CHKERR moab.get_connectivity(tets, verts, true);
    Range ref_edges;
    CHKERR moab.get_adjacencies(verts, 1, true, ref_edges,
                                moab::Interface::UNION);

    CHKERR refiner->add_vertices_in_the_middle_of_edges(ref_edges, bit);
    CHKERR refiner->refine_TET(vOlume, bit, false, verb);

    CHKERR update_range(fixed_edges);
    CHKERR update_range(&vOlume);
    CHKERR m_field.getInterface<MeshsetsManager>()
        ->updateAllMeshsetsByEntitiesChildren(bit);

    Range bit_tets;
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit, BitRefLevel().set(), MBTET, bit_tets);
    vOlume = intersect(vOlume, bit_tets);

    Range bit_edges;
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit, BitRefLevel().set(), MBEDGE, bit_edges);
    if (fixed_edges)
      *fixed_edges = intersect(*fixed_edges, bit_edges);

    MoFEMFunctionReturn(0);
  };

  for (int ll = init_bit_level; ll != init_bit_level + surf_levels; ++ll) {
    CHKERR findLevelSetVolumes(verb, debug);
    CHKERR refine(BitRefLevel().set(ll + 1),
                  unite(cutSurfaceVolumes, cutFrontVolumes));
  }

  for (int ll = init_bit_level + surf_levels;
       ll != init_bit_level + surf_levels + front_levels; ++ll) {
    CHKERR findLevelSetVolumes(verb, debug);
    CHKERR refine(BitRefLevel().set(ll + 1), cutFrontVolumes);
  }

  CHKERR findLevelSetVolumes(verb, debug);

  if (debug)
    CHKERR SaveData(m_field.get_moab())("refinedVolume.vtk", vOlume);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::findEdgesToCut(Range vol, int verb,
                                                const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  edgesToCut.clear();
  cutEdges.clear();

  Tag th_signed_dist;
  CHKERR moab.tag_get_handle("DIST_SURFACE_NORMAL_SIGNED", th_signed_dist);

  auto get_tag_edge_dist = [&](auto th, auto conn) {
    std::array<double, 2> r;
    CHKERR moab.tag_get_data(th, conn, 2, r.data());
    return r;
  };

  auto get_conn = [&](const auto e) {
    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(e, conn, num_nodes, true);
    return conn;
  };

  auto get_adj = [&moab](const Range r, const int dim) {
    Range a;
    if (dim)
      CHKERR moab.get_adjacencies(r, dim, false, a, moab::Interface::UNION);
    else
      CHKERR moab.get_connectivity(r, a, true);
    return a;
  };

  auto get_range = [](std::vector<EntityHandle> v) {
    Range r;
    r.insert_list(v.begin(), v.end());
    return r;
  };

  auto vol_edges = get_adj(vol, 1);

  aveLength = 0;
  maxLength = 0;
  int nb_ave_length = 0;
  for (auto e : vol_edges) {

    auto conn = get_conn(e);
    auto dist = get_tag_edge_dist(th_signed_dist, conn);
    const auto dist_max = std::max(dist[0], dist[1]);
    const auto dist_min = std::min(dist[0], dist[1]);
    const auto dot = dist[0] * dist[1];

    VectorDouble6 coords(6);
    CHKERR moab.get_coords(conn, 2, &coords[0]);
    VectorAdaptor n0 = getVectorAdaptor(&coords[0], 3);
    VectorAdaptor n1 = getVectorAdaptor(&coords[3], 3);
    VectorDouble3 ray = n1 - n0;
    const double ray_length = norm_2(ray);
    ray /= ray_length;

    if (

        (dot < 0 && dist_max > std::numeric_limits<double>::epsilon()) ||
        (std::abs(dist_min) < std::numeric_limits<double>::epsilon() &&
         dist_max > std::numeric_limits<double>::epsilon())

    ) {

      // Edges is on two sides of the surface
      const double s = dist[0] / (dist[0] - dist[1]);
      const double dist = s * ray_length;

      auto add_edge = [&](auto dist) {
        edgesToCut[e].dIst = dist;
        edgesToCut[e].lEngth = ray_length;
        edgesToCut[e].unitRayDir = ray;
        edgesToCut[e].rayPoint = n0;
        cutEdges.insert(e);

        aveLength += norm_2(ray);
        maxLength = fmax(maxLength, norm_2(ray));
        ++nb_ave_length;
      };

      add_edge(dist);
    }
  }
  aveLength /= nb_ave_length;

  if (debug)
    CHKERR SaveData(m_field.get_moab())("cut_edges.vtk", cutEdges);

  if (debug)
    CHKERR SaveData(m_field.get_moab())("cut_edges_vol.vtk",
                                        get_adj(cutEdges, 3));

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::projectZeroDistanceEnts(
    Range *fixed_edges, Range *corner_nodes, const double geometry_tol,
    const double close_tol, const int verb, const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  Skinner skin(&moab);
  MoFEMFunctionBegin;

  auto get_ent_adj = [&moab](const EntityHandle v, const int dim) {
    Range a;
    if (dim)
      CHKERR moab.get_adjacencies(&v, 1, dim, true, a);
    return a;
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

  auto get_conn = [&](auto e) {
    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(e, conn, num_nodes, true);
    return conn;
  };

  auto get_range = [](std::vector<EntityHandle> v) {
    Range r;
    r.insert_list(v.begin(), v.end());
    return r;
  };

  auto get_coords = [&](const auto v) {
    VectorDouble3 coords(3);
    CHKERR moab.get_coords(&v, 1, &coords[0]);
    return coords;
  };

  Tag th_normal_dist;
  CHKERR moab.tag_get_handle("DIST_SURFACE_NORMAL_VECTOR", th_normal_dist);

  auto get_normal_dits = [&](const EntityHandle v) {
    VectorDouble3 dist_vec(3);
    CHKERR moab.tag_get_data(th_normal_dist, &v, 1, &*dist_vec.begin());
    return dist_vec;
  };

  auto get_prj_point = [&](const EntityHandle v, const Range edges,
                           const bool geom_dist, const double geometry_tol) {
    auto get_tuple = [](const EntityHandle e, const double dist,
                        const double l) { return std::make_tuple(e, dist, l); };

    std::tuple<EntityHandle, double, double> min_tuple{0, 0, 0};

    for (auto e : edges) {
      auto eit = edgesToCut.find(e);
      if (eit != edgesToCut.end()) {

        auto get_dist = [&](auto eit) {
          int num_nodes;
          const EntityHandle *conn;
          CHKERR moab.get_connectivity(eit->first, conn, num_nodes, true);
          if (conn[0] == v)
            return eit->second.dIst;
          else if (conn[1] == v)
            return (eit->second.lEngth - eit->second.dIst);
          else
            THROW_MESSAGE("Data inconsistency");
        };

        const auto d = get_dist(eit);
        if (std::get<0>(min_tuple) == 0) {
          min_tuple = get_tuple(e, d, eit->second.lEngth);
        } else if (std::get<1>(min_tuple) > d) {
          min_tuple = get_tuple(e, d, eit->second.lEngth);
        }
      }
    }

    const auto geom_dist_vec = get_normal_dits(v);
    const auto geom_tol = norm_2(geom_dist_vec);
    const auto l = std::get<2>(min_tuple);

    if (geom_tol < l * geometry_tol) {

      return std::make_pair(get_coords(v), l);
    }

    if (geom_dist) {

      return std::make_pair(VectorDouble3(get_coords(v) + geom_dist_vec), l);

    } else {

      const auto &d = edgesToCut.at(std::get<0>(min_tuple));
      return std::make_pair(VectorDouble3(d.rayPoint + d.dIst * d.unitRayDir),
                            l);
    }
  };

  auto get_in_range = [](auto v, auto &r) { return (r.find(v) != r.end()); };

  auto project_nodes = [&](auto nodes_to_check) {
    auto get_fix_e = [](auto fixed_edges) {
      if (fixed_edges)
        return *fixed_edges;
      else
        return Range();
    };

    const auto fix_e = get_fix_e(fixed_edges);
    const auto skin_e = get_adj(unite(get_skin(vOlume), constrainSurface), 1);
    const auto cut_fix = intersect(cutEdges, fix_e);
    const auto cut_skin = intersect(cutEdges, skin_e);

    Range corner_n;
    if (corner_nodes)
      corner_n = intersect(*corner_nodes, nodes_to_check);
    auto fixe_n = intersect(get_adj(fix_e, 0), nodes_to_check);
    auto skin_n = intersect(get_adj(skin_e, 0), nodes_to_check);

    std::vector<std::pair<EntityHandle, TreeData>> vertices_on_cut_edges;
    vertices_on_cut_edges.reserve(nodes_to_check.size());

    auto add_zero_vertex = [&](const EntityHandle v, auto &&new_pos) {
      auto coords = get_coords(v);
      auto ray = new_pos - coords;
      auto l = norm_2(ray);
      TreeData data;
      data.dIst = l;
      if (l > std::numeric_limits<double>::epsilon()) {
        data.unitRayDir = ray / l;
      } else {
        data.dIst = 0;
        data.unitRayDir = ray;
      }
      data.rayPoint = coords;
      return std::make_pair(v, data);
    };

    auto intersect_v = [&](const auto v, const auto r) {
      return intersect(r, get_ent_adj(v, 1));
    };

    for (auto v : nodes_to_check) {
      if (get_in_range(v, corner_n)) {
        vertices_on_cut_edges.push_back(add_zero_vertex(v, get_coords(v)));
      } else if (get_in_range(v, fixe_n)) {

        const auto e = intersect_v(v, cut_fix);
        if (!e.empty()) {
          vertices_on_cut_edges.push_back(add_zero_vertex(
              v, get_prj_point(v, e, false, geometry_tol).first));

        } else {
          auto b = intersect_v(v, cutEdges);
          if (!b.empty()) {
            auto p = get_prj_point(v, b, true, geometry_tol);
            if (norm_2(get_coords(v) - p.first) < geometry_tol * p.second)
              vertices_on_cut_edges.push_back(add_zero_vertex(v, p.first));
          }
        }

      } else if (get_in_range(v, skin_n)) {

        const auto e = intersect_v(v, cut_skin);
        if (!e.empty()) {
          vertices_on_cut_edges.push_back(add_zero_vertex(
              v, get_prj_point(v, e, false, geometry_tol).first));

        } else {
          auto b = intersect_v(v, cutEdges);
          if (!b.empty()) {
            auto p = get_prj_point(v, b, true, geometry_tol);
            if (norm_2(get_coords(v) - p.first) < geometry_tol * p.second)
              vertices_on_cut_edges.push_back(add_zero_vertex(v, p.first));
          }
        }

      } else {
        const auto e = intersect_v(v, cutEdges);
        if (!e.empty())
          vertices_on_cut_edges.push_back(
              add_zero_vertex(v, get_prj_point(v, e, false, 0).first));
      }
    }

    auto get_distances = [&](auto &data) {
      Tag th;
      CHKERR moab.tag_get_handle("DIST_SURFACE_NORMAL_SIGNED", th);
      std::vector<EntityHandle> verts;
      verts.reserve(data.size());
      for (auto d : data)
        verts.emplace_back(d.first);
      std::vector<double> distances(verts.size());
      CHKERR moab.tag_get_data(th, &*verts.begin(), verts.size(),
                               &*distances.begin());
      std::map<EntityHandle, double> dist_map;
      for (size_t i = 0; i != verts.size(); ++i)
        dist_map[verts[i]] = distances[i];
      return dist_map;
    };

    auto dist_map = get_distances(vertices_on_cut_edges);

    auto cmp = [&dist_map](const auto &a, const auto &b) {
      const double da = dist_map[a.first];
      const double db = dist_map[b.first];
      if (da * db < 0)
        return dist_map[a.first] < dist_map[b.first];
      else
        return std::abs(dist_map[a.first]) < std::abs(dist_map[b.first]);
    };

    std::sort(vertices_on_cut_edges.begin(), vertices_on_cut_edges.end(), cmp);

    auto cmp_abs = [&dist_map](const auto &a, const auto &b) {
      return std::abs(dist_map[a.first]) < std::abs(dist_map[b.first]);
    };
    auto max_dist = std::max_element(vertices_on_cut_edges.begin(),
                                     vertices_on_cut_edges.end(), cmp_abs)
                        ->second.dIst;

    decltype(vertices_on_cut_edges) vertices_on_cut_edges_tol;
    vertices_on_cut_edges_tol.reserve(vertices_on_cut_edges.size());
    for (auto &t : vertices_on_cut_edges)
      if (t.second.dIst < close_tol * max_dist)
        vertices_on_cut_edges_tol.emplace_back(t);

    return vertices_on_cut_edges_tol;
  };

  auto get_min_quality =
      [&](const Range &adj_tets,
          const map<EntityHandle, TreeData> &vertices_on_cut_edges,
          const std::pair<EntityHandle, TreeData> *test_data_ptr) {
        double min_q = 1;
        for (auto t : adj_tets) {
          int num_nodes;
          const EntityHandle *conn;
          CHKERR m_field.get_moab().get_connectivity(t, conn, num_nodes, true);
          VectorDouble12 coords(12);
          CHKERR moab.get_coords(conn, num_nodes, &coords[0]);

          auto set_new_coord = [](auto d) {
            return d.rayPoint + d.dIst * d.unitRayDir;
          };

          for (auto n : {0, 1, 2, 3}) {
            auto n_coords = getVectorAdaptor(&coords[3 * n], 3);
            if (test_data_ptr) {
              if (test_data_ptr->first == conn[n]) {
                noalias(n_coords) = set_new_coord(test_data_ptr->second);
              }
            }
            auto mit = vertices_on_cut_edges.find(conn[n]);
            if (mit != vertices_on_cut_edges.end()) {
              noalias(n_coords) = set_new_coord(mit->second);
            }
          }
          min_q = std::min(min_q, Tools::volumeLengthQuality(&coords[0]));
        }
        return min_q;
      };

  auto get_zero_distance_verts = [&](const auto &&vertices_on_cut_edges) {
    verticesOnCutEdges.clear();
    std::vector<EntityHandle> zero_dist_vec;
    zero_dist_vec.reserve(vertices_on_cut_edges.size());
    for (auto t : vertices_on_cut_edges) {

      auto adj_tet = intersect(vOlume, get_ent_adj(t.first, 3));
      const double q0 = get_min_quality(adj_tet, verticesOnCutEdges, nullptr);
      const double q = get_min_quality(adj_tet, verticesOnCutEdges, &t);
      if (std::isnormal(q) && (q / q0) > projectEntitiesQualityTrashold) {
        verticesOnCutEdges[t.first] = t.second;
        zero_dist_vec.push_back(t.first);
      }
    }
    return zero_dist_vec;
  };

  auto vol_cut_ents = intersect(vOlume, get_adj(cutEdges, 3));

  auto get_zero_distant_ents = [&](auto bridge_ents, const int dim,
                                   const int bridge_dim) {
    auto ents =
        intersect(get_adj(bridge_ents, dim), get_adj(vol_cut_ents, dim));
    auto ents_to_remove =
        get_adj(subtract(get_adj(ents, bridge_dim), bridge_ents), dim);
    return subtract(ents, ents_to_remove);
  };

  zeroDistanceVerts =
      get_range(get_zero_distance_verts(project_nodes(get_adj(cutEdges, 0))));
  zeroDistanceEnts = subtract(get_zero_distant_ents(zeroDistanceVerts, 2, 0),
                              get_skin(vOlume));

  if (debug)
    CHKERR SaveData(moab)("zero_distance_verts.vtk", zeroDistanceVerts);
  if (debug)
    CHKERR SaveData(moab)("zero_distance_ents.vtk", zeroDistanceEnts);

  for (auto f : zeroDistanceVerts) {
    auto adj_edges = get_ent_adj(f, 1);
    for (auto e : adj_edges) {
      cutEdges.erase(e);
      edgesToCut.erase(e);
    }
  }

  if (debug)
    CHKERR SaveData(moab)("cut_edges_to_cut.vtk", cutEdges);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::cutEdgesInMiddle(const BitRefLevel bit,
                                                  Range &cut_vols,
                                                  Range &cut_surf,
                                                  Range &cut_verts,
                                                  const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MeshRefinement *refiner;
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBegin;

  if (cutEdges.size() != edgesToCut.size())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Data inconsistency");

  auto refine_mesh = [&]() {
    MoFEMFunctionBegin;
    CHKERR m_field.getInterface(refiner);
    CHKERR m_field.get_ref_ents(&ref_ents_ptr);
    CHKERR refiner->add_vertices_in_the_middle_of_edges(cutEdges, bit);
    CHKERR refiner->refine_TET(vOlume, bit, false, QUIET,
                               debug ? &cutEdges : NULL);
    MoFEMFunctionReturn(0);
  };

  CHKERR refine_mesh();

  if (debug)
    if (cutEdges.size() != edgesToCut.size())
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Data inconsistency");

  cut_vols.clear();
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, BitRefLevel().set(), MBTET, cut_vols);
  cut_surf.clear();
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, bit, MBTRI, cut_surf);

  if (debug)
    CHKERR SaveData(moab)("cut_surf_from_bit.vtk", cut_surf);

  // Find new vertices on cut edges
  cut_verts.clear();
  cut_verts.merge(zeroDistanceVerts);

  for (auto &m : edgesToCut) {
    auto vit = ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>().find(
        boost::make_tuple(MBVERTEX, m.first));
    if (vit ==
        ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>().end()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "No vertex on cut edges, that make no sense");
    }
    const boost::shared_ptr<RefEntity> &ref_ent = *vit;
    if ((ref_ent->getBitRefLevel() & bit).any()) {
      EntityHandle vert = ref_ent->getRefEnt();
      cut_verts.insert(vert);
      verticesOnCutEdges[vert] = m.second;
    } else {
      SETERRQ1(
          PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
          "Vertex has wrong bit ref level %s",
          boost::lexical_cast<std::string>(ref_ent->getBitRefLevel()).c_str());
    }
  }

  // Add zero distance entities faces
  Range tets_skin;
  Skinner skin(&moab);
  CHKERR skin.find_skin(0, cut_vols, false, tets_skin);
  tets_skin.merge(constrainSurface);

  // At that point cut_surf has all newly created faces, now take all
  // nodes on those faces and subtract nodes on cut edges. Faces adjacent to
  // nodes which left are not part of surface.
  cut_surf.merge(zeroDistanceEnts.subset_by_type(MBTRI));
  Range diff_verts;
  CHKERR moab.get_connectivity(cut_surf, diff_verts, true);
  diff_verts = subtract(diff_verts, cut_verts);

  Range subtract_faces;
  CHKERR moab.get_adjacencies(diff_verts, 2, false, subtract_faces,
                              moab::Interface::UNION);
  cut_surf = subtract(cut_surf, unite(subtract_faces, tets_skin));

  cut_verts.clear();
  CHKERR moab.get_connectivity(cut_surf, cut_verts, true);

  // Check non-mainfolds
  auto check_for_non_minfold = [&]() {
    MoFEMFunctionBegin;
    Range surf_edges;
    CHKERR moab.get_adjacencies(cut_surf, 1, false, surf_edges,
                                moab::Interface::UNION);
    for (auto e : surf_edges) {

      Range faces;
      CHKERR moab.get_adjacencies(&e, 1, 2, false, faces);
      faces = intersect(faces, cut_surf);
      if (faces.size() > 2) {

        bool resolved = false;

        // Check for haning node
        Range nodes;
        CHKERR moab.get_connectivity(faces, nodes, true);
        for (auto n : nodes) {
          Range adj_faces;
          CHKERR moab.get_adjacencies(&n, 1, 2, false, adj_faces);
          adj_faces = intersect(adj_faces, cut_surf);
          if (adj_faces.size() == 1) {
            cut_surf.erase(adj_faces[0]);
            resolved = true;
          }
        }

        // Check for two edges minfold
        Range adj_edges;
        CHKERR moab.get_adjacencies(faces, 1, false, adj_edges,
                                    moab::Interface::UNION);
        adj_edges = intersect(adj_edges, surf_edges);
        adj_edges.erase(e);
        for (auto other_e : adj_edges) {
          Range other_faces;
          CHKERR moab.get_adjacencies(&other_e, 1, 2, false, other_faces);
          other_faces = intersect(other_faces, cut_surf);
          if (other_faces.size() > 2) {
            other_faces = intersect(other_faces, faces);
            cut_surf = subtract(cut_surf, other_faces);
            resolved = true;
          }
        }

        if (!resolved && !debug)
          SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
                  "Non-mainfold surface");

        cut_verts.clear();
        CHKERR moab.get_connectivity(cut_surf, cut_verts, true);
      }
    }
    MoFEMFunctionReturn(0);
  };

  CHKERR check_for_non_minfold();

  if (debug)
    CHKERR SaveData(moab)("cut_surf.vtk", cut_surf);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::moveMidNodesOnCutEdges(Tag th) {
  MoFEMFunctionBeginHot;

  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  // Range out_side_vertices;
  for (auto m : verticesOnCutEdges) {
    double dist = m.second.dIst;
    VectorDouble3 new_coors = m.second.rayPoint + dist * m.second.unitRayDir;
    if (th)
      CHKERR moab.tag_set_data(th, &m.first, 1, &new_coors[0]);
    else
      CHKERR moab.set_coords(&m.first, 1, &new_coors[0]);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::moveMidNodesOnTrimmedEdges(Tag th) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  for (auto &v : verticesOnTrimEdges) {
    double dist = v.second.dIst;
    VectorDouble3 new_coors = v.second.rayPoint + dist * v.second.unitRayDir;
    if (th)
      CHKERR moab.tag_set_data(th, &v.first, 1, &new_coors[0]);
    else
      CHKERR moab.set_coords(&v.first, 1, &new_coors[0]);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::findEdgesToTrim(Range *fixed_edges,
                                                 Range *corner_nodes, Tag th,
                                                 const double tol_geometry,
                                                 const double tol_trim_close,
                                                 const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  // takes body skin
  Skinner skin(&moab);
  Range tets_skin;
  CHKERR skin.find_skin(0, cutNewVolumes, false, tets_skin);
  tets_skin.merge(constrainSurface);

  // vertices on the skin
  Range tets_skin_verts;
  CHKERR moab.get_connectivity(tets_skin, tets_skin_verts, true);
  // edges on the skin
  Range tets_skin_edges;
  CHKERR moab.get_adjacencies(tets_skin, 1, false, tets_skin_edges,
                              moab::Interface::UNION);
  // get edges on new surface
  Range cut_surface_edges;
  CHKERR moab.get_adjacencies(cutNewSurfaces, 1, false, cut_surface_edges,
                              moab::Interface::UNION);
  Range cut_surface_verts;
  CHKERR moab.get_connectivity(cutNewSurfaces, cut_surface_verts, true);

  Range corners;
  if (corner_nodes)
    corners = *corner_nodes;

  Range fix_edges;
  if (fixed_edges)
    fix_edges = *fixed_edges;

  Range fixed_edges_verts;
  if (fixed_edges)
    CHKERR moab.get_connectivity(*fixed_edges, fixed_edges_verts, true);

  Range surface_skin;
  if (fRont.empty())
    CHKERR skin.find_skin(0, sUrface, false, surface_skin);
  else
    surface_skin = fRont;

  auto get_point_coords = [&](EntityHandle v) {
    VectorDouble3 point(3);
    if (th)
      CHKERR moab.tag_get_data(th, &v, 1, &point[0]);
    else
      CHKERR moab.get_coords(&v, 1, &point[0]);
    return point;
  };

  struct VertMap {
    const double d;
    const EntityHandle v;
    const EntityHandle e;
    VertMap(const double d, const EntityHandle v, const EntityHandle e)
        : d(d), v(v), e(e) {}
  };

  typedef multi_index_container<
      VertMap,
      indexed_by<
          ordered_non_unique<member<VertMap, const double, &VertMap::d>>,
          ordered_non_unique<member<VertMap, const EntityHandle, &VertMap::v>>,
          ordered_non_unique<member<VertMap, const EntityHandle, &VertMap::e>>

          >>
      VertMapMultiIndex;

  VertMapMultiIndex verts_map;

  auto add_vert = [&](const EntityHandle v, const EntityHandle e,
                      const double dist) {
    verts_map.insert(VertMap(dist, v, e));
  };

  // clear data ranges
  trimEdges.clear();
  edgesToTrim.clear();
  verticesOnTrimEdges.clear();
  trimNewVertices.clear();

  auto cut_this_edge = [&](const EntityHandle e, const double length, auto &ray,
                           auto &ray_point) {
    trimEdges.insert(e);
    edgesToTrim[e].dIst = 1;
    edgesToTrim[e].lEngth = length;
    edgesToTrim[e].unitRayDir.resize(3, false);
    edgesToTrim[e].rayPoint.resize(3, false);
    for (int ii = 0; ii != 3; ++ii)
      edgesToTrim[e].unitRayDir[ii] = ray(ii);
    for (int ii = 0; ii != 3; ++ii)
      edgesToTrim[e].rayPoint[ii] = ray_point(ii);
  };

  FTensor::Index<'i', 3> i;
  int num_nodes;

  auto get_tag_data = [&](auto th, auto conn) {
    FTensor::Tensor1<double, 3> t;
    CHKERR moab.tag_get_data(th, &conn, 1, &t(0));
    return t;
  };

  double max_edge_length = 0;

  /// Project front entities on on the cut surface plane
  if (!fRont.empty()) {
    // Calculate distances
    treeSurfPtr = boost::shared_ptr<OrientedBoxTreeTool>(
        new OrientedBoxTreeTool(&moab, "ROOTSETSURF", true));
    CHKERR treeSurfPtr->build(cutNewSurfaces, rootSetSurf);

    for (auto s : surface_skin) {

      auto project_on_surface = [&](auto &t) {
        FTensor::Tensor1<double, 3> t_p;

        EntityHandle facets_out;
        CHKERR treeSurfPtr->closest_to_location(&t(0), rootSetSurf, &t_p(0),
                                                facets_out);

        FTensor::Tensor1<double, 3> t_normal;
        CHKERR m_field.getInterface<Tools>()->getTriNormal(facets_out,
                                                           &t_normal(0));
        t_normal(i) /= sqrt(t_normal(i) * t_normal(i));
        const double dot = t_normal(i) * (t_p(i) - t(i));
        t_p(i) = t(i) + dot * t_normal(i);

        return t_p;
      };

      const EntityHandle *conn;
      CHKERR moab.get_connectivity(s, conn, num_nodes, true);

      VectorDouble6 coords_front(6);
      CHKERR moab.get_coords(conn, num_nodes, &coords_front[0]);

      FTensor::Tensor1<double *, 3> t_f0(&coords_front[0], &coords_front[1],
                                         &coords_front[2]);
      FTensor::Tensor1<double *, 3> t_f1(&coords_front[3], &coords_front[4],
                                         &coords_front[5]);

      auto t_p_f0 = project_on_surface(t_f0);
      auto t_p_f1 = project_on_surface(t_f1);

      CHKERR moab.set_coords(&conn[0], 1, &t_p_f0(0));
      CHKERR moab.set_coords(&conn[1], 1, &t_p_f1(0));
    }
  }

  if (debug)
    CHKERR SaveData(moab)("surface_skin_to_trim.vtk", surface_skin);

  CHKERR createFrontLevelSets(cutNewSurfaces, th, QUIET, debug);
  Tag th_dist_front_vec;
  CHKERR moab.tag_get_handle("DIST_FRONT_VECTOR", th_dist_front_vec);

  auto create_tag = [&](const std::string name, const int dim) {
    Tag th;
    rval = moab.tag_get_handle(name.c_str(), th);
    if (rval == MB_SUCCESS)
      return th;
    std::vector<double> def_val(dim, 0);
    CHKERR moab.tag_get_handle(name.c_str(), dim, MB_TYPE_DOUBLE, th,
                               MB_TAG_CREAT | MB_TAG_SPARSE, &*def_val.begin());

    return th;
  };

  // iterate over edges on cut surface
  for (auto e : cut_surface_edges) {

    auto get_conn_edge = [&moab](auto e) {
      // Get edge connectivity and coords
      const EntityHandle *conn_edge;
      int num_nodes;
      CHKERR moab.get_connectivity(e, conn_edge, num_nodes, true);
      return conn_edge;
    };

    auto get_coords_edge = [&moab](auto conn_edge) {
      std::array<double,6> coords_edge;
      CHKERR moab.get_coords(conn_edge, 2, coords_edge.data());
      return coords_edge;
    };

    auto get_ftensor_coords = [](const double *ptr) {
      return FTensor::Tensor1<double, 3>{ptr[0], ptr[1], ptr[2]};
    };

    auto conn_edge = get_conn_edge(e);
    auto coords_edge = get_coords_edge(conn_edge);

    auto t_dist0 = get_tag_data(th_dist_front_vec, conn_edge[0]);
    auto t_dist1 = get_tag_data(th_dist_front_vec, conn_edge[1]);
    auto t_e0 = get_ftensor_coords(&coords_edge[0]);
    auto t_e1 = get_ftensor_coords(&coords_edge[3]);

    FTensor::Tensor1<double, 3> t_edge_delta;
    t_edge_delta(i) = t_e1(i) - t_e0(i);
    const double edge_length2 = t_edge_delta(i) * t_edge_delta(i);
    const double edge_length = sqrt(edge_length2);
    if (edge_length == 0)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Zero edge length");

    max_edge_length = std::max(max_edge_length, edge_length);
    const double dot = t_dist0(i) * t_dist1(i);
    const double dot_direction0 = t_dist0(i) * t_edge_delta(i);
    const double dot_direction1 = t_dist1(i) * t_edge_delta(i);

    auto add_edge = [&](auto t_cut) {
      FTensor::Tensor1<double, 3> t_ray;
      if (t_cut < 0.5) {
        t_ray(i) = t_cut * t_edge_delta(i);
        add_vert(conn_edge[0], e, fabs(t_cut));
        add_vert(conn_edge[1], e, fabs(t_cut - 1));
        cut_this_edge(e, edge_length, t_ray, t_e0);
      } else {
        FTensor::Tensor1<double, 3> t_edge_point;
        t_edge_point(i) = t_e0(i) + t_cut * t_edge_delta(i);
        t_ray(i) = t_edge_point(i) - t_e1(i);
        add_vert(conn_edge[0], e, fabs(t_cut));
        add_vert(conn_edge[1], e, fabs(t_cut - 1));
        cut_this_edge(e, edge_length, t_ray, t_e1);
      }
    };

    if (dot < tol_geometry * edge_length &&
        dot_direction0 > -std::numeric_limits<double>::epsilon() &&
        dot_direction1 < std::numeric_limits<double>::epsilon()) {

      const double s0 = t_dist0(i) * t_edge_delta(i) / edge_length;
      const double s1 = t_dist1(i) * t_edge_delta(i) / edge_length;
      const double t_cut = s0 / (s0 - s1);

      add_edge(t_cut);

    } else if (std::min(dot_direction0, dot_direction1) < edge_length) {

      for (auto f : surface_skin) {

        auto conn_front = get_conn_edge(f);
        auto coords_front = get_coords_edge(conn_front);
        auto t_f0 = get_ftensor_coords(&coords_front[0]);
        auto t_f1 = get_ftensor_coords(&coords_front[3]);

        double tvw;
        double tlk;

        auto res = Tools::minDistanceFromSegments(&t_e0(0), &t_e1(0), &t_f0(0),
                                                  &t_f1(0), &tvw, &tlk);

        if (res != Tools::NO_SOLUTION) {

          auto check = [](auto v) {
            return v > -std::numeric_limits<double>::epsilon() &&
                   (v - 1) < std::numeric_limits<double>::epsilon();
          };

          if ( check(tvw) && check(tlk) ) {
            add_edge(tvw);
          }

        }

      }
    }
  }


  if (debug)
    CHKERR SaveData(moab)("edges_potentially_to_trim.vtk", cut_surface_edges);

  if (debug)
    CHKERR SaveData(moab)("edges_selected_to_trim.vtk", trimEdges);

  auto get_quality_change = [&](const Range &adj_tets, const EntityHandle &v,
                                const VectorDouble3 &new_pos) {
    double q0 = 1;
    CHKERR m_field.getInterface<Tools>()->minTetsQuality(adj_tets, q0, th);

    double min_q = 1;
    for (auto t : adj_tets) {
      int num_nodes;
      const EntityHandle *conn;
      CHKERR m_field.get_moab().get_connectivity(t, conn, num_nodes, true);
      VectorDouble12 coords(12);
      if (th)
        CHKERR moab.tag_get_data(th, conn, num_nodes, &coords[0]);
      else
        CHKERR moab.get_coords(conn, num_nodes, &coords[0]);

      for (int n = 0; n != 4; ++n) {
        auto n_coords = getVectorAdaptor(&coords[3 * n], 3);
        if (conn[n] == v) {
          noalias(n_coords) = new_pos;
        } else {
          auto m = verticesOnTrimEdges.find(conn[n]);
          if (m != verticesOnTrimEdges.end())
            noalias(n_coords) =
                m->second.rayPoint + m->second.dIst * m->second.unitRayDir;
        }
      }

      const double q = Tools::volumeLengthQuality(&coords[0]);
      if (!std::isnormal(q)) {
        min_q = -2;
        break;
      }
      min_q = std::min(min_q, q);
    }
    return min_q / q0;
  };

  Range checked_verts;
  auto add_trim_vert = [&](const EntityHandle v, const EntityHandle e) {
    if (!(trimNewVertices.find(v) != trimNewVertices.end())) {
      auto &r = edgesToTrim.at(e);
      VectorDouble3 ray_point = get_point_coords(v);
      VectorDouble3 new_pos = r.rayPoint + r.dIst * r.unitRayDir;
      VectorDouble3 unit_ray_dir = ray_point - new_pos;
      Range adj_tets;
      CHKERR moab.get_adjacencies(&v, 1, 3, false, adj_tets);
      adj_tets = intersect(adj_tets, cutNewVolumes);
      double q = get_quality_change(adj_tets, v, new_pos);
      if (q > projectEntitiesQualityTrashold) {
        VectorDouble3 ray_dir = new_pos - ray_point;
        double dist = norm_2(ray_dir);
        verticesOnTrimEdges[v].dIst = 1;
        verticesOnTrimEdges[v].lEngth = dist;
        verticesOnTrimEdges[v].unitRayDir = ray_dir;
        verticesOnTrimEdges[v].rayPoint = ray_point;
        trimNewVertices.insert(v);
      }
      checked_verts.insert(v);
    }
  };

  auto add_no_move_trim = [&](const EntityHandle v, const EntityHandle e) {
    if (!(trimNewVertices.find(v) != trimNewVertices.end())) {
      auto &m = edgesToTrim.at(e);
      verticesOnTrimEdges[v] = m;
      verticesOnTrimEdges[v].rayPoint = get_point_coords(v);
      verticesOnTrimEdges[v].lEngth = 0;
      verticesOnTrimEdges[v].dIst = 0;
      trimNewVertices.insert(v);
      checked_verts.insert(v);
    }
  };

  // Iterate over all vertives close to surface front and check if those can
  // be moved

  // filer nodes which distance is in given tolerance
  auto hi = verts_map.get<0>().upper_bound(tol_trim_close);
  verts_map.get<0>().erase(hi, verts_map.end());

  auto remove_verts = [&](Range nodes) {
    for (auto n : nodes) {
      auto r = verts_map.get<1>().equal_range(n);
      verts_map.get<1>().erase(r.first, r.second);
    }
  };

  auto trim_verts = [&](const Range verts, const bool move) {
    VertMapMultiIndex verts_map_tmp;
    for (auto p = corners.pair_begin(); p != corners.pair_end(); ++p) {
      auto lo = verts_map.get<1>().lower_bound(p->first);
      auto hi = verts_map.get<1>().upper_bound(p->second);
      verts_map_tmp.insert(lo, hi);
    }
    if (move) {
      for (auto &m : verts_map_tmp.get<0>())
        add_trim_vert(m.v, m.e);
    } else {
      for (auto &m : verts_map_tmp.get<0>())
        add_no_move_trim(m.v, m.e);
    }
  };

  auto trim_edges = [&](const Range ents, const bool move) {
    VertMapMultiIndex verts_map_tmp;
    for (auto p = ents.pair_begin(); p != ents.pair_end(); ++p) {
      auto lo = verts_map.get<2>().lower_bound(p->first);
      auto hi = verts_map.get<2>().upper_bound(p->second);
      verts_map_tmp.insert(lo, hi);
      verts_map.get<2>().erase(lo, hi);
    }
    if (move) {
      for (auto &m : verts_map_tmp.get<0>())
        add_trim_vert(m.v, m.e);
    } else {
      for (auto &m : verts_map_tmp.get<0>())
        add_no_move_trim(m.v, m.e);
    }
  };

  auto intersect_self = [&](Range &a, const Range b) { a = intersect(a, b); };

  map<std::string, Range> range_maps;
  CHKERR skin.find_skin(0, cutNewSurfaces, false, range_maps["surface_skin"]);
  intersect_self(range_maps["surface_skin"], trimEdges);

  range_maps["fixed_edges_on_surface_skin"] =
      intersect(range_maps["surface_skin"], fix_edges);

  CHKERR moab.get_adjacencies(fixed_edges_verts, 1, false,
                              range_maps["fixed_edges_verts_edges"],
                              moab::Interface::UNION);
  intersect_self(range_maps["fixed_edges_verts_edges"], trimEdges);
  CHKERR moab.get_connectivity(
      range_maps["fixed_edges_verts_edges"],
      range_maps["fixed_edges_verts_edges_verts_on_the_skin"], false);
  intersect_self(range_maps["fixed_edges_verts_edges_verts_on_the_skin"],
                 tets_skin_verts);

  // do not move nodes at the corners
  trim_verts(corners, false);
  remove_verts(corners);
  trim_edges(range_maps["fixed_edges_on_surface_skin"], true);
  remove_verts(range_maps["fixed_edges_on_surface_skin_verts"]);
  trim_verts(range_maps["fixed_edges_verts_edges_verts_on_the_skin"], false);
  remove_verts(range_maps["fixed_edges_verts_edges_verts_on_the_skin"]);
  trim_edges(range_maps["surface_skin"], true);
  trim_verts(tets_skin_verts, false);
  remove_verts(tets_skin_verts);

  for (auto &m : verts_map.get<0>())
    add_trim_vert(m.v, m.e);

  auto add_zero_dist_vertices = [&](const auto &&verts, const double geom_tol) {
    MoFEMFunctionBegin;

    std::vector<double> dits_vec(3*verts.size());
    CHKERR moab.tag_get_data(th_dist_front_vec, verts, &*dits_vec.begin());
    auto t_dist = FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
        &dits_vec[0], &dits_vec[1], &dits_vec[2]);

    for (auto v : verts) {

      if (trimNewVertices.find(v) == trimNewVertices.end()) {

        auto get_tag_data = [&](auto th, auto conn) {
          FTensor::Tensor1<double, 3> t;
          CHKERR moab.tag_get_data(th, &conn, 1, &t(0));
          return t;
        };

        auto t_dist = get_tag_data(th_dist_front_vec, v);
        const double d = sqrt(t_dist(i) * t_dist(i));
        if (d < geom_tol * max_edge_length) {
          verticesOnTrimEdges[v].dIst = 0;
          verticesOnTrimEdges[v].lEngth = 0;
          verticesOnTrimEdges[v].unitRayDir.resize(3);
          verticesOnTrimEdges[v].unitRayDir.clear();
          verticesOnTrimEdges[v].rayPoint = get_point_coords(v);
          trimNewVertices.insert(v);

        }
      }

      ++t_dist;
    }

    MoFEMFunctionReturn(0);
  };

  CHKERR add_zero_dist_vertices(subtract(cut_surface_verts, checked_verts),
                                tol_geometry);

  for (auto m : verticesOnTrimEdges) {
    EntityHandle v = m.first;
    Range adj;
    CHKERR moab.get_adjacencies(&v, 1, 1, false, adj);
    for (auto e : adj) {
      edgesToTrim.erase(e);
      trimEdges.erase(e);
    }
  }

  if (debug)
    if (!trimNewVertices.empty())
      CHKERR SaveData(moab)("trim_close_vertices.vtk", trimNewVertices);

  if (debug)
    if (!trimEdges.empty())
      CHKERR SaveData(moab)("trim_edges.vtk", trimEdges);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::trimEdgesInTheMiddle(const BitRefLevel bit,
                                                      const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MeshRefinement *refiner;
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBegin;

  CHKERR m_field.getInterface(refiner);
  CHKERR m_field.get_ref_ents(&ref_ents_ptr);
  CHKERR refiner->add_vertices_in_the_middle_of_edges(trimEdges, bit);
  CHKERR refiner->refine_TET(cutNewVolumes, bit, false, QUIET,
                             debug ? &trimEdges : NULL);

  trimNewVolumes.clear();
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, BitRefLevel().set(), MBTET, trimNewVolumes);

  for (map<EntityHandle, TreeData>::iterator mit = edgesToTrim.begin();
       mit != edgesToTrim.end(); mit++) {
    auto vit = ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>().find(
        boost::make_tuple(MBVERTEX, mit->first));
    if (vit ==
        ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>().end())
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "No vertex on trim edges, that make no sense");

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

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::trimSurface(Range *fixed_edges,
                                             Range *corner_nodes,
                                             const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  Skinner skin(&moab);
  Range trim_tets_skin;
  CHKERR skin.find_skin(0, trimNewVolumes, false, trim_tets_skin);
  Range trim_tets_skin_edges;
  CHKERR moab.get_adjacencies(trim_tets_skin, 1, false, trim_tets_skin_edges,
                              moab::Interface::UNION);

  Range barrier_vertices(trimNewVertices);

  if (fixed_edges && trimFixedEdges) {

    // get all vertices on fixed edges and surface
    Range trim_surface_edges;
    CHKERR moab.get_adjacencies(trimNewSurfaces, 1, false, trim_surface_edges,
                                moab::Interface::UNION);
    Range fixed_edges_on_trim_surface;
    fixed_edges_on_trim_surface = intersect(trim_surface_edges, *fixed_edges);
    Range fixed_edges_on_trim_surface_verts;
    CHKERR moab.get_connectivity(fixed_edges_on_trim_surface,
                                 fixed_edges_on_trim_surface_verts, false);

    // get faces adjacent to barrier_vertices
    Range barrier_vertices_faces;
    CHKERR moab.get_adjacencies(barrier_vertices, 2, false,
                                barrier_vertices_faces, moab::Interface::UNION);
    barrier_vertices_faces = intersect(barrier_vertices_faces, trimNewSurfaces);

    // get vertices on fixed edges
    Range fixed_edges_vertices;
    CHKERR moab.get_connectivity(*fixed_edges, fixed_edges_vertices, false);
    fixed_edges_vertices = intersect(barrier_vertices, fixed_edges_vertices);
    fixed_edges_vertices =
        subtract(fixed_edges_vertices, fixed_edges_on_trim_surface_verts);
    if (corner_nodes)
      fixed_edges_vertices.merge(intersect(barrier_vertices, *corner_nodes));

    // get faces adjacent to vertices on fixed edges
    Range fixed_edges_faces;
    CHKERR moab.get_adjacencies(fixed_edges_vertices, 2, false,
                                fixed_edges_faces, moab::Interface::UNION);
    fixed_edges_faces = intersect(fixed_edges_faces, barrier_vertices_faces);

    if (debug && !fixed_edges_faces.empty())
      CHKERR SaveData(m_field.get_moab())("fixed_edges_faces.vtk",
                                          fixed_edges_faces);

    // get nodes on faces
    Range fixed_edges_faces_vertices;
    CHKERR moab.get_connectivity(fixed_edges_faces, fixed_edges_faces_vertices,
                                 false);
    barrier_vertices.merge(fixed_edges_faces_vertices);
  }

  auto remove_faces_on_skin = [&]() {
    MoFEMFunctionBegin;
    Range skin_faces = intersect(trimNewSurfaces, trim_tets_skin);
    if (!skin_faces.empty()) {
      Range add_to_barrier;
      CHKERR moab.get_connectivity(skin_faces, add_to_barrier, false);
      barrier_vertices.merge(add_to_barrier);
      for (auto f : skin_faces)
        trimNewSurfaces.erase(f);
    }
    MoFEMFunctionReturn(0);
  };

  auto get_trim_free_edges = [&]() {
    // get current surface skin
    Range trim_surf_skin;
    CHKERR skin.find_skin(0, trimNewSurfaces, false, trim_surf_skin);
    trim_surf_skin = subtract(trim_surf_skin, trim_tets_skin_edges);

    Range trim_surf_skin_verts;
    CHKERR moab.get_connectivity(trim_surf_skin, trim_surf_skin_verts, true);
    for (auto e : barrier_vertices)
      trim_surf_skin_verts.erase(e);

    Range free_edges;
    CHKERR moab.get_adjacencies(trim_surf_skin_verts, 1, false, free_edges,
                                moab::Interface::UNION);
    free_edges = intersect(free_edges, trim_surf_skin);

    return free_edges;
  };

  CHKERR remove_faces_on_skin();

  if (debug && !barrier_vertices.empty())
    CHKERR SaveData(m_field.get_moab())("barrier_vertices.vtk",
                                        barrier_vertices);

  int nn = 0;
  Range out_edges;
  while (!(out_edges = get_trim_free_edges()).empty()) {

    if (debug && !trimNewSurfaces.empty())
      CHKERR SaveData(m_field.get_moab())(
          "trimNewSurfaces_" + boost::lexical_cast<std::string>(nn) + ".vtk",
          trimNewSurfaces);

    if (debug && !out_edges.empty())
      CHKERR SaveData(m_field.get_moab())(
          "trimNewSurfacesOutsideVerts_" +
              boost::lexical_cast<std::string>(nn) + ".vtk",
          out_edges);

    Range outside_faces;
    CHKERR moab.get_adjacencies(out_edges, 2, false, outside_faces,
                                moab::Interface::UNION);
    trimNewSurfaces = subtract(trimNewSurfaces, outside_faces);
    ++nn;
  }

  if (debug && !trimNewSurfaces.empty())
    CHKERR SaveData(m_field.get_moab())(
        "trimNewSurfaces_" + boost::lexical_cast<std::string>(nn) + ".vtk",
        trimNewSurfaces);

  MoFEMFunctionReturn(0);
}

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
    CHKERR interface->findFacesWithThreeNodesOnInternalSurfaceSkin(
        meshset, split_bit, true, front_tris);
    CHKERR moab.delete_entities(&meshset, 1);
    ents = subtract(ents, front_tris);
  }
  // Remove entities on skin
  Range tets;
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      split_bit, BitRefLevel().set(), MBTET, tets);
  // Remove entities on skin
  Skinner skin(&moab);
  Range tets_skin;
  CHKERR skin.find_skin(0, tets, false, tets_skin);
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
  CHKERR m_field.getInterface(interface);
  EntityHandle meshset_volume;
  CHKERR moab.create_meshset(MESHSET_SET, meshset_volume);
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      split_bit, BitRefLevel().set(), MBTET, meshset_volume);
  EntityHandle meshset_surface;
  CHKERR moab.create_meshset(MESHSET_SET, meshset_surface);
  CHKERR moab.add_entities(meshset_surface, ents);
  CHKERR interface->getSides(meshset_surface, split_bit, true);
  CHKERR interface->splitSides(meshset_volume, bit, meshset_surface, true,
                               true);
  CHKERR moab.delete_entities(&meshset_volume, 1);
  CHKERR moab.delete_entities(&meshset_surface, 1);
  if (th) {
    Range prisms;
    ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit, BitRefLevel().set(), MBPRISM, prisms);
    for (Range::iterator pit = prisms.begin(); pit != prisms.end(); pit++) {
      int num_nodes;
      const EntityHandle *conn;
      CHKERR moab.get_connectivity(*pit, conn, num_nodes, true);
      MatrixDouble data(3, 3);
      CHKERR moab.tag_get_data(th, conn, 3, &data(0, 0));
      CHKERR moab.tag_set_data(th, &conn[3], 3, &data(0, 0));
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::mergeBadEdges(
    const int fraction_level, const Range &tets, const Range &surface,
    const Range &fixed_edges, const Range &corner_nodes, Range &edges_to_merge,
    Range &out_tets, Range &new_surf, Tag th, const bool update_meshsets,
    const BitRefLevel *bit_ptr, const bool debug) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  /**
   * \brief Merge nodes
   */
  struct MergeNodes {
    CoreInterface &mField;
    const bool onlyIfImproveQuality;
    Tag tH;
    bool updateMehsets;

    MergeNodes(CoreInterface &m_field, const bool only_if_improve_quality,
               Tag th, bool update_mehsets)
        : mField(m_field), onlyIfImproveQuality(only_if_improve_quality),
          tH(th), updateMehsets(update_mehsets) {
      mField.getInterface(nodeMergerPtr);
    }
    NodeMergerInterface *nodeMergerPtr;
    MoFEMErrorCode mergeNodes(int line_search, EntityHandle father,
                              EntityHandle mother, Range &proc_tets,
                              bool add_child = true) {
      moab::Interface &moab = mField.get_moab();
      MoFEMFunctionBegin;
      const EntityHandle conn[] = {father, mother};
      Range vert_tets;
      CHKERR moab.get_adjacencies(conn, 2, 3, false, vert_tets,
                                  moab::Interface::UNION);
      vert_tets = intersect(vert_tets, proc_tets);
      Range out_tets;
      CHKERR nodeMergerPtr->mergeNodes(father, mother, out_tets, &vert_tets,
                                       onlyIfImproveQuality, 0, line_search,
                                       tH);

      if (add_child && nodeMergerPtr->getSuccessMerge()) {

        Range::iterator lo, hi = proc_tets.begin();
        for (auto pt = vert_tets.pair_begin(); pt != vert_tets.pair_end();
             ++pt) {
          lo = proc_tets.lower_bound(hi, proc_tets.end(), pt->first);
          if (lo != proc_tets.end()) {
            hi = proc_tets.upper_bound(lo, proc_tets.end(), pt->second);
            proc_tets.erase(lo, hi);
          } else
            break;
        }
        proc_tets.merge(out_tets);

        auto &parent_child_map = nodeMergerPtr->getParentChildMap();

        struct ChangeChild {
          EntityHandle child;
          ChangeChild(const EntityHandle child) : child(child) {}
          void operator()(NodeMergerInterface::ParentChild &p) {
            p.cHild = child;
          }
        };

        std::vector<decltype(parentsChildMap.get<0>().begin())> it_vec;
        it_vec.reserve(parentsChildMap.size());

        for (auto &p : parent_child_map) {

          it_vec.clear();
          for (auto it = parentsChildMap.get<0>().equal_range(p.pArent);
               it.first != it.second; ++it.first)
            it_vec.emplace_back(it.first);

          for (auto it = parentsChildMap.get<1>().equal_range(p.pArent);
               it.first != it.second; ++it.first)
            it_vec.emplace_back(parentsChildMap.project<0>(it.first));

          for (auto &it : it_vec)
            parentsChildMap.modify(it, ChangeChild(p.cHild));
        }

        parentsChildMap.insert(parent_child_map.begin(),
                               parent_child_map.end());
      }
      MoFEMFunctionReturn(0);
    }

    MoFEMErrorCode updateRangeByChilds(Range &new_surf, Range &edges_to_merge,
                                       Range &not_merged_edges,
                                       bool add_child) {
      moab::Interface &moab = mField.get_moab();
      MoFEMFunctionBegin;
      if (add_child) {

        std::vector<EntityHandle> parents_ents_vec(parentsChildMap.size());
        for (auto &it : parentsChildMap)
          parents_ents_vec.emplace_back(it.pArent);
        Range parent_ents;
        parent_ents.insert_list(parents_ents_vec.begin(),
                                parents_ents_vec.end());

        Range surf_parent_ents = intersect(new_surf, parent_ents);
        new_surf = subtract(new_surf, surf_parent_ents);
        Range child_surf_ents;
        CHKERR updateRangeByChilds(parentsChildMap, surf_parent_ents,
                                   child_surf_ents);
        new_surf.merge(child_surf_ents);

        Range edges_to_merge_parent_ents =
            intersect(edges_to_merge, parent_ents);
        edges_to_merge = subtract(edges_to_merge, edges_to_merge_parent_ents);
        Range merged_child_edge_ents;
        CHKERR updateRangeByChilds(parentsChildMap, edges_to_merge_parent_ents,
                                   merged_child_edge_ents);

        Range not_merged_edges_child_ents =
            intersect(not_merged_edges, parent_ents);
        not_merged_edges =
            subtract(not_merged_edges, not_merged_edges_child_ents);
        Range not_merged_child_edge_ents;
        CHKERR updateRangeByChilds(parentsChildMap, not_merged_edges_child_ents,
                                   not_merged_child_edge_ents);

        merged_child_edge_ents =
            subtract(merged_child_edge_ents, not_merged_child_edge_ents);
        edges_to_merge.merge(merged_child_edge_ents);
        not_merged_edges.merge(not_merged_child_edge_ents);

        if (updateMehsets) {
          for (_IT_CUBITMESHSETS_FOR_LOOP_(
                   (*mField.getInterface<MeshsetsManager>()), cubit_it)) {
            EntityHandle cubit_meshset = cubit_it->meshset;
            Range meshset_ents;
            CHKERR moab.get_entities_by_handle(cubit_meshset, meshset_ents,
                                               true);
            Range child_ents;
            CHKERR updateRangeByChilds(parentsChildMap, meshset_ents,
                                       child_ents);
            CHKERR moab.add_entities(cubit_meshset, child_ents);
          }
        }
      }

      MoFEMFunctionReturn(0);
    };

  private:
    NodeMergerInterface::ParentChildMap parentsChildMap;
    std::vector<EntityHandle> childsVec;

    inline MoFEMErrorCode updateRangeByChilds(
        const NodeMergerInterface::ParentChildMap &parent_child_map,
        const Range &parents, Range &childs) {
      MoFEMFunctionBeginHot;
      childsVec.clear();
      childsVec.reserve(parents.size());
      for (auto pit = parents.pair_begin(); pit != parents.pair_end(); pit++) {
        auto it = parent_child_map.lower_bound(pit->first);
        if (it != parent_child_map.end()) {
          for (auto hi_it = parent_child_map.upper_bound(pit->second);
               it != hi_it; ++it)
            childsVec.emplace_back(it->cHild);
        }
      }
      childs.insert_list(childsVec.begin(), childsVec.end());
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
          maxLength(max_length) {
      gammaL = 1.;
      gammaQ = 3.;
      acceptedThrasholdMergeQuality = 0.5;
    }

    double
        gammaL; ///< Controls importance of length when ranking edges for merge
    double
        gammaQ; ///< Controls importance of quality when ranking edges for merge
    double acceptedThrasholdMergeQuality; ///< Do not merge quality if quality
                                          ///< above accepted thrashold

    MoFEMErrorCode operator()(const Range &tets, const Range &edges,
                              LengthMapData_multi_index &length_map,
                              double &ave) const {
      int num_nodes;
      const EntityHandle *conn;
      std::array<double, 12> coords;
      MoFEMFunctionBegin;
      VectorAdaptor s0(3, ublas::shallow_array_adaptor<double>(3, &coords[0]));
      VectorAdaptor s1(3, ublas::shallow_array_adaptor<double>(3, &coords[3]));
      VectorDouble3 delta(3);

      struct NodeQuality {
        EntityHandle node;
        double quality;
        NodeQuality(const EntityHandle node) : node(node), quality(1) {}
      };

      typedef multi_index_container<
          NodeQuality, indexed_by<

                           sequenced<>,

                           hashed_non_unique<tag<Ent_mi_tag>,
                                             member<NodeQuality, EntityHandle,
                                                    &NodeQuality::node>>

                           >>
          NodeQuality_sequence;

      NodeQuality_sequence node_quality_sequence;

      Range edges_nodes;
      CHKERR moab.get_connectivity(tets, edges_nodes, false);
      Range edges_tets;
      CHKERR moab.get_adjacencies(edges, 3, false, edges_tets,
                                  moab::Interface::UNION);
      edges_tets = intersect(edges_tets, tets);

      for (auto node : edges_nodes)
        node_quality_sequence.emplace_back(node);

      for (auto tet : edges_tets) {

        CHKERR moab.get_connectivity(tet, conn, num_nodes, true);
        if (tH)
          CHKERR moab.tag_get_data(tH, conn, num_nodes, coords.data());
        else
          CHKERR moab.get_coords(conn, num_nodes, coords.data());

        const double q = Tools::volumeLengthQuality(coords.data());

        for (auto n : {0, 1, 2, 3}) {
          auto it = node_quality_sequence.get<1>().find(conn[n]);
          if (it != node_quality_sequence.get<1>().end())
            const_cast<double &>(it->quality) = std::min(q, it->quality);
        }
      }

      for (auto edge : edges) {

        CHKERR moab.get_connectivity(edge, conn, num_nodes, true);

        if (tH)
          CHKERR moab.tag_get_data(tH, conn, num_nodes, coords.data());
        else
          CHKERR moab.get_coords(conn, num_nodes, coords.data());

        double q = 1;
        for (auto n : {0, 1}) {
          auto it = node_quality_sequence.get<1>().find(conn[n]);
          if (it != node_quality_sequence.get<1>().end())
            q = std::min(q, it->quality);
        }

        if (q < acceptedThrasholdMergeQuality) {
          noalias(delta) = (s0 - s1) / maxLength;
          double dot = inner_prod(delta, delta);
          double val = pow(q, gammaQ) * pow(dot, 0.5 * gammaL);
          length_map.insert(LengthMapData(val, q, edge));
        }
      }

      ave = 0;
      for (LengthMapData_multi_index::nth_index<0>::type::iterator mit =
               length_map.get<0>().begin();
           mit != length_map.get<0>().end(); mit++) {
        ave += mit->qUality;
      }
      ave /= length_map.size();
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
      FREE = 0,
      SKIN = 1 << 0,
      SURFACE = 1 << 1,
      SURFACE_SKIN = 1 << 2,
      FRONT_ENDS = 1 << 3,
      FIX_EDGES = 1 << 4,
      FIX_CORNERS = 1 << 5
    };

    typedef map<int, Range> SetsMap;

    MoFEMErrorCode classifyVerts(const Range &surface, const Range &tets,
                                 const Range &fixed_edges,
                                 const Range &corner_nodes,
                                 const Range &constrain_surface,
                                 SetsMap &sets_map) const {
      moab::Interface &moab(mField.get_moab());
      Skinner skin(&moab);
      MoFEMFunctionBegin;

      sets_map[FIX_CORNERS].merge(corner_nodes);
      Range fixed_verts;
      CHKERR moab.get_connectivity(fixed_edges, fixed_verts, true);
      sets_map[FIX_EDGES].swap(fixed_verts);

      Range tets_skin;
      CHKERR skin.find_skin(0, tets, false, tets_skin);
      Range tets_skin_edges;
      CHKERR moab.get_adjacencies(tets_skin, 1, false, tets_skin_edges,
                                  moab::Interface::UNION);
      tets_skin.merge(constrain_surface);

      // surface skin
      Range surface_skin;
      CHKERR skin.find_skin(0, surface, false, surface_skin);
      Range front_in_the_body;
      front_in_the_body = subtract(surface_skin, tets_skin_edges);
      Range front_ends;
      CHKERR skin.find_skin(0, front_in_the_body, false, front_ends);
      sets_map[FRONT_ENDS].swap(front_ends);

      Range surface_skin_verts;
      CHKERR moab.get_connectivity(surface_skin, surface_skin_verts, true);
      sets_map[SURFACE_SKIN].swap(surface_skin_verts);

      // surface
      Range surface_verts;
      CHKERR moab.get_connectivity(surface, surface_verts, true);
      sets_map[SURFACE].swap(surface_verts);

      // skin
      Range tets_skin_verts;
      CHKERR moab.get_connectivity(tets_skin, tets_skin_verts, true);
      sets_map[SKIN].swap(tets_skin_verts);

      Range tets_verts;
      CHKERR moab.get_connectivity(tets, tets_verts, true);
      sets_map[FREE].swap(tets_verts);

      MoFEMFunctionReturn(0);
    }

    MoFEMErrorCode getProcTets(const Range &tets, const Range &edges_to_merge,
                               Range &proc_tets) const {
      moab::Interface &moab(mField.get_moab());
      MoFEMFunctionBegin;
      Range edges_to_merge_verts;
      CHKERR moab.get_connectivity(edges_to_merge, edges_to_merge_verts, true);
      Range edges_to_merge_verts_tets;
      CHKERR moab.get_adjacencies(edges_to_merge_verts, 3, false,
                                  edges_to_merge_verts_tets,
                                  moab::Interface::UNION);
      edges_to_merge_verts_tets = intersect(edges_to_merge_verts_tets, tets);
      proc_tets.swap(edges_to_merge_verts_tets);
      MoFEMFunctionReturn(0);
    }

    MoFEMErrorCode removeBadEdges(const Range &surface, const Range &tets,
                                  const Range &fixed_edges,
                                  const Range &corner_nodes,
                                  Range &edges_to_merge,
                                  Range &not_merged_edges) {
      moab::Interface &moab(mField.get_moab());
      MoFEMFunctionBegin;

      // find skin
      Skinner skin(&moab);
      Range tets_skin;
      CHKERR skin.find_skin(0, tets, false, tets_skin);
      Range surface_skin;
      CHKERR skin.find_skin(0, surface, false, surface_skin);

      // end nodes
      Range tets_skin_edges;
      CHKERR moab.get_adjacencies(tets_skin, 1, false, tets_skin_edges,
                                  moab::Interface::UNION);
      Range surface_front;
      surface_front = subtract(surface_skin, tets_skin_edges);
      Range surface_front_nodes;
      CHKERR moab.get_connectivity(surface_front, surface_front_nodes, true);
      Range ends_nodes;
      CHKERR skin.find_skin(0, surface_front, false, ends_nodes);

      // remove bad merges

      // get surface and body skin verts
      Range surface_edges;
      CHKERR moab.get_adjacencies(surface, 1, false, surface_edges,
                                  moab::Interface::UNION);
      // get nodes on the surface
      Range surface_edges_verts;
      CHKERR moab.get_connectivity(surface_edges, surface_edges_verts, true);
      // get vertices on the body skin
      Range tets_skin_edges_verts;
      CHKERR moab.get_connectivity(tets_skin_edges, tets_skin_edges_verts,
                                   true);

      Range edges_to_remove;

      // remove edges self connected to body skin
      {
        Range ents_nodes_and_edges;
        ents_nodes_and_edges.merge(tets_skin_edges_verts);
        ents_nodes_and_edges.merge(tets_skin_edges);
        CHKERR removeSelfConectingEdges(ents_nodes_and_edges, edges_to_remove,
                                        0, false);
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
        CHKERR removeSelfConectingEdges(ents_nodes_and_edges, edges_to_remove,
                                        0, false);
      }
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      // remove edges adjacent corner_nodes execpt those on fixed edges
      Range fixed_edges_nodes;
      CHKERR moab.get_connectivity(fixed_edges, fixed_edges_nodes, true);
      {
        Range ents_nodes_and_edges;
        ents_nodes_and_edges.merge(fixed_edges_nodes);
        ents_nodes_and_edges.merge(ends_nodes);
        ents_nodes_and_edges.merge(corner_nodes);
        ents_nodes_and_edges.merge(fixed_edges);
        CHKERR removeSelfConectingEdges(ents_nodes_and_edges, edges_to_remove,
                                        0, false);
      }
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      // remove edges self connected to surface
      CHKERR removeSelfConectingEdges(surface_edges, edges_to_remove, 0, false);
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      // remove edges self contented on surface skin
      {
        Range ents_nodes_and_edges;
        ents_nodes_and_edges.merge(surface_skin);
        ents_nodes_and_edges.merge(fixed_edges_nodes);
        CHKERR removeSelfConectingEdges(ents_nodes_and_edges, edges_to_remove,
                                        0, false);
      }
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      // remove edges connecting crack front and fixed edges, except those short
      {
        Range ents_nodes_and_edges;
        ents_nodes_and_edges.merge(surface_skin.subset_by_type(MBEDGE));
        ents_nodes_and_edges.merge(fixed_edges.subset_by_type(MBEDGE));
        CHKERR removeSelfConectingEdges(ents_nodes_and_edges, edges_to_remove,
                                        0, false);
      }
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      // remove crack front nodes connected to the surface, except those short
      {
        Range ents_nodes_and_edges;
        ents_nodes_and_edges.merge(surface_front_nodes);
        ents_nodes_and_edges.merge(surface_front);
        ents_nodes_and_edges.merge(tets_skin_edges_verts);
        ents_nodes_and_edges.merge(tets_skin_edges);
        CHKERR removeSelfConectingEdges(ents_nodes_and_edges, edges_to_remove,
                                        tOL, false);
      }
      edges_to_merge = subtract(edges_to_merge, edges_to_remove);
      not_merged_edges.merge(edges_to_remove);

      MoFEMFunctionReturn(0);
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
      if (length > 0) {
        Range::iterator eit = ents_nodes_edges.begin();
        for (; eit != ents_nodes_edges.end();) {

          int num_nodes;
          const EntityHandle *conn;
          CHKERR moab.get_connectivity(*eit, conn, num_nodes, true);
          double coords[6];
          if (tH)
            CHKERR moab.tag_get_data(tH, conn, num_nodes, coords);
          else
            CHKERR moab.get_coords(conn, num_nodes, coords);

          auto get_edge_length = [coords]() {
            FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3> t_coords(
                &coords[0], &coords[1], &coords[2]);
            FTensor::Tensor1<double, 3> t_delta;
            FTensor::Index<'i', 3> i;
            t_delta(i) = t_coords(i);
            ++t_coords;
            t_delta(i) -= t_coords(i);
            return sqrt(t_delta(i) * t_delta(i));
          };

          if (get_edge_length() < tOL) {
            eit = ents_nodes_edges.erase(eit);
          } else {
            eit++;
          }
        }
      }
      edges_to_remove.swap(ents_nodes_edges);
      if (debug) {
        CHKERR SaveData(moab)("edges_to_remove.vtk", edges_to_remove);
      }
      MoFEMFunctionReturn(0);
    }
  };

  Range not_merged_edges;
  const double tol = 1e-3;
  CHKERR Toplogy(m_field, th, tol * aveLength)
      .removeBadEdges(surface, tets, fixed_edges, corner_nodes, edges_to_merge,
                      not_merged_edges);
  Toplogy::SetsMap sets_map;
  CHKERR Toplogy(m_field, th, tol * aveLength)
      .classifyVerts(surface, tets, fixed_edges, corner_nodes, constrainSurface,
                     sets_map);
  if (debug) {
    for (Toplogy::SetsMap::reverse_iterator sit = sets_map.rbegin();
         sit != sets_map.rend(); sit++) {
      std::string name = "classification_verts_" +
                         boost::lexical_cast<std::string>(sit->first) + ".vtk";
      if (!sit->second.empty())
        CHKERR SaveData(moab)(name, sit->second);
    }
  }
  Range proc_tets;
  CHKERR Toplogy(m_field, th, tol * aveLength)
      .getProcTets(tets, edges_to_merge, proc_tets);
  out_tets = subtract(tets, proc_tets);

  if (bit_ptr) {
    Range all_out_ents = out_tets;
    for (int dd = 2; dd >= 0; dd--) {
      CHKERR moab.get_adjacencies(out_tets, dd, false, all_out_ents,
                                  moab::Interface::UNION);
    }
    CHKERR m_field.getInterface<BitRefManager>()->addBitRefLevel(all_out_ents,
                                                                 *bit_ptr);
  }

  int nb_nodes_merged = 0;
  LengthMapData_multi_index length_map;
  new_surf = surface;

  auto save_merge_step = [&](const int pp, const Range collapsed_edges) {
    MoFEMFunctionBegin;
    Range adj_faces;
    CHKERR moab.get_adjacencies(proc_tets, 2, false, adj_faces,
                                moab::Interface::UNION);
    std::string name;
    name = "node_merge_step_" + boost::lexical_cast<std::string>(pp) + ".vtk";
    CHKERR SaveData(moab)(
        name, unite(intersect(new_surf, adj_faces), collapsed_edges));
    name =
        "edges_to_merge_step_" + boost::lexical_cast<std::string>(pp) + ".vtk";
    CHKERR SaveData(moab)(
        name, unite(intersect(new_surf, adj_faces), edges_to_merge));
    MoFEMFunctionReturn(0);
  };

  if (debug)
    CHKERR save_merge_step(0, Range());

  double ave0 = 0, ave = 0, min = 0, min_p = 0, min_pp;
  for (int pp = 0; pp != nbMaxMergingCycles; ++pp) {

    int nb_nodes_merged_p = nb_nodes_merged;
    length_map.clear();
    min_pp = min_p;
    min_p = min;
    CHKERR LengthMap(m_field, th, aveLength)(proc_tets, edges_to_merge,
                                             length_map, ave);
    min = length_map.get<2>().begin()->qUality;
    if (pp == 0) {
      ave0 = ave;
    }

    int nn = 0;
    Range collapsed_edges;
    MergeNodes merge_nodes(m_field, true, th, update_meshsets);

    for (auto mit = length_map.get<0>().begin();
         mit != length_map.get<0>().end(); mit++, nn++) {

      if (!mit->skip) {

        auto get_father_and_mother =
            [&](EntityHandle &father, EntityHandle &mother, int &line_search) {
              MoFEMFunctionBegin;
              int num_nodes;
              const EntityHandle *conn;
              CHKERR moab.get_connectivity(mit->eDge, conn, num_nodes, true);
              std::array<int, 2> conn_type = {0, 0};
              for (int nn = 0; nn != 2; nn++) {
                for (Toplogy::SetsMap::reverse_iterator sit = sets_map.rbegin();
                     sit != sets_map.rend(); sit++) {
                  if (sit->second.find(conn[nn]) != sit->second.end()) {
                    conn_type[nn] |= sit->first;
                  }
                }
              }
              int type_father, type_mother;
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
              if (type_father == type_mother) {
                line_search = lineSearchSteps;
              }
              MoFEMFunctionReturn(0);
            };

        int line_search = 0;
        EntityHandle father, mother;
        CHKERR get_father_and_mother(father, mother, line_search);
        CHKERR merge_nodes.mergeNodes(line_search, father, mother, proc_tets);
        if (m_field.getInterface<NodeMergerInterface>()->getSuccessMerge()) {
          const EntityHandle father_and_mother[] = {father, mother};
          Range adj_tets;
          CHKERR moab.get_adjacencies(father_and_mother, 1, 3, false, adj_tets);
          Range adj_tets_nodes;
          CHKERR moab.get_connectivity(adj_tets, adj_tets_nodes, true);
          Range adj_edges;
          CHKERR moab.get_adjacencies(adj_tets_nodes, 1, false, adj_edges,
                                      moab::Interface::UNION);
          for (auto ait : adj_edges) {
            auto miit = length_map.get<1>().find(ait);
            if (miit != length_map.get<1>().end())
              (const_cast<LengthMapData &>(*miit)).skip = true;
          }
          nb_nodes_merged++;
          collapsed_edges.insert(mit->eDge);
        }

        if (nn > static_cast<int>(length_map.size() / fraction_level))
          break;
        if (mit->qUality > ave)
          break;
      }
    }

    CHKERR merge_nodes.updateRangeByChilds(new_surf, edges_to_merge,
                                           not_merged_edges, true);

    PetscPrintf(m_field.get_comm(),
                "(%d) Number of nodes merged %d ave q %3.4e min q %3.4e\n", pp,
                nb_nodes_merged, ave, min);

    if (debug)
      CHKERR save_merge_step(pp + 1, collapsed_edges);

    if (nb_nodes_merged == nb_nodes_merged_p)
      break;
    if (min > 1e-2 && min == min_pp)
      break;
    if (min > ave0)
      break;

    Range adj_edges;
    CHKERR moab.get_adjacencies(proc_tets, 1, false, adj_edges,
                                moab::Interface::UNION);
    edges_to_merge = intersect(edges_to_merge, adj_edges);
    CHKERR Toplogy(m_field, th, tol * aveLength)
        .removeBadEdges(new_surf, proc_tets, fixed_edges, corner_nodes,
                        edges_to_merge, not_merged_edges);
  }

  if (bit_ptr)
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevel(proc_tets,
                                                                 *bit_ptr);

  out_tets.merge(proc_tets);
  Range adj_faces;
  CHKERR moab.get_adjacencies(out_tets, 2, false, adj_faces,
                              moab::Interface::UNION);
  new_surf = intersect(new_surf, adj_faces);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::mergeBadEdges(
    const int fraction_level, const BitRefLevel cut_bit,
    const BitRefLevel trim_bit, const BitRefLevel bit, const Range &surface,
    const Range &fixed_edges, const Range &corner_nodes, Tag th,
    const bool update_meshsets, const bool debug) {
  CoreInterface &m_field = cOre;
  MoFEMFunctionBegin;
  Range tets_level;
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      trim_bit, BitRefLevel().set(), MBTET, tets_level);

  Range edges_to_merge;
  CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByRefLevel(
      trim_bit, cut_bit | trim_bit, edges_to_merge);

  // get all entities not in database
  Range all_ents_not_in_database_before;
  CHKERR cOre.getInterface<BitRefManager>()->getAllEntitiesNotInDatabase(
      all_ents_not_in_database_before);

  edges_to_merge = edges_to_merge.subset_by_type(MBEDGE);
  if (debug)
    CHKERR SaveData(m_field.get_moab())("edges_to_merge.vtk", edges_to_merge);

  Range out_new_tets, out_new_surf;
  CHKERR mergeBadEdges(fraction_level, tets_level, surface, fixed_edges,
                       corner_nodes, edges_to_merge, out_new_tets, out_new_surf,
                       th, update_meshsets, &bit, debug);

  // get all entities not in database after merge
  Range all_ents_not_in_database_after;
  CHKERR cOre.getInterface<BitRefManager>()->getAllEntitiesNotInDatabase(
      all_ents_not_in_database_after);

  // delete hanging entities
  all_ents_not_in_database_after =
      subtract(all_ents_not_in_database_after, all_ents_not_in_database_before);
  Range meshsets;
  CHKERR m_field.get_moab().get_entities_by_type(0, MBENTITYSET, meshsets,
                                                 false);
  for (auto m : meshsets)
    CHKERR m_field.get_moab().remove_entities(m,
                                              all_ents_not_in_database_after);

  m_field.get_moab().delete_entities(all_ents_not_in_database_after);

  mergedVolumes.swap(out_new_tets);
  mergedSurfaces.swap(out_new_surf);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::setTagData(Tag th, const BitRefLevel bit) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  Range nodes;
  if (bit.none())
    CHKERR moab.get_entities_by_type(0, MBVERTEX, nodes);
  else
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit, BitRefLevel().set(), MBVERTEX, nodes);
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
  if (bit.none())
    CHKERR moab.get_entities_by_type(0, MBVERTEX, nodes);
  else
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit, mask, MBVERTEX, nodes);
  std::vector<double> coords(3 * nodes.size());
  CHKERR moab.tag_get_data(th, nodes, &coords[0]);
  CHKERR moab.set_coords(nodes, &coords[0]);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::saveCutEdges(const std::string prefix) {
  CoreInterface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  CHKERR SaveData(moab)(prefix + "out_vol.vtk", vOlume);
  CHKERR SaveData(moab)(prefix + "out_surface.vtk", sUrface);
  CHKERR SaveData(moab)(prefix + "out_cut_edges.vtk", cutEdges);
  CHKERR SaveData(moab)(prefix + "out_cut_new_volumes.vtk", cutNewVolumes);
  CHKERR SaveData(moab)(prefix + "out_cut_new_surfaces.vtk", cutNewSurfaces);
  CHKERR SaveData(moab)(prefix + "out_cut_zero_distance_ents.vtk",
                        zeroDistanceEnts);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CutMeshInterface::saveTrimEdges() {
  moab::Interface &moab = cOre.getInterface<CoreInterface>()->get_moab();
  MoFEMFunctionBegin;
  CHKERR SaveData(moab)("out_trim_surface.vtk", sUrface);
  CHKERR SaveData(moab)("out_trim_new_volumes.vtk", trimNewVolumes);
  CHKERR SaveData(moab)("out_trim_new_surfaces.vtk", trimNewSurfaces);
  CHKERR SaveData(moab)("out_trim_edges.vtk", trimEdges);
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM