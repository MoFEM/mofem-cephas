/**
 * @file PostProc.cpp
 * @brief Post processing elements and operators
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <MoFEM.hpp>

namespace MoFEM {

PostProcGenerateRefMeshBase::PostProcGenerateRefMeshBase()
    : hoNodes(PETSC_TRUE), defMaxLevel(0), optPrefix(""), countEle(0),
      countVertEle(0), nbVertices(0), nbEles(0) {}

MoFEMErrorCode PostProcGenerateRefMeshBase::getOptions() {
  MoFEMFunctionBegin;

  CHKERR PetscOptionsGetInt(optPrefix.c_str(), "-max_post_proc_ref_level",
                            &defMaxLevel, PETSC_NULL);
  CHKERR PetscOptionsGetBool(optPrefix.c_str(), "-max_post_ho_nodes", &hoNodes,
                             PETSC_NULL);

  if (defMaxLevel < 0)
    SETERRQ(PETSC_COMM_WORLD, MOFEM_INVALID_DATA,
            "Wrong parameter -max_post_proc_ref_level "
            "should be positive number");

  MoFEMFunctionReturn(0);
};

MoFEMErrorCode PostProcGenerateRefMesh<MBTET>::generateReferenceElementMesh() {
  MoFEMFunctionBegin;

  CHKERR getOptions();

  const int max_level = defMaxLevel;

  moab::Core core_ref;
  moab::Interface &moab_ref = core_ref;

  auto create_reference_element = [&moab_ref]() {
    MoFEMFunctionBegin;
    constexpr double base_coords[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    EntityHandle nodes[4];
    for (int nn = 0; nn < 4; nn++) {
      CHKERR
      moab_ref.create_vertex(&base_coords[3 * nn], nodes[nn]);
    }
    EntityHandle tet;
    CHKERR moab_ref.create_element(MBTET, nodes, 4, tet);
    MoFEMFunctionReturn(0);
  };

  MoFEM::CoreTmp<-1> m_core_ref(moab_ref, PETSC_COMM_SELF, -2);
  MoFEM::Interface &m_field_ref = m_core_ref;

  auto refine_ref_tetrahedron = [this, &m_field_ref, max_level]() {
    MoFEMFunctionBegin;
    // seed ref mofem database by setting bit ref level to reference
    // tetrahedron
    CHKERR
    m_field_ref.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, BitRefLevel().set(0));
    for (int ll = 0; ll != max_level; ++ll) {
      MOFEM_TAG_AND_LOG_C("WORLD", Sev::noisy, "PostProc", "Refine Level %d",
                          ll);
      Range edges;
      CHKERR m_field_ref.getInterface<BitRefManager>()
          ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(ll),
                                         BitRefLevel().set(), MBEDGE, edges);
      Range tets;
      CHKERR m_field_ref.getInterface<BitRefManager>()
          ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(ll),
                                         BitRefLevel(ll).set(), MBTET, tets);
      // refine mesh
      MeshRefinement *m_ref;
      CHKERR m_field_ref.getInterface(m_ref);
      CHKERR m_ref->addVerticesInTheMiddleOfEdges(edges,
                                                  BitRefLevel().set(ll + 1));
      CHKERR m_ref->refineTets(tets, BitRefLevel().set(ll + 1));
    }
    MoFEMFunctionReturn(0);
  };

  auto get_ref_gauss_pts_and_shape_functions = [this, max_level, &moab_ref,
                                                &m_field_ref]() {
    MoFEMFunctionBegin;
    for (int ll = 0; ll != max_level + 1; ++ll) {
      Range tets;
      CHKERR
      m_field_ref.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
          BitRefLevel().set(ll), BitRefLevel().set(ll), MBTET, tets);
      if (hoNodes) {
        EntityHandle meshset;
        CHKERR moab_ref.create_meshset(MESHSET_SET, meshset);
        CHKERR moab_ref.add_entities(meshset, tets);
        CHKERR moab_ref.convert_entities(meshset, true, false, false);
        CHKERR moab_ref.delete_entities(&meshset, 1);
      }
      Range elem_nodes;
      CHKERR moab_ref.get_connectivity(tets, elem_nodes, false);

      auto &gauss_pts = levelGaussPtsOnRefMesh[ll];
      gauss_pts.resize(elem_nodes.size(), 4, false);
      std::map<EntityHandle, int> little_map;
      Range::iterator nit = elem_nodes.begin();
      for (int gg = 0; nit != elem_nodes.end(); nit++, gg++) {
        CHKERR moab_ref.get_coords(&*nit, 1, &gauss_pts(gg, 0));
        little_map[*nit] = gg;
      }
      gauss_pts = trans(gauss_pts);

      auto &ref_tets = levelRef[ll];
      Range::iterator tit = tets.begin();
      for (int tt = 0; tit != tets.end(); ++tit, ++tt) {
        const EntityHandle *conn;
        int num_nodes;
        CHKERR moab_ref.get_connectivity(*tit, conn, num_nodes, false);
        if (tt == 0) {
          ref_tets.resize(tets.size(), num_nodes);
        }
        for (int nn = 0; nn != num_nodes; ++nn) {
          ref_tets(tt, nn) = little_map[conn[nn]];
        }
      }

      auto &shape_functions = levelShapeFunctions[ll];
      shape_functions.resize(elem_nodes.size(), 4);
      CHKERR ShapeMBTET(&*shape_functions.data().begin(), &gauss_pts(0, 0),
                        &gauss_pts(1, 0), &gauss_pts(2, 0), elem_nodes.size());
    }
    MoFEMFunctionReturn(0);
  };

  levelRef.resize(max_level + 1);
  levelGaussPtsOnRefMesh.resize(max_level + 1);
  levelShapeFunctions.resize(max_level + 1);

  CHKERR create_reference_element();
  CHKERR refine_ref_tetrahedron();
  CHKERR get_ref_gauss_pts_and_shape_functions();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PostProcGenerateRefMesh<MBHEX>::generateReferenceElementMesh() {
  MoFEMFunctionBegin;

  CHKERR getOptions();
#ifndef NDEBUG
  if (defMaxLevel > 0)
    MOFEM_LOG("WORLD", Sev::warning)
        << "Refinement for hexes is not implemented";
#endif

  moab::Core core_ref;
  moab::Interface &moab_ref = core_ref;

  auto create_reference_element = [&moab_ref]() {
    MoFEMFunctionBegin;
    constexpr double base_coords[] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                      0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    EntityHandle nodes[8];
    for (int nn = 0; nn < 8; nn++)
      CHKERR moab_ref.create_vertex(&base_coords[3 * nn], nodes[nn]);
    EntityHandle hex;
    CHKERR moab_ref.create_element(MBHEX, nodes, 8, hex);
    MoFEMFunctionReturn(0);
  };

  auto add_ho_nodes = [&]() {
    MoFEMFunctionBegin;
    Range hexes;
    CHKERR moab_ref.get_entities_by_type(0, MBHEX, hexes, true);
    EntityHandle meshset;
    CHKERR moab_ref.create_meshset(MESHSET_SET, meshset);
    CHKERR moab_ref.add_entities(meshset, hexes);
    CHKERR moab_ref.convert_entities(meshset, true, true, true);
    CHKERR moab_ref.delete_entities(&meshset, 1);
    MoFEMFunctionReturn(0);
  };

  auto set_gauss_pts = [&](std::map<EntityHandle, int> &little_map) {
    MoFEMFunctionBegin;
    Range hexes;
    CHKERR moab_ref.get_entities_by_type(0, MBHEX, hexes, true);
    Range hexes_nodes;
    CHKERR moab_ref.get_connectivity(hexes, hexes_nodes, false);
    auto &gauss_pts = levelGaussPtsOnRefMesh[0];
    gauss_pts.resize(hexes_nodes.size(), 4, false);
    size_t gg = 0;
    for (auto node : hexes_nodes) {
      CHKERR moab_ref.get_coords(&node, 1, &gauss_pts(gg, 0));
      little_map[node] = gg;
      ++gg;
    }
    gauss_pts = trans(gauss_pts);
    MoFEMFunctionReturn(0);
  };

  auto set_ref_hexes = [&](std::map<EntityHandle, int> &little_map) {
    MoFEMFunctionBegin;
    Range hexes;
    CHKERR moab_ref.get_entities_by_type(0, MBHEX, hexes, true);
    size_t hh = 0;
    auto &ref_hexes = levelRef[0];
    for (auto hex : hexes) {
      const EntityHandle *conn;
      int num_nodes;
      CHKERR moab_ref.get_connectivity(hex, conn, num_nodes, false);
      if (ref_hexes.size2() != num_nodes) {
        ref_hexes.resize(hexes.size(), num_nodes);
      }
      for (int nn = 0; nn != num_nodes; ++nn) {
        ref_hexes(hh, nn) = little_map[conn[nn]];
      }
      ++hh;
    }
    MoFEMFunctionReturn(0);
  };

  auto set_shape_functions = [&]() {
    MoFEMFunctionBegin;
    auto &gauss_pts = levelGaussPtsOnRefMesh[0];
    auto &shape_functions = levelShapeFunctions[0];
    const auto nb_gauss_pts = gauss_pts.size2();
    shape_functions.resize(nb_gauss_pts, 8);
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      const double ksi = gauss_pts(0, gg);
      const double zeta = gauss_pts(1, gg);
      const double eta = gauss_pts(2, gg);
      shape_functions(gg, 0) = N_MBHEX0(ksi, zeta, eta);
      shape_functions(gg, 1) = N_MBHEX1(ksi, zeta, eta);
      shape_functions(gg, 2) = N_MBHEX2(ksi, zeta, eta);
      shape_functions(gg, 3) = N_MBHEX3(ksi, zeta, eta);
      shape_functions(gg, 4) = N_MBHEX4(ksi, zeta, eta);
      shape_functions(gg, 5) = N_MBHEX5(ksi, zeta, eta);
      shape_functions(gg, 6) = N_MBHEX6(ksi, zeta, eta);
      shape_functions(gg, 7) = N_MBHEX7(ksi, zeta, eta);
    }
    MoFEMFunctionReturn(0);
  };

  levelRef.resize(1);
  levelGaussPtsOnRefMesh.resize(1);
  levelShapeFunctions.resize(1);

  CHKERR create_reference_element();
  if (hoNodes)
    CHKERR add_ho_nodes();
  std::map<EntityHandle, int> little_map;
  CHKERR set_gauss_pts(little_map);
  CHKERR set_ref_hexes(little_map);
  CHKERR set_shape_functions();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PostProcGenerateRefMesh<MBTRI>::generateReferenceElementMesh() {
  MoFEMFunctionBegin;

  CHKERR getOptions();
  const int max_level = defMaxLevel;

  moab::Core core_ref;
  moab::Interface &moab_ref = core_ref;

  auto create_reference_element = [&moab_ref]() {
    MoFEMFunctionBegin;
    constexpr double base_coords[] = {

        0, 0,
        0,

        1, 0,
        0,

        0, 1,
        0

    };
    EntityHandle nodes[3];
    for (int nn = 0; nn != 3; ++nn)
      CHKERR moab_ref.create_vertex(&base_coords[3 * nn], nodes[nn]);
    EntityHandle tri;
    CHKERR moab_ref.create_element(MBTRI, nodes, 3, tri);

    Range edges;
    CHKERR moab_ref.get_adjacencies(&tri, 1, 1, true, edges,
                                    moab::Interface::UNION);

    MoFEMFunctionReturn(0);
  };

  CHKERR create_reference_element();

  MoFEM::CoreTmp<-1> m_core_ref(moab_ref, PETSC_COMM_SELF, -2);
  MoFEM::Interface &m_field_ref = m_core_ref;

  auto refine_ref_triangles = [this, &m_field_ref, max_level]() {
    MoFEMFunctionBegin;
    // seed ref mofem database by setting bit ref level to reference
    // tetrahedron
    CHKERR
    m_field_ref.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 2, BitRefLevel().set(0));

    for (int ll = 0; ll != max_level; ++ll) {
      MOFEM_TAG_AND_LOG_C("WORLD", Sev::noisy, "PostProc", "Refine Level %d",
                          ll);
      Range edges;
      CHKERR m_field_ref.getInterface<BitRefManager>()
          ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(ll),
                                         BitRefLevel().set(), MBEDGE, edges);
      Range tris;
      CHKERR m_field_ref.getInterface<BitRefManager>()
          ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(ll),
                                         BitRefLevel(ll).set(), MBTRI, tris);
      // refine mesh
      auto m_ref = m_field_ref.getInterface<MeshRefinement>();
      CHKERR m_ref->addVerticesInTheMiddleOfEdges(edges,
                                                  BitRefLevel().set(ll + 1));
      CHKERR m_ref->refineTris(tris, BitRefLevel().set(ll + 1));
    }
    MoFEMFunctionReturn(0);
  };

  auto set_gauss_pts = [&](std::map<EntityHandle, int> &little_map) {
    MoFEMFunctionBegin;
    Range faces;
    CHKERR moab_ref.get_entities_by_type(0, MBTRI, faces, true);
    Range faces_nodes;
    CHKERR moab_ref.get_connectivity(faces, faces_nodes, false);
    auto &gauss_pts = levelGaussPtsOnRefMesh[0];
    gauss_pts.resize(faces_nodes.size(), 4, false);
    size_t gg = 0;
    for (auto node : faces_nodes) {
      CHKERR moab_ref.get_coords(&node, 1, &gauss_pts(gg, 0));
      little_map[node] = gg;
      ++gg;
    }
    gauss_pts = trans(gauss_pts);
    MoFEMFunctionReturn(0);
  };

  auto get_ref_gauss_pts_and_shape_functions = [this, max_level, &moab_ref,
                                                &m_field_ref]() {
    MoFEMFunctionBegin;
    for (int ll = 0; ll != max_level + 1; ++ll) {

      Range tris;
      CHKERR
      m_field_ref.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
          BitRefLevel().set(ll), BitRefLevel().set(ll), MBTRI, tris);

      if (hoNodes) {
        EntityHandle meshset;
        CHKERR moab_ref.create_meshset(MESHSET_SET, meshset);
        CHKERR moab_ref.add_entities(meshset, tris);
        CHKERR moab_ref.convert_entities(meshset, true, false, false);
        CHKERR moab_ref.delete_entities(&meshset, 1);
      }

      Range elem_nodes;
      CHKERR moab_ref.get_connectivity(tris, elem_nodes, false);

      auto &gauss_pts = levelGaussPtsOnRefMesh[ll];
      gauss_pts.resize(elem_nodes.size(), 3, false);
      std::map<EntityHandle, int> little_map;
      Range::iterator nit = elem_nodes.begin();
      for (int gg = 0; nit != elem_nodes.end(); nit++, gg++) {
        CHKERR moab_ref.get_coords(&*nit, 1, &gauss_pts(gg, 0));
        little_map[*nit] = gg;
      }
      gauss_pts = trans(gauss_pts);

      auto &ref_tris = levelRef[ll];
      Range::iterator tit = tris.begin();
      for (int tt = 0; tit != tris.end(); ++tit, ++tt) {
        const EntityHandle *conn;
        int num_nodes;
        CHKERR moab_ref.get_connectivity(*tit, conn, num_nodes, false);
        if (tt == 0) {
          ref_tris.resize(tris.size(), num_nodes);
        }
        for (int nn = 0; nn != num_nodes; ++nn) {
          ref_tris(tt, nn) = little_map[conn[nn]];
        }
      }

      auto &shape_functions = levelShapeFunctions[ll];
      shape_functions.resize(elem_nodes.size(), 3);
      CHKERR ShapeMBTRI(&*shape_functions.data().begin(), &gauss_pts(0, 0),
                        &gauss_pts(1, 0), elem_nodes.size());
    }
    MoFEMFunctionReturn(0);
  };

  levelRef.resize(max_level + 1);
  levelGaussPtsOnRefMesh.resize(max_level + 1);
  levelShapeFunctions.resize(max_level + 1);

  CHKERR refine_ref_triangles();
  CHKERR get_ref_gauss_pts_and_shape_functions();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PostProcGenerateRefMesh<MBQUAD>::generateReferenceElementMesh() {
  MoFEMFunctionBegin;

  CHKERR getOptions();
#ifndef NDEBUG
  if (defMaxLevel > 0)
    MOFEM_LOG("WORLD", Sev::warning)
        << "Refinement for quad is not implemented";
#endif

  moab::Core core_ref;
  moab::Interface &moab_ref = core_ref;

  auto create_reference_element = [&moab_ref]() {
    MoFEMFunctionBegin;
    constexpr double base_coords[] = {

        0, 0,
        0,

        1, 0,
        0,

        1, 1,
        0,

        0, 1,
        0

    };
    EntityHandle nodes[4];
    for (int nn = 0; nn < 4; nn++)
      CHKERR moab_ref.create_vertex(&base_coords[3 * nn], nodes[nn]);
    EntityHandle quad;
    CHKERR moab_ref.create_element(MBQUAD, nodes, 4, quad);
    MoFEMFunctionReturn(0);
  };

  auto add_ho_nodes = [&]() {
    MoFEMFunctionBegin;
    Range quads;
    CHKERR moab_ref.get_entities_by_type(0, MBQUAD, quads, true);
    EntityHandle meshset;
    CHKERR moab_ref.create_meshset(MESHSET_SET, meshset);
    CHKERR moab_ref.add_entities(meshset, quads);
    CHKERR moab_ref.convert_entities(meshset, true, true, true);
    CHKERR moab_ref.delete_entities(&meshset, 1);
    MoFEMFunctionReturn(0);
  };

  auto set_gauss_pts = [&](std::map<EntityHandle, int> &little_map) {
    MoFEMFunctionBegin;
    Range quads;
    CHKERR moab_ref.get_entities_by_type(0, MBQUAD, quads, true);
    Range quads_nodes;
    CHKERR moab_ref.get_connectivity(quads, quads_nodes, false);
    auto &gauss_pts = levelGaussPtsOnRefMesh[0];
    gauss_pts.resize(quads_nodes.size(), 4, false);
    size_t gg = 0;
    for (auto node : quads_nodes) {
      CHKERR moab_ref.get_coords(&node, 1, &gauss_pts(gg, 0));
      little_map[node] = gg;
      ++gg;
    }
    gauss_pts = trans(gauss_pts);
    MoFEMFunctionReturn(0);
  };

  auto set_ref_quads = [&](std::map<EntityHandle, int> &little_map) {
    MoFEMFunctionBegin;
    Range quads;
    CHKERR moab_ref.get_entities_by_type(0, MBQUAD, quads, true);
    size_t hh = 0;
    auto &ref_quads = levelRef[0];
    for (auto quad : quads) {
      const EntityHandle *conn;
      int num_nodes;
      CHKERR moab_ref.get_connectivity(quad, conn, num_nodes, false);
      if (ref_quads.size2() != num_nodes) {
        ref_quads.resize(quads.size(), num_nodes);
      }
      for (int nn = 0; nn != num_nodes; ++nn) {
        ref_quads(hh, nn) = little_map[conn[nn]];
      }
      ++hh;
    }
    MoFEMFunctionReturn(0);
  };

  auto set_shape_functions = [&]() {
    MoFEMFunctionBegin;
    auto &gauss_pts = levelGaussPtsOnRefMesh[0];
    auto &shape_functions = levelShapeFunctions[0];
    const auto nb_gauss_pts = gauss_pts.size2();
    shape_functions.resize(nb_gauss_pts, 4);
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      const double ksi = gauss_pts(0, gg);
      const double zeta = gauss_pts(1, gg);
      shape_functions(gg, 0) = N_MBQUAD0(ksi, zeta);
      shape_functions(gg, 1) = N_MBQUAD1(ksi, zeta);
      shape_functions(gg, 2) = N_MBQUAD2(ksi, zeta);
      shape_functions(gg, 3) = N_MBQUAD3(ksi, zeta);
    }
    MoFEMFunctionReturn(0);
  };

  levelRef.resize(1);
  levelGaussPtsOnRefMesh.resize(1);
  levelShapeFunctions.resize(1);

  CHKERR create_reference_element();
  if (hoNodes)
    CHKERR add_ho_nodes();
  std::map<EntityHandle, int> little_map;
  CHKERR set_gauss_pts(little_map);
  CHKERR set_ref_quads(little_map);
  CHKERR set_shape_functions();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PostProcGenerateRefMesh<MBEDGE>::generateReferenceElementMesh() {
  MoFEMFunctionBegin;

  CHKERR getOptions();
#ifndef NDEBUG
  if (defMaxLevel > 0)
    MOFEM_LOG("WORLD", Sev::warning)
        << "Refinement for edges is not implemented";
#endif

  auto set_gauss_pts = [&](std::map<EntityHandle, int> &little_map) {
    MoFEMFunctionBegin;

    int nb_nodes = 2;
    if (hoNodes)
      nb_nodes = 3;

    auto &gauss_pts = levelGaussPtsOnRefMesh[0];
    gauss_pts.resize(2, nb_nodes, false);
    gauss_pts.clear();

    int nn = 0;
    for (; nn != 2; ++nn) {
      gauss_pts(0, nn) = static_cast<double>(nn);
      little_map[nn] = nn;
    }

    if (nn < nb_nodes) {
      gauss_pts(0, nn) = 0.5;
      little_map[nn] = 2;
    }

    MoFEMFunctionReturn(0);
  };

  auto set_ref_edges = [&](std::map<EntityHandle, int> &little_map) {
    MoFEMFunctionBegin;

    int level = 0;
    int nb_edges = level + 1;

    int nb_nodes = 2;
    if (hoNodes)
      nb_nodes = 3;

    auto &ref_edges = levelRef[level];
    ref_edges.resize(nb_edges, nb_nodes, false);

    for (int ee = 0; ee != nb_edges; ++ee) {
      int nn = 0;
      for (; nn != 2; ++nn) {
        ref_edges(ee, nn) = nb_nodes * ee + nn;
      }
      if (nn < nb_nodes) {
        ref_edges(ee, nn) = nb_nodes * ee + 2;
      }
    }

    MoFEMFunctionReturn(0);
  };

  auto set_shape_functions = [&]() {
    MoFEMFunctionBegin;
    auto &gauss_pts = levelGaussPtsOnRefMesh[0];
    auto &shape_functions = levelShapeFunctions[0];
    const auto nb_gauss_pts = gauss_pts.size2();
    shape_functions.resize(nb_gauss_pts, 2);
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      const double ksi = gauss_pts(0, gg);
      shape_functions(gg, 0) = N_MBEDGE0(ksi);
      shape_functions(gg, 1) = N_MBEDGE1(ksi);
    }
    MoFEMFunctionReturn(0);
  };

  levelRef.resize(1);
  levelGaussPtsOnRefMesh.resize(1);
  levelShapeFunctions.resize(1);

  std::map<EntityHandle, int> little_map;
  CHKERR set_gauss_pts(little_map);
  CHKERR set_ref_edges(little_map);
  CHKERR set_shape_functions();

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM