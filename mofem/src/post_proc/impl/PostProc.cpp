/**
 * @file PostProc.cpp
 * @brief Post processing elements and operators
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <MoFEM.hpp>

namespace MoFEM {

struct PostProcGenerateRefMeshBase {

  std::vector<MatrixDouble> levelShapeFunctions;
  std::vector<MatrixDouble> levelGaussPtsOnRefMesh;
  std::vector<ublas::matrix<int>> levelRef;

  EntityHandle startingVertEleHandle;
  std::vector<double *> verticesOnEleArrays;
  EntityHandle startingEleHandle;
  EntityHandle *eleConn;

  int countEle;
  int countVertEle;

  int nbVertices;
  int nbEles;

  PostProcGenerateRefMeshBase();
  MoFEMErrorCode getOptions();

  virtual MoFEMErrorCode generateReferenceElementMesh() = 0;

  PetscBool hoNodes;
  int defMaxLevel;
  std::string optPrefix;
};

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

template <>
struct PostProcGenerateRefMesh<MBTET> : public PostProcGenerateRefMeshBase {

  using PostProcGenerateRefMeshBase::PostProcGenerateRefMeshBase;

  MoFEMErrorCode generateReferenceElementMesh() {
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
        MOFEM_TAG_AND_LOG_C("WORLD", Sev::verbose, "PostProc",
                            "Refine Level %d", ll);
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
                          &gauss_pts(1, 0), &gauss_pts(2, 0),
                          elem_nodes.size());
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
};

template <>
struct PostProcGenerateRefMesh<MBHEX> : public PostProcGenerateRefMeshBase {
  MoFEMErrorCode generateReferenceElementMesh() {
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
};

template <>
struct PostProcGenerateRefMesh<MBTRI> : public PostProcGenerateRefMeshBase {

  using PostProcGenerateRefMeshBase::PostProcGenerateRefMeshBase;

  MoFEMErrorCode generateReferenceElementMesh() {
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

    CHKERR moab_ref.write_file("aaa.vtk", "VTK", "");



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
        MOFEM_TAG_AND_LOG_C("WORLD", Sev::verbose, "PostProc",
                            "Refine Level %d", ll);
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
};

template <>
struct PostProcGenerateRefMesh<MBQUAD> : public PostProcGenerateRefMeshBase {
  MoFEMErrorCode generateReferenceElementMesh() {
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
};

template <>
struct PostProcGenerateRefMesh<MBEDGE> : public PostProcGenerateRefMeshBase {
  MoFEMErrorCode generateReferenceElementMesh() {
    MoFEMFunctionBegin;

    CHKERR getOptions();
#ifndef NDEBUG
    if (defMaxLevel > 0)
      MOFEM_LOG("WORLD", Sev::warning)
          << "Refinement for hexes is not implemented";
#endif

    auto set_gauss_pts = [&](std::map<EntityHandle, int> &little_map) {
      MoFEMFunctionBegin;

      int nb_nodes = 2;
      if(hoNodes)
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
      if(hoNodes)
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
};

template <typename E> int PostProcBrokenMeshInMoabBase<E>::getMaxLevel() const {
  auto get_element_max_dofs_order = [&]() {
    int max_order = 0;
    auto dofs_vec = E::getDataVectorDofsPtr();
    for (auto &dof : *dofs_vec) {
      const int dof_order = dof->getDofOrder();
      max_order = (max_order < dof_order) ? dof_order : max_order;
    };
    return max_order;
  };
  const auto dof_max_order = get_element_max_dofs_order();
  return (dof_max_order > 0) ? (dof_max_order - 1) / 2 : 0;
};

template <typename E>
PostProcBrokenMeshInMoabBase<E>::PostProcBrokenMeshInMoabBase(
    MoFEM::Interface &m_field)
    : E(m_field), postProcMesh(coreMesh) {}

template <typename E>
PostProcBrokenMeshInMoabBase<E>::~PostProcBrokenMeshInMoabBase() {
  ParallelComm *pcomm_post_proc_mesh =
      ParallelComm::get_pcomm(&postProcMesh, MYPCOMM_INDEX);
  if (pcomm_post_proc_mesh != NULL)
    delete pcomm_post_proc_mesh;
}

template <typename E> int PostProcBrokenMeshInMoabBase<E>::getRule(int order) {
  return -1;
};

template <typename E>
MoFEMErrorCode PostProcBrokenMeshInMoabBase<E>::setGaussPts(int order) {
  MoFEMFunctionBegin;

  auto type = type_from_handle(this->getFEEntityHandle());

  PostProcGenerateRefMeshPtr ref_ele;

  try {
    ref_ele = refElementsMap.at(type);
  } catch (const out_of_range &e) {
    SETERRQ1(
        PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
        "Generation of reference elements for type <%s> is not implemented",
        moab::CN::EntityTypeName(type));
  }

  auto set_gauss_pts = [&](auto &level_gauss_pts_on_ref_mesh, auto &level_ref,
                           auto &level_shape_functions,

                           auto start_vert_handle, auto start_ele_handle,
                           auto &verts_array, auto &conn, auto &ver_count,
                           auto &ele_count

                       ) {
    MoFEMFunctionBegin;

    size_t level = getMaxLevel();
    level = std::min(level, level_gauss_pts_on_ref_mesh.size() - 1);

    auto &level_ref_gauss_pts = level_gauss_pts_on_ref_mesh[level];
    auto &level_ref_ele = level_ref[level];
    auto &shape_functions = level_shape_functions[level];
    E::gaussPts.resize(level_ref_gauss_pts.size1(), level_ref_gauss_pts.size2(),
                       false);
    noalias(E::gaussPts) = level_ref_gauss_pts;

    const auto fe_ent = E::numeredEntFiniteElementPtr->getEnt();
    auto get_fe_coords = [&]() {
      const EntityHandle *conn;
      int num_nodes;
      CHK_MOAB_THROW(
          E::mField.get_moab().get_connectivity(fe_ent, conn, num_nodes, true),
          "error get connectivity");
      VectorDouble coords(num_nodes * 3);
      CHK_MOAB_THROW(
          E::mField.get_moab().get_coords(conn, num_nodes, &*coords.begin()),
          "error get coordinates");
      return coords;
    };

    auto coords = get_fe_coords();

    const int num_nodes = level_ref_gauss_pts.size2();
    mapGaussPts.resize(level_ref_gauss_pts.size2());

    FTensor::Index<'i', 3> i;
    FTensor::Tensor0<FTensor::PackPtr<double *, 1>> t_n(
        &*shape_functions.data().begin());
    FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3> t_coords(
        &verts_array[0][ver_count], &verts_array[1][ver_count],
        &verts_array[2][ver_count]);
    for (int gg = 0; gg != num_nodes; ++gg, ++ver_count) {

      mapGaussPts[gg] = start_vert_handle + ver_count;

      auto set_float_precision = [](const double x) {
        if (std::abs(x) < std::numeric_limits<float>::epsilon())
          return 0.;
        else
          return x;
      };

      t_coords(i) = 0;
      auto t_ele_coords = getFTensor1FromArray<3, 3>(coords);
      for (int nn = 0; nn != CN::VerticesPerEntity(type); ++nn) {
        t_coords(i) += t_n * t_ele_coords(i);
        ++t_ele_coords;
        ++t_n;
      }

      for (auto ii : {0, 1, 2})
        t_coords(ii) = set_float_precision(t_coords(ii));

      ++t_coords;
    }

    Tag th;
    int def_in_the_loop = -1;
    CHKERR postProcMesh.tag_get_handle("NB_IN_THE_LOOP", 1, MB_TYPE_INTEGER, th,
                                       MB_TAG_CREAT | MB_TAG_SPARSE,
                                       &def_in_the_loop);

    postProcElements.clear();
    const int num_el = level_ref_ele.size1();
    const int num_nodes_on_ele = level_ref_ele.size2();
    auto start_e = start_ele_handle + ele_count;
    postProcElements = Range(start_e, start_e + num_el - 1);
    for (auto tt = 0; tt != level_ref_ele.size1(); ++tt, ++ele_count) {
      for (int nn = 0; nn != num_nodes_on_ele; ++nn) {
        conn[num_nodes_on_ele * ele_count + nn] =
            mapGaussPts[level_ref_ele(tt, nn)];
      }
    }

    const int n_in_the_loop = E::nInTheLoop;
    CHKERR postProcMesh.tag_clear_data(th, postProcElements, &n_in_the_loop);

    MoFEMFunctionReturn(0);
  };

  CHKERR set_gauss_pts(

      ref_ele->levelGaussPtsOnRefMesh, ref_ele->levelRef,
      ref_ele->levelShapeFunctions,

      ref_ele->startingVertEleHandle, ref_ele->startingEleHandle,
      ref_ele->verticesOnEleArrays, ref_ele->eleConn, ref_ele->countVertEle,
      ref_ele->countEle

  );

  MoFEMFunctionReturn(0);
};

template <typename E>
MoFEMErrorCode PostProcBrokenMeshInMoabBase<E>::preProcess() {
  moab::Interface &moab = coreMesh;
  MoFEMFunctionBegin;

  ParallelComm *pcomm_post_proc_mesh =
      ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
  if (pcomm_post_proc_mesh != NULL)
    delete pcomm_post_proc_mesh;

  CHKERR postProcMesh.delete_mesh();

  auto get_ref_ele = [&](const EntityType type) {
    PostProcGenerateRefMeshPtr ref_ele_ptr;

    try {

      ref_ele_ptr = refElementsMap.at(type);

    } catch (const out_of_range &e) {

      switch (type) {
      case MBTET:
        ref_ele_ptr = boost::make_shared<PostProcGenerateRefMesh<MBTET>>();
        break;
      case MBHEX:
        ref_ele_ptr = boost::make_shared<PostProcGenerateRefMesh<MBHEX>>();
        break;
      case MBTRI:
        ref_ele_ptr = boost::make_shared<PostProcGenerateRefMesh<MBTRI>>();
        break;
      case MBQUAD:
        ref_ele_ptr = boost::make_shared<PostProcGenerateRefMesh<MBQUAD>>();
        break;
      case MBEDGE:
        ref_ele_ptr = boost::make_shared<PostProcGenerateRefMesh<MBEDGE>>();
        break;
      default:
        MOFEM_LOG("SELF", Sev::error)
            << "Generation of reference elements for type < "
            << moab::CN::EntityTypeName(type) << " > is not implemented";
        CHK_THROW_MESSAGE(MOFEM_NOT_IMPLEMENTED, "Element not implemented");
      }
    }

    CHK_THROW_MESSAGE(ref_ele_ptr->generateReferenceElementMesh(),
                      "Error when generating reference element");

    refElementsMap[type] = ref_ele_ptr;

    return ref_ele_ptr;
  };

  auto fe_ptr = this->problemPtr->numeredFiniteElementsPtr;

  auto miit =
      fe_ptr->template get<Composite_Name_And_Part_mi_tag>().lower_bound(
          boost::make_tuple(this->getFEName(), this->getLoFERank()));
  auto hi_miit =
      fe_ptr->template get<Composite_Name_And_Part_mi_tag>().upper_bound(
          boost::make_tuple(this->getFEName(), this->getHiFERank()));

  const int number_of_ents_in_the_loop = this->getLoopSize();
  if (std::distance(miit, hi_miit) != number_of_ents_in_the_loop) {
    SETERRQ(E::mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "Wrong size of indicies. Inconsistent size number of iterated "
            "elements iterated by problem and from range.");
  }

  for (auto &m : refElementsMap) {
    m.second->nbVertices = 0;
    m.second->nbEles = 0;
    m.second->countEle = 0;
    m.second->countVertEle = 0;
  }

  for (; miit != hi_miit; ++miit) {
    auto type = (*miit)->getEntType();
    auto ref_ele = get_ref_ele(type);

    // Set pointer to element. So that getDataVectorDofsPtr in getMaxLevel
    // can work
    E::numeredEntFiniteElementPtr = *miit;
    bool add = true;
    if (E::exeTestHook) {
      add = E::exeTestHook(this);
    }

    if (add) {
      size_t level = getMaxLevel();
      level = std::min(level, ref_ele->levelGaussPtsOnRefMesh.size() - 1);
      ref_ele->nbVertices += ref_ele->levelGaussPtsOnRefMesh[level].size2();
      ref_ele->nbEles += ref_ele->levelRef[level].size1();
    }
  }

  auto alloc_vertices_and_elements_on_post_proc_mesh = [&]() {
    MoFEMFunctionBegin;

    ReadUtilIface *iface;
    CHKERR postProcMesh.query_interface(iface);

    for (auto &m : refElementsMap) {
      if (m.second->nbEles) {
        CHKERR iface->get_node_coords(3, m.second->nbVertices, 0,
                                      m.second->startingVertEleHandle,
                                      m.second->verticesOnEleArrays);
        CHKERR iface->get_element_connect(
            m.second->nbEles, m.second->levelRef[0].size2(), m.first, 0,
            m.second->startingEleHandle, m.second->eleConn);

        m.second->countEle = 0;
        m.second->countVertEle = 0;
      }
    }

    MoFEMFunctionReturn(0);
  };

  CHKERR alloc_vertices_and_elements_on_post_proc_mesh();

  MoFEMFunctionReturn(0);
}

template <typename E>
MoFEMErrorCode PostProcBrokenMeshInMoabBase<E>::postProcess() {
  MoFEMFunctionBeginHot;

  auto update_elements = [&]() {
    ReadUtilIface *iface;
    CHKERR this->postProcMesh.query_interface(iface);
    MoFEMFunctionBegin;

    for (auto &m : refElementsMap) {
      if (m.second->nbEles) {
        MOFEM_TAG_AND_LOG("SELF", Sev::noisy, "PostProc")
            << "Update < " << moab::CN::EntityTypeName(m.first)
            << m.second->countEle;
        CHKERR iface->update_adjacencies(
            m.second->startingEleHandle, m.second->countEle,
            m.second->levelRef[0].size2(), m.second->eleConn);
      }
    }

    MoFEMFunctionReturn(0);
  };

  auto resolve_shared_ents = [&]() {
    MoFEMFunctionBegin;
    auto get_lower_dimension = [&]() {
      int min_dim = 3;
      for (auto &m : refElementsMap) {
        const int dim = moab::CN::Dimension(m.first);
        if (m.second->nbEles) {
          min_dim = std::min(dim, min_dim);
        }
      }
      return min_dim;
    };

    auto remove_obsolete_entities = [&](auto &&min_dim) {
      Range edges;
      CHKERR postProcMesh.get_entities_by_type(0, MBEDGE, edges, false);

      Range faces;
      CHKERR postProcMesh.get_entities_by_dimension(0, 2, faces, false);

      Range vols;
      CHKERR postProcMesh.get_entities_by_dimension(0, 3, vols, false);

      Range ents;

      if (min_dim > 1)
        CHKERR postProcMesh.delete_entities(edges);
      else
        ents.merge(edges);

      if (min_dim > 2)
        CHKERR postProcMesh.delete_entities(faces);
      else
        ents.merge(faces);

      ents.merge(vols);

      return ents;
    };

    auto ents = remove_obsolete_entities(get_lower_dimension());

    ParallelComm *pcomm_post_proc_mesh =
        ParallelComm::get_pcomm(&(postProcMesh), MYPCOMM_INDEX);
    if (pcomm_post_proc_mesh == NULL) {
      // wrapRefMeshComm =
      // boost::make_shared<WrapMPIComm>(T::mField.get_comm(), false);
      pcomm_post_proc_mesh = new ParallelComm(
          &(postProcMesh),
          PETSC_COMM_WORLD /*(T::wrapRefMeshComm)->get_comm()*/);
    }

    int rank = E::mField.get_comm_rank();
    CHKERR postProcMesh.tag_clear_data(pcomm_post_proc_mesh->part_tag(), ents,
                                       &rank);

    CHKERR pcomm_post_proc_mesh->resolve_shared_ents(0);

    MoFEMFunctionReturn(0);
  };

  CHKERR update_elements();
  CHKERR resolve_shared_ents();

  MoFEMFunctionReturnHot(0);
}

template <>
boost::shared_ptr<PostProcBrokenMeshInMoab<VolumeElementForcesAndSourcesCore>>
make_post_proc_fe_in_moab(MoFEM::Interface &m_field) {
  return boost::make_shared<
      PostProcBrokenMeshInMoab<VolumeElementForcesAndSourcesCore>>(m_field);
}

template <>
boost::shared_ptr<PostProcBrokenMeshInMoab<FaceElementForcesAndSourcesCore>>
make_post_proc_fe_in_moab(MoFEM::Interface &m_field) {
  return boost::make_shared<
      PostProcBrokenMeshInMoab<FaceElementForcesAndSourcesCore>>(m_field);
}

template <>
boost::shared_ptr<PostProcBrokenMeshInMoab<EdgeElementForcesAndSourcesCore>>
make_post_proc_fe_in_moab(MoFEM::Interface &m_field) {
  return boost::make_shared<
      PostProcBrokenMeshInMoab<EdgeElementForcesAndSourcesCore>>(m_field);
}

} // namespace MoFEM