/**
 * \file PoissonCutFEM.hpp
 * \example PoissonCutFEM.hpp
 *
 * Example of implementation for cutfem
 */

extern "C" {
#include <phg-quadrule/quad.h>
}

// Define name if it has not been defined yet
#ifndef __POISSONCUTFEM_HPP__
#define __POISSONCUTFEM_HPP__

// Namespace that contains necessary UDOs, will be included in the main program
namespace Poisson3DCutFEMOperators {

// Declare FTensor index for 2D problem
FTensor::Index<'i', SPACE_DIM> i;

struct HexMeshRefinement : public MeshRefinement {
  // reusing constructor from MeshRefinement
  using MeshRefinement::MeshRefinement;

  MoFEMErrorCode refineHexToTets(const Range &_hexes, const BitRefLevel &bit,
                                 SetEdgeBitsFun set_edge_bits, int verb,
                                 const bool debug) {
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    auto refined_ents_ptr = m_field.get_ref_ents();
    ReadUtilIface *read_util;
    MoFEMFunctionBegin;

    CHKERR m_field.get_moab().query_interface(read_util);

    Range ents_to_set_bit;

    auto get_parent_ents_view = [&](const auto type_child,
                                    const auto type_parent) {
      auto &ref_ents =
          refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>();
      auto range =
          ref_ents.equal_range(boost::make_tuple(type_child, type_parent));

      RefEntity_multiIndex_view_by_ordered_parent_entity ref_parent_ents_view;

      using I = decltype(range.first);

      boost::function<bool(I it)> tester = [&](I it) {
        return ((*it)->getBitRefLevel() & bit).any();
      };

      boost::function<MoFEMErrorCode(I f, I s)> inserter = [&](I f, I s) {
        ref_parent_ents_view.insert(f, s);
        return 0;
      };

      rangeInserter(range.first, range.second, tester, inserter);

      return ref_parent_ents_view;
    };

    auto ref_parent_ents_view = get_parent_ents_view(MBVERTEX, MBEDGE);

    // Only work on Hexahedrons
    auto hexes = _hexes.subset_by_type(MBHEX);

    std::vector<EntityHandle> parent_hexes_refined;
    std::vector<EntityHandle> vertices_of_created_tets;
    std::vector<BitRefEdges> parent_edges_bit_vec;
    std::vector<int> nb_new_tets_vec;
    std::vector<int> sub_type_vec;

    parent_hexes_refined.reserve(hexes.size());
    vertices_of_created_tets.reserve(
        6 * hexes.size()); // Maximum 6 tets per hex, FIXME:
    parent_edges_bit_vec.reserve(hexes.size());
    nb_new_tets_vec.reserve(hexes.size());
    sub_type_vec.reserve(hexes.size());

    // Loop over all hexes which are going to be refined
    for (auto hit : hexes) {
      // Get hexahedron connectivity
      const EntityHandle *conn;
      int num_nodes;
      CHKERR moab.get_connectivity(hit, conn, num_nodes, true);
      // std::cout << " The num nodes are.. " << num_nodes << "\n"; 
      // Initialize variables for hexahedron refinement
      BitRefEdges parent_edges_bit(0);
      // EntityHandle edge_new_nodes[12] = {0}; // 12 edges in a hexahedron

      if (true) {
        // Swap nodes forward
        EntityHandle _conn_[8];
        std::copy(&conn[0], &conn[8], &_conn_[0]);

        // Build connectivity for refined tetrahedrons
        EntityHandle new_tets_conns[6 * 4]; // Maximum 6 tets per hex FIXME: 28 for structured
        std::fill(&new_tets_conns[0], &new_tets_conns[6 * 4], no_handle);

        int sub_type = -1, nb_new_tets = 0;
        int is_unstructured_ref = 1;
        switch (is_unstructured_ref) {
        // Define cases for different numbers of split edges in a hexahedron
        case 1:
          sub_type = 0;
          hex_to_tet_type_1(_conn_, new_tets_conns);
          nb_new_tets = 6;
          break;
        // case 2: //FIXME:
        //   sub_type = 1;
        //   hex_to_tet_type_2(_conn_, new_tets_conns);
        //   nb_new_tets = 28;
        // break;

        default:
          SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "Impossible case");
        }
        // std::cout << " We are at a problem at line: " << __LINE__ << "\n"; 

        // If this element was not refined or was refined with different
        // patterns of split edges
        parent_hexes_refined.emplace_back(hit);
        for (int ht = 0; ht != nb_new_tets; ++ht) {
          for (auto nn : {0, 1, 2, 3}) // 4 vertices per tet
            vertices_of_created_tets.emplace_back(new_tets_conns[4 * ht + nn]);
        }
        parent_edges_bit_vec.emplace_back(parent_edges_bit);
        nb_new_tets_vec.emplace_back(nb_new_tets);
        sub_type_vec.emplace_back(sub_type);
        
        // std::cout << " We are at a problem at line: " << __LINE__ << "\n"; 
      } else {
        if (debug) {
          RefEntity_multiIndex::iterator hit_miit;
          hit_miit = refined_ents_ptr->find(hit);
          if (hit_miit == refined_ents_ptr->end())
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "Cannot find this hex");
        }
        ents_to_set_bit.insert(hit);
      }
    }
    //  std::cout << " We are at a problem at line: " << __LINE__ << "\n"; 
    // Create tetrahedrons
    EntityHandle start_e = 0;
    EntityHandle *conn_moab;
    const int nb_new_tets = vertices_of_created_tets.size() / 4;
    read_util->get_element_connect(nb_new_tets, 4, MBTET, 0, start_e,
                                   conn_moab);
    std::copy(vertices_of_created_tets.begin(), vertices_of_created_tets.end(),
              conn_moab);
    CHKERR read_util->update_adjacencies(start_e, nb_new_tets, 4, conn_moab);
    ents_to_set_bit.insert(start_e, start_e + nb_new_tets - 1);

    // Create adj entities
    Range ents;
    for (auto d : {1, 2}) {
      CHKERR moab.get_adjacencies(ents_to_set_bit, d, true, ents,
                                  moab::Interface::UNION);
    }
    // std::cout << " We are at a problem at line: " << __LINE__ << "\n"; 
    SetParent set_parent;

    // Set parents and adjacencies for the new tets
    for (int idx = 0; idx != parent_hexes_refined.size(); ++idx) {
      const EntityHandle hit = parent_hexes_refined[idx];

      const BitRefEdges &parent_edges_bit = parent_edges_bit_vec[idx];
      const int nb_new_tets = nb_new_tets_vec[idx];

      std::array<EntityHandle, 8> ref_tets;
      for (int tt = 0; tt != nb_new_tets; ++tt, ++start_e)
        ref_tets[tt] = start_e;

      if (nb_new_tets) {

        CHKERR moab.tag_clear_data(cOre.get_th_RefBitEdge(), &ref_tets[0],
                                   nb_new_tets, &parent_edges_bit);
        CHKERR moab.tag_clear_data(cOre.get_th_RefParentHandle(), &ref_tets[0],
                                   nb_new_tets, &hit);

        // hash map of nodes (RefEntity) by edges (EntityHandle)
        std::map<EntityHandle /*edge*/, EntityHandle /*node*/>
            map_ref_nodes_by_edges;
        // std::cout << " We are at a problem at line: " << __LINE__ << "\n"; 
        Range tit_edges;
        CHKERR moab.get_adjacencies(&hit, 1, 1, false, tit_edges);
        for (auto edge : tit_edges) {
          auto eit = ref_parent_ents_view.find(edge);
          if (eit != ref_parent_ents_view.end()) {
            map_ref_nodes_by_edges[(*eit)->getParentEnt()] =
                eit->get()->getEnt();
          }
        }

        // std::cout << " We are at a problem at line: " << __LINE__ << "\n"; 
        // std::cout << " We ar tit_edges.size() " << tit_edges.size() << "\n"; 
        // find parents for new edges and faces
        // get tet edges and faces
        Range tit_faces;
        CHKERR moab.get_adjacencies(&hit, 1, 2, false, tit_faces);
        // std::cout << " We ar tit_faces.size() " << tit_faces.size() << "\n"; 
        Range edges_nodes[12], faces_nodes[6];
        // for edges - add ref nodes
        // edges_nodes[ee] - contains all nodes on edge ee including mid nodes
        // if exist
        // std::cout << " We are at a problem at line: " << __LINE__ << "\n"; 
        Range::iterator eit = tit_edges.begin();
        for (int ee = 0; eit != tit_edges.end(); eit++, ee++) {
          CHKERR moab.get_connectivity(&*eit, 1, edges_nodes[ee], true);
          auto map_miit = map_ref_nodes_by_edges.find(*eit);
          if (map_miit != map_ref_nodes_by_edges.end()) {
            edges_nodes[ee].insert(map_miit->second);
          }
        }
        // std::cout << " We are at a problem at line: " << __LINE__ << "\n"; 
        // for faces - add ref nodes
        // faces_nodes[ff] - contains all nodes on face ff including mid nodes
        // if exist
        Range::iterator fit = tit_faces.begin();
        for (int ff = 0; fit != tit_faces.end(); fit++, ff++) {
          CHKERR moab.get_connectivity(&*fit, 1, faces_nodes[ff], true);
          // Get edges on face and loop over those edges to add mid-nodes to
          // range
          Range fit_edges;
          CHKERR moab.get_adjacencies(&*fit, 1, 1, false, fit_edges);
          for (Range::iterator eit2 = fit_edges.begin();
               eit2 != fit_edges.end(); eit2++) {
            auto map_miit = map_ref_nodes_by_edges.find(*eit2);
            if (map_miit != map_ref_nodes_by_edges.end()) {
              faces_nodes[ff].insert(map_miit->second);
            }
          }
        }
        // add ref nodes to tet
        // tet_nodes contains all nodes on tet including mid edge nodes
        Range tet_nodes;
        CHKERR moab.get_connectivity(&hit, 1, tet_nodes, true);
        for (std::map<EntityHandle, EntityHandle>::iterator map_miit =
                 map_ref_nodes_by_edges.begin();
             map_miit != map_ref_nodes_by_edges.end(); map_miit++) {
          tet_nodes.insert(map_miit->second);
        }
        Range ref_edges;
        // Get all all edges of refined tets
        CHKERR moab.get_adjacencies(ref_tets.data(), nb_new_tets, 1, false,
                                    ref_edges, moab::Interface::UNION);
        // Check for all ref edge and set parents
        for (Range::iterator reit = ref_edges.begin(); reit != ref_edges.end();
             reit++) {
          Range ref_edges_nodes;
          CHKERR moab.get_connectivity(&*reit, 1, ref_edges_nodes, true);

          // Check if ref edge is an coarse edge (loop over coarse tet edges)
          int ee = 0;
          for (; ee < 6; ee++) {
            // Two nodes are common (node[0],node[1],ref_node (if exist))
            // this tests if given edge is contained by edge of refined
            // tetrahedral
            if (intersect(edges_nodes[ee], ref_edges_nodes).size() == 2) {
              EntityHandle edge = tit_edges[ee];
              CHKERR set_parent(*reit, edge, refined_ents_ptr, cOre);
              break;
            }
          }
          if (ee < 6)
            continue; // this refined edge is contained by edge of tetrahedral
          // check if ref edge is in coarse face
          int ff = 0;
          for (; ff < 4; ff++) {
            // two nodes are common (node[0],node[1],ref_node (if exist))
            // this tests if given face is contained by face of  tetrahedral
            if (intersect(faces_nodes[ff], ref_edges_nodes).size() == 2) {
              EntityHandle face = tit_faces[ff];
              CHKERR set_parent(*reit, face, refined_ents_ptr, cOre);
              break;
            }
          }
          if (ff < 4)
            continue; // this refined edge is contained by face of tetrahedral

          // check if ref edge is in coarse tetrahedral (i.e. that is internal
          // edge of refined tetrahedral)
          if (intersect(tet_nodes, ref_edges_nodes).size() == 2) {
            CHKERR set_parent(*reit, hit, refined_ents_ptr, cOre);
            continue;
          }

          // Refined edge is not child of any edge, face or tetrahedral, this is
          // Impossible edge
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "impossible refined edge");
        }

        Range ref_faces;
        CHKERR moab.get_adjacencies(ref_tets.data(), nb_new_tets, 2, false,
                                    ref_faces, moab::Interface::UNION);
        Tag th_interface_side;
        const int def_side[] = {0};
        CHKERR moab.tag_get_handle("INTERFACE_SIDE", 1, MB_TYPE_INTEGER,
                                   th_interface_side,
                                   MB_TAG_CREAT | MB_TAG_SPARSE, def_side);
        // Check for all ref faces
        for (auto rfit : ref_faces) {
          Range ref_faces_nodes;
          CHKERR moab.get_connectivity(&rfit, 1, ref_faces_nodes, true);
          // Check if ref face is in coarse face
          int ff = 0;
          for (; ff < 4; ff++) {
            // Check if refined triangle is contained by face of tetrahedral
            if (intersect(faces_nodes[ff], ref_faces_nodes).size() == 3) {
              EntityHandle face = tit_faces[ff];
              CHKERR set_parent(rfit, face, refined_ents_ptr, cOre);
              int side = 0;
              // Set face side if it is on interface
              CHKERR moab.tag_get_data(th_interface_side, &face, 1, &side);
              CHKERR moab.tag_set_data(th_interface_side, &rfit, 1, &side);
              break;
            }
          }
          if (ff < 4)
            continue; // this face is contained by one of tetrahedrons
          // check if ref face is in coarse tetrahedral
          // this is ref face which is contained by tetrahedral volume
          if (intersect(tet_nodes, ref_faces_nodes).size() == 3) {
            CHKERR set_parent(rfit, hit, refined_ents_ptr, cOre);
            continue;
          }
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "impossible refined face");
        }
      }
    }
    // std::cout << " We are at a problem at line: " << __LINE__ << "\n"; 
    CHKERR set_parent(refined_ents_ptr);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevel(
        ents_to_set_bit, bit, true, verb, &ents);
    // std::cout << " We are at a problem at line: " << __LINE__ << "\n"; 
    MoFEMFunctionReturn(0);
  }

  MoFEMErrorCode refineHexToTets(const Range &_hexes, const BitRefLevel &bit,
                                 int verb = QUIET, const bool debug = false) {

    MoFEMFunctionBegin;

    auto set_edge_bits = [](moab::Interface &moab,
                            RefEntity_multiIndex_view_by_ordered_parent_entity
                                &ref_parent_ents_view,
                            EntityHandle hex, BitRefEdges &parent_edges_bit,
                            EntityHandle *edge_new_nodes, int *split_edges) {
      MoFEMFunctionBeginHot;

      for (int ee = 0; ee != 8; ++ee) {
        EntityHandle edge;
        CHKERR moab.side_element(hex, 1, ee, edge);
        auto eit = ref_parent_ents_view.find(edge);
        if (eit != ref_parent_ents_view.end()) {
          edge_new_nodes[ee] = (*eit)->getEnt();
          split_edges[parent_edges_bit.count()] = ee;
          parent_edges_bit.set(ee, 1);
        }
      }

      MoFEMFunctionReturnHot(0);
    };

    CHKERR refineHexToTets(_hexes, bit, set_edge_bits, verb, debug);

    MoFEMFunctionReturn(0);
  }
};

struct SetIntegrationOnActiveCells : public VolumeElementForcesAndSourcesCore {

  SetIntegrationOnActiveCells(MoFEM::Interface &m_field, const Range &active_cells, const Range &inside_cells)
       : VolumeElementForcesAndSourcesCore(m_field), activeCells(active_cells), insideCells(inside_cells) {}

  MoFEMErrorCode operator()(ForcesAndSourcesCore *fe_raw_ptr, int order_row,
                            int order_col, int order_data) {
    MoFEMFunctionBegin;

    constexpr bool debug = false;

    constexpr int numNodes = 4;
    constexpr int numEdges = 6;
    constexpr int refinementLevels = 3;

    auto fe_ptr = static_cast<VolumeElementForcesAndSourcesCore *>(fe_raw_ptr);
    auto &m_field = fe_raw_ptr->mField;
    auto &moab = fe_raw_ptr->mField.get_moab();
    auto fe_handle = fe_raw_ptr->getFEEntityHandle();
    auto bit_mng = m_field.getInterface<BitRefManager>();

    const int rule = 2 * order_data + 1; // FIXME:
    // gaussPts = fe_ptr->gaussPts; // FIXME:

    // print insideCells
    // std::cout << "Inside cells: " << insideCells << "\n";

  // auto calc_base_for_tet = [&]() {
  //   MoFEMFunctionBegin;
  //   // const size_t nb_gauss_pts = gaussPts.size2();
  //   const size_t nb_gauss_pts = QUAD_3D_TABLE[rule]->npoints;
  //   auto &gauss_pts = fe_ptr->gaussPts;
  //   gauss_pts.resize(4, nb_gauss_pts, false);

  //   auto &data = fe_ptr->getDataOnElementBySpaceArray()[H1];
  //   auto &base = data->dataOnEntities[MBVERTEX][0].getN(NOBASE);
  //   auto &diff_base = data->dataOnEntities[MBVERTEX][0].getDiffN(NOBASE);
  //   base.resize(nb_gauss_pts, 4, false);
  //   diff_base.resize(nb_gauss_pts, 12, false);
  //   CHKERR Tools::shapeFunMBTET(&*base.data().begin(), &gauss_pts(0, 0),
  //                               &gauss_pts(1, 0), &gauss_pts(2, 0), nb_gauss_pts);
  //   double *diff_shape_ptr = &*diff_base.data().begin();
  //   for (int gg = 0; gg != nb_gauss_pts; ++gg) {
  //     for (int nn = 0; nn != 4; ++nn) {
  //       for (int dd = 0; dd != 3; ++dd, ++diff_shape_ptr) {
  //         *diff_shape_ptr = Tools::diffShapeFunMBTET[3 * nn + dd];
  //       }
  //     }
  //   }
  //   MoFEMFunctionReturn(0);
  // };

    // auto set_integration_pts_for_tet = [&]() {
    //   MoFEMFunctionBegin;
    //   if (rule < QUAD_3D_TABLE_SIZE) {
    //     if (QUAD_3D_TABLE[rule]->dim != 3) {
    //       SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong dimension");
    //     }
    //     if (QUAD_3D_TABLE[rule]->order < rule) {
    //       SETERRQ2(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
    //                "wrong order %d != %d", QUAD_3D_TABLE[rule]->order, rule);
    //     }
    //     size_t nb_gauss_pts = QUAD_3D_TABLE[rule]->npoints;
    //     gaussPts.resize(4, nb_gauss_pts, false);
    //     cblas_dcopy(nb_gauss_pts, &QUAD_3D_TABLE[rule]->points[1], 4,
    //                 &gaussPts(0, 0), 1);
    //     cblas_dcopy(nb_gauss_pts, &QUAD_3D_TABLE[rule]->points[2], 4,
    //                 &gaussPts(1, 0), 1);
    //     cblas_dcopy(nb_gauss_pts, &QUAD_3D_TABLE[rule]->points[3], 4,
    //                 &gaussPts(2, 0), 1);
    //     cblas_dcopy(nb_gauss_pts, QUAD_3D_TABLE[rule]->weights, 1,
    //                 &gaussPts(3, 0), 1);

    //     // CHKERR calc_base_for_tet();

    //     auto &base = dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE);
    //     auto &diff_base = dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE);
    //     base.resize(nb_gauss_pts, 4, false);
    //     diff_base.resize(nb_gauss_pts, 12, false);
    //     double *shape_ptr = &*base.data().begin();
    //     cblas_dcopy(4 * nb_gauss_pts, QUAD_3D_TABLE[rule]->points, 1, shape_ptr,
    //                 1);
    //     double *diff_shape_ptr = &*diff_base.data().begin();
    //     for (int gg = 0; gg != nb_gauss_pts; ++gg) {
    //       for (int nn = 0; nn != 4; ++nn) {
    //         for (int dd = 0; dd != 3; ++dd, ++diff_shape_ptr) {
    //           *diff_shape_ptr = Tools::diffShapeFunMBTET[3 * nn + dd];
    //         }
    //       }
    //     }

    //   } else {
    //     SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
    //              "rule > quadrature order %d < %d", rule, QUAD_3D_TABLE_SIZE);
    //   }
    //   MoFEMFunctionReturn(0);
    // };

    auto set_base_quadrature = [&]() {
      MoFEMFunctionBegin;
      if (rule < QUAD_3D_TABLE_SIZE) {
        if (QUAD_3D_TABLE[rule]->dim != 3) {
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "wrong dimension");
        }
        if (QUAD_3D_TABLE[rule]->order < rule) {
          SETERRQ2(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                   "wrong order %d != %d", QUAD_3D_TABLE[rule]->order, rule);
        }
        const size_t nb_gauss_pts = QUAD_3D_TABLE[rule]->npoints;
        auto &gauss_pts = fe_ptr->gaussPts;
        gauss_pts.resize(4, nb_gauss_pts, false);
        cblas_dcopy(nb_gauss_pts, &QUAD_3D_TABLE[rule]->points[1], 4,
                    &gauss_pts(0, 0), 1);
        cblas_dcopy(nb_gauss_pts, &QUAD_3D_TABLE[rule]->points[2], 4,
                    &gauss_pts(1, 0), 1);
        cblas_dcopy(nb_gauss_pts, &QUAD_3D_TABLE[rule]->points[3], 4,
                    &gauss_pts(2, 0), 1);
        cblas_dcopy(nb_gauss_pts, QUAD_3D_TABLE[rule]->weights, 1,
                    &gauss_pts(3, 0), 1);
        // auto &data = fe_ptr->dataOnElement[H1];
        auto &data = fe_ptr->getDataOnElementBySpaceArray()[H1];
        // data->dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 4,
        //                                                       false);
        // double *shape_ptr =
        //     &*data->dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
        // cblas_dcopy(4 * nb_gauss_pts, QUAD_3D_TABLE[rule]->points, 1, shape_ptr,
        //             1);

        auto &base = data->dataOnEntities[MBVERTEX][0].getN(NOBASE);
        auto &diff_base = data->dataOnEntities[MBVERTEX][0].getDiffN(NOBASE);
        base.resize(nb_gauss_pts, 4, false);
        diff_base.resize(nb_gauss_pts, 12, false);
        CHKERR Tools::shapeFunMBTET(&*base.data().begin(), &gauss_pts(0, 0),
                                    &gauss_pts(1, 0), &gauss_pts(2, 0),
                                    nb_gauss_pts);
        double *diff_shape_ptr = &*diff_base.data().begin();
        for (int gg = 0; gg != nb_gauss_pts; ++gg) {
          for (int nn = 0; nn != 4; ++nn) {
            for (int dd = 0; dd != 3; ++dd, ++diff_shape_ptr) {
              *diff_shape_ptr = Tools::diffShapeFunMBTET[3 * nn + dd];
            }
          }
        }

      } else {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "rule > quadrature order %d < %d", rule, QUAD_3D_TABLE_SIZE);
      }
      MoFEMFunctionReturn(0);
    };

    auto set_integration_for_hex = [&]() {
      MoFEMFunctionBegin;
      CHKERR Tools::outerProductOfEdgeIntegrationPtsForHex(fe_ptr->gaussPts,
                                                           rule, rule, rule);
      MoFEMFunctionReturn(0);
    };

     auto calc_base_for_hex = [&](MatrixDouble &gaussPts) {
        MoFEMFunctionBegin;
        // gaussPts = fe_ptr->gaussPts;
        const size_t nb_gauss_pts = gaussPts.size2();
        auto &data = fe_ptr->getDataOnElementBySpaceArray()[H1];
        auto &base = data->dataOnEntities[MBVERTEX][0].getN(NOBASE);
        auto &diff_base = data->dataOnEntities[MBVERTEX][0].getDiffN(NOBASE);
        base.resize(nb_gauss_pts, 8, false);
        diff_base.resize(nb_gauss_pts, 24, false);
        for (int gg = 0; gg != nb_gauss_pts; ++gg) {
          const double ksi = gaussPts(0, gg);
          const double zeta = gaussPts(1, gg);
          const double eta = gaussPts(2, gg);
          base(gg, 0) = N_MBHEX0(ksi, zeta, eta);
          base(gg, 1) = N_MBHEX1(ksi, zeta, eta);
          base(gg, 2) = N_MBHEX2(ksi, zeta, eta);
          base(gg, 3) = N_MBHEX3(ksi, zeta, eta);
          base(gg, 4) = N_MBHEX4(ksi, zeta, eta);
          base(gg, 5) = N_MBHEX5(ksi, zeta, eta);
          base(gg, 6) = N_MBHEX6(ksi, zeta, eta);
          base(gg, 7) = N_MBHEX7(ksi, zeta, eta);
          diff_base(gg, 0 * 3 + 0) = diffN_MBHEX0x(zeta, eta);
          diff_base(gg, 1 * 3 + 0) = diffN_MBHEX1x(zeta, eta);
          diff_base(gg, 2 * 3 + 0) = diffN_MBHEX2x(zeta, eta);
          diff_base(gg, 3 * 3 + 0) = diffN_MBHEX3x(zeta, eta);
          diff_base(gg, 4 * 3 + 0) = diffN_MBHEX4x(zeta, eta);
          diff_base(gg, 5 * 3 + 0) = diffN_MBHEX5x(zeta, eta);
          diff_base(gg, 6 * 3 + 0) = diffN_MBHEX6x(zeta, eta);
          diff_base(gg, 7 * 3 + 0) = diffN_MBHEX7x(zeta, eta);
          diff_base(gg, 0 * 3 + 1) = diffN_MBHEX0y(ksi, eta);
          diff_base(gg, 1 * 3 + 1) = diffN_MBHEX1y(ksi, eta);
          diff_base(gg, 2 * 3 + 1) = diffN_MBHEX2y(ksi, eta);
          diff_base(gg, 3 * 3 + 1) = diffN_MBHEX3y(ksi, eta);
          diff_base(gg, 4 * 3 + 1) = diffN_MBHEX4y(ksi, eta);
          diff_base(gg, 5 * 3 + 1) = diffN_MBHEX5y(ksi, eta);
          diff_base(gg, 6 * 3 + 1) = diffN_MBHEX6y(ksi, eta);
          diff_base(gg, 7 * 3 + 1) = diffN_MBHEX7y(ksi, eta);
          diff_base(gg, 0 * 3 + 2) = diffN_MBHEX0z(ksi, zeta);
          diff_base(gg, 1 * 3 + 2) = diffN_MBHEX1z(ksi, zeta);
          diff_base(gg, 2 * 3 + 2) = diffN_MBHEX2z(ksi, zeta);
          diff_base(gg, 3 * 3 + 2) = diffN_MBHEX3z(ksi, zeta);
          diff_base(gg, 4 * 3 + 2) = diffN_MBHEX4z(ksi, zeta);
          diff_base(gg, 5 * 3 + 2) = diffN_MBHEX5z(ksi, zeta);
          diff_base(gg, 6 * 3 + 2) = diffN_MBHEX6z(ksi, zeta);
          diff_base(gg, 7 * 3 + 2) = diffN_MBHEX7z(ksi, zeta);
        }
        MoFEMFunctionReturn(0);
      };

     auto set_gauss_pts = [&](auto &ref_gauss_pts) {
       MoFEMFunctionBegin;
       fe_ptr->gaussPts.swap(ref_gauss_pts);
       const size_t nb_gauss_pts = fe_ptr->gaussPts.size2();
       auto &data = fe_ptr->getDataOnElementBySpaceArray()[H1];

      CHKERR calc_base_for_hex(fe_ptr->gaussPts);
      MoFEMFunctionReturnHot(0);
      //  CHKERR calc_base_for_hex(fe_ptr->gaussPts);
      //  CHKERR ShapeMBTET(shape_ptr, &fe_ptr->gaussPts(0, 0),
      //                    &fe_ptr->gaussPts(1, 0), &fe_ptr->gaussPts(2, 0),
      //                    nb_gauss_pts);
       // double diff_shape_fun[12];
       // CHKERR ShapeDiffMBTET(diff_shape_fun);

       // double *shape_ptr = 
       //     &*data->dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
       // double *diff_shape_fun =
       //     &*data->dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).data().begin();
       std::cout << " Setting shape funcs pts hex... \n";
       MoFEMFunctionReturnHot(0);

       // auto base = data->dataOnEntities[MBVERTEX][0].getBase();
       auto base = DEMKOWICZ_JACOBI_BASE;
       CHKERR HexPolynomialBase().getValue(
           fe_ptr->gaussPts,
           boost::shared_ptr<EntPolynomialBaseCtx>(new EntPolynomialBaseCtx(
               *data, "U", H1, CONTINUOUS, base, base)));

       MoFEMFunctionReturnHot(0);
       // CHKERR HexPolynomialBase().getValueH1(fe_ptr->gaussPts);
       // CHKERR ShapeMBTET(shape_ptr, &fe_ptr->gaussPts(0, 0),
       //                   &fe_ptr->gaussPts(1, 0), &fe_ptr->gaussPts(2, 0),
       //                   nb_gauss_pts);
       // std::cout << " Setting shape funcs pts2... \n";
       // double diff_shape_fun[12];
       // CHKERR ShapeDiffMBTET(diff_shape_fun);
       std::cout << " Setting shape funcs pts3... \n";

       MoFEMFunctionReturn(0);
     };

     // CHKERR set_base_quadrature();

     if (activeCells.find(fe_handle) != activeCells.end()) {
       std::cout << " Setting integration points for this element: "
                 << fe_handle << " " << type_from_handle(fe_handle) << "\n";
       auto refine_quadrature = [&]() {
         MoFEMFunctionBegin;

         std::function<Range(Range &)> get_all_children =
             [&](Range &ents) -> Range {
           Range all_children;
           CHKERR(bit_mng->updateRangeByChildren(ents, all_children));
           // if (!all_children.empty()) {
           //   ents.merge(get_all_children(all_children));
           // }

           Range children_level;
           for (auto &ent : all_children.subset_by_dimension(3)) {
             Range parent(ent, ent);
             children_level.merge(get_all_children(parent));
           }

           return children_level.empty() ? ents : children_level;
         };

         Range single_ent(fe_handle, fe_handle);

         auto all_refined_tets = get_all_children(single_ent);
         all_refined_tets = all_refined_tets.subset_by_dimension(3);
         // if (debug) {
         //   CHKERR save_range(moab, "ref_tets.vtk", tets);
         // }

         MatrixDouble ref_coords(all_refined_tets.size(), 12, false);
         std::cout << " Number of tets: " << all_refined_tets.size() << "\n";
         // std::cout << " all refined tets: " << all_refined_tets << "\n";
         int tt = 0;
         for (Range::iterator tit = all_refined_tets.begin();
              tit != all_refined_tets.end(); tit++, tt++) {
           int num_nodes;
           const EntityHandle *conn;
           CHKERR moab.get_connectivity(*tit, conn, num_nodes, true); // false
           CHKERR moab.get_coords(conn, num_nodes, &ref_coords(tt, 0));
         }
         std::cout << " We have the coords... \n";
         auto &data = fe_ptr->getDataOnElementBySpaceArray()[H1];
         // CHKERR calc_base_for_tet();
         // CHKERR set_base_quadrature();
         CHKERR set_integration_for_hex();
         // const size_t nb_gauss_pts = fe_ptr->gaussPts.size2();
         const size_t nb_gauss_pts = QUAD_3D_TABLE[rule]->npoints;
         CHKERR calc_base_for_hex(fe_ptr->gaussPts);

         // const size_t nb_gauss_pts = QUAD_3D_TABLE[rule]->npoints;
         MatrixDouble ref_gauss_pts(4, nb_gauss_pts * ref_coords.size1());
         MatrixDouble &shape_n = data->dataOnEntities[MBVERTEX][0].getN();
         std::cout << " Calculating shape funcs... \n";
         int gg = 0;
         for (size_t tt = 0; tt != ref_coords.size1(); tt++) {
           double *tet_coords = &ref_coords(tt, 0);
           double det = 6 * Tools::tetVolume(tet_coords);

           for (size_t ggg = 0; ggg != nb_gauss_pts; ++ggg, ++gg) {
             for (int dd = 0; dd != 3; dd++) {
               ref_gauss_pts(dd, gg) =
                   shape_n(ggg, 0) * tet_coords[3 * 0 + dd] +
                   shape_n(ggg, 1) * tet_coords[3 * 1 + dd] +
                   shape_n(ggg, 2) * tet_coords[3 * 2 + dd] +
                   shape_n(ggg, 3) * tet_coords[3 * 3 + dd];
             }
             ref_gauss_pts(3, gg) = fe_ptr->gaussPts(3, ggg) * det;
           }
         }

         std::cout << " Setting gauss pts... \n";
         CHKERR set_gauss_pts(ref_gauss_pts);

         std::cout << " Setting gauss pts2... \n";
         MoFEMFunctionReturn(0);
       };

       CHKERR refine_quadrature();
    } else if (insideCells.find(fe_handle) != insideCells.end()) {
      std::cout << " Setting STANDARD integration points for this element: "
                << fe_handle << " " << type_from_handle(fe_handle) << "\n";

      // CHKERR set_base_quadrature();
      // CHKERR setIntegrationPts();

      std::cout << " Calculating base for hex... \n";
      CHKERR set_integration_for_hex();
      CHKERR calc_base_for_hex(fe_ptr->gaussPts);
      std::cout << " Calculating base for hex2... \n";
    } else {
      // CHKERR set_integration_for_hex();
      // CHKERR calc_base_for_hex(fe_ptr->gaussPts);
      std::cout << " No integration points are set for this element: " << fe_handle << " " << type_from_handle(fe_handle) << "\n";
    }

    MoFEMFunctionReturn(0);
  }


private:
  struct Fe : public ForcesAndSourcesCore {
    using ForcesAndSourcesCore::dataOnElement;

  private:
    using ForcesAndSourcesCore::ForcesAndSourcesCore;
  };

  const Range &activeCells;
  const Range &insideCells;

  static inline std::map<long int, MatrixDouble> mapRefCoords;
};

enum ElementSide { LEFT_SIDE = 0, RIGHT_SIDE };
// data for skeleton computation
std::array<VectorInt, 2>
    indicesRowSideMap; ///< indices on rows for left hand-side
std::array<VectorInt, 2>
    indicesColSideMap; ///< indices on columns for left hand-side
std::array<MatrixDouble, 2> rowBaseSideMap; // base functions on rows
std::array<MatrixDouble, 2> colBaseSideMap; // base function  on columns
std::array<MatrixDouble, 2> rowDiffBaseSideMap; // derivative of base functions
std::array<MatrixDouble, 2> colDiffBaseSideMap; // derivative of base functions
std::array<double, 2> areaMap; // area/volume of elements on the side
std::array<int, 2> senseMap; // orientaton of local element edge/face in respect
                             // to global orientation of edge/face

Range activeCells;
Range insideCells;
Range activeSides;

struct SurfaceKDTree {

  MoFEM::Interface &mField;
  moab::Interface &moabEx;
  EntityHandle kdTreeRootMeshset;
  AdaptiveKDTree kdTree;

  SurfaceKDTree(MoFEM::Interface &m_field, moab::Interface &mb)
      : mField(m_field), moabEx(mb), kdTree(&mb) {}

  Range sKin; //< Skin mesh, i.e. mesh on surface
  Range sUrface;

  /**
   * \brief Take a skin from a mesh.
   *
   * @param  volume Volume elements in the mesh
   * @return        Error code
   */
  PetscErrorCode takeASkin(Range &volume) {
    PetscFunctionBegin;
    Skinner skin(&moabEx);
    CHKERR skin.find_skin(0, volume, false, sKin);
    
    // Only for debugging
    // EntityHandle meshset;
    // CHKERR moabEx.create_meshset(MESHSET_SET,meshset);
    // CHKERR moabEx.add_entities(meshset,sKin); CHKERRQ(rval);
    // CHKERR moabEx.write_file("skin.vtk","VTK","",&meshset,1);
    // CHKERRQ(rval);
    PetscFunctionReturn(0);
  }

  /**
   * \brief Build tree
   *
   * @return Error code
   */
  PetscErrorCode buildTree() {
    PetscFunctionBegin;
    CHKERR kdTree.build_tree(sKin, &kdTreeRootMeshset);
    
    PetscFunctionReturn(0);
  }

  PetscErrorCode buildTree(Range ents) {
    PetscFunctionBegin;
    CHKERR kdTree.build_tree(ents, &kdTreeRootMeshset);
    
    PetscFunctionReturn(0);
  }

  MoFEMErrorCode copySurface(const Range surface, moab::Interface &moab) {
    MoFEMFunctionBegin;
    sUrface.clear();
    std::map<EntityHandle, EntityHandle> verts_map;
    for (auto tit = surface.begin(); tit != surface.end(); tit++) {
      int num_nodes;
      const EntityHandle *conn;
      CHKERR moabEx.get_connectivity(*tit, conn, num_nodes, true);
      MatrixDouble coords(num_nodes, 3);

      CHKERR moabEx.get_coords(conn, num_nodes, &coords(0, 0));

      EntityHandle new_verts[num_nodes];
      for (int nn = 0; nn != num_nodes; nn++) {
        if (verts_map.find(conn[nn]) != verts_map.end()) {
          new_verts[nn] = verts_map[conn[nn]];
        } else {
          CHKERR moab.create_vertex(&coords(nn, 0), new_verts[nn]);
          verts_map[conn[nn]] = new_verts[nn];
        }
      }
      EntityHandle ele;
      CHKERR moab.create_element(MBTRI, new_verts, num_nodes, ele);
      sUrface.insert(ele);
    }
    // if (!save_mesh.empty())
    //   CHKERR SaveData(m_field.get_moab())(save_mesh, sUrface);
    MoFEMFunctionReturn(0);
  }

  /**
   * \brief Find point on surface which is closet to given point coordinates
   * @param  x  1st coordinate of point
   * @param  y  2nd coordinate of point
   * @param  z  3rd coordinate of point
   * @param  dx coordinate of distance vector
   * @param  dy coordinate of distance vector
   * @param  z  coordinate of distance vector
   * @return    Error code
   */
  PetscErrorCode findClosestPointToTheSurface(const double x, const double y,
                                              const double z, double &dx,
                                              double &dy, double &dz) {
    PetscFunctionBegin;
    const double coords[] = {x, y, z};
    double closest_point[3];
    EntityHandle triangle;
    CHKERR kdTree.closest_triangle(kdTreeRootMeshset, coords, closest_point,
                                   triangle);

    cblas_daxpy(3, -1, coords, 1, closest_point, 1);
    dx = closest_point[0];
    dy = closest_point[1];
    dz = closest_point[2];
    PetscFunctionReturn(0);
  }

  int isPointOutside(const double coords[3]) {
    std::vector<EntityHandle> triangles_out;
    vector<double> dist;
    double ray_dir[3] = {(double)rand(), (double)rand(), (double)rand()};
    const double tol = 1e-6; //
    CHKERR kdTree.ray_intersect_triangles(kdTreeRootMeshset, tol, ray_dir,
                                          coords, triangles_out, dist);

    if (triangles_out.size() % 2 == 0) {
      return 1;
    } else {
      return -1;
    }
  }

  MoFEMErrorCode getInsideCells(const Range &volume, Range &inside_cells) {
    MoFEMFunctionBegin;
    
    for (const auto &ent : volume) {
      const EntityHandle *conn;
      int num_nodes;
      CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
      int point_outside = 0;
      for (int nn = 0; nn != num_nodes; nn++) {
        double coords[3];
        CHKERR mField.get_moab().get_coords(&conn[nn], 1, coords);
        point_outside += isPointOutside(coords);
      }
      if (point_outside == num_nodes) 
        continue;
      inside_cells.insert(ent);
    }

    MoFEMFunctionReturn(0);
  }
  /**
   * \brief Set distance to vertices on mesh
   *
   * @param  field_name  Name of field for distances
   * @return            Error Code
   */
  PetscErrorCode setDistanceFromSurface(const std::string field_name) {
    PetscFunctionBegin;
    auto field_ents = mField.get_field_ents();
    auto field_bit_number = mField.get_field_bit_number(field_name);
    auto it = field_ents->lower_bound(
        FieldEntity::getLoBitNumberUId(field_bit_number));
    auto hi_it = field_ents->upper_bound(
        FieldEntity::getHiBitNumberUId(field_bit_number));
    // Get values at nodes first
    std::cout << " The field_bit_number is " << field_name << " "
              << field_bit_number << std::endl;
    for (; it != hi_it; it++) {
      // std::cout << "We enter the entities " << std::endl;
      EntityType type = it->get()->getEntType();
      if (type != MBVERTEX) {
        continue;
      }
      EntityHandle ent = it->get()->getEnt();
      double coords[3];
      CHKERR mField.get_moab().get_coords(&ent, 1, coords);
      int point_outside = isPointOutside(coords);

      double delta[3];
      CHKERR findClosestPointToTheSurface(coords[0], coords[1], coords[2],
                                          delta[0], delta[1], delta[2]);
      // std::cout << "the distance is " << delta[0] << " " << delta[1] << " " << delta[2] << std::endl;
      // Get vector of DOFs on vertex
      VectorAdaptor dofs = it->get()->getEntFieldData();
      // Set values
      dofs[0] = cblas_dnrm2(3, delta, 1) * point_outside;
    }
    // std::cout << "saved data on vertices at least" << std::endl;
    // Get values at edges and project those values on approximation base
    it = field_ents->lower_bound(
        FieldEntity::getLoBitNumberUId(field_bit_number));
    for (; it != hi_it; it++) {
      EntityType type = it->get()->getEntType();
      if (type != MBEDGE) {
        continue;
      }
      const EntityHandle *conn;
      int num_ndodes;
      EntityHandle ent = it->get()->getEnt();
      CHKERR mField.get_moab().get_connectivity(ent, conn, num_ndodes, true);
      
      double coords[3 * num_ndodes];
      CHKERR mField.get_moab().get_coords(conn, num_ndodes, coords);
      int point_outside = isPointOutside(coords);
      
      cblas_daxpy(3, 1, &coords[3], 1, coords, 1);
      cblas_dscal(3, 0.5, coords, 1);
      double delta[3];
      CHKERR findClosestPointToTheSurface(coords[0], coords[1], coords[2],
                                          delta[0], delta[1], delta[2]);

      auto conn_it0 = field_ents->find(
          FieldEntity::getLoLocalEntityBitNumber(field_bit_number, conn[0]));
      auto conn_it1 = field_ents->find(
          FieldEntity::getLoLocalEntityBitNumber(field_bit_number, conn[1]));

      if (conn_it0 == field_ents->end()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "entity not found");
      }
      if (conn_it1 == field_ents->end()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "entity not found");
      }

      VectorAdaptor vec_conn0 = conn_it0->get()->getEntFieldData();
      VectorAdaptor vec_conn1 = conn_it1->get()->getEntFieldData();
      // cerr << vec_conn0 << " " << vec_conn1 << endl;
      double ave_delta = 0.5 * (vec_conn0[0] + vec_conn1[0]);
      double edge_shape_function_val = 0.25;
      if (it->get()->getApproxBase() == AINSWORTH_LEGENDRE_BASE) {
      } else if (it->get()->getApproxBase() == AINSWORTH_LOBATTO_BASE) {
        edge_shape_function_val *= LOBATTO_PHI0(0);
      } 
      // else {
      //   SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Base not implemented");
      // }
      // Get vector of DOFs on vertex
      VectorAdaptor dofs = it->get()->getEntFieldData();
      // cerr << mid_val_dx << " " << ave_dx << endl;
      // Set values
      dofs[0] = point_outside *
          (cblas_dnrm2(3, delta, 1) - ave_delta) / edge_shape_function_val;
    }
    PetscFunctionReturn(0);
  }
};

/**
 * @brief Operator tp collect data from elements on the side of Edge/Face
 * 
 */
struct OpCalculateSideData : public FaceSideOp {

  OpCalculateSideData(std::string field_name, std::string col_field_name)
      : FaceSideOp(field_name, col_field_name, FaceSideOp::OPROWCOL) {

    std::fill(&doEntities[MBVERTEX], &doEntities[MBMAXTYPE], false);

    for (auto t = moab::CN::TypeDimensionMap[SPACE_DIM].first;
         t <= moab::CN::TypeDimensionMap[SPACE_DIM].second; ++t)
      doEntities[t] = true;
  }

  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type, EntData &row_data,
                        EntData &col_data) {
    MoFEMFunctionBeginHot;

    // Note: THat for L2 base data rows, and columns are the same, so operator
    // above can be simpler operator for the right hand side, and data can be
    // stored only for rows, since for columns data are the same. However for
    // complex multi-physics problems that not necessary would be a case. For
    // generality, we keep generic case.

    if ((CN::Dimension(row_type) == SPACE_DIM) &&
        (CN::Dimension(col_type) == SPACE_DIM)) {
      const auto nb_in_loop = getFEMethod()->nInTheLoop;
      indicesRowSideMap[nb_in_loop] = row_data.getIndices();
      indicesColSideMap[nb_in_loop] = col_data.getIndices();
      rowBaseSideMap[nb_in_loop] = row_data.getN();
      colBaseSideMap[nb_in_loop] = col_data.getN();
      rowDiffBaseSideMap[nb_in_loop] = row_data.getDiffN();
      colDiffBaseSideMap[nb_in_loop] = col_data.getDiffN();
      areaMap[nb_in_loop] = getMeasure();
      senseMap[nb_in_loop] = getSkeletonSense();
      if (!nb_in_loop) {
        indicesRowSideMap[1].clear();
        indicesColSideMap[1].clear();
        rowBaseSideMap[1].clear();
        colBaseSideMap[1].clear();
        rowDiffBaseSideMap[1].clear();
        colDiffBaseSideMap[1].clear();
        areaMap[1] = 0;
        senseMap[1] = 0;
      }
    } else {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Should not happen");
    }

    MoFEMFunctionReturnHot(0);
  }
};

template <typename T> inline auto get_ntensor(T &base_mat) {
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(
      &*base_mat.data().begin());
};

template <typename T> inline auto get_ntensor(T &base_mat, int gg, int bb) {
  double *ptr = &base_mat(gg, bb);
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(ptr);
};

template <typename T> inline auto get_diff_ntensor(T &base_mat) {
  double *ptr = &*base_mat.data().begin();
  return getFTensor1FromPtr<3>(ptr);
};

template <typename T>
inline auto get_diff_ntensor(T &base_mat, int gg, int bb) {
  double *ptr = &base_mat(gg, 3 * bb);
  return getFTensor1FromPtr<3>(ptr);
};


/**
 * @brief Operator the left hand side matrix 
 * 
 */
struct OpL2LhsGhostPenalty : public BoundaryEleOp {
public:

  /**
   * @brief Construct a new OpL2LhsGhostPenalty 
   * 
   * @param side_fe_ptr pointer to FE to evaluate side elements
   */
  OpL2LhsGhostPenalty(boost::shared_ptr<FaceSideEle> side_fe_ptr)
      : BoundaryEleOp(NOSPACE, BoundaryEleOp::OPSPACE), sideFEPtr(side_fe_ptr) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;

    // Collect data from side domian elements
    CHKERR loopSideVolumes("dFE", *sideFEPtr);
    // CHKERR loopSideFaces("dFE", *sideFEPtr);
    const auto in_the_loop =
        sideFEPtr->nInTheLoop; // return number of elements on the side

#ifndef NDEBUG
    const std::array<std::string, 2> ele_type_name = {"BOUNDARY", "SKELETON"};
    MOFEM_LOG("SELF", Sev::noisy)
        << "OpL2LhsPenalty inTheLoop " << ele_type_name[in_the_loop];
    MOFEM_LOG("SELF", Sev::noisy)
        << "OpL2LhsPenalty sense " << senseMap[0] << " " << senseMap[1];
#endif

    // calculate  penalty
    const double s = getMeasure() / (areaMap[0] + areaMap[1]);
    const double p = penalty * s;

    // get normal of the face or edge
    auto t_normal = getFTensor1Normal();
    auto inv_nrm = 1.0 / t_normal.l2();
    t_normal(i) *= inv_nrm;
    // t_normal.normalize();

    // get number of integration points on face
    const size_t nb_integration_pts = getGaussPts().size2();

    // beta paramter is zero, when penalty method is used, also, takes value 1,
    // when boundary edge/face is evaluated, and 2 when is skeleton edge/face.
    const double beta = static_cast<double>(nitsche) / (in_the_loop + 1);

    // iterate the sides rows
    for (auto s0 : {LEFT_SIDE, RIGHT_SIDE}) {

      // gent number of DOFs on the right side. 
      const auto nb_rows = indicesRowSideMap[s0].size();

      if (nb_rows) {

        // get orientation of the local element edge
        const auto sense_row = senseMap[s0];

        // iterate the side cols
        const auto nb_row_base_functions = rowBaseSideMap[s0].size2();
        for (auto s1 : {LEFT_SIDE, RIGHT_SIDE}) {

          // get orientation of the edge
          const auto sense_col = senseMap[s1];

          // number of dofs, for homogenous approximation this value is
          // constant.
          const auto nb_cols = indicesColSideMap[s1].size();

          // resize local element matrix
          locMat.resize(nb_rows, nb_cols, false);
          locMat.clear();

          // get base functions, and integration weights
          auto t_row_base = get_ntensor(rowBaseSideMap[s0]);
          auto t_diff_row_base = get_diff_ntensor(rowDiffBaseSideMap[s0]);
          auto t_w = getFTensor0IntegrationWeight();

          // iterate integration points on face/edge
          for (size_t gg = 0; gg != nb_integration_pts; ++gg) {

            // t_w is integration weight, and measure is element area, or
            // volume, depending if problem is in 2d/3d.
            const double alpha = getMeasure() * t_w;
            auto t_mat = locMat.data().begin();
            
            // iterate rows
            size_t rr = 0;
            for (; rr != nb_rows; ++rr) {

              // calculate tetting function 
              FTensor::Tensor1<double, SPACE_DIM> t_vn_plus;
              t_vn_plus(i) = beta * (phi * t_diff_row_base(i) / p);
              FTensor::Tensor1<double, SPACE_DIM> t_vn;
              t_vn(i) = t_row_base * t_normal(i) * sense_row - t_vn_plus(i);

              // get base functions on columns
              auto t_col_base = get_ntensor(colBaseSideMap[s1], gg, 0);
              auto t_diff_col_base =
                  get_diff_ntensor(colDiffBaseSideMap[s1], gg, 0);

              // iterate columns
              for (size_t cc = 0; cc != nb_cols; ++cc) {

                // calculate variance of tested function 
                FTensor::Tensor1<double, SPACE_DIM> t_un;
                t_un(i) = -p * (t_col_base * t_normal(i) * sense_col -
                                beta * t_diff_col_base(i) / p);

                // assemble matrix
                *t_mat -= alpha * (t_vn(i) * t_un(i));
                *t_mat -= alpha * (t_vn_plus(i) * (beta * t_diff_col_base(i)));

                // move to next column base and element of matrix
                ++t_col_base;
                ++t_diff_col_base;
                ++t_mat;
              }

              // move to next row base
              ++t_row_base;
              ++t_diff_row_base;
            }

            // this is obsolete for this particular example, we keep it for
            // generality. in case of multi-physics problems different fields can
            // chare hierarchical base but use different approx. order, so is
            // possible to have more base functions than DOFs on element.
            for (; rr < nb_row_base_functions; ++rr) {
              ++t_row_base;
              ++t_diff_row_base;
            }

            ++t_w;
          }

          // assemble system
          CHKERR ::MatSetValues(getKSPB(), indicesRowSideMap[s0].size(),
                                &*indicesRowSideMap[s0].begin(),
                                indicesColSideMap[s1].size(),
                                &*indicesColSideMap[s1].begin(),
                                &*locMat.data().begin(), ADD_VALUES);

          // stop of boundary element
          if (!in_the_loop)
            MoFEMFunctionReturnHot(0);
        }
      }
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<FaceSideEle>
      sideFEPtr; ///< pointer to element to get data on edge/face sides
  MatrixDouble locMat; ///< local operator matrix
};

}; // namespace Poisson3DCutFEMOperators

#endif //__POISSONCUTFEM_HPP__