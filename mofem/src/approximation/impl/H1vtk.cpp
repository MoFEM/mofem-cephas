/** \file Hcurl.cpp

  \brief Implementation of H-curl base

  Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by Mark Ainsworth and Joe Coyle and by Demkowicz
  Shape functions for MBTRI/MBTET and HCurl space

*/

#include <Includes.hpp>
#include <FTensor.hpp>
#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>
#include <definitions.h>
#include <Common.hpp>
#include <H1vtk.hpp>

using namespace MoFEM;

#ifdef GENERATE_VTK_WITH_H1_BASE

 #include <MoFEM.hpp>
using namespace MoFEM;
using namespace boost::numeric;

MoFEMErrorCode VTK_Demkowicz_H1_MBHEX(const string file_name) {
  MoFEMFunctionBegin;

  double base_coords[] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                          0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};

  moab::Core core_ref;
  moab::Interface &moab_ref = core_ref;

  EntityHandle nodes[8];
  for (int nn = 0; nn < 8; nn++) {
    CHKERR moab_ref.create_vertex(&base_coords[3 * nn], nodes[nn]);
  }
  EntityHandle hex;
  CHKERR moab_ref.create_element(MBHEX, nodes, 8, hex);

  MoFEM::Core m_core_ref(moab_ref, PETSC_COMM_SELF, -2);
  MoFEM::Interface &m_field_ref = m_core_ref;

  // CHKERR m_field_ref.getInterface<BitRefManager>()->setBitRefLevelByDim(
  //     0, 3, BitRefLevel().set(0));

  const int max_level = 3;
  // for (int ll = 0; ll != max_level; ll++) {
  //   Range edges;
  //   CHKERR m_field_ref.getInterface<BitRefManager>()
  //       ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(ll),
  //                                      BitRefLevel().set(), MBEDGE, edges);
  //   Range hexes;
  //   CHKERR m_field_ref.getInterface<BitRefManager>()
  //       ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(ll),
  //                                      BitRefLevel(ll).set(), MBHEX, hexes);
  //   // refine mesh
  //   MeshRefinement *m_ref;
  //   CHKERR m_field_ref.getInterface(m_ref);
  //   CHKERR m_ref->add_verices_in_the_middel_of_edges(edges,
  //                                                    BitRefLevel().set(ll + 1));
  //   CHKERR m_ref->refine_TET(hexes, BitRefLevel().set(ll + 1));
  // }

  Range hexes;
  CHKERR moab_ref.get_entities_by_type(0, MBHEX, hexes);
  // CHKERR m_field_ref.getInterface<BitRefManager>()
  //     ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(max_level),
  //                                    BitRefLevel().set(max_level), MBHEX,
  //                                    hexes);

  // Use 10 node hexes to print base
  if (1) {
    EntityHandle meshset;
    CHKERR moab_ref.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, meshset);
    CHKERR moab_ref.add_entities(meshset, hexes);
    CHKERR moab_ref.convert_entities(meshset, false, false, false);
    CHKERR moab_ref.delete_entities(&meshset, 1);
  }

  Range elem_nodes;
  rval = moab_ref.get_connectivity(hexes, elem_nodes, false);
  CHKERRG(rval);

  const int nb_gauss_pts = elem_nodes.size();
  MatrixDouble gauss_pts(nb_gauss_pts, 8);
  gauss_pts.clear();
  Range::iterator nit = elem_nodes.begin();
  for (int gg = 0; nit != elem_nodes.end(); nit++, gg++) {
    rval = moab_ref.get_coords(&*nit, 1, &gauss_pts(gg, 0));
    CHKERRG(rval);
  }
  gauss_pts = trans(gauss_pts);

  MatrixDouble shape_fun;
  shape_fun.resize(nb_gauss_pts, 8);
  CHKERR ShapeMBHEX(&*shape_fun.data().begin(), &gauss_pts(0, 0),
                    &gauss_pts(1, 0), &gauss_pts(2, 0), nb_gauss_pts);

  MatrixDouble diff_shape_fun;
  diff_shape_fun.resize(nb_gauss_pts, 24);
  CHKERR ShapeDiffMBHEX(&*diff_shape_fun.data().begin(), &gauss_pts(0, 0),
                        &gauss_pts(1, 0), &gauss_pts(2, 0), nb_gauss_pts);

  int edge_sense[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  const int order = 2;
  int edge_order[12] = {order, order, order, order, order, order, order, order, order, order, order, order};

  
  //begin
  // MatrixDouble edge_phi(6, 3 * NBEDGE_DEMKOWICZ_HCURL(order) * nb_gauss_pts);
  // MatrixDouble edge_diff_phi(6,
  //                            9 * NBEDGE_DEMKOWICZ_HCURL(order) * nb_gauss_pts);

  // edge_phi.clear();
  // edge_diff_phi.clear();

  // double *edge_phi_ptr[] = {&edge_phi(0, 0), &edge_phi(1, 0), &edge_phi(2, 0),
  //                           &edge_phi(3, 0), &edge_phi(4, 0), &edge_phi(5, 0)};
  // double *edge_diff_phi_ptr[] = {&edge_diff_phi(0, 0), &edge_diff_phi(1, 0),
  //                                &edge_diff_phi(2, 0), &edge_diff_phi(3, 0),
  //                                &edge_diff_phi(4, 0), &edge_diff_phi(5, 0)};

  // CHKERR H1_EdgeShapeFunctions_MBHEX(
  //     edge_sense, edge_order, &*shape_fun.data().begin(), diff_shape_fun,
  //     edge_phi_ptr, edge_diff_phi_ptr, nb_gauss_pts);
//end




// PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
//                                      double *diffL, const int dim) =
//       cTx->basePolynomialsType0;
  
    // edges
  
  double *h1_edge_n[12], *diff_h1_edge_n[12];
  //for (int ee = 0; ee != 12; ++ee) {

    int nb_dofs = NBEDGE_H1(order);

    MatrixDouble edge_phi(12, 3 * nb_dofs * nb_gauss_pts);
    MatrixDouble edge_diff_phi(12, 9 * nb_dofs * nb_gauss_pts);

    edge_phi.clear();
    edge_diff_phi.clear();

    double *edge_phi_ptr[] = {
        &edge_phi(0, 0), &edge_phi(1, 0), &edge_phi(2, 0),  &edge_phi(3, 0),
        &edge_phi(4, 0), &edge_phi(5, 0), &edge_phi(6, 0),  &edge_phi(7, 0),
        &edge_phi(8, 0), &edge_phi(9, 0), &edge_phi(10, 0), &edge_phi(11, 0)};

    double *edge_diff_phi_ptr[] = {
        &edge_diff_phi(0, 0), &edge_diff_phi(1, 0), &edge_diff_phi(2, 0),
        &edge_diff_phi(3, 0), &edge_diff_phi(4, 0), &edge_diff_phi(5, 0),
        &edge_diff_phi(6, 0), &edge_diff_phi(7, 0), &edge_diff_phi(8, 0),
        &edge_diff_phi(9, 0), &edge_diff_phi(10, 0), &edge_diff_phi(11, 0)};

    //}
    CHKERR H1_EdgeShapeFunctions_MBHEX(
        edge_sense, edge_order,
        &*shape_fun.data().begin(),
        &*diff_shape_fun.data().begin(),
        edge_phi_ptr, edge_diff_phi_ptr, nb_gauss_pts, Jacobi_polynomials);




  for (int ee = 0; ee != 12; ++ee) {
      for (int ll = 0; ll != NBEDGE_H1(edge_order[ee]); ++ll) {
      double def_val[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
      std::string tag_name = "E" + boost::lexical_cast<std::string>(ee) + "_" +
           boost::lexical_cast<std::string>(ll);
      Tag th;
      CHKERR moab_ref.tag_get_handle(
          tag_name.c_str(),
          3, MB_TYPE_DOUBLE, th, MB_TAG_CREAT | MB_TAG_SPARSE, def_val);

      int gg = 0;
      for (Range::iterator nit = elem_nodes.begin(); nit != elem_nodes.end();
           nit++, gg++) {
        int idx = 3 * NBEDGE_H1(edge_order[ee]) * gg + 3 * ll;
        CHKERR moab_ref.tag_set_data(th, &*nit, 1, &edge_phi(ee,idx));
      }
    }
  }

  // int faces_order[] = {order, order, order, order};
  // int faces_nodes[] = {0, 1, 3, 1, 2, 3, 0, 2, 3, 0, 1, 2};

  // MatrixDouble face_phi(4, 3 * NBFACETRI_DEMKOWICZ_HCURL(order) * nb_gauss_pts);
  // MatrixDouble face_diff_phi(4, 9 * NBFACETRI_DEMKOWICZ_HCURL(order) *
  //                                   nb_gauss_pts);
  // face_phi.clear();
  // face_diff_phi.clear();

  // double *face_phi_ptr[] = {&face_phi(0, 0), &face_phi(1, 0), &face_phi(2, 0),
  //                           &face_phi(3, 0)};
  // double *face_diff_phi_ptr[] = {&face_diff_phi(0, 0), &face_diff_phi(1, 0),
  //                                &face_diff_phi(2, 0), &face_diff_phi(3, 0)};

  // CHKERR Hcurl_Demkowicz_FaceBaseFunctions_MBTET(
  //     faces_nodes, faces_order, &*shape_fun.data().begin(), diff_shape_fun,
  //     face_phi_ptr, face_diff_phi_ptr, nb_gauss_pts);

  // for (int ff = 0; ff != 4; ++ff) {
  //   for (int ll = 0; ll != NBFACETRI_DEMKOWICZ_HCURL(order); ++ll) {
  //     double def_val[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  //     std::string tag_name = "F" + boost::lexical_cast<std::string>(ff) + "_" +
  //                            boost::lexical_cast<std::string>(ll);
  //     Tag th;
  //     CHKERR moab_ref.tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, th,
  //                                    MB_TAG_CREAT | MB_TAG_SPARSE, def_val);

  //     int gg = 0;
  //     for (Range::iterator nit = elem_nodes.begin(); nit != elem_nodes.end();
  //          nit++, gg++) {
  //       int idx = 3 * NBFACETRI_DEMKOWICZ_HCURL(order) * gg + 3 * ll;
  //       CHKERR moab_ref.tag_set_data(th, &*nit, 1, &face_phi(ff, idx));
  //     }
  //   }
  // }

  // MatrixDouble vol_phi(nb_gauss_pts, 3 * NBVOLUMETET_DEMKOWICZ_HCURL(order));
  // MatrixDouble diff_vol_phi(nb_gauss_pts,
  //                           9 * NBVOLUMETET_DEMKOWICZ_HCURL(order));

  // CHKERR Hcurl_Demkowicz_VolumeBaseFunctions_MBHEX(
  //   order, &*shape_fun.data().begin(), diff_shape_fun, &vol_phi(0,0),
  //   &diff_vol_phi(0,0), nb_gauss_pts);

  // for (int ll = 0; ll != NBVOLUMETET_DEMKOWICZ_HCURL(order); ++ll) {
  //   double def_val[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  //   std::string tag_name = "V_" + boost::lexical_cast<std::string>(ll);
  //   Tag th;
  //   CHKERR moab_ref.tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, th,
  //                                  MB_TAG_CREAT | MB_TAG_SPARSE, def_val);

  //   int gg = 0;
  //   for (Range::iterator nit = elem_nodes.begin(); nit != elem_nodes.end();
  //        nit++, gg++) {
  //     int idx = 3 * ll;
  //     CHKERR moab_ref.tag_set_data(th, &*nit, 1, &vol_phi(gg, idx));
  //   }
  // }

  EntityHandle meshset;
  CHKERR moab_ref.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, meshset);
  CHKERR moab_ref.add_entities(meshset, hexes);
  CHKERR moab_ref.write_file(file_name.c_str(), "VTK", "", &meshset, 1);

  MoFEMFunctionReturn(0);
}

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {
    // CHKERR VTK_Ainsworth_Hcurl_MBTET("out_curl_vtk_ainsworth_base_on_tet.vtk");
    CHKERR VTK_Demkowicz_H1_MBHEX("out_h1_vtk_demkowicz_base_on_hex.vtk");
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}

#endif // GENERATE_VTK_WITH_H1_BASE

