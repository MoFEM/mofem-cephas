/** \file CGGTonsorialBubbleBase.hpp

  \brief Implementation of tonsorial bubble base div(v) = 0.

  Implementation is based and motiveted by \cite cockburn2010new. This base
  is used to approximate stresses using Hdiv base with weakly enforced
  symmetry.

*/

#include <Includes.hpp>
#include <version.h>
#include <config.h>
#include <definitions.h>
#include <FTensor.hpp>
#include <Common.hpp>
#include <h1_hdiv_hcurl_l2.h>
#include <CGGTonsorialBubbleBase.hpp>

using namespace MoFEM;
using namespace FTensor;

#ifndef GENERATE_VTK_CGG_BUBBLE_BASE

MoFEMErrorCode MoFEM::CGG_BubbleBase_MBTET(const int p, const double *N,
                                    const double *diffN, const double *l2_base,
                                    const double *diff_l2_base,
                                    Tensor1<PackPtr<double *, 3>, 3> &t_phi,
                                    const int gdim) {
  MoFEMFunctionBeginHot;

  Index<'i', 3> i;
  Index<'j', 3> j;

  Tensor1<double, 3> t_diff_n[4];
  {
    Tensor1<PackPtr<const double *, 3>, 3> t_diff_n_tmp(
        &diffN[0], &diffN[1], &diffN[2]);
    for (int ii = 0; ii != 4; ++ii) {
      t_diff_n[ii](i) = t_diff_n_tmp(i);
      ++t_diff_n_tmp;
    }
  }

  MatrixBoundedArray<double, 18> w(3, 6);
  MatrixBoundedArray<double, 54> diff_w(9, 6);

  Tensor0<PackPtr<const double *, 1>> t_l2(l2_base);
  Tensor1<PackPtr<const double *, 3>, 3> t_diff_l2(
      &diff_l2_base[0], &diff_l2_base[1], &diff_l2_base[2]);
  FTensor::Tensor1<PackPtr<const double *, 4>, 4> t_N(&N[0], &N[1], &N[2],
                                                      &N[3]);

  for (int gg = 0; gg != gdim; ++gg) {

    auto get_t_w = [&]() {
      return Tensor1<PackPtr<double *, 1>, 3>(&w(0, 0), &w(1, 0), &w(2, 0));
    };

    auto get_t_diff_w = [&]() {
      return Tensor2<PackPtr<double *, 1>, 3, 3>(
          &diff_w(0, 0), &diff_w(1, 0), &diff_w(2, 0), &diff_w(3, 0),
          &diff_w(4, 0), &diff_w(5, 0), &diff_w(6, 0), &diff_w(7, 0),
          &diff_w(8, 0));
    };

    const int bi[3][2][4] = {

        {{1, 2, 3, 0}, {0, 1, 2, 3}},

        {{2, 3, 0, 1}, {1, 2, 3, 0}},

        {{3, 0, 1, 2}, {2, 3, 0, 1}}

    };

    auto calc_t_w = [&](const int ii,const int kk) {
      const int ii3 = bi[ii][kk][0];
      const int ii2 = bi[ii][kk][1];
      const int ii1 = bi[ii][kk][2];
      const int ii0 = bi[ii][kk][3];
      const double a = t_N(ii3) * t_N(ii2) * t_N(ii1);
      Tensor1<double, 3> t_w;
      Index<'i', 3> i;
      t_w(i) = a * t_diff_n[ii0](i);
      return t_w;
    };

    auto calc_t_diff_w = [&](const int ii,const int kk) {
      const int ii3 = bi[ii][kk][0];
      const int ii2 = bi[ii][kk][1];
      const int ii1 = bi[ii][kk][2];
      const int ii0 = bi[ii][kk][3];
      Index<'i', 3> i;
      Index<'j', 3> j;
      Tensor1<double, 3> t_diff_a;
      t_diff_a(i) = t_diff_n[ii3](i) * (t_N(ii2) * t_N(ii1));
      t_diff_a(i) += t_diff_n[ii2](i) * (t_N(ii3) * t_N(ii1));
      t_diff_a(i) += t_diff_n[ii1](i) * (t_N(ii3) * t_N(ii2));
      Tensor2<double, 3, 3> t_diff_w;
      t_diff_w(i, j) = t_diff_n[ii0](i) * t_diff_a(j);
      return t_diff_w;
    };

    auto calc_t_w_and_t_diff_w = [&]() {
      auto t_w = get_t_w();
      auto t_diff_w = get_t_diff_w();
      for (int ii = 0; ii != 3; ++ii) {
        auto t_w_ii = calc_t_w(ii,0);
        auto t_diff_w_ii = calc_t_diff_w(ii,0);
        auto t_w_jj = calc_t_w(ii,1);
        auto t_diff_w_jj = calc_t_diff_w(ii,1);
        t_w(i) = t_w_ii(i) - t_w_jj(i);
        t_diff_w(i, j) = t_diff_w_ii(i, j) - t_diff_w_jj(i, j);
        ++t_w;
        ++t_diff_w;
      }
    };
    calc_t_w_and_t_diff_w();

    for (int oo = 0; oo < NBVOLUMETET_L2(p - 2); ++oo) {

      auto t_w = get_t_w();
      auto t_diff_w = get_t_diff_w();

      for (int ii = 0; ii != 3; ++ii) {
        Tensor2<double, 3, 3> t_diff_pw;
        t_diff_pw(i, j) = t_l2 * t_diff_w(i, j) + t_w(i) * t_diff_l2(j);
        t_phi(0) = t_diff_pw(2, 1) - t_diff_pw(1, 2);
        t_phi(1) = t_diff_pw(0, 2) - t_diff_pw(2, 0);
        t_phi(2) = t_diff_pw(1, 0) - t_diff_pw(0, 1);
        ++t_w;
        ++t_diff_w;
        ++t_phi;
      }

      ++t_l2;
      ++t_diff_l2;
    }

    ++t_N;
  }

  MoFEMFunctionReturnHot(0);
}

#endif // GENERATE_VTK_CGG_BUBBLE_BASE

#ifdef GENERATE_VTK_CGG_BUBBLE_BASE

#include <MoFEM.hpp>
using namespace MoFEM;
using namespace boost::numeric;

MoFEMErrorCode VTK_cgg_bubble_base_MBTET(const string file_name) {
  MoFEMFunctionBegin;

  double base_coords[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};

  moab::Core core_ref;
  moab::Interface &moab_ref = core_ref;

  EntityHandle nodes[4];
  for (int nn = 0; nn < 4; nn++) {
    CHKERR moab_ref.create_vertex(&base_coords[3 * nn], nodes[nn]);
  }
  EntityHandle tet;
  CHKERR moab_ref.create_element(MBTET, nodes, 4, tet);

  MoFEM::Core m_core_ref(moab_ref, PETSC_COMM_SELF, -2);
  MoFEM::Interface &m_field_ref = m_core_ref;

  CHKERR m_field_ref.getInterface<BitRefManager>()->setBitRefLevelByDim(
      0, 3, BitRefLevel().set(0));

  const int max_level = 4;
  for (int ll = 0; ll != max_level; ll++) {
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
    CHKERR m_ref->add_verices_in_the_middel_of_edges(edges,
                                                     BitRefLevel().set(ll + 1));
    CHKERR m_ref->refine_TET(tets, BitRefLevel().set(ll + 1));
  }

  Range tets;
  CHKERR m_field_ref.getInterface<BitRefManager>()
      ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(max_level),
                                     BitRefLevel().set(max_level), MBTET, tets);

  // Use 10 node tets to print base
  if (1) {
    EntityHandle meshset;
    CHKERR moab_ref.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, meshset);
    CHKERR moab_ref.add_entities(meshset, tets);
    CHKERR moab_ref.convert_entities(meshset, true, false, false);
    CHKERR moab_ref.delete_entities(&meshset, 1);
  }

  Range elem_nodes;
  CHKERR moab_ref.get_connectivity(tets, elem_nodes, false);

  const int nb_gauss_pts = elem_nodes.size();
  MatrixDouble gauss_pts(nb_gauss_pts, 4);
  gauss_pts.clear();
  Range::iterator nit = elem_nodes.begin();
  for (int gg = 0; nit != elem_nodes.end(); nit++, gg++) {
    CHKERR moab_ref.get_coords(&*nit, 1, &gauss_pts(gg, 0));
  }
  gauss_pts = trans(gauss_pts);

  MatrixDouble shape_fun;
  shape_fun.resize(nb_gauss_pts, 4);
  CHKERR ShapeMBTET(&*shape_fun.data().begin(), &gauss_pts(0, 0),
                    &gauss_pts(1, 0), &gauss_pts(2, 0), nb_gauss_pts);
  double diff_shape_fun[12];
  CHKERR ShapeDiffMBTET(diff_shape_fun);

  int p = 1;
  MatrixDouble l2_phi(nb_gauss_pts, NBVOLUMETET_L2(p));
  MatrixDouble l2_diff_phi(nb_gauss_pts,3 * NBVOLUMETET_L2(p));
  CHKERR L2_Ainsworth_ShapeFunctions_MBTET(p, &shape_fun(0, 0), diff_shape_fun,
                                           &l2_phi(0,0), &l2_diff_phi(0, 0),
                                           nb_gauss_pts, Legendre_polynomials);

  MatrixDouble phi(nb_gauss_pts, 3 * NBVOLUMETET_CCG_BUBBLE(p));
  Tensor1<PackPtr<double *, 3>, 3> t_phi(&phi(0,0), &phi(0,1), &phi(0,2));
  CHKERR CGG_BubbleBase_MBTET(
      p+2, &shape_fun(0, 0), diff_shape_fun, &l2_phi(0,0), &l2_diff_phi(0, 0),
      t_phi, nb_gauss_pts);


  for (int ll = 0; ll != NBVOLUMETET_L2(p); ++ll) {

    double def_val[] = {0};
    std::string tag_name = "L2_" + boost::lexical_cast<std::string>(ll);
    Tag th;
    CHKERR moab_ref.tag_get_handle(tag_name.c_str(), 1, MB_TYPE_DOUBLE, th,
                                   MB_TAG_CREAT | MB_TAG_SPARSE, def_val);

    int gg = 0;
    for (Range::iterator nit = elem_nodes.begin(); nit != elem_nodes.end();
         nit++, gg++) {
      int idx = ll;
      double data[] = {l2_phi(gg, idx)};
      CHKERR moab_ref.tag_set_data(th, &*nit, 1, data);
    }

  }

  for (int ll = 0; ll != NBVOLUMETET_L2(p); ++ll) {

    double def_val[] = {0};
    std::string tag_name = "DiffL2_" + boost::lexical_cast<std::string>(ll);
    Tag th;
    CHKERR moab_ref.tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, th,
                                   MB_TAG_CREAT | MB_TAG_SPARSE, def_val);

    int gg = 0;
    for (Range::iterator nit = elem_nodes.begin(); nit != elem_nodes.end();
         nit++, gg++) {
      int idx = 3*ll;
      double data[] = {l2_diff_phi(gg, idx + 0), l2_diff_phi(gg, idx + 1),
                       l2_diff_phi(gg, idx + 2)};
      CHKERR moab_ref.tag_set_data(th, &*nit, 1, data);
    }

  }

  for (int ll = 0; ll != NBVOLUMETET_CCG_BUBBLE(p); ++ll) {

    double def_val[] = {0, 0, 0};
    std::string tag_name = "B_" + boost::lexical_cast<std::string>(ll);
    Tag th;
    CHKERR moab_ref.tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, th,
                                   MB_TAG_CREAT | MB_TAG_SPARSE, def_val);

    int gg = 0;
    for (Range::iterator nit = elem_nodes.begin(); nit != elem_nodes.end();
         nit++, gg++) {
      int idx = 3 * ll;
      double data[] = {phi(gg, idx + 0), phi(gg, idx + 1), phi(gg, idx + 2)};
      CHKERR moab_ref.tag_set_data(th, &*nit, 1, data);
    }

  }

  EntityHandle meshset;
  CHKERR moab_ref.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, meshset);
  CHKERR moab_ref.add_entities(meshset, tets);
  CHKERR moab_ref.write_file(file_name.c_str(), "VTK", "", &meshset, 1);

  MoFEMFunctionReturn(0);
}

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {
    CHKERR VTK_cgg_bubble_base_MBTET("out_curl_vtk_cgg_bubble_base_on_tet.vtk");
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}

#endif // GENERATE_VTK_CGG_BUBBLE_BASE