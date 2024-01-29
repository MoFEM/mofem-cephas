/** \file gauss_points_on_outer_product.cpp
  \example gauss_points_on_outer_product.cpp
  \brief Testing gauss points coordinates and weights for a quad face

*/

#include <MoFEM.hpp>

#ifdef __cplusplus
extern "C" {
#endif
#include <cblas.h>
#include <quad.h>
#ifdef __cplusplus
}
#endif

using namespace MoFEM;
static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);
  try {
    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    MoFEM::Core core(moab);

    auto test_quad = [&]() {
      MoFEMFunctionBegin;
      MatrixDouble pts_quad;
      int rule_ksi = 6;
      int rule_eta = 8;
      CHKERR Tools::outerProductOfEdgeIntegrationPtsForQuad(pts_quad, rule_ksi,
                                                            rule_eta);
      int nb_gauss_pts = pts_quad.size2();
      double sum_coords = 0, sum_gauss_pts = 0;
      for (int i = 0; i < nb_gauss_pts; ++i) {
        sum_coords += pts_quad(0, i) + pts_quad(1, i);
        sum_gauss_pts += pts_quad(2, i);
      }
      double eps = 1e-8;
      if (fabs(20.0 - sum_coords) > eps) {
        SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong result %3.4e", sum_coords);
      }
      if (fabs(1.0 - sum_gauss_pts) > eps) {
        SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong result %3.4e", sum_gauss_pts);
      }
      MoFEMFunctionReturn(0);
    };

    auto test_hex = [&]() {
      MoFEMFunctionBegin;
      MatrixDouble pts_quad;
      int rule_ksi = 4;
      int rule_eta = 5;
      int rule_zet = 6;
      CHKERR Tools::outerProductOfEdgeIntegrationPtsForHex(pts_quad, rule_ksi,
                                                           rule_eta, rule_zet);
      int nb_gauss_pts = pts_quad.size2();
      double sum_coords = 0, sum_gauss_pts = 0;
      for (int i = 0; i < nb_gauss_pts; ++i) {
        sum_coords += pts_quad(0, i) + pts_quad(1, i) + pts_quad(2, i);
        sum_gauss_pts += pts_quad(3, i);
      }
      double eps = 1e-8;
      if (fabs(54.0 - sum_coords) > eps) {
        SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong result %3.4e", sum_coords);
      }
      if (fabs(1.0 - sum_gauss_pts) > eps) {
        SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong result %3.4e", sum_gauss_pts);
      }
      MoFEMFunctionReturn(0);
    };

    auto test_refined_triangle = [&]() {
      MoFEMFunctionBegin;

      auto refine = Tools::refineTriangle(2);
      auto new_gauss_pts = Tools::refineTriangleIntegrationPts(4, refine);

      int new_nb_gauss_pts = new_gauss_pts.size2();
      double sum_coords = 0, sum_gauss_pts = 0;
      for (int i = 0; i < new_nb_gauss_pts; ++i) {
        sum_coords += new_gauss_pts(0, i) + new_gauss_pts(1, i);
        sum_gauss_pts += new_gauss_pts(2, i);
      }
      constexpr double eps = 1e-8;
      // if (fabs(54.0 - sum_coords) > eps) {
      //   SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
      //            "wrong result %3.4e", sum_coords);
      // }
      if (fabs(1.0 - sum_gauss_pts) > eps) {
        SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong result %3.4e", sum_gauss_pts);
      }


      MoFEMFunctionReturn(0);
    };

    CHKERR test_quad();
    CHKERR test_hex();
    CHKERR test_refined_triangle();
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
