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
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong result %3.4e", sum_coords);
      }
      if (fabs(1.0 - sum_gauss_pts) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
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
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong result %3.4e", sum_coords);
      }
      if (fabs(1.0 - sum_gauss_pts) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
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

      MOFEM_LOG("WORLD", Sev::verbose)
          << "sum_gauss_pts " << sum_gauss_pts << endl;
      MOFEM_LOG("WORLD", Sev::verbose) << "sum_coords " << sum_coords << endl;
      if (fabs(64.0 - sum_coords) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong result %3.4e", sum_coords);
      }
      if (fabs(1.0 - sum_gauss_pts) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong result %3.4e", sum_gauss_pts);
      }

      constexpr bool debug = false;

      if (debug) {

        auto [nodes, triangles, level_index] = refine;

        auto get_coords = [&](auto t) {
          auto [nodes, triangles, level_index] = refine;
          std::array<double, 9> ele_coords;
          for (auto n : {0, 1, 2}) {
            for (auto i : {0, 1}) {
              ele_coords[3 * n + i] = nodes[2 * triangles[6 * t + n] + i];
            }
            ele_coords[3 * n + 2] = 0;
          }
          return ele_coords;
        };

        moab::Core mb_instance;
        moab::Interface &moab = mb_instance;

        for (auto t = level_index[level_index.size() - 2];
             t != level_index.back(); ++t) {
          std::array<EntityHandle, 3> conn;
          auto ele_coords = get_coords(t);

          for (auto n : {0, 1, 2}) {
            CHKERR moab.create_vertex(&ele_coords[3 * n], conn[n]);
            MOFEM_LOG("WORLD", Sev::noisy)
                << ele_coords[3 * n] << " " << ele_coords[3 * n + 1] << " "
                << ele_coords[3 * n + 2];
          }

          MOFEM_LOG("WORLD", Sev::noisy)
              << triangles[6 * t + 0] << " " << triangles[6 * t + 1] << " "
              << triangles[6 * t + 2] << std::endl;

          EntityHandle tri;
          CHKERR moab.create_element(MBTRI, conn.data(), 3, tri);
        }

        for (int i = 0; i < new_nb_gauss_pts; ++i) {
          std::array<double, 3> coords{new_gauss_pts(0, i), new_gauss_pts(1, i),
                                       0};
          EntityHandle node;
          CHKERR moab.create_vertex(coords.data(), node);
        }

        CHKERR moab.write_file("gauss_refine.vtk", "VTK", "");
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
