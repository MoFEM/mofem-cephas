/**
 * \file bernstein_bezier_generate_base.cpp
 * \example bernstein_bezier_generate_base.cpp
 *
 * Genarte and check Bernstein-Bezier base. Test validates three properties
 * of BB polynomials from \cite ainsworth2011bernstein. 
 * 
 * 1) Multiplication
 * 2) Integration
 * 3) Derivative
 *
 */

/* This file is part of MoFEM.
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

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

    // Create MoFEM instance
    MoFEM::Core core(moab);

    constexpr int N = 5;

    // Edge

    MatrixInt edge_alpha(2 + NBEDGE_H1(N), 2);
    CHKERR BernsteinBezier::generateIndicesVertexEdge(N, &edge_alpha(0, 0));
    CHKERR BernsteinBezier::generateIndicesEdgeEdge(N, &edge_alpha(2, 0));
    std::cout << "edge alpha " << edge_alpha << std::endl;

    std::array<double, 6> edge_x_k = {0, 0, 0, 1, 0, 0};
    MatrixDouble edge_x_alpha(edge_alpha.size1(), 3);
    CHKERR BernsteinBezier::domainPoints3d(N, 2, edge_alpha.size1(),
                                           &edge_alpha(0, 0), edge_x_k.data(),
                                           &edge_x_alpha(0, 0));
    std::cout << "domain points " << edge_x_alpha << endl;

    const int M = 50;
    MatrixDouble edge_base(M, edge_alpha.size1());
    MatrixDouble edge_diff_base(M, edge_alpha.size1());

    auto calc_lambda_on_edge = [](int M) {
      MatrixDouble edge_lambda(M, 2);
      for (size_t i = 0; i != M; ++i) {
        double x = static_cast<double>(i) / (M - 1);
        edge_lambda(i, 0) = N_MBEDGE0(x);
        edge_lambda(i, 1) = N_MBEDGE1(x);
      }
      return edge_lambda;
    };
    auto edge_lambda = calc_lambda_on_edge(M);

    CHKERR BernsteinBezier::baseFunctionsEdge(
        N, M, edge_alpha.size1(), &edge_alpha(0, 0), &edge_lambda(0, 0),
        Tools::diffShapeFunMBEDGE.data(), &edge_base(0, 0),
        &edge_diff_base(0, 0));

    auto print_edge_base = [](auto M, auto &edge_base, auto &edge_diff_base) {
      MoFEMFunctionBegin;
      for (size_t i = 0; i != M; ++i) {
        double x = static_cast<double>(i) / (M - 1);
        std::cout << "edge " << x << " ";
        for (size_t j = 0; j != edge_base.size2(); ++j)
          std::cout << edge_base(i, j) << " ";
        std::cout << endl;
      }

      for (size_t i = 0; i != M; ++i) {
        double x = static_cast<double>(i) / (M - 1);
        std::cout << "diff_edge " << x << " ";
        for (size_t j = 0; j != edge_diff_base.size2(); ++j)
          std::cout << edge_diff_base(i, j) << " ";
        std::cout << endl;
      }
      MoFEMFunctionReturn(0);
    };
    CHKERR print_edge_base(M, edge_base, edge_diff_base);

    auto diag_n = [](int n) { return n * (n + 1) / 2; };

    auto binomial_alpha_beta = [](int dim, const int *alpha, const int *beta) {
      double f = boost::math::binomial_coefficient<double>(alpha[0] + beta[0],
                                                           alpha[0]);
      for (int d = 1; d != dim; ++d) {
        f *= boost::math::binomial_coefficient<double>(alpha[d] + beta[d],
                                                       alpha[d]);
      }
      return f;
    };

    auto check_property_one =
        [N, diag_n, binomial_alpha_beta](
            const int M, const auto &edge_alpha, const auto &edge_lambda,
            auto &edge_base,
            boost::function<MoFEMErrorCode(
                const int N, const int gdim, const int n_alpha,
                const int *alpha, const double *lambda,
                const double *grad_lambda, double *base, double *grad_base)>
                calc_base,
            bool debug) {
          MoFEMFunctionBegin;

          MatrixInt edge_alpha2(diag_n(edge_alpha.size1()), edge_alpha.size2());
          int k = 0;
          for (int i = 0; i != edge_alpha.size1(); ++i) {
            for (int j = i; j != edge_alpha.size1(); ++j, ++k) {
              for (int d = 0; d != edge_alpha.size2(); ++d) {
                edge_alpha2(k, d) = edge_alpha(i, d) + edge_alpha(j, d);
              }
            }
          }

          MatrixDouble edge_base2(edge_base.size1(), edge_alpha2.size1());
          CHKERR calc_base(N + N, edge_base.size1(), edge_alpha2.size1(),
                           &edge_alpha2(0, 0), &edge_lambda(0, 0),
                           Tools::diffShapeFunMBEDGE.data(), &edge_base2(0, 0),
                           nullptr);

          const double f0 = boost::math::binomial_coefficient<double>(N + N, N);
          for (int g = 0; g != edge_base.size1(); ++g) {
            int k = 0;
            for (size_t i = 0; i != edge_alpha.size1(); ++i) {
              for (size_t j = i; j != edge_alpha.size1(); ++j, ++k) {
                const double f = binomial_alpha_beta(
                    edge_alpha.size2(), &edge_alpha(i, 0), &edge_alpha(j, 0));
                const double b = f / f0;
                const double B_check = b * edge_base2(g, k);
                const double B_mult = edge_base(g, i) * edge_base(g, j);
                const double error = std::abs(B_check - B_mult);
                if (debug)
                  std::cout << "( " << k << " " << i << " " << j << " )"
                            << B_check << " " << B_mult << " " << error
                            << std::endl;
                constexpr double eps = 1e-12;
                if (error > eps)
                  SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                          "Property one is not working");
              }
            }
          }
          MoFEMFunctionReturn(0);
        };

    CHKERR check_property_one(M, edge_alpha, edge_lambda, edge_base,
                              BernsteinBezier::baseFunctionsEdge, false);

    auto check_property_two_on_edge = [N](auto &edge_alpha, bool debug) {
      MoFEMFunctionBegin;
      const int rule = N;
      int nb_gauss_pts = QUAD_1D_TABLE[rule]->npoints;
      MatrixDouble gauss_pts(2, nb_gauss_pts, false);
      cblas_dcopy(nb_gauss_pts, &QUAD_1D_TABLE[rule]->points[1], 2,
                  &gauss_pts(0, 0), 1);
      cblas_dcopy(nb_gauss_pts, QUAD_1D_TABLE[rule]->weights, 1,
                  &gauss_pts(1, 0), 1);
      MatrixDouble edge_lambda(nb_gauss_pts, 2);
      cblas_dcopy(2 * nb_gauss_pts, QUAD_1D_TABLE[rule]->points, 1,
                  &edge_lambda(0, 0), 1);
      MatrixDouble edge_base(nb_gauss_pts, edge_alpha.size1());
      CHKERR BernsteinBezier::baseFunctionsEdge(
          N, nb_gauss_pts, edge_alpha.size1(), &edge_alpha(0, 0),
          &edge_lambda(0, 0), Tools::diffShapeFunMBEDGE.data(),
          &edge_base(0, 0), nullptr);

      VectorDouble integral(edge_alpha.size1());
      integral.clear();

      for (int g = 0; g != nb_gauss_pts; ++g) {
        for (size_t i = 0; i != edge_alpha.size1(); ++i) {
          integral[i] += gauss_pts(1, g) * edge_base(g, i);
        }
      }

      const double check_integral =
          1 / boost::math::binomial_coefficient<double>(N + 1, 1);

      for (auto i : integral) {
        const double error = std::abs(i - check_integral);
        if (debug)
          std::cout << "edge integral " << i << " " << check_integral << " "
                    << error << endl;
        constexpr double eps = 1e-12;
        if (error > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "Property two on edge is not working");
      }

      MoFEMFunctionReturn(0);
    };

    CHKERR check_property_two_on_edge(edge_alpha, false);

    auto check_property_three_on_edge_derivatives =
        [N](const int M, const auto &edge_alpha, const auto &edge_lambda,
            const auto &edge_diff_base, bool debug) {
          MoFEMFunctionBegin;

          MatrixInt edge_alpha_diff(2 + NBEDGE_H1(N - 1), 2);
          CHKERR BernsteinBezier::generateIndicesVertexEdge(
              N - 1, &edge_alpha_diff(0, 0));
          CHKERR BernsteinBezier::generateIndicesEdgeEdge(
              N - 1, &edge_alpha_diff(2, 0));

          MatrixDouble c0(2 + NBEDGE_H1(N), 2 + NBEDGE_H1(N - 1));
          std::array<int, 2> diff0 = {1, 0};
          CHKERR BernsteinBezier::genrateDerivativeIndicesEdges(
              N, edge_alpha.size1(), &edge_alpha(0, 0), diff0.data(),
              edge_alpha_diff.size1(), &edge_alpha_diff(0, 0), &c0(0, 0));
          MatrixDouble c1(2 + NBEDGE_H1(N), 2 + NBEDGE_H1(N - 1));
          std::array<int, 2> diff1 = {0, 1};
          CHKERR BernsteinBezier::genrateDerivativeIndicesEdges(
              N, edge_alpha.size1(), &edge_alpha(0, 0), diff1.data(),
              edge_alpha_diff.size1(), &edge_alpha_diff(0, 0), &c1(0, 0));

          MatrixDouble edge_base_diff(M, edge_alpha_diff.size1());
          CHKERR BernsteinBezier::baseFunctionsEdge(
              N - 1, M, edge_alpha_diff.size1(), &edge_alpha_diff(0, 0),
              &edge_lambda(0, 0), Tools::diffShapeFunMBEDGE.data(),
              &edge_base_diff(0, 0), nullptr);

          for (size_t i = 0; i != M; ++i) {
            for (size_t j = 0; j != edge_diff_base.size2(); ++j) {

              double check_diff_base = 0;
              for (int k = 0; k != edge_alpha_diff.size1(); k++) {
                check_diff_base += c0(j, k) * edge_base_diff(i, k) *
                                   Tools::diffShapeFunMBEDGE[0];
                check_diff_base += c1(j, k) * edge_base_diff(i, k) *
                                   Tools::diffShapeFunMBEDGE[1];
              }
              const double error =
                  std::abs(check_diff_base - edge_diff_base(i, j));

              if (debug)
                std::cout << "edge_diff_base " << check_diff_base << " "
                          << edge_diff_base(i, j) << " " << error << std::endl;
              constexpr double eps = 1e-12;
              if (error > eps)
                SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                        "Property three on edge is not working");
            }

            if (debug)
              std::cout << endl;
          }

          MoFEMFunctionReturn(0);
        };

    CHKERR check_property_three_on_edge_derivatives(M, edge_alpha, edge_lambda,
                                                    edge_diff_base, false);

    // Triangle

    MatrixInt face_alpha(3 + 3 * NBEDGE_H1(N) + NBFACETRI_H1(N), 3);
    CHKERR BernsteinBezier::generateIndicesVertexTri(N, &face_alpha(0, 0));
    std::array<int, 3> face_edge_n{N, N, N};
    std::array<int *, 3> face_edge_ptr{&face_alpha(3, 0),
                                       &face_alpha(3 + NBEDGE_H1(N), 0),
                                       &face_alpha(3 + 2 * NBEDGE_H1(N), 0)};
    CHKERR BernsteinBezier::generateIndicesEdgeTri(face_edge_n.data(),
                                                   face_edge_ptr.data());
    CHKERR BernsteinBezier::generateIndicesTriTri(
        N, &face_alpha(3 + 3 * NBEDGE_H1(N), 0));
    // std::cout << "face alpha " << face_alpha << std::endl;

    std::array<double, 9> face_x_k = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    MatrixDouble face_x_alpha(face_alpha.size1(), 3);
    CHKERR BernsteinBezier::domainPoints3d(N, 3, face_alpha.size1(),
                                           &face_alpha(0, 0), face_x_k.data(),
                                           &face_x_alpha(0, 0));
    // std::cout << "domain points " << face_x_alpha << std::endl << std::endl;

    auto calc_lambda_on_face = [](auto &face_x_alpha) {
      MatrixDouble face_lambda(face_x_alpha.size1(), 3);
      for (size_t i = 0; i != face_x_alpha.size1(); ++i) {
        face_lambda(i, 0) = 1 - face_x_alpha(i, 0) - face_x_alpha(i, 1);
        face_lambda(i, 1) = face_x_alpha(i, 0);
        face_lambda(i, 2) = face_x_alpha(i, 1);
      }
      return face_lambda;
    };
    auto face_lambda = calc_lambda_on_face(face_x_alpha);

    MatrixDouble face_base(face_x_alpha.size1(), face_alpha.size1());
    MatrixDouble face_diff_base(face_x_alpha.size1(), 2 * face_alpha.size1());
    CHKERR BernsteinBezier::baseFunctionsTri(
        N, face_x_alpha.size1(), face_alpha.size1(), &face_alpha(0, 0),
        &face_lambda(0, 0), Tools::diffShapeFunMBTRI.data(), &face_base(0, 0),
        &face_diff_base(0, 0));
    // std::cout << "face base " << face_base << std::endl;

    CHKERR check_property_one(face_x_alpha.size1(), face_alpha, face_lambda,
                              face_base, BernsteinBezier::baseFunctionsTri,
                              false);

    auto check_property_two_on_face = [N](auto &face_alpha, bool debug) {
      MoFEMFunctionBegin;
      const int rule = N;
      const size_t nb_gauss_pts = QUAD_2D_TABLE[rule]->npoints;
      MatrixDouble gauss_pts(3, nb_gauss_pts, false);
      cblas_dcopy(nb_gauss_pts, &QUAD_2D_TABLE[rule]->points[1], 3,
                  &gauss_pts(0, 0), 1);
      cblas_dcopy(nb_gauss_pts, &QUAD_2D_TABLE[rule]->points[2], 3,
                  &gauss_pts(1, 0), 1);
      cblas_dcopy(nb_gauss_pts, QUAD_2D_TABLE[rule]->weights, 1,
                  &gauss_pts(2, 0), 1);
      MatrixDouble face_lambda(nb_gauss_pts, 3);
      cblas_dcopy(3 * nb_gauss_pts, QUAD_2D_TABLE[rule]->points, 1,
                  &face_lambda(0, 0), 1);
      MatrixDouble face_base(nb_gauss_pts, face_alpha.size1());
      CHKERR BernsteinBezier::baseFunctionsTri(
          N, nb_gauss_pts, face_alpha.size1(), &face_alpha(0, 0),
          &face_lambda(0, 0), Tools::diffShapeFunMBTRI.data(), &face_base(0, 0),
          nullptr);

      VectorDouble integral(face_alpha.size1());
      integral.clear();

      for (int g = 0; g != nb_gauss_pts; ++g) {
        for (size_t i = 0; i != face_alpha.size1(); ++i) {
          integral[i] += gauss_pts(2, g) * face_base(g, i) / 2.;
        }
      }

      const double check_integral =
          0.5 / boost::math::binomial_coefficient<double>(N + 2, 2);

      for (auto i : integral) {
        const double error = std::abs(i - check_integral);
        if (debug)
          std::cout << "edge integral " << i << " " << check_integral << " "
                    << error << endl;
        constexpr double eps = 1e-12;
        if (error > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "Property two on edge is not working");
      }

      MoFEMFunctionReturn(0);
    };
    CHKERR check_property_two_on_face(face_alpha, false);

    auto check_property_three_on_face_derivatives =
        [N](const int M, const auto &face_alpha, const auto &face_lambda,
            const auto &face_diff_base, bool debug) {
          MoFEMFunctionBegin;

          const int nb_edge = NBEDGE_H1(N - 1);
          const int nb_face = NBFACETRI_H1(N - 1);

          MatrixInt face_alpha_diff(3 + 3 * nb_edge + nb_face, 3);
          CHKERR BernsteinBezier::generateIndicesVertexTri(
              N - 1, &face_alpha_diff(0, 0));

          const std::array<int, 3> edge_n{N - 1, N - 1, N - 1};
          std::array<int *, 3> edge_ptr{&face_alpha_diff(3, 0),
                                        &face_alpha_diff(3 + nb_edge, 0),
                                        &face_alpha_diff(3 + 2 * nb_edge, 0)};

          CHKERR BernsteinBezier::generateIndicesEdgeTri(edge_n.data(),
                                                         edge_ptr.data());
          CHKERR BernsteinBezier::generateIndicesTriTri(
              N - 1, &face_alpha_diff(3 + 3 * nb_edge, 0));

          MatrixDouble c0(face_alpha.size1(), face_alpha_diff.size1());
          MatrixDouble c1(face_alpha.size1(), face_alpha_diff.size1());
          MatrixDouble c2(face_alpha.size1(), face_alpha_diff.size1());

          std::array<int, 3> diff0 = {1, 0, 0};
          CHKERR BernsteinBezier::genrateDerivativeIndicesTri(
              N, face_alpha.size1(), &face_alpha(0, 0), diff0.data(),
              face_alpha_diff.size1(), &face_alpha_diff(0, 0), &c0(0, 0));
          std::array<int, 3> diff1 = {0, 1, 0};
          CHKERR BernsteinBezier::genrateDerivativeIndicesTri(
              N, face_alpha.size1(), &face_alpha(0, 0), diff1.data(),
              face_alpha_diff.size1(), &face_alpha_diff(0, 0), &c1(0, 0));
          std::array<int, 3> diff2 = {0, 0, 1};
          CHKERR BernsteinBezier::genrateDerivativeIndicesTri(
              N, face_alpha.size1(), &face_alpha(0, 0), diff2.data(),
              face_alpha_diff.size1(), &face_alpha_diff(0, 0), &c2(0, 0));

          MatrixDouble face_base_diff(M, face_alpha_diff.size1());
          CHKERR BernsteinBezier::baseFunctionsTri(
              N - 1, M, face_alpha_diff.size1(), &face_alpha_diff(0, 0),
              &face_lambda(0, 0), Tools::diffShapeFunMBTRI.data(),
              &face_base_diff(0, 0), nullptr);

          for (size_t i = 0; i != M; ++i) {
            for (size_t j = 0; j != face_diff_base.size2() / 2; ++j) {

              VectorDouble3 check_diff_base(2);
              check_diff_base.clear();
              for (int k = 0; k != face_alpha_diff.size1(); ++k) {
                for (int d = 0; d != 2; ++d) {
                  check_diff_base[d] += c0(j, k) * face_base_diff(i, k) *
                                        Tools::diffShapeFunMBTRI[2 * 0 + d];
                  check_diff_base[d] += c1(j, k) * face_base_diff(i, k) *
                                        Tools::diffShapeFunMBTRI[2 * 1 + d];
                  check_diff_base[d] += c2(j, k) * face_base_diff(i, k) *
                                        Tools::diffShapeFunMBTRI[2 * 2 + d];
                }
              }
              auto diff_base = getVectorAdaptor(&face_diff_base(i, 2 * j), 2);
              const double error = norm_2(check_diff_base - diff_base);

              if (debug)
                std::cout << "face_diff_base " << check_diff_base << " "
                          << diff_base << " " << error << std::endl;
              constexpr double eps = 1e-12;
              if (error > eps)
                SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                        "Property three on face is not working");
            }

            if (debug)
              std::cout << endl;
          }

          MoFEMFunctionReturn(0);
        };

    CHKERR check_property_three_on_face_derivatives(
        face_x_alpha.size1(), face_alpha, face_lambda, face_diff_base, false);

    // Tetrahedron

    MatrixInt tet_alpha_vec(4, 4);
    CHKERR BernsteinBezier::generateIndicesVertexTet(N, &tet_alpha_vec(0, 0));

    const int nb_dofs_on_egde = NBEDGE_H1(N);
    std::array<MatrixInt, 6> tet_alpha_edge{
        MatrixInt(nb_dofs_on_egde, 4), MatrixInt(nb_dofs_on_egde, 4),
        MatrixInt(nb_dofs_on_egde, 4), MatrixInt(nb_dofs_on_egde, 4),
        MatrixInt(nb_dofs_on_egde, 4), MatrixInt(nb_dofs_on_egde, 4)};
    std::array<int, 6> tet_edge_n{N, N, N, N, N, N};
    std::array<int *, 6> tet_edge_ptr{
        &tet_alpha_edge[0](0, 0), &tet_alpha_edge[1](0, 0),
        &tet_alpha_edge[2](0, 0), &tet_alpha_edge[3](0, 0),
        &tet_alpha_edge[4](0, 0), &tet_alpha_edge[5](0, 0)};
    CHKERR BernsteinBezier::generateIndicesEdgeTet(tet_edge_n.data(),
                                                   tet_edge_ptr.data());

    const int nb_dofs_on_tri = NBFACETRI_H1(N);
    std::array<MatrixInt, 4> tet_alpha_face{
        MatrixInt(nb_dofs_on_tri, 4), MatrixInt(nb_dofs_on_tri, 4),
        MatrixInt(nb_dofs_on_tri, 4), MatrixInt(nb_dofs_on_tri, 4)};
    std::array<int, 4> tet_face_n{N, N, N, N};
    std::array<int *, 4> tet_face_ptr{
        &tet_alpha_face[0](0, 0), &tet_alpha_face[1](0, 0),
        &tet_alpha_face[2](0, 0), &tet_alpha_face[3](0, 0)};
    CHKERR BernsteinBezier::generateIndicesTriTet(tet_face_n.data(),
                                                  tet_face_ptr.data());

    const int nb_dofs_on_tet = NBVOLUMETET_H1(N);
    MatrixInt tet_alpha_tet(nb_dofs_on_tet, 4);
    CHKERR BernsteinBezier::generateIndicesTetTet(N, &tet_alpha_tet(0, 0));

    std::vector<MatrixInt *> alpha_tet;
    alpha_tet.push_back(&tet_alpha_vec);
    for (int e = 0; e != 6; ++e)
      alpha_tet.push_back(&tet_alpha_edge[e]);
    for (int f = 0; f != 4; ++f)
      alpha_tet.push_back(&tet_alpha_face[f]);
    alpha_tet.push_back(&tet_alpha_tet);

    auto create_tet_mesh = [&](auto &moab_ref, Range &tets) {
      MoFEMFunctionBegin;

      std::array<double, 12> base_coords{0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
      EntityHandle nodes[4];
      for (int nn = 0; nn < 4; nn++)
        CHKERR moab_ref.create_vertex(&base_coords[3 * nn], nodes[nn]);
      EntityHandle tet;
      CHKERR moab_ref.create_element(MBTET, nodes, 4, tet);

      MoFEM::CoreTmp<-1> m_core_ref(moab_ref, PETSC_COMM_SELF, -2);
      MoFEM::Interface &m_field_ref = m_core_ref;
      CHKERR m_field_ref.getInterface<BitRefManager>()->setBitRefLevelByDim(
          0, 3, BitRefLevel().set(0));
      const int max_level = 3;
      for (int ll = 0; ll != max_level; ll++) {
        Range edges;
        CHKERR m_field_ref.getInterface<BitRefManager>()
            ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(ll),
                                           BitRefLevel().set(), MBEDGE, edges);
        Range tets;
        CHKERR m_field_ref.getInterface<BitRefManager>()
            ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(ll),
                                           BitRefLevel(ll).set(), MBTET, tets);
        MeshRefinement *m_ref;
        CHKERR m_field_ref.getInterface(m_ref);
        CHKERR m_ref->add_vertices_in_the_middle_of_edges(
            edges, BitRefLevel().set(ll + 1));
        CHKERR m_ref->refine_TET(tets, BitRefLevel().set(ll + 1));
      }

      CHKERR m_field_ref.getInterface<BitRefManager>()
          ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(max_level),
                                         BitRefLevel().set(max_level), MBTET,
                                         tets);

      // Use 10 node tets to print base
      if (1) {
        EntityHandle meshset;
        CHKERR moab_ref.create_meshset(MESHSET_SET, meshset);
        CHKERR moab_ref.add_entities(meshset, tets);
        CHKERR moab_ref.convert_entities(meshset, true, false, false);
        CHKERR moab_ref.delete_entities(&meshset, 1);
      }

      MoFEMFunctionReturn(0);
    };

    Range tets;
    CHKERR create_tet_mesh(moab, tets);
    Range elem_nodes;
    CHKERR moab.get_connectivity(tets, elem_nodes, false);
    MatrixDouble coords(elem_nodes.size(), 3);
    CHKERR moab.get_coords(elem_nodes, &coords(0, 0));

    auto calc_lambda_on_tet = [](auto &coords) {
      MatrixDouble lambda(coords.size1(), 4);
      for (size_t i = 0; i != coords.size1(); ++i) {
        lambda(i, 0) = 1 - coords(i, 0) - coords(i, 1) - coords(i, 2);
        lambda(i, 1) = coords(i, 0);
        lambda(i, 2) = coords(i, 1);
        lambda(i, 3) = coords(i, 2);
      }
      return lambda;
    };
    auto lambda = calc_lambda_on_tet(coords);
    // cerr << lambda << endl;

    int nn = 0;
    for (auto alpha_ptr : alpha_tet) {

      MatrixDouble base(coords.size1(), alpha_ptr->size1());
      MatrixDouble diff_base(coords.size1(), 3 * alpha_ptr->size1());
      CHKERR BernsteinBezier::baseFunctionsTet(
          N, coords.size1(), alpha_ptr->size1(), &(*alpha_ptr)(0, 0),
          &lambda(0, 0), Tools::diffShapeFunMBTET.data(), &base(0, 0),
          &diff_base(0, 0));

      CHKERR check_property_one(coords.size1(), *alpha_ptr, lambda, base,
                                BernsteinBezier::baseFunctionsTet, false);

      // cerr << *alpha_ptr << endl;
      // cerr << base << endl;
      MatrixDouble trans_base = trans(base);
      for (int j = 0; j != base.size2(); ++j) {
        Tag th;
        double def_val[] = {0};
        CHKERR moab.tag_get_handle(
            ("base_" + boost::lexical_cast<std::string>(nn) + "_" +
             boost::lexical_cast<std::string>(j))
                .c_str(),
            1, MB_TYPE_DOUBLE, th, MB_TAG_CREAT | MB_TAG_DENSE, def_val);
        CHKERR moab.tag_set_data(th, elem_nodes, &trans_base(j, 0));
      }
      ++nn;
    }

    EntityHandle meshset;
    CHKERR moab.create_meshset(MESHSET_SET, meshset);
    CHKERR moab.add_entities(meshset, tets);
    CHKERR moab.write_file("bb.vtk", "VTK", "", &meshset, 1);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}