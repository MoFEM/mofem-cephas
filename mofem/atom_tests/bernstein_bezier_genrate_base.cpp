/**
 * \file bernstein_bezier_genrate_base.cpp
 * \example bernstein_bezier_genrate_base.cpp
 *
 * Genarte Bernstein-Bezier base
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

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
    // Create MoFEM instance
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    const int N = 6;

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
    cerr << edge_diff_base << endl;

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

    auto check_property_one = [&](const int M, auto &edge_alpha,
                                 auto &edge_base) {
      MoFEMFunctionBegin;

      MatrixInt edge_alpha2(diag_n(edge_alpha.size1()), 2);
      int k = 0;
      for (int i = 0; i != edge_alpha.size1(); ++i) {
        for (int j = i; j != edge_alpha.size1(); ++j, ++k) {
          for (int d = 0; d != edge_alpha.size2(); ++d) {
            edge_alpha2(k, d) = edge_alpha(i, d) + edge_alpha(j, d);
          }
        }
      }

      MatrixDouble edge_base2(M, edge_alpha2.size1());
      MatrixDouble edge_diff_base2(M, edge_alpha2.size1());
      CHKERR BernsteinBezier::baseFunctionsEdge(
          M + M, M, edge_alpha2.size1(), &edge_alpha2(0, 0), &edge_lambda(0, 0),
          Tools::diffShapeFunMBEDGE.data(), &edge_base2(0, 0),
          &edge_diff_base2(0, 0));

      const double f0 = boost::math::binomial_coefficient<double>(N + N, N);
      for (int g = 0; g != M; ++g) {
        int k = 0;
        for (size_t i = 0; i != edge_alpha.size1(); ++i) {
          for (size_t j = i; j != edge_alpha.size1(); ++j, ++k) {
            const double f =
                binomial_alpha_beta(2, &edge_alpha(i, 0), &edge_alpha(j, 0));
            const double b = f / f0;
            const double B_check = b * edge_base2(g, k);
            const double B_mult = edge_base(g, i) * edge_base(g, j);
            const double error = std::abs(B_check - B_mult);
            std::cout << "( " << k << " " << i << " " << j << " )" << B_check
                      << " " << B_mult << " " << error << std::endl;
          }
        }
      }
      MoFEMFunctionReturn(0);
    };

    // CHKERR check_property_one(M, edge_alpha, edge_base);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}