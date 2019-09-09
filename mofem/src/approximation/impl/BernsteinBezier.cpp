/** \file BernsteinBezier.cpp
\brief Bernstein-Bezier polynomials for H1 space

Implementation based on \cite ainsworth2011bernstein

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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>.
 */

template <int D, int S>
FTensor::Tensor1<FTensor::PackPtr<double *, S>, D>
BernsteinBezier::getFTensor1(double *x) {
  static_assert((D != 3 && S != 3), "not implemented");
  return FTensor::Tensor1<FTensor::PackPtr<double *, S>, D>();
}

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
BernsteinBezier::getFTensor1(double *x) {
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(x, &x[1], &x[2]);
}

template <int D, int Side>
MoFEMErrorCode BernsteinBezier::generateIndicesVertex(const int N, int *alpha) {
  MoFEMFunctionBeginHot;
  alpha[(D + 1) * Side + Side] = N;
  MoFEMFunctionReturnHot(0);
}

template <int D, int Side>
MoFEMErrorCode BernsteinBezier::generateIndicesEdgeOnSimplex(const int N,
                                                             int *alpha) {
  constexpr int edge_nodes[6][2] = {{0, 1}, {1, 2}, {2, 0},
                                    {0, 3}, {1, 3}, {2, 3}};
  MoFEMFunctionBeginHot;
  std::fill(alpha, &alpha[(D + 1) * NBEDGE_H1(N)], 0);
  for (int n = N - 1; n != 0; --n) {

    // + o o o o +

    alpha[edge_nodes[Side][0]] = n;
    alpha[edge_nodes[Side][1]] = N - n;
    alpha += D + 1;
  }
  MoFEMFunctionReturnHot(0);
}

template <int D, int Side>
MoFEMErrorCode BernsteinBezier::generateIndicesFaceOnSimplex(const int N,
                                                             int *alpha) {
  constexpr int tri_nodes[4][3] = {{0, 1, 3}, {1, 2, 3}, {0, 2, 3}, {0, 1, 2}};
  MoFEMFunctionBeginHot;
  std::fill(alpha, &alpha[(D + 1) * NBFACETRI_H1(N)], 0);
  for (int n0 = 1; n0 != N - 1; ++n0) {
    for (int n1 = 1; n1 != N - n0; ++n1) {

      // +
      // + +
      // + o +
      // + o o +
      // + o o o +
      // + + + + + +

      alpha[tri_nodes[Side][0]] = n0;
      alpha[tri_nodes[Side][1]] = n1;
      alpha[tri_nodes[Side][1]] = N - n0 - n1;

      alpha += D + 1;
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::generateIndicesVolumeOnSimplex(const int N,
                                                               int *alpha) {
  MoFEMFunctionBeginHot;
  std::fill(alpha, &alpha[4 * NBVOLUMETET_H1(N)], 0);
  for (int n0 = 1; n0 != N - 1; ++n0) {
    for (int n1 = 1; n1 != N - n0; ++n1) {
      for (int n2 = 1; n2 != N - n0 - n1; ++n2) {
        alpha[0] = n0;
        alpha[1] = n1;
        alpha[2] = n2;
        alpha[3] = N - n0 - n1 - n2;
        alpha += 4;
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::generateIndicesVertexEdge(const int N,
                                                          int *alpha) {
  MoFEMFunctionBeginHot;
  CHKERR generateIndicesVertex<1, 0>(N, alpha);
  CHKERR generateIndicesVertex<1, 1>(N, alpha);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::generateIndicesEdgeEdge(const int N,
                                                        int *alpha) {
  return generateIndicesEdgeOnSimplex<1, 0>(N, alpha);
}

template <int D>
MoFEMErrorCode BernsteinBezier::domainPoints(const int N, const int n_x,
                                             const int n_alpha, int *alpha,
                                             double *x_k, double *x_alpha) {
  FTensor::Index<'i', D> i;
  MoFEMFunctionBeginHot;

  auto t_x_alpha = getFTensor1<D, D>(x_alpha);
  for (int n0 = 0; n0 != n_alpha; ++n0) {
    t_x_alpha(i) = 0;
    auto t_x_k = getFTensor1<D, D>(x_k);
    for (int n1 = 0; n1 != n_x; ++n1) {
      t_x_alpha(i) += static_cast<double>(*alpha) * t_x_k(i);
      ++t_x_k;
      ++alpha;
    }
    t_x_alpha(i) /= static_cast<double>(N);
    ++t_x_alpha;
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::domainPoints3d(const int N, const int n_x,
                                               const int n_alpha, int *alpha,
                                               double *x_k, double *x_alpha) {
  return domainPoints<3>(N, n_x, n_alpha, alpha, x_k, x_alpha);
}

// template <int D>
// MoFEMErrorCode BernsteinBezier::baseFunctions(const int N, int *alpha,
//                                               double *lambda,
//                                               double *grad_lambda, double
//                                               *base, double *grad_base) {
//   FTensor::Index<'i', D> i;
//   MoFEMFunctionBeginHot;

//   auto t_base_grad = getFTensor1<D, D>(grad_base);

//   for (int n = 0; n != N; ++n) {
//     auto t_lambda_grad = getFTensor1<D, D>(grad_lambda);
//     for (int d = 0; d != D + 1; ++d) {
//       const double b = boost::math::binomial_coefficient<double>(N, *alpha);
//       const double a = b * pow((*lambda), (*alpha) - 1);
//       *base = a * (*lambda);
//       t_base_grad(i) = (*alpha) * a * t_lambda_grad(i);
//       ++lambda;
//       ++alpha;
//       ++t_lambda_grad;
//       ++base;
//       ++t_base_grad;
//     }
//   }

//   MoFEMFunctionReturnHot(0);
// }

// MoFEMErrorCode BernsteinBezier::nodeDomainPointsOnEdge3d(const int N,
//                                                          double *x_k,
//                                                          double *x_alpha) {
//   MoFEMFunctionBegin;
//   std::array<int, 2> alpha;
//   CHKERR generateIndicesVertex<1>(2, alpha.data());
//   CHKERR domainPoints<3>(N, 2, 2, alpha.data(), x_k, x_alpha);
//   MoFEMFunctionReturn(0);
// }

// MoFEMErrorCode BernsteinBezier::nodeBaseFunctionsOnEdge(const int N,
//                                                         double *lambda,
//                                                         double *grad_lambda,
//                                                         double *base,
//                                                         double *grad_base) {
//   MoFEMFunctionBegin;
//   std::array<int, 2> alpha;
//   CHKERR generateIndicesVertex<1>(N, alpha.data());
//   CHKERR baseFunctions<1>(N, alpha.data(), lambda, grad_lambda, base,
//                           grad_base);
//   MoFEMFunctionReturn(0);
// }

// MoFEMErrorCode BernsteinBezier::nodeBaseFunctionsOnTriangle(const int N,
//                                                             double *lambda,
//                                                             double
//                                                             *grad_lambda,
//                                                             double *base,
//                                                             double
//                                                             *grad_base) {
//   MoFEMFunctionBegin;
//   std::array<int, 3> alpha;
//   CHKERR generateIndicesVertex<2>(N, alpha.data());
//   CHKERR baseFunctions<2>(N, alpha.data(), lambda, grad_lambda, base,
//                           grad_base);
//   MoFEMFunctionReturn(0);
// }

// MoFEMErrorCode BernsteinBezier::nodeBaseFunctionsOnTetrahedron(
//     const int N, double *lambda, double *grad_lambda, double *base,
//     double *grad_base) {
//   MoFEMFunctionBegin;
//   std::array<int, 4> alpha;
//   CHKERR generateIndicesVertex<3>(N, alpha.data());
//   CHKERR baseFunctions<3>(N, alpha.data(), lambda, grad_lambda, base,
//                           grad_base);
//   MoFEMFunctionReturn(0);
// }

// MoFEMErrorCode BernsteinBezier::edgeBaseFunctionsOnEdge(const int N,
//                                                         double *lambda,
//                                                         double *grad_lambda,
//                                                         double *base,
//                                                         double *grad_base) {
//   MoFEMFunctionBegin;
//   VectorInt alpha(NBEDGE_H1(N));
//   CHKERR generateIndicesEdgeOnSimplex<1, 0>(N, &*alpha.data().begin());
//   CHKERR baseFunctions<1>(N, &*alpha.data().begin(), lambda, grad_lambda,
//   base,
//                           grad_base);
//   MoFEMFunctionReturn(0);
// }