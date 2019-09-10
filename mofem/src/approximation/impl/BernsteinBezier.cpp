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
  static_assert(

      !(D == 3 && S == 3),

      "not implemented");
  return FTensor::Tensor1<FTensor::PackPtr<double *, S>, D>();
}

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
BernsteinBezier::getFTensor1<3, 3>(double *x) {
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
MoFEMErrorCode BernsteinBezier::genrateDerivativeIndices(
    const int N, const int n_alpha, const int *alpha, const int *diff,
    const int n_alpha_diff, const int *alpha_diff, double *c) {
  MoFEMFunctionBeginHot;

  std::fill(c, &c[n_alpha * n_alpha_diff], 0);
  int abs_diff = diff[0];
  for (int d = 1; d != D + 1; ++d)
    abs_diff += diff[d];

  const double b = boost::math::factorial<double>(N) /
                   boost::math::factorial<double>(N - abs_diff);

  for (int i = 0; i != n_alpha; ++i) {
    for (int j = 0; j != n_alpha_diff; ++j) {
      int d = 0;
      for (; d != D + 1; ++d)
        if (alpha_diff[(D + 1) * j + d] + diff[d] != alpha[(D + 1) * i + d])
          break;
      if (d == D + 1) {
        c[i * n_alpha_diff + j] = b;
        break;
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::genrateDerivativeIndicesEdges(
    const int N, const int n_alpha, const int *alpha, const int *diff,
    const int n_alpha_diff, const int *alpha_diff, double *c) {
  return genrateDerivativeIndices<1>(N, n_alpha, alpha, diff, n_alpha_diff,
                                     alpha_diff, c);
}

template <int D>
MoFEMErrorCode
BernsteinBezier::domainPoints(const int N, const int n_x, const int n_alpha,
                              const int *alpha, const double *x_k,
                              double *x_alpha) {
  FTensor::Index<'i', D> i;
  MoFEMFunctionBeginHot;

  auto t_x_alpha = getFTensor1<D, D>(x_alpha);
  for (int n0 = 0; n0 != n_alpha; ++n0) {
    t_x_alpha(i) = 0;
    auto t_x_k = getFTensor1<D, D>(const_cast<double *>(x_k));
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
                                               const int n_alpha,
                                               const int *alpha,
                                               const double *x_k,
                                               double *x_alpha) {
  return domainPoints<3>(N, n_x, n_alpha, alpha, x_k, x_alpha);
}

template <int D>
MoFEMErrorCode
BernsteinBezier::baseFunctions(const int N, const int gdim, const int n_alpha,
                               const int *alpha, const double *lambda,
                               const double *grad_lambda, double *base,
                               double *grad_base) {
  MoFEMFunctionBeginHot;

  const int *const alpha0 = alpha;
  const double fN = boost::math::factorial<double>(N);
  std::array<double, D + 1> terms, diff_terms;
  const double *const grad_lambda0 = grad_lambda;

  for (int g = 0; g != gdim; ++g) {

    const double *const lambda0 = lambda;

    for (int n0 = 0; n0 != n_alpha; ++n0) {

      grad_lambda = grad_lambda0;

      double f = boost::math::factorial<double>(*alpha);
      terms[0] = pow(*lambda, (*alpha));
      if (*alpha > 0)
        diff_terms[0] = (*alpha) * pow(*lambda, (*alpha) - 1);
      else
        diff_terms[0] = 0;
      *base = terms[0];
      ++alpha;
      ++lambda;

      for (int n1 = 1; n1 < D + 1; ++n1) {
        f *= boost::math::factorial<double>(*alpha);
        terms[n1] = pow(*lambda, (*alpha));
        if (*alpha > 0)
          diff_terms[n1] = (*alpha) * pow(*lambda, (*alpha) - 1);
        else
          diff_terms[n1] = 0;
        *base *= terms[n1];
        ++alpha;
        ++lambda;
      }

      double z = diff_terms[0];
      for (int n2 = 1; n2 != D + 1; ++n2)
        z *= terms[n2];
      for (int d = 0; d != D; ++d) {
        grad_base[d] = z * (*grad_lambda);
        ++grad_lambda;
      }

      for (int n1 = 1; n1 < D + 1; ++n1) {
        z = diff_terms[n1];
        int n2 = 0;
        for (; n2 != n1; ++n2)
          z *= terms[n2];
        n2++;
        for (; n2 < D + 1; ++n2)
          z *= terms[n2];
        for (int d = 0; d != D; ++d) {
          grad_base[d] += z * (*grad_lambda);
          ++grad_lambda;
        }
      }

      const double b = fN / f;
      *base *= b;
      for (int d = 0; d != D; ++d)
        grad_base[d] *= b;

      ++base;
      grad_base += D;

      lambda = lambda0;
    }

    alpha = alpha0;
    lambda += D + 1;
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::baseFunctionsEdge(
    const int N, const int gdim, const int n_alpha, const int *alpha,
    const double *lambda, const double *grad_lambda, double *base,
    double *grad_base) {
  return baseFunctions<1>(N, gdim, n_alpha, alpha, lambda, grad_lambda, base,
                          grad_base);
}
