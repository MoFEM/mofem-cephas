/** \file BernsteinBezier.cpp
\brief Bernstein-Bezier polynomials for H1 space

Implementation based on \cite ainsworth2011bernstein

*/

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
MoFEMErrorCode BernsteinBezier::generateIndicesTriOnSimplex(const int N,
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
      alpha[tri_nodes[Side][2]] = N - n0 - n1;

      alpha += D + 1;
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::generateIndicesTetOnSimplex(const int N,
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

MoFEMErrorCode BernsteinBezier::generateIndicesVertexTri(const int N,
                                                         int *alpha) {
  MoFEMFunctionBeginHot;
  CHKERR generateIndicesVertex<2, 0>(N, alpha);
  CHKERR generateIndicesVertex<2, 1>(N, alpha);
  CHKERR generateIndicesVertex<2, 2>(N, alpha);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::generateIndicesEdgeTri(const int N[],
                                                       int *alpha[]) {
  MoFEMFunctionBeginHot;
  CHKERR generateIndicesEdgeOnSimplex<2, 0>(N[0], alpha[0]);
  CHKERR generateIndicesEdgeOnSimplex<2, 1>(N[1], alpha[1]);
  CHKERR generateIndicesEdgeOnSimplex<2, 2>(N[2], alpha[2]);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::generateIndicesEdgeTri(const int side,
                                                       const int N,
                                                       int *alpha) {
  switch (side) {
  case 0:
    return generateIndicesEdgeOnSimplex<2, 0>(N, alpha);
  case 1:
    return generateIndicesEdgeOnSimplex<2, 1>(N, alpha);
  case 2:
    return generateIndicesEdgeOnSimplex<2, 2>(N, alpha);
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Wrong side number, on triangle are edges from 0 to 2");
  }
}

MoFEMErrorCode BernsteinBezier::generateIndicesTriTri(const int N, int *alpha) {
  return generateIndicesTriOnSimplex<2, 3>(N, alpha);
}

MoFEMErrorCode BernsteinBezier::generateIndicesVertexTet(const int N,
                                                         int *alpha) {
  MoFEMFunctionBeginHot;
  CHKERR generateIndicesVertex<3, 0>(N, alpha);
  CHKERR generateIndicesVertex<3, 1>(N, alpha);
  CHKERR generateIndicesVertex<3, 2>(N, alpha);
  CHKERR generateIndicesVertex<3, 3>(N, alpha);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::generateIndicesEdgeTet(const int N[],
                                                       int *alpha[]) {
  MoFEMFunctionBeginHot;
  CHKERR generateIndicesEdgeOnSimplex<3, 0>(N[0], alpha[0]);
  CHKERR generateIndicesEdgeOnSimplex<3, 1>(N[1], alpha[1]);
  CHKERR generateIndicesEdgeOnSimplex<3, 2>(N[2], alpha[2]);
  CHKERR generateIndicesEdgeOnSimplex<3, 3>(N[3], alpha[3]);
  CHKERR generateIndicesEdgeOnSimplex<3, 4>(N[4], alpha[4]);
  CHKERR generateIndicesEdgeOnSimplex<3, 5>(N[5], alpha[5]);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::generateIndicesEdgeTet(const int side,
                                                       const int N,
                                                       int *alpha) {
  MoFEMFunctionBeginHot;
  switch (side) {
  case 0:
    return generateIndicesEdgeOnSimplex<3, 0>(N, alpha);
  case 1:
    return generateIndicesEdgeOnSimplex<3, 1>(N, alpha);
  case 2:
    return generateIndicesEdgeOnSimplex<3, 2>(N, alpha);
  case 3:
    return generateIndicesEdgeOnSimplex<3, 3>(N, alpha);
  case 4:
    return generateIndicesEdgeOnSimplex<3, 4>(N, alpha);
  case 5:
    return generateIndicesEdgeOnSimplex<3, 5>(N, alpha);
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Wrong side number, on tetrahedron are edges from 0 to 5");
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::generateIndicesTriTet(const int N[],
                                                      int *alpha[]) {
  MoFEMFunctionBeginHot;
  CHKERR generateIndicesTriOnSimplex<3, 0>(N[0], alpha[0]);
  CHKERR generateIndicesTriOnSimplex<3, 1>(N[1], alpha[1]);
  CHKERR generateIndicesTriOnSimplex<3, 2>(N[2], alpha[2]);
  CHKERR generateIndicesTriOnSimplex<3, 3>(N[3], alpha[3]);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::generateIndicesTriTet(const int side,
                                                      const int N, int *alpha) {
  MoFEMFunctionBeginHot;
  switch (side) {
  case 0:
    return generateIndicesTriOnSimplex<3, 0>(N, alpha);
  case 1:
    return generateIndicesTriOnSimplex<3, 1>(N, alpha);
  case 2:
    return generateIndicesTriOnSimplex<3, 2>(N, alpha);
  case 3:
    return generateIndicesTriOnSimplex<3, 3>(N, alpha);
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Wrong side number, on tetrahedron are from from 0 to 3");
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::generateIndicesTetTet(const int N, int *alpha) {
  return generateIndicesTetOnSimplex(N, alpha);
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

MoFEMErrorCode BernsteinBezier::genrateDerivativeIndicesTri(
    const int N, const int n_alpha, const int *alpha, const int *diff,
    const int n_alpha_diff, const int *alpha_diff, double *c) {
  return genrateDerivativeIndices<2>(N, n_alpha, alpha, diff, n_alpha_diff,
                                     alpha_diff, c);
}

MoFEMErrorCode BernsteinBezier::genrateDerivativeIndicesTet(
    const int N, const int n_alpha, const int *alpha, const int *diff,
    const int n_alpha_diff, const int *alpha_diff, double *c) {
  return genrateDerivativeIndices<3>(N, n_alpha, alpha, diff, n_alpha_diff,
                                     alpha_diff, c);
}

template <int D>
MoFEMErrorCode
BernsteinBezier::domainPoints(const int N, const int n_x, const int n_alpha,
                              const int *alpha, const double *x_k,
                              double *x_alpha) {
  FTensor::Index<'i', D> i;
  MoFEMFunctionBeginHot;

  auto t_x_alpha = FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
      x_alpha, &x_alpha[1], &x_alpha[2]);
  for (int n0 = 0; n0 != n_alpha; ++n0) {
    t_x_alpha(i) = 0;
    auto t_x_k = FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3>(
        x_k, &x_k[1], &x_k[2]);
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

template <int D, bool GRAD_BASE>
MoFEMErrorCode
BernsteinBezier::baseFunctions(const int N, const int gdim, const int n_alpha,
                               const int *alpha, const double *lambda,
                               const double *grad_lambda, double *base,
                               double *grad_base) {
  MoFEMFunctionBeginHot;

  const int *const alpha0 = alpha;
  const double fN = boost::math::factorial<double>(N);
  constexpr int MAX_ALPHA = 12;
  int max_alpha = *std::max_element(alpha, &alpha[(D + 1) * n_alpha]);
  if(max_alpha > MAX_ALPHA)
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Is assumed maximal order not to be bigger than %d", MAX_ALPHA);
  std::array<double, (D + 1) * (MAX_ALPHA + 1)> pow_alpha;
  std::array<double, (MAX_ALPHA + 1)> factorial_alpha;
  std::array<double, D + 1> terms, diff_terms;
  const double *const grad_lambda0 = grad_lambda;

  factorial_alpha[0] = 1;
  if (max_alpha >= 1)
    factorial_alpha[1] = 1;
  if (max_alpha >= 2)
    for (int a = 2; a <= max_alpha; ++a)
      factorial_alpha[a] = factorial_alpha[a - 1] * a;

  for (int g = 0; g != gdim; ++g) {

    for (int n = 0; n != D + 1; ++n, ++lambda) {
      const size_t shift = (MAX_ALPHA + 1) * n;
      double *pow_alpha_ptr = &pow_alpha[shift];
      *pow_alpha_ptr = 1;

      if (max_alpha >= 1) {
        ++pow_alpha_ptr;
        *pow_alpha_ptr = *lambda;
      }

      if (max_alpha >= 2)
        for (int a = 2; a <= max_alpha; ++a) {
          const double p = (*pow_alpha_ptr) * (*lambda);
          ++pow_alpha_ptr;
          *pow_alpha_ptr = p;
        }
    }

    for (int n0 = 0; n0 != n_alpha; ++n0) {

      grad_lambda = grad_lambda0;

      double f = factorial_alpha[(*alpha)];
      double *terms_ptr = terms.data();
      double *diff_terms_ptr = diff_terms.data();
      *terms_ptr = pow_alpha[(*alpha)];
      if (GRAD_BASE) {
        if (*alpha > 0)
          *diff_terms_ptr = (*alpha) * pow_alpha[(*alpha) - 1];
        else
          *diff_terms_ptr = 0;
      }
      *base = *terms_ptr;
      ++alpha;
      ++terms_ptr;
      ++diff_terms_ptr;

      for (int n1 = 1; n1 < D + 1;
           ++n1, ++alpha, ++terms_ptr, ++diff_terms_ptr) {
        f *= factorial_alpha[(*alpha)];
        const size_t shift = (MAX_ALPHA + 1) * n1;
        *terms_ptr = pow_alpha[shift + (*alpha)];
        if (GRAD_BASE) {
          if (*alpha > 0)
            *diff_terms_ptr = (*alpha) * pow_alpha[shift + (*alpha) - 1];
          else
            *diff_terms_ptr = 0;
        }
        *base *= *terms_ptr;
      }

      const double b = fN / f;
      *base *= b;
      ++base;

      if (GRAD_BASE) {
        double *terms_ptr = terms.data();
        double *diff_terms_ptr = diff_terms.data();
        double z = *diff_terms_ptr;
        ++terms_ptr;
        ++diff_terms_ptr;
        for (int n2 = 1; n2 != D + 1; ++n2, ++terms_ptr) 
          z *= *terms_ptr;
        
        for (int d = 0; d != D; ++d, ++grad_lambda) 
          grad_base[d] = z * (*grad_lambda);

        for (int n1 = 1; n1 < D + 1; ++n1) {
          z = *diff_terms_ptr;

          int n2 = 0;
          for (terms_ptr = terms.data(); n2 != n1; ++n2, ++terms_ptr)
            z *= *terms_ptr;

          ++n2;
          ++terms_ptr;
          
          for (; n2 < D + 1; ++n2, ++terms_ptr)
            z *= *terms_ptr;

          for (int d = 0; d != D; ++d, ++grad_lambda) 
            grad_base[d] += z * (*grad_lambda);

          ++diff_terms_ptr;
        }
        for (int d = 0; d != D; ++d)
          grad_base[d] *= b;
        grad_base += D;
      }

    }

    alpha = alpha0;
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BernsteinBezier::baseFunctionsEdge(
    const int N, const int gdim, const int n_alpha, const int *alpha,
    const double *lambda, const double *grad_lambda, double *base,
    double *grad_base) {
  if (grad_base)
    return baseFunctions<1, true>(N, gdim, n_alpha, alpha, lambda, grad_lambda,
                                  base, grad_base);
  else
    return baseFunctions<1, false>(N, gdim, n_alpha, alpha, lambda, grad_lambda,
                                   base, grad_base);
}

MoFEMErrorCode BernsteinBezier::baseFunctionsTri(
    const int N, const int gdim, const int n_alpha, const int *alpha,
    const double *lambda, const double *grad_lambda, double *base,
    double *grad_base) {
  if (grad_base)
    return baseFunctions<2, true>(N, gdim, n_alpha, alpha, lambda, grad_lambda,
                                  base, grad_base);
  else
    return baseFunctions<2, false>(N, gdim, n_alpha, alpha, lambda, grad_lambda,
                                   base, grad_base);
}

MoFEMErrorCode BernsteinBezier::baseFunctionsTet(
    const int N, const int gdim, const int n_alpha, const int *alpha,
    const double *lambda, const double *grad_lambda, double *base,
    double *grad_base) {
  if (grad_base)
    return baseFunctions<3, true>(N, gdim, n_alpha, alpha, lambda, grad_lambda,
                                  base, grad_base);
  else
    return baseFunctions<3, false>(N, gdim, n_alpha, alpha, lambda, grad_lambda,
                                   base, grad_base);
}
