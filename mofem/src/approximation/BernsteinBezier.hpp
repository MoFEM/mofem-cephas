/** \file BernsteinBezier.hpp
\brief Bernstein-Bezier polynomials for H1 space

*/

#ifndef _H1_BERNSTEIN_BEZIER__HPP_
#define _H1_BERNSTEIN_BEZIER__HPP_

namespace MoFEM {

struct BernsteinBezier {

private:
  template <int D, int S>
  static inline FTensor::Tensor1<FTensor::PackPtr<double *, S>, D>
  getFTensor1(double *x);

  template <int D, int Side>
  static MoFEMErrorCode generateIndicesVertex(const int N, int *alpha);

  template <int D, int Side>
  static MoFEMErrorCode generateIndicesEdgeOnSimplex(const int N, int *alpha);

  template <int D, int Side>
  static MoFEMErrorCode generateIndicesFaceOnSimplex(const int N, int *alpha);

  static MoFEMErrorCode generateIndicesVolumeOnSimplex(const int N, int *alpha);

  template <int D>
  static MoFEMErrorCode domainPoints(const int N, int *alpha, double *x_k,
                                     double *x_alpha);

  template <int D>
  static MoFEMErrorCode baseFunctions(const int N, int *alpha, double *lambda,
                                     double *grad_lambda, double *base,
                                     double *grad_base);
};

} // namespace MoFEM

#endif // _H1_BERNSTEIN_BEZIER__HPP_