/** \file BernsteinBezier.hpp
\brief Bernstein-Bezier polynomials for H1 space

Implementation based on \cite ainsworth2011bernstein

*/

#ifndef _H1_BERNSTEIN_BEZIER__HPP_
#define _H1_BERNSTEIN_BEZIER__HPP_

namespace MoFEM {

struct BernsteinBezier {

  static MoFEMErrorCode generateIndicesVertexEdge(const int N, int *alpha);
  static MoFEMErrorCode generateIndicesEdgeEdge(const int N, int *alpha);

  static MoFEMErrorCode domainPoints3d(const int N, const int n_x,
                                       const int n_alpha, int *alpha,
                                       double *x_k, double *x_alpha);

  // static MoFEMErrorCode nodeDomainPointsOnEdge3d(const int N, double *x_k,
  //                                                double *x_alpha);

  // static MoFEMErrorCode nodeBaseFunctionsOnEdge(const int N, double *lambda,
  //                                               double *grad_lambda,
  //                                               double *base,
  //                                               double *grad_base);

  // static MoFEMErrorCode nodeBaseFunctionsOnTriangle(const int N, double *lambda,
  //                                                   double *grad_lambda,
  //                                                   double *base,
  //                                                   double *grad_base);

  // static MoFEMErrorCode nodeBaseFunctionsOnTetrahedron(const int N,
  //                                                      double *lambda,
  //                                                      double *grad_lambda,
  //                                                      double *base,
  //                                                      double *grad_base);

  // static MoFEMErrorCode edgeBaseFunctionsOnEdge(const int N, double *lambda,
  //                                               double *grad_lambda,
  //                                               double *base,
  //                                               double *grad_base);


private:
  template <int D, int S>
  static inline FTensor::Tensor1<FTensor::PackPtr<double *, S>, D>
  getFTensor1(double *x);

  template <int D, int Side>
  inline static MoFEMErrorCode generateIndicesVertex(const int N, int *alpha);

  template <int D, int Side>
  inline static MoFEMErrorCode generateIndicesEdgeOnSimplex(const int N,
                                                            int *alpha);

  template <int D, int Side>
  inline static MoFEMErrorCode generateIndicesFaceOnSimplex(const int N,
                                                            int *alpha);

  static MoFEMErrorCode generateIndicesVolumeOnSimplex(const int N, int *alpha);

  template <int D>
  inline static MoFEMErrorCode domainPoints(const int N, const int n_x,
                                            const int n_alpha, int *alpha,
                                            double *x_k, double *x_alpha);

  // template <int D>
  // inline static MoFEMErrorCode
  // baseFunctions(const int N, int *alpha, double *lambda, double *grad_lambda,
  //               double *base, double *grad_base);
};

} // namespace MoFEM

#endif // _H1_BERNSTEIN_BEZIER__HPP_