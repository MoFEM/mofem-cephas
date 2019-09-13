/** \file BernsteinBezier.hpp
\brief Bernstein-Bezier polynomials for H1 space

Implementation based on \cite ainsworth2011bernstein

*/

#ifndef _H1_BERNSTEIN_BEZIER__HPP_
#define _H1_BERNSTEIN_BEZIER__HPP_

namespace MoFEM {

/**
 * @brief Evaluating BB polynomial
 *
 * Low level class for evaluation of Bernstein-Bezier polynomials and associated
 * tools useful for fast intergartion.
 *
 */
struct BernsteinBezier {

  /** \name Edge BB functions */

  /**@{*/

  static MoFEMErrorCode generateIndicesVertexEdge(const int N, int *alpha);
  static MoFEMErrorCode generateIndicesEdgeEdge(const int N, int *alpha);
  static MoFEMErrorCode baseFunctionsEdge(const int N, const int gdim,
                                          const int n_alpha, const int *alpha,
                                          const double *lambda,
                                          const double *grad_lambda,
                                          double *base, double *grad_base);
  static MoFEMErrorCode genrateDerivativeIndicesEdges(
      const int N, const int n_alpha, const int *alpha, const int *diff,
      const int n_alpha_diff, const int *alpha_diff, double *c);

  /**@}*/

  /** \name Triangle BB functions */

  /**@{*/

  static MoFEMErrorCode generateIndicesVertexTri(const int N, int *alpha);
  static MoFEMErrorCode generateIndicesEdgeTri(const int N[], int *alpha[]);
  static MoFEMErrorCode generateIndicesTriTri(const int N, int *alpha);
  static MoFEMErrorCode baseFunctionsTri(const int N, const int gdim,
                                         const int n_alpha, const int *alpha,
                                         const double *lambda,
                                         const double *grad_lambda,
                                         double *base, double *grad_base);
  static MoFEMErrorCode
  genrateDerivativeIndicesTri(const int N, const int n_alpha, const int *alpha,
                              const int *diff, const int n_alpha_diff,
                              const int *alpha_diff, double *c);
  
  /**@}*/

  /** \name Tetrahedron BB functions */

  static MoFEMErrorCode generateIndicesVertexTet(const int N, int *alpha);
  static MoFEMErrorCode generateIndicesEdgeTet(const int N[], int *alpha[]);
  static MoFEMErrorCode generateIndicesTriTet(const int N[], int *alpha[]);
  static MoFEMErrorCode generateIndicesTetTet(const int N, int *alpha);
  static MoFEMErrorCode baseFunctionsTet(const int N, const int gdim,
                                         const int n_alpha, const int *alpha,
                                         const double *lambda,
                                         const double *grad_lambda,
                                         double *base, double *grad_base);
  static MoFEMErrorCode
  genrateDerivativeIndicesTet(const int N, const int n_alpha, const int *alpha,
                              const int *diff, const int n_alpha_diff,
                              const int *alpha_diff, double *c);

  /**@}*/
 
  /**
   * @brief Genrate BB points in 3d
   * 
   * @param N 
   * @param n_x 
   * @param n_alpha 
   * @param alpha 
   * @param x_k 
   * @param x_alpha 
   * @return MoFEMErrorCode 
   */
  static MoFEMErrorCode domainPoints3d(const int N, const int n_x,
                                       const int n_alpha, const int *alpha,
                                       const double *x_k, double *x_alpha);


  private:
  template <int D, int S>
  static inline FTensor::Tensor1<FTensor::PackPtr<double *, S>, D>
  getFTensor1(double *x);

  /**
   * @brief Generate BB indices on vertices
   *
   * @tparam D
   * @tparam Side
   * @param N
   * @param alpha
   * @return MoFEMErrorCode
   */
  template <int D, int Side>
  inline static MoFEMErrorCode generateIndicesVertex(const int N, int *alpha);

  /**
   * @brief Genarte BB incices od simplex edges
   *
   * @tparam D
   * @tparam Side
   * @param N
   * @param alpha
   * @return MoFEMErrorCode
   */
  template <int D, int Side>
  inline static MoFEMErrorCode generateIndicesEdgeOnSimplex(const int N,
                                                            int *alpha);
  /**
   * @brief Generate BB indices on simples triangles
   *
   * @tparam D
   * @tparam Side
   * @param N
   * @param alpha
   * @return MoFEMErrorCode
   */
  template <int D, int Side>
  inline static MoFEMErrorCode generateIndicesTriOnSimplex(const int N,
                                                           int *alpha);
  /**
   * @brief Genarte BB indices on simplex


   \f[
   \frac{\partial^{|\mathbf{v}|} \phi^n_{\pmb\alpha}}{\partial
   \pmb\lambda^\mathbf{v}}= \mathbf{c}_{\pmb\alpha\pmb\beta}
   \phi^{n-|\mathbf{v}|}_{\pmb\beta} 
   \f]

   * @param N
   * @param alpha
   * @return MoFEMErrorCode
   */
  static MoFEMErrorCode generateIndicesTetOnSimplex(const int N, int *alpha);

  /**
   * @brief Brief calculate coefficients for directive of base functions
   *
   * \f[
   *
   * \f]
   *
   * @tparam D
   * @param N
   * @param n_alpha
   * @param alpha
   * @param diff
   * @param n_alpha_diff
   * @param alpha_diff
   * @param c
   * @return MoFEMErrorCode
   */
  template <int D>
  static MoFEMErrorCode
  genrateDerivativeIndices(const int N, const int n_alpha, const int *alpha,
                           const int *diff, const int n_alpha_diff,
                           const int *alpha_diff, double *c);

  /**
   * @brief Evaluate coordinates of BB points
   *
   * @tparam D dimension
   * @param N polynomail order
   * @param n_x
   * @param n_alpha
   * @param alpha
   * @param x_k cartesian coordinates of simplex nodes
   * @param x_alpha BB points
   * @return MoFEMErrorCode
   */
  template <int D>
  inline static MoFEMErrorCode domainPoints(const int N, const int n_x,
                                            const int n_alpha, const int *alpha,
                                            const double *x_k, double *x_alpha);

  /**
   * @brief BB base function
   *
   * @tparam D dimension
   * @tparam GRAD_BASE true to calculate base
   * @param N polynomial order
   * @param gdim number of evaluation points
   * @param n_alpha number BB base functions
   * @param alpha BB coefficients
   * @param lambda barycentric coefficients at evaluating points
   * @param grad_lambda gradients of barycentric over spatial coordinates
   * @param base base function
   * @param grad_base vector base functions
   * @return MoFEMErrorCode
   */
  template <int D, bool GRAD_BASE>
  inline static MoFEMErrorCode
  baseFunctions(const int N, const int gdim, const int n_alpha,
                const int *alpha, const double *lambda,
                const double *grad_lambda, double *base, double *grad_base);
};

} // namespace MoFEM

#endif // _H1_BERNSTEIN_BEZIER__HPP_