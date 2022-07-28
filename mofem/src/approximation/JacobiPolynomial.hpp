/** \file JacobiPolynomial.hpp
\brief Implementation of Legendre polynomial

*/



#ifndef __JACOBIPOLYNOMIALS_HPP__
#define __JACOBIPOLYNOMIALS_HPP__

namespace MoFEM {

/**
 * \brief Class used to give arguments to Legendre base functions
 * \ingroup mofem_base_functions
 */
struct JacobiPolynomialCtx : public BaseFunctionCtx {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  int P;
  double *diffX;
  double *diffT;
  int dIm;

  double aLpha;

  boost::shared_ptr<MatrixDouble> baseFunPtr;
  boost::shared_ptr<MatrixDouble> baseDiffFunPtr;

  PetscErrorCode (*basePolynomialsType1)(int p, double alpha, double x,
                                         double t, double *diff_x,
                                         double *diff_t, double *L,
                                         double *diffL, const int dim);

  JacobiPolynomialCtx(int p, double *diff_x, double *diff_t, int dim,
                      double alpha,
                      boost::shared_ptr<MatrixDouble> &base_fun_ptr,
                      boost::shared_ptr<MatrixDouble> &base_diff_fun_ptr)
      : P(p), diffX(diff_x), diffT(diff_t), dIm(dim), aLpha(alpha),
        baseFunPtr(base_fun_ptr), baseDiffFunPtr(base_diff_fun_ptr),
        basePolynomialsType1(Jacobi_polynomials) {}
  ~JacobiPolynomialCtx() {}
};

/**
 * \brief Calculating Legendre base functions
 * \ingroup mofem_base_functions
 */
struct JacobiPolynomial : public BaseFunction {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  JacobiPolynomial() {}
  ~JacobiPolynomial() {}

  MoFEMErrorCode getValue(MatrixDouble &pts_x, MatrixDouble &pts_t,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);
};

struct IntegratedJacobiPolynomialCtx : public JacobiPolynomialCtx {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  IntegratedJacobiPolynomialCtx(
      int p, double *diff_x, double *diff_t, int dim, double alpha,
      boost::shared_ptr<MatrixDouble> &base_fun_ptr,
      boost::shared_ptr<MatrixDouble> &base_diff_fun_ptr)
      : JacobiPolynomialCtx(p, diff_x, diff_t, dim, alpha, base_fun_ptr,
                            base_diff_fun_ptr) {
    basePolynomialsType1 = IntegratedJacobi_polynomials;
  }
  ~IntegratedJacobiPolynomialCtx() {}
};

struct IntegratedJacobiPolynomial : public JacobiPolynomial {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  IntegratedJacobiPolynomial() : JacobiPolynomial() {}
  ~IntegratedJacobiPolynomial() {}

  MoFEMErrorCode getValue(MatrixDouble &pts_x, MatrixDouble &pts_t,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);
};

} // namespace MoFEM

#endif //__JACOBIPOLYNOMIALS_HPP__
