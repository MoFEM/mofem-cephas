/** \file LegendrePolynomial.hpp
\brief Implementation of Legendre polynomial

*/



#ifndef __LEGENDREPOLYNOMIALS_HPP__
#define __LEGENDREPOLYNOMIALS_HPP__

namespace MoFEM {

/**
 * \brief Class used to give arguments to Legendre base functions
 * \ingroup mofem_base_functions
 */
struct LegendrePolynomialCtx : public BaseFunctionCtx {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  int P;
  double *diffS;
  int dIm;
  boost::shared_ptr<MatrixDouble> baseFunPtr;
  boost::shared_ptr<MatrixDouble> baseDiffFunPtr;

  PetscErrorCode (*basePolynomialsType0)(int p, double s, double *diff_s,
                                         double *L, double *diffL,
                                         const int dim);

  LegendrePolynomialCtx(int p, double *diff_s, int dim,
                        boost::shared_ptr<MatrixDouble> base_fun_ptr,
                        boost::shared_ptr<MatrixDouble> base_diff_fun_ptr)
      : P(p), diffS(diff_s), dIm(dim), baseFunPtr(base_fun_ptr),
        baseDiffFunPtr(base_diff_fun_ptr),
        basePolynomialsType0(Legendre_polynomials) {}
  ~LegendrePolynomialCtx() {}
};

/**
 * \brief Calculating Legendre base functions
 * \ingroup mofem_base_functions
 */
struct LegendrePolynomial : public BaseFunction {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  LegendrePolynomial() {}
  ~LegendrePolynomial() {}

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);
};

} // namespace MoFEM

#endif //__LEGENDREPOLYNOMIALS_HPP__
