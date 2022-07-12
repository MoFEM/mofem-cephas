/** \file JacobiPolynomial.hpp
\brief Implementation of Legendre polynomial

*/

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
