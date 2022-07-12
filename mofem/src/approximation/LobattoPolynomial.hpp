/** \file LobattoPolynomial.hpp
\brief Implementation of Lobatto polynomial

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

#ifndef __LOBATTOPOLYNOMIALS_HPP__
#define __LOBATTOPOLYNOMIALS_HPP__

namespace MoFEM {

/**
 * \brief Class used to give arguments to Lobatto base functions
 * \ingroup mofem_base_functions
 */
struct LobattoPolynomialCtx : public LegendrePolynomialCtx {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  LobattoPolynomialCtx(int p, double *diff_s, int dim,
                       boost::shared_ptr<MatrixDouble> base_fun_ptr,
                       boost::shared_ptr<MatrixDouble> base_diff_fun_ptr)
      : LegendrePolynomialCtx(p, diff_s, dim, base_fun_ptr, base_diff_fun_ptr) {
    basePolynomialsType0 = Lobatto_polynomials;
  }
  ~LobattoPolynomialCtx() {}
};

/**
 * \brief Calculating Lobatto base functions
 * \ingroup mofem_base_functions
 */
struct LobattoPolynomial : public LegendrePolynomial {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  LobattoPolynomial() {}
  ~LobattoPolynomial() {}

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);
};

/**
 * \brief Class used to give arguments to Kernel Lobatto base functions
 * \ingroup mofem_base_functions
 */
struct KernelLobattoPolynomialCtx : public LegendrePolynomialCtx {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  KernelLobattoPolynomialCtx(int p, double *diff_s, int dim,
                             boost::shared_ptr<MatrixDouble> base_fun_ptr,
                             boost::shared_ptr<MatrixDouble> base_diff_fun_ptr)
      : LegendrePolynomialCtx(p, diff_s, dim, base_fun_ptr, base_diff_fun_ptr) {
    basePolynomialsType0 = LobattoKernel_polynomials;
  }
  ~KernelLobattoPolynomialCtx() {}
};

/**
 * \brief Calculating Lobatto base functions
 * \ingroup mofem_base_functions
 */
struct KernelLobattoPolynomial : public LegendrePolynomial {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  KernelLobattoPolynomial() {}
  ~KernelLobattoPolynomial() {}

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);
};

} // namespace MoFEM

#endif
