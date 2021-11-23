/** \file LobattoPolynomial.hpp
\brief Implementation of Lobatto polynomial

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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

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
