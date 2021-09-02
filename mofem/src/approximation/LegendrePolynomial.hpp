/** \file LegendrePolynomial.hpp
\brief Implementation of Legendre polynomial

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
