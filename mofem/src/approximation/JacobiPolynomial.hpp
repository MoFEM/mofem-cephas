/** \file JacobiPolynomial.hpp
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

#ifndef __JACOBIPOLYNOMIALS_HPP__
#define __JACOBIPOLYNOMIALS_HPP__

namespace MoFEM {

  static const MOFEMuuid IDD_JACOBI_BASE_FUNCTION =
  MOFEMuuid(BitIntefaceId(JACOBI_BASE_FUNCTION_INTERFACE));

  /**
   * \brief Class used to give arguments to Legendre base functions
   * \ingroup mofem_base_functions
   */
  struct JacobiPolynomialCtx: public BaseFunctionCtx {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    int P;
    double *diffX;
    double *diffT;
    int dIm;

    double aLpha;

    boost::shared_ptr<ublas::matrix<double> > baseFunPtr;
    boost::shared_ptr<ublas::matrix<double> > baseDiffFunPtr;

    PetscErrorCode (*basePolynomialsType1)(
      int p,double alpha,
      double x,double t,
      double *diff_x,double *diff_t,
      double *L,double *diffL,
      const int dim
    );

    JacobiPolynomialCtx(
      int p,
      double *diff_x,
      double *diff_t,
      int dim,
      double alpha,
      boost::shared_ptr<ublas::matrix<double> > base_fun_ptr,
      boost::shared_ptr<ublas::matrix<double> > base_diff_fun_ptr
    ):
    P(p),
    diffX(diff_x),
    diffT(diff_t),
    dIm(dim),
    aLpha(alpha),
    baseFunPtr(base_fun_ptr),
    baseDiffFunPtr(base_diff_fun_ptr),
    basePolynomialsType1(Jacobi_polynomials) {
    }
    ~JacobiPolynomialCtx() {}

  };

  /**
   * \brief Calculating Legendre base functions
   * \ingroup mofem_base_functions
   */
  struct JacobiPolynomial: public BaseFunction {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    JacobiPolynomial() {}
    ~JacobiPolynomial() {}

    PetscErrorCode getValue(
      ublas::matrix<double> &pts_x,
      ublas::matrix<double> &pts_t,
      boost::shared_ptr<BaseFunctionCtx> ctx_ptr
    );

  };

}

#endif //__JACOBIPOLYNOMIALS_HPP__
