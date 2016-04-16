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

  static const int LEGENDRE_BASE_FUNCTION_INTERFACE = 1<<1;
  static const MOFEMuuid IDD_LEGENDRE_BASE_FUNCTION = MOFEMuuid(BitIntefaceId(LEGENDRE_BASE_FUNCTION_INTERFACE));

  struct LegendrePolynomialCtx: public BaseFunctionCtx {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    int P;
    double *diffS;
    int dIm;

    PetscErrorCode (*base_polynomials)(
      int p,double s,double *diff_s,double *L,double *diffL,const int dim
    );

    LegendrePolynomialCtx(int p,double *diff_s,int dim):
    P(p),
    diffS(diff_s),
    dIm(dim),
    base_polynomials(Legendre_polynomials) {
    }
    ~LegendrePolynomialCtx() {}

  };

  struct LegendrePolynomial: public BaseFunction {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    LegendrePolynomial() {}
    ~LegendrePolynomial() {}

    PetscErrorCode getValue(
      ublas::matrix<double> &pTs,
      boost::shared_ptr<ublas::matrix<double> > baseFunPtr,
      boost::shared_ptr<ublas::matrix<double> > baseDiffFunPtr,
      boost::shared_ptr<BaseFunctionCtx> ctxPtr
    );

  };

}


#endif //__LEGENDREPOLYNOMIALS_HPP__
