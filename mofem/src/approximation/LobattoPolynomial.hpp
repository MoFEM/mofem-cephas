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

  static const int LOBATTO_BASE_FUNCTION_INTERFACE = 1<<2;
  static const MOFEMuuid IDD_LOBATTO_BASE_FUNCTION = MOFEMuuid(BitIntefaceId(LEGENDRE_BASE_FUNCTION_INTERFACE));

  struct LobattoPolynomialCtx: public LegendrePolynomialCtx {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    LobattoPolynomialCtx(int p,double *diff_s,int dim):
    LegendrePolynomialCtx(p,diff_s,dim) {
      base_polynomials = Lobatto_polynomials;
    }
    ~LobattoPolynomialCtx() {}

  };

  struct LobattoPolynomial: public LegendrePolynomial {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    LobattoPolynomial() {}
    ~LobattoPolynomial() {}

    PetscErrorCode getValue(
      ublas::matrix<double> &pTs,
      boost::shared_ptr<ublas::matrix<double> > baseFunPtr,
      boost::shared_ptr<ublas::matrix<double> > baseDiffFunPtr,
      boost::shared_ptr<BaseFunctionCtx> ctxPtr
    );

  };

}

#endif
