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

  static const MOFEMuuid IDD_LOBATTO_BASE_FUNCTION = MOFEMuuid(BitIntefaceId(LOBATTO_BASE_FUNCTION_INTERFACE));
  static const MOFEMuuid IDD_KERNEL_BASE_FUNCTION = MOFEMuuid(BitIntefaceId(KERNEL_BASE_FUNCTION_INTERFACE));

  /**
   * \brief Class used to give arguments to Lobatto base functions
   * \ingroup mofem_base_functions
   */
  struct LobattoPolynomialCtx: public LegendrePolynomialCtx {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    LobattoPolynomialCtx(
      int p,
      double *diff_s,
      int dim,
      boost::shared_ptr<ublas::matrix<double> > base_fun_ptr,
      boost::shared_ptr<ublas::matrix<double> > base_diff_fun_ptr
    ):
    LegendrePolynomialCtx(p,diff_s,dim,base_fun_ptr,base_diff_fun_ptr) {
      basePolynomials = Lobatto_polynomials;
    }
    ~LobattoPolynomialCtx() {}

  };

  /**
   * \brief Calculating Lobatto base functions
   * \ingroup mofem_base_functions
   */
  struct LobattoPolynomial: public LegendrePolynomial {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    LobattoPolynomial() {}
    ~LobattoPolynomial() {}

    PetscErrorCode getValue(
      ublas::matrix<double> &pts,
      boost::shared_ptr<BaseFunctionCtx> ctx_ptr
    );

  };

  /**
   * \brief Class used to give arguments to Kernel Lobatto base functions
   * \ingroup mofem_base_functions
   */
  struct KernelLobattoPolynomialCtx: public LegendrePolynomialCtx {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    KernelLobattoPolynomialCtx(
      int p,
      double *diff_s,
      int dim,
      boost::shared_ptr<ublas::matrix<double> > base_fun_ptr,
      boost::shared_ptr<ublas::matrix<double> > base_diff_fun_ptr
    ):
    LegendrePolynomialCtx(p,diff_s,dim,base_fun_ptr,base_diff_fun_ptr) {
      basePolynomials = LobattoKernel_polynomials;
    }
    ~KernelLobattoPolynomialCtx() {}

  };

  /**
   * \brief Calculating Lobatto base functions
   * \ingroup mofem_base_functions
   */
  struct KernelLobattoPolynomial: public LegendrePolynomial {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    KernelLobattoPolynomial() {}
    ~KernelLobattoPolynomial() {}

    PetscErrorCode getValue(
      ublas::matrix<double> &pts,
      boost::shared_ptr<BaseFunctionCtx> ctx_ptr
    );

  };
  
}

#endif
