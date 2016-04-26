/** \file TetPolynomialBase.hpp
\brief Implementation of Ainsworth-Cole H1 base on tetrahedral
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

#ifndef __EDGEPOLYNOMIALBASE_HPP__
#define __EDGEPOLYNOMIALBASE_HPP__

namespace MoFEM {

  /**
   * \brief Calculate base functions on tetrahedral
   * 
   * \ingroup mofem_base_functions
   */
  struct EdgePolynomialBase: public BaseFunction {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    EdgePolynomialBase();
    ~EdgePolynomialBase();

    PetscErrorCode getValue(
      ublas::matrix<double> &pts,
      boost::shared_ptr<BaseFunctionCtx> ctx_ptr
    );

  private:

    EntPolynomialBaseCtx *cTx;

    ublas::vector<double> L,diffL;

    PetscErrorCode getValueH1(ublas::matrix<double> &pts);

    PetscErrorCode getValueL2(ublas::matrix<double> &pts);

    PetscErrorCode getValueHdiv(ublas::matrix<double> &pts);

    PetscErrorCode getValueHCurl(ublas::matrix<double> &pts);


  };


}

#endif //__EDGEPOLYNOMIALBASE_HPP__
