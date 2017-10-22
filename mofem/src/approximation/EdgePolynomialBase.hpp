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

    PetscErrorCode query_interface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface) const;

    EdgePolynomialBase();
    ~EdgePolynomialBase();

    PetscErrorCode getValue(
      MatrixDouble &pts,
      boost::shared_ptr<BaseFunctionCtx> ctx_ptr
    );

  private:

    EntPolynomialBaseCtx *cTx;

    VectorDouble L,diffL;

    PetscErrorCode getValueH1(MatrixDouble &pts);

    PetscErrorCode getValueL2(MatrixDouble &pts);

    PetscErrorCode getValueHdiv(MatrixDouble &pts);

    PetscErrorCode getValueHCurl(MatrixDouble &pts);


  };


}

#endif //__EDGEPOLYNOMIALBASE_HPP__
