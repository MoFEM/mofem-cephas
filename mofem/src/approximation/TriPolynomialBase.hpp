/** \file TriPolynomialBase.hpp
\brief Implementation of Ainsworth-Cole, Demkowicz  and other bases H1 base on triangle

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

#ifndef __H1TRIPOLYNOMIAL_HPP__
#define __H1TRIPOLYNOMIAL_HPP__

namespace MoFEM {

  /**
   * \brief Calculate base functions on triangle
   *
   * \ingroup mofem_base_functions
   */
  struct TriPolynomialBase: public BaseFunction {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    TriPolynomialBase();
    ~TriPolynomialBase();

    PetscErrorCode getValue(
      MatrixDouble &pts,
      boost::shared_ptr<BaseFunctionCtx> ctx_ptr
    );

  private:

    EntPolynomialBaseCtx *cTx;

    PetscErrorCode getValueH1(
      MatrixDouble &pts
    );

    PetscErrorCode getValueL2(
      MatrixDouble &pts
    );

    ublas::matrix<MatrixDouble > N_face_edge;
    ublas::vector<MatrixDouble > N_face_bubble;
    ublas::matrix<MatrixDouble > diffN_face_edge;
    ublas::vector<MatrixDouble > diffN_face_bubble;

    PetscErrorCode getValueHdiv(
      MatrixDouble &pts
    );

    PetscErrorCode getValueHCurl(
      MatrixDouble &pts
    );

  private:

    PetscErrorCode getValueHdivAinsworthBase(
      MatrixDouble &pts
    );

    PetscErrorCode getValueHdivDemkowiczBase(
      MatrixDouble &pts
    );

  };


}

#endif //__H1TRIPOLYNOMIAL_HPP__
