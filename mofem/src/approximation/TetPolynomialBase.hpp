/** \file TetPolynomialBase.hpp
\brief Implementation of Ainsworth-Coyle / Demkowicz or any other H1 base on tetrahedral

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

#ifndef __TETPOLYNOMIALBASE_HPP__
#define __TETPOLYNOMIALBASE_HPP__

namespace MoFEM {

  /**
   * \brief Calculate base functions on tetrahedral
   *
   * \ingroup mofem_base_functions
   */
  struct TetPolynomialBase: public BaseFunction {

    PetscErrorCode query_interface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface) const;

    TetPolynomialBase();
    ~TetPolynomialBase();

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
    ublas::vector<MatrixDouble > N_volume_edge;
    ublas::vector<MatrixDouble > N_volume_face;
    MatrixDouble N_volume_bubble;

    ublas::matrix<MatrixDouble > diffN_face_edge;
    ublas::vector<MatrixDouble > diffN_face_bubble;
    ublas::vector<MatrixDouble > diffN_volume_edge;
    ublas::vector<MatrixDouble > diffN_volume_face;
    MatrixDouble diffN_volume_bubble;


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

#endif //__TETPOLYNOMIALBASE_HPP__
