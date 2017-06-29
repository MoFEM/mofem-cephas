/** \file FatPrismPolynomialBase.hpp
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

#ifndef __FATPRISMPOLYNOMIALBASE_HPP__
#define __FATPRISMPOLYNOMIALBASE_HPP__

namespace MoFEM {

  static const MOFEMuuid IDD_FATPRISM_BASE_FUNCTION = MOFEMuuid(BitIntefaceId(FATPRISM_BASE_FUNCTION_INTERFACE));

  /**
  * \brief Class used to pass element data to calculate base functions on fat prism
  *
  * \ingroup mofem_base_functions
  * FIXME: Need moab and mofem finite element structure to work (that not perfect)
  */
  struct FatPrismPolynomialBaseCtx: public EntPolynomialBaseCtx {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    DataForcesAndSurcesCore& dataTrianglesOnly;
    DataForcesAndSurcesCore& dataTroughThickness;

    MatrixDouble& gaussPtsTrianglesOnly;
    MatrixDouble& gaussPtsThroughThickness;

    moab::Interface &mOab;
    const NumeredEntFiniteElement *fePtr;

    FatPrismPolynomialBaseCtx(
      DataForcesAndSurcesCore &data,
      DataForcesAndSurcesCore &data_triangles_only,
      DataForcesAndSurcesCore &data_trough_thickness,
      MatrixDouble& gauss_pts_triangles_only,
      MatrixDouble& gauss_pts_through_thickness,
      moab::Interface &moab,
      const NumeredEntFiniteElement *fe_ptr,
      const FieldSpace space,
      const FieldApproximationBase base,
      const FieldApproximationBase copy_node_base = LASTBASE
    );

    ~FatPrismPolynomialBaseCtx();

  };


  /**
  * \brief Calculate base functions on tetrahedral
  * \ingroup mofem_base_functions
  * FIXME: Need moab and mofem finite element structure to work (that not perfect)
  */
  struct FatPrismPolynomialBase: public BaseFunction {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    FatPrismPolynomialBase();
    ~FatPrismPolynomialBase();

    PetscErrorCode getValue(
      MatrixDouble &pts,boost::shared_ptr<BaseFunctionCtx> ctx_ptr
    );

  private:

    FatPrismPolynomialBaseCtx *cTx;

    PetscErrorCode getValueH1TrianglesOnly();

    PetscErrorCode getValueH1ThroughThickness();

    PetscErrorCode getValueH1(MatrixDouble &pts);

    PetscErrorCode getValueL2(MatrixDouble &pts);

    PetscErrorCode getValueHdiv(MatrixDouble &pts);

    PetscErrorCode getValueHCurl(MatrixDouble &pts);

    int faceNodes[2][3];
    MatrixDouble N;
    MatrixDouble diffN;

  };

}

#endif //__FATPRISMPOLYNOMIALBASE_HPP__
