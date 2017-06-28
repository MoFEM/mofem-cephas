/** \file FlatPrismPolynomialBase.hpp
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

#ifndef __FLATPRISMPOLYNOMIALBASE_HPP__
#define __FLATPRISMPOLYNOMIALBASE_HPP__

namespace MoFEM {

  static const MOFEMuuid IDD_FLATPRISM_BASE_FUNCTION = MOFEMuuid(BitIntefaceId(FLATPRISM_BASE_FUNCTION_INTERFACE));

  /**
  * \brief Class used to pass element data to calculate base functions on flat prism
  *
  * \ingroup mofem_base_functions
  * FIXME: Need moab and mofem finite element structure to work (that not perfect)
  */
  struct FlatPrismPolynomialBaseCtx: public EntPolynomialBaseCtx {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    moab::Interface &mOab;
    const NumeredEntFiniteElement *fePtr;

    FlatPrismPolynomialBaseCtx(
      DataForcesAndSurcesCore &data,
      moab::Interface &moab,
      const NumeredEntFiniteElement *fe_ptr,
      const FieldSpace space,
      const FieldApproximationBase base,
      const FieldApproximationBase copy_node_base = LASTBASE
    );

    ~FlatPrismPolynomialBaseCtx();

  };


  /**
  * \brief Calculate base functions on tetrahedral
  * \ingroup mofem_base_functions
  * FIXME: Need moab and mofem finite element structure to work (that not perfect)
  */
  struct FlatPrismPolynomialBase: public BaseFunction {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    FlatPrismPolynomialBase();
    ~FlatPrismPolynomialBase();

    PetscErrorCode getValue(
      MatrixDouble &pts,boost::shared_ptr<BaseFunctionCtx> ctx_ptr
    );

  private:

    FlatPrismPolynomialBaseCtx *cTx;

    PetscErrorCode getValueH1(MatrixDouble &pts);

    PetscErrorCode getValueL2(MatrixDouble &pts);

    PetscErrorCode getValueHdiv(MatrixDouble &pts);

    PetscErrorCode getValueHCurl(MatrixDouble &pts);

    int numNodes;
    const EntityHandle *connPrism;
    const EntityHandle *connFace3;
    const EntityHandle *connFace4;
    int faceNodes[2][3];
    MatrixDouble N;
    MatrixDouble diffN;

  };

}

#endif //__FLATPRISMPOLYNOMIALBASE_HPP__
