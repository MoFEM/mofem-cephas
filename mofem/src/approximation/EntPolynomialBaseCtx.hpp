/** \file EntPolynomialBaseCtx.hpp
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

#ifndef __ENTPOLYNOMIALCTX_HPP__
#define __ENTPOLYNOMIALCTX_HPP__

namespace MoFEM {

  static const MOFEMuuid IDD_TET_BASE_FUNCTION = MOFEMuuid(BitIntefaceId(TET_BASE_FUNCTION_INTERFACE));
  static const MOFEMuuid IDD_TRI_BASE_FUNCTION = MOFEMuuid(BitIntefaceId(TRI_BASE_FUNCTION_INTERFACE));
  static const MOFEMuuid IDD_EDGE_BASE_FUNCTION = MOFEMuuid(BitIntefaceId(EDGE_BASE_FUNCTION_INTERFACE));

  /**
   * \brief Class used to pass element data to calculate base functions on tet,triangle,edge
   *
   * \ingroup mofem_base_functions
   */
  struct EntPolynomialBaseCtx: public BaseFunctionCtx {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    PetscErrorCode (*basePolynomialsType0)(
      int p,double s,double *diff_s,double *L,double *diffL,const int dim
    );

    PetscErrorCode (*basePolynomialsType1)(
      int p,double alpha,
      double x,double t,
      double *diff_x,double *diff_t,
      double *L,double *diffL,
      const int dim
    );

    DataForcesAndSourcesCore &dAta;
    const FieldSpace sPace;
    const FieldApproximationBase bAse;
    const FieldApproximationBase copyNodeBase;
    const FEMethod *fePtr;

    EntPolynomialBaseCtx(
      DataForcesAndSourcesCore &data,
      const FieldSpace space,
      const FieldApproximationBase base,
      const FieldApproximationBase copy_node_base = LASTBASE
    );

    ~EntPolynomialBaseCtx();

  protected:

    PetscErrorCode setBase();

  };

}

#endif //__ENTPOLYNOMIALCTX_HPP__
