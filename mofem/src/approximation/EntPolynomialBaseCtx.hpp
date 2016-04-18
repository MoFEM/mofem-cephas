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

#ifndef __H1ENTPOLYNOMIALCTX_HPP__
#define __H1ENTPOLYNOMIALCTX_HPP__

namespace MoFEM {

  static const MOFEMuuid IDD_TET_BASE_FUNCTION = MOFEMuuid(BitIntefaceId(TET_BASE_FUNCTION_INTERFACE));
  static const MOFEMuuid IDD_TRI_BASE_FUNCTION = MOFEMuuid(BitIntefaceId(TRI_BASE_FUNCTION_INTERFACE));

  struct EntPolynomialBaseCtx: public BaseFunctionCtx {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    PetscErrorCode (*basePolynomials)(
      int p,double s,double *diff_s,double *L,double *diffL,const int dim
    );
    DataForcesAndSurcesCore &dAta;
    const FieldSpace sPace;
    const FieldApproximationBase bAse;
    const FieldApproximationBase copyNodeBase;

    EntPolynomialBaseCtx(
      DataForcesAndSurcesCore &data,
      const FieldSpace space,
      const FieldApproximationBase base,
      const FieldApproximationBase copy_node_base = LASTBASE
    );
    ~EntPolynomialBaseCtx();

  };

}

#endif //__H1ENTPOLYNOMIALCTX_HPP__
