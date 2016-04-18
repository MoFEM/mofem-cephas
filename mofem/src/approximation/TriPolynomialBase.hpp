/** \file TriPolynomialBase.hpp
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

#ifndef __H1TRIPOLYNOMIAL_HPP__
#define __H1TRIPOLYNOMIAL_HPP__

namespace MoFEM {

  static const MOFEMuuid IDD_H1TRI_BASE_FUNCTION = MOFEMuuid(BitIntefaceId(H1TRI_BASE_FUNCTION_INTERFACE));

  struct TriPolynomialBaseCtx: public BaseFunctionCtx {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    PetscErrorCode (*basePolynomials)(
      int p,double s,double *diff_s,double *L,double *diffL,const int dim
    );
    DataForcesAndSurcesCore &dAta;
    const FieldSpace sPace;
    const FieldApproximationBase bAse;
    const FieldApproximationBase copyNodeBase;

    TriPolynomialBaseCtx(
      DataForcesAndSurcesCore &data,
      const FieldSpace space,
      const FieldApproximationBase base,
      const FieldApproximationBase copy_node_base = LASTBASE
    );
    ~TriPolynomialBaseCtx();

  };

  struct TriPolynomialBase: public BaseFunction {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    TriPolynomialBase();
    ~TriPolynomialBase();

    PetscErrorCode getValue(
      ublas::matrix<double> &pts,
      boost::shared_ptr<BaseFunctionCtx> ctx_ptr
    );

  private:

    TriPolynomialBaseCtx *cTx;

    PetscErrorCode getValueH1(
      ublas::matrix<double> &pts
    );

    PetscErrorCode getValueL2(
      ublas::matrix<double> &pts
    );

    ublas::matrix<MatrixDouble > N_face_edge;
    ublas::vector<MatrixDouble > N_face_bubble;
    ublas::matrix<MatrixDouble > diffN_face_edge;
    ublas::vector<MatrixDouble > diffN_face_bubble;

    PetscErrorCode getValueHdiv(
      ublas::matrix<double> &pts
    );

    PetscErrorCode getValueHCurl(
      ublas::matrix<double> &pts
    );


  };


}

#endif //__H1TRIPOLYNOMIAL_HPP__
