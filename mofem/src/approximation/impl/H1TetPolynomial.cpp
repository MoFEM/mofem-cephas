/** \file H1TetPolynomial.cpp
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

#include <version.h>
#include <config.h>
#include <definitions.h>
#include <Includes.hpp>

#include <base_functions.h>
#include <fem_tools.h>
#include <h1_hdiv_hcurl_l2.h>
#include <Common.hpp>
#include <UnknownInterface.hpp>
using namespace MoFEM;

#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <DataStructures.hpp>

#include <BaseFunction.hpp>
#include <H1TetPolynomial.hpp>

PetscErrorCode H1TetPolynomialCtx::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_H1TET_BASE_FUNCTION) {
    *iface = dynamic_cast<H1TetPolynomialCtx*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = BaseFunctionCtx::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode H1TetPolynomial::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_H1TET_BASE_FUNCTION) {
    *iface = dynamic_cast<H1TetPolynomial*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = BaseFunction::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

H1TetPolynomialCtx::H1TetPolynomialCtx(
  DataForcesAndSurcesCore &data,
  const FieldSpace space,
  const FieldApproximationBase base
):
dAta(data),
sPace(space),
bAse(base) {
  switch(bAse) {
    case AINSWORTH_COLE_BASE:
    basePolynomials = Legendre_polynomials;
    break;
    case LOBATTO_BASE:
    basePolynomials = Lobatto_polynomials;
    break;
    default:
    THROW_MESSAGE("Not implemented for this base")
  }
}

H1TetPolynomialCtx::~H1TetPolynomialCtx() {
}

H1TetPolynomial::~H1TetPolynomial() {}
H1TetPolynomial::H1TetPolynomial() {}

PetscErrorCode H1TetPolynomial::getValueH1(ublas::matrix<double> &pts) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  DataForcesAndSurcesCore data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(
    int p,double s,double *diff_s,double *L,double *diffL,const int dim
  ) = cTx->basePolynomials;

  int nb_gauss_pts = pts.size2();
  if(!nb_gauss_pts) {
    PetscFunctionReturn(0);
  }

  if(pts.size1()<3) {
    SETERRQ(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "Wrong dimension of pts, should be at least 3 rows with coordinates"
    );
  }

  if(data.dataOnEntities[MBVERTEX][0].getN(base).size2()!=4) {
    data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts,4,false);
    ierr = ShapeMBTET(
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &pts(0,0),
      &pts(1,0),
      &pts(2,0),
      nb_gauss_pts
    ); CHKERRQ(ierr);
  }
  if(data.dataOnEntities[MBVERTEX][0].getN(base).size1()!=nb_gauss_pts) {
    SETERRQ1(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "Base functions or nodes has wrong number of integration points for base %s",
      ApproximationBaseNames[base]
    );
  }

  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(4,3,false);
  ierr = ShapeDiffMBTET(&*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin()); CHKERRQ(ierr);

  int _sense_[6],_order_[6];
  if(data.spacesOnEntities[MBEDGE].test(H1)) {
    //edges
    if(data.dataOnEntities[MBEDGE].size()!=6) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    double *_H1edgeN_[6],*_diffH1edgeN_[6];
    for(int ee = 0;ee<6;ee++) {
      if(data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      _sense_[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      _order_[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();
      int nb_dofs = NBEDGE_H1_AINSWORTH_COLE(data.dataOnEntities[MBEDGE][ee].getDataOrder());
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts,nb_dofs,false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts,3*nb_dofs,false);
      _H1edgeN_[ee] = &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      _diffH1edgeN_[ee] = &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }
    ierr = H1_EdgeShapeFunctions_MBTET(
      _sense_,_order_,
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      _H1edgeN_,_diffH1edgeN_,nb_gauss_pts,base_polynomials
    ); CHKERRQ(ierr);
  }

  if(data.spacesOnEntities[MBTRI].test(H1)) {

    //faces
    if(data.dataOnEntities[MBTRI].size()!=4) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    double *_H1faceN_[4],*_diffH1faceN_[4];
    for(int ff = 0;ff<4;ff++) {
      if(data.dataOnEntities[MBTRI][ff].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      int nb_dofs = NBFACETRI_H1_AINSWORTH_COLE(data.dataOnEntities[MBTRI][ff].getDataOrder());
      _order_[ff] = data.dataOnEntities[MBTRI][ff].getDataOrder();
      data.dataOnEntities[MBTRI][ff].getN(base).resize(nb_gauss_pts,nb_dofs,false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(nb_gauss_pts,3*nb_dofs,false);
      _H1faceN_[ff] = &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin();
      _diffH1faceN_[ff] = &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin();
    }
    if(data.facesNodes.size1() != 4) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    if(data.facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    ierr = H1_FaceShapeFunctions_MBTET(
      &*data.facesNodes.data().begin(),
      _order_,
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      _H1faceN_,
      _diffH1faceN_,
      nb_gauss_pts,
      base_polynomials
    ); CHKERRQ(ierr);

  }

  if(data.spacesOnEntities[MBTET].test(H1)) {

    //volume
    data.dataOnEntities[MBTET][0].getBase() = base;
    int order = data.dataOnEntities[MBTET][0].getDataOrder();
    int nb_vol_dofs = NBVOLUMETET_H1_AINSWORTH_COLE(order);
    data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts,nb_vol_dofs,false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts,3*nb_vol_dofs,false);
    ierr = H1_VolumeShapeFunctions_MBTET(
      data.dataOnEntities[MBTET][0].getDataOrder(),
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      &*data.dataOnEntities[MBTET][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin(),
      nb_gauss_pts,
      base_polynomials
    ); CHKERRQ(ierr);

  }

  PetscFunctionReturn(0);
}

PetscErrorCode H1TetPolynomial::getValue(
  ublas::matrix<double> &pts,
  boost::shared_ptr<BaseFunctionCtx> ctx_ptr
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  MoFEM::UnknownInterface *iface;
  ierr = ctx_ptr->queryInterface(IDD_H1TET_BASE_FUNCTION,&iface); CHKERRQ(ierr);
  cTx = reinterpret_cast<H1TetPolynomialCtx*>(iface);

  switch (cTx->sPace) {
    case H1:
    ierr = getValueH1(pts); CHKERRQ(ierr);
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Not yet implemented");
  }

  PetscFunctionReturn(0);
}
