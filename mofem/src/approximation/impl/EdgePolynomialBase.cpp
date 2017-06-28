/** \file EdgePolynomialBase.cpp
\brief Implementation of Ainsworth-Cole H1 base on edge
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

#include <FTensor.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <DataStructures.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <LoopMethods.hpp>

#include <BaseFunction.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <EdgePolynomialBase.hpp>

#include <Hcurl.hpp>

PetscErrorCode EdgePolynomialBase::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_EDGE_BASE_FUNCTION) {
    *iface = static_cast<EdgePolynomialBase*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = BaseFunction::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

EdgePolynomialBase::~EdgePolynomialBase() {}
EdgePolynomialBase::EdgePolynomialBase() {}

PetscErrorCode EdgePolynomialBase::getValue(
  MatrixDouble &pts,
  boost::shared_ptr<BaseFunctionCtx> ctx_ptr
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  MoFEM::UnknownInterface *iface;
  ierr = ctx_ptr->queryInterface(IDD_TET_BASE_FUNCTION,&iface); CHKERRQ(ierr);
  cTx = reinterpret_cast<EntPolynomialBaseCtx*>(iface);

  int nb_gauss_pts = pts.size2();
  if(!nb_gauss_pts) {
    PetscFunctionReturn(0);
  }

  if(pts.size1()<1) {
    SETERRQ(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "Wrong dimension of pts, should be at least 3 rows with coordinates"
    );
  }

  const FieldApproximationBase base = cTx->bAse;
  DataForcesAndSurcesCore& data = cTx->dAta;
  if(cTx->copyNodeBase==LASTBASE) {
    data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts,2,false);
    ierr = ShapeMBEDGE(
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &pts(0,0),
      nb_gauss_pts
    ); CHKERRQ(ierr);
  } else {
    data.dataOnEntities[MBVERTEX][0].getNSharedPtr(base) = data.dataOnEntities[MBVERTEX][0].getNSharedPtr(cTx->copyNodeBase);
  }
  if(data.dataOnEntities[MBVERTEX][0].getN(base).size1()!=(unsigned int)nb_gauss_pts) {
    SETERRQ1(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "Base functions or nodes has wrong number of integration points for base %s",
      ApproximationBaseNames[base]
    );
  }
  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(2,1,false);
  ierr = ShapeDiffMBEDGE(
    &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin()
  ); CHKERRQ(ierr);

  switch (cTx->sPace) {
    case H1:
    ierr = getValueH1(pts); CHKERRQ(ierr);
    break;
    case HDIV:
    ierr = getValueHdiv(pts); CHKERRQ(ierr);
    break;
    case HCURL:
    ierr = getValueHCurl(pts); CHKERRQ(ierr);
    break;
    case L2:
    ierr = getValueL2(pts); CHKERRQ(ierr);
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Not yet implemented");
  }

  PetscFunctionReturn(0);
}

PetscErrorCode EdgePolynomialBase::getValueH1(MatrixDouble &pts) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  DataForcesAndSurcesCore& data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(
    int p,double s,double *diff_s,double *L,double *diffL,const int dim
  ) = cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();


  //std::cerr << data.dataOnEntities[MBVERTEX][0].getN(base) << std::endl;
  //std::cerr << data.dataOnEntities[MBVERTEX][0].getDiffN(base) << std::endl;
  //
  // std::cerr << pts << std::endl;

  const int side_number = 0;
  int sense = data.dataOnEntities[MBEDGE][side_number].getSense();
  int order = data.dataOnEntities[MBEDGE][side_number].getDataOrder();
  data.dataOnEntities[MBEDGE][side_number].getN(base).resize(nb_gauss_pts,NBEDGE_H1(order),false);
  data.dataOnEntities[MBEDGE][side_number].getDiffN(base).resize(nb_gauss_pts,NBEDGE_H1(order),false);

  data.dataOnEntities[MBEDGE][side_number].getN(base).clear();
  data.dataOnEntities[MBEDGE][side_number].getDiffN(base).clear();

  L.resize(NBEDGE_H1(order),false);
  diffL.resize(NBEDGE_H1(order),false);

  // std::cerr << data.dataOnEntities[MBVERTEX][0].getN(base) << std::endl;

  if(data.dataOnEntities[MBEDGE][side_number].getDataOrder()>1) {

    double diff_s = 2.; // s = s(xi), ds/dxi = 2., because change of basis
    for(int gg = 0;gg<nb_gauss_pts;gg++) {

      double s = 2*pts(0,gg)-1; // makes form -1..1
      if(!sense) {
        s *= -1;
        diff_s *= -1;
      }

      // calculate Legendre polynomials at integration points
      ierr = base_polynomials(
        NBEDGE_H1(order)-1,s,&diff_s,&*L.data().begin(),&*diffL.data().begin(),1
      ); CHKERRQ(ierr);

      // std::cerr << "s " << s << " " << L << std::endl;

      for(unsigned int pp = 0;pp<data.dataOnEntities[MBEDGE][side_number].getN(base).size2();pp++) {

        // Calculate edge shape functions N0*N1*L(p), where N0 and N1 are nodal shape functions
        double v = data.dataOnEntities[MBVERTEX][0].getN(base)(gg,0)*data.dataOnEntities[MBVERTEX][0].getN(base)(gg,1);
        data.dataOnEntities[MBEDGE][side_number].getN(base)(gg,pp) = v*L(pp);

        // Calculate derivative edge shape functions
        // dN/dksi = dN0/dxi*N1*L + N0*dN1/ksi*L + N0*N1*dL/dxi
        data.dataOnEntities[MBEDGE][side_number].getDiffN(base)(gg,pp) =
        ((+1.)*data.dataOnEntities[MBVERTEX][0].getN(base)(gg,1)
        +
        data.dataOnEntities[MBVERTEX][0].getN(base)(gg,0)*(-1.))*L(pp)
        +
        v*diffL(pp);

      }
    }
  }

  // std::cerr << data.dataOnEntities[MBEDGE][0].getN(base) << std::endl;
  //std::cerr << data.dataOnEntities[MBEDGE][0].getDiffN(base) << std::endl;

  PetscFunctionReturn(0);
}

PetscErrorCode EdgePolynomialBase::getValueL2(MatrixDouble &pts) {
  PetscFunctionBegin;
  SETERRQ(
    PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
    "Make no sense, unless problem is 1d (1d not implemented yet)"
  );
  PetscFunctionReturn(0);
}

PetscErrorCode EdgePolynomialBase::getValueHdiv(MatrixDouble &pts) {
  SETERRQ(
    PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
    "Make no sense, unless problem is 2d (2d not implemented yet)"
  );
}

PetscErrorCode EdgePolynomialBase::getValueHCurl(MatrixDouble &pts) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  DataForcesAndSurcesCore& data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(
    int p,double s,double *diff_s,double *L,double *diffL,const int dim
  ) = cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  if(data.spacesOnEntities[MBEDGE].test(HCURL)) {
    //edges
    if(data.dataOnEntities[MBEDGE].size()!=1) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    int sense = data.dataOnEntities[MBEDGE][0].getSense();
    int order = data.dataOnEntities[MBEDGE][0].getDataOrder();
    int nb_dofs = NBEDGE_HCURL(data.dataOnEntities[MBEDGE][0].getDataOrder());
    data.dataOnEntities[MBEDGE][0].getN(base).resize(nb_gauss_pts,3*nb_dofs,false);
    data.dataOnEntities[MBEDGE][0].getDiffN(base).resize(nb_gauss_pts,0,false);
    // cerr << data.dataOnEntities[MBVERTEX][0].getDiffN(base) << endl;
    ierr = Hcurl_EdgeBaseFunctions_MBTET_ON_EDGE(
      sense,
      order,
      &data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
      &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
      &*data.dataOnEntities[MBEDGE][0].getN(base).data().begin(),
      NULL,
      nb_gauss_pts,
      base_polynomials
    ); CHKERRQ(ierr);
  } else {
    data.dataOnEntities[MBEDGE][0].getN(base).resize(nb_gauss_pts,0,false);
    data.dataOnEntities[MBEDGE][0].getDiffN(base).resize(nb_gauss_pts,0,false);
  }

  PetscFunctionReturn(0);
}
