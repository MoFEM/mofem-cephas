/** \file FlatPrismPolynomialBase.cpp
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
#include <FEMMultiIndices.hpp>
#include <DataStructures.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <LoopMethods.hpp>

#include <BaseFunction.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <FlatPrismPolynomialBase.hpp>

PetscErrorCode FlatPrismPolynomialBaseCtx::queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(
    uuid == IDD_FLATPRISM_BASE_FUNCTION
  ) {
    *iface = static_cast<FlatPrismPolynomialBaseCtx*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = EntPolynomialBaseCtx::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

FlatPrismPolynomialBaseCtx::FlatPrismPolynomialBaseCtx(
  DataForcesAndSurcesCore &data,
  moab::Interface &moab,
  const NumeredEntFiniteElement *fe_ptr,
  const FieldSpace space,
  const FieldApproximationBase base,
  const FieldApproximationBase copy_node_base
):
EntPolynomialBaseCtx(data,space,base,copy_node_base),
mOab(moab),
fePtr(fe_ptr) {
  PetscErrorCode ierr;
  ierr = setBase(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}
FlatPrismPolynomialBaseCtx::~FlatPrismPolynomialBaseCtx() {
}

PetscErrorCode FlatPrismPolynomialBase::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_FLATPRISM_BASE_FUNCTION) {
    *iface = static_cast<FlatPrismPolynomialBase*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = BaseFunction::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

FlatPrismPolynomialBase::~FlatPrismPolynomialBase() {}
FlatPrismPolynomialBase::FlatPrismPolynomialBase() {}

PetscErrorCode FlatPrismPolynomialBase::getValue(
  ublas::matrix<double> &pts,
  boost::shared_ptr<BaseFunctionCtx> ctx_ptr
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  MoFEM::UnknownInterface *iface;
  ierr = ctx_ptr->queryInterface(IDD_FLATPRISM_BASE_FUNCTION,&iface); CHKERRQ(ierr);
  cTx = reinterpret_cast<FlatPrismPolynomialBaseCtx*>(iface);
  if(!cTx->fePtr) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
      "Pointer to element should be given "
      "when EntPolynomialBaseCtx is constructed "
      "(use different constructor)"
    );
  }

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
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Not implemented");
  } else {
    data.dataOnEntities[MBVERTEX][0].getNSharedPtr(base) = data.dataOnEntities[MBVERTEX][0].getNSharedPtr(cTx->copyNodeBase);
  }
  data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts,6,false);
  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(nb_gauss_pts,12,false);
  if(data.dataOnEntities[MBVERTEX][0].getN(base).size1()!=(unsigned int)nb_gauss_pts) {
    SETERRQ1(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "Base functions or nodes has wrong number of integration points for base %s",
      ApproximationBaseNames[base]
    );
  }

  data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts,6,false);
  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(nb_gauss_pts,12,false);
  N.resize(nb_gauss_pts,3,false);
  diffN.resize(3,2,false);
  ierr = ShapeMBTRI(&*N.data().begin(),&pts(0,0),&pts(1,0),nb_gauss_pts); CHKERRQ(ierr);
  ierr = ShapeDiffMBTRI(&*diffN.data().begin()); CHKERRQ(ierr);

  // This is needed to have proper order of nodes on faces
  MoABErrorCode rval;
  rval = cTx->mOab.get_connectivity(cTx->fePtr->getEnt(),connPrism,numNodes,true); CHKERRQ_MOAB(rval);
  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(cTx->fePtr->get_side_number_table());
  SideNumber_multiIndex::nth_index<1>::type::iterator siit3 = side_table.get<1>().find(boost::make_tuple(MBTRI,3));
  SideNumber_multiIndex::nth_index<1>::type::iterator siit4 = side_table.get<1>().find(boost::make_tuple(MBTRI,4));
  if(siit3==side_table.get<1>().end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  if(siit4==side_table.get<1>().end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  rval = cTx->mOab.get_connectivity(siit3->get()->ent,connFace3,numNodes,true); CHKERRQ_MOAB(rval);
  rval = cTx->mOab.get_connectivity(siit4->get()->ent,connFace4,numNodes,true); CHKERRQ_MOAB(rval);

  for(int nn = 0;nn<3;nn++) {
    faceNodes[0][nn] = std::distance(connPrism,std::find(connPrism,connPrism+3,connFace3[nn]));
    faceNodes[1][nn] = std::distance(connPrism+3,std::find(connPrism+3,connPrism+6,connFace4[nn]));
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      double val = N(gg,nn);
      double val_x = diffN(nn,0);
      double val_y = diffN(nn,1);
      data.dataOnEntities[MBVERTEX][0].getN(base)(gg,nn) = val;
      data.dataOnEntities[MBVERTEX][0].getDiffN(base)(gg,2*nn+0) = val_x;
      data.dataOnEntities[MBVERTEX][0].getDiffN(base)(gg,2*nn+1) = val_y;
      data.dataOnEntities[MBVERTEX][0].getN(base)(gg,3+nn) = val;
      data.dataOnEntities[MBVERTEX][0].getDiffN(base)(gg,6+2*nn+0) = val_x;
      data.dataOnEntities[MBVERTEX][0].getDiffN(base)(gg,6+2*nn+1) = val_y;
    }
  }
  // for(int nn = 0;nn<3;nn++) {
  //   if(faceNodes[0][nn]!=faceNodes[1][nn]) {
  //     SETERRQ(
  //       PETSC_COMM_SELF,
  //       MOFEM_DATA_INCONSISTENCY,
  //       "Node order different on both faces"
  //     );
  //   }
  // }

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

PetscErrorCode FlatPrismPolynomialBase::getValueH1(ublas::matrix<double> &pts) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  DataForcesAndSurcesCore& data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(
    int p,double s,double *diff_s,double *L,double *diffL,const int dim
  ) = cTx->basePolynomials;

  int nb_gauss_pts = pts.size2();

  //edges
  if(data.dataOnEntities[MBEDGE].size()!=9) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  int valid_edges[] = { 1,1,1, 0,0,0, 1,1,1 };
  int sense[9],order[9];
  double *H1edgeN[9],*diffH1edgeN[9];
  if((data.spacesOnEntities[MBEDGE]).test(H1)) {
    for(int ee = 0;ee<9;ee++) {
      if(!valid_edges[ee]) continue;
      if(data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      order[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();
      int nb_dofs = NBEDGE_H1_AINSWORTH_COLE(data.dataOnEntities[MBEDGE][ee].getDataOrder());
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts,nb_dofs,false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts,2*nb_dofs,false);
      H1edgeN[ee] = &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      diffH1edgeN[ee] = &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }
    //shape functions on face 3
    ierr = H1_EdgeShapeFunctions_MBTRI(
      &sense[0],
      &order[0],
      &*N.data().begin(),
      &*diffN.data().begin(),
      &H1edgeN[0],
      &diffH1edgeN[0],
      nb_gauss_pts,base_polynomials
    ); CHKERRQ(ierr);
    //shape functions on face 4
    ierr = H1_EdgeShapeFunctions_MBTRI(
      &sense[6],
      &order[6],
      &*N.data().begin(),
      &*diffN.data().begin(),
      &H1edgeN[6],
      &diffH1edgeN[6],
      nb_gauss_pts,base_polynomials
    ); CHKERRQ(ierr);
  }

  //face
  if(data.dataOnEntities[MBTRI].size()!=5) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  if((data.spacesOnEntities[MBTRI]).test(H1)) {
    for(int ff = 3;ff<=4;ff++) {
      int nb_dofs = NBFACETRI_H1_AINSWORTH_COLE(data.dataOnEntities[MBTRI][ff].getDataOrder());
      data.dataOnEntities[MBTRI][ff].getN(base).resize(nb_gauss_pts,nb_dofs,false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(nb_gauss_pts,2*nb_dofs,false);
      ierr = H1_FaceShapeFunctions_MBTRI(
        faceNodes[ff-3],
        data.dataOnEntities[MBTRI][ff].getDataOrder(),
        &*N.data().begin(),
        &*diffN.data().begin(),
        &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin(),
        &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin(),
        nb_gauss_pts,base_polynomials
      ); CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode FlatPrismPolynomialBase::getValueL2(ublas::matrix<double> &pts) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Not yet implemented");
  PetscFunctionReturn(0);
}

PetscErrorCode FlatPrismPolynomialBase::getValueHdiv(ublas::matrix<double> &pts) {
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Not yet implemented");
}

PetscErrorCode FlatPrismPolynomialBase::getValueHCurl(ublas::matrix<double> &pts) {
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Not yet implemented");
}
