/** \file TetPolynomialBase.cpp
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
#include <EntPolynomialBaseCtx.hpp>
#include <TetPolynomialBase.hpp>

PetscErrorCode TetPolynomialBase::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_TET_BASE_FUNCTION) {
    *iface = dynamic_cast<TetPolynomialBase*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = BaseFunction::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

TetPolynomialBase::~TetPolynomialBase() {}
TetPolynomialBase::TetPolynomialBase() {}

PetscErrorCode TetPolynomialBase::getValueH1(ublas::matrix<double> &pts) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  DataForcesAndSurcesCore& data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(
    int p,double s,double *diff_s,double *L,double *diffL,const int dim
  ) = cTx->basePolynomials;

  int nb_gauss_pts = pts.size2();

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

    // cerr << "Aaaaa2\n";
    // cerr << data.dataOnEntities[MBVERTEX][0].getN(base) << endl;
    // cerr << ApproximationBaseNames[base] << endl;

  }

  if(data.spacesOnEntities[MBTET].test(H1)) {

    //volume
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

PetscErrorCode TetPolynomialBase::getValueL2(
  ublas::matrix<double> &pts
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  DataForcesAndSurcesCore& data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(
    int p,double s,double *diff_s,double *L,double *diffL,const int dim
  ) = cTx->basePolynomials;

  int nb_gauss_pts = pts.size2();

  data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts,NBVOLUMETET_L2_AINSWORTH_COLE(data.dataOnEntities[MBTET][0].getDataOrder()),false);
  data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts,3*NBVOLUMETET_L2_AINSWORTH_COLE(data.dataOnEntities[MBTET][0].getDataOrder()),false);

  ierr = L2_ShapeFunctions_MBTET(
    data.dataOnEntities[MBTET][0].getDataOrder(),
    &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
    &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
    &*data.dataOnEntities[MBTET][0].getN(base).data().begin(),
    &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin(),
    nb_gauss_pts
  ); CHKERRQ(ierr);


  PetscFunctionReturn(0);
}

PetscErrorCode TetPolynomialBase::getValueHdiv(
  ublas::matrix<double> &pts
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  DataForcesAndSurcesCore& data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(
    int p,double s,double *diff_s,double *L,double *diffL,const int dim
  ) = cTx->basePolynomials;

  int nb_gauss_pts = pts.size2();

  //face shape functions

  double *phi_f_e[4][3];
  double *phi_f[4];
  double *diff_phi_f_e[4][3];
  double *diff_phi_f[4];

  N_face_edge.resize(4,3,false);
  N_face_bubble.resize(4,false);
  diffN_face_edge.resize(4,3,false);
  diffN_face_bubble.resize(4,false);

  int faces_order[4];
  for(int ff = 0;ff<4;ff++) {
    if(data.dataOnEntities[MBTRI][ff].getSense() == 0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    faces_order[ff] = data.dataOnEntities[MBTRI][ff].getDataOrder();
    //three edges on face
    for(int ee = 0;ee<3;ee++) {
      N_face_edge(ff,ee).resize(nb_gauss_pts,3*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(faces_order[ff]),false);
      diffN_face_edge(ff,ee).resize(nb_gauss_pts,9*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(faces_order[ff]),false);
      phi_f_e[ff][ee] = &((N_face_edge(ff,ee))(0,0));
      diff_phi_f_e[ff][ee] = &((diffN_face_edge(ff,ee))(0,0));
    }
    N_face_bubble[ff].resize(nb_gauss_pts,3*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(faces_order[ff]),false);
    diffN_face_bubble[ff].resize(nb_gauss_pts,9*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(faces_order[ff]),false);
    phi_f[ff] = &*(N_face_bubble[ff].data().begin());
    diff_phi_f[ff] = &*(diffN_face_bubble[ff].data().begin());
  }

  ierr = Hdiv_EdgeFaceShapeFunctions_MBTET(
    &data.facesNodes(0,0),faces_order,
    &data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
    phi_f_e,diff_phi_f_e,nb_gauss_pts,
    base_polynomials
  ); CHKERRQ(ierr);

  ierr = Hdiv_FaceBubbleShapeFunctions_MBTET(
    &data.facesNodes(0,0),faces_order,
    &data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
    phi_f,diff_phi_f,nb_gauss_pts,
    base_polynomials
  ); CHKERRQ(ierr);

  //volume shape functions

  double *phi_v_e[6];
  double *phi_v_f[4];
  double *phi_v;
  double *diff_phi_v_e[6];
  double *diff_phi_v_f[4];
  double *diff_phi_v;

  int volume_order = data.dataOnEntities[MBTET][0].getDataOrder();
  double coords[] = { 0,0,0, 1,0,0, 0,1,0, 0,0,1 };

  N_volume_edge.resize(6,false);
  diffN_volume_edge.resize(6,false);
  for(int ee = 0;ee<6;ee++) {
    N_volume_edge[ee].resize(nb_gauss_pts,3*NBVOLUMETET_EDGE_HDIV_AINSWORTH_COLE(volume_order),false);
    diffN_volume_edge[ee].resize(nb_gauss_pts,9*NBVOLUMETET_EDGE_HDIV_AINSWORTH_COLE(volume_order),false);
    phi_v_e[ee] = &*(N_volume_edge[ee].data().begin());
    diff_phi_v_e[ee] = &*(diffN_volume_edge[ee].data().begin());
  }
  ierr = Hdiv_EdgeBasedVolumeShapeFunctions_MBTET(
    volume_order,coords,&data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
    phi_v_e,diff_phi_v_e,nb_gauss_pts,
    base_polynomials
  ); CHKERRQ(ierr);

  N_volume_face.resize(4,false);
  diffN_volume_face.resize(4,false);
  for(int ff = 0;ff<4;ff++) {
    N_volume_face[ff].resize(nb_gauss_pts,3*NBVOLUMETET_FACE_HDIV_AINSWORTH_COLE(volume_order),false);
    diffN_volume_face[ff].resize(nb_gauss_pts,9*NBVOLUMETET_FACE_HDIV_AINSWORTH_COLE(volume_order),false);
    phi_v_f[ff] = &*(N_volume_face[ff].data().begin());
    diff_phi_v_f[ff] = &*(diffN_volume_face[ff].data().begin());
  }
  ierr = Hdiv_FaceBasedVolumeShapeFunctions_MBTET(
    volume_order,coords,&data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
    phi_v_f,diff_phi_v_f,nb_gauss_pts,
    base_polynomials
  ); CHKERRQ(ierr);

  N_volume_bubble.resize(nb_gauss_pts,3*NBVOLUMETET_VOLUME_HDIV_AINSWORTH_COLE(volume_order),false);
  diffN_volume_bubble.resize(nb_gauss_pts,9*NBVOLUMETET_VOLUME_HDIV_AINSWORTH_COLE(volume_order),false);
  phi_v = &*(N_volume_bubble.data().begin());
  diff_phi_v = &*(diffN_volume_bubble.data().begin());
  ierr = Hdiv_VolumeBubbleShapeFunctions_MBTET(
    volume_order,coords,&data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
    phi_v,diff_phi_v,nb_gauss_pts,
    base_polynomials
  ); CHKERRQ(ierr);

  // Set shape functions into data strucrure Shape functions hast to be put
  // in arrays in order which guarantee hierarhical series of digrees of
  // freedom, i.e. in other words dofs form sub-entities has to be group
  // by order.

  //faces
  if(data.dataOnEntities[MBTRI].size()!=4) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  for(int ff = 0;ff<4;ff++) {
    data.dataOnEntities[MBTRI][ff].getHdivN(base).resize(nb_gauss_pts,3*NBFACETRI_HDIV_AINSWORTH_COLE(faces_order[ff]),false);
    data.dataOnEntities[MBTRI][ff].getDiffHdivN(base).resize(nb_gauss_pts,9*NBFACETRI_HDIV_AINSWORTH_COLE(faces_order[ff]),false);
    int col = 0,diff_col = 0;
    for(int oo = 0;oo<faces_order[ff];oo++) {
      for(int ee = 0;ee<3;ee++) {
        //values
        for(int dd = 3*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(oo);dd<3*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(oo+1);dd++,col++) {
          for(int gg = 0;gg!=nb_gauss_pts;gg++) {
            data.dataOnEntities[MBTRI][ff].getHdivN(base)(gg,col) = N_face_edge(ff,ee)(gg,dd);
          }
        }
        //direvatives
        for(int dd = 9*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(oo);dd<9*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(oo+1);dd++,diff_col++) {
          for(int gg = 0;gg!=nb_gauss_pts;gg++) {
            data.dataOnEntities[MBTRI][ff].getDiffHdivN(base)(gg,diff_col) = diffN_face_edge(ff,ee)(gg,dd);
          }
        }
      }
      //values
      for(int dd = 3*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(oo);dd<3*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(oo+1);dd++,col++) {
        for(int gg = 0;gg!=nb_gauss_pts;gg++) {
          data.dataOnEntities[MBTRI][ff].getHdivN(base)(gg,col) = N_face_bubble[ff](gg,dd);
        }
      }
      //direvatives
      for(int dd = 9*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(oo);dd<9*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(oo+1);dd++,diff_col++) {
        for(int gg = 0;gg!=nb_gauss_pts;gg++) {
          data.dataOnEntities[MBTRI][ff].getDiffHdivN(base)(gg,diff_col) = diffN_face_bubble[ff](gg,dd);
        }
      }
    }
  }

  //volume
  int col = 0,diff_col = 0;
  data.dataOnEntities[MBTET][0].getHdivN(base).resize(nb_gauss_pts,3*NBVOLUMETET_HDIV_AINSWORTH_COLE(volume_order),false);
  data.dataOnEntities[MBTET][0].getDiffHdivN(base).resize(nb_gauss_pts,9*NBVOLUMETET_HDIV_AINSWORTH_COLE(volume_order),false);
  for(int oo = 0;oo<volume_order;oo++) {
    for(int ee = 0;ee<6;ee++) {
      //values
      for(int dd = 3*NBVOLUMETET_EDGE_HDIV_AINSWORTH_COLE(oo);dd<3*NBVOLUMETET_EDGE_HDIV_AINSWORTH_COLE(oo+1);dd++,col++) {
        for(int gg = 0;gg!=nb_gauss_pts;gg++) {
          data.dataOnEntities[MBTET][0].getHdivN(base)(gg,col) = N_volume_edge[ee](gg,dd);
        }
      }
      //direvatives
      for(int dd = 9*NBVOLUMETET_EDGE_HDIV_AINSWORTH_COLE(oo);dd<9*NBVOLUMETET_EDGE_HDIV_AINSWORTH_COLE(oo+1);dd++,diff_col++) {
        for(int gg = 0;gg!=nb_gauss_pts;gg++) {
          data.dataOnEntities[MBTET][0].getDiffHdivN(base)(gg,diff_col) = diffN_volume_edge[ee](gg,dd);
        }
      }
    }
    for(int ff = 0;ff<4;ff++) {
      //values
      for(int dd = 3*NBVOLUMETET_FACE_HDIV_AINSWORTH_COLE(oo);dd<3*NBVOLUMETET_FACE_HDIV_AINSWORTH_COLE(oo+1);dd++,col++) {
        for(int gg = 0;gg!=nb_gauss_pts;gg++) {
          data.dataOnEntities[MBTET][0].getHdivN(base)(gg,col) = N_volume_face[ff](gg,dd);
        }
      }
      //direvatives
      for(int dd = 9*NBVOLUMETET_FACE_HDIV_AINSWORTH_COLE(oo);dd<9*NBVOLUMETET_FACE_HDIV_AINSWORTH_COLE(oo+1);dd++,diff_col++) {
        for(int gg = 0;gg!=nb_gauss_pts;gg++) {
          data.dataOnEntities[MBTET][0].getDiffHdivN(base)(gg,diff_col) = diffN_volume_face[ff](gg,dd);
        }
      }
    }
    //values
    for(int dd = 3*NBVOLUMETET_VOLUME_HDIV_AINSWORTH_COLE(oo);dd<3*NBVOLUMETET_VOLUME_HDIV_AINSWORTH_COLE(oo+1);dd++,col++) {
      for(int gg = 0;gg!=nb_gauss_pts;gg++) {
        data.dataOnEntities[MBTET][0].getHdivN(base)(gg,col) = N_volume_bubble(gg,dd);
      }
    }
    //direvatives
    for(int dd = 9*NBVOLUMETET_VOLUME_HDIV_AINSWORTH_COLE(oo);dd<9*NBVOLUMETET_VOLUME_HDIV_AINSWORTH_COLE(oo+1);dd++,diff_col++) {
      for(int gg = 0;gg!=nb_gauss_pts;gg++) {
        data.dataOnEntities[MBTET][0].getDiffHdivN(base)(gg,diff_col) = diffN_volume_bubble(gg,dd);
      }
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode TetPolynomialBase::getValueHCurl(
  ublas::matrix<double> &pts
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  DataForcesAndSurcesCore& data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(
    int p,double s,double *diff_s,double *L,double *diffL,const int dim
  ) = cTx->basePolynomials;

  int nb_gauss_pts = pts.size2();

  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Not yet implemented (You can do it)");

  PetscFunctionReturn(0);
}

PetscErrorCode TetPolynomialBase::getValue(
  ublas::matrix<double> &pts,
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

  if(pts.size1()<3) {
    SETERRQ(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "Wrong dimension of pts, should be at least 3 rows with coordinates"
    );
  }

  const FieldApproximationBase base = cTx->bAse;
  DataForcesAndSurcesCore& data = cTx->dAta;
  if(cTx->copyNodeBase==LASTBASE) {
    data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts,4,false);
    ierr = ShapeMBTET(
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &pts(0,0),
      &pts(1,0),
      &pts(2,0),
      nb_gauss_pts
    ); CHKERRQ(ierr);
  } else {
    data.dataOnEntities[MBVERTEX][0].getNSharedPtr(base) = data.dataOnEntities[MBVERTEX][0].getNSharedPtr(cTx->copyNodeBase);
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
  // if(cTx->sPace==H1) {
  //   MatrixDouble diffN(nb_gauss_pts,12);
  //   for(int gg = 0;gg<nb_gauss_pts;gg++) {
  //     for(int nn = 0;nn<4;nn++) {
  //       for(int dd = 0;dd<3;dd++) {
  //         diffN(gg,nn*3+dd) = data.dataOnEntities[MBVERTEX][0].getDiffN(base)(nn,dd);
  //       }
  //     }
  //   }
  //   data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(diffN.size1(),diffN.size2(),false);
  //   data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().swap(diffN.data());
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
