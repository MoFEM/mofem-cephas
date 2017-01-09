/** \file TetPolynomialBase.cpp
\brief Implementation of hierarchical bases on tetrahedral

A l2, h1, h-div and h-curl spaces are implemented.

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
#include <TetPolynomialBase.hpp>

#include <Hcurl.hpp>
#include <Hdiv.hpp>

PetscErrorCode TetPolynomialBase::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_TET_BASE_FUNCTION) {
    *iface = static_cast<TetPolynomialBase*>(this);
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
  ) = cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  int sense[6],order[6];
  if(data.spacesOnEntities[MBEDGE].test(H1)) {
    //edges
    if(data.dataOnEntities[MBEDGE].size()!=6) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    double *h1_edge_n[6],*diff_h1_egde_n[6];
    for(int ee = 0;ee<6;ee++) {
      if(data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      order[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();
      int nb_dofs = NBEDGE_H1(data.dataOnEntities[MBEDGE][ee].getDataOrder());
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts,nb_dofs,false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts,3*nb_dofs,false);
      h1_edge_n[ee] = &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      diff_h1_egde_n[ee] = &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }
    ierr = H1_EdgeShapeFunctions_MBTET(
      sense,order,
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      h1_edge_n,diff_h1_egde_n,nb_gauss_pts,base_polynomials
    ); CHKERRQ(ierr);
  }

  if(data.spacesOnEntities[MBTRI].test(H1)) {

    //faces
    if(data.dataOnEntities[MBTRI].size()!=4) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    double *h1_face_n[4],*diff_h1_face_n[4];
    for(int ff = 0;ff<4;ff++) {
      if(data.dataOnEntities[MBTRI][ff].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      int nb_dofs = NBFACETRI_H1(data.dataOnEntities[MBTRI][ff].getDataOrder());
      order[ff] = data.dataOnEntities[MBTRI][ff].getDataOrder();
      data.dataOnEntities[MBTRI][ff].getN(base).resize(nb_gauss_pts,nb_dofs,false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(nb_gauss_pts,3*nb_dofs,false);
      h1_face_n[ff] = &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin();
      diff_h1_face_n[ff] = &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin();
    }
    if(data.facesNodes.size1() != 4) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    if(data.facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    ierr = H1_FaceShapeFunctions_MBTET(
      &*data.facesNodes.data().begin(),
      order,
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      h1_face_n,
      diff_h1_face_n,
      nb_gauss_pts,
      base_polynomials
    ); CHKERRQ(ierr);

    // std::cerr << "Aaaaa2\n";
    // std::cerr << data.dataOnEntities[MBVERTEX][0].getN(base) << std::endl;
    // std::cerr << ApproximationBaseNames[base] << std::endl;

  }

  if(data.spacesOnEntities[MBTET].test(H1)) {

    //volume
    int order = data.dataOnEntities[MBTET][0].getDataOrder();
    int nb_vol_dofs = NBVOLUMETET_H1(order);
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
  ) = cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  data.dataOnEntities[MBTET][0].getN(base).resize(
    nb_gauss_pts,NBVOLUMETET_L2(data.dataOnEntities[MBTET][0].getDataOrder()),false
  );
  data.dataOnEntities[MBTET][0].getDiffN(base).resize(
    nb_gauss_pts,3*NBVOLUMETET_L2(data.dataOnEntities[MBTET][0].getDataOrder()),false
  );

  ierr = L2_ShapeFunctions_MBTET(
    data.dataOnEntities[MBTET][0].getDataOrder(),
    &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
    &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
    &*data.dataOnEntities[MBTET][0].getN(base).data().begin(),
    &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin(),
    nb_gauss_pts,
    base_polynomials
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
  ) = cTx->basePolynomialsType0;

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
      N_face_edge(ff,ee).resize(nb_gauss_pts,3*NBFACETRI_EDGE_HDIV(faces_order[ff]),false);
      diffN_face_edge(ff,ee).resize(nb_gauss_pts,9*NBFACETRI_EDGE_HDIV(faces_order[ff]),false);
      phi_f_e[ff][ee] = &((N_face_edge(ff,ee))(0,0));
      diff_phi_f_e[ff][ee] = &((diffN_face_edge(ff,ee))(0,0));
    }
    N_face_bubble[ff].resize(nb_gauss_pts,3*NBFACETRI_FACE_HDIV(faces_order[ff]),false);
    diffN_face_bubble[ff].resize(nb_gauss_pts,9*NBFACETRI_FACE_HDIV(faces_order[ff]),false);
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

  N_volume_edge.resize(6,false);
  diffN_volume_edge.resize(6,false);
  for(int ee = 0;ee<6;ee++) {
    N_volume_edge[ee].resize(nb_gauss_pts,3*NBVOLUMETET_EDGE_HDIV(volume_order),false);
    diffN_volume_edge[ee].resize(nb_gauss_pts,9*NBVOLUMETET_EDGE_HDIV(volume_order),false);
    phi_v_e[ee] = &*(N_volume_edge[ee].data().begin());
    diff_phi_v_e[ee] = &*(diffN_volume_edge[ee].data().begin());
  }
  ierr = Hdiv_EdgeBasedVolumeShapeFunctions_MBTET(
    volume_order,&data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
    phi_v_e,diff_phi_v_e,nb_gauss_pts,
    base_polynomials
  ); CHKERRQ(ierr);

  N_volume_face.resize(4,false);
  diffN_volume_face.resize(4,false);
  for(int ff = 0;ff<4;ff++) {
    N_volume_face[ff].resize(nb_gauss_pts,3*NBVOLUMETET_FACE_HDIV(volume_order),false);
    diffN_volume_face[ff].resize(nb_gauss_pts,9*NBVOLUMETET_FACE_HDIV(volume_order),false);
    phi_v_f[ff] = &*(N_volume_face[ff].data().begin());
    diff_phi_v_f[ff] = &*(diffN_volume_face[ff].data().begin());
  }
  ierr = Hdiv_FaceBasedVolumeShapeFunctions_MBTET(
    volume_order,&data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
    phi_v_f,diff_phi_v_f,nb_gauss_pts,
    base_polynomials
  ); CHKERRQ(ierr);

  N_volume_bubble.resize(nb_gauss_pts,3*NBVOLUMETET_VOLUME_HDIV(volume_order),false);
  diffN_volume_bubble.resize(nb_gauss_pts,9*NBVOLUMETET_VOLUME_HDIV(volume_order),false);
  phi_v = &*(N_volume_bubble.data().begin());
  diff_phi_v = &*(diffN_volume_bubble.data().begin());
  ierr = Hdiv_VolumeBubbleShapeFunctions_MBTET(
    volume_order,&data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
    phi_v,diff_phi_v,nb_gauss_pts,
    base_polynomials
  ); CHKERRQ(ierr);

  // Set shape functions into data structure Shape functions hast to be put
  // in arrays in order which guarantee hierarchical series of degrees of
  // freedom, i.e. in other words dofs form sub-entities has to be group
  // by order.

  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;

  //faces
  if(data.dataOnEntities[MBTRI].size()!=4) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  for(int ff = 0;ff!=4;ff++) {
    data.dataOnEntities[MBTRI][ff].getHdivN(base).resize(nb_gauss_pts,3*NBFACETRI_HDIV(faces_order[ff]),false);
    data.dataOnEntities[MBTRI][ff].getDiffHdivN(base).resize(nb_gauss_pts,9*NBFACETRI_HDIV(faces_order[ff]),false);
    if(NBFACETRI_HDIV(faces_order[ff])==0) continue;
    // face
    double *base_ptr = &*data.dataOnEntities[MBTRI][ff].getHdivN(base).data().begin();
    FTensor::Tensor1<double*,3> t_base(base_ptr,&base_ptr[HDIV1],&base_ptr[HDIV2],3);
    double *diff_base_ptr = &*data.dataOnEntities[MBTRI][ff].getDiffHdivN(base).data().begin();
    FTensor::Tensor2<double*,3,3> t_diff_base(
      &diff_base_ptr[HDIV0_0],&diff_base_ptr[HDIV0_1],&diff_base_ptr[HDIV0_2],
      &diff_base_ptr[HDIV1_0],&diff_base_ptr[HDIV1_1],&diff_base_ptr[HDIV1_2],
      &diff_base_ptr[HDIV2_0],&diff_base_ptr[HDIV2_1],&diff_base_ptr[HDIV2_2],9
    );
    // face-face
    boost::shared_ptr<FTensor::Tensor1<double*,3> > t_base_f;
    boost::shared_ptr<FTensor::Tensor2<double*,3,3> > t_diff_base_f;
    if(NBFACETRI_FACE_HDIV(faces_order[ff])>0) {
      base_ptr = phi_f[ff];
      t_base_f = boost::shared_ptr<FTensor::Tensor1<double*,3> >(
        new FTensor::Tensor1<double*,3>
        (base_ptr,&base_ptr[HDIV1],&base_ptr[HDIV2],3)
      );
      diff_base_ptr = diff_phi_f[ff];
      t_diff_base_f = boost::shared_ptr<FTensor::Tensor2<double*,3,3> >(
        new FTensor::Tensor2<double*,3,3>
        (
          &diff_base_ptr[HDIV0_0],&diff_base_ptr[HDIV0_1],&diff_base_ptr[HDIV0_2],
          &diff_base_ptr[HDIV1_0],&diff_base_ptr[HDIV1_1],&diff_base_ptr[HDIV1_2],
          &diff_base_ptr[HDIV2_0],&diff_base_ptr[HDIV2_1],&diff_base_ptr[HDIV2_2],9
        )
      );
    }
    // edge-face
    base_ptr = phi_f_e[ff][0];
    FTensor::Tensor1<double*,3> t_base_f_e0(base_ptr,&base_ptr[HDIV1],&base_ptr[HDIV2],3);
    diff_base_ptr = diff_phi_f_e[ff][0];
    FTensor::Tensor2<double*,3,3> t_diff_base_f_e0(
      &diff_base_ptr[HDIV0_0],&diff_base_ptr[HDIV0_1],&diff_base_ptr[HDIV0_2],
      &diff_base_ptr[HDIV1_0],&diff_base_ptr[HDIV1_1],&diff_base_ptr[HDIV1_2],
      &diff_base_ptr[HDIV2_0],&diff_base_ptr[HDIV2_1],&diff_base_ptr[HDIV2_2],9
    );
    base_ptr = phi_f_e[ff][1];
    FTensor::Tensor1<double*,3> t_base_f_e1(base_ptr,&base_ptr[HDIV1],&base_ptr[HDIV2],3);
    diff_base_ptr = diff_phi_f_e[ff][1];
    FTensor::Tensor2<double*,3,3> t_diff_base_f_e1(
      &diff_base_ptr[HDIV0_0],&diff_base_ptr[HDIV0_1],&diff_base_ptr[HDIV0_2],
      &diff_base_ptr[HDIV1_0],&diff_base_ptr[HDIV1_1],&diff_base_ptr[HDIV1_2],
      &diff_base_ptr[HDIV2_0],&diff_base_ptr[HDIV2_1],&diff_base_ptr[HDIV2_2],9
    );
    base_ptr = phi_f_e[ff][2];
    FTensor::Tensor1<double*,3> t_base_f_e2(base_ptr,&base_ptr[HDIV1],&base_ptr[HDIV2],3);
    diff_base_ptr = diff_phi_f_e[ff][2];
    FTensor::Tensor2<double*,3,3> t_diff_base_f_e2(
      &diff_base_ptr[HDIV0_0],&diff_base_ptr[HDIV0_1],&diff_base_ptr[HDIV0_2],
      &diff_base_ptr[HDIV1_0],&diff_base_ptr[HDIV1_1],&diff_base_ptr[HDIV1_2],
      &diff_base_ptr[HDIV2_0],&diff_base_ptr[HDIV2_1],&diff_base_ptr[HDIV2_2],9
    );
    for(int gg = 0;gg!=nb_gauss_pts;gg++) {
      for(int oo = 0;oo!=faces_order[ff];oo++) {
        for(int dd = NBFACETRI_EDGE_HDIV(oo);dd!=NBFACETRI_EDGE_HDIV(oo+1);dd++) {
          t_base(i) = t_base_f_e0(i);
          ++t_base;
          ++t_base_f_e0;
          t_diff_base(i,j) = t_diff_base_f_e0(i,j);
          ++t_diff_base;
          ++t_diff_base_f_e0;
          t_base(i) = t_base_f_e1(i);
          ++t_base;
          ++t_base_f_e1;
          t_diff_base(i,j) = t_diff_base_f_e1(i,j);
          ++t_diff_base;
          ++t_diff_base_f_e1;
          t_base(i) = t_base_f_e2(i);
          ++t_base;
          ++t_base_f_e2;
          t_diff_base(i,j) = t_diff_base_f_e2(i,j);
          ++t_diff_base;
          ++t_diff_base_f_e2;
        }
        for(int dd = NBFACETRI_FACE_HDIV(oo);dd!=NBFACETRI_FACE_HDIV(oo+1);dd++) {
          t_base(i) = (*t_base_f)(i);
          ++t_base;
          ++(*t_base_f);
          t_diff_base(i,j) = (*t_diff_base_f)(i,j);
          ++t_diff_base;
          ++(*t_diff_base_f);
        }
      }
    }
  }

  //volume
  data.dataOnEntities[MBTET][0].getHdivN(base).resize(nb_gauss_pts,3*NBVOLUMETET_HDIV(volume_order),false);
  data.dataOnEntities[MBTET][0].getDiffHdivN(base).resize(nb_gauss_pts,9*NBVOLUMETET_HDIV(volume_order),false);
  if(NBVOLUMETET_HDIV(volume_order)>0) {
    double *base_ptr = &*data.dataOnEntities[MBTET][0].getHdivN(base).data().begin();
    FTensor::Tensor1<double*,3> t_base(base_ptr,&base_ptr[HDIV1],&base_ptr[HDIV2],3);
    double *diff_base_ptr = &*data.dataOnEntities[MBTET][0].getDiffHdivN(base).data().begin();
    FTensor::Tensor2<double*,3,3> t_diff_base(
      &diff_base_ptr[HDIV0_0],&diff_base_ptr[HDIV0_1],&diff_base_ptr[HDIV0_2],
      &diff_base_ptr[HDIV1_0],&diff_base_ptr[HDIV1_1],&diff_base_ptr[HDIV1_2],
      &diff_base_ptr[HDIV2_0],&diff_base_ptr[HDIV2_1],&diff_base_ptr[HDIV2_2],9
    );
    // edges
    std::vector<FTensor::Tensor1<double*,3> > t_base_v_e;
    t_base_v_e.reserve(6);
    std::vector<FTensor::Tensor2<double*,3,3> > t_diff_base_v_e;
    t_diff_base_v_e.reserve(6);
    for(int ee = 0;ee!=6;ee++) {
      base_ptr = phi_v_e[ee];
      diff_base_ptr = diff_phi_v_e[ee];
      t_base_v_e.push_back(
        FTensor::Tensor1<double*,3>(base_ptr,&base_ptr[HDIV1],&base_ptr[HDIV2],3)
      );
      t_diff_base_v_e.push_back(
        FTensor::Tensor2<double*,3,3>(
          &diff_base_ptr[HDIV0_0],&diff_base_ptr[HDIV0_1],&diff_base_ptr[HDIV0_2],
          &diff_base_ptr[HDIV1_0],&diff_base_ptr[HDIV1_1],&diff_base_ptr[HDIV1_2],
          &diff_base_ptr[HDIV2_0],&diff_base_ptr[HDIV2_1],&diff_base_ptr[HDIV2_2],9
        )
      );
    }
    // faces
    std::vector<FTensor::Tensor1<double*,3> > t_base_v_f;
    t_base_v_f.reserve(4);
    std::vector<FTensor::Tensor2<double*,3,3> > t_diff_base_v_f;
    t_diff_base_v_f.reserve(4);
    if(NBVOLUMETET_FACE_HDIV(volume_order)>0) {
      for(int ff = 0;ff!=4;ff++) {
        base_ptr = phi_v_f[ff];
        diff_base_ptr = diff_phi_v_f[ff];
        t_base_v_f.push_back(
          FTensor::Tensor1<double*,3>(base_ptr,&base_ptr[HDIV1],&base_ptr[HDIV2],3)
        );
        t_diff_base_v_f.push_back(
          FTensor::Tensor2<double*,3,3>(
            &diff_base_ptr[HDIV0_0],&diff_base_ptr[HDIV0_1],&diff_base_ptr[HDIV0_2],
            &diff_base_ptr[HDIV1_0],&diff_base_ptr[HDIV1_1],&diff_base_ptr[HDIV1_2],
            &diff_base_ptr[HDIV2_0],&diff_base_ptr[HDIV2_1],&diff_base_ptr[HDIV2_2],9
          )
        );
      }
    }
    boost::shared_ptr<FTensor::Tensor1<double*,3> > t_base_v;
    boost::shared_ptr<FTensor::Tensor2<double*,3,3> > t_diff_base_v;
    if(NBVOLUMETET_VOLUME_HDIV(volume_order)>0) {
      base_ptr = phi_v;
      t_base_v = boost::shared_ptr<FTensor::Tensor1<double*,3> >(
        new FTensor::Tensor1<double*,3>
        (base_ptr,&base_ptr[HDIV1],&base_ptr[HDIV2],3)
      );
      diff_base_ptr = diff_phi_v;
      t_diff_base_v = boost::shared_ptr<FTensor::Tensor2<double*,3,3> >(
        new FTensor::Tensor2<double*,3,3>
        (
          &diff_base_ptr[HDIV0_0],&diff_base_ptr[HDIV0_1],&diff_base_ptr[HDIV0_2],
          &diff_base_ptr[HDIV1_0],&diff_base_ptr[HDIV1_1],&diff_base_ptr[HDIV1_2],
          &diff_base_ptr[HDIV2_0],&diff_base_ptr[HDIV2_1],&diff_base_ptr[HDIV2_2],9
        )
      );
    }
    for(int gg = 0;gg!=nb_gauss_pts;gg++) {
      for(int oo = 0;oo<volume_order;oo++) {
        for(int dd = NBVOLUMETET_EDGE_HDIV(oo);dd<NBVOLUMETET_EDGE_HDIV(oo+1);dd++) {
          for(int ee = 0;ee<6;ee++) {
            t_base(i) = t_base_v_e[ee](i);
            ++t_base;
            ++t_base_v_e[ee];
            t_diff_base(i,j) = t_diff_base_v_e[ee](i,j);
            ++t_diff_base;
            ++t_diff_base_v_e[ee];
          }
        }
        for(int dd = NBVOLUMETET_FACE_HDIV(oo);dd<NBVOLUMETET_FACE_HDIV(oo+1);dd++) {
          for(int ff = 0;ff<4;ff++) {
            t_base(i) = t_base_v_f[ff](i);
            ++t_base;
            ++t_base_v_f[ff];
            t_diff_base(i,j) = t_diff_base_v_f[ff](i,j);
            ++t_diff_base;
            ++t_diff_base_v_f[ff];
          }
        }
        for(int dd = NBVOLUMETET_VOLUME_HDIV(oo);dd<NBVOLUMETET_VOLUME_HDIV(oo+1);dd++) {
          t_base(i) = (*t_base_v)(i);
          ++t_base;
          ++(*t_base_v);
          t_diff_base(i,j) = (*t_diff_base_v)(i,j);
          ++t_diff_base;
          ++(*t_diff_base_v);
        }
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

  try {

  DataForcesAndSurcesCore& data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(
    int p,double s,double *diff_s,double *L,double *diffL,const int dim
  ) = cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  // edges
  if(data.spacesOnEntities[MBEDGE].test(HCURL)) {
    int sense[6],order[6];
    if(data.dataOnEntities[MBEDGE].size()!=6) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    double *hcurl_edge_n[6],*diff_hcurl_edge_n[6];
    for(int ee = 0;ee!=6;ee++) {
      if(data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      order[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();
      int nb_dofs = NBEDGE_HCURL(data.dataOnEntities[MBEDGE][ee].getDataOrder());
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts,3*nb_dofs,false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts,9*nb_dofs,false);
      hcurl_edge_n[ee] = &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      diff_hcurl_edge_n[ee] = &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }
    ierr = Hcurl_EdgeBaseFunctions_MBTET(
      sense,
      order,
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      hcurl_edge_n,
      diff_hcurl_edge_n,
      nb_gauss_pts,
      base_polynomials
    ); CHKERRQ(ierr);
  } else {
    for(int ee = 0;ee!=6;ee++) {
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts,0,false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts,0,false);
    }
  }

  // triangles
  if(data.spacesOnEntities[MBTRI].test(HCURL)) {
    int order[4];
    //faces
    if(data.dataOnEntities[MBTRI].size()!=4) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    double *hcurl_base_n[4],*diff_hcurl_base_n[4];
    for(int ff = 0;ff!=4;ff++) {
      if(data.dataOnEntities[MBTRI][ff].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      order[ff] = data.dataOnEntities[MBTRI][ff].getDataOrder();
      int nb_dofs = NBFACETRI_HCURL(order[ff] );
      data.dataOnEntities[MBTRI][ff].getN(base).resize(nb_gauss_pts,3*nb_dofs,false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(nb_gauss_pts,9*nb_dofs,false);
      hcurl_base_n[ff] = &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin();
      diff_hcurl_base_n[ff] = &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin();
    }
    if(data.facesNodes.size1() != 4) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    if(data.facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    ierr = Hcurl_FaceFunctions_MBTET(
      &*data.facesNodes.data().begin(),
      order,
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      hcurl_base_n,
      diff_hcurl_base_n,
      nb_gauss_pts,
      base_polynomials
    ); CHKERRQ(ierr);
  } else {
    for(int ff = 0;ff!=4;ff++) {
      data.dataOnEntities[MBTRI][ff].getN(base).resize(nb_gauss_pts,0,false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(nb_gauss_pts,0,false);
    }
  }

  if(data.spacesOnEntities[MBTET].test(HCURL)) {

    //volume
    int order = data.dataOnEntities[MBTET][0].getDataOrder();
    int nb_vol_dofs = NBVOLUMETET_HCURL(order);
    data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts,3*nb_vol_dofs,false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts,9*nb_vol_dofs,false);
    // cerr << data.dataOnEntities[MBVERTEX][0].getDiffN(base) << endl;
    ierr = Hcurl_VolumeFunctions_MBTET(
      data.dataOnEntities[MBTET][0].getDataOrder(),
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      &*data.dataOnEntities[MBTET][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin(),
      nb_gauss_pts,
      base_polynomials
    ); CHKERRQ(ierr);

  } else {
    data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts,0,false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts,0,false);
  }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

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
  if(data.dataOnEntities[MBVERTEX][0].getN(base).size1()!=(unsigned int)nb_gauss_pts) {
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
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Unknown space");
  }

  PetscFunctionReturn(0);
}
