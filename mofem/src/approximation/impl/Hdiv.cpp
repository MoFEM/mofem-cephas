/** \file Hdiv.cpp

  \brief Implementation of H-curl base

  Based on Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by Mark Ainsworth and Joe Coyle
  Shape functions for MBTRI/MBTET and HCurl space

*/

#include <petscsys.h>
#include <FTensor.hpp>
#include <h1_hdiv_hcurl_l2.h>
#include <Hdiv.hpp>
#include <definitions.h>

using namespace MoFEM;

PetscErrorCode MoFEM::Hdiv_EdgeFaceShapeFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *phi_f_e[4][3],double *diff_phi_f_e[4][3],int gdim,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  for(int ff = 0;ff<4;ff++) {
    if(diff_phi_f_e!=NULL) {
      ierr = Hdiv_EdgeFaceShapeFunctions_MBTET_ON_FACE(
        &faces_nodes[3*ff],p[ff],N,diffN,phi_f_e[ff],diff_phi_f_e[ff],gdim,4,base_polynomials
      ); CHKERRQ(ierr);
    } else {
      ierr = Hdiv_EdgeFaceShapeFunctions_MBTET_ON_FACE(
        &faces_nodes[3*ff],p[ff],N,diffN,phi_f_e[ff],NULL,gdim,4,base_polynomials
      ); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MoFEM::Hdiv_EdgeFaceShapeFunctions_MBTET_ON_FACE(
  int *faces_nodes,int p,double *N,double *diffN,double *phi_f_e[3],double *diff_phi_f_e[3],int gdim,int nb,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {

  const int face_edges_nodes[3][2] = { {0,1}, {1,2}, {2,0} };
  const int face_oposite_edges_node[] = { 2, 0, 1 };
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(p<1) PetscFunctionReturn(0);

  FTensor::Tensor1<double,3> t_edge_cross[3];
  FTensor::Tensor1<double,3> t_node_diff_ksi[4];
  FTensor::Tensor1<double,3> t_diff_ksi0i[3];
  if(diffN!=NULL) {
    t_node_diff_ksi[0] = FTensor::Tensor1<double,3>(diffN[0],diffN[ 1],diffN[ 2]);
    t_node_diff_ksi[1] = FTensor::Tensor1<double,3>(diffN[3],diffN[ 4],diffN[ 5]);
    t_node_diff_ksi[2] = FTensor::Tensor1<double,3>(diffN[6],diffN[ 7],diffN[ 8]);
    t_node_diff_ksi[3] = FTensor::Tensor1<double,3>(diffN[9],diffN[10],diffN[11]);
    for(int ee = 0;ee<3;ee++) {
      const int n0 = faces_nodes[face_edges_nodes[ee][0]];
      const int n1 = faces_nodes[face_edges_nodes[ee][1]];
      t_diff_ksi0i[ee](i) =
      t_node_diff_ksi[n1](i)-t_node_diff_ksi[n0](i);
      t_edge_cross[ee](0) =
      t_node_diff_ksi[n0](1)*t_node_diff_ksi[n1](2)-
      t_node_diff_ksi[n0](2)*t_node_diff_ksi[n1](1);
      t_edge_cross[ee](1) =
      t_node_diff_ksi[n0](2)*t_node_diff_ksi[n1](0)-
      t_node_diff_ksi[n0](0)*t_node_diff_ksi[n1](2);
      t_edge_cross[ee](2) =
      t_node_diff_ksi[n0](0)*t_node_diff_ksi[n1](1)-
      t_node_diff_ksi[n0](1)*t_node_diff_ksi[n1](0);
    }
  } else {
    for(int ee = 0;ee<3;ee++) {
      t_edge_cross[ee](0) = 1;
      t_edge_cross[ee](1) = 0;
      t_edge_cross[ee](2) = 0;
    }
  }
  double psi_l[p+1],diff_psi_l[3*(p+1)];
  boost::shared_ptr<FTensor::Tensor2<double*,3,3> > t_diff_phi_f_e_ptr;

  for(int ee = 0;ee!=3;ee++) {
    const int i0 = faces_nodes[face_edges_nodes[ee][0]];
    const int i1 = faces_nodes[face_edges_nodes[ee][1]];
    const int iO = faces_nodes[face_oposite_edges_node[ee]];
    FTensor::Tensor1<double*,3> t_psi_f_e(
      &phi_f_e[ee][0],&phi_f_e[ee][1],&phi_f_e[ee][2],3
    );
    if(diff_phi_f_e) {
      t_diff_phi_f_e_ptr = boost::shared_ptr<FTensor::Tensor2<double*,3,3> >(
        new FTensor::Tensor2<double*,3,3>(
          &diff_phi_f_e[ee][HDIV0_0],&diff_phi_f_e[ee][HDIV0_1],&diff_phi_f_e[ee][HDIV0_2],
          &diff_phi_f_e[ee][HDIV1_0],&diff_phi_f_e[ee][HDIV1_1],&diff_phi_f_e[ee][HDIV1_2],
          &diff_phi_f_e[ee][HDIV2_0],&diff_phi_f_e[ee][HDIV2_1],&diff_phi_f_e[ee][HDIV2_2],9
        )
      );
    }
    for(int ii = 0;ii!=gdim;ii++) {
      const int node_shift = ii*nb;
      const double n0 = N[node_shift+i0];
      const double n1 = N[node_shift+i1];
      const double lambda = N[node_shift+iO];
      const double ksi0i = n1 - n0;
      if(diff_phi_f_e) {
        ierr = base_polynomials(p,ksi0i,&t_diff_ksi0i[ee](0),psi_l,diff_psi_l,3); CHKERRQ(ierr);
      } else {
        ierr = base_polynomials(p,ksi0i,NULL,psi_l,NULL,3); CHKERRQ(ierr);
      }
      FTensor::Tensor0<double*> t_psi_l(&psi_l[0]);
      FTensor::Tensor1<double*,3> t_diff_psi_l(
        &diff_psi_l[0],&diff_psi_l[p+1],&diff_psi_l[2*p+2],1
      );
      for(int l = 0;l<=p-1;l++) {
        t_psi_f_e(i) = lambda*t_psi_l*t_edge_cross[ee](i);
        if(t_diff_phi_f_e_ptr) {
          (*t_diff_phi_f_e_ptr)(i,j) =
          (t_node_diff_ksi[iO](j)*t_psi_l+lambda*t_diff_psi_l(j))*t_edge_cross[ee](i);
          ++t_diff_psi_l;
          ++(*t_diff_phi_f_e_ptr);
        }
        ++t_psi_f_e;
        ++t_psi_l;
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MoFEM::Hdiv_FaceBubbleShapeFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *phi_f[],double *diff_phi_f[],int gdim,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  for(int ff = 0;ff<4;ff++) {
    double *diff;
    if(diff_phi_f!=NULL) {
      diff = diff_phi_f[ff];
    } else {
      diff = NULL;
    }
    ierr = Hdiv_FaceBubbleShapeFunctions_MBTET_ON_FACE(
      &faces_nodes[3*ff],p[ff],N,diffN,phi_f[ff],diff,gdim,4,base_polynomials
    ); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MoFEM::Hdiv_FaceBubbleShapeFunctions_MBTET_ON_FACE(
  int *face_nodes,int p,double *N,double *diffN,double *phi_f,double *diff_phi_f,int gdim,int nb,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(p<3) PetscFunctionReturn(0);

  const int vert_i = face_nodes[1];
  const int vert_j = face_nodes[2];
  const int i0 = face_nodes[0];
  FTensor::Tensor1<double,3> t_cross;
  FTensor::Tensor1<double,3> t_node_diff_ksi[4];
  FTensor::Tensor1<double,3> t_diff_ksi0i;
  FTensor::Tensor1<double,3> t_diff_ksi0j;

  if(diffN) {
    t_node_diff_ksi[0] = FTensor::Tensor1<double,3>(diffN[0],diffN[ 1],diffN[ 2]);
    t_node_diff_ksi[1] = FTensor::Tensor1<double,3>(diffN[3],diffN[ 4],diffN[ 5]);
    t_node_diff_ksi[2] = FTensor::Tensor1<double,3>(diffN[6],diffN[ 7],diffN[ 8]);
    t_node_diff_ksi[3] = FTensor::Tensor1<double,3>(diffN[9],diffN[10],diffN[11]);
    t_diff_ksi0i(i) =
    t_node_diff_ksi[vert_i](i)-t_node_diff_ksi[i0](i);
    t_diff_ksi0j(i) =
    t_node_diff_ksi[vert_j](i)-t_node_diff_ksi[i0](i);
    t_cross(0) =
    t_node_diff_ksi[vert_i](1)*t_node_diff_ksi[vert_j](2)-
    t_node_diff_ksi[vert_i](2)*t_node_diff_ksi[vert_j](1);
    t_cross(1) =
    t_node_diff_ksi[vert_i](2)*t_node_diff_ksi[vert_j](0)-
    t_node_diff_ksi[vert_i](0)*t_node_diff_ksi[vert_j](2);
    t_cross(2) =
    t_node_diff_ksi[vert_i](0)*t_node_diff_ksi[vert_j](1)-
    t_node_diff_ksi[vert_i](1)*t_node_diff_ksi[vert_j](0);
  } else {
    t_cross(0) = 1;
    t_cross(1) = 0;
    t_cross(2) = 0;
  }

  double psi_l[p+1],diff_psi_l[3*(p+1)];
  double psi_m[p+1],diff_psi_m[3*(p+1)];
  FTensor::Tensor1<double,3> t_diff_beta_0ij;

  FTensor::Tensor1<double*,3> t_psi_f(
    &phi_f[0],&phi_f[1],&phi_f[2],3
  );

  //FIXME: can be socped_ptr
  boost::shared_ptr<FTensor::Tensor2<double*,3,3> > t_diff_phi_f_ptr;
  if(diff_phi_f) {
    t_diff_phi_f_ptr = boost::shared_ptr<FTensor::Tensor2<double*,3,3> >(
      new FTensor::Tensor2<double*,3,3>(
        &diff_phi_f[HDIV0_0],&diff_phi_f[HDIV0_1],&diff_phi_f[HDIV0_2],
        &diff_phi_f[HDIV1_0],&diff_phi_f[HDIV1_1],&diff_phi_f[HDIV1_2],
        &diff_phi_f[HDIV2_0],&diff_phi_f[HDIV2_1],&diff_phi_f[HDIV2_2],9
      )
    );
  }

  for(int ii = 0;ii<gdim;ii++) {
    int node_shift = ii*nb;
    const double ni = N[ node_shift+vert_i ];
    const double nj = N[ node_shift+vert_j ];
    const double n0 = N[ node_shift+i0 ];
    const double ksi0i = ni - n0;
    const double ksi0j = nj - n0;
    double beta_0ij =  n0*ni*nj;
    if(diff_phi_f) {
      t_diff_beta_0ij(i) =
      (ni*nj)*t_node_diff_ksi[i0](i)+
      (n0*nj)*t_node_diff_ksi[vert_i](i)+
      (n0*ni)*t_node_diff_ksi[vert_j](i);
      ierr = base_polynomials(p,ksi0i,&t_diff_ksi0i(0),psi_l,diff_psi_l,3); CHKERRQ(ierr);
      ierr = base_polynomials(p,ksi0j,&t_diff_ksi0j(0),psi_m,diff_psi_m,3); CHKERRQ(ierr);
    } else {
      ierr = base_polynomials(p,ksi0i,NULL,psi_l,NULL,3); CHKERRQ(ierr);
      ierr = base_polynomials(p,ksi0j,NULL,psi_m,NULL,3); CHKERRQ(ierr);
    }

    int jj = 0;
    int oo = 0;
    for(;oo<=p-3;oo++) {
      FTensor::Tensor0<double*> t_psi_l(&psi_l[0]);
      FTensor::Tensor1<double*,3> t_diff_psi_l(
        diff_psi_l,&diff_psi_l[p+1],&diff_psi_l[2*p+2],1
      );
      for(int l = 0;l<=oo;l++) {
        int m = oo - l;
        if(m>=0) {
          FTensor::Tensor1<double,3> t_diff_psi_m(
            diff_psi_m[m],diff_psi_m[p+1+m],diff_psi_m[2*p+2+m]
          );
          t_psi_f(i) = (beta_0ij*t_psi_l*psi_m[m])*t_cross(i);
          if(diff_phi_f) {
            (*t_diff_phi_f_ptr)(i,j) = (
              (t_psi_l*psi_m[m])*t_diff_beta_0ij(j)+
              (beta_0ij*psi_m[m])*t_diff_psi_l(j)+
              (beta_0ij*t_psi_l)*t_diff_psi_m(j)
            )*t_cross(i);
            ++(*t_diff_phi_f_ptr);
          }
          ++t_psi_f;
          jj++;
        }
      }
      ++t_psi_l;
    }
    if(jj!=NBFACETRI_FACE_HDIV(p)) {
      SETERRQ2(
        PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong order %d != %d",
        jj,NBFACETRI_FACE_HDIV(p)
      );
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MoFEM::Hdiv_EdgeBasedVolumeShapeFunctions_MBTET(
  int p,
  double *N,
  double *diffN,
  double *phi_v_e[6],
  double *diff_phi_v_e[6],
  int gdim,
  PetscErrorCode (*base_polynomials)(
    int p,double s,double *diff_s,double *L,double *diffL,const int dim
  )
) {
  PetscErrorCode ierr;
  const int edges_nodes[6][2] = {
    {0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}
  };

  PetscFunctionBegin;
  if(p<2) PetscFunctionReturn(0);
  if(diffN == NULL) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }

  FTensor::Tensor1<double,3> t_coords[4] = {
    FTensor::Tensor1<double,3>(0,0,0),
    FTensor::Tensor1<double,3>(1,0,0),
    FTensor::Tensor1<double,3>(0,1,0),
    FTensor::Tensor1<double,3>(0,0,1)
  };
  FTensor::Tensor1<double*,3> t_node_diff_ksi[4] = {
    FTensor::Tensor1<double*,3>(&diffN[0],&diffN[ 1],&diffN[ 2]),
    FTensor::Tensor1<double*,3>(&diffN[3],&diffN[ 4],&diffN[ 5]),
    FTensor::Tensor1<double*,3>(&diffN[6],&diffN[ 7],&diffN[ 8]),
    FTensor::Tensor1<double*,3>(&diffN[9],&diffN[10],&diffN[11])
  };

  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;

  FTensor::Tensor1<double,3> t_tou_e;
  FTensor::Tensor1<double,3> t_diff_ksi0i;
  FTensor::Tensor1<double,3> t_diff_beta_e;

  double psi_l[p+1];
  double diff_psi_l[3*(p+1)];

  for(int ee = 0;ee!=6;ee++) {
    t_tou_e(i) =
    t_coords[edges_nodes[ee][1]](i)-t_coords[edges_nodes[ee][0]](i);
    t_diff_ksi0i(i) =
    t_node_diff_ksi[edges_nodes[ee][1]](i)-t_node_diff_ksi[edges_nodes[ee][0]](i);
    FTensor::Tensor1<double*,3> t_psi_v_e(
      &phi_v_e[ee][0],&phi_v_e[ee][1],&phi_v_e[ee][2],3
    );
    FTensor::Tensor2<double*,3,3> t_diff_phi_v_e(
      &diff_phi_v_e[ee][HDIV0_0],&diff_phi_v_e[ee][HDIV0_1],&diff_phi_v_e[ee][HDIV0_2],
      &diff_phi_v_e[ee][HDIV1_0],&diff_phi_v_e[ee][HDIV1_1],&diff_phi_v_e[ee][HDIV1_2],
      &diff_phi_v_e[ee][HDIV2_0],&diff_phi_v_e[ee][HDIV2_1],&diff_phi_v_e[ee][HDIV2_2],9
    );
    for(int ii = 0;ii!=gdim;ii++) {
      const int node_shift = ii*4;
      const double ni = N[node_shift+edges_nodes[ee][1]];
      const double n0 = N[node_shift+edges_nodes[ee][0]];
      const double beta_e = ni*n0;
      const double ksi0i = ni-n0;
      if(diff_phi_v_e) {
        t_diff_beta_e(i) =
        ni*t_node_diff_ksi[edges_nodes[ee][0]](i)+
        t_node_diff_ksi[edges_nodes[ee][1]](i)*n0;
        ierr = base_polynomials(p,ksi0i,&t_diff_ksi0i(0),psi_l,diff_psi_l,3); CHKERRQ(ierr);
      } else {
        ierr = base_polynomials(p,ksi0i,NULL,psi_l,NULL,3); CHKERRQ(ierr);
      }
      FTensor::Tensor0<double*> t_psi_l(&psi_l[0]);
      FTensor::Tensor1<double*,3> t_diff_psi_l(
        &diff_psi_l[0],&diff_psi_l[p+1],&diff_psi_l[2*p+2],1
      );
      for(int l = 0;l<=p-2;l++) {
        t_psi_v_e(i) = (beta_e*t_psi_l)*t_tou_e(i);
        ++t_psi_v_e;
        if(diff_phi_v_e) {
          t_diff_phi_v_e(i,j) =
          (t_diff_beta_e(j)*t_psi_l+beta_e*t_diff_psi_l(j))*
          t_tou_e(i);
          ++t_diff_phi_v_e;
          ++t_diff_psi_l;
        }
        ++t_psi_l;
      }
    }
  }

  PetscFunctionReturn(0);
}


// PetscErrorCode MoFEM::Hdiv_VolumeBubbleShapeFunctions_MBTET(
//   int p,double *N,double *diffN,double *phi_v,double *diff_phi_v,int gdim,
//   PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
// ) {
//   PetscErrorCode ierr;
//   PetscFunctionBegin;
//   if(p<4) PetscFunctionReturn(0);
//
//   FTensor::Tensor1<double,3> t_coords[4] = {
//     FTensor::Tensor1<double,3>(0,0,0),
//     FTensor::Tensor1<double,3>(1,0,0),
//     FTensor::Tensor1<double,3>(0,1,0),
//     FTensor::Tensor1<double,3>(0,0,1)
//   };
//   FTensor::Tensor1<double*,3> t_node_diff_ksi[4] = {
//     FTensor::Tensor1<double*,3>(&diffN[0],&diffN[ 1],&diffN[ 2]),
//     FTensor::Tensor1<double*,3>(&diffN[3],&diffN[ 4],&diffN[ 5]),
//     FTensor::Tensor1<double*,3>(&diffN[6],&diffN[ 7],&diffN[ 8]),
//     FTensor::Tensor1<double*,3>(&diffN[9],&diffN[10],&diffN[11])
//   };
//
//   FTensor::Index<'i',3> i;
//   FTensor::Index<'j',3> j;
//
//   FTensor::Tensor1<double,3> t_diff_ksi0i;
//   FTensor::Tensor1<double,3> t_diff_ksi0j;
//   FTensor::Tensor1<double,3> t_diff_ksi0k;
//
//   t_diff_ksi0i(i) = t_node_diff_ksi[1](i)-t_node_diff_ksi[0](i);
//   t_diff_ksi0j(i) = t_node_diff_ksi[2](i)-t_node_diff_ksi[0](i);
//   t_diff_ksi0k(i) = t_node_diff_ksi[3](i)-t_node_diff_ksi[0](i);
//
//   double psi_l[p+1];
//   double diff_psi_l[3*(p+1)];
//   double psi_m[p+1];
//   double diff_psi_m[3*(p+1)];
//   double psi_n[p+1];
//   double diff_psi_n[3*(p+1
//
//   FTensor::Tensor1<double,3> t_pgi_v(phi_v,&phi_v[HDIV1],&phi_v[HDIV2],3);
//
//   FTensor::Tensor1<double,3> t_diff_beta_v;
//   for(int ii = 0;ii<gdim;ii++) {
//     const int node_shift = ii*4;
//     const double n0 = N[0];
//     const double ni = N[1];
//     const double nj = N[2];
//     const double nk = N[3];
//     const double ksi0i =ni-n0;
//     const double ksi0j =nj-n0;
//     const double ksi0k =nk-n0;
//     const double beta_v = n0*n1*n2*n3;
//     if(diff_phi_v != NULL) {
//       t_diff_beta_v(i) =
//       (n1*n2*n3)*t_node_diff_ksi[0](i)+
//       (n0*n2*n3)*t_node_diff_ksi[1](i)+
//       (n0*n1*n3)*t_node_diff_ksi[2](i)+
//       (n0*n1*n2)*t_node_diff_ksi[3](i);
//       ierr = base_polynomials(p,ksi_0i,&t_diff_ksi0i(0),psi_l,diff_psi_l,3); CHKERRQ(ierr);
//       ierr = base_polynomials(p,ksi_0j,&t_diff_ksi0j(0),psi_m,diff_psi_m,3); CHKERRQ(ierr);
//       ierr = base_polynomials(p,ksi_0k,&t_diff_ksi0k(0),psi_n,diff_psi_n,3); CHKERRQ(ierr);
//     } else {
//       ierr = base_polynomials(p,ksi_0i,NULL,psi_l,NULL,3); CHKERRQ(ierr);
//       ierr = base_polynomials(p,ksi_0j,NULL,psi_m,NULL,3); CHKERRQ(ierr);
//       ierr = base_polynomials(p,ksi_0k,NULL,psi_n,NULL,3); CHKERRQ(ierr);
//     }
//
//     FTensor::Tensor0<double*> t_psi_l(&psi_l[0]);
//     FTensor::Tensor1<double*,3> t_diff_psi_l(
//       &diff_psi_l[0],&diff_psi_l[p+1],&diff_psi_l[2*p+2],1
//     );
//     FTensor::Tensor0<double*> t_psi_m(&psi_l[0]);
//     FTensor::Tensor1<double*,3> t_diff_psi_m(
//       &diff_psi_l[0],&diff_psi_l[p+1],&diff_psi_l[2*p+2],1
//     );
//     FTensor::Tensor0<double*> t_psi_n(&psi_l[0]);
//     FTensor::Tensor1<double*,3> t_diff_psi_n(
//       &diff_psi_l[0],&diff_psi_l[p+1],&diff_psi_l[2*p+2],1
//     );
//
//
//     int jj = 0;
//     for(int oo = 0;oo<=p-4;oo++) {
//       for(int l = 0;l<=oo;l++) {
//         for(m = 0;(l+m)<=oo;m++) {
//           int n = oo - l - m;
//           if(n>=0) {
//           }
//         }
//       }
//     }
//
//     if(3*jj!=NBVOLUMETET_VOLUME_HDIV(p)) {
//       SETERRQ2(
//         PETSC_COMM_SELF,
//         MOFEM_DATA_INCONSISTENCY,
//         "wrong order %d != %d",jj,
//         NBVOLUMETET_VOLUME_HDIV(p)
//       );
//     }
//
//   }
//
//
//   // double coords[] = { 0,0,0, 1,0,0, 0,1,0, 0,0,1 };
//   // double ed[3][3];
//   // int nn = 0;
//   // for(;nn<3;nn++) {
//   //   cblas_dcopy(3,&coords[3*(nn+1)],1,ed[nn],1);
//   //   cblas_daxpy(3,-1,&coords[0],1,ed[nn],1);
//   //   // double nrm2 = cblas_dnrm2(3,ed[nn],1);
//   //   // cblas_dscal(3,1./nrm2,ed[nn],1);
//   // }
//   // int ii = 0;
//   // for(;ii<GDIM;ii++) {
//   //   int node_shift = ii*4;
//   //   double Beta_0ijk =
//   //   N[ node_shift + 0]*N[ node_shift + 1]*N[ node_shift + 2]*N[ node_shift + 3];
//   //   double ksi_0i = N[ node_shift+1 ] - N[ node_shift+0 ];
//   //   double ksi_0j = N[ node_shift+2 ] - N[ node_shift+0 ];
//   //   double ksi_0k = N[ node_shift+3 ] - N[ node_shift+0 ];
//   //   double diff_Beta_0ijk[3] = {0,0,0};
//   //   double diff_ksi_0i[3],diff_ksi_0j[3],diff_ksi_0k[3];
//   //   double Psi_l[p+1],Psi_m[p+1],Psi_n[p+1];
//   //   double diff_Psi_l[3*(p+1)],diff_Psi_m[3*(p+1)],diff_Psi_n[3*(p+1)];
//   //   if(diffPHI_v != NULL) {
//   //     int dd = 0;
//   //     for(;dd<3;dd++) {
//   //       diff_Beta_0ijk[dd] =
//   //       diffN[ 3*0+dd ]*N[ node_shift + 1]*N[ node_shift + 2]*N[ node_shift + 3]+
//   //       N[ node_shift + 0]*diffN[ 3*1+dd ]*N[ node_shift + 2]*N[ node_shift + 3]+
//   //       N[ node_shift + 0]*N[ node_shift + 1]*diffN[ 3*2+dd ]*N[ node_shift + 3]+
//   //       N[ node_shift + 0]*N[ node_shift + 1]*N[ node_shift + 2]*diffN[ 3*3+dd ];
//   //       diff_ksi_0i[dd] = diffN[ 3*1+dd ] - diffN[ 3*0+dd ];
//   //       diff_ksi_0j[dd] = diffN[ 3*2+dd ] - diffN[ 3*0+dd ];
//   //       diff_ksi_0k[dd] = diffN[ 3*3+dd ] - diffN[ 3*0+dd ];
//   //     }
//   //     ierr = base_polynomials(p,ksi_0i,diff_ksi_0i,Psi_l,diff_Psi_l,3); CHKERRQ(ierr);
//   //     ierr = base_polynomials(p,ksi_0j,diff_ksi_0j,Psi_m,diff_Psi_m,3); CHKERRQ(ierr);
//   //     ierr = base_polynomials(p,ksi_0k,diff_ksi_0k,Psi_n,diff_Psi_n,3); CHKERRQ(ierr);
//   //   } else {
//   //     ierr = base_polynomials(p,ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
//   //     ierr = base_polynomials(p,ksi_0j,NULL,Psi_m,NULL,3); CHKERRQ(ierr);
//   //     ierr = base_polynomials(p,ksi_0k,NULL,Psi_n,NULL,3); CHKERRQ(ierr);
//   //   }
//   //   int shift = ii*NBVOLUMETET_VOLUME_HDIV(p);
//   //   int jj = 0;
//   //   int oo = 0;
//   //   for(;oo<=p-4;oo++) {
//   //     int l = 0;
//   //     for(;l<=oo;l++) {
//   //       int m = 0;
//   //       for(;(l+m)<=oo;m++) {
//   //         int n = oo - l - m;
//   //         if(n>=0) {
//   //           double s = Beta_0ijk*Psi_l[l]*Psi_m[m]*Psi_n[n];
//   //           int kk = 0;
//   //           for(;kk<3;kk++) {
//   //             PHI_v[3*shift + 3*3*jj + 3*0 + kk] = s*ed[0][kk];
//   //             PHI_v[3*shift + 3*3*jj + 3*1 + kk] = s*ed[1][kk];
//   //             PHI_v[3*shift + 3*3*jj + 3*2 + kk] = s*ed[2][kk];
//   //           }
//   //           if(diffPHI_v!=NULL) {
//   //             int dd = 0;
//   //             for(;dd<3;dd++) {
//   //               double diff =
//   //               diff_Beta_0ijk[dd]*Psi_l[l]*Psi_m[m]*Psi_n[n]+
//   //               Beta_0ijk*diff_Psi_l[dd*(p+1)+l]*Psi_m[m]*Psi_n[n]+
//   //               Beta_0ijk*Psi_l[l]*diff_Psi_m[dd*(p+1)+m]*Psi_n[n]+
//   //               Beta_0ijk*Psi_l[l]*Psi_m[m]*diff_Psi_n[dd*(p+1)+n];
//   //               int kk = 0;
//   //               for(;kk<3;kk++) {
//   //                 diffPHI_v[9*shift + 3*9*jj + 9*0 + 3*dd + kk] = diff*ed[0][kk];
//   //                 diffPHI_v[9*shift + 3*9*jj + 9*1 + 3*dd + kk] = diff*ed[1][kk];
//   //                 diffPHI_v[9*shift + 3*9*jj + 9*2 + 3*dd + kk] = diff*ed[2][kk];
//   //               }
//   //             }
//   //           }
//   //           jj++;
//   //         }
//   //       }
//   //     }
//   //   }
//   //   if(3*jj!=NBVOLUMETET_VOLUME_HDIV(p)) SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong order %d != %d",jj,NBVOLUMETET_VOLUME_HDIV(p));
//   // }
//   PetscFunctionReturn(0);
// }
