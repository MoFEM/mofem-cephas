/** \file hdiv.c

  Based on Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by Mark Ainsworth and Joe Coyle
  Shape functions for MBTRI/MBTET and HDiv space

*/

/**
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#include <petscsys.h>
#include <cblas.h>

#include <definitions.h>
#include <fem_tools.h>
#include <h1_hdiv_hcurl_l2.h>

PetscErrorCode Spin(double *spinOmega,double *vecOmega);

PetscErrorCode Hdiv_EdgeFaceShapeFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *PHI_f_e[4][3],double *diffPHI_f_e[4][3],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  int ff = 0;
  for(;ff<4;ff++) {
    if(diffPHI_f_e!=NULL) {
      ierr = Hdiv_EdgeFaceShapeFunctions_MBTET_ON_FACE(&faces_nodes[3*ff],p[ff],N,diffN,PHI_f_e[ff],diffPHI_f_e[ff],GDIM,4,base_polynomials); CHKERRQ(ierr);
    } else {
      ierr = Hdiv_EdgeFaceShapeFunctions_MBTET_ON_FACE(&faces_nodes[3*ff],p[ff],N,diffN,PHI_f_e[ff],NULL,GDIM,4,base_polynomials); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Hdiv_FaceBubbleShapeFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *PHI_f[],double *diffPHI_f[],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  int ff = 0;
  for(;ff<4;ff++) {
    double *diff = NULL;
    if(diffPHI_f!=NULL) diff = diffPHI_f[ff];
    ierr = Hdiv_FaceBubbleShapeFunctions_MBTET_ON_FACE(&faces_nodes[3*ff],p[ff],N,diffN,PHI_f[ff],diff,GDIM,4,base_polynomials); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Hdiv_EdgeBasedVolumeShapeFunctions_MBTET(
  int p,double *coords,double *N,double *diffN,double *PHI_v_e[6],double *diffPHI_v_e[6],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;
  if(p<2) PetscFunctionReturn(0);
  PetscErrorCode ierr;
  const int edges_nodes[] = { 0,1, 1,2, 2,0, 0,3, 1,3, 2,3 };
  double tau_e[6][3];
  int ee = 0;
  for(;ee<6;ee++) {
    cblas_dcopy(3,&coords[3*edges_nodes[2*ee+1]],1,tau_e[ee],1);
    cblas_daxpy(3,-1,&coords[3*edges_nodes[2*ee+0]],1,tau_e[ee],1);
    double nrm2 = cblas_dnrm2(3,tau_e[ee],1);
    cblas_dscal(3,1./nrm2,tau_e[ee],1);
  }
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*4;
    int shift = ii*NBVOLUMETET_EDGE_HDIV(p);
    ee = 0;
    for(;ee<6;ee++) {
      double Beta_e = N[ node_shift+edges_nodes[2*ee+1] ]*N[ node_shift+edges_nodes[2*ee+0] ];
      double ksi_0i = N[ node_shift+edges_nodes[2*ee+1] ]-N[ node_shift+edges_nodes[2*ee+0] ];
      double diff_Beta_e[3];
      double diff_ksi_0i[3];
      double Psi_l[p+1],diff_Psi_l[3*(p+1)];
      if(diffPHI_v_e!=NULL) {
        if(diffN == NULL) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        int dd = 0;
        for(;dd<3;dd++) {
          diff_Beta_e[dd] =
          diffN[3*edges_nodes[2*ee+1]+dd]*N[node_shift+edges_nodes[2*ee+0]]
          + N[node_shift+edges_nodes[2*ee+1]]*diffN[3*edges_nodes[2*ee+0]+dd];
          diff_ksi_0i[dd] =
          diffN[3*edges_nodes[2*ee+1]+dd]-diffN[3*edges_nodes[2*ee+0]+dd];
        }
        ierr = base_polynomials(p,ksi_0i,diff_ksi_0i,Psi_l,diff_Psi_l,3); CHKERRQ(ierr);
      } else {
        ierr = base_polynomials(p,ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
      }
      int l = 0;
      for(;l<=p-2;l++) {
        cblas_dcopy(3,tau_e[ee],1,&(PHI_v_e[ee])[3*shift+3*l],1);
        cblas_dscal(3,Beta_e*Psi_l[l],&(PHI_v_e[ee])[3*shift+3*l],1);
        if(diffPHI_v_e!=NULL) {
          int dd = 0;
          for(;dd<3;dd++) {
            cblas_dcopy(3,tau_e[ee],1,&(diffPHI_v_e[ee])[9*shift+9*l+3*dd],1);
            double diff = diff_Beta_e[dd]*Psi_l[l] + Beta_e*diff_Psi_l[dd*(p+1)+l];
            cblas_dscal(3,diff,&(diffPHI_v_e[ee])[9*shift+9*l+3*dd],1);
          }
        }
      }
      if(l!=NBVOLUMETET_EDGE_HDIV(p)) SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong order %d != %d",l,NBVOLUMETET_FACE_HDIV(p));
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Hdiv_FaceBasedVolumeShapeFunctions_MBTET(
  int p,double *coords,double *N,double *diffN,double *PHI_v_f[],double *diffPHI_v_f[],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;
  if(p<3) PetscFunctionReturn(0);
  PetscErrorCode ierr;
  const int faces_nodes[] = { 0,1,3, 1,2,3, 0,2,3, 0,1,2 };
  double tau_0i[4][3],tau_0j[4][3];
  int ff = 0;
  for(;ff<4;ff++) {
    int idx_node0 = faces_nodes[3*ff+0];
    int idx_node1 = faces_nodes[3*ff+1];
    int idx_node2 = faces_nodes[3*ff+2];
    cblas_dcopy(3,&coords[3*idx_node1],1,tau_0i[ff],1);
    cblas_daxpy(3,-1,&coords[3*idx_node0],1,tau_0i[ff],1);
    double nrm2_0i = cblas_dnrm2(3,tau_0i[ff],1);
    cblas_dscal(3,1./nrm2_0i,tau_0i[ff],1);
    cblas_dcopy(3,&coords[3*idx_node2],1,tau_0j[ff],1);
    cblas_daxpy(3,-1,&coords[3*idx_node0],1,tau_0j[ff],1);
    double nrm2_0j = cblas_dnrm2(3,tau_0j[ff],1);
    cblas_dscal(3,1./nrm2_0j,tau_0j[ff],1);
  }
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*4;
    ff = 0;
    for(;ff<4;ff++) {
      double ksi_0i = N[ node_shift+faces_nodes[3*ff+1] ] - N[ node_shift+faces_nodes[3*ff+0] ];
      double ksi_0j = N[ node_shift+faces_nodes[3*ff+2] ] - N[ node_shift+faces_nodes[3*ff+0] ];
      double diff_ksi_0i[3],diff_ksi_0j[3];
      if(diffPHI_v_f!=NULL) {
        if(diffN == NULL) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        int dd = 0;
        for(;dd<3;dd++) {
          diff_ksi_0i[dd] = diffN[ 3*faces_nodes[3*ff+1] + dd ] - diffN[ 3*faces_nodes[3*ff+0] + dd ];
          diff_ksi_0j[dd] = diffN[ 3*faces_nodes[3*ff+2] + dd ] - diffN[ 3*faces_nodes[3*ff+0] + dd ];
        }
      }
      double Psi_l[ p+1 ],Psi_m[ p+1 ];
      double diff_Psi_l[ 3*(p+1) ],diff_Psi_m[ 3*(p+1) ];
      if(diffPHI_v_f != NULL) {
        ierr = base_polynomials(p,ksi_0i,diff_ksi_0i,Psi_l,diff_Psi_l,3); CHKERRQ(ierr);
        ierr = base_polynomials(p,ksi_0j,diff_ksi_0j,Psi_m,diff_Psi_m,3); CHKERRQ(ierr);
      } else {
        ierr = base_polynomials(p,ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
        ierr = base_polynomials(p,ksi_0j,NULL,Psi_m,NULL,3); CHKERRQ(ierr);
      }
      double Beta_0ij =
      N[node_shift+faces_nodes[3*ff+0]]*N[node_shift+faces_nodes[3*ff+1]]*N[node_shift+faces_nodes[3*ff+2]];
      double diff_Beta_0ij[3];
      if(diffPHI_v_f!=NULL) {
        int dd = 0;
        for(;dd<3;dd++) {
          diff_Beta_0ij[dd] =
          diffN[3*faces_nodes[3*ff+0]+dd]*N[node_shift+faces_nodes[3*ff+1]]*N[node_shift+faces_nodes[3*ff+2]]+
          N[node_shift+faces_nodes[3*ff+0]]*diffN[3*faces_nodes[3*ff+1]+dd]*N[node_shift+faces_nodes[3*ff+2]]+
          N[node_shift+faces_nodes[3*ff+0]]*N[node_shift+faces_nodes[3*ff+1]]*diffN[3*faces_nodes[3*ff+2]+dd];
        }
      }
      int shift = ii*NBVOLUMETET_FACE_HDIV(p);
      int jj = 0;
      int oo = 0;
      for(;oo<=p-3;oo++) {
        int l = 0;
        for(;l<=oo;l++) {
          int m = oo - l;
          if(m>=0) {
            double scale = Beta_0ij*Psi_l[l]*Psi_m[m];
            cblas_dcopy(3,tau_0i[ff],1,&(PHI_v_f[ff])[3*shift + 3*jj],1);
            cblas_dscal(3,scale,&(PHI_v_f[ff])[3*shift + 3*jj],1);
            double diff[3];
            if(diffPHI_v_f!=NULL) {
              int dd = 0;
              for(;dd<3;dd++) {
                cblas_dcopy(3,tau_0i[ff],1,&(diffPHI_v_f[ff])[9*shift + 9*jj + 3*dd],1);
                diff[dd] =
                diff_Beta_0ij[dd]*Psi_l[l]*Psi_m[m]+
                Beta_0ij*diff_Psi_l[dd*(p+1)+l]*Psi_m[m]+
                Beta_0ij*Psi_l[l]*diff_Psi_m[dd*(p+1)+m];
                cblas_dscal(3,diff[dd],&(diffPHI_v_f[ff])[9*shift + 9*jj + 3*dd],1);
              }
            }
            jj++;
            cblas_dcopy(3,tau_0j[ff],1,&(PHI_v_f[ff])[3*shift + 3*jj],1);
            cblas_dscal(3,scale,&(PHI_v_f[ff])[3*shift + 3*jj],1);
            if(diffPHI_v_f!=NULL) {
              int dd = 0;
              for(;dd<3;dd++) {
                cblas_dcopy(3,tau_0j[ff],1,&(diffPHI_v_f[ff])[9*shift + 9*jj + 3*dd],1);
                cblas_dscal(3,diff[dd],&(diffPHI_v_f[ff])[9*shift + 9*jj + 3*dd],1);
              }
            }
            jj++;
          }
        }
      }
      if(jj!=NBVOLUMETET_FACE_HDIV(p)) SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong order %d != %d",jj,NBVOLUMETET_FACE_HDIV(p));
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Hdiv_VolumeBubbleShapeFunctions_MBTET(
  int p,double *coords,double *N,double *diffN,double *PHI_v,double *diffPHI_v,int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;
  if(p<4) PetscFunctionReturn(0);
  PetscErrorCode ierr;
  double ed[3][3];
  int nn = 0;
  for(;nn<3;nn++) {
    cblas_dcopy(3,&coords[3*(nn+1)],1,ed[nn],1);
    cblas_daxpy(3,-1,&coords[0],1,ed[nn],1);
    double nrm2 = cblas_dnrm2(3,ed[nn],1);
    cblas_dscal(3,1./nrm2,ed[nn],1);
  }
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*4;
    double Beta_0ijk =
    N[ node_shift + 0]*N[ node_shift + 1]*N[ node_shift + 2]*N[ node_shift + 3];
    double ksi_0i = N[ node_shift+1 ] - N[ node_shift+0 ];
    double ksi_0j = N[ node_shift+2 ] - N[ node_shift+0 ];
    double ksi_0k = N[ node_shift+3 ] - N[ node_shift+0 ];
    double diff_Beta_0ijk[3] = {0,0,0};
    double diff_ksi_0i[3],diff_ksi_0j[3],diff_ksi_0k[3];
    double Psi_l[p+1],Psi_m[p+1],Psi_n[p+1];
    double diff_Psi_l[3*(p+1)],diff_Psi_m[3*(p+1)],diff_Psi_n[3*(p+1)];
    if(diffPHI_v != NULL) {
      int dd = 0;
      for(;dd<3;dd++) {
        diff_Beta_0ijk[dd] =
        diffN[ 3*0+dd ]*N[ node_shift + 1]*N[ node_shift + 2]*N[ node_shift + 3]+
        N[ node_shift + 0]*diffN[ 3*1+dd ]*N[ node_shift + 2]*N[ node_shift + 3]+
        N[ node_shift + 0]*N[ node_shift + 1]*diffN[ 3*2+dd ]*N[ node_shift + 3]+
        N[ node_shift + 0]*N[ node_shift + 1]*N[ node_shift + 2]*diffN[ 3*3+dd ];
        diff_ksi_0i[dd] = diffN[ 3*1+dd ] - diffN[ 3*0+dd ];
        diff_ksi_0j[dd] = diffN[ 3*2+dd ] - diffN[ 3*0+dd ];
        diff_ksi_0k[dd] = diffN[ 3*3+dd ] - diffN[ 3*0+dd ];
      }
      ierr = base_polynomials(p,ksi_0i,diff_ksi_0i,Psi_l,diff_Psi_l,3); CHKERRQ(ierr);
      ierr = base_polynomials(p,ksi_0j,diff_ksi_0j,Psi_m,diff_Psi_m,3); CHKERRQ(ierr);
      ierr = base_polynomials(p,ksi_0k,diff_ksi_0k,Psi_n,diff_Psi_n,3); CHKERRQ(ierr);
    } else {
      ierr = base_polynomials(p,ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
      ierr = base_polynomials(p,ksi_0j,NULL,Psi_m,NULL,3); CHKERRQ(ierr);
      ierr = base_polynomials(p,ksi_0k,NULL,Psi_n,NULL,3); CHKERRQ(ierr);
    }
    int shift = ii*NBVOLUMETET_VOLUME_HDIV(p);
    int jj = 0;
    int oo = 0;
    for(;oo<=p-4;oo++) {
      int l = 0;
      for(;l<=oo;l++) {
        int m = 0;
        for(;(l+m)<=oo;m++) {
          int n = oo - l - m;
          if(n>=0) {
            double s = Beta_0ijk*Psi_l[l]*Psi_m[m]*Psi_n[n];
            int kk = 0;
            for(;kk<3;kk++) {
              PHI_v[3*shift + 3*3*jj + 3*0 + kk] = s*ed[0][kk];
              PHI_v[3*shift + 3*3*jj + 3*1 + kk] = s*ed[1][kk];
              PHI_v[3*shift + 3*3*jj + 3*2 + kk] = s*ed[2][kk];
            }
            if(diffPHI_v!=NULL) {
              int dd = 0;
              for(;dd<3;dd++) {
                double diff =
                diff_Beta_0ijk[dd]*Psi_l[l]*Psi_m[m]*Psi_n[n]+
                Beta_0ijk*diff_Psi_l[dd*(p+1)+l]*Psi_m[m]*Psi_n[n]+
                Beta_0ijk*Psi_l[l]*diff_Psi_m[dd*(p+1)+m]*Psi_n[n]+
                Beta_0ijk*Psi_l[l]*Psi_m[m]*diff_Psi_n[dd*(p+1)+n];
                int kk = 0;
                for(;kk<3;kk++) {
                  diffPHI_v[9*shift + 3*9*jj + 9*0 + 3*dd + kk] = diff*ed[0][kk];
                  diffPHI_v[9*shift + 3*9*jj + 9*1 + 3*dd + kk] = diff*ed[1][kk];
                  diffPHI_v[9*shift + 3*9*jj + 9*2 + 3*dd + kk] = diff*ed[2][kk];
                }
              }
            }
            jj++;
          }
        }
      }
    }
    if(3*jj!=NBVOLUMETET_VOLUME_HDIV(p)) SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong order %d != %d",jj,NBVOLUMETET_VOLUME_HDIV(p));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Hdiv_EdgeFaceShapeFunctions_MBTET_ON_FACE(
  int *faces_nodes,int p,double *N,double *diffN,double *PHI_f_e[3],double *diffPHI_f_e[3],int GDIM,int NB,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;
  if(p<1) PetscFunctionReturn(0);
  PetscErrorCode ierr;
  const int face_edges_nodes[] = { 0,1, 1,2, 2,0 };
  const int face_oposite_edges_node[] = { 2, 0, 1 };
  double Phi_f_e[3][3];
  if(diffN!=NULL) {
    int ee = 0;
    for(;ee<3;ee++) {
      double _Spin_[9];
      int n0_idx = faces_nodes[face_edges_nodes[2*ee+0]];
      int n1_idx = faces_nodes[face_edges_nodes[2*ee+1]];
      ierr = Spin(_Spin_,&diffN[3*n0_idx]); CHKERRQ(ierr);
      cblas_dgemv(
        CblasRowMajor,CblasNoTrans,
        3,3,1.,_Spin_,3,&diffN[3*n1_idx],1,0,Phi_f_e[ee],1
      );
    }
  } else {
    int ee = 0;
    for(;ee<3;ee++) {
      Phi_f_e[ee][0] = 1;
      Phi_f_e[ee][1] = 0;
      Phi_f_e[ee][2] = 0;
    }
  }
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*NB;
    int shift = ii*NBFACETRI_EDGE_HDIV(p);
    int ee = 0;
    for(;ee<3;ee++) {
      int n0_idx = faces_nodes[face_edges_nodes[2*ee+0]];
      int n1_idx = faces_nodes[face_edges_nodes[2*ee+1]];
      double ksi_0i = N[node_shift+n1_idx]-N[node_shift+n0_idx];
      double Psi_l[p+1],diff_Psi_l[3*(p+1)];
      double diff_ksi_0i[3];
      if(diffPHI_f_e != 0) {
        if(diffN == NULL) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        int dd = 0;
        for(;dd<3;dd++) {
          diff_ksi_0i[dd] = diffN[3*n1_idx+dd] - diffN[3*n0_idx+dd];
        }
        ierr = base_polynomials(p,ksi_0i,diff_ksi_0i,Psi_l,diff_Psi_l,3); CHKERRQ(ierr);
      } else {
        ierr = base_polynomials(p,ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
      }
      int nOposite_idx = faces_nodes[face_oposite_edges_node[ee]];
      double lambda = N[node_shift+nOposite_idx];
      double diff_lambda[3];
      if(diffPHI_f_e!=NULL) {
        int dd = 0;
        for(;dd<3;dd++) {
          diff_lambda[dd] = diffN[3*nOposite_idx+dd];
        }
      }
      int l = 0;
      for(;l<=p-1;l++) {
        int idx = 3*shift+3*l;
        cblas_dcopy(3,Phi_f_e[ee],1,&(PHI_f_e[ee])[idx],1);
        cblas_dscal(3,lambda*Psi_l[l],&(PHI_f_e[ee])[idx],1);
        if(diffPHI_f_e!=NULL) {
          int diff_idx = 9*shift+9*l;
          int dd = 0;
          for(;dd<3;dd++) {
            cblas_dcopy(3,Phi_f_e[ee],1,&(diffPHI_f_e[ee])[diff_idx+3*dd],1);
            double diff =
            diff_lambda[dd]*Psi_l[l] + lambda*diff_Psi_l[dd*(p+1)+l];
            cblas_dscal(3,diff,&(diffPHI_f_e[ee])[diff_idx+3*dd],1);
          }
        }
      }
      if(l!=NBFACETRI_EDGE_HDIV(p)) SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong order %d != %d",l,NBFACETRI_EDGE_HDIV(p));
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Hdiv_FaceBubbleShapeFunctions_MBTET_ON_FACE(
  int *faces_nodes,int p,double *N,double *diffN,double *PHI_f,double *diffPHI_f,int GDIM,int NB,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;
  if(p<3) PetscFunctionReturn(0);
  PetscErrorCode ierr;
  double Phi_f[3];
  if(diffN!=NULL) {
    int vert_i = faces_nodes[1];
    int vert_j = faces_nodes[2];
    double spin[9];
    ierr = Spin(spin,&diffN[3*vert_i]); CHKERRQ(ierr);
    cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1.,spin,3,&diffN[3*vert_j],1,0,Phi_f,1);
  } else {
    Phi_f[0] = 1;
    Phi_f[1] = 0;
    Phi_f[2] = 0;
  }
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*NB;
    double ksi_0i = N[ node_shift+faces_nodes[1] ] - N[ node_shift+faces_nodes[0] ];
    double ksi_0j = N[ node_shift+faces_nodes[2] ] - N[ node_shift+faces_nodes[0] ];
    double Psi_l[p+1],Psi_m[p+1];
    double diff_Psi_l[3*(p+1)],diff_Psi_m[3*(p+1)];
    if(diffPHI_f!=NULL) {
      if(diffN == NULL) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      double diff_ksi_0i[3];
      double diff_ksi_0j[3];
      int dd = 0;
      for(;dd<3;dd++) {
        diff_ksi_0i[dd] = diffN[3*faces_nodes[1] + dd] - diffN[3*faces_nodes[0] + dd];
        diff_ksi_0j[dd] = diffN[3*faces_nodes[2] + dd] - diffN[3*faces_nodes[0] + dd];
      }
      ierr = base_polynomials(p,ksi_0i,diff_ksi_0i,Psi_l,diff_Psi_l,3); CHKERRQ(ierr);
      ierr = base_polynomials(p,ksi_0j,diff_ksi_0j,Psi_m,diff_Psi_m,3); CHKERRQ(ierr);
    } else {
      ierr = base_polynomials(p,ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
      ierr = base_polynomials(p,ksi_0j,NULL,Psi_m,NULL,3); CHKERRQ(ierr);
    }
    double Beta_0ij =
    N[node_shift+faces_nodes[0]]*N[node_shift+faces_nodes[1]]*N[node_shift+faces_nodes[2]];
    double diff_Beta_0ij[3];
    if(diffPHI_f!=NULL) {
      int dd = 0;
      for(;dd<3;dd++) {
        diff_Beta_0ij[dd] =
        diffN[3*faces_nodes[0]+dd]*N[node_shift+faces_nodes[1]]*N[node_shift+faces_nodes[2]] +
        N[node_shift+faces_nodes[0]]*diffN[3*faces_nodes[1]+dd]*N[node_shift+faces_nodes[2]] +
        N[node_shift+faces_nodes[0]]*N[node_shift+faces_nodes[1]]*diffN[3*faces_nodes[2]+dd];
      }
    }
    int shift = ii*NBFACETRI_FACE_HDIV(p);
    int jj = 0;
    int oo = 0;
    for(;oo<=p-3;oo++) {
      int l = 0;
      for(;l<=oo;l++) {
        int m = 0;
        m = oo - l;
        if(m>=0) {
          double *phi_f = &(PHI_f)[3*shift+3*jj];
          cblas_dcopy(3,Phi_f,1,phi_f,1);
          cblas_dscal(3,Beta_0ij*Psi_l[l]*Psi_m[m],phi_f,1);
          if(diffPHI_f!=NULL) {
            int dd = 0;
            for(;dd<3;dd++) {
              double *diff_phi_f = &(diffPHI_f)[9*shift+9*jj+3*dd];
              cblas_dcopy(3,Phi_f,1,diff_phi_f,1);
              double diff =
              diff_Beta_0ij[dd]*Psi_l[l]*Psi_m[m]+
              Beta_0ij*diff_Psi_l[dd*(p+1)+l]*Psi_m[m]+
              Beta_0ij*Psi_l[l]*diff_Psi_m[dd*(p+1)+m];
              cblas_dscal(3,diff,diff_phi_f,1);
            }
          }
          jj++;
        }
      }
    }
    if(jj!=NBFACETRI_FACE_HDIV(p)) SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong order %d != %d",jj,NBFACETRI_FACE_HDIV(p));
  }
  PetscFunctionReturn(0);
}
