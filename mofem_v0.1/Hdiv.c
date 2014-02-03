/**  
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

// based on Hierarchic Finite Element Bases on Unstructured Tetrahedral
// Meshes, by Mark Ainsworth and Joe Coyle

#include<FEM.h>

#include<H1HdivHcurlL2.h>
#include<strings.h>
#include<assert.h>

PetscErrorCode Hdiv_EdgeFaceShapeFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *PHI_f_e[4][3],int GDIM) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  const int face_edges_nodes[] = { 0,1, 1,2, 2,0 };
  const int face_oposite_edges_node[] = { 2, 0, 1 };
  double Phi_f_e[4][3][3];
  int ff = 0;
  for(;ff<4;ff++) {
    int ee = 0;
    for(;ee<3;ee++) {
      double _Spin_[9];
      int n0_idx = faces_nodes[3*ff+face_edges_nodes[2*ee+0]];
      int n1_idx = faces_nodes[3*ff+face_edges_nodes[2*ee+1]];
      ierr = Spin(_Spin_,&diffN[3*n0_idx]); CHKERRQ(ierr);
      cblas_dgemv(CblasRowMajor,CblasNoTrans,
	3,3,1.,_Spin_,3,&diffN[3*n1_idx],1,0,Phi_f_e[ff][ee],1);
    }
  }
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*4;
    ff = 0;
    for(;ff<4;ff++) {
      if(p[ff]<1) continue;
      int shift = ii*NBFACE_EDGE_Hdiv(p[ff]); 
      int ee = 0;
      for(;ee<3;ee++) {
	int n0_idx = faces_nodes[3*ff+face_edges_nodes[2*ee+0]];
	int n1_idx = faces_nodes[3*ff+face_edges_nodes[2*ee+1]];
	double ksi_0i = N[node_shift+n1_idx]-N[node_shift+n0_idx];
	double Psi_l[p[ff]+1];
	ierr = Lagrange_basis(p[ff],ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
	int nOposite_idx = faces_nodes[3*ff+face_oposite_edges_node[ee]];
	double lambda = N[node_shift+nOposite_idx];
	int l = 0;
	for(;l<=p[ff]-1;l++) {
	  int idx = 3*shift+3*l;
	  cblas_dcopy(3,Phi_f_e[ff][ee],1,&(PHI_f_e[ff][ee])[idx],1);
	  cblas_dscal(3,lambda*Psi_l[l],&(PHI_f_e[ff][ee])[idx],1);
	  //int dd = 0;
	  //for(;dd<3;dd++) (PHI_f_e[ff][ee])[idx+dd] = ((ii+1)*10000)+(ff*1000)+(ee*100)+(dd*10)+l;
	}
	if(l!=NBFACE_EDGE_Hdiv(p[ff])) SETERRQ2(PETSC_COMM_SELF,1,"wrong order %d != %d",l,NBFACE_EDGE_Hdiv(p[ff]));
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Hdiv_FaceBubbleShapeFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *PHI_f[],int GDIM) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  double Phi_f[4][3];
  int ff = 0;
  for(;ff<4;ff++) {
    int vert_i = faces_nodes[3*ff+1];
    int vert_j = faces_nodes[3*ff+2];
    double _Spin_[9];
    ierr = Spin(_Spin_,&diffN[3*vert_i]); CHKERRQ(ierr);
    cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1.,_Spin_,3,&diffN[3*vert_j],1,0,Phi_f[ff],1);
  }
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*4;
    ff = 0;
    for(;ff<4;ff++) {
      if(p[ff]<3) continue;
      double ksi_0i = N[ node_shift+faces_nodes[3*ff+1] ] - N[ node_shift+faces_nodes[3*ff+0] ];
      double ksi_0j = N[ node_shift+faces_nodes[3*ff+2] ] - N[ node_shift+faces_nodes[3*ff+0] ];
      double Psi_l[p[ff]+1],Psi_m[p[ff]+1];
      ierr = Lagrange_basis(p[ff],ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
      ierr = Lagrange_basis(p[ff],ksi_0j,NULL,Psi_m,NULL,3); CHKERRQ(ierr);
      double Beta_0ij = 
	N[node_shift+faces_nodes[3*ff+0]]*N[node_shift+faces_nodes[3*ff+1]]*N[node_shift+faces_nodes[3*ff+2]];
      int shift = ii*NBFACE_FACE_Hdiv(p[ff]); 
      int jj = 0;
      int oo = 0;
      for(;oo<=p[ff]-3;oo++) {
	int l = 0;
	for(;l<=oo;l++) {
	  int m = 0;
	  m = oo - l;
	  if(m>=0) {
	    double *phi_f = &(PHI_f[ff])[3*shift+3*jj];
	    cblas_dcopy(3,Phi_f[ff],1,phi_f,1);
	    cblas_dscal(3,Beta_0ij*Psi_l[l]*Psi_m[m],phi_f,1);
	    jj++;
	  } 
	}
      }
      if(jj!=NBFACE_FACE_Hdiv(p[ff])) SETERRQ2(PETSC_COMM_SELF,1,"wrong order %d != %d",jj,NBFACE_FACE_Hdiv(p[ff]));
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Hdiv_EdgeBasedVolumeShapeFunctions_MBTET(
  int p,double *coords,double *N,double *PHI_v_e[6],int GDIM) {
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
    int shift = ii*NBVOLUME_EDGE_Hdiv(p);
    ee = 0;
    for(;ee<6;ee++) {
      double Beta_e = N[ node_shift+edges_nodes[2*ee+1] ]*N[ node_shift+edges_nodes[2*ee+0] ];
      double ksi_0i = N[ node_shift+edges_nodes[2*ee+1] ]-N[ node_shift+edges_nodes[2*ee+0] ];
      double Psi_l[p+1];
      ierr = Lagrange_basis(p,ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
      int l = 0;
      for(;l<=p-2;l++) {
	cblas_dcopy(3,tau_e[ee],1,&(PHI_v_e[ee])[3*shift+3*l],1);
	cblas_dscal(3,Beta_e*Psi_l[l],&(PHI_v_e[ee])[3*shift+3*l],1);
      }
      if(l!=NBVOLUME_EDGE_Hdiv(p)) SETERRQ2(PETSC_COMM_SELF,1,"wrong order %d != %d",l,NBVOLUME_FACE_Hdiv(p));
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Hdiv_FaceBasedVolumeShapeFunctions_MBTET(
  int p,double *coords,double *N,double *PHI_v_f[],int GDIM) {
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
      double Psi_l[ p+1 ],Psi_m[ p+1 ];
      ierr = Lagrange_basis(p,ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
      ierr = Lagrange_basis(p,ksi_0j,NULL,Psi_m,NULL,3); CHKERRQ(ierr);
      double Beta_0ij = 
	N[node_shift+faces_nodes[3*ff+0]]*N[node_shift+faces_nodes[3*ff+1]]*N[node_shift+faces_nodes[3*ff+2]];
      int shift = ii*NBVOLUME_FACE_Hdiv(p);
      int jj = 0;
      int oo = 0;
      for(;oo<=p-3;oo++) {
	int l = 0;
	for(;l<=oo;l++) {
	  int m = oo - l;
	  if(m>=0) {
	    cblas_dcopy(3,tau_0i[ff],1,&(PHI_v_f[ff])[3*shift + 3*jj],1);
	    cblas_dscal(3,Beta_0ij*Psi_l[l]*Psi_m[m],&(PHI_v_f[ff])[3*shift + 3*jj],1);
	    jj++;
	    cblas_dcopy(3,tau_0j[ff],1,&(PHI_v_f[ff])[3*shift + 3*jj],1);
	    cblas_dscal(3,Beta_0ij*Psi_l[l]*Psi_m[m],&(PHI_v_f[ff])[3*shift + 3*jj],1);
	    jj++;
	  }
	}
      }
      if(jj!=NBVOLUME_FACE_Hdiv(p)) SETERRQ2(PETSC_COMM_SELF,1,"wrong order %d != %d",jj,NBVOLUME_FACE_Hdiv(p));
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Hdiv_VolumeBubbleShapeFunctions_MBTET(
  int p,double *coords,double *N,double *PHI_v,int GDIM) {
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
    double Psi_l[p+1],Psi_m[p+1],Psi_n[p+1];
    ierr = Lagrange_basis(p,ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
    ierr = Lagrange_basis(p,ksi_0j,NULL,Psi_m,NULL,3); CHKERRQ(ierr);
    ierr = Lagrange_basis(p,ksi_0k,NULL,Psi_n,NULL,3); CHKERRQ(ierr);
    int shift = ii*NBVOLUME_VOLUME_Hdiv(p);
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
	    PHI_v[3*shift + 3*3*jj + 3*0 + 0] = s*ed[0][0]; 
	    PHI_v[3*shift + 3*3*jj + 3*0 + 1] = s*ed[0][1]; 
	    PHI_v[3*shift + 3*3*jj + 3*0 + 2] = s*ed[0][2];
	    //
	    PHI_v[3*shift + 3*3*jj + 3*1 + 0] = s*ed[1][0]; 
	    PHI_v[3*shift + 3*3*jj + 3*1 + 1] = s*ed[1][1]; 
	    PHI_v[3*shift + 3*3*jj + 3*1 + 2] = s*ed[1][2];
	    //
	    PHI_v[3*shift + 3*3*jj + 3*2 + 0] = s*ed[2][0]; 
	    PHI_v[3*shift + 3*3*jj + 3*2 + 1] = s*ed[2][1]; 
	    PHI_v[3*shift + 3*3*jj + 3*2 + 2] = s*ed[2][2];
	    jj++;
	  }
	}
      }
    }
    if(3*jj!=NBVOLUME_VOLUME_Hdiv(p)) SETERRQ2(PETSC_COMM_SELF,1,"wrong order %d != %d",jj,NBVOLUME_VOLUME_Hdiv(p));
  }
  PetscFunctionReturn(0);
}
