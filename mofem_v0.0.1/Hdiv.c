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
// Shape functions for MBTRI and H1 approximation

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
  double Phi_f_e[4*3*3];
  int ff = 0;
  for(;ff<4;ff++) {
    int ee = 0;
    for(;ee<3;ee++) {
      double _Spin_[9];
      ierr = Spin(_Spin_,&diffN[3*faces_nodes[2*ff+face_edges_nodes[2*ee+0]]]); CHKERRQ(ierr);
      cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1.,_Spin_,3,&diffN[3*faces_nodes[2*ff+face_edges_nodes[2*ee+0]]],1,
	0,&(Phi_f_e[ff*3*3+ee*3]),1);
    }
  }
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*4;
    ff = 0;
    for(;ff<4;ff++) {
      int ee = 0;
      for(;ee<3;ee++) {
	double ksi_0i = 
	  (N[node_shift+faces_nodes[3*ff+face_edges_nodes[2*ee+1]]] 
	  - N[node_shift+faces_nodes[3*ff+face_edges_nodes[2*ee+0]]]);
	double Psi_l[p[ff]];
	ierr = Lagrange_basis(p[ff],ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
	double lambda;
	lambda = N[node_shift+faces_nodes[3*ff+face_oposite_edges_node[ee]]];
	int jj = 0;
	cblas_dcopy(3,&Phi_f_e[ff*3*3+ee*3],1,&(PHI_f_e[ff][ee])[3*jj],1);
	cblas_dscal(3,lambda,&(PHI_f_e[ff][ee])[3*jj],1);
	jj++;
	int l = 0;
	for(;l<p[ff];l++) {
	  cblas_dcopy(3,&Phi_f_e[ff*3*3+ee*3],1,&(PHI_f_e[ff][ee])[3*jj],1);
	  cblas_dscal(3,lambda*Psi_l[l],&(PHI_f_e[ff][ee])[3*jj],1);
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Hdiv_FaceBubbleShapeFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *PHI_f[],int GDIM) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  double Phi_f[4*3];
  int ff = 0;
  for(;ff<4;ff++) {
    int vert_i = faces_nodes[3*ff+1];
    int vert_j = faces_nodes[3*ff+2];
    double _Spin_[9];
    ierr = Spin(_Spin_,&diffN[3*vert_i]); CHKERRQ(ierr);
    cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1.,_Spin_,3,&diffN[3*vert_j],1,0,&(Phi_f[3*ff]),1);
  }
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*4;
    ff = 0;
    for(;ff<4;ff++) {
      double ksi_0i = N[ node_shift+faces_nodes[3*ff+1] ] - N[ node_shift+faces_nodes[3*ff+0] ];
      double ksi_0j = N[ node_shift+faces_nodes[3*ff+2] ] - N[ node_shift+faces_nodes[3*ff+0] ];
      double Psi_l[p[ff]-3],Psi_m[p[ff]-3];
      ierr = Lagrange_basis(p[ff]-2,ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
      ierr = Lagrange_basis(p[ff]-2,ksi_0j,NULL,Psi_m,NULL,3); CHKERRQ(ierr);
      double Beta_0ij = 
	N[node_shift+faces_nodes[3*ff+0]]*N[node_shift+faces_nodes[3*ff+1]]*N[node_shift+faces_nodes[3*ff+2]];
      int jj = 0;
      int l = 0;
      for(;l<p[ff]-3;l++) {
	int m = 0;
	for(;(l+m)<p[ff]-3;m++) {
	  cblas_dcopy(3,Phi_f,1,&(PHI_f[ff])[3*jj],1);
	  cblas_dscal(3,Beta_0ij*Psi_l[l]*Psi_m[m],&(PHI_f[ff])[3*jj],1);
	  jj++;
	} 
      }   
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Hdiv_EdgeBasedVolumeShapeFunctions_MBTET(
  int *sense,int *p,double *coords,double *N,double *PHI_v_e[6],int GDIM) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  const int edges_nodes[2*6] = { 0,1, 1,2, 2,0, 0,3, 1,3, 2,3 };
  double tau_e[6*3];
  int ee = 0;
  for(;ee<6;ee++) {
    cblas_dcopy(3,&coords[edges_nodes[2*ee+1]],1,&tau_e[3*ee],1);
    cblas_daxpy(3,-1,&coords[edges_nodes[2*ee+0]],1,&tau_e[3*ee],1);
    cblas_dscal(3,sense[ee],&tau_e[3*ee],1);
  } 
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*4;
    ee = 0;
    for(;ee<6;ee++) {
      double Beta_e = N[ node_shift+edges_nodes[2*ee+0] ]*N[ node_shift+edges_nodes[2*ee+1] ];
      double ksi_0i = 
	(N[ node_shift+edges_nodes[3*ee+1] ] - N[ node_shift+edges_nodes[3*ee+0] ])*sense[ee];
      double Psi_l[ p[ee]-2 ];
      ierr = Lagrange_basis(p[ee]-2,ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
      int l = 0;
      for(;l<p[ee]-2;l++) {
	cblas_dcopy(3,&tau_e[ee*3],1,&(PHI_v_e[ee])[3*l],1);
	cblas_dscal(3,Beta_e*Psi_l[l],&(PHI_v_e[ee])[3*l],1);
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Hdiv_FaceBasedVolumeShapeFunctions_MBTET(
  int *faces_nodes,int *p,double *coords,double *N,double *PHI_v_f[],int GDIM) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  double tau_0i[4*3],tau_0j[4*3];
  int ff = 0;
  for(;ff<4;ff++) {
    cblas_dcopy(3,&coords[faces_nodes[3*ff+1]],1,&tau_0i[3*ff],1);
    cblas_daxpy(3,-1,&coords[faces_nodes[3*ff+0]],1,&tau_0i[3*ff],1);
    cblas_dcopy(3,&coords[faces_nodes[3*ff+2]],1,&tau_0j[3*ff],1);
    cblas_daxpy(3,-1,&coords[faces_nodes[3*ff+0]],1,&tau_0j[3*ff],1);
  }
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*4;
    ff = 0;
    for(;ff<4;ff++) {
      double ksi_0i = N[ node_shift+faces_nodes[3*ff+1] ] - N[ node_shift+faces_nodes[3*ff+0] ];
      double ksi_0j = N[ node_shift+faces_nodes[3*ff+2] ] - N[ node_shift+faces_nodes[3*ff+0] ];
      double Psi_l[ p[ff]-3 ],Psi_m[ p[ff]-3 ];
      ierr = Lagrange_basis(p[ff]-2,ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
      ierr = Lagrange_basis(p[ff]-2,ksi_0j,NULL,Psi_m,NULL,3); CHKERRQ(ierr);
      double Beta_0ij = 
	N[node_shift+faces_nodes[3*ff+0]]*N[node_shift+faces_nodes[3*ff+1]]*N[node_shift+faces_nodes[3*ff+2]];
      int jj = 0;
      int l = 0;
      for(;l<p[ff]-3;l++) {
	int m = 0;
	for(;(l+m)<p[ff]-3;m++) {
	  cblas_dcopy(3,&tau_0i[3*ff],1,&(PHI_v_f[ff])[3*2*jj + 0],1);
	  cblas_dscal(3,Beta_0ij*Psi_l[l]*Psi_m[m],&(PHI_v_f[ff])[3*2*jj + 0],1);
	  cblas_dcopy(3,&tau_0j[3*ff],1,&(PHI_v_f[ff])[3*2*jj + 3],1);
	  cblas_dscal(3,Beta_0ij*Psi_l[l]*Psi_m[m],&(PHI_v_f[ff])[3*2*jj + 3],1);
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Hdiv_VolumeBubbleShapeFunctions_MBTET(
  int *sense,int p,double *N,double *PHI_v,int GDIM) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*4;
    double Beta_0ijk = 
      N[ node_shift + 0]*N[ node_shift + 1]*N[ node_shift + 2]*N[ node_shift + 3];
    double ksi_0i = N[ node_shift+1 ] - N[ node_shift+0 ];
    double ksi_0j = N[ node_shift+2 ] - N[ node_shift+0 ];
    double ksi_0k = N[ node_shift+3 ] - N[ node_shift+0 ];
    double Psi_l[p-4],Psi_m[p-4],Psi_n[p-4];
    ierr = Lagrange_basis(p-4,ksi_0i,NULL,Psi_l,NULL,3); CHKERRQ(ierr);
    ierr = Lagrange_basis(p-4,ksi_0j,NULL,Psi_m,NULL,3); CHKERRQ(ierr);
    ierr = Lagrange_basis(p-4,ksi_0k,NULL,Psi_n,NULL,3); CHKERRQ(ierr);
    int jj = 0;
    int l = 0;
    for(;l<p-4;l++) {
      int m = 0;
      for(;(l+m)<p-4;m++) {
	int  n = 0;
	for(;(l+m+n)<p-4;n++) {
	  double s = Beta_0ijk*Psi_l[0]*Psi_m[0]*Psi_n[0];
	  PHI_v[3*3*jj + 3*0 + 0] = s; PHI_v[3*3*jj + 3*0 + 1] = 0; PHI_v[3*3*jj + 3*0 + 2] = 0;
	  PHI_v[3*3*jj + 3*0 + 0] = 0; PHI_v[3*3*jj + 3*1 + 1] = s; PHI_v[3*3*jj + 3*0 + 2] = 0;
	  PHI_v[3*3*jj + 3*0 + 0] = 0; PHI_v[3*3*jj + 3*0 + 1] = 0; PHI_v[3*3*jj + 3*2 + 2] = s;
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
