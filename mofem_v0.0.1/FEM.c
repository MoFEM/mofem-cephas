/* Copyright (C) 2009, Lukasz Kaczmarczyk (likask AT civil.gla.ac.uk)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <string.h>

#include<FEM.h>

#define EPS 1.e-12
void print_mat(double *M,int m,int n) {
  int ii,jj;
  for(ii = 0; ii<m; ii++) {
    jj = 0;
    for(; jj<n; jj++) {
      double val = M[ii*n+jj];
      if( fabs(val) < EPS ) val = 0; 
      printf("%+8.7e ",val);
    }; printf("\n");
  }; 
  fflush(stdout);
}
void print_mat_sym_upper(double *M,int m,int n) {
  int ii,jj;
  for(ii = 0; ii<m; ii++) {
    jj = 0;
    for(; jj<n; jj++) {
      double val = M[ii*n+jj];
      if( jj<ii ) val = M[jj*n+ii];
      if( fabs(val) < EPS ) val = 0; 
      printf("%+8.7e ",val);
    }; printf("\n");
  }; 
  fflush(stdout);
}
void print_mat_complex(__CLPK_doublecomplex *M,int m,int n) {
  int ii,jj;
  for(ii = 0; ii<m; ii++) {
    jj = 0;
    for(; jj<n; jj++) {
      double val_re = M[ii*n+jj].r;
      double val_im = M[ii*n+jj].i;
      printf("%+8.7e%+8.7ei ",val_re,val_im);
    }; printf("\n");
  }; 
  fflush(stdout);
}

double Shape_detJac(double *Jac) {
  double detJac;
  __CLPK_integer IPIV[4];
  __CLPK_integer info = lapack_dgetrf(3,3,Jac,3,IPIV);
  if(info !=0) return -1;
  int i = 0,j = 0;
  detJac = 1.;
  for(; i<3; i++) {
    detJac *= Jac[3*i+i];
    if( IPIV[i] != i+1 ) j++;
  }
  if ( (j - ( j/2 )*2) != 0 ) 
    detJac = - detJac;
  return detJac;
}
PetscErrorCode Shape_invJac(double *Jac) {
  PetscFunctionBegin;
  __CLPK_integer IPIV[4];
  __CLPK_doublereal WORK[3];
  __CLPK_integer LWORK = 3;
  __CLPK_integer info;
  info = lapack_dgetrf(3,3,Jac,3,IPIV);
  if(info != 0) SETERRQ1(PETSC_COMM_SELF,1,"info = %d",info);
  info = lapack_dgetri(3,Jac,3,IPIV,WORK,LWORK);
  if(info != 0) SETERRQ1(PETSC_COMM_SELF,1,"info = %d",info);
  PetscFunctionReturn(0);
}

//MBTRI
#define N_MBTRI0(x, y) ( 1.-x-y )
#define N_MBTRI1(x, y) ( x )
#define N_MBTRI2(x, y) ( y )
#define diffN_MBTRI0x ( -1. )
#define diffN_MBTRI0y ( -1. )
#define diffN_MBTRI1x ( 1 )
#define diffN_MBTRI1y ( 0 )
#define diffN_MBTRI2x ( 0 )
#define diffN_MBTRI2y ( 1 )
void ShapeMBTRI_GAUSS(double *N,const double *X,const double *Y,const int G_DIM) {
  int ii = 0;
  for(; ii<G_DIM; ii++) {
    double x = X[ii],y = Y[ii];
    N[3*ii+0] = N_MBTRI0(x,y);
    N[3*ii+1] = N_MBTRI1(x,y);
    N[3*ii+2] = N_MBTRI2(x,y);
  }
}
void ShapeDiffMBTRI(double *diffN) {
  diffN[0] = diffN_MBTRI0x; diffN[1] = diffN_MBTRI0y;
  diffN[2] = diffN_MBTRI1x; diffN[3] = diffN_MBTRI1y;
  diffN[4] = diffN_MBTRI2x; diffN[5] = diffN_MBTRI2y;
}
void ShapeFaceNormalMBTRI(double *diffN,const double *coords,double *normal) {
  double diffX_x,diffX_y,diffX_z;
  double diffY_x,diffY_y,diffY_z;
  diffX_x = 0.;diffX_y = 0.;diffX_z = 0.;
  diffY_x = 0.;diffY_y = 0.;diffY_z = 0.;
  int ii = 0;
  for(; ii<3; ii++) {
    diffX_x += coords[3*ii + 0]*diffN[2*ii+0];
    diffX_y += coords[3*ii + 1]*diffN[2*ii+0];
    diffX_z += coords[3*ii + 2]*diffN[2*ii+0];
    diffY_x += coords[3*ii + 0]*diffN[2*ii+1];
    diffY_y += coords[3*ii + 1]*diffN[2*ii+1];
    diffY_z += coords[3*ii + 2]*diffN[2*ii+1];
  }
  normal[0] = diffX_y*diffY_z - diffX_z*diffY_y;
  normal[1] = diffX_z*diffY_x - diffX_x*diffY_z;
  normal[2] = diffX_x*diffY_y - diffX_y*diffY_x;
  cblas_dscal(3,0.5,normal,1);
}
void ShapeJacMBTRI(double *diffN,const double *coords,double *Jac) {
  int ii,jj,kk;
  bzero(Jac,sizeof(double)*9);
  for(ii = 0; ii<3; ii++) 	//shape func.
    for(jj = 0; jj<3; jj++) 	//space
      for(kk = 0; kk<3; kk++) 	//direvative of shape func.
	Jac[ jj*3+kk ] += 
	diffN[ ii*3+kk ]*coords[ ii*3+jj ];
}
void ShapeDiffMBTRIinvJ(double *diffN,double *invJac,double *diffNinvJac) {
  int ii = 0;
  for(;ii<3; ii++) {
   cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,invJac,3,&diffN[ii*3],1,0.,&diffNinvJac[ii*3],1);
  }
}

//MBTET
#define N_MBTET0(x, y, z) ( 1.-x-y-z )
#define N_MBTET1(x, y, z) ( x )
#define N_MBTET2(x, y, z) ( y )
#define N_MBTET3(x, y, z) ( z )
#define diffN_MBTET0x ( -1. )
#define diffN_MBTET0y ( -1. )
#define diffN_MBTET0z ( -1. )
#define diffN_MBTET1x ( 1 )
#define diffN_MBTET1y ( 0 )
#define diffN_MBTET1z ( 0 )
#define diffN_MBTET2x ( 0 )
#define diffN_MBTET2y ( 1 )
#define diffN_MBTET2z ( 0 )
#define diffN_MBTET3x ( 0 )
#define diffN_MBTET3y ( 0 )
#define diffN_MBTET3z ( 1 )

PetscErrorCode ShapeJacMBTET(double *diffN,const double *coords,double *Jac) {
  PetscFunctionBegin;
  int ii,jj,kk;
  bzero(Jac,sizeof(double)*9);
  for(ii = 0; ii<4; ii++) 	//shape func.
    for(jj = 0; jj<3; jj++) 	//space
      for(kk = 0; kk<3; kk++) 	//direvative of shape func.
	Jac[ jj*3+kk ] += 
	diffN[ ii*3+kk ]*coords[ ii*3+jj ];
  PetscFunctionReturn(0);
}
double Shape_intVolumeMBTET(double *diffN,const double *coords) {
  double Jac[9];
  ShapeJacMBTET(diffN,coords,Jac);
  double detJac = Shape_detJac(Jac);
  //printf("detJac = +%6.4e\n",detJac);
  //print_mat(Jac,3,3);
  return detJac*G_TET_W1[0]/6.;
}
void ShapeMBTET(double *N,const double *G_X,const double *G_Y,const double *G_Z,int DIM) {
  int ii = 0;
  for(; ii<DIM; ii++) {
    double x = G_X[ii],y = G_Y[ii],z = G_Z[ii];
    N[4*ii+0] = N_MBTET0(x,y,z);
    N[4*ii+1] = N_MBTET1(x,y,z);
    N[4*ii+2] = N_MBTET2(x,y,z);
    N[4*ii+3] = N_MBTET3(x,y,z);
  }
}
void ShapeDiffMBTET(double *diffN) {
  diffN[0] = diffN_MBTET0x; diffN[1] = diffN_MBTET0y; diffN[2] = diffN_MBTET0z;
  diffN[3] = diffN_MBTET1x; diffN[4] = diffN_MBTET1y; diffN[5] = diffN_MBTET1z;
  diffN[6] = diffN_MBTET2x; diffN[7] = diffN_MBTET2y; diffN[8] = diffN_MBTET2z;
  diffN[9] = diffN_MBTET3x; diffN[10] = diffN_MBTET3y; diffN[11] = diffN_MBTET3z;
}
PetscErrorCode ShapeMBTET_inverse(double *N,double *diffN,const double *elem_coords,const double *glob_coords,double *loc_coords) {
  PetscFunctionBegin;
  double A[3*3];  
  double R[3];  
  int IPIV[3];
  //COL MAJOR
  //X
  A[0+3*0] = cblas_ddot(4,&diffN[0*3+0],3,&elem_coords[0*3+0],3);
  A[0+3*1] = cblas_ddot(4,&diffN[0*3+1],3,&elem_coords[0*3+0],3);
  A[0+3*2] = cblas_ddot(4,&diffN[0*3+2],3,&elem_coords[0*3+0],3);
  R[0] = glob_coords[0] - cblas_ddot(4,&N[0],1,&elem_coords[0*3+0],3);
  //printf("A\n[ %3.2f %3.2f %3.2f ] %3.2f \n",A[0*3],A[1*3],A[2*3],R[0]);
  //Y 
  A[1+3*0] = cblas_ddot(4,&diffN[0*3+0],3,&elem_coords[0*3+1],3);
  A[1+3*1] = cblas_ddot(4,&diffN[0*3+1],3,&elem_coords[0*3+1],3);
  A[1+3*2] = cblas_ddot(4,&diffN[0*3+2],3,&elem_coords[0*3+1],3);
  R[1] = glob_coords[1] - cblas_ddot(4,&N[0],1,&elem_coords[0*3+1],3);
  //printf("[ %3.2f %3.2f %3.2f ] %3.2f \n",A[1+3*0],A[1+3*1],A[1+3*2],R[1]);
  //Z
  A[2+3*0] = cblas_ddot(4,&diffN[0*3+0],3,&elem_coords[0*3+2],3);
  A[2+3*1] = cblas_ddot(4,&diffN[0*3+1],3,&elem_coords[0*3+2],3);
  A[2+3*2] = cblas_ddot(4,&diffN[0*3+2],3,&elem_coords[0*3+2],3);
  R[2] = glob_coords[2] - cblas_ddot(4,&N[0],1,&elem_coords[0*3+2],3);
  //printf("[ %3.2f %3.2f %3.2f ] %3.2f \n",A[2+3*0],A[2+3*1],A[2+3*2],R[1]);
  int info = lapack_dgesv(3,1,&A[0],3,(__CLPK_integer*)IPIV,R,3);
  if(info != 0) SETERRQ1(PETSC_COMM_SELF,1,"info == %d",info);
  //assert( info == 0 );
  cblas_dcopy(3,R,1,loc_coords,1);
  PetscFunctionReturn(0);
}
PetscErrorCode ShapeDiffMBTETinvJ(double *diffN,double *invJac,double *diffNinvJac) {
  PetscFunctionBegin;
  int ii = 0;
  for(;ii<4; ii++) {
   cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,invJac,3,&diffN[ii*3],1,0.,&diffNinvJac[ii*3],1);
  }
  PetscFunctionReturn(0);
}
void GradientOfDeformation(double *diffN,double *dofs,double *F) {
  int col,row = 0;
  for(;row<3;row++)
  for(col = 0;col<3;col++) {
    F[3*row+col] = cblas_ddot(4,&diffN[col],3,&dofs[row],3); }
}

// Approximation
PetscErrorCode Lagrange_basis(int p,double s,double *diff_s,double *L,double *diffL,const int dim) {
  PetscFunctionBegin;
  if(dim < 2) SETERRQ(PETSC_COMM_SELF,1,"dim < 2");
  if(dim > 3) SETERRQ(PETSC_COMM_SELF,1,"dim > 3");
  //assert(dim >= 2);
  //assert(dim <= 3);
  if(p==0) PetscFunctionReturn(0);
  L[0] = 1;
  diffL[0*p+0] = 0;
  diffL[1*p+0] = 0;
  if(dim == 3) diffL[2*p+0] = 0;
  if(p==1) PetscFunctionReturn(0);
  L[1] = s;
  diffL[0*p+1] = diff_s[0];
  diffL[1*p+1] = diff_s[1];
  if(dim == 3) diffL[2*p+1] = diff_s[2];
  if(p==2) PetscFunctionReturn(0);
  int l = 2;
  for(;l<p;l++) {
    double A = ( (2*(double)l-1)/((double)l) );
    double B = ( (double)(l-1)/((double)l) );
    L[l] = A*s*L[l-1] - B*L[l-2]; 
    diffL[0*p+l] = A*(s*diffL[0*p+l-1] + diff_s[0]*L[l-1]) - B*diffL[0*p+l-2]; 
    diffL[1*p+l] = A*(s*diffL[1*p+l-1] + diff_s[1]*L[l-1]) - B*diffL[1*p+l-2]; 
    if(dim == 2) continue;
    diffL[2*p+l] = A*(s*diffL[2*p+l-1] + diff_s[2]*L[l-1]) - B*diffL[2*p+l-2]; 
  }
  PetscFunctionReturn(0);
}

//ALL COMPLEX FROM NOW
void ShapeDiffMBTETinvJ_complex(double *diffN,__CLPK_doublecomplex *invJac,__CLPK_doublecomplex *diffNinvJac,const enum CBLAS_TRANSPOSE Trans) {
  __CLPK_doublecomplex tmp1 = {1.,0.},tmp2 = {0.,0.};
  int ii = 0,jj;
  for(;ii<4; ii++) {
    __CLPK_doublecomplex tmp3[3];
    for(jj = 0;jj<3;jj++) {
      tmp3[jj].r = diffN[ii*3+jj];
      tmp3[jj].i = 0; }
    cblas_zgemv(CblasRowMajor,Trans,3,3,&tmp1,invJac,3,tmp3,1,&tmp2,&diffNinvJac[ii*3],1); }
}
void ShapeFaceNormalMBTRI_complex(double *diffN,__CLPK_doublecomplex *xcoords,__CLPK_doublecomplex *xnormal) {
  double complex diffX_x,diffX_y,diffX_z;
  double complex diffY_x,diffY_y,diffY_z;
  diffX_x = diffX_y = diffX_z = 0.;
  diffY_x = diffY_y = diffY_z = 0.;
  int ii;
  for(ii = 0; ii<3; ii++) {
    diffX_x += (xcoords[3*ii + 0].r+I*xcoords[3*ii + 0].i)*diffN[2*ii+0];
    diffX_y += (xcoords[3*ii + 1].r+I*xcoords[3*ii + 1].i)*diffN[2*ii+0];
    diffX_z += (xcoords[3*ii + 2].r+I*xcoords[3*ii + 2].i)*diffN[2*ii+0];
    diffY_x += (xcoords[3*ii + 0].r+I*xcoords[3*ii + 0].i)*diffN[2*ii+1];
    diffY_y += (xcoords[3*ii + 1].r+I*xcoords[3*ii + 1].i)*diffN[2*ii+1];
    diffY_z += (xcoords[3*ii + 2].r+I*xcoords[3*ii + 2].i)*diffN[2*ii+1];
  }
  double complex tmp;
  tmp = diffX_y*diffY_z - diffX_z*diffY_y;
  xnormal[0].r = creal(tmp);
  xnormal[0].i = cimag(tmp);
  tmp = diffX_z*diffY_x - diffX_x*diffY_z;
  xnormal[1].r = creal(tmp);
  xnormal[1].i = cimag(tmp);
  tmp = diffX_x*diffY_y - diffX_y*diffY_x;
  xnormal[2].r = creal(tmp);
  xnormal[2].i = cimag(tmp);
}
void MakeComplexTensor(double *reA,double *imA,__CLPK_doublecomplex *xA) {
  int ii = 0,jj;
  for(;ii<3;ii++) {
    for(jj=0;jj<3;jj++) {
      xA[3*ii+jj].r = reA[3*ii+jj];
      xA[3*ii+jj].i = imA[3*ii+jj];
    }
  }
}
PetscErrorCode InvertComplexGradient(__CLPK_doublecomplex *xF) {
  PetscFunctionBegin;
  __CLPK_integer IPIV[4];
  __CLPK_doublecomplex WORK[4];
  __CLPK_integer LWORK = 4;
  __CLPK_integer info;
  info = lapack_zgetrf(3,3,xF,3,IPIV);
  if(info == 0) SETERRQ(PETSC_COMM_SELF,1,"info == 0");
  //assert(info == 0);
  info = lapack_zgetri(3,xF,3,IPIV,WORK,LWORK);
  if(info == 0) SETERRQ(PETSC_COMM_SELF,1,"info == 0");
  //assert(info == 0);
  PetscFunctionReturn(0);
}
PetscErrorCode InvertComplexSymmMatrix3by3(__CLPK_doublecomplex *xC) {
  PetscFunctionBegin;
  __CLPK_integer info;
  info = lapack_zpotrf('L',3,xC,3);
  if(info == 0) SETERRQ(PETSC_COMM_SELF,1,"info == 0");
  //assert(info == 0);
  info = lapack_zpotri('L',3,xC,3);
  if(info == 0) SETERRQ(PETSC_COMM_SELF,1,"info == 0");
  //assert(info == 0);
  PetscFunctionReturn(0);
}
PetscErrorCode DeterminantComplexGradient(__CLPK_doublecomplex *xF,__CLPK_doublecomplex *det_xF) {
  PetscFunctionBegin;
  __CLPK_integer IPIV[4];
  if(lapack_zgetrf(3,3,xF,3,IPIV) != 0) SETERRQ(PETSC_COMM_SELF,1,"lapack_zgetrf(3,3,xF,3,IPIV) != 0");
  double complex det = 1;
  int i = 0,j = 0;
  for(; i<3; i++) {
    det *= xF[3*i+i].r + I*xF[3*i+i].i;
    if( IPIV[i] != i+1 ) j++;
  }
  if ( (j - ( j/2 )*2) != 0 ) 
    det = - det;
  (*det_xF).r = creal(det);
  (*det_xF).i = cimag(det);
  PetscFunctionReturn(0);
}


