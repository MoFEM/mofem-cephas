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

#include <fem_tools.h>
#include <gm_rule.h>

#include <h1_hdiv_hcurl_l2.h>

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
PetscErrorCode Grundmann_Moeller_integration_points_1D_EDGE(int rule,double *G_TRI_X,double *G_TRI_W) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  int dim_num=1;
  int point;
  int point_num;
  double *w;
  double *x;
	
//  GM_RULE_SET determines the weights and abscissas
//  pof a Grundmann-Moeller quadrature rule for
//  the DIM_NUM dimensional simplex,
//  using a rule of in index RULE,
//	  which will have degree of exactness 2*RULE+1.
	
//  printf ( "  Here we use DIM_NUM = %d\n", dim_num  );
//  printf ( "  RULE = %d\n", rule );
//  printf ( "  DEGREE = %d\n", 2 * rule + 1 );
	
  point_num = gm_rule_size ( rule, dim_num );
	
  ierr = PetscMalloc(point_num*sizeof(double),&w); CHKERRQ(ierr);
  ierr = PetscMalloc(dim_num*point_num*sizeof(double),&x); CHKERRQ(ierr);
	
  gm_rule_set ( rule, dim_num, point_num, w, x );
	
  for( point = 0; point < point_num; point++ ){
    G_TRI_X[point] = x[0+point*dim_num];
    G_TRI_W[point] = w[point];
  }

  ierr = PetscFree(w); CHKERRQ(ierr);
  ierr = PetscFree(x); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode Grundmann_Moeller_integration_points_2D_TRI(int rule,double *G_TRI_X,double *G_TRI_Y,double *G_TRI_W){
  PetscFunctionBegin;

  PetscErrorCode ierr;

  int dim_num=2;
  int point;
  int point_num;
  double *w;
  double *x;
	
//  GM_RULE_SET determines the weights and abscissas
//  pof a Grundmann-Moeller quadrature rule for
//  the DIM_NUM dimensional simplex,
//  using a rule of in index RULE,
//	  which will have degree of exactness 2*RULE+1.
	
//  printf ( "  Here we use DIM_NUM = %d\n", dim_num  );
//  printf ( "  RULE = %d\n", rule );
//  printf ( "  DEGREE = %d\n", 2 * rule + 1 );
	
  point_num = gm_rule_size ( rule, dim_num );
	
  ierr = PetscMalloc(point_num*sizeof(double),&w); CHKERRQ(ierr);
  ierr = PetscMalloc(dim_num*point_num*sizeof(double),&x); CHKERRQ(ierr);
	
  gm_rule_set ( rule, dim_num, point_num, w, x );
	
  for( point = 0; point < point_num; point++ ){
      G_TRI_X[point] = x[0+point*dim_num];
      G_TRI_Y[point] = x[1+point*dim_num];
      G_TRI_W[point] = w[point];
  }

  ierr = PetscFree(w); CHKERRQ(ierr);
  ierr = PetscFree(x); CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

PetscErrorCode Grundmann_Moeller_integration_points_3D_TET(int rule,double *G_TET_X,double *G_TET_Y,double *G_TET_Z, double *G_TET_W){
  PetscFunctionBegin;
	
  PetscErrorCode ierr;

  int dim_num=3;
  int point;
  int point_num;
  double *w;
  double *x;
	
//	printf ( "  Here we use DIM_NUM = %d\n", dim_num  );
//	printf ( "  RULE = %d\n", rule );
//	printf ( "  DEGREE = %d\n", 2 * rule + 1 );
		
  point_num = gm_rule_size ( rule, dim_num );

  ierr = PetscMalloc(point_num*sizeof(double),&w); CHKERRQ(ierr);
  ierr = PetscMalloc(dim_num*point_num*sizeof(double),&x); CHKERRQ(ierr);

  gm_rule_set ( rule, dim_num, point_num, w, x );
  for(point = 0; point < point_num; point++ ){
      G_TET_X[point] = x[0+point*dim_num];
      G_TET_Y[point] = x[1+point*dim_num];
      G_TET_Z[point] = x[2+point*dim_num];
      G_TET_W[point] = w[point];
  }

  ierr = PetscFree(w); CHKERRQ(ierr);
  ierr = PetscFree(x); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode ShapeMBTRI(double *N,const double *X,const double *Y,const int G_DIM) {
  PetscFunctionBegin;
  int ii = 0;
  for(; ii<G_DIM; ii++) {
    double x = X[ii],y = Y[ii];
    N[3*ii+0] = N_MBTRI0(x,y);
    N[3*ii+1] = N_MBTRI1(x,y);
    N[3*ii+2] = N_MBTRI2(x,y);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ShapeDiffMBTRI(double *diffN) {
  PetscFunctionBegin;
  diffN[0] = diffN_MBTRI0x; diffN[1] = diffN_MBTRI0y;
  diffN[2] = diffN_MBTRI1x; diffN[3] = diffN_MBTRI1y;
  diffN[4] = diffN_MBTRI2x; diffN[5] = diffN_MBTRI2y;
  PetscFunctionReturn(0);
}
PetscErrorCode ShapeFaceBaseMBTRI(
  double *diffN,const double *coords,
  double *normal,double *s1,double *s2) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  double diffX_ksi[3];
  double diffX_eta[3];
  int ii = 0;
  for(; ii<3; ii++) {
    diffX_ksi[ii] = cblas_ddot(3,&coords[ii],3,&diffN[0],2);
    diffX_eta[ii] = cblas_ddot(3,&coords[ii],3,&diffN[1],2);
  }
  if(s1 != NULL) {
    cblas_dcopy(3,diffX_ksi,1,s1,1);
  }
  if(s2 != NULL) {
    cblas_dcopy(3,diffX_eta,1,s2,1);
  }
  double Spin_diffX_ksi[9];
  ierr = Spin(Spin_diffX_ksi,diffX_ksi); CHKERRQ(ierr);
  cblas_dgemv(
    CblasRowMajor,CblasNoTrans,3,3,
    1.,Spin_diffX_ksi,3,diffX_eta,1,0.,normal,1);
  PetscFunctionReturn(0);
}

PetscErrorCode ShapeFaceNormalMBTRI(
  double *diffN,const double *coords,double *normal) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = ShapeFaceBaseMBTRI(diffN,coords,normal,NULL,NULL);  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode ShapeFaceDiffNormal_MBTRI(double *diffN,const double *coords,double *diff_normal) {
  PetscFunctionBegin;
  // N = Spin(dX/dksi)*dX/deta = -Spin(dX/deta)*dX/dksi
  PetscErrorCode ierr;
  double diffX_ksi[3];
  double diffX_eta[3];
  int ii = 0;
  for(; ii<3; ii++) {
    diffX_ksi[ii] = cblas_ddot(3,&coords[ii],3,&diffN[0],2);
    diffX_eta[ii] = cblas_ddot(3,&coords[ii],3,&diffN[1],2);
  }
  double Spin_diffX_ksi[9];
  ierr = Spin(Spin_diffX_ksi,diffX_ksi); CHKERRQ(ierr);
  double Spin_diffX_eta[9];
  ierr = Spin(Spin_diffX_eta,diffX_eta); CHKERRQ(ierr);
  double B_ksi[3*9];
  bzero(B_ksi,3*9*sizeof(double));
  double B_eta[3*9];
  bzero(B_eta,3*9*sizeof(double));
  //B_ksi[] = [
  // diffN[2*0+0],	0,	0,	diffN[2*1+0],	0,	0,	diffN[2*2+0],	0, 	0
  // 0,		diffN[2*0+0],	0,	0,	diffN[2*1+0],	0,	0,	diffN[2*2+0],	0
  // 0,		0,	diffM[2*0+0],	0,	0,	diffN[2*1+0],	0,	0,	diffN[2*2+0]	
  //]
  //B_eta[] = [
  // diffN[2*0+1],	0,	0,	diffN[2*1+1],	0,	0,	diffN[2*2+1],	0, 	0
  // 0,		diffN[2*0+1],	0,	0,	diffN[2*1+1],	0,	0,	diffN[2*2+1],	0
  // 0,		0,	diffM[2*0+1],	0,	0,	diffN[2*1+1],	0,	0,	diffN[2*2+1]	
  //]
  ii = 0;
  for(;ii<3;ii++) {
    cblas_dcopy(3,&diffN[0],2,&B_ksi[ii*9+ii],3);
    cblas_dcopy(3,&diffN[1],2,&B_eta[ii*9+ii],3);
  }
  cblas_dgemm(
    CblasRowMajor,CblasNoTrans,CblasNoTrans,3,9,3,
      +1.,Spin_diffX_ksi,3,B_eta,9,0.,diff_normal,9);
  cblas_dgemm(
    CblasRowMajor,CblasNoTrans,CblasNoTrans,3,9,3,
      -1.,Spin_diffX_eta,3,B_ksi,9,1.,diff_normal,9);
  PetscFunctionReturn(0);
}

//MBTET
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
PetscErrorCode ShapeMBTET(double *N,const double *G_X,const double *G_Y,const double *G_Z,int DIM) {
  PetscFunctionBegin;
  int ii = 0;
  for(; ii<DIM; ii++) {
    double x = G_X[ii],y = G_Y[ii],z = G_Z[ii];
    N[4*ii+0] = N_MBTET0(x,y,z);
    N[4*ii+1] = N_MBTET1(x,y,z);
    N[4*ii+2] = N_MBTET2(x,y,z);
    N[4*ii+3] = N_MBTET3(x,y,z);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ShapeDiffMBTET(double *diffN) {
  PetscFunctionBegin;
  diffN[0] = diffN_MBTET0x; diffN[1] = diffN_MBTET0y; diffN[2] = diffN_MBTET0z;
  diffN[3] = diffN_MBTET1x; diffN[4] = diffN_MBTET1y; diffN[5] = diffN_MBTET1z;
  diffN[6] = diffN_MBTET2x; diffN[7] = diffN_MBTET2y; diffN[8] = diffN_MBTET2z;
  diffN[9] = diffN_MBTET3x; diffN[10] = diffN_MBTET3y; diffN[11] = diffN_MBTET3z;
  PetscFunctionReturn(0);
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
PetscErrorCode GradientOfDeformation(double *diffN,double *dofs,double *F) {
  PetscFunctionBegin;
  int col,row = 0;
  for(;row<3;row++)
  for(col = 0;col<3;col++) {
    F[3*row+col] = cblas_ddot(4,&diffN[col],3,&dofs[row],3); 
  }
  PetscFunctionReturn(0);
}

// Approximation
PetscErrorCode Lagrange_basis(int p,double s,double *diff_s,double *L,double *diffL,const int dim) {
  PetscFunctionBegin;
  if(dim < 1) SETERRQ(PETSC_COMM_SELF,1,"dim < 1");
  if(dim > 3) SETERRQ(PETSC_COMM_SELF,1,"dim > 3");
  if(p<0) SETERRQ(PETSC_COMM_SELF,1,"p < 0");
  L[0] = 1;
  if(diffL!=NULL) {
    diffL[0*(p+1)+0] = 0;
    if(dim >= 2) {
      diffL[1*(p+1)+0] = 0;
      if(dim == 3) diffL[2*(p+1)+0] = 0;
    }
  }
  if(p==0) PetscFunctionReturn(0);
  L[1] = s;
  if(diffL != NULL) {
    if(diff_s == NULL) {
      SETERRQ(PETSC_COMM_SELF,1,"diff_s == NULL");
    }
    diffL[0*(p+1)+1] = diff_s[0];
    if(dim >= 2) {
      diffL[1*(p+1)+1] = diff_s[1];
      if(dim == 3) diffL[2*(p+1)+1] = diff_s[2];
    }
  }
  if(p==1) PetscFunctionReturn(0);
  int l = 1;
  for(;l<p;l++) {
    double A = ( (2*(double)l+1)/((double)l+1) );
    double B = ( (double)l/((double)l+1) );
    L[l+1] = A*s*L[l] - B*L[l-1]; 
    if(diffL!=NULL) {
      if(diff_s==NULL) {
	SETERRQ(PETSC_COMM_SELF,1,"diff_s == NULL");
      }
      diffL[0*(p+1)+l+1] = A*(s*diffL[0*(p+1)+l] + diff_s[0]*L[l]) - B*diffL[0*(p+1)+l-1]; 
      if(dim >= 2) {
	diffL[1*(p+1)+l+1] = A*(s*diffL[1*(p+1)+l] + diff_s[1]*L[l]) - B*diffL[1*(p+1)+l-1]; 
	if(dim == 3) diffL[2*(p+1)+l+1] = A*(s*diffL[2*(p+1)+l] + diff_s[2]*L[l]) - B*diffL[2*(p+1)+l-1];
      }
    }
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
PetscErrorCode ShapeFaceNormalMBTRI_complex(double *diffN,__CLPK_doublecomplex *xcoords,__CLPK_doublecomplex *xnormal) {
  PetscFunctionBegin;
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
  PetscFunctionReturn(0);
}
PetscErrorCode MakeComplexTensor(double *reA,double *imA,__CLPK_doublecomplex *xA) {
  PetscFunctionBegin;
  int ii = 0,jj;
  for(;ii<3;ii++) {
    for(jj=0;jj<3;jj++) {
      xA[3*ii+jj].r = reA[3*ii+jj];
      xA[3*ii+jj].i = imA[3*ii+jj];
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode InvertComplexGradient(__CLPK_doublecomplex *xF) {
  PetscFunctionBegin;
  __CLPK_integer IPIV[4];
  __CLPK_doublecomplex WORK[4];
  __CLPK_integer LWORK = 4;
  __CLPK_integer info;
  info = lapack_zgetrf(3,3,xF,3,IPIV);
  if(info != 0) SETERRQ(PETSC_COMM_SELF,1,"info == 0");
  info = lapack_zgetri(3,xF,3,IPIV,WORK,LWORK);
  if(info != 0) SETERRQ(PETSC_COMM_SELF,1,"info == 0");
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
  if(lapack_zgetrf(3,3,xF,3,IPIV) != 0) {
    SETERRQ(PETSC_COMM_SELF,1,"lapack_zgetrf(3,3,xF,3,IPIV) != 0");
  }
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
PetscErrorCode Spin(double *spinOmega,double *vecOmega) {
  PetscFunctionBegin;
  bzero(spinOmega,9*sizeof(double));
  spinOmega[0*3+1] = -vecOmega[2];
  spinOmega[0*3+2] = +vecOmega[1];
  spinOmega[1*3+0] = +vecOmega[2];
  spinOmega[1*3+2] = -vecOmega[0];
  spinOmega[2*3+0] = -vecOmega[1];
  spinOmega[2*3+1] = +vecOmega[0]; 
  PetscFunctionReturn(0);
}
PetscErrorCode make_complex_matrix(double *reA,double *imA,__CLPK_doublecomplex *xA) {
  PetscFunctionBegin;
  int ii = 0,jj;
  for(;ii<3;ii++) {
    for(jj=0;jj<3;jj++) {
      xA[3*ii+jj].r = reA[3*ii+jj];
      xA[3*ii+jj].i = imA[3*ii+jj];
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Normal_hierarchical(
  int order_approx,int *order_edge_approx,
  int order,int *order_edge,
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *dofs,double *dofs_edge[],double *dofs_face,
  double *idofs,double *idofs_edge[],double *idofs_face,
    __CLPK_doublecomplex *xnormal,
    __CLPK_doublecomplex *xs1,
    __CLPK_doublecomplex *xs2,
    int gg) {
  PetscFunctionBegin;
  int nn,ee,dd;
  //node
  double complex diffX_x_node,diffX_y_node,diffX_z_node;
  double complex diffY_x_node,diffY_y_node,diffY_z_node;
  diffX_x_node = 0.;diffX_y_node = 0.;diffX_z_node = 0.;
  diffY_x_node = 0.;diffY_y_node = 0.;diffY_z_node = 0.;
  if(dofs!=NULL || idofs != NULL) {
    nn = 0;
    for(; nn<3; nn++) {
      if(dofs!=NULL) {
	diffX_x_node += dofs[3*nn + 0]*diffN[2*nn+0];
	diffX_y_node += dofs[3*nn + 1]*diffN[2*nn+0];
	diffX_z_node += dofs[3*nn + 2]*diffN[2*nn+0];
	diffY_x_node += dofs[3*nn + 0]*diffN[2*nn+1];
	diffY_y_node += dofs[3*nn + 1]*diffN[2*nn+1];
	diffY_z_node += dofs[3*nn + 2]*diffN[2*nn+1]; 
      }
      if(idofs!=NULL) {
	diffX_x_node += I*idofs[3*nn + 0]*diffN[2*nn+0];
	diffX_y_node += I*idofs[3*nn + 1]*diffN[2*nn+0];
	diffX_z_node += I*idofs[3*nn + 2]*diffN[2*nn+0];
	diffY_x_node += I*idofs[3*nn + 0]*diffN[2*nn+1];
	diffY_y_node += I*idofs[3*nn + 1]*diffN[2*nn+1];
	diffY_z_node += I*idofs[3*nn + 2]*diffN[2*nn+1]; 
      }
    }
  }
  double complex diffX_x,diffX_y,diffX_z;
  double complex diffY_x,diffY_y,diffY_z;
  diffX_x = diffX_x_node;diffX_y = diffX_y_node;diffX_z = diffX_z_node;
  diffY_x = diffY_x_node;diffY_y = diffY_y_node;diffY_z = diffY_z_node;
  if(dofs_face!=NULL ||idofs_face!=NULL) {
    int nb_dofs_face = NBFACE_H1(order);
    int nb_dofs_approx_face = NBFACE_H1(order_approx);
    if(nb_dofs_face>0) {
      if(dofs_face!=NULL) {
	diffX_x += cblas_ddot(nb_dofs_face,&dofs_face[0],3,&diffN_face[gg*2*nb_dofs_approx_face+0],2);
	diffX_y += cblas_ddot(nb_dofs_face,&dofs_face[1],3,&diffN_face[gg*2*nb_dofs_approx_face+0],2);
	diffX_z += cblas_ddot(nb_dofs_face,&dofs_face[2],3,&diffN_face[gg*2*nb_dofs_approx_face+0],2);
	diffY_x += cblas_ddot(nb_dofs_face,&dofs_face[0],3,&diffN_face[gg*2*nb_dofs_approx_face+1],2);
	diffY_y += cblas_ddot(nb_dofs_face,&dofs_face[1],3,&diffN_face[gg*2*nb_dofs_approx_face+1],2);
	diffY_z += cblas_ddot(nb_dofs_face,&dofs_face[2],3,&diffN_face[gg*2*nb_dofs_approx_face+1],2); 
      }
      if(idofs_face!=NULL) {
	diffX_x += I*cblas_ddot(nb_dofs_face,&idofs_face[0],3,&diffN_face[gg*2*nb_dofs_approx_face+0],2);
	diffX_y += I*cblas_ddot(nb_dofs_face,&idofs_face[1],3,&diffN_face[gg*2*nb_dofs_approx_face+0],2);
	diffX_z += I*cblas_ddot(nb_dofs_face,&idofs_face[2],3,&diffN_face[gg*2*nb_dofs_approx_face+0],2);
	diffY_x += I*cblas_ddot(nb_dofs_face,&idofs_face[0],3,&diffN_face[gg*2*nb_dofs_approx_face+1],2);
	diffY_y += I*cblas_ddot(nb_dofs_face,&idofs_face[1],3,&diffN_face[gg*2*nb_dofs_approx_face+1],2);
	diffY_z += I*cblas_ddot(nb_dofs_face,&idofs_face[2],3,&diffN_face[gg*2*nb_dofs_approx_face+1],2); 
      }
    }
  }
  ee = 0;
  if(dofs_edge!=NULL || idofs_edge!=NULL) {
    for(;ee<3;ee++) {
      int nb_dofs_edge = NBEDGE_H1(order_edge[ee]);
      int nb_dofs_approx_edge = NBEDGE_H1(order_edge_approx[ee]);
      if(nb_dofs_edge>0) {
	if(dofs_edge!=NULL) {
	  if(dofs_edge[ee]!=NULL) {
	    diffX_x += cblas_ddot(nb_dofs_edge,&(dofs_edge[ee])[0],3,&(diffN_edge[ee])[gg*2*nb_dofs_approx_edge+0],2);
	    diffX_y += cblas_ddot(nb_dofs_edge,&(dofs_edge[ee])[1],3,&(diffN_edge[ee])[gg*2*nb_dofs_approx_edge+0],2);
	    diffX_z += cblas_ddot(nb_dofs_edge,&(dofs_edge[ee])[2],3,&(diffN_edge[ee])[gg*2*nb_dofs_approx_edge+0],2);
	    diffY_x += cblas_ddot(nb_dofs_edge,&(dofs_edge[ee])[0],3,&(diffN_edge[ee])[gg*2*nb_dofs_approx_edge+1],2);
	    diffY_y += cblas_ddot(nb_dofs_edge,&(dofs_edge[ee])[1],3,&(diffN_edge[ee])[gg*2*nb_dofs_approx_edge+1],2);
	    diffY_z += cblas_ddot(nb_dofs_edge,&(dofs_edge[ee])[2],3,&(diffN_edge[ee])[gg*2*nb_dofs_approx_edge+1],2); 
	  }
	}
	if(idofs_edge!=NULL) {
	  if(idofs_edge[ee]==NULL) continue;
	  diffX_x += I*cblas_ddot(nb_dofs_edge,&(idofs_edge[ee])[0],3,&(diffN_edge[ee])[gg*2*nb_dofs_approx_edge+0],2);
	  diffX_y += I*cblas_ddot(nb_dofs_edge,&(idofs_edge[ee])[1],3,&(diffN_edge[ee])[gg*2*nb_dofs_approx_edge+0],2);
	  diffX_z += I*cblas_ddot(nb_dofs_edge,&(idofs_edge[ee])[2],3,&(diffN_edge[ee])[gg*2*nb_dofs_approx_edge+0],2);
	  diffY_x += I*cblas_ddot(nb_dofs_edge,&(idofs_edge[ee])[0],3,&(diffN_edge[ee])[gg*2*nb_dofs_approx_edge+1],2);
	  diffY_y += I*cblas_ddot(nb_dofs_edge,&(idofs_edge[ee])[1],3,&(diffN_edge[ee])[gg*2*nb_dofs_approx_edge+1],2);
	  diffY_z += I*cblas_ddot(nb_dofs_edge,&(idofs_edge[ee])[2],3,&(diffN_edge[ee])[gg*2*nb_dofs_approx_edge+1],2); 
	}
      }
    }
  }
  double complex normal[3];
  normal[0] = diffX_y*diffY_z - diffX_z*diffY_y;
  normal[1] = diffX_z*diffY_x - diffX_x*diffY_z;
  normal[2] = diffX_x*diffY_y - diffX_y*diffY_x;
  dd = 0;
  for(;dd<3;dd++) {
    xnormal[dd].r = creal(normal[dd]);
    xnormal[dd].i = cimag(normal[dd]);
  }
  if(xs1!=NULL) {
    xs1[0].r = creal(diffX_x); xs1[0].i = cimag(diffX_x);
    xs1[1].r = creal(diffX_y); xs1[1].i = cimag(diffX_y);
    xs1[2].r = creal(diffX_z); xs1[2].i = cimag(diffX_z);
  }
  if(xs2!=NULL) {
    xs2[0].r = creal(diffY_x); xs2[0].i = cimag(diffY_x);
    xs2[1].r = creal(diffY_y); xs2[1].i = cimag(diffY_y);
    xs2[2].r = creal(diffY_z); xs2[2].i = cimag(diffY_z);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Base_scale(
  __CLPK_doublecomplex *xnormal,__CLPK_doublecomplex *xs1,__CLPK_doublecomplex *xs2) {
  PetscFunctionBegin;
  complex double xnrm2_normal = csqrt( 
      cpow(xnormal[0].r+I*xnormal[0].i,2) +
      cpow(xnormal[1].r+I*xnormal[1].i,2) +
      cpow(xnormal[2].r+I*xnormal[2].i,2) );
  int dd = 0;
  for(;dd<3;dd++) {
      complex double s1 = (xs1[dd].r+I*xs1[dd].i)*xnrm2_normal;
      complex double s2 = (xs2[dd].r+I*xs2[dd].i)*xnrm2_normal;
      xs1[dd].r = creal(s1); 
      xs1[dd].i = cimag(s1);
      xs2[dd].r = creal(s2);
      xs2[dd].i = cimag(s2);
    }
  PetscFunctionReturn(0);
}

//MBEDGE
PetscErrorCode ShapeMBEDGE(double *N,const double *G_X,int DIM) {
  PetscFunctionBegin;
  int ii = 0;
  for(; ii<DIM; ii++) {
    double x = G_X[ii];
    N[2*ii+0] = N_MBEDGE0(x);
    N[2*ii+1] = N_MBEDGE1(x);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ShapeDiffMBEDGE(double *diffN) {
  PetscFunctionBegin;
  diffN[0] = diffN_MBEDGE0;
  diffN[1] = diffN_MBEDGE1;
  PetscFunctionReturn(0);
}

//FIXME: NOT PROPERLY TESTED YET
//HO 
//MBTRIQ
#define N_MBTRIQ0(x, y) ( (1.-x-y)*(2*(1.-x-y)-1.) )
#define N_MBTRIQ1(x, y) ( x*(2.*x-1.) )
#define N_MBTRIQ2(x, y) ( y*(2.*y-1.) )
#define N_MBTRIQ3(x, y) ( 4.*(1.-x-y)*x )
#define N_MBTRIQ4(x, y) ( 4.*x*y )
#define N_MBTRIQ5(x, y) ( 4.*(1.-x-y)*y )
#define diffN_MBTRIQ0x(x, y) ( x+y-3.*(1.-x-y) )
#define diffN_MBTRIQ0y(x, y) ( x+y-3.*(1.-x-y) )
#define diffN_MBTRIQ1x(x, y) ( -1.+4.*x )
#define diffN_MBTRIQ1y(x, y) ( 0. )
#define diffN_MBTRIQ2x(x, y) ( 0. )
#define diffN_MBTRIQ2y(x, y) ( -1.+4.*y )
#define diffN_MBTRIQ3x(x, y) ( 4.*((1.-x-y)-x) )
#define diffN_MBTRIQ3y(x, y) ( -4.*x )
#define diffN_MBTRIQ4x(x, y) ( 4.*y )
#define diffN_MBTRIQ4y(x, y) ( 4.*x )
#define diffN_MBTRIQ5x(x, y) ( -4.*y )
#define diffN_MBTRIQ5y(x, y) ( 4.*((1.-x-y)-y) )
PetscErrorCode ShapeMBTRIQ_GAUSS(double *N,const double *X,const double *Y,const int G_DIM) {
  PetscFunctionBegin;
  int ii = 0;
  for(; ii<G_DIM; ii++) {
    double x = X[ii],y = Y[ii];
    N[6*ii+0] = N_MBTRIQ0(x,y);
    N[6*ii+1] = N_MBTRIQ1(x,y);
    N[6*ii+2] = N_MBTRIQ2(x,y);
    N[6*ii+3] = N_MBTRIQ3(x,y);
    N[6*ii+4] = N_MBTRIQ4(x,y);
    N[6*ii+5] = N_MBTRIQ5(x,y);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ShapeMBTRIQ(double *N,const double x,const double y) {
  PetscFunctionBegin;
  N[0] = N_MBTRIQ0(x,y);
  N[1] = N_MBTRIQ1(x,y);
  N[2] = N_MBTRIQ2(x,y);
  N[3] = N_MBTRIQ3(x,y);
  N[4] = N_MBTRIQ4(x,y);
  N[5] = N_MBTRIQ5(x,y);
  PetscFunctionReturn(0);
}
PetscErrorCode ShapeDiffMBTRIQ(double *diffN,const double x,const double y) {
  PetscFunctionBegin;
  diffN[0] = diffN_MBTRIQ0x(x,y);
  diffN[1] = diffN_MBTRIQ0y(x,y);
  diffN[2] = diffN_MBTRIQ1x(x,y);
  diffN[3] = diffN_MBTRIQ1y(x,y);
  diffN[4] = diffN_MBTRIQ2x(x,y);
  diffN[5] = diffN_MBTRIQ2y(x,y);
  diffN[6] = diffN_MBTRIQ3x(x,y);
  diffN[7] = diffN_MBTRIQ3y(x,y);
  diffN[8] = diffN_MBTRIQ4x(x,y);
  diffN[9] = diffN_MBTRIQ4y(x,y);
  diffN[10] = diffN_MBTRIQ5x(x,y);
  diffN[11] = diffN_MBTRIQ5y(x,y);
  PetscFunctionReturn(0);
}

//MBTETQ (JULIEN WORK)
#define N_MBTETQ0(x, y, z) ( (2.*(1.-x-y-z)-1.)*(1.-x-y-z) )
#define N_MBTETQ1(x, y, z) ( (2.*x-1.)*x )
#define N_MBTETQ2(x, y, z) ( (2.*y-1.)*y )
#define N_MBTETQ3(x, y, z) ( (2.*z-1.)*z )
#define N_MBTETQ4(x, y, z) ( 4.*(1.-x-y-z)*x )
#define N_MBTETQ5(x, y, z) ( 4.*x*y )
#define N_MBTETQ6(x, y, z) ( 4.*(1.-x-y-z)*y )
#define N_MBTETQ7(x, y, z) ( 4.*(1.-x-y-z)*z )
#define N_MBTETQ8(x, y, z) ( 4.*x*z )
#define N_MBTETQ9(x, y, z) ( 4.*y*z )
#define diffN_MBTETQ0x(x, y, z) ( -3.+4.*x+4.*y+4.*z )
#define diffN_MBTETQ0y(x, y, z) ( -3.+4.*x+4.*y+4.*z )
#define diffN_MBTETQ0z(x, y, z) ( -3.+4.*x+4.*y+4.*z )
#define diffN_MBTETQ1x(x, y, z) ( 4.*x-1. )
#define diffN_MBTETQ1y(x, y, z) ( 0. )
#define diffN_MBTETQ1z(x, y, z) ( 0. )
#define diffN_MBTETQ2x(x, y, z) ( 0. )
#define diffN_MBTETQ2y(x, y, z) ( 4.*y-1. )
#define diffN_MBTETQ2z(x, y, z) ( 0. )
#define diffN_MBTETQ3x(x, y, z) ( 0. )
#define diffN_MBTETQ3y(x, y, z) ( 0. )
#define diffN_MBTETQ3z(x, y, z) ( 4.*z-1. )
#define diffN_MBTETQ4x(x, y, z) ( -8.*x+4.-4.*y-4.*z )
#define diffN_MBTETQ4y(x, y, z) ( -4.*x )
#define diffN_MBTETQ4z(x, y, z) ( -4.*x )
#define diffN_MBTETQ5x(x, y, z) ( 4.*y )
#define diffN_MBTETQ5y(x, y, z) ( 4.*x)
#define diffN_MBTETQ5z(x, y, z) ( 0. )
#define diffN_MBTETQ6x(x, y, z) ( -4.*y )
#define diffN_MBTETQ6y(x, y, z) ( -8.*y+4.-4.*x-4.*z )
#define diffN_MBTETQ6z(x, y, z) ( -4.*y )
#define diffN_MBTETQ7x(x, y, z) ( -4.*z )
#define diffN_MBTETQ7y(x, y, z) ( -4.*z )
#define diffN_MBTETQ7z(x, y, z) ( -8.*z+4.-4.*x-4.*y )
#define diffN_MBTETQ8x(x, y, z) ( 4.*z )
#define diffN_MBTETQ8y(x, y, z) ( 0. )
#define diffN_MBTETQ8z(x, y, z) ( 4.*x )
#define diffN_MBTETQ9x(x, y, z) ( 0. )
#define diffN_MBTETQ9y(x, y, z) ( 4.*z )
#define diffN_MBTETQ9z(x, y, z) ( 4.*y )
PetscErrorCode ShapeMBTETQ(double *N,const double x,const double y,const double z) {
  PetscFunctionBegin;
  N[0] = N_MBTETQ0(x,y,z);
  N[1] = N_MBTETQ1(x,y,z);
  N[2] = N_MBTETQ2(x,y,z);
  N[3] = N_MBTETQ3(x,y,z);
  N[4] = N_MBTETQ4(x,y,z);
  N[5] = N_MBTETQ5(x,y,z);
  N[6] = N_MBTETQ6(x,y,z);
  N[7] = N_MBTETQ7(x,y,z);
  N[8] = N_MBTETQ8(x,y,z);
  N[9] = N_MBTETQ9(x,y,z);
  PetscFunctionReturn(0);
}
PetscErrorCode ShapeDiffMBTETQ(double *diffN,const double x,const double y,const double z) {
  PetscFunctionBegin;
  diffN[0] = diffN_MBTETQ0x(x,y,z);
  diffN[1] = diffN_MBTETQ0y(x,y,z);
  diffN[2] = diffN_MBTETQ0z(x,y,z);
  diffN[3] = diffN_MBTETQ1x(x,y,z);
  diffN[4] = diffN_MBTETQ1y(x,y,z);
  diffN[5] = diffN_MBTETQ1z(x,y,z);
  diffN[6] = diffN_MBTETQ2x(x,y,z);
  diffN[7] = diffN_MBTETQ2y(x,y,z);
  diffN[8] = diffN_MBTETQ2z(x,y,z);
  diffN[9] = diffN_MBTETQ3x(x,y,z);
  diffN[10] = diffN_MBTETQ3y(x,y,z);
  diffN[11] = diffN_MBTETQ3z(x,y,z);
  diffN[12] = diffN_MBTETQ4x(x,y,z);
  diffN[13] = diffN_MBTETQ4y(x,y,z);
  diffN[14] = diffN_MBTETQ4z(x,y,z);
  diffN[15] = diffN_MBTETQ5x(x,y,z);
  diffN[16] = diffN_MBTETQ5y(x,y,z);
  diffN[17] = diffN_MBTETQ5z(x,y,z);
  diffN[18] = diffN_MBTETQ6x(x,y,z);
  diffN[19] = diffN_MBTETQ6y(x,y,z);
  diffN[20] = diffN_MBTETQ6z(x,y,z);
  diffN[21] = diffN_MBTETQ7x(x,y,z);
  diffN[22] = diffN_MBTETQ7y(x,y,z);
  diffN[23] = diffN_MBTETQ7z(x,y,z);
  diffN[24] = diffN_MBTETQ8x(x,y,z);
  diffN[25] = diffN_MBTETQ8y(x,y,z);
  diffN[26] = diffN_MBTETQ8z(x,y,z);
  diffN[27] = diffN_MBTETQ9x(x,y,z);
  diffN[28] = diffN_MBTETQ9y(x,y,z);
  diffN[29] = diffN_MBTETQ9z(x,y,z);
  PetscFunctionReturn(0);
}
PetscErrorCode ShapeMBTETQ_GAUSS(double *N,const double *X,const double *Y,const double *Z,const int G_DIM) {
  PetscFunctionBegin;
  int ii = 0;
  for(; ii<G_DIM; ii++) {
    double x = X[ii],y = Y[ii],z = Z[ii];
    N[10*ii+0] = N_MBTETQ0(x,y,z);
    N[10*ii+1] = N_MBTETQ1(x,y,z);
    N[10*ii+2] = N_MBTETQ2(x,y,z);
    N[10*ii+3] = N_MBTETQ3(x,y,z);
    N[10*ii+4] = N_MBTETQ4(x,y,z);
    N[10*ii+5] = N_MBTETQ5(x,y,z);
    N[10*ii+6] = N_MBTETQ6(x,y,z);
    N[10*ii+7] = N_MBTETQ7(x,y,z);
    N[10*ii+8] = N_MBTETQ8(x,y,z);
    N[10*ii+9] = N_MBTETQ9(x,y,z);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ShapeDiffMBTETQ_GAUSS(double *diffN,const double *X,const double *Y,const double *Z,const int G_DIM) {
  PetscFunctionBegin;
  int ii = 0;
  for(; ii<G_DIM; ii++) {
    double x = X[ii],y = Y[ii],z = Z[ii];
    diffN[30*ii+0] = diffN_MBTETQ0x(x,y,z); diffN[30*ii+1] = diffN_MBTETQ0y(x,y,z); diffN[30*ii+2] = diffN_MBTETQ0z(x,y,z);
    diffN[30*ii+3] = diffN_MBTETQ1x(x,y,z); diffN[30*ii+4] = diffN_MBTETQ1y(x,y,z); diffN[30*ii+5] = diffN_MBTETQ1z(x,y,z);
    diffN[30*ii+6] = diffN_MBTETQ2x(x,y,z); diffN[30*ii+7] = diffN_MBTETQ2y(x,y,z); diffN[30*ii+8] = diffN_MBTETQ2z(x,y,z);
    diffN[30*ii+9] = diffN_MBTETQ3x(x,y,z); diffN[30*ii+10] = diffN_MBTETQ3y(x,y,z); diffN[30*ii+11] = diffN_MBTETQ3z(x,y,z);
    diffN[30*ii+12] = diffN_MBTETQ4x(x,y,z); diffN[30*ii+13] = diffN_MBTETQ4y(x,y,z); diffN[30*ii+14] = diffN_MBTETQ4z(x,y,z);
    diffN[30*ii+15] = diffN_MBTETQ5x(x,y,z); diffN[30*ii+16] = diffN_MBTETQ5y(x,y,z); diffN[30*ii+17] = diffN_MBTETQ5z(x,y,z);
    diffN[30*ii+18] = diffN_MBTETQ6x(x,y,z); diffN[30*ii+19] = diffN_MBTETQ6y(x,y,z); diffN[30*ii+20] = diffN_MBTETQ6z(x,y,z);
    diffN[30*ii+21] = diffN_MBTETQ7x(x,y,z); diffN[30*ii+22] = diffN_MBTETQ7y(x,y,z); diffN[30*ii+23] = diffN_MBTETQ7z(x,y,z);
    diffN[30*ii+24] = diffN_MBTETQ8x(x,y,z); diffN[30*ii+25] = diffN_MBTETQ8y(x,y,z); diffN[30*ii+26] = diffN_MBTETQ8z(x,y,z);
    diffN[30*ii+27] = diffN_MBTETQ9x(x,y,z); diffN[30*ii+28] = diffN_MBTETQ9y(x,y,z); diffN[30*ii+29] = diffN_MBTETQ9z(x,y,z);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ShapeJacMBTETQ(const double *diffN,const double *coords,double *Jac) {
  PetscFunctionBegin;
  int ii,jj,kk;
  bzero(Jac,sizeof(double)*9);
  for(ii = 0; ii<10; ii++) 	//shape func.
    for(jj = 0; jj<3; jj++) 	//space
      for(kk = 0; kk<3; kk++) 	//direvative of shape func.
	Jac[ jj*3+kk ] += 
	  diffN[ ii*3+kk ]*coords[ ii*3+jj ];
  PetscFunctionReturn(0);
}
PetscErrorCode ShapeMBTETQ_detJac_at_Gauss_Points(double *detJac_at_Gauss_Points,const double *diffN,const double *coords,int G_DIM) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  double Jac[9];
  int ii = 0;
  for(; ii<G_DIM; ii++) {
    ierr = ShapeJacMBTETQ(&diffN[30*ii],coords,Jac); CHKERRQ(ierr);
    detJac_at_Gauss_Points[ii] = Shape_detJac(Jac);
  }
  PetscFunctionReturn(0);
}
double Shape_intVolumeMBTETQ(const double *diffN,const double *coords,int G_DIM,double *G_TET_W) {
  PetscErrorCode ierr;
  int ii = 0;
  double vol = 0;
  double detJac_at_Gauss_Points[G_DIM];
  ierr = ShapeMBTETQ_detJac_at_Gauss_Points(detJac_at_Gauss_Points,diffN,coords,G_DIM); CHKERRQ(ierr);
  for(; ii<G_DIM; ii++) {
    vol += G_TET_W[ii]*(detJac_at_Gauss_Points[ii])/6;
  }
  return vol;
}
/*PetscErrorCode ShapeMBTETQ_inverse(ShapeData *data,const double *elem_coords,const double *glob_coords,double *loc_coords,const double eps) {
  PetscFunctionBegin;
  double A[3*3];  
  double R[3];  
  double *N = (*data).N;
  double *diffN = (*data).diffN;
  int IPIV[3];
  float NORM_dR = 1000.;
  float NORM_R0;
  ShapeMBTETQ(data,0.1,0.1,0.1);
  ShapeDiffMBTETQ(data,0.1,0.1,0.1);
  R[0] = glob_coords[0] - cblas_ddot(10,&N[0],1,&elem_coords[0],3);
  R[1] = glob_coords[1] - cblas_ddot(10,&N[0],1,&elem_coords[1],3);
  R[2] = glob_coords[2] - cblas_ddot(10,&N[0],1,&elem_coords[2],3);
  NORM_R0 = cblas_dnrm2(3,&R[0],1);
  while ( (NORM_dR/NORM_R0)>eps){
    //COL MAJOR
    //X
    A[0+3*0] = cblas_ddot(10,&diffN[0*3+0],3,&elem_coords[0],3);
    A[0+3*1] = cblas_ddot(10,&diffN[0*3+1],3,&elem_coords[0],3);
    A[0+3*2] = cblas_ddot(10,&diffN[0*3+2],3,&elem_coords[0],3);
    R[0] = glob_coords[0] - cblas_ddot(10,&N[0],1,&elem_coords[0],3);
    //Y 
    A[1+3*0] = cblas_ddot(10,&diffN[0*3+0],3,&elem_coords[1],3);
    A[1+3*1] = cblas_ddot(10,&diffN[0*3+1],3,&elem_coords[1],3);
    A[1+3*2] = cblas_ddot(10,&diffN[0*3+2],3,&elem_coords[1],3);
    R[1] = glob_coords[1] - cblas_ddot(10,&N[0],1,&elem_coords[1],3);
    //Z
    A[2+3*0] = cblas_ddot(10,&diffN[0*3+0],3,&elem_coords[0*3+2],3);
    A[2+3*1] = cblas_ddot(10,&diffN[0*3+1],3,&elem_coords[0*3+2],3);
    A[2+3*2] = cblas_ddot(10,&diffN[0*3+2],3,&elem_coords[0*3+2],3);
    R[2] = glob_coords[2] - cblas_ddot(10,&N[0],1,&elem_coords[2],3);
    assert( lapack_dgesv(3,1,&A[0],3,(__CLPK_integer*)IPIV,R,3) == 0 );
    cblas_daxpy(3,1.,R,1,loc_coords,1);
    NORM_dR = cblas_dnrm2(3,&R[0],1);
    ShapeMBTETQ(data,loc_coords[0],loc_coords[1],loc_coords[2]);
    ShapeDiffMBTETQ(data,loc_coords[0],loc_coords[1],loc_coords[2]);
 }
 PetscFunctionReturn(0);
}
*/
