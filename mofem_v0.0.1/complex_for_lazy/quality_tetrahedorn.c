/* Copyright (C) 2009, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
 */

/* This file is part of mofem.
 * mofem is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * mofem is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with mofem. If not, see <http://www.gnu.org/licenses/>. */

static PetscErrorCode ierr;

PetscErrorCode calculate_lrms(double *dofs_egdes_X,double *lrms) {
  PetscFunctionBegin;
  *lrms = 0;
  int ee = 0;
  //loop over edges
  for(;ee<6;ee++) { 
    double edge_relative[3];
    bzero(edge_relative,3*sizeof(double));
    int jj = 0;
    for(;jj<3;jj++) {
      edge_relative[jj] = dofs_egdes_X[6*ee+jj];
      edge_relative[jj] -= dofs_egdes_X[6*ee+3+jj];  
      *lrms += pow(edge_relative[jj],2.); 
    }
  }
  *lrms = sqrt( (1./6.)*(*lrms) ); 
  PetscFunctionReturn(0);
}
PetscErrorCode calculate_push_edge_relative(double *edge_coords,__CLPK_doublecomplex *xH,__CLPK_doublecomplex *dofs_egdes_X,__CLPK_doublecomplex *xlrms) {
  PetscFunctionBegin;
  bzero(dofs_egdes_X,3*6*sizeof(__CLPK_doublecomplex));
  double complex lrms = 0; 
  int ee = 0;
  //loop over edges
  for(;ee<6;ee++) { 
    int jj = 0;
    if(xH!=NULL) {
      __CLPK_doublecomplex edge_relative[3];
      bzero(edge_relative,3*sizeof(__CLPK_doublecomplex));
      for(;jj<3;jj++) {
	edge_relative[jj].r = edge_coords[6*ee+jj];
	edge_relative[jj].r -= edge_coords[6*ee+3+jj]; }
      __CLPK_doublecomplex tmp1 = {1.,0.},tmp2 = {0.,0.};
      cblas_zgemv(CblasRowMajor,CblasNoTrans,3,3,&tmp1,xH,3,edge_relative,1,&tmp2,&dofs_egdes_X[3*ee],1); 
    }
    else {
      for(;jj<3;jj++) {
	dofs_egdes_X[3*ee+jj].r = edge_coords[6*ee+jj];
	dofs_egdes_X[3*ee+jj].r -= edge_coords[6*ee+3+jj];
      }}
    jj = 0;
    for(;jj<3;jj++) {
      lrms += cpow(dofs_egdes_X[3*ee+jj].r+I*dofs_egdes_X[3*ee+jj].i,2.); 
    } 
  }
  lrms = csqrt( (1./6.)*lrms ); 
  (*xlrms).r = creal(lrms);
  (*xlrms).i = cimag(lrms); 
  PetscFunctionReturn(0);
}
PetscErrorCode calculate_xQ(
  __CLPK_doublecomplex *xlrms,__CLPK_doublecomplex *inv_xH/*gradient of defomeation*/,
  __CLPK_doublecomplex *dofs_egdes_X_1,__CLPK_doublecomplex *dofs_egdes_X_2,__CLPK_doublecomplex *xQ) {
  PetscFunctionBegin;
  //add transpose inv_xH
  int nn = 0;
  for(;nn<3;nn++) {	
    int mm = 0;
    for(;mm<3;mm++) {
      xQ[nn*3+mm].r = inv_xH[mm*3+nn].r;
      xQ[nn*3+mm].i = inv_xH[mm*3+nn].i; 
  }}
  //loop over edges
  double complex tmp0 = -1./(2.*cpow(((*xlrms).r+I*(*xlrms).i),2.)); 
  int ee = 0;
  for(;ee<6;ee++) {
    int nn = 0;
    for(;nn<3;nn++) {
      complex double tmp1 = dofs_egdes_X_1[3*ee+nn].r+I*dofs_egdes_X_1[3*ee+nn].i;
      int mm = 0;
      for(;mm<3;mm++) {
	complex double tmp2 = tmp0*tmp1*(dofs_egdes_X_2[3*ee+mm].r+I*dofs_egdes_X_2[3*ee+mm].i);
	xQ[mm*3+nn].r += creal(tmp2);
	xQ[mm*3+nn].i += cimag(tmp2); }
    }}
  PetscFunctionReturn(0);
}
static int qual_ver = 1;
void set_qual_ver(int ver) { qual_ver = ver; }
int get_qual_ver() { return qual_ver; }

#define QUALITY_VOLUME_LENGTH \
  __CLPK_doublecomplex xlrms; \
  __CLPK_doublecomplex dofs_egdes_X[3*6]; \
  ierr = calculate_push_edge_relative(coords_edges,xH,dofs_egdes_X,&xlrms); CHKERRQ(ierr); \
  /* calculeate xQ */ \
  __CLPK_doublecomplex xQ[9]; \
  ierr = calculate_xQ(&xlrms,inv_xH,dofs_egdes_Chi,dofs_egdes_X,xQ); CHKERRQ(ierr); \
  /* some useful varibles */ \
  __CLPK_doublecomplex xV; \
  xV.r = det_xH.r*V; \
  xV.i = det_xH.i*V; \
  complex double complex_b = (det_xH.r+I*det_xH.i)/cpow((xlrms.r+I*xlrms.i)/xlrms0.r,3); \
  __CLPK_doublecomplex xb = { creal(complex_b), cimag(complex_b) }; \
  complex double complex_q = 6.*sqrt(2.)*(xV.r+I*xV.i)/cpow(xlrms.r+I*xlrms.i,3.); \
  __CLPK_doublecomplex xq = { creal(complex_q), cimag(complex_q) }; \
  complex double complex_grad; \
  if( qual_ver == 0 ) { \
    /* barrier and quality gradient change */ \
    complex_grad = (xb.r+I*xb.i)/(1.-gamma)-1./((xb.r+I*xb.i)-gamma); \
  } \
  if( qual_ver == 1 ) { \
    /* barrier and quality gradient */ \
    complex_grad = (xq.r+I*xq.i)/(1-gamma)-1./((xq.r+I*xq.i)-gamma); \
  } \
  if( qual_ver == 2 ) { \
    /* quality gradient */ \
    complex_grad = xq.r+I*xq.i; \
  } \
  __CLPK_doublecomplex xgrad = { creal(complex_grad), cimag(complex_grad) }; \
  cblas_zscal(9,&xgrad,xQ,1); 

PetscErrorCode quality_volume_length_F(double V,double *alpha2,double gamma,double *diffN,
  double *coords_edges,double *dofs_X,double *dofs_x,double *dofs_iX,double *dofs_ix,double *quality0,double *quality,double *b,
  double *F,double *iF) {
  PetscFunctionBegin;
  double N4[4*4]; 
  ierr = ShapeMBTET(N4,G_TET_X4,G_TET_Y4,G_TET_Z4,4); CHKERRQ(ierr);
  double H[9];
  ierr = GradientOfDeformation(diffN,dofs_X,H);  CHKERRQ(ierr);
  double iH[9];
  if(dofs_iX != NULL) {
    ierr = GradientOfDeformation(diffN,dofs_iX,iH);    CHKERRQ(ierr);
  } else {
    bzero(iH,9*sizeof(double));
  } 
  __CLPK_doublecomplex xH[9];
  ierr = MakeComplexTensor(H,iH,xH); CHKERRQ(ierr);
  __CLPK_doublecomplex det_xH;
  ierr = DeterminantComplexGradient(xH,&det_xH); CHKERRQ(ierr);
  ierr = MakeComplexTensor(H,iH,xH); CHKERRQ(ierr);
  __CLPK_doublecomplex inv_xH[9];
  cblas_zcopy(9,xH,1,inv_xH,1);
  ierr = InvertComplexGradient(inv_xH); CHKERRQ(ierr);
  __CLPK_doublecomplex xlrms0; 
  __CLPK_doublecomplex dofs_egdes_Chi[3*6]; 
  ierr = calculate_push_edge_relative(coords_edges,NULL,dofs_egdes_Chi,&xlrms0);  CHKERRQ(ierr);
  QUALITY_VOLUME_LENGTH
  //printf("%12.10e %12.10e %12.10e\n",xlrms.r-xlrms0.r,xb.r-1,det_xH.r-1);
  //print_mat_complex(xQ,3,3); 
  //print_mat_complex(xH,3,3); 
  //print_mat_complex(inv_xH,3,3); 
  *quality0 = 6.*sqrt(2.)*(V)/pow(xlrms0.r,3); 
  *quality = xq.r;
  *b = xb.r;
  if( F==NULL ) { 
    PetscFunctionReturn(0);
  }
 
  double reQ[9];
  TakeRe(xQ,reQ);
  double imQ[9];
  TakeIm(xQ,imQ);
  int gg = 0;
  for(;gg<4;gg++) {
    double alpha2_val = cblas_ddot(4,&N4[4*gg],1,alpha2,1);
    int node = 0;
    for(;node<4;node++) {
      if(F!=NULL) {
	F[3*node + 0] += G_TET_W4[gg]*alpha2_val*cblas_ddot(3,&diffN[node*3+0],1,&reQ[0],1);
	F[3*node + 1] += G_TET_W4[gg]*alpha2_val*cblas_ddot(3,&diffN[node*3+0],1,&reQ[3],1);
	F[3*node + 2] += G_TET_W4[gg]*alpha2_val*cblas_ddot(3,&diffN[node*3+0],1,&reQ[6],1); 
      } 
      if(iF!=NULL) {
	iF[3*node + 0] += G_TET_W4[gg]*alpha2_val*cblas_ddot(3,&diffN[node*3+0],1,&imQ[0],1);
	iF[3*node + 1] += G_TET_W4[gg]*alpha2_val*cblas_ddot(3,&diffN[node*3+0],1,&imQ[3],1);
	iF[3*node + 2] += G_TET_W4[gg]*alpha2_val*cblas_ddot(3,&diffN[node*3+0],1,&imQ[6],1); 
      }
  }}
  PetscFunctionReturn(0);
}
int quality_volume_length_K(double eps,double V,double *alpha2,double gamma,double *diffN,double *coords_edges,double *dofs_X,double *dofs_x,double *K,double *Koff) {
  double N4[4*4]; 
  ierr = ShapeMBTET(N4,G_TET_X4,G_TET_Y4,G_TET_Z4,4); CHKERRQ(ierr);
  double H[9];
  ierr = GradientOfDeformation(diffN,dofs_X,H);  CHKERRQ(ierr);
  double ZERO[9];
  bzero(ZERO,sizeof(double)*9);
  __CLPK_doublecomplex xlrms0; 
  __CLPK_doublecomplex dofs_egdes_Chi[3*6]; 
  ierr = calculate_push_edge_relative(coords_edges,NULL,dofs_egdes_Chi,&xlrms0);  CHKERRQ(ierr);
  double _idofs_X[12],_iH[9];
  int dd = 0;
  for(;dd<12;dd++) {
    bzero(_idofs_X,sizeof(double)*12);
    _idofs_X[dd] = eps;
    ierr = GradientOfDeformation(diffN,_idofs_X,_iH);  CHKERRQ(ierr);
    __CLPK_doublecomplex xH[9];
    ierr = MakeComplexTensor(H,_iH,xH); CHKERRQ(ierr);
    __CLPK_doublecomplex det_xH;
    ierr = DeterminantComplexGradient(xH,&det_xH);   CHKERRQ(ierr);
    ierr = MakeComplexTensor(H,_iH,xH); CHKERRQ(ierr);
    __CLPK_doublecomplex inv_xH[9];
    cblas_zcopy(9,xH,1,inv_xH,1);
    ierr = InvertComplexGradient(inv_xH); CHKERRQ(ierr);
    QUALITY_VOLUME_LENGTH
    double imQ[9];
    TakeIm(xQ,imQ);
    cblas_dscal(9,1./eps,imQ,1);
    int gg = 0;
    for(;gg<4;gg++) {
      double alpha2_val = cblas_ddot(4,&N4[4*gg],1,alpha2,1);
      int node = 0;
      for(;node<4;node++) {
	K[3*12*node + 0*12 + dd] += G_TET_W4[gg]*alpha2_val*cblas_ddot(3,&diffN[node*3+0],1,&imQ[0],1);
	K[3*12*node + 1*12 + dd] += G_TET_W4[gg]*alpha2_val*cblas_ddot(3,&diffN[node*3+0],1,&imQ[3],1);
	K[3*12*node + 2*12 + dd] += G_TET_W4[gg]*alpha2_val*cblas_ddot(3,&diffN[node*3+0],1,&imQ[6],1); 
  }}}
  return 0;
}


