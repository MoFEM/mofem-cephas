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

PetscErrorCode SpatialGradientOfDeformation(__CLPK_doublecomplex *xh,__CLPK_doublecomplex *inv_xH,__CLPK_doublecomplex *xF) {
  PetscFunctionBegin;
  __CLPK_doublecomplex tmp1 = {1.,0.},tmp2 = {0.,0.};
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,3,3,&tmp1,xh,3,inv_xH,3,&tmp2,xF,3);
  PetscFunctionReturn(0);
}
PetscErrorCode CauchyGreenDeformation(__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC) {
  PetscFunctionBegin;
  //Note: is symmetric
  bzero(xC,sizeof(__CLPK_doublecomplex)*9);
  __CLPK_doublecomplex tmp1 = {1,0},tmp2 = {0,0};
  cblas_zsyrk(CblasRowMajor,CblasUpper,CblasTrans,3,3,&tmp1,xF,3,&tmp2,xC,3);
  int ii=0,jj;
  for(;ii<3;ii++) {
    for(jj=0;jj<ii;jj++) {
      xC[ii*3+jj].r = xC[jj*3+ii].r; 
      xC[ii*3+jj].i = xC[jj*3+ii].i; 
  }}
  PetscFunctionReturn(0);
}
PetscErrorCode TakeIm(__CLPK_doublecomplex *xA,double *imA) {
  PetscFunctionBegin;
  int jj,ii = 0;
  for(;ii<3;ii++) {
    for(jj=0;jj<3;jj++) {
      imA[ii*3+jj] = xA[ii*3+jj].i;
    } } 
  PetscFunctionReturn(0);
}
PetscErrorCode TakeRe(__CLPK_doublecomplex *xA,double *reA) {
  PetscFunctionBegin;
  int jj,ii = 0;
  for(;ii<3;ii++) {
    for(jj=0;jj<3;jj++) {
      reA[ii*3+jj] = xA[ii*3+jj].r; 
    } 
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode PiolaKrihoff1_PullBack(__CLPK_doublecomplex *det_xH,__CLPK_doublecomplex *inv_xH,__CLPK_doublecomplex *xP,__CLPK_doublecomplex *xP_PullBack) {
  PetscFunctionBegin;
  __CLPK_doublecomplex tmp2 = {0,0};
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,3,3,3,det_xH,xP,3,inv_xH,3,&tmp2,xP_PullBack,3);
  PetscFunctionReturn(0);
}
PetscErrorCode ElshebyStress_PullBack(__CLPK_doublecomplex *det_xH,__CLPK_doublecomplex *inv_xH,__CLPK_doublecomplex *xStress,__CLPK_doublecomplex *xStress_PullBack) {
  PetscFunctionBegin;
  __CLPK_doublecomplex tmp2 = {0,0};
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,3,3,3,det_xH,xStress,3,inv_xH,3,&tmp2,xStress_PullBack,3);
  PetscFunctionReturn(0);
}

#define COMP_STRESSES \
  ierr = SpatialGradientOfDeformation(xh,inv_xH,xF); CHKERRQ(ierr); \
  if(dofs_T!=NULL) { \
    double temperature = 0; \
    temperature = cblas_ddot(4,&N[4*gg],1,dofs_T,1); \
    /*fprintf(stdout,"temp %f %f %f %f %f\n",temperature,dofs_T[0],dofs_T[1],dofs_T[2],dofs_T[3]);*/ \
    __CLPK_doublecomplex tmp1 = {1.,0.},tmp2 = {0.,0.}; \
    __CLPK_doublecomplex xT = { temperature, 0 }; \
    __CLPK_doublecomplex inv_xF_themp[9]; \
    ierr = ThermalDeformationGradient(thermal_expansion,thermal_load_factor,i_thermal_load_factor,xT,inv_xF_themp); CHKERRQ(ierr); \
    ierr = InvertComplexGradient(inv_xF_themp); CHKERRQ(ierr); \
    /*{ print_mat_complex(inv_xF_themp,3,3); fprintf(stdout,"inv_xF_themp\n"); }*/ \
    /*print_mat_complex(xF,3,3); fprintf(stdout,"xF\n");*/ \
    cblas_zcopy(9,xF,1,xF_tmp,1); \
    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,3,3,&tmp1,xF_tmp,3,inv_xF_themp,3,&tmp2,xF,3); \
    /*print_mat_complex(xF,3,3); fprintf(stdout,"xF\n");*/ \
  } \
  cblas_zcopy(9,xF,1,xF_tmp,1); \
  ierr = DeterminantComplexGradient(xF_tmp,&det_xF); CHKERRQ(ierr); \
  ierr = CauchyGreenDeformation(xF,xC); CHKERRQ(ierr); \
  ierr = StrainEnergy(lambda,mu,xF,xC,&det_xF,&xPsi,matctx); CHKERRQ(ierr); \
  ierr = PiolaKirhoiff2(lambda,mu,xF,xC,&det_xF,xS,matctx); CHKERRQ(ierr); \
  ierr = PiolaKirhoiff1(lambda,mu,xF,xS,xP); CHKERRQ(ierr); \
  ierr = PiolaKrihoff1_PullBack(&det_xH,inv_xH,xP,xP_PullBack); \
  ierr = ElshebyStress(&xPsi,xF,xP,xSigma); \
  ierr = ElshebyStress_PullBack(&det_xH,inv_xH,xSigma,xSigma_PullBack);

PetscErrorCode HierarhicalDeformationGradient(
  int *order_max_edge,int *order_max_face,int order_max_volume,
  int *order_edge,int *order_face,int order_volume,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_edge[],double *dofs_face[],double *dofs_volume,
  int gg,double *GRAD) {
  PetscFunctionBegin;
  int ee = 0;
  for(;ee<6;ee++) {
    if(NBEDGE_H1(order_edge[ee])==0) continue; 
    double edge_grad[9]; 
    bzero(edge_grad,9*sizeof(double)); 
    double *diff = &((diffN_edge[ee])[gg*3*NBEDGE_H1(order_max_edge[ee])]); 
    ierr = H1_EdgeGradientOfDeformation_hierachical(order_edge[ee],diff,dofs_edge[ee],edge_grad); CHKERRQ(ierr); 
    cblas_daxpy(9,1,edge_grad,1,GRAD,1); } 
  int ff = 0; 
  for(;ff<4;ff++) { 
    if(NBFACE_H1(order_face[ff])==0) continue; 
    double face_grad[9]; 
    bzero(face_grad,9*sizeof(double)); 
    double *diff = &((diffN_face[ff])[gg*3*NBFACE_H1(order_max_face[ff])]); 
    ierr = H1_FaceGradientOfDeformation_hierachical(order_face[ff],diff,dofs_face[ff],face_grad); CHKERRQ(ierr); 
    cblas_daxpy(9,1,face_grad,1,GRAD,1); } 
  if(NBVOLUME_H1(order_volume)>0) { 
    double volume_grad[9]; 
    bzero(volume_grad,9*sizeof(double)); 
    double *diff = &((diffN_volume)[gg*3*NBVOLUME_H1(order_max_volume)]); 
    ierr = H1_VolumeGradientOfDeformation_hierachical(order_volume,diff,dofs_volume,volume_grad); CHKERRQ(ierr); 
    cblas_daxpy(9,1,volume_grad,1,GRAD,1); } 
  PetscFunctionReturn(0);
}

PetscErrorCode Calculate_Stresses_at_GaussPoint(
  int *order_max_edge,int *order_max_face,int order_max_volume,
  int *order_X_edge,int *order_X_face,int order_X_volume,
  int *order_x_edge,int *order_x_face,int order_x_volume,
  double alpha,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X_node,double *dofs_X_edge[],double *dofs_X_face[],double *dofs_X_volume,
  double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  //temperature
  double thermal_expansion,double thermal_load_factor,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //
  double *F,double *Piola1Stress,double *CauhyStress,double *EshelbyStress,double *Psi,double *J,double *themp,
  int gg) {
  PetscFunctionBegin;

  double ZERO[9];
  bzero(ZERO,sizeof(double)*9);

  __CLPK_doublecomplex xh[9],xH[9],inv_xH[9],xF[9],xF_tmp[9],xC[9],xS[9],xP[9],xSigma[9],xCauchyStress[9],xPsi,det_xF,det_xH;
  
  double H[9];
  ierr = GradientOfDeformation(diffN,dofs_X_node,H);  CHKERRQ(ierr);
  ierr = HierarhicalDeformationGradient(
    order_max_edge,order_max_face,order_max_volume,
    order_X_edge,order_X_face,order_X_volume,
    diffN,diffN_edge,diffN_face,diffN_volume,
    dofs_X_edge,dofs_X_face,dofs_X_volume,
    gg,H); CHKERRQ(ierr); 
  ierr = MakeComplexTensor(H,ZERO,xH); CHKERRQ(ierr);
  cblas_zcopy(9,xH,1,inv_xH,1);
  //print_mat_complex(xH,3,3); fprintf(stdout,"xH\n"); 
  ierr = DeterminantComplexGradient(xH,&det_xH); CHKERRQ(ierr);
  ierr = InvertComplexGradient(inv_xH); CHKERRQ(ierr);

  double h[9];
  ierr = GradientOfDeformation(diffN,dofs_x_node,h);  CHKERRQ(ierr);
  ierr = HierarhicalDeformationGradient(
    order_max_edge,order_max_face,order_max_volume,
    order_x_edge,order_x_face,order_x_volume,
    diffN,diffN_edge,diffN_face,diffN_volume,
    dofs_x_edge,dofs_x_face,dofs_x_volume,
    gg,h); CHKERRQ(ierr); 
  ierr = MakeComplexTensor(h,ZERO,xh); CHKERRQ(ierr);
  //print_mat_complex(xh,3,3); fprintf(stdout,"xh\n"); 
  ierr = SpatialGradientOfDeformation(xh,inv_xH,xF); CHKERRQ(ierr); 
  //print_mat_complex(xF,3,3); fprintf(stdout,"xF\n"); 

  //temperature
  *themp = 0;
  if(dofs_T!=NULL) {
    *themp = cblas_ddot(4,&N[4*gg],1,dofs_T,1);
    __CLPK_doublecomplex tmp1 = {1.,0.},tmp2 = {0.,0.};
    __CLPK_doublecomplex xT = { *themp, 0 };
    __CLPK_doublecomplex inv_xF_themp[9];
    ierr = ThermalDeformationGradient(thermal_expansion,thermal_load_factor,0,xT,inv_xF_themp); CHKERRQ(ierr);
    ierr = InvertComplexGradient(inv_xF_themp); CHKERRQ(ierr);
    cblas_zcopy(9,xF,1,xF_tmp,1);
    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,3,3,&tmp1,xF_tmp,3,inv_xF_themp,3,&tmp2,xF,3);
  }

  //print_mat_complex(xF,3,3); fprintf(stdout,"xF themp\n"); 
  cblas_zcopy(9,xF,1,xF_tmp,1);
  //print_mat_complex(xF,3,3); fprintf(stdout,"xF tmp\n"); 
  ierr = DeterminantComplexGradient(xF_tmp,&det_xF); CHKERRQ(ierr); 
  //printf("det_xF %6.4e+%6.4ei\n", det_xF.r,det_xF.i);
  //printf("\n\n\n");

  ierr = CauchyGreenDeformation(xF,xC); CHKERRQ(ierr); 
  ierr = StrainEnergy(lambda,mu,xF,xC,&det_xF,&xPsi,matctx); CHKERRQ(ierr); 

  ierr = PiolaKirhoiff2(lambda,mu,xF,xC,&det_xF,xS,matctx); CHKERRQ(ierr); 
  ierr = PiolaKirhoiff1(lambda,mu,xF,xS,xP); CHKERRQ(ierr); 
  ierr = ElshebyStress(&xPsi,xF,xP,xSigma); 
  ierr = CauchyStress(xF,&det_xF,xP,xCauchyStress); CHKERRQ(ierr);

  *Psi = xPsi.r;
  *J = det_xF.r;

  TakeRe(xF,F);
  TakeRe(xP,Piola1Stress);
  TakeRe(xSigma,EshelbyStress);
  TakeRe(xCauchyStress,CauhyStress);

  PetscFunctionReturn(0);
}
PetscErrorCode bzero_Fint(
  int *order_edge,int *order_face,int order_volume,
  double *Fint,double *Fint_edge[],double *Fint_face[],double *Fint_volume) {
  PetscFunctionBegin;
  if(Fint!=NULL) bzero(Fint,12*sizeof(double));
  if(Fint_edge!=NULL) {
    int ee = 0;
    for(;ee<6;ee++) {
      if(Fint_edge[ee]==NULL) continue;
      bzero(Fint_edge[ee],3*NBEDGE_H1(order_edge[ee])*sizeof(double));
    }
  }
  if(Fint_face!=NULL) {
    int ff = 0;
    for(;ff<4;ff++) {
      if(Fint_face[ff]==NULL) continue;
      bzero(Fint_face[ff],3*NBFACE_H1(order_face[ff])*sizeof(double));
    }
  }
  if(Fint_volume!=NULL) {
    bzero(Fint_volume,3*NBVOLUME_H1(order_volume)*sizeof(double));
  }  
  PetscFunctionReturn(0);
}
PetscErrorCode Fint_Hh_hierarchical(
  int *order_max_edge,int *order_max_face,int order_max_volume,
  int *order_X_edge,int *order_X_face,int order_X_volume,
  int *order_x_edge,int *order_x_face,int order_x_volume,
  double alpha,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X_node,double *dofs_X_edge[],double *dofs_X_face[],double *dofs_X_volume,
  double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  double *dofs_iX,double *dofs_ix_node,
  //temperature
  double thermal_expansion,double thermal_load_factor,double i_thermal_load_factor,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //
  double *Fint_H,double *Fint_H_edge[],double *Fint_H_face[],double *Fint_H_volume,
  double *Fint_h,double *Fint_h_edge[],double *Fint_h_face[],double *Fint_h_volume,
  double *Fint_iH,double *Fint_iH_edge[],double *Fint_iH_face[],double *Fint_iH_volume,
  double *Fint_ih,double *Fint_ih_edge[],double *Fint_ih_face[],double *Fint_ih_volume,
  int G_DIM,const double *G_W) {
  PetscFunctionBegin;
  double ZERO[9];
  bzero(ZERO,sizeof(double)*9);
  __CLPK_doublecomplex xh[9],xH[9],inv_xH[9],xF[9],xF_tmp[9],xC[9],xS[9],xP[9],xP_PullBack[9],xSigma[9],xSigma_PullBack[9],xPsi,det_xF,det_xH;
  double reP[9],reSigma[9],imP[9],imSigma[9];
  ierr = bzero_Fint(order_x_edge,order_x_face,order_x_volume,
    Fint_h,Fint_h_edge,Fint_h_face,Fint_h_volume); CHKERRQ(ierr);
  ierr = bzero_Fint(order_x_edge,order_x_face,order_x_volume,
    Fint_ih,Fint_ih_edge,Fint_h_face,Fint_h_volume); CHKERRQ(ierr);
  ierr = bzero_Fint(order_X_edge,order_X_face,order_X_volume,
    Fint_H,Fint_H_edge,Fint_H_face,Fint_H_volume); CHKERRQ(ierr);
  ierr = bzero_Fint(order_X_edge,order_X_face,order_X_volume,
    Fint_iH,Fint_iH_edge,Fint_iH_face,Fint_iH_volume); CHKERRQ(ierr);
  int gg = 0;
  for(;gg<G_DIM;gg++) {
    //material
    double H[9],iH[9];
    ierr = GradientOfDeformation(diffN,dofs_X_node,H);  CHKERRQ(ierr);
    ierr = HierarhicalDeformationGradient(
      order_max_edge,order_max_face,order_max_volume,
      order_X_edge,order_X_face,order_X_volume,
      diffN,diffN_edge,diffN_face,diffN_volume,
      dofs_X_edge,dofs_X_face,dofs_X_volume,
      gg,H); CHKERRQ(ierr);
    if(dofs_iX == NULL) {
      ierr = MakeComplexTensor(H,ZERO,xH); CHKERRQ(ierr);
    } else {
      ierr = GradientOfDeformation(diffN,dofs_iX,iH); CHKERRQ(ierr);
      ierr = MakeComplexTensor(H,iH,xH); CHKERRQ(ierr);
    }
    cblas_zcopy(9,xH,1,inv_xH,1);
    //print_mat_complex(xH,3,3); fprintf(stdout,"xH\n"); 
    ierr = DeterminantComplexGradient(xH,&det_xH); CHKERRQ(ierr);
    ierr = InvertComplexGradient(inv_xH); CHKERRQ(ierr);
    //spatial
    double h[9],ih[9];
    ierr = GradientOfDeformation(diffN,dofs_x_node,h);  CHKERRQ(ierr);
    if(dofs_ix_node != NULL) {
      ierr = GradientOfDeformation(diffN,dofs_ix_node,ih); CHKERRQ(ierr);
    } else {
      bzero(ih,9*sizeof(double));
    }
    ierr = HierarhicalDeformationGradient(
      order_max_edge,order_max_face,order_max_volume,
      order_x_edge,order_x_face,order_x_volume,
      diffN,diffN_edge,diffN_face,diffN_volume,
      dofs_x_edge,dofs_x_face,dofs_x_volume,
      gg,h); CHKERRQ(ierr); 
    ierr = MakeComplexTensor(h,ih,xh); CHKERRQ(ierr);
    //print_mat_complex(xh,3,3); fprintf(stdout,"xh\n"); 
    COMP_STRESSES
    //print_mat_complex(xF,3,3); fprintf(stdout,"xF\n"); 
    TakeRe(xP_PullBack,reP);
    TakeRe(xSigma_PullBack,reSigma);
    TakeIm(xP_PullBack,imP);
    TakeIm(xSigma_PullBack,imSigma);
    int node = 0;
    for(;node<4;node++) {
      if(Fint_H!=NULL) {
	Fint_H[3*node + 0] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&reSigma[0],1);
	Fint_H[3*node + 1] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&reSigma[3],1);
	Fint_H[3*node + 2] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&reSigma[6],1); }
      if(Fint_h!=NULL) {
	Fint_h[3*node + 0] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&reP[0],1);
	Fint_h[3*node + 1] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&reP[3],1);
	Fint_h[3*node + 2] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&reP[6],1); }
      if(Fint_iH!=NULL) {
	Fint_iH[3*node + 0] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[0],1);
	Fint_iH[3*node + 1] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[3],1);
	Fint_iH[3*node + 2] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[6],1); }
      if(Fint_ih!=NULL) {
	Fint_ih[3*node + 0] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[0],1);
	Fint_ih[3*node + 1] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[3],1);
	Fint_ih[3*node + 2] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[6],1); }
    }
    int ee = 0;
    for(;ee<6;ee++) {
      int pp = 0;
      for(;pp<NBEDGE_H1(order_x_edge[ee]);pp++) {
	double *diff = &((diffN_edge[ee])[gg*3*NBEDGE_H1(order_max_edge[ee])+3*pp]);
	if(Fint_h_edge!=NULL)
	if(Fint_h_edge[ee]!=NULL) {
	  (Fint_h_edge[ee])[3*pp + 0] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reP[0],1);
	  (Fint_h_edge[ee])[3*pp + 1] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reP[3],1);
	  (Fint_h_edge[ee])[3*pp + 2] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reP[6],1); }
	if(Fint_ih_edge!=NULL)
	if(Fint_ih_edge[ee]!=NULL) {
	  (Fint_ih_edge[ee])[3*pp + 0] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	  (Fint_ih_edge[ee])[3*pp + 1] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	  (Fint_ih_edge[ee])[3*pp + 2] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }
      }
      for(;pp<NBEDGE_H1(order_X_edge[ee]);pp++) {
	double *diff = &((diffN_edge[ee])[gg*3*NBEDGE_H1(order_max_edge[ee])+3*pp]);
	if(Fint_H_edge!=NULL)
	if(Fint_H_edge[ee]!=NULL) {
	  (Fint_H_edge[ee])[3*pp + 0] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reSigma[0],1);
	  (Fint_H_edge[ee])[3*pp + 1] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reSigma[3],1);
	  (Fint_H_edge[ee])[3*pp + 2] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reSigma[6],1); }
	if(Fint_iH_edge!=NULL)
	if(Fint_iH_edge[ee]!=NULL) {
	  (Fint_iH_edge[ee])[3*pp + 0] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imSigma[0],1);
	  (Fint_iH_edge[ee])[3*pp + 1] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imSigma[3],1);
	  (Fint_iH_edge[ee])[3*pp + 2] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imSigma[6],1); }
      }
    } 
    int ff = 0;
    for(;ff<4;ff++) {
      int pp = 0;
      for(;pp<NBFACE_H1(order_x_face[ff]);pp++) {
	double *diff = &((diffN_face[ff])[gg*3*NBFACE_H1(order_max_face[ff])+3*pp]);
	if(Fint_h_face!=NULL) 
	if(Fint_h_face[ff]!=NULL) {
	  (Fint_h_face[ff])[3*pp + 0] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reP[0],1);
	  (Fint_h_face[ff])[3*pp + 1] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reP[3],1);
	  (Fint_h_face[ff])[3*pp + 2] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reP[6],1); }
	if(Fint_ih_face!=NULL) 
	if(Fint_ih_face[ff]!=NULL) {
	  (Fint_ih_face[ff])[3*pp + 0] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	  (Fint_ih_face[ff])[3*pp + 1] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	  (Fint_ih_face[ff])[3*pp + 2] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }
      }
      for(;pp<NBFACE_H1(order_X_face[ff]);pp++) {
	double *diff = &((diffN_face[ff])[gg*3*NBFACE_H1(order_max_face[ff])+3*pp]);
	if(Fint_H_face!=NULL) 
	if(Fint_H_face[ff]!=NULL) {
	  (Fint_H_face[ff])[3*pp + 0] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reSigma[0],1);
	  (Fint_H_face[ff])[3*pp + 1] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reSigma[3],1);
	  (Fint_H_face[ff])[3*pp + 2] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reSigma[6],1); }
	if(Fint_iH_face!=NULL) 
	if(Fint_iH_face[ff]!=NULL) {
	  (Fint_iH_face[ff])[3*pp + 0] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imSigma[0],1);
	  (Fint_iH_face[ff])[3*pp + 1] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imSigma[3],1);
	  (Fint_iH_face[ff])[3*pp + 2] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imSigma[6],1); }
      }
    }
    int pp = 0;
    for(;pp<NBVOLUME_H1(order_x_volume);pp++) {
      double *diff = &((diffN_volume)[gg*3*NBVOLUME_H1(order_max_volume)+3*pp]);
      if(Fint_h_volume!=NULL) {
	(Fint_h_volume)[3*pp + 0] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reP[0],1);
	(Fint_h_volume)[3*pp + 1] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reP[3],1);
	(Fint_h_volume)[3*pp + 2] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reP[6],1); }
      if(Fint_ih_volume!=NULL) {
	(Fint_ih_volume)[3*pp + 0] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	(Fint_ih_volume)[3*pp + 1] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	(Fint_ih_volume)[3*pp + 2] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }
    }
    for(;pp<NBVOLUME_H1(order_X_volume);pp++) {
      double *diff = &((diffN_volume)[gg*3*NBVOLUME_H1(order_max_volume)+3*pp]);
      if(Fint_H_volume!=NULL) {
	(Fint_H_volume)[3*pp + 0] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reSigma[0],1);
	(Fint_H_volume)[3*pp + 1] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reSigma[3],1);
	(Fint_H_volume)[3*pp + 2] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&reSigma[6],1); }
      if(Fint_iH_volume!=NULL) {
	(Fint_iH_volume)[3*pp + 0] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imSigma[0],1);
	(Fint_iH_volume)[3*pp + 1] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imSigma[3],1);
	(Fint_iH_volume)[3*pp + 2] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imSigma[6],1); }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Tangent_HH_hierachical(
  int *order_max_edge,int *order_max_face,int order_max_volume,
  int *order_X_edge,int *order_X_face,int order_X_volume,
  int *order_x_edge,int *order_x_face,int order_x_volume,
  double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X_node,double *dofs_X_edge[],double *dofs_X_face[],double *dofs_X_volume,
  double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  //temperature
  double thermal_expansion,double thermal_load_factor,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //
  double *K,double *Koff,double *Koff_edge[6],double *Koff_face[4],double *Koff_volume,int G_DIM,const double *G_W) {
  PetscFunctionBegin;
  double i_thermal_load_factor = 0;
  double ZERO[9];
  bzero(ZERO,sizeof(double)*9);
  __CLPK_doublecomplex xh[9],xH[9],inv_xH[9],xF[9],xF_tmp[9],xC[9],xS[9],xP[9],xP_PullBack[9],xSigma[9],xSigma_PullBack[9],det_xF,det_xH,xPsi;
  double _idofs_X_node[12],_iH[9],imP[9],imSigma[9];
  if(K!=NULL) bzero(K,12*12*sizeof(double));
  if(Koff!=NULL) bzero(Koff,12*12*sizeof(double));
  //zero edges
  int ee = 0;	
  if(Koff_edge!=NULL) {
    for(;ee<6;ee++) {
      if(NBEDGE_H1(order_x_edge[ee])==0) continue;
      int nb = 3*NBEDGE_H1(order_x_edge[ee]);
      if(Koff_edge[ee]!=NULL) bzero(Koff_edge[ee],nb*12*sizeof(double)); 
    }
  }
  //zero faces
  int ff = 0;	
  if(Koff_face!=NULL) {
    for(;ff<4;ff++) {
      if(NBFACE_H1(order_x_face[ff])==0) continue;
      int nb = 3*NBFACE_H1(order_x_face[ff]);
      if(Koff_face[ff]!=NULL) bzero(Koff_face[ff],nb*12*sizeof(double)); 
    }
  }
  //zero volume
  if(Koff_volume!=NULL) {
    if(NBVOLUME_H1(order_x_volume)!=0) {
      int nb = 3*NBVOLUME_H1(order_x_volume);
      if(Koff_volume!=NULL) bzero(Koff_volume,nb*12*sizeof(double)); 
    }
  }
  int dd = 0;
  for(;dd<12;dd++) {
    bzero(_idofs_X_node,sizeof(double)*12);
    _idofs_X_node[dd] = eps;
    int gg = 0;
    for(;gg<G_DIM;gg++) {
      double H[9];
      ierr = GradientOfDeformation(diffN,dofs_X_node,H);  CHKERRQ(ierr);
      ierr = HierarhicalDeformationGradient(
	order_max_edge,order_max_face,order_max_volume,
	order_X_edge,order_X_face,order_X_volume,
	diffN,diffN_edge,diffN_face,diffN_volume,
	dofs_X_edge,dofs_X_face,dofs_X_volume,
	gg,H); CHKERRQ(ierr);
      ierr = GradientOfDeformation(diffN,_idofs_X_node,_iH);  CHKERRQ(ierr);
      ierr = MakeComplexTensor(H,_iH,xH);  CHKERRQ(ierr);
      cblas_zcopy(9,xH,1,inv_xH,1);
      ierr = DeterminantComplexGradient(xH,&det_xH); CHKERRQ(ierr);
      ierr = InvertComplexGradient(inv_xH); CHKERRQ(ierr);
      double h[9];
      ierr = GradientOfDeformation(diffN,dofs_x_node,h);  CHKERRQ(ierr);
      ierr = HierarhicalDeformationGradient(
	order_max_edge,order_max_face,order_max_volume,
	order_x_edge,order_x_face,order_x_volume,
	diffN,diffN_edge,diffN_face,diffN_volume,
	dofs_x_edge,dofs_x_face,dofs_x_volume,
	gg,h); CHKERRQ(ierr); 
      ierr = MakeComplexTensor(h,ZERO,xh);  CHKERRQ(ierr);
      COMP_STRESSES
      TakeIm(xP_PullBack,imP);
      cblas_dscal(9,1./eps,imP,1);
      TakeIm(xSigma_PullBack,imSigma);
      cblas_dscal(9,1./eps,imSigma,1);
      int node = 0;
      for(;node<4;node++) {
	if(K!=NULL) {
	  K[3*12*node + 0*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[0],1);
	  K[3*12*node + 1*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[3],1);
	  K[3*12*node + 2*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[6],1); }
	if(Koff!=NULL) {
	  Koff[3*12*node + 0*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[0],1);
	  Koff[3*12*node + 1*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[3],1);
	  Koff[3*12*node + 2*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[6],1); }}
      ee = 0;
      if(Koff_edge!=NULL) {
	for(;ee<6;ee++) {		
	  int pp = 0;
	  for(;pp<NBEDGE_H1(order_x_edge[ee]);pp++) {
	    double *diff = &((diffN_edge[ee])[gg*3*NBEDGE_H1(order_max_edge[ee])+3*pp]);
	    if(Koff_edge[ee]!=NULL) {
	      (Koff_edge[ee])[3*pp*12 + 0*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	      (Koff_edge[ee])[3*pp*12 + 1*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	      (Koff_edge[ee])[3*pp*12 + 2*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }}}
      }
      ff = 0;
      if(Koff_face!=NULL) {
	for(;ff<4;ff++) {		
	  int pp = 0;
	  for(;pp<NBFACE_H1(order_x_face[ff]);pp++) {
	    double *diff = &((diffN_face[ff])[gg*3*NBFACE_H1(order_max_face[ff])+3*pp]);
	    if(Koff_face[ff]!=NULL) {
	      (Koff_face[ff])[3*pp*12 + 0*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	      (Koff_face[ff])[3*pp*12 + 1*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	      (Koff_face[ff])[3*pp*12 + 2*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }}}
      }
      if(Koff_volume!=NULL) {
	int pp = 0;
	for(;pp<NBVOLUME_H1(order_x_volume);pp++) {
	  double *diff = &((diffN_volume)[gg*3*NBVOLUME_H1(order_max_volume)+3*pp]);
	  if(Koff_volume!=NULL) {
	    (Koff_volume)[3*pp*12 + 0*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	    (Koff_volume)[3*pp*12 + 1*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	    (Koff_volume)[3*pp*12 + 2*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }}
      }
  }}
  PetscFunctionReturn(0);
}
PetscErrorCode Tangent_hh_hierachical(
  int *order_max_edge,int *order_max_face,int order_max_volume,
  int *order_X_edge,int *order_X_face,int order_X_volume,
  int *order_x_edge,int *order_x_face,int order_x_volume,
  double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X_node,double *dofs_X_edge[],double *dofs_X_face[],double *dofs_X_volume,
  double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  //temperature
  double thermal_expansion,double thermal_load_factor,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //
  double *K,double *Koff,double *K_edge[6],double *K_face[4],double *K_volume,int G_DIM,const double *G_W) {
  PetscFunctionBegin;
  double i_thermal_load_factor = 0;
  double ZERO[9];
  bzero(ZERO,sizeof(double)*9);
  __CLPK_doublecomplex xh[9],xH[9],inv_xH[9],xF[9],xF_tmp[9],xC[9],xS[9],xP[9],xP_PullBack[9],xSigma[9],xSigma_PullBack[9],det_xF,det_xH,xPsi;
  double _idofs_x[12],_ih[9],imP[9],imSigma[9];
  if(K!=NULL) bzero(K,12*12*sizeof(double));
  if(Koff!=NULL) bzero(Koff,12*12*sizeof(double));
  //zero edges
  int ee = 0;	
  for(;ee<6;ee++) {
    if(NBEDGE_H1(order_x_edge[ee])==0) continue;
    int nb = 3*NBEDGE_H1(order_x_edge[ee]);
    if(K_edge[ee]!=NULL) bzero(K_edge[ee],12*nb*sizeof(double)); }
  //zero faces
  int ff = 0;	
  for(;ff<4;ff++) {
    if(NBFACE_H1(order_x_face[ff])==0) continue;
    int nb = 3*NBFACE_H1(order_x_face[ff]);
    if(K_face[ff]!=NULL) bzero(K_face[ff],nb*12*sizeof(double)); }
  //zero volume
  if(NBVOLUME_H1(order_x_volume)!=0) {
    int nb = 3*NBVOLUME_H1(order_x_volume);
    if(K_volume!=NULL) bzero(K_volume,nb*12*sizeof(double)); 
  }
  int gg = 0;
  for(;gg<G_DIM;gg++) {
    double H[9];
    ierr = GradientOfDeformation(diffN,dofs_X_node,H);  CHKERRQ(ierr);
    ierr = HierarhicalDeformationGradient(
	order_max_edge,order_max_face,order_max_volume,
	order_X_edge,order_X_face,order_X_volume,
	diffN,diffN_edge,diffN_face,diffN_volume,
	dofs_X_edge,dofs_X_face,dofs_X_volume,
	gg,H); CHKERRQ(ierr); 
    ierr = MakeComplexTensor(H,ZERO,xH);  CHKERRQ(ierr);
    cblas_zcopy(9,xH,1,inv_xH,1);
    ierr = DeterminantComplexGradient(xH,&det_xH); CHKERRQ(ierr);
    ierr = InvertComplexGradient(inv_xH); CHKERRQ(ierr);
    double h[9];
    ierr = GradientOfDeformation(diffN,dofs_x_node,h);  CHKERRQ(ierr);
    ierr = HierarhicalDeformationGradient(
	order_max_edge,order_max_face,order_max_volume,
	order_x_edge,order_x_face,order_x_volume,
	diffN,diffN_edge,diffN_face,diffN_volume,
	dofs_x_edge,dofs_x_face,dofs_x_volume,
	gg,h); CHKERRQ(ierr); 
    int dd = 0;
    for(;dd<12;dd++) {
      bzero(_idofs_x,sizeof(double)*12);
      _idofs_x[dd] = eps;
      ierr = GradientOfDeformation(diffN,_idofs_x,_ih);  CHKERRQ(ierr);
      ierr = MakeComplexTensor(h,_ih,xh);  CHKERRQ(ierr);
      COMP_STRESSES
      TakeIm(xP_PullBack,imP);
      cblas_dscal(9,1./eps,imP,1);
      TakeIm(xSigma_PullBack,imSigma);
      cblas_dscal(9,1./eps,imSigma,1);
      int node = 0;
      for(;node<4;node++) {
	if(K!=NULL) {
	  K[3*12*node + 0*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[0],1);
	  K[3*12*node + 1*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[3],1);
	  K[3*12*node + 2*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[6],1); }
	if(Koff!=NULL) {
	  Koff[3*12*node + 0*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[0],1);
	  Koff[3*12*node + 1*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[3],1);
	  Koff[3*12*node + 2*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[6],1); }}
      ee = 0;
      for(;ee<6;ee++) {
        int pp = 0;
        for(;pp<NBEDGE_H1(order_x_edge[ee]);pp++) {
	  double *diff = &((diffN_edge[ee])[gg*3*NBEDGE_H1(order_max_edge[ee])+3*pp]);
	  if(K_edge[ee]!=NULL) {
	    (K_edge[ee])[3*pp*12 + 0*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	    (K_edge[ee])[3*pp*12 + 1*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	    (K_edge[ee])[3*pp*12 + 2*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }}}
      ff = 0;
      for(;ff<4;ff++) {
        int pp = 0;
        for(;pp<NBFACE_H1(order_x_face[ff]);pp++) {
	  double *diff = &((diffN_face[ff])[gg*3*NBFACE_H1(order_max_face[ff])+3*pp]);
	  if(K_face[ff]!=NULL) {
	    (K_face[ff])[3*pp*12 + 0*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	    (K_face[ff])[3*pp*12 + 1*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	    (K_face[ff])[3*pp*12 + 2*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }}}
      int pp = 0;
      for(;pp<NBVOLUME_H1(order_x_volume);pp++) {
	double *diff = &((diffN_volume)[gg*3*NBVOLUME_H1(order_max_volume)+3*pp]);
	if(K_volume!=NULL) {
	  (K_volume)[3*pp*12 + 0*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
  	  (K_volume)[3*pp*12 + 1*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	  (K_volume)[3*pp*12 + 2*12 + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }}
  }}
  PetscFunctionReturn(0);
}
PetscErrorCode Tangent_hh_hierachical_edge(
  int *order_max_edge,int *order_max_face,int order_max_volume,
  int *order_X_edge,int *order_X_face,int order_X_volume,
  int *order_x_edge,int *order_x_face,int order_x_volume,
  double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X_node,double *dofs_X_edge[],double *dofs_X_face[],double *dofs_X_volume,
  double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  //temperature
  double thermal_expansion,double thermal_load_factor,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //
  double *K[6],double *Koff[6],
  double *K_edge[6][6],double *K_face[4][6],double *K_volume[6],
  int G_DIM,const double *G_W) {
  PetscFunctionBegin;
  double i_thermal_load_factor = 0;
  double ZERO[9];
  bzero(ZERO,sizeof(double)*9);
  __CLPK_doublecomplex xh[9],xH[9],inv_xH[9],xF[9],xF_tmp[9],xC[9],xS[9],xP[9],xP_PullBack[9],xSigma[9],xSigma_PullBack[9],det_xF,det_xH,xPsi;
  double _ih[9],imP[9],imSigma[9];
  int ee = 0;	
  for(;ee<6;ee++) {
    if(NBEDGE_H1(order_x_edge[ee])==0) continue;
    int nb = 3*NBEDGE_H1(order_x_edge[ee]);
    if(K!=NULL) if(K[ee]!=NULL) bzero(K[ee],12*nb*sizeof(double));
    if(Koff!=NULL) if(Koff[ee]!=NULL) bzero(Koff[ee],12*nb*sizeof(double));
    int EE = 0;
    for(;EE<6;EE++) {
      int nb2 = 3*NBEDGE_H1(order_x_edge[EE]);
      if(K_edge[EE][ee]!=NULL) bzero(K_edge[EE][ee],nb2*nb*sizeof(double)); }
    int FF = 0;
    for(;FF<4;FF++) {
      int nb2 = 3*NBFACE_H1(order_x_face[FF]);
      if(K_face[FF][ee]!=NULL) bzero(K_face[FF][ee],nb2*nb*sizeof(double)); 
    }
    int nb2 = 3*NBVOLUME_H1(order_x_volume);
    if(K_volume[ee]!=NULL) bzero(K_volume[ee],nb*nb2*sizeof(double)); 
  }
  int EE = 0;
  for(;EE<6;EE++) {
    int nb_edge_dofs = 3*NBEDGE_H1(order_x_edge[EE]);
    if(nb_edge_dofs == 0) continue;
    double _idofs_x[nb_edge_dofs];
    int gg = 0;
    for(;gg<G_DIM;gg++) {
      double H[9];
      ierr = GradientOfDeformation(diffN,dofs_X_node,H);  CHKERRQ(ierr);
      ierr = HierarhicalDeformationGradient(
	order_max_edge,order_max_face,order_max_volume,
	order_X_edge,order_X_face,order_X_volume,
	diffN,diffN_edge,diffN_face,diffN_volume,
	dofs_X_edge,dofs_X_face,dofs_X_volume,
	gg,H); CHKERRQ(ierr); 
      ierr = MakeComplexTensor(H,ZERO,xH);  CHKERRQ(ierr);
      cblas_zcopy(9,xH,1,inv_xH,1);
      ierr = DeterminantComplexGradient(xH,&det_xH); CHKERRQ(ierr);
      ierr = InvertComplexGradient(inv_xH); CHKERRQ(ierr);
      double h[9];
      ierr = GradientOfDeformation(diffN,dofs_x_node,h);  CHKERRQ(ierr);
      ierr = HierarhicalDeformationGradient(
	order_max_edge,order_max_face,order_max_volume,
	order_x_edge,order_x_face,order_x_volume,
	diffN,diffN_edge,diffN_face,diffN_volume,
	dofs_x_edge,dofs_x_face,dofs_x_volume,
	gg,h); CHKERRQ(ierr); 
      int dd = 0;
      for(;dd<nb_edge_dofs;dd++) {
	bzero(_idofs_x,sizeof(double)*nb_edge_dofs);
	_idofs_x[dd] = eps;   
	double *diff_edge = &(diffN_edge[EE])[gg*3*NBEDGE_H1(order_x_edge[EE])];
	H1_EdgeGradientOfDeformation_hierachical(order_x_edge[EE],diff_edge,_idofs_x,_ih); 
        ierr = MakeComplexTensor(h,_ih,xh);  CHKERRQ(ierr);
        COMP_STRESSES
        TakeIm(xP_PullBack,imP);
        cblas_dscal(9,1./eps,imP,1);
        TakeIm(xSigma_PullBack,imSigma);
        cblas_dscal(9,1./eps,imSigma,1);
        int node = 0;
        for(;node<4;node++) {
	  if(K!=NULL) if(K[EE]!=NULL) {
	    (K[EE])[3*node*nb_edge_dofs + 0*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[0],1);
	    (K[EE])[3*node*nb_edge_dofs + 1*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[3],1);
	    (K[EE])[3*node*nb_edge_dofs + 2*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[6],1); }
	  if(Koff!=NULL) if(Koff[EE]!=NULL) {
	    (Koff[EE])[3*node*nb_edge_dofs + 0*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[0],1);
	    (Koff[EE])[3*node*nb_edge_dofs + 1*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[3],1);
	    (Koff[EE])[3*node*nb_edge_dofs + 2*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[6],1); }}
        ee = 0;
        for(;ee<6;ee++) {
          int pp = 0;
          for(;pp<NBEDGE_H1(order_x_edge[ee]);pp++) {
	    double *diff = &((diffN_edge[ee])[gg*3*NBEDGE_H1(order_max_edge[ee])+3*pp]);
	    if(K_edge[ee][EE]!=NULL) {
	      (K_edge[ee][EE])[3*pp*nb_edge_dofs + 0*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	      (K_edge[ee][EE])[3*pp*nb_edge_dofs + 1*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1); 
	      (K_edge[ee][EE])[3*pp*nb_edge_dofs + 2*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); } }}
        int ff = 0;
        for(;ff<4;ff++) {
          int pp = 0;
          for(;pp<NBFACE_H1(order_x_face[ff]);pp++) {
	    double *diff = &((diffN_face[ff])[gg*3*NBFACE_H1(order_max_face[ff])+3*pp]);
	    if(K_face[ff][EE]!=NULL) {
	      (K_face[ff][EE])[3*pp*nb_edge_dofs + 0*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	      (K_face[ff][EE])[3*pp*nb_edge_dofs + 1*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1); 
	      (K_face[ff][EE])[3*pp*nb_edge_dofs + 2*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); } }}
	int pp = 0;
	for(;pp<NBVOLUME_H1(order_x_volume);pp++) {
	  double *diff = &((diffN_volume)[gg*3*NBVOLUME_H1(order_max_volume)+3*pp]);
	  if(K_volume!=NULL) {
	    (K_volume[EE])[3*pp*nb_edge_dofs + 0*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	    (K_volume[EE])[3*pp*nb_edge_dofs + 1*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	    (K_volume[EE])[3*pp*nb_edge_dofs + 2*nb_edge_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }}
    }}}
  PetscFunctionReturn(0);
}
PetscErrorCode Tangent_hh_hierachical_face(
  int *order_max_edge,int *order_max_face,int order_max_volume,
  int *order_X_edge,int *order_X_face,int order_X_volume,
  int *order_x_edge,int *order_x_face,int order_x_volume,
  double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X_node,double *dofs_X_edge[],double *dofs_X_face[],double *dofs_X_volume,
  double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  //temperature
  double thermal_expansion,double thermal_load_factor,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //
  double *K[4],double *Koff[4],
  double *K_edge[6][4],double *K_face[4][4],double *K_volume[4],
  int G_DIM,const double *G_W) {
  PetscFunctionBegin;
  double i_thermal_load_factor = 0;
  double ZERO[9];
  bzero(ZERO,sizeof(double)*9);
  __CLPK_doublecomplex xh[9],xH[9],inv_xH[9],xF[9],xF_tmp[9],xC[9],xS[9],xP[9],xP_PullBack[9],xSigma[9],xSigma_PullBack[9],det_xF,det_xH,xPsi;
  double _ih[9],imP[9],imSigma[9];
  int ff = 0;	
  for(;ff<4;ff++) {
    if(NBFACE_H1(order_x_face[ff])==0) continue;
    int nb = 3*NBFACE_H1(order_x_face[ff]);
    if(K[ff]!=NULL) bzero(K[ff],12*nb*sizeof(double));
    if(Koff!=NULL) {
      if(Koff[ff]!=NULL) bzero(Koff[ff],12*nb*sizeof(double));
    }
    int EE = 0;
    for(;EE<6;EE++) {
      int nb2 = 3*NBEDGE_H1(order_x_edge[EE]);
      if(K_edge[EE][ff]!=NULL) bzero(K_edge[EE][ff],nb2*nb*sizeof(double)); }
    int FF = 0;
    for(;FF<4;FF++) {
      int nb2 = 3*NBFACE_H1(order_x_face[FF]);
      if(K_face[FF][ff]!=NULL) bzero(K_face[FF][ff],nb2*nb*sizeof(double)); 
    }
    int nb2 = 3*NBVOLUME_H1(order_x_volume);
    if(K_volume[ff]!=NULL) bzero(K_volume[ff],nb*nb2*sizeof(double)); 
  }
  int FF = 0;
  for(;FF<4;FF++) {
    int nb_face_dofs = 3*NBFACE_H1(order_x_face[FF]);
    if(nb_face_dofs == 0) continue;
    double _idofs_x[nb_face_dofs];
    int dd = 0;
    for(;dd<nb_face_dofs;dd++) {
      bzero(_idofs_x,sizeof(double)*nb_face_dofs);
      _idofs_x[dd] = eps;   
      int gg = 0;
      for(;gg<G_DIM;gg++) {
	double H[9];
	ierr = GradientOfDeformation(diffN,dofs_X_node,H);  CHKERRQ(ierr);
	ierr = HierarhicalDeformationGradient(
	  order_max_edge,order_max_face,order_max_volume,
	  order_X_edge,order_X_face,order_X_volume,
	  diffN,diffN_edge,diffN_face,diffN_volume,
	  dofs_X_edge,dofs_X_face,dofs_X_volume,
	  gg,H); CHKERRQ(ierr); 
	ierr = MakeComplexTensor(H,ZERO,xH);  CHKERRQ(ierr);
	cblas_zcopy(9,xH,1,inv_xH,1);
	ierr = DeterminantComplexGradient(xH,&det_xH); CHKERRQ(ierr);
	ierr = InvertComplexGradient(inv_xH); CHKERRQ(ierr);
	double h[9];
	ierr = GradientOfDeformation(diffN,dofs_x_node,h);  CHKERRQ(ierr);
	ierr = HierarhicalDeformationGradient(
	  order_max_edge,order_max_face,order_max_volume,
	  order_x_edge,order_x_face,order_x_volume,
	  diffN,diffN_edge,diffN_face,diffN_volume,
	  dofs_x_edge,dofs_x_face,dofs_x_volume,
	  gg,h); CHKERRQ(ierr); 
	double *diff_face = &(diffN_face[FF])[gg*3*NBFACE_H1(order_x_face[FF])];
	H1_FaceGradientOfDeformation_hierachical(order_x_face[FF],diff_face,_idofs_x,_ih); 
        ierr = MakeComplexTensor(h,_ih,xh);  CHKERRQ(ierr);
        COMP_STRESSES
        TakeIm(xP_PullBack,imP);
        cblas_dscal(9,1./eps,imP,1);
        TakeIm(xSigma_PullBack,imSigma);
        cblas_dscal(9,1./eps,imSigma,1);
        int node = 0;
        for(;node<4;node++) {
	  if(K[FF]!=NULL) {
	    (K[FF])[3*node*nb_face_dofs + 0*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[0],1);
	    (K[FF])[3*node*nb_face_dofs + 1*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[3],1);
	    (K[FF])[3*node*nb_face_dofs + 2*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[6],1); }
	  if(Koff!=NULL) {
	  if(Koff[FF]!=NULL) {
	    (Koff[FF])[3*node*nb_face_dofs + 0*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[0],1);
	    (Koff[FF])[3*node*nb_face_dofs + 1*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[3],1);
	    (Koff[FF])[3*node*nb_face_dofs + 2*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[6],1); }}
	}
        int ee = 0;
        for(;ee<6;ee++) {
          int pp = 0;
          for(;pp<NBEDGE_H1(order_x_edge[ee]);pp++) {
	    double *diff = &((diffN_edge[ee])[gg*3*NBEDGE_H1(order_max_edge[ee])+3*pp+0]);
	    if(K_edge[ee][FF]!=NULL) {
	      (K_edge[ee][FF])[3*pp*nb_face_dofs + 0*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	      (K_edge[ee][FF])[3*pp*nb_face_dofs + 1*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1); 
	      (K_edge[ee][FF])[3*pp*nb_face_dofs + 2*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); } }}
        ff = 0;
        for(;ff<4;ff++) {
          int pp = 0;
          for(;pp<NBFACE_H1(order_x_face[ff]);pp++) {
	    double *diff = &((diffN_face[ff])[gg*3*NBFACE_H1(order_max_face[ff])+3*pp+0]);
	    if(K_face[ff][FF]!=NULL) {
	      (K_face[ff][FF])[3*pp*nb_face_dofs + 0*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	      (K_face[ff][FF])[3*pp*nb_face_dofs + 1*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1); 
	      (K_face[ff][FF])[3*pp*nb_face_dofs + 2*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); } }}
	int pp = 0;
	for(;pp<NBVOLUME_H1(order_x_volume);pp++) {
	  double *diff = &((diffN_volume)[gg*3*NBVOLUME_H1(order_max_volume)+3*pp]);
	  if(K_volume!=NULL) {
	    (K_volume[FF])[3*pp*nb_face_dofs + 0*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	    (K_volume[FF])[3*pp*nb_face_dofs + 1*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	    (K_volume[FF])[3*pp*nb_face_dofs + 2*nb_face_dofs + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }}
    }}}
  PetscFunctionReturn(0);
}
PetscErrorCode Tangent_hh_hierachical_volume(
  int *order_max_edge,int *order_max_face,int order_max_volume,
  int *order_X_edge,int *order_X_face,int order_X_volume,
  int *order_x_edge,int *order_x_face,int order_x_volume,
  double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X_node,double *dofs_X_edge[],double *dofs_X_face[],double *dofs_X_volume,
  double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  //temperature
  double thermal_expansion,double thermal_load_factor,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //
  double *K,double *Koff,double *K_edge[6],double *K_face[4],double *K_volume,int G_DIM,const double *G_W) {
  PetscFunctionBegin;
  double i_thermal_load_factor = 0;
  double ZERO[9];
  bzero(ZERO,sizeof(double)*9);
  __CLPK_doublecomplex xh[9],xH[9],inv_xH[9],xF[9],xF_tmp[9],xC[9],xS[9],xP[9],xP_PullBack[9],xSigma[9],xSigma_PullBack[9],det_xF,det_xH,xPsi;
  double _ih[9],imP[9],imSigma[9];
  int nb_dofs_volume = 3*NBVOLUME_H1(order_x_volume);
  if(K!=NULL) bzero(K,nb_dofs_volume*12*sizeof(double));
  if(Koff!=NULL) bzero(Koff,nb_dofs_volume*12*sizeof(double));
  //zero edges
  int ee = 0;	
  for(;ee<6;ee++) {
    if(NBEDGE_H1(order_x_edge[ee])==0) continue;
    int nb2 = 3*NBEDGE_H1(order_x_edge[ee]);
    if(K_edge[ee]!=NULL) bzero(K_edge[ee],nb_dofs_volume*nb2*sizeof(double)); 
  }
  //zero faces
  int ff = 0;	
  for(;ff<4;ff++) {
    if(NBFACE_H1(order_x_face[ff])==0) continue;
    int nb2 = 3*NBFACE_H1(order_x_face[ff]);
    if(K_face[ff]!=NULL) bzero(K_face[ff],nb_dofs_volume*nb2*sizeof(double)); 
  }
  //zero volume
  if(NBVOLUME_H1(order_x_volume)!=0) {
    if(K_volume!=NULL) bzero(K_volume,nb_dofs_volume*nb_dofs_volume*sizeof(double)); 
  }
  double _idofs_x[nb_dofs_volume];
  int gg = 0;
  for(;gg<G_DIM;gg++) {
    double H[9];
    ierr = GradientOfDeformation(diffN,dofs_X_node,H);  CHKERRQ(ierr);
    ierr = HierarhicalDeformationGradient(
	order_max_edge,order_max_face,order_max_volume,
	order_X_edge,order_X_face,order_X_volume,
	diffN,diffN_edge,diffN_face,diffN_volume,
	dofs_X_edge,dofs_X_face,dofs_X_volume,
	gg,H); CHKERRQ(ierr); 
    ierr = MakeComplexTensor(H,ZERO,xH);  CHKERRQ(ierr);
    cblas_zcopy(9,xH,1,inv_xH,1);
    ierr = DeterminantComplexGradient(xH,&det_xH); CHKERRQ(ierr);
    ierr = InvertComplexGradient(inv_xH); CHKERRQ(ierr);
    double h[9];
    ierr = GradientOfDeformation(diffN,dofs_x_node,h);  CHKERRQ(ierr);
    ierr = HierarhicalDeformationGradient(
	order_max_edge,order_max_face,order_max_volume,
	order_x_edge,order_x_face,order_x_volume,
	diffN,diffN_edge,diffN_face,diffN_volume,
	dofs_x_edge,dofs_x_face,dofs_x_volume,
	gg,h); CHKERRQ(ierr); 
    int dd = 0;
    for(;dd<nb_dofs_volume;dd++) {
      bzero(_idofs_x,sizeof(double)*nb_dofs_volume);
      _idofs_x[dd] = eps;
      double *diff_volume = &((diffN_volume)[gg*3*NBVOLUME_H1(order_x_volume)]);
      H1_VolumeGradientOfDeformation_hierachical(order_x_volume,diff_volume,_idofs_x,_ih); 
      ierr = MakeComplexTensor(h,_ih,xh);  CHKERRQ(ierr);
      COMP_STRESSES
      TakeIm(xP_PullBack,imP);
      cblas_dscal(9,1./eps,imP,1);
      TakeIm(xSigma_PullBack,imSigma);
      cblas_dscal(9,1./eps,imSigma,1);
      int node = 0;
      for(;node<4;node++) {
	if(K!=NULL) {
	  K[3*node*nb_dofs_volume + 0*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[0],1);
	  K[3*node*nb_dofs_volume + 1*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[3],1);
	  K[3*node*nb_dofs_volume + 2*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imP[6],1); }
	if(Koff!=NULL) {
	  Koff[3*node*nb_dofs_volume + 0*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[0],1);
	  Koff[3*node*nb_dofs_volume + 1*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[3],1);
	  Koff[3*node*nb_dofs_volume + 2*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,&diffN[node*3+0],1,&imSigma[6],1); }}
      ee = 0;
      for(;ee<6;ee++) {
        int pp = 0;
        for(;pp<NBEDGE_H1(order_x_edge[ee]);pp++) {
	  double *diff = &((diffN_edge[ee])[gg*3*NBEDGE_H1(order_max_edge[ee])+3*pp]);
	  if(K_edge[ee]!=NULL) {
	    (K_edge[ee])[3*pp*nb_dofs_volume + 0*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	    (K_edge[ee])[3*pp*nb_dofs_volume + 1*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	    (K_edge[ee])[3*pp*nb_dofs_volume + 2*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }}}
      ff = 0;
      for(;ff<4;ff++) {
        int pp = 0;
        for(;pp<NBFACE_H1(order_x_face[ff]);pp++) {
	  double *diff = &((diffN_face[ff])[gg*3*NBFACE_H1(order_max_face[ff])+3*pp]);
	  if(K_face[ff]!=NULL) {
	    (K_face[ff])[3*pp*nb_dofs_volume + 0*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
	    (K_face[ff])[3*pp*nb_dofs_volume + 1*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	    (K_face[ff])[3*pp*nb_dofs_volume + 2*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }}}
      int pp = 0;
      for(;pp<NBVOLUME_H1(order_x_volume);pp++) {
	double *diff = &((diffN_volume)[gg*3*NBVOLUME_H1(order_max_volume)+3*pp]);
	if(K_volume!=NULL) {
	  (K_volume)[3*pp*nb_dofs_volume + 0*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[0],1);
  	  (K_volume)[3*pp*nb_dofs_volume + 1*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[3],1);
	  (K_volume)[3*pp*nb_dofs_volume + 2*nb_dofs_volume + dd] += alpha*G_W[gg]*cblas_ddot(3,diff,1,&imP[6],1); }}
  }}
  PetscFunctionReturn(0);
}
//External forces
PetscErrorCode Traction_hierarchical(int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *t,double *t_edge[],double *t_face,
  double *traction,int gg) {
  PetscFunctionBegin;
  int dd,ee;
  for(dd = 0;dd<3;dd++) traction[dd] = cblas_ddot(3,&N[gg*3],1,&t[dd],3);
  if(t_face!=NULL) {
    int nb_dofs_face = NBFACE_H1(order);
    if(nb_dofs_face>0) {
      for(dd = 0;dd<3;dd++) traction[dd] += cblas_ddot(nb_dofs_face,&N_face[gg*nb_dofs_face],1,&t_face[dd],3);
    }
  }
  if(t_edge!=NULL) {
    ee = 0;
    for(;ee<3;ee++) {
      if(t_edge[ee] == NULL) continue;
      int nb_dofs_edge = NBEDGE_H1(order_edge[ee]);
      if(nb_dofs_edge>0) {
	for(dd = 0;dd<3;dd++) {
	  traction[dd] += cblas_ddot(nb_dofs_edge,&(N_edge[ee][gg*nb_dofs_edge]),1,&(t_edge[ee][dd]),3);
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Fext_h_hierarchical(int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *t,double *t_edge[],double *t_face,
  double *dofs_x,double *dofs_x_edge[],double *dofs_x_face,
  double *idofs_x,double *idofs_x_edge[],double *idofs_x_face,
  double *Fext,double *Fext_edge[],double *Fext_face,
  double *iFext,double *iFext_edge[],double *iFext_face,
  int g_dim,const double *g_w) {
  PetscFunctionBegin;
  int dd,nn,ee,gg;
  if(Fext!=NULL) bzero(Fext,9*sizeof(double));
  if(iFext!=NULL) bzero(iFext,9*sizeof(double));
  ee = 0;
  for(;ee<3;ee++) {
    int nb_dofs_edge = NBEDGE_H1(order_edge[ee]);
    if(nb_dofs_edge==0) continue;
    if(Fext_edge!=NULL) bzero(Fext_edge[ee],3*nb_dofs_edge*sizeof(double));
    if(iFext_edge!=NULL) bzero(iFext_edge[ee],3*nb_dofs_edge*sizeof(double));
  } 
  int nb_dofs_face = NBFACE_H1(order);
  if(nb_dofs_face!=0) {
    if(Fext_face!=NULL) bzero(Fext_face,3*nb_dofs_face*sizeof(double));
    if(iFext_face!=NULL) bzero(iFext_face,3*nb_dofs_face*sizeof(double));
  }
  gg = 0;
  for(;gg<g_dim;gg++) {
    double traction[3] = {0,0,0};
    ierr = Traction_hierarchical(order,order_edge,N,N_face,N_edge,t,t_edge,t_face,traction,gg); CHKERRQ(ierr);
    __CLPK_doublecomplex xnormal[3],xs1[3],xs2[3];
    ierr = Normal_hierarchical(
      //FIXME: order of triacions approximation could be diffrernt for field approx
      order,order_edge, 
      order,order_edge,
      diffN,diffN_face,diffN_edge,
      dofs_x,dofs_x_edge,dofs_x_face,idofs_x,
      idofs_x_edge,idofs_x_face,
      xnormal,xs1,xs2,gg); CHKERRQ(ierr);
    ierr = Base_scale(xnormal,xs1,xs2); CHKERRQ(ierr);
    double normal_real[3],s1_real[3],s2_real[3];
    double normal_imag[3],s1_imag[3],s2_imag[3];
    for(dd = 0;dd<3;dd++) {
      normal_real[dd] = 0.5*xnormal[dd].r;
      normal_imag[dd] = 0.5*xnormal[dd].i;
      s1_real[dd] = 0.5*xs1[dd].r;
      s1_imag[dd] = 0.5*xs1[dd].i;
      s2_real[dd] = 0.5*xs2[dd].r;
      s2_imag[dd] = 0.5*xs2[dd].i;
    }
    nn = 0;
    for(;nn<3;nn++) {
      if(Fext!=NULL) 
	for(dd = 0;dd<3;dd++) {
	  /*fprintf(stderr,"%d %f %f %f %f %f %f\n",
	    gg,g_w[gg],N[3*gg+nn],normal_real[dd],
	    traction[0],traction[1],traction[2]);*/
	  Fext[3*nn+dd] += g_w[gg]*N[3*gg+nn]*normal_real[dd]*traction[2];
	  Fext[3*nn+dd] += g_w[gg]*N[3*gg+nn]*s1_real[dd]*traction[0];
	  Fext[3*nn+dd] += g_w[gg]*N[3*gg+nn]*s2_real[dd]*traction[1];
	}
      if(iFext!=NULL) 
	for(dd = 0;dd<3;dd++) {
	  iFext[3*nn+dd] += g_w[gg]*N[3*gg+nn]*normal_imag[dd]*traction[2];
	  iFext[3*nn+dd] += g_w[gg]*N[3*gg+nn]*s1_imag[dd]*traction[0];
	  iFext[3*nn+dd] += g_w[gg]*N[3*gg+nn]*s2_imag[dd]*traction[1];
	}
    }
    ee = 0;
    for(;ee<3;ee++) {
      int nb_dofs_edge = NBEDGE_H1(order_edge[ee]);
      if(nb_dofs_edge == 0) continue;
      int nn = 0;
      for(;nn<nb_dofs_edge;nn++) {
	if(Fext_edge!=NULL) 
	  for(dd = 0;dd<3;dd++) {
	    Fext_edge[ee][3*nn+dd] += g_w[gg]*N_edge[ee][gg*nb_dofs_edge+nn]*normal_real[dd]*traction[2];
	    Fext_edge[ee][3*nn+dd] += g_w[gg]*N_edge[ee][gg*nb_dofs_edge+nn]*s1_real[dd]*traction[0];
	    Fext_edge[ee][3*nn+dd] += g_w[gg]*N_edge[ee][gg*nb_dofs_edge+nn]*s2_real[dd]*traction[1];
	  }
	if(iFext_edge!=NULL) {
	  for(dd = 0;dd<3;dd++) {
	    iFext_edge[ee][3*nn+dd] += g_w[gg]*N_edge[ee][gg*nb_dofs_edge+nn]*normal_imag[dd]*traction[2];
	    iFext_edge[ee][3*nn+dd] += g_w[gg]*N_edge[ee][gg*nb_dofs_edge+nn]*s1_imag[dd]*traction[0];
	    iFext_edge[ee][3*nn+dd] += g_w[gg]*N_edge[ee][gg*nb_dofs_edge+nn]*s2_imag[dd]*traction[1];
	  }
	}
      }
    }
    if(nb_dofs_face!=0) {
      nn = 0;
      for(;nn<nb_dofs_face;nn++) {
	if(Fext_face!=NULL) 
	  for(dd = 0;dd<3;dd++) {
	    Fext_face[3*nn+dd] += g_w[gg]*N_face[gg*nb_dofs_face+nn]*normal_real[dd]*traction[2];
	    Fext_face[3*nn+dd] += g_w[gg]*N_face[gg*nb_dofs_face+nn]*s1_real[dd]*traction[0];
	    Fext_face[3*nn+dd] += g_w[gg]*N_face[gg*nb_dofs_face+nn]*s2_real[dd]*traction[1];
	  }
	if(iFext_face!=NULL) 
	  for(dd = 0;dd<3;dd++) {
	    iFext_face[3*nn+dd] += g_w[gg]*N_face[gg*nb_dofs_face+nn]*normal_imag[dd]*traction[2];
	    iFext_face[3*nn+dd] += g_w[gg]*N_face[gg*nb_dofs_face+nn]*s1_imag[dd]*traction[0];
	    iFext_face[3*nn+dd] += g_w[gg]*N_face[gg*nb_dofs_face+nn]*s2_imag[dd]*traction[1];
	  }
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode KExt_hh_hierarchical(double eps,int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *t,double *t_edge[],double *t_face,
  double *dofs_x,double *dofs_x_edge[],double *dofs_x_face,
  double *KExt_hh,double* KExt_edgeh[],double *KExt_faceh,
  int g_dim,const double *g_w) {
  PetscFunctionBegin;
  int gg,dd,ii,nn,ee;
  bzero(KExt_hh,9*9*sizeof(double));
  if(KExt_edgeh!=NULL) {
    for(ee = 0;ee<3;ee++) {
      int nb_dofs_edge = NBEDGE_H1(order_edge[ee]);
      bzero(KExt_edgeh[ee],9*3*nb_dofs_edge*sizeof(double));
    }
  }
  int nb_dofs_face = NBFACE_H1(order);
  if(KExt_faceh!=NULL) {
    bzero(KExt_faceh,9*3*nb_dofs_face*sizeof(double));
  }
  for(gg = 0;gg<g_dim;gg++) {
    double traction[3] = {0,0,0};
    ierr = Traction_hierarchical(order,order_edge,N,N_face,N_edge,t,t_edge,t_face,traction,gg); CHKERRQ(ierr);
    //
    __CLPK_doublecomplex xnormal[3],xs1[3],xs2[3];
    double idofs_x[9];
    for(ii = 0;ii<9;ii++) {
      bzero(idofs_x,9*sizeof(double));
      idofs_x[ii] = eps;
      ierr = Normal_hierarchical(
	order,order_edge, //FIXME
	order,order_edge,
	diffN,diffN_face,diffN_edge,
	dofs_x,dofs_x_edge,dofs_x_face,
	idofs_x,NULL,NULL,
	xnormal,xs1,xs2,gg); CHKERRQ(ierr);
      ierr = Base_scale(xnormal,xs1,xs2); CHKERRQ(ierr);
      double normal_imag[3],s1_imag[3],s2_imag[3];
      for(dd = 0;dd<3;dd++) {
	normal_imag[dd] = 0.5*xnormal[dd].i/eps;
	s1_imag[dd] = 0.5*xs1[dd].i/eps;
	s2_imag[dd] = 0.5*xs2[dd].i/eps;
      }
      nn = 0;
      for(;nn<3;nn++) {
	for(dd = 0;dd<3;dd++) {
	  KExt_hh[ii+9*3*nn+9*dd] += g_w[gg]*N[3*gg+nn]*normal_imag[dd]*traction[2];
	  KExt_hh[ii+9*3*nn+9*dd] += g_w[gg]*N[3*gg+nn]*s1_imag[dd]*traction[0];
	  KExt_hh[ii+9*3*nn+9*dd] += g_w[gg]*N[3*gg+nn]*s2_imag[dd]*traction[1];
	}
      }
      if(KExt_edgeh!=NULL) {
	for(ee = 0;ee<3;ee++) {
	  int nb_dofs_edge = NBEDGE_H1(order_edge[ee]);
	  for(nn = 0;nn<nb_dofs_edge;nn++) {
	    for(dd = 0;dd<3;dd++) {
	      KExt_edgeh[ee][ii+9*3*nn + 9*dd] += g_w[gg]*N_edge[ee][nb_dofs_edge*gg+nn]*normal_imag[dd]*traction[2];
	      KExt_edgeh[ee][ii+9*3*nn + 9*dd] += g_w[gg]*N_edge[ee][nb_dofs_edge*gg+nn]*s1_imag[dd]*traction[0];
	      KExt_edgeh[ee][ii+9*3*nn + 9*dd] += g_w[gg]*N_edge[ee][nb_dofs_edge*gg+nn]*s2_imag[dd]*traction[1];
	    }
	  }
	}
      }
      if(KExt_faceh!=NULL) {
	for(nn = 0;nn<nb_dofs_face;nn++) {
	  for(dd = 0;dd<3;dd++) {
	    KExt_faceh[ii+3*9*nn+9*dd] += g_w[gg]*N_face[nb_dofs_face*gg+nn]*normal_imag[dd]*traction[2];
	    KExt_faceh[ii+3*9*nn+9*dd] += g_w[gg]*N_face[nb_dofs_face*gg+nn]*s1_imag[dd]*traction[0];
	    KExt_faceh[ii+3*9*nn+9*dd] += g_w[gg]*N_face[nb_dofs_face*gg+nn]*s2_imag[dd]*traction[1];
	  }
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode KExt_hh_hierarchical_edge(double eps,int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *t,double *t_edge[],double *t_face,
  double *dofs_x,double *dofs_x_edge[],double *dofs_x_face,
  double* KExt_hedge[3],double* KExt_edgeedge[3][3],double *KExt_faceedge[3],
  int g_dim,const double *g_w) {
  PetscFunctionBegin;
  int gg,dd,ii,nn,ee,EE;
  int nb_dofs_face = NBFACE_H1(order);
  for(EE=0;EE<3;EE++) {
    int nb_dofs_edge_EE = NBEDGE_H1(order_edge[EE]);
    bzero(KExt_hedge[EE],9*3*nb_dofs_edge_EE*sizeof(double));
    if(KExt_edgeedge!=NULL) {
      for(ee = 0;ee<3;ee++) {
	int nb_dofs_edge = NBEDGE_H1(order_edge[ee]);
	bzero(KExt_edgeedge[EE][ee],3*nb_dofs_edge_EE*3*nb_dofs_edge*sizeof(double));
      }
    }
    if(KExt_faceedge!=NULL) {
      bzero(KExt_faceedge[EE],3*nb_dofs_edge_EE*3*nb_dofs_face*sizeof(double));
    }
  }
  for(gg = 0;gg<g_dim;gg++) {
    double traction[3]  = {0,0,0};
    ierr = Traction_hierarchical(order,order_edge,N,N_face,N_edge,t,t_edge,t_face,traction,gg); CHKERRQ(ierr);
    for(EE = 0;EE<3;EE++) {
      int nb_dofs_edge_EE = NBEDGE_H1(order_edge[EE]);
      double* idofs_x_edge[3] = { NULL, NULL, NULL };
      double idofs_x_edge_EE[3*nb_dofs_edge_EE];
      idofs_x_edge[EE] = idofs_x_edge_EE;
      for(ii = 0;ii<3*nb_dofs_edge_EE;ii++) {
        bzero(idofs_x_edge_EE,3*nb_dofs_edge_EE*sizeof(double));
        idofs_x_edge_EE[ii] = eps;
	__CLPK_doublecomplex xnormal[3],xs1[3],xs2[3];
	ierr = Normal_hierarchical(
	  order,order_edge, //FIXME
	  order,order_edge,
	  diffN,diffN_face,diffN_edge,
	  dofs_x,dofs_x_edge,dofs_x_face,NULL,idofs_x_edge,NULL,
	  xnormal,xs1,xs2,gg); CHKERRQ(ierr);
	ierr = Base_scale(xnormal,xs1,xs2); CHKERRQ(ierr);
	double normal_imag[3],s1_imag[3],s2_imag[3];
	for(dd = 0;dd<3;dd++) {
	  normal_imag[dd] = 0.5*xnormal[dd].i/eps;
	  s1_imag[dd] = 0.5*xs1[dd].i/eps;
	  s2_imag[dd] = 0.5*xs2[dd].i/eps;
	}
        for(nn = 0;nn<3;nn++) {
	  for(dd = 0;dd<3;dd++) {
	    KExt_hedge[EE][ii + 3*nb_dofs_edge_EE*3*nn + 3*nb_dofs_edge_EE*dd] += g_w[gg]*N[3*gg+nn]*normal_imag[dd]*traction[2];
	    KExt_hedge[EE][ii + 3*nb_dofs_edge_EE*3*nn + 3*nb_dofs_edge_EE*dd] += g_w[gg]*N[3*gg+nn]*s1_imag[dd]*traction[0];
	    KExt_hedge[EE][ii + 3*nb_dofs_edge_EE*3*nn + 3*nb_dofs_edge_EE*dd] += g_w[gg]*N[3*gg+nn]*s2_imag[dd]*traction[1];
	  }
        }
        if(KExt_edgeedge!=NULL) {
	 for(ee = 0;ee<3;ee++) {
	  int nb_dofs_edge = NBEDGE_H1(order_edge[ee]);
  	  for(nn = 0;nn<nb_dofs_edge;nn++) {
  	    for(dd = 0;dd<3;dd++) {
	      KExt_edgeedge[EE][ee][ii + 3*nb_dofs_edge_EE*3*nn + 3*nb_dofs_edge_EE*dd] += g_w[gg]*N_edge[ee][nb_dofs_edge*gg+nn]*normal_imag[dd]*traction[2];
	      KExt_edgeedge[EE][ee][ii + 3*nb_dofs_edge_EE*3*nn + 3*nb_dofs_edge_EE*dd] += g_w[gg]*N_edge[ee][nb_dofs_edge*gg+nn]*s1_imag[dd]*traction[0];
	      KExt_edgeedge[EE][ee][ii + 3*nb_dofs_edge_EE*3*nn + 3*nb_dofs_edge_EE*dd] += g_w[gg]*N_edge[ee][nb_dofs_edge*gg+nn]*s2_imag[dd]*traction[1];
	    }
  	  }
	 }
        }
        if(KExt_faceedge!=NULL) {
	  for(nn = 0;nn<nb_dofs_face;nn++) {
	    for(dd = 0;dd<3;dd++) {
	      KExt_faceedge[EE][ii + 3*nb_dofs_edge_EE*3*nn + 3*nb_dofs_edge_EE*dd] += g_w[gg]*N_face[nb_dofs_face*gg+nn]*normal_imag[dd]*traction[2];
	      KExt_faceedge[EE][ii + 3*nb_dofs_edge_EE*3*nn + 3*nb_dofs_edge_EE*dd] += g_w[gg]*N_face[nb_dofs_face*gg+nn]*s1_imag[dd]*traction[0];
	      KExt_faceedge[EE][ii + 3*nb_dofs_edge_EE*3*nn + 3*nb_dofs_edge_EE*dd] += g_w[gg]*N_face[nb_dofs_face*gg+nn]*s2_imag[dd]*traction[1];
	    }
	  }
        }
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode KExt_hh_hierarchical_face(double eps,int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *t,double *t_edge[],double *t_face,
  double *dofs_x,double *dofs_x_edge[],double *dofs_x_face,
  double *KExt_hface,double *KExt_edgeface[3],double *KExt_faceface,
  int g_dim,const double *g_w) {
  PetscFunctionBegin;
  int gg,dd,ii,nn,ee;
  int nb_dofs_face = NBFACE_H1(order);
  bzero(KExt_hface,9*3*nb_dofs_face*sizeof(double));
  if(KExt_edgeface!=NULL) {	
    for(ee = 0;ee<3;ee++) {
      int nb_dofs_edge = NBEDGE_H1(order_edge[ee]);
      bzero(KExt_edgeface[ee],3*nb_dofs_face*3*nb_dofs_edge*sizeof(double));
    }
  }
  if(KExt_faceface!=NULL) {
    bzero(KExt_faceface,3*nb_dofs_face*3*nb_dofs_face*sizeof(double));
  }
  for(gg = 0;gg<g_dim;gg++) {
    double traction[3] = {0,0,0}; 
    ierr = Traction_hierarchical(order,order_edge,N,N_face,N_edge,t,t_edge,t_face,traction,gg); CHKERRQ(ierr);
    double idofs_x_face[3*nb_dofs_face];
    for(ii = 0;ii<3*nb_dofs_face;ii++) {
      bzero(idofs_x_face,3*nb_dofs_face*sizeof(double));
      idofs_x_face[ii] = eps;
      __CLPK_doublecomplex xnormal[3],xs1[3],xs2[3];
      ierr = Normal_hierarchical(
	order,order_edge, //FIXME
	order,order_edge,
	diffN,diffN_face,diffN_edge,
	dofs_x,dofs_x_edge,dofs_x_face,NULL,NULL,idofs_x_face,
	xnormal,xs1,xs2,gg); CHKERRQ(ierr);
      ierr = Base_scale(xnormal,xs1,xs2); CHKERRQ(ierr);
      double normal_imag[3],s1_imag[3],s2_imag[3];
      for(dd = 0;dd<3;dd++) {
	normal_imag[dd] = 0.5*xnormal[dd].i/eps;
	s1_imag[dd] = 0.5*xs1[dd].i/eps;
	s2_imag[dd] = 0.5*xs2[dd].i/eps;
      }
      for(nn = 0;nn<3;nn++) {
	for(dd = 0;dd<3;dd++) {
	  KExt_hface[ii + 3*nb_dofs_face*3*nn + 3*nb_dofs_face*dd] += g_w[gg]*N[3*gg+nn]*normal_imag[dd]*traction[2];
	  KExt_hface[ii + 3*nb_dofs_face*3*nn + 3*nb_dofs_face*dd] += g_w[gg]*N[3*gg+nn]*s1_imag[dd]*traction[0];
	  KExt_hface[ii + 3*nb_dofs_face*3*nn + 3*nb_dofs_face*dd] += g_w[gg]*N[3*gg+nn]*s2_imag[dd]*traction[1];
	}
      }
      if(KExt_edgeface!=NULL) {
	for(ee = 0;ee<3;ee++) {
	  int nb_dofs_edge = NBEDGE_H1(order_edge[ee]);
	  for(nn = 0;nn<nb_dofs_edge;nn++) {
	    for(dd = 0;dd<3;dd++) {
	      KExt_edgeface[ee][ii + 3*nb_dofs_face*3*nn + 3*nb_dofs_face*dd] += g_w[gg]*N_edge[ee][nb_dofs_edge*gg+nn]*normal_imag[dd]*traction[2];
	      KExt_edgeface[ee][ii + 3*nb_dofs_face*3*nn + 3*nb_dofs_face*dd] += g_w[gg]*N_edge[ee][nb_dofs_edge*gg+nn]*s1_imag[dd]*traction[0];
	      KExt_edgeface[ee][ii + 3*nb_dofs_face*3*nn + 3*nb_dofs_face*dd] += g_w[gg]*N_edge[ee][nb_dofs_edge*gg+nn]*s2_imag[dd]*traction[1];
	    }
	  }
	}
      }
      if(KExt_faceface!=NULL) {
	for(nn = 0;nn<nb_dofs_face;nn++) {
	  for(dd = 0;dd<3;dd++) {
	    KExt_faceface[ii + 3*nb_dofs_face*3*nn + 3*nb_dofs_face*dd] += g_w[gg]*N_face[nb_dofs_face*gg+nn]*normal_imag[dd]*traction[2];
	    KExt_faceface[ii + 3*nb_dofs_face*3*nn + 3*nb_dofs_face*dd] += g_w[gg]*N_face[nb_dofs_face*gg+nn]*s1_imag[dd]*traction[0];
	    KExt_faceface[ii + 3*nb_dofs_face*3*nn + 3*nb_dofs_face*dd] += g_w[gg]*N_face[nb_dofs_face*gg+nn]*s2_imag[dd]*traction[1];
	  }
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
//Material
PetscErrorCode Fext_H_hierarchical(
  int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *t,double *t_edge[],double *t_face,
  double *dofs_X,double *dofs_X_edge[],double *dofs_X_face,
  double *idofs_X,
  double *Fext,double *iFext,int g_dim,const double *g_w) {
  PetscFunctionBegin;
  int dd,nn,gg;
  if(Fext!=NULL) bzero(Fext,9*sizeof(double));
  if(iFext!=NULL) bzero(iFext,9*sizeof(double));
  gg = 0;
  for(;gg<g_dim;gg++) {
    double traction[3] = {0,0,0};
    ierr = Traction_hierarchical(order,order_edge,N,N_face,N_edge,t,t_edge,t_face,traction,gg); CHKERRQ(ierr);
    __CLPK_doublecomplex xnormal[3],xs1[3],xs2[3];
    ierr = Normal_hierarchical(
      order,order_edge, //FIXME
      order,order_edge,
      diffN,diffN_face,diffN_edge,
      dofs_X,dofs_X_edge,dofs_X_face,
      idofs_X,NULL,NULL,
      xnormal,xs1,xs2,gg); CHKERRQ(ierr);
    ierr = Base_scale(xnormal,xs1,xs2); CHKERRQ(ierr);
    double normal_real[3],s1_real[3],s2_real[3];
    double normal_imag[3],s1_imag[3],s2_imag[3];
    for(dd = 0;dd<3;dd++) {
      normal_real[dd] = 0.5*xnormal[dd].r;
      normal_imag[dd] = 0.5*xnormal[dd].i;
      s1_real[dd] = 0.5*xs1[dd].r;
      s1_imag[dd] = 0.5*xs1[dd].i;
      s2_real[dd] = 0.5*xs2[dd].r;
      s2_imag[dd] = 0.5*xs2[dd].i;
    }
    nn = 0;
    for(;nn<3;nn++) {
      if(Fext!=NULL) 
	for(dd = 0;dd<3;dd++) {
	  Fext[3*nn+dd] += g_w[gg]*N[3*gg+nn]*normal_real[dd]*traction[2];
	  Fext[3*nn+dd] += g_w[gg]*N[3*gg+nn]*s1_real[dd]*traction[0];
	  Fext[3*nn+dd] += g_w[gg]*N[3*gg+nn]*s2_real[dd]*traction[1];
	}
      if(iFext!=NULL) 
	for(dd = 0;dd<3;dd++) {
	  iFext[3*nn+dd] += g_w[gg]*N[3*gg+nn]*normal_imag[dd]*traction[2];
	  iFext[3*nn+dd] += g_w[gg]*N[3*gg+nn]*s1_imag[dd]*traction[0];
	  iFext[3*nn+dd] += g_w[gg]*N[3*gg+nn]*s2_imag[dd]*traction[1];
	}
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode KExt_HH_hierarchical(
  double eps,int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *t,double *t_edge[],double *t_face,
  double *dofs_X,double *dofs_X_edge[],double *dofs_X_face,
  double *KExt_HH,int g_dim,const double *g_w) {
  PetscFunctionBegin;
  int gg,dd,ii,nn;
  bzero(KExt_HH,9*9*sizeof(double));
  for(gg = 0;gg<g_dim;gg++) {
    double traction[3] = {0,0,0};
    ierr = Traction_hierarchical(order,order_edge,N,N_face,N_edge,t,t_edge,t_face,traction,gg); CHKERRQ(ierr);
    //
    __CLPK_doublecomplex xnormal[3],xs1[3],xs2[3];
    double idofs_X[9];
    for(ii = 0;ii<9;ii++) {
      bzero(idofs_X,9*sizeof(double));
      idofs_X[ii] = eps;
      ierr = Normal_hierarchical(
	order,order_edge, //FIXME
	order,order_edge,
	diffN,diffN_face,diffN_edge,
	dofs_X,dofs_X_edge,dofs_X_face,
	idofs_X,NULL,NULL,
	xnormal,xs1,xs2,gg); CHKERRQ(ierr);
      ierr = Base_scale(xnormal,xs1,xs2); CHKERRQ(ierr);
      double normal_imag[3],s1_imag[3],s2_imag[3];
      for(dd = 0;dd<3;dd++) {
	normal_imag[dd] = 0.5*xnormal[dd].i/eps;
	s1_imag[dd] = 0.5*xs1[dd].i/eps;
	s2_imag[dd] = 0.5*xs2[dd].i/eps;
      }
      nn = 0;
      for(;nn<3;nn++) {
	for(dd = 0;dd<3;dd++) {
	  KExt_HH[ii+9*3*nn+9*dd] += g_w[gg]*N[3*gg+nn]*normal_imag[dd]*traction[2];
	  KExt_HH[ii+9*3*nn+9*dd] += g_w[gg]*N[3*gg+nn]*s1_imag[dd]*traction[0];
	  KExt_HH[ii+9*3*nn+9*dd] += g_w[gg]*N[3*gg+nn]*s2_imag[dd]*traction[1];
	}
      }
    }
  }
  PetscFunctionReturn(0);
}

