/* Copyright (C) 2009, Lukasz Kaczmarczyk (likask AT civil.gla.ac.uk)
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


static enum phisical_equation_volume ph_eq_vol = hooke;
void set_PhysicalEquationNumber(enum phisical_equation_volume eq) {
  ph_eq_vol = eq;
}
enum phisical_equation_volume get_PhysicalEquationNumber() {
  return ph_eq_vol;
}

static enum thremal_deformation_equation themp_eq_deformation = linear_expanison;
void set_ThermalDeformationEquationNumber(enum thremal_deformation_equation eq) {
  themp_eq_deformation = eq;
}
enum thremal_deformation_equation get_ThermalDeformationEquationNumber() {
  return themp_eq_deformation;
}

static PetscErrorCode ierr;

//Phusical Equations
static PetscErrorCode StrainEnergy_Hooke(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xPsi,void *ctx);
static PetscErrorCode PiolaKirhoiff2_Hooke(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xS,void *ctx);
static PetscErrorCode StrainEnergy_Kirchhoff(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xPsi,void *ctx);
static PetscErrorCode PiolaKirhoiff2_Kirchhoff(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xS,void *ctx);
static PetscErrorCode StrainEnergy_NeoHookean(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xPsi,void *ctx);
static PetscErrorCode PiolaKirhoiff2_NeoHookean(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xS,void *ctx);
static PetscErrorCode StrainEnergy_EberleinHolzapfel1(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xPsi,void *ctx);
static PetscErrorCode PiolaKirhoiff2_EberleinHolzapfel1(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xS,void *ctx);
//
static PetscErrorCode (*tab_func_Strain_energy[])(double,double,__CLPK_doublecomplex*,__CLPK_doublecomplex*,__CLPK_doublecomplex*,__CLPK_doublecomplex*,void*) = {
  StrainEnergy_Hooke,StrainEnergy_Kirchhoff,StrainEnergy_NeoHookean,StrainEnergy_EberleinHolzapfel1
};
static PetscErrorCode (*tab_func_PiolaKirhoiff2[])(double lambda,double mu,__CLPK_doublecomplex*,__CLPK_doublecomplex*,__CLPK_doublecomplex*,__CLPK_doublecomplex*,void*) = {
  PiolaKirhoiff2_Hooke,PiolaKirhoiff2_Kirchhoff,PiolaKirhoiff2_NeoHookean,PiolaKirhoiff2_EberleinHolzapfel1
};
//
static PetscErrorCode ThermalDeformationGradient_linear_expanison(const double alpha,
  const double lambda,const double i_lambda,__CLPK_doublecomplex xT,__CLPK_doublecomplex *xF) {
  PetscFunctionBegin;
  bzero(xF,sizeof(__CLPK_doublecomplex)*9);
  double complex streach = 1 + alpha*( (lambda+I*i_lambda)*(xT.r+I*xT.i) );
  int dd = 0;
  for(;dd<3;dd++) {
    xF[3*dd+dd].r = creal(streach);
    xF[3*dd+dd].i = cimag(streach);
  }
  PetscFunctionReturn(0);
}
static PetscErrorCode ThermalDeformationGradient_linear_expansion_true_volume(const double alpha,
  const double lambda,const double i_lambda,__CLPK_doublecomplex xT,__CLPK_doublecomplex *xF) {
  PetscFunctionBegin;
  bzero(xF,sizeof(__CLPK_doublecomplex)*9);
  double complex streach = cexp( alpha*(lambda+I*i_lambda)*(xT.r+I*xT.i) );
  int dd = 0;
  for(;dd<3;dd++) {
    xF[3*dd+dd].r = creal(streach);
    xF[3*dd+dd].i = cimag(streach);
  }
  PetscFunctionReturn(0);
}
static PetscErrorCode (*tab_func_ThermalDeformationGradient[])(
  const double alpha,const double lambda,const double i_lambda,__CLPK_doublecomplex xT,__CLPK_doublecomplex *xF)  = { 
  ThermalDeformationGradient_linear_expanison, ThermalDeformationGradient_linear_expansion_true_volume
};
PetscErrorCode ThermalDeformationGradient(const double alpha,const double lambda,const double i_lambda,__CLPK_doublecomplex xT,__CLPK_doublecomplex *xF) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = tab_func_ThermalDeformationGradient[themp_eq_deformation](alpha,lambda,i_lambda,xT,xF); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
//
PetscErrorCode StrainEnergy(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xPsi,void *ctx) {
  PetscFunctionBegin;
  ierr = tab_func_Strain_energy[ph_eq_vol](lambda,mu,xF,xC,xJ,xPsi,ctx); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode PiolaKirhoiff2(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xS,void *ctx) {
  PetscFunctionBegin;
  ierr = tab_func_PiolaKirhoiff2[ph_eq_vol](lambda,mu,xF,xC,xJ,xS,ctx); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
//*************************//
PetscErrorCode CauchyStress(__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xP,__CLPK_doublecomplex *xCauchyStress) {
  PetscFunctionBegin;
  if(ph_eq_vol == hooke) {
    cblas_zcopy(9,xP,1,xCauchyStress,1);
    PetscFunctionReturn(0);
  } else {
    __CLPK_doublecomplex tmp2 = {0,0};
    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,3,3,3,xJ,xP,3,xF,3,&tmp2,xCauchyStress,3);
    PetscFunctionReturn(0);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode PiolaKirhoiff1(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xS,__CLPK_doublecomplex *xP) {
  PetscFunctionBegin;
  if(ph_eq_vol == hooke) {
    cblas_zcopy(9,xS,1,xP,1);
    PetscFunctionReturn(0);
  } else {
    __CLPK_doublecomplex tmp1 = {1,0},tmp2 = {0,0};
    cblas_zsymm(CblasRowMajor,CblasRight,CblasUpper,3,3,&tmp1,xS,3,xF,3,&tmp2,xP,3);
    PetscFunctionReturn(0);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ElshebyStress(__CLPK_doublecomplex *xPsi,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xP,__CLPK_doublecomplex *xSigma) {
  PetscFunctionBegin;
  __CLPK_doublecomplex tmp1 = {-1,0},tmp2 = {0,0};
  cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,3,3,3,&tmp1,xF,3,xP,3,&tmp2,xSigma,3);
  xSigma[3*0+0].r += (*xPsi).r;
  xSigma[3*1+1].r += (*xPsi).r;
  xSigma[3*2+2].r += (*xPsi).r;
  xSigma[3*0+0].i += (*xPsi).i;
  xSigma[3*1+1].i += (*xPsi).i;
  xSigma[3*2+2].i += (*xPsi).i;
  PetscFunctionReturn(0);
}
//*************************//
//Hooke
static PetscErrorCode StrainEnergy_Hooke(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xPsi,void *ctx) {
  PetscFunctionBegin;
  double complex Psi = 0,tr_epsilon = 0;
  int ii = 0;
  for(;ii<3;ii++) {
    int jj = 0;
    for(;jj<3;jj++) {
      double complex epsilon = xF[ii*3+jj].r + I*xF[ii*3+jj].i; 
      if(ii == jj) {
	epsilon -= 1;
	tr_epsilon += epsilon;
      } else {
	epsilon = 0.5*(epsilon+xF[jj*3+ii].r + I*xF[jj*3+ii].i);
      }
      Psi += mu*epsilon*epsilon; 
  }}
  Psi += 0.5*lambda*tr_epsilon*tr_epsilon;
  (*xPsi).r = creal(Psi);
  (*xPsi).i = cimag(Psi);
  PetscFunctionReturn(0);
}
static PetscErrorCode PiolaKirhoiff2_Hooke(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xS,void *ctx) {
  PetscFunctionBegin;
  bzero(xS,sizeof(__CLPK_doublecomplex)*9);
  double complex tr_epsilon = 0;
  int ii = 0;
  for(;ii<3;ii++) {
    int jj = 0;
    for(;jj<3;jj++) {
      double complex epsilon = xF[ii*3+jj].r + I*xF[ii*3+jj].i; 
      if(ii == jj) {
	epsilon -= 1;
	tr_epsilon += epsilon;
      } else {
	epsilon = 0.5*(epsilon+xF[jj*3+ii].r + I*xF[jj*3+ii].i);
      }
      double complex stress = 2*mu*epsilon;
      xS[3*ii+jj].r = creal(stress);
      xS[3*ii+jj].i = cimag(stress);
    }
  }
  ii = 0;
  for(;ii<3;ii++) {
    xS[3*ii+ii].r += lambda*creal(tr_epsilon);
    xS[3*ii+ii].i += lambda*cimag(tr_epsilon);
  }
  PetscFunctionReturn(0);
}
//Kirchhoff
static PetscErrorCode StrainEnergy_Kirchhoff(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xPsi,void *ctx) {
  PetscFunctionBegin;
  double complex Psi = 0,tr_E = 0;
  int ii = 0;
  for(;ii<3;ii++) {
    int jj = 0;
    for(;jj<3;jj++) {
      double complex E = 0.5*(xC[ii*3+jj].r + I*xC[ii*3+jj].i); 
      if(ii == jj) {
	E -= 0.5;
	tr_E += E;
      }
      Psi += mu*E*E; 
  }}
  Psi += 0.5*lambda*tr_E*tr_E;
  (*xPsi).r = creal(Psi);
  (*xPsi).i = cimag(Psi);
  PetscFunctionReturn(0);
}
static PetscErrorCode PiolaKirhoiff2_Kirchhoff(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xS,void *ctx) {
  PetscFunctionBegin;
  bzero(xS,sizeof(__CLPK_doublecomplex)*9);
  double complex tr_E = 0;
  int ii = 0;
  for(;ii<3;ii++) {
    int jj = 0;
    for(;jj<3;jj++) {
      double complex E = 0.5*(xC[ii*3+jj].r + I*xC[ii*3+jj].i); 
      if(ii == jj) {
	E -= 0.5;
	tr_E += E;
      }
      double complex stress = 2*mu*E;
      xS[3*ii+jj].r = creal(stress);
      xS[3*ii+jj].i = cimag(stress);
    }
  }
  ii = 0;
  for(;ii<3;ii++) {
    xS[3*ii+ii].r += lambda*creal(tr_E);
    xS[3*ii+ii].i += lambda*cimag(tr_E);
  }
  PetscFunctionReturn(0);
}
//Neo-Hookean form Bonet Book
//Wc = (mu/2)*(Ic-3) + (lambda/2)*(ln(J))^2
static PetscErrorCode StrainEnergy_NeoHookean(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xPsi,void *ctx) {
  PetscFunctionBegin;
  int ii = 0;
  double complex Ic = 0;
  for(;ii<3;ii++) {
    Ic += xC[ii*3+ii].r + I*xC[ii*3+ii].i; }
  double complex J = (*xJ).r + I*(*xJ).i;
  double complex logJ = clog(J);
  double complex Psi = (mu/2.)*(Ic-3.) - mu*logJ + (lambda/2.)*cpow(logJ,2.);
  (*xPsi).r = creal(Psi);
  (*xPsi).i = cimag(Psi);
  PetscFunctionReturn(0);
}
//S= 2*dPsi/dC = mu*(I-invC) + lambda*ln(J)*invC
static PetscErrorCode PiolaKirhoiff2_NeoHookean(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xS,void *ctx) {
  PetscFunctionBegin;
  __CLPK_doublecomplex inv_xC[9];
  cblas_zcopy(9,xC,1,inv_xC,1);
  InvertComplexGradient(inv_xC);
  bzero(xS,sizeof(__CLPK_doublecomplex)*9);
  int ii = 0;
  for(;ii<3;ii++) xS[ii*3+ii].r = mu;
  __CLPK_doublecomplex tmp1 = {-mu,0};
  cblas_zaxpy(9,&tmp1,inv_xC,1,xS,1);
  double complex J = (*xJ).r + I*(*xJ).i;
  double complex logJ = clog(J);
  double complex xtmp2 = lambda*logJ;
  __CLPK_doublecomplex tmp2 = {creal(xtmp2),cimag(xtmp2)};
  cblas_zaxpy(9,&tmp2,inv_xC,1,xS,1);
  PetscFunctionReturn(0);
}
//Fibres, by R. Eberlein, G.A. Holzapfel, and C.A.J. Schulze-Bauer. 
//An anisotropic model for annulus tissue and enhanced finite element analyses of intact lumbar disc bodies. 
//Computer Methods in Biomechanics and Biomedical Engineering, 4(3):209–229, 2001. ISSN 1025-5842.
static PetscErrorCode preProc_EberleinHolzapfel1(__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,void *ctx,
  __CLPK_doublecomplex *xA1,__CLPK_doublecomplex *xA2,__CLPK_doublecomplex *xI1,__CLPK_doublecomplex *xI2) {
  PetscFunctionBegin;
  ctx_EberleinHolzapfel1 *my_ctx = (ctx_EberleinHolzapfel1*)ctx;
  double A1[9],A2[9];
  bzero(A1,sizeof(double)*9);
  bzero(A2,sizeof(double)*9);
  cblas_dger(CblasRowMajor,3,3,1.,(*my_ctx).fibre_vector_a1,1,(*my_ctx).fibre_vector_a1,1,A1,3);
  cblas_dger(CblasRowMajor,3,3,1.,(*my_ctx).fibre_vector_a2,1,(*my_ctx).fibre_vector_a2,1,A2,3);
  double ZERO[9];
  bzero(ZERO,sizeof(double)*9);
  MakeComplexTensor(A1,ZERO,xA1);   
  MakeComplexTensor(A2,ZERO,xA2);
  (*xI1).r = 0;
  (*xI1).i = 0;
  (*xI2).r = 0;
  (*xI2).i = 0;
  double complex pow_xJ = cpow((*xJ).r+I*(*xJ).i,2./3.);
  int ii = 0;
  for(;ii<3;ii++) {
    int jj = 0;
    for(;jj<3;jj++) {
      double complex cx = pow_xJ*(xC[ii*3+jj].r+I*xC[ii*3+jj].i);
      double complex t1 = cx*(xA1[ii*3+jj].r+I*xA1[ii*3+jj].i);
      (*xI1).r += creal(t1);
      (*xI1).i += cimag(t1);
      double complex t2 = cx*(xA2[ii*3+jj].r+I*xA2[ii*3+jj].i);
      (*xI2).r += creal(t2);
      (*xI2).i += cimag(t2);
  }}
  PetscFunctionReturn(0);
}
static PetscErrorCode StrainEnergy_EberleinHolzapfel1(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xPsi,void *ctx) {
  PetscFunctionBegin;
  ctx_EberleinHolzapfel1 *my_ctx = (ctx_EberleinHolzapfel1*)ctx;
  ierr = tab_func_Strain_energy[(*my_ctx).eq_solid](lambda,mu,xF,xC,xJ,xPsi,ctx); CHKERRQ(ierr);
  __CLPK_doublecomplex xA1[9],xA2[9],xI1,xI2;
  preProc_EberleinHolzapfel1(xF,xC,xJ,ctx,xA1,xA2,&xI1,&xI2);
  double k1 = (*my_ctx).k1;
  double k2 = (*my_ctx).k2;
  double complex Psi_f = 0;
  if(xI1.r>1) {
    Psi_f += (k1/(2*k2))*(cexp(k2*cpow(xI1.r+I*xI1.i-1,2))-1.);
  }
  if(xI2.r>1) {
    Psi_f += (k1/(2*k2))*(cexp(k2*cpow(xI2.r+I*xI2.i-1,2))-1.);
  }
  (*xPsi).r += creal(Psi_f);
  (*xPsi).i += cimag(Psi_f);
  PetscFunctionReturn(0);
}
static PetscErrorCode PiolaKirhoiff2_EberleinHolzapfel1(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xS,void *ctx) {
  PetscFunctionBegin;
  ctx_EberleinHolzapfel1 *my_ctx = (ctx_EberleinHolzapfel1*)ctx;
  ierr = tab_func_PiolaKirhoiff2[(*my_ctx).eq_solid](lambda,mu,xF,xC,xJ,xS,ctx); CHKERRQ(ierr);
  __CLPK_doublecomplex xA1[9],xA2[9],xI1,xI2;
  preProc_EberleinHolzapfel1(xF,xC,xJ,ctx,xA1,xA2,&xI1,&xI2);
  double k1 = (*my_ctx).k1;
  double k2 = (*my_ctx).k2;
  if(xI1.r>1) {
    double complex c = 2*k1*cexp(k2*cpow(xI1.r+I*xI1.i-1,2))*(xI1.r+I*xI1.i-1);
   __CLPK_doublecomplex xc = { creal(c), cimag(c) }; 
    cblas_zaxpy(9,&xc,xA1,1,xS,1);
  }
  if(xI2.r>1) {
    double complex c = 2*k1*cexp(k2*cpow(xI2.r+I*xI2.i-1,2))*(xI2.r+I*xI2.i-1);
   __CLPK_doublecomplex xc = { creal(c), cimag(c) }; 
    cblas_zaxpy(9,&xc,xA2,1,xS,1);
  }
  PetscFunctionReturn(0);
}



