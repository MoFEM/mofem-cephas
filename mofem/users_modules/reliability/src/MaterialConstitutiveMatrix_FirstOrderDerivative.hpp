/* Copyright (C) 2014, 
 *   Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
 *   Xiao-Yi Zhou (xiaoyi.zhou AT newcastle.ac.uk)
 * --------------------------------------------------------------
 * This routine calculates the first-order partial derivative of constitutive 
 * matrix, D_r, with respect to Young's modulus in z-direction, E_z, Young's
 * modulus p-direction, E_p, Poisson's ratio in p-direction, NU_p, Poisson's 
 * ratio in z-direction, NU_pz, and shear modulus in z-direction, G_zp, for
 * transversely isotropic material, which is usually for fibre/inclusion in
 * composite material, in the principal material coordinate system.
 *
 * HISTORY
 *
 * 2014.09.08 (first version)
 *
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

#ifndef __MATERIALCONSTITUTIVEMATRIX_FIRSTORDERDERIVATIVE_HPP__
#define __MATERIALCONSTITUTIVEMATRIX_FIRSTORDERDERIVATIVE_HPP__

#include <boost/numeric/ublas/symmetric.hpp>

namespace MoFEM {
  
// =============================================================================
//
//  TRANSVERSELY ISOTROPIC MATERIAL
//
// =============================================================================  
  struct TransverseIsotropicStiffnessMatrix_FirstOrderDerivative {
    double nu_p, nu_pz, E_p, E_z, G_zp;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rPoissonP;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rPoissonPZ;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rYoungP;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rYoungZ;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rShearZP;
    /***************************************************************************
     *
     * With repect to Poisson's ration in p-direction
     * 
     **************************************************************************/
      virtual PetscErrorCode D_r_PoissonP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      // cout<<"Hello from D_r_ElasticFEMethodTransIso"<<endl;
      StiffnessMatrix_rPoissonP.resize(6);
      StiffnessMatrix_rPoissonP.clear();
  
      StiffnessMatrix_rPoissonP(0,0) = pow(E_p,3)/(2*pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2))-E_p/(2*pow(nu_p+1,2));
      StiffnessMatrix_rPoissonP(0,1) = pow(E_p,3)/(2*pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2))+E_p/(2*pow(nu_p+1,2));
      StiffnessMatrix_rPoissonP(0,2) = (pow(E_p,2)*E_z*nu_pz)/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);

      //StiffnessMatrix_rPoissonP(1,0) = StiffnessMatrix_rPoissonP(0,1);
      StiffnessMatrix_rPoissonP(1,1) = StiffnessMatrix_rPoissonP(0,0);
      StiffnessMatrix_rPoissonP(1,2) = StiffnessMatrix_rPoissonP(0,2);

      //StiffnessMatrix_rPoissonP(2,0) = StiffnessMatrix_rPoissonP(0,2);
      //StiffnessMatrix_rPoissonP(2,1) = StiffnessMatrix_rPoissonP(0,2);
      StiffnessMatrix_rPoissonP(2,2) = (2*E_p*pow(E_z*nu_pz,2))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);

      StiffnessMatrix_rPoissonP(3,3) = -E_p/(2*pow((nu_p + 1),2));

      StiffnessMatrix_rPoissonP(4,4) = 0;

      StiffnessMatrix_rPoissonP(5,5) = 0;
      PetscFunctionReturn(0);
    }

    /***************************************************************************
     *
     * With repect to Poisson's ration in z-direction, nu_pz
     * 
     **************************************************************************/
      virtual PetscErrorCode D_r_PoissonPZ(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from TransverseIsotropicStiffnessMatrix_rPoissonPZ "<<endl;
      StiffnessMatrix_rPoissonPZ.resize(6);
      StiffnessMatrix_rPoissonPZ.clear();
  
      StiffnessMatrix_rPoissonPZ(0,0) = (2*E_p*E_p*E_z*nu_pz)/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);
      StiffnessMatrix_rPoissonPZ(0,1) = StiffnessMatrix_rPoissonPZ(0,0);
      StiffnessMatrix_rPoissonPZ(0,2) = (E_p*E_z*(2*E_z*pow(nu_pz,2)+E_p-E_p*nu_p))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);

      //StiffnessMatrix_rPoissonPZ(1,0) = StiffnessMatrix_rPoissonPZ(0,1);
      StiffnessMatrix_rPoissonPZ(1,1) = StiffnessMatrix_rPoissonPZ(0,0);
      StiffnessMatrix_rPoissonPZ(1,2) = StiffnessMatrix_rPoissonPZ(0,2);

      //StiffnessMatrix_rPoissonPZ(2,0) = StiffnessMatrix_rPoissonPZ(0,2);
      //StiffnessMatrix_rPoissonPZ(2,1) = StiffnessMatrix_rPoissonPZ(0,2);
      StiffnessMatrix_rPoissonPZ(2,2) = (4*E_p*E_z*E_z*nu_pz*(1-nu_p))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);

      StiffnessMatrix_rPoissonPZ(3,3) = 0;

      StiffnessMatrix_rPoissonPZ(4,4) = 0;

      StiffnessMatrix_rPoissonPZ(5,5) = 0;
      PetscFunctionReturn(0);
    }

    /***************************************************************************
     *
     * With repect to Young's modulus in p-direction, E_p
     * 
     **************************************************************************/
      virtual PetscErrorCode D_r_YoungP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
//      cout<<"Hello from TransverseIsotropicStiffnessMatrix_rPoissonP "<<endl;
      StiffnessMatrix_rYoungP.resize(6);
      StiffnessMatrix_rYoungP.clear();
  
      StiffnessMatrix_rYoungP(0,0) = (2*pow(E_z,2)*pow(nu_pz,4))/((nu_p-1)*pow((2*E_z*nu_pz*nu_pz+E_p*(nu_p-1)),2))-1/(pow(nu_p,2)-1);
      StiffnessMatrix_rYoungP(0,1) = (2*pow(E_z,2)*pow(nu_pz,4))/((nu_p-1)*pow((2*E_z*nu_pz*nu_pz+E_p*(nu_p-1)),2))-nu_p/(pow(nu_p,2)-1);
      StiffnessMatrix_rYoungP(0,2) = -(2*E_z*E_z*pow(nu_pz,3))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);

      //StiffnessMatrix_rYoungP(1,0) = StiffnessMatrix_rYoungP(0,1);
      StiffnessMatrix_rYoungP(1,1) = StiffnessMatrix_rYoungP(0,0);
      StiffnessMatrix_rYoungP(1,2) = StiffnessMatrix_rYoungP(0,2);

      //StiffnessMatrix_rYoungP(2,0) = StiffnessMatrix_rYoungP(0,2);
      //StiffnessMatrix_rYoungP(2,1) = StiffnessMatrix_rYoungP(0,2);
      StiffnessMatrix_rYoungP(2,2) = (2*pow(E_z*nu_pz,2)*(nu_p-1))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);

      StiffnessMatrix_rYoungP(3,3) = 1/(2*(nu_p + 1));

      StiffnessMatrix_rYoungP(4,4) = 0;

      StiffnessMatrix_rYoungP(5,5) = 0;
      PetscFunctionReturn(0);
    }

    /***************************************************************************
     *
     * With repect to Young's modulus in z-direction, E_z
     * 
     **************************************************************************/
      virtual PetscErrorCode D_r_YoungZ(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      // cout<<"Hello from TransverseIsotropicStiffnessMatrix_rYoungZ "<<endl;
      StiffnessMatrix_rYoungZ.resize(6);
      StiffnessMatrix_rYoungZ.clear();
  
      StiffnessMatrix_rYoungZ(0,0) = pow(E_p*nu_pz,2)/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);
      StiffnessMatrix_rYoungZ(0,1) = pow(E_p*nu_pz,2)/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);
      StiffnessMatrix_rYoungZ(0,2) = (pow(E_p,2)*nu_pz*(1-nu_p))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);

      //StiffnessMatrix_rYoungZ(1,0) = StiffnessMatrix_rYoungZ(0,0);
      StiffnessMatrix_rYoungZ(1,1) = StiffnessMatrix_rYoungZ(0,1);
      StiffnessMatrix_rYoungZ(1,2) = StiffnessMatrix_rYoungZ(0,2);

      //StiffnessMatrix_rYoungZ(2,0) = StiffnessMatrix_rYoungZ(0,2);
      //StiffnessMatrix_rYoungZ(2,1) = StiffnessMatrix_rYoungZ(0,2);
      StiffnessMatrix_rYoungZ(2,2) = pow(E_p*(nu_p-1),2)/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);

      StiffnessMatrix_rYoungZ(3,3) = 0;

      StiffnessMatrix_rYoungZ(4,4) = 0;

      StiffnessMatrix_rYoungZ(5,5) = 0;
      //cout<<StiffnessMatrix_rYoungZ<<endl;
      //cout<<nu_p<<"\t"<<nu_pz<<"\t"<<E_p<<"\t"<<E_z<<"\t"<<G_zp<<"\n";
      PetscFunctionReturn(0);
    }

    /***************************************************************************
     *
     * With repect to Shear modulus in z-direction, G_zp
     * 
     **************************************************************************/
      virtual PetscErrorCode D_r_ShearZP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from TransverseIsotropicStiffnessMatrix_rShearZP"<<endl;
      StiffnessMatrix_rShearZP.resize(6);
      StiffnessMatrix_rShearZP.clear();

      StiffnessMatrix_rShearZP(4,4) = 1;

      StiffnessMatrix_rShearZP(5,5) = 1;
      PetscFunctionReturn(0);
    }
    
  };
// =============================================================================
//
//  ISOTROPIC MATERIAL
//
// =============================================================================
  struct IsotropicStiffnessMatrix_FirstOrderDerivative {
    double young, nu;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rPoisson;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rYoung;
    /***************************************************************************
     *
     * With repect to Young's modulus
     * 
     **************************************************************************/
    virtual PetscErrorCode D_r_Young(double young, double nu){
      PetscFunctionBegin;
      StiffnessMatrix_rYoung.resize(6);
      StiffnessMatrix_rYoung.clear();

      double D00,D01,D33,constt;
      constt = 1/((1+nu)*(1-2*nu));
      
      D00 = constt*(1-nu);
      D01 = constt*nu;
      D33 = constt*(1-2*nu)/2;
      
      StiffnessMatrix_rYoung(0,0) = D00;  
      StiffnessMatrix_rYoung(0,1) = D01;  
      StiffnessMatrix_rYoung(0,2) = D01;

      //StiffnessMatrix_rYoung(1,0) = D01;  
      StiffnessMatrix_rYoung(1,1) = D00;  
      StiffnessMatrix_rYoung(1,2) = D01;

      //StiffnessMatrix_rYoung(2,0) = D01;  
      //StiffnessMatrix_rYoung(2,1) = D01;   
      StiffnessMatrix_rYoung(2,2) = D00;

      StiffnessMatrix_rYoung(3,3) = D33;

      StiffnessMatrix_rYoung(4,4) = D33;

      StiffnessMatrix_rYoung(5,5) = D33;

      PetscFunctionReturn(0);
    }
    /***************************************************************************
     *
     * With repect to Poisson's ratio
     * 
     **************************************************************************/
    virtual PetscErrorCode D_r_Poisson(double young, double nu){
      PetscFunctionBegin;
      StiffnessMatrix_rPoisson.resize(6);
      StiffnessMatrix_rPoisson.clear();

      double D00,D01,D33;
      
      D00 = -(2*young*nu*(nu - 2))/pow((2*nu*nu + nu - 1),2);
      D01 = (young*(2*nu*nu + 1))/pow((2*nu*nu + nu - 1),2);
      D33 = -young/(2*pow((nu + 1),2));
      
      StiffnessMatrix_rPoisson(0,0) = D00;  
      StiffnessMatrix_rPoisson(0,1) = D01;  
      StiffnessMatrix_rPoisson(0,2) = D01;

      //StiffnessMatrix_rPoisson(1,0) = D01;  
      StiffnessMatrix_rPoisson(1,1) = D00;  
      StiffnessMatrix_rPoisson(1,2) = D01;

      //StiffnessMatrix_rPoisson(2,0) = D01;  
      //StiffnessMatrix_rPoisson(2,1) = D01;   
      StiffnessMatrix_rPoisson(2,2) = D00;

      StiffnessMatrix_rPoisson(3,3) = D33;

      StiffnessMatrix_rPoisson(4,4) = D33;

      StiffnessMatrix_rPoisson(5,5) = D33;
      PetscFunctionReturn(0);
    }
};
  
  // ===========================================================================
  //
  // TWO-PHASE YARN - TRANSVERSELY ISOTROPIC MATERIAL
  //
  // ===========================================================================
  struct YarnStiffnessMatrix_Geom_FirstOrderDerivative{
	// fibre waviness
	double dI1, dI3, dI5, dI6, dI8;
//	double ampltidue, length; // waviness parameters: waviness amplitude & periodic length
//	// misalignment angle
//	double theta;
	double dmn, dm2, dn2, dm2n2, dmn3, dm3n, dn4, dm4; 
	// fibre volume fraction
//	double vf;
//	double kf, lf, mf, nf, pf; // fibre related Hill's moduli
//	double km, lm, mm, nm, pm; // matrix related Hill's moduli
	double dkc, dmc, dpc, dlc, dnc; // fibre volume fraction
    
    /***************************************************************************
     *
     * With repect to amplitude of waviness
     *
     **************************************************************************/
    virtual PetscErrorCode D_r_Amplitude(double amplitude, double length){
      PetscFunctionBegin;
      
      double alpha;
      alpha = 2*M_PI*amplitude/length;
      
      dI1 = -(pow(alpha,3) + 4*alpha)/(2*pow(alpha*alpha + 1,2.5))*(2*M_PI/length);
      dI3 = -(pow(alpha,3) - 2*alpha)/(2*pow(alpha*alpha + 1,2.5))*(2*M_PI/length);
      dI5 = 3*pow(alpha,3)/(2*pow(alpha*alpha + 1,2.5))*(2*M_PI/length);
      dI6 = -alpha/pow(alpha*alpha + 1,1.5)*(2*M_PI/length);
      dI8 =  alpha/pow(alpha*alpha + 1,1.5)*(2*M_PI/length);
			
      PetscFunctionReturn(0);
    }
    
    /***************************************************************************
     *
     * With repect to periodic length
     *
     **************************************************************************/
    virtual PetscErrorCode D_r_Length(double amplitude, double length){
      PetscFunctionBegin;
      
      double alpha;
      alpha = 2*M_PI*amplitude/length;
      
      dI1 = -(pow(alpha,3) + 4*alpha)/(2*pow(alpha*alpha + 1,2.5))*(-2*M_PI*amplitude/length/length);
      dI3 = -(pow(alpha,3) - 2*alpha)/(2*pow(alpha*alpha + 1,2.5))*(-2*M_PI*amplitude/length/length);
      dI5 = 3*pow(alpha,3)/(2*pow(alpha*alpha + 1,2.5))*(-2*M_PI*amplitude/length/length);
      dI6 = -alpha/pow(alpha*alpha + 1,1.5)*(-2*M_PI*amplitude/length/length);
      dI8 =  alpha/pow(alpha*alpha + 1,1.5)*(-2*M_PI*amplitude/length/length);
			
      PetscFunctionReturn(0);
    }
 
    /***************************************************************************
     *
     * With repect to misalignment angle
     *
     **************************************************************************/  
    virtual PetscErrorCode D_r_Angle(double theta) {
      PetscFunctionBegin;
			
      dmn = -sin(theta)*sin(theta)+cos(theta)*cos(theta); // mn
      dm2 = -2*cos(theta)*sin(theta); // m2
      dn2 =  2*sin(theta)*cos(theta);
      dmn3 = -pow(sin(theta),4) + 3*pow(cos(theta),2)*pow(sin(theta),2); // mn3
      dm3n =  pow(cos(theta),4) - 3*pow(cos(theta),2)*pow(sin(theta),2); // nm3
      dm2n2 = -2*cos(theta)*pow(sin(theta),3) + 2*sin(theta)*pow(cos(theta),3); // m2n2
      dm4 = -4*pow(cos(theta),3)*sin(theta); // m4
      dn4 = 4*pow(sin(theta),3)*cos(theta);
			
      PetscFunctionReturn(0);
	}
      
	/***************************************************************************
	 *
	 * With repect to misalignment angle
	 *
	 **************************************************************************/
	virtual PetscErrorCode D_r_Fraction(double vf, 
										double kf, double mf, double pf, double lf, double nf,
										double km, double mm, double pm, double lm, double nm) {
      PetscFunctionBegin;
			
      dkc = ((kf + mm)*(km + mm)*(kf - km))/pow(kf + mm - vf*kf + vf*km,2);
      dmc = (km*mm*(mf - mm))/(km*mm + (vf*mm - mf*(vf - 1))*(km + 2*mm)) 
			+ ((km + 2*mm)*(mf - mm)*(km*mm*(mf*vf - mm*(vf - 1)) + mf*mm*(km + 2*mm)))
			/pow((km*mm + (vf*mm - mf*(vf - 1))*(km + 2*mm)),2);
      dpc = (2*pm*(pf*pf - pm*pm))/pow(pf + pm - vf*pf + vf*pm,2);
	  dlc = -((km + mm)*(kf*lm + lm*mm - lf*pf - lf*pm))/pow(pf + pm + km*vf + mm*vf - pf*vf - pm*vf,2);
	  dnc = nf - nm - ((lf - lm)*(lf - lm - (lf*(km + mm) - lm*(kf + mm))
		  /(vf*(km + mm) - (pf + pm)*(vf - 1)) 
		  + ((lf*vf*(km + mm) - lm*(kf + mm)*(vf - 1))*(km + mm - pf - pm))
		  /pow((vf*(km + mm) - (pf + pm)*(vf - 1)),2)))/(kf - km);
 	
	  PetscFunctionReturn(0);
	}
	
    /***************************************************************************
     *
     * With repect to Young's modulus in z-direction, E_z
     * 
     **************************************************************************/
	virtual PetscErrorCode D_r_YoungZ(double nu_p, double nu_pz, double E_p, double E_z, double G_zp,
	                                  double Em, double NUm, double vf){
		
      PetscFunctionBegin;
      
	  dkc = (4*vf*pow(E_p*Em*nu_pz*(NUm - 1),2))
	        /pow((E_p*E_p*NUm + E_p*E_p*vf - E_p*E_p + 2*pow(E_p*NUm,2) - E_p*Em 
			- 2*pow(E_p*NUm,2)*vf + 2*E_p*Em*NUm + E_p*Em*nu_p - E_p*Em*vf 
			+ 2*E_z*Em*nu_pz*nu_pz - E_p*E_p*NUm*vf - 4*E_z*Em*NUm*nu_pz*nu_pz 
			+ 2*E_z*Em*nu_pz*nu_pz*vf - 2*E_p*Em*NUm*nu_p + E_p*Em*nu_p*vf),2);
      dmc = 0;
	  dpc = 0;
	  dlc = (2*E_p*E_p*Em*nu_pz*(vf + NUm*nu_pz - NUm*vf - nu_p*vf + NUm*nu_p*vf - NUm*nu_pz*vf))
	        /(pow((2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p),2)*(Em + 2*G_zp + Em*vf 
			- 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf));
	  dnc = (E_p*vf*(nu_p - 1))/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p) 
	        - (4*nu_pz*nu_pz*((Em*NUm)/(2*NUm*NUm + NUm - 1) 
			- (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
			*((2*((Em*NUm*(vf - 1)*(E_p*E_p*NUm + E_p*E_p + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
			/(2*(2*NUm - 1)*pow((NUm + 1),2)*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			+ (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*nu_pz*nu_pz - 2*E_p + 2*E_p*nu_p)
			*(2*NUm*NUm + NUm - 1)))*(2*NUm*NUm + NUm - 1))
			/(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf) 
			- (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
			+ (E_p*E_z*nu_pz*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))*(2*NUm*NUm + NUm - 1))
			/(E_p*E_p*NUm - E_p*E_p + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz) 
			- (2*E_p*E_z*nu_pz*nu_pz*vf*(nu_p - 1))/pow((E_p*nu_p - E_p + 2*E_z*nu_pz*nu_pz),2)
            - (4*Em*nu_pz*nu_pz*((Em*NUm)/(2*NUm*NUm + NUm - 1) 
			- (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))*((2*((Em*NUm*(vf - 1)
			*(E_p*E_p*NUm + E_p*E_p + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
			/(2*(2*NUm - 1)*pow((NUm + 1),2)*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			+ (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*nu_pz*nu_pz - 2*E_p + 2*E_p*nu_p)
			*(2*NUm*NUm + NUm - 1)))*(2*NUm*NUm + NUm - 1))
			/(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm 
			- 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf) 
			- (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
			+ (E_p*E_z*nu_pz*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
			*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)*(2*NUm*NUm + NUm - 1))
			/pow((E_p*E_p*NUm - E_p*E_p + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz),2) 
			+ (2*E_p*E_p*nu_pz*(nu_p - 1)*((2*((Em*NUm*(vf - 1)
			*(E_p*E_p*NUm + E_p*E_p + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
			/(2*(2*NUm - 1)*pow((NUm + 1),2)*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			+ (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*nu_pz*nu_pz - 2*E_p + 2*E_p*nu_p)
			*(2*NUm*NUm + NUm - 1)))*(2*NUm*NUm + NUm - 1))
			/(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm 
			- 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf) 
			- (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
			+ (E_p*E_z*nu_pz*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
			*(2*NUm*NUm + NUm - 1))/((2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)
			*(E_p*E_p*NUm - E_p*E_p + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz)) 
			+ (2*E_p*E_p*nu_pz*((Em*NUm)/(2*NUm*NUm + NUm - 1) 
			- (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))*(vf - 1)
			*(2*NUm*NUm + NUm - 1)*(Em*vf - 2*G_zp*vf + 2*Em*NUm*nu_pz 
			+ 2*G_zp*NUm*vf - Em*nu_p*vf + 2*G_zp*nu_p*vf + 4*G_zp*NUm*NUm*vf 
			- 4*G_zp*NUm*NUm*nu_p*vf - 2*G_zp*NUm*nu_p*vf))
			/((2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)*(E_p*E_p*NUm - E_p*E_p 
			+ 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz)
		   *(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm 
		   - 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf));
	  
      PetscFunctionReturn(0);
    }
    /***************************************************************************
     *
     * With repect to transverse modulus, E_p
     * 
     **************************************************************************/	
	virtual PetscErrorCode D_r_YoungP(double nu_p, double nu_pz, double E_p, 
	                                  double E_z, double G_zp,
	                                  double Em, double NUm, double vf){
		
      PetscFunctionBegin;
	  
	  dkc = -(2*E_p*vf*pow(Em*(NUm - 1),2)*(4*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
	         /pow((E_p*E_p*NUm + E_p*E_p*vf - E_p*E_p + 2*E_p*E_p*NUm*NUm - E_p*Em 
			 - 2*E_p*E_p*NUm*NUm*vf + 2*E_p*Em*NUm + E_p*Em*nu_p - E_p*Em*vf 
			 + 2*E_z*Em*nu_pz*nu_pz - E_p*E_p*NUm*vf - 4*E_z*Em*NUm*nu_pz*nu_pz 
			 + 2*E_z*Em*nu_pz*nu_pz*vf - 2*E_p*Em*NUm*nu_p + E_p*Em*nu_p*vf),2);
	  dmc = (8*vf*pow(Em*(NUm - 1),2)*(nu_p + 1))
	        /pow((3*E_p + Em + Em*nu_p - 3*E_p*vf + 3*Em*vf - 4*E_p*NUm*NUm 
			- E_p*NUm + E_p*NUm*vf - 4*Em*NUm*vf + 3*Em*nu_p*vf + 4*E_p*NUm*NUm*vf 
			- 4*Em*NUm*nu_p*vf),2);
	  dpc = 0;
	  dlc = (Em*(E_p*E_p*NUm - 4*E_z*E_z*pow(nu_pz,3)*vf - E_p*E_p*NUm*nu_p 
	        - E_p*E_p*NUm*vf - 4*E_p*E_z*NUm*nu_pz*nu_pz + E_p*E_p*NUm*nu_p*vf 
			+ 4*E_z*E_z*NUm*pow(nu_pz,3)*vf + 4*E_p*E_z*NUm*nu_pz*nu_pz*vf))
			/(pow((2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p),2)
			*(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm 
			- 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf));
	  dnc = (E_z*vf*(nu_p - 1))/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p) 
	        - (2*((Em*NUm)/(2*NUm*NUm + NUm - 1) 
			- (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
			*(nu_p - 1)*((2*((Em*NUm*(vf - 1)*(E_p*E_p*NUm + E_p*E_p + E_p*Em 
			- E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
			/(2*(2*NUm - 1)*pow((NUm + 1),2)*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			+ (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*nu_pz*nu_pz - 2*E_p + 2*E_p*nu_p)
			*(2*NUm*NUm + NUm - 1)))*(2*NUm*NUm + NUm - 1))
			/(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm 
			- 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf) 
			- (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
			+ (E_p*E_z*nu_pz*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
			*(2*NUm*NUm + NUm - 1))/(E_p*E_p*NUm - E_p*E_p + 2*E_p*E_p*NUm*NUm 
			+ E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz) 
			- (E_p*E_z*vf*pow((nu_p - 1),2))/pow((E_p*nu_p - E_p + 2*E_z*nu_pz*nu_pz),2) 
			+ (4*E_z*E_z*pow(nu_pz,3)*((2*((Em*NUm*(vf - 1)*(E_p*E_p*NUm 
			+ E_p*E_p + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
			/(2*(2*NUm - 1)*pow((NUm + 1),2)*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			+ (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*nu_pz*nu_pz - 2*E_p + 2*E_p*nu_p)
			*(2*NUm*NUm + NUm - 1)))*(2*NUm*NUm + NUm - 1))
			/(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm 
			- 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf) 
			- (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
			+ (E_p*E_z*nu_pz*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))*(2*NUm*NUm + NUm - 1))
			/((2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)*(E_p*E_p*NUm - E_p*E_p 
			+ 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz)) 
			+ (2*((Em*NUm)/(2*NUm*NUm + NUm - 1) - (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz 
			- E_p + E_p*nu_p))*((2*((Em*NUm*(vf - 1)*(E_p*E_p *NUm + E_p*E_p  
			+ E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
			/(2*(2*NUm - 1)*pow((NUm + 1),2)*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			+ (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*nu_pz*nu_pz - 2*E_p + 2*E_p*nu_p)
			*(2*NUm*NUm + NUm - 1)))*(2*NUm*NUm + NUm - 1))
			/(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm 
			- 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf) 
			- (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
			+ (E_p*E_z*nu_pz*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
			*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)*(2*NUm*NUm + NUm - 1)
			*(4*E_p*NUm*NUm + 2*E_p*NUm - 2*E_p + Em - Em*nu_p))
			/pow((E_p*E_p*NUm - E_p*E_p + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p 
			- 2*E_z*Em*nu_pz*nu_pz),2) - (2*((Em*NUm)/(2*NUm*NUm + NUm - 1) 
			- (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))*(vf - 1)
			*(2*NUm*NUm + NUm - 1)*(E_p*E_p*Em*NUm*nu_p - E_p*E_p*Em*NUm 
			+ 2*E_z*E_p*Em*pow(nu_pz,3)*vf - 4*E_z*E_z*G_zp*pow(nu_pz,3)*vf 
			+ 8*E_z*E_z*G_zp*NUm*NUm*pow(nu_pz,3)*vf + 4*E_p*E_z*Em*NUm*nu_pz*nu_pz 
			+ 4*E_z*E_z*G_zp*NUm*pow(nu_pz,3)*vf))/((2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)
			*(E_p*E_p*NUm - E_p*E_p + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p 
			- 2*E_z*Em*nu_pz*nu_pz)*(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm 
			- 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf));

      PetscFunctionReturn(0);
    }
	
    /***************************************************************************
     *
     * With repect to transverse Poisson's ratio, nu_p
     * 
     **************************************************************************/	
	virtual PetscErrorCode D_r_PoissonP(double nu_p, double nu_pz, double E_p, 
	                                  double E_z, double G_zp,
	                                  double Em, double NUm, double vf){
		
      PetscFunctionBegin;
	  
	  dkc = (2*pow(E_p,3)*Em*Em*vf*(NUm - 1)*(NUm - 1))
	       /pow((E_p*E_p*NUm + E_p*E_p*vf - E_p*E_p + 2*E_p*E_p*NUm*NUm - E_p*Em 
		   - 2*E_p*E_p*NUm*NUm*vf + 2*E_p*Em*NUm + E_p*Em*nu_p - E_p*Em*vf 
		   + 2*E_z*Em*nu_pz*nu_pz - E_p*E_p*NUm*vf - 4*E_z*Em*NUm*nu_pz*nu_pz 
		   + 2*E_z*Em*nu_pz*nu_pz*vf - 2*E_p*Em*NUm*nu_p + E_p*Em*nu_p*vf),2);
	  dmc = -(8*E_p*Em*Em*vf*(NUm - 1)*(NUm - 1))
	        /pow((3*E_p + Em + Em*nu_p - 3*E_p*vf + 3*Em*vf - 4*E_p*NUm*NUm - E_p*NUm 
			+ E_p*NUm*vf - 4*Em*NUm*vf + 3*Em*nu_p*vf + 4*E_p*NUm*NUm*vf - 4*Em*NUm*nu_p*vf),2);
	  dpc = 0;
	  dlc = (E_p*E_p*Em*(E_p*NUm - E_p*NUm*vf + 2*E_z*nu_pz*vf - 2*E_z*NUm*nu_pz*vf))
	        /(pow((2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p),2)
			*(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm 
			- 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf));
	  dnc = (E_p*E_z*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p) 
		   - (E_p*E_p*E_z*vf*(nu_p - 1))/pow((E_p*nu_p - E_p + 2*E_z*nu_pz*nu_pz),2)
		   - (2*E_p*((Em*NUm)/(2*NUm*NUm + NUm - 1) 
		   - (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
		   *((2*((Em*NUm*(vf - 1)*(E_p*E_p*NUm + E_p*E_p + E_p*Em 
		   - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
		   /(2*(2*NUm - 1)*(NUm + 1)*(NUm + 1)*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
		   + (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*nu_pz*nu_pz 
		   - 2*E_p + 2*E_p*nu_p)*(2*NUm*NUm + NUm - 1)))
		   *(2*NUm*NUm + NUm - 1))/(Em + 2*G_zp + Em*vf - 2*G_zp*vf 
		   - 4*G_zp*NUm*NUm - 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf 
		   + 4*G_zp*NUm*NUm*vf) - (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
		   + (E_p*E_z*nu_pz*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
		   *(2*NUm*NUm + NUm - 1))/(E_p*E_p*NUm - E_p*E_p 
		   + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz) 
		   - (2*E_p*E_p*E_z*nu_pz*((2*((Em*NUm*(vf - 1)
		   *(E_p*E_p*NUm + E_p*E_p + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
		   /(2*(2*NUm - 1)*(NUm + 1)*(NUm + 1)*(2*E_z*nu_pz*nu_pz - E_p 
		   + E_p*nu_p)) + (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))
		   /((4*E_z*nu_pz*nu_pz - 2*E_p + 2*E_p*nu_p)
		   *(2*NUm*NUm + NUm - 1)))*(2*NUm*NUm + NUm - 1))
		   /(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm 
		   - 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf) 
		   - (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
		   + (E_p*E_z*nu_pz*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
		   *(2*NUm*NUm + NUm - 1))/((2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)
		   *(E_p*E_p*NUm - E_p*E_p + 2*E_p*E_p*NUm*NUm + E_p*Em 
		   - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz)) 
		   - (2*E_p*Em*((Em*NUm)/(2*NUm*NUm + NUm - 1) 
		   - (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
		   *((2*((Em*NUm*(vf - 1)*(E_p*E_p*NUm + E_p*E_p + E_p*Em 
		   - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
		   /(2*(2*NUm - 1)*(NUm + 1)*(NUm + 1)*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
		   + (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*nu_pz*nu_pz - 2*E_p + 2*E_p*nu_p)
		   *(2*NUm*NUm + NUm - 1)))*(2*NUm*NUm + NUm - 1))/(Em + 2*G_zp 
		   + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm - 2*G_zp*NUm 
		   + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf) 
		   - (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
		   + (E_p*E_z*nu_pz*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
		   *(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)*(2*NUm*NUm + NUm - 1))
		   /pow((E_p*E_p*NUm - E_p*E_p + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p 
		   - 2*E_z*Em*nu_pz*nu_pz),2) + (2*E_p*E_p*((Em*NUm)/(2*NUm*NUm + NUm - 1) 
		   - (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
		   *(vf - 1)*(2*NUm*NUm + NUm - 1)*(E_p*Em*NUm + E_z*Em*nu_pz*vf 
		   - 2*E_z*G_zp*nu_pz*vf + 2*E_z*G_zp*NUm*nu_pz*vf + 4*E_z*G_zp*NUm*NUm*nu_pz*vf))
		   /((2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)*(E_p*E_p*NUm - E_p*E_p 
		   + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz)
		   *(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm 
		   - 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf));
 
      PetscFunctionReturn(0);
    }	
	
    /***************************************************************************
     *
     * With repect to Young's modulus in z-direction, E_z
     * 
     **************************************************************************/
	virtual PetscErrorCode D_r_PoissonPZ(double nu_p, double nu_pz, double E_p, double E_z, double G_zp,
	                                  double Em, double NUm, double vf){
		
	  PetscFunctionBegin;
      
	  dkc = (8*E_p*E_p*E_z*Em*Em*nu_pz*vf*(NUm - 1)*(NUm - 1))
	       /pow((E_p*E_p*NUm + E_p*E_p*vf - E_p*E_p + 2*E_p*E_p*NUm*NUm 
		   - E_p*Em - 2*E_p*E_p*NUm*NUm*vf + 2*E_p*Em*NUm + E_p*Em*nu_p 
		   - E_p*Em*vf + 2*E_z*Em*nu_pz*nu_pz - E_p*E_p*NUm*vf 
		   - 4*E_z*Em*NUm*nu_pz*nu_pz + 2*E_z*Em*nu_pz*nu_pz*vf 
		   - 2*E_p*Em*NUm*nu_p + E_p*Em*nu_p*vf),2);
	  dmc = 0;	
      dpc = 0;  
      dlc = (2*E_p*E_z*Em*(E_p*vf + 2*E_p*NUm*nu_pz - E_p*NUm*vf 
		    - E_p*nu_p*vf + 2*E_z*nu_pz*nu_pz*vf - 2*E_z*NUm*nu_pz*nu_pz*vf 
			+ E_p*NUm*nu_p*vf - 2*E_p*NUm*nu_pz*vf))
			/(pow((2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p),2)
			*(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm 
			- 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf)); 
      dnc = (2*E_p*E_z*((Em*NUm)/(2*NUm*NUm + NUm - 1) 
		   - (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))*(vf - 1)
		   *(2*NUm*NUm + NUm - 1)*(E_p*Em*vf - 2*E_p*G_zp*vf + 4*E_p*G_zp*NUm*NUm*vf 
		   + 2*E_z*Em*nu_pz*nu_pz*vf - 4*E_z*G_zp*nu_pz*nu_pz*vf 
		   + 4*E_p*Em*NUm*nu_pz + 2*E_p*G_zp*NUm*vf - E_p*Em*nu_p*vf 
		   + 2*E_p*G_zp*nu_p*vf - 2*E_p*G_zp*NUm*nu_p*vf 
		   - 4*E_p*G_zp*NUm*NUm*nu_p*vf + 4*E_z*G_zp*NUm*nu_pz*nu_pz*vf 
		   + 8*E_z*G_zp*NUm*NUm*nu_pz*nu_pz*vf))
		   /((2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)*(E_p*E_p*NUm - E_p*E_p 
		   + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz)
		   *(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm- 2*Em*NUm 
		   - 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf)) 
		   - (4*E_p*E_z*E_z*nu_pz*vf*(nu_p - 1))/pow((E_p*nu_p - E_p + 2*E_z*nu_pz*nu_pz),2) 
		   - (2*E_p*E_z*((2*((Em*NUm*(vf - 1)*(E_p*E_p*NUm + E_p*E_p + E_p*Em 
		   - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))/(2*(2*NUm - 1)*(NUm + 1)*(NUm + 1)
		   *(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) + (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))
		   /((4*E_z*nu_pz*nu_pz - 2*E_p + 2*E_p*nu_p)*(2*NUm*NUm + NUm - 1)))
		   *(2*NUm*NUm + NUm - 1))/(Em + 2*G_zp + Em*vf - 2*G_zp*vf 
		   - 4*G_zp*NUm*NUm - 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf 
		   + 4*G_zp*NUm*NUm*vf) - (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
		   + (E_p*E_z*nu_pz*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
		   *(2*E_z*nu_pz*nu_pz + E_p - E_p*nu_p)*(2*NUm*NUm + NUm - 1))
		   /((2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)*(E_p*E_p*NUm - E_p*E_p 
		   + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz)) 
		   - (8*E_z*Em*nu_pz*((Em*NUm)/(2*NUm*NUm + NUm - 1) 
		   - (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
		   *((2*((Em*NUm*(vf - 1)*(E_p*E_p*NUm + E_p*E_p + E_p*Em 
		   - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))/(2*(2*NUm - 1)*(NUm + 1)*(NUm + 1)
		   *(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) + (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))
		   /((4*E_z*nu_pz*nu_pz - 2*E_p + 2*E_p*nu_p)*(2*NUm*NUm + NUm - 1)))
		   *(2*NUm*NUm + NUm - 1))/(Em + 2*G_zp + Em*vf - 2*G_zp*vf 
		   - 4*G_zp*NUm*NUm - 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf 
		   + 4*G_zp*NUm*NUm*vf) - (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
		   + (E_p*E_z*nu_pz*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
		   *(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)*(2*NUm*NUm + NUm - 1))
		   /pow((E_p*E_p*NUm - E_p*E_p + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p 
		   - 2*E_z*Em*nu_pz*nu_pz),2) - (8*E_z*nu_pz*((Em*NUm)/(2*NUm*NUm + NUm - 1) 
		   - (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
		   *((2*((Em*NUm*(vf - 1)*(E_p*E_p*NUm + E_p*E_p + E_p*Em 
		   - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))/(2*(2*NUm - 1)*(NUm + 1)*(NUm + 1)
		   *(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
		   + (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*nu_pz*nu_pz - 2*E_p 
		   + 2*E_p*nu_p)*(2*NUm*NUm + NUm - 1)))*(2*NUm*NUm + NUm - 1))
		   /(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm 
		   - 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf) 
		   - (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) + (E_p*E_z*nu_pz*vf)
		   /(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))*(2*NUm*NUm + NUm - 1))
		   /(E_p*E_p*NUm - E_p*E_p + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p 
		   - 2*E_z*Em*nu_pz*nu_pz);
	  
	  PetscFunctionReturn(0);
    }
    /***************************************************************************
     *
     * With repect to Shear modulus in z-direction, G_zp
     * 
     **************************************************************************/	
	virtual PetscErrorCode D_r_ShearZP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp,
	                                  double Em, double NUm, double vf){
		
	  PetscFunctionBegin;
      
	  dkc = 0;
	  dmc = 0;
	  dpc = (4*Em*Em*vf)/pow((Em + 2*G_zp + Em*vf - 2*G_zp*vf + 2*G_zp*NUm - 2*G_zp*NUm*vf),2);
	  dlc = -(((Em*NUm*(1/((2*nu_p - 2)/E_p + (4*E_z*nu_pz*nu_pz)/(E_p*E_p)) 
	        - Em/(2*(NUm + 1)))*(vf - 1))/((2*NUm - 1)*(NUm + 1)) 
			+ (2*E_z*nu_pz*vf*(Em/(NUm + 1) 
			- (Em*NUm)/((2*NUm - 1)*(NUm + 1))))/(E_p*((2*nu_p - 2)/E_p 
			+ (4*E_z*nu_pz*nu_pz)/E_p*E_p)))*(vf - 1))/pow((vf*(Em/(NUm + 1) 
			- (Em*NUm)/((2*NUm - 1)*(NUm + 1))) - (G_zp + Em/(2*NUm + 2))*(vf - 1)),2);
	  dnc = (8*((Em*NUm*(vf - 1)*(E_p*E_p*NUm + E_p*E_p + E_p*Em - E_p*Em*nu_p 
	       - 2*E_z*Em*nu_pz*nu_pz))/(2*(2*NUm - 1)*(NUm + 1)*(NUm + 1)
		   *(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
		   + (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*nu_pz*nu_pz - 2*E_p 
		   + 2*E_p*nu_p)*(2*NUm*NUm + NUm - 1)))*((Em*NUm)/(2*NUm*NUm + NUm - 1) 
		   - (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))*(vf - 1)
		   *(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)*pow((2*NUm*NUm + NUm - 1),3))
		   /((E_p*E_p*NUm - E_p*E_p + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p 
		   - 2*E_z*Em*nu_pz*nu_pz)*pow((Em + 2*G_zp + Em*vf - 2*G_zp*vf 
		   - 4*G_zp*NUm*NUm - 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf 
		   + 4*G_zp*NUm*NUm*vf),2));
	  
	  PetscFunctionReturn(0);
    }

    /***************************************************************************
     *
     * With repect to Young's modulus in z-direction, E_z
     * 
     **************************************************************************/
	virtual PetscErrorCode D_r_Young(double nu_p, double nu_pz, double E_p, double E_z, double G_zp,
	                                  double Em, double NUm, double vf){
		
	  PetscFunctionBegin;
      
	  dkc = (E_p*E_p*NUm + E_p*E_p*vf + E_p*E_p + 2*E_p*Em 
	        - 2*E_p*E_p*NUm*NUm*vf - 2*E_p*Em*nu_p - 2*E_p*Em*vf 
			- 4*E_z*Em*nu_pz*nu_pz - E_p*E_p*NUm*vf + 4*E_z*Em*nu_pz*nu_pz*vf 
			+ 2*E_p*Em*nu_p*vf)/(4*(2*NUm - 1)*(NUm + 1)*(NUm + 1)
			*(((vf - 1)*(E_p*E_p*NUm + E_p*E_p + E_p*Em - E_p*Em*nu_p 
			- 2*E_z*Em*nu_pz*nu_pz))/(2*(NUm + 1)*(2*E_z*nu_pz*nu_pz - E_p 
			+ E_p*nu_p)) + (Em*vf*(NUm - 1))/(2*NUm*NUm + NUm - 1))
			*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			+ (Em*(vf - 2*NUm + 1)*(E_p*E_p*NUm + E_p*E_p*vf + E_p*E_p 
			+ E_p*Em - 2*E_p*E_p*NUm*NUm*vf - E_p*Em*nu_p - E_p*Em*vf 
			- 2*E_z*Em*nu_pz*nu_pz - E_p*E_p*NUm*vf + 2*E_z*Em*nu_pz*nu_pz*vf 
			+ E_p*Em*nu_p*vf))/(4*(2*NUm - 1)*(NUm + 1)*(NUm + 1)*pow((((vf - 1)
			*(E_p*E_p*NUm + E_p*E_p + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
			/(2*(NUm + 1)*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			+ (Em*vf*(NUm - 1))/(2*NUm*NUm + NUm - 1)),2)*(4*NUm*NUm + 2*NUm - 2)
			*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p));
     dmc = -((vf - 1)*(2*Em*Em*nu_p - 6*E_p*E_p*NUm + 3*E_p*E_p*vf + 3*Em*Em*vf 
		    + 9*E_p*E_p + Em*Em - 23*E_p*E_p*NUm*NUm + 8*E_p*E_p*pow(NUm,3) 
			+ 16*E_p*E_p*pow(NUm,4) + pow(Em*nu_p,2) + 6*E_p*Em - 5*pow(E_p*NUm,2)*vf 
			- 4*E_p*E_p*pow(NUm,3)*vf - 2*E_p*Em*NUm + 3*pow(Em*nu_p,2)*vf 
			+ 6*E_p*Em*nu_p - 6*E_p*Em*vf - 8*E_p*Em*NUm*NUm + 2*pow(E_p*NUm,2)*vf 
			- 4*Em*Em*NUm*vf + 6*Em*Em*nu_p*vf - 8*E_p*Em*NUm*NUm*nu_p 
			+ 8*E_p*Em*NUm*NUm*vf - 8*Em*Em*NUm*nu_p*vf - 4*NUm*pow(Em*nu_p,2)*vf 
			- 2*E_p*Em*NUm*nu_p + 2*E_p*Em*NUm*vf - 6*E_p*Em*nu_p*vf 
			+ 2*E_p*Em*NUm*nu_p*vf + 8*E_p*Em*NUm*NUm*nu_p*vf))
			/(2*(NUm + 1)*pow((3*E_p + Em + Em*nu_p - 3*E_p*vf + 3*Em*vf 
			- 4*E_p*NUm*NUm - E_p*NUm + E_p*NUm*vf - 4*Em*NUm*vf 
			+ 3*Em*nu_p*vf + 4*E_p*NUm*NUm*vf - 4*Em*NUm*nu_p*vf),2));
	  dpc = ((Em/(2*pow((NUm + 1),2)) + G_zp/(2*(NUm + 1)))*(vf - 1) 
		  - (G_zp*vf)/(NUm + 1))/((G_zp + Em/(2*(NUm + 1)))*(vf - 1) 
		  - (Em*vf)/(NUm + 1)) - (((Em*Em/(4*(NUm + 1)*(NUm + 1)) 
		  + (Em*G_zp)/(2*(NUm + 1)))*(vf - 1) - (Em*G_zp*vf)/(NUm + 1))
		  *((vf - 1)/(2*(NUm + 1)) - vf/(NUm + 1)))/pow(((G_zp + Em/(2*NUm + 2))*(vf - 1) - (Em*vf)/(NUm + 1)),2);
	  dlc = - ((NUm*(1/((2*nu_p - 2)/E_p + (4*E_z*nu_pz*nu_pz)/(E_p*E_p)) 
	        - Em/(2*(NUm + 1)))*(vf - 1))/((2*NUm - 1)*(NUm + 1)) 
			- (Em*NUm*(vf - 1))/(2*(2*NUm - 1)*(NUm + 1)*(NUm + 1)) 
			+ (2*E_z*nu_pz*vf*(1/(NUm + 1) - NUm/((2*NUm - 1)*(NUm + 1))))
			/(E_p*((2*nu_p - 2)/E_p + (4*E_z*nu_pz*nu_pz)/(E_p*E_p))))/(vf*(Em/(NUm + 1) 
			- (Em*NUm)/((2*NUm - 1)*(NUm + 1))) - (G_zp + Em/(2*(NUm + 1)))*(vf - 1)) 
			- (((vf - 1)/(2*(NUm + 1)) - vf*(1/(NUm + 1) 
			- NUm/((2*NUm - 1)*(NUm + 1))))*((Em*NUm*(1/((2*nu_p - 2)/E_p 
			+ (4*E_z*nu_pz*nu_pz)/(E_p*E_p)) - Em/(2*(NUm + 1)))*(vf - 1))
			/((2*NUm - 1)*(NUm + 1)) + (2*E_z*nu_pz*vf*(Em/(NUm + 1) 
			- (Em*NUm)/((2*NUm - 1)*(NUm + 1))))/(E_p*((2*nu_p - 2)/E_p 
			+ (4*E_z*nu_pz*nu_pz)/(E_p*E_p)))))/pow((vf*(Em/(NUm + 1) 
			- (Em*NUm)/((2*NUm - 1)*(NUm + 1))) 
			- (G_zp + Em/(2*NUm + 2))*(vf - 1)),2);
      dnc = (2*((Em*NUm)/(2*NUm*NUm + NUm - 1) - (E_p*E_z*nu_pz)
	        /(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))*((NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
			- (2*((NUm*(vf - 1)*(E_p*E_p*NUm + E_p*E_p + E_p*Em - E_p*Em*nu_p 
			- 2*E_z*Em*nu_pz*nu_pz))/(2*(2*NUm - 1)*(NUm + 1)*(NUm + 1)*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			- (Em*NUm*(vf - 1))/(2*(2*NUm - 1)*(NUm + 1)*(NUm + 1)) 
			+ (2*E_p*E_z*nu_pz*vf*(NUm - 1))/((4*E_z*nu_pz*nu_pz - 2*E_p + 2*E_p*nu_p)
			*(2*NUm*NUm + NUm - 1)))*(2*NUm*NUm + NUm - 1))/(Em + 2*G_zp + Em*vf 
			- 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf 
			+ 4*G_zp*NUm*NUm*vf) + (2*((Em*NUm*(vf - 1)*(E_p*E_p*NUm + E_p*E_p + E_p*Em 
			- E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))/(2*(2*NUm - 1)*(NUm + 1)*(NUm + 1)
			*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			+ (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*nu_pz*nu_pz - 2*E_p + 2*E_p*nu_p)
			*(2*NUm*NUm + NUm - 1)))*(vf - 2*NUm + 1)*(2*NUm*NUm + NUm - 1))
			/pow((Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm 
			- 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf),2))
			*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)*(2*NUm*NUm + NUm - 1))
			/(E_p*E_p*NUm - E_p*E_p + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p 
			- 2*E_z*Em*nu_pz*nu_pz) - (2*NUm*((2*((Em*NUm*(vf - 1)*(E_p*E_p*NUm 
			+ E_p*E_p + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
			/(2*(2*NUm - 1)*(NUm + 1)*(NUm + 1)*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			+ (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*nu_pz*nu_pz - 2*E_p 
			+ 2*E_p*nu_p)*(2*NUm*NUm + NUm - 1)))*(2*NUm*NUm + NUm - 1))
			/(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm 
			- 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf) 
			- (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
			+ (E_p*E_z*nu_pz*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
			*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))/(E_p*E_p*NUm - E_p*E_p 
			+ 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz) 
			- ((NUm - 1)*(vf - 1))/(2*NUm*NUm + NUm - 1) 
			- (2*((Em*NUm)/(2*NUm*NUm + NUm - 1) 
			- (E_p*E_z*nu_pz)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
			*((2*((Em*NUm*(vf - 1)*(E_p*E_p*NUm + E_p*E_p + E_p*Em 
			- E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))/(2*(2*NUm - 1)*(NUm + 1)*(NUm + 1)
			*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			+ (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))
			/((4*E_z*nu_pz*nu_pz - 2*E_p + 2*E_p*nu_p)
			*(2*NUm*NUm + NUm - 1)))*(2*NUm*NUm + NUm - 1))
			/(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*NUm*NUm - 2*Em*NUm 
			- 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*NUm*NUm*vf) 
			- (Em*NUm*(vf - 1))/(2*NUm*NUm + NUm - 1) 
			+ (E_p*E_z*nu_pz*vf)/(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p))
			*pow((2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p),2)*(2*NUm*NUm + NUm - 1))
			/pow((E_p*E_p*NUm - E_p*E_p + 2*E_p*E_p*NUm*NUm + E_p*Em - E_p*Em*nu_p 
			- 2*E_z*Em*nu_pz*nu_pz),2);

	  PetscFunctionReturn(0);
    }	 
    /***************************************************************************
     *
     * With repect to Young's modulus in z-direction, E_z
     * 
     **************************************************************************/
	virtual PetscErrorCode D_r_Poisson(double nu_p, double nu_pz, double E_p, double E_z, double G_zp,
	                                  double Em, double NUm, double vf){
		
	  PetscFunctionBegin;
      
	  dkc = - (Em*(5*E_p*E_p*NUm - E_p*E_p*vf + E_p*E_p + 4*E_p*E_p*NUm*NUm 
	        - 4*E_p*E_p*pow(NUm,3)*vf + 6*E_p*Em*NUm + 3*E_p*E_p*NUm*vf 
			- 12*E_z*Em*NUm*nu_pz*nu_pz - 6*E_p*Em*NUm*nu_p - 6*E_p*Em*NUm*vf 
			+ 6*E_p*Em*NUm*nu_p*vf + 12*E_z*Em*NUm*nu_pz*nu_pz*vf))
			/(4*(2*NUm - 1)*(2*NUm - 1)*pow((NUm + 1),3)*(((vf - 1)*(E_p*E_p*NUm 
			+ E_p*E_p + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
			/(2*(NUm + 1)*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			+ (Em*vf*(NUm - 1))/(2*NUm*NUm + NUm - 1))
			*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
			- (Em*Em*(4*NUm + vf + 4*NUm*vf - 4*NUm*NUm - 1)
			*(E_p*E_p*NUm + E_p*E_p*vf + E_p*E_p + E_p*Em 
			- 2*E_p*E_p*NUm*NUm*vf - E_p*Em*nu_p - E_p*Em*vf 
			- 2*E_z*Em*nu_pz*nu_pz - E_p*E_p*NUm*vf + 2*E_z*Em*nu_pz*nu_pz*vf 
			+ E_p*Em*nu_p*vf))/(8*(2*NUm - 1)*(NUm + 1)*(NUm + 1)*pow((((vf - 1)
			*(E_p*E_p*NUm + E_p*E_p + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
		   /(2*(NUm + 1)*(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)) 
		   + (Em*vf*(NUm - 1))/(2*NUm*NUm + NUm - 1)),2)
		   *(2*E_z*nu_pz*nu_pz - E_p + E_p*nu_p)*pow((2*NUm*NUm + NUm - 1),2));
     dmc = -(Em*(vf - 1)*(6*E_p*E_p*NUm - 2*Em*Em*nu_p + E_p*E_p*vf 
		    + Em*Em*vf - 9*E_p*E_p - Em*Em + 23*E_p*E_p*NUm*NUm - 8*E_p*E_p*pow(NUm,3) 
			- 16*E_p*E_p*pow(NUm,4) - Em*Em*nu_p*nu_p - 6*E_p*Em + 17*E_p*E_p*NUm*NUm*vf 
			+ 8*E_p*E_p*pow(NUm,3)*vf + 2*E_p*Em*NUm + Em*Em*nu_p*nu_p*vf 
			- 6*E_p*Em*nu_p - 2*E_p*Em*vf + 8*E_p*Em*NUm*NUm 
			+ 10*E_p*E_p*NUm*vf + 8*Em*Em*NUm*vf + 2*Em*Em*nu_p*vf 
			+ 8*E_p*Em*NUm*NUm*nu_p - 16*E_p*Em*NUm*NUm*vf 
			+ 16*Em*Em*NUm*nu_p*vf + 8*Em*Em*NUm*nu_p*nu_p*vf 
			+ 2*E_p*Em*NUm*nu_p - 18*E_p*Em*NUm*vf - 2*E_p*Em*nu_p*vf 
			- 18*E_p*Em*NUm*nu_p*vf - 16*E_p*Em*NUm*NUm*nu_p*vf))
			/(2*(NUm + 1)*(NUm + 1)*pow((3*E_p + Em + Em*nu_p - 3*E_p*vf 
			+ 3*Em*vf - 4*E_p*NUm*NUm - E_p*NUm + E_p*NUm*vf 
			- 4*Em*NUm*vf + 3*Em*nu_p*vf + 4*E_p*NUm*NUm*vf 
			- 4*Em*NUm*nu_p*vf),2));
      dpc = - ((Em*Em/(2*pow((NUm + 1),3)) + (Em*G_zp)/(2*(NUm + 1)*(NUm + 1)))*(vf - 1) 
	        - (Em*G_zp*vf)/pow(NUm + 1,2))/((G_zp + Em/(2*(NUm + 1)))*(vf - 1) 
			- (Em*vf)/(NUm + 1)) - (((Em*vf)/pow(NUm + 1,2) 
			- (Em*(vf - 1))/(2*(NUm + 1)*(NUm + 1)))*((Em*Em/(4*(NUm + 1)*(NUm + 1)) 
			+ (Em*G_zp)/(2*(NUm + 1)))*(vf - 1) 
			- (Em*G_zp*vf)/(NUm + 1)))/pow(((G_zp + Em/(2*NUm + 2))*(vf - 1) 
			- (Em*vf)/(NUm + 1)),2);
	  dlc = ((2*E_z*nu_pz*vf*(Em/pow((NUm + 1),2) + Em/((2*NUm - 1)*(NUm + 1)) 
	        - (Em*NUm)/((2*NUm - 1)*(NUm + 1)*(NUm + 1)) 
			 - (2*Em*NUm)/((2*NUm - 1)*(2*NUm - 1)*(NUm + 1))))/(E_p*((2*nu_p - 2)/E_p 
			 + (4*E_z*nu_pz*nu_pz)/pow(E_p,2))) - (Em*Em*NUm*(vf - 1))
			 /(2*(2*NUm - 1)*pow(NUm + 1,3)) - (Em*(1/((2*nu_p - 2)/E_p 
			 + (4*E_z*nu_pz*nu_pz)/pow(E_p,2)) - Em/(2*(NUm + 1)))*(vf - 1))
			 /((2*NUm - 1)*(NUm + 1)) + (Em*NUm*(1/((2*nu_p - 2)/E_p 
			 + (4*E_z*nu_pz*nu_pz)/pow(E_p,2)) - Em/(2*(NUm + 1)))*(vf - 1))
			 /((2*NUm - 1)*pow(NUm + 1,2)) + (2*Em*NUm*(1/((2*nu_p - 2)/E_p 
			 + (4*E_z*nu_pz*nu_pz)/pow(E_p,2)) - Em/(2*(NUm + 1)))*(vf - 1))
			 /((2*NUm - 1)*(2*NUm - 1)*(NUm + 1)))/(vf*(Em/(NUm + 1) 
			 - (Em*NUm)/((2*NUm - 1)*(NUm + 1))) 
			 - (G_zp + Em/(2*(NUm + 1)))*(vf - 1)) 
			 - (((Em*NUm*(1/((2*nu_p - 2)/E_p + (4*E_z*nu_pz*nu_pz)/pow(E_p,2)) 
			 - Em/(2*(NUm + 1)))*(vf - 1))/((2*NUm - 1)*(NUm + 1)) 
			 + (2*E_z*nu_pz*vf*(Em/(NUm + 1) 
			 - (Em*NUm)/((2*NUm - 1)*(NUm + 1))))/(E_p*((2*nu_p - 2)/E_p 
			 + (4*E_z*nu_pz*nu_pz)/pow(E_p,2))))*(vf*(Em/pow(NUm + 1,2) 
			 + Em/((2*NUm - 1)*(NUm + 1)) - (Em*NUm)/((2*NUm - 1)*pow(NUm + 1,2)) 
			 - (2*Em*NUm)/((2*NUm - 1)*(2*NUm - 1)*(NUm + 1))) 
			 - (Em*(vf - 1))/(2*pow(NUm + 1,2))))/pow((vf*(Em/(NUm + 1) 
			 - (Em*NUm)/((2*NUm - 1)*(NUm + 1))) 
			 - (G_zp + Em/(2*NUm + 2))*(vf - 1)),2);		
		dnc = (Em*(4*NUm + 1)*(NUm - 1)*(vf - 1))/pow((NUm + 2*pow(NUm,2) - 1),2) 
		     - (2*((Em*NUm)/(2*pow(NUm,2) + NUm - 1) 
			 - (E_p*E_z*nu_pz)/(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p))
			 *(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p)*(2*pow(NUm,2) + NUm - 1)
			 *((2*(4*NUm + 1)*((Em*NUm*(vf - 1)*(pow(E_p,2)*NUm + pow(E_p,2) 
			 + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*nu_pz*nu_pz))
			 /(2*(2*NUm - 1)*pow(NUm+1,2)*(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p)) 
			 + (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*pow(nu_pz,2)- 2*E_p + 2*E_p*nu_p)
			 *(2*pow(NUm,2) + NUm - 1))))/(Em + 2*G_zp + Em*vf - 2*G_zp*vf 
			 - 4*G_zp*pow(NUm,2) - 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf 
			 + 4*G_zp*pow(NUm,2)*vf) - (Em*(vf - 1))/(2*pow(NUm,2) + NUm - 1) 
			 + (2*((Em*NUm*(vf - 1)*(pow(E_p,2)*NUm + pow(E_p,2) + E_p*Em 
			 - E_p*Em*nu_p - 2*E_z*Em*pow(nu_pz,2)))/(2*(2*NUm - 1)*pow(NUm+1,2)
			 *(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p)) 
			 + (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*pow(nu_pz,2)
			 - 2*E_p + 2*E_p*nu_p)*(2*pow(NUm,2) + NUm - 1)))
			 *(2*pow(NUm,2) + NUm - 1)*(2*Em + 2*G_zp - 2*G_zp*vf 
			 + 8*G_zp*NUm - 8*G_zp*NUm*vf))/pow((Em + 2*G_zp + Em*vf 
			 - 2*G_zp*vf - 4*G_zp*pow(NUm,2) - 2*Em*NUm - 2*G_zp*NUm 
			 + 2*G_zp*NUm*vf + 4*G_zp*pow(NUm,2)*vf),2)
			 + (Em*NUm*(4*NUm + 1)*(vf - 1))/pow((NUm + 2*pow(NUm,2) - 1),2) 
			 + (Em*(2*pow(NUm,2) + NUm - 1)*(pow(E_p,2)*NUm - pow(E_p,2)*vf 
			 + pow(E_p,2) + 2*pow(E_p,2)*pow(NUm,2) + 2*pow(E_p,2)*pow(NUm,3) + E_p*Em 
			 - 2*pow(E_p,2)*pow(NUm,2)*vf - 2*pow(E_p,2)*pow(NUm,3)*vf - E_p*Em*NUm 
			 - E_p*Em*nu_p - E_p*Em*vf + 4*E_p*Em*pow(NUm,2) 
			 - 2*E_z*Em*pow(nu_pz,2)- pow(E_p,2)*NUm*vf - 4*E_p*Em*pow(NUm,2)*nu_p 
			 + 2*E_z*Em*NUm*pow(nu_pz,2)- 4*E_p*Em*pow(NUm,2)*vf 
			 + 2*E_z*Em*pow(nu_pz,2)*vf - 8*E_z*Em*pow(NUm,2)*pow(nu_pz,2)
			 + E_p*Em*NUm*nu_p + E_p*Em*NUm*vf + E_p*Em*nu_p*vf 
			 + 8*E_p*E_z*NUm*nu_pz*vf - E_p*Em*NUm*nu_p*vf 
			 + 4*E_p*E_z*pow(NUm,2)*nu_pz*vf - 4*E_p*E_z*pow(NUm,3)*nu_pz*vf 
			 + 4*E_p*Em*pow(NUm,2)*nu_p*vf - 2*E_z*Em*NUm*pow(nu_pz,2)*vf 
			 + 8*E_z*Em*pow(NUm,2)*pow(nu_pz,2)*vf))/(pow((2*NUm - 1),2)
			 *pow(NUm + 1,3)*(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p)
			 *(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*pow(NUm,2) 
			 - 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*pow(NUm,2)*vf))))
			 /(pow(E_p,2)*NUm - pow(E_p,2) + 2*pow(E_p,2)*pow(NUm,2) + E_p*Em 
			 - E_p*Em*nu_p - 2*E_z*Em*pow(nu_pz,2)) 
			 - (Em*(vf - 1))/(2*pow(NUm,2) + NUm - 1) 
			 - (2*(4*NUm + 1)*((Em*NUm)/(2*pow(NUm,2) + NUm - 1) 
			 - (E_p*E_z*nu_pz)/(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p))
			 *((2*((Em*NUm*(vf - 1)*(pow(E_p,2)*NUm + pow(E_p,2) + E_p*Em 
			 - E_p*Em*nu_p - 2*E_z*Em*pow(nu_pz,2)))/(2*(2*NUm - 1)
			 *pow(NUm+1,2)*(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p)) 
			 + (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*pow(nu_pz,2)
			 - 2*E_p + 2*E_p*nu_p)*(2*pow(NUm,2) + NUm - 1)))
			 *(2*pow(NUm,2) + NUm - 1))/(Em + 2*G_zp + Em*vf 
			 - 2*G_zp*vf - 4*G_zp*pow(NUm,2) - 2*Em*NUm - 2*G_zp*NUm 
			 + 2*G_zp*NUm*vf + 4*G_zp*pow(NUm,2)*vf) 
			 - (Em*NUm*(vf - 1))/(2*pow(NUm,2) + NUm - 1) 
			 + (E_p*E_z*nu_pz*vf)/(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p))
			 *(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p))/(pow(E_p,2)*NUm - pow(E_p,2) 
			 + 2*pow(E_p,2)*pow(NUm,2) + E_p*Em - E_p*Em*nu_p 
			 - 2*E_z*Em*pow(nu_pz,2)) + (2*Em*(2*pow(NUm,2) + 1)
			 *((2*((Em*NUm*(vf - 1)*(pow(E_p,2)*NUm + pow(E_p,2) 
			 + E_p*Em - E_p*Em*nu_p - 2*E_z*Em*pow(nu_pz,2)))
			 /(2*(2*NUm - 1)*pow(NUm+1,2)*(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p)) 
			 + (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*pow(nu_pz,2)- 2*E_p + 2*E_p*nu_p)
			 *(2*pow(NUm,2) + NUm - 1)))*(2*pow(NUm,2) + NUm - 1))
			 /(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*pow(NUm,2) 
			 - 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*pow(NUm,2)*vf) 
			 - (Em*NUm*(vf - 1))/(2*pow(NUm,2) + NUm - 1) 
			 + (E_p*E_z*nu_pz*vf)/(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p))
			 *(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p))/((2*pow(NUm,2) + NUm - 1)
			 *(pow(E_p,2)*NUm - pow(E_p,2) + 2*pow(E_p,2)*pow(NUm,2) + E_p*Em 
			 - E_p*Em*nu_p - 2*E_z*Em*pow(nu_pz,2))) 
			 + (2*pow(E_p,2)*(4*NUm + 1)*((Em*NUm)/(2*pow(NUm,2) + NUm - 1) 
			 - (E_p*E_z*nu_pz)/(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p))
			 *((2*((Em*NUm*(vf - 1)*(pow(E_p,2)*NUm + pow(E_p,2) + E_p*Em 
			 - E_p*Em*nu_p - 2*E_z*Em*pow(nu_pz,2)))/(2*(2*NUm - 1)
			 *pow(NUm+1,2)*(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p)) 
			 + (2*E_p*E_z*Em*nu_pz*vf*(NUm - 1))/((4*E_z*pow(nu_pz,2)
			 - 2*E_p + 2*E_p*nu_p)*(2*pow(NUm,2) + NUm - 1)))*(2*pow(NUm,2) + NUm - 1))
			 /(Em + 2*G_zp + Em*vf - 2*G_zp*vf - 4*G_zp*pow(NUm,2) 
			 - 2*Em*NUm - 2*G_zp*NUm + 2*G_zp*NUm*vf + 4*G_zp*pow(NUm,2)*vf) 
			 - (Em*NUm*(vf - 1))/(2*pow(NUm,2) + NUm - 1) 
			 + (E_p*E_z*nu_pz*vf)/(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p))
			 *(2*E_z*pow(nu_pz,2)- E_p + E_p*nu_p)*(2*pow(NUm,2) + NUm - 1))
			 /pow((pow(E_p,2)*NUm - pow(E_p,2) + 2*pow(E_p,2)*pow(NUm,2) + E_p*Em 
			 - E_p*Em*nu_p - 2*E_z*Em*pow(nu_pz,2)),2);
	 
	  PetscFunctionReturn(0);
    }		
  };
}
#endif //__D_R_ELASTICFEMETHODTRANSISO_HPP__

