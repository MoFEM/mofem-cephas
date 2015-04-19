/* Copyright (C) 2014, 
 *   Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
 *   Xiao-Yi Zhou (xiaoyi.zhou AT newcastle.ac.uk)
 * --------------------------------------------------------------
 * This routine calculates the second-order partial derivative of constitutive 
 * matrix, D_rs, with respect to Young's modulus in z-direction, E_z, Young's
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

#ifndef __MATERIALCONSTITUTIVEMATRIX_SECONDORDERDERIVATIVE_HPP__
#define __MATERIALCONSTITUTIVEMATRIX_SECONDORDERDERIVATIVE_HPP__


#include <boost/numeric/ublas/symmetric.hpp>

namespace MoFEM {
  
// =============================================================================
//
//  TRANSVERSELY ISOTROPIC MATERIAL
//
// =============================================================================
  struct TransverseIsotropicStiffnessMatrix_SecondOrderDerivative {
    double nu_p, nu_pz, E_p, E_z, G_zp;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsPoissonP;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsPoissonPPoissonPZ;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsPoissonPYoungP;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsPoissonPYoungZ;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsPoissonPShearZP;
	
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsPoissonPZ;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsPoissonPZYoungP;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsPoissonPZYoungZ;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsPoissonPZShearZP;
	
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsYoungP;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsYoungPYoungZ;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsYoungPShearZP;
	
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsYoungZ;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsYoungZShearZP;
	
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsShearZP;
	
    /***************************************************************************
     *
     * With repect to Poisson's ration in p-direction
     * 
     **************************************************************************/
      virtual PetscErrorCode D_rs_PoissonP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from D_rs_ElasticFEMethodTransIso"<<endl;
      StiffnessMatrix_rsPoissonP.resize(6);
      StiffnessMatrix_rsPoissonP.clear();
  
      StiffnessMatrix_rsPoissonP(0,0) = -(2*E_p*(-E_z*nu_pz*nu_pz+E_p)*(3*E_p*E_p*nu_p*nu_p+E_p*E_p+6*E_p*E_z*nu_p*nu_pz*nu_pz-2*E_p*E_z*nu_pz*nu_pz+4*E_z*E_z*pow(nu_pz,4)))/(pow((nu_p+1),3)*pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3));
      StiffnessMatrix_rsPoissonP(0,1) = -(2*E_p*(pow(E_p,3)*pow(nu_p,3)+3*pow(E_p,3)*nu_p+3*E_p*E_p*E_z*nu_p*nu_p*nu_pz*nu_pz-6*E_p*E_p*E_z*nu_p*nu_pz*nu_pz+3*E_p*E_p*E_z*nu_pz*nu_pz+6*E_p*E_z*E_z*nu_p*pow(nu_pz,4)-6*E_p*E_z*E_z*pow(nu_pz,4)+4*pow(E_z,3)*pow(nu_pz,6)))/(pow((nu_p+1),3)*pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3));
      StiffnessMatrix_rsPoissonP(0,2) = -(2*pow(E_p,3)*E_z*nu_pz)/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);

      //StiffnessMatrix_rsPoissonP(1,0) = StiffnessMatrix_rsPoissonP(0,1);
      StiffnessMatrix_rsPoissonP(1,1) = StiffnessMatrix_rsPoissonP(0,0);
      StiffnessMatrix_rsPoissonP(1,2) = StiffnessMatrix_rsPoissonP(0,2);

      //StiffnessMatrix_rsPoissonP(2,0) = StiffnessMatrix_rsPoissonP(0,2);
      //StiffnessMatrix_rsPoissonP(2,1) = StiffnessMatrix_rsPoissonP(0,2);
      StiffnessMatrix_rsPoissonP(2,2) = -(4*E_p*E_p*E_z*E_z*nu_pz*nu_pz)/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);

      StiffnessMatrix_rsPoissonP(3,3) = E_p/pow((nu_p + 1),3);

      StiffnessMatrix_rsPoissonP(4,4) = 0;

      StiffnessMatrix_rsPoissonP(5,5) = 0;
      PetscFunctionReturn(0);
    }
	
   /**************************************************************************************************
     *
     * With repect to Poisson's ration in p-direction, nu_p, and Poisson's ratio in z-direction, nu_pz
     * 
     ************************************************************************************************/
      virtual PetscErrorCode D_rs_PoissonPPoissonPZ(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from D_rs_ElasticFEMethodTransIso"<<endl;
      StiffnessMatrix_rsPoissonPPoissonPZ.resize(6);
      StiffnessMatrix_rsPoissonPPoissonPZ.clear();
  
      StiffnessMatrix_rsPoissonPPoissonPZ(0,0) = -(4*pow(E_p,3)*E_z*nu_pz)/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);
      StiffnessMatrix_rsPoissonPPoissonPZ(0,1) = StiffnessMatrix_rsPoissonPPoissonPZ(0,0);
      StiffnessMatrix_rsPoissonPPoissonPZ(0,2) = -(pow(E_p,2)*E_z*(6*E_z*pow(nu_pz,2) + E_p - E_p*nu_p))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);

      //StiffnessMatrix_rsPoissonPPoissonPZ(1,0) = StiffnessMatrix_rsPoissonPPoissonPZ(0,1);
      StiffnessMatrix_rsPoissonPPoissonPZ(1,1) = StiffnessMatrix_rsPoissonPPoissonPZ(0,0);
      StiffnessMatrix_rsPoissonPPoissonPZ(1,2) = StiffnessMatrix_rsPoissonPPoissonPZ(0,2);

      //StiffnessMatrix_rsPoissonPPoissonPZ(2,0) = StiffnessMatrix_rsPoissonPPoissonPZ(0,2);
      //StiffnessMatrix_rsPoissonPPoissonPZ(2,1) = StiffnessMatrix_rsPoissonPPoissonPZ(0,2);
      StiffnessMatrix_rsPoissonPPoissonPZ(2,2) = -(8*E_p*pow(E_z,3)*pow(nu_pz,3) - 4*pow(E_p,2)*pow(E_z,2)*nu_pz*(nu_p - 1))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);

      StiffnessMatrix_rsPoissonPPoissonPZ(3,3) = 0;

      StiffnessMatrix_rsPoissonPPoissonPZ(4,4) = 0;

      StiffnessMatrix_rsPoissonPPoissonPZ(5,5) = 0;
      PetscFunctionReturn(0);
    }	
	
    /**************************************************************************************************
     *
     * With repect to Poisson's ration in p-direction, nu_p, and Young's modulus in p-direction, E_p
     * 
     ************************************************************************************************/
      virtual PetscErrorCode D_rs_PoissonPYoungP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from D_rs_ElasticFEMethodTransIso"<<endl;
      StiffnessMatrix_rsPoissonPYoungP.resize(6);
      StiffnessMatrix_rsPoissonPYoungP.clear();
  
      StiffnessMatrix_rsPoissonPYoungP(0,0) = -(nu_p*(2*pow(E_p,3) - 12*pow(E_p*nu_pz,2)*E_z + 6*E_p*pow(E_z,2)*pow(nu_pz,4)) - 2*pow(E_p,3)*pow(nu_p,2) + 4*pow(E_z,3)*pow(nu_pz,6) - 6*E_p*pow(E_z,2)*pow(nu_pz,4))/(pow((nu_p + 1),2)*pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3));
      StiffnessMatrix_rsPoissonPYoungP(0,1) =  (nu_p*(pow(E_p,3) + 6*E_p*pow(E_z,2)*pow(nu_pz,4)) - pow(E_p,3) - pow(nu_p,2)*(pow(E_p,3) - 6*pow(E_p,2)*E_z*pow(nu_pz,2)) + pow(E_p*nu_p,3) + 4*pow(E_z,3)*pow(nu_pz,6) + 6*E_z*pow(E_p*nu_pz,2) - 6*E_p*pow(E_z,2)*pow(nu_pz,4))/(pow((nu_p + 1),2)*pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3));
      StiffnessMatrix_rsPoissonPYoungP(0,2) =  (4*E_p*pow(E_z,2)*pow(nu_pz,3))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);

      //StiffnessMatrix_rsPoissonPYoungP(1,0) = StiffnessMatrix_rsPoissonPYoungP(0,1);
      StiffnessMatrix_rsPoissonPYoungP(1,1) = StiffnessMatrix_rsPoissonPYoungP(0,0);
      StiffnessMatrix_rsPoissonPYoungP(1,2) = StiffnessMatrix_rsPoissonPYoungP(0,2);

      //StiffnessMatrix_rsPoissonPYoungP(2,0) = StiffnessMatrix_rsPoissonPYoungP(0,2);
      //StiffnessMatrix_rsPoissonPYoungP(2,1) = StiffnessMatrix_rsPoissonPYoungP(0,2);
      StiffnessMatrix_rsPoissonPYoungP(2,2) = (4*pow(E_z,3)*pow(nu_pz,4) - 2*E_p*pow(E_z*nu_pz,2)*(nu_p - 1))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);

      StiffnessMatrix_rsPoissonPYoungP(3,3) = -1/(2*pow((nu_p + 1),2));

      StiffnessMatrix_rsPoissonPYoungP(4,4) = 0;

      StiffnessMatrix_rsPoissonPYoungP(5,5) = 0;
      PetscFunctionReturn(0);
    }
	
    /**************************************************************************************************
     *
     * With repect to Poisson's ration in p-direction, nu_p, and Young's modulus in z-direction, E_z
     * 
     ************************************************************************************************/
      virtual PetscErrorCode D_rs_PoissonPYoungZ(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from D_rs_ElasticFEMethodTransIso"<<endl;
      StiffnessMatrix_rsPoissonPYoungZ.resize(6);
      StiffnessMatrix_rsPoissonPYoungZ.clear();
  
      StiffnessMatrix_rsPoissonPYoungZ(0,0) = -(2*pow(E_p,3)*pow(nu_pz,2))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);
      StiffnessMatrix_rsPoissonPYoungZ(0,1) =  StiffnessMatrix_rsPoissonPYoungZ(0,0);
      StiffnessMatrix_rsPoissonPYoungZ(0,2) =  -(pow(E_p,2)*nu_pz*(2*E_z*pow(nu_pz,2) + E_p - E_p*nu_p))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);

      //StiffnessMatrix_rsPoissonPYoungZ(1,0) = StiffnessMatrix_rsPoissonPYoungZ(0,1);
      StiffnessMatrix_rsPoissonPYoungZ(1,1) = StiffnessMatrix_rsPoissonPYoungZ(0,0);
      StiffnessMatrix_rsPoissonPYoungZ(1,2) = StiffnessMatrix_rsPoissonPYoungZ(0,2);

      //StiffnessMatrix_rsPoissonPYoungZ(2,0) = StiffnessMatrix_rsPoissonPYoungZ(0,2);
      //StiffnessMatrix_rsPoissonPYoungZ(2,1) = StiffnessMatrix_rsPoissonPYoungZ(0,2);
      StiffnessMatrix_rsPoissonPYoungZ(2,2) = (4*pow(E_p*nu_pz,2)*E_z*(nu_p - 1))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);

      StiffnessMatrix_rsPoissonPYoungZ(3,3) = 0;

      StiffnessMatrix_rsPoissonPYoungZ(4,4) = 0;

      StiffnessMatrix_rsPoissonPYoungZ(5,5) = 0;
      PetscFunctionReturn(0);
    }
	
    /**************************************************************************************************
     *
     * With repect to Poisson's ration in p-direction, nu_p, and Shear modulus in z-direction, G_zp
     * 
     ************************************************************************************************/
      virtual PetscErrorCode D_rs_PoissonPShearZP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from D_rs_ElasticFEMethodTransIso"<<endl;
      StiffnessMatrix_rsPoissonPShearZP.resize(6);
      StiffnessMatrix_rsPoissonPShearZP.clear();
      PetscFunctionReturn(0);
    }

    /***************************************************************************
     *
     * With repect to Poisson's ration in z-direction, nu_pz
     * 
     **************************************************************************/
      virtual PetscErrorCode D_rs_PoissonPZ(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from TransverseIsotropicStiffnessMatrix_rsPoissonPZ "<<endl;
      StiffnessMatrix_rsPoissonPZ.resize(6);
      StiffnessMatrix_rsPoissonPZ.clear();
      
      StiffnessMatrix_rsPoissonPZ(0,0) = -(2*E_p*E_p*E_z*(6*E_z*nu_pz*nu_pz+E_p-E_p*nu_p))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);
      StiffnessMatrix_rsPoissonPZ(0,1) = StiffnessMatrix_rsPoissonPZ(0,0);
      StiffnessMatrix_rsPoissonPZ(0,2) = -(8*E_p*pow(E_z,3)*pow(nu_pz,3)-4*pow(E_p,2)*pow(E_z,2)*nu_pz*(3*nu_p-3))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);

      //StiffnessMatrix_rsPoissonPZ(1,0) = StiffnessMatrix_rsPoissonPZ(0,1);
      StiffnessMatrix_rsPoissonPZ(1,1) = StiffnessMatrix_rsPoissonPZ(0,0);
      StiffnessMatrix_rsPoissonPZ(1,2) = StiffnessMatrix_rsPoissonPZ(0,2);

      //StiffnessMatrix_rsPoissonPZ(2,0) = StiffnessMatrix_rsPoissonPZ(0,2);
      //StiffnessMatrix_rsPoissonPZ(2,1) = StiffnessMatrix_rsPoissonPZ(0,2);
      StiffnessMatrix_rsPoissonPZ(2,2) = -(4*pow(E_p,2)*pow(E_z,2)*pow((nu_p-1),2)-24*E_p*pow(E_z,3)*pow(nu_pz,2)*(nu_p-1))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);

      StiffnessMatrix_rsPoissonPZ(3,3) = 0;

      StiffnessMatrix_rsPoissonPZ(4,4) = 0;

      StiffnessMatrix_rsPoissonPZ(5,5) = 0;
      PetscFunctionReturn(0);
    }

    /**************************************************************************************************
     *
     * With repect to Poisson's ration in z-direction, nu_z, and Young's modulus in p-direction, E_p
     * 
     ************************************************************************************************/
      virtual PetscErrorCode D_rs_PoissonPZYoungP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from D_rs_ElasticFEMethodTransIso"<<endl;
      StiffnessMatrix_rsPoissonPZYoungP.resize(6);
      StiffnessMatrix_rsPoissonPZYoungP.clear();
  
      StiffnessMatrix_rsPoissonPZYoungP(0,0) = (8*E_p*pow(E_z,2)*pow(nu_pz,3))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);
      StiffnessMatrix_rsPoissonPZYoungP(0,1) =  StiffnessMatrix_rsPoissonPZYoungP(0,0);
      StiffnessMatrix_rsPoissonPZYoungP(0,2) =  (4*pow(E_z,3)*pow(nu_pz,4) - 2*E_p*pow(E_z*nu_pz,2)*(3*nu_p - 3))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);

      //StiffnessMatrix_rsPoissonPZYoungP(1,0) = StiffnessMatrix_rsPoissonPZYoungP(0,1);
      StiffnessMatrix_rsPoissonPZYoungP(1,1) = StiffnessMatrix_rsPoissonPZYoungP(0,0);
      StiffnessMatrix_rsPoissonPZYoungP(1,2) = StiffnessMatrix_rsPoissonPZYoungP(0,2);

      //StiffnessMatrix_rsPoissonPZYoungP(2,0) = StiffnessMatrix_rsPoissonPZYoungP(0,2);
      //StiffnessMatrix_rsPoissonPZYoungP(2,1) = StiffnessMatrix_rsPoissonPZYoungP(0,2);
      StiffnessMatrix_rsPoissonPZYoungP(2,2) = -(8*pow(E_z*nu_pz,3)*(nu_p - 1) - 4*E_p*pow(E_z,2)*nu_pz*pow((nu_p - 1),2))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);

      StiffnessMatrix_rsPoissonPZYoungP(3,3) = 0;

      StiffnessMatrix_rsPoissonPZYoungP(4,4) = 0;

      StiffnessMatrix_rsPoissonPZYoungP(5,5) = 0;
      PetscFunctionReturn(0);
    }
	
    /**************************************************************************************************
     *
     * With repect to Poisson's ration in z-direction, nu_z, and Young's modulus in z-direction, E_z
     * 
     ************************************************************************************************/
      virtual PetscErrorCode D_rs_PoissonPZYoungZ(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from D_rs_ElasticFEMethodTransIso"<<endl;
      StiffnessMatrix_rsPoissonPZYoungZ.resize(6);
      StiffnessMatrix_rsPoissonPZYoungZ.clear();
  
      StiffnessMatrix_rsPoissonPZYoungZ(0,0) =  -(2*pow(E_p,2)*nu_pz*(2*E_z*pow(nu_pz,2) + E_p - E_p*nu_p))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);
      StiffnessMatrix_rsPoissonPZYoungZ(0,1) =  StiffnessMatrix_rsPoissonPZYoungZ(0,0);
      StiffnessMatrix_rsPoissonPZYoungZ(0,2) =  (pow(E_p,2)*(nu_p - 1)*(6*E_z*pow(nu_pz,2) + E_p - E_p*nu_p))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);

      //StiffnessMatrix_rsPoissonPZYoungZ(1,0) = StiffnessMatrix_rsPoissonPZYoungZ(0,1);
      StiffnessMatrix_rsPoissonPZYoungZ(1,1) = StiffnessMatrix_rsPoissonPZYoungZ(0,0);
      StiffnessMatrix_rsPoissonPZYoungZ(1,2) = StiffnessMatrix_rsPoissonPZYoungZ(0,2);

      //StiffnessMatrix_rsPoissonPZYoungZ(2,0) = StiffnessMatrix_rsPoissonPZYoungZ(0,2);
      //StiffnessMatrix_rsPoissonPZYoungZ(2,1) = StiffnessMatrix_rsPoissonPZYoungZ(0,2);
      StiffnessMatrix_rsPoissonPZYoungZ(2,2) = -(8*E_z*nu_pz*pow(E_p*(nu_p - 1),2))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);

      StiffnessMatrix_rsPoissonPZYoungZ(3,3) = 0;

      StiffnessMatrix_rsPoissonPZYoungZ(4,4) = 0;

      StiffnessMatrix_rsPoissonPZYoungZ(5,5) = 0;
      PetscFunctionReturn(0);
    }
	
    /**************************************************************************************************
     *
     * With repect to Poisson's ration in z-direction, nu_pz, and Shear modulus in z-direction, G_zp
     * 
     ************************************************************************************************/
      virtual PetscErrorCode D_rs_PoissonPZShearZP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from D_rs_ElasticFEMethodTransIso"<<endl;
      StiffnessMatrix_rsPoissonPZShearZP.resize(6);
      StiffnessMatrix_rsPoissonPZShearZP.clear();
      PetscFunctionReturn(0);
    }
	
    /***************************************************************************
     *
     * With repect to Young's modulus in p-direction, E_p
     * 
     **************************************************************************/
      virtual PetscErrorCode D_rs_YoungP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
//      cout<<"Hello from TransverseIsotropicStiffnessMatrix_rPoissonP "<<endl;
      StiffnessMatrix_rsYoungP.resize(6);
      StiffnessMatrix_rsYoungP.clear();
  
      StiffnessMatrix_rsYoungP(0,0) = -(4*E_z*E_z*pow(nu_pz,4))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);
      StiffnessMatrix_rsYoungP(0,1) = -(4*E_z*E_z*pow(nu_pz,4))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);
      StiffnessMatrix_rsYoungP(0,2) = (4*E_z*E_z*pow(nu_pz,3)*(nu_p-1))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);

      //StiffnessMatrix_rsYoungP(1,0) = StiffnessMatrix_rsYoungP(0,1);
      StiffnessMatrix_rsYoungP(1,1) = StiffnessMatrix_rsYoungP(0,0);
      StiffnessMatrix_rsYoungP(1,2) = StiffnessMatrix_rsYoungP(0,2);

      //StiffnessMatrix_rsYoungP(2,0) = StiffnessMatrix_rsYoungP(0,2);
      //StiffnessMatrix_rsYoungP(2,1) = StiffnessMatrix_rsYoungP(0,2);
      StiffnessMatrix_rsYoungP(2,2) = -(4*E_z*E_z*nu_pz*nu_pz*pow((nu_p-1),2))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);

      StiffnessMatrix_rsYoungP(3,3) = 0;

      StiffnessMatrix_rsYoungP(4,4) = 0;

      StiffnessMatrix_rsYoungP(5,5) = 0;
      PetscFunctionReturn(0);
    }
	
    /**************************************************************************************************
     *
     * With repect to Young's modulus in p-direction, E_p, and Young's modulus in z-direction, E_z
     * 
     ************************************************************************************************/
      virtual PetscErrorCode D_rs_YoungPYoungZ(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from D_rs_ElasticFEMethodTransIso"<<endl;
      StiffnessMatrix_rsYoungPYoungZ.resize(6);
      StiffnessMatrix_rsYoungPYoungZ.clear();
  
      StiffnessMatrix_rsYoungPYoungZ(0,0) =  (4*E_p*E_z*pow(nu_pz,4))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);
      StiffnessMatrix_rsYoungPYoungZ(0,1) =  StiffnessMatrix_rsYoungPYoungZ(0,0);
      StiffnessMatrix_rsYoungPYoungZ(0,2) =  -(4*E_p*E_z*pow(nu_pz,3)*(nu_p - 1))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);

      //StiffnessMatrix_rsYoungPYoungZ(1,0) = StiffnessMatrix_rsYoungPYoungZ(0,1);
      StiffnessMatrix_rsYoungPYoungZ(1,1) = StiffnessMatrix_rsYoungPYoungZ(0,0);
      StiffnessMatrix_rsYoungPYoungZ(1,2) = StiffnessMatrix_rsYoungPYoungZ(0,2);

      //StiffnessMatrix_rsYoungPYoungZ(2,0) = StiffnessMatrix_rsYoungPYoungZ(0,2);
      //StiffnessMatrix_rsYoungPYoungZ(2,1) = StiffnessMatrix_rsYoungPYoungZ(0,2);
      StiffnessMatrix_rsYoungPYoungZ(2,2) = (4*E_p*E_z*pow(nu_pz*(nu_p - 1),2))/pow((2*E_z*pow(nu_pz,2) - E_p + E_p*nu_p),3);

      StiffnessMatrix_rsYoungPYoungZ(3,3) = 0;

      StiffnessMatrix_rsYoungPYoungZ(4,4) = 0;

      StiffnessMatrix_rsYoungPYoungZ(5,5) = 0;
      PetscFunctionReturn(0);
    }
	
     /**************************************************************************************************
     *
     * With repect to Young's modulus in p-direction, E_p, and Shear modulus in z-direction, G_zp
     * 
     ************************************************************************************************/
      virtual PetscErrorCode D_rs_YoungPShearZP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from D_rs_ElasticFEMethodTransIso"<<endl;
      StiffnessMatrix_rsYoungPShearZP.resize(6);
      StiffnessMatrix_rsYoungPShearZP.clear();
      PetscFunctionReturn(0);
    }
	
    /***************************************************************************
     *
     * With repect to Young's modulus in z-direction, E_z
     * 
     **************************************************************************/
      virtual PetscErrorCode D_rs_YoungZ(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from TransverseIsotropicStiffnessMatrix_rsYoungZ "<<endl;
      StiffnessMatrix_rsYoungZ.resize(6);
      StiffnessMatrix_rsYoungZ.clear();
  
      StiffnessMatrix_rsYoungZ(0,0) = -(4*pow(E_p,2)*pow(nu_pz,4))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);
      StiffnessMatrix_rsYoungZ(0,1) = -(4*pow(E_p,2)*pow(nu_pz,4))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);
      StiffnessMatrix_rsYoungZ(0,2) = (4*pow(E_p,2)*pow(nu_pz,3)*(nu_p-1))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);

      //StiffnessMatrix_rsYoungZ(1,0) = StiffnessMatrix_rsYoungZ(0,0);
      StiffnessMatrix_rsYoungZ(1,1) = StiffnessMatrix_rsYoungZ(0,1);
      StiffnessMatrix_rsYoungZ(1,2) = StiffnessMatrix_rsYoungZ(0,2);

      //StiffnessMatrix_rsYoungZ(2,0) = StiffnessMatrix_rsYoungZ(0,2);
      //StiffnessMatrix_rsYoungZ(2,1) = StiffnessMatrix_rsYoungZ(0,2);
      StiffnessMatrix_rsYoungZ(2,2) = -(4*pow(E_p,2)*pow(nu_pz,2)*pow(nu_p-1,2))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);

      StiffnessMatrix_rsYoungZ(3,3) = 0;

      StiffnessMatrix_rsYoungZ(4,4) = 0;

      StiffnessMatrix_rsYoungZ(5,5) = 0;
      PetscFunctionReturn(0);
    }
	
    /**************************************************************************************************
     *
     * With repect to Young's modulus in z-direction, E_z, and Shear modulus in z-direction, G_zp
     * 
     ************************************************************************************************/
      virtual PetscErrorCode D_rs_YoungZShearZP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      PetscFunctionBegin;
      //cout<<"Hello from D_rs_ElasticFEMethodTransIso"<<endl;
      StiffnessMatrix_rsYoungZShearZP.resize(6);
      StiffnessMatrix_rsYoungZShearZP.clear();
      PetscFunctionReturn(0);
    }

    /***************************************************************************
     *
     * With repect to Shear modulus in z-direction, G_zp
     * 
     **************************************************************************/
      virtual PetscErrorCode D_rs_ShearZP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
      // cout<<"Hello from TransverseIsotropicStiffnessMatrix_rsShearZP"<<endl;
      PetscFunctionBegin;
      StiffnessMatrix_rsShearZP.resize(6);
      StiffnessMatrix_rsShearZP.clear();
      PetscFunctionReturn(0);
    }
    
  };

// =============================================================================
//
//  ISOTROPIC MATERIAL
//
// =============================================================================
  struct IsotropicStiffnessMatrix_SecondOrderDerivative {
    double young, nu;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsPoisson;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsYoung;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsYoungPoisson;
    /***************************************************************************
     *
     * With repect to Young's modulus, E
     * 
     **************************************************************************/
    virtual PetscErrorCode D_rs_Young(double young, double nu){
      PetscFunctionBegin;
      StiffnessMatrix_rsYoung.resize(6);
      StiffnessMatrix_rsYoung.clear();

      PetscFunctionReturn(0);
    }

    /***************************************************************************
     *
     * With repect to Poisson's ratio, NU
     * 
     **************************************************************************/
    virtual PetscErrorCode D_rs_Poisson(double young, double nu){
      PetscFunctionBegin;
      StiffnessMatrix_rsPoisson.resize(6);
      StiffnessMatrix_rsPoisson.clear();

      double D_rs00, D_rs01, D_rs33;

      D_rs00 = -(4*young*(- 2*nu*nu*nu + 6*nu*nu + 1))/pow((2*nu*nu + nu - 1),3);
      D_rs01 = -(2*young*(4*nu*nu*nu + 6*nu + 1))/pow((2*nu*nu + nu - 1),3);
      D_rs33 = young/pow((nu + 1),3);

      StiffnessMatrix_rsPoisson(0,0) = D_rs00;  
      StiffnessMatrix_rsPoisson(0,1) = D_rs01;   
      StiffnessMatrix_rsPoisson(0,2) = D_rs01;

      //StiffnessMatrix_rsPoisson(1,0) = D_rs01;  
      StiffnessMatrix_rsPoisson(1,1) = D_rs00;   
      StiffnessMatrix_rsPoisson(1,2) = D_rs01;

      //StiffnessMatrix_rsPoisson(2,0) = D_rs01;  
      //StiffnessMatrix_rsPoisson(2,1) = D_rs01;   
      StiffnessMatrix_rsPoisson(2,2) = D_rs00;

      StiffnessMatrix_rsPoisson(3,3) = D_rs33;

      StiffnessMatrix_rsPoisson(4,4) = D_rs33;

      StiffnessMatrix_rsPoisson(5,5) = D_rs33;

      PetscFunctionReturn(0);
    }

    /***************************************************************************
     *
     * With repect to Young's modulus and Poisson's ratio, NU
     * 
     **************************************************************************/
    virtual PetscErrorCode D_rs_YoungPoisson(double young, double nu){
      PetscFunctionBegin;
      StiffnessMatrix_rsYoungPoisson.resize(6);
      StiffnessMatrix_rsYoungPoisson.clear();

      double D_rs00, D_rs01, D_rs33;

      D_rs00 = -(2*nu*(nu - 2))/pow((2*nu*nu + nu - 1),2);
      D_rs01 = (2*nu*nu + 1)/pow((2*nu*nu + nu - 1),2);
      D_rs33 = -1/(2*pow((nu + 1),2));

      StiffnessMatrix_rsYoungPoisson(0,0) = D_rs00;  
      StiffnessMatrix_rsYoungPoisson(0,1) = D_rs01;   
      StiffnessMatrix_rsYoungPoisson(0,2) = D_rs01;

      //StiffnessMatrix_rsYoungPoisson(1,0) = D_rs01;  
      StiffnessMatrix_rsYoungPoisson(1,1) = D_rs00;   
      StiffnessMatrix_rsYoungPoisson(1,2) = D_rs01;

      //StiffnessMatrix_rsYoungPoisson(2,0) = D_rs01;  
      //StiffnessMatrix_rsYoungPoisson(2,1) = D_rs01;   
      StiffnessMatrix_rsYoungPoisson(2,2) = D_rs00;

      StiffnessMatrix_rsYoungPoisson(3,3) = D_rs33;

      StiffnessMatrix_rsYoungPoisson(4,4) = D_rs33;

      StiffnessMatrix_rsYoungPoisson(5,5) = D_rs33;

      PetscFunctionReturn(0);
    }
};
  
  // ===========================================================================
  //
  // TWO-PHASE YARN - TRANSVERSELY ISOTROPIC MATERIAL
  //
  // ===========================================================================
  struct YarnStiffnessMatrix_Geom_SecondOrderDerivative {
		// fibre waviness
		double ddI1, ddI3, ddI5, ddI6, ddI8;
		double ampltidue, length; // waviness parameters: waviness amplitude & periodic length
		// fibre misalignment
		double theta; // misalignment angle
		double ddmn, ddm2, ddn2, ddmn3, ddm3n, ddm2n2, ddm4, ddn4;
		// fibre volume fraction
		double vf;
		double kf, lf, mf, nf, pf; // fibre related Hill's moduli
		double km, lm, mm, nm, pm; // matrix related Hill's moduli
		double ddkc, ddmc, ddpc, ddlc, ddnc; // fibre volume fraction 
    
    /***************************************************************************
     *
     * With repect to amplitude of waviness
     *
     **************************************************************************/
    virtual PetscErrorCode D_rs_AmpAmp(double amplitude, double length) {
      
      PetscFunctionBegin;
      
      double alpha;
	    alpha = 2*M_PI*amplitude/length;
      
      ddI1 = (2*pow(alpha,4) + 13*pow(alpha,2) - 4)/(2*pow(alpha*alpha + 1,3.5))
             *pow(2*M_PI/length,2);
      ddI3 = (2*pow(alpha,4) - 11*pow(alpha,2) + 2)/(2*pow(alpha*alpha + 1,3.5))
             *pow(2*M_PI/length,2);
      ddI5 = (-6*pow(alpha,4) + 9*pow(alpha,2))/(2*pow(alpha*alpha + 1,3.5))
             *pow(2*M_PI/length,2);
      ddI6 = (2*pow(alpha,2) - 1)/pow(alpha*alpha + 1,2.5)*pow(2*M_PI/length,2);
      ddI8 = -(2*pow(alpha,2) - 1)/pow(alpha*alpha + 1,2.5)*pow(2*M_PI/length,2);
      
      PetscFunctionReturn(0);
    }
    
    /***************************************************************************
     *
     * With repect to periodic length
     *
     **************************************************************************/
    virtual PetscErrorCode D_rs_LenLen(double amplitude, double length) {
      
      PetscFunctionBegin;
	  
	    double alpha;
	    alpha = 2*M_PI*amplitude/length;
      
      ddI1 = (2*pow(alpha,4) + 13*pow(alpha,2) - 4)/(2*pow(alpha*alpha + 1,3.5))
             *pow(2*M_PI*amplitude/length/length,2) +
             -(pow(alpha,3) + 4*alpha)/(2*pow(alpha*alpha + 1,2.5))
             *(4*M_PI*amplitude/pow(length,3));
      ddI3 = (2*pow(alpha,4) - 11*pow(alpha,2) + 1)/(2*pow(alpha*alpha + 1,3.5))
            *pow(2*M_PI*amplitude/length/length,2) +
            -(pow(alpha,3) - 2*alpha)/(2*pow(alpha*alpha + 1,2.5))
            *(4*M_PI*amplitude/pow(length,3));
      ddI5 = (-6*pow(alpha,4) + 9*pow(alpha,2))/(2*pow(alpha*alpha + 1,3.5))
            *pow(2*M_PI*amplitude/length/length,2) +
            (3*pow(alpha,3))/(2*pow(alpha*alpha + 1,2.5))
            *(4*M_PI*amplitude/pow(length,3));
      ddI6 = (2*pow(alpha,2) - 1)/(pow(alpha*alpha + 1,2.5))
             *pow(2*M_PI*amplitude/length/length,2) +
             (-alpha)/(pow(alpha*alpha + 1,1.5))
             *(4*M_PI*amplitude/pow(length,3));
      ddI8 = -(2*pow(alpha,2) - 1)/(pow(alpha*alpha + 1,2.5))
             *pow(2*M_PI*amplitude/length/length,2) +
             (alpha)/(pow(alpha*alpha + 1,1.5))
             *(4*M_PI*amplitude/pow(length,3));
      
      PetscFunctionReturn(0);
      
    }
    
    /***************************************************************************
     *
     * With repect to periodic length
     *
     **************************************************************************/
    virtual PetscErrorCode D_rs_AmpLen(double amplitude, double length) {
      
      PetscFunctionBegin;

	    double alpha;
	    alpha = 2*M_PI*amplitude/length;
      
      ddI1 = (2*pow(alpha,4) + 13*pow(alpha,2) - 4)/(2*pow(alpha*alpha + 1,3.5))*(-2*M_PI*alpha/length/length) +
             -(pow(alpha,3) + 4*alpha)/(2*pow(alpha*alpha + 1,2.5))*(-2*M_PI/pow(length,2));
      ddI3 = (2*pow(alpha,4) - 11*pow(alpha,2) + 1)/(2*pow(alpha*alpha + 1,3.5))*(-2*M_PI*alpha/length/length) +
             -(pow(alpha,3) - 2*alpha)/(2*pow(alpha*alpha + 1,2.5))*(-2*M_PI/pow(length,2));
      ddI5 = (-6*pow(alpha,4) + 9*pow(alpha,2))/(2*pow(alpha*alpha + 1,3.5))*(-2*M_PI*alpha/length/length) +
             (3*pow(alpha,3))/(2*pow(alpha*alpha + 1,2.5))*(-2*M_PI/pow(length,2));
      ddI6 = (2*pow(alpha,2) - 1)/(pow(alpha*alpha + 1,2.5))*(-2*M_PI*alpha/length/length) +
             2*M_PI*alpha/(pow(alpha*alpha + 1,1.5))/pow(length,2);
      ddI8 = -(2*pow(alpha,2) - 1)/(pow(alpha*alpha + 1,2.5))*(-2*M_PI*alpha/length/length) -
             (2*M_PI*alpha)/(pow(alpha*alpha + 1,1.5))/pow(length,2);
      
      PetscFunctionReturn(0);
    }
		
		/***************************************************************************
     *
     * With repect to misalignment angle
     *
     **************************************************************************/
		virtual PetscErrorCode D_rs_Misalignment(double theta) {
			PetscFunctionBegin;
			
			ddmn   = -4*cos(theta)*sin(theta); // mn
			ddm2   = 2*(pow(sin(theta),2) - pow(cos(theta),2)); // m2
			ddn2   = 2*(pow(cos(theta),2) - pow(sin(theta),2)); // n2
			ddmn3  = -10*pow(sin(theta),3)*cos(theta)+6*pow(cos(theta),3)*sin(theta); // mn3
			ddm3n  = 6*cos(theta)*pow(sin(theta),3)-10*pow(cos(theta),3)*sin(theta); // m3n
			ddm2n2 = 2*pow(sin(theta),4)-12*pow(sin(theta),2)*pow(cos(theta),2)+2*pow(cos(theta),4); // m2n2
			ddm4   = 4*pow(cos(theta),4)-12*pow(sin(theta),2)*pow(cos(theta),2); // m4
			ddn4   = 12*pow(sin(theta),2)*pow(cos(theta),2) - 4*pow(sin(theta),4); // n4
			
			PetscFunctionReturn(0);
		}
		
		/***************************************************************************
		 *
		 * With repect to misalignment angle
		 *
		 **************************************************************************/
		 virtual PetscErrorCode D_rs_Fraction(double vf, 
																					 double kf, double mf, double pf, double lf, double nf,
																					 double km, double mm, double pm, double lm, double nm) {
				PetscFunctionBegin;
			
				ddkc = (2*(kf + mm)*(km + mm)*pow(kf - km,2))/pow(kf + mm - vf*kf + vf*km,3);
				ddmc = (2*(km*mm*(vf*mf - mm*(vf - 1)) + mf*mm*(km + 2*mm))*pow(km+2*mm,2)*pow(mf-mm,2))
				       /pow((km*mm + (vf*mm - mf*(vf - 1))*(km + 2*mm)),3) 
							 + (2*km*mm*(km + 2*mm)*pow(mf - mm,2))/pow((km*mm + (vf*mm - mf*(vf - 1))*(km + 2*mm)),2);
				ddpc = (4*pm*(pf + pm)*pow((pf - pm),2))/pow((pf + pm - vf*pf + vf*pm),3);
				ddlc = (2*(kf + mm)*(km + mm)*(kf - km)*(lf - lm))/pow((kf + mm - vf*kf + vf*km),3);
				ddnc = (2*(kf + mm)*(km + mm)*pow((lf - lm),2))/pow((kf + mm - vf*kf + vf*km),3);
				
			  PetscFunctionReturn(0);
      }
  };
}
#endif //__D_RS_ELASTICFEMETHODTRANSISO_HPP__

