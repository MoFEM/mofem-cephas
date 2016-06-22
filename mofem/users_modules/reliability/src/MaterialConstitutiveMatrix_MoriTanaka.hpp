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

#ifndef __MATERIALCONSTITUTIVEMATRIX_MORITANAKA_HPP__
#define __MATERIALCONSTITUTIVEMATRIX_MORITANAKA_HPP__

#include <boost/numeric/ublas/symmetric.hpp>

namespace MoFEM {
  
// =============================================================================
//
//  CONSTITUTIVE MATRIX OF UD COMPOSITE
//  RULE OF MIXTURE
//
// =============================================================================  
  struct TransverseIsotropicStiffnessMatrix_MoriTanaka {
    
    //double nu_23, nu_21, E_1, E_2, G_12;
    
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_NUpf;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_NUpzf;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_Epf;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_Ezf;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_Gzpf;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_Em;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_NUm;
    
    double nu_p, nu_pz, E_p, E_z, G_zp; // engineering parameters for trans-iso
    double lambda, mu;                 // lame parameters for isotropic material
    double vf;                                          // fibre volume fraction
    
    virtual PetscErrorCode  EEP_MT(double nu_p, double nu_pz, double E_p, double E_z, double G_zp,
                                   double lambda, double mu,
                                   double vf){
      PetscFunctionBegin;
      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();
      double k_t, l_t, n_t, m_t, p_t; // Hill's moduli for transversely isotropic material
      double k_i, l_i, n_i, m_i, p_i; // Hill's moduli for isotropic material
      
      // Calculate Hill's moduli for tran-iso material using its engineering parameters
      double nu_zp=(nu_pz*E_z)/E_p;
      double G_p = E_p/(2*(1+nu_p));
      
      k_t = 1/(2*(1-nu_p)/E_p-4*nu_zp*nu_zp/E_z);
      l_t = 2*k_t*nu_zp;
      n_t = E_z+l_t*l_t/k_t;
      m_t = G_p;
      p_t = G_zp; // need to be checked G_21 = G_12 [?]
      
      // Calculate Hill's moduli for isotropic material using its engineering parameters
      double K_value = lambda+2*mu/3;
      double G_value = mu;
      
      k_i = K_value + G_value/3;
      l_i = K_value - 2*G_value/3;
      n_i = K_value + 4*G_value/3;
      m_i = G_value;
      p_i = G_value;
      
      /*
       * Calculate Hill's modulus for two-phase composite with fibre of
       * transversely isotropic material and matrix of isotropic material using
       * Mori-Tanaka asymptotic method based on Dvorak derived simple formula
       */
      
      double k_c, l_c, n_c, m_c, p_c;
      double k_f, l_f, n_f, m_f, p_f;            // Hill's moduli for fibre yarn
      double k_m, l_m, n_m, m_m, p_m;                // Hill's moduli for matrix
      double vm;                                       // matrix volume fraction
      k_f=k_t; l_f=l_t; n_f=n_t; m_f=m_t; p_f=p_t;
      k_m=k_i; l_m=l_i; n_m=n_i; m_m=m_i; p_m=p_i;
      
      //vf=0.6;
      vm = 1-vf;
      //cout<<"Volume fraction: \t"<<vf<<endl;
      
      k_c = (vf*k_f*(k_m+m_m)+vm*k_m*(k_f+m_m))/(vf*(k_m+m_m)+vm*(k_f+m_m));
      m_c = (m_f*m_m*(k_m+2*m_m)+k_m*m_m*(vf*m_f+vm*m_m))/(k_m*m_m+(k_m+2*m_m)*(vf*m_m+vm*m_f));
      p_c = (2*vf*p_f*p_m+vm*(p_f*p_m+p_m*p_m))/(2*vf*p_m+vm*(p_f+p_m));
      l_c = (vf*l_f*(k_m+m_m)+vm*l_m*(k_f+m_m))/(vf*(k_m+m_m)+vm*(p_f+p_m));
      n_c = vf*n_f+vm*n_m+(l_c-vf*l_f-vm*l_m)*(l_f-l_m)/(k_f-k_m);
      
      // case 1: fibre direction in x-axis
      /*StiffnessMatrix(0,0) = n_c;
       StiffnessMatrix(0,1) = StiffnessMatrix(0,2) = l_c;
       StiffnessMatrix(1,1) = StiffnessMatrix(2,2) = k_c + m_c;
       StiffnessMatrix(1,2) = k_c - m_c;
       StiffnessMatrix(3,3) = m_c;
       StiffnessMatrix(4,4) = StiffnessMatrix(5,5) = p_c;*/
      // case 2: fibre direction in z-axis
      StiffnessMatrix(0,0) = StiffnessMatrix(1,1) = k_c + m_c;
      StiffnessMatrix(0,1) = k_c - m_c;
      StiffnessMatrix(0,2) = StiffnessMatrix(1,2) = l_c;
      StiffnessMatrix(2,2) = n_c;
      StiffnessMatrix(3,3) = m_c;
      StiffnessMatrix(4,4) = StiffnessMatrix(5,5) = p_c;
      // cout<<"Yarn C matrix \t"<<StiffnessMatrix<<endl;
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix);
      
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
      double dkc, dpc, dnc, dmc, dlc;
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
      
      // fibre direction in z-axis
      StiffnessMatrix_r_Ezf.resize(6); StiffnessMatrix_r_Ezf.clear();
      StiffnessMatrix_r_Ezf(0,0) = StiffnessMatrix_r_Ezf(1,1) = dkc + dmc;
      StiffnessMatrix_r_Ezf(0,1) = dkc - dmc;
      StiffnessMatrix_r_Ezf(0,2) = StiffnessMatrix_r_Ezf(1,2) = dlc;
      StiffnessMatrix_r_Ezf(2,2) = dnc;
      StiffnessMatrix_r_Ezf(3,3) = dmc;
      StiffnessMatrix_r_Ezf(4,4) = StiffnessMatrix_r_Ezf(5,5) = dpc;
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_Ezf);
      
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
      double dkc, dmc, dpc, dlc, dnc;
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
      
      // fibre direction in z-axis
      StiffnessMatrix_r_Epf.resize(6); StiffnessMatrix_r_Epf.clear();
      StiffnessMatrix_r_Epf(0,0) = StiffnessMatrix_r_Epf(1,1) = dkc + dmc;
      StiffnessMatrix_r_Epf(0,1) = dkc - dmc;
      StiffnessMatrix_r_Epf(0,2) = StiffnessMatrix_r_Epf(1,2) = dlc;
      StiffnessMatrix_r_Epf(2,2) = dnc;
      StiffnessMatrix_r_Epf(3,3) = dmc;
      StiffnessMatrix_r_Epf(4,4) = StiffnessMatrix_r_Epf(5,5) = dpc;
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_Epf);
      
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
      double dkc, dmc, dpc, dnc, dlc;
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
      
      // fibre direction in z-axis
      StiffnessMatrix_r_NUpf.resize(6); StiffnessMatrix_r_NUpf.clear();
      StiffnessMatrix_r_NUpf(0,0) = StiffnessMatrix_r_NUpf(1,1) = dkc + dmc;
      StiffnessMatrix_r_NUpf(0,1) = dkc - dmc;
      StiffnessMatrix_r_NUpf(0,2) = StiffnessMatrix_r_NUpf(1,2) = dlc;
      StiffnessMatrix_r_NUpf(2,2) = dnc;
      StiffnessMatrix_r_NUpf(3,3) = dmc;
      StiffnessMatrix_r_NUpf(4,4) = StiffnessMatrix_r_NUpf(5,5) = dpc;
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_NUpf);
      
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
      double dkc, dmc, dpc, dlc,dnc;
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
      
      // fibre direction in z-axis
      StiffnessMatrix_r_NUpzf.resize(6); StiffnessMatrix_r_NUpzf.clear();
      StiffnessMatrix_r_NUpzf(0,0) = StiffnessMatrix_r_NUpzf(1,1) = dkc + dmc;
      StiffnessMatrix_r_NUpzf(0,1) = dkc - dmc;
      StiffnessMatrix_r_NUpzf(0,2) = StiffnessMatrix_r_NUpzf(1,2) = dlc;
      StiffnessMatrix_r_NUpzf(2,2) = dnc;
      StiffnessMatrix_r_NUpzf(3,3) = dmc;
      StiffnessMatrix_r_NUpzf(4,4) = StiffnessMatrix_r_NUpzf(5,5) = dpc;
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_NUpzf);
      
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
      double dkc, dmc, dpc, dlc, dnc;
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
      
      // fibre direction in z-axis
      StiffnessMatrix_r_Gzpf.resize(6); StiffnessMatrix_r_Gzpf.clear();
      StiffnessMatrix_r_Gzpf(0,0) = StiffnessMatrix_r_Gzpf(1,1) = dkc + dmc;
      StiffnessMatrix_r_Gzpf(0,1) = dkc - dmc;
      StiffnessMatrix_r_Gzpf(0,2) = StiffnessMatrix_r_Gzpf(1,2) = dlc;
      StiffnessMatrix_r_Gzpf(2,2) = dnc;
      StiffnessMatrix_r_Gzpf(3,3) = dmc;
      StiffnessMatrix_r_Gzpf(4,4) = StiffnessMatrix_r_Gzpf(5,5) = dpc;
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_Gzpf);
      
      PetscFunctionReturn(0);
    }
    
    /***************************************************************************
     *
     * With repect to Young's modulus in z-direction, E_z
     *
     **************************************************************************/
    virtual PetscErrorCode D_r_YoungM(double nu_p, double nu_pz, double E_p, double E_z, double G_zp,
                                      double Em, double NUm, double vf){
      
      PetscFunctionBegin;
      double dkc, dmc, dpc, dlc, dnc;
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
      
      // fibre direction in z-axis
      StiffnessMatrix_r_Em.resize(6); StiffnessMatrix_r_Em.clear();
      StiffnessMatrix_r_Em(0,0) = StiffnessMatrix_r_Em(1,1) = dkc + dmc;
      StiffnessMatrix_r_Em(0,1) = dkc - dmc;
      StiffnessMatrix_r_Em(0,2) = StiffnessMatrix_r_Em(1,2) = dlc;
      StiffnessMatrix_r_Em(2,2) = dnc;
      StiffnessMatrix_r_Em(3,3) = dmc;
      StiffnessMatrix_r_Em(4,4) = StiffnessMatrix_r_Em(5,5) = dpc;
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_Em);
      
      PetscFunctionReturn(0);
    }
    /***************************************************************************
     *
     * With repect to Young's modulus in z-direction, E_z
     *
     **************************************************************************/
    virtual PetscErrorCode D_r_PoissonM(double nu_p, double nu_pz, double E_p, double E_z, double G_zp,
                                        double Em, double NUm, double vf){
      
      PetscFunctionBegin;
      double dkc, dmc, dpc, dlc, dnc;
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
      
      // fibre direction in z-axis
      StiffnessMatrix_r_NUm.resize(6); StiffnessMatrix_r_NUm.clear();
      StiffnessMatrix_r_NUm(0,0) = StiffnessMatrix_r_NUm(1,1) = dkc + dmc;
      StiffnessMatrix_r_NUm(0,1) = dkc - dmc;
      StiffnessMatrix_r_NUm(0,2) = StiffnessMatrix_r_NUm(1,2) = dlc;
      StiffnessMatrix_r_NUm(2,2) = dnc;
      StiffnessMatrix_r_NUm(3,3) = dmc;
      StiffnessMatrix_r_NUm(4,4) = StiffnessMatrix_r_NUm(5,5) = dpc;
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_NUm);
      
      PetscFunctionReturn(0);
    }
    
    /***************************************************************************
     *
     * Rotate fibre direction from z-axis to x-axis
     *
     * Different convention of material coordination system is used for
     *   transversly isotropic material. This is for fibre reinforced polymer
     *   composite.
     *
     *   x - fibre direction
     *   y - transverse direction
     *   z - through thickness direction
     *
     * Thus, a coordination transformation is thus conducted, which is to rotate
     *   about y with angle of 90 degree.
     *
     **************************************************************************/
    
    virtual PetscErrorCode FibreDirection_z2x(ublas::symmetric_matrix<FieldData,ublas::upper> &StiffnessMatrix) {
      //virtual PetscErrorCode FibreDirection_z2x(ublas::matrix<double> &StiffnessMatrix) {
      PetscFunctionBegin;
      
      double AxVector_PSFE[3] = {0.0,1.0,0.0};
      double AxAngle_PSFE[1] = {0.5*M_PI};
      double negAxAngle_PSFE[1]; negAxAngle_PSFE[0] = -AxAngle_PSFE[0];
      
      ublas::matrix<double> DummyMatrix_PSFE, DummyMatrix2_PSFE;
      DummyMatrix_PSFE = ublas::zero_matrix<FieldData>(6,6);
      DummyMatrix_PSFE = StiffnessMatrix;
      
      StressTransformation StressRotMat_PSFE(&AxVector_PSFE[0], AxAngle_PSFE[0]);
      StrainTransformation invStrainRotMat_PSFE(&AxVector_PSFE[0], negAxAngle_PSFE[0]);
      
      ublas::matrix<double> TrpMatrixStress_PSFE;
      TrpMatrixStress_PSFE = ublas::zero_matrix<FieldData>(6,6);
      TrpMatrixStress_PSFE = StressRotMat_PSFE.StressRotMat;
      
      ublas::matrix<double> TrpMatrixInvStrain_PSFE;
      TrpMatrixInvStrain_PSFE = ublas::zero_matrix<FieldData>(6,6);
      TrpMatrixInvStrain_PSFE = invStrainRotMat_PSFE.StrainRotMat;
      
      DummyMatrix2_PSFE = ublas::zero_matrix<FieldData>(6,6);
      ublas::matrix< FieldData > dummyA_PSFE = prod( DummyMatrix_PSFE , TrpMatrixInvStrain_PSFE );
      DummyMatrix2_PSFE = prod(TrpMatrixStress_PSFE, dummyA_PSFE);
      DummyMatrix_PSFE = ublas::zero_matrix<FieldData>(6,6);
      DummyMatrix_PSFE = DummyMatrix2_PSFE;
      StiffnessMatrix.clear(); StiffnessMatrix = DummyMatrix_PSFE;
      
      PetscFunctionReturn(0);
    }
  };
}
#endif //__D_R_ELASTICFEMETHODTRANSISO_HPP__

