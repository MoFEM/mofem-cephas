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

#ifndef __MATERIALCONSTITUTIVEMATRIX_RULEOFMIXTURE_HPP__
#define __MATERIALCONSTITUTIVEMATRIX_RULEOFMIXTURE_HPP__

#include <boost/numeric/ublas/symmetric.hpp>

namespace MoFEM {
  
// =============================================================================
//
//  CONSTITUTIVE MATRIX OF UD COMPOSITE
//  RULE OF MIXTURE
//
// =============================================================================  
  struct TransverseIsotropicStiffnessMatrix_RuleOfMixture {
    
    double nu_23, nu_21, E_1, E_2, G_12;
    
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_NUpf;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_NUpzf;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_Epf;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_Ezf;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_Gzpf;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_Em;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_r_NUm;
    
    virtual PetscErrorCode EEP_RM(double nu_p, double nu_pz, double E_p, double E_z,
                                  double G_zp, double E_m, double nu_m, double V_f) {
      
      PetscFunctionBegin;
      
      // Use the rule of mixture to estimate effective material properties
      double G_m, G_p, nu_zp;
      G_m    = E_m/(1 + nu_m)/2;
      G_p    = E_p/(1 + nu_p)/2;
      nu_zp  = E_z/E_p*nu_pz;
      
      double G_23, nu_12;
      E_1   = E_z*V_f + E_m*(1-V_f);
      E_2   = 1/(V_f/E_p + (1-V_f)/E_m);
      G_12  = 1/(V_f/G_zp + (1 - V_f)/G_m);
      G_23  = G_m/(1 - V_f*(1 - G_m/G_p));
      nu_23 = E_2/G_23/2 - 1;
      nu_12 = nu_zp*V_f + nu_m*(1-V_f);
      nu_21 = E_2/E_1*nu_12;
      
      // Use <TransverseIsotropicStiffnessMatrix> to calculate constitutive matrix
      TransverseIsotropicStiffnessMatrix ROM_Dmat(nu_23, nu_21, E_2, E_1, G_12);
      //cout<<"\nDmat from Rule of Mixture: "<<ROM_Dmat.StiffnessMatrix<<endl;
      // Rotate fibre direction from z-axis to x-axis
      StiffnessMatrix = ROM_Dmat.StiffnessMatrix;
      FibreDirection_z2x(StiffnessMatrix);
      //cout<<"\nAfter coordination transformation: "<<StiffnessMatrix<<endl;
      
      PetscFunctionReturn(0);
    }
    
    /***************************************************************************
     *
     * With repect to Poisson's ration in p-direction
     *
     **************************************************************************/
    virtual PetscErrorCode D_r_PoissonP(double nu_p, double nu_pz, double E_p,
                                        double E_z, double G_zp, double E_m,
                                        double nu_m, double V_f) {
      PetscFunctionBegin;
      
      TransverseIsotropicStiffnessMatrix_FirstOrderDerivative ROM_Dmat_r;
      ROM_Dmat_r.D_r_PoissonP(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_PoissonPZ(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_YoungP(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_YoungZ(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_ShearZP(nu_23, nu_21, E_2, E_1, G_12);
      
      StiffnessMatrix_r_NUpf.resize(6);
      StiffnessMatrix_r_NUpf.clear();
      
      double dE1dNUp, dE2dNUp, dNU21dNUp, dG12dNUp, dNU23dNUp;
      dE1dNUp = 0;
      dE2dNUp = 0;
      dNU21dNUp = 0;
      dG12dNUp = 0;
      
      double G_m, G_p, G_23;
      G_m    = E_m/(1 + nu_m)/2;
      G_p    = E_p/(1 + nu_p)/2;
      G_23   = G_m/(1 - V_f*(1 - G_m/G_p));
      double dG23dNUp, dGpdNUp;
      dGpdNUp = -E_p/2/pow((1 + nu_p),2);
      dG23dNUp = G_m/pow((1 - V_f*(1 - G_m/G_p)),2)*(-V_f*G_m/pow(G_p,2))*dGpdNUp;
      dNU23dNUp = -E_2/2/pow(G_23,2)*dG23dNUp;
      
      StiffnessMatrix_r_NUpf = ROM_Dmat_r.StiffnessMatrix_rPoissonP*dNU23dNUp;
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_NUpf);
      //cout<<"\nDmat wrt NUp: "<<StiffnessMatrix_r_NUpf<<endl;
      
      PetscFunctionReturn(0);
    }

    /***************************************************************************
     *
     * With repect to Poisson's ration in z-direction, nu_pz
     *
     **************************************************************************/
    virtual PetscErrorCode D_r_PoissonPZ(double nu_p, double nu_pz, double E_p,
                                         double E_z, double G_zp, double E_m,
                                         double nu_m, double V_f) {
      PetscFunctionBegin;
      //cout<<"Hello from TransverseIsotropicStiffnessMatrix_rPoissonPZ "<<endl;
      StiffnessMatrix_r_NUpzf.resize(6);
      StiffnessMatrix_r_NUpzf.clear();
      
      TransverseIsotropicStiffnessMatrix_FirstOrderDerivative ROM_Dmat_r;
      ROM_Dmat_r.D_r_PoissonP(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_PoissonPZ(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_YoungP(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_YoungZ(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_ShearZP(nu_23, nu_21, E_2, E_1, G_12);
      
      double dE1dNUpz, dE2dNUpz, dNU21dNUpz, dG12dNUpz, dNU23dNUpz;
      dE1dNUpz   = 0;
      dE2dNUpz   = 0;
      dNU21dNUpz = E_2/E_1*V_f*E_z/E_p;
      dG12dNUpz  = 0;
      dNU23dNUpz = 0;
      
      StiffnessMatrix_r_NUpzf = ROM_Dmat_r.StiffnessMatrix_rPoissonPZ*dNU21dNUpz;
      
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_NUpzf);
      //cout<<"\nDmat wrt NUpz: "<<StiffnessMatrix_r_NUpzf<<endl;
      
      PetscFunctionReturn(0);
    }

    /***************************************************************************
     *
     * With repect to Young's modulus in p-direction, E_p
     *
     **************************************************************************/
    virtual PetscErrorCode D_r_YoungP(double nu_p, double nu_pz, double E_p,
                                      double E_z, double G_zp, double E_m,
                                      double nu_m, double V_f) {
      PetscFunctionBegin;
      StiffnessMatrix_r_Epf.resize(6);
      StiffnessMatrix_r_Epf.clear();
      
      TransverseIsotropicStiffnessMatrix_FirstOrderDerivative ROM_Dmat_r;
      ROM_Dmat_r.D_r_PoissonP(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_PoissonPZ(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_YoungP(nu_23, nu_21, E_2, E_1, G_12);
      
      double dE2dEp, dNU21dEp, dNU23dEp;
      double nu_zp = E_z/E_p*nu_pz;
      dE2dEp   = -1/pow((V_f/E_p + (1-V_f)/E_m),2)*(-V_f/pow(E_p,2));
      
      double nu_12, dNU12dEp;
      nu_12 = nu_zp*V_f + nu_m*(1-V_f);
      dNU12dEp = -V_f*E_z*nu_pz/pow(E_p,2);
      dNU21dEp = dE2dEp/E_1*nu_12 + E_2/E_1*dNU12dEp;
      
      double G_23, G_p;
      double G_m = E_m/(1 + nu_m)/2;
      G_p     = E_p/(1 + nu_p)/2;
      G_23   = G_m/(1 - V_f*(1 - G_m/G_p));
      double dG23dEp, dGpdEp;
      dGpdEp = 1/(1 + nu_p)/2;
      dG23dEp = -G_m/pow((1 - V_f*(1 - G_m/G_p)),2)*(-V_f*G_m/pow(G_p,2)*dGpdEp);
      dNU23dEp =  0.5*(dE2dEp*G_23 - E_2*dG23dEp)/pow(G_23,2);
  
      StiffnessMatrix_r_Epf = ROM_Dmat_r.StiffnessMatrix_rYoungP*dE2dEp +
                              ROM_Dmat_r.StiffnessMatrix_rPoissonPZ*dNU21dEp +
                              ROM_Dmat_r.StiffnessMatrix_rPoissonP*dNU23dEp;
      
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_Epf);
      //cout<<"\nDmat wrt Ep: "<<StiffnessMatrix_r_Epf<<endl;
      
      PetscFunctionReturn(0);
    }

    /***************************************************************************
     *
     * With repect to Young's modulus in z-direction, E_z
     *
     **************************************************************************/
    virtual PetscErrorCode D_r_YoungZ(double nu_p, double nu_pz, double E_p,
                                      double E_z, double G_zp, double E_m,
                                      double nu_m, double V_f) {
      PetscFunctionBegin;
      // cout<<"Hello from TransverseIsotropicStiffnessMatrix_rYoungZ "<<endl;
      StiffnessMatrix_r_Ezf.resize(6);
      StiffnessMatrix_r_Ezf.clear();
      
      TransverseIsotropicStiffnessMatrix_FirstOrderDerivative ROM_Dmat_r;
      ROM_Dmat_r.D_r_PoissonPZ(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_YoungZ(nu_23, nu_21, E_2, E_1, G_12);
      
      double dE1dEz, dNU21dEz;
      dE1dEz = V_f;
      double nu_12, nu_zp;
      nu_zp  = E_z/E_p*nu_pz;
      nu_12 = nu_zp*V_f + nu_m*(1-V_f);
      double dNU12dEz;
      dNU12dEz = V_f/E_p*nu_pz;
      dNU21dEz = -E_2/pow(E_1,2)*dE1dEz*nu_12 + E_2/E_1*dNU12dEz;
      
      StiffnessMatrix_r_Ezf = ROM_Dmat_r.StiffnessMatrix_rYoungZ*dE1dEz +
                              ROM_Dmat_r.StiffnessMatrix_rPoissonPZ*dNU21dEz;
      
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_Ezf);
      //cout<<"\nDmat wrt Ez: "<<StiffnessMatrix_r_Ezf<<endl;
      PetscFunctionReturn(0);
    }
    
    /***************************************************************************
     *
     * With repect to Shear modulus in z-direction, G_zp
     *
     **************************************************************************/
    virtual PetscErrorCode D_r_ShearZP(double nu_p, double nu_pz, double E_p,
                                       double E_z, double G_zp, double E_m,
                                       double nu_m, double V_f) {
      PetscFunctionBegin;
      //cout<<"Hello from TransverseIsotropicStiffnessMatrix_rShearZP"<<endl;
      StiffnessMatrix_r_Gzpf.resize(6);
      StiffnessMatrix_r_Gzpf.clear();
      
      TransverseIsotropicStiffnessMatrix_FirstOrderDerivative ROM_Dmat_r;
      ROM_Dmat_r.D_r_ShearZP(nu_23, nu_21, E_2, E_1, G_12);
      
      double dG12dGzp;
      double G_m;
      G_m    = E_m/(1 + nu_m)/2;
      dG12dGzp = -1/pow((V_f/G_zp + (1 - V_f)/G_m), 2)*(-V_f/pow(G_zp,2));
      
      StiffnessMatrix_r_Gzpf = ROM_Dmat_r.StiffnessMatrix_rShearZP*dG12dGzp;
      
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_Gzpf);
      //cout<<"\nDmat wrt Gzp: "<<StiffnessMatrix_r_Gzpf<<endl;
      
      PetscFunctionReturn(0);
    }

    /***************************************************************************
     *
     * With repect to Young's modulus
     *
     **************************************************************************/
    virtual PetscErrorCode D_r_YoungM(double nu_p, double nu_pz, double E_p,
                                      double E_z, double G_zp, double E_m,
                                      double nu_m, double V_f) {
      PetscFunctionBegin;
      StiffnessMatrix_r_Em.resize(6);
      StiffnessMatrix_r_Em.clear();
      
      TransverseIsotropicStiffnessMatrix_FirstOrderDerivative ROM_Dmat_r;
      ROM_Dmat_r.D_r_PoissonP(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_PoissonPZ(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_YoungP(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_YoungZ(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_ShearZP(nu_23, nu_21, E_2, E_1, G_12);
      
      double dE1dEm, dE2dEm, dNU21dEm, dG12dEm, dNU23dEm;
      dE1dEm = (1-V_f);
      dE2dEm = -1/pow((V_f/E_p + (1-V_f)/E_m), 2)*(-1*(1-V_f)/pow(E_m,2));
      double G_23, G_m, G_p;
      G_m    = E_m/(1 + nu_m)/2;
      G_p     = E_p/(1 + nu_p)/2;
      G_23   = G_m/(1 - V_f*(1 - G_m/G_p));
      
      double dGmdEm, dG23dEm;
      dGmdEm   = 1/(1 + nu_m)/2;
      dG23dEm  = (dGmdEm*(1 - V_f*(1 - G_m/G_p)) - G_m*V_f*dGmdEm/G_p)/pow((1 - V_f*(1 - G_m/G_p)), 2);
      dNU23dEm =  0.5*(dE2dEm*G_23 - E_2*dG23dEm)/pow(G_23,2);
      
      double nu_12, nu_zp;
      nu_zp  = E_z/E_p*nu_pz;
      nu_12 = nu_zp*V_f + nu_m*(1-V_f);
      dNU21dEm = nu_12*(dE2dEm*E_1 - E_2*dE1dEm)/pow(E_1,2);
      
      dG12dEm  = 1/pow((V_f/G_zp + (1 - V_f)/G_m),2)*(1-V_f)/pow(G_m,2)*dGmdEm;
      
      StiffnessMatrix_r_Em = ROM_Dmat_r.StiffnessMatrix_rPoissonP*dNU23dEm +
      ROM_Dmat_r.StiffnessMatrix_rPoissonPZ*dNU21dEm +
      ROM_Dmat_r.StiffnessMatrix_rYoungP*dE2dEm +
      ROM_Dmat_r.StiffnessMatrix_rYoungZ*dE1dEm +
      ROM_Dmat_r.StiffnessMatrix_rShearZP*dG12dEm;
      
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_Em);
      //cout<<"\nDmat wrt Em: "<<StiffnessMatrix_r_Em<<endl;
      
      PetscFunctionReturn(0);
    }
    
    /***************************************************************************
     *
     * With repect to Poisson's ratio
     *
     **************************************************************************/
    virtual PetscErrorCode D_r_PoissonM(double nu_p, double nu_pz, double E_p,
                                        double E_z, double G_zp, double E_m,
                                        double nu_m, double V_f) {
      PetscFunctionBegin;
      StiffnessMatrix_r_NUm.resize(6);
      StiffnessMatrix_r_NUm.clear();
      
      TransverseIsotropicStiffnessMatrix_FirstOrderDerivative ROM_Dmat_r;
      ROM_Dmat_r.D_r_PoissonP(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_PoissonPZ(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_YoungP(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_YoungZ(nu_23, nu_21, E_2, E_1, G_12);
      ROM_Dmat_r.D_r_ShearZP(nu_23, nu_21, E_2, E_1, G_12);
      
      double dNU21dNUm, dG12dNUm, dNU23dNUm;
      dNU21dNUm = E_2*(1-V_f)/E_1;
      
      double G_m, dGmdNUm;
      G_m    = E_m/(1 + nu_m)/2;
      dGmdNUm = -2*E_m/pow(2*(1 + nu_m),2);
      dG12dNUm  = 1/pow((V_f/G_zp + (1 - V_f)/G_m),2)*(1 - V_f)/pow(G_m,2)*dGmdNUm;
      
      double G_p, G_23, dG23dNUm;
      G_p   = E_p/(1 + nu_p)/2;
      G_23  = G_m/(1 - V_f*(1 - G_m/G_p));
      dG23dNUm   =  (dGmdNUm*(1 - V_f*(1 - G_m/G_p)) - G_m*V_f*dGmdNUm/G_p)/pow((1 - V_f*(1 - G_m/G_p)),2) ;
      dNU23dNUm = -E_2/2/pow(G_23,2)*dG23dNUm;
      
      StiffnessMatrix_r_NUm =  ROM_Dmat_r.StiffnessMatrix_rPoissonP*dNU23dNUm +
      ROM_Dmat_r.StiffnessMatrix_rPoissonPZ*dNU21dNUm +
      ROM_Dmat_r.StiffnessMatrix_rShearZP*dG12dNUm;
      
      // Rotate fibre direction from z-axis to x-axis
      FibreDirection_z2x(StiffnessMatrix_r_NUm);
      //cout<<"\nDmat wrt NUm: "<<StiffnessMatrix_r_NUm<<endl;
      
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

