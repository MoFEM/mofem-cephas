/* Copyright (C) 2014, 
 *   Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
 *   Xiao-Yi Zhou (xiaoyi.zhou AT newcastle.ac.uk)
 * --------------------------------------------------------------
 * This routine calculates the second-order partial derivative of right-hand 
 * side, Rhs, with respect to considered random variable of material parameter
 * in the perturbation-based stochastic finite element analysis for two-phase 
 * (isotropic matrix + isotropic/transversely-isotropic reinforcement/inclusion
 * /fibre) composite materials.
 *
 * HISTORY
 *
 * 2014.09.12 (first version)
 *
 * REFERENCES
 * 1. Kleiber M. and Hien T. D. (1992) The stochastic finite element method - 
 *      Basic perturbation technique and computer implementation. John Wiley & 
 *      Sons.
 * 2. Kaczamarczyk L., Pearce C. J. and Bicanic N. (2008) Scale transition and 
 *      enforcement of RVE boundary conditions in second-order computational
 *      homogenization. International Journal for Numerical Methods in 
 *      Engineering, 74(3) p506-522. 
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

#ifndef __TRANS_ISO_RHS_RS_PSFEM_HPP__
#define __TRANS_ISO_RHS_RS_PSFEM_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
extern "C" {
#include <gm_rule.h>
}

namespace MoFEM {

  struct Trans_Iso_Rhs_rs_PSFEM: public TranIsotropicFibreDirRotElasticFEMethod {
    
    double young, pois; // young's modulus and poisson's ratio for isotropic material 
    Vec ddF;
    const string second_field;
    const string ix_first_randvar, ix_second_randvar; // index for considered random variables 
    const string material_type;     // Type of material: isotropic or transversly isotropic
    const string material_function; // Function of material: matrix or reinforcement/inclusion/fibre
    
    Trans_Iso_Rhs_rs_PSFEM(FieldInterface& _mField,Mat &_Aij,Vec _D,Vec _F,const string& _second_field,const string& _ix_first_randvar, const string& _ix_second_randvar, const string& _material_type, const string& _material_function):
    TranIsotropicFibreDirRotElasticFEMethod(_mField,_Aij,_D,_F),ddF(_F),second_field(_second_field),ix_first_randvar(_ix_first_randvar),ix_second_randvar(_ix_second_randvar),material_type(_material_type),material_function(_material_function){};

    
    
// =============================================================================
//
// Calculate material constitutive matrix in global coordinate system for
// transversely isotropic material
//
// =============================================================================
    virtual PetscErrorCode calculateD_r_Trans(double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp) {
      PetscFunctionBegin;
      ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix_r;
      ConstitutiveMatrix_r.resize(6);
      ConstitutiveMatrix_r.clear();

      /*************************************************************************
       *
       * Get the constitutive matrix
       *
       ************************************************************************/ 
      TransverseIsotropicStiffnessMatrix_FirstOrderDerivative mymat;
      if (ix_first_randvar.compare(0,8,"PoissonP") == 0){
         // cout<<"first variable poisson p \t";
         ierr = mymat.D_r_PoissonP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_r = mymat.StiffnessMatrix_rPoissonP;
      }
      else if (ix_first_randvar.compare(0,8,"PoissonZ") == 0){
         // cout<<"first variable poisson z \t";
         ierr = mymat.D_r_PoissonPZ(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_r = mymat.StiffnessMatrix_rPoissonPZ;
      }
      else if (ix_first_randvar.compare(0,6,"YoungP") == 0){
         // cout<<"first variable Young p \t";
         ierr = mymat.D_r_YoungP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_r = mymat.StiffnessMatrix_rYoungP;
      }
      else if (ix_first_randvar.compare(0,6,"YoungZ") == 0){
         // cout<<"first variable Young z \t";
         ierr = mymat.D_r_YoungZ(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_r = mymat.StiffnessMatrix_rYoungZ;
      }
      else if (ix_first_randvar.compare(0,7,"ShearZP") == 0){
         // cout<<"first variable Shear z \t";
         ierr = mymat.D_r_ShearZP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_r = mymat.StiffnessMatrix_rShearZP;
      }
      else {
           cout<<"Invalid input of random variable"<<endl;
      }

      /*************************************************************************
       * 
       * Constitutive matrix transformation
       * For textile, the constitutive matrix of fibre needs to transform to the
       * global coordinate system using element direction information obtained 
       * through 'potential flow' calculation
       *
       ************************************************************************/      
      // Rotating the Stiffness matrix according a set of axes of rotations and their respective angle
      // cout<<material_function<<"\t";
      D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
      
      if (material_function.compare(0,13,"reinforcement")==0){
         //cout<<"Reinforcement \t";
      	 vector< ublas::matrix< double > > normalized_phi;
      	 normalized_phi.resize(coords_at_Gauss_nodes.size());
         ierr = ComputeFibreDirection(normalized_phi); CHKERRQ(ierr);
      
      	 for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
        
            int noOfRotations = 1; //Number of Rotations
        
            double zVec[3]={0.0,0.0,1.0};
            double AxVector[3]={normalized_phi[gg](0,1)*zVec[2]-normalized_phi[gg](0,2)*zVec[1] , normalized_phi[gg](0,2)*zVec[0]-normalized_phi[gg](0,0)*zVec[2] , normalized_phi[gg](0,0)*zVec[1]-normalized_phi[gg](0,1)*zVec[0]};
            double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))/(sqrt(pow(normalized_phi[gg](0,0),2)+pow(normalized_phi[gg](0,1),2)+pow(normalized_phi[gg](0,2),2)))*(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};
        
            double negAxAngle[noOfRotations];
            for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];
        
            ublas::matrix<double> DummyMatrix,DummyMatrix2;
            DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
            DummyMatrix = ConstitutiveMatrix_r;
        
            ///Rotating Stiffness over a number of axis/angle rotations
            for (int aa=0; aa<noOfRotations; aa++) {
          
                StressTransformation StressRotMat(&AxVector[3*aa], AxAngle[aa]);
                StrainTransformation invStrainRotMat(&AxVector[3*aa], negAxAngle[aa]);
          
                ublas::matrix<double> TrpMatrixStress;
                TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
                TrpMatrixStress=StressRotMat.StressRotMat;
          
                ublas::matrix<double> TrpMatrixInvStrain;
                TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6);
                TrpMatrixInvStrain=invStrainRotMat.StrainRotMat;
          
                DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
                ublas::matrix< FieldData > dummyA = prod( DummyMatrix , TrpMatrixInvStrain );
                DummyMatrix2 = prod(TrpMatrixStress,dummyA);
                DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
                DummyMatrix = DummyMatrix2;
            }
        
            D_At_GaussPoint[gg].resize(6,6);
            D_At_GaussPoint[gg].clear();
            D_At_GaussPoint[gg] = DummyMatrix;
         }
      }
      else if (material_function.compare(0,6,"matrix")==0){
         cout<<"Matrix \t";
         for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
            D_At_GaussPoint[gg].resize(6,6);
            D_At_GaussPoint[gg].clear();
            D_At_GaussPoint[gg] = ConstitutiveMatrix_r;
         }
      }
      else {
         cout<<"Undefined material function!"<<endl;
      }
      
      PetscFunctionReturn(0);
    }
    //**************************************************************************

    virtual PetscErrorCode calculateD_rs_Trans(double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp) {
      PetscFunctionBegin;
      
      ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix_rs;
      ConstitutiveMatrix_rs.resize(6);
      ConstitutiveMatrix_rs.clear();
      /*************************************************************************
       *
       * Get the constitutive matrix
       *
       ************************************************************************/
      TransverseIsotropicStiffnessMatrix_SecondOrderDerivative mymat;
      if (ix_first_randvar.compare(0,8,"PoissonP") == 0 && ix_second_randvar.compare(0,8,"PoissonP") == 0){
         // cout<<"second variable Poisson p \t";
         ierr = mymat.D_rs_PoissonP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsPoissonP;
      }
	  else if (ix_first_randvar.compare(0,8,"PoissonP") == 0 && ix_second_randvar.compare(0,8,"PoissonZ") == 0){
         // cout<<"second variable Poisson z \t";
         ierr = mymat.D_rs_PoissonPPoissonPZ(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsPoissonPPoissonPZ;
      }
	  else if (ix_first_randvar.compare(0,8,"PoissonP") == 0 && ix_second_randvar.compare(0,6,"YoungP") == 0){
         // cout<<"second variable Poisson z \t";
         ierr = mymat.D_rs_PoissonPYoungP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsPoissonPYoungP;
      }
	  else if (ix_first_randvar.compare(0,8,"PoissonP") == 0 && ix_second_randvar.compare(0,6,"YoungZ") == 0){
         // cout<<"second variable Poisson z \t";
         ierr = mymat.D_rs_PoissonPYoungZ(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsPoissonPYoungZ;
      }
	  else if (ix_first_randvar.compare(0,8,"PoissonP") == 0 && ix_second_randvar.compare(0,7,"ShearZP") == 0){
         // cout<<"second variable Poisson z \t";
         ierr = mymat.D_rs_PoissonPShearZP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsPoissonPShearZP;
      }
      else if (ix_first_randvar.compare(0,8,"PoissonZ") == 0 && ix_second_randvar.compare(0,8,"PoissonZ") == 0){
         // cout<<"second variable Poisson z \t";
         ierr = mymat.D_rs_PoissonPZ(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsPoissonPZ;
      }
	  else if (ix_first_randvar.compare(0,8,"PoissonZ") == 0 && ix_second_randvar.compare(0,6,"YoungP") == 0){
         // cout<<"second variable Poisson z \t";
         ierr = mymat.D_rs_PoissonPZYoungP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsPoissonPZYoungP;
      }
	  else if (ix_first_randvar.compare(0,8,"PoissonZ") == 0 && ix_second_randvar.compare(0,6,"YoungZ") == 0){
         // cout<<"second variable Poisson z \t";
         ierr = mymat.D_rs_PoissonPZYoungZ(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsPoissonPZYoungZ;
      }
	  else if (ix_first_randvar.compare(0,8,"PoissonZ") == 0 && ix_second_randvar.compare(0,7,"ShearZP") == 0){
         // cout<<"second variable Poisson z \t";
         ierr = mymat.D_rs_PoissonPZShearZP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsPoissonPZShearZP;
      }
      else if (ix_first_randvar.compare(0,6,"YoungP") == 0 && ix_second_randvar.compare(0,6,"YoungP") == 0){
         // cout<<"second variable Young p \t";
         ierr = mymat.D_rs_YoungP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsYoungP;
      }
	  else if (ix_first_randvar.compare(0,6,"YoungP") == 0 && ix_second_randvar.compare(0,6,"YoungZ") == 0){
         // cout<<"second variable Young p \t";
         ierr = mymat.D_rs_YoungPYoungZ(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsYoungPYoungZ;
      }
	  else if (ix_first_randvar.compare(0,6,"YoungP") == 0 && ix_second_randvar.compare(0,7,"ShearZP") == 0){
         // cout<<"second variable Young p \t";
         ierr = mymat.D_rs_YoungPShearZP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsYoungPShearZP;
      }
      else if (ix_first_randvar.compare(0,6,"YoungZ") == 0 && ix_second_randvar.compare(0,6,"YoungZ") == 0){
         // cout<<"second variable Young Z \t";
         ierr = mymat.D_rs_YoungZ(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsYoungZ;
      }
	  else if (ix_first_randvar.compare(0,6,"YoungZ") == 0 && ix_second_randvar.compare(0,7,"ShearZP") == 0){
         // cout<<"second variable Young p \t";
         ierr = mymat.D_rs_YoungZShearZP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsYoungZShearZP;
      }
      else if (ix_first_randvar.compare(0,7,"ShearZP") == 0 && ix_second_randvar.compare(0,7,"ShearZP") == 0){
         // cout<<"second variable Shear zp \t";
         ierr = mymat.D_rs_ShearZP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsShearZP;
      }
	  else if (ix_second_randvar.compare(0,6,"YoungM") == 0 || ix_second_randvar.compare(0,8,"PoissonM") == 0){
         // cout<<"second variable Poisson z \t";
         ConstitutiveMatrix_rs.clear();
      }
      else {
           cout<<"Invalid input of random variable"<<endl;
      }
      
      /*************************************************************************
       * 
       * Constitutive matrix transformation
       * For textile, the constitutive matrix of fibre needs to transform to the
       * global coordinate system using element direction information obtained 
       * through 'potential flow' calculation
       *
       ************************************************************************/      
      // Rotating the Stiffness matrix according a set of axes of rotations and their respective angle
      // cout<<material_function<<"\t";
      D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
      
      if (material_function.compare(0,13,"reinforcement")==0){
         //cout<<"Reinforcement \t";
      	 vector< ublas::matrix< double > > normalized_phi;
      	 normalized_phi.resize(coords_at_Gauss_nodes.size());
         ierr = ComputeFibreDirection(normalized_phi); CHKERRQ(ierr);
      
      	 for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
        
            int noOfRotations = 1; //Number of Rotations
        
            double zVec[3]={0.0,0.0,1.0};
            double AxVector[3]={normalized_phi[gg](0,1)*zVec[2]-normalized_phi[gg](0,2)*zVec[1] , normalized_phi[gg](0,2)*zVec[0]-normalized_phi[gg](0,0)*zVec[2] , normalized_phi[gg](0,0)*zVec[1]-normalized_phi[gg](0,1)*zVec[0]};
            double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))/(sqrt(pow(normalized_phi[gg](0,0),2)+pow(normalized_phi[gg](0,1),2)+pow(normalized_phi[gg](0,2),2)))*(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};
        
            double negAxAngle[noOfRotations];
            for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];
        
            ublas::matrix<double> DummyMatrix,DummyMatrix2;
            DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
            DummyMatrix = ConstitutiveMatrix_rs;
        
            ///Rotating Stiffness over a number of axis/angle rotations
            for (int aa=0; aa<noOfRotations; aa++) {
          
                StressTransformation StressRotMat(&AxVector[3*aa], AxAngle[aa]);
                StrainTransformation invStrainRotMat(&AxVector[3*aa], negAxAngle[aa]);
          
                ublas::matrix<double> TrpMatrixStress;
                TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
                TrpMatrixStress=StressRotMat.StressRotMat;
          
                ublas::matrix<double> TrpMatrixInvStrain;
                TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6);
                TrpMatrixInvStrain=invStrainRotMat.StrainRotMat;
          
                DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
                ublas::matrix< FieldData > dummyA = prod( DummyMatrix , TrpMatrixInvStrain );
                DummyMatrix2 = prod(TrpMatrixStress,dummyA);
                DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
                DummyMatrix = DummyMatrix2;
            }
        
            D_At_GaussPoint[gg].resize(6,6);
            D_At_GaussPoint[gg].clear();
            D_At_GaussPoint[gg] = DummyMatrix;
         }
      }
      else if (material_function.compare(0,6,"matrix")==0){
         cout<<"Matrix \t";
         for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
            D_At_GaussPoint[gg].resize(6,6);
            D_At_GaussPoint[gg].clear();
            D_At_GaussPoint[gg] = ConstitutiveMatrix_rs;
         }
      }
      else {
         cout<<"Undefined material function!"<<endl;
      }
    
    PetscFunctionReturn(0);
  }

// =============================================================================
//
// Calculate material constitutive matrix in global coordinate system for
// isotropic material
//
// =============================================================================
    virtual PetscErrorCode calculateD_r_Iso(double _young, double _nu) {
      PetscFunctionBegin;
      ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix_r;
      ConstitutiveMatrix_r.resize(6);
      ConstitutiveMatrix_r.clear();

      /*************************************************************************
       *
       * Get the constitutive matrix
       *
       ************************************************************************/ 
      IsotropicStiffnessMatrix_FirstOrderDerivative mymat;
      if (ix_first_randvar.compare(0,5,"Young") == 0 ){
         //cout<<"Young's modulus of isotropic material"<<endl;
         ierr = mymat.D_r_Young(_young,_nu); CHKERRQ(ierr);
         ConstitutiveMatrix_r = mymat.StiffnessMatrix_rYoung;
      }
      else if (ix_first_randvar.compare(0,7,"Poisson") == 0){
         //cout<<"Poisson's ratio of isotropic material"<<endl;
         ierr = mymat.D_r_Poisson(_young,_nu); CHKERRQ(ierr);
         ConstitutiveMatrix_r = mymat.StiffnessMatrix_rPoisson;
      }
      else {
         cout<<"Invalid input of random variable"<<endl;
      }

      /*************************************************************************
       * 
       * Constitutive matrix transformation
       * For textile, the constitutive matrix of fibre needs to transform to the
       * global coordinate system using element direction information obtained 
       * through 'potential flow' calculation
       *
       ************************************************************************/      
      // Rotating the Stiffness matrix according a set of axes of rotations and their respective angle
      // cout<<material_function<<"\t";
      D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
      
      if (material_function.compare(0,13,"reinforcement")==0){
         //cout<<"Reinforcement \t";
      	 vector< ublas::matrix< double > > normalized_phi;
      	 normalized_phi.resize(coords_at_Gauss_nodes.size());
         ierr = ComputeFibreDirection(normalized_phi); CHKERRQ(ierr);
      
      	 for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
        
            int noOfRotations = 1; //Number of Rotations
        
            double zVec[3]={0.0,0.0,1.0};
            double AxVector[3]={normalized_phi[gg](0,1)*zVec[2]-normalized_phi[gg](0,2)*zVec[1] , normalized_phi[gg](0,2)*zVec[0]-normalized_phi[gg](0,0)*zVec[2] , normalized_phi[gg](0,0)*zVec[1]-normalized_phi[gg](0,1)*zVec[0]};
            double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))/(sqrt(pow(normalized_phi[gg](0,0),2)+pow(normalized_phi[gg](0,1),2)+pow(normalized_phi[gg](0,2),2)))*(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};
        
            double negAxAngle[noOfRotations];
            for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];
        
            ublas::matrix<double> DummyMatrix,DummyMatrix2;
            DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
            DummyMatrix = ConstitutiveMatrix_r;
        
            ///Rotating Stiffness over a number of axis/angle rotations
            for (int aa=0; aa<noOfRotations; aa++) {
          
                StressTransformation StressRotMat(&AxVector[3*aa], AxAngle[aa]);
                StrainTransformation invStrainRotMat(&AxVector[3*aa], negAxAngle[aa]);
          
                ublas::matrix<double> TrpMatrixStress;
                TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
                TrpMatrixStress=StressRotMat.StressRotMat;
          
                ublas::matrix<double> TrpMatrixInvStrain;
                TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6);
                TrpMatrixInvStrain=invStrainRotMat.StrainRotMat;
          
                DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
                ublas::matrix< FieldData > dummyA = prod( DummyMatrix , TrpMatrixInvStrain );
                DummyMatrix2 = prod(TrpMatrixStress,dummyA);
                DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
                DummyMatrix = DummyMatrix2;
            }
        
            D_At_GaussPoint[gg].resize(6,6);
            D_At_GaussPoint[gg].clear();
            D_At_GaussPoint[gg] = DummyMatrix;
         }
      }
      else if (material_function.compare(0,6,"matrix")==0){
         // cout<<"Matrix \t";
         for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
            D_At_GaussPoint[gg].resize(6,6);
            D_At_GaussPoint[gg].clear();
            D_At_GaussPoint[gg] = ConstitutiveMatrix_r;
         }
      }
      else {
         cout<<"Undefined material function!"<<endl;
      }
      
      PetscFunctionReturn(0);
    }
    //**************************************************************************

    virtual PetscErrorCode calculateD_rs_Iso(double _young, double _nu) {
      PetscFunctionBegin;
      
      ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix_rs;
      ConstitutiveMatrix_rs.resize(6);
      ConstitutiveMatrix_rs.clear();
      /*************************************************************************
       *
       * Get the constitutive matrix
       *
       ************************************************************************/
      IsotropicStiffnessMatrix_SecondOrderDerivative mymat;
      if (ix_first_randvar.compare(0,7,"Poisson") == 0 && ix_second_randvar.compare(0,7,"Poisson") == 0){
         // cout<<"second variable Poisson p \t";
         ierr = mymat.D_rs_Poisson(_young,_nu); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsPoisson;
      }
      else if (ix_first_randvar.compare(0,5,"Young") == 0 && ix_second_randvar.compare(0,5,"Young") == 0){
         // cout<<"second variable Young p \t";
         ierr = mymat.D_rs_Young(_young,_nu); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsYoung;
      }
      else if (ix_first_randvar.compare(0,5,"Young") == 0 && ix_second_randvar.compare(0,7,"Poisson") == 0){
         // cout<<"second variable Young Z \t";
         ierr = mymat.D_rs_YoungPoisson(_young,_nu); CHKERRQ(ierr);
         ConstitutiveMatrix_rs = mymat.StiffnessMatrix_rsYoungPoisson;
      }
      else {
           cout<<"Invalid input of random variable"<<endl;
      }
      
      /*************************************************************************
       * 
       * Constitutive matrix transformation
       * For textile, the constitutive matrix of fibre needs to transform to the
       * global coordinate system using element direction information obtained 
       * through 'potential flow' calculation
       *
       ************************************************************************/      
      // Rotating the Stiffness matrix according a set of axes of rotations and their respective angle
      // cout<<material_function<<"\t";
      D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
      
      if (material_function.compare(0,13,"reinforcement")==0){
         //cout<<"Reinforcement \t";
      	 vector< ublas::matrix< double > > normalized_phi;
      	 normalized_phi.resize(coords_at_Gauss_nodes.size());
         ierr = ComputeFibreDirection(normalized_phi); CHKERRQ(ierr);
      
      	 for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
        
            int noOfRotations = 1; //Number of Rotations
        
            double zVec[3]={0.0,0.0,1.0};
            double AxVector[3]={normalized_phi[gg](0,1)*zVec[2]-normalized_phi[gg](0,2)*zVec[1] , normalized_phi[gg](0,2)*zVec[0]-normalized_phi[gg](0,0)*zVec[2] , normalized_phi[gg](0,0)*zVec[1]-normalized_phi[gg](0,1)*zVec[0]};
            double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))/(sqrt(pow(normalized_phi[gg](0,0),2)+pow(normalized_phi[gg](0,1),2)+pow(normalized_phi[gg](0,2),2)))*(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};
        
            double negAxAngle[noOfRotations];
            for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];
        
            ublas::matrix<double> DummyMatrix,DummyMatrix2;
            DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
            DummyMatrix = ConstitutiveMatrix_rs;
        
            ///Rotating Stiffness over a number of axis/angle rotations
            for (int aa=0; aa<noOfRotations; aa++) {
          
                StressTransformation StressRotMat(&AxVector[3*aa], AxAngle[aa]);
                StrainTransformation invStrainRotMat(&AxVector[3*aa], negAxAngle[aa]);
          
                ublas::matrix<double> TrpMatrixStress;
                TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
                TrpMatrixStress=StressRotMat.StressRotMat;
          
                ublas::matrix<double> TrpMatrixInvStrain;
                TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6);
                TrpMatrixInvStrain=invStrainRotMat.StrainRotMat;
          
                DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
                ublas::matrix< FieldData > dummyA = prod( DummyMatrix , TrpMatrixInvStrain );
                DummyMatrix2 = prod(TrpMatrixStress,dummyA);
                DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
                DummyMatrix = DummyMatrix2;
            }
        
            D_At_GaussPoint[gg].resize(6,6);
            D_At_GaussPoint[gg].clear();
            D_At_GaussPoint[gg] = DummyMatrix;
         }
      }
      else if (material_function.compare(0,6,"matrix")==0){
         // cout<<"Matrix \t";
         for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
            D_At_GaussPoint[gg].resize(6,6);
            D_At_GaussPoint[gg].clear();
            D_At_GaussPoint[gg] = ConstitutiveMatrix_rs;
         }
      }
      else {
         cout<<"Undefined material function!"<<endl;
      }
    
    PetscFunctionReturn(0);
  }

// =============================================================================
//
// Calculate element stiffness matrix
//
// =============================================================================
    ublas::matrix<ublas::matrix<FieldData> > K_r;
    virtual PetscErrorCode StiffnessK_r() {
      PetscFunctionBegin;
      /**********************************
       *
       * Get material constitutive matrix
       *
       *********************************/
      if (material_type.compare(0,5,"trans") == 0){
         // cout<<"Transversely isotropic material \t";
         double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
         ierr = GetMatParameters(&_E_p,&_E_z,&_nu_p,&_nu_pz,&_G_zp); CHKERRQ(ierr);
         ierr = calculateD_r_Trans(_E_p,_E_z,_nu_p,_nu_pz,_G_zp); CHKERRQ(ierr);
//       cout<<"D_At_GaussPoint[0] "<< D_At_GaussPoint[0] <<endl;
      }
      else if (material_type.compare(0,3,"iso") == 0){
         // cout<<"Isotropic material \t";
         double _young,_pois;
  	 for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)){
            string name = it->get_name();
            if (name.compare(0,13,"MAT_ELASTIC_1") == 0){
               Mat_Elastic mydata;
               ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
               //cout << mydata;
              _young = mydata.data.Young;
              _pois  = mydata.data.Poisson;
            }
         }
         ierr = calculateD_r_Iso(_young,_pois); CHKERRQ(ierr);
         //cout<<"_young "<< _young <<endl;
         //cout<<"_pois "<< _pois <<endl;
      }
      else {
         cout<<"Only isotropic or transversely isotropic material is considered currently!"<<endl;
      }

//      cout<<"D_At_GaussPoint[0] "<< D_At_GaussPoint[0] <<endl;
      K_r.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        for(int gg = 0;gg<g_dim;gg++) {
          ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
          double w = V*G_TET_W[gg];
          if(detH.size()>0) {
            w *= detH[gg];
          }
          BD.resize(6,row_Mat.size2());
          cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
                      BD.size1(),BD.size2(),
                      w,&*D_At_GaussPoint[gg].data().begin(),D_At_GaussPoint[gg].size2(),
                      &*row_Mat.data().begin(),row_Mat.size2(),
                      0.,&*BD.data().begin(),BD.size2());
          for(int cc = 0;cc<col_mat;cc++) {
            if(ColGlob[cc].size()==0) continue;
            ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
            if(gg == 0) {
              K_r(rr,cc).resize(BD.size2(),col_Mat.size2());
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          0.,&*K_r(rr,cc).data().begin(),K_r(rr,cc).size2());
            } else {
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          1.,&*K_r(rr,cc).data().begin(),K_r(rr,cc).size2());
            }
          }
        }
      }
      PetscFunctionReturn(0);
    }

    
    ublas::matrix<ublas::matrix<FieldData> > K_rs;
    virtual PetscErrorCode StiffnessK_rs() {
      PetscFunctionBegin;
      /**********************************
       *
       * Get material constitutive matrix
       *
       *********************************/
      if (material_type.compare(0,5,"trans") == 0){
         // cout<<"Transversely isotropic material \t";
         double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
         ierr = GetMatParameters(&_E_p,&_E_z,&_nu_p,&_nu_pz,&_G_zp); CHKERRQ(ierr);
         ierr = calculateD_rs_Trans(_E_p,_E_z,_nu_p,_nu_pz,_G_zp); CHKERRQ(ierr);
//       cout<<"D_At_GaussPoint[0] "<< D_At_GaussPoint[0] <<endl;
      }
      else if (material_type.compare(0,3,"iso") == 0){
         // cout<<"Isotropic material \t";
         double _young,_pois;
  	 for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)){
            string name = it->get_name();
            if (name.compare(0,13,"MAT_ELASTIC_1") == 0){
               Mat_Elastic mydata;
               ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
              _young = mydata.data.Young;
              _pois  = mydata.data.Poisson;
            }
         }
         ierr = calculateD_rs_Iso(_young,_pois); CHKERRQ(ierr);
         //cout<<"_young "<< _young <<endl;
         //cout<<"_pois "<< _pois <<endl;
      }
      else {
         cout<<"Only isotropic or transversely isotropic material is considered currently!"<<endl;
      }
//      cout<<"D_At_GaussPoint[0] "<< D_At_GaussPoint[0] <<endl;
      K_rs.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        for(int gg = 0;gg<g_dim;gg++) {
          ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
          double w = V*G_TET_W[gg];
          if(detH.size()>0) {
            w *= detH[gg];
          }
          BD.resize(6,row_Mat.size2());
          cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
                      BD.size1(),BD.size2(),
                      w,&*D_At_GaussPoint[gg].data().begin(),D_At_GaussPoint[gg].size2(),
                      &*row_Mat.data().begin(),row_Mat.size2(),
                      0.,&*BD.data().begin(),BD.size2());
          for(int cc = 0;cc<col_mat;cc++) {
            if(ColGlob[cc].size()==0) continue;
            ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
            if(gg == 0) {
              K_rs(rr,cc).resize(BD.size2(),col_Mat.size2());
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          0.,&*K_rs(rr,cc).data().begin(),K_rs(rr,cc).size2());
            } else {
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          1.,&*K_rs(rr,cc).data().begin(),K_rs(rr,cc).size2());
            }
          }
        }
      }
      PetscFunctionReturn(0);
      
    }


// =============================================================================
//
// Calculate element force using as right-hand side of finite element 
// equilibrium equation
//
// F: calculate the second-order partial derivative of "external force" which
//    is referred as right-hand side in the algebriac equation [K][U_rs] = [F_rs]
//    [F_rs] = - [K_rs][U] - 2[K_r][U_s]
// =============================================================================
  vector<ublas::vector<FieldData> > f_el_rs; // element force
  virtual PetscErrorCode Rhs() {
    PetscFunctionBegin;
    //cout<<" Rhs() "<<endl;
    ierr = StiffnessK_r(); CHKERRQ(ierr);    // get K_r
    ierr = StiffnessK_rs(); CHKERRQ(ierr);   // get K_rs

    // displacements for nodes in each element and,
    // first-order derivative of displacements for nodes in each element
    vector<ublas::vector<FieldData> > D_elm;
    vector<ublas::vector<FieldData> > D_elm_r;
//     cout<<"col_mat = "<< col_mat << endl;
    D_elm.resize(col_mat);
    D_elm_r.resize(col_mat);
    
    int col_mat1 = 0;  //only nodes (1st order)
    ierr = GetDataVector("DISPLACEMENT",D_elm[col_mat1]); CHKERRQ(ierr);
    ierr = GetDataVector(second_field,D_elm_r[col_mat1]); CHKERRQ(ierr);
//    cout<<"D_elm[col_mat] = "<< D_elm[col_mat1] << endl;
    col_mat1++;
    
    for(int ee=0; ee<6; ee++) { //edges
      if(ColGlob[col_mat1].size()!=0) {
        ierr = GetDataVector("DISPLACEMENT",MBEDGE,D_elm[col_mat1],ee); CHKERRQ(ierr);
        ierr = GetDataVector(second_field,MBEDGE,D_elm_r[col_mat1],ee); CHKERRQ(ierr);
//          cout<<"Edges D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
        col_mat1++;
      }
    }
    
    for(int ff=0; ff<4; ff++) { //faces
      if(ColGlob[col_mat1].size()!=0) {
        ierr = GetDataVector("DISPLACEMENT",MBTRI,D_elm[col_mat1],ff); CHKERRQ(ierr);
        ierr = GetDataVector(second_field,MBTRI,D_elm_r[col_mat1],ff); CHKERRQ(ierr);
//          cout<<"Faces D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
        col_mat1++;
      }
    }
    
    if(ColGlob[col_mat1].size()!=0) { // volumes
      ierr = GetDataVector("DISPLACEMENT",MBTET,D_elm[col_mat1]); CHKERRQ(ierr);
      ierr = GetDataVector(second_field,MBTET,D_elm_r[col_mat1]); CHKERRQ(ierr);
//    cout<<"Faces D_elm[col_mat] = "<< D_elm[col_mat1] << endl;
    }
    
    // calculate element nodal forces, f_el_rs
    f_el_rs.resize(row_mat);
    for(int rr = 0;rr<row_mat;rr++) {
      if(RowGlob[rr].size()==0) continue;
      int rr_start=0;
      for(int cc = 0;cc<col_mat;cc++) {
        if(ColGlob[cc].size()==0) continue;
//          cout<<"rr "<<rr<<endl;
//          cout<<"cc "<<cc<<endl;
        if(rr_start == 0) {
//            cout<<"K_r(rr,cc) "<<K_r(rr,cc)<<endl;
//            cout<<"D_elm_r[cc] "<<D_elm_r[cc]<<endl;
          f_el_rs[rr] = - prod( K_rs(rr,cc), D_elm[cc] )
          - 2*prod(K_r(rr,cc), D_elm_r[cc]);
          rr_start++;
        } else {
          f_el_rs[rr] -= prod( K_rs(rr,cc), D_elm[cc] )
          + 2*prod(K_r(rr,cc), D_elm_r[cc]);
        }
      }
//      cout<<"f_el_rs[rr] "<<f_el_rs[rr]<<endl;
    }
    
    // assemble the obtained element nodal forces, fe_rs, into the
    // global nodal force vector, ddF,
    for(int rr = 0;rr<row_mat;rr++) {
      if(RowGlob[rr].size()==0) continue;
      if(RowGlob[rr].size()!=f_el_rs[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      ierr = VecSetValues(ddF,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_el_rs[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

// =============================================================================
//
// 
//
// =============================================================================
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

// =============================================================================
//
// 
//
// =============================================================================
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      
      ierr = Get_g_NTET(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
//      cout<<"Hi from K_rsPoissonP_ElasticFEMethodTransIso "<<endl;
      ierr = GetMatrices(); CHKERRQ(ierr);
      ierr = Rhs(); CHKERRQ(ierr);
//      std::string wait;
//      std::cin >> wait;
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
   
    
  };  

}

#endif //__TRANS_ISO_RHS_RS_PSFEM_HPP__
