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
    const string zeroth_field;
    const string first_field;
    const string ix_first_randvar, ix_second_randvar; // index for considered random variables 
    const string material_type;     // Type of material: isotropic or transversly isotropic
    const string material_function; // Function of material: matrix or reinforcement/inclusion/fibre
    
    Trans_Iso_Rhs_rs_PSFEM(FieldInterface& _mField,
                           Mat &_Aij,
                           Vec _D,
                           Vec _F,
                           const string & _zeroth_field,
                           const string& _first_field,
                           const string& _ix_first_randvar,
                           const string& _ix_second_randvar,
                           const string& _material_type,
                           const string& _material_function):
    TranIsotropicFibreDirRotElasticFEMethod(_mField,_Aij,_D,_F,_zeroth_field),ddF(_F),zeroth_field(_zeroth_field),first_field(_first_field),ix_first_randvar(_ix_first_randvar),ix_second_randvar(_ix_second_randvar),material_type(_material_type),material_function(_material_function){};

    
    
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
//      else {
//           cout<<"Invalid input of random variable"<<endl;
//      }

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
//      else {
//           cout<<"Invalid input of random variable"<<endl;
//      }
      
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
//      else {
//         cout<<"Invalid input of random variable"<<endl;
//      }

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
//      else {
//           cout<<"Invalid input of random variable"<<endl;
//      }
      
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
    ierr = GetDataVector(zeroth_field,D_elm[col_mat1]); CHKERRQ(ierr);
    ierr = GetDataVector(first_field,D_elm_r[col_mat1]); CHKERRQ(ierr);
//    cout<<"D_elm[col_mat] = "<< D_elm[col_mat1] << endl;
    col_mat1++;
    
    for(int ee=0; ee<6; ee++) { //edges
      if(ColGlob[col_mat1].size()!=0) {
        ierr = GetDataVector(zeroth_field,MBEDGE,D_elm[col_mat1],ee); CHKERRQ(ierr);
        ierr = GetDataVector(first_field,MBEDGE,D_elm_r[col_mat1],ee); CHKERRQ(ierr);
//          cout<<"Edges D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
        col_mat1++;
      }
    }
    
    for(int ff=0; ff<4; ff++) { //faces
      if(ColGlob[col_mat1].size()!=0) {
        ierr = GetDataVector(zeroth_field,MBTRI,D_elm[col_mat1],ff); CHKERRQ(ierr);
        ierr = GetDataVector(first_field,MBTRI,D_elm_r[col_mat1],ff); CHKERRQ(ierr);
//          cout<<"Faces D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
        col_mat1++;
      }
    }
    
    if(ColGlob[col_mat1].size()!=0) { // volumes
      ierr = GetDataVector(zeroth_field,MBTET,D_elm[col_mat1]); CHKERRQ(ierr);
      ierr = GetDataVector(first_field,MBTET,D_elm_r[col_mat1]); CHKERRQ(ierr);
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


/*****************************************************************************
   *
   * Second-order derivative for right-hand-side with respect to:
   *  1. fibre waviness:
   *      Amp - waviness amplitude
   *      Len - waviness periodic length
   *  2. fibre misalignment:
   *      Theta - misalignment angle
   *  3. fibre volume fraction:
   *      vf - fibre volume fraction
   *
   ****************************************************************************/
  struct Trans_Iso_Geom_Rhs_rs_PSFEM: public TranIso_FibreWavinessElasticFEMethod {
    
    Vec ddF;
    string zeroth_field;
    const string first_field;
    const string ix_first_randvar, ix_second_randvar; // index for considered random variables 
    const string material_type;     // Type of material: isotropic or transversly isotropic
    const string material_function; // Function of material: matrix or reinforcement/inclusion/fibre
    
    Trans_Iso_Geom_Rhs_rs_PSFEM(FieldInterface& _mField,
                                Mat &_Aij,
                                Vec _D,
                                Vec _F,
                                const string& _zeroth_field,
                                const string& _first_field,
                                const string& _ix_first_randvar,
                                const string& _ix_second_randvar,
                                const string& _material_type,
                                const string& _material_function):
    TranIso_FibreWavinessElasticFEMethod(_mField,_Aij,_D,_F,_zeroth_field),ddF(_F),
                                         zeroth_field(_zeroth_field),
                                         first_field(_first_field),
				                                 ix_first_randvar(_ix_first_randvar),
                                         ix_second_randvar(_ix_second_randvar),
                                         material_type(_material_type),
											                   material_function(_material_function){};

// =============================================================================
//
// Calculate material constitutive matrix in global coordinate system for
// transversely isotropic material
//
// =============================================================================
    virtual PetscErrorCode calculateD_r_Waviness(double _E_p, double _E_z, double _nu_p,
                                             double _nu_pz, double _G_zp,
                                             double _lambda, double _mu,
                                             double _vf,
                                             double _theta_f, double _WavinessFactor) {
      PetscFunctionBegin;
	  
	  // -----------
      // 1. Get the original compliance matrix
	  ublas::symmetric_matrix<FieldData,ublas::upper> SMat;
      SMat.resize(6);
      SMat.clear();
      
      //TransverseIsotropicComplianceMatrix TranIsoMat_S(_nu_p,_nu_pz,_E_p,_E_z,_G_zp);
      YarnComplianceMatrix TranIsoMat_S(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,
                                        _lambda, _mu, _vf);
      SMat = TranIsoMat_S.ComplianceMatrix;
	  
	  // -----------
	  // 2. Get the wavniess parameters
	  double WavinessAmplitude, WavinessLength;
	  double dI1, dI3, dI5, dI6, dI8; // 1st-order derivative
	  //WavinessAmplitude = 1.19; // unit: mm
	  WavinessLength = 27.9;   // unit: mm
	  WavinessAmplitude = WavinessLength*_WavinessFactor;
	  
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix;
	  ComplianceMatrix.resize(6);
	  ComplianceMatrix.clear();
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix_r;
	  ComplianceMatrix_r.resize(6);
	  ComplianceMatrix_r.clear();
      
	  YarnStiffnessMatrix_Geom_FirstOrderDerivative mymat;
	  if (ix_first_randvar.compare(0,9,"Amplitude") == 0) {
		ierr = mymat.D_r_Amplitude(WavinessAmplitude,WavinessLength); CHKERRQ(ierr);
		dI1 = mymat.dI1;
		dI3 = mymat.dI3;
		dI5 = mymat.dI5;
		dI6 = mymat.dI6;
		dI8 = mymat.dI8;
	  }
	  else if (ix_first_randvar.compare(0,6,"Length") == 0) {
		ierr = mymat.D_r_Length(WavinessAmplitude,WavinessLength); CHKERRQ(ierr);
		dI1 = mymat.dI1;
		dI3 = mymat.dI3;
		dI5 = mymat.dI5;
		dI6 = mymat.dI6;
		dI8 = mymat.dI8;
	  }
	  else {
        cout<<"Invalid input of random variable"<<endl;
	  }
	  
	  // -----------
      // 3. calculate first-order derivative of transformed compliance matrix
	  // rotate about x-axis
//      ComplianceMatrix_r(0,0) = 0;
//      ComplianceMatrix_r(0,1) = SMat(0,1)*dI6 + SMat(0,2)*dI8;
//      ComplianceMatrix_r(0,2) = SMat(0,2)*dI6 + SMat(0,1)*dI8;
//      ComplianceMatrix_r(1,1) = SMat(0,0)*dI1 + 2*(SMat(0,0)-SMat(0,1)+SMat(0,2))*dI3 + SMat(2,2)*dI5;
//      ComplianceMatrix_r(1,2) = SMat(0,2)*dI1 + (SMat(2,2)-SMat(0,0)+2*SMat(0,1))*dI3 + SMat(0,2)*dI5;
//      ComplianceMatrix_r(2,2) = SMat(2,2)*dI1 + SMat(0,0)*dI5 + 2*(SMat(0,0)-SMat(0,1)+SMat(0,2))*dI3;
//      ComplianceMatrix_r(3,3) = 2*(SMat(0,0)-SMat(0,1))*(dI1-2*dI3+dI5)+4*(SMat(0,0)-2*SMat(0,2)+SMat(2,2))*dI3;
//      ComplianceMatrix_r(4,4) = 0;
//      ComplianceMatrix_r(5,5) = 0;
	  
	  // rotate about y-axis
	  ComplianceMatrix_r(0,0) = SMat(0,0)*dI1 + (2*SMat(0,2)+SMat(5,5))*dI3 + SMat(2,2)*dI5;
	  ComplianceMatrix_r(0,1) = SMat(0,1)*dI6 + SMat(0,2)*dI8;
	  ComplianceMatrix_r(0,2) = SMat(0,2)*(dI1+dI5) + (SMat(0,0)+SMat(2,2)-SMat(5,5))*dI3;
	  ComplianceMatrix_r(1,1) = 0;
	  ComplianceMatrix_r(1,2) = SMat(0,2)*dI6 + SMat(0,1)*dI8;
	  ComplianceMatrix_r(2,2) = SMat(2,2)*dI1 + SMat(0,0)*dI5 + (2*SMat(0,2)+SMat(5,5))*dI3;
	  ComplianceMatrix_r(3,3) = 2*(SMat(0,0)-SMat(0,1))*dI6 + SMat(5,5)*dI8;
	  ComplianceMatrix_r(4,4) = 2*(2*SMat(0,0)-4*SMat(0,2)+2*SMat(2,2)-SMat(5,5))*dI3 + SMat(5,5)*(dI1+dI5);
	  ComplianceMatrix_r(5,5) = SMat(5,5)*dI6 + 2*(SMat(0,0)-SMat(0,1))*dI8;
	  
	  // -----------
	  // 4. retrieve transformed stiffness matrix
	  // 4.1 compliance matrix
      double I1, I3, I5, I6, I8;
      double alpha_w;
      alpha_w = 2*M_PI*_WavinessFactor;
      
      I1 = (1+alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);        // m^44
      I3 = (alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);          // m^2 n^2
      I5 = 1-(1+3*alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);    // n^4
      I6 = 1/sqrt(1+alpha_w*alpha_w);                                 // m^2
      I8 = 1-1/sqrt(1+alpha_w*alpha_w);
	  
      ComplianceMatrix(0,0) = SMat(0,0)*I1 + (2*SMat(0,2)+SMat(5,5))*I3 + SMat(2,2)*I5;
      ComplianceMatrix(0,1) = SMat(0,1)*I6 + SMat(0,2)*I8;
      ComplianceMatrix(0,2) = SMat(0,2)*(I1+I5) + (SMat(0,0)+SMat(2,2)-SMat(5,5))*I3;
      ComplianceMatrix(1,1) = SMat(0,0);
      ComplianceMatrix(1,2) = SMat(0,2)*I6 + SMat(0,1)*I8;
      ComplianceMatrix(2,2) = SMat(2,2)*I1 + SMat(0,0)*I5 + (2*SMat(0,2)+SMat(5,5))*I3;
      ComplianceMatrix(3,3) = 2*(SMat(0,0)-SMat(0,1))*I6 + SMat(5,5)*I8;
      ComplianceMatrix(4,4) = 2*(2*SMat(0,0)-4*SMat(0,2)+2*SMat(2,2)-SMat(5,5))*I3 + SMat(5,5)*(I1+I5);
      ComplianceMatrix(5,5) = SMat(5,5)*I6 + 2*(SMat(0,0)-SMat(0,1))*I8;
	  
	  // 4.2 stiffness matrix
	  
	  ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix;
	  ConstitutiveMatrix.resize(6);
	  ConstitutiveMatrix.clear();
	  
      double SVal;
      double S11, S12, S13, S21, S22, S23, S31, S32, S33;
      S11 = ComplianceMatrix(0,0);
      S12 = ComplianceMatrix(0,1);
      S13 = ComplianceMatrix(0,2);
      
      S21 = ComplianceMatrix(1,0);
      S22 = ComplianceMatrix(1,1);
      S23 = ComplianceMatrix(1,2);
      
      S31 = ComplianceMatrix(2,0);
      S32 = ComplianceMatrix(2,1);
      S33 = ComplianceMatrix(2,2);
      
      SVal = S11*(S22 * S33 - S23 * S32) + S12*(S23 * S31 - S21 * S33)
              + S13*(S21 * S32 - S22 * S31);
      
      ConstitutiveMatrix(0,0) = (ComplianceMatrix(1,1)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(1,2),2))/SVal;
      ConstitutiveMatrix(0,1) = (ComplianceMatrix(0,2)*ComplianceMatrix(1,2) - ComplianceMatrix(0,1)*ComplianceMatrix(2,2))/SVal;
      ConstitutiveMatrix(0,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(1,2) - ComplianceMatrix(0,2)*ComplianceMatrix(1,1))/SVal;
      ConstitutiveMatrix(1,1) = (ComplianceMatrix(0,0)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(0,2),2))/SVal;
      ConstitutiveMatrix(1,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(0,2) - ComplianceMatrix(1,2)*ComplianceMatrix(0,0))/SVal;
      ConstitutiveMatrix(2,2) = (ComplianceMatrix(0,0)*ComplianceMatrix(1,1) - pow(ComplianceMatrix(0,1),2))/SVal;
      ConstitutiveMatrix(3,3) = 1/ComplianceMatrix(3,3);
      ConstitutiveMatrix(4,4) = 1/ComplianceMatrix(4,4);
      ConstitutiveMatrix(5,5) = 1/ComplianceMatrix(5,5);
	  
	  // 5. Calculate derivative of transformed stiffness matrix by using the
	  //     formula: dC/dx = - C [dS/dt] C
	  
	  ublas::matrix<double> ConstitutiveMatrix_r;
      ConstitutiveMatrix_r = ublas::zero_matrix<FieldData>(6,6);
	  
	  ublas::matrix<double> DummyMatrix3;
	  DummyMatrix3 = ublas::zero_matrix<FieldData>(6,6);
	  ublas::matrix< FieldData > dummyAA = prod(ComplianceMatrix_r,ConstitutiveMatrix);
	  DummyMatrix3 = prod(ConstitutiveMatrix, dummyAA);
	  ConstitutiveMatrix_r  = -1*DummyMatrix3;
      //cout<<"Transformed stiffness matrix "<<ConstitutiveMatrix_r<<endl;
	  
	  // -----------
      // 6. Rotate the constitutive matrix
      D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
      //cout<<material_function<<endl;
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
		   //cout<<DummyMatrix<<endl;
         }
      }
      else if (material_function.compare(0,6,"matrix")==0){
		//cout<<"Matrix \t";
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

    virtual PetscErrorCode calculateD_rs_Waviness(double _E_p, double _E_z, double _nu_p,
                                             double _nu_pz, double _G_zp,
                                             double _lambda, double _mu,
                                             double _vf,
                                             double _theta_f, double _WavinessFactor) {
      PetscFunctionBegin;
      
	  
	  // -----------
      // 1. Get the original compliance matrix
	  ublas::symmetric_matrix<FieldData,ublas::upper> SMat;
      SMat.resize(6);
      SMat.clear();
      
      //TransverseIsotropicComplianceMatrix TranIsoMat_S(_nu_p,_nu_pz,_E_p,_E_z,_G_zp);
      YarnComplianceMatrix TranIsoMat_S(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,
                                        _lambda, _mu, _vf);
      SMat = TranIsoMat_S.ComplianceMatrix;
	  
	  // -----------
	  // 2. Get the wavniess parameters
	  double WavinessAmplitude, WavinessLength;
	  double ddI1, ddI3, ddI5, ddI6, ddI8; // 2st-order partial derivative
	  double dI1_r, dI3_r, dI5_r, dI6_r, dI8_r;
	  double dI1_s, dI3_s, dI5_s, dI6_s, dI8_s;
	  //WavinessAmplitude = 1.19; // unit: mm
	  WavinessLength = 27.9;   // unit: mm
	  WavinessAmplitude = WavinessLength*_WavinessFactor;
	  
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix;
      ComplianceMatrix.resize(6);
      ComplianceMatrix.clear();
	  
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix_r;
      ComplianceMatrix_r.resize(6);
      ComplianceMatrix_r.clear();
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix_s;
      ComplianceMatrix_s.resize(6);
      ComplianceMatrix_s.clear();
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix_rs;
      ComplianceMatrix_rs.resize(6);
      ComplianceMatrix_rs.clear();
      
      YarnStiffnessMatrix_Geom_SecondOrderDerivative mymat;
      YarnStiffnessMatrix_Geom_FirstOrderDerivative mymat_1st;
      if (ix_first_randvar.compare(0,9,"Amplitude") == 0 && ix_second_randvar.compare(0,9,"Amplitude")==0) {
        ierr = mymat.D_rs_AmpAmp(WavinessAmplitude,WavinessLength); CHKERRQ(ierr);
		ddI1 = mymat.ddI1;
		ddI3 = mymat.ddI3;
		ddI5 = mymat.ddI5;
		ddI6 = mymat.ddI6;
		ddI8 = mymat.ddI8;
		//
		ierr = mymat_1st.D_r_Amplitude(WavinessAmplitude,WavinessLength); CHKERRQ(ierr);
		dI1_r = mymat_1st.dI1;
		dI3_r = mymat_1st.dI3;
		dI5_r = mymat_1st.dI5;
		dI6_r = mymat_1st.dI6;
		dI8_r = mymat_1st.dI8;
		//
		dI1_s = mymat_1st.dI1;
		dI3_s = mymat_1st.dI3;
		dI5_s = mymat_1st.dI5;
		dI6_s = mymat_1st.dI6;
		dI8_s = mymat_1st.dI8;
		
		//cout<<"ddI1 = "<<ddI1<<"\t ddI3 = "<<ddI3<<"\t ddI5 = "<<ddI5<<"\t ddI6 = "<<ddI6<<"\t ddI8 = "<<ddI8<<endl;
      }
      else if (ix_first_randvar.compare(0,6,"Length") == 0 && ix_second_randvar.compare(0,6,"Length")==0) {
		ierr = mymat.D_rs_LenLen(WavinessAmplitude,WavinessLength); CHKERRQ(ierr);
		ddI1 = mymat.ddI1;
		ddI3 = mymat.ddI3;
		ddI5 = mymat.ddI5;
		ddI6 = mymat.ddI6;
		ddI8 = mymat.ddI8;
		//
		ierr = mymat_1st.D_r_Length(WavinessAmplitude,WavinessLength); CHKERRQ(ierr);
		dI1_r = mymat_1st.dI1;
		dI3_r = mymat_1st.dI3;
		dI5_r = mymat_1st.dI5;
		dI6_r = mymat_1st.dI6;
		dI8_r = mymat_1st.dI8;
		//
		dI1_s = mymat_1st.dI1;
		dI3_s = mymat_1st.dI3;
		dI5_s = mymat_1st.dI5;
		dI6_s = mymat_1st.dI6;
		dI8_s = mymat_1st.dI8;
      }
      else if (ix_first_randvar.compare(0,9,"Amplitude") == 0 && ix_second_randvar.compare(0,6,"Length")==0) {
		ierr = mymat.D_rs_AmpLen(WavinessAmplitude,WavinessLength); CHKERRQ(ierr);
		ddI1 = mymat.ddI1;
		ddI3 = mymat.ddI3;
		ddI5 = mymat.ddI5;
		ddI6 = mymat.ddI6;
		ddI8 = mymat.ddI8;
		//
		ierr = mymat_1st.D_r_Amplitude(WavinessAmplitude,WavinessLength); CHKERRQ(ierr);
		dI1_r = mymat_1st.dI1;
		dI3_r = mymat_1st.dI3;
		dI5_r = mymat_1st.dI5;
		dI6_r = mymat_1st.dI6;
		dI8_r = mymat_1st.dI8;
		//
		ierr = mymat_1st.D_r_Length(WavinessAmplitude,WavinessLength); CHKERRQ(ierr);
		dI1_s = mymat_1st.dI1;
		dI3_s = mymat_1st.dI3;
		dI5_s = mymat_1st.dI5;
		dI6_s = mymat_1st.dI6;
		dI8_s = mymat_1st.dI8;
      }
      else {
        cout<<"Invalid input of random variable"<<endl;
      }
	  
	  // -----------
	  // 3. calculate the 1st- & 2nd-orderpartial derivative of the transformed compliance matrix
	  // rotate about x-axis
	  // second-order
//      ComplianceMatrix_rs(0,0) = 0;
//      ComplianceMatrix_rs(0,1) = SMat(0,1)*ddI6 + SMat(0,2)*ddI8;
//      ComplianceMatrix_rs(0,2) = SMat(0,2)*ddI6 + SMat(0,1)*ddI8;
//      ComplianceMatrix_rs(1,1) = SMat(0,0)*ddI1 + 2*(SMat(0,0)-SMat(0,1)+SMat(0,2))*ddI3 + SMat(2,2)*ddI5;
//      ComplianceMatrix_rs(1,2) = SMat(0,2)*ddI1 + (SMat(2,2)-SMat(0,0)+2*SMat(0,1))*ddI3 + SMat(0,2)*ddI5;
//      ComplianceMatrix_rs(2,2) = SMat(2,2)*ddI1 + SMat(0,0)*ddI5 + 2*(SMat(0,0)-SMat(0,1)+SMat(0,2))*ddI3;
//      ComplianceMatrix_rs(3,3) = 2*(SMat(0,0)-SMat(0,1))*(ddI1-2*ddI3+ddI5)+4*(SMat(0,0)-2*SMat(0,2)+SMat(2,2))*ddI3;
//      ComplianceMatrix_rs(4,4) = 0;
//      ComplianceMatrix_rs(5,5) = 0;
	  // first-order - 1st variable
//      ComplianceMatrix_r(0,0) = 0;
//      ComplianceMatrix_r(0,1) = SMat(0,1)*dI6_r + SMat(0,2)*dI8_r;
//      ComplianceMatrix_r(0,2) = SMat(0,2)*dI6_r + SMat(0,1)*dI8_r;
//      ComplianceMatrix_r(1,1) = SMat(0,0)*dI1_r + 2*(SMat(0,0)-SMat(0,1)+SMat(0,2))*dI3_r + SMat(2,2)*dI5_r;
//      ComplianceMatrix_r(1,2) = SMat(0,2)*dI1_r + (SMat(2,2)-SMat(0,0)+2*SMat(0,1))*dI3_r + SMat(0,2)*dI5_r;
//      ComplianceMatrix_r(2,2) = SMat(2,2)*dI1_r + SMat(0,0)*dI5_r + 2*(SMat(0,0)-SMat(0,1)+SMat(0,2))*dI3_r;
//      ComplianceMatrix_r(3,3) = 2*(SMat(0,0)-SMat(0,1))*(dI1_r-2*dI3_r+dI5_r)+4*(SMat(0,0)-2*SMat(0,2)+SMat(2,2))*dI3_r;
//      ComplianceMatrix_r(4,4) = 0;
//      ComplianceMatrix_r(5,5) = 0;	
	  // first-order - 2nd variable
//      ComplianceMatrix_s(0,0) = 0;
//      ComplianceMatrix_s(0,1) = SMat(0,1)*dI6_s + SMat(0,2)*dI8_s;
//      ComplianceMatrix_s(0,2) = SMat(0,2)*dI6_s + SMat(0,1)*dI8_s;
//      ComplianceMatrix_s(1,1) = SMat(0,0)*dI1_s + 2*(SMat(0,0)-SMat(0,1)+SMat(0,2))*dI3_s + SMat(2,2)*dI5_s;
//      ComplianceMatrix_s(1,2) = SMat(0,2)*dI1_s + (SMat(2,2)-SMat(0,0)+2*SMat(0,1))*dI3_s + SMat(0,2)*dI5_s;
//      ComplianceMatrix_s(2,2) = SMat(2,2)*dI1_s + SMat(0,0)*dI5_s + 2*(SMat(0,0)-SMat(0,1)+SMat(0,2))*dI3_s;
//      ComplianceMatrix_s(3,3) = 2*(SMat(0,0)-SMat(0,1))*(dI1_s-2*dI3_s+dI5_s)+4*(SMat(0,0)-2*SMat(0,2)+SMat(2,2))*dI3_s;
//      ComplianceMatrix_s(4,4) = 0;
//      ComplianceMatrix_s(5,5) = 0;
	  
	  // rotate about x-axis
	  // second-order
      ComplianceMatrix_rs(0,0) = SMat(0,0)*ddI1 + (2*SMat(0,2)+SMat(5,5))*ddI3 + SMat(2,2)*ddI5;
      ComplianceMatrix_rs(0,1) = SMat(0,1)*ddI6 + SMat(0,2)*ddI8;
      ComplianceMatrix_rs(0,2) = SMat(0,2)*(ddI1+ddI5) + (SMat(0,0)+SMat(2,2)-SMat(5,5))*ddI3;
      ComplianceMatrix_rs(1,1) = 0;
      ComplianceMatrix_rs(1,2) = SMat(0,2)*ddI6 + SMat(0,1)*ddI8;
      ComplianceMatrix_rs(2,2) = SMat(2,2)*ddI1 + SMat(0,0)*ddI5 + (2*SMat(0,2)+SMat(5,5))*ddI3;
      ComplianceMatrix_rs(3,3) = 2*(SMat(0,0)-SMat(0,1))*ddI6 + SMat(5,5)*ddI8;
      ComplianceMatrix_rs(4,4) = 2*(2*SMat(0,0)-4*SMat(0,2)+2*SMat(2,2)-SMat(5,5))*ddI3 + SMat(5,5)*(ddI1+ddI5);
      ComplianceMatrix_rs(5,5) = SMat(5,5)*ddI6 + 2*(SMat(0,0)-SMat(0,1))*ddI8;
	  // first-order - 1st variable
      ComplianceMatrix_r(0,0) = SMat(0,0)*dI1_r + (2*SMat(0,2)+SMat(5,5))*dI3_r + SMat(2,2)*dI5_r;
      ComplianceMatrix_r(0,1) = SMat(0,1)*dI6_r + SMat(0,2)*dI8_r;
      ComplianceMatrix_r(0,2) = SMat(0,2)*(dI1_r+dI5_r) + (SMat(0,0)+SMat(2,2)-SMat(5,5))*dI3_r;
      ComplianceMatrix_r(1,1) = 0;
      ComplianceMatrix_r(1,2) = SMat(0,2)*dI6_r + SMat(0,1)*dI8_r;
      ComplianceMatrix_r(2,2) = SMat(2,2)*dI1_r + SMat(0,0)*dI5_r + (2*SMat(0,2)+SMat(5,5))*dI3_r;
      ComplianceMatrix_r(3,3) = 2*(SMat(0,0)-SMat(0,1))*dI6_r + SMat(5,5)*dI8_r;
      ComplianceMatrix_r(4,4) = 2*(2*SMat(0,0)-4*SMat(0,2)+2*SMat(2,2)-SMat(5,5))*dI3_r + SMat(5,5)*(dI1_r+dI5_r);
      ComplianceMatrix_r(5,5) = SMat(5,5)*dI6_r + 2*(SMat(0,0)-SMat(0,1))*dI8_r;
	  // first-order - 2nd variable
      ComplianceMatrix_s(0,0) = SMat(0,0)*dI1_s + (2*SMat(0,2)+SMat(5,5))*dI3_s + SMat(2,2)*dI5_s;
      ComplianceMatrix_s(0,1) = SMat(0,1)*dI6_s + SMat(0,2)*dI8_s;
      ComplianceMatrix_s(0,2) = SMat(0,2)*(dI1_s+dI5_s) + (SMat(0,0)+SMat(2,2)-SMat(5,5))*dI3_s;
      ComplianceMatrix_s(1,1) = 0;
      ComplianceMatrix_s(1,2) = SMat(0,2)*dI6_s + SMat(0,1)*dI8_s;
      ComplianceMatrix_s(2,2) = SMat(2,2)*dI1_s + SMat(0,0)*dI5_s + (2*SMat(0,2)+SMat(5,5))*dI3_s;
      ComplianceMatrix_s(3,3) = 2*(SMat(0,0)-SMat(0,1))*dI6_s + SMat(5,5)*dI8_s;
      ComplianceMatrix_s(4,4) = 2*(2*SMat(0,0)-4*SMat(0,2)+2*SMat(2,2)-SMat(5,5))*dI3_s + SMat(5,5)*(dI1_s+dI5_s);
      ComplianceMatrix_s(5,5) = SMat(5,5)*dI6_s + 2*(SMat(0,0)-SMat(0,1))*dI8_s;
	  
	  // -----------
      // 4. retrieve transformed stiffness matrix
	  // 4.1 compliance matrix
      double I1, I3, I5, I6, I8;
      double alpha_w;
      alpha_w = 2*M_PI*_WavinessFactor;
      
      I1 = (1+alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);        // m^44
      I3 = (alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);          // m^2 n^2
      I5 = 1-(1+3*alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);    // n^4
      I6 = 1/sqrt(1+alpha_w*alpha_w);                                 // m^2
      I8 = 1-1/sqrt(1+alpha_w*alpha_w);
	  
      ComplianceMatrix(0,0) = SMat(0,0)*I1 + (2*SMat(0,2)+SMat(5,5))*I3 + SMat(2,2)*I5;
      ComplianceMatrix(0,1) = SMat(0,1)*I6 + SMat(0,2)*I8;
      ComplianceMatrix(0,2) = SMat(0,2)*(I1+I5) + (SMat(0,0)+SMat(2,2)-SMat(5,5))*I3;
      ComplianceMatrix(1,1) = SMat(0,0);
      ComplianceMatrix(1,2) = SMat(0,2)*I6 + SMat(0,1)*I8;
      ComplianceMatrix(2,2) = SMat(2,2)*I1 + SMat(0,0)*I5 + (2*SMat(0,2)+SMat(5,5))*I3;
      ComplianceMatrix(3,3) = 2*(SMat(0,0)-SMat(0,1))*I6 + SMat(5,5)*I8;
      ComplianceMatrix(4,4) = 2*(2*SMat(0,0)-4*SMat(0,2)+2*SMat(2,2)-SMat(5,5))*I3 + SMat(5,5)*(I1+I5);
      ComplianceMatrix(5,5) = SMat(5,5)*I6 + 2*(SMat(0,0)-SMat(0,1))*I8;
	  
	  // 4.2 stiffness matrix
	  
	  ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix;
	  ConstitutiveMatrix.resize(6);
	  ConstitutiveMatrix.clear();
	  
      double SVal;
      double S11, S12, S13, S21, S22, S23, S31, S32, S33;
      S11 = ComplianceMatrix(0,0);
      S12 = ComplianceMatrix(0,1);
      S13 = ComplianceMatrix(0,2);
      
      S21 = ComplianceMatrix(1,0);
      S22 = ComplianceMatrix(1,1);
      S23 = ComplianceMatrix(1,2);
      
      S31 = ComplianceMatrix(2,0);
      S32 = ComplianceMatrix(2,1);
      S33 = ComplianceMatrix(2,2);
      
      SVal = S11*(S22 * S33 - S23 * S32) + S12*(S23 * S31 - S21 * S33)
              + S13*(S21 * S32 - S22 * S31);
      
      ConstitutiveMatrix(0,0) = (ComplianceMatrix(1,1)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(1,2),2))/SVal;
      ConstitutiveMatrix(0,1) = (ComplianceMatrix(0,2)*ComplianceMatrix(1,2) - ComplianceMatrix(0,1)*ComplianceMatrix(2,2))/SVal;
      ConstitutiveMatrix(0,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(1,2) - ComplianceMatrix(0,2)*ComplianceMatrix(1,1))/SVal;
      ConstitutiveMatrix(1,1) = (ComplianceMatrix(0,0)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(0,2),2))/SVal;
      ConstitutiveMatrix(1,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(0,2) - ComplianceMatrix(1,2)*ComplianceMatrix(0,0))/SVal;
      ConstitutiveMatrix(2,2) = (ComplianceMatrix(0,0)*ComplianceMatrix(1,1) - pow(ComplianceMatrix(0,1),2))/SVal;
      ConstitutiveMatrix(3,3) = 1/ComplianceMatrix(3,3);
      ConstitutiveMatrix(4,4) = 1/ComplianceMatrix(4,4);
      ConstitutiveMatrix(5,5) = 1/ComplianceMatrix(5,5);
	  
	  // 5. Calculate derivative of transformed stiffness matrix by using the
	  //     formula: d2C/dxdy = C [dS/dy] C [dS/dx] C - C [d2S/dxdy]C  
	  //                       + C [dS/dx] C [dS/dy] C
	  ublas::matrix<double> ConstitutiveMatrix_rs;
	  ConstitutiveMatrix_rs = ublas::zero_matrix<FieldData>(6,6);
	  
	  ublas::matrix<double> DummyMatrix1,DummyMatrix2,DummyMatrix3,DummyMatrix4;
	  DummyMatrix1 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix3 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix4 = ublas::zero_matrix<FieldData>(6,6);
	  ublas::matrix<double> CMat_rs1,CMat_rs2,CMat_rs3;
	  CMat_rs1 = ublas::zero_matrix<FieldData>(6,6);
	  CMat_rs2 = ublas::zero_matrix<FieldData>(6,6);
	  CMat_rs3 = ublas::zero_matrix<FieldData>(6,6);
	  
	  DummyMatrix1 = prod(ComplianceMatrix_r,ConstitutiveMatrix);
	  DummyMatrix2 = prod(ConstitutiveMatrix, DummyMatrix1);
	  DummyMatrix3 = prod(ComplianceMatrix_s, DummyMatrix2);
	  DummyMatrix4 = prod(ConstitutiveMatrix, DummyMatrix3);
	  CMat_rs1 = DummyMatrix4; 
	  
	  DummyMatrix1.clear();DummyMatrix2.clear();
	  DummyMatrix1 = prod(ComplianceMatrix_rs,ConstitutiveMatrix);
	  DummyMatrix2 = prod(ConstitutiveMatrix, DummyMatrix1);
	  CMat_rs2 = DummyMatrix2; 
	
	  DummyMatrix1.clear();DummyMatrix2.clear();
	  DummyMatrix3.clear();DummyMatrix4.clear();
	  DummyMatrix1 = prod(ComplianceMatrix_s,ConstitutiveMatrix);
	  DummyMatrix2 = prod(ConstitutiveMatrix, DummyMatrix1);
	  DummyMatrix3 = prod(ComplianceMatrix_r, DummyMatrix2);
	  DummyMatrix4 = prod(ConstitutiveMatrix, DummyMatrix3);
	  CMat_rs3 = DummyMatrix4;
	  
	  ConstitutiveMatrix_rs  = CMat_rs1 - CMat_rs2 + CMat_rs3;
      //cout<<"Transformed stiffness matrix "<<ConstitutiveMatrix_r<<endl;
	  
	  // -----------
      // 5. Rotate the constitutive matrix
	  D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
      //cout<<material_function<<endl;
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
		   //cout<<DummyMatrix<<endl;
         }
      }
      else if (material_function.compare(0,6,"matrix")==0){
		//cout<<"Matrix \t";
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
											 

	// =========================================================================
	//
	// Calculate the second-order derivative of material constitutive matrix for
	//   yarn with respect to misalignment angle
	//
	// =========================================================================
	virtual PetscErrorCode calculateD_r_Misalignment(double _E_p, double _E_z, 
                                             double _nu_p, double _nu_pz, double _G_zp,
                                             double _lambda, double _mu,
                                             double _vf,
                                             double _theta_f, double _WavinessFactor) {
	  PetscFunctionBegin;
	  // -----------
	  // 1. Get the original stiffness matrix
	  ublas::symmetric_matrix<FieldData,ublas::upper> SMat;
	  SMat.resize(6);
	  SMat.clear();
      YarnStiffnessMatrix TranIsoMat(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,
                                     _lambda, _mu, _vf);
      SMat = TranIsoMat.StiffnessMatrix;
			
	  // -----------
	  // 2. Get parameter due to fibre misalignment
	  double MisalignmentAngle;
	  double dmn, dm2, dn2, dm2n2, dmn3, dm3n, dn4, dm4; // 1st-order derivative
	  MisalignmentAngle = 0.08;   // unit: rad
	  YarnStiffnessMatrix_Geom_FirstOrderDerivative mymat;
      ierr = mymat.D_r_Angle(MisalignmentAngle); CHKERRQ(ierr);
      dmn   = mymat.dmn;
      dm2   = mymat.dm2;
      dn2   = mymat.dn2;
      dm2n2 = mymat.dm2n2;
      dmn3  = mymat.dmn3;
      dm3n  = mymat.dm3n;
      dm4   = mymat.dm4;
      dn4   = mymat.dn4;
			
      // -----------
      // 3. Rotate stiffness matrix according to misalignment angle
      //    about y-axis
      ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix_r;
      ConstitutiveMatrix_r.resize(6);
      ConstitutiveMatrix_r.clear();
			
      ConstitutiveMatrix_r(0,0) = SMat(0,0)*dm4 + SMat(2,2)*dn4
									+ 2*(SMat(0,2)+2*SMat(5,5))*dm2n2;
      ConstitutiveMatrix_r(0,1) = SMat(0,1)*dm2 + SMat(0,2)*dn2;
      ConstitutiveMatrix_r(0,2) = SMat(0,2)*(dm4 + dn4)
									+ (SMat(0,0)+SMat(0,2)-4*SMat(5,5))*dm2n2;
      ConstitutiveMatrix_r(0,4) = dm3n*(SMat(0,2)-SMat(0,0)+2*SMat(5,5)) 															
									+ dmn3*(SMat(2,2)-SMat(0,2)-2*SMat(5,5));
      ConstitutiveMatrix_r(1,2) = SMat(0,2)*dm2 + SMat(0,1)*dn2;
      ConstitutiveMatrix_r(1,4) = (SMat(0,2) - SMat(0,1))*dmn;
      ConstitutiveMatrix_r(2,2) = SMat(2,2)*dm4 + SMat(0,0)*dn4
									+ 2*(SMat(0,2)+2*SMat(5,5))*dm2n2;
      ConstitutiveMatrix_r(2,4) = dm3n*(SMat(2,2)-SMat(0,2)-2*SMat(5,5)) 															
									+ dmn3*(SMat(0,2)-SMat(0,0)+2*SMat(5,5));
      ConstitutiveMatrix_r(3,3) = (SMat(0,0)-SMat(0,1))*dm2/2 + SMat(5,5)*dn2;
      ConstitutiveMatrix_r(3,5) = (SMat(0,0)-SMat(0,1))*dmn/2 - SMat(5,5)*dmn;
      ConstitutiveMatrix_r(4,4) = SMat(5,5)*(dm4+dn4)+
			                (SMat(0,0)-2*SMat(0,2)+SMat(2,2)-2*SMat(5,5))*dm2n2;
      ConstitutiveMatrix_r(5,5) = 	(SMat(0,0)-SMat(0,1))*dn2/2 + SMat(5,5)*dm2;						
			
      // -----------
      // 4. Rotate the constitutive matrix
      D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
       //cout<<material_function<<endl;
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
		   //cout<<DummyMatrix<<endl;
         }
      }
      else if (material_function.compare(0,6,"matrix")==0){
		//cout<<"Matrix \t";
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
		
	virtual PetscErrorCode calculateD_rs_Misalignment(double _E_p, double _E_z, 
                                             double _nu_p, double _nu_pz, double _G_zp,
                                             double _lambda, double _mu,
                                             double _vf,
                                             double _theta_f, double _WavinessFactor) {
      PetscFunctionBegin;
      // -----------
      // 1. Get the original stiffness matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> SMat;
      SMat.resize(6);
      SMat.clear();
      YarnStiffnessMatrix TranIsoMat(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,
                                     _lambda, _mu, _vf);
      SMat = TranIsoMat.StiffnessMatrix;
			
      // -----------
      // 2. Get parameter due to fibre misalignment
      double MisalignmentAngle;
      double ddmn, ddm2, ddn2, ddm2n2, ddmn3, ddm3n, ddn4, ddm4; // 1st-order derivative
      MisalignmentAngle = 0.08;   // unit: rad
      YarnStiffnessMatrix_Geom_SecondOrderDerivative mymat;
      ierr = mymat.D_rs_Misalignment(MisalignmentAngle); CHKERRQ(ierr);
      ddmn   = mymat.ddmn;
      ddm2   = mymat.ddm2;
      ddn2   = mymat.ddn2;
      ddm2n2 = mymat.ddm2n2;
      ddmn3  = mymat.ddmn3;
      ddm3n  = mymat.ddm3n;
      ddm4   = mymat.ddm4;
      ddn4   = mymat.ddn4;
			
      // -----------
      // 3. Rotate stiffness matrix according to misalignment angle
      //    about y-axis
      ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix_rs;
      ConstitutiveMatrix_rs.resize(6);
      ConstitutiveMatrix_rs.clear();
			
      ConstitutiveMatrix_rs(0,0) = SMat(0,0)*ddm4 + SMat(2,2)*ddn4
										+ 2*(SMat(0,2)+2*SMat(5,5))*ddm2n2;
      ConstitutiveMatrix_rs(0,1) = SMat(0,1)*ddm2 + SMat(0,2)*ddn2;
      ConstitutiveMatrix_rs(0,2) = SMat(0,2)*(ddm4 + ddn4)
										+ (SMat(0,0)+SMat(0,2)-4*SMat(5,5))*ddm2n2;
      ConstitutiveMatrix_rs(0,4) = ddm3n*(SMat(0,2)-SMat(0,0)+2*SMat(5,5)) 															
										+ ddmn3*(SMat(2,2)-SMat(0,2)-2*SMat(5,5));
      ConstitutiveMatrix_rs(1,2) = SMat(0,2)*ddm2 + SMat(0,1)*ddn2;
      ConstitutiveMatrix_rs(1,4) = (SMat(0,2) - SMat(0,1))*ddmn;
      ConstitutiveMatrix_rs(2,2) = SMat(2,2)*ddm4 + SMat(0,0)*ddn4
										+ 2*(SMat(0,2)+2*SMat(5,5))*ddm2n2;
      ConstitutiveMatrix_rs(2,4) = ddm3n*(SMat(2,2)-SMat(0,2)-2*SMat(5,5)) 															
										+ ddmn3*(SMat(0,2)-SMat(0,0)+2*SMat(5,5));
      ConstitutiveMatrix_rs(3,3) = (SMat(0,0)-SMat(0,1))*ddm2/2 + SMat(5,5)*ddn2;
      ConstitutiveMatrix_rs(3,5) = (SMat(0,0)-SMat(0,1))*ddmn/2 - SMat(5,5)*ddmn;
      ConstitutiveMatrix_rs(4,4) = SMat(5,5)*(ddm4+ddn4)+
			                (SMat(0,0)-2*SMat(0,2)+SMat(2,2)-2*SMat(5,5))*ddm2n2;
      ConstitutiveMatrix_rs(5,5) = 	(SMat(0,0)-SMat(0,1))*ddn2/2 + SMat(5,5)*ddm2;						
			
      // -----------
      // 4. Rotate the constitutive matrix
      D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
      //cout<<material_function<<endl;
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
		   //cout<<DummyMatrix<<endl;
         }
      }
      else if (material_function.compare(0,6,"matrix")==0){
		//cout<<"Matrix \t";
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
		
	
    // =========================================================================
	//
	// Calculate the second-order derivative of material constitutive matrix for
	//   yarn with respect to fibre volume fraction
	//
	// =========================================================================

	virtual PetscErrorCode calculateD_r_Fraction(double _E_p, double _E_z, 
                                             double _nu_p, double _nu_pz, double _G_zp,
                                             double _lambda, double _mu,
                                             double _vf,
                                             double _theta_f, double _WavinessFactor) {
      PetscFunctionBegin;
	  
      // -----------
	  // 1. Get Hill's moduli
	  //    1.1 tran-iso material
      double _nu_zp=(_nu_pz*_E_z)/_E_p;
      double _G_p = _E_p/(2*(1+_nu_p));
      double k_t, l_t, n_t, m_t, p_t;
      
      k_t = 1/(2*(1-_nu_p)/_E_p-4*_nu_zp*_nu_zp/_E_z);
      l_t = 2*k_t*_nu_zp;
      n_t = _E_z+l_t*l_t/k_t;
      m_t = _G_p;
      p_t = _G_zp; // need to be checked G_21 = G_12 [?]
      
      //    1.2 isotropic material
      double K_value = _lambda+2*_mu/3;
      double G_value = _mu;
      double k_i, l_i, n_i, m_i, p_i;
      
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
      
      double k_f, l_f, n_f, m_f, p_f;            // Hill's moduli for fibre yarn
      double k_m, l_m, n_m, m_m, p_m;                // Hill's moduli for matrix
      k_f=k_t; l_f=l_t; n_f=n_t; m_f=m_t; p_f=p_t;
      k_m=k_i; l_m=l_i; n_m=n_i; m_m=m_i; p_m=p_i;
			
	  // -----------
	  // 2. Calculate partial derivative of the effective stiffness/compliance matrix
	  //  2.1 get stiffness matrix
	  double dk_c, dm_c, dp_c, dl_c, dn_c; // 1st-order fibre volume fraction
	  YarnStiffnessMatrix_Geom_FirstOrderDerivative mymat;
	  ierr = mymat.D_r_Fraction(_vf,k_f, m_f, p_f, l_f, n_f,k_m, m_m, p_m, l_m, n_m); CHKERRQ(ierr);
	  dk_c = mymat.dkc;
	  dm_c = mymat.dmc;
	  dp_c = mymat.dpc;
	  dl_c = mymat.dlc;
	  dn_c = mymat.dnc;
			
	  ublas::symmetric_matrix<FieldData,ublas::upper> CMat_r;
	  CMat_r.resize(6);
	  CMat_r.clear();
	  
//	  // case 1: fibre direction in x-axis
//	  ConstitutiveMatrix_r(0,0) = dn_c;
//	  ConstitutiveMatrix_r(0,1) = ConstitutiveMatrix_r(0,2) = dl_c;
//	  ConstitutiveMatrix_r(1,1) = ConstitutiveMatrix_r(2,2) = dk_c + dm_c;
//	  ConstitutiveMatrix_r(1,2) = dk_c - dm_c;
//	  ConstitutiveMatrix_r(3,3) = dm_c;
//	  ConstitutiveMatrix_r(4,4) = ConstitutiveMatrix_r(5,5) = dp_c;
	  
	  // case 2: fibre direction in z-axis
	  CMat_r(0,0) = CMat_r(1,1) = dk_c + dm_c;
	  CMat_r(0,1) = dk_c - dm_c;
	  CMat_r(0,2) = CMat_r(1,2) = dl_c;
	  CMat_r(2,2) = dn_c;
	  CMat_r(3,3) = dm_c;
	  CMat_r(4,4) = CMat_r(5,5) = dp_c;
	  
	  // 2.2 get original compliance matrix
	  ublas::symmetric_matrix<FieldData,ublas::upper> SMat;
      SMat.resize(6);
      SMat.clear();
      ublas::symmetric_matrix<FieldData,ublas::upper> SMat_r;
      SMat_r.resize(6);
      SMat_r.clear();
	  
      YarnComplianceMatrix TranIsoMat_S(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,
                                       _lambda, _mu, _vf);
      SMat = TranIsoMat_S.ComplianceMatrix;
	  
	  // 2.3 calculate partial derivative of compliance matrix with respect to fibre volume fraction
	  //   using formula: dS/dx = - S [dC/dx] S
	  
	  ublas::matrix<double> DummyMatrix1;
	  DummyMatrix1 = ublas::zero_matrix<FieldData>(6,6);
	  ublas::matrix< FieldData > DummyMatrix2 = prod(CMat_r,SMat);
	  DummyMatrix1 = prod(SMat, DummyMatrix2);
	  SMat_r  = -1*DummyMatrix1;
	  
	  // 3. Calculate transfromed stiffess/compliance matrix due to fibre waviness
     //   3.1 Get waviness parameter
      ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix;
      ComplianceMatrix.resize(6);
      ComplianceMatrix.clear();
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix_r;
      ComplianceMatrix_r.resize(6);
      ComplianceMatrix_r.clear();
      
      double I1, I3, I5, I6, I8;
      double W_mgn, L_wave, alpha_w;
      // W_mgn  =  _amplitude;//1.19;   // unit: mm
      L_wave = 27.9 ;   // unit: mm
      W_mgn  = L_wave*_WavinessFactor;
      
      alpha_w = 2*M_PI*W_mgn/L_wave;
      
      I1 = (1+alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);        // m^44
      I3 = (alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);          // m^2 n^2
      I5 = 1-(1+3*alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);    // n^4
      I6 = 1/sqrt(1+alpha_w*alpha_w);                                 // m^2
      I8 = 1-1/sqrt(1+alpha_w*alpha_w);                               // n^2
      //cout<<"alpha ="<<alpha_w<<"\t I1 = "<<I1<<"\t I3 = "<<I3<<"\t I5 = "<<I5<<"\t I6 = "<<I6<<"\t I8 = "<<I8<<endl;

      // --------
      //   3.2 Calculate transformed compliance matrix
	  // rotate about y-axis
      ComplianceMatrix(0,0) = SMat(0,0)*I1 + (2*SMat(0,2)+SMat(5,5))*I3 + SMat(2,2)*I5;
      ComplianceMatrix(0,1) = SMat(0,1)*I6 + SMat(0,2)*I8;
      ComplianceMatrix(0,2) = SMat(0,2)*(I1+I5) + (SMat(0,0)+SMat(2,2)-SMat(5,5))*I3;
      ComplianceMatrix(1,1) = SMat(0,0);
      ComplianceMatrix(1,2) = SMat(0,2)*I6 + SMat(0,1)*I8;
      ComplianceMatrix(2,2) = SMat(2,2)*I1 + SMat(0,0)*I5 + (2*SMat(0,2)+SMat(5,5))*I3;
      ComplianceMatrix(3,3) = 2*(SMat(0,0)-SMat(0,1))*I6 + SMat(5,5)*I8;
      ComplianceMatrix(4,4) = 2*(2*SMat(0,0)-4*SMat(0,2)+2*SMat(2,2)-SMat(5,5))*I3 + SMat(5,5)*(I1+I5);
      ComplianceMatrix(5,5) = SMat(5,5)*I6 + 2*(SMat(0,0)-SMat(0,1))*I8;

      //   3.3 Calculate partial derivative of the transfomred compliance matrix	
	  
	  ComplianceMatrix_r(0,0) = SMat_r(0,0)*I1 + (2*SMat_r(0,2)+SMat_r(5,5))*I3 + SMat_r(2,2)*I5;
      ComplianceMatrix_r(0,1) = SMat_r(0,1)*I6 + SMat_r(0,2)*I8;
      ComplianceMatrix_r(0,2) = SMat_r(0,2)*(I1+I5) + (SMat_r(0,0)+SMat_r(2,2)-SMat_r(5,5))*I3;
      ComplianceMatrix_r(1,1) = SMat_r(0,0);
      ComplianceMatrix_r(1,2) = SMat_r(0,2)*I6 + SMat_r(0,1)*I8;
      ComplianceMatrix_r(2,2) = SMat_r(2,2)*I1 + SMat_r(0,0)*I5 + (2*SMat_r(0,2)+SMat_r(5,5))*I3;
      ComplianceMatrix_r(3,3) = 2*(SMat_r(0,0)-SMat_r(0,1))*I6 + SMat_r(5,5)*I8;
      ComplianceMatrix_r(4,4) = 2*(2*SMat_r(0,0)-4*SMat_r(0,2)+2*SMat_r(2,2)-SMat_r(5,5))*I3 + SMat_r(5,5)*(I1+I5);
      ComplianceMatrix_r(5,5) = SMat_r(5,5)*I6 + 2*(SMat_r(0,0)-SMat_r(0,1))*I8;
	  
	  //   3.4 Calculate transforemd stiffness matrix
	  ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix;
	  ConstitutiveMatrix.resize(6);
	  ConstitutiveMatrix.clear();
	  
      double SVal;
      double S11, S12, S13, S21, S22, S23, S31, S32, S33;
      S11 = ComplianceMatrix(0,0);
      S12 = ComplianceMatrix(0,1);
      S13 = ComplianceMatrix(0,2);
      
      S21 = ComplianceMatrix(1,0);
      S22 = ComplianceMatrix(1,1);
      S23 = ComplianceMatrix(1,2);
      
      S31 = ComplianceMatrix(2,0);
      S32 = ComplianceMatrix(2,1);
      S33 = ComplianceMatrix(2,2);
      
      SVal = S11*(S22 * S33 - S23 * S32) + S12*(S23 * S31 - S21 * S33)
              + S13*(S21 * S32 - S22 * S31);
      
      ConstitutiveMatrix(0,0) = (ComplianceMatrix(1,1)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(1,2),2))/SVal;
      ConstitutiveMatrix(0,1) = (ComplianceMatrix(0,2)*ComplianceMatrix(1,2) - ComplianceMatrix(0,1)*ComplianceMatrix(2,2))/SVal;
      ConstitutiveMatrix(0,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(1,2) - ComplianceMatrix(0,2)*ComplianceMatrix(1,1))/SVal;
      ConstitutiveMatrix(1,1) = (ComplianceMatrix(0,0)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(0,2),2))/SVal;
      ConstitutiveMatrix(1,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(0,2) - ComplianceMatrix(1,2)*ComplianceMatrix(0,0))/SVal;
      ConstitutiveMatrix(2,2) = (ComplianceMatrix(0,0)*ComplianceMatrix(1,1) - pow(ComplianceMatrix(0,1),2))/SVal;
      ConstitutiveMatrix(3,3) = 1/ComplianceMatrix(3,3);
      ConstitutiveMatrix(4,4) = 1/ComplianceMatrix(4,4);
      ConstitutiveMatrix(5,5) = 1/ComplianceMatrix(5,5);
	  
	  //  3.5. Calculate derivative of transformed stiffness matrix by using the
	  //     formula: dC/dx = - C [dS/dt] C
	  
	  ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix_r;
	  ConstitutiveMatrix_r.resize(6);
	  ConstitutiveMatrix_r.clear();
	  
	  ublas::matrix<double> DummyMatrix3;
	  DummyMatrix3 = ublas::zero_matrix<FieldData>(6,6);
	  ublas::matrix< FieldData > dummyAA = prod(ComplianceMatrix_r,ConstitutiveMatrix);
	  DummyMatrix3 = prod(ConstitutiveMatrix, dummyAA);
	  ConstitutiveMatrix_r  = -1*DummyMatrix3;
			
      // -----------
      // 4. Rotate the constitutive matrix
      D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
       //cout<<material_function<<endl;
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
		   //cout<<DummyMatrix<<endl;
         }
      }
      else if (material_function.compare(0,6,"matrix")==0){
		//cout<<"Matrix \t";
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
		
		
	virtual PetscErrorCode calculateD_rs_Fraction(double _E_p, double _E_z, 
                                             double _nu_p, double _nu_pz, double _G_zp,
                                             double _lambda, double _mu,
                                             double _vf,
                                             double _theta_f, double _WavinessFactor) {
      PetscFunctionBegin;
      // -----------
	  // 1. Get Hill's moduli
	  //    1.1 tran-iso material
      double _nu_zp=(_nu_pz*_E_z)/_E_p;
      double _G_p = _E_p/(2*(1+_nu_p));
      double k_t, l_t, n_t, m_t, p_t;
      
      k_t = 1/(2*(1-_nu_p)/_E_p-4*_nu_zp*_nu_zp/_E_z);
      l_t = 2*k_t*_nu_zp;
      n_t = _E_z+l_t*l_t/k_t;
      m_t = _G_p;
      p_t = _G_zp; // need to be checked G_21 = G_12 [?]
      
      //    1.2 isotropic material
      double K_value = _lambda+2*_mu/3;
      double G_value = _mu;
      double k_i, l_i, n_i, m_i, p_i;
      
      k_i = K_value + G_value/3;
      l_i = K_value - 2*G_value/3;
      n_i = K_value + 4*G_value/3;
      m_i = G_value;
      p_i = G_value;
	  
//	  ublas::symmetric_matrix<FieldData,ublas::upper> CMat_T;
//	  CMat_T.resize(6);
//	  CMat_T.clear();
//	  
//	  CMat_T(0,0) = CMat_T(1,1) = k_t + m_t;
//	  CMat_T(0,1) = k_t - m_t;
//	  CMat_T(0,2) = CMat_T(1,2) = l_t;
//	  CMat_T(2,2) = n_t;
//	  CMat_T(3,3) = m_t;
//	  CMat_T(4,4) = CMat_T(5,5) = p_t;
//	  
//	  cout<<"Transisotropic \t"<<CMat_T<<endl;
      
      /*
       * Calculate Hill's modulus for two-phase composite with fibre of
       * transversely isotropic material and matrix of isotropic material using
       * Mori-Tanaka asymptotic method based on Dvorak derived simple formula
       */
      
      double k_f, l_f, n_f, m_f, p_f;            // Hill's moduli for fibre yarn
      double k_m, l_m, n_m, m_m, p_m;                // Hill's moduli for matrix
      k_f=k_t; l_f=l_t; n_f=n_t; m_f=m_t; p_f=p_t;
      k_m=k_i; l_m=l_i; n_m=n_i; m_m=m_i; p_m=p_i;
			
      // -----------
      // 2. Get parameter due to fibre misalignment
	  double dk_c, dm_c, dp_c, dl_c, dn_c; // 1st-order fibre volume fraction
	  YarnStiffnessMatrix_Geom_FirstOrderDerivative mymat_1st;
	  ierr = mymat_1st.D_r_Fraction(_vf,k_f, m_f, p_f, l_f, n_f,k_m, m_m, p_m, l_m, n_m); CHKERRQ(ierr);
	  dk_c = mymat_1st.dkc;
	  dm_c = mymat_1st.dmc;
	  dp_c = mymat_1st.dpc;
	  dl_c = mymat_1st.dlc;
	  dn_c = mymat_1st.dnc;
	  
      double ddk_c, ddm_c, ddp_c, ddl_c, ddn_c; // 2nd-order fibre volume fraction
      YarnStiffnessMatrix_Geom_SecondOrderDerivative mymat;
      ierr = mymat.D_rs_Fraction(_vf,k_f, m_f, p_f, l_f, n_f,k_m, m_m, p_m, l_m, n_m); CHKERRQ(ierr);
      ddk_c = mymat.ddkc;
      ddm_c = mymat.ddmc;
      ddp_c = mymat.ddpc;
      ddl_c = mymat.ddlc;
      ddn_c = mymat.ddnc;
			
      ublas::symmetric_matrix<FieldData,ublas::upper> CMat_r;
	  CMat_r.resize(6);
	  CMat_r.clear();
	  
      ublas::symmetric_matrix<FieldData,ublas::upper> CMat_rs;
      CMat_rs.resize(6);
      CMat_rs.clear();
	  
//	  // case 1: fibre direction in x-axis
//	  ConstitutiveMatrix_rs(0,0) = ddn_c;
//	  ConstitutiveMatrix_rs(0,1) = ConstitutiveMatrix_rs(0,2) = ddl_c;
//	  ConstitutiveMatrix_rs(1,1) = ConstitutiveMatrix_rs(2,2) = ddk_c + ddm_c;
//	  ConstitutiveMatrix_rs(1,2) = ddk_c - ddm_c;
//	  ConstitutiveMatrix_rs(3,3) = ddm_c;
//	  ConstitutiveMatrix_rs(4,4) = ConstitutiveMatrix_r(5,5) = ddp_c;
	  
	  // case 2: fibre direction in z-axis
	  CMat_r(0,0) = CMat_r(1,1) = dk_c + dm_c;
	  CMat_r(0,1) = dk_c - dm_c;
	  CMat_r(0,2) = CMat_r(1,2) = dl_c;
	  CMat_r(2,2) = dn_c;
	  CMat_r(3,3) = dm_c;
	  CMat_r(4,4) = CMat_r(5,5) = dp_c;
	  
	  CMat_rs(0,0) = CMat_rs(1,1) = ddk_c + ddm_c;
	  CMat_rs(0,1) = ddk_c - ddm_c;
	  CMat_rs(0,2) = CMat_rs(1,2) = ddl_c;
	  CMat_rs(2,2) = ddn_c;
	  CMat_rs(3,3) = ddm_c;
	  CMat_rs(4,4) = CMat_rs(5,5) = ddp_c;	

      // 2.2 get original compliance matrix
	  ublas::symmetric_matrix<FieldData,ublas::upper> SMat;
      SMat.resize(6);
      SMat.clear();
      ublas::symmetric_matrix<FieldData,ublas::upper> SMat_r;
      SMat_r.resize(6);
      SMat_r.clear();
	  ublas::symmetric_matrix<FieldData,ublas::upper> SMat_rs;
      SMat_rs.resize(6);
      SMat_rs.clear();
	  
	  
      YarnComplianceMatrix TranIsoMat_S(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,
                                       _lambda, _mu, _vf);
      SMat = TranIsoMat_S.ComplianceMatrix;
	  
	  // 2.3 calculate partial derivative of compliance matrix with respect to fibre volume fraction
	  //   1st-order using formula: dS/dx = - S [dC/dx] S
	  
	  ublas::matrix<double> DummyMatrix1;
	  DummyMatrix1 = ublas::zero_matrix<FieldData>(6,6);
	  ublas::matrix< FieldData > DummyMatrix2 = prod(CMat_r,SMat);
	  DummyMatrix1 = prod(SMat, DummyMatrix2);
	  SMat_r  = -1*DummyMatrix1;  
	  
	  //   2nd-order using formula: d2C/dxdy = C [dS/dy] C [dS/dx] C - C [d2C/dxdy]C  
	  //                       + C [dS/dx] C [dS/dy] C
	  ublas::matrix<double> DummyMatrix5,DummyMatrix6,DummyMatrix7,DummyMatrix8;
	  DummyMatrix5 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix6 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix7 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix8 = ublas::zero_matrix<FieldData>(6,6);
	  ublas::matrix<double> SMat_rs1,SMat_rs2,SMat_rs3;
	  SMat_rs1 = ublas::zero_matrix<FieldData>(6,6);
	  SMat_rs2 = ublas::zero_matrix<FieldData>(6,6);
	  SMat_rs3 = ublas::zero_matrix<FieldData>(6,6);
	  
	  DummyMatrix5 = prod(CMat_r,SMat);
	  DummyMatrix6 = prod(SMat, DummyMatrix5);
	  DummyMatrix7 = prod(CMat_r, DummyMatrix6);
	  DummyMatrix8 = prod(SMat, DummyMatrix7);
	  SMat_rs1 = DummyMatrix8; 
	  
	  //DummyMatrix5.clear();DummyMatrix6.clear();
	  DummyMatrix5 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix6 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix5 = prod(CMat_rs,SMat);
	  DummyMatrix6 = prod(SMat, DummyMatrix5);
	  SMat_rs2 = DummyMatrix6; 
	
	  //DummyMatrix5.clear();DummyMatrix6.clear();
	  //DummyMatrix7.clear();DummyMatrix8.clear();
	  DummyMatrix5 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix6 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix7 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix8 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix5 = prod(CMat_r,SMat);
	  DummyMatrix6 = prod(SMat, DummyMatrix5);
	  DummyMatrix7 = prod(CMat_r, DummyMatrix6);
	  DummyMatrix8 = prod(SMat, DummyMatrix7);
	  SMat_rs3 = DummyMatrix8;
	  
	  SMat_rs  = SMat_rs1 - SMat_rs2 + SMat_rs3;
	  
	  // 3. Calculate transfromed stiffess/compliance matrix due to fibre waviness
     //   3.1 Get waviness parameter
      ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix;
      ComplianceMatrix.resize(6);
      ComplianceMatrix.clear();
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix_r;
      ComplianceMatrix_r.resize(6);
      ComplianceMatrix_r.clear();
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix_rs;
      ComplianceMatrix_rs.resize(6);
      ComplianceMatrix_rs.clear();
      
      double I1, I3, I5, I6, I8;
      double W_mgn, L_wave, alpha_w;
      // W_mgn  =  _amplitude;//1.19;   // unit: mm
      L_wave = 27.9 ;   // unit: mm
      W_mgn  = L_wave*_WavinessFactor;
      
      alpha_w = 2*M_PI*W_mgn/L_wave;
      
      I1 = (1+alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);        // m^44
      I3 = (alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);          // m^2 n^2
      I5 = 1-(1+3*alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);    // n^4
      I6 = 1/sqrt(1+alpha_w*alpha_w);                                 // m^2
      I8 = 1-1/sqrt(1+alpha_w*alpha_w);                               // n^2
      //cout<<"alpha ="<<alpha_w<<"\t I1 = "<<I1<<"\t I3 = "<<I3<<"\t I5 = "<<I5<<"\t I6 = "<<I6<<"\t I8 = "<<I8<<endl;

      // --------
      //   3.2 Calculate transformed compliance matrix
	  // rotate about y-axis
      ComplianceMatrix(0,0) = SMat(0,0)*I1 + (2*SMat(0,2)+SMat(5,5))*I3 + SMat(2,2)*I5;
      ComplianceMatrix(0,1) = SMat(0,1)*I6 + SMat(0,2)*I8;
      ComplianceMatrix(0,2) = SMat(0,2)*(I1+I5) + (SMat(0,0)+SMat(2,2)-SMat(5,5))*I3;
      ComplianceMatrix(1,1) = SMat(0,0);
      ComplianceMatrix(1,2) = SMat(0,2)*I6 + SMat(0,1)*I8;
      ComplianceMatrix(2,2) = SMat(2,2)*I1 + SMat(0,0)*I5 + (2*SMat(0,2)+SMat(5,5))*I3;
      ComplianceMatrix(3,3) = 2*(SMat(0,0)-SMat(0,1))*I6 + SMat(5,5)*I8;
      ComplianceMatrix(4,4) = 2*(2*SMat(0,0)-4*SMat(0,2)+2*SMat(2,2)-SMat(5,5))*I3 + SMat(5,5)*(I1+I5);
      ComplianceMatrix(5,5) = SMat(5,5)*I6 + 2*(SMat(0,0)-SMat(0,1))*I8;

      //   3.3 Calculate partial derivative of the transfomred compliance matrix	
	  
	  ComplianceMatrix_r(0,0) = SMat_r(0,0)*I1 + (2*SMat_r(0,2)+SMat_r(5,5))*I3 + SMat_r(2,2)*I5;
      ComplianceMatrix_r(0,1) = SMat_r(0,1)*I6 + SMat_r(0,2)*I8;
      ComplianceMatrix_r(0,2) = SMat_r(0,2)*(I1+I5) + (SMat_r(0,0)+SMat_r(2,2)-SMat_r(5,5))*I3;
      ComplianceMatrix_r(1,1) = SMat_r(0,0);
      ComplianceMatrix_r(1,2) = SMat_r(0,2)*I6 + SMat_r(0,1)*I8;
      ComplianceMatrix_r(2,2) = SMat_r(2,2)*I1 + SMat_r(0,0)*I5 + (2*SMat_r(0,2)+SMat_r(5,5))*I3;
      ComplianceMatrix_r(3,3) = 2*(SMat_r(0,0)-SMat_r(0,1))*I6 + SMat_r(5,5)*I8;
      ComplianceMatrix_r(4,4) = 2*(2*SMat_r(0,0)-4*SMat_r(0,2)+2*SMat_r(2,2)-SMat_r(5,5))*I3 + SMat_r(5,5)*(I1+I5);
      ComplianceMatrix_r(5,5) = SMat_r(5,5)*I6 + 2*(SMat_r(0,0)-SMat_r(0,1))*I8;
	  
	  //   3.3 Calculate partial derivative of the transfomred compliance matrix	
	  
	  ComplianceMatrix_rs(0,0) = SMat_rs(0,0)*I1 + (2*SMat_rs(0,2)+SMat_rs(5,5))*I3 + SMat_rs(2,2)*I5;
      ComplianceMatrix_rs(0,1) = SMat_rs(0,1)*I6 + SMat_rs(0,2)*I8;
      ComplianceMatrix_rs(0,2) = SMat_rs(0,2)*(I1+I5) + (SMat_rs(0,0)+SMat_rs(2,2)-SMat_rs(5,5))*I3;
      ComplianceMatrix_rs(1,1) = SMat_rs(0,0);
      ComplianceMatrix_rs(1,2) = SMat_rs(0,2)*I6 + SMat_rs(0,1)*I8;
      ComplianceMatrix_rs(2,2) = SMat_rs(2,2)*I1 + SMat_rs(0,0)*I5 + (2*SMat_rs(0,2)+SMat_rs(5,5))*I3;
      ComplianceMatrix_rs(3,3) = 2*(SMat_rs(0,0)-SMat_rs(0,1))*I6 + SMat_rs(5,5)*I8;
      ComplianceMatrix_rs(4,4) = 2*(2*SMat_rs(0,0)-4*SMat_rs(0,2)+2*SMat_rs(2,2)-SMat_rs(5,5))*I3 + SMat_rs(5,5)*(I1+I5);
      ComplianceMatrix_rs(5,5) = SMat_rs(5,5)*I6 + 2*(SMat_rs(0,0)-SMat_rs(0,1))*I8;
	  
	  //   3.4 Calculate transforemd stiffness matrix
	  ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix;
	  ConstitutiveMatrix.resize(6);
	  ConstitutiveMatrix.clear();
	  
	  ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix_rs;
	  ConstitutiveMatrix_rs.resize(6);
	  ConstitutiveMatrix_rs.clear();
	  
      double SVal;
      double S11, S12, S13, S21, S22, S23, S31, S32, S33;
      S11 = ComplianceMatrix(0,0);
      S12 = ComplianceMatrix(0,1);
      S13 = ComplianceMatrix(0,2);
      
      S21 = ComplianceMatrix(1,0);
      S22 = ComplianceMatrix(1,1);
      S23 = ComplianceMatrix(1,2);
      
      S31 = ComplianceMatrix(2,0);
      S32 = ComplianceMatrix(2,1);
      S33 = ComplianceMatrix(2,2);
      
      SVal = S11*(S22 * S33 - S23 * S32) + S12*(S23 * S31 - S21 * S33)
              + S13*(S21 * S32 - S22 * S31);
      
      ConstitutiveMatrix(0,0) = (ComplianceMatrix(1,1)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(1,2),2))/SVal;
      ConstitutiveMatrix(0,1) = (ComplianceMatrix(0,2)*ComplianceMatrix(1,2) - ComplianceMatrix(0,1)*ComplianceMatrix(2,2))/SVal;
      ConstitutiveMatrix(0,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(1,2) - ComplianceMatrix(0,2)*ComplianceMatrix(1,1))/SVal;
      ConstitutiveMatrix(1,1) = (ComplianceMatrix(0,0)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(0,2),2))/SVal;
      ConstitutiveMatrix(1,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(0,2) - ComplianceMatrix(1,2)*ComplianceMatrix(0,0))/SVal;
      ConstitutiveMatrix(2,2) = (ComplianceMatrix(0,0)*ComplianceMatrix(1,1) - pow(ComplianceMatrix(0,1),2))/SVal;
      ConstitutiveMatrix(3,3) = 1/ComplianceMatrix(3,3);
      ConstitutiveMatrix(4,4) = 1/ComplianceMatrix(4,4);
      ConstitutiveMatrix(5,5) = 1/ComplianceMatrix(5,5);
	  
	  // 5. Calculate derivative of transformed stiffness matrix by using the
	  //     formula: d2C/dxdy = C [dS/dy] C [dS/dx] C - C [d2C/dxdy]C  
	  //                       + C [dS/dx] C [dS/dy] C
	  
	  ublas::matrix<double> DummyMatrix3,DummyMatrix4;
	  DummyMatrix1 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix3 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix4 = ublas::zero_matrix<FieldData>(6,6);
	  ublas::matrix<double> CMat_rs1,CMat_rs2,CMat_rs3;
	  CMat_rs1 = ublas::zero_matrix<FieldData>(6,6);
	  CMat_rs2 = ublas::zero_matrix<FieldData>(6,6);
	  CMat_rs3 = ublas::zero_matrix<FieldData>(6,6);
	  
	  DummyMatrix1 = prod(ComplianceMatrix_r,ConstitutiveMatrix);
	  DummyMatrix2 = prod(ConstitutiveMatrix, DummyMatrix1);
	  DummyMatrix3 = prod(ComplianceMatrix_r, DummyMatrix2);
	  DummyMatrix4 = prod(ConstitutiveMatrix, DummyMatrix3);
	  CMat_rs1 = DummyMatrix4; 
	  
	  //DummyMatrix1.clear();DummyMatrix2.clear();
	  DummyMatrix1 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix1 = prod(ComplianceMatrix_rs,ConstitutiveMatrix);
	  DummyMatrix2 = prod(ConstitutiveMatrix, DummyMatrix1);
	  CMat_rs2 = DummyMatrix2; 
	
	  //DummyMatrix1.clear();DummyMatrix2.clear();
	  //DummyMatrix3.clear();DummyMatrix4.clear();
	  DummyMatrix1 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix3 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix4 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix1 = prod(ComplianceMatrix_r,ConstitutiveMatrix);
	  DummyMatrix2 = prod(ConstitutiveMatrix, DummyMatrix1);
	  DummyMatrix3 = prod(ComplianceMatrix_r, DummyMatrix2);
	  DummyMatrix4 = prod(ConstitutiveMatrix, DummyMatrix3);
	  CMat_rs3 = DummyMatrix4;
	  
	  ConstitutiveMatrix_rs  = CMat_rs1 - CMat_rs2 + CMat_rs3;
			
      // -----------
      // 4. Rotate the constitutive matrix
      D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
      //cout<<material_function<<endl;
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
		   //cout<<DummyMatrix<<endl;
         }
      }
      else if (material_function.compare(0,6,"matrix")==0){
		//cout<<"Matrix \t";
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
		
    // =========================================================================
	//
	// Calculate the second-order derivative of material constitutive matrix for
	//   yarn with respect to material properties
	//
	// =========================================================================
	virtual PetscErrorCode calculateD_r_Material(double _E_p, double _E_z, double _nu_p,
                                             double _nu_pz, double _G_zp,
                                             double _lambda, double _mu,
                                             double _vf,
                                             double _theta_f, double _WavinessFactor) {
	  PetscFunctionBegin;
	  
	  // -----------
	  // 1. Get the original compliance matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> SMat;
      SMat.resize(6);
      SMat.clear();
      
      YarnComplianceMatrix TranIsoMat_S(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,
                                       _lambda, _mu, _vf);
      SMat = TranIsoMat_S.ComplianceMatrix;
	  
	  //  Get waviness parameter
	  double I1, I3, I5, I6, I8;
      double W_mgn, L_wave, alpha_w;
      // W_mgn  =  _amplitude;//1.19;   // unit: mm
      L_wave = 27.9 ;   // unit: mm
      W_mgn  = L_wave*_WavinessFactor;
      
      alpha_w = 2*M_PI*W_mgn/L_wave;
      
      I1 = (1+alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);        // m^44
      I3 = (alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);          // m^2 n^2
      I5 = 1-(1+3*alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);    // n^4
      I6 = 1/sqrt(1+alpha_w*alpha_w);                                 // m^2
      I8 = 1-1/sqrt(1+alpha_w*alpha_w);                               // n^2
	  
	  //  Calculate partial derivative of original stiffness matrix
	  double dk_c, dm_c, dp_c, dl_c, dn_c;        // 1st-order Hill's modulus
	  YarnStiffnessMatrix_Geom_FirstOrderDerivative mymat;
	  double Em, NUm;
	  Em = _mu*(3*_lambda+2*_mu)/(_lambda+_mu);
      NUm = _lambda/2/(_lambda+_mu); 
      if (ix_first_randvar.compare(0,6,"YoungZ") == 0) {
        ierr = mymat.D_r_YoungZ(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ix_first_randvar.compare(0,6,"YoungP") == 0) {
		ierr = mymat.D_r_YoungP(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ix_first_randvar.compare(0,8,"PoissonP") == 0) {
		ierr = mymat.D_r_PoissonP(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ix_first_randvar.compare(0,8,"PoissonZ") == 0) {
		ierr = mymat.D_r_PoissonPZ(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ix_first_randvar.compare(0,6,"YoungM") == 0) {
		ierr = mymat.D_r_Young(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ix_first_randvar.compare(0,8,"PoissonM") == 0) {
		ierr = mymat.D_r_Poisson(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ix_first_randvar.compare(0,7,"ShearZP") == 0) {
		ierr = mymat.D_r_ShearZP(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else {
			cout<<"Invalid input of random variable"<<endl;
	  }
	  
	  dk_c = mymat.dkc;
	  dm_c = mymat.dmc;
	  dp_c = mymat.dpc;
	  dl_c = mymat.dlc;
	  dn_c = mymat.dnc;	  
	  
	  ublas::symmetric_matrix<FieldData,ublas::upper> CMat_r;
	  CMat_r.resize(6);
	  CMat_r.clear();
	  
//	  // case 1: fibre direction in x-axis
//	  ConstitutiveMatrix_r(0,0) = dn_c;
//	  ConstitutiveMatrix_r(0,1) = ConstitutiveMatrix_r(0,2) = dl_c;
//	  ConstitutiveMatrix_r(1,1) = ConstitutiveMatrix_r(2,2) = dk_c + dm_c;
//	  ConstitutiveMatrix_r(1,2) = dk_c - dm_c;
//	  ConstitutiveMatrix_r(3,3) = dm_c;
//	  ConstitutiveMatrix_r(4,4) = ConstitutiveMatrix_r(5,5) = dp_c;
	  
	  // case 2: fibre direction in z-axis
	  CMat_r(0,0) = CMat_r(1,1) = dk_c + dm_c;
	  CMat_r(0,1) = dk_c - dm_c;
	  CMat_r(0,2) = CMat_r(1,2) = dl_c;
	  CMat_r(2,2) = dn_c;
	  CMat_r(3,3) = dm_c;
	  CMat_r(4,4) = CMat_r(5,5) = dp_c;
	  
	  //  Calculate partial derivative of original compliance matrix using formula
	  //           dS/dx = - S [dC/dx] S
	  ublas::symmetric_matrix<FieldData,ublas::upper> SMat_r;
	  SMat_r.resize(6);
	  SMat_r.clear();
	  
	  ublas::matrix<double> DummyMatrix1;
	  DummyMatrix1 = ublas::zero_matrix<FieldData>(6,6);
	  ublas::matrix< FieldData > DummyMatrix2 = prod(CMat_r,SMat);
	  DummyMatrix1 = prod(SMat, DummyMatrix2);
	  SMat_r  = -1*DummyMatrix1;
	  
	  //  Calculate partial derivativ of transformed compliance matrix
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix_r;
      ComplianceMatrix_r.resize(6);
      ComplianceMatrix_r.clear();
	  
	  ComplianceMatrix_r(0,0) = SMat_r(0,0)*I1 + (2*SMat_r(0,2)+SMat_r(5,5))*I3 + SMat_r(2,2)*I5;
      ComplianceMatrix_r(0,1) = SMat_r(0,1)*I6 + SMat_r(0,2)*I8;
      ComplianceMatrix_r(0,2) = SMat_r(0,2)*(I1+I5) + (SMat_r(0,0)+SMat_r(2,2)-SMat_r(5,5))*I3;
      ComplianceMatrix_r(1,1) = SMat_r(0,0);
      ComplianceMatrix_r(1,2) = SMat_r(0,2)*I6 + SMat_r(0,1)*I8;
      ComplianceMatrix_r(2,2) = SMat_r(2,2)*I1 + SMat_r(0,0)*I5 + (2*SMat_r(0,2)+SMat_r(5,5))*I3;
      ComplianceMatrix_r(3,3) = 2*(SMat_r(0,0)-SMat_r(0,1))*I6 + SMat_r(5,5)*I8;
      ComplianceMatrix_r(4,4) = 2*(2*SMat_r(0,0)-4*SMat_r(0,2)+2*SMat_r(2,2)-SMat_r(5,5))*I3 + SMat_r(5,5)*(I1+I5);
      ComplianceMatrix_r(5,5) = SMat_r(5,5)*I6 + 2*(SMat_r(0,0)-SMat_r(0,1))*I8;
	  
	  // Convert partial derivative of transformed compliance to partial derivative
	  // of transfomred stiffness matrix using formula: dC/dx = - C [ds/dx] C
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix;
      ComplianceMatrix.resize(6);
      ComplianceMatrix.clear();
	  
	  // rotate about y-axis
      ComplianceMatrix(0,0) = SMat(0,0)*I1 + (2*SMat(0,2)+SMat(5,5))*I3 + SMat(2,2)*I5;
      ComplianceMatrix(0,1) = SMat(0,1)*I6 + SMat(0,2)*I8;
      ComplianceMatrix(0,2) = SMat(0,2)*(I1+I5) + (SMat(0,0)+SMat(2,2)-SMat(5,5))*I3;
      ComplianceMatrix(1,1) = SMat(0,0);
      ComplianceMatrix(1,2) = SMat(0,2)*I6 + SMat(0,1)*I8;
      ComplianceMatrix(2,2) = SMat(2,2)*I1 + SMat(0,0)*I5 + (2*SMat(0,2)+SMat(5,5))*I3;
      ComplianceMatrix(3,3) = 2*(SMat(0,0)-SMat(0,1))*I6 + SMat(5,5)*I8;
      ComplianceMatrix(4,4) = 2*(2*SMat(0,0)-4*SMat(0,2)+2*SMat(2,2)-SMat(5,5))*I3 + SMat(5,5)*(I1+I5);
      ComplianceMatrix(5,5) = SMat(5,5)*I6 + 2*(SMat(0,0)-SMat(0,1))*I8;
	  
	  //   3.4 Calculate transforemd stiffness matrix
	  ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix;
	  ConstitutiveMatrix.resize(6);
	  ConstitutiveMatrix.clear();
	  
      double SVal;
      double S11, S12, S13, S21, S22, S23, S31, S32, S33;
      S11 = ComplianceMatrix(0,0);
      S12 = ComplianceMatrix(0,1);
      S13 = ComplianceMatrix(0,2);
      
      S21 = ComplianceMatrix(1,0);
      S22 = ComplianceMatrix(1,1);
      S23 = ComplianceMatrix(1,2);
      
      S31 = ComplianceMatrix(2,0);
      S32 = ComplianceMatrix(2,1);
      S33 = ComplianceMatrix(2,2);
      
      SVal = S11*(S22 * S33 - S23 * S32) + S12*(S23 * S31 - S21 * S33)
              + S13*(S21 * S32 - S22 * S31);
      
      ConstitutiveMatrix(0,0) = (ComplianceMatrix(1,1)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(1,2),2))/SVal;
      ConstitutiveMatrix(0,1) = (ComplianceMatrix(0,2)*ComplianceMatrix(1,2) - ComplianceMatrix(0,1)*ComplianceMatrix(2,2))/SVal;
      ConstitutiveMatrix(0,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(1,2) - ComplianceMatrix(0,2)*ComplianceMatrix(1,1))/SVal;
      ConstitutiveMatrix(1,1) = (ComplianceMatrix(0,0)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(0,2),2))/SVal;
      ConstitutiveMatrix(1,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(0,2) - ComplianceMatrix(1,2)*ComplianceMatrix(0,0))/SVal;
      ConstitutiveMatrix(2,2) = (ComplianceMatrix(0,0)*ComplianceMatrix(1,1) - pow(ComplianceMatrix(0,1),2))/SVal;
      ConstitutiveMatrix(3,3) = 1/ComplianceMatrix(3,3);
      ConstitutiveMatrix(4,4) = 1/ComplianceMatrix(4,4);
      ConstitutiveMatrix(5,5) = 1/ComplianceMatrix(5,5);
	  
	  ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix_r;
	  ConstitutiveMatrix_r.resize(6);
	  ConstitutiveMatrix_r.clear();
	  
	  ublas::matrix<double> DummyMatrix3;
	  DummyMatrix3 = ublas::zero_matrix<FieldData>(6,6);
	  ublas::matrix< FieldData > dummyAA = prod(ComplianceMatrix_r,ConstitutiveMatrix);
	  DummyMatrix3 = prod(ConstitutiveMatrix, dummyAA);
	  ConstitutiveMatrix_r  = -1*DummyMatrix3;
	  
	  // -----------
	  // 4. Rotate the constitutive matrix
	  D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
	  //cout<<material_function<<endl;
	  if (material_function.compare(0,13,"reinforcement")==0){
				//cout<<"Reinforcement \t";
      	 vector< ublas::matrix< double > > normalized_phi;
      	 normalized_phi.resize(coords_at_Gauss_nodes.size());
        ierr = ComputeFibreDirection(normalized_phi); CHKERRQ(ierr);
        
      	 for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
           
           int noOfRotations = 1; //Number of Rotations
           
           double zVec[3]={0.0,0.0,1.0};
           double AxVector[3]={normalized_phi[gg](0,1)*zVec[2]-normalized_phi[gg](0,2)*zVec[1],
                               normalized_phi[gg](0,2)*zVec[0]-normalized_phi[gg](0,0)*zVec[2],
							    normalized_phi[gg](0,0)*zVec[1]-normalized_phi[gg](0,1)*zVec[0]};
           double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))
		                      /(sqrt(pow(normalized_phi[gg](0,0),2)+pow(normalized_phi[gg](0,1),2)
							  +pow(normalized_phi[gg](0,2),2)))
							  *(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};
           
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
		   //cout<<DummyMatrix<<endl;
         }
      }
      else if (material_function.compare(0,6,"matrix")==0){
		//cout<<"Matrix \t";
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
	
	virtual PetscErrorCode calculateD_rs_Material(double _E_p, double _E_z, double _nu_p,
                                             double _nu_pz, double _G_zp,
                                             double _lambda, double _mu,
                                             double _vf,
                                             double _theta_f, double _WavinessFactor) {
	  PetscFunctionBegin;
	  
	  // -----------
	  // 1. Get the original compliance matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> SMat;
      SMat.resize(6);
      SMat.clear();
      
      YarnComplianceMatrix TranIsoMat_S(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,
                                       _lambda, _mu, _vf);
      SMat = TranIsoMat_S.ComplianceMatrix;
	  
	  // 2. Calculate partial derivative of original stiffness matrix
	  //  2.1 Calculate 1st-order partial derivative of original stiffness matrix
	  double dk_c, dm_c, dp_c, dl_c, dn_c;        // 1st-order Hill's modulus
	  YarnStiffnessMatrix_Geom_FirstOrderDerivative mymat;
	  double Em, NUm;
	  Em = _mu*(3*_lambda+2*_mu)/(_lambda+_mu);
      NUm = _lambda/2/(_lambda+_mu); //cout<<"Modulus: \t"<<Em<<"Poisson: \t"<<NUm<<endl;
	  if (ix_first_randvar.compare(0,6,"YoungZ") == 0) {
        ierr = mymat.D_r_YoungZ(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ix_first_randvar.compare(0,6,"YoungP") == 0) {
		ierr = mymat.D_r_YoungP(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ix_first_randvar.compare(0,8,"PoissonP") == 0) {
		ierr = mymat.D_r_PoissonP(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ix_first_randvar.compare(0,8,"PoissonZ") == 0) {
		ierr = mymat.D_r_PoissonPZ(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ix_first_randvar.compare(0,6,"YoungM") == 0) {
		ierr = mymat.D_r_Young(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ix_first_randvar.compare(0,8,"PoissonM") == 0) {
		ierr = mymat.D_r_Poisson(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ix_first_randvar.compare(0,7,"ShearZP") == 0) {
		ierr = mymat.D_r_ShearZP(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else {
			cout<<"Invalid input of random variable"<<endl;
	  }
      
	  dk_c = mymat.dkc;
	  dm_c = mymat.dmc;
	  dp_c = mymat.dpc;
	  dl_c = mymat.dlc;
	  dn_c = mymat.dnc;	
	  
	  ublas::symmetric_matrix<FieldData,ublas::upper> CMat_r;
	  CMat_r.resize(6);
	  CMat_r.clear();
	  
//	  // case 1: fibre direction in x-axis
//	  ConstitutiveMatrix_r(0,0) = dn_c;
//	  ConstitutiveMatrix_r(0,1) = ConstitutiveMatrix_r(0,2) = dl_c;
//	  ConstitutiveMatrix_r(1,1) = ConstitutiveMatrix_r(2,2) = dk_c + dm_c;
//	  ConstitutiveMatrix_r(1,2) = dk_c - dm_c;
//	  ConstitutiveMatrix_r(3,3) = dm_c;
//	  ConstitutiveMatrix_r(4,4) = ConstitutiveMatrix_r(5,5) = dp_c;
	  
	  // case 2: fibre direction in z-axis
	  CMat_r(0,0) = CMat_r(1,1) = dk_c + dm_c;
	  CMat_r(0,1) = dk_c - dm_c;
	  CMat_r(0,2) = CMat_r(1,2) = dl_c;
	  CMat_r(2,2) = dn_c;
	  CMat_r(3,3) = dm_c;
	  CMat_r(4,4) = CMat_r(5,5) = dp_c;
	  
	  //  2.2 Calculate 2nd-order partial derivative of original stiffness matrix
	  double ddk_c, ddm_c, ddp_c, ddl_c, ddn_c;   // 2nd-order Hill's modulus
	  YarnStiffnessMatrix_Geom_SecondOrderDerivative mymat_2nd;
	  if ((ix_first_randvar.compare(0,6,"YoungZ") == 0) && (ix_second_randvar.compare(0,6,"YoungZ")==0)) {
        ierr = mymat_2nd.D_rs_YoungZ(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if ((ix_first_randvar.compare(0,6,"YoungP") == 0) && (ix_second_randvar.compare(0,6,"YoungP")==0)) {
		ierr = mymat_2nd.D_rs_YoungP(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if ((ix_first_randvar.compare(0,8,"PoissonP") == 0) && (ix_second_randvar.compare(0,8,"PoissonP")==0)) {
		ierr = mymat_2nd.D_rs_PoissonP(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if ((ix_first_randvar.compare(0,8,"PoissonZ") == 0) && (ix_second_randvar.compare(0,8,"PoissonZ")==0)) {
		ierr = mymat_2nd.D_rs_PoissonPZ(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if ((ix_first_randvar.compare(0,6,"YoungM") == 0) && (ix_second_randvar.compare(0,6,"YoungM")==0)) {
		ierr = mymat_2nd.D_rs_Young(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if ((ix_first_randvar.compare(0,8,"PoissonM") == 0) && (ix_second_randvar.compare(0,8,"PoissonM")==0)) {
		ierr = mymat_2nd.D_rs_Poisson(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if ((ix_first_randvar.compare(0,7,"ShearZP") == 0) && (ix_second_randvar.compare(0,7,"ShearZP")==0)) {
		ierr = mymat_2nd.D_rs_ShearZP(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else {
		cout<<"Invalid input of random variable"<<endl;
	  }
	  
	  ddk_c = mymat_2nd.ddkc;
	  ddm_c = mymat_2nd.ddmc;
	  ddp_c = mymat_2nd.ddpc;
	  ddl_c = mymat_2nd.ddlc;
	  ddn_c = mymat_2nd.ddnc;
	  
	  ublas::symmetric_matrix<FieldData,ublas::upper> CMat_rs;
	  CMat_rs.resize(6);
	  CMat_rs.clear();	
	  
//	  // case 1: fibre direction in x-axis
//	  CMat_rs(0,0) = ddn_c;
//	  CMat_rs(0,1) = CMat_rs(0,2) = ddl_c;
//	  CMat_rs(1,1) = CMat_rs(2,2) = ddk_c + ddm_c;
//	  CMat_rs(1,2) = ddk_c - ddm_c;
//	  CMat_rs(3,3) = ddm_c;
//	  CMat_rs(4,4) = CMat_rs(5,5) = ddp_c;
	  
	  // case 2: fibre direction in z-axis
	  CMat_rs(0,0) = CMat_rs(1,1) = ddk_c + ddm_c;
	  CMat_rs(0,1) = ddk_c - ddm_c;
	  CMat_rs(0,2) = CMat_rs(1,2) = ddl_c;
	  CMat_rs(2,2) = ddn_c;
	  CMat_rs(3,3) = ddm_c;
	  CMat_rs(4,4) = CMat_rs(5,5) = ddp_c;	
//	  if (ix_first_randvar.compare(0,8,"PoissonP") == 0) {
//	    cout<<"CMat_rs = \t"<<CMat_rs<<endl; }
	  
	  //  2.3 Calculate 1st-order partial derivative of original compliance 
	  //        matrix using formula dS/dx = - S [dC/dx] S
	  // ublas::symmetric_matrix<FieldData,ublas::upper> SMat_r;
	  // SMat_r.resize(6);
	  // SMat_r.clear();
	  ublas::matrix<double> SMat_r;
	  SMat_r = ublas::zero_matrix<FieldData>(6,6);
	  
	  ublas::matrix<double> DummyMatrix1;
	  DummyMatrix1 = ublas::zero_matrix<FieldData>(6,6);
	  ublas::matrix< FieldData > DummyMatrix2 = prod(CMat_r,SMat);
	  DummyMatrix1 = prod(SMat, DummyMatrix2);
	  SMat_r  = -1*DummyMatrix1;

//	  if (ix_first_randvar.compare(0,8,"PoissonP") == 0) {
//	    cout<<"SMat_r= \t"<<SMat_r<<endl; }
	  
	  //  2.4 Calculate 2nd-order partial derivative of original compliance 
	  //        matrix using formula: d2S/dxdy = S [dC/dy] S [dC/dx] S 
	  //                      - S [d2C/dxdy]S + S [dC/dx] S [dC/dy] S
	  // ublas::symmetric_matrix<FieldData,ublas::upper> SMat_rs;
	  // SMat_rs.resize(6);
	  // SMat_rs.clear();
	  ublas::matrix<double> SMat_rs;
	  SMat_rs = ublas::zero_matrix<FieldData>(6,6);
	  
	  ublas::matrix<double> DummyMatrix5,DummyMatrix6,DummyMatrix7,DummyMatrix8;
	  DummyMatrix5 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix6 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix7 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix8 = ublas::zero_matrix<FieldData>(6,6);
	  ublas::matrix<double> SMat_rs1,SMat_rs2,SMat_rs3;
	  SMat_rs1 = ublas::zero_matrix<FieldData>(6,6);
	  SMat_rs2 = ublas::zero_matrix<FieldData>(6,6);
	  SMat_rs3 = ublas::zero_matrix<FieldData>(6,6);
	  
	  DummyMatrix5 = prod(CMat_r,SMat);
	  DummyMatrix6 = prod(SMat, DummyMatrix5);
	  DummyMatrix7 = prod(CMat_r, DummyMatrix6);
	  DummyMatrix8 = prod(SMat, DummyMatrix7);
	  SMat_rs1 = DummyMatrix8; 
	  
//	 if (ix_first_randvar.compare(0,8,"PoissonP") == 0) {
//	    cout<<"SMat_rs1= \t"<<SMat_rs1<<endl; }
	  
	  //DummyMatrix5.clear();DummyMatrix6.clear();
	  DummyMatrix5 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix6 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix5 = prod(CMat_rs,SMat);
	  DummyMatrix6 = prod(SMat, DummyMatrix5);
	  SMat_rs2 = DummyMatrix6; 
	  
//	  if (ix_first_randvar.compare(0,8,"PoissonP") == 0) {
//	    cout<<"SMat_rs2= \t"<<SMat_rs2<<endl; }
	
	  //DummyMatrix5.clear();DummyMatrix6.clear();
	  //DummyMatrix7.clear();DummyMatrix8.clear();
	  DummyMatrix5 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix6 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix7 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix8 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix5 = prod(CMat_r,SMat);
	  DummyMatrix6 = prod(SMat, DummyMatrix5);
	  DummyMatrix7 = prod(CMat_r, DummyMatrix6);
	  DummyMatrix8 = prod(SMat, DummyMatrix7);
	  SMat_rs3 = DummyMatrix8;
	  
//	  if (ix_first_randvar.compare(0,8,"PoissonP") == 0) {
//	    cout<<"SMat_rs3= \t"<<SMat_rs3<<endl; }
	  
	  SMat_rs  = SMat_rs1 - SMat_rs2 + SMat_rs3;
	  
//	  if (ix_first_randvar.compare(0,8,"PoissonP") == 0) {
//	    cout<<"SMat_rs= \t"<<SMat_rs<<endl; }
	  
	  // 3. Calculated 2nd-order partial derivative of transformed constitutive matrix
	  //  3.1. Get waviness parameter
	  double I1, I3, I5, I6, I8;
      double W_mgn, L_wave, alpha_w;
      // W_mg<< = 27.9 ;   // unit: mm
	  L_wave = 27.9;
      W_mgn  = L_wave*_WavinessFactor;
      
      alpha_w = 2*M_PI*W_mgn/L_wave;
      
      I1 = (1+alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);        // m^44
      I3 = (alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);          // m^2 n^2
      I5 = 1-(1+3*alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);    // n^4
      I6 = 1/sqrt(1+alpha_w*alpha_w);                                 // m^2
      I8 = 1-1/sqrt(1+alpha_w*alpha_w);                               // n^2
	  
	  //  3.2 Calculate 1st-order partial derivative of transformed compliance matrix
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix_r;
	  ComplianceMatrix_r.resize(6);
	  ComplianceMatrix_r.clear();
	  
	  ComplianceMatrix_r(0,0) = SMat_r(0,0)*I1 + (2*SMat_r(0,2)+SMat_r(5,5))*I3 + SMat_r(2,2)*I5;
	  ComplianceMatrix_r(0,1) = SMat_r(0,1)*I6 + SMat_r(0,2)*I8;
	  ComplianceMatrix_r(0,2) = SMat_r(0,2)*(I1+I5) + (SMat_r(0,0)+SMat_r(2,2)-SMat_r(5,5))*I3;
	  ComplianceMatrix_r(1,1) = SMat_r(0,0);
	  ComplianceMatrix_r(1,2) = SMat_r(0,2)*I6 + SMat_r(0,1)*I8;
	  ComplianceMatrix_r(2,2) = SMat_r(2,2)*I1 + SMat_r(0,0)*I5 + (2*SMat_r(0,2)+SMat_r(5,5))*I3;
	  ComplianceMatrix_r(3,3) = 2*(SMat_r(0,0)-SMat_r(0,1))*I6 + SMat_r(5,5)*I8;
	  ComplianceMatrix_r(4,4) = 2*(2*SMat_r(0,0)-4*SMat_r(0,2)+2*SMat_r(2,2)-SMat_r(5,5))*I3 + SMat_r(5,5)*(I1+I5);
	  ComplianceMatrix_r(5,5) = SMat_r(5,5)*I6 + 2*(SMat_r(0,0)-SMat_r(0,1))*I8;
	  
	  // cout<<"Transformed SMat_r = \t"<<ComplianceMatrix_r<<endl;
	  
	  //  3.3 Calculate 2nd-order partial derivative of transformed compliance matrix
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix_rs;
	  ComplianceMatrix_rs.resize(6);
	  ComplianceMatrix_rs.clear();
	  
	  ComplianceMatrix_rs(0,0) = SMat_rs(0,0)*I1 + (2*SMat_rs(0,2)+SMat_rs(5,5))*I3 + SMat_rs(2,2)*I5;
	  ComplianceMatrix_rs(0,1) = SMat_rs(0,1)*I6 + SMat_rs(0,2)*I8;
	  ComplianceMatrix_rs(0,2) = SMat_rs(0,2)*(I1+I5) + (SMat_rs(0,0)+SMat_rs(2,2)-SMat_rs(5,5))*I3;
	  ComplianceMatrix_rs(1,1) = SMat_rs(0,0);
	  ComplianceMatrix_rs(1,2) = SMat_rs(0,2)*I6 + SMat_rs(0,1)*I8;
	  ComplianceMatrix_rs(2,2) = SMat_rs(2,2)*I1 + SMat_rs(0,0)*I5 + (2*SMat_rs(0,2)+SMat_rs(5,5))*I3;
	  ComplianceMatrix_rs(3,3) = 2*(SMat_rs(0,0)-SMat_rs(0,1))*I6 + SMat_rs(5,5)*I8;
	  ComplianceMatrix_rs(4,4) = 2*(2*SMat_rs(0,0)-4*SMat_rs(0,2)+2*SMat_rs(2,2)-SMat_rs(5,5))*I3 + SMat_rs(5,5)*(I1+I5);
	  ComplianceMatrix_rs(5,5) = SMat_rs(5,5)*I6 + 2*(SMat_rs(0,0)-SMat_rs(0,1))*I8;	  
	  
	  // cout<<"Transformed SMat_rs = \t"<<ComplianceMatrix_rs<<endl;
	  
	  // 3.4 Convert partial derivative of transformed compliance to partial derivative
	  // of transfomred stiffness matrix
	  ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix;
      ComplianceMatrix.resize(6);
      ComplianceMatrix.clear();
	  
	  // rotate about y-axis
      ComplianceMatrix(0,0) = SMat(0,0)*I1 + (2*SMat(0,2)+SMat(5,5))*I3 + SMat(2,2)*I5;
      ComplianceMatrix(0,1) = SMat(0,1)*I6 + SMat(0,2)*I8;
      ComplianceMatrix(0,2) = SMat(0,2)*(I1+I5) + (SMat(0,0)+SMat(2,2)-SMat(5,5))*I3;
      ComplianceMatrix(1,1) = SMat(0,0);
      ComplianceMatrix(1,2) = SMat(0,2)*I6 + SMat(0,1)*I8;
      ComplianceMatrix(2,2) = SMat(2,2)*I1 + SMat(0,0)*I5 + (2*SMat(0,2)+SMat(5,5))*I3;
      ComplianceMatrix(3,3) = 2*(SMat(0,0)-SMat(0,1))*I6 + SMat(5,5)*I8;
      ComplianceMatrix(4,4) = 2*(2*SMat(0,0)-4*SMat(0,2)+2*SMat(2,2)-SMat(5,5))*I3 + SMat(5,5)*(I1+I5);
      ComplianceMatrix(5,5) = SMat(5,5)*I6 + 2*(SMat(0,0)-SMat(0,1))*I8;
	  
	  // cout<<"Transformed SMat = \t"<<ComplianceMatrix<<endl;
	  
	  //   3.4.1 Calculate transforemd stiffness matrix
	  ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix;
	  ConstitutiveMatrix.resize(6);
	  ConstitutiveMatrix.clear();
	  
      double SVal;
      double S11, S12, S13, S21, S22, S23, S31, S32, S33;
      S11 = ComplianceMatrix(0,0);
      S12 = ComplianceMatrix(0,1);
      S13 = ComplianceMatrix(0,2);
      
      S21 = ComplianceMatrix(1,0);
      S22 = ComplianceMatrix(1,1);
      S23 = ComplianceMatrix(1,2);
      
      S31 = ComplianceMatrix(2,0);
      S32 = ComplianceMatrix(2,1);
      S33 = ComplianceMatrix(2,2);
      
      SVal = S11*(S22 * S33 - S23 * S32) + S12*(S23 * S31 - S21 * S33)
              + S13*(S21 * S32 - S22 * S31);
      
      ConstitutiveMatrix(0,0) = (ComplianceMatrix(1,1)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(1,2),2))/SVal;
      ConstitutiveMatrix(0,1) = (ComplianceMatrix(0,2)*ComplianceMatrix(1,2) - ComplianceMatrix(0,1)*ComplianceMatrix(2,2))/SVal;
      ConstitutiveMatrix(0,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(1,2) - ComplianceMatrix(0,2)*ComplianceMatrix(1,1))/SVal;
      ConstitutiveMatrix(1,1) = (ComplianceMatrix(0,0)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(0,2),2))/SVal;
      ConstitutiveMatrix(1,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(0,2) - ComplianceMatrix(1,2)*ComplianceMatrix(0,0))/SVal;
      ConstitutiveMatrix(2,2) = (ComplianceMatrix(0,0)*ComplianceMatrix(1,1) - pow(ComplianceMatrix(0,1),2))/SVal;
      ConstitutiveMatrix(3,3) = 1/ComplianceMatrix(3,3);
      ConstitutiveMatrix(4,4) = 1/ComplianceMatrix(4,4);
      ConstitutiveMatrix(5,5) = 1/ComplianceMatrix(5,5);
	  
	  // cout<<"Transformed CMat = \t"<<ConstitutiveMatrix<<endl;
//	  if (ix_first_randvar.compare(0,8,"PoissonP") == 0) {
//		cout<<"PI = \t"<<M
//		cout<<"I1 = \t"<<I1<<"\t I3 = \t"<<I3<<"\t I5 = \t"<<I5<<"\t I6 = \t"<<I6<<"\t I8 = \t"<<I8<<endl;
//		cout<<"original SMat = \t"<<SMat<<endl;
//		cout<<"transformed SMat = \t"<<ComplianceMatrix<<endl;
//	    cout<<"CMat= \t"<<ConstitutiveMatrix<<endl; }	  
	  
	 //   3.4.1 Calculate derivative of transformed stiffness matrix by using the
	  //         formula: d2C/dxdy = C [dS/dy] C [dS/dx] C - C [d2C/dxdy]C  
	  //                       + C [dS/dx] C [dS/dy] C
	  ublas::matrix<double> ConstitutiveMatrix_rs;
	  ConstitutiveMatrix_rs = ublas::zero_matrix<FieldData>(6,6);
	  
	  ublas::matrix<double> DummyMatrix3,DummyMatrix4;
	  DummyMatrix1 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix3 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix4 = ublas::zero_matrix<FieldData>(6,6);
	  ublas::matrix<double> CMat_rs1,CMat_rs2,CMat_rs3;
	  CMat_rs1 = ublas::zero_matrix<FieldData>(6,6);
	  CMat_rs2 = ublas::zero_matrix<FieldData>(6,6);
	  CMat_rs3 = ublas::zero_matrix<FieldData>(6,6);
	  
	  DummyMatrix1 = prod(ComplianceMatrix_r,ConstitutiveMatrix);
	  DummyMatrix2 = prod(ConstitutiveMatrix, DummyMatrix1);
	  DummyMatrix3 = prod(ComplianceMatrix_r, DummyMatrix2);
	  DummyMatrix4 = prod(ConstitutiveMatrix, DummyMatrix3);
	  CMat_rs1 = DummyMatrix4; 
//	  if (ix_first_randvar.compare(0,8,"PoissonP") == 0) {
//	    cout<<"CMat_rs1= \t"<<CMat_rs1<<endl; }	  
	  
	  //DummyMatrix1.clear();DummyMatrix2.clear();
	  DummyMatrix1 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix1 = prod(ComplianceMatrix_rs,ConstitutiveMatrix);
	  DummyMatrix2 = prod(ConstitutiveMatrix, DummyMatrix1);
	  CMat_rs2 = DummyMatrix2; 
//	  if (ix_first_randvar.compare(0,8,"PoissonP") == 0) {
//	    cout<<"CMat_rs2= \t"<<CMat_rs2<<endl; }	  
	
	  //DummyMatrix1.clear();DummyMatrix2.clear();
	  //DummyMatrix3.clear();DummyMatrix4.clear();
	  DummyMatrix1 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix3 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix4 = ublas::zero_matrix<FieldData>(6,6);
	  DummyMatrix1 = prod(ComplianceMatrix_r,ConstitutiveMatrix);
	  DummyMatrix2 = prod(ConstitutiveMatrix, DummyMatrix1);
	  DummyMatrix3 = prod(ComplianceMatrix_r, DummyMatrix2);
	  DummyMatrix4 = prod(ConstitutiveMatrix, DummyMatrix3);
	  CMat_rs3 = DummyMatrix4;
//	  if (ix_first_randvar.compare(0,8,"PoissonP") == 0) {
//	    cout<<"CMat_rs3= \t"<<CMat_rs3<<endl; }		  
	  
	  ConstitutiveMatrix_rs  = CMat_rs1 - CMat_rs2 + CMat_rs3;
	  
	  // cout<<"Transformed CMat_rs = \t"<<ConstitutiveMatrix_rs<<endl;

//	  if (ix_first_randvar.compare(0,8,"PoissonP") == 0) {
//	    cout<<"CMat_rs= \t"<<ConstitutiveMatrix_rs<<endl; }	  
	  
	  // -----------
	  // 4. Rotate the constitutive matrix
	  D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
	  //cout<<material_function<<endl;
	  if (material_function.compare(0,13,"reinforcement")==0){
				//cout<<"Reinforcement \t";
      	 vector< ublas::matrix< double > > normalized_phi;
      	 normalized_phi.resize(coords_at_Gauss_nodes.size());
        ierr = ComputeFibreDirection(normalized_phi); CHKERRQ(ierr);
        
      	 for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
           
           int noOfRotations = 1; //Number of Rotations
           
           double zVec[3]={0.0,0.0,1.0};
           double AxVector[3]={normalized_phi[gg](0,1)*zVec[2]-normalized_phi[gg](0,2)*zVec[1],
                               normalized_phi[gg](0,2)*zVec[0]-normalized_phi[gg](0,0)*zVec[2],
							    normalized_phi[gg](0,0)*zVec[1]-normalized_phi[gg](0,1)*zVec[0]};
           double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))
		                      /(sqrt(pow(normalized_phi[gg](0,0),2)+pow(normalized_phi[gg](0,1),2)
							  +pow(normalized_phi[gg](0,2),2)))
							  *(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};
           
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
		   //cout<<DummyMatrix<<endl;
         }
      }
      else if (material_function.compare(0,6,"matrix")==0){
		//cout<<"Matrix \t";
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
	  // -----------
      // 1. Get material parameters
      double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
      double _lambda, _mu;
      double _vf;
      double _theta_f;
      double _WavinessFactor;
      ierr = GetMatParameters(&_E_p, &_E_z, &_nu_p, &_nu_pz, &_G_zp,
                              &_lambda, &_mu, &_vf, &_theta_f, &_WavinessFactor); CHKERRQ(ierr);
      // -----------
      // 2. Calculate material constitutive matrix
      if ((ix_first_randvar.compare(0,9,"Amplitude") == 0) || (ix_first_randvar.compare(0,6,"Length") == 0)) {
        ierr = calculateD_r_Waviness(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                                 _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);
      }
      else if (ix_first_randvar.compare(0,5,"Angle") == 0) {
        ierr = calculateD_r_Misalignment(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                                         _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);
      }
      else if (ix_first_randvar.compare(0,8,"Fraction") == 0) {
        ierr = calculateD_r_Fraction(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                                     _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);
      }
	  else if ((ix_first_randvar.compare(0,6,"YoungZ") == 0) ||
			   (ix_first_randvar.compare(0,6,"YoungP") == 0) ||
			   (ix_first_randvar.compare(0,8,"PoissonP") == 0) ||
			   (ix_first_randvar.compare(0,8,"PoissonZ") == 0) ||
			   (ix_first_randvar.compare(0,7,"ShearZP") == 0) ||
			   (ix_first_randvar.compare(0,6,"YoungM") == 0) ||
			   (ix_first_randvar.compare(0,8,"PoissonM") == 0)) {
        ierr = calculateD_r_Material(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                                         _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);
      }
      else {
        //cout<<"Invalid input of random variable"<<endl;
      }
      
	  // -----------
      // 3. Calculate element stiffness matrix					   
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
	  // -----------
      // 1. Get material parameters
      double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
      double _lambda, _mu;
      double _vf;
      double _theta_f;
      double _WavinessFactor;
      ierr = GetMatParameters(&_E_p, &_E_z, &_nu_p, &_nu_pz, &_G_zp,
                              &_lambda, &_mu, &_vf, &_theta_f, &_WavinessFactor); CHKERRQ(ierr);
      // -----------
      // 2. Calculate material constitutive matrix
      
      if ((ix_first_randvar.compare(0,9,"Amplitude") == 0) || (ix_first_randvar.compare(0,6,"Length") == 0)) {
        ierr = calculateD_rs_Waviness(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                                 _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);
      }
      else if (ix_first_randvar.compare(0,5,"Angle") == 0) {
        ierr = calculateD_rs_Misalignment(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                                         _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);
      }
      else if (ix_first_randvar.compare(0,8,"Fraction") == 0) {
        ierr = calculateD_rs_Fraction(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                                     _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);
      }
	  else if ((ix_first_randvar.compare(0,6,"YoungZ") == 0)  ||
			   (ix_first_randvar.compare(0,6,"YoungP") == 0)   ||
			   (ix_first_randvar.compare(0,8,"PoissonP") == 0) ||
			   (ix_first_randvar.compare(0,8,"PoissonZ") == 0) ||
			   (ix_first_randvar.compare(0,7,"ShearZP") == 0)  ||
			   (ix_first_randvar.compare(0,6,"YoungM") == 0)   ||
			   (ix_first_randvar.compare(0,8,"PoissonM") == 0)) {
        ierr = calculateD_rs_Material(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
				                       _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);
      }
      else {
        //cout<<"Invalid input of random variable"<<endl;
      }
      
	  // -----------
      // 3. Calculate element stiffness matrix		
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
	
	// ======================================================================
    //
    // Calculate element force to establish RHS of FE equilibrium equation
    //
	// F: calculate the second-order partial derivative of "external force" 
	//    which is referred as right-hand side in the algebriac equation 
	//    K][U_rs] = [F_rs]
	//
	//    [F_rs] = - [K_rs][U] - 2[K_r][U_s]
    // ======================================================================
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
	  ierr = GetDataVector(zeroth_field,D_elm[col_mat1]); CHKERRQ(ierr);
	  ierr = GetDataVector(first_field,D_elm_r[col_mat1]); CHKERRQ(ierr);
	  //    cout<<"D_elm[col_mat] = "<< D_elm[col_mat1] << endl;
	  col_mat1++;
    
	  for(int ee=0; ee<6; ee++) { //edges
		if(ColGlob[col_mat1].size()!=0) {
		  ierr = GetDataVector(zeroth_field,MBEDGE,D_elm[col_mat1],ee); CHKERRQ(ierr);
		  ierr = GetDataVector(first_field,MBEDGE,D_elm_r[col_mat1],ee); CHKERRQ(ierr);
//          cout<<"Edges D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
		  col_mat1++;
		}
	  }
    
	  for(int ff=0; ff<4; ff++) { //faces
		if(ColGlob[col_mat1].size()!=0) {
		  ierr = GetDataVector(zeroth_field,MBTRI,D_elm[col_mat1],ff); CHKERRQ(ierr);
		  ierr = GetDataVector(first_field,MBTRI,D_elm_r[col_mat1],ff); CHKERRQ(ierr);
//          cout<<"Faces D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
		  col_mat1++;
		}
	  }
    
	  if(ColGlob[col_mat1].size()!=0) { // volumes
		ierr = GetDataVector(zeroth_field,MBTET,D_elm[col_mat1]); CHKERRQ(ierr);
		ierr = GetDataVector(first_field,MBTET,D_elm_r[col_mat1]); CHKERRQ(ierr);
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
	// ======================================================================
    //
    // Postprocessing & Operatating
    //
    // ======================================================================	
	PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }
								 
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
