/* Copyright (C) 2014, 
 *   Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
 *   Xiao-Yi Zhou (xiaoyi.zhou AT newcastle.ac.uk)
 * --------------------------------------------------------------
 * This routine calculates the first-order partial derivative of right-hand 
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

#ifndef __TRANS_ISO_RHS_R_PSFEM_DEGRADATION_HPP__
#define __TRANS_ISO_RHS_R_PSFEM_DEGRADATION_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
extern "C" {
#include <gm_rule.h>
}

namespace MoFEM {
  
  struct Trans_Iso_Rhs_r_PSFEM_Degradation: public Trans_Iso_Rhs_r_PSFEM {
    
    double wt_Gauss;

    Trans_Iso_Rhs_r_PSFEM_Degradation(FieldInterface& _mField,
                          Mat _Aij,
                          Vec _D,
                          Vec _F,
                          string _zeroth_field,
                          string _ixrandvar,
                          string _material_type,
                          string _material_function,
                          double _wt_Gauss):
    Trans_Iso_Rhs_r_PSFEM(_mField,_Aij,_D,_F,_zeroth_field, _ixrandvar,_material_type,_material_function),
    wt_Gauss(_wt_Gauss){};

// =============================================================================
//
// Calculate material constitutive matrix in global coordinate system for
// isotropic material
//
// ============================================================================= 
    virtual PetscErrorCode calculateD_r_Iso(double _young, double _nu) {
      PetscFunctionBegin;
      /*************************************************************************
       *
       * Get the constitutive matrix
       *
       ************************************************************************/
      ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix_r;
      ConstitutiveMatrix_r.resize(6);
      ConstitutiveMatrix_r.clear();
       
      IsotropicStiffnessMatrix_FirstOrderDerivative mymat;
      if (ixrandvar.compare(0,5,"Young") == 0 ){
         //cout<<"Young's modulus of isotropic material"<<endl;
         ierr = mymat.D_r_Young(_young,_nu); CHKERRQ(ierr);
         ConstitutiveMatrix_r = mymat.StiffnessMatrix_rYoung;
      }
      else if (ixrandvar.compare(0,7,"Poisson") == 0){
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
            D_At_GaussPoint[gg] = wt_Gauss*ConstitutiveMatrix_r;
         }
      }
      else {
         cout<<"Undefined material function!"<<endl;
      }
      PetscFunctionReturn(0);
    }
    
  };
  
}

#endif //__TRANS_ISO_RHS_R_PSFEM_DEGRADATION_HPP__
