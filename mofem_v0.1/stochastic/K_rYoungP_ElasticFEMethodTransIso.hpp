/* Copyright (C) 2013, Zahur Ullah <Zahur.Ullah@glasgow.ac.uk>
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

#ifndef __K_RYOUNGP_ELASTICFEMETHODTRANSISO_HPP__
#define __K_RYOUNGP_ELASTICFEMETHODTRANSISO_HPP__


#include <boost/numeric/ublas/symmetric.hpp>
#include "K_rPoissonP_ElasticFEMethodTransIso.hpp"

namespace MoFEM {
  
  
  struct TransverseIsotropicStiffnessMatrix_rYoungP {
    double nu_p, nu_pz, E_p, E_z, G_zp;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rYoungP;
    TransverseIsotropicStiffnessMatrix_rYoungP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
//      cout<<"Hello from TransverseIsotropicStiffnessMatrix_rPoissionP "<<endl;
      StiffnessMatrix_rYoungP.resize(6);
      StiffnessMatrix_rYoungP.clear();
  
      StiffnessMatrix_rYoungP(0,0)=StiffnessMatrix_rYoungP(1,1)=-(E_p*E_p*nu_p-E_p*E_p-2*E_z*E_z*pow(nu_pz,4)+4*E_p*E_z*nu_pz*nu_pz)/((nu_p+1)*pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2));
      StiffnessMatrix_rYoungP(2,2)=(2*E_z*E_z*nu_pz*nu_pz*(nu_p*nu_p-1))/((nu_p+1)*pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2));
      StiffnessMatrix_rYoungP(0,1)=StiffnessMatrix_rYoungP(1,0)=-(E_p*E_p*nu_p*nu_p-E_p*E_p*nu_p+4*E_p*E_z*nu_p*nu_pz*nu_pz+2*E_z*E_z*pow(nu_pz,4))/((nu_p+1)*pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2));
      StiffnessMatrix_rYoungP(0,2)=StiffnessMatrix_rYoungP(2,0)=StiffnessMatrix_rYoungP(1,2)=StiffnessMatrix_rYoungP(2,1)=-(2*E_z*E_z*pow(nu_pz,3))/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);
      StiffnessMatrix_rYoungP(3,3)=1/(2*(nu_p + 1));
      StiffnessMatrix_rYoungP(4,4)=StiffnessMatrix_rYoungP(5,5)=0;
//      cout<<"StiffnessMatrix_rYoungP "<<StiffnessMatrix_rYoungP<<endl<<endl;
//      std::string wait;
//      std::cin >> wait;
    }
    
  };

  struct K_rYoungP_ElasticFEMethodTransIso: public K_rPoissonP_ElasticFEMethodTransIso {
    
    K_rYoungP_ElasticFEMethodTransIso(FieldInterface& _mField,Mat &_Aij,Vec _D,Vec _F):
    K_rPoissonP_ElasticFEMethodTransIso(_mField,_Aij,_D,_F){};

    
    
    virtual PetscErrorCode calculateD_r(double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp) {
      PetscFunctionBegin;
      
      ///Get Stiffness Matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rYoungP;
      StiffnessMatrix_rYoungP.resize(6);
      StiffnessMatrix_rYoungP.clear();
      TransverseIsotropicStiffnessMatrix_rYoungP TranIsoMat(_nu_p,_nu_pz,_E_p,_E_z,_G_zp);
      StiffnessMatrix_rYoungP=TranIsoMat.StiffnessMatrix_rYoungP;
//      cout<<"StiffnessMatrix "<<StiffnessMatrix_rPoissionP <<endl;
      
      ///Rotating the Stiffness matrix according a set of axes of rotations and their respective angle
      D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
      
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
        DummyMatrix = StiffnessMatrix_rYoungP;
        
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
      
      PetscFunctionReturn(0);
    }
    
  };
}
#endif //__K_RYOUNGP_ELASTICFEMETHODTRANSISO_HPP__

