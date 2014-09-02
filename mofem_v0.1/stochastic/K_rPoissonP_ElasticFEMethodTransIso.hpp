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

#ifndef __K_RPOISSONP_ELASTICFEMETHODTRANSISO_HPP__
#define __K_RPOISSONP_ELASTICFEMETHODTRANSISO_HPP__


#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFEMethodTransIso.hpp"

namespace MoFEM {
  
  
  struct TransverseIsotropicStiffnessMatrix_rPoissionP {
    double nu_p, nu_pz, E_p, E_z, G_zp;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rPoissionP;
    TransverseIsotropicStiffnessMatrix_rPoissionP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
//      cout<<"Hello from TransverseIsotropicStiffnessMatrix_rPoissionP "<<endl;
      StiffnessMatrix_rPoissionP.resize(6);
      StiffnessMatrix_rPoissionP.clear();
  
      StiffnessMatrix_rPoissionP(0,0)=StiffnessMatrix_rPoissionP(1,1)=(2*E_p*(-E_z*nu_pz*nu_pz+E_p)*(E_z*nu_pz*nu_pz+E_p*nu_p))/(pow((nu_p+1),2)*pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2));
      StiffnessMatrix_rPoissionP(2,2)=(2*E_p*E_z*E_z*nu_pz*nu_pz)/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);
      StiffnessMatrix_rPoissionP(0,1)=StiffnessMatrix_rPoissionP(1,0)=(E_p*(E_p*E_p*nu_p*nu_p+E_p*E_p+2*E_p*E_z*nu_p*nu_pz*nu_pz-2*E_p*E_z*nu_pz*nu_pz+2*E_z*E_z*pow(nu_pz,4)))/(pow((nu_p+1),2)*pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2));
      StiffnessMatrix_rPoissionP(0,2)=StiffnessMatrix_rPoissionP(2,0)=StiffnessMatrix_rPoissionP(1,2)=StiffnessMatrix_rPoissionP(2,1)=(E_p*E_p*E_z*nu_pz)/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),2);
      StiffnessMatrix_rPoissionP(3,3)=-E_p/(2*pow((nu_p + 1),2));
      StiffnessMatrix_rPoissionP(4,4)=StiffnessMatrix_rPoissionP(5,5)=0;
    }
    
  };

  struct K_rPoissonP_ElasticFEMethodTransIso: public TranIsotropicFibreDirRotElasticFEMethod {
    
    Vec dF;
    K_rPoissonP_ElasticFEMethodTransIso(FieldInterface& _mField,Mat &_Aij,Vec _D,Vec _F):
    TranIsotropicFibreDirRotElasticFEMethod(_mField,_Aij,_D,_F),dF(_F){};

    
    
    virtual PetscErrorCode calculateD_r(double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp) {
      PetscFunctionBegin;
      
      ///Get Stiffness Matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rPoissionP;
      StiffnessMatrix_rPoissionP.resize(6);
      StiffnessMatrix_rPoissionP.clear();
      TransverseIsotropicStiffnessMatrix_rPoissionP TranIsoMat(_nu_p,_nu_pz,_E_p,_E_z,_G_zp);
      StiffnessMatrix_rPoissionP=TranIsoMat.StiffnessMatrix_rPoissionP;
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
        DummyMatrix = StiffnessMatrix_rPoissionP;
        
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

    
    
    
    PetscErrorCode Stiffness() {
      PetscFunctionBegin;
      
      double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
      ierr = GetMatParameters(&_E_p,&_E_z,&_nu_p,&_nu_pz,&_G_zp); CHKERRQ(ierr);
//      cout<<"_E_p "<< _E_p <<endl;
//      cout<<"_E_z "<< _E_z <<endl;
//      cout<<"_nu_p "<< _nu_p <<endl;
//      cout<<"_nu_pz "<< _nu_pz <<endl;
//      cout<<"_G_zp "<< _G_zp <<endl;
      ierr = calculateD_r(_E_p,_E_z,_nu_p,_nu_pz,_G_zp); CHKERRQ(ierr);
//      cout<<"D_At_GaussPoint[0] "<< D_At_GaussPoint[0] <<endl;
      K.resize(row_mat,col_mat);
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
          //ublas::noalias(BD) = prod( w*D,row_Mat );
          cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
                      BD.size1(),BD.size2(),
                      w,&*D_At_GaussPoint[gg].data().begin(),D_At_GaussPoint[gg].size2(),
                      &*row_Mat.data().begin(),row_Mat.size2(),
                      0.,&*BD.data().begin(),BD.size2());
          for(int cc = 0;cc<col_mat;cc++) {
            if(ColGlob[cc].size()==0) continue;
            ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
            if(gg == 0) {
              K(rr,cc).resize(BD.size2(),col_Mat.size2());
              //ublas::noalias(K(rr,cc)) = prod(trans(BD) , col_Mat ); // int BT*D*B
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          0.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
            } else {
              //ublas::noalias(K(rr,cc)) += prod(trans(BTD) , col_Mat ); // int BT*D*B
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          1.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
            }
          }
        }
      }
      PetscFunctionReturn(0);
      
    }

    vector<ublas::vector<FieldData> > f_re;
    virtual PetscErrorCode Rhs() {
      PetscFunctionBegin;
//      cout<<"Before stiffness matix "<<endl;
      ierr = Stiffness(); CHKERRQ(ierr);
      
      vector<ublas::vector<FieldData> > D_elm;
//       cout<<"col_mat = "<< col_mat << endl;
      
      D_elm.resize(col_mat);
      int col_mat1 = 0;  //nodes
      ierr = GetDataVector("DISPLACEMENT",D_elm[col_mat1]); CHKERRQ(ierr);
//       cout<<"D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
      col_mat1++;
      
      
      for(int ee=0; ee<6; ee++) { //edges
        if(ColGlob[col_mat1].size()!=0) {
          ierr = GetDataVector("DISPLACEMENT",MBEDGE,D_elm[col_mat1],ee); CHKERRQ(ierr);
//          cout<<"Edges D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
          col_mat1++;
        }
      }
      
      for(int ff=0; ff<4; ff++) { //faces
        if(ColGlob[col_mat1].size()!=0) {
          ierr = GetDataVector("DISPLACEMENT",MBTRI,D_elm[col_mat1],ff); CHKERRQ(ierr);
//          cout<<"Faces D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
          col_mat1++;
        }
      }
      
      if(ColGlob[col_mat1].size()!=0) {
        ierr = GetDataVector("DISPLACEMENT",MBTET,D_elm[col_mat1]); CHKERRQ(ierr);
//        cout<<"Faces D_elm[col_mat] = "<< D_elm[col_mat1] << endl;
      }
      
      f_re.resize(row_mat);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        
        int rr_start=0;
        for(int cc = 0;cc<col_mat;cc++) {
          if(ColGlob[cc].size()==0) continue;
//          cout<<"rr "<<rr<<endl;
//          cout<<"cc "<<cc<<endl;
          if(rr_start == 0) {
//            cout<<"K(rr,cc) "<<K(rr,cc)<<endl;
//            cout<<"D_elm[cc] "<<D_elm[cc]<<endl;
            f_re[rr] =  -prod(K(rr,cc),D_elm[cc]);
            rr_start++;
          } else {
            f_re[rr] -= prod(K(rr,cc),D_elm[cc]);
          }
        }
//        cout<<"f_re[rr] "<<f_re[rr]<<endl;
      }
      
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        if(RowGlob[rr].size()!=f_re[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        ierr = VecSetValues(dF,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_re[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      
      PetscFunctionReturn(0);
    }



    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      
      ierr = Get_g_NTET(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
//      cout<<"Hi from K_rPoissonP_ElasticFEMethodTransIso "<<endl;
      ierr = GetMatrices(); CHKERRQ(ierr);
      ierr = Rhs(); CHKERRQ(ierr);
//      std::string wait;
//      std::cin >> wait;
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
   
    
  };
}
#endif //__K_RPOISSONP_ELASTICFEMETHODTRANSISO_HPP__

