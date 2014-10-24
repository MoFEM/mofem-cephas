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

#ifndef __K_RSPOISSONP_ELASTICFEMETHODTRANSISO_HPP__
#define __K_RSPOISSONP_ELASTICFEMETHODTRANSISO_HPP__


#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFEMethodTransIso.hpp"

namespace MoFEM {
  
  
  struct TransverseIsotropicStiffnessMatrix_rsPoissionP {
    double nu_p, nu_pz, E_p, E_z, G_zp;
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsPoissionP;
    TransverseIsotropicStiffnessMatrix_rsPoissionP(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
//      cout<<"Hello from TransverseIsotropicStiffnessMatrix_rsPoissionP "<<endl;
      StiffnessMatrix_rsPoissionP.resize(6);
      StiffnessMatrix_rsPoissionP.clear();
      
      StiffnessMatrix_rsPoissionP(0,0)=StiffnessMatrix_rsPoissionP(1,1)=-(2*E_p*(-E_z*nu_pz*nu_pz+E_p)*(3*E_p*E_p*nu_p*nu_p+E_p*E_p+6*E_p*E_z*nu_p*nu_pz*nu_pz-2*E_p*E_z*nu_pz*nu_pz+4*E_z*E_z*pow(nu_pz,4)))/(pow((nu_p+1),3)*pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3));
      StiffnessMatrix_rsPoissionP(2,2)=-(4*E_p*E_p*E_z*E_z*nu_pz*nu_pz)/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);
      StiffnessMatrix_rsPoissionP(0,1)=StiffnessMatrix_rsPoissionP(1,0)=-(2*E_p*(pow(E_p,3)*pow(nu_p,3)+3*pow(E_p,3)*nu_p+3*E_p*E_p*E_z*nu_p*nu_p*nu_pz*nu_pz-6*E_p*E_p*E_z*nu_p*nu_pz*nu_pz+3*E_p*E_p*E_z*nu_pz*nu_pz+6*E_p*E_z*E_z*nu_p*pow(nu_pz,4)-6*E_p*E_z*E_z*pow(nu_pz,4)+4*pow(E_z,3)*pow(nu_pz,6)))/(pow((nu_p+1),3)*pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3));
      StiffnessMatrix_rsPoissionP(0,2)=StiffnessMatrix_rsPoissionP(2,0)=StiffnessMatrix_rsPoissionP(1,2)=StiffnessMatrix_rsPoissionP(2,1)=-(2*pow(E_p,3)*E_z*nu_pz)/pow((2*E_z*nu_pz*nu_pz-E_p+E_p*nu_p),3);
      StiffnessMatrix_rsPoissionP(3,3)=E_p/pow((nu_p + 1),3);
      StiffnessMatrix_rsPoissionP(4,4)=StiffnessMatrix_rsPoissionP(5,5)=0;
    }
  };

  

  struct K_rsPoissonP_ElasticFEMethodTransIso: public TranIsotropicFibreDirRotElasticFEMethod {
    
    Vec ddF;
    const string second_field;
    
    K_rsPoissonP_ElasticFEMethodTransIso(FieldInterface& _mField,Mat &_Aij,Vec _D,Vec _F,const string& _second_field):
    TranIsotropicFibreDirRotElasticFEMethod(_mField,_Aij,_D,_F),ddF(_F),second_field(_second_field){};

    
    
    //*********************************************************************************************************************
    virtual PetscErrorCode calculateD_r(double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp) {
      PetscFunctionBegin;
      ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rPoissionP;
      StiffnessMatrix_rPoissionP.resize(6);
      StiffnessMatrix_rPoissionP.clear();
      TransverseIsotropicStiffnessMatrix_rPoissionP TranIsoMat(_nu_p,_nu_pz,_E_p,_E_z,_G_zp);
      StiffnessMatrix_rPoissionP=TranIsoMat.StiffnessMatrix_rPoissionP;
//      cout<<"StiffnessMatrix_rPoissionP "<<StiffnessMatrix_rPoissionP <<endl;
      
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
    //*********************************************************************************************************************

    virtual PetscErrorCode calculateD_rs(double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp) {
      PetscFunctionBegin;
      
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix_rsPoissionP;
    StiffnessMatrix_rsPoissionP.resize(6);
    StiffnessMatrix_rsPoissionP.clear();
    TransverseIsotropicStiffnessMatrix_rsPoissionP TranIsoMat(_nu_p,_nu_pz,_E_p,_E_z,_G_zp);
    StiffnessMatrix_rsPoissionP=TranIsoMat.StiffnessMatrix_rsPoissionP;
//    cout<<"StiffnessMatrix_rsPoissionP "<<StiffnessMatrix_rsPoissionP <<endl;
    
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
      DummyMatrix = StiffnessMatrix_rsPoissionP;
      
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
    //*********************************************************************************************************************
    ublas::matrix<ublas::matrix<FieldData> > K_r;
    PetscErrorCode StiffnessK_r() {
      PetscFunctionBegin;
      
      double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
      ierr = GetMatParameters(&_E_p,&_E_z,&_nu_p,&_nu_pz,&_G_zp); CHKERRQ(ierr);
      ierr = calculateD_r(_E_p,_E_z,_nu_p,_nu_pz,_G_zp); CHKERRQ(ierr);
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
    //*********************************************************************************************************************
    
    ublas::matrix<ublas::matrix<FieldData> > K_rs;
    PetscErrorCode StiffnessK_rs() {
      PetscFunctionBegin;
      
      double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
      ierr = GetMatParameters(&_E_p,&_E_z,&_nu_p,&_nu_pz,&_G_zp); CHKERRQ(ierr);
      ierr = calculateD_rs(_E_p,_E_z,_nu_p,_nu_pz,_G_zp); CHKERRQ(ierr);
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

    //*********************************************************************************************************************

    
  //F: calculate the second-order partial derivative of "external force" which
  //   is referred as right-hand side in the algebriac equation [K][U_rs] = [F_rs]
  //   [F_rs] = - [K_rs][U] - 2[K_r][U_s]
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
#endif //__K_RSPOISSONP_ELASTICFEMETHODTRANSISO_HPP__

