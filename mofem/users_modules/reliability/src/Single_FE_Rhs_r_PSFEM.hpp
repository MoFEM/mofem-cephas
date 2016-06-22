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

#ifndef __SINGLE_FE_RHS_R_PSFEM_HPP__
#define __SINGLE_FE_RHS_R_PSFEM_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
extern "C" {
#include <gm_rule.h>
}

namespace MoFEM {
  
  struct Single_FE_Rhs_r_PSFEM: public TranIsotropicPlyElasticFEMethod{
 
    double young, pois; // young's modulus and poisson's ratio for isotropic material
    Vec dF;
    string ixrandvar; // index for considered random variable
    string material_type;     // Type of material: isotropic or transversly isotropic
    string zeroth_field;
    ublas::matrix<double> Dmat_r;
    ublas::matrix<double> Dmat_r_xyz;
    
    Single_FE_Rhs_r_PSFEM(FieldInterface& _mField,
                          Mat &_Aij,
                          Vec _D,
                          Vec _F,
                          int _noAA,
                          double *_AxVector,
                          double *_AxAngle,
                          string _zeroth_field,
                          string _ixrandvar,
                          string _material_type):
    TranIsotropicPlyElasticFEMethod(_mField,_Aij,_D,_F,_noAA,_AxVector,_AxAngle,_zeroth_field),
    dF(_F),zeroth_field(_zeroth_field),ixrandvar(_ixrandvar), material_type(_material_type){};

// =============================================================================
//
// Calculate material constitutive matrix in global coordinate system for
// transversely isotropic material
//
// =============================================================================
    virtual PetscErrorCode calculateD_r_Trans(double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp) {
      PetscFunctionBegin;

      /*************************************************************************
       *
       * Get the constitutive matrix
       *
       ************************************************************************/
      ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix_r;
      ConstitutiveMatrix_r.resize(6);
      ConstitutiveMatrix_r.clear();
       
      TransverseIsotropicStiffnessMatrix_FirstOrderDerivative mymat;
      if (ixrandvar.compare(0,8,"PoissonP") == 0){
         //cout<<"Poisson's ratio in p-direction"<<endl;
         ierr = mymat.D_r_PoissonP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_r = mymat.StiffnessMatrix_rPoissonP;
      }
      else if (ixrandvar.compare(0,8,"PoissonZ") == 0){
         //cout<<"Poisson's ratio in z-direction"<<endl;
         ierr = mymat.D_r_PoissonPZ(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_r = mymat.StiffnessMatrix_rPoissonPZ;
      }
      else if (ixrandvar.compare(0,6,"YoungP") == 0){
         //cout<<"Young's modulus in p-direction"<<endl;
         ierr = mymat.D_r_YoungP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_r = mymat.StiffnessMatrix_rYoungP;
      }
      else if (ixrandvar.compare(0,6,"YoungZ") == 0){
         //cout<<"Young's modulus in z-direction"<<endl;
         ierr = mymat.D_r_YoungZ(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_r = mymat.StiffnessMatrix_rYoungZ;
         //cout<<ConstitutiveMatrix_r<<"\t";
      }
      else if (ixrandvar.compare(0,7,"ShearZP") == 0){
         //cout<<"Shear modulus in z-direction"<<endl;
         ierr = mymat.D_r_ShearZP(_nu_p,_nu_pz,_E_p,_E_z,_G_zp); CHKERRQ(ierr);
         ConstitutiveMatrix_r = mymat.StiffnessMatrix_rShearZP;
      }
      else {
         cout<<"Invalid input of random variable"<<endl;
      }

      /*************************************************************************
       * 
       * Constitutive matrix transformation
       *
       ************************************************************************/
      
      double AxVector_PSFE[3] = {0.0,1.0,0.0};
      double AxAngle_PSFE[1] = {0.5*M_PI};
      double negAxAngle_PSFE[1]; negAxAngle_PSFE[0] = -AxAngle_PSFE[0];
      
      ublas::matrix<double> DummyMatrix_PSFE, DummyMatrix2_PSFE;
      DummyMatrix_PSFE = ublas::zero_matrix<FieldData>(6,6);
      DummyMatrix_PSFE = ConstitutiveMatrix_r;
      
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
      ConstitutiveMatrix_r.clear(); ConstitutiveMatrix_r = DummyMatrix_PSFE;
      
      
      // Rotating the Stiffness matrix according ply angle
      // It is realized by rotating about z-axis as z-direction is throughthickness
      int noOfRotations = noAA; //Number of Rotations
      double negAxAngle[noOfRotations];
      for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];
      
      ublas::matrix<double> DummyMatrix,DummyMatrix2;
      DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
      DummyMatrix = ConstitutiveMatrix_r;
      Dmat_r = ConstitutiveMatrix_r;
      
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
      
      D.resize(6,6);
      D.clear();
      D = DummyMatrix;
      Dmat_r_xyz = DummyMatrix;
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
      int noOfRotations = noAA; //Number of Rotations
      double negAxAngle[noOfRotations];
      for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];
      
      ublas::matrix<double> DummyMatrix,DummyMatrix2;
      DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
      DummyMatrix = ConstitutiveMatrix_r;
      Dmat_r = ConstitutiveMatrix_r;
      
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
      
      D.resize(6,6);
      D.clear();
      D = DummyMatrix;
      
      PetscFunctionReturn(0);
    }    

// =============================================================================
//
// Calculate element stiffness matrix
//
// =============================================================================
    virtual PetscErrorCode Stiffness() {
      PetscFunctionBegin;
      
      double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
      ierr = GetMatParameters(&_E_p,&_E_z,&_nu_p,&_nu_pz,&_G_zp); CHKERRQ(ierr);
      /**********************************
       *
       * Get material constitutive matrix
       *
       *********************************/
      if (material_type.compare(0,5,"trans") == 0){
         // cout<<"Transversely isotropic material \t";
         ierr = calculateD_r_Trans(_E_p,_E_z,_nu_p,_nu_pz,_G_zp); CHKERRQ(ierr);
//       cout<<"D_At_GaussPoint[0] "<< D_At_GaussPoint[0] <<endl;
      }
      else if (material_type.compare(0,3,"iso") == 0){
        // cout<<"Isotropic material \t";
        double _young,_pois;
        _young = _E_p;
        _pois  = _nu_p;
        ierr = calculateD_r_Iso(_young,_pois); CHKERRQ(ierr);
      }
      else {
         cout<<"Only isotropic or transversely isotropic material is considered currently!"<<endl;
      }

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
                      w,&*D.data().begin(),D.size2(),
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

// =============================================================================
//
// Calculate element force using as right-hand side of finite element 
// equilibrium equation
//
// =============================================================================
    vector<ublas::vector<FieldData> > f_re;
    virtual PetscErrorCode Rhs() {
      PetscFunctionBegin;
      ierr = Stiffness(); CHKERRQ(ierr);
      
      vector<ublas::vector<FieldData> > D_elm;
      
      D_elm.resize(col_mat);
      int col_mat1 = 0;  //nodes
      ierr = GetDataVector(zeroth_field,D_elm[col_mat1]); CHKERRQ(ierr);
      col_mat1++;
      
      
      for(int ee=0; ee<6; ee++) { //edges
        if(ColGlob[col_mat1].size()!=0) {
          ierr = GetDataVector(zeroth_field,MBEDGE,D_elm[col_mat1],ee); CHKERRQ(ierr);
          col_mat1++;
        }
      }
      
      for(int ff=0; ff<4; ff++) { //faces
        if(ColGlob[col_mat1].size()!=0) {
          ierr = GetDataVector(zeroth_field,MBTRI,D_elm[col_mat1],ff); CHKERRQ(ierr);
          col_mat1++;
        }
      }
      
      if(ColGlob[col_mat1].size()!=0) {
        ierr = GetDataVector(zeroth_field,MBTET,D_elm[col_mat1]); CHKERRQ(ierr);
      }
      
      f_re.resize(row_mat);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        
        int rr_start=0;
        for(int cc = 0;cc<col_mat;cc++) {
          if(ColGlob[cc].size()==0) continue;
          if(rr_start == 0) {
            f_re[rr] =  -prod(K(rr,cc),D_elm[cc]);
            rr_start++;
          } else {
            f_re[rr] -= prod(K(rr,cc),D_elm[cc]);
          }
        }
      }
      
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        if(RowGlob[rr].size()!=f_re[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        ierr = VecSetValues(dF,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_re[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
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

#endif //__TRANS_ISO_RHS_R_PSFEM_HPP__
