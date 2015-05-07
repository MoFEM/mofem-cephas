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

#ifndef __TRANS_ISO_RHS_R_PSFEM_HPP__
#define __TRANS_ISO_RHS_R_PSFEM_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
extern "C" {
#include <gm_rule.h>
}

namespace MoFEM {
  
  struct Trans_Iso_Rhs_r_PSFEM: public TranIsotropicFibreDirRotElasticFEMethod {
 
    double young, pois; // young's modulus and poisson's ratio for isotropic material
    Vec dF;
    const string ixrandvar; // index for considered random variable 
    const string material_type;     // Type of material: isotropic or transversly isotropic
    const string material_function; // Function of material: matrix or reinforcement/inclusion/fibre

    Trans_Iso_Rhs_r_PSFEM(FieldInterface& _mField,Mat &_Aij,Vec _D,Vec _F, const string& _ixrandvar, const string& _material_type, const string& _material_function):
    TranIsotropicFibreDirRotElasticFEMethod(_mField,_Aij,_D,_F),dF(_F), ixrandvar(_ixrandvar), material_type(_material_type), material_function(_material_function){};

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
            D_At_GaussPoint[gg] = ConstitutiveMatrix_r;
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
    virtual PetscErrorCode Stiffness() {
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
//       cout<<"_E_p "<< _E_p <<endl;
//       cout<<"_E_z "<< _E_z <<endl;
//       cout<<"_nu_p "<< _nu_p <<endl;
//       cout<<"_nu_pz "<< _nu_pz <<endl;
//       cout<<"_G_zp "<< _G_zp <<endl;
         ierr = calculateD_r_Trans(_E_p,_E_z,_nu_p,_nu_pz,_G_zp); CHKERRQ(ierr);
//       cout<<"D_At_GaussPoint[0] "<< D_At_GaussPoint[0] <<endl;
      }
      else if (material_type.compare(0,3,"iso") == 0){
         // cout<<"Isotropic material \t";
         double _young,_pois;
  	 for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)){
            //cout << endl << *it << endl;
    	    //Get block name
            string name = it->get_name();
            if (name.compare(0,13,"MAT_ELASTIC_1") == 0){
               Mat_Elastic mydata;
               ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
               //cout << mydata;
              _young = mydata.data.Young;
              _pois  = mydata.data.Poisson;
            }
         }

         //ierr = GetMatParameters_Iso(&_young,&_pois); CHKERRQ(ierr);
         ierr = calculateD_r_Iso(_young,_pois); CHKERRQ(ierr);
         //cout<<"_young "<< _young <<endl;
         //cout<<"_pois "<< _pois <<endl;
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

// =============================================================================
//
// Calculate element force using as right-hand side of finite element 
// equilibrium equation
//
// =============================================================================
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
  
  /*****************************************************************************
   *
   * First-order derivative for right-hand-side with respect to:
   *  1. fibre waviness:
   *      A - waviness amplitude
   *      l - waviness periodic length
   *  2. fibre misalignment:
   *      theta - misalignment angle
   *  3. fibre volume fraction:
   *      vf - fibre volume fraction
   *
   ****************************************************************************/
  
  
  struct Trans_Iso_Geom_Rhs_r_PSFEM: public TranIso_FibreWavinessElasticFEMethod {
		Vec dF;
		const string ixrandvar; // index for considered random variable 
    const string material_type;     // Type of material: isotropic or transversly isotropic
    const string material_function; // Function of material: matrix or reinforcement/inclusion/fibre
	
    Trans_Iso_Geom_Rhs_r_PSFEM(FieldInterface& _mField, Mat &_Aij, Vec _D,
                                      Vec _F, const string& _ixrandvar,
                                      const string& _material_type,
                                      const string& _material_function ):
    TranIso_FibreWavinessElasticFEMethod(_mField,_Aij,_D,_F),dF(_F), 
										     ixrandvar(_ixrandvar), 
											 material_type(_material_type), 
											 material_function(_material_function){}
    
    // =========================================================================
    //
    // Calculate the first-order derivative of material constitutive matrix for
    //   yarn
    //
    // =========================================================================
    virtual PetscErrorCode calculateD_r_Waviness(double _E_p, double _E_z, double _nu_p,
                                             double _nu_pz, double _G_zp,
                                             double _lambda, double _mu,
                                             double _vf,
                                             double _theta_f, double _WavinessFactor) {
      PetscFunctionBegin;
      
      ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix_r;
      ConstitutiveMatrix_r.resize(6);
      ConstitutiveMatrix_r.clear();
	  
      ublas::symmetric_matrix<FieldData,ublas::upper> ConstitutiveMatrix;
      ConstitutiveMatrix.resize(6);
      ConstitutiveMatrix.clear();
	  
	  // -----------
	  // 1. Get the original compliance matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> SMat;
      SMat.resize(6);
      SMat.clear();
      
      //TransverseIsotropicComplianceMatrix TranIsoMat_S(_nu_p,_nu_pz,_E_p,_E_z,_G_zp);
      YarnComplianceMatrix TranIsoMat_S(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,
                                       _lambda, _mu, _vf);
      SMat = TranIsoMat_S.ComplianceMatrix;
	  // cout<<"Original "<<SMat<<endl;
	  
      // 2. Get the wavniess parameters
      double WavinessAmplitude, WavinessLength;
      double dI1, dI3, dI5, dI6, dI8; // 1st-order derivative
      //WavinessAmplitude = 1.19; // unit: mm
      WavinessLength = 27.9;   // unit: mm
      WavinessAmplitude = WavinessLength*_WavinessFactor;
	  
      ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix_r;
      ComplianceMatrix_r.resize(6);
      ComplianceMatrix_r.clear();
	  
      ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix;
      ComplianceMatrix.resize(6);
      ComplianceMatrix.clear();
      
      YarnStiffnessMatrix_Geom_FirstOrderDerivative mymat;
      if (ixrandvar.compare(0,9,"Amplitude") == 0) {
        ierr = mymat.D_r_Amplitude(WavinessAmplitude,WavinessLength); CHKERRQ(ierr);
			dI1 = mymat.dI1;
			dI3 = mymat.dI3;
			dI5 = mymat.dI5;
			dI6 = mymat.dI6;
			dI8 = mymat.dI8;
		    //cout<<"dI1 = "<<dI1<<"\t dI3 = "<<dI3<<"\t dI5 = "<<dI5<<"\t dI6 = "<<dI6<<"\t dI8 = "<<dI8<<endl;
		}
		else if (ixrandvar.compare(0,6,"Length") == 0) {
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
		
	// =========================================================================
	//
	// Calculate the first-order derivative of material constitutive matrix for
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
    
	// =========================================================================
    //
    // Calculate the first-order derivative of material constitutive matrix
	//   for yarn with respect to fibre volume fraction
    //
    // =========================================================================
	virtual PetscErrorCode calculateD_r_Fraction(double _E_p, double _E_z, double _nu_p,
                                             double _nu_pz, double _G_zp,
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
	
	// =========================================================================
    //
    // Calculate the first-order derivative of material constitutive matrix
	//   for yarn with respect to fibre volume fraction
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
	  double dk_c, dm_c, dp_c, dl_c, dn_c; // 1st-order fibre volume fraction
	  YarnStiffnessMatrix_Geom_FirstOrderDerivative mymat;
	  double Em, NUm;
	  Em = _mu*(3*_lambda+2*_mu)/(_lambda+_mu);
      NUm = _lambda/2/(_lambda+_mu); 
      if (ixrandvar.compare(0,6,"YoungZ") == 0) {
        ierr = mymat.D_r_YoungZ(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
		}
	  else if (ixrandvar.compare(0,6,"YoungP") == 0) {
		ierr = mymat.D_r_YoungP(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ixrandvar.compare(0,8,"PoissonP") == 0) {
		ierr = mymat.D_r_PoissonP(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ixrandvar.compare(0,8,"PoissonZ") == 0) {
		ierr = mymat.D_r_PoissonPZ(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ixrandvar.compare(0,6,"YoungM") == 0) {
		ierr = mymat.D_r_Young(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ixrandvar.compare(0,8,"PoissonM") == 0) {
		ierr = mymat.D_r_Poisson(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,Em, NUm, _vf); CHKERRQ(ierr);
	  }
	  else if (ixrandvar.compare(0,7,"ShearZP") == 0) {
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
		
    // =========================================================================
    //
    // Calculate element stiffness matrix
    //
    // =========================================================================
    virtual PetscErrorCode Stiffness() {
      
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
      if ((ixrandvar.compare(0,9,"Amplitude") == 0) || (ixrandvar.compare(0,6,"Length") == 0)) {
        ierr = calculateD_r_Waviness(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                                 _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);
      }
      else if (ixrandvar.compare(0,5,"Angle") == 0) {
        ierr = calculateD_r_Misalignment(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                                         _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);
      }
      else if (ixrandvar.compare(0,8,"Fraction") == 0) {
        ierr = calculateD_r_Fraction(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                                     _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);
      }
	  else if ((ixrandvar.compare(0,6,"YoungZ") == 0) ||
			   (ixrandvar.compare(0,6,"YoungP") == 0) ||
			   (ixrandvar.compare(0,8,"PoissonP") == 0) ||
			   (ixrandvar.compare(0,8,"PoissonZ") == 0) ||
			   (ixrandvar.compare(0,7,"ShearZP") == 0) ||
			   (ixrandvar.compare(0,6,"YoungM") == 0) ||
			   (ixrandvar.compare(0,8,"PoissonM") == 0)) {
        ierr = calculateD_r_Material(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                                         _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);
      }
      else {
        //cout<<"Invalid input of random variable"<<endl;
      }
      
      
      // -----------
      // 3. Calculate element stiffness matrix
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
	  //cout<<"K matrix "<<K<<endl;
      
      PetscFunctionReturn(0);
    }
    
    // =========================================================================
    //
    // Calculate element force to establish RHS of FE equilibrium equation
    //
    // =========================================================================
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
    // =========================================================================
    //
    // Postprocessing & Operatating
    //
    // =========================================================================
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      
      ierr = Get_g_NTET(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      //cout<<"Hi from K_rPoissonP_ElasticFEMethodTransIso "<<endl;
      ierr = GetMatrices(); CHKERRQ(ierr);
      ierr = Rhs(); CHKERRQ(ierr);
	  //cout<<"Finish right hand side"<<endl;
      //      std::string wait;
      //      std::cin >> wait;
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  };
  
}

#endif //__TRANS_ISO_RHS_R_PSFEM_HPP__
