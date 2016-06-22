/* Copyright (C) 2015, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
 * --------------------------------------------------------------
 *
 * Description: Implementation of thermal stress, i.e. right hand side as result of thermal stresses
 *
 * This is not exactly procedure for linear elatic dynamics, since jacobian is
 * evaluated at every time step and snes procedure is involved. However it is
 * implemented like that, to test methodology for general nonlinear problem.
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

#ifndef __FE2_MACRO_SOLVER_HPP
#define __FE2_MACRO_SOLVER_HPP

namespace MoFEM {
  
  /*****************************************************************************
   *                                                                           *
   *                                                                           *
   *                         MULTIPLE LAYERS LAMINATE                          *
   *                                                                           *
   *                                                                           *
   ****************************************************************************/
  
  struct FE2_Macro_Solver_Laminate {
    
    //global variable Dmat
    ublas::matrix<double> Dmat;
    ublas::matrix<double> Dmat_r_Em;
    ublas::matrix<double> Dmat_r_NUm;
    ublas::matrix<double> Dmat_r_Ep;
    ublas::matrix<double> Dmat_r_Ez;
    ublas::matrix<double> Dmat_r_NUp;
    ublas::matrix<double> Dmat_r_NUpz;
    ublas::matrix<double> Dmat_r_Gzp;
    ublas::matrix<double> Dmat_r_Ef;
    ublas::matrix<double> Dmat_r_NUf;
    
    ublas::matrix<double> Dmat_1st_Ply;
    ublas::matrix<double> Dmat_1st_Ply_r_Em;
    ublas::matrix<double> Dmat_1st_Ply_r_NUm;
    ublas::matrix<double> Dmat_1st_Ply_r_Ep;
    ublas::matrix<double> Dmat_1st_Ply_r_Ez;
    ublas::matrix<double> Dmat_1st_Ply_r_NUp;
    ublas::matrix<double> Dmat_1st_Ply_r_NUpz;
    ublas::matrix<double> Dmat_1st_Ply_r_Gzp;
    ublas::matrix<double> Dmat_1st_Ply_r_Ef;
    ublas::matrix<double> Dmat_1st_Ply_r_NUf;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta_1;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta_2;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta_3;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta_4;
    
    ublas::matrix<double> Dmat_2nd_Ply;
    ublas::matrix<double> Dmat_2nd_Ply_r_Em;
    ublas::matrix<double> Dmat_2nd_Ply_r_NUm;
    ublas::matrix<double> Dmat_2nd_Ply_r_Ep;
    ublas::matrix<double> Dmat_2nd_Ply_r_Ez;
    ublas::matrix<double> Dmat_2nd_Ply_r_NUp;
    ublas::matrix<double> Dmat_2nd_Ply_r_NUpz;
    ublas::matrix<double> Dmat_2nd_Ply_r_Gzp;
    ublas::matrix<double> Dmat_2nd_Ply_r_Ef;
    ublas::matrix<double> Dmat_2nd_Ply_r_NUf;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta_1;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta_2;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta_3;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta_4;
    
    ublas::matrix<double> Dmat_3rd_Ply;
    ublas::matrix<double> Dmat_3rd_Ply_r_Em;
    ublas::matrix<double> Dmat_3rd_Ply_r_NUm;
    ublas::matrix<double> Dmat_3rd_Ply_r_Ep;
    ublas::matrix<double> Dmat_3rd_Ply_r_Ez;
    ublas::matrix<double> Dmat_3rd_Ply_r_NUp;
    ublas::matrix<double> Dmat_3rd_Ply_r_NUpz;
    ublas::matrix<double> Dmat_3rd_Ply_r_Gzp;
    ublas::matrix<double> Dmat_3rd_Ply_r_Ef;
    ublas::matrix<double> Dmat_3rd_Ply_r_NUf;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta_1;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta_2;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta_3;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta_4;
    
    ublas::matrix<double> Dmat_4th_Ply;
    ublas::matrix<double> Dmat_4th_Ply_r_Em;
    ublas::matrix<double> Dmat_4th_Ply_r_NUm;
    ublas::matrix<double> Dmat_4th_Ply_r_Ep;
    ublas::matrix<double> Dmat_4th_Ply_r_Ez;
    ublas::matrix<double> Dmat_4th_Ply_r_NUp;
    ublas::matrix<double> Dmat_4th_Ply_r_NUpz;
    ublas::matrix<double> Dmat_4th_Ply_r_Gzp;
    ublas::matrix<double> Dmat_4th_Ply_r_Ef;
    ublas::matrix<double> Dmat_4th_Ply_r_NUf;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta_1;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta_2;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta_3;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta_4;
    
    //------------------------------------------------------------------------------
    // To transform constitutive matrix
    //
    
    PetscErrorCode Dmat_Transformation_old(double theta,
                                       ublas::matrix<FieldData> Dmat_123,
                                       ublas::matrix<FieldData> &Dmat_xyz) {
      //virtual PetscErrorCode FibreDirection_z2x(ublas::matrix<double> &StiffnessMatrix) {
      PetscFunctionBegin;
      
      double AxVector_PSFE[3] = {0.0,0.0,1.0};
      double AxAngle_PSFE[1] = {theta};
      double negAxAngle_PSFE[1]; negAxAngle_PSFE[0] = -AxAngle_PSFE[0];
      
      ublas::matrix<double> DummyMatrix_PSFE, DummyMatrix2_PSFE;
      DummyMatrix_PSFE = ublas::zero_matrix<FieldData>(6,6);
      DummyMatrix_PSFE = Dmat_123;
      
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
      Dmat_xyz.clear(); Dmat_xyz = DummyMatrix_PSFE;
      
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode Dmat_Transformation(double theta, ublas::matrix<FieldData> Dmat_123, ublas::matrix<FieldData> &Dmat_xyz) {
      PetscFunctionBegin;
      
      double l1, l2, l3, m1, m2, m3, n1, n2, n3;
      l1=cos(theta);
      m1=sin(theta);
      n1=0.0;
      l2=-sin(theta);
      m2=cos(theta);
      n2=0.0;
      l3=0.0;
      m3=0.0;
      n3=1.0;
      
      ublas::matrix<FieldData> T_strain;    T_strain.resize(6,6); T_strain.clear();
      T_strain(0,0)=l1*l1;
      T_strain(0,1)=m1*m1;
      T_strain(0,2)=n1*n1;
      T_strain(0,3)=l1*m1;
      T_strain(0,4)=m1*n1;
      T_strain(0,5)=l1*n1;
      
      T_strain(1,0)=l2*l2;
      T_strain(1,1)=m2*m2;
      T_strain(1,2)=n2*n2;
      T_strain(1,3)=l2*m2;
      T_strain(1,4)=m2*n2;
      T_strain(1,5)=l2*n2;
      
      T_strain(2,0)=l3*l3;
      T_strain(2,1)=m3*m3;
      T_strain(2,2)=n3*n3;
      T_strain(2,3)=l3*m3;
      T_strain(2,4)=m3*n3;
      T_strain(2,5)=l3*n3;
      
      T_strain(3,0)=2*l1*l2;
      T_strain(3,1)=2*m1*m2;
      T_strain(3,2)=2*n1*n2;
      T_strain(3,3)=l1*m2+m1*l2;
      T_strain(3,4)=m1*n2+n1*m2;
      T_strain(3,5)=l1*n2+n1*l2;
      
      T_strain(4,0)=2*l2*l3;
      T_strain(4,1)=2*m2*m3;
      T_strain(4,2)=2*n2*n3;
      T_strain(4,3)=l2*m3+m2*l3;
      T_strain(4,4)=m2*n3+n2*m3;
      T_strain(4,5)=l2*n3+n2*l3;
      
      T_strain(5,0)=2*l1*l3;
      T_strain(5,1)=2*m1*m3;
      T_strain(5,2)=2*n1*n3;
      T_strain(5,3)=l1*m3+m1*l3;
      T_strain(5,4)=m1*n3+n1*m3;
      T_strain(5,5)=l1*n3+n1*l3;
      
      //cout<<"\n\nT_strain = "<<T_strain<<endl;
      
      ublas::matrix<FieldData> Mat1=prod(Dmat_123,T_strain);
      Dmat_xyz = prod(trans(T_strain), Mat1);
      
      
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode Dmat_Transformation_r_Theta(double theta, ublas::matrix<FieldData> Dmat_123, ublas::matrix<FieldData> &Dmat_xyz) {
      PetscFunctionBegin;
      
      double l1, l2, l3, m1, m2, m3, n1, n2, n3;
      double l1_r, l2_r, l3_r, m1_r, m2_r, m3_r, n1_r, n2_r, n3_r;
      l1   =  cos(theta); l1_r = -sin(theta);
      m1   =  sin(theta); m1_r =  cos(theta);
      n1   =  0.0;        n1_r =  0.0;
      l2   = -sin(theta); l2_r = -cos(theta);
      m2   =  cos(theta); m2_r = -sin(theta);
      n2   =  0.0;        n2_r =  0.0;
      l3   =  0.0;        l3_r =  0.0;
      m3   =  0.0;        m3_r =  0.0;
      n3   =  1.0;        n3_r =  0.0;
      
      ublas::matrix<FieldData> T_strain;    T_strain.resize(6,6); T_strain.clear();
      T_strain(0,0) = l1*l1;
      T_strain(0,1) = m1*m1;
      T_strain(0,2) = n1*n1;
      T_strain(0,3) = l1*m1;
      T_strain(0,4) = m1*n1;
      T_strain(0,5) = l1*n1;
      
      T_strain(1,0) = l2*l2;
      T_strain(1,1) = m2*m2;
      T_strain(1,2) = n2*n2;
      T_strain(1,3) = l2*m2;
      T_strain(1,4) = m2*n2;
      T_strain(1,5) = l2*n2;
      
      T_strain(2,0) = l3*l3;
      T_strain(2,1) = m3*m3;
      T_strain(2,2) = n3*n3;
      T_strain(2,3) = l3*m3;
      T_strain(2,4) = m3*n3;
      T_strain(2,5) = l3*n3;
      
      T_strain(3,0) = 2*l1*l2;
      T_strain(3,1) = 2*m1*m2;
      T_strain(3,2) = 2*n1*n2;
      T_strain(3,3) = l1*m2+m1*l2;
      T_strain(3,4) = m1*n2+n1*m2;
      T_strain(3,5) = l1*n2+n1*l2;
      
      T_strain(4,0) = 2*l2*l3;
      T_strain(4,1) = 2*m2*m3;
      T_strain(4,2) = 2*n2*n3;
      T_strain(4,3) = l2*m3+m2*l3;
      T_strain(4,4) = m2*n3+n2*m3;
      T_strain(4,5) = l2*n3+n2*l3;
      
      T_strain(5,0) = 2*l1*l3;
      T_strain(5,1) = 2*m1*m3;
      T_strain(5,2) = 2*n1*n3;
      T_strain(5,3) = l1*m3+m1*l3;
      T_strain(5,4) = m1*n3+n1*m3;
      T_strain(5,5) = l1*n3+n1*l3;
      
      ublas::matrix<FieldData> T_strain_r_Theta;    T_strain_r_Theta.resize(6,6); T_strain_r_Theta.clear();
      T_strain_r_Theta(0,0) = 2*l1*l1_r;
      T_strain_r_Theta(0,1) = 2*m1*m1_r;
      T_strain_r_Theta(0,2) = 2*n1*n1_r;
      T_strain_r_Theta(0,3) = l1_r*m1 + l1*m1_r;
      T_strain_r_Theta(0,4) = m1_r*n1 + m1*n1_r;
      T_strain_r_Theta(0,5) = l1_r*n1 + l1*n1_r;
      
      T_strain_r_Theta(1,0) = 2*l2*l2_r;
      T_strain_r_Theta(1,1) = 2*m2*m2_r;
      T_strain_r_Theta(1,2) = 2*n2*n2_r;
      T_strain_r_Theta(1,3) = l2_r*m2 + l2*m2_r;
      T_strain_r_Theta(1,4) = m2_r*n2 + m2*n2_r;
      T_strain_r_Theta(1,5) = l2_r*n2 + l2*n2_r;
      
      T_strain_r_Theta(2,0) = 2*l3*l3_r;
      T_strain_r_Theta(2,1) = 2*m3*m3_r;
      T_strain_r_Theta(2,2) = 2*n3*n3_r;
      T_strain_r_Theta(2,3) = l3_r*m3 + l3*m3_r;
      T_strain_r_Theta(2,4) = m3_r*n3 + m3*n3_r;
      T_strain_r_Theta(2,5) = l3_r*n3 + l3*n3_r;
      
      T_strain_r_Theta(3,0) = 2*(l1_r*l2 + l1*l2_r);
      T_strain_r_Theta(3,1) = 2*(m1_r*m2 + m1*m2_r);
      T_strain_r_Theta(3,2) = 2*(n1_r*n2 + n1*n2_r);
      T_strain_r_Theta(3,3) = l1_r*m2 + l1*m2_r + m1_r*l2 + m1*l2_r;
      T_strain_r_Theta(3,4) = m1_r*n2 + m1*n2_r + n1_r*m2 + n1*m2_r;
      T_strain_r_Theta(3,5) = l1_r*n2 + l1*n2_r + n1_r*l2 + n1*l2_r;
      
      T_strain_r_Theta(4,0) = 2*(l2_r*l3 + l2*l3_r);
      T_strain_r_Theta(4,1) = 2*(m2_r*m3 + m2*m3_r);
      T_strain_r_Theta(4,2) = 2*(n2_r*n3 + n2*n3_r);
      T_strain_r_Theta(4,3) = l2_r*m3 + l2*m3_r + m2_r*l3 + m2*l3_r;
      T_strain_r_Theta(4,4) = m2_r*n3 + m2*n3_r + n2_r*m3 + n2*m3_r;
      T_strain_r_Theta(4,5) = l2_r*n3 + l2*n3_r + n2_r*l3 + n2*l3_r;
      
      T_strain_r_Theta(5,0) = 2*(l1_r*l3 + l1*l3_r);
      T_strain_r_Theta(5,1) = 2*(m1_r*m3 + m1*m3_r);
      T_strain_r_Theta(5,2) = 2*(n1_r*n3 + n1*n3_r);
      T_strain_r_Theta(5,3) = l1_r*m3 + l1*m3_r + m1_r*l3 + m1*l3_r;
      T_strain_r_Theta(5,4) = m1_r*n3 + m1*n3_r + n1_r*m3 + n1*m3_r;
      T_strain_r_Theta(5,5) = l1_r*n3 + l1*n3_r + n1_r*l3 + n1*l3_r;
      
      T_strain_r_Theta = T_strain_r_Theta*(M_PI/180);
      
      //cout<<"\n\nT_strain = "<<T_strain<<endl;
      
      ublas::matrix<FieldData> Mat1 = prod(Dmat_123,T_strain);
      ublas::matrix<FieldData> Mat2 = prod(trans(T_strain_r_Theta), Mat1);
      
      ublas::matrix<FieldData> Mat3 = prod(Dmat_123,T_strain_r_Theta);
      ublas::matrix<FieldData> Mat4 = prod(trans(T_strain), Mat3);
      
      Dmat_xyz = Mat2 + Mat4;
      
      
      PetscFunctionReturn(0);
    }
    
    
    // =========================================================================
    //
    //  A.VI. SOLUTION PHASE:
    //        Caculate RVE constitutive matrix Dmat
    //        Computational Homogenization Method
    //
    // =========================================================================
    
    PetscErrorCode Micro_FE_Dmat(FieldInterface &m_field_RVE,
                                 int &nvars, int &nders,
                                 vector<string> &stochastic_fields,
                                 ublas::vector<double> matprop,
                                 int num_rvars,
                                 vector<string> vars_name) {
      
      PetscFunctionBegin;
      cout <<"Hi from Calculate_RVEDmat"<<endl;
      
      //ErrorCode rval;
      PetscErrorCode ierr;
      
      Dmat.resize(6,6); Dmat.clear();
      
      Dmat_r_Em.resize(6,6);   Dmat_r_Em.clear();
      Dmat_r_NUm.resize(6,6);  Dmat_r_NUm.clear();
      Dmat_r_Ep.resize(6,6);   Dmat_r_Ep.clear();
      Dmat_r_Ez.resize(6,6);   Dmat_r_Ez.clear();
      Dmat_r_NUp.resize(6,6);  Dmat_r_NUp.clear();
      Dmat_r_NUpz.resize(6,6); Dmat_r_NUpz.clear();
      Dmat_r_Gzp.resize(6,6);  Dmat_r_Gzp.clear();
      Dmat_r_Ef.resize(6,6);   Dmat_r_Ef.clear();
      Dmat_r_NUf.resize(6,6);  Dmat_r_NUf.clear();
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      Vec F1,F2,F3,F4,F5,F6,D;
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F6); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D); CHKERRQ(ierr);
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat Aij;
      ierr = m_field_RVE.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_RVE",&Aij); CHKERRQ(ierr);
      
      struct MyElasticFEMethod: public ElasticFEMethod {
        MyElasticFEMethod(FieldInterface& _m_field,
                          Mat& _Aij,Vec& _D,Vec& _F,double _lambda,double _mu, string _field_name = "DISPLACEMENT"):
        ElasticFEMethod(_m_field,_Aij,_D,_F,_lambda,_mu,_field_name) {};
        
        PetscErrorCode Fint(Vec F_int) {
          PetscFunctionBegin;
          ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
          for(int rr = 0;rr<row_mat;rr++) {
            if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
            if(RowGlob[rr].size()==0) continue;
            f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
            ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
          }
          PetscFunctionReturn(0);
        }
      };
      
      //Assemble F and Aij
      double YoungModulus = 3500;
      double PoissonRatio = 0.3;
      double YoungModulus_Fibre, PoissonRatio_Fibre;
      //double alpha;
      int field_rank=3;
      
      /*************************************************************************
       *
       *  2. Get the volume of RVE
       *
       ************************************************************************/
      double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
      Vec RVE_volume_Vec;
      ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field_RVE.get_moab(),MYPCOMM_INDEX);
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
      ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
      
      RVEVolume MyRVEVol(m_field_RVE,Aij,D,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
      RVEVolumeTrans MyRVEVolTrans(m_field_RVE,Aij,D,F1, RVE_volume_Vec);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVolTrans);  CHKERRQ(ierr);
      //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
      cout<<"Final RVE_volume = "<< RVE_volume <<endl;
      
      
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)) {
        cout << endl << *it << endl;
        
        //Get block name
        string name = it->get_name();
        // ---------------------------------
        //
        // Modify matrix material properties
        //
        // ---------------------------------
        if (name.compare(0,13,"MAT_ELASTIC_1") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_rvars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Em") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = matprop(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUm") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = matprop(ii-1);
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Young   = matprop(0);
           mydata.data.Poisson = matprop(1);
           */
          YoungModulus=mydata.data.Young;
          PoissonRatio=mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Matrix material:\n" << mydata;
        }
        
        // ---------------------------------
        //
        // Modify matrix fibre properties
        // Case 1: Isotropic material
        //
        // ---------------------------------
        if (name.compare(0,19,"MAT_FIBRE_ISOTROPIC") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_rvars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Ef") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = matprop(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUf") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = matprop(ii-1);
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Young   = matprop(0);
           mydata.data.Poisson = matprop(1);
           */
          YoungModulus_Fibre = mydata.data.Young;
          PoissonRatio_Fibre = mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre material:\n" << mydata;
        }
        
        // ---------------------------------
        //
        // Modify fibre material properties
        // Case 2: Transversely isotropic material
        //
        // ---------------------------------
        
        if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
          Mat_Elastic_TransIso mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int i=1;i<=num_rvars;i++) {
            ParameterName = vars_name[i];
            cout<<"the variable name is "<<vars_name[i]<<endl;
            if (ParameterName.compare(0,2,"Ez") == 0) {
              mydata.data.Youngz = matprop(i-1);
            }
            else if (ParameterName.compare(0,2,"Ep") == 0) {
              mydata.data.Youngp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUp") == 0) {
              mydata.data.Poissonp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUz") == 0) {
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"Gzp") == 0) {
              mydata.data.Shearzp = matprop(i-1);
            }
            else if (ParameterName.compare(0,2,"Ef") == 0) {
              mydata.data.Youngz = matprop(i-1);
              mydata.data.Youngp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUf") == 0) {
              mydata.data.Poissonp  = matprop(i-1);
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if ((ParameterName.compare(0,3,"NUf") == 0) || (ParameterName.compare(0,2,"Ef") == 0)) {
              mydata.data.Shearzp = mydata.data.Youngz/(2*(1+mydata.data.Poissonp));
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Poissonp  = matprop(2);
           mydata.data.Poissonpz = matprop(3);
           mydata.data.Youngp    = matprop(4);
           mydata.data.Youngz    = matprop(5);
           mydata.data.Shearzp   = matprop(6);
           */
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre material:\n" << mydata;
        }
      }
      
      ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
      applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
      
      MyElasticFEMethod my_fe(m_field_RVE,Aij,D,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),"DISP_RVE");
      TranIsotropicFibreDirRotElasticFEMethod_PSFE MyTIsotFE(m_field_RVE,Aij,D,F1,"DISP_RVE");
      
      ElasticFE_RVELagrange_Disp_Multi_Rhs MyFE_RVELagrange(m_field_RVE,Aij,D,F1,F2,F3,F4,F5,F6,applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank);
      
      cout<<"After ElasticFE_RVELagrange_Disp_Multi_Rhs "<<endl;
      
      ierr = VecZeroEntries(F1); CHKERRQ(ierr);
      ierr = VecZeroEntries(F2); CHKERRQ(ierr);
      ierr = VecZeroEntries(F3); CHKERRQ(ierr);
      ierr = VecZeroEntries(F4); CHKERRQ(ierr);
      ierr = VecZeroEntries(F5); CHKERRQ(ierr);
      ierr = VecZeroEntries(F6); CHKERRQ(ierr);
      
      ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
      
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyTIsotFE);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVELagrange);  CHKERRQ(ierr);
      
      
      ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F1); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F1); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F2); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F2); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F3); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F3); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F4); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F4); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F5); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F5); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F6); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F6); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  3. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ************************************************************************/
      
      //Solver
      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      ierr = KSPSetUp(solver); CHKERRQ(ierr);
      
      //create a vector for 6 components of homogenized stress
      Vec Stress_Homo, Stress_Homo_r;
      PetscScalar *avec;
      
      if(pcomm->rank()==0) {
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo_r);
      } else {
        int ghost[] = {0,1,2,3,4,5};
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo_r);
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 1: applied macro strain: [1 0 0 0 0 0]^T
      //------------------------------------------------------------------------
      
      cout<<"===============================================================\n";
      cout<<"        Applied strain [1 0 0 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      
      //solve for F1 and D1
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = KSPSolve(solver,F1,D); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // calculate homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_1(m_field_RVE,Aij,D,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_1);  CHKERRQ(ierr);
      VecGetArray(Stress_Homo, &avec);
      for (int ii=0; ii<6; ii++) {
        Dmat(ii,0)=*avec;
        avec++;
      }
      
      /*if (pcomm->rank()==0) {
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,0)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo,&avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=1; ii <= num_rvars; ii++) {
        Vec dF, dD;
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF); CHKERRQ(ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD); CHKERRQ(ierr);
        int idx_disp = 0; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        
        ierr = VecZeroEntries(dF); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
        
        if (VariableName.compare(0,2,"Em") == 0) { // due to Young's modulus of matrix - isotropic
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r(m_field_RVE,Aij,D,dF,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUm") == 0) { // due to Poisson's ratio in matrix - isotropic
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUp") == 0) { // due to Poisson's ratio in p-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUz") == 0) { // due to Poisson's ratio in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ep") == 0) { // due to Young's modulus in p-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ez") == 0) { // due to Young's modulus in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) { // due to shear modulus in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D,dF,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ef") == 0) { // due to Young's modulus of fibre - isotropic material
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ef(m_field_RVE,Aij,D,dF,"DISP_RVE","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ef);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUf") == 0) { // due to Poisson's ratio of fibre - isotropic material
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUf(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUf);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if (idx_disp == 1) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[0];
          //------------
          // a. Solving
          //------------
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if (idx_disp == 1) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r1(m_field_RVE,Aij,dD,dF,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r1);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            for(int irow=0; irow<6; irow++){
              cout.precision(15);
              //cout<<*avec_r<<endl;
              if (VariableName.compare(0,2,"Em") == 0)       { Dmat_r_Em(irow,0) = *avec_r; }
              else if (VariableName.compare(0,3,"NUm") == 0) { Dmat_r_NUm(irow,0) = *avec_r; }
              else if (VariableName.compare(0,3,"NUp") == 0) { Dmat_r_NUp(irow,0) = *avec_r; }
              else if (VariableName.compare(0,3,"NUz") == 0) { Dmat_r_NUpz(irow,0) = *avec_r; }
              else if (VariableName.compare(0,2,"Ep") == 0)  { Dmat_r_Ep(irow,0) = *avec_r; }
              else if (VariableName.compare(0,2,"Ez") == 0)  { Dmat_r_Ez(irow,0) = *avec_r; }
              else if (VariableName.compare(0,3,"Gzp") == 0) { Dmat_r_Gzp(irow,0) = *avec_r; }
              else if (VariableName.compare(0,2,"Ef") == 0)  { Dmat_r_Ef(irow,0) = *avec_r; }
              else if (VariableName.compare(0,3,"NUf") == 0) { Dmat_r_NUf(irow,0) = *avec_r; }

              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        //cout<< "\n\n";
        VariableName.clear();
        ierr = VecDestroy(&dD); CHKERRQ(ierr);
        ierr = VecDestroy(&dF); CHKERRQ(ierr);
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 2: applied macro strain: [0 1 0 0 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 1 0 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      
      // solve for F2 and D2
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = KSPSolve(solver,F2,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_2(m_field_RVE,Aij,D,F2,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_2);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,1)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,1)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo,&avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=1; ii <= num_rvars; ii++) {
        Vec dF, dD;
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF); CHKERRQ(ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD); CHKERRQ(ierr);
        int idx_disp = 0; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        
        ierr = VecZeroEntries(dF); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        
        if (VariableName.compare(0,2,"Em") == 0) { // due to Young's modulus of matrix - isotropic
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r(m_field_RVE,Aij,D,dF,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUm") == 0) { // due to Poisson's ratio of matrix - isotropic
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUp") == 0) { // due to Poisson's ratio in p-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUz") == 0) { // due to Poisson's ratio in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ep") == 0) { // due to Young's modulus in p-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ez") == 0) { // due to Young's modulus in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) { // due to shear modulus in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D,dF,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ef") == 0) { // due to Young's modulus of fibre - isotropic material
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ef(m_field_RVE,Aij,D,dF,"DISP_RVE","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ef);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUf") == 0) { // due to Poisson's ratio of fibre - isotropic material
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUf(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUf);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if (idx_disp == 1) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[0];
          //------------
          // a. Solving
          //------------
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if (idx_disp == 1) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r2(m_field_RVE,Aij,dD,dF,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r2);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            for(int irow=0; irow<6; irow++){
              cout.precision(15);
              //cout<<*avec_r<<endl;
              if (VariableName.compare(0,2,"Em") == 0)       { Dmat_r_Em(irow,1) = *avec_r; }
              else if (VariableName.compare(0,3,"NUm") == 0) { Dmat_r_NUm(irow,1) = *avec_r; }
              else if (VariableName.compare(0,3,"NUp") == 0) { Dmat_r_NUp(irow,1) = *avec_r; }
              else if (VariableName.compare(0,3,"NUz") == 0) { Dmat_r_NUpz(irow,1) = *avec_r; }
              else if (VariableName.compare(0,2,"Ep") == 0)  { Dmat_r_Ep(irow,1) = *avec_r; }
              else if (VariableName.compare(0,2,"Ez") == 0)  { Dmat_r_Ez(irow,1) = *avec_r; }
              else if (VariableName.compare(0,3,"Gzp") == 0) { Dmat_r_Gzp(irow,1) = *avec_r; }
              else if (VariableName.compare(0,2,"Ef") == 0)  { Dmat_r_Ef(irow,1) = *avec_r; }
              else if (VariableName.compare(0,3,"NUf") == 0) { Dmat_r_NUf(irow,1) = *avec_r; }
              
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        VariableName.clear();
        ierr = VecDestroy(&dD); CHKERRQ(ierr);
        ierr = VecDestroy(&dF); CHKERRQ(ierr);
        //cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 3: applied macro strain: [0 0 1 0 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 1 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F3 and D3
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = KSPSolve(solver,F3,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_3(m_field_RVE,Aij,D,F3,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_3);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,2)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,2)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=1; ii <= num_rvars; ii++) {
        Vec dF, dD;
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF); CHKERRQ(ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD); CHKERRQ(ierr);
        int idx_disp = 0; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];

        ierr = VecZeroEntries(dF); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        
        if (VariableName.compare(0,2,"Em") == 0) { // due to Young's modulus of matrix - isotropic
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r(m_field_RVE,Aij,D,dF,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUm") == 0) { // due to Poisson's ratio of matrix - isotropic
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUp") == 0) { // due to Poisson's ratio in p-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUz") == 0) { // due to Poisson's ratio in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ep") == 0) { // due to Young's modulus in p-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ez") == 0) { // due to Young's modulus in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) { // due to shear modulus in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D,dF,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ef") == 0) { // due to Young's modulus of fibre - isotropic material
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ef(m_field_RVE,Aij,D,dF,"DISP_RVE","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ef);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUf") == 0) { // due to Poisson's ratio of fibre - isotropic material
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUf(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUf);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if (idx_disp == 1) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[0];
          //------------
          // a. Solving
          //------------
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if (idx_disp==1) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r3(m_field_RVE,Aij,dD,dF,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r3);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              if (VariableName.compare(0,2,"Em") == 0)       { Dmat_r_Em(irow,2) = *avec_r; }
              else if (VariableName.compare(0,3,"NUm") == 0) { Dmat_r_NUm(irow,2) = *avec_r; }
              else if (VariableName.compare(0,3,"NUp") == 0) { Dmat_r_NUp(irow,2) = *avec_r; }
              else if (VariableName.compare(0,3,"NUz") == 0) { Dmat_r_NUpz(irow,2) = *avec_r; }
              else if (VariableName.compare(0,2,"Ep") == 0)  { Dmat_r_Ep(irow,2) = *avec_r; }
              else if (VariableName.compare(0,2,"Ez") == 0)  { Dmat_r_Ez(irow,2) = *avec_r; }
              else if (VariableName.compare(0,3,"Gzp") == 0) { Dmat_r_Gzp(irow,2) = *avec_r; }
              else if (VariableName.compare(0,2,"Ef") == 0)  { Dmat_r_Ef(irow,2) = *avec_r; }
              else if (VariableName.compare(0,3,"NUf") == 0) { Dmat_r_NUf(irow,2) = *avec_r; }
              
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        VariableName.clear();
        ierr = VecDestroy(&dD); CHKERRQ(ierr);
        ierr = VecDestroy(&dF); CHKERRQ(ierr);
        //cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 4: applied macro strain: [0 0 0 1 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 0 1 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F4 and D4
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = KSPSolve(solver,F4,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      
      // Extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_4(m_field_RVE,Aij,D,F4,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_4);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,3)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,3)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=1; ii <= num_rvars; ii++) {
        Vec dF, dD;
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF); CHKERRQ(ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD); CHKERRQ(ierr);
        int idx_disp = 0; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
      
        ierr = VecZeroEntries(dF); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        
        if (VariableName.compare(0,2,"Em") == 0) { // due to Young's modulus of matrix - isotropic
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r(m_field_RVE,Aij,D,dF,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUm") == 0) { // due to Poisson's ratio of matrix - isotropic
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUp") == 0) { // due to Poisson's ratio in p-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUz") == 0) { // due to Poisson's ratio in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ep") == 0) { // due to Young's modulus in p-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ez") == 0) { // due to Young's modulus in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) { // due to shear modulus in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D,dF,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ef") == 0) { // due to Young's modulus of fibre - isotropic material
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ef(m_field_RVE,Aij,D,dF,"DISP_RVE","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ef);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUf") == 0) { // due to Poisson's ratio of fibre - isotropic material
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUf(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUf);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if (idx_disp==1) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[0];
          //------------
          // a. Solving
          //------------
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if (idx_disp==1) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r4(m_field_RVE,Aij,dD,dF,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r4);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              if (VariableName.compare(0,2,"Em") == 0)       { Dmat_r_Em(irow,3) = *avec_r; }
              else if (VariableName.compare(0,3,"NUm") == 0) { Dmat_r_NUm(irow,3) = *avec_r; }
              else if (VariableName.compare(0,3,"NUp") == 0) { Dmat_r_NUp(irow,3) = *avec_r; }
              else if (VariableName.compare(0,3,"NUz") == 0) { Dmat_r_NUpz(irow,3) = *avec_r; }
              else if (VariableName.compare(0,2,"Ep") == 0)  { Dmat_r_Ep(irow,3) = *avec_r; }
              else if (VariableName.compare(0,2,"Ez") == 0)  { Dmat_r_Ez(irow,3) = *avec_r; }
              else if (VariableName.compare(0,3,"Gzp") == 0) { Dmat_r_Gzp(irow,3) = *avec_r; }
              else if (VariableName.compare(0,2,"Ef") == 0)  { Dmat_r_Ef(irow,3) = *avec_r; }
              else if (VariableName.compare(0,3,"NUf") == 0) { Dmat_r_NUf(irow,3) = *avec_r; }
    
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        VariableName.clear();
        ierr = VecDestroy(&dD); CHKERRQ(ierr);
        ierr = VecDestroy(&dF); CHKERRQ(ierr);
        //cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 5: applied macro strain: [0 0 0 0 1 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 0 0 1 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F5 and D5
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = KSPSolve(solver,F5,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_5(m_field_RVE,Aij,D,F5,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_5);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,4)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,4)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=1; ii <= num_rvars; ii++) {
        Vec dF, dD;
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF); CHKERRQ(ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD); CHKERRQ(ierr);
        int idx_disp = 0; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        
        ierr = VecZeroEntries(dF); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        
        if (VariableName.compare(0,2,"Em") == 0) { // due to Young's modulus of matrix - isotropic
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r(m_field_RVE,Aij,D,dF,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUm") == 0) { // due to Poisson's ratio of matrix - isotropic
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUp") == 0) { // due to Poisson's ratio in p-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUz") == 0) { // due to Poisson's ratio in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ep") == 0) { // due to Young's modulus in p-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUz") == 0) { // due to Young's modulus in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) { // due to shear modulus in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D,dF,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ef") == 0) { // due to Young's modulus of fibre - isotropic material
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ef(m_field_RVE,Aij,D,dF,"DISP_RVE","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ef);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUf") == 0) { // due to Poisson's ratio of fibre - isotropic material
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUf(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUf);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if (idx_disp==1) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[0];
          //------------
          // a. Solving
          //------------
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if (idx_disp==1) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r5(m_field_RVE,Aij,dD,dF,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r5);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              if (VariableName.compare(0,2,"Em") == 0)       { Dmat_r_Em(irow,4) = *avec_r; }
              else if (VariableName.compare(0,3,"NUm") == 0) { Dmat_r_NUm(irow,4) = *avec_r; }
              else if (VariableName.compare(0,3,"NUp") == 0) { Dmat_r_NUp(irow,4) = *avec_r; }
              else if (VariableName.compare(0,3,"NUz") == 0) { Dmat_r_NUpz(irow,4) = *avec_r; }
              else if (VariableName.compare(0,2,"Ep") == 0)  { Dmat_r_Ep(irow,4) = *avec_r; }
              else if (VariableName.compare(0,2,"Ez") == 0)  { Dmat_r_Ez(irow,4) = *avec_r; }
              else if (VariableName.compare(0,3,"Gzp") == 0) { Dmat_r_Gzp(irow,4) = *avec_r; }
              else if (VariableName.compare(0,2,"Ef") == 0)  { Dmat_r_Ef(irow,4) = *avec_r; }
              else if (VariableName.compare(0,3,"NUf") == 0) { Dmat_r_NUf(irow,4) = *avec_r; }
              
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        VariableName.clear();
        ierr = VecDestroy(&dD); CHKERRQ(ierr);
        ierr = VecDestroy(&dF); CHKERRQ(ierr);
        //cout<< "\n\n";
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 6: applied macro strain: [0 0 0 0 0 1]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"         Applied strain [0 0 0 0 0 1]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F6 and D6
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = KSPSolve(solver,F6,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_6(m_field_RVE,Aij,D,F6,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_6);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,5)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,5)<<endl;
       }
       cout<< "\n\n";
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=1; ii <= num_rvars; ii++) {
        Vec dF, dD;
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF); CHKERRQ(ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD); CHKERRQ(ierr);
        int idx_disp = 0; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
  
        ierr = VecZeroEntries(dF); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        
        if (VariableName.compare(0,2,"Em") == 0) { // due to Young's modulus of matrix - isotropic
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r(m_field_RVE,Aij,D,dF,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUm") == 0) { // due to Poisson's ratio of matrix - isotropic
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUp") == 0) { // due to Poisson's ratio in p-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUz") == 0) { // due to Poisson's ratio in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ep") == 0) { // due to Young's modulus in p-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ez") == 0) { // due to Young's modulus in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) { // due to shear modulus in z-direction of fibre
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D,dF,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ef") == 0) { // due to Young's modulus of fibre - isotropic material
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ef(m_field_RVE,Aij,D,dF,"DISP_RVE","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ef);  CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUf") == 0) { // due to Poisson's ratio of fibre - isotropic material
          idx_disp = 1;
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUf(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUf);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if (idx_disp==1) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[0];
          //------------
          // a. Solving
          //------------
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if (idx_disp==1) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r6(m_field_RVE,Aij,dD,dF,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r6);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              if (VariableName.compare(0,2,"Em") == 0)       { Dmat_r_Em(irow,5) = *avec_r; }
              else if (VariableName.compare(0,3,"NUm") == 0) { Dmat_r_NUm(irow,5) = *avec_r; }
              else if (VariableName.compare(0,3,"NUp") == 0) { Dmat_r_NUp(irow,5) = *avec_r; }
              else if (VariableName.compare(0,3,"NUz") == 0) { Dmat_r_NUpz(irow,5) = *avec_r; }
              else if (VariableName.compare(0,2,"Ep") == 0)  { Dmat_r_Ep(irow,5) = *avec_r; }
              else if (VariableName.compare(0,2,"Ez") == 0)  { Dmat_r_Ez(irow,5) = *avec_r; }
              else if (VariableName.compare(0,3,"Gzp") == 0) { Dmat_r_Gzp(irow,5) = *avec_r; }
              else if (VariableName.compare(0,2,"Ef") == 0)  { Dmat_r_Ef(irow,5) = *avec_r; }
              else if (VariableName.compare(0,3,"NUf") == 0) { Dmat_r_NUf(irow,5) = *avec_r; }

              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        VariableName.clear();
        ierr = VecDestroy(&dD); CHKERRQ(ierr);
        ierr = VecDestroy(&dF); CHKERRQ(ierr);
        //cout<< "\n\n";
      }
      
      /*************************************************************************
       *
       *  4. FINISH
       *
       ************************************************************************/
      
      //Destroy matrices and vectors
      ierr = VecDestroy(&F1); CHKERRQ(ierr);
      ierr = VecDestroy(&F2); CHKERRQ(ierr);
      ierr = VecDestroy(&F3); CHKERRQ(ierr);
      ierr = VecDestroy(&F4); CHKERRQ(ierr);
      ierr = VecDestroy(&F5); CHKERRQ(ierr);
      ierr = VecDestroy(&F6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      
      ierr = VecDestroy(&Stress_Homo); CHKERRQ(ierr);
      ierr = VecDestroy(&Stress_Homo_r); CHKERRQ(ierr);
      ierr = VecDestroy(&RVE_volume_Vec); CHKERRQ(ierr);
      
      ierr = MatDestroy(&Aij); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      
      //cout<<"\n\n";
      //m_field_RVE.list_fields();
      
      
      PetscFunctionReturn(0);
    }
    
    // =========================================================================
    //
    //  A.VI. SOLUTION PHASE:
    //        Caculate RVE constitutive matrix Dmat
    //        Rule of Mixture
    //
    // =========================================================================
    
    PetscErrorCode Dmat_RuleMixture(FieldInterface &m_field_RVE,
                                    int &nvars, int &nders,
                                    vector<string> &stochastic_fields,
                                    ublas::vector<double> matprop,
                                    int num_rvars,
                                    vector<string> vars_name) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      
      Dmat.resize(6,6); Dmat.clear();
      
      Dmat_r_Em.resize(6,6);   Dmat_r_Em.clear();
      Dmat_r_NUm.resize(6,6);  Dmat_r_NUm.clear();
      Dmat_r_Ep.resize(6,6);   Dmat_r_Ep.clear();
      Dmat_r_Ez.resize(6,6);   Dmat_r_Ez.clear();
      Dmat_r_NUp.resize(6,6);  Dmat_r_NUp.clear();
      Dmat_r_NUpz.resize(6,6); Dmat_r_NUpz.clear();
      Dmat_r_Gzp.resize(6,6);  Dmat_r_Gzp.clear();
      Dmat_r_Ef.resize(6,6);   Dmat_r_Ef.clear();
      Dmat_r_NUf.resize(6,6);  Dmat_r_NUf.clear();
      
      double Em, NUm;
      double Ef, NUf;
      double Ep, Ez, NUp, NUpz, Gzp;
      double Vf;
      
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)) {
        cout << endl << *it << endl;
        
        //Get block name
        string name = it->get_name();
        // ---------------------------------
        //
        // Modify matrix material properties
        //
        // ---------------------------------
        if (name.compare(0,13,"MAT_ELASTIC_1") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_rvars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Em") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = matprop(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUm") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = matprop(ii-1);
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Young   = matprop(0);
           mydata.data.Poisson = matprop(1);
           */
          Em  = mydata.data.Young;
          NUm = mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Matrix material:\n" << mydata;
        }
        
        // ---------------------------------
        //
        // Modify matrix fibre properties
        // Case 1: Isotropic material
        //
        // ---------------------------------
        if (name.compare(0,19,"MAT_FIBRE_ISOTROPIC") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_rvars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Ef") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = matprop(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUf") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = matprop(ii-1);
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Young   = matprop(0);
           mydata.data.Poisson = matprop(1);
           */
          Ef  = mydata.data.Young;
          NUf = mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre material:\n" << mydata;
        }
        
        // ---------------------------------
        //
        // Modify fibre material properties
        // Case 2: Transversely isotropic material
        //
        // ---------------------------------
        
        if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
          Mat_Elastic_TransIso mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int i=1;i<=num_rvars;i++) {
            ParameterName = vars_name[i];
            cout<<"the variable name is "<<vars_name[i]<<endl;
            if (ParameterName.compare(0,2,"Ez") == 0) {
              mydata.data.Youngz = matprop(i-1);
            }
            else if (ParameterName.compare(0,2,"Ep") == 0) {
              mydata.data.Youngp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUp") == 0) {
              mydata.data.Poissonp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUz") == 0) {
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"Gzp") == 0) {
              mydata.data.Shearzp = matprop(i-1);
            }
            else if (ParameterName.compare(0,2,"Ef") == 0) {
              mydata.data.Youngz = matprop(i-1);
              mydata.data.Youngp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUf") == 0) {
              mydata.data.Poissonp  = matprop(i-1);
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if ((ParameterName.compare(0,3,"NUf") == 0) || (ParameterName.compare(0,2,"Ef") == 0)) {
              mydata.data.Shearzp = mydata.data.Youngz/(2*(1+mydata.data.Poissonp));
            }
            ParameterName.clear();
          }
          
          NUp  = mydata.data.Poissonp;
          NUpz = mydata.data.Poissonpz;
          Ep   = mydata.data.Youngp;
          Ez   = mydata.data.Youngz;
          Gzp  =  mydata.data.Shearzp;
           
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre material:\n" << mydata;
        }
      }
      
      // Use the rule of mixture to estimate effective material properties
      TransverseIsotropicStiffnessMatrix_RuleOfMixture RuleOfMixture;
      Vf = 0.6;
      RuleOfMixture.EEP_RM(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
      Dmat = RuleOfMixture.StiffnessMatrix;
      
      // Calculate the first-order derivative of the effective elastic matrix
      string ParameterName;
      for (int i=1;i<=num_rvars;i++) {
        ParameterName = vars_name[i];
        if (ParameterName.compare(0,3,"NUp") == 0) {
          RuleOfMixture.D_r_PoissonP(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_NUp = RuleOfMixture.StiffnessMatrix_r_NUpf;
        }
        else if (ParameterName.compare(0,3,"NUz") == 0) {
          RuleOfMixture.D_r_PoissonPZ(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_NUpz = RuleOfMixture.StiffnessMatrix_r_NUpzf;
        }
        else if (ParameterName.compare(0,2,"Ep") == 0) {
          RuleOfMixture.D_r_YoungP(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_Ep = RuleOfMixture.StiffnessMatrix_r_Epf;
        }
        else if (ParameterName.compare(0,2,"Ez") == 0) {
          RuleOfMixture.D_r_YoungZ(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_Ez = RuleOfMixture.StiffnessMatrix_r_Ezf;
        }
        else if (ParameterName.compare(0,3,"Gzp") == 0) {
          RuleOfMixture.D_r_ShearZP(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_Gzp = RuleOfMixture.StiffnessMatrix_r_Gzpf;
        }
        else if (ParameterName.compare(0,2,"Em") == 0) {
          RuleOfMixture.D_r_YoungM(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_Em = RuleOfMixture.StiffnessMatrix_r_Em;
        }
        else if (ParameterName.compare(0,3,"NUm") == 0) {
          RuleOfMixture.D_r_PoissonM(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_NUm = RuleOfMixture.StiffnessMatrix_r_NUm;
        }
      }
      
      PetscFunctionReturn(0);
    }

    // =========================================================================
    //
    //  A.VI. SOLUTION PHASE:
    //        Caculate RVE constitutive matrix Dmat
    //        Mori-Tanaka method
    //
    // =========================================================================
    
    PetscErrorCode Dmat_MoriTanaka(FieldInterface &m_field_RVE,
                                    int &nvars, int &nders,
                                    vector<string> &stochastic_fields,
                                    ublas::vector<double> matprop,
                                    int num_rvars,
                                    vector<string> vars_name) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      
      Dmat.resize(6,6); Dmat.clear();
      
      Dmat_r_Em.resize(6,6);   Dmat_r_Em.clear();
      Dmat_r_NUm.resize(6,6);  Dmat_r_NUm.clear();
      Dmat_r_Ep.resize(6,6);   Dmat_r_Ep.clear();
      Dmat_r_Ez.resize(6,6);   Dmat_r_Ez.clear();
      Dmat_r_NUp.resize(6,6);  Dmat_r_NUp.clear();
      Dmat_r_NUpz.resize(6,6); Dmat_r_NUpz.clear();
      Dmat_r_Gzp.resize(6,6);  Dmat_r_Gzp.clear();
      Dmat_r_Ef.resize(6,6);   Dmat_r_Ef.clear();
      Dmat_r_NUf.resize(6,6);  Dmat_r_NUf.clear();
      
      double Em, NUm;
      double Ef, NUf;
      double Ep, Ez, NUp, NUpz, Gzp;
      double Vf;
      
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)) {
        cout << endl << *it << endl;
        
        //Get block name
        string name = it->get_name();
        // ---------------------------------
        //
        // Modify matrix material properties
        //
        // ---------------------------------
        if (name.compare(0,13,"MAT_ELASTIC_1") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_rvars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Em") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = matprop(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUm") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = matprop(ii-1);
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Young   = matprop(0);
           mydata.data.Poisson = matprop(1);
           */
          Em  = mydata.data.Young;
          NUm = mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Matrix material:\n" << mydata;
        }
        
        // ---------------------------------
        //
        // Modify matrix fibre properties
        // Case 1: Isotropic material
        //
        // ---------------------------------
        if (name.compare(0,19,"MAT_FIBRE_ISOTROPIC") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_rvars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Ef") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = matprop(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUf") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = matprop(ii-1);
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Young   = matprop(0);
           mydata.data.Poisson = matprop(1);
           */
          Ef  = mydata.data.Young;
          NUf = mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre material:\n" << mydata;
        }
        
        // ---------------------------------
        //
        // Modify fibre material properties
        // Case 2: Transversely isotropic material
        //
        // ---------------------------------
        
        if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
          Mat_Elastic_TransIso mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int i=1;i<=num_rvars;i++) {
            ParameterName = vars_name[i];
            cout<<"the variable name is "<<vars_name[i]<<endl;
            if (ParameterName.compare(0,2,"Ez") == 0) {
              mydata.data.Youngz = matprop(i-1);
            }
            else if (ParameterName.compare(0,2,"Ep") == 0) {
              mydata.data.Youngp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUp") == 0) {
              mydata.data.Poissonp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUz") == 0) {
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"Gzp") == 0) {
              mydata.data.Shearzp = matprop(i-1);
            }
            else if (ParameterName.compare(0,2,"Ef") == 0) {
              mydata.data.Youngz = matprop(i-1);
              mydata.data.Youngp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUf") == 0) {
              mydata.data.Poissonp  = matprop(i-1);
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if ((ParameterName.compare(0,3,"NUf") == 0) || (ParameterName.compare(0,2,"Ef") == 0)) {
              mydata.data.Shearzp = mydata.data.Youngz/(2*(1+mydata.data.Poissonp));
            }
            ParameterName.clear();
          }
          
          NUp  = mydata.data.Poissonp;
          NUpz = mydata.data.Poissonpz;
          Ep   = mydata.data.Youngp;
          Ez   = mydata.data.Youngz;
          Gzp  =  mydata.data.Shearzp;
           
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre material:\n" << mydata;
        }
      }
      
      // Use the rule of mixture to estimate effective material properties
      TransverseIsotropicStiffnessMatrix_MoriTanaka Mori_Tanaka;
      Vf = 0.6;
	  double lambda, mu;
	  lambda = Em*NUm/((1+NUm)*(1-2*NUm));
	  mu       =  Em/(2*(1+NUm));
      Mori_Tanaka.EEP_MT(NUp,NUpz,Ep,Ez,Gzp,lambda,mu,Vf);
      Dmat = Mori_Tanaka.StiffnessMatrix;
      
      // Calculate the first-order derivative of the effective elastic matrix
      string ParameterName;
      for (int i=1;i<=num_rvars;i++) {
        ParameterName = vars_name[i];
        if (ParameterName.compare(0,3,"NUp") == 0) {
          Mori_Tanaka.D_r_PoissonP(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_NUp = Mori_Tanaka.StiffnessMatrix_r_NUpf;
        }
        else if (ParameterName.compare(0,3,"NUz") == 0) {
          Mori_Tanaka.D_r_PoissonPZ(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_NUpz = Mori_Tanaka.StiffnessMatrix_r_NUpzf;
        }
        else if (ParameterName.compare(0,2,"Ep") == 0) {
          Mori_Tanaka.D_r_YoungP(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_Ep = Mori_Tanaka.StiffnessMatrix_r_Epf;
        }
        else if (ParameterName.compare(0,2,"Ez") == 0) {
          Mori_Tanaka.D_r_YoungZ(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_Ez = Mori_Tanaka.StiffnessMatrix_r_Ezf;
        }
        else if (ParameterName.compare(0,3,"Gzp") == 0) {
          Mori_Tanaka.D_r_ShearZP(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_Gzp = Mori_Tanaka.StiffnessMatrix_r_Gzpf;
        }
        else if (ParameterName.compare(0,2,"Em") == 0) {
          Mori_Tanaka.D_r_YoungM(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_Em = Mori_Tanaka.StiffnessMatrix_r_Em;
        }
        else if (ParameterName.compare(0,3,"NUm") == 0) {
          Mori_Tanaka.D_r_PoissonM(NUp,NUpz,Ep,Ez,Gzp,Em,NUm,Vf);
          Dmat_r_NUm = Mori_Tanaka.StiffnessMatrix_r_NUm;
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    // =========================================================================
    //
    //  B.VI. SOLUTION PHASE:
    //        Solve Macroscale FE equation
    //
    // =========================================================================
    
    PetscErrorCode Macro_FE_REL(FieldInterface &m_field_Macro,
                                int &nvars, int &nders,
                                vector<string> &stochastic_fields,
                                ublas::vector<double> TheVariables,
                                int num_rvars,
                                vector<string> vars_name,
                                ublas::vector<double> PlyAngle,
                                int ExaminedPly,
                                PetscInt NO_Layers,
                                ublas::vector<ublas::matrix <double> > &TheStress) {
      
      PetscFunctionBegin;
      
      //ErrorCode rval;
      PetscErrorCode ierr;
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      //create matrices
      
      Vec F, D;
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&F); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&D); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat A;
      ierr = m_field_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_MACRO",&A); CHKERRQ(ierr);
      
      struct MyElasticFEMethod_Macro: public FE2_ElasticFEMethod {
        MyElasticFEMethod_Macro(FieldInterface& _m_field_Macro,Mat _A,Vec _D,Vec& _F, ublas::matrix<FieldData> _Dmat,string _field_name):
        FE2_ElasticFEMethod(_m_field_Macro,_A,_D,_F, _Dmat, _field_name) {};
        
        virtual PetscErrorCode RhsAndLhs() {
          PetscFunctionBegin;
          
          ierr = Lhs(); CHKERRQ(ierr);
          
          PetscFunctionReturn(0);
        }
      };
      
      
      Projection10NodeCoordsOnField ent_method_material_Macro(m_field_Macro,"MESH_NODE_POSITIONS");
      ierr = m_field_Macro.loop_dofs("MESH_NODE_POSITIONS",ent_method_material_Macro); CHKERRQ(ierr);
      
      //Assemble F and A
      DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field_Macro,"DISP_MACRO",A,D,F);
      
      
      /*****************************************************************************
       *
       * Read the saved Dmat mechancial (from the computational homgenisaiton of the 0deg RVE)
       *
       ****************************************************************************/
      
      double theta;
      
      /*Dmat.clear();
      Dmat(0,0) = 1.2894e5; Dmat(0,1) = 0.0525e5; Dmat(0,2) = 0.0525e5;
      Dmat(1,0) = 0.0525e5; Dmat(1,1) = 0.1331e5; Dmat(1,2) = 0.0545e5;
      Dmat(2,0) = 0.0525e5; Dmat(2,1) = 0.0545e5; Dmat(2,2) = 0.1331e5;
      Dmat(3,3) = 0.0660e5; Dmat(4,4) = 0.0393e5; Dmat(5,5) = 0.0660e5;*/
      // First-layer
      theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
      Dmat_1st_Ply.resize(6,6);   Dmat_1st_Ply.clear();cout<<"\nOriginal Dmat: "<<Dmat<<endl;
      ierr = Dmat_Transformation_old(theta, Dmat, Dmat_1st_Ply); CHKERRQ(ierr);cout<<"\Transformmed Dmat: "<<Dmat_1st_Ply<<endl;
      ierr = Dmat_Transformation(theta, Dmat, Dmat_1st_Ply); CHKERRQ(ierr);cout<<"\Transformmed Dmat: "<<Dmat_1st_Ply<<endl;
      // Second-layer
      if (NO_Layers > 1) {
        theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
        Dmat_2nd_Ply.resize(6,6);   Dmat_2nd_Ply.clear();
        ierr = Dmat_Transformation(theta, Dmat, Dmat_2nd_Ply); CHKERRQ(ierr);//cout<<"\n\nLayer 2: \t"<<PlyAngle(1)<<"\n"<<Dmat_2nd_Ply<<endl;
      }
      // Third-layer
      if (NO_Layers > 2) {
        theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
        Dmat_3rd_Ply.resize(6,6);   Dmat_3rd_Ply.clear();
        ierr = Dmat_Transformation(theta, Dmat, Dmat_3rd_Ply); CHKERRQ(ierr);//cout<<"\n\nLayer 3: \t"<<PlyAngle(2)<<"\n"<<Dmat_3rd_Ply<<endl;
      }
      // Fourth-layer
      if (NO_Layers > 3) {
        theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
        Dmat_4th_Ply.resize(6,6);   Dmat_4th_Ply.clear();
        ierr = Dmat_Transformation(theta, Dmat, Dmat_4th_Ply); CHKERRQ(ierr);//cout<<"\n\nLayer 4: \t"<<PlyAngle(3)<<"\n"<<Dmat_4th_Ply<<endl;
      }
      
      MyElasticFEMethod_Macro my_fe_1st_ply (m_field_Macro,A,D,F,Dmat_1st_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_2nd_ply (m_field_Macro,A,D,F,Dmat_2nd_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_3rd_ply (m_field_Macro,A,D,F,Dmat_3rd_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_4th_ply (m_field_Macro,A,D,F,Dmat_4th_Ply,"DISP_MACRO");
      
      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = MatZeroEntries(A); CHKERRQ(ierr);
      
      //ierr = m_field_Macro.set_global_VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      //preproc
      ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
      //loop elems
      //PetscBarrier(PETSC_NULL);
      // First layer
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe_1st_ply);  CHKERRQ(ierr);
      // Second layer
      if (NO_Layers > 1) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe_2nd_ply);  CHKERRQ(ierr);
      }
      // Third layer
      if (NO_Layers > 2) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe_3rd_ply);  CHKERRQ(ierr);
      }
      // Fourth layer
      if (NO_Layers > 3) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe_4th_ply);  CHKERRQ(ierr);
      }
      
      
      //forces and preassures on surface
      boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
      
      // Check whether force is considered as random variable or not
      int idx_force = 100;
      for (int ii = 1; ii<=num_rvars; ii++) {
        string VariableName;
        VariableName = vars_name[ii];
        if (VariableName.compare(0,5,"force") == 0) {
          idx_force = ii;cout<<"The force is: "<<TheVariables(ii-1)<<endl;
        }
      }
      
      
      if (idx_force==100) {
        MetaNeummanForces Zeroth_FE;
        ierr = Zeroth_FE.setNeumannFiniteElementOperators(m_field_Macro,neumann_forces,F,"DISP_MACRO"); CHKERRQ(ierr);
      } else {
        MyMetaNeummanForces Zeroth_FE;
        ierr = Zeroth_FE.setNeumannFiniteElementOperators(m_field_Macro,
                                                          neumann_forces,F,
                                                          "DISP_MACRO",
                                                          "MESH_NODE_POSITIONS",
                                                          TheVariables(idx_force-1)); CHKERRQ(ierr);
      }
      
      boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
      for(;mit!=neumann_forces.end();mit++) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
      }
      
      //postproc
      ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
      
      //set matrix possitives define and symetric for cholesky and icc preceonditionser
      ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
      
      ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
      
      //VecCopy(ElemForce,F);
      
      
      /*************************************************************************
       *
       *  2. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ************************************************************************/
      //Solver
      KSP solver_Macro;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver_Macro); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver_Macro,A,A); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver_Macro); CHKERRQ(ierr);
      ierr = KSPSetUp(solver_Macro); CHKERRQ(ierr);
      
      //MatView(A,PETSC_VIEWER_STDOUT_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
      
      // elastic analys
      ierr = KSPSolve(solver_Macro,F,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      
      //Save data on mesh
      ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      //VecView(D,PETSC_VIEWER_STDOUT_WORLD);
      //VecView(F,PETSC_VIEWER_STDOUT_WORLD);
      
      cout<<"Solving the zeroth-order equation is finish. \n";
      
      ublas::matrix<double> StressGP(3,3);
      double Angle_REL;
      ublas::matrix<double> Dmat_REL; Dmat_REL.clear();
      Angle_REL = PlyAngle(ExaminedPly - 1)*(M_PI/180.0);
      
      switch (ExaminedPly) {
        case 1: { Dmat_REL = Dmat_1st_Ply; break; }
        case 2: { Dmat_REL = Dmat_2nd_Ply; break; }
        case 3: { Dmat_REL = Dmat_3rd_Ply; break; }
        case 4: { Dmat_REL = Dmat_4th_Ply; break; }
      }
      
      FE2_PostProcStressForReliability_Zeroth Calc_Stress(m_field_Macro,"DISP_MACRO",Dmat_REL);
      //cout<<"Dmat_REL: "<<Dmat_REL<<endl;
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress); CHKERRQ(ierr);
      //cout<<"\nStress "<<Calc_Stress.StressGP<<endl;
      StressGP.clear(); REL_Stress_Transformation(Angle_REL, Calc_Stress.StressGP, StressGP);
      TheStress(0).resize(3,3); TheStress(0).clear(); TheStress(0) = StressGP;
      //cout<<"Stress at GP in 123 in side FE2: "<<StressGP<<endl;
      
      /*************************************************************************
       *
       *  3. SOLVE THE FIRST-ORDER AND THE SECOND-ORDER FE EQUILIBRIUM EQUATION
       *     1st order-[K][U_r] = -[K_r][U}
       *     2nd order-[K][U_rs] = -[K_rs][U]-2[K_r][U_s]
       *
       ************************************************************************/
      int idx_disp; // index of displacement field in the field name vector <stochastic_fields>
      string VariableName;
      ublas::matrix<double> Dmat_REL_r;
      ublas::matrix<double> StressGP_r(3,3);
      for (int ii=1; ii<=num_rvars; ii++) {
        
        Vec dF, dD;
        ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&dF); CHKERRQ(ierr);
        ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&dD); CHKERRQ(ierr);
        idx_disp = 0; VariableName = vars_name[ii];
        Dmat_REL_r.resize(6,6); Dmat_REL_r.clear();
        
        ierr = VecZeroEntries(dD); CHKERRQ(ierr);
        ierr = VecZeroEntries(dF); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        
        /***********************************************************************
         *
         * 3.1. Case 1: Material properties are  treated as random variables
         *
         **********************************************************************/
       
        if (VariableName.compare(0,2,"Em") == 0) {// due to Young's modulus of matrix (Em)
          idx_disp = 1;
          // Examined ply
          ierr = Dmat_Transformation(Angle_REL, Dmat_r_Em, Dmat_REL_r); CHKERRQ(ierr);
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_Em.resize(6,6);   Dmat_1st_Ply_r_Em.clear();
          ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_1st_Ply_r_Em); CHKERRQ(ierr);
          // Second-layer
          if (NO_Layers > 1) {
            theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_2nd_Ply_r_Em.resize(6,6);   Dmat_2nd_Ply_r_Em.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_2nd_Ply_r_Em); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 2) {
            theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_3rd_Ply_r_Em.resize(6,6);   Dmat_3rd_Ply_r_Em.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_3rd_Ply_r_Em); CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_4th_Ply_r_Em.resize(6,6);   Dmat_4th_Ply_r_Em.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_4th_Ply_r_Em); CHKERRQ(ierr);
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_r_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Em,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Em,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Em,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Em,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_1st_Ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_2nd_Ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_3rd_Ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_4th_Ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r); CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUm") == 0) { // due to Poisson's ratio of matrix (NUm)
          idx_disp = 1;
          // Examined ply
          ierr = Dmat_Transformation(Angle_REL, Dmat_r_NUm, Dmat_REL_r); CHKERRQ(ierr);
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_NUm.resize(6,6);   Dmat_1st_Ply_r_NUm.clear();
          ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_1st_Ply_r_NUm); CHKERRQ(ierr);
          // Second-layer
          if (NO_Layers > 1) {
            theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_2nd_Ply_r_NUm.resize(6,6);   Dmat_2nd_Ply_r_NUm.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_2nd_Ply_r_NUm); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 2) {
            theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_3rd_Ply_r_NUm.resize(6,6);   Dmat_3rd_Ply_r_NUm.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_3rd_Ply_r_NUm); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 3) {
            theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_4th_Ply_r_NUm.resize(6,6);   Dmat_4th_Ply_r_NUm.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_4th_Ply_r_NUm); CHKERRQ(ierr);
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_NUm,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_NUm,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_NUm,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_NUm,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUm(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUm_1st_Ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_NUm_2nd_Ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_NUm_3rd_Ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_NUm_4th_Ply);  CHKERRQ(ierr);
          }
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {// due to transversal Poisson's ratio of fibre (NUp)
          idx_disp = 1;
          // Examined ply
          ierr = Dmat_Transformation(Angle_REL, Dmat_r_NUp, Dmat_REL_r); CHKERRQ(ierr);
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_NUp.resize(6,6);   Dmat_1st_Ply_r_NUp.clear();
          ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_1st_Ply_r_NUp); CHKERRQ(ierr);
          // Second-layer
          if (NO_Layers > 1) {
            theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_2nd_Ply_r_NUp.resize(6,6);   Dmat_2nd_Ply_r_NUp.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_2nd_Ply_r_NUp); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 2) {
            theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_3rd_Ply_r_NUp.resize(6,6);   Dmat_3rd_Ply_r_NUp.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_3rd_Ply_r_NUp); CHKERRQ(ierr);
          }
          // Fourth-layer
          if (NO_Layers > 3) {
            theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_4th_Ply_r_NUp.resize(6,6);   Dmat_4th_Ply_r_NUp.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_4th_Ply_r_NUp); CHKERRQ(ierr);
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_NUp,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_NUp,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_NUp,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_NUp,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUp(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUp_1st_Ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_NUp_2nd_Ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_NUp_3rd_Ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_NUp_4th_Ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUz") == 0) {// due to axial Poisson's ratio of fibre (NUpz)
          idx_disp = 1;
          // Examined ply
          ierr = Dmat_Transformation(Angle_REL, Dmat_r_NUpz, Dmat_REL_r); CHKERRQ(ierr);
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_NUpz.resize(6,6);   Dmat_1st_Ply_r_NUpz.clear();
          ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_1st_Ply_r_NUpz); CHKERRQ(ierr);
          // Second-layer
          if (NO_Layers > 1) {
            theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_2nd_Ply_r_NUpz.resize(6,6);   Dmat_2nd_Ply_r_NUpz.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_2nd_Ply_r_NUpz); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 2) {
            theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_3rd_Ply_r_NUpz.resize(6,6);   Dmat_3rd_Ply_r_NUpz.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_3rd_Ply_r_NUpz); CHKERRQ(ierr);
          }
          // Fourth-layer
          if (NO_Layers > 3) {
            theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_4th_Ply_r_NUpz.resize(6,6);   Dmat_4th_Ply_r_NUpz.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_4th_Ply_r_NUpz); CHKERRQ(ierr);
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_NUpz,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_NUpz,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_NUpz,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_NUpz,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUpz(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUpz_1st_Ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_NUpz_2nd_Ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_NUpz_3rd_Ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_NUpz_4th_Ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {// due to transversal modulus of fibre (Ep)
          idx_disp = 1;
          // Examined ply
          ierr = Dmat_Transformation(Angle_REL, Dmat_r_Ep, Dmat_REL_r); CHKERRQ(ierr);
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_Ep.resize(6,6);   Dmat_1st_Ply_r_Ep.clear();
          ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_1st_Ply_r_Ep); CHKERRQ(ierr);
          // Second-layer
          if (NO_Layers > 1) {
            theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_2nd_Ply_r_Ep.resize(6,6);   Dmat_2nd_Ply_r_Ep.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_2nd_Ply_r_Ep); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 2) {
            theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_3rd_Ply_r_Ep.resize(6,6);   Dmat_3rd_Ply_r_Ep.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_3rd_Ply_r_Ep); CHKERRQ(ierr);
          }
          // Fourth-layer
          if (NO_Layers > 3) {
            theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_4th_Ply_r_Ep.resize(6,6);   Dmat_4th_Ply_r_Ep.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_4th_Ply_r_Ep); CHKERRQ(ierr);
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_r_Ep_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Ep,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Ep_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Ep,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Ep_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Ep,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Ep_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Ep,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ep(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Ep_1st_Ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Ep_2nd_Ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Ep_3rd_Ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Ep_4th_Ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {// due to axial modulus of fibre (Ez)
          idx_disp = 1;cout<<"\nDmat_r_Ez: "<<Dmat_r_Ez<<endl;
          // Examined ply
          ierr = Dmat_Transformation(Angle_REL, Dmat_r_Ez, Dmat_REL_r); CHKERRQ(ierr);
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_Ez.resize(6,6);   Dmat_1st_Ply_r_Ez.clear();
          ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_1st_Ply_r_Ez); CHKERRQ(ierr);
          // Second-layer
          if (NO_Layers > 1) {
            theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_2nd_Ply_r_Ez.resize(6,6);   Dmat_2nd_Ply_r_Ez.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_2nd_Ply_r_Ez); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 2) {
            theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_3rd_Ply_r_Ez.resize(6,6);   Dmat_3rd_Ply_r_Ez.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_3rd_Ply_r_Ez); CHKERRQ(ierr);
          }
          // Fourth-layer
          if (NO_Layers > 3) {
            theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_4th_Ply_r_Ez.resize(6,6);   Dmat_4th_Ply_r_Ez.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_4th_Ply_r_Ez); CHKERRQ(ierr);
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_r_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Ez,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Ez,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Ez,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Ez,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_1st_Ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_2nd_Ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_3rd_Ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_4th_Ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r); CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {// due to shear modulus of fibre (Gzp)
          idx_disp = 1;
          // Examined ply
          ierr = Dmat_Transformation(Angle_REL, Dmat_r_Gzp, Dmat_REL_r); CHKERRQ(ierr);
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_Gzp.resize(6,6);   Dmat_1st_Ply_r_Gzp.clear();
          ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_1st_Ply_r_Gzp); CHKERRQ(ierr);
          // Second-layer
          if (NO_Layers > 1) {
            theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_2nd_Ply_r_Gzp.resize(6,6);   Dmat_2nd_Ply_r_Gzp.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_2nd_Ply_r_Gzp); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 2) {
            theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_3rd_Ply_r_Gzp.resize(6,6);   Dmat_3rd_Ply_r_Gzp.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_3rd_Ply_r_Gzp); CHKERRQ(ierr);
          }
          // Fourth-layer
          if (NO_Layers > 3) {
            theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_4th_Ply_r_Gzp.resize(6,6);   Dmat_4th_Ply_r_Gzp.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_4th_Ply_r_Gzp); CHKERRQ(ierr);
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Gzp,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Gzp,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Gzp,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Gzp,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Gzp(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Gzp_1st_Ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Gzp_2nd_Ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Gzp_3rd_Ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Gzp_4th_Ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,2,"Ef") == 0) {// due to Young's modulus of fibre (Ef) - isotropic
          idx_disp = 1;
          // Examined ply
          ierr = Dmat_Transformation(Angle_REL, Dmat_r_Ef, Dmat_REL_r); CHKERRQ(ierr);
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_Ef.resize(6,6);   Dmat_1st_Ply_r_Ef.clear();
          ierr = Dmat_Transformation(theta, Dmat_r_Ef, Dmat_1st_Ply_r_Ef); CHKERRQ(ierr);
          // Second-layer
          if (NO_Layers > 1) {
            theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_2nd_Ply_r_Ef.resize(6,6);   Dmat_2nd_Ply_r_Ef.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ef, Dmat_2nd_Ply_r_Ef); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 2) {
            theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_3rd_Ply_r_Ef.resize(6,6);   Dmat_3rd_Ply_r_Ef.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ef, Dmat_3rd_Ply_r_Ef); CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_4th_Ply_r_Ef.resize(6,6);   Dmat_4th_Ply_r_Ef.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ef, Dmat_4th_Ply_r_Ef); CHKERRQ(ierr);
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_r_Ef_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Ef,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Ef_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Ef,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Ef_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Ef,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Ef_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Ef,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ef(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ef); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Ef_1st_Ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Ef_2nd_Ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Ef_3rd_Ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Ef_4th_Ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ef); CHKERRQ(ierr);
        }
        else if (VariableName.compare(0,3,"NUf") == 0) { // due to Poisson's ratio of fibre (NUf) - isotropic
          idx_disp = 1;
          // Examined ply
          ierr = Dmat_Transformation(Angle_REL, Dmat_r_NUf, Dmat_REL_r); CHKERRQ(ierr);
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_NUf.resize(6,6);   Dmat_1st_Ply_r_NUf.clear();
          ierr = Dmat_Transformation(theta, Dmat_r_NUf, Dmat_1st_Ply_r_NUf); CHKERRQ(ierr);
          // Second-layer
          if (NO_Layers > 1) {
            theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_2nd_Ply_r_NUf.resize(6,6);   Dmat_2nd_Ply_r_NUf.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUf, Dmat_2nd_Ply_r_NUf); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 2) {
            theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_3rd_Ply_r_NUf.resize(6,6);   Dmat_3rd_Ply_r_NUf.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUf, Dmat_3rd_Ply_r_NUf); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 3) {
            theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_4th_Ply_r_NUf.resize(6,6);   Dmat_4th_Ply_r_NUf.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUf, Dmat_4th_Ply_r_NUf); CHKERRQ(ierr);
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUf_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_NUf,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUf_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_NUf,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUf_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_NUf,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_NUf_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_NUf,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUf(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUf); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUf_1st_Ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_NUf_2nd_Ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_NUf_3rd_Ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_NUf_4th_Ply);  CHKERRQ(ierr);
          }
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUf); CHKERRQ(ierr);
        }
        /***********************************************************************
         *
         * 3.2. Case 2: Applied forces are treated as random variables
         *
         **********************************************************************/
        else if (VariableName.compare(0,5,"force") == 0) {
          idx_disp = 1;
          
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = MatZeroEntries(A); CHKERRQ(ierr);
          
          // Establish an object of elastic FE method
          MyElasticFEMethod_Macro my_fe_1st_Ply_r_F (m_field_Macro,A,D,F,Dmat_1st_Ply,"DISP_MACRO");
          MyElasticFEMethod_Macro my_fe_2nd_Ply_r_F (m_field_Macro,A,D,F,Dmat_2nd_Ply,"DISP_MACRO");
          MyElasticFEMethod_Macro my_fe_3rd_Ply_r_F (m_field_Macro,A,D,F,Dmat_3rd_Ply,"DISP_MACRO");
          MyElasticFEMethod_Macro my_fe_4th_Ply_r_F (m_field_Macro,A,D,F,Dmat_4th_Ply,"DISP_MACRO");
          
          //preproc
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
          
          // Calculate applied forces and preassures on surface and
          // assemble force vector
          boost::ptr_map<string,NeummanForcesSurface> my_neumann_forces;
          MyMetaNeummanForces_r_PSFEM First_FE;
          ierr = First_FE.addNeumannBCElements(m_field_Macro,"DISP_MACRO"); CHKERRQ(ierr);
          ierr = First_FE.setNeumannFiniteElementOperators(m_field_Macro,my_neumann_forces,dF,"DISP_MACRO"); CHKERRQ(ierr);
          boost::ptr_map<string,NeummanForcesSurface>::iterator mitt = my_neumann_forces.begin();
          for(;mitt!=my_neumann_forces.end();mitt++) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO",mitt->first,mitt->second->getLoopFe()); CHKERRQ(ierr);
          }
          
          // Assemble stiffness matrix
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe_1st_Ply_r_F);     CHKERRQ(ierr);
          // Third layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe_2nd_Ply_r_F);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe_3rd_Ply_r_F);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe_4th_Ply_r_F);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
          
          ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
          
        }
        
        /***********************************************************************
         *
         * 3.3. Case 3: Ply orientation angle is treated as random variable
         *
         **********************************************************************/
        else if (VariableName.compare(0,11,"orientation") == 0) {
          idx_disp = 1;
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_Theta.resize(6,6);   Dmat_1st_Ply_r_Theta.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_1st_Ply_r_Theta); CHKERRQ(ierr);
          // Second-layer
          if (NO_Layers > 1) {
            theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_2nd_Ply_r_Theta.resize(6,6);   Dmat_2nd_Ply_r_Theta.clear();
            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_2nd_Ply_r_Theta); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 2) {
            theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_3rd_Ply_r_Theta.resize(6,6);   Dmat_3rd_Ply_r_Theta.clear();
            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_3rd_Ply_r_Theta); CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_4th_Ply_r_Theta.resize(6,6);   Dmat_4th_Ply_r_Theta.clear();
            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_4th_Ply_r_Theta); CHKERRQ(ierr);
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Theta_1st_Ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Theta_2nd_Ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Theta_3rd_Ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Theta_4th_Ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta); CHKERRQ(ierr);
          
        }
        else if (VariableName.compare(0,6,"theta1") == 0) { // first layer angle is a random variable
          idx_disp = 1;
          // Examined ply
          if (ExaminedPly == 1) {
            Dmat_REL_r = Dmat_1st_Ply_r_Theta_1;
          }
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_Theta_1.resize(6,6);   Dmat_1st_Ply_r_Theta_1.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_1st_Ply_r_Theta_1); CHKERRQ(ierr);
          //cout<<"\n\nThe first ply angle is "<<PlyAngle(0)<<endl;
          //cout<<"\n\n Dmat_1st_Ply_r_Theta_1 is "<<Dmat_1st_Ply_r_Theta_1<<endl;
          // Second-layer
          if (NO_Layers > 1) {
            Dmat_2nd_Ply_r_Theta_1.resize(6,6);   Dmat_2nd_Ply_r_Theta_1.clear();
          }
          // Third-layer
          if (NO_Layers > 2) {
            Dmat_3rd_Ply_r_Theta_1.resize(6,6);   Dmat_3rd_Ply_r_Theta_1.clear();
          }
          // Fourth layer
          if (NO_Layers > 3) {
            Dmat_4th_Ply_r_Theta_1.resize(6,6);   Dmat_4th_Ply_r_Theta_1.clear();
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta_1,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta_1,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta_1,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta_1,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_1(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_1); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_1);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_1);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_1);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_1);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_1); CHKERRQ(ierr);
          
          cout<<"Solving the first order derivative of Theta_1"<<endl;
        }
        else if (VariableName.compare(0,6,"theta2") == 0) { // second layer angle is a random variable
          idx_disp = 1;
          // Examined ply
          if (ExaminedPly == 2) {
            Dmat_REL_r = Dmat_2nd_Ply_r_Theta_2;
          }
          // First-layer
          Dmat_1st_Ply_r_Theta_2.resize(6,6);   Dmat_1st_Ply_r_Theta_2.clear();
          // Second-layer
          theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_2nd_Ply_r_Theta_2.resize(6,6);   Dmat_2nd_Ply_r_Theta_2.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_2nd_Ply_r_Theta_2); CHKERRQ(ierr);
          // Third-layer
          if (NO_Layers > 2) {
            Dmat_3rd_Ply_r_Theta_2.resize(6,6);   Dmat_3rd_Ply_r_Theta_2.clear();
          }
          // Fourth layer
          if (NO_Layers > 3) {
            Dmat_4th_Ply_r_Theta_2.resize(6,6);   Dmat_4th_Ply_r_Theta_2.clear();
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta_2,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta_2,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta_2,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta_2,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_2(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_2); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_2);  CHKERRQ(ierr);
          // Second layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_2);  CHKERRQ(ierr);
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_2);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_2);  CHKERRQ(ierr);
          }
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_2); CHKERRQ(ierr);
    
          cout<<"Solving the first order derivative of Theta_2"<<endl;
        }
        else if (VariableName.compare(0,6,"theta3") == 0) { // third layer angle is a random variable
          idx_disp = 1;
          // Examined ply
          if (ExaminedPly == 3) {
            Dmat_REL_r = Dmat_3rd_Ply_r_Theta_3;
          }
          // First-layer
          Dmat_1st_Ply_r_Theta_3.resize(6,6);   Dmat_1st_Ply_r_Theta_3.clear();
          // Second-layer
          Dmat_2nd_Ply_r_Theta_3.resize(6,6);   Dmat_2nd_Ply_r_Theta_3.clear();
          // Third-layer
          theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_3rd_Ply_r_Theta_3.resize(6,6);   Dmat_3rd_Ply_r_Theta_3.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_3rd_Ply_r_Theta_3); CHKERRQ(ierr);
          // Fourth layer
          if (NO_Layers > 3) {
            Dmat_4th_Ply_r_Theta_3.resize(6,6);   Dmat_4th_Ply_r_Theta_3.clear();
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta_3,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta_3,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta_3,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta_3,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_3(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_3); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_3);  CHKERRQ(ierr);
          // Second layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_3);  CHKERRQ(ierr);
          // Third layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_3);  CHKERRQ(ierr);
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_3);  CHKERRQ(ierr);
          }
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_3); CHKERRQ(ierr);
          
          cout<<"Solving the first order derivative of Theta_3"<<endl;
        }
        else if (VariableName.compare(0,6,"theta4") == 0) { // fourth layer angle is a random variable
          idx_disp = 1;
          // Examined ply
          if (ExaminedPly == 3) {
            Dmat_REL_r = Dmat_4th_Ply_r_Theta_4;
          }
          // First-layer
          Dmat_1st_Ply_r_Theta_4.resize(6,6);   Dmat_1st_Ply_r_Theta_4.clear();
          // Second-layer
          Dmat_2nd_Ply_r_Theta_4.resize(6,6);   Dmat_2nd_Ply_r_Theta_4.clear();
          // Third-layer
          Dmat_3rd_Ply_r_Theta_4.resize(6,6);   Dmat_3rd_Ply_r_Theta_4.clear();
          // Fourth layer
          theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_4th_Ply_r_Theta_4.resize(6,6);   Dmat_4th_Ply_r_Theta_4.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_4th_Ply_r_Theta_4); CHKERRQ(ierr);
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta_4,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta_4,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta_4,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta_4,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_4(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_4); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_4);  CHKERRQ(ierr);
          // Second layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_4);  CHKERRQ(ierr);
          // Third layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_4);  CHKERRQ(ierr);
          // Fourth layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_4);  CHKERRQ(ierr);
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_4); CHKERRQ(ierr);
          
          cout<<"Solving the first order derivative of Theta_4"<<endl;
        }
        
        // post-processing
        if (idx_disp==1) {
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);//ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          
          // Calculate stress for examined element
          
          FE2_PostProcStressForReliability_First Calc_Stress_r(m_field_Macro,"DISP_MACRO","DISP_MACRO_r",Dmat_REL,Dmat_REL_r);
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_r);  CHKERRQ(ierr);
          StressGP_r.clear(); REL_Stress_Transformation(Angle_REL, Calc_Stress_r.StressGP_r, StressGP_r);
          TheStress(ii).resize(3,3); TheStress(ii).clear(); TheStress(ii) = StressGP_r;
        }
        
        VariableName.clear();
        
        ierr = VecDestroy(&dF); CHKERRQ(ierr);
        ierr = VecDestroy(&dD); CHKERRQ(ierr);
      }
      
      /***************************************************************************
       *
       *  4. FINISH
       *
       **************************************************************************/
      
      //Destroy matrices
      ierr = VecDestroy(&F); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = MatDestroy(&A); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver_Macro); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    //~FE2_Macro_Solver_Laminate() {};
  };
  
  //typedef struct _FE2_Macro_Solver_Laminate *FE2_Macro_Solver_Laminate;
}

#endif //__FE2_MACRO_SOLVER_HPP