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

#ifndef __FE2_MACRO_SOLVER_MCS_HPP
#define __FE2_MACRO_SOLVER_MCS_HPP

namespace MoFEM {
  struct FE2_Macro_Solver_MCS {
    ublas::matrix<FieldData> Dmat;
    ublas::matrix<double> Dmat_1st_Ply;
    ublas::matrix<double> Dmat_2nd_Ply;
    ublas::matrix<double> Dmat_3rd_Ply;
    ublas::matrix<double> Dmat_4th_Ply;
    
    //------------------------------------------------------------------------------
    // To transform constitutive matrix
    //
    
    virtual PetscErrorCode Dmat_Transformation(double theta, ublas::matrix<FieldData> Dmat_123, ublas::matrix<FieldData> &Dmat_xyz) {
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
    
    virtual PetscErrorCode Dmat_Transformation_r_Theta(double theta, ublas::matrix<FieldData> Dmat_123, ublas::matrix<FieldData> &Dmat_xyz) {
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
    //
    // =========================================================================
    virtual PetscErrorCode Calculate_RVEDmat (FieldInterface &m_field_RVE,
                                              ublas::vector<double> vars_value,
                                              int num_rvars,
                                              vector<string> vars_name) {
      PetscFunctionBegin;
      ErrorCode rval;
      PetscErrorCode ierr;
      
      Dmat.resize(6,6); Dmat.clear();
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      Vec F1,F2,F3,F4,F5,F6,D1,D2,D3,D4,D5,D6;
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D6); CHKERRQ(ierr);
      
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
      
      
      double YoungModulus = 1;
      double PoissonRatio = 0.0;
      double alpha;
      int field_rank=3;
      
      /*****************************************************************************
       *
       *  2. Get the volume of RVE
       *
       ****************************************************************************/
      double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
      Vec RVE_volume_Vec;
      ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field_RVE.get_moab(),MYPCOMM_INDEX);
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
      ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
      
      RVEVolume MyRVEVol(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
      RVEVolumeTrans MyRVEVolTrans(m_field_RVE,Aij,D1,F1, RVE_volume_Vec);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVolTrans);  CHKERRQ(ierr);
      //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
      cout<<"Final RVE_volume = "<< RVE_volume <<endl;
      
      
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)) {
        cout << endl << *it << endl;
        //
        // Get block name
        //
        string name = it->get_name();
        //
        // Modify material properties of matrix
        //
        if (name.compare(0,13,"MAT_ELASTIC_1") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_rvars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Em") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = vars_value(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUm") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = vars_value(ii-1);
            }
            ParameterName.clear();
          }
          
          YoungModulus=mydata.data.Young;
          PoissonRatio=mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Matrix material:\n" << mydata;
        }
        //
        // Update material properties of inclusion/fibre reinforcement
        //
        if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
          Mat_Elastic_TransIso mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int i=1;i<=num_rvars;i++) {
            ParameterName = vars_name[i];
            cout<<"the variable name is "<<vars_name[i]<<endl;
            if (ParameterName.compare(0,2,"Ez") == 0) {
              mydata.data.Youngz = vars_value(i-1);
            }
            else if (ParameterName.compare(0,2,"Ep") == 0) {
              mydata.data.Youngp = vars_value(i-1);
            }
            else if (ParameterName.compare(0,3,"NUp") == 0) {
              mydata.data.Poissonp = vars_value(i-1);
            }
            else if (ParameterName.compare(0,4,"NUpz") == 0) {
              mydata.data.Poissonpz = vars_value(i-1);
            }
            else if (ParameterName.compare(0,3,"Gzp") == 0) {
              mydata.data.Shearzp = vars_value(i-1);
            }
            ParameterName.clear();
          }
          
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre material:\n" << mydata;
        }
      }
      
      cout<<YoungModulus<<"\t"<<PoissonRatio<<endl;
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
      applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
      
      MyElasticFEMethod my_fe(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),"DISP_RVE");
      TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(m_field_RVE,Aij,D1,F1,"DISP_RVE");
      ElasticFE_RVELagrange_Disp_Multi_Rhs MyFE_RVELagrange(m_field_RVE,Aij,D1,F1,F2,F3,F4,F5,F6,applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank);
      
      cout<<"After ElasticFE_RVELagrange_Disp_Multi_Rhs = "<<endl;
      
      ierr = VecZeroEntries(F1); CHKERRQ(ierr);
      ierr = VecZeroEntries(F2); CHKERRQ(ierr);
      ierr = VecZeroEntries(F3); CHKERRQ(ierr);
      ierr = VecZeroEntries(F4); CHKERRQ(ierr);
      ierr = VecZeroEntries(F5); CHKERRQ(ierr);
      ierr = VecZeroEntries(F6); CHKERRQ(ierr);
      ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
      ierr = VecZeroEntries(D1); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
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
      
      
      /*****************************************************************************
       *
       *  3. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ****************************************************************************/
      //Solver
      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      ierr = KSPSetUp(solver); CHKERRQ(ierr);
      
      //create a vector for 6 components of homogenized stress
      Vec Stress_Homo;
      PetscScalar *avec;
      if(pcomm->rank()==0) {
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo);
      } else {
        int ghost[] = {0,1,2,3,4,5};
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo);
        
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 1: applied macro strain: [1 0 0 0 0 0]^T
      //------------------------------------------------------------------------
      
      cout<<"===============================================================\n";
      cout<<"        Applied strain [1 0 0 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //solve for F1 and D1
      ierr = KSPSolve(solver,F1,D1); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // calculate homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_1(m_field_RVE,Aij,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
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
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 2: applied macro strain: [0 1 0 0 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 1 0 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      // solve for F2 and D2
      ierr = KSPSolve(solver,F2,D2); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_2(m_field_RVE,Aij,D2,F2,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
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
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 3: applied macro strain: [0 0 1 0 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 1 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //solve for F3 and D3
      ierr = KSPSolve(solver,F3,D3); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_3(m_field_RVE,Aij,D3,F3,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
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
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 4: applied macro strain: [0 0 0 1 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 0 1 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //solve for F4 and D4
      ierr = KSPSolve(solver,F4,D4); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      
      // Extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_4(m_field_RVE,Aij,D4,F4,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
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
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 5: applied macro strain: [0 0 0 0 1 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 0 0 1 0]^T\n";
      cout<<"===============================================================\n";
      
      //solve for F5 and D5
      ierr = KSPSolve(solver,F5,D5); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_5(m_field_RVE,Aij,D5,F5,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
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
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 6: applied macro strain: [0 0 0 0 0 1]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"         Applied strain [0 0 0 0 0 1]^T\n";
      cout<<"===============================================================\n";
      
      //solve for F6 and D6
      ierr = KSPSolve(solver,F6,D6); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_6(m_field_RVE,Aij,D6,F6,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
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
      
      ierr = VecDestroy(&D1); CHKERRQ(ierr);
      ierr = VecDestroy(&D2); CHKERRQ(ierr);
      ierr = VecDestroy(&D3); CHKERRQ(ierr);
      ierr = VecDestroy(&D4); CHKERRQ(ierr);
      ierr = VecDestroy(&D5); CHKERRQ(ierr);
      ierr = VecDestroy(&D6); CHKERRQ(ierr);
      
      ierr = MatDestroy(&Aij); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    
    // =========================================================================
    //
    //  A.VI. SOLUTION PHASE:
    //        Caculate RVE constitutive matrix Dmat
    //
    // =========================================================================
    virtual PetscErrorCode Calculate_RVEDmat_Periodic(FieldInterface &m_field_RVE,
                                                      ublas::vector<double> vars_value,
                                                      int num_rvars,
                                                      vector<string> vars_name) {
      PetscFunctionBegin;
      ErrorCode rval;
      PetscErrorCode ierr;
      
      Dmat.resize(6,6); Dmat.clear();
      
      /*****************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ****************************************************************************/
      
      vector<Vec> F(6);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F[0]); CHKERRQ(ierr);
      for(int ii = 1;ii<6;ii++) {
        ierr = VecDuplicate(F[0],&F[ii]); CHKERRQ(ierr);
      }
      
      Vec D;
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
      double alpha;
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
      
      RVEVolume MyRVEVol(m_field_RVE,Aij,D,F[0],LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
      RVEVolumeTrans MyRVEVolTrans(m_field_RVE,Aij,D,F[0], RVE_volume_Vec);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVolTrans);  CHKERRQ(ierr);
      //ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
      cout<<"Final RVE_volume = "<< RVE_volume <<endl;
      
      
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)) {
        cout << endl << *it << endl;
        //
        // Get block name
        //
        string name = it->get_name();
        //
        // Modify material properties of matrix
        //
        if (name.compare(0,13,"MAT_ELASTIC_1") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_rvars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Em") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = vars_value(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUm") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = vars_value(ii-1);
            }
            ParameterName.clear();
          }
          
          YoungModulus=mydata.data.Young;
          PoissonRatio=mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Matrix material:\n" << mydata;
        }
        //
        // Update material properties of inclusion/fibre reinforcement
        //
        if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
          Mat_Elastic_TransIso mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int i=1;i<=num_rvars;i++) {
            ParameterName = vars_name[i];
            cout<<"the variable name is "<<vars_name[i]<<endl;
            if (ParameterName.compare(0,2,"Ez") == 0) {
              mydata.data.Youngz = vars_value(i-1);
            }
            else if (ParameterName.compare(0,2,"Ep") == 0) {
              mydata.data.Youngp = vars_value(i-1);
            }
            else if (ParameterName.compare(0,3,"NUp") == 0) {
              mydata.data.Poissonp = vars_value(i-1);
            }
            else if (ParameterName.compare(0,4,"NUpz") == 0) {
              mydata.data.Poissonpz = vars_value(i-1);
            }
            else if (ParameterName.compare(0,3,"Gzp") == 0) {
              mydata.data.Shearzp = vars_value(i-1);
            }
            ParameterName.clear();
          }
          
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre material:\n" << mydata;
        }
      }
      
      cout<<YoungModulus<<"\t"<<PoissonRatio<<endl;
      
      ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
      applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();

      MyElasticFEMethod MyFE(m_field_RVE,Aij,D,F[0],LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),"DISP_RVE");
      TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(m_field_RVE,Aij,D,F[0],"DISP_RVE");
      
      
      ElasticFE_RVELagrange_Periodic_Multi_Rhs MyFE_RVELagrangePeriodic(m_field_RVE,Aij,D,F[0],F[1],F[2],F[3],F[4],F[5],applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ElasticFE_RVELagrange_RigidBodyTranslation MyFE_RVELagrangeRigidBodyTrans(m_field_RVE,Aij,D,F[0],applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank,"Lagrange_mul_disp_rigid_trans");
      
      for(int ii = 0; ii<6; ii++) {
        ierr = VecZeroEntries(F[ii]); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(F[ii],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(F[ii],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      }
      ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyFE);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyTIsotFE);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVELagrangePeriodic);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE_rigid_trans",MyFE_RVELagrangeRigidBodyTrans);  CHKERRQ(ierr);
      
      for(int ii = 0; ii<6; ii++) {
        ierr = VecGhostUpdateBegin(F[ii],ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(F[ii],ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(F[ii]); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F[ii]); CHKERRQ(ierr);
      }
      
      ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      
      
      /*****************************************************************************
       *
       *  3. SOLVE THE ZEROTH-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ****************************************************************************/
      
      //----------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //----------------------------------------------------------------------------
      //Solver
      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      ierr = KSPSetUp(solver); CHKERRQ(ierr);
      
      //create a vector for 6 components of homogenized stress
      Vec Stress_Homo;
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
      
      for(int ics = 0; ics<6; ics++) {
        
        cout<<"===============================================================\n";
        switch (ics) {
          case 0:
            cout<<"        Applied strain [1 0 0 0 0 0]^T\n"; break;
          case 1:
            cout<<"        Applied strain [0 1 0 0 0 0]^T\n"; break;
          case 2:
            cout<<"        Applied strain [0 0 1 0 0 0]^T\n"; break;
          case 3:
            cout<<"        Applied strain [0 0 0 1 0 0]^T\n"; break;
          case 4:
            cout<<"        Applied strain [0 0 0 0 1 0]^T\n"; break;
          case 5:
            cout<<"        Applied strain [0 0 0 0 0 1]^T\n"; break;
        }
        cout<<"===============================================================\n";
        
        ierr = VecZeroEntries(D); CHKERRQ(ierr);
        ierr = KSPSolve(solver,F[ics],D); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        
        //----------------------------------------------------------------------------
        // 3.2 Calculating zeroth-order homogenized stress using volume averaging theorem
        //----------------------------------------------------------------------------
        
        ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
        ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic(m_field_RVE,Aij,D,F[ics],&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressPeriodic);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec;
          VecGetArray(Stress_Homo, &avec);
          //cout<< "\nStress_Homo =\n";
          for(int kk=0; kk<6; kk++) {
            Dmat(kk,ics) = *avec;
            //cout.precision(15);
            //cout <<*avec<<endl;
            avec++;
          }
          VecRestoreArray(Stress_Homo,&avec);
        }
      }
      
      /*************************************************************************
       *
       *  4. FINISH
       *
       ************************************************************************/
      
      //Destroy matrices and vectors
      for(int ii = 0;ii<6;ii++) {
        ierr = VecDestroy(&F[ii]); CHKERRQ(ierr);
      }
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = MatDestroy(&Aij); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      ierr = VecDestroy(&RVE_volume_Vec); CHKERRQ(ierr);
      ierr = VecDestroy(&Stress_Homo); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    
    // =========================================================================
    //
    //  B.VI. SOLUTION PHASE:
    //        Solve Macroscale FE equation
    //
    // =========================================================================
    virtual PetscErrorCode Macro_FE_Solver(FieldInterface &m_field_Macro,
                                           ublas::vector<double> vars_value,
                                           int num_rvars,
                                           vector<string> vars_name) {
      PetscFunctionBegin;
      
      ErrorCode rval;
      PetscErrorCode ierr;
      
      /*****************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ****************************************************************************/
      //create matrices
      Vec F,D;
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
      MyElasticFEMethod_Macro my_fe_Macro(m_field_Macro,A,D,F,Dmat,"DISP_MACRO");
      
      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = MatZeroEntries(A); CHKERRQ(ierr);
      
      ierr = m_field_Macro.set_global_VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      //preproc
      ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
      //loop elems
      //PetscBarrier(PETSC_NULL);
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe_Macro);  CHKERRQ(ierr);
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe_Macro);  CHKERRQ(ierr);
      
      //forces and preassures on surface
      boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
      // Check whether force is considered as random variable or not
      int idx_force = 100;
      for (int ii = 1; ii<=num_rvars; ii++) {
        string VariableName;
        VariableName = vars_name[ii];
        if (VariableName.compare(0,5,"force") == 0) {
          idx_force = ii;cout<<"The force is: "<<vars_value(ii-1)<<endl;
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
                                                          vars_value(idx_force-1)); CHKERRQ(ierr);
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
      
      
      // ===========================================================================
      //
      //  B.VIII. FINISH
      //
      // ===========================================================================
      //Destroy matrices
      ierr = VecDestroy(&F); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = MatDestroy(&A); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver_Macro); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    // =========================================================================
    //
    //  B.VI. SOLUTION PHASE:
    //        Solve Macroscale FE equation
    //
    // =========================================================================
    
    virtual PetscErrorCode Macro_FE_Laminate(FieldInterface &m_field_Macro,
                                        ublas::vector<double> TheVariables,
                                        int num_rvars,
                                        vector<string> vars_name,
                                        ublas::vector<double> PlyAngle,
                                        PetscInt NO_Layers) {
      PetscFunctionBegin;
      
      ErrorCode rval;
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
      Dmat_1st_Ply.resize(6,6);   Dmat_1st_Ply.clear();
      ierr = Dmat_Transformation(theta, Dmat, Dmat_1st_Ply); CHKERRQ(ierr);
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
      
      MyElasticFEMethod_Macro my_fe_1st_ply(m_field_Macro,A,D,F,Dmat_1st_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_2nd_ply(m_field_Macro,A,D,F,Dmat_2nd_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_3rd_ply(m_field_Macro,A,D,F,Dmat_3rd_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_4th_ply(m_field_Macro,A,D,F,Dmat_4th_Ply,"DISP_MACRO");
      
      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = MatZeroEntries(A); CHKERRQ(ierr);
      
      ierr = m_field_Macro.set_global_VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
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
    
  };
  
}

#endif //__FE2_MACRO_SOLVER_MCS_HPP