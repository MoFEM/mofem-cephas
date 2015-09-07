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
  struct FE2_Macro_Solver {
    
    //global variable Dmat
    ublas::matrix<double> Dmat;
    ublas::matrix<double> Dmat_r_Em;
    ublas::matrix<double> Dmat_r_NUm;
    ublas::matrix<double> Dmat_r_Ep;
    ublas::matrix<double> Dmat_r_Ez;
    ublas::matrix<double> Dmat_r_NUp;
    ublas::matrix<double> Dmat_r_NUpz;
    ublas::matrix<double> Dmat_r_Gzp;
    
    ublas::matrix<double> Dmat_rs_EmEm;
    ublas::matrix<double> Dmat_rs_NUmNUm;
    ublas::matrix<double> Dmat_rs_EpEp;
    ublas::matrix<double> Dmat_rs_EzEz;
    ublas::matrix<double> Dmat_rs_NUpNUp;
    ublas::matrix<double> Dmat_rs_NUpzNUpz;
    ublas::matrix<double> Dmat_rs_GzpGzp;
    
    
    // =========================================================================
    //
    //  A.VI. SOLUTION PHASE:
    //        Caculate RVE constitutive matrix Dmat
    //
    // =========================================================================
    
    virtual PetscErrorCode Calculate_RVEDmat(FieldInterface &m_field_RVE,
                                     int &nvars, int &nders,
                                     vector<string> &stochastic_fields,
                                     ublas::vector<double> matprop,
                                     int num_rvars,
                                     vector<string> vars_name) {
      PetscFunctionBegin;
      cout <<"Hi from Calculate_RVEDmat"<<endl;
      
      ErrorCode rval;
      PetscErrorCode ierr;
      
      Dmat.resize(6,6); Dmat.clear();
      
      Dmat_r_Em.resize(6,6);   Dmat_r_Em.clear();
      Dmat_r_NUm.resize(6,6);   Dmat_r_NUm.clear();
      Dmat_r_Ep.resize(6,6);   Dmat_r_Ep.clear();
      Dmat_r_Ez.resize(6,6);   Dmat_r_Ez.clear();
      Dmat_r_NUp.resize(6,6);  Dmat_r_NUp.clear();
      Dmat_r_NUpz.resize(6,6); Dmat_r_NUpz.clear();
      Dmat_r_Gzp.resize(6,6);  Dmat_r_Gzp.clear();
      
      Dmat_rs_EmEm.resize(6,6);     Dmat_rs_EmEm.clear();
      Dmat_rs_NUmNUm.resize(6,6);     Dmat_rs_NUmNUm.clear();
      Dmat_rs_EpEp.resize(6,6);     Dmat_rs_EpEp.clear();
      Dmat_rs_EzEz.resize(6,6);     Dmat_rs_EzEz.clear();
      Dmat_rs_NUpNUp.resize(6,6);   Dmat_rs_NUpNUp.clear();
      Dmat_rs_NUpzNUpz.resize(6,6); Dmat_rs_NUpzNUpz.clear();
      Dmat_rs_GzpGzp.resize(6,6);   Dmat_rs_GzpGzp.clear();
      
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      Vec F1,F2,F3,F4,F5,F6,D1,D2,D3,D4,D5,D6;
      Vec dF1,dF2,dF3,dF4,dF5,dF6;
      Vec ddF1,ddF2,ddF3,ddF4,ddF5,ddF6;
      Vec dD1,dD2,dD3,dD4,dD5,dD6;
      Vec ddD1,ddD2,ddD3,ddD4,ddD5,ddD6;
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD6); CHKERRQ(ierr);
      
      
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
      
      RVEVolume MyRVEVol(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
      RVEVolumeTrans MyRVEVolTrans(m_field_RVE,Aij,D1,F1, RVE_volume_Vec);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVolTrans);  CHKERRQ(ierr);
      //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
      cout<<"Final RVE_volume = "<< RVE_volume <<endl;
      
      
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)) {
        cout << endl << *it << endl;
        
        //Get block name
        string name = it->get_name();
        
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
            else if (ParameterName.compare(0,4,"NUpz") == 0) {
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"Gzp") == 0) {
              mydata.data.Shearzp = matprop(i-1);
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
      ierr = VecZeroEntries(D2); CHKERRQ(ierr);
      ierr = VecZeroEntries(D3); CHKERRQ(ierr);
      ierr = VecZeroEntries(D4); CHKERRQ(ierr);
      ierr = VecZeroEntries(D5); CHKERRQ(ierr);
      ierr = VecZeroEntries(D6); CHKERRQ(ierr);
      
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
      Vec Stress_Homo, Stress_Homo_r, Stress_Homo_rs;
      PetscScalar *avec;
      
      if(pcomm->rank()==0) {
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo_r);
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo_rs);
      } else {
        int ghost[] = {0,1,2,3,4,5};
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo_r);
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo_rs);
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
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      
      for(int ii=1; ii <= num_rvars; ii++) {
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }

        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D1,dF1,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D1,dF1,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D1,dF1,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D1,dF1,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D1,dF1,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
                 
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF1); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF1); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF1,dD1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF1); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF1); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF1,ddD1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r1(m_field_RVE,Aij,dD1,dF1,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r1);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++){
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,0) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,0) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,0) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,0) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,0) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,0) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,0) = *avec_r; break;
                case 7:
                  Dmat_rs_EmEm(irow,0) = *avec_r; break;
                case 8:
                  Dmat_rs_NUmNUm(irow,0) = *avec_r; break;
                case 9:
                  Dmat_rs_NUpNUp(irow,0) = *avec_r; break;
                case 10:
                  Dmat_rs_NUpzNUpz(irow,0) = *avec_r; break;
                case 11:
                  Dmat_rs_EpEp(irow,0) = *avec_r; break;
                case 12:
                  Dmat_rs_EzEz(irow,0) = *avec_r; break;
                case 13:
                  Dmat_rs_GzpGzp(irow,0) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        cout<< "\n\n";
        
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
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=1; ii <= num_rvars; ii++) {
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }

        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D2,dF2,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D2,dF2,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D2,dF2,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D2,dF2,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D2,dF2,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF2); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF2); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF2,dD2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF2); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF2); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF2,ddD2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r2(m_field_RVE,Aij,dD2,dF2,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r2);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,1) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,1) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,1) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,1) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,1) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,1) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,1) = *avec_r; break;
                case 7:
                  Dmat_rs_EmEm(irow,1) = *avec_r; break;
                case 8:
                  Dmat_rs_NUmNUm(irow,1) = *avec_r; break;
                case 9:
                  Dmat_rs_NUpNUp(irow,1) = *avec_r; break;
                case 10:
                  Dmat_rs_NUpzNUpz(irow,1) = *avec_r; break;
                case 11:
                  Dmat_rs_EpEp(irow,1) = *avec_r; break;
                case 12:
                  Dmat_rs_EzEz(irow,1) = *avec_r; break;
                case 13:
                  Dmat_rs_GzpGzp(irow,1) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        
        cout<< "\n\n";
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
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=1; ii <= num_rvars; ii++) {
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        VariableName.clear();
        
          
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D3,dF3,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D3,dF3,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D3,dF3,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D3,dF3,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D3,dF3,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF3); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF3); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF3,dD3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF3); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF3); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF3,ddD3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r3(m_field_RVE,Aij,dD3,dF3,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r3);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,2) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,2) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,2) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,2) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,2) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,2) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,2) = *avec_r; break;
                case 7:
                  Dmat_rs_EmEm(irow,2) = *avec_r; break;
                case 8:
                  Dmat_rs_NUmNUm(irow,2) = *avec_r; break;
                case 9:
                  Dmat_rs_NUpNUp(irow,2) = *avec_r; break;
                case 10:
                  Dmat_rs_NUpzNUpz(irow,2) = *avec_r; break;
                case 11:
                  Dmat_rs_EpEp(irow,2) = *avec_r; break;
                case 12:
                  Dmat_rs_EzEz(irow,2) = *avec_r; break;
                case 13:
                  Dmat_rs_GzpGzp(irow,2) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        
        cout<< "\n\n";
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
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=1; ii <= num_rvars; ii++) {
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D4,dF4,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D4,dF4,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D4,dF4,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D4,dF4,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D4,dF4,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF4); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF4); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF4,dD4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF4); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF4); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF4,ddD4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r4(m_field_RVE,Aij,dD4,dF4,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r4);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,3) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,3) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,3) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,3) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,3) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,3) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,3) = *avec_r; break;
                case 7:
                  Dmat_rs_EmEm(irow,3) = *avec_r; break;
                case 8:
                  Dmat_rs_NUmNUm(irow,3) = *avec_r; break;
                case 9:
                  Dmat_rs_NUpNUp(irow,3) = *avec_r; break;
                case 10:
                  Dmat_rs_NUpzNUpz(irow,3) = *avec_r; break;
                case 11:
                  Dmat_rs_EpEp(irow,3) = *avec_r; break;
                case 12:
                  Dmat_rs_EzEz(irow,3) = *avec_r; break;
                case 13:
                  Dmat_rs_GzpGzp(irow,3) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        
        cout<< "\n\n";
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
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=1; ii <= num_rvars; ii++) {
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D5,dF5,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D5,dF5,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D5,dF5,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D5,dF5,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D5,dF5,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF5); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF5); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF5,dD5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF5); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF5); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF5,ddD5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r5(m_field_RVE,Aij,dD5,dF5,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r5);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,4) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,4) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,4) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,4) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,4) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,4) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,4) = *avec_r; break;
                case 7:
                  Dmat_rs_EmEm(irow,4) = *avec_r; break;
                case 8:
                  Dmat_rs_NUmNUm(irow,4) = *avec_r; break;
                case 9:
                  Dmat_rs_NUpNUp(irow,4) = *avec_r; break;
                case 10:
                  Dmat_rs_NUpzNUpz(irow,4) = *avec_r; break;
                case 11:
                  Dmat_rs_EpEp(irow,4) = *avec_r; break;
                case 12:
                  Dmat_rs_EzEz(irow,4) = *avec_r; break;
                case 13:
                  Dmat_rs_GzpGzp(irow,4) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        cout<< "\n\n";
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
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=1; ii <= num_rvars; ii++) {
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D6,dF6,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D6,dF6,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D6,dF6,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D6,dF6,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D6,dF6,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF6); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF6); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF6,dD6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF6); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF6); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF6,ddD6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r6(m_field_RVE,Aij,dD6,dF6,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r6);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,5) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,5) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,5) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,5) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,5) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,5) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,5) = *avec_r; break;
                case 7:
                  Dmat_rs_EmEm(irow,5) = *avec_r; break;
                case 8:
                  Dmat_rs_NUmNUm(irow,5) = *avec_r; break;
                case 9:
                  Dmat_rs_NUpNUp(irow,5) = *avec_r; break;
                case 10:
                  Dmat_rs_NUpzNUpz(irow,5) = *avec_r; break;
                case 11:
                  Dmat_rs_EpEp(irow,5) = *avec_r; break;
                case 12:
                  Dmat_rs_EzEz(irow,5) = *avec_r; break;
                case 13:
                  Dmat_rs_GzpGzp(irow,5) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        cout<< "\n\n";
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
      
      ierr = VecDestroy(&dF1); CHKERRQ(ierr);
      ierr = VecDestroy(&dF2); CHKERRQ(ierr);
      ierr = VecDestroy(&dF3); CHKERRQ(ierr);
      ierr = VecDestroy(&dF4); CHKERRQ(ierr);
      ierr = VecDestroy(&dF5); CHKERRQ(ierr);
      ierr = VecDestroy(&dF6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&ddF1); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF2); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF3); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF4); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF5); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&D1); CHKERRQ(ierr);
      ierr = VecDestroy(&D2); CHKERRQ(ierr);
      ierr = VecDestroy(&D3); CHKERRQ(ierr);
      ierr = VecDestroy(&D4); CHKERRQ(ierr);
      ierr = VecDestroy(&D5); CHKERRQ(ierr);
      ierr = VecDestroy(&D6); CHKERRQ(ierr);
      
      
      ierr = VecDestroy(&dD1); CHKERRQ(ierr);
      ierr = VecDestroy(&dD2); CHKERRQ(ierr);
      ierr = VecDestroy(&dD3); CHKERRQ(ierr);
      ierr = VecDestroy(&dD4); CHKERRQ(ierr);
      ierr = VecDestroy(&dD5); CHKERRQ(ierr);
      ierr = VecDestroy(&dD6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&ddD1); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD2); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD3); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD4); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD5); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD6); CHKERRQ(ierr);
      
      ierr = MatDestroy(&Aij); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    // =========================================================================
    //
    //  B.VI. SOLUTION PHASE:
    //        Solve Macroscale FE equation
    //
    // =========================================================================
    
    virtual PetscErrorCode Macro_FE_Solver(FieldInterface &m_field_Macro,
                                           int &nvars, int &nders,
                                           vector<string> &stochastic_fields,
                                           ublas::vector<double> TheVariables,
                                           int num_rvars,
                                           vector<string> vars_name) {
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
      Vec dF, dD;
      Vec ddF, ddD;
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&F); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&D); CHKERRQ(ierr);
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&dF); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&dD); CHKERRQ(ierr);
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&ddF); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&ddD); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat A;
      ierr = m_field_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_MACRO",&A); CHKERRQ(ierr);
      
      //Matrix View
      //MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
      //std::string wait;
      //std::cin >> wait;
      
      
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
      
      
      /*************************************************************************
       *
       *  3. SOLVE THE FIRST-ORDER AND THE SECOND-ORDER FE EQUILIBRIUM EQUATION
       *     1st order-[K][U_r] = -[K_r][U}
       *     2nd order-[K][U_rs] = -[K_rs][U]-2[K_r][U_s]
       *
       ************************************************************************/
      for (int ii=1; ii<=num_rvars; ii++) {
        
        /***********************************************************************
         *
         * 3.1. Case 1: Material properties are  treated as random variables
         *
         **********************************************************************/
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        else if (VariableName.compare(0,5,"force") == 0) {
          idx_disp = 80;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddD); CHKERRQ(ierr);
          ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        
        switch (idx_disp) {
          case 0: {// due to Young's modulus of matrix (Em)
            FE2_Rhs_r_PSFEM my_fe2_k_r_Em(m_field_Macro,A,dD,dF,Dmat_r_Em,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Em(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_Em);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_Em);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
            break;
          }
          case 1: { // due to Poisson's ratio of matrix (NUm)
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUm(m_field_Macro,A,D,dF,Dmat_r_NUm,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUm(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_NUm);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_NUm);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
            break;
          }
          case 2: {// due to transversal Poisson's ratio of fibre (NUp)
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUp(m_field_Macro,A,D,dF,Dmat_r_NUp,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_NUp);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_NUp);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
            break;
          }
          case 3: {// due to axial Poisson's ratio of fibre (NUpz)
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz(m_field_Macro,A,D,dF,Dmat_r_NUpz,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUpz(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_NUpz);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_NUpz);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
            break;
          }
          case 4: {// due to transversal modulus of fibre (Ep)
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ep(m_field_Macro,A,D,dF,Dmat_r_Ep,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ep(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_Ep);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_Ep);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
            break;
          }
          case 5: {// due to axial modulus of fibre (Ez)
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ez(m_field_Macro,A,dD,dF,Dmat_r_Ez,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ez(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_Ez);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_Ez);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
            break;
          }
          case 6: {// due to shear modulus of fibre (Gzp)
            FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp(m_field_Macro,A,D,dF,Dmat_r_Gzp,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Gzp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_Gzp);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_Gzp);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
            break;
          }
          case 7: {// 2nd order due to Em & Em
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EmEm(m_field_Macro,A,D,ddF,Dmat_r_Em,"DISP_MACRO",Dmat_rs_EmEm,"DISP_MACRO_r_Em");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_Em(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Em); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_EmEm);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_EmEm);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Em); CHKERRQ(ierr);
            break;
          }
          case 8: {// 2nd order due to NUm & NUm
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUmNUm(m_field_Macro,A,D,ddF,Dmat_r_NUm,"DISP_MACRO",Dmat_rs_NUmNUm,"DISP_MACRO_r_NUm");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_NUm(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUm); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_NUmNUm);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_NUmNUm);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUm); CHKERRQ(ierr);
            break;
          }
          case 9: {// 2nd order due to NUp & NUp
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUpNUp(m_field_Macro,A,D,ddF,Dmat_r_NUp,"DISP_MACRO",Dmat_rs_NUpNUp,"DISP_MACRO_r_NUp");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_NUp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_NUpNUp);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_NUpNUp);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUp); CHKERRQ(ierr);
            break;
          }
          case 10: {// 2nd order due to NUpz & NUpz
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUpzNUpz(m_field_Macro,A,D,ddF,Dmat_r_NUpz,"DISP_MACRO",Dmat_rs_NUpzNUpz,"DISP_MACRO_r_NUpz");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_NUpz(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUpz); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_NUpzNUpz);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_NUpzNUpz);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUpz); CHKERRQ(ierr);
            break;
          }
          case 11: {// 2nd order due to Ep & Ep
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EpEp(m_field_Macro,A,D,ddF,Dmat_r_Ep,"DISP_MACRO",Dmat_rs_EpEp,"DISP_MACRO_r_Ep");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_Ep(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Ep); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_EpEp);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_EpEp);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Ep); CHKERRQ(ierr);
            break;
          }
          case 12: {// 2nd order due to Ez & Ez
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EzEz(m_field_Macro,A,D,ddF,Dmat_r_Ez,"DISP_MACRO",Dmat_rs_EzEz,"DISP_MACRO_r_Ez");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_Ez(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Ez); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_EzEz);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_EzEz);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Ez); CHKERRQ(ierr);
            break;
          }
          case 13: {// 2nd order due to Gzp & Gzp
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_GzpGzp(m_field_Macro,A,D,ddF,Dmat_r_Gzp,"DISP_MACRO",Dmat_rs_GzpGzp,"DISP_MACRO_r_Gzp");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_Gzp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Gzp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_GzpGzp);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_GzpGzp);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Gzp); CHKERRQ(ierr);
            break;
          }
        }
        // post-processing
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ss_field.str(""); ss_field.clear();
          ss_field << "DISP_MACRO" << stochastic_fields[idx_disp];
          if (idx_disp<nvars) {
            
            ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
            
            //cout<<"First order derivative of dD"<<endl;
            //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
            
            ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);//ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            cout<<"Solving the first-order equation "<<ss_field.str().c_str()<<" is finish. \n";
            
            //cout<<"First order derivative of F"<<endl;
            //ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
            
          }
          else {
            ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
            
            ierr = KSPSolve(solver_Macro,ddF,ddD); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",ss_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            cout<<"Solving the second-order equation "<<ss_field.str().c_str()<<" is finish. \n";
          }
          //ierr = KSPReset(solver_Macro); CHKERRQ(ierr);
        }
        
        /***********************************************************************
         *
         * 3.2. Case 2: Applied forces are treated as random variables
         *
         **********************************************************************/
        
        if (idx_disp == 80) {
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = MatZeroEntries(A); CHKERRQ(ierr);
          
          // Establish an object of elastic FE method
          MyElasticFEMethod_Macro my_fe_Macro_r_F(m_field_Macro,A,dD,dF,Dmat,"DISP_MACRO");
          
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
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe_Macro_r_F);  CHKERRQ(ierr);
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe_Macro_r_F);  CHKERRQ(ierr);
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
          
          ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
          
          // Insert value into the force vector
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          //cout<<"First order derivative of F"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          // Solve the FE equation
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_F",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          //ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        }
        
      }
      
    /***************************************************************************
     *
     *  4. FINISH
     *
     **************************************************************************/
      
      //Destroy matrices
      ierr = VecDestroy(&F); CHKERRQ(ierr);
      ierr = VecDestroy(&dF); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = VecDestroy(&dD); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD); CHKERRQ(ierr);
      ierr = MatDestroy(&A); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver_Macro); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    
    
    
    // =========================================================================
    //
    //  A.VI. SOLUTION PHASE:
    //        Caculate RVE constitutive matrix Dmat
    //
    // =========================================================================
    
    virtual PetscErrorCode Calculate_RVEDmat_PSFE(FieldInterface &m_field_RVE,
                                             int &nvars, int &nders,
                                             vector<string> &stochastic_fields) {
      PetscFunctionBegin;
      cout <<"Hi from Calculate_RVEDmat"<<endl;
      
      ErrorCode rval;
      PetscErrorCode ierr;
      
      Dmat.resize(6,6); Dmat.clear();
      
      Dmat_r_Em.resize(6,6);   Dmat_r_Em.clear();
      Dmat_r_NUm.resize(6,6);   Dmat_r_NUm.clear();
      Dmat_r_Ep.resize(6,6);   Dmat_r_Ep.clear();
      Dmat_r_Ez.resize(6,6);   Dmat_r_Ez.clear();
      Dmat_r_NUp.resize(6,6);  Dmat_r_NUp.clear();
      Dmat_r_NUpz.resize(6,6); Dmat_r_NUpz.clear();
      Dmat_r_Gzp.resize(6,6);  Dmat_r_Gzp.clear();
      
      Dmat_rs_EmEm.resize(6,6);     Dmat_rs_EmEm.clear();
      Dmat_rs_NUmNUm.resize(6,6);     Dmat_rs_NUmNUm.clear();
      Dmat_rs_EpEp.resize(6,6);     Dmat_rs_EpEp.clear();
      Dmat_rs_EzEz.resize(6,6);     Dmat_rs_EzEz.clear();
      Dmat_rs_NUpNUp.resize(6,6);   Dmat_rs_NUpNUp.clear();
      Dmat_rs_NUpzNUpz.resize(6,6); Dmat_rs_NUpzNUpz.clear();
      Dmat_rs_GzpGzp.resize(6,6);   Dmat_rs_GzpGzp.clear();
      
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      Vec F1,F2,F3,F4,F5,F6,D1,D2,D3,D4,D5,D6;
      Vec dF1,dF2,dF3,dF4,dF5,dF6;
      Vec ddF1,ddF2,ddF3,ddF4,ddF5,ddF6;
      Vec dD1,dD2,dD3,dD4,dD5,dD6;
      Vec ddD1,ddD2,ddD3,ddD4,ddD5,ddD6;
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD6); CHKERRQ(ierr);
      
      
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
      
      RVEVolume MyRVEVol(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
      RVEVolumeTrans MyRVEVolTrans(m_field_RVE,Aij,D1,F1, RVE_volume_Vec);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVolTrans);  CHKERRQ(ierr);
      //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
      cout<<"Final RVE_volume = "<< RVE_volume <<endl;
      
      
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
      ierr = VecZeroEntries(D2); CHKERRQ(ierr);
      ierr = VecZeroEntries(D3); CHKERRQ(ierr);
      ierr = VecZeroEntries(D4); CHKERRQ(ierr);
      ierr = VecZeroEntries(D5); CHKERRQ(ierr);
      ierr = VecZeroEntries(D6); CHKERRQ(ierr);
      
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
      Vec Stress_Homo, Stress_Homo_r, Stress_Homo_rs;
      PetscScalar *avec;
      
      if(pcomm->rank()==0) {
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo_r);
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo_rs);
      } else {
        int ghost[] = {0,1,2,3,4,5};
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo_r);
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo_rs);
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
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D1,dF1,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D1,dF1,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D1,dF1,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D1,dF1,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D1,dF1,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF1); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF1); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF1,dD1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF1); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF1); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF1,ddD1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r1(m_field_RVE,Aij,dD1,dF1,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r1);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++){
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,0) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,0) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,0) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,0) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,0) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,0) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,0) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,0) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,0) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,0) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,0) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,0) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,0) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,0) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
        
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
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D2,dF2,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D2,dF2,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D2,dF2,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D2,dF2,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D2,dF2,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF2); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF2); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF2,dD2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF2); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF2); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF2,ddD2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r2(m_field_RVE,Aij,dD2,dF2,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r2);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,1) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,1) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,1) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,1) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,1) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,1) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,1) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,1) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,1) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,1) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,1) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,1) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,1) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,1) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
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
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D3,dF3,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D3,dF3,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D3,dF3,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D3,dF3,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D3,dF3,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF3); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF3); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF3,dD3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF3); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF3); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF3,ddD3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r3(m_field_RVE,Aij,dD3,dF3,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r3);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,2) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,2) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,2) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,2) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,2) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,2) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,2) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,2) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,2) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,2) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,2) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,2) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,2) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,2) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
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
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D4,dF4,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D4,dF4,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D4,dF4,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D4,dF4,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D4,dF4,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF4); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF4); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF4,dD4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF4); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF4); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF4,ddD4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r4(m_field_RVE,Aij,dD4,dF4,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r4);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,3) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,3) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,3) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,3) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,3) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,3) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,3) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,3) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,3) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,3) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,3) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,3) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,3) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,3) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
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
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D5,dF5,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D5,dF5,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D5,dF5,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D5,dF5,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D5,dF5,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF5); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF5); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF5,dD5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF5); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF5); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF5,ddD5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r5(m_field_RVE,Aij,dD5,dF5,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r5);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,4) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,4) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,4) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,4) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,4) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,4) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,4) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,4) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,4) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,4) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,4) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,4) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,4) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,4) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
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
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D6,dF6,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D6,dF6,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D6,dF6,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D6,dF6,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D6,dF6,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF6); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF6); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF6,dD6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF6); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF6); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF6,ddD6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r6(m_field_RVE,Aij,dD6,dF6,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r6);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,5) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,5) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,5) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,5) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,5) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,5) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,5) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,5) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,5) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,5) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,5) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,5) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,5) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,5) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
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
      
      ierr = VecDestroy(&dF1); CHKERRQ(ierr);
      ierr = VecDestroy(&dF2); CHKERRQ(ierr);
      ierr = VecDestroy(&dF3); CHKERRQ(ierr);
      ierr = VecDestroy(&dF4); CHKERRQ(ierr);
      ierr = VecDestroy(&dF5); CHKERRQ(ierr);
      ierr = VecDestroy(&dF6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&ddF1); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF2); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF3); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF4); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF5); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&D1); CHKERRQ(ierr);
      ierr = VecDestroy(&D2); CHKERRQ(ierr);
      ierr = VecDestroy(&D3); CHKERRQ(ierr);
      ierr = VecDestroy(&D4); CHKERRQ(ierr);
      ierr = VecDestroy(&D5); CHKERRQ(ierr);
      ierr = VecDestroy(&D6); CHKERRQ(ierr);
      
      
      ierr = VecDestroy(&dD1); CHKERRQ(ierr);
      ierr = VecDestroy(&dD2); CHKERRQ(ierr);
      ierr = VecDestroy(&dD3); CHKERRQ(ierr);
      ierr = VecDestroy(&dD4); CHKERRQ(ierr);
      ierr = VecDestroy(&dD5); CHKERRQ(ierr);
      ierr = VecDestroy(&dD6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&ddD1); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD2); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD3); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD4); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD5); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD6); CHKERRQ(ierr);
      
      ierr = MatDestroy(&Aij); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    
  };
  
}

#endif //__FE2_MACRO_SOLVER_HPP