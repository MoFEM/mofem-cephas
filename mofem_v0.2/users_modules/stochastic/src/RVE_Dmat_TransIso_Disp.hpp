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

#ifndef __RVE_DMAT_TRANSISO_DISP_HPP
#define __RVE_DMAT_TRANSISO_DISP_HPP

namespace MoFEM {
  
  struct RVE_Dmat_TransIso_Disp {
    
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
    
    //
    // Caculate RVE constitutive matrix Dmat
    //
    PetscErrorCode Calculate_RVEDmat(FieldInterface &m_field_RVE,
                                     int &nvars, int &nders,
                                     vector<string> &stochastic_fields) {
      PetscFunctionBegin;
      cout <<"Hi from Calculate_RVE_Dmat_iso_Disp1 "<<endl;
      
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
      
      
      /*****************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ****************************************************************************/
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
      
      
      /*****************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ****************************************************************************/
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
      double YoungModulus = 1;
      double PoissonRatio = 0.0;
      double alpha;
      int field_rank=3;
      
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)) {
        cout << endl << *it << endl;
        
        //Get block name
        string name = it->get_name();
        if (name.compare(0,13,"MAT_ELASTIC_1") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << mydata;
          YoungModulus=mydata.data.Young;
          PoissonRatio=mydata.data.Poisson;
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
        
      if (pcomm->rank()==0) {
        cout<< "\nStress_Homo = \n\n";
        for(int ii=0; ii<6; ii++){
          cout <<Dmat(ii,0)<<endl;
        }
      }
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
          
          cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++){
            cout.precision(15);
            cout<<*avec_r<<endl;
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
      
      if(pcomm->rank()==0){
        cout<< "\nStress_Homo = \n\n";
        for(int ii=0; ii<6; ii++){
          cout <<Dmat(ii,1)<<endl;
        }
      }
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
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D1,dF2,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D1,dF2,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D1,dF2,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D1,dF2,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D1,dF2,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D1,dF2,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D1,dF2,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D1,ddF2,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D1,ddF2,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D1,ddF2,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D1,ddF2,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D1,ddF2,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D1,ddF2,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D1,ddF2,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
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
          
          cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            cout<<*avec_r<<endl;
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
      
      if(pcomm->rank()==0){
        cout<< "\nStress_Homo = \n\n";
        for(int ii=0; ii<6; ii++){
          cout <<Dmat(ii,2)<<endl;
        }
      }
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
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D1,dF3,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D1,dF3,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D1,dF3,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D1,dF3,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D1,dF3,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D1,dF3,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D1,dF3,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D1,ddF3,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D1,ddF3,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D1,ddF3,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D1,ddF3,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D1,ddF3,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D1,ddF3,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D1,ddF3,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
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
          
          cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            cout<<*avec_r<<endl;
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
      
      if(pcomm->rank()==0){
        cout<< "\nStress_Homo = \n\n";
        for(int ii=0; ii<6; ii++){
          cout <<Dmat(ii,3)<<endl;
        }
      }
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
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D1,dF4,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D1,dF4,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D1,dF4,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D1,dF4,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D1,dF4,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D1,dF4,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D1,dF4,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D1,ddF4,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D1,ddF4,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D1,ddF4,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D1,ddF4,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D1,ddF4,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D1,ddF4,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D1,ddF4,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
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
          
          cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            cout<<*avec_r<<endl;
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
      
      if(pcomm->rank()==0){
        cout<< "\nStress_Homo = \n\n";
        for(int ii=0; ii<6; ii++){
          cout <<Dmat(ii,4)<<endl;
        }
      }
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
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D1,dF5,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D1,dF5,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D1,dF5,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D1,dF5,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D1,dF5,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D1,dF5,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D1,dF5,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D1,ddF5,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D1,ddF5,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D1,ddF5,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D1,ddF5,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D1,ddF5,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D1,ddF5,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D1,ddF5,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
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
          
          cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            cout<<*avec_r<<endl;
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
      
      if(pcomm->rank()==0){
        cout<< "\nStress_Homo = \n\n";
        for(int ii=0; ii<6; ii++){
          cout <<Dmat(ii,5)<<endl;
        }
        cout<< "\n\n";
      }
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
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D1,dF6,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D1,dF6,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D1,dF6,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D1,dF6,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D1,dF6,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D1,dF6,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D1,dF6,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D1,ddF6,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D1,ddF6,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D1,ddF6,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D1,ddF6,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D1,ddF6,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D1,ddF6,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D1,ddF6,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
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
          
          cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            cout<<*avec_r<<endl;
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
      
      
      /*****************************************************************************
       *  3.2. SOLVE THE FIRST-ORDER AND THE SECOND-ORDER FE EQUILIBRIUM EQUATION
       *     1st order-[K][U_r] = -[K_r][U}
       *     2nd order-[K][U_rs] = -[K_rs][U]-2[K_r][U_s]
       *
       ****************************************************************************/
      
      
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

#endif //__CALCULATE_RVE_DMAT_ISO_DISP_HPP