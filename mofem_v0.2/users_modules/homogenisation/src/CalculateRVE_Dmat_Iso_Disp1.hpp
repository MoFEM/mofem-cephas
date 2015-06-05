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

#ifndef __CALCULATE_RVE_DMAT_ISO_DISP1_HPP
#define __CALCULATE_RVE_DMAT_ISO_DISP1_HPP

namespace MoFEM {

  struct Calculate_RVE_Dmat_iso_Disp1 {
    
    
    //global variable Dmat
    ublas::matrix<double> Dmat;

    
    
    
    PetscErrorCode Calculate_RVEDmat(FieldInterface &mField_RVE) {
      PetscFunctionBegin;
      cout <<"Hi from Calculate_RVE_Dmat_iso_Disp1 "<<endl;

      ErrorCode rval;
      PetscErrorCode ierr;

      Vec F1,F2,F3,F4,F5,F6,D1,D2,D3,D4,D5,D6;
      ierr = mField_RVE.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F1); CHKERRQ(ierr);
      ierr = mField_RVE.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F2); CHKERRQ(ierr);
      ierr = mField_RVE.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F3); CHKERRQ(ierr);
      ierr = mField_RVE.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F4); CHKERRQ(ierr);
      ierr = mField_RVE.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F5); CHKERRQ(ierr);
      ierr = mField_RVE.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F6); CHKERRQ(ierr);

      ierr = mField_RVE.VecCreateGhost("ELASTIC_MECHANICS",COL,&D1); CHKERRQ(ierr);
      ierr = mField_RVE.VecCreateGhost("ELASTIC_MECHANICS",COL,&D2); CHKERRQ(ierr);
      ierr = mField_RVE.VecCreateGhost("ELASTIC_MECHANICS",COL,&D3); CHKERRQ(ierr);
      ierr = mField_RVE.VecCreateGhost("ELASTIC_MECHANICS",COL,&D4); CHKERRQ(ierr);
      ierr = mField_RVE.VecCreateGhost("ELASTIC_MECHANICS",COL,&D5); CHKERRQ(ierr);
      ierr = mField_RVE.VecCreateGhost("ELASTIC_MECHANICS",COL,&D6); CHKERRQ(ierr);
      
      
      //create matrices (here F, D and Aij are matrices for the full problem)
      Mat Aij;
      ierr = mField_RVE.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

      struct MyElasticFEMethod: public ElasticFEMethod {
        MyElasticFEMethod(FieldInterface& _m_field,Mat _Aij,Vec _D,Vec& _F,double _lambda,double _mu):
        ElasticFEMethod(_m_field,_Aij,_D,_F,_lambda,_mu) {};

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
      const double young_modulus = 1;
      const double poisson_ratio = 0.0;
      int field_rank=3;

      ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
      applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();

      MyElasticFEMethod my_fe(mField_RVE,Aij,D1,F1,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
      ElasticFE_RVELagrange_Disp_Multi_Rhs MyFE_RVELagrange(mField_RVE,Aij,D1,F1,F2,F3,F4,F5,F6,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank);

      cout<<"After ElasticFE_RVELagrange_Disp_Multi_Rhs = "<<endl;
      ierr = VecZeroEntries(F1); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = VecZeroEntries(F2); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = VecZeroEntries(F3); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = VecZeroEntries(F4); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = VecZeroEntries(F5); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = VecZeroEntries(F6); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
      
      //  cout<<"Before  my_fe = "<<endl;
      ierr = mField_RVE.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",my_fe);  CHKERRQ(ierr);
      //  cout<<"Before  MyFE_RVELagrange = "<<endl;
      ierr = mField_RVE.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVELagrange);  CHKERRQ(ierr);
      //  cout<<"After  MyFE_RVELagrange = "<<endl;

      ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

      ierr = VecGhostUpdateBegin(F1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F1); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F1); CHKERRQ(ierr);

      ierr = VecGhostUpdateBegin(F2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F2); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F2); CHKERRQ(ierr);

      ierr = VecGhostUpdateBegin(F3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F3); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F3); CHKERRQ(ierr);

      ierr = VecGhostUpdateBegin(F4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F4); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F4); CHKERRQ(ierr);

      ierr = VecGhostUpdateBegin(F5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F5); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F5); CHKERRQ(ierr);

      ierr = VecGhostUpdateBegin(F6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F6); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F6); CHKERRQ(ierr);

      //=============================================================================================================
      //Calculation of RVE volume for homogenised stress calculaiton
      //=============================================================================================================
      double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
      Vec RVE_volume_Vec;
      ParallelComm* pcomm = ParallelComm::get_pcomm(&mField_RVE.get_moab(),MYPCOMM_INDEX);
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
      ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);

      RVEVolume MyRVEVol(mField_RVE,Aij,D1,F1,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio), RVE_volume_Vec);
      ierr = mField_RVE.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyRVEVol);  CHKERRQ(ierr);
      //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
      cout<<"Final RVE_volume = "<< RVE_volume <<endl;
      //=============================================================================================================

      //Solver
      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      ierr = KSPSetUp(solver); CHKERRQ(ierr);

      //solve for F1 and D1
      ierr = KSPSolve(solver,F1,D1); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = mField_RVE.set_global_ghost_vector("ELASTIC_MECHANICS",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

      //=============================================================================================================
      // homogenised stress for strian [1 0 0 0 0 0]^T
      //=============================================================================================================

//      ublas::matrix<FieldData> Dmat;
      Dmat.resize(6,6); Dmat.clear();

      //create a vector for 6 components of homogenized stress
      Vec Stress_Homo;
      if(pcomm->rank()==0) {
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo);
      } else {
        int ghost[] = {0,1,2,3,4,5};
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo);

      }


      {
        ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);

        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(mField_RVE,Aij,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank);


        ierr = mField_RVE.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
        VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
        VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);

    //    if(pcomm->rank() == 0) cout<< " Stress_Homo =  "<<endl;
    //    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

        PetscScalar *avec;
        VecGetArray(Stress_Homo, &avec);
        for(int ii=0; ii<6; ii++){
          Dmat(ii,0)=*avec;
          avec++;
        }
        
        if(pcomm->rank()==0){
          cout<< "\nStress_Homo = \n\n";
          for(int ii=0; ii<6; ii++){
            cout <<Dmat(ii,0)<<endl;
          }
        }
        
      }
    //=============================================================================================================


    //solve for F2 and D2
    ierr = KSPSolve(solver,F2,D2); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = mField_RVE.set_global_ghost_vector("ELASTIC_MECHANICS",ROW,D2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

    //=============================================================================================================
    // homogenised stress for strian [0 1 0 0 0 0]^T
    //=============================================================================================================

    {
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
  //    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(mField_RVE,Aij,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
      ierr = mField_RVE.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);

      VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
      VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
  //    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

      //    if(pcomm->rank() == 0) cout<< " Stress_Homo =  "<<endl;
  //    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      PetscScalar *avec;
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
    }
      //=============================================================================================================
    //solve for F3 and D3
    ierr = KSPSolve(solver,F3,D3); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = mField_RVE.set_global_ghost_vector("ELASTIC_MECHANICS",ROW,D3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    //=============================================================================================================
    // homogenised stress for strian [0 0 1 0 0 0]^T
    //=============================================================================================================
    {
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(mField_RVE,Aij,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
      ierr = mField_RVE.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);

      VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
      VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
      PetscScalar *avec;
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
    }
  //  //=============================================================================================================
    //solve for F4 and D4
    ierr = KSPSolve(solver,F4,D4); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = mField_RVE.set_global_ghost_vector("ELASTIC_MECHANICS",ROW,D4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    //=============================================================================================================
    // homogenised stress for strian [0 0 0 1 0 0]^T
    //=============================================================================================================
    {
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(mField_RVE,Aij,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
      ierr = mField_RVE.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);

      VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
      VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
      PetscScalar *avec;
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
    }
    //=============================================================================================================
    //solve for F5 and D5
    ierr = KSPSolve(solver,F5,D5); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = mField_RVE.set_global_ghost_vector("ELASTIC_MECHANICS",ROW,D5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  //  //=============================================================================================================
  //  // homogenised stress for strian [0 0 0 0 1 0]^T
  //  //=============================================================================================================

    {
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(mField_RVE,Aij,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
      ierr = mField_RVE.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);

      VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
      VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
      PetscScalar *avec;
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
    }
    //=============================================================================================================
    //solve for F6 and D6
    ierr = KSPSolve(solver,F6,D6); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = mField_RVE.set_global_ghost_vector("ELASTIC_MECHANICS",ROW,D6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

    //=============================================================================================================
    // homogenised stress for strian [0 0 0 0 0 1]^T
    //=============================================================================================================
    {
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(mField_RVE,Aij,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
      ierr = mField_RVE.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
      
      VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
      VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
      PetscScalar *avec;
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
        cout<< "Dmat = "<<Dmat<<endl;
      }
    }
    
    
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

  
  };
}

#endif //__CALCULATE_RVE_DMAT_ISO_DISP_HPP



