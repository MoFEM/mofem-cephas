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

#ifndef __CALCULATE_RVE_DMAT_HPP
#define __CALCULATE_RVE_DMAT_HPP

namespace MoFEM {

  struct Calculate_RVE_Dmat {
    
    struct MyVolumeFE: public TetElementForcesAndSourcesCore {
      MyVolumeFE(FieldInterface &_mField): TetElementForcesAndSourcesCore(_mField) {}
      int getRule(int order) { return order-1; };
    };
    MyVolumeFE feRhs; ///< cauclate right hand side for tetrahedral elements
    MyVolumeFE& getLoopFeRhs() { return feRhs; } ///< get rhs volume element

    FieldInterface &mField;
    Calculate_RVE_Dmat(FieldInterface &m_field):
    feRhs(m_field),
    mField(m_field) {}

    
    struct CommonData {
      ublas::vector<double> wtAtGaussPts;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > Dmat_RVE;
    };
    CommonData commonData;

    
    
    template<typename OP>
    struct OpGetFieldAtGaussPts: public OP::UserDataOperator {
      
      ublas::vector<double> &fieldAtGaussPts;
      OpGetFieldAtGaussPts(const string field_name,ublas::vector<double> &field_at_gauss_pts):
      OP::UserDataOperator(field_name),
      fieldAtGaussPts(field_at_gauss_pts) {}
      
      /** \brief operator calculating temperature and rate of temperature
       *
       * temperature temperature or rate of temperature is calculated multiplyingshape functions by degrees of freedom
       */
      PetscErrorCode doWork(
                            int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        try {
//          cout<<"form OpGetFieldAtGaussPts "<<endl;
          
          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          int nb_dofs = data.getFieldData().size();
          int nb_gauss_pts = data.getN().size1();
          
          //initialize
          fieldAtGaussPts.resize(nb_gauss_pts);
          if(type == MBVERTEX) {
            //loop over shape functions on entities allways start from
            //vertices, so if nodal shape functions are processed, vector of
            //field values is zeroad at initialization
            fill(fieldAtGaussPts.begin(),fieldAtGaussPts.end(),0);
          }
          
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            fieldAtGaussPts[gg] += inner_prod(data.getN(gg,nb_dofs),data.getFieldData());
//            cout<<"fieldAtGaussPts[gg] "<<fieldAtGaussPts[gg] <<endl;
          }
          
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
        
        PetscFunctionReturn(0);
      }
      
    };

    
    struct OpGetWtAtGaussPts: public OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore> {
      OpGetWtAtGaussPts(const string wt_field_name,CommonData &common_data):
      OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore>(wt_field_name,common_data.wtAtGaussPts) {}
    };

  
  

    
    
    struct OpCalculate_RVEDmat: public TetElementForcesAndSourcesCore::UserDataOperator {
      
      PetscErrorCode ierr;
      FieldInterface &m_field_RVE;
      
      CommonData &commonData;
      OpCalculate_RVEDmat(FieldInterface &m_field_RVE, const string field_name,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      m_field_RVE(m_field_RVE),commonData(common_data){
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F2); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F3); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F4); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F5); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F6); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&D1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&D2); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&D3); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&D4); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&D5); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&D6); CHKERRABORT(PETSC_COMM_WORLD,ierr);

        ierr = m_field_RVE.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_RVE",&A); CHKERRABORT(PETSC_COMM_WORLD,ierr);
      }
      
      ~OpCalculate_RVEDmat(){
        ierr = VecDestroy(&F1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = VecDestroy(&F2); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = VecDestroy(&F3); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = VecDestroy(&F4); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = VecDestroy(&F5); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = VecDestroy(&F6); CHKERRABORT(PETSC_COMM_WORLD,ierr);

        ierr = VecDestroy(&D1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = VecDestroy(&D2); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = VecDestroy(&D3); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = VecDestroy(&D4); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = VecDestroy(&D5); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = VecDestroy(&D6); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        
        ierr = MatDestroy(&A); CHKERRABORT(PETSC_COMM_WORLD,ierr);
      }
      
    
      struct MyElasticFEMethod: public ElasticFEMethod {
        MyElasticFEMethod(FieldInterface& _m_field,Mat _Aij,Vec _D,Vec& _F,double _lambda,double _mu, string _field_name):
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

      Vec F1,F2,F3,F4,F5,F6,D1,D2,D3,D4,D5,D6;
      Mat A;

      
      
      PetscErrorCode doWork(
                            int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        try {
          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          int nb_gauss_pts = data.getN().size1();
          
          EntityHandle fe_ent = getMoFEMFEPtr()->get_ent(); //handle of finite element
          cout<<"fe_ent "<<fe_ent <<endl;
          
          
          if(type == MBVERTEX) {  //the doWork loop is 1+6+4+1 times but we want to loop over Guass points only onece
            const double young_modulus = 1;
            const double poisson_ratio = 0.0;
            MyElasticFEMethod my_fe(m_field_RVE,A,D1,F1,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_RVE");
            
            ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
            int field_rank=3; // it is mechanical problem
            applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
            ElasticFE_RVELagrange_Disp_Multi_Rhs MyFE_RVELagrange(m_field_RVE,A,D1,F1,F2,F3,F4,F5,F6,applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank);

            //=============================================================================================================
            //Calculation of RVE volume
            //=============================================================================================================
            double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
            Vec RVE_volume_Vec;
            ParallelComm* pcomm_RVE = ParallelComm::get_pcomm(&m_field_RVE.get_moab(),MYPCOMM_INDEX);
            ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm_RVE->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
            ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
            RVEVolume MyRVEVol(m_field_RVE,A,D1,F1,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio), RVE_volume_Vec);
            //=============================================================================================================

            for(int gg = 0;gg<nb_gauss_pts;gg++) {
              cout<<"gg Start =  "<<gg <<endl;
              
              
              ierr = VecZeroEntries(F1); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(F1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(F1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              
              ierr = VecZeroEntries(F2); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(F2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(F2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              
              ierr = VecZeroEntries(F1); CHKERRQ(ierr);
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
              
              ierr = MatZeroEntries(A); CHKERRQ(ierr);

              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe);  CHKERRQ(ierr);
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVELagrange);  CHKERRQ(ierr);
              
              ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
              ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
              
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

              if(gg==0){//do it for only one Gauss point as it is same for all others
                //=============================================================================================================
                ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
                //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
                ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
                cout<<"Final RVE_volume = "<< RVE_volume <<endl;
                //=============================================================================================================
              }
              
              //Solver
              KSP solver;
              ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
              ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
              ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
              ierr = KSPSetUp(solver); CHKERRQ(ierr);

              //solve for F1 and D1
              ierr = KSPSolve(solver,F1,D1); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(D1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(D1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

              //=============================================================================================================
              // homogenised stress for strian [1 0 0 0 0 0]^T
              //=============================================================================================================
              
              ublas::matrix<FieldData> Dmat;
              Dmat.resize(6,6); Dmat.clear();
              
              //create a vector for 6 components of homogenized stress
              Vec Stress_Homo;
              if(pcomm_RVE->rank()==0) {
                VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo);
              } else {
                int ghost[] = {0,1,2,3,4,5};
                VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo);
                
              }
              
              
              {
                ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
                
                ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
                
                
                ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
                VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
                VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
                
                PetscScalar *avec;
                VecGetArray(Stress_Homo, &avec);
                for(int ii=0; ii<6; ii++){
                  Dmat(ii,0)=*avec;
                  avec++;
                }
                
                if(pcomm_RVE->rank()==0){
                  cout<< "\nStress_Homo = \n\n";
                  for(int ii=0; ii<6; ii++){
                    cout <<Dmat(ii,0)<<endl;
                  }
                }
                
              }
              //=============================================================================================================

              cout<<"gg End =  "<<gg <<endl;

              string wait;
              cin >> wait;
            
              
              
            }
          }
          
          
          
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        PetscFunctionReturn(0);
      }
      
    };
    
    
    
    
    
    PetscErrorCode setRVE_DmatRhsOperators(FieldInterface &m_field_RVE, string field_name,string wt_field_name) {
      PetscFunctionBegin;
      //first calculate wt at each gauss point
      feRhs.get_op_to_do_Rhs().push_back(new OpGetWtAtGaussPts(wt_field_name,commonData));
      //At each gauss point run RVE with its own mesh
      feRhs.get_op_to_do_Rhs().push_back(new OpCalculate_RVEDmat(m_field_RVE,field_name,commonData));

      
      
      
//      map<int,BlockData>::iterator sit = setOfBlocks.begin();
//      feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(field_name,commonData));
//      for(;sit!=setOfBlocks.end();sit++) {
//        //add finite element
//        feRhs.get_op_to_do_Rhs().push_back(new OpThermalRhs(field_name,F,sit->second,commonData));
//      }
      PetscFunctionReturn(0);
    }

    
    
    
  };
  
}

#endif //__CALCULATE_RVE_DMAT_HPP



