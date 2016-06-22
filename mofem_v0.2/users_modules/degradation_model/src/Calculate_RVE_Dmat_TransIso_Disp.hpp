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

#ifndef __CALCULATE_RVE_DMAT_TRANSISO_DISP_HPP
#define __CALCULATE_RVE_DMAT_TRANSISO_DISP_HPP

namespace MoFEM {

  struct Calculate_RVE_Dmat_TransIso_Disp {
    
    struct MyVolumeFE: public VolumeElementForcesAndSourcesCore {
      MyVolumeFE(FieldInterface &_mField): VolumeElementForcesAndSourcesCore(_mField) {}
      
      
      int getRule(int order) { return -1; }; //with -1 this function will not work
      
      //This is the same funciton as used in the elastic element to make sure we use equal number of gauss point for calculation
      //of Dmat and subsequently in the elastic element.
      
      PetscErrorCode setGaussPts(int order) {
        //ublas::matrix<double> gaussPts;
        //gausPts.resize(nb_gauss_pts,4);
        //X,Y,Z,W
        order = 1;
        for(_IT_GET_FEDATA_BY_NAME_DOFS_FOR_LOOP_(this,"DISP_MACRO",dof)) {
          order = max(order,dof->get_max_order());
        }
        int rule = max(0,order-1);
        if( 2*rule + 1 < 2*(order-1) ) {
          SETERRQ2(PETSC_COMM_SELF,1,"wrong rule %d %d",order,rule);
        }
        int nb_gauss_pts = gm_rule_size(rule,3);
        if(gaussPts.size2() == (unsigned int)nb_gauss_pts) {
          PetscFunctionReturn(0);
        }
        gaussPts.resize(4,nb_gauss_pts);
        ierr = Grundmann_Moeller_integration_points_3D_TET(rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),&gaussPts(3,0)); CHKERRQ(ierr);
      }
    };
    MyVolumeFE feRhs; ///< cauclate right hand side for tetrahedral elements
    MyVolumeFE& getLoopFeRhs() { return feRhs; } ///< get rhs volume element

    FieldInterface &mField;
    Calculate_RVE_Dmat_TransIso_Disp(FieldInterface &m_field):
    feRhs(m_field),
    mField(m_field) {}

    
    struct CommonData {
      ublas::vector<double> wtAtGaussPts;
      //ublas::vector<double> tempAtGaussPts;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > Dmat_RVE;
    };
    CommonData commonData;
    
    PetscErrorCode addElasticElements(const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      ErrorCode rval;
      ierr = mField.add_finite_element("ELASTIC_FE_MACRO",MF_ZERO); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row("ELASTIC_FE_MACRO",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col("ELASTIC_FE_MACRO",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("ELASTIC_FE_MACRO",field_name); CHKERRQ(ierr);
      if(mField.check_field(mesh_nodals_positions)) {
        ierr = mField.modify_finite_element_add_field_data("ELASTIC_FE_MACRO",mesh_nodals_positions); CHKERRQ(ierr);
      }
      
      // loop over all blocksets and get data which name is MAT_ELASTICSET
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
        Range tEts;
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,tEts,true); CHKERR_PETSC(rval);
        ierr = mField.add_ents_to_finite_element_by_TETs(tEts,"ELASTIC_FE_MACRO"); CHKERRQ(ierr);
      }
      
      PetscFunctionReturn(0);
    }
    
    template<typename OP>
    struct OpGetFieldAtGaussPts: public OP::UserDataOperator {
      
      ublas::vector<double> &fieldAtGaussPts;
      OpGetFieldAtGaussPts(const string field_name,ublas::vector<double> &field_at_gauss_pts):
      OP::UserDataOperator(field_name,OP::UserDataOperator::OPROW),
      fieldAtGaussPts(field_at_gauss_pts) {}
      
      /** \brief operator calculating temperature and rate of temperature
       *
       * temperature temperature or rate of temperature is calculated multiplyingshape functions by degrees of freedom
       */
      PetscErrorCode doWork(
                            int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        try {
          //cout<<"form OpGetFieldAtGaussPts "<<endl;
          
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
            //cout<<"fieldAtGaussPts[gg] "<<fieldAtGaussPts[gg] <<endl;
          }
          
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
        
        PetscFunctionReturn(0);
      }
      
    };
    
    struct OpGetWtAtGaussPts: public OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore> {
      OpGetWtAtGaussPts(const string wt_field_name,CommonData &common_data):
      OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore>(wt_field_name,common_data.wtAtGaussPts) {}
    };
    
//    struct OpGetTempAtGaussPts: public OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore> {
//      OpGetTempAtGaussPts(const string temp_field_name,CommonData &common_data):
//      OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore>(temp_field_name,common_data.tempAtGaussPts) {}
//    };

  
    struct OpCalculate_RVEDmat: public VolumeElementForcesAndSourcesCore::UserDataOperator {
      
      PetscErrorCode ierr;
      FieldInterface &m_field_RVE;
      
      CommonData &commonData;
      OpCalculate_RVEDmat(FieldInterface &m_field_RVE, const string field_name,CommonData &common_data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name, UserDataOperator::OPROW),
      m_field_RVE(m_field_RVE),commonData(common_data){
      }
      
      ~OpCalculate_RVEDmat(){
      }


      PetscErrorCode doWork(
                            int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        try {
          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          int nb_gauss_pts = data.getN().size1();
//
          EntityHandle fe_ent = getMoFEMFEPtr()->get_ent(); //handle of finite element
//          cout<<"fe_ent "<<fe_ent <<endl;
          commonData.Dmat_RVE[fe_ent].resize(nb_gauss_pts);
//          map<EntityHandle, ublas::vector<ublas::matrix<double> > > Dmat_RVE;

          if(type == MBVERTEX) {  //the doWork loop is 1+6+4+1 times but we want to loop over Guass points only onece
            ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
            int field_rank=3; // it is mechanical problem
            applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();

            cout<<"nb_gauss_pts "<<nb_gauss_pts<<endl;
            
            for(int gg = 0;gg<nb_gauss_pts;gg++) {
              Vec F1,F2,F3,F4,F5,F6,D1;
              Mat A;
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F2); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F3); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F4); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F5); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F6); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&D1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              
              ierr = m_field_RVE.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_RVE",&A); CHKERRABORT(PETSC_COMM_WORLD,ierr);

              //=============================================================================================================
              //Calculation of RVE volume
              //=============================================================================================================
              double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
              Vec RVE_volume_Vec;
              ParallelComm* pcomm_RVE = ParallelComm::get_pcomm(&m_field_RVE.get_moab(),MYPCOMM_INDEX);
              //            cout<<" pcomm_RVE->size() = "<<pcomm_RVE->size()<<endl;
              ierr = VecCreateMPI(PETSC_COMM_SELF, 1, pcomm_RVE->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
              ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
              RVEVolume MyRVEVol(m_field_RVE,A,D1,F1,0.0,0.0, RVE_volume_Vec);
              
              //=============================================================================================================
              cout<<"\n\n";
              cout<<"gg Start =  "<<gg <<endl;
              cout<<"wt Start =  "<<commonData.wtAtGaussPts(gg) <<endl;
              //cout<<"temp Start =  "<<commonData.tempAtGaussPts(gg) <<endl;
              
              //We don't need to calculate internal forces for RVE, as ElasticFEMethod is used to assemble A matirx only
              //so noo need to create MyElasticFEMethod here
              ElasticFEMethod_Matrix my_fe_marix(m_field_RVE,A,D1,F1,0.0,0.0,commonData.wtAtGaussPts(gg),"DISP_RVE");
              TranIsotropicFibreDirRotElasticFEMethod my_fe_transiso(m_field_RVE,A,D1,F1,"DISP_RVE");
              ElasticFE_RVELagrange_Disp_Multi_Rhs MyFE_RVELagrange(m_field_RVE,A,D1,F1,F2,F3,F4,F5,F6,applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank);

//              cout<<"commonData.wtAtGaussPts(gg) =  "<<commonData.wtAtGaussPts(gg) <<endl;

              ierr = VecZeroEntries(F1); CHKERRQ(ierr);
              ierr = VecZeroEntries(F2); CHKERRQ(ierr);
              ierr = VecZeroEntries(F3); CHKERRQ(ierr);
              ierr = VecZeroEntries(F4); CHKERRQ(ierr);
              ierr = VecZeroEntries(F5); CHKERRQ(ierr);
              ierr = VecZeroEntries(F6); CHKERRQ(ierr);
              ierr = MatZeroEntries(A); CHKERRQ(ierr);
              ierr = VecZeroEntries(D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_marix);  CHKERRQ(ierr);
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_transiso);  CHKERRQ(ierr);
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVELagrange);  CHKERRQ(ierr);

              ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
              ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
              
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

              if(gg==0){//do it for only one Gauss point as it is same for all others
                //=============================================================================================================
                ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
                ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVol);  CHKERRQ(ierr);

                //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
                ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
//                cout<<"Final RVE_volume = "<< RVE_volume <<endl;
//                cout<<"Element Number fe_ent"<<fe_ent<<endl;
//                string wait;
//                cin >>wait;
                //=============================================================================================================
              }
              
              //Solver
              KSP solver;
              ierr = KSPCreate(PETSC_COMM_SELF,&solver); CHKERRQ(ierr);
              ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
              ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
              ierr = KSPSetUp(solver); CHKERRQ(ierr);
//              ierr = VecView(D1,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

              //=============================================================================================================
              // homogenised stress for strian [1 0 0 0 0 0]^T
              //=============================================================================================================
              //solve for F1
              ierr = KSPSolve(solver,F1,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//              ierr = VecView(D1,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

              ublas::matrix<FieldData> Dmat;
              Dmat.resize(6,6); Dmat.clear();
              
              //create a vector for 6 components of homogenized stress
              Vec Stress_Homo;
              if(pcomm_RVE->rank()==0) {
                VecCreateGhost(PETSC_COMM_SELF,6,6,0,PETSC_NULL,&Stress_Homo);
              } else {
                int ghost[] = {0,1,2,3,4,5};
                VecCreateGhost(PETSC_COMM_SELF,0,6,6,ghost,&Stress_Homo);
                
              }
              

              {
                ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
                
                ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);

                
                ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
                PetscScalar *avec;
                VecGetArray(Stress_Homo, &avec);
                for(int ii=0; ii<6; ii++){
                  Dmat(ii,0)=*avec;
                  avec++;
                }
//                if(pcomm_RVE->rank()==0){
//                  cout<< "\nStress_Homo = \n\n";
//                  for(int ii=0; ii<6; ii++){
//                    cout <<Dmat(ii,0)<<endl;
//                  }
//                }
                
              }
              
              
              //=============================================================================================================
              // homogenised stress for strian [0 1 0 0 0 0]^T
              //=============================================================================================================
              //solve for F2
              ierr = KSPSolve(solver,F2,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              
              {
                ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
                //    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
                
                ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
                ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
                
                //    if(pcomm->rank() == 0) cout<< " Stress_Homo =  "<<endl;
                //    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
                PetscScalar *avec;
                VecGetArray(Stress_Homo, &avec);
                for(int ii=0; ii<6; ii++){
                  Dmat(ii,1)=*avec;
                  avec++;
                }
                
//                if(pcomm_RVE->rank()==0){
//                  cout<< "\nStress_Homo = \n\n";
//                  for(int ii=0; ii<6; ii++){
//                    cout <<Dmat(ii,1)<<endl;
//                  }
//                }
              }

              //=============================================================================================================
              // homogenised stress for strian [0 0 1 0 0 0]^T
              //=============================================================================================================
              //solve for F3
              ierr = KSPSolve(solver,F3,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              
              {
                ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
                
                ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
                ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
                
                //    if(pcomm->rank() == 0) cout<< " Stress_Homo =  "<<endl;
                //    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
                PetscScalar *avec;
                VecGetArray(Stress_Homo, &avec);
                for(int ii=0; ii<6; ii++){
                  Dmat(ii,2)=*avec;
                  avec++;
                }
                
//                if(pcomm_RVE->rank()==0){
//                  cout<< "\nStress_Homo = \n\n";
//                  for(int ii=0; ii<6; ii++){
//                    cout <<Dmat(ii,2)<<endl;
//                  }
//                }
              }
              //=============================================================================================================
              // homogenised stress for strian [0 0 0 1 0 0]^T
              //=============================================================================================================
              //solve for F4
              ierr = KSPSolve(solver,F4,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              
              
              {
                ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
                ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
                ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
                
                //    if(pcomm->rank() == 0) cout<< " Stress_Homo =  "<<endl;
                //    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
                PetscScalar *avec;
                VecGetArray(Stress_Homo, &avec);
                for(int ii=0; ii<6; ii++){
                  Dmat(ii,3)=*avec;
                  avec++;
                }
                
//                if(pcomm_RVE->rank()==0){
//                  cout<< "\nStress_Homo = \n\n";
//                  for(int ii=0; ii<6; ii++){
//                    cout <<Dmat(ii,3)<<endl;
//                  }
//                }
              }
              //=============================================================================================================
              // homogenised stress for strian [0 0 0 0 1 0]^T
              //=============================================================================================================
              //solve for F5
              ierr = KSPSolve(solver,F5,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              
              {
                ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
                
                ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
                ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
                
                //    if(pcomm->rank() == 0) cout<< " Stress_Homo =  "<<endl;
                //    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
                PetscScalar *avec;
                VecGetArray(Stress_Homo, &avec);
                for(int ii=0; ii<6; ii++){
                  Dmat(ii,4)=*avec;
                  avec++;
                }
                
//                if(pcomm_RVE->rank()==0){
//                  cout<< "\nStress_Homo = \n\n";
//                  for(int ii=0; ii<6; ii++){
//                    cout <<Dmat(ii,4)<<endl;
//                  }
//                }
              }
              //=============================================================================================================
              // homogenised stress for strian [0 0 0 0 0 1]^T
              //=============================================================================================================
              //solve for F6
              ierr = KSPSolve(solver,F6,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              
              {
                ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
                ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
                ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
                
                //    if(pcomm->rank() == 0) cout<< " Stress_Homo =  "<<endl;
                //    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
                PetscScalar *avec;
                VecGetArray(Stress_Homo, &avec);
                for(int ii=0; ii<6; ii++){
                  Dmat(ii,5)=*avec;
                  avec++;
                }
//                if(pcomm_RVE->rank()==0){
//                  cout<< "\nStress_Homo = \n\n";
//                  for(int ii=0; ii<6; ii++){
//                    cout <<Dmat(ii,5)<<endl;
//                  }
//                }
              }
              //=============================================================================================================
             
//              if(pcomm_RVE->rank()==0){
//                cout<< "\nStress_Homo = \n\n";
//                for(int ii=0; ii<6; ii++){
//                  for(int jj=0; jj<6; jj++){
//                    cout <<Dmat(ii,jj)<<"    ";
//                  }
//                  cout<<endl;
//                }
//              }
//              
              cout<<"fe_ent "<<fe_ent << "gg  "<<gg<<endl;
              cout<<"Dmat = "<<Dmat<<endl;

              commonData.Dmat_RVE[fe_ent](gg).resize(6,6);
              commonData.Dmat_RVE[fe_ent](gg)=Dmat;
              
              ierr = VecDestroy(&F1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&F2); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&F3); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&F4); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&F5); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&F6); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&D1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&Stress_Homo); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&RVE_volume_Vec); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              
              ierr = MatDestroy(&A); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = KSPDestroy(&solver); CHKERRQ(ierr);
              
//              string wait;
//              cin>>wait;
            }
          }
          
//          cout<<"fe_ent "<<fe_ent <<endl;
//          cout<<"nb_gauss_pts "<<nb_gauss_pts <<endl;
//          for(int gg = 0;gg<nb_gauss_pts;gg++) {
//            cout<<"Gauss Number =  "<<gg <<endl;
//            cout<<"commonData.Dmat_RVE[fe_ent](gg) =  "<<commonData.Dmat_RVE[fe_ent](gg) <<endl;
//          }
 
          
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        PetscFunctionReturn(0);
      }
      
    };
    
    PetscErrorCode setRVE_DmatRhsOperators(FieldInterface &m_field_RVE, string field_name,string wt_field_name){//,string temp_field_name) {
      PetscFunctionBegin;
      //first calculate wt at each gauss point
      feRhs.getOpPtrVector().push_back(new OpGetWtAtGaussPts(wt_field_name,commonData));
      //feRhs.getOpPtrVector().push_back(new OpGetTempAtGaussPts(temp_field_name,commonData));
      //At each gauss point run RVE with its own mesh
      feRhs.getOpPtrVector().push_back(new OpCalculate_RVEDmat(m_field_RVE,field_name,commonData));

      PetscFunctionReturn(0);
    }

  
  };
  
}

#endif //__CALCULATE_RVE_DMAT_TRANSISO_DISP_HPP



