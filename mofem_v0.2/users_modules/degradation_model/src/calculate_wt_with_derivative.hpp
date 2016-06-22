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

#ifndef __CALCULATE_WT_WITH_DERIVATIVE_HPP
#define __CALCULATE_WT_WITH_DERIVATIVE_HPP

namespace MoFEM {
  
  
  struct Calculate_wt_with_derivative {
    
    struct MyVolumeFE: public VolumeElementForcesAndSourcesCore {
      MyVolumeFE(FieldInterface &_mField): VolumeElementForcesAndSourcesCore(_mField) {}
      int getRule(int order) { return order-1; };
    };
    
    
    MyVolumeFE feRhsDegradation; ///< cauclate right hand side for tetrahedral elements
    MyVolumeFE& getLoopFeRhs() { return feRhsDegradation; } ///< get rhs volume element
    
    MyVolumeFE feLhsDegradation; //< calculate left hand side for tetrahedral elements
    MyVolumeFE& getLoopFeLhs() { return feLhsDegradation; } ///< get lhs volume element

    
    FieldInterface &mField;
    Calculate_wt_with_derivative(FieldInterface &m_field): feRhsDegradation(m_field), feLhsDegradation(m_field), mField(m_field) {}
    
    struct BlockData {
      double saturation_Conc;
      BlockData(): saturation_Conc(100) {}
      Range tEts;
    };
    map<int,BlockData> setOfBlocks;
    
    struct CommonData {
      ublas::vector<double> temperatureAtGaussPts;
      ublas::vector<double> concentrationAtGaussPts;
      ublas::vector<double> fieldAtGaussPts;
      ublas::vector<double> fieldRateAtGaussPts;
      ublas::vector<double> field_r_BetaAtGaussPts;
    };
    CommonData commonData;
    
    
    
  struct LoadTimeSeries: public FEMethod {
    FieldInterface& mField;
    const string tempSeries;
    const string moisSeries;
    LoadTimeSeries(FieldInterface& _mField, const string temp_series, const string mois_series):mField(_mField), tempSeries(temp_series), moisSeries(mois_series) {}

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      SeriesRecorder *recorder_ptr;
      ierr = mField.query_interface(recorder_ptr); CHKERRQ(ierr);

//      cout << "ts_step "<<ts_step<< endl;

      ierr = recorder_ptr->load_series_data(tempSeries,ts_step); CHKERRQ(ierr);
      ierr = recorder_ptr->load_series_data(moisSeries,ts_step); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }
  };

    
    
    /** \brief this calas is to control time stepping
     * \infroup mofem_thermal_elem
     *
     * It is used to save data for temerature rate vectot to MoFEM field.
     */
    struct UpdateAndControl: public FEMethod {
      
      FieldInterface& mField;
      TS tS;
      const string tempName;
      const string rateName;
      int jacobianLag;
      UpdateAndControl(FieldInterface& _mField,TS _ts,
                       const string temp_name,
                       const string rate_name): mField(_mField),tS(_ts),
      tempName(temp_name),rateName(rate_name),jacobianLag(-1) {}
      
      PetscErrorCode preProcess() {
        PetscFunctionBegin;
        PetscErrorCode ierr;
        ierr = mField.set_other_local_ghost_vector(problemPtr,tempName,rateName,
                                                   ROW,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        PetscFunctionReturn(0);
      }
      
      PetscErrorCode postProcess() {
        PetscFunctionBegin;
        PetscErrorCode ierr;
        SNES snes;
        ierr = TSGetSNES(tS,&snes); CHKERRQ(ierr);
        ierr = SNESSetLagJacobian(snes,jacobianLag); CHKERRQ(ierr);
        PetscFunctionReturn(0);
      }
      
    };

    
    /** \brief TS monitore it records temperature at time steps
     * \infroup mofem_thermal_elem
     */
    struct TimeSeriesMonitor: public FEMethod {
      
      FieldInterface &mField;
      const string seriesName;
      const string tempName;
      BitRefLevel mask;
      
      TimeSeriesMonitor(FieldInterface &m_field,const string series_name,const string temp_name):
      mField(m_field),seriesName(series_name),tempName(temp_name) {
        mask.set();
      }
      
      PetscErrorCode postProcess() {
        PetscFunctionBegin;
        PetscErrorCode ierr;
        
        ierr = mField.set_global_ghost_vector(problemPtr,ROW,ts_u,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        
        BitRefLevel proble_bit_level = problemPtr->get_BitRefLevel();
        
        SeriesRecorder *recorder_ptr;
        ierr = mField.query_interface(recorder_ptr); CHKERRQ(ierr);
        ierr = recorder_ptr->record_begin(seriesName); CHKERRQ(ierr);
        ierr = recorder_ptr->record_field(seriesName,tempName,proble_bit_level,mask); CHKERRQ(ierr);
        ierr = recorder_ptr->record_end(seriesName); CHKERRQ(ierr);
        
        PetscFunctionReturn(0);
      }
      
    };

    
    /** \brief opearator to caulate tempereature  and rate of temperature at Gauss points
     * \infroup mofem_thermal_elem
     */
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
          }
//          cout << "data.getFieldData() "<<data.getFieldData()<< endl;
//          cout << "data.getN(gg,nb_dofs) "<<data.getN()<< endl;

        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
        
        PetscFunctionReturn(0);
      }
      
    };
    
    /** \brief operator to calculate tempereature at Gauss pts
     * \infroup mofem_Calculate_wt
     */
    struct OpGetTempAtGaussPts: public OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore> {
      OpGetTempAtGaussPts(const string thermal_field_name,CommonData &common_data):
      OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore>(thermal_field_name,common_data.temperatureAtGaussPts) {}
    };
    
    /** \brief operator to calculate moisture concentration at Gauss pts
     * \infroup mofem_Calculate_wt
     */
    struct OpGetConcAtGaussPts: public OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore> {
      OpGetConcAtGaussPts(const string conc_field_name,CommonData &common_data):
      OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore>(conc_field_name,common_data.concentrationAtGaussPts) {}
    };

    
    /** \brief operator to calculate Wt at Gauss pts
     * \infroup mofem_Calculate_wt
     */
    struct OpGetTetWtAtGaussPts: public OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore> {
      OpGetTetWtAtGaussPts(const string field_name,CommonData &common_data):
      OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore>(field_name,common_data.fieldAtGaussPts) {}
    };
   
    
    /** \brief operator to calculate Wt rate at Gauss pts
     * \infroup mofem_Calculate_wt
     */
    struct OpGetTetWtRateAtGaussPts: public OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore> {
      OpGetTetWtRateAtGaussPts(const string rate_name,CommonData &common_data):
      OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore>(rate_name,common_data.fieldRateAtGaussPts) {}
    };
    
    /** \brief operator to calculate Wt rate at Gauss pts
     * \infroup mofem_Calculate_wt
     */
    struct OpGetTetWt_r_BetaAtGaussPts: public OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore> {
      OpGetTetWt_r_BetaAtGaussPts(const string r_beta_name,CommonData &common_data):
      OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore>(r_beta_name,common_data.field_r_BetaAtGaussPts) {}
    };
    
    
    /** \biref operator to calculate right hand side of heat conductivity terms
     * \infroup mofem_thermal_elem
     */
    struct OpGetDegradationRhs: public VolumeElementForcesAndSourcesCore::UserDataOperator {
      
      BlockData &dAta;
      CommonData &commonData;
      bool useTsF;
      OpGetDegradationRhs(const string field_name,BlockData &data,CommonData &common_data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name, UserDataOperator::OPROW),
      dAta(data),commonData(common_data),useTsF(true) {}
      
      ublas::vector<double> Nf;
      
      /** \brief calculate thermal conductivity matrix
       *
       * F = int diffN^T [qw_dot + C*beta*log(1-T/Tg)qw] dOmega^3
       *
       */
      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        
        if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
          PetscFunctionReturn(0);
        }
        
        try {

          if(data.getIndices().size()==0) PetscFunctionReturn(0);
          if(dAta.tEts.find(getMoFEMFEPtr()->get_ent())==dAta.tEts.end()) PetscFunctionReturn(0);

          PetscErrorCode ierr;

          int nb_row_dofs = data.getIndices().size();
          Nf.resize(nb_row_dofs);
          bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));

//          cout << "data.getIndices() "<<data.getIndices()<< endl;
//          cout << "data.getN() "<<data.getN()<< endl;
          
          for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
              double val = getVolume()*getGaussPts()(3,gg);
            if(getHoGaussPtsDetJac().size()>0) {
              val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
            }
            double T=commonData.temperatureAtGaussPts(gg);
            double C=commonData.concentrationAtGaussPts(gg);
            double qw=commonData.fieldAtGaussPts(gg);
            double qw_dot=commonData.fieldRateAtGaussPts(gg);
            
//            cout << "T "<<T<< endl;
//            cout << "C "<<C<< endl;
//            cout << "qw "<<qw<< endl;
//            cout << "qw_dot "<<qw_dot<< endl;
            
            double beta=-0.001682;
            double Tg=126;

//            cout << "dAta.refTemperature "<<dAta.saturation_Conc<< endl;
            C=C/dAta.saturation_Conc;  //relative concentration

            T=T+273.15; Tg=Tg+273.15;  //convert T to Kelvin as the degradation model works for kelvin
            double Nf1=qw_dot+C*beta*log(1-T/Tg)*qw;
            
//            cout << "nb_row_dofs "<<nb_row_dofs<< endl;
//            cout << "data.getN(gg,nb_row_dofs) "<<data.getN(gg,nb_row_dofs)<< endl;
            Nf+=val*Nf1*data.getN(gg,nb_row_dofs);
          }
//          cout << "Nf "<<Nf<< endl;
            ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),
                                &data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
          } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        PetscFunctionReturn(0);
      }
      
    };

    

    
    
    /** \biref operator to calculate left hand side of miisture capacity terms
     * \infroup mofem_thermal_elem
     */
    struct OpGetDegradationLhs: public VolumeElementForcesAndSourcesCore::UserDataOperator {
      
      BlockData &dAta;
      CommonData &commonData;
      OpGetDegradationLhs(const string field_name,BlockData &data,CommonData &common_data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name, UserDataOperator::OPROWCOL),
      dAta(data),commonData(common_data) {}
      
      ublas::matrix<double> M,transM;
      
      /** \brief calculate heat capacity matrix
       *
       * M = int N^T  N dOmega
       *
       */
      PetscErrorCode doWork(
                            int row_side,int col_side,
                            EntityType row_type,EntityType col_type,
                            DataForcesAndSurcesCore::EntData &row_data,
                            DataForcesAndSurcesCore::EntData &col_data) {
        PetscFunctionBegin;
        //        cout<<"OpMoistureCapacityLhs  start "<<endl;
        
        
        try {
          if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
          if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
          
          int nb_row = row_data.getN().size2();
          int nb_col = col_data.getN().size2();
          
          //          cout<<"nb_row  =  "<<nb_row<<endl;
          //          cout<<"nb_col  =  "<<nb_col<<endl;
          //          //  std::string wait;
          //          //  std::cin >> wait;
          
          M.resize(nb_row,nb_col);
          bzero(&*M.data().begin(),nb_row*nb_col*sizeof(double));
          
          for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
            double val = getVolume()*getGaussPts()(3,gg);
            if(getHoGaussPtsDetJac().size()>0) {
              val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
            }

            double T=commonData.temperatureAtGaussPts(gg);
            double C=commonData.concentrationAtGaussPts(gg);
            double beta=-0.001682;
            double Tg=126;
            
            C=C/dAta.saturation_Conc;  //relative concentration
            T=T+273.15; Tg=Tg+273.15;  //convert T to Kelvin as the degradation model works for kelvin
            
            val *=C*beta*log(1-T/Tg)+getFEMethod()->ts_a;

            //cblas
            //double *N_row,*N_col;
            //N_row = &row_data.getN()(gg,0);
            //N_col = &col_data.getN()(gg,0);
            //cblas_dger(CblasRowMajor,
            //  nb_row,nb_col,val,N_row,1,N_col,1,&M(0,0),nb_col);
            
            //ublas
            noalias(M) += val*outer_prod( row_data.getN(gg,nb_row),col_data.getN(gg,nb_col) );
            
          }
          
//          cout<<"M  =  "<<M<<endl;

          PetscErrorCode ierr;
          ierr = MatSetValues(
                              (getFEMethod()->ts_B),
                              nb_row,&row_data.getIndices()[0],
                              nb_col,&col_data.getIndices()[0],
                              &M(0,0),ADD_VALUES); CHKERRQ(ierr);
          if(row_side != col_side || row_type != col_type) {
            transM.resize(nb_col,nb_row);
            noalias(transM) = trans(M);
            ierr = MatSetValues(
                                (getFEMethod()->ts_B),
                                nb_col,&col_data.getIndices()[0],
                                nb_row,&row_data.getIndices()[0],
                                &transM(0,0),ADD_VALUES); CHKERRQ(ierr);
          }
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        
        //        cout<<"OpMoistureCapacityLhs  end "<<endl;
        
        PetscFunctionReturn(0);
      }
      
    };

    
    
    
    PetscErrorCode addWtElement(const string problem_name,const string fe_name,
                                const string field_name,
                                const string rate_name,
                                const string thermal_field_name,
                                const string conc_field_name,
                                const string field_r_beta_name,
                                const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {

      PetscFunctionBegin;
      
      if(mField.check_field(thermal_field_name)) {
        if(mField.check_field(conc_field_name)) {
          
          PetscErrorCode ierr;
          ErrorCode rval;
          
          ierr = mField.add_finite_element(fe_name,MF_ZERO); CHKERRQ(ierr);
          ierr = mField.modify_finite_element_add_field_row(fe_name,field_name); CHKERRQ(ierr);
          ierr = mField.modify_finite_element_add_field_col(fe_name,field_name); CHKERRQ(ierr);
          ierr = mField.modify_finite_element_add_field_data(fe_name,field_name); CHKERRQ(ierr);
          ierr = mField.modify_finite_element_add_field_data(fe_name,rate_name); CHKERRQ(ierr);
          ierr = mField.modify_finite_element_add_field_data(fe_name,thermal_field_name); CHKERRQ(ierr);
          ierr = mField.modify_finite_element_add_field_data(fe_name,conc_field_name); CHKERRQ(ierr);
          ierr = mField.modify_finite_element_add_field_data(fe_name,field_r_beta_name); CHKERRQ(ierr);
          
          if(mField.check_field(mesh_nodals_positions)) {
            ierr = mField.modify_finite_element_add_field_data(fe_name,mesh_nodals_positions); CHKERRQ(ierr);
          }
          ierr = mField.modify_problem_add_finite_element(problem_name,fe_name); CHKERRQ(ierr);
          
          for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
            if(it->get_name().compare(0,17,"MAT_RVE_DIFFUSION") == 0){
              rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
              ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,fe_name); CHKERRQ(ierr);
              //reading saturation concentraiton from the command line (can be input from the CUBIT BLOCKSETS)
              double saturation_conc;
              PetscBool flg;
              ierr = PetscOptionsGetReal(PETSC_NULL,"-my_saturation_conc",&saturation_conc,&flg); CHKERRQ(ierr);
              if(flg == PETSC_TRUE) {
                PetscPrintf(PETSC_COMM_WORLD,"set saturation concentration %3.2f\n",saturation_conc);
                setOfBlocks[it->get_msId()].saturation_Conc = saturation_conc;
              }

            }
          }
        }
      }
      PetscFunctionReturn(0);
    }
    
    
    /** \brief set up operators for unsedy heat flux problem
     * \infroup mofem_thermal_elem
     */
    PetscErrorCode setTimeSteppingProblem(TsCtx &ts_ctx,string fe_name,
                                          string field_name,
                                          string rate_name,
                                          string thermal_field_name,
                                          string conc_field_name,
                                          string field_r_beta_name,
                                          const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;
      
      {
        map<int,BlockData>::iterator sit = setOfBlocks.begin();
        for(;sit!=setOfBlocks.end();sit++) {
//           cout<<"HI before setTimeSteppingProblem() "<<endl;
          //add finite element
          feLhsDegradation.getOpPtrVector().push_back(new OpGetDegradationLhs(field_name,sit->second,commonData));
          
          feRhsDegradation.getOpPtrVector().push_back(new OpGetTempAtGaussPts(thermal_field_name, commonData));
          feRhsDegradation.getOpPtrVector().push_back(new OpGetConcAtGaussPts(conc_field_name,    commonData));
          feRhsDegradation.getOpPtrVector().push_back(new OpGetTetWtAtGaussPts(field_name,commonData));
          feRhsDegradation.getOpPtrVector().push_back(new OpGetTetWtRateAtGaussPts(rate_name,commonData));
          feRhsDegradation.getOpPtrVector().push_back(new OpGetTetWt_r_BetaAtGaussPts(field_r_beta_name,commonData));
          feRhsDegradation.getOpPtrVector().push_back(new OpGetDegradationRhs(field_name, sit->second,  commonData));

//          cout<<"HI After setTimeSteppingProblem() "<<endl;
        }
      }
      {
        bool ho_geometry = false;
        if(mField.check_field(mesh_nodals_positions)) {
          ho_geometry = true;
        }
      }
      
      //rhs
      TsCtx::loops_to_do_type& loops_to_do_Rhs = ts_ctx.get_loops_to_do_IFunction();
      loops_to_do_Rhs.push_back(TsCtx::loop_pair_type(fe_name,&feRhsDegradation));

      //lhs
      TsCtx::loops_to_do_type& loops_to_do_Mat = ts_ctx.get_loops_to_do_IJacobian();
      loops_to_do_Mat.push_back(TsCtx::loop_pair_type(fe_name,&feLhsDegradation));
      
      PetscFunctionReturn(0);
    }
};
  
}

#endif //__CALCULATE_WT_WITH_DERIVATIVE_HPP



