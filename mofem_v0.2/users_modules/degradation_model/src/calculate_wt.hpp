/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
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

#ifndef __CALCULATE_WT_HPP
#define __CALCULATE_WT_HPP

namespace MoFEM {
  
  
  struct Calculate_wt {
    
    struct MyVolumeFE: public TetElementForcesAndSourcesCore {
      MyVolumeFE(FieldInterface &_mField): TetElementForcesAndSourcesCore(_mField) {}
      int getRule(int order) { return order-1; };
    };
    
    
    MyVolumeFE feDegradation;
    MyVolumeFE& getLoopDegradation() { return feDegradation; }

    FieldInterface &mField;
    Calculate_wt(FieldInterface &m_field): feDegradation(m_field), mField(m_field) {}
    
    struct BlockData {
//      double youngModulus;
//      double poissonRatio;
//      double thermalExpansion;
      double saturation_Conc;
      BlockData(): saturation_Conc(100) {}
      Range tEts;
    };
    map<int,BlockData> setOfBlocks;
    
    struct CommonData {
      ublas::vector<double> temperatureAtGaussPts;
      ublas::vector<double> concentrationAtGaussPts;
    };
    CommonData commonData;
    
    
    
  struct LoadTimeSeries: public FEMethod {
    FieldInterface& mField;
    const string tempField;
    const string moisField;
    LoadTimeSeries(FieldInterface& _mField, const string temp_field, const string mois_field):mField(_mField), tempField(temp_field), moisField(mois_field) {}

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      SeriesRecorder *recorder_ptr;
      ierr = mField.query_interface(recorder_ptr); CHKERRQ(ierr);

      ierr = recorder_ptr->load_series_data(tempField,ts_step); CHKERRQ(ierr);
      ierr = recorder_ptr->load_series_data(moisField,ts_step); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }
  };

    
    
    struct OpGetTempAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {
      
      CommonData &commonData;
      int verb;
      OpGetTempAtGaussPts(const string thermal_field_name, CommonData &common_data,int _verb = 0):
      TetElementForcesAndSourcesCore::UserDataOperator(thermal_field_name),
      commonData(common_data),verb(_verb) {}
      
      //loop over EntityType for nodes(1), edges(6), faces(4) and volume(1)
      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        try {
//          cout<<"from OpGetTempConcAtGaussPts ================================ "<<endl;
//          cout<<" data.getFieldData() = " <<data.getFieldData()<<endl;
//          std::string wait;
//          std::cin >> wait;
          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          int nb_dofs = data.getFieldData().size();
          int nb_gauss_pts = data.getN().size1();
          //initialize
          commonData.temperatureAtGaussPts.resize(nb_gauss_pts);
          if(type == MBVERTEX) { //initialize data only once for nodes
            fill(commonData.temperatureAtGaussPts.begin(),commonData.temperatureAtGaussPts.end(),0);
          }
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            commonData.temperatureAtGaussPts[gg] += inner_prod(data.getN(gg,nb_dofs),data.getFieldData());
          }
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        PetscFunctionReturn(0);
      }
      
    };
    
    
    struct OpGetConcAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {
      CommonData &commonData;
      int verb;
      OpGetConcAtGaussPts(const string conc_field_name, CommonData &common_data,int _verb = 0):
      TetElementForcesAndSourcesCore::UserDataOperator(conc_field_name),
      commonData(common_data),verb(_verb) {}
      
      //loop over EntityType for nodes(1), edges(6), faces(4) and volume(1)
      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        try {
//          cout<<"from OpGetConcAtGaussPts ================================ "<<endl;
//          cout<<" data.getFieldData() = " <<data.getFieldData()<<endl;
//          std::string wait;
//          std::cin >> wait;
          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          int nb_dofs = data.getFieldData().size();
          int nb_gauss_pts = data.getN().size1();
          //initialize
          commonData.concentrationAtGaussPts.resize(nb_gauss_pts);
          if(type == MBVERTEX) { //initialize data only once for nodes
            fill(commonData.concentrationAtGaussPts.begin(),commonData.concentrationAtGaussPts.end(),0);
          }
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            commonData.concentrationAtGaussPts[gg] += inner_prod(data.getN(gg,nb_dofs),data.getFieldData());
          }
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        PetscFunctionReturn(0);
      }
      
    };

    
    
    struct OpGetWtAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {
      Vec F;
      BlockData &dAta;
      CommonData &commonData;
      int verb;
      OpGetWtAtGaussPts(const string field_name,BlockData &data,CommonData &common_data,int _verb = 0):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name)
      ,dAta(data),commonData(common_data),verb(_verb) { }
      
      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        
        try {
//          cout<<"from OpGetWtAtGaussPts ================================ "<<endl;
//          cout<<"data.getIndices().size() =  "<<data.getIndices().size()<<endl;

          if(data.getIndices().size()==0) PetscFunctionReturn(0);
          if(dAta.tEts.find(getMoFEMFEPtr()->get_ent())==dAta.tEts.end()) PetscFunctionReturn(0);
          
          PetscErrorCode ierr;

          const FENumeredDofMoFEMEntity *dof_ptr;
          ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
          int rank = dof_ptr->get_max_rank();
          if(rank != 3) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
          }
          
          unsigned int nb_dofs = data.getIndices().size();
          if(nb_dofs % rank != 0) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
          }
          if(data.getN().size2() < nb_dofs/rank) {
            SETERRQ3(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,
                     "data inconsistency N.size2 %d nb_dof %d rank %d",data.getN().size2(),nb_dofs,rank);
          }
          
          if(verb > 0) {
            if(type == MBVERTEX) {
              cout << commonData.temperatureAtGaussPts << endl;
//              cout << "thermal expansion " << dAta.thermalExpansion << endl;
            }
          }
          
          

          for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
            
            //eps_thermal = (phi-phi_ref)*alpha
            //sig_thernal = - (E/1-2mu) * eps_thermal
            //var_eps = [ diff_N[0], diffN[1], diffN[2] ]
            
//            if(dAta.refTemperature != dAta.refTemperature) {
//              SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA ,"invalid data");
//            }
            
            
            double temp = commonData.temperatureAtGaussPts[gg];
            double conc = commonData.concentrationAtGaussPts[gg];
            cout<<"temp "<<temp<<endl;
            cout<<"conc "<<conc<<endl;
            cout<<"dAta.saturation_Conc "<<dAta.saturation_Conc<<endl;
            conc=conc/dAta.saturation_Conc;
            
            
            
            std::string wait;
            std::cin >> wait;

//            double phi = (commonData.temperatureAtGaussPts[gg]-dAta.refTemperature);
//            double eps_thermal = dAta.thermalExpansion*phi;
//            double sig_thermal = -eps_thermal*(dAta.youngModulus/(1.-2*dAta.poissonRatio));
//            
//            double val = sig_thermal*getVolume()*getGaussPts()(3,gg);
//            
//            double *diff_N;
//            diff_N = &data.getDiffN()(gg,0);
//            cblas_daxpy(nb_dofs,val,diff_N,1,&Nf[0],1);
//            
          }
//
//          /*for(unsigned int ii = 0;ii<Nf.size();ii++) {
//           if(Nf[ii] != Nf[ii]) {
//           SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA ,"invalid data");
//           }
//           }*/
//          
//          ierr = VecSetValues(F,data.getIndices().size(),
//                              &data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
//          
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
        
        PetscFunctionReturn(0);
      }
      
    };
    
    
    
    PetscErrorCode addWtElement(const string problem_name,const string fe_name,const string field_name,const string thermal_field_name,
                                           const string conc_field_name, const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;
      
      if(mField.check_field(thermal_field_name)) {
        
        PetscErrorCode ierr;
        ErrorCode rval;
        
        ierr = mField.add_finite_element(fe_name,MF_ZERO); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_row(fe_name,field_name); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_col(fe_name,field_name); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_data(fe_name,field_name); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_data(fe_name,thermal_field_name); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_data(fe_name,conc_field_name); CHKERRQ(ierr);

        if(mField.check_field(mesh_nodals_positions)) {
          ierr = mField.modify_finite_element_add_field_data(fe_name,mesh_nodals_positions); CHKERRQ(ierr);
        }
        ierr = mField.modify_problem_add_finite_element(problem_name,fe_name); CHKERRQ(ierr);
        
        for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_MOISTURESET,it)) {
//          Mat_Elastic mydata;
//          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
//          setOfBlocks[it->get_msId()].youngModulus = mydata.data.Young;
//          setOfBlocks[it->get_msId()].poissonRatio = mydata.data.Poisson;
//          setOfBlocks[it->get_msId()].thermalExpansion = mydata.data.ThermalExpansion;
          rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
          ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,fe_name); CHKERRQ(ierr);
          
          double saturation_conc;
          PetscBool flg;
          ierr = PetscOptionsGetReal(PETSC_NULL,"-my_saturation_conc",&saturation_conc,&flg); CHKERRQ(ierr);
          if(flg == PETSC_TRUE) {
            PetscPrintf(PETSC_COMM_WORLD,"set saturation concentration %3.2f\n",saturation_conc);
            setOfBlocks[it->get_msId()].saturation_Conc = saturation_conc;
          }
        }
        
      }
      
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode setCalculateWtOperators(string field_name, string thermal_field_name,string conc_field_name,Vec &F,int verb = 0) {
      PetscFunctionBegin;
      if(mField.check_field(thermal_field_name)) {
        if(mField.check_field(conc_field_name)) {
          map<int,BlockData>::iterator sit = setOfBlocks.begin();
          for(;sit!=setOfBlocks.end();sit++) {
            cout<<"HI before feDegradation.get_op_to_do_Rhs() "<<endl;
            feDegradation.get_op_to_do_Rhs().push_back(new OpGetTempAtGaussPts(thermal_field_name, commonData,verb));
            feDegradation.get_op_to_do_Rhs().push_back(new OpGetConcAtGaussPts(conc_field_name,    commonData,verb));
            feDegradation.get_op_to_do_Rhs().push_back(new OpGetWtAtGaussPts(field_name,  sit->second,  commonData,verb));

//            feThermalStressRhs.get_op_to_do_Rhs().push_back(new OpThermalStressRhs(field_name,F,sit->second,commonData,verb));

            
            cout<<"HI After feDegradation.get_op_to_do_Rhs() "<<endl;
          }
          
        }
      }
      PetscFunctionReturn(0);
    }
    
    
  };
  
}

#endif //__CALCULATE_WT_HPP



