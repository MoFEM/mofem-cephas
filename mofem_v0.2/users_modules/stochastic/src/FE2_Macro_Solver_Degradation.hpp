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

#ifndef __FE2_MACRO_SOLVER_DEGRADATION_HPP
#define __FE2_MACRO_SOLVER_DEGRADATION_HPP

namespace MoFEM {
  
  /*****************************************************************************
   *                                                                           *
   *                                                                           *
   *                   CALCULATE RVE Dmat AND ITS DERIVATIVES                  *
   *                                                                           *
   *                                                                           *
   ****************************************************************************/
  
  struct Get_Dmat_4_REL {
    
    struct MyVolumeFE: public VolumeElementForcesAndSourcesCore {
      MyVolumeFE(FieldInterface &_mField): VolumeElementForcesAndSourcesCore(_mField) {}
      
      
      int getRule(int order) { return -1; }; //with -1 this function will not work
      
      // This is the same funciton as used in the elastic element to make sure
      //   we use equal number of gauss point for calculation of Dmat and
      //   subsequently in the elastic element.
      
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
    Get_Dmat_4_REL(FieldInterface &m_field):
    feRhs(m_field),
    mField(m_field) {}
    
    struct OpGet_RVEDmat: public VolumeElementForcesAndSourcesCore::UserDataOperator {
      FieldInterface &m_Field;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > &The_Dmat;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > &The_Dmat_r_Em;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > &The_Dmat_r_NUm;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > &The_Dmat_r_Ep;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > &The_Dmat_r_Ez;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > &The_Dmat_r_NUp;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > &The_Dmat_r_NUz;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > &The_Dmat_r_Gzp;
      
      ublas::matrix<double> &Output_Dmat;
      ublas::matrix<double> &Output_Dmat_r_Em;
      ublas::matrix<double> &Output_Dmat_r_NUm;
      ublas::matrix<double> &Output_Dmat_r_Ep;
      ublas::matrix<double> &Output_Dmat_r_Ez;
      ublas::matrix<double> &Output_Dmat_r_NUp;
      ublas::matrix<double> &Output_Dmat_r_NUz;
      ublas::matrix<double> &Output_Dmat_r_Gzp;
      
      OpGet_RVEDmat(FieldInterface &m_Field, const string field_name,
                    map<EntityHandle, ublas::vector<ublas::matrix<double> > > &_The_Dmat,
                    map<EntityHandle, ublas::vector<ublas::matrix<double> > > &_The_Dmat_r_Em,
                    map<EntityHandle, ublas::vector<ublas::matrix<double> > > &_The_Dmat_r_NUm,
                    map<EntityHandle, ublas::vector<ublas::matrix<double> > > &_The_Dmat_r_Ep,
                    map<EntityHandle, ublas::vector<ublas::matrix<double> > > &_The_Dmat_r_Ez,
                    map<EntityHandle, ublas::vector<ublas::matrix<double> > > &_The_Dmat_r_NUp,
                    map<EntityHandle, ublas::vector<ublas::matrix<double> > > &_The_Dmat_r_NUz,
                    map<EntityHandle, ublas::vector<ublas::matrix<double> > > &_The_Dmat_r_Gzp,
                    ublas::matrix<double> &_Output_Dmat,
                    ublas::matrix<double> &_Output_Dmat_r_Em,
                    ublas::matrix<double> &_Output_Dmat_r_NUm,
                    ublas::matrix<double> &_Output_Dmat_r_Ep,
                    ublas::matrix<double> &_Output_Dmat_r_Ez,
                    ublas::matrix<double> &_Output_Dmat_r_NUp,
                    ublas::matrix<double> &_Output_Dmat_r_NUz,
                    ublas::matrix<double> &_Output_Dmat_r_Gzp):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name, UserDataOperator::OPROW),
      m_Field(m_Field),
      The_Dmat(_The_Dmat),The_Dmat_r_Em(_The_Dmat_r_Em),The_Dmat_r_NUm(_The_Dmat_r_NUm),
      The_Dmat_r_Ep(_The_Dmat_r_Ep),The_Dmat_r_Ez(_The_Dmat_r_Ez),
      The_Dmat_r_NUp(_The_Dmat_r_NUp),The_Dmat_r_NUz(_The_Dmat_r_NUz),The_Dmat_r_Gzp(_The_Dmat_r_Gzp),
      Output_Dmat(_Output_Dmat),Output_Dmat_r_Em(_Output_Dmat_r_Em),Output_Dmat_r_NUm(_Output_Dmat_r_NUm),
      Output_Dmat_r_Ep(_Output_Dmat_r_Ep),Output_Dmat_r_Ez(_Output_Dmat_r_Ez),
      Output_Dmat_r_NUp(_Output_Dmat_r_NUp),Output_Dmat_r_NUz(_Output_Dmat_r_NUz),Output_Dmat_r_Gzp(_Output_Dmat_r_Gzp){}
  
      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        
        try {
          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          int nb_gauss_pts = data.getN().size1();
          EntityHandle fe_ent = getMoFEMFEPtr()->get_ent(); // handle of finite element
          // cout<<"\nFrom OpGet_RVEDmat fe_ent "<<fe_ent<<endl;
          
          if(type == MBVERTEX) {
            for(int gg = 0; gg < nb_gauss_pts; gg++) {
              // cout<<"The Dmat at GP "<<gg<<" is ";
              Output_Dmat       = The_Dmat[fe_ent](gg);
              Output_Dmat_r_Em  = The_Dmat_r_Em[fe_ent](gg);
              Output_Dmat_r_NUm = The_Dmat_r_NUm[fe_ent](gg);
              Output_Dmat_r_Ez  = The_Dmat_r_Ez[fe_ent](gg);
              Output_Dmat_r_Ep  = The_Dmat_r_Ep[fe_ent](gg);
              Output_Dmat_r_NUp = The_Dmat_r_NUp[fe_ent](gg);
              Output_Dmat_r_NUz = The_Dmat_r_NUz[fe_ent](gg);
              Output_Dmat_r_Gzp = The_Dmat_r_Gzp[fe_ent](gg);
              // cout<<Output_Dmat<<endl;
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
    
    PetscErrorCode setMacro_DmatRhsOperators(FieldInterface &m_field_Macro,
                                             string field_name,
                                             map<EntityHandle, ublas::vector<ublas::matrix<double> > > &Dmat_RVE,
                                             map<EntityHandle, ublas::vector<ublas::matrix<double> > > &Dmat_RVE_r_Em,
                                             map<EntityHandle, ublas::vector<ublas::matrix<double> > > &Dmat_RVE_r_NUm,
                                             map<EntityHandle, ublas::vector<ublas::matrix<double> > > &Dmat_RVE_r_Ep,
                                             map<EntityHandle, ublas::vector<ublas::matrix<double> > > &Dmat_RVE_r_Ez,
                                             map<EntityHandle, ublas::vector<ublas::matrix<double> > > &Dmat_RVE_r_NUp,
                                             map<EntityHandle, ublas::vector<ublas::matrix<double> > > &Dmat_RVE_r_NUz,
                                             map<EntityHandle, ublas::vector<ublas::matrix<double> > > &Dmat_RVE_r_Gzp,
                                             ublas::matrix<double> &Output_Dmat,
                                             ublas::matrix<double> &Output_Dmat_r_Em,
                                             ublas::matrix<double> &Output_Dmat_r_NUm,
                                             ublas::matrix<double> &Output_Dmat_r_Ep,
                                             ublas::matrix<double> &Output_Dmat_r_Ez,
                                             ublas::matrix<double> &Output_Dmat_r_NUp,
                                             ublas::matrix<double> &Output_Dmat_r_NUz,
                                             ublas::matrix<double> &Output_Dmat_r_Gzp) {
      PetscFunctionBegin;
      
      feRhs.getOpPtrVector().push_back(new OpGet_RVEDmat(m_field_Macro,field_name,
                                                         Dmat_RVE,Dmat_RVE_r_Em,Dmat_RVE_r_NUm,
                                                         Dmat_RVE_r_Ep,Dmat_RVE_r_Ez,
                                                         Dmat_RVE_r_NUp,Dmat_RVE_r_NUz,Dmat_RVE_r_Gzp,
                                                         Output_Dmat,Output_Dmat_r_Em,Output_Dmat_r_NUm,
                                                         Output_Dmat_r_Ep,Output_Dmat_r_Ez,
                                                         Output_Dmat_r_NUp,Output_Dmat_r_NUz,Output_Dmat_r_Gzp));
      
      PetscFunctionReturn(0);
    }
    
    
  };
  
  struct FE2_RVE_Dmat_Degradation_Disp {
    
    struct MyVolumeFE: public VolumeElementForcesAndSourcesCore {
      MyVolumeFE(FieldInterface &_mField): VolumeElementForcesAndSourcesCore(_mField) {}
      
      
      int getRule(int order) { return -1; }; //with -1 this function will not work
      
      // This is the same funciton as used in the elastic element to make sure
      //   we use equal number of gauss point for calculation of Dmat and
      //   subsequently in the elastic element.
      
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
    FE2_RVE_Dmat_Degradation_Disp(FieldInterface &m_field):
    feRhs(m_field),
    mField(m_field) {}
    
    
    struct CommonData {
      ublas::vector<double> wtAtGaussPts;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > Dmat_RVE;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > Dmat_RVE_r_Em;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > Dmat_RVE_r_NUm;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > Dmat_RVE_r_Ep;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > Dmat_RVE_r_Ez;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > Dmat_RVE_r_NUp;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > Dmat_RVE_r_NUz;
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > Dmat_RVE_r_Gzp;
    };
    CommonData commonData;
    
    
    
    template<typename OP>
    struct OpGetFieldAtGaussPts: public OP::UserDataOperator {
      
      ublas::vector<double> &fieldAtGaussPts;
      OpGetFieldAtGaussPts(const string field_name,ublas::vector<double> &field_at_gauss_pts):
      OP::UserDataOperator(field_name,OP::UserDataOperator::OPROW),
      fieldAtGaussPts(field_at_gauss_pts) {}
      
      /** \brief operator calculating temperature and rate of temperature
       *
       *  temperature or rate of temperature is calculated by multiplying shape
       *  functions by degrees of freedom
       */
      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        try {
          // cout<<"form OpGetFieldAtGaussPts "<<endl;
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
          
          for(int gg = 0; gg < nb_gauss_pts; gg++) {
            fieldAtGaussPts[gg] += inner_prod(data.getN(gg,nb_dofs),data.getFieldData());
            // cout<<"fieldAtGaussPts[gg] "<<fieldAtGaussPts[gg] <<endl;
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
    
    
    struct OpCalculate_RVEDmat: public VolumeElementForcesAndSourcesCore::UserDataOperator {
      
      PetscErrorCode ierr;
      
      FieldInterface &m_field_RVE;
      int nvars;
      int nders;
      int num_rvars;
      ublas::vector<double> matprop;
      vector<string> stochastic_fields;
      vector<string> vars_name;
      CommonData &commonData;
      
      OpCalculate_RVEDmat(FieldInterface &m_field_RVE, const string field_name,
                          CommonData &common_data, int _nvars, int _nders,
                          vector<string> _stochastic_fields,
                          ublas::vector<double> _matprop,
                          int _num_rvars,
                          vector<string> _vars_name):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name, UserDataOperator::OPROW),
      m_field_RVE(m_field_RVE),
      commonData(common_data),
      nvars(_nvars),
      nders(_nders),
      stochastic_fields(_stochastic_fields),
      matprop(_matprop),
      num_rvars(_num_rvars),
      vars_name(_vars_name){}
      
      ~OpCalculate_RVEDmat(){}
      
      
      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        try {
          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          int nb_gauss_pts = data.getN().size1();
          EntityHandle fe_ent = getMoFEMFEPtr()->get_ent(); //handle of finite element
          cout<<"\nElement: "<<fe_ent<<endl;
          commonData.Dmat_RVE[fe_ent].resize(nb_gauss_pts);
          commonData.Dmat_RVE_r_Em[fe_ent].resize(nb_gauss_pts);
          commonData.Dmat_RVE_r_NUm[fe_ent].resize(nb_gauss_pts);
          commonData.Dmat_RVE_r_Ep[fe_ent].resize(nb_gauss_pts);
          commonData.Dmat_RVE_r_Ez[fe_ent].resize(nb_gauss_pts);
          commonData.Dmat_RVE_r_NUp[fe_ent].resize(nb_gauss_pts);
          commonData.Dmat_RVE_r_NUz[fe_ent].resize(nb_gauss_pts);
          commonData.Dmat_RVE_r_Gzp[fe_ent].resize(nb_gauss_pts);
          // map<EntityHandle, ublas::vector<ublas::matrix<double> > > Dmat_RVE;
          
          if(type == MBVERTEX) {  //the doWork loop is 1+6+4+1 times but we want to loop over Guass points only onece
            ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
            int field_rank=3; // it is mechanical problem
            applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
            
            //cout<<"nb_gauss_pts "<<nb_gauss_pts<<endl;
            
            for(int gg = 0; gg < nb_gauss_pts; gg++) {
              Vec F1,F2,F3,F4,F5,F6;
              Vec D1;
              Vec dF1,dF2,dF3,dF4,dF5,dF6;
              Vec dD1;
              Vec ddF1,ddF2,ddF3,ddF4,ddF5,ddF6;
              Vec ddD1;
              Mat A;
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F2); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F3); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F4); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F5); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F6); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF2); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF3); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF4); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF5); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF6); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF2); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF3); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF4); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF5); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF6); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&D1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD1); CHKERRABORT(PETSC_COMM_WORLD, ierr);
              
              ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD1); CHKERRABORT(PETSC_COMM_WORLD, ierr);
              
              ierr = m_field_RVE.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_RVE",&A); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              
              /*****************************************************************
               *
               *  2. Get the volume of RVE
               *
               ****************************************************************/
              double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
              Vec RVE_volume_Vec;
              ParallelComm* pcomm_RVE = ParallelComm::get_pcomm(&m_field_RVE.get_moab(),MYPCOMM_INDEX);
              //            cout<<" pcomm_RVE->size() = "<<pcomm_RVE->size()<<endl;
              ierr = VecCreateMPI(PETSC_COMM_SELF, 1, pcomm_RVE->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
              ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
              RVEVolume MyRVEVol(m_field_RVE,A,D1,F1,0.0,0.0, RVE_volume_Vec);
              
              /*****************************************************************
               *
               * 3. Update values of material properties
               *
               ****************************************************************/
              
              for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)) {
                //cout << endl << *it << endl;
                
                //Get block name
                string name = it->get_name();
                // -------------------------------------
                //
                // 3.1 Modify matrix material properties
                //
                // -------------------------------------
                if (name.compare(0,13,"MAT_ELASTIC_1") == 0) {
                  Mat_Elastic mydata;
                  ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
                  
                  string ParameterName;
                   for (int ii = 1;ii <= num_rvars;ii++) {
                   ParameterName = vars_name[ii];
                   
                   if (ParameterName.compare(0,2,"Em") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
                   mydata.data.Young = matprop(ii-1);
                   }
                   else if (ParameterName.compare(0,3,"NUm") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
                   mydata.data.Poisson = matprop(ii-1);
                   }
                   ParameterName.clear();
                   }
                   
                   // YoungModulus = mydata.data.Young;
                   // PoissonRatio = mydata.data.Poisson;
                   ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
                  //cout <<"\n\nMatrix material:" << mydata;
                  //cout<<"This is matrix material properties"<<endl;
                }
                
                // ------------------------------------
                //
                // 3.2 Modify fibre material properties
                //
                // ------------------------------------
                if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
                  Mat_Elastic_TransIso mydata;
                  ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
                  
                  string ParameterName;//cout<<"\nthe size of num_vars is "<<vars_name.size()<<"\t"<<num_rvars<<endl;
                  for (int ii = 1; ii <= num_rvars; ii++) {
                    ParameterName = vars_name[ii];
                    //cout<<"the variable name is "<<vars_name[ii]<<endl;
                    if (ParameterName.compare(0,2,"Ez") == 0) {
                      mydata.data.Youngz = matprop(ii-1);
                    }
                    else if (ParameterName.compare(0,2,"Ep") == 0) {
                      mydata.data.Youngp = matprop(ii-1);
                    }
                    else if (ParameterName.compare(0,3,"NUp") == 0) {
                      mydata.data.Poissonp = matprop(ii-1);
                    }
                    else if (ParameterName.compare(0,3,"NUz") == 0) {
                      mydata.data.Poissonpz = matprop(ii-1);
                    }
                    else if (ParameterName.compare(0,3,"Gzp") == 0) {
                      mydata.data.Shearzp = matprop(ii-1);
                    }
                    else if (ParameterName.compare(0,2,"Ef") == 0) {
                      mydata.data.Youngz = matprop(ii-1);
                      mydata.data.Youngp = matprop(ii-1);
                    }
                    else if (ParameterName.compare(0,3,"NUf") == 0) {
                      mydata.data.Poissonp  = matprop(ii-1);
                      mydata.data.Poissonpz = matprop(ii-1);
                    }
                    else if ((ParameterName.compare(0,3,"NUf") == 0) || (ParameterName.compare(0,2,"Ef") == 0)) {
                      mydata.data.Shearzp = mydata.data.Youngz/(2*(1+mydata.data.Poissonp));
                    }
                    ParameterName.clear();
                  }
                   
                  ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
                  //cout<<"\n\nFibre material:" << mydata;
                  //cout<<"This is fibre material properties"<<endl;
                }
              }
              
              //=============================================================================================================
              
              //              cout<<"gg Start =  "<<gg <<endl;
              //We don't need to calculate internal forces for RVE, as ElasticFEMethod is used to assemble A matirx only
              //so noo need to create MyElasticFEMethod here
              ElasticFEMethod_Matrix my_fe_marix(m_field_RVE,A,D1,F1,0.0,0.0,commonData.wtAtGaussPts(gg),"DISP_RVE");cout<<"The wt is "<<commonData.wtAtGaussPts(gg)<<endl;
              TranIsotropicFibreDirRotElasticFEMethod my_fe_transiso(m_field_RVE,A,D1,F1,"DISP_RVE");
              ElasticFE_RVELagrange_Disp_Multi_Rhs MyFE_RVELagrange(m_field_RVE,A,D1,F1,F2,F3,F4,F5,F6,applied_strain,
                                                                    "DISP_RVE","Lagrange_mul_disp",field_rank);
              
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
                ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
                ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
                ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
              }
              
              
              /*************************************************************************
               *
               *  3. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
               *     [K][U] = [F]
               *
               ************************************************************************/
              //Solver
              KSP solver_RVE;
              ierr = KSPCreate(PETSC_COMM_SELF,&solver_RVE); CHKERRQ(ierr);
              ierr = KSPSetOperators(solver_RVE,A,A); CHKERRQ(ierr);
              ierr = KSPSetFromOptions(solver_RVE); CHKERRQ(ierr);
              ierr = KSPSetUp(solver_RVE); CHKERRQ(ierr);
              //ierr = VecView(D1,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
              
              ublas::matrix<FieldData> Dmat;
              Dmat.resize(6,6); Dmat.clear();
              
              ublas::matrix<FieldData> Dmat_r_Em, Dmat_r_NUm;
              ublas::matrix<FieldData> Dmat_r_Ep, Dmat_r_Ez, Dmat_r_NUp, Dmat_r_NUz, Dmat_r_Gzp;
              Dmat_r_Em.resize(6,6);    Dmat_r_Em.clear();
              Dmat_r_NUm.resize(6,6);   Dmat_r_NUm.clear();
              
              Dmat_r_Ep.resize(6,6);    Dmat_r_Ep.clear();
              Dmat_r_Ez.resize(6,6);    Dmat_r_Ez.clear();
              Dmat_r_NUp.resize(6,6);   Dmat_r_NUp.clear();
              Dmat_r_NUz.resize(6,6);   Dmat_r_NUz.clear();
              Dmat_r_Gzp.resize(6,6);   Dmat_r_Gzp.clear();
              
              //create a vector for 6 components of homogenized stress
              Vec Stress_Homo, Stress_Homo_r;
              if(pcomm_RVE->rank()==0) {
                VecCreateGhost(PETSC_COMM_SELF,6,6,0,PETSC_NULL,&Stress_Homo);
                VecCreateGhost(PETSC_COMM_SELF,6,6,0,PETSC_NULL,&Stress_Homo_r);
              } else {
                int ghost[] = {0,1,2,3,4,5};
                VecCreateGhost(PETSC_COMM_SELF,0,6,6,ghost,&Stress_Homo);
                VecCreateGhost(PETSC_COMM_SELF,0,6,6,ghost,&Stress_Homo_r);
              }
              
              //----------------------------------------------------------------
              // 3.1 Solving the equation to get nodal displacement
              //     case 1: applied macro strain: [1 0 0 0 0 0]^T
              //----------------------------------------------------------------
              /*cout<<"=======================================================\n";
              cout<<"        Applied strain [1 0 0 0 0 0]^T\n";
              cout<<"=======================================================\n";*/
              
              //-------------
              // ZEROTH-ORDER
              //-------------
              ierr = KSPSolve(solver_RVE,F1,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              //ierr = VecView(D1,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
              
              ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
              
              ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_1(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
              
              // calculate homogenized stress
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_1);  CHKERRQ(ierr);
              PetscScalar *avec;
              VecGetArray(Stress_Homo, &avec);
              for(int ii=0; ii<6; ii++){
                Dmat(ii,0)=*avec;
                avec++;
              }
              VecRestoreArray(Stress_Homo,&avec);
              
              //---------------------
              // FIRST & SECOND-ORDER
              //---------------------
              for(int ii = 1; ii <= num_rvars; ii++) {
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
                else if (VariableName.compare(0,3,"NUz") == 0) {
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
                else if (VariableName.compare(0,2,"Ef") == 0) {
                  idx_disp = 7;//cout<<"\n\nThe random variable is Ef"<<endl;
                }
                else if (VariableName.compare(0,3,"NUf") == 0) {
                  idx_disp = 8;//cout<<"\n\nThe random variable is NUf"<<endl;
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
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Em(m_field_RVE,A,D1,dF1,"DISP_RVE","Young","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
                }
                else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUm(m_field_RVE,A,D1,dF1,"DISP_RVE","Poisson","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
                }
                else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUp(m_field_RVE,A,D1,dF1,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
                }
                else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUpz(m_field_RVE,A,D1,dF1,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
                }
                else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ep(m_field_RVE,A,D1,dF1,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
                }
                else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ez(m_field_RVE,A,D1,dF1,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
                }
                else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Gzp(m_field_RVE,A,D1,dF1,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
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
                  
                  ierr = KSPSolve(solver_RVE,dF1,dD1); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  //cout<<"\n\nFirst order case"<<endl;
                }
                //------------
                // b. Calculating first-order homogenized stress
                //------------
                if ((idx_disp>=0) && (idx_disp<nders)) {
                  ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
                  
                  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r1(m_field_RVE,A,dD1,dF1,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r1);  CHKERRQ(ierr);
                  
                  if(pcomm_RVE->rank()==0) {
                    PetscScalar    *avec_r;
                    VecGetArray(Stress_Homo_r, &avec_r);
                    
                    //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
                    //cout<< "\n"<<ss_field<<" = \n\n";
                    for(int irow = 0; irow < 6; irow++){
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
                          Dmat_r_NUz(irow,0) = *avec_r; break;
                        case 4:
                          Dmat_r_Ep(irow,0) = *avec_r; break;
                        case 5:
                          Dmat_r_Ez(irow,0) = *avec_r; break;
                        case 6:
                          Dmat_r_Gzp(irow,0) = *avec_r; break;
                          /*case 7:
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
                           Dmat_rs_GzpGzp(irow,0) = *avec_r; break;*/
                      }
                      // write result to output file
                      //TheFile<<setprecision(15)<<*avec_r<<'\n';
                      avec_r++;
                    }
                    VecRestoreArray(Stress_Homo_r, &avec_r);
                  }
                }
              }
              
              //----------------------------------------------------------------
              // 3.1 Solving the equation to get nodal displacement
              //     case 2: applied macro strain: [0 1 0 0 0 0]^T
              //----------------------------------------------------------------
              /*cout<<"=======================================================\n";
              cout<<"        Applied strain [0 1 0 0 0 0]^T\n";
              cout<<"=======================================================\n";*/
              
              //-------------
              // ZEROTH-ORDER
              //-------------
              ierr = KSPSolve(solver_RVE,F2,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              
              ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
              //    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
              ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_2(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_2);  CHKERRQ(ierr);
              
              VecGetArray(Stress_Homo, &avec);
              for(int ii = 0; ii < 6; ii++) {
                Dmat(ii,1) = *avec;
                avec++;
              }
              VecRestoreArray(Stress_Homo,&avec);
              
              //---------------------
              // FIRST & SECOND-ORDER
              //---------------------
              for(int ii = 1; ii <= num_rvars; ii++) {
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
                else if (VariableName.compare(0,3,"NUz") == 0) {
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
                else if (VariableName.compare(0,2,"Ef") == 0) {
                  idx_disp = 7;
                }
                else if (VariableName.compare(0,3,"NUf") == 0) {
                  idx_disp = 8;
                }
                VariableName.clear();
                
                if ((idx_disp >=0) && (idx_disp < nvars)) {
                  ierr = VecZeroEntries(dF2); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                } else if ((idx_disp >= nvars) && (idx_disp < nders)) {
                  ierr = VecZeroEntries(ddF2); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                }
                
                if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Em(m_field_RVE,A,D1,dF2,"DISP_RVE","Young","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
                }
                else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUm(m_field_RVE,A,D1,dF2,"DISP_RVE","Poisson","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
                }
                else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUp(m_field_RVE,A,D1,dF2,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
                }
                else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUpz(m_field_RVE,A,D1,dF2,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
                }
                else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ep(m_field_RVE,A,D1,dF2,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
                }
                else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ez(m_field_RVE,A,D1,dF2,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
                }
                else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Gzp(m_field_RVE,A,D1,dF2,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
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
                  
                  ierr = KSPSolve(solver_RVE,dF2,dD1); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  //cout<<"\n\nFirst order case"<<endl;
                }
                //------------
                // b. Calculating first-order homogenized stress
                //------------
                if ((idx_disp>=0) && (idx_disp<nders)) {
                  ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
                  
                  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r2(m_field_RVE,A,dD1,dF2,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r2);  CHKERRQ(ierr);
                  
                  if(pcomm_RVE->rank()==0) {
                    PetscScalar    *avec_r;
                    VecGetArray(Stress_Homo_r, &avec_r);
                    
                    //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
                    //cout<< "\n"<<ss_field<<" = \n\n";
                    for(int irow = 0; irow < 6; irow++){
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
                          Dmat_r_NUz(irow,1) = *avec_r; break;
                        case 4:
                          Dmat_r_Ep(irow,1) = *avec_r; break;
                        case 5:
                          Dmat_r_Ez(irow,1) = *avec_r; break;
                        case 6:
                          Dmat_r_Gzp(irow,1) = *avec_r; break;
                          /*case 7:
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
                           Dmat_rs_GzpGzp(irow,0) = *avec_r; break;*/
                      }
                      avec_r++;
                    }
                    VecRestoreArray(Stress_Homo_r, &avec_r);
                  }
                }
              }
              
              //----------------------------------------------------------------
              // 3.1 Solving the equation to get nodal displacement
              //     case 3: applied macro strain: [0 0 1 0 0 0]^T
              //----------------------------------------------------------------
              /*cout<<"=======================================================\n";
              cout<<"        Applied strain [0 0 1 0 0 0]^T\n";
              cout<<"=======================================================\n";*/
              
              //-------------
              // ZEROTH-ORDER
              //-------------
              ierr = KSPSolve(solver_RVE,F3,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
              
              ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_3(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_3);  CHKERRQ(ierr);
              
              VecGetArray(Stress_Homo, &avec);
              for(int ii = 0; ii < 6; ii++){
                Dmat(ii,2)=*avec;
                avec++;
              }
              VecRestoreArray(Stress_Homo, &avec);
              
              
              //---------------------
              // FIRST & SECOND-ORDER
              //---------------------
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
                else if (VariableName.compare(0,3,"NUz") == 0) {
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
                else if (VariableName.compare(0,2,"Ef") == 0) {
                  idx_disp = 7;
                }
                else if (VariableName.compare(0,3,"NUf") == 0) {
                  idx_disp = 8;
                }
                VariableName.clear();
                
                if ((idx_disp >=0) && (idx_disp < nvars)) {
                  ierr = VecZeroEntries(dF3); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                } else if ((idx_disp >= nvars) && (idx_disp < nders)) {
                  ierr = VecZeroEntries(ddF3); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                }
                if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Em(m_field_RVE,A,D1,dF3,"DISP_RVE","Young","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
                }
                else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUm(m_field_RVE,A,D1,dF3,"DISP_RVE","Poisson","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
                }
                else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUp(m_field_RVE,A,D1,dF3,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
                }
                else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUpz(m_field_RVE,A,D1,dF3,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
                }
                else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ep(m_field_RVE,A,D1,dF3,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
                }
                else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ez(m_field_RVE,A,D1,dF3,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
                }
                else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Gzp(m_field_RVE,A,D1,dF3,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
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
                  
                  ierr = KSPSolve(solver_RVE,dF3,dD1); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  //cout<<"\n\nFirst order case"<<endl;
                }
                //------------
                // b. Calculating first-order homogenized stress
                //------------
                if ((idx_disp>=0) && (idx_disp<nders)) {
                  ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
                  
                  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r3(m_field_RVE,A,dD1,dF3,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r3);  CHKERRQ(ierr);
                  
                  if(pcomm_RVE->rank()==0) {
                    PetscScalar    *avec_r;
                    VecGetArray(Stress_Homo_r, &avec_r);
                    
                    //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
                    //cout<< "\n"<<ss_field<<" = \n\n";
                    for(int irow = 0; irow < 6; irow++){
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
                          Dmat_r_NUz(irow,2) = *avec_r; break;
                        case 4:
                          Dmat_r_Ep(irow,2) = *avec_r; break;
                        case 5:
                          Dmat_r_Ez(irow,2) = *avec_r; break;
                        case 6:
                          Dmat_r_Gzp(irow,2) = *avec_r; break;
                          /*case 7:
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
                           Dmat_rs_GzpGzp(irow,0) = *avec_r; break;*/
                      }
                      avec_r++;
                    }
                    VecRestoreArray(Stress_Homo_r, &avec_r);
                  }
                }
              }
              
              //----------------------------------------------------------------
              // 3.1 Solving the equation to get nodal displacement
              //     case 4: applied macro strain: [0 0 0 1 0 0]^T
              //----------------------------------------------------------------
              /*cout<<"=======================================================\n";
              cout<<"        Applied strain [0 0 0 1 0 0]^T\n";
              cout<<"=======================================================\n";*/
              
              //-------------
              // ZEROTH-ORDER
              //-------------
              ierr = KSPSolve(solver_RVE,F4,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
              ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_4(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_4);  CHKERRQ(ierr);
              
              VecGetArray(Stress_Homo, &avec);
              for(int ii = 0; ii < 6; ii++){
                Dmat(ii,3)=*avec;
                avec++;
              }
              VecRestoreArray(Stress_Homo, &avec);
              
              
              //---------------------
              // FIRST & SECOND-ORDER
              //---------------------
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
                else if (VariableName.compare(0,3,"NUz") == 0) {
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
                else if (VariableName.compare(0,2,"Ef") == 0) {
                  idx_disp = 7;
                }
                else if (VariableName.compare(0,3,"NUf") == 0) {
                  idx_disp = 8;
                }
                VariableName.clear();
                
                if ((idx_disp >=0) && (idx_disp < nvars)) {
                  ierr = VecZeroEntries(dF4); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                } else if ((idx_disp >= nvars) && (idx_disp < nders)) {
                  ierr = VecZeroEntries(ddF4); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                }
                if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Em(m_field_RVE,A,D1,dF4,"DISP_RVE","Young","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
                }
                else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUm(m_field_RVE,A,D1,dF4,"DISP_RVE","Poisson","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
                }
                else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUp(m_field_RVE,A,D1,dF4,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
                }
                else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUpz(m_field_RVE,A,D1,dF4,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
                }
                else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ep(m_field_RVE,A,D1,dF4,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
                }
                else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ez(m_field_RVE,A,D1,dF4,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
                }
                else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Gzp(m_field_RVE,A,D1,dF4,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
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
                  
                  ierr = KSPSolve(solver_RVE,dF4,dD1); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  //cout<<"\n\nFirst order case"<<endl;
                }
                //------------
                // b. Calculating first-order homogenized stress
                //------------
                if ((idx_disp>=0) && (idx_disp<nders)) {
                  ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
                  
                  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r4(m_field_RVE,A,dD1,dF4,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r4);  CHKERRQ(ierr);
                  
                  if(pcomm_RVE->rank()==0) {
                    PetscScalar    *avec_r;
                    VecGetArray(Stress_Homo_r, &avec_r);
                    
                    //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
                    //cout<< "\n"<<ss_field<<" = \n\n";
                    for(int irow = 0; irow < 6; irow++){
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
                          Dmat_r_NUz(irow,3) = *avec_r; break;
                        case 4:
                          Dmat_r_Ep(irow,3) = *avec_r; break;
                        case 5:
                          Dmat_r_Ez(irow,3) = *avec_r; break;
                        case 6:
                          Dmat_r_Gzp(irow,3) = *avec_r; break;
                          /*case 7:
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
                           Dmat_rs_GzpGzp(irow,0) = *avec_r; break;*/
                      }
                      avec_r++;
                    }
                    VecRestoreArray(Stress_Homo_r, &avec_r);
                  }
                }
              }
              
              //----------------------------------------------------------------
              // 3.1 Solving the equation to get nodal displacement
              //     case 5: applied macro strain: [0 0 0 0 1 0]^T
              //----------------------------------------------------------------
              /*cout<<"=======================================================\n";
              cout<<"        Applied strain [0 0 0 0 1 0]^T\n";
              cout<<"=======================================================\n";*/
              
              //-------------
              // ZEROTH-ORDER
              //-------------
              ierr = KSPSolve(solver_RVE,F5,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              
              ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
              
              ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_5(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_5);  CHKERRQ(ierr);
              
              VecGetArray(Stress_Homo, &avec);
              for(int ii = 0; ii < 6; ii++){
                Dmat(ii,4) = *avec;
                avec++;
              }
              VecRestoreArray(Stress_Homo, &avec);
              
              
              //---------------------
              // FIRST & SECOND-ORDER
              //---------------------
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
                else if (VariableName.compare(0,3,"NUz") == 0) {
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
                else if (VariableName.compare(0,2,"Ef") == 0) {
                  idx_disp = 7;
                }
                else if (VariableName.compare(0,3,"NUf") == 0) {
                  idx_disp = 8;
                }
                VariableName.clear();
                
                if ((idx_disp >=0) && (idx_disp < nvars)) {
                  ierr = VecZeroEntries(dF5); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                } else if ((idx_disp >= nvars) && (idx_disp < nders)) {
                  ierr = VecZeroEntries(ddF5); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                }
                if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Em(m_field_RVE,A,D1,dF5,"DISP_RVE","Young","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
                }
                else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUm(m_field_RVE,A,D1,dF5,"DISP_RVE","Poisson","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
                }
                else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUp(m_field_RVE,A,D1,dF5,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
                }
                else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUpz(m_field_RVE,A,D1,dF5,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
                }
                else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ep(m_field_RVE,A,D1,dF5,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
                }
                else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ez(m_field_RVE,A,D1,dF5,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
                }
                else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Gzp(m_field_RVE,A,D1,dF5,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
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
                  
                  ierr = KSPSolve(solver_RVE,dF5,dD1); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  //cout<<"\n\nFirst order case"<<endl;
                }
                //------------
                // b. Calculating first-order homogenized stress
                //------------
                if ((idx_disp>=0) && (idx_disp<nders)) {
                  ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
                  
                  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r5(m_field_RVE,A,dD1,dF5,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r5);  CHKERRQ(ierr);
                  
                  if(pcomm_RVE->rank()==0) {
                    PetscScalar    *avec_r;
                    VecGetArray(Stress_Homo_r, &avec_r);
                    
                    //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
                    //cout<< "\n"<<ss_field<<" = \n\n";
                    for(int irow = 0; irow < 6; irow++){
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
                          Dmat_r_NUz(irow,4) = *avec_r; break;
                        case 4:
                          Dmat_r_Ep(irow,4) = *avec_r; break;
                        case 5:
                          Dmat_r_Ez(irow,4) = *avec_r; break;
                        case 6:
                          Dmat_r_Gzp(irow,4) = *avec_r; break;
                          /*case 7:
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
                           Dmat_rs_GzpGzp(irow,0) = *avec_r; break;*/
                      }
                      avec_r++;
                    }
                    VecRestoreArray(Stress_Homo_r, &avec_r);
                  }
                }
              }
              
              
              //----------------------------------------------------------------
              // 3.1 Solving the equation to get nodal displacement
              //     case 6: applied macro strain: [0 0 0 0 0 1]^T
              //----------------------------------------------------------------
              /*cout<<"=======================================================\n";
              cout<<"         Applied strain [0 0 0 0 0 1]^T\n";
              cout<<"=======================================================\n";*/
              
              //-------------
              // ZEROTH-ORDER
              //-------------
              ierr = KSPSolve(solver_RVE,F6,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              
              ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
              ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_6(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_6);  CHKERRQ(ierr);
              
              VecGetArray(Stress_Homo, &avec);
              for(int ii = 0; ii < 6; ii++){
                Dmat(ii,5)=*avec;
                avec++;
              }
              VecRestoreArray(Stress_Homo, &avec);
              
              
              //---------------------
              // FIRST & SECOND-ORDER
              //---------------------
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
                else if (VariableName.compare(0,3,"NUz") == 0) {
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
                else if (VariableName.compare(0,2,"Ef") == 0) {
                  idx_disp = 7;
                }
                else if (VariableName.compare(0,3,"NUf") == 0) {
                  idx_disp = 8;
                }
                VariableName.clear();
                
                if ((idx_disp >=0) && (idx_disp < nvars)) {
                  ierr = VecZeroEntries(dF6); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                } else if ((idx_disp >= nvars) && (idx_disp < nders)) {
                  ierr = VecZeroEntries(ddF6); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                }
                if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Em(m_field_RVE,A,D1,dF6,"DISP_RVE","Young","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
                }
                else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUm(m_field_RVE,A,D1,dF6,"DISP_RVE","Poisson","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
                }
                else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUp(m_field_RVE,A,D1,dF6,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
                }
                else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUpz(m_field_RVE,A,D1,dF6,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
                }
                else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ep(m_field_RVE,A,D1,dF6,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
                }
                else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ez(m_field_RVE,A,D1,dF6,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
                }
                else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Gzp(m_field_RVE,A,D1,dF6,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
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
                  
                  ierr = KSPSolve(solver_RVE,dF6,dD1); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  //cout<<"\n\nFirst order case"<<endl;
                }
                //------------
                // b. Calculating first-order homogenized stress
                //------------
                if ((idx_disp>=0) && (idx_disp<nders)) {
                  ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
                  
                  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r6(m_field_RVE,A,dD1,dF6,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r6);  CHKERRQ(ierr);
                  
                  if(pcomm_RVE->rank()==0) {
                    PetscScalar    *avec_r;
                    VecGetArray(Stress_Homo_r, &avec_r);
                    
                    //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
                    //cout<< "\n"<<ss_field<<" = \n\n";
                    for(int irow = 0; irow < 6; irow++){
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
                          Dmat_r_NUz(irow,5) = *avec_r; break;
                        case 4:
                          Dmat_r_Ep(irow,5) = *avec_r; break;
                        case 5:
                          Dmat_r_Ez(irow,5) = *avec_r; break;
                        case 6:
                          Dmat_r_Gzp(irow,5) = *avec_r; break;
                          /*case 7:
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
                           Dmat_rs_GzpGzp(irow,0) = *avec_r; break;*/
                      }
                      avec_r++;
                    }
                    VecRestoreArray(Stress_Homo_r, &avec_r);
                  }
                }
              }
              
              // -
              
              //cout<<"fe_ent "<<fe_ent << "gg  "<<gg<<endl;
              //cout<<"\nDmat = "<<Dmat<<endl;
              /*cout<<"\nDmat_r_Em   = "<<Dmat_r_Em<<endl;
               cout<<"\nDmat_r_NUm  = "<<Dmat_r_NUm<<endl;
               cout<<"\nDmat_r_Ep   = "<<Dmat_r_Ep<<endl;
               cout<<"\nDmat_r_Ez   = "<<Dmat_r_Ez<<endl;
               cout<<"\nDmat_r_Gzp  = "<<Dmat_r_Gzp<<endl;
               cout<<"\nDmat_r_NUp  = "<<Dmat_r_NUp<<endl;
               cout<<"\nDmat_r_NUz  = "<<Dmat_r_NUz<<endl;*/
              
              commonData.Dmat_RVE[fe_ent](gg).resize(6,6);
              commonData.Dmat_RVE[fe_ent](gg) = Dmat;
              
              commonData.Dmat_RVE_r_Em[fe_ent](gg).resize(6,6);
              commonData.Dmat_RVE_r_Em[fe_ent](gg) = Dmat_r_Em;
              
              commonData.Dmat_RVE_r_NUm[fe_ent](gg).resize(6,6);
              commonData.Dmat_RVE_r_NUm[fe_ent](gg) = Dmat_r_NUm;
              
              commonData.Dmat_RVE_r_Ep[fe_ent](gg).resize(6,6);
              commonData.Dmat_RVE_r_Ep[fe_ent](gg) = Dmat_r_Ep;
              
              commonData.Dmat_RVE_r_Ez[fe_ent](gg).resize(6,6);
              commonData.Dmat_RVE_r_Ez[fe_ent](gg) = Dmat_r_Ez;
              
              commonData.Dmat_RVE_r_NUp[fe_ent](gg).resize(6,6);
              commonData.Dmat_RVE_r_NUp[fe_ent](gg) = Dmat_r_NUp;
              
              commonData.Dmat_RVE_r_NUz[fe_ent](gg).resize(6,6);
              commonData.Dmat_RVE_r_NUz[fe_ent](gg) = Dmat_r_NUz;
              
              commonData.Dmat_RVE_r_Gzp[fe_ent](gg).resize(6,6);
              commonData.Dmat_RVE_r_Gzp[fe_ent](gg) = Dmat_r_Gzp;
              
              ierr = VecDestroy(&F1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&F2); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&F3); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&F4); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&F5); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&F6); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              
              ierr = VecDestroy(&dF1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&dF2); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&dF3); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&dF4); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&dF5); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&dF6); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              
              ierr = VecDestroy(&D1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&dD1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              
              ierr = VecDestroy(&Stress_Homo);    CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&Stress_Homo_r);  CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = VecDestroy(&RVE_volume_Vec); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              
              ierr = MatDestroy(&A); CHKERRABORT(PETSC_COMM_WORLD,ierr);
              ierr = KSPDestroy(&solver_RVE); CHKERRQ(ierr);
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
    
    PetscErrorCode setRVE_DmatRhsOperators(FieldInterface &m_field_RVE,
                                           string field_name,
                                           string wt_field_name,
                                           int nvars,int nders,
                                           vector<string> stochasticfields,
                                           ublas::vector<double> matprop,
                                           int num_rvars,
                                           vector<string> vars_name) {
      PetscFunctionBegin;
      //first calculate wt at each gauss point
      feRhs.getOpPtrVector().push_back(new OpGetWtAtGaussPts(wt_field_name,commonData));
      //At each gauss point run RVE with its own mesh
      feRhs.getOpPtrVector().push_back(new OpCalculate_RVEDmat(m_field_RVE,
                                                               field_name,
                                                               commonData,
                                                               nvars,nders,
                                                               stochasticfields,
                                                               matprop,
                                                               num_rvars,vars_name));
      
      PetscFunctionReturn(0);
    }
    
  };

  
  /*****************************************************************************
   *                                                                           *
   *                                                                           *
   *                         MULTIPLE LAYERS LAMINATE                          *
   *                                                                           *
   *                                                                           *
   ****************************************************************************/
  
  struct FE2_Macro_Laminate_Degradation {
    
    map<EntityHandle, ublas::vector<ublas::matrix<double> > > Macro_GP_Dmat;
    map<EntityHandle, ublas::vector<ublas::matrix<double> > > Macro_GP_Dmat_r_Em;
    map<EntityHandle, ublas::vector<ublas::matrix<double> > > Macro_GP_Dmat_r_NUm;
    map<EntityHandle, ublas::vector<ublas::matrix<double> > > Macro_GP_Dmat_r_Ep;
    map<EntityHandle, ublas::vector<ublas::matrix<double> > > Macro_GP_Dmat_r_Ez;
    map<EntityHandle, ublas::vector<ublas::matrix<double> > > Macro_GP_Dmat_r_NUp;
    map<EntityHandle, ublas::vector<ublas::matrix<double> > > Macro_GP_Dmat_r_NUz;
    map<EntityHandle, ublas::vector<ublas::matrix<double> > > Macro_GP_Dmat_r_Gzp;
    
    ublas::matrix<double> REL_Dmat;
    ublas::matrix<double> REL_Dmat_r_Em;
    ublas::matrix<double> REL_Dmat_r_NUm;
    ublas::matrix<double> REL_Dmat_r_Ep;
    ublas::matrix<double> REL_Dmat_r_Ez;
    ublas::matrix<double> REL_Dmat_r_NUp;
    ublas::matrix<double> REL_Dmat_r_NUz;
    ublas::matrix<double> REL_Dmat_r_Gzp;
    
    // =========================================================================
    //
    //  B.VI. SOLUTION PHASE:
    //        Solve Macroscale FE equation
    //
    // =========================================================================
    virtual PetscErrorCode Calc_Dmat_at_Macro_GP(FieldInterface &m_field_Macro,
                                                 FieldInterface &m_field_RVE,
                                                 int &nvars, int &nders,
                                                 vector<string> &stochastic_fields,
                                                 ublas::vector<double> TheVariables,
                                                 int num_rvars,
                                                 vector<string> vars_name,
                                                 ublas::vector<double> PlyAngle,
                                                 PetscInt NO_Layers,
                                                 int TimePoint) {
      PetscFunctionBegin;
      
      ErrorCode rval;
      PetscErrorCode ierr;
      Get_Dmat_4_REL   Obj_Get_Dmat_4_REL(m_field_Macro);
      FE2_RVE_Dmat_Degradation_Disp Calc_RVE_Dmat(m_field_Macro);
      
      ierr = Calc_RVE_Dmat.setRVE_DmatRhsOperators(m_field_RVE, "DISP_MACRO","Wt",nvars,nders,stochastic_fields,TheVariables,num_rvars,vars_name); CHKERRQ(ierr);
      ierr = Obj_Get_Dmat_4_REL.setMacro_DmatRhsOperators(m_field_Macro, "DISP_MACRO",
                                                          Calc_RVE_Dmat.commonData.Dmat_RVE,
                                                          Calc_RVE_Dmat.commonData.Dmat_RVE_r_Em,
                                                          Calc_RVE_Dmat.commonData.Dmat_RVE_r_NUm,
                                                          Calc_RVE_Dmat.commonData.Dmat_RVE_r_Ep,
                                                          Calc_RVE_Dmat.commonData.Dmat_RVE_r_Ez,
                                                          Calc_RVE_Dmat.commonData.Dmat_RVE_r_NUp,
                                                          Calc_RVE_Dmat.commonData.Dmat_RVE_r_NUz,
                                                          Calc_RVE_Dmat.commonData.Dmat_RVE_r_Gzp,
                                                          REL_Dmat,
                                                          REL_Dmat_r_Em,
                                                          REL_Dmat_r_NUm,
                                                          REL_Dmat_r_Ep,
                                                          REL_Dmat_r_Ez,
                                                          REL_Dmat_r_NUp,
                                                          REL_Dmat_r_NUz,
                                                          REL_Dmat_r_Gzp); CHKERRQ(ierr);
      
      SeriesRecorder *recorder_ptr;
      ierr = m_field_Macro.query_interface(recorder_ptr); CHKERRQ(ierr);
      int count = 0;
      if( recorder_ptr->check_series("Wt_SERIES") ) {
        cout<<"============== Wt_SERIES exists =============== "<<endl;
        for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder_ptr,"Wt_SERIES",sit)) {
          cout<<"\n**************************************";
          cout<<"\n*      Time step: "<<count;
          cout<<"\n**************************************"<<endl;
          //if(count%8==0) {
          if (count == TimePoint) {
            
            PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());
            ierr = recorder_ptr->load_series_data("Wt_SERIES",sit->get_step_number()); CHKERRQ(ierr);
            
            cout<<"\n";
            cout<<"======================================\n";
            cout<<"      Start of RVE calculation        \n";
            cout<<"======================================"<<endl;
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",Calc_RVE_Dmat.getLoopFeRhs()); CHKERRQ(ierr);
            
            cout<<"\n";
            cout<<"======================================\n";
            cout<<"        End of RVE calculation        \n";
            cout<<"======================================"<<endl;
            
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Obj_Get_Dmat_4_REL.getLoopFeRhs()); CHKERRQ(ierr);
            
            PetscPrintf(PETSC_COMM_WORLD,"End of step %d\n",sit->get_step_number());
          }
          count++;
        }
      }
      
      Macro_GP_Dmat       = Calc_RVE_Dmat.commonData.Dmat_RVE;
      Macro_GP_Dmat_r_Em  = Calc_RVE_Dmat.commonData.Dmat_RVE_r_Em;
      Macro_GP_Dmat_r_NUm = Calc_RVE_Dmat.commonData.Dmat_RVE_r_NUm;
      Macro_GP_Dmat_r_Ep  = Calc_RVE_Dmat.commonData.Dmat_RVE_r_Ep;
      Macro_GP_Dmat_r_Ez  = Calc_RVE_Dmat.commonData.Dmat_RVE_r_Ez;
      Macro_GP_Dmat_r_NUp = Calc_RVE_Dmat.commonData.Dmat_RVE_r_NUp;
      Macro_GP_Dmat_r_NUz = Calc_RVE_Dmat.commonData.Dmat_RVE_r_NUz;
      Macro_GP_Dmat_r_Gzp = Calc_RVE_Dmat.commonData.Dmat_RVE_r_Gzp;
      //REL_Dmat = Obj_Get_Dmat_4_REL.Output_Dmat;
      
      PetscFunctionReturn(0);
    }
    
    virtual PetscErrorCode Macro_FE_REL(FieldInterface &m_field_Macro,
                                        FieldInterface &m_field_RVE,
                                        int &nvars, int &nders,
                                        vector<string> &stochastic_fields,
                                        ublas::vector<double> TheVariables,
                                        int num_rvars,
                                        vector<string> vars_name,
                                        ublas::vector<double> PlyAngle,
                                        PetscInt NO_Layers,
                                        int TimePoint) {
      PetscFunctionBegin;
      
      ErrorCode rval;
      PetscErrorCode ierr;
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      //create matrices
      
      Vec F, dF;
      Vec D, dD;
      Mat A;
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&F); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&dF); CHKERRQ(ierr);
      
      ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
      ierr = VecDuplicate(dF,&dD); CHKERRQ(ierr);
      
      ierr = m_field_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_MACRO",&A); CHKERRQ(ierr);
      
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      
      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      
      //External forces, This vector is assemble only once (as this is not function of Dmat)
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
      
      // ierr = MetaNeummanForces::setNeumannFiniteElementOperators(m_field_Macro,neumann_forces,F,"DISP_MACRO"); CHKERRQ(ierr);
      
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
      
      //  ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      Vec Fint;
      ierr = VecDuplicate(F,&Fint); CHKERRQ(ierr);
      
      
      PostPocOnRefinedMesh post_proc(m_field_Macro);
      ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
      ierr = post_proc.addFieldValuesPostProc("DISP_MACRO"); CHKERRQ(ierr);
      
      //read time series and do thermo elastci analysis
      SeriesRecorder *recorder_ptr;
      ierr = m_field_Macro.query_interface(recorder_ptr); CHKERRQ(ierr);
      
      // Here we use ElasticFEMethod_Dmat_input, so will multiply Fint with -1
      DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field_Macro,"DISP_MACRO",A,D,F);
      
      // Here pointer to object is used instead of object, because the pointer
      //   will be destroyed at the end of code before PetscFinalize to make
      //   sure all its internal matrices and vectores are destroyed.
      
      ElasticFEMethod_Dmat_input my_fe(m_field_Macro,A,D,Fint,0.0,0.0,Macro_GP_Dmat,"DISP_MACRO");
      
      //preproc
      ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
      
      
      // We need to assemble matrix A and internal force vector Fint at each
      // time step as these depends on Dmat, which will change at each time step
      ierr = MatZeroEntries(A); CHKERRQ(ierr);
      
      ierr = VecZeroEntries(Fint); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(Fint,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(Fint,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      
      ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      //      ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      
      
      // ------------------------
      // calculate Dmat for all Gauss points in the macro-mesh
      // ------------------------
      /*cout<<"\n\nSet 01"<<endl;
       ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",Calc_RVE_Dmat.getLoopFeRhs()); CHKERRQ(ierr);*/
      
      ierr = VecScale(Fint,-1); CHKERRQ(ierr); //Multiply Fint with -1 (Fint=-Fint)
      ierr = VecAXPY(Fint,1,F); CHKERRQ(ierr); //Fint=Fint+F
      
      //      cin>>wait;
      //      map<EntityHandle, ublas::vector<ublas::matrix<double> > >::iterator mit = calculate_rve_dmat.commonData.Dmat_RVE.begin();
      //      for(;mit!=calculate_rve_dmat.commonData.Dmat_RVE.end();mit++) {
      //        cerr << mit->first << " " << mit->second << endl;
      //      }
      cout<<"\n\nSet 02"<<endl;
      //loop over macro elemnts to assemble A matrix and Fint vector
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe);  CHKERRQ(ierr);
      
      my_dirichlet_bc.snes_B = A;
      my_dirichlet_bc.snes_x = D;
      my_dirichlet_bc.snes_f = Fint;
      cout<<"\n\nSet 03"<<endl;
      //postproc
      ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
      cout<<"\n\nSet 04"<<endl;
      //Solver
      KSP solver_Macro;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver_Macro); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver_Macro,A,A); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver_Macro); CHKERRQ(ierr);
      ierr = KSPSetUp(solver_Macro); CHKERRQ(ierr);
      
      ierr = KSPSolve(solver_Macro,Fint,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      
      //Save data on mesh
      ierr = m_field_Macro.set_local_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",post_proc); CHKERRQ(ierr);
      /*ostringstream o1;
      o1 << "FE2_out_" << sit->step_number << ".h5m";
      rval = post_proc.postProcMesh.write_file(o1.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);*/
      
      /*if(count==100){
       // save the solution file for subsequent strain calculation analysis
       ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
       rval = moab_Macro.write_file("FE2_solution_100.h5m"); CHKERR_PETSC(rval);
       }*/
      
      cout<<"\n\nThe Zeroth-order problem is done!"<<endl;
      /***********************************************************************
       *
       *  3. SOLVE THE FIRST- AND SECOND-ORDER FE EQUILIBRIUM EQUATION
       *     1st order: [K][U_r]  = -[K_r][U}
       *     2nd order: [K][U_rs] = -[K_rs][U]-2[K_r][U_s]
       *
       **********************************************************************/
      
      for (int irv = 1; irv <= num_rvars; irv++) {
        
        /***********************************************************************
         *
         * 3.1. Case 1: Material properties are  treated as random variables
         *
         **********************************************************************/
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[irv];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,3,"NUz") == 0) {
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
        else if (VariableName.compare(0,2,"Ef") == 0) {
          idx_disp = 7;
        }
        else if (VariableName.compare(0,3,"NUf") == 0) {
          idx_disp = 8;
        }
        else if (VariableName.compare(0,5,"force") == 0) {
          idx_disp = 80;
        }
        else if (VariableName.compare(0,11,"orientation") == 0) {
          idx_disp = 90;
        }
        else if (VariableName.compare(0,6,"theta1") == 0) {
          idx_disp = 91;
        }
        else if (VariableName.compare(0,6,"theta2") == 0) {
          idx_disp = 92;
        }
        else if (VariableName.compare(0,6,"theta3") == 0) {
          idx_disp = 93;
        }
        else if (VariableName.compare(0,6,"theta4") == 0) {
          idx_disp = 94;
        }
        VariableName.clear();
        
        // initiation
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        switch (idx_disp) {
          case 0: { cout<<"\n\nWith respect to Em"<<endl;// w.r.t. - Em
            FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_Em(m_field_Macro,A,dD,dF,Macro_GP_Dmat_r_Em,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Em(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Em);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
            break;
          }
          case 1: { cout<<"\n\nWith respect to NUm"<<endl;// w.r.t. - NUm
            FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_NUm(m_field_Macro,A,D,dF,Macro_GP_Dmat_r_NUm,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUm(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUm);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
            break;
          }
          case 2: { cout<<"\n\nWith respect to NUp"<<endl;// w.r.t. - NUp
            FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_NUp(m_field_Macro,A,dD,dF,Macro_GP_Dmat_r_NUp,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUp);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
            break;
          }
          case 3: { cout<<"\n\nWith respect to NUz"<<endl;// w.r.t. - NUz
            FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_NUz(m_field_Macro,A,dD,dF,Macro_GP_Dmat_r_NUz,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUz(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUz); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUz);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUz); CHKERRQ(ierr);
            break;
          }
          case 4: { cout<<"\n\nWith respect to Ep"<<endl;// w.r.t. - Ep
            FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_Ep(m_field_Macro,A,dD,dF,Macro_GP_Dmat_r_Ep,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ep(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Ep);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
            break;
          }
          case 5: { cout<<"\n\nWith respect to Ez"<<endl;// w.r.t. - Ez
            FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_Ez(m_field_Macro,A,dD,dF,Macro_GP_Dmat_r_Ez,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ez(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Ez);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
            break;
          }
          case 6: { cout<<"\n\nWith respect to Gzp"<<endl;// w.r.t. - Gzp
            FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_Gzp(m_field_Macro,A,dD,dF,Macro_GP_Dmat_r_Gzp,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Gzp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Gzp);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
            break;
          }
        }
        
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ss_field.str(""); ss_field.clear();
          ss_field << "DISP_MACRO" << stochastic_fields[idx_disp];
          if (idx_disp<nvars) {
            ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
            
            ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);//ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          }
        }
        
        /***********************************************************************
         *
         * 3.2. Case 2: Applied forces are treated as random variables
         *
         **********************************************************************/
        
        if (idx_disp == 80) {
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ElasticFEMethod_Dmat_input my_fe_1st_Ply_r_F(m_field_Macro,A,D,F,0.0,0.0,Macro_GP_Dmat,"DISP_MACRO");
          
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
          
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe_1st_Ply_r_F);     CHKERRQ(ierr);
          
          /*my_dirichlet_bc.snes_B = A;
          my_dirichlet_bc.snes_x = D;
          my_dirichlet_bc.snes_f = Fint;*/
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
          
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
        }
        
        
      }
      cout<<"\n\nThe First-order problem is done!"<<endl;
      
      // --------------------------------------------
      ierr = KSPDestroy(&solver_Macro); CHKERRQ(ierr);
      
      //        string wait;
      //        cin>>wait;
      
      /***************************************************************************
       *
       *  4. FINISH
       *
       **************************************************************************/
      //Destroy matrices
      ierr = VecDestroy(&F); CHKERRQ(ierr);
      ierr = VecDestroy(&Fint); CHKERRQ(ierr);
      ierr = VecDestroy(&dF); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = VecDestroy(&dD); CHKERRQ(ierr);
      ierr = MatDestroy(&A); CHKERRQ(ierr);
      //ierr = PetscFinalize(); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    virtual PetscErrorCode Macro_FE_REL2(FieldInterface &m_field_Macro,
                                        FieldInterface &m_field_RVE,
                                        int &nvars, int &nders,
                                        vector<string> &stochastic_fields,
                                        ublas::vector<double> TheVariables,
                                        int num_rvars,
                                        vector<string> vars_name,
                                        ublas::vector<double> PlyAngle,
                                        PetscInt NO_Layers,
                                        int TimePoint) {
      PetscFunctionBegin;
      
      ErrorCode rval;
      PetscErrorCode ierr;
      
      FE2_RVE_Dmat_Degradation_Disp Calc_RVE_Dmat(m_field_Macro);
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      //create matrices
      
      Vec F, dF;
      Vec D, dD;
      Mat A;
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&F); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&dF); CHKERRQ(ierr);
      
      ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
      ierr = VecDuplicate(dF,&dD); CHKERRQ(ierr);
      
      ierr = m_field_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_MACRO",&A); CHKERRQ(ierr);
      
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      
      //External forces, This vector is assemble only once (as this is not function of Dmat)
      boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
      ierr = MetaNeummanForces::setNeumannFiniteElementOperators(m_field_Macro,neumann_forces,F,"DISP_MACRO"); CHKERRQ(ierr);
      boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
      for(;mit!=neumann_forces.end();mit++) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
      }
      
      //  ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      
      ierr = Calc_RVE_Dmat.setRVE_DmatRhsOperators(m_field_RVE, "DISP_MACRO","Wt",nvars,nders,stochastic_fields,TheVariables,num_rvars,vars_name); CHKERRQ(ierr);
      Vec Fint;
      ierr = VecDuplicate(F,&Fint); CHKERRQ(ierr);
      
      
      PostPocOnRefinedMesh post_proc(m_field_Macro);
      ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
      ierr = post_proc.addFieldValuesPostProc("DISP_MACRO"); CHKERRQ(ierr);
      
      //read time series and do thermo elastci analysis
      SeriesRecorder *recorder_ptr;
      ierr = m_field_Macro.query_interface(recorder_ptr); CHKERRQ(ierr);
      int count = 0;
      if( recorder_ptr->check_series("Wt_SERIES") ) {
        cout<<"============== Wt_SERIES exists =============== "<<endl;
        for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder_ptr,"Wt_SERIES",sit)) {
          cout<<"\n**************************************";
          cout<<"\n*      Time step: "<<count;
          cout<<"\n**************************************"<<endl;
          //if(count%8==0) {
          if (count == TimePoint) {
            
            PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());
            ierr = recorder_ptr->load_series_data("Wt_SERIES",sit->get_step_number()); CHKERRQ(ierr);
            
            cout<<"\n\nSet 01"<<endl;
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",Calc_RVE_Dmat.getLoopFeRhs()); CHKERRQ(ierr);
            
            // Here we use ElasticFEMethod_Dmat_input, so will multiply Fint with -1
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field_Macro,"DISP_MACRO",A,D,F);
            
            // Here pointer to object is used instead of object, because the pointer
            //   will be destroyed at the end of code before PetscFinalize to make
            //   sure all its internal matrices and vectores are destroyed.
            
            ElasticFEMethod_Dmat_input my_fe(m_field_Macro,A,D,Fint,0.0,0.0,Calc_RVE_Dmat.commonData.Dmat_RVE,"DISP_MACRO");
            
            //preproc
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
            
            
            // We need to assemble matrix A and internal force vector Fint at each
            // time step as these depends on Dmat, which will change at each time step
            ierr = MatZeroEntries(A); CHKERRQ(ierr);
            ierr = VecZeroEntries(Fint); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(Fint,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(Fint,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            
            ierr = VecZeroEntries(D); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            
            ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            //      ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
            
            
            // ------------------------
            // calculate Dmat for all Gauss points in the macro-mesh
            // ------------------------
            /*cout<<"\n\nSet 01"<<endl;
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",Calc_RVE_Dmat.getLoopFeRhs()); CHKERRQ(ierr);*/
            
            ierr = VecScale(Fint,-1); CHKERRQ(ierr); //Multiply Fint with -1 (Fint=-Fint)
            ierr = VecAXPY(Fint,1,F); CHKERRQ(ierr); //Fint=Fint+F
            
            //      cin>>wait;
            //      map<EntityHandle, ublas::vector<ublas::matrix<double> > >::iterator mit = calculate_rve_dmat.commonData.Dmat_RVE.begin();
            //      for(;mit!=calculate_rve_dmat.commonData.Dmat_RVE.end();mit++) {
            //        cerr << mit->first << " " << mit->second << endl;
            //      }
            cout<<"\n\nSet 02"<<endl;
            //loop over macro elemnts to assemble A matrix and Fint vector
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe);  CHKERRQ(ierr);
            
            my_dirichlet_bc.snes_B = A;
            my_dirichlet_bc.snes_x = D;
            my_dirichlet_bc.snes_f = Fint;
            cout<<"\n\nSet 03"<<endl;
            //postproc
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
            cout<<"\n\nSet 04"<<endl;
            //Solver
            KSP solver_Macro;
            ierr = KSPCreate(PETSC_COMM_WORLD,&solver_Macro); CHKERRQ(ierr);
            ierr = KSPSetOperators(solver_Macro,A,A); CHKERRQ(ierr);
            ierr = KSPSetFromOptions(solver_Macro); CHKERRQ(ierr);
            ierr = KSPSetUp(solver_Macro); CHKERRQ(ierr);
            
            ierr = KSPSolve(solver_Macro,Fint,D); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            
            //Save data on mesh
            ierr = m_field_Macro.set_local_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",post_proc); CHKERRQ(ierr);
            ostringstream o1;
            o1 << "FE2_out_" << sit->step_number << ".h5m";
            rval = post_proc.postProcMesh.write_file(o1.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
            
            /*if(count==100){
              // save the solution file for subsequent strain calculation analysis
              ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              rval = moab_Macro.write_file("FE2_solution_100.h5m"); CHKERR_PETSC(rval);
            }*/
            
            cout<<"\n\nThe Zeroth-order problem is done!"<<endl;
            /***********************************************************************
             *
             *  3. SOLVE THE FIRST- AND SECOND-ORDER FE EQUILIBRIUM EQUATION
             *     1st order: [K][U_r]  = -[K_r][U}
             *     2nd order: [K][U_rs] = -[K_rs][U]-2[K_r][U_s]
             *
             **********************************************************************/
            
            for (int irv = 1; irv <= num_rvars; irv++) {
              
              int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
              string VariableName = vars_name[irv];
              if (VariableName.compare(0,2,"Em") == 0) {
                idx_disp = 0;
              }
              else if (VariableName.compare(0,3,"NUm") == 0) {
                idx_disp = 1;
              }
              else if (VariableName.compare(0,3,"NUp") == 0) {
                idx_disp = 2;
              }
              else if (VariableName.compare(0,3,"NUz") == 0) {
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
              else if (VariableName.compare(0,2,"Ef") == 0) {
                idx_disp = 7;
              }
              else if (VariableName.compare(0,3,"NUf") == 0) {
                idx_disp = 8;
              }
              else if (VariableName.compare(0,5,"force") == 0) {
                idx_disp = 80;
              }
              else if (VariableName.compare(0,11,"orientation") == 0) {
                idx_disp = 90;
              }
              else if (VariableName.compare(0,6,"theta1") == 0) {
                idx_disp = 91;
              }
              else if (VariableName.compare(0,6,"theta2") == 0) {
                idx_disp = 92;
              }
              else if (VariableName.compare(0,6,"theta3") == 0) {
                idx_disp = 93;
              }
              else if (VariableName.compare(0,6,"theta4") == 0) {
                idx_disp = 94;
              }
              VariableName.clear();
              
              // initiation
              if ((idx_disp >=0) && (idx_disp < nvars)) {
                ierr = VecZeroEntries(dD); CHKERRQ(ierr);
                ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                ierr = VecZeroEntries(dF); CHKERRQ(ierr);
                ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              }
              
              switch (idx_disp) {
                case 0: { cout<<"\n\nWith respect to Em"<<endl;// w.r.t. - Em
                  FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_Em(m_field_Macro,A,dD,dF,Calc_RVE_Dmat.commonData.Dmat_RVE_r_Em,"DISP_MACRO");
                  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Em(m_field_Macro,"DISP_MACRO",A,dD,dF);
                  ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
                  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Em);  CHKERRQ(ierr);
                  ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
                  break;
                }
                case 1: { cout<<"\n\nWith respect to NUm"<<endl;// w.r.t. - NUm
                  FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_NUm(m_field_Macro,A,D,dF,Calc_RVE_Dmat.commonData.Dmat_RVE_r_NUm,"DISP_MACRO");
                  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUm(m_field_Macro,"DISP_MACRO",A,dD,dF);
                  ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
                  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUm);  CHKERRQ(ierr);
                  ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
                  break;
                }
                case 2: { cout<<"\n\nWith respect to NUp"<<endl;// w.r.t. - NUp
                  FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_NUp(m_field_Macro,A,dD,dF,Calc_RVE_Dmat.commonData.Dmat_RVE_r_NUp,"DISP_MACRO");
                  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUp(m_field_Macro,"DISP_MACRO",A,dD,dF);
                  ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
                  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUp);  CHKERRQ(ierr);
                  ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
                  break;
                }
                case 3: { cout<<"\n\nWith respect to NUz"<<endl;// w.r.t. - NUz
                  FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_NUz(m_field_Macro,A,dD,dF,Calc_RVE_Dmat.commonData.Dmat_RVE_r_NUz,"DISP_MACRO");
                  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUz(m_field_Macro,"DISP_MACRO",A,dD,dF);
                  ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUz); CHKERRQ(ierr);
                  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUz);  CHKERRQ(ierr);
                  ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUz); CHKERRQ(ierr);
                  break;
                }
                case 4: { cout<<"\n\nWith respect to Ep"<<endl;// w.r.t. - Ep
                  FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_Ep(m_field_Macro,A,dD,dF,Calc_RVE_Dmat.commonData.Dmat_RVE_r_Ep,"DISP_MACRO");
                  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ep(m_field_Macro,"DISP_MACRO",A,dD,dF);
                  ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
                  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Ep);  CHKERRQ(ierr);
                  ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
                  break;
                }
                case 5: { cout<<"\n\nWith respect to Ez"<<endl;// w.r.t. - Ez
                  FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_Ez(m_field_Macro,A,dD,dF,Calc_RVE_Dmat.commonData.Dmat_RVE_r_Ez,"DISP_MACRO");
                  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ez(m_field_Macro,"DISP_MACRO",A,dD,dF);
                  ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
                  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Ez);  CHKERRQ(ierr);
                  ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
                  break;
                }
                case 6: { cout<<"\n\nWith respect to Gzp"<<endl;// w.r.t. - Gzp
                  FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_Gzp(m_field_Macro,A,dD,dF,Calc_RVE_Dmat.commonData.Dmat_RVE_r_Gzp,"DISP_MACRO");
                  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Gzp(m_field_Macro,"DISP_MACRO",A,dD,dF);
                  ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
                  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Gzp);  CHKERRQ(ierr);
                  ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
                  break;
                }
              }
              
              ostringstream ss_field;
              if ((idx_disp >=0) && (idx_disp < nvars)) {
                ss_field.str(""); ss_field.clear();
                ss_field << "DISP_MACRO" << stochastic_fields[idx_disp];
                if (idx_disp<nvars) {
                  ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
                  ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
                  
                  ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);//ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                }
              }
            }
            //cout<<"\n\nThe First-order problem is done!"<<endl;
            
            // --------------------------------------------
            ierr = KSPDestroy(&solver_Macro); CHKERRQ(ierr);
            PetscPrintf(PETSC_COMM_WORLD,"End of step %d\n",sit->get_step_number());
            //        string wait;
            //        cin>>wait;
          }
          count++;
        }
      }
      
      /***************************************************************************
       *
       *  4. FINISH
       *
       **************************************************************************/
      //Destroy matrices
      ierr = VecDestroy(&F); CHKERRQ(ierr);
      ierr = VecDestroy(&Fint); CHKERRQ(ierr);
      ierr = VecDestroy(&dF); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = VecDestroy(&dD); CHKERRQ(ierr);
      ierr = MatDestroy(&A); CHKERRQ(ierr);
      //ierr = PetscFinalize(); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
  };
}

#endif //__FE2_MACRO_SOLVER_HPP