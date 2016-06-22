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

#ifndef __CALCULATE_RVE_DMAT_R_TRANSISO_DISP_HPP
#define __CALCULATE_RVE_DMAT_R_TRANSISO_DISP_HPP

namespace MoFEM {

  struct Calculate_RVE_Dmat_r_TransIso_Disp {
    
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
    Calculate_RVE_Dmat_r_TransIso_Disp(FieldInterface &m_field):
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
    
    PetscErrorCode addElasticElements(const string field_name,
                                      vector<string> stochastic_fields,
                                      int num_vars,
                                      const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      ErrorCode rval;
      ierr = mField.add_finite_element("ELASTIC_FE_MACRO",MF_ZERO); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row("ELASTIC_FE_MACRO",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col("ELASTIC_FE_MACRO",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("ELASTIC_FE_MACRO",field_name); CHKERRQ(ierr);
      
      // First-order field
      for (int ii = 0; ii < num_vars; ii++) {
        ostringstream ss_field;
        ss_field << field_name << stochastic_fields[ii];
        ierr = mField.modify_finite_element_add_field_data("ELASTIC_FE_MACRO",ss_field.str().c_str()); CHKERRQ(ierr);
      }
      
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
    
    PetscErrorCode addElasticElements_Laminate(vector<string> FE_Group_Name,
                                               int NO_Layers,
                                               const string field_name,
                                               vector<string> stochastic_fields,
                                               int num_vars,
                                               const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
                                      
      PetscFunctionBegin;
      PetscErrorCode ierr;
      ErrorCode rval;
      for (int ife = 0; ife <= NO_Layers; ife++) {
        ierr = mField.add_finite_element(FE_Group_Name[ife],MF_ZERO); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_row(FE_Group_Name[ife],field_name); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_col(FE_Group_Name[ife],field_name); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_data(FE_Group_Name[ife],field_name); CHKERRQ(ierr);
        
        ierr = mField.modify_finite_element_add_field_data(FE_Group_Name[ife],"Wt"); CHKERRQ(ierr);
        
        // First-order field
        for (int ii = 0; ii < num_vars; ii++) {
          ostringstream ss_field;
          ss_field << field_name << stochastic_fields[ii];
          ierr = mField.modify_finite_element_add_field_data(FE_Group_Name[ife],ss_field.str().c_str()); CHKERRQ(ierr);
        }
        
        
        if(mField.check_field(mesh_nodals_positions)) {
          ierr = mField.modify_finite_element_add_field_data(FE_Group_Name[ife],mesh_nodals_positions); CHKERRQ(ierr);
        }
        
      }
      
     /* 
      // loop over all blocksets and get data which name is MAT_ELASTICSET
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
        Range tEts;
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,tEts,true); CHKERR_PETSC(rval);
        ierr = mField.add_ents_to_finite_element_by_TETs(tEts,"ELASTIC_FE_MACRO"); CHKERRQ(ierr);
      }
      */
      
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
      int numders;
      int numvars;
      ublas::vector<double> matprop;
      vector<string> stochastic_fields;
      CommonData &commonData;
      OpCalculate_RVEDmat(FieldInterface &m_field_RVE, const string field_name,
                          CommonData &common_data,int num_ders,int num_vars,
                          vector<string> _stochastic_fields):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name, UserDataOperator::OPROW),
                                                          m_field_RVE(m_field_RVE),
                                                          commonData(common_data),
                                                          numders(num_ders),
                                                          numvars(num_vars),
                                                          stochastic_fields(_stochastic_fields){}
      
      ~OpCalculate_RVEDmat(){
      }


      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        try {
          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          int nb_gauss_pts = data.getN().size1();
          EntityHandle fe_ent = getMoFEMFEPtr()->get_ent(); //handle of finite element
          // cout<<"fe_ent "<<fe_ent <<endl;
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

            cout<<"nb_gauss_pts "<<nb_gauss_pts<<endl;
            
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
              
              //=============================================================================================================

//              cout<<"gg Start =  "<<gg <<endl;
              //We don't need to calculate internal forces for RVE, as ElasticFEMethod is used to assemble A matirx only
              //so noo need to create MyElasticFEMethod here
              ElasticFEMethod_Matrix my_fe_marix(m_field_RVE,A,D1,F1,0.0,0.0,commonData.wtAtGaussPts(gg),"DISP_RVE");
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
              for(int ii = 0; ii < numders; ii++) {
                if (ii < numvars) {
                  ierr = VecZeroEntries(dF1); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                } else {
                  ierr = VecZeroEntries(ddF1); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                }
                if (ii == 0) { // due to Young's modulus of matrix - isotropic
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Em(m_field_RVE,A,D1,dF1,"DISP_RVE","Young","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
                }
                else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUm(m_field_RVE,A,D1,dF1,"DISP_RVE","Poisson","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
                }
                else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUp(m_field_RVE,A,D1,dF1,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
                }
                else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUpz(m_field_RVE,A,D1,dF1,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
                }
                else if (ii == 4) { // due to Young's modulus in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ep(m_field_RVE,A,D1,dF1,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
                }
                else if (ii == 5) { // due to Young's modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ez(m_field_RVE,A,D1,dF1,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
                }
                else if (ii == 6) { // due to shear modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Gzp(m_field_RVE,A,D1,dF1,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
                }
                
                //
                ostringstream ss_field;
                ss_field << "DISP_RVE" << stochastic_fields[ii];
                if (ii < numvars){ // solution for first-order problem
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
                    switch (ii) {
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
              
              //------------------------------------------------------------------------
              // 3.1 Solving the equation to get nodal displacement
              //     case 2: applied macro strain: [0 1 0 0 0 0]^T
              //------------------------------------------------------------------------
              ierr = KSPSolve(solver_RVE,F2,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              
              ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
              //    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
              ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_2(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_2);  CHKERRQ(ierr);
              
              VecGetArray(Stress_Homo, &avec);
              for(int ii=0; ii<6; ii++){
                Dmat(ii,1)=*avec;
                avec++;
              }
              VecRestoreArray(Stress_Homo,&avec);

              //---------------------
              // FIRST & SECOND-ORDER
              //---------------------
              for(int ii = 0; ii < numders; ii++) {
                if (ii < numvars) {
                  ierr = VecZeroEntries(dF2); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                } else {
                  ierr = VecZeroEntries(ddF2); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                }
                if (ii == 0) { // due to Young's modulus of matrix - isotropic
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Em(m_field_RVE,A,D1,dF2,"DISP_RVE","Young","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
                }
                else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUm(m_field_RVE,A,D1,dF2,"DISP_RVE","Poisson","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
                }
                else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUp(m_field_RVE,A,D1,dF2,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
                }
                else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUpz(m_field_RVE,A,D1,dF2,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
                }
                else if (ii == 4) { // due to Young's modulus in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ep(m_field_RVE,A,D1,dF2,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
                }
                else if (ii == 5) { // due to Young's modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ez(m_field_RVE,A,D1,dF2,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
                }
                else if (ii == 6) { // due to shear modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Gzp(m_field_RVE,A,D1,dF2,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
                }
                
                //
                ostringstream ss_field;
                ss_field << "DISP_RVE" << stochastic_fields[ii];
                if (ii < numvars){ // solution for first-order problem
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
                    switch (ii) {
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
              
              //----------------------------------------------------------------
              // 3.1 Solving the equation to get nodal displacement
              //     case 3: applied macro strain: [0 0 1 0 0 0]^T
              //----------------------------------------------------------------
              ierr = KSPSolve(solver_RVE,F3,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              
              //{
              
              //}
              ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
              
              ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_3(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_3);  CHKERRQ(ierr);
              
              VecGetArray(Stress_Homo, &avec);
              for(int ii=0; ii<6; ii++){
                Dmat(ii,2)=*avec;
                avec++;
              }
              VecRestoreArray(Stress_Homo, &avec);
              
              
              //---------------------
              // FIRST & SECOND-ORDER
              //---------------------
              for(int ii = 0; ii < numders; ii++) {
                if (ii < numvars) {
                  ierr = VecZeroEntries(dF3); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                } else {
                  ierr = VecZeroEntries(ddF3); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                }
                if (ii == 0) { // due to Young's modulus of matrix - isotropic
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Em(m_field_RVE,A,D1,dF3,"DISP_RVE","Young","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
                }
                else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUm(m_field_RVE,A,D1,dF3,"DISP_RVE","Poisson","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
                }
                else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUp(m_field_RVE,A,D1,dF3,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
                }
                else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUpz(m_field_RVE,A,D1,dF3,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
                }
                else if (ii == 4) { // due to Young's modulus in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ep(m_field_RVE,A,D1,dF3,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
                }
                else if (ii == 5) { // due to Young's modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ez(m_field_RVE,A,D1,dF3,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
                }
                else if (ii == 6) { // due to shear modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Gzp(m_field_RVE,A,D1,dF3,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
                }
                
                //
                ostringstream ss_field;
                ss_field << "DISP_RVE" << stochastic_fields[ii];
                if (ii < numvars){ // solution for first-order problem
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
                    switch (ii) {
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
              
              //----------------------------------------------------------------
              // 3.1 Solving the equation to get nodal displacement
              //     case 4: applied macro strain: [0 0 0 1 0 0]^T
              //----------------------------------------------------------------
              ierr = KSPSolve(solver_RVE,F4,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
              ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_4(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_4);  CHKERRQ(ierr);
              
              VecGetArray(Stress_Homo, &avec);
              for(int ii=0; ii<6; ii++){
                Dmat(ii,3)=*avec;
                avec++;
              }
              VecRestoreArray(Stress_Homo, &avec);

              
              //---------------------
              // FIRST & SECOND-ORDER
              //---------------------
              for(int ii = 0; ii < numders; ii++) {
                if (ii < numvars) {
                  ierr = VecZeroEntries(dF4); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                } else {
                  ierr = VecZeroEntries(ddF4); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                }
                if (ii == 0) { // due to Young's modulus of matrix - isotropic
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Em(m_field_RVE,A,D1,dF4,"DISP_RVE","Young","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
                }
                else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUm(m_field_RVE,A,D1,dF4,"DISP_RVE","Poisson","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
                }
                else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUp(m_field_RVE,A,D1,dF4,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
                }
                else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUpz(m_field_RVE,A,D1,dF4,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
                }
                else if (ii == 4) { // due to Young's modulus in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ep(m_field_RVE,A,D1,dF4,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
                }
                else if (ii == 5) { // due to Young's modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ez(m_field_RVE,A,D1,dF4,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
                }
                else if (ii == 6) { // due to shear modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Gzp(m_field_RVE,A,D1,dF4,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
                }
                
                //
                ostringstream ss_field;
                ss_field << "DISP_RVE" << stochastic_fields[ii];
                if (ii < numvars){ // solution for first-order problem
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
                    switch (ii) {
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
              
              //----------------------------------------------------------------
              // 3.1 Solving the equation to get nodal displacement
              //     case 5: applied macro strain: [0 0 0 0 1 0]^T
              //----------------------------------------------------------------
              ierr = KSPSolve(solver_RVE,F5,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              
              ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
              
              ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_5(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_5);  CHKERRQ(ierr);
              
              VecGetArray(Stress_Homo, &avec);
              for(int ii=0; ii<6; ii++){
                Dmat(ii,4)=*avec;
                avec++;
              }
              VecRestoreArray(Stress_Homo, &avec);
              
              
              //---------------------
              // FIRST & SECOND-ORDER
              //---------------------
              for(int ii = 0; ii < numders; ii++) {
                if (ii < numvars) {
                  ierr = VecZeroEntries(dF5); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                } else {
                  ierr = VecZeroEntries(ddF5); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                }
                if (ii == 0) { // due to Young's modulus of matrix - isotropic
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Em(m_field_RVE,A,D1,dF5,"DISP_RVE","Young","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
                }
                else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUm(m_field_RVE,A,D1,dF5,"DISP_RVE","Poisson","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
                }
                else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUp(m_field_RVE,A,D1,dF5,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
                }
                else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUpz(m_field_RVE,A,D1,dF5,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
                }
                else if (ii == 4) { // due to Young's modulus in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ep(m_field_RVE,A,D1,dF5,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
                }
                else if (ii == 5) { // due to Young's modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ez(m_field_RVE,A,D1,dF5,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
                }
                else if (ii == 6) { // due to shear modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Gzp(m_field_RVE,A,D1,dF5,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
                }
                
                //
                ostringstream ss_field;
                ss_field << "DISP_RVE" << stochastic_fields[ii];
                if (ii < numvars){ // solution for first-order problem
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
                    switch (ii) {
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
              
              
              //----------------------------------------------------------------
              // 3.1 Solving the equation to get nodal displacement
              //     case 6: applied macro strain: [0 0 0 0 0 1]^T
              //----------------------------------------------------------------
              ierr = KSPSolve(solver_RVE,F6,D1); CHKERRQ(ierr);
              ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              
              ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
              ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_6(m_field_RVE,A,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
              ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_6);  CHKERRQ(ierr);
              
              VecGetArray(Stress_Homo, &avec);
              for(int ii=0; ii<6; ii++){
                Dmat(ii,5)=*avec;
                avec++;
              }
              VecRestoreArray(Stress_Homo, &avec);

              
              //---------------------
              // FIRST & SECOND-ORDER
              //---------------------
              for(int ii = 0; ii < numders; ii++) {
                if (ii < numvars) {
                  ierr = VecZeroEntries(dF6); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                } else {
                  ierr = VecZeroEntries(ddF6); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                }
                if (ii == 0) { // due to Young's modulus of matrix - isotropic
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Em(m_field_RVE,A,D1,dF6,"DISP_RVE","Young","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
                }
                else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUm(m_field_RVE,A,D1,dF6,"DISP_RVE","Poisson","isotropic","matrix",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
                }
                else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUp(m_field_RVE,A,D1,dF6,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
                }
                else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_NUpz(m_field_RVE,A,D1,dF6,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
                }
                else if (ii == 4) { // due to Young's modulus in p-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ep(m_field_RVE,A,D1,dF6,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
                }
                else if (ii == 5) { // due to Young's modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Ez(m_field_RVE,A,D1,dF6,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
                }
                else if (ii == 6) { // due to shear modulus in z-direction of fibre
                  Trans_Iso_Rhs_r_PSFEM_Degradation my_fe_k_r_Gzp(m_field_RVE,A,D1,dF6,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement",commonData.wtAtGaussPts(gg));
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
                }
                
                //
                ostringstream ss_field;
                ss_field << "DISP_RVE" << stochastic_fields[ii];
                if (ii < numvars){ // solution for first-order problem
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
                    switch (ii) {
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
              
              // -
           
              //cout<<"fe_ent "<<fe_ent << "gg  "<<gg<<endl;
              cout<<"\nDmat = "<<Dmat<<endl;
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
                                           
                                           int numders,
                                           int numvars,
                                           vector<string> stochasticfields) {
      PetscFunctionBegin;
      //first calculate wt at each gauss point
      feRhs.getOpPtrVector().push_back(new OpGetWtAtGaussPts(wt_field_name,commonData));
      //At each gauss point run RVE with its own mesh
      feRhs.getOpPtrVector().push_back(new OpCalculate_RVEDmat(m_field_RVE,field_name,commonData,numders,numvars,stochasticfields));

      PetscFunctionReturn(0);
    }

  
  };
  
}

#endif //__CALCULATE_RVE_DMAT_R_TRANSISO_DISP_HPP



