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
   *                         MULTIPLE LAYERS LAMINATE                          *
   *                                                                           *
   *                                                                           *
   ****************************************************************************/
  
  struct FE2_Macro_Solver_Laminate_Degradation {
    
    //global variable Dmat
    ublas::matrix<double> Dmat;
    ublas::matrix<double> Dmat_r_Em;
    ublas::matrix<double> Dmat_r_NUm;
    ublas::matrix<double> Dmat_r_Ep;
    ublas::matrix<double> Dmat_r_Ez;
    ublas::matrix<double> Dmat_r_NUp;
    ublas::matrix<double> Dmat_r_NUpz;
    ublas::matrix<double> Dmat_r_Gzp;
    ublas::matrix<double> Dmat_r_Ef;
    ublas::matrix<double> Dmat_r_NUf;
    
    ublas::matrix<double> Dmat_rs_EmEm;
    ublas::matrix<double> Dmat_rs_NUmNUm;
    ublas::matrix<double> Dmat_rs_EpEp;
    ublas::matrix<double> Dmat_rs_EzEz;
    ublas::matrix<double> Dmat_rs_NUpNUp;
    ublas::matrix<double> Dmat_rs_NUpzNUpz;
    ublas::matrix<double> Dmat_rs_GzpGzp;
    ublas::matrix<double> Dmat_rs_EfEf;
    ublas::matrix<double> Dmat_rs_NUfNUf;
    
    ublas::matrix<double> Dmat_1st_Ply;
    ublas::matrix<double> Dmat_1st_Ply_r_Em;
    ublas::matrix<double> Dmat_1st_Ply_r_NUm;
    ublas::matrix<double> Dmat_1st_Ply_r_Ep;
    ublas::matrix<double> Dmat_1st_Ply_r_Ez;
    ublas::matrix<double> Dmat_1st_Ply_r_NUp;
    ublas::matrix<double> Dmat_1st_Ply_r_NUpz;
    ublas::matrix<double> Dmat_1st_Ply_r_Gzp;
    ublas::matrix<double> Dmat_1st_Ply_r_Ef;
    ublas::matrix<double> Dmat_1st_Ply_r_NUf;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta_1;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta_2;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta_3;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta_4;
    
    ublas::matrix<double> Dmat_2nd_Ply;
    ublas::matrix<double> Dmat_2nd_Ply_r_Em;
    ublas::matrix<double> Dmat_2nd_Ply_r_NUm;
    ublas::matrix<double> Dmat_2nd_Ply_r_Ep;
    ublas::matrix<double> Dmat_2nd_Ply_r_Ez;
    ublas::matrix<double> Dmat_2nd_Ply_r_NUp;
    ublas::matrix<double> Dmat_2nd_Ply_r_NUpz;
    ublas::matrix<double> Dmat_2nd_Ply_r_Gzp;
    ublas::matrix<double> Dmat_2nd_Ply_r_Ef;
    ublas::matrix<double> Dmat_2nd_Ply_r_NUf;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta_1;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta_2;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta_3;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta_4;
    
    ublas::matrix<double> Dmat_3rd_Ply;
    ublas::matrix<double> Dmat_3rd_Ply_r_Em;
    ublas::matrix<double> Dmat_3rd_Ply_r_NUm;
    ublas::matrix<double> Dmat_3rd_Ply_r_Ep;
    ublas::matrix<double> Dmat_3rd_Ply_r_Ez;
    ublas::matrix<double> Dmat_3rd_Ply_r_NUp;
    ublas::matrix<double> Dmat_3rd_Ply_r_NUpz;
    ublas::matrix<double> Dmat_3rd_Ply_r_Gzp;
    ublas::matrix<double> Dmat_3rd_Ply_r_Ef;
    ublas::matrix<double> Dmat_3rd_Ply_r_NUf;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta_1;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta_2;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta_3;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta_4;
    
    ublas::matrix<double> Dmat_4th_Ply;
    ublas::matrix<double> Dmat_4th_Ply_r_Em;
    ublas::matrix<double> Dmat_4th_Ply_r_NUm;
    ublas::matrix<double> Dmat_4th_Ply_r_Ep;
    ublas::matrix<double> Dmat_4th_Ply_r_Ez;
    ublas::matrix<double> Dmat_4th_Ply_r_NUp;
    ublas::matrix<double> Dmat_4th_Ply_r_NUpz;
    ublas::matrix<double> Dmat_4th_Ply_r_Gzp;
    ublas::matrix<double> Dmat_4th_Ply_r_Ef;
    ublas::matrix<double> Dmat_4th_Ply_r_NUf;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta_1;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta_2;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta_3;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta_4;
    
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
    

    
    // =========================================================================
    //
    //  B.VI. SOLUTION PHASE:
    //        Solve Macroscale FE equation
    //
    // =========================================================================
    
    virtual PetscErrorCode Macro_FE_REL(FieldInterface &m_field_Macro,
                                        int &nvars, int &nders,
                                        vector<string> &stochastic_fields,
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
      Vec dF, dD;
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&F); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&D); CHKERRQ(ierr);
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&dF); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&dD); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat A;
      ierr = m_field_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_MACRO",&A); CHKERRQ(ierr);
      
      
      Projection10NodeCoordsOnField ent_method_material_Macro(m_field_Macro,"MESH_NODE_POSITIONS");
      ierr = m_field_Macro.loop_dofs("MESH_NODE_POSITIONS",ent_method_material_Macro); CHKERRQ(ierr);
      
      //Assemble F and A
      DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field_Macro,"DISP_MACRO",A,D,F);
      
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
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
      
      ierr = Calc_RVE_Dmat.setRVE_DmatRhsOperators(m_field_RVE, "DISP_MACRO","Wt",nders,nvars,stochastic_fields); CHKERRQ(ierr);
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
          cout<<"\n**************************************"<<endl;
          cout<<"\n*      Time step: "<<count<<endl;
          cout<<"\n**************************************"<<endl;
          if(count%8==0) {
            //if(count == 0) {
            
            PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());
            ierr = recorder_ptr->load_series_data("Wt_SERIES",sit->get_step_number()); CHKERRQ(ierr);
            
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
            
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",Calc_RVE_Dmat.getLoopFeRhs()); CHKERRQ(ierr);
            
            ierr = VecScale(Fint,-1); CHKERRQ(ierr); //Multiply Fint with -1 (Fint=-Fint)
            ierr = VecAXPY(Fint,1,F); CHKERRQ(ierr); //Fint=Fint+F
            
            //      cin>>wait;
            //      map<EntityHandle, ublas::vector<ublas::matrix<double> > >::iterator mit = calculate_rve_dmat.commonData.Dmat_RVE.begin();
            //      for(;mit!=calculate_rve_dmat.commonData.Dmat_RVE.end();mit++) {
            //        cerr << mit->first << " " << mit->second << endl;
            //      }
            
            //loop over macro elemnts to assemble A matrix and Fint vector
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe);  CHKERRQ(ierr);
            
            my_dirichlet_bc.snes_B = A;
            my_dirichlet_bc.snes_x = D;
            my_dirichlet_bc.snes_f = Fint;
            
            //postproc
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
            
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
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",post_proc); CHKERRQ(ierr);
            ostringstream o1;
            o1 << "FE2_out_" << sit->step_number << ".h5m";
            rval = post_proc.postProcMesh.write_file(o1.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
            
            if(count==100){
              // save the solution file for subsequent strain calculation analysis
              ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              rval = moab_Macro.write_file("FE2_solution_100.h5m"); CHKERR_PETSC(rval);
            }
            
            cout<<"\n\nThe Zeroth-order problem is done!"<<endl;
            /***********************************************************************
             *
             *  3. SOLVE THE FIRST- AND SECOND-ORDER FE EQUILIBRIUM EQUATION
             *     1st order: [K][U_r]  = -[K_r][U}
             *     2nd order: [K][U_rs] = -[K_rs][U]-2[K_r][U_s]
             *
             **********************************************************************/
            
            for (int irv = 0; irv < nvars; irv++) {
              // initiation
              ierr = VecZeroEntries(dD); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecZeroEntries(dF); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              
              switch (irv) {
                case 0: { cout<<"\n\nWith respect to Em"<<endl;// w.r.t. - Em
                  FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_Em(m_field_Macro,A,dD,dF,Calc_RVE_Dmat.commonData.Dmat_RVE_r_Em,"DISP_MACRO");
                  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Em(m_field_Macro,"DISP_MACRO",A,dD,dF);
                  ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
                  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_Em);  CHKERRQ(ierr);
                  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_Em);  CHKERRQ(ierr);
                  ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
                  break;
                }
                case 1: { cout<<"\n\nWith respect to NUm"<<endl;// w.r.t. - NUm
                  FE2_Rhs_r_PSFEM_Degradation my_fe2_k_r_NUm(m_field_Macro,A,D,dF,Calc_RVE_Dmat.commonData.Dmat_RVE_r_NUm,"DISP_MACRO");
                  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUm(m_field_Macro,"DISP_MACRO",A,dD,dF);
                  ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
                  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_NUm);  CHKERRQ(ierr);
                  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_NUm);  CHKERRQ(ierr);
                  ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
                  break;
                }
              }
              
              ostringstream ss_field;
              ss_field.str(""); ss_field.clear();
              ss_field << "DISP_MACRO" << stochastic_fields[irv];
              
              ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
              ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
              
              ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);//ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            }
            cout<<"\n\nThe First-order problem is done!"<<endl;
            
            // --------------------------------------------
            ierr = KSPDestroy(&solver_Macro); CHKERRQ(ierr);
            PetscPrintf(PETSC_COMM_WORLD,"End of step %d\n",sit->get_step_number());
            //        string wait;
            //        cin>>wait;
          }
          count++;
        }
      }
      
      ierr = VecDestroy(&F); CHKERRQ(ierr);
      ierr = VecDestroy(&Fint); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = MatDestroy(&A); CHKERRQ(ierr);
      ierr = PetscFinalize(); CHKERRQ(ierr);
      
      
      
      
      
      
      
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
      
      MyElasticFEMethod_Macro my_fe_1st_ply (m_field_Macro,A,D,F,Dmat_1st_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_2nd_ply (m_field_Macro,A,D,F,Dmat_2nd_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_3rd_ply (m_field_Macro,A,D,F,Dmat_3rd_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_4th_ply (m_field_Macro,A,D,F,Dmat_4th_Ply,"DISP_MACRO");
      
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
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_Em.resize(6,6);   Dmat_1st_Ply_r_Em.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_1st_Ply_r_Em); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_Em.resize(6,6);   Dmat_2nd_Ply_r_Em.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_2nd_Ply_r_Em); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_Em.resize(6,6);   Dmat_3rd_Ply_r_Em.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_3rd_Ply_r_Em); CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_Em.resize(6,6);   Dmat_4th_Ply_r_Em.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_4th_Ply_r_Em); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_Em_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Em,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Em_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Em,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Em_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Em,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Em_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Em,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Em(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Em_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Em_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Em_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Em_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
            break;
          }
          case 1: { // due to Poisson's ratio of matrix (NUm)
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_NUm.resize(6,6);   Dmat_1st_Ply_r_NUm.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_1st_Ply_r_NUm); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_NUm.resize(6,6);   Dmat_2nd_Ply_r_NUm.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_2nd_Ply_r_NUm); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_NUm.resize(6,6);   Dmat_3rd_Ply_r_NUm.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_3rd_Ply_r_NUm); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_NUm.resize(6,6);   Dmat_4th_Ply_r_NUm.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_4th_Ply_r_NUm); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_NUm,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_NUm,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_NUm,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_NUm,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUm(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUm_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_NUm_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_NUm_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_NUm_4th_Ply);  CHKERRQ(ierr);
            }
            
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
            break;
          }
          case 2: {// due to transversal Poisson's ratio of fibre (NUp)
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_NUp.resize(6,6);   Dmat_1st_Ply_r_NUp.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_1st_Ply_r_NUp); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_NUp.resize(6,6);   Dmat_2nd_Ply_r_NUp.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_2nd_Ply_r_NUp); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_NUp.resize(6,6);   Dmat_3rd_Ply_r_NUp.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_3rd_Ply_r_NUp); CHKERRQ(ierr);
            }
            // Fourth-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_NUp.resize(6,6);   Dmat_4th_Ply_r_NUp.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_4th_Ply_r_NUp); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_NUp,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_NUp,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_NUp,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_NUp,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUp_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_NUp_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_NUp_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_NUp_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
            break;
          }
          case 3: {// due to axial Poisson's ratio of fibre (NUpz)
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_NUpz.resize(6,6);   Dmat_1st_Ply_r_NUpz.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_1st_Ply_r_NUpz); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_NUpz.resize(6,6);   Dmat_2nd_Ply_r_NUpz.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_2nd_Ply_r_NUpz); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_NUpz.resize(6,6);   Dmat_3rd_Ply_r_NUpz.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_3rd_Ply_r_NUpz); CHKERRQ(ierr);
            }
            // Fourth-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_NUpz.resize(6,6);   Dmat_4th_Ply_r_NUpz.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_4th_Ply_r_NUpz); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_NUpz,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_NUpz,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_NUpz,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_NUpz,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUpz(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUpz_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_NUpz_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_NUpz_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_NUpz_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
            break;
          }
          case 4: {// due to transversal modulus of fibre (Ep)
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_Ep.resize(6,6);   Dmat_1st_Ply_r_Ep.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_1st_Ply_r_Ep); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_Ep.resize(6,6);   Dmat_2nd_Ply_r_Ep.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_2nd_Ply_r_Ep); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_Ep.resize(6,6);   Dmat_3rd_Ply_r_Ep.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_3rd_Ply_r_Ep); CHKERRQ(ierr);
            }
            // Fourth-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_Ep.resize(6,6);   Dmat_4th_Ply_r_Ep.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_4th_Ply_r_Ep); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ep_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Ep,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ep_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Ep,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ep_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Ep,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ep_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Ep,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ep(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Ep_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Ep_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Ep_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Ep_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
            break;
          }
          case 5: {// due to axial modulus of fibre (Ez)
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_Ez.resize(6,6);   Dmat_1st_Ply_r_Ez.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_1st_Ply_r_Ez); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_Ez.resize(6,6);   Dmat_2nd_Ply_r_Ez.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_2nd_Ply_r_Ez); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_Ez.resize(6,6);   Dmat_3rd_Ply_r_Ez.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_3rd_Ply_r_Ez); CHKERRQ(ierr);
            }
            // Fourth-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_Ez.resize(6,6);   Dmat_4th_Ply_r_Ez.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_4th_Ply_r_Ez); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ez_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Ez,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ez_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Ez,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ez_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Ez,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ez_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Ez,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ez(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Ez_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Ez_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Ez_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Ez_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
            break;
          }
          case 6: {// due to shear modulus of fibre (Gzp)
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_Gzp.resize(6,6);   Dmat_1st_Ply_r_Gzp.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_1st_Ply_r_Gzp); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_Gzp.resize(6,6);   Dmat_2nd_Ply_r_Gzp.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_2nd_Ply_r_Gzp); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_Gzp.resize(6,6);   Dmat_3rd_Ply_r_Gzp.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_3rd_Ply_r_Gzp); CHKERRQ(ierr);
            }
            // Fourth-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_Gzp.resize(6,6);   Dmat_4th_Ply_r_Gzp.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_4th_Ply_r_Gzp); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Gzp,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Gzp,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Gzp,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Gzp,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Gzp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Gzp_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Gzp_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Gzp_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Gzp_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
            break;
          }
          case 7: {// due to Young's modulus of fibre (Ef) - isotropic
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_Ef.resize(6,6);   Dmat_1st_Ply_r_Ef.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ef, Dmat_1st_Ply_r_Ef); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_Ef.resize(6,6);   Dmat_2nd_Ply_r_Ef.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ef, Dmat_2nd_Ply_r_Ef); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_Ef.resize(6,6);   Dmat_3rd_Ply_r_Ef.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ef, Dmat_3rd_Ply_r_Ef); CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_Ef.resize(6,6);   Dmat_4th_Ply_r_Ef.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ef, Dmat_4th_Ply_r_Ef); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ef_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Ef,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ef_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Ef,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ef_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Ef,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ef_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Ef,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ef(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ef); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Ef_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Ef_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Ef_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Ef_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ef); CHKERRQ(ierr);
            break;
          }
          case 8: { // due to Poisson's ratio of fibre (NUf) - isotropic
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_NUf.resize(6,6);   Dmat_1st_Ply_r_NUf.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUf, Dmat_1st_Ply_r_NUf); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_NUf.resize(6,6);   Dmat_2nd_Ply_r_NUf.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUf, Dmat_2nd_Ply_r_NUf); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_NUf.resize(6,6);   Dmat_3rd_Ply_r_NUf.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUf, Dmat_3rd_Ply_r_NUf); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_NUf.resize(6,6);   Dmat_4th_Ply_r_NUf.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUf, Dmat_4th_Ply_r_NUf); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUf_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_NUf,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUf_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_NUf,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUf_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_NUf,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUf_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_NUf,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUf(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUf); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUf_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_NUf_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_NUf_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_NUf_4th_Ply);  CHKERRQ(ierr);
            }
            
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUf); CHKERRQ(ierr);
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
          MyElasticFEMethod_Macro my_fe_1st_Ply_r_F (m_field_Macro,A,D,F,Dmat_1st_Ply,"DISP_MACRO");
          MyElasticFEMethod_Macro my_fe_2nd_Ply_r_F (m_field_Macro,A,D,F,Dmat_2nd_Ply,"DISP_MACRO");
          MyElasticFEMethod_Macro my_fe_3rd_Ply_r_F (m_field_Macro,A,D,F,Dmat_3rd_Ply,"DISP_MACRO");
          MyElasticFEMethod_Macro my_fe_4th_Ply_r_F (m_field_Macro,A,D,F,Dmat_4th_Ply,"DISP_MACRO");
          
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
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe_1st_Ply_r_F);     CHKERRQ(ierr);
          // Third layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe_2nd_Ply_r_F);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe_3rd_Ply_r_F);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe_4th_Ply_r_F);  CHKERRQ(ierr);
          }
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
        
        /***********************************************************************
         *
         * 3.3. Case 3: Ply orientation angle is treated as random variable
         *
         **********************************************************************/
        if (idx_disp == 90) {
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_Theta.resize(6,6);   Dmat_1st_Ply_r_Theta.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_1st_Ply_r_Theta); CHKERRQ(ierr);
          // Second-layer
          if (NO_Layers > 1) {
            theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_2nd_Ply_r_Theta.resize(6,6);   Dmat_2nd_Ply_r_Theta.clear();
            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_2nd_Ply_r_Theta); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 2) {
            theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_3rd_Ply_r_Theta.resize(6,6);   Dmat_3rd_Ply_r_Theta.clear();
            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_3rd_Ply_r_Theta); CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_4th_Ply_r_Theta.resize(6,6);   Dmat_4th_Ply_r_Theta.clear();
            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_4th_Ply_r_Theta); CHKERRQ(ierr);
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Theta_1st_Ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Theta_2nd_Ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Theta_3rd_Ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Theta_4th_Ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if (idx_disp == 91) { // first layer angle is a random variable
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          //cout<<"\n\nThe angle case Theta_1 starts from here!!!"<<endl;
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_Theta_1.resize(6,6);   Dmat_1st_Ply_r_Theta_1.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_1st_Ply_r_Theta_1); CHKERRQ(ierr);
          //cout<<"\n\nThe first ply angle is "<<PlyAngle(0)<<endl;
          //cout<<"\n\n Dmat_1st_Ply_r_Theta_1 is "<<Dmat_1st_Ply_r_Theta_1<<endl;
          // Second-layer
          if (NO_Layers > 1) {
            Dmat_2nd_Ply_r_Theta_1.resize(6,6);   Dmat_2nd_Ply_r_Theta_1.clear();
          }
          // Third-layer
          if (NO_Layers > 2) {
            Dmat_3rd_Ply_r_Theta_1.resize(6,6);   Dmat_3rd_Ply_r_Theta_1.clear();
          }
          // Fourth layer
          if (NO_Layers > 3) {
            Dmat_4th_Ply_r_Theta_1.resize(6,6);   Dmat_4th_Ply_r_Theta_1.clear();
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta_1,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta_1,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta_1,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta_1,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_1(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_1); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_1);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_1);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_1);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_1);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_1); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta_1st_Ply",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"Solving the first order derivative of Theta_1"<<endl;
        }
        else if (idx_disp == 92) { // second layer angle is a random variable
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          //cout<<"\n\nThe angle case Theta_2 starts from here!!!"<<endl;
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          Dmat_1st_Ply_r_Theta_2.resize(6,6);   Dmat_1st_Ply_r_Theta_2.clear();
          // Second-layer
          theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_2nd_Ply_r_Theta_2.resize(6,6);   Dmat_2nd_Ply_r_Theta_2.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_2nd_Ply_r_Theta_2); CHKERRQ(ierr);
          // Third-layer
          if (NO_Layers > 2) {
            Dmat_3rd_Ply_r_Theta_2.resize(6,6);   Dmat_3rd_Ply_r_Theta_2.clear();
          }
          // Fourth layer
          if (NO_Layers > 3) {
            Dmat_4th_Ply_r_Theta_2.resize(6,6);   Dmat_4th_Ply_r_Theta_2.clear();
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta_2,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta_2,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta_2,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta_2,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_2(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_2); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_2);  CHKERRQ(ierr);
          // Second layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_2);  CHKERRQ(ierr);
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_2);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_2);  CHKERRQ(ierr);
          }
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_2); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta_2nd_Ply",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"Solving the first order derivative of Theta_2"<<endl;
        }
        else if (idx_disp == 93) { // third layer angle is a random variable
          // Initiate the involved parameters to zero
          //cout<<"\n\nThe angle case Theta_3 starts from here!!!"<<endl;
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          Dmat_1st_Ply_r_Theta_3.resize(6,6);   Dmat_1st_Ply_r_Theta_3.clear();
          // Second-layer
          Dmat_2nd_Ply_r_Theta_3.resize(6,6);   Dmat_2nd_Ply_r_Theta_3.clear();
          // Third-layer
          theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_3rd_Ply_r_Theta_3.resize(6,6);   Dmat_3rd_Ply_r_Theta_3.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_3rd_Ply_r_Theta_3); CHKERRQ(ierr);
          // Fourth layer
          if (NO_Layers > 3) {
            Dmat_4th_Ply_r_Theta_3.resize(6,6);   Dmat_4th_Ply_r_Theta_3.clear();
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta_3,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta_3,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta_3,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta_3,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_3(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_3); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_3);  CHKERRQ(ierr);
          // Second layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_3);  CHKERRQ(ierr);
          // Third layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_3);  CHKERRQ(ierr);
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_3);  CHKERRQ(ierr);
          }
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_3); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta_3rd_Ply",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"Solving the first order derivative of Theta_3"<<endl;
        }
        else if (idx_disp == 94) { // fourth layer angle is a random variable
          // Initiate the involved parameters to zero
          // cout<<"\n\nThe angle case Theta_4 starts from here!!!"<<endl;
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          Dmat_1st_Ply_r_Theta_4.resize(6,6);   Dmat_1st_Ply_r_Theta_4.clear();
          // Second-layer
          Dmat_2nd_Ply_r_Theta_4.resize(6,6);   Dmat_2nd_Ply_r_Theta_4.clear();
          // Third-layer
          Dmat_3rd_Ply_r_Theta_4.resize(6,6);   Dmat_3rd_Ply_r_Theta_4.clear();
          // Fourth layer
          theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_4th_Ply_r_Theta_4.resize(6,6);   Dmat_4th_Ply_r_Theta_4.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_4th_Ply_r_Theta_4); CHKERRQ(ierr);
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta_4,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta_4,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta_4,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta_4,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_4(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_4); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_4);  CHKERRQ(ierr);
          // Second layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_4);  CHKERRQ(ierr);
          // Third layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_4);  CHKERRQ(ierr);
          // Fourth layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_4);  CHKERRQ(ierr);
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_4); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta_4th_Ply",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"Solving the first order derivative of Theta_4"<<endl;
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
  };
}

#endif //__FE2_MACRO_SOLVER_HPP