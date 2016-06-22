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

#ifndef __SINGLE_FE_SOLVER_LAMINATE_HPP
#define __SINGLE_FE_SOLVER_LAMINATE_HPP

namespace MoFEM {
  
  /*****************************************************************************
   *                                                                           *
   *                                                                           *
   *                         MULTIPLE LAYERS LAMINATE                          *
   *                                                                           *
   *                                                                           *
   ****************************************************************************/
  
  struct Single_FE_Solver_Laminate {
    
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
    
    PetscErrorCode Dmat_Transformation(double theta, ublas::matrix<FieldData> Dmat_123, ublas::matrix<FieldData> &Dmat_xyz) {
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
    
    PetscErrorCode Dmat_Transformation_r_Theta(double theta, ublas::matrix<FieldData> Dmat_123, ublas::matrix<FieldData> &Dmat_xyz) {
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
    //  B.VI. SOLUTION PHASE:
    //        Solve Macroscale FE equation
    //
    // =========================================================================
    
    PetscErrorCode Macro_FE_REL(FieldInterface &m_field_Macro,
                                int &nvars, int &nders,
                                vector<string> &stochastic_fields,
                                ublas::vector<double> TheVariables,
                                int num_rvars,
                                vector<string> vars_name,
                                ublas::vector<double> PlyAngle,
                                int ExaminedPly,
                                PetscInt NO_Layers,
                                ublas::vector<ublas::matrix <double> > &TheStress) {
      
      PetscFunctionBegin;
      
      //ErrorCode rval;
      PetscErrorCode ierr;
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      //create matrices
      
      Vec F, D;
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&F); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&D); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat A;
      ierr = m_field_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_MACRO",&A); CHKERRQ(ierr);
      
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
      
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_Macro,BLOCKSET,it)) {
        cout << endl << *it << endl;
        
        //Get block name
        string name = it->get_name();
//        // ---------------------------------
//        //
//        // Case 1: Isotropic material
//        //
//        // ---------------------------------
//        if (name.compare(0,13,"MAT_ISOTROPIC") == 0) {
//          Mat_Elastic mydata;
//          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
//          
//          string ParameterName;
//          for (int ii=1;ii<=num_rvars;ii++) {
//            ParameterName = vars_name[ii];
//            
//            if (ParameterName.compare(0,2,"Ef") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
//              mydata.data.Young = matprop(ii-1);
//            }
//            else if (ParameterName.compare(0,3,"NUf") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
//              mydata.data.Poisson = matprop(ii-1);
//            }
//            ParameterName.clear();
//          }
//          
//          /*
//           mydata.data.Young   = matprop(0);
//           mydata.data.Poisson = matprop(1);
//           */
//          YoungModulus_Fibre = mydata.data.Young;
//          PoissonRatio_Fibre = mydata.data.Poisson;
//          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
//          cout << "Fibre material:\n" << mydata;
//        }
        
        // ---------------------------------
        //
        // Modify fibre material properties
        // Case 2: Transversely isotropic material
        //
        // ---------------------------------
        
        if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
          Mat_Elastic_TransIso mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int i=1;i<=num_rvars;i++) {
            ParameterName = vars_name[i];
            cout<<"the variable name is "<<vars_name[i]<<endl;
            if (ParameterName.compare(0,2,"Ez") == 0) {
              mydata.data.Youngz = TheVariables(i-1);
            }
            else if (ParameterName.compare(0,2,"Ep") == 0) {
              mydata.data.Youngp = TheVariables(i-1);
            }
            else if (ParameterName.compare(0,3,"NUp") == 0) {
              mydata.data.Poissonp = TheVariables(i-1);
            }
            else if (ParameterName.compare(0,3,"NUz") == 0) {
              mydata.data.Poissonpz = TheVariables(i-1);
            }
            else if (ParameterName.compare(0,3,"Gzp") == 0) {
              mydata.data.Shearzp = TheVariables(i-1);
            }
            else if (ParameterName.compare(0,2,"Ef") == 0) {
              mydata.data.Youngz = TheVariables(i-1);
              mydata.data.Youngp = TheVariables(i-1);
            }
            else if (ParameterName.compare(0,3,"NUf") == 0) {
              mydata.data.Poissonp  = TheVariables(i-1);
              mydata.data.Poissonpz = TheVariables(i-1);
            }
            else if ((ParameterName.compare(0,3,"NUf") == 0) || (ParameterName.compare(0,2,"Ef") == 0)) {
              mydata.data.Shearzp = mydata.data.Youngz/(2*(1+mydata.data.Poissonp));
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
      
      Projection10NodeCoordsOnField ent_method_material_Macro(m_field_Macro,"MESH_NODE_POSITIONS");
      ierr = m_field_Macro.loop_dofs("MESH_NODE_POSITIONS",ent_method_material_Macro); CHKERRQ(ierr);
      
      //Assemble F and A
      DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field_Macro,"DISP_MACRO",A,D,F);
      
      /*****************************************************************************
       *
       * Read the saved Dmat mechancial (from the computational homgenisaiton of the 0deg RVE)
       *
       ****************************************************************************/
      
      double theta;
      double AxVector[3] = {0,0,1.0};
      double AxAngle_1st[1] = {0}, AxAngle_2nd[1] = {0}, AxAngle_3rd[1] = {0}, AxAngle_4th[1] = {0};
      
      // First-layer
      AxAngle_1st[0] = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
      // Second-layer
      if (NO_Layers > 1) {
        AxAngle_2nd[0] = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
      }
      // Third-layer
      if (NO_Layers > 2) {
        AxAngle_3rd[0] = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
      }
      // Fourth-layer
      if (NO_Layers > 3) {
        AxAngle_4th[0] = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
      }
      
      TranIsotropicPlyElasticFEMethod my_fe_1st_ply(m_field_Macro,A,D,F,1,AxVector,AxAngle_1st,"DISP_MACRO");
      TranIsotropicPlyElasticFEMethod my_fe_2nd_ply(m_field_Macro,A,D,F,1,AxVector,AxAngle_2nd,"DISP_MACRO");
      TranIsotropicPlyElasticFEMethod my_fe_3rd_ply(m_field_Macro,A,D,F,1,AxVector,AxAngle_3rd,"DISP_MACRO");
      TranIsotropicPlyElasticFEMethod my_fe_4th_ply(m_field_Macro,A,D,F,1,AxVector,AxAngle_4th,"DISP_MACRO");
      
      
      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = MatZeroEntries(A); CHKERRQ(ierr);
      
      ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

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
      //MatView(A,PETSC_VIEWER_STDOUT_WORLD);
      
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
      
      ublas::matrix<double> StressGP(3,3);
      double Angle_REL;
      ublas::matrix<double> Dmat_REL; Dmat_REL.clear();
      ublas::matrix<double> TrpMatrixStress; TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
      
//      Angle_REL = -1*PlyAngle(0)*(M_PI/180.0);
//      cout<<"\nNumber of layers: "<<NO_Layers<<endl;
//      cout<<"\n First layer"<<endl;
//      cout<<"\n Angle: "<<AxAngle_1st[0]<<endl;
//      cout<<"\nDmat in 123: "<<my_fe_1st_ply.Dmat<<endl;
//      cout<<"\nDmat in xyz: "<<my_fe_1st_ply.Dmat_xyz<<endl;
//      ierr = Dmat_Transformation(Angle_REL, my_fe_1st_ply.Dmat, Dmat_REL); CHKERRQ(ierr);
//      cout<<"\nDmat in xyz: "<<Dmat_REL<<endl;
//      
//      cout<<"\n Second layer"<<endl;
//      cout<<"\n Angle: "<<AxAngle_2nd[0]<<endl;
//      Angle_REL = -1*PlyAngle(1)*(M_PI/180.0);
//      cout<<"\nDmat in 123: "<<my_fe_2nd_ply.Dmat<<endl;
//      cout<<"\nDmat in xyz: "<<my_fe_2nd_ply.Dmat_xyz<<endl;
//      ierr = Dmat_Transformation(Angle_REL, my_fe_2nd_ply.Dmat, Dmat_REL); CHKERRQ(ierr);
//      cout<<"\nDmat in xyz: "<<Dmat_REL<<endl;
//      
//      cout<<"\n Third layer"<<endl;
//      cout<<"\n Angle: "<<AxAngle_3rd[0]<<endl;
//      Angle_REL = -1*PlyAngle(2)*(M_PI/180.0);
//      cout<<"\nDmat in 123: "<<my_fe_3rd_ply.Dmat<<endl;
//      cout<<"\nDmat in xyz: "<<my_fe_3rd_ply.Dmat_xyz<<endl;
//      ierr = Dmat_Transformation(Angle_REL, my_fe_3rd_ply.Dmat, Dmat_REL); CHKERRQ(ierr);
//      cout<<"\nDmat in xyz: "<<Dmat_REL<<endl;
//      
//      cout<<"\n Fourth layer"<<endl;
//      cout<<"\n Angle: "<<AxAngle_4th[0]<<endl;
//      Angle_REL = -1*PlyAngle(3)*(M_PI/180.0);
//      cout<<"\nDmat in 123: "<<my_fe_4th_ply.Dmat<<endl;
//      cout<<"\nDmat in xyz: "<<my_fe_4th_ply.Dmat_xyz<<endl;
//      ierr = Dmat_Transformation(Angle_REL, my_fe_4th_ply.Dmat, Dmat_REL); CHKERRQ(ierr);
//      cout<<"\nDmat in xyz: "<<Dmat_REL<<endl;
      
      Angle_REL = PlyAngle(ExaminedPly - 1)*(M_PI/180.0);
      
      switch (ExaminedPly) {
        case 1: { Dmat_REL = my_fe_1st_ply.Dmat_xyz; TrpMatrixStress = my_fe_1st_ply.T_Stress; break; }
        case 2: { Dmat_REL = my_fe_2nd_ply.Dmat_xyz; TrpMatrixStress = my_fe_2nd_ply.T_Stress; break; }
        case 3: { Dmat_REL = my_fe_3rd_ply.Dmat_xyz; TrpMatrixStress = my_fe_3rd_ply.T_Stress; break; }
        case 4: { Dmat_REL = my_fe_4th_ply.Dmat_xyz; TrpMatrixStress = my_fe_4th_ply.T_Stress; break; }
      }
      cout<<"\n The Dmat "<<Dmat_REL<<endl;
      FE2_PostProcStressForReliability_Zeroth Calc_Stress(m_field_Macro,"DISP_MACRO",Dmat_REL);
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress); CHKERRQ(ierr);
//      StressGP.clear(); REL_Stress_Transformation(Angle_REL, Calc_Stress.StressGP, StressGP);
//      TheStress(0).resize(3,3); TheStress(0).clear(); TheStress(0) = StressGP;
//      cout<<"\nStress at GP in 123 in side FE2: "<<setprecision(15)<<StressGP<<endl;
//      cout<<"\nStress at GP in xyz in side FE2: "<<setprecision(15)<<Calc_Stress.StressGP<<endl;
//      
      ublas::vector<double> Stress_VectorNotation_xyz(6); Stress_VectorNotation_xyz.clear();
      Stress_VectorNotation_xyz(0) = Calc_Stress.StressGP(0,0);
      Stress_VectorNotation_xyz(1) = Calc_Stress.StressGP(1,1);
      Stress_VectorNotation_xyz(2) = Calc_Stress.StressGP(2,2);
      Stress_VectorNotation_xyz(3) = Calc_Stress.StressGP(0,1);
      Stress_VectorNotation_xyz(4) = Calc_Stress.StressGP(1,2);
      Stress_VectorNotation_xyz(5) = Calc_Stress.StressGP(2,0);
      
      ublas::vector<double> Stress_VectorNotation_123 = prod(TrpMatrixStress, Stress_VectorNotation_xyz);
      ublas::matrix<double> Stress_123(3,3); Stress_123.clear();
      Stress_123(0,0) = Stress_VectorNotation_123(0);
      Stress_123(1,1) = Stress_VectorNotation_123(1);
      Stress_123(2,2) = Stress_VectorNotation_123(2);
      Stress_123(0,1) = Stress_123(1,0) = Stress_VectorNotation_123(3);
      Stress_123(1,2) = Stress_123(2,1) = Stress_VectorNotation_123(4);
      Stress_123(2,0) = Stress_123(0,2) = Stress_VectorNotation_123(5);
      
      TheStress(0).resize(3,3); TheStress(0).clear(); TheStress(0) = Stress_123;
      cout<<"\nStress at GP in 123 in side FE2: "<<setprecision(15)<<Stress_123<<endl;
      

      /*************************************************************************
       *
       *  3. SOLVE THE FIRST-ORDER AND THE SECOND-ORDER FE EQUILIBRIUM EQUATION
       *     1st order-[K][U_r] = -[K_r][U}
       *     2nd order-[K][U_rs] = -[K_rs][U]-2[K_r][U_s]
       *
       ************************************************************************/
      int idx_disp; // index of displacement field in the field name vector <stochastic_fields>
      string VariableName;
      ublas::matrix<double> Dmat_REL_r;
      ublas::matrix<double> StressGP_r(3,3);
      for (int ii=1; ii<=num_rvars; ii++) {
        
        Vec dF, dD;
        ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&dF); CHKERRQ(ierr);
        ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&dD); CHKERRQ(ierr);
        idx_disp = 0; VariableName = vars_name[ii];
        Dmat_REL_r.resize(6,6); Dmat_REL_r.clear();
        
        ierr = VecZeroEntries(dD); CHKERRQ(ierr);
        ierr = VecZeroEntries(dF); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

        /***********************************************************************
         *
         * 3.1. Case 1: Material properties are  treated as random variables
         *
         **********************************************************************/
       
        if (VariableName.compare(0,2,"Ez") == 0) {// due to Young's modulus of matrix (Em)
          idx_disp = 1;
          
          Single_FE_Rhs_r_PSFEM my_fe_r_1st_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_1st,"DISP_MACRO","YoungZ","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_2nd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_2nd,"DISP_MACRO","YoungZ","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_3rd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_3rd,"DISP_MACRO","YoungZ","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_4th_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_4th,"DISP_MACRO","YoungZ","trans");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe_r_1st_ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe_r_2nd_ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe_r_3rd_ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe_r_4th_ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r); CHKERRQ(ierr);
          // Examined ply
          switch (ExaminedPly) {
            case 1: { Dmat_REL_r = my_fe_r_1st_ply.Dmat_r_xyz; break; }
            case 2: { Dmat_REL_r = my_fe_r_2nd_ply.Dmat_r_xyz; break; }
            case 3: { Dmat_REL_r = my_fe_r_3rd_ply.Dmat_r_xyz; break; }
            case 4: { Dmat_REL_r = my_fe_r_4th_ply.Dmat_r_xyz; break; }
          }
        }
        else if (VariableName.compare(0,3,"NUp") == 0) { // due to Poisson's ratio of matrix (NUm)
          idx_disp = 1;
          
          Single_FE_Rhs_r_PSFEM my_fe_r_1st_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_1st,"DISP_MACRO","PoissonP","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_2nd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_2nd,"DISP_MACRO","PoissonP","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_3rd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_3rd,"DISP_MACRO","PoissonP","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_4th_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_4th,"DISP_MACRO","PoissonP","trans");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUm(m_field_Macro,"DISP_MACRO",A,dD,dF);
          
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe_r_1st_ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe_r_2nd_ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe_r_3rd_ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe_r_4th_ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
          
          // Examined ply
  
          switch (ExaminedPly) {
            case 1: { Dmat_REL_r = my_fe_r_1st_ply.Dmat_r_xyz; break; }
            case 2: { Dmat_REL_r = my_fe_r_2nd_ply.Dmat_r_xyz; break; }
            case 3: { Dmat_REL_r = my_fe_r_3rd_ply.Dmat_r_xyz; break; }
            case 4: { Dmat_REL_r = my_fe_r_4th_ply.Dmat_r_xyz; break; }
          }
        }
        else if (VariableName.compare(0,3,"NUz") == 0) {// due to axial Poisson's ratio of fibre (NUpz)
          idx_disp = 1;
          
          Single_FE_Rhs_r_PSFEM my_fe_r_1st_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_1st,"DISP_MACRO","PoissonZ","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_2nd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_2nd,"DISP_MACRO","PoissonZ","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_3rd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_3rd,"DISP_MACRO","PoissonZ","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_4th_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_4th,"DISP_MACRO","PoissonZ","trans");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUpz(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe_r_1st_ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe_r_2nd_ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe_r_3rd_ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe_r_4th_ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
          
          // Examined ply
          switch (ExaminedPly) {
            case 1: { Dmat_REL_r = my_fe_r_1st_ply.Dmat_r_xyz; break; }
            case 2: { Dmat_REL_r = my_fe_r_2nd_ply.Dmat_r_xyz; break; }
            case 3: { Dmat_REL_r = my_fe_r_3rd_ply.Dmat_r_xyz; break; }
            case 4: { Dmat_REL_r = my_fe_r_4th_ply.Dmat_r_xyz; break; }
          }
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {// due to transversal modulus of fibre (Ep)
          idx_disp = 1;
          
          Single_FE_Rhs_r_PSFEM my_fe_r_1st_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_1st,"DISP_MACRO","YoungP","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_2nd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_2nd,"DISP_MACRO","YoungP","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_3rd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_3rd,"DISP_MACRO","YoungP","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_4th_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_4th,"DISP_MACRO","YoungP","trans");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ep(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe_r_1st_ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe_r_2nd_ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe_r_3rd_ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe_r_4th_ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
          
          // Examined ply
          switch (ExaminedPly) {
            case 1: { Dmat_REL_r = my_fe_r_1st_ply.Dmat_r_xyz; break; }
            case 2: { Dmat_REL_r = my_fe_r_2nd_ply.Dmat_r_xyz; break; }
            case 3: { Dmat_REL_r = my_fe_r_3rd_ply.Dmat_r_xyz; break; }
            case 4: { Dmat_REL_r = my_fe_r_4th_ply.Dmat_r_xyz; break; }
          }
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {// due to shear modulus of fibre (Gzp)
          idx_disp = 1;
          
          Single_FE_Rhs_r_PSFEM my_fe_r_1st_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_1st,"DISP_MACRO","ShearZP","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_2nd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_2nd,"DISP_MACRO","ShearZP","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_3rd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_3rd,"DISP_MACRO","ShearZP","trans");
          Single_FE_Rhs_r_PSFEM my_fe_r_4th_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_4th,"DISP_MACRO","ShearZP","trans");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Gzp(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe_r_1st_ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe_r_2nd_ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe_r_3rd_ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe_r_4th_ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
          
          // Examined ply
          switch (ExaminedPly) {
            case 1: { Dmat_REL_r = my_fe_r_1st_ply.Dmat_r_xyz; break; }
            case 2: { Dmat_REL_r = my_fe_r_2nd_ply.Dmat_r_xyz; break; }
            case 3: { Dmat_REL_r = my_fe_r_3rd_ply.Dmat_r_xyz; break; }
            case 4: { Dmat_REL_r = my_fe_r_4th_ply.Dmat_r_xyz; break; }
          }
        }
        else if (VariableName.compare(0,2,"Ef") == 0) {// due to Young's modulus of fibre (Ef) - isotropic
          idx_disp = 1;
          
          Single_FE_Rhs_r_PSFEM my_fe_r_1st_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_1st,"DISP_MACRO","Young","iso");
          Single_FE_Rhs_r_PSFEM my_fe_r_2nd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_2nd,"DISP_MACRO","Young","iso");
          Single_FE_Rhs_r_PSFEM my_fe_r_3rd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_3rd,"DISP_MACRO","Young","iso");
          Single_FE_Rhs_r_PSFEM my_fe_r_4th_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_4th,"DISP_MACRO","Young","iso");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ef(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ef); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe_r_1st_ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe_r_2nd_ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe_r_3rd_ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe_r_4th_ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ef); CHKERRQ(ierr);
          
          // Examined ply
          switch (ExaminedPly) {
            case 1: { Dmat_REL_r = my_fe_r_1st_ply.Dmat_r_xyz; break; }
            case 2: { Dmat_REL_r = my_fe_r_2nd_ply.Dmat_r_xyz; break; }
            case 3: { Dmat_REL_r = my_fe_r_3rd_ply.Dmat_r_xyz; break; }
            case 4: { Dmat_REL_r = my_fe_r_4th_ply.Dmat_r_xyz; break; }
          }
        }
        else if (VariableName.compare(0,3,"NUf") == 0) { // due to Poisson's ratio of fibre (NUf) - isotropic
          idx_disp = 1;
          
          Single_FE_Rhs_r_PSFEM my_fe_r_1st_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_1st,"DISP_MACRO","Poisson","iso");
          Single_FE_Rhs_r_PSFEM my_fe_r_2nd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_2nd,"DISP_MACRO","Poisson","iso");
          Single_FE_Rhs_r_PSFEM my_fe_r_3rd_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_3rd,"DISP_MACRO","Poisson","iso");
          Single_FE_Rhs_r_PSFEM my_fe_r_4th_ply(m_field_Macro,A,dD,dF,1,AxVector,AxAngle_4th,"DISP_MACRO","Poisson","iso");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUf(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUf); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe_r_1st_ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe_r_2nd_ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe_r_3rd_ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe_r_4th_ply);  CHKERRQ(ierr);
          }
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUf); CHKERRQ(ierr);
          
          // Examined ply
          switch (ExaminedPly) {
            case 1: { Dmat_REL_r = my_fe_r_1st_ply.Dmat_r_xyz; break; }
            case 2: { Dmat_REL_r = my_fe_r_2nd_ply.Dmat_r_xyz; break; }
            case 3: { Dmat_REL_r = my_fe_r_3rd_ply.Dmat_r_xyz; break; }
            case 4: { Dmat_REL_r = my_fe_r_4th_ply.Dmat_r_xyz; break; }
          }
        }
        /***********************************************************************
         *
         * 3.2. Case 2: Applied forces are treated as random variables
         *
         **********************************************************************/
        else if (VariableName.compare(0,5,"force") == 0) {
          idx_disp = 2;
          
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = MatZeroEntries(A); CHKERRQ(ierr);
          
          // Establish an object of elastic FE method
          TranIsotropicPlyElasticFEMethod my_fe_r_1st_ply(m_field_Macro,A,D,F,1,AxVector,AxAngle_1st,"DISP_MACRO");
          TranIsotropicPlyElasticFEMethod my_fe_r_2nd_ply(m_field_Macro,A,D,F,1,AxVector,AxAngle_2nd,"DISP_MACRO");
          TranIsotropicPlyElasticFEMethod my_fe_r_3rd_ply(m_field_Macro,A,D,F,1,AxVector,AxAngle_3rd,"DISP_MACRO");
          TranIsotropicPlyElasticFEMethod my_fe_r_4th_ply(m_field_Macro,A,D,F,1,AxVector,AxAngle_4th,"DISP_MACRO");
          
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
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe_r_1st_ply);     CHKERRQ(ierr);
          // Third layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe_r_2nd_ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe_r_3rd_ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe_r_4th_ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
          
          ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
          
        }
//        /***********************************************************************
//         *
//         * 3.3. Case 3: Ply orientation angle is treated as random variable
//         *
//         **********************************************************************/
//        else if (VariableName.compare(0,11,"orientation") == 0) {
//          idx_disp = 1;
//          // First-layer
//          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
//          Dmat_1st_Ply_r_Theta.resize(6,6);   Dmat_1st_Ply_r_Theta.clear();
//          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_1st_Ply_r_Theta); CHKERRQ(ierr);
//          // Second-layer
//          if (NO_Layers > 1) {
//            theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
//            Dmat_2nd_Ply_r_Theta.resize(6,6);   Dmat_2nd_Ply_r_Theta.clear();
//            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_2nd_Ply_r_Theta); CHKERRQ(ierr);
//          }
//          // Third-layer
//          if (NO_Layers > 2) {
//            theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
//            Dmat_3rd_Ply_r_Theta.resize(6,6);   Dmat_3rd_Ply_r_Theta.clear();
//            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_3rd_Ply_r_Theta); CHKERRQ(ierr);
//          }
//          // Fourth layer
//          if (NO_Layers > 3) {
//            theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
//            Dmat_4th_Ply_r_Theta.resize(6,6);   Dmat_4th_Ply_r_Theta.clear();
//            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_4th_Ply_r_Theta); CHKERRQ(ierr);
//          }
//          
//          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta,"DISP_MACRO");
//          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta,"DISP_MACRO");
//          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta,"DISP_MACRO");
//          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta,"DISP_MACRO");
//          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta(m_field_Macro,"DISP_MACRO",A,dD,dF);
//          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta); CHKERRQ(ierr);
//          // First layer
//          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Theta_1st_Ply);  CHKERRQ(ierr);
//          // Second layer
//          if (NO_Layers > 1) {
//            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Theta_2nd_Ply);  CHKERRQ(ierr);
//          }
//          // Third layer
//          if (NO_Layers > 2) {
//            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Theta_3rd_Ply);  CHKERRQ(ierr);
//          }
//          // Fourth layer
//          if (NO_Layers > 3) {
//            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Theta_4th_Ply);  CHKERRQ(ierr);
//          }
//          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta); CHKERRQ(ierr);
//          
//        }
        
        // post-processing
        ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
        
        ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);//ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        
        // Calculate stress for examined element
        if (idx_disp == 2) {
          Dmat_REL_r.clear();
        }

        FE2_PostProcStressForReliability_First Calc_Stress_r(m_field_Macro,"DISP_MACRO","DISP_MACRO_r",Dmat_REL,Dmat_REL_r);
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_r);  CHKERRQ(ierr);
//        StressGP_r.clear(); REL_Stress_Transformation(Angle_REL, Calc_Stress_r.StressGP_r, StressGP_r);
//        cout<<"Stress_r for "<<VariableName<<" is: "<<setprecision(15)<<StressGP_r<<endl;
//        TheStress(ii).resize(3,3); TheStress(ii).clear(); TheStress(ii) = StressGP_r;
        
        ublas::vector<double> Stress_VectorNotation_xyz(6); Stress_VectorNotation_xyz.clear();
        Stress_VectorNotation_xyz(0) = Calc_Stress_r.StressGP_r(0,0);
        Stress_VectorNotation_xyz(1) = Calc_Stress_r.StressGP_r(1,1);
        Stress_VectorNotation_xyz(2) = Calc_Stress_r.StressGP_r(2,2);
        Stress_VectorNotation_xyz(3) = Calc_Stress_r.StressGP_r(0,1);
        Stress_VectorNotation_xyz(4) = Calc_Stress_r.StressGP_r(1,2);
        Stress_VectorNotation_xyz(5) = Calc_Stress_r.StressGP_r(2,0);
        
        ublas::vector<double> Stress_VectorNotation_123 = prod(TrpMatrixStress, Stress_VectorNotation_xyz);
        ublas::matrix<double> Stress_123(3,3); Stress_123.clear();
        Stress_123(0,0) = Stress_VectorNotation_123(0);
        Stress_123(1,1) = Stress_VectorNotation_123(1);
        Stress_123(2,2) = Stress_VectorNotation_123(2);
        Stress_123(0,1) = Stress_123(1,0) = Stress_VectorNotation_123(3);
        Stress_123(1,2) = Stress_123(2,1) = Stress_VectorNotation_123(4);
        Stress_123(2,0) = Stress_123(0,2) = Stress_VectorNotation_123(5);
        
        TheStress(ii).resize(3,3); TheStress(ii).clear(); TheStress(ii) = Stress_123;
        cout<<"\nStress at GP in 123 in side FE2: "<<setprecision(15)<<Stress_123<<endl;
        
        // Reset parameters
        VariableName.clear();
        ierr = VecDestroy(&dF); CHKERRQ(ierr);
        ierr = VecDestroy(&dD); CHKERRQ(ierr);
      }

      /***************************************************************************
       *
       *  4. FINISH
       *
       **************************************************************************/
      
      //Destroy matrices
      ierr = VecDestroy(&F); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = MatDestroy(&A); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver_Macro); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }

  };
  
}

#endif //__SINGLE_FE_SOLVER_LAMINATE_HPP