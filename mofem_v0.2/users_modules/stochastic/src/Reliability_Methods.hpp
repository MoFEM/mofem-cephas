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

#ifndef __FE2_RELIABILITY_METHODS_HPP
#define __FE2_RELIABILITY_METHODS_HPP

namespace MoFEM {
  
  struct Reliability_Methods {
    /***************************************************************************
     *                                                                         *
     *                                                                         *
     *             FIRST-ORDER RELIABILITY METHOD - FORM                       *
     *                                                                         *
     *                                                                         *
     **************************************************************************/
    virtual PetscErrorCode FORM(FieldInterface &m_field_RVE,
                                FieldInterface &m_field_Macro,
                                int &nvars, int &nders,
                                vector<string> &stochastic_fields,
                                PetscInt NO_Layers,
                                Stochastic_Model probdata,
                                Reliability_Options ReliabOpt,
                                int FailureCriterion,
                                Reliability_Results &MsFE_Reliability_Results) {
      PetscFunctionBegin;
      
      //////////////////////////////////////////////////////////////////////////////
      //                                                                          //
      //        STEP 2: use probability transformation to obtain y in             //
      //                standard normal distribution space                        //
      //                y = norminv(F(x))                                         //
      //                                                                          //
      //////////////////////////////////////////////////////////////////////////////
      NatafTransformation my_nataf_transformation;
      ierr = my_nataf_transformation.ModCorrMat_Empirical(probdata.num_vars,
                                                          probdata.marg,
                                                          probdata.correlation,
                                                          probdata.mod_correlation); CHKERRQ(ierr);
      
      //
      // Perform Cholesky decomposition for the modified correlation matrix
      //  A = LL'
      //
      
      ublas::triangular_matrix<double, ublas::lower> Lmat(probdata.num_vars,probdata.num_vars);
      cholesky_decompose(probdata.mod_correlation,Lmat);
      probdata.Lo.resize(probdata.num_vars,probdata.num_vars);
      probdata.Lo = Lmat;
      
      // Compute the inverse of Lo
      probdata.inv_Lo.resize(probdata.num_vars,probdata.num_vars);
      bool singular = false;
      probdata.inv_Lo = gjinverse(probdata.Lo,singular);
      
      
      //////////////////////////////////////////////////////////////////////////////
      //                                                                          //
      //        STEP 3: select an initial checking point, x                       //
      //                                                                          //
      //////////////////////////////////////////////////////////////////////////////
      
      ublas::vector<double> x(probdata.num_vars);
      for (int i=0; i<probdata.num_vars; i++) {
        x(i) = probdata.marg(i,3);
      }
      
      ublas::vector<double> u;
      ierr = my_nataf_transformation.x_to_u(x,probdata.num_vars,probdata.marg,probdata.inv_Lo,u); CHKERRQ(ierr);
      
      //////////////////////////////////////////////////////////////////////////////
      //                                                                          //
      //        STEP 4: Perform iterative loop to find design point               //
      //                                                                          //
      //////////////////////////////////////////////////////////////////////////////
      
      // Declare matrix for stroring RVE constitutive matrix & its derivatives
      ublas::matrix<double> Dmat(6,6);           Dmat.clear();
      ublas::matrix<double> Dmat_r_Em(6,6);      Dmat_r_Em.clear();
      ublas::matrix<double> Dmat_r_NUm(6,6);     Dmat_r_NUm.clear();
      ublas::matrix<double> Dmat_r_Ep(6,6);      Dmat_r_Ep.clear();
      ublas::matrix<double> Dmat_r_Ez(6,6);      Dmat_r_Ez.clear();
      ublas::matrix<double> Dmat_r_NUp(6,6);     Dmat_r_NUp.clear();
      ublas::matrix<double> Dmat_r_NUpz(6,6);    Dmat_r_NUpz.clear();
      ublas::matrix<double> Dmat_r_Gzp(6,6);     Dmat_r_Gzp.clear();
      ublas::matrix<double> Dmat_r_Ef(6,6);      Dmat_r_NUpz.clear();
      ublas::matrix<double> Dmat_r_NUf(6,6);     Dmat_r_NUf.clear();
      ublas::matrix<double> Dmat_r_F(6,6);       Dmat_r_F.clear();
      ublas::matrix<double> Dmat_r_Theta(6,6);   Dmat_r_Theta.clear();
      ublas::matrix<double> Dmat_r_Theta_1(6,6); Dmat_r_Theta_1.clear();
      ublas::matrix<double> Dmat_r_Theta_2(6,6); Dmat_r_Theta_2.clear();
      ublas::matrix<double> Dmat_r_Theta_3(6,6); Dmat_r_Theta_3.clear();
      ublas::matrix<double> Dmat_r_Theta_4(6,6); Dmat_r_Theta_4.clear();
      
      // Declate matrix for storing stress matrix & its derivatives
      ublas::matrix<double> StressGP_Global(3,3);    StressGP_Global.clear();
      ublas::matrix<double> StressGP(3,3);           StressGP.clear();
      ublas::matrix<double> StressGP_r_Em(3,3);      StressGP_r_Em.clear();
      ublas::matrix<double> StressGP_r_NUm(3,3);     StressGP_r_NUm.clear();
      ublas::matrix<double> StressGP_r_Ep(3,3);      StressGP_r_Ep.clear();
      ublas::matrix<double> StressGP_r_Ez(3,3);      StressGP_r_Ez.clear();
      ublas::matrix<double> StressGP_r_NUp(3,3);     StressGP_r_NUp.clear();
      ublas::matrix<double> StressGP_r_NUpz(3,3);    StressGP_r_NUpz.clear();
      ublas::matrix<double> StressGP_r_Gzp(3,3);     StressGP_r_Gzp.clear();
      ublas::matrix<double> StressGP_r_Ef(3,3);      StressGP_r_Ef.clear();
      ublas::matrix<double> StressGP_r_NUf(3,3);     StressGP_r_NUf.clear();
      ublas::matrix<double> StressGP_r_F(3,3);       StressGP_r_F.clear();
      ublas::matrix<double> StressGP_r_Theta(3,3);   StressGP_r_Theta.clear();
      ublas::matrix<double> StressGP_r_Theta_1(3,3); StressGP_r_Theta_1.clear();
      ublas::matrix<double> StressGP_r_Theta_2(3,3); StressGP_r_Theta_2.clear();
      ublas::matrix<double> StressGP_r_Theta_3(3,3); StressGP_r_Theta_3.clear();
      ublas::matrix<double> StressGP_r_Theta_4(3,3); StressGP_r_Theta_4.clear();
      
      int    echo_flag = ReliabOpt.echo_flag;
      double e1        = ReliabOpt.e1;
      double e2        = ReliabOpt.e2;
      int    istep_max = ReliabOpt.istep_max;
      double step_code = ReliabOpt.step_code;
      int    beta_flag = ReliabOpt.Recorded_beta;
      
      //
      // Set parameters for the iterative loop
      //
      int    istep = 1;     // Initialize iterative counter
      int    conv_flag = 0; // Convergence is achieved when this flag is set to 1
      
      double val_G, val_G0;                                         // LSF for given inputs
      double step_size;                                             // Step size
      ublas::vector<double> grad_g(probdata.num_vars); // Gradient of LSF in x space
      ublas::vector<double> grad_G(probdata.num_vars); // Gradient of LSF in u space
      ublas::vector<double> alpha(probdata.num_vars);  // Direction cosine vector
      ublas::matrix<double> dudx(probdata.num_vars,probdata.num_vars);
      ublas::matrix<double> inv_dudx(probdata.num_vars,probdata.num_vars);
      ublas::vector<double> u_dir(probdata.num_vars);  // Direction
      ublas::vector<double> u_new(probdata.num_vars);  // New trial of checking point
      
      double theta_angle;
      ublas::vector<double> PlyAngle_new;
      PlyAngle_new = probdata.PlyAngle;
      
      //
      // Start iteration
      //
      
      clock_t rel_t1, rel_t2;
      double rel_calc_time;
      rel_t1 = clock();
      
      FE2_Macro_Solver_Laminate Solve_FE2_Problem;
      LSF_Composite_Lamina TheLSF;
      
      ofstream BetaFile;
      if (beta_flag == 1) {
        BetaFile.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//04_ReliabilityAnalysis//Result_Beta.txt",ofstream::out);
      }
      
      do {
        
        cout<<"\n\n*************************************************\n*\n";
        cout<<"*    This is "<<istep<<" step!"<<"\n*\n";
        cout<<"*************************************************\n";
        
        if (echo_flag) {
          cout<<"-------------------------------- \n";
          cout<<"Now carrying out iteration number: \t"<<istep<<endl;
        }
        
        
        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //        STEP 5: calculate Jacobian, J = dy/dx                           //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        
        //
        // Transformation from u to x space
        //
        double detj = 1.0;
        x.resize(probdata.num_vars); x.clear();cout<<"u value: "<<u<<endl;
        ierr = my_nataf_transformation.u_to_x(u,probdata.num_vars,probdata.marg,probdata.Lo,x,detj); CHKERRQ(ierr);
        
        // ----------------
        // Update ply angle
        // -----------------
        for (unsigned i=1; i<=x.size();i++) {
          if (probdata.NameVars[i].compare(0,11,"orientation") == 0) {
            cout<<"\nAngle "<<x(i-1)<<endl;
            for (unsigned j=0; j<probdata.PlyAngle.size(); j++) {
              cout << "The original ply angle "<<probdata.PlyAngle(j)<<"\t delta x: "<<x(i-1)<<endl;
              PlyAngle_new(j) = probdata.PlyAngle(j) + x(i-1);
              cout << "The modified ply angle "<<PlyAngle_new(j)<<endl;
            }
          }
          else if (probdata.NameVars[i].compare(0,6,"theta1") == 0) {
            PlyAngle_new(0) = x(i-1);
          }
          else if (probdata.NameVars[i].compare(0,6,"theta2") == 0) {
            PlyAngle_new(1) = x(i-1);
          }
          else if (probdata.NameVars[i].compare(0,6,"theta3") == 0) {
            PlyAngle_new(2) = x(i-1);
          }
          else if (probdata.NameVars[i].compare(0,6,"theta4") == 0) {
            PlyAngle_new(3) = x(i-1);
          }
        }
        //cout << "\n\nThe modified ply angle "<<PlyAngle_new<<endl;
        
        // Jacobian
        dudx.resize(probdata.num_vars,probdata.num_vars); dudx.clear();
        ierr = my_nataf_transformation.Jacobian_u_x(x,u,probdata.num_vars,probdata.marg,probdata.Lo,probdata.inv_Lo,dudx); CHKERRQ(ierr);
        
        inv_dudx.resize(probdata.num_vars,probdata.num_vars); inv_dudx.clear();
        inv_dudx = gjinverse(dudx,singular);
        
        //
        // Evaluate limit-state function and its gradient
        //
        
        ierr = Solve_FE2_Problem.Micro_FE_Dmat(m_field_RVE, nvars, nders,
                                               stochastic_fields, x, probdata.num_vars,
                                               probdata.NameVars); CHKERRQ(ierr);
        ierr = Solve_FE2_Problem.Macro_FE_REL(m_field_Macro, nvars, nders,
                                              stochastic_fields, x, probdata.num_vars,
                                              probdata.NameVars,PlyAngle_new,//probdata.PlyAngle,
                                              NO_Layers); CHKERRQ(ierr);
        
        switch (probdata.ExaminedPly) {
          case 1:{cout<<"\n\nThe 1st layer is under examination.\n\n";
            //--------------------------
            // D matrix for the 1st ply
            //--------------------------
            Dmat             = Solve_FE2_Problem.Dmat_1st_Ply; cout<<"Dmat: "<<Dmat<<"\n\n";
            Dmat_r_Em        = Solve_FE2_Problem.Dmat_1st_Ply_r_Em; //cout<<"Dmat_r_Em: "<<Dmat_r_Em<<"\n\n";
            Dmat_r_NUm       = Solve_FE2_Problem.Dmat_1st_Ply_r_NUm;
            Dmat_r_Ep        = Solve_FE2_Problem.Dmat_1st_Ply_r_Ep;
            Dmat_r_Ez        = Solve_FE2_Problem.Dmat_1st_Ply_r_Ez; //cout<<"Dmat_r_Ez: "<<Dmat_r_Ez<<"\n\n";
            Dmat_r_NUp       = Solve_FE2_Problem.Dmat_1st_Ply_r_NUp;
            Dmat_r_NUpz      = Solve_FE2_Problem.Dmat_1st_Ply_r_NUpz;
            Dmat_r_Gzp       = Solve_FE2_Problem.Dmat_1st_Ply_r_Gzp;
            Dmat_r_Ef        = Solve_FE2_Problem.Dmat_1st_Ply_r_Ef;
            Dmat_r_NUf       = Solve_FE2_Problem.Dmat_1st_Ply_r_NUf;
            Dmat_r_Theta     = Solve_FE2_Problem.Dmat_1st_Ply_r_Theta;
            Dmat_r_Theta_1   = Solve_FE2_Problem.Dmat_1st_Ply_r_Theta_1;
            Dmat_r_Theta_2   = Solve_FE2_Problem.Dmat_1st_Ply_r_Theta_2;
            Dmat_r_Theta_3   = Solve_FE2_Problem.Dmat_1st_Ply_r_Theta_3;
            Dmat_r_Theta_4   = Solve_FE2_Problem.Dmat_1st_Ply_r_Theta_4;
            break;
          }
          case 2: {cout<<"\n\nThe 2nd layer is under examination.\n\n";
            //--------------------------
            // D matrix for the 2nd ply
            //--------------------------
            Dmat             = Solve_FE2_Problem.Dmat_2nd_Ply; //cout<<"Dmat: "<<Dmat<<"\n\n";
            Dmat_r_Em        = Solve_FE2_Problem.Dmat_2nd_Ply_r_Em; //cout<<"Dmat_r_Em: "<<Dmat_r_Em<<"\n\n";
            Dmat_r_NUm       = Solve_FE2_Problem.Dmat_2nd_Ply_r_NUm;
            Dmat_r_Ep        = Solve_FE2_Problem.Dmat_2nd_Ply_r_Ep;
            Dmat_r_Ez        = Solve_FE2_Problem.Dmat_2nd_Ply_r_Ez; //cout<<"Dmat_r_Ez: "<<Dmat_r_Ez<<"\n\n";
            Dmat_r_NUp       = Solve_FE2_Problem.Dmat_2nd_Ply_r_NUp;
            Dmat_r_NUpz      = Solve_FE2_Problem.Dmat_2nd_Ply_r_NUpz;
            Dmat_r_Gzp       = Solve_FE2_Problem.Dmat_2nd_Ply_r_Gzp;
            Dmat_r_Ef        = Solve_FE2_Problem.Dmat_2nd_Ply_r_Ef;
            Dmat_r_NUf       = Solve_FE2_Problem.Dmat_2nd_Ply_r_NUf;
            Dmat_r_Theta     = Solve_FE2_Problem.Dmat_2nd_Ply_r_Theta;
            Dmat_r_Theta_1   = Solve_FE2_Problem.Dmat_2nd_Ply_r_Theta_1;
            Dmat_r_Theta_2   = Solve_FE2_Problem.Dmat_2nd_Ply_r_Theta_2;
            Dmat_r_Theta_3   = Solve_FE2_Problem.Dmat_2nd_Ply_r_Theta_3;
            Dmat_r_Theta_4   = Solve_FE2_Problem.Dmat_2nd_Ply_r_Theta_4;
            break;
          }
          case 3: {cout<<"\n\nThe 3rd layer is under examination.\n\n";
            //--------------------------
            // D matrix for the 3rd ply
            //--------------------------
            Dmat             = Solve_FE2_Problem.Dmat_3rd_Ply; //cout<<"Dmat: "<<Dmat<<"\n\n";
            Dmat_r_Em        = Solve_FE2_Problem.Dmat_3rd_Ply_r_Em; //cout<<"Dmat_r_Em: "<<Dmat_r_Em<<"\n\n";
            Dmat_r_NUm       = Solve_FE2_Problem.Dmat_3rd_Ply_r_NUm;
            Dmat_r_Ep        = Solve_FE2_Problem.Dmat_3rd_Ply_r_Ep;
            Dmat_r_Ez        = Solve_FE2_Problem.Dmat_3rd_Ply_r_Ez; //cout<<"Dmat_r_Ez: "<<Dmat_r_Ez<<"\n\n";
            Dmat_r_NUp       = Solve_FE2_Problem.Dmat_3rd_Ply_r_NUp;
            Dmat_r_NUpz      = Solve_FE2_Problem.Dmat_3rd_Ply_r_NUpz;
            Dmat_r_Gzp       = Solve_FE2_Problem.Dmat_3rd_Ply_r_Gzp;
            Dmat_r_Ef        = Solve_FE2_Problem.Dmat_3rd_Ply_r_Ef;
            Dmat_r_NUf       = Solve_FE2_Problem.Dmat_3rd_Ply_r_NUf;
            Dmat_r_Theta     = Solve_FE2_Problem.Dmat_3rd_Ply_r_Theta;
            Dmat_r_Theta_1   = Solve_FE2_Problem.Dmat_3rd_Ply_r_Theta_1;
            Dmat_r_Theta_2   = Solve_FE2_Problem.Dmat_3rd_Ply_r_Theta_2;
            Dmat_r_Theta_3   = Solve_FE2_Problem.Dmat_3rd_Ply_r_Theta_3;
            Dmat_r_Theta_4   = Solve_FE2_Problem.Dmat_3rd_Ply_r_Theta_4;
            break;
          }
          case 4: {cout<<"\n\nThe 4th layer is under examination.\n\n";
            //--------------------------
            // D matrix for the 4th ply
            //--------------------------
            Dmat             = Solve_FE2_Problem.Dmat_4th_Ply; //cout<<"Dmat: "<<Dmat<<"\n\n";
            Dmat_r_Em        = Solve_FE2_Problem.Dmat_4th_Ply_r_Em; //cout<<"Dmat_r_Em: "<<Dmat_r_Em<<"\n\n";
            Dmat_r_NUm       = Solve_FE2_Problem.Dmat_4th_Ply_r_NUm;
            Dmat_r_Ep        = Solve_FE2_Problem.Dmat_4th_Ply_r_Ep;
            Dmat_r_Ez        = Solve_FE2_Problem.Dmat_4th_Ply_r_Ez; //cout<<"Dmat_r_Ez: "<<Dmat_r_Ez<<"\n\n";
            Dmat_r_NUp       = Solve_FE2_Problem.Dmat_4th_Ply_r_NUp;
            Dmat_r_NUpz      = Solve_FE2_Problem.Dmat_4th_Ply_r_NUpz;
            Dmat_r_Gzp       = Solve_FE2_Problem.Dmat_4th_Ply_r_Gzp;
            Dmat_r_Ef        = Solve_FE2_Problem.Dmat_4th_Ply_r_Ef;
            Dmat_r_NUf       = Solve_FE2_Problem.Dmat_4th_Ply_r_NUf;
            Dmat_r_Theta     = Solve_FE2_Problem.Dmat_4th_Ply_r_Theta;
            Dmat_r_Theta_1   = Solve_FE2_Problem.Dmat_4th_Ply_r_Theta_1;
            Dmat_r_Theta_2   = Solve_FE2_Problem.Dmat_4th_Ply_r_Theta_2;
            Dmat_r_Theta_3   = Solve_FE2_Problem.Dmat_4th_Ply_r_Theta_3;
            Dmat_r_Theta_4   = Solve_FE2_Problem.Dmat_4th_Ply_r_Theta_4;
            break;
          }
        }
        
        Dmat_r_F.resize(6,6); Dmat_r_F.clear();
        
        theta_angle = PlyAngle_new(probdata.ExaminedPly - 1)*(M_PI/180.0);
        
        // Calculate the zeroth-order stress at Gauss points for specific element(s)
        FE2_PostProcStressForReliability_Zeroth Calc_Stress(m_field_Macro,"DISP_MACRO",Dmat);
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress); CHKERRQ(ierr);
        StressGP.clear(); REL_Stress_Transformation(theta_angle, Calc_Stress.StressGP, StressGP);
        //cout<<"Stress at GP in xyz: "<<Calc_Stress.StressGP<<endl;
        //cout<<"Stress at GP in 123: "<<StressGP<<endl;
        
        // Calculate the first-order partial derivative stress at Gauss points for specific element(s)
        for (unsigned i=1; i<=x.size(); i++) {
          cout<<"The variable name is "<<probdata.NameVars[i]<<endl;
          if (probdata.NameVars[i].compare(0,2,"Em") == 0) {
            // w.r.t. Em
            FE2_PostProcStressForReliability_First Calc_Stress_Em(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Em",Dmat,Dmat_r_Em);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Em);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_Em.StressGP_r;
            StressGP_r_Em.clear(); REL_Stress_Transformation(theta_angle, StressGP_Global, StressGP_r_Em);
            cout<<"Stress_r_Em at GP: "<<StressGP_r_Em<<endl;
          }
          else if (probdata.NameVars[i].compare(0,3,"NUm") == 0) {
            // w.r.t. NUm
            FE2_PostProcStressForReliability_First Calc_Stress_NUm(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_NUm",Dmat,Dmat_r_NUm);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_NUm);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_NUm.StressGP_r;
            StressGP_r_NUm.clear(); REL_Stress_Transformation(theta_angle, StressGP_Global, StressGP_r_NUm);
            cout<<"Stress_r_NUm at GP: "<<StressGP_r_NUm<<endl;
          }
          else if (probdata.NameVars[i].compare(0,3,"NUp") == 0) {
            // w.r.t. NUp
            FE2_PostProcStressForReliability_First Calc_Stress_NUp(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_NUp",Dmat,Dmat_r_NUp);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_NUp);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_NUp.StressGP_r;
            StressGP_r_NUp.clear(); REL_Stress_Transformation(theta_angle, StressGP_Global, StressGP_r_NUp);
            cout<<"Stress_r_NUp at GP: "<<StressGP_r_NUp<<endl;
          }
          else if (probdata.NameVars[i].compare(0,3,"NUz") == 0) {
            // w.r.t. NUpz
            FE2_PostProcStressForReliability_First Calc_Stress_NUpz(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_NUpz",Dmat,Dmat_r_NUpz);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_NUpz);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_NUpz.StressGP_r;
            StressGP_r_NUpz.clear(); REL_Stress_Transformation(theta_angle, StressGP_Global, StressGP_r_NUpz);
            cout<<"Stress_r_NUpz at GP: "<<StressGP_r_NUpz<<endl;
          }
          else if (probdata.NameVars[i].compare(0,2,"Ep") == 0) {
            // w.r.t. Ep
            FE2_PostProcStressForReliability_First Calc_Stress_Ep(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Ep",Dmat,Dmat_r_Ep);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Ep);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_Ep.StressGP_r;
            StressGP_r_Ep.clear(); REL_Stress_Transformation(theta_angle, StressGP_Global, StressGP_r_Ep);
            cout<<"Stress_r_Ep at GP: "<<StressGP_r_Ep<<endl;
          }
          else if (probdata.NameVars[i].compare(0,2,"Ez") == 0) {
            // w.r.t. Ez
            FE2_PostProcStressForReliability_First Calc_Stress_Ez(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Ez",Dmat,Dmat_r_Ez);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Ez);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_Ez.StressGP_r;
            StressGP_r_Ez.clear(); REL_Stress_Transformation(theta_angle, StressGP_Global, StressGP_r_Ez);
            cout<<"Stress_r_Ez at GP: "<<StressGP_r_Ez<<endl;
          }
          else if (probdata.NameVars[i].compare(0,3,"Gzp") == 0) {
            // w.r.t. Gzp
            FE2_PostProcStressForReliability_First Calc_Stress_Gzp(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Gzp",Dmat,Dmat_r_Gzp);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Gzp);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_Gzp.StressGP_r;
            StressGP_r_Gzp.clear(); REL_Stress_Transformation(theta_angle, StressGP_Global, StressGP_r_Gzp);
            cout<<"Stress_r_Gzp at GP: "<<StressGP_r_Gzp<<endl;
          }
          else if (probdata.NameVars[i].compare(0,2,"Ef") == 0) {
            // w.r.t. Ef - fibre with isotropic material
            FE2_PostProcStressForReliability_First Calc_Stress_Ef(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Ef",Dmat,Dmat_r_Ef);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Ef);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_Ef.StressGP_r;
            StressGP_r_Ef.clear(); REL_Stress_Transformation(theta_angle, StressGP_Global, StressGP_r_Ef);
            cout<<"Stress_r_Ef at GP: "<<StressGP_r_Ef<<endl;
          }
          else if (probdata.NameVars[i].compare(0,3,"NUf") == 0) {
            // w.r.t. NUf - fibre with isotropic material
            FE2_PostProcStressForReliability_First Calc_Stress_NUf(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_NUf",Dmat,Dmat_r_NUf);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_NUf);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_NUf.StressGP_r;
            StressGP_r_NUf.clear(); REL_Stress_Transformation(theta_angle, StressGP_Global, StressGP_r_NUf);
            cout<<"Stress_r_NUf at GP: "<<StressGP_r_NUf<<endl;
          }
          else if (probdata.NameVars[i].compare(0,5,"force") == 0) {
            // w.r.t. F
            // with respect to F
            FE2_PostProcStressForReliability_First Calc_Stress_F(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_F",Dmat,Dmat_r_F);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_F);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_F.StressGP_r;
            StressGP_r_F.clear(); REL_Stress_Transformation(theta_angle,StressGP_Global,StressGP_r_F);
            cout<<"Stress_r_F at GP: "<<StressGP_r_F<<endl;
          }
          else if (probdata.NameVars[i].compare(0,11,"orientation") == 0) {
            // w.r.t. ply angle
            FE2_PostProcStressForReliability_First Calc_Stress_Theta(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Theta",Dmat,Dmat_r_Theta);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Theta);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_Theta.StressGP_r;
            StressGP_r_Theta.clear(); REL_Stress_Transformation(theta_angle,StressGP_Global,StressGP_r_Theta);
            //cout<<"StressGP_r_Theta_xyz: "<<Calc_Stress_Theta.StressGP_r<<endl;
            cout<<"StressGP_r_Theta_123: "<<StressGP_r_Theta<<endl;
          }
          else if (probdata.NameVars[i].compare(0,6,"theta1") == 0) {
            // w.r.t. ply angle
            // 1st-layer
            FE2_PostProcStressForReliability_First Calc_Stress_Theta_1(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Theta_1st_Ply",Dmat,Dmat_r_Theta_1);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Theta_1);  CHKERRQ(ierr);
            StressGP_r_Theta_1.clear();
            if (probdata.ExaminedPly == 1) {
              REL_Stress_Transformation_Theta(theta_angle,
                                              Calc_Stress.StressGP,
                                              Calc_Stress_Theta_1.StressGP_r,
                                              StressGP_r_Theta_1);
            } else {
              REL_Stress_Transformation(theta_angle,Calc_Stress_Theta_1.StressGP_r,
                                        StressGP_r_Theta_1);
            }
          }
          else if (probdata.NameVars[i].compare(0,6,"theta2") == 0) {
            // w.r.t. ply angle
            // 2nd-layer
            FE2_PostProcStressForReliability_First Calc_Stress_Theta_2(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Theta_2nd_Ply",Dmat,Dmat_r_Theta_2);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Theta_2);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_Theta_2.StressGP_r;
            StressGP_r_Theta_2.clear();
            if (probdata.ExaminedPly == 2) {
              REL_Stress_Transformation_Theta(theta_angle,
                                              Calc_Stress.StressGP,
                                              Calc_Stress_Theta_2.StressGP_r,
                                              StressGP_r_Theta_2);
            } else {
              REL_Stress_Transformation(theta_angle,StressGP_Global,StressGP_r_Theta_2);
            }
          }
          else if (probdata.NameVars[i].compare(0,6,"theta3") == 0) {
            // w.r.t. ply angle
            // 3rd-layer
            FE2_PostProcStressForReliability_First Calc_Stress_Theta_3(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Theta_3rd_Ply",Dmat,Dmat_r_Theta_3);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Theta_3);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_Theta_3.StressGP_r;
            if (probdata.ExaminedPly == 3) {
              StressGP_r_Theta_3.clear(); REL_Stress_Transformation_Theta(theta_angle,
                                                                          Calc_Stress.StressGP,
                                                                          Calc_Stress_Theta_3.StressGP_r,
                                                                          StressGP_r_Theta_3);
            } else {
              StressGP_r_Theta_3.clear(); REL_Stress_Transformation(theta_angle,StressGP_Global,StressGP_r_Theta_3);
            }
          }
          else if (probdata.NameVars[i].compare(0,6,"theta4") == 0) {
            // w.r.t. ply angle
            // 4th-layer
            FE2_PostProcStressForReliability_First Calc_Stress_Theta_4(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Theta_4th_Ply",Dmat,Dmat_r_Theta_4);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Theta_4);  CHKERRQ(ierr);
            StressGP_Global.clear(); StressGP_Global = Calc_Stress_Theta_4.StressGP_r;
            if (probdata.ExaminedPly == 4) {
              StressGP_r_Theta_4.clear(); REL_Stress_Transformation_Theta(theta_angle,
                                                                          Calc_Stress.StressGP,
                                                                          Calc_Stress_Theta_4.StressGP_r,
                                                                          StressGP_r_Theta_4);
            } else {
              StressGP_r_Theta_4.clear(); REL_Stress_Transformation(theta_angle,StressGP_Global,StressGP_r_Theta_4);
            }
          }
        }
        
        // Evaluate LSF and its gradient
        grad_g.resize(probdata.num_vars); grad_g.clear();
        
        switch (FailureCriterion) {
          case 12013: {
             //
             // Maximum stress theory - fibre failure
             //
             MsFE_Reliability_Results.NameOfFailureCriterion = "Maximum stress theory - Fibre failure";
             ierr = TheLSF.gfun_ply_MS_LD(x,probdata.NameVars,probdata.MatStrength,
                                          StressGP,
                                          StressGP_r_Em,StressGP_r_NUm,
                                          StressGP_r_NUp,StressGP_r_NUpz,
                                          StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                          StressGP_r_Ef,StressGP_r_NUf,
                                          StressGP_r_F,StressGP_r_Theta,
                                          StressGP_r_Theta_1,StressGP_r_Theta_2,
                                          StressGP_r_Theta_3,StressGP_r_Theta_4,
                                          val_G,grad_g); CHKERRQ(ierr);
              break;
          }
          case 22013: {
             //
             // Maximum stress theory - matrix failure
             //
             MsFE_Reliability_Results.NameOfFailureCriterion = "Maximum stress theory - Matrix failure";
             ierr = TheLSF.gfun_ply_MS_TD(x,probdata.NameVars,probdata.MatStrength,
                                          StressGP,
                                          StressGP_r_Em,StressGP_r_NUm,
                                          StressGP_r_NUp,StressGP_r_NUpz,
                                          StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                          StressGP_r_Ef,StressGP_r_NUf,
                                          StressGP_r_F,StressGP_r_Theta,
                                          StressGP_r_Theta_1,StressGP_r_Theta_2,
                                          StressGP_r_Theta_3,StressGP_r_Theta_4,
                                          val_G,grad_g); CHKERRQ(ierr);
             break;
          }
          case 22014: {
             //
             // Maximum stress theory - shear failure
             //
             MsFE_Reliability_Results.NameOfFailureCriterion = "Maximum stress theory - Shear failure";
             ierr = TheLSF.gfun_ply_MS_Shear(x,probdata.NameVars,probdata.MatStrength,
                                             StressGP,
                                             StressGP_r_Em,StressGP_r_NUm,
                                             StressGP_r_NUp,StressGP_r_NUpz,
                                             StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                             StressGP_r_Ef,StressGP_r_NUf,
                                             StressGP_r_F,StressGP_r_Theta,
                                             StressGP_r_Theta_1,StressGP_r_Theta_2,
                                             StressGP_r_Theta_3,StressGP_r_Theta_4,
                                             val_G,grad_g); CHKERRQ(ierr);
             break;
          }
          case 13033: {
             //
             // Hashin failure theory - fibre failure
             //
             MsFE_Reliability_Results.NameOfFailureCriterion = "Hashin failure theory - Fibre failure";
             ierr = TheLSF.gfun_ply_HF(x,probdata.NameVars,probdata.MatStrength,
                                       StressGP,
                                       StressGP_r_Em,StressGP_r_NUm,
                                       StressGP_r_NUp,StressGP_r_NUpz,
                                       StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                       StressGP_r_Ef,StressGP_r_NUf,
                                       StressGP_r_F,StressGP_r_Theta,
                                       StressGP_r_Theta_1,StressGP_r_Theta_2,
                                       StressGP_r_Theta_3,StressGP_r_Theta_4,
                                       val_G,grad_g); CHKERRQ(ierr);
             break;
          }
          case 23033: {
             //
             // Hashin failure theory - matrix failure
             //
             MsFE_Reliability_Results.NameOfFailureCriterion = "Hashin failure theory - Matrix failure";
             ierr = TheLSF.gfun_ply_HM(x,probdata.NameVars,probdata.MatStrength,
                                       StressGP,
                                       StressGP_r_Em,StressGP_r_NUm,
                                       StressGP_r_NUp,StressGP_r_NUpz,
                                       StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                       StressGP_r_Ef,StressGP_r_NUf,
                                       StressGP_r_F,StressGP_r_Theta,
                                       StressGP_r_Theta_1,StressGP_r_Theta_2,
                                       StressGP_r_Theta_3,StressGP_r_Theta_4,
                                       val_G,grad_g); CHKERRQ(ierr);
             break;
          }
          case 42050: {
             //
             // Tsai-Wu failure criteria
             //
             MsFE_Reliability_Results.NameOfFailureCriterion = "Tsai-Wu - 2D stress state";
             ierr = TheLSF.gfun_ply_Tsai_Wu_2D(x,probdata.NameVars,probdata.MatStrength,
                                               StressGP,
                                               StressGP_r_Em,StressGP_r_NUm,
                                               StressGP_r_NUp,StressGP_r_NUpz,
                                               StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                               StressGP_r_Ef,StressGP_r_NUf,
                                               StressGP_r_F,StressGP_r_Theta,
                                               StressGP_r_Theta_1,StressGP_r_Theta_2,
                                               StressGP_r_Theta_3,StressGP_r_Theta_4,
                                               val_G,grad_g); CHKERRQ(ierr);
             break;
          }
          case 43050: {
            //
            // Tsai-Wu failure criteria
            //
            MsFE_Reliability_Results.NameOfFailureCriterion = "Tsai-Wu - 3D stress state";
            ierr = TheLSF.gfun_ply_Tsai_Wu(x,probdata.NameVars,probdata.MatStrength,
                                           StressGP,
                                           StressGP_r_Em,StressGP_r_NUm,
                                           StressGP_r_NUp,StressGP_r_NUpz,
                                           StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                           StressGP_r_Ef,StressGP_r_NUf,
                                           StressGP_r_F,StressGP_r_Theta,
                                           StressGP_r_Theta_1,StressGP_r_Theta_2,
                                           StressGP_r_Theta_3,StressGP_r_Theta_4,
                                           val_G,grad_g); CHKERRQ(ierr);
            break;
          }
          case 44050: {
             //
             // Tsai-Wu failure criteria
             //
             MsFE_Reliability_Results.NameOfFailureCriterion = "Tsai-Wu-Christensen - 3D stress state";
             ierr = TheLSF.gfun_ply_Tsai_Wu_Christensen(x,probdata.NameVars,probdata.MatStrength,
                                                        StressGP,
                                                        StressGP_r_Em,StressGP_r_NUm,
                                                        StressGP_r_NUp,StressGP_r_NUpz,
                                                        StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                                        StressGP_r_Ef,StressGP_r_NUf,
                                                        StressGP_r_F,StressGP_r_Theta,
                                                        StressGP_r_Theta_1,StressGP_r_Theta_2,
                                                        StressGP_r_Theta_3,StressGP_r_Theta_4,
                                                        val_G,grad_g); CHKERRQ(ierr);
            break;
          }
          case 42060: {
             //
             // Tsai-Hill failure criteria
             //
             MsFE_Reliability_Results.NameOfFailureCriterion = "Tsai-Hill - 2D stress-state";
             ierr = TheLSF.gfun_ply_Tsai_Hill_2D(x,probdata.NameVars,probdata.MatStrength,
                                                 StressGP,
                                                 StressGP_r_Em,StressGP_r_NUm,
                                                 StressGP_r_NUp,StressGP_r_NUpz,
                                                 StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                                 StressGP_r_Ef,StressGP_r_NUf,
                                                 StressGP_r_F,StressGP_r_Theta,
                                                 StressGP_r_Theta_1,StressGP_r_Theta_2,
                                                 StressGP_r_Theta_3,StressGP_r_Theta_4,
                                                 val_G,grad_g); CHKERRQ(ierr);
            break;
          }
          case 43060: {
             //
             // Tsai-Hill failure criteria
             //
             MsFE_Reliability_Results.NameOfFailureCriterion = "Tsai-Hill - 3D stress-state";
             ierr = TheLSF.gfun_ply_Tsai_Hill(x,probdata.NameVars,probdata.MatStrength,
                                              StressGP,
                                              StressGP_r_Em,StressGP_r_NUm,
                                              StressGP_r_NUp,StressGP_r_NUpz,
                                              StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                              StressGP_r_Ef,StressGP_r_NUf,
                                              StressGP_r_F,StressGP_r_Theta,
                                              StressGP_r_Theta_1,StressGP_r_Theta_2,
                                              StressGP_r_Theta_3,StressGP_r_Theta_4,
                                              val_G,grad_g); CHKERRQ(ierr);
            break;
          }
          case 13073: {
             //
             // Richard Christen: Fibre controlled failure
             //
             MsFE_Reliability_Results.NameOfFailureCriterion = "Christensen - Fibre controlled failure";
             ierr = TheLSF.gfun_ply_RCF(x,probdata.NameVars,probdata.MatStrength,
                                        StressGP,
                                        StressGP_r_Em,StressGP_r_NUm,
                                        StressGP_r_NUp,StressGP_r_NUpz,
                                        StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                        StressGP_r_Ef,StressGP_r_NUf,
                                        StressGP_r_F,StressGP_r_Theta,
                                        StressGP_r_Theta_1,StressGP_r_Theta_2,
                                        StressGP_r_Theta_3,StressGP_r_Theta_4,
                                        val_G,grad_g); CHKERRQ(ierr);
            break;
          }
          case 23073: {
             //
             // Richard Christen: Matrix controlled failure
             //
             MsFE_Reliability_Results.NameOfFailureCriterion = "Christensen - Matrix controlled failure";
             ierr = TheLSF.gfun_ply_RCM(x,probdata.NameVars,probdata.MatStrength,
                                        StressGP,
                                        StressGP_r_Em,StressGP_r_NUm,
                                        StressGP_r_NUp,StressGP_r_NUpz,
                                        StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                        StressGP_r_Ef,StressGP_r_NUf,
                                        StressGP_r_F,StressGP_r_Theta,
                                        StressGP_r_Theta_1,StressGP_r_Theta_2,
                                        StressGP_r_Theta_3,StressGP_r_Theta_4,
                                        val_G,grad_g); CHKERRQ(ierr);
            break;
          }
          case 42080: {
             //
             // Hoffman failure theory
             //
             MsFE_Reliability_Results.NameOfFailureCriterion = "Hoffman failure theory - 2D stress states";
             ierr = TheLSF.gfun_ply_Hoffman_2D(x,probdata.NameVars,probdata.MatStrength,
                                               StressGP,
                                               StressGP_r_Em,StressGP_r_NUm,
                                               StressGP_r_NUp,StressGP_r_NUpz,
                                               StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                               StressGP_r_Ef,StressGP_r_NUf,
                                               StressGP_r_F,StressGP_r_Theta,
                                               StressGP_r_Theta_1,StressGP_r_Theta_2,
                                               StressGP_r_Theta_3,StressGP_r_Theta_4,
                                               val_G,grad_g); CHKERRQ(ierr);
            break;
          }
          case 43080: {
             //
             // Hoffman failure theory
             //
             MsFE_Reliability_Results.NameOfFailureCriterion = "Hoffman failure theory - 3D stress states";
             ierr = TheLSF.gfun_ply_Hoffman(x,probdata.NameVars,probdata.MatStrength,
                                            StressGP,
                                            StressGP_r_Em,StressGP_r_NUm,
                                            StressGP_r_NUp,StressGP_r_NUpz,
                                            StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                                            StressGP_r_Ef,StressGP_r_NUf,
                                            StressGP_r_F,StressGP_r_Theta,
                                            StressGP_r_Theta_1,StressGP_r_Theta_2,
                                            StressGP_r_Theta_3,StressGP_r_Theta_4,
                                            val_G,grad_g); CHKERRQ(ierr);
            break;
          }
            //default: {}
        }
        
        grad_G.resize(probdata.num_vars); grad_G.clear();
        grad_G = prod(grad_g,inv_dudx);
        
        cout<<"LSF value is \t"<<val_G<<endl;
        cout<<"Gradient of LSF is \t"<<grad_g<<"\t"<<grad_G<<endl;
        
        //
        // Set scale parameter G0 and inform about structural response
        //
        if (istep == 1) {
          val_G0 = val_G;
          if (echo_flag) {
            cout<<"Value of limit-state function in the first step: \t"<<val_G<<endl;
          }
        }
        
        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //        STEP 6: compute the direction cosine                            //
        //                a) dg/dy = inv(J) dG/dx                                 //
        //                b) alpha_i = dg/dy_i / norm(dg/dy)                      //
        //                c) reliability index estimate: beta = sqrt(norm(y))     //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        
        //
        // Compute direction cosine alpha vector
        //
        alpha.resize(probdata.num_vars); alpha.clear();
        alpha = -grad_G/norm_2(grad_G);
        cout<<"Direction cosine "<<alpha<<endl;
        
        //
        // Check convergence
        //
        if (((abs(val_G/val_G0)<e1) && (norm_2(u-inner_prod(alpha,u)*alpha)<e2)) || (istep == istep_max)) {
          conv_flag = 1;
        }
        
        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //        STEP 7: convergence check and compute new trial point:          //
        //                7.1 convergence check                                   //
        //                    (a) design point                                    //
        //                    (b) estimate of reliability index                   //
        //                7.2 new trial point                                     //
        //                           y_n = -alpha[beta + g/norm(dg/dy)]           //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        
        //
        // Take a step if convergence is not achieved
        //
        if (conv_flag == 0) {
          // Determine search direction
          u_dir.resize(probdata.num_vars); u_dir.clear();
          search_dir(val_G,grad_G,u,u_dir);
          // Determine step size
          if (step_code != 0) {
            step_size = step_code;
          }
          else {
            cout<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!";
            cout<<"\nArmijo rule will be used to determine step size for setting new trial point, \n";
            cout<<"but it isn't available yet in current version!!! \n";
            cout<<"Re-assign a nonzero value to step_code.\n";
            cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
            conv_flag = 1;
          }
          
          cout<<"\n\nReliability index estimate at "<<istep<<" is: "<<((val_G/norm_2(grad_G)) + inner_prod(alpha, u));
          cout<<"\t"<<norm_2(u)<<"\t"<<inner_prod(alpha,u)<<endl;
          //cout<<"design point: "<<u<<endl;
          
          // Write reliability index in file
          if (beta_flag == 1) {
            BetaFile<<setprecision(15)<<inner_prod(alpha,u)<<"\t"<<val_G;
            for (int i=0;i<probdata.num_vars;i++) {
              BetaFile<<"\t"<<x(i);
            }
            BetaFile<<"\n";
          }
          
          
          // Determin new trial point
          u_new.resize(probdata.num_vars); u_new.clear();
          u_new = u + step_size*u_dir;  // when step_size is 1, it is HLRF search algorithm
          // Prepare for a new round in the loop
          u.resize(probdata.num_vars); u.clear();
          u = u_new;
          istep = istep + 1;
        }
        
        rel_t2 = clock();
        rel_calc_time  = (double)(rel_t2 - rel_t1)/CLOCKS_PER_SEC;
        rel_t1 = rel_t2;
        cout<<"Elapsed time at this step is: "<<rel_calc_time<<" seconds.\n";
        
      } while (conv_flag == 0);
      
      // Close beta value writting file
      if (beta_flag == 1) { BetaFile.close(); }
      
      //////////////////////////////////////////////////////////////////////////////
      //                                                                          //
      //        STEP 8: Calculate the estimate of reliability index               //
      //                                                                          //
      //////////////////////////////////////////////////////////////////////////////
      
      MsFE_Reliability_Results.beta_FORM = inner_prod(alpha,u);
      MsFE_Reliability_Results.DesignPoint = u;
      
      cout<<u<<endl;
      
      using boost::math::normal_distribution;
      normal_distribution<> snorm(0,1);
      MsFE_Reliability_Results.prob_failure_FORM = cdf(snorm,-MsFE_Reliability_Results.beta_FORM);
      
      //////////////////////////////////////////////////////////////////////////////
      //                                                                          //
      //                               FINISH                                     //
      //                                                                          //
      //////////////////////////////////////////////////////////////////////////////
      
      
      cout<<"\n\n*************************************************\n*\n";
      if (istep == istep_max) {
        cout<<"*  The maximum number of iteration is reached.\n";
      } else {
        cout<<"*  The number of iterations is: "<<istep<<".\n";
      }
      
      cout<<"*  Reliability method is: FORM"<<endl;
      cout<<"*  Optimal reliability index (beta) is: "<<MsFE_Reliability_Results.beta_FORM<<endl;
      cout<<"*  Optimal failure probability is:      "<<setprecision(15)<<MsFE_Reliability_Results.prob_failure_FORM<<endl;
      cout<<"*  The failure criterion used is:       "<<MsFE_Reliability_Results.NameOfFailureCriterion<<endl;
      
      PetscFunctionReturn(0);
    }
    
    
    /***************************************************************************
     *                                                                         *
     *                                                                         *
     *             MONTE CARLO IMPORTANCE SAMPLING - MCIS                      *
     *                                                                         *
     *                                                                         *
     **************************************************************************/
    
    virtual PetscErrorCode Crude_MCS(FieldInterface &m_field_RVE,
                                     FieldInterface &m_field_Macro,
                                     PetscInt NO_Layers,
                                     Stochastic_Model probdata,
                                     Reliability_Options ReliabOpt,
                                     int FailureCriterion,
                                     Reliability_Results &MsFE_Reliability_Results) {
      PetscFunctionBegin;
      
      double val_G_TW;
      /*double val_G_MS_LD, val_G_MS_TD, val_G_MS_Shear;
      double val_G_HF, val_G_HM;
      double val_G_TW_2D, val_G_TW;
      double val_G_TH_2D, val_G_TH;
      double val_G_RCF, val_G_RCM;
      double val_G_Hoffman_2D, val_G_Hoffman;*/
      
      int no_fail = 0;
      int no_mcs = 120000;
      ublas::vector<double> x;
      ublas::matrix<double> StressGP;
      ublas::matrix<double> StressGP_1st, StressGP_2nd, StressGP_3rd, StressGP_4th;
      
      double theta_angle;
      ublas::vector<double> PlyAngle_new;
      PlyAngle_new = probdata.PlyAngle;
      
      //
      // Start loop
      //
      
      clock_t rel_t1, rel_t2;
      double rel_calc_time;
      rel_t1 = clock();
      
      ublas::matrix<double> Dmat;
      ublas::matrix<double> Dmat_1st, Dmat_2nd, Dmat_3rd, Dmat_4th;
      
      MCS_RNG sampling_values;
      FE2_Macro_Solver_MCS Solve_FE2_Problem;
      LimitStateFunction_MCS TheLSF;
      
      for (int imcs=1; imcs<=no_mcs; imcs++) {
        cout<<"\n\n*************************************************\n*\n";
        cout<<"*    This is "<<imcs<<" step!"<<"\n*\n";
        cout<<"*************************************************\n";
        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //        STEP 2: generate samples                                        //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        ierr = sampling_values.myRNG(probdata.num_vars,probdata.marg,x); CHKERRQ(ierr);
        cout<<"\nThe random variables: "<<x<<endl;
        
        // ----------------
        // Update ply angle
        // -----------------
        for (unsigned i=1; i<=x.size();i++) {
          if (probdata.NameVars[i].compare(0,11,"orientation") == 0) {
            cout<<"\nAngle "<<x(i-1)<<endl;
            for (unsigned j=0; j<probdata.PlyAngle.size(); j++) {
              cout << "The original ply angle "<<probdata.PlyAngle(j)<<"\t delta x: "<<x(i-1)<<endl;
              PlyAngle_new(j) = probdata.PlyAngle(j) + x(i-1);
              cout << "The modified ply angle "<<PlyAngle_new(j)<<endl;
            }
          }
          else if (probdata.NameVars[i].compare(0,6,"theta1") == 0) {
            PlyAngle_new(0) = x(i-1);
          }
          else if (probdata.NameVars[i].compare(0,6,"theta2") == 0) {
            PlyAngle_new(1) = x(i-1);
          }
          else if (probdata.NameVars[i].compare(0,6,"theta3") == 0) {
            PlyAngle_new(2) = x(i-1);
          }
          else if (probdata.NameVars[i].compare(0,6,"theta4") == 0) {
            PlyAngle_new(3) = x(i-1);
          }
        }
        
        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //        STEP 3: call FE program to calculate structural response        //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        
        //
        // Evaluate limit-state function and its gradient
        //
        ierr = Solve_FE2_Problem.Calculate_RVEDmat(m_field_RVE,x,probdata.num_vars,
                                                   probdata.NameVars); CHKERRQ(ierr);
        ierr = Solve_FE2_Problem.Macro_FE_Laminate(m_field_Macro,x,probdata.num_vars,
                                                   probdata.NameVars,
                                                   PlyAngle_new,
                                                   NO_Layers); CHKERRQ(ierr);
        //
        // Get constitutive matrix or D matrix at Gauss points
        //
        switch (probdata.ExaminedPly) {
          case 0:{cout<<"\n\nAll plies are under examination.\n\n";
            //--------------------------
            // D matrix for all plies
            //--------------------------
            Dmat_1st = Solve_FE2_Problem.Dmat_1st_Ply;
            Dmat_2nd = Solve_FE2_Problem.Dmat_2nd_Ply;
            Dmat_3rd = Solve_FE2_Problem.Dmat_3rd_Ply;
            Dmat_4th = Solve_FE2_Problem.Dmat_4th_Ply;
            break;
          }
          case 1:{cout<<"\n\nThe 1st layer is under examination.\n\n";
            //--------------------------
            // D matrix for the 1st ply
            //--------------------------
            Dmat = Solve_FE2_Problem.Dmat_1st_Ply;
            break;
          }
          case 2: {cout<<"\n\nThe 2nd layer is under examination.\n\n";
            //--------------------------
            // D matrix for the 2nd ply
            //--------------------------
            Dmat = Solve_FE2_Problem.Dmat_2nd_Ply;
            break;
          }
          case 3: {cout<<"\n\nThe 3rd layer is under examination.\n\n";
            //--------------------------
            // D matrix for the 3rd ply
            //--------------------------
            Dmat = Solve_FE2_Problem.Dmat_3rd_Ply;
            break;
          }
          case 4: {cout<<"\n\nThe 4th layer is under examination.\n\n";
            //--------------------------
            // D matrix for the 4th ply
            //--------------------------
            Dmat = Solve_FE2_Problem.Dmat_4th_Ply;
            break;
          }
        }
        
        //
        // Calculate the stress at Gauss points for specific element(s)
        //
        if (probdata.ExaminedPly == 0) {
          FE2_PostProcStressForReliability_Zeroth Calc_Stress_1st(m_field_Macro,"DISP_MACRO",Dmat_1st);
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_1st",Calc_Stress_1st);  CHKERRQ(ierr);
          
          FE2_PostProcStressForReliability_Zeroth Calc_Stress_2nd(m_field_Macro,"DISP_MACRO",Dmat_2nd);
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_2nd",Calc_Stress_2nd);  CHKERRQ(ierr);
          
          FE2_PostProcStressForReliability_Zeroth Calc_Stress_3rd(m_field_Macro,"DISP_MACRO",Dmat_3rd);
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_3rd",Calc_Stress_3rd);  CHKERRQ(ierr);
          
          FE2_PostProcStressForReliability_Zeroth Calc_Stress_4th(m_field_Macro,"DISP_MACRO",Dmat_4th);
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_4th",Calc_Stress_4th);  CHKERRQ(ierr);
          
          theta_angle = PlyAngle_new(0)*(M_PI/180.0);
          StressGP_1st.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_1st.StressGP,StressGP_1st);
          cout<<"1st ply stress at GP: "<<StressGP_1st<<endl;
          theta_angle = PlyAngle_new(1)*(M_PI/180.0);
          StressGP_2nd.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_2nd.StressGP,StressGP_2nd);
          cout<<"2nd ply stress at GP: "<<StressGP_2nd<<endl;
          theta_angle = PlyAngle_new(2)*(M_PI/180.0);
          StressGP_3rd.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_3rd.StressGP,StressGP_3rd);
          cout<<"3rd ply stress at GP: "<<StressGP_3rd<<endl;
          theta_angle = PlyAngle_new(3)*(M_PI/180.0);
          StressGP_4th.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_4th.StressGP,StressGP_4th);
          cout<<"4th ply stress at GP: "<<StressGP_4th<<endl;
        }
        else {
          FE2_PostProcStressForReliability_Zeroth Calc_Stress(m_field_Macro,"DISP_MACRO",Dmat);
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress);  CHKERRQ(ierr);
          
          theta_angle = PlyAngle_new(0)*(M_PI/180.0);
          StressGP.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress.StressGP,StressGP);
          cout<<"1st ply stress at GP: "<<StressGP<<endl;
        }
        
        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //        STEP 4: evaluate limit state function                           //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        
        switch (FailureCriterion) {
          case 43050: {
            //
            // Tsai-Wu failure criteria
            //
            MsFE_Reliability_Results.NameOfFailureCriterion = "Tsai-Wu";
            ierr = TheLSF.gfun_ply_Tsai_Wu(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_TW); CHKERRQ(ierr);
            
            if (val_G_TW<0) { no_fail++;}
            break;
          }
        }
        
        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //        STEP 4: output limit state function                             //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        
        rel_t2 = clock();
        rel_calc_time  = (double)(rel_t2 - rel_t1)/CLOCKS_PER_SEC;
        rel_t1 = rel_t2;
        cout<<"Elapsed time at this step is: "<<rel_calc_time<<" seconds.\n";
      }
      
      //////////////////////////////////////////////////////////////////////////////
      //                                                                          //
      //        STEP 5: Calculate the probability of failure and corresponding    //
      //                estimate of reliability index                             //
      //                                                                          //
      //////////////////////////////////////////////////////////////////////////////
      
      MsFE_Reliability_Results.prob_failure_MCS = no_fail/no_mcs;
      using boost::math::normal_distribution;
      normal_distribution<> snorm(0,1);
      if (no_fail == 0) {
        MsFE_Reliability_Results.beta_MCS = 0.0;
      } else {
        MsFE_Reliability_Results.beta_MCS = quantile(snorm,MsFE_Reliability_Results.prob_failure_MCS);
      }
      
      //////////////////////////////////////////////////////////////////////////////
      //                                                                          //
      //                               FINISH                                     //
      //                                                                          //
      //////////////////////////////////////////////////////////////////////////////
      
      cout<<"\n\n*************************************************\n*\n";
      //if (istep == istep_max) {
      //  cout<<"*  The maximum number of iteration is reached.\n";
      //} else {
      //  cout<<"*  The number of iterations is: "<<istep<<".\n";
      //}
      cout<<"*  Reliability method is Crude Monte Carlo simulation."<<endl;
      cout<<"*  Reliability index (beta) is:    "<<MsFE_Reliability_Results.beta_MCS<<endl;
      cout<<"*  Probability of failure (pf) is: "<<setprecision(15)<<MsFE_Reliability_Results.prob_failure_MCS<<endl;
      cout<<"*  Number of failure events is:    "<<no_fail<<endl;
      cout<<"*  The failure criterion used is:  "<<MsFE_Reliability_Results.NameOfFailureCriterion<<endl;
      
      PetscFunctionReturn(0);
    }
    
    
    /***************************************************************************
     *                                                                         *
     *                                                                         *
     *             MONTE CARLO IMPORTANCE SAMPLING - MCIS                      *
     *                                                                         *
     *                                                                         *
     **************************************************************************/

    virtual PetscErrorCode MCIS(FieldInterface &m_field_RVE,
                                FieldInterface &m_field_Macro,
                                PetscInt NO_Layers,
                                Stochastic_Model probdata,
                                Reliability_Options ReliabOpt,
                                int FailureCriterion,
                                Reliability_Results &MsFE_Reliability_Results) {
      PetscFunctionBegin;
      
      cout<<"\n\nHello from Importance sampling!\n\n"<<endl;
      
      double val_G_MS_LD, val_G_MS_TD, val_G_MS_Shear;
      double val_G_HF, val_G_HM;
      double val_G_TW;
      double val_G_TH;
      double val_G_RCF, val_G_RCM;
      double val_G_Hoffman;
      
      int no_mcs = 120000;
      ublas::vector<double> x;
      ublas::vector<double> u, u_temp;
      ublas::matrix<double> StressGP;
      ublas::matrix<double> StressGP_1st, StressGP_2nd, StressGP_3rd, StressGP_4th;
      
      double theta_angle;
      ublas::vector<double> PlyAngle_new;
      PlyAngle_new = probdata.PlyAngle;
      
      //
      // Start loop
      //
      
      clock_t rel_t1, rel_t2;
      double rel_calc_time;
      rel_t1 = clock();
      
      ublas::matrix<double> Dmat;
      ublas::matrix<double> Dmat_1st, Dmat_2nd, Dmat_3rd, Dmat_4th;
      

      FE2_Macro_Solver_MCS Solve_FE2_Problem;
      LimitStateFunction_MCS TheLSF;

      NatafTransformation my_nataf_transformation;
      ierr = my_nataf_transformation.ModCorrMat_Empirical(probdata.num_vars,
                                                          probdata.marg,
                                                          probdata.correlation,
                                                          probdata.mod_correlation); CHKERRQ(ierr);
      //
      // Perform Cholesky decomposition for the modified correlation matrix
      //  A = LL'
      //
      ublas::triangular_matrix<double, ublas::lower> Lmat(probdata.num_vars,probdata.num_vars);
      cholesky_decompose(probdata.mod_correlation,Lmat);
      probdata.Lo.resize(probdata.num_vars,probdata.num_vars);
      probdata.Lo = Lmat;
      
      //
      // Compute the inverse of Lo
      //
      probdata.inv_Lo.resize(probdata.num_vars,probdata.num_vars);
      bool singular = false;
      probdata.inv_Lo = gjinverse(probdata.Lo,singular);
      
      double detj = 1.0;
      double pf = 0.0;
      double vpf = 0.0;

      MsFE_Reliability_Results.prob_failure_MCIS = MsFE_Reliability_Results.prob_failure_FORM;

      boost::mt19937 rng(time(0));
      boost::random::normal_distribution<> my_norm(0,1);
      
      ofstream MCIS_File;
      MCIS_File.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//04_ReliabilityAnalysis//Result_MCIS.txt",ofstream::out);
      
      for (int imcs=1; imcs<=no_mcs; imcs++) {
        
        cout<<"\n\n*************************************************\n*\n";
        cout<<"*    This is "<<imcs<<" step!"<<"\n*\n";
        cout<<"*************************************************\n";
        
        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //        STEP 2: generate samples                                        //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        
        //
        // Generate random points in standard normal space
        //
        u.resize(probdata.num_vars); u.clear();
        for (int irv = 0; irv<probdata.num_vars; irv++) {
          boost::variate_generator <boost::mt19937&,boost::random::normal_distribution<> > normrnd(rng,my_norm);
          u(irv) = normrnd();cout<<irv+1<<"-th rv is "<<u(irv)<<endl;
        }
        cout<<u<<endl;
        //
        // Move points to
        //
        u_temp.resize(probdata.num_vars); u_temp.clear();
        u_temp = u + MsFE_Reliability_Results.DesignPoint;
        //
        // Transform random number from standard normal space to original space
        //
        x.resize(probdata.num_vars); x.clear();
        ierr = my_nataf_transformation.u_to_x(u_temp,probdata.num_vars,probdata.marg,probdata.Lo,x,detj); CHKERRQ(ierr);
        cout<<"\nRandom vaiables in U-space: "<<u<<endl;
        cout<<"\nRandom vaiables in X-space: "<<x<<endl;
        
        // ----------------
        // Update ply angle
        // -----------------
        for (unsigned i=1; i<=x.size();i++) {
          if (probdata.NameVars[i].compare(0,11,"orientation") == 0) {
            cout<<"\nAngle "<<x(i-1)<<endl;
            for (unsigned j=0; j<probdata.PlyAngle.size(); j++) {
              cout << "The original ply angle "<<probdata.PlyAngle(j)<<"\t delta x: "<<x(i-1)<<endl;
              PlyAngle_new(j) = probdata.PlyAngle(j) + x(i-1);
              cout << "The modified ply angle "<<PlyAngle_new(j)<<endl;
            }
          }
          else if (probdata.NameVars[i].compare(0,6,"theta1") == 0) {
            PlyAngle_new(0) = x(i-1);
          }
          else if (probdata.NameVars[i].compare(0,6,"theta2") == 0) {
            PlyAngle_new(1) = x(i-1);
          }
          else if (probdata.NameVars[i].compare(0,6,"theta3") == 0) {
            PlyAngle_new(2) = x(i-1);
          }
          else if (probdata.NameVars[i].compare(0,6,"theta4") == 0) {
            PlyAngle_new(3) = x(i-1);
          }
        }
        
        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //        STEP 3: call FE program to calculate structural response        //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        
        //
        // Evaluate limit-state function and its gradient
        //
        ierr = Solve_FE2_Problem.Calculate_RVEDmat(m_field_RVE,x,probdata.num_vars,
                                                   probdata.NameVars); CHKERRQ(ierr);
        ierr = Solve_FE2_Problem.Macro_FE_Laminate(m_field_Macro,x,probdata.num_vars,
                                                   probdata.NameVars,
                                                   PlyAngle_new,
                                                   NO_Layers); CHKERRQ(ierr);
        //
        // Get constitutive matrix or D matrix at Gauss points
        //
        switch (probdata.ExaminedPly) {
          case 0:{cout<<"\n\nAll plies are under examination.\n\n";
            //--------------------------
            // D matrix for all plies
            //--------------------------
            Dmat_1st = Solve_FE2_Problem.Dmat_1st_Ply;
            Dmat_2nd = Solve_FE2_Problem.Dmat_2nd_Ply;
            Dmat_3rd = Solve_FE2_Problem.Dmat_3rd_Ply;
            Dmat_4th = Solve_FE2_Problem.Dmat_4th_Ply;
            break;
          }
          case 1:{cout<<"\n\nThe 1st layer is under examination.\n\n";
            //--------------------------
            // D matrix for the 1st ply
            //--------------------------
            Dmat = Solve_FE2_Problem.Dmat_1st_Ply;
            break;
          }
          case 2: {cout<<"\n\nThe 2nd layer is under examination.\n\n";
            //--------------------------
            // D matrix for the 2nd ply
            //--------------------------
            Dmat = Solve_FE2_Problem.Dmat_2nd_Ply;
            break;
          }
          case 3: {cout<<"\n\nThe 3rd layer is under examination.\n\n";
            //--------------------------
            // D matrix for the 3rd ply
            //--------------------------
            Dmat = Solve_FE2_Problem.Dmat_3rd_Ply;
            break;
          }
          case 4: {cout<<"\n\nThe 4th layer is under examination.\n\n";
            //--------------------------
            // D matrix for the 4th ply
            //--------------------------
            Dmat = Solve_FE2_Problem.Dmat_4th_Ply;
            break;
          }
        }
        
        //
        // Calculate the stress at Gauss points for specific element(s)
        //
        if (probdata.ExaminedPly == 0) {
          FE2_PostProcStressForReliability_Zeroth Calc_Stress_1st(m_field_Macro,"DISP_MACRO",Dmat_1st);
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_1st",Calc_Stress_1st);  CHKERRQ(ierr);
          
          FE2_PostProcStressForReliability_Zeroth Calc_Stress_2nd(m_field_Macro,"DISP_MACRO",Dmat_2nd);
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_2nd",Calc_Stress_2nd);  CHKERRQ(ierr);
          
          FE2_PostProcStressForReliability_Zeroth Calc_Stress_3rd(m_field_Macro,"DISP_MACRO",Dmat_3rd);
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_3rd",Calc_Stress_3rd);  CHKERRQ(ierr);
          
          FE2_PostProcStressForReliability_Zeroth Calc_Stress_4th(m_field_Macro,"DISP_MACRO",Dmat_4th);
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_4th",Calc_Stress_4th);  CHKERRQ(ierr);
          
          theta_angle = PlyAngle_new(0)*(M_PI/180.0);
          StressGP_1st.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_1st.StressGP,StressGP_1st);
          cout<<"1st ply stress at GP: "<<StressGP_1st<<endl;
          theta_angle = PlyAngle_new(1)*(M_PI/180.0);
          StressGP_2nd.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_2nd.StressGP,StressGP_2nd);
          cout<<"2nd ply stress at GP: "<<StressGP_2nd<<endl;
          theta_angle = PlyAngle_new(2)*(M_PI/180.0);
          StressGP_3rd.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_3rd.StressGP,StressGP_3rd);
          cout<<"3rd ply stress at GP: "<<StressGP_3rd<<endl;
          theta_angle = PlyAngle_new(3)*(M_PI/180.0);
          StressGP_4th.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_4th.StressGP,StressGP_4th);
          cout<<"4th ply stress at GP: "<<StressGP_4th<<endl;
        }
        else {
          FE2_PostProcStressForReliability_Zeroth Calc_Stress(m_field_Macro,"DISP_MACRO",Dmat);
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress);  CHKERRQ(ierr);
          
          theta_angle = PlyAngle_new(probdata.ExaminedPly-1)*(M_PI/180.0);
          StressGP.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress.StressGP,StressGP);
          cout<<"1st ply stress at GP: "<<StressGP<<endl;
        }
        
        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //        STEP 4: evaluate limit state function                           //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        
        switch (FailureCriterion) {
          case 12013: {
            //
            // Maximum stress theory - fibre failure
            //
            MsFE_Reliability_Results.NameOfFailureCriterion = "Maximum stress theory - Fibre failure";
            ierr = TheLSF.gfun_ply_MS_LD(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_MS_LD); CHKERRQ(ierr);
            cout<<"LSF value is: "<<val_G_MS_LD<<endl;
            Update_Pf_MCIS(u, val_G_MS_LD, pf, vpf, imcs, MsFE_Reliability_Results);
            break;
          }
          case 22013: {
            //
            // Maximum stress theory - matrix failure
            //
            MsFE_Reliability_Results.NameOfFailureCriterion = "Maximum stress theory - Matrix failure";
            ierr = TheLSF.gfun_ply_MS_TD(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_MS_TD); CHKERRQ(ierr);
            cout<<"LSF value is: "<<val_G_MS_TD<<endl;
            Update_Pf_MCIS(u, val_G_MS_TD, pf, vpf, imcs, MsFE_Reliability_Results);
            break;
          }
          case 22014: {
            //
            // Maximum stress theory - shear failure
            //
            MsFE_Reliability_Results.NameOfFailureCriterion = "Maximum stress theory - Shear failure";
            ierr = TheLSF.gfun_ply_MS_Shear(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_MS_Shear); CHKERRQ(ierr);
            cout<<"LSF value is: "<<val_G_MS_Shear<<endl;
            Update_Pf_MCIS(u, val_G_MS_Shear, pf, vpf, imcs, MsFE_Reliability_Results);
            break;
          }
          case 13033: {
            //
            // Hashin failure theory - fibre failure
            //
            MsFE_Reliability_Results.NameOfFailureCriterion = "Hashin failure theory - Fibre failure";
            ierr = TheLSF.gfun_ply_HF(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_HF); CHKERRQ(ierr);
            cout<<"LSF value is: "<<val_G_HF<<endl;
            Update_Pf_MCIS(u, val_G_HF, pf, vpf, imcs, MsFE_Reliability_Results);
            break;
          }
          case 23033: {
            //
            // Hashin failure theory - matrix failure
            //
            MsFE_Reliability_Results.NameOfFailureCriterion = "Hashin failure theory - Matrix failure";
            ierr = TheLSF.gfun_ply_HM(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_HM); CHKERRQ(ierr);
            cout<<"LSF value is: "<<val_G_HM<<endl;
            Update_Pf_MCIS(u, val_G_HM, pf, vpf, imcs, MsFE_Reliability_Results);
            break;
          }
          case 43050: {
            //
            // Tsai-Wu failure criteria
            //
            MsFE_Reliability_Results.NameOfFailureCriterion = "Tsai-Wu";
            ierr = TheLSF.gfun_ply_Tsai_Wu(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_TW); CHKERRQ(ierr);
            cout<<"LSF value is: "<<val_G_TW<<endl;
            Update_Pf_MCIS(u, val_G_TW, pf, vpf, imcs, MsFE_Reliability_Results);
            break;
          }
          case 43060: {
            //
            // Tsai-Hill failure criteria
            //
            MsFE_Reliability_Results.NameOfFailureCriterion = "Tsai-Hill - 3D stress-state";
            ierr = TheLSF.gfun_ply_Tsai_Hill(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_TH); CHKERRQ(ierr);
            cout<<"LSF value is: "<<val_G_TH<<endl;
            Update_Pf_MCIS(u, val_G_TH, pf, vpf, imcs, MsFE_Reliability_Results);
            break;
          }
          case 13073: {
            //
            // Richard Christen: Fibre controlled failure
            //
            MsFE_Reliability_Results.NameOfFailureCriterion = "Christensen - Fibre controlled failure";
            ierr = TheLSF.gfun_ply_RCF(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_RCF); CHKERRQ(ierr);
            cout<<"LSF value is: "<<val_G_RCF<<endl;
            Update_Pf_MCIS(u, val_G_RCF, pf, vpf, imcs, MsFE_Reliability_Results);
            break;
          }
          case 23073: {
            //
            // Richard Christen: Matrix controlled failure
            //
            MsFE_Reliability_Results.NameOfFailureCriterion = "Christensen - Matrix controlled failure";
            ierr = TheLSF.gfun_ply_RCM(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_RCM); CHKERRQ(ierr);
            cout<<"LSF value is: "<<val_G_RCM<<endl;
            Update_Pf_MCIS(u, val_G_RCM, pf, vpf, imcs, MsFE_Reliability_Results);
            break;
          }
          case 43080: {
            //
            // Hoffman failure theory
            //
            MsFE_Reliability_Results.NameOfFailureCriterion = "Hoffman failure theory - 3D stress states";
            ierr = TheLSF.gfun_ply_Hoffman(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_Hoffman); CHKERRQ(ierr);
            cout<<"LSF value is: "<<val_G_Hoffman<<endl;
            Update_Pf_MCIS(u, val_G_Hoffman, pf, vpf, imcs, MsFE_Reliability_Results);
            break;
          }
        }
        
        cout<<"\n\n";
        cout<<"Failure probability from MCIS: "<<MsFE_Reliability_Results.prob_failure_MCIS<<endl;
        cout<<"Failure probability from FORM: "<<MsFE_Reliability_Results.prob_failure_FORM<<endl;
        
        MCIS_File<<setprecision(15)<<MsFE_Reliability_Results.prob_failure_MCIS<<"\t";
        MCIS_File<<setprecision(15)<<MsFE_Reliability_Results.beta_MCIS<<"\t";
        MCIS_File<<setprecision(15)<<MsFE_Reliability_Results.prob_failure_FORM<<"\t";
        MCIS_File<<setprecision(15)<<MsFE_Reliability_Results.beta_FORM<<"\n";
        
        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //        STEP 4: output limit state function                             //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        
        rel_t2 = clock();
        rel_calc_time  = (double)(rel_t2 - rel_t1)/CLOCKS_PER_SEC;
        rel_t1 = rel_t2;
        cout<<"Elapsed time at this step is: "<<rel_calc_time<<" seconds.\n";
      }
      
      //////////////////////////////////////////////////////////////////////////////
      //                                                                          //
      //        STEP 5: Calculate the probability of failure and corresponding    //
      //                estimate of reliability index                             //
      //                                                                          //
      //////////////////////////////////////////////////////////////////////////////
      
      
      
      //////////////////////////////////////////////////////////////////////////////
      //                                                                          //
      //                               FINISH                                     //
      //                                                                          //
      //////////////////////////////////////////////////////////////////////////////
      
      MCIS_File.close();
                     
      cout<<"\n\n*************************************************\n*\n";
      cout<<"*  Reliability method is Importance Sampling."<<endl;
      cout<<"*  Reliability index (beta) is:    "<<MsFE_Reliability_Results.beta_MCIS<<endl;
      cout<<"*  Probability of failure (pf) is: "<<setprecision(15)<<MsFE_Reliability_Results.prob_failure_MCIS<<endl;
      cout<<"*  The failure criterion used is:  "<<MsFE_Reliability_Results.NameOfFailureCriterion<<endl;
      
      PetscFunctionReturn(0);
    }
    
    
    /***************************************************************************
     *                                                                         *
     *                                                                         *
     *             MONTE CARLO IMPORTANCE SAMPLING - MCIS                      *
     *                                                                         *
     *                                                                         *
     **************************************************************************/
    
    virtual PetscErrorCode Multi_MCIS(FieldInterface &m_field_RVE,
                                      FieldInterface &m_field_Macro,
                                      PetscInt NO_Layers,
                                      Stochastic_Model probdata,
                                      Reliability_Options ReliabOpt,
                                      ublas::vector<int> FC,
                                      vector<Reliability_Results> Obj_Rel_Res) {
      PetscFunctionBegin;
      
      cout<<"\n\nHello from Importance sampling!\n\n"<<endl;
      
      double val_G_MS_LD, val_G_MS_TD, val_G_MS_Shear;
      double val_G_HF, val_G_HM;
      double val_G_TW;
      double val_G_TH;
      double val_G_RCF, val_G_RCM;
      double val_G_Hoffman;
      
      int no_mcs = 120000;
      ublas::vector<double> x;
      ublas::vector<double> u, u_temp;
      ublas::matrix<double> StressGP;
      ublas::matrix<double> StressGP_1st, StressGP_2nd, StressGP_3rd, StressGP_4th;
      
      double theta_angle;
      ublas::vector<double> PlyAngle_new;
      PlyAngle_new = probdata.PlyAngle;
      
      //
      // Start loop
      //
      
      clock_t rel_t1, rel_t2;
      double rel_calc_time;
      rel_t1 = clock();
      
      ublas::matrix<double> Dmat;
      ublas::matrix<double> Dmat_1st, Dmat_2nd, Dmat_3rd, Dmat_4th;
      
      
      FE2_Macro_Solver_MCS Solve_FE2_Problem;
      LimitStateFunction_MCS TheLSF;
      
      NatafTransformation my_nataf_transformation;
      ierr = my_nataf_transformation.ModCorrMat_Empirical(probdata.num_vars,
                                                          probdata.marg,
                                                          probdata.correlation,
                                                          probdata.mod_correlation); CHKERRQ(ierr);
      //
      // Perform Cholesky decomposition for the modified correlation matrix
      //  A = LL'
      //
      ublas::triangular_matrix<double, ublas::lower> Lmat(probdata.num_vars,probdata.num_vars);
      cholesky_decompose(probdata.mod_correlation,Lmat);
      probdata.Lo.resize(probdata.num_vars,probdata.num_vars);
      probdata.Lo = Lmat;
      
      //
      // Compute the inverse of Lo
      //
      probdata.inv_Lo.resize(probdata.num_vars,probdata.num_vars);
      bool singular = false;
      probdata.inv_Lo = gjinverse(probdata.Lo,singular);
      
      ublas::vector<double> detj(FC.size());
      ublas::vector<double> pf(FC.size());
      ublas::vector<double> vpf(FC.size());
      for (unsigned i=0; i<FC.size(); i++) {
        detj(i) = 1.0;
        pf(i)   = 0.0;
        vpf(i)  = 0.0;
      }
      
      // Obj_Rel_Res.prob_failure_MCIS = Obj_Rel_Res.prob_failure_FORM;
      
      boost::mt19937 rng(time(0));
      boost::random::normal_distribution<> my_norm(0,1);
      
      ofstream MCIS_File;
      MCIS_File.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//04_ReliabilityAnalysis//Result_Multi_MCIS.txt",ofstream::out);
      
      for (int imcs=1; imcs<=no_mcs; imcs++) {
        
        cout<<"\n\n*************************************************\n*\n";
        cout<<"*    This is "<<imcs<<" step!"<<"\n*\n";
        cout<<"*************************************************\n";
        
        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //        STEP 2: generate samples                                        //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        
        //
        // Generate random points in standard normal space
        //
        u.resize(probdata.num_vars); u.clear();
        for (int irv = 0; irv<probdata.num_vars; irv++) {
          boost::variate_generator <boost::mt19937&,boost::random::normal_distribution<> > normrnd(rng,my_norm);
          u(irv) = normrnd();cout<<irv+1<<"-th rv is "<<u(irv)<<endl;
        }
        cout<<u<<endl;
        
        for (unsigned ifc = 1; ifc<=FC.size();ifc++) {
          //
          // Move points w.r.t. design point
          //
          u_temp.resize(probdata.num_vars); u_temp.clear();
          u_temp = u + Obj_Rel_Res[ifc-1].DesignPoint;
          //
          // Transform random number from standard normal space to original space
          //
          x.resize(probdata.num_vars); x.clear();
          ierr = my_nataf_transformation.u_to_x(u_temp,probdata.num_vars,probdata.marg,probdata.Lo,x,detj(ifc-1)); CHKERRQ(ierr);
          cout<<"\nRandom vaiables in U-space: "<<u<<endl;
          cout<<"\nRandom vaiables in X-space: "<<x<<endl;
          // ----------------
          // Update ply angle
          // ----------------
          for (unsigned i=1; i<=x.size();i++) {
            if (probdata.NameVars[i].compare(0,11,"orientation") == 0) {
              cout<<"\nAngle "<<x(i-1)<<endl;
              for (unsigned j=0; j<probdata.PlyAngle.size(); j++) {
                cout << "The original ply angle "<<probdata.PlyAngle(j)<<"\t delta x: "<<x(i-1)<<endl;
                PlyAngle_new(j) = probdata.PlyAngle(j) + x(i-1);
                cout << "The modified ply angle "<<PlyAngle_new(j)<<endl;
              }
            }
            else if (probdata.NameVars[i].compare(0,6,"theta1") == 0) {
              PlyAngle_new(0) = x(i-1);
            }
            else if (probdata.NameVars[i].compare(0,6,"theta2") == 0) {
              PlyAngle_new(1) = x(i-1);
            }
            else if (probdata.NameVars[i].compare(0,6,"theta3") == 0) {
              PlyAngle_new(2) = x(i-1);
            }
            else if (probdata.NameVars[i].compare(0,6,"theta4") == 0) {
              PlyAngle_new(3) = x(i-1);
            }
          }
          
          ////////////////////////////////////////////////////////////////////////////
          //                                                                        //
          //        STEP 3: call FE program to calculate structural response        //
          //                                                                        //
          ////////////////////////////////////////////////////////////////////////////

          //
          // Evaluate limit-state function and its gradient
          //
          ierr = Solve_FE2_Problem.Calculate_RVEDmat(m_field_RVE,x,probdata.num_vars,
                                                     probdata.NameVars); CHKERRQ(ierr);
          ierr = Solve_FE2_Problem.Macro_FE_Laminate(m_field_Macro,x,probdata.num_vars,
                                                     probdata.NameVars,
                                                     PlyAngle_new,
                                                     NO_Layers); CHKERRQ(ierr);
          //
          // Get constitutive matrix or D matrix at Gauss points
          //
          switch (probdata.ExaminedPly) {
            case 0:{cout<<"\n\nAll plies are under examination.\n\n";
              //--------------------------
              // D matrix for all plies
              //--------------------------
              Dmat_1st = Solve_FE2_Problem.Dmat_1st_Ply;
              Dmat_2nd = Solve_FE2_Problem.Dmat_2nd_Ply;
              Dmat_3rd = Solve_FE2_Problem.Dmat_3rd_Ply;
              Dmat_4th = Solve_FE2_Problem.Dmat_4th_Ply;
              break;
            }
            case 1:{cout<<"\n\nThe 1st layer is under examination.\n\n";
              //--------------------------
              // D matrix for the 1st ply
              //--------------------------
              Dmat = Solve_FE2_Problem.Dmat_1st_Ply;
              break;
            }
            case 2: {cout<<"\n\nThe 2nd layer is under examination.\n\n";
              //--------------------------
              // D matrix for the 2nd ply
              //--------------------------
              Dmat = Solve_FE2_Problem.Dmat_2nd_Ply;
              break;
            }
            case 3: {cout<<"\n\nThe 3rd layer is under examination.\n\n";
              //--------------------------
              // D matrix for the 3rd ply
              //--------------------------
              Dmat = Solve_FE2_Problem.Dmat_3rd_Ply;
              break;
            }
            case 4: {cout<<"\n\nThe 4th layer is under examination.\n\n";
              //--------------------------
              // D matrix for the 4th ply
              //--------------------------
              Dmat = Solve_FE2_Problem.Dmat_4th_Ply;
              break;
            }
          }
          //
          // Calculate the stress at Gauss points for specific element(s)
          //
          if (probdata.ExaminedPly == 0) {
            FE2_PostProcStressForReliability_Zeroth Calc_Stress_1st(m_field_Macro,"DISP_MACRO",Dmat_1st);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_1st",Calc_Stress_1st);  CHKERRQ(ierr);
            
            FE2_PostProcStressForReliability_Zeroth Calc_Stress_2nd(m_field_Macro,"DISP_MACRO",Dmat_2nd);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_2nd",Calc_Stress_2nd);  CHKERRQ(ierr);
            
            FE2_PostProcStressForReliability_Zeroth Calc_Stress_3rd(m_field_Macro,"DISP_MACRO",Dmat_3rd);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_3rd",Calc_Stress_3rd);  CHKERRQ(ierr);
            
            FE2_PostProcStressForReliability_Zeroth Calc_Stress_4th(m_field_Macro,"DISP_MACRO",Dmat_4th);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_4th",Calc_Stress_4th);  CHKERRQ(ierr);
            
            theta_angle = PlyAngle_new(0)*(M_PI/180.0);
            StressGP_1st.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_1st.StressGP,StressGP_1st);
            cout<<"1st ply stress at GP: "<<StressGP_1st<<endl;
            theta_angle = PlyAngle_new(1)*(M_PI/180.0);
            StressGP_2nd.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_2nd.StressGP,StressGP_2nd);
            cout<<"2nd ply stress at GP: "<<StressGP_2nd<<endl;
            theta_angle = PlyAngle_new(2)*(M_PI/180.0);
            StressGP_3rd.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_3rd.StressGP,StressGP_3rd);
            cout<<"3rd ply stress at GP: "<<StressGP_3rd<<endl;
            theta_angle = PlyAngle_new(3)*(M_PI/180.0);
            StressGP_4th.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_4th.StressGP,StressGP_4th);
            cout<<"4th ply stress at GP: "<<StressGP_4th<<endl;
          }
          else {
            FE2_PostProcStressForReliability_Zeroth Calc_Stress(m_field_Macro,"DISP_MACRO",Dmat);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress);  CHKERRQ(ierr);
            
            theta_angle = PlyAngle_new(probdata.ExaminedPly-1)*(M_PI/180.0);
            StressGP.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress.StressGP,StressGP);
            cout<<"1st ply stress at GP: "<<StressGP<<endl;
          }
          
          ////////////////////////////////////////////////////////////////////////////
          //                                                                        //
          //        STEP 4: evaluate limit state function                           //
          //                                                                        //
          ////////////////////////////////////////////////////////////////////////////
          
          switch (FC(ifc-1)) {
            case 12013: {
              //
              // Maximum stress theory - fibre failure
              //
              Obj_Rel_Res[ifc-1].NameOfFailureCriterion = "Maximum stress theory - Fibre failure";
              ierr = TheLSF.gfun_ply_MS_LD(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_MS_LD); CHKERRQ(ierr);
              cout<<"LSF value is: "<<val_G_MS_LD<<endl;
              Update_Pf_MCIS(u, val_G_MS_LD, pf(ifc-1), vpf(ifc-1), imcs, Obj_Rel_Res[ifc-1]);
              break;
            }
            case 22013: {
              //
              // Maximum stress theory - matrix failure
              //
              Obj_Rel_Res[ifc-1].NameOfFailureCriterion = "Maximum stress theory - Matrix failure";
              ierr = TheLSF.gfun_ply_MS_TD(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_MS_TD); CHKERRQ(ierr);
              cout<<"LSF value is: "<<val_G_MS_TD<<endl;
              Update_Pf_MCIS(u, val_G_MS_TD, pf(ifc-1), vpf(ifc-1), imcs, Obj_Rel_Res[ifc-1]);
              break;
            }
            case 22014: {
              //
              // Maximum stress theory - shear failure
              //
              Obj_Rel_Res[ifc-1].NameOfFailureCriterion = "Maximum stress theory - Shear failure";
              ierr = TheLSF.gfun_ply_MS_Shear(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_MS_Shear); CHKERRQ(ierr);
              cout<<"LSF value is: "<<val_G_MS_Shear<<endl;
              Update_Pf_MCIS(u, val_G_MS_Shear, pf(ifc-1), vpf(ifc-1), imcs, Obj_Rel_Res[ifc-1]);
              break;
            }
            case 13033: {
              //
              // Hashin failure theory - fibre failure
              //
              Obj_Rel_Res[ifc-1].NameOfFailureCriterion = "Hashin failure theory - Fibre failure";
              ierr = TheLSF.gfun_ply_HF(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_HF); CHKERRQ(ierr);
              cout<<"LSF value is: "<<val_G_HF<<endl;
              Update_Pf_MCIS(u, val_G_HF, pf(ifc-1), vpf(ifc-1), imcs, Obj_Rel_Res[ifc-1]);
              break;
            }
            case 23033: {
              //
              // Hashin failure theory - matrix failure
              //
              Obj_Rel_Res[ifc-1].NameOfFailureCriterion = "Hashin failure theory - Matrix failure";
              ierr = TheLSF.gfun_ply_HM(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_HM); CHKERRQ(ierr);
              cout<<"LSF value is: "<<val_G_HM<<endl;
              Update_Pf_MCIS(u, val_G_HM, pf(ifc-1), vpf(ifc-1), imcs, Obj_Rel_Res[ifc-1]);
              break;
            }
            case 43050: {
              //
              // Tsai-Wu failure criteria
              //
              Obj_Rel_Res[ifc-1].NameOfFailureCriterion = "Tsai-Wu";
              ierr = TheLSF.gfun_ply_Tsai_Wu(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_TW); CHKERRQ(ierr);
              cout<<"LSF value is: "<<val_G_TW<<endl;
              Update_Pf_MCIS(u, val_G_TW, pf(ifc-1), vpf(ifc-1), imcs, Obj_Rel_Res[ifc-1]);
              break;
            }
            case 43060: {
              //
              // Tsai-Hill failure criteria
              //
              Obj_Rel_Res[ifc-1].NameOfFailureCriterion = "Tsai-Hill - 3D stress-state";
              ierr = TheLSF.gfun_ply_Tsai_Hill(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_TH); CHKERRQ(ierr);
              cout<<"LSF value is: "<<val_G_TH<<endl;
              Update_Pf_MCIS(u, val_G_TH, pf(ifc-1), vpf(ifc-1), imcs, Obj_Rel_Res[ifc-1]);
              break;
            }
            case 13073: {
              //
              // Richard Christen: Fibre controlled failure
              //
              Obj_Rel_Res[ifc-1].NameOfFailureCriterion = "Christensen - Fibre controlled failure";
              ierr = TheLSF.gfun_ply_RCF(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_RCF); CHKERRQ(ierr);
              cout<<"LSF value is: "<<val_G_RCF<<endl;
              Update_Pf_MCIS(u, val_G_RCF, pf(ifc-1), vpf(ifc-1), imcs, Obj_Rel_Res[ifc-1]);
              break;
            }
            case 23073: {
              //
              // Richard Christen: Matrix controlled failure
              //
              Obj_Rel_Res[ifc-1].NameOfFailureCriterion = "Christensen - Matrix controlled failure";
              ierr = TheLSF.gfun_ply_RCM(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_RCM); CHKERRQ(ierr);
              cout<<"LSF value is: "<<val_G_RCM<<endl;
              Update_Pf_MCIS(u, val_G_RCM, pf(ifc-1), vpf(ifc-1), imcs, Obj_Rel_Res[ifc-1]);
              break;
            }
            case 43080: {
              //
              // Hoffman failure theory
              //
              Obj_Rel_Res[ifc-1].NameOfFailureCriterion = "Hoffman failure theory - 3D stress states";
              ierr = TheLSF.gfun_ply_Hoffman(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_Hoffman); CHKERRQ(ierr);
              cout<<"LSF value is: "<<val_G_Hoffman<<endl;
              Update_Pf_MCIS(u, val_G_Hoffman, pf(ifc-1), vpf(ifc-1), imcs, Obj_Rel_Res[ifc-1]);
              break;
            }
          }
          
          ////////////////////////////////////////////////////////////////////////////
          //                                                                        //
          //        STEP 5: write results in text file                              //
          //                display results                                         //
          //                                                                        //
          ////////////////////////////////////////////////////////////////////////////
          
          cout<<"\n\n";
          cout<<"/////////////////////////////////////////////////////////////////////////////////"<<endl;
          cout<<"//"<<endl;
          cout<<"//   The failure criterion used is: "<<Obj_Rel_Res[ifc-1].NameOfFailureCriterion<<endl;
          cout<<"//   Failure probability from MCIS: "<<Obj_Rel_Res[ifc-1].prob_failure_MCIS<<endl;
          cout<<"//   Failure probability from FORM: "<<Obj_Rel_Res[ifc-1].prob_failure_FORM<<endl;
          cout<<"//"<<endl;
          cout<<"/////////////////////////////////////////////////////////////////////////////////"<<endl;
          cout<<"\n\n";
          
          MCIS_File<<setprecision(15)<<Obj_Rel_Res[ifc-1].prob_failure_MCIS<<"\t";
          MCIS_File<<setprecision(15)<<Obj_Rel_Res[ifc-1].beta_MCIS<<"\t";
          MCIS_File<<setprecision(15)<<Obj_Rel_Res[ifc-1].prob_failure_FORM<<"\t";
          MCIS_File<<setprecision(15)<<Obj_Rel_Res[ifc-1].beta_FORM<<"\t";
        }
        
        MCIS_File<<"\n";
        
        rel_t2 = clock();
        rel_calc_time  = (double)(rel_t2 - rel_t1)/CLOCKS_PER_SEC;
        rel_t1 = rel_t2;
        cout<<"Elapsed time at this step is: "<<rel_calc_time<<" seconds.\n";
      }
      
      //////////////////////////////////////////////////////////////////////////////
      //                                                                          //
      //        STEP 5: Calculate the probability of failure and corresponding    //
      //                estimate of reliability index                             //
      //                                                                          //
      //////////////////////////////////////////////////////////////////////////////
      
      
      
      //////////////////////////////////////////////////////////////////////////////
      //                                                                          //
      //                               FINISH                                     //
      //                                                                          //
      //////////////////////////////////////////////////////////////////////////////
      
      MCIS_File.close();
      
      cout<<"\n\n*************************************************\n*\n";
      cout<<"*  Reliability method is Importance Sampling."<<endl;
      //cout<<"*  Reliability index (beta) is:    "<<MsFE_Reliability_Results.beta_MCIS<<endl;
      //cout<<"*  Probability of failure (pf) is: "<<setprecision(15)<<MsFE_Reliability_Results.prob_failure_MCIS<<endl;
      //cout<<"*  The failure criterion used is:  "<<MsFE_Reliability_Results.NameOfFailureCriterion<<endl;
      
      PetscFunctionReturn(0);
    }
    
    //==========================================================================
    //
    // Function: Updating probability of failure for Importance Sampling
    //
    //==========================================================================
    virtual PetscErrorCode Update_Pf_MCIS(ublas::vector<double> u, double g_val,
                                          double &pf, double &vpf,
                                          int iter_prob_MCIS,
                                          Reliability_Results &MsFE_Reliability_Results) {
      PetscFunctionBegin;
      
      ublas::vector<double> u_star = MsFE_Reliability_Results.DesignPoint;
      double beta_FORM = MsFE_Reliability_Results.beta_FORM;
      
      if (g_val<0) {
        pf = pf + exp(-inner_prod(u,u_star) - pow(beta_FORM,2)/2);
        vpf = vpf + pow(exp(-inner_prod(u,u_star) - pow(beta_FORM,2)/2),2);
      }
      cout<<"Sub-pf is: "<<pf<<endl;
      MsFE_Reliability_Results.prob_failure_MCIS = max(pf/iter_prob_MCIS,1e-20);
      //vpf1 = 1/iter_prob_MCIS*(vpf/iter_prob_MCIS - Pf1^2);
      
      using boost::math::normal_distribution;
      normal_distribution<> snorm(0,1);
      MsFE_Reliability_Results.beta_MCIS = -1*quantile(snorm,MsFE_Reliability_Results.prob_failure_MCIS);
      
      PetscFunctionReturn(0);
    }
  };
}

#endif //__FE2_RELIABILITY_METHODS_HPP