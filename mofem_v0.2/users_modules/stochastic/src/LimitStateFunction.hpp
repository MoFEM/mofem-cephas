/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * LIMITSTATEFUNCTION: This structure defines limit state function (LSF) used to
 *   conduct structural reliability analysis. In general, it includes two
 *   classes of LSF, one is for structures with traditional materials, and 
 *   another is for fibre reinforced polymer (FRP) based strucutres. For FRPs, 
 *   three groups of LSF are considered according to failure criteria of 
 *   composites. There are three categories of failure criteria for chioce.
 *   (1) The so called micromechanics level takes the individual fibres and the
 *       seprate matrix phase in between them as the size scale for homogeneity.
 *   (2) The next level up is the aligned fibre, lamina level, which then is much
 *       larger than the size of the individual filament or fibre.
 *   (3) Finally, at yet a still much larger scale, the homogenization could be
 *       taken at the laminate level, involving the stacking of various lamina
 *       in various directions.
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

#ifndef __LIMITSTATEFUNCTION_HPP__
#define __LIMITSTATEFUNCTION_HPP__

namespace MoFEM {
  struct LimitStateFunction_Composite_Microscale {
  };
  
  struct LimitStateFunction_Composite_Lamina {
  };
  
  struct LimitStateFunction_Composite_Laminate {
  };
  
  struct LimitStateFunction_MCS {
    
    //------------------------------------------------------------------------------
    // TSAI-WU
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Wu(ublas::vector<double> x,//ublas::matrix<double> PlyStrength,
                                       vector<string> vars_name,
                                       ublas::matrix<double> MatStrength,
                                       ublas::matrix<double> StressGP,
                                       double &val_lsf) {
      PetscFunctionBegin;
      
      ErrorCode rval;
      PetscErrorCode ierr;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      double F1, F2, F3, F11, F22, F12, F44, F66;

      X_T  = MatStrength(0,0);
      X_C  = MatStrength(1,0);
      Y_T  = MatStrength(2,0);
      Y_C  = MatStrength(3,0);
      S_12 = MatStrength(4,0);
      S_23 = MatStrength(5,0);
      
      // Update the values if the parameters are considered as random variables.
      for (int i=1; i<=x.size();i++) {
        if (vars_name[i].compare(0,2,"XT") == 0) {
          X_T = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          X_C = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          Y_T = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          Y_C = x(i-1);
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          S_12 = x(i-1);
        }
        else if (vars_name[i].compare(0,3,"S23") == 0) {
          S_23 = x(i-1);
        }
      }
      
      F1 = (1/X_T - 1/X_C);
      F2 = (1/Y_T - 1/Y_C);
      F11 = 1/(X_T*X_C);
      F22 = 1/(Y_T*Y_C);
      F12 = -1/sqrt(X_T*X_C*Y_T*Y_C); // Empirical suggestion: Mises-Hencky criterion
      F44 = 1/(S_23*S_23);
      F66 = 1/(S_12*S_12);
      
      
      val_lsf = 1- (  F1*StressGP(0,0)
                    + F2*(StressGP(1,1) + StressGP(2,2))
                    + F11*StressGP(0,0)*StressGP(0,0)
                    + F22*pow((StressGP(1,1) + StressGP(2,2)),2)
                    + F44*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2))
                    + F66*(StressGP(0,1)*StressGP(0,1) + StressGP(2,0)*StressGP(2,0))
                    + F12*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
      
      PetscFunctionReturn(0);
      
    }
    
    
    //------------------------------------------------------------------------------
    // Tsai-Hill
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Hill(ublas::vector<double> x,
                                              vector<string> vars_name,
                                              ublas::matrix<double> PlyStrength,
                                              ublas::matrix<double> StressGP,
                                              double &val_lsf) {
      PetscFunctionBegin;
      
      ErrorCode rval;
      PetscErrorCode ierr;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double X,Y;
      
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = PlyStrength(0,0);
      X_C  = PlyStrength(1,0);
      Y_T  = PlyStrength(2,0);
      Y_C  = PlyStrength(3,0);
      S_12 = PlyStrength(4,0);
      
      
      // Update the values if the parameters are considered as random variables.
      for (int i=1; i<=x.size();i++) {
        if (vars_name[i].compare(0,2,"XT") == 0) {
          X_T = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          X_C = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          Y_T = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          Y_C = x(i-1);
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          S_12 = x(i-1);
        }
      }
      
      // Evaluate the limit state function
      if (StressGP(0,0)>0) {
        X = X_T;
      } else {
        X = X_C;
      }
      if (StressGP(1,1)>0) {
        Y = Y_T;
      } else {
        Y = Y_C;
      }
      
      val_lsf = 1- (StressGP(0,0)*StressGP(0,0)/(X*X)
                    + StressGP(1,1)*StressGP(1,1)/(Y*Y)
                    + StressGP(0,1)*StressGP(0,1)/(S_12*S_12)
                    - StressGP(0,0)*StressGP(1,1)/(X*X));
      
      PetscFunctionReturn(0);
    }
    
  };
  
  struct LimitStateFunction {
    
    //------------------------------------------------------------------------------
    // To compute limit state function and its gradient w.r.t. the given trail
    // checking point
    //
    virtual PetscErrorCode gfun(ublas::vector<double> x,
                                double &val_lsf,
                                ublas::vector<double> &grad_lsf) {
      /*
       * Define limit state function
       *   Here, using the limit state function of Example 4.5 in
       *   [Choi et al., 2007] Reliability-based structural design
       *
       *   g = E*I - 78.12*P
       *     where E = Young's modulus, denoted as x(0)
       *           I = area moment of the cross-section, denoted as x(1)
       *           P = applied load, denoted as x(2)
       */
      
      // Evaluate the limit state function
      val_lsf = x(0)*x(1)-78.12*x(2);
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      // w.r.t. the 1st random variable
      grad_lsf(0) = x(1);
      // w.r.t. the 2nd random variable
      grad_lsf(1) = x(0);
      // w.r.t. the 3rd random variable
      grad_lsf(2) = -78.12;
    }
    
    //------------------------------------------------------------------------------
    // CHRISTENSEN: MATRIX CONTROLLED FAILURE
    //
    
    virtual PetscErrorCode gfun_composite(ublas::matrix<double> PlyStrength,
                                          ublas::matrix<double> StressGP,
                                          ublas::matrix<double> StressGP_r_Em,
                                          ublas::matrix<double> StressGP_r_NUm,
                                          ublas::matrix<double> StressGP_r_NUp,
                                          ublas::matrix<double> StressGP_r_NUpz,
                                          ublas::matrix<double> StressGP_r_Ep,
                                          ublas::matrix<double> StressGP_r_Ez,
                                          ublas::matrix<double> StressGP_r_Gzp,
                                          double &val_lsf,
                                          ublas::vector<double> &grad_lsf) {
      
      double T_11; // tensile strength in the fibre direction
      double T_22; // tensile strength in the transverse direction
      double C_11; // compressive strength in the fibre direction
      double C_22; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      T_11 = PlyStrength(0,0);
      C_11 = PlyStrength(1,0);
      T_22 = PlyStrength(2,0);
      C_22 = PlyStrength(3,0);
      S_12 = PlyStrength(4,0);
      S_23 = PlyStrength(5,0);
    
      
      // matrix controlled failure
      
      StressGP = StressGP;
      // Evaluate the limit state function
      val_lsf = 1- ((1/T_22 - 1/C_22)*(StressGP(1,1) + StressGP(2,2))
                + 1/(T_22*C_22)*(StressGP(1,1) + StressGP(2,2))*(StressGP(1,1) + StressGP(2,2))
                + 1/(S_23*S_23)*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2))
                + 1/(S_12*S_12)*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2)));
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      
      // w.r.t. Em
      dsdx.clear(); dsdx = StressGP_r_Em;
      grad_lsf(0) = (1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                    + 1/(T_22*C_22)*2*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                    + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2)-StressGP(1,1)*dsdx(2,2))
                    + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(2,0));
      
      // w.r.t. NUm
      dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
      grad_lsf(1) = (1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
      + 1/(T_22*C_22)*2*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
      + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2)-StressGP(1,1)*dsdx(2,2))
      + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(2,0));
      
      // w.r.t. NUp
      dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
      grad_lsf(2) = (1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
      + 1/(T_22*C_22)*2*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
      + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2)-StressGP(1,1)*dsdx(2,2))
      + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(2,0));
      
      // w.r.t. NUpz
      dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
      grad_lsf(3) = (1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
      + 1/(T_22*C_22)*2*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
      + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2)-StressGP(1,1)*dsdx(2,2))
      + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(2,0));
      
      // w.r.t. Ep
      dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
      grad_lsf(4) = (1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
      + 1/(T_22*C_22)*2*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
      + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2)-StressGP(1,1)*dsdx(2,2))
      + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(2,0));
      
      // w.r.t. Ez
      dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
      grad_lsf(5) = (1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
      + 1/(T_22*C_22)*2*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
      + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2)-StressGP(1,1)*dsdx(2,2))
      + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(2,0));
      
      // w.r.t. Gzp
      dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
      grad_lsf(6) = (1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
      + 1/(T_22*C_22)*2*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
      + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2)-StressGP(1,1)*dsdx(2,2))
      + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(2,0));
      
      grad_lsf = -1*grad_lsf;
    }
    
    
    //------------------------------------------------------------------------------
    // TSAI-WU
    //
    
    virtual PetscErrorCode gfun_ply_TW(ublas::vector<double> x,//ublas::matrix<double> PlyStrength,
                                       vector<string> vars_name,
                                       ublas::matrix<double> MatStrength,
                                       ublas::matrix<double> StressGP,
                                       ublas::matrix<double> StressGP_r_Em,
                                       ublas::matrix<double> StressGP_r_NUm,
                                       ublas::matrix<double> StressGP_r_NUp,
                                       ublas::matrix<double> StressGP_r_NUpz,
                                       ublas::matrix<double> StressGP_r_Ep,
                                       ublas::matrix<double> StressGP_r_Ez,
                                       ublas::matrix<double> StressGP_r_Gzp,
                                       ublas::matrix<double> StressGP_r_F,
                                       double &val_lsf,
                                       ublas::vector<double> &grad_lsf) {
      PetscFunctionBegin;
      
      ErrorCode rval;
      PetscErrorCode ierr;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      double F1, F2, F3, F11, F22, F12, F44, F66;
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      /*
       X_T  = x(7);
      X_C  = x(8);
      Y_T  = x(9);
      Y_C  = x(10);
      S_12 = x(11);
      S_23 = x(11);
       */
      X_T  = MatStrength(0,0);
      X_C  = MatStrength(1,0);
      Y_T  = MatStrength(2,0);
      Y_C  = MatStrength(3,0);cout<<"The Y_T is: "<<Y_T<<endl;
      S_12 = MatStrength(4,0);
      S_23 = MatStrength(5,0);
      
      // Update the values if the parameters are considered as random variables.
      for (int i=1; i<=x.size();i++) {
        if (vars_name[i].compare(0,2,"XT") == 0) {
          X_T = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          X_C = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          Y_T = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          Y_C = x(i-1);
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          S_12 = x(i-1);
        }
        else if (vars_name[i].compare(0,3,"S23") == 0) {
          S_23 = x(i-1);
        }
      }
      
      F1 = (1/X_T - 1/X_C);
      F2 = (1/Y_T - 1/Y_C);
      F11 = 1/(X_T*X_C);
      F22 = 1/(Y_T*Y_C);
      F12 = -1/sqrt(X_T*X_C*Y_T*Y_C); // Empirical suggestion: Mises-Hencky criterion
      F44 = 1/(S_23*S_23);
      F66 = 1/(S_12*S_12);
      
      
      val_lsf = 1- (  F1*StressGP(0,0)
                    + F2*(StressGP(1,1) + StressGP(2,2))
                    + F11*StressGP(0,0)*StressGP(0,0)
                    + F22*pow((StressGP(1,1) + StressGP(2,2)),2)
                    + F44*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2))
                    + F66*(StressGP(0,1)*StressGP(0,1) + StressGP(2,0)*StressGP(2,0))
                    + F12*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      for (int i=1; i<=x.size(); i++) {
        cout<<"The variable name is "<<vars_name[i]<<endl;
        if (vars_name[i].compare(0,2,"Em") == 0) {
          // w.r.t. Em
          dsdx.clear(); dsdx = StressGP_r_Em;
          grad_lsf(i-1) = - (  F1*StressGP_r_Em(0,0)
                             + F2*(StressGP_r_Em(1,1) + StressGP_r_Em(2,2))
                             + 2*F11*StressGP(0,0)*StressGP_r_Em(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(StressGP_r_Em(1,1) + StressGP_r_Em(2,2))
                             + F44*(2*StressGP(1,2)*StressGP_r_Em(1,2)
                                    - StressGP_r_Em(1,1)*StressGP(2,2)
                                    - StressGP(1,1)*StressGP_r_Em(2,2))
                             + F66*(2*StressGP(0,1)*StressGP_r_Em(0,1) + 2*StressGP(2,0)*StressGP_r_Em(2,0))
                             + F12*StressGP_r_Em(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + F12*StressGP(0,0)*(StressGP_r_Em(1,1)+StressGP_r_Em(2,2)));
        }
        else if (vars_name[i].compare(0,3,"NUm") == 0) {
          // w.r.t. NUm
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
          grad_lsf(i-1) = - (  F1*StressGP_r_NUm(0,0)
                           + F2*(StressGP_r_NUm(1,1) + StressGP_r_NUm(2,2))
                           + 2*F11*StressGP(0,0)*StressGP_r_NUm(0,0)
                           + 2*F22*(StressGP(1,1) + StressGP(2,2))*(StressGP_r_NUm(1,1) + StressGP_r_NUm(2,2))
                           + F44*(2*StressGP(1,2)*StressGP_r_NUm(1,2)
                                  - StressGP_r_NUm(1,1)*StressGP(2,2)
                                  - StressGP(1,1)*StressGP_r_NUm(2,2))
                           + F66*(2*StressGP(0,1)*StressGP_r_NUm(0,1) + 2*StressGP(2,0)*StressGP_r_NUm(2,0))
                           + F12*StressGP_r_NUm(0,0)*(StressGP(1,1)+StressGP(2,2))
                           + F12*StressGP(0,0)*(StressGP_r_NUm(1,1)+StressGP_r_NUm(2,2)));
        }
        else if (vars_name[i].compare(0,3,"NUp") == 0) {
          // w.r.t. NUp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
          grad_lsf(i-1) = - (  F1*StressGP_r_NUp(0,0)
                           + F2*(StressGP_r_NUp(1,1) + StressGP_r_NUp(2,2))
                           + 2*F11*StressGP(0,0)*StressGP_r_NUp(0,0)
                           + 2*F22*(StressGP(1,1) + StressGP(2,2))*(StressGP_r_NUp(1,1) + StressGP_r_NUp(2,2))
                           + F44*(2*StressGP(1,2)*StressGP_r_NUp(1,2)
                                  - StressGP_r_NUp(1,1)*StressGP(2,2)
                                  - StressGP(1,1)*StressGP_r_NUp(2,2))
                           + F66*(2*StressGP(0,1)*StressGP_r_NUp(0,1) + 2*StressGP(2,0)*StressGP_r_NUp(2,0))
                           + F12*StressGP_r_NUp(0,0)*(StressGP(1,1)+StressGP(2,2))
                           + F12*StressGP(0,0)*(StressGP_r_NUp(1,1)+StressGP_r_NUp(2,2)));
        }
        else if (vars_name[i].compare(0,4,"NUpz") == 0) {
          // w.r.t. NUpz
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
          grad_lsf(i-1) = - (  F1*StressGP_r_NUpz(0,0)
                           + F2*(StressGP_r_NUpz(1,1) + StressGP_r_NUpz(2,2))
                           + 2*F11*StressGP(0,0)*StressGP_r_NUpz(0,0)
                           + 2*F22*(StressGP(1,1) + StressGP(2,2))*(StressGP_r_NUpz(1,1) + StressGP_r_NUpz(2,2))
                           + F44*(2*StressGP(1,2)*StressGP_r_NUpz(1,2)
                                  - StressGP_r_NUpz(1,1)*StressGP(2,2)
                                  - StressGP(1,1)*StressGP_r_NUpz(2,2))
                           + F66*(2*StressGP(0,1)*StressGP_r_NUpz(0,1) + 2*StressGP(2,0)*StressGP_r_NUpz(2,0))
                           + F12*StressGP_r_NUpz(0,0)*(StressGP(1,1)+StressGP(2,2))
                           + F12*StressGP(0,0)*(StressGP_r_NUpz(1,1)+StressGP_r_NUpz(2,2)));
          
        }
        else if (vars_name[i].compare(0,2,"Ep") == 0) {
          // w.r.t. Ep
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
          grad_lsf(i-1) = - (  F1*StressGP_r_Ep(0,0)
                           + F2*(StressGP_r_Ep(1,1) + StressGP_r_Ep(2,2))
                           + 2*F11*StressGP(0,0)*StressGP_r_Ep(0,0)
                           + 2*F22*(StressGP(1,1) + StressGP(2,2))*(StressGP_r_Ep(1,1) + StressGP_r_Ep(2,2))
                           + F44*(2*StressGP(1,2)*StressGP_r_Ep(1,2)
                                  - StressGP_r_Ep(1,1)*StressGP(2,2)
                                  - StressGP(1,1)*StressGP_r_Ep(2,2))
                           + F66*(2*StressGP(0,1)*StressGP_r_Ep(0,1) + 2*StressGP(2,0)*StressGP_r_Ep(2,0))
                           + F12*StressGP_r_Ep(0,0)*(StressGP(1,1)+StressGP(2,2))
                           + F12*StressGP(0,0)*(StressGP_r_Ep(1,1)+StressGP_r_Ep(2,2)));
        }
        else if (vars_name[i].compare(0,2,"Ez") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
          grad_lsf(i-1) = - (  F1*StressGP_r_Ez(0,0)
                           + F2*(StressGP_r_Ez(1,1) + StressGP_r_Ez(2,2))
                           + 2*F11*StressGP(0,0)*StressGP_r_Ez(0,0)
                           + 2*F22*(StressGP(1,1) + StressGP(2,2))*(StressGP_r_Ez(1,1) + StressGP_r_Ez(2,2))
                           + F44*(2*StressGP(1,2)*StressGP_r_Ez(1,2)
                                  - StressGP_r_Ez(1,1)*StressGP(2,2)
                                  - StressGP(1,1)*StressGP_r_Ez(2,2))
                           + F66*(2*StressGP(0,1)*StressGP_r_Ez(0,1) + 2*StressGP(2,0)*StressGP_r_Ez(2,0))
                           + F12*StressGP_r_Ez(0,0)*(StressGP(1,1)+StressGP(2,2))
                           + F12*StressGP(0,0)*(StressGP_r_Ez(1,1)+StressGP_r_Ez(2,2)));
        }
        else if (vars_name[i].compare(0,3,"Gzp") == 0) {
          // w.r.t. Gzp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
          grad_lsf(i-1) = - (  F1*StressGP_r_Gzp(0,0)
                           + F2*(StressGP_r_Gzp(1,1) + StressGP_r_Gzp(2,2))
                           + 2*F11*StressGP(0,0)*StressGP_r_Gzp(0,0)
                           + 2*F22*(StressGP(1,1) + StressGP(2,2))*(StressGP_r_Gzp(1,1) + StressGP_r_Gzp(2,2))
                           + F44*(2*StressGP(1,2)*StressGP_r_Gzp(1,2)
                                  - StressGP_r_Gzp(1,1)*StressGP(2,2)
                                  - StressGP(1,1)*StressGP_r_Gzp(2,2))
                           + F66*(2*StressGP(0,1)*StressGP_r_Gzp(0,1) + 2*StressGP(2,0)*StressGP_r_Gzp(2,0))
                           + F12*StressGP_r_Gzp(0,0)*(StressGP(1,1)+StressGP(2,2))
                           + F12*StressGP(0,0)*(StressGP_r_Gzp(1,1)+StressGP_r_Gzp(2,2)));
        }
        else if (vars_name[i].compare(0,5,"force") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
          grad_lsf(i-1) = - (  F1*StressGP_r_F(0,0)
                             + F2*(StressGP_r_F(1,1) + StressGP_r_F(2,2))
                             + 2*F11*StressGP(0,0)*StressGP_r_F(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(StressGP_r_F(1,1) + StressGP_r_F(2,2))
                             + F44*(2*StressGP(1,2)*StressGP_r_F(1,2)
                                    - StressGP_r_F(1,1)*StressGP(2,2)
                                    - StressGP(1,1)*StressGP_r_F(2,2))
                             + F66*(2*StressGP(0,1)*StressGP_r_F(0,1) + 2*StressGP(2,0)*StressGP_r_F(2,0))
                             + F12*StressGP_r_F(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + F12*StressGP(0,0)*(StressGP_r_F(1,1)+StressGP_r_F(2,2)));
        }
        else if (vars_name[i].compare(0,2,"XT") == 0) {
          // w.r.t X_T
          double F12_XT = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_C*Y_T*Y_C;
          double F1_XT = -1/(X_T*X_T);
          double F11_XT = -1/(X_T*X_T*X_C);
          grad_lsf(i-1) = - (  F1_XT*StressGP(0,0)
                           + F11_XT*StressGP(0,0)*StressGP(0,0)
                           + 2*F12_XT*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          // w.r.t X_C
          double F1_XC = 1/(X_C*X_C);
          double F11_XC = -1/(X_T*X_C*X_C);
          double F12_XC = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*Y_T*Y_C;
          grad_lsf(i-1) = - (  F1_XC*StressGP(0,0)
                           + F11_XC*StressGP(0,0)*StressGP(0,0)
                           + 2*F12_XC*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          // w.r.t Y_T
          double F2_YT = -1/(Y_T*Y_T);
          double F22_YT = -1/(Y_T*Y_T*Y_C);
          double F12_YT = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*X_C*Y_C;
          grad_lsf(i-1) = - ( F2_YT*(StressGP(1,1) + StressGP(2,2))
                           + F22_YT*pow((StressGP(1,1) + StressGP(2,2)),2)
                           + 2*F12_YT*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          // w.r.t Y_C
          double F2_YC = 1/(Y_C*Y_C);
          double F22_YC = -1/(Y_T*Y_C*Y_C);
          double F12_YC = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*X_C*Y_T;
          grad_lsf(i-1) = - ( F2_YC*(StressGP(1,1) + StressGP(2,2))
                            + F22_YC*pow((StressGP(1,1) + StressGP(2,2)),2)
                            + 2*F12_YC*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          // w.r.t S_12
          grad_lsf(i-1) = - (-2*StressGP(0,1)*StressGP(0,1)/pow(S_12,3));
        }
      
      }
      
      PetscFunctionReturn(0);
      
    }
    
    //--------------------------------------------------------------------------
    // TSAI-WU: reduction version for the case of an orthotropic lamina under
    //          plane stress conditions
    //
    
    virtual PetscErrorCode gfun_ply_TW_PlaneStress(ublas::vector<double> x,//ublas::matrix<double> PlyStrength,
                                                         ublas::matrix<double> StressGP,
                                                         ublas::matrix<double> StressGP_r_Em,
                                                         ublas::matrix<double> StressGP_r_NUm,
                                                         ublas::matrix<double> StressGP_r_NUp,
                                                         ublas::matrix<double> StressGP_r_NUpz,
                                                         ublas::matrix<double> StressGP_r_Ep,
                                                         ublas::matrix<double> StressGP_r_Ez,
                                                         ublas::matrix<double> StressGP_r_Gzp,
                                                         double &val_lsf,
                                                         ublas::vector<double> &grad_lsf) {
      PetscFunctionBegin;
      
      ErrorCode rval;
      PetscErrorCode ierr;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      double F1, F2, F3, F11, F22, F12, F44, F66;
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = x(7);
      X_C  = x(8);
      Y_T  = x(9);
      Y_C  = x(10);
      S_12 = x(11);
      S_23 = x(11);
      
      F1 = (1/X_T - 1/X_C);
      F2 = (1/Y_T - 1/Y_C);
      F11 = 1/(X_T*X_C);
      F22 = 1/(Y_T*Y_C);
      F12 = -1/2*sqrt(1/(X_T*X_C*Y_T*Y_C)); // Empirical suggestion: Mises-Hencky criterion
      F44 = 1/(S_23*S_23);
      F66 = 1/(S_12*S_12);
      
      //StressGP = StressGP;
      // Evaluate the limit state function
      val_lsf = 1- (  F1*StressGP(0,0)
                    + F2*StressGP(1,1)
                    + F11*StressGP(0,0)*StressGP(0,0)
                    + F22*StressGP(1,1)*StressGP(1,1)
                    + F66*StressGP(0,1)*StressGP(0,1)
                    + 2*F12*StressGP(0,0)*StressGP(1,1));
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      
      // w.r.t. Em
      dsdx.clear(); dsdx = StressGP_r_Em;
      grad_lsf(0) = - (  F1*StressGP_r_Em(0,0)
                       + F2*StressGP_r_Em(1,1)
                       + 2*F11*StressGP(0,0)*StressGP_r_Em(0,0)
                       + 2*F22*StressGP(1,1)*StressGP_r_Em(1,1)
                       + 2*F66*StressGP(0,1)*StressGP_r_Em(0,1)
                       + 2*F12*StressGP_r_Em(0,0)*StressGP(1,1)
                       + 2*F12*StressGP(0,0)*StressGP_r_Em(1,1));
      
      // w.r.t. NUm
      dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
      grad_lsf(1) = - (  F1*StressGP_r_NUm(0,0)
                       + F2*StressGP_r_NUm(1,1)
                       + 2*F11*StressGP(0,0)*StressGP_r_NUm(0,0)
                       + 2*F22*StressGP(1,1)*StressGP_r_NUm(1,1)
                       + 2*F66*StressGP(0,1)*StressGP_r_NUm(0,1)
                       + 2*F12*StressGP_r_NUm(0,0)*StressGP(1,1)
                       + 2*F12*StressGP(0,0)*StressGP_r_NUm(1,1));
      
      
      // w.r.t. NUp
      dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
      grad_lsf(2) = - (  F1*StressGP_r_NUp(0,0)
                       + F2*StressGP_r_NUp(1,1)
                       + 2*F11*StressGP(0,0)*StressGP_r_NUp(0,0)
                       + 2*F22*StressGP(1,1)*StressGP_r_NUp(1,1)
                       + 2*F66*StressGP(0,1)*StressGP_r_NUp(0,1)
                       + 2*F12*StressGP_r_NUp(0,0)*StressGP(1,1)
                       + 2*F12*StressGP(0,0)*StressGP_r_NUp(1,1));
      
      // w.r.t. NUpz
      dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
      grad_lsf(3) = - (  F1*StressGP_r_NUp(0,0)
                       + F2*StressGP_r_NUp(1,1)
                       + 2*F11*StressGP(0,0)*StressGP_r_NUp(0,0)
                       + 2*F22*StressGP(1,1)*StressGP_r_NUp(1,1)
                       + 2*F66*StressGP(0,1)*StressGP_r_NUp(0,1)
                       + 2*F12*StressGP_r_NUp(0,0)*StressGP(1,1)
                       + 2*F12*StressGP(0,0)*StressGP_r_NUp(1,1));
      
      // w.r.t. Ep
      dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
      grad_lsf(4) = - (  F1*StressGP_r_Ep(0,0)
                       + F2*StressGP_r_Ep(1,1)
                       + 2*F11*StressGP(0,0)*StressGP_r_Ep(0,0)
                       + 2*F22*StressGP(1,1)*StressGP_r_Ep(1,1)
                       + 2*F66*StressGP(0,1)*StressGP_r_Ep(0,1)
                       + 2*F12*StressGP_r_Ep(0,0)*StressGP(1,1)
                       + 2*F12*StressGP(0,0)*StressGP_r_Ep(1,1));
      
      // w.r.t. Ez
      dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
      grad_lsf(5) = - (  F1*StressGP_r_Ez(0,0)
                       + F2*StressGP_r_Ez(1,1)
                       + 2*F11*StressGP(0,0)*StressGP_r_Ez(0,0)
                       + 2*F22*StressGP(1,1)*StressGP_r_Ez(1,1)
                       + 2*F66*StressGP(0,1)*StressGP_r_Ez(0,1)
                       + 2*F12*StressGP_r_Ez(0,0)*StressGP(1,1)
                       + 2*F12*StressGP(0,0)*StressGP_r_Ez(1,1));
      
      // w.r.t. Gzp
      dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
      grad_lsf(6) = - (  F1*StressGP_r_Gzp(0,0)
                       + F2*StressGP_r_Gzp(1,1)
                       + 2*F11*StressGP(0,0)*StressGP_r_Gzp(0,0)
                       + 2*F22*StressGP(1,1)*StressGP_r_Gzp(1,1)
                       + 2*F66*StressGP(0,1)*StressGP_r_Gzp(0,1)
                       + 2*F12*StressGP_r_Gzp(0,0)*StressGP(1,1)
                       + 2*F12*StressGP(0,0)*StressGP_r_Gzp(1,1));
      
      // w.r.t X_T
      double F12_XT = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_C*Y_T*Y_C;
      double F1_XT = -1/(X_T*X_T);
      double F11_XT = -1/(X_T*X_T*X_C);
      grad_lsf(7) = - (  F1_XT*StressGP(0,0)
                       + F11_XT*StressGP(0,0)*StressGP(0,0)
                       + 2*F12_XT*StressGP(0,0)*StressGP(1,1));
      // w.r.t X_C
      double F1_XC = 1/(X_C*X_C);
      double F11_XC = -1/(X_T*X_C*X_C);
      double F12_XC = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*Y_T*Y_C;
      grad_lsf(8) = - (  F1_XC*StressGP(0,0)
                       + F11_XC*StressGP(0,0)*StressGP(0,0)
                       + 2*F12_XC*StressGP(0,0)*StressGP(1,1));
      // w.r.t Y_T
      double F2_YT = -1/(Y_T*Y_T);
      double F22_YT = -1/(Y_T*Y_T*Y_C);
      double F12_YT = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*X_C*Y_C;
      grad_lsf(9) = - (  F2_YT*StressGP(1,1)
                       + F22_YT*StressGP(1,1)*StressGP(1,1)
                       + 2*F12_YT*StressGP(0,0)*StressGP(1,1));
      // w.r.t Y_C
      double F2_YC = 1/(Y_C*Y_C);
      double F22_YC = -1/(Y_T*Y_C*Y_C);
      double F12_YC = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*X_C*Y_T;
      grad_lsf(10) = - ( F2_YC*StressGP(1,1)
                        + F22_YC*StressGP(1,1)*StressGP(1,1)
                        + 2*F12_YC*StressGP(0,0)*StressGP(1,1));
      // w.r.t S_12
      grad_lsf(11) = - (-2*StressGP(0,1)*StressGP(0,1)/pow(S_12,3));
      
      PetscFunctionReturn(0);
    }
    
    
    //------------------------------------------------------------------------------
    // Tsai-Hill
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Hill(ublas::vector<double> x,
                                              vector<string> vars_name,
                                              ublas::matrix<double> PlyStrength,
                                              ublas::matrix<double> StressGP,
                                              ublas::matrix<double> StressGP_r_Em,
                                              ublas::matrix<double> StressGP_r_NUm,
                                              ublas::matrix<double> StressGP_r_NUp,
                                              ublas::matrix<double> StressGP_r_NUpz,
                                              ublas::matrix<double> StressGP_r_Ep,
                                              ublas::matrix<double> StressGP_r_Ez,
                                              ublas::matrix<double> StressGP_r_Gzp,
                                              ublas::matrix<double> StressGP_r_F,
                                              double &val_lsf,
                                              ublas::vector<double> &grad_lsf) {
      
      PetscFunctionBegin;
      
      ErrorCode rval;
      PetscErrorCode ierr;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double X,Y;
      
      
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = PlyStrength(0,0);
      X_C  = PlyStrength(1,0);
      Y_T  = PlyStrength(2,0);
      Y_C  = PlyStrength(3,0);
      S_12 = PlyStrength(4,0);
      
      
      // Update the values if the parameters are considered as random variables.
      for (int i=1; i<=x.size();i++) {
        if (vars_name[i].compare(0,2,"XT") == 0) {
          X_T = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          X_C = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          Y_T = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          Y_C = x(i-1);
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          S_12 = x(i-1);
        }
      }
      
      // Evaluate the limit state function
      if (StressGP(0,0)>0) {
        X = X_T;
      } else {
        X = X_C;
      }
      if (StressGP(1,1)>0) {
        Y = Y_T;
      } else {
        Y = Y_C;
      }
      
      val_lsf = 1- (StressGP(0,0)*StressGP(0,0)/(X*X)
                    + StressGP(1,1)*StressGP(1,1)/(Y*Y)
                    + StressGP(0,1)*StressGP(0,1)/(S_12*S_12)
                    - StressGP(0,0)*StressGP(1,1)/(X*X));
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      for (int i=1; i<=x.size(); i++) {
        cout<<"The variable name is "<<vars_name[i]<<endl;
        if (vars_name[i].compare(0,2,"Em") == 0) {
          // w.r.t. Em
          dsdx.clear(); dsdx = StressGP_r_Em;
          grad_lsf(i-1) = - (2*StressGP(0,0)*StressGP_r_Em(0,0)/(X*X)
                           + 2*StressGP(1,1)*StressGP_r_Em(1,1)/(Y*Y)
                           + 2*StressGP(0,1)*StressGP_r_Em(0,1)/(S_12*S_12)
                           - StressGP_r_Em(0,0)*StressGP(1,1)/(X*X)
                           - StressGP(0,0)*StressGP_r_Em(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,3,"NUm") == 0) {
          // w.r.t. NUm
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
          grad_lsf(i-1) = - (2*StressGP(0,0)*StressGP_r_NUm(0,0)/(X*X)
                           + 2*StressGP(1,1)*StressGP_r_NUm(1,1)/(Y*Y)
                           + 2*StressGP(0,1)*StressGP_r_NUm(0,1)/(S_12*S_12)
                           - StressGP_r_NUm(0,0)*StressGP(1,1)/(X*X)
                           - StressGP(0,0)*StressGP_r_NUm(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,3,"NUp") == 0) {
          // w.r.t. NUp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
          grad_lsf(i-1) = - (2*StressGP(0,0)*StressGP_r_NUp(0,0)/(X*X)
                           + 2*StressGP(1,1)*StressGP_r_NUp(1,1)/(Y*Y)
                           + 2*StressGP(0,1)*StressGP_r_NUp(0,1)/(S_12*S_12)
                           - StressGP_r_NUp(0,0)*StressGP(1,1)/(X*X)
                           - StressGP(0,0)*StressGP_r_NUp(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,4,"NUpz") == 0) {
          // w.r.t. NUpz
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
          grad_lsf(i-1) = - (2*StressGP(0,0)*StressGP_r_NUpz(0,0)/(X*X)
                           + 2*StressGP(1,1)*StressGP_r_NUpz(1,1)/(Y*Y)
                           + 2*StressGP(0,1)*StressGP_r_NUpz(0,1)/(S_12*S_12)
                           - StressGP_r_NUpz(0,0)*StressGP(1,1)/(X*X)
                           - StressGP(0,0)*StressGP_r_NUpz(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,2,"Ep") == 0) {
          // w.r.t. Ep
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
          grad_lsf(i-1) = - (2*StressGP(0,0)*StressGP_r_Ep(0,0)/(X*X)
                           + 2*StressGP(1,1)*StressGP_r_Ep(1,1)/(Y*Y)
                           + 2*StressGP(0,1)*StressGP_r_Ep(0,1)/(S_12*S_12)
                           - StressGP_r_Ep(0,0)*StressGP(1,1)/(X*X)
                           - StressGP(0,0)*StressGP_r_Ep(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,2,"Ez") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
          grad_lsf(i-1) = - (2*StressGP(0,0)*StressGP_r_Ez(0,0)/(X*X)
                           + 2*StressGP(1,1)*StressGP_r_Ez(1,1)/(Y*Y)
                           + 2*StressGP(0,1)*StressGP_r_Ez(0,1)/(S_12*S_12)
                           - StressGP_r_Ez(0,0)*StressGP(1,1)/(X*X)
                           - StressGP(0,0)*StressGP_r_Ez(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,2,"Ez") == 0) {
          // w.r.t. Gzp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
          grad_lsf(i-1) = - (2*StressGP(0,0)*StressGP_r_Gzp(0,0)/(X*X)
                           + 2*StressGP(1,1)*StressGP_r_Gzp(1,1)/(Y*Y)
                           + 2*StressGP(0,1)*StressGP_r_Gzp(0,1)/(S_12*S_12)
                           - StressGP_r_Gzp(0,0)*StressGP(1,1)/(X*X)
                           - StressGP(0,0)*StressGP_r_Gzp(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,2,"XT") == 0) {
          // w.r.t. XT
          if (StressGP(0,0)>0) {
            grad_lsf(i-1) =  2*StressGP(0,0)*StressGP(0,0)/(X*X*X)
                            -2*StressGP(0,0)*StressGP(1,1)/(X*X*X);
          } else {
            grad_lsf(i-1) = 0;
          }
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          // w.r.t. XC
          if (StressGP(0,0)<0) {
            grad_lsf(i-1) =  2*StressGP(0,0)*StressGP(0,0)/(X*X*X)
                            -2*StressGP(0,0)*StressGP(1,1)/(X*X*X);
          } else {
            grad_lsf(i-1) = 0;
          }
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          // w.r.t. YT
          if (StressGP(1,1)>0) {
            grad_lsf(i-1) = 2*StressGP(1,1)*StressGP(1,1)/(Y*Y*Y);
          } else {
            grad_lsf(i-1) = 0;
          }
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          // w.r.t. YC
          if (StressGP(1,1)<0) {
            grad_lsf(i-1) = 2*StressGP(1,1)*StressGP(1,1)/(Y*Y*Y);
          } else {
            grad_lsf(i-1) = 0;
          }
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          // w.r.t. S12
          grad_lsf(i-1) = 2*StressGP(0,1)*StressGP(0,1)/(S_12*S_12*S_12);
        }
        else if (vars_name[i].compare(0,5,"force") == 0) {
          // w.r.t. Em
          dsdx.clear(); dsdx = StressGP_r_F;
          grad_lsf(i-1) = - (2*StressGP(0,0)*StressGP_r_F(0,0)/(X*X)
                             + 2*StressGP(1,1)*StressGP_r_F(1,1)/(Y*Y)
                             + 2*StressGP(0,1)*StressGP_r_F(0,1)/(S_12*S_12)
                             - StressGP_r_F(0,0)*StressGP(1,1)/(X*X)
                             - StressGP(0,0)*StressGP_r_F(1,1)/(X*X));
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Hashin criterion
    //
    
  };
  
}

#endif //__LIMITSTATEFUNCTION_HPP__
