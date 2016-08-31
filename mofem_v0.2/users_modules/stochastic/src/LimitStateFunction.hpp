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
  // -------------------------------------------------------------------------
  //
  // Code convention: A B C D E
  // A: Failure object
  //    1. Fibre failure
  //    2. Matrix failure
  //    3. Interface
  //    4. Ply failure
  
  // B: required stress state
  //   2. two-dimension
  //   3. three-dimension
  // CD: Failure criterion
  //   01. Maximum stress
  //   02. Maximum strain
  //   03. Hashin
  //   04. Rotem
  //   05. Tsai-Wu
  //   06. Tsai-Hill
  //   07. Christensen
  //   08. Hoffman
  // E: Failure mode
  //    0: General
  //    1. Tension
  //    2. Compression
  //    3. Tension and compression
  //    4. Shear
  //    5. Delamination
  
  struct LimitStateFunction_MCS {
    
    //------------------------------------------------------------------------------
    // Code: 12013
    // MAXIMUM STRESS - Longitudinal direction
    //
    
    virtual PetscErrorCode gfun_ply_MS_LD(ublas::vector<double> x,
                                          vector<string> vars_name,
                                          ublas::vector<double> PlyStrength,
                                          ublas::matrix<double> StressGP,
                                          double &val_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      S_23 = PlyStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      if (StressGP(0,0)>0) {                         // Fibre failure in tension
        // Evaluate the limite state function
        val_lsf = X_T - StressGP(0,0);
      }
      else {                                     // Fibre failure in compression
        // Evaluate the limit state function
        val_lsf = X_C + StressGP(0,0);
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 22013
    // MAXIMUM STRESS - Transversal direction
    //
    
    virtual PetscErrorCode gfun_ply_MS_TD(ublas::vector<double> x,//ublas::matrix<double> PlyStrength,
                                          vector<string> vars_name,
                                          ublas::vector<double> PlyStrength,
                                          ublas::matrix<double> StressGP,
                                          double &val_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      ublas::matrix<double> dsdx(3,3);
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      if (StressGP(1,1)>0) {
        // Evaluate limit state function
        val_lsf = Y_T - StressGP(1,1);
      }
      else {
        // Evaluate limit state function
        val_lsf = Y_C + StressGP(1,1);
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 22014
    // MAXIMUM STRESS - Shear failure
    //
    
    virtual PetscErrorCode gfun_ply_MS_Shear(ublas::vector<double> x,
                                             vector<string> vars_name,
                                             ublas::vector<double> PlyStrength,
                                             ublas::matrix<double> StressGP,
                                             double &val_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      if (StressGP(0,1)>0) {
        // Evaluate limit state function
        val_lsf = S_12 - StressGP(0,1);
      }
      else {
        // Evaluate limit state function
        val_lsf = S_12 + StressGP(0,1);
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 13033
    // Hashin failure theory - Fibre Mode
    //
    virtual PetscErrorCode gfun_ply_HF(ublas::vector<double> x,
                                       vector<string> vars_name,
                                       ublas::vector<double> PlyStrength,
                                       ublas::matrix<double> StressGP,
                                       double &val_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      S_23 = PlyStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      if (StressGP(0,0)>0) { // Tensile fibre mode
        // Evaluate the limit state function
        val_lsf = 1 - (1/(X_T*X_T)*pow(StressGP(0,0),2)
                       + 1/(S_12*S_12)*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2)));
      }
      else { // Compressive fibre mode
        // Evaluate the limit state function
        val_lsf = 1 - (StressGP(0,0)*StressGP(0,0))/(X_C*X_C);
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 23033
    // Hashin failure theory - Matrix Mode
    //
    virtual PetscErrorCode gfun_ply_HM(ublas::vector<double> x,
                                       vector<string> vars_name,
                                       ublas::vector<double> PlyStrength,
                                       ublas::matrix<double> StressGP,
                                       double &val_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      S_23 = PlyStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      // Calculate the transverse shear strength using Christensen's formula
      double eta;
      eta = 0.2519;
      S_23 = eta*Y_T*Y_C;
      
      //cout<<"\n\nThe value is "<<(StressGP(1,1) + StressGP(2,2))<<endl;
      
      if ((StressGP(1,1) + StressGP(2,2))>0) { // Tensile matrix mode
        // Evaluate the limit state function
        val_lsf = 1 - (1/(Y_T*Y_T)*pow((StressGP(1,1)+StressGP(2,2)),2)
                       + 1/(S_23*S_23)*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2))
                       + 1/(S_12*S_12)*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2)));
      }
      else { // Compressive matrix mode
        // Evaluate the limit state function
        val_lsf = 1 - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(StressGP(1,1) + StressGP(2,2))
                       + 1/(4*S_23*S_23)*pow((StressGP(1,1) + StressGP(2,2)),2)
                       + 1/(S_23*S_23)*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2))
                       + 1/(S_12*S_12)*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2)));
      }
      
      PetscFunctionReturn(0);
    }
    
    
    //------------------------------------------------------------------------------
    // TSAI-WU
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Wu_0(ublas::vector<double> x,//ublas::matrix<double> PlyStrength,
                                            vector<string> vars_name,
                                            ublas::vector<double> PlyStrength,
                                            ublas::matrix<double> StressGP,
                                            double &val_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      double F1, F2, F11, F22, F12, F44, F66;
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      S_23 = PlyStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
    
    //--------------------------------------------------------------------------
    // Code: 42050
    // TSAI-WU failure theory
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Wu_2D(ublas::vector<double> x,
                                               vector<string> vars_name,
                                               ublas::vector<double> PlyStrength,
                                               ublas::matrix<double> StressGP,
                                               double &val_lsf) {
      
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      double F1, F2, F11, F22, F12, F66;
      F1 = (1/X_T - 1/X_C);
      F2 = (1/Y_T - 1/Y_C);
      F11 = 1/(X_T*X_C);
      F22 = 1/(Y_T*Y_C);
      F12 = -1/(2*sqrt(X_T*X_C*Y_T*Y_C)); // Empirical suggestion: Mises-Hencky criterion
      F66 = 1/(S_12*S_12);
      
      //StressGP = StressGP;
      // Evaluate the limit state function
      val_lsf = 1- (  F1*StressGP(0,0) + F2*StressGP(1,1)
                    + F11*StressGP(0,0)*StressGP(0,0)
                    + F22*StressGP(1,1)*StressGP(1,1)
                    + F66*StressGP(0,1)*StressGP(0,1)
                    + 2*F12*StressGP(0,0)*StressGP(1,1));
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 43050
    // TSAI-WU
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Wu(ublas::vector<double> x,
                                            vector<string> vars_name,
                                            ublas::vector<double> MatStrength,
                                            ublas::matrix<double> StressGP,
                                            double &val_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      double F1, F2, F11, F22, F44, F66, F12, F23;
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = MatStrength(0);
      X_C  = MatStrength(1);
      Y_T  = MatStrength(2);
      Y_C  = MatStrength(3);//cout<<"The Y_T is: "<<Y_T<<endl;
      S_12 = MatStrength(4);
      S_23 = MatStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      // Calculate the transverse shear strength using Christensen's formula
      double eta;
      eta = 0.2519; // CFRP
      //eta = 0.2910; // GFRP
      S_23 = eta*Y_T*Y_C;
      
      F1 = (1/X_T - 1/X_C);
      F2 = (1/Y_T - 1/Y_C);
      F11 = 1/(X_T*X_C);
      F22 = 1/(Y_T*Y_C);
      F12 = -1/(2*sqrt(X_T*X_C*Y_T*Y_C)); // Empirical suggestion: Mises-Hencky criterion
      F23 = -1/(2*Y_T*Y_C);
      F44 = 1/(S_23*S_23);
      F66 = 1/(S_12*S_12);
      
      val_lsf = 1- (  F1*StressGP(0,0)
                    + F2*(StressGP(1,1) + StressGP(2,2))
                    + F11*StressGP(0,0)*StressGP(0,0)
                    + F22*(StressGP(1,1)*StressGP(1,1) + StressGP(2,2)*StressGP(2,2))
                    + F44*StressGP(1,2)*StressGP(1,2)
                    + F66*(StressGP(0,1)*StressGP(0,1) + StressGP(2,0)*StressGP(2,0))
                    + 2*F12*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2))
                    + 2*F23*StressGP(1,1)*StressGP(2,2));
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 44050
    // TSAI-WU
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Wu_Christensen(ublas::vector<double> x,
                                                        vector<string> vars_name,
                                                        ublas::vector<double> MatStrength,
                                                        ublas::matrix<double> StressGP,
                                                        double &val_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      double F1, F2, F11, F22, F44, F66, F12;
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = MatStrength(0);
      X_C  = MatStrength(1);
      Y_T  = MatStrength(2);
      Y_C  = MatStrength(3);//cout<<"The Y_T is: "<<Y_T<<endl;
      S_12 = MatStrength(4);
      S_23 = MatStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      // Calculate the transverse shear strength using Christensen's formula
      double eta;
      eta = 0.2519;
      S_23 = eta*Y_T*Y_C;
      
      F1 = (1/X_T - 1/X_C);
      F2 = (1/Y_T - 1/Y_C);
      F11 = 1/(X_T*X_C);
      F22 = 1/(Y_T*Y_C);
      F12 = -1/(2*sqrt(X_T*X_C*Y_T*Y_C)); // Empirical suggestion: Mises-Hencky criterion
      F44 = 1/(S_23*S_23);
      F66 = 1/(S_12*S_12);
      
      val_lsf = 1- (  F1*StressGP(0,0)
                    + F2*(StressGP(1,1) + StressGP(2,2))
                    + F11*StressGP(0,0)*StressGP(0,0)
                    + F22*pow((StressGP(1,1) + StressGP(2,2)),2)
                    + F44*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2))
                    + F66*(StressGP(0,1)*StressGP(0,1) + StressGP(2,0)*StressGP(2,0))
                    + 2*F12*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
      
      PetscFunctionReturn(0);
    }
    
    
    //------------------------------------------------------------------------------
    // Code: 42060
    // Tsai-Hill failure theory
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Hill_2D(ublas::vector<double> x,
                                                 vector<string> vars_name,
                                                 ublas::vector<double> PlyStrength,
                                                 ublas::matrix<double> StressGP,
                                                 double &val_lsf) {
      
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double X,Y;
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
    
    //------------------------------------------------------------------------------
    // Code: 43060
    // Tsai-Hill failure theory
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Hill(ublas::vector<double> x,
                                              vector<string> vars_name,
                                              ublas::vector<double> PlyStrength,
                                              ublas::matrix<double> StressGP,
                                              double &val_lsf) {
      
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      double X,Y;
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      S_23 = PlyStrength(5);
      
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      // Calculate the transverse shear strength using Christensen's formula
      double eta;
      eta = 0.2519;
      S_23 = eta*Y_T*Y_C;
      
      val_lsf = 1- (  StressGP(0,0)*StressGP(0,0)/(X*X)
                    + pow((StressGP(1,1) - StressGP(2,2)),2)/(Y*Y)
                    + (StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2))/(S_12*S_12)
                    + StressGP(1,2)*StressGP(1,2)/(S_23*S_23)
                    - (StressGP(0,0)*StressGP(1,1) + StressGP(0,0)*StressGP(2,2) - StressGP(1,1)*StressGP(2,2))/(X*X));
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 13073
    // RICHARD CHRISTENSEN: FIBRE CONTROLLED FAILURE
    //
    
    virtual PetscErrorCode gfun_ply_RCF(ublas::vector<double> x,
                                        vector<string> vars_name,
                                        ublas::vector<double> PlyStrength,
                                        ublas::matrix<double> StressGP,
                                        double &val_lsf) {
      PetscFunctionBegin;
      
      double T_11; // tensile strength in the fibre direction
      double T_22; // tensile strength in the transverse direction
      double C_11; // compressive strength in the fibre direction
      double C_22; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength

      T_11 = PlyStrength(0);
      C_11 = PlyStrength(1);
      T_22 = PlyStrength(2);
      C_22 = PlyStrength(3);
      S_12 = PlyStrength(4);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
        if (vars_name[i].compare(0,2,"XT") == 0) {
          T_11 = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          C_11 = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          T_22 = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          C_22 = x(i-1);
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          S_12 = x(i-1);
        }
      }
      
      
      // Evaluate the limit state function
      val_lsf = 1 - ((1/T_11 - 1/C_11)*StressGP(0,0)
                     + 1/(T_11*C_11)*StressGP(0,0)*StressGP(0,0));
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 23073
    // RICHARD CHRISTENSEN: MATRIX CONTROLLED FAILURE
    //
    
    virtual PetscErrorCode gfun_ply_RCM(ublas::vector<double> x,
                                        vector<string> vars_name,
                                        ublas::vector<double> PlyStrength,
                                        ublas::matrix<double> StressGP,
                                        double &val_lsf) {
      PetscFunctionBegin;
      
      double T_11; // tensile strength in the fibre direction
      double T_22; // tensile strength in the transverse direction
      double C_11; // compressive strength in the fibre direction
      double C_22; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength

      T_11 = PlyStrength(0);
      C_11 = PlyStrength(1);
      T_22 = PlyStrength(2);
      C_22 = PlyStrength(3);
      S_12 = PlyStrength(4);
      S_23 = PlyStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
        if (vars_name[i].compare(0,2,"XT") == 0) {
          T_11 = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          C_11 = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          T_22 = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          C_22 = x(i-1);
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          S_12 = x(i-1);
        }
        else if (vars_name[i].compare(0,3,"S23") == 0) {
          S_23 = x(i-1);
        }
      }
      
      // Calculate the transverse shear strength using Christensen's formula
      double eta;
      eta = 0.2519;
      S_23 = eta*T_22*C_22;
      
      // Evaluate the limit state function
      val_lsf = 1 - ((1/T_22 - 1/C_22)*(StressGP(1,1) + StressGP(2,2))
                     + 1/(T_22*C_22)*pow((StressGP(1,1)+ StressGP(2,2)),2)
                     + 1/(S_12*S_12)*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2))
                     + 1/(S_23*S_23)*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2)));
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 42080
    // Hoffman failure theory
    
    virtual PetscErrorCode gfun_ply_Hoffman_2D(ublas::vector<double> x,
                                               vector<string> vars_name,
                                               ublas::vector<double> PlyStrength,
                                               ublas::matrix<double> StressGP,
                                               double &val_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      val_lsf = 1 - ( StressGP(0,0)*(1/X_T - 1/X_C)
                     + StressGP(1,1)*(1/Y_T - 1/Y_C)
                     + StressGP(0,0)*StressGP(0,0)/(X_T*X_C)
                     + StressGP(1,1)*StressGP(1,1)/(Y_T*Y_C)
                     + StressGP(0,1)*StressGP(0,1)/(S_12*S_12)
                     - StressGP(0,0)*StressGP(1,1)/(X_T*X_C));
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 43080
    // Hoffman failure theory
    
    virtual PetscErrorCode gfun_ply_Hoffman(ublas::vector<double> x,
                                            vector<string> vars_name,
                                            ublas::vector<double> PlyStrength,
                                            ublas::matrix<double> StressGP,
                                            double &val_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      S_23 = PlyStrength(5);
      
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      // Calculate the transverse shear strength using Christensen's formula
      double eta;
      eta = 0.2519;
      S_23 = eta*Y_T*Y_C;
      
      // Evaluate the limit state function
      val_lsf = 1 - (  StressGP(0,0)*(1/X_T - 1/X_C)
                     + (StressGP(1,1) + StressGP(2,2))*(1/Y_T - 1/Y_C)
                     + StressGP(0,0)*StressGP(0,0)/(X_T*X_C)
                     + pow((StressGP(1,1) - StressGP(2,2)),2)/(Y_T*Y_C)
                     + (StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2))/(S_12*S_12)
                     + StressGP(1,2)*StressGP(1,2)/(S_23*S_23)
                     - (StressGP(0,0)*StressGP(1,1) + StressGP(0,0)*StressGP(2,2) - StressGP(1,1)*StressGP(2,2))/(X_T*X_C));
      
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
  };
  
  struct LSF_Composite_Lamina {
    
    //------------------------------------------------------------------------------
    // Code: 12013
    // MAXIMUM STRESS - Longitudinal direction
    //
    
    virtual PetscErrorCode gfun_ply_MS_LD(ublas::vector<double> x,//ublas::matrix<double> PlyStrength,
                                          vector<string> vars_name,
                                          ublas::vector<double> PlyStrength,
                                          ublas::matrix<double> StressGP,
                                          ublas::matrix<double> StressGP_r_Em,
                                          ublas::matrix<double> StressGP_r_NUm,
                                          ublas::matrix<double> StressGP_r_NUp,
                                          ublas::matrix<double> StressGP_r_NUpz,
                                          ublas::matrix<double> StressGP_r_Ep,
                                          ublas::matrix<double> StressGP_r_Ez,
                                          ublas::matrix<double> StressGP_r_Gzp,
                                          ublas::matrix<double> StressGP_r_Ef,
                                          ublas::matrix<double> StressGP_r_NUf,
                                          ublas::matrix<double> StressGP_r_F,
                                          ublas::matrix<double> StressGP_r_Theta,
                                          ublas::matrix<double> StressGP_r_Theta_1,
                                          ublas::matrix<double> StressGP_r_Theta_2,
                                          ublas::matrix<double> StressGP_r_Theta_3,
                                          ublas::matrix<double> StressGP_r_Theta_4,
                                          double &val_lsf,
                                          ublas::vector<double> &grad_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      ublas::matrix<double> dsdx(3,3);
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      S_23 = PlyStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      if (StressGP(0,0)>0) {                         // Fibre failure in tension
        // Evaluate the limite state function
        val_lsf = X_T - StressGP(0,0);
        
        // Evaluate the partial derovative of the LSF w.r.t. basic variables
        for (unsigned i=1; i<=x.size(); i++) {
          cout<<"The variable name is "<<vars_name[i]<<endl;
          if (vars_name[i].compare(0,2,"Em") == 0) {
            // w.r.t. Em
            dsdx.clear(); dsdx = StressGP_r_Em;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,3,"NUm") == 0) {
            // w.r.t. NUm
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,3,"NUp") == 0) {
            // w.r.t. NUp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,3,"NUz") == 0) {
            // w.r.t. NUpz
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,2,"Ep") == 0) {
            // w.r.t. Ep
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,2,"Ez") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,3,"Gzp") == 0) {
            // w.r.t. Gzp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,2,"Ef") == 0) {
            // w.r.t. Ef - fibre with isotropic material
            dsdx.clear(); dsdx = StressGP_r_Ef;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,3,"NUf") == 0) {
            // w.r.t. NUf - fibre with isotropic material
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,5,"force") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,11,"orientation") == 0) {
            // w.r.t. ply angle
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,6,"theta1") == 0) {
            // w.r.t. ply angle
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,6,"theta2") == 0) {
            // w.r.t. ply angle
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,6,"theta3") == 0) {
            // w.r.t. ply angle
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,6,"theta4") == 0) {
            // w.r.t. ply angle
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
            grad_lsf(i-1) = - dsdx(0,0);
          }
          else if (vars_name[i].compare(0,2,"XT") == 0) {
            // w.r.t X_T
            grad_lsf(i-1) = 1;
          }
          else if (vars_name[i].compare(0,2,"XC") == 0) {
            // w.r.t X_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YT") == 0) {
            // w.r.t Y_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YC") == 0) {
            // w.r.t Y_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,3,"S12") == 0) {
            // w.r.t S_12
            grad_lsf(i-1) = 0;
          }
        }
      }
      else {                                     // Fibre failure in compression
        // Evaluate the limit state function
        val_lsf = X_C + StressGP(0,0);
        // Evaluate the partial derovative of the LSF w.r.t. basic variables
        for (unsigned i=1; i<=x.size(); i++) {
          cout<<"The variable name is "<<vars_name[i]<<endl;
          if (vars_name[i].compare(0,2,"Em") == 0) {
            // w.r.t. Em
            dsdx.clear(); dsdx = StressGP_r_Em;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,3,"NUm") == 0) {
            // w.r.t. NUm
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,3,"NUp") == 0) {
            // w.r.t. NUp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,3,"NUz") == 0) {
            // w.r.t. NUpz
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,2,"Ep") == 0) {
            // w.r.t. Ep
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,2,"Ez") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,3,"Gzp") == 0) {
            // w.r.t. Gzp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,2,"Ef") == 0) {
            // w.r.t. Ef - fibre with isotropic material
            dsdx.clear(); dsdx = StressGP_r_Ef;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,3,"NUf") == 0) {
            // w.r.t. NUf - fibre with isotropic materila
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,5,"force") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,11,"orientation") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,6,"theta1") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,6,"theta2") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,6,"theta3") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,6,"theta4") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
            grad_lsf(i-1) = dsdx(0,0);
          }
          else if (vars_name[i].compare(0,2,"XT") == 0) {
            // w.r.t X_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"XC") == 0) {
            // w.r.t X_C
            grad_lsf(i-1) = 1;
          }
          else if (vars_name[i].compare(0,2,"YT") == 0) {
            // w.r.t Y_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YC") == 0) {
            // w.r.t Y_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,3,"S12") == 0) {
            // w.r.t S_12
            grad_lsf(i-1) = 0;
          }
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 22013
    // MAXIMUM STRESS - Transversal direction
    //
    
    virtual PetscErrorCode gfun_ply_MS_TD(ublas::vector<double> x,//ublas::matrix<double> PlyStrength,
                                          vector<string> vars_name,
                                          ublas::vector<double> PlyStrength,
                                          ublas::matrix<double> StressGP,
                                          ublas::matrix<double> StressGP_r_Em,
                                          ublas::matrix<double> StressGP_r_NUm,
                                          ublas::matrix<double> StressGP_r_NUp,
                                          ublas::matrix<double> StressGP_r_NUpz,
                                          ublas::matrix<double> StressGP_r_Ep,
                                          ublas::matrix<double> StressGP_r_Ez,
                                          ublas::matrix<double> StressGP_r_Gzp,
                                          ublas::matrix<double> StressGP_r_Ef,
                                          ublas::matrix<double> StressGP_r_NUf,
                                          ublas::matrix<double> StressGP_r_F,
                                          ublas::matrix<double> StressGP_r_Theta,
                                          ublas::matrix<double> StressGP_r_Theta_1,
                                          ublas::matrix<double> StressGP_r_Theta_2,
                                          ublas::matrix<double> StressGP_r_Theta_3,
                                          ublas::matrix<double> StressGP_r_Theta_4,
                                          double &val_lsf,
                                          ublas::vector<double> &grad_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      ublas::matrix<double> dsdx(3,3);
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      if (StressGP(1,1)>0) {
        // Evaluate limit state function
        val_lsf = Y_T - StressGP(1,1);
        // Evaluate the partial derovative of the LSF w.r.t. basic variables
        for (unsigned i=1; i<=x.size(); i++) {
          cout<<"The variable name is "<<vars_name[i]<<endl;
          if (vars_name[i].compare(0,2,"Em") == 0) {
            // w.r.t. Em
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Em;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,3,"NUm") == 0) {
            // w.r.t. NUm
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,3,"NUp") == 0) {
            // w.r.t. NUp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,3,"NUz") == 0) {
            // w.r.t. NUpz
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,2,"Ep") == 0) {
            // w.r.t. Ep
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,2,"Ez") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,3,"Gzp") == 0) {
            // w.r.t. Gzp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,2,"Ef") == 0) {
            // w.r.t. Ef - fibre with isotropic material
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ef;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,3,"NUf") == 0) {
            // w.r.t. NUf - fibre with isotropic material
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,5,"force") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,11,"orientation") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,6,"theta1") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,6,"theta2") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,6,"theta3") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,6,"theta4") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
            grad_lsf(i-1) = - dsdx(1,1);
          }
          else if (vars_name[i].compare(0,2,"XT") == 0) {
            // w.r.t X_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"XC") == 0) {
            // w.r.t X_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YT") == 0) {
            // w.r.t Y_T
            grad_lsf(i-1) = 1;
          }
          else if (vars_name[i].compare(0,2,"YC") == 0) {
            // w.r.t Y_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,3,"S12") == 0) {
            // w.r.t S_12
            grad_lsf(i-1) = 0;
          }
        }
      }
      else {
        // Evaluate limit state function
        val_lsf = Y_C + StressGP(1,1);
        // Evaluate the partial derovative of the LSF w.r.t. basic variables
        for (unsigned i=1; i<=x.size(); i++) {
          cout<<"The variable name is "<<vars_name[i]<<endl;
          if (vars_name[i].compare(0,2,"Em") == 0) {
            // w.r.t. Em
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Em;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,3,"NUm") == 0) {
            // w.r.t. NUm
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,3,"NUp") == 0) {
            // w.r.t. NUp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,3,"NUz") == 0) {
            // w.r.t. NUpz
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,2,"Ep") == 0) {
            // w.r.t. Ep
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,2,"Ez") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,3,"Gzp") == 0) {
            // w.r.t. Gzp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,2,"Ef") == 0) {
            // w.r.t. Ef - fibre with isotropic material
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ef;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,3,"NUf") == 0) {
            // w.r.t. NUf - fibre with isotropic material
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,5,"force") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,11,"orientation") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,6,"theta1") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,6,"theta2") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,6,"theta3") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,6,"theta4") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
            grad_lsf(i-1) = dsdx(1,1);
          }
          else if (vars_name[i].compare(0,2,"XT") == 0) {
            // w.r.t X_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"XC") == 0) {
            // w.r.t X_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YT") == 0) {
            // w.r.t Y_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YC") == 0) {
            // w.r.t Y_C
            grad_lsf(i-1) = 1;
          }
          else if (vars_name[i].compare(0,3,"S12") == 0) {
            // w.r.t S_12
            grad_lsf(i-1) = 0;
          }
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 22014
    // MAXIMUM STRESS - Shear failure
    //
    
    virtual PetscErrorCode gfun_ply_MS_Shear(ublas::vector<double> x,//ublas::matrix<double> PlyStrength,
                                             vector<string> vars_name,
                                             ublas::vector<double> PlyStrength,
                                             ublas::matrix<double> StressGP,
                                             ublas::matrix<double> StressGP_r_Em,
                                             ublas::matrix<double> StressGP_r_NUm,
                                             ublas::matrix<double> StressGP_r_NUp,
                                             ublas::matrix<double> StressGP_r_NUpz,
                                             ublas::matrix<double> StressGP_r_Ep,
                                             ublas::matrix<double> StressGP_r_Ez,
                                             ublas::matrix<double> StressGP_r_Gzp,
                                             ublas::matrix<double> StressGP_r_Ef,
                                             ublas::matrix<double> StressGP_r_NUf,
                                             ublas::matrix<double> StressGP_r_F,
                                             ublas::matrix<double> StressGP_r_Theta,
                                             ublas::matrix<double> StressGP_r_Theta_1,
                                             ublas::matrix<double> StressGP_r_Theta_2,
                                             ublas::matrix<double> StressGP_r_Theta_3,
                                             ublas::matrix<double> StressGP_r_Theta_4,
                                             double &val_lsf,
                                             ublas::vector<double> &grad_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      ublas::matrix<double> dsdx(3,3);
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      if (StressGP(0,1)>0) {
        // Evaluate limit state function
        val_lsf = S_12 - StressGP(0,1);
        // Evaluate the partial derovative of the LSF w.r.t. basic variables
        for (unsigned i=1; i<=x.size(); i++) {
          cout<<"The variable name is "<<vars_name[i]<<endl;
          if (vars_name[i].compare(0,2,"Em") == 0) {
            // w.r.t. Em
            dsdx.clear(); dsdx = StressGP_r_Em;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,3,"NUm") == 0) {
            // w.r.t. NUm
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,3,"NUp") == 0) {
            // w.r.t. NUp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,3,"NUz") == 0) {
            // w.r.t. NUpz
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,2,"Ep") == 0) {
            // w.r.t. Ep
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,2,"Ez") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,3,"Gzp") == 0) {
            // w.r.t. Gzp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,2,"Ef") == 0) {
            // w.r.t. Ef - fibre with isotropic material
            dsdx.clear(); dsdx = StressGP_r_Ef;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,3,"NUf") == 0) {
            // w.r.t. NUf - fibre with isotropic material
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,5,"force") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,11,"orientation") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,6,"theta1") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,6,"theta2") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,6,"theta3") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,6,"theta4") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
            grad_lsf(i-1) = - dsdx(0,1);
          }
          else if (vars_name[i].compare(0,2,"XT") == 0) {
            // w.r.t X_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"XC") == 0) {
            // w.r.t X_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YT") == 0) {
            // w.r.t Y_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YC") == 0) {
            // w.r.t Y_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,3,"S12") == 0) {
            // w.r.t S_12
            grad_lsf(i-1) = 1;
          }
        }
      }
      else {
        // Evaluate limit state function
        val_lsf = S_12 + StressGP(0,1);
        // Evaluate the partial derovative of the LSF w.r.t. basic variables
        for (unsigned i=1; i<=x.size(); i++) {
          cout<<"The variable name is "<<vars_name[i]<<endl;
          if (vars_name[i].compare(0,2,"Em") == 0) {
            // w.r.t. Em
            dsdx.clear(); dsdx = StressGP_r_Em;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,3,"NUm") == 0) {
            // w.r.t. NUm
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,3,"NUp") == 0) {
            // w.r.t. NUp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,3,"NUz") == 0) {
            // w.r.t. NUpz
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,2,"Ep") == 0) {
            // w.r.t. Ep
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,2,"Ez") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,3,"Gzp") == 0) {
            // w.r.t. Gzp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,2,"Ef") == 0) {
            // w.r.t. Ef - fibre with isotropic material
            dsdx.clear(); dsdx = StressGP_r_Ef;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,3,"NUf") == 0) {
            // w.r.t. NUf - fibre with isotropic mateirla
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,5,"force") == 0) {
            // w.r.t. F
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,11,"orientation") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,6,"theta1") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,6,"theta2") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,6,"theta3") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,6,"theta4") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
            grad_lsf(i-1) = dsdx(0,1);
          }
          else if (vars_name[i].compare(0,2,"XT") == 0) {
            // w.r.t X_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"XC") == 0) {
            // w.r.t X_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YT") == 0) {
            // w.r.t Y_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YC") == 0) {
            // w.r.t Y_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,3,"S12") == 0) {
            // w.r.t S_12
            grad_lsf(i-1) = 1;
          }
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 13033
    // Hashin failure theory - Fibre Mode
    //
    virtual PetscErrorCode gfun_ply_HF(ublas::vector<double> x,
                                       vector<string> vars_name,
                                       ublas::vector<double> PlyStrength,
                                       ublas::matrix<double> StressGP,
                                       ublas::matrix<double> StressGP_r_Em,
                                       ublas::matrix<double> StressGP_r_NUm,
                                       ublas::matrix<double> StressGP_r_NUp,
                                       ublas::matrix<double> StressGP_r_NUpz,
                                       ublas::matrix<double> StressGP_r_Ep,
                                       ublas::matrix<double> StressGP_r_Ez,
                                       ublas::matrix<double> StressGP_r_Gzp,
                                       ublas::matrix<double> StressGP_r_Ef,
                                       ublas::matrix<double> StressGP_r_NUf,
                                       ublas::matrix<double> StressGP_r_F,
                                       ublas::matrix<double> StressGP_r_Theta,
                                       ublas::matrix<double> StressGP_r_Theta_1,
                                       ublas::matrix<double> StressGP_r_Theta_2,
                                       ublas::matrix<double> StressGP_r_Theta_3,
                                       ublas::matrix<double> StressGP_r_Theta_4,
                                       double &val_lsf,
                                       ublas::vector<double> &grad_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      ublas::matrix<double> dsdx(3,3);
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      S_23 = PlyStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      cout<<"\n\nThe value is "<<(StressGP(1,1) + StressGP(2,2))<<endl;
      
      if (StressGP(0,0)>0) { // Tensile fibre mode
        // Evaluate the limit state function
        val_lsf = 1 - (1/(X_T*X_T)*pow(StressGP(0,0),2)
                       + 1/(S_12*S_12)*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2)));
        // Evaluate the partial derovative of the LSF w.r.t. basic variables
        for (unsigned i=1; i<=x.size(); i++) {
          cout<<"The variable name is "<<vars_name[i]<<endl;
          if (vars_name[i].compare(0,2,"Em") == 0) {
            // w.r.t. Em
            dsdx.clear(); dsdx = StressGP_r_Em;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"NUm") == 0) {
            // w.r.t. NUm
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"NUp") == 0) {
            // w.r.t. NUp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"NUz") == 0) {
            // w.r.t. NUpz
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,2,"Ep") == 0) {
            // w.r.t. Ep
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,2,"Ez") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"Gzp") == 0) {
            // w.r.t. Gzp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,2,"Ef") == 0) {
            // w.r.t. Ef - fibre with isotropic material
            dsdx.clear(); dsdx = StressGP_r_Ef;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"NUf") == 0) {
            // w.r.t. NUf - fibre with isotropic material
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,5,"force") == 0) {
            // w.r.t. F
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,11,"orientation") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,6,"theta1") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,6,"theta2") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,6,"theta3") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,6,"theta4") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
            grad_lsf(i-1) = - (2/(X_T*X_T)*StressGP(0,0)*dsdx(0,0)
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,2,"XT") == 0) {
            // w.r.t X_T
            grad_lsf(i-1) = 2/(X_T*X_T*X_T)*pow(StressGP(0,0),2);
          }
          else if (vars_name[i].compare(0,2,"XC") == 0) {
            // w.r.t X_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YT") == 0) {
            // w.r.t Y_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YC") == 0) {
            // w.r.t Y_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,3,"S12") == 0) {
            // w.r.t S_12
            grad_lsf(i-1) = 2/(S_12*S_12*S_12)*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2));
          }
        }
      }
      else { // Compressive fibre mode
        // Evaluate the limit state function
        val_lsf = 1 - (StressGP(0,0)*StressGP(0,0))/(X_C*X_C);
        // Evaluate the partial derovative of the LSF w.r.t. basic variables
        for (unsigned i=1; i<=x.size(); i++) {
          cout<<"The variable name is "<<vars_name[i]<<endl;
          if (vars_name[i].compare(0,2,"Em") == 0) {
            // w.r.t. Em
            dsdx.clear(); dsdx = StressGP_r_Em;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,3,"NUm") == 0) {
            // w.r.t. NUm
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,3,"NUp") == 0) {
            // w.r.t. NUp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,3,"NUz") == 0) {
            // w.r.t. NUpz
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,2,"Ep") == 0) {
            // w.r.t. Ep
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,2,"Ez") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,3,"Gzp") == 0) {
            // w.r.t. Gzp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,2,"Ef") == 0) {
            // w.r.t. Ef - fibre with isotropic material
            dsdx.clear(); dsdx = StressGP_r_Ef;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,3,"NUf") == 0) {
            // w.r.t. NUf - fibre with isotropci material
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,5,"force") == 0) {
            // w.r.t. F
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,11,"orientation") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,6,"theta1") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,6,"theta2") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,6,"theta3") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,6,"theta4") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
            grad_lsf(i-1) = - 2*StressGP(0,0)*dsdx(0,0)/(X_C*X_C);
          }
          else if (vars_name[i].compare(0,2,"XT") == 0) {
            // w.r.t X_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"XC") == 0) {
            // w.r.t X_C
            grad_lsf(i-1) = 2*(StressGP(0,0)*StressGP(0,0))/(X_C*X_C*X_C);
          }
          else if (vars_name[i].compare(0,2,"YT") == 0) {
            // w.r.t Y_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YC") == 0) {
            // w.r.t Y_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,3,"S12") == 0) {
            // w.r.t S_12
            grad_lsf(i-1) = 0;
          }
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 23033
    // Hashin failure theory - Matrix Mode
    //
    virtual PetscErrorCode gfun_ply_HM(ublas::vector<double> x,
                                       vector<string> vars_name,
                                       ublas::vector<double> PlyStrength,
                                       ublas::matrix<double> StressGP,
                                       ublas::matrix<double> StressGP_r_Em,
                                       ublas::matrix<double> StressGP_r_NUm,
                                       ublas::matrix<double> StressGP_r_NUp,
                                       ublas::matrix<double> StressGP_r_NUpz,
                                       ublas::matrix<double> StressGP_r_Ep,
                                       ublas::matrix<double> StressGP_r_Ez,
                                       ublas::matrix<double> StressGP_r_Gzp,
                                       ublas::matrix<double> StressGP_r_Ef,
                                       ublas::matrix<double> StressGP_r_NUf,
                                       ublas::matrix<double> StressGP_r_F,
                                       ublas::matrix<double> StressGP_r_Theta,
                                       ublas::matrix<double> StressGP_r_Theta_1,
                                       ublas::matrix<double> StressGP_r_Theta_2,
                                       ublas::matrix<double> StressGP_r_Theta_3,
                                       ublas::matrix<double> StressGP_r_Theta_4,
                                       double &val_lsf,
                                       ublas::vector<double> &grad_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      ublas::matrix<double> dsdx(3,3);
      
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      S_23 = PlyStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      // Calculate the transverse shear strength using Christensen's formula
      double eta;
      eta = 0.2519;
      S_23 = eta*Y_T*Y_C;
      
      cout<<"\n\nThe value is "<<(StressGP(1,1) + StressGP(2,2))<<endl;
      
      if ((StressGP(1,1) + StressGP(2,2))>0) { // Tensile matrix mode
        // Evaluate the limit state function
        val_lsf = 1 - (1/(Y_T*Y_T)*pow((StressGP(1,1)+StressGP(2,2)),2)
                       + 1/(S_23*S_23)*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2))
                       + 1/(S_12*S_12)*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2)));
        // Evaluate the partial derovative of the LSF w.r.t. basic variables
        for (unsigned i=1; i<=x.size(); i++) {
          cout<<"The variable name is "<<vars_name[i]<<endl;
          if (vars_name[i].compare(0,2,"Em") == 0) {
            // w.r.t. Em
            dsdx.clear(); dsdx = StressGP_r_Em;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"NUm") == 0) {
            // w.r.t. NUm
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"NUp") == 0) {
            // w.r.t. NUp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"NUz") == 0) {
            // w.r.t. NUpz
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,2,"Ep") == 0) {
            // w.r.t. Ep
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,2,"Ez") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"Gzp") == 0) {
            // w.r.t. Gzp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          if (vars_name[i].compare(0,2,"Ef") == 0) {
            // w.r.t. Ef - fibre with isotropic material
            dsdx.clear(); dsdx = StressGP_r_Ef;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"NUf") == 0) {
            // w.r.t. NUf - fibre with isotropic material
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,5,"force") == 0) {
            // w.r.t. F
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,11,"orientation") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,6,"theta1") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,6,"theta2") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,6,"theta3") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,6,"theta4") == 0) {
            // w.r.t. Theta
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
            grad_lsf(i-1) = - (2/(Y_T*Y_T)*(StressGP(1,1)+StressGP(2,2))*(dsdx(1,1)+dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,2,"XT") == 0) {
            // w.r.t X_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"XC") == 0) {
            // w.r.t X_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YT") == 0) {
            // w.r.t Y_T
            // grad_lsf(i-1) = 2/(Y_T*Y_T*Y_T)*pow((StressGP(1,1)+StressGP(2,2)),2);
            grad_lsf(i-1) = 2/(Y_T*Y_T*Y_T)*pow((StressGP(1,1)+StressGP(2,2)),2)
                           + 2*eta*Y_C/(S_23*S_23*S_23)*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2));
          }
          else if (vars_name[i].compare(0,2,"YC") == 0) {
            // w.r.t Y_C
            //grad_lsf(i-1) = 0;
            grad_lsf(i-1) = 2*eta*Y_T/(S_23*S_23*S_23)*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2));
          }
          else if (vars_name[i].compare(0,3,"S12") == 0) {
            // w.r.t S_12
            grad_lsf(i-1) = 2/(S_12*S_12*S_12)*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2));
          }
        }
      }
      else { // Compressive matrix mode
        // Evaluate the limit state function
        val_lsf = 1 - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(StressGP(1,1) + StressGP(2,2))
                       + 1/(4*S_23*S_23)*pow((StressGP(1,1) + StressGP(2,2)),2)
                       + 1/(S_23*S_23)*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2))
                       + 1/(S_12*S_12)*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2)));
        // Evaluate the partial derovative of the LSF w.r.t. basic variables
        for (unsigned i=1; i<=x.size(); i++) {
          cout<<"The variable name is "<<vars_name[i]<<endl;
          if (vars_name[i].compare(0,2,"Em") == 0) {
            // w.r.t. Em
            dsdx.clear(); dsdx = StressGP_r_Em;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"NUm") == 0) {
            // w.r.t. NUm
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"NUp") == 0) {
            // w.r.t. NUp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"NUz") == 0) {
            // w.r.t. NUpz
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,2,"Ep") == 0) {
            // w.r.t. Ep
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,2,"Ez") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"Gzp") == 0) {
            // w.r.t. Gzp
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,2,"Ef") == 0) {
            // w.r.t. Ef - fibre with isotropic material
            dsdx.clear(); dsdx = StressGP_r_Ef;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,3,"NUf") == 0) {
            // w.r.t. NUf - fibre with isotropic material
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,5,"force") == 0) {
            // w.r.t. F
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,11,"orientation") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,6,"theta1") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,6,"theta2") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,6,"theta3") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,6,"theta4") == 0) {
            // w.r.t. Ez
            dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
            grad_lsf(i-1) = - (1/Y_C*(pow(Y_C/2/S_23,2)-1)*(dsdx(1,1) + dsdx(2,2))
                               + 1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                               + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                               + 1/(S_12*S_12)*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(0,2)*dsdx(0,2)));
          }
          else if (vars_name[i].compare(0,2,"XT") == 0) {
            // w.r.t X_T
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"XC") == 0) {
            // w.r.t X_C
            grad_lsf(i-1) = 0;
          }
          else if (vars_name[i].compare(0,2,"YT") == 0) {
            // w.r.t Y_T
            //grad_lsf(i-1) = 0;
            grad_lsf(i-1) = - (- eta*Y_C*Y_C/(2*pow(S_23,3))*(StressGP(1,1) + StressGP(2,2))
                               - eta*Y_C/(2*S_23*S_23*S_23)*pow((StressGP(1,1) + StressGP(2,2)),2)
                               - 2*eta*Y_C/(S_23*S_23*S_23)*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2)));
          }
          else if (vars_name[i].compare(0,2,"YC") == 0) {
            // w.r.t Y_C
            //grad_lsf(i-1) = - (-2/(Y_C*Y_C)*(pow(Y_C/2/S_23,2)-1)*(StressGP(1,1) + StressGP(2,2))
            //                   +1/(2*S_23*S_23)*(StressGP(1,1) + StressGP(2,2)));
            grad_lsf(i-1) = - (-1/(Y_C*Y_C)*(pow(1/(2*eta*Y_T),2)-1)*(StressGP(1,1) + StressGP(2,2))
                               - eta*Y_T/(2*S_23*S_23*S_23)*pow((StressGP(1,1) + StressGP(2,2)),2)
                               - 2*eta*Y_T/(S_23*S_23*S_23)*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2)));
                   
          }
          else if (vars_name[i].compare(0,3,"S12") == 0) {
            // w.r.t S_12
            grad_lsf(i-1) = 2/(S_12*S_12*S_12)*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2));
          }
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    //--------------------------------------------------------------------------
    // Code: 42050
    // TSAI-WU failure theory
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Wu_2D(ublas::vector<double> x,
                                               vector<string> vars_name,
                                               ublas::vector<double> PlyStrength,
                                               ublas::matrix<double> StressGP,
                                               ublas::matrix<double> StressGP_r_Em,
                                               ublas::matrix<double> StressGP_r_NUm,
                                               ublas::matrix<double> StressGP_r_NUp,
                                               ublas::matrix<double> StressGP_r_NUpz,
                                               ublas::matrix<double> StressGP_r_Ep,
                                               ublas::matrix<double> StressGP_r_Ez,
                                               ublas::matrix<double> StressGP_r_Gzp,
                                               ublas::matrix<double> StressGP_r_Ef,
                                               ublas::matrix<double> StressGP_r_NUf,
                                               ublas::matrix<double> StressGP_r_F,
                                               ublas::matrix<double> StressGP_r_Theta,
                                               ublas::matrix<double> StressGP_r_Theta_1,
                                               ublas::matrix<double> StressGP_r_Theta_2,
                                               ublas::matrix<double> StressGP_r_Theta_3,
                                               ublas::matrix<double> StressGP_r_Theta_4,
                                               double &val_lsf,
                                               ublas::vector<double> &grad_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      double F1, F2, F11, F22, F12, F66;
      F1 = (1/X_T - 1/X_C);
      F2 = (1/Y_T - 1/Y_C);
      F11 = 1/(X_T*X_C);
      F22 = 1/(Y_T*Y_C);
      F12 = -1/(2*sqrt(X_T*X_C*Y_T*Y_C)); // Empirical suggestion: Mises-Hencky criterion
      F66 = 1/(S_12*S_12);
      
      //StressGP = StressGP;
      // Evaluate the limit state function
      val_lsf = 1- (  F1*StressGP(0,0) + F2*StressGP(1,1)
                    + F11*StressGP(0,0)*StressGP(0,0)
                    + F22*StressGP(1,1)*StressGP(1,1)
                    + F66*StressGP(0,1)*StressGP(0,1)
                    + 2*F12*StressGP(0,0)*StressGP(1,1));
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      for (unsigned i=1; i<=x.size(); i++) {
        cout<<"The variable name is "<<vars_name[i]<<endl;
        if (vars_name[i].compare(0,2,"Em") == 0) {
          // w.r.t. Em
          dsdx.clear(); dsdx = StressGP_r_Em;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,3,"NUm") == 0) {
          // w.r.t. NUm
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,3,"NUp") == 0) {
          // w.r.t. NUp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,3,"NUz") == 0) {
          // w.r.t. NUpz
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,2,"Ep") == 0) {
          // w.r.t. Ep
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,2,"Ez") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,3,"Gzp") == 0) {
          // w.r.t. Gzp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,2,"Ef") == 0) {
          // w.r.t. Ef - fibre with isotropic material
          dsdx.clear(); dsdx = StressGP_r_Ef;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,3,"NUf") == 0) {
          // w.r.t. NUf - fibre with isotropic material
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,2,"XT") == 0) {
          // w.r.t. XT
          double F12_XT = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_C*Y_T*Y_C;
          double F1_XT = -1/(X_T*X_T);
          double F11_XT = -1/(X_T*X_T*X_C);
          grad_lsf(i-1) = - (  F1_XT*StressGP(0,0)
                             + F11_XT*StressGP(0,0)*StressGP(0,0)
                             + 2*F12_XT*StressGP(0,0)*StressGP(1,1));
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          // w.r.t. XC
          double F1_XC = 1/(X_C*X_C);
          double F11_XC = -1/(X_T*X_C*X_C);
          double F12_XC = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*Y_T*Y_C;
          grad_lsf(i-1) = - (  F1_XC*StressGP(0,0)
                             + F11_XC*StressGP(0,0)*StressGP(0,0)
                             + 2*F12_XC*StressGP(0,0)*StressGP(1,1));
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          // w.r.t. YT
          double F2_YT = -1/(Y_T*Y_T);
          double F22_YT = -1/(Y_T*Y_T*Y_C);
          double F12_YT = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*X_C*Y_C;
          grad_lsf(i-1) = - (  F2_YT*StressGP(1,1)
                             + F22_YT*StressGP(1,1)*StressGP(1,1)
                             + 2*F12_YT*StressGP(0,0)*StressGP(1,1));
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          // w.r.t. YC
          double F2_YC = 1/(Y_C*Y_C);
          double F22_YC = -1/(Y_T*Y_C*Y_C);
          double F12_YC = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*X_C*Y_T;
          grad_lsf(i-1) = - (  F2_YC*StressGP(1,1)
                             + F22_YC*StressGP(1,1)*StressGP(1,1)
                             + 2*F12_YC*StressGP(0,0)*StressGP(1,1));
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          // w.r.t. S12
          grad_lsf(i-1) = - (-2*StressGP(0,1)*StressGP(0,1)/pow(S_12,3));
        }
        else if (vars_name[i].compare(0,5,"force") == 0) {
          // w.r.t. F
          dsdx.clear(); dsdx = StressGP_r_F;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,11,"orientation") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,6,"theta1") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_1;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,6,"theta2") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_2;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,6,"theta3") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_3;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
        else if (vars_name[i].compare(0,6,"theta4") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_4;
          grad_lsf(i-1) = - (  F1*dsdx(0,0) + F2*dsdx(1,1)
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*StressGP(1,1)*dsdx(1,1)
                             + 2*F66*StressGP(0,1)*dsdx(0,1)
                             + 2*F12*dsdx(0,0)*StressGP(1,1)
                             + 2*F12*StressGP(0,0)*dsdx(1,1));
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 43050
    // TSAI-WU
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Wu(ublas::vector<double> x,
                                            vector<string> vars_name,
                                            ublas::vector<double> MatStrength,
                                            ublas::matrix<double> StressGP,
                                            ublas::matrix<double> StressGP_r_Em,
                                            ublas::matrix<double> StressGP_r_NUm,
                                            ublas::matrix<double> StressGP_r_NUp,
                                            ublas::matrix<double> StressGP_r_NUpz,
                                            ublas::matrix<double> StressGP_r_Ep,
                                            ublas::matrix<double> StressGP_r_Ez,
                                            ublas::matrix<double> StressGP_r_Gzp,
                                            ublas::matrix<double> StressGP_r_Ef,
                                            ublas::matrix<double> StressGP_r_NUf,
                                            ublas::matrix<double> StressGP_r_F,
                                            ublas::matrix<double> StressGP_r_Theta,
                                            ublas::matrix<double> StressGP_r_Theta_1,
                                            ublas::matrix<double> StressGP_r_Theta_2,
                                            ublas::matrix<double> StressGP_r_Theta_3,
                                            ublas::matrix<double> StressGP_r_Theta_4,
                                            double &val_lsf,
                                            ublas::vector<double> &grad_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      double F1, F2, F11, F22, F44, F66, F12, F23;
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = MatStrength(0);
      X_C  = MatStrength(1);
      Y_T  = MatStrength(2);
      Y_C  = MatStrength(3);//cout<<"The Y_T is: "<<Y_T<<endl;
      S_12 = MatStrength(4);
      S_23 = MatStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      // Calculate the transverse shear strength using Christensen's formula
      double eta;
      eta = 0.2519; // CFRP
      //eta = 0.2910; // GFRP
      S_23 = eta*Y_T*Y_C;
      
      F1 = (1/X_T - 1/X_C);
      F2 = (1/Y_T - 1/Y_C);
      F11 = 1/(X_T*X_C);
      F22 = 1/(Y_T*Y_C);
      F12 = -1/(2*sqrt(X_T*X_C*Y_T*Y_C)); // Empirical suggestion: Mises-Hencky criterion
      F23 = -1/(2*Y_T*Y_C);
      F44 = 1/(S_23*S_23);
      F66 = 1/(S_12*S_12);
      
      val_lsf = 1- (  F1*StressGP(0,0)
                    + F2*(StressGP(1,1) + StressGP(2,2))
                    + F11*StressGP(0,0)*StressGP(0,0)
                    + F22*(StressGP(1,1)*StressGP(1,1) + StressGP(2,2)*StressGP(2,2))
                    + F44*StressGP(1,2)*StressGP(1,2)
                    + F66*(StressGP(0,1)*StressGP(0,1) + StressGP(2,0)*StressGP(2,0))
                    + 2*F12*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2))
                    + 2*F23*StressGP(1,1)*StressGP(2,2));
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      for (unsigned i=1; i<=x.size(); i++) {
        cout<<"The variable name is "<<vars_name[i]<<endl;
        if (vars_name[i].compare(0,2,"Em") == 0) {
          // w.r.t. Em
          dsdx.clear(); dsdx = StressGP_r_Em;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
        else if (vars_name[i].compare(0,3,"NUm") == 0) {
          // w.r.t. NUm
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
        else if (vars_name[i].compare(0,3,"NUp") == 0) {
          // w.r.t. NUp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
        else if (vars_name[i].compare(0,3,"NUz") == 0) {
          // w.r.t. NUpz
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
          
        }
        else if (vars_name[i].compare(0,2,"Ep") == 0) {
          // w.r.t. Ep
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
        else if (vars_name[i].compare(0,2,"Ez") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
        else if (vars_name[i].compare(0,3,"Gzp") == 0) {
          // w.r.t. Gzp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
        else if (vars_name[i].compare(0,2,"Ef") == 0) {
          // w.r.t. Em
          dsdx.clear(); dsdx = StressGP_r_Ef;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
        else if (vars_name[i].compare(0,3,"NUf") == 0) {
          // w.r.t. NUm
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
        else if (vars_name[i].compare(0,5,"force") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
        else if (vars_name[i].compare(0,11,"orientation") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
        else if (vars_name[i].compare(0,6,"theta1") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
        else if (vars_name[i].compare(0,6,"theta2") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
        else if (vars_name[i].compare(0,6,"theta3") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
        else if (vars_name[i].compare(0,6,"theta4") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                             + 2*F44*StressGP(1,2)*dsdx(1,2)
                             + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                             + 2*F23*dsdx(1,1)*StressGP(2,2)
                             + 2*F23*StressGP(1,1)*dsdx(2,2));
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
          double F44_YT = -2*eta*Y_C/(S_23*S_23*S_23);
          double F23_YT = 1/(2*Y_C*Y_T*Y_T);
          grad_lsf(i-1) = - (  F2_YT*(StressGP(1,1) + StressGP(2,2))
                             + F22_YT*(StressGP(1,1)*StressGP(1,1) + StressGP(2,2)*StressGP(2,2))
                             + F44_YT*StressGP(1,2)*StressGP(1,2)
                             + 2*F12_YT*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F23_YT*StressGP(1,1)*StressGP(2,2));
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          // w.r.t Y_C
          double F2_YC = 1/(Y_C*Y_C);
          double F22_YC = -1/(Y_T*Y_C*Y_C);
          double F12_YC = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*X_C*Y_T;
          double F44_YC = -2*eta*Y_T/(S_23*S_23*S_23);
          double F23_YC = 1/(2*Y_C*Y_C*Y_T);
          grad_lsf(i-1) = - (  F2_YC*(StressGP(1,1) + StressGP(2,2))
                             + F22_YC*(StressGP(1,1)*StressGP(1,1) + StressGP(2,2)*StressGP(2,2))
                             + F44_YC*StressGP(1,2)*StressGP(1,2)
                             + 2*F12_YC*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F23_YC*StressGP(1,1)*StressGP(2,2));
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          // w.r.t S_12
          grad_lsf(i-1) = 2*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2))/pow(S_12,3);
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 43050
    // TSAI-WU
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Wu_New(ublas::vector<double> x,
                                                vector<string> vars_name,
                                                ublas::vector<double> MatStrength,
                                                ublas::matrix<double> StressGP,
                                                ublas::vector<ublas::matrix<double> > StressGP_r,
                                                ublas::matrix<ublas::matrix<double> > StressGP_rs,
                                                double &val_lsf,
                                                ublas::vector<double> &grad_lsf,
                                                ublas::matrix<double> &hess_lsf,
                                                int num_ply_vars,
                                                int PSFE_order) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      double F1, F2, F11, F22, F44, F66, F12, F23;
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = MatStrength(0);
      X_C  = MatStrength(1);
      Y_T  = MatStrength(2);
      Y_C  = MatStrength(3);//cout<<"The Y_T is: "<<Y_T<<endl;
      S_12 = MatStrength(4);
      S_23 = MatStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      // Calculate the transverse shear strength using Christensen's formula
      double eta;
      eta = 0.2519; // CFRP
      //eta = 0.2910; // GFRP
      S_23 = eta*Y_T*Y_C;
      
      F1 = (1/X_T - 1/X_C);
      F2 = (1/Y_T - 1/Y_C);
      F11 = 1/(X_T*X_C);
      F22 = 1/(Y_T*Y_C);
      F12 = -1/(2*sqrt(X_T*X_C*Y_T*Y_C)); // Empirical suggestion: Mises-Hencky criterion
      F23 = -1/(2*Y_T*Y_C);
      F44 = 1/(S_23*S_23);
      F66 = 1/(S_12*S_12);
      
      val_lsf = 1- (  F1*StressGP(0,0)
                    + F2*(StressGP(1,1) + StressGP(2,2))
                    + F11*StressGP(0,0)*StressGP(0,0)
                    + F22*(StressGP(1,1)*StressGP(1,1) + StressGP(2,2)*StressGP(2,2))
                    + F44*StressGP(1,2)*StressGP(1,2)
                    + F66*(StressGP(0,1)*StressGP(0,1) + StressGP(2,0)*StressGP(2,0))
                    + 2*F12*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2))
                    + 2*F23*StressGP(1,1)*StressGP(2,2));
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      
      // ============================================
      //
      // Calculate the first order derivative of LSF
      //
      // ============================================
      for (unsigned i=0; i<x.size(); i++) {
        cout<<"The variable name is "<<vars_name[i+1]<<endl;
        if (vars_name[i+1].compare(0,2,"XT") == 0) {
          // w.r.t X_T
          double F12_XT = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_C*Y_T*Y_C;
          double F1_XT = -1/(X_T*X_T);
          double F11_XT = -1/(X_T*X_T*X_C);
          grad_lsf(i) = - (  F1_XT*StressGP(0,0)
                             + F11_XT*StressGP(0,0)*StressGP(0,0)
                             + 2*F12_XT*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
          
        }
        else if (vars_name[i+1].compare(0,2,"XC") == 0) {
          // w.r.t X_C
          double F1_XC = 1/(X_C*X_C);
          double F11_XC = -1/(X_T*X_C*X_C);
          double F12_XC = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*Y_T*Y_C;
          grad_lsf(i) = - (  F1_XC*StressGP(0,0)
                             + F11_XC*StressGP(0,0)*StressGP(0,0)
                             + 2*F12_XC*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
        }
        else if (vars_name[i+1].compare(0,2,"YT") == 0) {
          // w.r.t Y_T
          double F2_YT = -1/(Y_T*Y_T);
          double F22_YT = -1/(Y_T*Y_T*Y_C);
          double F12_YT = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*X_C*Y_C;
          double F44_YT = -2*eta*Y_C/(S_23*S_23*S_23);
          double F23_YT = 1/(2*Y_C*Y_T*Y_T);
          grad_lsf(i) = - (  F2_YT*(StressGP(1,1) + StressGP(2,2))
                             + F22_YT*(StressGP(1,1)*StressGP(1,1) + StressGP(2,2)*StressGP(2,2))
                             + F44_YT*StressGP(1,2)*StressGP(1,2)
                             + 2*F12_YT*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F23_YT*StressGP(1,1)*StressGP(2,2));
        }
        else if (vars_name[i+1].compare(0,2,"YC") == 0) {
          // w.r.t Y_C
          double F2_YC = 1/(Y_C*Y_C);
          double F22_YC = -1/(Y_T*Y_C*Y_C);
          double F12_YC = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*X_C*Y_T;
          double F44_YC = -2*eta*Y_T/(S_23*S_23*S_23);
          double F23_YC = 1/(2*Y_C*Y_C*Y_T);
          grad_lsf(i) = - (  F2_YC*(StressGP(1,1) + StressGP(2,2))
                             + F22_YC*(StressGP(1,1)*StressGP(1,1) + StressGP(2,2)*StressGP(2,2))
                             + F44_YC*StressGP(1,2)*StressGP(1,2)
                             + 2*F12_YC*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F23_YC*StressGP(1,1)*StressGP(2,2));
        }
        else if (vars_name[i+1].compare(0,3,"S12") == 0) {
          // w.r.t S_12
          grad_lsf(i) = 2*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2))/pow(S_12,3);
        }
        else {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r(i);
          grad_lsf(i) = - (  F1*dsdx(0,0)
                           + F2*(dsdx(1,1) + dsdx(2,2))
                           + 2*F11*StressGP(0,0)*dsdx(0,0)
                           + 2*F22*(StressGP(1,1)*dsdx(1,1) + StressGP(2,2)*dsdx(2,2))
                           + 2*F44*StressGP(1,2)*dsdx(1,2)
                           + 2*F66*(StressGP(0,1)*dsdx(0,1) + StressGP(2,0)*dsdx(2,0))
                           + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                           + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2))
                           + 2*F23*dsdx(1,1)*StressGP(2,2)
                           + 2*F23*StressGP(1,1)*dsdx(2,2));
        }
      }
      
      // ============================================
      //
      // Calculate the second order derivative of LSF
      //
      // ============================================
      
      if (PSFE_order == 2) {
        ublas::matrix<double> dsdx_r(3,3);
        ublas::matrix<double> dsdx_s(3,3);
        ublas::matrix<double> dsdx_rs(3,3);
        //      dsdx_rs.clear();
        
        for (unsigned ivar=0; ivar<x.size(); ivar++) {
          for (unsigned jvar=0; jvar<x.size(); jvar++) {
            dsdx_r.clear(); dsdx_s.clear(); dsdx_rs.clear();
            //
            if ((ivar<num_ply_vars) && (jvar<num_ply_vars)) {
              dsdx_r = StressGP_r(ivar); dsdx_s = StressGP_r(jvar);
              dsdx_rs = StressGP_rs(ivar,jvar);
              
              hess_lsf(ivar, jvar) = - (  F1*dsdx_rs(0,0)
                                 + F2*(dsdx_rs(1,1) + dsdx_rs(2,2))
                                 + 2*F11*(dsdx_s(0,0)*dsdx_r(0,0) + StressGP(0,0)*dsdx_rs(0,0))
                                 + 2*F22*(dsdx_s(1,1)*dsdx_r(1,1) + StressGP(1,1)*dsdx_rs(1,1) + dsdx_s(2,2)*dsdx_r(2,2) + StressGP(2,2)*dsdx_rs(2,2))
                                 + 2*F44*(dsdx_s(1,2)*dsdx_r(1,2) + StressGP(1,2)*dsdx_rs(1,2))
                                 + 2*F66*(dsdx_s(0,1)*dsdx_r(0,1) + StressGP(0,1)*dsdx_rs(0,1) + dsdx_s(2,0)*dsdx_r(2,0) + StressGP(2,0)*dsdx_rs(2,0))
                                 + 2*F12*(dsdx_rs(0,0)*(StressGP(1,1)+StressGP(2,2)) + dsdx_r(0,0)*(dsdx_s(1,1)+dsdx_s(2,2)))
                                 + 2*F12*(dsdx_s(0,0)*(dsdx_r(1,1)+dsdx_r(2,2)) + StressGP(0,0)*(dsdx_rs(1,1)+dsdx_rs(2,2)))
                                 + 2*F23*(dsdx_rs(1,1)*StressGP(2,2) + dsdx_r(1,1)*dsdx_s(2,2))
                                 + 2*F23*(dsdx_s(1,1)*dsdx_r(2,2) + StressGP(1,1)*dsdx_rs(2,2)));
              hess_lsf(jvar, ivar) = hess_lsf(ivar, jvar);
            }
            else if ((ivar<num_ply_vars) && (vars_name[jvar+1].compare(0,2,"XT")==0)) {
              // XT Em
              double F12_XT = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_C*Y_T*Y_C;
              double F1_XT = -1/(X_T*X_T);
              double F11_XT = -1/(X_T*X_T*X_C);
              
              hess_lsf(ivar, jvar) = - (  F1_XT*dsdx_r(0,0)
                                        + 2*F11_XT*StressGP(0,0)*dsdx_r(0,0)
                                        + 2*F12_XT*(dsdx_r(0,0)*(StressGP(1,1)+StressGP(2,2))
                                        + StressGP(0,0)*(dsdx_r(1,1)+dsdx_r(2,2))));
              hess_lsf(jvar, ivar) = hess_lsf(ivar, jvar);
            }
            else if ((vars_name[ivar+1].compare(0,2,"XT")==0) && (vars_name[jvar+1].compare(0,2,"XT")==0)) {
              // XT XT
              double F12_XTXT = -3/8*pow(X_T*X_C*Y_T*Y_C,-2.5)*pow(X_C*Y_T*Y_C,2);
              double F1_XTXT  = 2/(X_T*X_T*X_T);
              double F11_XTXT = 2/(X_T*X_T*X_T*X_C);
              hess_lsf(ivar, jvar) = - (  F1_XTXT*StressGP(0,0)
                                 + F11_XTXT*StressGP(0,0)*StressGP(0,0)
                                 + 2*F12_XTXT*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
            }
            else if ((ivar<num_ply_vars) && (vars_name[jvar+1].compare(0,2,"XC")==0)) {
              // X_C Em
              double F1_XC = 1/(X_C*X_C);
              double F11_XC = -1/(X_T*X_C*X_C);
              double F12_XC = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*Y_T*Y_C;

              hess_lsf(ivar, jvar) = - (  F1_XC*dsdx_r(0,0)
                                 + 2*F11_XC*StressGP(0,0)*dsdx_r(0,0)
                                 + 2*F12_XC*(dsdx_r(0,0)*(StressGP(1,1)+StressGP(2,2))
                                             + StressGP(0,0)*(dsdx_r(1,1)+dsdx_r(2,2))));
              hess_lsf(jvar, ivar) = hess_lsf(ivar, jvar);
            }
            else if ((vars_name[ivar+1].compare(0,2,"XC")==0) && (vars_name[jvar+1].compare(0,2,"XC")==0)) {
              // X_C X_C
              double F1_XCXC  = -2/(X_C*X_C*X_C);
              double F11_XCXC = 2/(X_T*X_C*X_C*X_C);
              double F12_XCXC = -3/8*pow(X_T*X_C*Y_T*Y_C,-2.5)*pow(X_T*Y_T*Y_C,2);
              hess_lsf(ivar, jvar) = - (  F1_XCXC*StressGP(0,0)
                                 + F11_XCXC*StressGP(0,0)*StressGP(0,0)
                                 + 2*F12_XCXC*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
            }
            else if ((ivar<num_ply_vars) && (vars_name[jvar+1].compare(0,2,"YT")==0)) {
              // Y_T Em
              double F2_YT = -1/(Y_T*Y_T);
              double F22_YT = -1/(Y_T*Y_T*Y_C);
              double F12_YT = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*X_C*Y_C;
              double F44_YT = -2*eta*Y_C/(S_23*S_23*S_23);
              double F23_YT = 1/(2*Y_C*Y_T*Y_T);
              
              hess_lsf(ivar, jvar) = - (  F2_YT*(dsdx_r(1,1) + dsdx_r(2,2))
                                        + 2*F22_YT*(StressGP(1,1)*dsdx_r(1,1) + StressGP(2,2)*dsdx_r(2,2))
                                        + 2*F44_YT*StressGP(1,2)*dsdx_r(1,2)
                                        + 2*F12_YT*(dsdx_r(0,0)*(StressGP(1,1)+StressGP(2,2)) + StressGP(0,0)*(dsdx_r(1,1)+dsdx_r(2,2)))
                                        + 2*F23_YT*(dsdx_r(1,1)*StressGP(2,2) + StressGP(1,1)*dsdx_r(2,2)));
              hess_lsf(jvar, ivar) = hess_lsf(ivar, jvar);
            }
            else if ((vars_name[ivar+1].compare(0,2,"YT")==0) && (vars_name[jvar+1].compare(0,2,"YT")==0)) {
              // Y_T Y_T
              double F2_YTYT  = 2/(Y_T*Y_T*Y_T);
              double F22_YTYT = 2/(Y_T*Y_T*Y_T*Y_C);
              double F12_YTYT = -3/8*pow(X_T*X_C*Y_T*Y_C,-2.5)*pow(X_T*X_C*Y_C,2);
              double F44_YTYT = 6*pow(eta*Y_C,2)/pow(S_23,4);
              double F23_YTYT = -2/(2*Y_C*Y_T*Y_T*Y_T);
              
              hess_lsf(ivar, jvar) = - (  F2_YTYT*(StressGP(1,1) + StressGP(2,2))
                                        + F22_YTYT*(StressGP(1,1)*StressGP(1,1) + StressGP(2,2)*StressGP(2,2))
                                        + F44_YTYT*StressGP(1,2)*StressGP(1,2)
                                        + 2*F12_YTYT*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2))
                                        + 2*F23_YTYT*StressGP(1,1)*StressGP(2,2));
            }
            else if ((ivar<num_ply_vars) && (vars_name[jvar+1].compare(0,2,"YC")==0)) {
              // Y_C EM
              double F2_YC = 1/(Y_C*Y_C);
              double F22_YC = -1/(Y_T*Y_C*Y_C);
              double F12_YC = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*X_C*Y_T;
              double F44_YC = -2*eta*Y_T/(S_23*S_23*S_23);
              double F23_YC = 1/(2*Y_C*Y_C*Y_T);
              
              hess_lsf(ivar, jvar) = - (  F2_YC*(dsdx_r(1,1) + dsdx_s(2,2))
                                        + 2*F22_YC*(StressGP(1,1)*dsdx_r(1,1) + StressGP(2,2)*dsdx_r(2,2))
                                        + 2*F44_YC*StressGP(1,2)*dsdx_r(1,2)
                                        + 2*F12_YC*(dsdx_r(0,0)*(StressGP(1,1)+StressGP(2,2)) + StressGP(0,0)*(dsdx_r(1,1)+dsdx_r(2,2)))
                                        + 2*F23_YC*(dsdx_r(1,1)*StressGP(2,2) + StressGP(1,1)*dsdx_r(2,2)));
              hess_lsf(jvar, ivar) = hess_lsf(ivar, jvar);
            }
            else if ((vars_name[ivar+1].compare(0,2,"YC")==0) && (vars_name[jvar+1].compare(0,2,"YC")==0)) {
              // Y_C Y_C
              double F2_YCYC  = -2/(Y_C*Y_C*Y_C);
              double F22_YCYC = 2/(Y_T*Y_C*Y_C*Y_C);
              double F12_YCYC = -3/8*pow(X_T*X_C*Y_T*Y_C,-2.5)*pow(X_T*X_C*Y_T,2);
              double F44_YCYC = 6*pow(eta*Y_T,2)/pow(S_23,4);
              double F23_YCYC = -2/(2*Y_C*Y_C*Y_C*Y_T);
              
              hess_lsf(ivar, jvar) = - (  F2_YCYC*(StressGP(1,1) + StressGP(2,2))
                                        + F22_YCYC*(StressGP(1,1)*StressGP(1,1) + StressGP(2,2)*StressGP(2,2))
                                        + F44_YCYC*StressGP(1,2)*StressGP(1,2)
                                        + 2*F12_YCYC*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2))
                                        + 2*F23_YCYC*StressGP(1,1)*StressGP(2,2));
            }
            else if ((ivar<num_ply_vars) && (vars_name[jvar+1].compare(0,3,"S12")==0)) {
              // S12 Em
              hess_lsf(ivar, jvar) = 4*(StressGP(0,1)*dsdx_r(0,1) + StressGP(0,2)*dsdx_r(0,2))/pow(S_12,3);
              hess_lsf(jvar, ivar) = hess_lsf(ivar, jvar);
            }
            else if ((vars_name[ivar+1].compare(0,3,"S12")==0) && (vars_name[jvar+1].compare(0,2,"S12")==0)) {
              // S12 S12
              hess_lsf(ivar, jvar) = -6*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2))/pow(S_12,4);
            }
          }
        }
      }
      
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 44050
    // TSAI-WU
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Wu_Christensen(ublas::vector<double> x,
                                                        vector<string> vars_name,
                                                        ublas::vector<double> MatStrength,
                                                        ublas::matrix<double> StressGP,
                                                        ublas::matrix<double> StressGP_r_Em,
                                                        ublas::matrix<double> StressGP_r_NUm,
                                                        ublas::matrix<double> StressGP_r_NUp,
                                                        ublas::matrix<double> StressGP_r_NUpz,
                                                        ublas::matrix<double> StressGP_r_Ep,
                                                        ublas::matrix<double> StressGP_r_Ez,
                                                        ublas::matrix<double> StressGP_r_Gzp,
                                                        ublas::matrix<double> StressGP_r_Ef,
                                                        ublas::matrix<double> StressGP_r_NUf,
                                                        ublas::matrix<double> StressGP_r_F,
                                                        ublas::matrix<double> StressGP_r_Theta,
                                                        ublas::matrix<double> StressGP_r_Theta_1,
                                                        ublas::matrix<double> StressGP_r_Theta_2,
                                                        ublas::matrix<double> StressGP_r_Theta_3,
                                                        ublas::matrix<double> StressGP_r_Theta_4,
                                                        double &val_lsf,
                                                        ublas::vector<double> &grad_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      double F1, F2, F11, F22, F44, F66, F12;
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = MatStrength(0);
      X_C  = MatStrength(1);
      Y_T  = MatStrength(2);
      Y_C  = MatStrength(3);//cout<<"The Y_T is: "<<Y_T<<endl;
      S_12 = MatStrength(4);
      S_23 = MatStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      // Calculate the transverse shear strength using Christensen's formula
      double eta;
      eta = 0.2519;
      S_23 = eta*Y_T*Y_C;
      
      F1 = (1/X_T - 1/X_C);
      F2 = (1/Y_T - 1/Y_C);
      F11 = 1/(X_T*X_C);
      F22 = 1/(Y_T*Y_C);
      F12 = -1/(2*sqrt(X_T*X_C*Y_T*Y_C)); // Empirical suggestion: Mises-Hencky criterion
      F44 = 1/(S_23*S_23);
      F66 = 1/(S_12*S_12);
      
      val_lsf = 1- (  F1*StressGP(0,0)
                    + F2*(StressGP(1,1) + StressGP(2,2))
                    + F11*StressGP(0,0)*StressGP(0,0)
                    + F22*pow((StressGP(1,1) + StressGP(2,2)),2)
                    + F44*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2))
                    + F66*(StressGP(0,1)*StressGP(0,1) + StressGP(2,0)*StressGP(2,0))
                    + 2*F12*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      for (unsigned i=1; i<=x.size(); i++) {
        cout<<"The variable name is "<<vars_name[i]<<endl;
        if (vars_name[i].compare(0,2,"Em") == 0) {
          // w.r.t. Em
          dsdx.clear(); dsdx = StressGP_r_Em;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,3,"NUm") == 0) {
          // w.r.t. NUm
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,3,"NUp") == 0) {
          // w.r.t. NUp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,3,"NUz") == 0) {
          // w.r.t. NUpz
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
          
        }
        else if (vars_name[i].compare(0,2,"Ep") == 0) {
          // w.r.t. Ep
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,2,"Ez") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,3,"Gzp") == 0) {
          // w.r.t. Gzp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,2,"Ef") == 0) {
          // w.r.t. Ef - fibre with isotropic material
          dsdx.clear(); dsdx = StressGP_r_Ef;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,3,"NUf") == 0) {
          // w.r.t. NUf - fibre with isotropic material
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,5,"force") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,11,"orientation") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,6,"theta1") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,6,"theta2") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,6,"theta3") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,6,"theta4") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
          grad_lsf(i-1) = - (  F1*dsdx(0,0)
                             + F2*(dsdx(1,1) + dsdx(2,2))
                             + 2*F11*StressGP(0,0)*dsdx(0,0)
                             + 2*F22*(StressGP(1,1) + StressGP(2,2))*(dsdx(1,1) + dsdx(2,2))
                             + F44*(2*StressGP(1,2)*dsdx(1,2) - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))
                             + F66*(2*StressGP(0,1)*dsdx(0,1) + 2*StressGP(2,0)*dsdx(2,0))
                             + 2*F12*dsdx(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + 2*F12*StressGP(0,0)*(dsdx(1,1)+dsdx(2,2)));
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
          //grad_lsf(i-1) = - ( F2_YT*(StressGP(1,1) + StressGP(2,2))
          //                   + F22_YT*pow((StressGP(1,1) + StressGP(2,2)),2)
          //                   + 2*F12_YT*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
          double F44_YT = -2*eta*Y_C/(S_23*S_23*S_23);
          grad_lsf(i-1) = - ( F2_YT*(StressGP(1,1) + StressGP(2,2))
                             + F22_YT*pow((StressGP(1,1) + StressGP(2,2)),2)
                             + 2*F12_YT*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + F44_YT*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2)));
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          // w.r.t Y_C
          double F2_YC = 1/(Y_C*Y_C);
          double F22_YC = -1/(Y_T*Y_C*Y_C);
          double F12_YC = 1/4*pow(X_T*X_C*Y_T*Y_C,-1.5)*X_T*X_C*Y_T;
          //grad_lsf(i-1) = - ( F2_YC*(StressGP(1,1) + StressGP(2,2))
          //                   + F22_YC*pow((StressGP(1,1) + StressGP(2,2)),2)
          //                   + 2*F12_YC*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2)));
          double F44_YC = -2*eta*Y_T/(S_23*S_23*S_23);
          grad_lsf(i-1) = - ( F2_YC*(StressGP(1,1) + StressGP(2,2))
                             + F22_YC*pow((StressGP(1,1) + StressGP(2,2)),2)
                             + 2*F12_YC*StressGP(0,0)*(StressGP(1,1)+StressGP(2,2))
                             + F44_YC*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2)));
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          // w.r.t S_12
          grad_lsf(i-1) = 2*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2))/pow(S_12,3);
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 42060
    // Tsai-Hill failure theory
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Hill_2D(ublas::vector<double> x,
                                                 vector<string> vars_name,
                                                 ublas::vector<double> PlyStrength,
                                                 ublas::matrix<double> StressGP,
                                                 ublas::matrix<double> StressGP_r_Em,
                                                 ublas::matrix<double> StressGP_r_NUm,
                                                 ublas::matrix<double> StressGP_r_NUp,
                                                 ublas::matrix<double> StressGP_r_NUpz,
                                                 ublas::matrix<double> StressGP_r_Ep,
                                                 ublas::matrix<double> StressGP_r_Ez,
                                                 ublas::matrix<double> StressGP_r_Gzp,
                                                 ublas::matrix<double> StressGP_r_Ef,
                                                 ublas::matrix<double> StressGP_r_NUf,
                                                 ublas::matrix<double> StressGP_r_F,
                                                 ublas::matrix<double> StressGP_r_Theta,
                                                 ublas::matrix<double> StressGP_r_Theta_1,
                                                 ublas::matrix<double> StressGP_r_Theta_2,
                                                 ublas::matrix<double> StressGP_r_Theta_3,
                                                 ublas::matrix<double> StressGP_r_Theta_4,
                                                 double &val_lsf,
                                                 ublas::vector<double> &grad_lsf) {
      
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double X,Y;
      
      
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      for (unsigned i=1; i<=x.size(); i++) {
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
        else if (vars_name[i].compare(0,3,"NUz") == 0) {
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
        else if (vars_name[i].compare(0,3,"Gzp") == 0) {
          // w.r.t. Gzp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
          grad_lsf(i-1) = - (2*StressGP(0,0)*StressGP_r_Gzp(0,0)/(X*X)
                             + 2*StressGP(1,1)*StressGP_r_Gzp(1,1)/(Y*Y)
                             + 2*StressGP(0,1)*StressGP_r_Gzp(0,1)/(S_12*S_12)
                             - StressGP_r_Gzp(0,0)*StressGP(1,1)/(X*X)
                             - StressGP(0,0)*StressGP_r_Gzp(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,2,"Ef") == 0) {
          // w.r.t. Ef - fibre with isotropic material
          dsdx.clear(); dsdx = StressGP_r_Ef;
          grad_lsf(i-1) = - (2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y*Y)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X*X)
                             - StressGP(0,0)*dsdx(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,3,"NUf") == 0) {
          // w.r.t. NUf - fibre with isotropic material
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
          grad_lsf(i-1) = - (2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y*Y)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X*X)
                             - StressGP(0,0)*dsdx(1,1)/(X*X));
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
          // w.r.t. F
          dsdx.clear(); dsdx = StressGP_r_F;
          grad_lsf(i-1) = - (2*StressGP(0,0)*StressGP_r_F(0,0)/(X*X)
                             + 2*StressGP(1,1)*StressGP_r_F(1,1)/(Y*Y)
                             + 2*StressGP(0,1)*StressGP_r_F(0,1)/(S_12*S_12)
                             - StressGP_r_F(0,0)*StressGP(1,1)/(X*X)
                             - StressGP(0,0)*StressGP_r_F(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,11,"orientation") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta;
          grad_lsf(i-1) = - (2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y*Y)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X*X)
                             - StressGP(0,0)*dsdx(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,6,"theta1") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_1;
          grad_lsf(i-1) = - (2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y*Y)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X*X)
                             - StressGP(0,0)*dsdx(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,6,"theta2") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_2;
          grad_lsf(i-1) = - (2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y*Y)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X*X)
                             - StressGP(0,0)*dsdx(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,6,"theta3") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_3;
          grad_lsf(i-1) = - (2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y*Y)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X*X)
                             - StressGP(0,0)*dsdx(1,1)/(X*X));
        }
        else if (vars_name[i].compare(0,6,"theta4") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_4;
          grad_lsf(i-1) = - (2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y*Y)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X*X)
                             - StressGP(0,0)*dsdx(1,1)/(X*X));
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 43060
    // Tsai-Hill failure theory
    //
    
    virtual PetscErrorCode gfun_ply_Tsai_Hill(ublas::vector<double> x,
                                              vector<string> vars_name,
                                              ublas::vector<double> PlyStrength,
                                              ublas::matrix<double> StressGP,
                                              ublas::matrix<double> StressGP_r_Em,
                                              ublas::matrix<double> StressGP_r_NUm,
                                              ublas::matrix<double> StressGP_r_NUp,
                                              ublas::matrix<double> StressGP_r_NUpz,
                                              ublas::matrix<double> StressGP_r_Ep,
                                              ublas::matrix<double> StressGP_r_Ez,
                                              ublas::matrix<double> StressGP_r_Gzp,
                                              ublas::matrix<double> StressGP_r_Ef,
                                              ublas::matrix<double> StressGP_r_NUf,
                                              ublas::matrix<double> StressGP_r_F,
                                              ublas::matrix<double> StressGP_r_Theta,
                                              ublas::matrix<double> StressGP_r_Theta_1,
                                              ublas::matrix<double> StressGP_r_Theta_2,
                                              ublas::matrix<double> StressGP_r_Theta_3,
                                              ublas::matrix<double> StressGP_r_Theta_4,
                                              double &val_lsf,
                                              ublas::vector<double> &grad_lsf) {
      
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      double X,Y;
      
      
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      S_23 = PlyStrength(5);
      
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      // Calculate the transverse shear strength using Christensen's formula
      double eta;
      eta = 0.2519;
      S_23 = eta*Y_T*Y_C;
      
      val_lsf = 1- (  StressGP(0,0)*StressGP(0,0)/(X*X)
                    + pow((StressGP(1,1) - StressGP(2,2)),2)/(Y*Y)
                    + (StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2))/(S_12*S_12)
                    + StressGP(1,2)*StressGP(1,2)/(S_23*S_23)
                    - (StressGP(0,0)*StressGP(1,1) + StressGP(0,0)*StressGP(2,2) - StressGP(1,1)*StressGP(2,2))/(X*X));
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      for (unsigned i=1; i<=x.size(); i++) {
        cout<<"The variable name is "<<vars_name[i]<<endl;
        if (vars_name[i].compare(0,2,"Em") == 0) {
          // w.r.t. Em
          dsdx.clear(); dsdx = StressGP_r_Em;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,3,"NUm") == 0) {
          // w.r.t. NUm
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,3,"NUp") == 0) {
          // w.r.t. NUp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,3,"NUz") == 0) {
          // w.r.t. NUpz
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,2,"Ep") == 0) {
          // w.r.t. Ep
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,2,"Ez") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,3,"Gzp") == 0) {
          // w.r.t. Gzp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,2,"Ef") == 0) {
          // w.r.t. Ef - fibre with isotropic material
          dsdx.clear(); dsdx = StressGP_r_Ef;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,3,"NUf") == 0) {
          // w.r.t. NUf - fibre with isotropic material
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,5,"force") == 0) {
          // w.r.t. Force
          dsdx.clear(); dsdx = StressGP_r_F;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,11,"orientation") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,6,"theta1") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_1;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,6,"theta2") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_2;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,6,"theta3") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_3;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,6,"theta4") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_4;
          grad_lsf(i-1) = - (  2*StressGP(0,0)*dsdx(0,0)/(X*X)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y*Y)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1)
                                + dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2)
                                - dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X*X));
        }
        else if (vars_name[i].compare(0,2,"XT") == 0) {
          // w.r.t. XT
          if (StressGP(0,0)>0) {
            grad_lsf(i-1) =  2*StressGP(0,0)*StressGP(0,0)/(X*X*X)
                            -2*(StressGP(0,0)*StressGP(1,1) + StressGP(0,0)*StressGP(2,2) - StressGP(1,1)*StressGP(2,2))/(X*X*X);
          } else {
            grad_lsf(i-1) = 0;
          }
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          // w.r.t. XC
          if (StressGP(0,0)<0) {
            grad_lsf(i-1) =  2*StressGP(0,0)*StressGP(0,0)/(X*X*X)
                            -2*(StressGP(0,0)*StressGP(1,1) + StressGP(0,0)*StressGP(2,2) - StressGP(1,1)*StressGP(2,2))/(X*X*X);
          } else {
            grad_lsf(i-1) = 0;
          }
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          // w.r.t. YT
          if (StressGP(1,1)>0) {
            grad_lsf(i-1) = (2*pow((StressGP(1,1) - StressGP(2,2)),2)/(Y*Y*Y)
                            +2*eta*Y_C*StressGP(1,2)*StressGP(1,2)/(S_23*S_23*S_23));
          } else {
            grad_lsf(i-1) = 0;
          }
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          // w.r.t. YC
          if (StressGP(1,1)<0) {
            grad_lsf(i-1) = (2*pow((StressGP(1,1) - StressGP(2,2)),2)/(Y*Y*Y)
                             +2*eta*Y_T*StressGP(1,2)*StressGP(1,2)/(S_23*S_23*S_23));
          } else {
            grad_lsf(i-1) = 0;
          }
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          // w.r.t. S12
          grad_lsf(i-1) = 2*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2))/(S_12*S_12*S_12);
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 13073
    // RICHARD CHRISTENSEN: FIBRE CONTROLLED FAILURE
    //
    
    virtual PetscErrorCode gfun_ply_RCF(ublas::vector<double> x,
                                        vector<string> vars_name,
                                        ublas::vector<double> PlyStrength,
                                        ublas::matrix<double> StressGP,
                                        ublas::matrix<double> StressGP_r_Em,
                                        ublas::matrix<double> StressGP_r_NUm,
                                        ublas::matrix<double> StressGP_r_NUp,
                                        ublas::matrix<double> StressGP_r_NUpz,
                                        ublas::matrix<double> StressGP_r_Ep,
                                        ublas::matrix<double> StressGP_r_Ez,
                                        ublas::matrix<double> StressGP_r_Gzp,
                                        ublas::matrix<double> StressGP_r_Ef,
                                        ublas::matrix<double> StressGP_r_NUf,
                                        ublas::matrix<double> StressGP_r_F,
                                        ublas::matrix<double> StressGP_r_Theta,
                                        ublas::matrix<double> StressGP_r_Theta_1,
                                        ublas::matrix<double> StressGP_r_Theta_2,
                                        ublas::matrix<double> StressGP_r_Theta_3,
                                        ublas::matrix<double> StressGP_r_Theta_4,
                                        double &val_lsf,
                                        ublas::vector<double> &grad_lsf) {
      
      PetscFunctionBegin;
      
      double T_11; // tensile strength in the fibre direction
      double T_22; // tensile strength in the transverse direction
      double C_11; // compressive strength in the fibre direction
      double C_22; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      T_11 = PlyStrength(0);
      C_11 = PlyStrength(1);
      T_22 = PlyStrength(2);
      C_22 = PlyStrength(3);
      S_12 = PlyStrength(4);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
        if (vars_name[i].compare(0,2,"XT") == 0) {
          T_11 = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          C_11 = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          T_22 = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          C_22 = x(i-1);
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          S_12 = x(i-1);
        }
      }
      
      
      // Evaluate the limit state function
      val_lsf = 1 - ((1/T_11 - 1/C_11)*StressGP(0,0)
                     + 1/(T_11*C_11)*StressGP(0,0)*StressGP(0,0));
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      for (unsigned i=1; i<=x.size(); i++) {
        cout<<"The variable name is "<<vars_name[i]<<endl;
        if (vars_name[i].compare(0,2,"Em") == 0) {
          // w.r.t. Em
          dsdx.clear(); dsdx = StressGP_r_Em;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,3,"NUm") == 0) {
          // w.r.t. NUm
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,3,"NUp") == 0) {
          // w.r.t. NUp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,3,"NUz") == 0) {
          // w.r.t. NUpz
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,2,"Ep") == 0) {
          // w.r.t. Ep
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,2,"Ez") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,3,"Gzp") == 0) {
          // w.r.t. Gzp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,2,"Ef") == 0) {
          // w.r.t. Ef - fibre with isotropic material
          dsdx.clear(); dsdx = StressGP_r_Ef;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,3,"NUf") == 0) {
          // w.r.t. NUf - fibre with isotropic material
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,5,"force") == 0) {
          // w.r.t. F
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,11,"orientation") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,6,"theta1") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,6,"theta2") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,6,"theta3") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,6,"theta4") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
          grad_lsf(i-1) = - ((1/T_11 - 1/C_11)*dsdx(0,0)
                             + 2/(T_11*C_11)*StressGP(0,0)*dsdx(0,0));
        }
        else if (vars_name[i].compare(0,2,"XT") == 0) {
          // w.r.t X_T
          grad_lsf(i-1) =  1/(T_11*T_11)*StressGP(0,0)
                         + 1/(T_11*T_11*C_11)*StressGP(0,0)*StressGP(0,0);
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          // w.r.t X_C
          grad_lsf(i-1) = - 1/(C_11*C_11)*StressGP(0,0)
                          + 1/(T_11*C_11*C_11)*StressGP(0,0)*StressGP(0,0);
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          // w.r.t Y_T
          grad_lsf(i-1) = 0;
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          // w.r.t Y_C
          grad_lsf(i-1) = 0;
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          // w.r.t S_12
          grad_lsf(i-1) = 0;
        }
        
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 23073
    // RICHARD CHRISTENSEN: MATRIX CONTROLLED FAILURE
    //
    
    virtual PetscErrorCode gfun_ply_RCM(ublas::vector<double> x,
                                        vector<string> vars_name,
                                        ublas::vector<double> PlyStrength,
                                        ublas::matrix<double> StressGP,
                                        ublas::matrix<double> StressGP_r_Em,
                                        ublas::matrix<double> StressGP_r_NUm,
                                        ublas::matrix<double> StressGP_r_NUp,
                                        ublas::matrix<double> StressGP_r_NUpz,
                                        ublas::matrix<double> StressGP_r_Ep,
                                        ublas::matrix<double> StressGP_r_Ez,
                                        ublas::matrix<double> StressGP_r_Gzp,
                                        ublas::matrix<double> StressGP_r_Ef,
                                        ublas::matrix<double> StressGP_r_NUf,
                                        ublas::matrix<double> StressGP_r_F,
                                        ublas::matrix<double> StressGP_r_Theta,
                                        ublas::matrix<double> StressGP_r_Theta_1,
                                        ublas::matrix<double> StressGP_r_Theta_2,
                                        ublas::matrix<double> StressGP_r_Theta_3,
                                        ublas::matrix<double> StressGP_r_Theta_4,
                                        double &val_lsf,
                                        ublas::vector<double> &grad_lsf) {
      
      PetscFunctionBegin;
      
      double T_11; // tensile strength in the fibre direction
      double T_22; // tensile strength in the transverse direction
      double C_11; // compressive strength in the fibre direction
      double C_22; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      T_11 = PlyStrength(0);
      C_11 = PlyStrength(1);
      T_22 = PlyStrength(2);
      C_22 = PlyStrength(3);
      S_12 = PlyStrength(4);
      S_23 = PlyStrength(5);
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
        if (vars_name[i].compare(0,2,"XT") == 0) {
          T_11 = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          C_11 = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          T_22 = x(i-1);
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          C_22 = x(i-1);
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          S_12 = x(i-1);
        }
        else if (vars_name[i].compare(0,3,"S23") == 0) {
          S_23 = x(i-1);
        }
      }
      
      // Calculate the transverse shear strength using Christensen's formula
      double eta;
      eta = 0.2519;
      S_23 = eta*T_22*C_22;
    
      // Evaluate the limit state function
      val_lsf = 1 - ((1/T_22 - 1/C_22)*(StressGP(1,1) + StressGP(2,2))
                     + 1/(T_22*C_22)*pow((StressGP(1,1)+ StressGP(2,2)),2)
                     + 1/(S_12*S_12)*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2))
                     + 1/(S_23*S_23)*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2)));
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      for (unsigned i=1; i<=x.size(); i++) {
        cout<<"The variable name is "<<vars_name[i]<<endl;
        if (vars_name[i].compare(0,2,"Em") == 0) {
          // w.r.t. Em
          dsdx.clear(); dsdx = StressGP_r_Em;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,3,"NUm") == 0) {
          // w.r.t. NUm
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,3,"NUp") == 0) {
          // w.r.t. NUp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,3,"NUz") == 0) {
          // w.r.t. NUpz
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,2,"Ep") == 0) {
          // w.r.t. Ep
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,2,"Ez") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,3,"Gzp") == 0) {
          // w.r.t. Gzp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,2,"Ef") == 0) {
          // w.r.t. Ef - fibre with isotropic material
          dsdx.clear(); dsdx = StressGP_r_Ef;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,3,"NUf") == 0) {
          // w.r.t. NUf - fibre with isotropic material
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,5,"force") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_F;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,11,"orientation") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,6,"theta1") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_1;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,6,"theta2") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_2;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,6,"theta3") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_3;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,6,"theta4") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Theta_4;
          grad_lsf(i-1) = - ((1/T_22 - 1/C_22)*(dsdx(1,1) + dsdx(2,2))
                             + 2/(T_22*C_22)*(StressGP(1,1)+ StressGP(2,2))*(dsdx(1,1)+ dsdx(2,2))
                             + 2/(S_12*S_12)*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))
                             + 1/(S_23*S_23)*(2*StressGP(1,2)*dsdx(1,2)
                                              - dsdx(1,1)*StressGP(2,2)
                                              - StressGP(1,1)*dsdx(2,2)));
        }
        else if (vars_name[i].compare(0,2,"XT") == 0) {
          // w.r.t X_T
          grad_lsf(i-1) = 0;
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          // w.r.t X_C
          grad_lsf(i-1) = 0;
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          // w.r.t Y_T
          grad_lsf(i-1) = - (- 1/(T_22*T_22)*(StressGP(1,1) + StressGP(2,2))
                             - 1/(T_22*C_22*T_22)*pow((StressGP(1,1)+ StressGP(2,2)),2)
                             - 2*eta*C_22/(S_23*S_23*S_23)*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2)));
          
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          // w.r.t Y_C
          grad_lsf(i-1) = - (1/(C_22*C_22)*(StressGP(1,1) + StressGP(2,2))
                             - 1/(T_22*C_22*C_22)*pow((StressGP(1,1)+ StressGP(2,2)),2)
                             - 2*eta*T_22/(S_23*S_23*S_23)*(StressGP(1,2)*StressGP(1,2) - StressGP(1,1)*StressGP(2,2)));
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          // w.r.t S_12
          grad_lsf(i-1) = 2/(S_12*S_12*S_12)*(StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2));
        }
        
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 42080
    // Hoffman failure theory
    
    virtual PetscErrorCode gfun_ply_Hoffman_2D(ublas::vector<double> x,
                                               vector<string> vars_name,
                                               ublas::vector<double> PlyStrength,
                                               ublas::matrix<double> StressGP,
                                               ublas::matrix<double> StressGP_r_Em,
                                               ublas::matrix<double> StressGP_r_NUm,
                                               ublas::matrix<double> StressGP_r_NUp,
                                               ublas::matrix<double> StressGP_r_NUpz,
                                               ublas::matrix<double> StressGP_r_Ep,
                                               ublas::matrix<double> StressGP_r_Ez,
                                               ublas::matrix<double> StressGP_r_Gzp,
                                               ublas::matrix<double> StressGP_r_Ef,
                                               ublas::matrix<double> StressGP_r_NUf,
                                               ublas::matrix<double> StressGP_r_F,
                                               ublas::matrix<double> StressGP_r_Theta,
                                               ublas::matrix<double> StressGP_r_Theta_1,
                                               ublas::matrix<double> StressGP_r_Theta_2,
                                               ublas::matrix<double> StressGP_r_Theta_3,
                                               ublas::matrix<double> StressGP_r_Theta_4,
                                               double &val_lsf,
                                               ublas::vector<double> &grad_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      val_lsf = 1 - ( StressGP(0,0)*(1/X_T - 1/X_C)
                     + StressGP(1,1)*(1/Y_T - 1/Y_C)
                     + StressGP(0,0)*StressGP(0,0)/(X_T*X_C)
                     + StressGP(1,1)*StressGP(1,1)/(Y_T*Y_C)
                     + StressGP(0,1)*StressGP(0,1)/(S_12*S_12)
                     - StressGP(0,0)*StressGP(1,1)/(X_T*X_C));
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      for (unsigned i=1; i<=x.size(); i++) {
        cout<<"The variable name is "<<vars_name[i]<<endl;
        if (vars_name[i].compare(0,2,"Em") == 0) {
          // w.r.t. Em
          dsdx.clear(); dsdx = StressGP_r_Em;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,3,"NUm") == 0) {
          // w.r.t. NUm
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,3,"NUp") == 0) {
          // w.r.t. NUp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,3,"NUz") == 0) {
          // w.r.t. NUpz
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,2,"Ep") == 0) {
          // w.r.t. Ep
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,2,"Ez") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,3,"Gzp") == 0) {
          // w.r.t. Gzp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,2,"Ef") == 0) {
          // w.r.t. Ef - fibre with isotropic material
          dsdx.clear(); dsdx = StressGP_r_Ef;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,3,"NUf") == 0) {
          // w.r.t. NUf - fibre with isotropic material
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,2,"XT") == 0) {
          // w.r.t. XT
          grad_lsf(i-1) =  - (  StressGP(0,0)*(-1/(X_T*X_T))
                              - StressGP(0,0)*StressGP(0,0)/(X_T*X_T*X_C)
                              + StressGP(0,0)*StressGP(1,1)/(X_T*X_T*X_C));
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          // w.r.t. XC
          grad_lsf(i-1) =  - ( StressGP(0,0)*( 1/(X_C*X_C))
                              - StressGP(0,0)*StressGP(0,0)/(X_T*X_C*X_C)
                              + StressGP(0,0)*StressGP(1,1)/(X_T*X_C*X_C));
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          // w.r.t. YT
          grad_lsf(i-1) = - (  StressGP(1,1)*(-1/(Y_T*Y_T))
                             - StressGP(1,1)*StressGP(1,1)/(Y_T*Y_T*Y_C));
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          // w.r.t. YC
          grad_lsf(i-1) = - (  StressGP(1,1)*( 1/(Y_C*Y_C))
                             - StressGP(1,1)*StressGP(1,1)/(Y_T*Y_C*Y_C));
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          // w.r.t. S12
          grad_lsf(i-1) = 2*StressGP(0,1)*StressGP(0,1)/(S_12*S_12*S_12);
        }
        else if (vars_name[i].compare(0,5,"force") == 0) {
          // w.r.t. Em
          dsdx.clear(); dsdx = StressGP_r_F;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,11,"orientation") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,6,"theta1") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_1;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,6,"theta2") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_2;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,6,"theta3") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_3;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,6,"theta4") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_4;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + dsdx(1,1)*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*StressGP(1,1)*dsdx(1,1)/(Y_T*Y_C)
                             + 2*StressGP(0,1)*dsdx(0,1)/(S_12*S_12)
                             - dsdx(0,0)*StressGP(1,1)/(X_T*X_C)
                             - StressGP(0,0)*dsdx(1,1)/(X_T*X_C));
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    //------------------------------------------------------------------------------
    // Code: 43080
    // Hoffman failure theory
    
    virtual PetscErrorCode gfun_ply_Hoffman(ublas::vector<double> x,
                                            vector<string> vars_name,
                                            ublas::vector<double> PlyStrength,
                                            ublas::matrix<double> StressGP,
                                            ublas::matrix<double> StressGP_r_Em,
                                            ublas::matrix<double> StressGP_r_NUm,
                                            ublas::matrix<double> StressGP_r_NUp,
                                            ublas::matrix<double> StressGP_r_NUpz,
                                            ublas::matrix<double> StressGP_r_Ep,
                                            ublas::matrix<double> StressGP_r_Ez,
                                            ublas::matrix<double> StressGP_r_Gzp,
                                            ublas::matrix<double> StressGP_r_Ef,
                                            ublas::matrix<double> StressGP_r_NUf,
                                            ublas::matrix<double> StressGP_r_F,
                                            ublas::matrix<double> StressGP_r_Theta,
                                            ublas::matrix<double> StressGP_r_Theta_1,
                                            ublas::matrix<double> StressGP_r_Theta_2,
                                            ublas::matrix<double> StressGP_r_Theta_3,
                                            ublas::matrix<double> StressGP_r_Theta_4,
                                            double &val_lsf,
                                            ublas::vector<double> &grad_lsf) {
      PetscFunctionBegin;
      
      double X_T; // tensile strength in the fibre direction
      double X_C; // tensile strength in the transverse direction
      double Y_T; // compressive strength in the fibre direction
      double Y_C; // compressive strength in the transverse direction
      double S_12; // Longitudinal or axial shear strength
      double S_23; // Transverse shear strength
      
      ublas::matrix<double> dsdx(3,3);
      //cout<<"\nPly strength: "<<PlyStrength<<endl;
      X_T  = PlyStrength(0);
      X_C  = PlyStrength(1);
      Y_T  = PlyStrength(2);
      Y_C  = PlyStrength(3);
      S_12 = PlyStrength(4);
      S_23 = PlyStrength(5);
      
      
      // Update the values if the parameters are considered as random variables.
      for (unsigned i=1; i<=x.size();i++) {
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
      
      // Calculate the transverse shear strength using Christensen's formula
      double eta;
      eta = 0.2519;
      S_23 = eta*Y_T*Y_C;
      
      // Evaluate the limit state function
      val_lsf = 1 - (  StressGP(0,0)*(1/X_T - 1/X_C)
                     + (StressGP(1,1) + StressGP(2,2))*(1/Y_T - 1/Y_C)
                     + StressGP(0,0)*StressGP(0,0)/(X_T*X_C)
                     + pow((StressGP(1,1) - StressGP(2,2)),2)/(Y_T*Y_C)
                     + (StressGP(0,1)*StressGP(0,1) + StressGP(0,2)*StressGP(0,2))/(S_12*S_12)
                     + StressGP(1,2)*StressGP(1,2)/(S_23*S_23)
                     - (StressGP(0,0)*StressGP(1,1) + StressGP(0,0)*StressGP(2,2) - StressGP(1,1)*StressGP(2,2))/(X_T*X_C));
      
      // Evaluate the partial derovative of the LSF w.r.t. basic variables
      for (unsigned i=1; i<=x.size(); i++) {
        cout<<"The variable name is "<<vars_name[i]<<endl;
        if (vars_name[i].compare(0,2,"Em") == 0) {
          // w.r.t. Em
          dsdx.clear(); dsdx = StressGP_r_Em;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,3,"NUm") == 0) {
          // w.r.t. NUm
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUm;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,3,"NUp") == 0) {
          // w.r.t. NUp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUp;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,3,"NUz") == 0) {
          // w.r.t. NUpz
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUpz;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,2,"Ep") == 0) {
          // w.r.t. Ep
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ep;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,2,"Ez") == 0) {
          // w.r.t. Ez
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Ez;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,3,"Gzp") == 0) {
          // w.r.t. Gzp
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_Gzp;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,2,"Ef") == 0) {
          // w.r.t. Ef - fibre with isotropic material
          dsdx.clear(); dsdx = StressGP_r_Ef;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,3,"NUf") == 0) {
          // w.r.t. NUf - fibre with isotropic material
          dsdx.clear(); dsdx.resize(3,3); dsdx = StressGP_r_NUf;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,5,"force") == 0) {
          // w.r.t. Force
          dsdx.clear(); dsdx = StressGP_r_F;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,11,"orientation") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,6,"theta1") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_1;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,6,"theta2") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_2;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,6,"theta3") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_3;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,6,"theta4") == 0) {
          // w.r.t. Theta
          dsdx.clear(); dsdx = StressGP_r_Theta_4;
          grad_lsf(i-1) = - (  dsdx(0,0)*(1/X_T - 1/X_C)
                             + (dsdx(1,1) + dsdx(2,2))*(1/Y_T - 1/Y_C)
                             + 2*StressGP(0,0)*dsdx(0,0)/(X_T*X_C)
                             + 2*(StressGP(1,1) - StressGP(2,2))*(dsdx(1,1) - dsdx(2,2))/(Y_T*Y_C)
                             + 2*(StressGP(0,1)*dsdx(0,1) + StressGP(0,2)*dsdx(0,2))/(S_12*S_12)
                             + 2*StressGP(1,2)*dsdx(1,2)/(S_23*S_23)
                             - (dsdx(0,0)*StressGP(1,1) + StressGP(0,0)*dsdx(1,1) +
                                dsdx(0,0)*StressGP(2,2) + StressGP(0,0)*dsdx(2,2) -
                                dsdx(1,1)*StressGP(2,2) - StressGP(1,1)*dsdx(2,2))/(X_T*X_C));
        }
        else if (vars_name[i].compare(0,2,"XT") == 0) {
          // w.r.t. XT
          grad_lsf(i-1) =  - (- StressGP(0,0)/(X_T*X_T)
                              - StressGP(0,0)*StressGP(0,0)/(X_T*X_T*X_C)
                              + (StressGP(0,0)*StressGP(1,1) +
                                 StressGP(0,0)*StressGP(2,2) -
                                 StressGP(1,1)*StressGP(2,2))/(X_T*X_T*X_C));
        }
        else if (vars_name[i].compare(0,2,"XC") == 0) {
          // w.r.t. XC
          grad_lsf(i-1) =  - (  StressGP(0,0)/(X_C*X_C)
                              - StressGP(0,0)*StressGP(0,0)/(X_T*X_C*X_C)
                              + (StressGP(0,0)*StressGP(1,1) +
                                 StressGP(0,0)*StressGP(2,2) -
                                 StressGP(1,1)*StressGP(2,2))/(X_T*X_C*X_C));
        }
        else if (vars_name[i].compare(0,2,"YT") == 0) {
          // w.r.t. YT
          grad_lsf(i-1) = - (- (StressGP(1,1) + StressGP(2,2))/(Y_T*Y_T)
                             - pow((StressGP(1,1) - StressGP(2,2)),2)/(Y_T*Y_T*Y_C)
                             - 2*eta*Y_C*StressGP(1,2)*StressGP(1,2)/(S_23*S_23*S_23));
        }
        else if (vars_name[i].compare(0,2,"YC") == 0) {
          // w.r.t. YC
          grad_lsf(i-1) = - (  (StressGP(1,1) + StressGP(2,2))*( 1/(Y_C*Y_C))
                             - pow((StressGP(1,1) - StressGP(2,2)),2)/(Y_T*Y_C*Y_C)
                             - 2*eta*Y_T*StressGP(1,2)*StressGP(1,2)/(S_23*S_23*S_23));
        }
        else if (vars_name[i].compare(0,3,"S12") == 0) {
          // w.r.t. S12
          grad_lsf(i-1) = 2*StressGP(0,1)*StressGP(0,1)/(S_12*S_12*S_12);
        }
      }
      
      PetscFunctionReturn(0);
    }
    
  };
  
}

#endif //__LIMITSTATEFUNCTION_HPP__
