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

#ifndef __RELIABILITY_FUNCTIONS_HPP
#define __RELIABILITY_FUNCTIONS_HPP

namespace MoFEM {
  
  //======================================================
  
  //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  //                                                                            \\
  //              FUNCTIONS FOR STRUCTURAL RELIABILITY ANLYSIS                  \\
  //                                                                            \\
  //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  
  //------------------------------------------------------------------------------
  // Construct data structure <Stochastic_Model> for collecting data representing
  //   statistical information of inputs including probability distribution,
  //   correlation matrix and etc.
  
  struct Stochastic_Model {
    int    num_vars;                       // Number of variables
    double dist_type;                      // Distribution type index
    double transf_type;                    // Type of joint distribution
    double R0_method;                      // Method for computation of the modified Nataf correlation matrix
    int flag_sens;                         // Flag for computation of sensitivities w.r.t. parameters
    int ExaminedPly;
    vector<string> NameVars;               // Name of random variables
    string HomoMethod;                     // Homogenization method
    ublas::matrix<double> correlation;     // Correlation matrix
    ublas::matrix<double> marg;            // Marginal distribution for each random variable
    ublas::matrix<double> mod_correlation; // modified correlation matrix
    ublas::matrix<double> Lo;              // Chelosky decomposition
    ublas::matrix<double> inv_Lo;          // inverse of matrix Lo
    ublas::vector<double> MatStrength;     // Material strength
    ublas::vector<double> PlyAngle;        // Angle of orientation
  };
  
  //------------------------------------------------------------------------------
  // Construct data structure <Reliability_Options> to define calculation options
  //
  
  struct Reliability_Options {
    int echo_flag;            // Program interactive mode, 0: silent mode
    // FORM analysis options
    int istep_max;            // Maximum number of interations allowed in the search algorithm
    double e1;                // Tolerance on how close design point is to limit-state surface
    double e2;                // Tolerance on how accurately the gradient points towards the origin
    double step_code;         // 0: step size by Armijo rule, otherwise: given value is the step size
    int Recorded_u;           // 0: u-vector not recorded at all iterations, 1: u-vector recorded at all iterations
    int Recorded_x;           // 0: x-vector not recorded at all iterations, 1: x-vector recorded at all iterations
    int Recorded_beta;        // 0: beta not recorded at all iterations, 1: recorded at all iterations
    string grad_G;            // "PSFEM": perturbation based SFEM, "DDM": direct differentiation, 'ADM': automatic differentiation
  };
  
  //------------------------------------------------------------------------------
  // Construct data structure <LSF_Options> to define limit-state function options
  //
  
  struct LSF_Options {
    string evaluator; // Type of limit-state function evaluator,
    //   "basic", the LSF is defined by means of an analytical expression
    //   "failure", the LSF is defined by using failure criterion
    int flag_sens;    // Flag for sensitivity computation w.r.t. specified parameters of the LSF, 0: no, 1: yes
  };
  
  
  //------------------------------------------------------------------------------
  // To determine search direction
  //
  
  void search_dir(double val_G,
                  ublas::vector<double> grad_G,
                  ublas::vector<double> u,
                  ublas::vector<double> &u_dir) {
    // Determin direction cosine vector
    boost::numeric::ublas::vector<double> alpha;
    alpha = -grad_G/norm_2(grad_G);
    
    // Compute direction
    u_dir = ((val_G/norm_2(grad_G)) + inner_prod(alpha, u))*alpha - u;
  }
  
  //------------------------------------------------------------------------------
  // To generate sample
  //
  void REL_Stress_Transformation(double theta, ublas::matrix<double> Stress_xyz, ublas::matrix<double> &Stress_123) {
    
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
    
    ublas::vector<double> Stress_VectorNotation_xyz(6);
    Stress_VectorNotation_xyz(0) = Stress_xyz(0,0);
    Stress_VectorNotation_xyz(1) = Stress_xyz(1,1);
    Stress_VectorNotation_xyz(2) = Stress_xyz(2,2);
    Stress_VectorNotation_xyz(3) = Stress_xyz(0,1);
    Stress_VectorNotation_xyz(4) = Stress_xyz(1,2);
    Stress_VectorNotation_xyz(5) = Stress_xyz(2,0);
    
//    Stress_VectorNotation_xyz(3) = Stress_xyz(1,2);
//    Stress_VectorNotation_xyz(4) = Stress_xyz(0,2);
//    Stress_VectorNotation_xyz(5) = Stress_xyz(0,1);
    
    ublas::matrix<FieldData> T_strain_inv_T;    T_strain_inv_T.resize(6,6);
    T_strain_inv_T(0,0) =   l1*l1;
    T_strain_inv_T(0,1) =   m1*m1;
    T_strain_inv_T(0,2) =   n1*n1;
    T_strain_inv_T(0,3) = 2*l1*m1;
    T_strain_inv_T(0,4) = 2*m1*n1;
    T_strain_inv_T(0,5) = 2*l1*n1;
    
    T_strain_inv_T(1,0) =   l2*l2;
    T_strain_inv_T(1,1) =   m2*m2;
    T_strain_inv_T(1,2) =   n2*n2;
    T_strain_inv_T(1,3) = 2*l2*m2;
    T_strain_inv_T(1,4) = 2*m2*n2;
    T_strain_inv_T(1,5) = 2*l2*n2;
    
    T_strain_inv_T(2,0) =   l3*l3;
    T_strain_inv_T(2,1) =   m3*m3;
    T_strain_inv_T(2,2) =   n3*n3;
    T_strain_inv_T(2,3) = 2*l3*m3;
    T_strain_inv_T(2,4) = 2*m3*n3;
    T_strain_inv_T(2,5) = 2*l3*n3;
    
    T_strain_inv_T(3,0) = l1*l2;
    T_strain_inv_T(3,1) = m1*m2;
    T_strain_inv_T(3,2) = n1*n2;
    T_strain_inv_T(3,3) = l1*m2+m1*l2;
    T_strain_inv_T(3,4) = m1*n2+n1*m2;
    T_strain_inv_T(3,5) = l1*n2+n1*l2;
    
    T_strain_inv_T(4,0) = l2*l3;
    T_strain_inv_T(4,1) = m2*m3;
    T_strain_inv_T(4,2) = n2*n3;
    T_strain_inv_T(4,3) = l2*m3+m2*l3;
    T_strain_inv_T(4,4) = m2*n3+n2*m3;
    T_strain_inv_T(4,5) = l2*n3+n2*l3;
    
    T_strain_inv_T(5,0) = l1*l3;
    T_strain_inv_T(5,1) = m1*m3;
    T_strain_inv_T(5,2) = n1*n3;
    T_strain_inv_T(5,3) = l1*m3+m1*l3;
    T_strain_inv_T(5,4) = m1*n3+n1*m3;
    T_strain_inv_T(5,5) = l1*n3+n1*l3;
    
    ublas::vector<double> Stress_VectorNotation_123 = prod(T_strain_inv_T, Stress_VectorNotation_xyz);
    
    Stress_123.resize(3,3); Stress_123.clear();
    Stress_123(0,0) = Stress_VectorNotation_123(0);
    Stress_123(1,1) = Stress_VectorNotation_123(1);
    Stress_123(2,2) = Stress_VectorNotation_123(2);
    Stress_123(0,1) = Stress_123(1,0) = Stress_VectorNotation_123(3);
    Stress_123(1,2) = Stress_123(2,1) = Stress_VectorNotation_123(4);
    Stress_123(2,0) = Stress_123(0,2) = Stress_VectorNotation_123(5);

//    Stress_123(0,1) = Stress_123(1,0) = Stress_VectorNotation_123(5);
//    Stress_123(1,2) = Stress_123(2,1) = Stress_VectorNotation_123(3);
//    Stress_123(2,0) = Stress_123(0,2) = Stress_VectorNotation_123(4);
    
  }
  
  
  //------------------------------------------------------------------------------
  // Stress transformation
  //
  
  void REL_Stress_Transformation_Theta(double theta, ublas::matrix<double> Stress_xyz,
                                       ublas::matrix<double> Stress_xyz_r,
                                       ublas::matrix<double> &Stress_123_r) {
    
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
    
    ublas::vector<double> Stress_VectorNotation_xyz(6);
    Stress_VectorNotation_xyz(0) = Stress_xyz(0,0);
    Stress_VectorNotation_xyz(1) = Stress_xyz(1,1);
    Stress_VectorNotation_xyz(2) = Stress_xyz(2,2);
    Stress_VectorNotation_xyz(3) = Stress_xyz(0,1);
    Stress_VectorNotation_xyz(4) = Stress_xyz(1,2);
    Stress_VectorNotation_xyz(5) = Stress_xyz(2,0);
    
    ublas::vector<double> Stress_VectorNotation_xyz_r(6);
    Stress_VectorNotation_xyz_r(0) = Stress_xyz_r(0,0);
    Stress_VectorNotation_xyz_r(1) = Stress_xyz_r(1,1);
    Stress_VectorNotation_xyz_r(2) = Stress_xyz_r(2,2);
    Stress_VectorNotation_xyz_r(3) = Stress_xyz_r(0,1);
    Stress_VectorNotation_xyz_r(4) = Stress_xyz_r(1,2);
    Stress_VectorNotation_xyz_r(5) = Stress_xyz_r(2,0);
    
    ublas::matrix<FieldData> T_strain_inv_T;    T_strain_inv_T.resize(6,6);
    T_strain_inv_T(0,0) =   l1*l1;
    T_strain_inv_T(0,1) =   m1*m1;
    T_strain_inv_T(0,2) =   n1*n1;
    T_strain_inv_T(0,3) = 2*l1*m1;
    T_strain_inv_T(0,4) = 2*m1*n1;
    T_strain_inv_T(0,5) = 2*l1*n1;
    
    T_strain_inv_T(1,0) =   l2*l2;
    T_strain_inv_T(1,1) =   m2*m2;
    T_strain_inv_T(1,2) =   n2*n2;
    T_strain_inv_T(1,3) = 2*l2*m2;
    T_strain_inv_T(1,4) = 2*m2*n2;
    T_strain_inv_T(1,5) = 2*l2*n2;
    
    T_strain_inv_T(2,0) =   l3*l3;
    T_strain_inv_T(2,1) =   m3*m3;
    T_strain_inv_T(2,2) =   n3*n3;
    T_strain_inv_T(2,3) = 2*l3*m3;
    T_strain_inv_T(2,4) = 2*m3*n3;
    T_strain_inv_T(2,5) = 2*l3*n3;
    
    T_strain_inv_T(3,0) = l1*l2;
    T_strain_inv_T(3,1) = m1*m2;
    T_strain_inv_T(3,2) = n1*n2;
    T_strain_inv_T(3,3) = l1*m2+m1*l2;
    T_strain_inv_T(3,4) = m1*n2+n1*m2;
    T_strain_inv_T(3,5) = l1*n2+n1*l2;
    
    T_strain_inv_T(4,0) = l2*l3;
    T_strain_inv_T(4,1) = m2*m3;
    T_strain_inv_T(4,2) = n2*n3;
    T_strain_inv_T(4,3) = l2*m3+m2*l3;
    T_strain_inv_T(4,4) = m2*n3+n2*m3;
    T_strain_inv_T(4,5) = l2*n3+n2*l3;
    
    T_strain_inv_T(5,0) = l1*l3;
    T_strain_inv_T(5,1) = m1*m3;
    T_strain_inv_T(5,2) = n1*n3;
    T_strain_inv_T(5,3) = l1*m3+m1*l3;
    T_strain_inv_T(5,4) = m1*n3+n1*m3;
    T_strain_inv_T(5,5) = l1*n3+n1*l3;
    
    ublas::matrix<FieldData> T_strain_inv_T_r_Theta;    T_strain_inv_T_r_Theta.resize(6,6);
    T_strain_inv_T_r_Theta(0,0) =   2*l1*l1_r;
    T_strain_inv_T_r_Theta(0,1) =   2*m1*m1_r;
    T_strain_inv_T_r_Theta(0,2) =   2*n1*n1_r;
    T_strain_inv_T_r_Theta(0,3) = 2*(l1_r*m1 + l1*m1_r);
    T_strain_inv_T_r_Theta(0,4) = 2*(m1_r*n1 + m1*n1_r);
    T_strain_inv_T_r_Theta(0,5) = 2*(l1_r*n1 + l1*n1_r);
    
    T_strain_inv_T_r_Theta(1,0) =   2*l2*l2_r;
    T_strain_inv_T_r_Theta(1,1) =   2*m2*m2_r;
    T_strain_inv_T_r_Theta(1,2) =   2*n2*n2_r;
    T_strain_inv_T_r_Theta(1,3) = 2*(l2_r*m2 + l2*m2_r);
    T_strain_inv_T_r_Theta(1,4) = 2*(m2_r*n2 + m2*n2_r);
    T_strain_inv_T_r_Theta(1,5) = 2*(l2_r*n2 + l2*n2_r);
    
    T_strain_inv_T_r_Theta(2,0) =   2*l3*l3_r;
    T_strain_inv_T_r_Theta(2,1) =   2*m3*m3_r;
    T_strain_inv_T_r_Theta(2,2) =   2*n3*n3_r;
    T_strain_inv_T_r_Theta(2,3) = 2*(l3_r*m3 + l3*m3_r);
    T_strain_inv_T_r_Theta(2,4) = 2*(m3_r*n3 + m3*n3_r);
    T_strain_inv_T_r_Theta(2,5) = 2*(l3_r*n3 + l3*n3_r);
    
    T_strain_inv_T_r_Theta(3,0) = l1_r*l2 + l1*l2_r;
    T_strain_inv_T_r_Theta(3,1) = m1_r*m2 + m1*m2_r;
    T_strain_inv_T_r_Theta(3,2) = n1_r*n2 + n1*n2_r;
    T_strain_inv_T_r_Theta(3,3) = l1_r*m2 + l1*m2_r + m1_r*l2 + m1*l2_r;
    T_strain_inv_T_r_Theta(3,4) = m1_r*n2 + m1*n2_r + n1_r*m2 + n1*m2_r;
    T_strain_inv_T_r_Theta(3,5) = l1_r*n2 + l1*n2_r + n1_r*l2 + n1*l2_r;
    
    T_strain_inv_T_r_Theta(4,0) = l2_r*l3 + l2*l3_r;
    T_strain_inv_T_r_Theta(4,1) = m2_r*m3 + m2*m3_r;
    T_strain_inv_T_r_Theta(4,2) = n2_r*n3 + n2*n3_r;
    T_strain_inv_T_r_Theta(4,3) = l2_r*m3 + l2*m3_r + m2_r*l3 + m2*l3_r;
    T_strain_inv_T_r_Theta(4,4) = m2_r*n3 + m2*n3_r + n2_r*m3 + n2*m3_r;
    T_strain_inv_T_r_Theta(4,5) = l2_r*n3 + l2*n3_r + n2_r*l3 + n2*l3_r;
    
    T_strain_inv_T_r_Theta(5,0) = l1_r*l3 + l1*l3_r;
    T_strain_inv_T_r_Theta(5,1) = m1_r*m3 + m1*m3_r;
    T_strain_inv_T_r_Theta(5,2) = n1_r*n3 + n1*n3_r;
    T_strain_inv_T_r_Theta(5,3) = l1_r*m3 + l1*m3_r + m1_r*l3 + m1*l3_r;
    T_strain_inv_T_r_Theta(5,4) = m1_r*n3 + m1*n3_r + n1_r*m3 + n1*m3_r;
    T_strain_inv_T_r_Theta(5,5) = l1_r*n3 + l1*n3_r + n1_r*l3 + n1*l3_r;
    
    T_strain_inv_T_r_Theta = T_strain_inv_T_r_Theta*(M_PI/180);
    
    ublas::vector<double> StressVector_r_1 = prod(T_strain_inv_T, Stress_VectorNotation_xyz_r);
    ublas::vector<double> StressVector_r_2 = prod(T_strain_inv_T_r_Theta, Stress_VectorNotation_xyz);
    ublas::vector<double> Stress_VectorNotation_123_r = StressVector_r_1 + StressVector_r_2;
    
    Stress_123_r.resize(3,3); Stress_123_r.clear();
    Stress_123_r(0,0) = Stress_VectorNotation_123_r(0);
    Stress_123_r(1,1) = Stress_VectorNotation_123_r(1);
    Stress_123_r(2,2) = Stress_VectorNotation_123_r(2);
    Stress_123_r(0,1) = Stress_123_r(1,0) = Stress_VectorNotation_123_r(3);
    Stress_123_r(1,2) = Stress_123_r(2,1) = Stress_VectorNotation_123_r(4);
    Stress_123_r(2,0) = Stress_123_r(0,2) = Stress_VectorNotation_123_r(5);
    
  }
  
}

#endif //__RELIABILITY_FUNCTIONS_HPP