//
//  main.cpp
//  Reliability_Analysis
//
//  Created by Xiaoyi Zhou on 22/05/2015.
//  Copyright (c) 2015 Xiaoyi Zhou. All rights reserved.
//
// This code adopts the open-source Matlab toolbox - FERUM 4.1 (Finite Element
//   Reliability Using Matlab) devloped by Jean-Marc BOURINET at the IFMA (Institut
//   Français de Mécanique Avancée) in Clermont-Ferrand, France, which is a new
//   version of FERUM originally developed and maintained by Terje Haukaas,
//   Armen Der Kiureghian and other contributors at the University of California
//   at Berkeley initiated in 1999.
//
// FERUM 4.1 is available at http://www.ifma.fr/Recherche/laboratoires_recherche/FERUM
// FERUM 3.0 is available at http://www.ce.berkeley.edu/projects/ferum/
//
// References:
//   [1] Merchers, R. (1999) Structural reliability analysis and prediction,
//       2nd edition, Wiley.
//   [2] Ditlevsen, O. (2005) Structural reliability methods, electronic version
//   [3] Choi, et al. (2007) Reliability-based structural design, Springer
//   [4] Bourinet, J.-M. (2009) FERUM 4.0 User's Guide
//

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//              PROCEDURE OF SECOND-ORDER RELIABILITY METHOD                  //
//                                                                            //
//  1. Read FORM results including:                                           //
//       design point in x and u spaces,                                      //
//       values of LSF and its first-order derivative;                        //
//       cosine direction or alpha                                            //
//  2. Compute the second-order derivative of LSF                             //
//  3. Construct an orthogonal matrix P based on alpha                        //
//  4. Compute the matrix A = PHP'/||G'|| and form A11 by deleting the last   //
//       row and column of A                                                  //
//  5. Compute the eigenvalues of A11 to get curvatures, k_i, i = 1, ..., n-1 //
//  6. Calculate SORM approximation of faiure probability by empirical formula//
//       Breitung                                                             //
//       Tvedt                                                                //
//       Curve fitting [Der Kiureghian and De Stefano]                        //
//       Point fitting [Der Kiureghian et al.]                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

struct Stochastic_Model {
  int    num_vars;                       // Number of variables
  int    num_mat_vars;                   // Number of variables for material properties
  double dist_type;                      // Distribution type index
  double transf_type;                    // Type of joint distribution
  double R0_method;                      // Method for computation of the modified Nataf correlation matrix
  int flag_sens;                         // Flag for computation of sensitivities w.r.t. parameters
  int ExaminedPly;
  int AnalysisType;                      // Analysis type 10: FORM 20: SORM 30: MCS 31: MCIS
  vector<string> NameVars;               // Name of random variables
  vector<string> NameMatVars;            // Name of random variables for material properties
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