/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

#include <MoFEM.hpp>
using namespace MoFEM;

#include <DirichletBC.hpp>
#include <Projection10NodeCoordsOnField.hpp>

#include <petsctime.h>

#include <SurfacePressure.hpp>
#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcOnRefMesh.hpp>
#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

#include <ElasticFEMethod.hpp>
#include "ElasticFEMethodTransIso.hpp"

#include "ElasticFE_RVELagrange_Disp.hpp"
#include "ElasticFE_RVELagrange_Disp_Multi_Rhs.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Disp.hpp"
#include "RVEVolume.hpp"

using namespace ObosleteUsersModules;
#include "MaterialConstitutiveMatrix_FirstOrderDerivative.hpp"
#include "MaterialConstitutiveMatrix_SecondOrderDerivative.hpp"
#include "Trans_Iso_Rhs_r_PSFEM.hpp"
#include "Trans_Iso_Rhs_rs_PSFEM.hpp"

#include <FE2_ElasticFEMethod.hpp>

#include <FE2_Rhs_r_PSFEM.hpp>
#include <FE2_Rhs_rs_PSFEM.hpp>

#include <Reliability_SurfacePressure.hpp>

#include <FE2_Macro_Solver.hpp>

#include <FE2_PostProcStressForReliability.hpp>

using namespace boost::numeric;

//======================================================
// Declaration for reliability analysis
extern "C" {
#include <gm_rule.h>
#include <ltqnorm.h>
}

#include <iostream>
#include <fstream>
#include <ctime>
//#include <vector>
#include <new>
#include <ctype.h>

#include <boost/random.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions/gamma.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <cholesky.hpp>
#include <MatrixInverse.hpp> // download from http://proteowizard.sourceforge.net/dox/_matrix_inverse_8hpp.html

#include <ImportProbData.hpp>
#include <NatafTransformation.hpp>
#include <LimitStateFunction.hpp>
#include <Reliability_Input.hpp>
#include <SORM.hpp>

//======================================================

/*******************************************************************************
 *                                                                             *
 *             FUNCTIONS FOR STRUCTURAL RELIABILITY ANLYSIS                    *
 *                                                                             *
/*******************************************************************************/

//------------------------------------------------------------------------------
// Construct data structure <Stochastic_Model> for collecting data representing
//   statistical information of inputs including probability distribution,
//   correlation matrix and etc.

//struct Stochastic_Model {
//  int    num_vars;                       // Number of variables
//  int    num_mat_vars;                   // Number of variables for material properties
//  double dist_type;                      // Distribution type index
//  double transf_type;                    // Type of joint distribution
//  double R0_method;                      // Method for computation of the modified Nataf correlation matrix
//  int flag_sens;                         // Flag for computation of sensitivities w.r.t. parameters
//  int ExaminedPly;
//  int AnalysisType;                       // Analysis type 10: FORM 20: SORM 30: MCS 31: MCIS
//  vector<string> NameVars;               // Name of random variables
//  vector<string> NameMatVars;            // Name of random variables for material properties
//  ublas::matrix<double> correlation;     // Correlation matrix
//  ublas::matrix<double> marg;            // Marginal distribution for each random variable
//  ublas::matrix<double> mod_correlation; // modified correlation matrix
//  ublas::matrix<double> Lo;              // Chelosky decomposition
//  ublas::matrix<double> inv_Lo;          // inverse of matrix Lo
//  ublas::vector<double> MatStrength;     // Material strength
//  ublas::vector<double> PlyAngle;        // Angle of orientation
//};

//------------------------------------------------------------------------------
// Construct data structure <Reliability_Options> to define calculation options
//

//struct Reliability_Options {
//  int echo_flag;            // Program interactive mode, 0: silent mode
//  // FORM analysis options
//  int istep_max;            // Maximum number of interations allowed in the search algorithm
//  double e1;                // Tolerance on how close design point is to limit-state surface
//  double e2;                // Tolerance on how accurately the gradient points towards the origin
//  double step_code;         // 0: step size by Armijo rule, otherwise: given value is the step size
//  int Recorded_u;           // 0: u-vector not recorded at all iterations, 1: u-vector recorded at all iterations
//  int Recorded_x;           // 0: x-vector not recorded at all iterations, 1: x-vector recorded at all iterations
//  int Recorded_beta;        // 0: beta not recorded at all iterations, 1: recorded at all iterations
//  string grad_G;            // "PSFEM": perturbation based SFEM, "DDM": direct differentiation, 'ADM': automatic differentiation
//};

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


ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

//------------------------------------------------------------------------------
// To transform stress
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
  
}

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


int main(int argc, char *argv[]) {
  
  clock_t start_time, finish_time;
  double total_time;
  start_time = clock();
  
  const char* args_1st[] = {
    "_r_1",  "_r_2",  "_r_3",  "_r_4",  "_r_5",                                              // 1st order
    "_r_6",  "_r_7",  "_r_8",  "_r_9",  "_r_10",
    "_r_11", "_r_12", "_r_13", "_r_14", "_r_15"
  };
  const char* args_2nd[] = {
    "_rs_r1s1",
    "_rs_r1s2",  "_rs_r2s2",
    "_rs_r1s3",  "_rs_r2s3",  "_rs_r3s3",
    "_rs_r1s4",  "_rs_r2s4",  "_rs_r3s4",  "_rs_r4s4",
    "_rs_r1s5",  "_rs_r2s5",  "_rs_r3s5",  "_rs_r4s5",  "_rs_r5s5",
    "_rs_r1s6",  "_rs_r2s6",  "_rs_r3s6",  "_rs_r4s6",  "_rs_r5s6",  "_rs_r6s6",
    "_rs_r1s7",  "_rs_r2s7",  "_rs_r3s7",  "_rs_r4s7",  "_rs_r5s7",  "_rs_r6s7",  "_rs_r7s7",
    "_rs_r1s8",  "_rs_r2s8",  "_rs_r3s8",  "_rs_r4s8",  "_rs_r5s8",  "_rs_r6s8",  "_rs_r7s8",  "_rs_r8s8",
    "_rs_r1s9",  "_rs_r2s9",  "_rs_r3s9",  "_rs_r4s9",  "_rs_r5s9",  "_rs_r6s9",  "_rs_r7s9",  "_rs_r8s9",  "_rs_r9s9",
    "_rs_r1s10", "_rs_r2s10", "_rs_r3s10", "_rs_r4s10", "_rs_r5s10", "_rs_r6s10", "_rs_r7s10", "_rs_r8s10", "_rs_r9s10", "_rs_r10s10",
    "_rs_r1s11", "_rs_r2s11", "_rs_r3s11", "_rs_r4s11", "_rs_r5s11", "_rs_r6s11", "_rs_r7s11", "_rs_r8s11", "_rs_r9s11", "_rs_r10s11", "_rs_r11s11",
    "_rs_r1s12", "_rs_r2s12", "_rs_r3s12", "_rs_r4s12", "_rs_r5s12", "_rs_r6s12", "_rs_r7s12", "_rs_r8s12", "_rs_r9s12", "_rs_r10s12", "_rs_r11s12", "_rs_r12s12",
    "_rs_r1s13", "_rs_r2s13", "_rs_r3s13", "_rs_r4s13", "_rs_r5s13", "_rs_r6s13", "_rs_r7s13", "_rs_r8s13", "_rs_r9s13", "_rs_r10s13", "_rs_r11s13", "_rs_r12s13", "_rs_r13s13",
    "_rs_r1s14", "_rs_r2s14", "_rs_r3s14", "_rs_r4s14", "_rs_r5s14", "_rs_r6s14", "_rs_r7s14", "_rs_r8s14", "_rs_r9s14", "_rs_r10s14", "_rs_r11s14", "_rs_r12s13", "_rs_r13s14",  "_rs_r14s14",
    "_rs_r1s15", "_rs_r2s15", "_rs_r3s15", "_rs_r4s15", "_rs_r5s15", "_rs_r6s15", "_rs_r7s15", "_rs_r8s15", "_rs_r9s15", "_rs_r10s15", "_rs_r11s15", "_rs_r12s13", "_rs_r13s15",  "_rs_r14s15", "_rs_r15s15",
  }; // 2nd order
  
  //============================================================================
  //
  //  A. Micro (RVE) Problem
  //
  //============================================================================
  
  PetscInitialize(&argc,&argv,(char *)0,help);
  
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  moab::Core mb_instance_RVE;
  Interface& moab_RVE = mb_instance_RVE;
  
  
  // ===========================================================================
  //
  //  A.I. READ MESH DATA AND FEA COMPUTATION PARAMETERS FROM FILE
  //
  // ===========================================================================
  
  // ==
  Stochastic_Model probdata;
  
  // Import data from textile file
  ImportProbData readprobdata;
  
  ierr = readprobdata.ProbdataFileIn(); CHKERRQ(ierr);
  probdata.marg         = readprobdata.MargProb;
  probdata.correlation  = readprobdata.CorrMat;
  probdata.num_vars     = readprobdata.NumVars;
  probdata.MatStrength  = readprobdata.MatStrength;
  probdata.NameVars     = readprobdata.NameVars;
  probdata.PlyAngle     = readprobdata.PlyAngle;
  probdata.ExaminedPly  = readprobdata.ExaminedLayer;
  probdata.AnalysisType = readprobdata.AnalysisType;
  probdata.num_mat_vars = readprobdata.NumMatVars;
  probdata.NameMatVars  = readprobdata.NameMatVars;
  
  int nvars_rve_mat = 0;
  int nvars_ply_mat = 0;
  int nders_rve_mat = 0;
  int nders_ply_mat = 0;
  
  vector<string> name_rve_mat; name_rve_mat.clear();
  vector<string> name_ply_mat; name_ply_mat.clear();
  
  
  // ---------------------------------
  //
  // Identify random variables for
  //   microscale material properties
  //
  // ---------------------------------
  cout<<probdata.num_vars<<endl;
  for (int i = 1; i<=probdata.num_vars; i++) {
    // The variable names is
    cout<<"The variable is "<<probdata.NameVars[i]<<endl;
    if (probdata.NameVars[i].compare(0,2,"Em") == 0) {
      nvars_rve_mat++; nders_rve_mat = nders_rve_mat + nvars_rve_mat;
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_rve_mat.push_back(probdata.NameVars[i]);
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,3,"NUm") == 0) {
      nvars_rve_mat++; nders_rve_mat = nders_rve_mat + nvars_rve_mat;
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_rve_mat.push_back(probdata.NameVars[i]);
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,3,"NUz") == 0) {
      nvars_rve_mat++; nders_rve_mat = nders_rve_mat + nvars_rve_mat;
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_rve_mat.push_back(probdata.NameVars[i]);
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,4,"NUpz") == 0) {
      nvars_rve_mat++; nders_rve_mat = nders_rve_mat + nvars_rve_mat;
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_rve_mat.push_back(probdata.NameVars[i]);
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,2,"Ez") == 0) {
      nvars_rve_mat++; nders_rve_mat = nders_rve_mat + nvars_rve_mat;
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_rve_mat.push_back(probdata.NameVars[i]);
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,2,"Ep") == 0) {
      nvars_rve_mat++; nders_rve_mat = nders_rve_mat + nvars_rve_mat;
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_rve_mat.push_back(probdata.NameVars[i]);
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,3,"Gzp") == 0) {
      nvars_rve_mat++; nders_rve_mat = nders_rve_mat + nvars_rve_mat;
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_rve_mat.push_back(probdata.NameVars[i]);
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,2,"Ef") == 0) {
      nvars_rve_mat++; nders_rve_mat = nders_rve_mat + nvars_rve_mat;
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_rve_mat.push_back(probdata.NameVars[i]);
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,3,"NUf") == 0) {
      nvars_rve_mat++; nders_rve_mat = nders_rve_mat + nvars_rve_mat;
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_rve_mat.push_back(probdata.NameVars[i]);
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,11,"orientation") == 0) {
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,6,"theta1") == 0) {
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,6,"theta2") == 0) {
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,6,"theta3") == 0) {
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,6,"theta4") == 0) {
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
    else if (probdata.NameVars[i].compare(0,5,"force") == 0) {
      nvars_ply_mat++; nders_ply_mat = nders_ply_mat + nvars_ply_mat;
      name_ply_mat.push_back(probdata.NameVars[i]);
    }
  }
  
  nders_rve_mat = nders_rve_mat + nvars_rve_mat;
  nders_ply_mat = nders_ply_mat + nvars_ply_mat;
  
  cout<<"The number of variables for RVE material properties is "<<nvars_rve_mat<<"\t"<<nders_rve_mat<<endl;
  cout<<"The number of variables for ply material properties is "<<nvars_ply_mat<<"\t"<<nders_ply_mat<<endl;
  
  const char* args_rve[nders_rve_mat];
  for (int i = 0; i<nders_rve_mat; i++) {
    if (i<nvars_rve_mat) {
      args_rve[i] = args_1st[i];
    } else {
      args_rve[i] = args_2nd[i-nvars_rve_mat];
    }
    cout<<"The symbol is "<<args_rve[i]<<endl;
  }
  
  const char* args_ply[nders_ply_mat];
  for (int i = 0; i<nders_ply_mat; i++) {
    if (i<nvars_ply_mat) {
      args_ply[i] = args_1st[i];
    } else {
      args_ply[i] = args_2nd[i-nvars_ply_mat];
    }
    cout<<"The symbol is "<<args_ply[i]<<endl;
  }
  
//  int nvars = 9;    // number of variables
//  int nders = 54;   // number of partial derivatives (first- and second- order)
//  nvars = nvars_rve_mat;
//  nders = nders_rve_mat;
//  vector<string> stochastic_fields(args_rve, args_rve + nders_rve_mat);
  
  vector<string> stochastic_fields_rve(args_rve, args_rve + nders_rve_mat);
  vector<string> stochastic_fields_ply(args_ply, args_ply + nders_ply_mat);
  
  /*****************************************************************************
   *
   * Read parameters from line command
   *
   ****************************************************************************/
  
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_RVE",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_RVE (MESH FILE NEEDED)");
  }
  
  PetscInt order_RVE;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order_RVE",&order_RVE,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order_RVE = 1;
  }
cout<<"stop 1"<<endl;
  
  /*****************************************************************************
   *
   * Transfer mesh data to MOAB database
   *
   ****************************************************************************/
  
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab_RVE.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab_RVE,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab_RVE,PETSC_COMM_WORLD);
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core_RVE(moab_RVE);
  FieldInterface& m_field_RVE = core_RVE;
  
 cout<<"stop 2"<<endl;
  /*****************************************************************************
   *
   * Get fibre direction information from the potential-flow calculation
   *
   ****************************************************************************/
  
  Tag th_phi;
  //    double def_val  = 0;
  rval = moab_RVE.tag_get_handle("PHI",1,MB_TYPE_DOUBLE,th_phi); CHKERR_PETSC(rval);
  
  Tag th_meshset_info;
  int def_meshset_info[2] = {0,0};
  rval = moab_RVE.tag_get_handle("MESHSET_INFO",2,MB_TYPE_INTEGER,th_meshset_info,MB_TAG_CREAT|MB_TAG_SPARSE,&def_meshset_info);
  
  int meshset_data[2];
  EntityHandle root = moab_RVE.get_root_set();
  rval = moab_RVE.tag_get_data(th_meshset_info,&root,1,meshset_data); CHKERR_PETSC(rval);
  
  vector<BitRefLevel> bit_levels;
  bit_levels.push_back(BitRefLevel().set(meshset_data[0]-1));
  
  ierr = m_field_RVE.build_fields(); CHKERRQ(ierr);
  ierr = m_field_RVE.build_finite_elements(); CHKERRQ(ierr);
  ierr = m_field_RVE.build_adjacencies(bit_levels.back()); CHKERRQ(ierr);
  ierr = m_field_RVE.build_problems(); CHKERRQ(ierr);
cout<<"stop 3"<<endl;
  /*****************************************************************************
   *
   * Group element into various mesh-set
   * meshset_level0: all element
   * meshset_Matrix: element set for matrix
   * meshset_Inclusion: element set for inclusion/fibre
   *
   ****************************************************************************/
  
  EntityHandle out_meshset;
  rval = moab_RVE.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
  //    ierr = m_field_RVE.get_problem_finite_elements_entities("POTENTIAL_PROBLEM","POTENTIAL_ELEM",out_meshset); CHKERRQ(ierr);
  ierr = m_field_RVE.get_entities_by_ref_level(bit_levels.back(),BitRefLevel().set(),out_meshset); CHKERRQ(ierr);
  Range LatestRefinedTets;
  rval = moab_RVE.get_entities_by_type(out_meshset, MBTET,LatestRefinedTets,true); CHKERR_PETSC(rval);
  
  Range LatestRefinedPrisms;
  rval = moab_RVE.get_entities_by_type(out_meshset, MBPRISM,LatestRefinedPrisms,true); CHKERR_PETSC(rval);
  
  cout<<"No of Prisms/Interfaces = "<<LatestRefinedPrisms.size()<<endl;
  
  BitRefLevel problem_bit_level_RVE = bit_levels.back();
  
  EntityHandle meshset_Matrix, meshset_Fibre;
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_Matrix); CHKERR_PETSC(rval);
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_Fibre); CHKERR_PETSC(rval);
  
  ///Getting No. of Fibres to be used for Potential Flow Problem
  int noOfFibres=0;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET|UNKNOWNCUBITNAME,it)) {
    
    std::size_t found=it->get_name().find("PotentialFlow");
    if (found==std::string::npos) continue;
    noOfFibres += 1;
  }
  cout<<"No. of Fibres for Potential Flow : "<<noOfFibres<<endl;
  
  //  vector<int> fibreList(noOfFibres,0);
  //  for (int aa=0; aa<noOfFibres; aa++) {
  //    fibreList[aa] = aa + 1;
  //  }
  
  vector<Range> RangeFibre(noOfFibres);
  vector<EntityHandle> fibre_meshset(noOfFibres);
  
  for (int ii=0; ii<noOfFibres; ii++) {
    ostringstream sss;
    sss << "POTENTIAL_ELEM" << ii+1;
    for(_IT_GET_FES_BY_NAME_FOR_LOOP_(m_field_RVE, sss.str().c_str() ,it)){
      RangeFibre[ii].insert(it->get_ent());
      rval = moab_RVE.create_meshset(MESHSET_SET,fibre_meshset[ii]); CHKERR_PETSC(rval);
      rval = moab_RVE.add_entities(fibre_meshset[ii],RangeFibre[ii]); CHKERR_PETSC(rval);
      rval = moab_RVE.unite_meshset(meshset_Fibre,fibre_meshset[ii]); CHKERR_PETSC(rval);
    }
  }
  
  //rval = moab_RVE.write_file("meshset_Fibre.vtk","VTK","",&meshset_Fibre,1); CHKERR_PETSC(rval);
  
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)){
    
    if(it->get_name() == "MAT_ELASTIC_1") {
      Range TetsInBlock;
      rval = moab_RVE.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
      Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);
      
      cout<<"=============  TetsInBlock  "<< TetsInBlock.size() <<endl;
      
      rval = moab_RVE.add_entities(meshset_Matrix,block_rope_bit_level);CHKERR_PETSC(rval);
      
    }
  }
  ierr = m_field_RVE.seed_finite_elements(meshset_Matrix); CHKERRQ(ierr);
  
  Range prims_on_problem_bit_level_RVE;
  ierr = m_field_RVE.get_entities_by_type_and_ref_level(problem_bit_level_RVE,BitRefLevel().set(),MBPRISM,prims_on_problem_bit_level_RVE); CHKERRQ(ierr);
  //to create meshset from range
  EntityHandle meshset_prims_on_problem_bit_level_RVE;
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_prims_on_problem_bit_level_RVE); CHKERR_PETSC(rval);
  rval = moab_RVE.add_entities(meshset_prims_on_problem_bit_level_RVE,prims_on_problem_bit_level_RVE); CHKERR_PETSC(rval);
  ierr = m_field_RVE.seed_ref_level_MESHSET(meshset_prims_on_problem_bit_level_RVE,BitRefLevel().set()); CHKERRQ(ierr);

  
  // ===========================================================================
  //
  // A.II. DEFINE PROBLEM
  //
  // ===========================================================================
  
  /*****************************************************************************
   *
   * Add field
   *  (1) Deterministic fields
   *  (2) Stochastic fields
   *       (total 14 field for 1st and 2nd order stochastic PSFEM)
   *
   ****************************************************************************/
  // Deterministic fields
  int field_rank=3;
  ierr = m_field_RVE.add_field("DISP_RVE",H1,field_rank); CHKERRQ(ierr);
  ierr = m_field_RVE.add_field("Lagrange_mul_disp",H1,field_rank); CHKERRQ(ierr);
  
  // Stochastic fields for perturbation method
  for(int ii=0; ii < nders_rve_mat; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields_rve[ii];
    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_RVE.add_field(ss_field.str().c_str(),H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  }
  
  
  /*****************************************************************************
   *
   * Create finite element for the defined fields
   *
   ****************************************************************************/
  ierr = m_field_RVE.add_finite_element("ELASTIC_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.add_finite_element("TRAN_ISO_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.add_finite_element("Lagrange_FE"); CHKERRQ(ierr);
  
  
  /*****************************************************************************
   *
   * set field data which finite element use
   * set field row which finite element use
   *
   * [K][U] = [F]                                           Zeroth-order problem
   * [K][D_r U] = - [D_r K][U]                               First-order problem
   * [K][H_rs U] = - [H_rs K][U] - 2[D_r K][D_s U]          Second-order problem
   *
   ****************************************************************************/
  //Define rows/cols and element data
  ierr = m_field_RVE.modify_finite_element_add_field_row("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  
  // Stochastic
  //for(int ii=0; ii < nvars; ii++ ) {
  for(int ii=0; ii < nders_rve_mat; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields_rve[ii];
    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_RVE.modify_finite_element_add_field_data("ELASTIC_FE_RVE",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  //FE Transverse Isotropic
  ierr = m_field_RVE.modify_finite_element_add_field_row("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE","POTENTIAL_FIELD"); CHKERRQ(ierr);
  
  
  //for(int ii=0; ii < nvars; ii++ ) {
  for(int ii=0; ii < nders_rve_mat; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields_rve[ii];
    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  //C and CT
  //======================================================================================================
  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = m_field_RVE.modify_finite_element_add_field_row("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  
  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = m_field_RVE.modify_finite_element_add_field_col("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_row("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  
  //data
  ierr = m_field_RVE.modify_finite_element_add_field_data("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  //======================================================================================================
  
  
  //define problems
  ierr = m_field_RVE.add_problem("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  
  
  //set finite elements for problem
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","Lagrange_FE"); CHKERRQ(ierr);
  
  
  //set refinment level for problem
  ierr = m_field_RVE.modify_problem_ref_level_add_bit("ELASTIC_PROBLEM_RVE",problem_bit_level_RVE); CHKERRQ(ierr); // problem_bit_level
  
  
  // ===========================================================================
  //
  // A.III. DECLARE PROBLEM
  //
  // ===========================================================================
  
  /*****************************************************************************
   *
   * Add entitities (by tets) to the field
   *
   ****************************************************************************/
  ierr = m_field_RVE.add_ents_to_field_by_TETs(0,"DISP_RVE"); CHKERRQ(ierr);
  for(int ii=0; ii < nders_rve_mat; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields_rve[ii];
    ierr = m_field_RVE.add_ents_to_field_by_TETs(0,ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  
  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  ierr = m_field_RVE.add_ents_to_finite_element_by_TETs(meshset_Matrix,"ELASTIC_FE_RVE",true); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_finite_element_by_TETs(meshset_Fibre,"TRAN_ISO_FE_RVE",true); CHKERRQ(ierr);
  
  Range SurfacesFaces;
  ierr = m_field_RVE.get_cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_FE"); CHKERRQ(ierr);
  
  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab_RVE.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
  rval = moab_RVE.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = m_field_RVE.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Lagrange_mul_disp",2); CHKERRQ(ierr);
  
  /*****************************************************************************
   *
   * Set applied order
   * See reference for detals:
   *   Ainsworth M. and Coyle J. (2003) Hierarchic finite element bases on
   unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/
  ierr = m_field_RVE.set_field_order(0,MBTET,   "DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBTRI,   "DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,  "DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"DISP_RVE",1); CHKERRQ(ierr);
  
  ierr = m_field_RVE.set_field_order(0,MBTRI,   "Lagrange_mul_disp",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,  "Lagrange_mul_disp",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);
  
  int order_st = order_RVE;
  for(int ii=0; ii < nders_rve_mat; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields_rve[ii];
    //    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_RVE.set_field_order(0,MBTET,   ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_RVE.set_field_order(0,MBTRI,   ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_RVE.set_field_order(0,MBEDGE,  ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_RVE.set_field_order(0,MBVERTEX,ss_field.str().c_str(),1);        CHKERRQ(ierr);
  }
  
  // ===========================================================================
  //
  //  A.IV. BUILD DATABASE
  //
  // ===========================================================================
  
  //build field
  ierr = m_field_RVE.build_fields(); CHKERRQ(ierr);
  
  //build finite elemnts
  ierr = m_field_RVE.build_finite_elements(); CHKERRQ(ierr);
  
  //build adjacencies
  ierr = m_field_RVE.build_adjacencies(problem_bit_level_RVE); CHKERRQ(ierr);
  
  //build problem
  ierr = m_field_RVE.build_problems(); CHKERRQ(ierr);
  
  // ===========================================================================
  //
  //  A.V. MESH PARTITION
  //
  // ===========================================================================
  
  //partition
  ierr = m_field_RVE.partition_problem("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.partition_finite_elements("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field_RVE.partition_ghost_dofs("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  
  //print bcs
  ierr = m_field_RVE.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field_RVE.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = m_field_RVE.print_cubit_materials_set(); CHKERRQ(ierr);
  
  
  
  //============================================================================
  //
  //  B. Macro Problem
  //
  //============================================================================
  
  // ===========================================================================
  //
  // B. I. READ MESH DATA AND FINITE ELEMENT ANALYSIS CONTROL PARAMETERS FROM FILE
  //
  // ===========================================================================
  moab::Core mb_instance_Macro;
  Interface& moab_Macro = mb_instance_Macro;
  
  //ParallelComm* pcomm_Macro = ParallelComm::get_pcomm(&moab_Macro,MYPCOMM_INDEX);
  //if(pcomm == NULL) pcomm_Macro =  new ParallelComm(&moab_Macro,PETSC_COMM_WORLD);
  
  /*****************************************************************************
   *
   * Read parameters from line command
   *
   ****************************************************************************/
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_macro",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_macro (MESH FILE NEEDED)");
  }
  
  PetscInt order_Macro;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order_Macro",&order_Macro,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order_Macro = 1;
  }
  
  PetscInt NO_Layers;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_NO_Layers",&NO_Layers,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    NO_Layers = 1;
  }
  cout<<"\n\nNumber of layers: "<<NO_Layers<<endl;
  
  /*****************************************************************************
   *
   * Transfer mesh data to MOAB database
   *
   ****************************************************************************/
  rval = moab_Macro.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core_Macro(moab_Macro);
  FieldInterface& m_field_Macro = core_Macro;

  //ref meshset ref level 0
  ierr = m_field_Macro.seed_ref_level_3D(0,0); CHKERRQ(ierr);

  // stl::bitset see for more details
  BitRefLevel bit_level0_Macro;
  bit_level0_Macro.set(0);

  EntityHandle meshset_level0;
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field_Macro.seed_ref_level_3D(0,bit_level0_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.get_entities_by_ref_level(bit_level0_Macro,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  
  Range TetsInBlock_1st_Ply, TetsInBlock_2nd_Ply, TetsInBlock_3rd_Ply, TetsInBlock_4th_Ply, TetsInBlock_Reliability;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_Macro,BLOCKSET,it)){
    if(it->get_name() == "MAT_ELASTIC_First") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_1st_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_Second") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_2nd_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_Third") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_3rd_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_Fourth") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_4th_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "RELIABILITY") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_Reliability,true); CHKERR_PETSC(rval);
    }
  }
  
  
  // ===========================================================================
  //
  // B.II. DEFINE PROBLEM
  //
  // ===========================================================================
  

  /*****************************************************************************
   *
   * Add stochastic field
   * (total 14 field for 1st and 2nd order stochastic PSFEM)
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_field("DISP_MACRO",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field_Macro.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  
  // Stochastic fields for perturbation methods at macroscale
  for(int ii=0; ii < nders_ply_mat; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields_ply[ii];
    ierr = m_field_Macro.add_field(ss_field.str().c_str(),H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  }
  
  /*****************************************************************************
   *
   * Create finite element for the defined fields
   *
   ****************************************************************************/
  // First layer
  ierr = m_field_Macro.add_finite_element("ELASTIC_1st_Ply"); CHKERRQ(ierr);
  // Second layer
  if (NO_Layers > 1) {
    ierr = m_field_Macro.add_finite_element("ELASTIC_2nd_Ply"); CHKERRQ(ierr);
  }
  // Third layer
  if (NO_Layers > 2) {
    ierr = m_field_Macro.add_finite_element("ELASTIC_3rd_Ply"); CHKERRQ(ierr);
  }
  // Fourth layer
  if (NO_Layers > 3) {
    ierr = m_field_Macro.add_finite_element("ELASTIC_4th_Ply"); CHKERRQ(ierr);
  }
  
  ierr = m_field_Macro.add_finite_element("ELASTIC_FE_MACRO_REL"); CHKERRQ(ierr);
 

  /*****************************************************************************
   *
   * set field data which finite element use
   *
   ****************************************************************************/
  // First layer
  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_1st_Ply","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_1st_Ply","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  for(int ii=0; ii < nders_ply_mat; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields_ply[ii];
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  // Second layer
  if (NO_Layers > 1) {
    //Define rows/cols and element data
    ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_2nd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_2nd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    
    for(int ii=0; ii < nders_ply_mat; ii++ ) {
      ostringstream ss_field;
      ss_field << "DISP_MACRO" << stochastic_fields_ply[ii];
      ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply",ss_field.str().c_str()); CHKERRQ(ierr);
    }
  }
  
  // Third layer
  if (NO_Layers > 2) {
    //Define rows/cols and element data
    ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_3rd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_3rd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    
    for(int ii=0; ii < nders_ply_mat; ii++ ) {
      ostringstream ss_field;
      ss_field << "DISP_MACRO" << stochastic_fields_ply[ii];
      ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply",ss_field.str().c_str()); CHKERRQ(ierr);
    }
  }
  
  // Fourth layer
  if (NO_Layers > 3) {
    //Define rows/cols and element data
    ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_4th_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_4th_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    
    for(int ii=0; ii < nders_ply_mat; ii++ ) {
      ostringstream ss_field;
      ss_field << "DISP_MACRO" << stochastic_fields_ply[ii];
      ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply",ss_field.str().c_str()); CHKERRQ(ierr);
    }
  }

  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_FE_MACRO_REL","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_FE_MACRO_REL","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  for(int ii=0; ii < nders_ply_mat; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields_ply[ii];
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  
  //define problems
  ierr = m_field_Macro.add_problem("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  
  //set finite elements for problem
  // First layer
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply"); CHKERRQ(ierr);
  // Second layer
  if (NO_Layers > 1) {
    ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply"); CHKERRQ(ierr);
  }
  // Third layer
  if (NO_Layers > 2) {
    ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply"); CHKERRQ(ierr);
  }
  // Fourth layer
  if (NO_Layers > 3) {
    ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply"); CHKERRQ(ierr);
  }
  
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL"); CHKERRQ(ierr);
  
  //set refinment level for problem
  ierr = m_field_Macro.modify_problem_ref_level_add_bit("ELASTIC_PROBLEM_MACRO",bit_level0_Macro); CHKERRQ(ierr);
  
  // ===========================================================================
  //
  // B.III. DECLARE PROBLEM
  //
  // ===========================================================================

  /*****************************************************************************
   *
   * Add entitities (by tets) to the field
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  for(int ii=0; ii < nders_ply_mat; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields_ply[ii];
    ierr = m_field_Macro.add_ents_to_field_by_TETs(0,ss_field.str().c_str()); CHKERRQ(ierr);
  }

  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  // First layer
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_1st_Ply, "ELASTIC_1st_Ply"); CHKERRQ(ierr);
  // Second layer
  if (NO_Layers > 1) {
    ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_2nd_Ply,"ELASTIC_2nd_Ply"); CHKERRQ(ierr);
  }
  // Third layer
  if (NO_Layers > 2) {
    ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_3rd_Ply,"ELASTIC_3rd_Ply"); CHKERRQ(ierr);
  }
  // Fourth layer
  if (NO_Layers > 3) {
    ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_4th_Ply,"ELASTIC_4th_Ply"); CHKERRQ(ierr);
  }
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_Reliability,"ELASTIC_FE_MACRO_REL"); CHKERRQ(ierr);

  /*****************************************************************************
   *
   * Set applied order
   * See reference for detals:
   *   Ainsworth M. and Coyle J. (2003) Hierarchic finite element bases on
   *      unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/
  ierr = m_field_Macro.set_field_order(0,MBTET,   "DISP_MACRO",order_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,   "DISP_MACRO",order_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,  "DISP_MACRO",order_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"DISP_MACRO",1); CHKERRQ(ierr);
  
  order_st = order_Macro;
  for(int ii=0; ii < nders_ply_mat; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields_ply[ii];
    ierr = m_field_Macro.set_field_order(0,MBTET,   ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBTRI,   ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBEDGE,  ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBVERTEX,ss_field.str().c_str(),1); CHKERRQ(ierr);
  }
  
  //
  ierr = m_field_Macro.set_field_order(0,MBTET,   "MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,   "MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,  "MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
 
 //***************************************************************************
  ierr = MetaNeummanForces::addNeumannBCElements(m_field_Macro,"DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","FORCE_FE"); CHKERRQ(ierr);
  

  // ===========================================================================
  //
  //  B.IV. BUILD DATABASE
  //
  // ===========================================================================
  
  //build field
  ierr = m_field_Macro.build_fields(); CHKERRQ(ierr);
  
  //build finite elemnts
  ierr = m_field_Macro.build_finite_elements(); CHKERRQ(ierr);
  
  //build adjacencies
  ierr = m_field_Macro.build_adjacencies(bit_level0_Macro); CHKERRQ(ierr);
  
  //build problem
  ierr = m_field_Macro.build_problems(); CHKERRQ(ierr);
  

  // ===========================================================================
  //
  //  B.V. MESH PARTITION
  //
  // ===========================================================================
  
  //partition
  ierr = m_field_Macro.partition_problem("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.partition_finite_elements("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.partition_ghost_dofs("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  
  //m_field_Macro.list_dofs_by_field_name("DISP_MACORO",true);
  
  //print bcs
  ierr = m_field_Macro.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field_Macro.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = m_field_Macro.print_cubit_materials_set(); CHKERRQ(ierr);

 
  // ===========================================================================
  //
  //  C. RELIABILITY ANLYSIS
  //
  // ===========================================================================
  
  cout<<"\n\n";
  cout<<"///////////////////////////////////////////////////////////////////\n";
  cout<<"//                                                               //\n";
  cout<<"//           Reliability calculation starts from here!           //\n";
  cout<<"//                                                               //\n";
  cout<<"/////////////////////////////////////////////////////////////////\n\n";
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 1: read inputs from data file                                //
  //                a) probability data                                       //
  //                 .1) marginal distribution type                           //
  //                 .2) parameters                                           //
  //                 .3) correlation matrix                                   //
  //                b) limit state function                                   //
  //                 .1) indicator of how limit-state function is given       //
  //                     1:                                                   //
  //                c) analysis option                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  double beta;                    // Reliability index
  
  /*
   *  Set reliability calculation options
   */
  Reliability_Options ReliabOpt;
  
  ReliabOpt.echo_flag     = 1;       // Program interactive mode, 0: silent mode
  // For FORM
  ReliabOpt.istep_max     = 2000;    // Maximum number of interation allowed in the search algorithm
  ReliabOpt.e1            = 0.001;   // Tolerance on how close design point is to limit-state surface
  ReliabOpt.e2            = 0.001;   // Tolerance on how accurately the gradient points towards the origin
  ReliabOpt.step_code     = 1.0;       // 0: step size by Armijo rule, otherwise: given value is the step size
  ReliabOpt.Recorded_u    = 1;       // 0: u-vector not recorded at all iterations, 1: u-vector recorded at all iterations
  ReliabOpt.Recorded_x    = 1;       // 0: x-vector not recorded at all iterations, 1: x-vector recorded at all iterations
  ReliabOpt.Recorded_beta = 1;
  ReliabOpt.grad_G        = "PSFEM"; // "PSFEM": perturbation, "DDM": direct differentiation, 'ADM': automatic differentiation
  
  int    echo_flag = ReliabOpt.echo_flag;
  double e1        = ReliabOpt.e1;
  double e2        = ReliabOpt.e2;
  int    istep_max = ReliabOpt.istep_max;
  double step_code = ReliabOpt.step_code;
  int    beta_flag = ReliabOpt.Recorded_beta;
  
  /*
   *  Read inputs' statistical properties from file to insert into <probdata>
   */
  /*Stochastic_Model probdata;
  
  // Import data from textile file
  // ImportProbData readprobdata;
  
  ierr = readprobdata.ProbdataFileIn(); CHKERRQ(ierr);
  probdata.marg         = readprobdata.MargProb;
  probdata.correlation  = readprobdata.CorrMat;
  probdata.num_vars     = readprobdata.NumVars;
  probdata.MatStrength  = readprobdata.MatStrength;
  probdata.NameVars     = readprobdata.NameVars;
  probdata.PlyAngle     = readprobdata.PlyAngle;
  probdata.ExaminedPly  = readprobdata.ExaminedLayer;
  probdata.AnalysisType = readprobdata.AnalysisType;
  probdata.num_mat_vars = readprobdata.NumMatVars;
  probdata.NameMatVars  = readprobdata.NameMatVars;
   */
  
  int FailureCriterion; // 1: Tsai-Wu, 2: Tsai-Hill
  string NameOfFailureCriterion;
  FailureCriterion = readprobdata.FailureCriterion;
  
  step_code = readprobdata.SearchStep;
  
  string var_name;
  
  int PSFE_order = 1;
  
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

  /*
   * Perform Cholesky decomposition for the modified correlation matrix
   *    A = LL'
   */
  ublas::triangular_matrix<double, ublas::lower> Lmat(probdata.num_vars,probdata.num_vars);
  cholesky_decompose(probdata.mod_correlation,Lmat);
  probdata.Lo.resize(probdata.num_vars,probdata.num_vars);
  probdata.Lo = Lmat;
  
  //cout<<"Cholesky decomposed matrix: "<<Lmat<<endl;
  
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
  
  // ===========================================================================
  //
  //  C. SOLVING FE EQUATION
  //
  // ===========================================================================
  // Declare matrix for stroring RVE constitutive matrix & its derivatives
  ublas::matrix<double> Dmat;
  ublas::vector<ublas::matrix<double> > Dmat_r;
  ublas::matrix<ublas::matrix<double> > Dmat_rs;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 4: Perform iterative loop to find design point               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  /*
   *  Set parameters for the iterative loop
   */
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
  
  ublas::matrix<double> Hess_g;                 // Hessian matrix of LSF in x space
  grad_g.resize(probdata.num_vars); grad_g.clear();
  Hess_g.resize(probdata.num_vars,probdata.num_vars); Hess_g.clear();
  
  //
  ublas::matrix<double> StressGP_Global(3,3);    StressGP_Global.clear();
  ublas::matrix<double> StressGP(3,3);           StressGP.clear();
  ublas::vector<ublas::matrix<double> > StressGP_r(nvars_ply_mat);
  ublas::matrix<ublas::matrix<double> > StressGP_rs(nvars_ply_mat,nvars_ply_mat);
  
  for (int ivar = 0; ivar<nvars_ply_mat; ivar++) {
    StressGP_r(ivar).resize(3,3); StressGP_r(ivar).clear();
    for (int jvar = 0; jvar<nvars_ply_mat; jvar++) {
      StressGP_rs(ivar,jvar).resize(3,3); StressGP_rs(ivar,jvar).clear();
    }
  }
  
  double theta_angle;
  ublas::vector<double> PlyAngle_new;
  PlyAngle_new = probdata.PlyAngle;
  
  /*
   *  Start iteration
   */
  
  clock_t rel_t1, rel_t2;
  double rel_calc_time;
  rel_t1 = clock();
  
  FE2_Macro_Solver_Laminate Solve_FE2_Problem;
  LSF_Composite_Lamina TheLSF;
  //cout<<beta_flag<<endl;
  
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
    for (int i=1; i<=x.size();i++) {
      if (probdata.NameVars[i].compare(0,11,"orientation") == 0) {
        cout<<"\nAngle "<<x(i-1)<<endl;
        for (int j=0; j<probdata.PlyAngle.size(); j++) {
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
    // Micro-scale FE
    ierr = Solve_FE2_Problem.RVE_Dmat_Disp(m_field_RVE,
                                           stochastic_fields_rve,
                                           x,
                                           probdata.num_vars,
                                           probdata.NameVars,
                                           PSFE_order); CHKERRQ(ierr);
    // Macro-scale FE
    ierr = Solve_FE2_Problem.Macro_FE_REL_FSORM(m_field_Macro,
                                                nvars_rve_mat,
                                                nvars_ply_mat,
                                                name_rve_mat,
                                                name_ply_mat,
                                                stochastic_fields_ply,
                                                x,PlyAngle_new,//probdata.PlyAngle,
                                                NO_Layers,
                                                PSFE_order); CHKERRQ(ierr);
    
    Dmat.resize(6,6);
    Dmat_r.resize(nvars_ply_mat);
    Dmat_rs.resize(nvars_ply_mat, nvars_ply_mat);
    
    switch (probdata.ExaminedPly) {
      case 1:{cout<<"\n\nThe 1st layer is under examination.\n\n";
        //--------------------------
        // D matrix for the 1st ply
        //--------------------------
        
        for (int i=0; i<nvars_ply_mat; i++) {
          Dmat_r(i) = Solve_FE2_Problem.Ply_1st_Dmat_r(i);
          for (int j=0; j<nvars_ply_mat; j++) {
            Dmat_rs(i,j) = Solve_FE2_Problem.Ply_1st_Dmat_rs(i,j);
          }
        }
        Dmat = Solve_FE2_Problem.Dmat_1st_Ply; //cout<<"Dmat: "<<Dmat<<"\n\n";
        break;
      }
      case 2: {cout<<"\n\nThe 2nd layer is under examination.\n\n";
        //--------------------------
        // D matrix for the 2nd ply
        //--------------------------
        for (int i=0; i<nvars_ply_mat; i++) {
          Dmat_r(i) = Solve_FE2_Problem.Ply_2nd_Dmat_r(i);
          for (int j=0; j<nvars_ply_mat; j++) {
            Dmat_rs(i,j) = Solve_FE2_Problem.Ply_2nd_Dmat_rs(i,j);
          }
        }
        
        Dmat = Solve_FE2_Problem.Dmat_2nd_Ply; //cout<<"Dmat: "<<Dmat<<"\n\n";
        break;
      }
      case 3: {cout<<"\n\nThe 3rd layer is under examination.\n\n";
        //--------------------------
        // D matrix for the 3rd ply
        //--------------------------
        for (int i=0; i<nvars_ply_mat; i++) {
          Dmat_r(i) = Solve_FE2_Problem.Ply_3rd_Dmat_r(i);
          for (int j=0; j<nvars_ply_mat; j++) {
            Dmat_rs(i,j) = Solve_FE2_Problem.Ply_3rd_Dmat_rs(i,j);
          }
        }
        
        Dmat = Solve_FE2_Problem.Dmat_3rd_Ply; //cout<<"Dmat: "<<Dmat<<"\n\n";
        break;
      }
      case 4: {cout<<"\n\nThe 4th layer is under examination.\n\n";
        //--------------------------
        // D matrix for the 4th ply
        //--------------------------
        for (int i=0; i<nvars_ply_mat; i++) {
          Dmat_r(i) = Solve_FE2_Problem.Ply_4th_Dmat_r(i);
          for (int j=0; j<nvars_ply_mat; j++) {
            Dmat_rs(i,j) = Solve_FE2_Problem.Ply_4th_Dmat_rs(i,j);
          }
        }
        
        Dmat = Solve_FE2_Problem.Dmat_4th_Ply; //cout<<"Dmat: "<<Dmat<<"\n\n";
        break;
      }
    }
    
    //theta_angle = probdata.PlyAngle(probdata.ExaminedPly - 1)*(M_PI/180.0); cout<<"\n the angle "<<theta_angle<<endl;
    theta_angle = PlyAngle_new(probdata.ExaminedPly - 1)*(M_PI/180.0);
    
    // ==============================================
    //
    // Recover the stress at Gauss points for specific element(s)
    //
    // ==============================================
    
    FE2_PostProcStressForReliability_Zeroth Calc_Stress(m_field_Macro,"DISP_MACRO",Dmat);
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress); CHKERRQ(ierr);
    StressGP.clear(); REL_Stress_Transformation(theta_angle, Calc_Stress.StressGP, StressGP);
    //cout<<"Stress at GP in xyz: "<<Calc_Stress.StressGP<<endl;
    cout<<"Stress at GP in 123: "<<StressGP<<endl;

    // ==============================================
    //
    // Recover the first-order derivatives of stress
    //
    // ==============================================
    
    for (unsigned i=0; i<nvars_ply_mat; i++) {
      var_name.clear(); var_name = probdata.NameVars[i+1];
      
      ostringstream ss_field;
      ss_field.str(""); ss_field.clear();
      ss_field << "DISP_MACRO" << stochastic_fields_ply[i];
      cout<<"The first-order field is "<< ss_field.str().c_str() <<endl;
      
      FE2_PostProcStressForReliability_First Calc_Stress_r(m_field_Macro,"DISP_MACRO",ss_field.str().c_str(),Dmat,Dmat_r(i));
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_r);  CHKERRQ(ierr);
      StressGP_r(i).resize(3,3); StressGP_r(i).clear();
      //cout<<"The variable name is "<<probdata.NameVars[i]<<endl;
      
      if (var_name.compare(0,6,"theta1") == 0) {
        // w.r.t. ply angle for the 1st-layer
        if (probdata.ExaminedPly == 1) {
          REL_Stress_Transformation_Theta(theta_angle, Calc_Stress.StressGP, Calc_Stress_r.StressGP_r, StressGP_r(i));
        } else {
          REL_Stress_Transformation(theta_angle, Calc_Stress_r.StressGP_r, StressGP_r(i));
        }
      }
      else if (var_name.compare(0,6,"theta2") == 0) {
        // w.r.t. ply angle for the 2nd-layer
        if (probdata.ExaminedPly == 2) {
          REL_Stress_Transformation_Theta(theta_angle, Calc_Stress.StressGP, Calc_Stress_r.StressGP_r, StressGP_r(i));
        } else {
          REL_Stress_Transformation(theta_angle, Calc_Stress_r.StressGP_r,StressGP_r(i));
        }
      }
      else if (var_name.compare(0,6,"theta3") == 0) {
        // w.r.t. ply angle for the 3rd-layer
        if (probdata.ExaminedPly == 3) {
          REL_Stress_Transformation_Theta(theta_angle, Calc_Stress.StressGP, Calc_Stress_r.StressGP_r, StressGP_r(i));
        } else {
          REL_Stress_Transformation(theta_angle, Calc_Stress_r.StressGP_r,StressGP_r(i));
        }
      }
      else if (var_name.compare(0,6,"theta4") == 0) {
        // w.r.t. ply angle the 4th-layer
        if (probdata.ExaminedPly == 4) {
          REL_Stress_Transformation_Theta(theta_angle, Calc_Stress.StressGP, Calc_Stress_r.StressGP_r, StressGP_r(i));
        } else {
          REL_Stress_Transformation(theta_angle, Calc_Stress_r.StressGP_r, StressGP_r(i));
        }
      }
      else {
        REL_Stress_Transformation(theta_angle, Calc_Stress_r.StressGP_r, StressGP_r(i));
      }
      cout<<StressGP_r(i)<<endl;
      var_name.clear();
    }
    
    // Evaluate LSF and its gradient
    grad_g.resize(probdata.num_vars); grad_g.clear();
    //ierr = TheLSF.gfun(x,val_G,grad_g); CHKERRQ(ierr);
    
    switch (FailureCriterion) {
      case 43050: {
        //
        // Tsai-Wu failure criteria
        //
        NameOfFailureCriterion = "Tsai-Wu - 3D stress state";
        ierr = TheLSF.gfun_ply_Tsai_Wu_New(x, probdata.NameVars, probdata.MatStrength,
                                           StressGP, StressGP_r, StressGP_rs,
                                           val_G, grad_g, Hess_g,
                                           nvars_ply_mat, PSFE_order); CHKERRQ(ierr);
        break;
      }
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
    //if (istep>12) {
    //  conv_flag = 1;
    //}
  } while (conv_flag == 0);
  
  
  // Close beta value writting file
  if (beta_flag == 1) { BetaFile.close(); }
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 8: Calculate the estimate of reliability index               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  beta = inner_prod(alpha,u);
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 9: Adjust the estimate of reliability index by using         //
  //                the Second Order Reliability Method (SORM)                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  PSFE_order = probdata.AnalysisType;
  
  if (PSFE_order==2) {
    cout<<"\n\n";
    cout<<"/////////////////////////////////////////////////////////////////\n";
    cout<<"//                                                             //\n";
    cout<<"// The Second Order Reliability Method starts from here!       //\n";
    cout<<"//                                                             //\n";
    cout<<"/////////////////////////////////////////////////////////////////\n";
    
    // =======================================
    //
    // 9.1 Conducting analysis for RVE problem
    //
    // =======================================
    ierr = Solve_FE2_Problem.RVE_Dmat_Disp(m_field_RVE,
                                           stochastic_fields_rve,
                                           x,
                                           probdata.num_vars,
                                           probdata.NameVars,
                                           PSFE_order); CHKERRQ(ierr);
  
    // =======================================
    //
    // 9.2 Conducting analysis for Macro problem
    //
    // =======================================
    ierr = Solve_FE2_Problem.Macro_FE_REL_FSORM(m_field_Macro,
                                                nvars_rve_mat,
                                                nvars_ply_mat,
                                                name_rve_mat,
                                                name_ply_mat,
                                                stochastic_fields_ply,
                                                x,PlyAngle_new,
                                                NO_Layers,
                                                PSFE_order); CHKERRQ(ierr);
    
    // =======================================
    //
    // 9.3 Recover the 2nd-order derivatives of stress
    //
    // =======================================
    
    int sub_nvars = 0;
    int var_pos;
    for (int ivar = 0; ivar<nvars_ply_mat; ivar++) {
      // Get the first-order field
      ostringstream first_field_r;
      first_field_r.str(""); first_field_r.clear();
      first_field_r << "DISP_MACRO" << stochastic_fields_ply[ivar];
      
      for (int jvar = ivar; jvar<nvars_ply_mat; jvar++) {
        // Get the first-order derivative of field w.r.t. - j-th variable
        ostringstream first_field_s;
        first_field_s.str(""); first_field_s.clear();
        first_field_s << "DISP_MACRO" << stochastic_fields_ply[jvar];
        
        //
        if (jvar == ivar) {
          sub_nvars = sub_nvars + (jvar + 1);
          var_pos = sub_nvars;
        } else {
          var_pos = var_pos + jvar;
        }
        
        // Get the second-order field
        ostringstream second_field;
        second_field.str(""); second_field.clear();
        second_field << "DISP_MACRO" << stochastic_fields_ply[var_pos + nvars_ply_mat -1];
        cout<<"The second-order field is "<<second_field.str().c_str()<<endl;
        
        //
        FE2_PostProcStressForReliability_Second Calc_Stress_2nd(m_field_Macro,
                                                                "DISP_MACRO",
                                                                first_field_r.str().c_str(),
                                                                first_field_s.str().c_str(),
                                                                second_field.str().c_str(),
                                                                Dmat,
                                                                Dmat_r(ivar),Dmat_r(jvar),
                                                                Dmat_rs(ivar,jvar));
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_2nd);  CHKERRQ(ierr);
        //ublas::matrix<double> istress_rs; istress_rs.clear();
        REL_Stress_Transformation(theta_angle, Calc_Stress_2nd.StressGP_rs, StressGP_rs(ivar,jvar));
        cout<<"Stress_rs at GP: "<<StressGP_rs(ivar,jvar)<<endl;
      }
    }
    
    // =======================================
    //
    // 9.4 Calculate Hessian matrix of LSF in x-space
    //
    // =======================================
    
    ierr = TheLSF.gfun_ply_Tsai_Wu_New(x, probdata.NameVars, probdata.MatStrength,
                                       StressGP, StressGP_r, StressGP_rs,
                                       val_G, grad_g, Hess_g,
                                       nvars_ply_mat, PSFE_order); CHKERRQ(ierr);
    
    // =======================================
    //
    // 9.5 Construct orthogonal matrix
    //
    // =======================================
    SORM theSORM;
    
    cout<<"\nStart to calculate orthogonal matrix"<<endl;
    ublas::matrix<double> Qmatrix;
    theSORM.orthonormal_matrix(alpha,Qmatrix);
    cout<<"\nThe Q matrix: "<<Qmatrix<<endl;

    ublas::vector<double> dxdu;
    ublas::matrix<double> ddxddu;
    cout<<"\n The x: "<<x<<endl;
    cout<<"\n The u: "<<u<<endl;
    theSORM.d2x_dudu(x, u, probdata, dxdu, ddxddu);

    // =======================================
    //
    // 9.6 Transform Hessian matrix from x-space to u-space
    //
    // =======================================
    ublas::matrix<double> Hess_G;
    theSORM.Hessian_Matrix(dxdu, grad_g, Hess_g, Hess_G, ddxddu);

    ublas::matrix<double> A_Matrix;
    ublas::matrix<double> Temp_A_Matrix;
    Temp_A_Matrix = prod(Qmatrix, Hess_G);
    A_Matrix = prod(Temp_A_Matrix, trans(Qmatrix))/norm_2(grad_G);

    cout<<"\nHessian matrix: "<<Hess_G<<endl;
    cout<<"\nA matrix: "<<A_Matrix<<endl;

    ublas::matrix<double> New_A_Matrix;
    int Size_A; Size_A = A_Matrix.size1() - 1;
    New_A_Matrix.resize(Size_A,Size_A); New_A_Matrix.clear();
    for (int i = 0; i<Size_A; i++) {
      for (int j = 0; j<Size_A; j++) {
        New_A_Matrix(i,j) = A_Matrix(i,j);
      }
    }

    //LAPACK - eigenvalues and vectors. Applied twice for initial creates memory space
    ublas::matrix<double> eigen_vectors = New_A_Matrix;
    ublas::vector<double> kappa(Size_A);

    int lda, info, lwork = -1; lda = Size_A;
    double wkopt;
    info = lapack_dsyev('N','U',Size_A,&(eigen_vectors.data()[0]),lda,&(kappa.data()[0]),&wkopt,lwork);
    if(info != 0) SETERRQ1(PETSC_COMM_SELF,1,"is something wrong with lapack_dsyev info = %d",info);
    lwork = (int)wkopt;
    double work[lwork];
    info = lapack_dsyev('V','U',Size_A,&(eigen_vectors.data()[0]),lda,&(kappa.data()[0]),work,lwork);
    if(info != 0) SETERRQ1(PETSC_COMM_SELF,1,"is something wrong with lapack_dsyev info = %d",info);
    
    
    using boost::math::normal_distribution;
    normal_distribution<> snorm(0,1);
    double pf_Breitung = 0.0;
    for (int i = 0; i<Size_A; i++) {
      pf_Breitung = pf_Breitung + 1/sqrt(1 + kappa(i)*beta);
    }
    pf_Breitung = cdf(snorm,-beta)*pf_Breitung;
    double beta_Breitung = 0.0;
    //beta_Breitung = -quantile(snorm,pf_Breitung);
    
    cout<<"\nThe eigen values are: "<<kappa<<endl;
    cout<<"\nThe Breitung probability of failure is: "<<pf_Breitung<<endl;
    //cout<<"\nThe Breitung reliability index is: "<<beta_Breitung<<endl;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //                               FINISH                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  finish_time = clock();
  total_time  = (double)(finish_time - start_time)/CLOCKS_PER_SEC;
  
  cout<<"\n\n*************************************************\n*\n";
  
  
  if (istep == istep_max) {
    cout<<"*  The maximum number of iteration is reached.\n";
  } else {
    cout<<"*  The number of iterations is: "<<istep<<".\n";
  }
  
  
  cout<<"*  Optimal reliability index (beta) is: "<<beta<<endl;
  cout<<"*  The failure criterion used is: "<<NameOfFailureCriterion<<endl;
  cout<<"*  Elapsed time is "<<total_time<<" seconds.\n";
  cout<<"*  The program finishes !!! \n*\n";
  cout<<"*************************************************"<<endl;
  
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
  
}
