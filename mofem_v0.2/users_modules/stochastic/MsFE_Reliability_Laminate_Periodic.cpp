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

using namespace ObosleteUsersModules;

#include "ElasticFE_RVELagrange_Disp.hpp"
#include "ElasticFE_RVELagrange_Disp_Multi_Rhs.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Disp.hpp"

#include "ElasticFE_RVELagrange_Periodic.hpp"
#include "ElasticFE_RVELagrange_RigidBodyTranslation.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Periodic.hpp"

#include "RVEVolume.hpp"

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

#include <MCS_RNG.hpp>
#include <FE2_Macro_Solver_MCS.hpp>
//-----------------------
#include <cholesky.hpp>
#include <MatrixInverse.hpp> // download from http://proteowizard.sourceforge.net/dox/_matrix_inverse_8hpp.html

#include <ImportProbData.hpp>
#include <NatafTransformation.hpp>
#include <LimitStateFunction.hpp>

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

struct Stochastic_Model {
  int    num_vars;                       // Number of variables
  double dist_type;                      // Distribution type index
  double transf_type;                    // Type of joint distribution
  double R0_method;                      // Method for computation of the modified Nataf correlation matrix
  int flag_sens;                         // Flag for computation of sensitivities w.r.t. parameters
  int ExaminedPly;
  int AnalysisType;
  vector<string> NameVars;               // Name of random variables
  ublas::matrix<double> correlation;     // Correlation matrix
  ublas::matrix<double> marg;            // Marginal distribution for each random variable
  ublas::matrix<double> mod_correlation; // modified correlation matrix
  ublas::matrix<double> Lo;              // Chelosky decomposition
  ublas::matrix<double> inv_Lo;          // inverse of matrix Lo
  ublas::vector<double> MatStrength;     // Material strength
  ublas::vector<double> PlyAngle;        // Angle of orientation
};

struct Reliability_Results {
  double beta_FORM;                      // Reliability index from the FORM
  double beta_MCS;                       // Reliability index from the Crude MCS
  double beta_MCIS;                      // Reliability index from the Importance Sampling
  double prob_failure_FORM;              // Failure probability from the FORM
  double prob_failure_MCS;               // Failure probability from the Crudee MCS
  double prob_failure_MCIS;              // Failure probability from the Importance Sampling
  string NameOfFailureCriterion;         // Failure crition
  ublas::vector<double> DesignPoint;     // Design point
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


ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

const char* args[] = {
  "_r_Em", "_r_NUm",
  "_r_NUp", "_r_NUpz", "_r_Ep", "_r_Ez", "_r_Gzp",
  "_r_Ef", "_r_NUf",                                                // 1st order
  "_rs_EmEm", "_rs_NUmNUm",
  "_rs_NUpNUp", "_rs_NUpzNUpz", "_rs_EpEp", "_rs_EzEz", "_rs_GzpGzp",
  "_rs_EfEf", "_rs_NUfNUf"                                          // 2nd order
};

int nvars = 9;    // number of variables
int nders = 18;   // number of partial derivatives (firsr- and second- order)
vector<string> stochastic_fields(args, args + 18);


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

//==============================================================================
//
// Define class and multindex container to store data for triangles on the
// boundary of the RVE (it cannot be defined within main)
//
//==============================================================================

struct Face_CenPos_Handle {
  double xcoord, ycoord, zcoord;
  const EntityHandle  Tri_Hand;
  Face_CenPos_Handle(double _xcoord, double _ycoord,  double _zcoord,  const EntityHandle _Tri_Hand):xcoord(_xcoord),
  ycoord(_ycoord), zcoord(_zcoord), Tri_Hand(_Tri_Hand) {}
};

struct xcoord_tag {}; //tags to used in the multindex container
struct ycoord_tag {};
struct zcoord_tag {};
struct Tri_Hand_tag {};
struct Composite_xyzcoord {};

typedef multi_index_container<
Face_CenPos_Handle,
indexed_by<
ordered_non_unique<
tag<xcoord_tag>, member<Face_CenPos_Handle,double,&Face_CenPos_Handle::xcoord> >,

ordered_non_unique<
tag<ycoord_tag>, member<Face_CenPos_Handle,double,&Face_CenPos_Handle::ycoord> >,

ordered_non_unique<
tag<zcoord_tag>, member<Face_CenPos_Handle,double,&Face_CenPos_Handle::zcoord> >,

ordered_unique<
tag<Tri_Hand_tag>, member<Face_CenPos_Handle,const EntityHandle,&Face_CenPos_Handle::Tri_Hand> >,

ordered_unique<
tag<Composite_xyzcoord>,
composite_key<
Face_CenPos_Handle,
member<Face_CenPos_Handle,double,&Face_CenPos_Handle::xcoord>,
member<Face_CenPos_Handle,double,&Face_CenPos_Handle::ycoord>,
member<Face_CenPos_Handle,double,&Face_CenPos_Handle::zcoord> > >
> > Face_CenPos_Handle_multiIndex;

#include <Reliability_Methods_Periodic.hpp>

int main(int argc, char *argv[]) {
  
  clock_t start_time, finish_time;
  double total_time;
  start_time = clock();
  
  
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
  ierr = m_field_RVE.add_field("Lagrange_mul_disp",H1,field_rank); CHKERRQ(ierr);  //For lagrange multipliers to control the periodic motion
  ierr = m_field_RVE.add_field("Lagrange_mul_disp_rigid_trans",NOFIELD,3); CHKERRQ(ierr);  //To control the rigid body motion (3 Translations)
  
  // Stochastic fields for perturbation method
  for(int ii=0; ii < nders; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields[ii];
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
  ierr = m_field_RVE.add_finite_element("Lagrange_FE_rigid_trans"); CHKERRQ(ierr); //For rigid body control
  
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
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields[ii];
    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_RVE.modify_finite_element_add_field_data("ELASTIC_FE_RVE",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  //FE Transverse Isotropic
  ierr = m_field_RVE.modify_finite_element_add_field_row("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE","POTENTIAL_FIELD"); CHKERRQ(ierr);
  
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields[ii];
    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE",ss_field.str().c_str()); CHKERRQ(ierr);
  }

  //======================================================================================================
  // Define rows/cols and element data for C and CT (for lagrange multipliers)
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
  // Define rows/cols and element data for C1 and C1T (for lagrange multipliers to contol the rigid body motion)
  //======================================================================================================
  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = m_field_RVE.modify_finite_element_add_field_row("Lagrange_FE_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("Lagrange_FE_rigid_trans","DISP_RVE"); CHKERRQ(ierr);
  
  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = m_field_RVE.modify_finite_element_add_field_col("Lagrange_FE_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_row("Lagrange_FE_rigid_trans","DISP_RVE"); CHKERRQ(ierr);
  
  //data
  ierr = m_field_RVE.modify_finite_element_add_field_data("Lagrange_FE_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("Lagrange_FE_rigid_trans","DISP_RVE"); CHKERRQ(ierr);
  //======================================================================================================
  
  //define problems
  ierr = m_field_RVE.add_problem("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  
  //set finite elements for problem
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","Lagrange_FE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","Lagrange_FE_rigid_trans"); CHKERRQ(ierr);
  
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
  for(int ii=0; ii < nders; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields[ii];
    ierr = m_field_RVE.add_ents_to_field_by_TETs(0,ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  ierr = m_field_RVE.add_ents_to_finite_element_by_TETs(meshset_Matrix,"ELASTIC_FE_RVE",true); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_finite_element_by_TETs(meshset_Fibre,"TRAN_ISO_FE_RVE",true); CHKERRQ(ierr);
  
  // Add finite element to lagrange element for rigid body translation
  Range Tris_NewWholeMesh, Tri_OldNewSurf, SurfacesFaces;
  ierr = m_field_RVE.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBTRI,Tris_NewWholeMesh); CHKERRQ(ierr);
  ierr = m_field_RVE.get_cubit_msId_entities_by_dimension(103,SIDESET,2,Tri_OldNewSurf,true); CHKERRQ(ierr);
  SurfacesFaces = intersect(Tris_NewWholeMesh,Tri_OldNewSurf);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  
  // Create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab_RVE.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
  rval = moab_RVE.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = m_field_RVE.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);

  ierr = m_field_RVE.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_FE_rigid_trans"); CHKERRQ(ierr);

  //=======================================================================================================
  // Add Periodic Prisms Between Triangles on -ve and +ve faces to implement periodic boundary conditions
  //======================================================================================================= 
  //if (choise_value == HOMOBCPERIODIC) {
  //Populating the Multi-index container with -ve triangles
    Range Tri_OldNewSurfNeg, SurTrisNeg;
    ierr = m_field_RVE.get_cubit_msId_entities_by_dimension(101,SIDESET,2,Tri_OldNewSurfNeg,true); CHKERRQ(ierr);
    SurTrisNeg = intersect(Tris_NewWholeMesh,Tri_OldNewSurfNeg);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET 101 = %d\n",SurTrisNeg.size()); CHKERRQ(ierr);

    Face_CenPos_Handle_multiIndex Face_CenPos_Handle_varNeg, Face_CenPos_Handle_varPos;
    double TriCen[3], coords_Tri[9];

    double roundfact=1000.0;
    for(Range::iterator it = SurTrisNeg.begin(); it!=SurTrisNeg.end();  it++) {
        // cout<<"count1 ="<<count1<<endl;
        const EntityHandle* conn_face;  int num_nodes_Tri;

        //get nodes attached to the face
        rval = moab_RVE.get_connectivity(*it,conn_face,num_nodes_Tri,true); CHKERR_PETSC(rval);
        //get nodal coordinates
        rval = moab_RVE.get_coords(conn_face,num_nodes_Tri,coords_Tri); CHKERR_PETSC(rval);

       //Find triangle centriod
        TriCen[0]= (coords_Tri[0]+coords_Tri[3]+coords_Tri[6])/3.0;
        TriCen[1]= (coords_Tri[1]+coords_Tri[4]+coords_Tri[7])/3.0;
        TriCen[2]= (coords_Tri[2]+coords_Tri[5]+coords_Tri[8])/3.0;

        //round values to 3 disimal places
        if(TriCen[0]>=0) TriCen[0]=double(int(TriCen[0]*roundfact+0.5))/roundfact;  else TriCen[0]=double(int(TriCen[0]*roundfact-0.5))/roundfact; //-ve and +ve value
        if(TriCen[1]>=0) TriCen[1]=double(int(TriCen[1]*roundfact+0.5))/roundfact;  else TriCen[1]=double(int(TriCen[1]*roundfact-0.5))/roundfact;
        if(TriCen[2]>=0) TriCen[2]=double(int(TriCen[2]*roundfact+0.5))/roundfact;  else TriCen[2]=double(int(TriCen[2]*roundfact-0.5))/roundfact;
         // cout<<"   TriCen[0]= "<<TriCen[0] << "   TriCen[1]= "<< TriCen[1] << "   TriCen[2]= "<< TriCen[2] <<endl;
        //fill the multi-index container with centriod coordinates and triangle handles
        Face_CenPos_Handle_varNeg.insert(Face_CenPos_Handle(TriCen[0], TriCen[1], TriCen[2], *it));
    }

    //Populating the Multi-index container with +ve triangles
    Range Tri_OldNewSurfPos, SurTrisPos;
    ierr = m_field_RVE.get_cubit_msId_entities_by_dimension(102,SIDESET,2,Tri_OldNewSurfPos,true); CHKERRQ(ierr);
    SurTrisPos = intersect(Tris_NewWholeMesh,Tri_OldNewSurfPos);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET 102 = %d\n",SurTrisPos.size()); CHKERRQ(ierr);

    for(Range::iterator it = SurTrisPos.begin(); it!=SurTrisPos.end();  it++) {
        const EntityHandle* conn_face;  int num_nodes_Tri;

        //get nodes attached to the face
        rval = moab_RVE.get_connectivity(*it,conn_face,num_nodes_Tri,true); CHKERR_PETSC(rval);
        //get nodal coordinates
        rval = moab_RVE.get_coords(conn_face,num_nodes_Tri,coords_Tri); CHKERR_PETSC(rval);

        //Find triangle centriod
        TriCen[0]= (coords_Tri[0]+coords_Tri[3]+coords_Tri[6])/3.0;
        TriCen[1]= (coords_Tri[1]+coords_Tri[4]+coords_Tri[7])/3.0;
        TriCen[2]= (coords_Tri[2]+coords_Tri[5]+coords_Tri[8])/3.0;

        //round values to 3 disimal places
        if(TriCen[0]>=0) TriCen[0]=double(int(TriCen[0]*roundfact+0.5))/roundfact;  else TriCen[0]=double(int(TriCen[0]*roundfact-0.5))/roundfact;
        if(TriCen[1]>=0) TriCen[1]=double(int(TriCen[1]*roundfact+0.5))/roundfact;  else TriCen[1]=double(int(TriCen[1]*roundfact-0.5))/roundfact;
        if(TriCen[2]>=0) TriCen[2]=double(int(TriCen[2]*roundfact+0.5))/roundfact;  else TriCen[2]=double(int(TriCen[2]*roundfact-0.5))/roundfact;
       // cout<<"TriCen[0]= "<<TriCen[0] << "   TriCen[1]= "<< TriCen[1] << "   TriCen[2]= "<< TriCen[2] <<endl;

        //fill the multi-index container with centriod coordinates and triangle handles
        Face_CenPos_Handle_varPos.insert(Face_CenPos_Handle(TriCen[0], TriCen[1], TriCen[2], *it));
    }

    //Find minimum and maximum X, Y and Z coordinates of the RVE (using multi-index container)
    double XcoordMin, YcoordMin, ZcoordMin, XcoordMax, YcoordMax, ZcoordMax;
    typedef Face_CenPos_Handle_multiIndex::index<xcoord_tag>::type::iterator Tri_Xcoord_iterator;
    typedef Face_CenPos_Handle_multiIndex::index<ycoord_tag>::type::iterator Tri_Ycoord_iterator;
    typedef Face_CenPos_Handle_multiIndex::index<zcoord_tag>::type::iterator Tri_Zcoord_iterator;
    Tri_Xcoord_iterator XcoordMin_it, XcoordMax_it;
    Tri_Ycoord_iterator YcoordMin_it, YcoordMax_it;
    Tri_Zcoord_iterator ZcoordMin_it, ZcoordMax_it;

    //XcoordMax_it-- because .end() will point iterator after the data range but .begin() will point the iteratore to the first value of range
    XcoordMin_it=Face_CenPos_Handle_varNeg.get<xcoord_tag>().begin();                  XcoordMin=XcoordMin_it->xcoord;
    XcoordMax_it=Face_CenPos_Handle_varPos.get<xcoord_tag>().end();    XcoordMax_it--; XcoordMax=XcoordMax_it->xcoord;
    YcoordMin_it=Face_CenPos_Handle_varNeg.get<ycoord_tag>().begin();                  YcoordMin=YcoordMin_it->ycoord;
    YcoordMax_it=Face_CenPos_Handle_varPos.get<ycoord_tag>().end();    YcoordMax_it--; YcoordMax=YcoordMax_it->ycoord;
    ZcoordMin_it=Face_CenPos_Handle_varNeg.get<zcoord_tag>().begin();                  ZcoordMin=ZcoordMin_it->zcoord;
    ZcoordMax_it=Face_CenPos_Handle_varPos.get<zcoord_tag>().end();    ZcoordMax_it--; ZcoordMax=ZcoordMax_it->zcoord;

    cout<<"XcoordMin "<<XcoordMin << "      XcoordMax "<<XcoordMax <<endl;
    cout<<"YcoordMin "<<YcoordMin << "      YcoordMax "<<YcoordMax <<endl;
    cout<<"ZcoordMin "<<ZcoordMin << "      ZcoordMax "<<ZcoordMax <<endl;

    //Creating Prisims between triangles on -ve and +ve faces
    typedef Face_CenPos_Handle_multiIndex::index<Tri_Hand_tag>::type::iterator Tri_Hand_iterator;
    Tri_Hand_iterator Tri_Neg;
    typedef Face_CenPos_Handle_multiIndex::index<Composite_xyzcoord>::type::iterator xyzcoord_iterator;
    xyzcoord_iterator Tri_Pos;
    Range PrismRange;
    double XPos, YPos, ZPos;
    //int count=0;

    // loop over -ve triangles to create prisims elemenet between +ve and -ve triangles
    // count1=1;
  for(Range::iterator it = SurTrisNeg.begin(); it!=SurTrisNeg.end();  it++) {
    //        cout<<"count1 ="<<count1<<endl;  count1++;
    Tri_Neg=Face_CenPos_Handle_varNeg.get<Tri_Hand_tag>().find(*it);
    //corresponding +ve triangle
    if(Tri_Neg->xcoord==XcoordMin){XPos=XcoordMax;         YPos=Tri_Neg->ycoord;  ZPos=Tri_Neg->zcoord;};
    if(Tri_Neg->ycoord==YcoordMin){XPos=Tri_Neg->xcoord;   YPos=YcoordMax;        ZPos=Tri_Neg->zcoord;};
    if(Tri_Neg->zcoord==ZcoordMin){XPos=Tri_Neg->xcoord;   YPos=Tri_Neg->ycoord;  ZPos=ZcoordMax;      };
    Tri_Pos=Face_CenPos_Handle_varPos.get<Composite_xyzcoord>().find(boost::make_tuple(XPos, YPos, ZPos));
    EntityHandle PrismNodes[6];
    vector<EntityHandle> TriNodesNeg, TriNodesPos;
    double CoordNodeNeg[9], CoordNodePos[9];
    rval = moab_RVE.get_connectivity(&(Tri_Neg->Tri_Hand),1,TriNodesNeg,true); CHKERR_PETSC(rval);
    rval = moab_RVE.get_connectivity(&(Tri_Pos->Tri_Hand),1,TriNodesPos,true); CHKERR_PETSC(rval);
    rval = moab_RVE.get_coords(&TriNodesNeg[0],3,CoordNodeNeg);  CHKERR_THROW(rval);
    rval = moab_RVE.get_coords(&TriNodesPos[0],3,CoordNodePos);  CHKERR_THROW(rval);
    for(int ii=0; ii<3; ii++){
      PrismNodes[ii]=TriNodesNeg[ii];
    }
    
    //Match exact nodes to each other to avoide the problem of twisted prisms
    double XNodeNeg, YNodeNeg, ZNodeNeg, XNodePos, YNodePos, ZNodePos;
    for(int ii=0; ii<3; ii++){
      if(Tri_Neg->xcoord==XcoordMin){XNodeNeg=XcoordMax;          YNodeNeg=CoordNodeNeg[3*ii+1];   ZNodeNeg=CoordNodeNeg[3*ii+2];};
      if(Tri_Neg->ycoord==YcoordMin){XNodeNeg=CoordNodeNeg[3*ii]; YNodeNeg=YcoordMax;              ZNodeNeg=CoordNodeNeg[3*ii+2];};
      if(Tri_Neg->zcoord==ZcoordMin){XNodeNeg=CoordNodeNeg[3*ii]; YNodeNeg=CoordNodeNeg[3*ii+1];   ZNodeNeg=ZcoordMax;};
      for(int jj=0; jj<3; jj++){
        XNodePos=CoordNodePos[3*jj]; YNodePos=CoordNodePos[3*jj+1]; ZNodePos=CoordNodePos[3*jj+2];
        
        if(XNodeNeg==XNodePos  &&  YNodeNeg==YNodePos  &&  ZNodeNeg==ZNodePos){
          PrismNodes[3+ii]=TriNodesPos[jj];
          break;
        }
      }
    }
    //prism nodes and their coordinates
    double CoordNodesPrisms[18];
    rval = moab_RVE.get_coords(PrismNodes,6,CoordNodesPrisms);  CHKERR_THROW(rval);
    EntityHandle PeriodicPrism;
    rval = moab_RVE.create_element(MBPRISM,PrismNodes,6,PeriodicPrism); CHKERR_PETSC(rval);
    PrismRange.insert(PeriodicPrism);
  }
//  }
  //cout<<"\n Constructing Prisms"<<endl;
  //cout<<"PrismRange "<<PrismRange<<endl;
  //Saving prisms in interface.vtk
  EntityHandle out_meshset1;
  rval = moab_RVE.create_meshset(MESHSET_SET,out_meshset1); CHKERR_PETSC(rval);
  rval = moab_RVE.add_entities(out_meshset1,PrismRange); CHKERR_PETSC(rval);
  rval = moab_RVE.write_file("Prisms.vtk","VTK","",&out_meshset1,1); CHKERR_PETSC(rval);
  cout << "Prisms.vtk output" <<endl;
  
  //Adding Prisims to Element Lagrange_elem (to loop over these prisims)
  EntityHandle PrismRangeMeshset;
  rval = moab_RVE.create_meshset(MESHSET_SET,PrismRangeMeshset); CHKERR_PETSC(rval);
  rval = moab_RVE.add_entities(PrismRangeMeshset,PrismRange); CHKERR_PETSC(rval);
  //    cout << PrismRange <<endl;
  ierr = m_field_RVE.seed_ref_level_3D(PrismRangeMeshset,problem_bit_level_RVE); CHKERRQ(ierr);
  //    mField.seed_finite_elements(PrismRange);
  ierr = m_field_RVE.get_entities_by_ref_level(bit_levels.back(),BitRefLevel().set(),out_meshset); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_finite_element_by_PRISMs(PrismRange, "Lagrange_FE"); CHKERRQ(ierr);
  
  //Adding only -ve surfaces to the field Lagrange_mul_disp (in periodic boundary conditions size of C (3M/2 x 3N))
  //to create meshset from range  SurTrisNeg
  EntityHandle SurTrisNegMeshset;
  rval = moab_RVE.create_meshset(MESHSET_SET,SurTrisNegMeshset); CHKERR_PETSC(rval);
  rval = moab_RVE.add_entities(SurTrisNegMeshset,SurTrisNeg); CHKERR_PETSC(rval);
  ierr = m_field_RVE.add_ents_to_field_by_TRIs(SurTrisNegMeshset,"Lagrange_mul_disp",2); CHKERRQ(ierr);
  
  /*****************************************************************************
   *
   * Set applied order
   * See reference for detals:
   *   Ainsworth M. and Coyle J. (2003) Hierarchic finite element bases on
   unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/
  ierr = m_field_RVE.set_field_order(0,MBTET,"DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBTRI,"DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,"DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"DISP_RVE",1); CHKERRQ(ierr);
  
  ierr = m_field_RVE.set_field_order(0,MBTRI,"Lagrange_mul_disp",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);
  
  int order_st = order_RVE;
  for(int ii=0; ii < nders; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields[ii];
    //    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_RVE.set_field_order(0,MBTET,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_RVE.set_field_order(0,MBTRI,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_RVE.set_field_order(0,MBEDGE,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_RVE.set_field_order(0,MBVERTEX,ss_field.str().c_str(),1); CHKERRQ(ierr);
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
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.add_field(ss_field.str().c_str(),H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  }
  ierr = m_field_Macro.add_field("DISP_MACRO_r_F",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field_Macro.add_field("DISP_MACRO_r_Theta",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field_Macro.add_field("DISP_MACRO_r_Theta_1st_Ply",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field_Macro.add_field("DISP_MACRO_r_Theta_2nd_Ply",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field_Macro.add_field("DISP_MACRO_r_Theta_3rd_Ply",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field_Macro.add_field("DISP_MACRO_r_Theta_4th_Ply",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  
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
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO_r_F"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO_r_Theta"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO_r_Theta_1st_Ply"); CHKERRQ(ierr);
  
  
  // Second layer
  if (NO_Layers > 1) {
    //Define rows/cols and element data
    ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_2nd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_2nd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO_r_Theta_2nd_Ply"); CHKERRQ(ierr);
    
    for(int ii=0; ii < nvars; ii++ ) {
      ostringstream ss_field;
      ss_field << "DISP_MACRO" << stochastic_fields[ii];
      ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply",ss_field.str().c_str()); CHKERRQ(ierr);
    }
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO_r_F"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO_r_Theta"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO_r_Theta_1st_Ply"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO_r_Theta_2nd_Ply"); CHKERRQ(ierr);
  }
  // Third layer
  if (NO_Layers > 2) {
    //Define rows/cols and element data
    ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_3rd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_3rd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO_r_Theta_3rd_Ply"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO_r_Theta_3rd_Ply"); CHKERRQ(ierr);
    
    for(int ii=0; ii < nvars; ii++ ) {
      ostringstream ss_field;
      ss_field << "DISP_MACRO" << stochastic_fields[ii];
      ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply",ss_field.str().c_str()); CHKERRQ(ierr);
    }
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO_r_F"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO_r_Theta"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO_r_Theta_1st_Ply"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO_r_Theta_2nd_Ply"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO_r_Theta_3rd_Ply"); CHKERRQ(ierr);
  }
  // Fourth layer
  if (NO_Layers > 3) {
    //Define rows/cols and element data
    ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_4th_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_4th_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO_r_Theta_4th_Ply"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO_r_Theta_4th_Ply"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO_r_Theta_4th_Ply"); CHKERRQ(ierr);
    
    for(int ii=0; ii < nvars; ii++ ) {
      ostringstream ss_field;
      ss_field << "DISP_MACRO" << stochastic_fields[ii];
      ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply",ss_field.str().c_str()); CHKERRQ(ierr);
    }
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO_r_F"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO_r_Theta"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO_r_Theta_1st_Ply"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO_r_Theta_2nd_Ply"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO_r_Theta_3rd_Ply"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO_r_Theta_4th_Ply"); CHKERRQ(ierr);
  }

  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_FE_MACRO_REL","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_FE_MACRO_REL","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO_r_F"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO_r_Theta"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO_r_Theta_1st_Ply"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO_r_Theta_2nd_Ply"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO_r_Theta_3rd_Ply"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO_r_Theta_4th_Ply"); CHKERRQ(ierr);
  
  
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
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.add_ents_to_field_by_TETs(0,ss_field.str().c_str()); CHKERRQ(ierr);
  }
  ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"DISP_MACRO_r_F"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"DISP_MACRO_r_Theta"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"DISP_MACRO_r_Theta_1st_Ply"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"DISP_MACRO_r_Theta_2nd_Ply"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"DISP_MACRO_r_Theta_3rd_Ply"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"DISP_MACRO_r_Theta_4th_Ply"); CHKERRQ(ierr);
  

  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_1st_Ply, "ELASTIC_1st_Ply"); CHKERRQ(ierr);
  if (NO_Layers > 1) {
    ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_2nd_Ply,"ELASTIC_2nd_Ply"); CHKERRQ(ierr);
  }
  if (NO_Layers > 2) {
    ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_3rd_Ply,"ELASTIC_3rd_Ply"); CHKERRQ(ierr);
  }
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
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.set_field_order(0,MBTET,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBTRI,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBEDGE,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBVERTEX,ss_field.str().c_str(),1); CHKERRQ(ierr);
  }
  // Applied force
  ierr = m_field_Macro.set_field_order(0,MBTET,   "DISP_MACRO_r_F",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,   "DISP_MACRO_r_F",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,  "DISP_MACRO_r_F",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"DISP_MACRO_r_F",1); CHKERRQ(ierr);
  // Ply orietation
  ierr = m_field_Macro.set_field_order(0,MBTET,   "DISP_MACRO_r_Theta",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,   "DISP_MACRO_r_Theta",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,  "DISP_MACRO_r_Theta",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"DISP_MACRO_r_Theta",1); CHKERRQ(ierr);
  
  ierr = m_field_Macro.set_field_order(0,MBTET,   "DISP_MACRO_r_Theta_1st_Ply",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,   "DISP_MACRO_r_Theta_1st_Ply",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,  "DISP_MACRO_r_Theta_1st_Ply",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"DISP_MACRO_r_Theta_1st_Ply",1); CHKERRQ(ierr);
  
  ierr = m_field_Macro.set_field_order(0,MBTET,   "DISP_MACRO_r_Theta_2nd_Ply",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,   "DISP_MACRO_r_Theta_2nd_Ply",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,  "DISP_MACRO_r_Theta_2nd_Ply",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"DISP_MACRO_r_Theta_2nd_Ply",1); CHKERRQ(ierr);
  
  ierr = m_field_Macro.set_field_order(0,MBTET,   "DISP_MACRO_r_Theta_3rd_Ply",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,   "DISP_MACRO_r_Theta_3rd_Ply",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,  "DISP_MACRO_r_Theta_3rd_Ply",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"DISP_MACRO_r_Theta_3rd_Ply",1); CHKERRQ(ierr);
  
  ierr = m_field_Macro.set_field_order(0,MBTET,   "DISP_MACRO_r_Theta_4th_Ply",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,   "DISP_MACRO_r_Theta_4th_Ply",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,  "DISP_MACRO_r_Theta_4th_Ply",order_st); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"DISP_MACRO_r_Theta_4th_Ply",1); CHKERRQ(ierr);
  
  //
  ierr = m_field_Macro.set_field_order(0,MBTET,   "MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,   "MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,  "MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
 
  // load
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
  ReliabOpt.step_code     = 1;       // 0: step size by Armijo rule, otherwise: given value is the step size
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
  Stochastic_Model probdata;
  
  // Import data from textile file
  ImportProbData readprobdata;
  
  ierr = readprobdata.ProbdataFileIn(); CHKERRQ(ierr);
  probdata.marg        = readprobdata.MargProb;
  probdata.correlation = readprobdata.CorrMat;
  probdata.num_vars    = readprobdata.NumVars;
  probdata.MatStrength = readprobdata.MatStrength;
  probdata.NameVars    = readprobdata.NameVars;
  probdata.PlyAngle    = readprobdata.PlyAngle;
  probdata.ExaminedPly = readprobdata.ExaminedLayer;
  probdata.AnalysisType = readprobdata.AnalysisType;
  
  int FailureCriterion; // 1: Tsai-Wu, 2: Tsai-Hill
  string NameOfFailureCriterion;
  FailureCriterion = readprobdata.FailureCriterion;
  Reliability_Results MsFE_Reliability_Results;
  step_code = readprobdata.SearchStep;
  
  switch (probdata.AnalysisType) {
    case 1: {
      cout<<"\n\n";
      cout<<"====================================="<<endl;
      cout<<"   First-order reliability method"<<endl;
      cout<<"====================================="<<endl;
      cout<<"\n\n";
      
      Reliability_Methods_Periodic theREL;
      theREL.FORM(m_field_RVE,m_field_Macro,nvars,nders,stochastic_fields,
                  NO_Layers,probdata,ReliabOpt,FailureCriterion,
                  MsFE_Reliability_Results);
      break;
    }
//    case 2: {
//      cout<<"\n\n";
//      cout<<"====================================="<<endl;
//      cout<<"   Crude Monte Carlo simulation"<<endl;
//      cout<<"====================================="<<endl;
//      cout<<"\n\n";
//      
//      Reliability_Methods_Periodic theREL;
//      theREL.Crude_MCS(m_field_RVE,m_field_Macro,NO_Layers,probdata,ReliabOpt,FailureCriterion,
//                       MsFE_Reliability_Results);
//      break;
//    }
    case 3: {
      cout<<"\n\n";
      cout<<"====================================="<<endl;
      cout<<"   Importance sampling method"<<endl;
      cout<<"====================================="<<endl;
      cout<<"\n\n";
      
      Reliability_Methods_Periodic theREL;
      //
      // Run FORM to get initial design point
      //
      theREL.FORM(m_field_RVE,m_field_Macro,nvars,nders,stochastic_fields,
                  NO_Layers,probdata,ReliabOpt,FailureCriterion,
                  MsFE_Reliability_Results);
      // MCIS
      theREL.MCIS(m_field_RVE,m_field_Macro,NO_Layers,probdata,ReliabOpt,FailureCriterion,
                  MsFE_Reliability_Results);
      
      break;
    }
      
  }
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //                               FINISH                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  finish_time = clock();
  total_time  = (double)(finish_time - start_time)/CLOCKS_PER_SEC;
  
  cout<<"*  Elapsed time is "<<total_time<<" seconds.\n";
  cout<<"*  The program finishes !!! \n*\n";
  cout<<"*************************************************"<<endl;
  
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
  
}
