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
// To compute limit state function and its gradient w.r.t. the given trail
// checking point
//

void gfun(Stochastic_Model probdata,
          ublas::vector<double> x,
          double &val_lsf,
          ublas::vector<double> &grad_lsf,
          ublas::matrix<double> &Hess_lsf) {
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
  int nvars;
  nvars = probdata.num_vars;
  grad_lsf.resize(nvars); grad_lsf.clear();
  Hess_lsf.resize(nvars,nvars); Hess_lsf.clear();
  
  /*
  // ================================
  //
  // Example 1
  //
  // ================================
  // Evaluate the limit state function
  val_lsf = x(0)*x(1)-78.12*x(2);
  
  // Evaluate the partial derovative of the LSF w.r.t. basic variables
  // w.r.t. the 1st random variable
  grad_lsf(0) = x(1);
  // w.r.t. the 2nd random variable
  grad_lsf(1) = x(0);
  // w.r.t. the 3rd random variable
  grad_lsf(2) = -78.12;
  
  // Evaluate the 2nd-order partial derivatives of the LSF w.r.t. basic variables
  // w.r.t. the 1st random variable
  Hess_lsf(0,0) = 0;
  Hess_lsf(0,1) = 1;
  Hess_lsf(0,2) = 0;
  // w.r.t. the 2nd random variable
  Hess_lsf(1,0) = 1;
  Hess_lsf(1,1) = 0;
  Hess_lsf(1,2) = 0;
  // w.r.t. the 3rd random variable
  Hess_lsf(2,0) = 0;
  Hess_lsf(2,1) = 0;
  Hess_lsf(2,2) = 0;
  */
  
  /*
  val_lsf = x(0)*x(1)-600*x(2);
  
  grad_lsf(0) = x(1);
  grad_lsf(1) = x(0);
  grad_lsf(2) = -600;*/
  
  
  // ================================
  //
  // Example 3
  //
  // ================================
  // Evaluate the limit state function
//  val_lsf = pow(x(0),4) + 2*pow(x(1),4) - 20;
//  
//  // Evaluate the 1st-order partial derivatives of the LSF w.r.t. basic variables
//  // w.r.t. the 1st random variable
//  grad_lsf(0) = 4*pow(x(0),3);
//  // w.r.t. the 2nd random variable
//  grad_lsf(1) = 8*pow(x(1),3);
//  
//  // Evaluate the 2nd-order partial derivatives of the LSF w.r.t. basic variables
//  // w.r.t. the 1st random variable
//  Hess_lsf(0,0) = 12*pow(x(0),2);
//  // w.r.t. the 2nd random variable
//  Hess_lsf(1,1) = 24*pow(x(1),2);
  

  /**/
  // ================================
  //
  // Example 4
  //
  // ================================
  // Evaluate the limit state function
  val_lsf = x(2) - sqrt(300*pow(x(0),2) + 1.92*pow(x(1),2));
  
  // Evaluate the 1st-order partial derivatives of the LSF w.r.t. basic variables
  // w.r.t. the 1st random variable
  grad_lsf(0) = - pow(300*pow(x(0),2) + 1.92*pow(x(1),2),-0.5)*300*x(0);
  // w.r.t. the 2nd random variable
  grad_lsf(1) = - pow(300*pow(x(0),2) + 1.92*pow(x(1),2),-0.5)*1.92*x(1);
  // w.r.t. the 3rd random variable
  grad_lsf(2) = 1;
  
  // Evaluate the 2nd-order partial derivatives of the LSF w.r.t. basic variables
  // w.r.t. the 1st random variable
  double a;
  a = pow(300*pow(x(0),2) + 1.92*pow(x(1),2),-1.5)*576;
  Hess_lsf(0,0) = -pow(x(1),2);
  Hess_lsf(0,1) = x(0)*x(1);
  Hess_lsf(0,2) = 0;
  // w.r.t. the 2nd random variable
  Hess_lsf(1,0) = x(0)*x(1);
  Hess_lsf(1,1) = -pow(x(0),2);
  Hess_lsf(1,2) = 0;
  // w.r.t. the 3rd random variable
  Hess_lsf(2,0) = 0;
  Hess_lsf(2,1) = 0;
  Hess_lsf(2,2) = 0;
  
  Hess_lsf = a*Hess_lsf;
  /**/
  cout<<"\n The Hessian matrix: "<<Hess_lsf<<endl;
  
}


ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

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
  
  /*****************************************************************************
   *
   * Read parameters from line command
   *
   ****************************************************************************/
  
//  PetscBool flg = PETSC_TRUE;
//  char mesh_file_name[255];
//  
//  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_RVE",mesh_file_name,255,&flg); CHKERRQ(ierr);
//  if(flg != PETSC_TRUE) {
//    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_RVE (MESH FILE NEEDED)");
//  }
//  
//  PetscInt order_RVE;
//  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order_RVE",&order_RVE,&flg); CHKERRQ(ierr);
//  if(flg != PETSC_TRUE) {
//    order_RVE = 1;
//  }
  
 
  // ===========================================================================
  //
  //  C. RELIABILITY ANLYSIS
  //
  // ===========================================================================
  
  cout<<"\n\n";
  cout<<"///////////////////////////////////////////////////////////////////\n";
  cout<<"//                                                               //\n  ";
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
  ReliabOpt.step_code     = 0.05;       // 0: step size by Armijo rule, otherwise: given value is the step size
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
  cout<<"Search step: "<<step_code<<endl;
  
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
    x.resize(probdata.num_vars); x.clear();
    ierr = my_nataf_transformation.u_to_x(u,probdata.num_vars,probdata.marg,probdata.Lo,x,detj); CHKERRQ(ierr);
    
    cout<<"x is "<<x<<endl;
    cout<<"u is "<<u<<endl;
    
    // Jacobian
    dudx.resize(probdata.num_vars,probdata.num_vars); dudx.clear();
    ierr = my_nataf_transformation.Jacobian_u_x(x,u,probdata.num_vars,probdata.marg,probdata.Lo,probdata.inv_Lo,dudx); CHKERRQ(ierr);
    
    inv_dudx.resize(probdata.num_vars,probdata.num_vars); inv_dudx.clear();
    inv_dudx = gjinverse(dudx,singular);
    
    // Evaluate LSF and its gradient
    grad_g.resize(probdata.num_vars); grad_g.clear();
    gfun(probdata, x,val_G, grad_g, Hess_g);
    
//    //ierr = TheLSF.gfun(x,val_G,grad_g); CHKERRQ(ierr);
//    
//    switch (FailureCriterion) {
//      case 43050: {
//        //
//        // Tsai-Wu failure criteria
//        //
//        NameOfFailureCriterion = "Tsai-Wu - 3D stress state";
//        ierr = TheLSF.gfun_ply_Tsai_Wu_New(x, probdata.NameVars, probdata.MatStrength,
//                                           StressGP, StressGP_r, StressGP_rs,
//                                           val_G, grad_g, Hess_g,
//                                           nvars_ply_mat, PSFE_order); CHKERRQ(ierr);
//        break;
//      }
//    }
//    
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
    if (istep>8) {
      //conv_flag = 1;
    }
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
  
  double beta_Breitung;
  double pf_Breitung;
  
  if (PSFE_order==2) {
    cout<<"\n\n";
    cout<<"/////////////////////////////////////////////////////////////////\n";
    cout<<"//                                                             //\n";
    cout<<"// The Second Order Reliability Method starts from here!       //\n";
    cout<<"//                                                             //\n";
    cout<<"/////////////////////////////////////////////////////////////////\n";
            
    // =======================================
    //
    // 9.4 Calculate Hessian matrix of LSF in x-space
    //
    // =======================================
    
//    ierr = TheLSF.gfun_ply_Tsai_Wu_New(x, probdata.NameVars, probdata.MatStrength,
//                                       StressGP, StressGP_r, StressGP_rs,
//                                       val_G, grad_g, Hess_g,
//                                       nvars_ply_mat, PSFE_order); CHKERRQ(ierr);
    
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
    pf_Breitung = cdf(snorm,-beta);
    for (int i = 0; i<Size_A; i++) {
      pf_Breitung = pf_Breitung/sqrt(1 + kappa(i)*beta);
    }
    beta_Breitung = -quantile(snorm,pf_Breitung);
    
    cout<<"\nThe eigen values are: "<<kappa<<endl;
  }

  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //                               FINISH                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  finish_time = clock();
  total_time  = (double)(finish_time - start_time)/CLOCKS_PER_SEC;
  
  cout<<"\n\n************************************************************\n*\n";
  
  
  if (istep == istep_max) {
    cout<<"*  The maximum number of iteration is reached.\n";
  } else {
    cout<<"*  The number of iterations is: "<<istep<<".\n";
  }
  
  
    cout<<"*  Optimal reliability index (beta) is:    "<<beta<<endl;
  if (PSFE_order==2) {
    cout<<"*  The Breitung reliability index is:      "<<beta_Breitung<<endl;
    cout<<"*  The Breitung probability of failure is: "<<pf_Breitung<<endl;
  }
  cout<<"*  Elapsed time is "<<total_time<<" seconds.\n";
  cout<<"*  The program finishes !!! \n*\n";
  cout<<"************************************************************"<<endl;
  
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
  
}
