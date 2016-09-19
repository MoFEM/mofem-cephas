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

extern "C" {
#include <gm_rule.h>
#include <ltqnorm.h>
}


#include <MoFEM.hpp>
using namespace MoFEM;

#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
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

using namespace boost::numeric;
//using namespace MoFEM;
using namespace std;

static char help[] = "...\n\n";


//------------------------------------------------------------------------------
// Construct data structure <Stochastic_Model> for collecting data representing
//   statistical information of inputs including probability distribution,
//   correlation matrix and etc.

struct Stochastic_Model {
  int nvars;                             // Number of variables
  double dist_type;                      // Distribution type index
  double transf_type;                    // Type of joint distribution
  double R0_method;                      // Method for computation of the modified Nataf correlation matrix
  int flag_sens;                         // Flag for computation of sensitivities w.r.t. parameters
  ublas::matrix<double> correlation;     // Correlation matrix
  ublas::matrix<double> marg;            // Marginal distribution for each random variable
  ublas::matrix<double> mod_correlation; // modified correlation matrix
  ublas::matrix<double> Lo;              // Chelosky decomposition
  ublas::matrix<double> inv_Lo;          // inverse of matrix Lo
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
// To count number of a specific character or a substring in a string;
//

void str_cnt(char *mstr,char substr,int &cnt) {
  int j=0;
  for (int i=0;i<strlen(mstr);i++)
  {
    if (*(mstr+i)==substr)
    {
      j++;
    }
  }
  cnt=j;
}

//------------------------------------------------------------------------------
// To determine the position of a specific charater in a string
//

void str_pos(char *mstr,char substr,int *pos) {
  int j=1;
  for (int i=0;i<strlen(mstr);i++) {
    if (*(mstr+i)==substr) {
      //cout<<*(mstr+i)<<'\t'<<i<<'\n';
      pos[j]=i;
      j++;
    }
  }
}


//------------------------------------------------------------------------------
// To transform x from original space to standard normal distribution probability
// space using Nataf transformation
//

double x_to_u(ublas::vector<double> x,
              Stochastic_Model probdata,
              ublas::vector<double> &u) {
  int nvars;
  nvars = probdata.nvars;
  ublas::vector<double> z(nvars);
  u.resize(nvars);
  
  using boost::math::normal_distribution;
  normal_distribution<> snorm(0,1);
  double imean, istd, iloc, ishape, iscale, ilambda, ilower, iupper;
  int dist_type;
  
  for (int i=0;i<nvars;i++) {
    dist_type = probdata.marg(i,0);
    switch (dist_type) { // distribution type
      case 1: {
        /*
         * distribution type 1: Normal
         *           parameter: location - mean
         *                      scale    - standard deviation
         */
        imean = probdata.marg(i,4);
        istd  = probdata.marg(i,5);
        normal_distribution<> mynorm(imean,istd);
        z(i) = quantile(snorm,cdf(mynorm,x(i)));
        break;
      }
      case 2: {
        /*
         * distribution type 2: lognormal
         *           parameter: location - mean of logrithm rv
         *                      scale    - std of logrithm rv
         */
        using boost::math::lognormal_distribution;
        iloc   = probdata.marg(i,4);
        iscale = probdata.marg(i,5);
        lognormal_distribution<> my_logn(iloc,iscale);
        z(i) = quantile(snorm,cdf(my_logn,x(i)));
        break;
      }
      case 3: {
        /*
         * distribution type 3: exponential
         *           parameter: lambda = 1/mu
         */
        using boost::math::exponential_distribution;
        ilambda = 1/probdata.marg(i,1);
        exponential_distribution<> my_exp(ilambda);
        z(i) = quantile(snorm,cdf(my_exp,x(i)));
        break;
      }
      case 4: {
        /*
         * distribution type 4: Extreme value distribution - I (Gumbel)
         *           parameter: shape
         *                      location
         *                      scale
         */
        using boost::math::extreme_value_distribution;
        iloc   = probdata.marg(i,4);
        iscale = probdata.marg(i,5);
        extreme_value_distribution<> my_EVI(iloc,iscale);
        z(i) = quantile(snorm,cdf(my_EVI,x(i)));
        break;
      }
      case 5: {
        /*
         * distribution type 5: Extreme value distribution - III (Weibull)
         *           parameter: shape
         *                      scale
         */
        using boost::math::weibull_distribution;
        iloc   = probdata.marg(i,4);
        iscale = probdata.marg(i,5);
        weibull_distribution<> my_wbl(iloc,iscale);
        z(i) = quantile(snorm,cdf(my_wbl,x(i)));
        break;
      }
      case 6: {
        /*
         * distribution type 6: Uniform distribution
         *           parameter: shape
         *                      scale
         */
        using boost::math::uniform_distribution;
        ilower = probdata.marg(i,4);
        iupper = probdata.marg(i,5);
        uniform_distribution<> my_unif(ilower,iupper);
        z(i) = quantile(snorm,cdf(my_unif,x(i)));
        break;
      }
      case 7: {
        /*
         * distribution type 7: Gamma distribution
         *           parameter: shape
         *                      scale
         */
        using boost::math::gamma_distribution;
        ishape = probdata.marg(i,4);
        iscale = probdata.marg(i,5);
        gamma_distribution<> my_gamma(ishape,iscale);
        z(i) = quantile(snorm,cdf(my_gamma,x(i)));
        break;
      }
      default:
        cout<<"The distribution type index should be value between 0 and 7!"<<endl;
        break;
    }

    cout<<i<<"th variable \t"<<z(i)<<endl;
  }

  u = prod(probdata.inv_Lo,z);
  
}


//------------------------------------------------------------------------------
// To transform u from standard normal distribution probability
// space to original space to obtain x using Nataf transformation
//

double u_to_x(ublas::vector<double> u,
              Stochastic_Model probdata,
              ublas::vector<double> &x,
              double &detj) {
  int nvars;
  nvars = probdata.nvars;
  x.resize(nvars);
  ublas::vector<double> z(nvars);
  z = prod(probdata.Lo,u);
  
  using boost::math::normal_distribution;
  normal_distribution<> snorm(0,1);
  double imean, istd, iloc, ishape, iscale, ilambda, ilower, iupper;
  int dist_type;
  detj = 1.0;
  
  for (int i=0;i<nvars;i++) {
    dist_type = probdata.marg(i,0);
    switch (dist_type) { // distribution type
      case 1: {
        /*
         * distribution type 1: Normal
         *           parameter: location - mean
         *                      scale    - standard deviation
         */
        imean = probdata.marg(i,4);
        istd  = probdata.marg(i,5);
        normal_distribution<> mynorm(imean,istd);
        x(i) = quantile(mynorm,cdf(snorm,z(i)));
        detj = detj*pdf(mynorm,x(i))/pdf(snorm,u(i));
        break;
      }
      case 2: {
        /*
         * distribution type 2: lognormal
         *           parameter: location - mean of logrithm rv
         *                      scale    - std of logrithm rv
         */
        using boost::math::lognormal_distribution;
        iloc   = probdata.marg(i,4);
        iscale = probdata.marg(i,5);
        lognormal_distribution<> my_logn(iloc,iscale);
        x(i) = quantile(my_logn,cdf(snorm,z(i)));
        detj = detj*pdf(my_logn,x(i))/pdf(snorm,u(i));
        break;
      }
      case 3: {
        /*
         * distribution type 3: exponential
         *           parameter: lambda = 1/mu
         */
        using boost::math::exponential_distribution;
        ilambda = 1/probdata.marg(i,1);
        exponential_distribution<> my_exp(ilambda);
        x(i) = quantile(my_exp,cdf(snorm,z(i)));
        detj = detj*pdf(my_exp,x(i))/pdf(snorm,u(i));
        break;
      }
      case 4: {
        /*
         * distribution type 4: Extreme value distribution - I (Gumbel)
         *           parameter: shape
         *                      location
         *                      scale
         */
        using boost::math::extreme_value_distribution;
        iloc   = probdata.marg(i,4);
        iscale = probdata.marg(i,5);
        extreme_value_distribution<> my_EVI(iloc,iscale);
        x(i) = quantile(my_EVI,cdf(snorm,z(i)));
        detj = detj*pdf(my_EVI,x(i))/pdf(snorm,u(i));
        break;
      }
      case 5: {
        /*
         * distribution type 5: Extreme value distribution - III (Weibull)
         *           parameter: shape
         *                      scale
         */
        using boost::math::weibull_distribution;
        iloc   = probdata.marg(i,4);
        iscale = probdata.marg(i,5);
        weibull_distribution<> my_wbl(iloc,iscale);
        x(i) = quantile(my_wbl,cdf(snorm,z(i)));
        detj = detj*pdf(my_wbl,x(i))/pdf(snorm,u(i));
        break;
      }
      case 6: {
        /*
         * distribution type 6: Uniform distribution
         *           parameter: shape
         *                      scale
         */
        using boost::math::uniform_distribution;
        ilower = probdata.marg(i,4);
        iupper = probdata.marg(i,5);
        uniform_distribution<> my_unif(ilower,iupper);
        x(i) = quantile(my_unif,cdf(snorm,z(i)));
        detj = detj*pdf(my_unif,x(i))/pdf(snorm,u(i));
        break;
      }
      case 7: {
        /*
         * distribution type 7: Gamma distribution
         *           parameter: shape
         *                      scale
         */
        using boost::math::gamma_distribution;
        ishape = probdata.marg(i,4);
        iscale = probdata.marg(i,5);
        gamma_distribution<> my_gamma(ishape,iscale);
        x(i) = quantile(my_gamma,cdf(snorm,z(i)));
        detj = detj*pdf(my_gamma,x(i))/pdf(snorm,u(i));
        break;
      }
      default:
        cout<<"The distribution type index should be value between 0 and 7!"<<endl;
        break;
    }
    
    cout<<i<<"th variable \t"<<x(i)<<"\t"<<u(i)<<endl;
  }
  
}

//------------------------------------------------------------------------------
// To compute modified correlation matrix by using Nataf transformation
//

void modify_correlation_mat(Stochastic_Model &probdata) {
  /*
   * [Liu and Der Kiureghian, 1986] provided an empirical formula to calculate
   *   transformed correlation coefficient.
   *       R = rho_mod/rho_init
   *       R = a + b*Vi + c*Vi^2 + d*rho_init + e*rho_init^2 + f*rho_init*Vi
   *             + g*Vj + k*rho_init*Vj + l*Vi*Vj
   *       where a, b, c, d, e, f, g, k and l are coefficients, which are given
   *       below accoording to the reference.
   *
   * [Melchers, 1999] Structural reliability analysis and prediction
   * [Ditlevsen, 2005] Structural reliability methods
   *
   */
  vector<vector<vector<double> > > coef;
  
  // set up size for 3D array 'coef'
  coef.resize(10);
  for (int i = 0; i < 10; i++) {
    coef[i].resize(10);
    for (int j = 0; j < 10; j++) {
      coef[i][j].resize(10);
    }
  }
  // Initialization the 3D arrayv 'coef' with zeros
  for (int i = 0; i<10; i++) {
    for (int j = 0; j<10; j++) {
      for (int k = 0; k<10; k++) {
        coef[i][j][k] = 0.0;
      }
    }
  }
  
  // Assign values to coef Xj Xi
  // .1 Case: Xj of Normal and Xi of ?
  //   .1 Normal to Normal
  coef[0][0][0] = 1.0;
  //   .2 Normal to Lognormal
  coef[0][1][0] = 0.0; // Value will be given by a specific formula
  //   .3 Normal to Exponential
  coef[0][2][0] = 1.107;
  //   .4 Normal to Gumbel
  coef[0][3][0] = 1.031;
  //   .5 Normal to Weibull
  coef[0][4][0] =  1.031;
  coef[0][4][1] = -0.195;
  coef[0][4][2] =  0.0328;
  //   .6 Normal to Uniform
  coef[0][5][0] = 1.023;
  
  // .2 Case: Xj of Lognormal to Xi of ?
  //   .1 Lognormal to Normal
  coef[1][0][0] = 0.0; // Value will be given by a specific formula
  //   .2 Lognormal to Lognormal
  coef[1][1][0] = 0.0; // Value will be given by a specific formula
  //   .3 Lognormal to Exponential
  coef[1][2][0] =  1.098;
  coef[1][2][1] =  0.019;
  coef[1][2][2] =  0.0303;
  coef[1][2][3] =  0.003;
  coef[1][2][4] =  0.025;
  coef[1][2][5] = -0.437;
  //   .4 Lognormal to Gumbel
  coef[1][3][0] =  1.029;
  coef[1][3][1] =  0.014;
  coef[1][3][2] =  0.233;
  coef[1][3][3] =  0.001;
  coef[1][3][4] =  0.004;
  coef[1][3][5] = -0.197;
  //   .5 Lognormal to Weibull
  coef[1][4][0] =  1.031;
  coef[1][4][1] =  0.011;
  coef[1][4][2] =  0.220;
  coef[1][4][3] =  0.052;
  coef[1][4][4] =  0.002;
  coef[1][4][5] =  0.005;
  coef[1][4][6] = -0.210;
  coef[1][4][7] =  0.350;
  coef[1][4][8] = -0.174;
  coef[1][4][9] =  0.009;
  //   .6 Lognormal to Uniform
  coef[1][5][0] =  1.019;
  coef[1][5][1] =  0.014;
  coef[1][5][2] =  0.249;
  coef[1][5][3] =  0.0;
  coef[1][5][4] =  0.010;
  
  // .3 Case: Xj of Exponential to Xi of ?
  //   .1 Exponential to Normal
  coef[2][0][0] = 1.107;
  //   .1 Exponential to Lognormal
  coef[2][1][0] =  1.098;
  coef[2][1][1] =  0.019;
  coef[2][1][2] =  0.303;
  coef[2][1][3] =  0.003;
  coef[2][1][4] =  0.025;
  coef[2][1][5] = -0.437;
  //   .2 Exponential to Exponetial
  coef[2][2][0] =  1.229;
  coef[2][2][1] =  0.0;
  coef[2][2][2] =  0.0;
  coef[2][2][3] = -0.367;
  coef[2][2][4] =  0.153;
  //   .3 Exponential to Gumbel
  coef[2][3][0] =  1.142;
  coef[2][3][1] =  0.0;
  coef[2][3][2] =  0.0;
  coef[2][3][3] = -0.154;
  coef[2][3][4] =  0.031;
  //   .4 Exponential to Weibull
  //      !!! difference exists between Ditlevsen & Melchers
  //          the values used here follow those given in Ditlevsen's book.
  coef[2][4][0] =  1.147;
  coef[2][4][1] = -0.271;
  coef[2][4][2] =  0.459;
  coef[2][4][3] =  0.145;
  coef[2][4][4] =  0.010;
  coef[2][4][5] = -0.467;
  //   .5 Exponential to Uniform
  coef[2][5][0] =  1.133;
  coef[2][5][1] =  0.0;
  coef[2][5][2] =  0.0;
  coef[2][5][3] =  0.0;
  coef[2][5][4] =  0.029;
  
  // .4 Case: Xj of Gumbel to Xi of ?
  //   .1 Gumbel to Normal
  coef[3][0][0] =  1.031;
  //   .2 Gumbel to Lognormal
  coef[3][1][0] =  1.029;
  coef[3][1][1] =  0.014;
  coef[3][1][2] =  0.233;
  coef[3][1][3] =  0.001;
  coef[3][1][4] =  0.004;
  coef[3][1][5] = -0.197;
  //   .3 Gumbel to Exponential
  coef[3][2][0] =  1.142;
  coef[3][2][1] =  0.0;
  coef[3][2][2] =  0.0;
  coef[3][2][3] = -0.154;
  coef[3][2][4] =  0.031;
  //   .4 Gumbel to Gumbel
  coef[3][3][0] =  1.064;
  coef[3][3][1] =  0.0;
  coef[3][3][2] =  0.0;
  coef[3][3][3] = -0.069;
  coef[3][3][4] =  0.005;
  //   .5 Gumbel to Weibull
  coef[3][4][0] =  1.064;
  coef[3][4][1] = -0.210;
  coef[3][4][2] =  0.356;
  coef[3][4][3] =  0.065;
  coef[3][4][4] =  0.003;
  coef[3][4][5] = -0.211;
  //   .6 Gumbel to Uniform
  coef[3][5][0] =  1.055;
  coef[3][5][1] =  0.0;
  coef[3][5][2] =  0.0;
  coef[3][5][3] =  0.0;
  coef[3][5][4] =  0.015;
  
  // .5 Case: Xj of Weibull to Xi of ?
  //   .1 Weibull to Normal
  coef[4][0][0] =  1.031;
  coef[4][0][1] = -0.195;
  coef[4][0][2] =  0.328;
  //   .2 Weibull to Lognormal
  //      !!! difference exists between Ditlevsen & Melchers
  //          the values used here follow those given in Ditlevsen's book.
  coef[4][1][0] =  1.031;
  coef[4][1][1] =  0.011;
  coef[4][1][2] =  0.220;
  coef[4][1][3] =  0.052;
  coef[4][1][4] =  0.002;
  coef[4][1][5] =  0.005;
  coef[4][1][6] = -0.210;
  coef[4][1][7] =  0.350;
  coef[4][1][8] = -0.174;
  coef[4][1][9] =  0.009;
  //   .3 Weibull to Exponential
  //      !!! difference exists between Ditlevsen & Melchers
  //          the values used here follow those given in Ditlevsen's book.
  coef[4][2][0] =  1.147;
  coef[4][2][1] = -0.271;
  coef[4][2][2] =  0.459;
  coef[4][2][3] =  0.145;
  coef[4][2][4] =  0.010;
  coef[4][2][5] = -0.467;
  //   .4 Weibull to Gumbel
  coef[4][3][0] =  1.064;
  coef[4][3][1] = -0.210;
  coef[4][3][2] =  0.356;
  coef[4][3][3] =  0.065;
  coef[4][3][4] =  0.003;
  coef[4][3][5] = -0.211;
  //   .5 Weibull to Weibull
  coef[4][4][0] =  1.063;
  coef[4][4][1] = -0.200;
  coef[4][4][2] =  0.337;
  coef[4][4][3] = -0.004;
  coef[4][4][4] = -0.001;
  coef[4][4][5] =  0.007;
  coef[4][4][6] = -0.200;
  coef[4][4][7] =  0.337;
  coef[4][4][8] =  0.007;
  coef[4][4][9] = -0.007;
  //   .6 Weibull to Uniform
  coef[4][5][0] =  1.061;
  coef[4][5][1] = -0.237;
  coef[4][5][2] =  0.379;
  coef[4][5][3] =  0.0;
  coef[4][5][4] = -0.005;
  
  // .6 Case: Xj of Uniform to Xi of ?
  //   .1 Uniform to Normal
  coef[5][0][0] =  1.023;
  //   .2 Uniform to Lognormal
  coef[5][1][0] =  1.019;
  coef[5][1][1] =  0.014;
  coef[5][1][2] =  0.249;
  coef[5][1][3] =  0.0;
  coef[5][1][4] =  0.01;
  //   .3 Uniform to Exponential
  coef[5][2][0] =  1.133;
  coef[5][2][1] =  0.0;
  coef[5][2][2] =  0.0;
  coef[5][2][3] =  0.0;
  coef[5][2][4] =  0.029;
  //   .4 Uniform to Gumbel
  coef[5][3][0] =  1.055;
  coef[5][3][1] =  0.0;
  coef[5][3][2] =  0.0;
  coef[5][3][3] =  0.0;
  coef[5][3][4] =  0.015;
  //   .5 Uniform to Weibull
  coef[5][4][0] =  1.061;
  coef[5][4][1] = -0.237;
  coef[5][4][2] =  0.379;
  coef[5][4][3] =  0.0;
  coef[5][4][4] = -0.005;
  //   .6 Uniform to Uniform
  coef[5][5][0] =  1.047;
  coef[5][5][1] =  0.0;
  coef[5][5][2] =  0.0;
  coef[5][5][3] =  0.0;
  coef[5][5][4] = -0.047;
  
  // Define the size of the modified correlation matrix
  probdata.mod_correlation.resize(probdata.nvars,probdata.nvars);
  
  int    idisttype, jdisttype; // distribution type index for ith & jth variable
  double R;                    // ratio of correlation coefficient
  double Vi, Vj;               // coefficient of variation of ith & jth variables
  double rho_ji;               // correlation coefficient of jth & ith variables
  double rcoef[10];            // coefficients for R value calculation function
  
  for (int ivar = 0;ivar<probdata.nvars;ivar++) {
    idisttype = probdata.marg(ivar,0);
    
    for (int jvar = ivar+1; jvar<probdata.nvars; jvar++) {
      jdisttype = probdata.marg(jvar,0);
      //cout<<"ith variable = "<<idisttype<<"\t jth variable = "<<jdisttype<<endl;
      
      if ((idisttype == 1) && (jdisttype == 1)) {
        // cout<<"Both are normal distribution"<<endl;
        probdata.mod_correlation(jvar,ivar) = probdata.correlation(jvar,ivar);
        
      }
      else if ((idisttype == 1) && (jdisttype == 2)) {
        // cout<<"Normal and Lognormal"<<endl;
        Vj = probdata.marg(jvar,2);
        R = Vj/sqrt(log(1 + Vj*Vj));
        probdata.mod_correlation(jvar,ivar) = R * probdata.correlation(jvar,ivar);
        
      }
      else if ((idisttype == 2) && (jdisttype == 1)) {
        // cout<<"Lognormal and Normal"<<endl;
        Vj = probdata.marg(jvar,2);
        R = Vj/sqrt(log(1 + Vj*Vj));
        probdata.mod_correlation(jvar,ivar) = R * probdata.correlation(jvar,ivar);
        
      }
      else if ((idisttype == 2) && (jdisttype == 2)) {
        // cout<<"Lognormal and Lognormal"<<endl;
        Vi = probdata.marg(ivar,2);
        Vj = probdata.marg(jvar,2);
        rho_ji = probdata.correlation(jvar,ivar);
        R = log(1 + rho_ji * Vi * Vj) / (sqrt(log(1 + Vj * Vj) * log(1 + Vi * Vi)));
        probdata.mod_correlation(jvar,ivar) = R * probdata.correlation(jvar,ivar);
        
      }
      else {
        // cout<<"Other cases"<<endl;
        
        // Get coefficients from coefficient matrix
        for (int k = 0; k<10; k++) {
          rcoef[k] = coef[jvar][ivar][k];
        }
        //
        Vi = probdata.marg(ivar,2);
        Vj = probdata.marg(jvar,2);
        rho_ji = probdata.correlation(jvar,ivar);
        // calculate ratio R
        R =   rcoef[0]                                        // 1
        + rcoef[6] * Vi + rcoef[7] * Vi * Vi              // Vi & Vi*Vi
        + rcoef[1] * Vj + rcoef[2] * Vj * Vj              // Vj & Vj*Vj
        + rcoef[9] * Vi * Vj                              // Vi*Vj
        + rcoef[3] * rho_ji + rcoef[4] * rho_ji * rho_ji  // rho & rho*rho
        + rcoef[8] * rho_ji * Vi + rcoef[5] * rho_ji * Vj;// rho*Vi & rho*Vj
        probdata.mod_correlation(jvar,ivar) = R * probdata.correlation(jvar,ivar);
        
      }
      //cout<<R<<'\t'<<probdata.mod_correlation(jvar,ivar)<<endl;
    }
  }
  
  // Fill the rest elements of the modified correlation matrix
  for (int irow = 0; irow<probdata.nvars; irow++) {
    for (int jcol = 0; jcol<probdata.nvars; jcol++) {
      if (jcol == irow) {
        probdata.mod_correlation(jcol,irow) = 1.0;
      }
      if ( jcol > irow ) {
        probdata.mod_correlation(irow,jcol) = probdata.mod_correlation(jcol,irow);
      }
      //cout<<probdata.correlation(irow,jcol)<<'\t';
      //cout<<probdata.mod_correlation(irow,jcol)<<'\t';
    }
    //cout<<endl;
  }
}

//------------------------------------------------------------------------------
// To compute Jacobian
//

void Jacobian_u_x(ublas::vector<double> x,
                  ublas::vector<double> u,
                  Stochastic_Model probdata,
                  ublas::matrix<double> &dudx) {
  int nvars;
  nvars = probdata.nvars;
  double imean, istd, iloc, ishape, iscale, ilambda, ilower, iupper;
  int dist_type;
  
  ublas::vector<double> z;
  z = prod(probdata.Lo, u);
  //z = u;
  dudx.resize(nvars,nvars);
  dudx.clear();
  
  ublas::matrix<double> dzdx(nvars,nvars);
  dzdx.clear();
  
  using boost::math::normal_distribution;
  normal_distribution<> snorm(0,1);
  
  for (int i=0; i<nvars; i++) {
    dist_type = probdata.marg(i,0);
    switch (dist_type) {
      case 1: { // Normal distribution
        imean = probdata.marg(i,4);
        istd  = probdata.marg(i,5);
        normal_distribution<> mynorm(imean,istd);
        dzdx(i, i) = pdf(mynorm,x(i))/pdf(snorm,z(i));
        break;
      }
      case 2: { // Lognormal distribution
        using boost::math::lognormal_distribution;
        iloc   = probdata.marg(i,4);
        iscale = probdata.marg(i,5);
        lognormal_distribution<> my_logn(iloc,iscale);
        dzdx(i, i) = pdf(my_logn,x(i))/pdf(snorm,z(i));
        break;
      }
      case 3: { // Exponential distribution
        using boost::math::exponential_distribution;
        ilambda = 1/probdata.marg(i,1);
        exponential_distribution<> my_exp(ilambda);
        dzdx(i, i) = pdf(my_exp,x(i))/pdf(snorm,z(i));
        break;
      }
      case 4: { // Gumbel distribution
        using boost::math::extreme_value_distribution;
        iloc   = probdata.marg(i,4);
        iscale = probdata.marg(i,5);
        extreme_value_distribution<> my_EVI(iloc,iscale);
        dzdx(i, i) = pdf(my_EVI,x(i))/pdf(snorm,z(i));
        break;
      }
      case 5: { // Weibull distribution
        using boost::math::weibull_distribution;
        iloc   = probdata.marg(i,4);
        iscale = probdata.marg(i,5);
        weibull_distribution<> my_wbl(iloc,iscale);
        dzdx(i, i) = pdf(my_wbl,x(i))/pdf(snorm,z(i));
        break;
      }
      case 6: { // Uniform distribution
        using boost::math::uniform_distribution;
        ilower = probdata.marg(i,4);
        iupper = probdata.marg(i,5);
        uniform_distribution<> my_unif(ilower,iupper);
        dzdx(i, i) = pdf(my_unif,x(i))/pdf(snorm,z(i));
        break;
      }
      case 7: { // Gamma distribution
        using boost::math::gamma_distribution;
        ishape = probdata.marg(i,4);
        iscale = probdata.marg(i,5);
        gamma_distribution<> my_gamma(ishape,iscale);
        dzdx(i, i) = pdf(my_gamma,x(i))/pdf(snorm,z(i));
        break;
      }
        
      default:
        break;
    }
  }
  dudx = prod(probdata.inv_Lo, dzdx);
}

//------------------------------------------------------------------------------
// To compute limit state function and its gradient w.r.t. the given trail
// checking point
//

void gfun(Stochastic_Model probdata,
          ublas::vector<double> x,
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
  
  
  // ================================
  //
  // Example 1
  //
  // ================================
  /*
  // Evaluate the limit state function
  val_lsf = x(0)*x(1)-78.12*x(2);
  
  // Evaluate the partial derovative of the LSF w.r.t. basic variables
  // w.r.t. the 1st random variable
  grad_lsf(0) = x(1);
  // w.r.t. the 2nd random variable
  grad_lsf(1) = x(0);
  // w.r.t. the 3rd random variable
  grad_lsf(2) = -78.12;
  */
  
  /*
  val_lsf = x(0)*x(1)-600*x(2);
  
  grad_lsf(0) = x(1);
  grad_lsf(1) = x(0);
  grad_lsf(2) = -600;
  */
  
  
  // ================================
  //
  // Example 3
  //
  // ================================
  // Evaluate the limit state function
  val_lsf = pow(x(0),4) + 2*pow(x(1),4) - 20;
  
  // Evaluate the 1st-order partial derivatives of the LSF w.r.t. basic variables
  // w.r.t. the 1st random variable
  grad_lsf(0) = 4*pow(x(0),3);
  // w.r.t. the 2nd random variable
  grad_lsf(1) = 8*pow(x(1),3);
  
  // Evaluate the 2nd-order partial derivatives of the LSF w.r.t. basic variables
  // w.r.t. the 1st random variable
  //Hess_lsf(0,0) = 12*pow(x(0),2);
  // w.r.t. the 2nd random variable
  //Hess_lsf(1,1) = 24*pow(x(1),2);
  
  
  /*
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
  Hess_lsf(0,0) = a*(-pow(x(1),2));
  Hess_lsf(0,1) = a*(-x(0)*x(1));
  Hess_lsf(0,2) = 0;
  // w.r.t. the 2nd random variable
  Hess_lsf(1,0) = a*(-x(0)*x(1));
  Hess_lsf(1,1) = a*(-pow(x(0),2));
  Hess_lsf(1,2) = 0;
  // w.r.t. the 3rd random variable
  Hess_lsf(2,0) = 0;
  Hess_lsf(2,1) = 0;
  Hess_lsf(2,2) = 0;
  */
  
  
}

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




/*******************************************************************************
 *                                                                             *
 *             MAIN CODE FOR STRUCTURAL RELIABILITY ANLYSIS                    *
 *                                                                             *
/*******************************************************************************/



int main(int argc, char * argv[]) {
  
  clock_t start_time, finish_time;
  double total_time;
  start_time = clock();
  
  ErrorCode rval;
  PetscErrorCode ierr;
  PetscInitialize(&argc,&argv,(char *)0,help);
  PetscBool flg = PETSC_TRUE;
  const char *option;
  PetscInt order;
  
//  mt19937 rng;                         // using pseudo-random generator: mt19937
//  double myrnd;
//  normal_distribution<> normdist(0,1);        // setting normal distribution rnd
//  variate_generator <mt19937&,normal_distribution<> > normrnd(rng,normdist);
//  // combination of distribution and generator
//  myrnd = normrnd(); cout<<myrnd<<endl;
//  double value;
//  value = ltqnorm(myrnd);
//  cout<<value<<endl;
  
  
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
  
  ReliabOpt.echo_flag  = 1;       // Program interactive mode, 0: silent mode
  // FORM analysis options
  ReliabOpt.istep_max  = 2000;     // Maximum number of interation allowed in the search algorithm
  ReliabOpt.e1         = 0.001;   // Tolerance on how close design point is to limit-state surface
  ReliabOpt.e2         = 0.001;   // Tolerance on how accurately the gradient points towards the origin
  ReliabOpt.step_code  = 0.025;//0.025;   // 0: step size by Armijo rule, otherwise: given value is the step size
  ReliabOpt.Recorded_u = 1;       // 0: u-vector not recorded at all iterations, 1: u-vector recorded at all iterations
  ReliabOpt.Recorded_x = 1;       // 0: x-vector not recorded at all iterations, 1: x-vector recorded at all iterations
  ReliabOpt.grad_G     = "PSFEM"; // "PSFEM": perturbation based SFEM, "DDM": direct differentiation, 'ADM': automatic differentiation
  
  
  int    echo_flag = ReliabOpt.echo_flag;
  double e1        = ReliabOpt.e1;
  double e2        = ReliabOpt.e2;
  int    istep_max = ReliabOpt.istep_max;
  double step_code = ReliabOpt.step_code;
  
  /*
   *  Read inputs' statistical properties from file to insert into <probdata>
   */
  Stochastic_Model probdata;
  cout<<"\n\nRead file"<<endl;
  ifstream ProbDataFile("/Users/nxz6/Dropbox/DURACOMP_Cal/009_MoFEM/04_ReliabilityAnalysis/Input_probdata_Example01.txt");
  cout<<"\n Open the file"<<endl;
  
  char   buffer[256];
  string stringbuf;
  string substringbuf;
  string datatype;
  int    cnt;
  int    *pos;
  int    MAR_IX, COR_IX;
  MAR_IX = 0; COR_IX = 0;
  while (!ProbDataFile.eof()){
    ProbDataFile.getline(buffer,100);
    if (strlen(buffer)>0) {
      //stringbuf = (string)buffer;
      //substringbuf = stringbuf.substr(0,1);
      //cout<<substringbuf<<'\t'<<atof(substringbuf.c_str())<<endl;
      if (isdigit(buffer[0]) == 0) {
        stringbuf = (string)buffer;
        if (stringbuf.compare(0,3,"NUM") == 0) {
          // cout<<"Next line is data for number of variables"<<endl;
          datatype = "NUMBER";
        }
        else if (stringbuf.compare(0,3,"COR") == 0) {
          // cout<<"Next line is data for correlation matrix"<<endl;
          datatype = "CORRELATION";
        }
        else if (stringbuf.compare(0,3,"MAR") == 0) {
          // cout<<"Next line is data for marginal distribution"<<endl;
          datatype = "MARGINAL";
        }
        stringbuf.clear();
      }
      else {
        str_cnt(buffer,',',cnt);
        pos = new int[cnt];
        str_pos(buffer,',',pos);
        stringbuf = (string)buffer;
        
        if (datatype.compare(0,3,"NUM") == 0) { // number of variables
          substringbuf = stringbuf.substr(0,pos[1]);
          probdata.nvars = atoi(substringbuf.c_str());
        }
        else if (datatype.compare(0,3,"MAR") == 0) { // marginal distribution
          // Declaration
          if (MAR_IX ==0) {
            probdata.marg.resize(probdata.nvars,10);
          }
          MAR_IX ++;
          
          // Insert data into probdata
          for (int i=1; i<=cnt;i++) {
            if (i==1) {
              substringbuf = stringbuf.substr(0,pos[i]);
            }
            else {
              substringbuf = stringbuf.substr(pos[i-1]+1,(pos[i]-pos[i-1])-1);
            }
            probdata.marg(MAR_IX-1,i-1) = atof(substringbuf.c_str());
          }
        }
        else if (datatype.compare(0,3,"COR") == 0) { // correlation matrix
          // Declaration
          if (COR_IX ==0) {
            probdata.correlation.resize(probdata.nvars,probdata.nvars);
          }
          COR_IX ++;
          
          // Insert data into probdata
          for (int i=1; i<=cnt;i++) {
            if (i == 1) {
              substringbuf = stringbuf.substr(0,pos[i]);
            }
            else {
              substringbuf = stringbuf.substr(pos[i-1]+1,(pos[i]-pos[i-1])-1);
            }
            probdata.correlation(COR_IX-1, i-1) = atof(substringbuf.c_str());
          }
        }
        // free dynamically allocated memory
        delete [] pos;
        // cout<<datatype<<'\t'<<cnt<<'\t'<<buffer<<endl;
      }
    }
  }
  

  cout<<"number of variables: "<<probdata.nvars<<endl;
  //cout<<"distribution type: "<<probdata.marg[0][0]<<endl;
  //cout<<"correlation coefficient: "<<probdata.correlation[1][1]<<endl;

  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 2: use probability transformation to obtain y in             //
  //                standard normal distribution space                        //
  //                y = norminv(F(x))                                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  modify_correlation_mat(probdata);
  
  /*
   * Perform Cholesky decomposition for the modified correlation matrix
   *    A = LL'
   */
  ublas::triangular_matrix<double, ublas::lower> Lmat(probdata.nvars,probdata.nvars);
  cholesky_decompose(probdata.mod_correlation,Lmat);
  probdata.Lo.resize(probdata.nvars,probdata.nvars);
  probdata.Lo = Lmat;
  
  //cout<<"Cholesky decomposed matrix: "<<Lmat<<endl;
  
  
  // Compute the inverse of Lo
  probdata.inv_Lo.resize(probdata.nvars,probdata.nvars);
  bool singular = false;
  probdata.inv_Lo = gjinverse(probdata.Lo,singular);
  //cout<<"Lmat \t"<<probdata.Lo<<endl;
  //cout<<"inverse of Lmat \t"<<probdata.inv_Lo<<endl;
  
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 3: select an initial checking point, x                       //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////

  ublas::vector<double> x(probdata.nvars);
  for (int i=0; i<probdata.nvars; i++) {
    x(i) = probdata.marg(i,3);
  }
  double detj;
  ublas::vector<double> u(probdata.nvars);
  x_to_u(x,probdata,u);
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 4: Perform iterative loop to find design point               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  /*
   *  Set parameters for the iterative loop
   */
  int    istep = 1;     // Initialize iterative counter
  int    conv_flag = 0; // Convergence is achieved when this flag is set to
  
  double val_G, val_G0;                         // LSF for given inputs
  double step_size;                             // Step size
  ublas::vector<double> grad_g(probdata.nvars); // Gradient of LSF in x space
  ublas::vector<double> grad_G(probdata.nvars); // Gradient of LSF in u space
  ublas::vector<double> alpha(probdata.nvars);  // Direction cosine vector
  ublas::matrix<double> dudx(probdata.nvars,probdata.nvars);
  ublas::matrix<double> inv_dudx(probdata.nvars,probdata.nvars);
  ublas::vector<double> u_dir(probdata.nvars);  // Direction
  ublas::vector<double> u_new(probdata.nvars);  // New trial of checking point
  
  /*
   *  Start iteration
   */
  
  do {
    
    if (echo_flag) {
      cout<<"-------------------------------- \n";
      cout<<"Now carrying out iteration number: \t"<<istep<<endl;
    }
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //        STEP 5: calculate Jacobian, J = dy/dx                           //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    /*
     * Transformation from u to x space
     */
    x.resize(probdata.nvars); x.clear();
    u_to_x(u,probdata,x,detj);
    
    // Jacobian
    dudx.resize(probdata.nvars,probdata.nvars); dudx.clear();
    Jacobian_u_x(x,u,probdata,dudx);
    
    inv_dudx.resize(probdata.nvars,probdata.nvars);
    inv_dudx.clear();
    inv_dudx = gjinverse(dudx,singular);
    
    /*
     * Evaluate limit-state function and its gradient
     */
    grad_g.resize(probdata.nvars); grad_g.clear();
    gfun(probdata,x,val_G,grad_g);
    
    grad_G.resize(probdata.nvars); grad_G.clear();
    cout<<dudx<<"\t"<<inv_dudx<<endl;
    grad_G = prod(grad_g,inv_dudx);
    
    cout<<"LSF value is \t"<<val_G<<endl;
    cout<<"Gradient of LSF at normal space and original space are \t"<<grad_g<<"\t"<<grad_G<<endl;
    
    /*
     *  Set scale parameter G0 and inform about structural response
     */
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
    /*
     * Compute direction cosine alpha vector
     */
    alpha.resize(probdata.nvars); alpha.clear();
    alpha = -grad_G/norm_2(grad_G);
    cout<<"Direction cosine "<<alpha<<endl;
    
    /*
     * Check convergence
     */
    if (((abs(val_G/val_G0)<e1) && (norm_2(u-inner_prod(alpha,u)*alpha))) || (istep == istep_max)) {
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
    
    /*
     * Take a step if convergence is not achieved
     */
    if (conv_flag == 0) {
      // Determine search direction
      search_dir(val_G,grad_G,u,u_dir);
      // Determine step size
      if (step_code != 0) {
        step_size = step_code;
      }
      else {
        cout<<"\n!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!";
        cout<<"\nArmijo rule will be used to determine step size for setting new trial point, \n";
        cout<<"but it isn't available yet in current version!!! \n";
        cout<<"Re-assign a nonzero value to step_code.\n";
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        conv_flag = 1;
      }
      
      // Determin new trial point
      u_new.resize(probdata.nvars); u_new.clear();
      u_new = u + step_size*u_dir;
      
      cout<<"\n\nReliability index estimate at "<<istep<<" is: \t"<<((val_G/norm_2(grad_G)) + inner_prod(alpha, u))<<endl;
      
      // Prepare for a new round in the loop
      u.resize(probdata.nvars); u.clear();
      u = u_new;
      istep = istep + 1;
    }
    
    
  } while (conv_flag == 0); // end of while
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 8: Calculate the estimate of reliability index               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////

  beta = inner_prod(alpha,u);
  
  
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
  }
  
  cout<<"*  Optimal reliability index (beta) is: "<<beta<<endl;
  cout<<"*  Elapsed time is "<<total_time<<" seconds.\n";
  cout<<"*  The program finishes !!! \n*\n";
  cout<<"*************************************************"<<endl;
  
  

  ierr = PetscFinalize(); CHKERRQ(ierr);
}
