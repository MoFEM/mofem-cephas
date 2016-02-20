/* \brief ModifyCorrelationMatrix
 *
 * Edited and modified by Xiaoyi Zhou.
 *
 * This routine iss used to modify correlation matrix when transforming non-gaussian
 *   random variables to Gaussian.
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

struct NatafTransformation {
  
  
  //------------------------------------------------------------------------------
  // To transform x from original space to standard normal distribution probability
  // space using Nataf transformation
  //
  
  virtual PetscErrorCode x_to_u(ublas::vector<double> x,
                                int num_vars,
                                ublas::matrix<double> MargProb,
                                ublas::matrix<double> inv_Lo,
                                ublas::vector<double> &u) {
    PetscFunctionBegin;
    
    ublas::vector<double> z(num_vars);
    u.resize(num_vars);u.clear();
    
    using boost::math::normal_distribution;
    normal_distribution<> snorm(0,1);
    double imean, istd, iloc, ishape, iscale, ilambda, ilower, iupper;
    int dist_type;
    
    for (unsigned i=0; i < num_vars; i++) {
      dist_type = MargProb(i,0);
      switch (dist_type) { // distribution type
        case 1: {
          /*
           * distribution type 1: Normal
           *           parameter: location - mean
           *                      scale    - standard deviation
           */
          imean = MargProb(i,4);
          istd  = MargProb(i,5);
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
          iloc   = MargProb(i,4);
          iscale = MargProb(i,5);
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
          ilambda = 1/MargProb(i,1);
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
          iloc   = MargProb(i,4);
          iscale = MargProb(i,5);
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
          iloc   = MargProb(i,4);
          iscale = MargProb(i,5);
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
          ilower = MargProb(i,4);
          iupper = MargProb(i,5);
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
          ishape = MargProb(i,4);
          iscale = MargProb(i,5);
          gamma_distribution<> my_gamma(ishape,iscale);
          z(i) = quantile(snorm,cdf(my_gamma,x(i)));
          break;
        }
        default:
          cout<<"The distribution type index should be value between 0 and 7!"<<endl;
          break;
      }
      
      /*if (i == 0) { cout<<"Value of 1st variable \t"<<z(i)<<endl; }
      else if (i == 1) { cout<<"Value of 2nd variable \t"<<z(i)<<endl; }
      else if (i == 2) { cout<<"Value of 3rd variable \t"<<z(i)<<endl; }
      else { cout<<"Value of "<<i+1<<"th variable \t"<<z(i)<<endl;}*/
    }
    u = prod(inv_Lo,z); // u = Lo^{-1} z
    
    PetscFunctionReturn(0);
  }
  
  //------------------------------------------------------------------------------
  // To transform u from standard normal distribution probability
  // space to original space to obtain x using Nataf transformation
  //
  
  virtual PetscErrorCode u_to_x(ublas::vector<double> u,
                                int num_vars,
                                ublas::matrix<double> MargProb,
                                ublas::matrix<double> Lo,
                                ublas::vector<double> &x,
                                double &detj) {
    PetscFunctionBegin;
    
    x.resize(num_vars);
    ublas::vector<double> z(num_vars);
    z = prod(Lo,u);
    
    using boost::math::normal_distribution;
    normal_distribution<> snorm(0,1);
    double imean, istd, iloc, ishape, iscale, ilambda, ilower, iupper;
    int dist_type;
    detj = 1.0;
    
    for (int i=0;i<num_vars;i++) {
      dist_type = MargProb(i,0);
      switch (dist_type) { // distribution type
        case 1: {
          /*
           * distribution type 1: Normal
           *           parameter: location - mean
           *                      scale    - standard deviation
           */
          imean = MargProb(i,4);
          istd  = MargProb(i,5);
          normal_distribution<> mynorm(imean,istd);
          x(i) = quantile(mynorm,cdf(snorm,z(i)));
          detj = detj*cdf(mynorm,x(i))/cdf(snorm,u(i));
          break;
        }
        case 2: {
          /*
           * distribution type 2: lognormal
           *           parameter: location - mean of logrithm rv
           *                      scale    - std of logrithm rv
           */
          using boost::math::lognormal_distribution;
          iloc   = MargProb(i,4);
          iscale = MargProb(i,5);
          lognormal_distribution<> my_logn(iloc,iscale);
          x(i) = quantile(my_logn,cdf(snorm,z(i)));
          detj = detj*cdf(my_logn,x(i))/cdf(snorm,u(i));
          break;
        }
        case 3: {
          /*
           * distribution type 3: exponential
           *           parameter: lambda = 1/mu
           */
          using boost::math::exponential_distribution;
          ilambda = 1/MargProb(i,1);
          exponential_distribution<> my_exp(ilambda);
          x(i) = quantile(my_exp,cdf(snorm,z(i)));
          detj = detj*cdf(my_exp,x(i))/cdf(snorm,u(i));
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
          iloc   = MargProb(i,4);
          iscale = MargProb(i,5);
          extreme_value_distribution<> my_EVI(iloc,iscale);
          x(i) = quantile(my_EVI,cdf(snorm,z(i)));
          detj = detj*cdf(my_EVI,x(i))/cdf(snorm,u(i));
          break;
        }
        case 5: {
          /*
           * distribution type 5: Extreme value distribution - III (Weibull)
           *           parameter: shape
           *                      scale
           */
          using boost::math::weibull_distribution;
          iloc   = MargProb(i,4);
          iscale = MargProb(i,5);
          weibull_distribution<> my_wbl(iloc,iscale);
          x(i) = quantile(my_wbl,cdf(snorm,z(i)));
          detj = detj*cdf(my_wbl,x(i))/cdf(snorm,u(i));
          break;
        }
        case 6: {
          /*
           * distribution type 6: Uniform distribution
           *           parameter: shape
           *                      scale
           */
          using boost::math::uniform_distribution;
          ilower = MargProb(i,4);
          iupper = MargProb(i,5);
          uniform_distribution<> my_unif(ilower,iupper);
          x(i) = quantile(my_unif,cdf(snorm,z(i)));
          detj = detj*cdf(my_unif,x(i))/cdf(snorm,u(i));
          break;
        }
        case 7: {
          /*
           * distribution type 7: Gamma distribution
           *           parameter: shape
           *                      scale
           */
          using boost::math::gamma_distribution;
          ishape = MargProb(i,4);
          iscale = MargProb(i,5);
          gamma_distribution<> my_gamma(ishape,iscale);
          x(i) = quantile(my_gamma,cdf(snorm,z(i)));
          detj = detj*cdf(my_gamma,x(i))/cdf(snorm,u(i));
          break;
        }
        default:
          cout<<"The distribution type index should be value between 0 and 7!"<<endl;
          break;
      }
      
      /*if (i == 0) { cout<<"Value of 1st variable \t"<<x(i)<<"\t"<<u(i)<<endl; }
      else if (i == 1) { cout<<"Value of 2nd variable \t"<<x(i)<<"\t"<<u(i)<<endl; }
      else if (i == 2) { cout<<"Value of 3rd variable \t"<<x(i)<<"\t"<<u(i)<<endl; }
      else { cout<<"Value of "<<i+1<<"th variable \t"<<x(i)<<"\t"<<u(i)<<endl;}*/
    }
    
    PetscFunctionReturn(0);
  }
  
  //------------------------------------------------------------------------------
  // To compute Jacobian
  //
  
  virtual PetscErrorCode Jacobian_u_x(ublas::vector<double> x,
                                      ublas::vector<double> u,
                                      int num_vars,
                                      ublas::matrix<double> MargProb,
                                      ublas::matrix<double> Lo,
                                      ublas::matrix<double> inv_Lo,
                                      ublas::matrix<double> &dudx) {
    PetscFunctionBegin;
    
    double imean, istd, iloc, ishape, iscale, ilambda, ilower, iupper;
    int dist_type;
    
    ublas::vector<double> z;
    z = prod(Lo, u);
    dudx.resize(num_vars,num_vars); dudx.clear();
    
    ublas::matrix<double> dzdx(num_vars,num_vars); dzdx.clear();
    
    using boost::math::normal_distribution;
    normal_distribution<> snorm(0,1);
    
    for (int i=0; i<num_vars; i++) {
      dist_type = MargProb(i,0);
      switch (dist_type) {
        case 1: { // Normal distribution
          imean = MargProb(i,4);
          istd  = MargProb(i,5);
          normal_distribution<> mynorm(imean,istd);
          dzdx(i, i) = pdf(mynorm,x(i))/pdf(snorm,z(i));
          break;
        }
        case 2: { // Lognormal distribution
          using boost::math::lognormal_distribution;
          iloc   = MargProb(i,4);
          iscale = MargProb(i,5);
          lognormal_distribution<> my_logn(iloc,iscale);
          dzdx(i, i) = pdf(my_logn,x(i))/pdf(snorm,z(i));
          break;
        }
        case 3: { // Exponential distribution
          using boost::math::exponential_distribution;
          ilambda = 1/MargProb(i,1);
          exponential_distribution<> my_exp(ilambda);
          dzdx(i, i) = pdf(my_exp,x(i))/pdf(snorm,z(i));
          break;
        }
        case 4: { // Gumbel distribution
          using boost::math::extreme_value_distribution;
          iloc   = MargProb(i,4);
          iscale = MargProb(i,5);
          extreme_value_distribution<> my_EVI(iloc,iscale);
          dzdx(i, i) = pdf(my_EVI,x(i))/pdf(snorm,z(i));
          break;
        }
        case 5: { // Weibull distribution
          using boost::math::weibull_distribution;
          iloc   = MargProb(i,4);
          iscale = MargProb(i,5);
          weibull_distribution<> my_wbl(iloc,iscale);
          dzdx(i, i) = pdf(my_wbl,x(i))/pdf(snorm,z(i));
          break;
        }
        case 6: { // Uniform distribution
          using boost::math::uniform_distribution;
          ilower = MargProb(i,4);
          iupper = MargProb(i,5);
          uniform_distribution<> my_unif(ilower,iupper);
          dzdx(i, i) = pdf(my_unif,x(i))/pdf(snorm,z(i));
          break;
        }
        case 7: { // Gamma distribution
          using boost::math::gamma_distribution;
          ishape = MargProb(i,4);
          iscale = MargProb(i,5);
          gamma_distribution<> my_gamma(ishape,iscale);
          dzdx(i, i) = pdf(my_gamma,x(i))/pdf(snorm,z(i));
          break;
        }
          
        default:
          break;
      }
    }
    
    dudx = prod(inv_Lo, dzdx);
    
    PetscFunctionReturn(0);
  }
  
  //------------------------------------------------------------------------------
  // To compute modified correlation matrix by using Nataf transformation
  //
  
  virtual PetscErrorCode ModCorrMat_Empirical(int num_vars,
                                              ublas::matrix<double> MargProb,
                                              ublas::matrix<double> CorrMat,
                                              ublas::matrix<double> &ModifiedCorrMat) {
    PetscFunctionBegin;
    
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
    ModifiedCorrMat.resize(num_vars,num_vars);
    
    int    idisttype, jdisttype; // distribution type index for ith & jth variable
    double R;                    // ratio of correlation coefficient
    double Vi, Vj;               // coefficient of variation of ith & jth variables
    double rho_ji;               // correlation coefficient of jth & ith variables
    double rcoef[10];            // coefficients for R value calculation function
    
    for (int ivar = 0; ivar < num_vars; ivar++) {
      idisttype = MargProb(ivar,0);
      
      for (int jvar = ivar+1; jvar<num_vars; jvar++) {
        jdisttype = MargProb(jvar,0);
        //cout<<"ith variable = "<<idisttype<<"\t jth variable = "<<jdisttype<<endl;
        
        if ((idisttype == 1) && (jdisttype == 1)) {
          // cout<<"Both are normal distribution"<<endl;
          ModifiedCorrMat(jvar,ivar) = CorrMat(jvar,ivar);
          
        }
        else if ((idisttype == 1) && (jdisttype == 2)) {
          // cout<<"Normal and Lognormal"<<endl;
          Vj = MargProb(jvar,2);
          R = Vj/sqrt(log(1 + Vj*Vj));
          ModifiedCorrMat(jvar,ivar) = R * CorrMat(jvar,ivar);
          
        }
        else if ((idisttype == 2) && (jdisttype == 1)) {
          // cout<<"Lognormal and Normal"<<endl;
          Vj = MargProb(jvar,2);
          R = Vj/sqrt(log(1 + Vj*Vj));
          ModifiedCorrMat(jvar,ivar) = R * CorrMat(jvar,ivar);
          
        }
        else if ((idisttype == 2) && (jdisttype == 2)) {
          // cout<<"Lognormal and Lognormal"<<endl;
          Vi = MargProb(ivar,2);
          Vj = MargProb(jvar,2);
          rho_ji = CorrMat(jvar,ivar);
          R = log(1 + rho_ji * Vi * Vj) / (sqrt(log(1 + Vj * Vj) * log(1 + Vi * Vi)));
          ModifiedCorrMat(jvar,ivar) = R * CorrMat(jvar,ivar);
          
        }
        else {
          //cout<<"Other cases"<<endl;
          //cout<<"Number of variables "<<num_vars<<endl;
          // Get coefficients from coefficient matrix
          for (int k = 0; k<10; k++) {
            rcoef[k] = coef[jdisttype][idisttype][k];
          }
          //
          Vi = MargProb(ivar,2);
          Vj = MargProb(jvar,2);
          rho_ji = CorrMat(jvar,ivar);
          // calculate ratio R
          R =   rcoef[0]                                        // 1
          + rcoef[6] * Vi + rcoef[7] * Vi * Vi              // Vi & Vi*Vi
          + rcoef[1] * Vj + rcoef[2] * Vj * Vj              // Vj & Vj*Vj
          + rcoef[9] * Vi * Vj                              // Vi*Vj
          + rcoef[3] * rho_ji + rcoef[4] * rho_ji * rho_ji  // rho & rho*rho
          + rcoef[8] * rho_ji * Vi + rcoef[5] * rho_ji * Vj;// rho*Vi & rho*Vj
          ModifiedCorrMat(jvar,ivar) = R * CorrMat(jvar,ivar);
          
        }
        //cout<<R<<'\t'<<ModifiedCorrMat(jvar,ivar)<<endl;
      }
    }
    
    // Fill the rest elements of the modified correlation matrix
    for (int irow = 0; irow<num_vars; irow++) {
      for (int jcol = 0; jcol<num_vars; jcol++) {
        if (jcol == irow) {
          ModifiedCorrMat(jcol,irow) = 1.0;
        }
        if ( jcol > irow ) {
          ModifiedCorrMat(irow,jcol) = ModifiedCorrMat(jcol,irow);
        }
        //cout<<CorrMat(irow,jcol)<<'\t';
        //cout<<ModifiedCorrMat(irow,jcol)<<'\t';
      }
      //cout<<endl;
    }
    //cout<<"\n\nModified correlation matrix: "<<ModifiedCorrMat<<endl;
    PetscFunctionReturn(0);
  }
  
};

