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

struct MCS_RNG {
  
  
  //------------------------------------------------------------------------------
  // To generate sample
  //
  
  virtual PetscErrorCode myRNG(int num_vars,
                               ublas::matrix<double> MargProb,
                               ublas::vector<double> &x) {
    PetscFunctionBegin;
    
    ErrorCode rval;
    PetscErrorCode ierr;
    
    x.resize(num_vars);x.clear();
    
    boost::mt19937 rng(time(0));        // using pseudo-random generator: mt19937
    
    double imean, istd, iloc, ishape, iscale, ilambda, ilower, iupper;
    int dist_type;
    double myrnd;
    
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
          boost::random::normal_distribution<> my_norm(imean,istd);
          boost::variate_generator <boost::mt19937&,boost::random::normal_distribution<> > normrnd(rng,my_norm);
          x(i) = normrnd();
          break;
        }
        case 2: {
          /*
           * distribution type 2: lognormal
           *           parameter: location - mean of logrithm rv
           *                      scale    - std of logrithm rv
           */
          iloc   = MargProb(i,4);
          iscale = MargProb(i,5);
          boost::random::lognormal_distribution<> my_logn(iloc,iscale);
          boost::variate_generator <boost::mt19937&,boost::random::lognormal_distribution<> > lognrnd(rng,my_logn);
          x(i) = lognrnd();
          break;
        }
        case 3: {
          /*
           * distribution type 3: exponential
           *           parameter: lambda = 1/mu
           */
          ilambda = 1/MargProb(i,1);
          boost::random::exponential_distribution<> my_exp(ilambda);
          boost::variate_generator <boost::mt19937&,boost::random::exponential_distribution<> > exprnd(rng,my_exp);
          x(i) = exprnd();
          break;
        }
        case 4: {
          /*
           * distribution type 4: Extreme value distribution - I (Gumbel)
           *           parameter: shape
           *                      location
           *                      scale
           */
          iloc   = MargProb(i,4);
          iscale = MargProb(i,5);
          boost::random::extreme_value_distribution<> my_EVI(iloc,iscale);
          boost::variate_generator <boost::mt19937&,boost::random::extreme_value_distribution<> > evrnd(rng,my_EVI);
          x(i) = evrnd();
          break;
        }
        case 5: {
          /*
           * distribution type 5: Extreme value distribution - III (Weibull)
           *           parameter: shape
           *                      scale
           */
          iloc   = MargProb(i,4);
          iscale = MargProb(i,5);
          boost::random::weibull_distribution<> my_wbl(iloc,iscale);
          boost::variate_generator <boost::mt19937&,boost::random::weibull_distribution<> > wblrnd(rng,my_wbl);
          x(i) = wblrnd();
          break;
        }
        case 6: {
          /*
           * distribution type 6: Uniform distribution
           *           parameter: shape
           *                      scale
           */
          ilower = MargProb(i,4);
          iupper = MargProb(i,5);
          boost::uniform_real<> my_unif(ilower,iupper);
          boost::variate_generator <boost::mt19937&,boost::uniform_real<> > unifrnd(rng,my_unif);
          x(i) = unifrnd();
          break;
        }
        case 7: {
          /*
           * distribution type 7: Gamma distribution
           *           parameter: shape
           *                      scale
           */
          ishape = MargProb(i,4);
          iscale = MargProb(i,5);
          boost::random::gamma_distribution<> my_gamma(ishape,iscale);
          boost::variate_generator <boost::mt19937&,boost::random::gamma_distribution<> > gamrnd(rng,my_gamma);
          x(i) = gamrnd();
          break;
        }
        default:
          cout<<"The distribution type index should be value between 0 and 7!"<<endl;
          break;
      }
      
      cout<<i<<"th variable \t"<<x(i)<<endl;
    }
    
    PetscFunctionReturn(0);
  }
  
};

