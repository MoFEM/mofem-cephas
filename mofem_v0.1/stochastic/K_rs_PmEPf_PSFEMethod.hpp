/* Copyright (C) 2014, Xiao-Yi Zhou (xiaoyi.zhou@ncl.ac.uk)
 * --------------------------------------------------------------
 * This routine calculates the second-order partial derivative of stiffness 
 * matrix, K_rs, with respect to Poisson's ratio of matrix and Young's modulus 
 * or Poisson's ratio of inclusion/fibre for two-phase material for the 
 * stochastic finite element. 
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

#ifndef __K_RS_PMEPF_PSFEMETHOD_HPP__
#define __K_RS_PMEPF_PSFEMETHOD_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "K_rsPoissonFEMethod.hpp"

extern "C" {
#include <gm_rule.h>
}

namespace MoFEM {

struct K_rs_PmEPf_PSFEMethod: public K_rs_PSFEMethod {

    K_rs_PmEPf_PSFEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young,double _pois,const string& _second_field):
      K_rs_PSFEMethod( _mField,_Aij,_X,_F,_young,_pois,_second_field){
        
    }
    //F: calculate the first-order derivative of constitutive matrix  
    //   with respect to Poisson's ratio of matrix
    virtual PetscErrorCode calculateD_r(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear();

      double D_r00,D_r01,D_r33;
      D_r00 = -(2*young*nu*(nu - 2))/pow((2*nu*nu + nu - 1),2);
      D_r01 = (young*(2*nu*nu + 1))/pow((2*nu*nu + nu - 1),2);
      D_r33 = -young/(2*pow((nu + 1),2));
      
      D(0,0)=D_r00;  D(0,1)=D_r01;   D(0,2)=D_r01;
      D(1,0)=D_r01;  D(1,1)=D_r00;   D(1,2)=D_r01;
      D(2,0)=D_r01;  D(2,1)=D_r01;   D(2,2)=D_r00;
      D(3,3)=D_r33;
      D(4,4)=D_r33;
      D(5,5)=D_r33;
      // cout<<"D_r = "<<D;
      PetscFunctionReturn(0);
    }

    //F: calculate the second-order derivative of constitutive matrix 
    //   with respect to Poisson's ratio of matrix and Young's modulus of fibre
    virtual PetscErrorCode calculateD_rs(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear(); // Initiate D_rs to zero

      // cout<<"D_rs = "<<D;

      PetscFunctionReturn(0);
    }
  
};

    
}

#endif //__K_RS_PMEPF_PSFEMETHOD_HPP__
