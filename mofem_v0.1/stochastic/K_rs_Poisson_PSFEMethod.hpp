/* Copyright (C) 2014, Xiao-Yi Zhou (xiaoyi.zhou@ncl.ac.uk)
 * --------------------------------------------------------------
 * This routine calculates the second-order derivative of stiffness matrix,
 * K_rs, with respect to Poisson's ratio for isotropic material for the 
 * stochastic finite element
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

#ifndef __K_RS_POISSON_PSFEMETHOD_HPP__
#define __K_RS_POISSON_PSFEMETHOD_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
extern "C" {
#include <gm_rule.h>
}

namespace MoFEM {

struct K_rs_Poisson_PSFEMethod: public K_rs_PSFEMethod {

    K_rs_Poisson_PSFEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young,double _pois,const string& _second_field):
    K_rs_PSFEMethod( _mField,_Aij,_X,_F,_young,_pois,_second_field){
      
    }
    //F: calculate the first-order derivative of constitutive matrix  
    //   with respect to Poisson's ratio
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
    //   with respect to Poisson's ratio
    virtual PetscErrorCode calculateD_rs(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear(); // Initiate D_rs to zero

      double D_rs00, D_rs01, D_rs33;

      D_rs00 = -(4*young*(- 2*nu*nu*nu + 6*nu*nu + 1))/pow((2*nu*nu + nu - 1),3);
      D_rs01 = -(2*young*(4*nu*nu*nu + 6*nu + 1))/pow((2*nu*nu + nu - 1),3);
      D_rs33 = young/pow((nu + 1),3);
      
      
       // Set nonzero values
      D(0,0) = D_rs00; D(0,1) = D_rs01;  D(0,2) = D_rs01; 
      D(1,0) = D_rs01; D(1,1) = D_rs00;  D(1,2) = D_rs01; 
      D(2,0) = D_rs01; D(2,1) = D_rs01;  D(2,2) = D_rs00;
      D(3,3) = D_rs33;
      D(4,4) = D_rs33;
      D(5,5) = D_rs33;

      // cout<<"D_rs = "<<D;

      PetscFunctionReturn(0);
    }
    
};

    
}

#endif //__K_RS_POISSON_PSFEMETHOD_HPP__
