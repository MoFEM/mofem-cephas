/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah@glasgow.ac.uk)
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

#ifndef __K_RYOUNGSPOISSONFEMETHOD_HPP__
#define __K_RYOUNGSPOISSONFEMETHOD_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
extern "C" {
#include <gm_rule.h>
}

namespace MoFEM {

struct K_rYoungPoissonFEMethod: public K_rsPoissonFEMethod {

    K_rYoungPoissonFEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young,double _pois,const string& _second_field):
      K_rsPoissonFEMethod( _mField,_Aij,_X,_F,_young,_pois,_second_field){
        
    }
    //F: calculate the first-order derivative of constitutive matrix  
    //   with respect to Poisson's ratio
    virtual PetscErrorCode calculateD_r(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear();
      double D_r00,D_r01,D_r33,constt;
      constt=1/((1+nu)*(1-2*nu));

      D_r00 = constt*(1-nu);
      D_r01 = constt*nu;
      D_r33 = constt*(1-2*nu)/2;
      
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

      D_rs00 = -(2*nu*(nu - 2))/pow((2*nu*nu + nu - 1),2);
      D_rs01 = (2*nu*nu + 1)/pow((2*nu*nu + nu - 1),2);
      D_rs33 = -1/(2*pow((nu + 1),2));

      
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

#endif //__K_RSPOISSONFEMETHOD_HPP__
