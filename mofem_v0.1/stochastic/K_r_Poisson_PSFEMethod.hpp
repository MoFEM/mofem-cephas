/* Copyright (C) 2014,
 *   Zahur Ullah (Zahur.Ullah@glasgow.ac.uk)
 *   Xiao-Yi Zhou (xiaoyi.zhou AT newcastle.ac.uk)             
 * --------------------------------------------------------------
 * This routine calculates the right-hand side of the 1st-order finite-element 
 * equilibrium equation, Rhs = [K_r][U] for implementing 
 * perturbation-based stochastic finite element (PSFE).
 *
 * HISTORY
 *
 * 2014.08.29 (first version)
 *
 * REFERENCES
 * 1. Kleiber M. and Hien T. D. (1992) The stochastic finite element method - 
 *      Basic perturbation technique and computer implementation. John Wiley & 
 *      Sons.
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

#ifndef __K_R_POISSON_PSFEMETHOD_HPP__
#define __K_R_POISSON_PSFEMETHOD_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
extern "C" {
#include <gm_rule.h>
}

namespace MoFEM {
  
  struct K_r_Poisson_PSFEMethod: public K_r_PSFEMethod {

    K_r_Poisson_PSFEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young,double _pois):
    K_r_PSFEMethod( _mField,_Aij,_X,_F,_young,_pois) {
      
    }
    
    
    virtual PetscErrorCode calculateD_r(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear();
      double D00,D01,D33;
      
      D00=-(2*young*nu*(nu - 2))/pow((2*nu*nu + nu - 1),2);
      D01=(young*(2*nu*nu + 1))/pow((2*nu*nu + nu - 1),2);
      D33=-young/(2*pow((nu + 1),2));
      
      D(0,0)=D00;  D(0,1)=D01;  D(0,2)=D01;
      D(1,0)=D01;  D(1,1)=D00;  D(1,2)=D01;
      D(2,0)=D01;  D(2,1)=D01;   D(2,2)=D00;
      D(3,3)=D33;
      D(4,4)=D33;
      D(5,5)=D33;
      //      cout<<"D = "<<D;
      PetscFunctionReturn(0);
    }

        
  };
  
  
}

#endif //__K_R_POISSON_PSFEMETHOD_HPP__
