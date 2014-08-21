/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah@glasgow.ac.uk)
 * --------------------------------------------------------------
 * Implemnetation of the stochastic finite element componenet K_r
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

#ifndef __K_RYOUNGFEMETHOD_HPP__
#define __K_RYOUNGFEMETHOD_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "K_rPoissonFEMethod.hpp"

extern "C" {
#include <gm_rule.h>
}

namespace MoFEM {
  
  struct K_rYoungFEMethod: public K_rPoissonFEMethod {
    
    K_rYoungFEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young,double _pois):
    K_rPoissonFEMethod( _mField,_Aij,_X,_F,_young,_pois) {
      
    }
    
    PetscErrorCode calculateD(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear();
      
      double D00,D01,D33,constt;
      constt=1/((1+nu)*(1-2*nu));
      
      D00=constt*(1-nu);
      D01=constt*nu;
      D33=constt*(1-2*nu)/2;
      
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

#endif //__K_RYOUNGFEMETHOD_HPP__
