/* Copyright (C) 2015, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#ifndef __ELASTICFEMETHOD_MATRIX_HPP__
#define __ELASTICFEMETHOD_MATRIX_HPP__

namespace ObosleteUsersModules {
  
  struct ElasticFEMethod_Matrix: public ElasticFEMethod {
    
    double wt_Gauss;
    ElasticFEMethod_Matrix(FieldInterface& _mField,Mat _Aij,Vec _X,Vec _F,double _lambda,double _mu, double _wt_Gauss, string _field_name = "DISPLACEMENT"):
    ElasticFEMethod(_mField,_Aij,_X,_F,_lambda,_mu,_field_name), wt_Gauss(_wt_Gauss) {}
    
    virtual PetscErrorCode calculateD(double _lambda,double _mu) {
      PetscFunctionBegin;
      D = _lambda*D_lambda + _mu*D_mu;
//      cout<<" D before Degradation "<<D<<endl;
      wt_Gauss=1.0;
      D=wt_Gauss*D;
//      cout<<" wt_Gauss "<<wt_Gauss<<endl;
//      cout<<" D After Degradation "<<D<<endl;
      PetscFunctionReturn(0);
    }
    
  };
  
  
}

#endif //__ELASTICFEMETHOD_MATRIX_HPP__
