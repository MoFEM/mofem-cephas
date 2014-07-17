/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#ifndef __ELASTICFE_RVELAGRANGE_TRACTION_HPP__
#define __ELASTICFE_RVELAGRANGE_TRACTION_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFE_RVELagrange_Disp.hpp"

namespace MoFEM {
  
  struct ElasticFE_RVELagrange_Traction: public ElasticFE_RVELagrange_Disp {
    
    ElasticFE_RVELagrange_Traction(
                                   FieldInterface& _mField,Mat &_Aij,Vec &_D,Vec& _F,ublas::vector<FieldData> _applied_strain):
    ElasticFE_RVELagrange_Disp(_mField,_Aij, _D, _F, _applied_strain){};
    
    double coords_face[9];
    double area;
    virtual PetscErrorCode GetN_and_Indices();
    
    ublas::vector<ublas::matrix<FieldData> > H_mat;
    virtual PetscErrorCode Get_H_mat();
    
    //Calculate and assemble NT x N matrix
    //    ublas::matrix<ublas::matrix<FieldData> > NTN;
    virtual PetscErrorCode Stiffness();
    virtual PetscErrorCode Lhs();
    
    ublas::matrix<FieldData> X_mat, nodes_coord, gauss_coord;
    ublas::vector<ublas::matrix<FieldData> > D_mat;
    ublas::vector<FieldData> f; //f.resize(9);
    
    //Calculate the right hand side vector, i.e. f=D_mat * applied_strain and assemble it into the global force vector F
    virtual PetscErrorCode Rhs();
    
    PetscErrorCode operator()();
    
  };
}

#endif //__ElasticFE_RVELagrange_Traction_HPP__
