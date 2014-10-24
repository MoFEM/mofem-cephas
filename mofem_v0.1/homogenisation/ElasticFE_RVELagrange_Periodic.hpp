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

#ifndef __ELASTICFE_RVELAGRANGE_PERIODIC_HPP__
#define __ELASTICFE_RVELAGRANGE_PERIODIC_HPP__

#include "ElasticFE_RVELagrange_Disp.hpp"

namespace MoFEM {
  
  struct ElasticFE_RVELagrange_Periodic: public ElasticFE_RVELagrange_Disp {
    
    
    ElasticFE_RVELagrange_Periodic(
                                   FieldInterface& _mField,Mat &_Aij,Vec &_D,Vec& _F, ublas::vector<FieldData>_applied_strain, const string& _field_main, const string& _field_lagrange, int _rank_field):
    ElasticFE_RVELagrange_Disp(_mField,_Aij, _D, _F, _applied_strain, _field_main, _field_lagrange, _rank_field){};
    
    vector<vector<vector<DofIdx> > > RowGlob;  //The outer vector is of size 2 (one for -ve triangles and one for +ve triangles)
    vector<vector<vector<DofIdx> > > ColGlob;
    EntityHandle prism_periodic;
    double coords_prism[18];
    const EntityHandle* conn_Prism;
    int  num_nodes1;
    double coords_face[2][9];
    double area;  //area is the same for two triangular faces of the prism so don't need area[2]
    
    virtual PetscErrorCode GetN_and_Indices();
    
    ublas::vector<ublas::vector<ublas::matrix<FieldData> > > H_mat;
    virtual PetscErrorCode Get_H_mat();
    
    //Calculate and assemble NT x N matrix
    ublas::vector<ublas::matrix<ublas::matrix<FieldData> > > NTN;  //The outer most vector is of size(2) for negative and positive triangles
    
    virtual PetscErrorCode Stiffness();
    virtual PetscErrorCode Lhs();

    //Calculate the right hand side vector, i.e. f=D_mat * applied_strain and assemble it into the global force vector F
    virtual PetscErrorCode Rhs();
    PetscErrorCode operator()();
    
    
  };
  
  
}

#endif //__ElasticFE_RVELagrange_Periodic__
