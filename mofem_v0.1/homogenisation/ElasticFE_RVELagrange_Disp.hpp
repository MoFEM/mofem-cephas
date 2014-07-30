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

#ifndef __ELASTICFE_RVELAGRANGE_DISP_HPP__
#define __ELASTICFE_RVELAGRANGE_DISP_HPP__

#include "FEMethod_UpLevelStudent.hpp"
#include <moab/ParallelComm.hpp>

namespace MoFEM {
  
  struct ElasticFE_RVELagrange_Disp: public FEMethod_UpLevelStudent {
    
    FieldInterface& mField;
    Mat Aij;
    Vec F;
    ublas::vector<FieldData> applied_strain;
    const string field_main;
    const string field_lagrange;
    int rank_field;
    
    ElasticFE_RVELagrange_Disp(FieldInterface& _mField,Mat &_Aij,Vec &_D,Vec& _F,ublas::vector<FieldData> _applied_strain,
                               const string& _field_main, const string& _field_lagrange, int _rank_field);
    
    
    ErrorCode rval;
    ParallelComm* pcomm;
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;
    
    vector<double> g_NTRI;
    ublas::matrix<double> g_NTRI_mat;
    const double* G_W_TRI;
    
    PetscErrorCode preProcess();
    
    int row_mat,col_mat,g_TRI_dim;
    vector<vector<DofIdx> > RowGlob;
    vector<vector<DofIdx> > ColGlob;
    vector<vector<ublas::matrix<double> > > rowNMatrices;

    PetscErrorCode postProcess();
    
    //Find indices and calculate shape function  and arrange in required form
    double coords_face[9];
    double area;
    virtual PetscErrorCode GetN_and_Indices();
    
    //Calculate H matrix
    ublas::vector<ublas::matrix<FieldData> > H_mat;
    virtual PetscErrorCode Get_H_mat();
    
    //Calculate and assemble NT x N matrix
    ublas::matrix<ublas::matrix<FieldData> > NTN;
    virtual PetscErrorCode Stiffness();
    virtual PetscErrorCode Lhs();
    
    //Calculate the right hand side vector, i.e. f=D_max * applied_strain and assemble it into the global force vector F
    ublas::matrix<FieldData> X_mat, nodes_coord, gauss_coord;
    ublas::vector<ublas::matrix<FieldData> > D_mat;
    ublas::vector<FieldData> f; //f.resize(9);
    virtual PetscErrorCode Rhs();
    
    //Loop over all the elements
    PetscErrorCode operator()();
    
  };
  
}

#endif
