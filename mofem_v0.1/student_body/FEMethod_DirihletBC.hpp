/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
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

#ifndef __MOABFEMETHOD_DIRICHLETBC_HPP__
#define __MOABFEMETHOD_DIRICHLETBC_HPP__

#include "CoreDataStructures.hpp"
#include "FieldInterface.hpp"
#include "FEMethod_LowLevelStudent.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

/** 
 * \brief The student user interface for Dirichlet boundary conditions
 * 
*/
struct BaseDirichletBC {

  BaseDirichletBC();

  virtual PetscErrorCode SetDirichletBC_to_ElementIndiciesRow(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirichletBC);
  virtual PetscErrorCode SetDirichletBC_to_ElementIndiciesCol(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirichletBC);
  virtual PetscErrorCode SetDirichletBC_to_ElementIndicies(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirichletBC);
  virtual PetscErrorCode SetDirichletBC_to_ElementIndiciesFace(
    vector<DofIdx>& DirichletBC,vector<DofIdx>& FaceNodeGlobalDofs,vector<vector<DofIdx> > &FaceEdgeGlobalDofs,vector<DofIdx> &FaceGlobalDofs);
  virtual PetscErrorCode SetDirichletBC_to_FieldData(FieldInterface::FEMethod *fe_method_ptr,Vec D);
  virtual PetscErrorCode SetDirichletBC_to_MatrixDiagonal(FieldInterface::FEMethod *fe_method_ptr,Mat Aij);
  virtual PetscErrorCode SetDirichletBC_to_RHS(FieldInterface::FEMethod *fe_method_ptr,Vec F);

};

struct CubitDisplacementDirichletBC: public BaseDirichletBC {
  FieldInterface& mField;
  string problem_name;  
  string field_name;

  CubitDisplacementDirichletBC(FieldInterface& _mField,const string _problem_name,const string _field_name); 

  PetscErrorCode ierr;
  ErrorCode rval;

  map<int,Range> bc_map[3];
  map<int,double> bc_map_val[3];

  PetscErrorCode Init();

  PetscErrorCode SetDirichletBC_to_ElementIndiciesRow(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirichletBC);
  PetscErrorCode SetDirichletBC_to_ElementIndiciesCol(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirichletBC);
  PetscErrorCode SetDirichletBC_to_ElementIndicies(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirichletBC);
  PetscErrorCode SetDirichletBC_to_ElementIndiciesFace(
    vector<DofIdx>& DirichletBC,vector<DofIdx> &FaceNodeIndices, vector<vector<DofIdx> > &FaceEdgeIndices, vector<DofIdx> &FaceIndices);
  PetscErrorCode SetDirichletBC_to_MatrixDiagonal(FieldInterface::FEMethod *fe_method_ptr,Mat Aij);
  PetscErrorCode SetDirichletBC_to_RHS(FieldInterface::FEMethod *fe_method_ptr,Vec F);
  PetscErrorCode SetDirichletBC_to_FieldData(FieldInterface::FEMethod *fe_method_ptr,Vec D);

};
    
struct CubitDisplacementDirihletBC_MatZeroRowsColumns: public CubitDisplacementDirihletBC {

  
    
};

}

#endif //__MOABFEMETHOD_DIRICHLETBC_HPP__
