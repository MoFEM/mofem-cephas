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

#ifndef __MOABFEMETHOD_DIRIHLETBC_HPP__
#define __MOABFEMETHOD_DIRIHLETBC_HPP__

#include "Core_dataStructures.hpp"
#include "moabField.hpp"
#include "moabFEMethod_LowLevelStudent.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

/** 
 * \brief The student user interface for Dirihlet boundary conditions
 * 
*/
struct BaseDirihletBC {

  BaseDirihletBC() {};

  virtual PetscErrorCode SetDirihletBC_to_ElementIndicies(
    moabField::FEMethod *fe_method_ptr,string field_name,
    vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"sorry.. you need to tell me what to do");
    NOT_USED(fe_method_ptr);
    NOT_USED(field_name);
    NOT_USED(RowGlobDofs);
    NOT_USED(ColGlobDofs);
    NOT_USED(DirihletBC);
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode SetDirihletBC_to_ElementIndiciesFace(
    vector<DofIdx>& DirihletBC,vector<DofIdx>& FaceNodeGlobalDofs,vector<vector<DofIdx> > &FaceEdgeGlobalDofs,vector<DofIdx> &FaceGlobalDofs) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"sorry.. you need to tell me what to do");
    NOT_USED(FaceNodeGlobalDofs);
    NOT_USED(FaceEdgeGlobalDofs);
    NOT_USED(FaceGlobalDofs);
    NOT_USED(DirihletBC);
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode SetDirihletBC_to_MatrixDiagonal(
    moabField::FEMethod *fe_method_ptr,Mat Aij) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"sorry.. you need to tell me what to do");
    NOT_USED(fe_method_ptr);
    NOT_USED(Aij);
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode SetDirihletBC_to_RHS(moabField::FEMethod *fe_method_ptr) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"sorry.. you need to tell me what to do");
    NOT_USED(fe_method_ptr);
    PetscFunctionReturn(0);
  }

};

}

#endif //__MOABFEMETHOD_DIRIHLETBC_HPP__
