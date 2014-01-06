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

#ifndef __MOABSURFACECONSTRAINS_HPP__
#define __MOABSURFACECONSTRAINS_HPP__

#include "FieldInterface.hpp"
#include "DirihletBC.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

struct C_SURFACE_FEMethod:public FieldInterface::FEMethod {
  ErrorCode rval;
  PetscErrorCode ierr;
  Interface& moab;
  BaseDirihletBC *dirihlet_bc_method_ptr;
  Mat C;
  string lambda_field_name;
  bool updated;

  vector<double> diffNTRI;
  vector<double> g_NTRI3;
  const double *G_TRI_W;

  C_SURFACE_FEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,Mat _C,string _lambda_field_name,int _verbose = 0);
  C_SURFACE_FEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,Mat _C,int _verbose = 0);
  void run_in_constructor();

  PetscErrorCode preProcess();

  ublas::bounded_matrix<double,3,9 > C_MAT_ELEM,iC_MAT_ELEM;
  ublas::bounded_matrix<double,9,3 > CT_MAT_ELEM;
  ublas::bounded_matrix<double,3,9 > dC_MAT_ELEM;
  ublas::bounded_matrix<double,9,9 > dCT_MAT_ELEM;
  ublas::vector<double,ublas::bounded_array<double,3> > ig_VEC_ELEM;
  ublas::vector<double,ublas::bounded_array<double,9> > if_VEC_ELEM;

  ublas::vector<DofIdx,ublas::bounded_array<DofIdx,3> > lambda_global_row_indices,lambda_global_col_indices;
  vector<DofIdx> DirihletBC;
  vector<vector<DofIdx> > dof_global_row_indices,dof_global_col_indices;
  ublas::vector<double,ublas::bounded_array<double,3> > ent_lambda_data;
  ublas::vector<double,ublas::bounded_array<double,9> > ent_dofs_data,ent_idofs_data;
  ublas::vector<double,ublas::bounded_array<double,9> > coords;

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = operator()(false,false); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  };

  PetscErrorCode cOnstrain(double *dofs_X,double *dofs_iX,double *C,double *iC);
  virtual PetscErrorCode iNtegrate(bool transpose,bool nonlinear);
  virtual PetscErrorCode aSsemble(bool transpose,bool nonlinear);
  virtual PetscErrorCode operator()(bool transpose,bool nonlinear);
  PetscErrorCode postProcess();

};

struct g_SURFACE_FEMethod: public C_SURFACE_FEMethod {
  Vec g;
  g_SURFACE_FEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,Vec _g,string _lambda_field_name,int _verbose = 0); 
  g_SURFACE_FEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,Vec _g,int _verbose = 0); 

  ublas::vector<double,ublas::bounded_array<double,3> > g_VEC_ELEM;
  ublas::vector<double,ublas::bounded_array<double,9> > f_VEC_ELEM;
  PetscErrorCode iNtegrate(bool transpose,bool nonlinear);
  PetscErrorCode aSsemble(bool transpose,bool nonlinear);
};

struct C_FEMethod_ForSnes: public FieldInterface::FEMethod {

  PetscErrorCode setElemData(FieldInterface::FEMethod &e) {
    PetscFunctionBegin;
    //copy all elem that to other FEMethod class
    e.problem_ptr = problem_ptr;
    e.moabfields = moabfields;
    e.ents_moabfield = ents_moabfield;
    e.dofs_moabfield = dofs_moabfield;
    e.finite_elements = finite_elements;
    e.finite_elements_moabents = finite_elements_moabents;
    e.fem_adjacencies = fem_adjacencies;
    e.fe_name = fe_name;
    e.fe_ptr = fe_ptr;
    e.data_multiIndex = data_multiIndex;
    e.row_multiIndex = row_multiIndex;
    e.col_multiIndex = col_multiIndex;
    PetscFunctionReturn(0);
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

};

struct C_SURFACE_FEMethod_ForSnes: public C_FEMethod_ForSnes {

  FieldInterface& mField;
  C_SURFACE_FEMethod MatMethod;
  g_SURFACE_FEMethod FMethod;
  bool nonlinear;

  C_SURFACE_FEMethod_ForSnes(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,int _verbose = 0):
   mField(_mField),
    MatMethod(_mField.get_moab(),_dirihlet_bc_method_ptr,PETSC_NULL),
    FMethod(_mField.get_moab(),_dirihlet_bc_method_ptr,snes_f),nonlinear(false) {}

  C_SURFACE_FEMethod_ForSnes(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,string _lambda_field_name,int _verbose = 0):
   mField(_mField),
    MatMethod(_mField.get_moab(),_dirihlet_bc_method_ptr,PETSC_NULL,_lambda_field_name),
    FMethod(_mField.get_moab(),_dirihlet_bc_method_ptr,snes_f,_lambda_field_name),nonlinear(false) {}

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	FMethod.g = snes_f;
	ierr = setElemData(FMethod); CHKERRQ(ierr);
	ierr = FMethod.operator()(true,nonlinear); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	MatMethod.C = *snes_B;
	ierr = setElemData(MatMethod); CHKERRQ(ierr);
	ierr = MatMethod.operator()(true,nonlinear); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatMethod.postProcess(); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }


};

}

#endif //__MOABSURFACECONSTRAINS_HPP__
