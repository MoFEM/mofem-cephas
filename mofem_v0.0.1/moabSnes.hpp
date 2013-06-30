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

#ifndef __MOABSNES_HPP__
#define __MOABSNES_HPP__

#include "moabField.hpp"
#include <petsc.h>
#include <petscmat.h>
#include <petscsnes.h>


namespace MoFEM {

struct moabSnesCtx {

  ErrorCode rval;
  PetscErrorCode ierr;

  moabField &mField;
  Interface &moab;
  ParallelComm* pcomm;

  string problem_name;
  typedef vector< pair<string,moabField::FEMethod*> > loops_to_do_type;
  loops_to_do_type loops_to_do;

  moabSnesCtx(moabField &_mField,const string &_problem_name): 
    mField(_mField),moab(_mField.get_moab()),
    problem_name(_problem_name) {
    pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);

  }

  const moabField& get_mField() const { return mField; }
  const Interface& get_moab() const { return moab; }

  friend PetscErrorCode SnesFunc(SNES snes,Vec x,Vec f,moabSnesCtx *);

};

PetscErrorCode SnesRhs(SNES snes,Vec x,Vec f,moabSnesCtx *ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  moabSnesCtx::loops_to_do_type::iterator lit = ctx->loops_to_do.begin();
  for(;lit!=ctx->loops_to_do.end();lit++) {
    ierr = lit->second->set_snes(snes); CHKERRQ(ierr);
    ierr = lit->second->set_x(x); CHKERRQ(ierr);
    ierr = lit->second->set_f(f); CHKERRQ(ierr);
    ierr = ctx->mField.loop_finite_elements(ctx->problem_name,lit->first,*(lit->second));  CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode SnesMat(SNES snes,Vec x,Mat *A,Mat *B,MatStructure *flag,moabSnesCtx *ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  moabSnesCtx::loops_to_do_type::iterator lit = ctx->loops_to_do.begin();
  for(;lit!=ctx->loops_to_do.end();lit++) {
    ierr = lit->second->set_snes(snes); CHKERRQ(ierr);
    ierr = lit->second->set_x(x); CHKERRQ(ierr);
    ierr = lit->second->set_A(A); CHKERRQ(ierr);
    ierr = lit->second->set_B(B); CHKERRQ(ierr);
    ierr = lit->second->set_flag(flag); CHKERRQ(ierr);
    ierr = ctx->mField.loop_finite_elements(ctx->problem_name,lit->first,*(lit->second));  CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

}

#endif // __MOABSNES_HPP__
