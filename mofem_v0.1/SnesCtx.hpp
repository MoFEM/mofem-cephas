/** \file SnesCtx.hpp 
 * \brief Context for PETSc SNES, i.e. nonlinear spolver
 */

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

namespace MoFEM {

struct SnesCtx {

  ErrorCode rval;
  PetscErrorCode ierr;

  FieldInterface &mField;
  Interface &moab;

  string problem_name;

  typedef pair<string,FieldInterface::FEMethod*> loop_pair_type;
  typedef vector<loop_pair_type > loops_to_do_type;
  loops_to_do_type loops_to_do_Mat;
  loops_to_do_type loops_to_do_Rhs;

  typedef vector<FieldInterface::BasicMethod*> basic_method_to_do;
  basic_method_to_do preProcess_Mat;
  basic_method_to_do postProcess_Mat;
  basic_method_to_do preProcess_Rhs;
  basic_method_to_do postProcess_Rhs;

  PetscLogEvent USER_EVENT_SnesRhs;
  PetscLogEvent USER_EVENT_SnesMat;

  SnesCtx(FieldInterface &_mField,const string &_problem_name): 
    mField(_mField),moab(_mField.get_moab()),problem_name(_problem_name) {
    PetscLogEventRegister("LoopSnesRhs",0,&USER_EVENT_SnesRhs);
    PetscLogEventRegister("LoopSnesMat",0,&USER_EVENT_SnesMat);
  }

  const FieldInterface& get_mField() const { return mField; }
  const Interface& get_moab() const { return moab; }
  loops_to_do_type& get_loops_to_do_Mat() { return loops_to_do_Mat; }
  loops_to_do_type& get_loops_to_do_Rhs() { return loops_to_do_Rhs; }
  basic_method_to_do& get_preProcess_to_do_Rhs() { return preProcess_Rhs; }
  basic_method_to_do& get_postProcess_to_do_Rhs() { return postProcess_Rhs; }
  basic_method_to_do& get_preProcess_to_do_Mat() { return preProcess_Mat; }
  basic_method_to_do& get_postProcess_to_do_Mat() { return postProcess_Mat; }

  friend PetscErrorCode SnesRhs(SNES snes,Vec x,Vec f,void *ctx);
  friend PetscErrorCode SnesFunc(SNES snes,Vec x,Vec f,SnesCtx *);
  friend PetscErrorCode SnesMat(SNES snes,Vec x,Mat *A,Mat *B,MatStructure *flag,void *ctx);

};

PetscErrorCode SnesRhs(SNES snes,Vec x,Vec f,void *ctx);
PetscErrorCode SnesFunc(SNES snes,Vec x,Vec f,SnesCtx *);
PetscErrorCode SnesMat(SNES snes,Vec x,Mat *A,Mat *B,MatStructure *flag,void *ctx);

}

#endif // __MOABSNES_HPP__
