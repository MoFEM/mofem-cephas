/** \file SnesCtx.hpp
 * \brief Context for PETSc SNES, i.e. nonlinear solver
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

#ifndef __SNESCTX_HPP__
#define __SNESCTX_HPP__

namespace MoFEM {

/** \brief Interface for nonlinear (SNES) solver
  * \ingroup petsc_context_struture
  */
struct SnesCtx {

  ErrorCode rval;
  PetscErrorCode ierr;

  MoFEM::Interface &mField;
  moab::Interface &moab;

  std::string problemName;
  bool zeroPreCondMatrixB;

  typedef std::pair<std::string,FEMethod*> loop_pair_type;
  typedef std::vector<loop_pair_type > loops_to_do_type;
  loops_to_do_type loops_to_do_Mat;
  loops_to_do_type loops_to_do_Rhs;

  typedef std::vector<BasicMethod*> basic_method_to_do;
  basic_method_to_do preProcess_Mat;
  basic_method_to_do postProcess_Mat;
  basic_method_to_do preProcess_Rhs;
  basic_method_to_do postProcess_Rhs;

  PetscLogEvent USER_EVENT_SnesRhs;
  PetscLogEvent USER_EVENT_SnesMat;

  SnesCtx(Interface &m_field,const std::string &problem_name):
    mField(m_field),
    moab(m_field.get_moab()),
    problemName(problem_name),
    zeroPreCondMatrixB(true) {
    PetscLogEventRegister("LoopSNESRhs",0,&USER_EVENT_SnesRhs);
    PetscLogEventRegister("LoopSNESMat",0,&USER_EVENT_SnesMat);
  }
  virtual ~SnesCtx() {}

  const MoFEM::Interface& getm_field() const { return mField; }
  const moab::Interface& get_moab() const { return moab; }

  loops_to_do_type& get_loops_to_do_Mat() { return loops_to_do_Mat; }
  loops_to_do_type& get_loops_to_do_Rhs() { return loops_to_do_Rhs; }
  basic_method_to_do& get_preProcess_to_do_Rhs() { return preProcess_Rhs; }
  basic_method_to_do& get_postProcess_to_do_Rhs() { return postProcess_Rhs; }
  basic_method_to_do& get_preProcess_to_do_Mat() { return preProcess_Mat; }
  basic_method_to_do& get_postProcess_to_do_Mat() { return postProcess_Mat; }

  friend PetscErrorCode SnesRhs(SNES snes,Vec x,Vec f,void *ctx);
  friend PetscErrorCode SnesMat(SNES snes,Vec x,Mat A,Mat B,void *ctx);

};

/**
 * \brief This is MoFEM implementation for the right side evaluation in SNES solver
 *
 * For more information pleas look to PETSc manual, i.e. SNESSetFunction
 * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESSetFunction.html>
 *
 * @param  snes SNES solver
 * @param  x    Solution vector at current iteration
 * @param  f    The right hand side vector
 * @param  ctx  Pointer to context thata, i.e. SnesCtx
 * @return      Error code
 */
PetscErrorCode SnesRhs(SNES snes,Vec x,Vec f,void *ctx);

PetscErrorCode SnesMat(SNES snes,Vec x,Mat A,Mat B,void *ctx);

}

#endif // __SNESCTX_HPP__
