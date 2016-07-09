/** \file KspCtx.hpp
 * \brief Context for PETSc KSP, i.e. nonlinear solver
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

#ifndef __KSPCTX_HPP__
#define __KSPCTX_HPP__

namespace MoFEM {

/** \brief Interface for linear (KSP) solver
  * \ingroup petsc_context_struture
  */
struct KspCtx {

  ErrorCode rval;
  PetscErrorCode ierr;

  MoFEM::Interface &mField;
  moab::Interface &moab;

  std::string problemName;   ///< Problem name

  typedef std::pair<std::string,FEMethod*> loop_pair_type;
  typedef std::vector<loop_pair_type > loops_to_do_type;
  loops_to_do_type loops_to_do_Mat;
  loops_to_do_type loops_to_do_Rhs;

  typedef std::vector<BasicMethod*> basic_method_to_do;
  basic_method_to_do preProcess_Mat;
  basic_method_to_do postProcess_Mat;
  basic_method_to_do preProcess_Rhs;
  basic_method_to_do postProcess_Rhs;

  PetscLogEvent USER_EVENT_KspRhs;
  PetscLogEvent USER_EVENT_KspMat;

  KspCtx(MoFEM::Interface &m_field,const std::string &_problem_name):
    mField(m_field),
    moab(m_field.get_moab()),
    problemName(_problem_name) {
    PetscLogEventRegister("LoopKSPRhs",0,&USER_EVENT_KspRhs);
    PetscLogEventRegister("LoopKSPMat",0,&USER_EVENT_KspMat);
  }
  virtual ~KspCtx() {}

  const MoFEM::Interface& get_mField() const { return mField; }
  const moab::Interface& get_moab() const { return moab; }

  loops_to_do_type& get_loops_to_do_Mat() { return loops_to_do_Mat; }
  loops_to_do_type& get_loops_to_do_Rhs() { return loops_to_do_Rhs; }
  basic_method_to_do& get_preProcess_to_do_Rhs() { return preProcess_Rhs; }
  basic_method_to_do& get_postProcess_to_do_Rhs() { return postProcess_Rhs; }
  basic_method_to_do& get_preProcess_to_do_Mat() { return preProcess_Mat; }
  basic_method_to_do& get_postProcess_to_do_Mat() { return postProcess_Mat; }

  friend PetscErrorCode KspRhs(KSP ksp,Vec f,void *ctx);
  friend PetscErrorCode KspMat(KSP ksp,Mat A,Mat B,void *ctx);

};

PetscErrorCode KspRhs(KSP ksp,Vec f,void *ctx);
PetscErrorCode KspMat(KSP ksp,Mat A,Mat B,void *ctx);

}

/***************************************************************************//**
 * \defgroup petsc_context_struture Solver context structures
 * \brief Context structures used to exchange information between PETSc and MoFEM
 *
 * \ingroup mofem
 ******************************************************************************/




#endif // __KSPCTX_HPP__
