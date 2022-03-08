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
 * \ingroup mofem_petsc_solvers
 */
struct KspCtx {

  MoFEM::Interface &mField;
  moab::Interface &moab;

  std::string problemName; ///< Problem name
  MoFEMTypes bH;           ///< If set to MF_EXIST check if element exist

  typedef MoFEM::PairNameFEMethodPtr PairNameFEMethodPtr;
  typedef MoFEM::FEMethodsSequence FEMethodsSequence;
  typedef MoFEM::BasicMethodsSequence BasicMethodsSequence;

  FEMethodsSequence loops_to_do_Mat; ///< Sequence of finite elements instances
                                     ///< assembling tangent matrix
  FEMethodsSequence loops_to_do_Rhs;   ///< Sequence of finite elements
                                       ///< instances assembling residual vector
  BasicMethodsSequence preProcess_Mat; ///< Sequence of methods run before
                                       ///< tangent matrix is assembled
  BasicMethodsSequence postProcess_Mat; ///< Sequence of methods run after
                                        ///< tangent matrix is assembled
  BasicMethodsSequence preProcess_Rhs;  ///< Sequence of methods run before
                                        ///< residual is assembled
  BasicMethodsSequence postProcess_Rhs; ///< Sequence of methods run after
                                        ///< residual is assembled

  KspCtx(MoFEM::Interface &m_field, const std::string &_problem_name)
      : mField(m_field), moab(m_field.get_moab()), problemName(_problem_name),
        bH(MF_EXIST) {
    PetscLogEventRegister("LoopKSPRhs", 0, &MOFEM_EVENT_KspRhs);
    PetscLogEventRegister("LoopKSPMat", 0, &MOFEM_EVENT_KspMat);
  }
  virtual ~KspCtx() = default;

  /**
   * @return return reference to vector with FEMethod to calculate matrix
   */
  FEMethodsSequence &get_loops_to_do_Mat() { return loops_to_do_Mat; }

  /**
   * @return return vector to vector with FEMethod to vector
   */
  FEMethodsSequence &get_loops_to_do_Rhs() { return loops_to_do_Rhs; }

  /**
   * The sequence of BasicMethod is executed before residual is calculated. It
   * can be used to setup data structures, e.g. zero global variable which is
   * integrated in domain, e.g. for calculation of strain energy.
   *
   * @return reference to BasicMethod for preprocessing
   */
  BasicMethodsSequence &get_preProcess_to_do_Rhs() { return preProcess_Rhs; }

  /**
   * The sequence of BasicMethod is executed after residual is calculated. It
   * can be used to setup data structures, e.g. aggregate data from processors
   * or to apply essential boundary conditions.
   *
   * @return reference to BasicMethod for postprocessing
   */
  BasicMethodsSequence &get_postProcess_to_do_Rhs() { return postProcess_Rhs; }

  /**
   * @return reference to BasicMethod for preprocessing
   */
  BasicMethodsSequence &get_preProcess_to_do_Mat() { return preProcess_Mat; }

  /**
   * The sequence of BasicMethod is executed after tangent matrix is calculated.
   * It can be used to setup data structures, e.g. aggregate data from
   * processors or to apply essential boundary conditions.
   *
   * @return reference to BasicMethod for postprocessing
   */
  BasicMethodsSequence &get_postProcess_to_do_Mat() { return postProcess_Mat; }

  friend PetscErrorCode KspRhs(KSP ksp, Vec f, void *ctx);
  friend PetscErrorCode KspMat(KSP ksp, Mat A, Mat B, void *ctx);

  /**
   * @brief Clear loops 
   * 
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode clearLoops();

private:
  PetscLogEvent MOFEM_EVENT_KspRhs;
  PetscLogEvent MOFEM_EVENT_KspMat;

  boost::movelib::unique_ptr<bool> vecAssembleSwitch;
  boost::movelib::unique_ptr<bool> matAssembleSwitch;
};

/**
 * \brief Run over elements in the lists
 * @param  ksp KSP solver
 * @param  f   the right hand side vector
 * @param  ctx data context, i.e. KspCtx
 * @return     error code
 */
PetscErrorCode KspRhs(KSP ksp, Vec f, void *ctx);

/**
 * \brief Run over elenents in the list
 * @param  ksp KSP solver
 * @param  A   matrix
 * @param  B   Preconditioner matrix
 * @param  ctx data context, i.e. KspCtx
 * @return     error code
 */
PetscErrorCode KspMat(KSP ksp, Mat A, Mat B, void *ctx);

} // namespace MoFEM

#endif // __KSPCTX_HPP__
