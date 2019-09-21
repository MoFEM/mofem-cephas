/** \file TsCtx.hpp
 * \brief Context for PETSc Time Stepping
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

#ifndef __TSCTX_HPP__
#define __TSCTX_HPP__

namespace MoFEM {

/** \brief Interface for Time Stepping (TS) solver
 * \ingroup petsc_context_struture
 */
struct TsCtx {

  MoFEM::Interface &mField;
  moab::Interface &moab;

  std::string problemName;
  MoFEMTypes bH; ///< If set to MF_EXIST check if element exist

  /// \deprecated use PairNameFEMethodPtr
  DEPRECATED typedef MoFEM::PairNameFEMethodPtr loop_pair_type;

  /// \deprecated use FEMethodsSequence
  DEPRECATED typedef MoFEM::FEMethodsSequence loops_to_do_type;

  /// \deprecated use BasicMethodsSequence
  DEPRECATED typedef MoFEM::BasicMethodsSequence basic_method_to_do;

  typedef MoFEM::PairNameFEMethodPtr PairNameFEMethodPtr;
  typedef MoFEM::FEMethodsSequence FEMethodsSequence;
  typedef MoFEM::BasicMethodsSequence BasicMethodsSequence;

  FEMethodsSequence loops_to_do_IJacobian;
  FEMethodsSequence loops_to_do_IFunction;
  FEMethodsSequence loops_to_do_Monitor;
  FEMethodsSequence loops_to_do_RHSJacobian;
  FEMethodsSequence loops_to_do_RHSFunction;
  BasicMethodsSequence preProcess_IJacobian;
  BasicMethodsSequence postProcess_IJacobian;
  BasicMethodsSequence preProcess_IFunction;
  BasicMethodsSequence postProcess_IFunction;
  BasicMethodsSequence preProcess_Monitor;
  BasicMethodsSequence postProcess_Monitor;
  BasicMethodsSequence preProcess_RHSJacobian;
  BasicMethodsSequence preProcess_RHSFunction;
  BasicMethodsSequence postProcess_RHSJacobian;
  BasicMethodsSequence postProcess_RHSFunction;

  bool zeroMatrix;

  TsCtx(MoFEM::Interface &m_field, const std::string &problem_name)
      : mField(m_field), moab(m_field.get_moab()), problemName(problem_name),
        bH(MF_EXIST), zeroMatrix(true) {
    PetscLogEventRegister("LoopTsIFunction", 0, &MOFEM_EVENT_TsCtxIFunction);
    PetscLogEventRegister("LoopTsIJacobian", 0, &MOFEM_EVENT_TsCtxIJacobian);
    PetscLogEventRegister("LoopTsRHSFunction", 0,
                          &MOFEM_EVENT_TsCtxRHSFunction);
    PetscLogEventRegister("LoopTsRHSJacobian", 0,
                          &MOFEM_EVENT_TsCtxRHSJacobian);
    PetscLogEventRegister("LoopTsMonitor", 0, &MOFEM_EVENT_TsCtxMonitor);
    PetscLogEventRegister("LoopTsI2Function", 0, &MOFEM_EVENT_TsCtxI2Function);
    PetscLogEventRegister("LoopTsI2Jacobian", 0, &MOFEM_EVENT_TsCtxI2Jacobian);
  }

  /**
   * @brief Get the loops to do IFunction object
   *
   * It is sequence of finite elements used to evaluate the right hand side of
   * implicit time solver.
   *
   * @return FEMethodsSequence&
   */
  FEMethodsSequence &get_loops_to_do_IFunction() {
    return loops_to_do_IFunction;
  }

  /**
   * @brief Get the loops to do RHSFunction object
   *
   * It is sequence of finite elements used to evaluate the right hand side of
   * implicit time solver.
   *
   * @return FEMethodsSequence&
   */
  FEMethodsSequence &get_loops_to_do_RHSFunction() {
    return loops_to_do_RHSFunction;
  }

  /**
   * @brief Get the loops to do IJacobian object
   *
   * It is sequence of finite elements used to evalite the left hand sie of
   * implimcit time solver.
   *
   * @return FEMethodsSequence&
   */
  FEMethodsSequence &get_loops_to_do_IJacobian() {
    return loops_to_do_IJacobian;
  }

  /**
   * @brief Get the loops to do RHSJacobian object
   *
   * It is sequence of finite elements used to evalite the left hand sie of
   * implimcit time solver.
   *
   * @return FEMethodsSequence&
   */
  FEMethodsSequence &get_loops_to_do_RHSJacobian() {
    return loops_to_do_RHSJacobian;
  }

  /**
   * @brief Get the loops to do Monitor object
   *
   * It is sequence used to monitor solution of time solver.
   *
   * @return FEMethodsSequence&
   */
  FEMethodsSequence &get_loops_to_do_Monitor() { return loops_to_do_Monitor; }

  /**
   * @brief Get the preProcess to do IFunction object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &get_preProcess_to_do_IFunction() {
    return preProcess_IFunction;
  }

  /**
   * @brief Get the postProcess to do IFunction object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &get_postProcess_to_do_IFunction() {
    return postProcess_IFunction;
  }

  /**
   * @brief Get the preProcess to do IJacobian object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &get_preProcess_to_do_IJacobian() {
    return preProcess_IJacobian;
  }

  /**
   * @brief Get the postProcess to do IJacobian object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &get_postProcess_to_do_IJacobian() {
    return postProcess_IJacobian;
  }

  /**
   * @brief Get the preProcess to do Monitor object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &get_preProcess_to_do_Monitor() {
    return preProcess_Monitor;
  }

  /**
   * @brief Get the postProcess to do Monitor object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &get_postProcess_to_do_Monitor() {
    return postProcess_Monitor;
  }

  /**
   * @brief Get the preProcess to do RHSJacobian object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &get_preProcess_to_do_RHSJacobian() {
    return preProcess_RHSJacobian;
  }

  /**
   * @brief Get the postProcess to do RHSJacobian object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &get_postProcess_to_do_RHSJacobian() {
    return postProcess_RHSJacobian;
  }

  /**
   * @brief Get the preProcess to do RHSFunction object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &get_preProcess_to_do_RHSFunction() {
    return preProcess_RHSFunction;
  }

  /**
   * @brief Get the postProcess to do RHSFunction object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &get_postProcess_to_do_RHSFunction() {
    return postProcess_RHSFunction;
  }

  friend PetscErrorCode TsSetIFunction(TS ts, PetscReal t, Vec u, Vec u_t,
                                       Vec F, void *ctx);
  friend PetscErrorCode TsSetIJacobian(TS ts, PetscReal t, Vec u, Vec U_t,
                                       PetscReal a, Mat A, Mat B, void *ctx);
  friend PetscErrorCode TsMonitorSet(TS ts, PetscInt step, PetscReal t, Vec u,
                                     void *ctx);
  friend PetscErrorCode TsSetRHSFunction(TS ts, PetscReal t, Vec u, Vec F,
                                         void *ctx);
  friend PetscErrorCode TsSetRHSJacobian(TS ts, PetscReal t, Vec u, Mat A,
                                         Mat B, void *ctx);

  friend PetscErrorCode TsSetI2Function(TS ts, PetscReal t, Vec U, Vec U_t,
                                        Vec U_tt, Vec F, void *ctx);

  friend PetscErrorCode TsSetI2Jacobian(TS ts, PetscReal t, Vec U, Vec U_t,
                                        Vec U_tt, PetscReal v, PetscReal a,
                                        Mat J, Mat P, void *ctx);

private:
  PetscLogEvent MOFEM_EVENT_TsCtxRHSFunction;
  PetscLogEvent MOFEM_EVENT_TsCtxRHSJacobian;
  PetscLogEvent MOFEM_EVENT_TsCtxIFunction;
  PetscLogEvent MOFEM_EVENT_TsCtxIJacobian;
  PetscLogEvent MOFEM_EVENT_TsCtxMonitor;
  PetscLogEvent MOFEM_EVENT_TsCtxI2Function;
  PetscLogEvent MOFEM_EVENT_TsCtxI2Jacobian;

  boost::movelib::unique_ptr<bool> vecAssembleSwitch;
  boost::movelib::unique_ptr<bool> matAssembleSwitch;
};

/**
 * @brief Set IFunction for TS solver
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetIFunction.html>See
 * petsc for details</a>
 *
 * @param ts
 * @param t
 * @param u
 * @param u_t
 * @param F
 * @param ctx
 * @return PetscErrorCode
 */
PetscErrorCode TsSetIFunction(TS ts, PetscReal t, Vec u, Vec u_t, Vec F,
                              void *ctx);

/**
 * @brief Set function evaluating jacobina in TS solver
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-3.1/docs/manualpages/TS/TSSetIJacobian.html>See
 * PETSc for details</a>
 *
 * @param ts
 * @param t
 * @param u
 * @param u_t
 * @param a
 * @param A
 * @param B
 * @param ctx
 * @return PetscErrorCode
 */
PetscErrorCode TsSetIJacobian(TS ts, PetscReal t, Vec u, Vec u_t, PetscReal a,
                              Mat A, Mat B, void *ctx);

/**
 * @brief Set monitor for TS solver
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSMonitorSet.html>See
 * PETSc for details</a>
 *
 * @param ts
 * @param step
 * @param t
 * @param u
 * @param ctx
 * @return PetscErrorCode
 */
PetscErrorCode TsMonitorSet(TS ts, PetscInt step, PetscReal t, Vec u,
                            void *ctx);

/// \deprecate Do not use, change to TsSetIFunction
DEPRECATED inline PetscErrorCode f_TSSetIFunction(TS ts, PetscReal t, Vec u,
                                                  Vec u_t, Vec F, void *ctx) {
  return TsSetIFunction(ts, t, u, u_t, F, ctx);
}

/// \deprecated Do not use, change to TsSetIJacobian
DEPRECATED inline PetscErrorCode f_TSSetIJacobian(TS ts, PetscReal t, Vec u,
                                                  Vec u_t, PetscReal a, Mat A,
                                                  Mat B, void *ctx) {
  return TsSetIJacobian(ts, t, u, u_t, a, A, B, ctx);
}

/// \deprecated Do not use, change to TsMonitorSet
DEPRECATED inline PetscErrorCode f_TSMonitorSet(TS ts, PetscInt step,
                                                PetscReal t, Vec u, void *ctx) {
  return TsMonitorSet(ts, step, t, u, ctx);
}

/**
 * @brief TS solver function
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-3.11/docs/manualpages/TS/TSSetRHSFunction.html#TSSetRHSFunction>See
 * PETSc for details</a>
 *
 * @param ts
 * @param t
 * @param u
 * @param F
 * @param ctx
 * @return PetscErrorCode
 */
PetscErrorCode TsSetRHSFunction(TS ts, PetscReal t, Vec u, Vec F, void *ctx);

/**
 * @brief TS solver function
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-3.11/docs/manualpages/TS/TSSetRHSJacobian.html#TSSetRHSJacobian>See
 * PETSc for details</a>
 *
 * @param ts
 * @param t
 * @param u
 * @param A
 * @param B
 * @param ctx
 * @return PetscErrorCode
 */
PetscErrorCode TsSetRHSJacobian(TS ts, PetscReal t, Vec u, Mat A, Mat B,
                                void *ctx);

/**
 * @brief Calculation Jaconian for second order PDE in time
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetI2Jacobian.html>PETSc
 * for details</a>
 *
 * @param ts
 * @param J
 * @param P
 * @param jac
 * @param ctx
 * @return PetscErrorCode
 */
PetscErrorCode TsSetI2Jacobian(TS ts, PetscReal t, Vec u, Vec u_t, Vec u_tt,
                               PetscReal v, PetscReal a, Mat J, Mat P,
                               void *ctx);

/**
 * @brief Calculation the right hand side for second order PDE in time
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetI2Function.html>PETSc
 * for details</a>
 *
 * @param ts
 * @param F
 * @param fun
 * @param ctx
 * @return PetscErrorCode
 */
PetscErrorCode TsSetI2Function(TS ts, PetscReal t, Vec u, Vec u_t, Vec u_tt,
                               Vec F, void *ctx);

} // namespace MoFEM

#endif // __TSCTX_HPP__
