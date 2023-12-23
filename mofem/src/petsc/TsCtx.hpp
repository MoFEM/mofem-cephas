/** \file TsCtx.hpp
 * \brief Context for PETSc Time Stepping
 */

#ifndef __TSCTX_HPP__
#define __TSCTX_HPP__

#include <petsc/private/tsimpl.h>

#define TSADAPTMOFEM "TSMoFEMAdapt"

namespace MoFEM {

/** \brief Interface for Time Stepping (TS) solver
 * \ingroup mofem_petsc_solvers
 */
struct TsCtx {

  MoFEM::Interface &mField;
  moab::Interface &moab;

  std::string problemName;
  MoFEMTypes bH; ///< If set to MF_EXIST check if element exist

  typedef MoFEM::PairNameFEMethodPtr PairNameFEMethodPtr;
  typedef MoFEM::FEMethodsSequence FEMethodsSequence;
  typedef MoFEM::BasicMethodsSequence BasicMethodsSequence;

  FEMethodsSequence loopsIJacobian;
  FEMethodsSequence loopsIFunction;
  FEMethodsSequence loopsMonitor;
  FEMethodsSequence loopsRHSJacobian;
  FEMethodsSequence loopsRHSFunction;
  BasicMethodsSequence preProcessIJacobian;
  BasicMethodsSequence postProcessIJacobian;
  BasicMethodsSequence preProcessIFunction;
  BasicMethodsSequence postProcessIFunction;
  BasicMethodsSequence preProcessMonitor;
  BasicMethodsSequence postProcessMonitor;
  BasicMethodsSequence preProcessRHSJacobian;
  BasicMethodsSequence preProcessRHSFunction;
  BasicMethodsSequence postProcessRHSJacobian;
  BasicMethodsSequence postProcessRHSFunction;

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

  virtual ~TsCtx() = default;

  /**
   * @brief Get the loops to do IFunction object
   *
   * It is sequence of finite elements used to evaluate the right hand side of
   * implicit time solver.
   *
   * @return FEMethodsSequence&
   */
  FEMethodsSequence &getLoopsIFunction() { return loopsIFunction; }

  /**
   * @brief Get the loops to do RHSFunction object
   *
   * It is sequence of finite elements used to evaluate the right hand side of
   * implicit time solver.
   *
   * @return FEMethodsSequence&
   */
  FEMethodsSequence &getLoopsRHSFunction() { return loopsRHSFunction; }

  /**
   * @brief Get the loops to do IJacobian object
   *
   * It is sequence of finite elements used to evalite the left hand sie of
   * implicit time solver.
   *
   * @return FEMethodsSequence&
   */
  FEMethodsSequence &getLoopsIJacobian() { return loopsIJacobian; }

  /**
   * @brief Get the loops to do RHSJacobian object
   *
   * It is sequence of finite elements used to evalite the left hand sie of
   * implicit time solver.
   *
   * @return FEMethodsSequence&
   */
  FEMethodsSequence &getLoopsRHSJacobian() { return loopsRHSJacobian; }

  /**
   * @brief Get the loops to do Monitor object
   *
   * It is sequence used to monitor solution of time solver.
   *
   * @return FEMethodsSequence&
   */
  FEMethodsSequence &getLoopsMonitor() { return loopsMonitor; }

  /**
   * @brief Get the preProcess to do IFunction object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &getPreProcessIFunction() { return preProcessIFunction; }

  /**
   * @brief Get the postProcess to do IFunction object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &getPostProcessIFunction() {
    return postProcessIFunction;
  }

  /**
   * @brief Get the preProcess to do IJacobian object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &getPreProcessIJacobian() { return preProcessIJacobian; }

  /**
   * @brief Get the postProcess to do IJacobian object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &getPostProcessIJacobian() {
    return postProcessIJacobian;
  }

  /**
   * @brief Get the preProcess to do Monitor object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &getPreProcessMonitor() { return preProcessMonitor; }

  /**
   * @brief Get the postProcess to do Monitor object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &getPostProcessMonitor() { return postProcessMonitor; }

  /**
   * @brief Get the preProcess to do RHSJacobian object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &getPreProcessRHSJacobian() {
    return preProcessRHSJacobian;
  }

  /**
   * @brief Get the postProcess to do RHSJacobian object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &getPostProcessRHSJacobian() {
    return postProcessRHSJacobian;
  }

  /**
   * @brief Get the preProcess to do RHSFunction object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &getPreProcessRHSFunction() {
    return preProcessRHSFunction;
  }

  /**
   * @brief Get the postProcess to do RHSFunction object
   *
   * @return BasicMethodsSequence&
   */
  BasicMethodsSequence &getPostProcessRHSFunction() {
    return postProcessRHSFunction;
  }

  /**
   * @brief Clear loops
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode clearLoops();

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
 * @brief Set function evaluating jacobian in TS solver
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
 * @brief Calculation Jacobian for second order PDE in time
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetI2Jacobian.html>See
 * PETSc for details</a>
 *
 * @param ts
 * @param t time at step/stage being solved
 * @param u state vectora
 * @param u_t time derivative of state vector
 * @param u_tt second time derivative of state vector
 * @param a shift for u_t
 * @param aa shift for u_tt
 * @param A Jacobian of G(U) = F(t,U,W+v*U,W'+a*U), equivalent to dF/dU +
 * v*dF/dU_t + a*dF/dU_tt
 * @param B preconditioning matrix for J, may be same as J
 * @param ctx TsCtx context for matrix evaluation routine
 * @return PetscErrorCode
 */
PetscErrorCode TsSetI2Jacobian(TS ts, PetscReal t, Vec u, Vec u_t, Vec u_tt,
                               PetscReal a, PetscReal aa, Mat A, Mat B,
                               void *ctx);

/**
 * @brief Calculation the right hand side for second order PDE in time
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetI2Function.html>PETSc
 * for details</a>
 *
 * @param ts
 * @param t
 * @param u
 * @param u_t
 * @param u_tt
 * @param F
 * @param ctx
 * @return PetscErrorCode
 */
PetscErrorCode TsSetI2Function(TS ts, PetscReal t, Vec u, Vec u_t, Vec u_tt,
                               Vec F, void *ctx);

/** \brief Custom TSAdaptivity in MoFEM
 * 
 * \code
 *  CHKERR TSAdaptRegister(TSADAPTMOFEM, TSAdaptCreateMoFEM);
    TSAdapt adapt;
    CHKERR TSGetAdapt(solver, &adapt);
    \endcode
 *
 */
struct TSAdaptMoFEM {

  TSAdaptMoFEM() : alpha(0.75), gamma(0.5), desiredIt(6) {
    CHKERR PetscOptionsGetScalar("", "-ts_mofem_adapt_alpha", &alpha,
                                 PETSC_NULL);
    CHKERR PetscOptionsGetScalar("", "-ts_mofem_adapt_gamma", &gamma,
                                 PETSC_NULL);
    CHKERR PetscOptionsGetInt("", "-ts_mofem_adapt_desired_it", &desiredIt,
                              PETSC_NULL);
    MOFEM_LOG_TAG("WORLD", "TSMoFEMAdapt");
    MOFEM_LOG("WORLD", Sev::inform)
        << "TS adaptivity: alpha = " << alpha << ", gamma = " << gamma
        << ", desiredIt = " << desiredIt;
  }

  double alpha;
  double gamma;
  int desiredIt;
};

static PetscErrorCode TSAdaptChooseMoFEM(TSAdapt adapt, TS ts, PetscReal h,
                                         PetscInt *next_sc, PetscReal *next_h,
                                         PetscBool *accept, PetscReal *wlte,
                                         PetscReal *wltea, PetscReal *wlter);
static PetscErrorCode TSAdaptResetMoFEM(TSAdapt adapt);

static PetscErrorCode TSAdaptDestroyMoFEM(TSAdapt adapt);

/**
 * @brief Craete MOFEM adapt 
 *
 * \code
 *  CHKERR TSAdaptRegister(TSADAPTMOFEM, TSAdaptCreateMoFEM);
 *   TSAdapt adapt;
 *  CHKERR TSGetAdapt(solver, &adapt);
 * \endcode
 *
 * @param adapt
 * @return PetscErrorCode
 */
PetscErrorCode TSAdaptCreateMoFEM(TSAdapt adapt);

} // namespace MoFEM

#endif // __TSCTX_HPP__
