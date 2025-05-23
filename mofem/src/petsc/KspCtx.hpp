/** \file KspCtx.hpp
 * \brief Context for PETSc KSP, i.e. nonlinear solver
 */

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
  FEMethodsSequence loops_to_do_Rhs; ///< Sequence of finite elements
                                     ///< instances assembling residual vector
  BasicMethodsSequence preProcess_Mat;  ///< Sequence of methods run before
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
  FEMethodsSequence &getSetOperators() { return loops_to_do_Mat; }

  /**
   * @return return vector to vector with FEMethod to vector
   */
  FEMethodsSequence &getComputeRhs() { return loops_to_do_Rhs; }

  /**
   * The sequence of BasicMethod is executed before residual is calculated. It
   * can be used to setup data structures, e.g. zero global variable which is
   * integrated in domain, e.g. for calculation of strain energy.
   *
   * @return reference to BasicMethod for preprocessing
   */
  BasicMethodsSequence &getPreProcComputeRhs() { return preProcess_Rhs; }

  /**
   * The sequence of BasicMethod is executed after residual is calculated. It
   * can be used to setup data structures, e.g. aggregate data from processors
   * or to apply essential boundary conditions.
   *
   * @return reference to BasicMethod for postprocessing
   */
  BasicMethodsSequence &getPostProcComputeRhs() { return postProcess_Rhs; }

  /**
   * @return reference to BasicMethod for preprocessing
   */
  BasicMethodsSequence &getPreProcSetOperators() { return preProcess_Mat; }

  /**
   * The sequence of BasicMethod is executed after tangent matrix is calculated.
   * It can be used to setup data structures, e.g. aggregate data from
   * processors or to apply essential boundary conditions.
   *
   * @return reference to BasicMethod for postprocessing
   */
  BasicMethodsSequence &getPostProcSetOperators() { return postProcess_Mat; }

  friend PetscErrorCode KspRhs(KSP ksp, Vec f, void *ctx);
  friend PetscErrorCode KspMat(KSP ksp, Mat A, Mat B, void *ctx);

  /**
   * @brief Clear loops
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode clearLoops();

  /** @deprecated use getSetOperator */
  DEPRECATED FEMethodsSequence &get_loops_to_do_Mat() {
    return getSetOperators();
  }

  /** @deprecated use getComputeRhs */
  DEPRECATED FEMethodsSequence &get_loops_to_do_Rhs() {
    return getComputeRhs();
  }

  /** @deprecated use getPreProcComputeRhs */
  DEPRECATED BasicMethodsSequence &get_preProcess_to_do_Rhs() {
    return getPreProcComputeRhs();
  }

  /** @deprecated use getPostProcComputeRhs */
  BasicMethodsSequence &get_postProcess_to_do_Rhs() {
    return getPostProcComputeRhs();
  }

  /** @deprecated use getPreProcSetOperators */
  DEPRECATED BasicMethodsSequence &get_preProcess_to_do_Mat() {
    return getPreProcSetOperators();
  }

  /** @deprecated use  getPostProcSetOperators */
  DEPRECATED BasicMethodsSequence &get_postProcess_to_do_Mat() {
    return getPostProcSetOperators();
  }

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
 * \brief Run over elements in the list
 * @param  ksp KSP solver
 * @param  A   matrix
 * @param  B   Preconditioned matrix
 * @param  ctx data context, i.e. KspCtx
 * @return     error code
 */
PetscErrorCode KspMat(KSP ksp, Mat A, Mat B, void *ctx);

} // namespace MoFEM

#endif // __KSPCTX_HPP__
