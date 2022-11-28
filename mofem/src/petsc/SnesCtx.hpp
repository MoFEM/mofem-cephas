/** \file SnesCtx.hpp
 * \brief Context for PETSc SNES, i.e. nonlinear solver
 */

#ifndef __SNESCTX_HPP__
#define __SNESCTX_HPP__

namespace MoFEM {

/** \brief Interface for nonlinear (SNES) solver
 * \ingroup mofem_petsc_solvers
 */
struct SnesCtx {

  MoFEM::Interface &mField; ///< database Interface
  moab::Interface &moab;    ///< moab Interface

  std::string problemName; ///< problem name
  MoFEMTypes
      bH; ///< If set to MF_EXIST check if element exist, default MF_EXIST
  bool zeroPreCondMatrixB; ///< If true zero matrix, otherwise user need to do
                           ///< it, default true
  MatAssemblyType typeOfAssembly; ///< type of assembly at the end
  bool vErify;                    ///< If true verify vector

  typedef MoFEM::PairNameFEMethodPtr PairNameFEMethodPtr;
  typedef MoFEM::FEMethodsSequence FEMethodsSequence;
  typedef MoFEM::BasicMethodsSequence BasicMethodsSequence;

  FEMethodsSequence loops_to_do_Mat; ///< Sequence of finite elements instances
                                     ///< assembling tangent matrix
  FEMethodsSequence loops_to_do_Rhs; ///< Sequence of finite elements instances
                                     ///< assembling residual vector
  BasicMethodsSequence preProcess_Mat;  ///< Sequence of methods run before
                                        ///< tangent matrix is assembled
  BasicMethodsSequence postProcess_Mat; ///< Sequence of methods run after
                                        ///< tangent matrix is assembled
  BasicMethodsSequence
      preProcess_Rhs; ///< Sequence of methods run before residual is assembled
  BasicMethodsSequence
      postProcess_Rhs; ///< Sequence of methods run after residual is assembled

  SnesCtx(Interface &m_field, const std::string &problem_name)
      : mField(m_field), moab(m_field.get_moab()), problemName(problem_name),
        bH(MF_EXIST), zeroPreCondMatrixB(true),
        typeOfAssembly(MAT_FINAL_ASSEMBLY), vErify(false) {
    PetscLogEventRegister("LoopSNESRhs", 0, &MOFEM_EVENT_SnesRhs);
    PetscLogEventRegister("LoopSNESMat", 0, &MOFEM_EVENT_SnesMat);
    if (!LogManager::checkIfChannelExist("SNES_WORLD")) {
      auto core_log = logging::core::get();
      core_log->add_sink(
          LogManager::createSink(LogManager::getStrmWorld(), "SNES_WORLD"));
      LogManager::setLog("SNES_WORLD");
      MOFEM_LOG_TAG("SNES_WORLD", "SNES");
    }
  }

  virtual ~SnesCtx() = default;

  /**
   * @return return reference to vector with FEMethod to calculate tangent
   * matrix
   */
  FEMethodsSequence &get_loops_to_do_Mat() { return loops_to_do_Mat; }

  /**
   * @return return vector to vector with FEMethod to calculate residual
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

  /**
   * \brief Copy sequences from other SNES contex
   * @param  snes_ctx SNES contex from which Sequence is copied from
   * @return          error code
   */
  MoFEMErrorCode copyLoops(const SnesCtx &snes_ctx);

  /**
   * @brief Clear loops
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode clearLoops();

  friend PetscErrorCode SnesRhs(SNES snes, Vec x, Vec f, void *ctx);
  friend PetscErrorCode SnesMat(SNES snes, Vec x, Mat A, Mat B, void *ctx);

  friend MoFEMErrorCode SNESMoFEMSetAssmblyType(SNES snes,
                                                MatAssemblyType type);
  friend MoFEMErrorCode SnesMoFEMSetBehavior(SNES snes, MoFEMTypes bh);

private:
  boost::movelib::unique_ptr<bool> vecAssembleSwitch;
  boost::movelib::unique_ptr<bool> matAssembleSwitch;
  PetscLogEvent MOFEM_EVENT_SnesRhs; ///< Log events to assemble residual
  PetscLogEvent MOFEM_EVENT_SnesMat; ///< Log events to assemble tangent matrix
};

/**
 * \brief This is MoFEM implementation for the right hand side (residual vector)
 * evaluation in SNES solver
 *
 * For more information pleas look to PETSc manual, i.e. SNESSetFunction
 * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESSetFunction.html>
 *
 * @param  snes SNES solver
 * @param  x    Solution vector at current iteration
 * @param  f    The right hand side vector
 * @param  ctx  Pointer to context i.e. SnesCtx
 * @return      Error code
 */
PetscErrorCode SnesRhs(SNES snes, Vec x, Vec f, void *ctx);

/**
 * \brief This is MoFEM implementation for the left hand side (tangent matrix)
 * evaluation in SNES solver
 *
 * For more information pleas look to PETSc manual, i.e. SNESSetJacobian
 * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESSetJacobian.html#SNESSetJacobian>
 *
 * @param  snes SNES solver
 * @param  x    Solution vector at current iteration
 * @param  A    Tangent matrix
 * @param  B    Preconditioner tangent matrix
 * @param  ctx  Pointer to context i.e. SnesCtx
 * @return      Error code
 */
PetscErrorCode SnesMat(SNES snes, Vec x, Mat A, Mat B, void *ctx);

/**
 * \brief Set assembly type at the end of SnesMat
 *
 * \note Note that tangent matrix need have to have final assembly, you would
 * use flush assembly in special case that you call SnesMat form other function
 * set to SNESSetJacobian
 *
 * @param  snes
 * @param  type type of assembly, either MAT_FLUSH_ASSEMBLY or
 * MAT_FINAL_ASSEMBLY
 * @return      error code
 */
MoFEMErrorCode SnesMoFEMSetAssemblyType(SNES snes, MatAssemblyType type);

/**
 * \brief Set behavior if finite element in sequence does not exist
 * @param  snes
 * @param  bh   If set to MF_EXIST check if element exist, default MF_EXIST.
 * Otherwise set MF_ZERO
 * @return      error code
 */
MoFEMErrorCode SnesMoFEMSetBehavior(SNES snes, MoFEMTypes bh);

/**
 * @brief Sens monitor printing residual field by field
 *
 */
MoFEMErrorCode MoFEMSNESMonitorFields(SNES snes, PetscInt its, PetscReal fgnorm,
                                      SnesCtx *snes_ctx);

} // namespace MoFEM

#endif // __SNESCTX_HPP__
