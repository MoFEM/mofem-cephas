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

    MoFEM::Interface &mField;   ///< database Interface
    moab::Interface &moab;      ///< moab Interface

    std::string problemName;    ///< problem name
    MoFEMTypes bH;              ///< If set to MF_EXIST check if element exist, default MF_EXIST
    bool zeroPreCondMatrixB;    ///< If true zero matrix, otherwise user need to do it, default true
    MatAssemblyType typeOfAssembly; ///< type of assembly at the end

    /// \deprecated use PairNameFEMethodPtr
    DEPRECATED  typedef MoFEM::PairNameFEMethodPtr loop_pair_type;

    /// \deprecated use FEMethodsSequence
    DEPRECATED typedef MoFEM::FEMethodsSequence loops_to_do_type;

    /// \deprecated use BasicMethodsSequence
    DEPRECATED typedef MoFEM::BasicMethodsSequence basic_method_to_do;

    typedef MoFEM::PairNameFEMethodPtr PairNameFEMethodPtr;
    typedef MoFEM::FEMethodsSequence FEMethodsSequence;
    typedef MoFEM::BasicMethodsSequence BasicMethodsSequence;

    FEMethodsSequence loops_to_do_Mat;    ///< Sequence of finite elements instances assembiling tangent matrix
    FEMethodsSequence loops_to_do_Rhs;    ///< Sequence of finite elements instances assembiling residual vector
    BasicMethodsSequence preProcess_Mat;  ///< Sequence of methods run before tangent matrix is assembled
    BasicMethodsSequence postProcess_Mat; ///< Sequence of methods run after tangent matrix is assembled
    BasicMethodsSequence preProcess_Rhs;  ///< Sequence of methods run before residual is assembled
    BasicMethodsSequence postProcess_Rhs; ///< Sequence of methods run after residual is assembled

    /**
     * \brief Copy seqences from other SNES contex
     * @param  snes_ctx SNES contex from which Sequence is copied from
     * @return          error code
     */
    PetscErrorCode copyLoops(const SnesCtx &snes_ctx) {
      MoFEMFunctionBeginHot;
      loops_to_do_Mat = snes_ctx.loops_to_do_Mat;
      loops_to_do_Rhs = snes_ctx.loops_to_do_Rhs;
      preProcess_Mat = snes_ctx.preProcess_Mat;
      postProcess_Mat = snes_ctx.postProcess_Mat;
      preProcess_Rhs = snes_ctx.preProcess_Rhs;
      postProcess_Rhs = snes_ctx.postProcess_Rhs;
      MoFEMFunctionReturnHot(0);
    }

    PetscLogEvent USER_EVENT_SnesRhs; ///< Log events to assemble residual
    PetscLogEvent USER_EVENT_SnesMat; ///< Log events to assemble tangent matrix

    SnesCtx(Interface &m_field,const std::string &problem_name):
    mField(m_field),
    moab(m_field.get_moab()),
    problemName(problem_name),
    bH(MF_EXIST),
    zeroPreCondMatrixB(true),
    typeOfAssembly(MAT_FINAL_ASSEMBLY) {
      PetscLogEventRegister("LoopSNESRhs",0,&USER_EVENT_SnesRhs);
      PetscLogEventRegister("LoopSNESMat",0,&USER_EVENT_SnesMat);
    }

    virtual ~SnesCtx() {
    }

    /**
    * @return return reference to vector with FEMethod to calculate tangent matrix
    */
    FEMethodsSequence& get_loops_to_do_Mat() { return loops_to_do_Mat; }

    /**
    * @return return vector to vector with FEMethod to calculate residual
    */
    FEMethodsSequence& get_loops_to_do_Rhs() { return loops_to_do_Rhs; }

    /**
    * The sequence of BasicMethod is executed before residual is calculated. It can be
    * used to setup data structures, e.g. zero global variable which is integrated in
    * domain, e.g. for calculation of strain energy.
    *
    * @return reference to BasicMethod for preprocessing
    */
    BasicMethodsSequence& get_preProcess_to_do_Rhs() { return preProcess_Rhs; }

    /**
    * The sequence of BasicMethod is executed after residual is calculated. It can be
    * used to setup data structures, e.g. aggregate data from processors or to apply
    * essential boundary conditions.
    *
    * @return reference to BasicMethod for postprocessing
    */
    BasicMethodsSequence& get_postProcess_to_do_Rhs() { return postProcess_Rhs; }

    /**
    * @return reference to BasicMethod for preprocessing
    */
    BasicMethodsSequence& get_preProcess_to_do_Mat() { return preProcess_Mat; }

    /**
    * The sequence of BasicMethod is executed after tangent matrix is calculated. It can be
    * used to setup data structures, e.g. aggregate data from processors or to apply
    * essential boundary conditions.
    *
    * @return reference to BasicMethod for postprocessing
    */
    BasicMethodsSequence& get_postProcess_to_do_Mat() { return postProcess_Mat; }

    friend PetscErrorCode SnesRhs(SNES snes,Vec x,Vec f,void *ctx);
    friend PetscErrorCode SnesMat(SNES snes,Vec x,Mat A,Mat B,void *ctx);

    friend PetscErrorCode SNESMoFEMSetAssmblyType(SNES snes,MatAssemblyType type);
    friend PetscErrorCode SNESMoFEMSetBehavior(SNES snes,MoFEMTypes bh);
  };

  /**
  * \brief This is MoFEM implementation for the right hand side (residual vector) evaluation in SNES solver
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

  /**
  * \brief This is MoFEM implementation for the left hand side (tangent matrix) evaluation in SNES solver
  *
  * For more information pleas look to PETSc manual, i.e. SNESSetJacobian
  * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESSetJacobian.html#SNESSetJacobian>
  *
  * @param  snes SNES solver
  * @param  x    Solution vector at current iteration
  * @param  A    Tangent matrix
  * @param  B    Preconditioner tangent matrix
  * @param  ctx  Pointer to context thata, i.e. SnesCtx
  * @return      Error code
  */
  PetscErrorCode SnesMat(SNES snes,Vec x,Mat A,Mat B,void *ctx);

  /**
   * \brief Set assembly type at the end of SnesMat
   *
   * \note Note that tangent matrix need have to have final assembly, you would
   * use flush assembly in special case that you call SnesMat form other function
   * set to SNESSetJacobian
   *
   * @param  snes
   * @param  type type of assembly, either MAT_FLUSH_ASSEMBLY or MAT_FINAL_ASSEMBLY
   * @return      error code
   */
  PetscErrorCode SNESMoFEMSetAssmblyType(SNES snes,MatAssemblyType type);

  /**
   * \brief Set behavior if finite element in sequence does not exist
   * @param  snes
   * @param  bh   If set to MF_EXIST check if element exist, default MF_EXIST. Otherwise set MF_ZERO
   * @return      error code
   */
  PetscErrorCode SNESMoFEMSetBehavior(SNES snes,MoFEMTypes bh);


}

#endif // __SNESCTX_HPP__
