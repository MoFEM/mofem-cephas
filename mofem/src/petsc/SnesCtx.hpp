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
    

    MoFEM::Interface &mField;   ///< database Interface
    moab::Interface &moab;      ///< moab Interface

    std::string problemName;    ///< problem name
    MoFEMTypes bH;              ///< If set to MF_EXIST check if element exist, default MF_EXIST
    bool zeroPreCondMatrixB;    ///< If true zero matrix, otherwise user need to do it, default true

    struct LoopPairType: public std::pair<std::string,FEMethod*> {
      LoopPairType(std::string name,FEMethod *ptr):
      std::pair<std::string,FEMethod*>(name,ptr) {}
      LoopPairType(std::string name,boost::shared_ptr<FEMethod> ptr):
      std::pair<std::string,FEMethod*>(name,ptr.get()),
      fePtr(ptr) {}
      virtual ~LoopPairType() {}
    private:
      boost::shared_ptr<FEMethod> fePtr;
    };
    typedef LoopPairType loop_pair_type;
    typedef std::vector<loop_pair_type > loops_to_do_type;

    loops_to_do_type loops_to_do_Mat;
    loops_to_do_type loops_to_do_Rhs;

    struct BasicMethodPtr {
      BasicMethodPtr(BasicMethod *ptr):
      rawPtr(ptr) {}
      BasicMethodPtr(boost::shared_ptr<BasicMethod> ptr):
      rawPtr(ptr.get()),
      bmPtr(ptr) {}
      BasicMethodPtr(boost::shared_ptr<FEMethod> ptr):
      rawPtr(ptr.get()),
      bmPtr(ptr) {}
      inline BasicMethod& operator*() const { return *rawPtr; };
      inline BasicMethod* operator->() const { return rawPtr; }
    private:
      BasicMethod* rawPtr;
      boost::shared_ptr<BasicMethod> bmPtr;
    };
    typedef std::vector<BasicMethodPtr> basic_method_to_do;

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
    bH(MF_EXIST),
    zeroPreCondMatrixB(true) {
      PetscLogEventRegister("LoopSNESRhs",0,&USER_EVENT_SnesRhs);
      PetscLogEventRegister("LoopSNESMat",0,&USER_EVENT_SnesMat);
    }

    virtual ~SnesCtx() {}

    /**
    * @return return reference to vector with FEMethod to calculate tangent matrix
    */
    loops_to_do_type& get_loops_to_do_Mat() { return loops_to_do_Mat; }

    /**
    * @return return vector to vector with FEMethod to calculate residual
    */
    loops_to_do_type& get_loops_to_do_Rhs() { return loops_to_do_Rhs; }

    /**
    * The sequence of BasicMethod is executed before residual is calculated. It can be
    * used to setup data structures, e.g. zero global variable which is integrated in
    * domain, e.g. for calculation of strain energy.
    *
    * @return reference to BasicMethod for preprocessing
    */
    basic_method_to_do& get_preProcess_to_do_Rhs() { return preProcess_Rhs; }

    /**
    * The sequence of BasicMethod is executed after residual is calculated. It can be
    * used to setup data structures, e.g. aggregate data from processors or to apply
    * essential boundary conditions.
    *
    * @return reference to BasicMethod for postprocessing
    */
    basic_method_to_do& get_postProcess_to_do_Rhs() { return postProcess_Rhs; }

    /**
    * @return reference to BasicMethod for preprocessing
    */
    basic_method_to_do& get_preProcess_to_do_Mat() { return preProcess_Mat; }

    /**
    * The sequence of BasicMethod is executed after tangent matrix is calculated. It can be
    * used to setup data structures, e.g. aggregate data from processors or to apply
    * essential boundary conditions.
    *
    * @return reference to BasicMethod for postprocessing
    */
    basic_method_to_do& get_postProcess_to_do_Mat() { return postProcess_Mat; }

    friend PetscErrorCode SnesRhs(SNES snes,Vec x,Vec f,void *ctx);
    friend PetscErrorCode SnesMat(SNES snes,Vec x,Mat A,Mat B,void *ctx);

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
  * For more information pleas look to PETSc manual, i.e. SNESSetFunction
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

}

#endif // __SNESCTX_HPP__
