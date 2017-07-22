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

  ErrorCode rval;
  

  MoFEM::Interface &mField;
  moab::Interface &moab;

  std::string problemName;
  MoFEMTypes bH; ///< If set to MF_EXIST check if element exist

  struct LoopPairType: public std::pair<std::string,FEMethod*> {
    LoopPairType(std::string name,FEMethod *ptr):
    std::pair<std::string,FEMethod*>(name,ptr) {}
    LoopPairType(std::string name,boost::shared_ptr<FEMethod> ptr):
    std::pair<std::string,FEMethod*>(name,ptr.get()),
    fePtr(ptr) {
    }
    virtual ~LoopPairType() {}
  private:
    boost::shared_ptr<FEMethod> fePtr;
  };
  typedef LoopPairType loop_pair_type;

  typedef std::vector<loop_pair_type > loops_to_do_type;
  loops_to_do_type loops_to_do_IJacobian;
  loops_to_do_type loops_to_do_IFunction;
  loops_to_do_type loops_to_do_Monitor;

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

  basic_method_to_do preProcess_IJacobian;
  basic_method_to_do postProcess_IJacobian;
  basic_method_to_do preProcess_IFunction;
  basic_method_to_do postProcess_IFunction;
  basic_method_to_do preProcess_Monitor;
  basic_method_to_do postProcess_Monitor;

  PetscLogEvent USER_EVENT_TsCtxRHSFunction;
  PetscLogEvent USER_EVENT_TsCtxRHSJacobian;
  PetscLogEvent USER_EVENT_TsCtxIFunction;
  PetscLogEvent USER_EVENT_TsCtxIJacobian;
  PetscLogEvent USER_EVENT_TsCtxMonitor;

  bool zeroMatrix;
  TsCtx(MoFEM::Interface &m_field,const std::string &problem_name):
    mField(m_field),
    moab(m_field.get_moab()),
    problemName(problem_name),
    bH(MF_EXIST),
    zeroMatrix(true) {
    PetscLogEventRegister("LoopTsIFunction",0,&USER_EVENT_TsCtxIFunction);
    PetscLogEventRegister("LoopTsIJacobian",0,&USER_EVENT_TsCtxIJacobian);
    PetscLogEventRegister("LoopTsRHSFunction",0,&USER_EVENT_TsCtxRHSFunction);
    PetscLogEventRegister("LoopTsRHSJacobian",0,&USER_EVENT_TsCtxRHSJacobian);
    PetscLogEventRegister("LoopTsMonitor",0,&USER_EVENT_TsCtxMonitor);
  }

  loops_to_do_type& get_loops_to_do_IFunction() { return loops_to_do_IFunction; }
  loops_to_do_type& get_loops_to_do_IJacobian() { return loops_to_do_IJacobian; }
  loops_to_do_type& get_loops_to_do_Monitor() { return loops_to_do_Monitor; }

  basic_method_to_do& get_preProcess_to_do_IFunction() { return preProcess_IFunction; }
  basic_method_to_do& get_postProcess_to_do_IFunction() { return postProcess_IFunction; }
  basic_method_to_do& get_preProcess_to_do_IJacobian() { return preProcess_IJacobian; }
  basic_method_to_do& get_postProcess_to_do_IJacobian() { return postProcess_IJacobian; }
  basic_method_to_do& get_preProcess_to_do_Monitor() { return preProcess_Monitor; }
  basic_method_to_do& get_postProcess_to_do_Monitor() { return postProcess_Monitor; }

  friend PetscErrorCode f_TSSetIFunction(TS ts,PetscReal t,Vec u,Vec u_t,Vec F,void *ctx);
  friend PetscErrorCode f_TSSetIJacobian(TS ts,PetscReal t,Vec u,Vec U_t,PetscReal a,Mat A,Mat B,void *ctx);
  friend PetscErrorCode f_TSMonitorSet(TS ts,PetscInt step,PetscReal t,Vec u,void *ctx);

};

PetscErrorCode f_TSSetIFunction(TS ts,PetscReal t,Vec u,Vec u_t,Vec F,void *ctx);
PetscErrorCode f_TSSetIJacobian(TS ts,PetscReal t,Vec u,Vec u_t,PetscReal a,Mat A,Mat B,void *ctx);
PetscErrorCode f_TSMonitorSet(TS ts,PetscInt step,PetscReal t,Vec u,void *ctx);


}

#endif // __TSCTX_HPP__
