/** \brieg DMMoFEM.hpp
  \brief Discrete manager interface for MoFEM
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


#ifndef __DMMMOFEM_H
#define __DMMMOFEM_H

#define DM_NO_ELEMENT "DMNONEFE"

/**
  * \brief Register MoFEM problem
  * \ingroup dm
  */
PetscErrorCode DMRegister_MoFEM(const char sname[]);

/**
  * \brief Must be called by user to set MoFEM data structures
  * \ingroup dm
  */
PetscErrorCode DMMoFEMCreateMoFEM(DM dm,MoFEM::FieldInterface *m_field_ptr,const char problem_name[],const MoFEM::BitRefLevel &bit_level);

/**
  * \brief Get pointer to problem data structure
  * \ingroup dm
  */
PetscErrorCode DMMoFEMGetProblemPtr(DM dm,const MoFEM::MoFEMProblem **problem_ptr);

/**
  * \brief set squared problem
  * \ingroup dm

  It if true is assumed that matrix has the same indexing on rows and
  columns. This reduces interprocessor communication.

  */
PetscErrorCode DMMoFEMSetSquareProblem(DM dm,PetscBool square_problem);

/**
  * \brief get squared problem
  * \ingroup dm

  It if true is assumed that matrix has the same indexing on rows and
  columns. This reduces interprocessor communication.

  */
PetscErrorCode DMMoFEMGetSquareProblem(DM dm,PetscBool *square_problem);

/**
 * \brief Resolve shared entities
 * @param  dm      dm
 * @param  fe_name finite element for which shared entities are resolved
 * @return         error code
 *

 * This allows for tag reduction or tag exchange, f.e.

 \code
 ierr = DMMoFEMGetSquareProblem(dm,"SHELL_ELEMENT"); CHKERRQ(ierr);
 Tag th;
 rval = mField.get_moab().tag_get_handle("ADAPT_ORDER",th); CHKERRQ_MOAB(rval);
 ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
 // rval = pcomm->reduce_tags(th,MPI_SUM,prisms);
 rval = pcomm->exchange_tags(th,prisms);
 \endcode

 * \ingroup dm
 */
PetscErrorCode DMMoFEMResolveSharedEntities(DM dm,const char fe_name[]);

/**
 * \brief Get finite elements layout in the problem
 *
 * In layout is stored information how many elements is on each processor, for
 * more information look int petsc documentation
 * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/IS/PetscLayoutCreate.html#PetscLayoutCreate>
 *
 * @param  dm     discrete manager for this problem
 * @param  fe_name finite element name
 * @param  layout pointer to layout, for created layout user takes responsibility for destroying it.
 * @return        error code
 *
 * \ingroup dm
 */
PetscErrorCode DMMoFEMGetProblemFiniteElementLayout(DM dm,const char fe_name[],PetscLayout *layout);

/**
  * \brief add element to dm
  * \ingroup dm
  */
PetscErrorCode DMMoFEMAddElement(DM dm,const char fe_name[]);

/**
  * \brief unset element from dm
  * \ingroup dm
  */
PetscErrorCode DMMoFEMUnSetElement(DM dm,const char fe_name[]);

/**
  * \brief set local (or ghosted) vector values on mesh for partition only
  * \ingroup dm

  * \param l vector
  * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
  * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
  *
  * SCATTER_REVERSE set data to field entities from V vector.
  *
  * SCATTER_FORWARD set vector V from data field entities

  */
PetscErrorCode DMoFEMMeshToLocalVector(DM dm,Vec l,InsertMode mode,ScatterMode scatter_mode);

/**
  * \brief set ghosted vector values on all existing mesh entities
  * \ingroup dm

  * \param g vector
  * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
  * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
  *
  * SCATTER_REVERSE set data to field entities from V vector.
  *
  * SCATTER_FORWARD set vector V from data field entities

  */
PetscErrorCode DMoFEMMeshToGlobalVector(DM dm,Vec g,InsertMode mode,ScatterMode scatter_mode);

/**
  * \brief execute finite element method for each element in dm (problem)
  * \ingroup dm
  */
PetscErrorCode DMoFEMPreProcessFiniteElements(DM dm,MoFEM::FEMethod *method);

/**
  * \brief execute finite element method for each element in dm (problem)
  * \ingroup dm
  */
PetscErrorCode DMoFEMPostProcessFiniteElements(DM dm,MoFEM::FEMethod *method);

/**
 * \brief Executes FEMethod for finite elements in DM
 * @param  dm       MoFEM discrete manager
 * @param  fe_name  name of finite element
 * @param  method   pointer to \ref MoFEM::FEMethod
 * @param  low_rank lowest rank of processor
 * @param  up_rank  upper run of processor
 * @return          Error code
 * \ingroup dm
 */
PetscErrorCode DMoFEMLoopFiniteElementsUpAndLowRank(DM dm,const char fe_name[],MoFEM::FEMethod *method,int low_rank,int up_rank);

/**
 * \brief Executes FEMethod for finite elements in DM
 * @param  dm      MoFEM discrete manager
 * @param  fe_name name of element
 * @param  method  pointer to \ref MOFEM::FEMethod
 * @return         Error code
 * \ingroup dm
 */
PetscErrorCode DMoFEMLoopFiniteElements(DM dm,const char fe_name[],MoFEM::FEMethod *method);

/**
  * \brief execute method for dofs on field in problem
  * \ingroup dm
  */
PetscErrorCode DMoFEMLoopDofs(DM dm,const char field_name[],MoFEM::EntMethod *method);

// /**
//  * \brief Set compute operator for KSP solver via sub-matrix and IS
//  *
//  * @param  dm   DM
//  * @return      error code
//  *
//  * \ingroup dm
//  */
// PetscErrorCode DMMoFEMKSPSetComputeOperatorsViaSubMatrixbByIs(DM dm);

/**
  * \brief set KSP right hand side evaluation function
  * \ingroup dm
  */
PetscErrorCode DMMoFEMKSPSetComputeRHS(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only);

/**
 * \brief Set KSP opetators and push mofem finite element methods
 *
 * @param  dm        DM
 * @param  fe_name   finite element name
 * @param  method    method on the element (executed for each element in the problem which given name)
 * @param  pre_only  method for pre-process before element method
 * @param  post_only method for post-process after element method
 * @return           error code
 *
 * \ingroup dm
 */
PetscErrorCode DMMoFEMKSPSetComputeOperators(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only);

/**
  * \brief set SNES residual evaluation function
  * \ingroup dm
  */
PetscErrorCode DMMoFEMSNESSetFunction(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only);

/**
  * \brief set SNES Jacobian evaluation function
  * \ingroup dm
  */
PetscErrorCode DMMoFEMSNESSetJacobian(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only);

/**
  * \brief set TS implicit function evaluation function
  * \ingroup dm
  */
PetscErrorCode DMMoFEMTSSetIFunction(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only);

/**
  * \brief set TS Jacobian evaluation function
  * \ingroup dm
  */
PetscErrorCode DMMoFEMTSSetIJacobian(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only);

#ifdef __KSPCTX_HPP__

/**
  * \brief get MoFEM::KspCtx data structure
  * \ingroup dm
  */
PetscErrorCode DMMoFEMGetKspCtx(DM dm,MoFEM::KspCtx **ksp_ctx);

#endif //__KSPCTX_HPP__

#ifdef __SNESCTX_HPP__

/**
  * \brief get MoFEM::SnesCtx data structure
  * \ingroup dm
  */
PetscErrorCode DMMoFEMGetSnesCtx(DM dm,MoFEM::SnesCtx **snes_ctx);

/**
  * \brief Set MoFEM::SnesCtx data structure
  * \ingroup dm

  It take over pointer, do not delete it, DM will destroy pointer
  when is destroyed.

  */
PetscErrorCode DMMoFEMSetSnesCtx(DM dm,MoFEM::SnesCtx * const snes_ctx);

#endif //__SNESCTX_HPP__

#ifdef __TSCTX_HPP__

/**
  * \brief get MoFEM::TsCtx data structure
  * \ingroup dm
  */
PetscErrorCode DMMoFEMGetTsCtx(DM dm,MoFEM::TsCtx **ts_ctx);

/**
  * \brief Set MoFEM::TsCtx data structure
  * \ingroup dm

  It take over pointer, do not delete it, DM will destroy pointer
  when is destroyed.

  */
PetscErrorCode DMMoFEMSetTsCtx(DM dm,MoFEM::TsCtx * const ts_ctx);

#endif //__TSCTX_HPP__

/** sets if read mesh is partitioned
  * \ingroup dm
  */
PetscErrorCode DMMoFEMSetIsPartitioned(DM dm,PetscBool is_partitioned);

/** get if read mesh is partitioned
  * \ingroup dm
  */
PetscErrorCode DMMoFEMGetIsPartitioned(DM dm,PetscBool *is_partitioned);

/**
 * \brief Set operators for MoFEM dm
 * @param  dm
 * @return  error code
 * \ingroup dm
 */
PetscErrorCode DMSetOperators_MoFEM(DM dm);

/**
  * \brief Create dm data structure with MoFEM data structure
  * \ingroup dm
  */
PetscErrorCode DMCreate_MoFEM(DM dm);

/**
  * \brief Destroys dm with MoFEM data structure
  * \ingroup dm
  */
PetscErrorCode DMDestroy_MoFEM(DM dm);

/**
 * \brief DMShellSetCreateGlobalVector
 * \ingroup dm
 *
 * sets the routine to create a global vector
 * associated with the shell DM
 */
PetscErrorCode DMCreateGlobalVector_MoFEM(DM dm,Vec *g);

/**
 * \brief DMShellSetCreateLocalVector
 * \ingroup dm
 *
 * sets the routine to create a local vector
 * associated with the shell DM
 */
PetscErrorCode DMCreateLocalVector_MoFEM(DM dm,Vec *l);

/**
  * DMShellSetCreateMatrix
  * \ingroup dm
  *
  * sets the routine to create a matrix associated with the shell DM
  */
PetscErrorCode DMCreateMatrix_MoFEM(DM dm,Mat *M);

/**
  * Set options for MoFEM DM
  * \ingroup dm
  */
#if PETSC_VERSION_GE(3,5,3)
    PetscErrorCode DMSetFromOptions_MoFEM(PetscOptions *PetscOptionsObject,DM dm);
#else
    PetscErrorCode DMSetFromOptions_MoFEM(DM dm);
#endif

/**
  * sets up the MoFEM structures inside a DM object
  * \ingroup dm
  */
PetscErrorCode DMSetUp_MoFEM(DM dm);

/**
  * destroy the MoFEM structure
  * \ingroup dm
  */
PetscErrorCode DMDestroy_MoFEM(DM dm);

/**
  * DMShellSetGlobalToLocal
  * \ingroup dm
  *
  * the routine that begins the global to local scatter
  */
PetscErrorCode DMGlobalToLocalBegin_MoFEM(DM dm,Vec,InsertMode,Vec);

/**
  * DMShellSetGlobalToLocal
  * \ingroup dm
  *
  * the routine that begins the global to local scatter
  */
PetscErrorCode DMGlobalToLocalEnd_MoFEM(DM dm,Vec,InsertMode,Vec);

/**
  * DMShellSetLocalToGlobal
  * \ingroup dm
  *
  * the routine that begins the local to global scatter
  */
PetscErrorCode DMLocalToGlobalBegin_MoFEM(DM,Vec,InsertMode,Vec);

/**
  * DMShellSetLocalToGlobal
  * \ingroup dm
  *
  * the routine that ends the local to global scatter
  */
PetscErrorCode DMLocalToGlobalBegin_MoFEM(DM,Vec,InsertMode,Vec);

/**
  * DMShellSetLocalToGlobal
  * \ingroup dm
  *
  * the routine that ends the local to global scatter
  */
PetscErrorCode DMLocalToGlobalEnd_MoFEM(DM,Vec,InsertMode,Vec);


namespace MoFEM {

  static const int DMCTX_INTERFACE = 1<<0;
  static const MOFEMuuid IDD_DMCTX = MOFEMuuid(BitIntefaceId(DMCTX_INTERFACE));

  /**
   * \brief PETSc  Discrete Manager data structure
   *
   * This structure should not be accessed or modified by user. Is not available
   * from outside MoFEM DM manager. However user can inherit dat class and
   * add data for additional functionality.
   *
   * This is part of implementation for PETSc interface, see more details in
   * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/index.html>
   *
   */
  struct DMCtx: public MoFEM::UnknownInterface {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,UnknownInterface** iface);

    FieldInterface *mField_ptr; 		///< MoFEM interface
    PetscBool isProblemBuild;      ///< True if problem is build
    string problemName;			        ///< Problem name

    KspCtx *kspCtx;			  ///< data structure KSP
    SnesCtx *snesCtx;			///< data structure SNES
    TsCtx	*tsCtx;				   ///< data structure for TS solver

    //options
    PetscBool isPartitioned;		///< true if read mesh is on parts
    PetscBool isSquareMatrix;		///< true if rows equals to cols
    PetscInt verbosity;			    ///< verbosity

    int rAnk,sIze;

    //pointer to data structures
    const MoFEMProblem *problemPtr;	  ///< pinter to problem data structure

    DMCtx();
    virtual ~DMCtx();

    int referenceNumber;

  };


}

#endif //__DMMMOFEM_H

/***************************************************************************//**
 * \defgroup dm Discreet manager
 * \ingroup mofem
 ******************************************************************************/
