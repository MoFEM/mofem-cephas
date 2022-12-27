/** \file DMMoFEM.hpp
  \brief Discrete manager interface for MoFEM
  */



#ifndef __DMMMOFEM_H
#define __DMMMOFEM_H

#define DM_NO_ELEMENT "DMNONEFE"

namespace MoFEM {

/**
 * \brief Register MoFEM problem
 * \ingroup dm
 */
PetscErrorCode DMRegister_MoFEM(const char sname[]);

/**
 * \brief Must be called by user to set MoFEM data structures
 * \ingroup dm
 */
PetscErrorCode DMMoFEMCreateMoFEM(
    DM dm, MoFEM::Interface *m_field_ptr, const char problem_name[],
    const MoFEM::BitRefLevel bit_level,
    const MoFEM::BitRefLevel bit_mask = MoFEM::BitRefLevel().set());

/**
 * @brief Duplicate internal data struture
 * 
 * @param dm 
 * @param dm_duplicate 
 * @return PetscErrorCode 
 */
PetscErrorCode DMMoFEMDuplicateDMCtx(DM dm, DM dm_duplicate);

/**
 * \brief Must be called by user to set Sub DM MoFEM data structures
 * \ingroup dm
 */
PetscErrorCode DMMoFEMCreateSubDM(DM subdm, DM dm, const char problem_name[]);

/**
 * \brief Get pointer to MoFEM::Interface
 * @param  dm          Distributed mesh manager
 * @param  m_field_ptr Pointer to pointer of field interface
 * @return             Error code
 * \ingroup dm
 */
PetscErrorCode DMoFEMGetInterfacePtr(DM dm, MoFEM::Interface **m_field_ptr);

/**
 * \brief Get pointer to problem data structure
 * \ingroup dm
 */
PetscErrorCode DMMoFEMGetProblemPtr(DM dm, const MoFEM::Problem **problem_ptr);

/**
 * If this is set to PETSC_TRUE problem is deleted with DM
 * @param  dm             the DM object
 * @param  destroy        if PETSC_TRUE problem is destroyed
 * @return                error code
 */
PetscErrorCode DMMoFEMSetDestroyProblem(DM dm, PetscBool destroy_problem);

/**
 * Get if problem will be destroyed with DM
 * @param  dm             the DM object
 * @param  destroy        return if PETSC_TRUE problem is destroyed
 * @return                error code
 */
PetscErrorCode DMMoFEMGetDestroyProblem(DM dm, PetscBool *destroy_problem);

/**
  * \brief set squared problem
  * \ingroup dm

  It if true is assumed that matrix has the same indexing on rows and
  columns. This reduces interprocessor communication.

  */
PetscErrorCode DMMoFEMSetSquareProblem(DM dm, PetscBool square_problem);

/**
  * \brief get squared problem
  * \ingroup dm

  It if true is assumed that matrix has the same indexing on rows and
  columns. This reduces interprocessor communication.

  */
PetscErrorCode DMMoFEMGetSquareProblem(DM dm, PetscBool *square_problem);

/**
 * \brief Resolve shared entities
 *
 * @param  dm      dm
 * @param  fe_name finite element for which shared entities are resolved
 * @return         error code
 *
 * \note This function is valid for parallel algebra and serial mesh. It
 * should be run collectively, i.e. on all processors.
 *
 * This allows for tag reduction or tag exchange, f.e.
 *
 * \code
 * CHKERR DMMoFEMResolveSharedFiniteElements(dm,"SHELL_ELEMENT");
 * Tag th;
 * CHKERR mField.get_moab().tag_get_handle("ADAPT_ORDER",th);
 * ParallelComm* pcomm =
 * ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
 * // CHKERR pcomm->reduce_tags(th,MPI_SUM,prisms);
 * CHKERR pcomm->exchange_tags(th,prisms);
 * \endcode
 *
 * \ingroup dm
 */
PetscErrorCode DMMoFEMResolveSharedFiniteElements(DM dm, const char fe_name[]);

/**
 * @deprecated Use DMMoFEMResolveSharedFiniteElements
 */
DEPRECATED PetscErrorCode DMMoFEMResolveSharedEntities(DM dm,
                                                       const char fe_name[]);

/**
 * \brief Get finite elements layout in the problem
 *
 * In layout is stored information how many elements is on each processor, for
 * more information look int petsc documentation
 * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/IS/PetscLayoutCreate.html#PetscLayoutCreate>
 *
 * @param  dm     discrete manager for this problem
 * @param  fe_name finite element name
 * @param  layout pointer to layout, for created layout user takes
 * responsibility for destroying it.
 * @return        error code
 *
 * \ingroup dm
 */
PetscErrorCode DMMoFEMGetProblemFiniteElementLayout(DM dm, const char fe_name[],
                                                    PetscLayout *layout);

/**
 * \brief add element to dm
 * \ingroup dm
 *
 * \note add_file is a collective, should be executed on all processors.
 * Otherwise could lead to deadlock.
 *
 */
PetscErrorCode DMMoFEMAddElement(DM dm, std::string fe_name);

/**
 * \brief add element to dm
 * \ingroup dm
 *
 * \note add_file is a collective, should be executed on all processors.
 * Otherwise could lead to deadlock.
 *
 */
PetscErrorCode DMMoFEMAddElement(DM dm, std::vector<std::string> fe_name);

/**
 * \brief unset element from dm
 * \ingroup dm
 */
PetscErrorCode DMMoFEMUnSetElement(DM dm, std::string fe_name);

/**
  * \brief set local (or ghosted) vector values on mesh for partition only
  * \ingroup dm

  * \param l vector
  * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
  * \param scatter_mode see petsc manual for ScatterMode (The available modes
  are: SCATTER_FORWARD or SCATTER_REVERSE)
  *
  * SCATTER_REVERSE set data to field entities from V vector.
  *
  * SCATTER_FORWARD set vector V from data field entities

  */
PetscErrorCode DMoFEMMeshToLocalVector(DM dm, Vec l, InsertMode mode,
                                       ScatterMode scatter_mode);

/**
  * \brief set ghosted vector values on all existing mesh entities
  * \ingroup dm

  * \param g vector
  * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
  * \param scatter_mode see petsc manual for ScatterMode (The available modes
  are: SCATTER_FORWARD or SCATTER_REVERSE)
  *
  * SCATTER_REVERSE set data to field entities from V vector.
  *
  * SCATTER_FORWARD set vector V from data field entities

  */
PetscErrorCode DMoFEMMeshToGlobalVector(DM dm, Vec g, InsertMode mode,
                                        ScatterMode scatter_mode);

/**
 * \brief execute finite element method for each element in dm (problem)
 * \ingroup dm
 */
PetscErrorCode DMoFEMPreProcessFiniteElements(DM dm, MoFEM::FEMethod *method);

/**
 * \brief execute finite element method for each element in dm (problem)
 * \ingroup dm
 */
PetscErrorCode DMoFEMPostProcessFiniteElements(DM dm, MoFEM::FEMethod *method);

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
PetscErrorCode DMoFEMLoopFiniteElementsUpAndLowRank(
    DM dm, const char fe_name[], MoFEM::FEMethod *method, int low_rank,
    int up_rank, CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr());

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
PetscErrorCode DMoFEMLoopFiniteElementsUpAndLowRank(
    DM dm, const std::string fe_name, boost::shared_ptr<MoFEM::FEMethod> method,
    int low_rank, int up_rank,
    CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr());

/**
 * \brief Executes FEMethod for finite elements in DM
 * @param  dm      MoFEM discrete manager
 * @param  fe_name name of element
 * @param  method  pointer to \ref MOFEM::FEMethod
 * @return         Error code
 * \ingroup dm
 */
PetscErrorCode
DMoFEMLoopFiniteElements(DM dm, const char fe_name[], MoFEM::FEMethod *method,
                         CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr());

/**
 * \brief Executes FEMethod for finite elements in DM
 * @param  dm      MoFEM discrete manager
 * @param  fe_name name of element
 * @param  method  pointer to \ref MOFEM::FEMethod
 * @return         Error code
 * \ingroup dm
 */
PetscErrorCode
DMoFEMLoopFiniteElements(DM dm, const std::string fe_name,
                         boost::shared_ptr<MoFEM::FEMethod> method,
                         CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr());

/**
 * \brief execute method for dofs on field in problem
 * \ingroup dm
 */
PetscErrorCode DMoFEMLoopDofs(DM dm, const char field_name[],
                              MoFEM::DofMethod *method);

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
PetscErrorCode DMMoFEMKSPSetComputeRHS(DM dm, const char fe_name[],
                                       MoFEM::FEMethod *method,
                                       MoFEM::BasicMethod *pre_only,
                                       MoFEM::BasicMethod *post_only);

/**
 * \brief set KSP right hand side evaluation function
 * \ingroup dm
 */
PetscErrorCode
DMMoFEMKSPSetComputeRHS(DM dm, const std::string fe_name,
                        boost::shared_ptr<MoFEM::FEMethod> method,
                        boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                        boost::shared_ptr<MoFEM::BasicMethod> post_only);

/**
 * \brief Set KSP operators and push mofem finite element methods
 *
 * @param  dm        DM
 * @param  fe_name   finite element name
 * @param  method    method on the element (executed for each element in the
 * problem which given name)
 * @param  pre_only  method for pre-process before element method
 * @param  post_only method for post-process after element method
 * @return           error code
 *
 * \ingroup dm
 */
PetscErrorCode DMMoFEMKSPSetComputeOperators(DM dm, const char fe_name[],
                                             MoFEM::FEMethod *method,
                                             MoFEM::BasicMethod *pre_only,
                                             MoFEM::BasicMethod *post_only);

/**
 * \brief Set KSP operators and push mofem finite element methods
 *
 * @param  dm        DM
 * @param  fe_name   finite element name
 * @param  method    method on the element (executed for each element in the
 * problem which given name)
 * @param  pre_only  method for pre-process before element method
 * @param  post_only method for post-process after element method
 * @return           error code
 *
 * \ingroup dm
 */
PetscErrorCode
DMMoFEMKSPSetComputeOperators(DM dm, const std::string fe_name,
                              boost::shared_ptr<MoFEM::FEMethod> method,
                              boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                              boost::shared_ptr<MoFEM::BasicMethod> post_only);

/**
 * \brief set SNES residual evaluation function
 * \ingroup dm
 */
PetscErrorCode DMMoFEMSNESSetFunction(DM dm, const char fe_name[],
                                      MoFEM::FEMethod *method,
                                      MoFEM::BasicMethod *pre_only,
                                      MoFEM::BasicMethod *post_only);

/**
 * \brief set SNES residual evaluation function
 * \ingroup dm
 */
PetscErrorCode
DMMoFEMSNESSetFunction(DM dm, const std::string fe_name,
                       boost::shared_ptr<MoFEM::FEMethod> method,
                       boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                       boost::shared_ptr<MoFEM::BasicMethod> post_only);

/**
 * \brief set SNES Jacobian evaluation function
 * \ingroup dm
 */
PetscErrorCode DMMoFEMSNESSetJacobian(DM dm, const char fe_name[],
                                      MoFEM::FEMethod *method,
                                      MoFEM::BasicMethod *pre_only,
                                      MoFEM::BasicMethod *post_only);

/**
 * \brief set SNES Jacobian evaluation function
 * \ingroup dm
 */
PetscErrorCode
DMMoFEMSNESSetJacobian(DM dm, const std::string fe_name,
                       boost::shared_ptr<MoFEM::FEMethod> method,
                       boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                       boost::shared_ptr<MoFEM::BasicMethod> post_only);

/**
 * \brief set TS implicit function evaluation function
 * \ingroup dm
 */
PetscErrorCode DMMoFEMTSSetIFunction(DM dm, const char fe_name[],
                                     MoFEM::FEMethod *method,
                                     MoFEM::BasicMethod *pre_only,
                                     MoFEM::BasicMethod *post_only);

/**
 * \brief set TS implicit function evaluation function
 * \ingroup dm
 */
PetscErrorCode
DMMoFEMTSSetIFunction(DM dm, const std::string fe_name,
                      boost::shared_ptr<MoFEM::FEMethod> method,
                      boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                      boost::shared_ptr<MoFEM::BasicMethod> post_only);

/**
 * \brief set TS Jacobian evaluation function
 * \ingroup dm
 */
PetscErrorCode
DMMoFEMTSSetIJacobian(DM dm, const std::string fe_name,
                      boost::shared_ptr<MoFEM::FEMethod> method,
                      boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                      boost::shared_ptr<MoFEM::BasicMethod> post_only);

/**
 * \brief set TS Jacobian evaluation function
 * \ingroup dm
 */
PetscErrorCode DMMoFEMTSSetIJacobian(DM dm, const char fe_name[],
                                     MoFEM::FEMethod *method,
                                     MoFEM::BasicMethod *pre_only,
                                     MoFEM::BasicMethod *post_only);

/**
 * @brief set TS the right hand side function
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetRHSFunction.html#TSSetRHSFunction>See
 * petsc documentation</a>
 *
 * @param dm
 * @param fe_name
 * @param method
 * @param pre_only
 * @param post_only
 * @return PetscErrorCode
 */
PetscErrorCode
DMMoFEMTSSetRHSFunction(DM dm, const std::string fe_name,
                        boost::shared_ptr<MoFEM::FEMethod> method,
                        boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                        boost::shared_ptr<MoFEM::BasicMethod> post_only);

PetscErrorCode DMMoFEMTSSetRHSFunction(DM dm, const char fe_name[],
                                       MoFEM::FEMethod *method,
                                       MoFEM::BasicMethod *pre_only,
                                       MoFEM::BasicMethod *post_only);

/**
 * @brief set TS the right hand side jacobian
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetRHSJacobian.html>See
 * petsc documentation</a>
 *
 * @param dm
 * @param fe_name
 * @param method
 * @param pre_only
 * @param post_only
 * @return PetscErrorCode
 */
PetscErrorCode
DMMoFEMTSSetRHSJacobian(DM dm, const std::string fe_name,
                        boost::shared_ptr<MoFEM::FEMethod> method,
                        boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                        boost::shared_ptr<MoFEM::BasicMethod> post_only);

PetscErrorCode DMMoFEMTSSetRHSJacobian(DM dm, const char fe_name[],
                                       MoFEM::FEMethod *method,
                                       MoFEM::BasicMethod *pre_only,
                                       MoFEM::BasicMethod *post_only);

/**
 * \brief set TS implicit function evaluation function
 * \ingroup dm
 */
PetscErrorCode
DMMoFEMTSSetI2Function(DM dm, const std::string fe_name,
                       boost::shared_ptr<MoFEM::FEMethod> method,
                       boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                       boost::shared_ptr<MoFEM::BasicMethod> post_only);
/**
 * \brief set TS implicit function evaluation function
 * \ingroup dm
 */
PetscErrorCode DMMoFEMTSSetI2Function(DM dm, const char fe_name[],
                                      MoFEM::FEMethod *method,
                                      MoFEM::BasicMethod *pre_only,
                                      MoFEM::BasicMethod *post_only);

/**
 * \brief set TS Jacobian evaluation function
 * \ingroup dm
 */
PetscErrorCode
DMMoFEMTSSetI2Jacobian(DM dm, const std::string fe_name,
                       boost::shared_ptr<MoFEM::FEMethod> method,
                       boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                       boost::shared_ptr<MoFEM::BasicMethod> post_only);
/**
 * \brief set TS Jacobian evaluation function
 * \ingroup dm
 */
PetscErrorCode
DMMoFEMTSSetI2Jacobian(DM dm, const char fe_name[],
                                     MoFEM::FEMethod *method,
                                     MoFEM::BasicMethod *pre_only,
                                     MoFEM::BasicMethod *post_only);

/**
 * @brief Set Monitor To TS solver
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSMonitorSet.html>See
 * PETSc documentaton here</a>
 *
 * @param dm
 * @param ts time solver
 * @param fe_name
 * @param method
 * @param pre_only
 * @param post_only
 * @return PetscErrorCod
 */
PetscErrorCode
DMMoFEMTSSetMonitor(DM dm, TS ts, const std::string fe_name,
                    boost::shared_ptr<MoFEM::FEMethod> method,
                    boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                    boost::shared_ptr<MoFEM::BasicMethod> post_only);

/**
 * @brief Set Monitor To TS solver
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSMonitorSet.html>See
 * PETSc documentaton here</a>
 *
 * @param dm
 * @param ts time solver
 * @param fe_name
 * @param method
 * @param pre_only
 * @param post_only
 * @return PetscErrorCod
 */
PetscErrorCode DMMoFEMTSSetMonitor(DM dm, TS ts, const char fe_name[],
                                   MoFEM::FEMethod *method,
                                   MoFEM::BasicMethod *pre_only,
                                   MoFEM::BasicMethod *post_only);

/**
 * \brief get MoFEM::KspCtx data structure
 * \ingroup dm
 */
PetscErrorCode DMMoFEMGetKspCtx(DM dm, MoFEM::KspCtx **ksp_ctx);

/**
 * \brief get MoFEM::KspCtx data structure
 * \ingroup dm
 */
PetscErrorCode
DMMoFEMGetKspCtx(DM dm, const boost::shared_ptr<MoFEM::KspCtx> &ksp_ctx);

/**
 * \brief set MoFEM::KspCtx data structure
 * \ingroup dm
 */
PetscErrorCode DMMoFEMSetKspCtx(DM dm,
                                boost::shared_ptr<MoFEM::KspCtx> ksp_ctx);

/**
 * \brief get MoFEM::SnesCtx data structure
 * \ingroup dm
 */
PetscErrorCode DMMoFEMGetSnesCtx(DM dm, MoFEM::SnesCtx **snes_ctx);

/**
 * \brief get MoFEM::SnesCtx data structure
 * \ingroup dm
 */
PetscErrorCode
DMMoFEMGetSnesCtx(DM dm, const boost::shared_ptr<MoFEM::SnesCtx> &snes_ctx);

/**
  * \brief Set MoFEM::SnesCtx data structure
  * \ingroup dm

  */
PetscErrorCode DMMoFEMSetSnesCtx(DM dm,
                                 boost::shared_ptr<MoFEM::SnesCtx> snes_ctx);

/**
 * \brief get MoFEM::TsCtx data structure
 * \ingroup dm
 */
PetscErrorCode DMMoFEMGetTsCtx(DM dm, MoFEM::TsCtx **ts_ctx);

/**
 * \brief get MoFEM::TsCtx data structure
 * \ingroup dm
 */
PetscErrorCode DMMoFEMGetTsCtx(DM dm,
                               const boost::shared_ptr<MoFEM::TsCtx> &ts_ctx);

/**
  * \brief Set MoFEM::TsCtx data structure
  * \ingroup dm

  It take over pointer, do not delete it, DM will destroy pointer
  when is destroyed.

  */
PetscErrorCode DMMoFEMSetTsCtx(DM dm, boost::shared_ptr<MoFEM::TsCtx> ts_ctx);

/** sets if read mesh is partitioned
 * \ingroup dm
 */
PetscErrorCode DMMoFEMSetIsPartitioned(DM dm, PetscBool is_partitioned);

/** get if read mesh is partitioned
 * \ingroup dm
 */
PetscErrorCode DMMoFEMGetIsPartitioned(DM dm, PetscBool *is_partitioned);

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
PetscErrorCode DMCreateGlobalVector_MoFEM(DM dm, Vec *g);

/**
 * \brief DMShellSetCreateGlobalVector
 * \ingroup dm
 *
 * sets the routine to create a global vector
 * associated with the shell DM
 */
PetscErrorCode DMCreateGlobalVector_MoFEM(DM dm, SmartPetscObj<Vec> &g_ptr);

/**
 * \brief DMShellSetCreateLocalVector
 * \ingroup dm
 *
 * sets the routine to create a local vector
 * associated with the shell DM
 */
PetscErrorCode DMCreateLocalVector_MoFEM(DM dm, Vec *l);

/**
 * DMShellSetCreateMatrix
 * \ingroup dm
 *
 * sets the routine to create a matrix associated with the shell DM
 */
PetscErrorCode DMCreateMatrix_MoFEM(DM dm, Mat *M);

/**
 * DMShellSetCreateMatrix
 * \ingroup dm
 *
 * sets the routine to create a matrix associated with the shell DM
 */
PetscErrorCode DMCreateMatrix_MoFEM(DM dm, SmartPetscObj<Mat> &M);

/**
 * Set options for MoFEM DM
 * \ingroup dm
 */
#if PETSC_VERSION_GE(3, 7, 0)
PetscErrorCode DMSetFromOptions_MoFEM(PetscOptionItems *PetscOptionsObject,
                                      DM dm);
#elif PETSC_VERSION_GE(3, 5, 3)
PetscErrorCode DMSetFromOptions_MoFEM(PetscOptions *PetscOptionsObject, DM dm);
#else
PetscErrorCode DMSetFromOptions_MoFEM(DM dm);
#endif

/**
 * sets up the MoFEM structures inside a DM object
 * \ingroup dm
 */
PetscErrorCode DMSetUp_MoFEM(DM dm);

/**
 * Sets up the MoFEM structures inside a DM object for sub dm
 * \ingroup dm
 */
PetscErrorCode DMSubDMSetUp_MoFEM(DM subdm);

/**
 * Add field to sub dm problem on rows
 * \ingroup dm
 */
PetscErrorCode DMMoFEMAddSubFieldRow(DM dm, const char field_name[],
                                     EntityType lo_type = MBVERTEX,
                                     EntityType hi_type = MBMAXTYPE);

/**
 * Add field to sub dm problem on columns
 * \ingroup dm
 */
PetscErrorCode DMMoFEMAddSubFieldCol(DM dm, const char field_name[],
                                     EntityType lo_type = MBVERTEX,
                                     EntityType hi_type = MBMAXTYPE);

/**
 * Return true if this DM is sub problem
  \ingroup dm
 * @param  dm            the DM object
 * @param  is_subproblem true if subproblem
 * @return               error code
 */
PetscErrorCode DMMoFEMGetIsSubDM(DM dm, PetscBool *is_sub_dm);

/**
 * \brief get sub problem is
 * @param  dm has to be created with DMMoFEMSetSquareProblem
 * @param  is return is on the row
 * @return    error code
 *
 * Returns IS with global indices of the DM used to create SubDM
 *
 */
PetscErrorCode DMMoFEMGetSubRowIS(DM dm, IS *is);

/**
 * \brief get sub problem is
 * @param  dm has to be created with DMMoFEMSetSquareProblem
 * @param  is return is on the row
 * @return    error code
 *
 * Returns IS with global indices of the DM used to create SubDM
 *
 */
PetscErrorCode DMMoFEMGetSubColIS(DM dm, IS *is);

/**
 * \brief Add problem to composite DM on row
  \ingroup dm
 *
 * This create block on row with DOFs from problem of given name
 *
 * @param  dm       the DM object
 * @param  prb_name add problem name
 * @return          error code
 */
PetscErrorCode DMMoFEMAddRowCompositeProblem(DM dm, const char prb_name[]);

/**
 * \brief Add problem to composite DM on col
 * \ingroup dm
 *
 * This create block on col with DOFs from problem of given name
 *
 * @param  dm       the DM object
 * @param  prb_name add problem name
 * @return          error code
 */
PetscErrorCode DMMoFEMAddColCompositeProblem(DM dm, const char prb_name[]);

/**
 * \brief Get if this DM is composite DM
 * \ingroup dm
 *
 * @param  dm         the DM object
 * @param  is_comp_dm return true if composite problem here
 * @return            error code
 */
PetscErrorCode DMMoFEMGetIsCompDM(DM dm, PetscBool *is_comp_dm);

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
PetscErrorCode DMGlobalToLocalBegin_MoFEM(DM dm, Vec, InsertMode, Vec);

/**
 * DMShellSetGlobalToLocal
 * \ingroup dm
 *
 * the routine that begins the global to local scatter
 */
PetscErrorCode DMGlobalToLocalEnd_MoFEM(DM dm, Vec, InsertMode, Vec);

/**
 * DMShellSetLocalToGlobal
 * \ingroup dm
 *
 * the routine that begins the local to global scatter
 */
PetscErrorCode DMLocalToGlobalBegin_MoFEM(DM, Vec, InsertMode, Vec);

/**
 * DMShellSetLocalToGlobal
 * \ingroup dm
 *
 * the routine that ends the local to global scatter
 */
PetscErrorCode DMLocalToGlobalBegin_MoFEM(DM, Vec, InsertMode, Vec);

/**
 * DMShellSetLocalToGlobal
 * \ingroup dm
 *
 * the routine that ends the local to global scatter
 */
PetscErrorCode DMLocalToGlobalEnd_MoFEM(DM, Vec, InsertMode, Vec);

/**
 * Creates a set of IS objects with the global indices of dofs for each field
 * @param  dm         The number of fields (or NULL if not requested)
 *
 * Output:
 * @param  numFields  The number of fields (or NULL if not requested)
 * @param  fieldNames The name for each field (or NULL if not requested)
 * @param  fields     The global indices for each field (or NULL if not
 requested)
 *
 * @return            error code

 * \note The user is responsible for freeing all requested arrays. In
 particular,
 * every entry of names should be freed with PetscFree(), every entry of fields
 * should be destroyed with ISDestroy(), and both arrays should be freed with
 * PetscFree().

  \ingroup dm

 */
PetscErrorCode DMCreateFieldIS_MoFEM(DM dm, PetscInt *numFields,
                                     char ***fieldNames, IS **fields);

/**
 * \brief get field is in the problem
 * @param  dm         the DM objects
 * @param  rc         ROW or COL (no difference is problem is squared)
 * @param  field_name name of the field
 * @param  is         returned the IS object
 * @return            error code
 *
 * \code
 * IS is;
 * ierr = DMMoFEMGetFieldIS(dm,ROW,"DISP",&is_disp); CHKERRG(ierr);
 * \endcode
 *
 *
  \ingroup dm
 */
PetscErrorCode DMMoFEMGetFieldIS(DM dm, RowColData rc, const char field_name[],
                                 IS *is);

/**
 * @brief Set verbosity level
 *
 * @param dm
 * @param verb see VERBOSITY_LEVELS for list of the levels
 * @return PetscErrorCode
 */
PetscErrorCode DMMoFEMSetVerbosity(DM dm, const int verb);

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
 * \ingroup dm
 *
 */
struct DMCtx : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  Interface *mField_ptr;    ///< MoFEM interface
  PetscBool isProblemBuild; ///< True if problem is build
  std::string problemName;  ///< Problem name

  // Options
  PetscBool isPartitioned;  ///< true if read mesh is on parts
  PetscBool isSquareMatrix; ///< true if rows equals to cols

  int rAnk, sIze;

  // pointer to data structures
  const Problem *problemPtr; ///< pinter to problem data structure

  // sub problem
  PetscBool isSubDM;
  std::vector<std::string> rowFields;
  std::vector<std::string> colFields;
  const Problem *problemMainOfSubPtr; ///< pinter to main problem to sub-problem

  PetscBool isCompDM;
  std::vector<std::string> rowCompPrb;
  std::vector<std::string> colCompPrb;
  boost::shared_ptr<std::map<std::string, std::pair<EntityType, EntityType>>>
      mapTypeRow;
  boost::shared_ptr<std::map<std::string, std::pair<EntityType, EntityType>>>
      mapTypeCol;

  PetscBool destroyProblem; ///< If true destroy problem with DM

  DMCtx();
  virtual ~DMCtx() = default;

  int verbosity; ///< verbosity
  int referenceNumber;

  boost::shared_ptr<KspCtx> kspCtx;   ///< data structure KSP
  boost::shared_ptr<SnesCtx> snesCtx; ///< data structure SNES
  boost::shared_ptr<TsCtx> tsCtx;     ///< data structure for TS solver
};

/**
 * @brief get problem pointer from DM
 * 
 */
inline auto getProblemPtr(DM dm) {
  const MoFEM::Problem *problem_ptr;
  CHK_THROW_MESSAGE(DMMoFEMGetProblemPtr(dm, &problem_ptr),
                    "Get cot get problem ptr from DM");
  return problem_ptr;
};

/**
 * @brief Get smart matrix from DM
 * \ingroup dm
 * 
 */
inline auto smartCreateDMMatrix(DM dm) {
  SmartPetscObj<Mat> a;
  ierr = DMCreateMatrix_MoFEM(dm, a);
  CHKERRABORT(getCommFromPetscObject(reinterpret_cast<PetscObject>(dm)), ierr);
  return a;
};

/**
 * @brief Get smart vector from DM
 * \ingroup dm
 * 
 */
inline auto smartCreateDMVector(DM dm) {
  SmartPetscObj<Vec> v;
  ierr = DMCreateGlobalVector_MoFEM(dm, v);
  CHKERRABORT(getCommFromPetscObject(reinterpret_cast<PetscObject>(dm)), ierr);
  return v;
};

/**
 * @brief Get KSP context data structure used by DM
 * 
 */
inline auto smartGetDMKspCtx(DM dm) {
  boost::shared_ptr<MoFEM::KspCtx> ksp_ctx;
  ierr = DMMoFEMGetKspCtx(dm, ksp_ctx);
  CHKERRABORT(getCommFromPetscObject(reinterpret_cast<PetscObject>(dm)), ierr);
  return ksp_ctx;
};

/**
 * @brief Get SNES context data structure used by DM
 * 
 */
inline auto smartGetDMSnesCtx(DM dm) {
  boost::shared_ptr<MoFEM::SnesCtx> snes_ctx;
  ierr = DMMoFEMGetSnesCtx(dm, snes_ctx);
  CHKERRABORT(getCommFromPetscObject(reinterpret_cast<PetscObject>(dm)), ierr);
  return snes_ctx;
};

/**
 * @brief Get TS context data structure used by DM
 * 
 */
inline auto smartGetDMTsCtx(DM dm) {
  boost::shared_ptr<MoFEM::TsCtx> ts_ctx;
  ierr = DMMoFEMGetTsCtx(dm, ts_ctx);
  CHKERRABORT(getCommFromPetscObject(reinterpret_cast<PetscObject>(dm)), ierr);
  return ts_ctx;
};


} // namespace MoFEM

#endif //__DMMMOFEM_H

/**
 * \defgroup dm Distributed mesh manager
 * \brief Implementation of PETSc DM, managing interactions between mesh data
 *structures and vectors and matrices
 *
 * DM objects are used to manage communication between the algebraic structures
 *in PETSc (Vec and Mat) and mesh data structures in PDE-based (or other)
 * simulations.
 *
 * DM is abstract interface, here is it particular implementation for MoFEM
 *code.
 *
 * \ingroup mofem
 **/
