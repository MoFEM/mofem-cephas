/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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
  * \brief set squared problem
  * \ingroup dm

  It if trure is assumed that matrix has the same indexing on rows and
  collumns. This reduces interprocessor communication.

  */
PetscErrorCode DMMoFEMSetSquareProblem(DM dm,PetscBool square_problem);

/** 
  * \brief get squared problem
  * \ingroup dm

  It if trure is assumed that matrix has the same indexing on rows and
  collumns. This reduces interprocessor communication.

  */
PetscErrorCode DMMoFEMGetSquareProblem(DM dm,PetscBool *square_problem);

/** 
  * \brief add element to dm
  * \ingroup dm
  */
PetscErrorCode DMMoFEMAddElement(DM dm,const char fe_name[]);

/** 
  * \brief set local (or ghosted) vector values on mesh for partition only
  * \ingroup dm
  */
PetscErrorCode DMoFEMMeshToLocalVector(DM dm,Vec l,InsertMode mode,ScatterMode scatter_mode);

/** 
  * \brief set ghosted vector values on all existing mesh entities
  * \ingroup dm
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
  * \brief execute finite element method for each element in dm (problem)
  * \ingroup dm
  */
PetscErrorCode DMoFEMLoopFiniteElements(DM dm,const char fe_name[],MoFEM::FEMethod *method);

/** 
  * \brief execute method for dofs on field in problem 
  * \ingroup dm
  */
PetscErrorCode DMoFEMLoopDofs(DM dm,const char field_name[],MoFEM::EntMethod *method);

/**
  * \brief set KSP right hand side evaluation function
  * \ingroup dm
  */
PetscErrorCode DMMoFEMKSPSetComputeRHS(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only);

/**
  * \brief set KSP matrix evaluation function
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

#ifdef __MOABKSP_HPP__

/**
  * \brief get MoFEM::KspCtx data structure
  * \ingroup dm 
  */
PetscErrorCode DMMoFEMGetKspCtx(DM dm,MoFEM::KspCtx **ksp_ctx); 

#endif

#ifdef __MOABSNES_HPP__

/**
  * \brief get MoFEM::SnesCtx data structure
  * \ingroup dm 
  */
PetscErrorCode DMMoFEMGetSnesCtx(DM dm,MoFEM::SnesCtx **snes_ctx); 

#endif 

#ifdef __MOABTS_HPP__

/**
  * \brief get MoFEM::TsCtx data structure
  * \ingroup dm 
  */
PetscErrorCode DMMoFEMGetTsCtx(DM dm,MoFEM::TsCtx **ts_ctx); 

#endif

/** sets if read mesh is partitioned
  * \ingroup dm
  */
PetscErrorCode DMMoFEMSetIsPartitioned(DM dm,PetscBool is_partitioned);

/** get if read mesh is partitioned
  * \ingroup dm
  */
PetscErrorCode DMMoFEMGetIsPartitioned(DM dm,PetscBool *is_partitioned);

/** 
  * \brief Create dm data structure with MoFEM data structure
  * \ingroup dm
  */
PetscErrorCode DMCreate_MoFEM(DM dm);

/** 
  * \brief Destroys dm with MoFEM data structure
  * \ingroup dm
  */
PetscErrorCode DMDestroym_MoFEM(DM dm);

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
  * destroy the MoFEM strutcture 
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

#endif //__DMMMOFEM_H

/***************************************************************************//**
 * \defgroup dm MoFEM discreat manager
 ******************************************************************************/


