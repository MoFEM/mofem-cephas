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

/** 
  * \brief Register MoFEM problem
  * \ingroup dm
  */
PetscErrorCode DMRegister_MoFEM(const char sname[]);

/** 
  * \brief Must be called by user to set MoFEM data structures
  * \ingroup dm
  */
PetscErrorCode DMMoFEMCreateMoFEM(MoFEM::FieldInterface *m_field_ptr,const char problem_name[],const MoFEM::BitRefLevel &bit_level,DM dm);

/** 
  * \brief add element to dm
  * \ingroup dm
  */
PetscErrorCode DMMoFEMAddElement(const char fe_name[],DM dm);

/** 
  * \brief set local (or ghosted) vector values on mesh for partition only
  * \ingroup dm
  */
PetscErrorCode DMoFEMMeshToLocalVector(Vec l,InsertMode mode,ScatterMode scatter_mode,DM dm);

/** 
  * \brief set ghosted vector values on all existing mesh entities
  * \ingroup dm
  */
PetscErrorCode DMoFEMMeshToGlobalVector(Vec g,InsertMode mode,ScatterMode scatter_mode,DM dm);

/** 
  * \brief execute finite element method for each element in dm (problem)
  * \ingroup dm
  */
PetscErrorCode DMoFEMPreProcessFiniteElements(MoFEM::FEMethod *method,DM dm);

/** 
  * \brief execute finite element method for each element in dm (problem)
  * \ingroup dm
  */
PetscErrorCode DMoFEMPostProcessFiniteElements(MoFEM::FEMethod *method,DM dm);

/** 
  * \brief execute finite element method for each element in dm (problem)
  * \ingroup dm
  */
PetscErrorCode DMoFEMLoopFiniteElements(const char fe_name[],MoFEM::FEMethod *method,DM dm);

/** 
  * \brief execute method for dofs on field in problem 
  * \ingroup dm
  */
PetscErrorCode DMoFEMLoopDofs(const char field_name[],MoFEM::EntMethod *method,DM dm);

/**
  * \brief set KSP right hand side evaluation function
  * \ingroup dm
  */
PetscErrorCode DMMoFEMKSPSetComputeRHS(const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only);

/**
  * \brief set KSP matrix evaluation function
  * \ingroup dm
  */
PetscErrorCode DMMoFEMKSPSetComputeOperators(const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only);




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
PetscErrorCode DMCreateGlobalVector_MoFEM(DM dm,Vec *globV);

/** 
 * \brief DMShellSetCreateLocalVector 
 * \ingroup dm
 *
 * sets the routine to create a local vector
 * associated with the shell DM
 */
PetscErrorCode DMCreateLocalVector_MoFEM(DM dm,Vec *locV);

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
PetscErrorCode DMSetFromOptions_MoFEM(DM dm);

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
 * \defgroup dm MoFem discreat manager
 * \ingroup mofem
 ******************************************************************************/


