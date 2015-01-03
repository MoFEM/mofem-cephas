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

namespace MoFEM {

struct DMCtx {

  FieldInterface &mField; 	//< MoFEM interface
  string problemName; 		//< problen name

  /// constructor
  DMCtx(FieldInterface &m_field,string problem_name); 

  // destructor
  ~DMCtx(); 

  const FieldInterface& get_mField() const { return mField; }
  const Interface& get_moab() const { return mField.get_moab(); }

  PetscLogEvent USER_EVENT_CreateGlobalVector;
  PetscLogEvent USER_EVENT_CreateMatrix;

  friend PetscErrorCode DMCreate_MoFEM(DM dm);
  friend PetscErrorCode DMDestroym_MoFEM(DM dm);
  friend PetscErrorCode DMCreateGlobalVector_MoFEM(DM dm,Vec *globV);
  friend PetscErrorCode DMCreateLocalVector_MoFEM(DM dm,Vec *locV);
  friend PetscErrorCode DMCreateMatrix_MoFEM(DM dm,Mat *M);
  friend PetscErrorCode DMSetUp_MoFEM(DM dm); 
  friend PetscErrorCode DMSetFromOptions_MoFEM(DM dm);
  friend PetscErrorCode DMGlobalToLocalBegin_MoFEM(DM dm,Vec,InsertMode,Vec);
  friend PetscErrorCode DMGlobalToLocalEnd_MoFEM(DM dm,Vec,InsertMode,Vec);
  friend PetscErrorCode DMLocalToGlobalBegin_MoFEM(DM,Vec,InsertMode,Vec);
  friend PetscErrorCode DMLocalToGlobalEnd_MoFEM(DM,Vec,InsertMode,Vec);

};

PetscErrorCode DMCreate_MoFEM(DM dm);
PetscErrorCode DMDestroym_MoFEM(DM dm);

/** 
 * \brief DMShellSetCreateGlobalVector 
 *
 * sets the routine to create a global vector
 * associated with the shell DM
 */
PetscErrorCode DMCreateGlobalVector_MoFEM(DM dm,Vec *globV);

/** 
 * \brief DMShellSetCreateLocalVector 
 *
 * sets the routine to create a local vector
 * associated with the shell DM
 */
PetscErrorCode DMCreateLocalVector_MoFEM(DM dm,Vec *locV);

/**
  * DMShellSetCreateMatrix
  *
  * sets the routine to create a matrix associated with the shell DM
  */
PetscErrorCode DMCreateMatrix_MoFEM(DM dm,Mat *M);

/**
  * Set options for MoFEM DM
  */
PetscErrorCode DMSetFromOptions_MoFEM(DM dm);

/** 
  * sets up the MoFEM structures inside a DM object
  */
PetscErrorCode DMSetUp_MoFEM(DM dm);

/** 
  * destroy the MoFEM strutcture 
  */
PetscErrorCode DMDestroy_MoFEM(DM dm);

/**
  * DMShellSetGlobalToLocal
  *
  * the routine that begins the global to local scatter
  */
PetscErrorCode DMGlobalToLocalBegin_MoFEM(DM dm,Vec,InsertMode,Vec);

/**
  * DMShellSetGlobalToLocal
  *
  * the routine that begins the global to local scatter
  */
PetscErrorCode DMGlobalToLocalEnd_MoFEM(DM dm,Vec,InsertMode,Vec);

/**
  * DMShellSetLocalToGlobal
  *
  * the routine that begins the local to global scatter
  */
PetscErrorCode DMLocalToGlobalBegin_MoFEM(DM,Vec,InsertMode,Vec);

/**
  * DMShellSetLocalToGlobal
  *
  * the routine that ends the local to global scatter
  */
PetscErrorCode DMLocalToGlobalBegin_MoFEM(DM,Vec,InsertMode,Vec);

/**
  * DMShellSetLocalToGlobal
  *
  * the routine that ends the local to global scatter
  */
PetscErrorCode DMLocalToGlobalEnd_MoFEM(DM,Vec,InsertMode,Vec);

}

#endif //__DMMMOFEM_H

