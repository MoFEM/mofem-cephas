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

//petsc
#include <petscsys.h>
#include <petscvec.h> 
#include <petscmat.h> 
#include <petscsnes.h> 
#include <petscts.h> 

#include <petsc-private/dmimpl.h> /*I  "petscdm.h"   I*/
#include <petsc-private/vecimpl.h> /*I  "petscdm.h"   I*/

//moab
#include <moab/ParallelComm.hpp>

//mofem
#include <definitions.h>
#include <Common.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>

#include <DMMoFEM.hpp>

#include <cblas.h>

namespace MoFEM {

struct DMCtx {

  FieldInterface *mField_ptr; 	//< MoFEM interface
  string problemName;

  //options
  PetscBool isPartitioned;	//< true if read mesh is on parts
  PetscInt verbosity;		//< verbosity

  //global control
  static PetscBool isProblemsBuild;

  DMCtx(); 
  // destructor
  ~DMCtx(); 

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

PetscBool DMCtx::isProblemsBuild = PETSC_FALSE;

DMCtx::DMCtx(): 
  mField_ptr(PETSC_NULL),
  isPartitioned(PETSC_FALSE),
  verbosity(0)
  {}
}

using namespace MoFEM;

PetscErrorCode DMRegister_MoFEM(const char sname[],DM *dmb) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DMRegister(sname,DMCreate_MoFEM); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMMoFEMCreateMoFEM(FieldInterface *m_field_ptr,const char problem_name[],DM dmb) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)(dmb->data);
  dm_field->mField_ptr = m_field_ptr;
  dm_field->problemName = string(problem_name);
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreate_MoFEM(DM dm) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;

  ierr = PetscNewLog(dm,(DMCtx**)&dm->data); CHKERRQ(ierr);

  dm->ops->createglobalvector       = DMCreateGlobalVector_MoFEM;
  dm->ops->createlocalvector        = DMCreateLocalVector_MoFEM;
  dm->ops->creatematrix             = DMCreateMatrix_MoFEM;
  dm->ops->setup                    = DMSetUp_MoFEM;
  dm->ops->destroy                  = DMDestroy_MoFEM;
  dm->ops->setfromoptions           = DMSetFromOptions_MoFEM;
  dm->ops->globaltolocalbegin       = DMGlobalToLocalBegin_MoFEM;
  dm->ops->globaltolocalend         = DMGlobalToLocalEnd_MoFEM;
  dm->ops->localtoglobalbegin       = DMLocalToGlobalBegin_MoFEM;
  dm->ops->localtoglobalend         = DMLocalToGlobalEnd_MoFEM;

  PetscFunctionReturn(0);
}

PetscErrorCode DMDestroym_MoFEM(DM dm) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateGlobalVector_MoFEM(DM dm,Vec *globV) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  FieldInterface *m_field_ptr = dm_field->mField_ptr;
  string &problem_name = dm_field->problemName;
  ierr = m_field_ptr->VecCreateGhost(problem_name,ROW,globV); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateLocalVector_MoFEM(DM dm,Vec *locV) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  FieldInterface *m_field_ptr = dm_field->mField_ptr;
  string &problem_name = dm_field->problemName;
  ierr = m_field_ptr->VecCreateSeq(problem_name,ROW,locV); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateMatrix_MoFEM(DM dm,Mat *M) {
  PetscErrorCode ierr;	
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  FieldInterface *m_field_ptr = dm_field->mField_ptr;
  string &problem_name = dm_field->problemName;
  ierr = m_field_ptr->MatCreateMPIAIJWithArrays(problem_name,M); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMSetFromOptions_MoFEM(DM dm) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = PetscOptionsHead("DMMoFEM Options");CHKERRQ(ierr);
  ierr = PetscOptionsBool("-dm_is_partitioned","set if mesh is partitioned (works which native MOAB file formata, i.e. h5m","DMSetUp",dm_field->isPartitioned,&dm_field->isPartitioned,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMSetUp_MoFEM(DM dm) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  FieldInterface *m_field_ptr = dm_field->mField_ptr;
  string &problem_name = dm_field->problemName;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field_ptr->get_moab(),MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&m_field_ptr->get_moab(),PETSC_COMM_WORLD);
  if(dm_field->isPartitioned) {
    if(!dm_field->isPartitioned) {
      ierr = m_field_ptr->build_partitioned_problems(PETSC_COMM_WORLD,1); CHKERRQ(ierr);
      dm_field->isProblemsBuild = PETSC_TRUE;
    }
    ierr = m_field_ptr->partition_finite_elements(problem_name,true,0,pcomm->size(),1); CHKERRQ(ierr);
  } else {
    if(!dm_field->isPartitioned) {
      ierr = m_field_ptr->build_problems(); CHKERRQ(ierr);
      dm_field->isProblemsBuild = PETSC_TRUE;
    }
    ierr = m_field_ptr->partition_problem(problem_name); CHKERRQ(ierr);
    ierr = m_field_ptr->partition_finite_elements(problem_name); CHKERRQ(ierr);
  }
  ierr = m_field_ptr->partition_ghost_dofs(problem_name); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMGlobalToLocalBegin_MoFEM(DM dm,Vec g,InsertMode mode,Vec l) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;

  DMCtx *dm_field = (DMCtx*)dm->data;
  FieldInterface *m_field_ptr = dm_field->mField_ptr;
  string &problem_name = dm_field->problemName;
  const MoFEMProblem *problem_ptr;
  ierr = m_field_ptr->get_problem(problem_name,&problem_ptr); CHKERRQ(ierr);
  int nb_dofs = problem_ptr->get_nb_local_dofs_row();
  int nb_ghost = problem_ptr->get_nb_ghost_dofs_row();

  double *array_loc,*array_glob;
  ierr = VecGetArray(l,&array_loc); CHKERRQ(ierr);
  ierr = VecGetArray(g,&array_glob); CHKERRQ(ierr);
  switch (mode) {
    case INSERT_VALUES:
      cblas_dcopy(nb_dofs+nb_ghost,array_glob,1,array_loc,1);
    break;
    case ADD_VALUES:
      cblas_daxpy(nb_dofs+nb_ghost,1,array_glob,1,array_loc,1);
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  cblas_dcopy(nb_ghost,&array_glob[nb_dofs],1,&array_loc[nb_dofs],1);
  ierr = VecRestoreArray(l,&array_loc); CHKERRQ(ierr);
  ierr = VecRestoreArray(g,&array_glob); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMGlobalToLocalEnd_MoFEM(DM dm,Vec g,InsertMode mode,Vec l) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}


PetscErrorCode DMLocalToGlobalBegin_MoFEM(DM dm,Vec l,InsertMode mode,Vec g) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  
  DMCtx *dm_field = (DMCtx*)dm->data;
  FieldInterface *m_field_ptr = dm_field->mField_ptr;
  string &problem_name = dm_field->problemName;
  const MoFEMProblem *problem_ptr;
  ierr = m_field_ptr->get_problem(problem_name,&problem_ptr); CHKERRQ(ierr);
  int nb_dofs = problem_ptr->get_nb_local_dofs_row();
  int nb_ghost = problem_ptr->get_nb_local_dofs_row();

  double *array_loc,*array_glob;
  ierr = VecGetArray(l,&array_loc); CHKERRQ(ierr);
  ierr = VecGetArray(g,&array_glob); CHKERRQ(ierr);
  switch (mode) {
    case INSERT_VALUES:
      cblas_dcopy(nb_dofs+nb_ghost,array_loc,1,array_glob,1);
    break;
    case ADD_VALUES:
      cblas_daxpy(nb_dofs+nb_ghost,1,array_loc,1,array_glob,1);
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  ierr = VecRestoreArray(l,&array_loc); CHKERRQ(ierr);
  ierr = VecRestoreArray(g,&array_glob); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(g,mode,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(g,mode,SCATTER_REVERSE); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DMLocalToGlobalEnd_MoFEM(DM,Vec l,InsertMode mode,Vec g) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecGhostUpdateBegin(g,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(g,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}






