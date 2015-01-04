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

namespace MoFEM {

struct DMCtx {

  FieldInterface *mField_ptr; 	//< MoFEM interface
  string problemName;

  //options
  PetscBool is_partitioned;	//< true if read mesh is on parts

  //global control
  static PetscBool is_problems_build;

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

PetscBool DMCtx::is_problems_build = PETSC_FALSE;

DMCtx::DMCtx(): 
  mField_ptr(PETSC_NULL),
  is_partitioned(PETSC_FALSE)
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

  ierr = PetscNewLog(dm,(DMCtx**)&dm->data);CHKERRQ(ierr);

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
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"DM function not implenented into MoFEM");
  PetscFunctionReturn(0);
}

PetscErrorCode DMMoFEMSetCreateGlobalVector_MoFEM(DM dm,Vec *globV) {
  PetscFunctionBegin;
  DMCtx *dmmofem;
  dmmofem = (DMCtx*)(dm)->data;


  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"DM function not implenented into MoFEM");
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateLocalVector_MoFEM(DM dm,Vec *locV) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"DM function not implenented into MoFEM");
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateMatrix_MoFEM(DM dm,Mat *M) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"DM function not implenented into MoFEM");
  PetscFunctionReturn(0);
}

PetscErrorCode DMSetFromOptions_MoFEM(DM dm) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = PetscOptionsHead("DMMoFEM Options");CHKERRQ(ierr);
  ierr  = PetscOptionsBool("-dm_is_partitioned","set if mesh is partitioned (works which native MOAB file formata, i.e. h5m","DMSetUp",dm_field->is_partitioned,&dm_field->is_partitioned,NULL);CHKERRQ(ierr);
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
  if(dm_field->is_partitioned) {
    if(!dm_field->is_partitioned) {
      ierr = m_field_ptr->build_partitioned_problems(PETSC_COMM_WORLD,1); CHKERRQ(ierr);
      dm_field->is_partitioned = PETSC_TRUE;
    }
    ierr = m_field_ptr->partition_finite_elements(problem_name,true,0,pcomm->size(),1); CHKERRQ(ierr);
  } else {
    if(!dm_field->is_partitioned) {
      ierr = m_field_ptr->build_problems(); CHKERRQ(ierr);
      dm_field->is_partitioned = PETSC_TRUE;
    }
    ierr = m_field_ptr->partition_problem(problem_name); CHKERRQ(ierr);
    ierr = m_field_ptr->partition_finite_elements(problem_name); CHKERRQ(ierr);
  }
  ierr = m_field_ptr->partition_ghost_dofs(problem_name); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMGlobalToLocalBegin_MoFEM(DM dm,Vec,InsertMode,Vec) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"DM function not implenented into MoFEM");
  PetscFunctionReturn(0);
}

PetscErrorCode DMGlobalToLocalEnd_MoFEM(DM dm,Vec,InsertMode,Vec) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"DM function not implenented into MoFEM");
  PetscFunctionReturn(0);
}

PetscErrorCode DMMoFEMSetGlobalToLocalVecScatter(DM dm, VecScatter gtol) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"DM function not implenented into MoFEM");
  PetscFunctionReturn(0);
}

PetscErrorCode DMLocalToGlobalBegin_MoFEM(DM,Vec,InsertMode,Vec) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"DM function not implenented into MoFEM");
  PetscFunctionReturn(0);
}

PetscErrorCode DMLocalToGlobalEnd_MoFEM(DM,Vec,InsertMode,Vec) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"DM function not implenented into MoFEM");
  PetscFunctionReturn(0);
}






