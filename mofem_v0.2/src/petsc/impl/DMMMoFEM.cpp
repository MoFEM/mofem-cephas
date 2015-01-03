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

#include <petscsys.h>
#include <petscvec.h> 
#include <petscmat.h> 
#include <petscsnes.h> 
#include <petscts.h> 

#include <definitions.h>
#include <Common.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>
#include <DMMoFEM.hpp>

#include <petsc-private/dmimpl.h> /*I  "petscdm.h"   I*/
#include <petsc-private/vecimpl.h> /*I  "petscdm.h"   I*/

namespace MoFEM {

DMCtx::DMCtx(FieldInterface &m_field,string problem_name): 
  mField(m_field), problemName(problem_name) {
  PetscLogEventRegister("DMMoFEMCreateGlobalVector",0,&USER_EVENT_CreateGlobalVector);
  PetscLogEventRegister("DMMMoFEMCreateMatrix",0,&USER_EVENT_CreateMatrix);

  string sname = "MoFEM_"+problem_name;
  DMRegister(sname.c_str(),DMCreate_MoFEM);

}

PetscErrorCode DMCreate_MoFEM(DM dm) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

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
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"DM function not implenented into MoFEM");
  PetscFunctionReturn(0);
}

PetscErrorCode DMSetUp_MoFEM(DM dm) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"DM function not implenented into MoFEM");
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

}





