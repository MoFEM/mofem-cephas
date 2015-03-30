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

#include <KspCtx.hpp>
#include <SnesCtx.hpp>
#include <TsCtx.hpp>

#include <cblas.h>

namespace MoFEM {

struct DMCtx {

  FieldInterface *mField_ptr; 		//< MoFEM interface
  string problemName;			//< problem name
  
  KspCtx *kspCtx;			//< data structure KSP
  SnesCtx *snesCtx;			//< data structure SNES
  TsCtx	*tsCtx;				//< data structure for TS solver

  //options
  PetscBool isPartitioned;		//< true if read mesh is on parts
  PetscInt verbosity;			//< verbosity

  int rAnk,sIze;

  //global control
  static PetscBool isProblemsBuild;

  //pouinter to data structures
  const MoFEMProblem *problemPtr;	//< pinter to problem data struture

  DMCtx(); 
  virtual ~DMCtx();

  friend PetscErrorCode DMCreate_MoFEM(DM dm);
  friend PetscErrorCode DMDestroy_MoFEM(DM dm);
  friend PetscErrorCode DMCreateGlobalVector_MoFEM(DM dm,Vec *globV);
  friend PetscErrorCode DMCreateLocalVector_MoFEM(DM dm,Vec *locV);
  friend PetscErrorCode DMCreateMatrix_MoFEM(DM dm,Mat *M);
  friend PetscErrorCode DMSetUp_MoFEM(DM dm); 
  #if PETSC_VERSION_GE(3,5,3)
    friend PetscErrorCode DMSetFromOptions_MoFEM(PetscOptions *PetscOptionsObject,DM dm);
  #else 
    friend PetscErrorCode DMSetFromOptions_MoFEM(DM dm);
  #endif
  friend PetscErrorCode DMGlobalToLocalBegin_MoFEM(DM dm,Vec,InsertMode,Vec);
  friend PetscErrorCode DMGlobalToLocalEnd_MoFEM(DM dm,Vec,InsertMode,Vec);
  friend PetscErrorCode DMLocalToGlobalBegin_MoFEM(DM,Vec,InsertMode,Vec);
  friend PetscErrorCode DMLocalToGlobalEnd_MoFEM(DM,Vec,InsertMode,Vec);

};

PetscBool DMCtx::isProblemsBuild = PETSC_FALSE;

DMCtx::DMCtx(): 
  mField_ptr(PETSC_NULL),
  kspCtx(NULL),snesCtx(NULL),tsCtx(NULL),
  isPartitioned(PETSC_FALSE),
  verbosity(0) {
}
DMCtx::~DMCtx() {
  delete kspCtx;
  delete snesCtx;
  delete tsCtx;
}


}

using namespace MoFEM;

PetscErrorCode DMRegister_MoFEM(const char sname[]) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DMRegister(sname,DMCreate_MoFEM); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreate_MoFEM(DM dm) {
  //PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;

  dm->data = new DMCtx();

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

PetscErrorCode DMDestroy_MoFEM(DM dm) {
  //PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  delete (DMCtx*)dm->data;
  PetscFunctionReturn(0);
}

PetscErrorCode DMMoFEMCreateMoFEM(DM dm,MoFEM::FieldInterface *m_field_ptr,const char problem_name[],const MoFEM::BitRefLevel &bit_level) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(!dm->data) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"data structure for MoFEM not yet created");
  }
  if(!m_field_ptr) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"DM function not implemented into MoFEM");
  } 
  dm_field->mField_ptr = m_field_ptr;
  dm_field->problemName = problem_name;
  ierr = dm_field->mField_ptr->add_problem(dm_field->problemName,MF_ZERO); CHKERRQ(ierr);
  ierr = dm_field->mField_ptr->modify_problem_ref_level_add_bit(dm_field->problemName,bit_level); CHKERRQ(ierr);
  dm_field->kspCtx = new KspCtx(*m_field_ptr,problem_name);
  dm_field->snesCtx = new SnesCtx(*m_field_ptr,problem_name);
  dm_field->tsCtx = new TsCtx(*m_field_ptr,problem_name);

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dm,&comm); CHKERRQ(ierr);
  int result = 0;
  MPI_Comm_compare(comm,m_field_ptr->get_comm(),&result);
  //cerr << result << " " << MPI_IDENT << " " << MPI_CONGRUENT << " " << MPI_SIMILAR << " " << MPI_UNEQUAL << endl;
  if(result > MPI_CONGRUENT) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"MoFEM and DM using different communicators");
  }
  MPI_Comm_size(comm,&dm_field->sIze);
  MPI_Comm_rank(comm,&dm_field->rAnk);

  PetscFunctionReturn(0);
}

PetscErrorCode DMMoFEMAddElement(DM dm,const char fe_name[]) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->modify_problem_add_finite_element(dm_field->problemName,fe_name); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMoFEMMeshToLocalVector(DM dm,Vec l,InsertMode mode,ScatterMode scatter_mode) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->set_local_ghost_vector(dm_field->problemPtr,ROW,l,mode,scatter_mode); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMoFEMMeshToGlobalVector(DM dm,Vec g,InsertMode mode,ScatterMode scatter_mode) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->set_global_ghost_vector(dm_field->problemPtr,ROW,g,mode,scatter_mode); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMoFEMPreProcessFiniteElements(DM dm,MoFEM::FEMethod *method) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->problem_basic_method_preProcess(dm_field->problemPtr,*method); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMoFEMPostProcessFiniteElements(DM dm,MoFEM::FEMethod *method) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->problem_basic_method_postProcess(dm_field->problemPtr,*method); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMoFEMLoopFiniteElements(DM dm,const char fe_name[],MoFEM::FEMethod *method) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->loop_finite_elements(dm_field->problemPtr,fe_name,*method,dm_field->rAnk,dm_field->rAnk); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMoFEMLoopDofs(DM dm,const char field_name[],MoFEM::EntMethod *method) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->loop_dofs(dm_field->problemPtr,field_name,ROW,*method,dm_field->rAnk,dm_field->rAnk); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMMoFEMKSPSetComputeRHS(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(pre_only!=NULL) {
    dm_field->kspCtx->get_preProcess_to_do_Rhs().push_back(pre_only);
  }
  if(method!=NULL) {
    dm_field->kspCtx->get_loops_to_do_Rhs().push_back(KspCtx::loop_pair_type(fe_name,method));
  }
  if(post_only!=NULL) {
    dm_field->kspCtx->get_postProcess_to_do_Rhs().push_back(post_only);
  }
  ierr = DMKSPSetComputeRHS(dm,KspRhs,dm_field->kspCtx); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode DMMoFEMKSPSetComputeOperators(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(pre_only!=NULL) {
    dm_field->kspCtx->get_preProcess_to_do_Mat().push_back(pre_only);
  }
  if(method!=NULL) {
    dm_field->kspCtx->get_loops_to_do_Mat().push_back(KspCtx::loop_pair_type(fe_name,method));
  }
  if(post_only!=NULL) {
    dm_field->kspCtx->get_postProcess_to_do_Mat().push_back(post_only);
  }
  ierr = DMKSPSetComputeOperators(dm,KspMat,dm_field->kspCtx); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMMoFEMSNESSetFunction(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(pre_only!=NULL) {
    dm_field->snesCtx->get_preProcess_to_do_Rhs().push_back(pre_only);
  }
  if(method!=NULL) {
    dm_field->snesCtx->get_loops_to_do_Rhs().push_back(KspCtx::loop_pair_type(fe_name,method));
  }
  if(post_only!=NULL) {
    dm_field->snesCtx->get_postProcess_to_do_Rhs().push_back(post_only);
  }
  ierr = DMSNESSetFunction(dm,SnesRhs,dm_field->snesCtx); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMMoFEMSNESSetJacobian(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(pre_only!=NULL) {
    dm_field->snesCtx->get_preProcess_to_do_Rhs().push_back(pre_only);
  }
  if(method!=NULL) {
    dm_field->snesCtx->get_loops_to_do_Rhs().push_back(KspCtx::loop_pair_type(fe_name,method));
  }
  if(post_only!=NULL) {
    dm_field->snesCtx->get_postProcess_to_do_Rhs().push_back(post_only);
  }
  ierr = DMSNESSetJacobian(dm,SnesMat,dm_field->snesCtx); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMMoFEMTSSetIFunction(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(pre_only!=NULL) {
    dm_field->tsCtx->get_preProcess_to_do_IFunction().push_back(pre_only);
  }
  if(method!=NULL) {
    dm_field->tsCtx->get_loops_to_do_IFunction().push_back(KspCtx::loop_pair_type(fe_name,method));
  }
  if(post_only!=NULL) {
    dm_field->tsCtx->get_postProcess_to_do_IFunction().push_back(post_only);
  }
  ierr = DMTSSetIFunction(dm,f_TSSetIFunction,dm_field->tsCtx); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMMoFEMTSSetIJacobian(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(pre_only!=NULL) {
    dm_field->tsCtx->get_preProcess_to_do_IJacobian().push_back(pre_only);
  }
  if(method!=NULL) {
    dm_field->tsCtx->get_loops_to_do_IJacobian().push_back(KspCtx::loop_pair_type(fe_name,method));
  }
  if(post_only!=NULL) {
    dm_field->tsCtx->get_postProcess_to_do_IJacobian().push_back(post_only);
  }
  ierr = DMTSSetIJacobian(dm,f_TSSetIJacobian,dm_field->tsCtx); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMMoFEMGetKspCtx(DM dm,MoFEM::KspCtx **ksp_ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  *ksp_ctx = dm_field->kspCtx;
  PetscFunctionReturn(0);
}

PetscErrorCode DMMoFEMGetSnesCtx(DM dm,MoFEM::SnesCtx **snes_ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  *snes_ctx = dm_field->snesCtx;
  PetscFunctionReturn(0);
}

/** get if read mesh is partitioned
  * \ingroup dm
  */
PetscErrorCode DMMoFEMSetIsPartitioned(DM dm,PetscBool is_partitioned) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  dm_field->isPartitioned = is_partitioned;
  PetscFunctionReturn(0);
}

/** get if read mesh is partitioned
  * \ingroup dm
  */
PetscErrorCode DMMoFEMGetIsPartitioned(DM dm,PetscBool *is_partitioned) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  *is_partitioned = dm_field->isPartitioned;
  PetscFunctionReturn(0);
}

PetscErrorCode DMMoFEMGetTsCtx(DM dm,MoFEM::TsCtx **ts_ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  *ts_ctx = dm_field->tsCtx;
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateGlobalVector_MoFEM(DM dm,Vec *g) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->VecCreateGhost(dm_field->problemName,ROW,g); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateLocalVector_MoFEM(DM dm,Vec *l) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->VecCreateSeq(dm_field->problemName,ROW,l); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateMatrix_MoFEM(DM dm,Mat *M) {
  PetscErrorCode ierr;	
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->MatCreateMPIAIJWithArrays(dm_field->problemName,M); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

  #if PETSC_VERSION_GE(3,5,3)
    PetscErrorCode DMSetFromOptions_MoFEM(PetscOptions *PetscOptionsObject,DM dm) {
  #else 
    PetscErrorCode DMSetFromOptions_MoFEM(DM dm) {
  #endif 
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  #if PETSC_VERSION_GE(3,5,3)
    ierr = PetscOptionsHead(PetscOptionsObject,"DMMoFEM Options");CHKERRQ(ierr);
  #else 
    ierr = PetscOptionsHead("DMMoFEM Options");CHKERRQ(ierr);
  #endif
  ierr = PetscOptionsBool(
    "-dm_is_partitioned",
    "set if mesh is partitioned (works which native MOAB file formata, i.e. h5m","DMSetUp",
    dm_field->isPartitioned,&dm_field->isPartitioned,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMSetUp_MoFEM(DM dm) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(dm_field->isPartitioned) {
    ierr = dm_field->mField_ptr->build_partitioned_problems(); CHKERRQ(ierr);
    dm_field->isProblemsBuild = PETSC_TRUE;
    ierr = dm_field->mField_ptr->partition_finite_elements(dm_field->problemName,true,0,dm_field->sIze,1); CHKERRQ(ierr);
  } else {
    ierr = dm_field->mField_ptr->build_problems(); CHKERRQ(ierr);
    dm_field->isProblemsBuild = PETSC_TRUE;
    ierr = dm_field->mField_ptr->partition_problem(dm_field->problemName); CHKERRQ(ierr);
    ierr = dm_field->mField_ptr->partition_finite_elements(dm_field->problemName); CHKERRQ(ierr);
  }
  ierr = dm_field->mField_ptr->partition_ghost_dofs(dm_field->problemName); CHKERRQ(ierr);
  // dmmofem struture
  ierr = dm_field->mField_ptr->get_problem(dm_field->problemName,&dm_field->problemPtr); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMGlobalToLocalBegin_MoFEM(DM dm,Vec g,InsertMode mode,Vec l) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  
  ierr = VecGhostUpdateBegin(g,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DMGlobalToLocalEnd_MoFEM(DM dm,Vec g,InsertMode mode,Vec l) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;

  ierr = VecGhostUpdateEnd(g,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  DMCtx *dm_field = (DMCtx*)dm->data;
  int nb_dofs = dm_field->problemPtr->get_nb_local_dofs_row();
  int nb_ghost = dm_field->problemPtr->get_nb_ghost_dofs_row();

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
  ierr = VecRestoreArray(l,&array_loc); CHKERRQ(ierr);
  ierr = VecRestoreArray(g,&array_glob); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode DMLocalToGlobalBegin_MoFEM(DM dm,Vec l,InsertMode mode,Vec g) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  
  DMCtx *dm_field = (DMCtx*)dm->data;
  int nb_dofs = dm_field->problemPtr->get_nb_local_dofs_row();
  int nb_ghost = dm_field->problemPtr->get_nb_ghost_dofs_row();

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

  PetscFunctionReturn(0);
}

PetscErrorCode DMLocalToGlobalEnd_MoFEM(DM,Vec l,InsertMode mode,Vec g) {
  //PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}






