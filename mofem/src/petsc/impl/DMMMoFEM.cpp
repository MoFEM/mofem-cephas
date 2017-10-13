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

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <ProblemsManager.hpp>
#include <ISManager.hpp>
#include <VecManager.hpp>
#include <Core.hpp>

#include <AuxPETSc.hpp>
#include <KspCtx.hpp>
#include <SnesCtx.hpp>
#include <TsCtx.hpp>

// #undef PETSC_VERSION_RELEASE
// #define PETSC_VERSION_RELEASE 1

#if PETSC_VERSION_GE(3,6,0)
  #include <petsc/private/dmimpl.h> /*I  "petscdm.h"   I*/
  // #include <petsc/private/vecimpl.h> /*I  "petscdm.h"   I*/
#else
  #include <petsc-private/dmimpl.h> /*I  "petscdm.h"   I*/
  #include <petsc-private/vecimpl.h> /*I  "petscdm.h"   I*/
#endif

#include <DMMoFEM.hpp>

using namespace MoFEM;

DMCtx::DMCtx():
mField_ptr(PETSC_NULL),
isProblemBuild(PETSC_FALSE),
isPartitioned(PETSC_FALSE),
isSquareMatrix(PETSC_TRUE),
isSubDM(PETSC_FALSE),
isCompDM(PETSC_FALSE),
destroyProblem(PETSC_FALSE),
verbosity(0),
referenceNumber(0) {}

DMCtx::~DMCtx() {
  // cerr << "Snes " << snesCtx.use_count() << endl;
  // cerr << "Destroy DMCtx" << endl;
}

PetscErrorCode DMCtx::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if(uuid == IDD_DMCTX) {
    *iface = dynamic_cast<DMCtx*>(this);
    MoFEMFunctionReturnHot(0);
  }
  if(uuid == IDD_MOFEMUnknown) {
    *iface = dynamic_cast<UnknownInterface*>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMRegister_MoFEM(const char sname[]) {
  MoFEMFunctionBeginHot;
  ierr = DMRegister(sname,DMCreate_MoFEM); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMSetOperators_MoFEM(DM dm) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;

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
  dm->ops->createfieldis            = DMCreateFieldIS_MoFEM;

  // Default matrix type
  ierr = DMSetMatType(dm,MATMPIAIJ); CHKERRQ(ierr);

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMCreate_MoFEM(DM dm) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  dm->data = new DMCtx();
  ierr = DMSetOperators_MoFEM(dm); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMDestroy_MoFEM(DM dm) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  DMCtx *dm_field = (DMCtx*)dm->data;
  MoFEMFunctionBeginHot;
  if(((DMCtx*)dm->data)->referenceNumber==0) {
    if(dm_field->destroyProblem) {
      if(dm_field->mField_ptr->check_problem(dm_field->problemName)) {
        dm_field->mField_ptr->delete_problem(dm_field->problemName);
      } // else problem has to be deleted by the user
    }
    // cerr << "Destroy " << dm_field->problemName << endl;
    delete ((DMCtx*)dm->data);
  } else {
    // cerr << "Derefrence " << dm_field->problemName << " "  << ((DMCtx*)dm->data)->referenceNumber << endl;
    (((DMCtx*)dm->data)->referenceNumber)--;
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMCreateMoFEM(
  DM dm,
  MoFEM::Interface *m_field_ptr,
  const char problem_name[],
  const MoFEM::BitRefLevel bit_level,
  const MoFEM::BitRefLevel bit_mask
) {
  MoFEMFunctionBeginHot;

  DMCtx *dm_field = (DMCtx*)dm->data;
  if(!dm->data) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"data structure for MoFEM not yet created");
  }
  if(!m_field_ptr) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"DM function not implemented into MoFEM");
  }
  dm_field->mField_ptr = m_field_ptr;
  dm_field->problemName = problem_name;
  if(!m_field_ptr->check_problem(dm_field->problemName)) {
    // problem is not defined, declare problem internally set bool to destroyProblem
    // problem with DM
    dm_field->destroyProblem = PETSC_TRUE;
    ierr = dm_field->mField_ptr->add_problem(dm_field->problemName); CHKERRQ(ierr);
  } else {
    dm_field->destroyProblem = PETSC_FALSE;
  }
  ierr = dm_field->mField_ptr->modify_problem_ref_level_add_bit(
    dm_field->problemName,bit_level
  ); CHKERRQ(ierr);
  ierr = dm_field->mField_ptr->modify_problem_mask_ref_level_set_bit(
    dm_field->problemName,bit_mask
  ); CHKERRQ(ierr);
  dm_field->kspCtx = boost::shared_ptr<KspCtx>(new KspCtx(*m_field_ptr,problem_name));
  dm_field->snesCtx = boost::shared_ptr<SnesCtx>(new SnesCtx(*m_field_ptr,problem_name));
  dm_field->tsCtx = boost::shared_ptr<TsCtx>(new TsCtx(*m_field_ptr,problem_name));

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dm,&comm); CHKERRQ(ierr);
  int result = 0;
  MPI_Comm_compare(comm,m_field_ptr->get_comm(),&result);
  //std::cerr << result << " " << MPI_IDENT << " " << MPI_CONGRUENT << " " << MPI_SIMILAR << " " << MPI_UNEQUAL << std::endl;
  if(result > MPI_CONGRUENT) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"MoFEM and DM using different communicators");
  }
  MPI_Comm_size(comm,&dm_field->sIze);
  MPI_Comm_rank(comm,&dm_field->rAnk);

  // problem structure
  ierr = dm_field->mField_ptr->get_problem(dm_field->problemName,&dm_field->problemPtr); CHKERRQ(ierr);

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMCreateSubDM(DM subdm,DM dm,const char problem_name[]) {

  MoFEMFunctionBeginHot;

  DMCtx *dm_field = (DMCtx*)dm->data;
  if(!dm->data) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"data structure for MoFEM not yet created");
  }
  ierr = DMMoFEMCreateMoFEM(
    subdm,dm_field->mField_ptr,problem_name,dm_field->problemPtr->getBitRefLevel()
  ); CHKERRQ(ierr);

  DMCtx *subdm_field = (DMCtx*)subdm->data;
  subdm_field->isSubDM = PETSC_TRUE;
  subdm_field->problemMainOfSubPtr = dm_field->problemPtr;
  subdm_field->isPartitioned = dm_field->isPartitioned;
  subdm_field->isSquareMatrix = PETSC_FALSE;
  subdm->ops->setup = DMSubDMSetUp_MoFEM;

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMAddSubFieldRow(DM dm,const char field_name[]) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(!dm->data) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"data structure for MoFEM not yet created");
  }
  if(!dm_field->isSubDM) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"this is not sub-dm");
  }
  dm_field->rowFields.push_back(field_name);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMAddSubFieldCol(DM dm,const char field_name[]) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(!dm->data) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"data structure for MoFEM not yet created");
  }
  if(!dm_field->isSubDM) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"this is not sub-dm");
  }
  dm_field->colFields.push_back(field_name);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetIsSubDM(DM dm,PetscBool *is_sub_dm) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  *is_sub_dm = dm_field->isSubDM;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetSubRowIS(DM dm,IS *is) {

  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(dm_field->isSubDM!=PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"This DM is not created as a SubDM");
  }
  if(dm_field->isProblemBuild!=PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Problem is not build");
  }
  boost::shared_ptr<Problem::SubProblemData> sub_data = dm_field->problemPtr->getSubData();
  ierr = sub_data->getRowIs(is); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetSubColIS(DM dm,IS *is) {

  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(dm_field->isSubDM!=PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"This DM is not created as a SubDM");
  }
  if(dm_field->isProblemBuild!=PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Problem is not build");
  }
  boost::shared_ptr<Problem::SubProblemData> sub_data = dm_field->problemPtr->getSubData();
  ierr = sub_data->getColIs(is); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMAddRowCompositeProblem(DM dm,const char prb_name[]) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(!dm->data) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"data structure for MoFEM not yet created");
  }
  if(!dm_field->isCompDM) {
    dm_field->isCompDM = PETSC_TRUE;
  }
  dm_field->rowCompPrb.push_back(prb_name);
  if(dm_field->isSquareMatrix) {
    dm_field->colCompPrb.push_back(prb_name);
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMAddColCompositeProblem(DM dm,const char prb_name[]) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(!dm->data) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"data structure for MoFEM not yet created");
  }
  if(!dm_field->isCompDM) {
    dm_field->isCompDM = PETSC_TRUE;
  }
  if(dm_field->isSquareMatrix) {
    SETERRQ(
      PETSC_COMM_SELF,
      MOFEM_INVALID_DATA,
      "No need to add problem on column when problem block structurally symmetric"
    );
  }
  dm_field->colCompPrb.push_back(prb_name);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetIsCompDM(DM dm,PetscBool *is_comp_dm) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  *is_comp_dm = dm_field->isCompDM;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMGetInterfacePtr(DM dm,const MoFEM::Interface **m_field_ptr) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(!dm->data) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"data structure for MoFEM not yet created");
  }
  *m_field_ptr = dm_field->mField_ptr;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetProblemPtr(DM dm,const MoFEM::Problem **problem_ptr) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(!dm->data) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"data structure for MoFEM not yet created");
  }
  *problem_ptr = dm_field->problemPtr;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMSetDestroyProblem(DM dm,PetscBool destroy_problem)  {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  dm_field->destroyProblem = destroy_problem;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetDestroyProblem(DM dm,PetscBool *destroy_problem)  {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  *destroy_problem = dm_field->destroyProblem;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMSetSquareProblem(DM dm,PetscBool square_problem) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  dm_field->isSquareMatrix = square_problem;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMResolveSharedEntities(DM dm,const char fe_name[]) {

  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->resolve_shared_ents(dm_field->problemPtr,fe_name); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetProblemFiniteElementLayout(DM dm,const char fe_name[],PetscLayout *layout) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dm,&comm); CHKERRQ(ierr);
  ierr = dm_field->problemPtr->getNumberOfElementsByNameAndPart(comm,fe_name,layout); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetSquareProblem(DM dm,PetscBool *square_problem) {
  MoFEMFunctionBeginHot;
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  *square_problem = dm_field->isSquareMatrix;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMAddElement(DM dm,const char fe_name[]) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->modify_problem_add_finite_element(dm_field->problemName,fe_name); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMUnSetElement(DM dm,const char fe_name[]) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->modify_problem_unset_finite_element(dm_field->problemName,fe_name); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMMeshToLocalVector(DM dm,Vec l,InsertMode mode,ScatterMode scatter_mode) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->query_interface<VecManager>()->setLocalGhostVector(
    dm_field->problemPtr,COL,l,mode,scatter_mode
  ); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMMeshToGlobalVector(DM dm,Vec g,InsertMode mode,ScatterMode scatter_mode) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->query_interface<VecManager>()->setGlobalGhostVector(
    dm_field->problemPtr,COL,g,mode,scatter_mode
  ); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMPreProcessFiniteElements(DM dm,MoFEM::FEMethod *method) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->problem_basic_method_preProcess(dm_field->problemPtr,*method); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMPostProcessFiniteElements(DM dm,MoFEM::FEMethod *method) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->problem_basic_method_postProcess(dm_field->problemPtr,*method); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMLoopFiniteElementsUpAndLowRank(DM dm,const char fe_name[],MoFEM::FEMethod *method,int low_rank,int up_rank) {
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->loop_finite_elements(dm_field->problemPtr,fe_name,*method,low_rank,up_rank); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMLoopFiniteElementsUpAndLowRank(
  DM dm,const std::string& fe_name,boost::shared_ptr<MoFEM::FEMethod> method,int low_rank,int up_rank
) {
  return DMoFEMLoopFiniteElementsUpAndLowRank(dm,fe_name.c_str(),method.get(),low_rank,up_rank);
}

PetscErrorCode DMoFEMLoopFiniteElements(DM dm,const char fe_name[],MoFEM::FEMethod *method) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = DMoFEMLoopFiniteElementsUpAndLowRank(dm,fe_name,method,dm_field->rAnk,dm_field->rAnk); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMLoopFiniteElements(DM dm,const std::string& fe_name,boost::shared_ptr<MoFEM::FEMethod> method) {
  return DMoFEMLoopFiniteElements(dm,fe_name.c_str(),method.get());
}

PetscErrorCode DMoFEMLoopDofs(DM dm,const char field_name[],MoFEM::EntMethod *method) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->loop_dofs(dm_field->problemPtr,field_name,COL,*method,dm_field->rAnk,dm_field->rAnk); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

template<class S,class T>
static PetscErrorCode DMMoFEMKSPSetComputeRHS(DM dm,S fe_name,T method,T pre_only,T post_only) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(pre_only) {
    dm_field->kspCtx->get_preProcess_to_do_Rhs().push_back(pre_only);
  }
  if(method) {
    dm_field->kspCtx->get_loops_to_do_Rhs().push_back(PairNameFEMethodPtr(fe_name,method));
  }
  if(post_only) {
    dm_field->kspCtx->get_postProcess_to_do_Rhs().push_back(post_only);
  }
  ierr = DMKSPSetComputeRHS(dm,KspRhs,dm_field->kspCtx.get()); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMKSPSetComputeRHS(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only) {
  return DMMoFEMKSPSetComputeRHS<const char *,MoFEM::FEMethod*>(dm,fe_name,method,pre_only,post_only);
}

PetscErrorCode DMMoFEMKSPSetComputeRHS(
  DM dm,
  const std::string& fe_name,
  boost::shared_ptr<MoFEM::FEMethod> method,
  boost::shared_ptr<MoFEM::FEMethod> pre_only,
  boost::shared_ptr<MoFEM::FEMethod> post_only
) {
  return DMMoFEMKSPSetComputeRHS<const std::string&,boost::shared_ptr<MoFEM::FEMethod> >(
    dm,fe_name,method,pre_only,post_only
  );
}

template<class S,class T>
static PetscErrorCode DMMoFEMKSPSetComputeOperators(DM dm,S fe_name,T method,T pre_only,T post_only) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(pre_only) {
    dm_field->kspCtx->get_preProcess_to_do_Mat().push_back(pre_only);
  }
  if(method) {
    dm_field->kspCtx->get_loops_to_do_Mat().push_back(PairNameFEMethodPtr(fe_name,method));
  }
  if(post_only) {
    dm_field->kspCtx->get_postProcess_to_do_Mat().push_back(post_only);
  }
  ierr = DMKSPSetComputeOperators(dm,KspMat,dm_field->kspCtx.get()); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMKSPSetComputeOperators(
  DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only
) {
  return DMMoFEMKSPSetComputeOperators<const char *,MoFEM::FEMethod*>(dm,fe_name,method,pre_only,post_only);
}

PetscErrorCode DMMoFEMKSPSetComputeOperators(
  DM dm,
  const std::string& fe_name,
  boost::shared_ptr<MoFEM::FEMethod> method,
  boost::shared_ptr<MoFEM::FEMethod> pre_only,
  boost::shared_ptr<MoFEM::FEMethod> post_only
) {
  return DMMoFEMKSPSetComputeOperators<const std::string&,boost::shared_ptr<MoFEM::FEMethod> >(
    dm,fe_name,method,pre_only,post_only
  );
}

template<class S,class T>
static PetscErrorCode DMMoFEMSNESSetFunction(DM dm,S fe_name,T method,T pre_only,T post_only) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(pre_only) {
    dm_field->snesCtx->get_preProcess_to_do_Rhs().push_back(pre_only);
  }
  if(method) {
    dm_field->snesCtx->get_loops_to_do_Rhs().push_back(PairNameFEMethodPtr(fe_name,method));
  }
  if(post_only) {
    dm_field->snesCtx->get_postProcess_to_do_Rhs().push_back(post_only);
  }
  ierr = DMSNESSetFunction(dm,SnesRhs,dm_field->snesCtx.get()); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMSNESSetFunction(
  DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only
) {
  return DMMoFEMSNESSetFunction<const char *,MoFEM::FEMethod*>(dm,fe_name,method,pre_only,post_only);
}

PetscErrorCode DMMoFEMSNESSetFunction(
  DM dm,const std::string& fe_name,
  boost::shared_ptr<MoFEM::FEMethod> method,
  boost::shared_ptr<MoFEM::FEMethod> pre_only,
  boost::shared_ptr<MoFEM::FEMethod> post_only
) {
  return DMMoFEMSNESSetFunction<const std::string&,boost::shared_ptr<MoFEM::FEMethod> >(
    dm,fe_name,method,pre_only,post_only
  );
}

template<class S,class T>
static PetscErrorCode DMMoFEMSNESSetJacobian(DM dm,S fe_name,T method,T pre_only,T post_only) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(pre_only) {
    dm_field->snesCtx->get_preProcess_to_do_Mat().push_back(pre_only);
  }
  if(method) {
    dm_field->snesCtx->get_loops_to_do_Mat().push_back(PairNameFEMethodPtr(fe_name,method));
  }
  if(post_only) {
    dm_field->snesCtx->get_postProcess_to_do_Mat().push_back(post_only);
  }
  ierr = DMSNESSetJacobian(dm,SnesMat,dm_field->snesCtx.get()); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMSNESSetJacobian(
  DM dm,const char fe_name[],
  MoFEM::FEMethod *method,
  MoFEM::FEMethod *pre_only,
  MoFEM::FEMethod *post_only
) {
  return DMMoFEMSNESSetJacobian<const char *,MoFEM::FEMethod*>(dm,fe_name,method,pre_only,post_only);
}

PetscErrorCode DMMoFEMSNESSetJacobian(
  DM dm,const std::string& fe_name,
  boost::shared_ptr<MoFEM::FEMethod> method,
  boost::shared_ptr<MoFEM::FEMethod> pre_only,
  boost::shared_ptr<MoFEM::FEMethod> post_only
) {
  return DMMoFEMSNESSetJacobian<const std::string&,boost::shared_ptr<MoFEM::FEMethod> >(
    dm,fe_name,method,pre_only,post_only
  );
}

template<class S,class T>
static PetscErrorCode DMMoFEMTSSetIFunction(DM dm,S fe_name,T method,T pre_only,T post_only) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(pre_only) {
    dm_field->tsCtx->get_preProcess_to_do_IFunction().push_back(pre_only);
  }
  if(method) {
    dm_field->tsCtx->get_loops_to_do_IFunction().push_back(PairNameFEMethodPtr(fe_name,method));
  }
  if(post_only) {
    dm_field->tsCtx->get_postProcess_to_do_IFunction().push_back(post_only);
  }
  ierr = DMTSSetIFunction(dm,f_TSSetIFunction,dm_field->tsCtx.get()); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMTSSetIFunction(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only) {
  return DMMoFEMTSSetIFunction<const char *,MoFEM::FEMethod*>(dm,fe_name,method,pre_only,post_only);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMTSSetIFunction(
  DM dm,const std::string& fe_name,
  boost::shared_ptr<MoFEM::FEMethod> method,
  boost::shared_ptr<MoFEM::FEMethod> pre_only,
  boost::shared_ptr<MoFEM::FEMethod> post_only
) {
  return DMMoFEMTSSetIFunction<const std::string,boost::shared_ptr<MoFEM::FEMethod> >(
    dm,fe_name,method,pre_only,post_only
  );
  MoFEMFunctionReturnHot(0);
}

template<class S,class T>
static PetscErrorCode DMMoFEMTSSetIJacobian(DM dm,S fe_name,T method,T pre_only,T post_only) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(pre_only) {
    dm_field->tsCtx->get_preProcess_to_do_IJacobian().push_back(pre_only);
  }
  if(method) {
    dm_field->tsCtx->get_loops_to_do_IJacobian().push_back(PairNameFEMethodPtr(fe_name,method));
  }
  if(post_only) {
    dm_field->tsCtx->get_postProcess_to_do_IJacobian().push_back(post_only);
  }
  ierr = DMTSSetIJacobian(dm,f_TSSetIJacobian,dm_field->tsCtx.get()); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMTSSetIJacobian(DM dm,const char fe_name[],MoFEM::FEMethod *method,MoFEM::FEMethod *pre_only,MoFEM::FEMethod *post_only) {
  return DMMoFEMTSSetIJacobian<const char*,FEMethod*>(dm,fe_name,method,pre_only,post_only);
}

PetscErrorCode DMMoFEMTSSetIJacobian(
  DM dm,const std::string& fe_name,
  boost::shared_ptr<MoFEM::FEMethod> method,
  boost::shared_ptr<MoFEM::FEMethod> pre_only,
  boost::shared_ptr<MoFEM::FEMethod> post_only
) {
  return DMMoFEMTSSetIJacobian<const std::string&,boost::shared_ptr<MoFEM::FEMethod> >(
    dm,fe_name,method,pre_only,post_only
  );
}

PetscErrorCode DMMoFEMGetKspCtx(DM dm,MoFEM::KspCtx **ksp_ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  *ksp_ctx = dm_field->kspCtx.get();
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetKspCtx(DM dm,const boost::shared_ptr<MoFEM::KspCtx>& ksp_ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  const_cast<boost::shared_ptr<MoFEM::KspCtx>& >(ksp_ctx) = dm_field->kspCtx;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMSetKspCtx(DM dm,boost::shared_ptr<MoFEM::KspCtx>& ksp_ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  dm_field->kspCtx = ksp_ctx;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetSnesCtx(DM dm,MoFEM::SnesCtx **snes_ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  *snes_ctx = dm_field->snesCtx.get();
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetSnesCtx(DM dm,const boost::shared_ptr<MoFEM::SnesCtx>& snes_ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  const_cast<boost::shared_ptr<MoFEM::SnesCtx>& >(snes_ctx) = dm_field->snesCtx;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMSetSnesCtx(DM dm,boost::shared_ptr<MoFEM::SnesCtx> &snes_ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  dm_field->snesCtx = snes_ctx;
  MoFEMFunctionReturnHot(0);
}

/** get if read mesh is partitioned
  * \ingroup dm
  */
PetscErrorCode DMMoFEMSetIsPartitioned(DM dm,PetscBool is_partitioned) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  dm_field->isPartitioned = is_partitioned;
  MoFEMFunctionReturnHot(0);
}

/** get if read mesh is partitioned
  * \ingroup dm
  */
PetscErrorCode DMMoFEMGetIsPartitioned(DM dm,PetscBool *is_partitioned) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  *is_partitioned = dm_field->isPartitioned;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetTsCtx(DM dm,MoFEM::TsCtx **ts_ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  *ts_ctx = dm_field->tsCtx.get();
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetTsCtx(DM dm,const boost::shared_ptr<MoFEM::TsCtx>& ts_ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  const_cast<boost::shared_ptr<MoFEM::TsCtx>& >(ts_ctx) = dm_field->tsCtx;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMSetTsCtx(DM dm,boost::shared_ptr<MoFEM::TsCtx>& ts_ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  dm_field->tsCtx = ts_ctx;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMCreateGlobalVector_MoFEM(DM dm,Vec *g) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->query_interface<VecManager>()->vecCreateGhost(
    dm_field->problemName,COL,g
  ); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMCreateLocalVector_MoFEM(DM dm,Vec *l) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->query_interface<VecManager>()->vecCreateSeq(
    dm_field->problemName,COL,l
  ); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMCreateMatrix_MoFEM(DM dm,Mat *M) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  if(strcmp(dm->mattype,MATMPIAIJ)==0) {
    ierr = dm_field->mField_ptr->MatCreateMPIAIJWithArrays(dm_field->problemName,M); CHKERRQ(ierr);
  } else if(strcmp(dm->mattype,MATAIJ)==0) {
    PetscInt *i;
    PetscInt *j;
    PetscScalar *v;
    #if PETSC_VERSION_GE(3,7,0)
      ierr = dm_field->mField_ptr->MatCreateSeqAIJWithArrays(
        dm_field->problemName,M,&i,&j,&v
      ); CHKERRQ(ierr);
      ierr = MatConvert(*M,MATAIJ,MAT_INPLACE_MATRIX,M); CHKERRQ(ierr);
    #else
      Mat N;
      ierr = dm_field->mField_ptr->MatCreateSeqAIJWithArrays(
        dm_field->problemName,&N,&i,&j,&v
      ); CHKERRQ(ierr);
      ierr = MatConvert(N,MATAIJ,MAT_INITIAL_MATRIX,M); CHKERRQ(ierr);
      ierr = MatDestroy(&N); CHKERRQ(ierr);
    #endif
    ierr = PetscFree(i); CHKERRQ(ierr);
    ierr = PetscFree(j); CHKERRQ(ierr);
    ierr = PetscFree(v); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Matrix type not implemented");
  }
  MoFEMFunctionReturnHot(0);
}

#if PETSC_VERSION_GE(3,7,0)
PetscErrorCode DMSetFromOptions_MoFEM(PetscOptionItems *PetscOptionsObject,DM dm) {
#elif PETSC_VERSION_GE(3,5,3)
PetscErrorCode DMSetFromOptions_MoFEM(PetscOptions *PetscOptionsObject,DM dm) {
#else
PetscErrorCode DMSetFromOptions_MoFEM(DM dm) {
#endif

  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  #if PETSC_VERSION_GE(3,5,3)
    ierr = PetscOptionsHead(PetscOptionsObject,"DMMoFEM Options");CHKERRQ(ierr);
  #else
    ierr = PetscOptionsHead("DMMoFEM Options");CHKERRQ(ierr);
  #endif
  ierr = PetscOptionsBool(
    "-dm_is_partitioned",
    "set if mesh is partitioned (works which native MOAB file format, i.e. h5m","DMSetUp",
    dm_field->isPartitioned,&dm_field->isPartitioned,NULL
  ); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMSetUp_MoFEM(DM dm) {

  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ProblemsManager *prb_mng_ptr;
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->query_interface(prb_mng_ptr); CHKERRQ(ierr);

  if(dm_field->isCompDM) {
    // It is composite probelm
    ierr = prb_mng_ptr->buildCompsedProblem(
      dm_field->problemName,
      dm_field->rowCompPrb,
      dm_field->colCompPrb,
      dm_field->isSquareMatrix == PETSC_TRUE
    ); CHKERRQ(ierr);
  } else {
    if(dm_field->isPartitioned) {
      ierr = prb_mng_ptr->buildProblemOnDistributedMesh(
        dm_field->problemName,dm_field->isSquareMatrix == PETSC_TRUE
      ); CHKERRQ(ierr);
    } else {
      ierr = prb_mng_ptr->buildProblem(
        dm_field->problemName,dm_field->isSquareMatrix == PETSC_TRUE
      ); CHKERRQ(ierr);
      ierr = prb_mng_ptr->partitionProblem(dm_field->problemName); CHKERRQ(ierr);
    }
  }

  // Partition finite elements
  if(dm_field->isPartitioned) {
    ierr = prb_mng_ptr->partitionFiniteElements(
      dm_field->problemName,true,0,dm_field->sIze,1
    ); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionGhostDofsOnDistributedMesh(
      dm_field->problemName
    ); CHKERRQ(ierr);
  } else {
    ierr = prb_mng_ptr->partitionFiniteElements(
      dm_field->problemName
    ); CHKERRQ(ierr);
    // Get ghost DOFs
    ierr = prb_mng_ptr->partitionGhostDofs(dm_field->problemName); CHKERRQ(ierr);
  }

  // Set flag that problem is build and partitioned
  dm_field->isProblemBuild = PETSC_TRUE;

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMSubDMSetUp_MoFEM(DM subdm) {
  PetscValidHeaderSpecific(subdm,DM_CLASSID,1);
  ProblemsManager *prb_mng_ptr;
  MoFEMFunctionBeginHot;

  DMCtx *subdm_field = (DMCtx*)subdm->data;

  // build sub dm problem
  ierr = subdm_field->mField_ptr->query_interface(prb_mng_ptr); CHKERRQ(ierr);
  ierr = prb_mng_ptr->buildSubProblem(
    subdm_field->problemName,
    subdm_field->rowFields,
    subdm_field->colFields,
    subdm_field->problemMainOfSubPtr->getName(),
    subdm_field->isSquareMatrix == PETSC_TRUE
  ); CHKERRQ(ierr);

  // partition problem
  subdm_field->isPartitioned = subdm_field->isPartitioned;
  if(subdm_field->isPartitioned) {
    ierr = prb_mng_ptr->partitionFiniteElements(
      subdm_field->problemName,true,0,subdm_field->sIze,1
    ); CHKERRQ(ierr);
    // set ghost nodes
    ierr = prb_mng_ptr->partitionGhostDofsOnDistributedMesh(
      subdm_field->problemName
    ); CHKERRQ(ierr);
  } else {
    ierr = prb_mng_ptr->partitionFiniteElements(
      subdm_field->problemName
    ); CHKERRQ(ierr);
    // set ghost nodes
    ierr = prb_mng_ptr->partitionGhostDofs(subdm_field->problemName); CHKERRQ(ierr);
  }

  subdm_field->isProblemBuild = PETSC_TRUE;

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMGlobalToLocalBegin_MoFEM(DM dm,Vec g,InsertMode mode,Vec l) {
  MoFEMFunctionBeginHot;


  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;

  ierr = VecGhostUpdateBegin(g,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMGlobalToLocalEnd_MoFEM(DM dm,Vec g,InsertMode mode,Vec l) {
  MoFEMFunctionBeginHot;

  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;

  ierr = VecGhostUpdateEnd(g,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  DMCtx *dm_field = (DMCtx*)dm->data;
  int nb_dofs = dm_field->problemPtr->getNbLocalDofsRow();
  int nb_ghost = dm_field->problemPtr->getNbGhostDofsRow();

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
  MoFEMFunctionReturnHot(0);
}


PetscErrorCode DMLocalToGlobalBegin_MoFEM(DM dm,Vec l,InsertMode mode,Vec g) {

  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;

  DMCtx *dm_field = (DMCtx*)dm->data;
  int nb_dofs = dm_field->problemPtr->getNbLocalDofsRow();
  int nb_ghost = dm_field->problemPtr->getNbGhostDofsRow();

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

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMLocalToGlobalEnd_MoFEM(DM,Vec l,InsertMode mode,Vec g) {
  //
  MoFEMFunctionBeginHot;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMCreateFieldIS_MoFEM(
  DM dm,PetscInt *numFields,char ***fieldNames, IS **fields
) {

  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;

  if (numFields) {
    *numFields = 0;
  }
  if (fieldNames) {
    *fieldNames = NULL;
  }
  if (fields) {
    *fields = NULL;
  }

  DMCtx *dm_field = (DMCtx*)dm->data;
  const Field_multiIndex *fields_ptr;
  ierr = dm_field->mField_ptr->get_fields(&fields_ptr); CHKERRQ(ierr);
  Field_multiIndex::iterator fit,hi_fit;
  fit = fields_ptr->begin();
  hi_fit = fields_ptr->end();
  *numFields = std::distance(fit,hi_fit);

  if(fieldNames) {
    PetscMalloc1(*numFields,fieldNames);
  }
  if(fields) {
    PetscMalloc1(*numFields,fields);
  }

  for(int f = 0;fit!=hi_fit;fit++,f++) {
    if(fieldNames) {
      PetscStrallocpy(fit->get()->getName().c_str(),(char**)&(*fieldNames)[f]);
    }
    if(fields) {
      ierr = dm_field->mField_ptr->query_interface<ISManager>()->isCreateProblemFieldAndRank(
        dm_field->problemPtr->getName(),
        ROW,fit->get()->getName(),0,fit->get()->getNbOfCoeffs(),
        &(*fields)[f]
      ); CHKERRQ(ierr);
    }
  }

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetFieldIS(DM dm,RowColData rc,const char field_name[],IS *is) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx*)dm->data;
  ierr = dm_field->mField_ptr->query_interface<ISManager>()->isCreateProblemFieldAndRank(
    dm_field->problemPtr->getName(),
    ROW,field_name,0,1000,
    is
  ); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}
