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

// #undef PETSC_VERSION_RELEASE
// #define PETSC_VERSION_RELEASE 1

#if PETSC_VERSION_GE(3, 6, 0)
#include <petsc/private/dmimpl.h> /*I  "petscdm.h"   I*/
// #include <petsc/private/vecimpl.h> /*I  "petscdm.h"   I*/
#else
#include <petsc-private/dmimpl.h>  /*I  "petscdm.h"   I*/
#include <petsc-private/vecimpl.h> /*I  "petscdm.h"   I*/
#endif

#include <DMMoFEM.hpp>

namespace MoFEM {

DMCtx::DMCtx()
    : mField_ptr(PETSC_NULL), isProblemBuild(PETSC_FALSE),
      isPartitioned(PETSC_FALSE), isSquareMatrix(PETSC_TRUE),
      isSubDM(PETSC_FALSE), isCompDM(PETSC_FALSE), destroyProblem(PETSC_FALSE),
      verbosity(VERBOSE), referenceNumber(0) {

  if (!LogManager::checkIfChannelExist("DMWORLD")) {
    auto core_log = logging::core::get();
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "DMWORLD"));
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSync(), "DMSYNC"));
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSelf(), "DMSELF"));
    LogManager::setLog("DMWORLD");
    LogManager::setLog("DMSYNC");
    LogManager::setLog("DMSELF");
    MOFEM_LOG_TAG("DMWORLD", "DM");
    MOFEM_LOG_TAG("DMSYNC", "DM");
    MOFEM_LOG_TAG("DMSELF", "DM");
  }
}

MoFEMErrorCode DMCtx::query_interface(boost::typeindex::type_index type_index,
                                      UnknownInterface **iface) const {
  *iface = const_cast<DMCtx *>(this);
  return 0;
}

PetscErrorCode DMRegister_MoFEM(const char sname[]) {
  MoFEMFunctionBegin;
  CHKERR DMRegister(sname, DMCreate_MoFEM);
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMSetOperators_MoFEM(DM dm) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;

  dm->ops->createglobalvector = DMCreateGlobalVector_MoFEM;
  dm->ops->createlocalvector = DMCreateLocalVector_MoFEM;
  dm->ops->creatematrix = DMCreateMatrix_MoFEM;
  dm->ops->setup = DMSetUp_MoFEM;
  dm->ops->destroy = DMDestroy_MoFEM;
  dm->ops->setfromoptions = DMSetFromOptions_MoFEM;
  dm->ops->globaltolocalbegin = DMGlobalToLocalBegin_MoFEM;
  dm->ops->globaltolocalend = DMGlobalToLocalEnd_MoFEM;
  dm->ops->localtoglobalbegin = DMLocalToGlobalBegin_MoFEM;
  dm->ops->localtoglobalend = DMLocalToGlobalEnd_MoFEM;
  dm->ops->createfieldis = DMCreateFieldIS_MoFEM;

  // Default matrix type
  CHKERR DMSetMatType(dm, MATMPIAIJ);

  MoFEMFunctionReturn(0);
}

PetscErrorCode DMCreate_MoFEM(DM dm) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  dm->data = new DMCtx();
  CHKERR DMSetOperators_MoFEM(dm);
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMDestroy_MoFEM(DM dm) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  MoFEMFunctionBeginHot;

  MPI_Comm comm;
  CHKERR PetscObjectGetComm((PetscObject)dm, &comm);

  int result;
  MPI_Comm_compare(comm, PETSC_COMM_SELF, &result);
  if (result == MPI_IDENT)
    MOFEM_LOG("DMSELF", Sev::noisy)
        << "MoFEM DM destroy for problem " << dm_field->problemName
        << " referenceNumber " << dm_field->referenceNumber;
  else
    MOFEM_LOG("DMWORLD", Sev::noisy)
        << "MoFEM DM destroy for problem " << dm_field->problemName
        << " referenceNumber " << dm_field->referenceNumber;

  if (dm_field->referenceNumber == 0) {
    if (dm_field->destroyProblem) {

      if (dm_field->mField_ptr->check_problem(dm_field->problemName)) {
        dm_field->mField_ptr->delete_problem(dm_field->problemName);
      } // else problem has to be deleted by the user
    }

    delete static_cast<DMCtx *>(dm->data);

  } else
    --dm_field->referenceNumber;

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMCreateMoFEM(DM dm, MoFEM::Interface *m_field_ptr,
                                  const char problem_name[],
                                  const MoFEM::BitRefLevel bit_level,
                                  const MoFEM::BitRefLevel bit_mask) {
  MoFEMFunctionBegin;

  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (!dm->data) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "data structure for MoFEM not yet created");
  }
  if (!m_field_ptr) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "DM function not implemented into MoFEM");
  }
  dm_field->mField_ptr = m_field_ptr;
  dm_field->problemName = problem_name;
  if (!m_field_ptr->check_problem(dm_field->problemName)) {
    // problem is not defined, declare problem internally set bool to
    // destroyProblem problem with DM
    dm_field->destroyProblem = PETSC_TRUE;
    CHKERR dm_field->mField_ptr->add_problem(dm_field->problemName, MF_EXCL,
                                             dm_field->verbosity);
  } else {
    dm_field->destroyProblem = PETSC_FALSE;
  }
  CHKERR dm_field->mField_ptr->modify_problem_ref_level_add_bit(
      dm_field->problemName, bit_level);
  CHKERR dm_field->mField_ptr->modify_problem_mask_ref_level_add_bit(
      dm_field->problemName, bit_mask);
  dm_field->kspCtx =
      boost::shared_ptr<KspCtx>(new KspCtx(*m_field_ptr, problem_name));
  dm_field->snesCtx =
      boost::shared_ptr<SnesCtx>(new SnesCtx(*m_field_ptr, problem_name));
  dm_field->tsCtx =
      boost::shared_ptr<TsCtx>(new TsCtx(*m_field_ptr, problem_name));

  MPI_Comm comm;
  CHKERR PetscObjectGetComm((PetscObject)dm, &comm);
  int result = 0;
  MPI_Comm_compare(comm, m_field_ptr->get_comm(), &result);
  if (result > MPI_CONGRUENT) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "MoFEM and DM using different communicators");
  }
  MPI_Comm_size(comm, &dm_field->sIze);
  MPI_Comm_rank(comm, &dm_field->rAnk);

  // problem structure
  CHKERR dm_field->mField_ptr->get_problem(dm_field->problemName,
                                           &dm_field->problemPtr);

  MPI_Comm_compare(comm, PETSC_COMM_SELF, &result);
  if (result == MPI_IDENT)
    MOFEM_LOG("DMSELF", Sev::noisy)
        << "MoFEM DM created for problem " << dm_field->problemName;
  else
    MOFEM_LOG("DMWORLD", Sev::noisy)
        << "MoFEM DM created for problem " << dm_field->problemName;

  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMDuplicateDMCtx(DM dm, DM dm_duplicate) {
  MoFEMFunctionBegin;

  auto *dm_field = static_cast<DMCtx *>(dm->data);
  if (!dm->data)
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "data structure for MoFEM not yet created");

  if (static_cast<DMCtx *>(dm_duplicate->data)->referenceNumber == 0)
    delete static_cast<DMCtx *>(dm_duplicate->data);

  dm_duplicate->data = dm->data;
  ++(static_cast<DMCtx *>(dm_duplicate->data)->referenceNumber);

  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMCreateSubDM(DM subdm, DM dm, const char problem_name[]) {
  MoFEMFunctionBegin;

  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (!dm->data) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "data structure for MoFEM not yet created");
  }
  CHKERR DMMoFEMCreateMoFEM(subdm, dm_field->mField_ptr, problem_name,
                            dm_field->problemPtr->getBitRefLevel());

  DMCtx *subdm_field = (DMCtx *)subdm->data;
  subdm_field->isSubDM = PETSC_TRUE;
  subdm_field->problemMainOfSubPtr = dm_field->problemPtr;
  subdm_field->isPartitioned = dm_field->isPartitioned;
  subdm_field->isSquareMatrix = PETSC_FALSE;
  subdm->ops->setup = DMSubDMSetUp_MoFEM;

  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMAddSubFieldRow(DM dm, const char field_name[],
                                     EntityType lo_type, EntityType hi_type) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (!dm->data) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "data structure for MoFEM not yet created");
  }
  if (!dm_field->isSubDM) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "this is not sub-dm");
  }
  dm_field->rowFields.push_back(field_name);
  if (lo_type != MBVERTEX || hi_type != MBMAXTYPE) {
    if (!dm_field->mapTypeRow)
      dm_field->mapTypeRow = boost::make_shared<
          std::map<std::string, std::pair<EntityType, EntityType>>>();
    (*dm_field->mapTypeRow)[field_name] =
        std::pair<EntityType, EntityType>(lo_type, hi_type);
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMAddSubFieldCol(DM dm, const char field_name[],
                                     EntityType lo_type, EntityType hi_type) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (!dm->data) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "data structure for MoFEM not yet created");
  }
  if (!dm_field->isSubDM) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "this is not sub-dm");
  }
  dm_field->colFields.push_back(field_name);
  if (lo_type != MBVERTEX || hi_type != MBMAXTYPE) {
    if (!dm_field->mapTypeCol)
      dm_field->mapTypeCol = boost::make_shared<
          std::map<std::string, std::pair<EntityType, EntityType>>>();
    (*dm_field->mapTypeCol)[field_name] =
        std::pair<EntityType, EntityType>(lo_type, hi_type);
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetIsSubDM(DM dm, PetscBool *is_sub_dm) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  *is_sub_dm = dm_field->isSubDM;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetSubRowIS(DM dm, IS *is) {
  MoFEMFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (dm_field->isSubDM != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "This DM is not created as a SubDM");
  }
  if (dm_field->isProblemBuild != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Problem is not build");
  }
  boost::shared_ptr<Problem::SubProblemData> sub_data =
      dm_field->problemPtr->getSubData();
  CHKERR sub_data->getRowIs(is);
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMGetSubColIS(DM dm, IS *is) {
  MoFEMFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (dm_field->isSubDM != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "This DM is not created as a SubDM");
  }
  if (dm_field->isProblemBuild != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Problem is not build");
  }
  boost::shared_ptr<Problem::SubProblemData> sub_data =
      dm_field->problemPtr->getSubData();
  CHKERR sub_data->getColIs(is);
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMAddRowCompositeProblem(DM dm, const char prb_name[]) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (!dm->data) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "data structure for MoFEM not yet created");
  }
  if (!dm_field->isCompDM) {
    dm_field->isCompDM = PETSC_TRUE;
  }
  dm_field->rowCompPrb.push_back(prb_name);
  if (dm_field->isSquareMatrix) {
    dm_field->colCompPrb.push_back(prb_name);
  }
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMAddColCompositeProblem(DM dm, const char prb_name[]) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (!dm->data) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "data structure for MoFEM not yet created");
  }
  if (!dm_field->isCompDM) {
    dm_field->isCompDM = PETSC_TRUE;
  }
  if (dm_field->isSquareMatrix) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
            "No need to add problem on column when problem block structurally "
            "symmetric");
  }
  dm_field->colCompPrb.push_back(prb_name);
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMGetIsCompDM(DM dm, PetscBool *is_comp_dm) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  *is_comp_dm = dm_field->isCompDM;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMGetInterfacePtr(DM dm, MoFEM::Interface **m_field_ptr) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (!dm->data) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "data structure for MoFEM not yet created");
  }
  *m_field_ptr = dm_field->mField_ptr;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetProblemPtr(DM dm, const MoFEM::Problem **problem_ptr) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (!dm->data) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "data structure for MoFEM not yet created");
  }
  *problem_ptr = dm_field->problemPtr;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMSetDestroyProblem(DM dm, PetscBool destroy_problem) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  dm_field->destroyProblem = destroy_problem;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetDestroyProblem(DM dm, PetscBool *destroy_problem) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  *destroy_problem = dm_field->destroyProblem;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMSetSquareProblem(DM dm, PetscBool square_problem) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  dm_field->isSquareMatrix = square_problem;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMResolveSharedFiniteElements(DM dm, const char fe_name[]) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  CHKERR dm_field->mField_ptr->getInterface<CommInterface>()
      ->resolveSharedFiniteElements(dm_field->problemPtr, fe_name);
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMResolveSharedEntities(DM dm, const char fe_name[]) {
  return DMMoFEMResolveSharedFiniteElements(dm, fe_name);
}

PetscErrorCode DMMoFEMGetProblemFiniteElementLayout(DM dm, const char fe_name[],
                                                    PetscLayout *layout) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);

  MPI_Comm comm;
  CHKERR PetscObjectGetComm((PetscObject)dm, &comm);
  CHKERR dm_field->problemPtr->getNumberOfElementsByNameAndPart(comm, fe_name,
                                                                layout);
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMGetSquareProblem(DM dm, PetscBool *square_problem) {
  MoFEMFunctionBeginHot;
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  *square_problem = dm_field->isSquareMatrix;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMAddElement(DM dm, const char fe_name[]) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  ierr = dm_field->mField_ptr->modify_problem_add_finite_element(
      dm_field->problemName, fe_name);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMUnSetElement(DM dm, const char fe_name[]) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  ierr = dm_field->mField_ptr->modify_problem_unset_finite_element(
      dm_field->problemName, fe_name);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMMeshToLocalVector(DM dm, Vec l, InsertMode mode,
                                       ScatterMode scatter_mode) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  ierr = dm_field->mField_ptr->getInterface<VecManager>()->setLocalGhostVector(
      dm_field->problemPtr, COL, l, mode, scatter_mode);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMMeshToGlobalVector(DM dm, Vec g, InsertMode mode,
                                        ScatterMode scatter_mode) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  ierr = dm_field->mField_ptr->getInterface<VecManager>()->setGlobalGhostVector(
      dm_field->problemPtr, COL, g, mode, scatter_mode);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMPreProcessFiniteElements(DM dm, MoFEM::FEMethod *method) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  ierr = dm_field->mField_ptr->problem_basic_method_preProcess(
      dm_field->problemPtr, *method);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMPostProcessFiniteElements(DM dm, MoFEM::FEMethod *method) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  ierr = dm_field->mField_ptr->problem_basic_method_postProcess(
      dm_field->problemPtr, *method);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMoFEMLoopFiniteElementsUpAndLowRank(DM dm, const char fe_name[],
                                                    MoFEM::FEMethod *method,
                                                    int low_rank, int up_rank) {
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  ierr = dm_field->mField_ptr->loop_finite_elements(
      dm_field->problemPtr, fe_name, *method, low_rank, up_rank);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode
DMoFEMLoopFiniteElementsUpAndLowRank(DM dm, const std::string fe_name,
                                     boost::shared_ptr<MoFEM::FEMethod> method,
                                     int low_rank, int up_rank) {
  return DMoFEMLoopFiniteElementsUpAndLowRank(dm, fe_name.c_str(), method.get(),
                                              low_rank, up_rank);
}

PetscErrorCode DMoFEMLoopFiniteElements(DM dm, const char fe_name[],
                                        MoFEM::FEMethod *method) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  ierr = DMoFEMLoopFiniteElementsUpAndLowRank(dm, fe_name, method,
                                              dm_field->rAnk, dm_field->rAnk);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode
DMoFEMLoopFiniteElements(DM dm, const std::string fe_name,
                         boost::shared_ptr<MoFEM::FEMethod> method) {
  return DMoFEMLoopFiniteElements(dm, fe_name.c_str(), method.get());
}

PetscErrorCode DMoFEMLoopDofs(DM dm, const char field_name[],
                              MoFEM::DofMethod *method) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  ierr =
      dm_field->mField_ptr->loop_dofs(dm_field->problemPtr, field_name, COL,
                                      *method, dm_field->rAnk, dm_field->rAnk);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

template <class S, class T0, class T1, class T2>
static PetscErrorCode DMMoFEMKSPSetComputeRHS(DM dm, S fe_name, T0 method,
                                              T1 pre_only, T2 post_only) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (pre_only) {
    dm_field->kspCtx->get_preProcess_to_do_Rhs().push_back(pre_only);
  }
  if (method) {
    dm_field->kspCtx->get_loops_to_do_Rhs().push_back(
        PairNameFEMethodPtr(fe_name, method));
  }
  if (post_only) {
    dm_field->kspCtx->get_postProcess_to_do_Rhs().push_back(post_only);
  }
  CHKERR DMKSPSetComputeRHS(dm, KspRhs, dm_field->kspCtx.get());
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMKSPSetComputeRHS(DM dm, const char fe_name[],
                                       MoFEM::FEMethod *method,
                                       MoFEM::BasicMethod *pre_only,
                                       MoFEM::BasicMethod *post_only) {
  return DMMoFEMKSPSetComputeRHS<const char *, MoFEM::FEMethod *,
                                 MoFEM::BasicMethod *, MoFEM::BasicMethod *>(
      dm, fe_name, method, pre_only, post_only);
}

PetscErrorCode
DMMoFEMKSPSetComputeRHS(DM dm, const std::string fe_name,
                        boost::shared_ptr<MoFEM::FEMethod> method,
                        boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                        boost::shared_ptr<MoFEM::BasicMethod> post_only) {
  return DMMoFEMKSPSetComputeRHS<const std::string,
                                 boost::shared_ptr<MoFEM::FEMethod>,
                                 boost::shared_ptr<MoFEM::BasicMethod>,
                                 boost::shared_ptr<MoFEM::BasicMethod>>(
      dm, fe_name, method, pre_only, post_only);
}

template <class S, class T0, class T1, class T2>
static PetscErrorCode DMMoFEMKSPSetComputeOperators(DM dm, S fe_name, T0 method,
                                                    T1 pre_only, T2 post_only) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (pre_only) {
    dm_field->kspCtx->get_preProcess_to_do_Mat().push_back(pre_only);
  }
  if (method) {
    dm_field->kspCtx->get_loops_to_do_Mat().push_back(
        PairNameFEMethodPtr(fe_name, method));
  }
  if (post_only) {
    dm_field->kspCtx->get_postProcess_to_do_Mat().push_back(post_only);
  }
  CHKERR DMKSPSetComputeOperators(dm, KspMat, dm_field->kspCtx.get());
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMKSPSetComputeOperators(DM dm, const char fe_name[],
                                             MoFEM::FEMethod *method,
                                             MoFEM::BasicMethod *pre_only,
                                             MoFEM::BasicMethod *post_only) {
  return DMMoFEMKSPSetComputeOperators<const char *, MoFEM::FEMethod *,
                                       MoFEM::BasicMethod *,
                                       MoFEM::BasicMethod *>(
      dm, fe_name, method, pre_only, post_only);
}

PetscErrorCode
DMMoFEMKSPSetComputeOperators(DM dm, const std::string fe_name,
                              boost::shared_ptr<MoFEM::FEMethod> method,
                              boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                              boost::shared_ptr<MoFEM::BasicMethod> post_only) {
  return DMMoFEMKSPSetComputeOperators<const std::string,
                                       boost::shared_ptr<MoFEM::FEMethod>>(
      dm, fe_name, method, pre_only, post_only);
}

template <class S, class T0, class T1, class T2>
static PetscErrorCode DMMoFEMSNESSetFunction(DM dm, S fe_name, T0 method,
                                             T1 pre_only, T2 post_only) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (pre_only) {
    dm_field->snesCtx->get_preProcess_to_do_Rhs().push_back(pre_only);
  }
  if (method) {
    dm_field->snesCtx->get_loops_to_do_Rhs().push_back(
        PairNameFEMethodPtr(fe_name, method));
  }
  if (post_only) {
    dm_field->snesCtx->get_postProcess_to_do_Rhs().push_back(post_only);
  }
  CHKERR DMSNESSetFunction(dm, SnesRhs, dm_field->snesCtx.get());
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMSNESSetFunction(DM dm, const char fe_name[],
                                      MoFEM::FEMethod *method,
                                      MoFEM::BasicMethod *pre_only,
                                      MoFEM::BasicMethod *post_only) {
  return DMMoFEMSNESSetFunction<const char *, MoFEM::FEMethod *,
                                MoFEM::BasicMethod *, MoFEM::BasicMethod *>(
      dm, fe_name, method, pre_only, post_only);
}

PetscErrorCode
DMMoFEMSNESSetFunction(DM dm, const std::string fe_name,
                       boost::shared_ptr<MoFEM::FEMethod> method,
                       boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                       boost::shared_ptr<MoFEM::BasicMethod> post_only) {
  return DMMoFEMSNESSetFunction<const std::string,
                                boost::shared_ptr<MoFEM::FEMethod>,
                                boost::shared_ptr<MoFEM::BasicMethod>,
                                boost::shared_ptr<MoFEM::BasicMethod>>(
      dm, fe_name, method, pre_only, post_only);
}

template <class S, class T0, class T1, class T2>
static PetscErrorCode DMMoFEMSNESSetJacobian(DM dm, S fe_name, T0 method,
                                             T1 pre_only, T2 post_only) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (pre_only) {
    dm_field->snesCtx->get_preProcess_to_do_Mat().push_back(pre_only);
  }
  if (method) {
    dm_field->snesCtx->get_loops_to_do_Mat().push_back(
        PairNameFEMethodPtr(fe_name, method));
  }
  if (post_only) {
    dm_field->snesCtx->get_postProcess_to_do_Mat().push_back(post_only);
  }
  CHKERR DMSNESSetJacobian(dm, SnesMat, dm_field->snesCtx.get());
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMSNESSetJacobian(DM dm, const char fe_name[],
                                      MoFEM::FEMethod *method,
                                      MoFEM::BasicMethod *pre_only,
                                      MoFEM::BasicMethod *post_only) {
  return DMMoFEMSNESSetJacobian<const char *, MoFEM::FEMethod *,
                                MoFEM::BasicMethod *, MoFEM::BasicMethod *>(
      dm, fe_name, method, pre_only, post_only);
}

PetscErrorCode
DMMoFEMSNESSetJacobian(DM dm, const std::string fe_name,
                       boost::shared_ptr<MoFEM::FEMethod> method,
                       boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                       boost::shared_ptr<MoFEM::BasicMethod> post_only) {
  return DMMoFEMSNESSetJacobian<const std::string,
                                boost::shared_ptr<MoFEM::FEMethod>,
                                boost::shared_ptr<MoFEM::BasicMethod>,
                                boost::shared_ptr<MoFEM::BasicMethod>>(
      dm, fe_name, method, pre_only, post_only);
}

template <class S, class T0, class T1, class T2>
static PetscErrorCode DMMoFEMTSSetIFunction(DM dm, S fe_name, T0 method,
                                            T1 pre_only, T2 post_only) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (pre_only) {
    dm_field->tsCtx->get_preProcess_to_do_IFunction().push_back(pre_only);
  }
  if (method) {
    dm_field->tsCtx->get_loops_to_do_IFunction().push_back(
        PairNameFEMethodPtr(fe_name, method));
  }
  if (post_only) {
    dm_field->tsCtx->get_postProcess_to_do_IFunction().push_back(post_only);
  }
  CHKERR DMTSSetIFunction(dm, TsSetIFunction, dm_field->tsCtx.get());
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMTSSetIFunction(DM dm, const char fe_name[],
                                     MoFEM::FEMethod *method,
                                     MoFEM::BasicMethod *pre_only,
                                     MoFEM::BasicMethod *post_only) {
  return DMMoFEMTSSetIFunction<const char *, MoFEM::FEMethod *,
                               MoFEM::BasicMethod *, MoFEM::BasicMethod *>(
      dm, fe_name, method, pre_only, post_only);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode
DMMoFEMTSSetIFunction(DM dm, const std::string fe_name,
                      boost::shared_ptr<MoFEM::FEMethod> method,
                      boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                      boost::shared_ptr<MoFEM::BasicMethod> post_only) {
  return DMMoFEMTSSetIFunction<const std::string,
                               boost::shared_ptr<MoFEM::FEMethod>,
                               boost::shared_ptr<MoFEM::BasicMethod>,
                               boost::shared_ptr<MoFEM::BasicMethod>>(
      dm, fe_name, method, pre_only, post_only);
  MoFEMFunctionReturnHot(0);
}

template <class S, class T0, class T1, class T2>
static PetscErrorCode DMMoFEMTSSetIJacobian(DM dm, S fe_name, T0 method,
                                            T1 pre_only, T2 post_only) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (pre_only) {
    dm_field->tsCtx->get_preProcess_to_do_IJacobian().push_back(pre_only);
  }
  if (method) {
    dm_field->tsCtx->get_loops_to_do_IJacobian().push_back(
        PairNameFEMethodPtr(fe_name, method));
  }
  if (post_only) {
    dm_field->tsCtx->get_postProcess_to_do_IJacobian().push_back(post_only);
  }
  CHKERR DMTSSetIJacobian(dm, TsSetIJacobian, dm_field->tsCtx.get());
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMTSSetIJacobian(DM dm, const char fe_name[],
                                     MoFEM::FEMethod *method,
                                     MoFEM::BasicMethod *pre_only,
                                     MoFEM::BasicMethod *post_only) {
  return DMMoFEMTSSetIJacobian<const char *, FEMethod *, MoFEM::BasicMethod *,
                               MoFEM::BasicMethod *>(dm, fe_name, method,
                                                     pre_only, post_only);
}

PetscErrorCode
DMMoFEMTSSetIJacobian(DM dm, const std::string fe_name,
                      boost::shared_ptr<MoFEM::FEMethod> method,
                      boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                      boost::shared_ptr<MoFEM::BasicMethod> post_only) {
  return DMMoFEMTSSetIJacobian<const std::string,
                               boost::shared_ptr<MoFEM::FEMethod>,
                               boost::shared_ptr<MoFEM::BasicMethod>,
                               boost::shared_ptr<MoFEM::BasicMethod>>(
      dm, fe_name, method, pre_only, post_only);
}

template <class S, class T0, class T1, class T2>
static PetscErrorCode DMMoFEMTSSetRHSFunction(DM dm, S fe_name, T0 method,
                                              T1 pre_only, T2 post_only) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (pre_only)
    dm_field->tsCtx->get_preProcess_to_do_RHSFunction().push_back(pre_only);
  if (method)
    dm_field->tsCtx->get_loops_to_do_RHSFunction().push_back(
        PairNameFEMethodPtr(fe_name, method));
  if (post_only)
    dm_field->tsCtx->get_postProcess_to_do_RHSFunction().push_back(post_only);
  CHKERR DMTSSetRHSFunction(dm, TsSetRHSFunction, dm_field->tsCtx.get());
  MoFEMFunctionReturn(0);
}

PetscErrorCode
DMMoFEMTSSetRHSFunction(DM dm, const std::string fe_name,
                        boost::shared_ptr<MoFEM::FEMethod> method,
                        boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                        boost::shared_ptr<MoFEM::BasicMethod> post_only) {
  return DMMoFEMTSSetRHSFunction<const std::string,
                                 boost::shared_ptr<MoFEM::FEMethod>,
                                 boost::shared_ptr<MoFEM::BasicMethod>,
                                 boost::shared_ptr<MoFEM::BasicMethod>>(
      dm, fe_name, method, pre_only, post_only);
  MoFEMFunctionReturnHot(0);
}

template <class S, class T0, class T1, class T2>
static PetscErrorCode DMMoFEMTSSetRHSJacobian(DM dm, S fe_name, T0 method,
                                              T1 pre_only, T2 post_only) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (pre_only)
    dm_field->tsCtx->get_preProcess_to_do_RHSFunction().push_back(pre_only);
  if (method)
    dm_field->tsCtx->get_loops_to_do_RHSFunction().push_back(
        PairNameFEMethodPtr(fe_name, method));
  if (post_only)
    dm_field->tsCtx->get_postProcess_to_do_RHSFunction().push_back(post_only);
  CHKERR DMTSSetRHSJacobian(dm, TsSetRHSJacobian, dm_field->tsCtx.get());
  MoFEMFunctionReturn(0);
}

PetscErrorCode
DMMoFEMTSSetRHSJacobian(DM dm, const std::string fe_name,
                        boost::shared_ptr<MoFEM::FEMethod> method,
                        boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                        boost::shared_ptr<MoFEM::BasicMethod> post_only) {
  return DMMoFEMTSSetRHSJacobian<const std::string,
                                 boost::shared_ptr<MoFEM::FEMethod>,
                                 boost::shared_ptr<MoFEM::BasicMethod>,
                                 boost::shared_ptr<MoFEM::BasicMethod>>(
      dm, fe_name, method, pre_only, post_only);
  MoFEMFunctionReturnHot(0);
}

template <class S, class T0, class T1, class T2>
static PetscErrorCode DMMoFEMTSSetI2Function(DM dm, S fe_name, T0 method,
                                             T1 pre_only, T2 post_only) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (pre_only) {
    dm_field->tsCtx->get_preProcess_to_do_IFunction().push_back(pre_only);
  }
  if (method) {
    dm_field->tsCtx->get_loops_to_do_IFunction().push_back(
        PairNameFEMethodPtr(fe_name, method));
  }
  if (post_only) {
    dm_field->tsCtx->get_postProcess_to_do_IFunction().push_back(post_only);
  }
  CHKERR DMTSSetI2Function(dm, TsSetI2Function, dm_field->tsCtx.get());
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMTSSetI2Function(DM dm, const char fe_name[],
                                      MoFEM::FEMethod *method,
                                      MoFEM::BasicMethod *pre_only,
                                      MoFEM::BasicMethod *post_only) {
  return DMMoFEMTSSetI2Function<const char *, MoFEM::FEMethod *,
                                MoFEM::BasicMethod *, MoFEM::BasicMethod *>(
      dm, fe_name, method, pre_only, post_only);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode
DMMoFEMTSSetI2Function(DM dm, const std::string fe_name,
                       boost::shared_ptr<MoFEM::FEMethod> method,
                       boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                       boost::shared_ptr<MoFEM::BasicMethod> post_only) {
  return DMMoFEMTSSetI2Function<const std::string,
                                boost::shared_ptr<MoFEM::FEMethod>,
                                boost::shared_ptr<MoFEM::BasicMethod>,
                                boost::shared_ptr<MoFEM::BasicMethod>>(
      dm, fe_name, method, pre_only, post_only);
  MoFEMFunctionReturnHot(0);
}

template <class S, class T0, class T1, class T2>
static PetscErrorCode DMMoFEMTSSetI2Jacobian(DM dm, S fe_name, T0 method,
                                             T1 pre_only, T2 post_only) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (pre_only) {
    dm_field->tsCtx->get_preProcess_to_do_IJacobian().push_back(pre_only);
  }
  if (method) {
    dm_field->tsCtx->get_loops_to_do_IJacobian().push_back(
        PairNameFEMethodPtr(fe_name, method));
  }
  if (post_only) {
    dm_field->tsCtx->get_postProcess_to_do_IJacobian().push_back(post_only);
  }
  CHKERR DMTSSetI2Jacobian(dm, TsSetI2Jacobian, dm_field->tsCtx.get());
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMTSSetI2Jacobian(DM dm, const char fe_name[],
                                      MoFEM::FEMethod *method,
                                      MoFEM::BasicMethod *pre_only,
                                      MoFEM::BasicMethod *post_only) {
  return DMMoFEMTSSetI2Jacobian<const char *, FEMethod *, MoFEM::BasicMethod *,
                                MoFEM::BasicMethod *>(dm, fe_name, method,
                                                      pre_only, post_only);
}

PetscErrorCode
DMMoFEMTSSetI2Jacobian(DM dm, const std::string fe_name,
                       boost::shared_ptr<MoFEM::FEMethod> method,
                       boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                       boost::shared_ptr<MoFEM::BasicMethod> post_only) {
  return DMMoFEMTSSetI2Jacobian<const std::string,
                                boost::shared_ptr<MoFEM::FEMethod>,
                                boost::shared_ptr<MoFEM::BasicMethod>,
                                boost::shared_ptr<MoFEM::BasicMethod>>(
      dm, fe_name, method, pre_only, post_only);
}

template <class S, class T0, class T1, class T2>
static PetscErrorCode DMMoFEMTSSetMonitor(DM dm, TS ts, S fe_name, T0 method,
                                          T1 pre_only, T2 post_only) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (pre_only)
    dm_field->tsCtx->get_preProcess_to_do_Monitor().push_back(pre_only);
  if (method)
    dm_field->tsCtx->get_loops_to_do_Monitor().push_back(
        PairNameFEMethodPtr(fe_name, method));
  if (post_only)
    dm_field->tsCtx->get_postProcess_to_do_Monitor().push_back(post_only);
  CHKERR TSMonitorSet(ts, TsMonitorSet, dm_field->tsCtx.get(), PETSC_NULL);
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMTSSetMonitor(DM dm, TS ts, const char fe_name[],
                                   MoFEM::FEMethod *method,
                                   MoFEM::BasicMethod *pre_only,
                                   MoFEM::BasicMethod *post_only) {
  return DMMoFEMTSSetMonitor<const char *, MoFEM::FEMethod *,
                             MoFEM::BasicMethod *, MoFEM::BasicMethod *>(
      dm, ts, fe_name, method, pre_only, post_only);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode
DMMoFEMTSSetMonitor(DM dm, TS ts, const std::string fe_name,
                    boost::shared_ptr<MoFEM::FEMethod> method,
                    boost::shared_ptr<MoFEM::BasicMethod> pre_only,
                    boost::shared_ptr<MoFEM::BasicMethod> post_only) {
  return DMMoFEMTSSetMonitor<const std::string,
                             boost::shared_ptr<MoFEM::FEMethod>,
                             boost::shared_ptr<MoFEM::BasicMethod>,
                             boost::shared_ptr<MoFEM::BasicMethod>>(
      dm, ts, fe_name, method, pre_only, post_only);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetKspCtx(DM dm, MoFEM::KspCtx **ksp_ctx) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  *ksp_ctx = dm_field->kspCtx.get();
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode
DMMoFEMGetKspCtx(DM dm, const boost::shared_ptr<MoFEM::KspCtx> &ksp_ctx) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  const_cast<boost::shared_ptr<MoFEM::KspCtx> &>(ksp_ctx) = dm_field->kspCtx;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMSetKspCtx(DM dm,
                                boost::shared_ptr<MoFEM::KspCtx> ksp_ctx) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  dm_field->kspCtx = ksp_ctx;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetSnesCtx(DM dm, MoFEM::SnesCtx **snes_ctx) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = (DMCtx *)dm->data;
  *snes_ctx = dm_field->snesCtx.get();
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode
DMMoFEMGetSnesCtx(DM dm, const boost::shared_ptr<MoFEM::SnesCtx> &snes_ctx) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  const_cast<boost::shared_ptr<MoFEM::SnesCtx> &>(snes_ctx) = dm_field->snesCtx;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMSetSnesCtx(DM dm,
                                 boost::shared_ptr<MoFEM::SnesCtx> snes_ctx) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  dm_field->snesCtx = snes_ctx;
  MoFEMFunctionReturnHot(0);
}

/** get if read mesh is partitioned
 * \ingroup dm
 */
PetscErrorCode DMMoFEMSetIsPartitioned(DM dm, PetscBool is_partitioned) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  dm_field->isPartitioned = is_partitioned;
  MoFEMFunctionReturnHot(0);
}

/** get if read mesh is partitioned
 * \ingroup dm
 */
PetscErrorCode DMMoFEMGetIsPartitioned(DM dm, PetscBool *is_partitioned) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  *is_partitioned = dm_field->isPartitioned;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetTsCtx(DM dm, MoFEM::TsCtx **ts_ctx) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  *ts_ctx = dm_field->tsCtx.get();
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMGetTsCtx(DM dm,
                               const boost::shared_ptr<MoFEM::TsCtx> &ts_ctx) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  const_cast<boost::shared_ptr<MoFEM::TsCtx> &>(ts_ctx) = dm_field->tsCtx;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMMoFEMSetTsCtx(DM dm, boost::shared_ptr<MoFEM::TsCtx> ts_ctx) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  dm_field->tsCtx = ts_ctx;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMCreateGlobalVector_MoFEM(DM dm, Vec *g) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  CHKERR dm_field->mField_ptr->getInterface<VecManager>()->vecCreateGhost(
      dm_field->problemName, COL, g);
  CHKERR VecSetDM(*g, dm);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMCreateGlobalVector_MoFEM(DM dm, SmartPetscObj<Vec> &g_ptr) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  CHKERR dm_field->mField_ptr->getInterface<VecManager>()->vecCreateGhost(
      dm_field->problemName, COL, g_ptr);
  CHKERR VecSetDM(g_ptr, dm);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMCreateLocalVector_MoFEM(DM dm, Vec *l) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  CHKERR dm_field->mField_ptr->getInterface<VecManager>()->vecCreateSeq(
      dm_field->problemName, COL, l);
  CHKERR VecSetDM(*l, dm);
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMCreateMatrix_MoFEM(DM dm, Mat *M) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (strcmp(dm->mattype, MATMPIAIJ) == 0) {
    CHKERR dm_field->mField_ptr->getInterface<MatrixManager>()
        ->createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>(dm_field->problemName,
                                                        M);
  } else if (strcmp(dm->mattype, MATAIJ) == 0) {
    CHKERR dm_field->mField_ptr->getInterface<MatrixManager>()
        ->createSeqAIJWithArrays<PetscLocalIdx_mi_tag>(dm_field->problemName,
                                                       M);
  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Matrix type not implemented");
  }
  CHKERR MatSetDM(*M, dm);
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMCreateMatrix_MoFEM(DM dm, SmartPetscObj<Mat> &M) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  if (strcmp(dm->mattype, MATMPIAIJ) == 0) {
    CHKERR dm_field->mField_ptr->getInterface<MatrixManager>()
        ->createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>(dm_field->problemName,
                                                        M);
  } else if (strcmp(dm->mattype, MATAIJ) == 0) {
    CHKERR dm_field->mField_ptr->getInterface<MatrixManager>()
        ->createSeqAIJWithArrays<PetscLocalIdx_mi_tag>(dm_field->problemName,
                                                       M);
  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Matrix type not implemented");
  }
  CHKERR MatSetDM(M, dm);
  MoFEMFunctionReturn(0);
}

#if PETSC_VERSION_GE(3, 7, 0)
PetscErrorCode DMSetFromOptions_MoFEM(PetscOptionItems *PetscOptionsObject,
                                      DM dm) {
#elif PETSC_VERSION_GE(3, 5, 3)
PetscErrorCode DMSetFromOptions_MoFEM(PetscOptions *PetscOptionsObject, DM dm) {
#else
PetscErrorCode DMSetFromOptions_MoFEM(DM dm) {
#endif

  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
#if PETSC_VERSION_GE(3, 5, 3)
  ierr = PetscOptionsHead(PetscOptionsObject, "DMMoFEM Options");
  CHKERRG(ierr);
#else
  ierr = PetscOptionsHead("DMMoFEM Options");
  CHKERRG(ierr);
#endif
  ierr = PetscOptionsBool("-dm_is_partitioned",
                          "set if mesh is partitioned (works which native MOAB "
                          "file format, i.e. h5m",
                          "DMSetUp", dm_field->isPartitioned,
                          &dm_field->isPartitioned, NULL);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMSetUp_MoFEM(DM dm) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  ProblemsManager *prb_mng_ptr;
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  CHKERR dm_field->mField_ptr->getInterface(prb_mng_ptr);

  if (dm_field->isCompDM) {
    // It is composite probelm
    CHKERR prb_mng_ptr->buildCompsedProblem(
        dm_field->problemName, dm_field->rowCompPrb, dm_field->colCompPrb,
        dm_field->isSquareMatrix == PETSC_TRUE, dm_field->verbosity);
  } else {
    if (dm_field->isPartitioned) {
      CHKERR prb_mng_ptr->buildProblemOnDistributedMesh(
          dm_field->problemName, dm_field->isSquareMatrix == PETSC_TRUE,
          dm_field->verbosity);
    } else {
      CHKERR prb_mng_ptr->buildProblem(dm_field->problemName,
                                       dm_field->isSquareMatrix == PETSC_TRUE,
                                       dm_field->verbosity);
      CHKERR prb_mng_ptr->partitionProblem(dm_field->problemName,
                                           dm_field->verbosity);
    }
  }

  // Partition finite elements
  if (dm_field->isPartitioned) {
    CHKERR prb_mng_ptr->partitionFiniteElements(
        dm_field->problemName, true, 0, dm_field->sIze, dm_field->verbosity);
    CHKERR prb_mng_ptr->partitionGhostDofsOnDistributedMesh(
        dm_field->problemName, dm_field->verbosity);
  } else {
    // partition finite elemnets
    CHKERR prb_mng_ptr->partitionFiniteElements(dm_field->problemName, false,
                                                -1, -1, dm_field->verbosity);
    // Get ghost DOFs
    CHKERR prb_mng_ptr->partitionGhostDofs(dm_field->problemName,
                                           dm_field->verbosity);
  }

  // Set flag that problem is build and partitioned
  dm_field->isProblemBuild = PETSC_TRUE;

  MoFEMFunctionReturn(0);
}

PetscErrorCode DMSubDMSetUp_MoFEM(DM subdm) {
  PetscValidHeaderSpecific(subdm, DM_CLASSID, 1);
  ProblemsManager *prb_mng_ptr;
  MoFEMFunctionBegin;

  DMCtx *subdm_field = static_cast<DMCtx *>(subdm->data);

  // build sub dm problem
  CHKERR subdm_field->mField_ptr->getInterface(prb_mng_ptr);

  map<std::string, std::pair<EntityType, EntityType>> *entity_map_row = nullptr;
  map<std::string, std::pair<EntityType, EntityType>> *entity_map_col = nullptr;

  if (subdm_field->mapTypeRow)
    entity_map_row = subdm_field->mapTypeRow.get();
  if (subdm_field->mapTypeCol)
    entity_map_row = subdm_field->mapTypeCol.get();

  CHKERR prb_mng_ptr->buildSubProblem(
      subdm_field->problemName, subdm_field->rowFields, subdm_field->colFields,
      subdm_field->problemMainOfSubPtr->getName(),
      subdm_field->isSquareMatrix == PETSC_TRUE, entity_map_row, entity_map_col,
      subdm_field->verbosity);

  // partition problem
  subdm_field->isPartitioned = subdm_field->isPartitioned;
  if (subdm_field->isPartitioned) {
    CHKERR prb_mng_ptr->partitionFiniteElements(subdm_field->problemName, true,
                                                0, subdm_field->sIze,
                                                subdm_field->verbosity);
    // set ghost nodes
    CHKERR prb_mng_ptr->partitionGhostDofsOnDistributedMesh(
        subdm_field->problemName, subdm_field->verbosity);
  } else {
    // partition finite elements
    CHKERR prb_mng_ptr->partitionFiniteElements(subdm_field->problemName, false,
                                                -1, -1, subdm_field->verbosity);
    // set ghost nodes
    CHKERR prb_mng_ptr->partitionGhostDofs(subdm_field->problemName,
                                           subdm_field->verbosity);
  }

  subdm_field->isProblemBuild = PETSC_TRUE;

  MoFEMFunctionReturn(0);
}

PetscErrorCode DMGlobalToLocalBegin_MoFEM(DM dm, Vec g, InsertMode mode,
                                          Vec l) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  CHKERR VecGhostUpdateBegin(g, INSERT_VALUES, SCATTER_FORWARD);
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMGlobalToLocalEnd_MoFEM(DM dm, Vec g, InsertMode mode, Vec l) {
  MoFEMFunctionBeginHot;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;

  CHKERR VecGhostUpdateEnd(g, INSERT_VALUES, SCATTER_FORWARD);

  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  int nb_dofs = dm_field->problemPtr->getNbLocalDofsRow();
  int nb_ghost = dm_field->problemPtr->getNbGhostDofsRow();

  double *array_loc, *array_glob;
  CHKERR VecGetArray(l, &array_loc);
  CHKERR VecGetArray(g, &array_glob);
  switch (mode) {
  case INSERT_VALUES:
    cblas_dcopy(nb_dofs + nb_ghost, array_glob, 1, array_loc, 1);
    break;
  case ADD_VALUES:
    cblas_daxpy(nb_dofs + nb_ghost, 1, array_glob, 1, array_loc, 1);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  CHKERR VecRestoreArray(l, &array_loc);
  CHKERR VecRestoreArray(g, &array_glob);
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMLocalToGlobalBegin_MoFEM(DM dm, Vec l, InsertMode mode,
                                          Vec g) {

  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;

  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  int nb_dofs = dm_field->problemPtr->getNbLocalDofsRow();
  int nb_ghost = dm_field->problemPtr->getNbGhostDofsRow();

  double *array_loc, *array_glob;
  CHKERR VecGetArray(l, &array_loc);
  CHKERR VecGetArray(g, &array_glob);
  switch (mode) {
  case INSERT_VALUES:
    cblas_dcopy(nb_dofs + nb_ghost, array_loc, 1, array_glob, 1);
    break;
  case ADD_VALUES:
    cblas_daxpy(nb_dofs + nb_ghost, 1, array_loc, 1, array_glob, 1);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  CHKERR VecRestoreArray(l, &array_loc);
  CHKERR VecRestoreArray(g, &array_glob);

  MoFEMFunctionReturn(0);
}

PetscErrorCode DMLocalToGlobalEnd_MoFEM(DM, Vec l, InsertMode mode, Vec g) {
  //
  MoFEMFunctionBeginHot;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode DMCreateFieldIS_MoFEM(DM dm, PetscInt *numFields,
                                     char ***fieldNames, IS **fields) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;

  if (numFields) {
    *numFields = 0;
  }
  if (fieldNames) {
    *fieldNames = NULL;
  }
  if (fields) {
    *fields = NULL;
  }

  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  auto fields_ptr = dm_field->mField_ptr->get_fields();
  Field_multiIndex::iterator fit, hi_fit;
  fit = fields_ptr->begin();
  hi_fit = fields_ptr->end();
  *numFields = std::distance(fit, hi_fit);

  if (fieldNames) {
    CHKERR PetscMalloc1(*numFields, fieldNames);
  }
  if (fields) {
    CHKERR PetscMalloc1(*numFields, fields);
  }

  for (int f = 0; fit != hi_fit; fit++, f++) {
    if (fieldNames) {
      CHKERR PetscStrallocpy(fit->get()->getName().c_str(),
                             (char **)&(*fieldNames)[f]);
    }
    if (fields) {
      CHKERR dm_field->mField_ptr->getInterface<ISManager>()
          ->isCreateProblemFieldAndRank(
              dm_field->problemPtr->getName(), ROW, fit->get()->getName(), 0,
              fit->get()->getNbOfCoeffs(), &(*fields)[f]);
    }
  }

  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMGetFieldIS(DM dm, RowColData rc, const char field_name[],
                                 IS *is) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  CHKERR dm_field->mField_ptr->getInterface<ISManager>()
      ->isCreateProblemFieldAndRank(dm_field->problemPtr->getName(), ROW,
                                    field_name, 0, 1000, is);
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMMoFEMSetVerbosity(DM dm, const int verb) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBeginHot;
  DMCtx *dm_field = static_cast<DMCtx *>(dm->data);
  dm_field->verbosity = verb;
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM