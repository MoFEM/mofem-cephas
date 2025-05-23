/** \file VecManager.cpp
 * \brief Implementation of VecManager
 */

namespace MoFEM {

MoFEMErrorCode
VecManager::query_interface(boost::typeindex::type_index type_index,
                            UnknownInterface **iface) const {
  *iface = const_cast<VecManager *>(this);
  return 0;
}

VecManager::VecManager(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)), dEbug(false) {
  if (!LogManager::checkIfChannelExist("VECWORLD")) {
    auto core_log = logging::core::get();

    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "VECWORLD"));
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSync(), "VECSYNC"));
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSelf(), "VECSELF"));

    LogManager::setLog("VECWORLD");
    LogManager::setLog("VECSYNC");
    LogManager::setLog("VECSELF");

    MOFEM_LOG_TAG("VECWORLD", "VecMng");
    MOFEM_LOG_TAG("VECSYNC", "VecMng");
    MOFEM_LOG_TAG("VECSELF", "VecMng");
  }

  MOFEM_LOG("VECWORLD", Sev::noisy) << "Vector manager created";
}

MoFEMErrorCode VecManager::vecCreateSeq(const std::string name, RowColData rc,
                                        Vec *V) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(name, &problem_ptr);
  DofIdx nb_local_dofs, nb_ghost_dofs;
  switch (rc) {
  case ROW:
    nb_local_dofs = problem_ptr->getNbLocalDofsRow();
    nb_ghost_dofs = problem_ptr->getNbGhostDofsRow();
    break;
  case COL:
    nb_local_dofs = problem_ptr->getNbLocalDofsCol();
    nb_ghost_dofs = problem_ptr->getNbGhostDofsCol();
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }
  CHKERR ::VecCreateSeq(PETSC_COMM_SELF, nb_local_dofs + nb_ghost_dofs, V);
  CHKERR ::VecSetOption(
      *V, VEC_IGNORE_NEGATIVE_INDICES,
      PETSC_TRUE); // As default in MOFEM we will assume that this is always on
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VecManager::vecCreateGhost(const std::string name, RowColData rc,
                                          Vec *V) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(name, &problem_ptr);
  DofIdx nb_dofs, nb_local_dofs, nb_ghost_dofs;
  NumeredDofEntityByLocalIdx *dofs;
  switch (rc) {
  case ROW:
    nb_dofs = problem_ptr->getNbDofsRow();
    nb_local_dofs = problem_ptr->getNbLocalDofsRow();
    nb_ghost_dofs = problem_ptr->getNbGhostDofsRow();
    dofs = const_cast<NumeredDofEntityByLocalIdx *>(
        &problem_ptr->numeredRowDofsPtr->get<PetscLocalIdx_mi_tag>());
    break;
  case COL:
    nb_dofs = problem_ptr->getNbDofsCol();
    nb_local_dofs = problem_ptr->getNbLocalDofsCol();
    nb_ghost_dofs = problem_ptr->getNbGhostDofsCol();
    dofs = const_cast<NumeredDofEntityByLocalIdx *>(
        &problem_ptr->numeredColDofsPtr->get<PetscLocalIdx_mi_tag>());
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  // get ghost dofs
  auto miit = dofs->lower_bound(nb_local_dofs);
  auto hi_miit = dofs->upper_bound(nb_local_dofs + nb_ghost_dofs);
  int count = std::distance(miit, hi_miit);
  if (count != nb_ghost_dofs) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  std::vector<DofIdx> ghost_idx;
  ghost_idx.reserve(count);
  for (; miit != hi_miit; ++miit) {
    if ((*miit)->getActive()) {
      ghost_idx.emplace_back((*miit)->petscGloablDofIdx);

#ifndef NDEBUG
      if (PetscUnlikely(ghost_idx.back() > nb_dofs || ghost_idx.back() < 0)) {
        MOFEM_LOG_FUNCTION();
        MOFEM_LOG_ATTRIBUTES("VECSELF",
                             LogManager::BitLineID | LogManager::BitScope);
        MOFEM_LOG("VECSELF", Sev::error) << "Problem with DOF " << **miit;
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Wrong ghost index");
      }
#endif
    }
  }
  CHKERR ::VecCreateGhost(m_field.get_comm(), nb_local_dofs, nb_dofs,
                          ghost_idx.size(), &*ghost_idx.begin(), V);
  CHKERR ::VecSetOption(
      *V, VEC_IGNORE_NEGATIVE_INDICES,
      PETSC_TRUE); // As default in MOFEM we will assume that this is always on
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VecManager::vecCreateGhost(const std::string name, RowColData rc,
                                          SmartPetscObj<Vec> &v_ptr) const {
  MoFEMFunctionBegin;
  Vec vec;
  CHKERR vecCreateGhost(name, rc, &vec);
  v_ptr.reset(vec, false);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VecManager::vecScatterCreate(
    Vec xin, const std::string x_problem, const std::string x_field_name,
    RowColData x_rc, Vec yin, const std::string y_problem,
    const std::string y_field_name, RowColData y_rc, VecScatter *newctx) const {
  const MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  std::vector<int> idx(0), idy(0);
  CHKERR m_field.getInterface<ISManager>()
      ->isCreateFromProblemFieldToOtherProblemField(
          x_problem, x_field_name, x_rc, y_problem, y_field_name, y_rc, idx,
          idy);
  IS ix, iy;
  CHKERR ISCreateGeneral(PETSC_COMM_WORLD, idx.size(), &idx[0],
                         PETSC_USE_POINTER, &ix);
  CHKERR ISCreateGeneral(PETSC_COMM_WORLD, idy.size(), &idy[0],
                         PETSC_USE_POINTER, &iy);
  CHKERR ::VecScatterCreate(xin, ix, yin, iy, newctx);
  CHKERR ISDestroy(&ix);
  CHKERR ISDestroy(&iy);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VecManager::vecScatterCreate(
    Vec xin, const std::string x_problem, RowColData x_rc, Vec yin,
    const std::string y_problem, RowColData y_rc, VecScatter *newctx) const {
  const MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  std::vector<int> idx(0), idy(0);
  CHKERR m_field.getInterface<ISManager>()->isCreateFromProblemToOtherProblem(
      x_problem, x_rc, y_problem, y_rc, idx, idy);
  IS ix, iy;
  CHKERR ISCreateGeneral(PETSC_COMM_WORLD, idx.size(), &idx[0],
                         PETSC_USE_POINTER, &ix);
  CHKERR ISCreateGeneral(PETSC_COMM_WORLD, idy.size(), &idy[0],
                         PETSC_USE_POINTER, &iy);
  if (dEbug) {
    ISView(ix, PETSC_VIEWER_STDOUT_WORLD);
    ISView(iy, PETSC_VIEWER_STDOUT_WORLD);
  }
  CHKERR ::VecScatterCreate(xin, ix, yin, iy, newctx);
  CHKERR ISDestroy(&ix);
  CHKERR ISDestroy(&iy);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VecManager::vecScatterCreate(Vec xin, const std::string x_problem,
                             const std::string x_field_name, RowColData x_rc,
                             Vec yin, const std::string y_problem,
                             const std::string y_field_name, RowColData y_rc,
                             SmartPetscObj<VecScatter> &smart_newctx) const {
  MoFEMFunctionBegin;
  VecScatter newctx;
  CHKERR vecScatterCreate(xin, x_problem, x_field_name, x_rc, yin, y_problem,
                          y_field_name, y_rc, &newctx);
  smart_newctx = SmartPetscObj<VecScatter>(newctx);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VecManager::vecScatterCreate(Vec xin, const std::string x_problem,
                             RowColData x_rc, Vec yin,
                             const std::string y_problem, RowColData y_rc,
                             SmartPetscObj<VecScatter> &smart_newctx) const {
  MoFEMFunctionBegin;
  VecScatter newctx;
  CHKERR vecScatterCreate(xin, x_problem, x_rc, yin, y_problem, y_rc, &newctx);
  smart_newctx = SmartPetscObj<VecScatter>(newctx);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VecManager::setLocalGhostVector(const Problem *problem_ptr,
                                               RowColData rc, Vec V,
                                               InsertMode mode,
                                               ScatterMode scatter_mode) const {
  MoFEMFunctionBegin;
  NumeredDofEntityByLocalIdx *dofs;
  DofIdx nb_local_dofs, nb_ghost_dofs;
  switch (rc) {
  case ROW:
    nb_local_dofs = problem_ptr->getNbLocalDofsRow();
    nb_ghost_dofs = problem_ptr->getNbGhostDofsRow();
    dofs = const_cast<NumeredDofEntityByLocalIdx *>(
        &problem_ptr->numeredRowDofsPtr->get<PetscLocalIdx_mi_tag>());
    break;
  case COL:
    nb_local_dofs = problem_ptr->getNbLocalDofsCol();
    nb_ghost_dofs = problem_ptr->getNbGhostDofsCol();
    dofs = const_cast<NumeredDofEntityByLocalIdx *>(
        &problem_ptr->numeredColDofsPtr->get<PetscLocalIdx_mi_tag>());
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  Vec Vlocal;
  CHKERR VecGhostGetLocalForm(V, &Vlocal);
  int size;
  CHKERR VecGetLocalSize(V, &size);
  if (size != nb_local_dofs) {
    SETERRQ3(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "data inconsistency: check ghost vector, problem <%s> with nb. of "
             "local nodes %d != %d",
             problem_ptr->getName().c_str(), size, nb_local_dofs);
  }
  CHKERR VecGetLocalSize(Vlocal, &size);
  if (size != nb_local_dofs + nb_ghost_dofs) {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "data inconsistency: check ghost vector, problem with nb. of "
             "ghost nodes %d != ",
             size, nb_local_dofs + nb_ghost_dofs);
  }
  auto miit = dofs->lower_bound(0);
  auto hi_miit = dofs->upper_bound(nb_local_dofs + nb_ghost_dofs);
  DofIdx ii = 0;
  switch (scatter_mode) {
  case SCATTER_FORWARD: {
    PetscScalar *array;
    CHKERR VecGetArray(Vlocal, &array);
    switch (mode) {
    case INSERT_VALUES:
      for (; miit != hi_miit; ++miit, ++ii)
        array[ii] = (*miit)->getFieldData();
      break;
    case ADD_VALUES:
      for (; miit != hi_miit; ++miit, ++ii)
        array[ii] += (*miit)->getFieldData();
      break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
    }
    CHKERR VecRestoreArray(Vlocal, &array);
  }; break;
  case SCATTER_REVERSE: {
    const PetscScalar *array;
    CHKERR VecGetArrayRead(Vlocal, &array);
    switch (mode) {
    case INSERT_VALUES:
      for (; miit != hi_miit; ++miit, ++ii) {
        (*miit)->getFieldData() = array[ii];
      }
      break;
    case ADD_VALUES:
      for (; miit != hi_miit; ++miit, ++ii)
        (*miit)->getFieldData() += array[ii];
      break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
    }
    CHKERR VecRestoreArrayRead(Vlocal, &array);
  }; break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  CHKERR VecGhostRestoreLocalForm(V, &Vlocal);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VecManager::setLocalGhostVector(const std::string name,
                                               RowColData rc, Vec V,
                                               InsertMode mode,
                                               ScatterMode scatter_mode) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(name, &problem_ptr);
  CHKERR setLocalGhostVector(problem_ptr, rc, V, mode, scatter_mode);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VecManager::setGlobalGhostVector(const Problem *problem_ptr, RowColData rc,
                                 Vec V, InsertMode mode,
                                 ScatterMode scatter_mode) const {
  MoFEMFunctionBegin;
  typedef NumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type
      DofsByGlobalIdx;
  DofsByGlobalIdx *dofs;
  DofIdx nb_dofs;
  switch (rc) {
  case ROW:
    nb_dofs = problem_ptr->getNbDofsRow();
    dofs = const_cast<DofsByGlobalIdx *>(
        &problem_ptr->numeredRowDofsPtr->get<PetscGlobalIdx_mi_tag>());
    break;
  case COL:
    nb_dofs = problem_ptr->getNbDofsCol();
    dofs = const_cast<DofsByGlobalIdx *>(
        &problem_ptr->numeredColDofsPtr->get<PetscGlobalIdx_mi_tag>());
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Function works only for ROWs and COLs");
  }
  DofsByGlobalIdx::iterator miit = dofs->lower_bound(0);
  DofsByGlobalIdx::iterator hi_miit = dofs->upper_bound(nb_dofs);
  auto comm = PetscObjectComm((PetscObject)V);
  switch (scatter_mode) {
  case SCATTER_REVERSE: {
    VecScatter ctx;
    Vec V_glob;
    CHKERR VecScatterCreateToAll(V, &ctx, &V_glob);
    CHKERR VecScatterBegin(ctx, V, V_glob, INSERT_VALUES, SCATTER_FORWARD);
    CHKERR VecScatterEnd(ctx, V, V_glob, INSERT_VALUES, SCATTER_FORWARD);
    int size;
    CHKERR VecGetSize(V_glob, &size);
    if (size != nb_dofs) {
      SETERRQ(comm, MOFEM_DATA_INCONSISTENCY,
              "Size of vector is inconsistent with size of problem. You could "
              "use wrong vector with wrong problem, or you created vector "
              "before you remove DOFs from problem.");
    }
    PetscScalar *array;
    CHKERR VecGetArray(V_glob, &array);
    switch (mode) {
    case INSERT_VALUES:
      for (; miit != hi_miit; miit++)
        (*miit)->getFieldData() = array[(*miit)->getPetscGlobalDofIdx()];
      break;
    case ADD_VALUES:
      for (; miit != hi_miit; miit++)
        (*miit)->getFieldData() += array[(*miit)->getPetscGlobalDofIdx()];
      break;
    default:
      SETERRQ(comm, MOFEM_NOT_IMPLEMENTED, "not implemented");
    }
    CHKERR VecRestoreArray(V_glob, &array);
    CHKERR VecScatterDestroy(&ctx);
    CHKERR VecDestroy(&V_glob);
    break;
  }
  default:
    SETERRQ(comm, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VecManager::setGlobalGhostVector(const std::string name, RowColData rc, Vec V,
                                 InsertMode mode,
                                 ScatterMode scatter_mode) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(name, &problem_ptr);
  CHKERR setGlobalGhostVector(problem_ptr, rc, V, mode, scatter_mode);
  MoFEMFunctionReturn(0);
}

template <int MODE> struct SetOtherLocalGhostVector {
  template <typename A0, typename A1, typename A2, typename A3, typename A4>
  inline MoFEMErrorCode operator()(A0 dofs_ptr, A1 array, A2 miit, A3 hi_miit,
                                   A4 &cpy_bit_number,
                                   const std::string cpy_field_name) {
    MoFEMFunctionBegin;
    for (; miit != hi_miit; miit++) {
      if (miit->get()->getHasLocalIndex()) {
        const auto uid = DofEntity::getUniqueIdCalculate(
            (*miit)->getEntDofIdx(), FieldEntity::getLocalUniqueIdCalculate(
                                         cpy_bit_number, (*miit)->getEnt()));
        auto diiiit = dofs_ptr->template get<Unique_mi_tag>().find(uid);
        if (diiiit == dofs_ptr->template get<Unique_mi_tag>().end()) {
          MOFEM_LOG_FUNCTION();
          MOFEM_LOG_ATTRIBUTES("VECSELF",
                               LogManager::BitLineID | LogManager::BitScope);
          MOFEM_LOG("VECSELF", Sev::error)
              << "Problem finding DOFs in the copy field";
          MOFEM_LOG("VECSELF", Sev::error) << "Field DOF: " << (**miit);
          MOFEM_LOG("VECSELF", Sev::error)
              << "Copy field name: " << cpy_field_name;
          SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                  "Automatic creation of entity and dof not implemented");
        }
        if constexpr (MODE == INSERT_VALUES)
          (*diiiit)->getFieldData() = array[(*miit)->getPetscLocalDofIdx()];
        else if constexpr (MODE == ADD_VALUES)
          (*diiiit)->getFieldData() += array[(*miit)->getPetscLocalDofIdx()];
      }
    }
    MoFEMFunctionReturn(0);
  }
};

MoFEMErrorCode VecManager::setOtherLocalGhostVector(
    const Problem *problem_ptr, const std::string field_name,
    const std::string cpy_field_name, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode) const {
  const MoFEM::Interface &m_field = cOre;
  auto fields_ptr = m_field.get_fields();
  auto dofs_ptr = m_field.get_dofs();
  MoFEMFunctionBegin;
  using DofsByUId = NumeredDofEntity_multiIndex::index<Unique_mi_tag>::type;
  DofsByUId *dofs;

  switch (rc) {
  case ROW:
    dofs = const_cast<DofsByUId *>(
        &problem_ptr->numeredRowDofsPtr->get<Unique_mi_tag>());
    break;
  case COL:
    dofs = const_cast<DofsByUId *>(
        &problem_ptr->numeredColDofsPtr->get<Unique_mi_tag>());
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }

  auto cpy_fit = fields_ptr->get<FieldName_mi_tag>().find(cpy_field_name);
  if (cpy_fit == fields_ptr->get<FieldName_mi_tag>().end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
             "cpy field < %s > not found, (top tip: check spelling)",
             cpy_field_name.c_str());
  const auto cpy_bit_number = (*cpy_fit)->getBitNumber();

  auto miit = dofs->lower_bound(
      FieldEntity::getLoBitNumberUId(m_field.get_field_bit_number(field_name)));
  if (miit == dofs->end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
             "cpy field < %s > not found, (top tip: check spelling)",
             field_name.c_str());
  }
  auto hi_miit = dofs->upper_bound(
      FieldEntity::getHiBitNumberUId(m_field.get_field_bit_number(field_name)));

  if ((*miit)->getSpace() != (*cpy_fit)->getSpace()) {
    SETERRQ4(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "fields have to have same space (%s) %s != (%s) %s",
             (*miit)->getName().c_str(), FieldSpaceNames[(*miit)->getSpace()],
             cpy_field_name.c_str(), FieldSpaceNames[(*cpy_fit)->getSpace()]);
  }
  if ((*miit)->getNbOfCoeffs() != (*cpy_fit)->getNbOfCoeffs()) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "fields have to have same rank");
  }

  switch (scatter_mode) {
  case SCATTER_REVERSE: {

    PetscScalar *array;
    CHKERR VecGetArray(V, &array);
    if (mode == INSERT_VALUES)
      CHKERR SetOtherLocalGhostVector<INSERT_VALUES>()(
          dofs_ptr, array, miit, hi_miit, cpy_bit_number, cpy_field_name);
    else if (mode == ADD_VALUES)
      CHKERR SetOtherLocalGhostVector<ADD_VALUES>()(
          dofs_ptr, array, miit, hi_miit, cpy_bit_number, cpy_field_name);
    else
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
              "Operation mode not implemented");
    CHKERR VecRestoreArray(V, &array);
  } break;
  case SCATTER_FORWARD: {
    for (; miit != hi_miit; miit++) {
      if (!miit->get()->getHasLocalIndex())
        continue;
      const auto uid = DofEntity::getUniqueIdCalculate(
          (*miit)->getEntDofIdx(), FieldEntity::getLocalUniqueIdCalculate(
                                       cpy_bit_number, (*miit)->getEnt()));
      auto diiiit = dofs_ptr->get<Unique_mi_tag>().find(uid);
      if (diiiit == dofs_ptr->get<Unique_mi_tag>().end()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "no data to fill the vector (top tip: you want scatter forward "
                "or scatter reverse?)");
      }
      CHKERR VecSetValue(V, (*miit)->getPetscGlobalDofIdx(),
                         (*diiiit)->getFieldData(), mode);
    }
    CHKERR VecAssemblyBegin(V);
    CHKERR VecAssemblyEnd(V);
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "not implemented");
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VecManager::setOtherLocalGhostVector(
    const std::string name, const std::string field_name,
    const std::string cpy_field_name, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(name, &problem_ptr);
  CHKERR setOtherLocalGhostVector(problem_ptr, field_name, cpy_field_name, rc,
                                  V, mode, scatter_mode);
  MoFEMFunctionReturn(0);
}

template <int MODE> struct SetOtherGlobalGhostVector {
  template <typename A0, typename A1, typename A2, typename A3, typename A4>
  inline MoFEMErrorCode operator()(A0 dofs_ptr, A1 array, A2 miit, A3 hi_miit,
                                   A4 &cpy_bit_number) {
    MoFEMFunctionBegin;
    for (; miit != hi_miit; miit++) {
      const auto uid = DofEntity::getUniqueIdCalculate(
          (*miit)->getEntDofIdx(), FieldEntity::getLocalUniqueIdCalculate(
                                       cpy_bit_number, (*miit)->getEnt()));
      auto diiiit = dofs_ptr->template get<Unique_mi_tag>().find(uid);
#ifndef NDEBUG
      if (diiiit == dofs_ptr->template get<Unique_mi_tag>().end()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                "Automatic creation of entity and dof not implemented");
      }
#endif
      if constexpr (MODE == INSERT_VALUES)
        (*diiiit)->getFieldData() = array[(*miit)->getPetscGlobalDofIdx()];
      else if constexpr (MODE == ADD_VALUES)
        (*diiiit)->getFieldData() += array[(*miit)->getPetscGlobalDofIdx()];
    }
    MoFEMFunctionReturn(0);
  }
};

MoFEMErrorCode VecManager::setOtherGlobalGhostVector(
    const Problem *problem_ptr, const std::string field_name,
    const std::string cpy_field_name, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode) const {
  const MoFEM::Interface &m_field = cOre;
  auto fields_ptr = m_field.get_fields();
  auto dofs_ptr = m_field.get_dofs();
  MoFEMFunctionBegin;
  NumeredDofEntityByUId *dofs;
  DofIdx nb_dofs;
  switch (rc) {
  case ROW:
    nb_dofs = problem_ptr->getNbDofsRow();
    dofs = const_cast<NumeredDofEntityByUId *>(
        &problem_ptr->numeredRowDofsPtr->get<Unique_mi_tag>());
    break;
  case COL:
    nb_dofs = problem_ptr->getNbDofsCol();
    dofs = const_cast<NumeredDofEntityByUId *>(
        &problem_ptr->numeredColDofsPtr->get<Unique_mi_tag>());
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "not implemented");
  }
  auto cpy_fit = fields_ptr->get<FieldName_mi_tag>().find(cpy_field_name);
  if (cpy_fit == fields_ptr->get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
             "cpy field < %s > not found, (top tip: check spelling)",
             cpy_field_name.c_str());
  }
  const auto cpy_bit_number = (*cpy_fit)->getBitNumber();

  auto miit = dofs->lower_bound(
      FieldEntity::getLoBitNumberUId(m_field.get_field_bit_number(field_name)));
  if (miit == dofs->end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
             "problem field < %s > not found, (top tip: check spelling)",
             field_name.c_str());
  }
  auto hi_miit = dofs->upper_bound(
      FieldEntity::getHiBitNumberUId(m_field.get_field_bit_number(field_name)));
  if ((*miit)->getSpace() != (*cpy_fit)->getSpace()) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "fields have to have same space");
  }
  if ((*miit)->getNbOfCoeffs() != (*cpy_fit)->getNbOfCoeffs()) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "fields have to have same rank");
  }
  switch (scatter_mode) {
  case SCATTER_REVERSE: {
    Vec V_glob;
    VecScatter ctx;
    CHKERR VecScatterCreateToAll(V, &ctx, &V_glob);
    CHKERR VecScatterBegin(ctx, V, V_glob, INSERT_VALUES, SCATTER_FORWARD);
    CHKERR VecScatterEnd(ctx, V, V_glob, INSERT_VALUES, SCATTER_FORWARD);
    int size;
    CHKERR VecGetSize(V_glob, &size);
    if (size != nb_dofs)
      SETERRQ(
          PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
          "data inconsistency: nb. of dofs and declared nb. dofs in database");
    PetscScalar *array;
    CHKERR VecGetArray(V_glob, &array);

    if (mode == INSERT_VALUES)
      CHKERR SetOtherGlobalGhostVector<INSERT_VALUES>()(
          dofs_ptr, array, miit, hi_miit, cpy_bit_number);
    else if (mode == ADD_VALUES)
      CHKERR SetOtherGlobalGhostVector<ADD_VALUES>()(dofs_ptr, array, miit,
                                                     hi_miit, cpy_bit_number);
    else
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
              "Operation mode not implemented");

    CHKERR VecRestoreArray(V_glob, &array);
    CHKERR VecDestroy(&V_glob);
    CHKERR VecScatterDestroy(&ctx);
  } break;
  case SCATTER_FORWARD: {
    for (; miit != hi_miit; miit++) {
      const auto uid = DofEntity::getUniqueIdCalculate(
          (*miit)->getEntDofIdx(), FieldEntity::getLocalUniqueIdCalculate(
                                       cpy_bit_number, (*miit)->getEnt()));
      auto diiiit = dofs_ptr->get<Unique_mi_tag>().find(uid);
      if (diiiit == dofs_ptr->get<Unique_mi_tag>().end()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "No data to fill the vector (top tip: you want scatter forward "
                "or scatter reverse?)");
      }
      CHKERR VecSetValue(V, (*miit)->getPetscGlobalDofIdx(),
                         (*diiiit)->getFieldData(), mode);
    }
    CHKERR VecAssemblyBegin(V);
    CHKERR VecAssemblyEnd(V);
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VecManager::setOtherGlobalGhostVector(
    const std::string name, const std::string field_name,
    const std::string cpy_field_name, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(name, &problem_ptr);
  CHKERR setOtherGlobalGhostVector(problem_ptr, field_name, cpy_field_name, rc,
                                   V, mode, scatter_mode);
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
