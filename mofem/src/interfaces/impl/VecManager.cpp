/** \file VecManager.cpp
 * \brief Implementation of VecManager
 *
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

namespace MoFEM {

MoFEMErrorCode VecManager::query_interface(const MOFEMuuid &uuid,
                                           UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMVEC) {
    *iface = const_cast<VecManager *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

VecManager::VecManager(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)), dEbug(false) {}
VecManager::~VecManager() {}

MoFEMErrorCode VecManager::vecCreateSeq(const std::string &name, RowColData rc,
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
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VecManager::vecCreateGhost(const std::string &name,
                                          RowColData rc, Vec *V) const {
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
        &problem_ptr->numeredDofsRows->get<PetscLocalIdx_mi_tag>());
    break;
  case COL:
    nb_dofs = problem_ptr->getNbDofsCol();
    nb_local_dofs = problem_ptr->getNbLocalDofsCol();
    nb_ghost_dofs = problem_ptr->getNbGhostDofsCol();
    dofs = const_cast<NumeredDofEntityByLocalIdx *>(
        &problem_ptr->numeredDofsCols->get<PetscLocalIdx_mi_tag>());
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
  std::vector<DofIdx> ghost_idx(count);
  for (auto vit = ghost_idx.begin(); miit != hi_miit; ++miit, ++vit) {
    *vit = (*miit)->petscGloablDofIdx;
  }
  CHKERR ::VecCreateGhost(PETSC_COMM_WORLD, nb_local_dofs, nb_dofs,
                          nb_ghost_dofs, &ghost_idx[0], V);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VecManager::vecScatterCreate(Vec xin, const std::string &x_problem,
                             const std::string &x_field_name, RowColData x_rc,
                             Vec yin, const std::string &y_problem,
                             const std::string &y_field_name, RowColData y_rc,
                             VecScatter *newctx) const {
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
    Vec xin, const std::string &x_problem, RowColData x_rc, Vec yin,
    const std::string &y_problem, RowColData y_rc, VecScatter *newctx) const {
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
        &problem_ptr->numeredDofsRows->get<PetscLocalIdx_mi_tag>());
    break;
  case COL:
    nb_local_dofs = problem_ptr->getNbLocalDofsCol();
    nb_ghost_dofs = problem_ptr->getNbGhostDofsCol();
    dofs = const_cast<NumeredDofEntityByLocalIdx *>(
        &problem_ptr->numeredDofsCols->get<PetscLocalIdx_mi_tag>());
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  Vec Vlocal;
  CHKERR VecGhostGetLocalForm(V, &Vlocal);
  int size;
  CHKERR VecGetLocalSize(V, &size);
  if (size != nb_local_dofs) {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "data inconsistency: check ghost vector, problem with nb. of "
             "local nodes %d != %d",
             size, nb_local_dofs);
  }
  CHKERR VecGetLocalSize(Vlocal, &size);
  if (size != nb_local_dofs + nb_ghost_dofs) {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "data inconsistency: check ghost vector, problem with nb. of "
             "ghost nodes %d != ",
             size, nb_local_dofs + nb_ghost_dofs);
  }
  NumeredDofEntityByLocalIdx::iterator miit = dofs->lower_bound(0);
  NumeredDofEntityByLocalIdx::iterator hi_miit =
      dofs->upper_bound(nb_local_dofs + nb_ghost_dofs);
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
        // std::cerr << *miit << std::endl;
        // std::cerr << array[ii] << std::endl;
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

MoFEMErrorCode VecManager::setLocalGhostVector(const std::string &name,
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
        &problem_ptr->numeredDofsRows->get<PetscGlobalIdx_mi_tag>());
    break;
  case COL:
    nb_dofs = problem_ptr->getNbDofsCol();
    dofs = const_cast<DofsByGlobalIdx *>(
        &problem_ptr->numeredDofsCols->get<PetscGlobalIdx_mi_tag>());
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  DofsByGlobalIdx::iterator miit = dofs->lower_bound(0);
  DofsByGlobalIdx::iterator hi_miit = dofs->upper_bound(nb_dofs);
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
      SETERRQ(
          PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
          "data inconsistency: nb. of dofs and declared nb. dofs in database");
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
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
    }
    CHKERR VecRestoreArray(V_glob, &array);
    CHKERR VecScatterDestroy(&ctx);
    CHKERR VecDestroy(&V_glob);
    break;
  }
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
VecManager::setGlobalGhostVector(const std::string &name, RowColData rc, Vec V,
                                 InsertMode mode,
                                 ScatterMode scatter_mode) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(name, &problem_ptr);
  CHKERR setGlobalGhostVector(problem_ptr, rc, V, mode, scatter_mode);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VecManager::setOtherLocalGhostVector(
    const Problem *problem_ptr, const std::string &field_name,
    const std::string &cpy_field_name, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode) const {
  const MoFEM::Interface &m_field = cOre;
  const Field_multiIndex *fields_ptr;
  const DofEntity_multiIndex *dofs_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_fields(&fields_ptr);
  CHKERR m_field.get_dofs(&dofs_ptr);
  typedef NumeredDofEntity_multiIndex::index<FieldName_mi_tag>::type DofsByName;
  DofsByName *dofs;
  switch (rc) {
  case ROW:
    dofs = const_cast<DofsByName *>(
        &problem_ptr->numeredDofsRows->get<FieldName_mi_tag>());
    break;
  case COL:
    dofs = const_cast<DofsByName *>(
        &problem_ptr->numeredDofsCols->get<FieldName_mi_tag>());
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator cpy_fit =
      fields_ptr->get<FieldName_mi_tag>().find(cpy_field_name);
  if (cpy_fit == fields_ptr->get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
             "cpy field < %s > not found, (top tip: check spelling)",
             cpy_field_name.c_str());
  }
  DofsByName::iterator miit = dofs->lower_bound(field_name);
  if (miit == dofs->end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
             "cpy field < %s > not found, (top tip: check spelling)",
             field_name.c_str());
  }
  DofsByName::iterator hi_miit = dofs->upper_bound(field_name);
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
    bool alpha = true;
    switch (mode) {
    case INSERT_VALUES:
      break;
    case ADD_VALUES:
      alpha = false;
      break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "not implemented");
    }
    PetscScalar *array;
    VecGetArray(V, &array);
    for (; miit != hi_miit; miit++) {
      if (!miit->get()->getHasLocalIndex())
        continue;
      DofEntity_multiIndex::index<
          Composite_Name_And_Ent_And_EntDofIdx_mi_tag>::type::iterator diiiit;
      diiiit =
          dofs_ptr->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().find(
              boost::make_tuple(cpy_field_name, (*miit)->getEnt(),
                                (*miit)->getEntDofIdx()));
      if (diiiit ==
          dofs_ptr->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().end()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
                "equivalent dof does not exist, use "
                "set_other_global_ghost_vector to create dofs entries");
      }
      if (alpha) {
        (*diiiit)->getFieldData() = array[(*miit)->getPetscLocalDofIdx()];
      } else {
        (*diiiit)->getFieldData() += array[(*miit)->getPetscLocalDofIdx()];
      }
    }
    CHKERR VecRestoreArray(V, &array);
  } break;
  case SCATTER_FORWARD: {
    for (; miit != hi_miit; miit++) {
      if (!miit->get()->getHasLocalIndex())
        continue;
      DofEntity_multiIndex::index<
          Composite_Name_And_Ent_And_EntDofIdx_mi_tag>::type::iterator diiiit;
      diiiit =
          dofs_ptr->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().find(
              boost::make_tuple(cpy_field_name, (*miit)->getEnt(),
                                (*miit)->getEntDofIdx()));
      if (diiiit ==
          dofs_ptr->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().end()) {
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
    const std::string &name, const std::string &field_name,
    const std::string &cpy_field_name, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(name, &problem_ptr);
  CHKERR setOtherLocalGhostVector(problem_ptr, field_name, cpy_field_name, rc,
                                  V, mode, scatter_mode);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VecManager::setOtherGlobalGhostVector(
    const Problem *problem_ptr, const std::string &field_name,
    const std::string &cpy_field_name, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode) const {
  const MoFEM::Interface &m_field = cOre;
  const Field_multiIndex *fields_ptr;
  const DofEntity_multiIndex *dofs_ptr;
  const FieldEntity_multiIndex *field_ents;
  MoFEMFunctionBegin;
  CHKERR m_field.get_fields(&fields_ptr);
  CHKERR m_field.get_dofs(&dofs_ptr);
  CHKERR m_field.get_field_ents(&field_ents);
  typedef NumeredDofEntityByFieldName DofsByName;
  DofsByName *dofs;
  DofIdx nb_dofs;
  switch (rc) {
  case ROW:
    nb_dofs = problem_ptr->getNbDofsRow();
    dofs = const_cast<DofsByName *>(
        &problem_ptr->numeredDofsRows->get<FieldName_mi_tag>());
    break;
  case COL:
    nb_dofs = problem_ptr->getNbDofsCol();
    dofs = const_cast<DofsByName *>(
        &problem_ptr->numeredDofsCols->get<FieldName_mi_tag>());
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "not implemented");
  }
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator cpy_fit =
      fields_ptr->get<FieldName_mi_tag>().find(cpy_field_name);
  if (cpy_fit == fields_ptr->get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
             "cpy field < %s > not found, (top tip: check spelling)",
             cpy_field_name.c_str());
  }
  DofsByName::iterator miit = dofs->lower_bound(field_name);
  if (miit == dofs->end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
             "problem field < %s > not found, (top tip: check spelling)",
             field_name.c_str());
  }
  DofsByName::iterator hi_miit = dofs->upper_bound(field_name);
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
    VecGetArray(V_glob, &array);
    bool alpha = true;
    switch (mode) {
    case INSERT_VALUES:
      break;
    case ADD_VALUES:
      alpha = false;
      break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "not implemented");
    }
    for (; miit != hi_miit; miit++) {
      if ((*miit)->getPetscGlobalDofIdx() >= size) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency: nb. of dofs and declared nb. dofs in "
                "database");
      }
      DofEntity_multiIndex::index<
          Composite_Name_And_Ent_And_EntDofIdx_mi_tag>::type::iterator diiiit;
      diiiit =
          dofs_ptr->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().find(
              boost::make_tuple(cpy_field_name, (*miit)->getEnt(),
                                (*miit)->getEntDofIdx()));
      if (diiiit ==
          dofs_ptr->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().end()) {
        EntityHandle ent = (*miit)->getEnt();
        rval = const_cast<moab::Interface &>(m_field.get_moab())
                   .add_entities((*cpy_fit)->getMeshset(), &ent, 1);
        CHKERRQ_MOAB(rval);
        // create field moabent
        ApproximationOrder order = (*miit)->getMaxOrder();
        std::pair<FieldEntity_multiIndex::iterator, bool> p_e_miit;
        try {
          boost::shared_ptr<FieldEntity> moabent(
              new FieldEntity(*cpy_fit, (*miit)->getRefEntityPtr()));
          p_e_miit =
              const_cast<FieldEntity_multiIndex *>(field_ents)->insert(moabent);
        } catch (const std::exception &ex) {
          std::ostringstream ss;
          ss << "throw in method: " << ex.what() << " at line " << __LINE__
             << " in file " << __FILE__ << std::endl;
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, ss.str().c_str());
        }
        if ((*p_e_miit.first)->getMaxOrder() < order) {
          bool success =
              const_cast<FieldEntity_multiIndex *>(field_ents)
                  ->modify(p_e_miit.first, FieldEntity_change_order(order));
          if (!success)
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
        }
        // create field moabdof
        DofEntityByNameAndEnt::iterator hi_diit, diit;
        diit = dofs_ptr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(
            boost::make_tuple(field_name, (*miit)->getEnt()));
        hi_diit = dofs_ptr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(
            boost::make_tuple(field_name, (*miit)->getEnt()));
        for (; diit != hi_diit; diit++) {
          boost::shared_ptr<DofEntity> mdof =
              boost::shared_ptr<DofEntity>(new DofEntity(
                  *(p_e_miit.first), (*diit)->getDofOrder(),
                  (*diit)->getDofCoeffIdx(), (*diit)->getEntDofIdx()));
          std::pair<DofEntity_multiIndex::iterator, bool> cpy_p_diit;
          cpy_p_diit =
              const_cast<DofEntity_multiIndex *>(dofs_ptr)->insert(mdof);
          if (cpy_p_diit.second) {
            bool success = const_cast<DofEntity_multiIndex *>(dofs_ptr)->modify(
                cpy_p_diit.first, DofEntity_active_change(true));
            if (!success)
              SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                      "modification unsuccessful");
          }
        }
        diiiit =
            dofs_ptr->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().find(
                boost::make_tuple(cpy_field_name, (*miit)->getEnt(),
                                  (*miit)->getEntDofIdx()));
        if (diiiit ==
            dofs_ptr->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().end())
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "data inconsistency");
      }
      if (alpha)
        (*diiiit)->getFieldData() = 0;
      (*diiiit)->getFieldData() += array[(*miit)->getPetscGlobalDofIdx()];
      if (dEbug) {
        std::ostringstream ss;
        ss << *(*diiiit) << "set " << array[(*miit)->getPetscGlobalDofIdx()]
           << std::endl;
        PetscPrintf(PETSC_COMM_WORLD, ss.str().c_str());
      }
    }
    CHKERR VecRestoreArray(V_glob, &array);
    CHKERR VecDestroy(&V_glob);
    CHKERR VecScatterDestroy(&ctx);
  } break;
  case SCATTER_FORWARD: {
    for (; miit != hi_miit; miit++) {
      DofEntity_multiIndex::index<
          Composite_Name_And_Ent_And_EntDofIdx_mi_tag>::type::iterator diiiit;
      diiiit =
          dofs_ptr->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().find(
              boost::make_tuple(cpy_field_name, (*miit)->getEnt(),
                                (*miit)->getEntDofIdx()));
      if (diiiit ==
          dofs_ptr->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().end()) {
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
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VecManager::setOtherGlobalGhostVector(
    const std::string &name, const std::string &field_name,
    const std::string &cpy_field_name, RowColData rc, Vec V, InsertMode mode,
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
