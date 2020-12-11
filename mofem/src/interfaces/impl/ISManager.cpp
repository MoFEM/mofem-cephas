/** \file ISManager.cpp
 * \brief IS creating
 * \ingroup mofem_is_managers
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
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

MoFEMErrorCode ISManager::query_interface(const MOFEMuuid &uuid,
                                          UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMISManager) {
    *iface = const_cast<ISManager *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

ISManager::ISManager(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)), dEbug(false) {}
ISManager::~ISManager() {}

MoFEMErrorCode ISManager::sectionCreate(const std::string &problem_name,
                                        PetscSection *s,
                                        const RowColData row_col) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  auto fields_ptr = m_field.get_fields();
  auto fe_ptr = m_field.get_finite_elements();
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(problem_name, &problem_ptr);
  boost::shared_ptr<NumeredDofEntity_multiIndex> dofs;
  BitFieldId fields_ids;
  switch (row_col) {
  case ROW:
    dofs = problem_ptr->numeredRowDofsPtr;
    for (FiniteElement_multiIndex::iterator fit = fe_ptr->begin();
         fit != fe_ptr->end(); fit++) {
      if ((fit->get()->getId() & problem_ptr->getBitFEId()).any()) {
        fields_ids |= fit->get()->getBitFieldIdRow();
      }
    }
    break;
  case COL:
    dofs = problem_ptr->numeredColDofsPtr;
    for (FiniteElement_multiIndex::iterator fit = fe_ptr->begin();
         fit != fe_ptr->end(); fit++) {
      if ((fit->get()->getId() & problem_ptr->getBitFEId()).any()) {
        fields_ids |= fit->get()->getBitFieldIdCol();
      }
    }
    break;
  default:
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
            "Has to be ROW or COLUMN");
  }
  // get fields names on the problem
  map<std::string, std::pair<int, int>> fields_map;
  {
    int field = 0;
    for (Field_multiIndex::iterator fit = fields_ptr->begin();
         fit != fields_ptr->end(); fit++) {
      if ((fit->get()->getId() & fields_ids).any()) {
        fields_map[fit->get()->getName()].first = field++;
        fields_map[fit->get()->getName()].second = fit->get()->getNbOfCoeffs();
      }
    }
  }
  const int proc = m_field.get_comm_rank();
  CHKERR PetscSectionCreate(PETSC_COMM_WORLD, s);
  CHKERR PetscSectionSetNumFields(*s, fields_map.size());
  for (map<std::string, std::pair<int, int>>::iterator mit = fields_map.begin();
       mit != fields_map.end(); mit++) {
    CHKERR PetscSectionSetFieldName(*s, mit->second.first, mit->first.c_str());
    CHKERR PetscSectionSetFieldComponents(*s, mit->second.first,
                                          mit->second.second);
  }
  // determine number of points
  int nb_charts = 0;
  {
    auto dit = dofs->begin();
    auto hi_dit = dofs->end();
    for (; dit != hi_dit;) {
      if (static_cast<int>(dit->get()->getPart()) == proc) {
        const auto &ent_uid = dit->get()->getEntLocalUniqueId();
        while (dit != hi_dit && dit->get()->getEntLocalUniqueId() == ent_uid) {
          ++dit;
        }
        ++nb_charts;
      } else {
        ++dit;
      }
    }
  }
  // get layout, i.e. chart
  PetscLayout layout;
  CHKERR PetscLayoutCreate(PETSC_COMM_WORLD, &layout);
  CHKERR PetscLayoutSetBlockSize(layout, 1);
  CHKERR PetscLayoutSetLocalSize(layout, nb_charts);
  CHKERR PetscLayoutSetUp(layout);
  int rstart, rend;
  CHKERR PetscLayoutGetRange(layout, &rstart, &rend);
  CHKERR PetscLayoutDestroy(&layout);
  CHKERR PetscSectionSetChart(*s, rstart, rend);

  // loop of all dofs
  {
    auto dit = dofs->begin();
    auto hi_dit = dofs->end();
    int point = rstart;
    for (; dit != hi_dit;) {
      if (static_cast<int>(dit->get()->getPart()) == proc) {

        const auto &field_name = dit->get()->getName();

        int dd = 0;
        const auto &ent_uid = dit->get()->getEntLocalUniqueId();
        while (dit != hi_dit && dit->get()->getEntLocalUniqueId() == ent_uid) {
          ++dd;
          ++dit;
        }

        if (fields_map.find(field_name) == fields_map.end()) {
          MOFEM_LOG_C("SELF", Sev::warning, "Warning: Field %s not found",
                      dit->get()->getName().c_str());
        } else {
          CHKERR PetscSectionAddDof(*s, point, dd);
          int field = fields_map.at(field_name).first;
          CHKERR PetscSectionSetFieldDof(*s, point, field, dd);
        }

        ++point;

      } else {
        ++dit;
      }

    }
  }
  // cerr << "done " << proc << endl;
  CHKERR PetscSectionSetUp(*s);
  // cerr << "end " << proc << endl;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateProblemOrder(const std::string &problem,
                                               RowColData rc, int min_order,
                                               int max_order, IS *is) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(problem, &problem_ptr);

  typedef multi_index_container<
      boost::shared_ptr<NumeredDofEntity>,

      indexed_by<

          sequenced<>,

          ordered_non_unique<
              tag<Order_mi_tag>,
              const_mem_fun<NumeredDofEntity::interface_type_DofEntity,
                            ApproximationOrder, &NumeredDofEntity::getDofOrder>>

          >>
      NumeredDofEntity_order_view_multiIndex;

  const int rank = m_field.get_comm_rank();

  NumeredDofEntity_order_view_multiIndex dofs_part_view;
  auto insert_part_range = [&dofs_part_view, rank](auto &dofs) {
    dofs_part_view.insert(dofs_part_view.end(), dofs.lower_bound(rank),
                          dofs.upper_bound(rank));
  };

  switch (rc) {
  case ROW:
    insert_part_range(problem_ptr->numeredRowDofsPtr->get<Part_mi_tag>());
    break;
  case COL:
    insert_part_range(problem_ptr->numeredColDofsPtr->get<Part_mi_tag>());
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }

  auto lo = dofs_part_view.get<Order_mi_tag>().lower_bound(min_order);
  auto hi = dofs_part_view.get<Order_mi_tag>().upper_bound(max_order);
  const int size = std::distance(lo, hi);
  int *id;
  CHKERR PetscMalloc(size * sizeof(int), &id);
  int *id_it = id;
  for (; lo != hi; ++lo, ++id_it)
    *id_it = (*lo)->getPetscGlobalDofIdx();
  sort(id, &id[size]);

  CHKERR ISCreateGeneral(PETSC_COMM_WORLD, size, id, PETSC_OWN_POINTER, is);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateProblemFieldAndRank(
    const std::string &problem, RowColData rc, const std::string &field,
    int min_coeff_idx, int max_coeff_idx, IS *is, Range *ents_ptr) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(problem, &problem_ptr);
  const auto bit_number = m_field.get_field_bit_number(field);

  typedef NumeredDofEntity_multiIndex::index<Unique_mi_tag>::type DofsByUId;
  DofsByUId::iterator it, hi_it;
  int nb_loc_dofs;
  switch (rc) {
  case ROW:
    nb_loc_dofs = problem_ptr->getNbLocalDofsRow();
    it = problem_ptr->numeredRowDofsPtr->get<Unique_mi_tag>().lower_bound(
        FieldEntity::getLoBitNumberUId(bit_number));
    hi_it = problem_ptr->numeredRowDofsPtr->get<Unique_mi_tag>().upper_bound(
        FieldEntity::getHiBitNumberUId(bit_number));
    break;
  case COL:
    nb_loc_dofs = problem_ptr->getNbLocalDofsCol();
    it = problem_ptr->numeredColDofsPtr->get<Unique_mi_tag>().lower_bound(
        FieldEntity::getLoBitNumberUId(bit_number));
    hi_it = problem_ptr->numeredColDofsPtr->get<Unique_mi_tag>().upper_bound(
        FieldEntity::getHiBitNumberUId(bit_number));
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }



  std::vector<int> idx_vec;
  idx_vec.reserve(std::distance(it, hi_it));
  for (; it != hi_it; ++it) {

    auto true_if_dof_on_entity = [&]() {
      if (ents_ptr) {
        return ents_ptr->find((*it)->getEnt()) != ents_ptr->end();
      } else {
        return true;
      }
    };

    auto check = [&]() {
      const auto ceff_idx = (*it)->getDofCoeffIdx();
      if (

          (*it)->getPetscLocalDofIdx() >= nb_loc_dofs ||

          ceff_idx < min_coeff_idx || ceff_idx > max_coeff_idx

      )
        return false;
      else
        return true;
    };

    if (check()) {
      if (true_if_dof_on_entity()) {
        idx_vec.emplace_back((*it)->getPetscGlobalDofIdx());
      }
    }
  }

  int *id;
  CHKERR PetscMalloc(idx_vec.size() * sizeof(int), &id);
  std::copy(idx_vec.begin(), idx_vec.end(), id);
  CHKERR ISCreateGeneral(PETSC_COMM_WORLD, idx_vec.size(), id,
                         PETSC_OWN_POINTER, is);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateProblemFieldAndRank(
    const std::string &problem, RowColData rc, const std::string &field,
    int min_coeff_idx, int max_coeff_idx, SmartPetscObj<IS> &smart_is,
    Range *ents_ptr) const {
  MoFEMFunctionBegin;
  IS is;
  CHKERR isCreateProblemFieldAndRank(problem, rc, field, min_coeff_idx,
                                     max_coeff_idx, &is, ents_ptr);
  smart_is = SmartPetscObj<IS>(is);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateProblemFieldAndEntityType(
    const std::string &problem, RowColData rc, const std::string &field,
    EntityType low_type, EntityType hi_type, int min_coeff_idx,
    int max_coeff_idx, IS *is, Range *ents_ptr) const {
  const MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  EntityHandle field_meshset = m_field.get_field_meshset(field);
  Range ents;
  for (; low_type <= hi_type; ++low_type)
    CHKERR m_field.get_moab().get_entities_by_type(field_meshset, low_type,
                                                   ents, true);
  if (ents_ptr)
    ents = intersect(ents, *ents_ptr);
  CHKERR isCreateProblemFieldAndRank(problem, rc, field, min_coeff_idx,
                                     max_coeff_idx, is, &ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem, const std::string &x_field_name,
    RowColData x_rc, const std::string &y_problem,
    const std::string &y_field_name, RowColData y_rc, std::vector<int> &idx,
    std::vector<int> &idy) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *px_ptr;
  const Problem *py_ptr;
  MoFEMFunctionBegin;

  CHKERR m_field.get_problem(x_problem, &px_ptr);
  CHKERR m_field.get_problem(y_problem, &py_ptr);

  typedef multi_index_container<
      boost::shared_ptr<NumeredDofEntity>,

      indexed_by<

          sequenced<>,

          ordered_non_unique<
              tag<Composite_Ent_And_EntDofIdx_mi_tag>,
              composite_key<
                  NumeredDofEntity,
                  const_mem_fun<NumeredDofEntity::interface_type_DofEntity,
                                EntityHandle, &NumeredDofEntity::getEnt>,
                  const_mem_fun<NumeredDofEntity::interface_type_DofEntity,
                                DofIdx, &NumeredDofEntity::getEntDofIdx>>>

          >>
      NumeredDofEntity_view_multiIndex;

  NumeredDofEntity_view_multiIndex dofs_view;

  auto x_bit_number = m_field.get_field_bit_number(x_field_name);

  switch (x_rc) {
  case ROW:
    dofs_view.insert(dofs_view.end(),
                     px_ptr->numeredRowDofsPtr->get<Unique_mi_tag>().lower_bound(
                         FieldEntity::getLoBitNumberUId(x_bit_number)),
                     px_ptr->numeredRowDofsPtr->get<Unique_mi_tag>().upper_bound(
                         FieldEntity::getHiBitNumberUId(x_bit_number)));
    break;
  case COL:
    dofs_view.insert(dofs_view.end(),
                     px_ptr->numeredColDofsPtr->get<Unique_mi_tag>().lower_bound(
                         FieldEntity::getLoBitNumberUId(x_bit_number)),
                     px_ptr->numeredColDofsPtr->get<Unique_mi_tag>().upper_bound(
                         FieldEntity::getHiBitNumberUId(x_bit_number)));
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "only makes sense for ROWS and COLS");
  }

  decltype(py_ptr->numeredRowDofsPtr) dofs_ptr;
  switch (y_rc) {
  case ROW:
    dofs_ptr = py_ptr->numeredRowDofsPtr;
    break;
  case COL:
    dofs_ptr = py_ptr->numeredColDofsPtr;
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "only makes sense for ROWS and COLS");
  }

  std::map<int, int> global_dofs_map;
  const auto y_bit_number = m_field.get_field_bit_number(y_field_name);
  auto dit = dofs_ptr->get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId(y_bit_number));
  auto hi_dit = dofs_ptr->get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId(y_bit_number));
  const auto rank = m_field.get_comm_rank();
  for (; dit != hi_dit; ++dit) {
    if ((*dit)->getPart() == rank) {
      auto x_dit = dofs_view.get<Composite_Ent_And_EntDofIdx_mi_tag>().find(
          boost::make_tuple((*dit)->getEnt(), (*dit)->getEntDofIdx()));
      if (x_dit != dofs_view.get<Composite_Ent_And_EntDofIdx_mi_tag>().end()) {
        global_dofs_map[(*x_dit)->getPetscGlobalDofIdx()] =
            (*dit)->getPetscGlobalDofIdx();
      }
    }
  }

  idx.resize(global_dofs_map.size());
  idy.resize(global_dofs_map.size());
  {
    auto ix = idx.begin();
    auto iy = idy.begin();
    for (auto mit = global_dofs_map.begin(); mit != global_dofs_map.end();
         mit++, ix++, iy++) {
      *ix = mit->first;
      *iy = mit->second;
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem, const std::string &x_field_name,
    RowColData x_rc, const std::string &y_problem,
    const std::string &y_field_name, RowColData y_rc, IS *ix, IS *iy) const {
  MoFEMFunctionBegin;
  std::vector<int> idx(0), idy(0);
  CHKERR isCreateFromProblemFieldToOtherProblemField(
      x_problem, x_field_name, x_rc, y_problem, y_field_name, y_rc, idx, idy);
  if (ix != PETSC_NULL) {
    CHKERR ISCreateGeneral(PETSC_COMM_WORLD, idx.size(), &idx[0],
                           PETSC_COPY_VALUES, ix);
  }
  CHKERR ISCreateGeneral(PETSC_COMM_WORLD, idy.size(), &idy[0],
                         PETSC_COPY_VALUES, iy);
  if (dEbug) {
    ISView(*ix, PETSC_VIEWER_STDOUT_WORLD);
    ISView(*iy, PETSC_VIEWER_STDOUT_WORLD);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateFromProblemToOtherProblem(
    const std::string &x_problem, RowColData x_rc, const std::string &y_problem,
    RowColData y_rc, std::vector<int> &idx, std::vector<int> &idy) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *px_ptr;
  const Problem *py_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(x_problem, &px_ptr);
  CHKERR m_field.get_problem(y_problem, &py_ptr);
  NumeredDofEntityByLocalIdx::iterator y_dit, hi_y_dit;
  switch (y_rc) {
  case ROW:
    y_dit = py_ptr->numeredRowDofsPtr->get<PetscLocalIdx_mi_tag>().lower_bound(0);
    hi_y_dit = py_ptr->numeredRowDofsPtr->get<PetscLocalIdx_mi_tag>().lower_bound(
        py_ptr->getNbLocalDofsRow()); // should be lower
    break;
  case COL:
    y_dit = py_ptr->numeredColDofsPtr->get<PetscLocalIdx_mi_tag>().lower_bound(0);
    hi_y_dit = py_ptr->numeredColDofsPtr->get<PetscLocalIdx_mi_tag>().lower_bound(
        py_ptr->getNbLocalDofsCol()); // should be lower
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  const NumeredDofEntityByUId *x_numered_dofs_by_uid;
  switch (x_rc) {
  case ROW:
    x_numered_dofs_by_uid = &(px_ptr->numeredRowDofsPtr->get<Unique_mi_tag>());
    break;
  case COL:
    x_numered_dofs_by_uid = &(px_ptr->numeredColDofsPtr->get<Unique_mi_tag>());
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  for (; y_dit != hi_y_dit; y_dit++) {
    if ((*y_dit)->getPart() != (unsigned int)m_field.get_comm_rank()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    NumeredDofEntityByUId::iterator x_dit;
    x_dit = x_numered_dofs_by_uid->find((*y_dit)->getLocalUniqueId());
    if (x_dit == x_numered_dofs_by_uid->end())
      continue;
    idx.push_back((*x_dit)->getPetscGlobalDofIdx());
    idy.push_back((*y_dit)->getPetscGlobalDofIdx());
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateFromProblemToOtherProblem(
    const std::string &x_problem, RowColData x_rc, const std::string &y_problem,
    RowColData y_rc, IS *ix, IS *iy) const {
  MoFEMFunctionBegin;
  std::vector<int> idx(0), idy(0);
  CHKERR isCreateFromProblemToOtherProblem(x_problem, x_rc, y_problem, y_rc,
                                           idx, idy);
  CHKERR ISCreateGeneral(PETSC_COMM_WORLD, idx.size(), &idx[0],
                         PETSC_COPY_VALUES, ix);
  CHKERR ISCreateGeneral(PETSC_COMM_WORLD, idy.size(), &idy[0],
                         PETSC_COPY_VALUES, iy);
  if (dEbug) {
    ISView(*ix, PETSC_VIEWER_STDOUT_WORLD);
    ISView(*iy, PETSC_VIEWER_STDOUT_WORLD);
  }
  MoFEMFunctionReturn(0);
}
} // namespace MoFEM
