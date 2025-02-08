/** \file ISManager.cpp
 * \brief IS creating
 * \ingroup mofem_is_managers
 */


namespace MoFEM {

MoFEMErrorCode
ISManager::query_interface(boost::typeindex::type_index type_index,
                           UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = const_cast<ISManager *>(this);
  MoFEMFunctionReturnHot(0);
}

ISManager::ISManager(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)), dEbug(false) {}

MoFEMErrorCode ISManager::sectionCreate(const std::string problem_name,
                                        PetscSection *s,
                                        const RowColData row_col) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr = m_field.get_problem(problem_name);
  auto fields_ptr = m_field.get_fields();
  auto fe_ptr = m_field.get_finite_elements();
  MoFEMFunctionBegin;
  boost::shared_ptr<NumeredDofEntity_multiIndex> dofs;
  BitFieldId fields_ids(0);
  switch (row_col) {
  case ROW:
    dofs = problem_ptr->numeredRowDofsPtr;
    for (auto fit = fe_ptr->begin(); fit != fe_ptr->end(); fit++) {
      if ((fit->get()->getId() & problem_ptr->getBitFEId()).any()) {
        fields_ids |= fit->get()->getBitFieldIdRow();
      }
    }
    break;
  case COL:
    dofs = problem_ptr->numeredColDofsPtr;
    for (auto fit = fe_ptr->begin(); fit != fe_ptr->end(); fit++) {
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
    for (auto fit = fields_ptr->begin(); fit != fields_ptr->end(); fit++) {
      if ((fit->get()->getId() & fields_ids).any()) {
        fields_map[fit->get()->getName()].first = field++;
        fields_map[fit->get()->getName()].second = fit->get()->getNbOfCoeffs();
      }
    }
  }
  const int proc = m_field.get_comm_rank();
  CHKERR PetscSectionCreate(PETSC_COMM_WORLD, s);
  CHKERR PetscSectionSetNumFields(*s, fields_map.size());
  for (auto mit = fields_map.begin(); mit != fields_map.end(); mit++) {
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
          const DofIdx loc_idx = dit->get()->getPetscLocalDofIdx();
          if (loc_idx >= 0)
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

SmartPetscObj<PetscSection>
ISManager::sectionCreate(const std::string problem_name,
                         const RowColData row_col) const {

  PetscSection s;
  CHK_THROW_MESSAGE(sectionCreate(problem_name, &s, row_col),
                    "Section not created");
  return SmartPetscObj<PetscSection>(s, false);
}

MoFEMErrorCode ISManager::isCreateProblem(const std::string problem_name,
                                          RowColData rc, IS *is) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(problem_name, &problem_ptr);

  const int rank = m_field.get_comm_rank();

  decltype(problem_ptr->numeredRowDofsPtr) dofs_ptr;

  switch (rc) {
  case ROW:
    dofs_ptr = problem_ptr->numeredRowDofsPtr;
    break;
  case COL:
    dofs_ptr = problem_ptr->numeredColDofsPtr;
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }

  auto lo = dofs_ptr->get<Part_mi_tag>().lower_bound(rank);
  auto hi = dofs_ptr->get<Part_mi_tag>().upper_bound(rank);
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

MoFEMErrorCode ISManager::isCreateProblem(const std::string problem_name,
                                          RowColData rc,
                                          SmartPetscObj<IS> &is) const {
  MoFEMFunctionBegin;
  IS raw_is;
  CHKERR isCreateProblem(problem_name, rc, &raw_is);
  is = SmartPetscObj<IS>(raw_is);
  MoFEMFunctionReturn(0);
}

SmartPetscObj<IS> ISManager::isCreateProblem(const std::string problem_name,
                                             RowColData rc) const {
  SmartPetscObj<IS> is;
  CHK_THROW_MESSAGE(isCreateProblem(problem_name, rc, is), "IS not created");
  return is;
}

MoFEMErrorCode ISManager::isCreateProblemOrder(const std::string problem_name,
                                               RowColData rc, int min_order,
                                               int max_order, IS *is) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(problem_name, &problem_ptr);

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

  NumeredDofEntity_order_view_multiIndex dofs_by_order;
  auto insert_part_range = [&dofs_by_order, rank](auto &dofs) {
    dofs_by_order.insert(dofs_by_order.end(), dofs.lower_bound(rank),
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

  auto lo = dofs_by_order.get<Order_mi_tag>().lower_bound(min_order);
  auto hi = dofs_by_order.get<Order_mi_tag>().upper_bound(max_order);
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

MoFEMErrorCode ISManager::isCreateProblemOrder(const std::string problem_name,
                                               RowColData rc, int min_order,
                                               int max_order,
                                               SmartPetscObj<IS> &is) const {
  MoFEMFunctionBegin;
  IS raw_is;
  CHKERR isCreateProblemOrder(problem_name, rc, min_order, max_order, &raw_is);
  is = SmartPetscObj<IS>(raw_is);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateProblemFieldAndRank(
    const std::string problem_name, RowColData rc, const std::string field,
    int min_coeff_idx, int max_coeff_idx, IS *is, Range *ents_ptr) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(problem_name, &problem_ptr);
  const auto bit_number = m_field.get_field_bit_number(field);

  using DofsByUId = NumeredDofEntity_multiIndex::index<Unique_mi_tag>::type;
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
      const auto coeff_idx = (*it)->getDofCoeffIdx();
      if (

          (*it)->getPetscLocalDofIdx() >= nb_loc_dofs ||

          coeff_idx < min_coeff_idx || coeff_idx > max_coeff_idx

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
  CHKERR ISCreateGeneral(m_field.get_comm(), idx_vec.size(), id,
                         PETSC_OWN_POINTER, is);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateProblemFieldAndRank(
    const std::string problem_name, RowColData rc, const std::string field,
    int min_coeff_idx, int max_coeff_idx, SmartPetscObj<IS> &smart_is,
    Range *ents_ptr) const {
  MoFEMFunctionBegin;
  IS is;
  CHKERR isCreateProblemFieldAndRank(problem_name, rc, field, min_coeff_idx,
                                     max_coeff_idx, &is, ents_ptr);
  smart_is = SmartPetscObj<IS>(is);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateProblemBrokenFieldAndRank(
    const std::vector<boost::weak_ptr<NumeredDofEntity>> &dofs_vec,
    SmartPetscObj<IS> &smart_is) const {
  const MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;

  std::vector<int> idx_vec;
  idx_vec.reserve(dofs_vec.size());
  for (auto &dof : dofs_vec) {
    if (auto d = dof.lock()) {
      idx_vec.emplace_back(d->getPetscGlobalDofIdx());
    }
  }

  IS is_raw;
  int *id;
  CHKERR PetscMalloc(idx_vec.size() * sizeof(int), &id);
  std::copy(idx_vec.begin(), idx_vec.end(), id);
  CHKERR ISCreateGeneral(m_field.get_comm(), idx_vec.size(), id,
                         PETSC_OWN_POINTER, &is_raw);

  smart_is = SmartPetscObj<IS>(is_raw);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateProblemFieldAndRankLocal(
    const std::string problem_name, RowColData rc, const std::string field,
    int min_coeff_idx, int max_coeff_idx, IS *is, Range *ents_ptr) const {
  const MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(problem_name, &problem_ptr);
  const auto bit_number = m_field.get_field_bit_number(field);

  auto get_low_hi_uid = [&]() {
    return std::make_pair(FieldEntity::getLoBitNumberUId(bit_number),
                          FieldEntity::getHiBitNumberUId(bit_number));
  };

  auto get_low_hi_uid_by_entities = [&](auto f, auto s) {
    return std::make_pair(DofEntity::getLoFieldEntityUId(bit_number, f),
                          DofEntity::getHiFieldEntityUId(bit_number, s));
  };

  auto get_low_hi = [&](auto lo_uid, auto hi_uid) {
    using DofsByUId = NumeredDofEntity_multiIndex::index<Unique_mi_tag>::type;
    DofsByUId::iterator it, hi_it;
    switch (rc) {
    case ROW:
      it = problem_ptr->numeredRowDofsPtr->get<Unique_mi_tag>().lower_bound(
          lo_uid);
      hi_it = problem_ptr->numeredRowDofsPtr->get<Unique_mi_tag>().upper_bound(
          hi_uid);
      break;
    case COL:
      it = problem_ptr->numeredColDofsPtr->get<Unique_mi_tag>().lower_bound(
          lo_uid);
      hi_it = problem_ptr->numeredColDofsPtr->get<Unique_mi_tag>().upper_bound(
          hi_uid);
      break;
    default:
      THROW_MESSAGE("not implemented");
    }
    return std::make_pair(it, hi_it);
  };

  auto check = [&](auto it) {
    const auto coeff_idx = (*it)->getDofCoeffIdx();
    if (

        coeff_idx < min_coeff_idx || coeff_idx > max_coeff_idx

    )
      return false;
    else
      return true;
  };

  auto emplace_indices = [&](auto it, auto hi_it, auto &idx_vec) {
    for (; it != hi_it; ++it) {
      if (check(it))
        idx_vec.emplace_back((*it)->getPetscLocalDofIdx());
    }
  };

  auto [lo_uid, hi_uid] = get_low_hi_uid();
  auto [lo, hi] = get_low_hi(lo_uid, hi_uid);
  std::vector<int> idx_vec;
  idx_vec.reserve(std::distance(lo, hi));

  if (ents_ptr) {
    for (auto pit = ents_ptr->const_pair_begin();
         pit != ents_ptr->const_pair_end(); ++pit) {
      auto [lo_uid, hi_uid] =
          get_low_hi_uid_by_entities(pit->first, pit->second);
      auto [lo, hi] = get_low_hi(lo_uid, hi_uid);
      emplace_indices(lo, hi, idx_vec);
    }
  } else {
    emplace_indices(lo, hi, idx_vec);
  }

  int *id;
  CHKERR PetscMalloc(idx_vec.size() * sizeof(int), &id);
  std::copy(idx_vec.begin(), idx_vec.end(), id);
  CHKERR ISCreateGeneral(PETSC_COMM_SELF, idx_vec.size(), id, PETSC_OWN_POINTER,
                         is);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateProblemFieldAndRankLocal(
    const std::string problem_name, RowColData rc, const std::string field,
    int min_coeff_idx, int max_coeff_idx, SmartPetscObj<IS> &smart_is,
    Range *ents_ptr) const {
  MoFEMFunctionBegin;
  IS is;
  CHKERR isCreateProblemFieldAndRankLocal(problem_name, rc, field,
                                          min_coeff_idx, max_coeff_idx, &is,
                                          ents_ptr);
  smart_is = SmartPetscObj<IS>(is);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateProblemFieldAndEntityType(
    const std::string problem_name, RowColData rc, const std::string field,
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
  CHKERR isCreateProblemFieldAndRank(problem_name, rc, field, min_coeff_idx,
                                     max_coeff_idx, is, &ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateFromProblemFieldToOtherProblemField(
    const std::string x_problem, const std::string x_field_name,
    RowColData x_rc, const std::string y_problem,
    const std::string y_field_name, RowColData y_rc, std::vector<int> &idx,
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
    dofs_view.insert(
        dofs_view.end(),
        px_ptr->numeredRowDofsPtr->get<Unique_mi_tag>().lower_bound(
            FieldEntity::getLoBitNumberUId(x_bit_number)),
        px_ptr->numeredRowDofsPtr->get<Unique_mi_tag>().upper_bound(
            FieldEntity::getHiBitNumberUId(x_bit_number)));
    break;
  case COL:
    dofs_view.insert(
        dofs_view.end(),
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
    const std::string x_problem, const std::string x_field_name,
    RowColData x_rc, const std::string y_problem,
    const std::string y_field_name, RowColData y_rc, IS *ix, IS *iy) const {
  MoFEMFunctionBegin;
  const MoFEM::Interface &m_field = cOre;
  std::vector<int> idx(0), idy(0);
  CHKERR isCreateFromProblemFieldToOtherProblemField(
      x_problem, x_field_name, x_rc, y_problem, y_field_name, y_rc, idx, idy);
  if (ix != PETSC_NULL) {
    CHKERR ISCreateGeneral(m_field.get_comm(), idx.size(), &idx[0],
                           PETSC_COPY_VALUES, ix);
  }
  CHKERR ISCreateGeneral(m_field.get_comm(), idy.size(), &idy[0],
                         PETSC_COPY_VALUES, iy);
  if (dEbug) {
    ISView(*ix, PETSC_VIEWER_STDOUT_WORLD);
    ISView(*iy, PETSC_VIEWER_STDOUT_WORLD);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ISManager::isCreateFromProblemToOtherProblem(
    const std::string x_problem, RowColData x_rc, const std::string y_problem,
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
    y_dit =
        py_ptr->numeredRowDofsPtr->get<PetscLocalIdx_mi_tag>().lower_bound(0);
    hi_y_dit =
        py_ptr->numeredRowDofsPtr->get<PetscLocalIdx_mi_tag>().lower_bound(
            py_ptr->getNbLocalDofsRow()); // should be lower
    break;
  case COL:
    y_dit =
        py_ptr->numeredColDofsPtr->get<PetscLocalIdx_mi_tag>().lower_bound(0);
    hi_y_dit =
        py_ptr->numeredColDofsPtr->get<PetscLocalIdx_mi_tag>().lower_bound(
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
    const std::string x_problem, RowColData x_rc, const std::string y_problem,
    RowColData y_rc, IS *ix, IS *iy) const {
  const MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  std::vector<int> idx(0), idy(0);
  CHKERR isCreateFromProblemToOtherProblem(x_problem, x_rc, y_problem, y_rc,
                                           idx, idy);
  CHKERR ISCreateGeneral(m_field.get_comm(), idx.size(), &idx[0],
                         PETSC_COPY_VALUES, ix);
  CHKERR ISCreateGeneral(m_field.get_comm(), idy.size(), &idy[0],
                         PETSC_COPY_VALUES, iy);
  if (dEbug) {
    ISView(*ix, PETSC_VIEWER_STDOUT_WORLD);
    ISView(*iy, PETSC_VIEWER_STDOUT_WORLD);
  }
  MoFEMFunctionReturn(0);
}
} // namespace MoFEM
