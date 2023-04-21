/** \file ProblemsManager.cpp
 * \brief Managing complexities for problem
 * \ingroup mofem_problems_manager
 */

namespace MoFEM {

#define ProblemManagerFunctionBegin                                            \
  MoFEMFunctionBegin;                                                          \
  MOFEM_LOG_CHANNEL("WORLD");                                                  \
  MOFEM_LOG_CHANNEL("SYNC");                                                   \
  MOFEM_LOG_FUNCTION();                                                        \
  MOFEM_LOG_TAG("SYNC", "ProblemsManager");                                    \
  MOFEM_LOG_TAG("WORLD", "ProblemsManager")

struct IdxDataType {
  IdxDataType(const UId uid, const int dof) {
    bcopy(&uid, dAta, 4 * sizeof(int));
    dAta[4] = dof;
  }

private:
  int dAta[5];
};

struct IdxDataTypePtr {
  IdxDataTypePtr(const int *ptr) : pTr(ptr) {}
  inline int getDofIdx() const {
    int global_dof = pTr[4];
    return global_dof;
  }
  inline UId getUId() const {
    unsigned int b0, b1, b2, b3;
    bcopy(&pTr[0], &b0, sizeof(int));
    bcopy(&pTr[1], &b1, sizeof(int));
    bcopy(&pTr[2], &b2, sizeof(int));
    bcopy(&pTr[3], &b3, sizeof(int));
    UId uid = static_cast<UId>(b0) | static_cast<UId>(b1) << 8 * sizeof(int) |
              static_cast<UId>(b2) << 16 * sizeof(int) |
              static_cast<UId>(b3) << 24 * sizeof(int);
    return uid;
  }

private:
  const int *pTr;
};

MoFEMErrorCode
ProblemsManager::query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const {
  *iface = const_cast<ProblemsManager *>(this);
  return 0;
}

ProblemsManager::ProblemsManager(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)),
      buildProblemFromFields(PETSC_FALSE),
      synchroniseProblemEntities(PETSC_FALSE) {
  PetscLogEventRegister("ProblemsManager", 0, &MOFEM_EVENT_ProblemsManager);
}

MoFEMErrorCode ProblemsManager::getOptions() {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBeginHot;
  CHKERR PetscOptionsBegin(m_field.get_comm(), "", "Problem manager", "none");
  {
    CHKERR PetscOptionsBool(
        "-problem_build_from_fields",
        "Add DOFs to problem directly from fields not through DOFs on elements",
        "", buildProblemFromFields, &buildProblemFromFields, NULL);
  }
  ierr = PetscOptionsEnd();
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ProblemsManager::partitionMesh(
    const Range &ents, const int dim, const int adj_dim, const int n_parts,
    Tag *th_vertex_weights, Tag *th_edge_weights, Tag *th_part_weights,
    int verb, const bool debug) {
  return static_cast<MoFEM::Interface &>(cOre)
      .getInterface<CommInterface>()
      ->partitionMesh(ents, dim, adj_dim, n_parts, th_vertex_weights,
                      th_edge_weights, th_part_weights, verb, debug);
}

MoFEMErrorCode ProblemsManager::buildProblem(const std::string name,
                                             const bool square_matrix,
                                             int verb) {

  MoFEM::Interface &m_field = cOre;
  ProblemManagerFunctionBegin;
  if (!(cOre.getBuildMoFEM() & (1 << 0)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "fields not build");
  if (!(cOre.getBuildMoFEM() & (1 << 1)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "FEs not build");
  if (!(cOre.getBuildMoFEM() & (1 << 2)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "adjacencies not build");
  const Problem *problem_ptr;
  CHKERR m_field.get_problem(name, &problem_ptr);
  CHKERR buildProblem(const_cast<Problem *>(problem_ptr), square_matrix, verb);
  cOre.getBuildMoFEM() |= 1 << 3; // It is assumed that user who uses this
                                  // function knows what he is doing
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::buildProblem(Problem *problem_ptr,
                                             const bool square_matrix,
                                             int verb) {
  MoFEM::Interface &m_field = cOre;
  auto dofs_field_ptr = m_field.get_dofs();
  auto ents_field_ptr = m_field.get_field_ents();
  auto adjacencies_ptr = m_field.get_ents_elements_adjacency();
  ProblemManagerFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_ProblemsManager, 0, 0, 0, 0);

  // Note: Only allowed changes on problem_ptr structure which not influence
  // multi-index.

  if (problem_ptr->getBitRefLevel().none()) {
    SETERRQ1(PETSC_COMM_SELF, 1, "problem <%s> refinement level not set",
             problem_ptr->getName().c_str());
  }
  CHKERR m_field.clear_problem(problem_ptr->getName());

  // zero finite elements
  problem_ptr->numeredFiniteElementsPtr->clear();

  DofEntity_multiIndex_active_view dofs_rows, dofs_cols;
  auto make_rows_and_cols_view = [&](auto &dofs_rows, auto &dofs_cols) {
    MoFEMFunctionBeginHot;
    for (auto it = ents_field_ptr->begin(); it != ents_field_ptr->end(); ++it) {

      const auto uid = (*it)->getLocalUniqueId();

      auto r = adjacencies_ptr->get<Unique_mi_tag>().equal_range(uid);
      for (auto lo = r.first; lo != r.second; ++lo) {

        if ((lo->getBitFEId() & problem_ptr->getBitFEId()).any()) {
          std::array<bool, 2> row_col = {false, false};

          const BitRefLevel &prb_bit = problem_ptr->getBitRefLevel();
          const BitRefLevel &prb_mask = problem_ptr->getBitRefLevelMask();
          const BitRefLevel &fe_bit = lo->entFePtr->getBitRefLevel();

          // if entity is not problem refinement level
          if ((fe_bit & prb_mask) != fe_bit)
            continue;
          if ((fe_bit & prb_bit).none())
            continue;

          auto add_to_view = [&](auto &nb_dofs, auto &view, auto rc) {
            auto dit = nb_dofs->lower_bound(uid);
            decltype(dit) hi_dit;
            if (dit != nb_dofs->end()) {
              hi_dit = nb_dofs->upper_bound(
                  uid | static_cast<UId>(MAX_DOFS_ON_ENTITY - 1));
              view.insert(dit, hi_dit);
              row_col[rc] = true;
            }
          };

          if ((lo->entFePtr->getBitFieldIdRow() & (*it)->getId()).any())
            add_to_view(dofs_field_ptr, dofs_rows, ROW);

          if (!square_matrix)
            if ((lo->entFePtr->getBitFieldIdCol() & (*it)->getId()).any())
              add_to_view(dofs_field_ptr, dofs_cols, COL);

          if (square_matrix && row_col[ROW])
            break;
          else if (row_col[ROW] && row_col[COL])
            break;
        }
      }
    }
    MoFEMFunctionReturnHot(0);
  };

  CHKERR make_rows_and_cols_view(dofs_rows, dofs_cols);

  // Add row dofs to problem
  {
    // zero rows
    problem_ptr->nbDofsRow = 0;
    problem_ptr->nbLocDofsRow = 0;
    problem_ptr->nbGhostDofsRow = 0;

    // add dofs for rows
    DofEntity_multiIndex_active_view::nth_index<0>::type::iterator miit,
        hi_miit;
    hi_miit = dofs_rows.get<0>().end();

    int count_dofs = dofs_rows.get<1>().count(true);
    boost::shared_ptr<std::vector<NumeredDofEntity>> dofs_array =
        boost::shared_ptr<std::vector<NumeredDofEntity>>(
            new std::vector<NumeredDofEntity>());
    problem_ptr->getRowDofsSequence()->push_back(dofs_array);
    dofs_array->reserve(count_dofs);
    miit = dofs_rows.get<0>().begin();
    for (; miit != hi_miit; miit++) {
      if ((*miit)->getActive()) {
        dofs_array->emplace_back(*miit);
        dofs_array->back().dofIdx = (problem_ptr->nbDofsRow)++;
      }
    }
    auto hint = problem_ptr->numeredRowDofsPtr->end();
    for (auto &v : *dofs_array) {
      hint = problem_ptr->numeredRowDofsPtr->emplace_hint(hint, dofs_array, &v);
    }
  }

  // Add col dofs to problem
  if (!square_matrix) {
    // zero cols
    problem_ptr->nbDofsCol = 0;
    problem_ptr->nbLocDofsCol = 0;
    problem_ptr->nbGhostDofsCol = 0;

    // add dofs for cols
    DofEntity_multiIndex_active_view::nth_index<0>::type::iterator miit,
        hi_miit;
    hi_miit = dofs_cols.get<0>().end();

    int count_dofs = 0;
    miit = dofs_cols.get<0>().begin();
    for (; miit != hi_miit; miit++) {
      if (!(*miit)->getActive()) {
        continue;
      }
      count_dofs++;
    }

    boost::shared_ptr<std::vector<NumeredDofEntity>> dofs_array =
        boost::shared_ptr<std::vector<NumeredDofEntity>>(
            new std::vector<NumeredDofEntity>());
    problem_ptr->getColDofsSequence()->push_back(dofs_array);
    dofs_array->reserve(count_dofs);
    miit = dofs_cols.get<0>().begin();
    for (; miit != hi_miit; miit++) {
      if (!(*miit)->getActive()) {
        continue;
      }
      dofs_array->emplace_back(*miit);
      dofs_array->back().dofIdx = problem_ptr->nbDofsCol++;
    }
    auto hint = problem_ptr->numeredColDofsPtr->end();
    for (auto &v : *dofs_array) {
      hint = problem_ptr->numeredColDofsPtr->emplace_hint(hint, dofs_array, &v);
    }
  } else {
    problem_ptr->numeredColDofsPtr = problem_ptr->numeredRowDofsPtr;
    problem_ptr->nbLocDofsCol = problem_ptr->nbLocDofsRow;
    problem_ptr->nbDofsCol = problem_ptr->nbDofsRow;
  }

  // job done, some debugging and postprocessing
  if (verb >= VERBOSE) {
    MOFEM_LOG("SYNC", Sev::verbose)
        << problem_ptr->getName() << " Nb. local dofs "
        << problem_ptr->numeredRowDofsPtr->size() << " by "
        << problem_ptr->numeredColDofsPtr->size();
    MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::verbose);
  }

  if (verb >= NOISY) {
    MOFEM_LOG("SYNC", Sev::noisy)
        << "FEs row dofs " << *problem_ptr << " Nb. row dof "
        << problem_ptr->getNbDofsRow();
    for (auto &miit : *problem_ptr->numeredRowDofsPtr)
      MOFEM_LOG("SYNC", Sev::noisy) << *miit;

    MOFEM_LOG("SYNC", Sev::noisy)
        << "FEs col dofs " << *problem_ptr << " Nb. col dof "
        << problem_ptr->getNbDofsCol();
    for (auto &miit : *problem_ptr->numeredColDofsPtr)
      MOFEM_LOG("SYNC", Sev::noisy) << *miit;
    MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::noisy);
  }

  cOre.getBuildMoFEM() |= Core::BUILD_PROBLEM; // It is assumed that user who
                                               // uses this function knows
                                               // what he is doing

  PetscLogEventEnd(MOFEM_EVENT_ProblemsManager, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::buildProblemOnDistributedMesh(
    const std::string name, const bool square_matrix, int verb) {
  MoFEM::Interface &m_field = cOre;
  ProblemManagerFunctionBegin;

  if (!((cOre.getBuildMoFEM()) & Core::BUILD_FIELD))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "fields not build");
  if (!((cOre.getBuildMoFEM()) & Core::BUILD_FE))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "FEs not build");
  if (!((cOre.getBuildMoFEM()) & Core::BUILD_ADJ))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "adjacencies not build");

  const Problem *problem_ptr;
  CHKERR m_field.get_problem(name, &problem_ptr);
  CHKERR buildProblemOnDistributedMesh(const_cast<Problem *>(problem_ptr),
                                       square_matrix, verb);

  cOre.getBuildMoFEM() |= Core::BUILD_PROBLEM;
  cOre.getBuildMoFEM() |= Core::PARTITION_PROBLEM;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::buildProblemOnDistributedMesh(
    Problem *problem_ptr, const bool square_matrix, int verb) {
  MoFEM::Interface &m_field = cOre;
  auto fields_ptr = m_field.get_fields();
  auto fe_ptr = m_field.get_finite_elements();
  auto fe_ent_ptr = m_field.get_ents_finite_elements();
  auto ents_field_ptr = m_field.get_field_ents();
  auto dofs_field_ptr = m_field.get_dofs();
  ProblemManagerFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_ProblemsManager, 0, 0, 0, 0);

  // clear data structures
  CHKERR m_field.clear_problem(problem_ptr->getName());

  CHKERR getOptions();

  if (problem_ptr->getBitRefLevel().none())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
             "problem <%s> refinement level not set",
             problem_ptr->getName().c_str());

  int loop_size = 2;
  if (square_matrix) {
    loop_size = 1;
    problem_ptr->numeredColDofsPtr = problem_ptr->numeredRowDofsPtr;
  } else if (problem_ptr->numeredColDofsPtr == problem_ptr->numeredRowDofsPtr) {
    problem_ptr->numeredColDofsPtr =
        boost::shared_ptr<NumeredDofEntity_multiIndex>(
            new NumeredDofEntity_multiIndex());
  }

  const BitRefLevel &prb_bit = problem_ptr->getBitRefLevel();
  const BitRefLevel &prb_mask = problem_ptr->getBitRefLevelMask();

  // // get rows and cols dofs view based on data on elements
  DofEntity_multiIndex_active_view dofs_rows, dofs_cols;

  // Add DOFs to problem by visiting all elements and adding DOFs from
  // elements to the problem
  if (buildProblemFromFields == PETSC_FALSE) {

    auto make_rows_and_cols_view = [&](auto &dofs_rows, auto &dofs_cols) {
      auto ents_field_ptr = m_field.get_field_ents();
      auto adjacencies_ptr = m_field.get_ents_elements_adjacency();
      MoFEMFunctionBeginHot;
      for (auto it = ents_field_ptr->begin(); it != ents_field_ptr->end();
           ++it) {

        const auto uid = (*it)->getLocalUniqueId();

        auto r = adjacencies_ptr->get<Unique_mi_tag>().equal_range(uid);
        for (auto lo = r.first; lo != r.second; ++lo) {

          if ((lo->getBitFEId() & problem_ptr->getBitFEId()).any()) {
            std::array<bool, 2> row_col = {false, false};

            const BitRefLevel &fe_bit = lo->entFePtr->getBitRefLevel();

            // if entity is not problem refinement level
            if ((fe_bit & prb_bit).any() && (fe_bit & prb_mask) == fe_bit) {

              auto add_to_view = [&](auto dofs, auto &view, auto rc) {
                auto dit = dofs->lower_bound(uid);
                decltype(dit) hi_dit;
                if (dit != dofs->end()) {
                  hi_dit = dofs->upper_bound(
                      uid | static_cast<UId>(MAX_DOFS_ON_ENTITY - 1));
                  view.insert(dit, hi_dit);
                  row_col[rc] = true;
                }
              };

              if ((lo->entFePtr->getBitFieldIdRow() & (*it)->getId()).any())
                add_to_view(dofs_field_ptr, dofs_rows, ROW);

              if (!square_matrix)
                if ((lo->entFePtr->getBitFieldIdCol() & (*it)->getId()).any())
                  add_to_view(dofs_field_ptr, dofs_cols, COL);

              if (square_matrix && row_col[ROW])
                break;
              else if (row_col[ROW] && row_col[COL])
                break;
            }
          }
        }
      }
      MoFEMFunctionReturnHot(0);
    };

    CHKERR make_rows_and_cols_view(dofs_rows, dofs_cols);
  }

  // Add DOFS to the problem by searching all the filedes, and adding to
  // problem owned or shared DOFs
  if (buildProblemFromFields == PETSC_TRUE) {
    // Get fields IDs on elements
    BitFieldId fields_ids_row, fields_ids_col;
    for (auto fit = fe_ptr->begin(); fit != fe_ptr->end(); fit++) {
      if ((fit->get()->getId() & problem_ptr->getBitFEId()).any()) {
        fields_ids_row |= fit->get()->getBitFieldIdRow();
        fields_ids_col |= fit->get()->getBitFieldIdCol();
      }
    }
    // Get fields DOFs
    for (auto fit = fields_ptr->begin(); fit != fields_ptr->end(); fit++) {
      if ((fit->get()->getId() & (fields_ids_row | fields_ids_col)).any()) {

        auto dit = dofs_field_ptr->get<Unique_mi_tag>().lower_bound(
            FieldEntity::getLoBitNumberUId((*fit)->getBitNumber()));
        auto hi_dit = dofs_field_ptr->get<Unique_mi_tag>().upper_bound(
            FieldEntity::getHiBitNumberUId((*fit)->getBitNumber()));

        for (; dit != hi_dit; dit++) {

          const int owner_proc = dit->get()->getOwnerProc();
          if (owner_proc != m_field.get_comm_rank()) {
            const unsigned char pstatus = dit->get()->getPStatus();
            if (pstatus == 0) {
              continue;
            }
          }

          const auto &dof_bit = (*dit)->getBitRefLevel();
          // if entity is not problem refinement level
          if ((dof_bit & prb_bit).any() && (dof_bit & prb_mask) == dof_bit) {

            if ((fit->get()->getId() & fields_ids_row).any()) {
              dofs_rows.insert(*dit);
            }
            if (!square_matrix) {
              if ((fit->get()->getId() & fields_ids_col).any()) {
                dofs_cols.insert(*dit);
              }
            }

          }
        }
      }
    }
  }

  if (synchroniseProblemEntities) {
    // Get fields IDs on elements
    BitFieldId fields_ids_row, fields_ids_col;
    BitFieldId *fields_ids[2] = {&fields_ids_row, &fields_ids_col};
    for (FiniteElement_multiIndex::iterator fit = fe_ptr->begin();
         fit != fe_ptr->end(); fit++) {
      if ((fit->get()->getId() & problem_ptr->getBitFEId()).any()) {
        fields_ids_row |= fit->get()->getBitFieldIdRow();
        fields_ids_col |= fit->get()->getBitFieldIdCol();
      }
    }

    DofEntity_multiIndex_active_view *dofs_ptr[] = {&dofs_rows, &dofs_cols};
    for (int ss = 0; ss != ((square_matrix) ? 1 : 2); ++ss) {
      DofEntity_multiIndex_active_view::nth_index<1>::type::iterator miit,
          hi_miit;
      miit = dofs_ptr[ss]->get<1>().lower_bound(1);
      hi_miit = dofs_ptr[ss]->get<1>().upper_bound(1);
      Range ents_to_synchronise;
      for (; miit != hi_miit; ++miit) {
        if (miit->get()->getEntDofIdx() != 0)
          continue;
        ents_to_synchronise.insert(miit->get()->getEnt());
      }
      Range tmp_ents = ents_to_synchronise;
      CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(
          ents_to_synchronise, nullptr, verb);
      ents_to_synchronise = subtract(ents_to_synchronise, tmp_ents);
      for (auto fit = fields_ptr->begin(); fit != fields_ptr->end(); fit++) {
        if ((fit->get()->getId() & *fields_ids[ss]).any()) {
          const auto bit_number = (*fit)->getBitNumber();
          for (Range::pair_iterator pit = ents_to_synchronise.pair_begin();
               pit != ents_to_synchronise.pair_end(); ++pit) {
            const auto f = pit->first;
            const auto s = pit->second;
            const auto lo_uid =
                FieldEntity::getLocalUniqueIdCalculate(bit_number, f);
            const auto hi_uid =
                FieldEntity::getLocalUniqueIdCalculate(bit_number, s);

            auto dit = dofs_field_ptr->get<Unique_mi_tag>().lower_bound(lo_uid);
            auto hi_dit =
                dofs_field_ptr->get<Unique_mi_tag>().upper_bound(hi_uid);

            dofs_ptr[ss]->insert(dit, hi_dit);
          }
        }
      }
    }
  }

  // add dofs for rows and cols and set ownership
  DofEntity_multiIndex_active_view *dofs_ptr[] = {&dofs_rows, &dofs_cols};
  boost::shared_ptr<NumeredDofEntity_multiIndex> numered_dofs_ptr[] = {
      problem_ptr->numeredRowDofsPtr, problem_ptr->numeredColDofsPtr};
  int *nbdof_ptr[] = {&problem_ptr->nbDofsRow, &problem_ptr->nbDofsCol};
  int *local_nbdof_ptr[] = {&problem_ptr->nbLocDofsRow,
                            &problem_ptr->nbLocDofsCol};
  int *ghost_nbdof_ptr[] = {&problem_ptr->nbGhostDofsRow,
                            &problem_ptr->nbGhostDofsCol};
  for (int ss = 0; ss < 2; ss++) {
    *(nbdof_ptr[ss]) = 0;
    *(local_nbdof_ptr[ss]) = 0;
    *(ghost_nbdof_ptr[ss]) = 0;
  }

  // Loop over dofs on rows and columns and add to multi-indices in dofs
  // problem structure,  set partition for each dof
  int nb_local_dofs[] = {0, 0};
  for (int ss = 0; ss < loop_size; ss++) {
    DofEntity_multiIndex_active_view::nth_index<1>::type::iterator miit,
        hi_miit;
    miit = dofs_ptr[ss]->get<1>().lower_bound(1);
    hi_miit = dofs_ptr[ss]->get<1>().upper_bound(1);
    for (; miit != hi_miit; miit++) {
      int owner_proc = (*miit)->getOwnerProc();
      if (owner_proc == m_field.get_comm_rank()) {
        nb_local_dofs[ss]++;
      }
    }
  }

  // get layout
  int start_ranges[2], end_ranges[2];
  for (int ss = 0; ss != loop_size; ss++) {
    PetscLayout layout;
    CHKERR PetscLayoutCreate(m_field.get_comm(), &layout);
    CHKERR PetscLayoutSetBlockSize(layout, 1);
    CHKERR PetscLayoutSetLocalSize(layout, nb_local_dofs[ss]);
    CHKERR PetscLayoutSetUp(layout);
    CHKERR PetscLayoutGetSize(layout, &*nbdof_ptr[ss]);
    CHKERR PetscLayoutGetRange(layout, &start_ranges[ss], &end_ranges[ss]);
    CHKERR PetscLayoutDestroy(&layout);
  }
  if (square_matrix) {
    nbdof_ptr[1] = nbdof_ptr[0];
    nb_local_dofs[1] = nb_local_dofs[0];
    start_ranges[1] = start_ranges[0];
    end_ranges[1] = end_ranges[0];
  }

  // if(sizeof(UId) != SIZEOFUID) {
  //   SETERRQ2(
  //     PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
  //     "check size of UId, size of UId is %u != %u",
  //     sizeof(UId),SIZEOFUID
  //   );
  // }

  // set local and global indices on own dofs
  const size_t idx_data_type_size = sizeof(IdxDataType);
  const size_t data_block_size = idx_data_type_size / sizeof(int);

  if (sizeof(IdxDataType) % sizeof(int)) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  std::vector<std::vector<IdxDataType>> ids_data_packed_rows(
      m_field.get_comm_size()),
      ids_data_packed_cols(m_field.get_comm_size());

  // Loop over dofs on this processor and prepare those dofs to send on
  // another proc
  for (int ss = 0; ss != loop_size; ++ss) {

    DofEntity_multiIndex_active_view::nth_index<0>::type::iterator miit,
        hi_miit;
    hi_miit = dofs_ptr[ss]->get<0>().end();

    boost::shared_ptr<std::vector<NumeredDofEntity>> dofs_array =
        boost::shared_ptr<std::vector<NumeredDofEntity>>(
            new std::vector<NumeredDofEntity>());
    int nb_dofs_to_add = 0;
    miit = dofs_ptr[ss]->get<0>().begin();
    for (; miit != hi_miit; ++miit) {
      // Only set global idx for dofs on this processor part
      if (!(miit->get()->getActive()))
        continue;
      ++nb_dofs_to_add;
    }
    dofs_array->reserve(nb_dofs_to_add);
    if (ss == 0) {
      problem_ptr->getRowDofsSequence()->push_back(dofs_array);
    } else {
      problem_ptr->getColDofsSequence()->push_back(dofs_array);
    }

    int &local_idx = *local_nbdof_ptr[ss];
    miit = dofs_ptr[ss]->get<0>().begin();
    for (; miit != hi_miit; ++miit) {

      // Only set global idx for dofs on this processor part
      if (!(miit->get()->getActive()))
        continue;

      dofs_array->emplace_back(*miit);

      int owner_proc = dofs_array->back().getOwnerProc();
      if (owner_proc < 0) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }

      if (owner_proc != m_field.get_comm_rank()) {
        dofs_array->back().pArt = owner_proc;
        dofs_array->back().dofIdx = -1;
        dofs_array->back().petscGloablDofIdx = -1;
        dofs_array->back().petscLocalDofIdx = -1;
      } else {

        // Set part and indexes
        int glob_idx = start_ranges[ss] + local_idx;
        dofs_array->back().pArt = owner_proc;
        dofs_array->back().dofIdx = glob_idx;
        dofs_array->back().petscGloablDofIdx = glob_idx;
        dofs_array->back().petscLocalDofIdx = local_idx;
        local_idx++;

        unsigned char pstatus = dofs_array->back().getPStatus();
        // check id dof is shared, if that is a case global idx added to data
        // structure and is sended to other processors
        if (pstatus > 0) {

          for (int proc = 0;
               proc < MAX_SHARING_PROCS &&
               -1 != dofs_array->back().getSharingProcsPtr()[proc];
               proc++) {
            // make it different for rows and columns when needed
            if (ss == 0) {
              ids_data_packed_rows[dofs_array->back()
                                       .getSharingProcsPtr()[proc]]
                  .emplace_back(dofs_array->back().getGlobalUniqueId(),
                                glob_idx);
            } else {
              ids_data_packed_cols[dofs_array->back()
                                       .getSharingProcsPtr()[proc]]
                  .emplace_back(dofs_array->back().getGlobalUniqueId(),
                                glob_idx);
            }
            if (!(pstatus & PSTATUS_MULTISHARED)) {
              break;
            }
          }
        }
      }
    }

    auto hint = numered_dofs_ptr[ss]->end();
    for (auto &v : *dofs_array)
      hint = numered_dofs_ptr[ss]->emplace_hint(hint, dofs_array, &v);
  }
  if (square_matrix) {
    local_nbdof_ptr[1] = local_nbdof_ptr[0];
  }

  int nsends_rows = 0, nsends_cols = 0;
  // Non zero lengths[i] represent a message to i of length lengths[i].
  std::vector<int> lengths_rows(m_field.get_comm_size()),
      lengths_cols(m_field.get_comm_size());
  lengths_rows.clear();
  lengths_cols.clear();
  for (int proc = 0; proc < m_field.get_comm_size(); proc++) {
    lengths_rows[proc] = ids_data_packed_rows[proc].size() * data_block_size;
    lengths_cols[proc] = ids_data_packed_cols[proc].size() * data_block_size;
    if (!ids_data_packed_rows[proc].empty())
      nsends_rows++;
    if (!ids_data_packed_cols[proc].empty())
      nsends_cols++;
  }

  MPI_Status *status;
  CHKERR PetscMalloc1(m_field.get_comm_size(), &status);

  // Do rows
  int nrecvs_rows;    // number of messages received
  int *onodes_rows;   // list of node-ids from which messages are expected
  int *olengths_rows; // corresponding message lengths
  int **rbuf_row;     // must bee freed by user

  // make sure it is a PETSc comm
  MPI_Comm comm;
  CHKERR PetscCommDuplicate(m_field.get_comm(), &comm, NULL);

  {

    // rows

    // Computes the number of messages a node expects to receive
    CHKERR PetscGatherNumberOfMessages(comm, NULL, &lengths_rows[0],
                                       &nrecvs_rows);
    // std::cerr << nrecvs_rows << std::endl;

    // Computes info about messages that a MPI-node will receive, including
    // (from-id,length) pairs for each message.
    CHKERR PetscGatherMessageLengths(comm, nsends_rows, nrecvs_rows,
                                     &lengths_rows[0], &onodes_rows,
                                     &olengths_rows);

    // Gets a unique new tag from a PETSc communicator. All processors that
    // share the communicator MUST call this routine EXACTLY the same number
    // of times. This tag should only be used with the current objects
    // communicator; do NOT use it with any other MPI communicator.
    int tag_row;
    CHKERR PetscCommGetNewTag(comm, &tag_row);

    // Allocate a buffer sufficient to hold messages of size specified in
    // olengths. And post Irecvs on these buffers using node info from onodes
    MPI_Request *r_waits_row; // must bee freed by user
    // rbuf has a pointers to messeges. It has size of of nrecvs (number of
    // messages) +1. In the first index a block is allocated,
    // such that rbuf[i] = rbuf[i-1]+olengths[i-1].

    CHKERR PetscPostIrecvInt(comm, tag_row, nrecvs_rows, onodes_rows,
                             olengths_rows, &rbuf_row, &r_waits_row);
    CHKERR PetscFree(onodes_rows);

    MPI_Request *s_waits_row; // status of sens messages
    CHKERR PetscMalloc1(nsends_rows, &s_waits_row);

    // Send messeges
    for (int proc = 0, kk = 0; proc < m_field.get_comm_size(); proc++) {
      if (!lengths_rows[proc])
        continue; // no message to send to this proc
      CHKERR MPI_Isend(&(ids_data_packed_rows[proc])[0], // buffer to send
                       lengths_rows[proc],               // message length
                       MPIU_INT, proc,                   // to proc
                       tag_row, comm, s_waits_row + kk);
      kk++;
    }

    if (nrecvs_rows) {
      CHKERR MPI_Waitall(nrecvs_rows, r_waits_row, status);
    }
    if (nsends_rows) {
      CHKERR MPI_Waitall(nsends_rows, s_waits_row, status);
    }

    CHKERR PetscFree(r_waits_row);
    CHKERR PetscFree(s_waits_row);
  }

  // cols
  int nrecvs_cols = nrecvs_rows;
  int *olengths_cols = olengths_rows;
  PetscInt **rbuf_col = rbuf_row;
  if (!square_matrix) {

    // Computes the number of messages a node expects to receive
    CHKERR PetscGatherNumberOfMessages(comm, NULL, &lengths_cols[0],
                                       &nrecvs_cols);

    // Computes info about messages that a MPI-node will receive, including
    // (from-id,length) pairs for each message.
    int *onodes_cols;
    CHKERR PetscGatherMessageLengths(comm, nsends_cols, nrecvs_cols,
                                     &lengths_cols[0], &onodes_cols,
                                     &olengths_cols);

    // Gets a unique new tag from a PETSc communicator.
    int tag_col;
    CHKERR PetscCommGetNewTag(comm, &tag_col);

    // Allocate a buffer sufficient to hold messages of size specified in
    // olengths. And post Irecvs on these buffers using node info from onodes
    MPI_Request *r_waits_col; // must bee freed by user
    CHKERR PetscPostIrecvInt(comm, tag_col, nrecvs_cols, onodes_cols,
                             olengths_cols, &rbuf_col, &r_waits_col);
    CHKERR PetscFree(onodes_cols);

    MPI_Request *s_waits_col; // status of sens messages
    CHKERR PetscMalloc1(nsends_cols, &s_waits_col);

    // Send messages
    for (int proc = 0, kk = 0; proc < m_field.get_comm_size(); proc++) {
      if (!lengths_cols[proc])
        continue; // no message to send to this proc
      CHKERR MPI_Isend(&(ids_data_packed_cols[proc])[0], // buffer to send
                       lengths_cols[proc],               // message length
                       MPIU_INT, proc,                   // to proc
                       tag_col, comm, s_waits_col + kk);
      kk++;
    }

    if (nrecvs_cols) {
      CHKERR MPI_Waitall(nrecvs_cols, r_waits_col, status);
    }
    if (nsends_cols) {
      CHKERR MPI_Waitall(nsends_cols, s_waits_col, status);
    }

    CHKERR PetscFree(r_waits_col);
    CHKERR PetscFree(s_waits_col);
  }

  CHKERR PetscCommDestroy(&comm);
  CHKERR PetscFree(status);

  DofEntity_multiIndex_global_uid_view dofs_glob_uid_view;
  auto hint = dofs_glob_uid_view.begin();
  for (auto dof : *m_field.get_dofs())
    dofs_glob_uid_view.emplace_hint(hint, dof);

  // set values received from other processors
  for (int ss = 0; ss != loop_size; ++ss) {

    int nrecvs;
    int *olengths;
    int **data_procs;
    if (ss == 0) {
      nrecvs = nrecvs_rows;
      olengths = olengths_rows;
      data_procs = rbuf_row;
    } else {
      nrecvs = nrecvs_cols;
      olengths = olengths_cols;
      data_procs = rbuf_col;
    }

    UId uid;
    for (int kk = 0; kk != nrecvs; ++kk) {
      int len = olengths[kk];
      int *data_from_proc = data_procs[kk];
      for (int dd = 0; dd < len; dd += data_block_size) {
        uid = IdxDataTypePtr(&data_from_proc[dd]).getUId();
        auto ddit = dofs_glob_uid_view.find(uid);

        if (PetscUnlikely(ddit == dofs_glob_uid_view.end())) {

#ifndef NDEBUG
          MOFEM_LOG("SELF", Sev::error)
              << "DOF is shared or multishared between processors. For example "
                 "if order of field on given entity is set in inconsistently, "
                 "has different value on two processor, error such as this is "
                 "triggered";

          MOFEM_LOG("SELF", Sev::error) << "UId " << uid << " is not found";
          const auto owner_proc = FieldEntity::getOwnerFromUniqueId(uid);
          MOFEM_LOG("SELF", Sev::error)
              << "Problematic UId owner proc is " << owner_proc;
          const auto uid_handle = FieldEntity::getHandleFromUniqueId(uid);
          MOFEM_LOG("SELF", Sev::error)
              << "Problematic UId entity owning processor handle is "
              << uid_handle << " entity type "
              << moab::CN::EntityTypeName(type_from_handle(uid_handle));
          const auto uid_bit_number =
              FieldEntity::getFieldBitNumberFromUniqueId(uid);
          MOFEM_LOG("SELF", Sev::error)
              << "Problematic UId field is "
              << m_field.get_field_name(
                     field_bit_from_bit_number(uid_bit_number));

          FieldEntity_multiIndex_global_uid_view field_view;
          field_view.insert(ents_field_ptr->begin(), ents_field_ptr->end());
          auto fe_it = field_view.find(FieldEntity::getGlobalUniqueIdCalculate(
              owner_proc, uid_bit_number, uid_handle));
          if (fe_it == field_view.end()) {
            MOFEM_LOG("SELF", Sev::error)
                << "Also, no field entity in database for given global UId";
          } else {
            MOFEM_LOG("SELF", Sev::error) << "Field entity in databse exist "
                                             "(but have no DOF wih give UId";
            MOFEM_LOG("SELF", Sev::error) << **fe_it;

            // Save file with missing entity
            auto error_file_name =
                "error_with_missing_entity_" +
                boost::lexical_cast<std::string>(m_field.get_comm_rank()) +
                ".vtk";
            MOFEM_LOG("SELF", Sev::error)
                << "Look to file < " << error_file_name
                << " > it contains entity with missing DOF.";

            auto tmp_msh = get_temp_meshset_ptr(m_field.get_moab());
            const auto local_fe_ent = (*fe_it)->getEnt();
            CHKERR m_field.get_moab().add_entities(*tmp_msh, &local_fe_ent, 1);
            CHKERR m_field.get_moab().write_file(error_file_name.c_str(), "VTK",
                                                 "", tmp_msh->get_ptr(), 1);
          }
#endif

          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "DOF with global UId not found (Compile code in Debug to "
                  "learn more about problem");
        }

        auto dit = numered_dofs_ptr[ss]->find((*ddit)->getLocalUniqueId());

        if (dit != numered_dofs_ptr[ss]->end()) {

          int global_idx = IdxDataTypePtr(&data_from_proc[dd]).getDofIdx();
#ifndef NDEBUG
          if (PetscUnlikely(global_idx < 0))
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "received negative dof");
#endif
          bool success;
          success = numered_dofs_ptr[ss]->modify(
              dit, NumeredDofEntity_mofem_index_change(global_idx));

#ifndef NDEBUG
          if (PetscUnlikely(!success))
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
#endif
          success = numered_dofs_ptr[ss]->modify(
              dit, NumeredDofEntity_part_and_glob_idx_change((*dit)->getPart(),
                                                             global_idx));
#ifndef NDEBUG
          if (PetscUnlikely(!success))
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
#endif
        };

#ifndef NDEBUG
        if (PetscUnlikely(ddit->get()->getPStatus() == 0)) {

          // Dof is shared on this processor, however there is no element
          // which have this dof. If DOF is not shared and received from other
          // processor, but not marked as a shared on other that means that is
          // data inconstancy and error should be thrown.

          std::ostringstream zz;
          zz << **ddit << std::endl;
          SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                   "data inconsistency, dofs is not shared, but received "
                   "from other proc\n"
                   "%s",
                   zz.str().c_str());
        }
#endif
      }
    }
  }

  if (square_matrix) {
    (problem_ptr->nbDofsCol) = (problem_ptr->nbDofsRow);
    (problem_ptr->nbLocDofsCol) = (problem_ptr->nbLocDofsRow);
  }

  CHKERR PetscFree(olengths_rows);
  CHKERR PetscFree(rbuf_row[0]);
  CHKERR PetscFree(rbuf_row);
  if (!square_matrix) {
    CHKERR PetscFree(olengths_cols);
    CHKERR PetscFree(rbuf_col[0]);
    CHKERR PetscFree(rbuf_col);
  }

  if (square_matrix) {
    if (numered_dofs_ptr[0]->size() != numered_dofs_ptr[1]->size()) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "data inconsistency for square_matrix %d!=%d",
               numered_dofs_ptr[0]->size(), numered_dofs_ptr[1]->size());
    }
    if (problem_ptr->numeredRowDofsPtr != problem_ptr->numeredColDofsPtr) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "data inconsistency for square_matrix");
    }
  }

  CHKERR printPartitionedProblem(problem_ptr, verb);
  CHKERR debugPartitionedProblem(problem_ptr, verb);

  PetscLogEventEnd(MOFEM_EVENT_ProblemsManager, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::buildSubProblem(
    const std::string out_name,

    const std::vector<std::string> &fields_row,
    const std::vector<std::string> &fields_col,

    const std::string main_problem, const bool square_matrix,

    const map<std::string, boost::shared_ptr<Range>> *entityMapRow,
    const map<std::string, boost::shared_ptr<Range>> *entityMapCol,

    int verb) {
  MoFEM::Interface &m_field = cOre;
  auto problems_ptr = m_field.get_problems();
  ProblemManagerFunctionBegin;

  CHKERR m_field.clear_problem(out_name);

  // get reference to all problems
  using ProblemByName = decltype(problems_ptr->get<Problem_mi_tag>());
  auto &problems_by_name =
      const_cast<ProblemByName &>(problems_ptr->get<Problem_mi_tag>());

  // get iterators to out problem, i.e. build problem
  auto out_problem_it = problems_by_name.find(out_name);
  if (out_problem_it == problems_by_name.end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "subproblem with name < %s > not defined (top tip check spelling)",
             out_name.c_str());
  }
  // get iterator to main problem, i.e. out problem is subproblem of main
  // problem
  auto main_problem_it = problems_by_name.find(main_problem);
  if (main_problem_it == problems_by_name.end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "problem of subproblem with name < %s > not defined (top tip "
             "check spelling)",
             main_problem.c_str());
  }

  // get dofs for row & columns for out problem,
  boost::shared_ptr<NumeredDofEntity_multiIndex> out_problem_dofs[] = {
      out_problem_it->numeredRowDofsPtr, out_problem_it->numeredColDofsPtr};
  // get dofs for row & columns for main problem
  boost::shared_ptr<NumeredDofEntity_multiIndex> main_problem_dofs[] = {
      main_problem_it->numeredRowDofsPtr, main_problem_it->numeredColDofsPtr};
  // get local indices counter
  int *nb_local_dofs[] = {&out_problem_it->nbLocDofsRow,
                          &out_problem_it->nbLocDofsCol};
  // get global indices counter
  int *nb_dofs[] = {&out_problem_it->nbDofsRow, &out_problem_it->nbDofsCol};

  // set number of ghost nodes to zero
  {
    out_problem_it->nbGhostDofsRow = 0;
    out_problem_it->nbGhostDofsCol = 0;
  }

  // put rows & columns field names in array
  std::vector<std::string> fields[] = {fields_row, fields_col};
  const map<std::string, boost::shared_ptr<Range>> *entityMap[] = {
      entityMapRow, entityMapCol};

  // make data structure fos sub-problem data
  out_problem_it->subProblemData =
      boost::make_shared<Problem::SubProblemData>();

  // Loop over rows and columns
  for (int ss = 0; ss != (square_matrix ? 1 : 2); ++ss) {

    // reset dofs and columns counters
    (*nb_local_dofs[ss]) = 0;
    (*nb_dofs[ss]) = 0;
    // clear arrays
    out_problem_dofs[ss]->clear();

    // If DOFs are cleared clear finite elements too.
    out_problem_it->numeredFiniteElementsPtr->clear();

    // get dofs by field name and insert them in out problem multi-indices
    for (auto field : fields[ss]) {

      // Following reserve memory in sequences, only two allocations are here,
      // once for array of objects, next for array of shared pointers

      // aliased sequence of pointer is killed with element
      boost::shared_ptr<std::vector<NumeredDofEntity>> dofs_array =
          boost::make_shared<std::vector<NumeredDofEntity>>();
      // reserve memory for field  dofs
      if (!ss)
        out_problem_it->getRowDofsSequence()->emplace_back(dofs_array);
      else
        out_problem_it->getColDofsSequence()->emplace_back(dofs_array);

      // create elements objects
      auto bit_number = m_field.get_field_bit_number(field);

      auto add_dit_to_dofs_array = [&](auto &dit) {
        if (dit->get()->getPetscGlobalDofIdx() >= 0)
          dofs_array->emplace_back(
              dit->get()->getDofEntityPtr(), dit->get()->getPetscGlobalDofIdx(),
              dit->get()->getPetscGlobalDofIdx(),
              dit->get()->getPetscLocalDofIdx(), dit->get()->getPart());
      };

      auto get_dafult_dof_range = [&]() {
        auto dit = main_problem_dofs[ss]->get<Unique_mi_tag>().lower_bound(
            FieldEntity::getLoBitNumberUId(bit_number));
        auto hi_dit = main_problem_dofs[ss]->get<Unique_mi_tag>().upper_bound(
            FieldEntity::getHiBitNumberUId(bit_number));
        return std::make_pair(dit, hi_dit);
      };

      if (entityMap[ss]) {
        auto mit = entityMap[ss]->find(field);

        if (mit != entityMap[ss]->end()) {
          for (auto p = mit->second->pair_begin(); p != mit->second->pair_end();
               ++p) {
            const auto lo_ent = p->first;
            const auto hi_ent = p->second;
            auto dit = main_problem_dofs[ss]->get<Unique_mi_tag>().lower_bound(
                DofEntity::getLoFieldEntityUId(bit_number, lo_ent));
            auto hi_dit =
                main_problem_dofs[ss]->get<Unique_mi_tag>().upper_bound(
                    DofEntity::getHiFieldEntityUId(bit_number, hi_ent));
            dofs_array->reserve(std::distance(dit, hi_dit));
            for (; dit != hi_dit; dit++) {
              add_dit_to_dofs_array(dit);
            }
          }
        } else {
          auto [dit, hi_dit] = get_dafult_dof_range();
          dofs_array->reserve(std::distance(dit, hi_dit));
          for (; dit != hi_dit; dit++)
            add_dit_to_dofs_array(dit);
        }
      } else {
        auto [dit, hi_dit] = get_dafult_dof_range();
        dofs_array->reserve(std::distance(dit, hi_dit));
        for (; dit != hi_dit; dit++)
          add_dit_to_dofs_array(dit);
      }

      // fill multi-index
      auto hint = out_problem_dofs[ss]->end();
      for (auto &v : *dofs_array)
        hint = out_problem_dofs[ss]->emplace_hint(hint, dofs_array, &v);
    }
    // Set local indexes
    {
      auto dit = out_problem_dofs[ss]->get<Idx_mi_tag>().begin();
      auto hi_dit = out_problem_dofs[ss]->get<Idx_mi_tag>().end();
      for (; dit != hi_dit; dit++) {
        int idx = -1; // if dof is not part of partition, set local index to -1
        if (dit->get()->getPart() == (unsigned int)m_field.get_comm_rank()) {
          idx = (*nb_local_dofs[ss])++;
        }
        bool success = out_problem_dofs[ss]->modify(
            out_problem_dofs[ss]->project<0>(dit),
            NumeredDofEntity_local_idx_change(idx));
        if (!success) {
          SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                  "operation unsuccessful");
        }
      };
    }
    // Set global indexes, compress global indices
    {
      auto dit =
          out_problem_dofs[ss]->get<PetscLocalIdx_mi_tag>().lower_bound(0);
      auto hi_dit =
          out_problem_dofs[ss]->get<PetscLocalIdx_mi_tag>().upper_bound(
              out_problem_dofs[ss]->size());
      const int nb = std::distance(dit, hi_dit);
      // get main problem global indices
      std::vector<int> main_indices(nb);
      for (auto it = main_indices.begin(); dit != hi_dit; dit++, it++) {
        *it = dit->get()->getPetscGlobalDofIdx();
      }
      // create is with global dofs
      IS is;
      CHKERR ISCreateGeneral(m_field.get_comm(), nb, &*main_indices.begin(),
                             PETSC_USE_POINTER, &is);
      // create map form main problem global indices to out problem global
      // indices
      AO ao;
      CHKERR AOCreateMappingIS(is, PETSC_NULL, &ao);
      if (ss == 0) {
        IS is_dup;
        CHKERR ISDuplicate(is, &is_dup);
        out_problem_it->getSubData()->rowIs = SmartPetscObj<IS>(is_dup, false);
        out_problem_it->getSubData()->rowMap = SmartPetscObj<AO>(ao, true);
      } else {
        IS is_dup;
        CHKERR ISDuplicate(is, &is_dup);
        out_problem_it->getSubData()->colIs = SmartPetscObj<IS>(is_dup, false);
        out_problem_it->getSubData()->colMap = SmartPetscObj<AO>(ao, true);
      }
      CHKERR AOApplicationToPetscIS(ao, is);
      // set global number of DOFs
      CHKERR ISGetSize(is, nb_dofs[ss]);
      CHKERR ISDestroy(&is);
      // set out problem global indices after applying map
      dit = out_problem_dofs[ss]->get<PetscLocalIdx_mi_tag>().lower_bound(0);
      for (std::vector<int>::iterator it = main_indices.begin(); dit != hi_dit;
           dit++, it++) {
        bool success = out_problem_dofs[ss]->modify(
            out_problem_dofs[ss]->project<0>(dit),
            NumeredDofEntity_part_and_all_indices_change(
                dit->get()->getPart(), *it, *it,
                dit->get()->getPetscLocalDofIdx()));
        if (!success) {
          SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                  "operation unsuccessful");
        }
      }
      // set global indices to nodes not on this part
      {
        NumeredDofEntityByLocalIdx::iterator dit =
            out_problem_dofs[ss]->get<PetscLocalIdx_mi_tag>().lower_bound(-1);
        NumeredDofEntityByLocalIdx::iterator hi_dit =
            out_problem_dofs[ss]->get<PetscLocalIdx_mi_tag>().upper_bound(-1);
        const int nb = std::distance(dit, hi_dit);
        std::vector<int> main_indices_non_local(nb);
        for (auto it = main_indices_non_local.begin(); dit != hi_dit;
             dit++, it++) {
          *it = dit->get()->getPetscGlobalDofIdx();
        }
        IS is;
        CHKERR ISCreateGeneral(m_field.get_comm(), nb,
                               &*main_indices_non_local.begin(),
                               PETSC_USE_POINTER, &is);
        CHKERR AOApplicationToPetscIS(ao, is);
        CHKERR ISDestroy(&is);
        dit = out_problem_dofs[ss]->get<PetscLocalIdx_mi_tag>().lower_bound(-1);
        for (auto it = main_indices_non_local.begin(); dit != hi_dit;
             dit++, it++) {
          bool success = out_problem_dofs[ss]->modify(
              out_problem_dofs[ss]->project<0>(dit),
              NumeredDofEntity_part_and_all_indices_change(
                  dit->get()->getPart(), dit->get()->getDofIdx(), *it,
                  dit->get()->getPetscLocalDofIdx()));
          if (!success) {
            SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                    "operation unsuccessful");
          }
        }
      }
      CHKERR AODestroy(&ao);
    }
  }

  if (square_matrix) {
    out_problem_it->numeredColDofsPtr = out_problem_it->numeredRowDofsPtr;
    out_problem_it->nbLocDofsCol = out_problem_it->nbLocDofsRow;
    out_problem_it->nbDofsCol = out_problem_it->nbDofsRow;
    out_problem_it->getSubData()->colIs = out_problem_it->getSubData()->rowIs;
    out_problem_it->getSubData()->colMap = out_problem_it->getSubData()->rowMap;
  }

  CHKERR printPartitionedProblem(&*out_problem_it, verb);
  CHKERR debugPartitionedProblem(&*out_problem_it, verb);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::buildComposedProblem(
    const std::string out_name, const std::vector<std::string> add_row_problems,
    const std::vector<std::string> add_col_problems, const bool square_matrix,
    int verb) {
  if (!(cOre.getBuildMoFEM() & Core::BUILD_FIELD))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "fields not build");
  if (!(cOre.getBuildMoFEM() & Core::BUILD_FE))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "FEs not build");
  if (!(cOre.getBuildMoFEM() & Core::BUILD_ADJ))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "adjacencies not build");
  MoFEM::Interface &m_field = cOre;
  auto problems_ptr = m_field.get_problems();
  ProblemManagerFunctionBegin;

  CHKERR m_field.clear_problem(out_name);
  // get reference to all problems
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemByName;
  ProblemByName &problems_by_name =
      const_cast<ProblemByName &>(problems_ptr->get<Problem_mi_tag>());

  // Get iterators to out problem, i.e. build problem
  ProblemByName::iterator out_problem_it = problems_by_name.find(out_name);
  if (out_problem_it == problems_by_name.end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "problem with name < %s > not defined (top tip check spelling)",
             out_name.c_str());
  }
  // Make data structure for composed-problem data
  out_problem_it->composedProblemsData =
      boost::make_shared<ComposedProblemsData>();
  boost::shared_ptr<ComposedProblemsData> cmp_prb_data =
      out_problem_it->getComposedProblemsData();

  const std::vector<std::string> *add_prb[] = {&add_row_problems,
                                               &add_col_problems};
  std::vector<const Problem *> *add_prb_ptr[] = {&cmp_prb_data->rowProblemsAdd,
                                                 &cmp_prb_data->colProblemsAdd};
  std::vector<SmartPetscObj<IS>> *add_prb_is[] = {&cmp_prb_data->rowIs,
                                                  &cmp_prb_data->colIs};

  // Get local indices counter
  int *nb_local_dofs[] = {&out_problem_it->nbLocDofsRow,
                          &out_problem_it->nbLocDofsCol};
  // Get global indices counter
  int *nb_dofs[] = {&out_problem_it->nbDofsRow, &out_problem_it->nbDofsCol};

  // Set number of ghost nodes to zero
  {
    out_problem_it->nbDofsRow = 0;
    out_problem_it->nbDofsCol = 0;
    out_problem_it->nbLocDofsRow = 0;
    out_problem_it->nbLocDofsCol = 0;
    out_problem_it->nbGhostDofsRow = 0;
    out_problem_it->nbGhostDofsCol = 0;
  }
  int nb_dofs_reserve[] = {0, 0};

  // Loop over rows and columns in the main problem and sub-problems
  for (int ss = 0; ss != ((square_matrix) ? 1 : 2); ss++) {
    add_prb_ptr[ss]->reserve(add_prb[ss]->size());
    add_prb_is[ss]->reserve(add_prb[ss]->size());
    for (std::vector<std::string>::const_iterator vit = add_prb[ss]->begin();
         vit != add_prb[ss]->end(); vit++) {
      ProblemByName::iterator prb_it = problems_by_name.find(*vit);
      if (prb_it == problems_by_name.end()) {
        SETERRQ1(
            PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "problem with name < %s > not defined (top tip check spelling)",
            vit->c_str());
      }
      add_prb_ptr[ss]->push_back(&*prb_it);
      // set number of dofs on rows and columns
      if (ss == 0) {
        // row
        *nb_dofs[ss] += add_prb_ptr[ss]->back()->getNbDofsRow();
        *nb_local_dofs[ss] += add_prb_ptr[ss]->back()->getNbLocalDofsRow();
        nb_dofs_reserve[ss] +=
            add_prb_ptr[ss]->back()->numeredRowDofsPtr->size();
      } else {
        // column
        *nb_dofs[ss] += add_prb_ptr[ss]->back()->getNbDofsCol();
        *nb_local_dofs[ss] += add_prb_ptr[ss]->back()->getNbLocalDofsCol();
        nb_dofs_reserve[ss] +=
            add_prb_ptr[ss]->back()->numeredColDofsPtr->size();
      }
    }
  }
  // if squre problem, rows and columns are the same
  if (square_matrix) {
    add_prb_ptr[1]->reserve(add_prb_ptr[0]->size());
    add_prb_is[1]->reserve(add_prb_ptr[0]->size());
    out_problem_it->numeredColDofsPtr = out_problem_it->numeredRowDofsPtr;
    *nb_dofs[1] = *nb_dofs[0];
    *nb_local_dofs[1] = *nb_local_dofs[0];
  }

  // reserve memory for dofs
  boost::shared_ptr<std::vector<NumeredDofEntity>> dofs_array[2];
  // Reserve memory
  for (int ss = 0; ss != ((square_matrix) ? 1 : 2); ss++) {
    dofs_array[ss] = boost::make_shared<std::vector<NumeredDofEntity>>();
    dofs_array[ss]->reserve(nb_dofs_reserve[ss]);
    if (!ss)
      out_problem_it->getRowDofsSequence()->emplace_back(dofs_array[ss]);
    else
      out_problem_it->getColDofsSequence()->emplace_back(dofs_array[ss]);
  }

  // Push back DOFs
  for (int ss = 0; ss != ((square_matrix) ? 1 : 2); ss++) {
    NumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator
        dit,
        hi_dit;
    int shift_glob = 0;
    int shift_loc = 0;
    for (unsigned int pp = 0; pp != add_prb_ptr[ss]->size(); pp++) {
      PetscInt *dofs_out_idx_ptr;
      int nb_local_dofs = (*add_prb_ptr[ss])[pp]->getNbLocalDofsRow();
      CHKERR PetscMalloc(nb_local_dofs * sizeof(int), &dofs_out_idx_ptr);
      if (ss == 0) {
        dit = (*add_prb_ptr[ss])[pp]
                  ->numeredRowDofsPtr->get<PetscGlobalIdx_mi_tag>()
                  .begin();
        hi_dit = (*add_prb_ptr[ss])[pp]
                     ->numeredRowDofsPtr->get<PetscGlobalIdx_mi_tag>()
                     .end();
      } else {
        dit = (*add_prb_ptr[ss])[pp]
                  ->numeredColDofsPtr->get<PetscGlobalIdx_mi_tag>()
                  .begin();
        hi_dit = (*add_prb_ptr[ss])[pp]
                     ->numeredColDofsPtr->get<PetscGlobalIdx_mi_tag>()
                     .end();
      }
      int is_nb = 0;
      for (; dit != hi_dit; dit++) {
        const BitRefLevel &prb_bit = out_problem_it->getBitRefLevel();
        const BitRefLevel &prb_mask = out_problem_it->getBitRefLevelMask();
        const BitRefLevel &dof_bit = dit->get()->getBitRefLevel();
        if ((dof_bit & prb_bit).none() || ((dof_bit & prb_mask) != dof_bit))
          continue;
        const int rank = m_field.get_comm_rank();
        const int part = dit->get()->getPart();
        const int glob_idx = shift_glob + dit->get()->getPetscGlobalDofIdx();
        const int loc_idx =
            (part == rank) ? (shift_loc + dit->get()->getPetscLocalDofIdx())
                           : -1;
        dofs_array[ss]->emplace_back(dit->get()->getDofEntityPtr(), glob_idx,
                                     glob_idx, loc_idx, part);
        if (part == rank) {
          dofs_out_idx_ptr[is_nb++] = glob_idx;
        }
      }
      if (is_nb > nb_local_dofs) {
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "Data inconsistency");
      }
      IS is;
      CHKERR ISCreateGeneral(m_field.get_comm(), is_nb, dofs_out_idx_ptr,
                             PETSC_OWN_POINTER, &is);
      auto smart_is = SmartPetscObj<IS>(is);
      (*add_prb_is[ss]).push_back(smart_is);
      if (ss == 0) {
        shift_glob += (*add_prb_ptr[ss])[pp]->getNbDofsRow();
        shift_loc += (*add_prb_ptr[ss])[pp]->getNbLocalDofsRow();
      } else {
        shift_glob += (*add_prb_ptr[ss])[pp]->getNbDofsCol();
        shift_loc += (*add_prb_ptr[ss])[pp]->getNbLocalDofsCol();
      }
      if (square_matrix) {
        (*add_prb_ptr[1]).push_back((*add_prb_ptr[0])[pp]);
        (*add_prb_is[1]).push_back(smart_is);
      }
    }
  }

  if ((*add_prb_is[1]).size() != (*add_prb_is[0]).size()) {
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY, "Data inconsistency");
  }

  // Insert DOFs to problem multi-index
  for (int ss = 0; ss != ((square_matrix) ? 1 : 2); ss++) {
    auto hint = (ss == 0) ? out_problem_it->numeredRowDofsPtr->end()
                          : out_problem_it->numeredColDofsPtr->end();
    for (auto &v : *dofs_array[ss])
      hint = (ss == 0) ? out_problem_it->numeredRowDofsPtr->emplace_hint(
                             hint, dofs_array[ss], &v)
                       : out_problem_it->numeredColDofsPtr->emplace_hint(
                             hint, dofs_array[ss], &v);
  }

  // Compress DOFs
  *nb_dofs[0] = 0;
  *nb_dofs[1] = 0;
  *nb_local_dofs[0] = 0;
  *nb_local_dofs[1] = 0;
  for (int ss = 0; ss != ((square_matrix) ? 1 : 2); ss++) {

    boost::shared_ptr<NumeredDofEntity_multiIndex> dofs_ptr;
    if (ss == 0) {
      dofs_ptr = out_problem_it->numeredRowDofsPtr;
    } else {
      dofs_ptr = out_problem_it->numeredColDofsPtr;
    }
    NumeredDofEntityByUId::iterator dit, hi_dit;
    dit = dofs_ptr->get<Unique_mi_tag>().begin();
    hi_dit = dofs_ptr->get<Unique_mi_tag>().end();
    std::vector<int> idx;
    idx.reserve(std::distance(dit, hi_dit));
    // set dofs in order entity and dof number on entity
    for (; dit != hi_dit; dit++) {
      if (dit->get()->getPart() == (unsigned int)m_field.get_comm_rank()) {
        bool success = dofs_ptr->get<Unique_mi_tag>().modify(
            dit, NumeredDofEntity_local_idx_change((*nb_local_dofs[ss])++));
        if (!success) {
          SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
        }
        idx.push_back(dit->get()->getPetscGlobalDofIdx());
      } else {
        if (dit->get()->getPetscLocalDofIdx() != -1) {
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "local index should be negative");
        }
      }
    }
    if (square_matrix) {
      *nb_local_dofs[1] = *nb_local_dofs[0];
    }

    // set new dofs mapping
    IS is;
    CHKERR ISCreateGeneral(m_field.get_comm(), idx.size(), &*idx.begin(),
                           PETSC_USE_POINTER, &is);
    CHKERR ISGetSize(is, nb_dofs[ss]);
    if (square_matrix) {
      *nb_dofs[1] = *nb_dofs[0];
    }

    AO ao;
    CHKERR AOCreateMappingIS(is, PETSC_NULL, &ao);
    for (unsigned int pp = 0; pp != (*add_prb_is[ss]).size(); pp++)
      CHKERR AOApplicationToPetscIS(ao, (*add_prb_is[ss])[pp]);

    // Set DOFs numeration
    {
      std::vector<int> idx_new;
      idx_new.reserve(dofs_ptr->size());
      for (NumeredDofEntityByUId::iterator dit =
               dofs_ptr->get<Unique_mi_tag>().begin();
           dit != dofs_ptr->get<Unique_mi_tag>().end(); dit++) {
        idx_new.push_back(dit->get()->getPetscGlobalDofIdx());
      }
      // set new global dofs numeration
      IS is_new;
      CHKERR ISCreateGeneral(m_field.get_comm(), idx_new.size(),
                             &*idx_new.begin(), PETSC_USE_POINTER, &is_new);
      CHKERR AOApplicationToPetscIS(ao, is_new);
      // set global indices to multi-index
      std::vector<int>::iterator vit = idx_new.begin();
      for (NumeredDofEntityByUId::iterator dit =
               dofs_ptr->get<Unique_mi_tag>().begin();
           dit != dofs_ptr->get<Unique_mi_tag>().end(); dit++) {
        bool success =
            dofs_ptr->modify(dit, NumeredDofEntity_part_and_glob_idx_change(
                                      dit->get()->getPart(), *(vit++)));
        if (!success) {
          SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
        }
      }
      CHKERR ISDestroy(&is_new);
    }
    CHKERR ISDestroy(&is);
    CHKERR AODestroy(&ao);
  }

  CHKERR printPartitionedProblem(&*out_problem_it, verb);
  CHKERR debugPartitionedProblem(&*out_problem_it, verb);

  // Inidcate that porble has been build
  cOre.getBuildMoFEM() |= Core::BUILD_PROBLEM;
  cOre.getBuildMoFEM() |= Core::PARTITION_PROBLEM;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::partitionSimpleProblem(const std::string name,
                                                       int verb) {

  MoFEM::Interface &m_field = cOre;
  auto problems_ptr = m_field.get_problems();
  ProblemManagerFunctionBegin;

  if (!(cOre.getBuildMoFEM() & Core::BUILD_FIELD))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "fields not build");
  if (!(cOre.getBuildMoFEM() & Core::BUILD_FE))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "FEs not build");
  if (!(cOre.getBuildMoFEM() & Core::BUILD_ADJ))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "adjacencies not build");
  if (!(cOre.getBuildMoFEM() & Core::BUILD_PROBLEM))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "adjacencies not build");
  MOFEM_LOG("WORLD", Sev::verbose) << "Simple partition problem " << name;

  // find p_miit
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemByName;
  ProblemByName &problems_set =
      const_cast<ProblemByName &>(problems_ptr->get<Problem_mi_tag>());
  ProblemByName::iterator p_miit = problems_set.find(name);
  if (p_miit == problems_set.end()) {
    SETERRQ1(PETSC_COMM_SELF, 1,
             "problem < %s > is not found (top tip: check spelling)",
             name.c_str());
  }
  typedef boost::multi_index::index<NumeredDofEntity_multiIndex,
                                    Idx_mi_tag>::type NumeredDofEntitysByIdx;
  NumeredDofEntitysByIdx &dofs_row_by_idx =
      p_miit->numeredRowDofsPtr->get<Idx_mi_tag>();
  NumeredDofEntitysByIdx &dofs_col_by_idx =
      p_miit->numeredColDofsPtr->get<Idx_mi_tag>();
  boost::multi_index::index<NumeredDofEntity_multiIndex,
                            Idx_mi_tag>::type::iterator miit_row,
      hi_miit_row;
  boost::multi_index::index<NumeredDofEntity_multiIndex,
                            Idx_mi_tag>::type::iterator miit_col,
      hi_miit_col;
  DofIdx &nb_row_local_dofs = p_miit->nbLocDofsRow;
  DofIdx &nb_row_ghost_dofs = p_miit->nbGhostDofsRow;
  nb_row_local_dofs = 0;
  nb_row_ghost_dofs = 0;
  DofIdx &nb_col_local_dofs = p_miit->nbLocDofsCol;
  DofIdx &nb_col_ghost_dofs = p_miit->nbGhostDofsCol;
  nb_col_local_dofs = 0;
  nb_col_ghost_dofs = 0;

  bool square_matrix = false;
  if (p_miit->numeredRowDofsPtr == p_miit->numeredColDofsPtr) {
    square_matrix = true;
  }

  // get row range of local indices
  PetscLayout layout_row;
  const int *ranges_row;

  DofIdx nb_dofs_row = dofs_row_by_idx.size();
  CHKERR PetscLayoutCreate(m_field.get_comm(), &layout_row);
  CHKERR PetscLayoutSetBlockSize(layout_row, 1);
  CHKERR PetscLayoutSetSize(layout_row, nb_dofs_row);
  CHKERR PetscLayoutSetUp(layout_row);
  CHKERR PetscLayoutGetRanges(layout_row, &ranges_row);
  // get col range of local indices
  PetscLayout layout_col;
  const int *ranges_col;
  if (!square_matrix) {
    DofIdx nb_dofs_col = dofs_col_by_idx.size();
    CHKERR PetscLayoutCreate(m_field.get_comm(), &layout_col);
    CHKERR PetscLayoutSetBlockSize(layout_col, 1);
    CHKERR PetscLayoutSetSize(layout_col, nb_dofs_col);
    CHKERR PetscLayoutSetUp(layout_col);
    CHKERR PetscLayoutGetRanges(layout_col, &ranges_col);
  }
  for (unsigned int part = 0; part < (unsigned int)m_field.get_comm_size();
       part++) {
    miit_row = dofs_row_by_idx.lower_bound(ranges_row[part]);
    hi_miit_row = dofs_row_by_idx.lower_bound(ranges_row[part + 1]);
    if (std::distance(miit_row, hi_miit_row) !=
        ranges_row[part + 1] - ranges_row[part]) {
      SETERRQ4(
          PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
          "data inconsistency, std::distance(miit_row,hi_miit_row) != rend - "
          "rstart (%d != %d - %d = %d) ",
          std::distance(miit_row, hi_miit_row), ranges_row[part + 1],
          ranges_row[part], ranges_row[part + 1] - ranges_row[part]);
    }
    // loop rows
    for (; miit_row != hi_miit_row; miit_row++) {
      bool success = dofs_row_by_idx.modify(
          miit_row,
          NumeredDofEntity_part_and_glob_idx_change(part, (*miit_row)->dofIdx));
      if (!success)
        SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                "modification unsuccessful");
      if (part == (unsigned int)m_field.get_comm_rank()) {
        success = dofs_row_by_idx.modify(
            miit_row, NumeredDofEntity_local_idx_change(nb_row_local_dofs++));
        if (!success)
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
      }
    }
    if (!square_matrix) {
      miit_col = dofs_col_by_idx.lower_bound(ranges_col[part]);
      hi_miit_col = dofs_col_by_idx.lower_bound(ranges_col[part + 1]);
      if (std::distance(miit_col, hi_miit_col) !=
          ranges_col[part + 1] - ranges_col[part]) {
        SETERRQ4(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
                 "data inconsistency, std::distance(miit_col,hi_miit_col) != "
                 "rend - "
                 "rstart (%d != %d - %d = %d) ",
                 std::distance(miit_col, hi_miit_col), ranges_col[part + 1],
                 ranges_col[part], ranges_col[part + 1] - ranges_col[part]);
      }
      // loop cols
      for (; miit_col != hi_miit_col; miit_col++) {
        bool success = dofs_col_by_idx.modify(
            miit_col, NumeredDofEntity_part_and_glob_idx_change(
                          part, (*miit_col)->dofIdx));
        if (!success)
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
        if (part == (unsigned int)m_field.get_comm_rank()) {
          success = dofs_col_by_idx.modify(
              miit_col, NumeredDofEntity_local_idx_change(nb_col_local_dofs++));
          if (!success)
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
        }
      }
    }
  }
  CHKERR PetscLayoutDestroy(&layout_row);
  if (!square_matrix) {
    CHKERR PetscLayoutDestroy(&layout_col);
  }
  if (square_matrix) {
    nb_col_local_dofs = nb_row_local_dofs;
    nb_col_ghost_dofs = nb_row_ghost_dofs;
  }
  CHKERR printPartitionedProblem(&*p_miit, verb);
  cOre.getBuildMoFEM() |= Core::PARTITION_PROBLEM;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::partitionProblem(const std::string name,
                                                 int verb) {
  MoFEM::Interface &m_field = cOre;
  auto problems_ptr = m_field.get_problems();
  ProblemManagerFunctionBegin;

  MOFEM_LOG("WORLD", Sev::noisy) << "Partition problem " << name;

  using NumeredDofEntitysByIdx =
      NumeredDofEntity_multiIndex::index<Idx_mi_tag>::type;
  using ProblemsByName = Problem_multiIndex::index<Problem_mi_tag>::type;

  // Find problem pointer by name
  auto &problems_set =
      const_cast<ProblemsByName &>(problems_ptr->get<Problem_mi_tag>());
  auto p_miit = problems_set.find(name);
  if (p_miit == problems_set.end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
             "problem with name %s not defined (top tip check spelling)",
             name.c_str());
  }
  int nb_dofs_row = p_miit->getNbDofsRow();

  if (m_field.get_comm_size() != 1) {

    if (!(cOre.getBuildMoFEM() & (1 << 0)))
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "fields not build");
    if (!(cOre.getBuildMoFEM() & (1 << 1)))
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "FEs not build");
    if (!(cOre.getBuildMoFEM() & (1 << 2)))
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "entFEAdjacencies not build");
    if (!(cOre.getBuildMoFEM() & (1 << 3)))
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Problems not build");

    Mat Adj;
    CHKERR m_field.getInterface<MatrixManager>()
        ->createMPIAdjWithArrays<Idx_mi_tag>(name, &Adj, verb);

    int m, n;
    CHKERR MatGetSize(Adj, &m, &n);
    if (verb > VERY_VERBOSE)
      MatView(Adj, PETSC_VIEWER_STDOUT_WORLD);

    // partitioning
    MatPartitioning part;
    IS is;
    CHKERR MatPartitioningCreate(m_field.get_comm(), &part);
    CHKERR MatPartitioningSetAdjacency(part, Adj);
    CHKERR MatPartitioningSetFromOptions(part);
    CHKERR MatPartitioningSetNParts(part, m_field.get_comm_size());
    CHKERR MatPartitioningApply(part, &is);
    if (verb > VERY_VERBOSE)
      ISView(is, PETSC_VIEWER_STDOUT_WORLD);

    // gather
    IS is_gather, is_num, is_gather_num;
    CHKERR ISAllGather(is, &is_gather);
    CHKERR ISPartitioningToNumbering(is, &is_num);
    CHKERR ISAllGather(is_num, &is_gather_num);
    const int *part_number, *petsc_idx;
    CHKERR ISGetIndices(is_gather, &part_number);
    CHKERR ISGetIndices(is_gather_num, &petsc_idx);
    int size_is_num, size_is_gather;
    CHKERR ISGetSize(is_gather, &size_is_gather);
    if (size_is_gather != (int)nb_dofs_row)
      SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
               "data inconsistency %d != %d", size_is_gather, nb_dofs_row);

    CHKERR ISGetSize(is_num, &size_is_num);
    if (size_is_num != (int)nb_dofs_row)
      SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
               "data inconsistency %d != %d", size_is_num, nb_dofs_row);

    bool square_matrix = false;
    if (p_miit->numeredRowDofsPtr == p_miit->numeredColDofsPtr)
      square_matrix = true;

    // if (!square_matrix) {
    //   // FIXME: This is for back compatibility, if deprecate interface
    //   function
    //   // build interfaces is removed, this part of the code will be obsolete
    //   auto mit_row = p_miit->numeredRowDofsPtr->get<Idx_mi_tag>().begin();
    //   auto hi_mit_row = p_miit->numeredRowDofsPtr->get<Idx_mi_tag>().end();
    //   auto mit_col = p_miit->numeredColDofsPtr->get<Idx_mi_tag>().begin();
    //   auto hi_mit_col = p_miit->numeredColDofsPtr->get<Idx_mi_tag>().end();
    //   for (; mit_row != hi_mit_row; mit_row++, mit_col++) {
    //     if (mit_col == hi_mit_col) {
    //       SETERRQ(
    //           PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
    //           "check finite element definition, nb. of rows is not equal to "
    //           "number for columns");
    //     }
    //     if (mit_row->get()->getLocalUniqueId() !=
    //         mit_col->get()->getLocalUniqueId()) {
    //       SETERRQ(
    //           PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
    //           "check finite element definition, nb. of rows is not equal to "
    //           "number for columns");
    //     }
    //   }
    // }

    auto number_dofs = [&](auto &dofs_idx, auto &counter) {
      MoFEMFunctionBegin;
      for (auto miit_dofs_row = dofs_idx.begin();
           miit_dofs_row != dofs_idx.end(); miit_dofs_row++) {
        const int part = part_number[(*miit_dofs_row)->dofIdx];
        if (part == (unsigned int)m_field.get_comm_rank()) {
          const bool success = dofs_idx.modify(
              miit_dofs_row,
              NumeredDofEntity_part_and_indices_change(
                  part, petsc_idx[(*miit_dofs_row)->dofIdx], counter++));
          if (!success) {
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
          }
        } else {
          const bool success = dofs_idx.modify(
              miit_dofs_row, NumeredDofEntity_part_and_glob_idx_change(
                                 part, petsc_idx[(*miit_dofs_row)->dofIdx]));
          if (!success) {
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
          }
        }
      }
      MoFEMFunctionReturn(0);
    };

    // Set petsc global indices
    auto &dofs_row_by_idx_no_const = const_cast<NumeredDofEntitysByIdx &>(
        p_miit->numeredRowDofsPtr->get<Idx_mi_tag>());
    int &nb_row_local_dofs = p_miit->nbLocDofsRow;
    int &nb_row_ghost_dofs = p_miit->nbGhostDofsRow;
    nb_row_local_dofs = 0;
    nb_row_ghost_dofs = 0;

    CHKERR number_dofs(dofs_row_by_idx_no_const, nb_row_local_dofs);

    int &nb_col_local_dofs = p_miit->nbLocDofsCol;
    int &nb_col_ghost_dofs = p_miit->nbGhostDofsCol;
    if (square_matrix) {
      nb_col_local_dofs = nb_row_local_dofs;
      nb_col_ghost_dofs = nb_row_ghost_dofs;
    } else {
      NumeredDofEntitysByIdx &dofs_col_by_idx_no_const =
          const_cast<NumeredDofEntitysByIdx &>(
              p_miit->numeredColDofsPtr->get<Idx_mi_tag>());
      nb_col_local_dofs = 0;
      nb_col_ghost_dofs = 0;
      CHKERR number_dofs(dofs_col_by_idx_no_const, nb_col_local_dofs);
    }

    CHKERR ISRestoreIndices(is_gather, &part_number);
    CHKERR ISRestoreIndices(is_gather_num, &petsc_idx);
    CHKERR ISDestroy(&is_num);
    CHKERR ISDestroy(&is_gather_num);
    CHKERR ISDestroy(&is_gather);
    CHKERR ISDestroy(&is);
    CHKERR MatPartitioningDestroy(&part);
    CHKERR MatDestroy(&Adj);
    CHKERR printPartitionedProblem(&*p_miit, verb);
  } else {

    auto number_dofs = [&](auto &dof_idx, auto &counter) {
      MoFEMFunctionBeginHot;
      for (auto miit_dofs_row = dof_idx.begin(); miit_dofs_row != dof_idx.end();
           miit_dofs_row++) {
        const bool success = dof_idx.modify(
            miit_dofs_row,
            NumeredDofEntity_part_and_indices_change(0, counter, counter));
        ++counter;
        if (!success) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
        }
      }
      MoFEMFunctionReturnHot(0);
    };

    auto &dofs_row_by_idx_no_const = const_cast<NumeredDofEntitysByIdx &>(
        p_miit->numeredRowDofsPtr->get<Idx_mi_tag>());
    int &nb_row_local_dofs = p_miit->nbLocDofsRow;
    int &nb_row_ghost_dofs = p_miit->nbGhostDofsRow;
    nb_row_local_dofs = 0;
    nb_row_ghost_dofs = 0;

    CHKERR number_dofs(dofs_row_by_idx_no_const, nb_row_local_dofs);

    bool square_matrix = false;
    if (p_miit->numeredRowDofsPtr == p_miit->numeredColDofsPtr)
      square_matrix = true;

    int &nb_col_local_dofs = p_miit->nbLocDofsCol;
    int &nb_col_ghost_dofs = p_miit->nbGhostDofsCol;
    if (square_matrix) {
      nb_col_local_dofs = nb_row_local_dofs;
      nb_col_ghost_dofs = nb_row_ghost_dofs;
    } else {
      NumeredDofEntitysByIdx &dofs_col_by_idx_no_const =
          const_cast<NumeredDofEntitysByIdx &>(
              p_miit->numeredColDofsPtr->get<Idx_mi_tag>());
      nb_col_local_dofs = 0;
      nb_col_ghost_dofs = 0;
      CHKERR number_dofs(dofs_col_by_idx_no_const, nb_col_local_dofs);
    }
  }

  cOre.getBuildMoFEM() |= 1 << 4;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::inheritPartition(
    const std::string name, const std::string problem_for_rows, bool copy_rows,
    const std::string problem_for_cols, bool copy_cols, int verb) {
  MoFEM::Interface &m_field = cOre;
  auto problems_ptr = m_field.get_problems();
  ProblemManagerFunctionBegin;

  if (!(cOre.getBuildMoFEM() & Core::BUILD_PROBLEM))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY, "pRoblems not build");

  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemByName;

  // find p_miit
  ProblemByName &problems_by_name =
      const_cast<ProblemByName &>(problems_ptr->get<Problem_mi_tag>());
  ProblemByName::iterator p_miit = problems_by_name.find(name);
  if (p_miit == problems_by_name.end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_NOT_FOUND,
             "problem with name < %s > not defined (top tip check spelling)",
             name.c_str());
  }
  if (verb > QUIET)
    MOFEM_LOG("WORLD", Sev::inform)
        << p_miit->getName() << " from rows of " << problem_for_rows
        << " and columns of " << problem_for_cols;

  // find p_miit_row
  ProblemByName::iterator p_miit_row = problems_by_name.find(problem_for_rows);
  if (p_miit_row == problems_by_name.end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "problem with name < %s > not defined (top tip check spelling)",
             problem_for_rows.c_str());
  }
  const boost::shared_ptr<NumeredDofEntity_multiIndex> dofs_row =
      p_miit_row->numeredRowDofsPtr;

  // find p_mit_col
  ProblemByName::iterator p_miit_col = problems_by_name.find(problem_for_cols);
  if (p_miit_col == problems_by_name.end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "problem with name < %s > not defined (top tip check spelling)",
             problem_for_cols.c_str());
  }
  boost::shared_ptr<NumeredDofEntity_multiIndex> dofs_col =
      p_miit_col->numeredColDofsPtr;

  bool copy[] = {copy_rows, copy_cols};
  boost::shared_ptr<NumeredDofEntity_multiIndex> composed_dofs[] = {
      p_miit->numeredRowDofsPtr, p_miit->numeredColDofsPtr};

  int *nb_local_dofs[] = {&p_miit->nbLocDofsRow, &p_miit->nbLocDofsCol};
  int *nb_dofs[] = {&p_miit->nbDofsRow, &p_miit->nbDofsCol};
  boost::shared_ptr<NumeredDofEntity_multiIndex> copied_dofs[] = {dofs_row,
                                                                  dofs_col};

  for (int ss = 0; ss < 2; ss++) {

    // build indices
    *nb_local_dofs[ss] = 0;
    if (!copy[ss]) {

      // only copy indices which are belong to some elements if this problem
      std::vector<int> is_local, is_new;

      NumeredDofEntityByUId &dofs_by_uid =
          copied_dofs[ss]->get<Unique_mi_tag>();
      for (NumeredDofEntity_multiIndex::iterator dit =
               composed_dofs[ss]->begin();
           dit != composed_dofs[ss]->end(); dit++) {
        NumeredDofEntityByUId::iterator diit =
            dofs_by_uid.find((*dit)->getLocalUniqueId());
        if (diit == dofs_by_uid.end()) {
          SETERRQ(
              m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "data inconsistency, could not find dof in composite problem");
        }
        int part_number = (*diit)->getPart(); // get part number
        int petsc_global_dof = (*diit)->getPetscGlobalDofIdx();
        bool success;
        success = composed_dofs[ss]->modify(
            dit, NumeredDofEntity_part_and_glob_idx_change(part_number,
                                                           petsc_global_dof));
        if (!success) {
          SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
        }
        if ((*dit)->getPart() == (unsigned int)m_field.get_comm_rank()) {
          success = composed_dofs[ss]->modify(
              dit, NumeredDofEntity_local_idx_change((*nb_local_dofs[ss])++));
          if (!success) {
            SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
          }
          is_local.push_back(petsc_global_dof);
        }
      }

      AO ao;
      CHKERR AOCreateMapping(m_field.get_comm(), is_local.size(), &is_local[0],
                             NULL, &ao);

      // apply local to global mapping
      is_local.resize(0);
      for (NumeredDofEntity_multiIndex::iterator dit =
               composed_dofs[ss]->begin();
           dit != composed_dofs[ss]->end(); dit++) {
        is_local.push_back((*dit)->getPetscGlobalDofIdx());
      }
      CHKERR AOPetscToApplication(ao, is_local.size(), &is_local[0]);
      int idx2 = 0;
      for (NumeredDofEntity_multiIndex::iterator dit =
               composed_dofs[ss]->begin();
           dit != composed_dofs[ss]->end(); dit++) {
        int part_number = (*dit)->getPart(); // get part number
        int petsc_global_dof = is_local[idx2++];
        bool success;
        success = composed_dofs[ss]->modify(
            dit, NumeredDofEntity_part_and_glob_idx_change(part_number,
                                                           petsc_global_dof));
        if (!success) {
          SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
        }
      }

      CHKERR AODestroy(&ao);

    } else {

      for (NumeredDofEntity_multiIndex::iterator dit = copied_dofs[ss]->begin();
           dit != copied_dofs[ss]->end(); dit++) {
        std::pair<NumeredDofEntity_multiIndex::iterator, bool> p;
        p = composed_dofs[ss]->insert(boost::shared_ptr<NumeredDofEntity>(
            new NumeredDofEntity((*dit)->getDofEntityPtr())));
        if (p.second) {
          (*nb_dofs[ss])++;
        }
        int dof_idx = (*dit)->getDofIdx();
        int part_number = (*dit)->getPart(); // get part number
        int petsc_global_dof = (*dit)->getPetscGlobalDofIdx();
        if (part_number == (unsigned int)m_field.get_comm_rank()) {
          const bool success = composed_dofs[ss]->modify(
              p.first, NumeredDofEntity_part_and_all_indices_change(
                           part_number, dof_idx, petsc_global_dof,
                           (*nb_local_dofs[ss])++));
          if (!success) {
            SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
          }
        } else {
          const bool success = composed_dofs[ss]->modify(
              p.first, NumeredDofEntity_part_and_mofem_glob_idx_change(
                           part_number, dof_idx, petsc_global_dof));
          if (!success) {
            SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
          }
        }
      }
    }
  }

  CHKERR printPartitionedProblem(&*p_miit, verb);
  CHKERR debugPartitionedProblem(&*p_miit, verb);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ProblemsManager::printPartitionedProblem(const Problem *problem_ptr, int verb) {
  MoFEM::Interface &m_field = cOre;
  ProblemManagerFunctionBegin;

  if (verb > QUIET) {

    MOFEM_LOG("SYNC", Sev::inform)
        << problem_ptr->getName() << " Nb. local dof "
        << problem_ptr->getNbLocalDofsRow() << " by "
        << problem_ptr->getNbLocalDofsCol() << " nb global dofs "
        << problem_ptr->getNbDofsRow() << " by " << problem_ptr->getNbDofsCol();

    MOFEM_LOG_SYNCHRONISE(m_field.get_comm())
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ProblemsManager::debugPartitionedProblem(const Problem *problem_ptr, int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBeginHot;

  auto save_ent = [](moab::Interface &moab, const std::string name,
                       const EntityHandle ent) {
    MoFEMFunctionBegin;
    EntityHandle out_meshset;
    CHKERR moab.create_meshset(MESHSET_SET, out_meshset);
    CHKERR moab.add_entities(out_meshset, &ent, 1);
    CHKERR moab.write_file(name.c_str(), "VTK", "", &out_meshset, 1);
    CHKERR moab.delete_entities(&out_meshset, 1);
    MoFEMFunctionReturn(0);
  };

  if (debug > 0) {

    typedef NumeredDofEntity_multiIndex::index<Idx_mi_tag>::type
        NumeredDofEntitysByIdx;
    NumeredDofEntitysByIdx::iterator dit, hi_dit;
    const NumeredDofEntitysByIdx *numered_dofs_ptr[] = {
        &(problem_ptr->numeredRowDofsPtr->get<Idx_mi_tag>()),
        &(problem_ptr->numeredColDofsPtr->get<Idx_mi_tag>())};

    int *nbdof_ptr[] = {&problem_ptr->nbDofsRow, &problem_ptr->nbDofsCol};
    int *local_nbdof_ptr[] = {&problem_ptr->nbLocDofsRow,
                              &problem_ptr->nbLocDofsCol};

    for (int ss = 0; ss < 2; ss++) {

      dit = numered_dofs_ptr[ss]->begin();
      hi_dit = numered_dofs_ptr[ss]->end();
      for (; dit != hi_dit; dit++) {
        if ((*dit)->getPart() == (unsigned int)m_field.get_comm_rank()) {
          if ((*dit)->getPetscLocalDofIdx() < 0) {
            std::ostringstream zz;
            zz << "rank " << m_field.get_comm_rank() << " " << **dit;
            SETERRQ2(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE,
                     "local dof index for %d (0-row, 1-col) not set, i.e. has "
                     "negative value\n %s",
                     ss, zz.str().c_str());
          }
          if ((*dit)->getPetscLocalDofIdx() >= *local_nbdof_ptr[ss]) {
            std::ostringstream zz;
            zz << "rank " << m_field.get_comm_rank() << " " << **dit;
            SETERRQ2(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE,
                     "local dofs for %d (0-row, 1-col) out of range\n %s", ss,
                     zz.str().c_str());
          }
        } else {
          if ((*dit)->getPetscGlobalDofIdx() < 0) {

            const EntityHandle ent = (*dit)->getEnt();
            CHKERR save_ent(
                m_field.get_moab(),
                "debug_part" +
                    boost::lexical_cast<std::string>(m_field.get_comm_rank()) +
                    "_negative_global_index.vtk",
                ent);

            std::ostringstream zz;
            zz << "rank " << m_field.get_comm_rank() << " "
               << dit->get()->getBitRefLevel() << " " << **dit;
            SETERRQ2(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE,
                     "global dof index for %d (0-row, 1-col) row not set, i.e. "
                     "has negative value\n %s",
                     ss, zz.str().c_str());
          }
          if ((*dit)->getPetscGlobalDofIdx() >= *nbdof_ptr[ss]) {
            std::ostringstream zz;
            zz << "rank " << m_field.get_comm_rank() << " nb_dofs "
               << *nbdof_ptr[ss] << " " << **dit;
            SETERRQ2(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE,
                     "global dofs for %d (0-row, 1-col) out of range\n %s", ss,
                     zz.str().c_str());
          }
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ProblemsManager::partitionFiniteElements(const std::string name,
                                                        bool part_from_moab,
                                                        int low_proc,
                                                        int hi_proc, int verb) {
  MoFEM::Interface &m_field = cOre;
  auto problems_ptr = m_field.get_problems();
  auto fe_ent_ptr = m_field.get_ents_finite_elements();
  ProblemManagerFunctionBegin;

  if (!(cOre.getBuildMoFEM() & Core::BUILD_FIELD))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY, "fields not build");

  if (!(cOre.getBuildMoFEM() & Core::BUILD_FE))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY, "FEs not build");

  if (!(cOre.getBuildMoFEM() & Core::BUILD_ADJ))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "adjacencies not build");

  if (!(cOre.getBuildMoFEM() & Core::BUILD_PROBLEM))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY, "problem not build");

  if (!(cOre.getBuildMoFEM() & Core::PARTITION_PROBLEM))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "problem not partitioned");

  if (low_proc == -1)
    low_proc = m_field.get_comm_rank();
  if (hi_proc == -1)
    hi_proc = m_field.get_comm_rank();

  // Find pointer to problem of given name
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemByName;
  auto &problems =
      const_cast<ProblemByName &>(problems_ptr->get<Problem_mi_tag>());
  ProblemByName::iterator p_miit = problems.find(name);
  if (p_miit == problems.end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_NOT_FOUND,
             "problem < %s > not found (top tip: check spelling)",
             name.c_str());
  }

  // Get reference on finite elements multi-index on the problem
  NumeredEntFiniteElement_multiIndex &problem_finite_elements =
      *p_miit->numeredFiniteElementsPtr;

  // Clear all elements and data, build it again
  problem_finite_elements.clear();

  // Check if dofs and columns are the same, i.e. structurally symmetric
  // problem
  bool do_cols_prob = true;
  if (p_miit->numeredRowDofsPtr == p_miit->numeredColDofsPtr) {
    do_cols_prob = false;
  }

  auto get_good_elems = [&]() {
    auto good_elems = std::vector<decltype(fe_ent_ptr->begin())>();
    good_elems.reserve(fe_ent_ptr->size());

    const auto prb_bit = p_miit->getBitRefLevel();
    const auto prb_mask = p_miit->getBitRefLevelMask();

    // Loop over all elements in database and if right element is there add it
    // to problem finite element multi-index
    for (auto efit = fe_ent_ptr->begin(); efit != fe_ent_ptr->end(); ++efit) {

      // if element is not part of problem
      if (((*efit)->getId() & p_miit->getBitFEId()).any()) {

        const auto &fe_bit = (*efit)->getBitRefLevel();

        // if entity is not problem refinement level
        if ((fe_bit & prb_mask) == fe_bit && (fe_bit & prb_bit).any())
          good_elems.emplace_back(efit);
      }
    }

    return good_elems;
  };

  auto good_elems = get_good_elems();

  auto numbered_good_elems_ptr =
      boost::make_shared<std::vector<NumeredEntFiniteElement>>();
  numbered_good_elems_ptr->reserve(good_elems.size());
  for (auto &efit : good_elems)
    numbered_good_elems_ptr->emplace_back(NumeredEntFiniteElement(*efit));

  if (!do_cols_prob) {
    for (auto &fe : *numbered_good_elems_ptr) {
      if (fe.sPtr->getRowFieldEntsPtr() == fe.sPtr->getColFieldEntsPtr()) {
        fe.getColFieldEntsPtr() = fe.getRowFieldEntsPtr();
      }
    }
  }

  if (part_from_moab) {
    for (auto &fe : *numbered_good_elems_ptr) {
      // if partition is taken from moab partition
      int proc = fe.getPartProc();
      if (proc == -1 && fe.getEntType() == MBVERTEX)
        proc = fe.getOwnerProc();
      fe.part = proc;
    }
  }

  for (auto &fe : *numbered_good_elems_ptr) {

    NumeredDofEntity_multiIndex_uid_view_ordered rows_view;
    CHKERR fe.sPtr->getRowDofView(*(p_miit->numeredRowDofsPtr), rows_view);

    if (!part_from_moab) {
      std::vector<int> parts(m_field.get_comm_size(), 0);
      for (auto &dof_ptr : rows_view)
        parts[dof_ptr->pArt]++;
      std::vector<int>::iterator pos = max_element(parts.begin(), parts.end());
      const auto max_part = std::distance(parts.begin(), pos);
      fe.part = max_part;
    }
  }

  for (auto &fe : *numbered_good_elems_ptr) {

    auto check_fields_and_dofs = [&]() {
      if (!part_from_moab) {
        if (fe.getBitFieldIdRow().none() && m_field.get_comm_size() == 0) {
          MOFEM_LOG("WORLD", Sev::warning)
              << "At least one field has to be added to element row to "
                 "determine partition of finite element. Check element " +
                     boost::lexical_cast<std::string>(fe.getName());
        }
      }

      return true;
    };

    if (check_fields_and_dofs()) {
      // Add element to the problem
      auto p = problem_finite_elements.insert(
          boost::shared_ptr<NumeredEntFiniteElement>(numbered_good_elems_ptr,
                                                     &fe));
      if (!p.second)
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "element is there");
    }
  }

  if (verb >= VERBOSE) {
    auto elements_on_rank =
        problem_finite_elements.get<Part_mi_tag>().equal_range(
            m_field.get_comm_rank());
    MOFEM_LOG("SYNC", Sev::verbose)
        << p_miit->getName() << " nb. elems "
        << std::distance(elements_on_rank.first, elements_on_rank.second);
    auto fe_ptr = m_field.get_finite_elements();
    for (auto &fe : *fe_ptr) {
      auto e_range =
          problem_finite_elements.get<Composite_Name_And_Part_mi_tag>()
              .equal_range(
                  boost::make_tuple(fe->getName(), m_field.get_comm_rank()));
      MOFEM_LOG("SYNC", Sev::noisy)
          << "Element " << fe->getName() << " nb. elems "
          << std::distance(e_range.first, e_range.second);
    }

    MOFEM_LOG_SYNCHRONISE(m_field.get_comm());
  }

  cOre.getBuildMoFEM() |= Core::PARTITION_FE;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::partitionGhostDofs(const std::string name,
                                                   int verb) {
  MoFEM::Interface &m_field = cOre;
  auto problems_ptr = m_field.get_problems();
  ProblemManagerFunctionBegin;

  if (!(cOre.getBuildMoFEM() & Core::PARTITION_PROBLEM))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "partition of problem not build");
  if (!(cOre.getBuildMoFEM() & Core::PARTITION_FE))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "partitions finite elements not build");

  // get problem pointer
  auto p_miit = problems_ptr->get<Problem_mi_tag>().find(name);
  if (p_miit == problems_ptr->get<Problem_mi_tag>().end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "Problem %s not fond",
             name.c_str());

  // get reference to number of ghost dofs
  int &nb_row_ghost_dofs = p_miit->nbGhostDofsRow;
  int &nb_col_ghost_dofs = p_miit->nbGhostDofsCol;
  nb_row_ghost_dofs = 0;
  nb_col_ghost_dofs = 0;

  // do work if more than one processor
  if (m_field.get_comm_size() > 1) {

    NumeredDofEntity_multiIndex_uid_view_ordered ghost_idx_row_view,
        ghost_idx_col_view;

    // get elements on this partition
    auto fe_range =
        p_miit->numeredFiniteElementsPtr->get<Part_mi_tag>().equal_range(
            m_field.get_comm_rank());

    // get dofs on elements which are not part of this partition

    struct Inserter {
      using Vec = std::vector<boost::shared_ptr<NumeredDofEntity>>;
      using It = Vec::iterator;
      It operator()(Vec &dofs_view, It &hint,
                    boost::shared_ptr<NumeredDofEntity> &&dof) {
        dofs_view.emplace_back(dof);
        return dofs_view.end();
      }
    };

    // rows
    std::vector<boost::shared_ptr<NumeredDofEntity>> fe_vec_view;
    auto hint_r = ghost_idx_row_view.begin();
    for (auto fe_ptr = fe_range.first; fe_ptr != fe_range.second; ++fe_ptr) {

      fe_vec_view.clear();
      CHKERR EntFiniteElement::getDofView((*fe_ptr)->getRowFieldEnts(),
                                          *(p_miit->getNumeredRowDofsPtr()),
                                          fe_vec_view, Inserter());

      for (auto &dof_ptr : fe_vec_view) {
        if (dof_ptr->getPart() != (unsigned int)m_field.get_comm_rank()) {
          hint_r = ghost_idx_row_view.emplace_hint(hint_r, dof_ptr);
        }
      }
    }

    // columns
    if (p_miit->numeredColDofsPtr == p_miit->numeredRowDofsPtr) {

      auto hint_c = ghost_idx_col_view.begin();
      for (auto fe_ptr = fe_range.first; fe_ptr != fe_range.second; ++fe_ptr) {

        fe_vec_view.clear();
        CHKERR EntFiniteElement::getDofView((*fe_ptr)->getColFieldEnts(),
                                            *(p_miit->getNumeredColDofsPtr()),
                                            fe_vec_view, Inserter());

        for (auto &dof_ptr : fe_vec_view) {
          if (dof_ptr->getPart() != (unsigned int)m_field.get_comm_rank()) {
            hint_c = ghost_idx_col_view.emplace_hint(hint_c, dof_ptr);
          }
        }
      }
    }

    int *nb_ghost_dofs[2] = {&nb_row_ghost_dofs, &nb_col_ghost_dofs};
    int nb_local_dofs[2] = {p_miit->nbLocDofsRow, p_miit->nbLocDofsCol};

    NumeredDofEntity_multiIndex_uid_view_ordered *ghost_idx_view[2] = {
        &ghost_idx_row_view, &ghost_idx_col_view};
    NumeredDofEntityByUId *dof_by_uid_no_const[2] = {
        &p_miit->numeredRowDofsPtr->get<Unique_mi_tag>(),
        &p_miit->numeredColDofsPtr->get<Unique_mi_tag>()};

    int loop_size = 2;
    if (p_miit->numeredColDofsPtr == p_miit->numeredRowDofsPtr) {
      loop_size = 1;
    }

    // set local ghost dofs indices
    for (int ss = 0; ss != loop_size; ++ss) {
      for (auto &gid : *ghost_idx_view[ss]) {
        NumeredDofEntityByUId::iterator dof =
            dof_by_uid_no_const[ss]->find(gid->getLocalUniqueId());
        if (PetscUnlikely((*dof)->petscLocalDofIdx != (DofIdx)-1))
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "inconsistent data, ghost dof already set");
        bool success = dof_by_uid_no_const[ss]->modify(
            dof, NumeredDofEntity_local_idx_change(nb_local_dofs[ss]++));
        if (PetscUnlikely(!success))
          SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
        (*nb_ghost_dofs[ss])++;
      }
    }
    if (loop_size == 1) {
      (*nb_ghost_dofs[1]) = (*nb_ghost_dofs[0]);
    }
  }

  if (verb > QUIET) {
    MOFEM_LOG("SYNC", Sev::inform)
        << " FEs ghost dofs on problem " << p_miit->getName()
        << " Nb. ghost dof " << p_miit->getNbGhostDofsRow() << " by "
        << p_miit->getNbGhostDofsCol() << " Nb. local dof "
        << p_miit->getNbLocalDofsCol() << " by " << p_miit->getNbLocalDofsCol();

    MOFEM_LOG_SYNCHRONISE(m_field.get_comm())
  }

  cOre.getBuildMoFEM() |= Core::PARTITION_GHOST_DOFS;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ProblemsManager::partitionGhostDofsOnDistributedMesh(const std::string name,
                                                     int verb) {
  MoFEM::Interface &m_field = cOre;
  auto problems_ptr = m_field.get_problems();
  ProblemManagerFunctionBegin;

  if (!(cOre.getBuildMoFEM() & Core::PARTITION_PROBLEM))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "partition of problem not build");
  if (!(cOre.getBuildMoFEM() & Core::PARTITION_FE))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "partitions finite elements not build");

  // get problem pointer
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  // find p_miit
  ProblemsByName &problems_set =
      const_cast<ProblemsByName &>(problems_ptr->get<Problem_mi_tag>());
  ProblemsByName::iterator p_miit = problems_set.find(name);

  // get reference to number of ghost dofs
  // get number of local dofs
  DofIdx *nb_ghost_dofs[2] = {&(p_miit->nbGhostDofsRow),
                              &(p_miit->nbGhostDofsCol)};
  DofIdx nb_local_dofs[2] = {p_miit->nbLocDofsRow, p_miit->nbLocDofsCol};
  for (int ss = 0; ss != 2; ++ss) {
    (*nb_ghost_dofs[ss]) = 0;
  }

  // do work if more than one processor
  if (m_field.get_comm_size() > 1) {
    // determine if rows on columns are different from dofs on rows
    int loop_size = 2;
    if (p_miit->numeredColDofsPtr == p_miit->numeredRowDofsPtr) {
      loop_size = 1;
    }

    typedef decltype(p_miit->numeredRowDofsPtr) NumbDofTypeSharedPtr;
    NumbDofTypeSharedPtr numered_dofs[] = {p_miit->numeredRowDofsPtr,
                                           p_miit->numeredColDofsPtr};

    // iterate over dofs on rows and dofs on columns
    for (int ss = 0; ss != loop_size; ++ss) {

      // create dofs view by uid
      auto r = numered_dofs[ss]->get<PetscLocalIdx_mi_tag>().equal_range(-1);

      std::vector<NumeredDofEntity_multiIndex::iterator> ghost_idx_view;
      ghost_idx_view.reserve(std::distance(r.first, r.second));
      for (; r.first != r.second; ++r.first)
        ghost_idx_view.emplace_back(numered_dofs[ss]->project<0>(r.first));

      auto cmp = [](auto a, auto b) {
        return (*a)->getLocalUniqueId() < (*b)->getLocalUniqueId();
      };
      sort(ghost_idx_view.begin(), ghost_idx_view.end(), cmp);

      // iterate over dofs which have negative local index
      for (auto gid_it : ghost_idx_view) {
        bool success = numered_dofs[ss]->modify(
            gid_it, NumeredDofEntity_local_idx_change((nb_local_dofs[ss])++));
        if (!success)
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
        ++(*nb_ghost_dofs[ss]);
      }
    }
    if (loop_size == 1) {
      (*nb_ghost_dofs[1]) = (*nb_ghost_dofs[0]);
    }
  }

  if (verb > QUIET) {
    MOFEM_LOG("SYNC", Sev::inform)
        << " FEs ghost dofs on problem " << p_miit->getName()
        << " Nb. ghost dof " << p_miit->getNbGhostDofsRow() << " by "
        << p_miit->getNbGhostDofsCol() << " Nb. local dof "
        << p_miit->getNbLocalDofsCol() << " by " << p_miit->getNbLocalDofsCol();

    MOFEM_LOG_SYNCHRONISE(m_field.get_comm())
  }

  cOre.getBuildMoFEM() |= Core::PARTITION_GHOST_DOFS;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::getFEMeshset(const std::string prb_name,
                                             const std::string &fe_name,
                                             EntityHandle *meshset) const {
  MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  const FiniteElement_multiIndex *fes_ptr;
  ProblemManagerFunctionBegin;

  CHKERR m_field.get_moab().create_meshset(MESHSET_SET, *meshset);
  CHKERR m_field.get_problem(prb_name, &problem_ptr);
  CHKERR m_field.get_finite_elements(&fes_ptr);

  auto fe_miit = fes_ptr->get<FiniteElement_name_mi_tag>().find(fe_name);
  if (fe_miit != fes_ptr->get<FiniteElement_name_mi_tag>().end()) {
    auto fit =
        problem_ptr->numeredFiniteElementsPtr->get<Unique_mi_tag>().lower_bound(
            EntFiniteElement::getLocalUniqueIdCalculate(
                0, (*fe_miit)->getFEUId()));
    auto hi_fe_it =
        problem_ptr->numeredFiniteElementsPtr->get<Unique_mi_tag>().upper_bound(
            EntFiniteElement::getLocalUniqueIdCalculate(
                get_id_for_max_type<MBENTITYSET>(), (*fe_miit)->getFEUId()));
    std::vector<EntityHandle> fe_vec;
    fe_vec.reserve(std::distance(fit, hi_fe_it));
    for (; fit != hi_fe_it; fit++)
      fe_vec.push_back(fit->get()->getEnt());
    CHKERR m_field.get_moab().add_entities(*meshset, &*fe_vec.begin(),
                                           fe_vec.size());
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ProblemsManager::getProblemElementsLayout(const std::string name,
                                          const std::string &fe_name,
                                          PetscLayout *layout) const {
  MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  ProblemManagerFunctionBegin;
  CHKERR m_field.get_problem(name, &problem_ptr);
  CHKERR problem_ptr->getNumberOfElementsByNameAndPart(PETSC_COMM_WORLD,
                                                       fe_name, layout);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::removeDofsOnEntities(
    const std::string problem_name, const std::string field_name,
    const Range ents, const int lo_coeff, const int hi_coeff,
    const int lo_order, const int hi_order, int verb, const bool debug) {

  MoFEM::Interface &m_field = cOre;
  ProblemManagerFunctionBegin;

  const Problem *prb_ptr;
  CHKERR m_field.get_problem(problem_name, &prb_ptr);

  decltype(prb_ptr->numeredRowDofsPtr) numered_dofs[2] = {
      prb_ptr->numeredRowDofsPtr, nullptr};
  if (prb_ptr->numeredRowDofsPtr != prb_ptr->numeredColDofsPtr)
    numered_dofs[1] = prb_ptr->numeredColDofsPtr;

  int *nbdof_ptr[] = {&prb_ptr->nbDofsRow, &prb_ptr->nbDofsCol};
  int *local_nbdof_ptr[] = {&prb_ptr->nbLocDofsRow, &prb_ptr->nbLocDofsCol};
  int *ghost_nbdof_ptr[] = {&prb_ptr->nbGhostDofsRow, &prb_ptr->nbGhostDofsCol};

  const int nb_init_row_dofs = prb_ptr->getNbDofsRow();
  const int nb_init_col_dofs = prb_ptr->getNbDofsCol();
  const int nb_init_loc_row_dofs = prb_ptr->getNbLocalDofsRow();
  const int nb_init_loc_col_dofs = prb_ptr->getNbLocalDofsCol();
  const int nb_init_ghost_row_dofs = prb_ptr->getNbGhostDofsRow();
  const int nb_init_ghost_col_dofs = prb_ptr->getNbGhostDofsCol();

  for (int s = 0; s != 2; ++s)
    if (numered_dofs[s]) {

      typedef multi_index_container<

          NumeredDofEntity_multiIndex::iterator, indexed_by<sequenced<>>

          >
          NumeredDofEntity_it_view_multiIndex;

      const auto bit_number = m_field.get_field_bit_number(field_name);
      NumeredDofEntity_it_view_multiIndex dofs_it_view;

      // Set -1 to global and local dofs indices
      for (auto pit = ents.const_pair_begin(); pit != ents.const_pair_end();
           ++pit) {
        auto lo = numered_dofs[s]->get<Unique_mi_tag>().lower_bound(
            DofEntity::getLoFieldEntityUId(bit_number, pit->first));
        auto hi = numered_dofs[s]->get<Unique_mi_tag>().upper_bound(
            DofEntity::getHiFieldEntityUId(bit_number, pit->second));

        for (; lo != hi; ++lo)
          if ((*lo)->getDofCoeffIdx() >= lo_coeff &&
              (*lo)->getDofCoeffIdx() <= hi_coeff &&
              (*lo)->getDofOrder() >= lo_order &&
              (*lo)->getDofOrder() <= hi_order)
            dofs_it_view.emplace_back(numered_dofs[s]->project<0>(lo));
      }

      if (verb > QUIET) {
        for (auto &dof : dofs_it_view)
          MOFEM_LOG("SYNC", Sev::noisy) << **dof;
        MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::noisy);
      }

      // create weak view
      std::vector<boost::weak_ptr<NumeredDofEntity>> dofs_weak_view;
      dofs_weak_view.reserve(dofs_it_view.size());
      for (auto dit : dofs_it_view)
        dofs_weak_view.push_back(*dit);

      if (verb >= NOISY)
        MOFEM_LOG_C("SYNC", Sev::noisy,
                    "Number of DOFs in multi-index %d and to delete %d\n",
                    numered_dofs[s]->size(), dofs_it_view.size());

      // erase dofs from problem
      for (auto weak_dit : dofs_weak_view)
        if (auto dit = weak_dit.lock()) {
          numered_dofs[s]->erase(dit->getLocalUniqueId());
        }

      if (verb >= NOISY)
        MOFEM_LOG_C("SYNC", Sev::noisy,
                    "Number of DOFs in multi-index after delete %d\n",
                    numered_dofs[s]->size());

      // get current number of ghost dofs
      int nb_local_dofs = 0;
      int nb_ghost_dofs = 0;
      for (auto dit = numered_dofs[s]->get<PetscLocalIdx_mi_tag>().begin();
           dit != numered_dofs[s]->get<PetscLocalIdx_mi_tag>().end(); ++dit) {
        if ((*dit)->getPetscLocalDofIdx() >= 0 &&
            (*dit)->getPetscLocalDofIdx() < *(local_nbdof_ptr[s]))
          ++nb_local_dofs;
        else if ((*dit)->getPetscLocalDofIdx() >= *(local_nbdof_ptr[s]))
          ++nb_ghost_dofs;
        else
          SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "Impossible case. You could run problem on no distributed "
                   "mesh. That is not implemented. Dof local index is %d",
                   (*dit)->getPetscLocalDofIdx());
      }

      // get indices
      auto get_indices_by_tag = [&](auto tag, auto &indices, bool only_local) {
        const int nb_dofs = numered_dofs[s]->size();
        indices.clear();
        indices.reserve(nb_dofs);
        for (auto dit = numered_dofs[s]->get<decltype(tag)>().begin();
             dit != numered_dofs[s]->get<decltype(tag)>().end(); ++dit) {
          bool add = true;
          if (only_local) {
            if ((*dit)->getPetscLocalDofIdx() < 0 ||
                (*dit)->getPetscLocalDofIdx() >= *(local_nbdof_ptr[s])) {
              add = false;
            }
          }
          if (add)
            indices.push_back(decltype(tag)::get_index(dit));
        }
      };

      auto get_indices_by_uid = [&](auto tag, auto &indices) {
        const int nb_dofs = numered_dofs[s]->size();
        indices.clear();
        indices.reserve(nb_dofs);
        for (auto dit = numered_dofs[s]->begin(); dit != numered_dofs[s]->end();
             ++dit)
          indices.push_back(decltype(tag)::get_index(dit));
      };

      auto get_sub_ao = [&](auto sub_data) {
        if (s == 0) {
          return sub_data->getSmartRowMap();
        } else {
          return sub_data->getSmartColMap();
        }
      };

      auto set_sub_is_and_ao = [&s, &prb_ptr](auto sub_data, auto is, auto ao) {
        if (s == 0) {
          sub_data->rowIs = is;
          sub_data->rowMap = ao;
        } else {
          sub_data->colIs = is;
          sub_data->colMap = ao;
        }
      };

      auto apply_symmetry = [&s, &prb_ptr](auto sub_data) {
        if (s == 0) {
          if (prb_ptr->numeredRowDofsPtr == prb_ptr->numeredColDofsPtr) {
            sub_data->colIs = sub_data->getSmartRowIs();
            sub_data->colMap = sub_data->getSmartRowMap();
          }
        }
      };

      auto concatenate_dofs = [&](auto tag, auto &indices,
                                  const auto local_only) {
        MoFEMFunctionBegin;
        get_indices_by_tag(tag, indices, local_only);

        SmartPetscObj<AO> ao;
        // Create AO from app indices (i.e. old), to pestc indices (new after
        // remove)
        if (local_only)
          ao = createAOMapping(m_field.get_comm(), indices.size(),
                               &*indices.begin(), PETSC_NULL);
        else
          ao = createAOMapping(PETSC_COMM_SELF, indices.size(),
                               &*indices.begin(), PETSC_NULL);

        // Set mapping to sub dm data
        if (local_only) {
          if (auto sub_data = prb_ptr->getSubData()) {
            // create is and then map it to main problem of sub-problem
            auto sub_is = createISGeneral(m_field.get_comm(), indices.size(),
                                          &*indices.begin(), PETSC_COPY_VALUES);
            // get old app, i.e. oroginal befor sub indices, and ao, from app,
            // to petsc sub indices.
            auto sub_ao = get_sub_ao(sub_data);
            CHKERR AOPetscToApplicationIS(sub_ao, sub_is);
            sub_ao = createAOMappingIS(sub_is, PETSC_NULL);
            // set new sub ao
            set_sub_is_and_ao(sub_data, sub_is, sub_ao);
            apply_symmetry(sub_data);
          } else {
            // create sub data
            prb_ptr->getSubData() =
                boost::make_shared<Problem::SubProblemData>();
            auto sub_is = createISGeneral(m_field.get_comm(), indices.size(),
                                          &*indices.begin(), PETSC_COPY_VALUES);
            // set sub is ao
            set_sub_is_and_ao(prb_ptr->getSubData(), sub_is, ao);
            apply_symmetry(prb_ptr->getSubData());
          }
        }

        get_indices_by_uid(tag, indices);
        CHKERR AOApplicationToPetsc(ao, indices.size(), &*indices.begin());

        MoFEMFunctionReturn(0);
      };

      // set indices index
      auto set_concatenated_indices = [&]() {
        std::vector<int> global_indices;
        std::vector<int> local_indices;
        MoFEMFunctionBegin;
        CHKERR concatenate_dofs(PetscGlobalIdx_mi_tag(), global_indices, true);
        CHKERR concatenate_dofs(PetscLocalIdx_mi_tag(), local_indices, false);
        auto gi = global_indices.begin();
        auto li = local_indices.begin();
        for (auto dit = numered_dofs[s]->begin(); dit != numered_dofs[s]->end();
             ++dit) {
          auto mod = NumeredDofEntity_part_and_all_indices_change(
              (*dit)->getPart(), (*dit)->getDofIdx(), *gi, *li);
          bool success = numered_dofs[s]->modify(dit, mod);
          if (!success)
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "can not set negative indices");
          ++gi;
          ++li;
        }
        MoFEMFunctionReturn(0);
      };
      CHKERR set_concatenated_indices();

      MPI_Allreduce(&nb_local_dofs, nbdof_ptr[s], 1, MPI_INT, MPI_SUM,
                    m_field.get_comm());
      *(local_nbdof_ptr[s]) = nb_local_dofs;
      *(ghost_nbdof_ptr[s]) = nb_ghost_dofs;

      if (debug)
        for (auto dof : (*numered_dofs[s])) {
          if (dof->getPetscGlobalDofIdx() < 0) {
            SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "Negative global idx");
          }
          if (dof->getPetscLocalDofIdx() < 0) {
            SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "Negative local idx");
          }
        }

    } else {

      *(nbdof_ptr[1]) = *(nbdof_ptr[0]);
      *(local_nbdof_ptr[1]) = *(local_nbdof_ptr[0]);
      *(ghost_nbdof_ptr[1]) = *(ghost_nbdof_ptr[0]);
    }

  if (verb > QUIET) {
    MOFEM_LOG_C(
        "WORLD", Sev::inform,
        "Removed DOFs from problem %s dofs [%d / %d (before %d / %d) global]",
        prb_ptr->getName().c_str(), prb_ptr->getNbDofsRow(),
        prb_ptr->getNbDofsCol(), nb_init_row_dofs, nb_init_col_dofs);
    MOFEM_LOG_C("SYNC", Sev::verbose,
                "Removed DOFs from problem %s dofs [ %d / %d  "
                "(before %d / %d) local, %d / %d (before %d / %d)]",
                prb_ptr->getName().c_str(), prb_ptr->getNbLocalDofsRow(),
                prb_ptr->getNbLocalDofsCol(), nb_init_loc_row_dofs,
                nb_init_loc_row_dofs, prb_ptr->getNbGhostDofsRow(),
                prb_ptr->getNbGhostDofsCol(), nb_init_ghost_row_dofs,
                nb_init_ghost_col_dofs);
    MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::verbose);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::removeDofsOnEntitiesNotDistributed(
    const std::string problem_name, const std::string field_name,
    const Range ents, const int lo_coeff, const int hi_coeff,
    const int lo_order, const int hi_order, int verb, const bool debug) {

  MoFEM::Interface &m_field = cOre;
  ProblemManagerFunctionBegin;

  const Problem *prb_ptr;
  CHKERR m_field.get_problem(problem_name, &prb_ptr);

  decltype(prb_ptr->numeredRowDofsPtr) numered_dofs[2] = {
      prb_ptr->numeredRowDofsPtr, nullptr};
  if (prb_ptr->numeredRowDofsPtr != prb_ptr->numeredColDofsPtr)
    numered_dofs[1] = prb_ptr->numeredColDofsPtr;

  int *nbdof_ptr[] = {&prb_ptr->nbDofsRow, &prb_ptr->nbDofsCol};
  int *local_nbdof_ptr[] = {&prb_ptr->nbLocDofsRow, &prb_ptr->nbLocDofsCol};
  int *ghost_nbdof_ptr[] = {&prb_ptr->nbGhostDofsRow, &prb_ptr->nbGhostDofsCol};

  const int nb_init_row_dofs = prb_ptr->getNbDofsRow();
  const int nb_init_col_dofs = prb_ptr->getNbDofsCol();
  const int nb_init_loc_row_dofs = prb_ptr->getNbLocalDofsRow();
  const int nb_init_loc_col_dofs = prb_ptr->getNbLocalDofsCol();
  const int nb_init_ghost_row_dofs = prb_ptr->getNbGhostDofsRow();
  const int nb_init_ghost_col_dofs = prb_ptr->getNbGhostDofsCol();

  const std::array<int, 2> nb_init_dofs = {nb_init_row_dofs, nb_init_col_dofs};

  for (int s = 0; s != 2; ++s)
    if (numered_dofs[s]) {

      typedef multi_index_container<

          NumeredDofEntity_multiIndex::iterator, indexed_by<sequenced<>>

          >
          NumeredDofEntity_it_view_multiIndex;

      const auto bit_number = m_field.get_field_bit_number(field_name);
      NumeredDofEntity_it_view_multiIndex dofs_it_view;

      // Set -1 to global and local dofs indices
      for (auto pit = ents.const_pair_begin(); pit != ents.const_pair_end();
           ++pit) {
        auto lo = numered_dofs[s]->get<Unique_mi_tag>().lower_bound(
            DofEntity::getLoFieldEntityUId(bit_number, pit->first));
        auto hi = numered_dofs[s]->get<Unique_mi_tag>().upper_bound(
            DofEntity::getHiFieldEntityUId(bit_number, pit->second));

        for (; lo != hi; ++lo)
          if ((*lo)->getDofCoeffIdx() >= lo_coeff &&
              (*lo)->getDofCoeffIdx() <= hi_coeff &&
              (*lo)->getDofOrder() >= lo_order &&
              (*lo)->getDofOrder() <= hi_order)
            dofs_it_view.emplace_back(numered_dofs[s]->project<0>(lo));
      }

      if (verb > QUIET) {
        for (auto &dof : dofs_it_view)
          MOFEM_LOG("SYNC", Sev::noisy) << **dof;
      }

      // set negative index
      auto mod = NumeredDofEntity_part_and_all_indices_change(-1, -1, -1, -1);
      for (auto dit : dofs_it_view) {
        bool success = numered_dofs[s]->modify(dit, mod);
        if (!success)
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "can not set negative indices");
      }

      // create weak view
      std::vector<boost::weak_ptr<NumeredDofEntity>> dosf_weak_view;
      dosf_weak_view.reserve(dofs_it_view.size());
      for (auto dit : dofs_it_view)
        dosf_weak_view.push_back(*dit);

      if (verb >= NOISY)
        MOFEM_LOG_C("SYNC", Sev::noisy,
                    "Number of DOFs in multi-index %d and to delete %d\n",
                    numered_dofs[s]->size(), dofs_it_view.size());

      // erase dofs from problem
      for (auto weak_dit : dosf_weak_view)
        if (auto dit = weak_dit.lock()) {
          numered_dofs[s]->erase(dit->getLocalUniqueId());
        }

      if (verb >= NOISY)
        MOFEM_LOG_C("SYNC", Sev::noisy,
                    "Number of DOFs in multi-index after delete %d\n",
                    numered_dofs[s]->size());

      // get current number of ghost dofs
      int nb_global_dof = 0;
      int nb_local_dofs = 0;
      int nb_ghost_dofs = 0;

      for (auto dit = numered_dofs[s]->begin(); dit != numered_dofs[s]->end();
           ++dit) {

        if ((*dit)->getDofIdx() >= 0) {

          if ((*dit)->getPetscLocalDofIdx() >= 0 &&
              (*dit)->getPetscLocalDofIdx() < *(local_nbdof_ptr[s]))
            ++nb_local_dofs;
          else if ((*dit)->getPetscLocalDofIdx() >= *(local_nbdof_ptr[s]))
            ++nb_ghost_dofs;

          ++nb_global_dof;
        }
      }

      if (debug) {
        MPI_Allreduce(&nb_local_dofs, nbdof_ptr[s], 1, MPI_INT, MPI_SUM,
                      m_field.get_comm());
        if (*(nbdof_ptr[s]) != nb_global_dof)
          SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "Number of local DOFs do not add up %d != %d",
                   *(nbdof_ptr[s]), nb_global_dof);
      }

      *(nbdof_ptr[s]) = nb_global_dof;
      *(local_nbdof_ptr[s]) = nb_local_dofs;
      *(ghost_nbdof_ptr[s]) = nb_ghost_dofs;

      // get indices
      auto get_indices_by_tag = [&](auto tag) {
        std::vector<int> indices;
        indices.resize(nb_init_dofs[s], -1);
        for (auto dit = numered_dofs[s]->get<Idx_mi_tag>().lower_bound(0);
             dit != numered_dofs[s]->get<Idx_mi_tag>().end(); ++dit) {
          indices[(*dit)->getDofIdx()] = decltype(tag)::get_index(dit);
        }
        return indices;
      };

      auto renumber = [&](auto tag, auto &indices) {
        MoFEMFunctionBegin;
        int idx = 0;
        for (auto dit = numered_dofs[s]->get<decltype(tag)>().lower_bound(0);
             dit != numered_dofs[s]->get<decltype(tag)>().end(); ++dit) {
          indices[(*dit)->getDofIdx()] = idx++;
        }
        MoFEMFunctionReturn(0);
      };

      auto get_sub_ao = [&](auto sub_data) {
        if (s == 0) {
          return sub_data->getSmartRowMap();
        } else {
          return sub_data->getSmartColMap();
        }
      };

      auto set_sub_is_and_ao = [&s, &prb_ptr](auto sub_data, auto is, auto ao) {
        if (s == 0) {
          sub_data->rowIs = is;
          sub_data->rowMap = ao;
        } else {
          sub_data->colIs = is;
          sub_data->colMap = ao;
        }
      };

      auto apply_symmetry = [&s, &prb_ptr](auto sub_data) {
        if (s == 0) {
          if (prb_ptr->numeredRowDofsPtr == prb_ptr->numeredColDofsPtr) {
            sub_data->colIs = sub_data->getSmartRowIs();
            sub_data->colMap = sub_data->getSmartRowMap();
          }
        }
      };

      auto set_sub_data = [&](auto &indices) {
        MoFEMFunctionBegin;
        if (auto sub_data = prb_ptr->getSubData()) {
          // create is and then map it to main problem of sub-problem
          auto sub_is = createISGeneral(m_field.get_comm(), indices.size(),
                                        &*indices.begin(), PETSC_COPY_VALUES);
          // get old app, i.e. oroginal befor sub indices, and ao, from
          // app, to petsc sub indices.
          auto sub_ao = get_sub_ao(sub_data);
          CHKERR AOPetscToApplicationIS(sub_ao, sub_is);
          sub_ao = createAOMappingIS(sub_is, PETSC_NULL);
          // set new sub ao
          set_sub_is_and_ao(sub_data, sub_is, sub_ao);
          apply_symmetry(sub_data);
        } else {
          prb_ptr->getSubData() = boost::make_shared<Problem::SubProblemData>();
          auto sub_is = createISGeneral(m_field.get_comm(), indices.size(),
                                        &*indices.begin(), PETSC_COPY_VALUES);
          auto sub_ao = createAOMappingIS(sub_is, PETSC_NULL);
          // set sub is ao
          set_sub_is_and_ao(prb_ptr->getSubData(), sub_is, sub_ao);
          apply_symmetry(prb_ptr->getSubData());
        }
        MoFEMFunctionReturn(0);
      };

      auto global_indices = get_indices_by_tag(PetscGlobalIdx_mi_tag());
      auto local_indices = get_indices_by_tag(PetscLocalIdx_mi_tag());
      CHKERR set_sub_data(global_indices);
      CHKERR renumber(PetscGlobalIdx_mi_tag(), global_indices);
      CHKERR renumber(PetscLocalIdx_mi_tag(), local_indices);

      int i = 0;    
      for (auto dit = numered_dofs[s]->begin(); dit != numered_dofs[s]->end();
           ++dit) {
        auto idx = (*dit)->getDofIdx();
        if (idx >= 0) {
          auto mod = NumeredDofEntity_part_and_all_indices_change(
              (*dit)->getPart(), i++, global_indices[idx], local_indices[idx]);
          bool success = numered_dofs[s]->modify(dit, mod);
          if (!success)
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "can not set negative indices");
        } else {
          auto mod = NumeredDofEntity_part_and_all_indices_change(
              (*dit)->getPart(), -1, -1, -1);
          bool success = numered_dofs[s]->modify(dit, mod);
          if (!success)
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "can not set negative indices");
 
        }
      };

      if (debug) {
        for (auto dof : (*numered_dofs[s])) {
          if (dof->getDofIdx() >= 0 && dof->getPetscGlobalDofIdx() < 0) {
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "Negative global idx");
          }
        }

      }

    } else {

      *(nbdof_ptr[1]) = *(nbdof_ptr[0]);
      *(local_nbdof_ptr[1]) = *(local_nbdof_ptr[0]);
      *(ghost_nbdof_ptr[1]) = *(ghost_nbdof_ptr[0]);
    }

  if (verb >= NOISY)
    MOFEM_LOG_SYNCHRONISE(m_field.get_comm());

  if (verb > QUIET) {
    MOFEM_LOG_C(
        "WORLD", Sev::inform,
        "Removed DOFs from problem %s dofs [%d / %d (before %d / %d) global]",
        prb_ptr->getName().c_str(), prb_ptr->getNbDofsRow(),
        prb_ptr->getNbDofsCol(), nb_init_row_dofs, nb_init_col_dofs);
    MOFEM_LOG_C("SYNC", Sev::verbose,
                "Removed DOFs from problem %s dofs [ %d / %d  "
                "(before %d / %d) local, %d / %d (before %d / %d)]",
                prb_ptr->getName().c_str(), prb_ptr->getNbLocalDofsRow(),
                prb_ptr->getNbLocalDofsCol(), nb_init_loc_row_dofs,
                nb_init_loc_row_dofs, prb_ptr->getNbGhostDofsRow(),
                prb_ptr->getNbGhostDofsCol(), nb_init_ghost_row_dofs,
                nb_init_ghost_col_dofs);
    MOFEM_LOG_SYNCHRONISE(m_field.get_comm());
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::removeDofsOnEntities(
    const std::string problem_name, const std::string field_name,
    const BitRefLevel bit_ref_level, const BitRefLevel bit_ref_mask,
    Range *ents_ptr, const int lo_coeff, const int hi_coeff, const int lo_order,
    const int hi_order, int verb, const bool debug) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;

  auto bit_manager = m_field.getInterface<BitRefManager>();

  Range ents;
  if (ents_ptr) {
    ents = *ents_ptr;
    CHKERR bit_manager->filterEntitiesByRefLevel(bit_ref_level, bit_ref_mask,
                                                 ents, verb);
  } else {
    CHKERR bit_manager->getEntitiesByRefLevel(bit_ref_level, bit_ref_mask, ents,
                                              verb);
  }

  CHKERR removeDofsOnEntities(problem_name, field_name, ents, lo_coeff,
                              hi_coeff, lo_order, hi_order, verb, debug);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::removeDofsOnEntitiesNotDistributed(
    const std::string problem_name, const std::string field_name,
    const BitRefLevel bit_ref_level, const BitRefLevel bit_ref_mask,
    Range *ents_ptr, const int lo_coeff, const int hi_coeff, const int lo_order,
    const int hi_order, int verb, const bool debug) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;

  auto bit_manager = m_field.getInterface<BitRefManager>();

  Range ents;
  if (ents_ptr) {
    ents = *ents_ptr;
    CHKERR bit_manager->filterEntitiesByRefLevel(bit_ref_level, bit_ref_mask,
                                                 ents, verb);
  } else {
    CHKERR bit_manager->getEntitiesByRefLevel(bit_ref_level, bit_ref_mask, ents,
                                              verb);
  }

  CHKERR removeDofsOnEntitiesNotDistributed(problem_name, field_name, ents,
                                            lo_coeff, hi_coeff, lo_order,
                                            hi_order, verb, debug);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ProblemsManager::markDofs(const std::string problem_name, RowColData rc,
                          const enum MarkOP op, const Range ents,
                          std::vector<unsigned char> &marker) const {

  Interface &m_field = cOre;
  const Problem *problem_ptr;
  ProblemManagerFunctionBegin;
  CHKERR m_field.get_problem(problem_name, &problem_ptr);
  boost::shared_ptr<NumeredDofEntity_multiIndex> dofs;
  switch (rc) {
  case ROW:
    dofs = problem_ptr->getNumeredRowDofsPtr();
    break;
  case COL:
    dofs = problem_ptr->getNumeredColDofsPtr();
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "Should be row or column");
  }
  marker.resize(dofs->size(), 0);
  std::vector<unsigned char> marker_tmp;

  switch (op) {
  case MarkOP::OR:
    for (auto p = ents.pair_begin(); p != ents.pair_end(); ++p) {
      auto lo = dofs->get<Ent_mi_tag>().lower_bound(p->first);
      auto hi = dofs->get<Ent_mi_tag>().upper_bound(p->second);
      for (; lo != hi; ++lo)
        marker[(*lo)->getPetscLocalDofIdx()] |= 1;
    }
    break;
  case MarkOP::AND:
    marker_tmp.resize(dofs->size(), 0);
    for (auto p = ents.pair_begin(); p != ents.pair_end(); ++p) {
      auto lo = dofs->get<Ent_mi_tag>().lower_bound(p->first);
      auto hi = dofs->get<Ent_mi_tag>().upper_bound(p->second);
      for (; lo != hi; ++lo)
        marker_tmp[(*lo)->getPetscLocalDofIdx()] = 1;
    }
    for (int i = 0; i != marker.size(); ++i) {
      marker[i] &= marker_tmp[i];
    }
    break;
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::modifyMarkDofs(
    const std::string problem_name, RowColData rc, const std::string field_name,
    const int lo, const int hi, const enum ProblemsManager::MarkOP op,
    const unsigned char c, std::vector<unsigned char> &marker) const {

  Interface &m_field = cOre;
  const Problem *problem_ptr;
  ProblemManagerFunctionBegin;
  CHKERR m_field.get_problem(problem_name, &problem_ptr);
  boost::shared_ptr<NumeredDofEntity_multiIndex> dofs;
  switch (rc) {
  case ROW:
    dofs = problem_ptr->getNumeredRowDofsPtr();
    break;
  case COL:
    dofs = problem_ptr->getNumeredColDofsPtr();
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "Should be row or column");
  }
  marker.resize(dofs->size(), 0);

  auto dof_lo = dofs->get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId(m_field.get_field_bit_number(field_name)));
  auto dof_hi = dofs->get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId(m_field.get_field_bit_number(field_name)));

  auto marker_ref = [marker](auto &it) -> unsigned int & {
    return marker[(*it)->getPetscLocalDofIdx()];
  };

  switch (op) {
  case MarkOP::OR:
    for (; dof_lo != dof_hi; ++dof_lo)
      if ((*dof_lo)->getDofCoeffIdx() >= lo &&
          (*dof_lo)->getDofCoeffIdx() <= hi)
        marker[(*dof_lo)->getPetscLocalDofIdx()] |= c;
    break;
  case MarkOP::AND:
    for (; dof_lo != dof_hi; ++dof_lo)
      if ((*dof_lo)->getDofCoeffIdx() >= lo &&
          (*dof_lo)->getDofCoeffIdx() <= hi)
        marker[(*dof_lo)->getPetscLocalDofIdx()] &= c;
    break;
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ProblemsManager::addFieldToEmptyFieldBlocks(const std::string problem_name,
                                            const std::string row_field,
                                            const std::string col_field) const {

  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  const auto problem_ptr = m_field.get_problem(problem_name);
  auto get_field_id = [&](const std::string field_name) {
    return m_field.get_field_structure(field_name)->getId();
  };
  const auto row_id = get_field_id(row_field);
  const auto col_id = get_field_id(col_field);

  problem_ptr->addFieldToEmptyFieldBlocks(EmptyFieldBlocks(row_id, col_id));

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
