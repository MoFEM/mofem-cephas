/** \file ProblemsManager.cpp
 * \brief Managing complexities for problem
 * \ingroup mofem_problems_manager
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

#ifdef PARMETIS

MoFEMErrorCode MatPartitioningApply_Parmetis_MoFEM(MatPartitioning part,
                                                   IS *partitioning);

#endif // PARMETIS

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
ProblemsManager::query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMProblemsManager) {
    *iface = const_cast<ProblemsManager *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
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
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  MOFEM_LOG_CHANNEL("WORLD");
  MOFEM_LOG_TAG("WORLD", "Partition mesh");
  MOFEM_LOG_CHANNEL("SYNC");
  MOFEM_LOG_TAG("SYNC", "Partition mesh");

  // get layout
  int rstart, rend, nb_elems;
  {
    PetscLayout layout;
    CHKERR PetscLayoutCreate(m_field.get_comm(), &layout);
    CHKERR PetscLayoutSetBlockSize(layout, 1);
    CHKERR PetscLayoutSetSize(layout, ents.size());
    CHKERR PetscLayoutSetUp(layout);
    CHKERR PetscLayoutGetSize(layout, &nb_elems);
    CHKERR PetscLayoutGetRange(layout, &rstart, &rend);
    CHKERR PetscLayoutDestroy(&layout);
    if (verb >= VERBOSE) {
      MOFEM_LOG("SYNC", LogManager::SeverityLevel::inform)
          << "Finite elements in problem: row lower " << rstart << " row upper "
          << rend << " nb. elems " << nb_elems << " ( " << ents.size() << " )";
      MOFEM_LOG_SYNCHORMISE(PETSC_COMM_WORLD)
    }
  }

  std::vector<EntityHandle> weight_ents;
  weight_ents.reserve(rend - rstart + 1);

  struct AdjBridge {
    EntityHandle ent;
    std::vector<int> adj;
    AdjBridge(const EntityHandle ent, std::vector<int> &adj):
      ent(ent), adj(adj) {}
  };

  typedef multi_index_container<
      AdjBridge,
      indexed_by<

          hashed_unique<member<AdjBridge, EntityHandle, &AdjBridge::ent>>

          >>
      AdjBridgeMap;

  Range all_dim_ents;
  CHKERR m_field.get_moab().get_adjacencies(ents, adj_dim, true, all_dim_ents,
                                            moab::Interface::UNION);

  AdjBridgeMap adj_bridge_map;
  auto hint = adj_bridge_map.begin();
  std::vector<int> adj;
  for (auto ent : all_dim_ents) {
    Range adj_ents;
    CHKERR m_field.get_moab().get_adjacencies(&ent, 1, dim, false, adj_ents);
    adj_ents = intersect(adj_ents, ents);
    adj.clear();
    adj.reserve(adj_ents.size());
    for (auto a : adj_ents)
      adj.emplace_back(ents.index(a));
    hint = adj_bridge_map.emplace_hint(hint, ent, adj);
  }

  int *_i;
  int *_j;
  {
    const int nb_loc_elements = rend - rstart;
    std::vector<int> i(nb_loc_elements + 1, 0), j;
    {
      std::vector<int> row_adj;
      Range::iterator fe_it;
      int ii, jj;
      size_t max_row_size;
      for (

          fe_it = ents.begin(), ii = 0, jj = 0, max_row_size = 0;

          fe_it != ents.end(); ++fe_it, ++ii) {

        if (ii < rstart)
          continue;
        if (ii >= rend)
          break;

        if (m_field.get_moab().type_from_handle(*fe_it) == MBENTITYSET) {
          SETERRQ(
              PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
              "not yet implemented, don't know what to do for meshset element");
        } else {

          Range dim_ents;
          CHKERR m_field.get_moab().get_adjacencies(&*fe_it, 1, adj_dim, false,
                                                    dim_ents);
          dim_ents = intersect(dim_ents, all_dim_ents);

          row_adj.clear();
          for (auto e : dim_ents) {
            auto adj_it = adj_bridge_map.find(e);
            if (adj_it != adj_bridge_map.end()) {
              
              for (const auto idx : adj_it->adj)
                row_adj.push_back(idx);

            } else
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                      "Entity not found");
          }

          std::sort(row_adj.begin(), row_adj.end());
          auto end = std::unique(row_adj.begin(), row_adj.end());

          size_t row_size = std::distance(row_adj.begin(), end);
          max_row_size = std::max(max_row_size, row_size);
          if (j.capacity() < (nb_loc_elements - jj) * max_row_size)
            j.reserve(nb_loc_elements * max_row_size);

          i[jj] = j.size();
          auto diag = ents.index(*fe_it);
          for (auto it = row_adj.begin(); it != end; ++it)
            if (*it != diag) 
              j.push_back(*it);
        }

        ++jj;

        if (th_vertex_weights != NULL)
          weight_ents.push_back(*fe_it);

      }

      i[jj] = j.size();
    }

    CHKERR PetscMalloc(i.size() * sizeof(int), &_i);
    CHKERR PetscMalloc(j.size() * sizeof(int), &_j);
    copy(i.begin(), i.end(), _i);
    copy(j.begin(), j.end(), _j);

  }

  // get weights
  int *vertex_weights = NULL;
  if (th_vertex_weights != NULL) {
    CHKERR PetscMalloc(weight_ents.size() * sizeof(int), &vertex_weights);
    CHKERR m_field.get_moab().tag_get_data(*th_vertex_weights,
                                           &*weight_ents.begin(),
                                           weight_ents.size(), vertex_weights);
  }

  {
    Mat Adj;
    // Adjacency matrix used to partition problems, f.e. METIS
    CHKERR MatCreateMPIAdj(m_field.get_comm(), rend - rstart, nb_elems, _i, _j,
                           PETSC_NULL, &Adj);
    CHKERR MatSetOption(Adj, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE);

    if (debug) {
      Mat A;
      CHKERR MatConvert(Adj, MATMPIAIJ, MAT_INITIAL_MATRIX, &A);
      CHKERR MatView(A, PETSC_VIEWER_DRAW_WORLD);
      std::string wait;
      std::cin >> wait;
      CHKERR MatDestroy(&A);
    }

    // run pets to do partitioning
    MOFEM_LOG("WORLD", LogManager::SeverityLevel::verbose) << "Start";

    MatPartitioning part;
    IS is;
    CHKERR MatPartitioningCreate(m_field.get_comm(), &part);

    CHKERR MatPartitioningSetAdjacency(part, Adj);
    CHKERR MatPartitioningSetFromOptions(part);
    CHKERR MatPartitioningSetNParts(part, n_parts);
    if (th_vertex_weights != NULL) {
      CHKERR MatPartitioningSetVertexWeights(part, vertex_weights);
    }
    PetscBool same;
    PetscObjectTypeCompare((PetscObject)part, MATPARTITIONINGPARMETIS, &same);
    if (same) {
#ifdef PARMETIS
      CHKERR MatPartitioningApply_Parmetis_MoFEM(part, &is);
#endif
    } else {
      CHKERR MatPartitioningApply(part, &is);
    }

    MOFEM_LOG("WORLD", LogManager::SeverityLevel::verbose) << "End";

    // gather
    IS is_gather, is_num, is_gather_num;
    CHKERR ISAllGather(is, &is_gather);
    CHKERR ISPartitioningToNumbering(is, &is_num);
    CHKERR ISAllGather(is_num, &is_gather_num);

    const int *part_number, *gids;
    CHKERR ISGetIndices(is_gather, &part_number);
    CHKERR ISGetIndices(is_gather_num, &gids);

    // set partition tag and gid tag to entities
    ParallelComm *pcomm = ParallelComm::get_pcomm(
        &m_field.get_moab(), m_field.get_basic_entity_data_ptr()->pcommID);
    Tag gid_tag;
    Tag part_tag = pcomm->part_tag();
    {
      const int zero = 0;
      CHKERR m_field.get_moab().tag_get_handle(
          GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid_tag,
          MB_TAG_DENSE | MB_TAG_CREAT, &zero);
      // get any sets already with this tag, and clear them
      CHKERR m_field.get_moab().tag_set_data(part_tag, ents, part_number);
      // rval = moab.tag_set_data(gid_tag,ents,&gids[0]); CHKERRQ_MOAB(rval);
      // std::vector<int> add_one(ents.size());
      // for(int ii = 0;ii<ents.size();ii++) {
      //   add_one[ii] = gids[ii]+1;
      // }
      // rval = moab.tag_set_data(gid_tag,ents,&add_one[0]); CHKERRQ_MOAB(rval);
    }

    std::map<int, Range> parts_ents;
    {
      // get entities on each part
      Range::iterator eit = ents.begin();
      for (int ii = 0; eit != ents.end(); eit++, ii++) {
        parts_ents[part_number[ii]].insert(*eit);
      }
      Range tagged_sets;
      CHKERR m_field.get_moab().get_entities_by_type_and_tag(
          0, MBENTITYSET, &part_tag, NULL, 1, tagged_sets,
          moab::Interface::UNION);
      // if(!tagged_sets.empty()) {
      //   CHKERR m_field.get_moab().delete_entities(tagged_sets);
      //   tagged_sets.clear();
      // }
      if (!tagged_sets.empty()) {
        CHKERR m_field.get_moab().tag_delete_data(part_tag, tagged_sets);
      }
      if (n_parts > (int)tagged_sets.size()) {
        // too few partition sets - create missing ones
        int num_new = n_parts - tagged_sets.size();
        for (int i = 0; i < num_new; i++) {
          EntityHandle new_set;
          CHKERR m_field.get_moab().create_meshset(
              MESHSET_SET | MESHSET_TRACK_OWNER, new_set);
          tagged_sets.insert(new_set);
        }
      } else if (n_parts < (int)tagged_sets.size()) {
        // too many partition sets - delete extras
        int num_del = tagged_sets.size() - n_parts;
        for (int i = 0; i < num_del; i++) {
          EntityHandle old_set = tagged_sets.pop_back();
          CHKERR m_field.get_moab().delete_entities(&old_set, 1);
        }
      }
      // write a tag to those sets denoting they're partition sets, with a value
      // of the proc number
      std::vector<int> dum_ids(n_parts);
      for (int i = 0; i < n_parts; i++)
        dum_ids[i] = i;
      CHKERR m_field.get_moab().tag_set_data(part_tag, tagged_sets,
                                             &*dum_ids.begin());
      CHKERR m_field.get_moab().clear_meshset(tagged_sets);

      // get lower dimension entities on each part
      for (int pp = 0; pp != n_parts; pp++) {
        Range dim_ents = parts_ents[pp].subset_by_dimension(dim);
        for (int dd = dim - 1; dd != -1; dd--) {
          Range adj_ents;
          if (dim > 0) {
            CHKERR m_field.get_moab().get_adjacencies(
                dim_ents, dd, true, adj_ents, moab::Interface::UNION);
          } else {
            CHKERR m_field.get_moab().get_connectivity(dim_ents, adj_ents,
                                                       true);
          }
          parts_ents[pp].merge(adj_ents);
          // std::cerr << pp << " add " << parts_ents[pp].size() << std::endl;
        }
      }
      for (int pp = 1; pp != n_parts; pp++) {
        for (int ppp = 0; ppp != pp; ppp++) {
          // std::cerr << pp << "<-" << ppp << " " << parts_ents[pp].size() << "
          // " << parts_ents[ppp].size();
          parts_ents[pp] = subtract(parts_ents[pp], parts_ents[ppp]);
          // std::cerr << " " << parts_ents[pp].size() << std::endl;
        }
      }
      if (debug) {
        for (int rr = 0; rr != m_field.get_comm_size(); rr++) {
          ostringstream ss;
          ss << "out_part_" << rr << ".vtk";
          EntityHandle meshset;
          CHKERR m_field.get_moab().create_meshset(MESHSET_SET, meshset);
          CHKERR m_field.get_moab().add_entities(meshset, parts_ents[rr]);
          CHKERR m_field.get_moab().write_file(ss.str().c_str(), "VTK", "",
                                               &meshset, 1);
          CHKERR m_field.get_moab().delete_entities(&meshset, 1);
        }
      }
      for (int pp = 0; pp != n_parts; pp++) {
        CHKERR m_field.get_moab().add_entities(tagged_sets[pp], parts_ents[pp]);
      }
      // for(int rr = 0;rr!=m_field.get_comm_size();rr++) {
      //   ostringstream ss;
      //   ss << "out_part_meshsets_" << rr << ".vtk";
      //   EntityHandle meshset = tagged_sets[rr];
      //   rval =
      //   m_field.get_moab().write_file(ss.str().c_str(),"VTK","",&meshset,1);
      //   CHKERRQ_MOAB(rval);
      // }

      // set gid to lower dimension entities
      for (int dd = 0; dd <= dim; dd++) {
        int gid = 1; // moab indexing from 1
        for (int pp = 0; pp != n_parts; pp++) {
          Range dim_ents = parts_ents[pp].subset_by_dimension(dd);
          // std::cerr << dim_ents.size() << " " << dd  << " " << pp <<
          // std::endl;
          for (Range::iterator eit = dim_ents.begin(); eit != dim_ents.end();
               eit++) {
            if (dd > 0) {
              CHKERR m_field.get_moab().tag_set_data(part_tag, &*eit, 1, &pp);
            }
            CHKERR m_field.get_moab().tag_set_data(gid_tag, &*eit, 1, &gid);
            gid++;
          }
        }
      }
    }

    CHKERR ISRestoreIndices(is_gather, &part_number);
    CHKERR ISRestoreIndices(is_gather_num, &gids);
    CHKERR ISDestroy(&is_num);
    CHKERR ISDestroy(&is_gather_num);
    CHKERR ISDestroy(&is_gather);
    CHKERR ISDestroy(&is);
    CHKERR MatPartitioningDestroy(&part);
    CHKERR MatDestroy(&Adj);
  }

  cOre.getBuildMoFEM() |= Core::PARTITION_MESH;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::buildProblem(const std::string name,
                                             const bool square_matrix,
                                             int verb) {

  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
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
  const EntFiniteElement_multiIndex *fe_ent_ptr;
  const DofEntity_multiIndex *dofs_field_ptr;
  MoFEMFunctionBegin;
  MOFEM_LOG_CHANNEL("SYNC");
  MOFEM_LOG_TAG("SYNC", "buildProblem");
  PetscLogEventBegin(MOFEM_EVENT_ProblemsManager, 0, 0, 0, 0);

  // Note: Only allowed changes on problem_ptr structure which not influence
  // multi-index.

  if (problem_ptr->getBitRefLevel().none()) {
    SETERRQ1(PETSC_COMM_SELF, 1, "problem <%s> refinement level not set",
             problem_ptr->getName().c_str());
  }
  CHKERR m_field.clear_problem(problem_ptr->getName());
  CHKERR m_field.get_ents_finite_elements(&fe_ent_ptr);
  CHKERR m_field.get_dofs(&dofs_field_ptr);

  // zero finite elements
  problem_ptr->numeredFiniteElements->clear();

  DofEntity_multiIndex_active_view dofs_rows, dofs_cols;
  {
    EntFiniteElement_multiIndex::iterator miit = fe_ent_ptr->begin();
    EntFiniteElement_multiIndex::iterator hi_miit = fe_ent_ptr->end();
    // iterate all finite element entities in database
    for (; miit != hi_miit; miit++) {
      // if element is in problem
      if (((*miit)->getId() & problem_ptr->getBitFEId()).any()) {
        BitRefLevel prb_bit = problem_ptr->getBitRefLevel();
        BitRefLevel prb_mask = problem_ptr->getMaskBitRefLevel();
        BitRefLevel fe_bit = (*miit)->getBitRefLevel();
        // if entity is not problem refinement level
        if ((fe_bit & prb_mask) != fe_bit)
          continue;
        if ((fe_bit & prb_bit) != prb_bit)
          continue;
        // get dof uids for rows and columns
        CHKERR(*miit)->getRowDofView(*dofs_field_ptr, dofs_rows);
        if (!square_matrix) {
          CHKERR(*miit)->getColDofView(*dofs_field_ptr, dofs_cols);
        }
      }
    }
  }

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
    auto hint = problem_ptr->numeredDofsRows->end();
    for (auto &v : *dofs_array) {
      hint = problem_ptr->numeredDofsRows->emplace_hint(hint, dofs_array, &v);
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
    auto hint = problem_ptr->numeredDofsCols->end();
    for (auto &v : *dofs_array) {
      hint = problem_ptr->numeredDofsCols->emplace_hint(hint, dofs_array, &v);
    }
  } else {
    problem_ptr->numeredDofsCols = problem_ptr->numeredDofsRows;
    problem_ptr->nbLocDofsCol = problem_ptr->nbLocDofsRow;
    problem_ptr->nbDofsCol = problem_ptr->nbDofsRow;
  }

  // job done, some debugging and postprocessing
  if (verb >= VERBOSE) {
    MOFEM_LOG("SYNC", LogManager::SeverityLevel::verbose)
        << problem_ptr->getName() << " Nb. local dofs "
        << problem_ptr->numeredDofsRows->size() << " by "
        << problem_ptr->numeredDofsCols->size();
    MOFEM_LOG_SYNCHORMISE(PETSC_COMM_WORLD);
  }

  if (verb >= NOISY) {
    MOFEM_LOG("SYNC", LogManager::SeverityLevel::noisy)
        << "FEs data for problem " << *problem_ptr;
    for(auto &miit : *fe_ent_ptr)
      MOFEM_LOG("SYNC", LogManager::SeverityLevel::noisy) << *miit;

    MOFEM_LOG("SYNC", LogManager::SeverityLevel::noisy)
        << "FEs row dofs " << *problem_ptr << " Nb. row dof "
        << problem_ptr->getNbDofsRow();
    for (auto &miit : *problem_ptr->numeredDofsRows)
      MOFEM_LOG("SYNC", LogManager::SeverityLevel::noisy) << *miit;

    MOFEM_LOG("SYNC", LogManager::SeverityLevel::noisy)
        << "FEs col dofs " << *problem_ptr << " Nb. col dof "
        << problem_ptr->getNbDofsCol();
    for (auto &miit : *problem_ptr->numeredDofsCols)
      MOFEM_LOG("SYNC", LogManager::SeverityLevel::noisy) << *miit;
    MOFEM_LOG_SYNCHORMISE(PETSC_COMM_WORLD);
  }

  cOre.getBuildMoFEM() |= Core::BUILD_PROBLEM; // It is assumed that user who
                                               // uses this function knows what
                                               // he is doing

  PetscLogEventEnd(MOFEM_EVENT_ProblemsManager, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::buildProblemOnDistributedMesh(
    const std::string name, const bool square_matrix, int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;

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
  const Field_multiIndex *fields_ptr;
  const FiniteElement_multiIndex *fe_ptr;
  const EntFiniteElement_multiIndex *fe_ent_ptr;
  const DofEntity_multiIndex *dofs_field_ptr;
  MoFEMFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_ProblemsManager, 0, 0, 0, 0);

  // clear data structures
  CHKERR m_field.clear_problem(problem_ptr->getName());

  CHKERR getOptions();
  CHKERR m_field.get_fields(&fields_ptr);
  CHKERR m_field.get_finite_elements(&fe_ptr);
  CHKERR m_field.get_ents_finite_elements(&fe_ent_ptr);
  CHKERR m_field.get_dofs(&dofs_field_ptr);

  if (problem_ptr->getBitRefLevel().none()) {
    SETERRQ1(PETSC_COMM_SELF, 1, "problem <%s> refinement level not set",
             problem_ptr->getName().c_str());
  }

  int loop_size = 2;
  if (square_matrix) {
    loop_size = 1;
    problem_ptr->numeredDofsCols = problem_ptr->numeredDofsRows;
  } else if (problem_ptr->numeredDofsCols == problem_ptr->numeredDofsRows) {
    problem_ptr->numeredDofsCols =
        boost::shared_ptr<NumeredDofEntity_multiIndex>(
            new NumeredDofEntity_multiIndex());
  }

  const BitRefLevel prb_bit = problem_ptr->getBitRefLevel();
  const BitRefLevel prb_mask = problem_ptr->getMaskBitRefLevel();

  // // get rows and cols dofs view based on data on elements
  DofEntity_multiIndex_active_view dofs_rows, dofs_cols;

  // Add DOFs to problem by visiting all elements and adding DOFs from
  // elements to the problem
  if (buildProblemFromFields == PETSC_FALSE) {
    // fe_miit iterator for finite elements
    EntFiniteElement_multiIndex::iterator fe_miit = fe_ent_ptr->begin();
    EntFiniteElement_multiIndex::iterator hi_fe_miit = fe_ent_ptr->end();
    // iterate all finite elements entities in database
    for (; fe_miit != hi_fe_miit; fe_miit++) {
      // if element is in problem
      if (((*fe_miit)->getId() & problem_ptr->getBitFEId()).any()) {

        const BitRefLevel fe_bit = (*fe_miit)->getBitRefLevel();
        // if entity is not problem refinement level
        if ((fe_bit & prb_mask) != fe_bit)
          continue;
        if ((fe_bit & prb_bit) != prb_bit)
          continue;

        // get dof uids for rows and columns
        CHKERR(*fe_miit)->getRowDofView(*dofs_field_ptr, dofs_rows);
        if (!square_matrix) {
          CHKERR(*fe_miit)->getColDofView(*dofs_field_ptr, dofs_cols);
        }
      }
    }
  }

  // Add DOFS to the proble by searching all the filedes, and adding to problem
  // owned or shared DOFs
  if (buildProblemFromFields == PETSC_TRUE) {
    // Get fields IDs on elements
    BitFieldId fields_ids_row, fields_ids_col;
    for (FiniteElement_multiIndex::iterator fit = fe_ptr->begin();
         fit != fe_ptr->end(); fit++) {
      if ((fit->get()->getId() & problem_ptr->getBitFEId()).any()) {
        fields_ids_row |= fit->get()->getBitFieldIdRow();
        fields_ids_col |= fit->get()->getBitFieldIdCol();
      }
    }
    // Get fields DOFs
    for (Field_multiIndex::iterator fit = fields_ptr->begin();
         fit != fields_ptr->end(); fit++) {
      if ((fit->get()->getId() & (fields_ids_row | fields_ids_col)).any()) {
        for (DofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit =
                 dofs_field_ptr->get<FieldName_mi_tag>().lower_bound(
                     fit->get()->getName());
             dit != dofs_field_ptr->get<FieldName_mi_tag>().upper_bound(
                        fit->get()->getName());
             dit++) {
          const int owner_proc = dit->get()->getOwnerProc();
          if (owner_proc != m_field.get_comm_rank()) {
            const unsigned char pstatus = dit->get()->getPStatus();
            if (pstatus == 0) {
              continue;
            }
          }
          const BitRefLevel dof_bit = (*dit)->getBitRefLevel();
          // if entity is not problem refinement level
          if ((dof_bit & prb_mask) != dof_bit)
            continue;
          if ((dof_bit & prb_bit) != prb_bit)
            continue;
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
          ents_to_synchronise, verb);
      ents_to_synchronise = subtract(ents_to_synchronise, tmp_ents);
      for (Field_multiIndex::iterator fit = fields_ptr->begin();
           fit != fields_ptr->end(); fit++) {
        if ((fit->get()->getId() & *fields_ids[ss]).any()) {
          for (Range::pair_iterator pit = ents_to_synchronise.pair_begin();
               pit != ents_to_synchronise.pair_end(); ++pit) {
            const EntityHandle f = pit->first;
            const EntityHandle s = pit->second;
            DofEntity_multiIndex::index<
                Composite_Name_And_Ent_mi_tag>::type::iterator dit,
                hi_dit;
            dit = dofs_field_ptr->get<Composite_Name_And_Ent_mi_tag>()
                      .lower_bound(boost::make_tuple(fit->get()->getName(), f));
            hi_dit =
                dofs_field_ptr->get<Composite_Name_And_Ent_mi_tag>()
                    .upper_bound(boost::make_tuple(fit->get()->getName(), s));
            dofs_ptr[ss]->insert(dit, hi_dit);
          }
        }
      }
    }
  }

  // add dofs for rows and cols and set ownership
  DofEntity_multiIndex_active_view *dofs_ptr[] = {&dofs_rows, &dofs_cols};
  boost::shared_ptr<NumeredDofEntity_multiIndex> numered_dofs_ptr[] = {
      problem_ptr->numeredDofsRows, problem_ptr->numeredDofsCols};
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

  // Loop over dofs on rows and columns and add to multi-indices in dofs problem
  // structure,  set partition for each dof
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

  // Loop over dofs on this processor and prepare those dofs to send on another
  // proc
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

  {

    // make sure it is a PETSc comm
    CHKERR PetscCommDuplicate(m_field.get_comm(), &m_field.get_comm(), NULL);

    // rows

    // Computes the number of messages a node expects to receive
    CHKERR PetscGatherNumberOfMessages(m_field.get_comm(), NULL,
                                       &lengths_rows[0], &nrecvs_rows);
    // std::cerr << nrecvs_rows << std::endl;

    // Computes info about messages that a MPI-node will receive, including
    // (from-id,length) pairs for each message.
    CHKERR PetscGatherMessageLengths(m_field.get_comm(), nsends_rows,
                                     nrecvs_rows, &lengths_rows[0],
                                     &onodes_rows, &olengths_rows);

    // Gets a unique new tag from a PETSc communicator. All processors that
    // share the communicator MUST call this routine EXACTLY the same number of
    // times. This tag should only be used with the current objects
    // communicator; do NOT use it with any other MPI communicator.
    int tag_row;
    CHKERR PetscCommGetNewTag(m_field.get_comm(), &tag_row);

    // Allocate a buffer sufficient to hold messages of size specified in
    // olengths. And post Irecvs on these buffers using node info from onodes
    MPI_Request *r_waits_row; // must bee freed by user
    // rbuf has a pointers to messeges. It has size of of nrecvs (number of
    // messages) +1. In the first index a block is allocated,
    // such that rbuf[i] = rbuf[i-1]+olengths[i-1].

    CHKERR PetscPostIrecvInt(m_field.get_comm(), tag_row, nrecvs_rows,
                             onodes_rows, olengths_rows, &rbuf_row,
                             &r_waits_row);
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
                       tag_row, m_field.get_comm(), s_waits_row + kk);
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
    CHKERR PetscGatherNumberOfMessages(m_field.get_comm(), NULL,
                                       &lengths_cols[0], &nrecvs_cols);

    // Computes info about messages that a MPI-node will receive, including
    // (from-id,length) pairs for each message.
    int *onodes_cols;
    CHKERR PetscGatherMessageLengths(m_field.get_comm(), nsends_cols,
                                     nrecvs_cols, &lengths_cols[0],
                                     &onodes_cols, &olengths_cols);

    // Gets a unique new tag from a PETSc communicator.
    int tag_col;
    CHKERR PetscCommGetNewTag(m_field.get_comm(), &tag_col);

    // Allocate a buffer sufficient to hold messages of size specified in
    // olengths. And post Irecvs on these buffers using node info from onodes
    MPI_Request *r_waits_col; // must bee freed by user
    CHKERR PetscPostIrecvInt(m_field.get_comm(), tag_col, nrecvs_cols,
                             onodes_cols, olengths_cols, &rbuf_col,
                             &r_waits_col);
    CHKERR PetscFree(onodes_cols);

    MPI_Request *s_waits_col; // status of sens messages
    CHKERR PetscMalloc1(nsends_cols, &s_waits_col);

    // Send messeges
    for (int proc = 0, kk = 0; proc < m_field.get_comm_size(); proc++) {
      if (!lengths_cols[proc])
        continue; // no message to send to this proc
      CHKERR MPI_Isend(&(ids_data_packed_cols[proc])[0], // buffer to send
                       lengths_cols[proc],               // message length
                       MPIU_INT, proc,                   // to proc
                       tag_col, m_field.get_comm(), s_waits_col + kk);
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

  CHKERR PetscFree(status);

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

    const DofEntity_multiIndex *dofs_ptr;
    CHKERR m_field.get_dofs(&dofs_ptr);

    UId uid;
    NumeredDofEntity_multiIndex::iterator dit;
    for (int kk = 0; kk != nrecvs; ++kk) {
      int len = olengths[kk];
      int *data_from_proc = data_procs[kk];
      for (int dd = 0; dd < len; dd += data_block_size) {
        uid = IdxDataTypePtr(&data_from_proc[dd]).getUId();
        dit = numered_dofs_ptr[ss]->find(uid);
        if (dit == numered_dofs_ptr[ss]->end()) {
          // Dof is shared to this processor, however there is no element which
          // have this dof
          // continue;
          DofEntity_multiIndex::iterator ddit = dofs_ptr->find(uid);
          if (ddit != dofs_ptr->end()) {
            unsigned char pstatus = ddit->get()->getPStatus();
            if (pstatus > 0) {
              continue;
            } else {
              std::ostringstream zz;
              zz << **ddit << std::endl;
              SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                       "data inconsistency, dofs is not shared, but received "
                       "from other proc\n"
                       "%s",
                       zz.str().c_str());
            }
          } else {
            std::ostringstream zz;
            zz << uid << std::endl;
            SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                     "no such dof %s in mofem database", zz.str().c_str());
          }
          std::ostringstream zz;
          zz << uid << std::endl;
          SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                   "dof %s not found", zz.str().c_str());
        }
        int global_idx = IdxDataTypePtr(&data_from_proc[dd]).getDofIdx();
        if (global_idx < 0) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "received negative dof");
        }
        bool success;
        success = numered_dofs_ptr[ss]->modify(
            dit, NumeredDofEntity_mofem_index_change(global_idx));
        if (!success)
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
        success = numered_dofs_ptr[ss]->modify(
            dit, NumeredDofEntity_part_and_glob_idx_change((*dit)->getPart(),
                                                           global_idx));
        if (!success)
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
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
    if (problem_ptr->numeredDofsRows != problem_ptr->numeredDofsCols) {
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
    const std::string out_name, const std::vector<std::string> &fields_row,
    const std::vector<std::string> &fields_col, const std::string main_problem,
    const bool square_matrix,
    const map<std::string, std::pair<EntityType, EntityType>> *entityMapRow,
    const map<std::string, std::pair<EntityType, EntityType>> *entityMapCol,
    int verb) {
  MoFEM::Interface &m_field = cOre;
  const Problem_multiIndex *problems_ptr;
  MoFEMFunctionBegin;

  CHKERR m_field.clear_problem(out_name);
  CHKERR m_field.get_problems(&problems_ptr);

  // get reference to all problems
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemByName;
  ProblemByName &problems_by_name =
      const_cast<ProblemByName &>(problems_ptr->get<Problem_mi_tag>());

  // get iterators to out problem, i.e. build problem
  ProblemByName::iterator out_problem_it = problems_by_name.find(out_name);
  if (out_problem_it == problems_by_name.end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "subproblem with name < %s > not defined (top tip check spelling)",
             out_name.c_str());
  }
  // get iterator to main problem, i.e. out problem is subproblem of main
  // problem
  ProblemByName::iterator main_problem_it = problems_by_name.find(main_problem);
  if (main_problem_it == problems_by_name.end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "problem of subproblem with name < %s > not defined (top tip "
             "check spelling)",
             main_problem.c_str());
  }

  // get dofs for row & columns for out problem,
  boost::shared_ptr<NumeredDofEntity_multiIndex> out_problem_dofs[] = {
      out_problem_it->numeredDofsRows, out_problem_it->numeredDofsCols};
  // get dofs for row & columns for main problem
  boost::shared_ptr<NumeredDofEntity_multiIndex> main_problem_dofs[] = {
      main_problem_it->numeredDofsRows, main_problem_it->numeredDofsCols};
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
  const map<std::string, std::pair<EntityType, EntityType>> *entityMap[] = {
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
    out_problem_it->numeredFiniteElements->clear();

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
      auto dit =
          main_problem_dofs[ss]->get<FieldName_mi_tag>().lower_bound(field);
      auto hi_dit =
          main_problem_dofs[ss]->get<FieldName_mi_tag>().upper_bound(field);

      auto add_dit_to_dofs_array = [&](auto &dit) {
        dofs_array->emplace_back(
            dit->get()->getDofEntityPtr(), dit->get()->getPetscGlobalDofIdx(),
            dit->get()->getPetscGlobalDofIdx(),
            dit->get()->getPetscLocalDofIdx(), dit->get()->getPart());
      };

      if (entityMap[ss]) {
        auto mit = entityMap[ss]->find(field);
        if (mit != entityMap[ss]->end()) {
          EntityType lo_type = mit->second.first;
          EntityType hi_type = mit->second.second;
          int count = 0;
          for (auto diit = dit; diit != hi_dit; ++diit) {
            EntityType ent_type = (*diit)->getEntType();
            if (ent_type >= lo_type && ent_type <= hi_type)
              ++count;
          }
          dofs_array->reserve(count);
          for (; dit != hi_dit; ++dit) {
            EntityType ent_type = (*dit)->getEntType();
            if (ent_type >= lo_type && ent_type <= hi_type)
              add_dit_to_dofs_array(dit);
          }
        } else {
          dofs_array->reserve(std::distance(dit, hi_dit));
          for (; dit != hi_dit; dit++)
            add_dit_to_dofs_array(dit);
        }
      } else {
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
        CHKERR ISDuplicate(is, &(out_problem_it->getSubData()->rowIs));
        // CHKERR ISSort(out_problem_it->getSubData()->rowIs);
        out_problem_it->getSubData()->rowMap = ao;
        CHKERR PetscObjectReference((PetscObject)ao);
      } else {
        CHKERR ISDuplicate(is, &(out_problem_it->getSubData()->colIs));
        // CHKERR ISSort(out_problem_it->getSubData()->colIs);
        out_problem_it->getSubData()->colMap = ao;
        CHKERR PetscObjectReference((PetscObject)ao);
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
    out_problem_it->numeredDofsCols = out_problem_it->numeredDofsRows;
    out_problem_it->nbLocDofsCol = out_problem_it->nbLocDofsRow;
    out_problem_it->nbDofsCol = out_problem_it->nbDofsRow;
    out_problem_it->getSubData()->colIs = out_problem_it->getSubData()->rowIs;
    out_problem_it->getSubData()->colMap = out_problem_it->getSubData()->rowMap;
    CHKERR PetscObjectReference(
        (PetscObject)out_problem_it->getSubData()->rowIs);
    CHKERR PetscObjectReference(
        (PetscObject)out_problem_it->getSubData()->rowMap);
  }

  CHKERR printPartitionedProblem(&*out_problem_it, verb);
  CHKERR debugPartitionedProblem(&*out_problem_it, verb);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::buildCompsedProblem(
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
  const Problem_multiIndex *problems_ptr;
  MoFEMFunctionBegin;

  CHKERR m_field.clear_problem(out_name);
  CHKERR m_field.get_problems(&problems_ptr);
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
  std::vector<IS> *add_prb_is[] = {&cmp_prb_data->rowIs, &cmp_prb_data->colIs};

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
    // cerr << "SS " << ss << endl;
    // cerr << add_prb[ss]->size() << endl;
    // cerr << add_prb_ptr[ss]->size() << endl;
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
        nb_dofs_reserve[ss] += add_prb_ptr[ss]->back()->numeredDofsRows->size();
      } else {
        // column
        *nb_dofs[ss] += add_prb_ptr[ss]->back()->getNbDofsCol();
        *nb_local_dofs[ss] += add_prb_ptr[ss]->back()->getNbLocalDofsCol();
        nb_dofs_reserve[ss] += add_prb_ptr[ss]->back()->numeredDofsCols->size();
      }
    }
  }
  // if squre problem, rows and columns are the same
  if (square_matrix) {
    add_prb_ptr[1]->reserve(add_prb_ptr[0]->size());
    add_prb_is[1]->reserve(add_prb_ptr[0]->size());
    out_problem_it->numeredDofsCols = out_problem_it->numeredDofsRows;
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
                  ->numeredDofsRows->get<PetscGlobalIdx_mi_tag>()
                  .begin();
        hi_dit = (*add_prb_ptr[ss])[pp]
                     ->numeredDofsRows->get<PetscGlobalIdx_mi_tag>()
                     .end();
      } else {
        dit = (*add_prb_ptr[ss])[pp]
                  ->numeredDofsCols->get<PetscGlobalIdx_mi_tag>()
                  .begin();
        hi_dit = (*add_prb_ptr[ss])[pp]
                     ->numeredDofsCols->get<PetscGlobalIdx_mi_tag>()
                     .end();
      }
      int is_nb = 0;
      for (; dit != hi_dit; dit++) {
        BitRefLevel prb_bit = out_problem_it->getBitRefLevel();
        BitRefLevel prb_mask = out_problem_it->getMaskBitRefLevel();
        BitRefLevel dof_bit = dit->get()->getBitRefLevel();
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
      (*add_prb_is[ss]).push_back(is);
      if (ss == 0) {
        shift_glob += (*add_prb_ptr[ss])[pp]->getNbDofsRow();
        shift_loc += (*add_prb_ptr[ss])[pp]->getNbLocalDofsRow();
      } else {
        shift_glob += (*add_prb_ptr[ss])[pp]->getNbDofsCol();
        shift_loc += (*add_prb_ptr[ss])[pp]->getNbLocalDofsCol();
      }
      if (square_matrix) {
        (*add_prb_ptr[1]).push_back((*add_prb_ptr[0])[pp]);
        (*add_prb_is[1]).push_back(is);
        CHKERR PetscObjectReference((PetscObject)is);
      }
    }
  }

  if ((*add_prb_is[1]).size() != (*add_prb_is[0]).size()) {
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY, "Data inconsistency");
  }

  // Insert DOFs to problem multi-index
  for (int ss = 0; ss != ((square_matrix) ? 1 : 2); ss++) {
    auto hint = (ss == 0) ? out_problem_it->numeredDofsRows->end()
                          : out_problem_it->numeredDofsCols->end();
    for (auto &v : *dofs_array[ss])
      hint = (ss == 0) ? out_problem_it->numeredDofsRows->emplace_hint(
                             hint, dofs_array[ss], &v)
                       : out_problem_it->numeredDofsCols->emplace_hint(
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
      dofs_ptr = out_problem_it->numeredDofsRows;
    } else {
      dofs_ptr = out_problem_it->numeredDofsCols;
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
  const Problem_multiIndex *problems_ptr;
  MoFEMFunctionBegin;
  MOFEM_LOG_CHANNEL("WORLD");

  if (!(cOre.getBuildMoFEM() & Core::BUILD_FIELD))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "fields not build");
  if (!(cOre.getBuildMoFEM() & Core::BUILD_FE))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "FEs not build");
  if (!(cOre.getBuildMoFEM() & Core::BUILD_ADJ))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "adjacencies not build");
  if (!(cOre.getBuildMoFEM() & Core::BUILD_PROBLEM))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "adjacencies not build");
  MOFEM_LOG("WORLD", LogManager::SeverityLevel::verbose)
      << "Simple partition problem " << name;

  CHKERR m_field.get_problems(&problems_ptr);
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
      p_miit->numeredDofsRows->get<Idx_mi_tag>();
  NumeredDofEntitysByIdx &dofs_col_by_idx =
      p_miit->numeredDofsCols->get<Idx_mi_tag>();
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
  if (p_miit->numeredDofsRows == p_miit->numeredDofsCols) {
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
        SETERRQ4(
            PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
            "data inconsistency, std::distance(miit_col,hi_miit_col) != rend - "
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
  const Problem_multiIndex *problems_ptr;
  MoFEMFunctionBegin;
  MOFEM_LOG_CHANNEL("WORLD");
  MOFEM_LOG_TAG("WORLD", "partitionProblem");
  MOFEM_LOG("WORLD", LogManager::SeverityLevel::noisy) << name;


  if (!(cOre.getBuildMoFEM() & (1 << 0)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "fields not build");
  if (!(cOre.getBuildMoFEM() & (1 << 1)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "FEs not build");
  if (!(cOre.getBuildMoFEM() & (1 << 2)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "entFEAdjacencies not build");
  if (!(cOre.getBuildMoFEM() & (1 << 3)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Problems not build");


  typedef NumeredDofEntity_multiIndex::index<Idx_mi_tag>::type
      NumeredDofEntitysByIdx;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;

  // Find problem pointer by name
  CHKERR m_field.get_problems(&problems_ptr);
  auto &problems_set =
      const_cast<ProblemsByName &>(problems_ptr->get<Problem_mi_tag>());
  auto p_miit = problems_set.find(name);
  if (p_miit == problems_set.end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
             "problem with name %s not defined (top tip check spelling)",
             name.c_str());
  }
  int nb_dofs_row = p_miit->getNbDofsRow();

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
    SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "data inconsistency %d != %d",
             size_is_gather, nb_dofs_row);

  CHKERR ISGetSize(is_num, &size_is_num);
  if (size_is_num != (int)nb_dofs_row)
    SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "data inconsistency %d != %d",
             size_is_num, nb_dofs_row);

  bool square_matrix = false;
  if (p_miit->numeredDofsRows == p_miit->numeredDofsCols)
    square_matrix = true;

  if (!square_matrix) {
    // FIXME: This is for back compatibility, if deprecate interface function
    // build interfaces is removed, this part of the code will be obsolete
    auto mit_row = p_miit->numeredDofsRows->get<Idx_mi_tag>().begin();
    auto hi_mit_row = p_miit->numeredDofsRows->get<Idx_mi_tag>().end();
    auto mit_col = p_miit->numeredDofsCols->get<Idx_mi_tag>().begin();
    auto hi_mit_col = p_miit->numeredDofsCols->get<Idx_mi_tag>().end();
    for (; mit_row != hi_mit_row; mit_row++, mit_col++) {
      if (mit_col == hi_mit_col) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "check finite element definition, nb. of rows is not equal to "
                "number for columns");
      }
      if (mit_row->get()->getGlobalUniqueId() !=
          mit_col->get()->getGlobalUniqueId()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "check finite element definition, nb. of rows is not equal to "
                "number for columns");
      }
    }
  }

  // Set petsc global indices
  auto &dofs_row_by_idx_no_const = const_cast<NumeredDofEntitysByIdx &>(
      p_miit->numeredDofsRows->get<Idx_mi_tag>());
  int &nb_row_local_dofs = p_miit->nbLocDofsRow;
  int &nb_row_ghost_dofs = p_miit->nbGhostDofsRow;
  nb_row_local_dofs = 0;
  nb_row_ghost_dofs = 0;

  for (auto miit_dofs_row = dofs_row_by_idx_no_const.begin();
       miit_dofs_row != dofs_row_by_idx_no_const.end(); miit_dofs_row++) {
    const int part = part_number[(*miit_dofs_row)->dofIdx];
    if (part == (unsigned int)m_field.get_comm_rank()) {
      const bool success = dofs_row_by_idx_no_const.modify(
          miit_dofs_row,
          NumeredDofEntity_part_and_indices_change(
              part, petsc_idx[(*miit_dofs_row)->dofIdx], nb_row_local_dofs++));
      if (!success) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                "modification unsuccessful");
      }
    } else {
      const bool success = dofs_row_by_idx_no_const.modify(
          miit_dofs_row, NumeredDofEntity_part_and_glob_idx_change(
                             part, petsc_idx[(*miit_dofs_row)->dofIdx]));
      if (!success) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                "modification unsuccessful");
      }
    }
  }

  int &nb_col_local_dofs = p_miit->nbLocDofsCol;
  int &nb_col_ghost_dofs = p_miit->nbGhostDofsCol;
  if (square_matrix) {
    nb_col_local_dofs = nb_row_local_dofs;
    nb_col_ghost_dofs = nb_row_ghost_dofs;
  } else {
    NumeredDofEntitysByIdx &dofs_col_by_idx_no_const =
        const_cast<NumeredDofEntitysByIdx &>(
            p_miit->numeredDofsCols->get<Idx_mi_tag>());
    nb_col_local_dofs = 0;
    nb_col_ghost_dofs = 0;
    for (auto miit_dofs_col = dofs_col_by_idx_no_const.begin();
         miit_dofs_col != dofs_col_by_idx_no_const.end(); miit_dofs_col++) {
      const int part = part_number[(*miit_dofs_col)->dofIdx];
      if (part == (unsigned int)m_field.get_comm_rank()) {
        const bool success = dofs_col_by_idx_no_const.modify(
            miit_dofs_col, NumeredDofEntity_part_and_indices_change(
                               part, petsc_idx[(*miit_dofs_col)->dofIdx],
                               nb_col_local_dofs++));
        if (!success) {
          SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
        }
      } else {
        const bool success = dofs_col_by_idx_no_const.modify(
            miit_dofs_col, NumeredDofEntity_part_and_glob_idx_change(
                               part, petsc_idx[(*miit_dofs_col)->dofIdx]));
        if (!success) {
          SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
        }
      }
    }
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

  cOre.getBuildMoFEM() |= 1 << 4;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::inheritPartition(
    const std::string name, const std::string problem_for_rows, bool copy_rows,
    const std::string problem_for_cols, bool copy_cols, int verb) {

  MoFEM::Interface &m_field = cOre;
  const Problem_multiIndex *problems_ptr;
  MoFEMFunctionBegin;

  if (!(cOre.getBuildMoFEM() & Core::BUILD_PROBLEM))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY, "pRoblems not build");

  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemByName;

  // find p_miit
  CHKERR m_field.get_problems(&problems_ptr);
  ProblemByName &problems_by_name =
      const_cast<ProblemByName &>(problems_ptr->get<Problem_mi_tag>());
  ProblemByName::iterator p_miit = problems_by_name.find(name);
  if (p_miit == problems_by_name.end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_NOT_FOUND,
             "problem with name < %s > not defined (top tip check spelling)",
             name.c_str());
  }
  if (verb > QUIET)
    PetscPrintf(m_field.get_comm(),
                "Compose problem %s from rows of %s and columns of %s\n",
                p_miit->getName().c_str(), problem_for_rows.c_str(),
                problem_for_cols.c_str());

  // find p_miit_row
  ProblemByName::iterator p_miit_row = problems_by_name.find(problem_for_rows);
  if (p_miit_row == problems_by_name.end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "problem with name < %s > not defined (top tip check spelling)",
             problem_for_rows.c_str());
  }
  const boost::shared_ptr<NumeredDofEntity_multiIndex> dofs_row =
      p_miit_row->numeredDofsRows;

  // find p_mit_col
  ProblemByName::iterator p_miit_col = problems_by_name.find(problem_for_cols);
  if (p_miit_col == problems_by_name.end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "problem with name < %s > not defined (top tip check spelling)",
             problem_for_cols.c_str());
  }
  boost::shared_ptr<NumeredDofEntity_multiIndex> dofs_col =
      p_miit_col->numeredDofsCols;

  bool copy[] = {copy_rows, copy_cols};
  boost::shared_ptr<NumeredDofEntity_multiIndex> composed_dofs[] = {
      p_miit->numeredDofsRows, p_miit->numeredDofsCols};

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
            dofs_by_uid.find((*dit)->getGlobalUniqueId());
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
  MoFEMFunctionBeginHot;
  MOFEM_LOG_CHANNEL("SYNC");
  MOFEM_LOG_TAG("SYNC", "partitionedProblem");

  if (verb > QUIET) {

    MOFEM_LOG("SYNC", LogManager::SeverityLevel::inform)
        << problem_ptr->getName() << " Nb. local dof "
        << problem_ptr->getNbLocalDofsRow() << " by "
        << problem_ptr->getNbLocalDofsCol() << " nb global dofs "
        << problem_ptr->getNbDofsRow() << " by " << problem_ptr->getNbDofsCol();

    MOFEM_LOG_SYNCHORMISE(PETSC_COMM_WORLD)
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
ProblemsManager::debugPartitionedProblem(const Problem *problem_ptr, int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBeginHot;
  if (debug > 0) {

    typedef NumeredDofEntity_multiIndex::index<Idx_mi_tag>::type
        NumeredDofEntitysByIdx;
    NumeredDofEntitysByIdx::iterator dit, hi_dit;
    const NumeredDofEntitysByIdx *numered_dofs_ptr[] = {
        &(problem_ptr->numeredDofsRows->get<Idx_mi_tag>()),
        &(problem_ptr->numeredDofsCols->get<Idx_mi_tag>())};

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
            SETERRQ2(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
                     "local dof index for %d (0-row, 1-col) not set, i.e. has "
                     "negative value\n %s",
                     ss, zz.str().c_str());
          }
          if ((*dit)->getPetscLocalDofIdx() >= *local_nbdof_ptr[ss]) {
            std::ostringstream zz;
            zz << "rank " << m_field.get_comm_rank() << " " << **dit;
            SETERRQ2(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
                     "local dofs for %d (0-row, 1-col) out of range\n %s", ss,
                     zz.str().c_str());
          }
        } else {
          if ((*dit)->getPetscGlobalDofIdx() < 0) {
            std::ostringstream zz;
            zz << "rank " << m_field.get_comm_rank() << " "
               << dit->get()->getBitRefLevel() << " " << **dit;
            SETERRQ2(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
                     "global dof index for %d (0-row, 1-col) row not set, i.e. "
                     "has negative value\n %s",
                     ss, zz.str().c_str());
          }
          if ((*dit)->getPetscGlobalDofIdx() >= *nbdof_ptr[ss]) {
            std::ostringstream zz;
            zz << "rank " << m_field.get_comm_rank() << " nb_dofs "
               << *nbdof_ptr[ss] << " " << **dit;
            SETERRQ2(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
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
  const Problem_multiIndex *problems_ptr;
  const EntFiniteElement_multiIndex *fe_ent_ptr;
  MoFEMFunctionBegin;
  MOFEM_LOG_CHANNEL("WORLD");
  MOFEM_LOG_CHANNEL("SYNC");

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

  // Get pointer to problem data struture
  CHKERR m_field.get_problems(&problems_ptr);

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
      *p_miit->numeredFiniteElements;

  // Clear all elements and data, build it again
  problem_finite_elements.clear();

  // Check if dofs and columns are the same, i.e. structurally symmetric problem
  bool do_cols_prob = true;
  if (p_miit->numeredDofsRows == p_miit->numeredDofsCols) {
    do_cols_prob = false;
  }

  // Loop over all elements in database and if right element is there add it
  // to problem finite element multi-index
  CHKERR m_field.get_ents_finite_elements(&fe_ent_ptr);
  for (auto efit = fe_ent_ptr->begin(); efit != fe_ent_ptr->end(); efit++) {

    // if element is not part of problem
    if (((*efit)->getId() & p_miit->getBitFEId()).none())
      continue;

    BitRefLevel prb_bit = p_miit->getBitRefLevel();
    BitRefLevel prb_mask = p_miit->getMaskBitRefLevel();
    BitRefLevel fe_bit = (*efit)->getBitRefLevel();
    // if entity is not problem refinement level
    if ((fe_bit & prb_mask) != fe_bit)
      continue;
    if ((fe_bit & prb_bit) != prb_bit)
      continue;

    // create element
    boost::shared_ptr<NumeredEntFiniteElement> numered_fe(
        new NumeredEntFiniteElement(*efit));

    // check if rows and columns are the same on this element
    bool do_cols_fe = true;
    if ((numered_fe->sPtr->row_field_ents_view ==
         numered_fe->sPtr->col_field_ents_view) &&
        !do_cols_prob) {
      do_cols_fe = false;
      numered_fe->cols_dofs = numered_fe->rows_dofs;
    } else {
      // different dofs on rows and columns
      numered_fe->cols_dofs = boost::shared_ptr<FENumeredDofEntity_multiIndex>(
          new FENumeredDofEntity_multiIndex());
    }
    // get pointer to dofs multi-index on rows and columns
    auto rows_dofs = numered_fe->rows_dofs;
    auto cols_dofs = numered_fe->cols_dofs;
    // clear multi-indices
    rows_dofs->clear();
    if (do_cols_fe) {
      cols_dofs->clear();
    }
    NumeredDofEntity_multiIndex_uid_view_ordered rows_view;
    NumeredDofEntity_multiIndex_uid_view_ordered cols_view;

    // set partition to the element
    {
      if (part_from_moab) {
        // if partition is taken from moab partition
        int proc = (*efit)->getPartProc();
        if (proc == -1 && (*efit)->getEntType() == MBVERTEX) {
          proc = (*efit)->getOwnerProc();
        }
        NumeredEntFiniteElement_change_part(proc).operator()(numered_fe);
        
      } else {

        // Count partition of the dofs in row, the larges dofs with given
        // partition is used to set partition of the element
        CHKERR(*efit)->getRowDofView(*(p_miit->numeredDofsRows), rows_view,
                                     moab::Interface::UNION);
        std::vector<int> parts(m_field.get_comm_size(), 0);
        for(auto &dof_ptr : rows_view)
          parts[dof_ptr->pArt]++;
        std::vector<int>::iterator pos =
            max_element(parts.begin(), parts.end());
        unsigned int max_part = std::distance(parts.begin(), pos);
        NumeredEntFiniteElement_change_part(max_part).operator()(numered_fe);
      }
    }

    // set dofs on rows and columns (if are different)
    if ((numered_fe->getPart() >= (unsigned int)low_proc) &&
        (numered_fe->getPart() <= (unsigned int)hi_proc)) {

      std::array<NumeredDofEntity_multiIndex_uid_view_ordered *, 2> dofs_view{
          &rows_view, &cols_view};
      std::array<FENumeredDofEntity_multiIndex *, 2> fe_dofs{rows_dofs.get(),
                                                             cols_dofs.get()};

      for (int ss = 0; ss != (do_cols_fe ? 2 : 1); ss++) {

        if (ss == 0) {
          if (part_from_moab) {
            // get row_view
            CHKERR(*efit)->getRowDofView(*(p_miit->numeredDofsRows),
                                         *dofs_view[ss],
                                         moab::Interface::UNION);
          }
        } else {
          // get cols_views
          CHKERR(*efit)->getColDofView(*(p_miit->numeredDofsCols),
                                       *dofs_view[ss], moab::Interface::UNION);
        }

        // Following reserve memory in sequences, only two allocations are here,
        // once for array of objects, next for array of shared pointers

        // reserve memory for field  dofs
        boost::shared_ptr<std::vector<FENumeredDofEntity>> dofs_array(
            new std::vector<FENumeredDofEntity>());

        if (!ss) {
          numered_fe->getRowDofsSequence() = dofs_array;
          if (!do_cols_fe)
            numered_fe->getColDofsSequence() = dofs_array;
        } else
          numered_fe->getColDofsSequence() = dofs_array;

        auto vit = dofs_view[ss]->begin();
        auto hi_vit = dofs_view[ss]->end();

        dofs_array->reserve(std::distance(vit, hi_vit));

        // create elements objects
        for (; vit != hi_vit; vit++) {
          boost::shared_ptr<SideNumber> side_number_ptr;
          side_number_ptr = (*efit)->getSideNumberPtr((*vit)->getEnt());
          dofs_array->emplace_back(side_number_ptr, *vit);
        }

        // finally add DoFS to multi-indices
        auto hint = fe_dofs[ss]->end();
        for (auto &v : *dofs_array)
          hint = fe_dofs[ss]->emplace_hint(hint, dofs_array, &v);
      }
    }

    auto check_fields_and_dofs = [part_from_moab,
                                  &m_field](const auto &numered_fe) {
      auto &fe = *(numered_fe->sPtr);


      if (!part_from_moab) {

        if (fe.getBitFieldIdRow().none() && m_field.get_comm_size() == 0) {
          MOFEM_LOG("WORLD", LogManager::SeverityLevel::warning)
              << "At least one field has to be added to element row to "
                 "determine partition of finite element. Check element " +
                     boost::lexical_cast<std::string>(fe.getName());
        }
      }

      // Adding elements if row or column has DOFs, or there is no field set to
      // rows and columns. The second case would be used by elements performing
      // tasks which do not assemble matrices or vectors, but evaluate fields or
      // modify base functions.

      return (!fe.row_field_ents_view->empty() ||
              !fe.col_field_ents_view->empty())

             ||

             (fe.getBitFieldIdRow().none() || fe.getBitFieldIdCol().none());
    };

    if (check_fields_and_dofs(numered_fe)) {
      // Add element to the problem
      auto p = problem_finite_elements.insert(numered_fe);
      if (!p.second)
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "element is there");

      if (verb >= VERY_VERBOSE) {
        MOFEM_LOG("SYNC", LogManager::SeverityLevel::very_verbose) << *p_miit;
        MOFEM_LOG("SYNC", LogManager::SeverityLevel::very_verbose) << *p.first;
        for (auto &dof_ptr : (*p.first)->rows_dofs->get<Unique_mi_tag>())
          MOFEM_LOG("SYNC", LogManager::SeverityLevel::noisy)
              << "row: " << *dof_ptr;
        for (auto &dof_ptr : (*p.first)->cols_dofs->get<Unique_mi_tag>())
          MOFEM_LOG("SYNC", LogManager::SeverityLevel::noisy)
              << "col: " << *dof_ptr;
      }
    }

  }

  if (verb >= VERBOSE) {
    auto elements_on_rank =
        problem_finite_elements.get<Part_mi_tag>().equal_range(
            m_field.get_comm_rank());
    MOFEM_LOG("SYNC", LogManager::SeverityLevel::verbose)
        << *p_miit << " Nb. elems "
        << std::distance(elements_on_rank.first, elements_on_rank.second);
    MOFEM_LOG_SYNCHORMISE(PETSC_COMM_WORLD);
  }

  cOre.getBuildMoFEM() |= Core::PARTITION_FE;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::partitionGhostDofs(const std::string name,
                                                   int verb) {
  MoFEM::Interface &m_field = cOre;
  const Problem_multiIndex *problems_ptr;
  MoFEMFunctionBegin;

  if (!(cOre.getBuildMoFEM() & Core::PARTITION_PROBLEM))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "partition of problem not build");
  if (!(cOre.getBuildMoFEM() & Core::PARTITION_FE))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "partitions finite elements not build");

  // find p_miit
  CHKERR m_field.get_problems(&problems_ptr);

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
        p_miit->numeredFiniteElements->get<Part_mi_tag>().equal_range(
            m_field.get_comm_rank());

    // get dofs on elements which are not part of this partition

    // rows
    auto hint_r = ghost_idx_row_view.begin();
    for (auto fe_ptr = fe_range.first; fe_ptr != fe_range.second; ++fe_ptr) {
      for (auto &dof_ptr : *(*fe_ptr)->rows_dofs) {
        if (dof_ptr->getPart() != (unsigned int)m_field.get_comm_rank()) {
          hint_r = ghost_idx_row_view.emplace_hint(
              hint_r, dof_ptr->getNumeredDofEntityPtr());
        }
      }
    }

    // columns
    if (p_miit->numeredDofsCols == p_miit->numeredDofsRows) {
      auto hint_c = ghost_idx_col_view.begin();
      for (auto fe_ptr = fe_range.first; fe_ptr != fe_range.second; ++fe_ptr) {
        for (auto &dof_ptr : *(*fe_range.first)->cols_dofs) {
          if (dof_ptr->getPart() != (unsigned int)m_field.get_comm_rank()) {
            hint_c = ghost_idx_col_view.emplace_hint(
                hint_c, dof_ptr->getNumeredDofEntityPtr());
          }
        }
      }
    }

    int *nb_ghost_dofs[2] = {&nb_row_ghost_dofs, &nb_col_ghost_dofs};
    int nb_local_dofs[2] = {p_miit->nbLocDofsRow, p_miit->nbLocDofsCol};

    NumeredDofEntity_multiIndex_uid_view_ordered *ghost_idx_view[2] = {
        &ghost_idx_row_view, &ghost_idx_col_view};
    NumeredDofEntityByUId *dof_by_uid_no_const[2] = {
        &p_miit->numeredDofsRows->get<Unique_mi_tag>(),
        &p_miit->numeredDofsCols->get<Unique_mi_tag>()};

    int loop_size = 2;
    if (p_miit->numeredDofsCols == p_miit->numeredDofsRows) {
      loop_size = 1;
    }

    // set local ghost dofs indices
    for (int ss = 0; ss != loop_size; ++ss) {
      for (auto &gid : *ghost_idx_view[ss]) {
        NumeredDofEntityByUId::iterator dof =
            dof_by_uid_no_const[ss]->find(gid->getGlobalUniqueId());
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
    std::ostringstream ss;
    ss << "partition_ghost_col_dofs: rank = " << m_field.get_comm_rank()
       << " FEs col ghost dofs " << *p_miit << " Nb. col ghost dof "
       << p_miit->getNbGhostDofsCol() << " Nb. local dof "
       << p_miit->getNbLocalDofsCol() << std::endl;
    ss << "partition_ghost_row_dofs: rank = " << m_field.get_comm_rank()
       << " FEs row ghost dofs " << *p_miit << " Nb. row ghost dof "
       << p_miit->getNbGhostDofsRow() << " Nb. local dof "
       << p_miit->getNbLocalDofsRow() << std::endl;
    if (verb > VERBOSE) {
      NumeredDofEntity_multiIndex::iterator miit_dd_col =
          p_miit->numeredDofsCols->begin();
      for (; miit_dd_col != p_miit->numeredDofsCols->end(); miit_dd_col++) {
        if ((*miit_dd_col)->pArt == (unsigned int)m_field.get_comm_rank())
          continue;
        if ((*miit_dd_col)->petscLocalDofIdx == (DofIdx)-1)
          continue;
        ss << *(*miit_dd_col) << std::endl;
      }
    }
    PetscSynchronizedPrintf(m_field.get_comm(), ss.str().c_str());
    PetscSynchronizedFlush(m_field.get_comm(), PETSC_STDOUT);
  }

  cOre.getBuildMoFEM() |= Core::PARTITION_GHOST_DOFS;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ProblemsManager::partitionGhostDofsOnDistributedMesh(const std::string name,
                                                     int verb) {
  MoFEM::Interface &m_field = cOre;
  const Problem_multiIndex *problems_ptr;
  MoFEMFunctionBegin;

  if (!(cOre.getBuildMoFEM() & Core::PARTITION_PROBLEM))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "partition of problem not build");
  if (!(cOre.getBuildMoFEM() & Core::PARTITION_FE))
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "partitions finite elements not build");

  // get problem pointer
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  // find p_miit
  CHKERR m_field.get_problems(&problems_ptr);
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
    if (p_miit->numeredDofsCols == p_miit->numeredDofsRows) {
      loop_size = 1;
    }

    typedef decltype(p_miit->numeredDofsRows) NumbDofTypeSharedPtr;
    NumbDofTypeSharedPtr numered_dofs[] = {p_miit->numeredDofsRows,
                                           p_miit->numeredDofsCols};

    // interate over dofs on rows and dofs on columns
    for (int ss = 0; ss != loop_size; ++ss) {

      // create dofs view by uid
      auto r = numered_dofs[ss]->get<PetscLocalIdx_mi_tag>().equal_range(-1);

      std::vector<NumeredDofEntity_multiIndex::iterator> ghost_idx_view;
      ghost_idx_view.reserve(std::distance(r.first, r.second));
      for (; r.first != r.second; ++r.first)
        ghost_idx_view.emplace_back(numered_dofs[ss]->project<0>(r.first));

      auto cmp = [](auto a, auto b) {
        return (*a)->getGlobalUniqueId() < (*b)->getGlobalUniqueId();
      };
      sort(ghost_idx_view.begin(), ghost_idx_view.end(), cmp);

      // intare over dofs which have negative local index
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
    std::ostringstream ss;
    ss << "partition_ghost_row_dofs: rank = " << m_field.get_comm_rank()
       << " FEs row ghost dofs " << *p_miit << " Nb. row ghost dof "
       << p_miit->getNbGhostDofsRow() << " Nb. local dof "
       << p_miit->getNbLocalDofsRow() << std::endl;
    ss << "partition_ghost_col_dofs: rank = " << m_field.get_comm_rank()
       << " FEs col ghost dofs " << *p_miit << " Nb. col ghost dof "
       << p_miit->getNbGhostDofsCol() << " Nb. local dof "
       << p_miit->getNbLocalDofsCol() << std::endl;
    if (verb > VERBOSE) {
      NumeredDofEntity_multiIndex::iterator miit_dd_col =
          p_miit->numeredDofsCols->begin();
      for (; miit_dd_col != p_miit->numeredDofsCols->end(); miit_dd_col++) {
        if ((*miit_dd_col)->pArt == (unsigned int)m_field.get_comm_rank())
          continue;
        if ((*miit_dd_col)->petscLocalDofIdx == (DofIdx)-1)
          continue;
        ss << *(*miit_dd_col) << std::endl;
      }
    }
    PetscSynchronizedPrintf(m_field.get_comm(), ss.str().c_str());
    PetscSynchronizedFlush(m_field.get_comm(), PETSC_STDOUT);
  }

  cOre.getBuildMoFEM() |= Core::PARTITION_GHOST_DOFS;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::getFEMeshset(const std::string prb_name,
                                             const std::string fe_name,
                                             EntityHandle *meshset) const {
  MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_moab().create_meshset(MESHSET_SET, *meshset);
  CHKERR m_field.get_problem(prb_name, &problem_ptr);
  auto fit =
      problem_ptr->numeredFiniteElements->get<FiniteElement_name_mi_tag>()
          .lower_bound(fe_name);
  auto hi_fe_it =
      problem_ptr->numeredFiniteElements->get<FiniteElement_name_mi_tag>()
          .upper_bound(fe_name);
  std::vector<EntityHandle> fe_vec;
  fe_vec.reserve(std::distance(fit, hi_fe_it));
  for (; fit != hi_fe_it; fit++)
    fe_vec.push_back(fit->get()->getEnt());
  CHKERR m_field.get_moab().add_entities(*meshset, &*fe_vec.begin(),
                                         fe_vec.size());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ProblemsManager::getProblemElementsLayout(const std::string name,
                                          const std::string fe_name,
                                          PetscLayout *layout) const {
  MoFEM::Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(name, &problem_ptr);
  CHKERR problem_ptr->getNumberOfElementsByNameAndPart(PETSC_COMM_WORLD,
                                                       fe_name, layout);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::removeDofsOnEntities(
    const std::string problem_name, const std::string field_name,
    const Range ents, const int lo_coeff, const int hi_coeff, int verb,
    const bool debug) {

  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;

  const Problem *prb_ptr;
  CHKERR m_field.get_problem(problem_name, &prb_ptr);

  decltype(prb_ptr->numeredDofsRows) numered_dofs[2] = {
      prb_ptr->numeredDofsRows, nullptr};
  if (prb_ptr->numeredDofsRows != prb_ptr->numeredDofsCols)
    numered_dofs[1] = prb_ptr->numeredDofsCols;

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

      NumeredDofEntity_it_view_multiIndex dofs_it_view;

      // Set -1 to global and local dofs indices
      for (auto pit = ents.const_pair_begin(); pit != ents.const_pair_end();
           ++pit) {
        auto lo =
            numered_dofs[s]
                ->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>()
                .lower_bound(boost::make_tuple(field_name, pit->first, 0));
        auto hi = numered_dofs[s]
                      ->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>()
                      .lower_bound(boost::make_tuple(field_name, pit->second,
                                                     MAX_DOFS_ON_ENTITY));
        for (; lo != hi; ++lo)
          if ((*lo)->getDofCoeffIdx() >= lo_coeff &&
              (*lo)->getDofCoeffIdx() <= hi_coeff)
            dofs_it_view.emplace_back(numered_dofs[s]->project<0>(lo));
      }

      if (verb >= VERY_NOISY) {
        for (auto &dof : dofs_it_view) {
          std::ostringstream ss;
          ss << **dof;
          PetscSynchronizedPrintf(m_field.get_comm(), "%s\n", ss.str().c_str());
        }
        PetscSynchronizedFlush(m_field.get_comm(), PETSC_STDOUT);
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
        PetscSynchronizedPrintf(
            m_field.get_comm(),
            "Number of DOFs in multi-index %d and to delete %d\n",
            numered_dofs[s]->size(), dofs_it_view.size());

      // erase dofs from problem
      for (auto weak_dit : dosf_weak_view)
        if (auto dit = weak_dit.lock()) {
          numered_dofs[s]->erase(dit->getGlobalUniqueId());
        }

      if (verb >= NOISY)
        PetscSynchronizedPrintf(
            m_field.get_comm(),
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
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Imposible case");
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

      auto concatenate_dofs = [&](auto tag, auto &indices,
                                  const auto local_only) {
        MoFEMFunctionBegin;
        get_indices_by_tag(tag, indices, local_only);

        AO ao;
        if (local_only)
          CHKERR AOCreateMapping(m_field.get_comm(), indices.size(),
                                 &*indices.begin(), PETSC_NULL, &ao);
        else
          CHKERR AOCreateMapping(PETSC_COMM_SELF, indices.size(),
                                 &*indices.begin(), PETSC_NULL, &ao);

        get_indices_by_uid(tag, indices);
        CHKERR AOApplicationToPetsc(ao, indices.size(), &*indices.begin());
        CHKERR AODestroy(&ao);
        MoFEMFunctionReturn(0);
      };

      // set indices index
      auto set_concatinated_indices = [&]() {
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
      CHKERR set_concatinated_indices();

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
    PetscSynchronizedPrintf(
        m_field.get_comm(),
        "removed ents on rank %d from problem %s dofs [ %d / %d  (before %d / "
        "%d) local, %d / %d (before %d / %d) "
        "ghost, %d / %d (before %d / %d) global]\n",
        m_field.get_comm_rank(), prb_ptr->getName().c_str(),
        prb_ptr->getNbLocalDofsRow(), prb_ptr->getNbLocalDofsCol(),
        nb_init_row_dofs, nb_init_col_dofs, prb_ptr->getNbGhostDofsRow(),
        prb_ptr->getNbGhostDofsCol(), nb_init_ghost_row_dofs,
        nb_init_ghost_col_dofs, prb_ptr->getNbDofsRow(),
        prb_ptr->getNbDofsCol(), nb_init_loc_row_dofs, nb_init_loc_col_dofs);
    PetscSynchronizedFlush(m_field.get_comm(), PETSC_STDOUT);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ProblemsManager::markDofs(const std::string problem_name,
                                         RowColData rc, const Range ents,
                                         std::vector<bool> &marker) {

  Interface &m_field = cOre;
  const Problem *problem_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_problem(problem_name, &problem_ptr);
  boost::shared_ptr<NumeredDofEntity_multiIndex> dofs;
  switch (rc) {
  case ROW:
    dofs = problem_ptr->getNumeredDofsRows();
    break;
  case COL:
    dofs = problem_ptr->getNumeredDofsCols();
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "Should be row or column");
  }
  marker.resize(dofs->size());
  marker.clear();
  for(auto p = ents.pair_begin();p!=ents.pair_end(); ++p) {
    auto lo = dofs->get<Ent_mi_tag>().lower_bound(p->first);
    auto hi = dofs->get<Ent_mi_tag>().upper_bound(p->second);
    for (; lo != hi; ++lo)
      marker[(*lo)->getPetscLocalDofIdx()] = true;
  }
  MoFEMFunctionReturn(0);
}

} // MOFEM namespace
