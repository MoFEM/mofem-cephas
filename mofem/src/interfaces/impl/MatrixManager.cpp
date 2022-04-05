/**
 * \brief Create adjacent matrices using different indices

 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>

*/

#define MatrixManagerFunctionBegin                                             \
  MoFEMFunctionBegin;                                                          \
  MOFEM_LOG_CHANNEL("WORLD");                                                  \
  MOFEM_LOG_CHANNEL("SYNC");                                                   \
  MOFEM_LOG_FUNCTION();                                                        \
  MOFEM_LOG_TAG("SYNC", "MatrixManager");                                      \
  MOFEM_LOG_TAG("WORLD", "MatrixManager")

namespace MoFEM {

/** \brief Create compressed matrix

  \note Only function class members are allowed in this class. NO VARIABLES.

  \todo It is obsolete implementation, code should be moved to interface
  MatrixManager.

  \todo While matrix is created is assumed that all entities on element are
  adjacent to each other, in some cases, this is created denser matrices than it
  should be. Some kind of filtering can be added.

  \todo Creation of the block matrices

  \todo Some efficiency improvemnt are possible


  */
struct CreateRowComressedADJMatrix : public Core {

  CreateRowComressedADJMatrix(moab::Interface &moab,
                              MPI_Comm comm = PETSC_COMM_WORLD, int verbose = 1)
      : Core(moab, comm, verbose) {}

  typedef FieldEntityEntFiniteElementAdjacencyMap_multiIndex::index<
      Unique_mi_tag>::type AdjByEnt;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  typedef NumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type
      DofByGlobalPetscIndex;

  /** \brief Create matrix adjacencies

    Depending on TAG type, which some structure used to number dofs, matrix is
    partitioned using part number stored in multi-index, or is partitioned on
    parts based only on index number.

    See: Idx_mi_tag  PetscGlobalIdx_mi_tag and PetscLocalIdx_mi_tag

    */
  template <typename TAG>
  MoFEMErrorCode
  createMatArrays(ProblemsByName::iterator p_miit, const MatType type,
                  std::vector<PetscInt> &i, std::vector<PetscInt> &j,
                  const bool no_diagonals = true, int verb = QUIET) const;

  /** \brief Get element adjacencies
   */
  template <typename TAG>
  MoFEMErrorCode getEntityAdjacenies(
      ProblemsByName::iterator p_miit,
      typename boost::multi_index::index<NumeredDofEntity_multiIndex,
                                         TAG>::type::iterator mit_row,
      boost::shared_ptr<FieldEntity> mofem_ent_ptr,
      std::vector<int> &dofs_col_view, int verb) const;
};

template <typename TAG>
MoFEMErrorCode CreateRowComressedADJMatrix::getEntityAdjacenies(
    ProblemsByName::iterator p_miit,
    typename boost::multi_index::index<NumeredDofEntity_multiIndex,
                                       TAG>::type::iterator mit_row,
    boost::shared_ptr<FieldEntity> mofem_ent_ptr,
    std::vector<int> &dofs_col_view, int verb) const {
  MoFEMFunctionBegin;

  BitRefLevel prb_bit = p_miit->getBitRefLevel();
  BitRefLevel prb_mask = p_miit->getBitRefLevelMask();

  const EmptyFieldBlocks &empty_field_blocks = p_miit->getEmptyFieldBlocks();
  const DofIdx nb_dofs_col = p_miit->getNbDofsCol();

  dofs_col_view.clear();
  for (auto r = entFEAdjacencies.get<Unique_mi_tag>().equal_range(
           mofem_ent_ptr->getLocalUniqueId());
       r.first != r.second; ++r.first) {

    if (r.first->byWhat & BYROW) {

      if ((r.first->entFePtr->getId() & p_miit->getBitFEId()).none()) {
        // if element is not part of problem
        continue;
      }

      const BitRefLevel &fe_bit = r.first->entFePtr->getBitRefLevel();
      // if entity is not problem refinement level
      if ((fe_bit & prb_mask) != fe_bit)
        continue;
      if ((fe_bit & prb_bit).none())
        continue;

      const BitRefLevel &dof_bit = (*mit_row)->getBitRefLevel();
      const bool empty_row_block =
          (empty_field_blocks.first & (*mit_row)->getId()).none();

      // if entity is not problem refinement level
      if ((fe_bit & dof_bit).any()) {

        for (auto &it : r.first->entFePtr->getColFieldEnts()) {
          if (auto e = it.lock()) {
            if (empty_row_block ||
                (empty_field_blocks.second & e->getId()).none()) {

              if (auto cache = e->entityCacheColDofs.lock()) {
                const auto lo = cache->loHi[0];
                const auto hi = cache->loHi[1];
                for (auto vit = lo; vit != hi; ++vit) {

                  const int idx = TAG::get_index(vit);
                  if (PetscLikely(idx >= 0))
                    dofs_col_view.push_back(idx);

#ifndef NDEBUG
                  if (PetscUnlikely(idx >= nb_dofs_col)) {
                    MOFEM_LOG("SELF", Sev::error)
                        << "Problem with dof: " << std::endl
                        << "Rank " << rAnk << " : " << *(*vit);
                    SETERRQ(
                        mofemComm, PETSC_ERR_ARG_SIZ,
                        "Index of dof larger than number of DOFs in column");
                  }
#endif
                }

              } else {
                SETERRQ(mofemComm, MOFEM_DATA_INCONSISTENCY, "Cache not set");
              }
            }
          }
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

template <typename TAG>
MoFEMErrorCode CreateRowComressedADJMatrix::createMatArrays(
    ProblemsByName::iterator p_miit, const MatType type,
    std::vector<PetscInt> &i, std::vector<PetscInt> &j, const bool no_diagonals,
    int verb) const {
  MatrixManagerFunctionBegin;

  PetscLogEventBegin(MOFEM_EVENT_createMat, 0, 0, 0, 0);

  auto cache = boost::make_shared<std::vector<EntityCacheNumeredDofs>>(
      entsFields.size());

  size_t idx = 0;
  for (auto it = entsFields.begin(); it != entsFields.end(); ++it, ++idx) {

    const auto uid = (*it)->getLocalUniqueId();
    auto r = entFEAdjacencies.get<Unique_mi_tag>().equal_range(uid);
    for (auto lo = r.first; lo != r.second; ++lo) {

      if ((lo->getBitFEId() & p_miit->getBitFEId()).any()) {

        const BitRefLevel &prb_bit = p_miit->getBitRefLevel();
        const BitRefLevel &prb_mask = p_miit->getBitRefLevelMask();
        const BitRefLevel &fe_bit = lo->entFePtr->getBitRefLevel();

        // if entity is not problem refinement level
        if ((fe_bit & prb_mask) != fe_bit)
          continue;
        if ((fe_bit & prb_bit).none())
          continue;

        auto dit = p_miit->numeredColDofsPtr->lower_bound(uid);

        decltype(dit) hi_dit;
        if (dit != p_miit->numeredColDofsPtr->end())
          hi_dit = p_miit->numeredColDofsPtr->upper_bound(
              uid | static_cast<UId>(MAX_DOFS_ON_ENTITY - 1));
        else
          hi_dit = dit;

        (*it)->entityCacheColDofs =
            boost::shared_ptr<EntityCacheNumeredDofs>(cache, &((*cache)[idx]));
        (*cache)[idx].loHi = {dit, hi_dit};

        break;
      }
    }
  }

  using NumeredDofEntitysByIdx =
      typename boost::multi_index::index<NumeredDofEntity_multiIndex,
                                         TAG>::type;

  // Get multi-indices for rows and columns
  const NumeredDofEntitysByIdx &dofs_row_by_idx =
      p_miit->numeredRowDofsPtr->get<TAG>();
  int nb_dofs_row = p_miit->getNbDofsRow();
  if (nb_dofs_row == 0) {
    SETERRQ1(mofemComm, MOFEM_DATA_INCONSISTENCY, "problem <%s> has zero rows",
             p_miit->getName().c_str());
  }

  // Get adjacencies form other processors
  std::map<int, std::vector<int>> adjacent_dofs_on_other_parts;

  // If not partitioned set petsc layout for matrix. If partitioned need to get
  // adjacencies form other parts. Note if algebra is only partitioned no need
  // to collect adjacencies form other entities. Those are already on mesh
  // which is assumed that is on each processor the same.
  typename boost::multi_index::index<NumeredDofEntity_multiIndex,
                                     TAG>::type::iterator miit_row,
      hi_miit_row;

  if (TAG::IamNotPartitioned) {

    // Get range of local indices
    PetscLayout layout;
    CHKERR PetscLayoutCreate(mofemComm, &layout);
    CHKERR PetscLayoutSetBlockSize(layout, 1);
    CHKERR PetscLayoutSetSize(layout, nb_dofs_row);
    CHKERR PetscLayoutSetUp(layout);
    PetscInt rstart, rend;
    CHKERR PetscLayoutGetRange(layout, &rstart, &rend);
    CHKERR PetscLayoutDestroy(&layout);

    if (verb >= VERBOSE) {
      MOFEM_LOG("SYNC", Sev::noisy)
          << "row lower " << rstart << " row upper " << rend;
      MOFEM_LOG_SYNCHRONISE(mofemComm)
    }

    miit_row = dofs_row_by_idx.lower_bound(rstart);
    hi_miit_row = dofs_row_by_idx.lower_bound(rend);
    if (std::distance(miit_row, hi_miit_row) != rend - rstart) {
      SETERRQ4(
          mofemComm, PETSC_ERR_ARG_SIZ,
          "data inconsistency, std::distance(miit_row,hi_miit_row) != rend - "
          "rstart (%d != %d - %d = %d) ",
          std::distance(miit_row, hi_miit_row), rend, rstart, rend - rstart);
    }

  } else {

    MPI_Comm comm;
    // // Make sure it is a PETSc mofemComm
    CHKERR PetscCommDuplicate(get_comm(), &comm, NULL);

    // get adjacent nodes on other partitions
    std::vector<std::vector<int>> dofs_vec(sIze);

    boost::shared_ptr<FieldEntity> mofem_ent_ptr;
    std::vector<int> dofs_col_view;

    typename boost::multi_index::index<NumeredDofEntity_multiIndex,
                                       TAG>::type::iterator mit_row,
        hi_mit_row;

    mit_row = dofs_row_by_idx.begin();
    hi_mit_row = dofs_row_by_idx.end();
    for (; mit_row != hi_mit_row; mit_row++) {

      // Shared or multishared row and not owned. Those data should be send to
      // other side.

      // Get entity adjacencies, no need to repeat that operation for dofs when
      // are on the same entity. For simplicity is assumed that those sheered
      // the same adjacencies.
      unsigned char pstatus = (*mit_row)->getPStatus();
      if ((pstatus & PSTATUS_NOT_OWNED) &&
          (pstatus & (PSTATUS_SHARED | PSTATUS_MULTISHARED))) {

        bool get_adj_col = true;
        if (mofem_ent_ptr) {
          if (mofem_ent_ptr->getLocalUniqueId() ==
              (*mit_row)->getFieldEntityPtr()->getLocalUniqueId()) {
            get_adj_col = false;
          }
        }

        if (get_adj_col) {
          // Get entity adjacencies
          mofem_ent_ptr = (*mit_row)->getFieldEntityPtr();
          CHKERR getEntityAdjacenies<TAG>(p_miit, mit_row, mofem_ent_ptr,
                                          dofs_col_view, verb);
          // Sort, unique and resize dofs_col_view
          {
            std::sort(dofs_col_view.begin(), dofs_col_view.end());
            std::vector<int>::iterator new_end =
                std::unique(dofs_col_view.begin(), dofs_col_view.end());
            int new_size = std::distance(dofs_col_view.begin(), new_end);
            dofs_col_view.resize(new_size);
          }
          // Add that row. Patterns is that first index is row index, second is
          // size of adjacencies after that follows column adjacencies.
          int owner = (*mit_row)->getOwnerProc();
          dofs_vec[owner].emplace_back(TAG::get_index(mit_row)); // row index
          dofs_vec[owner].emplace_back(
              dofs_col_view.size()); // nb. of column adjacencies
          // add adjacent cools
          dofs_vec[owner].insert(dofs_vec[owner].end(), dofs_col_view.begin(),
                                 dofs_col_view.end());
        }
      }
    }

    int nsends = 0;                         // number of messages to send
    std::vector<int> dofs_vec_length(sIze); // length of the message to proc
    for (int proc = 0; proc < sIze; proc++) {

      if (!dofs_vec[proc].empty()) {

        dofs_vec_length[proc] = dofs_vec[proc].size();
        nsends++;

      } else {

        dofs_vec_length[proc] = 0;
      }
    }

    std::vector<MPI_Status> status(sIze);

    // Computes the number of messages a node expects to receive
    int nrecvs; // number of messages received
    CHKERR PetscGatherNumberOfMessages(comm, NULL, &dofs_vec_length[0],
                                       &nrecvs);

    // Computes info about messages that a MPI-node will receive, including
    // (from-id,length) pairs for each message.
    int *onodes;   // list of node-ids from which messages are expected
    int *olengths; // corresponding message lengths
    CHKERR PetscGatherMessageLengths(comm, nsends, nrecvs, &dofs_vec_length[0],
                                     &onodes, &olengths);

    // Gets a unique new tag from a PETSc communicator.
    int tag;
    CHKERR PetscCommGetNewTag(comm, &tag);

    // Allocate a buffer sufficient to hold messages of size specified in
    // olengths. And post Irecvs on these buffers using node info from onodes
    int **rbuf;           // must bee freed by user
    MPI_Request *r_waits; // must bee freed by user

    // rbuf has a pointers to messages. It has size of of nrecvs (number of
    // messages) +1. In the first index a block is allocated,
    // such that rbuf[i] = rbuf[i-1]+olengths[i-1].
    CHKERR PetscPostIrecvInt(comm, tag, nrecvs, onodes, olengths, &rbuf,
                             &r_waits);

    MPI_Request *s_waits; // status of sens messages
    CHKERR PetscMalloc1(nsends, &s_waits);

    // Send messages
    for (int proc = 0, kk = 0; proc < sIze; proc++) {
      if (!dofs_vec_length[proc])
        continue;                             // no message to send to this proc
      CHKERR MPI_Isend(&(dofs_vec[proc])[0],  // buffer to send
                       dofs_vec_length[proc], // message length
                       MPIU_INT, proc,        // to proc
                       tag, comm, s_waits + kk);
      kk++;
    }

    // Wait for received
    if (nrecvs) {
      CHKERR MPI_Waitall(nrecvs, r_waits, &status[0]);
    }
    // Wait for send messages
    if (nsends) {
      CHKERR MPI_Waitall(nsends, s_waits, &status[0]);
    }

    for (int kk = 0; kk < nrecvs; kk++) {

      int len = olengths[kk];
      int *data_from_proc = rbuf[kk];

      for (int ii = 0; ii < len;) {

        int row_idx = data_from_proc[ii++];     // get row number
        int nb_adj_dofs = data_from_proc[ii++]; // get nb. of adjacent dofs

        if (debug) {

          DofByGlobalPetscIndex::iterator dit;
          dit = p_miit->numeredRowDofsPtr->get<PetscGlobalIdx_mi_tag>().find(
              row_idx);
          if (dit ==
              p_miit->numeredRowDofsPtr->get<PetscGlobalIdx_mi_tag>().end()) {
            SETERRQ1(get_comm(), MOFEM_DATA_INCONSISTENCY,
                     "dof %d can not be found in problem", row_idx);
          }
        }

        for (int jj = 0; jj < nb_adj_dofs; jj++) {
          adjacent_dofs_on_other_parts[row_idx].push_back(data_from_proc[ii++]);
        }
      }
    }

    // Cleaning
    CHKERR PetscFree(s_waits);
    CHKERR PetscFree(rbuf[0]);
    CHKERR PetscFree(rbuf);
    CHKERR PetscFree(r_waits);
    CHKERR PetscFree(onodes);
    CHKERR PetscFree(olengths);

    miit_row = dofs_row_by_idx.begin();
    hi_miit_row = dofs_row_by_idx.end();

    CHKERR PetscCommDestroy(&comm);
  }

  boost::shared_ptr<FieldEntity> mofem_ent_ptr;
  int row_last_evaluated_idx = -1;

  std::vector<int> dofs_vec;
  std::vector<int> dofs_col_view;

  // loop local rows
  int nb_loc_row_from_iterators = 0;
  unsigned int rows_to_fill = p_miit->getNbLocalDofsRow();
  i.reserve(rows_to_fill + 1);
  for (; miit_row != hi_miit_row; miit_row++) {

    if (!TAG::IamNotPartitioned) {
      if (static_cast<int>((*miit_row)->getPart()) != rAnk)
        continue;
    }
    // This is only for cross-check if everything is ok
    nb_loc_row_from_iterators++;

    // add next row to compressed matrix
    i.push_back(j.size());

    // Get entity adjacencies, no need to repeat that operation for dofs when
    // are on the same entity. For simplicity is assumed that those share the
    // same adjacencies.
    if ((!mofem_ent_ptr)
            ? 1
            : (mofem_ent_ptr->getLocalUniqueId() !=
               (*miit_row)->getFieldEntityPtr()->getLocalUniqueId())) {

      // get entity adjacencies
      mofem_ent_ptr = (*miit_row)->getFieldEntityPtr();
      CHKERR getEntityAdjacenies<TAG>(p_miit, miit_row, mofem_ent_ptr,
                                      dofs_col_view, verb);
      row_last_evaluated_idx = TAG::get_index(miit_row);

      dofs_vec.resize(0);
      // insert dofs_col_view
      dofs_vec.insert(dofs_vec.end(), dofs_col_view.begin(),
                      dofs_col_view.end());

      unsigned char pstatus = (*miit_row)->getPStatus();
      if (pstatus > 0) {
        std::map<int, std::vector<int>>::iterator mit;
        mit = adjacent_dofs_on_other_parts.find(row_last_evaluated_idx);
        if (mit == adjacent_dofs_on_other_parts.end()) {
          // NOTE: Dof can adjacent to other part but no elements are there
          // which use that dof std::cerr << *miit_row << std::endl; SETERRQ1(
          //   get_comm(),MOFEM_DATA_INCONSISTENCY,
          //   "data inconsistency row_last_evaluated_idx = %d",
          //   row_last_evaluated_idx
          // );
        } else {
          dofs_vec.insert(dofs_vec.end(), mit->second.begin(),
                          mit->second.end());
        }
      }

      // sort and make unique
      sort(dofs_vec.begin(), dofs_vec.end());
      std::vector<int>::iterator new_end =
          unique(dofs_vec.begin(), dofs_vec.end());
      int new_size = std::distance(dofs_vec.begin(), new_end);
      dofs_vec.resize(new_size);
    }

    // Try to be smart reserving memory
    if (j.capacity() < j.size() + dofs_vec.size()) {

      // TODO: [CORE-55] Improve algorithm estimating size of compressed matrix
      unsigned int nb_nonzero = j.size() + dofs_vec.size();
      unsigned int average_row_fill = nb_nonzero / i.size() + 1;
      j.reserve(rows_to_fill * average_row_fill);
    }

    auto hi_diit = dofs_vec.end();
    for (auto diit = dofs_vec.begin(); diit != hi_diit; diit++) {

      if (no_diagonals) {
        if (*diit == TAG::get_index(miit_row)) {
          continue;
        }
      }
      j.push_back(*diit);
    }
  }

  // build adj matrix
  i.push_back(j.size());

  if (strcmp(type, MATMPIADJ) == 0) {

    // Adjacency matrix used to partition problems, f.e. METIS
    if (i.size() - 1 != (unsigned int)nb_loc_row_from_iterators) {
      SETERRQ(get_comm(), PETSC_ERR_ARG_SIZ, "data inconsistency");
    }

  } else if (strcmp(type, MATMPIAIJ) == 0) {

    // Compressed MPIADJ matrix
    if (i.size() - 1 != (unsigned int)nb_loc_row_from_iterators) {
      SETERRQ(get_comm(), PETSC_ERR_ARG_SIZ, "data inconsistency");
    }
    PetscInt nb_local_dofs_row = p_miit->getNbLocalDofsRow();
    if ((unsigned int)nb_local_dofs_row != i.size() - 1) {
      SETERRQ(get_comm(), PETSC_ERR_ARG_SIZ, "data inconsistency");
    }

  } else if (strcmp(type, MATAIJ) == 0) {

    // Sequential compressed ADJ matrix
    if (i.size() - 1 != (unsigned int)nb_loc_row_from_iterators) {
      SETERRQ(get_comm(), PETSC_ERR_ARG_SIZ, "data inconsistency");
    }
    PetscInt nb_local_dofs_row = p_miit->getNbLocalDofsRow();
    if ((unsigned int)nb_local_dofs_row != i.size() - 1) {
      SETERRQ(get_comm(), PETSC_ERR_ARG_SIZ, "data inconsistency");
    }

  } else {

    SETERRQ(get_comm(), PETSC_ERR_ARG_NULL, "not implemented");
  }

  PetscLogEventEnd(MOFEM_EVENT_createMat, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
MatrixManager::query_interface(boost::typeindex::type_index type_index,
                               UnknownInterface **iface) const {
  *iface = const_cast<MatrixManager *>(this);
  return 0;
}

MatrixManager::MatrixManager(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)) {
  PetscLogEventRegister("MatrixManagerCreateMPIAIJ", 0,
                        &MOFEM_EVENT_createMPIAIJ);
  PetscLogEventRegister("MatrixManagerCreateMPIAIJWithArrays", 0,
                        &MOFEM_EVENT_createMPIAIJWithArrays);
  PetscLogEventRegister("MatrixManagerCreateMPIAdjWithArrays", 0,
                        &MOFEM_EVENT_createMPIAdjWithArrays);
  PetscLogEventRegister("MatrixManagerCreateSeqAIJWithArrays", 0,
                        &MOFEM_EVENT_createSeqAIJWithArrays);
  PetscLogEventRegister("MatrixManagerCheckMPIAIJWithArraysMatrixFillIn", 0,
                        &MOFEM_EVENT_checkMatrixFillIn);
}

template <>
MoFEMErrorCode MatrixManager::createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>(
    const std::string name, Mat *Aij, int verb) {
  MoFEM::CoreInterface &m_field = cOre;
  CreateRowComressedADJMatrix *core_ptr =
      static_cast<CreateRowComressedADJMatrix *>(&cOre);
  MoFEMFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_createMPIAIJWithArrays, 0, 0, 0, 0);

  auto problems_ptr = m_field.get_problems();
  auto &prb = problems_ptr->get<Problem_mi_tag>();
  auto p_miit = prb.find(name);
  if (p_miit == prb.end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_NOT_FOUND,
             "problem < %s > is not found (top tip: check spelling)",
             name.c_str());
  }

  std::vector<int> i_vec, j_vec;
  j_vec.reserve(10000);
  CHKERR core_ptr->createMatArrays<PetscGlobalIdx_mi_tag>(
      p_miit, MATMPIAIJ, i_vec, j_vec, false, verb);

  int nb_row_dofs = p_miit->getNbDofsRow();
  int nb_col_dofs = p_miit->getNbDofsCol();
  int nb_local_dofs_row = p_miit->getNbLocalDofsRow();
  int nb_local_dofs_col = p_miit->getNbLocalDofsCol();

  CHKERR ::MatCreateMPIAIJWithArrays(
      m_field.get_comm(), nb_local_dofs_row, nb_local_dofs_col, nb_row_dofs,
      nb_col_dofs, &*i_vec.begin(), &*j_vec.begin(), PETSC_NULL, Aij);

  PetscLogEventEnd(MOFEM_EVENT_createMPIAIJWithArrays, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode
MatrixManager::createMPIAIJ<PetscGlobalIdx_mi_tag>(const std::string name,
                                                   Mat *Aij, int verb) {
  MoFEM::CoreInterface &m_field = cOre;
  CreateRowComressedADJMatrix *core_ptr =
      static_cast<CreateRowComressedADJMatrix *>(&cOre);
  MoFEMFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_createMPIAIJ, 0, 0, 0, 0);

  auto problems_ptr = m_field.get_problems();
  auto &prb = problems_ptr->get<Problem_mi_tag>();
  auto p_miit = prb.find(name);
  if (p_miit == prb.end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_NOT_FOUND,
             "problem < %s > is not found (top tip: check spelling)",
             name.c_str());
  }

  std::vector<int> i_vec, j_vec;
  j_vec.reserve(10000);
  CHKERR core_ptr->createMatArrays<PetscGlobalIdx_mi_tag>(
      p_miit, MATMPIAIJ, i_vec, j_vec, false, verb);

  int nb_row_dofs = p_miit->getNbDofsRow();
  int nb_col_dofs = p_miit->getNbDofsCol();
  int nb_local_dofs_row = p_miit->getNbLocalDofsRow();
  int nb_local_dofs_col = p_miit->getNbLocalDofsCol();

  auto get_layout = [&]() {
    int start_ranges, end_ranges;
    PetscLayout layout;
    CHKERR PetscLayoutCreate(m_field.get_comm(), &layout);
    CHKERR PetscLayoutSetBlockSize(layout, 1);
    CHKERR PetscLayoutSetLocalSize(layout, nb_local_dofs_col);
    CHKERR PetscLayoutSetUp(layout);
    CHKERR PetscLayoutGetRange(layout, &start_ranges, &end_ranges);
    CHKERR PetscLayoutDestroy(&layout);
    return std::make_pair(start_ranges, end_ranges);
  };

  auto get_nnz = [&](auto &d_nnz, auto &o_nnz) {
    MoFEMFunctionBeginHot;
    auto layout = get_layout();
    int j = 0;
    for (int i = 0; i != nb_local_dofs_row; ++i) {
      for (; j != i_vec[i + 1]; ++j) {
        if (j_vec[j] < layout.second && j_vec[j] >= layout.first)
          ++(d_nnz[i]);
        else
          ++(o_nnz[i]);
      }
    }
    MoFEMFunctionReturnHot(0);
  };

  std::vector<int> d_nnz(nb_local_dofs_row, 0), o_nnz(nb_local_dofs_row, 0);
  CHKERR get_nnz(d_nnz, o_nnz);

  CHKERR MatCreate(m_field.get_comm(), Aij);
  CHKERR MatSetSizes(*Aij, nb_local_dofs_row, nb_local_dofs_col, nb_row_dofs,
                     nb_col_dofs);
  CHKERR MatSetType(*Aij, MATMPIAIJ);
  CHKERR MatMPIAIJSetPreallocation(*Aij, 0, &*d_nnz.begin(), 0,
                                   &*o_nnz.begin());

  PetscLogEventEnd(MOFEM_EVENT_createMPIAIJ, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode
MatrixManager::createMPIAdjWithArrays<Idx_mi_tag>(const std::string name,
                                                  Mat *Adj, int verb) {
  MoFEM::CoreInterface &m_field = cOre;
  CreateRowComressedADJMatrix *core_ptr =
      static_cast<CreateRowComressedADJMatrix *>(&cOre);
  MoFEMFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_createMPIAdjWithArrays, 0, 0, 0, 0);

  auto problems_ptr = m_field.get_problems();
  auto &prb = problems_ptr->get<Problem_mi_tag>();
  auto p_miit = prb.find(name);
  if (p_miit == prb.end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_NOT_FOUND,
             "problem < %s > is not found (top tip: check spelling)",
             name.c_str());
  }

  std::vector<int> i_vec, j_vec;
  j_vec.reserve(10000);
  CHKERR core_ptr->createMatArrays<Idx_mi_tag>(p_miit, MATMPIADJ, i_vec, j_vec,
                                               true, verb);
  int *_i, *_j;
  CHKERR PetscMalloc(i_vec.size() * sizeof(int), &_i);
  CHKERR PetscMalloc(j_vec.size() * sizeof(int), &_j);
  copy(i_vec.begin(), i_vec.end(), _i);
  copy(j_vec.begin(), j_vec.end(), _j);

  int nb_col_dofs = p_miit->getNbDofsCol();
  CHKERR MatCreateMPIAdj(m_field.get_comm(), i_vec.size() - 1, nb_col_dofs, _i,
                         _j, PETSC_NULL, Adj);
  CHKERR MatSetOption(*Adj, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE);

  PetscLogEventEnd(MOFEM_EVENT_createMPIAdjWithArrays, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode MatrixManager::createSeqAIJWithArrays<PetscLocalIdx_mi_tag>(
    const std::string name, Mat *Aij, int verb) {
  MoFEM::CoreInterface &m_field = cOre;
  CreateRowComressedADJMatrix *core_ptr =
      static_cast<CreateRowComressedADJMatrix *>(&cOre);
  MoFEMFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_createMPIAIJWithArrays, 0, 0, 0, 0);

  auto problems_ptr = m_field.get_problems();
  auto &prb = problems_ptr->get<Problem_mi_tag>();
  auto p_miit = prb.find(name);
  if (p_miit == prb.end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_NOT_FOUND,
             "problem < %s > is not found (top tip: check spelling)",
             name.c_str());
  }

  std::vector<int> i_vec, j_vec;
  j_vec.reserve(10000);
  CHKERR core_ptr->createMatArrays<PetscGlobalIdx_mi_tag>(p_miit, MATAIJ, i_vec,
                                                          j_vec, false, verb);

  int nb_local_dofs_row = p_miit->getNbLocalDofsRow();
  int nb_local_dofs_col = p_miit->getNbLocalDofsCol();

  double *_a;
  CHKERR PetscMalloc(j_vec.size() * sizeof(double), &_a);

  Mat tmpMat;
  CHKERR ::MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, nb_local_dofs_row,
                                     nb_local_dofs_col, &*i_vec.begin(),
                                     &*j_vec.begin(), _a, &tmpMat);
  CHKERR MatDuplicate(tmpMat, MAT_SHARE_NONZERO_PATTERN, Aij);
  CHKERR MatDestroy(&tmpMat);

  CHKERR PetscFree(_a);

  PetscLogEventEnd(MOFEM_EVENT_createMPIAIJWithArrays, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MatrixManager::checkMatrixFillIn(const std::string problem_name,
                                                int row_print, int col_print,
                                                Mat A, int verb) {
  MoFEM::CoreInterface &m_field = cOre;
  MatrixManagerFunctionBegin;

  PetscLogEventBegin(MOFEM_EVENT_checkMatrixFillIn, 0, 0, 0, 0);

  struct TestMatrixFillIn : public FEMethod {
    CoreInterface *mFieldPtr;

    Mat A;

    int rowPrint, colPrint;

    TestMatrixFillIn(CoreInterface *m_field_ptr, Mat a, int row_print,
                     int col_print)
        : mFieldPtr(m_field_ptr), A(a), rowPrint(row_print),
          colPrint(col_print){};

    MoFEMErrorCode preProcess() {
      MoFEMFunctionBeginHot;
      MoFEMFunctionReturnHot(0);
    }

    MoFEMErrorCode operator()() {
      MoFEMFunctionBegin;

      if (refinedFiniteElementsPtr->find(
              numeredEntFiniteElementPtr->getEnt()) ==
          refinedFiniteElementsPtr->end()) {
        SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }

      auto row_dofs = getRowDofsPtr();
      auto col_dofs = getColDofsPtr();

      for (auto cit = col_dofs->begin(); cit != col_dofs->end(); cit++) {

        if (refinedEntitiesPtr->find((*cit)->getEnt()) ==
            refinedEntitiesPtr->end()) {
          SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "data inconsistency");
        }
        if (!(*cit)->getActive()) {
          SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "data inconsistency");
        }

        FieldEntityEntFiniteElementAdjacencyMap_multiIndex::index<
            Composite_Unique_mi_tag>::type::iterator ait;
        ait = adjacenciesPtr->get<Composite_Unique_mi_tag>().find(
            boost::make_tuple((*cit)->getFieldEntityPtr()->getLocalUniqueId(),
                              numeredEntFiniteElementPtr->getLocalUniqueId()));
        if (ait == adjacenciesPtr->end()) {
          SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "adjacencies data inconsistency");
        } else {
          UId uid = ait->getEntUniqueId();
          if (entitiesPtr->find(uid) == entitiesPtr->end()) {
            SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency");
          }
          if (dofsPtr->find((*cit)->getLocalUniqueId()) == dofsPtr->end()) {
            SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency");
          }
        }

        if ((*cit)->getEntType() != MBVERTEX) {

          auto range =
              col_dofs->get<Ent_mi_tag>().equal_range((*cit)->getEnt());
          int nb_dofs_on_ent = std::distance(range.first, range.second);

          int max_order = (*cit)->getMaxOrder();
          if ((*cit)->getNbOfCoeffs() * (*cit)->getOrderNbDofs(max_order) !=
              nb_dofs_on_ent) {
            MOFEM_LOG("SELF", Sev::warning)
                << "Warning: Number of Dofs in Col diffrent than number "
                   "of dofs for given entity order "
                << (*cit)->getNbOfCoeffs() * (*cit)->getOrderNbDofs(max_order)
                << " " << nb_dofs_on_ent;
          }
        }
      }

      for (auto rit = row_dofs->begin(); rit != row_dofs->end(); rit++) {

        if (refinedEntitiesPtr->find((*rit)->getEnt()) ==
            refinedEntitiesPtr->end()) {
          SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "data inconsistency");
        }
        if (!(*rit)->getActive()) {
          SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "data inconsistency");
        }

        FieldEntityEntFiniteElementAdjacencyMap_multiIndex::index<
            Composite_Unique_mi_tag>::type::iterator ait;
        ait = adjacenciesPtr->get<Composite_Unique_mi_tag>().find(
            boost::make_tuple((*rit)->getFieldEntityPtr()->getLocalUniqueId(),
                              numeredEntFiniteElementPtr->getLocalUniqueId()));
        if (ait == adjacenciesPtr->end()) {
          MOFEM_LOG_ATTRIBUTES("SYNC", LogManager::BitScope);
          MOFEM_LOG("SELF", Sev::error) << *(*rit);
          MOFEM_LOG("SELF", Sev::error) << *(*rit);
          MOFEM_LOG("SELF", Sev::error) << *numeredEntFiniteElementPtr;
          MOFEM_LOG("SELF", Sev::error) << "dof: " << (*rit)->getBitRefLevel();
          MOFEM_LOG("SELF", Sev::error)
              << "fe: " << numeredEntFiniteElementPtr->getBitRefLevel();
          MOFEM_LOG("SELF", Sev::error)
              << "problem: " << problemPtr->getBitRefLevel();
          MOFEM_LOG("SELF", Sev::error)
              << "problem mask: " << problemPtr->getBitRefLevelMask();
          SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "adjacencies data inconsistency");
        } else {
          UId uid = ait->getEntUniqueId();
          if (entitiesPtr->find(uid) == entitiesPtr->end()) {
            SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency");
          }
          if (dofsPtr->find((*rit)->getLocalUniqueId()) == dofsPtr->end()) {
            SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency");
          }
        }
        int row = (*rit)->getPetscGlobalDofIdx();

        auto col_dofs = getColDofsPtr();
        for (auto cit = col_dofs->begin(); cit != col_dofs->end(); cit++) {

          int col = (*cit)->getPetscGlobalDofIdx();

          if (row == rowPrint && col == colPrint) {
            MOFEM_LOG("SELF", Sev::noisy) << "fe:\n"
                                          << *numeredEntFiniteElementPtr;
            MOFEM_LOG("SELF", Sev::noisy) << "row:\n" << *(*rit);
            MOFEM_LOG("SELF", Sev::noisy) << "col:\n" << *(*cit);
            MOFEM_LOG("SELF", Sev::noisy)
                << "fe:\n"
                << numeredEntFiniteElementPtr->getBitRefLevel();
            MOFEM_LOG("SELF", Sev::noisy) << "row:\n"
                                          << (*rit)->getBitRefLevel();
            MOFEM_LOG("SELF", Sev::noisy) << "col:\n"
                                          << (*cit)->getBitRefLevel();
          }

          CHKERR MatSetValue(A, row, col, 1, INSERT_VALUES);
        }

        if ((*rit)->getEntType() != MBVERTEX) {

          auto range =
              row_dofs->get<Ent_mi_tag>().equal_range((*rit)->getEnt());
          int nb_dofs_on_ent = std::distance(range.first, range.second);

          int max_order = (*rit)->getMaxOrder();
          if ((*rit)->getNbOfCoeffs() * (*rit)->getOrderNbDofs(max_order) !=
              nb_dofs_on_ent) {
            MOFEM_LOG("SELF", Sev::warning)
                << "Warning: Number of Dofs in Row diffrent than number "
                   "of dofs for given entity order "
                << (*rit)->getNbOfCoeffs() * (*rit)->getOrderNbDofs(max_order)
                << " " << nb_dofs_on_ent;
          }
        }
      }

      MoFEMFunctionReturn(0);
    }

    MoFEMErrorCode postProcess() {
      MoFEMFunctionBegin;

      // cerr << mFieldPtr->get_comm_rank() << endl;
      CHKERR MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY);
      CHKERR MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY);

      MoFEMFunctionReturn(0);
    }
  };

  // create matrix
  CHKERR MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

  if (verb >= VERY_VERBOSE) {
    MatView(A, PETSC_VIEWER_STDOUT_WORLD);
  }

  if (verb >= NOISY) {
    MatView(A, PETSC_VIEWER_DRAW_WORLD);
    std::string wait;
    std::cin >> wait;
  }

  TestMatrixFillIn method(&m_field, A, row_print, col_print);

  // get problem
  auto problems_ptr = m_field.get_problems();
  auto &prb_set = problems_ptr->get<Problem_mi_tag>();
  auto p_miit = prb_set.find(problem_name);
  if (p_miit == prb_set.end())
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "problem < %s > not found (top tip: check spelling)",
             problem_name.c_str());
  MOFEM_LOG_C("WORLD", Sev::inform, "check problem < %s >",
              problem_name.c_str());

  // loop all elements in problem and check if assemble is without error
  auto fe_ptr = m_field.get_finite_elements();
  for (auto &fe : *fe_ptr) {
    MOFEM_LOG_C("WORLD", Sev::verbose, "\tcheck element %s",
                fe->getName().c_str());
    CHKERR m_field.loop_finite_elements(problem_name, fe->getName(), method,
                                        nullptr, MF_EXIST,
                                        CacheTupleSharedPtr(), verb);
  }

  CHKERR MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  CHKERR MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  PetscLogEventEnd(MOFEM_EVENT_checkMatrixFillIn, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode
MatrixManager::checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>(
    const std::string problem_name, int row_print, int col_print, int verb) {
  MoFEM::CoreInterface &m_field = cOre;
  MatrixManagerFunctionBegin;
  // create matrix
  SmartPetscObj<Mat> A;
  CHKERR createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>(problem_name, A, verb);
  CHKERR MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  CHKERR checkMatrixFillIn(problem_name, row_print, col_print, A, verb);
  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode MatrixManager::checkMPIAIJMatrixFillIn<PetscGlobalIdx_mi_tag>(
    const std::string problem_name, int row_print, int col_print, int verb) {
  MoFEM::CoreInterface &m_field = cOre;
  MatrixManagerFunctionBegin;
  // create matrix
  SmartPetscObj<Mat> A;
  CHKERR createMPIAIJ<PetscGlobalIdx_mi_tag>(problem_name, A, verb);
  CHKERR MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  CHKERR checkMatrixFillIn(problem_name, row_print, col_print, A, verb);
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
