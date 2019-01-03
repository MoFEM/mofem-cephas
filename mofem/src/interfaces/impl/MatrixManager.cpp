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
  MoFEMErrorCode createMatArrays(ProblemsByName::iterator p_miit,
                                 const MatType type, std::vector<PetscInt> &i,
                                 std::vector<PetscInt> &j,
                                 const bool no_diagonals = true, int verb = -1);

  template <typename TAG>
  MoFEMErrorCode createMat(const std::string &name, Mat *M, const MatType type,
                           PetscInt **_i, PetscInt **_j, PetscScalar **_v,
                           const bool no_diagonals = true, int verb = -1);

  /** \brief Get element adjacencies
   */
  template <typename TAG>
  MoFEMErrorCode getEntityAdjacenies(
      ProblemsByName::iterator p_miit,
      typename boost::multi_index::index<NumeredDofEntity_multiIndex,
                                         TAG>::type::iterator mit_row,
      boost::shared_ptr<FieldEntity> mofem_ent_ptr,
      std::vector<int> &dofs_col_view, int verb);

  MoFEMErrorCode buildFECol(ProblemsByName::iterator p_miit,
                            boost::shared_ptr<EntFiniteElement> ent_fe_ptr,
                            bool do_cols_prob,
                            boost::shared_ptr<NumeredEntFiniteElement> &fe_ptr);
};

MoFEMErrorCode CreateRowComressedADJMatrix::buildFECol(
    ProblemsByName::iterator p_miit,
    boost::shared_ptr<EntFiniteElement> ent_fe_ptr, bool do_cols_prob,
    boost::shared_ptr<NumeredEntFiniteElement> &fe_ptr) {
  MoFEMFunctionBegin;

  if (!ent_fe_ptr) {
    SETERRQ(cOmm, MOFEM_DATA_INCONSISTENCY,
            "Pointer to EntFiniteElement not given");
  }
  // if element is not part of problem
  if ((ent_fe_ptr->getId() & p_miit->getBitFEId()).none())
    MoFEMFunctionReturnHot(0);

  BitRefLevel prb_bit = p_miit->getBitRefLevel();
  BitRefLevel prb_mask = p_miit->getMaskBitRefLevel();
  BitRefLevel fe_bit = ent_fe_ptr->getBitRefLevel();
  // if entity is not problem refinement level
  if ((fe_bit & prb_mask) != fe_bit)
    MoFEMFunctionReturnHot(0);
  if ((fe_bit & prb_bit) != prb_bit)
    MoFEMFunctionReturnHot(0);

  NumeredEntFiniteElement_multiIndex::iterator fe_it =
      p_miit->numeredFiniteElements.find(ent_fe_ptr->getGlobalUniqueId());

  // Create element if is not there
  if (fe_it == p_miit->numeredFiniteElements.end()) {
    std::pair<NumeredEntFiniteElement_multiIndex::iterator, bool> p;
    p = p_miit->numeredFiniteElements.insert(
        boost::make_shared<NumeredEntFiniteElement>(ent_fe_ptr));
    if (!p.second) {
      SETERRQ(cOmm, MOFEM_DATA_INCONSISTENCY,
              "Data inconsistency, this element should be created");
    }
    fe_it = p.first;
  }
  fe_ptr = *fe_it;

  if (fe_ptr) {

    // Build DOFs on columns
    if (fe_ptr->cols_dofs->empty()) {

      // Get dofs on columns
      NumeredDofEntity_multiIndex_uid_view_ordered cols_view;
      CHKERR fe_ptr->getEntFiniteElement()->getColDofView(
          *(p_miit->numeredDofsCols), cols_view, moab::Interface::UNION);

      // Reserve memory for field  dofs
      boost::shared_ptr<std::vector<FENumeredDofEntity>> dofs_array =
          boost::make_shared<std::vector<FENumeredDofEntity>>();
      fe_ptr->getColDofsSequence() = dofs_array;
      dofs_array->reserve(cols_view.size());
      // Reserve memory for shared pointers now
      std::vector<boost::shared_ptr<FENumeredDofEntity>> dofs_shared_array;
      dofs_shared_array.reserve(dofs_array->size());
      // Create dofs objects
      for (NumeredDofEntity_multiIndex_uid_view_ordered::iterator it =
               cols_view.begin();
           it != cols_view.end(); it++) {
        if (!*it) {
          SETERRQ(cOmm, MOFEM_DATA_INCONSISTENCY, "Null pointer to dof");
        }
        boost::shared_ptr<SideNumber> side_number_ptr;
        side_number_ptr = fe_ptr->getSideNumberPtr(it->get()->getEnt());
        dofs_array->push_back(FENumeredDofEntity(side_number_ptr, *it));
        dofs_shared_array.push_back(boost::shared_ptr<FENumeredDofEntity>(
            dofs_array, &dofs_array->back()));
      }
      // Finally add DoFS to multi-indices
      fe_ptr->cols_dofs->insert(dofs_shared_array.begin(),
                                dofs_shared_array.end());
    }

  } else {
    SETERRQ(cOmm, MOFEM_DATA_INCONSISTENCY,
            "At that point ptr to finite element should be well known");
  }

  MoFEMFunctionReturn(0);
}

template <typename TAG>
MoFEMErrorCode CreateRowComressedADJMatrix::getEntityAdjacenies(
    ProblemsByName::iterator p_miit,
    typename boost::multi_index::index<NumeredDofEntity_multiIndex,
                                       TAG>::type::iterator mit_row,
    boost::shared_ptr<FieldEntity> mofem_ent_ptr,
    std::vector<int> &dofs_col_view, int verb) {
  MoFEMFunctionBegin;

  // check if dofs and columns are the same, i.e. structurally symmetric problem
  bool do_cols_prob = true;
  if (p_miit->numeredDofsRows == p_miit->numeredDofsCols) {
    do_cols_prob = false;
  }

  AdjByEnt::iterator adj_miit, hi_adj_miit;
  adj_miit = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(
      mofem_ent_ptr->getGlobalUniqueId());
  hi_adj_miit = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(
      mofem_ent_ptr->getGlobalUniqueId());

  dofs_col_view.clear();
  for (; adj_miit != hi_adj_miit; adj_miit++) {

    if (adj_miit->byWhat & BYROW) {

      if ((adj_miit->entFePtr->getId() & p_miit->getBitFEId()).none()) {
        // if element is not part of problem
        continue;
      }

      BitRefLevel prb_bit = p_miit->getBitRefLevel();
      BitRefLevel prb_mask = p_miit->getMaskBitRefLevel();
      BitRefLevel fe_bit = adj_miit->entFePtr->getBitRefLevel();
      // if entity is not problem refinement level
      if ((fe_bit & prb_mask) != fe_bit)
        continue;
      if ((fe_bit & prb_bit) != prb_bit)
        continue;
      BitRefLevel dof_bit = mit_row->get()->getBitRefLevel();
      // if entity is not problem refinement level
      if ((fe_bit & dof_bit).none())
        continue;

      boost::shared_ptr<NumeredEntFiniteElement> fe_ptr;
      // get element, if element is not in database build columns dofs
      CHKERR buildFECol(p_miit, adj_miit->entFePtr, do_cols_prob, fe_ptr);

      if (fe_ptr) {
        for (FENumeredDofEntity_multiIndex::iterator vit =
                 fe_ptr.get()->cols_dofs->begin();
             vit != fe_ptr.get()->cols_dofs->end(); vit++) {
          const int idx = TAG::get_index(vit);
          dofs_col_view.push_back(idx);
          // Only check if index is correct
          if (idx < 0) {
            std::ostringstream zz;
            zz << "rank " << rAnk << " ";
            zz << *(*vit) << std::endl;
            SETERRQ(cOmm, PETSC_ERR_ARG_SIZ, zz.str().c_str());
          }
          if (idx >= p_miit->getNbDofsCol()) {
            std::ostringstream zz;
            zz << "rank " << rAnk << " ";
            zz << *(*vit) << std::endl;
            SETERRQ(cOmm, PETSC_ERR_ARG_SIZ, zz.str().c_str());
          }
        }
        if (verb >= NOISY) {
          std::stringstream ss;
          ss << "rank " << rAnk << ":  numeredDofsCols" << std::endl;
          FENumeredDofEntity_multiIndex::iterator dit, hi_dit;
          dit = fe_ptr.get()->cols_dofs->begin();
          hi_dit = fe_ptr.get()->cols_dofs->end();
          for (; dit != hi_dit; dit++) {
            ss << "\t" << **dit << std::endl;
          }
          PetscSynchronizedPrintf(cOmm, "%s", ss.str().c_str());
        }
      } else {
        SETERRQ(cOmm, MOFEM_DATA_INCONSISTENCY,
                "Element should be here, otherwise matrix will have missing "
                "elements");
      }
    }
  }

  MoFEMFunctionReturn(0);
}

template <typename TAG>
MoFEMErrorCode CreateRowComressedADJMatrix::createMatArrays(
    ProblemsByName::iterator p_miit, const MatType type,
    std::vector<PetscInt> &i, std::vector<PetscInt> &j, const bool no_diagonals,
    int verb) {
  MoFEMFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_createMat, 0, 0, 0, 0);

  typedef
      typename boost::multi_index::index<NumeredDofEntity_multiIndex, TAG>::type
          NumeredDofEntitysByIdx;

  // Get multi-indices for rows and columns
  const NumeredDofEntitysByIdx &dofs_row_by_idx =
      p_miit->numeredDofsRows->get<TAG>();
  int nb_dofs_row = p_miit->getNbDofsRow();
  if (nb_dofs_row == 0) {
    SETERRQ1(cOmm, MOFEM_DATA_INCONSISTENCY, "problem <%s> has zero rows",
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
    CHKERR PetscLayoutCreate(cOmm, &layout);
    CHKERR PetscLayoutSetBlockSize(layout, 1);
    CHKERR PetscLayoutSetSize(layout, nb_dofs_row);
    CHKERR PetscLayoutSetUp(layout);
    PetscInt rstart, rend;
    CHKERR PetscLayoutGetRange(layout, &rstart, &rend);
    CHKERR PetscLayoutDestroy(&layout);
    if (verb >= VERBOSE) {
      PetscSynchronizedPrintf(cOmm, "\tcreate_Mat: row lower %d row upper %d\n",
                              rstart, rend);
      PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
    }
    miit_row = dofs_row_by_idx.lower_bound(rstart);
    hi_miit_row = dofs_row_by_idx.lower_bound(rend);
    if (std::distance(miit_row, hi_miit_row) != rend - rstart) {
      SETERRQ4(cOmm, PETSC_ERR_ARG_SIZ,
               "data inconsistency, std::distance(miit_row,hi_miit_row) != rend - "
               "rstart (%d != %d - %d = %d) ",
               std::distance(miit_row, hi_miit_row), rend, rstart, rend - rstart);
    }

  } else {

    // Make sure it is a PETSc cOmm
    CHKERR PetscCommDuplicate(cOmm, &cOmm, NULL);

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
          if (mofem_ent_ptr->getGlobalUniqueId() ==
              (*mit_row)->getFieldEntityPtr()->getGlobalUniqueId()) {
            get_adj_col = false;
          }
        }

        if (get_adj_col) {
          // Get entity adjacencies
          mofem_ent_ptr = (*mit_row)->getFieldEntityPtr();
          CHKERR getEntityAdjacenies<TAG>(p_miit, mit_row, mofem_ent_ptr,
                                          dofs_col_view, verb);
          // Sort, uniqe and resize dofs_col_view
          {
            sort(dofs_col_view.begin(), dofs_col_view.end());
            std::vector<int>::iterator new_end =
                unique(dofs_col_view.begin(), dofs_col_view.end());
            int new_size = std::distance(dofs_col_view.begin(), new_end);
            dofs_col_view.resize(new_size);
          }
          // Add that row. Patterns is that first index is row index, second is
          // size of adjacencies after that follows column adjacencies.
          int owner = (*mit_row)->getOwnerProc();
          dofs_vec[owner].push_back(TAG::get_index(mit_row)); // row index
          dofs_vec[owner].push_back(
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
    CHKERR PetscGatherNumberOfMessages(cOmm, NULL, &dofs_vec_length[0],
                                       &nrecvs);

    // Computes info about messages that a MPI-node will receive, including
    // (from-id,length) pairs for each message.
    int *onodes;   // list of node-ids from which messages are expected
    int *olengths; // corresponding message lengths
    CHKERR PetscGatherMessageLengths(cOmm, nsends, nrecvs, &dofs_vec_length[0],
                                     &onodes, &olengths);

    // Gets a unique new tag from a PETSc communicator.
    int tag;
    CHKERR PetscCommGetNewTag(cOmm, &tag);

    // Allocate a buffer sufficient to hold messages of size specified in
    // olengths. And post Irecvs on these buffers using node info from onodes
    int **rbuf;           // must bee freed by user
    MPI_Request *r_waits; // must bee freed by user

    // rbuf has a pointers to messages. It has size of of nrecvs (number of
    // messages) +1. In the first index a block is allocated,
    // such that rbuf[i] = rbuf[i-1]+olengths[i-1].
    CHKERR PetscPostIrecvInt(cOmm, tag, nrecvs, onodes, olengths, &rbuf,
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
                       tag, cOmm, s_waits + kk);
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
          dit = p_miit->numeredDofsRows->get<PetscGlobalIdx_mi_tag>().find(
              row_idx);
          if (dit ==
              p_miit->numeredDofsRows->get<PetscGlobalIdx_mi_tag>().end()) {
            SETERRQ1(cOmm, MOFEM_DATA_INCONSISTENCY,
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
    // if(strcmp(type,MATMPIADJ)==0) {
    // int idx = TAG::get_index(miit_row);
    // if((*dofs_col_by_idx.find(idx))->getGlobalUniqueId()!=(*miit_row)->getGlobalUniqueId())
    // {
    //   SETERRQ(cOmm,PETSC_ERR_ARG_SIZ,"data inconsistency");
    // }
    // }

    // Get entity adjacencies, no need to repeat that operation for dofs when
    // are on the same entity. For simplicity is assumed that those share the
    // same adjacencies.
    if ((!mofem_ent_ptr)
            ? 1
            : (mofem_ent_ptr->getGlobalUniqueId() !=
               (*miit_row)->getFieldEntityPtr()->getGlobalUniqueId())) {

      if (verb >= NOISY) {
        std::stringstream ss;
        ss << "rank " << rAnk << ": row " << **miit_row << std::endl;
        PetscSynchronizedPrintf(cOmm, "%s", ss.str().c_str());
      }

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
          //   cOmm,MOFEM_DATA_INCONSISTENCY,
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
      if (verb >= NOISY) {
        std::stringstream ss;
        ss << "rank " << rAnk << ": dofs_vec for " << *mofem_ent_ptr
           << std::endl;
        PetscSynchronizedPrintf(cOmm, "%s", ss.str().c_str());
      }
    }

    // Try to be smart reserving memory
    if (j.capacity() < j.size() + dofs_vec.size()) {

      unsigned int nb_nonzero = j.size() + dofs_vec.size();
      unsigned int average_row_fill =
          nb_nonzero / i.size() + nb_nonzero % i.size();
      if (j.capacity() < rows_to_fill * average_row_fill) {
        j.reserve(rows_to_fill * average_row_fill);
      }
    }

    // add indices to compressed matrix
    if (verb >= VERY_VERBOSE) {
      PetscSynchronizedPrintf(cOmm, "rank %d: ", rAnk);
    }
    std::vector<int>::iterator diit, hi_diit;
    diit = dofs_vec.begin();
    hi_diit = dofs_vec.end();
    for (; diit != hi_diit; diit++) {

      if (no_diagonals) {
        if (*diit == TAG::get_index(miit_row)) {
          continue;
        }
      }
      j.push_back(*diit);

      if (verb >= VERY_VERBOSE) {
        PetscSynchronizedPrintf(cOmm, "%d ", *diit);
      }
    }
    if (verb >= VERY_VERBOSE) {
      PetscSynchronizedPrintf(cOmm, "\n", *diit);
    }
  }

  if (verb >= VERY_VERBOSE) {
    PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
  }

  // build adj matrix
  i.push_back(j.size());

  if (strcmp(type, MATMPIADJ) == 0) {

    // Adjacency matrix used to partition problems, f.e. METIS
    if (i.size() - 1 != (unsigned int)nb_loc_row_from_iterators) {
      SETERRQ(cOmm, PETSC_ERR_ARG_SIZ, "data inconsistency");
    }

  } else if (strcmp(type, MATMPIAIJ) == 0) {

    // Compressed MPIADJ matrix
    if (i.size() - 1 != (unsigned int)nb_loc_row_from_iterators) {
      SETERRQ(cOmm, PETSC_ERR_ARG_SIZ, "data inconsistency");
    }
    PetscInt nb_local_dofs_row = p_miit->getNbLocalDofsRow();
    if ((unsigned int)nb_local_dofs_row != i.size() - 1) {
      SETERRQ(cOmm, PETSC_ERR_ARG_SIZ, "data inconsistency");
    }

  } else if (strcmp(type, MATAIJ) == 0) {

    // Sequential compressed ADJ matrix
    if (i.size() - 1 != (unsigned int)nb_loc_row_from_iterators) {
      SETERRQ(cOmm, PETSC_ERR_ARG_SIZ, "data inconsistency");
    }
    PetscInt nb_local_dofs_row = p_miit->getNbLocalDofsRow();
    if ((unsigned int)nb_local_dofs_row != i.size() - 1) {
      SETERRQ(cOmm, PETSC_ERR_ARG_SIZ, "data inconsistency");
    }

  } else {

    SETERRQ(cOmm, PETSC_ERR_ARG_NULL, "not implemented");
  }

  PetscLogEventEnd(MOFEM_EVENT_createMat, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

template <typename TAG>
MoFEMErrorCode CreateRowComressedADJMatrix::createMat(
    const std::string &name, Mat *M, const MatType type, PetscInt **_i,
    PetscInt **_j, PetscScalar **_v, const bool no_diagonals, int verb) {
  MoFEMFunctionBegin;

  ProblemsByName &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = pRoblems_set.find(name);
  if (p_miit == pRoblems_set.end()) {
    SETERRQ1(cOmm, MOFEM_NOT_FOUND,
             "problem < %s > is not found (top tip: check spelling)",
             name.c_str());
  }

  std::vector<PetscInt> i, j;
  CHKERR createMatArrays<TAG>(p_miit, type, i, j, no_diagonals, verb);

  CHKERR PetscMalloc(i.size() * sizeof(PetscInt), _i);
  CHKERR PetscMalloc(j.size() * sizeof(PetscInt), _j);
  copy(i.begin(), i.end(), *_i);
  copy(j.begin(), j.end(), *_j);

  PetscInt nb_row_dofs = p_miit->getNbDofsRow();
  PetscInt nb_col_dofs = p_miit->getNbDofsCol();

  if (strcmp(type, MATMPIADJ) == 0) {

    // Adjacency matrix used to partition problems, f.e. METIS
    CHKERR MatCreateMPIAdj(cOmm, i.size() - 1, nb_col_dofs, *_i, *_j,
                           PETSC_NULL, M);
    CHKERR MatSetOption(*M, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE);

  } else if (strcmp(type, MATMPIAIJ) == 0) {
    if (_v != PETSC_NULL) {
      CHKERR PetscMalloc(j.size() * sizeof(PetscScalar), _v);
    }
    PetscScalar *v = (_v != PETSC_NULL) ? *_v : PETSC_NULL;
    // Compressed MPIADJ matrix
    PetscInt nb_local_dofs_row = p_miit->getNbLocalDofsRow();
    PetscInt nb_local_dofs_col = p_miit->getNbLocalDofsCol();
    CHKERR ::MatCreateMPIAIJWithArrays(cOmm, nb_local_dofs_row,
                                       nb_local_dofs_col, nb_row_dofs,
                                       nb_col_dofs, *_i, *_j, v, M);

  } else if (strcmp(type, MATAIJ) == 0) {
    if (_v != PETSC_NULL) {
      CHKERR PetscMalloc(j.size() * sizeof(PetscScalar), _v);
    }
    PetscScalar *v = (_v != PETSC_NULL) ? *_v : PETSC_NULL;
    // Sequential compressed AIJ matrix
    PetscInt nb_local_dofs_row = p_miit->getNbLocalDofsRow();
    PetscInt nb_local_dofs_col = p_miit->getNbLocalDofsCol();
    CHKERR ::MatCreateSeqAIJWithArrays(cOmm, nb_local_dofs_row,
                                       nb_local_dofs_col, *_i, *_j, v, M);

  } else {

    SETERRQ(cOmm, PETSC_ERR_ARG_NULL, "not implemented");
  }
  // MatView(*M,PETSC_VIEWER_STDOUT_WORLD);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::MatCreateMPIAIJWithArrays(const std::string &name,
                                               Mat *Aij, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  int *_i, *_j;
  CreateRowComressedADJMatrix *core_ptr =
      static_cast<CreateRowComressedADJMatrix *>(const_cast<Core *>(this));
  CHKERR core_ptr->createMat<PetscGlobalIdx_mi_tag>(
      name, Aij, MATMPIAIJ, &_i, &_j, PETSC_NULL, false, verb);
  CHKERR PetscFree(_i);
  CHKERR PetscFree(_j);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::MatCreateMPIAdj_with_Idx_mi_tag(const std::string &name,
                                                     Mat *Adj, int verb) {
  MoFEMFunctionBegin;
  int *i, *j;
  if (verb >= VERY_VERBOSE) {
    PetscPrintf(cOmm, "\tCreate Adj matrix\n");
  }
  CreateRowComressedADJMatrix *core_ptr =
      static_cast<CreateRowComressedADJMatrix *>(const_cast<Core *>(this));
  CHKERR core_ptr->createMat<Idx_mi_tag>(name, Adj, MATMPIADJ, &i, &j,
                                         PETSC_NULL, true, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::MatCreateSeqAIJWithArrays(const std::string &name,
                                               Mat *Aij, PetscInt **i,
                                               PetscInt **j, PetscScalar **v,
                                               int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CreateRowComressedADJMatrix *core_ptr =
      static_cast<CreateRowComressedADJMatrix *>(const_cast<Core *>(this));
  CHKERR core_ptr->createMat<PetscLocalIdx_mi_tag>(name, Aij, MATAIJ, i, j, v,
                                                   false, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::partition_check_matrix_fill_in(const std::string &problem_name,
                                     int row_print, int col_print, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;

  struct TestMatrixFillIn : public FEMethod {
    Interface *mFieldPtr;

    Mat A;

    int rowPrint, colPrint;

    TestMatrixFillIn(Interface *m_field_ptr, Mat a, int row_print,
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

      for (FENumeredDofEntity_multiIndex::iterator cit = colPtr->begin();
           cit != colPtr->end(); cit++) {

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
            boost::make_tuple((*cit)->getFieldEntityPtr()->getGlobalUniqueId(),
                              numeredEntFiniteElementPtr->getGlobalUniqueId()));
        if (ait == adjacenciesPtr->end()) {
          SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "adjacencies data inconsistency");
        } else {
          UId uid = ait->getEntUniqueId();
          if (entitiesPtr->find(uid) == entitiesPtr->end()) {
            SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency");
          }
          if (dofsPtr->find((*cit)->getGlobalUniqueId()) == dofsPtr->end()) {
            SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency");
          }
        }

        if ((*cit)->getEntType() != MBVERTEX) {

          FENumeredDofEntity_multiIndex::index<
              Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,
              hi_dit;
          dit = colPtr->get<Composite_Name_Type_And_Side_Number_mi_tag>()
                    .lower_bound(boost::make_tuple(
                        (*cit)->getName(), (*cit)->getEntType(),
                        (*cit)->sideNumberPtr->side_number));
          hi_dit = colPtr->get<Composite_Name_Type_And_Side_Number_mi_tag>()
                       .upper_bound(boost::make_tuple(
                           (*cit)->getName(), (*cit)->getEntType(),
                           (*cit)->sideNumberPtr->side_number));
          int nb_dofs_on_ent = std::distance(dit, hi_dit);

          int max_order = (*cit)->getMaxOrder();
          if ((*cit)->getNbOfCoeffs() * (*cit)->getOrderNbDofs(max_order) !=
              nb_dofs_on_ent) {
            std::cerr << "Warning: Number of Dofs in Col diffrent than number "
                         "of dofs for given entity order "
                      << (*cit)->getNbOfCoeffs() *
                             (*cit)->getOrderNbDofs(max_order)
                      << " " << nb_dofs_on_ent << std::endl;
          }
        }
      }

      FENumeredDofEntity_multiIndex::iterator rit = rowPtr->begin();
      for (; rit != rowPtr->end(); rit++) {

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
            boost::make_tuple((*rit)->getFieldEntityPtr()->getGlobalUniqueId(),
                              numeredEntFiniteElementPtr->getGlobalUniqueId()));
        if (ait == adjacenciesPtr->end()) {
          std::ostringstream ss;
          ss << *(*rit) << std::endl;
          ss << *numeredEntFiniteElementPtr << std::endl;
          ss << "dof: " << (*rit)->getBitRefLevel() << std::endl;
          ss << "fe: " << numeredEntFiniteElementPtr->getBitRefLevel()
             << std::endl;
          ss << "problem: " << problemPtr->getBitRefLevel() << std::endl;
          ss << "problem mask: " << problemPtr->getMaskBitRefLevel()
             << std::endl;
          PetscPrintf(mFieldPtr->get_comm(), "%s", ss.str().c_str());
          SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "adjacencies data inconsistency");
        } else {
          UId uid = ait->getEntUniqueId();
          if (entitiesPtr->find(uid) == entitiesPtr->end()) {
            SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency");
          }
          if (dofsPtr->find((*rit)->getGlobalUniqueId()) == dofsPtr->end()) {
            SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency");
          }
        }
        int row = (*rit)->getPetscGlobalDofIdx();

        FENumeredDofEntity_multiIndex::iterator cit = colPtr->begin();
        for (; cit != colPtr->end(); cit++) {

          int col = (*cit)->getPetscGlobalDofIdx();

          if (row == rowPrint && col == colPrint) {

            std::ostringstream ss;
            ss << "fe:\n" << *numeredEntFiniteElementPtr << std::endl;
            ss << "row:\n" << *(*rit) << std::endl;
            ss << "col:\n" << *(*cit) << std::endl;

            ss << "fe:\n"
               << numeredEntFiniteElementPtr->getBitRefLevel() << std::endl;
            ss << "row:\n" << (*rit)->getBitRefLevel() << std::endl;
            ss << "col:\n" << (*cit)->getBitRefLevel() << std::endl;

            std::cerr << ss.str() << std::endl;

            // PetscPrintf(mFieldPtr->get_comm(),"%s\n",ss.str().c_str());
          }

          CHKERR MatSetValue(A, row, col, 1, INSERT_VALUES);
        }

        if ((*rit)->getEntType() != MBVERTEX) {

          FENumeredDofEntity_multiIndex::index<
              Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,
              hi_dit;
          dit = rowPtr->get<Composite_Name_Type_And_Side_Number_mi_tag>()
                    .lower_bound(boost::make_tuple(
                        (*rit)->getName(), (*rit)->getEntType(),
                        (*rit)->sideNumberPtr->side_number));
          hi_dit = rowPtr->get<Composite_Name_Type_And_Side_Number_mi_tag>()
                       .upper_bound(boost::make_tuple(
                           (*rit)->getName(), (*rit)->getEntType(),
                           (*rit)->sideNumberPtr->side_number));
          int nb_dofs_on_ent = std::distance(dit, hi_dit);

          int max_order = (*rit)->getMaxOrder();
          if ((*rit)->getNbOfCoeffs() * (*rit)->getOrderNbDofs(max_order) !=
              nb_dofs_on_ent) {
            std::cerr << "Warning: Number of Dofs in Row diffrent than number "
                         "of dofs for given entity order "
                      << (*rit)->getNbOfCoeffs() *
                             (*rit)->getOrderNbDofs(max_order)
                      << " " << nb_dofs_on_ent << std::endl;
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

  Mat A;
  CHKERR MatCreateMPIAIJWithArrays(problem_name, &A);
  CHKERR MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

  if (verb >= VERY_VERBOSE) {
    MatView(A, PETSC_VIEWER_STDOUT_WORLD);
  }

  if (verb >= NOISY) {
    MatView(A, PETSC_VIEWER_DRAW_WORLD);
    std::string wait;
    std::cin >> wait;
  }

  TestMatrixFillIn method(this, A, row_print, col_print);

  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  // find p_miit
  ProblemsByName &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = pRoblems_set.find(problem_name);
  if (p_miit == pRoblems_set.end()) {
    SETERRQ1(cOmm, 1, "problem < %s > not found (top tip: check spelling)",
             problem_name.c_str());
  }
  if (verb >= VERBOSE) {
    PetscPrintf(cOmm, "check problem < %s >\n", problem_name.c_str());
  }

  // Loop all elements in problem and check if assemble is without error
  FiniteElement_multiIndex::iterator fe = finiteElements.begin();
  FiniteElement_multiIndex::iterator hi_fe = finiteElements.end();
  for (; fe != hi_fe; fe++) {

    if (verb >= VERBOSE) {
      PetscPrintf(cOmm, "\tcheck element %s\n", (*fe)->getName().c_str());
    }

    CHKERR loop_finite_elements(problem_name, (*fe)->getName(), method,
                                MF_EXIST, verb);
  }

  CHKERR MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  CHKERR MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  CHKERR MatDestroy(&A);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
MatrixManager::query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMMatrixManager) {
    *iface = const_cast<MatrixManager *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

MatrixManager::MatrixManager(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)) {
  PetscLogEventRegister("MatrixManagerCreateMPIAIJWithArrays", 0,
                        &MOFEM_EVENT_createMPIAIJWithArrays);
  PetscLogEventRegister("MatrixManagerCreateMPIAdjWithArrays", 0,
                        &MOFEM_EVENT_createMPIAdjWithArrays);
  PetscLogEventRegister("MatrixManagerCreateSeqAIJWithArrays", 0,
                        &MOFEM_EVENT_createSeqAIJWithArrays);
  PetscLogEventRegister("MatrixManagerCheckMPIAIJWithArraysMatrixFillIn", 0,
                        &MOFEM_EVENT_checkMPIAIJWithArraysMatrixFillIn);
}
MatrixManager::~MatrixManager() {}

template <>
MoFEMErrorCode MatrixManager::createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>(
    const std::string &name, Mat *Aij, PetscInt **i, PetscInt **j,
    PetscScalar **v, int verb) {
  MoFEMFunctionBegin;

  int *_i, *_j;
  double *_v = v ? *v : PETSC_NULL;
  CreateRowComressedADJMatrix *core_ptr =
      static_cast<CreateRowComressedADJMatrix *>(&cOre);
  CHKERR core_ptr->createMat<PetscGlobalIdx_mi_tag>(name, Aij, MATMPIAIJ, &_i,
                                                    &_j, &_v, false, verb);

  if (i)
    *i = _i;
  else
    CHKERR PetscFree(_i);

  if (j)
    *j = _j;
  else
    CHKERR PetscFree(_j);

  if (v)
    *v = _v;

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode MatrixManager::createMPIAdjWithArrays<Idx_mi_tag>(
    const std::string &name, Mat *Adj, PetscInt **i, PetscInt **j,
    PetscScalar **v, int verb) {
  MoFEMFunctionBegin;

  int *_i, *_j;
  double *_v = v ? *v : PETSC_NULL;
  CreateRowComressedADJMatrix *core_ptr =
      static_cast<CreateRowComressedADJMatrix *>(&cOre);
  CHKERR core_ptr->createMat<Idx_mi_tag>(name, Adj, MATMPIADJ, &_i, &_j, &_v,
                                         true, verb);

  if (v)
    *v = _v;

  if (i)
    *i = _i;

  if (j)
    *j = _j;

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode MatrixManager::createSeqAIJWithArrays<PetscLocalIdx_mi_tag>(
    const std::string &name, Mat *Aij, PetscInt **i, PetscInt **j,
    PetscScalar **v, int verb) {
  MoFEMFunctionBegin;

  int *_i, *_j;
  double *_v = v ? *v : PETSC_NULL;
  CreateRowComressedADJMatrix *core_ptr =
      static_cast<CreateRowComressedADJMatrix *>(&cOre);
  CHKERR core_ptr->createMat<PetscLocalIdx_mi_tag>(name, Aij, MATAIJ, i, j, v,
                                                   false, verb);

  if (i)
    *i = _i;
  else
    CHKERR PetscFree(_i);

  if (j)
    *j = _j;
  else
    CHKERR PetscFree(_j);

  if (v)
    *v = _v;

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode
MatrixManager::checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>(
    const std::string &problem_name, int row_print, int col_print,
    int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;

  struct TestMatrixFillIn : public FEMethod {
    Interface *mFieldPtr;

    Mat A;

    int rowPrint, colPrint;

    TestMatrixFillIn(Interface *m_field_ptr, Mat a, int row_print,
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

      for (FENumeredDofEntity_multiIndex::iterator cit = colPtr->begin();
           cit != colPtr->end(); cit++) {

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
            boost::make_tuple((*cit)->getFieldEntityPtr()->getGlobalUniqueId(),
                              numeredEntFiniteElementPtr->getGlobalUniqueId()));
        if (ait == adjacenciesPtr->end()) {
          SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "adjacencies data inconsistency");
        } else {
          UId uid = ait->getEntUniqueId();
          if (entitiesPtr->find(uid) == entitiesPtr->end()) {
            SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency");
          }
          if (dofsPtr->find((*cit)->getGlobalUniqueId()) == dofsPtr->end()) {
            SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency");
          }
        }

        if ((*cit)->getEntType() != MBVERTEX) {

          FENumeredDofEntity_multiIndex::index<
              Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,
              hi_dit;
          dit = colPtr->get<Composite_Name_Type_And_Side_Number_mi_tag>()
                    .lower_bound(boost::make_tuple(
                        (*cit)->getName(), (*cit)->getEntType(),
                        (*cit)->sideNumberPtr->side_number));
          hi_dit = colPtr->get<Composite_Name_Type_And_Side_Number_mi_tag>()
                       .upper_bound(boost::make_tuple(
                           (*cit)->getName(), (*cit)->getEntType(),
                           (*cit)->sideNumberPtr->side_number));
          int nb_dofs_on_ent = std::distance(dit, hi_dit);

          int max_order = (*cit)->getMaxOrder();
          if ((*cit)->getNbOfCoeffs() * (*cit)->getOrderNbDofs(max_order) !=
              nb_dofs_on_ent) {
            std::cerr << "Warning: Number of Dofs in Col diffrent than number "
                         "of dofs for given entity order "
                      << (*cit)->getNbOfCoeffs() *
                             (*cit)->getOrderNbDofs(max_order)
                      << " " << nb_dofs_on_ent << std::endl;
          }
        }
      }

      FENumeredDofEntity_multiIndex::iterator rit = rowPtr->begin();
      for (; rit != rowPtr->end(); rit++) {

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
            boost::make_tuple((*rit)->getFieldEntityPtr()->getGlobalUniqueId(),
                              numeredEntFiniteElementPtr->getGlobalUniqueId()));
        if (ait == adjacenciesPtr->end()) {
          std::ostringstream ss;
          ss << *(*rit) << std::endl;
          ss << *numeredEntFiniteElementPtr << std::endl;
          ss << "dof: " << (*rit)->getBitRefLevel() << std::endl;
          ss << "fe: " << numeredEntFiniteElementPtr->getBitRefLevel()
             << std::endl;
          ss << "problem: " << problemPtr->getBitRefLevel() << std::endl;
          ss << "problem mask: " << problemPtr->getMaskBitRefLevel()
             << std::endl;
          PetscPrintf(mFieldPtr->get_comm(), "%s", ss.str().c_str());
          SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "adjacencies data inconsistency");
        } else {
          UId uid = ait->getEntUniqueId();
          if (entitiesPtr->find(uid) == entitiesPtr->end()) {
            SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency");
          }
          if (dofsPtr->find((*rit)->getGlobalUniqueId()) == dofsPtr->end()) {
            SETERRQ(mFieldPtr->get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency");
          }
        }
        int row = (*rit)->getPetscGlobalDofIdx();

        FENumeredDofEntity_multiIndex::iterator cit = colPtr->begin();
        for (; cit != colPtr->end(); cit++) {

          int col = (*cit)->getPetscGlobalDofIdx();

          if (row == rowPrint && col == colPrint) {

            std::ostringstream ss;
            ss << "fe:\n" << *numeredEntFiniteElementPtr << std::endl;
            ss << "row:\n" << *(*rit) << std::endl;
            ss << "col:\n" << *(*cit) << std::endl;

            ss << "fe:\n"
               << numeredEntFiniteElementPtr->getBitRefLevel() << std::endl;
            ss << "row:\n" << (*rit)->getBitRefLevel() << std::endl;
            ss << "col:\n" << (*cit)->getBitRefLevel() << std::endl;

            std::cerr << ss.str() << std::endl;

            // PetscPrintf(mFieldPtr->get_comm(),"%s\n",ss.str().c_str());
          }

          CHKERR MatSetValue(A, row, col, 1, INSERT_VALUES);
        }

        if ((*rit)->getEntType() != MBVERTEX) {

          FENumeredDofEntity_multiIndex::index<
              Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,
              hi_dit;
          dit = rowPtr->get<Composite_Name_Type_And_Side_Number_mi_tag>()
                    .lower_bound(boost::make_tuple(
                        (*rit)->getName(), (*rit)->getEntType(),
                        (*rit)->sideNumberPtr->side_number));
          hi_dit = rowPtr->get<Composite_Name_Type_And_Side_Number_mi_tag>()
                       .upper_bound(boost::make_tuple(
                           (*rit)->getName(), (*rit)->getEntType(),
                           (*rit)->sideNumberPtr->side_number));
          int nb_dofs_on_ent = std::distance(dit, hi_dit);

          int max_order = (*rit)->getMaxOrder();
          if ((*rit)->getNbOfCoeffs() * (*rit)->getOrderNbDofs(max_order) !=
              nb_dofs_on_ent) {
            std::cerr << "Warning: Number of Dofs in Row diffrent than number "
                         "of dofs for given entity order "
                      << (*rit)->getNbOfCoeffs() *
                             (*rit)->getOrderNbDofs(max_order)
                      << " " << nb_dofs_on_ent << std::endl;
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
  Mat A;
  CHKERR createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>(
      problem_name, &A, PETSC_NULL, PETSC_NULL, PETSC_NULL, verb);
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
  const Problem_multiIndex *problems_ptr;
  CHKERR m_field.get_problems(&problems_ptr);
  auto &prb_set = problems_ptr->get<Problem_mi_tag>();
  auto p_miit = prb_set.find(problem_name);
  if (p_miit == prb_set.end())
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "problem < %s > not found (top tip: check spelling)",
             problem_name.c_str());

  if (verb >= VERBOSE)
    PetscPrintf(m_field.get_comm(), "check problem < %s >\n",
                problem_name.c_str());

  // loop all elements in problem and check if assemble is without error
  const FiniteElement_multiIndex *fe_ptr;
  CHKERR m_field.get_finite_elements(&fe_ptr);
  for (auto &fe : *fe_ptr) {
    if (verb >= VERBOSE)
      PetscPrintf(m_field.get_comm(), "\tcheck element %s\n",
                  fe->getName().c_str());
    CHKERR m_field.loop_finite_elements(problem_name, fe->getName(), method,
                                        MF_EXIST, verb);
  }

  CHKERR MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  CHKERR MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  CHKERR MatDestroy(&A);

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
