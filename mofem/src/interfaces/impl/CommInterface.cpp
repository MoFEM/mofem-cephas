/** \file CommInterface.cpp
 * \brief Functions for interprocessor communications
 * \mofem_comm
 */

namespace MoFEM {

#ifdef PARMETIS

MoFEMErrorCode MatPartitioningApply_Parmetis_MoFEM(MatPartitioning part,
                                                   IS *partitioning);

#endif // PARMETIS

MoFEMErrorCode
CommInterface::query_interface(boost::typeindex::type_index type_index,
                               UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = const_cast<CommInterface *>(this);
  MoFEMFunctionReturnHot(0);
}

CommInterface::CommInterface(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)), dEbug(false) {}

MoFEMErrorCode CommInterface::synchroniseEntities(
    Range &ents, std::map<int, Range> *received_ents, int verb) {
  MoFEM::Interface &m_field = cOre;
  ParallelComm *pcomm = ParallelComm::get_pcomm(
      &m_field.get_moab(), m_field.get_basic_entity_data_ptr()->pcommID);
  MoFEMFunctionBegin;

  auto get_pstatus = [&](const auto ent) {
    unsigned char pstatus;
    CHK_MOAB_THROW(m_field.get_moab().tag_get_data(pcomm->pstatus_tag(), &ent,
                                                   1, &pstatus),
                   "can not get pstatus");
    return pstatus;
  };

  auto get_sharing_procs = [&](const auto ent, const auto pstatus) {
    std::vector<int> sharing_procs(MAX_SHARING_PROCS, -1);
    if (pstatus & PSTATUS_MULTISHARED) {
      // entity is multi shared
      CHK_MOAB_THROW(m_field.get_moab().tag_get_data(
                         pcomm->sharedps_tag(), &ent, 1, &sharing_procs[0]),
                     "can not ger sharing_procs_ptr");
    } else if (pstatus & PSTATUS_SHARED) {
      // shared
      CHK_MOAB_THROW(m_field.get_moab().tag_get_data(pcomm->sharedp_tag(), &ent,
                                                     1, &sharing_procs[0]),
                     "can not get sharing proc");
    }
    return sharing_procs;
  };

  auto get_sharing_handles = [&](const auto ent, const auto pstatus) {
    std::vector<EntityHandle> sharing_handles(MAX_SHARING_PROCS, 0);
    if (pstatus & PSTATUS_MULTISHARED) {
      // entity is multi shared
      CHK_MOAB_THROW(m_field.get_moab().tag_get_data(
                         pcomm->sharedhs_tag(), &ent, 1, &sharing_handles[0]),
                     "get shared handles");
    } else if (pstatus & PSTATUS_SHARED) {
      // shared
      CHK_MOAB_THROW(m_field.get_moab().tag_get_data(pcomm->sharedh_tag(), &ent,
                                                     1, &sharing_handles[0]),
                     "get sharing handle");
    }
    return sharing_handles;
  };

  // make a buffer entities to send
  std::vector<std::vector<EntityHandle>> sbuffer(m_field.get_comm_size());

  for (auto ent : ents) {

    auto pstatus = get_pstatus(ent);
    if (pstatus) {
      auto sharing_procs = get_sharing_procs(ent, pstatus);
      auto sharing_handles = get_sharing_handles(ent, pstatus);

      if (verb >= NOISY) {
        MOFEM_LOG("SYNC", Sev::noisy) << "pstatus " << std::bitset<8>(pstatus);
      }

      for (int proc = 0; proc < MAX_SHARING_PROCS && -1 != sharing_procs[proc];
           proc++) {
        if (sharing_procs[proc] == -1)
          SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE,
                  "sharing processor not set");

        if (sharing_procs[proc] == m_field.get_comm_rank())
          continue;

        const auto handle_on_sharing_proc = sharing_handles[proc];
        // add entity to send, handle on the other side
        sbuffer[sharing_procs[proc]].push_back(handle_on_sharing_proc);
        if (verb >= NOISY) {
          MOFEM_LOG_C("SYNC", Sev::noisy, "send %lu (%lu) to %d at %d\n", ent,
                      handle_on_sharing_proc, sharing_procs[proc],
                      m_field.get_comm_rank());
        }

        if (!(pstatus & PSTATUS_MULTISHARED))
          break;
      }
    }
  }
  if (verb >= NOISY) {
    MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::noisy);
  }

  int nsends = 0; // number of messages to send
  std::vector<int> sbuffer_lengths(
      m_field.get_comm_size()); // length of the message to proc
  const size_t block_size =
      sizeof(EntityHandle) /
      sizeof(int); // FIXME check if that works on given architecture!
  for (int proc = 0; proc < m_field.get_comm_size(); proc++) {

    if (!sbuffer[proc].empty()) {

      sbuffer_lengths[proc] = sbuffer[proc].size() * block_size;
      nsends++;

    } else {

      sbuffer_lengths[proc] = 0;
    }
  }

  // Make sure it is a PETSc m_field.get_comm()
  MPI_Comm comm;
  CHKERR PetscCommDuplicate(m_field.get_comm(), &comm, NULL);

  std::vector<MPI_Status> status(m_field.get_comm_size());

  // Computes the number of messages a node expects to receive
  int nrecvs; // number of messages received
  CHKERR PetscGatherNumberOfMessages(comm, NULL, &sbuffer_lengths[0], &nrecvs);

  // Computes info about messages that a MPI-node will receive, including
  // (from-id,length) pairs for each message.
  int *onodes;   // list of node-ids from which messages are expected
  int *olengths; // corresponding message lengths
  CHKERR PetscGatherMessageLengths(comm, nsends, nrecvs, &sbuffer_lengths[0],
                                   &onodes, &olengths);

  // Gets a unique new tag from a PETSc communicator. All processors that share
  // the communicator MUST call this routine EXACTLY the same number of times.
  // This tag should only be used with the current objects communicator; do NOT
  // use it with any other MPI communicator.
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
  for (int proc = 0, kk = 0; proc < m_field.get_comm_size(); proc++) {
    if (!sbuffer_lengths[proc])
      continue;                             // no message to send to this proc
    CHKERR MPI_Isend(&(sbuffer[proc])[0],   // buffer to send
                     sbuffer_lengths[proc], // message length
                     MPIU_INT, proc,        // to proc
                     tag, comm, s_waits + kk);
    kk++;
  }

  // Wait for received
  if (nrecvs)
    CHKERR MPI_Waitall(nrecvs, r_waits, &status[0]);

  // Wait for send messages
  if (nsends)
    CHKERR MPI_Waitall(nsends, s_waits, &status[0]);

  if (verb >= VERBOSE) {
    MOFEM_LOG_C("SYNC", Sev::verbose, "Rank %d nb. before ents %u\n",
                m_field.get_comm_rank(), ents.size());
    MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::verbose);
  }

  // synchronise range
  for (int kk = 0; kk < nrecvs; kk++) {

    int len = olengths[kk];
    int *data_from_proc = rbuf[kk];

    for (int ee = 0; ee < len;) {
      EntityHandle ent;
      bcopy(&data_from_proc[ee], &ent, sizeof(EntityHandle));
      ents.insert(ent);

      if (received_ents) {
        (*received_ents)[onodes[kk]].insert(ent);
      }

      ee += block_size;

      if (verb >= VERBOSE) {
        MOFEM_LOG_C("SYNC", Sev::noisy, "received %lu from %d at %d\n", ent,
                    onodes[kk], m_field.get_comm_rank());
      }
    }
  }
  if (verb >= VERBOSE) {
    MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::verbose);
  }

  if (verb >= VERBOSE) {
    MOFEM_LOG_C("SYNC", Sev::verbose, "Rank %d nb. after ents %u",
                m_field.get_comm_rank(), ents.size());
    MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::verbose);
  }

  // Cleaning
  CHKERR PetscFree(s_waits);
  CHKERR PetscFree(rbuf[0]);
  CHKERR PetscFree(rbuf);
  CHKERR PetscFree(r_waits);
  CHKERR PetscFree(onodes);
  CHKERR PetscFree(olengths);
  CHKERR PetscCommDestroy(&comm);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CommInterface::synchroniseEntities(Range &ents, int verb) {
  return synchroniseEntities(ents, nullptr, verb);
}

MoFEMErrorCode CommInterface::synchroniseFieldEntities(const std::string name,
                                                       int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  EntityHandle idm = m_field.get_field_meshset(name);
  Range ents;
  CHKERR m_field.get_moab().get_entities_by_handle(idm, ents, false);
  CHKERR synchroniseEntities(ents, nullptr, verb);
  CHKERR m_field.get_moab().add_entities(idm, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CommInterface::resolveParentEntities(const Range &ents,
                                                    int verb) {
  MoFEM::Interface &m_field = cOre;
  ParallelComm *pcomm = ParallelComm::get_pcomm(
      &m_field.get_moab(), m_field.get_basic_entity_data_ptr()->pcommID);

  MoFEMFunctionBegin;

  Range shared = ents;
  CHKERR pcomm->filter_pstatus(shared, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1,
                               nullptr);
  CHKERR pcomm->filter_pstatus(shared, PSTATUS_SHARED | PSTATUS_MULTISHARED,
                               PSTATUS_OR, -1, nullptr);

  auto th_RefParentHandle = cOre.get_th_RefParentHandle();
  auto th_RefBitLevel = cOre.get_th_RefBitLevel();

  auto get_pstatus = [&](const auto ent) {
    unsigned char pstatus;
    CHK_MOAB_THROW(m_field.get_moab().tag_get_data(pcomm->pstatus_tag(), &ent,
                                                   1, &pstatus),
                   "can not get pstatus");
    return pstatus;
  };

  auto get_sharing_procs = [&](const auto ent, const auto pstatus) {
    std::vector<int> sharing_procs(MAX_SHARING_PROCS, -1);
    if (pstatus & PSTATUS_MULTISHARED) {
      // entity is multi shared
      CHK_MOAB_THROW(m_field.get_moab().tag_get_data(
                         pcomm->sharedps_tag(), &ent, 1, &sharing_procs[0]),
                     "can not ger sharing_procs_ptr");
    } else if (pstatus & PSTATUS_SHARED) {
      // shared
      CHK_MOAB_THROW(m_field.get_moab().tag_get_data(pcomm->sharedp_tag(), &ent,
                                                     1, &sharing_procs[0]),
                     "can not get sharing proc");
    }
    return sharing_procs;
  };

  auto get_sharing_handles = [&](const auto ent, const auto pstatus) {
    std::vector<EntityHandle> sharing_handles(MAX_SHARING_PROCS, 0);
    if (pstatus & PSTATUS_MULTISHARED) {
      // entity is multi shared
      CHK_MOAB_THROW(m_field.get_moab().tag_get_data(
                         pcomm->sharedhs_tag(), &ent, 1, &sharing_handles[0]),
                     "get shared handles");
    } else if (pstatus & PSTATUS_SHARED) {
      // shared
      CHK_MOAB_THROW(m_field.get_moab().tag_get_data(pcomm->sharedh_tag(), &ent,
                                                     1, &sharing_handles[0]),
                     "get sharing handle");
    }
    return sharing_handles;
  };

  auto get_parent_and_bit = [&](const auto ent) {
    EntityHandle parent;
    CHK_MOAB_THROW(
        m_field.get_moab().tag_get_data(th_RefParentHandle, &ent, 1, &parent),
        "get parent");
    BitRefLevel bit;
    CHK_MOAB_THROW(
        m_field.get_moab().tag_get_data(th_RefBitLevel, &ent, 1, &bit),
        "get parent");
    return std::make_pair(parent, bit);
  };

  auto set_parent = [&](const auto ent, const auto parent) {
    return m_field.get_moab().tag_set_data(th_RefParentHandle, &ent, 1,
                                           &parent);
  };

  auto set_bit = [&](const auto ent, const auto bit) {
    return m_field.get_moab().tag_set_data(th_RefBitLevel, &ent, 1, &bit);
  };

  // make a buffer
  std::vector<std::vector<unsigned long long>> sbuffer(m_field.get_comm_size());

  for (auto ent : shared) {

    auto pstatus = get_pstatus(ent);
    auto sharing_procs = get_sharing_procs(ent, pstatus);
    auto sharing_handles = get_sharing_handles(ent, pstatus);
    auto [parent, bit] = get_parent_and_bit(ent);

    if (verb >= NOISY)
      MOFEM_LOG("SYNC", Sev::noisy) << "pstatus " << std::bitset<8>(pstatus);

    if (parent) {
      auto pstatus_parent = get_pstatus(parent);
      auto sharing_procs_parent = get_sharing_procs(parent, pstatus_parent);
      auto sharing_handles_parent = get_sharing_handles(parent, pstatus_parent);

      for (int proc = 0; proc < MAX_SHARING_PROCS && -1 != sharing_procs[proc];
           proc++) {
        if (sharing_procs[proc] == -1)
          SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE,
                  "sharing processor not set");

        if (sharing_procs[proc] != m_field.get_comm_rank()) {

          auto it = std::find(sharing_procs_parent.begin(),
                              sharing_procs_parent.end(), sharing_procs[proc]);
          if (it == sharing_procs_parent.end()) {
            SETERRQ1(
                PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Sharing proc for parent entity can not be found proc = %u",
                sharing_procs[proc]);
          }

          auto handle_on_sharing_proc = sharing_handles[proc];
          auto parent_handle_on_sharing_proc =
              sharing_handles_parent[std::distance(sharing_procs_parent.begin(),
                                                   it)];
          sbuffer[sharing_procs[proc]].push_back(handle_on_sharing_proc);
          sbuffer[sharing_procs[proc]].push_back(parent_handle_on_sharing_proc);
          try {
            sbuffer[sharing_procs[proc]].push_back(bit.to_ullong());
          } catch (std::exception &ex) {
            MOFEM_LOG("SELF", Sev::warning) << ex.what();
            MOFEM_LOG("SELF", Sev::warning)
                << "On " << ent << " "
                << moab::CN::EntityTypeName(type_from_handle(ent));
            MOFEM_LOG("SELF", Sev::warning) << "For bit ref " << bit;
          }
          if (verb >= NOISY)
            MOFEM_LOG_C("SYNC", Sev::noisy, "send %lu (%lu) to %d at %d\n", ent,
                        handle_on_sharing_proc, sharing_procs[proc],
                        m_field.get_comm_rank());

          if (!(pstatus & PSTATUS_MULTISHARED))
            break;
        }
      }
    } else {
      for (int proc = 0; proc < MAX_SHARING_PROCS && -1 != sharing_procs[proc];
           proc++) {
        if (sharing_procs[proc] == -1)
          SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE,
                  "sharing processor not set");

        if (sharing_procs[proc] != m_field.get_comm_rank()) {
          auto handle_on_sharing_proc = sharing_handles[proc];
          sbuffer[sharing_procs[proc]].push_back(handle_on_sharing_proc);
          sbuffer[sharing_procs[proc]].push_back(parent);

          try {
            sbuffer[sharing_procs[proc]].push_back(bit.to_ullong());
          } catch (std::exception &ex) {
            MOFEM_LOG("SELF", Sev::warning) << ex.what();
            MOFEM_LOG("SELF", Sev::warning)
                << "On " << ent << " "
                << moab::CN::EntityTypeName(type_from_handle(ent));
            MOFEM_LOG("SELF", Sev::warning) << "For bit ref " << bit;
          }

          if (verb >= NOISY)
            MOFEM_LOG_C("SYNC", Sev::noisy, "send %lu (%lu) to %d at %d\n", ent,
                        handle_on_sharing_proc, sharing_procs[proc],
                        m_field.get_comm_rank());

          if (!(pstatus & PSTATUS_MULTISHARED))
            break;
        }
      }
    }
  }

  int nsends = 0; // number of messages to send
  std::vector<int> sbuffer_lengths(
      m_field.get_comm_size()); // length of the message to proc

  const size_t block_size = sizeof(unsigned long long) / sizeof(int);
  for (int proc = 0; proc < m_field.get_comm_size(); proc++) {

    if (!sbuffer[proc].empty()) {

      sbuffer_lengths[proc] = sbuffer[proc].size() * block_size;
      nsends++;

    } else {

      sbuffer_lengths[proc] = 0;
    }
  }

  // Make sure it is a PETSc m_field.get_comm()
  MPI_Comm comm;
  CHKERR PetscCommDuplicate(m_field.get_comm(), &comm, NULL);

  std::vector<MPI_Status> status(m_field.get_comm_size());

  // Computes the number of messages a node expects to receive
  int nrecvs; // number of messages received
  CHKERR PetscGatherNumberOfMessages(comm, NULL, &sbuffer_lengths[0], &nrecvs);

  // Computes info about messages that a MPI-node will receive, including
  // (from-id,length) pairs for each message.
  int *onodes;   // list of node-ids from which messages are expected
  int *olengths; // corresponding message lengths
  CHKERR PetscGatherMessageLengths(comm, nsends, nrecvs, &sbuffer_lengths[0],
                                   &onodes, &olengths);

  // Gets a unique new tag from a PETSc communicator. All processors that share
  // the communicator MUST call this routine EXACTLY the same number of times.
  // This tag should only be used with the current objects communicator; do NOT
  // use it with any other MPI communicator.
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
  for (int proc = 0, kk = 0; proc < m_field.get_comm_size(); proc++) {
    if (!sbuffer_lengths[proc])
      continue;                             // no message to send to this proc
    CHKERR MPI_Isend(&(sbuffer[proc])[0],   // buffer to send
                     sbuffer_lengths[proc], // message length
                     MPIU_INT, proc,        // to proc
                     tag, comm, s_waits + kk);
    kk++;
  }

  // Wait for received
  if (nrecvs)
    CHKERR MPI_Waitall(nrecvs, r_waits, &status[0]);

  // Wait for send messages
  if (nsends)
    CHKERR MPI_Waitall(nsends, s_waits, &status[0]);

  if (verb >= VERBOSE) {
    MOFEM_LOG_C("SYNC", Sev::verbose,
                "Rank %d nb. shared to synchronise parents ents %u\n",
                m_field.get_comm_rank(), shared.size());
  }

  // synchronise range
  for (int kk = 0; kk < nrecvs; kk++) {

    int len = olengths[kk];
    int *data_from_proc = rbuf[kk];

    for (int ee = 0; ee < len;) {
      EntityHandle ent;
      EntityHandle parent;
      unsigned long long uulong_bit;
      bcopy(&data_from_proc[ee], &ent, sizeof(EntityHandle));
      ee += block_size;
      bcopy(&data_from_proc[ee], &parent, sizeof(EntityHandle));
      ee += block_size;
      bcopy(&data_from_proc[ee], &uulong_bit, sizeof(unsigned long long));
      ee += block_size;

      CHKERR set_parent(ent, parent);
      CHKERR set_bit(ent, BitRefLevel(uulong_bit));

      if (verb >= VERBOSE) {
        MOFEM_LOG_C("SYNC", Sev::noisy, "received %lu (%lu) from %d at %d\n",
                    ent, parent, onodes[kk], m_field.get_comm_rank());
        MOFEM_LOG("SYNC", Sev::noisy) << "Bit " << BitRefLevel(uulong_bit);
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
  CHKERR PetscCommDestroy(&comm);

  if (verb >= VERBOSE)
    MOFEM_LOG_SYNCHRONISE(m_field.get_comm());

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CommInterface::resolveSharedFiniteElements(
    const Problem *problem_ptr, const std::string &fe_name, int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  ParallelComm *pcomm = ParallelComm::get_pcomm(
      &m_field.get_moab(), m_field.get_basic_entity_data_ptr()->pcommID);
  std::vector<int> shprocs(MAX_SHARING_PROCS, 0);
  std::vector<EntityHandle> shhandles(MAX_SHARING_PROCS, 0);
  Range ents;
  Tag th_gid = m_field.get_moab().globalId_tag();
  PetscLayout layout;
  CHKERR problem_ptr->getNumberOfElementsByNameAndPart(m_field.get_comm(),
                                                       fe_name, &layout);
  int gid, last_gid;
  CHKERR PetscLayoutGetRange(layout, &gid, &last_gid);
  CHKERR PetscLayoutDestroy(&layout);

  const FiniteElement_multiIndex *fes_ptr;
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
    for (; fit != hi_fe_it; ++fit) {

      auto ent = (*fit)->getEnt();
      auto part = (*fit)->getPart();

      ents.insert(ent);
      CHKERR m_field.get_moab().tag_set_data(pcomm->part_tag(), &ent, 1, &part);
      if (part == pcomm->rank()) {
        CHKERR m_field.get_moab().tag_set_data(th_gid, &ent, 1, &gid);
        gid++;
      }
      shprocs.clear();
      shhandles.clear();

      if (pcomm->size() > 1) {

        unsigned char pstatus = 0;
        if (pcomm->rank() != part) {
          pstatus = PSTATUS_NOT_OWNED;
          pstatus |= PSTATUS_GHOST;
        }

        if (pcomm->size() > 2) {
          pstatus |= PSTATUS_SHARED;
          pstatus |= PSTATUS_MULTISHARED;
        } else {
          pstatus |= PSTATUS_SHARED;
        }

        size_t rrr = 0;
        for (size_t rr = 0; rr < pcomm->size(); ++rr) {
          if (rr != pcomm->rank()) {
            shhandles[rrr] = ent;
            shprocs[rrr] = rr;
            ++rrr;
          }
        }
        for (; rrr != pcomm->size(); ++rrr)
          shprocs[rrr] = -1;

        if (pstatus & PSTATUS_SHARED) {
          CHKERR m_field.get_moab().tag_set_data(pcomm->sharedp_tag(), &ent, 1,
                                                 &shprocs[0]);
          CHKERR m_field.get_moab().tag_set_data(pcomm->sharedh_tag(), &ent, 1,
                                                 &shhandles[0]);
        }

        if (pstatus & PSTATUS_MULTISHARED) {
          CHKERR m_field.get_moab().tag_set_data(pcomm->sharedps_tag(), &ent, 1,
                                                 &shprocs[0]);
          CHKERR m_field.get_moab().tag_set_data(pcomm->sharedhs_tag(), &ent, 1,
                                                 &shhandles[0]);
        }
        CHKERR m_field.get_moab().tag_set_data(pcomm->pstatus_tag(), &ent, 1,
                                               &pstatus);
      }
    }
  }

  CHKERR pcomm->exchange_tags(th_gid, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CommInterface::resolveSharedFiniteElements(
    const std::string name, const std::string &fe_name, int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  const Problem *problem_ptr;
  CHKERR m_field.get_problem(name, &problem_ptr);
  CHKERR resolveSharedFiniteElements(problem_ptr, fe_name, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
CommInterface::makeEntitiesMultishared(const EntityHandle *entities,
                                       const int num_entities,
                                       const int owner_proc, int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;

  if (m_field.get_comm_size() > 1) {

    ParallelComm *pcomm = ParallelComm::get_pcomm(
        &m_field.get_moab(), m_field.get_basic_entity_data_ptr()->pcommID);

    Range all_ents_range;
    all_ents_range.insert_list(entities, entities + num_entities);

    auto get_tag = [&]() {
      const int negone = -1;
      Tag th;
      m_field.get_moab().tag_get_handle("TMP_GLOBAL_ID_TAG_NAME", 1,
                                        MB_TYPE_INTEGER, th,
                                        MB_TAG_CREAT | MB_TAG_SPARSE, &negone);
      return th;
    };

    auto delete_tag = [&](auto &&th_gid) {
      MoFEMFunctionBegin;
      CHKERR m_field.get_moab().tag_delete(th_gid);
      MoFEMFunctionReturn(0);
    };

    auto resolve_shared_ents = [&](auto &&th_gid, auto &all_ents_range) {
      auto set_gid = [&](auto &th_gid) {
        std::vector<int> gids(num_entities);
        for (size_t g = 0; g != all_ents_range.size(); ++g)
          gids[g] = g + 1;
        CHKERR m_field.get_moab().tag_set_data(th_gid, all_ents_range,
                                               &*gids.begin());

        return &th_gid;
      };

      auto get_skin_ents = [&](auto &all_ents_range) {
        std::array<Range, 4> proc_ents_skin;
        proc_ents_skin[3] = all_ents_range.subset_by_dimension(3);
        proc_ents_skin[2] = all_ents_range.subset_by_dimension(2);
        proc_ents_skin[1] = all_ents_range.subset_by_dimension(1);
        proc_ents_skin[0] = all_ents_range.subset_by_dimension(0);
        return proc_ents_skin;
      };

      auto resolve_dim = [&](auto &all_ents_range) {
        for (int resolve_dim = 3; resolve_dim >= 0; --resolve_dim) {
          if (all_ents_range.num_of_dimension(resolve_dim))
            return resolve_dim;
        }
        return -1;
      };

      auto get_proc_ent = [&](auto &all_ents_range) {
        Range proc_ent;
        if (m_field.get_comm_rank() == owner_proc)
          proc_ent = all_ents_range;
        return proc_ent;
      };

      auto resolve_shared_ents = [&](auto &&proc_ents, auto &&skin_ents) {
        return pcomm->resolve_shared_ents(
            0, proc_ents, resolve_dim(all_ents_range),
            resolve_dim(all_ents_range), skin_ents.data(), set_gid(th_gid));
      };

      CHKERR resolve_shared_ents(get_proc_ent(all_ents_range),
                                 get_skin_ents(all_ents_range));

      return th_gid;
    };

    CHKERR delete_tag(resolve_shared_ents(get_tag(), all_ents_range));

    if (verb >= NOISY) {

      auto print_owner = [&](const EntityHandle e) {
        MoFEMFunctionBegin;
        int moab_owner_proc;
        EntityHandle moab_owner_handle;
        CHKERR pcomm->get_owner_handle(e, moab_owner_proc, moab_owner_handle);

        unsigned char pstatus = 0;

        CHKERR m_field.get_moab().tag_get_data(pcomm->pstatus_tag(), &e, 1,
                                               &pstatus);

        std::vector<int> shprocs(MAX_SHARING_PROCS, 0);
        std::vector<EntityHandle> shhandles(MAX_SHARING_PROCS, 0);

        CHKERR m_field.get_moab().tag_get_data(pcomm->sharedp_tag(), &e, 1,
                                               &shprocs[0]);
        CHKERR m_field.get_moab().tag_get_data(pcomm->sharedh_tag(), &e, 1,
                                               &shhandles[0]);
        if (pstatus & PSTATUS_MULTISHARED) {
          CHKERR m_field.get_moab().tag_get_data(pcomm->sharedps_tag(), &e, 1,
                                                 &shprocs[0]);
          CHKERR m_field.get_moab().tag_get_data(pcomm->sharedhs_tag(), &e, 1,
                                                 &shhandles[0]);
        }

        std::ostringstream ss;

        ss << "Rank " << m_field.get_comm_rank() << " ";
        if (!(pstatus & PSTATUS_NOT_OWNED))
          ss << "OWNER ";
        if (pstatus & PSTATUS_SHARED)
          ss << "PSTATUS_SHARED ";
        if (pstatus & PSTATUS_MULTISHARED)
          ss << "PSTATUS_MULTISHARED ";

        ss << "owner " << moab_owner_proc << " (" << owner_proc << ") ";

        ss << "shprocs: ";
        for (size_t r = 0; r != m_field.get_comm_size() + 1; ++r)
          ss << shprocs[r] << " ";

        ss << "shhandles: ";
        for (size_t r = 0; r != m_field.get_comm_size() + 1; ++r)
          ss << shhandles[r] << " ";

        ss << std::endl;
        MOFEM_LOG("SYNC", Sev::noisy) << ss.str();
        MOFEM_LOG_SYNCHRONISE(m_field.get_comm());

        MoFEMFunctionReturn(0);
      };

      for (auto e : all_ents_range)
        CHKERR print_owner(e);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CommInterface::makeEntitiesMultishared(Range &entities,
                                                      const int owner_proc,
                                                      int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  if (m_field.get_comm_size() > 1) {
    const int num_ents = entities.size();
    std::vector<EntityHandle> vec_ents(num_ents);
    std::copy(entities.begin(), entities.end(), vec_ents.begin());
    CHKERR makeEntitiesMultishared(&*vec_ents.begin(), num_ents, owner_proc,
                                   verb);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
CommInterface::makeFieldEntitiesMultishared(const std::string field_name,
                                            const int owner_proc, int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  if (m_field.get_comm_size() > 1) {
    EntityHandle field_meshset = m_field.get_field_meshset(field_name);
    std::vector<EntityHandle> field_ents;
    CHKERR m_field.get_moab().get_entities_by_handle(field_meshset, field_ents,
                                                     true);
    CHKERR makeEntitiesMultishared(&*field_ents.begin(), field_ents.size(),
                                   owner_proc, verb);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CommInterface::exchangeFieldData(const std::string field_name,
                                                int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  if (m_field.get_comm_size() > 1) {

    Range exchange_ents_data_verts, exchange_ents_data;

    auto *field_ents = m_field.get_field_ents();
    auto field_bit_number = m_field.get_field_bit_number(field_name);
    auto lo = field_ents->get<Unique_mi_tag>().lower_bound(
        FieldEntity::getLoBitNumberUId(field_bit_number));
    auto hi = field_ents->get<Unique_mi_tag>().lower_bound(
        FieldEntity::getHiBitNumberUId(field_bit_number));

    for (auto it = lo; it != hi; ++it)
      if (

          ((*it)->getPStatus()) &&

          (*it)->getNbDofsOnEnt()

      ) {
        if ((*it)->getEntType() == MBVERTEX)
          exchange_ents_data_verts.insert((*it)->getEnt());
        else
          exchange_ents_data.insert((*it)->getEnt());
      }

    auto field_ptr = m_field.get_field_structure(field_name);
    ParallelComm *pcomm = ParallelComm::get_pcomm(
        &m_field.get_moab(), m_field.get_basic_entity_data_ptr()->pcommID);

    auto exchange = [&](const Range &ents, Tag th) {
      MoFEMFunctionBegin;
      if (!ents.empty()) {
        std::vector<Tag> tags;
        tags.push_back(th);
        CHKERR pcomm->exchange_tags(tags, tags, ents);
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR exchange(exchange_ents_data_verts, field_ptr->th_FieldDataVerts);
    CHKERR exchange(exchange_ents_data, field_ptr->th_FieldData);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
CommInterface::partitionMesh(const Range &ents, const int dim,
                             const int adj_dim, const int n_parts,
                             Tag *th_vertex_weights, Tag *th_edge_weights,
                             Tag *th_part_weights, int verb, const bool debug) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;

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
      MOFEM_LOG("SYNC", Sev::inform)
          << "Finite elements in problem: row lower " << rstart << " row upper "
          << rend << " nb. elems " << nb_elems << " ( " << ents.size() << " )";
      MOFEM_LOG_SYNCHRONISE(m_field.get_comm())
    }
  }

  std::vector<EntityHandle> weight_ents;
  weight_ents.reserve(rend - rstart + 1);

  struct AdjBridge {
    EntityHandle ent;
    std::vector<int> adj;
    AdjBridge(const EntityHandle ent, std::vector<int> &adj)
        : ent(ent), adj(adj) {}
  };

  typedef multi_index_container<
      AdjBridge,
      indexed_by<

          hashed_unique<member<AdjBridge, EntityHandle, &AdjBridge::ent>>

          >>
      AdjBridgeMap;

  auto get_it = [&](auto i) {
    auto it = ents.begin();
    for (; i > 0; --i) {
      if (it == ents.end())
        break;
      ++it;
    }
    return it;
  };

  Range proc_ents;
  proc_ents.insert(get_it(rstart), get_it(rend));
  if (proc_ents.size() != rend - rstart)
    SETERRQ2(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
             "Wrong number of elements in range %d != %d", proc_ents.size(),
             rend - rstart);

  Range all_dim_ents;
  CHKERR m_field.get_moab().get_adjacencies(
      proc_ents, adj_dim, true, all_dim_ents, moab::Interface::UNION);

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
      for (fe_it = proc_ents.begin(), ii = rstart, jj = 0, max_row_size = 0;
           fe_it != proc_ents.end(); ++fe_it, ++ii) {

        if (type_from_handle(*fe_it) == MBENTITYSET) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                  "not yet implemented, don't know what to do for meshset "
                  "element");
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
    MOFEM_LOG("WORLD", Sev::verbose) << "Start";

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

    MOFEM_LOG("WORLD", Sev::verbose) << "End";

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
    Tag part_tag = pcomm->part_tag();
    CHKERR m_field.get_moab().tag_set_data(part_tag, ents, part_number);
    Tag gid_tag = m_field.get_moab().globalId_tag();

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
      if (!tagged_sets.empty())
        CHKERR m_field.get_moab().tag_delete_data(part_tag, tagged_sets);

      if (n_parts > (int)tagged_sets.size()) {
        // too few partition sets - create missing ones
        int num_new = n_parts - tagged_sets.size();
        for (int i = 0; i < num_new; i++) {
          EntityHandle new_set;
          CHKERR m_field.get_moab().create_meshset(MESHSET_SET, new_set);
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
      // write a tag to those sets denoting they're partition sets, with a
      // value of the proc number
      std::vector<int> dum_ids(n_parts);
      for (int i = 0; i < n_parts; i++)
        dum_ids[i] = i;
      CHKERR m_field.get_moab().tag_set_data(part_tag, tagged_sets,
                                             &*dum_ids.begin());
      CHKERR m_field.get_moab().clear_meshset(tagged_sets);

      // get lower dimension entities on each part
      for (int pp = 0; pp != n_parts; pp++) {
        Range dim_ents = parts_ents[pp].subset_by_dimension(dim);
        for (int dd = dim - 1; dd >= 0; dd--) {
          Range adj_ents;
          CHKERR m_field.get_moab().get_adjacencies(
              dim_ents, dd, false, adj_ents, moab::Interface::UNION);
          parts_ents[pp].merge(adj_ents);
        }
      }
      for (int pp = 1; pp != n_parts; pp++) {
        for (int ppp = 0; ppp != pp; ppp++) {
          parts_ents[pp] = subtract(parts_ents[pp], parts_ents[ppp]);
        }
      }

      for (int pp = 0; pp != n_parts; pp++) {
        CHKERR m_field.get_moab().add_entities(tagged_sets[pp], parts_ents[pp]);
      }

      auto set_part = [&]() {
        MoFEMFunctionBegin;
        for (EntityType t = MBEDGE; t != MBENTITYSET; ++t) {
          for (int pp = 0; pp != n_parts; pp++) {
            Range type_ents = parts_ents[pp].subset_by_type(t);
            CHKERR m_field.get_moab().tag_clear_data(part_tag, type_ents, &pp);
          }
        }
        MoFEMFunctionReturn(0);
      };

      auto set_gid = [&]() {
        MoFEMFunctionBegin;
        for (EntityType t = MBVERTEX; t != MBENTITYSET; ++t) {

          void *ptr;
          int count;

          int gid = 1; // moab indexing from 1a
          for (int pp = 0; pp != n_parts; pp++) {
            Range type_ents = parts_ents[pp].subset_by_type(t);

            auto eit = type_ents.begin();
            for (; eit != type_ents.end();) {
              CHKERR m_field.get_moab().tag_iterate(
                  gid_tag, eit, type_ents.end(), count, ptr);
              auto gid_tag_ptr = static_cast<int *>(ptr);
              for (; count > 0; --count) {
                *gid_tag_ptr = gid;
                ++eit;
                ++gid;
                ++gid_tag_ptr;
              }
            }
          }
        }

        MoFEMFunctionReturn(0);
      };

      CHKERR set_part();
      CHKERR set_gid();
    }

    if (debug) {
      if (m_field.get_comm_rank() == 0) {
        for (int rr = 0; rr != n_parts; rr++) {
          ostringstream ss;
          ss << "out_part_" << rr << ".vtk";
          MOFEM_LOG("SELF", Sev::inform) << "Save debug part mesh " << ss.str();
          EntityHandle meshset;
          CHKERR m_field.get_moab().create_meshset(MESHSET_SET, meshset);
          CHKERR m_field.get_moab().add_entities(meshset, parts_ents[rr]);
          CHKERR m_field.get_moab().write_file(ss.str().c_str(), "VTK", "",
                                               &meshset, 1);
          CHKERR m_field.get_moab().delete_entities(&meshset, 1);
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

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CommInterface::loadFileRootProcAllRestDistributed(
    moab::Interface &moab, const char *file_name, int dim,
    LoadFileFun proc_skin_fun, const char *options) {
  MoFEMFunctionBegin;

  CHKERR moab.load_file(file_name, 0, options);
  ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
  if (pcomm == NULL)
    pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
  if (pcomm->size() == 1)
    MoFEMFunctionReturnHot(0);

  Skinner skin(&moab);
  MOFEM_LOG_CHANNEL("WORLD");
  MOFEM_LOG_CHANNEL("SYNC");

  auto print_range_on_procs = [&](const Range &range, std::string name) {
    MOFEM_LOG("SYNC", sev) << name << " on proc [" << pcomm->rank()
                           << "] : " << range.size();
    MOFEM_LOG_SEVERITY_SYNC(PETSC_COMM_WORLD, sev);
  };

  auto save_range_to_file = [&](const Range range, std::string name = "part") {
    if (!debug)
      return;
    int rr = pcomm->rank();
    ostringstream ss;
    ss << "out_" << name << "_" << rr << ".vtk";
    MOFEM_LOG("SYNC", sev) << "Save debug part mesh " << ss.str();
    EntityHandle meshset;
    CHKERR moab.create_meshset(MESHSET_SET, meshset);
    CHKERR moab.add_entities(meshset, range);
    if (!range.empty())
      CHKERR moab.write_file(ss.str().c_str(), "VTK", "", &meshset, 1);
    MOFEM_LOG("SYNC", sev) << "Empty range";
    CHKERR moab.delete_entities(&meshset, 1);
    MOFEM_LOG_SEVERITY_SYNC(PETSC_COMM_WORLD, sev);
  };

  auto get_skin = [&](auto e) {
    Range s;
    CHKERR skin.find_skin(0, e.subset_by_dimension(dim), false, s);
    return s;
  };

  auto get_adj = [&](auto ents, auto dim) {
    Range adj;
    CHKERR moab.get_adjacencies(ents, dim, false, adj, moab::Interface::UNION);
    return adj;
  };

  Range all_ents;
  CHKERR moab.get_entities_by_handle(0, all_ents, false);
  save_range_to_file(all_ents, "all_ents");

  Tag part_tag = pcomm->part_tag();
  Range tagged_sets, proc_ents, off_proc_ents;

  CHKERR moab.get_entities_by_type_and_tag(0, MBENTITYSET, &part_tag, NULL, 1,
                                           tagged_sets, moab::Interface::UNION);
  print_range_on_procs(tagged_sets, "tagged_sets");
  for (auto &mit : tagged_sets) {
    int part;
    CHKERR moab.tag_get_data(part_tag, &mit, 1, &part);
    if (part == pcomm->rank()) {
      CHKERR moab.get_entities_by_handle(mit, proc_ents, true);
    } else {
      CHKERR moab.get_entities_by_handle(mit, off_proc_ents, true);
    }
  }

  print_range_on_procs(proc_ents, "proc_ents");
  save_range_to_file(proc_ents, "proc_ents");

  auto number_ents = [&]() {
    auto gid_tag = moab.globalId_tag();
    void *ptr;
    int count;

    for (auto d = 1; d != dim - 1; ++d) {
      int gid = 1;
      Range ents;
      CHKERR moab.get_entities_by_dimension(0, d, ents, true);
      auto eit = ents.begin();
      for (; eit != ents.end();) {
        CHKERR moab.tag_iterate(gid_tag, eit, ents.end(), count, ptr);
        auto gid_tag_ptr = static_cast<int *>(ptr);
        for (; count > 0; --count) {
          *gid_tag_ptr = gid;
          ++eit;
          ++gid;
          ++gid_tag_ptr;
        }
      }
    }
  };

  auto get_proc_ents_skin = [&]() {
    std::array<Range, 4>
        proc_ents_skin; // stores only entities shared with other processors
    auto all_skin = get_skin(all_ents.subset_by_dimension(dim));
    auto proc_skin = get_skin(proc_ents.subset_by_dimension(dim));
    if (dim == 3) {
      proc_ents_skin[2] = subtract(proc_skin, all_skin);       // faces
      proc_ents_skin[1] = get_adj(proc_ents_skin[dim - 1], 1); // edges
    }
    if (dim == 2) {
      proc_ents_skin[1] = subtract(proc_skin, all_skin); // edges
    }
    CHKERR moab.get_connectivity(proc_ents_skin[dim - 1], proc_ents_skin[0],
                                 false); // vertices
    return proc_ents_skin;
  };

  auto add_post_proc_skin = [&](auto &&proc_ents_skin) {
    CubitMeshSet_multiIndex cubit_meshsets_index; ///< cubit meshsets

    Tag nsTag, ssTag, nsTag_data, ssTag_data, bhTag, bhTag_header;
    CHKERR MeshsetsManager::getTags(moab, nsTag, ssTag, nsTag_data, ssTag_data,
                                    bhTag, bhTag_header, QUIET);

    Range meshsets;
    CHKERR moab.get_entities_by_type(0, MBENTITYSET, meshsets, false);
    for (auto m : meshsets) {
      // check if meshset is cubit meshset
      CubitMeshSets block(moab, m);
      if ((block.cubitBcType & CubitBCType(NODESET | SIDESET | BLOCKSET))
              .any()) {
        auto p = cubit_meshsets_index.insert(block);
        if (!p.second) {
          // blockset/nodeset/sideset set exist, could be created on other
          // processor.
          Range ents;
          CHKERR moab.get_entities_by_handle(m, ents, true);
          CHKERR moab.add_entities(p.first->getMeshset(), ents);
          CHKERR moab.delete_entities(&m, 1);
        }
      }
    }

    for (auto &m : cubit_meshsets_index) {
      MOFEM_LOG("WORLD", sev) << "LoadFile read " << m;
    }
    MOFEM_LOG_SEVERITY_SYNC(PETSC_COMM_WORLD, sev);

    std::vector<const CubitMeshSets *> vec_ptr;
    for (auto &m : cubit_meshsets_index) {
      vec_ptr.push_back(&m);
    }
    MeshsetsManager::sortMeshsets(vec_ptr);

    return proc_skin_fun(std::move(proc_ents_skin), std::move(vec_ptr));
  };

  auto remove_off_proc_ents = [&](auto &all_ents, auto &off_proc_ents,
                                  auto &proc_ents_skin) {
    auto to_remove = off_proc_ents;
    for (auto d = dim - 1; d >= 0; --d) {
      to_remove = subtract(to_remove, proc_ents_skin[d]);
    }

    Range meshsets;
    CHKERR moab.get_entities_by_type(0, MBENTITYSET, meshsets, true);
    for (auto m : meshsets) {
      CHKERR moab.remove_entities(m, to_remove);
    }
    for (int dd = dim; dd > 0; --dd) {
      CHKERR moab.delete_entities(to_remove.subset_by_dimension(dd));
    }
  };

  number_ents(); // number higher order entities
  auto proc_ents_skin = add_post_proc_skin(get_proc_ents_skin());
  for (auto d = dim - 1; d >= 0; --d) {
    save_range_to_file(proc_ents_skin[d], "proc_ents_skin" + std::to_string(d));
  }

  if (pcomm->rank()) {
    remove_off_proc_ents(all_ents, off_proc_ents, proc_ents_skin);
  }

  for (auto d = dim - 1; d >= 0; --d) {
    print_range_on_procs(proc_ents_skin[d],
                         "proc_ents_skin by dim" + std::to_string(d));
  }

  CHKERR pcomm->resolve_shared_ents(0, proc_ents, dim, -1,
                                    proc_ents_skin.data());

  if (debug) {

    auto filter_owners = [&](auto &&skin) {
      Range owner_ents;
      CHKERR pcomm->filter_pstatus(skin, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1,
                                   &owner_ents);
      return owner_ents;
    };

    // Create a tag
    int def_val = -1;
    Tag tag_handle;
    CHKERR moab.tag_get_handle("TestGather", 1, MB_TYPE_INTEGER, tag_handle,
                               MB_TAG_DENSE | MB_TAG_CREAT, &def_val);
    Range gather_ents;
    if (pcomm->rank() > 0) {
      gather_ents = proc_ents.subset_by_dimension(3);
      std::vector<int> vals(gather_ents.size());
      std::for_each(vals.begin(), vals.end(), [](auto &v) { v = std::rand(); });
      CHKERR moab.tag_set_data(tag_handle, gather_ents, vals.data());
    } else {
      gather_ents = off_proc_ents.subset_by_dimension(3);
    }
    auto gid_tag = moab.globalId_tag();
    CHKERR pcomm->gather_data(gather_ents, tag_handle, gid_tag);
  }

  Range all_ents_after;
  CHKERR moab.get_entities_by_handle(0, all_ents_after);
  save_range_to_file(all_ents_after.subset_by_dimension(3),
                     "all_ents_part_after");

  MoFEMFunctionReturn(0);
}

Range CommInterface::getPartEntities(moab::Interface &moab, int part) {

  ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
  if (pcomm == NULL)
    pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
  if (pcomm->size() == 1) {
    Range ents;
    CHK_THROW_MESSAGE(moab.get_entities_by_handle(0, ents, false),
                      "get_entities_by_handle failed");
    return subtract(ents, ents.subset_by_type(MBENTITYSET));
  }

  Range proc_ents;

  auto get_proc_ents = [&]() {
    MoFEMFunctionBegin;

    Tag part_tag = pcomm->part_tag();

    Range tagged_sets;
    CHKERR moab.get_entities_by_type_and_tag(0, MBENTITYSET, &part_tag, NULL, 1,
                                             tagged_sets,
                                             moab::Interface::UNION);
    for (auto &mit : tagged_sets) {
      int part;
      CHKERR moab.tag_get_data(part_tag, &mit, 1, &part);
      if (part == pcomm->rank()) {
        CHKERR moab.get_entities_by_handle(mit, proc_ents, true);
        MoFEMFunctionReturnHot(0);
      }
    }
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Part not found in partitioned mesh");
    MoFEMFunctionReturn(0);
  };

  CHK_THROW_MESSAGE(get_proc_ents(), "get_proc_ents failed");

  return proc_ents;
}

} // namespace MoFEM