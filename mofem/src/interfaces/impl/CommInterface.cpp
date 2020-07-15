/** \file CommInterface.cpp
 * \brief Functions for interprocessor communications
 * \mofem_comm
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

MoFEMErrorCode CommInterface::query_interface(const MOFEMuuid &uuid,
                                              UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMComm) {
    *iface = const_cast<CommInterface *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

CommInterface::CommInterface(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)), dEbug(false) {}
CommInterface::~CommInterface() {}

MoFEMErrorCode CommInterface::synchroniseEntities(Range &ents, int verb) {
  MoFEM::Interface &m_field = cOre;
  auto ref_ents_ptr = m_field.get_ref_ents();
  MoFEMFunctionBegin;

  // make a buffer
  std::vector<std::vector<EntityHandle>> sbuffer(m_field.get_comm_size());

  Range::iterator eit = ents.begin();
  for (; eit != ents.end(); eit++) {

    auto meit = ref_ents_ptr->get<Ent_mi_tag>().find(*eit);
    if (meit == ref_ents_ptr->get<Ent_mi_tag>().end()) {
      continue;
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "rank %d entity %lu not exist on database, local entity can not "
               "be found for this owner",
               m_field.get_comm_rank(), *eit);
    }

    unsigned char pstatus = (*meit)->getPStatus();

    if (pstatus == 0)
      continue;

    if (verb >= NOISY) {
      std::ostringstream zz;
      zz << "pstatus " << std::bitset<8>(pstatus) << " ";
      PetscSynchronizedPrintf(m_field.get_comm(), "%s", zz.str().c_str());
    }

    for (int proc = 0;
         proc < MAX_SHARING_PROCS && -1 != (*meit)->getSharingProcsPtr()[proc];
         proc++) {
      if ((*meit)->getSharingProcsPtr()[proc] == -1)
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
                "sharing processor not set");

      if ((*meit)->getSharingProcsPtr()[proc] == m_field.get_comm_rank())
        continue;

      EntityHandle handle_on_sharing_proc =
          (*meit)->getSharingHandlersPtr()[proc];
      sbuffer[(*meit)->getSharingProcsPtr()[proc]].push_back(
          handle_on_sharing_proc);
      if (verb >= NOISY)
        PetscSynchronizedPrintf(
            m_field.get_comm(), "send %lu (%lu) to %d at %d\n",
            (*meit)->getEnt(), handle_on_sharing_proc,
            (*meit)->getSharingProcsPtr()[proc], m_field.get_comm_rank());

      if (!(pstatus & PSTATUS_MULTISHARED))
        break;
    }
  }

  int nsends = 0; // number of messages to send
  std::vector<int> sbuffer_lengths(
      m_field.get_comm_size()); // length of the message to proc
  const size_t block_size = sizeof(EntityHandle) / sizeof(int);
  for (int proc = 0; proc < m_field.get_comm_size(); proc++) {

    if (!sbuffer[proc].empty()) {

      sbuffer_lengths[proc] = sbuffer[proc].size() * block_size;
      nsends++;

    } else {

      sbuffer_lengths[proc] = 0;
    }
  }

  // Make sure it is a PETSc m_field.get_comm()
  CHKERR PetscCommDuplicate(m_field.get_comm(), &m_field.get_comm(), NULL);

  std::vector<MPI_Status> status(m_field.get_comm_size());

  // Computes the number of messages a node expects to receive
  int nrecvs; // number of messages received
  CHKERR PetscGatherNumberOfMessages(m_field.get_comm(), NULL,
                                     &sbuffer_lengths[0], &nrecvs);

  // Computes info about messages that a MPI-node will receive, including
  // (from-id,length) pairs for each message.
  int *onodes;   // list of node-ids from which messages are expected
  int *olengths; // corresponding message lengths
  CHKERR PetscGatherMessageLengths(m_field.get_comm(), nsends, nrecvs,
                                   &sbuffer_lengths[0], &onodes, &olengths);

  // Gets a unique new tag from a PETSc communicator. All processors that share
  // the communicator MUST call this routine EXACTLY the same number of times.
  // This tag should only be used with the current objects communicator; do NOT
  // use it with any other MPI communicator.
  int tag;
  CHKERR PetscCommGetNewTag(m_field.get_comm(), &tag);

  // Allocate a buffer sufficient to hold messages of size specified in
  // olengths. And post Irecvs on these buffers using node info from onodes
  int **rbuf;           // must bee freed by user
  MPI_Request *r_waits; // must bee freed by user
  // rbuf has a pointers to messages. It has size of of nrecvs (number of
  // messages) +1. In the first index a block is allocated,
  // such that rbuf[i] = rbuf[i-1]+olengths[i-1].
  CHKERR PetscPostIrecvInt(m_field.get_comm(), tag, nrecvs, onodes, olengths,
                           &rbuf, &r_waits);

  MPI_Request *s_waits; // status of sens messages
  CHKERR PetscMalloc1(nsends, &s_waits);

  // Send messages
  for (int proc = 0, kk = 0; proc < m_field.get_comm_size(); proc++) {
    if (!sbuffer_lengths[proc])
      continue;                             // no message to send to this proc
    CHKERR MPI_Isend(&(sbuffer[proc])[0],   // buffer to send
                     sbuffer_lengths[proc], // message length
                     MPIU_INT, proc,        // to proc
                     tag, m_field.get_comm(), s_waits + kk);
    kk++;
  }

  // Wait for received
  if (nrecvs)
    CHKERR MPI_Waitall(nrecvs, r_waits, &status[0]);

  // Wait for send messages
  if (nsends)
    CHKERR MPI_Waitall(nsends, s_waits, &status[0]);

  if (verb >= VERY_VERBOSE) {
    PetscSynchronizedPrintf(m_field.get_comm(), "Rank %d nb. before ents %u\n",
                            m_field.get_comm_rank(), ents.size());
  }

  // synchronise range
  for (int kk = 0; kk < nrecvs; kk++) {

    int len = olengths[kk];
    int *data_from_proc = rbuf[kk];

    for (int ee = 0; ee < len; ee += block_size) {

      EntityHandle ent;
      bcopy(&data_from_proc[ee], &ent, sizeof(EntityHandle));
      auto meit = ref_ents_ptr->get<Ent_mi_tag>().find(ent);
      if (meit == ref_ents_ptr->get<Ent_mi_tag>().end())
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "rank %d entity %lu not exist on database, local entity can "
                 "not be found for this owner",
                 m_field.get_comm_rank(), ent);

      if (verb >= VERY_VERBOSE)
        PetscSynchronizedPrintf(
            m_field.get_comm(), "received %ul (%ul) from %d at %d\n",
            (*meit)->getEnt(), ent, onodes[kk], m_field.get_comm_rank());

      ents.insert((*meit)->getEnt());
    }
  }

  if (verb >= VERBOSE)
    PetscSynchronizedPrintf(m_field.get_comm(), "Rank %d nb. after ents %u\n",
                            m_field.get_comm_rank(), ents.size());

  // Cleaning
  CHKERR PetscFree(s_waits);
  CHKERR PetscFree(rbuf[0]);
  CHKERR PetscFree(rbuf);
  CHKERR PetscFree(r_waits);
  CHKERR PetscFree(onodes);
  CHKERR PetscFree(olengths);

  if (verb >= VERBOSE)
    PetscSynchronizedFlush(m_field.get_comm(), PETSC_STDOUT);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CommInterface::synchroniseFieldEntities(const std::string name,
                                                       int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  EntityHandle idm = m_field.get_field_meshset(name);
  Range ents;
  CHKERR m_field.get_moab().get_entities_by_handle(idm, ents, false);
  CHKERR synchroniseEntities(ents, verb);
  CHKERR m_field.get_moab().add_entities(idm, ents);
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
  Tag th_gid;
  const int zero = 0;
  CHKERR m_field.get_moab().tag_get_handle(GLOBAL_ID_TAG_NAME, 1,
                                           MB_TYPE_INTEGER, th_gid,
                                           MB_TAG_DENSE | MB_TAG_CREAT, &zero);
  PetscLayout layout;
  CHKERR problem_ptr->getNumberOfElementsByNameAndPart(m_field.get_comm(),
                                                       fe_name, &layout);
  int gid, last_gid;
  CHKERR PetscLayoutGetRange(layout, &gid, &last_gid);
  CHKERR PetscLayoutDestroy(&layout);
  for (_IT_NUMEREDFE_BY_NAME_FOR_LOOP_(problem_ptr, fe_name, fe_it)) {
    EntityHandle ent = (*fe_it)->getEnt();
    ents.insert(ent);
    unsigned int part = (*fe_it)->getPart();
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
  CHKERR pcomm->exchange_tags(th_gid, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CommInterface::resolveSharedFiniteElements(
    const std::string &name, const std::string &fe_name, int verb) {
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
      Tag th_gid;
      const int zero = 0;
      CHKERR m_field.get_moab().tag_get_handle(
          "TMP_GLOBAL_ID_TAG_NAME", 1, MB_TYPE_INTEGER, th_gid,
          MB_TAG_SPARSE | MB_TAG_CREAT, &zero);
      return th_gid;
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
        PetscSynchronizedPrintf(m_field.get_comm(), "%s", ss.str().c_str());
        PetscSynchronizedFlush(m_field.get_comm(), PETSC_STDOUT);

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

} // namespace MoFEM