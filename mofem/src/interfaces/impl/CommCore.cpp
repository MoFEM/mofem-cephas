/** \file CommCore.cpp
 * \brief Core functions interface for managing communication
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

MoFEMErrorCode Core::synchronise_entities(Range &ents, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;

  // make a buffer
  std::vector<std::vector<EntityHandle>> sbuffer(sIze);

  Range::iterator eit = ents.begin();
  for (; eit != ents.end(); eit++) {

    RefEntity_multiIndex::iterator meit;
    meit = refinedEntities.get<Ent_mi_tag>().find(*eit);
    if (meit == refinedEntities.get<Ent_mi_tag>().end()) {
      continue;
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "rank %d entity %lu not exist on database, local entity can not "
               "be found for this owner",
               rAnk, *eit);
    }

    unsigned char pstatus = (*meit)->getPStatus();

    if (pstatus == 0)
      continue;

    if (verb >= NOISY) {
      std::ostringstream zz;
      zz << "pstatus " << std::bitset<8>(pstatus) << " ";
      PetscSynchronizedPrintf(cOmm, "%s", zz.str().c_str());
    }

    for (int proc = 0;
         proc < MAX_SHARING_PROCS && -1 != (*meit)->getSharingProcsPtr()[proc];
         proc++) {
      if ((*meit)->getSharingProcsPtr()[proc] == -1)
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
                "sharing processor not set");

      if ((*meit)->getSharingProcsPtr()[proc] == rAnk)
        continue;

      EntityHandle handle_on_sharing_proc =
          (*meit)->getSharingHandlersPtr()[proc];
      sbuffer[(*meit)->getSharingProcsPtr()[proc]].push_back(
          handle_on_sharing_proc);
      if (verb >= NOISY)
        PetscSynchronizedPrintf(cOmm, "send %lu (%lu) to %d at %d\n",
                                (*meit)->getRefEnt(), handle_on_sharing_proc,
                                (*meit)->getSharingProcsPtr()[proc], rAnk);

      if (!(pstatus & PSTATUS_MULTISHARED))
        break;
    }
  }

  int nsends = 0;                         // number of messages to send
  std::vector<int> sbuffer_lengths(sIze); // length of the message to proc
  const size_t block_size = sizeof(EntityHandle) / sizeof(int);
  for (int proc = 0; proc < sIze; proc++) {

    if (!sbuffer[proc].empty()) {

      sbuffer_lengths[proc] = sbuffer[proc].size() * block_size;
      nsends++;

    } else {

      sbuffer_lengths[proc] = 0;
    }
  }

  // Make sure it is a PETSc cOmm
  CHKERR PetscCommDuplicate(cOmm, &cOmm, NULL);

  std::vector<MPI_Status> status(sIze);

  // Computes the number of messages a node expects to receive
  int nrecvs; // number of messages received
  CHKERR PetscGatherNumberOfMessages(cOmm, NULL, &sbuffer_lengths[0], &nrecvs);

  // Computes info about messages that a MPI-node will receive, including
  // (from-id,length) pairs for each message.
  int *onodes;   // list of node-ids from which messages are expected
  int *olengths; // corresponding message lengths
  CHKERR PetscGatherMessageLengths(cOmm, nsends, nrecvs, &sbuffer_lengths[0],
                                   &onodes, &olengths);

  // Gets a unique new tag from a PETSc communicator. All processors that share
  // the communicator MUST call this routine EXACTLY the same number of times.
  // This tag should only be used with the current objects communicator; do NOT
  // use it with any other MPI communicator.
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
    if (!sbuffer_lengths[proc])
      continue;                             // no message to send to this proc
    CHKERR MPI_Isend(&(sbuffer[proc])[0],   // buffer to send
                     sbuffer_lengths[proc], // message length
                     MPIU_INT, proc,        // to proc
                     tag, cOmm, s_waits + kk);
    kk++;
  }

  // Wait for received
  if (nrecvs)
    CHKERR MPI_Waitall(nrecvs, r_waits, &status[0]);

  // Wait for send messages
  if (nsends)
    CHKERR MPI_Waitall(nsends, s_waits, &status[0]);

  if (verb >= VERY_VERBOSE) {
    PetscSynchronizedPrintf(cOmm, "Rank %d nb. before ents %u\n", rAnk,
                            ents.size());
  }

  // synchronise range
  for (int kk = 0; kk < nrecvs; kk++) {

    int len = olengths[kk];
    int *data_from_proc = rbuf[kk];

    for (int ee = 0; ee < len; ee += block_size) {

      EntityHandle ent;
      bcopy(&data_from_proc[ee], &ent, sizeof(EntityHandle));
      RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator meit;
      meit = refinedEntities.get<Ent_mi_tag>().find(ent);
      if (meit == refinedEntities.get<Ent_mi_tag>().end())
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "rank %d entity %lu not exist on database, local entity can "
                 "not be found for this owner",
                 rAnk, ent);

      if (verb >= VERY_VERBOSE)
        PetscSynchronizedPrintf(cOmm, "received %ul (%ul) from %d at %d\n",
                                (*meit)->getRefEnt(), ent, onodes[kk], rAnk);

      ents.insert((*meit)->getRefEnt());
    }
  }

  if (verb >= VERBOSE)
    PetscSynchronizedPrintf(cOmm, "Rank %d nb. after ents %u\n", rAnk,
                            ents.size());

  // Cleaning
  CHKERR PetscFree(s_waits);
  CHKERR PetscFree(rbuf[0]);
  CHKERR PetscFree(rbuf);
  CHKERR PetscFree(r_waits);
  CHKERR PetscFree(onodes);
  CHKERR PetscFree(olengths);

  if (verb >= VERBOSE)
    PetscSynchronizedFlush(cOmm, PETSC_STDOUT);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::synchronise_field_entities(const BitFieldId id, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  EntityHandle idm;
  idm = get_field_meshset(id);
  Range ents;
  CHKERR get_moab().get_entities_by_handle(idm, ents, false);
  CHKERR synchronise_entities(ents, verb);
  CHKERR get_moab().add_entities(idm, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::synchronise_field_entities(const std::string &name,
                                                int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR synchronise_field_entities(getBitFieldId(name), verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::resolve_shared_finite_elements(const Problem *problem_ptr,
                                                    const std::string &fe_name,
                                                    int verb) {
  MoFEMFunctionBegin;
  ParallelComm *pcomm =
      ParallelComm::get_pcomm(&get_moab(), basicEntityDataPtr->pcommID);
  std::vector<int> shprocs(MAX_SHARING_PROCS, 0);
  std::vector<EntityHandle> shhandles(MAX_SHARING_PROCS, 0);
  Range ents;
  Tag th_gid;
  const int zero = 0;
  CHKERR get_moab().tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER,
                                   th_gid, MB_TAG_DENSE | MB_TAG_CREAT, &zero);
  PetscLayout layout;
  CHKERR problem_ptr->getNumberOfElementsByNameAndPart(get_comm(), fe_name,
                                                       &layout);
  int gid, last_gid;
  CHKERR PetscLayoutGetRange(layout, &gid, &last_gid);
  CHKERR PetscLayoutDestroy(&layout);
  for (_IT_NUMEREDFE_BY_NAME_FOR_LOOP_(problem_ptr, fe_name, fe_it)) {
    EntityHandle ent = (*fe_it)->getEnt();
    ents.insert(ent);
    unsigned int part = (*fe_it)->getPart();
    CHKERR get_moab().tag_set_data(pcomm->part_tag(), &ent, 1, &part);
    if (part == pcomm->rank()) {
      CHKERR get_moab().tag_set_data(th_gid, &ent, 1, &gid);
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
        CHKERR get_moab().tag_set_data(pcomm->sharedp_tag(), &ent, 1,
                                       &shprocs[0]);
        CHKERR get_moab().tag_set_data(pcomm->sharedh_tag(), &ent, 1,
                                       &shhandles[0]);
      }

      if (pstatus & PSTATUS_MULTISHARED) {
        CHKERR get_moab().tag_set_data(pcomm->sharedps_tag(), &ent, 1,
                                       &shprocs[0]);
        CHKERR get_moab().tag_set_data(pcomm->sharedhs_tag(), &ent, 1,
                                       &shhandles[0]);
      }
      CHKERR get_moab().tag_set_data(pcomm->pstatus_tag(), &ent, 1, &pstatus);
    }
  }
  CHKERR pcomm->exchange_tags(th_gid, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::resolve_shared_finite_elements(const std::string &name,
                                                    const std::string &fe_name,
                                                    int verb) {
  MoFEMFunctionBegin;
  // find p_miit
  auto p_miit = pRoblems.get<Problem_mi_tag>().find(name);
  if (p_miit == pRoblems.get<Problem_mi_tag>().end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
             "problem with name < %s > not defined (top tip check spelling)",
             name.c_str());
  CHKERR resolve_shared_finite_elements(&*p_miit, fe_name, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::make_entities_multishared(const EntityHandle *entities,
                                               const int num_entities,
                                               const int owner_proc, int verb) {
  MoFEMFunctionBegin;

  if (sIze > 1) {

    ParallelComm *pcomm =
        ParallelComm::get_pcomm(&get_moab(), basicEntityDataPtr->pcommID);
    const EntityHandle *ent = entities;
    const EntityHandle *const end = entities + num_entities;
    std::vector<EntityHandle> all_ents_vec(ent, end);

    auto print_vec = [&](auto name, auto &vec) {
      std::ostringstream ss;
      ss << "Proc " << rAnk << " ";
      ss << "vector " << name << " [ ";
      for (auto v : vec)
        ss << v << " ";
      ss << "]" << endl;
      PetscSynchronizedPrintf(cOmm, "%s", ss.str().c_str());
      PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
    };

    if (verb >= NOISY)
      print_vec("all_ents_vec", all_ents_vec);

    const size_t block_size = sizeof(EntityHandle) / sizeof(int);

    int tag;
    CHKERR PetscCommGetNewTag(cOmm, &tag);

    std::vector<EntityHandle> recv_ents_vec[sIze];
    std::vector<MPI_Request> r_waits(sIze);
    for (int proc = 0, kk = 0; proc < sIze; ++proc) {
      if (proc != rAnk) {
        recv_ents_vec[proc].resize(all_ents_vec.size());
        CHKERR MPI_Irecv(&*recv_ents_vec[proc].begin(), // buffer to receive
                         recv_ents_vec[proc].size() *
                             block_size, // message length
                         MPIU_INT, proc, // to proc
                         tag, cOmm, &r_waits[kk]);
        ++kk;
      }
    }

    std::vector<MPI_Request> s_waits(sIze); // status of send messages
    for (int proc = 0, kk = 0; proc < sIze; ++proc) {
      if (proc != rAnk) {
        CHKERR MPI_Isend(&*all_ents_vec.begin(),           // buffer to send
                         all_ents_vec.size() * block_size, // message length
                         MPIU_INT, proc,                   // to proc
                         tag, cOmm, &s_waits[kk]);
        ++kk;
      }
    }

    std::vector<MPI_Status> status(sIze);

    // Wait for received
    CHKERR MPI_Waitall(sIze - 1, &*r_waits.begin(), &*status.begin());

    // Wait for send messages
    CHKERR MPI_Waitall(sIze - 1, &*s_waits.begin(), &*status.begin());

    if (verb >= NOISY) {
      for (int proc = 0, kk = 0; proc < sIze; ++proc) {
        if (proc != rAnk) {
          print_vec("recv_ents_vec received from " +
                        boost::lexical_cast<std::string>(proc),
                    recv_ents_vec[proc]);
        }
      }
    }

    unsigned char pstatus = 0;

    if (rAnk != owner_proc)
      pstatus = PSTATUS_NOT_OWNED;

    pstatus |= PSTATUS_SHARED;
    if(sIze > 2)
      pstatus |= PSTATUS_MULTISHARED;

    if (pcomm->size() != sIze)
      SETERRQ(PETSC_COMM_WORLD, MOFEM_NOT_IMPLEMENTED,
              "Implementation that size of PComm and MoFEM are diffrent not "
              "implemented");

    std::vector<int> shprocs(MAX_SHARING_PROCS, 0);
    std::vector<EntityHandle> shhandles(MAX_SHARING_PROCS, 0);

    for (int e = 0; e != all_ents_vec.size(); ++e) {

      size_t rrr = 0;
      if (rAnk != owner_proc) {
        shhandles[rrr] = (recv_ents_vec[owner_proc])[e];
        shprocs[rrr] = owner_proc;
        ++rrr;
      }

      for (size_t rr = 0; rr < sIze; ++rr) {
        if (rr != owner_proc && rr != rAnk) {
          shhandles[rrr] = (recv_ents_vec[rr])[e];
          shprocs[rrr] = rr;
          ++rrr;
        }
      }

      for (; rrr != MAX_SHARING_PROCS; ++rrr)
        shprocs[rrr] = -1;

      CHKERR get_moab().tag_set_data(pcomm->sharedp_tag(), &all_ents_vec[e], 1,
                                     &shprocs[0]);
      CHKERR get_moab().tag_set_data(pcomm->sharedh_tag(), &all_ents_vec[e], 1,
                                     &shhandles[0]);
      if (pstatus & PSTATUS_MULTISHARED) {
        CHKERR get_moab().tag_set_data(pcomm->sharedps_tag(), &all_ents_vec[e],
                                       1, &shprocs[0]);
        CHKERR get_moab().tag_set_data(pcomm->sharedhs_tag(), &all_ents_vec[e],
                                       1, &shhandles[0]);
      }

      CHKERR get_moab().tag_set_data(pcomm->pstatus_tag(), &all_ents_vec[e], 1,
                                     &pstatus);

      auto print_owner = [&]() {
        MoFEMFunctionBegin;
        int moab_owner_proc;
        EntityHandle moab_owner_handle;
        CHKERR pcomm->get_owner_handle(all_ents_vec[e], moab_owner_proc,
                                       moab_owner_handle);

        std::ostringstream ss;

        ss << "Rank " << rAnk << " ";
        if (!(pstatus & PSTATUS_NOT_OWNED))
          ss << "OWNER ";
        if (pstatus & PSTATUS_SHARED)
          ss << "PSTATUS_SHARED ";
        else if (pstatus & PSTATUS_MULTISHARED)
          ss << "PSTATUS_MULTISHARED ";

        ss << "owner " << moab_owner_proc << " (" << owner_proc << ") ";

        ss << std::endl;
        PetscSynchronizedPrintf(cOmm, "%s", ss.str().c_str());
        PetscSynchronizedFlush(cOmm, PETSC_STDOUT);

        MoFEMFunctionReturn(0);
      };

      if (verb >= NOISY)
        CHKERR print_owner();

      auto check_owner = [&]() {
        MoFEMFunctionBegin;
        int moab_owner_proc;
        EntityHandle moab_owner_handle;
        CHKERR pcomm->get_owner_handle(all_ents_vec[e], moab_owner_proc,
                                       moab_owner_handle);
        if (owner_proc != moab_owner_proc)
          SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "%d != %d",
                   owner_proc, moab_owner_proc);
        if (rAnk != owner_proc)
          if (recv_ents_vec[owner_proc][e] != moab_owner_handle)
            SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "%lu != %lu",
                     recv_ents_vec[owner_proc][e], moab_owner_handle);
        MoFEMFunctionReturn(0);
      };

      CHKERR check_owner();
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::make_entities_multishared(Range &entities,
                                               const int owner_proc, int verb) {
  MoFEMFunctionBegin;
  if (sIze > 1) {
    const int num_ents = entities.size();
    std::vector<EntityHandle> vec_ents(num_ents);
    std::copy(entities.begin(), entities.end(), vec_ents.begin());
    CHKERR make_entities_multishared(&*vec_ents.begin(), num_ents, owner_proc,
                                     verb);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::make_field_entities_multishared(const std::string field_name,
                                      const int owner_proc, int verb) {
  MoFEMFunctionBegin;
  if (sIze > 1) {
    EntityHandle field_meshset = get_field_meshset(field_name);
    std::vector<EntityHandle> field_ents;
    CHKERR get_moab().get_entities_by_handle(field_meshset, field_ents, true);
    CHKERR make_entities_multishared(&*field_ents.begin(), field_ents.size(),
                                     owner_proc, verb);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::exchange_field_data(const std::string field_name,
                                         int verb) {
  MoFEMFunctionBegin;
  if (sIze > 1) {

    Range exchange_ents_data_verts, exchange_ents_data;

    for (auto it = entsFields.get<FieldName_mi_tag>().lower_bound(field_name);
         it != entsFields.get<FieldName_mi_tag>().upper_bound(field_name); ++it)
      if (

          ((*it)->getPStatus()) &&

          (*it)->getNbDofsOnEnt()

      ) {
        if ((*it)->getEntType() == MBVERTEX)
          exchange_ents_data_verts.insert((*it)->getRefEnt());
        else 
          exchange_ents_data.insert((*it)->getRefEnt());  
      }

    auto field_ptr = get_field_structure(field_name);
    ParallelComm *pcomm =
        ParallelComm::get_pcomm(&get_moab(), basicEntityDataPtr->pcommID);

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
