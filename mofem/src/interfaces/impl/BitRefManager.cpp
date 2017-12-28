/** \file BitRefManager.cpp
 * \brief Managing BitRefLevels
 * \mofem_bit_ref
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

MoFEMErrorCode BitRefManager::query_interface(const MOFEMuuid &uuid,
                                              UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMBitRefManager) {
    *iface = const_cast<BitRefManager *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

BitRefManager::BitRefManager(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)), dEbug(false) {}
BitRefManager::~BitRefManager() {}

MoFEMErrorCode BitRefManager::setBitRefLevel(const Range &ents,
                                             const BitRefLevel &bit,
                                             const bool only_tets,
                                             int verb) const {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ents_ptr;
  const RefElement_multiIndex *ref_fe_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.get_ref_ents(&ref_ents_ptr);
  CHKERRG(ierr);
  ierr = m_field.get_ref_finite_elements(&ref_fe_ptr);
  CHKERRG(ierr);
  Range seeded_ents;
  try {

    if (verb > 1) {
      PetscSynchronizedPrintf(m_field.get_comm(), "nb. entities for seed %d\n",
                              ents.size());
    }

    /// tool class with methods used more than twp times
    struct Tool {

      const BitRefLevel &bIt;                          ///< bit to set
      const RefEntity_multiIndex *refEntsPtr;          ///< access to database
      boost::shared_ptr<BasicEntityData> &baseEntData; ///< base entity data

      /// constrictor
      Tool(MoFEM::Interface &m_field, const BitRefLevel &bit,
           const RefEntity_multiIndex *ref_ents_ptr)
          : bIt(bit), refEntsPtr(ref_ents_ptr),
            baseEntData(m_field.get_basic_entity_data_ptr()) {}

      /// find entities and change entity bit if in database
      MoFEMErrorCode
      findEntsToAdd(EntityHandle f, EntityHandle s,
                    std::vector<EntityHandle> &seed_ents_vec,
                    std::vector<boost::shared_ptr<RefEntity> >
                        *shared_ref_ents_vec_for_fe = NULL) const {
        MoFEMFunctionBeginHot;
        RefEntity_multiIndex::iterator rit, hi_rit;
        // get lower bound of multi-index
        rit = refEntsPtr->lower_bound(f);
        if (rit == refEntsPtr->end()) {
          // all enties in range are added to database
          seed_ents_vec.reserve(s - f + 1);
          for (; f <= s; f++) {
            seed_ents_vec.push_back(f);
          }
        } else {
          // some entities from range are in database
          hi_rit = refEntsPtr->upper_bound(s);
          for (; f <= s; f++) {
            if (f == rit->get()->getRefEnt()) {
              // entity is in database, change bit level only
              bool success = const_cast<RefEntity_multiIndex *>(refEntsPtr)
                                 ->modify(rit, RefEntity_change_add_bit(bIt));
              if (!success) {
                SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                        "modification unsuccessful");
              }
              if (shared_ref_ents_vec_for_fe != NULL) {
                shared_ref_ents_vec_for_fe->push_back(*rit);
              }
              rit++; // move to next one
              if (rit == hi_rit) {
                // break loop, rest of the entities in range are not in database
                break;
              }
            } else {
              // this entities added to database
              seed_ents_vec.push_back(f);
            }
          }
          // add rest entities to vector of entities going to be added to
          // database
          for (; f <= s; f++) {
            seed_ents_vec.push_back(f);
          }
        }
        MoFEMFunctionReturnHot(0);
      }

      /// add entities to database
      MoFEMErrorCode
      addToDatabase(std::vector<EntityHandle> &seed_ents_vec,
                    std::vector<boost::shared_ptr<RefEntity> >
                        *shared_ref_ents_vec_for_fe = NULL) const {
        MoFEMFunctionBeginHot;
        // add entities to database
        boost::shared_ptr<std::vector<RefEntity> > ref_ents_vec =
            boost::make_shared<std::vector<RefEntity> >();
        ref_ents_vec->reserve(seed_ents_vec.size());
        // create ref entity instances
        for (std::vector<EntityHandle>::const_iterator vit =
                 seed_ents_vec.begin();
             vit != seed_ents_vec.end(); vit++) {
          ref_ents_vec->push_back(RefEntity(baseEntData, *vit));
          RefEntity_change_add_bit(bIt).operator()(ref_ents_vec->back());
        }
        std::vector<boost::shared_ptr<RefEntity> > shared_ref_ents_vec;
        shared_ref_ents_vec.reserve(ref_ents_vec->size());
        // create aliased shared pointers to ref entity instances
        for (std::vector<RefEntity>::iterator vit = ref_ents_vec->begin();
             vit != ref_ents_vec->end(); vit++) {
          shared_ref_ents_vec.push_back(
              boost::shared_ptr<RefEntity>(ref_ents_vec, &*vit));
        }
        if (shared_ref_ents_vec_for_fe) {
          shared_ref_ents_vec_for_fe->insert(shared_ref_ents_vec_for_fe->end(),
                                             shared_ref_ents_vec.begin(),
                                             shared_ref_ents_vec.end());
        }
        // add shared pointers to database
        const_cast<RefEntity_multiIndex *>(refEntsPtr)
            ->insert(shared_ref_ents_vec.begin(), shared_ref_ents_vec.end());
        MoFEMFunctionReturnHot(0);
      }
    };

    for (Range::const_pair_iterator pit = ents.pair_begin();
         pit != ents.pair_end(); pit++) {
      // get first and last element of range
      EntityHandle f = pit->first;
      EntityHandle s = pit->second;
      std::vector<EntityHandle>
          seed_ents_vec; // entities seeded not in database
      std::vector<boost::shared_ptr<RefEntity> > shared_ref_ents_vec_for_fe;
      // find ents to add
      ierr =
          Tool(m_field, bit, ref_ents_ptr)
              .findEntsToAdd(f, s, seed_ents_vec, &shared_ref_ents_vec_for_fe);
      CHKERRG(ierr);
      // add elements
      if (!seed_ents_vec.empty()) {
        ierr = Tool(m_field, bit, ref_ents_ptr)
                   .addToDatabase(seed_ents_vec, &shared_ref_ents_vec_for_fe);
        CHKERRG(ierr);
      }

      // create finite elements
      for (std::vector<boost::shared_ptr<RefEntity> >::iterator vit =
               shared_ref_ents_vec_for_fe.begin();
           vit != shared_ref_ents_vec_for_fe.end(); vit++) {
        std::pair<RefElement_multiIndex::iterator, bool> p_fe;
        switch ((*vit)->getEntType()) {
        case MBVERTEX:
          p_fe =
              const_cast<RefElement_multiIndex *>(ref_fe_ptr)
                  ->insert(ptrWrapperRefElement(boost::shared_ptr<RefElement>(
                      new RefElement_VERTEX(*vit))));
          seeded_ents.insert((*vit)->getRefEnt());
          break;
        case MBEDGE:
          p_fe =
              const_cast<RefElement_multiIndex *>(ref_fe_ptr)
                  ->insert(ptrWrapperRefElement(boost::shared_ptr<RefElement>(
                      new RefElement_EDGE(*vit))));
          seeded_ents.insert((*vit)->getRefEnt());
          break;
        case MBTRI:
          p_fe =
              const_cast<RefElement_multiIndex *>(ref_fe_ptr)
                  ->insert(ptrWrapperRefElement(
                      boost::shared_ptr<RefElement>(new RefElement_TRI(*vit))));
          seeded_ents.insert((*vit)->getRefEnt());
          break;
        case MBTET:
          p_fe =
              const_cast<RefElement_multiIndex *>(ref_fe_ptr)
                  ->insert(ptrWrapperRefElement(
                      boost::shared_ptr<RefElement>(new RefElement_TET(*vit))));
          seeded_ents.insert((*vit)->getRefEnt());
          break;
        case MBPRISM:
          p_fe =
              const_cast<RefElement_multiIndex *>(ref_fe_ptr)
                  ->insert(ptrWrapperRefElement(boost::shared_ptr<RefElement>(
                      new RefElement_PRISM(*vit))));
          if (!only_tets) {
            seeded_ents.insert((*vit)->getRefEnt());
          }
          break;
        case MBENTITYSET:
          p_fe =
              const_cast<RefElement_multiIndex *>(ref_fe_ptr)
                  ->insert(ptrWrapperRefElement(boost::shared_ptr<RefElement>(
                      new RefElement_MESHSET(*vit))));
          break;
        default:
          SETERRQ(m_field.get_comm(), MOFEM_NOT_IMPLEMENTED, "not implemented");
        }
      }
    }

    if (!seeded_ents.empty()) {
      int dim = m_field.get_moab().dimension_from_handle(seeded_ents[0]);
      for (int dd = 0; dd < dim; dd++) {
        Range ents;
        rval = m_field.get_moab().get_adjacencies(seeded_ents, dd, true, ents,
                                                  moab::Interface::UNION);
        CHKERRQ_MOAB(rval);
        if (dd == 2 && only_tets) {
          // currently only works with triangles
          ents = ents.subset_by_type(MBTRI);
        }
        std::vector<EntityHandle>
            seed_ents_vec; // entities seeded not in database
        for (Range::pair_iterator pit = ents.pair_begin();
             pit != ents.pair_end(); pit++) {
          seed_ents_vec.clear();
          // get first and last element of range
          EntityHandle f = pit->first;
          EntityHandle s = pit->second;
          ierr = Tool(m_field, bit, ref_ents_ptr)
                     .findEntsToAdd(f, s, seed_ents_vec);
          CHKERRG(ierr);
          if (!seed_ents_vec.empty()) {
            ierr =
                Tool(m_field, bit, ref_ents_ptr).addToDatabase(seed_ents_vec);
            CHKERRG(ierr);
          }
        }
      }
    }

  } catch (MoFEMException const &e) {
    SETERRQ(m_field.get_comm(), e.errorCode, e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::addToDatabaseBitRefLevelByType(
    const EntityType type, const BitRefLevel &bit, const bool only_tets,
    int verb) const {
  MoFEMFunctionBeginHot;
  Range ents;
  ierr = getEntitiesByTypeAndRefLevel(bit, BitRefLevel().set(), type, ents);
  CHKERRG(ierr);
  // Add bit ref level to database
  ierr = setBitRefLevel(ents, bit);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::addToDatabaseBitRefLevelByDim(
    const int dim, const BitRefLevel &bit, const bool only_tets,
    int verb) const {
  MoFEMFunctionBeginHot;
  Range ents;
  ierr = getEntitiesByDimAndRefLevel(bit, BitRefLevel().set(), dim, ents);
  CHKERRG(ierr);
  // Add bit ref level to database
  ierr = setBitRefLevel(ents, bit);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::setBitLevelToMeshset(const EntityHandle meshset,
                                                   const BitRefLevel &bit,
                                                   int verb) const {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ents_ptr;
  const RefElement_multiIndex *ref_fe_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.get_ref_ents(&ref_ents_ptr);
  CHKERRG(ierr);
  ierr = m_field.get_ref_finite_elements(&ref_fe_ptr);
  CHKERRG(ierr);
  std::pair<RefEntity_multiIndex::iterator, bool> p_ent =
      const_cast<RefEntity_multiIndex *>(ref_ents_ptr)
          ->insert(boost::shared_ptr<RefEntity>(
              new RefEntity(m_field.get_basic_entity_data_ptr(), meshset)));
  const_cast<RefEntity_multiIndex *>(ref_ents_ptr)
      ->modify(p_ent.first, RefEntity_change_add_bit(bit));
  ptrWrapperRefElement pack_fe(
      boost::shared_ptr<RefElement>(new RefElement_MESHSET(*p_ent.first)));
  std::pair<RefElement_multiIndex::iterator, bool> p_fe =
      const_cast<RefElement_multiIndex *>(ref_fe_ptr)->insert(pack_fe);
  if (verb > 0) {
    std::ostringstream ss;
    ss << "add meshset as ref_ent " << *(p_fe.first->getRefElement())
       << std::endl;
    PetscPrintf(m_field.get_comm(), ss.str().c_str());
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::setBitRefLevelByDim(const EntityHandle meshset,
                                                  const int dim,
                                                  const BitRefLevel &bit,
                                                  int verb) const {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBeginHot;
  Range ents;
  rval =
      m_field.get_moab().get_entities_by_dimension(meshset, dim, ents, false);
  CHKERRQ_MOAB(rval);
  ierr = setBitRefLevel(ents, bit, false, verb);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::setBitRefLevelByType(const EntityHandle meshset,
                                                   const EntityType type,
                                                   const BitRefLevel &bit,
                                                   int verb) const {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBeginHot;
  Range ents;
  rval = m_field.get_moab().get_entities_by_type(meshset, type, ents, false);
  CHKERRQ_MOAB(rval);
  ierr = setBitRefLevel(ents, bit, false, verb);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::addBitRefLevel(const Range &ents,
                                             const BitRefLevel &bit,
                                             int verb) const {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ent_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.get_ref_ents(&ref_ent_ptr);
  for (Range::const_pair_iterator pit = ents.const_pair_begin();
       pit != ents.const_pair_end(); pit++) {
    EntityHandle first = pit->first;
    EntityHandle second = pit->second;
    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator dit;
    dit = ref_ent_ptr->get<Ent_mi_tag>().lower_bound(first);
    if (dit == ref_ent_ptr->get<Ent_mi_tag>().end())
      continue;
    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator hi_dit;
    hi_dit = ref_ent_ptr->get<Ent_mi_tag>().upper_bound(second);
    for (; dit != hi_dit; dit++) {
      bool success = const_cast<RefEntity_multiIndex *>(ref_ent_ptr)
                         ->modify(ref_ent_ptr->project<0>(dit),
                                  RefEntity_change_add_bit(bit));
      if (!success) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                "operation unsuccessful");
      };
      if (verb >= VERY_NOISY) {
        cerr << **dit << endl;
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::addBitRefLevelByDim(const EntityHandle meshset,
                                                  const int dim,
                                                  const BitRefLevel &bit,
                                                  int verb) const {
  MoFEM::Interface& m_field = cOre;
  moab::Interface& moab = m_field.get_moab();
  Range ents, adj;
  MoFEMFunctionBeginHot;
  rval = moab.get_entities_by_dimension(meshset, dim, ents, true);
  CHKERRQ_MOAB(rval);
  for (int dd = dim - 1; dd >= 0; dd--) {
    rval = moab.get_adjacencies(ents, dd, false, adj, moab::Interface::UNION);
    CHKERRQ_MOAB(rval);
  }
  ents.merge(adj);
  if (verb == VERY_NOISY) {
    cerr << ents << endl;
  }
  ierr = addBitRefLevel(ents, bit, verb);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::setNthBitRefLevel(const Range &ents, const int n,
                                                const bool b, int verb) const {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ent_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.get_ref_ents(&ref_ent_ptr);
  for (Range::const_pair_iterator pit = ents.const_pair_begin();
       pit != ents.const_pair_end(); pit++) {
    EntityHandle first = pit->first;
    EntityHandle second = pit->second;
    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator dit;
    dit = ref_ent_ptr->get<Ent_mi_tag>().lower_bound(first);
    if (dit == ref_ent_ptr->get<Ent_mi_tag>().end())
      continue;
    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator hi_dit;
    hi_dit = ref_ent_ptr->get<Ent_mi_tag>().upper_bound(second);
    for (; dit != hi_dit; dit++) {
      bool success = const_cast<RefEntity_multiIndex *>(ref_ent_ptr)
                         ->modify(ref_ent_ptr->project<0>(dit),
                                  RefEntity_change_set_nth_bit(n, b));
      if (!success) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                "operation unsuccessful");
      };
      if (verb > 0) {
        cerr << **dit << endl;
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::setNthBitRefLevel(const int n, const bool b,
                                                int verb) const {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ent_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.get_ref_ents(&ref_ent_ptr);
  RefEntity_multiIndex::iterator dit, hi_dit;
  dit = ref_ent_ptr->begin();
  hi_dit = ref_ent_ptr->end();
  for (; dit != hi_dit; dit++) {
    bool success = const_cast<RefEntity_multiIndex *>(ref_ent_ptr)
                       ->modify(dit, RefEntity_change_set_nth_bit(n, b));
    if (!success) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
              "operation unsuccessful");
    };
    if (verb > 0) {
      cerr << **dit << endl;
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::shiftLeftBitRef(const int shift,
                                              const BitRefLevel mask,
                                              int verb) const {
  MoFEMFunctionBeginHot;
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::shiftRightBitRef(const int shift,
                                               const BitRefLevel mask,
                                               int verb) const {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ent_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_ref_ents(&ref_ent_ptr);
  for (int ii = 0; ii < shift; ii++) {
    // delete bits on the right which are shifted to zero
    BitRefLevel delete_bits = BitRefLevel().set(0) & mask;
    if (delete_bits.any()) {
      CHKERR m_field.delete_ents_by_bit_ref(delete_bits, delete_bits, true,
                                            verb);
    }
    for (RefEntity_multiIndex::iterator ent_it = ref_ent_ptr->begin();
         ent_it != ref_ent_ptr->end(); ent_it++) {
      if (verb > NOISY) {
        std::cerr << (*ent_it)->getBitRefLevel() << " : ";
      }
      bool success =
          const_cast<RefEntity_multiIndex *>(ref_ent_ptr)
              ->modify(ent_it, RefEntity_change_right_shift(1, mask));
      if (!success)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "inconsistency in data");
      if (verb >= VERY_NOISY) {
        std::cerr << (*ent_it)->getBitRefLevel() << std::endl;
      }
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::writeBitLevelByType(
    const BitRefLevel &bit, const BitRefLevel &mask, const EntityType type,
    const char *file_name, const char *file_type, const char *options) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBeginHot;
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET, meshset);
  CHKERRQ_MOAB(rval);
  ierr = getEntitiesByTypeAndRefLevel(bit, mask, type, meshset);
  CHKERRG(ierr);
  rval = moab.write_file(file_name, file_type, options, &meshset, 1);
  CHKERRQ_MOAB(rval);
  rval = moab.delete_entities(&meshset, 1);
  CHKERRQ_MOAB(rval);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::writeEntitiesNotInDatabase(
    const char *file_name, const char *file_type, const char *options) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBeginHot;
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET, meshset);
  CHKERRG(rval);
  Range ents;
  ierr = getAllEntitiesNotInDatabase(ents);
  CHKERRG(ierr);
  rval = moab.add_entities(meshset, ents);
  CHKERRG(rval);
  rval = moab.write_file(file_name, file_type, options, &meshset, 1);
  CHKERRG(rval);
  rval = moab.delete_entities(&meshset, 1);
  CHKERRG(rval);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::getEntitiesByTypeAndRefLevel(
    const BitRefLevel &bit, const BitRefLevel &mask, const EntityType type,
    const EntityHandle meshset, int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBeginHot;
  Range ents;
  ierr = getEntitiesByTypeAndRefLevel(bit, mask, type, ents, verb);
  CHKERRG(ierr);
  rval = moab.add_entities(meshset, ents);
  CHKERRG(rval);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::filterEntitiesByRefLevel(const BitRefLevel &bit,
                                                       const BitRefLevel &mask,
                                                       Range &ents,
                                                       int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBeginHot;

  const BitRefLevel *tag_bit;

  Range swap_ents;
  for (Range::pair_iterator p_eit = ents.pair_begin(); p_eit != ents.pair_end();
       ++p_eit) {
    EntityHandle f = p_eit->first;
    EntityHandle s = p_eit->second;

    EntityHandle ff = 0;
    EntityHandle ss = 0;

    for (; f <= s; ++f) {
      // Get entity bit ref level
      rval = moab.tag_get_by_ptr(cOre.get_th_RefBitLevel(), &f, 1,
                                 (const void **)(&tag_bit));
      if (PetscLikely(rval == MB_SUCCESS)) {
        if ((mask.any() && tag_bit->none()) ||
            ((*tag_bit) & mask) != (*tag_bit) || (((*tag_bit) & bit).none())) {
          // Entity not on BitRefLevel.
          if (ff != 0) {
            // Add sub-range
            swap_ents.insert(ff, ss);
            ff = ss = 0;
          }
        } else {
          if (ff == 0) {
            // Start new sub-range
            ff = ss = f;
          } else {
            // Add entity to sub-range
            ss = f;
          }
        }
      } else {
        if (ff != 0) {
          // Add sub-range
          swap_ents.insert(ff, ss);
        }
        // Start new sub-range
        ff = ss = 0;
      }
    }

    if (ff != 0) {
      swap_ents.insert(ff, ss);
    }

  }

  ents.swap(swap_ents);

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::getEntitiesByTypeAndRefLevel(
    const BitRefLevel &bit, const BitRefLevel &mask, const EntityType type,
    Range &ents, int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBeginHot;
  ierr = moab.get_entities_by_type(0, type, ents, false);
  CHKERRG(ierr);
  ierr = filterEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::getEntitiesByDimAndRefLevel(const BitRefLevel &bit,
                                           const BitRefLevel &mask,
                                           const int dim,
                                           const EntityHandle meshset,
                                           int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBeginHot;
  Range ents;
  ierr = getEntitiesByDimAndRefLevel(bit, mask, dim, ents, verb);
  CHKERRG(ierr);
  rval = moab.add_entities(meshset, ents);
  CHKERRQ_MOAB(rval);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::getEntitiesByDimAndRefLevel(
    const BitRefLevel &bit, const BitRefLevel &mask, const int dim,
    Range &ents, int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBeginHot;
  ierr = moab.get_entities_by_dimension(0, dim, ents, false);
  ierr = filterEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::getEntitiesByRefLevel(const BitRefLevel &bit,
                                                    const BitRefLevel &mask,
                                                    const EntityHandle meshset,
                                                    const int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  Range ents;
  CHKERR getEntitiesByRefLevel(bit, mask, ents,verb);
  CHKERR moab.add_entities(meshset, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::getEntitiesByRefLevel(const BitRefLevel &bit,
                                                    const BitRefLevel &mask,
                                                    Range &ents,
                                                    const int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  Range meshset_ents;
  CHKERR moab.get_entities_by_type(0, MBENTITYSET, meshset_ents, false);
  CHKERR moab.get_entities_by_handle(0, ents, false);
  ents.merge(meshset_ents);
  CHKERR filterEntitiesByRefLevel(bit, mask, ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::getEntitiesByParentType(const BitRefLevel &bit,
                                                      const BitRefLevel &mask,
                                                      const EntityType type,
                                                      Range &ents) const {
  MoFEM::Interface &m_field = cOre;
  // moab::Interface &moab = m_field.get_moab();
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.get_ref_ents(&ref_ents_ptr);
  CHKERRG(ierr);
  typedef RefEntity_multiIndex::index<ParentEntType_mi_tag>::type
      RefEntsByParentType;
  const RefEntsByParentType &ref_ents =
      ref_ents_ptr->get<ParentEntType_mi_tag>();
  RefEntsByParentType::iterator it, hi_it;
  it = ref_ents.lower_bound(type);
  hi_it = ref_ents.upper_bound(type);
  for (; it != hi_it; it++) {
    const BitRefLevel &ent_bit = it->get()->getBitRefLevel();
    if ((ent_bit & mask) == ent_bit && (ent_bit & bit).any()) {
      ents.insert(it->get()->getRefEnt());
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
BitRefManager::getAllEntitiesNotInDatabase(Range &ents) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  rval = moab.get_entities_by_handle(0,ents,false); CHKERRQ_MOAB(rval);
  ents = subtract(ents,ents.subset_by_type(MBENTITYSET));
  ierr = filterEntitiesNotInDatabase(ents); CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::filterEntitiesNotInDatabase(Range &ents) const {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.get_ref_ents(&ref_ents_ptr);
  CHKERRG(ierr);
  Range::iterator eit = ents.begin();
  for (; eit != ents.end(); ) {
    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator rit;
    rit = ref_ents_ptr->get<Ent_mi_tag>().find(*eit);
    if (rit != ref_ents_ptr->get<Ent_mi_tag>().end()) {
      eit = ents.erase(eit);
    } else {
      eit++;
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
BitRefManager::getAdjacenciesEquality(const EntityHandle from_entity,
                                      const int to_dimension,
                                      Range &adj_entities) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBeginHot;
  RefEntity from_ref_entity(m_field.get_basic_entity_data_ptr(), from_entity);
  // std::cerr << "from:\n";
  // std::cerr << from_ref_entity << std::endl;
  rval =
      moab.get_adjacencies(&from_entity, 1, to_dimension, false, adj_entities);
  CHKERRQ_MOAB(rval);
  std::vector<BitRefLevel> bit_levels(adj_entities.size());
  rval = moab.tag_get_data(cOre.get_th_RefBitLevel(), adj_entities,
                           &*bit_levels.begin());
  CHKERRQ_MOAB(rval);
  std::vector<BitRefLevel>::iterator b_it = bit_levels.begin();
  Range::iterator eit = adj_entities.begin();
  // std::cerr << "to:\n";
  for (; eit != adj_entities.end(); b_it++) {
    // RefEntity adj_entity(moab,*eit);
    // std::cerr << "\t" << adj_entity << std::endl;
    if (from_ref_entity.getBitRefLevel() !=
        *b_it /*adj_entity.getBitRefLevel()*/) {
      eit = adj_entities.erase(eit);
    } else {
      eit++;
    }
  }
  if (b_it != bit_levels.end()) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Data inconsistency");
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::getAdjacenciesAny(const EntityHandle from_entity,
                                                const int to_dimension,
                                                Range &adj_entities) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBeginHot;
  RefEntity from_ref_entity(m_field.get_basic_entity_data_ptr(), from_entity);
  // std::cerr << "from:\n";
  // std::cerr << from_ref_entity << std::endl;
  rval =
      moab.get_adjacencies(&from_entity, 1, to_dimension, false, adj_entities);
  CHKERRQ_MOAB(rval);
  std::vector<BitRefLevel> bit_levels(adj_entities.size());
  rval = moab.tag_get_data(cOre.get_th_RefBitLevel(), adj_entities,
                           &*bit_levels.begin());
  std::vector<BitRefLevel>::iterator b_it = bit_levels.begin();
  Range::iterator eit = adj_entities.begin();
  // std::cerr << "to:\n";
  for (; eit != adj_entities.end(); b_it++) {
    // RefEntity adj_entity(moab,*eit);
    // std::cerr << "\t" << adj_entity << std::endl;
    if (!(from_ref_entity.getBitRefLevel() & (*b_it))
             .any() /*adj_entity.getBitRefLevel()).any()*/) {
      eit = adj_entities.erase(eit);
    } else {
      eit++;
    }
  }
  if (b_it != bit_levels.end()) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Data inconsistency");
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::getAdjacencies(
    const Problem *problem_ptr, const EntityHandle *from_entities,
    const int num_entities, const int to_dimension, Range &adj_entities,
    const int operation_type, const int verb) const {
  MoFEMFunctionBeginHot;
  BitRefLevel bit = problem_ptr->getBitRefLevel();
  ierr = getAdjacencies(bit, from_entities, num_entities, to_dimension,
                        adj_entities, operation_type);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::getAdjacencies(
    const BitRefLevel &bit, const EntityHandle *from_entities,
    const int num_entities, const int to_dimension, Range &adj_entities,
    const int operation_type, const int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBeginHot;
  if (verb > 0) {
    std::ostringstream ss;
    ss << "from: " << bit << std::endl << "to: " << std::endl;
    PetscPrintf(m_field.get_comm(), ss.str().c_str());
  }
  rval = moab.get_adjacencies(from_entities, num_entities, to_dimension, false,
                              adj_entities, operation_type);
  CHKERRQ_MOAB(rval);
  std::vector<BitRefLevel> bit_levels(adj_entities.size());
  rval = moab.tag_get_data(cOre.get_th_RefBitLevel(), adj_entities,
                           &*bit_levels.begin());
  CHKERRQ_MOAB(rval);
  std::vector<BitRefLevel>::iterator b_it = bit_levels.begin();
  Range::iterator eit = adj_entities.begin();
  // std::cerr << "to:\n";
  for (; eit != adj_entities.end(); b_it++) {
    if (verb > 0) {
      RefEntity adj_entity(m_field.get_basic_entity_data_ptr(), *eit);
      std::ostringstream ss;
      ss << "\t" << adj_entity << std::endl;
      PetscPrintf(m_field.get_comm(), ss.str().c_str());
    }
    if (!((*b_it) & bit).any()) {
      eit = adj_entities.erase(eit);
    } else {
      eit++;
    }
  }
  if (b_it != bit_levels.end()) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Data inconsistency");
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::updateMeshsetByEntitiesChildren(
    const EntityHandle parent, const BitRefLevel &child_bit,
    const EntityHandle child, EntityType child_type, const bool recursive,
    int verb) {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.get_ref_ents(&ref_ents_ptr);
  CHKERRG(ierr);
  Range ents;
  rval = moab.get_entities_by_handle(parent, ents, recursive);
  if (rval != MB_SUCCESS) {
    std::cerr << parent << std::endl;
    std::cerr << moab.type_from_handle(parent) << " " << MBENTITYSET
              << std::endl;
  }
  CHKERRQ_MOAB(rval);

  typedef RefEntity_multiIndex::index<
      Composite_ParentEnt_And_EntType_mi_tag>::type RefEntsByComposite;
  RefEntsByComposite &ref_ents =
      const_cast<RefEntity_multiIndex *>(ref_ents_ptr)
          ->get<Composite_ParentEnt_And_EntType_mi_tag>();
  for (Range::iterator eit = ents.begin(); eit != ents.end(); eit++) {
    if (verb > 2) {
      std::ostringstream ss;
      ss << "ent " << *eit << std::endl;
      ;
      PetscPrintf(m_field.get_comm(), ss.str().c_str());
    }
    RefEntsByComposite::iterator miit =
        ref_ents.lower_bound(boost::make_tuple(*eit, child_type));
    RefEntsByComposite::iterator hi_miit =
        ref_ents.upper_bound(boost::make_tuple(*eit, child_type));
    for (; miit != hi_miit; miit++) {
      if (verb > 2) {
        std::ostringstream ss;
        ss << "any bit " << *miit << std::endl;
        ;
        PetscPrintf(m_field.get_comm(), ss.str().c_str());
      }
      if (((*miit)->getBitRefLevel() & child_bit).any()) {
        EntityHandle ref_ent = (*miit)->getRefEnt();
        if (ref_ent == *eit)
          continue;
        if (ref_ent == 0) {
          SETERRQ(m_field.get_comm(), MOFEM_IMPOSIBLE_CASE,
                  "this should not happen");
        }
        if (moab.type_from_handle(*eit) == MBENTITYSET) {
          SETERRQ(m_field.get_comm(), MOFEM_IMPOSIBLE_CASE,
                  "this should not happen");
        }
        rval = moab.add_entities(child, &ref_ent, 1);
        CHKERRQ_MOAB(rval);
        if (verb > 1) {
          std::ostringstream ss;
          ss << "good bit " << *miit << std::endl;
          PetscPrintf(m_field.get_comm(), ss.str().c_str());
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::updateFieldMeshsetByEntitiesChildren(
    const BitRefLevel &child_bit, int verb) {
  MoFEM::Interface &m_field = cOre;
  const Field_multiIndex *fields_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.get_fields(&fields_ptr);
  CHKERRG(ierr);
  Field_multiIndex::iterator fit = fields_ptr->begin();
  for (; fit != fields_ptr->end(); fit++) {
    EntityHandle meshset = (*fit)->getMeshset();
    ierr = updateMeshsetByEntitiesChildren(meshset, child_bit, meshset, MBTET,
                                           false, verb);
    CHKERRG(ierr);
    ierr = updateMeshsetByEntitiesChildren(meshset, child_bit, meshset, MBTRI,
                                           false, verb);
    CHKERRG(ierr);
    ierr = updateMeshsetByEntitiesChildren(meshset, child_bit, meshset, MBEDGE,
                                           false, verb);
    CHKERRG(ierr);
    ierr = updateMeshsetByEntitiesChildren(meshset, child_bit, meshset,
                                           MBVERTEX, false, verb);
    CHKERRG(ierr);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::updateFieldMeshsetByEntitiesChildren(
    const std::string name, const BitRefLevel &child_bit, int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBeginHot;
  try {
    EntityHandle meshset = m_field.get_field_structure(name)->getMeshset();
    ierr = updateMeshsetByEntitiesChildren(meshset, child_bit, meshset, MBTET,
                                           false, verb);
    CHKERRG(ierr);
    ierr = updateMeshsetByEntitiesChildren(meshset, child_bit, meshset, MBTRI,
                                           false, verb);
    CHKERRG(ierr);
    ierr = updateMeshsetByEntitiesChildren(meshset, child_bit, meshset, MBEDGE,
                                           false, verb);
    CHKERRG(ierr);
    ierr = updateMeshsetByEntitiesChildren(meshset, child_bit, meshset,
                                           MBVERTEX, false, verb);
    CHKERRG(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(m_field.get_comm(), e.errorCode, e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::updateFiniteElementMeshsetByEntitiesChildren(
    const std::string name, const BitRefLevel &child_bit,
    const EntityType fe_ent_type, int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBeginHot;
  try {
    EntityHandle meshset = m_field.get_finite_element_meshset(name);
    ierr = updateMeshsetByEntitiesChildren(meshset, child_bit, meshset,
                                           fe_ent_type, false, verb);
    CHKERRG(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(m_field.get_comm(), e.errorCode, e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::updateRange(const Range &parent_ents,
                                          Range &child_ents) {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.get_ref_ents(&ref_ents_ptr);
  CHKERRG(ierr);
  typedef RefEntity_multiIndex::index<Ent_Ent_mi_tag>::type RefEntsByParent;
  RefEntsByParent &ref_ents =
      const_cast<RefEntity_multiIndex *>(ref_ents_ptr)->get<Ent_Ent_mi_tag>();
  for (Range::iterator eit = parent_ents.begin(); eit != parent_ents.end();
       eit++) {
    RefEntsByParent::iterator miit = ref_ents.lower_bound(*eit);
    RefEntsByParent::iterator hi_miit = ref_ents.upper_bound(*eit);
    for (; miit != hi_miit; miit++) {
      EntityHandle ref_ent = (*miit)->getRefEnt();
      if (ref_ent == *eit)
        continue;
      if (ref_ent == 0) {
        SETERRQ(m_field.get_comm(), MOFEM_IMPOSIBLE_CASE,
                "this should not happen");
      }
      if (moab.type_from_handle(*eit) == MBENTITYSET) {
        SETERRQ(m_field.get_comm(), MOFEM_IMPOSIBLE_CASE,
                "this should not happen");
      }
      child_ents.insert(ref_ent);
    }
  }
  MoFEMFunctionReturnHot(0);
}
}
