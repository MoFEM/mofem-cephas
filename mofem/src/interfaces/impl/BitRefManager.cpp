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

/// tool class with methods used more than twp times
struct SetBitRefLevelTool {

  const BitRefLevel &bIt;                     ///< bit to set
  const RefEntity_multiIndex *refEntsPtr;     ///< access to ents database
  const RefElement_multiIndex *refElementPtr; ///< access to fe database

  boost::shared_ptr<BasicEntityData> &baseEntData; ///< base entity data

  /// constrictor
  SetBitRefLevelTool(MoFEM::Interface &m_field, const BitRefLevel &bit,
                     const RefEntity_multiIndex *ref_ents_ptr,
                     const RefElement_multiIndex *ref_element_ptr)
      : bIt(bit), refEntsPtr(ref_ents_ptr), refElementPtr(ref_element_ptr),
        baseEntData(m_field.get_basic_entity_data_ptr()) {}

  /// find entities and change entity bit if in database
  MoFEMErrorCode findEntsToAdd(EntityHandle f, EntityHandle s,
                               Range &seed_ents_range) const {
    MoFEMFunctionBeginHot;

    seed_ents_range.insert(f, s);
    // get lower bound of multi-index
    auto rit = refEntsPtr->lower_bound(f);
    if (rit == refEntsPtr->end()) {
      // all enties in range are added to database
      MoFEMFunctionReturnHot(0);
    } else {

      // some entities from range are in database
      auto hi_rit = refEntsPtr->upper_bound(s);

      Range to_erase;
      insertOrdered(to_erase, RefEntExtractor(), rit, hi_rit);
      if (bIt.any())
        for (; rit != hi_rit; ++rit)
          const_cast<BitRefLevel &>((*rit)->getBitRefLevel()) |= bIt;

      Range::iterator lo, hi = seed_ents_range.begin();
      for (auto pt = to_erase.pair_begin(); pt != to_erase.pair_end(); ++pt) {
        lo = seed_ents_range.lower_bound(hi, seed_ents_range.end(), pt->first);
        if (lo != seed_ents_range.end()) {
          hi = seed_ents_range.upper_bound(lo, seed_ents_range.end(),
                                           pt->second);
          seed_ents_range.erase(lo, hi);
        } else
          break;
      }

    }

    MoFEMFunctionReturnHot(0);
  }

  /// add entities to database
  MoFEMErrorCode addEntsToDatabase(const Range &seed_ents_range) const {
    MoFEMFunctionBeginHot;
    std::vector<boost::shared_ptr<RefEntity>> shared_ref_ents_vec;
    shared_ref_ents_vec.reserve(seed_ents_range.size());
    for (Range::const_pair_iterator pit = seed_ents_range.pair_begin();
         pit != seed_ents_range.pair_end(); pit++) {
      // add entities to database
      EntityHandle f = pit->first;
      EntityHandle s = pit->second;
      boost::shared_ptr<std::vector<RefEntity>> ref_ents_vec(
          new std::vector<RefEntity>());
      ref_ents_vec->reserve(s - f + 1);
      for (auto f : Range(f, s))
        ref_ents_vec->emplace_back(baseEntData, f);

      // Set bits to range
      if (bIt.any()) {
        boost::shared_ptr<std::vector<const void *>> bits_by_ptr(
            new std::vector<const void *>());
        bits_by_ptr->resize(s - f + 1);
        CHKERR baseEntData->moab.tag_get_by_ptr(
            baseEntData->th_RefBitLevel, Range(f, s), &*bits_by_ptr->begin());
        for (auto &v_bit_ptr : *bits_by_ptr)
          const_cast<BitRefLevel &>(
              *(static_cast<const BitRefLevel *>(v_bit_ptr))) |= bIt;
      }

      for (auto &re : *ref_ents_vec)
        shared_ref_ents_vec.emplace_back(ref_ents_vec, &re);
    }
    if (!shared_ref_ents_vec.empty()) {
      int s0 = refEntsPtr->size();
      const_cast<RefEntity_multiIndex *>(refEntsPtr)
          ->insert(shared_ref_ents_vec.begin(), shared_ref_ents_vec.end());
      if ((refEntsPtr->size() - s0) != shared_ref_ents_vec.size()) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Data inconsistency %d != %d", refEntsPtr->size() - s0,
                 shared_ref_ents_vec.size());
      }
    }
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode findElementsToAdd(EntityHandle f, EntityHandle s,
                                   Range &seed_fe_range) const {
    MoFEMFunctionBeginHot;
    seed_fe_range.insert(f, s);
    RefElement_multiIndex::iterator rit, hi_rit;
    // get lower bound of multi-index
    rit = refElementPtr->lower_bound(f);
    if (rit == refElementPtr->end()) {
      // all enties in range are added to database
      MoFEMFunctionReturnHot(0);
    } else {
      // some entities from range are in database
      hi_rit = refElementPtr->upper_bound(s);
      for (; rit != hi_rit; ++rit) {
        seed_fe_range.erase(rit->get()->getRefEnt());
      }
    }
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode addElementsToDatabase(Range &seed_fe_range) const {
    MoFEMFunctionBeginHot;
    std::vector<boost::shared_ptr<RefElement>> shared_ref_fe_vec;
    shared_ref_fe_vec.reserve(seed_fe_range.size());
    // create ref entity instances
    for (Range::const_pair_iterator pit = seed_fe_range.const_pair_begin();
         pit != seed_fe_range.const_pair_end(); ++pit) {
      RefEntity_multiIndex::iterator rit, hi_rit;
      rit = refEntsPtr->lower_bound(pit->first);
      hi_rit = refEntsPtr->upper_bound(pit->second);
      if (static_cast<int>(std::distance(rit, hi_rit)) !=
          static_cast<int>(pit->second - pit->first + 1)) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "data inconsistency %d != %d", std::distance(rit, hi_rit),
                 pit->second - pit->first + 1);
      }
      switch ((*rit)->getEntType()) {
      case MBVERTEX: {
        boost::shared_ptr<std::vector<RefElement_VERTEX>> ref_fe_vec =
            boost::make_shared<std::vector<RefElement_VERTEX>>();
        ref_fe_vec->reserve(pit->second - pit->first + 1);
        for (; rit != hi_rit; ++rit) {
          ref_fe_vec->push_back(RefElement_VERTEX(*rit));
          shared_ref_fe_vec.push_back(
              boost::shared_ptr<RefElement>(ref_fe_vec, &ref_fe_vec->back()));
        }
      } break;
      case MBEDGE: {
        boost::shared_ptr<std::vector<RefElement_EDGE>> ref_fe_vec =
            boost::make_shared<std::vector<RefElement_EDGE>>();
        ref_fe_vec->reserve(pit->second - pit->first + 1);
        for (; rit != hi_rit; ++rit) {
          ref_fe_vec->push_back(RefElement_EDGE(*rit));
          shared_ref_fe_vec.push_back(
              boost::shared_ptr<RefElement>(ref_fe_vec, &ref_fe_vec->back()));
        }
      } break;
      case MBTRI: {
        boost::shared_ptr<std::vector<RefElement_TRI>> ref_fe_vec =
            boost::make_shared<std::vector<RefElement_TRI>>();
        ref_fe_vec->reserve(pit->second - pit->first + 1);
        for (; rit != hi_rit; ++rit) {
          ref_fe_vec->push_back(RefElement_TRI(*rit));
          shared_ref_fe_vec.push_back(
              boost::shared_ptr<RefElement>(ref_fe_vec, &ref_fe_vec->back()));
        }
      } break;
      case MBTET: {
        boost::shared_ptr<std::vector<RefElement_TET>> ref_fe_vec =
            boost::make_shared<std::vector<RefElement_TET>>();
        ref_fe_vec->reserve(pit->second - pit->first + 1);
        for (; rit != hi_rit; ++rit) {
          ref_fe_vec->push_back(RefElement_TET(*rit));
          shared_ref_fe_vec.push_back(
              boost::shared_ptr<RefElement>(ref_fe_vec, &ref_fe_vec->back()));
        }
      } break;
      case MBPRISM: {
        boost::shared_ptr<std::vector<RefElement_PRISM>> ref_fe_vec =
            boost::make_shared<std::vector<RefElement_PRISM>>();
        ref_fe_vec->reserve(pit->second - pit->first + 1);
        for (; rit != hi_rit; ++rit) {
          ref_fe_vec->push_back(RefElement_PRISM(*rit));
          shared_ref_fe_vec.push_back(
              boost::shared_ptr<RefElement>(ref_fe_vec, &ref_fe_vec->back()));
        }
      } break;
      case MBENTITYSET: {
        boost::shared_ptr<std::vector<RefElement_MESHSET>> ref_fe_vec =
            boost::make_shared<std::vector<RefElement_MESHSET>>();
        ref_fe_vec->reserve(pit->second - pit->first + 1);
        for (; rit != hi_rit; ++rit) {
          ref_fe_vec->push_back(RefElement_MESHSET(*rit));
          shared_ref_fe_vec.push_back(
              boost::shared_ptr<RefElement>(ref_fe_vec, &ref_fe_vec->back()));
        }
      } break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
      }
    }
    // add shared pointers to database
    const_cast<RefElement_multiIndex *>(refElementPtr)
        ->insert(shared_ref_fe_vec.begin(), shared_ref_fe_vec.end());
    MoFEMFunctionReturnHot(0);
  }
};

MoFEMErrorCode BitRefManager::setBitRefLevel(const Range &ents,
                                             const BitRefLevel bit,
                                             const bool only_tets,
                                             int verb) const {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ents_ptr;
  const RefElement_multiIndex *ref_fe_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_ref_ents(&ref_ents_ptr);
  CHKERR m_field.get_ref_finite_elements(&ref_fe_ptr);

  if (verb > VERBOSE) 
    PetscSynchronizedPrintf(PETSC_COMM_SELF, "nb. entities to add %d\n",
                            ents.size());

  CHKERR setElementsBitRefLevel(ents, bit, verb);

  if (!ents.empty()) {
    for (int d = 3; d >= 1; --d) {
      Range dim_ents;
      if (only_tets && d == 3) {
        dim_ents = ents.subset_by_type(MBTET);
      } else {
        dim_ents = ents.subset_by_dimension(d);
      }
      if (!dim_ents.empty()) {
        for (int dd = 0; dd < d; ++dd) {
          Range adj_ents;
          CHKERR m_field.get_moab().get_adjacencies(
              dim_ents, dd, true, adj_ents, moab::Interface::UNION);
          for (Range::pair_iterator pit = adj_ents.pair_begin();
               pit != adj_ents.pair_end(); ++pit) {
            Range seed_ents_range;
            // get first and last element of range
            EntityHandle f = pit->first;
            EntityHandle s = pit->second;
            CHKERR SetBitRefLevelTool(m_field, bit, ref_ents_ptr, ref_fe_ptr)
                .findEntsToAdd(f, s, seed_ents_range);
            if (!seed_ents_range.empty()) {
              CHKERR SetBitRefLevelTool(m_field, bit, ref_ents_ptr, ref_fe_ptr)
                  .addEntsToDatabase(seed_ents_range);
            }
          }
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::setElementsBitRefLevel(const Range &ents,
                                                     const BitRefLevel bit,
                                                     int verb) const {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ents_ptr;
  const RefElement_multiIndex *ref_fe_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_ref_ents(&ref_ents_ptr);
  CHKERR m_field.get_ref_finite_elements(&ref_fe_ptr);

  for (Range::const_pair_iterator pit = ents.pair_begin();
       pit != ents.pair_end(); pit++) {
    // get first and last element of range
    EntityHandle f = pit->first;
    EntityHandle s = pit->second;
    Range seed_ents_range; // entities seeded not in database
    // find ents to add
    CHKERR SetBitRefLevelTool(m_field, bit, ref_ents_ptr, ref_fe_ptr)
        .findEntsToAdd(f, s, seed_ents_range);
    // add elements
    if (!seed_ents_range.empty()) {
      CHKERR SetBitRefLevelTool(m_field, bit, ref_ents_ptr, ref_fe_ptr)
          .addEntsToDatabase(seed_ents_range);
    }
    Range seed_fe_range;
    CHKERR SetBitRefLevelTool(m_field, bit, ref_ents_ptr, ref_fe_ptr)
        .findElementsToAdd(f, s, seed_fe_range);
    if (!seed_fe_range.empty()) {
      CHKERR SetBitRefLevelTool(m_field, bit, ref_ents_ptr, ref_fe_ptr)
          .addElementsToDatabase(seed_fe_range);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::setEntitiesBitRefLevel(const Range &ents,
                                                     const BitRefLevel bit,
                                                     int verb) const {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ents_ptr;
  const RefElement_multiIndex *ref_fe_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_ref_ents(&ref_ents_ptr);
  CHKERR m_field.get_ref_finite_elements(&ref_fe_ptr);

  for (Range::const_pair_iterator pit = ents.pair_begin();
       pit != ents.pair_end(); pit++) {
    // get first and last element of range
    EntityHandle f = pit->first;
    EntityHandle s = pit->second;
    Range seed_ents_range; // entities seeded not in database
    // find ents to add
    CHKERR SetBitRefLevelTool(m_field, bit, ref_ents_ptr, ref_fe_ptr)
        .findEntsToAdd(f, s, seed_ents_range);
    // add elements
    if (!seed_ents_range.empty()) {
      CHKERR SetBitRefLevelTool(m_field, bit, ref_ents_ptr, ref_fe_ptr)
          .addEntsToDatabase(seed_ents_range);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::addToDatabaseBitRefLevelByType(
    const EntityType type, const BitRefLevel bit, const BitRefLevel mask,
    int verb) const {
  MoFEMFunctionBegin;
  Range ents;
  CHKERR getEntitiesByTypeAndRefLevel(bit, mask, type, ents);
  CHKERR setBitRefLevel(ents, BitRefLevel());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::addToDatabaseBitRefLevelByDim(
    const int dim, const BitRefLevel bit, const BitRefLevel mask,
    int verb) const {
  MoFEMFunctionBegin;
  Range ents;
  CHKERR getEntitiesByDimAndRefLevel(bit, mask, dim, ents);
  CHKERR setBitRefLevel(ents, BitRefLevel());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::setBitLevelToMeshset(const EntityHandle meshset,
                                                   const BitRefLevel bit,
                                                   int verb) const {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ents_ptr;
  const RefElement_multiIndex *ref_fe_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_ref_ents(&ref_ents_ptr);
  CHKERR m_field.get_ref_finite_elements(&ref_fe_ptr);
  // Add ref entity
  std::pair<RefEntity_multiIndex::iterator, bool> p_ent =
      const_cast<RefEntity_multiIndex *>(ref_ents_ptr)
          ->insert(boost::shared_ptr<RefEntity>(
              new RefEntity(m_field.get_basic_entity_data_ptr(), meshset)));
  *(const_cast<RefEntity *>(p_ent.first->get())->getBitRefLevelPtr()) |= bit;
  // Add ref element
  boost::shared_ptr<RefElement> fe_ptr =
      boost::shared_ptr<RefElement>(new RefElement_MESHSET(*p_ent.first));
  std::pair<RefElement_multiIndex::iterator, bool> p_fe =
      const_cast<RefElement_multiIndex *>(ref_fe_ptr)->insert(fe_ptr);
  if (verb > 0) {
    std::ostringstream ss;
    ss << "add meshset as ref_ent " << **p_fe.first << std::endl;
    PetscPrintf(PETSC_COMM_SELF, ss.str().c_str());
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::setBitRefLevelByDim(const EntityHandle meshset,
                                                  const int dim,
                                                  const BitRefLevel bit,
                                                  int verb) const {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  Range ents;
  CHKERR m_field.get_moab().get_entities_by_dimension(meshset, dim, ents,
                                                      false);
  CHKERR setBitRefLevel(ents, bit, false, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::setBitRefLevelByType(const EntityHandle meshset,
                                                   const EntityType type,
                                                   const BitRefLevel bit,
                                                   int verb) const {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  Range ents;
  CHKERR m_field.get_moab().get_entities_by_type(meshset, type, ents, false);
  CHKERR setBitRefLevel(ents, bit, false, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::addBitRefLevel(const Range &ents,
                                             const BitRefLevel bit,
                                             int verb) const {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  std::vector<const BitRefLevel *> ents_bits_vec;
  CHKERR RefEntity::getBitRefLevel(m_field.get_moab(), ents, ents_bits_vec);
  for (auto it : ents_bits_vec)
    const_cast<BitRefLevel &>(*it) |= bit;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::addBitRefLevelByDim(const EntityHandle meshset,
                                                  const int dim,
                                                  const BitRefLevel bit,
                                                  int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  Range ents, adj;
  MoFEMFunctionBegin;
  CHKERR moab.get_entities_by_dimension(meshset, dim, ents, true);
  for (int dd = dim - 1; dd >= 0; dd--) {
    CHKERR moab.get_adjacencies(ents, dd, false, adj, moab::Interface::UNION);
  }
  ents.merge(adj);
  if (verb == VERY_NOISY) {
    cerr << ents << endl;
  }
  CHKERR addBitRefLevel(ents, bit, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::setNthBitRefLevel(const Range &ents, const int n,
                                                const bool b, int verb) const {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  std::vector<const BitRefLevel *> ents_bits_vec;
  CHKERR RefEntity::getBitRefLevel(m_field.get_moab(), ents, ents_bits_vec);
  for (std::vector<const BitRefLevel *>::iterator it = ents_bits_vec.begin();
       it != ents_bits_vec.end(); ++it) {
    const_cast<BitRefLevel &>(**it)[n] = b;
  }
  if (verb == VERY_NOISY) {
    cerr << ents << endl;
  }
  MoFEMFunctionReturn(0);
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
    (*const_cast<RefEntity *>(dit->get())->getBitRefLevelPtr())[n] = b;
    if (verb >= VERY_VERBOSE) {
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
  RefEntity_change_right_shift right_shift(1, mask);
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
      right_shift(const_cast<boost::shared_ptr<RefEntity> &>(*ent_it));
      if (verb >= VERY_NOISY) {
        std::cerr << (*ent_it)->getBitRefLevel() << std::endl;
      }
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::writeBitLevelByType(
    const BitRefLevel bit, const BitRefLevel mask, const EntityType type,
    const char *file_name, const char *file_type, const char *options,
    const bool check_for_empty) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  Range ents;
  CHKERR getEntitiesByTypeAndRefLevel(bit, mask, type, ents);
  if (check_for_empty && ents.empty())
    MoFEMFunctionReturnHot(0);
  EntityHandle meshset;
  CHKERR moab.create_meshset(MESHSET_SET, meshset);
  CHKERR moab.add_entities(meshset, ents);
  CHKERR moab.write_file(file_name, file_type, options, &meshset, 1);
  CHKERR moab.delete_entities(&meshset, 1);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::writeEntitiesNotInDatabase(
    const char *file_name, const char *file_type, const char *options,
    const bool check_for_empty) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  EntityHandle meshset;
  Range ents;
  CHKERR getAllEntitiesNotInDatabase(ents);
  if (check_for_empty && ents.empty())
    MoFEMFunctionReturnHot(0);
  CHKERR moab.create_meshset(MESHSET_SET, meshset);
  CHKERR moab.add_entities(meshset, ents);
  CHKERR moab.write_file(file_name, file_type, options, &meshset, 1);
  CHKERR moab.delete_entities(&meshset, 1);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::writeEntitiesAllBitLevelsByType(
    const BitRefLevel mask, const EntityType type, const char *file_name,
    const char *file_type, const char *options) {
  MoFEMFunctionBegin;
  for (int ll = 0; ll != BITREFLEVEL_SIZE; ++ll) {
    std::string name = boost::lexical_cast<std::string>(ll) + "_" + file_name;
    CHKERR writeBitLevelByType(BitRefLevel().set(ll), mask, type, name.c_str(),
                               file_type, options, true);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::getEntitiesByTypeAndRefLevel(
    const BitRefLevel bit, const BitRefLevel mask, const EntityType type,
    const EntityHandle meshset, int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  Range ents;
  CHKERR getEntitiesByTypeAndRefLevel(bit, mask, type, ents, verb);
  CHKERR moab.add_entities(meshset, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::filterEntitiesByRefLevel(const BitRefLevel bit,
                                                       const BitRefLevel mask,
                                                       Range &ents,
                                                       int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBeginHot;

  std::vector<EntityHandle> ents_vec;
  ents_vec.reserve(ents.size());

  std::vector<BitRefLevel *> tags_bits_ptr_vec(ents.size());

  Range swap_ents;
  auto hint = swap_ents.begin();

  for (Range::pair_iterator p_eit = ents.pair_begin(); p_eit != ents.pair_end();
       ++p_eit) {

    EntityHandle f = p_eit->first;
    const EntityHandle s = p_eit->second;

    // get bits on entities
    rval = moab.tag_get_by_ptr(cOre.get_th_RefBitLevel(), Range(f, s),
                             (const void **)(&*tags_bits_ptr_vec.begin()));

    if (rval == MB_SUCCESS) {

      auto bit_it = tags_bits_ptr_vec.begin();

      auto check = [&bit, &mask](const auto &entity_bit) -> bool {
        return

            (entity_bit & bit).any() &&

            ((entity_bit & mask) == entity_bit);
      };

      while (f != s + 1) {

        while (f != s + 1 && !check(**bit_it)) {
          ++bit_it;
          ++f;
        }

        if (f != s + 1) {

          const EntityHandle start = f;

          while (f != (s + 1) && check(**bit_it)) {
            ++bit_it;
            ++f;
          };

          hint = swap_ents.insert(hint, start, f - 1);
        }

      }
    }
  }

  ents.swap(swap_ents);

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::getEntitiesByTypeAndRefLevel(
    const BitRefLevel bit, const BitRefLevel mask, const EntityType type,
    Range &ents, int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  CHKERR moab.get_entities_by_type(0, type, ents, false);
  CHKERR filterEntitiesByRefLevel(bit, mask, ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::getEntitiesByDimAndRefLevel(
    const BitRefLevel bit, const BitRefLevel mask, const int dim,
    const EntityHandle meshset, int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  Range ents;
  CHKERR getEntitiesByDimAndRefLevel(bit, mask, dim, ents, verb);
  CHKERR moab.add_entities(meshset, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::getEntitiesByDimAndRefLevel(
    const BitRefLevel bit, const BitRefLevel mask, const int dim, Range &ents,
    int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  CHKERR moab.get_entities_by_dimension(0, dim, ents, false);
  CHKERR filterEntitiesByRefLevel(bit, mask, ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::getEntitiesByRefLevel(const BitRefLevel bit,
                                                    const BitRefLevel mask,
                                                    const EntityHandle meshset,
                                                    const int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  Range ents;
  CHKERR getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR moab.add_entities(meshset, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::getEntitiesByRefLevel(const BitRefLevel bit,
                                                    const BitRefLevel mask,
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

MoFEMErrorCode BitRefManager::getEntitiesByParentType(const BitRefLevel bit,
                                                      const BitRefLevel mask,
                                                      const EntityType type,
                                                      Range &ents,
                                                      const int verb) const {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_ref_ents(&ref_ents_ptr);
  auto &ref_ents = ref_ents_ptr->get<ParentEntType_mi_tag>();
  auto it = ref_ents.lower_bound(type);
  auto hi_it = ref_ents.upper_bound(type);
  std::vector<EntityHandle> ents_vec;
  ents_vec.reserve(std::distance(it, hi_it));
  for (; it != hi_it; it++) {
    const BitRefLevel &ent_bit = it->get()->getBitRefLevel();
    if ((ent_bit & mask) == ent_bit && (ent_bit & bit).any())
      ents_vec.emplace_back(it->get()->getRefEnt());
  }
  ents.insert_list(ents_vec.begin(), ents_vec.end());
  if (verb > NOISY)
    cerr << "getEntitiesByParentType: " << ents << endl;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::getAllEntitiesNotInDatabase(Range &ents) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  rval = moab.get_entities_by_handle(0, ents, false);
  CHKERRG(rval);
  ents = subtract(ents, ents.subset_by_type(MBENTITYSET));
  ierr = filterEntitiesNotInDatabase(ents);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::filterEntitiesNotInDatabase(Range &ents) const {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.get_ref_ents(&ref_ents_ptr);
  CHKERRG(ierr);
  Range::iterator eit = ents.begin();
  for (; eit != ents.end();) {
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
  MoFEMFunctionBegin;
  BitRefLevel bit_from_entity;
  CHKERR moab.tag_get_data(cOre.get_th_RefBitLevel(), &from_entity, 1,
                           &bit_from_entity);
  CHKERR moab.get_adjacencies(&from_entity, 1, to_dimension, false,
                              adj_entities);
  std::vector<BitRefLevel> bit_levels(adj_entities.size());
  CHKERR moab.tag_get_data(cOre.get_th_RefBitLevel(), adj_entities,
                           &*bit_levels.begin());
  auto b_it = bit_levels.begin();
  auto eit = adj_entities.begin();
  for (; eit != adj_entities.end(); b_it++) {
    if (bit_from_entity != *b_it) {
      eit = adj_entities.erase(eit);
    } else {
      eit++;
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::getAdjacenciesAny(const EntityHandle from_entity,
                                                const int to_dimension,
                                                Range &adj_entities) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBeginHot;
  BitRefLevel bit_from_entity;
  CHKERR moab.tag_get_data(cOre.get_th_RefBitLevel(), &from_entity, 1,
                           &bit_from_entity);
  CHKERR moab.get_adjacencies(&from_entity, 1, to_dimension, false,
                              adj_entities);
  std::vector<BitRefLevel> bit_levels(adj_entities.size());
  CHKERR moab.tag_get_data(cOre.get_th_RefBitLevel(), adj_entities,
                           &*bit_levels.begin());
  std::vector<BitRefLevel>::iterator b_it = bit_levels.begin();
  Range::iterator eit = adj_entities.begin();
  for (; eit != adj_entities.end(); b_it++) {
    if (!(bit_from_entity & (*b_it)).any()) {
      eit = adj_entities.erase(eit);
    } else {
      eit++;
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitRefManager::getAdjacencies(
    const Problem *problem_ptr, const EntityHandle *from_entities,
    const int num_entities, const int to_dimension, Range &adj_entities,
    const int operation_type, const int verb) const {
  MoFEMFunctionBegin;
  BitRefLevel bit = problem_ptr->getBitRefLevel();
  CHKERR getAdjacencies(bit, from_entities, num_entities, to_dimension,
                        adj_entities, operation_type);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::getAdjacencies(
    const BitRefLevel bit, const EntityHandle *from_entities,
    const int num_entities, const int to_dimension, Range &adj_entities,
    const int operation_type, const int verb) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  if (verb > QUIET) {
    std::ostringstream ss;
    ss << "from: " << bit << std::endl << "to: " << std::endl;
    PetscPrintf(PETSC_COMM_SELF, ss.str().c_str());
  }
  CHKERR moab.get_adjacencies(from_entities, num_entities, to_dimension, false,
                              adj_entities, operation_type);
  std::vector<BitRefLevel> bit_levels(adj_entities.size());
  CHKERR moab.tag_get_data(cOre.get_th_RefBitLevel(), adj_entities,
                           &*bit_levels.begin());
  std::vector<BitRefLevel>::iterator b_it = bit_levels.begin();
  // std::cerr << "to:\n";
  for (Range::iterator eit = adj_entities.begin(); eit != adj_entities.end();
       b_it++) {
    if (verb > VERBOSE) {
      RefEntity adj_entity(m_field.get_basic_entity_data_ptr(), *eit);
      std::ostringstream ss;
      ss << "\t" << adj_entity.getBitRefLevel() << " : " << adj_entity
         << std::endl;
      PetscPrintf(PETSC_COMM_SELF, ss.str().c_str());
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
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::updateMeshsetByEntitiesChildren(
    const EntityHandle parent, const BitRefLevel &parent_bit,
    const BitRefLevel &parent_mask, const BitRefLevel &child_bit,
    const BitRefLevel &child_mask, const EntityHandle child,
    EntityType child_type, const bool recursive, int verb) {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_ref_ents(&ref_ents_ptr);
  Range ents;
  CHKERR moab.get_entities_by_handle(parent, ents, recursive);
  CHKERR filterEntitiesByRefLevel(parent_bit, parent_mask, ents, verb);
  if (verb >= VERY_VERBOSE) {
    std::cerr << "Parnets:" << endl << parent << std::endl;
  }
  typedef RefEntity_multiIndex::index<
      Composite_ParentEnt_And_EntType_mi_tag>::type RefEntsByComposite;
  RefEntsByComposite &ref_ents =
      const_cast<RefEntity_multiIndex *>(ref_ents_ptr)
          ->get<Composite_ParentEnt_And_EntType_mi_tag>();
  Range children_ents;
  for (Range::pair_iterator pit = ents.pair_begin(); pit != ents.pair_end();
       ++pit) {
    const EntityHandle f = pit->first;
    const EntityHandle s = pit->second;
    auto lo_mit = ref_ents.lower_bound(boost::make_tuple(child_type, f));
    auto hi_mit = ref_ents.upper_bound(boost::make_tuple(child_type, s));
    insertOrdered(children_ents, RefEntExtractor(), lo_mit, hi_mit);
  }
  CHKERR filterEntitiesByRefLevel(child_bit, child_mask, children_ents, verb);
  if (verb >= VERY_VERBOSE) {
    std::cerr << "Children: " << endl << parent << std::endl;
  }
  CHKERR moab.add_entities(child, children_ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::updateMeshsetByEntitiesChildren(
    const EntityHandle parent, const BitRefLevel &child_bit,
    const EntityHandle child, EntityType child_type, const bool recursive,
    int verb) {
  MoFEMFunctionBegin;
  CHKERR updateMeshsetByEntitiesChildren(
      parent, BitRefLevel().set(), BitRefLevel().set(), child_bit, child_bit,
      child, child_type, recursive, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::updateFieldMeshsetByEntitiesChildren(
    const BitRefLevel &child_bit, int verb) {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  const Field_multiIndex *fields_ptr;
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_fields(&fields_ptr);
  CHKERR m_field.get_ref_ents(&ref_ents_ptr);

  for (auto &fit : (*fields_ptr)) {

    const EntityHandle meshset = fit->getMeshset();
    Range parent_ents;
    CHKERR moab.get_entities_by_handle(meshset, parent_ents, true);

    if (verb >= VERY_VERBOSE)
      std::cerr << "Parnets:" << endl << parent_ents << std::endl;

    Range children_ents;
    CHKERR updateRange(parent_ents, children_ents);
    CHKERR filterEntitiesByRefLevel(child_bit, BitRefLevel().set(),
                                    children_ents, verb);

    CHKERR moab.add_entities(meshset, children_ents);

    if (verb >= VERY_VERBOSE)
      std::cerr << "Children: " << endl << children_ents << std::endl;
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::updateFieldMeshsetByEntitiesChildren(
    const std::string name, const BitRefLevel &child_bit, int verb) {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  EntityHandle field_meshset = m_field.get_field_structure(name)->getMeshset();

  Range parent_ents;
  CHKERR moab.get_entities_by_handle(field_meshset, parent_ents, true);

  if (verb >= VERY_VERBOSE)
    std::cerr << "Parnets:" << endl << parent_ents << std::endl;

  Range children_ents;
  CHKERR updateRange(parent_ents, children_ents);
  CHKERR filterEntitiesByRefLevel(child_bit, BitRefLevel().set(), children_ents,
                                  verb);

  CHKERR moab.add_entities(field_meshset, children_ents);

  if (verb >= VERY_VERBOSE)
    std::cerr << "Children: " << endl << children_ents << std::endl;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::updateFiniteElementMeshsetByEntitiesChildren(
    const std::string name, const BitRefLevel &child_bit,
    const EntityType fe_ent_type, int verb) {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  EntityHandle meshset = m_field.get_finite_element_meshset(name);
  CHKERR updateMeshsetByEntitiesChildren(meshset, child_bit, meshset,
                                         fe_ent_type, false, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::updateRange(const Range &parent_ents,
                                          Range &child_ents) {
  MoFEM::Interface &m_field = cOre;
  const RefEntity_multiIndex *ref_ents_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_ref_ents(&ref_ents_ptr);
  auto &ref_ents =
      const_cast<RefEntity_multiIndex *>(ref_ents_ptr)->get<Ent_Ent_mi_tag>();
  std::vector<EntityHandle> child_ents_vec;
  child_ents_vec.reserve(ref_ents.size());
  for (Range::const_pair_iterator pit = parent_ents.const_pair_begin();
       pit != parent_ents.const_pair_end(); pit++) {
    auto it = ref_ents.lower_bound(pit->first);
    if (it != ref_ents.end()) {
      auto hi_it = ref_ents.upper_bound(pit->second);
      for (; it != hi_it; it++) {
        if (it->get()->getEntType() == MBENTITYSET) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
                  "this should not happen");
        }
        child_ents_vec.emplace_back((*it)->getRefEnt());
      }
    }
  }
  child_ents.insert_list(child_ents_vec.begin(), child_ents_vec.end());
  MoFEMFunctionReturn(0);
}
} // namespace MoFEM
