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
    : cOre(const_cast<MoFEM::Core &>(core)), dEbug(false) {

  if (!LogManager::checkIfChannelExist("BitRefSelf")) {
    auto core_log = logging::core::get();
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSelf(), "BitRefSelf"));
    LogManager::setLog("BitRefSelf");
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "BitRefWorld"));
    LogManager::setLog("BitRefWorld");
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSync(), "BitRefSync"));
    LogManager::setLog("BitRefSync");
    MOFEM_LOG_TAG("BitRefSelf", "BitRefManager");
    MOFEM_LOG_TAG("BitRefWorld", "BitRefManager");
    MOFEM_LOG_TAG("BitRefSync", "BitRefManager");
  }

  MOFEM_LOG_FUNCTION();
  MOFEM_LOG("BitRefWorld", Sev::noisy) << "BitRefManager interface created";
}

/// tool class with methods used more than twp times
struct SetBitRefLevelTool {

  MoFEM::Interface &mField;

  const BitRefLevel &bIt;                     ///< bit to set
  const RefEntity_multiIndex *refEntsPtr;     ///< access to ents database
  const RefElement_multiIndex *refElementPtr; ///< access to fe database

  boost::shared_ptr<BasicEntityData> &baseEntData; ///< base entity data

  /// constrictor
  SetBitRefLevelTool(MoFEM::Interface &m_field, const BitRefLevel &bit,
                     const RefEntity_multiIndex *ref_ents_ptr,
                     const RefElement_multiIndex *ref_element_ptr)
      : mField(m_field), bIt(bit), refEntsPtr(ref_ents_ptr),
        refElementPtr(ref_element_ptr),
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
  template <int N>
  MoFEMErrorCode addEntsToDatabaseImpl(const Range &seed_ents_range) const {
    MoFEMFunctionBeginHot;
    std::vector<boost::shared_ptr<RefEntity>> shared_ref_ents_vec;
    shared_ref_ents_vec.reserve(seed_ents_range.size());
    for (Range::const_pair_iterator pit = seed_ents_range.pair_begin();
         pit != seed_ents_range.pair_end(); pit++) {
      // add entities to database
      EntityHandle f = pit->first;
      EntityHandle s = pit->second;
      boost::shared_ptr<std::vector<RefEntityTmp<N>>> ref_ents_vec(
          new std::vector<RefEntityTmp<N>>());
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
        shared_ref_ents_vec.emplace_back(ref_ents_vec,
                                         static_cast<RefEntity *>(&re));
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

  MoFEMErrorCode addEntsToDatabase(const Range &seed_ents_range) const {
    MoFEMFunctionBeginHot;

    switch (mField.getValue()) {
      case -1:
        return addEntsToDatabaseImpl<-1>(seed_ents_range);
      case 0:
        return addEntsToDatabaseImpl<0>(seed_ents_range);
      case 1:
        return addEntsToDatabaseImpl<1>(seed_ents_range);
      default:
        SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                 "Core index can vary from -1 to %d", MAX_CORE_TMP);
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
        seed_fe_range.erase(rit->get()->getEnt());
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
      case MBTRI:
      case MBQUAD: {
        boost::shared_ptr<std::vector<RefElementFace>> ref_fe_vec =
            boost::make_shared<std::vector<RefElementFace>>();
        ref_fe_vec->reserve(pit->second - pit->first + 1);
        for (; rit != hi_rit; ++rit) {
          ref_fe_vec->push_back(RefElementFace(*rit));
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
  auto ref_ents_ptr = m_field.get_ref_ents();
  auto ref_fe_ptr = m_field.get_ref_finite_elements();
  MoFEMFunctionBegin;

  MOFEM_LOG_FUNCTION();
  MOFEM_LOG_C("BitRefSelf", Sev::noisy, "nb. entities to add %d", ents.size());

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
            if (!seed_ents_range.empty())
              CHKERR SetBitRefLevelTool(m_field, bit, ref_ents_ptr, ref_fe_ptr)
                  .addEntsToDatabase(seed_ents_range);
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
  auto ref_ents_ptr = m_field.get_ref_ents();
  auto ref_fe_ptr = m_field.get_ref_finite_elements();
  MoFEMFunctionBegin;

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
    if (!seed_ents_range.empty())
      CHKERR SetBitRefLevelTool(m_field, bit, ref_ents_ptr, ref_fe_ptr)
          .addEntsToDatabase(seed_ents_range);

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
  auto ref_ents_ptr = m_field.get_ref_ents();
  auto ref_fe_ptr = m_field.get_ref_finite_elements();
  MoFEMFunctionBegin;

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
    if (!seed_ents_range.empty())
      CHKERR SetBitRefLevelTool(m_field, bit, ref_ents_ptr, ref_fe_ptr)
          .addEntsToDatabase(seed_ents_range);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::setFieldEntitiesBitRefLevel(
    const std::string field_name, const BitRefLevel bit, int verb) const {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;
  EntityHandle field_meshset = m_field.get_field_meshset(field_name);
  Range field_ents;
  CHKERR m_field.get_moab().get_entities_by_handle(field_meshset, field_ents,
                                                   true);
  CHKERR setEntitiesBitRefLevel(field_ents, bit, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::addToDatabaseBitRefLevelByType(
    const EntityType type, const BitRefLevel bit, const BitRefLevel mask,
    int verb) const {
  MoFEMFunctionBegin;
  Range ents;
  CHKERR getEntitiesByTypeAndRefLevel(bit, mask, type, ents);
  CHKERR setBitRefLevel(ents, BitRefLevel(), false);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::addToDatabaseBitRefLevelByDim(
    const int dim, const BitRefLevel bit, const BitRefLevel mask,
    int verb) const {
  MoFEMFunctionBegin;
  Range ents;
  CHKERR getEntitiesByDimAndRefLevel(bit, mask, dim, ents);
  CHKERR setBitRefLevel(ents, BitRefLevel(), false);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::setBitLevelToMeshset(const EntityHandle meshset,
                                                   const BitRefLevel bit,
                                                   int verb) const {
  MoFEM::Interface &m_field = cOre;
  auto ref_ents_ptr = m_field.get_ref_ents();
  auto ref_fe_ptr = m_field.get_ref_finite_elements();
  MoFEMFunctionBegin;
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

  MOFEM_LOG_FUNCTION();
  MOFEM_LOG("BitRefSelf", Sev::noisy)
      << "Add meshset as ref_ent " << **p_fe.first;

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
    MOFEM_LOG_FUNCTION();
    MOFEM_LOG("BitRefSelf", Sev::noisy) << "Add add bit ref level by dim ";
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
  auto ref_ents_ptr = m_field.get_ref_ents();
  MoFEMFunctionBeginHot;
  RefEntity_multiIndex::iterator dit, hi_dit;
  dit = ref_ents_ptr->begin();
  hi_dit = ref_ents_ptr->end();
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
  auto ref_ents_ptr = m_field.get_ref_ents();
  MoFEMFunctionBegin;
  RefEntity_change_right_shift right_shift(1, mask);
  for (int ii = 0; ii < shift; ii++) {
    // delete bits on the right which are shifted to zero
    BitRefLevel delete_bits = BitRefLevel().set(0) & mask;
    if (delete_bits.any()) {
      CHKERR m_field.delete_ents_by_bit_ref(delete_bits, delete_bits, true,
                                            verb);
    }
    for (RefEntity_multiIndex::iterator ent_it = ref_ents_ptr->begin();
         ent_it != ref_ents_ptr->end(); ent_it++) {
      if (verb > NOISY) {
        MOFEM_LOG_FUNCTION();
        MOFEM_LOG("BitRefSelf", Sev::noisy)
            << (*ent_it)->getBitRefLevel() << " : ";
      }
      right_shift(const_cast<boost::shared_ptr<RefEntity> &>(*ent_it));
      if (verb >= VERY_NOISY) {
        MOFEM_LOG_FUNCTION();
        MOFEM_LOG("BitRefSelf", Sev::noisy) << (*ent_it)->getBitRefLevel();
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
  auto ref_ents_ptr = m_field.get_ref_ents();
  MoFEMFunctionBegin;
  auto &ref_ents = ref_ents_ptr->get<Ent_Ent_mi_tag>();
  auto it = ref_ents.lower_bound(get_id_for_min_type(type));
  auto hi_it = ref_ents.upper_bound(get_id_for_max_type(type));
  std::vector<EntityHandle> ents_vec;
  ents_vec.reserve(std::distance(it, hi_it));
  for (; it != hi_it; it++) {
    const BitRefLevel &ent_bit = it->get()->getBitRefLevel();
    if ((ent_bit & mask) == ent_bit && (ent_bit & bit).any())
      ents_vec.emplace_back(it->get()->getEnt());
  }
  ents.insert_list(ents_vec.begin(), ents_vec.end());
  if (verb > NOISY)
    cerr << "getEntitiesByParentType: " << ents << endl;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::getAllEntitiesNotInDatabase(Range &ents) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  CHKERR moab.get_entities_by_handle(0, ents, false);
  ents = subtract(ents, ents.subset_by_type(MBENTITYSET));
  CHKERR filterEntitiesNotInDatabase(ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitRefManager::filterEntitiesNotInDatabase(Range &ents) const {
  MoFEM::Interface &m_field = cOre;
  auto ref_ents_ptr = m_field.get_ref_ents();
  MoFEMFunctionBeginHot;
  auto eit = ents.begin();
  for (; eit != ents.end();) {
    auto rit = ref_ents_ptr->get<Ent_mi_tag>().find(*eit);
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

  if (verb > VERBOSE) {
    MOFEM_LOG_FUNCTION();
    MOFEM_LOG("BitRef", Sev::noisy) << "from: " << bit;
  }

  CHKERR moab.get_adjacencies(from_entities, num_entities, to_dimension, false,
                              adj_entities, operation_type);
  std::vector<BitRefLevel> bit_levels(adj_entities.size());
  CHKERR moab.tag_get_data(cOre.get_th_RefBitLevel(), adj_entities,
                           &*bit_levels.begin());
  std::vector<BitRefLevel>::iterator b_it = bit_levels.begin();
  for (Range::iterator eit = adj_entities.begin(); eit != adj_entities.end();
       b_it++) {

    if (verb > VERBOSE) {
      RefEntity adj_entity(m_field.get_basic_entity_data_ptr(), *eit);
      MOFEM_LOG_FUNCTION();
      MOFEM_LOG("BitRef", Sev::noisy)
          << "to: " << adj_entity.getBitRefLevel() << " : " << adj_entity;
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
  MoFEMFunctionBegin;
  Range parent_ents;
  CHKERR moab.get_entities_by_handle(parent, parent_ents, recursive);
  CHKERR filterEntitiesByRefLevel(parent_bit, parent_mask, parent_ents, verb);
  if (verb >= VERY_VERBOSE) {
    MOFEM_LOG_FUNCTION();
    MOFEM_LOG("BitRefSelf", Sev::noisy) << "Parnets: " << parent;
  }
  Range children_ents;
  CHKERR updateRange(parent_ents, children_ents);
  CHKERR filterEntitiesByRefLevel(child_bit, BitRefLevel().set(), children_ents,
                                  verb);
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
  auto fields_ptr = m_field.get_fields();
  auto ref_ents_ptr = m_field.get_ref_ents();
  MoFEMFunctionBegin;

  for (auto &fit : (*fields_ptr)) {

    const EntityHandle meshset = fit->getMeshset();
    Range parent_ents;
    CHKERR moab.get_entities_by_handle(meshset, parent_ents, true);

    if (verb >= VERY_VERBOSE) {
      MOFEM_LOG_FUNCTION();
      MOFEM_LOG("BitRefSelf", Sev::noisy) << "Parnets:" << std::endl
                                          << parent_ents << std::endl;
    }

    Range children_ents;
    CHKERR updateRange(parent_ents, children_ents);
    CHKERR filterEntitiesByRefLevel(child_bit, BitRefLevel().set(),
                                    children_ents, verb);

    CHKERR moab.add_entities(meshset, children_ents);

    if (verb >= VERY_VERBOSE) {
      MOFEM_LOG_FUNCTION();
      MOFEM_LOG("BitRefSelf", Sev::noisy) << "Children: " << std::endl
                                          << children_ents << std::endl;
    }
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

  if (verb >= VERBOSE) {
    MOFEM_LOG_FUNCTION();
    MOFEM_LOG("BitRefSelf", Sev::noisy) << "Parnets:" << endl
                                        << parent_ents << std::endl;
  }

  Range children_ents;
  CHKERR updateRange(parent_ents, children_ents);
  CHKERR filterEntitiesByRefLevel(child_bit, BitRefLevel().set(), children_ents,
                                  verb);

  CHKERR moab.add_entities(field_meshset, children_ents);

  if (verb >= VERBOSE) {
    MOFEM_LOG_FUNCTION();
    MOFEM_LOG("BitRefSelf", Sev::noisy) << "Children: " << endl
                                        << children_ents << std::endl;
  }

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
                                          Range &child_ents, MoFEMTypes bh) {
  MoFEM::Interface &m_field = cOre;
  auto ref_ents_ptr = m_field.get_ref_ents();
  MoFEMFunctionBegin;
  auto &ref_ents =
      const_cast<RefEntity_multiIndex *>(ref_ents_ptr)->get<Ent_Ent_mi_tag>();
  std::vector<EntityHandle> child_ents_vec;
  child_ents_vec.reserve(ref_ents.size());
  for (Range::const_pair_iterator pit = parent_ents.const_pair_begin();
       pit != parent_ents.const_pair_end(); pit++) {
    const auto f = pit->first;
    auto it = ref_ents.lower_bound(f);
    if (it != ref_ents.end()) {
      const auto s = pit->second;
      auto hi_it = ref_ents.upper_bound(s);
      if (bh == MF_EXIST) {
        if (std::distance(it, hi_it) != (s - f) + 1) {
          SETERRQ2(
              PETSC_COMM_SELF, MOFEM_NOT_FOUND,
              "Number of entities and enties parents is diffrent %d != %d ",
              std::distance(it, hi_it), (s - f) + 1);
        }
      }
      for (; it != hi_it; it++) {
        if (it->get()->getEntType() == MBENTITYSET) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
                  "This should not happen; Entity should not have part of the "
                  "meshset. It has no children.");
        }
        child_ents_vec.emplace_back((*it)->getEnt());
      }
    } else if (bh == MF_EXIST) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "Entities not found");
    }
  }
  child_ents.insert_list(child_ents_vec.begin(), child_ents_vec.end());
  MoFEMFunctionReturn(0);
}
} // namespace MoFEM
