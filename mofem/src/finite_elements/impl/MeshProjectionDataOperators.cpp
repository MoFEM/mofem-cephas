/** \file MeshProjectionDataOperators.cpp


*/

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

namespace MoFEM {

OpRunParent::OpRunParent(boost::shared_ptr<ForcesAndSourcesCore> parent_ele_ptr,
                         BitRefLevel bit_parent, BitRefLevel bit_parent_mask,
                         boost::shared_ptr<ForcesAndSourcesCore> this_ele_ptr,
                         BitRefLevel bit_this, BitRefLevel bit_this_mask,
                         int verb, Sev sev)
    : ForcesAndSourcesCore::UserDataOperator(
          NOSPACE, ForcesAndSourcesCore::UserDataOperator::OPSPACE),
      parentElePtr(parent_ele_ptr), bitParent(bit_parent),
      bitParentMask(bit_parent_mask), thisElePtr(this_ele_ptr),
      bitThis(bit_this), bitThisMask(bit_this_mask), verbosity(verb),
      severityLevel(sev) {}

MoFEMErrorCode OpRunParent::doWork(int side, EntityType type,
                                   EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  auto &bit = getFEMethod()->numeredEntFiniteElementPtr->getBitRefLevel();

  auto check = [&](auto &b, auto &m) {
    return ((bit & b).any()) && ((bit & m) == bit);
  };

  if (check(bitParent, bitParentMask)) {

    CHKERR loopParent(getFEName(), parentElePtr.get(), verbosity,
                      severityLevel);

  } else if (check(bitThis, bitThisMask)) {

    CHKERR loopThis(getFEName(), thisElePtr.get(), verbosity, severityLevel);
  }

  MoFEMFunctionReturn(0);
}

OpAddParentEntData::OpAddParentEntData(
    std::string field_name, OpType op_parent_type,
    boost::shared_ptr<ForcesAndSourcesCore> parent_ele_ptr,
    BitRefLevel bit_child, BitRefLevel bit_child_mask,
    BitRefLevel bit_parent_ent, BitRefLevel bit_parent_ent_mask, int verb,
    Sev sev)
    : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPSPACE),
      fieldName(field_name), opParentType(op_parent_type),
      parentElePtr(parent_ele_ptr), bitChild(bit_child),
      bitChildMask(bit_child_mask), bitParentEnt(bit_parent_ent),
      bitParentEntMask(bit_parent_ent_mask), verbosity(verb),
      severityLevel(sev) {}

MoFEMErrorCode OpAddParentEntData::doWork(int side, EntityType type,
                                          EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  auto check = [](auto &b, auto &m, auto &bit) {
    return ((bit & b).any()) && ((bit & m) == bit);
  };

  auto get_entities_field_data_ptr = [&](auto space) {
    switch (opParentType) {
    case OPROW:
      return getPtrFE()->getDataOnElementBySpaceArray()[space];
    case OPCOL:
      return getPtrFE()->getDerivedDataOnElement()[space];
    default:
      return boost::shared_ptr<EntitiesFieldData>();
    }
  };

  auto &bit_fe = getFEMethod()->numeredEntFiniteElementPtr->getBitRefLevel();
  if (check(bitChild, bitChildMask, bit_fe)) {

    MOFEM_LOG("SELF", severityLevel) << "Child FE bit: " << bit_fe;

    auto loop_parent_fe = [&]() {
      MoFEMFunctionBeginHot;

#ifndef NDEBUG
      if (verbosity >= VERBOSE) {
        MOFEM_LOG("SELF", severityLevel)
            << "Loop parent element in OpAddParentEntData";
      }
#endif

      // note live of op pointer is controlled by ptr_vec in in finite
      // element
      auto field_op =
          new ForcesAndSourcesCore::UserDataOperator(fieldName, opParentType);
      // that forces to run operator at last instance and collect data on
      // entities
      field_op->doWorkRhsHook = [&](DataOperator *op_ptr, int side,
                                    EntityType type,
                                    EntitiesFieldData::EntData &data) {
        MoFEMFunctionBegin;

        auto field_entities = data.getFieldEntities();
        if (field_entities.size()) {

          if (field_entities.size() == 1) {
            auto &bit_ent = field_entities[0]->getBitRefLevel();
            if (!check(bitParentEnt, bitParentEntMask, bit_ent))
              MoFEMFunctionReturnHot(0);
          }

#ifndef NDEBUG
          if (verbosity >= VERBOSE) {
            for (auto field_ent : field_entities) {
              MOFEM_LOG("SELF", severityLevel)
                  << "Parent FE bit: " << field_ent->getBitRefLevel() << " "
                  << *field_ent;
            }
          }
#endif

          auto space = data.getSpace();
          auto base = data.getBase();

          if (verbosity == VERBOSE) {
            MOFEM_LOG("SELF", severityLevel)
                << "Side/type: " << side << "/" << CN::EntityTypeName(type)
                << " op space/base: " << FieldSpaceNames[space] << "/"
                << ApproximationBaseNames[space];
          }

          auto entities_field_data_ptr = get_entities_field_data_ptr(space);
          if (!entities_field_data_ptr)
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "Only OPROW/OPCOL op type is allows");

          entities_field_data_ptr->dataOnEntities[MBENTITYSET].push_back(
              new EntitiesFieldData::EntData());
          auto &child_data =
              entities_field_data_ptr->dataOnEntities[MBENTITYSET].back();

          child_data.sPace = space;
          child_data.bAse = base;

          child_data.sEnse = data.getSense();
          child_data.oRder = data.getOrder();
          child_data.iNdices = data.getIndices();
          child_data.localIndices = data.getLocalIndices();
          child_data.dOfs = data.getFieldDofs();
          child_data.fieldEntities = field_entities;
          child_data.fieldData = data.getFieldData();

          int direvative = 0;
          for (auto &b : data.baseFunctionsAndBaseDerivatives) {
            if (b[base]) {
              child_data.baseFunctionsAndBaseDerivatives[direvative] = b;
            }
            ++direvative;
          }

          if (field_entities.size() > 1) {
            int dof = 0;
            for (auto &field_ent : field_entities) {
              auto &bit_ent = field_ent->getBitRefLevel();
              if (!check(bitParentEnt, bitParentEntMask, bit_ent)) {
                for (auto d = 0; d != field_ent->getNbDofsOnEnt(); ++d) {
                  child_data.iNdices[dof + d] = -1;
                  child_data.localIndices[dof + d] = -1;
                  child_data.fieldData[dof + d] = 0;
                }
              }
            }
          }
        }

        MoFEMFunctionReturn(0);
      };

      parentElePtr->getOpPtrVector().push_back(field_op);
      CHKERR loopParent(getFEName(), parentElePtr.get(), verbosity,
                        severityLevel);
                      
 #ifndef NDEBUG
      auto &parent_gauss_pts = parentElePtr->gaussPts;
      if (getGaussPts().size1() != parent_gauss_pts.size1()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Wrong number of weights");
      }
      if (getGaussPts().size2() != parent_gauss_pts.size2()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Wrong number of itegartion points");
      }
#endif

      // clean element from obsolete data operator
      parentElePtr->getOpPtrVector().pop_back();

      MoFEMFunctionReturnHot(0);
    };

    CHKERR loop_parent_fe();
  }

  MoFEMFunctionReturn(0);
}

OpRestoreEntData::OpRestoreEntData(FieldSpace space, OpType op_type, int verb,
                                   Sev sev)
    : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPSPACE), sPace(space),
      opType(op_type), verbosity(verb), severityLevel(sev) {}

MoFEMErrorCode OpRestoreEntData::doWork(int side, EntityType type,
                                        EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

#ifndef NDEBUG
  auto log_meshset_data = [&](auto &data_on_meshset) {
    MoFEMFunctionBegin;
    MOFEM_LOG("SELF", severityLevel)
        << "Nb meshset data " << data_on_meshset.size();
    int side = 0;
    for (auto &ent_data : data_on_meshset) {
      for (auto field_ent : ent_data.getFieldEntities()) {
        MOFEM_LOG("SELF", severityLevel)
            << "Side " << side << ": " << *field_ent;
      }
      ++side;
    }
    MoFEMFunctionReturn(0);
  };
#endif

  switch (opType) {
  case OPROW:
#ifndef NDEBUG
    CHKERR log_meshset_data(
        getPtrFE()->getDataOnElementBySpaceArray()[sPace]->dataOnEntities[MBENTITYSET]);
#endif
    getPtrFE()->getDataOnElementBySpaceArray()[sPace]->dataOnEntities[MBENTITYSET].clear();
    break;
  case OPCOL:
#ifndef NDEBUG
    CHKERR log_meshset_data(getPtrFE()
                                ->getDerivedDataOnElement()[sPace]
                                ->dataOnEntities[MBENTITYSET]);
#endif
    getPtrFE()
        ->getDerivedDataOnElement()[sPace]
        ->dataOnEntities[MBENTITYSET]
        .clear();
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Only OPROW/OPCOL op type is allows");
  }

  MoFEMFunctionReturn(0);
}

OpSwitchOffIndices::OpSwitchOffIndices(const std::string field_name,
                                       const OpType op_type,
                                       const BitRefLevel bit_ent,
                                       const BitRefLevel bit_mask_ent)
    : ForcesAndSourcesCore::UserDataOperator(field_name, op_type),
      bitEnt(bit_ent), bitEntMask(bit_mask_ent) {}

MoFEMErrorCode OpSwitchOffIndices::doWork(int side, EntityType type,
                                          EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  auto check = [&](auto &bit) {
    return ((bit & bitEnt).any()) && ((bit & bitEntMask) == bit);
  };

  size_t dd = 0;
  for (auto &dof : data.getFieldDofs()) {
    auto &bit = dof->getBitRefLevel();
    if (check(bit)) {
      data.getIndices()[dd] = -1;
      data.getLocalIndices()[dd] = -1;
      data.getFieldData()[dd] = 0;
    }
    ++dd;
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpAddParentEntData::markHangingSkin(
    MoFEM::Interface &m_field, const int dim, const BitRefLevel parent_bit,
    const BitRefLevel parent_mask, const BitRefLevel child_bit,
    const BitRefLevel child_mask, const BitRefLevel mark_bit,
    const bool resolve_shared, const std::string debug_file_name) {
  BitRefManager *bit_mng = m_field.getInterface<BitRefManager>();
  ParallelComm *pcomm =
      ParallelComm::get_pcomm(&m_field.get_moab(), MYPCOMM_INDEX);
  Skinner skin(&m_field.get_moab());
  MoFEMFunctionBegin;

  Range parent_level, child_only, child_level;
  CHKERR bit_mng->getEntitiesByDimAndRefLevel(parent_bit, parent_mask, dim,
                                              parent_level);
  CHKERR bit_mng->getEntitiesByDimAndRefLevel(child_bit, child_bit, dim,
                                              child_only);
  CHKERR bit_mng->getEntitiesByDimAndRefLevel(child_bit, child_mask, dim,
                                              child_level);
  auto parnets_on_child = subtract(child_level, child_only);

  Range parent_skin;
  CHKERR skin.find_skin(0, parent_level, false, parent_skin);
  if (resolve_shared)
    CHKERR pcomm->filter_pstatus(parent_skin,
                                 PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                 PSTATUS_NOT, -1, nullptr);
  Range child_skin;
  CHKERR skin.find_skin(0, child_level, false, child_skin);
  if (resolve_shared)
    CHKERR pcomm->filter_pstatus(child_skin,
                                 PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                 PSTATUS_NOT, -1, nullptr);

  Range parents_on_child_skin;
  CHKERR skin.find_skin(0, parnets_on_child, false, parents_on_child_skin);
  if (resolve_shared)
    CHKERR pcomm->filter_pstatus(parents_on_child_skin,
                                 PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                 PSTATUS_NOT, -1, nullptr);
  parents_on_child_skin = subtract(parents_on_child_skin, parent_skin);

  for (auto d = dim - 1; d >= 0; --d) {
    Range adj;
    CHKERR m_field.get_moab().get_adjacencies(parents_on_child_skin, d, false,
                                              adj, moab::Interface::UNION);
    parents_on_child_skin.merge(adj);
  }

  CHKERR bit_mng->addBitRefLevel(parents_on_child_skin, mark_bit);

  if (!debug_file_name.empty()) {
    auto proc_rank = m_field.get_comm_rank();
    auto name =
        boost::lexical_cast<std::string>(proc_rank) + "_" + debug_file_name;
    MOFEM_LOG("SELF", Sev::verbose) << "Save debug file: " << name;
    CHKERR bit_mng->writeBitLevel(mark_bit, parent_mask | mark_bit,
                                  name.c_str(), "VTK", "");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpAddParentEntData::markHangingSkinChildren(
    MoFEM::Interface &m_field, const BitRefLevel child_bit,
    const BitRefLevel child_mask, const BitRefLevel mark_bit,
    const BitRefLevel mark_mask, const std::string debug_file_name) {
  BitRefManager *bit_mng = m_field.getInterface<BitRefManager>();
  auto refined_ents_ptr = m_field.get_ref_ents();
  MoFEMFunctionBegin;

  Range markers;
  CHKERR bit_mng->getEntitiesByRefLevel(mark_bit, mark_mask, markers);

  std::vector<EntityHandle> handles;

  auto get_children = [&]() {
    Range children;

    for (auto p = markers.pair_begin(); p != markers.pair_end(); ++p) {
      const auto f = p->first;
      const auto s = p->second;

      auto &ref_ents_by_parents = refined_ents_ptr->get<Ent_Ent_mi_tag>();
      auto lo = ref_ents_by_parents.lower_bound(f);
      auto hi = ref_ents_by_parents.upper_bound(s);

      handles.resize(std::distance(lo, hi));
      for (; lo != hi; ++lo) {
        handles.push_back((*lo)->getEnt());
      }
      children.insert_list(handles.begin(), handles.end());

    }
 
    CHKERR bit_mng->filterEntitiesByRefLevel(child_bit, child_mask, children);

    return children;
  };

  auto children = get_children();

  CHKERR bit_mng->addBitRefLevel(children, mark_bit);

  if (!debug_file_name.empty()) {
    auto proc_rank = m_field.get_comm_rank();
    auto name =
        boost::lexical_cast<std::string>(proc_rank) + "_" + debug_file_name;
    MOFEM_LOG("SELF", Sev::verbose) << "Save debug file: " << name;
    CHKERR bit_mng->writeBitLevel(mark_bit, BitRefLevel().set() /*mark_bit | child_mask*/, name.c_str(),
                                  "VTK", "");
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM