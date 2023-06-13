/** \file MeshProjectionDataOperators.cpp


*/

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

  if (verbosity > QUIET) {
    MOFEM_LOG_CHANNEL("SELF");
    MOFEM_TAG_AND_LOG("SELF", severityLevel, "OpRunParent")
        << "FE bit " << bit
        << " check parent = " << check(bitParent, bitParentMask)
        << " check this " << check(bitThis, bitThisMask);
  }

  if (check(bitParent, bitParentMask)) {
    if (parentElePtr)
      CHKERR loopParent(getFEName(), parentElePtr.get(), verbosity,
                        severityLevel);

  } else if (check(bitThis, bitThisMask)) {
    if (thisElePtr)
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
    : ForcesAndSourcesCore::UserDataOperator(field_name, op_parent_type),
      fieldName(field_name), opParentType(op_parent_type),
      parentElePtr(parent_ele_ptr), bitChild(bit_child),
      bitChildMask(bit_child_mask), bitParentEnt(bit_parent_ent),
      bitParentEntMask(bit_parent_ent_mask), verbosity(verb),
      severityLevel(sev) {
  // Push op to collect data
  auto field_op =
      new ForcesAndSourcesCore::UserDataOperator(fieldName, opParentType);
  parentElePtr->getOpPtrVector().push_back(field_op);
}

OpAddParentEntData::OpAddParentEntData(
    FieldSpace space, OpType op_parent_type,
    boost::shared_ptr<ForcesAndSourcesCore> parent_ele_ptr,
    BitRefLevel bit_child, BitRefLevel bit_child_mask,
    BitRefLevel bit_parent_ent, BitRefLevel bit_parent_ent_mask, int verb,
    Sev sev)
    : ForcesAndSourcesCore::UserDataOperator(space, op_parent_type),
      fieldName(""), approxSpace(space), opParentType(op_parent_type),
      parentElePtr(parent_ele_ptr), bitChild(bit_child),
      bitChildMask(bit_child_mask), bitParentEnt(bit_parent_ent),
      bitParentEntMask(bit_parent_ent_mask), verbosity(verb),
      severityLevel(sev) {
  // Push op to collect data
  auto field_op =
      new ForcesAndSourcesCore::UserDataOperator(approxSpace, opParentType);
  parentElePtr->getOpPtrVector().push_back(field_op);
}

MoFEMErrorCode OpAddParentEntData::opRhs(EntitiesFieldData &entities_field_data,
                                         const bool error_if_no_base) {
  int count_meshset_sides = 0; // count number of data on parent element
  MoFEMFunctionBegin;

  auto check = [](auto &b, auto &m, auto &bit) {
    return ((bit & b).any()) && ((bit & m) == bit);
  };

  auto set_child_data_entity = [](auto &parent_data, auto &child_data) {
    MoFEMFunctionBeginHot;
    child_data.getEntDataBitRefLevel() = parent_data.getEntDataBitRefLevel();
    child_data.sPace = parent_data.getSpace();
    child_data.bAse = parent_data.getBase();
    child_data.sEnse = parent_data.getSense();
    child_data.oRder = parent_data.getOrder();
    child_data.iNdices.swap(parent_data.getIndices());
    child_data.localIndices.swap(parent_data.getLocalIndices());
    child_data.dOfs.swap(parent_data.getFieldDofs());
    child_data.fieldData.swap(parent_data.getFieldData());
    child_data.fieldEntities.swap(parent_data.getFieldEntities());
    MoFEMFunctionReturnHot(0);
  };

  auto set_child_data_vertex = [&](auto &parent_data, auto &child_data,
                                   int node, int num_nodes) {
    MoFEMFunctionBegin;

    child_data.getEntDataBitRefLevel().resize(1, false);
    child_data.getEntDataBitRefLevel()[0] =
        parent_data.getEntDataBitRefLevel()[node];
    child_data.sPace = parent_data.getSpace();
    child_data.bAse = parent_data.getBase();
    child_data.sEnse = parent_data.getSense();
    child_data.oRder = parent_data.getOrder();

    if (parent_data.dOfs.size()) {

      auto &field_entities = parent_data.getFieldEntities();
      auto &vertex_entity = field_entities[node];

      const auto nb_coeffs = vertex_entity->getNbOfCoeffs();

      child_data.fieldEntities.resize(1, false);
      child_data.fieldEntities[0] = vertex_entity;
#ifndef NDEBUG
      if (parent_data.dOfs.size() != nb_coeffs * num_nodes)
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Inconsistent number of DOFs and vertices %d != %d",
                 parent_data.dOfs.size(), nb_coeffs * num_nodes);
#endif

      // It could be a case that all DOFs on element are removed from the
      // problem, the size of indices vector is zero.
      if (parent_data.iNdices.size()) {

#ifndef NDEBUG
        if (parent_data.dOfs.size() != parent_data.iNdices.size())
          SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "Inconsistent size of DOFs %d != %d",
                   parent_data.dOfs.size(), parent_data.iNdices.size());
#endif

        child_data.iNdices.resize(nb_coeffs, false);
        child_data.localIndices.resize(nb_coeffs, false);
        child_data.dOfs.resize(nb_coeffs, false);
        child_data.fieldData.resize(nb_coeffs, false);

        int DD = 0;
        for (auto dd = node * nb_coeffs; dd != (node + 1) * nb_coeffs;
             ++dd, ++DD) {
          child_data.iNdices[DD] = parent_data.iNdices[dd];
          child_data.localIndices[DD] = parent_data.localIndices[dd];
          child_data.dOfs[DD] = parent_data.dOfs[dd];
          child_data.fieldData[DD] = parent_data.fieldData[dd];
        }
      } else {

        child_data.iNdices.clear();
        child_data.localIndices.clear();
        child_data.dOfs.resize(nb_coeffs, false);
        child_data.fieldData.resize(nb_coeffs, false);
        int DD = 0;
        for (auto dd = node * nb_coeffs; dd != (node + 1) * nb_coeffs;
             ++dd, ++DD) {
          child_data.dOfs[DD] = parent_data.dOfs[dd];
          child_data.fieldData[DD] = parent_data.fieldData[dd];
        }
        
      }
    }

    MoFEMFunctionReturn(0);
  };

  /**
   * @brief swap child base
   *
   * @todo add swap for Bernstein-Bezier base
   *
   */
  auto set_child_base_entity = [](auto &parent_data, auto &child_data) {
    MoFEMFunctionBeginHot;
    child_data.resetFieldDependentData();
    child_data.getEntDataBitRefLevel() = parent_data.getEntDataBitRefLevel();
    child_data.sPace = parent_data.getSpace();

#ifndef NDEBUG
    if (child_data.bAse == AINSWORTH_BERNSTEIN_BEZIER_BASE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Base not implemented");
#endif

    using BaseDerivatives = EntitiesFieldData::EntData::BaseDerivatives;
    for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
      for (auto derivative = 0; derivative != BaseDerivatives::LastDerivative;
           ++derivative) {
        auto parent_base =
            parent_data.getNSharedPtr(static_cast<FieldApproximationBase>(b),
                                      static_cast<BaseDerivatives>(derivative));
        if (parent_base) {
          auto &child_base = child_data.getNSharedPtr(
              static_cast<FieldApproximationBase>(b),
              static_cast<BaseDerivatives>(derivative));
          if (!child_base)
            child_base = boost::make_shared<MatrixDouble>();
          child_base->swap(*parent_base);
        }
      }
    }

    MoFEMFunctionReturnHot(0);
  };

  auto set_child_base_vertex = [](auto &parent_data, auto &child_data, int node,
                                  int num_nodes) {
    MoFEMFunctionBeginHot;
    child_data.resetFieldDependentData();
    child_data.getEntDataBitRefLevel().resize(1, false);
    child_data.getEntDataBitRefLevel()[0] =
        parent_data.getEntDataBitRefLevel()[node];
    child_data.sPace = parent_data.getSpace();

#ifndef NDEBUG
    if (child_data.bAse == AINSWORTH_BERNSTEIN_BEZIER_BASE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Base not implemented");
#endif

    using BaseDerivatives = EntitiesFieldData::EntData::BaseDerivatives;
    for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
      for (auto derivative = 0; derivative != BaseDerivatives::LastDerivative;
           ++derivative) {
        auto parent_base =
            parent_data.getNSharedPtr(static_cast<FieldApproximationBase>(b),
                                      static_cast<BaseDerivatives>(derivative));
        if (parent_base) {
          auto &child_base = child_data.getNSharedPtr(
              static_cast<FieldApproximationBase>(b),
              static_cast<BaseDerivatives>(derivative));

          const auto num_bases_per_node = parent_base->size2() / num_nodes;
          if (parent_base->size2() % num_nodes) {
            SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                     "Inconsistent nb of base functions and nodes mod(%d, %d)",
                     parent_base->size2(), num_nodes);
          }

          if (!child_base)
            child_base = boost::make_shared<MatrixDouble>(parent_base->size1(),
                                                          num_bases_per_node);
          else
            child_base->resize(parent_base->size1(), num_bases_per_node, false);

          for (auto gg = 0; gg != parent_base->size1(); ++gg) {
            int DD = 0;
            for (auto dd = node * num_bases_per_node;
                 dd != (node + 1) * num_bases_per_node; ++dd, ++DD) {
              (*child_base)(gg, DD) = (*parent_base)(gg, dd);
            }
          }
        }
      }
    }

    MoFEMFunctionReturnHot(0);
  };

  /**
   * @brief if DOF is on children set -1 to its repetition on parent
   *
   */
  auto switch_off_dofs_children = [&](auto &parent_ent_data, auto &child_data) {
    MoFEMFunctionBeginHot;
    for (auto i : parent_ent_data.getIndices()) {
      if (i >= 0) {
        for (auto &child_ent_data : child_data) {
          auto it = std::find(child_ent_data.getIndices().begin(),
                              child_ent_data.getIndices().end(), i);
          if (it != child_ent_data.getIndices().end()) {
            const auto dof_idx =
                std::distance(child_ent_data.getIndices().begin(), it);
            auto &bit_dof =
                child_ent_data.getFieldDofs()[dof_idx]->getBitRefLevel();
            if (check(bitParentEnt, bitParentEntMask, bit_dof)) {
              child_ent_data.getIndices()[dof_idx] = -1;
              child_ent_data.getLocalIndices()[dof_idx] = -1;
              child_ent_data.getFieldData()[dof_idx] = 0;
            }
          }
        }
      }
    }
    MoFEMFunctionReturnHot(0);
  };

  // that forces to run operator at last instance and collect data on
  // entities
  auto do_work_parent_hook = [&](DataOperator *op_ptr, int side,
                                 EntityType type,
                                 EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;

    auto up_op_ptr = static_cast<UserDataOperator *>(op_ptr);
    auto &field_entities = data.getFieldEntities();

#ifndef NDEBUG
    if (verbosity >= VERBOSE) {
      MOFEM_LOG("SELF", severityLevel)
          << "Add entity data to meshset "
          << "side/type: " << side << "/" << CN::EntityTypeName(type)
          << " op space/base: " << FieldSpaceNames[data.getSpace()] << "/"
          << ApproximationBaseNames[data.getBase()]
          << " OPSPACE: " << ((opParentType == OPSPACE) ? "Yes" : "No");
    }
    if (verbosity >= NOISY) {
      for (auto field_ent : field_entities) {
        MOFEM_LOG("SELF", severityLevel)
            << "Parent field entity bit: " << field_ent->getBitRefLevel() << " "
            << *field_ent;
      }
    }
#endif

    const auto nb_ents = field_entities.size();
    const auto nb_nodes = data.getEntDataBitRefLevel().size();

    // #ifndef NDEBUG
    if (type == MBVERTEX &&
        nb_nodes != up_op_ptr->getNumberOfNodesOnElement()) {
      SETERRQ2(
          PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
          "Inconsistency between bit ref levels and number of nodes %d != %d",
          nb_nodes, up_op_ptr->getNumberOfNodesOnElement());
    }
    // #endif

    auto get_child_meshset =
        [&](auto count_meshset_sides) -> MoFEM::EntitiesFieldData::EntData & {
      auto &data_on_meshset = entities_field_data.dataOnEntities[MBENTITYSET];
      if (data_on_meshset.size() < count_meshset_sides) {
        if (poolEntsVector.size()) {
          // transfer data from pool to data_on_meshset
          data_on_meshset.transfer(data_on_meshset.end(),
                                   poolEntsVector.begin(), poolEntsVector);
        } else {
          // resize, no data on pool
          entities_field_data.dataOnEntities[MBENTITYSET].resize(
              count_meshset_sides);
        }
      }
      return entities_field_data
          .dataOnEntities[MBENTITYSET][count_meshset_sides - 1];
    };

    if (opParentType == OPSPACE) {
      for (auto node = 0; node != nb_nodes; ++node) {
        auto &bit_ent = data.getEntDataBitRefLevel()[node];
        if (check(bitParentEnt, bitParentEntMask, bit_ent)) {
          ++count_meshset_sides;
          auto &child_data_meshset = get_child_meshset(count_meshset_sides);
          if (type == MBVERTEX) {
            const auto num_nodes = up_op_ptr->getNumberOfNodesOnElement();
            CHKERR set_child_base_vertex(data, child_data_meshset, node,
                                         num_nodes);
          } else {
            CHKERR set_child_base_entity(data, child_data_meshset);
          }
        }
      }

    } else {

      {
        boost::container::static_vector<EntityType, 128> ents_type;
        ents_type.reserve(field_entities.size());
        for (auto &field_entity : field_entities)
          ents_type.push_back(field_entity->getEntType());
        std::sort(ents_type.begin(), ents_type.end());
        auto end = std::unique(ents_type.begin(), ents_type.end());
        for (auto it_t = ents_type.begin(); it_t != end; ++it_t)
          CHKERR switch_off_dofs_children(
              data, entities_field_data.dataOnEntities[*it_t]);
        if (type == MBENTITYSET)
          CHKERR switch_off_dofs_children(
              data, entities_field_data.dataOnEntities[type]);
      }

      if (type == MBVERTEX &&
          nb_ents != up_op_ptr->getNumberOfNodesOnElement()) {
        SETERRQ2(
            PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Number of nodes is expected to much number of entities: %d != %d",
            nb_ents, up_op_ptr->getNumberOfNodesOnElement());
      } else if (type != MBVERTEX && nb_ents > 1) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Code can handle only one entity");
      }

      for (auto node = 0; node != nb_nodes; ++node) {
        auto &bit_ent = data.getEntDataBitRefLevel()[node];
        if (check(bitParentEnt, bitParentEntMask, bit_ent)) {
          ++count_meshset_sides;
          auto &child_data_meshset = get_child_meshset(count_meshset_sides);
          if (type == MBVERTEX) {
            const auto num_nodes = up_op_ptr->getNumberOfNodesOnElement();
            CHKERR set_child_data_vertex(data, child_data_meshset, node,
                                         num_nodes);
          } else {
            CHKERR set_child_data_entity(data, child_data_meshset);
          }
        }
      }
    }

    MoFEMFunctionReturn(0);
  };

  // iterate parents collect data
  auto &bit_fe = getFEMethod()->numeredEntFiniteElementPtr->getBitRefLevel();
  if (check(bitChild, bitChildMask, bit_fe)) { // check if FE is on right bit

    if (verbosity >= VERBOSE) {
      MOFEM_LOG("SELF", severityLevel) << "Child FE bit: " << bit_fe;
    }

    auto loop_parent_fe = [&]() {
      MoFEMFunctionBeginHot;

      // execute do_work_parent_hook on parent element, and collect data.
      parentElePtr->getOpPtrVector().back().doWorkRhsHook = do_work_parent_hook;
      CHKERR loopParent(getFEName(), parentElePtr.get(), verbosity,
                        severityLevel);

#ifndef NDEBUG
      if (parentElePtr->getLoopSize()) {
        auto &parent_gauss_pts = parentElePtr->gaussPts;
        if (getGaussPts().size1() != parent_gauss_pts.size1()) {
          MOFEM_LOG("SELF", Sev::error)
              << "Calling element: "
              << boost::typeindex::type_id_runtime(*parentElePtr).pretty_name();
          MOFEM_LOG("SELF", Sev::error) << getGaussPts();
          MOFEM_LOG("SELF", Sev::error) << parent_gauss_pts;
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Wrong number of weights");
        }
        if (getGaussPts().size2() != parent_gauss_pts.size2()) {
          MOFEM_LOG("SELF", Sev::error)
              << "Calling element: "
              << boost::typeindex::type_id_runtime(*parentElePtr).pretty_name();
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Wrong number of integration points");
        }
      }
#endif

      MoFEMFunctionReturnHot(0);
    };

    // iterate parent finite lents
    CHKERR loop_parent_fe();
  }

  if (verbosity >= VERBOSE) {
    MOFEM_LOG("SELF", severityLevel)
        << "Number of meshset entities "
        << entities_field_data.dataOnEntities[MBENTITYSET].size();
  }

  auto &data_on_meshset = entities_field_data.dataOnEntities[MBENTITYSET];
  auto it = data_on_meshset.begin();

  // transfer not used data on entities to pool, to be used later
  if (count_meshset_sides < data_on_meshset.size()) {
    for (auto s = 0; s != count_meshset_sides; ++s)
      ++it;
    poolEntsVector.transfer(poolEntsVector.end(), it, data_on_meshset.end(),
                            data_on_meshset);
  }

  auto set_up_derivative_ent_data = [&](auto &entities_field_data,
                                        auto &derivative_entities_field_data) {
    MoFEMFunctionBegin;

    using EntData = EntitiesFieldData::EntData;
    using DerivedEntData = DerivedEntitiesFieldData::DerivedEntData;

    auto &data_ptr = getPtrFE()->getDataOnElementBySpaceArray()[approxSpace];
    auto &ents_data = entities_field_data.dataOnEntities[MBENTITYSET];
    auto &dev_ents_data =
        derivative_entities_field_data.dataOnEntities[MBENTITYSET];
    dev_ents_data.clear();
    for (auto c = 0; c < ents_data.size(); ++c) {
      boost::shared_ptr<EntData> ent_data_ptr(data_ptr, &ents_data[c]);
      dev_ents_data.push_back(new DerivedEntData(ent_data_ptr));
    }
    MoFEMFunctionReturn(0);
  };

  if (opParentType == OPSPACE) {
    CHKERR set_up_derivative_ent_data(
        entities_field_data,
        *(getPtrFE()->getDerivedDataOnElementBySpaceArray()[approxSpace]));
  }

#ifndef NDEBUG
  auto &side_table =
      getPtrFE()->numeredEntFiniteElementPtr->getSideNumberTable();
  for (auto &data : entities_field_data.dataOnEntities[MBENTITYSET]) {
    for (auto dof : data.getFieldDofs()) {
      auto &bit_dof = dof->getBitRefLevel();
      if (check(bitParentEnt, bitParentEntMask, bit_dof)) {
        auto ent = dof->getEnt();
        if (side_table.find(ent) == side_table.end()) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Adjacency not found");
        }
      }
    }
  }
#endif

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM