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

#ifndef NDEBUG
  if (verbosity > QUIET) {
    MOFEM_LOG_CHANNEL("SELF");
    MOFEM_TAG_AND_LOG("SELF", severityLevel, "OpRunParent")
        << "FE bit " << bit
        << " check parent = " << check(bitParent, bitParentMask)
        << " check this " << check(bitThis, bitThisMask);
  }
#endif

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
  int count_meshset_sides = 0;
  MoFEMFunctionBegin;

  auto check = [](auto &b, auto &m, auto &bit) {
    return ((bit & b).any()) && ((bit & m) == bit);
  };

  auto set_child_data = [](auto &parent_data, auto &child_data) {
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

  /**
   * @brief swap child base
   *
   * @todo add swap for Bernstein-Bezier base
   *
   */
  auto set_child_base = [](auto &parent_data, auto &child_data) {
    MoFEMFunctionBeginHot;
    child_data.resetFieldDependentData();
    child_data.getEntDataBitRefLevel() = parent_data.getEntDataBitRefLevel();
    child_data.sPace = parent_data.getSpace();
    child_data.getEntDataBitRefLevel() = parent_data.getEntDataBitRefLevel();

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
          auto child_base = child_data.getNSharedPtr(
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

  auto switch_of_dofs_children = [](auto &parent_ent_data, auto &child_data) {
    MoFEMFunctionBeginHot;
    for (auto i : parent_ent_data.getIndices()) {
      if (i >= 0) {
        for (auto &child_ent_data : child_data) {
          auto it = std::find(child_ent_data.getIndices().begin(),
                              child_ent_data.getIndices().end(), i);
          if (it != child_ent_data.getIndices().end()) {
            const auto dof_idx =
                std::distance(child_ent_data.getIndices().begin(), it);
            child_ent_data.getIndices()[dof_idx] = -1;
            child_ent_data.getLocalIndices()[dof_idx] = -1;
            child_ent_data.getFieldData()[dof_idx] = 0;
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

    BitRefLevel &bit_ent = data.getEntDataBitRefLevel();

    // note all nodes from all added
    if (check(bitParentEnt, bitParentEntMask, bit_ent)) {
      ++count_meshset_sides;

      if (verbosity >= VERBOSE) {
        MOFEM_LOG("SELF", severityLevel)
            << "Add entity data to meshset "
            << "side/type: " << side << "/" << CN::EntityTypeName(type)
            << " op space/base: " << FieldSpaceNames[data.getSpace()] << "/"
            << ApproximationBaseNames[data.getBase()];
      }

      auto &data_on_meshset = entities_field_data.dataOnEntities[MBENTITYSET];
      if (data_on_meshset.size() < count_meshset_sides) {
        if (poolEntsVector.size()) {
          data_on_meshset.transfer(data_on_meshset.end(),
                                   poolEntsVector.begin(), poolEntsVector);
        } else {
          entities_field_data.dataOnEntities[MBENTITYSET].resize(
              count_meshset_sides);
        }
      }

      auto &child_data_meshset =
          entities_field_data
              .dataOnEntities[MBENTITYSET][count_meshset_sides - 1];

      if (opParentType == OPSPACE) {
        CHKERR set_child_base(data, child_data_meshset);
        child_data_meshset.resetFieldDependentData();

      } else {

        auto &field_entities = data.getFieldEntities();

        // note all nodes from all added
        if (field_entities.size() > 1) {
          int dof_idx = 0;
          for (auto dof : data.getFieldDofs()) {
            auto &bit_ent = dof->getBitRefLevel();
            if (!check(bitParentEnt, bitParentEntMask, bit_ent)) {
              data.getIndices()[dof_idx] = -1;
              data.getLocalIndices()[dof_idx] = -1;
              data.getFieldData()[dof_idx] = 0;
            }
            ++dof_idx;
          }
        }

#ifndef NDEBUG
        if (verbosity >= NOISY) {
          for (auto field_ent : field_entities) {
            MOFEM_LOG("SELF", severityLevel)
                << "Parent field entity bit: " << field_ent->getBitRefLevel()
                << " " << *field_ent;
          }
        }
#endif

        if (field_entities.size()) {
          boost::container::static_vector<EntityType, 128> ents_type;
          ents_type.reserve(field_entities.size());
          for (auto fe : field_entities)
            ents_type.push_back(fe->getEntType());
          std::sort(ents_type.begin(), ents_type.end());

          auto end = std::unique(ents_type.begin(), ents_type.end());
          for (auto it_t = ents_type.begin(); it_t != end; ++it_t)
            CHKERR switch_of_dofs_children(
                data, entities_field_data.dataOnEntities[*it_t]);
          if (type == MBENTITYSET)
            CHKERR switch_of_dofs_children(
                data, entities_field_data.dataOnEntities[type]);
        }

        CHKERR set_child_data(data, child_data_meshset);
      }
    }

    MoFEMFunctionReturn(0);
  };

  // iterate parents collect data
  auto &bit_fe = getFEMethod()->numeredEntFiniteElementPtr->getBitRefLevel();
  if (check(bitChild, bitChildMask, bit_fe)) {

    if (verbosity >= VERBOSE) {
      MOFEM_LOG("SELF", severityLevel) << "Child FE bit: " << bit_fe;
    }

    auto loop_parent_fe = [&]() {
      MoFEMFunctionBeginHot;

      parentElePtr->getOpPtrVector().back().doWorkRhsHook = do_work_parent_hook;
      CHKERR loopParent(getFEName(), parentElePtr.get(), verbosity,
                        severityLevel);

#ifndef NDEBUG
      auto &parent_gauss_pts = parentElePtr->gaussPts;
      if (getGaussPts().size1() != parent_gauss_pts.size1()) {
        MOFEM_LOG("SELF", Sev::error) << getGaussPts();
        MOFEM_LOG("SELF", Sev::error) << parent_gauss_pts;
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Wrong number of weights");
      }
      if (getGaussPts().size2() != parent_gauss_pts.size2()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Wrong number of integration points");
      }
#endif

      MoFEMFunctionReturnHot(0);
    };

    CHKERR loop_parent_fe();
  }

  if (verbosity >= VERBOSE) {
    MOFEM_LOG("SELF", severityLevel)
        << "Number of meshset entities "
        << entities_field_data.dataOnEntities[MBENTITYSET].size();
  }

  auto &data_on_meshset = entities_field_data.dataOnEntities[MBENTITYSET];
  auto it = data_on_meshset.begin();

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