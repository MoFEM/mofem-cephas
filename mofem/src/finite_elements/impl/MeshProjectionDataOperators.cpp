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
    : ForcesAndSourcesCore::UserDataOperator(field_name, op_parent_type),
      fieldName(field_name), opParentType(op_parent_type),
      
      
      parentElePtr(parent_ele_ptr), bitChild(bit_child),
      bitChildMask(bit_child_mask), bitParentEnt(bit_parent_ent),
      bitParentEntMask(bit_parent_ent_mask), verbosity(verb),
      severityLevel(sev) {}

MoFEMErrorCode OpAddParentEntData::opRhs(EntitiesFieldData &entities_field_data,
                                         const bool error_if_no_base) {
  int count_meshset_sides = 0;
  MoFEMFunctionBegin;

  auto check = [](auto &b, auto &m, auto &bit) {
    return ((bit & b).any()) && ((bit & m) == bit);
  };

  auto set_child_data = [](auto &parent_data, auto &child_data) {
    MoFEMFunctionBeginHot;
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
    child_data.sPace = parent_data.getSpace();
    child_data.getEntDataBitRefLevel() = parent_data.getEntDataBitRefLevel();

#ifndef NDEBUG
    if (child_data.bAse == AINSWORTH_BERNSTEIN_BEZIER_BASE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Base not implemented");
#endif

    using BaseDerivatives = EntitiesFieldData::EntData::BaseDerivatives;
    for (auto direvative = 0; direvative != BaseDerivatives::LastDerivative;
         ++direvative) {
      auto &parent_base = parent_data.getNSharedPtr(
          parent_data.getBase(), static_cast<BaseDerivatives>(direvative));
      if (parent_base) {
        auto &child_base = child_data.getNSharedPtr(
            parent_data.getBase(), static_cast<BaseDerivatives>(direvative));
        if (!child_base)
          child_base = boost::make_shared<MatrixDouble>();
        child_base->swap(*parent_base);
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

    auto field_entities = data.getFieldEntities();
    if (field_entities.size()) {

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

        if (data.getFieldEntities().size()) {
          // note: assumes that all entities on data are the same type. that
          // is good assumption, since we group entities on side by type.
          const auto ent_type = data.getFieldEntities()[0]->getEntType();
          CHKERR switch_of_dofs_children(
              data, entities_field_data.dataOnEntities[ent_type]);
          if (ent_type != type)
            CHKERR switch_of_dofs_children(
                data, entities_field_data.dataOnEntities[type]);
        }

        if (entities_field_data.dataOnEntities[MBENTITYSET].size() <
            count_meshset_sides)
          entities_field_data.dataOnEntities[MBENTITYSET].resize(
              count_meshset_sides);

        auto &child_data_meshset =
            entities_field_data
                .dataOnEntities[MBENTITYSET][count_meshset_sides - 1];

        CHKERR set_child_data(data, child_data_meshset);
        CHKERR set_child_base(data, child_data_meshset);
      }
    }

    MoFEMFunctionReturn(0);
  };

  auto &bit_fe = getFEMethod()->numeredEntFiniteElementPtr->getBitRefLevel();
  if (check(bitChild, bitChildMask, bit_fe)) {

    if (verbosity >= VERBOSE) {
      MOFEM_LOG("SELF", severityLevel) << "Child FE bit: " << bit_fe;
    }

    auto loop_parent_fe = [&]() {
      MoFEMFunctionBeginHot;

      // note live of op pointer is controlled by ptr_vec in in finite
      // element
      auto field_op =
          new ForcesAndSourcesCore::UserDataOperator(fieldName, opParentType);
      field_op->doWorkRhsHook = do_work_parent_hook;

      parentElePtr->getOpPtrVector().push_back(field_op);
      CHKERR loopParent(getFEName(), parentElePtr.get(), verbosity,
                        severityLevel);
      // clean element from obsolete data operator
      parentElePtr->getOpPtrVector().pop_back();

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

  entities_field_data.dataOnEntities[MBENTITYSET].resize(count_meshset_sides);

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM