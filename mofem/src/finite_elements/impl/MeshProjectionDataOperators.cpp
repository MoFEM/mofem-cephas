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
          NOSPACE, ForcesAndSourcesCore::UserDataOperator::OPLAST),
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
    : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPLAST),
      opParentType(op_parent_type), parentElePtr(parent_ele_ptr),
      bitChild(bit_child), bitChildMask(bit_child_mask),
      bitParentEnt(bit_parent_ent), bitParentEntMask(bit_parent_ent_mask),
      verbosity(verb), severityLevel(sev) {}

MoFEMErrorCode OpAddParentEntData::doWork(int side, EntityType type,
                                          EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  auto check = [](auto &b, auto &m, auto &bit) {
    return ((bit & b).any()) && ((bit & m) == bit);
  };

  auto get_entities_field_data_ptr = [&](auto space) {
    switch (opParentType) {
    case OPROW:
      return getPtrFE()->getDataOnElement()[space];
    case OPCOL:
      return getPtrFE()->getDerivedDataOnElement()[space];
    default:
      return boost::shared_ptr<EntitiesFieldData>();
    }
  };

  auto &bit_fe = getFEMethod()->numeredEntFiniteElementPtr->getBitRefLevel();
  if (check(bitChild, bitChildMask, bit_fe)) {

    auto loop_parent_fe = [&]() {
      MoFEMFunctionBeginHot;

      if (verbosity >= VERBOSE) {
        MOFEM_LOG("SELF", severityLevel)
            << "loop parent element in OpAddParentEntData";
      }

      // note live of op pointer is controlled by ptr_vec in in finite
      // element
      auto field_op = new ForcesAndSourcesCore::UserDataOperator(rowFieldName,
                                                                 opParentType);
      // that forces to run operator at last instance and collect data on
      // entities
      field_op->doWorkRhsHook = [&](DataOperator *op_ptr, int side,
                                    EntityType type,
                                    EntitiesFieldData::EntData &data) {
        MoFEMFunctionBegin;


        auto field_entities = data.getFieldEntities();
        if (field_entities.size() == 1) {
          auto &bit_ent = field_entities[0]->getBitRefLevel();
          if (!check(bitParentEnt, bitParentEntMask, bit_ent))
            MoFEMFunctionReturnHot(0);
        } 

        auto space = data.getSpace();
        auto base = data.getBase();

        auto entities_field_data_ptr = get_entities_field_data_ptr(space);
        if (entities_field_data_ptr)
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

        MoFEMFunctionReturn(0);
      };

      parentElePtr->getOpPtrVector().push_back(field_op);
      CHKERR loopParent(getFEName(), parentElePtr.get(), verbosity,
                        severityLevel);
      // clean element from obsolete data operator
      parentElePtr->getOpPtrVector().pop_back();

      MoFEMFunctionReturnHot(0);
    };

    CHKERR loop_parent_fe();
  }

  MoFEMFunctionReturn(0);
}

OpRestoreEntData::OpRestoreEntData(FieldSpace space, OpType op_type)
    : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPLAST), sPace(space),
      opType(op_type) {}

MoFEMErrorCode OpRestoreEntData::doWork(int side, EntityType type,
                                        EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  switch (opType) {
  case OPROW:
    getPtrFE()->getDataOnElement()[sPace]->dataOnEntities[MBENTITYSET].clear();
    break;
  case OPCOL:
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

} // namespace MoFEM