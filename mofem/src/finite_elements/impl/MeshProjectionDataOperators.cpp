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

// OpAddParentEntData::OpAddParentEntData(
//     std::string field_name, OpType op_parent_type,
//     boost::shared_ptr<ForcesAndSourcesCore> parent_ele_ptr,
//     BitRefLevel bit_child, BitRefLevel bit_child_mask,
//     BitRefLevel bit_parent_ent, BitRefLevel bit_parent_ent_mask, int verb,
//     Sev sev)
//     : ForcesAndSourcesCore::UserDataOperator(field_name, op_parent_type),
//       opParentType(op_parent_type), parentElePtr(parent_ele_ptr),
//       bitChild(bit_child), bitChildMask(bit_child_mask),
//       bitParentEnt(bit_parent_ent), bitParentEntMask(bit_parent_ent_mask),
//       verbosity(verb), severityLevel(sev) {
//   std::fill(doEntities.begin(), doEntities.end(), false);
//   doEntities[MBENTITYSET] = true;
// }

// MoFEMErrorCode OpAddParentEntData::doWork(int side, EntityType type,
//                                           EntitiesFieldData::EntData &data) {
//   MoFEMFunctionBegin;

//   // auto check = [](auto &b, auto &m, auto &bit) {
//   //   return ((bit & b).any()) && ((bit & m) == bit);
//   // };

//   // auto &bit_fe = getFEMethod()->numeredEntFiniteElementPtr->getBitRefLevel();
//   // if (check(bitChild, bitChildMask, bit_fe)) {

//   //   auto loop_parent_fe = [&]() {
//   //     MoFEMFunctionBeginHot;
//   //     if (verbosity >= VERBOSE) {
//   //       MOFEM_LOG("SELF", severityLevel)
//   //           << "loop parent element in OpAddParentEntData";
//   //     }
//   //     // note live of op pointer is controlled by ptr_vec in in finite
//   //     // element
//   //     auto field_op = new ForcesAndSourcesCore::UserDataOperator(rowFieldName,
//   //                                                                opParentType);
//   //     // that forces to run operator at last instance and collect data on
//   //     // entities
//   //     field_op->doWorkRhsHook =
//   //         [&](DataOperator *op_ptr, int side, EntityType type,
//   //             EntitiesFieldData::EntData &data) { return 0; };
//   //     parentElePtr->getOpPtrVector().push_back(field_op);
//   //     CHKERR loopParent(getFEName(), parentElePtr.get(), verbosity,
//   //                       severityLevel);
//   //     // clean element from obsolete data operator
//   //     parentElePtr->getOpPtrVector().pop_back();
//   //     MoFEMFunctionReturnHot(0);
//   //   };

//   //   auto

//   //   auto base = data.getBase();
//   //   auto space = data.getSpace();
//   //   auto entities_field_data = getPtrFE()->getDataOnElement()[space];

//   //   auto ent_finite_element = parentElePtr->getFEEntityHandle();
//   //   auto parent_entities_field_data = parentElePtr->getDataOnElement()[space];
//   //   auto &side_parent_table =
//   //       parentElePtr->numeredEntFiniteElementPtr->getSideNumberTable();

//   //   if (side == 0 && type == MBVERTEX) {

//   //     auto ent_data = getPtrFE()->getDataOnElement();
//   //     auto derived_data_on_element = getDerivedDataOnElement();

//   //     if (op.getOpType() & UserDataOperator::OPROW ||
//   //         op.getOpType() & UserDataOperator::OPROWCOL) {

//   //       for (auto &entities_field_data_ptr : ent_data) {
//   //         if(entities_field_data_ptr) {
//   //           for(auto &ent_data : *entities_field_data_ptr) {
//   //             if(ent_data->getFieldEntities().size() == 1) {
//   //               auto ent_data_bit = ent_data->getBitRefLevel();
//   //               if(check(bitParentEnt,bitParentEntMask, ent_data_bit)) {
//   //                 ent_data[space][MBENTITYSET].push_back(ent_data);
//   //               }
//   //             } else if(ent_data->getFieldEntities().size() > 1) {
//   //               ent_data[space][MBENTITYSET].push_back(ent_data);
//   //               auto &back_end_data = ent_data[space][MBENTITYSET].back();
//   //               int dd = 0;
//   //               for(auto dof : ent_data->getFieldDofs()) {
//   //                 auto &dof_bit = dof->getBitRefLevel();
//   //                 if (!check(bitParentEnt, bitParentEntMask, dof_bit)) {
//   //                   back_ent_data->getFieldData()[dd] = 0;
//   //                   back_ent_data->getIndices()[dd] = 1;
//   //                 }
//   //                 ++dd;
//   //               }
//   //             } 
//   //           }
//   //         }
//   //       }
//   //     }

//   //     if (op.getOpType() & UserDataOperator::OPCOL ||
//   //         op.getOpType() & UserDataOperator::OPROWCOL) {
//   //     }

//   //     if (parent_entities_field_data) {
//   //       for (auto &type_data : *parent_entities_field_data) {
//   //         for (auto &side_data : type_data) {
//   //           for (auto field_ent : side_data->getFieldEntities()) {
//   //             auto &bit_parent_ent = field_ent->getBitRefLevel();
//   //             auto parent_type = field_ent->getEntType();
//   //             if (check(bitParentEnt, bitParentEntMask, bit_parent_ent)) {


//   //             }
//   //           }
//   //         }
//   //       }
//   //     }
//   //   }

//   //   // auto base_swap = [&]() {
//   //   //   MoFEMFunctionBeginHot;
//   //   //   if (base == AINSWORTH_BERNSTEIN_BEZIER_BASE) {
//   //   //     CHKERR entities_field_data->baseSwap(rowFieldName,
//   //   //                                          AINSWORTH_BERNSTEIN_BEZIER_BASE);
//   //   //     CHKERR parent_entities_field_data->baseSwap(
//   //   //         rowFieldName, AINSWORTH_BERNSTEIN_BEZIER_BASE);
//   //   //   }
//   //   //   MoFEMFunctionReturnHot(0);
//   //   // };
//   // }

//   MoFEMFunctionReturn(0);
// }

// OpRetoreEntData::OpRetoreEntData()
//     : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPLAST) {}

// MoFEMErrorCode OpRetoreEntData::doWork(int side, EntityType type,
//                                        EntitiesFieldData::EntData &data) {
//   MoFEMFunctionBegin;
//   MoFEMFunctionReturn(0);
// }

} // namespace MoFEM