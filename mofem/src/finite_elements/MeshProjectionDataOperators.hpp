/** \file MeshProjectionDataOperators.hpp
  * \brief Mesh projection operators

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

#ifndef __MESH_PROJECTION_DATA_OPERATORS_HPP__
#define __MESH_PROJECTION_DATA_OPERATORS_HPP__

namespace MoFEM {

/**
 * @brief Operator to execute finite element instance on parent element.
 * This operator is typically used to project field from parent to child, or
 * vice versa. It enables to evaluate filed data of parent entity on chile
 * entity integration points.
 */
struct OpRunParent : public ForcesAndSourcesCore::UserDataOperator {

  /**
   * @brief Construct a new Op Run Parent object
   *
   * @note Finite element instance usually has to be class which has overloaded
   * method from projevting integration points from child tp parent.
   *
   * @note Typically parent_ele_ptr and bit_this_mask is the same instance
   *
   * @param parent_ele_ptr finite element instance executed on parent entity
   * @param bit_parent bit of parent entity
   * @param bit_parent_mask mask of parent entity
   * @param this_ele_ptr "this" element instance
   * @param bit_this bit of entity on which "this" finite element is executed
   * @param bit_this_mask mask of entity on which "this" finite element instance
   * is executed
   * @param verb verbosity level
   * @param sev logging severity level
   */
  OpRunParent(boost::shared_ptr<ForcesAndSourcesCore> parent_ele_ptr,
              BitRefLevel bit_parent, BitRefLevel bit_parent_mask,
              boost::shared_ptr<ForcesAndSourcesCore> this_ele_ptr,
              BitRefLevel bit_this, BitRefLevel bit_this_mask, int verb = QUIET,
              Sev sev = Sev::noisy);

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  boost::shared_ptr<ForcesAndSourcesCore> parentElePtr;
  boost::shared_ptr<ForcesAndSourcesCore> thisElePtr;
  const BitRefLevel bitParent;
  const BitRefLevel bitParentMask;
  const BitRefLevel bitThis;
  const BitRefLevel bitThisMask;
  int verbosity;
  Sev severityLevel;
};

/**
 * @brief Operator to project base functions from parent entity
 *
 * This operator project base functions, field data (i.e. indices, field values
 * of dofs, etc.), into parent element. Operator can be called as a hierarchy to
 * get access to information on lower refinement levels.
 *
 */
struct OpAddParentEntData : public ForcesAndSourcesCore::UserDataOperator {

  /**
   * @brief Construct a new Op Add Parent Ent Data object
   * 
   * @param field_name field name DOFs projected from parent
   * @param op_parent_type  type of user data operator
   * @param parent_ele_ptr parent finite element instance
   * @param bit_child bit of child entity
   * @param bit_child_mask bit mask of child
   * @param bit_parent_ent bit of parent entity
   * @param bit_parent_ent_mask bit mask of parent
   * @param verb verbosity level
   * @param sev severity level for logging
   */
  OpAddParentEntData(std::string field_name, OpType op_parent_type,
                     boost::shared_ptr<ForcesAndSourcesCore> parent_ele_ptr,
                     BitRefLevel bit_child, BitRefLevel bit_child_mask,
                     BitRefLevel bit_parent_ent,
                     BitRefLevel bit_parent_ent_mask, int verb = QUIET,
                     Sev sev = Sev::noisy);

  MoFEMErrorCode opRhs(EntitiesFieldData &data,
                       const bool error_if_no_base = false);

private:
  std::string fieldName;
  OpType opParentType;
  boost::shared_ptr<ForcesAndSourcesCore> parentElePtr;
  const BitRefLevel bitChild;
  const BitRefLevel bitChildMask;
  const BitRefLevel bitParentEnt;
  const BitRefLevel bitParentEntMask;
  int verbosity;
  Sev severityLevel;
};

} // namespace MoFEM

#endif //__MESH_PROJECTION_DATA_OPERATORS_HPP__