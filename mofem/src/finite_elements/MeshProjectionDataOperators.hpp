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

struct OpRunParent : public ForcesAndSourcesCore::UserDataOperator {
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

struct OpAddParentEntData : public ForcesAndSourcesCore::UserDataOperator {

  OpAddParentEntData(std::string field_name, OpType op_parent_type,
                     boost::shared_ptr<ForcesAndSourcesCore> parent_ele_ptr,
                     BitRefLevel bit_child, BitRefLevel bit_child_mask,
                     BitRefLevel bit_parent_ent,
                     BitRefLevel bit_parent_ent_mask, int verb = QUIET,
                     Sev sev = Sev::noisy);

  MoFEMErrorCode opRhs(EntitiesFieldData &data,
                       const bool error_if_no_base = false);

  static MoFEMErrorCode
  markHangingSkinParents(MoFEM::Interface &m_field, const int dim,
                  const BitRefLevel parent_bit, const BitRefLevel parent_mask,
                  const BitRefLevel child_bit, const BitRefLevel child_mask,
                  const BitRefLevel mark_bit, const bool resolve_shared,
                  const std::string debug_file_name = "");

  static MoFEMErrorCode markHangingSkinChildren(
      MoFEM::Interface &m_field, const BitRefLevel child_bit,
      const BitRefLevel child_mask, const BitRefLevel mark_bit,
      const BitRefLevel mark_mask, const std::string debug_file_name = "");
      

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