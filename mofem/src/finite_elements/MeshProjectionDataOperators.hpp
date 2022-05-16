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

  /**
   * @brief Construct a new Op Add Parent Ent Data object
   *
   * @param space field space
   * @param op_parent_type  type of user data operator
   * @param parent_ele_ptr parent finite element instance
   * @param bit_child bit of child entity
   * @param bit_child_mask bit mask of child
   * @param bit_parent_ent bit of parent entity
   * @param bit_parent_ent_mask bit mask of parent
   * @param verb verbosity level
   * @param sev severity level for logging
   */
  OpAddParentEntData(FieldSpace space, OpType op_parent_type,
                     boost::shared_ptr<ForcesAndSourcesCore> parent_ele_ptr,
                     BitRefLevel bit_child, BitRefLevel bit_child_mask,
                     BitRefLevel bit_parent_ent,
                     BitRefLevel bit_parent_ent_mask, int verb = QUIET,
                     Sev sev = Sev::noisy);

  MoFEMErrorCode opRhs(EntitiesFieldData &data,
                       const bool error_if_no_base = false);

private:
  std::string fieldName;
  FieldSpace approxSpace;
  OpType opParentType;
  boost::shared_ptr<ForcesAndSourcesCore> parentElePtr;
  const BitRefLevel bitChild;
  const BitRefLevel bitChildMask;
  const BitRefLevel bitParentEnt;
  const BitRefLevel bitParentEntMask;
  int verbosity;
  Sev severityLevel;
};

template <int DIM> struct ParentFiniteElementAdjacencyFunction {

  ParentFiniteElementAdjacencyFunction(BitRefLevel bit_parent,
                                       BitRefLevel bit_parent_mask,
                                       BitRefLevel bit_ent,
                                       BitRefLevel bit_ent_mask)
      : bitParent(bit_parent), bitParentMask(bit_parent_mask), bitEnt(bit_ent),
        bitEntMask(bit_ent_mask) {}

  /**
   * @brief Function setting adjacencies to DOFs of parent element
   *
   * @note elements form child, see dofs from parent, so DOFs localted on
   * adjacencies of parent entity has adjacent to dofs of chiold.
   *
   * @tparam DIM dimension of the element entity
   * @param moab
   * @param field
   * @param fe
   * @param adjacency
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode operator()(moab::Interface &moab, const Field &field,
                            const EntFiniteElement &fe,
                            std::vector<EntityHandle> &adjacency) {
    MoFEMFunctionBegin;

    static_assert(DIM >= 0 && DIM <= 3, "DIM is out of scope");

    auto check = [](auto &b, auto &m, auto &bit) {
      return ((bit & b).any()) && ((bit & m) == bit);
    };

    adjTmp.clear();

    if (field.getSpace() != NOFIELD) {

      auto basic_entity_data_ptr = fe.getBasicDataPtr();
      auto th_parent_handle = basic_entity_data_ptr->th_RefParentHandle;
      auto th_bit_level = basic_entity_data_ptr->th_RefBitLevel;

      using GetParent = boost::function<MoFEMErrorCode(
          EntityHandle fe, std::vector<EntityHandle> & parents)>;

      /**
       * @brief this function os called recursively, until all stack of parents
       * is found.
       *
       */
      GetParent get_parent = [&](EntityHandle fe,
                                 std::vector<EntityHandle> &parents) {
        MoFEMFunctionBegin;
        EntityHandle fe_parent;

        CHKERR moab.tag_get_data(th_parent_handle, &fe, 1, &fe_parent);
        auto parent_type = type_from_handle(fe_parent);
        auto back_type = type_from_handle(fe);
        BitRefLevel bit_parent;
        CHKERR moab.tag_get_data(th_bit_level, &fe_parent, 1, &bit_parent);
        if (check(bitParent, bitParentMask, bit_parent)) {
          if (fe_parent != 0 && fe_parent != fe && parent_type == back_type) {
            parents.push_back(fe_parent);
            CHKERR get_parent(parents.back(), parents);
          }
        }
        MoFEMFunctionReturn(0);
      };

      std::vector<EntityHandle> parents;
      parents.reserve(BITREFLEVEL_SIZE);

      CHKERR get_parent(fe.getEnt(), parents);

      for (auto fe_ent : parents) {
        switch (field.getSpace()) {
        case H1:
          CHKERR moab.get_adjacencies(&fe_ent, 1, 0, false, adjacency,
                                      moab::Interface::UNION);
        case HCURL:
          if constexpr (DIM >= 2)
            CHKERR moab.get_adjacencies(&fe_ent, 1, 1, false, adjacency,
                                        moab::Interface::UNION);
        case HDIV:
          if constexpr (DIM == 3)
            CHKERR moab.get_adjacencies(&fe_ent, 1, 2, false, adjacency,
                                        moab::Interface::UNION);
        case L2:
          adjacency.push_back(fe_ent);
          break;
        default:
          SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                  "this field is not implemented for face finite element");
        }
      }

      if (adjacency.size()) {

        std::sort(adjacency.begin(), adjacency.end());
        auto it = std::unique(adjacency.begin(), adjacency.end());
        adjacency.resize(std::distance(adjacency.begin(), it));
        bitLevels.resize(adjacency.size());
        CHKERR moab.tag_get_data(th_bit_level, &*adjacency.begin(),
                                 adjacency.size(), &*bitLevels.begin());


        adjTmp.reserve(adjacency.size());
        for (int i = 0; i != adjacency.size(); ++i) {
          const auto &bit = bitLevels[i];
          if (check(bitEnt, bitEntMask, bit)) {
            adjTmp.push_back(adjacency[i]);
          }
        }

      }
    }

    adjacency.clear();

    if constexpr (DIM == 3)
      CHKERR DefaultElementAdjacency::defaultVolume(moab, field, fe, adjacency);
    if constexpr (DIM == 2)
      CHKERR DefaultElementAdjacency::defaultFace(moab, field, fe, adjacency);
    else if constexpr (DIM == 1)
      CHKERR DefaultElementAdjacency::defaultEdge(moab, field, fe, adjacency);
    else if constexpr (DIM == 0)
      CHKERR DefaultElementAdjacency::defaultVertex(moab, field, fe, adjacency);

    adjacency.insert(adjacency.end(), adjTmp.begin(), adjTmp.end());

    std::sort(adjacency.begin(), adjacency.end());
    auto it = std::unique(adjacency.begin(), adjacency.end());
    adjacency.resize(std::distance(adjacency.begin(), it));

    for (auto e : adjacency) {
      auto side_table = fe.getSideNumberTable();
      if (side_table.find(e) == side_table.end())
        const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
            .insert(boost::shared_ptr<SideNumber>(new SideNumber(e, -1, 0, 0)));
    }

    MoFEMFunctionReturn(0);
  }

private:
  BitRefLevel bitParent;
  BitRefLevel bitParentMask;
  BitRefLevel bitEnt;
  BitRefLevel bitEntMask;
  std::vector<EntityHandle> adjTmp;
  std::vector<BitRefLevel> bitLevels;
};

} // namespace MoFEM

#endif //__MESH_PROJECTION_DATA_OPERATORS_HPP__