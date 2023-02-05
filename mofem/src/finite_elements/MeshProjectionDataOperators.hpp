/** \file MeshProjectionDataOperators.hpp
  * \brief Mesh projection operators

*/

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
   * method from projecting integration points from child tp parent.
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
 * @brief Operator to project base functions from parent entity to child
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

  boost::ptr_deque<EntitiesFieldData::EntData> poolEntsVector;
};

/**
 * @brief Create adjacency to parent elements.
 *
 * That class is used during entity finite element construction.
 *
 * @tparam DIM dimension of parent element
 */
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
   * @note elements form child, see dofs from parent, so DOFs located on
   * adjacencies of parent entity has adjacent to dofs of child.
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

    adjacency.clear();

    if (field.getSpace() != NOFIELD) {

      auto basic_entity_data_ptr = fe.getBasicDataPtr();
      auto th_parent_handle = basic_entity_data_ptr->th_RefParentHandle;
      auto th_bit_level = basic_entity_data_ptr->th_RefBitLevel;

      std::vector<EntityHandle> parents;
      parents.reserve(BITREFLEVEL_SIZE);

      CHKERR getParent(fe.getEnt(), parents, moab, th_parent_handle,
                       th_bit_level);

      CHKERR getParentsAdjacencies(field, moab, parents, adjacency);
    }

    adjTmp.clear();
    CHKERR getDefaultAdjacencies(moab, field, fe, adjTmp);
    adjacency.insert(adjacency.end(), adjTmp.begin(), adjTmp.end());

    std::sort(adjacency.begin(), adjacency.end());
    auto it = std::unique(adjacency.begin(), adjacency.end());
    adjacency.resize(std::distance(adjacency.begin(), it));

    for (auto e : adjacency) {
      auto &side_table = fe.getSideNumberTable();
      if (side_table.find(e) == side_table.end())
        const_cast<SideNumber_multiIndex &>(side_table)
            .insert(boost::shared_ptr<SideNumber>(new SideNumber(e, -1, 0, 0)));
    }

    MoFEMFunctionReturn(0);
  }

protected:
  MoFEMErrorCode getParent(EntityHandle fe, std::vector<EntityHandle> &parents,
                           moab::Interface &moab, Tag th_parent_handle,
                           Tag th_bit_level) {
    MoFEMFunctionBegin;

    auto check = [](auto &b, auto &m, auto &bit) {
      return ((bit & b).any()) && ((bit & m) == bit);
    };

    BitRefLevel bit_fe;
    CHKERR moab.tag_get_data(th_bit_level, &fe, 1, &bit_fe);
    if (check(bitEnt, bitEntMask, bit_fe)) {

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
          if ((fe_parent != 0) && (fe_parent != fe) &&
              (parent_type == back_type)) {
            parents.push_back(fe_parent);
            CHKERR get_parent(parents.back(), parents);
          }
        }
        MoFEMFunctionReturn(0);
      };

      CHKERR get_parent(fe, parents);
    }
    MoFEMFunctionReturn(0);
  }

  MoFEMErrorCode getParentsAdjacencies(const Field &field,
                                       moab::Interface &moab,
                                       std::vector<EntityHandle> &parents,
                                       std::vector<EntityHandle> &adjacency) {
    MoFEMFunctionBegin;

    switch (field.getSpace()) {
    case H1:
      for (auto fe_ent : parents)
        CHKERR moab.get_connectivity(&*parents.begin(), parents.size(),
                                     adjacency, true);
    case HCURL:
      if constexpr (DIM >= 1)
        for (auto fe_ent : parents)
          CHKERR moab.get_adjacencies(&fe_ent, 1, 1, false, adjacency,
                                      moab::Interface::UNION);
    case HDIV:
      if constexpr (DIM == 2)
        for (auto fe_ent : parents)
          CHKERR moab.get_adjacencies(&fe_ent, 1, 2, false, adjacency,
                                      moab::Interface::UNION);
    case L2:
      for (auto fe_ent : parents)
        adjacency.push_back(fe_ent);
      break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
              "this field is not implemented for face finite element");
    }

    if (adjacency.size()) {
      std::sort(adjacency.begin(), adjacency.end());
      auto it = std::unique(adjacency.begin(), adjacency.end());
      adjacency.resize(std::distance(adjacency.begin(), it));
    }

    MoFEMFunctionReturn(0);
  }

  MoFEMErrorCode getDefaultAdjacencies(moab::Interface &moab,
                                       const Field &field,
                                       const EntFiniteElement &fe,
                                       std::vector<EntityHandle> &adjacency) {
    MoFEMFunctionBegin;
    if constexpr (DIM == 3)
      CHKERR DefaultElementAdjacency::defaultVolume(moab, field, fe, adjacency);
    if constexpr (DIM == 2)
      CHKERR DefaultElementAdjacency::defaultFace(moab, field, fe, adjacency);
    else if constexpr (DIM == 1)
      CHKERR DefaultElementAdjacency::defaultEdge(moab, field, fe, adjacency);
    else if constexpr (DIM == 0)
      CHKERR DefaultElementAdjacency::defaultVertex(moab, field, fe, adjacency);
    MoFEMFunctionReturn(0);
  };

  BitRefLevel bitParent;
  BitRefLevel bitParentMask;
  BitRefLevel bitEnt;
  BitRefLevel bitEntMask;
  std::vector<EntityHandle> adjTmp;
};

/**
 * @brief Create adjacency to parent skeleton elements.
 *
 * That class is used during entity finite element construction.
 *
 * @tparam DIM dimension of parent element
 */
template <int DIM>
struct ParentFiniteElementAdjacencyFunctionSkeleton
    : public ParentFiniteElementAdjacencyFunction<DIM> {

  ParentFiniteElementAdjacencyFunctionSkeleton(BitRefLevel bit_parent,
                                               BitRefLevel bit_parent_mask,
                                               BitRefLevel bit_ent,
                                               BitRefLevel bit_ent_mask)
      : ParentFiniteElementAdjacencyFunction<DIM>(bit_parent, bit_parent_mask,
                                                  bit_ent, bit_ent_mask) {}

  /**
   * @brief Function setting adjacencies to DOFs of parent element
   *
   * @note elements form child, see dofs from parent, so DOFs located on
   * adjacencies of parent entity has adjacent to dofs of child.
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

    adjacency.clear();
    CHKERR this->getDefaultAdjacencies(moab, field, fe, adjacency);

    if (field.getSpace() != NOFIELD) {

      const auto fe_ent = fe.getEnt();
      brideAdjacencyEdge.clear();
      CHKERR moab.get_adjacencies(&fe_ent, 1, DIM + 1, false,
                                  brideAdjacencyEdge);

      std::vector<EntityHandle> parents;

      if (this->bitParent.any()) {
        auto basic_entity_data_ptr = fe.getBasicDataPtr();
        auto th_parent_handle = basic_entity_data_ptr->th_RefParentHandle;
        auto th_bit_level = basic_entity_data_ptr->th_RefBitLevel;
        parents.reserve(BITREFLEVEL_SIZE);
        for (auto bridge_fe : brideAdjacencyEdge) {
          CHKERR this->getParent(bridge_fe, parents, moab, th_parent_handle,
                                 th_bit_level);
        };
        parents.insert(parents.end(), brideAdjacencyEdge.begin(),
                       brideAdjacencyEdge.end());
      } else {
        parents.swap(brideAdjacencyEdge);
      }

      CHKERR this->getParentsAdjacencies(field, moab, parents, adjacency);

      std::sort(adjacency.begin(), adjacency.end());
      auto it = std::unique(adjacency.begin(), adjacency.end());
      adjacency.resize(std::distance(adjacency.begin(), it));

      for (auto e : adjacency) {
        auto &side_table = fe.getSideNumberTable();
        if (side_table.find(e) == side_table.end())
          const_cast<SideNumber_multiIndex &>(side_table)
              .insert(
                  boost::shared_ptr<SideNumber>(new SideNumber(e, -1, 0, 0)));
      }
    }

    MoFEMFunctionReturn(0);
  }

protected:
  std::vector<EntityHandle> brideAdjacencyEdge;
};

} // namespace MoFEM

#endif //__MESH_PROJECTION_DATA_OPERATORS_HPP__