/** \file ForcesAndSourcesCore.cpp

\brief Implementation of Elements on Entities for Forces and Sources
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

#ifdef __cplusplus
extern "C" {
#endif
#include <cblas.h>
#include <lapack_wrap.h>
// #include <gm_rule.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

ForcesAndSourcesCore::ForcesAndSourcesCore(Interface &m_field)
    :

      mField(m_field), getRuleHook(0), setRuleHook(0),
      dataOnElement{

          nullptr,
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET), // NOFIELD
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET), // H1
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET), // HCURL
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET), // HDIV
          boost::make_shared<DataForcesAndSourcesCore>(MBENTITYSET)  // L2

      },
      derivedDataOnElement{

          nullptr,
          boost::make_shared<DerivedDataForcesAndSourcesCore>(
              dataOnElement[NOFIELD]), // NOFIELD
          boost::make_shared<DerivedDataForcesAndSourcesCore>(
              dataOnElement[H1]), // H1
          boost::make_shared<DerivedDataForcesAndSourcesCore>(
              dataOnElement[HCURL]), // HCURL
          boost::make_shared<DerivedDataForcesAndSourcesCore>(
              dataOnElement[HDIV]), // HDIV
          boost::make_shared<DerivedDataForcesAndSourcesCore>(
              dataOnElement[L2]) // L2

      },
      dataNoField(*dataOnElement[NOFIELD].get()),
      dataH1(*dataOnElement[H1].get()), dataHcurl(*dataOnElement[HCURL].get()),
      dataHdiv(*dataOnElement[HDIV].get()), dataL2(*dataOnElement[L2].get()),
      lastEvaluatedElementEntityType(MBMAXTYPE), sidePtrFE(nullptr) {}

// ** Sense **

MoFEMErrorCode ForcesAndSourcesCore::getEntitySense(
    const EntityType type,
    boost::ptr_vector<DataForcesAndSourcesCore::EntData> &data) const {
  MoFEMFunctionBegin;

  auto &side_table = numeredEntFiniteElementPtr->getSideNumberTable().get<0>();
  auto sit = side_table.lower_bound(get_id_for_min_type(type));
  if (sit != side_table.end()) {
    auto hi_sit = side_table.upper_bound(get_id_for_max_type(type));
    for (; sit != hi_sit; ++sit) {
      const int side_number = (*sit)->side_number;
      if (side_number >= 0) {
        const int brother_side_number = (*sit)->brother_side_number;
        const int sense = (*sit)->sense;

        data[side_number].getSense() = sense;
        if (brother_side_number != -1)
          data[brother_side_number].getSense() = sense;
      }
    }
  }
  MoFEMFunctionReturn(0);
}

// ** Order **

template <typename ENTMULTIINDEX>
static inline int getMaxOrder(const ENTMULTIINDEX &multi_index) {
  int max_order = 0;
  for (auto ent_field_weak_ptr : multi_index)
    if (auto e = ent_field_weak_ptr.lock()) {
      const int order = e->getMaxOrder();
      max_order = (max_order < order) ? order : max_order;
    }
  return max_order;
}

int ForcesAndSourcesCore::getMaxDataOrder() const {
  int max_order = 0;
  for (auto e : getDataFieldEnts()) {
    const int order = e->getMaxOrder();
    max_order = (max_order < order) ? order : max_order;
  }
  return max_order;
}

int ForcesAndSourcesCore::getMaxRowOrder() const {
  return getMaxOrder(getRowFieldEnts());
}

int ForcesAndSourcesCore::getMaxColOrder() const {
  return getMaxOrder(getColFieldEnts());
}

MoFEMErrorCode ForcesAndSourcesCore::getEntityDataOrder(
    const EntityType type, const FieldSpace space,
    boost::ptr_vector<DataForcesAndSourcesCore::EntData> &data) const {
  MoFEMFunctionBegin;

  auto set_order = [&]() {
    MoFEMFunctionBeginHot;
    auto &side_table = numeredEntFiniteElementPtr->getSideNumberTable();

    for (unsigned int s = 0; s != data.size(); ++s)
      data[s].getDataOrder() = 0;

    auto &fields_ents =
        getDataFieldEnts().get<Composite_EntType_and_Space_mi_tag>();

    for (auto r = fields_ents.equal_range(boost::make_tuple(type, space));
         r.first != r.second; ++r.first) {

      auto &e = **r.first;

      auto sit = side_table.find(e.getEnt());
      if (sit != side_table.end()) {

        auto &side = *sit;
        const int side_number = side->side_number;
        if (side_number >= 0) {
          ApproximationOrder ent_order = e.getMaxOrder();
          auto &dat = data[side_number];
          dat.getDataOrder() =
              dat.getDataOrder() > ent_order ? dat.getDataOrder() : ent_order;
        }
      } else
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Entity on side of the element not found");
    }
    MoFEMFunctionReturnHot(0);
  };

  auto set_order_on_brother = [&]() {
    MoFEMFunctionBeginHot;
    auto &side_table =
        numeredEntFiniteElementPtr->getSideNumberTable().get<0>();
    auto sit = side_table.lower_bound(get_id_for_min_type(type));
    if (sit != side_table.end()) {
      auto hi_sit = side_table.upper_bound(get_id_for_max_type(type));
      for (; sit != hi_sit; ++sit) {
        const int brother_side_number = (*sit)->brother_side_number;
        if (brother_side_number != -1) {
          const int side_number = (*sit)->side_number;
          data[brother_side_number].getDataOrder() =
              data[side_number].getDataOrder();
        }
      }
    }
    MoFEMFunctionReturnHot(0);
  };

  CHKERR set_order();
  CHKERR set_order_on_brother();

  MoFEMFunctionReturn(0);
}

// ** Indices **

MoFEMErrorCode ForcesAndSourcesCore::getNodesIndices(
    const boost::string_ref field_name, FENumeredDofEntity_multiIndex &dofs,
    VectorInt &nodes_indices, VectorInt &local_nodes_indices) const {
  MoFEMFunctionBegin;

  auto &dofs_by_uid = dofs.get<Unique_mi_tag>();
  auto bit_number = mField.get_field_bit_number(&field_name[0]);
  auto dit = dofs_by_uid.lower_bound(DofEntity::getLoFieldEntityUId(
      bit_number, get_id_for_min_type<MBVERTEX>()));
  auto hi_dit = dofs_by_uid.upper_bound(DofEntity::getHiFieldEntityUId(
      bit_number, get_id_for_max_type<MBVERTEX>()));

  if (dit != hi_dit) {

    int num_nodes;
    CHKERR getNumberOfNodes(num_nodes);
    int max_nb_dofs = 0;
    const int nb_dofs_on_vert = (*dit)->getNbOfCoeffs();
    max_nb_dofs = nb_dofs_on_vert * num_nodes;
    nodes_indices.resize(max_nb_dofs, false);
    local_nodes_indices.resize(max_nb_dofs, false);
    if (std::distance(dit, hi_dit) != max_nb_dofs) {
      std::fill(nodes_indices.begin(), nodes_indices.end(), -1);
      std::fill(local_nodes_indices.begin(), local_nodes_indices.end(), -1);
    }

    for (; dit != hi_dit; dit++) {
      auto &dof = **dit;
      const int idx = dof.getPetscGlobalDofIdx();
      const int local_idx = dof.getPetscLocalDofIdx();
      const int side_number = dof.getSideNumberPtr()->side_number;
      const int pos = side_number * nb_dofs_on_vert + dof.getDofCoeffIdx();
      nodes_indices[pos] = idx;
      local_nodes_indices[pos] = local_idx;
      const int brother_side_number =
          (*dit)->getSideNumberPtr()->brother_side_number;
      if (brother_side_number != -1) {
        const int elem_idx =
            brother_side_number * nb_dofs_on_vert + (*dit)->getDofCoeffIdx();
        nodes_indices[elem_idx] = idx;
        local_nodes_indices[elem_idx] = local_idx;
      }
    }

  } else {
    nodes_indices.resize(0, false);
    local_nodes_indices.resize(0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getRowNodesIndices(DataForcesAndSourcesCore &data,
                                         const std::string &field_name) const {
  return getNodesIndices(field_name,
                         const_cast<FENumeredDofEntity_multiIndex &>(
                             numeredEntFiniteElementPtr->getRowDofs()),
                         data.dataOnEntities[MBVERTEX][0].getIndices(),
                         data.dataOnEntities[MBVERTEX][0].getLocalIndices());
}

MoFEMErrorCode
ForcesAndSourcesCore::getColNodesIndices(DataForcesAndSourcesCore &data,
                                         const std::string &field_name) const {
  return getNodesIndices(field_name,
                         const_cast<FENumeredDofEntity_multiIndex &>(
                             numeredEntFiniteElementPtr->getColDofs()),
                         data.dataOnEntities[MBVERTEX][0].getIndices(),
                         data.dataOnEntities[MBVERTEX][0].getLocalIndices());
}

MoFEMErrorCode ForcesAndSourcesCore::getEntityIndices(
    DataForcesAndSourcesCore &data, const std::string &field_name,
    FENumeredDofEntity_multiIndex &dofs, const EntityType type_lo,
    const EntityType type_hi) const {
  MoFEMFunctionBegin;

  for (EntityType t = type_lo; t != type_hi; ++t) {
    for (auto &dat : data.dataOnEntities[t]) {
      dat.getIndices().resize(0, false);
      dat.getLocalIndices().resize(0, false);
    }
  }

  auto &dofs_by_uid = dofs.get<Unique_mi_tag>();
  auto bit_number = mField.get_field_bit_number(field_name);

  auto dit = dofs_by_uid.lower_bound(
      DofEntity::getLoFieldEntityUId(bit_number, get_id_for_min_type(type_lo)));
  if (dit == dofs_by_uid.end())
    MoFEMFunctionReturnHot(0);
  auto hi_dit = dofs_by_uid.upper_bound(
      DofEntity::getHiFieldEntityUId(bit_number, get_id_for_max_type(type_hi)));

  for (; dit != hi_dit; ++dit) {
    auto &dof = **dit;
    const EntityType type = dof.getEntType();
    const int side = dof.getSideNumberPtr()->side_number;

    if (side >= 0) {

      auto &dat = data.dataOnEntities[type][side];
      const int nb_dofs_on_ent = dof.getNbDofsOnEnt();
      if (nb_dofs_on_ent) {
        const int brother_side = dof.getSideNumberPtr()->brother_side_number;
        auto &ent_field_indices = dat.getIndices();
        auto &ent_field_local_indices = dat.getLocalIndices();
        if (ent_field_indices.empty()) {
          ent_field_indices.resize(nb_dofs_on_ent, false);
          ent_field_local_indices.resize(nb_dofs_on_ent, false);
          std::fill(ent_field_indices.data().begin(),
                    ent_field_indices.data().end(), -1);
          std::fill(ent_field_local_indices.data().begin(),
                    ent_field_local_indices.data().end(), -1);
        }
        const int idx = dof.getEntDofIdx();
        ent_field_indices[idx] = dof.getPetscGlobalDofIdx();
        ent_field_local_indices[idx] = dof.getPetscLocalDofIdx();
        if (brother_side != -1) {
          auto &dat_brother = data.dataOnEntities[type][brother_side];
          dat_brother.getIndices().resize(nb_dofs_on_ent, false);
          dat_brother.getLocalIndices().resize(nb_dofs_on_ent, false);
          dat_brother.getIndices()[idx] = dat.getIndices()[idx];
          dat_brother.getLocalIndices()[idx] = dat.getLocalIndices()[idx];
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getNoFieldIndices(const std::string &field_name,
                                        FENumeredDofEntity_multiIndex &dofs,
                                        VectorInt &indices) const {
  MoFEMFunctionBeginHot;
  auto field_it = fieldsPtr->get<FieldName_mi_tag>().find(field_name);
  auto dit = dofs.get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId((*field_it)->getBitNumber()));
  auto hi_dit = dofs.get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId((*field_it)->getBitNumber()));
  indices.resize(std::distance(dit, hi_dit));
  for (; dit != hi_dit; dit++) {
    int idx = (*dit)->getPetscGlobalDofIdx();
    indices[(*dit)->getDofCoeffIdx()] = idx;
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ForcesAndSourcesCore::getNoFieldRowIndices(
    DataForcesAndSourcesCore &data, const std::string &field_name) const {
  MoFEMFunctionBegin;
  if (data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  CHKERR getNoFieldIndices(field_name,
                           const_cast<FENumeredDofEntity_multiIndex &>(
                               numeredEntFiniteElementPtr->getRowDofs()),
                           data.dataOnEntities[MBENTITYSET][0].getIndices());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ForcesAndSourcesCore::getNoFieldColIndices(
    DataForcesAndSourcesCore &data, const std::string &field_name) const {
  MoFEMFunctionBegin;
  if (data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  CHKERR getNoFieldIndices(field_name,
                           const_cast<FENumeredDofEntity_multiIndex &>(
                               numeredEntFiniteElementPtr->getColDofs()),
                           data.dataOnEntities[MBENTITYSET][0].getIndices());
  MoFEMFunctionReturn(0);
}

// ** Indices from problem **

MoFEMErrorCode ForcesAndSourcesCore::getProblemNodesIndices(
    const std::string &field_name, const NumeredDofEntity_multiIndex &dofs,
    VectorInt &nodes_indices) const {
  MoFEMFunctionBeginHot;

  const Field *field_struture = mField.get_field_structure(field_name);
  if (field_struture->getSpace() == H1) {

    int num_nodes;
    CHKERR getNumberOfNodes(num_nodes);
    nodes_indices.resize(field_struture->getNbOfCoeffs() * num_nodes, false);
    std::fill(nodes_indices.begin(), nodes_indices.end(), -1);

    auto &side_table = const_cast<SideNumber_multiIndex &>(
        numeredEntFiniteElementPtr->getSideNumberTable());
    auto siit = side_table.get<1>().lower_bound(boost::make_tuple(MBVERTEX, 0));
    auto hi_siit =
        side_table.get<1>().lower_bound(boost::make_tuple(MBVERTEX, 10000));

    int nn = 0;
    for (; siit != hi_siit; siit++, nn++) {

      if (siit->get()->side_number == -1)
        continue;

      auto bit_number = mField.get_field_bit_number(field_name);
      const EntityHandle ent = siit->get()->ent;
      auto dit = dofs.get<Unique_mi_tag>().lower_bound(
          FieldEntity::getLoLocalEntityBitNumber(bit_number, ent));
      auto hi_dit = dofs.get<Unique_mi_tag>().upper_bound(
          FieldEntity::getHiLocalEntityBitNumber(bit_number, ent));
      for (; dit != hi_dit; dit++) {
        nodes_indices[siit->get()->side_number * (*dit)->getNbOfCoeffs() +
                      (*dit)->getDofCoeffIdx()] =
            (*dit)->getPetscGlobalDofIdx();
      }
    }
  } else {
    nodes_indices.resize(0, false);
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ForcesAndSourcesCore::getProblemTypeIndices(
    const std::string &field_name, const NumeredDofEntity_multiIndex &dofs,
    EntityType type, int side_number, VectorInt &indices) const {
  MoFEMFunctionBeginHot;

  indices.resize(0);

  auto &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());
  auto siit =
      side_table.get<1>().lower_bound(boost::make_tuple(type, side_number));
  auto hi_siit =
      side_table.get<1>().upper_bound(boost::make_tuple(type, side_number));

  for (; siit != hi_siit; siit++) {

    if (siit->get()->side_number == -1)
      continue;

    const EntityHandle ent = siit->get()->ent;
    auto bit_number = mField.get_field_bit_number(field_name);
    auto dit = dofs.get<Unique_mi_tag>().lower_bound(
        FieldEntity::getLoLocalEntityBitNumber(bit_number, ent));
    auto hi_dit = dofs.get<Unique_mi_tag>().upper_bound(
        FieldEntity::getHiLocalEntityBitNumber(bit_number, ent));
    indices.resize(std::distance(dit, hi_dit));
    for (; dit != hi_dit; dit++) {

      indices[(*dit)->getEntDofIdx()] = (*dit)->getPetscGlobalDofIdx();
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ForcesAndSourcesCore::getProblemNodesRowIndices(
    const std::string &field_name, VectorInt &nodes_indices) const {
  return getProblemNodesIndices(field_name, *(problemPtr->numeredDofsRows),
                                nodes_indices);
}

MoFEMErrorCode
ForcesAndSourcesCore::getProblemTypeRowIndices(const std::string &field_name,
                                               EntityType type, int side_number,
                                               VectorInt &indices) const {
  return getProblemTypeIndices(field_name, *(problemPtr->numeredDofsRows), type,
                               side_number, indices);
}

MoFEMErrorCode ForcesAndSourcesCore::getProblemNodesColIndices(
    const std::string &field_name, VectorInt &nodes_indices) const {
  return getProblemNodesIndices(field_name, *(problemPtr->numeredDofsCols),
                                nodes_indices);
}

MoFEMErrorCode
ForcesAndSourcesCore::getProblemTypeColIndices(const std::string &field_name,
                                               EntityType type, int side_number,
                                               VectorInt &indices) const {
  return getProblemTypeIndices(field_name, *(problemPtr->numeredDofsCols), type,
                               side_number, indices);
}

// ** Data **

MoFEMErrorCode
ForcesAndSourcesCore::getNodesFieldData(DataForcesAndSourcesCore &data,
                                        const std::string &field_name) const {

  auto get_nodes_field_data =
      [&](FEDofEntity_multiIndex &dofs, VectorDouble &nodes_data,
          VectorDofs &nodes_dofs, FieldSpace &space,
          FieldApproximationBase &base, VectorInt &bb_node_order) {
        MoFEMFunctionBegin;

        auto &dofs_by_uid = dofs.get<Unique_mi_tag>();
        auto bit_number = mField.get_field_bit_number(field_name);
        auto dit = dofs_by_uid.lower_bound(DofEntity::getLoFieldEntityUId(
            bit_number, get_id_for_min_type<MBVERTEX>()));
        decltype(dit) hi_dit;
        if (dit == dofs_by_uid.end())
          hi_dit = dit;
        else
          hi_dit = dofs_by_uid.upper_bound(DofEntity::getHiFieldEntityUId(
              bit_number, get_id_for_max_type<MBVERTEX>()));

        if (dit != hi_dit) {
          auto &first_dof = **dit;
          space = first_dof.getSpace();
          base = first_dof.getApproxBase();
          int num_nodes;
          CHKERR getNumberOfNodes(num_nodes);
          bb_node_order.resize(num_nodes, false);
          bb_node_order.clear();
          const int nb_dof_idx = first_dof.getNbOfCoeffs();
          const int max_nb_dofs = nb_dof_idx * num_nodes;
          nodes_data.resize(max_nb_dofs, false);
          nodes_dofs.resize(max_nb_dofs, false);
          nodes_data.clear();

          std::vector<boost::weak_ptr<FEDofEntity>> brother_dofs_vec;
          for (; dit != hi_dit;) {
            const auto &dof_ptr = *dit;
            const auto &dof = *dof_ptr;
            const auto &sn = *dof.getSideNumberPtr();
            const int side_number = sn.side_number;
            const int brother_side_number = sn.brother_side_number;
            if (brother_side_number != -1)
              brother_dofs_vec.emplace_back(dof_ptr);

            bb_node_order[side_number] = dof.getMaxOrder();
            int pos = side_number * nb_dof_idx;
            auto ent_filed_data_vec = dof.getEntFieldData();
            for (int ii = 0; ii != nb_dof_idx; ++ii) {
              nodes_data[pos] = ent_filed_data_vec[ii];
              nodes_dofs[pos] = (*dit).get();
              ++pos;
              ++dit;
            }
          }

          for (auto &dof_ptr : brother_dofs_vec) {
            if (const auto d = dof_ptr.lock()) {
              const auto &sn = d->getSideNumberPtr();
              const int side_number = sn->side_number;
              const int brother_side_number = sn->brother_side_number;
              bb_node_order[brother_side_number] = bb_node_order[side_number];
              int pos = side_number * nb_dof_idx;
              int brother_pos = brother_side_number * nb_dof_idx;
              for (int ii = 0; ii != nb_dof_idx; ++ii) {
                nodes_data[brother_pos] = nodes_data[pos];
                nodes_dofs[brother_pos] = nodes_dofs[pos];
                ++pos;
                ++brother_pos;
              }
            }
          }

        } else {
          nodes_data.resize(0, false);
          nodes_dofs.resize(0, false);
        }

        MoFEMFunctionReturn(0);
      };

  return get_nodes_field_data(
      const_cast<FEDofEntity_multiIndex &>(getDataDofs()),
      data.dataOnEntities[MBVERTEX][0].getFieldData(),
      data.dataOnEntities[MBVERTEX][0].getFieldDofs(),
      data.dataOnEntities[MBVERTEX][0].getSpace(),
      data.dataOnEntities[MBVERTEX][0].getBase(),
      data.dataOnEntities[MBVERTEX][0].getBBNodeOrder());
}

MoFEMErrorCode ForcesAndSourcesCore::getEntityFieldData(
    DataForcesAndSourcesCore &data, const std::string &field_name,
    const EntityType type_lo, const EntityType type_hi) const {
  MoFEMFunctionBegin;
  for (EntityType t = type_lo; t != type_hi; ++t) {
    for (auto &dat : data.dataOnEntities[t]) {
      dat.getDataOrder() = 0;
      dat.getBase() = NOBASE;
      dat.getSpace() = NOSPACE;
      dat.getFieldData().resize(0, false);
      dat.getFieldDofs().resize(0, false);
    }
  }

  auto &dofs = const_cast<FEDofEntity_multiIndex &>(
      numeredEntFiniteElementPtr->getDataDofs());
  auto bit_number = mField.get_field_bit_number(field_name);
  auto &dofs_by_uid = dofs.get<Unique_mi_tag>();
  auto dit = dofs_by_uid.lower_bound(
      DofEntity::getLoFieldEntityUId(bit_number, get_id_for_min_type(type_lo)));
  if (dit == dofs_by_uid.end())
    MoFEMFunctionReturnHot(0);
  auto hi_dit = dofs_by_uid.upper_bound(
      DofEntity::getHiFieldEntityUId(bit_number, get_id_for_max_type(type_hi)));

  std::vector<boost::weak_ptr<FEDofEntity>> brother_dofs_vec;
  for (; dit != hi_dit;) {

    auto &dof = **dit;
    const int nb_dofs_on_ent = dof.getNbDofsOnEnt();
    if (nb_dofs_on_ent) {

      const EntityType type = dof.getEntType();
      const int side = dof.getSideNumberPtr()->side_number;
      if (side >= 0) {

        auto &dat = data.dataOnEntities[type][side];
        auto &ent_field_dofs = dat.getFieldDofs();
        auto &ent_field_data = dat.getFieldData();
        const int brother_side = dof.getSideNumberPtr()->brother_side_number;
        if (brother_side != -1)
          brother_dofs_vec.emplace_back(*dit);

        if (ent_field_data.empty()) {

          dat.getBase() = dof.getApproxBase();
          dat.getSpace() = dof.getSpace();
          const int ent_order = dof.getMaxOrder();
          dat.getDataOrder() =
              dat.getDataOrder() > ent_order ? dat.getDataOrder() : ent_order;
          ent_field_data.resize(nb_dofs_on_ent, false);
          noalias(ent_field_data) = dof.getEntFieldData();
          ent_field_dofs.resize(nb_dofs_on_ent, false);
          for (int ii = 0; ii != nb_dofs_on_ent; ++ii) {
            ent_field_dofs[(*dit)->getEntDofIdx()] = (*dit).get();
            ++dit;
          }
        }

      } else {

        for (int ii = 0; ii != nb_dofs_on_ent; ++ii)
          ++dit;
      }
    }
  }

  for (auto &dof_ptr : brother_dofs_vec) {
    if (auto d = dof_ptr.lock()) {
      const EntityType type = d->getEntType();
      const int side = d->getSideNumberPtr()->side_number;
      const int brother_side = d->getSideNumberPtr()->brother_side_number;
      auto &dat = data.dataOnEntities[type][side];
      auto &dat_brother = data.dataOnEntities[type][brother_side];
      dat_brother.getBase() = dat.getBase();
      dat_brother.getSpace() = dat.getSpace();
      dat_brother.getDataOrder() = dat.getDataOrder();
      dat_brother.getFieldData() = dat.getFieldData();
      dat_brother.getFieldDofs() = dat.getFieldDofs();
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ForcesAndSourcesCore::getNoFieldFieldData(
    const boost::string_ref field_name, FEDofEntity_multiIndex &dofs,
    VectorDouble &ent_field_data, VectorDofs &ent_field_dofs) const {
  MoFEMFunctionBeginHot;
  auto field_it = fieldsPtr->get<FieldName_mi_tag>().find(field_name);

  auto dit = dofs.get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId((*field_it)->getBitNumber()));
  auto hi_dit = dofs.get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId((*field_it)->getBitNumber()));

  int size = std::distance(dit, hi_dit);
  ent_field_data.resize(size, false);
  ent_field_dofs.resize(size, false);
  for (; dit != hi_dit; dit++) {
    int idx = (*dit)->getDofCoeffIdx();
    ent_field_data[idx] = (*dit)->getFieldData();
    ent_field_dofs[idx] = (*dit).get();
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ForcesAndSourcesCore::getNoFieldFieldData(
    DataForcesAndSourcesCore &data, const boost::string_ref field_name) const {
  MoFEMFunctionBegin;
  if (data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  CHKERR getNoFieldFieldData(
      field_name,
      const_cast<FEDofEntity_multiIndex &>(
          numeredEntFiniteElementPtr->getDataDofs()),
      data.dataOnEntities[MBENTITYSET][0].getFieldData(),
      data.dataOnEntities[MBENTITYSET][0].getFieldDofs());
  MoFEMFunctionReturn(0);
}

// ** Face **

MoFEMErrorCode
ForcesAndSourcesCore::getFaceTriNodes(DataForcesAndSourcesCore &data) const {
  MoFEMFunctionBegin;
  // PetscAttachDebugger();
  data.facesNodes.resize(4, 3, false);
  auto &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());
  auto siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI, 0));
  auto hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI, 4));
  if (std::distance(siit, hi_siit) != 4) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Should be 4 triangles on tet, side_table not initialized");
  }
  const int canonical_face_sense_p1[4][3] = {
      {0, 1, 3},
      {1, 2, 3},
      {0, 3, 2} /**/,
      {0, 2, 1} /**/}; // second index is offset (positive sense)
  const int canonical_face_sense_m1[4][3] = {
      {0, 3, 1},
      {1, 3, 2},
      {0, 2, 3},
      {0, 1, 2}}; // second index is offset (negative sense
  for (; siit != hi_siit; siit++) {
    const boost::shared_ptr<SideNumber> side = *siit;
    int face_conn[3] = {-1, -1, -1};
    if (side->offset == 0) {
      face_conn[0] = side->sense == 1
                         ? canonical_face_sense_p1[(int)side->side_number][0]
                         : canonical_face_sense_m1[(int)side->side_number][0];
      face_conn[1] = side->sense == 1
                         ? canonical_face_sense_p1[(int)side->side_number][1]
                         : canonical_face_sense_m1[(int)side->side_number][1];
      face_conn[2] = side->sense == 1
                         ? canonical_face_sense_p1[(int)side->side_number][2]
                         : canonical_face_sense_m1[(int)side->side_number][2];
    }
    if (side->offset == 1) {
      face_conn[0] =
          side->sense == 1
              ? canonical_face_sense_p1[(int)side->side_number][1]
              : canonical_face_sense_m1[(int)side->side_number][2] /**/;
      face_conn[1] = side->sense == 1
                         ? canonical_face_sense_p1[(int)side->side_number][2]
                         : canonical_face_sense_m1[(int)side->side_number][0];
      face_conn[2] = side->sense == 1
                         ? canonical_face_sense_p1[(int)side->side_number][0]
                         : canonical_face_sense_m1[(int)side->side_number][1];
    }
    if (side->offset == 2) {
      face_conn[0] =
          side->sense == 1
              ? canonical_face_sense_p1[(int)side->side_number][2]
              : canonical_face_sense_m1[(int)side->side_number][1] /**/;
      face_conn[1] = side->sense == 1
                         ? canonical_face_sense_p1[(int)side->side_number][0]
                         : canonical_face_sense_m1[(int)side->side_number][2];
      face_conn[2] = side->sense == 1
                         ? canonical_face_sense_p1[(int)side->side_number][1]
                         : canonical_face_sense_m1[(int)side->side_number][0];
    }
    for (int nn = 0; nn < 3; nn++)
      data.facesNodes(side->side_number, nn) = face_conn[nn];
    {
      const EntityHandle *conn_tet;
      int num_nodes_tet;
      EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
      CHKERR mField.get_moab().get_connectivity(ent, conn_tet, num_nodes_tet,
                                                true);
      if (num_nodes_tet != 4)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      int num_nodes_face;
      const EntityHandle *conn_face;
      CHKERR mField.get_moab().get_connectivity(side->ent, conn_face,
                                                num_nodes_face, true);
      if (num_nodes_face != 3)
        SETERRQ(PETSC_COMM_SELF, 1, "data inconsistency");
      if (conn_face[0] != conn_tet[data.facesNodes(side->side_number, 0)])
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      if (conn_face[1] != conn_tet[data.facesNodes(side->side_number, 1)])
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      if (conn_face[2] != conn_tet[data.facesNodes(side->side_number, 2)])
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
    }
  }
  MoFEMFunctionReturn(0);
}

// ** Space and Base **

MoFEMErrorCode ForcesAndSourcesCore::getSpacesAndBaseOnEntities(
    DataForcesAndSourcesCore &data) const {
  MoFEMFunctionBeginHot;

  if (nInTheLoop == 0) {
    data.sPace.reset();
    data.bAse.reset();
    for (EntityType t = MBVERTEX; t != MBMAXTYPE; ++t) {
      data.spacesOnEntities[t].reset();
      data.basesOnEntities[t].reset();
    }
    for (int s = 0; s != LASTSPACE; ++s) {
      data.basesOnSpaces[s].reset();
    }
  }

  if (getDataFieldEntsPtr())
    for (auto e : getDataFieldEnts()) {
      // get data from entity
      const EntityType type = e->getEntType();
      const FieldSpace space = e->getSpace();
      const FieldApproximationBase approx = e->getApproxBase();

      // set data
      data.sPace.set(space);
      data.bAse.set(approx);
      data.spacesOnEntities[type].set(space);
      data.basesOnEntities[type].set(approx);
      data.basesOnSpaces[space].set(approx);
    }
  else
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "data fields ents not allocated on element");

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ForcesAndSourcesCore::calHierarchicalBaseFunctionsOnElement(
    const FieldApproximationBase b) {
  MoFEMFunctionBegin;
  if (dataOnElement[H1]->bAse.test(b)) {
    switch (static_cast<FieldApproximationBase>(b)) {
    case NOBASE:
      break;
    case AINSWORTH_BERNSTEIN_BEZIER_BASE:
      // BERNSTEIN_BEZIER_BASE is not hierarchical base
      break;
    case AINSWORTH_LEGENDRE_BASE:
    case AINSWORTH_LOBATTO_BASE:
    case DEMKOWICZ_JACOBI_BASE:
      if (!getElementPolynomialBase())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Functions genrating approximation base not defined");

      for (int space = H1; space != LASTSPACE; ++space) {
        if (dataOnElement[H1]->sPace.test(space) &&
            dataOnElement[H1]->bAse.test(b) &&
            dataOnElement[H1]->basesOnSpaces[space].test(b)) {
          CHKERR getElementPolynomialBase()->getValue(
              gaussPts,
              boost::make_shared<EntPolynomialBaseCtx>(
                  *dataOnElement[space], static_cast<FieldSpace>(space),
                  static_cast<FieldApproximationBase>(b), NOBASE));
        }
      }
      break;
    case USER_BASE:
      if (!getUserPolynomialBase())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Functions genrating user approximation base not defined");

      for (int space = H1; space != LASTSPACE; ++space)
        if (dataOnElement[H1]->sPace.test(space) &&
            dataOnElement[H1]->bAse.test(b) &&
            dataOnElement[H1]->basesOnSpaces[space].test(b)) {
          CHKERR getUserPolynomialBase()->getValue(
              gaussPts,
              boost::make_shared<EntPolynomialBaseCtx>(
                  *dataOnElement[space], static_cast<FieldSpace>(space),
                  static_cast<FieldApproximationBase>(b), NOBASE));
        }
      break;
    default:
      SETERRQ1(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "Base <%s> not yet implemented",
               ApproximationBaseNames[static_cast<FieldApproximationBase>(b)]);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ForcesAndSourcesCore::calHierarchicalBaseFunctionsOnElement() {
  MoFEMFunctionBegin;
  /// Use the some node base. Node base is usually used for construction other
  /// bases.
  for (int space = HCURL; space != LASTSPACE; ++space) {
    dataOnElement[space]->dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) =
        dataOnElement[H1]->dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
  }
  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
    CHKERR calHierarchicalBaseFunctionsOnElement(
        static_cast<FieldApproximationBase>(b));
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::calBernsteinBezierBaseFunctionsOnElement() {
  MoFEMFunctionBegin;

  auto get_nodal_base_data = [&](DataForcesAndSourcesCore &data,
                                 const std::string &field_name) {
    MoFEMFunctionBegin;
    auto &space = data.dataOnEntities[MBVERTEX][0].getSpace();
    auto &base = data.dataOnEntities[MBVERTEX][0].getBase();
    auto &bb_node_order = data.dataOnEntities[MBVERTEX][0].getBBNodeOrder();

    auto bit_number = mField.get_field_bit_number(field_name);
    auto &dofs_by_uid = getDataDofs().get<Unique_mi_tag>();
    auto dit = dofs_by_uid.lower_bound(DofEntity::getLoFieldEntityUId(
        bit_number, get_id_for_min_type<MBVERTEX>()));
    if (dit == dofs_by_uid.end())
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "No nodal dofs on element");
    auto hi_dit = dofs_by_uid.upper_bound(DofEntity::getHiFieldEntityUId(
        bit_number, get_id_for_max_type<MBVERTEX>()));

    if (dit != hi_dit) {
      auto &first_dof = **dit;
      space = first_dof.getSpace();
      base = first_dof.getApproxBase();
      int num_nodes;
      CHKERR getNumberOfNodes(num_nodes);
      bb_node_order.resize(num_nodes, false);
      bb_node_order.clear();
      const int nb_dof_idx = first_dof.getNbOfCoeffs();

      std::vector<boost::weak_ptr<FEDofEntity>> brother_dofs_vec;
      for (; dit != hi_dit;) {
        const auto &dof_ptr = *dit;
        const auto &dof = *dof_ptr;
        const auto &sn = *dof.getSideNumberPtr();
        const int side_number = sn.side_number;
        const int brother_side_number = sn.brother_side_number;
        if (brother_side_number != -1)
          brother_dofs_vec.emplace_back(dof_ptr);
        bb_node_order[side_number] = dof.getMaxOrder();
        for (int ii = 0; ii != nb_dof_idx; ++ii)
          ++dit;
      }

      for (auto &dof_ptr : brother_dofs_vec) {
        if (const auto d = dof_ptr.lock()) {
          const auto &sn = d->getSideNumberPtr();
          const int side_number = sn->side_number;
          const int brother_side_number = sn->brother_side_number;
          bb_node_order[brother_side_number] = bb_node_order[side_number];
        }
      }
    }
    MoFEMFunctionReturn(0);
  };

  auto get_entity_base_data = [&](DataForcesAndSourcesCore &data,
                                  const std::string &field_name,
                                  const EntityType type_lo,
                                  const EntityType type_hi) {
    MoFEMFunctionBegin;
    for (EntityType t = type_lo; t != type_hi; ++t) {
      for (auto &dat : data.dataOnEntities[t]) {
        dat.getDataOrder() = 0;
        dat.getBase() = NOBASE;
        dat.getSpace() = NOSPACE;
        dat.getFieldData().resize(0, false);
        dat.getFieldDofs().resize(0, false);
      }
    }

    auto &dofs = const_cast<FEDofEntity_multiIndex &>(
        numeredEntFiniteElementPtr->getDataDofs());
    auto &dofs_by_uid = dofs.get<Unique_mi_tag>();
    auto bit_number = mField.get_field_bit_number(field_name);
    auto dit = dofs_by_uid.lower_bound(DofEntity::getLoFieldEntityUId(
        bit_number, get_id_for_min_type(type_lo)));
    if (dit == dofs_by_uid.end())
      MoFEMFunctionReturnHot(0);
    auto hi_dit = dofs_by_uid.upper_bound(DofEntity::getHiFieldEntityUId(
        bit_number, get_id_for_max_type(type_hi)));

    std::vector<boost::weak_ptr<FEDofEntity>> brother_dofs_vec;
    for (; dit != hi_dit;) {

      auto &dof = **dit;
      const int nb_dofs_on_ent = dof.getNbDofsOnEnt();
      if (nb_dofs_on_ent) {
        const EntityType type = dof.getEntType();
        const int side = dof.getSideNumberPtr()->side_number;
        auto &dat = data.dataOnEntities[type][side];

        const int brother_side = dof.getSideNumberPtr()->brother_side_number;
        if (brother_side != -1)
          brother_dofs_vec.emplace_back(*dit);

        dat.getBase() = dof.getApproxBase();
        dat.getSpace() = dof.getSpace();
        const int ent_order = dof.getMaxOrder();
        dat.getDataOrder() =
            dat.getDataOrder() > ent_order ? dat.getDataOrder() : ent_order;
        for (int ii = 0; ii != nb_dofs_on_ent; ++ii) {
          ++dit;
        }
      }
    }

    for (auto &dof_ptr : brother_dofs_vec) {
      if (auto d = dof_ptr.lock()) {
        const EntityType type = d->getEntType();
        const int side = d->getSideNumberPtr()->side_number;
        const int brother_side = d->getSideNumberPtr()->brother_side_number;
        auto &dat = data.dataOnEntities[type][side];
        auto &dat_brother = data.dataOnEntities[type][brother_side];
        dat_brother.getBase() = dat.getBase();
        dat_brother.getSpace() = dat.getSpace();
        dat_brother.getDataOrder() = dat.getDataOrder();
      }
    }
    MoFEMFunctionReturn(0);
  };

  for (auto &e : getDataFieldEnts()) {
    if (e->getApproxBase() == AINSWORTH_BERNSTEIN_BEZIER_BASE) {
      auto space = e->getSpace();
      for (EntityType t = MBVERTEX; t != MBPOLYHEDRON; ++t) {
        for (auto &dat : (*dataOnElement[space]).dataOnEntities[t]) {
          for (auto &ptr : dat.getBBAlphaIndicesByOrderArray())
            ptr.reset();
          for (auto &ptr : dat.getBBNByOrderArray())
            ptr.reset();
          for (auto &ptr : dat.getBBDiffNByOrderArray())
            ptr.reset();
        }
      }
    }
  }

  for (auto &e : getDataFieldEnts()) {
    if (e->getApproxBase() == AINSWORTH_BERNSTEIN_BEZIER_BASE) {
      auto field_name = e->getName();
      auto space = e->getSpace();
      CHKERR get_nodal_base_data(*dataOnElement[space], field_name);
      CHKERR get_entity_base_data(*dataOnElement[space], field_name, MBEDGE,
                                  MBPOLYHEDRON);
      CHKERR getElementPolynomialBase()->getValue(
          gaussPts,
          boost::make_shared<EntPolynomialBaseCtx>(
              *dataOnElement[space], field_name, static_cast<FieldSpace>(space),
              AINSWORTH_BERNSTEIN_BEZIER_BASE, NOBASE));
    }
  }

  MoFEMFunctionReturn(0);
};

MoFEMErrorCode ForcesAndSourcesCore::createDataOnElement() {
  MoFEMFunctionBegin;

  const EntityType type = numeredEntFiniteElementPtr->getEntType();
  if (type == lastEvaluatedElementEntityType)
    MoFEMFunctionReturnHot(0);

  // Data on elements for proper spaces
  for (int space = H1; space != LASTSPACE; ++space) {
    dataOnElement[space]->setElementType(type);
    derivedDataOnElement[space]->setElementType(type);
  }

  lastEvaluatedElementEntityType = type;

  MoFEMFunctionReturn(0);
}

#define FUNCTION_NAME_WITH_OP_NAME(OP)                                         \
  std::ostringstream ss;                                                       \
  ss << "(Calling user data operator "                                         \
     << boost::typeindex::type_id_runtime(OP).pretty_name() << " rowField "    \
     << (OP).rowFieldName << " colField " << (OP).colFieldName << ") "

#define CATCH_OP_ERRORS(OP)                                                    \
  catch (MoFEMExceptionInitial const &ex) {                                    \
    FUNCTION_NAME_WITH_OP_NAME(OP) << PETSC_FUNCTION_NAME;                     \
    return PetscError(PETSC_COMM_SELF, ex.lINE, ss.str().c_str(), __FILE__,    \
                      ex.errorCode, PETSC_ERROR_INITIAL, ex.what());           \
  }                                                                            \
  catch (MoFEMExceptionRepeat const &ex) {                                     \
    FUNCTION_NAME_WITH_OP_NAME(OP) << PETSC_FUNCTION_NAME;                     \
    return PetscError(PETSC_COMM_SELF, ex.lINE, ss.str().c_str(), __FILE__,    \
                      ex.errorCode, PETSC_ERROR_REPEAT, " ");                  \
  }                                                                            \
  catch (MoFEMException const &ex) {                                           \
    FUNCTION_NAME_WITH_OP_NAME(OP) << ex.errorMessage;                         \
    SETERRQ(PETSC_COMM_SELF, ex.errorCode, ss.str().c_str());                  \
  }                                                                            \
  catch (std::exception const &ex) {                                           \
    std::string message("Error: " + std::string(ex.what()) + " at " +          \
                        boost::lexical_cast<std::string>(__LINE__) + " : " +   \
                        std::string(__FILE__) + " in " +                       \
                        std::string(PETSC_FUNCTION_NAME));                     \
    FUNCTION_NAME_WITH_OP_NAME(OP) << message;                                 \
    SETERRQ(PETSC_COMM_SELF, MOFEM_STD_EXCEPTION_THROW, ss.str().c_str());     \
  }

MoFEMErrorCode ForcesAndSourcesCore::loopOverOperators() {
  MoFEMFunctionBegin;

  const EntityType type = numeredEntFiniteElementPtr->getEntType();
  const UserDataOperator::OpType types[2] = {UserDataOperator::OPROW,
                                             UserDataOperator::OPCOL};
  std::vector<std::string> last_eval_field_name(2);

  boost::ptr_vector<UserDataOperator>::iterator oit, hi_oit;
  oit = opPtrVector.begin();
  hi_oit = opPtrVector.end();

  for (; oit != hi_oit; oit++) {

    try {

      CHKERR oit->setPtrFE(this);

      if (oit->opType == UserDataOperator::OPLAST) {

        // Set field
        switch (oit->sPace) {
        case NOSPACE:
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Unknown space");
        case NOFIELD:
        case H1:
        case HCURL:
        case HDIV:
        case L2:
          break;
        default:
          SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                   "Not implemented for this space", oit->sPace);
        }

        // Reseat all data which all field dependent
        dataOnElement[oit->sPace]->resetFieldDependentData();
        last_eval_field_name[0] = "";

        // Run operator
        try {
          CHKERR oit->opRhs(*dataOnElement[oit->sPace], false);
        }
        CATCH_OP_ERRORS(*oit);

      } else {

        boost::shared_ptr<DataForcesAndSourcesCore> op_data[2];
        std::array<bool, 2> base_swap;
        std::array<std::pair<std::string, FieldApproximationBase>, 2>
            base_swap_data;
        auto swap_bases = [&]() {
          MoFEMFunctionBeginHot;
          for (size_t ss = 0; ss != 2; ++ss)
            if (base_swap[ss])
              CHKERR op_data[ss]->baseSwap(base_swap_data[ss].first,
                                           base_swap_data[ss].second);
          MoFEMFunctionReturnHot(0);
        };

        for (size_t ss = 0; ss != 2; ss++) {

          const std::string field_name =
              !ss ? oit->rowFieldName : oit->colFieldName;
          if (field_name.empty()) {
            SETERRQ2(
                PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "No field name in operator %d (0-row, 1-column) in operator %s",
                ss,
                (boost::typeindex::type_id_runtime(*oit).pretty_name())
                    .c_str());
          }
          const Field *field_struture = mField.get_field_structure(field_name);
          const BitFieldId data_id = field_struture->getId();
          const FieldSpace space = field_struture->getSpace();
          const FieldApproximationBase base = field_struture->getApproxBase();
          op_data[ss] =
              !ss ? dataOnElement[space] : derivedDataOnElement[space];

          switch (base) {
          case AINSWORTH_BERNSTEIN_BEZIER_BASE:
            base_swap_data[ss] = std::pair<std::string, FieldApproximationBase>(
                field_name, AINSWORTH_BERNSTEIN_BEZIER_BASE);
            base_swap[ss] = true;
            break;
          default:
            base_swap[ss] = false;
          };

          if ((oit->getNumeredEntFiniteElementPtr()->getBitFieldIdData() &
               data_id)
                  .none()) {
            SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                     "no data field < %s > on finite element < %s >",
                     field_name.c_str(), feName.c_str());
          }

          if (oit->getOpType() & types[ss] ||
              oit->getOpType() & UserDataOperator::OPROWCOL) {

            switch (space) {
            case NOSPACE:
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                      "unknown space");
              break;
            case NOFIELD:
            case H1:
            case HCURL:
            case HDIV:
            case L2:
              break;
            default:
              SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                       "Not implemented for this space", space);
            }

            if (last_eval_field_name[ss] != field_name) {

              CHKERR getEntityFieldData(*op_data[ss], field_name, MBEDGE);
              if (!ss)
                CHKERR getEntityRowIndices(*op_data[ss], field_name, MBEDGE);
              else
                CHKERR getEntityColIndices(*op_data[ss], field_name, MBEDGE);

              switch (space) {
              case NOSPACE:
                SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                        "unknown space");
                break;
              case H1:
                if (!ss)
                  CHKERR getRowNodesIndices(*op_data[ss], field_name);
                else
                  CHKERR getColNodesIndices(*op_data[ss], field_name);
                CHKERR getNodesFieldData(*op_data[ss], field_name);
                break;
              case HCURL:
              case HDIV:
                break;
              case L2:
                switch (type) {
                case MBVERTEX:
                  CHKERR getNodesFieldData(*op_data[ss], field_name);
                  break;
                default:
                  break;
                }
                break;
              case NOFIELD:
                if (!getNinTheLoop()) {
                  // NOFIELD data are the same for each element, can be
                  // retrieved only once
                  if (!ss) {
                    CHKERR getNoFieldRowIndices(*op_data[ss], field_name);
                  } else {
                    CHKERR getNoFieldColIndices(*op_data[ss], field_name);
                  }
                  CHKERR getNoFieldFieldData(*op_data[ss], field_name);
                }
                break;
              default:
                SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                         "Not implemented for this space", space);
              }
              last_eval_field_name[ss] = field_name;
            }
          }
        }

        CHKERR swap_bases();

        if (oit->getOpType() & UserDataOperator::OPROW) {
          try {
            CHKERR oit->opRhs(*op_data[0], false);
          }
          CATCH_OP_ERRORS(*oit);
        }

        if (oit->getOpType() & UserDataOperator::OPCOL) {
          try {
            CHKERR oit->opRhs(*op_data[1], false);
          }
          CATCH_OP_ERRORS(*oit);
        }

        if (oit->getOpType() & UserDataOperator::OPROWCOL) {
          try {
            CHKERR oit->opLhs(*op_data[0], *op_data[1]);
          }
          CATCH_OP_ERRORS(*oit);
        }

        CHKERR swap_bases();
      }
    }
    CATCH_OP_ERRORS(*oit);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ForcesAndSourcesCore::UserDataOperator::getProblemRowIndices(
    const std::string field_name, const EntityType type, const int side,
    VectorInt &indices) const {
  MoFEMFunctionBegin;
  if (ptrFE == NULL)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");

  switch (type) {
  case MBVERTEX:
    CHKERR ptrFE->getProblemNodesRowIndices(field_name, indices);
    break;
  default:
    CHKERR ptrFE->getProblemTypeRowIndices(field_name, type, side, indices);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ForcesAndSourcesCore::UserDataOperator::getProblemColIndices(
    const std::string field_name, const EntityType type, const int side,
    VectorInt &indices) const {
  MoFEMFunctionBegin;
  if (ptrFE == NULL)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");

  switch (type) {
  case MBVERTEX:
    CHKERR ptrFE->getProblemNodesColIndices(field_name, indices);
    break;
  default:
    CHKERR ptrFE->getProblemTypeColIndices(field_name, type, side, indices);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::setSideFEPtr(const ForcesAndSourcesCore *side_fe_ptr) {
  MoFEMFunctionBeginHot;
  sidePtrFE = const_cast<ForcesAndSourcesCore *>(side_fe_ptr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ForcesAndSourcesCore::UserDataOperator::loopSide(
    const string &fe_name, ForcesAndSourcesCore *side_fe, const size_t side_dim,
    const EntityHandle ent_for_side) {
  MoFEMFunctionBegin;
  const EntityHandle ent = ent_for_side ? ent_for_side : getFEEntityHandle();

  const Problem *problem_ptr = getFEMethod()->problemPtr;
  Range adjacent_ents;
  CHKERR ptrFE->mField.getInterface<BitRefManager>()->getAdjacenciesAny(
      ent, side_dim, adjacent_ents);
  typedef NumeredEntFiniteElement_multiIndex::index<
      Composite_Name_And_Ent_mi_tag>::type FEByComposite;
  FEByComposite &numered_fe =
      problem_ptr->numeredFiniteElements->get<Composite_Name_And_Ent_mi_tag>();

  side_fe->feName = fe_name;

  CHKERR side_fe->setSideFEPtr(ptrFE);
  CHKERR side_fe->copyBasicMethod(*getFEMethod());
  CHKERR side_fe->copyKsp(*getFEMethod());
  CHKERR side_fe->copySnes(*getFEMethod());
  CHKERR side_fe->copyTs(*getFEMethod());

  CHKERR side_fe->preProcess();

  int nn = 0;
  side_fe->loopSize = adjacent_ents.size();
  for (Range::iterator vit = adjacent_ents.begin(); vit != adjacent_ents.end();
       vit++) {
    FEByComposite::iterator miit =
        numered_fe.find(boost::make_tuple(fe_name, *vit));
    if (miit != numered_fe.end()) {
      side_fe->nInTheLoop = nn++;
      side_fe->numeredEntFiniteElementPtr = *miit;
      CHKERR (*side_fe)();
    }
  }

  CHKERR side_fe->postProcess();

  MoFEMFunctionReturn(0);
}

int ForcesAndSourcesCore::getRule(int order_row, int order_col,
                                  int order_data) {
  return getRuleHook ? getRuleHook(order_row, order_col, order_data)
                     : getRule(order_data);
}

MoFEMErrorCode ForcesAndSourcesCore::setGaussPts(int order_row, int order_col,
                                                 int order_data) {
  return setRuleHook ? setRuleHook(order_row, order_col, order_data)
                     : setGaussPts(order_data);
}

int ForcesAndSourcesCore::getRule(int order) { return 2 * order; }

/** \deprecated setGaussPts(int row_order, int col_order, int data order);
 */
MoFEMErrorCode ForcesAndSourcesCore::setGaussPts(int order) {
  MoFEMFunctionBeginHot;
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Sorry, not implemented");
  MoFEMFunctionReturnHot(0);
}

ForcesAndSourcesCore::UserDataOperator::UserDataOperator(const FieldSpace space,
                                                         const char type,
                                                         const bool symm)
    : DataOperator(symm), opType(type), sPace(space), ptrFE(nullptr) {}

ForcesAndSourcesCore::UserDataOperator::UserDataOperator(
    const std::string &field_name, const char type, const bool symm)
    : DataOperator(symm), opType(type), rowFieldName(field_name),
      colFieldName(field_name), sPace(LASTSPACE), ptrFE(nullptr) {}

ForcesAndSourcesCore::UserDataOperator::UserDataOperator(
    const std::string &row_field_name, const std::string &col_field_name,
    const char type, const bool symm)
    : DataOperator(symm), opType(type), rowFieldName(row_field_name),
      colFieldName(col_field_name), sPace(LASTSPACE), ptrFE(nullptr) {}

MoFEMErrorCode ForcesAndSourcesCore::preProcess() {
  MoFEMFunctionBeginHot;
  if (preProcessHook) {
    ierr = preProcessHook();
    CHKERRG(ierr);
  }
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode ForcesAndSourcesCore::operator()() {
  MoFEMFunctionBeginHot;
  if (operatorHook) {
    ierr = operatorHook();
    CHKERRG(ierr);
  }
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode ForcesAndSourcesCore::postProcess() {
  MoFEMFunctionBeginHot;
  if (postProcessHook) {
    ierr = postProcessHook();
    CHKERRG(ierr);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::UserDataOperator::setPtrFE(ForcesAndSourcesCore *ptr) {
  MoFEMFunctionBeginHot;
  if (!(ptrFE = dynamic_cast<ForcesAndSourcesCore *>(ptr)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "User operator and finite element do not work together");
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM
