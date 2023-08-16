/** \file ForcesAndSourcesCore.cpp

\brief Implementation of Elements on Entities for Forces and Sources
*/

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

static auto cmp_uid_lo(const boost::weak_ptr<FieldEntity> &a, const UId &b) {
  if (auto a_ptr = a.lock()) {
    if (a_ptr->getLocalUniqueId() < b)
      return true;
    else
      return false;
  } else {
    return false;
  }
}

static auto cmp_uid_hi(const UId &b, const boost::weak_ptr<FieldEntity> &a) {
  if (auto a_ptr = a.lock()) {
    if (b < a_ptr->getLocalUniqueId())
      return true;
    else
      return false;
  } else {
    return true;
  }
}

ForcesAndSourcesCore::ForcesAndSourcesCore(Interface &m_field)
    :

      mField(m_field), getRuleHook(0), setRuleHook(0),
      dataOnElement{

          boost::make_shared<EntitiesFieldData>(MBENTITYSET), // NOSPACE,
          boost::make_shared<EntitiesFieldData>(MBENTITYSET), // NOFIELD
          boost::make_shared<EntitiesFieldData>(MBENTITYSET), // H1
          boost::make_shared<EntitiesFieldData>(MBENTITYSET), // HCURL
          boost::make_shared<EntitiesFieldData>(MBENTITYSET), // HDIV
          boost::make_shared<EntitiesFieldData>(MBENTITYSET)  // L2

      },
      derivedDataOnElement{

          nullptr,
          boost::make_shared<DerivedEntitiesFieldData>(
              dataOnElement[NOFIELD]), // NOFIELD
          boost::make_shared<DerivedEntitiesFieldData>(dataOnElement[H1]), // H1
          boost::make_shared<DerivedEntitiesFieldData>(
              dataOnElement[HCURL]), // HCURL
          boost::make_shared<DerivedEntitiesFieldData>(
              dataOnElement[HDIV]),  // HDIV
          boost::make_shared<DerivedEntitiesFieldData>(dataOnElement[L2]) // L2

      },
      dataNoField(*dataOnElement[NOFIELD].get()),
      dataH1(*dataOnElement[H1].get()), dataHcurl(*dataOnElement[HCURL].get()),
      dataHdiv(*dataOnElement[HDIV].get()), dataL2(*dataOnElement[L2].get()),
      lastEvaluatedElementEntityType(MBMAXTYPE), sidePtrFE(nullptr),
      refinePtrFE(nullptr) {

  dataOnElement[NOSPACE]->dataOnEntities[MBENTITYSET].push_back(
      new EntitiesFieldData::EntData());

  dataOnElement[NOFIELD]->dataOnEntities[MBENTITYSET].push_back(
      new EntitiesFieldData::EntData());
  derivedDataOnElement[NOFIELD]->dataOnEntities[MBENTITYSET].push_back(
      new EntitiesFieldData::EntData());
}

// ** Sense **

MoFEMErrorCode ForcesAndSourcesCore::getEntitySense(
    const EntityType type,
    boost::ptr_vector<EntitiesFieldData::EntData> &data) const {
  MoFEMFunctionBegin;

  auto &side_table = numeredEntFiniteElementPtr->getSideNumberTable().get<0>();
  auto sit = side_table.lower_bound(get_id_for_min_type(type));
  if (sit != side_table.end()) {
    auto hi_sit = side_table.upper_bound(get_id_for_max_type(type));
    for (; sit != hi_sit; ++sit) {
      if (const auto side_number = (*sit)->side_number; side_number >= 0) {
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
    if (auto ptr = e.lock()) {
      const int order = ptr->getMaxOrder();
      max_order = (max_order < order) ? order : max_order;
    }
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
    boost::ptr_vector<EntitiesFieldData::EntData> &data) const {
  MoFEMFunctionBegin;

  auto set_order = [&]() {
    MoFEMFunctionBeginHot;
    auto &side_table = numeredEntFiniteElementPtr->getSideNumberTable();

    for (unsigned int s = 0; s != data.size(); ++s)
      data[s].getOrder() = 0;

    const FieldEntity_vector_view &data_field_ent = getDataFieldEnts();

    for (auto r = fieldsPtr->get<BitFieldId_space_mi_tag>().equal_range(space);
         r.first != r.second; ++r.first) {

      const auto field_bit_number = (*r.first)->getBitNumber();
      const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
          field_bit_number, get_id_for_min_type(type));
      auto lo = std::lower_bound(data_field_ent.begin(), data_field_ent.end(),
                                 lo_uid, cmp_uid_lo);
      if (lo != data_field_ent.end()) {
        const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
            field_bit_number, get_id_for_max_type(type));
        auto hi =
            std::upper_bound(lo, data_field_ent.end(), hi_uid, cmp_uid_hi);
        for (; lo != hi; ++lo) {

          if (auto ptr = lo->lock()) {

            auto &e = *ptr;
            auto sit = side_table.find(e.getEnt());
            if (sit != side_table.end()) {
              auto &side = *sit;
              if (const auto side_number = side->side_number;
                  side_number >= 0) {
                ApproximationOrder ent_order = e.getMaxOrder();
                auto &dat = data[side_number];
                dat.getOrder() =
                    dat.getOrder() > ent_order ? dat.getOrder() : ent_order;
              }
            } else
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                      "Entity on side of the element not found");
          }
        }
      }
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
          data[brother_side_number].getOrder() = data[side_number].getOrder();
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

template <typename EXTRACTOR>
MoFEMErrorCode ForcesAndSourcesCore::getNodesIndices(
    const int bit_number, FieldEntity_vector_view &ents_field,
    VectorInt &nodes_indices, VectorInt &local_nodes_indices,
    EXTRACTOR &&extractor) const {
  MoFEMFunctionBegin;

  auto field_it = fieldsPtr->get<BitFieldId_mi_tag>().find(
      BitFieldId().set(bit_number - 1));
  if (field_it != fieldsPtr->get<BitFieldId_mi_tag>().end()) {

#ifndef NDEBUG
    if ((*field_it)->getBitNumber() != bit_number)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Wrong bit number");
#endif
    const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
        bit_number, get_id_for_min_type<MBVERTEX>());
    auto lo = std::lower_bound(ents_field.begin(), ents_field.end(), lo_uid,
                               cmp_uid_lo);
    if (lo != ents_field.end()) {
      const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
          bit_number, get_id_for_max_type<MBVERTEX>());
      auto hi = std::upper_bound(lo, ents_field.end(), hi_uid, cmp_uid_hi);

      const int num_nodes = getNumberOfNodes();
      const int nb_dofs_on_vert = (*field_it)->getNbOfCoeffs();
      const int max_nb_dofs = nb_dofs_on_vert * num_nodes;

      int nb_dofs = 0;
      for (auto it = lo; it != hi; ++it) {
        if (auto e = it->lock()) {
          if (auto cache = extractor(e).lock()) {
            if (cache->loHi[0] != cache->loHi[1]) {
              nb_dofs += std::distance(cache->loHi[0], cache->loHi[1]);
              break;
            }
          }
        }
      }

      if (nb_dofs) {
        nodes_indices.resize(max_nb_dofs, false);
        local_nodes_indices.resize(max_nb_dofs, false);
      } else {
        nodes_indices.resize(0, false);
        local_nodes_indices.resize(0, false);
      }

      if (nb_dofs != max_nb_dofs) {
        std::fill(nodes_indices.begin(), nodes_indices.end(), -1);
        std::fill(local_nodes_indices.begin(), local_nodes_indices.end(), -1);
      }

      for (auto it = lo; it != hi; ++it) {
        if (auto e = it->lock()) {
          auto side_ptr = e->getSideNumberPtr();
          if (const auto side_number = side_ptr->side_number;
              side_number >= 0) {
            const auto brother_side_number = side_ptr->brother_side_number;
            if (auto cache = extractor(e).lock()) {
              for (auto dit = cache->loHi[0]; dit != cache->loHi[1]; ++dit) {
                auto &dof = **dit;
                const int idx = dof.getPetscGlobalDofIdx();
                const int local_idx = dof.getPetscLocalDofIdx();
                const int pos =
                    side_number * nb_dofs_on_vert + dof.getDofCoeffIdx();
                nodes_indices[pos] = idx;
                local_nodes_indices[pos] = local_idx;
                if (brother_side_number != -1) {
                  const int elem_idx = brother_side_number * nb_dofs_on_vert +
                                       (*dit)->getDofCoeffIdx();
                  nodes_indices[elem_idx] = idx;
                  local_nodes_indices[elem_idx] = local_idx;
                }
              }
            }
          }
        }
      }
    } else {
      nodes_indices.resize(0, false);
      local_nodes_indices.resize(0, false);
    }

  } else {
    nodes_indices.resize(0, false);
    local_nodes_indices.resize(0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getRowNodesIndices(EntitiesFieldData &data,
                                         const int bit_number) const {

  struct Extractor {
    boost::weak_ptr<EntityCacheNumeredDofs>
    operator()(boost::shared_ptr<FieldEntity> &e) {
      return e->entityCacheRowDofs;
    }
  };

  return getNodesIndices(bit_number, getRowFieldEnts(),
                         data.dataOnEntities[MBVERTEX][0].getIndices(),
                         data.dataOnEntities[MBVERTEX][0].getLocalIndices(),
                         Extractor());
}

MoFEMErrorCode
ForcesAndSourcesCore::getColNodesIndices(EntitiesFieldData &data,
                                         const int bit_number) const {

  struct Extractor {
    boost::weak_ptr<EntityCacheNumeredDofs>
    operator()(boost::shared_ptr<FieldEntity> &e) {
      return e->entityCacheColDofs;
    }
  };

  return getNodesIndices(bit_number, getColFieldEnts(),
                         data.dataOnEntities[MBVERTEX][0].getIndices(),
                         data.dataOnEntities[MBVERTEX][0].getLocalIndices(),
                         Extractor());
}

template <typename EXTRACTOR>
MoFEMErrorCode ForcesAndSourcesCore::getEntityIndices(
    EntitiesFieldData &data, const int bit_number,
    FieldEntity_vector_view &ents_field, const EntityType type_lo,
    const EntityType type_hi, EXTRACTOR &&extractor) const {
  MoFEMFunctionBegin;

  for (EntityType t = type_lo; t != type_hi; ++t) {
    for (auto &dat : data.dataOnEntities[t]) {
      dat.getIndices().resize(0, false);
      dat.getLocalIndices().resize(0, false);
    }
  }

  const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
      bit_number, get_id_for_min_type(type_lo));
  auto lo = std::lower_bound(ents_field.begin(), ents_field.end(), lo_uid,
                             cmp_uid_lo);
  if (lo != ents_field.end()) {
    const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
        bit_number, get_id_for_max_type(type_hi));
    auto hi = std::upper_bound(lo, ents_field.end(), hi_uid, cmp_uid_hi);
    if (lo != hi) {

      std::vector<boost::weak_ptr<FieldEntity>> brother_ents_vec;

      for (auto it = lo; it != hi; ++it)
        if (auto e = it->lock()) {

          const EntityType type = e->getEntType();
          auto side_ptr = e->getSideNumberPtr();
          if (const auto side = side_ptr->side_number; side >= 0) {
            const auto nb_dofs_on_ent = e->getNbDofsOnEnt();
            const auto brother_side = side_ptr->brother_side_number;
            auto &dat = data.dataOnEntities[type][side];
            auto &ent_field_indices = dat.getIndices();
            auto &ent_field_local_indices = dat.getLocalIndices();

            ent_field_indices.resize(nb_dofs_on_ent, false);
            ent_field_local_indices.resize(nb_dofs_on_ent, false);
            std::fill(ent_field_indices.data().begin(),
                      ent_field_indices.data().end(), -1);
            std::fill(ent_field_local_indices.data().begin(),
                      ent_field_local_indices.data().end(), -1);

            if (auto cache = extractor(e).lock()) {
              for (auto dit = cache->loHi[0]; dit != cache->loHi[1]; ++dit) {
                const int idx = (*dit)->getEntDofIdx();
                ent_field_indices[idx] = (*dit)->getPetscGlobalDofIdx();
                ent_field_local_indices[idx] = (*dit)->getPetscLocalDofIdx();
              }
            }

            if (brother_side != -1) {
              auto &dat_brother = data.dataOnEntities[type][brother_side];
              dat_brother.getIndices().resize(nb_dofs_on_ent, false);
              dat_brother.getLocalIndices().resize(nb_dofs_on_ent, false);
              noalias(dat_brother.getIndices()) = dat.getIndices();
              noalias(dat_brother.getLocalIndices()) = dat.getLocalIndices();
            }
          }
        }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ForcesAndSourcesCore::getEntityRowIndices(
    EntitiesFieldData &data, const int bit_number, const EntityType type_lo,
    const EntityType type_hi) const {

  struct Extractor {
    boost::weak_ptr<EntityCacheNumeredDofs>
    operator()(boost::shared_ptr<FieldEntity> &e) {
      return e->entityCacheRowDofs;
    }
  };

  return getEntityIndices(data, bit_number, getRowFieldEnts(), type_lo, type_hi,
                          Extractor());
}

MoFEMErrorCode ForcesAndSourcesCore::getEntityColIndices(
    EntitiesFieldData &data, const int bit_number, const EntityType type_lo,
    const EntityType type_hi) const {

  struct Extractor {
    boost::weak_ptr<EntityCacheNumeredDofs>
    operator()(boost::shared_ptr<FieldEntity> &e) {
      return e->entityCacheColDofs;
    }
  };

  return getEntityIndices(data, bit_number, getColFieldEnts(), type_lo, type_hi,
                          Extractor());
}

MoFEMErrorCode ForcesAndSourcesCore::getNoFieldIndices(
    const int bit_number, boost::shared_ptr<FENumeredDofEntity_multiIndex> dofs,
    VectorInt &indices) const {
  MoFEMFunctionBeginHot;
  auto dit = dofs->get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId(bit_number));
  auto hi_dit = dofs->get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId(bit_number));
  indices.resize(std::distance(dit, hi_dit));
  for (; dit != hi_dit; dit++) {
    int idx = (*dit)->getPetscGlobalDofIdx();
    indices[(*dit)->getDofCoeffIdx()] = idx;
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getNoFieldRowIndices(EntitiesFieldData &data,
                                           const int bit_number) const {
  MoFEMFunctionBegin;
#ifndef NDEBUG
  if (data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "data.dataOnEntities[MBENTITYSET] is empty");
  }
#endif
  CHKERR getNoFieldIndices(bit_number, getRowDofsPtr(),
                           data.dataOnEntities[MBENTITYSET][0].getIndices());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getNoFieldColIndices(EntitiesFieldData &data,
                                           const int bit_number) const {
  MoFEMFunctionBegin;
#ifndef NDEBUG
  if (data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "data.dataOnEntities[MBENTITYSET] is empty");
  }
#endif
  CHKERR getNoFieldIndices(bit_number, getColDofsPtr(),
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

    const int num_nodes = getNumberOfNodes();
    nodes_indices.resize(field_struture->getNbOfCoeffs() * num_nodes, false);
    std::fill(nodes_indices.begin(), nodes_indices.end(), -1);

    auto &side_table = const_cast<SideNumber_multiIndex &>(
        numeredEntFiniteElementPtr->getSideNumberTable());
    auto siit = side_table.get<1>().lower_bound(boost::make_tuple(MBVERTEX, 0));
    auto hi_siit =
        side_table.get<1>().lower_bound(boost::make_tuple(MBVERTEX, 10000));

    int nn = 0;
    for (; siit != hi_siit; siit++, nn++) {
      if (siit->get()->side_number >= 0) {
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
    if (siit->get()->side_number >= 0) {

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
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ForcesAndSourcesCore::getProblemNodesRowIndices(
    const std::string &field_name, VectorInt &nodes_indices) const {
  return getProblemNodesIndices(field_name, *(problemPtr->numeredRowDofsPtr),
                                nodes_indices);
}

MoFEMErrorCode
ForcesAndSourcesCore::getProblemTypeRowIndices(const std::string &field_name,
                                               EntityType type, int side_number,
                                               VectorInt &indices) const {
  return getProblemTypeIndices(field_name, *(problemPtr->numeredRowDofsPtr),
                               type, side_number, indices);
}

MoFEMErrorCode ForcesAndSourcesCore::getProblemNodesColIndices(
    const std::string &field_name, VectorInt &nodes_indices) const {
  return getProblemNodesIndices(field_name, *(problemPtr->numeredColDofsPtr),
                                nodes_indices);
}

MoFEMErrorCode
ForcesAndSourcesCore::getProblemTypeColIndices(const std::string &field_name,
                                               EntityType type, int side_number,
                                               VectorInt &indices) const {
  return getProblemTypeIndices(field_name, *(problemPtr->numeredColDofsPtr),
                               type, side_number, indices);
}

// ** Data **

MoFEMErrorCode ForcesAndSourcesCore::getBitRefLevelOnData() {
  MoFEMFunctionBegin;

  for (auto &data : dataOnElement) {
    if (data) {
      for (auto &dat : data->dataOnEntities) {
        for (auto &ent_dat : dat) {
          ent_dat.getEntDataBitRefLevel().clear();
        }
      }
    }
  }

  auto &field_ents = getDataFieldEnts();
  for (auto it : field_ents) {
    if (auto e = it.lock()) {
      const FieldSpace space = e->getSpace();
      if (space > NOFIELD) {
        const EntityType type = e->getEntType();
        const signed char side = e->getSideNumberPtr()->side_number;
        if (side >= 0) {
          if (auto &data = dataOnElement[space]) {
            if (type == MBVERTEX) {
              auto &dat = data->dataOnEntities[type][0];
              dat.getEntDataBitRefLevel().resize(getNumberOfNodes(), false);
              dat.getEntDataBitRefLevel()[side] = e->getBitRefLevel();
            } else {
              auto &dat = data->dataOnEntities[type][side];
              dat.getEntDataBitRefLevel().resize(1, false);
              dat.getEntDataBitRefLevel()[0] = e->getBitRefLevel();
            }
          }
        }
      } else {
        if (auto &data = dataOnElement[NOFIELD]) {
          auto &dat = data->dataOnEntities[MBENTITYSET][0];
          dat.getEntDataBitRefLevel().resize(1, false);
          dat.getEntDataBitRefLevel()[0] = e->getBitRefLevel();
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
};

MoFEMErrorCode
ForcesAndSourcesCore::getNodesFieldData(EntitiesFieldData &data,
                                        const int bit_number) const {

  auto get_nodes_field_data = [&](VectorDouble &nodes_data,
                                  VectorFieldEntities &field_entities,
                                  VectorDofs &nodes_dofs, FieldSpace &space,
                                  FieldApproximationBase &base,
                                  VectorInt &bb_node_order) {
    MoFEMFunctionBegin;

    nodes_data.resize(0, false);
    nodes_dofs.resize(0, false);
    field_entities.resize(0, false);

    auto field_it =
        fieldsPtr->get<BitFieldId_mi_tag>().find(BitFEId().set(bit_number - 1));
    if (field_it != fieldsPtr->get<BitFieldId_mi_tag>().end()) {

#ifndef NDEBUG
      if ((*field_it)->getBitNumber() != bit_number)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Wrong bit number");
#endif
      const int nb_dofs_on_vert = (*field_it)->getNbOfCoeffs();
      space = (*field_it)->getSpace();
      base = (*field_it)->getApproxBase();

      auto &field_ents = getDataFieldEnts();
      const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
          bit_number, get_id_for_min_type<MBVERTEX>());
      auto lo = std::lower_bound(field_ents.begin(), field_ents.end(), lo_uid,
                                 cmp_uid_lo);
      if (lo != field_ents.end()) {
        const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
            bit_number, get_id_for_max_type<MBVERTEX>());
        auto hi = std::upper_bound(lo, field_ents.end(), hi_uid, cmp_uid_hi);
        if (lo != hi) {

          int nb_dofs = 0;
          for (auto it = lo; it != hi; ++it) {
            if (auto e = it->lock()) {
              if (auto cache = e->entityCacheDataDofs.lock()) {
                if (cache->loHi[0] != cache->loHi[1]) {
                  nb_dofs += std::distance(cache->loHi[0], cache->loHi[1]);
                  break;
                }
              }
            }
          }

          if (nb_dofs) {

            const int num_nodes = getNumberOfNodes();
            bb_node_order.resize(num_nodes, false);
            bb_node_order.clear();
            const int max_nb_dofs = nb_dofs_on_vert * num_nodes;
            nodes_data.resize(max_nb_dofs, false);
            nodes_dofs.resize(max_nb_dofs, false);
            field_entities.resize(num_nodes, false);
            std::fill(nodes_data.begin(), nodes_data.end(), 0);
            std::fill(nodes_dofs.begin(), nodes_dofs.end(), nullptr);
            std::fill(field_entities.begin(), field_entities.end(), nullptr);

            std::vector<boost::weak_ptr<FieldEntity>> brother_ents_vec;

            for (auto it = lo; it != hi; ++it) {
              if (auto e = it->lock()) {
                const auto &sn = e->getSideNumberPtr();
                // Some field entities on skeleton can have negative side
                // number
                if (const auto side_number = sn->side_number;
                    side_number >= 0) {
                  const int brother_side_number = sn->brother_side_number;

                  field_entities[side_number] = e.get();
                  if (brother_side_number != -1) {
                    brother_ents_vec.emplace_back(e);
                    field_entities[side_number] = field_entities[side_number];
                  }

                  bb_node_order[side_number] = e->getMaxOrder();
                  int pos = side_number * nb_dofs_on_vert;
                  auto ent_filed_data_vec = e->getEntFieldData();
                  if (auto cache = e->entityCacheDataDofs.lock()) {
                    for (auto dit = cache->loHi[0]; dit != cache->loHi[1];
                         ++dit) {
                      const auto dof_idx = (*dit)->getEntDofIdx();
                      nodes_data[pos + dof_idx] = ent_filed_data_vec[dof_idx];
                      nodes_dofs[pos + dof_idx] =
                          reinterpret_cast<FEDofEntity *>((*dit).get());
                    }
                  }
                }
              }
            }

            for (auto &it : brother_ents_vec) {
              if (const auto e = it.lock()) {
                const auto &sn = e->getSideNumberPtr();
                const int side_number = sn->side_number;
                const int brother_side_number = sn->brother_side_number;
                bb_node_order[brother_side_number] = bb_node_order[side_number];
                int pos = side_number * nb_dofs_on_vert;
                int brother_pos = brother_side_number * nb_dofs_on_vert;
                for (int ii = 0; ii != nb_dofs_on_vert; ++ii) {
                  nodes_data[brother_pos] = nodes_data[pos];
                  nodes_dofs[brother_pos] = nodes_dofs[pos];
                  ++pos;
                  ++brother_pos;
                }
              }
            }
          }
        }
      }
    }

    MoFEMFunctionReturn(0);
  };

  return get_nodes_field_data(
      data.dataOnEntities[MBVERTEX][0].getFieldData(),
      data.dataOnEntities[MBVERTEX][0].getFieldEntities(),
      data.dataOnEntities[MBVERTEX][0].getFieldDofs(),
      data.dataOnEntities[MBVERTEX][0].getSpace(),
      data.dataOnEntities[MBVERTEX][0].getBase(),
      data.dataOnEntities[MBVERTEX][0].getBBNodeOrder());
}

MoFEMErrorCode ForcesAndSourcesCore::getEntityFieldData(
    EntitiesFieldData &data, const int bit_number, const EntityType type_lo,
    const EntityType type_hi) const {
  MoFEMFunctionBegin;
  for (EntityType t = type_lo; t != type_hi; ++t) {
    for (auto &dat : data.dataOnEntities[t]) {
      dat.getOrder() = 0;
      dat.getBase() = NOBASE;
      dat.getSpace() = NOSPACE;
      dat.getFieldData().resize(0, false);
      dat.getFieldDofs().resize(0, false);
      dat.getFieldEntities().resize(0, false);
    }
  }

  auto &field_ents = getDataFieldEnts();
  const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
      bit_number, get_id_for_min_type(type_lo));
  auto lo = std::lower_bound(field_ents.begin(), field_ents.end(), lo_uid,
                             cmp_uid_lo);
  if (lo != field_ents.end()) {
    const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
        bit_number, get_id_for_max_type(type_hi));
    auto hi = std::upper_bound(lo, field_ents.end(), hi_uid, cmp_uid_hi);
    if (lo != hi) {

      std::vector<boost::weak_ptr<FieldEntity>> brother_ents_vec;

      for (auto it = lo; it != hi; ++it)
        if (auto e = it->lock()) {
          auto side_ptr = e->getSideNumberPtr();
          if (const auto side = side_ptr->side_number; side >= 0) {
            const EntityType type = e->getEntType();
            auto &dat = data.dataOnEntities[type][side];
            auto &ent_field = dat.getFieldEntities();
            auto &ent_field_dofs = dat.getFieldDofs();
            auto &ent_field_data = dat.getFieldData();

            const int brother_side = side_ptr->brother_side_number;
            if (brother_side != -1)
              brother_ents_vec.emplace_back(e);

            dat.getBase() = e->getApproxBase();
            dat.getSpace() = e->getSpace();
            const int ent_order = e->getMaxOrder();
            dat.getOrder() =
                dat.getOrder() > ent_order ? dat.getOrder() : ent_order;

            auto ent_data = e->getEntFieldData();
            ent_field_data.resize(ent_data.size(), false);
            noalias(ent_field_data) = ent_data;
            ent_field_dofs.resize(ent_data.size(), false);
            std::fill(ent_field_dofs.begin(), ent_field_dofs.end(), nullptr);
            ent_field.resize(1, false);
            ent_field[0] = e.get();
            if (auto cache = e->entityCacheDataDofs.lock()) {
              for (auto dit = cache->loHi[0]; dit != cache->loHi[1]; ++dit) {
                ent_field_dofs[(*dit)->getEntDofIdx()] =
                    reinterpret_cast<FEDofEntity *>((*dit).get());
              }
            }
          }
        }

      for (auto &it : brother_ents_vec) {
        if (const auto e = it.lock()) {
          const EntityType type = e->getEntType();
          const int side = e->getSideNumberPtr()->side_number;
          const int brother_side = e->getSideNumberPtr()->brother_side_number;
          auto &dat = data.dataOnEntities[type][side];
          auto &dat_brother = data.dataOnEntities[type][brother_side];
          dat_brother.getBase() = dat.getBase();
          dat_brother.getSpace() = dat.getSpace();
          dat_brother.getOrder() = dat.getOrder();
          dat_brother.getFieldData() = dat.getFieldData();
          dat_brother.getFieldDofs() = dat.getFieldDofs();
          dat_brother.getFieldEntities() = dat.getFieldEntities();
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ForcesAndSourcesCore::getNoFieldFieldData(
    const int bit_number, VectorDouble &ent_field_data,
    VectorDofs &ent_field_dofs, VectorFieldEntities &ent_field) const {

  MoFEMFunctionBeginHot;

  ent_field_data.resize(0, false);
  ent_field_dofs.resize(0, false);
  ent_field.resize(0, false);

  auto &field_ents = getDataFieldEnts();
  const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
      bit_number, get_id_for_min_type<MBVERTEX>());
  auto lo = std::lower_bound(field_ents.begin(), field_ents.end(), lo_uid,
                             cmp_uid_lo);
  if (lo != field_ents.end()) {

    ent_field.resize(field_ents.size(), false);
    std::fill(ent_field.begin(), ent_field.end(), nullptr);

    const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
        bit_number, get_id_for_max_type<MBENTITYSET>());
    auto hi = std::upper_bound(lo, field_ents.end(), hi_uid, cmp_uid_hi);
    if (lo != hi) {

      int side = 0;
      for (auto it = lo; it != hi; ++it, ++side)
        if (auto e = it->lock()) {

          const auto size = e->getNbDofsOnEnt();
          ent_field_data.resize(size, false);
          ent_field_dofs.resize(size, false);
          ent_field[side] = e.get();
          noalias(ent_field_data) = e->getEntFieldData();

          if (auto cache = e->entityCacheDataDofs.lock()) {
            for (auto dit = cache->loHi[0]; dit != cache->loHi[1]; ++dit) {
              ent_field_dofs[(*dit)->getEntDofIdx()] =
                  reinterpret_cast<FEDofEntity *>((*dit).get());
            }
          }
        }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getNoFieldFieldData(EntitiesFieldData &data,
                                          const int bit_number) const {
  MoFEMFunctionBegin;
  if (data.dataOnEntities[MBENTITYSET].size() == 0)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "No space to insert data");

  CHKERR getNoFieldFieldData(
      bit_number, data.dataOnEntities[MBENTITYSET][0].getFieldData(),
      data.dataOnEntities[MBENTITYSET][0].getFieldDofs(),
      data.dataOnEntities[MBENTITYSET][0].getFieldEntities());
  MoFEMFunctionReturn(0);
}

// ** Face **

MoFEMErrorCode
ForcesAndSourcesCore::getFaceNodes(EntitiesFieldData &data) const {
  MoFEMFunctionBegin;
  auto &face_nodes = data.facesNodes;
  auto &face_nodes_order = data.facesNodesOrder;
  auto &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());
  const auto ent = numeredEntFiniteElementPtr->getEnt();
  const auto type = numeredEntFiniteElementPtr->getEntType();
  const auto nb_faces = CN::NumSubEntities(type, 2);
  const EntityHandle *conn_ele;
  int num_nodes_ele;
  CHKERR mField.get_moab().get_connectivity(ent, conn_ele, num_nodes_ele, true);
  auto side_ptr_it = side_table.get<1>().lower_bound(
      boost::make_tuple(CN::TypeDimensionMap[2].first, 0));
  auto hi_side_ptr_it = side_table.get<1>().upper_bound(
      boost::make_tuple(CN::TypeDimensionMap[2].second, 100));

  for (; side_ptr_it != hi_side_ptr_it; ++side_ptr_it) {
    const auto side = (*side_ptr_it)->side_number;
    const auto sense = (*side_ptr_it)->sense;
    const auto offset = (*side_ptr_it)->offset;

    EntityType face_type;
    int nb_nodes_face;
    auto face_indices =
        CN::SubEntityVertexIndices(type, 2, side, face_type, nb_nodes_face);
    face_nodes.resize(nb_faces, nb_nodes_face);
    face_nodes_order.resize(nb_faces, nb_nodes_face);

    if (sense == 1)
      for (int n = 0; n != nb_nodes_face; ++n)
        face_nodes_order(side, n) = (n + offset) % nb_nodes_face;
    else
      for (int n = 0; n != nb_nodes_face; ++n)
        face_nodes_order(side, n) =
            (nb_nodes_face - (n - offset) % nb_nodes_face) % nb_nodes_face;

    for (int n = 0; n != nb_nodes_face; ++n)
      face_nodes(side, n) = face_indices[face_nodes_order(side, n)];

#ifndef NDEBUG
    const auto face_ent = (*side_ptr_it)->ent;
    auto check = [&]() {
      MoFEMFunctionBegin;
      const EntityHandle *conn_face;
      // int nb_nodes_face;
      CHKERR mField.get_moab().get_connectivity(face_ent, conn_face,
                                                nb_nodes_face, true);
      face_nodes.resize(nb_faces, nb_nodes_face);
      for (int nn = 0; nn != nb_nodes_face; ++nn) {
        if (face_nodes(side, nn) !=
            std::distance(
                conn_ele,
                std::find(conn_ele, &conn_ele[num_nodes_ele], conn_face[nn]))) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Wrong face numeration");
        }
      }
      MoFEMFunctionReturn(0);
    };
    CHKERR check();
#endif
  }

  MoFEMFunctionReturn(0);
}

// ** Space and Base **

MoFEMErrorCode ForcesAndSourcesCore::getSpacesAndBaseOnEntities(
    EntitiesFieldData &data) const {
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

  auto fe_type = numeredEntFiniteElementPtr->getEntType();

  if (getDataFieldEntsPtr())
    for (auto e : getDataFieldEnts()) {
      if (auto ptr = e.lock()) {
        // get data from entity
        const EntityType type = ptr->getEntType();
        if (DefaultElementAdjacency::getDefTypeMap(fe_type, type)) {
          const FieldSpace space = ptr->getSpace();
          const FieldApproximationBase approx = ptr->getApproxBase();
          // set data
          data.sPace.set(space);
          data.bAse.set(approx);
          data.spacesOnEntities[type].set(space);
          data.basesOnEntities[type].set(approx);
          data.basesOnSpaces[space].set(approx);
        }
      }
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
                "Functions generating user approximation base not defined");

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
    dataOnElement[space]->dataOnEntities[MBVERTEX][0].getDiffNSharedPtr(
        NOBASE) =
        dataOnElement[H1]->dataOnEntities[MBVERTEX][0].getDiffNSharedPtr(
            NOBASE);
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

  const auto ele_type = numeredEntFiniteElementPtr->getEntType();

  auto get_nodal_base_data = [&](EntitiesFieldData &data, auto field_ptr) {
    MoFEMFunctionBegin;
    auto &space = data.dataOnEntities[MBVERTEX][0].getSpace();
    auto &base = data.dataOnEntities[MBVERTEX][0].getBase();
    auto &bb_node_order = data.dataOnEntities[MBVERTEX][0].getBBNodeOrder();

    auto &field_ents = getDataFieldEnts();
    auto bit_number = field_ptr->getBitNumber();
    const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
        bit_number, get_id_for_min_type<MBVERTEX>());
    auto lo = std::lower_bound(field_ents.begin(), field_ents.end(), lo_uid,
                               cmp_uid_lo);
    if (lo != field_ents.end()) {
      const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
          bit_number, get_id_for_max_type<MBVERTEX>());
      auto hi = std::upper_bound(lo, field_ents.end(), hi_uid, cmp_uid_hi);
      if (lo != hi) {

        for (auto it = lo; it != hi; ++it)
          if (auto first_e = it->lock()) {
            space = first_e->getSpace();
            base = first_e->getApproxBase();
            const int num_nodes = getNumberOfNodes();
            bb_node_order.resize(num_nodes, false);
            bb_node_order.clear();

            std::vector<boost::weak_ptr<FieldEntity>> brother_ents_vec;

            for (; it != hi; ++it) {
              if (auto e = it->lock()) {
                const auto &sn = e->getSideNumberPtr();
                const int side_number = sn->side_number;
                const int brother_side_number = sn->brother_side_number;
                if (brother_side_number != -1)
                  brother_ents_vec.emplace_back(e);
                bb_node_order[side_number] = e->getMaxOrder();
              }
            }

            for (auto &it : brother_ents_vec) {
              if (const auto e = it.lock()) {
                const auto &sn = e->getSideNumberPtr();
                const int side_number = sn->side_number;
                const int brother_side_number = sn->brother_side_number;
                bb_node_order[brother_side_number] = bb_node_order[side_number];
              }
            }

            break;
          }
      }
    }

    MoFEMFunctionReturn(0);
  };

  auto get_entity_base_data = [&](EntitiesFieldData &data, auto field_ptr,
                                  const EntityType type_lo,
                                  const EntityType type_hi) {
    MoFEMFunctionBegin;
    for (EntityType t = MBEDGE; t != MBPOLYHEDRON; ++t) {
      for (auto &dat : data.dataOnEntities[t]) {
        dat.getOrder() = 0;
        dat.getBase() = NOBASE;
        dat.getSpace() = NOSPACE;
        dat.getFieldData().resize(0, false);
        dat.getFieldDofs().resize(0, false);
      }
    }

    auto &field_ents = getDataFieldEnts();
    auto bit_number = field_ptr->getBitNumber();
    const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
        bit_number, get_id_for_min_type(type_lo));
    auto lo = std::lower_bound(field_ents.begin(), field_ents.end(), lo_uid,
                               cmp_uid_lo);
    if (lo != field_ents.end()) {
      const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
          bit_number, get_id_for_max_type(type_hi));
      auto hi = std::upper_bound(lo, field_ents.end(), hi_uid, cmp_uid_hi);

      std::vector<boost::weak_ptr<FieldEntity>> brother_ents_vec;
      for (; lo != hi; ++lo) {
        if (auto e = lo->lock()) {
          if (auto cache = e->entityCacheDataDofs.lock()) {
            if (cache->loHi[0] != cache->loHi[1]) {
              if (const auto side = e->getSideNumberPtr()->side_number;
                  side >= 0) {
                const EntityType type = e->getEntType();
                auto &dat = data.dataOnEntities[type][side];
                const int brother_side =
                    e->getSideNumberPtr()->brother_side_number;
                if (brother_side != -1)
                  brother_ents_vec.emplace_back(e);
                dat.getBase() = e->getApproxBase();
                dat.getSpace() = e->getSpace();
                const auto ent_order = e->getMaxOrder();
                dat.getOrder() =
                    dat.getOrder() > ent_order ? dat.getOrder() : ent_order;
              }
            }
          }
        }
      }

      for (auto &ent_ptr : brother_ents_vec) {
        if (auto e = ent_ptr.lock()) {
          const EntityType type = e->getEntType();
          const int side = e->getSideNumberPtr()->side_number;
          const int brother_side = e->getSideNumberPtr()->brother_side_number;
          auto &dat = data.dataOnEntities[type][side];
          auto &dat_brother = data.dataOnEntities[type][brother_side];
          dat_brother.getBase() = dat.getBase();
          dat_brother.getSpace() = dat.getSpace();
          dat_brother.getOrder() = dat.getOrder();
        }
      }
    }
    MoFEMFunctionReturn(0);
  };

  for (auto &e : getDataFieldEnts()) {
    if (auto ent_data_ptr = e.lock()) {
      if (ent_data_ptr->getApproxBase() == AINSWORTH_BERNSTEIN_BEZIER_BASE) {
        auto space = ent_data_ptr->getSpace();
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
  }

  auto check_space = [&](const auto space) {
    switch (space) {
    case H1:
      for (auto t = MBVERTEX; t <= ele_type; ++t) {
        if (dataOnElement[H1]->spacesOnEntities[t].test(H1))
          return true;
      }
      return false;
    case HCURL:
      for (auto t = MBEDGE; t <= ele_type; ++t) {
        if (dataOnElement[HCURL]->spacesOnEntities[t].test(HCURL))
          return true;
      }
      return false;
    case HDIV:
      for (auto t = MBTRI; t <= ele_type; ++t) {
        if (dataOnElement[HDIV]->spacesOnEntities[t].test(HDIV))
          return true;
      }
      return false;
    case L2:
      return dataOnElement[L2]->spacesOnEntities[ele_type].test(L2);
      break;
    default:
      THROW_MESSAGE("Not implemented");
    }
  };

  std::set<string> fields_list;
  for (auto &e : getDataFieldEnts()) {
    if (auto ent_data_ptr = e.lock()) {
      if (ent_data_ptr->getApproxBase() == AINSWORTH_BERNSTEIN_BEZIER_BASE) {
        auto field_name = ent_data_ptr->getName();
        if (fields_list.find(field_name) == fields_list.end()) {
          auto field_ptr = ent_data_ptr->getFieldRawPtr();
          auto space = ent_data_ptr->getSpace();
          CHKERR get_nodal_base_data(*dataOnElement[space], field_ptr);
          CHKERR get_entity_base_data(*dataOnElement[space], field_ptr, MBEDGE,
                                      ele_type);
          if (check_space(space)) {
            CHKERR getElementPolynomialBase()->getValue(
                gaussPts, boost::make_shared<EntPolynomialBaseCtx>(
                              *dataOnElement[space], field_name,
                              static_cast<FieldSpace>(space),
                              AINSWORTH_BERNSTEIN_BEZIER_BASE, NOBASE));
            fields_list.insert(field_name);
          }
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
};

MoFEMErrorCode ForcesAndSourcesCore::createDataOnElement(EntityType type) {
  MoFEMFunctionBegin;

  // Data on elements for proper spaces
  for (int space = H1; space != LASTSPACE; ++space) {
    dataOnElement[space]->setElementType(type);
    derivedDataOnElement[space]->setElementType(type);
  }

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

  using UDO = UserDataOperator;

  std::array<std::string, 2> field_name;
  std::array<const Field *, 3> field_struture;
  std::array<int, 2>
      field_id; // note the this is field bit number (nor field bit)
  std::array<FieldSpace, 2> space;
  std::array<FieldApproximationBase, 2> base;

  constexpr std::array<UDO::OpType, 2> types = {UDO::OPROW, UDO::OPCOL};
  std::array<int, 2> last_eval_field_id = {0, 0};

  std::array<boost::shared_ptr<EntitiesFieldData>, 2> op_data;

  auto swap_bases = [&](auto &op) {
    MoFEMFunctionBeginHot;
    for (size_t ss = 0; ss != 2; ++ss) {
      if (op.getOpType() & types[ss] || op.getOpType() & UDO::OPROWCOL) {
        switch (base[ss]) {
        case AINSWORTH_BERNSTEIN_BEZIER_BASE:
          CHKERR op_data[ss]->baseSwap(field_name[ss],
                                       AINSWORTH_BERNSTEIN_BEZIER_BASE);
        default:
          break;
        }
      }
    }

    MoFEMFunctionReturnHot(0);
  };

  const EntityType type = numeredEntFiniteElementPtr->getEntType();

  // evaluate entity data only, no field specific data provided or known
  auto evaluate_op_space = [&](auto &op) {
    MoFEMFunctionBeginHot;

    // reseat all data which all field dependent
    dataOnElement[op.sPace]->resetFieldDependentData();
    std::fill(last_eval_field_id.begin(), last_eval_field_id.end(), 0);

    switch (op.sPace) {
    case NOSPACE:
      try {
        CHKERR op.doWork(
            0, MBENTITYSET,
            dataOnElement[NOSPACE]->dataOnEntities[MBENTITYSET][0]);
      }
      CATCH_OP_ERRORS(op);
      break;
    case NOFIELD:
    case H1:
    case HCURL:
    case HDIV:
    case L2:
      try {

        for (EntityType t = MBVERTEX; t != MBMAXTYPE; ++t) {
          for (auto &e : dataOnElement[op.sPace]->dataOnEntities[t]) {
            e.getSpace() = op.sPace;
          }
        }

        CHKERR op.opRhs(*dataOnElement[op.sPace], false);
      }
      CATCH_OP_ERRORS(op);
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
               "Not implemented for this space", op.sPace);
    }

    MoFEMFunctionReturnHot(0);
  };

  // set entity data
  auto set_op_entityties_data = [&](auto ss, auto &op) {
    MoFEMFunctionBeginHot;

#ifndef NDEBUG
    if ((op.getNumeredEntFiniteElementPtr()->getBitFieldIdData() &
         mField.get_field_id(field_name[ss]))
            .none()) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "no data field < %s > on finite element < %s >",
               field_name[ss].c_str(), getFEName().c_str());
    }
#endif

    op_data[ss] =
        !ss ? dataOnElement[space[ss]] : derivedDataOnElement[space[ss]];

    for (auto &data : op_data[ss]->dataOnEntities[MBENTITYSET]) {
      CHKERR data.resetFieldDependentData();
    }

    auto get_data_for_nodes = [&]() {
      MoFEMFunctionBegin;
      if (!ss)
        CHKERR getRowNodesIndices(*op_data[ss], field_id[ss]);
      else
        CHKERR getColNodesIndices(*op_data[ss], field_id[ss]);
      CHKERR getNodesFieldData(*op_data[ss], field_id[ss]);
      MoFEMFunctionReturn(0);
    };

    // get data on entities but not nodes
    auto get_data_for_entities = [&]() {
      MoFEMFunctionBegin;
      CHKERR getEntityFieldData(*op_data[ss], field_id[ss], MBEDGE);
      if (!ss)
        CHKERR getEntityRowIndices(*op_data[ss], field_id[ss], MBEDGE);
      else
        CHKERR getEntityColIndices(*op_data[ss], field_id[ss], MBEDGE);
      MoFEMFunctionReturn(0);
    };

    auto get_data_for_meshset = [&]() {
      MoFEMFunctionBegin;
      if (!ss) {
        CHKERR getNoFieldRowIndices(*op_data[ss], field_id[ss]);
      } else {
        CHKERR getNoFieldColIndices(*op_data[ss], field_id[ss]);
      }
      CHKERR getNoFieldFieldData(*op_data[ss], field_id[ss]);
      MoFEMFunctionReturn(0);
    };

    switch (space[ss]) {
    case H1:
      CHKERR get_data_for_nodes();
    case HCURL:
    case HDIV:
      CHKERR get_data_for_entities();
      break;
    case L2:
      switch (type) {
      case MBVERTEX:
        CHKERR get_data_for_nodes();
        break;
      default:
        CHKERR get_data_for_entities();
      }
      break;
    case NOFIELD:
      // if (!getNinTheLoop()) {
      // NOFIELD data are the same for each element, can be
      // retrieved only once
      CHKERR get_data_for_meshset();
      // }
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
               "not implemented for this space < %s >",
               FieldSpaceNames[space[ss]]);
    }
    MoFEMFunctionReturnHot(0);
  };

  // evalate operators with field data
  auto evaluate_op_for_fields = [&](auto &op) {
    MoFEMFunctionBeginHot;

    if (op.getOpType() & UDO::OPROW) {
      try {
        CHKERR op.opRhs(*op_data[0], false);
      }
      CATCH_OP_ERRORS(op);
    }

    if (op.getOpType() & UDO::OPCOL) {
      try {
        CHKERR op.opRhs(*op_data[1], false);
      }
      CATCH_OP_ERRORS(op);
    }

    if (op.getOpType() & UDO::OPROWCOL) {
      try {
        CHKERR op.opLhs(*op_data[0], *op_data[1]);
      }
      CATCH_OP_ERRORS(op);
    }
    MoFEMFunctionReturnHot(0);
  };

  // Collect bit ref level on all entities, fields and spaces
  CHKERR getBitRefLevelOnData();

  auto oit = opPtrVector.begin();
  auto hi_oit = opPtrVector.end();

  // interate over all operators on element
  for (; oit != hi_oit; oit++) {

    try {

      CHKERR oit->setPtrFE(this);

      if ((oit->opType & UDO::OPSPACE) == UDO::OPSPACE) {

        // run operator for space but specific field
        CHKERR evaluate_op_space(*oit);

      } else if (

          (oit->opType & (UDO::OPROW | UDO::OPCOL | UDO::OPROWCOL)) ==
          oit->opType

      ) {

        field_name[0] = oit->rowFieldName;
        field_name[1] = oit->colFieldName;

        for (size_t ss = 0; ss != 2; ss++) {

#ifndef NDEBUG
          if (field_name[ss].empty()) {
            SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                     "Not set Field name in operator %d (0-row, 1-column) in "
                     "operator %s",
                     ss,
                     (boost::typeindex::type_id_runtime(*oit).pretty_name())
                         .c_str());
          }
#endif

          field_struture[ss] = mField.get_field_structure(field_name[ss]);
          field_id[ss] = field_struture[ss]->getBitNumber();
          space[ss] = field_struture[ss]->getSpace();
          base[ss] = field_struture[ss]->getApproxBase();
        }

        // not that if field name do not change between operators, entity field
        // data are nor rebuild
        for (size_t ss = 0; ss != 2; ss++) {

          if (oit->getOpType() & types[ss] ||
              oit->getOpType() & UDO::OPROWCOL) {
            if (last_eval_field_id[ss] != field_id[ss]) {
              CHKERR set_op_entityties_data(ss, *oit);
              last_eval_field_id[ss] = field_id[ss];
            }
          }
        }

        CHKERR swap_bases(*oit);

        // run operator for given field or row, column or both
        CHKERR evaluate_op_for_fields(*oit);

        CHKERR swap_bases(*oit);

      } else {

        for (int i = 0; i != 5; ++i)
          if (oit->opType & (1 << i))
            MOFEM_LOG("SELF", Sev::error) << UDO::OpTypeNames[i];
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Impossible operator type");
      }
    }
    CATCH_OP_ERRORS(*oit);
  }
  MoFEMFunctionReturn(0);
}

const char *const ForcesAndSourcesCore::UserDataOperator::OpTypeNames[] = {
    "OPROW", " OPCOL", "OPROWCOL", "OPSPACE", "OPLAST"};

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

MoFEMErrorCode ForcesAndSourcesCore::setRefineFEPtr(
    const ForcesAndSourcesCore *refine_fe_ptr) {
  MoFEMFunctionBeginHot;
  refinePtrFE = const_cast<ForcesAndSourcesCore *>(refine_fe_ptr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ForcesAndSourcesCore::UserDataOperator::loopSide(
    const string &fe_name, ForcesAndSourcesCore *side_fe, const size_t side_dim,
    const EntityHandle ent_for_side, const int verb,
    const LogManager::SeverityLevel sev, AdjCache *adj_cache) {
  MoFEMFunctionBegin;

  const auto *problem_ptr = getFEMethod()->problemPtr;
  auto &numered_fe =
      problem_ptr->numeredFiniteElementsPtr->get<Unique_mi_tag>();

  auto fe_miit = ptrFE->mField.get_finite_elements()
                     ->get<FiniteElement_name_mi_tag>()
                     .find(fe_name);
  if (fe_miit != ptrFE->mField.get_finite_elements()
                     ->get<FiniteElement_name_mi_tag>()
                     .end()) {

    const auto ent = ent_for_side ? ent_for_side : getFEEntityHandle();

    side_fe->feName = fe_name;

    CHKERR side_fe->setSideFEPtr(ptrFE);
    CHKERR side_fe->copyBasicMethod(*getFEMethod());
    CHKERR side_fe->copyPetscData(*getFEMethod());
    CHKERR side_fe->copyKsp(*getFEMethod());
    CHKERR side_fe->copySnes(*getFEMethod());
    CHKERR side_fe->copyTs(*getFEMethod());

    side_fe->cacheWeakPtr = getFEMethod()->cacheWeakPtr;

    CHKERR side_fe->preProcess();

    std::vector<boost::weak_ptr<NumeredEntFiniteElement>> fe_vec;
    auto get_numered_fe_ptr = [&](auto &fe_uid, Range &&adjacent_ents)
        -> std::vector<boost::weak_ptr<NumeredEntFiniteElement>> & {
      fe_vec.reserve(adjacent_ents.size());
      for (auto fe_ent : adjacent_ents) {
        auto miit = numered_fe.find(
            EntFiniteElement::getLocalUniqueIdCalculate(fe_ent, fe_uid));
        if (miit != numered_fe.end()) {
          fe_vec.emplace_back(*miit);
        }
      }
      return fe_vec;
    };

    auto get_bit_entity_adjacency = [&]() {
      Range adjacent_ents;
      CHKERR ptrFE->mField.getInterface<BitRefManager>()->getAdjacenciesAny(
          ent, side_dim, adjacent_ents);
      return adjacent_ents;
    };

    auto get_adj = [&](auto &fe_uid)
        -> std::vector<boost::weak_ptr<NumeredEntFiniteElement>> & {
      if (adj_cache) {
        try {
          return (*adj_cache).at(ent);
        } catch (const std::out_of_range &) {
          return (*adj_cache)[ent] =
                     get_numered_fe_ptr(fe_uid, get_bit_entity_adjacency());
        }
      } else {
        return get_numered_fe_ptr(fe_uid, get_bit_entity_adjacency());
      }
    };

    auto adj = get_adj((*fe_miit)->getFEUId());

    int nn = 0;
    side_fe->loopSize = adj.size();
    for (auto fe_weak_ptr : adj) {
      if (auto fe_ptr = fe_weak_ptr.lock()) {
        if (verb >= VERBOSE)
          MOFEM_LOG("SELF", sev) << "Side finite element "
                                 << "(" << nn << "): " << *fe_ptr;
        side_fe->nInTheLoop = nn;
        side_fe->numeredEntFiniteElementPtr = fe_ptr;
        CHKERR (*side_fe)();
      }
      ++nn;
    }

    CHKERR side_fe->postProcess();
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ForcesAndSourcesCore::UserDataOperator::loopThis(
    const string &fe_name, ForcesAndSourcesCore *this_fe, const int verb,
    const LogManager::SeverityLevel sev) {
  MoFEMFunctionBegin;

  if (verb >= VERBOSE)
    MOFEM_LOG("SELF", sev) << "This finite element: "
                           << *getNumeredEntFiniteElementPtr();

  this_fe->feName = fe_name;

  CHKERR this_fe->setRefineFEPtr(ptrFE);
  CHKERR this_fe->copyBasicMethod(*getFEMethod());
  CHKERR this_fe->copyPetscData(*getFEMethod());
  CHKERR this_fe->copyKsp(*getFEMethod());
  CHKERR this_fe->copySnes(*getFEMethod());
  CHKERR this_fe->copyTs(*getFEMethod());

  this_fe->cacheWeakPtr = getFEMethod()->cacheWeakPtr;

  CHKERR this_fe->preProcess();

  this_fe->nInTheLoop = getNinTheLoop();
  this_fe->numeredEntFiniteElementPtr = getNumeredEntFiniteElementPtr();

  CHKERR (*this_fe)();

  CHKERR this_fe->postProcess();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ForcesAndSourcesCore::UserDataOperator::loopParent(
    const string &fe_name, ForcesAndSourcesCore *parent_fe, const int verb,
    const LogManager::SeverityLevel sev) {
  MoFEMFunctionBegin;

  auto &fes =
      ptrFE->mField.get_finite_elements()->get<FiniteElement_name_mi_tag>();
  auto fe_miit = fes.find(fe_name);
  if (fe_miit != fes.end()) {

    const auto *problem_ptr = getFEMethod()->problemPtr;
    auto &numered_fe =
        problem_ptr->numeredFiniteElementsPtr->get<Unique_mi_tag>();


    parent_fe->feName = fe_name;
    CHKERR parent_fe->setRefineFEPtr(ptrFE);
    CHKERR parent_fe->copyBasicMethod(*getFEMethod());
    CHKERR parent_fe->copyPetscData(*getFEMethod());
    CHKERR parent_fe->copyKsp(*getFEMethod());
    CHKERR parent_fe->copySnes(*getFEMethod());
    CHKERR parent_fe->copyTs(*getFEMethod());

    parent_fe->cacheWeakPtr = getFEMethod()->cacheWeakPtr;

    const auto parent_ent = getNumeredEntFiniteElementPtr()->getParentEnt();
    auto miit = numered_fe.find(EntFiniteElement::getLocalUniqueIdCalculate(
        parent_ent, (*fe_miit)->getFEUId()));
    if (miit != numered_fe.end()) {
      if (verb >= VERBOSE)
        MOFEM_LOG("SELF", sev) << "Parent finite element: " << **miit;
      parent_fe->loopSize = 1;
      parent_fe->nInTheLoop = 0;
      parent_fe->numeredEntFiniteElementPtr = *miit;
      CHKERR parent_fe->preProcess();
      CHKERR (*parent_fe)();
      CHKERR parent_fe->postProcess();
    } else {
      if (verb >= VERBOSE)
        MOFEM_LOG("SELF", sev) << "Parent finite element: no parent";
      parent_fe->loopSize = 0;
      parent_fe->nInTheLoop = 0;
      CHKERR parent_fe->preProcess();
      CHKERR parent_fe->postProcess();
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ForcesAndSourcesCore::UserDataOperator::loopChildren(
    const string &fe_name, ForcesAndSourcesCore *child_fe, const int verb,
    const LogManager::SeverityLevel sev) {
  MoFEMFunctionBegin;

  auto fe_miit = ptrFE->mField.get_finite_elements()
                     ->get<FiniteElement_name_mi_tag>()
                     .find(fe_name);
  if (fe_miit != ptrFE->mField.get_finite_elements()
                     ->get<FiniteElement_name_mi_tag>()
                     .end()) {

    const auto *problem_ptr = getFEMethod()->problemPtr;
    auto &ref_ents = *getPtrFE()->mField.get_ref_ents();
    auto &numered_fe =
        problem_ptr->numeredFiniteElementsPtr->get<Unique_mi_tag>();

    const auto parent_ent = getNumeredEntFiniteElementPtr()->getEnt();
    const auto parent_type = getNumeredEntFiniteElementPtr()->getEntType();
    auto range =
        ref_ents.get<Composite_ParentEnt_And_EntType_mi_tag>().equal_range(
            boost::make_tuple(parent_type, parent_ent));

    if (auto size = std::distance(range.first, range.second)) {

      std::vector<EntityHandle> childs_vec;
      childs_vec.reserve(size);
      for (; range.first != range.second; ++range.first)
        childs_vec.emplace_back((*range.first)->getEnt());

      Range childs;

      if ((childs_vec.back() - childs_vec.front() + 1) == size)
        childs = Range(childs_vec.front(), childs_vec.back());
      else
        childs.insert_list(childs_vec.begin(), childs_vec.end());

      child_fe->feName = fe_name;

      CHKERR child_fe->setRefineFEPtr(ptrFE);
      CHKERR child_fe->copyBasicMethod(*getFEMethod());
      CHKERR child_fe->copyPetscData(*getFEMethod());
      CHKERR child_fe->copyKsp(*getFEMethod());
      CHKERR child_fe->copySnes(*getFEMethod());
      CHKERR child_fe->copyTs(*getFEMethod());

      child_fe->cacheWeakPtr = getFEMethod()->cacheWeakPtr;

      CHKERR child_fe->preProcess();

      int nn = 0;
      child_fe->loopSize = size;

      for (auto p = childs.pair_begin(); p != childs.pair_end(); ++p) {

        auto miit =
            numered_fe.lower_bound(EntFiniteElement::getLocalUniqueIdCalculate(
                p->first, (*fe_miit)->getFEUId()));
        auto hi_miit =
            numered_fe.upper_bound(EntFiniteElement::getLocalUniqueIdCalculate(
                p->second, (*fe_miit)->getFEUId()));

        for (; miit != hi_miit; ++miit) {

          if (verb >= VERBOSE)
            MOFEM_LOG("SELF", sev) << "Child finite element "
                                   << "(" << nn << "): " << **miit;

          child_fe->nInTheLoop = nn++;
          child_fe->numeredEntFiniteElementPtr = *miit;
          CHKERR (*child_fe)();
        }
      }
    }

    CHKERR child_fe->postProcess();
  }

  MoFEMFunctionReturn(0);
}

int ForcesAndSourcesCore::getRule(int order_row, int order_col,
                                  int order_data) {
  return getRuleHook ? getRuleHook(order_row, order_col, order_data)
                     : getRule(order_data);
}

MoFEMErrorCode ForcesAndSourcesCore::setGaussPts(int order_row, int order_col,
                                                 int order_data) {
  return setRuleHook ? setRuleHook(this, order_row, order_col, order_data)
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
    const std::string field_name, const char type, const bool symm)
    : DataOperator(symm), opType(type), rowFieldName(field_name),
      colFieldName(field_name), sPace(LASTSPACE), ptrFE(nullptr) {}

ForcesAndSourcesCore::UserDataOperator::UserDataOperator(
    const std::string row_field_name, const std::string col_field_name,
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
  } else {
#ifndef NDEBUG
    MOFEM_LOG("SELF", Sev::warning)
        << "No method operator() overloaded on element entity on finite "
           "element <"
        << boost::typeindex::type_id_runtime(*this).pretty_name() << ">";
#endif
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
