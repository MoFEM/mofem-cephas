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

#include <Common.hpp>
#include <Includes.hpp>
#include <definitions.h>
#include <version.h>

#include <fem_tools.h>
#include <h1_hdiv_hcurl_l2.h>

#include <UnknownInterface.hpp>

#include <AdjacencyMultiIndices.hpp>
#include <BCData.hpp>
#include <BCMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <DofsMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <MaterialBlocks.hpp>
#include <ProblemsMultiIndices.hpp>
#include <SeriesMultiIndices.hpp>
#include <TagMultiIndices.hpp>

#include <Core.hpp>
#include <Interface.hpp>
#include <LoopMethods.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>

#include <BaseFunction.hpp>
#include <DataOperators.hpp>
#include <DataStructures.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <FTensor.hpp>

#include <ForcesAndSourcesCore.hpp>

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

      mField(m_field),
      dataOnElement{

          nullptr,
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)), // NOFIELD
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)), // H1
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)), // HCURL
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)), // HDIV
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DataForcesAndSourcesCore(MBENTITYSET)) // L2

      },
      derivedDataOnElement{

          nullptr,
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnElement[NOFIELD])),
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnElement[H1])),
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnElement[HCURL])),
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnElement[HDIV])),
          boost::shared_ptr<DataForcesAndSourcesCore>(
              new DerivedDataForcesAndSourcesCore(dataOnElement[L2]))

      },
      dataNoField(*dataOnElement[NOFIELD].get()),
      dataH1(*dataOnElement[H1].get()), dataHcurl(*dataOnElement[HCURL].get()),
      dataHdiv(*dataOnElement[HDIV].get()), dataL2(*dataOnElement[L2].get()),
      getRuleHook(0), lastEvaluatedElementEntityType(MBMAXTYPE) {}

MoFEMErrorCode ForcesAndSourcesCore::getNumberOfNodes(int &num_nodes) const {
  MoFEMFunctionBeginHot;

  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  switch (mField.get_moab().type_from_handle(ent)) {
  case MBVERTEX:
    num_nodes = 1;
    break;
  case MBEDGE:
    num_nodes = 2;
    break;
  case MBTRI:
    num_nodes = 3;
    break;
  case MBQUAD:
    num_nodes = 4;
    break;
  case MBTET:
    num_nodes = 4;
    break;
  case MBPRISM:
    num_nodes = 6;
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }

  MoFEMFunctionReturnHot(0);
}

// ** Sense **

MoFEMErrorCode ForcesAndSourcesCore::getSense(
    EntityType type,
    boost::ptr_vector<DataForcesAndSourcesCore::EntData> &data) const {
  MoFEMFunctionBegin;
  auto &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());
  if (PetscUnlikely(data.size() < side_table.get<2>().count(type))) {
    // prims has 9 edges, some of edges for "flat" prism are not active
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "wrong number of sides %u < %u", data.size(),
             side_table.get<2>().count(type));
  }
  const auto &st = side_table.get<2>();
  auto sit = st.lower_bound(type);
  auto hi_sit = st.upper_bound(type);
  for (; sit != hi_sit; sit++) {
    const auto &side = **sit;
    data[side.side_number].getSense() = sit->get()->sense;
    if (side.brother_side_number != -1) {
      if (PetscUnlikely(data.size() < (unsigned)side.brother_side_number)) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data struture too small to keep data about brother node");
      }
      data[side.brother_side_number].getSense() = side.sense;
    }
  }
  MoFEMFunctionReturn(0);
}

// ** Order **

template <typename DOFMULTIINDEX>
static int getMaxOrder(const DOFMULTIINDEX &dof_multi_index) {
  auto dit = dof_multi_index.begin();
  auto hi_dit = dof_multi_index.end();
  int max_order = 0;
  for (; dit != hi_dit; dit++) {
    if ((*dit)->getEntDofIdx() == 0) {
      const int dit_max_order = (*dit)->getMaxOrder();
      max_order = (max_order > dit_max_order) ? max_order : dit_max_order;
    }
  }
  return max_order;
}

int ForcesAndSourcesCore::getMaxDataOrder() const {
  return getMaxOrder(numeredEntFiniteElementPtr->getDataDofs());
}

int ForcesAndSourcesCore::getMaxRowOrder() const {
  return getMaxOrder(numeredEntFiniteElementPtr->getRowsDofs());
}

int ForcesAndSourcesCore::getMaxColOrder() const {
  return getMaxOrder(numeredEntFiniteElementPtr->getColsDofs());
}

MoFEMErrorCode ForcesAndSourcesCore::getDataOrder(
    const EntityType type, const FieldSpace space,
    boost::ptr_vector<DataForcesAndSourcesCore::EntData> &data) const {
  MoFEMFunctionBegin;
  SideNumber_multiIndex &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());
  if (PetscUnlikely(data.size() < side_table.get<2>().count(type))) {
    // prims has 9 edges, some of edges for "flat" prism are not active
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "data structure too small to keep data %d < %d", data.size(),
             side_table.get<2>().count(type));
  }
  for (unsigned int side = 0; side != data.size(); ++side) {
    data[side].getDataOrder() = 0;
  }
  auto &data_dofs = const_cast<FEDofEntity_multiIndex::index<
      Composite_EntType_and_Space_mi_tag>::type &>(
      numeredEntFiniteElementPtr->getDataDofs()
          .get<Composite_EntType_and_Space_mi_tag>());
  auto tuple = boost::make_tuple(type, space);
  auto dit = data_dofs.lower_bound(tuple);
  if (dit != data_dofs.end()) {
    auto hi_dit = data_dofs.upper_bound(tuple);
    for (; dit != hi_dit; dit++) {
      auto &dof = **dit;
      ApproximationOrder ent_order = dof.getMaxOrder();
      int side_number = dof.sideNumberPtr->side_number;
      if (PetscUnlikely(side_number < 0)) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "side number can not be negative");
      }
      auto &dat = data[side_number];
      dat.getDataOrder() =
          dat.getDataOrder() > ent_order ? dat.getDataOrder() : ent_order;
      if (dof.sideNumberPtr->brother_side_number != -1) {
        if (PetscUnlikely(
                data.size() <
                (unsigned int)(*dit)->sideNumberPtr->brother_side_number)) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "data size not big enough to keep brother data");
        }
        data[dof.sideNumberPtr->brother_side_number].getDataOrder() =
            dat.getDataOrder();
      }
    }
  }
  MoFEMFunctionReturn(0);
}

// ** Indices **

MoFEMErrorCode ForcesAndSourcesCore::getNodesIndices(
    const boost::string_ref field_name, FENumeredDofEntity_multiIndex &dofs,
    VectorInt &nodes_indices, VectorInt &local_nodes_indices) const {
  MoFEMFunctionBegin;
  auto &dofs_by_type = dofs.get<Composite_Name_And_Type_mi_tag>();
  auto tuple = boost::make_tuple(field_name, MBVERTEX);
  auto dit = dofs_by_type.lower_bound(tuple);
  auto hi_dit = dofs_by_type.upper_bound(tuple);

  int num_nodes;
  CHKERR getNumberOfNodes(num_nodes);
  int max_nb_dofs = 0;
  if (dit != hi_dit) {
    max_nb_dofs = (*dit)->getNbOfCoeffs() * num_nodes;
  }

  const int nb_dofs = std::distance(dit, hi_dit);
  if (nb_dofs != max_nb_dofs) {
    nodes_indices.resize(max_nb_dofs, false);
    local_nodes_indices.resize(max_nb_dofs, false);
    for (int dd = 0; dd < max_nb_dofs; dd++) {
      nodes_indices[dd] = -1;
      local_nodes_indices[dd] = -1;
    }
  } else {
    nodes_indices.resize(nb_dofs, false);
    local_nodes_indices.resize(nb_dofs, false);
  }

  for (; dit != hi_dit; dit++) {
    auto &dof = **dit;
    const int idx = dof.getPetscGlobalDofIdx();
    const int local_idx = dof.getPetscLocalDofIdx();
    const int side_number = dof.sideNumberPtr->side_number;
    const int pos = side_number * dof.getNbOfCoeffs() + dof.getDofCoeffIdx();
    nodes_indices[pos] = idx;
    local_nodes_indices[pos] = local_idx;
    const int brother_side_number = (*dit)->sideNumberPtr->brother_side_number;
    if (brother_side_number != -1) {
      if (nodes_indices.size() <
          (unsigned int)(brother_side_number * (*dit)->getNbOfCoeffs() +
                         (*dit)->getNbOfCoeffs())) {
        nodes_indices.resize(brother_side_number * (*dit)->getNbOfCoeffs() +
                             (*dit)->getNbOfCoeffs());
      }
      const int elem_idx = brother_side_number * (*dit)->getNbOfCoeffs() +
                           (*dit)->getDofCoeffIdx();
      nodes_indices[elem_idx] = idx;
      local_nodes_indices[elem_idx] = local_idx;
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getRowNodesIndices(DataForcesAndSourcesCore &data,
                                         const std::string &field_name) const {
  return getNodesIndices(field_name,
                         const_cast<FENumeredDofEntity_multiIndex &>(
                             numeredEntFiniteElementPtr->getRowsDofs()),
                         data.dataOnEntities[MBVERTEX][0].getIndices(),
                         data.dataOnEntities[MBVERTEX][0].getLocalIndices());
}

MoFEMErrorCode
ForcesAndSourcesCore::getColNodesIndices(DataForcesAndSourcesCore &data,
                                         const std::string &field_name) const {
  return getNodesIndices(field_name,
                         const_cast<FENumeredDofEntity_multiIndex &>(
                             numeredEntFiniteElementPtr->getColsDofs()),
                         data.dataOnEntities[MBVERTEX][0].getIndices(),
                         data.dataOnEntities[MBVERTEX][0].getLocalIndices());
}

MoFEMErrorCode ForcesAndSourcesCore::getTypeIndices(
    const boost::string_ref field_name, FENumeredDofEntity_multiIndex &dofs,
    EntityType type, int side_number, VectorInt &indices,
    VectorInt &local_indices) const {
  MoFEMFunctionBeginHot;
  auto &dofs_by_side = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>();
  auto tuple = boost::make_tuple(field_name, type, side_number);
  auto dit = dofs_by_side.lower_bound(tuple);
  if (dit == dofs_by_side.end()) {
    indices.resize(0);
    local_indices.resize(0);
    MoFEMFunctionReturnHot(0);
  } else {
    auto hi_dit = dofs_by_side.upper_bound(tuple);
    auto &first_dof = **dit;
    indices.resize(first_dof.getNbDofsOnEnt(), false);
    local_indices.resize(first_dof.getNbDofsOnEnt(), false);
    for (; dit != hi_dit; dit++) {
      auto &dof = **dit;
      const int idx = dof.getPetscGlobalDofIdx();
      const int elem_idx = dof.getEntDofIdx();
      indices[elem_idx] = idx;
      const int local_idx = dof.getPetscLocalDofIdx();
      local_indices[elem_idx] = local_idx;
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ForcesAndSourcesCore::getTypeIndices(
    const boost::string_ref field_name, FENumeredDofEntity_multiIndex &dofs,
    EntityType type,
    boost::ptr_vector<DataForcesAndSourcesCore::EntData> &data) const {
  MoFEMFunctionBegin;
  SideNumber_multiIndex &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());
  auto siit = side_table.get<2>().lower_bound(type);
  auto hi_siit = side_table.get<2>().upper_bound(type);
  const auto tuple = boost::make_tuple(field_name, type);
  auto dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(tuple);
  if (dit == dofs.get<Composite_Name_And_Type_mi_tag>().end()) {
    for (SideNumber_multiIndex::nth_index<2>::type::iterator siiit = siit;
         siiit != hi_siit; siiit++) {
      auto &dat = data[siiit->get()->side_number];
      dat.getIndices().resize(0, false);
      dat.getLocalIndices().resize(0, false);
    }
    MoFEMFunctionReturnHot(0);
  }
  for (auto siiit = siit; siiit != hi_siit; siiit++) {
    const int side_number = siiit->get()->side_number;
    data[side_number].semaphore = false;
  }
  auto hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(tuple);
  for (; dit != hi_dit; dit++) {
    const int side = dit->get()->sideNumberPtr->side_number;
    const int nb_dofs_on_ent = dit->get()->getNbDofsOnEnt();
    VectorInt &indices = data[side].getIndices();
    VectorInt &local_indices = data[side].getLocalIndices();
    if (!data[side].semaphore) {
      data[side].semaphore = true;
      indices.resize(nb_dofs_on_ent, false);
      local_indices.resize(nb_dofs_on_ent, false);
    }
    if (!nb_dofs_on_ent) {
      continue;
    }
    const int idx = dit->get()->getEntDofIdx();
    indices[idx] = dit->get()->getPetscGlobalDofIdx();
    local_indices[idx] = dit->get()->getPetscLocalDofIdx();
  }
  for (; siit != hi_siit; siit++) {
    const int side_number = siit->get()->side_number;
    if (!data[side_number].semaphore) {
      data[side_number].getIndices().resize(0, false);
      data[side_number].getLocalIndices().resize(0, false);
    }
    if (siit->get()->brother_side_number != -1) {
      CHKERR getTypeIndices(
          field_name, dofs, type, side_number,
          data[siit->get()->brother_side_number].getIndices(),
          data[siit->get()->brother_side_number].getLocalIndices());
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getNoFieldIndices(const std::string &field_name,
                                        FENumeredDofEntity_multiIndex &dofs,
                                        VectorInt &indices) const {
  MoFEMFunctionBeginHot;
  auto dit = dofs.get<FieldName_mi_tag>().lower_bound(field_name);
  auto hi_dit = dofs.get<FieldName_mi_tag>().upper_bound(field_name);
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
                               numeredEntFiniteElementPtr->getRowsDofs()),
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
                               numeredEntFiniteElementPtr->getColsDofs()),
                           data.dataOnEntities[MBENTITYSET][0].getIndices());
  MoFEMFunctionReturn(0);
}

// ** Indices from problem **

MoFEMErrorCode ForcesAndSourcesCore::getProblemNodesIndices(
    const std::string &field_name, const NumeredDofEntity_multiIndex &dofs,
    VectorInt &nodes_indices) const {
  MoFEMFunctionBeginHot;
  nodes_indices.resize(0);
  auto &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());
  auto siit = side_table.get<1>().lower_bound(boost::make_tuple(MBVERTEX, 0));
  auto hi_siit =
      side_table.get<1>().lower_bound(boost::make_tuple(MBVERTEX, 10000));

  int nn = 0;
  for (; siit != hi_siit; siit++, nn++) {

    if (siit->get()->side_number == -1)
      continue;

    const EntityHandle ent = siit->get()->ent;
    auto dit =
        dofs.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().lower_bound(
            boost::make_tuple(field_name, ent, 0));
    auto hi_dit =
        dofs.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().upper_bound(
            boost::make_tuple(field_name, ent, 10000)); /// very large number
    if (dit != hi_dit) {
      if (!nn) {
        nodes_indices.resize((*dit)->getNbOfCoeffs() *
                             std::distance(siit, hi_siit));
      }
      for (; dit != hi_dit; dit++) {
        nodes_indices[siit->get()->side_number * (*dit)->getNbOfCoeffs() +
                      (*dit)->getDofCoeffIdx()] =
            (*dit)->getPetscGlobalDofIdx();
      }
    }
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
    NumeredDofEntity_multiIndex::index<
        Composite_Name_And_Ent_And_EntDofIdx_mi_tag>::type::iterator dit,
        hi_dit;
    dit = dofs.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().lower_bound(
        boost::make_tuple(field_name, ent, 0));
    hi_dit =
        dofs.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().upper_bound(
            boost::make_tuple(field_name, ent, 10000)); /// very large number

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

MoFEMErrorCode ForcesAndSourcesCore::getNodesFieldData(
    const boost::string_ref field_name, FEDofEntity_multiIndex &dofs,
    VectorDouble &nodes_data, VectorDofs &nodes_dofs, FieldSpace &space,
    FieldApproximationBase &base) const {
  MoFEMFunctionBegin;
  auto tuple = boost::make_tuple(field_name, MBVERTEX);
  auto &dofs_by_name_and_type = dofs.get<Composite_Name_And_Type_mi_tag>();
  auto dit = dofs_by_name_and_type.lower_bound(tuple);
  auto hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(tuple);

  auto &first_dof = **dit;
  space = first_dof.getSpace();
  base = first_dof.getApproxBase();

  const int nb_dofs = std::distance(dit, hi_dit);
  int num_nodes;
  CHKERR getNumberOfNodes(num_nodes);
  const int max_nb_dofs = first_dof.getNbOfCoeffs() * num_nodes;

  if (nb_dofs != max_nb_dofs) {
    nodes_data.resize(max_nb_dofs, false);
    nodes_data.clear();
    nodes_dofs.resize(max_nb_dofs, false);
  } else {
    nodes_data.resize(nb_dofs, false);
    nodes_dofs.resize(nb_dofs, false);
  }

  for (; dit != hi_dit; dit++) {
    const auto &dof_ptr = *dit;
    const auto &dof = *dof_ptr;
    const auto &sn = *dof.sideNumberPtr;
    const int side_number = sn.side_number;
    if (PetscUnlikely(side_number == -1)) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Side number not set");
    }
    const int pos = side_number * dof.getNbOfCoeffs() + dof.getDofCoeffIdx();
    nodes_data[pos] = dof.getFieldData();
    nodes_dofs[pos] = dof_ptr;
    const int brother_side_number = sn.brother_side_number;
    if (brother_side_number != -1) {
      if (nodes_data.size() <
          (unsigned int)(brother_side_number * dof.getNbOfCoeffs() +
                         dof.getNbOfCoeffs())) {
        nodes_data.resize(brother_side_number * dof.getNbOfCoeffs() +
                          dof.getNbOfCoeffs());
        nodes_dofs.resize(brother_side_number * dof.getNbOfCoeffs() +
                          dof.getNbOfCoeffs());
      }
      int brother_pos =
          brother_side_number * dof.getNbOfCoeffs() + dof.getDofCoeffIdx();
      nodes_data[brother_pos] = nodes_data[pos];
      nodes_dofs[brother_pos] = dof_ptr;
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getNodesFieldData(DataForcesAndSourcesCore &data,
                                        const std::string &field_name) const {
  return getNodesFieldData(field_name,
                           const_cast<FEDofEntity_multiIndex &>(
                               numeredEntFiniteElementPtr->getDataDofs()),
                           data.dataOnEntities[MBVERTEX][0].getFieldData(),
                           data.dataOnEntities[MBVERTEX][0].getFieldDofs(),
                           data.dataOnEntities[MBVERTEX][0].getSpace(),
                           data.dataOnEntities[MBVERTEX][0].getBase());
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
  auto &dofs_by_type = dofs.get<Composite_Name_And_Type_mi_tag>();
  auto diit = dofs_by_type.lower_bound(boost::make_tuple(field_name, type_lo));
  if (diit == dofs_by_type.end())
    MoFEMFunctionReturnHot(0);
  auto hi_diit =
      dofs_by_type.lower_bound(boost::make_tuple(field_name, type_hi));
  FEDofEntity_multiIndex_uid_view dof_uid_view;
  dof_uid_view.insert(diit, hi_diit);
  auto dit = dof_uid_view.begin();
  auto hi_dit  = dof_uid_view.end();
  for (; dit != hi_dit;) {
    auto &dof = **dit;
    const EntityType type = dof.getEntType();
    const int side = dof.sideNumberPtr->side_number;
    auto &dat = data.dataOnEntities[type][side];
    const int nb_dofs_on_ent = dof.getNbDofsOnEnt();

    if (nb_dofs_on_ent) {

      dat.getBase() = dof.getApproxBase();
      dat.getSpace() = dof.getSpace();
      const int ent_order = dof.getMaxOrder();
      dat.getDataOrder() =
          dat.getDataOrder() > ent_order ? dat.getDataOrder() : ent_order;
      const int brother_side = dof.sideNumberPtr->brother_side_number;

      auto &ent_field_data = dat.getFieldData();
      auto &ent_field_dofs = dat.getFieldDofs();
      if (ent_field_data.empty()) {
        ent_field_data.resize(nb_dofs_on_ent, false);
        ent_field_dofs.resize(nb_dofs_on_ent, false);
      }
      const auto dof_ent_field_data = dof.getEntFieldData();
      for (int ii = 0; ii != nb_dofs_on_ent; ++ii) {
        ent_field_data[ii] = dof_ent_field_data[ii];
        ent_field_dofs[ii] = *dit;
        ++dit;
      }

      if (brother_side != -1) {
        auto &dat_brother = data.dataOnEntities[type][brother_side];
        dat_brother.getFieldData() = dat.getFieldData();
        dat_brother.getFieldDofs() = dat.getFieldDofs();
        dat_brother.getBase() = dat.getBase();
        dat_brother.getSpace() = dat.getSpace();
        dat_brother.getDataOrder() = dat.getDataOrder();
      }

    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ForcesAndSourcesCore::getNoFieldFieldData(
    const boost::string_ref field_name, FEDofEntity_multiIndex &dofs,
    VectorDouble &ent_field_data, VectorDofs &ent_field_dofs) const {
  MoFEMFunctionBeginHot;
  auto dit = dofs.get<FieldName_mi_tag>().lower_bound(field_name);
  auto hi_dit = dofs.get<FieldName_mi_tag>().upper_bound(field_name);
  int size = std::distance(dit, hi_dit);
  ent_field_data.resize(size, false);
  ent_field_dofs.resize(size, false);
  for (; dit != hi_dit; dit++) {
    int idx = (*dit)->getDofCoeffIdx();
    ent_field_data[idx] = (*dit)->getFieldData();
    ent_field_dofs[idx] = *dit;
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
  for (_IT_GET_FEDATA_DOFS_FOR_LOOP_(this, dof)) {
    if (dof->get()->getEntDofIdx() != 0)
      continue;
    const EntityType type = dof->get()->getEntType();
    const FieldSpace space = dof->get()->getSpace();
    const FieldApproximationBase approx = dof->get()->getApproxBase();
    data.sPace.set(space);
    data.bAse.set(approx);
    data.spacesOnEntities[type].set(space);
    data.basesOnEntities[type].set(approx);
    data.basesOnSpaces[space].set(approx);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ForcesAndSourcesCore::calculateBaseFunctionsOnElement(
    const FieldApproximationBase b) {
  MoFEMFunctionBegin;
  if (dataOnElement[H1]->bAse.test(b)) {
    switch (static_cast<FieldApproximationBase>(b)) {
    case NOBASE:
      break;
    case AINSWORTH_LEGENDRE_BASE:
    case AINSWORTH_LOBATTO_BASE:
    case DEMKOWICZ_JACOBI_BASE:
      if (!getElementPolynomialBase()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Functions genrating approximation base not defined");
      }
      for (int space = H1; space != LASTSPACE; ++space) {
        if (dataOnElement[H1]->sPace.test(space) &&
            dataOnElement[H1]->bAse.test(b) &&
            dataOnElement[H1]->basesOnSpaces[space].test(b)) {
          CHKERR getElementPolynomialBase()->getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                  *dataOnElement[space], static_cast<FieldSpace>(space),
                  static_cast<FieldApproximationBase>(b), NOBASE)));
        }
      }
      break;
    case USER_BASE:
      if (!getUserPolynomialBase()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Functions genrating user approximation base not defined");
      }
      for (int space = H1; space != LASTSPACE; ++space)
        if (dataOnElement[H1]->sPace.test(space) &&
            dataOnElement[H1]->bAse.test(b) &&
            dataOnElement[H1]->basesOnSpaces[space].test(b)) {
          CHKERR getUserPolynomialBase()->getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                  *dataOnElement[space], static_cast<FieldSpace>(space),
                  static_cast<FieldApproximationBase>(b), NOBASE)));
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

MoFEMErrorCode ForcesAndSourcesCore::calculateBaseFunctionsOnElement() {
  MoFEMFunctionBegin;
  /// Use the some node base. Node base is usually used for construction other
  /// bases.
  for (int space = HCURL; space != LASTSPACE; ++space) {
    dataOnElement[space]->dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) =
        dataOnElement[H1]->dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
  }
  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
    CHKERR calculateBaseFunctionsOnElement(
        static_cast<FieldApproximationBase>(b));
  }
  MoFEMFunctionReturn(0);
}

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
  const int dim = mField.get_moab().dimension_from_handle(
      numeredEntFiniteElementPtr->getEnt());
  const UserDataOperator::OpType types[2] = {UserDataOperator::OPROW,
                                             UserDataOperator::OPCOL};
  std::vector<std::string> last_eval_field_name(2);
  DataForcesAndSourcesCore *op_data[2];
  FieldSpace space[2];

  boost::ptr_vector<UserDataOperator>::iterator oit, hi_oit;
  oit = opPtrVector.begin();
  hi_oit = opPtrVector.end();

  for (; oit != hi_oit; oit++) {

    oit->setPtrFE(this);

    if (oit->opType == UserDataOperator::OPLAST) {

      // Set field
      switch (oit->sPace) {
      case NOSPACE:
        SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY, "unknown space");
      case NOFIELD:
      case H1:
      case HCURL:
      case HDIV:
      case L2:
        op_data[0] = dataOnElement[oit->sPace].get();
        break;
      default:
        SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                 "Not implemented for this space", oit->sPace);
      }

      // Reseat all data which all field dependent
      op_data[0]->resetFieldDependentData();
      last_eval_field_name[0] = "";
      last_eval_field_name[1] = "";

      // Run operator
      try {
        CHKERR oit->opRhs(*op_data[0], oit->doVertices, oit->doEdges,
                          oit->doQuads, oit->doTris, oit->doTets, false);
      }
      CATCH_OP_ERRORS(*oit);

    } else {

      for (int ss = 0; ss != 2; ss++) {

        std::string field_name = !ss ? oit->rowFieldName : oit->colFieldName;
        const Field *field_struture = mField.get_field_structure(field_name);
        BitFieldId data_id = field_struture->getId();

        if ((oit->getNumeredEntFiniteElementPtr()->getBitFieldIdData() &
             data_id)
                .none()) {
          SETERRQ2(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                   "no data field < %s > on finite element < %s >",
                   field_name.c_str(), feName.c_str());
        }

        if (oit->getOpType() & types[ss] ||
            oit->getOpType() & UserDataOperator::OPROWCOL) {

          space[ss] = field_struture->getSpace();
          switch (space[ss]) {
          case NOSPACE:
            SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "unknown space");
            break;
          case NOFIELD:
          case H1:
          case HCURL:
          case HDIV:
          case L2:
            op_data[ss] = !ss ? dataOnElement[space[ss]].get()
                              : derivedDataOnElement[space[ss]].get();
            break;
          default:
            SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                     "Not implemented for this space", space[ss]);
          }

          if (last_eval_field_name[ss] != field_name) {

            CHKERR getEntityFieldData(*op_data[ss], field_name, MBEDGE);

            switch (space[ss]) {
            case NOSPACE:
              SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                      "unknown space");
              break;
            case H1:
              if (!ss) {
                CHKERR getRowNodesIndices(*op_data[ss], field_name);
              } else {
                CHKERR getColNodesIndices(*op_data[ss], field_name);
              }
              CHKERR getNodesFieldData(*op_data[ss], field_name);
              if (dim == 0)
                break;
            case HCURL:
              if (!ss) {
                CHKERR getEntityRowIndices<MBEDGE>(*op_data[ss], field_name);
              } else {
                CHKERR getEntityColIndices<MBEDGE>(*op_data[ss], field_name);
              }
              if (dim == 1)
                break;
            case HDIV:
              if (!ss) {
                CHKERR getEntityRowIndices<MBTRI>(*op_data[ss], field_name);
                CHKERR getEntityRowIndices<MBQUAD>(*op_data[ss], field_name);
              } else {
                CHKERR getEntityColIndices<MBTRI>(*op_data[ss], field_name);
                CHKERR getEntityColIndices<MBQUAD>(*op_data[ss], field_name);
              }
              if (dim == 2)
                break;
            case L2:
              switch (type) {
              case MBPRISM:
                if (!ss) {
                  CHKERR getEntityRowIndices<MBPRISM>(*op_data[ss], field_name);
                } else {
                  CHKERR getEntityColIndices<MBPRISM>(*op_data[ss], field_name);
                }
                break;
              case MBTET:
                if (!ss) {
                  CHKERR getEntityRowIndices<MBTET>(*op_data[ss], field_name);
                } else {
                  CHKERR getEntityColIndices<MBTET>(*op_data[ss], field_name);
                }
                break;
              case MBTRI:
                if (!ss) {
                  CHKERR getEntityRowIndices<MBTRI>(*op_data[ss], field_name);
                } else {
                  CHKERR getEntityColIndices<MBTRI>(*op_data[ss], field_name);
                }
                break;
              case MBEDGE:
                if (!ss) {
                  CHKERR getEntityRowIndices<MBEDGE>(*op_data[ss], field_name);
                } else {
                  CHKERR getEntityColIndices<MBEDGE>(*op_data[ss], field_name);
                }
                break;
              case MBVERTEX:
                if (!ss) {
                  CHKERR getRowNodesIndices(*op_data[ss], field_name);
                } else {
                  CHKERR getColNodesIndices(*op_data[ss], field_name);
                }
                CHKERR getNodesFieldData(*op_data[ss], field_name);
                break;
              default:
                SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                        "not implemented");
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
                       "Not implemented for this space", space[ss]);
            }
            last_eval_field_name[ss] = field_name;
          }
        }
      }

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
          CHKERR oit->opLhs(*op_data[0], *op_data[1], oit->sYmm);
        }
        CATCH_OP_ERRORS(*oit);
      }
    }
  }
  MoFEMFunctionReturn(0);
}

// **** Data Operator ****

static MoFEMErrorCode
get_porblem_row_indices(const ForcesAndSourcesCore *fe_ptr,
                        const EntityType type, const int side,
                        const std::string field_name, VectorInt &indices) {
  MoFEMFunctionBegin;
  switch (type) {
  case MBVERTEX:
    CHKERR fe_ptr->getProblemNodesRowIndices(field_name, indices);
    break;
  default:
    CHKERR fe_ptr->getProblemTypeRowIndices(field_name, type, side, indices);
  }
  MoFEMFunctionReturn(0);
}

static MoFEMErrorCode
get_porblem_col_indices(const ForcesAndSourcesCore *fe_ptr,
                        const EntityType type, const int side,
                        const std::string field_name, VectorInt &indices) {
  MoFEMFunctionBegin;
  switch (type) {
  case MBVERTEX:
    CHKERR fe_ptr->getProblemNodesColIndices(field_name, indices);
    break;
  default:
    CHKERR fe_ptr->getProblemTypeColIndices(field_name, type, side, indices);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ForcesAndSourcesCore::UserDataOperator::getProblemRowIndices(
    const std::string field_name, const EntityType type, const int side,
    VectorInt &indices) const {
  MoFEMFunctionBegin;
  if (ptrFE == NULL) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  CHKERR get_porblem_row_indices(ptrFE, type, side, field_name, indices);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ForcesAndSourcesCore::UserDataOperator::getProblemColIndices(
    const std::string field_name, const EntityType type, const int side,
    VectorInt &indices) const {
  MoFEMFunctionBegin;
  if (ptrFE == NULL) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  CHKERR get_porblem_col_indices(ptrFE, type, side, field_name, indices);
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
