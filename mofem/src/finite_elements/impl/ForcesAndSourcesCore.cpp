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
    : mField(m_field), getRuleHook(0) {

  const EntityType type = MBENTITYSET;
  dataOnElement[NOFIELD] = boost::shared_ptr<DataForcesAndSourcesCore>(
      new DataForcesAndSourcesCore(MBENTITYSET));
  derivedDataOnElement[NOFIELD] = boost::shared_ptr<DataForcesAndSourcesCore>(
      new DataForcesAndSourcesCore(MBENTITYSET));

  // Data on elements for proper spaces
  for (int space = H1; space != LASTSPACE; ++space) {
    dataOnElement[space] = boost::shared_ptr<DataForcesAndSourcesCore>(
        new DataForcesAndSourcesCore(type));
    derivedDataOnElement[space] = boost::shared_ptr<DataForcesAndSourcesCore>(
        dataOnElement[space],
        new DerivedDataForcesAndSourcesCore(dataOnElement[space]));
  }

}

ForcesAndSourcesCore::~ForcesAndSourcesCore() {}

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
  if (data.size() < side_table.get<2>().count(type)) {
    // prims has 9 edges, some of edges for "flat" prism are not active
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "data inconsistency %u < %u", data.size(),
             side_table.get<2>().count(type));
  }
  auto siit = side_table.get<2>().lower_bound(type);
  auto hi_siit = side_table.get<2>().upper_bound(type);
  for (; siit != hi_siit; siit++) {
    data[siit->get()->side_number].getSense() = siit->get()->sense;
    if (siit->get()->brother_side_number != -1) {
      if (data.size() < (unsigned)siit->get()->brother_side_number) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
      data[siit->get()->brother_side_number].getSense() = siit->get()->sense;
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getEdgesSense(DataForcesAndSourcesCore &data) const {
  return getSense(MBEDGE, data.dataOnEntities[MBEDGE]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getTrisSense(DataForcesAndSourcesCore &data) const {
  return getSense(MBTRI, data.dataOnEntities[MBTRI]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getQuadSense(DataForcesAndSourcesCore &data) const {
  return getSense(MBQUAD, data.dataOnEntities[MBQUAD]);
}

// ** Order **

template <typename DOFMULTIINDEX>
static int getMaxOrder(const DOFMULTIINDEX &dof_multi_index) {
  int max_order = 0;
  typename DOFMULTIINDEX::iterator dit, hi_dit;
  dit = dof_multi_index.begin();
  hi_dit = dof_multi_index.end();
  for (; dit != hi_dit; dit++) {
    if ((*dit)->getEntDofIdx() != 0)
      continue;
    int dit_max_order = (*dit)->getMaxOrder();
    max_order = (max_order > dit_max_order) ? max_order : dit_max_order;
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
  MoFEMFunctionBeginHot;
  try {
    SideNumber_multiIndex &side_table = const_cast<SideNumber_multiIndex &>(
        numeredEntFiniteElementPtr->getSideNumberTable());
    if (data.size() < side_table.get<2>().count(type)) {
      // prims has 9 edges, some of edges for "flat" prism are not active
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "data inconsistency %d < %d", data.size(),
               side_table.get<2>().count(type));
    }
    for (unsigned int side = 0; side < data.size(); side++) {
      data[side].getDataOrder() = 0;
    }
    FEDofEntity_multiIndex::index<Composite_EntType_and_Space_mi_tag>::type
        &data_dofs = const_cast<FEDofEntity_multiIndex::index<
            Composite_EntType_and_Space_mi_tag>::type &>(
            numeredEntFiniteElementPtr->getDataDofs()
                .get<Composite_EntType_and_Space_mi_tag>());
    FEDofEntity_multiIndex::index<
        Composite_EntType_and_Space_mi_tag>::type::iterator dit,
        hi_dit;
    dit = data_dofs.lower_bound(boost::make_tuple(type, space));
    hi_dit = data_dofs.upper_bound(boost::make_tuple(type, space));
    for (; dit != hi_dit; dit++) {
      ApproximationOrder ent_order = (*dit)->getMaxOrder();
      int side_number = (*dit)->sideNumberPtr->side_number;
      if (side_number < 0) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
      data[side_number].getDataOrder() =
          data[side_number].getDataOrder() > ent_order
              ? data[side_number].getDataOrder()
              : ent_order;
      if ((*dit)->sideNumberPtr->brother_side_number != -1) {
        if (data.size() <
            (unsigned int)(*dit)->sideNumberPtr->brother_side_number) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "data inconsistency");
        }
        data[(*dit)->sideNumberPtr->brother_side_number].getDataOrder() =
            data[side_number].getDataOrder();
      }
    }
  } catch (std::exception &ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__
       << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF, MOFEM_STD_EXCEPTION_THROW, ss.str().c_str());
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getEdgesDataOrder(DataForcesAndSourcesCore &data,
                                        const FieldSpace space) const {
  return getDataOrder(MBEDGE, space, data.dataOnEntities[MBEDGE]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getTrisDataOrder(DataForcesAndSourcesCore &data,
                                       const FieldSpace space) const {
  return getDataOrder(MBTRI, space, data.dataOnEntities[MBTRI]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getQuadDataOrder(DataForcesAndSourcesCore &data,
                                       const FieldSpace space) const {
  return getDataOrder(MBQUAD, space, data.dataOnEntities[MBQUAD]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getTetDataOrder(DataForcesAndSourcesCore &data,
                                      const FieldSpace space) const {
  return getDataOrder(MBTET, space, data.dataOnEntities[MBTET]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getPrismDataOrder(DataForcesAndSourcesCore &data,
                                        const FieldSpace space) const {
  return getDataOrder(MBPRISM, space, data.dataOnEntities[MBPRISM]);
}

MoFEMErrorCode ForcesAndSourcesCore::getDataOrderSpaceAndBase(
    const std::string &field_name, const EntityType type,
    boost::ptr_vector<DataForcesAndSourcesCore::EntData> &data) const {
  auto &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());
  if (data.size() < side_table.get<2>().count(type)) {
    // prims has 9 edges, some of edges for "flat" prism are not active
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "data inconsistency %d < %d", data.size(),
             side_table.get<2>().count(type));
  }

  for (unsigned int side = 0; side != data.size(); ++side) {
    data[side].getDataOrder() = 0;
    data[side].getBase() = NOBASE;
    data[side].getSpace() = NOSPACE;
  }

  // get dofs by name, type and side
  const auto &data_dofs =
      numeredEntFiniteElementPtr->getDataDofs()
          .get<Composite_Name_Type_And_Side_Number_mi_tag>();
  auto dit = data_dofs.lower_bound(boost::make_tuple(field_name, type, 0));
  if (dit == data_dofs.end()) {
    MoFEMFunctionReturnHot(0);
  }
  auto hi_dit =
      data_dofs.upper_bound(boost::make_tuple(field_name, type, data.size()));

  for (; dit != hi_dit;) {

    // std::cerr << ApproximationBaseNames[dit->getApproxBase()] << std::endl;

    int side_number = (*dit)->sideNumberPtr->side_number;
    if (data[side_number].getDataOrder()) {
      const int nb_dofs_on_ent = (*dit)->getNbDofsOnEnt();
      for (int i = 0; i != nb_dofs_on_ent; ++i, ++dit) {
      }
      continue;
    }

    ApproximationOrder ent_order = (*dit)->getMaxOrder();
    if (PetscUnlikely(side_number < 0)) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Side number is not set to dof entity");
    }
    data[side_number].getBase() = (*dit)->getApproxBase();
    data[side_number].getSpace() = (*dit)->getSpace();
    data[side_number].getDataOrder() =
        data[side_number].getDataOrder() > ent_order
            ? data[side_number].getDataOrder()
            : ent_order;
    if ((*dit)->sideNumberPtr->brother_side_number != -1) {
      if (data.size() <
          (unsigned int)(*dit)->sideNumberPtr->brother_side_number) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
      data[(*dit)->sideNumberPtr->brother_side_number].getBase() =
          data[side_number].getBase();
      data[(*dit)->sideNumberPtr->brother_side_number].getSpace() =
          data[side_number].getSpace();
      data[(*dit)->sideNumberPtr->brother_side_number].getDataOrder() =
          data[side_number].getDataOrder();
    }

    const int nb_dofs_on_ent = (*dit)->getNbDofsOnEnt();
    for (int i = 0; i != nb_dofs_on_ent; ++i, ++dit) {
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ForcesAndSourcesCore::getEdgesDataOrderSpaceAndBase(
    DataForcesAndSourcesCore &data, const std::string &field_name) const {
  return getDataOrderSpaceAndBase(field_name, MBEDGE,
                                  data.dataOnEntities[MBEDGE]);
}

MoFEMErrorCode ForcesAndSourcesCore::getTrisDataOrderSpaceAndBase(
    DataForcesAndSourcesCore &data, const std::string &field_name) const {
  return getDataOrderSpaceAndBase(field_name, MBTRI,
                                  data.dataOnEntities[MBTRI]);
}

MoFEMErrorCode ForcesAndSourcesCore::getQuadDataOrderSpaceAndBase(
    DataForcesAndSourcesCore &data, const std::string &field_name) const {
  return getDataOrderSpaceAndBase(field_name, MBQUAD,
                                  data.dataOnEntities[MBQUAD]);
}

MoFEMErrorCode ForcesAndSourcesCore::getTetDataOrderSpaceAndBase(
    DataForcesAndSourcesCore &data, const std::string &field_name) const {
  return getDataOrderSpaceAndBase(field_name, MBTET,
                                  data.dataOnEntities[MBTET]);
}

MoFEMErrorCode ForcesAndSourcesCore::getPrismDataOrderSpaceAndBase(
    DataForcesAndSourcesCore &data, const std::string &field_name) const {
  return getDataOrderSpaceAndBase(field_name, MBPRISM,
                                  data.dataOnEntities[MBPRISM]);
}

// ** Indices **

MoFEMErrorCode ForcesAndSourcesCore::getNodesIndices(
    const boost::string_ref field_name, FENumeredDofEntity_multiIndex &dofs,
    VectorInt &nodes_indices, VectorInt &local_nodes_indices) const {
  MoFEMFunctionBegin;
  auto dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(
      boost::make_tuple(field_name, MBVERTEX));
  auto hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(
      boost::make_tuple(field_name, MBVERTEX));

  int num_nodes;
  CHKERR getNumberOfNodes(num_nodes);
  int max_nb_dofs = 0;
  if (dit != hi_dit) {
    max_nb_dofs = (*dit)->getNbOfCoeffs() * num_nodes;
  }

  if (std::distance(dit, hi_dit) != max_nb_dofs) {
    nodes_indices.resize(max_nb_dofs, false);
    local_nodes_indices.resize(max_nb_dofs, false);
    for (int dd = 0; dd < max_nb_dofs; dd++) {
      nodes_indices[dd] = -1;
      local_nodes_indices[dd] = -1;
    }
  } else {
    nodes_indices.resize(std::distance(dit, hi_dit), false);
    local_nodes_indices.resize(std::distance(dit, hi_dit), false);
  }

  for (; dit != hi_dit; dit++) {
    const int idx = (*dit)->getPetscGlobalDofIdx();
    const int local_idx = (*dit)->getPetscLocalDofIdx();
    const int side_number = (*dit)->sideNumberPtr->side_number;
    const int pos =
        side_number * (*dit)->getNbOfCoeffs() + (*dit)->getDofCoeffIdx();
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
  //
  MoFEMFunctionBeginHot;
  auto dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(
      boost::make_tuple(field_name, type, side_number));
  if (dit == dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().end()) {
    indices.resize(0);
    local_indices.resize(0);
    MoFEMFunctionReturnHot(0);
  }
  auto hi_dit =
      dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(
          boost::make_tuple(field_name, type, side_number));
  if (dit != hi_dit) {
    indices.resize((*dit)->getNbDofsOnEnt(), false);
    local_indices.resize((*dit)->getNbDofsOnEnt(), false);
    for (; dit != hi_dit; dit++) {
      int idx = (*dit)->getPetscGlobalDofIdx();
      int elem_idx = (*dit)->getEntDofIdx();
      indices[elem_idx] = idx;
      int local_idx = (*dit)->getPetscLocalDofIdx();
      local_indices[elem_idx] = local_idx;
    }
  } else {
    indices.resize(0);
    local_indices.resize(0);
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
      data[siiit->get()->side_number].getIndices().resize(0, false);
      data[siiit->get()->side_number].getLocalIndices().resize(0, false);
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
ForcesAndSourcesCore::getEdgesRowIndices(DataForcesAndSourcesCore &data,
                                         const std::string &field_name) const {

  return getTypeIndices(field_name,
                        const_cast<FENumeredDofEntity_multiIndex &>(
                            numeredEntFiniteElementPtr->getRowsDofs()),
                        MBEDGE, data.dataOnEntities[MBEDGE]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getEdgesColIndices(DataForcesAndSourcesCore &data,
                                         const std::string &field_name) const {

  return getTypeIndices(field_name,
                        const_cast<FENumeredDofEntity_multiIndex &>(
                            numeredEntFiniteElementPtr->getColsDofs()),
                        MBEDGE, data.dataOnEntities[MBEDGE]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getTrisRowIndices(DataForcesAndSourcesCore &data,
                                        const std::string &field_name) const {
  return getTypeIndices(field_name,
                        const_cast<FENumeredDofEntity_multiIndex &>(
                            numeredEntFiniteElementPtr->getRowsDofs()),
                        MBTRI, data.dataOnEntities[MBTRI]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getTrisColIndices(DataForcesAndSourcesCore &data,
                                        const std::string &field_name) const {

  return getTypeIndices(field_name,
                        const_cast<FENumeredDofEntity_multiIndex &>(
                            numeredEntFiniteElementPtr->getColsDofs()),
                        MBTRI, data.dataOnEntities[MBTRI]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getTetsRowIndices(DataForcesAndSourcesCore &data,
                                        const std::string &field_name) const {

  MoFEMFunctionBegin;
  if (data.dataOnEntities[MBTET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  CHKERR getTypeIndices(field_name,
                        const_cast<FENumeredDofEntity_multiIndex &>(
                            numeredEntFiniteElementPtr->getRowsDofs()),
                        MBTET, 0, data.dataOnEntities[MBTET][0].getIndices(),
                        data.dataOnEntities[MBTET][0].getLocalIndices());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getTetsColIndices(DataForcesAndSourcesCore &data,
                                        const std::string &field_name) const {

  MoFEMFunctionBegin;
  if (data.dataOnEntities[MBTET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  CHKERR getTypeIndices(field_name,
                        const_cast<FENumeredDofEntity_multiIndex &>(
                            numeredEntFiniteElementPtr->getColsDofs()),
                        MBTET, 0, data.dataOnEntities[MBTET][0].getIndices(),
                        data.dataOnEntities[MBTET][0].getLocalIndices());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getQuadRowIndices(DataForcesAndSourcesCore &data,
                                        const std::string &field_name) const {

  return getTypeIndices(field_name,
                        const_cast<FENumeredDofEntity_multiIndex &>(
                            numeredEntFiniteElementPtr->getRowsDofs()),
                        MBQUAD, data.dataOnEntities[MBQUAD]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getQuadColIndices(DataForcesAndSourcesCore &data,
                                        const std::string &field_name) const {

  return getTypeIndices(field_name,
                        const_cast<FENumeredDofEntity_multiIndex &>(
                            numeredEntFiniteElementPtr->getColsDofs()),
                        MBQUAD, data.dataOnEntities[MBQUAD]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getPrismRowIndices(DataForcesAndSourcesCore &data,
                                         const std::string &field_name) const {

  return getTypeIndices(field_name,
                        const_cast<FENumeredDofEntity_multiIndex &>(
                            numeredEntFiniteElementPtr->getRowsDofs()),
                        MBPRISM, data.dataOnEntities[MBPRISM]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getPrismColIndices(DataForcesAndSourcesCore &data,
                                         const std::string &field_name) const {

  return getTypeIndices(field_name,
                        const_cast<FENumeredDofEntity_multiIndex &>(
                            numeredEntFiniteElementPtr->getColsDofs()),
                        MBPRISM, data.dataOnEntities[MBPRISM]);
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
  auto dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(
      boost::make_tuple(field_name, MBVERTEX));
  auto hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(
      boost::make_tuple(field_name, MBVERTEX));

  int num_nodes;
  CHKERR getNumberOfNodes(num_nodes);
  int max_nb_dofs = 0;
  if (dit != hi_dit) {
    max_nb_dofs = (*dit)->getNbOfCoeffs() * num_nodes;
  }

  if (std::distance(dit, hi_dit) != max_nb_dofs) {
    nodes_data.resize(max_nb_dofs, false);
    nodes_data.clear();
    nodes_dofs.resize(max_nb_dofs, false);
  } else {
    int size = std::distance(dit, hi_dit);
    nodes_data.resize(size, false);
    nodes_dofs.resize(size, false);
  }

  if (dit != hi_dit) {
    space = (*dit)->getSpace();
    base = (*dit)->getApproxBase();
  }

  for (; dit != hi_dit; dit++) {
    FieldData val = (*dit)->getFieldData();
    int side_number = (*dit)->sideNumberPtr->side_number;
    if (side_number == -1) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    int pos = side_number * (*dit)->getNbOfCoeffs() + (*dit)->getDofCoeffIdx();
    nodes_data[pos] = val;
    nodes_dofs[pos] = *dit;
    int brother_side_number = (*dit)->sideNumberPtr->brother_side_number;
    if (brother_side_number != -1) {
      if (nodes_data.size() <
          (unsigned int)(brother_side_number * (*dit)->getNbOfCoeffs() +
                         (*dit)->getNbOfCoeffs())) {
        nodes_data.resize(brother_side_number * (*dit)->getNbOfCoeffs() +
                          (*dit)->getNbOfCoeffs());
        nodes_dofs.resize(brother_side_number * (*dit)->getNbOfCoeffs() +
                          (*dit)->getNbOfCoeffs());
      }
      int brother_pos = brother_side_number * (*dit)->getNbOfCoeffs() +
                        (*dit)->getDofCoeffIdx();
      nodes_data[brother_pos] = val;
      nodes_dofs[brother_pos] = *dit;
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

MoFEMErrorCode ForcesAndSourcesCore::getTypeFieldData(
    const boost::string_ref field_name, FEDofEntity_multiIndex &dofs,
    EntityType type, int side_number, VectorDouble &ent_field_data,
    VectorDofs &ent_field_dofs) const {
  //
  MoFEMFunctionBeginHot;
  auto dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(
      boost::make_tuple(field_name, type, side_number));
  if (dit == dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().end()) {
    ent_field_data.resize(0, false);
    ent_field_dofs.resize(0, false);
    MoFEMFunctionReturnHot(0);
  }
  auto hi_dit =
      dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(
          boost::make_tuple(field_name, type, side_number));
  if (dit != hi_dit) {
    ent_field_data.resize((*dit)->getNbDofsOnEnt(), false);
    ent_field_dofs.resize((*dit)->getNbDofsOnEnt(), false);
    for (; dit != hi_dit; dit++) {
      const FieldData val = (*dit)->getFieldData();
      const int idx = (*dit)->getEntDofIdx();
      ent_field_data[idx] = val;
      ent_field_dofs[idx] = *dit;
    }
  } else {
    ent_field_data.resize(0, false);
    ent_field_dofs.resize(0, false);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ForcesAndSourcesCore::getTypeFieldData(
    const boost::string_ref field_name, FEDofEntity_multiIndex &dofs,
    EntityType type,
    boost::ptr_vector<DataForcesAndSourcesCore::EntData> &data) const {
  MoFEMFunctionBegin;
  auto &side_table = const_cast<SideNumber_multiIndex &>(
      numeredEntFiniteElementPtr->getSideNumberTable());
  // if(data.size() < side_table.get<2>().count(type)) {
  //   SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  // }
  auto siit = side_table.get<2>().lower_bound(type);
  auto hi_siit = side_table.get<2>().upper_bound(type);
  auto tuple = boost::make_tuple(field_name, type);
  auto dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(tuple);
  if (dit == dofs.get<Composite_Name_And_Type_mi_tag>().end()) {
    for (auto siiit = siit; siiit != hi_siit; siiit++) {
      const int side_number = siiit->get()->side_number;
      data[side_number].getFieldData().resize(0, false);
      data[side_number].getFieldDofs().resize(0, false);
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
    const int nb_dofs_on_ent = (*dit)->getNbDofsOnEnt();
    auto &ent_field_data = data[side].getFieldData();
    auto &ent_field_dofs = data[side].getFieldDofs();
    if (!data[side].semaphore) {
      data[side].semaphore = true;
      ent_field_data.resize(nb_dofs_on_ent, false);
      ent_field_dofs.resize(nb_dofs_on_ent, false);
    }
    if (!nb_dofs_on_ent) {
      continue;
    }
    const int idx = dit->get()->getEntDofIdx();
    ent_field_data[idx] = dit->get()->getFieldData();
    ent_field_dofs[idx] = *dit;
  }
  for (; siit != hi_siit; siit++) {
    const int side_number = siit->get()->side_number;
    if (!data[side_number].semaphore) {
      data[side_number].getFieldData().resize(0, false);
      data[side_number].getFieldDofs().resize(0, false);
    }
    if (siit->get()->brother_side_number != -1) {
      if (data.size() < (unsigned int)siit->get()->brother_side_number) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
      CHKERR getTypeFieldData(
          field_name, dofs, type, side_number,
          data[siit->get()->brother_side_number].getFieldData(),
          data[siit->get()->brother_side_number].getFieldDofs());
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

MoFEMErrorCode
ForcesAndSourcesCore::getEdgesFieldData(DataForcesAndSourcesCore &data,
                                        const std::string &field_name) const {
  return getTypeFieldData(field_name,
                          const_cast<FEDofEntity_multiIndex &>(
                              numeredEntFiniteElementPtr->getDataDofs()),
                          MBEDGE, data.dataOnEntities[MBEDGE]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getTrisFieldData(DataForcesAndSourcesCore &data,
                                       const std::string &field_name) const {
  return getTypeFieldData(field_name,
                          const_cast<FEDofEntity_multiIndex &>(
                              numeredEntFiniteElementPtr->getDataDofs()),
                          MBTRI, data.dataOnEntities[MBTRI]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getQuadFieldData(DataForcesAndSourcesCore &data,
                                       const std::string &field_name) const {
  return getTypeFieldData(field_name,
                          const_cast<FEDofEntity_multiIndex &>(
                              numeredEntFiniteElementPtr->getDataDofs()),
                          MBQUAD, data.dataOnEntities[MBQUAD]);
}

MoFEMErrorCode
ForcesAndSourcesCore::getTetsFieldData(DataForcesAndSourcesCore &data,
                                       const std::string &field_name) const {
  MoFEMFunctionBegin;
  if (data.dataOnEntities[MBTET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  CHKERR getTypeFieldData(field_name,
                          const_cast<FEDofEntity_multiIndex &>(
                              numeredEntFiniteElementPtr->getDataDofs()),
                          MBTET, 0,
                          data.dataOnEntities[MBTET][0].getFieldData(),
                          data.dataOnEntities[MBTET][0].getFieldDofs());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
ForcesAndSourcesCore::getPrismFieldData(DataForcesAndSourcesCore &data,
                                        const std::string &field_name) const {
  MoFEMFunctionBegin;
  if (data.dataOnEntities[MBPRISM].size() == 0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  CHKERR getTypeFieldData(field_name,
                          const_cast<FEDofEntity_multiIndex &>(
                              numeredEntFiniteElementPtr->getDataDofs()),
                          MBPRISM, 0,
                          data.dataOnEntities[MBPRISM][0].getFieldData(),
                          data.dataOnEntities[MBPRISM][0].getFieldDofs());
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

MoFEMErrorCode ForcesAndSourcesCore::calculateBaseFunctionsOnElement() {
  MoFEMFunctionBegin;
  /// Use the some node base. Node base is usually used for construction other
  /// bases.

  for (int space = HCURL; space != LASTSPACE; ++space) {
    dataOnElement[space]->dataOnEntities[MBVERTEX][0].getNSharedPtr(
        NOBASE) =
        dataOnElement[H1]->dataOnEntities[MBVERTEX][0].getNSharedPtr(
            NOBASE);
  }

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
    if (dataOnElement[H1]->bAse.test(b)) {
      switch (ApproximationBaseArray[b]) {
      case NOBASE:
        break;
      case AINSWORTH_LEGENDRE_BASE:
      case AINSWORTH_LOBATTO_BASE:
      case DEMKOWICZ_JACOBI_BASE:
        for (int space = H1; space != LASTSPACE; ++space) {
          if (dataOnElement[H1]->sPace.test(space) &&
              dataOnElement[H1]->bAse.test(b) &&
              dataOnElement[H1]->basesOnSpaces[space].test(b)) {
            CHKERR getElementPolynomialBase()->getValue(
                gaussPts,
                boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                    *dataOnElement[space], static_cast<FieldSpace>(space),
                    ApproximationBaseArray[b], NOBASE)));
          }
        }
        break;
      case USER_BASE:
        for (int space = H1; space != LASTSPACE; ++space)
          if (dataOnElement[H1]->sPace.test(space) &&
              dataOnElement[H1]->bAse.test(b) &&
              dataOnElement[H1]->basesOnSpaces[space].test(b)) {
            CHKERR getUserPolynomialBase()->getValue(
                gaussPts,
                boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                    *dataOnElement[space], static_cast<FieldSpace>(space),
                    ApproximationBaseArray[b], NOBASE)));
          }
        break;
      default:
        SETERRQ1(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                 "Base <%s> not yet implemented",
                 ApproximationBaseNames[ApproximationBaseArray[b]]);
      }
    }
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
     << boost::typeindex::type_id_runtime(OP).pretty_name() << ") "

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
      case H1:
      case HCURL:
      case HDIV:
      case L2:
      case NOFIELD:
        op_data[0] = dataOnElement[oit->sPace].get();
        break;
      case LASTSPACE:
        SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY, "unknown space");
        break;
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
          case H1:
          case HCURL:
          case HDIV:
          case L2:
          case NOFIELD:
            op_data[ss] = !ss ? dataOnElement[space[ss]].get()
                              : derivedDataOnElement[space[ss]].get();
            break;
          case LASTSPACE:
            SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                    "unknown space");
            break;
          }

          if (last_eval_field_name[ss] != field_name) {

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
                CHKERR getEdgesRowIndices(*op_data[ss], field_name);
              } else {
                CHKERR getEdgesColIndices(*op_data[ss], field_name);
              }
              CHKERR getEdgesDataOrderSpaceAndBase(*op_data[ss], field_name);
              CHKERR getEdgesFieldData(*op_data[ss], field_name);
              if (dim == 1)
                break;
            case HDIV:
              if (!ss) {
                CHKERR getTrisRowIndices(*op_data[ss], field_name);
              } else {
                CHKERR getTrisColIndices(*op_data[ss], field_name);
              }
              CHKERR getTrisDataOrderSpaceAndBase(*op_data[ss], field_name);
              CHKERR getTrisFieldData(*op_data[ss], field_name);
              if (dim == 2)
                break;
            case L2:
              switch (type) {
              case MBPRISM:
                if (!ss) {
                  CHKERR getPrismRowIndices(*op_data[ss], field_name);
                } else {
                  CHKERR getPrismColIndices(*op_data[ss], field_name);
                }
                CHKERR getPrismDataOrderSpaceAndBase(*op_data[ss], field_name);
                CHKERR getPrismFieldData(*op_data[ss], field_name);
                break;
              case MBTET:
                if (!ss) {
                  CHKERR getTetsRowIndices(*op_data[ss], field_name);
                } else {
                  CHKERR getTetsColIndices(*op_data[ss], field_name);
                }
                CHKERR getTetDataOrderSpaceAndBase(*op_data[ss], field_name);
                CHKERR getTetsFieldData(*op_data[ss], field_name);
                break;
              case MBTRI:
                if (!ss) {
                  CHKERR getTrisRowIndices(*op_data[ss], field_name);
                } else {
                  CHKERR getTrisColIndices(*op_data[ss], field_name);
                }
                CHKERR getTrisDataOrderSpaceAndBase(*op_data[ss], field_name);
                CHKERR getTrisFieldData(*op_data[ss], field_name);
                break;
              case MBEDGE:
                if (!ss) {
                  CHKERR getEdgesRowIndices(*op_data[ss], field_name);
                } else {
                  CHKERR getEdgesColIndices(*op_data[ss], field_name);
                }
                CHKERR getEdgesDataOrderSpaceAndBase(*op_data[ss], field_name);
                CHKERR getEdgesFieldData(*op_data[ss], field_name);
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
            case LASTSPACE:
              SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                      "unknown space");
              break;
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
